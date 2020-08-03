#include <ATen/BatchedFallback.h>
#include <ATen/VmapTransforms.h>

namespace at {

// Given a linear index, return the actual index.
// Example: Given linear_idx = 3, sizes = [5, 2], we would return [1, 0]
static SmallVector<indexing::TensorIndex,kVmapStaticDimVecSize>
computeIndex(int64_t linear_idx, IntArrayRef sizes) {
  SmallVector<indexing::TensorIndex,kVmapStaticDimVecSize> result;
  result.reserve(sizes.size());
  for (auto it = sizes.rbegin(); it != sizes.rend(); it++) {
    auto remainder = linear_idx % *it;
    result.push_back(remainder);
    linear_idx -= remainder;
    linear_idx /= *it;
  }
  std::reverse(std::begin(result), std::end(result));
  return result;
}

static bool areAllReturnsTensors(const FunctionSchema& schema) {
  return std::all_of(
      schema.returns().begin(),
      schema.returns().end(),
      [] (const Argument& arg) { return arg.type() == TensorType::get(); });
}

static bool areAnyArgumentsTensorList(const FunctionSchema& schema) {
  return std::any_of(
      schema.arguments().begin(),
      schema.arguments().end(),
      [] (const Argument& arg) { return arg.type()->isSubtypeOf(ListType::ofTensors()); });
}

// The general flow of the algorithm is as follows.
// - First, we figure out which arguments are BatchedTensors and save them
//   to a vector. We also store a vector of which index of the arguments list
//   each BatchedTensor appears in. This will be useful for bookkeeping later.
// - Next, we apply the MultiBatchVmapTransform to all of the BatchedTensors.
//   This returns a vector of VmapPhysicalView that hold tensors that contain
//   all of the collective batch dimensions at the front of the tensors.
// - Then, we attempt to call `op` once per slice of the inputs. To do this,
//   we repeatedly we slice the input arguments (if they are BatchedTensors),
//   put the sliced (or a not-sliced) version of the input onto the stack, invoke
//   the operator, and then pop the results off the stack.
// - Each result obtained from the previous step is a slice of the total result,
//   so we stack those tensors together to form the final result.
void batchedTensorForLoopFallback(const c10::OperatorHandle& op, torch::jit::Stack* stack) {
  const auto& schema = op.schema();
  const auto num_returns = schema.returns().size();
  TORCH_CHECK(!schema.is_mutable() && !schema.hasAnyAliasInfo(),
              "Batching rule not implemented for ", schema, "; ",
              "the fallback path doesn't work on in-place or view ops.");
  TORCH_CHECK(areAllReturnsTensors(schema) && !areAnyArgumentsTensorList(schema),
              "Batching rule not implemented for ", schema, ". ",
              "We could not generate a fallback.");
  TORCH_CHECK(num_returns == 1,
              "Batching rule not implemented for ", schema, ". ",
              "We do not yet support operations with multiple returns.");
  TORCH_WARN("Batching rule not implemented for ", schema, " falling back "
             "to slow (for loop and stack) implementation");

  const auto num_arguments = schema.arguments().size();
  const auto arguments = torch::jit::last(stack, num_arguments);
  const auto arguments_begin = stack->size() - num_arguments;

  // Figure out which arguments are BatchedTensor. Save them to a vector.
  // For each BatchedTensor, also record what position of `arguments` they came from.
  SmallVector<Tensor,kVmapTransformStaticInputSize> batched_tensor_inputs;
  VmapDimVector batched_tensor_inputs_position;
  for (int64_t idx = 0; idx < arguments.size(); ++idx) {
    const auto& ivalue = arguments[idx];
    if (!ivalue.isTensor()) {
      continue;
    }
    const auto& tensor = ivalue.toTensor();
    if (!tensor.defined()) {
      continue;
    }
    const auto* batched = maybeGetBatched(tensor);
    if (!batched) {
      continue;
    }
    batched_tensor_inputs.push_back(tensor);
    batched_tensor_inputs_position.push_back(idx);
  }
  TORCH_INTERNAL_ASSERT(batched_tensor_inputs.size() > 0);

  // MultiBatchVmapTransform the BatchedTensor arguments. This returns
  // VmapPhysicalViews that contain all of the batch dimensions.
  const auto input_physical_views = MultiBatchVmapTransform::logicalToPhysical(
      batched_tensor_inputs);

  // Compute the total number of batches
  auto num_batch_dims = input_physical_views.front().numBatchDims();
  auto some_sizes = input_physical_views.front().tensor().sizes();
  auto batch_sizes = ArrayRef<int64_t>(some_sizes.begin(), some_sizes.begin() + num_batch_dims);
  auto num_batches = std::accumulate(
      batch_sizes.begin(),
      batch_sizes.end(),
      1,
      std::multiplies<int64_t>());

  // Strategy: For each batch, we are going to push slices (where applicable)
  // of the arguments onto `stack`, call `op`, and store the result in
  // `output_shards`.
  std::vector<Tensor> output_shards;
  output_shards.reserve(num_batches);

  for (int64_t linear_idx = 0; linear_idx < num_batches; ++linear_idx) {
    auto index = computeIndex(linear_idx, batch_sizes);
    auto batched_tensor_inputs_pos_iter = batched_tensor_inputs_position.begin();
    auto input_physical_views_iter = input_physical_views.begin();
    for (int64_t arg_idx = 0; arg_idx < num_arguments; ++arg_idx) {
      // We assume that torch::jit::Stack is backed by vector<IValue> for
      // simplicity. When that is not the case, this code should be updated.
      const auto& argument = (*stack)[arguments_begin + arg_idx];
      if (arg_idx != *batched_tensor_inputs_pos_iter) {
        // argument isn't a BatchedTensor
        torch::jit::push(stack, argument);
        continue;
      }
      // argument is a BatchedTensor
      const auto& physical_view_for_argument = *input_physical_views_iter;
      torch::jit::push(stack, physical_view_for_argument.tensor().index(index));
      batched_tensor_inputs_pos_iter++;
      input_physical_views_iter++;
    }

    op.callBoxed(stack);

    // We assume there is a single tensor return
    output_shards.emplace_back(torch::jit::pop(stack).toTensor());
  }

  // Stack the tensors together to form the result.
  auto flat_output = at::stack(output_shards);
  VmapDimVector output_sizes(batch_sizes);
  output_sizes.insert(
      output_sizes.end(),
      flat_output.sizes().begin() + 1,
      flat_output.sizes().end());
  torch::jit::drop(stack, num_arguments);
  torch::jit::push(
      stack,
      input_physical_views.front().newLogicalFromPhysical(flat_output.view(output_sizes)));
}

} // namespace at
