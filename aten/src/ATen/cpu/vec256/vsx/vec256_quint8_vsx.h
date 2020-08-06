#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>
#include <ATen/cpu/vec256/vsx/vsx_helpers.h>
#include <c10/util/quint8.h>
#include <array>

// This file defines Vec256<> for the quantized types.
//
//
// Currently, we simply use these classes as efficient converters between
// the quantized types and Vec256<float>, usually in bandwidth-bound cases
// where doing the arithmetic in full-precision is acceptable (e.g.
// elementwise operators).
//
//
// Conversions are as follows:
//  Vec256<quint8> -> 4x Vec256<float>
//
// The size of the returned float vector is specified by the special
// constexpr function float_num_vecs. The type of the value returned
// from dequantize (and expected as an argument to quantize) is
// specified by float_vec_return_type.
//
// When writing kernels with these vectors, it is expected that floating-
// point operations will be carried out in a loop over Vec256<T>::float_num_vecs
// iterations.

namespace at {
namespace vec256 {
namespace {

const __vshi mask_unsigned = vec_splats((short int)0xFF);
template <>
struct Vec256<c10::quint8> {
 private:
  union {
    struct {
      __vchar _vec0;
      __vchar _vec1;
    };
    struct {
      __vcharb _vecb0;
      __vcharb _vecb1;
    };

  } __attribute__((__may_alias__));

 public:
  Vec256() {}
  static constexpr int size() {
    return 32;
  }

  static constexpr size_t float_num_vecs() {
    return 4;
  }
  static constexpr int int_num_vecs() {
    return 4;
  }
  using float_vec_return_type = std::array<Vec256<float>, 4>;
  using int_vec_return_type = std::array<Vec256<c10::qint32>, 4>;
  using value_type = typename c10::quint8::underlying;
  using vec_internal_type = __vchar;
  using vec_internal_mask_type = __vcharb;
  // Broadcast constructor
  __inline_attrs Vec256(const c10::quint8& val)
      : _vec0(vec_splats(val.val_)), _vec1(vec_splats(val.val_)) {}

  __inline_attrs Vec256(const Vec256<c10::quint8>& other)
      : _vec0{other._vec0}, _vec1(other._vec1) {}

  __inline_attrs Vec256(__vchar v) : _vec0{v}, _vec1{v} {}
  __inline_attrs Vec256(__vcharb vmask) : _vecb0{vmask}, _vecb1{vmask} {}
  __inline_attrs Vec256(__vchar v1, __vchar v2) : _vec0{v1}, _vec1{v2} {}
  __inline_attrs Vec256(__vcharb v1, __vcharb v2) : _vecb0{v1}, _vecb1{v2} {}

  inline __inline_attrs const vec_internal_type& vec0() const {
    return _vec0;
  }
  inline __inline_attrs const vec_internal_type& vec1() const {
    return _vec1;
  }

  static __inline_attrs Vec256<c10::quint8> loadu(
      const void* ptr,
      int count = size()) {
    if (count == size()) {
      return {
          vec_vsx_ld(offset0, reinterpret_cast<const value_type*>(ptr)),
          vec_vsx_ld(offset16, reinterpret_cast<const value_type*>(ptr))};
    }
    __at_align32__ value_type tmp_values[size()];
    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));
    return {vec_vsx_ld(offset0, tmp_values), vec_vsx_ld(offset16, tmp_values)};
  }
  void __inline_attrs store(void* ptr, int count = size()) const {
    if (count == size()) {
      vec_vsx_st(_vec0, offset0, reinterpret_cast<value_type*>(ptr));
      vec_vsx_st(_vec1, offset16, reinterpret_cast<value_type*>(ptr));
    } else if (count > 0) {
      __at_align32__ value_type tmp_values[size()];
      vec_vsx_st(_vec0, offset0, tmp_values);
      vec_vsx_st(_vec1, offset16, tmp_values);
      std::memcpy(
          ptr, tmp_values, std::min(count, size()) * sizeof(value_type));
    }
  }

 public:
  float_vec_return_type __inline_attrs dequantize(
      Vec256<float> scale,
      Vec256<float> zero_point,
      Vec256<float> scale_zp_premul) const {
    // unpacking unsigned as signed
    __vshi vecshi0 = vec_unpackh((__vchari)_vec0);
    __vshi vecshi1 = vec_unpackl((__vchari)_vec0);

    __vshi vecshi2 = vec_unpackh((__vchari)_vec1);
    __vshi vecshi3 = vec_unpackl((__vchari)_vec1);

    // signed ->  unsigned
    vecshi0 = vec_and(vecshi0, mask_unsigned);
    vecshi1 = vec_and(vecshi1, mask_unsigned);

    vecshi2 = vec_and(vecshi2, mask_unsigned);
    vecshi3 = vec_and(vecshi3, mask_unsigned);

    __vi veci0 = vec_unpackh(vecshi0);
    __vi veci1 = vec_unpackl(vecshi0);

    __vi veci2 = vec_unpackh(vecshi1);
    __vi veci3 = vec_unpackl(vecshi1);

    __vi veci4 = vec_unpackh(vecshi2);
    __vi veci5 = vec_unpackl(vecshi2);

    __vi veci6 = vec_unpackh(vecshi3);
    __vi veci7 = vec_unpackl(vecshi3);

    __vf vecf0_0 = vec_float(veci0);
    __vf vecf1_0 = vec_float(veci1);

    __vf vecf0_1 = vec_float(veci2);
    __vf vecf1_1 = vec_float(veci3);

    __vf vecf0_2 = vec_float(veci4);
    __vf vecf1_2 = vec_float(veci5);

    __vf vecf0_3 = vec_float(veci6);
    __vf vecf1_3 = vec_float(veci7);
    __vf scale_vec0 = scale.vec0();
    __vf scale_vec1 = scale.vec1();
    __vf scale_zp_premul0 = scale_zp_premul.vec0();
    __vf scale_zp_premul1 = scale_zp_premul.vec1();
    return {
        Vec256<float>{
            vec_madd(scale_vec0, vecf0_0, scale_zp_premul0),
            vec_madd(scale_vec1, vecf1_0, scale_zp_premul1)},
        Vec256<float>{
            vec_madd(scale_vec0, vecf0_1, scale_zp_premul0),
            vec_madd(scale_vec1, vecf1_1, scale_zp_premul1)},
        Vec256<float>{
            vec_madd(scale_vec0, vecf0_2, scale_zp_premul0),
            vec_madd(scale_vec1, vecf1_2, scale_zp_premul1)},
        Vec256<float>{
            vec_madd(scale_vec0, vecf0_3, scale_zp_premul0),
            vec_madd(scale_vec1, vecf1_3, scale_zp_premul1)}};
  }

  static Vec256<c10::quint8> quantize(
      const float_vec_return_type& rhs,
      float scale,
      int32_t zero_point,
      float inverse_scale) {
    // constexpr int32_t min_val = std::numeric_limits<value_type>::min();
    // constexpr int32_t max_val = std::numeric_limits<value_type>::max();

    __vf vec_inverse = vec_splats(inverse_scale);
    __vf vec_zero_point = vec_splats((float)zero_point);
    // __vui vmin = vec_splats(min_val);
    // __vui vmax = vec_splats(max_val);
    Vec256<float> vf0 = rhs[0];
    Vec256<float> vf1 = rhs[1];
    Vec256<float> vf2 = rhs[2];
    Vec256<float> vf3 = rhs[3];
    __vf vecf0 = vf0.vec0();
    __vf vecf1 = vf0.vec1();
    __vf vecf2 = vf1.vec0();
    __vf vecf3 = vf1.vec1();

    __vf vecf4 = vf2.vec0();
    __vf vecf5 = vf2.vec1();
    __vf vecf6 = vf3.vec0();
    __vf vecf7 = vf3.vec1();

    vecf0 = vec_mul(vecf0, vec_inverse);
    vecf1 = vec_mul(vecf1, vec_inverse);
    vecf2 = vec_mul(vecf2, vec_inverse);
    vecf3 = vec_mul(vecf3, vec_inverse);

    vecf4 = vec_mul(vecf4, vec_inverse);
    vecf5 = vec_mul(vecf5, vec_inverse);
    vecf6 = vec_mul(vecf6, vec_inverse);
    vecf7 = vec_mul(vecf7, vec_inverse);

    vecf0 = vec_add(vec_rint(vecf0), vec_zero_point);
    vecf1 = vec_add(vec_rint(vecf1), vec_zero_point);
    vecf2 = vec_add(vec_rint(vecf2), vec_zero_point);
    vecf3 = vec_add(vec_rint(vecf3), vec_zero_point);

    vecf4 = vec_add(vec_rint(vecf4), vec_zero_point);
    vecf5 = vec_add(vec_rint(vecf5), vec_zero_point);
    vecf6 = vec_add(vec_rint(vecf6), vec_zero_point);
    vecf7 = vec_add(vec_rint(vecf7), vec_zero_point);

    __vi veci0 = vec_signed(vecf0);
    __vi veci1 = vec_signed(vecf1);
    __vi veci2 = vec_signed(vecf2);
    __vi veci3 = vec_signed(vecf3);

    __vi veci4 = vec_signed(vecf4);
    __vi veci5 = vec_signed(vecf5);
    __vi veci6 = vec_signed(vecf6);
    __vi veci7 = vec_signed(vecf7);

    __vshi vecshi0 = vec_packs(veci0, veci1);
    __vshi vecshi1 = vec_packs(veci2, veci3);
    __vshi vecshi2 = vec_packs(veci4, veci5);
    __vshi vecshi3 = vec_packs(veci6, veci7);

    __vchar vec0 = vec_packsu(vecshi0, vecshi1);
    __vchar vec1 = vec_packsu(vecshi2, vecshi3);

    return {vec0, vec1};
  }

  Vec256<c10::quint8> __inline_attrs relu(Vec256<c10::quint8> zero_point) const {
    return {vec_max(_vec0, zero_point._vec0), vec_max(_vec1, zero_point._vec1)};
  }

  Vec256<c10::quint8> __inline_attrs
  relu6(Vec256<c10::quint8> zero_point, Vec256<c10::quint8> q_six) const {
    __vchar max0 = vec_max(_vec0, zero_point._vec0);
    __vchar max1 = vec_max(_vec1, zero_point._vec1);
    return {vec_min(max0, q_six._vec0), vec_min(max1, q_six._vec1)};
  }

  int_vec_return_type widening_subtract(Vec256<c10::quint8> b) const {
    __vshi vecshi0 = vec_unpackh((__vchari)_vec0);
    __vshi vecBshi0 = vec_unpackh((__vchari)b._vec0);
    __vshi vecshi1 = vec_unpackl((__vchari)_vec0);
    __vshi vecBshi1 = vec_unpackl((__vchari)b._vec0);

    __vshi vecshi2 = vec_unpackh((__vchari)_vec1);
    __vshi vecBshi2 = vec_unpackh((__vchari)b._vec1);
    __vshi vecshi3 = vec_unpackl((__vchari)_vec1);
    __vshi vecBshi3 = vec_unpackl((__vchari)b._vec1);

    vecshi0 = vec_and(vecshi0, mask_unsigned);
    vecBshi0 = vec_and(vecBshi0, mask_unsigned);
    vecshi1 = vec_and(vecshi1, mask_unsigned);
    vecBshi1 = vec_and(vecBshi1, mask_unsigned);

    vecshi2 = vec_and(vecshi2, mask_unsigned);
    vecBshi2 = vec_and(vecBshi2, mask_unsigned);
    vecshi3 = vec_and(vecshi3, mask_unsigned);
    vecBshi3 = vec_and(vecBshi3, mask_unsigned);

    __vi veci0 = vec_unpackh(vecshi0);
    __vi vecBi0 = vec_unpackh(vecBshi0);
    __vi veci1 = vec_unpackl(vecshi0);
    __vi vecBi1 = vec_unpackl(vecBshi0);

    __vi veci2 = vec_unpackh(vecshi1);
    __vi vecBi2 = vec_unpackh(vecBshi1);
    __vi veci3 = vec_unpackl(vecshi1);
    __vi vecBi3 = vec_unpackl(vecBshi1);

    __vi veci4 = vec_unpackh(vecshi2);
    __vi vecBi4 = vec_unpackh(vecBshi2);
    __vi veci5 = vec_unpackl(vecshi2);
    __vi vecBi5 = vec_unpackl(vecBshi2);

    __vi veci6 = vec_unpackh(vecshi3);
    __vi vecBi6 = vec_unpackh(vecBshi3);
    __vi veci7 = vec_unpackl(vecshi3);
    __vi vecBi7 = vec_unpackl(vecBshi3);

    return {
        Vec256<c10::qint32>(veci0 - vecBi0, veci1 - vecBi1),
        Vec256<c10::qint32>(veci2 - vecBi2, veci3 - vecBi3),
        Vec256<c10::qint32>(veci4 - vecBi4, veci5 - vecBi5),
        Vec256<c10::qint32>(veci6 - vecBi6, veci7 - vecBi7)};
  }

  static Vec256<c10::quint8> requantize_from_int(
      const int_vec_return_type& inp,
      float multiplier,
      int32_t zero_point) {
    __vf vec_multiplier = vec_splats(multiplier);
    __vi vec_zero_point = vec_splats(zero_point);

    Vec256<c10::qint32> vi0 = inp[0];
    Vec256<c10::qint32> vi1 = inp[1];
    Vec256<c10::qint32> vi2 = inp[2];
    Vec256<c10::qint32> vi3 = inp[3];

    __vf vecf0 = vec_float(vi0.vec0());
    __vf vecf1 = vec_float(vi0.vec1());
    __vf vecf2 = vec_float(vi1.vec0());
    __vf vecf3 = vec_float(vi1.vec1());

    __vf vecf4 = vec_float(vi2.vec0());
    __vf vecf5 = vec_float(vi2.vec1());
    __vf vecf6 = vec_float(vi3.vec0());
    __vf vecf7 = vec_float(vi3.vec1());

    vecf0 = vec_mul(vecf0, vec_multiplier);
    vecf1 = vec_mul(vecf1, vec_multiplier);
    vecf2 = vec_mul(vecf2, vec_multiplier);
    vecf3 = vec_mul(vecf3, vec_multiplier);

    vecf4 = vec_mul(vecf4, vec_multiplier);
    vecf5 = vec_mul(vecf5, vec_multiplier);
    vecf6 = vec_mul(vecf6, vec_multiplier);
    vecf7 = vec_mul(vecf7, vec_multiplier);

    vecf0 = vec_rint(vecf0);
    vecf1 = vec_rint(vecf1);
    vecf2 = vec_rint(vecf2);
    vecf3 = vec_rint(vecf3);

    vecf4 = vec_rint(vecf4);
    vecf5 = vec_rint(vecf5);
    vecf6 = vec_rint(vecf6);
    vecf7 = vec_rint(vecf7);

    __vi veci0 = vec_signed(vecf0);
    __vi veci1 = vec_signed(vecf1);
    __vi veci2 = vec_signed(vecf2);
    __vi veci3 = vec_signed(vecf3);

    __vi veci4 = vec_signed(vecf4);
    __vi veci5 = vec_signed(vecf5);
    __vi veci6 = vec_signed(vecf6);
    __vi veci7 = vec_signed(vecf7); 
    
    veci0 = vec_add(veci0, vec_zero_point);
    veci1 = vec_add(veci1, vec_zero_point);
    veci2 = vec_add(veci2, vec_zero_point);
    veci3 = vec_add(veci3, vec_zero_point);

    veci4 = vec_add(veci4, vec_zero_point);
    veci5 = vec_add(veci5, vec_zero_point);
    veci6 = vec_add(veci6, vec_zero_point);
    veci7 = vec_add(veci7, vec_zero_point);
 
    __vshi vecshi0 = vec_packs(veci0, veci1);
    __vshi vecshi1 = vec_packs(veci2, veci3);
    __vshi vecshi2 = vec_packs(veci4, veci5);
    __vshi vecshi3 = vec_packs(veci6, veci7);  

    __vchar vec0 = vec_packsu(vecshi0, vecshi1);
    __vchar vec1 = vec_packsu(vecshi2, vecshi3);

    return {vec0, vec1};
  }

  void dump() const {
    value_type vals[size()];
    store((void*)vals);
    for (int i = 0; i < size(); ++i) {
      std::cout << (int)(vals[i]) << " ";
    }
    std::cout << std::endl;
  }

  DEFINE_MEMBER_OP(operator==, c10::quint8, vec_cmpeq)
  DEFINE_MEMBER_OP(operator!=, c10::quint8, vec_cmpne)
  DEFINE_MEMBER_OP(operator<, c10::quint8, vec_cmplt)
  DEFINE_MEMBER_OP(operator<=, c10::quint8, vec_cmple)
  DEFINE_MEMBER_OP(operator>, c10::quint8, vec_cmpgt)
  DEFINE_MEMBER_OP(operator>=, c10::quint8, vec_cmpge)
  DEFINE_MEMBER_OP(operator+, c10::quint8, vec_add)
  DEFINE_MEMBER_OP(operator-, c10::quint8, vec_sub)
  DEFINE_MEMBER_OP(operator*, c10::quint8, vec_mul)
  DEFINE_MEMBER_EMULATE_BINARY_OP(operator/, c10::quint8, /)
  DEFINE_MEMBER_OP(maximum, c10::quint8, vec_max)
  DEFINE_MEMBER_OP(minimum, c10::quint8, vec_min)
  DEFINE_MEMBER_OP(operator&, c10::quint8, vec_and)
  DEFINE_MEMBER_OP(operator|, c10::quint8, vec_or)
  DEFINE_MEMBER_OP(operator^, c10::quint8, vec_xor)
};

template <>
Vec256<c10::quint8> inline maximum(
    const Vec256<c10::quint8>& a,
    const Vec256<c10::quint8>& b) {
  return a.maximum(b);
}

template <>
Vec256<c10::quint8> inline minimum(
    const Vec256<c10::quint8>& a,
    const Vec256<c10::quint8>& b) {
  return a.minimum(b);
}

} // namespace
} // namespace vec256
} // namespace at