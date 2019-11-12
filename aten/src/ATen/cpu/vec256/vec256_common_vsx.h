#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>
#include <ATen/cpu/vec256/vec256_double_vsx.h>
#include <ATen/cpu/vec256/vec256_float_vsx.h>
#include <ATen/cpu/vec256/vec256_int_vsx.h>
#include <ATen/cpu/vec256/vec256_qint_vsx.h>

namespace at {
namespace vec256 {
// See Note [Acceptable use of anonymous namespace in header]

namespace {
#if defined(__VSX__)

#define DEFINE_BINARY_OP(op, operand_type, func)                      \
  template <>                                                         \
  Vec256<operand_type> inline __inline_attrs op(                      \
      const Vec256<operand_type>& a, const Vec256<operand_type>& b) { \
    return Vec256<operand_type>{func(a.vec0(), b.vec0()),             \
                                func(a.vec1(), b.vec1())};            \
  }

#define DEFINE_CLAMP_FUNCS(operand_type)                                \
  template <>                                                           \
  Vec256<operand_type> inline __inline_attrs clamp(                     \
      const Vec256<operand_type>& a,                                    \
      const Vec256<operand_type>& min,                                  \
      const Vec256<operand_type>& max) {                                \
    return Vec256<operand_type>{                                        \
        vec_min(max.vec0(), vec_max(a.vec0(), min.vec0())),             \
        vec_min(max.vec1(), vec_max(a.vec1(), min.vec1()))};            \
  }                                                                     \
  template <>                                                           \
  Vec256<operand_type> inline __inline_attrs clamp_min(                 \
      const Vec256<operand_type>& a, const Vec256<operand_type>& min) { \
    return Vec256<operand_type>{vec_max(a.vec0(), min.vec0()),          \
                                vec_max(a.vec1(), min.vec1())};         \
  }                                                                     \
  template <>                                                           \
  Vec256<operand_type> inline __inline_attrs clamp_max(                 \
      const Vec256<operand_type>& a, const Vec256<operand_type>& max) { \
    return Vec256<operand_type>{vec_min(a.vec0(), max.vec0()),          \
                                vec_min(a.vec1(), max.vec1())};         \
  }

#define DEFINE_EMULATE_BINARY_FUNCS(fname, operand_type, binary_op)       \
  template <>                                                             \
  Vec256<operand_type> inline __inline_attrs fname(                       \
      const Vec256<operand_type>& a, const Vec256<operand_type>& b) {     \
    __at_align32__ Vec256<operand_type>::value_type values_a[Vec256<operand_type>::size()];   \
    __at_align32__ Vec256<operand_type>::value_type values_b[Vec256<operand_type>::size()];   \
    a.store(values_a);                                                    \
    b.store(values_b);                                                    \
    \
	for (int i = 0; i != Vec256<operand_type>::size(); i++) {             \
      values_a[i]  = values_a[i] binary_op values_b[i];                   \
    }                                                                     \
    return Vec256<operand_type>::loadu(values_a);                         \
  }

DEFINE_BINARY_OP(operator+, double, vec_add)
DEFINE_BINARY_OP(operator-, double, vec_sub)
DEFINE_BINARY_OP(operator*, double, vec_mul)
DEFINE_BINARY_OP(operator/, double, vec_div)
DEFINE_BINARY_OP(operator&, double, vec_and)
DEFINE_BINARY_OP(operator|, double, vec_or)
DEFINE_BINARY_OP(operator^, double, vec_xor)
DEFINE_BINARY_OP(maximum, double, vec_max)
DEFINE_BINARY_OP(minimum, double, vec_min)

DEFINE_BINARY_OP(operator+, float, vec_add)
DEFINE_BINARY_OP(operator-, float, vec_sub)
DEFINE_BINARY_OP(operator*, float, vec_mul)
DEFINE_BINARY_OP(operator/, float, vec_div)
DEFINE_BINARY_OP(operator&, float, vec_and)
DEFINE_BINARY_OP(operator|, float, vec_or)
DEFINE_BINARY_OP(operator^, float, vec_xor)
DEFINE_BINARY_OP(maximum, float, vec_max)
DEFINE_BINARY_OP(minimum, float, vec_min)

DEFINE_BINARY_OP(operator+, int64_t, vec_add)
DEFINE_BINARY_OP(operator-, int64_t, vec_sub)
DEFINE_BINARY_OP(operator*, int64_t, vec_mul)
DEFINE_BINARY_OP(operator&, int64_t, vec_and)
DEFINE_BINARY_OP(operator|, int64_t, vec_or)
DEFINE_BINARY_OP(operator^, int64_t, vec_xor)
DEFINE_BINARY_OP(maximum, int64_t, vec_max)
DEFINE_BINARY_OP(minimum, int64_t, vec_min)

DEFINE_BINARY_OP(operator+, int32_t, vec_add)
DEFINE_BINARY_OP(operator-, int32_t, vec_sub)
DEFINE_BINARY_OP(operator*, int32_t, vec_mul)
DEFINE_BINARY_OP(operator&, int32_t, vec_and)
DEFINE_BINARY_OP(operator|, int32_t, vec_or)
DEFINE_BINARY_OP(operator^, int32_t, vec_xor)
DEFINE_BINARY_OP(maximum, int32_t, vec_max)
DEFINE_BINARY_OP(minimum, int32_t, vec_min)

DEFINE_BINARY_OP(operator+, int16_t, vec_add)
DEFINE_BINARY_OP(operator-, int16_t, vec_sub)
DEFINE_BINARY_OP(operator*, int16_t, vec_mul)
DEFINE_BINARY_OP(operator&, int16_t, vec_and)
DEFINE_BINARY_OP(operator|, int16_t, vec_or)
DEFINE_BINARY_OP(operator^, int16_t, vec_xor)
DEFINE_BINARY_OP(maximum, int16_t, vec_max)
DEFINE_BINARY_OP(minimum, int16_t, vec_min)

DEFINE_BINARY_OP(operator+, c10::qint32, vec_add)
DEFINE_BINARY_OP(operator-, c10::qint32, vec_sub)
DEFINE_BINARY_OP(operator*, c10::qint32, vec_mul)
DEFINE_BINARY_OP(operator&, c10::qint32, vec_and)
DEFINE_BINARY_OP(operator|, c10::qint32, vec_or)
DEFINE_BINARY_OP(operator^, c10::qint32, vec_xor)
DEFINE_BINARY_OP(maximum, c10::qint32, vec_max)
DEFINE_BINARY_OP(minimum, c10::qint32, vec_min)

DEFINE_BINARY_OP(operator+, c10::qint8, vec_add)
DEFINE_BINARY_OP(operator-, c10::qint8, vec_sub)
DEFINE_BINARY_OP(operator*, c10::qint8, vec_mul)
DEFINE_BINARY_OP(operator&, c10::qint8, vec_and)
DEFINE_BINARY_OP(operator|, c10::qint8, vec_or)
DEFINE_BINARY_OP(operator^, c10::qint8, vec_xor)
DEFINE_BINARY_OP(maximum, c10::qint8, vec_max)
DEFINE_BINARY_OP(minimum, c10::qint8, vec_min)

DEFINE_BINARY_OP(operator+, c10::quint8, vec_add)
DEFINE_BINARY_OP(operator-, c10::quint8, vec_sub)
DEFINE_BINARY_OP(operator*, c10::quint8, vec_mul)
DEFINE_BINARY_OP(operator&, c10::quint8, vec_and)
DEFINE_BINARY_OP(operator|, c10::quint8, vec_or)
DEFINE_BINARY_OP(operator^, c10::quint8, vec_xor)
DEFINE_BINARY_OP(maximum, c10::quint8, vec_max)
DEFINE_BINARY_OP(minimum, c10::quint8, vec_min)

DEFINE_CLAMP_FUNCS(c10::quint8)
DEFINE_CLAMP_FUNCS(c10::qint8)
DEFINE_CLAMP_FUNCS(c10::qint32)
DEFINE_CLAMP_FUNCS(int16_t)
DEFINE_CLAMP_FUNCS(int32_t)
DEFINE_CLAMP_FUNCS(int64_t)
DEFINE_CLAMP_FUNCS(float)
DEFINE_CLAMP_FUNCS(double)

DEFINE_EMULATE_BINARY_FUNCS(operator/, int32_t, /)
DEFINE_EMULATE_BINARY_FUNCS(operator/, int16_t, /)
DEFINE_EMULATE_BINARY_FUNCS(operator/, c10::qint32, /)
DEFINE_EMULATE_BINARY_FUNCS(operator/, c10::qint8, /)
DEFINE_EMULATE_BINARY_FUNCS(operator/, c10::quint8, /)
template <>
Vec256<int64_t> inline __inline_attrs operator/(
    const Vec256<int64_t>& a,
    const Vec256<int64_t>& b) {
  __vlli a_vec0 = a.vec0();
  __vlli a_vec1 = a.vec1();
  __vlli b_vec0 = b.vec0();
  __vlli b_vec1 = b.vec1();
  auto r1 = a_vec0[0] / b_vec0[0];
  auto r2 = a_vec0[1] / b_vec0[1];
  auto r3 = a_vec1[0] / b_vec1[0];
  auto r4 = a_vec1[1] / b_vec1[1];
  return Vec256<int64_t>{r1, r2, r3, r4};
}

template <>
Vec256<double> inline __inline_attrs fmadd(
    const Vec256<double>& a,
    const Vec256<double>& b,
    const Vec256<double>& c) {
  return Vec256<double>{vec_madd(a.vec0(), b.vec0(), c.vec0()),
                        vec_madd(a.vec1(), b.vec1(), c.vec1())};
}
template <>
Vec256<float> inline __inline_attrs fmadd(
    const Vec256<float>& a,
    const Vec256<float>& b,
    const Vec256<float>& c) {
  return Vec256<float>{vec_madd(a.vec0(), b.vec0(), c.vec0()),
                       vec_madd(a.vec1(), b.vec1(), c.vec1())};
}
template <>
Vec256<int64_t> inline __inline_attrs fmadd(
    const Vec256<int64_t>& a,
    const Vec256<int64_t>& b,
    const Vec256<int64_t>& c) {
  return Vec256<int64_t>{a.vec0() * b.vec0() + c.vec0(),
                         a.vec1() * b.vec1() + c.vec1()};
}
template <>
Vec256<int32_t> inline __inline_attrs fmadd(
    const Vec256<int32_t>& a,
    const Vec256<int32_t>& b,
    const Vec256<int32_t>& c) {
  return Vec256<int32_t>{a.vec0() * b.vec0() + c.vec0(),
                         a.vec1() * b.vec1() + c.vec1()};
}
template <>
Vec256<int16_t> inline __inline_attrs fmadd(
    const Vec256<int16_t>& a,
    const Vec256<int16_t>& b,
    const Vec256<int16_t>& c) {
  return Vec256<int16_t>{a.vec0() * b.vec0() + c.vec0(),
                         a.vec1() * b.vec1() + c.vec1()};
}

#define DEFINE_REINTERPRET_CAST_FUNCS(                                 \
    first_type, cast_type, cast_inner_vector_type)                     \
  template <>                                                          \
  inline __inline_attrs Vec256<cast_type> cast<cast_type, first_type>( \
      const Vec256<first_type>& src) {                                 \
    return Vec256<cast_type>{(cast_inner_vector_type)src.vec0(),       \
                             (cast_inner_vector_type)src.vec1()};      \
  }

#define DEFINE_REINTERPRET_CAST_TO_ALL_FUNCS(first_type)       \
  DEFINE_REINTERPRET_CAST_FUNCS(first_type, double, __vd)      \
  DEFINE_REINTERPRET_CAST_FUNCS(first_type, float, __vf) \
			DEFINE_REINTERPRET_CAST_FUNCS(         \
      first_type, int64_t, __vlli)                             \
      DEFINE_REINTERPRET_CAST_FUNCS(first_type, int32_t, __vi) \
          DEFINE_REINTERPRET_CAST_FUNCS(first_type, int16_t, __vshi)

DEFINE_REINTERPRET_CAST_TO_ALL_FUNCS(float)
DEFINE_REINTERPRET_CAST_TO_ALL_FUNCS(double)
DEFINE_REINTERPRET_CAST_TO_ALL_FUNCS(int64_t)
DEFINE_REINTERPRET_CAST_TO_ALL_FUNCS(int32_t)
DEFINE_REINTERPRET_CAST_TO_ALL_FUNCS(int16_t)

template <>
Vec256<int64_t> inline __inline_attrs convert_to_int_of_same_size<double>(
    const Vec256<double>& src) {
  return Vec256<int64_t>{vec_signed(src.vec0()), vec_signed(src.vec1())};
}

template <>
Vec256<int32_t> inline __inline_attrs convert_to_int_of_same_size<float>(
    const Vec256<float>& src) {
  return Vec256<int32_t>{vec_signed(src.vec0()), vec_signed(src.vec1())};
}

template <>
inline void convert(const int32_t* src, float* dst, int64_t n) {
  int64_t i;
  // int32_t and float have same size
#pragma unroll
  for (i = 0; i <= (n - Vec256<int32_t>::size());
       i += Vec256<int32_t>::size()) {
    const int32_t* src_a = src + i;
    float* dst_a = dst + i;
    __vi input_vec0 = vec_vsx_ld(offset0, reinterpret_cast<const __vi*>(src_a));
    __vi input_vec1 =
        vec_vsx_ld(offset16, reinterpret_cast<const __vi*>(src_a));
    vec_vsx_st(vec_float(input_vec0), offset0, reinterpret_cast<float*>(dst_a));
    vec_vsx_st(
        vec_float(input_vec1), offset16, reinterpret_cast<float*>(dst_a));
  }
#ifndef _MSC_VER
#pragma unroll
#endif
  for (; i < n; i++) {
    dst[i] = static_cast<float>(src[i]);
  }
}

template <>
inline void convert(const int64_t* src, double* dst, int64_t n) {
  int64_t i;
  // int64_t and double have same size

#pragma unroll
  for (i = 0; i <= (n - Vec256<int64_t>::size());
       i += Vec256<int64_t>::size()) {
    const int64_t* src_a = src + i;
    double* dst_a = dst + i;
    __vlli input_vec0 =
        vec_vsx_ld(offset0, reinterpret_cast<const __vlli*>(src_a));
    __vlli input_vec1 =
        vec_vsx_ld(offset16, reinterpret_cast<const __vlli*>(src_a));
    vec_vsx_st(
        vec_double(input_vec0), offset0, reinterpret_cast<double*>(dst_a));
    vec_vsx_st(
        vec_double(input_vec1), offset16, reinterpret_cast<double*>(dst_a));
  }
#pragma unroll
  for (; i < n; i++) {
    dst[i] = static_cast<double>(src[i]);
  }
}

#undef DEFINE_BINARY_OP
#undef DEFINE_CLAMP_FUNCS
#undef DEFINE_EMULATE_BINARY_FUNCS
#undef DEFINE_REINTERPRET_CAST_TO_ALL_FUNCS
#undef DEFINE_REINTERPRET_CAST_FUNCS

template <>
std::pair<Vec256<double>, Vec256<double>> inline interleave2<double>(
    const Vec256<double>& a,
    const Vec256<double>& b) {
  // inputs:
  //   a      = {a0, a1, a2, a3}
  //   b      = {b0, b1, b2, b3}

  __vd ab00 = vec_xxpermdi(a.vec0(), b.vec0(), 0);
  __vd ab11 = vec_xxpermdi(a.vec0(), b.vec0(), 3);
  __vd ab2_00 = vec_xxpermdi(a.vec1(), b.vec1(), 0);
  __vd ab2_11 = vec_xxpermdi(a.vec1(), b.vec1(), 3);
  //   return {a0, b0, a1, b1}
  //          {a2, b2, a3, b3}
  return std::make_pair(
      Vec256<double>{ab00, ab11}, Vec256<double>{ab2_00, ab2_11});
}

template <>
std::pair<Vec256<double>, Vec256<double>> inline deinterleave2<double>(
    const Vec256<double>& a,
    const Vec256<double>& b) {
  // inputs:
  //   a = {a0, b0, a1, b1}
  //   b = {a2, b2, a3, b3}
  __vd aa01 = vec_xxpermdi(a.vec0(), a.vec1(), 0);
  __vd aa23 = vec_xxpermdi(b.vec0(), b.vec1(), 0);

  __vd bb_01 = vec_xxpermdi(a.vec0(), a.vec1(), 3);
  __vd bb_23 = vec_xxpermdi(b.vec0(), b.vec1(), 3);

  // swap lanes:
  //   return {a0, a1, a2, a3}
  //          {b0, b1, b2, b3}
  return std::make_pair(
      Vec256<double>{aa01, aa23}, Vec256<double>{bb_01, bb_23});
}

template <>
std::pair<Vec256<float>, Vec256<float>> inline interleave2<float>(
    const Vec256<float>& a,
    const Vec256<float>& b) {
  // inputs:
  //   a = {a0, a1, a2, a3,, a4, a5, a6, a7}
  //   b = {b0, b1, b2, b3,, b4, b5, b6, b7}

  __vf ab0011 = vec_mergeh(a.vec0(), b.vec0());
  __vf ab2233 = vec_mergel(a.vec0(), b.vec0());

  __vf ab2_0011 = vec_mergeh(a.vec1(), b.vec1());
  __vf ab2_2233 = vec_mergel(a.vec1(), b.vec1());
  // group cols crossing lanes:
  //   return {a0, b0, a1, b1,, a2, b2, a3, b3}
  //          {a4, b4, a5, b5,, a6, b6, a7, b7}

  return std::make_pair(
      Vec256<float>{ab0011, ab2233}, Vec256<float>{ab2_0011, ab2_2233});
}

template <>
std::pair<Vec256<float>, Vec256<float>> inline deinterleave2<float>(
    const Vec256<float>& a,
    const Vec256<float>& b) {
  // inputs:
  //   a = {a0, b0, a1, b1,, a2, b2, a3, b3}
  //   b = {a4, b4, a5, b5,, a6, b6, a7, b7}

  // {a0,a2,b0,b2} {a1,a3,b1,b3}
  __vf a0a2b0b2 = vec_mergeh(a.vec0(), a.vec1());
  __vf a1a3b1b3 = vec_mergel(a.vec0(), a.vec1());

  __vf aa0123 = vec_mergeh(a0a2b0b2, a1a3b1b3);
  __vf bb0123 = vec_mergel(a0a2b0b2, a1a3b1b3);

  __vf a0a2b0b2_2 = vec_mergeh(b.vec0(), b.vec1());
  __vf a1a3b1b3_2 = vec_mergel(b.vec0(), b.vec1());

  __vf aa0123_2 = vec_mergeh(a0a2b0b2_2, a1a3b1b3_2);
  __vf bb0123_2 = vec_mergel(a0a2b0b2_2, a1a3b1b3_2);

  // it could be done with vec_perm ,too
  // swap lanes:
  //   return {a0, a1, a2, a3,, a4, a5, a6, a7}
  //          {b0, b1, b2, b3,, b4, b5, b6, b7}

  return std::make_pair(
      Vec256<float>{aa0123, aa0123_2}, Vec256<float>{bb0123, bb0123_2});
}

#endif
} // namespace
} // namespace vec256
} // namespace at