#pragma once
#include <ATen/cpu/vec256/intrinsics.h>

#include <cstdint>

#define DEFINE_MEMBER_OP(op, op_type, func)                              \
  Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>& other) \
      const {                                                            \
    return Vec256<op_type>{func(_vec0, other._vec0),                     \
                           func(_vec1, other._vec1)};                    \
  }

#define DEFINE_MEMBER_BITWISE_OP(op, op_type, func)                      \
  Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>& other) \
      const {                                                            \
    return Vec256<op_type>{func(_vecb0, other._vecb0),                   \
                           func(_vecb1, other._vecb1)};                  \
  }

#define DEFINE_MEMBER_TERNARY_OP(op, op_type, func)                   \
  Vec256<op_type> __inline_attrs op(const Vec256<op_type>& b,         \
                                    const Vec256<op_type>& c) const { \
    return Vec256<op_type>{func(_vec0, b._vec0, c._vec0),             \
                           func(_vec1, b._vec1, c._vec1)};            \
  }

#define DEFINE_MEMBER_EMULATE_BINARY_OP(op, op_type, binary_op)       \
  Vec256<op_type> __inline_attrs op(const Vec256<op_type>& b) const { \
    Vec256<op_type>::vec_internal_type ret_0;                         \
    Vec256<op_type>::vec_internal_type ret_1;                         \
    for (int i = 0; i < Vec256<op_type>::size() / 2; i++) {           \
      ret_0[i] = _vec0[i] binary_op b._vec0[i];                       \
      ret_1[i] = _vec1[i] binary_op b._vec1[i];                       \
    }                                                                 \
    return Vec256<op_type>{ret_0, ret_1};                             \
  }

template <typename T>
struct SplatType {
  using type = T;
};

#define DEFINE_MEMBER_OP_AND_ONE(op, op_type, func)                      \
  Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>& other) \
      const {                                                            \
    using vvtype = Vec256<op_type>::vec_internal_type;                   \
    const vvtype v_one = vec_splats(static_cast<op_type>(1.0));          \
    vvtype ret0 = (vvtype)func(_vec0, other._vec0);                      \
    vvtype ret1 = (vvtype)func(_vec1, other._vec1);                      \
    return Vec256<op_type>{vec_and(ret0, v_one), vec_and(ret1, v_one)};  \
  }

#define DEFINE_CLAMP_FUNCS(operand_type)                                \
  template <>                                                           \
  Vec256<operand_type> inline __inline_attrs clamp(                     \
      const Vec256<operand_type>& a, const Vec256<operand_type>& min,   \
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

#define DEFINE_REINTERPRET_CAST_FUNCS(first_type, cast_type,           \
                                      cast_inner_vector_type)          \
  template <>                                                          \
  inline __inline_attrs Vec256<cast_type> cast<cast_type, first_type>( \
      const Vec256<first_type>& src) {                                 \
    return Vec256<cast_type>{(cast_inner_vector_type)src.vec0(),       \
                             (cast_inner_vector_type)src.vec1()};      \
  }

#define DEFINE_REINTERPRET_CAST_TO_ALL_FUNCS(first_type)     \
  DEFINE_REINTERPRET_CAST_FUNCS(first_type, double, __vd)    \
  DEFINE_REINTERPRET_CAST_FUNCS(first_type, float, __vf)     \
  DEFINE_REINTERPRET_CAST_FUNCS(first_type, int64_t, __vlli) \
  DEFINE_REINTERPRET_CAST_FUNCS(first_type, int32_t, __vi)   \
  DEFINE_REINTERPRET_CAST_FUNCS(first_type, int16_t, __vshi)

// it can be used to emulate blend faster
constexpr int blendChoice(int mask, int half1 = 0xF, int half2 = 0xF0) {
  int none = 0;
  int both = half1 | half2;
  // clamp it between 0 and both
  mask = mask & both;
  // return  (a._vec0, a._vec1)
  if (mask == none) return 0;
  // return (b._vec0,b._vec1)
  else if (mask == both)
    return 1;
  // return  (b._vec0,a._vec1)
  else if (mask == half1)
    return 2;
  // return  (a._vec0,b._vec1)
  else if (mask == half2)
    return 3;
  // return  (*_vec0,a._vec1)
  else if (mask > 0 && mask < half1)
    return 4;
  // return  (*_vec0,b._vec1)
  else if ((mask & half2) == half2)
    return 5;
  // return (a._vec0,*_vec1)
  else if ((mask & half1) == 0 && mask > half1)
    return 6;
  // return (b._vec0,*_vec1)
  else if ((mask & half1) == half1 && mask > half1)
    return 7;
  // return (*_vec0,*_vec1)
  return 8;
}

// it can be used to emulate blend faster
constexpr int blendChoiceDbl(int mask) {
  // clamp it 0 and 0xF
  return blendChoice(mask, 0x3, 0xC);
}

constexpr __vib VsxMask1(int mask) {
  uint32_t g0 = (mask & 1) * 0xffffffff;
  uint32_t g1 = ((mask & 2) >> 1) * 0xffffffff;
  uint32_t g2 = ((mask & 4) >> 2) * 0xffffffff;
  uint32_t g3 = ((mask & 8) >> 3) * 0xffffffff;
  return (__vib){g0, g1, g2, g3};
}

constexpr __vib VsxMask2(int mask) {
  uint32_t mask2 = (mask & 0xFF) >> 4;
  return VsxMask1(mask2);
}

constexpr __vllb VsxDblMask1(int mask) {
  uint64_t g0 = (mask & 1) * 0xffffffffffffffff;
  uint64_t g1 = ((mask & 2) >> 1) * 0xffffffffffffffff;
  return (__vllb){g0, g1};
}

constexpr __vllb VsxDblMask2(int mask) {
  uint32_t mask2 = (mask & 0xF) >> 2;
  return VsxDblMask1(mask2);
}

constexpr int maskForComplex(int mask) {
  mask = mask & 0xF;
  int complex_mask = 0;
  if (mask & 1) complex_mask |= 3;
  if (mask & 2) complex_mask |= (3 << 2);
  if (mask & 4) complex_mask |= (3 << 4);
  if (mask & 8) complex_mask |= (3 << 6);
  return complex_mask;
}

constexpr int maskForComplexDbl(int mask) {
  mask = mask & 0x3;
  int complex_mask = 0;
  if (mask & 1) complex_mask |= 3;
  if (mask & 2) complex_mask |= (3 << 2);
  return complex_mask;
}

constexpr int blendChoiceComplex(int mask) {
  return blendChoice(maskForComplex(mask));
}

constexpr int blendChoiceComplexDbl(int mask) {
  return blendChoiceDbl(maskForComplexDbl(mask));
}

constexpr __vib VsxComplexMask1(int mask) {
  return VsxMask1(maskForComplex(mask));
}

constexpr __vib VsxComplexMask2(int mask) {
  uint32_t mask2 = (mask & 0xF) >> 2;
  return VsxMask1(maskForComplex(mask2));
}

constexpr __vllb VsxComplexDblMask1(int mask) { return VsxDblMask1(mask); }

constexpr __vllb VsxComplexDblMask2(int mask) {
  uint32_t mask2 = (mask & 0xF) >> 2;
  return VsxDblMask1(mask2);
}

// constants
namespace at {
namespace vec256 {
// See Note [Acceptable use of anonymous namespace in header]
namespace {
//#Constants
const __vchar mask_zero_bits = __vchar{128, 128, 128, 128, 128, 128, 128, 128,
                                128, 128, 128, 128, 96,  64,  32,  0};

const __vchar swap_mask =
    __vchar{4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11};

const __vi v0x7f = vec_splats(0x7f);
const __vi vi_0 = vec_splats((int)(0));
const __vi vi_1 = vec_splats((int)1);
const __vi vi_2 = vec_splats((int)2);
const __vi vi_4 = vec_splats((int)4);
const __vi vi_inv1 = vec_splats((int)~1);
const __vui vu_29 = vec_splats(29u);
const __vui vu_23 = vec_splats(23u);

const __vib inv_mant_mask = (__vib)vec_splats((unsigned int)~0xff800000);
const __vib sign_mask = (__vib)vec_splats((int)0x80000000);
const __vib real_mask = __vib{0xFFFFFFFF, 0x0, 0xFFFFFFFF, 0x0};
const __vib imag_mask = __vib{0x0, 0xFFFFFFFF, 0x0, 0xFFFFFFFF};
const __vib isign_mask = __vib{0x0, 0x80000000, 0x0, 0x80000000};
const __vib rsign_mask = __vib{0x80000000, 0x0, 0x80000000, 0x0};

const __vllb vd_imag_mask  = __vllb{0x0, 0xFFFFFFFFFFFFFFFF};
const __vllb vd_real_mask  = __vllb{0xFFFFFFFFFFFFFFFF, 0x0};
const __vllb vd_isign_mask = __vllb{0x0, 0x8000000000000000};
const __vllb vd_rsign_mask = __vllb{0x8000000000000000, 0x0};

const __vf zero = vec_splats(0.f);
const __vf half = vec_splats(0.5f);
const __vf one = vec_splats(1.f);
const __vf two = vec_splats(2.0f);
const __vf _4div_pi = vec_splats(1.27323954473516f);
const __vf v_inf = (__vf)vec_splats(0x7f800000u);
const __vf v_nan = (__vf)vec_splats(0x7fffffff);
const __vf log10e_inv = vec_splats(0.43429448190325176f);
const __vf log2e_inv = vec_splats(1.4426950408889634f);
const __vf log2eB_inv = vec_splats(1.442695036924675f);
const __vf cephes_SQRTHF = vec_splats(0.707106781186547524f);
const __vf coscof_p0 = vec_splats(2.443315711809948E-005f);
const __vf coscof_p1 = vec_splats(-1.388731625493765E-003f);
const __vf coscof_p2 = vec_splats(4.166664568298827E-002f);
const __vf exp_hi = vec_splats(104.f);
const __vf exp_lo = vec_splats(-104.f);
const __vf exp_p0 = vec_splats(0.000198527617612853646278381f);
const __vf exp_p1 = vec_splats((0.00139304355252534151077271f));
const __vf exp_p2 = vec_splats(0.00833336077630519866943359f);
const __vf exp_p3 = vec_splats(0.0416664853692054748535156f);
const __vf exp_p4 = vec_splats(0.166666671633720397949219f);
const __vf exp_p5 = vec_splats(0.5f);
const __vf log_p0 = vec_splats(7.0376836292E-2f);
const __vf log_p1 = vec_splats(-1.1514610310E-1f);
const __vf log_p2 = vec_splats(1.1676998740E-1f);
const __vf log_p3 = vec_splats(-1.2420140846E-1f);
const __vf log_p4 = vec_splats(+1.4249322787E-1f);
const __vf log_p5 = vec_splats(-1.6668057665E-1f);
const __vf log_p6 = vec_splats(+2.0000714765E-1f);
const __vf log_p7 = vec_splats(-2.4999993993E-1f);
const __vf log_p8 = vec_splats(+3.3333331174E-1f);
const __vf log_q1 = vec_splats(-2.12194440e-4f);
const __vf log_q2 = vec_splats(0.693359375f);
const __vf max_logf = vec_splats(88.02969187150841f);
const __vf max_numf = vec_splats(1.7014117331926442990585209174225846272e38f);
const __vf min_inf = (__vf)vec_splats(0xff800000u);
const __vf min_norm_pos = (__vf)vec_splats(0x0800000u);
const __vf minus_cephes_dp1 = vec_splats(-0.78515625f);
const __vf minus_cephes_dp2 = vec_splats(-2.4187564849853515625e-4f);
const __vf minus_cephes_dp3 = vec_splats(-3.77489497744594108e-8f);
const __vf negln2f_hi = vec_splats(-0.693145751953125f);
const __vf negln2f_lo = vec_splats(-1.428606765330187045e-06f);
const __vf p0 = vec_splats(2.03721912945E-4f);
const __vf p1 = vec_splats(8.33028376239E-3f);
const __vf p2 = vec_splats(1.66667160211E-1f);
const __vf sincof_p0 = vec_splats(-1.9515295891E-4f);
const __vf sincof_p1 = vec_splats(8.3321608736E-3f);
const __vf sincof_p2 = vec_splats(-1.6666654611E-1f);
const __vf tanh_0p625 = vec_splats(0.625f);
const __vf tanh_half_max = vec_splats(44.014845935754205f);
const __vf tanh_p0 = vec_splats(-5.70498872745E-3f);
const __vf tanh_p1 = vec_splats(2.06390887954E-2f);
const __vf tanh_p2 = vec_splats(-5.37397155531E-2f);
const __vf tanh_p3 = vec_splats(1.33314422036E-1f);
const __vf tanh_p4 = vec_splats(-3.33332819422E-1f);
const __vf vcheck = vec_splats((float)(1LL << 24));
const __vf imag_one = __vf{0.f, 1.f, 0.f, 1.f};
const __vf imag_half = __vf{0.f, 0.5f, 0.f, 0.5f};
const __vf sqrt2_2 = __vf{0.70710676908493042f, 0.70710676908493042,
                          0.70710676908493042, 0.70710676908493042};
const __vf pi_2 = __vf{M_PI / 2, 0.0, M_PI / 2, 0.0};

const __vd vd_one = vec_splats(1.0);
const __vd vd_zero = vec_splats(0.0);
const __vd vd_log10e_inv = vec_splats(0.43429448190325176);
const __vd vd_log2e_inv = vec_splats(1.4426950408889634);
const __vd vd_imag_one = __vd{0.0, 1.0};
const __vd vd_imag_half = __vd{0.0, 0.5};
const __vd vd_sqrt2_2 = __vd{0.70710678118654757, 0.70710678118654757};
const __vd vd_pi_2 = __vd{M_PI / 2.0, 0.0};

}  // namespace
}  // namespace vec256
}  // namespace at