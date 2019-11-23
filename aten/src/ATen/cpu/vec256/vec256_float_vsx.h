#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>

namespace at {
namespace vec256 {
// See Note [Acceptable use of anonymous namespace in header]

namespace {

#if defined(__VSX__)

template <>
class Vec256<float> {
 private:
  union {
    struct {
      __vf _vec0;
      __vf _vec1;
    };
    struct {
      __vib _vecb0;
      __vib _vecb1;
    };

  } __attribute__((__may_alias__));

 public:
  using value_type = float;
  using vec_internal_type = __vf;

  static constexpr int size() {
    return 8;
  }
  Vec256() {}
  __inline_attrs Vec256(__vf v1, __vf v2) : _vec0{v1}, _vec1{v2} {}
  __inline_attrs Vec256(__vib v1, __vib v2) : _vecb0{v1}, _vecb1{v2} {}
  __inline_attrs Vec256(float scalar)
      : _vec0{vec_splats(scalar)}, _vec1{vec_splats(scalar)} {}
  __inline_attrs Vec256(
      float scalar1,
      float scalar2,
      float scalar3,
      float scalar4,
      float scalar5,
      float scalar6,
      float scalar7,
      float scalar8)
      : _vec0{__vf{scalar1, scalar2, scalar3, scalar4}},
        _vec1{__vf{scalar5, scalar6, scalar7, scalar8}} {}
  inline __inline_attrs const vec_internal_type& vec0() const {
    return _vec0;
  }
  inline __inline_attrs const vec_internal_type& vec1() const {
    return _vec1;
  }

  template <int64_t mask>
  static c10::guts::enable_if_t<mask == 0, Vec256<float>> __inline_attrs
  blend(const Vec256<float>& a, const Vec256<float>& b) {
    // here I am using intel style mask number
    // std::cout << "0" << std::endl;
    return a;
  }

  template <int64_t mask>
  static c10::guts::enable_if_t<(mask & 255) == 255, Vec256<float>>
      __inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
    // here I am using intel style mask number
    // std::cout << "127" << std::endl;
    return b;
  }

  template <int64_t mask>
  static c10::guts::enable_if_t<mask == 15, Vec256<float>> __inline_attrs
  blend(const Vec256<float>& a, const Vec256<float>& b) {
    // here I am using intel style mask number
    // std::cout << "15" << std::endl;
    return Vec256<float>{b._vec0, a._vec1};
  }

  template <int64_t mask>
  static c10::guts::enable_if_t<(mask > 0 && mask < 15), Vec256<float>>
      __inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
    constexpr uint32_t g0 = (mask & 1) * 0xffffffff;
    constexpr uint32_t g1 = ((mask & 2) >> 1) * 0xffffffff;
    constexpr uint32_t g2 = ((mask & 4) >> 2) * 0xffffffff;
    constexpr uint32_t g3 = ((mask & 8) >> 3) * 0xffffffff;
    const __vib mask_1st = (__vib){g0, g1, g2, g3};
    // std::cout << "mask 1st " << std::endl;
    return Vec256<float>{(__vf)vec_sel(a._vec0, b._vec0, (__vib)mask_1st),
                         a._vec1};
  }

  template <int64_t mask>
  static c10::guts::enable_if_t<
      (mask > 15 && (mask & 255) != 255 && ((mask & 15) == 15)),
      Vec256<float>>
      __inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
    // here I am using intel style mask number
    constexpr uint32_t mask2 = (mask & 255) >> 4;
    constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
    constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
    constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
    constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;
    // std::cout << "mask 2nd 1st b " << std::endl;
    const __vib mask_2nd = (__vib){g0_2, g1_2, g2_2, g3_2};
    // generated masks
    return Vec256<float>{b._vec0,
                         (__vf)vec_sel(a._vec1, b._vec1, (__vib)mask_2nd)};
  }

  template <int64_t mask>
  static c10::guts::enable_if_t<
      (mask > 15 && ((mask & 255) != 255) && ((mask & 15) == 0)),
      Vec256<float>>
      __inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
    // here I am using intel style mask number
    constexpr uint32_t mask2 = (mask & 255) >> 4;
    constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
    constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
    constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
    constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;
    // std::cout << "mask only 2nd " << std::endl;
    const __vib mask_2nd = (__vib){g0_2, g1_2, g2_2, g3_2};
    // generated masks
    return Vec256<float>{a, (__vf)vec_sel(a._vec1, b._vec1, (__vib)mask_2nd)};
  }

  template <int64_t mask>
  static c10::guts::enable_if_t<
      (mask > 15 && ((mask & 255) != 255) && ((mask & 15) != 0) &&
       ((mask & 15) != 15)),
      Vec256<float>>
      __inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
    // here I am using intel style mask number
    constexpr uint32_t g0 = (mask & 1) * 0xffffffff;
    constexpr uint32_t g1 = ((mask & 2) >> 1) * 0xffffffff;
    constexpr uint32_t g2 = ((mask & 4) >> 2) * 0xffffffff;
    constexpr uint32_t g3 = ((mask & 8) >> 3) * 0xffffffff;
    constexpr uint32_t mask2 = (mask & 255) >> 4;
    constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
    constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
    constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
    constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;
    std::cout << "mask both " << std::endl;
    const __vib mask_1st = (__vib){g0, g1, g2, g3};
    const __vib mask_2nd = (__vib){g0_2, g1_2, g2_2, g3_2};
    // generated masks
    return Vec256<float>{(__vf)vec_sel(a._vec0, b._vec0, (__vib)mask_1st),
                         (__vf)vec_sel(a._vec1, b._vec1, (__vib)mask_2nd)};
  }

  static Vec256<float> __inline_attrs blendv(
      const Vec256<float>& a,
      const Vec256<float>& b,
      const Vec256<float>& mask) {
    // the mask used here returned by comparision of vec256
    // assuming this we can use the same mask directly with vec_sel
    // warning intel style mask will not work properly
    return Vec256<float>{vec_sel(a._vec0, b._vec0, mask._vecb0),
                         vec_sel(a._vec1, b._vec1, mask._vecb1)};
  }

  static Vec256<float> arange(float base = 0.f, float step = 1.f) {
    return Vec256<float>(
        base,
        base + step,
        base + 2 * step,
        base + 3 * step,
        base + 4 * step,
        base + 5 * step,
        base + 6 * step,
        base + 7 * step);
  }
  static Vec256<float> set(
      const Vec256<float>& a,
      const Vec256<float>& b,
      int64_t count = size()) {
    switch (count) {
      case 0:
        return a;
      case 1:
        return blend<1>(a, b);
      case 2:
        return blend<3>(a, b);
      case 3:
        return blend<7>(a, b);
      case 4:
        return blend<15>(a, b);
      case 5:
        return blend<31>(a, b);
      case 6:
        return blend<63>(a, b);
      case 7:
        return blend<127>(a, b);
    }

    return b;
  }
  static Vec256<float> __inline_attrs
  loadu(const void* ptr, int64_t count = size()) {
    if (count == size()) {
      return Vec256<float>{
          vec_vsx_ld(offset0, reinterpret_cast<const float*>(ptr)),
          vec_vsx_ld(offset16, reinterpret_cast<const float*>(ptr))};
    }

    __at_align32__ float tmp_values[size()];
    __vf* vtmp = reinterpret_cast<__vf*>(tmp_values);
    std::memcpy(
        tmp_values, reinterpret_cast<const float*>(ptr), count * sizeof(float));

    return Vec256<float>{vtmp[0], vtmp[1]};
  }
  void __inline_attrs store(void* ptr, int count = size()) const {
    if (count == size()) {
      vec_vsx_st(_vec0, offset0, reinterpret_cast<float*>(ptr));
      vec_vsx_st(_vec1, offset16, reinterpret_cast<float*>(ptr));
    } else if (count > 0) {
      float tmp_values[size()];
      __vf* vtmp = reinterpret_cast<__vf*>(tmp_values);
      vtmp[0] = _vec0;
      vtmp[1] = _vec1;
      std::memcpy(ptr, tmp_values, count * sizeof(float));
    }
  }
  const float& operator[](int idx) const = delete;
  float& operator[](int idx) = delete;

  Vec256<float> map(float (*f)(float)) const {
    return Vec256<float>{f(_vec0[0]),
                         f(_vec0[1]),
                         f(_vec0[2]),
                         f(_vec0[3]),
                         f(_vec1[0]),
                         f(_vec1[1]),
                         f(_vec1[2]),
                         f(_vec1[3])};
  }

  Vec256<float> mapbi(float (*f)(float, float), const Vec256<float>& other)
      const {
    return Vec256<float>{f(_vec0[0], other._vec0[0]),
                         f(_vec0[1], other._vec0[1]),
                         f(_vec0[2], other._vec0[2]),
                         f(_vec0[3], other._vec0[3]),
                         f(_vec1[0], other._vec1[0]),
                         f(_vec1[1], other._vec1[1]),
                         f(_vec1[2], other._vec1[2]),
                         f(_vec1[3], other._vec1[3])};
  }

  Vec256<float> __inline_attrs abs() const {
    return Vec256<float>{vec_abs(_vec0), vec_abs(_vec1)};
  }

  Vec256<float> __inline_attrs acos() const {
    return map(std::acos);
  }
  Vec256<float> __inline_attrs asin() const {
    return map(std::asin);
  }
  Vec256<float> atan() const {
    return map(std::atan);
  }
  Vec256<float> atan2(const Vec256<float>& exp) const {
    return mapbi(std::atan2, exp);
  }

  Vec256<float> lgamma() const {
    return map(std::lgamma);
  }
  Vec256<float> erf() const {
    return map(std::erf);
  }
  Vec256<float> erfc() const {
    return map(std::erfc);
  }

  Vec256<float> erfinv() const {
    return map(calc_erfinv);
  }

  Vec256<float> angle() const {
    return Vec256<float>{0};
  }
  Vec256<float> real() const {
    return *this;
  }
  Vec256<float> imag() const {
    return Vec256<float>{0};
  }
  Vec256<float> conj() const {
    return *this;
  }

  Vec256<float> __inline_attrs exp() const {
    // implementation logic from avx_mathfun with some modifications
    __vf x0, x1, fx0, fx1, z0, z1, y0, y1, pow2n0, pow2n1, tmp0, tmp1;
    __vi imm0, imm1;

    const __vf one = vec_splats(1.f);
    const __vf exp_hi = vec_splats(88.3762626647949f);
    const __vf exp_lo = vec_splats(-88.3762626647949f);
    const __vf cephes_LOG2EF = vec_splats(1.44269504088896341f);
    const __vf cephes_exp_C1 = vec_splats(0.693359375f);
    const __vf cephes_exp_C2 = vec_splats(-2.12194440e-4f);
    const __vf cephes_exp_p0 = vec_splats(1.9875691500E-4f);
    const __vf cephes_exp_p1 = vec_splats(1.3981999507E-3f);
    const __vf cephes_exp_p2 = vec_splats(8.3334519073E-3f);
    const __vf cephes_exp_p3 = vec_splats(4.1665795894E-2f);
    const __vf cephes_exp_p4 = vec_splats(1.6666665459E-1f);
    const __vf cephes_exp_p5 = vec_splats(5.0000001201E-1f);
    const __vui v23 = vec_splats(23u);
    const __vi v0x7f = vec_splats(0x7f);
    __vf zerov = {0.f, 0.f, 0.f, 0.f};
    x0 = vec_min(_vec0, exp_hi);
    x1 = vec_min(_vec1, exp_hi);

    x0 = vec_max(_vec0, exp_lo);
    x1 = vec_max(_vec1, exp_lo);

    /* express exp(x) as exp(g + n*log(2)) */
    fx0 = vec_mul(x0, cephes_LOG2EF);
    fx1 = vec_mul(x1, cephes_LOG2EF);
    fx0 = vec_round(fx0);
    fx1 = vec_round(fx1);
    tmp0 = vec_mul(fx0, cephes_exp_C1);
    tmp1 = vec_mul(fx1, cephes_exp_C1);

    z0 = vec_mul(fx0, cephes_exp_C2);
    z1 = vec_mul(fx1, cephes_exp_C2);

    x0 = vec_sub(x0, tmp0);
    x1 = vec_sub(x1, tmp1);

    x0 = vec_sub(x0, z0);
    x1 = vec_sub(x1, z1);

    z0 = vec_mul(x0, x0);
    z1 = vec_mul(x1, x1);

    y0 = vec_madd(cephes_exp_p0, x0, cephes_exp_p1);
    y1 = vec_madd(cephes_exp_p1, x1, cephes_exp_p1);

    y0 = vec_madd(y0, x0, cephes_exp_p2);
    y1 = vec_madd(y1, x1, cephes_exp_p2);

    y0 = vec_madd(y0, x0, cephes_exp_p3);
    y1 = vec_madd(y1, x1, cephes_exp_p3);

    y0 = vec_madd(y0, x0, cephes_exp_p4);
    y1 = vec_madd(y1, x1, cephes_exp_p4);

    y0 = vec_madd(y0, x0, cephes_exp_p5);
    y1 = vec_madd(y1, x1, cephes_exp_p5);

    y0 = vec_madd(y0, z0, x0);
    y1 = vec_madd(y1, z1, x1);

    y0 = vec_add(y0, one);
    y1 = vec_add(y1, one);
    /* build 2^n */
    imm0 = vec_signed(
        fx0); //__vi{ static_cast<int>(fx0[0]),static_cast<int>(fx0[1])
              //,static_cast<int>(fx0[2]) ,static_cast<int>(fx0[3]) };
    imm1 = vec_signed(fx1);
    imm0 = vec_add(imm0, v0x7f);
    imm1 = vec_add(imm1, v0x7f);
    imm0 = vec_sl(imm0, v23);
    imm1 = vec_sl(imm1, v23);
    // treat imm0 as float vector without conversion
    pow2n0 = (__vf)imm0;
    pow2n1 = (__vf)imm1;
    y0 = vec_mul(y0, pow2n0);
    y1 = vec_mul(y1, pow2n1);
    return Vec256<float>{y0, y1};
  }
  Vec256<float> expm1() const {
    const __vf one = vec_splats(1.f);
    auto ret = exp();
    return Vec256<float>{
        vec_sub(ret._vec0, one),
        vec_sub(ret._vec1, one),
    };
  }

  Vec256<float> __inline_attrs log() const {
    __vi imm0, imm1;
    __vf x0, x1, e0, e1, ee0, ee1, y0, y1, z0, z1, tmp0, tmp1;
    __vib mask0, mask1;
    const __vf one = vec_splats(1.f);
    const __vui min_norm_pos = vec_splats((unsigned int)0x00800000);
    const __vui inv_mant_mask = vec_splats((unsigned int)(~0x7f800000));
    const __vi v0x7f = vec_splats(0x7f);
    const __vui v23 = vec_splats(23u);
    const __vf _0p5 = vec_splats(0.5f);
    const __vf cephes_SQRTHF = vec_splats(0.707106781186547524f);
    const __vf cephes_log_p0 = vec_splats(7.0376836292E-2f);
    const __vf cephes_log_p1 = vec_splats(-1.1514610310E-1f);
    const __vf cephes_log_p2 = vec_splats(1.1676998740E-1f);
    const __vf cephes_log_p3 = vec_splats(-1.2420140846E-1f);
    const __vf cephes_log_p4 = vec_splats(+1.4249322787E-1f);
    const __vf cephes_log_p5 = vec_splats(-1.6668057665E-1f);
    const __vf cephes_log_p6 = vec_splats(+2.0000714765E-1f);
    const __vf cephes_log_p7 = vec_splats(-2.4999993993E-1f);
    const __vf cephes_log_p8 = vec_splats(+3.3333331174E-1f);
    const __vf cephes_log_q1 = vec_splats(-2.12194440e-4f);
    const __vf cephes_log_q2 = vec_splats(0.693359375f);
    __vf zero = {0.f, 0.f, 0.f, 0.f};

    __vib invalid_mask0 = (__vib)vec_cmple(_vec0, zero);
    __vib invalid_mask1 = (__vib)vec_cmple(_vec1, zero);

    x0 = vec_max(_vec0, (__vf)min_norm_pos); /* cut off denormalized stuff */
    x1 = vec_max(_vec1, (__vf)min_norm_pos); /* cut off denormalized stuff */

    imm0 = vec_sr(__vi(x0), v23);
    imm1 = vec_sr(__vi(x1), v23);

    /* keep only the fractional part */
    x0 = vec_and(x0, (__vib)inv_mant_mask);
    x1 = vec_and(x1, (__vib)inv_mant_mask);
    x0 = vec_or(x0, _0p5);
    x1 = vec_or(x1, _0p5);

    imm0 = vec_sub(imm0, v0x7f);
    imm1 = vec_sub(imm1, v0x7f);
    e0 = vec_float(imm0);
    e1 = vec_float(imm1);

    e0 = vec_add(e0, one);
    e1 = vec_add(e1, one);
    /* part2:
                              if( x < SQRTHF ) {
                                    e -= 1;
                                    x = x + x - 1.0;
                              } else { x = x - 1.0; }
                       */
    mask0 = (__vib)vec_cmplt(x0, cephes_SQRTHF);
    mask1 = (__vib)vec_cmplt(x1, cephes_SQRTHF);
    tmp0 = vec_and(x0, mask0);
    tmp1 = vec_and(x1, mask1);
    x0 = vec_sub(x0, one);
    x1 = vec_sub(x1, one);
    e0 = vec_sub(e0, vec_and(one, mask0));
    e1 = vec_sub(e1, vec_and(one, mask1));
    x0 = vec_add(x0, tmp0);
    x1 = vec_add(x1, tmp1);

    z0 = vec_mul(x0, x0);
    z1 = vec_mul(x1, x1);

    y0 = vec_madd(x0, cephes_log_p0, cephes_log_p1);
    y1 = vec_madd(x1, cephes_log_p0, cephes_log_p1);
    y0 = vec_madd(y0, x0, cephes_log_p2);
    y1 = vec_madd(y1, x1, cephes_log_p2);
    y0 = vec_madd(y0, x0, cephes_log_p3);
    y1 = vec_madd(y1, x1, cephes_log_p3);
    y0 = vec_madd(y0, x0, cephes_log_p4);
    y1 = vec_madd(y1, x1, cephes_log_p4);
    y0 = vec_madd(y0, x0, cephes_log_p5);
    y1 = vec_madd(y1, x1, cephes_log_p5);
    y0 = vec_madd(y0, x0, cephes_log_p6);
    y1 = vec_madd(y1, x1, cephes_log_p6);
    y0 = vec_madd(y0, x0, cephes_log_p7);
    y1 = vec_madd(y1, x1, cephes_log_p7);
    y0 = vec_madd(y0, x0, cephes_log_p8);
    y1 = vec_madd(y1, x1, cephes_log_p8);
    y0 = vec_mul(y0, x0);
    y1 = vec_mul(y1, x1);

    y0 = vec_mul(y0, z0);
    y1 = vec_mul(y1, z1);

    y0 = vec_madd(e0, cephes_log_q1, y0);
    y1 = vec_madd(e1, cephes_log_q1, y1);

    y0 = y0 - z0 * _0p5;
    y1 = y1 - z1 * _0p5;
    x0 = vec_add(x0, y0);
    x1 = vec_add(x1, y1);
    x0 = vec_madd(e0, cephes_log_q2, x0);
    x1 = vec_madd(e1, cephes_log_q2, x1);
    x0 = vec_or(x0, invalid_mask0); // negative arg will be NAN
    x1 = vec_or(x1, invalid_mask1);
    return Vec256<float>{x0, x1};
  }
  Vec256<float> __inline_attrs log10() const {
    const __vf log10e_inv = vec_splats(0.43429448190325176f);
    auto ret = log();
    return Vec256<float>{vec_mul(ret._vec0, log10e_inv),
                         vec_mul(ret._vec1, log10e_inv)};
  }
  Vec256<float> __inline_attrs log1p() const {
    const __vf one = vec_splats(1.f);
    Vec256<float> ret = {vec_add(_vec0, one), vec_add(_vec1, one)};

    return ret.log();
  }
  Vec256<float> __inline_attrs log2() const {
    const __vf log2e_inv = vec_splats(1.4426950408889634f);
    auto ret = log();
    return Vec256<float>{vec_mul(ret._vec0, log2e_inv),
                         vec_mul(ret._vec1, log2e_inv)};
  }
  Vec256<float> __inline_attrs ceil() const {
    return Vec256<float>{vec_ceil(_vec0), vec_ceil(_vec1)};
  }
  Vec256<float> __inline_attrs cos() const {
    const __vui v29 = vec_splats(29u);
    const __vib sign_mask = (__vib)vec_splats((unsigned int)0x80000000);
    const __vi vi_0 = vec_splats((int)(0));
    const __vi vi_1 = vec_splats((int)1);
    const __vi vi_inv1 = vec_splats((int)~1);
    const __vi vi_2 = vec_splats((int)2);
    const __vi vi_4 = vec_splats((int)4);

    const __vf _0p5 = vec_splats(0.5f);
    const __vf one = vec_splats(1.f);
    const __vf minus_cephes_dp1 = vec_splats(-0.78515625f);
    const __vf minus_cephes_dp2 = vec_splats(-2.4187564849853515625e-4f);
    const __vf minus_cephes_dp3 = vec_splats(-3.77489497744594108e-8f);

    const __vf sincof_p0 = vec_splats(-1.9515295891E-4f);
    const __vf sincof_p1 = vec_splats(8.3321608736E-3f);
    const __vf sincof_p2 = vec_splats(-1.6666654611E-1f);

    const __vf coscof_p0 = vec_splats(2.443315711809948E-005f);
    const __vf coscof_p1 = vec_splats(-1.388731625493765E-003f);
    const __vf coscof_p2 = vec_splats(4.166664568298827E-002f);
    const __vf _4div_pi = vec_splats(1.27323954473516f);

    __vf x0, x1, y0, y1, z0, z1, y0_2, y1_2;
    __vib sign_bit0, sign_bit1, poly_mask0, poly_mask1;
    __vi imm0_0, imm1_0, imm0_2, imm1_2;

    /* take the absolute value */
    x0 = vec_abs(_vec0);
    x1 = vec_abs(_vec1);
    /* extract the sign bit (upper one) */
    sign_bit0 = (__vib)vec_and(_vec0, sign_mask);
    sign_bit1 = (__vib)vec_and(_vec1, sign_mask);
    /* scale by 4/Pi */
    y0 = vec_mul(x0, _4div_pi);
    y1 = vec_mul(x1, _4div_pi);
    /* store the integer part of y in mm0 */
    imm0_2 = vec_signed(y0);
    imm1_2 = vec_signed(y1);
    /* j=(j+1) & (~1) (see the cephes sources) */

    imm0_2 = vec_add(imm0_2, vi_1);
    imm1_2 = vec_add(imm1_2, vi_1);
    imm0_2 = vec_and(imm0_2, (__vib)vi_inv1);
    imm1_2 = vec_and(imm1_2, (__vib)vi_inv1);
    y0 = vec_float(imm0_2);
    y1 = vec_float(imm1_2);
    imm0_2 = vec_sub(imm0_2, (__vi)vi_2);
    imm1_2 = vec_sub(imm1_2, (__vi)vi_2);
    /* get the swap sign flag */
    imm0_0 = vec_and(vec_nand(imm0_2, imm0_2), vi_4);
    imm1_0 = vec_and(vec_nand(imm1_2, imm1_2), vi_4);
    sign_bit0 = (__vib)vec_sl(imm0_0, v29);
    sign_bit1 = (__vib)vec_sl(imm1_0, v29);
    /* get the polynom selection mask
       there is one polynom for 0 <= x <= Pi/4
       and another one for Pi/4<x<=Pi/2

       Both branches will be computed.
    */
    imm0_2 = vec_and(imm0_2, (__vib)vi_2);
    imm1_2 = vec_and(imm1_2, (__vib)vi_2);
    poly_mask0 = (__vib)vec_cmpeq(imm0_2, (__vi)vi_0);
    poly_mask1 = (__vib)vec_cmpeq(imm1_2, (__vi)vi_0);

    /* The magic pass: "Extended precision modular arithmetic"
       x = ((x - y * DP1) - y * DP2) - y * DP3; */

    x0 = vec_madd(y0, minus_cephes_dp1, x0);
    x1 = vec_madd(y1, minus_cephes_dp1, x1);
    x0 = vec_madd(y0, minus_cephes_dp2, x0);
    x1 = vec_madd(y1, minus_cephes_dp2, x1);
    x0 = vec_madd(y0, minus_cephes_dp3, x0);
    x1 = vec_madd(y1, minus_cephes_dp3, x1);

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */

    z0 = vec_mul(x0, x0);
    z1 = vec_mul(x1, x1);

    y0 = vec_madd(coscof_p0, z0, coscof_p1);
    y1 = vec_madd(coscof_p0, z1, coscof_p1);
    y0 = vec_madd(y0, z0, coscof_p2);
    y1 = vec_madd(y1, z1, coscof_p2);
    y0 = vec_mul(y0, z0);
    y1 = vec_mul(y1, z1);
    y0 = vec_mul(y0, z0);
    y1 = vec_mul(y1, z1);
    y0 = y0 - z0 * _0p5;
    y1 = y1 - z1 * _0p5;
    y0 = vec_add(y0, one);
    y1 = vec_add(y1, one);

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

    y0_2 = vec_madd(sincof_p0, z0, sincof_p1);
    y1_2 = vec_madd(sincof_p0, z1, sincof_p1);
    y0_2 = vec_madd(y0_2, z0, sincof_p2);
    y1_2 = vec_madd(y1_2, z1, sincof_p2);
    y0_2 = vec_mul(y0_2, z0);
    y1_2 = vec_mul(y1_2, z1);
    y0_2 = vec_madd(y0_2, x0, x0);
    y1_2 = vec_madd(y1_2, x1, x1);

    /* select the correct result from the two polynoms */
    y0 = vec_sel(y0, y0_2, poly_mask0);
    y1 = vec_sel(y1, y1_2, poly_mask1);
    /* update the sign */
    y0 = vec_xor(y0, sign_bit0);
    y1 = vec_xor(y1, sign_bit1);

    return Vec256<float>{y0, y1};
  }
  Vec256<float> __inline_attrs cosh() const {
    // cosh = 1/2 * (e^x + e^-x)
    Vec256<float> ret = this->exp();
    Vec256<float> ret_recp = ret.reciprocal();
    const __vf half = vec_splats(0.5f);
    __vf v0 = vec_mul(half, vec_add(ret._vec0, ret_recp._vec0));
    __vf v1 = vec_mul(half, vec_add(ret._vec1, ret_recp._vec1));
    return Vec256<float>{v0, v1};
  }
  Vec256<float> __inline_attrs floor() const {
    return Vec256<float>{vec_floor(_vec0), vec_floor(_vec1)};
  }
  Vec256<float> __inline_attrs neg() const {
    return Vec256<float>{vec_neg(_vec0), vec_neg(_vec1)};
  }
  Vec256<float> __inline_attrs round() const {
    return Vec256<float>{vec_round(_vec0), vec_round(_vec1)};
  }
  Vec256<float> __inline_attrs sin() const {
    const __vui v29 = vec_splats(29u);
    const __vib sign_mask = (__vib)vec_splats((unsigned int)0x80000000);
    const __vui vuint_0 = vec_splats(0u);
    const __vi vi_1 = vec_splats((int)1);
    const __vui vuint_inv1 = vec_splats((unsigned int)~1);
    const __vui vuint_2 = vec_splats(2u);
    const __vui vuint_4 = vec_splats(4u);

    const __vf _0p5 = vec_splats(0.5f);
    const __vf one = vec_splats(1.f);
    const __vf minus_cephes_dp1 = vec_splats(-0.78515625f);
    const __vf minus_cephes_dp2 = vec_splats(-2.4187564849853515625e-4f);
    const __vf minus_cephes_dp3 = vec_splats(-3.77489497744594108e-8f);

    const __vf sincof_p0 = vec_splats(-1.9515295891E-4f);
    const __vf sincof_p1 = vec_splats(8.3321608736E-3f);
    const __vf sincof_p2 = vec_splats(-1.6666654611E-1f);

    const __vf coscof_p0 = vec_splats(2.443315711809948E-005f);
    const __vf coscof_p1 = vec_splats(-1.388731625493765E-003f);
    const __vf coscof_p2 = vec_splats(4.166664568298827E-002f);
    const __vf _4div_pi = vec_splats(1.27323954473516f);

    __vf x0, x1, y0, y1, z0, z1, y0_2, y1_2;
    __vib swap_sign_bit0, swap_sign_bit1, sign_bit0, sign_bit1, poly_mask0,
        poly_mask1;
    __vi imm0_0, imm1_0, imm0_2, imm1_2;

    /* take the absolute value */
    x0 = vec_abs(_vec0);
    x1 = vec_abs(_vec1);
    /* extract the sign bit (upper one) */
    sign_bit0 = (__vib)vec_and(_vec0, sign_mask);
    sign_bit1 = (__vib)vec_and(_vec1, sign_mask);
    /* scale by 4/Pi */
    y0 = vec_mul(x0, _4div_pi);
    y1 = vec_mul(x1, _4div_pi);
    /* store the integer part of y in mm0 */
    imm0_2 = vec_signed(y0);
    imm1_2 = vec_signed(y1);
    /* j=(j+1) & (~1) (see the cephes sources) */

    imm0_2 = vec_add(imm0_2, vi_1);
    imm1_2 = vec_add(imm1_2, vi_1);
    imm0_2 = vec_and(imm0_2, (__vib)vuint_inv1);
    imm1_2 = vec_and(imm1_2, (__vib)vuint_inv1);
    y0 = vec_float(imm0_2);
    y1 = vec_float(imm1_2);
    /* get the swap sign flag */
    imm0_0 = vec_and(imm0_2, (__vib)vuint_4);
    imm1_0 = vec_and(imm1_2, (__vib)vuint_4);
    swap_sign_bit0 = (__vib)vec_sl(imm0_0, v29);
    swap_sign_bit1 = (__vib)vec_sl(imm1_0, v29);
    /* get the polynom selection mask
       there is one polynom for 0 <= x <= Pi/4
       and another one for Pi/4<x<=Pi/2

       Both branches will be computed.
    */
    imm0_2 = vec_and(imm0_2, (__vib)vuint_2);
    imm1_2 = vec_and(imm1_2, (__vib)vuint_2);
    poly_mask0 = (__vib)vec_cmpeq(imm0_2, (__vi)vuint_0);
    poly_mask1 = (__vib)vec_cmpeq(imm1_2, (__vi)vuint_0);

    sign_bit0 = vec_xor(sign_bit0, swap_sign_bit0);
    sign_bit1 = vec_xor(sign_bit1, swap_sign_bit1);
    /* The magic pass: "Extended precision modular arithmetic"
       x = ((x - y * DP1) - y * DP2) - y * DP3; */

    x0 = vec_madd(y0, minus_cephes_dp1, x0);
    x1 = vec_madd(y1, minus_cephes_dp1, x1);
    x0 = vec_madd(y0, minus_cephes_dp2, x0);
    x1 = vec_madd(y1, minus_cephes_dp2, x1);
    x0 = vec_madd(y0, minus_cephes_dp3, x0);
    x1 = vec_madd(y1, minus_cephes_dp3, x1);

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */

    z0 = vec_mul(x0, x0);
    z1 = vec_mul(x1, x1);

    y0 = vec_madd(coscof_p0, z0, coscof_p1);
    y1 = vec_madd(coscof_p0, z1, coscof_p1);
    y0 = vec_madd(y0, z0, coscof_p2);
    y1 = vec_madd(y1, z1, coscof_p2);
    y0 = vec_mul(y0, z0);
    y1 = vec_mul(y1, z1);
    y0 = vec_mul(y0, z0);
    y1 = vec_mul(y1, z1);
    y0 = y0 - z0 * _0p5;
    y1 = y1 - z1 * _0p5;
    y0 = vec_add(y0, one);
    y1 = vec_add(y1, one);

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

    y0_2 = vec_madd(sincof_p0, z0, sincof_p1);
    y1_2 = vec_madd(sincof_p0, z1, sincof_p1);
    y0_2 = vec_madd(y0_2, z0, sincof_p2);
    y1_2 = vec_madd(y1_2, z1, sincof_p2);
    y0_2 = vec_mul(y0_2, z0);
    y1_2 = vec_mul(y1_2, z1);
    y0_2 = vec_madd(y0_2, x0, x0);
    y1_2 = vec_madd(y1_2, x1, x1);

    /* select the correct result from the two polynoms */
    y0 = vec_sel(y0, y0_2, poly_mask0);
    y1 = vec_sel(y1, y1_2, poly_mask1);
    /* update the sign */
    y0 = vec_xor(y0, sign_bit0);
    y1 = vec_xor(y1, sign_bit1);

    return Vec256<float>{y0, y1};
  }
  Vec256<float> __inline_attrs sinh() const {
    // cosh = 1/2 * (e^x - e^-x)
    const __vf half = vec_splats(0.5f);
    Vec256<float> ret = this->exp();
    Vec256<float> ret_recp = ret.reciprocal();
    __vf v0 = vec_mul(half, vec_sub(ret._vec0, ret_recp._vec0));
    __vf v1 = vec_mul(half, vec_sub(ret._vec1, ret_recp._vec1));
    return Vec256<float>{v0, v1};
  }
  Vec256<float> __inline_attrs tan() const {
    return map(std::tan);
  }
  Vec256<float> __inline_attrs tanh() const {
    // (e^(2x) -1 ) / (e^(2x) + 1);
    const __vf one = vec_splats(1.0f);
    const __vf two = vec_splats(2.0f);
    Vec256<float> vec_tmp{vec_mul(two, _vec0), vec_mul(two, _vec1)};
    Vec256<float> exp2x = vec_tmp.exp();
    return Vec256<float>{
        vec_div(vec_sub(exp2x._vec0, one), vec_add(exp2x._vec0, one)),
        vec_div(vec_sub(exp2x._vec1, one), vec_add(exp2x._vec1, one))};
  }
  Vec256<float> __inline_attrs trunc() const {
    return Vec256<float>{vec_trunc(_vec0), vec_trunc(_vec1)};
  }

  Vec256<float> __inline_attrs frac() const {
    return Vec256<float>{vec_sub(_vec0, vec_trunc(_vec0)),
                         vec_sub(_vec1, vec_trunc(_vec1))};
  }

  Vec256<float> __inline_attrs sqrt() const {
    return Vec256<float>{vec_sqrt(_vec0), vec_sqrt(_vec1)};
  }
  Vec256<float> __inline_attrs reciprocal() const {
    const __vf one = vec_splats(1.f);
    return Vec256<float>{
        vec_div(
            one, _vec0), // vec_re(_vec0) - is estimated calc so we dont use it
        vec_div(one, _vec1)};
  }
  Vec256<float> __inline_attrs rsqrt() const {
    const __vf one = vec_splats(1.f);
    return Vec256<float>{
        vec_div(one, vec_sqrt(_vec0)), // vec_rsqrt- is estimated rsqrt
        vec_div(one, vec_sqrt(_vec1))};
  }

  Vec256<float> __inline_attrs pow(const Vec256<float>& pow_exp) const {
    __vf tem0, tem1, out_0, out_1, ret0, ret1;
    __vf odd0, odd1;
    __vf floor_b0, floor_b1, out1_0, out1_1;
    __vib man0, man1, sign0, sign1;
    const __vui mask = vec_splats(0x80000000);
    const __vf one = vec_splats(1.f);
    const __vui nan = vec_splats(0xffffffff);

    const __vi uint_1 = vec_splats(1);
    const __vi uint_0 = vec_splats(0);
    __vf zero = {0.f, 0.f, 0.f, 0.f};

    sign0 = (__vib)vec_and((__vib)mask, _vec0);
    sign1 = (__vib)vec_and((__vib)mask, _vec1);

    // |b|
    auto abs_b = pow_exp.abs();

    /* using ln fuction */
    auto temp = abs();
    temp = temp.log();

    tem0 = vec_mul(temp._vec0, pow_exp._vec0);
    tem1 = vec_mul(temp._vec1, pow_exp._vec1);

    /* using exp fuction */
    auto temp1 = Vec256<float>{tem0, tem1};

    auto cc = temp1.exp();

    man0 = (__vib)vec_cmpeq(
        vec_and(vec_signed(pow_exp._vec0), uint_1), uint_0); // even or odd
    man1 = (__vib)vec_cmpeq(
        vec_and(vec_signed(pow_exp._vec1), uint_1), uint_0); // even or odd
    // if even then then pow result should be absolute

    odd0 = (__vf)vec_or(cc._vec0, sign0); // copy_sign
    odd1 = (__vf)vec_or(cc._vec1, sign1);

    out_0 = vec_sel(odd0, cc._vec0, man0);
    out_1 = vec_sel(odd1, cc._vec1, man1);

    floor_b0 = vec_floor(pow_exp._vec0);
    floor_b1 = vec_floor(pow_exp._vec1);

    /* x<0 and y != N, then NAN */
    man0 = (__vib)vec_and(
        vec_cmpne(pow_exp._vec0, floor_b0), vec_cmplt(_vec0, zero));
    man1 = (__vib)vec_and(
        vec_cmpne(pow_exp._vec1, floor_b1), vec_cmplt(_vec1, zero));

    out1_0 = vec_sel(out_0, (__vf)nan, man0);
    out1_1 = vec_sel(out_1, (__vf)nan, man1);
    /* y = 0 then 1 */
    man0 = (__vib)vec_cmpeq(abs_b._vec0, zero);
    man1 = (__vib)vec_cmpeq(abs_b._vec1, zero);
    ret0 = vec_sel(out1_0, one, man0);
    ret1 = vec_sel(out1_1, one, man1);
    return Vec256<float>{ret0, ret1};
  }

  Vec256<float> __inline_attrs operator==(const Vec256<float>& other) const {
    return Vec256<float>{(__vib)vec_cmpeq(_vec0, other._vec0),
                         (__vib)vec_cmpeq(_vec1, other._vec1)};
  }

  Vec256<float> __inline_attrs operator!=(const Vec256<float>& other) const {
    return Vec256<float>{(__vib)vec_cmpne(_vec0, other._vec0),
                         (__vib)vec_cmpne(_vec1, other._vec1)};
  }

  Vec256<float> __inline_attrs operator<(const Vec256<float>& other) const {
    return Vec256<float>{(__vib)vec_cmplt(_vec0, other._vec0),
                         (__vib)vec_cmplt(_vec1, other._vec1)};
  }

  Vec256<float> __inline_attrs operator<=(const Vec256<float>& other) const {
    return Vec256<float>{(__vib)vec_cmple(_vec0, other._vec0),
                         (__vib)vec_cmple(_vec1, other._vec1)};
  }

  Vec256<float> __inline_attrs operator>(const Vec256<float>& other) const {
    return Vec256<float>{(__vib)vec_cmpgt(_vec0, other._vec0),
                         (__vib)vec_cmpgt(_vec1, other._vec1)};
  }

  Vec256<float> __inline_attrs operator>=(const Vec256<float>& other) const {
    return Vec256<float>{(__vib)vec_cmpge(_vec0, other._vec0),
                         (__vib)vec_cmpge(_vec1, other._vec1)};
  }
};

#endif

} // namespace
} // namespace vec256
} // namespace at
