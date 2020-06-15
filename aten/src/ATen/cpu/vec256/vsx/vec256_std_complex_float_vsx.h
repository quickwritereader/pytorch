
#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>
#include <ATen/cpu/vec256/vsx/vsx_helpers.h>

namespace at {
namespace vec256 {
// See Note [Acceptable use of anonymous namespace in header]
namespace {
using StdComplexFlt = std::complex<float>;

template <>
class Vec256<StdComplexFlt> {
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
  using value_type = StdComplexFlt;
  using vec_internal_type = __vf;
  using vec_internal_mask_type = __vib;

  static constexpr int size() { return 4; }
  Vec256() {}

  __inline_attrs Vec256(__vf v) : _vec0{v}, _vec1{v} {}
  __inline_attrs Vec256(__vib vmask) : _vecb0{vmask}, _vecb1{vmask} {}
  __inline_attrs Vec256(__vf v1, __vf v2) : _vec0{v1}, _vec1{v2} {}
  __inline_attrs Vec256(__vib v1, __vib v2) : _vecb0{v1}, _vecb1{v2} {}

  Vec256(StdComplexFlt val) {
    float real_value = val.real();
    float imag_value = val.imag();
    _vec0 = __vf{real_value, imag_value, real_value, imag_value};
    _vec1 = __vf{real_value, imag_value, real_value, imag_value};
  }

  Vec256(StdComplexFlt val1, StdComplexFlt val2, StdComplexFlt val3, StdComplexFlt val4) {
    _vec0 = __vf{val1.real(), val1.imag(), val2.real(), val2.imag()};
    _vec1 = __vf{val3.real(), val3.imag(), val4.real(), val4.imag()};
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 0, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    return a;
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 1, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    return b;
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 2, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    return {b._vec0, a._vec1};
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 3, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    return {a._vec0, b._vec1};
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 4, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    const __vib mask_1st = VsxComplexMask1(mask);
    return {(__vf)vec_sel(a._vec0, b._vec0, mask_1st), a._vec1};
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 5, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    const __vib mask_1st = VsxComplexMask1(mask);
    return {(__vf)vec_sel(a._vec0, b._vec0, mask_1st), b._vec1};
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 6, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    const __vib mask_2nd = VsxComplexMask2(mask);
    // generated masks
    return {a._vec0, (__vf)vec_sel(a._vec1, b._vec1, mask_2nd)};
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 7, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    const __vib mask_2nd = VsxComplexMask2(mask);
    // generated masks
    return {b._vec0, (__vf)vec_sel(a._vec1, b._vec1, mask_2nd)};
  }

  template <uint64_t mask>
  static std::enable_if_t<blendChoiceComplex(mask) == 8, Vec256<StdComplexFlt>>
      __inline_attrs blend(const Vec256<StdComplexFlt>& a,
                           const Vec256<StdComplexFlt>& b) {
    const __vib mask_1st = VsxComplexMask1(mask);
    const __vib mask_2nd = VsxComplexMask2(mask);
    return {(__vf)vec_sel(a._vec0, b._vec0, mask_1st),
            (__vf)vec_sel(a._vec1, b._vec1, mask_2nd)};
  }

  template <int64_t mask>
  static Vec256<StdComplexFlt> __inline_attrs
  el_blend(const Vec256<StdComplexFlt>& a, const Vec256<StdComplexFlt>& b) {
    const __vib mask_1st = VsxMask1(mask);
    const __vib mask_2nd = VsxMask2(mask);
    return {(__vf)vec_sel(a._vec0, b._vec0, mask_1st),
            (__vf)vec_sel(a._vec1, b._vec1, mask_2nd)};
  }

  static Vec256<StdComplexFlt> blendv(const Vec256<StdComplexFlt>& a,
                                   const Vec256<StdComplexFlt>& b,
                                   const Vec256<StdComplexFlt>& mask) {
    // convert std::complex<V> index mask to V index mask: xy -> xxyy
    auto mask_complex = Vec256<StdComplexFlt>(vec_mergeh(mask._vec0, mask._vec0),
                                           vec_mergeh(mask._vec1, mask._vec1));
    // mask_complex.dump();
    return {
        vec_sel(a._vec0, b._vec0, mask_complex._vec0),
        vec_sel(a._vec1, b._vec1, mask_complex._vec1),
    };
  }

  static Vec256<StdComplexFlt> elwise_blendv(const Vec256<StdComplexFlt>& a,
                                          const Vec256<StdComplexFlt>& b,
                                          const Vec256<StdComplexFlt>& mask) {
    return {
        vec_sel(a._vec0, b._vec0, mask._vec0),
        vec_sel(a._vec1, b._vec1, mask._vec1),
    };
  }

  template <typename step_t>
  static Vec256<StdComplexFlt> arange(StdComplexFlt base = 0.,
                                   step_t step = static_cast<step_t>(1)) {
    return Vec256<StdComplexFlt>(base, base + step, base + StdComplexFlt(2) * step,
                              base + StdComplexFlt(3) * step);
  }
  static Vec256<StdComplexFlt> set(const Vec256<StdComplexFlt>& a,
                                const Vec256<StdComplexFlt>& b,
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
    }
    return b;
  }

  static Vec256<value_type> __inline_attrs loadu(const void* ptr,
                                                 int count = size()) {
    if (count == size()) {
      return {vec_vsx_ld(offset0, reinterpret_cast<const float*>(ptr)),
              vec_vsx_ld(offset16, reinterpret_cast<const float*>(ptr))};
    }

    __at_align32__ value_type tmp_values[size()];
    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

    return {vec_vsx_ld(offset0, reinterpret_cast<const float*>(tmp_values)),
            vec_vsx_ld(offset16, reinterpret_cast<const float*>(tmp_values))};
  }
  void __inline_attrs store(void* ptr, int count = size()) const {
    if (count == size()) {
      vec_vsx_st(_vec0, offset0, reinterpret_cast<float*>(ptr));
      vec_vsx_st(_vec1, offset16, reinterpret_cast<float*>(ptr));
    } else if (count > 0) {
      __at_align32__ value_type tmp_values[size()];
      vec_vsx_st(_vec0, offset0, reinterpret_cast<float*>(tmp_values));
      vec_vsx_st(_vec1, offset16, reinterpret_cast<float*>(tmp_values));
      std::memcpy(ptr, tmp_values,
                  std::min(count, size()) * sizeof(value_type));
    }
  }

  const StdComplexFlt& operator[](int idx) const = delete;
  StdComplexFlt& operator[](int idx) = delete;

  Vec256<StdComplexFlt> map(StdComplexFlt (*f)(const StdComplexFlt&)) const {
    __at_align32__ StdComplexFlt tmp[size()];
    store(tmp);
    for (int i = 0; i < size(); i++) {
      tmp[i] = f(tmp[i]);
    }
    return loadu(tmp);
  }

  static Vec256<StdComplexFlt> horizontal_add_permD8(Vec256<StdComplexFlt>& first,
                                                  Vec256<StdComplexFlt>& second) {
    // we will simulate it differently with 6 instructions total
    // lets permute second so that we can add it getting horizontall sums
    auto first_perm = first.el_swapped();    // 2perm
    auto second_perm = second.el_swapped();  // 2perm
    // summ
    auto first_ret = first + first_perm;     // 2add
    auto second_ret = second + second_perm;  // 2 add
    // now lets choose evens
    return el_mergee(first_ret, second_ret);  // 2 mergee's
  }

  static Vec256<StdComplexFlt> horizontal_sub_permD8(Vec256<StdComplexFlt>& first,
                                                  Vec256<StdComplexFlt>& second) {
    // we will simulate it differently with 6 instructions total
    // lets permute second so that we can add it getting horizontall sums
    auto first_perm = first.el_swapped();    // 2perm
    auto second_perm = second.el_swapped();  // 2perm
    // summ
    auto first_ret = first - first_perm;     // 2sub
    auto second_ret = second - second_perm;  // 2 sub
    // now lets choose evens
    return el_mergee(first_ret, second_ret);  // 2 mergee's
  }

  Vec256<StdComplexFlt> abs_2_() const {
    auto a = (*this).elwise_mult(*this);
    auto permuted = a.el_swapped();
    a = a + permuted;
    return a.el_mergee();
  }

  Vec256<StdComplexFlt> abs_() const {
    auto ret = abs_2_();
    return ret.elwise_sqrt();
  }

  Vec256<StdComplexFlt> abs() const { return abs_() & real_mask; }

  Vec256<StdComplexFlt> real_() const { return *this & real_mask; }
  Vec256<StdComplexFlt> real() const { return *this & real_mask; }
  Vec256<StdComplexFlt> imag_() const { return *this & imag_mask; }
  Vec256<StdComplexFlt> imag() const {
    // we can use swap_mask or sldwi
    auto ret = imag_();
    return {vec_sldw(ret._vec0, ret._vec0, 3),
            vec_sldw(ret._vec1, ret._vec1, 3)};
  }

  Vec256<StdComplexFlt> conj_() const { return *this ^ isign_mask; }
  Vec256<StdComplexFlt> conj() const { return *this ^ isign_mask; }

  Vec256<StdComplexFlt> log() const {
    // Most trigonomic ops use the log() op to improve complex number
    // performance.
    return map(std::log);
  }

  Vec256<StdComplexFlt> log2() const {
    // log2eB_inv
    auto ret = log();
    return ret.elwise_mult(log2e_inv);
  }
  Vec256<StdComplexFlt> log10() const {
    auto ret = log();
    return ret.elwise_mult(log10e_inv);
  }


  Vec256<StdComplexFlt> el_swapped() const {
    __vf v0 = vec_perm(_vec0, _vec0, swap_mask);
    __vf v1 = vec_perm(_vec1, _vec1, swap_mask);
    return {v0, v1};
  }

  Vec256<StdComplexFlt> el_mergee() const {
    // as mergee phased in , we can use vec_perm with mask
    return {vec_mergee(_vec0, _vec0), vec_mergee(_vec1, _vec1)};
  }

  static Vec256<StdComplexFlt> el_mergee(Vec256<StdComplexFlt>& first,
                                      Vec256<StdComplexFlt>& second) {
    // as mergee phased in , we can use vec_perm with mask
    return {vec_mergee(first._vec0, second._vec0),
            vec_mergee(first._vec1, second._vec1)};
  }

  Vec256<StdComplexFlt> angle_() const {
    // angle = atan2(b/a)
    // auto b_a = _mm256_permute_ps(values, 0xB1); // b        a
    // return Sleef_atan2f8_u10(values, b_a); // 90-angle angle
    auto ret = el_swapped();
    for (int i = 0; i < 4; i++) {
      ret._vec0[i] = std::atan2(_vec0[i], ret._vec0[i]);
      ret._vec1[i] = std::atan2(_vec1[i], ret._vec0[i]);
    }
    return ret;
  }

  Vec256<StdComplexFlt> angle() const {
    auto a = angle_().el_swapped();
    return a & real_mask;
  }


  Vec256<StdComplexFlt> sin() const { return map(std::sin); }
  Vec256<StdComplexFlt> sinh() const { return map(std::sinh); }
  Vec256<StdComplexFlt> cos() const { return map(std::cos); }
  Vec256<StdComplexFlt> cosh() const { return map(std::cosh); }
  Vec256<StdComplexFlt> ceil() const { return {vec_ceil(_vec0), vec_ceil(_vec1)}; }
  Vec256<StdComplexFlt> floor() const {
    return {vec_floor(_vec0), vec_floor(_vec1)};
  }
  Vec256<StdComplexFlt> neg() const {
    auto z = Vec256<StdComplexFlt>(zero);
    return z - *this;
  }
  Vec256<StdComplexFlt> round() const {
    return {vec_round(_vec0), vec_round(_vec1)};
  }
  Vec256<StdComplexFlt> tan() const { return map(std::tan); }
  Vec256<StdComplexFlt> tanh() const { return map(std::tanh); }
  Vec256<StdComplexFlt> trunc() const {
    return {vec_trunc(_vec0), vec_trunc(_vec1)};
  }

  Vec256<StdComplexFlt> elwise_sqrt() const {
    return {vec_sqrt(_vec0), vec_sqrt(_vec1)};
  }

  void dump() const {
    std::cout << _vec0[0] << "," << _vec0[1] << "," << _vec0[2] << ","
              << _vec0[3] << ",";
    std::cout << _vec1[0] << "," << _vec1[1] << "," << _vec1[2] << ","
              << _vec1[3] << std::endl;
  }

  Vec256<StdComplexFlt> sqrt() const {
    //   sqrt(a + bi)
    // = sqrt(2)/2 * [sqrt(sqrt(a**2 + b**2) + a) + sgn(b)*sqrt(sqrt(a**2 +
    // b**2) - a)i] = sqrt(2)/2 * [sqrt(abs() + a) + sgn(b)*sqrt(abs() - a)i]

    auto sign = *this & isign_mask;
    auto factor = sign | sqrt2_2;
    auto a_a = el_mergee();
    // a_a.dump();
    a_a = a_a ^ isign_mask;  // a -a
    auto res_re_im =
        (abs_() + a_a).elwise_sqrt();  // sqrt(abs + a) sqrt(abs - a)
    return factor.elwise_mult(res_re_im);
  }

  Vec256<StdComplexFlt> reciprocal() const {
    // re + im*i = (a + bi)  / (c + di)
    // re = (ac + bd)/abs_2() = c/abs_2()
    // im = (bc - ad)/abs_2() = d/abs_2()
    auto c_d = *this ^ isign_mask;  // c       -d
    auto abs = abs_2_();
    return c_d.elwise_div(abs);
  }

  Vec256<StdComplexFlt> rsqrt() const { return sqrt().reciprocal(); }

  Vec256<StdComplexFlt> pow(const Vec256<StdComplexFlt>& exp) const {
    __at_align32__ StdComplexFlt x_tmp[size()];
    __at_align32__ StdComplexFlt y_tmp[size()];
    store(x_tmp);
    exp.store(y_tmp);
    for (int i = 0; i < size(); i++) {
      x_tmp[i] = std::pow(x_tmp[i], y_tmp[i]);
    }
    return loadu(x_tmp);
  }

  Vec256<StdComplexFlt> atan() const {
    // atan(x) = i/2 * ln((i + z)/(i - z))
    auto ione = Vec256(imag_one);
    auto sum = ione + *this;
    auto sub = ione - *this;
    auto ln = (sum / sub).log();  // ln((i + z)/(i - z))
    return ln * imag_half;        // i/2*ln()
  }

  Vec256<StdComplexFlt> acos() const {
    // acos(x) = pi/2 - asin(x)
    return Vec256(pi_2) - asin();
  }

  Vec256<StdComplexFlt> inline operator*(const Vec256<StdComplexFlt>& b) const {
    //(a + bi)  * (c + di) = (ac - bd) + (ad + bc)i
    auto ac_bd = elwise_mult(b);
    auto d_c = b.el_swapped();
    d_c = d_c ^ isign_mask;
    auto ad_bc = elwise_mult(d_c);
    auto ret = horizontal_sub_permD8(ac_bd, ad_bc);
    return ret;
  }

  Vec256<StdComplexFlt> inline operator/(const Vec256<StdComplexFlt>& b) const {
    // re + im*i = (a + bi)  / (c + di)
    // re = (ac + bd)/abs_2()
    // im = (bc - ad)/abs_2()
    auto ac_bd = elwise_mult(b);
    auto d_c = b.el_swapped();
    d_c = d_c ^ rsign_mask;
    auto ad_bc = elwise_mult(d_c);
    auto abs_b = b.abs_2_();
    auto re_im = horizontal_add_permD8(ac_bd, ad_bc);
    auto ret = re_im.elwise_div(abs_b);
    // ret.dump();//
    return ret;
  }

  Vec256<StdComplexFlt> asin() const {
    // asin(x)
    // = -i*ln(iz + sqrt(1 -z^2))
    // = -i*ln((ai - b) + sqrt(1 - (a + bi)*(a + bi)))
    // = -i*ln((-b + ai) + sqrt(1 - (a**2 - b**2) - 2*abi))

#if 1
    auto conj = conj_();
    auto b_a = conj.el_swapped();
    auto ab = conj.elwise_mult(b_a);
    auto im = ab + ab;
    auto val_2 = (*this).elwise_mult(*this);
    auto val_2_swapped = val_2.el_swapped();
    auto re = horizontal_sub_permD8(val_2, val_2_swapped);
    re = Vec256<StdComplexFlt>(one) - re;
    auto root = el_blend<0xAA>(re, im).sqrt();
    auto ln = (b_a + root).log();
    return ln.el_swapped().conj();
#else
    return map(std::asin);
#endif
  }

  Vec256<StdComplexFlt> exp() const { 
    return map(std::exp);
  }



  Vec256<StdComplexFlt> eq(const Vec256<StdComplexFlt>& other) const {
    auto ret = (*this == other);
    return ret & one;
  }
  Vec256<StdComplexFlt> ne(const Vec256<StdComplexFlt>& other) const {
    auto ret = (*this != other);
    return ret & one;
  }

  Vec256<StdComplexFlt> atan2(const Vec256<StdComplexFlt>& b) const {
      AT_ERROR("not supported for complex numbers");
  }
  Vec256<StdComplexFlt> erf() const {
      AT_ERROR("not supported for complex numbers");
  }
  Vec256<StdComplexFlt> erfc() const {
      AT_ERROR("not supported for complex numbers");
  }

  Vec256<StdComplexFlt> log1p() const {
      AT_ERROR("not supported for complex numbers");
  }

  Vec256<StdComplexFlt> expm1() const {
      AT_ERROR("not supported for complex numbers");
  }

  Vec256<StdComplexFlt> operator<(const Vec256<StdComplexFlt>& other) const {
      TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexFlt> operator<=(const Vec256<StdComplexFlt>& other) const {
      TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexFlt> operator>(const Vec256<StdComplexFlt>& other) const {
      TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexFlt> operator>=(const Vec256<StdComplexFlt>& other) const {
      TORCH_CHECK(false, "not supported for complex numbers");
  }

  Vec256<StdComplexFlt> lt(const Vec256<StdComplexFlt>& other) const {
      TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexFlt> le(const Vec256<StdComplexFlt>& other) const {
      TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexFlt> gt(const Vec256<StdComplexFlt>& other) const {
      TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexFlt> ge(const Vec256<StdComplexFlt>& other) const {
      TORCH_CHECK(false, "not supported for complex numbers");
  }


  DEFINE_MEMBER_OP(operator==, StdComplexFlt, vec_cmpeq)
  DEFINE_MEMBER_OP(operator!=, StdComplexFlt, vec_cmpne)

  DEFINE_MEMBER_OP(operator+, StdComplexFlt, vec_add)
  DEFINE_MEMBER_OP(operator-, StdComplexFlt, vec_sub)
  DEFINE_MEMBER_OP(operator&, StdComplexFlt, vec_and)
  DEFINE_MEMBER_OP(operator|, StdComplexFlt, vec_or)
  DEFINE_MEMBER_OP(operator^, StdComplexFlt, vec_xor)
  // elelemtwise helpers
  DEFINE_MEMBER_OP(elwise_mult, StdComplexFlt, vec_mul)
  DEFINE_MEMBER_OP(elwise_div, StdComplexFlt, vec_div)
  DEFINE_MEMBER_OP(elwise_gt, StdComplexFlt, vec_cmpgt)
  DEFINE_MEMBER_OP(elwise_ge, StdComplexFlt, vec_cmpge)
  DEFINE_MEMBER_OP(elwise_lt, StdComplexFlt, vec_cmplt)
  DEFINE_MEMBER_OP(elwise_le, StdComplexFlt, vec_cmple)
};

template <>
Vec256<StdComplexFlt> inline maximum(const Vec256<StdComplexFlt>& a,
                                  const Vec256<StdComplexFlt>& b) {
  auto abs_a = a.abs_2_();
  auto abs_b = b.abs_2_();
  // auto mask = _mm256_cmp_ps(abs_a, abs_b, _CMP_LT_OQ);
  // auto max = _mm256_blendv_ps(a, b, mask);
  auto mask = abs_a.elwise_lt(abs_b);
  auto max = Vec256<StdComplexFlt>::elwise_blendv(a, b, mask);

  return max;
  // Exploit the fact that all-ones is a NaN.
  // auto isnan = _mm256_cmp_ps(abs_a, abs_b, _CMP_UNORD_Q);
  // return _mm256_or_ps(max, isnan);
}

template <>
Vec256<StdComplexFlt> inline minimum(const Vec256<StdComplexFlt>& a,
                                  const Vec256<StdComplexFlt>& b) {
  auto abs_a = a.abs_2_();
  auto abs_b = b.abs_2_();
  // auto mask = _mm256_cmp_ps(abs_a, abs_b, _CMP_GT_OQ);
  // auto min = _mm256_blendv_ps(a, b, mask);
  auto mask = abs_a.elwise_gt(abs_b);
  auto min = Vec256<StdComplexFlt>::elwise_blendv(a, b, mask);
  return min;
  // Exploit the fact that all-ones is a NaN.
  // auto isnan = _mm256_cmp_ps(abs_a, abs_b, _CMP_UNORD_Q);
  // return _mm256_or_ps(min, isnan);
}

template <>
Vec256<StdComplexFlt> inline clamp(const Vec256<StdComplexFlt>& a,
                                const Vec256<StdComplexFlt>& min,
                                const Vec256<StdComplexFlt>& max) {
  auto abs_a = a.abs_2_();
  auto abs_min = min.abs_2_();
  // auto max_mask = _mm256_cmp_ps(abs_a, abs_min, _CMP_LT_OQ);
  auto max_mask = abs_a.elwise_lt(abs_min);
  auto abs_max = max.abs_2_();
  // auto min_mask = _mm256_cmp_ps(abs_a, abs_max, _CMP_GT_OQ);
  auto min_mask = abs_a.elwise_lt(abs_max);
  // return _mm256_blendv_ps(_mm256_blendv_ps(a, min, max_mask), max, min_mask);
  auto f = Vec256<StdComplexFlt>::elwise_blendv(a, min, max_mask);
  return Vec256<StdComplexFlt>::elwise_blendv(f, max, min_mask);
}

template <>
Vec256<StdComplexFlt> inline clamp_min(const Vec256<StdComplexFlt>& a,
                                    const Vec256<StdComplexFlt>& min) {
  auto abs_a = a.abs_2_();
  auto abs_min = min.abs_2_();
  // auto max_mask = _mm256_cmp_ps(abs_a, abs_min, _CMP_LT_OQ);
  // return _mm256_blendv_ps(a, min, max_mask);
  auto max_mask = abs_a.elwise_lt(abs_min);
  return Vec256<StdComplexFlt>::elwise_blendv(a, min, max_mask);
}

template <>
Vec256<StdComplexFlt> inline clamp_max(const Vec256<StdComplexFlt>& a,
                                    const Vec256<StdComplexFlt>& max) {
  auto abs_a = a.abs_2_();
  auto abs_max = max.abs_2_();
  // auto min_mask = _mm256_cmp_ps(abs_a, abs_max, _CMP_GT_OQ);
  // return _mm256_blendv_ps(a, max, min_mask);
  auto min_mask = abs_a.elwise_lt(abs_max);
  return Vec256<StdComplexFlt>::elwise_blendv(a, max, min_mask);
}

}  // namespace
}  // namespace vec256
}  // namespace at
