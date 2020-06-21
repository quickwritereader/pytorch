#pragma once
#include <c10/util/complex_type.h>
#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>
#include <ATen/cpu/vec256/vsx/vsx_helpers.h>

namespace at {
namespace vec256 {
// See Note [Acceptable use of anonymous namespace in header]
namespace {
using StdComplexDbl = std::complex<double>;

template <>
class Vec256<StdComplexDbl> {
  union {
    struct {
      __vd _vec0;
      __vd _vec1;
    };
    struct {
      __vllb _vecb0;
      __vllb _vecb1;
    };

  } __attribute__((__may_alias__));

 public:
  using value_type = StdComplexDbl;
  using vec_internal_type = __vd;
  using vec_internal_mask_type = __vllb;
  static constexpr int size() { return 2; }
  Vec256() {}
  __inline_attrs Vec256(__vd v) : _vec0{v}, _vec1{v} {}
  __inline_attrs Vec256(__vllb vmask) : _vecb0{vmask}, _vecb1{vmask} {}
  __inline_attrs Vec256(__vd v1, __vd v2) : _vec0{v1}, _vec1{v2} {}
  __inline_attrs Vec256(__vllb v1, __vllb v2) : _vecb0{v1}, _vecb1{v2} {}

  Vec256(StdComplexDbl val) {
    double real_value = val.real();
    double imag_value = val.imag();
    _vec0 = __vd{real_value, imag_value};
    _vec1 = __vd{real_value, imag_value};
  }
  Vec256(StdComplexDbl val1, StdComplexDbl val2) {
    _vec0 = __vd{val1.real(), val1.imag()};
    _vec1 = __vd{val2.real(), val2.imag()};
  }

  inline __inline_attrs const vec_internal_type& vec0() const { return _vec0; }
  inline __inline_attrs const vec_internal_type& vec1() const { return _vec1; }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceStdComplexDbl(mask) == 0, Vec256<StdComplexDbl>>
      __inline_attrs blend(const Vec256<StdComplexDbl>& a,
                           const Vec256<StdComplexDbl>& b) {
    return a;
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceStdComplexDbl(mask) == 1, Vec256<StdComplexDbl>>
      __inline_attrs blend(const Vec256<StdComplexDbl>& a,
                           const Vec256<StdComplexDbl>& b) {
    return b;
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceStdComplexDbl(mask) == 2, Vec256<StdComplexDbl>>
      __inline_attrs blend(const Vec256<StdComplexDbl>& a,
                           const Vec256<StdComplexDbl>& b) {
    return {b._vec0, a._vec1};
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceStdComplexDbl(mask) == 3, Vec256<StdComplexDbl>>
      __inline_attrs blend(const Vec256<StdComplexDbl>& a,
                           const Vec256<StdComplexDbl>& b) {
    return {a._vec0, b._vec1};
  }

  template <int64_t mask>
      static Vec256<StdComplexDbl>__inline_attrs el_blend(const Vec256<StdComplexDbl>& a,
                                         const Vec256<StdComplexDbl>& b) {
    const __vllb mask_1st = VsxDblMask1(mask);
    const __vllb mask_2nd = VsxDblMask2(mask);
    return {(__vd)vec_sel(a._vec0, b._vec0, mask_1st),
            (__vd)vec_sel(a._vec1, b._vec1, mask_2nd)};
  }

  static Vec256<StdComplexDbl> blendv(const Vec256<StdComplexDbl>& a,
                                   const Vec256<StdComplexDbl>& b,
                                   const Vec256<StdComplexDbl>& mask) {
    // convert std::complex<V> index mask to V index mask: xy -> xxyy
    auto mask_complex =
        Vec256<StdComplexDbl>(vec_splat(mask._vec0, 0), vec_splat(mask._vec1, 0));
    return {vec_sel(a._vec0, b._vec0, mask_complex._vecb0),
            vec_sel(a._vec1, b._vec1, mask_complex._vecb1)};
  }

  static Vec256<StdComplexDbl> __inline_attrs
  elwise_blendv(const Vec256<StdComplexDbl>& a, const Vec256<StdComplexDbl>& b,
            const Vec256<StdComplexDbl>& mask) {
    return {vec_sel(a._vec0, b._vec0, mask._vecb0),
            vec_sel(a._vec1, b._vec1, mask._vecb1)};
  }
  template <typename step_t>
  static Vec256<StdComplexDbl> arange(StdComplexDbl base = 0.,
                                   step_t step = static_cast<step_t>(1)) {
    return Vec256<StdComplexDbl>(base, base + step);
  }
  static Vec256<StdComplexDbl> set(const Vec256<StdComplexDbl>& a,
                                const Vec256<StdComplexDbl>& b,
                                int64_t count = size()) {
    switch (count) {
      case 0:
        return a;
      case 1:
        return blend<1>(a, b);
    }
    return b;
  }

  static Vec256<value_type> __inline_attrs loadu(const void* ptr,
                                                 int count = size()) {
    if (count == size()) {
      return {vec_vsx_ld(offset0, reinterpret_cast<const double*>(ptr)),
              vec_vsx_ld(offset16, reinterpret_cast<const double*>(ptr))};
    }

    __at_align32__ value_type tmp_values[size()];
    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

    return {vec_vsx_ld(offset0, reinterpret_cast<const double*>(tmp_values)),
            vec_vsx_ld(offset16, reinterpret_cast<const double*>(tmp_values))};
  }
  void __inline_attrs store(void* ptr, int count = size()) const {
    if (count == size()) {
      vec_vsx_st(_vec0, offset0, reinterpret_cast<double*>(ptr));
      vec_vsx_st(_vec1, offset16, reinterpret_cast<double*>(ptr));
    } else if (count > 0) {
      __at_align32__ value_type tmp_values[size()];
      vec_vsx_st(_vec0, offset0, reinterpret_cast<double*>(tmp_values));
      vec_vsx_st(_vec1, offset16, reinterpret_cast<double*>(tmp_values));
      std::memcpy(ptr, tmp_values,
                  std::min(count, size()) * sizeof(value_type));
    }
  }

  const StdComplexDbl& operator[](int idx) const = delete;
  StdComplexDbl& operator[](int idx) = delete;
  Vec256<StdComplexDbl> map(StdComplexDbl (*f)(const StdComplexDbl&)) const {
    __at_align32__ StdComplexDbl tmp[size()];
    store(tmp);
    for (int i = 0; i < size(); i++) {
      tmp[i] = f(tmp[i]);
    }
    return loadu(tmp);
  }

  Vec256<StdComplexDbl> el_swapped() const {
    __vd v0 = vec_xxpermdi(_vec0, _vec0, 2);
    __vd v1 = vec_xxpermdi(_vec1, _vec1, 2);
    return {v0, v1};
  }

  Vec256<StdComplexDbl> get_evens() const {
    __vd v0 = vec_splat(_vec0, 0);
    __vd v1 = vec_splat(_vec1, 0);
    return {v0, v1};
  }

  static Vec256<StdComplexDbl> get_evens(Vec256<StdComplexDbl>& first,
                                      Vec256<StdComplexDbl>& second) {
    // as mergee phased in , we can use vec_perm with mask
    return {vec_mergeh(first._vec0, second._vec0),
            vec_mergeh(first._vec1, second._vec1)};
  }

  Vec256<StdComplexDbl> abs_2_() const {
    auto a = (*this).elwise_mult(*this);
    auto permuted = a.el_swapped();
    a = a + permuted;
    return a;
  }

  Vec256<StdComplexDbl> abs_() const {
    auto ret = abs_2_();
    return ret.elwise_sqrt();
  }

  Vec256<StdComplexDbl> abs() const { return abs_() & vd_real_mask; }

  Vec256<StdComplexDbl> angle_() const {
    // angle = atan2(b/a)
    // auto b_a = _mm256_permute_pd(values, 0x05);     // b        a
    // return Sleef_atan2d4_u10(values, b_a);          // 90-angle angle
    auto ret = el_swapped();
    for (int i = 0; i < 2; i++) {
      ret._vec0[i] = std::atan2(_vec0[i], ret._vec0[i]);
      ret._vec1[i] = std::atan2(_vec1[i], ret._vec0[i]);
    }
    return ret;
  }

  Vec256<StdComplexDbl> angle() const {
    auto a = angle_().el_swapped();
    return a & vd_real_mask;
  }

  Vec256<StdComplexDbl> real_() const { return *this & vd_real_mask; }
  Vec256<StdComplexDbl> real() const { return *this & vd_real_mask; }
  Vec256<StdComplexDbl> imag_() const { return *this & vd_imag_mask; }
  Vec256<StdComplexDbl> imag() const { return imag_().el_swapped(); }

  Vec256<StdComplexDbl> conj_() const { return *this ^ vd_isign_mask; }
  Vec256<StdComplexDbl> conj() const { return *this ^ vd_isign_mask; }

  Vec256<StdComplexDbl> log() const {
    // Most trigonomic ops use the log() op to improve complex number
    // performance.
    return map(std::log);
  }

  Vec256<StdComplexDbl> log2() const {
    // log2eB_inv
    auto ret = log();
    return ret.elwise_mult(vd_log2e_inv);
  }
  Vec256<StdComplexDbl> log10() const {
    auto ret = log();
    return ret.elwise_mult(vd_log10e_inv);
  }


  Vec256<StdComplexDbl> asin() const {
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
    auto re = horizontal_sub(val_2, val_2_swapped);
    re = Vec256<StdComplexDbl>(vd_one) - re;
    auto root = el_blend<0x0A>(re, im).sqrt();
    auto ln = (b_a + root).log();
    return ln.el_swapped().conj();
#else
    return map(std::asin);
#endif
  }

  Vec256<StdComplexDbl> acos() const {
    // acos(x) = pi/2 - asin(x)
    return Vec256(vd_pi_2) - asin();
  }

  Vec256<StdComplexDbl> atan() const {
    // atan(x) = i/2 * ln((i + z)/(i - z))
    auto ione = Vec256(vd_imag_one);
    auto sum = ione + *this;
    auto sub = ione - *this;
    auto ln = (sum / sub).log();  // ln((i + z)/(i - z))
    return ln * vd_imag_half;     // i/2*ln()
  }

  Vec256<StdComplexDbl> sin() const { return map(std::sin); }
  Vec256<StdComplexDbl> sinh() const { return map(std::sinh); }
  Vec256<StdComplexDbl> cos() const { return map(std::cos); }
  Vec256<StdComplexDbl> cosh() const { return map(std::cosh); }

  Vec256<StdComplexDbl> tan() const { return map(std::tan); }

  Vec256<StdComplexDbl> ceil() const { return {vec_ceil(_vec0), vec_ceil(_vec1)}; }
  Vec256<StdComplexDbl> floor() const {
    return {vec_floor(_vec0), vec_floor(_vec1)};
  }
  Vec256<StdComplexDbl> neg() const {
    auto z = Vec256<StdComplexDbl>(vd_zero);
    return z - *this;
  }
  Vec256<StdComplexDbl> round() const {
    return {vec_round(_vec0), vec_round(_vec1)};
  }

  Vec256<StdComplexDbl> trunc() const {
    return {vec_trunc(_vec0), vec_trunc(_vec1)};
  }

  Vec256<StdComplexDbl> elwise_sqrt() const {
    return {vec_sqrt(_vec0), vec_sqrt(_vec1)};
  }

  void dump() const {
    std::cout << _vec0[0] << "," << _vec0[1] << ",";
    std::cout << _vec1[0] << "," << _vec1[1] << std::endl;
  }

  Vec256<StdComplexDbl> sqrt() const {
    //   sqrt(a + bi)
    // = sqrt(2)/2 * [sqrt(sqrt(a**2 + b**2) + a) + sgn(b)*sqrt(sqrt(a**2 +
    // b**2) - a)i] = sqrt(2)/2 * [sqrt(abs() + a) + sgn(b)*sqrt(abs() - a)i]

    auto sign = *this & vd_isign_mask;
    auto factor = sign | vd_sqrt2_2;
    auto a_a = get_evens();
    // a_a.dump();
    a_a = a_a ^ vd_isign_mask;  // a -a
    auto res_re_im =
        (abs_() + a_a).elwise_sqrt();  // sqrt(abs + a) sqrt(abs - a)
    return factor.elwise_mult(res_re_im);
  }

  Vec256<StdComplexDbl> reciprocal() const {
    // re + im*i = (a + bi)  / (c + di)
    // re = (ac + bd)/abs_2() = c/abs_2()
    // im = (bc - ad)/abs_2() = d/abs_2()
    auto c_d = *this ^ vd_isign_mask;  // c       -d
    auto abs = abs_2_();
    return c_d.elwise_div(abs);
  }

  Vec256<StdComplexDbl> rsqrt() const { return sqrt().reciprocal(); }

  static Vec256<StdComplexDbl> horizontal_add(Vec256<StdComplexDbl>& first,
                                           Vec256<StdComplexDbl>& second) {
    auto first_perm = first.el_swapped();    // 2perm
    auto second_perm = second.el_swapped();  // 2perm
    // summ
    auto first_ret = first + first_perm;     // 2add
    auto second_ret = second + second_perm;  // 2 add
    // now lets choose evens
    return get_evens(first_ret, second_ret);  // 2 mergee's
  }

  static Vec256<StdComplexDbl> horizontal_sub(Vec256<StdComplexDbl>& first,
                                           Vec256<StdComplexDbl>& second) {
    // we will simulate it differently with 6 instructions total
    // lets permute second so that we can add it getting horizontall sums
    auto first_perm = first.el_swapped();    // 2perm
    auto second_perm = second.el_swapped();  // 2perm
    // summ
    auto first_ret = first - first_perm;     // 2sub
    auto second_ret = second - second_perm;  // 2 sub
    // now lets choose evens
    return get_evens(first_ret, second_ret);  // 2 mergee's
  }

  Vec256<StdComplexDbl> inline operator*(const Vec256<StdComplexDbl>& b) const {
    //(a + bi)  * (c + di) = (ac - bd) + (ad + bc)i
    auto ac_bd = elwise_mult(b);
    auto d_c = b.el_swapped();
    d_c = d_c ^ vd_isign_mask;
    auto ad_bc = elwise_mult(d_c);
    auto ret = horizontal_sub(ac_bd, ad_bc);
    return ret;
  }

  Vec256<StdComplexDbl> inline operator/(const Vec256<StdComplexDbl>& b) const {
    // re + im*i = (a + bi)  / (c + di)
    // re = (ac + bd)/abs_2()
    // im = (bc - ad)/abs_2()
    auto ac_bd = elwise_mult(b);
    auto d_c = b.el_swapped();
    d_c = d_c ^ vd_rsign_mask;
    auto ad_bc = elwise_mult(d_c);
    auto abs_b = b.abs_2_();
    auto re_im = horizontal_add(ac_bd, ad_bc);
    auto ret = re_im.elwise_div(abs_b);
    // ret.dump();//
    return ret;
  }

  Vec256<StdComplexDbl> exp() const { return map(std::exp); }

  Vec256<StdComplexDbl> pow(const Vec256<StdComplexDbl>& exp) const {
    __at_align32__ StdComplexDbl x_tmp[size()];
    __at_align32__ StdComplexDbl y_tmp[size()];
    store(x_tmp);
    exp.store(y_tmp);
    for (int i = 0; i < size(); i++) {
      x_tmp[i] = std::pow(x_tmp[i], y_tmp[i]);
    }
    return loadu(x_tmp);
  }

  Vec256<StdComplexDbl> log1p() const {
      AT_ERROR("not supported for complex numbers");
  }

  Vec256<StdComplexDbl> atan2(const Vec256<StdComplexDbl>& b) const {
      AT_ERROR("not supported for complex numbers");
  }
  Vec256<StdComplexDbl> erf() const {
      AT_ERROR("not supported for complex numbers");
  }
  Vec256<StdComplexDbl> erfc() const {
      AT_ERROR("not supported for complex numbers");
  }

  Vec256<StdComplexDbl> expm1() const {
      AT_ERROR("not supported for complex numbers");
  }

  Vec256<StdComplexDbl> operator<(const Vec256<StdComplexDbl>& other) const {
    TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexDbl> operator<=(const Vec256<StdComplexDbl>& other) const {
    TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexDbl> operator>(const Vec256<StdComplexDbl>& other) const {
    TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexDbl> operator>=(const Vec256<StdComplexDbl>& other) const {
    TORCH_CHECK(false, "not supported for complex numbers");
  }

  Vec256<StdComplexDbl> eq(const Vec256<StdComplexDbl>& other) const {
    auto ret = (*this == other);
    return ret & vd_one;
  }
  Vec256<StdComplexDbl> ne(const Vec256<StdComplexDbl>& other) const {
    auto ret = (*this != other);
    return ret & vd_one;
  }

  Vec256<StdComplexDbl> lt(const Vec256<StdComplexDbl>& other) const {
    TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexDbl> le(const Vec256<StdComplexDbl>& other) const {
    TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexDbl> gt(const Vec256<StdComplexDbl>& other) const {
    TORCH_CHECK(false, "not supported for complex numbers");
  }
  Vec256<StdComplexDbl> ge(const Vec256<StdComplexDbl>& other) const {
    TORCH_CHECK(false, "not supported for complex numbers");
  }

  DEFINE_MEMBER_OP(operator==, StdComplexDbl, vec_cmpeq)
  DEFINE_MEMBER_OP(operator!=, StdComplexDbl, vec_cmpne)

  DEFINE_MEMBER_OP(operator+, StdComplexDbl, vec_add)
  DEFINE_MEMBER_OP(operator-, StdComplexDbl, vec_sub)
  DEFINE_MEMBER_OP(operator&, StdComplexDbl, vec_and)
  DEFINE_MEMBER_OP(operator|, StdComplexDbl, vec_or)
  DEFINE_MEMBER_OP(operator^, StdComplexDbl, vec_xor)
  // elelemtwise helpers
  DEFINE_MEMBER_OP(elwise_mult, StdComplexDbl, vec_mul)
  DEFINE_MEMBER_OP(elwise_div, StdComplexDbl, vec_div)
  DEFINE_MEMBER_OP(elwise_gt, StdComplexDbl, vec_cmpgt)
  DEFINE_MEMBER_OP(elwise_ge, StdComplexDbl, vec_cmpge)
  DEFINE_MEMBER_OP(elwise_lt, StdComplexDbl, vec_cmplt)
  DEFINE_MEMBER_OP(elwise_le, StdComplexDbl, vec_cmple)
};

template <>
Vec256<StdComplexDbl> inline maximum(const Vec256<StdComplexDbl>& a,
                                  const Vec256<StdComplexDbl>& b) {
  auto abs_a = a.abs_2_();
  auto abs_b = b.abs_2_();
  // auto mask = _mm256_cmp_ps(abs_a, abs_b, _CMP_LT_OQ);
  // auto max = _mm256_blendv_ps(a, b, mask);
  auto mask = abs_a.elwise_lt(abs_b);
  auto max = Vec256<StdComplexDbl>::elwise_blendv(a, b, mask);

  return max;
  // Exploit the fact that all-ones is a NaN.
  // auto isnan = _mm256_cmp_ps(abs_a, abs_b, _CMP_UNORD_Q);
  // return _mm256_or_ps(max, isnan);
}

template <>
Vec256<StdComplexDbl> inline minimum(const Vec256<StdComplexDbl>& a,
                                  const Vec256<StdComplexDbl>& b) {
  auto abs_a = a.abs_2_();
  auto abs_b = b.abs_2_();
  // auto mask = _mm256_cmp_ps(abs_a, abs_b, _CMP_GT_OQ);
  // auto min = _mm256_blendv_ps(a, b, mask);
  auto mask = abs_a.elwise_gt(abs_b);
  auto min = Vec256<StdComplexDbl>::elwise_blendv(a, b, mask);
  return min;
  // Exploit the fact that all-ones is a NaN.
  // auto isnan = _mm256_cmp_ps(abs_a, abs_b, _CMP_UNORD_Q);
  // return _mm256_or_ps(min, isnan);
}

template <>
Vec256<StdComplexDbl> inline clamp(const Vec256<StdComplexDbl>& a,
                                const Vec256<StdComplexDbl>& min,
                                const Vec256<StdComplexDbl>& max) {
  auto abs_a = a.abs_2_();
  auto abs_min = min.abs_2_();
  // auto max_mask = _mm256_cmp_ps(abs_a, abs_min, _CMP_LT_OQ);
  auto max_mask = abs_a.elwise_lt(abs_min);
  auto abs_max = max.abs_2_();
  // auto min_mask = _mm256_cmp_ps(abs_a, abs_max, _CMP_GT_OQ);
  auto min_mask = abs_a.elwise_lt(abs_max);
  // return _mm256_blendv_ps(_mm256_blendv_ps(a, min, max_mask), max, min_mask);
  auto f = Vec256<StdComplexDbl>::elwise_blendv(a, min, max_mask);
  return Vec256<StdComplexDbl>::elwise_blendv(f, max, min_mask);
}

template <>
Vec256<StdComplexDbl> inline clamp_min(const Vec256<StdComplexDbl>& a,
                                    const Vec256<StdComplexDbl>& min) {
  auto abs_a = a.abs_2_();
  auto abs_min = min.abs_2_();
  // auto max_mask = _mm256_cmp_ps(abs_a, abs_min, _CMP_LT_OQ);
  // return _mm256_blendv_ps(a, min, max_mask);
  auto max_mask = abs_a.elwise_lt(abs_min);
  return Vec256<StdComplexDbl>::elwise_blendv(a, min, max_mask);
}

template <>
Vec256<StdComplexDbl> inline clamp_max(const Vec256<StdComplexDbl>& a,
                                    const Vec256<StdComplexDbl>& max) {
  auto abs_a = a.abs_2_();
  auto abs_max = max.abs_2_();
  // auto min_mask = _mm256_cmp_ps(abs_a, abs_max, _CMP_GT_OQ);
  // return _mm256_blendv_ps(a, max, min_mask);
  auto min_mask = abs_a.elwise_lt(abs_max);
  return Vec256<StdComplexDbl>::elwise_blendv(a, max, min_mask);
}

}  // namespace
}  // namespace vec256
}  // namespace at
