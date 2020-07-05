#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>
#include <ATen/cpu/vec256/vsx/vsx_helpers.h>
namespace at {
namespace vec256 {

namespace {


template <>
class Vec256<double> {
 private:
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
  using value_type = double;
  using vec_internal_type = __vd;
  using vec_internal_mask_type = __vllb;
  static constexpr int size() {
    return 4;
  }
  Vec256() {}
  __inline_attrs Vec256(__vd v) : _vec0{v}, _vec1{v} {}
  __inline_attrs Vec256(__vllb vmask) : _vecb0{vmask}, _vecb1{vmask} {}
  __inline_attrs Vec256(__vd v1, __vd v2) : _vec0{v1}, _vec1{v2} {}
  __inline_attrs Vec256(__vllb v1, __vllb v2) : _vecb0{v1}, _vecb1{v2} {}
  __inline_attrs Vec256(double scalar)
      : _vec0{vec_splats(scalar)}, _vec1{vec_splats(scalar)} {}
  __inline_attrs Vec256(
      double scalar1,
      double scalar2,
      double scalar3,
      double scalar4)
      : _vec0{__vd{scalar1, scalar2}}, _vec1{__vd{scalar3, scalar4}} {}
  inline __inline_attrs const vec_internal_type& vec0() const {
    return _vec0;
  }
  inline __inline_attrs const vec_internal_type& vec1() const {
    return _vec1;
  }

  int zero_mask() const {
    auto cmp = (*this == vd_zero);
    return (cmp._vecb0[0] & 1) | (cmp._vecb0[1] & 2) | (cmp._vecb1[0] & 4) |
        (cmp._vecb1[1] & 8);
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 0, Vec256<double>> __inline_attrs
      blend(const Vec256<double>& a, const Vec256<double>& b) {
      return a;
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 1, Vec256<double>> __inline_attrs
      blend(const Vec256<double>& a, const Vec256<double>& b) {
      return b;
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 2, Vec256<double>> __inline_attrs
      blend(const Vec256<double>& a, const Vec256<double>& b) {
      return { b._vec0, a._vec1 };
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 3, Vec256<double>> __inline_attrs
      blend(const Vec256<double>& a, const Vec256<double>& b) {
      return { a._vec0, b._vec1 };
  }
 

  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 4, Vec256<double>> __inline_attrs
      blend(const Vec256<double>& a, const Vec256<double>& b) {
      const __vllb mask_1st = VsxDblMask1(mask);
      return { (__vd)vec_sel(a._vec0, b._vec0, mask_1st), a._vec1 };
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 5, Vec256<double>> __inline_attrs
      blend(const Vec256<double>& a, const Vec256<double>& b) {
      const __vllb mask_1st = VsxDblMask1(mask);
      return { (__vd)vec_sel(a._vec0, b._vec0, mask_1st), b._vec1 };
  }


  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 6,
      Vec256<double>>
      __inline_attrs blend(const Vec256<double>& a, const Vec256<double>& b) {
      const __vllb mask_2nd = VsxDblMask2(mask);
      // generated masks
      return { a._vec0,
          (__vd)vec_sel(a._vec1, b._vec1, mask_2nd) };
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 7,
      Vec256<double>>
      __inline_attrs blend(const Vec256<double>& a, const Vec256<double>& b) {
      const __vllb mask_2nd = VsxDblMask2(mask);
      // generated masks
      return { b._vec0,
          (__vd)vec_sel(a._vec1, b._vec1, mask_2nd) };
  }

  template <int64_t mask>
  static std::enable_if_t<blendChoiceDbl(mask) == 8, Vec256<double>>
      __inline_attrs blend(const Vec256<double>& a, const Vec256<double>& b) {
      const __vllb mask_1st = VsxDblMask1(mask);
      const __vllb mask_2nd = VsxDblMask2(mask);
      return {
          (__vd)vec_sel(a._vec0, b._vec0, mask_1st),
          (__vd)vec_sel(a._vec1, b._vec1, mask_2nd) };
  }


  static Vec256<double> __inline_attrs blendv(
      const Vec256<double>& a,
      const Vec256<double>& b,
      const Vec256<double>& mask) {
    // the mask used here returned by comparision of vec256

    return {
        vec_sel(a._vec0, b._vec0, mask._vecb0),
        vec_sel(a._vec1, b._vec1, mask._vecb1)};
  }
  static Vec256<double> arange(double base = 0., double step = 1.) {
    return Vec256<double>(base, base + step, base + 2 * step, base + 3 * step);
  }

  static Vec256<double> __inline_attrs
  set(const Vec256<double>& a, const Vec256<double>& b, size_t count = size()) {
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
  static Vec256<value_type> __inline_attrs
  loadu(const void* ptr, int count = size()) {
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
  const double& operator[](int idx) const = delete;
  double& operator[](int idx) = delete;

  Vec256<double> map(double (*f)(double)) const {
    return {f(_vec0[0]), f(_vec0[1]), f(_vec1[0]), f(_vec1[1])};
  }

  Vec256<double> mapbi(double (*f)(double, double), const Vec256<double>& exp)
      const {
    return {
        f(_vec0[0], exp._vec0[0]),
        f(_vec0[1], exp._vec0[1]),
        f(_vec1[0], exp._vec1[0]),
        f(_vec1[1], exp._vec1[1])};
  }
  Vec256<double> __inline_attrs abs() const {
    return {vec_abs(_vec0), vec_abs(_vec1)};
  }

  Vec256<double> __inline_attrs acos() const {
    return map(std::acos);
  }
  Vec256<double> __inline_attrs asin() const {
    return map(std::asin);
  }
  Vec256<double> atan() const {
    return map(std::atan);
  }
  Vec256<double> atan2(const Vec256<double>& exp) const {
    return mapbi(std::atan2, exp);
  }
  Vec256<double> erf() const {
    return map(std::erf);
  }
  Vec256<double> erfc() const {
    return map(std::erfc);
  }
  Vec256<double> __inline_attrs exp() const {
    return map(std::exp);
  }
  Vec256<double> expm1() const {
    return  map(std::expm1);
  }

  Vec256<double> lgamma() const {
    return map(std::lgamma);
  }

  Vec256<double> erfinv() const {
    return map(calc_erfinv);
  }

  Vec256<double> angle() const {
    return Vec256<double>{0};
  }
  Vec256<double> real() const {
    return *this;
  }
  Vec256<double> imag() const {
    return Vec256<double>{0};
  }
  Vec256<double> conj() const {
    return *this;
  }

  Vec256<double> __inline_attrs log() const {
    return map(std::log);
  }
  Vec256<double> __inline_attrs log10() const {
    return map(std::log10);
  }
  Vec256<double> __inline_attrs log1p() const {
    return map(std::log1p);
  }
  Vec256<double> __inline_attrs log2() const {
    return map(std::log2);
  }
  Vec256<double> __inline_attrs ceil() const {
    return {vec_ceil(_vec0), vec_ceil(_vec1)};
  }
  Vec256<double> __inline_attrs cos() const {
    return map(std::cos);
  }
  Vec256<double> __inline_attrs cosh() const {
    return map(std::cosh);
  }
  Vec256<double> __inline_attrs floor() const {
    return {vec_floor(_vec0), vec_floor(_vec1)};
  }
  Vec256<double> __inline_attrs neg() const {
    return {vec_neg(_vec0), vec_neg(_vec1)};
  }
  Vec256<double> __inline_attrs round() const {
    return {vec_round(_vec0), vec_round(_vec1)};
  }
  Vec256<double> __inline_attrs sin() const {
    return map(std::sin);
  }
  Vec256<double> __inline_attrs sinh() const {
    return map(std::sinh);
  }
  Vec256<double> __inline_attrs tan() const {
    return map(std::tan);
  }
  Vec256<double> __inline_attrs tanh() const {
    return map(std::tanh);
  }
  Vec256<double> __inline_attrs trunc() const {
    return {vec_trunc(_vec0), vec_trunc(_vec1)};
  }

  Vec256<double> __inline_attrs frac() const {
    return *this - trunc();
  }

  Vec256<double> __inline_attrs sqrt() const {
    return {vec_sqrt(_vec0), vec_sqrt(_vec1)};
  }
  Vec256<double> __inline_attrs reciprocal() const { 
    return {
        vec_div(vd_one, _vec0), // vec_re(_vec0) is estimated one.
        vec_div(vd_one, _vec1)};
  }
  Vec256<double> __inline_attrs rsqrt() const {
    return sqrt().reciprocal();
  }

  Vec256<double> __inline_attrs pow(const Vec256<double>& exp) const {
    return mapbi(std::pow, exp);
  }
  Vec256<double> __inline_attrs fmod(const Vec256<double>& q) const {
    return mapbi(std::fmod, q);
  }

  DEFINE_MEMBER_OP(operator==, double, vec_cmpeq)
  DEFINE_MEMBER_OP(operator!=, double, vec_cmpne)
  DEFINE_MEMBER_OP(operator<, double, vec_cmplt)
  DEFINE_MEMBER_OP(operator<=, double, vec_cmple)
  DEFINE_MEMBER_OP(operator>, double, vec_cmpgt)
  DEFINE_MEMBER_OP(operator>=, double, vec_cmpge)
  DEFINE_MEMBER_OP_AND_ONE(eq, double, vec_cmpeq)
  DEFINE_MEMBER_OP_AND_ONE(ne, double, vec_cmpne)
  DEFINE_MEMBER_OP_AND_ONE(lt, double, vec_cmplt)
  DEFINE_MEMBER_OP_AND_ONE(le, double, vec_cmple)
  DEFINE_MEMBER_OP_AND_ONE(gt, double, vec_cmpgt)
  DEFINE_MEMBER_OP_AND_ONE(ge, double, vec_cmpge)
  DEFINE_MEMBER_OP(operator+, double, vec_add)
  DEFINE_MEMBER_OP(operator-, double, vec_sub)
  DEFINE_MEMBER_OP(operator*, double, vec_mul)
  DEFINE_MEMBER_OP(operator/, double, vec_div)
  DEFINE_MEMBER_OP(maximum, double, vec_max)
  DEFINE_MEMBER_OP(minimum, double, vec_min)
  DEFINE_MEMBER_OP(operator&, double, vec_and)
  DEFINE_MEMBER_OP(operator|, double, vec_or)
  DEFINE_MEMBER_OP(operator^, double, vec_xor)
  DEFINE_MEMBER_TERNARY_OP(madd, double, vec_madd)
};
template <>
Vec256<double> inline maximum(
    const Vec256<double>& a,
    const Vec256<double>& b) {
  return a.maximum(b);
}

template <>
Vec256<double> inline minimum(
    const Vec256<double>& a,
    const Vec256<double>& b) {
  return a.minimum(b);
}
} // namespace
} // namespace vec256
} // namespace at
