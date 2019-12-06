#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>

namespace at {
namespace vec256 {
// See Note [Acceptable use of anonymous namespace in header]

namespace {

#if defined(__VSX__)

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
  static constexpr size_t size() {
    return 4;
  }
  Vec256() {}
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

  template <uint64_t mask>
  static c10::guts::enable_if_t<mask == 0, Vec256<double>> __inline_attrs
  blend(const Vec256<double>& a, const Vec256<double>& b) {
    return a;
  }

  template <uint64_t mask>
  static c10::guts::enable_if_t<mask == 3, Vec256<double>> __inline_attrs
  blend(const Vec256<double>& a, const Vec256<double>& b) {
    return Vec256<double>{b._vec0, a._vec1};
  }

  template <uint64_t mask>
  static c10::guts::enable_if_t<(mask & 15) == 15, Vec256<double>> __inline_attrs
  blend(const Vec256<double>& a, const Vec256<double>& b) {
    return b;
  }

  template <uint64_t mask>
  static c10::guts::enable_if_t<(mask > 0 && mask < 3), Vec256<double>>
      __inline_attrs blend(const Vec256<double>& a, const Vec256<double>& b) {
    
    constexpr uint64_t g0 = (mask & 1) * 0xffffffffffffffff;
    constexpr uint64_t g1 = ((mask & 2) >> 1) * 0xffffffffffffffff;
    const __vllb mask_1st = (__vllb){g0, g1};
    return Vec256<double>{(__vd)vec_sel(a._vec0, b._vec0, (__vllb)mask_1st),
                          a._vec1};
  }

  template <uint64_t mask>
  static c10::guts::enable_if_t<(mask > 3) && (mask & 3) == 0, Vec256<double>>
      __inline_attrs blend(const Vec256<double>& a, const Vec256<double>& b) {
    
    constexpr uint64_t g0_2 = ((mask & 4) >> 2) * 0xffffffffffffffff;
    constexpr uint64_t g1_2 = ((mask & 8) >> 3) * 0xffffffffffffffff;

    const __vllb mask_2nd = (__vllb){g0_2, g1_2};
    return Vec256<double>{a._vec0,
                          (__vd)vec_sel(a._vec1, b._vec1, (__vllb)mask_2nd)};
  }

  template <uint64_t mask>
  static c10::guts::enable_if_t<
      (mask > 3) && (mask & 3) != 0 && (mask & 15) != 15,
      Vec256<double>>
      __inline_attrs blend(const Vec256<double>& a, const Vec256<double>& b) {
    
    constexpr uint64_t g0 = (mask & 1) * 0xffffffffffffffff;
    constexpr uint64_t g1 = ((mask & 2) >> 1) * 0xffffffffffffffff;
    constexpr uint64_t g0_2 = ((mask & 4) >> 2) * 0xffffffffffffffff;
    constexpr uint64_t g1_2 = ((mask & 8) >> 3) * 0xffffffffffffffff;

    const __vllb mask_1st = (__vllb){g0, g1};
    const __vllb mask_2nd = (__vllb){g0_2, g1_2};
    return Vec256<double>{(__vd)vec_sel(a._vec0, b._vec0, (__vllb)mask_1st),
                          (__vd)vec_sel(a._vec1, b._vec1, (__vllb)mask_2nd)};
  }

  static Vec256<double> __inline_attrs blendv(
      const Vec256<double>& a,
      const Vec256<double>& b,
      const Vec256<double>& mask) {
    // the mask used here returned by comparision of vec256

    return Vec256<double>{vec_sel(a._vec0, b._vec0, mask._vecb0),
                          vec_sel(a._vec1, b._vec1, mask._vecb1)};
  }
  static Vec256<double> arange(double base = 0., double step = 1.) {
    return Vec256<double>(base, base + step, base + 2 * step, base + 3 * step);
  }

  static Vec256<double> __inline_attrs
  set(const Vec256<double>& a,
      const Vec256<double>& b,
      size_t count = size()) {
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
      loadu(const void* ptr, size_t count = size()) {
      if (count == size()) {
          return Vec256<value_type>{
              vec_vsx_ld(offset0, reinterpret_cast<const value_type*>(ptr)),
                  vec_vsx_ld(offset16, reinterpret_cast<const value_type*>(ptr))};
      }

      __at_align32__ value_type tmp_values[size()];
      std::memcpy(tmp_values, ptr, std::min(count,size()) * sizeof(value_type));

      return Vec256<value_type>{
          vec_vsx_ld(offset0, tmp_values),
              vec_vsx_ld(offset16, tmp_values)};
  }
  void __inline_attrs store(void* ptr, size_t count = size()) const {
      if (count == size()) {
          vec_vsx_st(_vec0, offset0, reinterpret_cast<value_type*>(ptr));
          vec_vsx_st(_vec1, offset16, reinterpret_cast<value_type*>(ptr));
      }
      else if (count > 0) {
          __at_align32__ value_type tmp_values[size()];
          vec_vsx_st(_vec0, offset0, tmp_values);
          vec_vsx_st(_vec1, offset16, tmp_values);
          std::memcpy(ptr, tmp_values, std::min(count,size()) * sizeof(value_type));
      }
  }
  const double& operator[](int idx) const = delete;
  double& operator[](int idx) = delete;

  Vec256<double> map(double (*f)(double)) const {
    return Vec256<double>{f(_vec0[0]), f(_vec0[1]), f(_vec1[0]), f(_vec1[1])};
  }

  Vec256<double> mapbi(double (*f)(double, double), const Vec256<double>& exp)
      const {
    return Vec256<double>{f(_vec0[0], exp._vec0[0]),
                          f(_vec0[1], exp._vec0[1]),
                          f(_vec1[0], exp._vec1[0]),
                          f(_vec1[1], exp._vec1[1])};
  }
  Vec256<double> __inline_attrs abs() const {
    return Vec256<double>{vec_abs(_vec0), vec_abs(_vec1)};
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
    return map(std::expm1);
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
    return Vec256<double>{vec_ceil(_vec0), vec_ceil(_vec1)};
  }
  Vec256<double> __inline_attrs cos() const {
    return map(std::cos);
  }
  Vec256<double> __inline_attrs cosh() const {
    return map(std::cosh);
  }
  Vec256<double> __inline_attrs floor() const {
    return Vec256<double>{vec_floor(_vec0), vec_floor(_vec1)};
  }
  Vec256<double> __inline_attrs neg() const {
    return Vec256<double>{vec_neg(_vec0), vec_neg(_vec1)};
  }
  Vec256<double> __inline_attrs round() const {
    return Vec256<double>{vec_round(_vec0), vec_round(_vec1)};
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
    return Vec256<double>{vec_trunc(_vec0), vec_trunc(_vec1)};
  }

  Vec256<double> __inline_attrs frac() const {
    return Vec256<double>{vec_sub(_vec0, vec_trunc(_vec0)),
                          vec_sub(_vec1, vec_trunc(_vec1))};
  }

  Vec256<double> __inline_attrs sqrt() const {
    return Vec256<double>{vec_sqrt(_vec0), vec_sqrt(_vec1)};
  }
  Vec256<double> __inline_attrs reciprocal() const {
    const __vd one = vec_splats(1.0);
    return Vec256<double>{
        vec_div(one, _vec0), // vec_re(_vec0) is estimated one.
        vec_div(one, _vec1)};
  }
  Vec256<double> __inline_attrs rsqrt() const {
    const __vd one = vec_splats(1.0);
    return Vec256<double>{
        vec_div(one, vec_sqrt(_vec0)), // vec_rsqrt- is estimated rsqrt
        vec_div(one, vec_sqrt(_vec1))};
  }

  Vec256<double> __inline_attrs pow(const Vec256<double>& exp) const {
    return mapbi(std::pow, exp);
  }

  Vec256<double> __inline_attrs operator==(const Vec256<double>& other) const {
    return Vec256<double>{(__vllb)vec_cmpeq(_vec0, other._vec0),
                          (__vllb)vec_cmpeq(_vec1, other._vec1)};
  }

  Vec256<double> __inline_attrs operator!=(const Vec256<double>& other) const {
    return Vec256<double>{(__vllb)vec_cmpne(_vec0, other._vec0),
                          (__vllb)vec_cmpne(_vec1, other._vec1)};
  }

  Vec256<double> __inline_attrs operator<(const Vec256<double>& other) const {
    return Vec256<double>{(__vllb)vec_cmplt(_vec0, other._vec0),
                          (__vllb)vec_cmplt(_vec1, other._vec1)};
  }

  Vec256<double> __inline_attrs operator<=(const Vec256<double>& other) const {
    return Vec256<double>{(__vllb)vec_cmple(_vec0, other._vec0),
                          (__vllb)vec_cmple(_vec1, other._vec1)};
  }

  Vec256<double> __inline_attrs operator>(const Vec256<double>& other) const {
    return Vec256<double>{(__vllb)vec_cmpgt(_vec0, other._vec0),
                          (__vllb)vec_cmpgt(_vec1, other._vec1)};
  }

  Vec256<double> __inline_attrs operator>=(const Vec256<double>& other) const {
    return Vec256<double>{(__vllb)vec_cmpge(_vec0, other._vec0),
                          (__vllb)vec_cmpge(_vec1, other._vec1)};
  }
};

#endif

} // namespace
} // namespace vec256
} // namespace at
