#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>

namespace at
{
    namespace vec256
    {
        // See Note [Acceptable use of anonymous namespace in header]
        namespace
        {

#if defined(__VSX__)

            template <>
            class Vec256<int64_t>
            {
            private:
                union {
                    struct {
                        __vlli _vec0;
                        __vlli _vec1;
                    };
                    struct {
                        __vllb _vecb0;
                        __vllb _vecb1;
                    };

                } __attribute__((__may_alias__));

            public:
                using value_type = int64_t;
                using vec_internal_type = __vlli;
                static constexpr size_t size() {
                    return 4;
                }
                Vec256() {}
                __inline_attrs Vec256(__vlli v1, __vlli v2) : _vec0{ v1 }, _vec1{ v2 } {}
                __inline_attrs Vec256(__vllb v1, __vllb v2) : _vecb0{ v1 }, _vecb1{ v2 } {}
                __inline_attrs Vec256(int64_t scalar)
                    : _vec0{ vec_splats(scalar) }, _vec1{ vec_splats(scalar) } {}
                __inline_attrs Vec256(
                    int64_t scalar1,
                    int64_t scalar2,
                    int64_t scalar3,
                    int64_t scalar4)
                    : _vec0{ __vlli{scalar1, scalar2} }, _vec1{ __vlli{scalar3, scalar4} } {}

                inline __inline_attrs const vec_internal_type& vec0() const {
                    return _vec0;
                }
                inline __inline_attrs const vec_internal_type& vec1() const {
                    return _vec1;
                }


                template <uint64_t mask>
                static c10::guts::enable_if_t<mask == 0, Vec256<int64_t>> __inline_attrs
                    blend(const Vec256<int64_t>& a, const Vec256<int64_t>& b) {
                    return a;
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<mask == 3, Vec256<int64_t>> __inline_attrs
                    blend(const Vec256<int64_t>& a, const Vec256<int64_t>& b) {
                    return Vec256<int64_t> {b._vec0, a._vec1};
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<(mask & 15) == 15, Vec256<int64_t>>
                    __inline_attrs blend(const Vec256<int64_t>& a, const Vec256<int64_t>& b) {
                    return b;
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<(mask > 0 && mask < 3), Vec256<int64_t>>
                    __inline_attrs blend(const Vec256<int64_t>& a, const Vec256<int64_t>& b) {
                    
                    constexpr uint64_t g0 = (mask & 1) * 0xffffffffffffffff;
                    constexpr uint64_t g1 = ((mask & 2) >> 1) * 0xffffffffffffffff;
                    const __vllb mask_1st = (__vllb){
                        g0, g1
                    };
                    return Vec256<int64_t> {
                        (__vlli)vec_sel(a._vec0, b._vec0, (__vllb)mask_1st),
                        a._vec1
                    };
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<(mask > 3) && (mask & 3) == 0, Vec256<int64_t>>
                    __inline_attrs blend(const Vec256<int64_t>& a, const Vec256<int64_t>& b) {
                    
                    constexpr uint64_t g0_2 = ((mask & 4) >> 2) * 0xffffffffffffffff;
                    constexpr uint64_t g1_2 = ((mask & 8) >> 3) * 0xffffffffffffffff;

                    const __vllb mask_2nd = (__vllb){
                        g0_2, g1_2
                    };
                    return Vec256<int64_t> {
                        a._vec0,
                        (__vlli)vec_sel(a._vec1, b._vec1, (__vllb)mask_2nd)
                    };
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<
                    (mask > 3) && (mask & 3) != 0 && (mask & 15) != 15,
                    Vec256<int64_t>>
                    __inline_attrs blend(const Vec256<int64_t>& a, const Vec256<int64_t>& b) {
                    
                    constexpr uint64_t g0 = (mask & 1) * 0xffffffffffffffff;
                    constexpr uint64_t g1 = ((mask & 2) >> 1) * 0xffffffffffffffff;
                    constexpr uint64_t g0_2 = ((mask & 4) >> 2) * 0xffffffffffffffff;
                    constexpr uint64_t g1_2 = ((mask & 8) >> 3) * 0xffffffffffffffff;

                    const __vllb mask_1st = (__vllb){
                        g0, g1
                    };
                    const __vllb mask_2nd = (__vllb){
                        g0_2, g1_2
                    };
                    return Vec256<int64_t> {
                        (__vlli)vec_sel(a._vec0, b._vec0, (__vllb)mask_1st),
                        (__vlli)vec_sel(a._vec1, b._vec1, (__vllb)mask_2nd)
                    };
                }

                static Vec256<int64_t> __inline_attrs blendv(
                    const Vec256<int64_t>& a,
                    const Vec256<int64_t>& b,
                    const Vec256<int64_t>& mask) {
                    // the mask used here returned by comparision of vec256
                    
                    return Vec256<int64_t> {
                        vec_sel(a._vec0, b._vec0, mask._vecb0),
                        vec_sel(a._vec1, b._vec1, mask._vecb1)
                    };
                }
                static Vec256<int64_t> arange(int64_t base = 0., int64_t step = 1.) {
                    return Vec256<int64_t>(base, base + step, base + 2 * step, base + 3 * step);
                }

                static Vec256<int64_t> __inline_attrs
                    set(const Vec256<int64_t>& a,
                        const Vec256<int64_t>& b,
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
                        static_assert(sizeof(double) == sizeof(value_type));
                        const double* dptr = reinterpret_cast<const double*>(ptr);
                        return Vec256<value_type> {
                            //treat it as double load
                            (__vlli)vec_vsx_ld(offset0, dptr),
                                (__vlli)vec_vsx_ld(offset16, dptr)
                        };
                    }

                    double tmp_values[size()] __at_align32__;
                    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

                    return Vec256<value_type> {
                        (__vlli)vec_vsx_ld(offset0, tmp_values),
                            (__vlli)vec_vsx_ld(offset16, tmp_values)
                    };
                }
                void __inline_attrs store(void* ptr, size_t count = size()) const {
                    if (count == size()) {
                        double* dptr = reinterpret_cast<double*>(ptr);
                        vec_vsx_st((__vd)_vec0, offset0, dptr);
                        vec_vsx_st((__vd)_vec1, offset16, dptr);
                    }
                    else if (count > 0) {
                        __at_align32__ double tmp_values[size()];
                        vec_vsx_st((__vd)_vec0, offset0, tmp_values);
                        vec_vsx_st((__vd)_vec1, offset16, tmp_values);
                        std::memcpy(ptr, tmp_values, std::min(count, size()) * sizeof(value_type));
                    }
                }
                const int64_t& operator[](int idx) const = delete;
                int64_t& operator[](int idx) = delete;


                Vec256<int64_t> angle() const {
                    return Vec256<int64_t> {0};
                }
                Vec256<int64_t> real() const {
                    return *this;
                }
                Vec256<int64_t> imag() const {
                    return Vec256<int64_t> {0};
                }
                Vec256<int64_t> conj() const {
                    return *this;
                }

                Vec256<int64_t> __inline_attrs abs() const {
                    return Vec256<int64_t> {vec_abs(_vec0), vec_abs(_vec1)};
                }

                Vec256<int64_t> __inline_attrs neg() const {
                    return Vec256<int64_t> {vec_neg(_vec0), vec_neg(_vec1)};
                }

                Vec256<int64_t> __inline_attrs
                    operator==(const Vec256<int64_t>& other) const {
                    return Vec256<int64_t> {(__vllb)vec_cmpeq(_vec0, other._vec0),
                        (__vllb)vec_cmpeq(_vec1, other._vec1)
                    };
                }

                Vec256<int64_t> __inline_attrs
                    operator!=(const Vec256<int64_t>& other) const {
                    return Vec256<int64_t> {(__vllb)vec_cmpne(_vec0, other._vec0),
                        (__vllb)vec_cmpne(_vec1, other._vec1)
                    };
                }

                Vec256<int64_t> __inline_attrs operator<(const Vec256<int64_t>& other) const {
                    return Vec256<int64_t> {(__vllb)vec_cmplt(_vec0, other._vec0),
                        (__vllb)vec_cmplt(_vec1, other._vec1)
                    };
                }

                Vec256<int64_t> __inline_attrs
                    operator<=(const Vec256<int64_t>& other) const {
                    return Vec256<int64_t> {(__vllb)vec_cmple(_vec0, other._vec0),
                        (__vllb)vec_cmple(_vec1, other._vec1)
                    };
                }

                Vec256<int64_t> __inline_attrs operator>(const Vec256<int64_t>& other) const {
                    return Vec256<int64_t> {(__vllb)vec_cmpgt(_vec0, other._vec0),
                        (__vllb)vec_cmpgt(_vec1, other._vec1)
                    };
                }

                Vec256<int64_t> __inline_attrs
                    operator>=(const Vec256<int64_t>& other) const {
                    return Vec256<int64_t> {(__vllb)vec_cmpge(_vec0, other._vec0),
                        (__vllb)vec_cmpge(_vec1, other._vec1)
                    };
                }
            };

            template <>
            class Vec256<int32_t>
            {
            private:
                union {
                    struct {
                        __vi _vec0;
                        __vi _vec1;
                    };
                    struct {
                        __vib _vecb0;
                        __vib _vecb1;
                    };

                } __attribute__((__may_alias__));

            public:
                using value_type = int32_t;
                using vec_internal_type = __vi;
                static constexpr size_t size() {
                    return 8;
                }
                Vec256() {}
                __inline_attrs Vec256(__vi v1, __vi v2) : _vec0{ v1 }, _vec1{ v2 } {}
                __inline_attrs Vec256(__vib v1, __vib v2) : _vecb0{ v1 }, _vecb1{ v2 } {}
                __inline_attrs Vec256(int32_t scalar)
                    : _vec0{ vec_splats(scalar) }, _vec1{ vec_splats(scalar) } {}
                __inline_attrs Vec256(
                    int32_t scalar1,
                    int32_t scalar2,
                    int32_t scalar3,
                    int32_t scalar4,
                    int32_t scalar5,
                    int32_t scalar6,
                    int32_t scalar7,
                    int32_t scalar8)
                    : _vec0{ __vi{scalar1, scalar2, scalar3, scalar4} },
                    _vec1{ __vi{scalar5, scalar6, scalar7, scalar8} } {}
                inline __inline_attrs const vec_internal_type& vec0() const {
                    return _vec0;
                }
                inline __inline_attrs const vec_internal_type& vec1() const {
                    return _vec1;
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<mask == 0, Vec256<int32_t>> __inline_attrs
                    blend(const Vec256<int32_t>& a, const Vec256<int32_t>& b) {
                    

                    return a;
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<(mask & 255) == 255, Vec256<int32_t>>
                    __inline_attrs blend(const Vec256<int32_t>& a, const Vec256<int32_t>& b) {
                    

                    return b;
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<mask == 15, Vec256<int32_t>> __inline_attrs
                    blend(const Vec256<int32_t>& a, const Vec256<int32_t>& b) {
                    

                    return Vec256<int32_t> {b._vec0, a._vec1};
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<(mask > 0 && mask < 15), Vec256<int32_t>>
                    __inline_attrs blend(const Vec256<int32_t>& a, const Vec256<int32_t>& b) {
                    constexpr uint32_t g0 = (mask & 1) * 0xffffffff;
                    constexpr uint32_t g1 = ((mask & 2) >> 1) * 0xffffffff;
                    constexpr uint32_t g2 = ((mask & 4) >> 2) * 0xffffffff;
                    constexpr uint32_t g3 = ((mask & 8) >> 3) * 0xffffffff;
                    const __vib mask_1st = (__vib){
                        g0, g1, g2, g3
                    };

                    return Vec256<int32_t> {(__vi)vec_sel(a._vec0, b._vec0, (__vib)mask_1st),
                        a._vec1
                    };
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<
                    (mask > 15 && (mask & 255) != 255 && ((mask & 15) == 15)),
                    Vec256<int32_t>>
                    __inline_attrs blend(const Vec256<int32_t>& a, const Vec256<int32_t>& b) {
                    
                    constexpr uint32_t mask2 = (mask & 255) >> 4;
                    constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
                    constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
                    constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
                    constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;

                    const __vib mask_2nd = (__vib){
                        g0_2, g1_2, g2_2, g3_2
                    };
                    // generated masks
                    return Vec256<int32_t> {b._vec0,
                        (__vi)vec_sel(a._vec1, b._vec1, (__vib)mask_2nd)
                    };
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<
                    (mask > 15 && ((mask & 255) != 255) && ((mask & 15) == 0)),
                    Vec256<int32_t>>
                    __inline_attrs blend(const Vec256<int32_t>& a, const Vec256<int32_t>& b) {
                    
                    constexpr uint32_t mask2 = (mask & 255) >> 4;
                    constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
                    constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
                    constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
                    constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;

                    const __vib mask_2nd = (__vib){
                        g0_2, g1_2, g2_2, g3_2
                    };
                    // generated masks
                    return Vec256<int32_t> {a, (__vi)vec_sel(a._vec1, b._vec1, (__vib)mask_2nd)};
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<
                    (mask > 15 && ((mask & 255) != 255) && ((mask & 15) != 0) &&
                    ((mask & 15) != 15)),
                    Vec256<int32_t>>
                    __inline_attrs blend(const Vec256<int32_t>& a, const Vec256<int32_t>& b) {
                    
                    constexpr uint32_t g0 = (mask & 1) * 0xffffffff;
                    constexpr uint32_t g1 = ((mask & 2) >> 1) * 0xffffffff;
                    constexpr uint32_t g2 = ((mask & 4) >> 2) * 0xffffffff;
                    constexpr uint32_t g3 = ((mask & 8) >> 3) * 0xffffffff;
                    constexpr uint32_t mask2 = (mask & 255) >> 4;
                    constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
                    constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
                    constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
                    constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;

                    const __vib mask_1st = (__vib){
                        g0, g1, g2, g3
                    };
                    const __vib mask_2nd = (__vib){
                        g0_2, g1_2, g2_2, g3_2
                    };
                    // generated masks
                    return Vec256<int32_t> {(__vi)vec_sel(a._vec0, b._vec0, (__vib)mask_1st),
                        (__vi)vec_sel(a._vec1, b._vec1, (__vib)mask_2nd)
                    };
                }

                static Vec256<int32_t> __inline_attrs blendv(
                    const Vec256<int32_t>& a,
                    const Vec256<int32_t>& b,
                    const Vec256<int32_t>& mask) {
                    // the mask used here returned by comparision of vec256
                    // assuming this we can use the same mask directly with vec_sel
                    // warning intel style mask will not work properly
                    return Vec256<int32_t> {vec_sel(a._vec0, b._vec0, mask._vecb0),
                        vec_sel(a._vec1, b._vec1, mask._vecb1)
                    };
                }

                static Vec256<int32_t> arange(int32_t base = 0.f, int32_t step = 1.f) {
                    return Vec256<int32_t>(
                        base,
                        base + step,
                        base + 2 * step,
                        base + 3 * step,
                        base + 4 * step,
                        base + 5 * step,
                        base + 6 * step,
                        base + 7 * step);
                }
                static Vec256<int32_t> set(
                    const Vec256<int32_t>& a,
                    const Vec256<int32_t>& b,
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
                static Vec256<value_type> __inline_attrs
                    loadu(const void* ptr, size_t count = size()) {
                    if (count == size()) {
                        return Vec256<value_type> {
                            vec_vsx_ld(offset0, reinterpret_cast<const value_type*>(ptr)),
                                vec_vsx_ld(offset16, reinterpret_cast<const value_type*>(ptr))
                        };
                    }

                    __at_align32__ value_type tmp_values[size()];
                    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

                    return Vec256<value_type> {
                        vec_vsx_ld(offset0, tmp_values),
                            vec_vsx_ld(offset16, tmp_values)
                    };
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
                        std::memcpy(ptr, tmp_values, std::min(count, size()) * sizeof(value_type));
                    }
                }
                const int32_t& operator[](int idx) const = delete;
                int32_t& operator[](int idx) = delete;

                Vec256<int32_t> angle() const {
                    return Vec256<int32_t> {0};
                }
                Vec256<int32_t> real() const {
                    return *this;
                }
                Vec256<int32_t> imag() const {
                    return Vec256<int32_t> {0};
                }
                Vec256<int32_t> conj() const {
                    return *this;
                }

                Vec256<int32_t> __inline_attrs abs() const {
                    return Vec256<int32_t> {vec_abs(_vec0), vec_abs(_vec1)};
                }

                Vec256<int32_t> __inline_attrs neg() const {
                    return Vec256<int32_t> {vec_neg(_vec0), vec_neg(_vec1)};
                }

                Vec256<int32_t> __inline_attrs
                    operator==(const Vec256<int32_t>& other) const {
                    return Vec256<int32_t> {(__vib)vec_cmpeq(_vec0, other._vec0),
                        (__vib)vec_cmpeq(_vec1, other._vec1)
                    };
                }

                Vec256<int32_t> __inline_attrs
                    operator!=(const Vec256<int32_t>& other) const {
                    return Vec256<int32_t> {(__vib)vec_cmpne(_vec0, other._vec0),
                        (__vib)vec_cmpne(_vec1, other._vec1)
                    };
                }

                Vec256<int32_t> __inline_attrs operator<(const Vec256<int32_t>& other) const {
                    return Vec256<int32_t> {(__vib)vec_cmplt(_vec0, other._vec0),
                        (__vib)vec_cmplt(_vec1, other._vec1)
                    };
                }

                Vec256<int32_t> __inline_attrs
                    operator<=(const Vec256<int32_t>& other) const {
                    return Vec256<int32_t> {(__vib)vec_cmple(_vec0, other._vec0),
                        (__vib)vec_cmple(_vec1, other._vec1)
                    };
                }

                Vec256<int32_t> __inline_attrs operator>(const Vec256<int32_t>& other) const {
                    return Vec256<int32_t> {(__vib)vec_cmpgt(_vec0, other._vec0),
                        (__vib)vec_cmpgt(_vec1, other._vec1)
                    };
                }

                Vec256<int32_t> __inline_attrs
                    operator>=(const Vec256<int32_t>& other) const {
                    return Vec256<int32_t> {(__vib)vec_cmpge(_vec0, other._vec0),
                        (__vib)vec_cmpge(_vec1, other._vec1)
                    };
                }
            };

            template <>
            class Vec256<int16_t>
            {
            private:
                union {
                    struct {
                        __vshi _vec0;
                        __vshi _vec1;
                    };
                    struct {
                        __vshb _vecb0;
                        __vshb _vecb1;
                    };

                } __attribute__((__may_alias__));

            public:
                using value_type = int16_t;
                using vec_internal_type = __vshi;
                static constexpr size_t size() {
                    return 16;
                }
                Vec256() {}
                __inline_attrs Vec256(__vshi v1, __vshi v2) : _vec0{ v1 }, _vec1{ v2 } {}
                __inline_attrs Vec256(__vshb v1, __vshb v2) : _vecb0{ v1 }, _vecb1{ v2 } {}
                __inline_attrs Vec256(int16_t scalar)
                    : _vec0{ vec_splats(scalar) }, _vec1{ vec_splats(scalar) } {}

                __inline_attrs Vec256(
                    int16_t scalar1,
                    int16_t scalar2,
                    int16_t scalar3,
                    int16_t scalar4,
                    int16_t scalar5,
                    int16_t scalar6,
                    int16_t scalar7,
                    int16_t scalar8,
                    int16_t scalar9,
                    int16_t scalar10,
                    int16_t scalar11,
                    int16_t scalar12,
                    int16_t scalar13,
                    int16_t scalar14,
                    int16_t scalar15,
                    int16_t scalar16)
                    : _vec0{ __vshi{scalar1,
                                    scalar2,
                                    scalar3,
                                    scalar4,
                                    scalar5,
                                    scalar6,
                                    scalar7,
                                    scalar8
                                   }
                },
                    _vec1{ __vshi{scalar9,
                                  scalar10,
                                  scalar11,
                                  scalar12,
                                  scalar13,
                                  scalar14,
                                  scalar15,
                                  scalar16
                                 }
                } {}
                inline __inline_attrs const vec_internal_type& vec0() const {
                    return _vec0;
                }
                inline __inline_attrs const vec_internal_type& vec1() const {
                    return _vec1;
                }


                template <uint64_t mask>
                static c10::guts::enable_if_t<mask == 0, Vec256<int16_t>> __inline_attrs
                    blend(const Vec256<int16_t>& a, const Vec256<int16_t>& b) {
                    

                    return a;
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<(mask & 65535) == 65535, Vec256<int16_t>>
                    __inline_attrs blend(const Vec256<int16_t>& a, const Vec256<int16_t>& b) {
                    

                    return b;
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<mask == 255, Vec256<int16_t>> __inline_attrs
                    blend(const Vec256<int16_t>& a, const Vec256<int16_t>& b) {
                    

                    return Vec256<int16_t> {b._vec0, a._vec1};
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<(mask > 0 && mask < 255), Vec256<int16_t>>
                    __inline_attrs blend(const Vec256<int16_t>& a, const Vec256<int16_t>& b) {
                    constexpr int16_t g0 = (mask & 1) * 0xffff;
                    constexpr int16_t g1 = ((mask & 2) >> 1) * 0xffff;
                    constexpr int16_t g2 = ((mask & 4) >> 2) * 0xffff;
                    constexpr int16_t g3 = ((mask & 8) >> 3) * 0xffff;
                    constexpr int16_t g4 = ((mask & 16) >> 4) * 0xffff;
                    constexpr int16_t g5 = ((mask & 32) >> 5) * 0xffff;
                    constexpr int16_t g6 = ((mask & 64) >> 6) * 0xffff;
                    constexpr int16_t g7 = ((mask & 128) >> 7) * 0xffff;
                    const __vshi mask_1st = __vshi{ g0, g1, g2, g3, g4, g5, g6, g7 };

                    return Vec256<int16_t> {(__vshi)vec_sel(a._vec0, b._vec0, (__vshb)mask_1st),
                        a._vec1
                    };
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<
                    (mask > 255 && (mask & 65535) != 65535 && ((mask & 255) == 255)),
                    Vec256<int16_t>>
                    __inline_attrs blend(const Vec256<int16_t>& a, const Vec256<int16_t>& b) {
                    
                    constexpr int16_t g0_2 = (mask & 1) * 0xffff;
                    constexpr int16_t g1_2 = ((mask & 2) >> 1) * 0xffff;
                    constexpr int16_t g2_2 = ((mask & 4) >> 2) * 0xffff;
                    constexpr int16_t g3_2 = ((mask & 8) >> 3) * 0xffff;
                    constexpr int16_t g4_2 = ((mask & 16) >> 4) * 0xffff;
                    constexpr int16_t g5_2 = ((mask & 32) >> 5) * 0xffff;
                    constexpr int16_t g6_2 = ((mask & 64) >> 6) * 0xffff;
                    constexpr int16_t g7_2 = ((mask & 128) >> 7) * 0xffff;

                    const __vshi mask_2nd =
                        __vshi{ g0_2, g1_2, g2_2, g3_2, g4_2, g5_2, g6_2, g7_2 };
                    // generated masks
                    return Vec256<int16_t> {b._vec0,
                        (__vshi)vec_sel(a._vec1, b._vec1, (__vshb)mask_2nd)
                    };
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<
                    (mask > 255 && ((mask & 65535) != 65535) && ((mask & 255) == 0)),
                    Vec256<int16_t>>
                    __inline_attrs blend(const Vec256<int16_t>& a, const Vec256<int16_t>& b) {
                    
                    constexpr int16_t mask2 = (mask & 65535) >> 16;
                    constexpr int16_t g0_2 = (mask & 1) * 0xffff;
                    constexpr int16_t g1_2 = ((mask & 2) >> 1) * 0xffff;
                    constexpr int16_t g2_2 = ((mask & 4) >> 2) * 0xffff;
                    constexpr int16_t g3_2 = ((mask & 8) >> 3) * 0xffff;
                    constexpr int16_t g4_2 = ((mask & 16) >> 4) * 0xffff;
                    constexpr int16_t g5_2 = ((mask & 32) >> 5) * 0xffff;
                    constexpr int16_t g6_2 = ((mask & 64) >> 6) * 0xffff;
                    constexpr int16_t g7_2 = ((mask & 128) >> 7) * 0xffff;

                    const __vshi mask_2nd =
                        __vshi{ g0_2, g1_2, g2_2, g3_2, g4_2, g5_2, g6_2, g7_2 };
                    // generated masks
                    return Vec256<int16_t> {a,
                        (__vshi)vec_sel(a._vec1, b._vec1, (__vshb)mask_2nd)
                    };
                }

                template <uint64_t mask>
                static c10::guts::enable_if_t<
                    (mask > 255 && ((mask & 65535) != 65535) && ((mask & 255) != 0) &&
                    ((mask & 255) != 255)),
                    Vec256<int16_t>>
                    __inline_attrs blend(const Vec256<int16_t>& a, const Vec256<int16_t>& b) {
                    
                    constexpr int16_t g0 = (mask & 1) * 0xffff;
                    constexpr int16_t g1 = ((mask & 2) >> 1) * 0xffff;
                    constexpr int16_t g2 = ((mask & 4) >> 2) * 0xffff;
                    constexpr int16_t g3 = ((mask & 8) >> 3) * 0xffff;
                    constexpr int16_t g4 = ((mask & 16) >> 4) * 0xffff;
                    constexpr int16_t g5 = ((mask & 32) >> 5) * 0xffff;
                    constexpr int16_t g6 = ((mask & 64) >> 6) * 0xffff;
                    constexpr int16_t g7 = ((mask & 128) >> 7) * 0xffff;
                    constexpr int16_t mask2 = (mask & 65535) >> 16;
                    constexpr int16_t g0_2 = (mask & 1) * 0xffff;
                    constexpr int16_t g1_2 = ((mask & 2) >> 1) * 0xffff;
                    constexpr int16_t g2_2 = ((mask & 4) >> 2) * 0xffff;
                    constexpr int16_t g3_2 = ((mask & 8) >> 3) * 0xffff;
                    constexpr int16_t g4_2 = ((mask & 16) >> 4) * 0xffff;
                    constexpr int16_t g5_2 = ((mask & 32) >> 5) * 0xffff;
                    constexpr int16_t g6_2 = ((mask & 64) >> 6) * 0xffff;
                    constexpr int16_t g7_2 = ((mask & 128) >> 7) * 0xffff;

                    const __vshi mask_1st = __vshi{ g0, g1, g2, g3, g4, g5, g6, g7 };
                    const __vshi mask_2nd =
                        __vshi{ g0_2, g1_2, g2_2, g3_2, g4_2, g5_2, g6_2, g7_2 };
                    // generated masks
                    return Vec256<int16_t> {(__vshi)vec_sel(a._vec0, b._vec0, (__vshb)mask_1st),
                        (__vshi)vec_sel(a._vec1, b._vec1, (__vshb)mask_2nd)
                    };
                }

                static Vec256<int16_t> __inline_attrs blendv(
                    const Vec256<int16_t>& a,
                    const Vec256<int16_t>& b,
                    const Vec256<int16_t>& mask) {
                    // the mask used here returned by comparision of vec256
                    // assuming this we can use the same mask directly with vec_sel
                    // warning intel style mask will not work properly
                    return Vec256<int16_t> {vec_sel(a._vec0, b._vec0, mask._vecb0),
                        vec_sel(a._vec1, b._vec1, mask._vecb1)
                    };
                }

                static Vec256<int16_t> arange(int16_t base = 0, int16_t step = 1) {
                    return Vec256<int16_t>(
                        base,
                        base + step,
                        base + 2 * step,
                        base + 3 * step,
                        base + 4 * step,
                        base + 5 * step,
                        base + 6 * step,
                        base + 7 * step,
                        base + 8 * step,
                        base + 9 * step,
                        base + 10 * step,
                        base + 11 * step,
                        base + 12 * step,
                        base + 13 * step,
                        base + 14 * step,
                        base + 15 * step);
                }
                static Vec256<int16_t> set(
                    const Vec256<int16_t>& a,
                    const Vec256<int16_t>& b,
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
                    case 4:
                        return blend<15>(a, b);
                    case 5:
                        return blend<31>(a, b);
                    case 6:
                        return blend<63>(a, b);
                    case 7:
                        return blend<127>(a, b);
                    case 8:
                        return blend<255>(a, b);
                    case 9:
                        return blend<511>(a, b);
                    case 10:
                        return blend<1023>(a, b);
                    case 11:
                        return blend<2047>(a, b);
                    case 12:
                        return blend<4095>(a, b);
                    case 13:
                        return blend<8191>(a, b);
                    case 14:
                        return blend<16383>(a, b);
                    case 15:
                        return blend<32767>(a, b);
                    }
                    return b;
                }
                static Vec256<value_type> __inline_attrs
                    loadu(const void* ptr, size_t count = size()) {
                    if (count == size()) {
                        return Vec256<value_type> {
                            vec_vsx_ld(offset0, reinterpret_cast<const value_type*>(ptr)),
                                vec_vsx_ld(offset16, reinterpret_cast<const value_type*>(ptr))
                        };
                    }

                    __at_align32__ value_type tmp_values[size()];
                    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

                    return Vec256<value_type> {
                        vec_vsx_ld(offset0, tmp_values),
                            vec_vsx_ld(offset16, tmp_values)
                    };
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
                        std::memcpy(ptr, tmp_values, std::min(count, size()) * sizeof(value_type));
                    }
                }
                const int16_t& operator[](int idx) const = delete;
                int16_t& operator[](int idx) = delete;

                Vec256<int16_t> angle() const {
                    return Vec256<int16_t> {0};
                }
                Vec256<int16_t> real() const {
                    return *this;
                }
                Vec256<int16_t> imag() const {
                    return Vec256<int16_t> {0};
                }
                Vec256<int16_t> conj() const {
                    return *this;
                }

                Vec256<int16_t> __inline_attrs abs() const {
                    return Vec256<int16_t> {vec_abs(_vec0), vec_abs(_vec1)};
                }

                Vec256<int16_t> __inline_attrs neg() const {
                    return Vec256<int16_t> {vec_neg(_vec0), vec_neg(_vec1)};
                }

                Vec256<int16_t> __inline_attrs
                    operator==(const Vec256<int16_t>& other) const {
                    return Vec256<int16_t> {(__vshb)vec_cmpeq(_vec0, other._vec0),
                        (__vshb)vec_cmpeq(_vec1, other._vec1)
                    };
                }

                Vec256<int16_t> __inline_attrs
                    operator!=(const Vec256<int16_t>& other) const {
                    return Vec256<int16_t> {(__vshb)vec_cmpne(_vec0, other._vec0),
                        (__vshb)vec_cmpne(_vec1, other._vec1)
                    };
                }

                Vec256<int16_t> __inline_attrs operator<(const Vec256<int16_t>& other) const {
                    return Vec256<int16_t> {(__vshb)vec_cmplt(_vec0, other._vec0),
                        (__vshb)vec_cmplt(_vec1, other._vec1)
                    };
                }

                Vec256<int16_t> __inline_attrs
                    operator<=(const Vec256<int16_t>& other) const {
                    return Vec256<int16_t> {(__vshb)vec_cmple(_vec0, other._vec0),
                        (__vshb)vec_cmple(_vec1, other._vec1)
                    };
                }

                Vec256<int16_t> __inline_attrs operator>(const Vec256<int16_t>& other) const {
                    return Vec256<int16_t> {(__vshb)vec_cmpgt(_vec0, other._vec0),
                        (__vshb)vec_cmpgt(_vec1, other._vec1)
                    };
                }

                Vec256<int16_t> __inline_attrs
                    operator>=(const Vec256<int16_t>& other) const {
                    return Vec256<int16_t> {(__vshb)vec_cmpge(_vec0, other._vec0),
                        (__vshb)vec_cmpge(_vec1, other._vec1)
                    };
                }
            };
#endif

        } // namespace
    } // namespace vec256
} // namespace at