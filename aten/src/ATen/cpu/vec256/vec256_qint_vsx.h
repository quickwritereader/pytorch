#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>
#include <c10/util/qint32.h>
#include <c10/util/qint8.h>
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
//  Vec256<qint8> -> 4x Vec256<float>
//  Vec256<quint8> -> 4x Vec256<float>
//  Vec256<qint32> -> 1x Vec256<float>
//
// The size of the returned float vector is specified by the special
// constexpr function float_num_vecs. The type of the value returned
// from dequantize (and expected as an argument to quantize) is
// specified by float_vec_return_type.
//
// When writing kernels with these vectors, it is expected that floating-
// point operations will be carried out in a loop over Vec256<T>::float_num_vecs
// iterations.

namespace at
{
    namespace vec256
    {
        namespace
        {
#if defined(__VSX__)

            template <>
            struct Vec256<c10::qint8> {
            private:
                union {
                    struct {
                        __vchari _vec0;
                        __vchari _vec1;
                    };
                    struct {
                        __vcharb _vecb0;
                        __vcharb _vecb1;
                    };

                } __attribute__((__may_alias__));

                Vec256() {}

            public:
                static constexpr size_t size() {
                    return 32;
                }

                static constexpr size_t float_num_vecs() {
                    return 4;
                }

                using float_vec_return_type = std::array<Vec256<float>, 4>;
                using value_type = typename c10::qint8::underlying;
                using vec_internal_type = __vchari;
                // Broadcast constructor
                __inline_attrs Vec256(const c10::qint8& val)
                    : _vec0{ vec_splats(val.val_) }, _vec1{ vec_splats(val.val_) } {}

                __inline_attrs Vec256(const Vec256<c10::qint8>& other)
                    : _vec0{ other._vec0 }, _vec1(other._vec1) {}

                __inline_attrs Vec256(__vchari v1, __vchari v2) : _vec0{ v1 }, _vec1{ v2 } {}
                __inline_attrs Vec256(__vcharb v1, __vcharb v2) : _vecb0{ v1 }, _vecb1{ v2 } {}

                inline __inline_attrs const vec_internal_type& vec0() const {
                    return _vec0;
                }
                inline __inline_attrs const vec_internal_type& vec1() const {
                    return _vec1;
                }


                static __inline_attrs Vec256<c10::qint8> loadu(const void* ptr, size_t count = size()) {
                    if (count == size()) {
                        return Vec256<c10::qint8> {
                            vec_vsx_ld(offset0, reinterpret_cast<const __vchari*>(ptr)),
                                vec_vsx_ld(offset16, reinterpret_cast<const __vchari*>(ptr))
                        };
                    }
                    __at_align32__ value_type tmp_values[size()];
                    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

                    return Vec256<c10::qint8> {
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


            public:
                float_vec_return_type __inline_attrs dequantize(
                    Vec256<float> scale,
                    Vec256<float> zero_point,
                    Vec256<float> scale_zp_premul) const {
                    __vshi vecshi0 = vec_unpackh(_vec0);
                    __vshi vecshi1 = vec_unpackl(_vec0);

                    __vshi vecshi2 = vec_unpackh(_vec1);
                    __vshi vecshi3 = vec_unpackl(_vec1);

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
                            vec_madd(scale_vec1, vecf1_0, scale_zp_premul1)
                        },
                        Vec256<float>{
                            vec_madd(scale_vec0, vecf0_1, scale_zp_premul0),
                            vec_madd(scale_vec1, vecf1_1, scale_zp_premul1)
                        },
                        Vec256<float>{
                            vec_madd(scale_vec0, vecf0_2, scale_zp_premul0),
                            vec_madd(scale_vec1, vecf1_2, scale_zp_premul1)
                        },
                        Vec256<float>{
                            vec_madd(scale_vec0, vecf0_3, scale_zp_premul0),
                            vec_madd(scale_vec1, vecf1_3, scale_zp_premul1)
                        }
                    };
                }

                static Vec256<c10::qint8> quantize(
                    const float_vec_return_type& rhs,
                    float scale,
                    int32_t zero_point,
                    float inverse_scale) {
                    // constexpr int32_t min_val = std::numeric_limits<value_type>::min();
                    // constexpr int32_t max_val = std::numeric_limits<value_type>::max();

                    __vf inverse_scale_v = vec_splats(inverse_scale);
                    __vf vec_zero_point = vec_splats((float)zero_point);
                    // __vi vmin = vec_splats(min_val);
                    // __vi vmax = vec_splats(max_val);

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

                    vecf0 = vec_madd(vecf0, inverse_scale_v, vec_zero_point);
                    vecf1 = vec_madd(vecf1, inverse_scale_v, vec_zero_point);
                    vecf2 = vec_madd(vecf2, inverse_scale_v, vec_zero_point);
                    vecf3 = vec_madd(vecf3, inverse_scale_v, vec_zero_point);

                    vecf4 = vec_madd(vecf4, inverse_scale_v, vec_zero_point);
                    vecf5 = vec_madd(vecf5, inverse_scale_v, vec_zero_point);
                    vecf6 = vec_madd(vecf6, inverse_scale_v, vec_zero_point);
                    vecf7 = vec_madd(vecf7, inverse_scale_v, vec_zero_point);

                    __vi veci0 = vec_signed(vecf0);
                    __vi veci1 = vec_signed(vecf1);
                    __vi veci2 = vec_signed(vecf2);
                    __vi veci3 = vec_signed(vecf3);

                    __vi veci4 = vec_signed(vecf4);
                    __vi veci5 = vec_signed(vecf5);
                    __vi veci6 = vec_signed(vecf6);
                    __vi veci7 = vec_signed(vecf7);

                    // veci0 = vec_min(vmax, vec_max( vmin, vecf0)) ;
                    // veci1 = vec_min(vmax, vec_max( vmin, vecf1)) ;
                    // veci2 = vec_min(vmax, vec_max( vmin, vecf2)) ;
                    // veci3 = vec_min(vmax, vec_max( vmin, vecf3)) ;

                    // veci4 = vec_min(vmax, vec_max( vmin, vecf4)) ;
                    // veci5 = vec_min(vmax, vec_max( vmin, vecf5)) ;
                    // veci6 = vec_min(vmax, vec_max( vmin, vecf6)) ;
                    // veci7 = vec_min(vmax, vec_max( vmin, vecf7)) ;
                    // vec_packs CLAMP already
                    __vshi vecshi0 = vec_packs(veci0, veci1);
                    __vshi vecshi1 = vec_packs(veci2, veci3);
                    __vshi vecshi2 = vec_packs(veci4, veci5);
                    __vshi vecshi3 = vec_packs(veci6, veci7);

                    __vchari vec0 = vec_packs(vecshi0, vecshi1);
                    __vchari vec1 = vec_packs(vecshi2, vecshi3);

                    return Vec256<c10::qint8> {vec0, vec1};
                }

                Vec256<c10::qint8> __inline_attrs relu(Vec256<c10::qint8> zero_point) {
                    return Vec256<c10::qint8> {vec_max(_vec0, zero_point._vec0),
                        vec_max(_vec1, zero_point._vec1)
                    };
                }

                Vec256<c10::qint8> __inline_attrs
                    relu6(Vec256<c10::qint8> zero_point, Vec256<c10::qint8> q_six) {

                    return Vec256<c10::qint8> {
                        vec_min(vec_max(_vec0, zero_point._vec0), q_six._vec0),
                            vec_min(vec_max(_vec1, zero_point._vec1), q_six._vec1)
                    };
                }

                void dump() const {
                    value_type vals[size()];
                    store((void*)vals);
                    for (size_t i = 0; i < size(); ++i) {
                        std::cout << (int)(vals[i]) << " ";
                    }
                    std::cout << std::endl;
                }
            };

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

                Vec256() {}

            public:
                static constexpr size_t size() {
                    return 32;
                }

                static constexpr size_t float_num_vecs() {
                    return 4;
                }

                using float_vec_return_type = std::array<Vec256<float>, 4>;
                using value_type = typename c10::quint8::underlying;
                using vec_internal_type = __vchar;
                // Broadcast constructor
                __inline_attrs Vec256(const c10::quint8& val)
                    : _vec0(vec_splats(val.val_)), _vec1(vec_splats(val.val_)) {}

                __inline_attrs Vec256(const Vec256<c10::quint8>& other)
                    : _vec0{ other._vec0 }, _vec1(other._vec1) {}

                __inline_attrs Vec256(__vchar v1, __vchar v2) : _vec0{ v1 }, _vec1{ v2 } {}
                __inline_attrs Vec256(__vcharb v1, __vcharb v2) : _vecb0{ v1 }, _vecb1{ v2 } {}

                inline __inline_attrs const vec_internal_type& vec0() const {
                    return _vec0;
                }
                inline __inline_attrs const vec_internal_type& vec1() const {
                    return _vec1;
                }


                static __inline_attrs Vec256<c10::quint8> loadu(const void* ptr, size_t count = size()) {
                    if (count == size()) {
                        return Vec256<c10::quint8> {
                            vec_vsx_ld(offset0, reinterpret_cast<const value_type*>(ptr)),
                                vec_vsx_ld(offset16, reinterpret_cast<const value_type*>(ptr))
                        };
                    }
                    __at_align32__ value_type tmp_values[size()];
                    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

                    return Vec256<c10::quint8> {
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


            public:
                float_vec_return_type __inline_attrs dequantize(
                    Vec256<float> scale,
                    Vec256<float> zero_point,
                    Vec256<float> scale_zp_premul) const {
                    //unpacking unsigned as signed
                    __vshi vecshi0 = vec_unpackh((__vchari)_vec0);
                    __vshi vecshi1 = vec_unpackl((__vchari)_vec0);

                    __vshi vecshi2 = vec_unpackh((__vchari)_vec1);
                    __vshi vecshi3 = vec_unpackl((__vchari)_vec1);

                    __vui veci0 = (__vui)vec_unpackh(vecshi0);
                    __vui veci1 = (__vui)vec_unpackl(vecshi0);

                    __vui veci2 = (__vui)vec_unpackh(vecshi1);
                    __vui veci3 = (__vui)vec_unpackl(vecshi1);

                    __vui veci4 = (__vui)vec_unpackh(vecshi2);
                    __vui veci5 = (__vui)vec_unpackl(vecshi2);

                    __vui veci6 = (__vui)vec_unpackh(vecshi3);
                    __vui veci7 = (__vui)vec_unpackl(vecshi3);

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
                            vec_madd(scale_vec1, vecf1_0, scale_zp_premul1)
                        },
                        Vec256<float>{
                            vec_madd(scale_vec0, vecf0_1, scale_zp_premul0),
                            vec_madd(scale_vec1, vecf1_1, scale_zp_premul1)
                        },
                        Vec256<float>{
                            vec_madd(scale_vec0, vecf0_2, scale_zp_premul0),
                            vec_madd(scale_vec1, vecf1_2, scale_zp_premul1)
                        },
                        Vec256<float>{
                            vec_madd(scale_vec0, vecf0_3, scale_zp_premul0),
                            vec_madd(scale_vec1, vecf1_3, scale_zp_premul1)
                        }
                    };
                }

                static Vec256<c10::quint8> quantize(
                    const float_vec_return_type& rhs,
                    float scale,
                    int32_t zero_point,
                    float inverse_scale) {
                    // constexpr int32_t min_val = std::numeric_limits<value_type>::min();
                    // constexpr int32_t max_val = std::numeric_limits<value_type>::max();

                    __vf inverse_scale_v = vec_splats(inverse_scale);
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

                    vecf0 = vec_madd(vecf0, inverse_scale_v, vec_zero_point);
                    vecf1 = vec_madd(vecf1, inverse_scale_v, vec_zero_point);
                    vecf2 = vec_madd(vecf2, inverse_scale_v, vec_zero_point);
                    vecf3 = vec_madd(vecf3, inverse_scale_v, vec_zero_point);

                    vecf4 = vec_madd(vecf4, inverse_scale_v, vec_zero_point);
                    vecf5 = vec_madd(vecf5, inverse_scale_v, vec_zero_point);
                    vecf6 = vec_madd(vecf6, inverse_scale_v, vec_zero_point);
                    vecf7 = vec_madd(vecf7, inverse_scale_v, vec_zero_point);

                    __vui veci0 = vec_unsigned(vecf0);
                    __vui veci1 = vec_unsigned(vecf1);
                    __vui veci2 = vec_unsigned(vecf2);
                    __vui veci3 = vec_unsigned(vecf3);

                    __vui veci4 = vec_unsigned(vecf4);
                    __vui veci5 = vec_unsigned(vecf5);
                    __vui veci6 = vec_unsigned(vecf6);
                    __vui veci7 = vec_unsigned(vecf7);

                    // vec_packs CLAMP already
                    __vsh vecshi0 = vec_packsu(veci0, veci1);
                    __vsh vecshi1 = vec_packsu(veci2, veci3);
                    __vsh vecshi2 = vec_packsu(veci4, veci5);
                    __vsh vecshi3 = vec_packsu(veci6, veci7);

                    __vchar vec0 = vec_packsu(vecshi0, vecshi1);
                    __vchar vec1 = vec_packsu(vecshi2, vecshi3);

                    return Vec256<c10::quint8> {vec0, vec1};
                }

                Vec256<c10::quint8> __inline_attrs relu(Vec256<c10::quint8> zero_point) {
                    return Vec256<c10::quint8> {vec_max(_vec0, zero_point._vec0),
                        vec_max(_vec1, zero_point._vec1)
                    };
                }

                Vec256<c10::quint8> __inline_attrs
                    relu6(Vec256<c10::quint8> zero_point, Vec256<c10::quint8> q_six) {
                    return Vec256<c10::quint8> {
                        vec_min(vec_max(_vec0, zero_point._vec0), q_six._vec0),
                            vec_min(vec_max(_vec1, zero_point._vec1), q_six._vec1)
                    };
                }

                void dump() const {
                    value_type vals[size()];
                    store((void*)vals);
                    for (size_t i = 0; i < size(); ++i) {
                        std::cout << (int)(vals[i]) << " ";
                    }
                    std::cout << std::endl;
                }
            };

            template <>
            struct Vec256<c10::qint32> {
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

                Vec256() {}

            public:
                static constexpr size_t size() {
                    return 8;
                }

                static constexpr size_t float_num_vecs() {
                    return 1;
                }

                using float_vec_return_type = std::array<Vec256<float>, 1>;
                using value_type = c10::qint32::underlying;
                using vec_internal_type = __vi;

                __inline_attrs Vec256(__vi v1, __vi v2) : _vec0{ v1 }, _vec1{ v2 } {}
                __inline_attrs Vec256(__vib v1, __vib v2) : _vecb0{ v1 }, _vecb1{ v2 } {}

                Vec256(const c10::qint32& val)
                    : _vec0(vec_splats(val.val_)), _vec1(vec_splats(val.val_)) {}

                static Vec256<c10::qint32> __inline_attrs
                    loadu(const void* ptr, size_t count = size()) {
                    if (count == size()) {
                        return Vec256<c10::qint32> {
                            vec_vsx_ld(offset0, reinterpret_cast<const value_type*>(ptr)),
                                vec_vsx_ld(offset16, reinterpret_cast<const value_type*>(ptr))
                        };
                    }

                    __at_align32__ value_type tmp_values[size()];
                    std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

                    return Vec256<c10::qint32> {
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

                inline __inline_attrs const vec_internal_type& vec0() const {
                    return _vec0;
                }
                inline __inline_attrs const vec_internal_type& vec1() const {
                    return _vec1;
                }

                float_vec_return_type dequantize(
                    Vec256<float> scale,
                    Vec256<float> zero_point,
                    Vec256<float> scale_zp_premul) const {
                    __vf float_vals0 = vec_float(_vec0);
                    __vf float_vals1 = vec_float(_vec1);
                    __vf scale_vec0 = scale.vec0();
                    __vf scale_vec1 = scale.vec1();
                    __vf scale_zp_premul0 = scale_zp_premul.vec0();
                    __vf scale_zp_premul1 = scale_zp_premul.vec1();
                    return { Vec256<float>{
                            vec_madd(scale_vec0, float_vals0, scale_zp_premul0),
                            vec_madd(scale_vec1, float_vals1, scale_zp_premul1)
                        }
                    };
                }

                static Vec256<c10::qint32> quantize(
                    const float_vec_return_type& rhs,
                    float scale,
                    int32_t zero_point,
                    float inverse_scale) {
                    Vec256<c10::qint32> retval;

                    const __vi vmin = vec_splats(std::numeric_limits<value_type>::min());
                    const __vi vmax = vec_splats(std::numeric_limits<value_type>::max());
                    __vf inverse_scale_v = vec_splats(inverse_scale);
                    __vf vec_zero_point = vec_splats((float)zero_point);
                    Vec256<float> vf0 = rhs[0];
                    __vf vecf0 = vf0.vec0();
                    __vf vecf1 = vf0.vec1();

                    vecf0 = vec_madd(vecf0, inverse_scale_v, vec_zero_point);
                    vecf1 = vec_madd(vecf1, inverse_scale_v, vec_zero_point);

                    __vi veci0 = vec_min(vmax, vec_max(vmin, vec_signed(vecf0)));
                    __vi veci1 = vec_min(vmax, vec_max(vmin, vec_signed(vecf1)));

                    return Vec256<c10::qint32> {veci0, veci1};
                }

                Vec256<c10::qint32> relu(Vec256<c10::qint32> zero_point) {
                    return Vec256<c10::qint32> {vec_max(_vec0, zero_point._vec0),
                        vec_max(_vec1, zero_point._vec1)
                    };
                }

                Vec256<c10::qint32> relu6(
                    Vec256<c10::qint32> zero_point,
                    Vec256<c10::qint32> q_six) {
                    return Vec256<c10::qint32> {
                        vec_min(vec_max(_vec0, zero_point._vec0), q_six._vec0),
                            vec_min(vec_max(_vec1, zero_point._vec1), q_six._vec1)
                    };
                }

                void dump() const {
                    value_type vals[size()];
                    store((void*)vals);
                    for (size_t i = 0; i < size(); ++i) {
                        std::cout << (int)(vals[i]) << " ";
                    }
                    std::cout << std::endl;
                }
            };

#endif

        } // namespace
    } // namespace vec256
} // namespace at