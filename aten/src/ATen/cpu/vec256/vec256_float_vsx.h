#pragma once

#include <ATen/cpu/vec256/intrinsics.h>
#include <ATen/cpu/vec256/vec256_base.h>
#include <ATen/cpu/vec256/vsx_macros.h> 
namespace at {
	namespace vec256 {
		// See Note [Acceptable use of anonymous namespace in header]

		namespace {

#if defined(__VSX__)

//#Constants
			static const __vf zero = vec_splats(0.f);
			static const __vf half = vec_splats(0.5f);
			static const __vf one = vec_splats(1.f);
			static const __vf two = vec_splats(2.0f);
			static const __vf _4div_pi = vec_splats(1.27323954473516f);
			static const __vf v_inf = (__vf)vec_splats(0x7f800000u);
			static const __vf v_nan = (__vf)vec_splats(0x7fffffff);
			static const __vf log10e_inv = vec_splats(0.43429448190325176f);
			static const __vf log2e_inv = vec_splats(1.4426950408889634f);

			static const __vf cephes_SQRTHF = vec_splats(0.707106781186547524f);
			static const __vf coscof_p0 = vec_splats(2.443315711809948E-005f);
			static const __vf coscof_p1 = vec_splats(-1.388731625493765E-003f);
			static const __vf coscof_p2 = vec_splats(4.166664568298827E-002f);
			static const __vf exp_hi = vec_splats(104.f);
			static const __vf exp_lo = vec_splats(-104.f);
			static const __vf exp_p0 = vec_splats(0.000198527617612853646278381f);
			static const __vf exp_p1 = vec_splats((0.00139304355252534151077271f));
			static const __vf exp_p2 = vec_splats(0.00833336077630519866943359f);
			static const __vf exp_p3 = vec_splats(0.0416664853692054748535156f);
			static const __vf exp_p4 = vec_splats(0.166666671633720397949219f);
			static const __vf exp_p5 = vec_splats(0.5f);
			static const __vf log_p0 = vec_splats(7.0376836292E-2f);
			static const __vf log_p1 = vec_splats(-1.1514610310E-1f);
			static const __vf log_p2 = vec_splats(1.1676998740E-1f);
			static const __vf log_p3 = vec_splats(-1.2420140846E-1f);
			static const __vf log_p4 = vec_splats(+1.4249322787E-1f);
			static const __vf log_p5 = vec_splats(-1.6668057665E-1f);
			static const __vf log_p6 = vec_splats(+2.0000714765E-1f);
			static const __vf log_p7 = vec_splats(-2.4999993993E-1f);
			static const __vf log_p8 = vec_splats(+3.3333331174E-1f);
			static const __vf log_q1 = vec_splats(-2.12194440e-4f);
			static const __vf log_q2 = vec_splats(0.693359375f);
			static const __vf max_logf = vec_splats(88.02969187150841f);
			static const __vf max_numf = vec_splats(1.7014117331926442990585209174225846272e38f);
			static const __vf min_inf = (__vf)vec_splats(0xff800000u);
			static const __vf min_norm_pos = (__vf)vec_splats(0x0800000u);
			static const __vf minus_cephes_dp1 = vec_splats(-0.78515625f);
			static const __vf minus_cephes_dp2 = vec_splats(-2.4187564849853515625e-4f);
			static const __vf minus_cephes_dp3 = vec_splats(-3.77489497744594108e-8f);
			static const __vf negln2f_hi = vec_splats(-0.693145751953125f);
			static const __vf negln2f_lo = vec_splats(-1.428606765330187045e-06f);
			static const __vf p0 = vec_splats(2.03721912945E-4f);
			static const __vf p1 = vec_splats(8.33028376239E-3f);
			static const __vf p2 = vec_splats(1.66667160211E-1f);
			static const __vf sincof_p0 = vec_splats(-1.9515295891E-4f);
			static const __vf sincof_p1 = vec_splats(8.3321608736E-3f);
			static const __vf sincof_p2 = vec_splats(-1.6666654611E-1f);
			static const __vf tanh_0p625 = vec_splats(0.625f);
			static const __vf tanh_half_max = vec_splats(44.014845935754205f);
			static const __vf tanh_p0 = vec_splats(-5.70498872745E-3f);
			static const __vf tanh_p1 = vec_splats(2.06390887954E-2f);
			static const __vf tanh_p2 = vec_splats(-5.37397155531E-2f);
			static const __vf tanh_p3 = vec_splats(1.33314422036E-1f);
			static const __vf tanh_p4 = vec_splats(-3.33332819422E-1f);
			static const __vf vcheck = vec_splats((float)(1LL << 24));
			static const __vi v0x7f = vec_splats(0x7f);
			static const __vi vi_0 = vec_splats((int)(0));
			static const __vi vi_1 = vec_splats((int)1);
			static const __vi vi_2 = vec_splats((int)2);
			static const __vi vi_4 = vec_splats((int)4);
			static const __vi vi_inv1 = vec_splats((int)~1);
			static const __vib inv_mant_mask = (__vib)vec_splats((unsigned int)~0xff800000);
			static const __vib sign_mask = (__vib)vec_splats((int)0x80000000);
			static const __vui vu_29 = vec_splats(29u);
			static const __vui vu_23 = vec_splats(23u);


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
				using vec_internal_mask_type = __vib;

				static constexpr size_t size() {
					return 8;
				}
				Vec256() {}
				__inline_attrs Vec256(__vf v1, __vf v2) : _vec0{ v1 }, _vec1{ v2 } {}
				__inline_attrs Vec256(__vib v1, __vib v2) : _vecb0{ v1 }, _vecb1{ v2 } {}
				__inline_attrs Vec256(float scalar)
					: _vec0{ vec_splats(scalar) }, _vec1{ vec_splats(scalar) } {}
				__inline_attrs Vec256(
					float scalar1,
					float scalar2,
					float scalar3,
					float scalar4,
					float scalar5,
					float scalar6,
					float scalar7,
					float scalar8)
					: _vec0{ __vf{scalar1, scalar2, scalar3, scalar4} },
					_vec1{ __vf{scalar5, scalar6, scalar7, scalar8} } {}
				inline __inline_attrs const vec_internal_type& vec0() const {
					return _vec0;
				}
				inline __inline_attrs const vec_internal_type& vec1() const {
					return _vec1;
				}

				template <uint64_t mask>
				static std::endable_if_t<mask == 0, Vec256<float>> __inline_attrs
					blend(const Vec256<float>& a, const Vec256<float>& b) {
					return a;
				}

				template <uint64_t mask>
				static std::endable_if_t<(mask & 255) == 255, Vec256<float>>
					__inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
					return b;
				}

				template <uint64_t mask>
				static std::endable_if_t<mask == 15, Vec256<float>> __inline_attrs
					blend(const Vec256<float>& a, const Vec256<float>& b) {
					return Vec256<float>{b._vec0, a._vec1};
				}

				template <uint64_t mask>
				static std::endable_if_t<(mask > 0 && mask < 15), Vec256<float>>
					__inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
					constexpr uint32_t g0 = (mask & 1) * 0xffffffff;
					constexpr uint32_t g1 = ((mask & 2) >> 1) * 0xffffffff;
					constexpr uint32_t g2 = ((mask & 4) >> 2) * 0xffffffff;
					constexpr uint32_t g3 = ((mask & 8) >> 3) * 0xffffffff;
					const __vib mask_1st = (__vib){ g0, g1, g2, g3 };
					return Vec256<float>{(__vf)vec_sel(a._vec0, b._vec0, (__vib)mask_1st),
						a._vec1};
				}

				template <uint64_t mask>
				static std::endable_if_t<
					(mask > 15 && (mask & 255) != 255 && ((mask & 15) == 15)),
					Vec256<float>>
					__inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
					constexpr uint32_t mask2 = (mask & 255) >> 4;
					constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
					constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
					constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
					constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;
					const __vib mask_2nd = (__vib){ g0_2, g1_2, g2_2, g3_2 };
					// generated masks
					return Vec256<float>{b._vec0,
						(__vf)vec_sel(a._vec1, b._vec1, (__vib)mask_2nd)};
				}

				template <uint64_t mask>
				static std::endable_if_t<
					(mask > 15 && ((mask & 255) != 255) && ((mask & 15) == 0)),
					Vec256<float>>
					__inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {
					constexpr uint32_t mask2 = (mask & 255) >> 4;
					constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
					constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
					constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
					constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;
					const __vib mask_2nd = (__vib){ g0_2, g1_2, g2_2, g3_2 };
					// generated masks
					return Vec256<float>{a, (__vf)vec_sel(a._vec1, b._vec1, (__vib)mask_2nd)};
				}

				template <uint64_t mask>
				static std::endable_if_t<
					(mask > 15 && ((mask & 255) != 255) && ((mask & 15) != 0) &&
						((mask & 15) != 15)),
					Vec256<float>>
					__inline_attrs blend(const Vec256<float>& a, const Vec256<float>& b) {

					constexpr uint32_t g0 = (mask & 1) * 0xffffffff;
					constexpr uint32_t g1 = ((mask & 2) >> 1) * 0xffffffff;
					constexpr uint32_t g2 = ((mask & 4) >> 2) * 0xffffffff;
					constexpr uint32_t g3 = ((mask & 8) >> 3) * 0xffffffff;
					constexpr uint32_t mask2 = (mask & 255) >> 4;
					constexpr uint32_t g0_2 = (mask2 & 1) * 0xffffffff;
					constexpr uint32_t g1_2 = ((mask2 & 2) >> 1) * 0xffffffff;
					constexpr uint32_t g2_2 = ((mask2 & 4) >> 2) * 0xffffffff;
					constexpr uint32_t g3_2 = ((mask2 & 8) >> 3) * 0xffffffff;

					const __vib mask_1st = (__vib){ g0, g1, g2, g3 };
					const __vib mask_2nd = (__vib){ g0_2, g1_2, g2_2, g3_2 };
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
						return Vec256<float>{
							vec_vsx_ld(offset0, reinterpret_cast<const value_type*>(ptr)),
								vec_vsx_ld(offset16, reinterpret_cast<const value_type*>(ptr))};
					}

					__at_align32__  value_type tmp_values[size()];
					std::memcpy(tmp_values, ptr, std::min(count, size()) * sizeof(value_type));

					return Vec256<float>{
						vec_vsx_ld(offset0, tmp_values),
							vec_vsx_ld(offset16, tmp_values)};
				}
				void __inline_attrs store(void* ptr, size_t count = size()) const {
					if (count == size()) {
						vec_vsx_st(_vec0, offset0, reinterpret_cast<value_type*>(ptr));
						vec_vsx_st(_vec1, offset16, reinterpret_cast<value_type*>(ptr));
					}
					else if (count > 0) {
						__at_align32__  value_type tmp_values[size()];
						vec_vsx_st(_vec0, offset0, tmp_values);
						vec_vsx_st(_vec1, offset16, tmp_values);
						std::memcpy(ptr, tmp_values, std::min(count, size()) * sizeof(value_type));
					}
				}

				const float& operator[](int idx) const = delete;
				float& operator[](int idx) = delete;

				Vec256<float> map(float (*f)(float)) const {
					Vec256<float> ret;
					for (int i = 0; i < 4; i++) {
						ret._vec0[i] = f(_vec0[i]);
					}
					for (int i = 0; i < 4; i++) {
						ret._vec1[i] = f(_vec1[i]);
					}
					return ret;
				}

				Vec256<float> mapbi(float (*f)(float, float), const Vec256<float>& other) const {
					Vec256<float> ret;
					for (int i = 0; i < 4; i++) {
						ret._vec0[i] = f(_vec0[i],other._vec0[i]);
					}
					for (int i = 0; i < 4; i++) {
						ret._vec1[i] = f(_vec1[i],other._vec1[i]);
					}
					return ret;
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
					// implementation logic from avx_mathfun with some modifications from sleef
			 
					/* Express e**x = e**g 2**n
					  *   = e**g e**( n loge(2) )
					  *   = e**( g + n loge(2) )
					  */
					auto tmp_x = *this;
					auto fx = (tmp_x * log2e_inv).round();

					auto x = fmadd(fx, negln2f_hi, tmp_x);
					x = fmadd(fx, negln2f_lo, x);
					auto z = x * x;
					auto y = fmadd(exp_p0, x, exp_p1);
					y = fmadd(y, x, exp_p2);
					y = fmadd(y, x, exp_p3);
					y = fmadd(y, x, exp_p4);
					y = fmadd(y, x, exp_p5);
					y = fmadd(y, z, x) + one;

					/* vm_pow2n 2^n */
					__vi imm0 = vec_signed(fx._vec0);
					__vi imm1 = vec_signed(fx._vec1);
					//this pow2n logic is  from Sleef code
					__vi imm00 = imm0 >> 1; //>>1 
					__vi imm01 = imm1 >> 1;
					__vi imm10 = imm0 - imm00;
					__vi imm11 = imm1 - imm01;
					imm00 = (imm00 + v0x7f) << vu_23;
					imm01 = (imm01 + v0x7f) << vu_23;
					imm10 = (imm10 + v0x7f) << vu_23;
					imm11 = (imm11 + v0x7f) << vu_23;
					// treat imm as float vector without conversion

					y._vec0 = (y._vec0 * (__vf)imm00) * (__vf)imm10;
					y._vec1 = (y._vec1 * (__vf)imm01) * (__vf)imm11;
					//boundary check
					auto tmp = select(y, v_inf, exp_hi <= tmp_x);
					y = select(tmp, zero, tmp_x < exp_lo);

					return y;
				}
				Vec256<float> expm1() const {
					return exp() - one; 
				}

				Vec256<float> __inline_attrs log() const {
					auto temp = *this;
					auto invalid_mask = temp < zero;
					auto x = maximum(temp, min_norm_pos); /* cut off denormalized stuff */

					__vi imm0 = vec_sr(__vi(x._vec0), vu_23);
					__vi imm1 = vec_sr(__vi(x._vec1), vu_23);
					/* keep only the fractional part */
					x = x & inv_mant_mask;
					x = x | half;

					imm0 = imm0 - v0x7f;
					imm1 = imm1 - v0x7f;
					Vec256<float> ex;
					ex._vec0 = vec_float(imm0);
					ex._vec1 = vec_float(imm1);
					ex = ex + one;
					auto mask = x < cephes_SQRTHF;
					auto t = x & mask;
					x = x - one;
					ex = ex - (one & mask);
					x = x + t;
					auto z = x * x;
					auto y = fmadd(x, log_p0, log_p1);
					y = fmadd(y, x, log_p2);
					y = fmadd(y, x, log_p3);
					y = fmadd(y, x, log_p4);
					y = fmadd(y, x, log_p5);
					y = fmadd(y, x, log_p6);
					y = fmadd(y, x, log_p7);
					y = fmadd(y, x, log_p8);
					y = y * x * z;
					y = fmadd(ex, log_q1, y);
					y = y - z * half;
					x = x + y;
					x = fmadd(ex, log_q2, x);
					x = select(x, v_nan, invalid_mask); // negative arg will be NAN
					 //zero is -inf  
					x = select(x, min_inf, temp == zero);
					return x;
				}
				Vec256<float> __inline_attrs log10() const {
					return log() * log10e_inv;
				}
				Vec256<float> __inline_attrs log1p() const {
					return ((*this)+one).log();
				}
				Vec256<float> __inline_attrs log2() const {
					return log() * log2e_inv;
				}
				Vec256<float> __inline_attrs ceil() const {
					return Vec256<float>{vec_ceil(_vec0), vec_ceil(_vec1)};
				}
				Vec256<float> __inline_attrs cos() const {
					/* take the absolute value */
					auto x = abs();
					/* extract the sign bit (upper one) */
					auto sign_bit = (*this) & sign_mask;
					/* scale by 4/Pi */
					auto y = x * _4div_pi;
					/* store the integer part of y in mm0 */
					/* j=(j+1) & (~1) (see the cephes sources) */
					__vi imm0 = (vec_signed(y._vec0) + vi_1) & vi_inv1;
					__vi imm1 = (vec_signed(y._vec1) + vi_1) & vi_inv1;
					y._vec0 = vec_float(imm0);
					y._vec1 = vec_float(imm1);

					imm0 = imm0 - vi_2;
					imm1 = imm1 - vi_2;
					Vec256<float>  poly_mask;
					/* get the swap sign flag */
					__vi tmp0 = vec_and(vec_nand(imm0, imm0), vi_4);
					__vi tmp1 = vec_and(vec_nand(imm1, imm1), vi_4);
					sign_bit._vecb0 = (__vib)vec_sl(tmp0, vu_29);
					sign_bit._vecb1 = (__vib)vec_sl(tmp1, vu_29);
					/* get the polynom selection mask
					   there is one polynom for 0 <= x <= Pi/4
					   and another one for Pi/4<x<=Pi/2
					   Both branches will be computed.
					*/
					poly_mask._vecb0 = (__vib)vec_cmpeq((imm0 & vi_2), vi_0);
					poly_mask._vecb1 = (__vib)vec_cmpeq((imm1 & vi_2), vi_0);

					/* The magic pass: "Extended precision modular arithmetic"
					   x = ((x - y * DP1) - y * DP2) - y * DP3; */
					x = fmadd(y, minus_cephes_dp1, x);
					x = fmadd(y, minus_cephes_dp2, x);
					x = fmadd(y, minus_cephes_dp3, x);

					/* Evaluate the first polynom  (0 <= x <= Pi/4) */
					auto z = x * x;
					y = fmadd(coscof_p0, z, coscof_p1);
					y = fmadd(y, z, coscof_p2);
					y = y * z * z;
					y = y - z * half + one;

					/* Evaluate the second polynom  (Pi/4 <= x <= 0) */
					auto y_2 = fmadd(sincof_p0, z, sincof_p1);
					y_2 = fmadd(y_2, z, sincof_p2);
					y_2 = y_2 * z;
					y_2 = fmadd(y_2, x, x);

					/* select the correct result from the two polynoms */
					y = select(y, y_2, poly_mask);
					/* update the sign */
					y = y ^ sign_bit;

					return y;
				}
				Vec256<float> __inline_attrs cosh() const {
					// cosh = 1/2 * (e^x + e^-x) 
					__vf recp_vec0, recp_vec1, v0, v1;
					auto e_x = abs().exp();
					return (e_x + one / e_x) * half; 
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
					/* take the absolute value and xtract sign*/
					auto x = abs();
					auto sign_bit = (*this) & sign_mask;

					/* scale by 4/Pi */
					auto y = x * _4div_pi;
					/* store the integer part of y in mm0 */

					/* j=(j+1) & (~1) (see the cephes sources) */ 
					__vi imm0 = (vec_signed(y._vec0) + vi_1) & vi_inv1;
					__vi imm1 = (vec_signed(y._vec1) + vi_1) & vi_inv1; 
					y._vec0 = vec_float(imm0);
					y._vec1 = vec_float(imm1);
					/* get the swap sign flag */
					Vec256<float> swap_sign_bit, poly_mask;
					swap_sign_bit._vecb0 = (__vib)vec_sl(imm0 & vi_4, vu_29);
					swap_sign_bit._vecb1 = (__vib)vec_sl(imm1 & vi_4, vu_29);
					/* get the polynom selection mask
					   there is one polynom for 0 <= x <= Pi/4
					   and another one for Pi/4<x<=Pi/2

					   Both branches will be computed.
					*/
					poly_mask._vecb0 = vec_cmpeq((imm0 & vi_2), vi_0); 
					poly_mask._vecb1 = vec_cmpeq((imm1 & vi_2), vi_0);
					sign_bit = sign_bit ^ swap_sign_bit; //xor operation 

					/* The magic pass: "Extended precision modular arithmetic"
					   x = ((x - y * DP1) - y * DP2) - y * DP3; */
					x = fmadd(y, minus_cephes_dp1, x); 
					x = fmadd(y, minus_cephes_dp2, x);
					x = fmadd(y, minus_cephes_dp3, x);

					/* Evaluate the first polynom  (0 <= x <= Pi/4) */
					auto z = x * x; 
					y = fmadd(coscof_p0, z, coscof_p1); 
					y = fmadd(y, z, coscof_p2); 
					y = y * z * z;
					y = y - z * half + one;

					/* Evaluate the second polynom  (Pi/4 <= x <= 0) */
					auto y2 = fmadd(sincof_p0, z, sincof_p1);
					y2 = fmadd(y2, z, sincof_p2);
					y2 =  y2 * z;
					y2 = fmadd(y2, x, x);;

					/* select the correct result from the two polynoms */
					y = select(y, y2, poly_mask);
					y = y ^ sign_bit;

					return y;
				}
				Vec256<float> __inline_attrs sinh() const {

					auto temp_abs = abs();
					//get exponent
					auto ret = temp_abs.exp();
					auto recp = half / ret;
					auto v = half * ret - recp;
					/* extract the sign bit (upper one) */
					auto sign_bit = (*this) & sign_mask;
					v = v | sign_bit;
					auto z = temp_abs * temp_abs;
					auto y = p0 * z + p1;
					y = y * z + p2;
					y = y * z * temp_abs + temp_abs;
					//check and select 
					return select(y, v, temp_abs > one);
				}
				Vec256<float> __inline_attrs tan() const {
					return map(std::tan);
				}
				Vec256<float> __inline_attrs tanh() const {

					auto x = *this;
					auto vabs = abs();
					//get exponent
					auto exp2x = (vabs + vabs).exp();
					auto vv = one - two / (exp2x + one);
					/* extract the sign bit (upper one) */
					auto sign_bit = (*this) & sign_mask;
					auto z = vabs * vabs;
					auto y   = fmadd(tanh_p0, z, tanh_p1);
					auto tmp = fmadd(y, z, tanh_p2);
					y   = fmadd(tmp, z, tanh_p3);
					tmp = fmadd(y, z, tanh_p4);
					y = tmp * z;
					tmp = fmadd(y, x, x);
					//add sign
					vv = vv | sign_bit;
					//check and select
					auto sel_mask = vabs >= tanh_0p625;
					auto max_mask = vabs > tanh_half_max;
					auto max_ret = sign_bit ^ one;
					return select(select(tmp, vv, sel_mask), max_ret, max_mask);
				}
				Vec256<float> __inline_attrs trunc() const {
					return Vec256<float>{vec_trunc(_vec0), vec_trunc(_vec1)};
				}

				Vec256<float> __inline_attrs frac() const {
					return Vec256<float>{
						vec_sub(_vec0, vec_trunc(_vec0)),
							vec_sub(_vec1, vec_trunc(_vec1))};
				}

				Vec256<float> __inline_attrs sqrt() const {
					return Vec256<float>{vec_sqrt(_vec0), vec_sqrt(_vec1)};
				}
				Vec256<float> __inline_attrs reciprocal() const {
					return one / (*this);
				}
				Vec256<float> __inline_attrs rsqrt() const {
					return one / sqrt();
				}

				Vec256<float> __inline_attrs pow(const Vec256<float>& pow_exp) const {

					auto x = *this;
					auto sign_bit = (*this) & sign_mask;
					// |b|
					auto pow_exp_abs = pow_exp.abs();
					auto pow_exp_trunc = pow_exp.trunc();
					Vec256<float> odd_mask;
					odd_mask._vecb0 = (vec_signed(pow_exp._vec0) & vi_1) != vi_0;
					odd_mask._vecb1 = (vec_signed(pow_exp._vec1) & vi_1) != vi_0;
					/* using ln fuction */
					auto temp = (abs().log() * pow_exp).exp();

					// is odd or even check from Sleef

					auto is_int = (pow_exp == pow_exp_trunc) | (pow_exp_abs >= vcheck);
					auto is_odd = odd_mask &  is_int & (pow_exp_abs < vcheck);
					// if even then then pow result should be absolute 
					auto temp_sign = temp | sign_bit; // copy_sign
					auto out = select(temp, temp_sign, is_odd);
					/* x<0 and y != N, then NAN */
					auto tmp_m = (pow_exp.floor() != pow_exp) & (x < zero);
					auto out1 = select(out, v_nan, tmp_m);
					/* y = 0 then 1 */
					tmp_m = pow_exp_abs == zero;
					return select(out1, one, tmp_m);
				}

 

				    DEFINE_CMP_OP(operator ==, float, vec_cmpeq)
					DEFINE_CMP_OP(operator !=, float, vec_cmpne)
					DEFINE_CMP_OP(operator <, float, vec_cmplt)
					DEFINE_CMP_OP(operator <=, float, vec_cmple)
					DEFINE_CMP_OP(operator >, float, vec_cmpgt)
					DEFINE_CMP_OP(operator >=, float, vec_cmpge)
					DEFINE_FRIEND_BINARY_OP(operator +, float, vec_add)
					DEFINE_FRIEND_BINARY_OP(operator -, float, vec_sub)
					DEFINE_FRIEND_BINARY_OP(operator *, float, vec_mul)
					DEFINE_FRIEND_BINARY_OP(operator /, float, vec_div)
					DEFINE_FRIEND_BINARY_OP(maximum, float, vec_max)
					DEFINE_FRIEND_BINARY_OP(minimum, float, vec_min)
					DEFINE_FRIEND_BITWISE_OP(operator &, float, vec_and)
					DEFINE_FRIEND_BITWISE_OP(operator |, float, vec_or)
					DEFINE_FRIEND_BITWISE_OP(operator ^, float, vec_xor)
					DEFINE_FRIEND_TERNARY_OP(fmadd, float, vec_madd)
					DEFINE_FRIEND_TERNARY_OP_MASK(select, float, vec_sel)



			};

#endif

		} // namespace
	} // namespace vec256
} // namespace at
