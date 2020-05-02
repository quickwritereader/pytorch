#pragma once

#include <gtest/gtest.h>
#include <limits>
#include <iostream>
#include <exception>
#include <ATen/native/Math.h>
#include <ATen/cpu/vec256/vec256.h>
#include <vector>
#include <random>
#include <functional>


#define  RESOLVE_OVERLOAD(...) \
        [](auto&&...args)->decltype(auto) {return __VA_ARGS__(std::forward<decltype(args)>(args)...); }


template<typename T>
using Vec = typename at::vec256::Vec256<T>;
using vfloat = Vec<float>;
using vdouble = Vec<double>;
using vlong = Vec<int64_t>;
using vint = Vec<int32_t>;
using vshort = Vec<int16_t>;
using vqint8 = Vec<c10::qint8>;
using vquint8 = Vec<c10::quint8>;
using vqint = Vec<c10::qint32>;

template<typename T>
using ValueType = typename T::value_type;

template<class T, size_t N>
constexpr size_t size(T(&)[N]) {
	return N;
}

constexpr size_t UNARY = 1;
constexpr size_t BINARY = 2;
constexpr size_t TERNARY = 3;



template<size_t argc, typename T, typename Op>
typename std::enable_if_t<(argc == UNARY), Vec<T>>
call_op(Op op, T* ptr, size_t vsize) {
	return op(ptr);
}

template<size_t argc, typename T, typename Op>
typename std::enable_if_t<(argc == BINARY), Vec<T>>
call_op(Op op, T* ptr, size_t vsize) {
	return op(ptr, ptr + vsize);
}

template<size_t argc, typename T, typename Op >
typename std::enable_if_t<(argc == TERNARY), Vec<T>>
call_op(Op op, T* ptr, size_t vsize) {
	return op(ptr, ptr + vsize, ptr + vsize + vsize);
}

template< size_t argc, typename Filter, typename T>
typename std::enable_if_t<std::is_same<Filter, nullptr_t>::value, void>
call_filter(Filter  filter, T* ptr, size_t vsize) {
}

template< size_t argc, typename Filter, typename T>
typename std::enable_if_t<!std::is_same<Filter, nullptr_t>::value && (argc == UNARY), void>
call_filter(Filter  filter, T* ptr, size_t vsize) {
	return filter(ptr, vsize);
}

template< size_t argc, typename Filter, typename T>
typename std::enable_if_t<!std::is_same<Filter, nullptr_t>::value && (argc == BINARY), void>
call_filter(Filter filter, T* ptr, size_t vsize) {
	return filter(ptr, ptr + vsize, vsize);
}

template< size_t argc, typename Filter, typename T>
typename std::enable_if_t<!std::is_same<Filter, nullptr_t>::value && (argc == TERNARY), void>
call_filter(Filter  filter, T* ptr, size_t vsize) {
	return filter(ptr, ptr + vsize, ptr + vsize + vsize, vsize);
}

template<int N>
struct BitStr {
	using type = uintmax_t;
};

template<>
struct BitStr<8> {
	using type = uint64_t;
};

template<>
struct BitStr<4> {
	using type = uint32_t;
};

template<>
struct BitStr<2> {
	using type = uint16_t;
};

template<>
struct BitStr<1> {
	using type = uint8_t;
};

template<typename T>
using BitValueType = typename  BitStr<sizeof(ValueType<T>)>::type;

template <typename T>
struct DomainRange {
	T start; //start [
	T end;  // end is not included
	//one could use  nextafter for including his end case for tests
};


template <typename T>
struct SpecArg {
	std::vector<T> Args;
	T desiredReturn;
};

template <typename T>
struct CheckWithinDomains {
	//each argument takes domain Range
	std::vector<DomainRange<T>> ArgsDomain;
	//check with error tolerance
	bool CheckWithAcceptance;
	T AcceptedError;
};

template <typename T>
class TestCaseBuilder;

template <typename T>
class TestingCase {
public:
	friend class TestCaseBuilder<T>;
	static TestCaseBuilder<T> getBuilder() {
		return TestCaseBuilder<T>{};
	}

	bool checkDefaultSpecials() const {
		return defaultSpecials;
	}

	bool checkNansAndInfinities() const {
		return test_nan_inf;
	}

	size_t getTrialCount() const {
		return trials;
	}

	bool isBitwise() const {
		return bitwise;
	}
	const std::vector<CheckWithinDomains<T>>& getDomains() const {
		return domains;
	}

	const std::vector<SpecArg<T>>& getCustomSpecials() const {
		return customSpecialCheck;
	}

private:
	// if domains is empty we will test default
	std::vector<CheckWithinDomains<T>> domains;
	bool defaultSpecials = false;
	std::vector<SpecArg<T>> customSpecialCheck;
	bool test_nan_inf = false;
	bool bitwise = false; //test bitlevel 
	size_t trials = 1;
};


template <typename T>
class TestCaseBuilder {
private:
	TestingCase<T> t_case;
public:

	TestCaseBuilder<T>& set(bool bitwise, bool allow_specials, bool test_nan_inf) {
		t_case.bitwise = bitwise;
		t_case.test_nan_inf = test_nan_inf;
		t_case.defaultSpecials = allow_specials;
		return *this;
	}

	TestCaseBuilder<T>& setTrialCount(size_t trial_count) {
		t_case.trials = trial_count;
		return *this;
	}

	TestCaseBuilder<T>& addDomain(const CheckWithinDomains<T>& domainCheck) {
		t_case.domains.emplace_back(domainCheck);
		return *this;
	}

	TestCaseBuilder<T>& addSpecial(const SpecArg<T>& specialArg) {
		t_case.customSpecialCheck.emplace_back(specialArg);
		return *this;
	}

	TestCaseBuilder<T>& testNansAndInfinities() {
		t_case.test_nan_inf = true & std::is_floating_point<T>::value;
		return *this;
	}

	TestCaseBuilder<T>& checkDefaultSpecials() {
		t_case.defaultSpecials = true;
		return *this;
	}

	TestCaseBuilder<T>& compareBitwise() {
		t_case.bitwise = true;
		return *this;
	}

	operator TestingCase<T> && () {
		return std::move(t_case);
	}


};

template<class To, class From>
typename std::enable_if<
	(sizeof(To) == sizeof(From)) &&
	std::is_trivially_copyable<From>::value&&
	std::is_trivial<To>::value,
	// this implementation requires that To is trivially default constructible
	To>::type
	bit_cast(const From& src) noexcept
{
	To dst;
	std::memcpy(&dst, &src, sizeof(To));
	return dst;
}

template<class To, class T>
To bit_cast_ptr(T* p, size_t N = sizeof(To)) noexcept {
	unsigned char p1[sizeof(To)] = {};
	std::memcpy(p1, p, std::min(N, sizeof(To)));
	return bit_cast<To>(p1);
}
//turn off optimization for this to work

template<typename T>
std::enable_if_t<std::is_same<T, double>::value, bool>
check_both_nan(T x, T y) {
	return  __isnan(x) && __isnan(y);  //(std::fpclassify(x) == FP_NAN && std::fpclassify(y) == FP_NAN);
}
template<typename T>
std::enable_if_t<std::is_same<T, float>::value, bool>
check_both_nan(T x, T y) {
	return  __isnanf(x) && __isnanf(y);  //(std::fpclassify(x) == FP_NAN && std::fpclassify(y) == FP_NAN);
}
template<typename T>
std::enable_if_t<!std::is_floating_point<T>::value, bool>
check_both_nan(T x, T y) {
	return false;
}

template<typename T>
T local_abs(T x) {
	return x > 0 ? x : -x;
}

template<typename T>
T reciprocal(T x) {
	return 1 / x;
}

template<typename T>
T rsqrt(T x) {
	return 1 / std::sqrt(x);
}

template<typename T>
T sqrt(T x) {
	return std::sqrt(x);
}

template<class T>
T  maximum(const T& a, const T& b) {
	return (a > b) ? a : b;
}

template<class T>
T  minimum(const T& a, const T& b) {
	return (a < b) ? a : b;
}

template<class T>
T  clamp(const T& a, const T& min, const T& max) {
	return a < min ? min : (a > max ? max : a);
}

template<class T>
T  clamp_max(const T& a, const T& max) {
	return a > max ? max : a;
}

template<class T>
T  clamp_min(const T& a, const T& min) {
	return a < min ? min : a;
}

template<class VT, size_t N>
void copy_interleave(VT(&vals)[N], VT(&interleaved)[N])
{
	static_assert(N % 2 == 0, "should be even");
	auto ptr1 = vals;
	auto ptr2 = vals + N / 2;
	for (size_t i = 0; i < N; i += 2) {
		interleaved[i] = *ptr1++;
		interleaved[i + 1] = *ptr2++;
	}
}

template<typename T>
std::enable_if_t<std::is_floating_point<T>::value, bool>
is_zero(T val) {
	return val == 0 || (std::is_floating_point<T>::value && std::fpclassify(val) == FP_ZERO);
}

template<typename T>
std::enable_if_t<!std::is_floating_point<T>::value, bool>
is_zero(T val) {
	return val == 0;
}

template<typename T>
void filter_mult_overflow(T* first, T* second, size_t vsize) {
	if (std::is_integral<T>::value == false) return;
	for (size_t i = 0; i < vsize; i++) {
		if (!is_zero(second[i])) {
			T c = (std::numeric_limits<T>::max() - 1) / second[i];
			if (abs(first[i]) >= c) {
				//correct first;
				first[i] = c;
			}
		}//is_zero
	}
}

template<typename T>
void filter_div_ub(T* first, T* second, size_t vsize) {
	if (std::is_integral<T>::value == false) return;
	for (size_t i = 0; i < vsize; i++) {
		if (is_zero(second[i])) {
			second[i] = 1;
		}
		else if (first[i] == std::numeric_limits<T>::min() && second[i] == -1) {
			second[i] = 1;
		}
	}
}

template<typename T>
void filter_op(T* f, T* s, size_t vsize, bool minus) {
	T max = std::numeric_limits<T>::max();
	T min = std::numeric_limits<T>::min();
	for (size_t i = 0; i < vsize; i++) {
		T a = f[i];
		T b = s[i];
		if (minus) {
			if (b == min) b = min + 1;
			b = -b;
		}
		bool sgn1 = a > 0;
		bool sgn2 = b > 0;
		if (sgn1 == sgn2) {
			if (sgn1 && a > max - b) {
				f[i] = max - b;
			}
			else if (!sgn1 && a < min - b) {
				f[i] = min - b;
			}
		}
	}
}

template<typename T>
void filter_add_overflow(T* first, T* second, size_t vsize) {
	if (std::is_integral<T>::value == false) return;
	return filter_op(first, second, vsize, false);
}

template<typename T>
void filter_minus_overflow(T* first, T* second, size_t vsize) {
	if (std::is_integral<T>::value == false) return;
	return filter_op(first, second, vsize, true);
}

template<typename T>
void filter_clamp(T* first, T* second, T* third, size_t vsize) {
	for (size_t i = 0; i < vsize; i++) {
		if (third[i] < second[i]) {
			T tmp = second[i];
			second[i] = tmp;
			third[i] = tmp;
		}
	}
}

template<typename T>
void filter_zero(T* ptr, size_t vsize) {
	T v = 1;
	for (size_t i = 0; i < vsize; i++) {
		ptr[i] = is_zero(ptr[i]) ? v : ptr[i];
	}
}


template<typename T>
void filter_int_minimum(T* ptr, size_t vsize) {
	if (!std::is_integral<T>::value) return;
	for (size_t i = 0; i < vsize; i++) {
		if (ptr[i] == std::numeric_limits<T>::min()) {
			ptr[i] = 0;
		}
	}
}



template<typename T>
void AssertVec256(T expected, T actual, std::function<std::string(int index)> get_details = {}, bool bitwise = false, bool check_absError = false, ValueType<T> absError = {}) {
	CACHE_ALIGN ValueType<T> exp[T::size()];
	CACHE_ALIGN ValueType<T> act[T::size()];
	expected.store(exp);
	actual.store(act);
	if (bitwise) {
		using bit_rep = BitValueType<T>;
		for (size_t i = 0; i < T::size(); i++) {
			bit_rep b_exp = bit_cast<bit_rep>(exp[i]);
			bit_rep b_act = bit_cast<bit_rep>(act[i]);
			ASSERT_EQ(b_exp, b_act) << (get_details ? get_details(i) : "");
		}
	}
	else if (check_absError) {
		for (size_t i = 0; i < T::size(); i++) {
			ASSERT_NEAR(exp[i], act[i], absError);
		}
	}
	else if (std::is_same<ValueType<T>, float>::value) {
		for (size_t i = 0; i < T::size(); i++) {
			if (!(check_both_nan(exp[i], act[i]))) {
				ASSERT_FLOAT_EQ(exp[i], act[i]) << (get_details ? get_details(i) : "");
			}
		}
	}
	else if (std::is_same<ValueType<T>, double>::value) {
		for (size_t i = 0; i < T::size(); i++) {
			if (!(check_both_nan(exp[i], act[i]))) {
				ASSERT_DOUBLE_EQ(exp[i], act[i]) << (get_details ? get_details(i) : "");
			}
		}
	}
	else {
		for (size_t i = 0; i < T::size(); i++) {
			ASSERT_EQ(exp[i], act[i]) << (get_details ? get_details(i) : "");
		}
	}
} 

template<typename T>
std::enable_if_t<std::is_integral<T>::value, void>
initialize_values(T* vec, size_t size_v, T start = std::numeric_limits<T>::min(), T end = std::numeric_limits<T>::max()) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<T> dis(start, end);
	for (size_t i = 0; i < size_v; i++) {
		vec[i] = dis(gen);
	}
}

template<typename T>
std::enable_if_t<std::is_floating_point<T>::value, void>
initialize_values(T* vec, size_t size_v, T start = std::numeric_limits<T>::lowest(), T end = std::numeric_limits<T>::max()) {
	//pseudo uniform real values over entire range  is not possible.
	//and undefined too, so we will get from two sets [-min ,0) and [0,max)
	//randomly
	//round some part of floats
	size_t round_index_start = size_v / 10;
	size_t round_index_end = size_v * 3 / 10;
	std::random_device rd;
	std::mt19937 gen(rd());
	bool use_two = start < 0;
	std::uniform_int_distribution<int> change(0, 1);
	std::uniform_real_distribution<T> dis_negative(start, 0);
	std::uniform_real_distribution<T> dis_positive(0, end);
	bool sign = false;
	for (size_t i = 0; i < size_v; i++) {
		T a = {};
		if (use_two) {
			if (change(gen)) {
				a = dis_negative(gen);
			}
			else {
				a = dis_positive(gen);
			}
		}
		else {
			a = dis_positive(gen);
		}
		if (i >= round_index_start && i < round_index_end) {
			T b = std::round(a);
			//only if rounded value within range 
			if (b >= start && b < end) a = b;
		}
		vec[i] = a;
	}
}

template<typename T>
void initialize_default_specials(T* vec, size_t size_v) {
	CACHE_ALIGN T  vec_spec[] = {
		std::numeric_limits<T>::max(),
		std::numeric_limits<T>::max() - 1,
		std::numeric_limits<T>::max() - 4,
		0,11,-117,117,-71,71,
		std::numeric_limits<T>::max() - 14,
		std::numeric_limits<T>::max() - 29,
		std::numeric_limits<T>::max() - 49,
		std::numeric_limits<T>::max() - 75,
		std::numeric_limits<T>::max() - 99,
		0,1,2,3,4,5,87,-87,-93,93,
		std::numeric_limits<T>::lowest() + 1,
		std::numeric_limits<T>::lowest() + 32,
		std::numeric_limits<T>::lowest() + 66,
		std::numeric_limits<T>::min(),
		std::numeric_limits<T>::min() + 1,
		std::numeric_limits<T>::min() + 4,
		0,44,23,-23,111,
		std::numeric_limits<T>::min() + 14,
		std::numeric_limits<T>::min() + 29,
		std::numeric_limits<T>::min() + 49,
		std::numeric_limits<T>::min() + 75,
		std::numeric_limits<T>::min() + 99,
		0,0,0,0,std::numeric_limits<T>::lowest(),
		std::numeric_limits<T>::max() ,0
	};

	//first we will copy exactly
	for (size_t i = 0; i < std::min(size(vec_spec), size_v); i++) {
		vec[i] = vec_spec[i];
	}

	//next we will choose randomly
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<size_t> dis(0, size(vec_spec));

	for (size_t i = std::min(size(vec_spec), size_v); i < size_v; i++) {
		size_t ind = dis(gen);
		vec[i] = vec_spec[ind];
	}


}


template< typename T, size_t arg_count = 1, typename Op1, typename Op2, typename Filter = nullptr_t>
std::enable_if_t<(arg_count >= 1), void>
test_array(
	std::string test_name,
	Op1 expected_f,
	Op2 actual_f, const TestingCase<ValueType<T>>& test_case, Filter filter = {}) {

	using vec_type = T;
	using VT = ValueType<T>;
	constexpr size_t el_count = vec_type::size();
	constexpr size_t size_v = 2048 * el_count;
	constexpr size_t size_specials = (arg_count == 1) ? 32 * el_count : size_v;
	CACHE_ALIGN VT  vals[size_v * arg_count];
	VT* ptr = vals;
	//return detailed failure
	bool bitwise = test_case.isBitwise();
	std::function<std::string(int i)> detail = [&test_name, &ptr, bitwise, size_v](int i) {
		std::stringstream stream;
		using bit_rep = BitValueType<T>;
		stream << "::" << test_name << "\n";
		stream << "{ ";
		if (bitwise) stream << std::hex;
		for (size_t ii = 0; ii < arg_count; ii++) {
			auto local = ptr[ii * size_v + i];
			stream << (bitwise ? bit_cast<bit_rep>(local) : local);
			stream << ((ii < arg_count - 1) ? ", " : " ");
		}
		stream << "}";
		return stream.str();
	};


	VT default_start = std::is_floating_point<VT>::value ? std::numeric_limits<VT>::lowest() : std::numeric_limits<VT>::min();
	VT default_end = std::numeric_limits<VT>::max();
	for (const CheckWithinDomains<VT>& dmn : test_case.getDomains()) {

		size_t dmn_argc = dmn.ArgsDomain.size();
		for (size_t trial = 0; trial < test_case.getTrialCount(); trial++) {
			//lets random initialize  the whole array
			for (size_t i = 0; i < arg_count; i++) {
				VT start = i < dmn_argc ? dmn.ArgsDomain[i].start : default_start;
				VT end = i < dmn_argc ? dmn.ArgsDomain[i].end : default_end;
				ptr = vals + i * size_v;
				initialize_values(ptr, size_v, start, end);
			}
			ptr = vals;
			call_filter<arg_count>(filter, ptr, size_v);
			//test
			for (size_t i = 0; i < size_v; i += el_count) {
				ptr = vals + i;
				AssertVec256(call_op<arg_count>(expected_f, ptr, size_v), call_op<arg_count>(actual_f, ptr, size_v), detail, bitwise,dmn.CheckWithAcceptance,dmn.AcceptedError);
			}
		}//trial

	}
	if (test_case.checkDefaultSpecials()) {
		//fill the whole array with default specials
		for (size_t i = 0; i < arg_count; i++) {
			ptr = vals + i * size_v;
			initialize_default_specials(ptr, size_v);
		}
		ptr = vals;
		call_filter<arg_count>(filter, ptr, size_v);
		//as special values are small count, we will test small range
		for (size_t i = 0; i < size_specials; i += el_count) {
			ptr = vals + i;
			AssertVec256(call_op<arg_count>(expected_f, ptr, size_v), call_op<arg_count>(actual_f, ptr, size_v), detail, bitwise);
		}
	}
	if (test_case.checkNansAndInfinities()) {
		for (size_t ii = 0; ii < el_count; ii++) {
			vals[ii] = std::numeric_limits<VT>::quiet_NaN();
		}
		ptr = vals;
		AssertVec256(call_op<arg_count>(expected_f, ptr, size_v), call_op<arg_count>(actual_f, ptr, size_v), detail, bitwise);
		//infinity
		for (size_t ii = 0; ii < el_count; ii++) {
			vals[ii] = std::numeric_limits<VT>::infinity();
		}
		ptr = vals;
		AssertVec256(call_op<arg_count>(expected_f, ptr, size_v), call_op<arg_count>(actual_f, ptr, size_v), detail, bitwise);
	}
}

template< typename T, size_t arg_count = 1, typename Op1, typename Op2, typename Filter = nullptr_t>
std::enable_if_t<(arg_count >= 1), void>
test_array(
	std::string test_name,
	Op1 expected_f,
	Op2 actual_f, bool bitwise = false,
	Filter filter = {}, bool allow_specials = false, bool test_nan_inf = false, size_t trials = 1) {
	TestingCase< ValueType<T>> test_case = TestingCase<ValueType<T>>::getBuilder()
		.set(bitwise, allow_specials, test_nan_inf)
		.addDomain(CheckWithinDomains<ValueType<T>>{})
		.setTrialCount(trials)
		;
	test_array<T, arg_count, Op1, Op2, Filter>(test_name, expected_f, actual_f, test_case, filter);
}


template<typename T, typename UnaryOp>
T func(UnaryOp call, ValueType<T>* vals) {
	CACHE_ALIGN ValueType<T> ret[T::size()];
	for (size_t i = 0; i < T::size(); i++) {
		ret[i] = call(vals[i]);
	}
	return T::loadu(ret);
}

template<typename T, typename BinaryOp>
T func_bi(BinaryOp call, ValueType<T>* f, ValueType<T>* s) {
	CACHE_ALIGN ValueType<T> ret[T::size()];
	for (size_t i = 0; i < T::size(); i++) {
		ret[i] = call(f[i], s[i]);
	}
	return T::loadu(ret);
}

template<typename T, typename BinaryOp>
T func_bitbi(BinaryOp call, ValueType<T>* f, ValueType<T>* s) {
	using bit_rep = BitValueType<T>;
	CACHE_ALIGN bit_rep ret[T::size()];
	for (size_t i = 0; i < T::size(); i++) {
		ret[i] = call(bit_cast<bit_rep>(f[i]), bit_cast<bit_rep>(s[i]));
	}
	//memory
	return T::loadu(ret);
}

template<typename T, typename BinaryOp>
T func_cmpbi(BinaryOp call, ValueType<T>* f, ValueType<T>* s) {
	using bit_rep = BitValueType<T>;
	constexpr bit_rep mask = std::numeric_limits<bit_rep>::max();
	CACHE_ALIGN bit_rep ret[T::size()];
	for (size_t i = 0; i < T::size(); i++) {
		ret[i] = call(f[i], s[i]) ? mask : 0;
	}
	//memory
	return T::loadu(ret);
}

template<typename T, typename TernaryOp>
T func_tri(TernaryOp call, ValueType<T>* f, ValueType<T>* s, ValueType<T>* th) {
	CACHE_ALIGN ValueType<T> ret[T::size()];
	for (size_t i = 0; i < T::size(); i++) {
		ret[i] = call(f[i], s[i], th[i]);
	}
	return T::loadu(ret);
}