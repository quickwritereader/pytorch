#include <ATen/test/Vec256Test.h>
namespace {

#if GTEST_HAS_TYPED_TEST

    template <typename T>
    class Memory : public ::testing::Test {};

    template <typename T>
    class Arithmetics : public ::testing::Test {};

    template <typename T>
    class Comparison : public ::testing::Test {};

    template <typename T>
    class Bitwise : public ::testing::Test {};

    template <typename T>
    class MinMax : public ::testing::Test {};

    template <typename T>
    class Interleave : public ::testing::Test {};

    template <typename T>
    class SignManipulation : public ::testing::Test {};

    template <typename T>
    class Rounding : public ::testing::Test {};

    template <typename T>
    class SqrtAndReciprocal : public ::testing::Test {};

    template <typename T>
    class SqrtAndReciprocalReal : public ::testing::Test {};

    template <typename T>
    class Trigonometric : public ::testing::Test {};

    template <typename T>
    class ErrorFunctions : public ::testing::Test {};

    template <typename T>
    class Exponents : public ::testing::Test {};

    template <typename T>
    class Hyperbolic : public ::testing::Test {};

    template <typename T>
    class InverseTrigonometric : public ::testing::Test {};

    template <typename T>
    class InverseTrigonometricReal : public ::testing::Test {};

    template <typename T>
    class LGamma : public ::testing::Test {};

    template <typename T>
    class Logarithm : public ::testing::Test {};

    template <typename T>
    class LogarithmReals : public ::testing::Test {};

    template <typename T>
    class Pow : public ::testing::Test {};

    template <typename T>
    class BitwiseFloatsAdditional : public ::testing::Test {};

    template <typename T>
    class BitwiseFloatsAdditional2 : public ::testing::Test {};

    template <typename T>
    class ComplexTests : public ::testing::Test {};

    using RealFloatTestedTypes = ::testing::Types<vfloat, vdouble>;
    using FloatTestedTypes = ::testing::Types<vfloat, vdouble, vcomplex, vcomplexDbl>;
    using ALLTestedTypes = ::testing::Types<vfloat, vdouble, vcomplex, vlong, vint, vshort,
        vqint8, vquint8, vqint>;
    using QuantTestedTypes = ::testing::Types<vqint8, vquint8, vqint>;
    using RealFloatIntTestedTypes =
        ::testing::Types<vfloat, vdouble, vlong, vint, vshort>;
    using FloatIntTestedTypes =
        ::testing::Types<vfloat, vdouble, vcomplex, vcomplexDbl, vlong, vint, vshort>;
    using SingleFloat = ::testing::Types<vfloat>;
    using ComplexTypes = ::testing::Types<  vcomplex, vcomplexDbl>;

    TYPED_TEST_CASE(Memory, ALLTestedTypes);

    TYPED_TEST_CASE(Arithmetics, FloatIntTestedTypes);

    TYPED_TEST_CASE(Comparison, RealFloatIntTestedTypes);

    TYPED_TEST_CASE(Bitwise, FloatIntTestedTypes);

    TYPED_TEST_CASE(MinMax, RealFloatIntTestedTypes);

    TYPED_TEST_CASE(Interleave, RealFloatIntTestedTypes);

    TYPED_TEST_CASE(SignManipulation, FloatIntTestedTypes);

    TYPED_TEST_CASE(Rounding, RealFloatTestedTypes);

    TYPED_TEST_CASE(SqrtAndReciprocal, FloatTestedTypes);

    TYPED_TEST_CASE(SqrtAndReciprocalReal, RealFloatTestedTypes);

    TYPED_TEST_CASE(Trigonometric, RealFloatTestedTypes);

    TYPED_TEST_CASE(ErrorFunctions, RealFloatTestedTypes);

    TYPED_TEST_CASE(Exponents, RealFloatTestedTypes);

    TYPED_TEST_CASE(Hyperbolic, RealFloatTestedTypes);

    TYPED_TEST_CASE(InverseTrigonometricReal, RealFloatTestedTypes);

    TYPED_TEST_CASE(InverseTrigonometric, FloatTestedTypes);

    TYPED_TEST_CASE(LGamma, RealFloatTestedTypes);

    TYPED_TEST_CASE(Logarithm, FloatTestedTypes);
    TYPED_TEST_CASE(LogarithmReals, RealFloatTestedTypes);

    TYPED_TEST_CASE(Pow, RealFloatTestedTypes);

    TYPED_TEST_CASE(BitwiseFloatsAdditional, RealFloatTestedTypes);
    TYPED_TEST_CASE(BitwiseFloatsAdditional2, FloatTestedTypes);

    TYPED_TEST(Memory, UnAlignedLoadStore) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        constexpr size_t b_size = vec_type::size() * sizeof(VT);
        CACHE_ALIGN unsigned char ref_storage[128 * b_size];
        CACHE_ALIGN unsigned char storage[128 * b_size];
        // fill with gibberish
        for (auto& x : ref_storage) {
            x = std::rand() % 255;
        }
        // test counted load stores
#ifdef __VSX__
        for (int i = 1; i < 2 * vec_type::size(); i++) {
            vec_type v = vec_type::loadu(ref_storage, i);
            v.store(storage);
            size_t count = std::min(i * sizeof(VT), b_size);
            bool cmp = (std::memcmp(ref_storage, storage, count) == 0);
            ASSERT_TRUE(cmp) << "failure count: " << i;
            // clear storage
            std::memset(storage, 0, b_size);
        }
#endif
        // testing unaligned load store
        for (size_t offset = 0; offset < b_size; offset += 1) {
            unsigned char* p1 = ref_storage + offset;
            unsigned char* p2 = storage + offset;
            for (; p1 + b_size <= std::end(ref_storage); p1 += b_size, p2 += b_size) {
                vec_type v = vec_type::loadu(p1);
                v.store(p2);
            }
            size_t written = p1 - ref_storage - offset;
            bool cmp =
                (std::memcmp(ref_storage + offset, storage + offset, written) == 0);
            ASSERT_TRUE(cmp) << "failure at unaligned offset: " << offset;
            // clear storage
            std::memset(storage, 0, sizeof storage);
        }
    }



    TYPED_TEST(SignManipulation, Absolute) {
        using vec_type = TypeParam; 
        test_unary<vec_type>(
            "absolute", RESOLVE_OVERLOAD(std::abs),
            [](vec_type v) { return v.abs(); }, false,
            RESOLVE_OVERLOAD(filter_int_minimum), false);
    }

    TYPED_TEST(SignManipulation, Negate) {
        using vec_type = TypeParam;
        // negate overflows for minimum on int and long
        test_unary<vec_type>(
            "negate", std::negate<ValueType<vec_type>>(),
            [](vec_type v) { return v.neg(); }, false,
            RESOLVE_OVERLOAD(filter_int_minimum), false);
    }


    TYPED_TEST(Rounding, Round) { 
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "round", RESOLVE_OVERLOAD(local_round),
            [](vec_type v) { return v.round(); });
    }

    TYPED_TEST(Rounding, Ceil) {
        using vec_type = TypeParam;
        using UVT = UvalueType<TypeParam>; 
        test_unary<vec_type>(
            "ceil", RESOLVE_OVERLOAD(std::ceil),
            [](vec_type v) { return v.ceil(); });
    }

    TYPED_TEST(Rounding, Floor) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "floor", RESOLVE_OVERLOAD(std::floor),
            [](vec_type v) { return v.floor(); });
    }

    TYPED_TEST(Rounding, Trunc) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "trunc", RESOLVE_OVERLOAD(std::trunc),
            [](vec_type v) { return v.trunc(); });
    }

    TYPED_TEST(SqrtAndReciprocal, Sqrt) {
        //pytorch complex sqrt precision differs from std 
        //we will check for error under 0.001 
        using vec_type = TypeParam;
        using UVT = UvalueType<TypeParam>;
        auto test_case =
            TestingCase<UVT>::getBuilder()
            .addDomain(CheckWithinDomains<UVT>{ {{-100, 100}}, true, 1.e-3f})
            .setTrialCount(200);
        test_unary<vec_type>(
            "sqrt", RESOLVE_OVERLOAD(std::sqrt),
            [](vec_type v) { return v.sqrt(); }, test_case);
    }

    TYPED_TEST(SqrtAndReciprocalReal, Sqrt) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "sqrt", RESOLVE_OVERLOAD(std::sqrt),
            [](vec_type v) { return v.sqrt(); }, false, {}, true);
    }

    TYPED_TEST(SqrtAndReciprocalReal, RSqrt) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "rsqrt", rsqrt<ValueType<vec_type>>,
            [](vec_type v) { return v.rsqrt(); }, false,
            RESOLVE_OVERLOAD(filter_zero));
    }

    TYPED_TEST(SqrtAndReciprocalReal, Reciprocal) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "reciprocal",
            reciprocal<ValueType<vec_type>>,
            [](vec_type v) { return v.reciprocal(); }, false,
            RESOLVE_OVERLOAD(filter_zero));
    }

    TYPED_TEST(Trigonometric, Sin) {
        using vec_type = TypeParam;
        using UVT = UvalueType<TypeParam>;
        auto test_case =
            TestingCase<UVT>::getBuilder()
            .addDomain(CheckWithinDomains<UVT>{ { {-4096, 4096}}, true, 1.2e-7f})
            .addDomain(CheckWithinDomains<UVT>{ { {-8192, 8192}}, true, 3.0e-7f})
            .setTrialCount(8000);
        test_unary<vec_type>(
            "sin",
            RESOLVE_OVERLOAD(std::sin),
            [](vec_type v) { return v.sin(); }, test_case);
    }

    TYPED_TEST(Trigonometric, Cos) {
        using vec_type = TypeParam;
        using UVT = UvalueType<TypeParam>;
        auto test_case =
            TestingCase<UVT>::getBuilder()
            .addDomain(CheckWithinDomains<UVT>{ { {-4096, 4096}}, true, 1.2e-7f})
            .addDomain(CheckWithinDomains<UVT>{ { {-8192, 8192}}, true, 3.0e-7f})
            .setTrialCount(8000);
        test_unary<vec_type>(
            "cos",
            RESOLVE_OVERLOAD(std::cos),
            [](vec_type v) { return v.cos(); }, test_case);
    }

    TYPED_TEST(Trigonometric, Tan) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "tan",
            RESOLVE_OVERLOAD(std::tan),
            [](vec_type v) { return v.tan(); });
    }

    TYPED_TEST(Hyperbolic, Tanh) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "tanH", RESOLVE_OVERLOAD(std::tanh),
            [](vec_type v) { return v.tanh(); });
    }

    TYPED_TEST(Hyperbolic, Sinh) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "sinh", RESOLVE_OVERLOAD(std::sinh),
            [](vec_type v) { return v.sinh(); });
    }

    TYPED_TEST(Hyperbolic, Cosh) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "cosh", RESOLVE_OVERLOAD(std::cosh),
            [](vec_type v) { return v.cosh(); });
    }

    TYPED_TEST(InverseTrigonometric, Asin) {
        bool checkRelativeErr = is_complex<ValueType<TypeParam>>();
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "asin", RESOLVE_OVERLOAD(std::asin),
            [](vec_type v) { return v.asin(); }, false, {}, checkRelativeErr);
    }

    TYPED_TEST(InverseTrigonometric, ACos) {
        bool checkRelativeErr = is_complex<ValueType<TypeParam>>();
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "acos", RESOLVE_OVERLOAD(std::acos),
            [](vec_type v) { return v.acos(); }, false, {}, checkRelativeErr);
    }

    TYPED_TEST(InverseTrigonometric, ATan) {
        bool checkRelativeErr = is_complex<ValueType<TypeParam>>();
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "atan", RESOLVE_OVERLOAD(std::atan),
            [](vec_type v) { return v.atan(); }, false, {}, checkRelativeErr);
    }

    TYPED_TEST(Logarithm, Log) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "log", RESOLVE_OVERLOAD(std::log),
            [](const vec_type& v) { return v.log(); });
    }

    TYPED_TEST(Logarithm, Log2) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "log2", RESOLVE_OVERLOAD(local_log2),
            [](const vec_type& v) { return v.log2(); });
    }

    TYPED_TEST(Logarithm, Log10) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "log10", RESOLVE_OVERLOAD(std::log10),
            [](const vec_type& v) { return v.log10(); });
    }

    TYPED_TEST(LogarithmReals, Log1p) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "log1p", RESOLVE_OVERLOAD(std::log1p),
            [](const vec_type& v) { return v.log1p(); }, false, {},true);
    }

    TYPED_TEST(Exponents, Exp) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "exp", RESOLVE_OVERLOAD(std::exp),
            [](const vec_type& v) { return v.exp(); });
    }

    TYPED_TEST(Exponents, Expm1) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "expm1", RESOLVE_OVERLOAD(std::expm1), 
            [](const vec_type& v) { return v.expm1(); }, false, {}, true);
    }

    TYPED_TEST(ErrorFunctions, Erf) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "erf",
            RESOLVE_OVERLOAD(std::erf),
            [](const vec_type& v) { return v.erf(); });
    }

    TYPED_TEST(ErrorFunctions, Erfc) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "erfc", RESOLVE_OVERLOAD(std::erfc), 
            [](const vec_type& v) { return v.erfc(); });
    }

    TYPED_TEST(ErrorFunctions, Erfinv) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "erfinv", RESOLVE_OVERLOAD(calc_erfinv),
            [](const vec_type& v) { return v.erfinv(); });
    }

    TYPED_TEST(LGamma, LGamma) {
        using vec_type = TypeParam;
        test_unary<vec_type>(
            "lgamma", RESOLVE_OVERLOAD(std::lgamma),
            [](vec_type v) { return v.lgamma(); });
    }

    TYPED_TEST(InverseTrigonometricReal, ATan2) {
        using vec_type = TypeParam;
        test_binary<vec_type>(
            "atan2", RESOLVE_OVERLOAD(std::atan2),
            [](vec_type v0, vec_type v1) {
                return v0.atan2(v1);
            });
    }


    TYPED_TEST(Pow, Pow) {
        using vec_type = TypeParam;
        test_binary<vec_type>(
            "pow", RESOLVE_OVERLOAD(std::pow),
            [](vec_type v0, vec_type v1) { return v0.pow(v1); });
    }

    TYPED_TEST(Interleave, Interleave) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        constexpr size_t N = vec_type::size() * 2;
        CACHE_ALIGN VT vals[N];
        CACHE_ALIGN VT interleaved[N]; 
        ValueGen<VT> generator;
        for (VT& v : vals) {
            v = generator.get();
        }
        copy_interleave(vals, interleaved);
        auto a = vec_type::loadu(vals);
        auto b = vec_type::loadu(vals + N / 2);
        auto cc = interleave2(a, b);
        size_t start_i = 0;
        std::function<std::string(int i)> detail = [&start_i](int i) {
            std::stringstream stream;
            stream << "::interleave\nfail index " << start_i + i;
            return stream.str();
        };
        AssertVec256(std::get<0>(cc), vec_type::loadu(interleaved), detail, true);
        start_i = N / 2;
        AssertVec256(std::get<1>(cc), vec_type::loadu(interleaved + start_i), detail,
            true);
    }

    TYPED_TEST(Interleave, DeInterleave) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        constexpr size_t N = vec_type::size() * 2;
        CACHE_ALIGN VT vals[N];
        CACHE_ALIGN VT interleaved[N];
        ValueGen<VT> generator;
        for (VT& v : vals) {
            v = generator.get();
        }
        copy_interleave(vals, interleaved);
        // test interleaved with vals this time
        auto a = vec_type::loadu(interleaved);
        auto b = vec_type::loadu(interleaved + N / 2);
        auto cc = deinterleave2(a, b);
        size_t start_i = 0;
        std::function<std::string(int i)> detail = [&start_i](int i) {
            std::stringstream stream;
            stream << "::deinterleave\nfail index " << start_i + i;
            return stream.str();
        };
        AssertVec256(std::get<0>(cc), vec_type::loadu(vals), detail, true);
        start_i = N / 2;
        AssertVec256(std::get<1>(cc), vec_type::loadu(vals + start_i), detail, true);
    }

    TYPED_TEST(Arithmetics, Plus) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "plus", std::plus<VT>(),
            [](const vec_type &v0, const vec_type& v1) -> vec_type {
                return v0 + v1;
            },
                false, RESOLVE_OVERLOAD(filter_add_overflow), false);
    }

    TYPED_TEST(Arithmetics, Minus) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "minus", std::minus<VT>(),
            [](const vec_type& v0, const vec_type& v1) -> vec_type {
                return v0 - v1;
            },
            false, RESOLVE_OVERLOAD(filter_minus_overflow), false);
    }


    TYPED_TEST(Arithmetics, Multiplication) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        if (is_complex<VT>()) {
            using UVT = UvalueType<TypeParam>;
            auto test_case =
                TestingCase<UVT>::getBuilder()
                .addDomain(CheckWithinDomains<UVT>{ { DomainRange<UVT>{(UVT)-100, (UVT)100}, DomainRange<UVT>{(UVT)-100, (UVT)100}}, true, (UVT)(1.e-5) });
            test_binary<vec_type>(
                "mult",
                std::multiplies<VT>(),
                [](const vec_type& v0, const vec_type& v1) { return v0 * v1; },
                test_case, RESOLVE_OVERLOAD(filter_mult_overflow));
        }
        else {
            test_binary<vec_type>(
                "mult", std::multiplies<VT>(),
                [](const vec_type& v0, const vec_type& v1) { return v0 * v1; },
                false, RESOLVE_OVERLOAD(filter_mult_overflow), false);
        }
    }



    TYPED_TEST(Arithmetics, Division) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        //for complex we will use small range and absError against std implementation
        //because inside our implementation we are using the same type and multiplication easily can become inf
        // try for example std::complex<float>(1.7852e+38,1.65523e+38)/std::complex<float>(1.74044e+38,1.57524e+38)
        if (is_complex<VT>()) {
            using UVT = UvalueType<TypeParam>;
            auto test_case =
                TestingCase<UVT>::getBuilder()
                .addDomain(CheckWithinDomains<UVT>{ { DomainRange<UVT>{(UVT)-10, (UVT)10}, DomainRange<UVT>{(UVT)-10, (UVT)10}}, true,  (UVT)(1.e-5) });
            test_binary<vec_type>(
                "division",
                std::divides<VT>(),
                [](const vec_type& v0, const vec_type& v1) { return v0 / v1; },
                test_case, RESOLVE_OVERLOAD(filter_div_ub));
        }
        else {

            test_binary<vec_type>(
                "division",
                std::divides<VT>(),
                [](const vec_type& v0, const vec_type& v1) { return v0 / v1; },
                false, RESOLVE_OVERLOAD(filter_div_ub), false);
        }
    }



    TYPED_TEST(Bitwise, BitAnd) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "bit_and",  local_and<VT>,
            [](const vec_type& v0, const vec_type& v1) { return v0 & v1; }, true);
    }
    
    TYPED_TEST(Bitwise, BitOr) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "bit_or", local_or<VT>,
            [](const vec_type& v0, const vec_type& v1) { return v0 | v1; }, true);
    }

    TYPED_TEST(Bitwise, BitXor) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "bit_xor", local_xor<VT>,
            [](const vec_type& v0, const vec_type& v1) { return v0 ^ v1; }, true);
    }


    TYPED_TEST(Comparison, Equal) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "==",
            [](const VT& v1, const VT& v2) {return func_cmp(std::equal_to<VT>(), v1, v2); },
            [](const vec_type& v0, const vec_type& v1) { return v0 == v1; },
            true, {});
    }

    TYPED_TEST(Comparison, NotEqual) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "!=",
            [](const VT& v1, const VT& v2) {return func_cmp(std::not_equal_to<VT>(), v1, v2); },
            [](const vec_type& v0, const vec_type& v1) { return v0 != v1; },
            true, {});
    }

    TYPED_TEST(Comparison, Greater) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            ">",
            [](const VT& v1, const VT& v2) {return func_cmp(std::greater<VT>(), v1, v2); },
            [](const vec_type& v0, const vec_type& v1) { return v0 > v1; },
            true, {});
    }

    TYPED_TEST(Comparison, Less) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "<",
            [](const VT& v1, const VT& v2) {return func_cmp(std::less<VT>(), v1, v2); },
            [](const vec_type& v0, const vec_type& v1) { return v0 < v1; },
            true, {});
    }

    TYPED_TEST(Comparison, GreaterEqual) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            ">=",
            [](const VT& v1, const VT& v2) {return func_cmp(std::greater_equal<VT>(), v1, v2); },
            [](const vec_type& v0, const vec_type& v1) { return v0 >= v1; },
            true, {});
    }

    TYPED_TEST(Comparison, LessEqual) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "<=",
            [](const VT& v1, const VT& v2) {return func_cmp(std::less_equal<VT>(), v1, v2); },
            [](const vec_type& v0, const vec_type& v1) { return v0 <= v1; },
            true, {});
    }

    TYPED_TEST(MinMax, Minimum) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "minimum",
            minimum<VT>,
            [](const vec_type& v0, const vec_type& v1) {
                return minimum(v0, v1);
            });
    }

    TYPED_TEST(MinMax, Maximum) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "maximum",
            maximum<VT>,
            [](const vec_type& v0, const vec_type& v1) {
                return maximum(v0, v1);
            });
    }

    TYPED_TEST(MinMax, ClampMin) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "clamp min",
            clamp_min<VT>,
            [](const vec_type& v0, const vec_type& v1) {
                return clamp_min(v0, v1);
            });
    }

    TYPED_TEST(MinMax, ClampMax) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_binary<vec_type>(
            "clamp max",
            clamp_max<VT>,
            [](const vec_type& v0, const vec_type& v1) {
                return clamp_max(v0, v1);
            });
    }

 
    TYPED_TEST(MinMax, Clamp) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_ternary<vec_type>(
            "clamp", clamp<VT>,
            [](const vec_type& v0, const vec_type& v1, const vec_type& v2) {
                return clamp(v0, v1, v2);
            },
                false, RESOLVE_OVERLOAD(filter_clamp));
    }
     

    TYPED_TEST(BitwiseFloatsAdditional, ZeroMask) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        CACHE_ALIGN VT test_vals[vec_type::size()];
        //all sets will be within 0  2^(n-1)
        auto power_sets = 1 << (vec_type::size());
        for (int expected = 0; expected < power_sets; expected++) {
            // generate test_val based on expected
            for (int i = 0; i < vec_type::size(); ++i)
            {
                if (expected & (1 << i)) {
                    test_vals[i] = (VT)0;
                }
                else {
                    test_vals[i] = (VT)0.897;
                }
            }

            int actual = vec_type::loadu(test_vals).zero_mask();
            ASSERT_EQ(expected, actual) << std::hex << actual << ";" << expected;

        } //

    }

    template<typename vec_type, typename VT, int64_t mask>
    typename std::enable_if_t<(mask < 0 || mask> 255), void>
        test_blend(VT expected_val[vec_type::size()], VT a[vec_type::size()], VT b[vec_type::size()])
    {
    }

    template<typename vec_type, typename VT, int64_t mask>
    typename std::enable_if_t<(mask >= 0 && mask <= 255), void>
        test_blend(VT expected_val[vec_type::size()], VT a[vec_type::size()], VT b[vec_type::size()])
    {
        //generate expected_val
        int64_t m = mask;
        for (int64_t i = 0; i < vec_type::size(); i++) {
            if (m & 0x01) {
                expected_val[i] = b[i];
            }
            else {
                expected_val[i] = a[i];
            }
            m = m >> 1;
        }
        //test with blend
        auto vec_a = vec_type::loadu(a);
        auto vec_b = vec_type::loadu(b);
        auto expected = vec_type::loadu(expected_val);
        auto actual = vec_type::template blend<mask>(vec_a, vec_b);
        std::function<std::string(int i)> detail = [](int i) {
            std::stringstream stream;
            stream << "mask " << mask << " index:" << i;
            return stream.str();
        };
        AssertVec256(expected, actual, detail);
        test_blend<vec_type, VT, mask - 1>(expected_val, a, b);
    }


    template<typename T, int N>
    void blend_init(T(&a)[N], T(&b)[N]) {
        a[0] = (T)1.0;
        b[0] = a[0] + (T)N;
        for (int i = 1; i < N; i++) {
            a[i] = a[i - 1] + (T)(1.0);
            b[i] = b[i - 1] + (T)(1.0);
        }
    }

    template<>
    void blend_init<std::complex<float>, 4>(std::complex<float>(&a)[4], std::complex<float>(&b)[4]) {
        auto add = std::complex<float>(1., 100.);
        a[0] = std::complex<float>(1., 100.);
        b[0] = std::complex<float>(5., 1000.);
        for (int i = 1; i < 4; i++) {
            a[i] = a[i - 1] + add;
            b[i] = b[i - 1] + add;
        }
    }

    template<>
    void blend_init<std::complex<double>, 2>(std::complex<double>(&a)[2], std::complex<double>(&b)[2]) {
        auto add = std::complex<double>(1.0, 100.0);
        a[0] = std::complex<double>(1.0, 100.0);
        b[0] = std::complex<double>(3.0, 1000.0);
        a[1] = a[0] + add;
        b[1] = b[0] + add;
    }

    TYPED_TEST(BitwiseFloatsAdditional2, Blend) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        CACHE_ALIGN VT a[vec_type::size()];
        CACHE_ALIGN VT b[vec_type::size()];
        CACHE_ALIGN VT expected_val[vec_type::size()];
        blend_init(a, b);
        constexpr int64_t power_sets = 1 << (vec_type::size());
        test_blend<vec_type, VT, power_sets - 1>(expected_val, a, b);
    }




    TEST(ComplexTests, TestComplexFloatImagRealConj) {
        //vcomplex a = { std::complex<float>(1,2),std::complex<float>(3,4) ,std::complex<float>(5,6) ,std::complex<float>(6,7) };
        float aa[] = { 1.5488e-28,2.5488e-28,3.5488e-28,4.5488e-28,5.5488e-28,6.5488e-28,7.5488e-28,8.5488e-28 };
        float exp[] = { aa[0],0,aa[2],0,aa[4],0,aa[6],0 }; 
        float exp3[] = { aa[1],0,aa[3],0,aa[5],0,aa[7],0 };
        float exp4[] = { 1.5488e-28, -2.5488e-28,3.5488e-28,-4.5488e-28,5.5488e-28,-6.5488e-28,7.5488e-28,-8.5488e-28 };
        auto a = vcomplex::loadu(aa); 
        auto actual1 = a.real(); 
        auto actual3 = a.imag();
        auto actual4 = a.conj(); 
        std::cout << actual1 << std::endl; 
        std::cout << actual3 << std::endl;
        std::cout << actual4 << std::endl; 
        auto expected1 = vcomplex::loadu(exp); 
        auto expected3 = vcomplex::loadu(exp3);
        auto expected4 = vcomplex::loadu(exp4); 
        AssertVec256(expected1, actual1); 
        AssertVec256(expected3, actual3);
        AssertVec256(expected4, actual4);
    }

   

#endif
}  // namespace