#if 0
#include <Vec256Test.h>
namespace {

#if GTEST_HAS_TYPED_TEST

    template<typename T>
    class Memory : public ::testing::Test {
    };

    template<typename T>
    class Arithmetics : public ::testing::Test {
    };

    template<typename T>
    class Comparison : public ::testing::Test {
    };

    template<typename T>
    class Bitwise : public ::testing::Test {
    };

    template<typename T>
    class MinMax : public ::testing::Test {
    };

    template<typename T>
    class Interleave : public ::testing::Test {
    };

    template<typename T>
    class SignManipulation : public ::testing::Test {
    };

    template<typename T>
    class Rounding : public ::testing::Test {
    };

    template<typename T>
    class SqrtAndReciprocal : public ::testing::Test {
    };

    template<typename T>
    class Trigonometric : public ::testing::Test {
    };

    template<typename T>
    class ErrorFunctions : public ::testing::Test {
    };

    template<typename T>
    class Exponents : public ::testing::Test {
    };

    template<typename T>
    class Hyperbolic : public ::testing::Test {
    };

    template<typename T>
    class InverseTrigonometric : public ::testing::Test {
    };

    template<typename T>
    class LGamma : public ::testing::Test {
    };

    template<typename T>
    class Logarithm : public ::testing::Test {
    };

    template<typename T>
    class Pow : public ::testing::Test {
    };
    using FloatTestedTypes = ::testing::Types<vfloat, vdouble >;
    using ALLTestedTypes = ::testing::Types<vfloat, vdouble, vlong, vint, vshort, vqint8, vquint8, vqint  >;
    using FloatTestedTypes = ::testing::Types<vfloat, vdouble >;
    using QuantTestedTypes = ::testing::Types<vqint8, vquint8, vqint >;
    using FloatIntTestedTypes = ::testing::Types<vfloat, vdouble, vlong, vint, vshort >;
    using SingleFloat = ::testing::Types<vfloat>;

    TYPED_TEST_CASE(Memory, ALLTestedTypes);

    TYPED_TEST_CASE(Arithmetics, FloatIntTestedTypes);

    TYPED_TEST_CASE(Comparison, FloatIntTestedTypes);

    TYPED_TEST_CASE(Bitwise, FloatIntTestedTypes);

    TYPED_TEST_CASE(MinMax, FloatIntTestedTypes);

    TYPED_TEST_CASE(Interleave, FloatIntTestedTypes);

    TYPED_TEST_CASE(SignManipulation, FloatIntTestedTypes);

    TYPED_TEST_CASE(Rounding, FloatTestedTypes);

    TYPED_TEST_CASE(SqrtAndReciprocal, FloatTestedTypes);

    TYPED_TEST_CASE(Trigonometric, FloatTestedTypes);

    TYPED_TEST_CASE(ErrorFunctions, FloatTestedTypes);

    TYPED_TEST_CASE(Exponents, FloatTestedTypes);

    TYPED_TEST_CASE(Hyperbolic, FloatTestedTypes);

    TYPED_TEST_CASE(InverseTrigonometric, FloatTestedTypes);

    TYPED_TEST_CASE(LGamma, FloatTestedTypes);

    TYPED_TEST_CASE(Logarithm, FloatTestedTypes);

    TYPED_TEST_CASE(Pow, FloatTestedTypes);

    TYPED_TEST_CASE(Pow, FloatTestedTypes);




    TYPED_TEST(Memory, UnAlignedLoadStore) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        constexpr size_t b_size = vec_type::size() * sizeof(VT);
        CACHE_ALIGN unsigned char  ref_storage[128 * b_size];
        CACHE_ALIGN unsigned char  storage[128 * b_size];
        //fill with gibberish
        for (auto& x : ref_storage) {
            x = std::rand() % 255;
        }
        //test counted load stores
#ifdef __VSX__
        for (size_t i = 1; i < 2 * vec_type::size(); i++) {
            vec_type v = vec_type::loadu(ref_storage, i);
            v.store(storage);
            size_t count = std::min(i * sizeof(VT), b_size);
            bool cmp = (std::memcmp(ref_storage, storage, count) == 0);
            ASSERT_TRUE(cmp) << "failure count: " << i;
            //clear storage
            std::memset(storage, 0, b_size);
        }
#endif
        //testing unaligned load store
        for (size_t offset = 0; offset < b_size; offset += 1)
        {
            unsigned char* p1 = ref_storage + offset;
            unsigned char* p2 = storage + offset;
            for (; p1 + b_size <= std::end(ref_storage); p1 += b_size, p2 += b_size) {
                vec_type v = vec_type::loadu(p1);
                v.store(p2);
            }
            size_t written = p1 - ref_storage - offset;
            bool cmp = (std::memcmp(ref_storage + offset, storage + offset, written) == 0);
            ASSERT_TRUE(cmp) << "failure at unaligned offset: " << offset;
            //clear storage
            std::memset(storage, 0, sizeof storage);
        }
    }

    TYPED_TEST(Arithmetics, Plus) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array< vec_type, BINARY>("plus",
            [](VT* f, VT* s) ->vec_type {
                return func_bi<vec_type>(std::plus<VT>(), f, s);
            },
            [](VT* f, VT* s) ->vec_type {
                return vec_type::loadu(f) + vec_type::loadu(s);
            },
                false, RESOLVE_OVERLOAD(filter_add_overflow), false);
    }
#if 1

    TYPED_TEST(Arithmetics, Minus) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("minus",
            [](VT* f, VT* s) {
                return func_bi<vec_type>(std::minus<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) - vec_type::loadu(s);
            },
                false, RESOLVE_OVERLOAD(filter_minus_overflow), false);
    }

    TYPED_TEST(Arithmetics, Multiplication) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("mult",
            [](VT* f, VT* s) {
                return func_bi<vec_type>(std::multiplies<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) * vec_type::loadu(s);
            },
                false, RESOLVE_OVERLOAD(filter_mult_overflow), false);
    }

    TYPED_TEST(Arithmetics, Division) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("division",
            [](VT* f, VT* s) {
                return func_bi<vec_type>(std::divides<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) / vec_type::loadu(s);
            },
                false, RESOLVE_OVERLOAD(filter_div_ub), false);
    }

    TYPED_TEST(Bitwise, BitAnd) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        using bit_rep = BitValueType<TypeParam>;
        test_array<vec_type, BINARY>("bit_and",
            [](VT* f, VT* s) {
                return func_bitbi<vec_type>(std::bit_and<bit_rep>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) & vec_type::loadu(s);
            });
    }

    TYPED_TEST(Bitwise, BitOr) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        using bit_rep = BitValueType<TypeParam>;
        test_array<vec_type, BINARY>("bit_or",
            [](VT* f, VT* s) {
                return func_bitbi<vec_type>(std::bit_or<bit_rep>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) | vec_type::loadu(s);
            });
    }

    TYPED_TEST(Bitwise, BitXor) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        using bit_rep = BitValueType<TypeParam>;
        test_array<vec_type, BINARY>("bit_xor",
            [](VT* f, VT* s) {
                return func_bitbi<vec_type>(std::bit_xor<bit_rep>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) ^ vec_type::loadu(s);
            });
    }

    TYPED_TEST(Comparison, Equal) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("==",
            [](VT* f, VT* s) {
                return func_cmpbi<vec_type>(std::equal_to<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) == vec_type::loadu(s);
            }
            , true, {}, true);
    }

    TYPED_TEST(Comparison, NotEqual) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("!=",
            [](VT* f, VT* s) {
                return func_cmpbi<vec_type>(std::not_equal_to<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) != vec_type::loadu(s);
            }
            , true, {}, true);
    }

    TYPED_TEST(Comparison, Greater) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>(">",
            [](VT* f, VT* s) {
                return func_cmpbi<vec_type>(std::greater<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) > vec_type::loadu(s);
            }
            , true, {}, true);
    }

    TYPED_TEST(Comparison, Less) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("<",
            [](VT* f, VT* s) {
                return func_cmpbi<vec_type>(std::less<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) < vec_type::loadu(s);
            }
            , true, {}, true);
    }

    TYPED_TEST(Comparison, GreaterEqual) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>(">=",
            [](VT* f, VT* s) {
                return func_cmpbi<vec_type>(std::greater_equal<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) >= vec_type::loadu(s);
            }
            , true, {}, true);
    }

    TYPED_TEST(Comparison, LessEqual) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("<=",
            [](VT* f, VT* s) {
                return func_cmpbi<vec_type>(std::less_equal<VT>(), f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f) <= vec_type::loadu(s);
            }
            , true, {}, true);
    }

    TYPED_TEST(MinMax, Minimum) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("minimum",
            [](VT* f, VT* s) {
                return func_bi<vec_type>(minimum<VT>, f, s);
            },
            [](VT* f, VT* s) {
                return minimum(vec_type::loadu(f), vec_type::loadu(s));
            }
            , false, {}, true);
    }

    TYPED_TEST(MinMax, Maximum) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("maximum",
            [](VT* f, VT* s) {
                return func_bi<vec_type>(maximum<VT>, f, s);
            },
            [](VT* f, VT* s) {
                return maximum(vec_type::loadu(f), vec_type::loadu(s));
            }
            , false, {}, true);
    }

    TYPED_TEST(MinMax, ClampMin) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("clamp min",
            [](VT* f, VT* s) {
                return func_bi<vec_type>(clamp_min<VT>, f, s);
            },
            [](VT* f, VT* s) {
                return clamp_min(vec_type::loadu(f), vec_type::loadu(s));
            }
            , false, {}, true);
    }

    TYPED_TEST(MinMax, ClampMax) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, BINARY>("clamp max",
            [](VT* f, VT* s) {
                return func_bi<vec_type>(clamp_max<VT>, f, s);
            },
            [](VT* f, VT* s) {
                return clamp_max(vec_type::loadu(f), vec_type::loadu(s));
            }
            , false, {}, true);
    }

    TYPED_TEST(MinMax, Clamp) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type, TERNARY>("clamp",
            [](VT* f, VT* s, VT* th) {
                return func_tri<vec_type>(clamp<VT>, f, s, th);
            },
            [](VT* f, VT* s, VT* th) {
                return clamp(vec_type::loadu(f), vec_type::loadu(s), vec_type::loadu(th));
            }
            , false, RESOLVE_OVERLOAD(filter_clamp), true);
    }

    TYPED_TEST(SignManipulation, Absolute) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("absolute",
            [](VT* vals) {
                return func<vec_type>(local_abs<VT>, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).abs();
            }, false, RESOLVE_OVERLOAD(filter_int_minimum), true);
    }

    TYPED_TEST(SignManipulation, Negate) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        //negate overflows for minimum on int and long
        test_array<vec_type>("negate",
            [](VT* vals) {
                return func<vec_type>(std::negate<VT>(), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).neg();
            }
            , false, RESOLVE_OVERLOAD(filter_int_minimum), true);
    }

    TYPED_TEST(Interleave, Interleave) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        constexpr size_t N = vec_type::size() * 2;
        CACHE_ALIGN VT vals[N];
        CACHE_ALIGN VT interleaved[N];
        initialize_values(vals, size(vals));
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
        AssertVec256(std::get<1>(cc), vec_type::loadu(interleaved + start_i), detail, true);
    }

    TYPED_TEST(Interleave, DeInterleave) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        constexpr size_t N = vec_type::size() * 2;
        CACHE_ALIGN VT vals[N];
        CACHE_ALIGN VT interleaved[N];
        initialize_values(vals, size(vals));
        copy_interleave(vals, interleaved);
        //test interleaved with vals this time
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

    TYPED_TEST(Rounding, Round) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("round",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::round), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).round();
            });
    }

    TYPED_TEST(Rounding, Ceil) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("ceil",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::ceil), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).ceil();
            });
    }

    TYPED_TEST(Rounding, Floor) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("floor",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::floor), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).floor();
            });
    }

    TYPED_TEST(Rounding, Trunc) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("trunc",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::trunc), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).trunc();
            });
    }

    TYPED_TEST(SqrtAndReciprocal, Sqrt) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("sqrt",
            [](VT* vals) {
                return func<vec_type>(sqrt<VT>, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).sqrt();
            }, false, {}, true);
    }

    TYPED_TEST(SqrtAndReciprocal, RSqrt) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("rsqrt",
            [](VT* vals) {
                return func<vec_type>(rsqrt<VT>, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).rsqrt();
            }, false, RESOLVE_OVERLOAD(filter_zero), true);
    }

    TYPED_TEST(SqrtAndReciprocal, Reciprocal) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("reciprocal",
            [](VT* vals) {
                return func<vec_type>(reciprocal<VT>, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).reciprocal();
            }, false, RESOLVE_OVERLOAD(filter_zero), true);
    }

    TYPED_TEST(Trigonometric, Sin) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;

        auto test_case = TestingCase<VT>::getBuilder()
            .addDomain(CheckWithinDomains<VT>{ { {-4096, 4096}}, true, 1.2e-7f })
            .addDomain(CheckWithinDomains<VT>{ { {-8192, 8192}}, true, 3.0e-7f })
            .setTrialCount(200);
        test_array<vec_type>("sin",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::sin), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).sin();
            }, test_case);
    }

    TYPED_TEST(Trigonometric, Cos) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        auto test_case = TestingCase<VT>::getBuilder()
            .addDomain(CheckWithinDomains<VT>{ { {-4096, 4096}}, true, 1.2e-7f })
            .addDomain(CheckWithinDomains<VT>{ { {-8192, 8192}}, true, 3.0e-7f })
            .setTrialCount(200);
        test_array<vec_type>("cos",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::cos), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).cos();
            },test_case);
    }

    TYPED_TEST(Trigonometric, Tan) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("tan",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::tan), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).tan();
            });
    }

    TYPED_TEST(Hyperbolic, Tanh) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("tanH",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::tanh), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).tanh();
            });
    }

    TYPED_TEST(Hyperbolic, Sinh) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("sinh",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::sinh), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).sinh();
            });
    }

    TYPED_TEST(Hyperbolic, Cosh) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("cosh",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::cosh), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).cosh();
            });
    }

    TYPED_TEST(InverseTrigonometric, Asin) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("asin",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::asin), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).asin();
            });
    }

    TYPED_TEST(InverseTrigonometric, ACos) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("acos",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::acos), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).acos();
            });
    }

    TYPED_TEST(InverseTrigonometric, ATan) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT) = &std::atan;
        test_array<vec_type>("atan",
            [call](VT* vals) {
                return func<vec_type>(call, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).atan();
            });
    }

    TYPED_TEST(InverseTrigonometric, ATan2) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT, VT) = &std::atan2;
        test_array<vec_type, BINARY>("atan2",
            [call](VT* f, VT* s) {
                return func_bi<vec_type>(call, f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f).atan2(vec_type::loadu(s));
            });
    }

    TYPED_TEST(Logarithm, Log) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT) = &std::log;
        test_array<vec_type>("log",
            [call](VT* vals) {
                return func<vec_type>(call, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).log();
            });
    }

    TYPED_TEST(Logarithm, Log2) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT) = &std::log2;
        test_array<vec_type>("log2",
            [call](VT* vals) {
                return func<vec_type>(call, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).log2();
            });
    }

    TYPED_TEST(Logarithm, Log10) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT) = &std::log10;
        test_array<vec_type>("log10",
            [call](VT* vals) {
                return func<vec_type>(call, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).log10();
            });
    }

    TYPED_TEST(Logarithm, Log1p) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT) = &std::log1p;
        test_array<vec_type>("log1p",
            [call](VT* vals) {
                return func<vec_type>(call, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).log1p();
            });
    }

    TYPED_TEST(Exponents, Exp) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT) = &std::exp;
        test_array<vec_type>("exp",
            [call](VT* vals) {
                return func<vec_type>(call, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).exp();
            });
    }

    TYPED_TEST(Exponents, Expm1) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("expm1",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::expm1), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).expm1();
            });
    }

    TYPED_TEST(ErrorFunctions, Erf) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("erf",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::erf), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).erf();
            });
    }

    TYPED_TEST(ErrorFunctions, Erfc) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        test_array<vec_type>("erfc",
            [](VT* vals) {
                return func<vec_type>(RESOLVE_OVERLOAD(std::erfc), vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).erfc();
            });
    }

    TYPED_TEST(ErrorFunctions, Erfinv) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT) = calc_erfinv;
        test_array<vec_type>("erfinv",
            [call](VT* vals) {
                return func<vec_type>(call, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).erfinv();
            });
    }

    TYPED_TEST(Pow, Pow) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT, VT) = &std::pow;
        test_array<vec_type, BINARY>("pow",
            [call](VT* f, VT* s) {
                return func_bi<vec_type>(call, f, s);
            },
            [](VT* f, VT* s) {
                return vec_type::loadu(f).pow(vec_type::loadu(s));
            });
    }

    TYPED_TEST(LGamma, LGamma) {
        using vec_type = TypeParam;
        using VT = ValueType<TypeParam>;
        VT(*call)(VT) = &std::lgamma;
        test_array<vec_type>("lgamma",
            [call](VT* vals) {
                return func<vec_type>(call, vals);
            },
            [](VT* vals) {
                return vec_type::loadu(vals).lgamma();
            });
    }
#endif
#endif
}

#endif