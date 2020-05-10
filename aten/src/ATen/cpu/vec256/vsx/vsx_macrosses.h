#pragma once
#define DEFINE_MEMBER_OP(op, op_type, func)                              \
  Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>& other) \
      const {                                                            \
    return Vec256<op_type>{                                              \
        func(_vec0, other._vec0), func(_vec1, other._vec1)};             \
  }

#define DEFINE_MEMBER_BITWISE_OP(op, op_type, func)                      \
  Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>& other) \
      const {                                                            \
    return Vec256<op_type>{                                              \
        func(_vecb0, other._vecb0), func(_vecb1, other._vecb1)};         \
  } 


#define DEFINE_MEMBER_TERNARY_OP(op, op_type, func)                        \
   Vec256<op_type>  __inline_attrs op(                                     \
      const Vec256<op_type>& b,                                            \
      const Vec256<op_type>& c) {                                          \
    return Vec256<op_type>{                                                \
        func(_vec0, b._vec0, c._vec0), func(_vec1, b._vec1, c._vec1)}; \
  }                                                                         

#define DEFINE_MEMBER_EMULATE_BINARY_OP(op, op_type, binary_op)  \
   Vec256<op_type>  __inline_attrs op(          \
       const Vec256<op_type>& b) { \
    Vec256<op_type>::vec_internal_type ret_0;               \
    Vec256<op_type>::vec_internal_type ret_1;               \
    for (size_t i = 0; i < Vec256<op_type>::size() / 2; i++) { \
      ret_0[i] = _vec0[i] binary_op b._vec0[i];           \
      ret_1[i] = _vec1[i] binary_op b._vec1[i];           \
    }                                                       \
    return Vec256<op_type>{ret_0, ret_1};                   \
  }                                                         

#define DEFINE_MEMBER_OP_AND_ONE(op, op_type, func)                          \
  Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>& other)     \
      const {                                                                \
    using vvtype =Vec256<op_type>::vec_internal_type;                          \
    const vvtype v_one = vec_splats(static_cast<op_type>(1.0)); \
    vvtype ret0 = (vvtype)func(_vec0, other._vec0);      \
    vvtype ret1 = (vvtype)func(_vec1, other._vec1);      \
    return Vec256<op_type>{vec_and(ret0, v_one), vec_and(ret1, v_one)};      \
  }






