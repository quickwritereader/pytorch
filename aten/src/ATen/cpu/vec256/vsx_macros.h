#pragma once

#define DEFINE_CMP_OP(op, op_type, func)\
Vec256<op_type> inline __inline_attrs op(const Vec256<op_type> & other) const {\
  return Vec256<op_type> {\
    func(_vec0, other._vec0),\
    func(_vec1, other._vec1)\
  };\
}\
friend Vec256<op_type> __inline_attrs op(const Vec256<op_type>& a,const Vec256<op_type>::vec_internal_type& other_vec)   {\
  return Vec256<op_type> {\
    func(a._vec0, other_vec),\
    func(a._vec1, other_vec)\
  };\
} \
friend Vec256<op_type> __inline_attrs op(const Vec256<op_type>::vec_internal_type& a_vec, const Vec256<op_type>& other)   {\
  return Vec256<op_type> {\
    func(a_vec, other._vec0),\
    func(a_vec, other._vec1)\
  };\
}

#define DEFINE_FRIEND_BITWISE_OP(op, op_type, func) \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type> & b) { \
  return Vec256<op_type> { \
    func(a._vec0, b._vec0), \
    func(a._vec1, b._vec1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type>::vec_internal_type& vec) { \
  return Vec256<op_type> { \
    func(a._vec0, vec), \
    func(a._vec1, vec) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>::vec_internal_type& a_vec, \
  const Vec256<op_type> & b) { \
  return Vec256<op_type> { \
    func(a_vec, b._vec0), \
    func(a_vec, b._vec1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type>::vec_internal_mask_type&  vec) { \
  return Vec256<op_type> { \
    func(a._vecb0, vec), \
    func(a._vecb1, vec) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>::vec_internal_mask_type&  a_vec, \
  const Vec256<op_type> & b) { \
  return Vec256<op_type> { \
    func(a_vec, b._vecb0), \
    func(a_vec, b._vecb1) \
  }; \
}
#define DEFINE_FRIEND_BINARY_OP(op, op_type, func) \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type> & b) { \
  return Vec256<op_type> { \
    func(a._vec0, b._vec0), \
    func(a._vec1, b._vec1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type>::vec_internal_type& vec) { \
  return Vec256<op_type> { \
    func(a._vec0, vec), \
    func(a._vec1, vec) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op(const Vec256<op_type>::vec_internal_type& a_vec, \
  const Vec256<op_type> & b) { \
  return Vec256<op_type> { \
    func(a_vec, b._vec0), \
    func(a_vec, b._vec1) \
  }; \
}
#define DEFINE_FRIEND_TERNARY_OP_MASK(op, op_type, func) \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type> & b, \
      const Vec256<op_type> & mask) { \
  return Vec256<op_type> { \
    func(a._vec0, b._vec0, mask._vecb0), \
    func(a._vec1, b._vec1, mask._vecb1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type> ::vec_internal_type b_vec, \
      const Vec256<op_type> & mask) { \
  return Vec256<op_type> { \
    func(a._vec0, b_vec, mask._vecb0), \
    func(a._vec1, b_vec, mask._vecb1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> ::vec_internal_type a_vec, \
    const Vec256<op_type> & b, \
      const Vec256<op_type> & mask) { \
  return Vec256<op_type> { \
    func(a_vec, b._vec0, mask._vecb0), \
    func(a_vec, b._vec1, mask._vecb1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type> & b, \
      const Vec256<op_type>::vec_internal_mask_type&  mask) { \
  return Vec256<op_type> { \
    func(a._vec0, b._vec0, mask), \
    func(a._vec1, b._vec1, mask) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type> ::vec_internal_type b_vec, \
      const Vec256<op_type>::vec_internal_mask_type&  mask) { \
  return Vec256<op_type> { \
    func(a._vec0, b_vec, mask), \
    func(a._vec1, b_vec, mask) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> ::vec_internal_type a_vec, \
    const Vec256<op_type> & b, \
      const Vec256<op_type>::vec_internal_mask_type&  mask) { \
  return Vec256<op_type> { \
    func(a_vec, b._vec0, mask), \
    func(a_vec, b._vec1, mask) \
  }; \
}
#define DEFINE_FRIEND_TERNARY_OP(op, op_type, func) \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type> & b, \
      const Vec256<op_type> & c) { \
  return Vec256<op_type> { \
    func(a._vec0, b._vec0, c._vec0), \
    func(a._vec1, b._vec1, c._vec1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type>::vec_internal_type& b_vec, \
      const Vec256<op_type> & c) { \
  return Vec256<op_type> { \
    func(a._vec0, b_vec, c._vec0), \
    func(a._vec1, b_vec, c._vec1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type>::vec_internal_type& a_vec, \
    const Vec256<op_type> & b, \
      const Vec256<op_type> & c) { \
  return Vec256<op_type> { \
    func(a_vec, b._vec0, c._vec0), \
    func(a_vec, b._vec1, c._vec1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type>& a , \
    const Vec256<op_type> & b, \
      const Vec256<op_type>::vec_internal_type & c_vec) { \
  return Vec256<op_type> { \
    func(a._vec0, b._vec0, c_vec), \
    func(a._vec1, b._vec1, c_vec) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type>::vec_internal_type& a_vec, \
    const Vec256<op_type>::vec_internal_type& b_vec, \
      const Vec256<op_type> & c) { \
  return Vec256<op_type> { \
    func(a_vec, b_vec, c._vec0), \
    func(a_vec, b_vec, c._vec1) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type> & a, \
    const Vec256<op_type>::vec_internal_type& b_vec, \
      const Vec256<op_type>::vec_internal_type& c_vec) { \
  return Vec256<op_type> { \
    func(a._vec0, b_vec, c_vec), \
    func(a._vec1, b_vec, c_vec) \
  }; \
} \
friend Vec256<op_type> inline __inline_attrs op( \
  const Vec256<op_type>::vec_internal_type& a_vec, \
    const Vec256<op_type> & b, \
      const Vec256<op_type>::vec_internal_type& c_vec) { \
  return Vec256<op_type> { \
    func(a_vec, b._vec0, c_vec), \
    func(a_vec, b._vec1, c_vec) \
  }; \
}