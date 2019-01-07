#include <stdint.h>

//------------------------------------------------------------------------------

#ifndef __TYPEDEF_HPP
#define __TYPEDEF_HPP

//------------------------------------------------------------------------------

#define DIM 2
#define REAL_TYPE double
#define INDEX_TYPE uint32_t

//------------------------------------------------------------------------------

/// Typedef for reals
typedef REAL_TYPE real_t;

/// Typedef for integers
typedef INDEX_TYPE index_t;

//------------------------------------------------------------------------------

/// Template for array/vector types
template <typename _type, uint32_t _dim>
struct array_t {
  // Constructors
  // Empty; initialize array_t to 0
  array_t () { for (uint32_t i=0;i<_dim;++i) x[i] = 0; }

  // Value; initialize array_to to val
  array_t (const _type& v1) { for (uint32_t i=0;i<_dim;++i) x[i] = v1; }

  // Value 2: initialize first value to v1, all others to v2
  array_t (const _type& v1, const _type& v2) { x[0] = v1; for (uint32_t i=1;i<_dim;++i) x[i] = v2; }

  // Copy-constructor: field of values
  array_t (const _type (&cp)[_dim]) { for (uint32_t i=0;i<_dim;++i) x[i] = cp[i]; }
  // Copy-constructor: array_t
  array_t (const array_t<_type, _dim>& cp) { for (uint32_t i=0;i<_dim;++i) x[i] = cp.x[i]; }

  // access operator
  _type& operator [] (uint32_t i) { return x[i]; }
  const _type& operator [] (uint32_t i) const { return x[i]; }

  // store values and size
  _type x[_dim];
  static const uint32_t dim = _dim;
};
/// Typedef for d-dimensional array of reals
typedef array_t<real_t, DIM> multi_real_t;
/// Typedef for d-dimensional array of integer
typedef array_t<index_t, DIM> multi_index_t;

//------------------------------------------------------------------------------

/// Forward declaration of classes used
class Iterator;
class Geometry;
class Grid;
class Parameter;
class Solver;
class Compute;

#endif // __TYPEDEF_HPP
