#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __ITERATOR_HPP
#define __ITERATOR_HPP
//------------------------------------------------------------------------------
/** Iterator base class
 */
class Iterator {
public:
	/// Constructs a new Iterator depending on a geometry
	Iterator (const Geometry* geom);
	/// Constructs a new Iterator on a geometry with a defined starting value
	Iterator (const Geometry* geom, const index_t& value);

	///	Returns the current position value
	virtual const index_t&	Value	() const;
	/// Cast operator to convert Iterators to integers
	virtual operator const	index_t&() const;
	/// Returns the position coordinates
	virtual multi_index_t	Pos		() const;

	/// Sets the iterator to the first element
	virtual void	First	();
	/// Goes to the next element of the iterator, disables it if position is end
	virtual void	Next	();

	/// Checks if the iterator still has a valid value
	virtual bool	Valid () const;

	/// Returns an Iterator that is located left from this one.
        // if we are at the left boundary, the cell sees itself
	virtual Iterator	Left	() const;

	/// Returns an Iterator that is located right from this one
        // If we are at the right boundary, the cell sees itself
	virtual Iterator	Right	() const;

	/// Returns an Iterator that is located above this one
        // If we are at the upper domain boundary, the cell sees itself
	virtual Iterator	Top		() const;

	/// Returns an Iterator that is located below this one
        // If we are at the lower domain boundary, the cell sees itself
	virtual	Iterator	Down	() const;
	/// Goes to the second next element of the iterator, disables it if position is end
	virtual void DoubleNext();

protected:
	const Geometry*	_geom;
	index_t			_value;
	bool			_valid;
};

//------------------------------------------------------------------------------
/** Iterator for interior cells
 */
class InteriorIterator : public Iterator {
public:
	/// Construct a new InteriorIterator
        InteriorIterator (const Geometry* geom);

	/// Sets the iterator to the first element
	void First ();
	/// Goes to the next element of the iterator, disables it if position is end
	void Next  ();
	/// Goes to the second next element of the iterator, disables it if position is end
	void DoubleNext();
};

//------------------------------------------------------------------------------
/** Iterator for domain boundary cells.
 */
class BoundaryIterator : public Iterator {
public:
	/// Constructs a new BoundaryIterator
	BoundaryIterator (const Geometry* geom);

	/// Sets the boundary to iterate
	void SetBoundary (const index_t& boundary);

	/// Sets the iterator to the first element
	void First ();
	/// Goes to the next element of the iterator, disables it if position is end
	void Next  ();
private:
    index_t		_boundary;
};
//------------------------------------------------------------------------------
#endif // __ITERATOR_HPP
