#include "iterator.hpp"
#include "geometry.hpp"
//------------------------------------------------------------------------------
Iterator::Iterator (const Geometry* geom) : _geom(geom) {
	_valid = true;
	_value = 0;
}
//------------------------------------------------------------------------------
Iterator::Iterator (const Geometry* geom, const index_t& value) : _geom(geom) {
	_value = value;
	_valid = true;
	if (_value >= _geom->Size()[0]*_geom->Size()[1])
		_valid = false;
}
//------------------------------------------------------------------------------
const index_t& Iterator::Value () const {
	return _value;
}

void Iterator::DoubleNext(){
	_value+= 2;
	_valid = ((_geom->Size()[0])*(_geom->Size()[1]) > _value); //Check if the value is even valid
}

//------------------------------------------------------------------------------
Iterator::operator const index_t&() const {
	return _value;
}
//------------------------------------------------------------------------------
multi_index_t Iterator::Pos () const {
	multi_index_t pos;
	pos[0] = _value%_geom->Size()[0];
	pos[1] = _value/_geom->Size()[0];
	return pos;
}
//------------------------------------------------------------------------------
void Iterator::First () {
	_value = 0;
	_valid = true;
}
//------------------------------------------------------------------------------
void Iterator::Next () {
	if (!_valid) return;
	++_value;
	if (_value >= _geom->Size()[0]*_geom->Size()[1])
		_valid = false;
}
//------------------------------------------------------------------------------
bool Iterator::Valid () const {
	return _valid;
}
//------------------------------------------------------------------------------
Iterator Iterator::Left () const {
    if (_value%_geom->Size()[0] == 0) return *this;
    return Iterator(_geom,_value - 1);
}
//------------------------------------------------------------------------------
Iterator Iterator::Right () const {
	if (_value%_geom->Size()[0] == (_geom->Size()[0] - 1)) return *this;
    return Iterator(_geom,_value + 1);
}
//------------------------------------------------------------------------------
Iterator Iterator::Top () const {
	if (_value/_geom->Size()[0] == (_geom->Size()[1] - 1)) return *this;
    return Iterator(_geom,_value + _geom->Size()[0]);
}
//------------------------------------------------------------------------------
Iterator Iterator::Down () const {
	if (_value/_geom->Size()[0] == 0) return *this;
    return Iterator(_geom,_value - _geom->Size()[0]);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
InteriorIterator::InteriorIterator (const Geometry* geom) : Iterator(geom) {
}
//------------------------------------------------------------------------------
void InteriorIterator::First () {
	_value = _geom->Size()[0] + 1;
	_valid = true;
}
//------------------------------------------------------------------------------
void InteriorIterator::Next  () {
    if (!_valid) return;
    ++_value;
    while(_geom->get_cellType(_value) != CellType_t::typeFluid){
    ++_value;
    if (_value%_geom->Size()[0] >= (_geom->Size()[0] - 1)) _value += 2;
    if (_value >= _geom->Size()[0]*(_geom->Size()[1]-1)) _valid = false;
    }
}

void InteriorIterator::DoubleNext(){
	if (this->Pos()[0] >= (_geom->Size()[0] - 3)){
		_value+= 4;
	} else {
		_value+= 2;
	}
    if (_value >= _geom->Size()[0]*(_geom->Size()[1]-1)) _valid = false;
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BoundaryIterator::BoundaryIterator (const Geometry* geom) : Iterator(geom) {
	_boundary = 0;
}
//------------------------------------------------------------------------------
// If coup is true, it is an iterator specifically for coupling boundaries,
// which knows where the boundaries are located
void BoundaryIterator::SetBoundary (const index_t& boundary, bool coup) {
	_boundary = boundary%4;
	_valid = false;
	_coup = coup;
}
//------------------------------------------------------------------------------
void BoundaryIterator::First () {
	_valid = true;
	switch (_boundary) {
	case 0:
		if(_coup){
			_value = _geom->Origin()[0]+_geom->Origin()[1]*_geom->Size()[0];
		}else{
			_value = 0;
		}
		break;
	case 1:
		if(_coup){
			_value = _geom->Origin()[0]-1+(_geom->Origin()[1]+1)*_geom->Size()[0];
		}else{
			_value = 0;
		}
		break;
	case 2:
		if(_coup){
			_value = _geom->Origin()[0]+(_geom->Origin()[1]+1)*_geom->Size()[0]+_geom->Coup()*_geom->Size()[0];
		}else{
			_value = _geom->Size()[0]*(_geom->Size()[1] - 1);
		}
		break;
	case 3:
		if(_coup){
			_value = _geom->Origin()[0]+(_geom->Origin()[1]+1)*_geom->Size()[0]+_geom->Coup();
		}else{
			_value = _geom->Size()[0] - 1;
		}
		break;
	default:
		_boundary = 0;
		_value = 0;
		break;
	};
}

//------------------------------------------------------------------------------
void BoundaryIterator::Next () {
	if (!_valid) return;
	switch (_boundary) {
	case 0:
		++_value;
		if(_coup){
			if (_value >= _geom->Origin()[0]+_geom->Origin()[1]*_geom->Size()[0]+_geom->Coup()) _valid = false;
		}else{
			if (_value >= _geom->Size()[0]) _valid = false;
		}
		break;
	case 1:
		_value += _geom->Size()[0];
		if(_coup){
			if (_value > _geom->Origin()[0]+_geom->Origin()[1]*_geom->Size()[0]+_geom->Coup()*_geom->Size()[0]) _valid = false;
		}else{
			if (_value >= _geom->Size()[0]*_geom->Size()[1]) _valid = false;
		}
		break;
	case 2:
		++_value;
		if(_coup){
			if (_value > _geom->Origin()[0]+(_geom->Origin()[1]+1)*_geom->Size()[0]+_geom->Coup()*_geom->Size()[0]+_geom->Coup()) _valid = false;
		}else{
			if (_value >= _geom->Size()[0]*_geom->Size()[1]) _valid = false;
		}
		break;
	case 3:
		_value += _geom->Size()[0];
		if(_coup){
			if (_value > _geom->Origin()[0]+_geom->Origin()[1]*_geom->Size()[0] + _geom->Coup()*_geom->Size()[0]+ _geom->Coup()) _valid = false;
		}else{
			if (_value >= _geom->Size()[0]*_geom->Size()[1]) _valid = false;
		}
		break;
	default:
		++_value;
		if (_value >= _geom->Size()[0]) _valid = false;
		break;
	};
}
//------------------------------------------------------------------------------
