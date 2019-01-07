/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "communicator.hpp"
#include "parameter.hpp"
#include <stdio.h>
#include <string.h>


/*
 *
 * typeInH:
 * HHHHHHHHHHHHHHHHHHHHHHHHHHHH
 *  \     |     |     |      /
 *   \    v     v     v     /
 *     ---------------------
 *
 * typeInV:
 *  V \
 *  V >\
 *  V > |
 *  V > |
 *  V >/
 *  V /
 *
 * typeSlipH:
 * ----------------------
 *  ->  ->  ->  ->  ->
 *
 * typeSlipV:
 * |
 * | ^
 * | |
 * |
 * | ^
 * | |
 * |
 *
 * */

//------------------------------------------------------------------------------
void Geometry::UpdateCellDirichlet_U(Grid *u, const real_t &value,
                                     const Iterator &it) const {
  switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].fluid) {
  case cellW:
    u->Cell(it.Left()) = value;
    u->Cell(it) = value;
    break;
  case cellNW:
    u->Cell(it.Left()) = value;
  case cellN:
    u->Cell(it) = 2.0 * value - u->Cell(it.Top());
    break;
  case cellSW:
    u->Cell(it.Left()) = value;
  case cellS:
    u->Cell(it) = 2.0 * value - u->Cell(it.Down());
    break;
  default:
    u->Cell(it) = value;
    break;
  };
}
//------------------------------------------------------------------------------
void Geometry::UpdateCellDirichlet_V(Grid *v, const real_t &value,
                                     const Iterator &it) const {
    switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].fluid) {
    case cellS:
        v->Cell(it.Down()) = value;
        v->Cell(it) = value;
        break;
    case cellSE:
        v->Cell(it.Down()) = value;
    case cellE:
        v->Cell(it) = 2.0 * value - v->Cell(it.Right());
        break;
    case cellSW:
        v->Cell(it.Down()) = value;
    case cellW:
        v->Cell(it) = 2.0 * value - v->Cell(it.Left());
        break;
    default:
    	v->Cell(it) = value;
    	break;
    };
}
//------------------------------------------------------------------------------
void Geometry::UpdateCellDirichlet_T(Grid *grid, const real_t &value, const Iterator &it) const {
    switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].fluid) {
    case cellW:
        grid->Cell(it) = 2.0 * value - grid->Cell(it.Left());
        break;
    case cellN:
        grid->Cell(it) = 2.0 * value - grid->Cell(it.Top());
        break;
    case cellE:
        grid->Cell(it) = 2.0 * value - grid->Cell(it.Right());
        break;
    case cellS:
        grid->Cell(it) = 2.0 * value - grid->Cell(it.Down());
        break;
    default:
        break;
  };
}//------------------------------------------------------------------------------
void Geometry::UpdateCellNeumann(Grid *grid, const Iterator &it) const {
    switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].fluid) {
    case cellW:
        grid->Cell(it) = grid->Cell(it.Left());
        break;
    case cellNW:
        grid->Cell(it) = 0.5 * (grid->Cell(it.Left()) + grid->Cell(it.Top()));
        break;
    case cellN:
        grid->Cell(it) = grid->Cell(it.Top());
        break;
    case cellNE:
        grid->Cell(it) = 0.5 * (grid->Cell(it.Right()) + grid->Cell(it.Top()));
        break;
    case cellE:
        grid->Cell(it) = grid->Cell(it.Right());
        break;
    case cellSE:
        grid->Cell(it) = 0.5 * (grid->Cell(it.Right()) + grid->Cell(it.Down()));
        break;
    case cellS:
        grid->Cell(it) = grid->Cell(it.Down());
        break;
    case cellSW:
        grid->Cell(it) = 0.5 * (grid->Cell(it.Left()) + grid->Cell(it.Down()));
        break;
    default:
        grid->Cell(it) = grid->Cell(it.Left());
        break;
  };
}
//------------------------------------------------------------------------------
void Geometry::UpdateCellNeumann_P(Grid *grid, const Iterator &it) const {
  UpdateCellNeumann(grid, it);
}
//------------------------------------------------------------------------------
Geometry::Geometry() : _comm(NULL) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;

  // create boundary halo
  _size[0] += 2;
  _size[1] += 2;

  _bsize = _size;
  _blength = _length;

  _cell = NULL;
  _boffset = 0;
}
//------------------------------------------------------------------------------
Geometry::Geometry(const Communicator *comm) : _comm(comm) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;

  _bsize[0] = _size[0] / _comm->ThreadDim()[0] + 2;
  _bsize[1] = _size[1] / _comm->ThreadDim()[1] + 2;
  if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
    _bsize[0] += _size[0] % _comm->ThreadDim()[0];
  if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
    _bsize[1] += _size[1] % _comm->ThreadDim()[0];

  _blength[0] = _h[0] * (_bsize[0] - 2);
  _blength[1] = _h[1] * (_bsize[1] - 2);

  // create boundary halo
  _size[0] += 2;
  _size[1] += 2;

  _cell = NULL;
  _boffset = 0;
}
//------------------------------------------------------------------------------
Geometry::~Geometry() {
  if (_cell)
    delete[] _cell;
}
//------------------------------------------------------------------------------
void Geometry::Load(const char *file) {
  FILE *handle = fopen(file, "r");
  double inval[2];
  char name[200000];
  while (!feof(handle)) {
    if (!fscanf(handle, "%s =", name))
      continue;
    if (strcmp(name, "size") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _size[0] = inval[0];
        _size[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "length") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _length[0] = inval[0];
        _length[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "velocity") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _velocity[0] = inval[0];
        _velocity[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "pressure") == 0) {
      if (fscanf(handle, " %lf\n", &inval[0]))
        _pressure = inval[0];
      continue;
    }
    if (strcmp(name, "geometry") == 0) {
      if (!fscanf(handle, "%s\n", name))
        continue;
      if (strcmp(name, "free") == 0) {
        if (_cell)
          delete[] _cell;
        _cell = new Cell_t[_size[0] * _size[1]];
        bool parabolic = false;
        // Read stuff from file
        for (int y = _size[1]; y-- > 0;) {
          if (feof(handle)) {
            delete[] _cell;
            _cell = NULL;
            break;
          }
          for (int x = 0; x < _size[0]; ++x) {
            _cell[x + y * _size[0]].fluid = cellNone;
            _cell[x + y * _size[0]].factor = 1.0;
            switch (getc(handle)) {
            case '#':
              _cell[x + y * _size[0]].type = typeSolid;
              break;
            case 'I':
              _cell[x + y * _size[0]].type = typeIn;
              break;
            case 'O':
              _cell[x + y * _size[0]].type = typeOut;
              break;
            case '|':
              _cell[x + y * _size[0]].type = typeSlipV;
              break;
            case 'T':
              _cell[x + y * _size[0]].type = typeHot;
              break;
            case 'C':
              _cell[x + y * _size[0]].type = typeCold;
              break;
            case 'W':
              _cell[x + y * _size[0]].type = typeInsul;
              break;
            case '-':
              _cell[x + y * _size[0]].type = typeSlipH;
              break;
            case 'H':
              _cell[x + y * _size[0]].type = typeInH;
              parabolic = true;
              break;
            case 'V':
              _cell[x + y * _size[0]].type = typeInV;
              parabolic = true;
              break;
            default:
              if (x == 0 || x == _size[0] - 1 || y == 0 || y == _size[1] - 1)
                _cell[x + y * _size[0]].type = typeSolid;
              else
                _cell[x + y * _size[0]].type = typeFluid;
              break;
            };
          }
          if (!fscanf(handle, "\n"))
            continue;
        }
        if (!_cell)
          break;
        // Process it
        for (int y = 0; y < _size[1]; ++y) {
          for (int x = 0; x < _size[0]; ++x) {
            int check = 0;
            if (_cell[x + y * _size[0]].type == typeFluid)
              continue;
            if (x < _size[0] - 1 &&
                _cell[x + 1 + y * _size[0]].type == typeFluid)
              check |= 8;
            if (x > 0 && _cell[x - 1 + y * _size[0]].type == typeFluid)
              check |= 2;
            if (y < _size[1] - 1 &&
                _cell[x + (y + 1) * _size[0]].type == typeFluid)
              check |= 1;
            if (y > 0 && _cell[x + (y - 1) * _size[0]].type == typeFluid)
              check |= 4;
            switch (check) {
            case 5:
            case 7:
            case 10:
            case 11:
            case 13:
            case 14:
            case 15:
              _cell[x + y * _size[0]].type = typeFluid;
              _cell[x + y * _size[0]].fluid = cellNone;
              if (x > 0)
                x--;
              if (y > 0)
                y--;
              break;
            case 1:
              _cell[x + y * _size[0]].fluid = cellN;
              break;
            case 2:
              _cell[x + y * _size[0]].fluid = cellW;
              break;
            case 3:
              _cell[x + y * _size[0]].fluid = cellNW;
              break;
            case 4:
              _cell[x + y * _size[0]].fluid = cellS;
              break;
            case 6:
              _cell[x + y * _size[0]].fluid = cellSW;
              break;
            case 8:
              _cell[x + y * _size[0]].fluid = cellE;
              break;
            case 9:
              _cell[x + y * _size[0]].fluid = cellNE;
              break;
            case 12:
              _cell[x + y * _size[0]].fluid = cellSE;
              break;
            };
          }
        }
        // Parabolic stuff
        if (parabolic) {
          for (int y = 0; y < _size[1]; ++y) {        // uint32_t durch int ersetzt
            for (int x = 0; x < _size[0]; ++x) {
              int32_t dist1 = 0;
              int32_t dist2 = 0;
              switch (_cell[x + y * _size[0]].type) {
              case typeInH:
                while (x - dist1 >= 0 &&
                       _cell[x - dist1 + y * _size[0]].type == typeInH)
                  ++dist1;
                while (x + dist2 < _size[0] &&
                       _cell[x + dist2 + y * _size[0]].type == typeInH)
                  ++dist2;
                _cell[x + y * _size[0]].factor =
                    4.0 * ((real_t)(dist1)-0.5) * ((real_t)(dist2)-0.5) /
                    (real_t)((dist1 + dist2 - 1) * (dist1 + dist2 - 1));
                break;
              case typeInV:
                while (y - dist1 >= 0 &&
                       _cell[x + (y - dist1) * _size[0]].type == typeInV)
                  ++dist1;
                while (y + dist2 < _size[1] &&
                       _cell[x + (dist2 + y) * _size[0]].type == typeInV)
                  ++dist2;
                _cell[x + y * _size[0]].factor =
                    4.0 * ((real_t)(dist1)-0.5) * ((real_t)(dist2)-0.5) /
                    (real_t)((dist1 + dist2 - 1) * (dist1 + dist2 - 1));
                break;
              default:
                break;
              };
            }
          }
        }
        _size[0] -= 2;
        _size[1] -= 2;
      }
    }
  }
  fclose(handle);
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];

  _blength = _length;

  _boffset = 0;

  if (_comm) {
    _bsize[0] = _size[0] / _comm->ThreadDim()[0] + 2;
    _bsize[1] = _size[1] / _comm->ThreadDim()[1] + 2;
    _boffset = (_bsize[0] - 2) * _comm->ThreadIdx()[0] +
               (_size[0] + 2) * (_bsize[1] - 2) * _comm->ThreadIdx()[1];
    if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
      _bsize[0] += _size[0] % _comm->ThreadDim()[0];
    if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
      _bsize[1] += _size[1] % _comm->ThreadDim()[0];

    _blength[0] = _h[0] * (_bsize[0] - 2);
    _blength[1] = _h[1] * (_bsize[1] - 2);
  }

  _size[0] += 2;
  _size[1] += 2;

  if (!_comm)
    _bsize = _size;
}
//------------------------------------------------------------------------------
const multi_index_t &Geometry::Size() const { return _bsize; }
//------------------------------------------------------------------------------
const multi_index_t &Geometry::TotalSize() const { return _size; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::Length() const { return _blength; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::TotalLength() const { return _length; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::Mesh() const { return _h; }
//------------------------------------------------------------------------------
void Geometry::Update_U(Grid *u) const {
  if (_cell) {
    Iterator it(this);
    for (it.First(); it.Valid(); it.Next()) {
      switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].type) {
      case typeHot:
      case typeCold:
      case typeInsul:
      case typeSlipV:
      case typeSolid:
        UpdateCellDirichlet_U(u, 0.0, it);
        break;
      case typeIn:
      case typeInH:
        UpdateCellDirichlet_U(u, _velocity[0], it);
        break;
      case typeSlipH:
      case typeOut:
        // std::cout << _size[0] << std::endl;
        UpdateCellNeumann(u, it);
        break;
      case typeInV:
        UpdateCellDirichlet_U(
            u,
            _velocity[0] *
                _cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].factor,
            it);
        break;
      default:
        break;
      };
    }
  } else {
    BoundaryIterator it(this);
    if (_comm && _comm->isBottom()) {
      it.SetBoundary(0);
      for (it.First(); it.Valid(); it.Next())
        u->Cell(it) = -u->Cell(it.Top());
    }
    if (_comm && _comm->isLeft()) {
      it.SetBoundary(1);
      for (it.First(); it.Valid(); it.Next())
        u->Cell(it) = 0.0;
    }
    if (_comm && _comm->isRight()) {
      it.SetBoundary(3);
      for (it.First(); it.Valid(); it.Next()) {
        u->Cell(it) = 0.0;
        u->Cell(it.Left()) = 0.0;
      }
    }
    if (_comm && _comm->isTop()) {
      it.SetBoundary(2);
      for (it.First(); it.Valid(); it.Next())
        u->Cell(it) = 2.0 * _velocity[0] - u->Cell(it.Down());
    }
  }
}

//------------------------------------------------------------------------------
void Geometry::Update_V(Grid *v) const {
  if (_cell) {
    Iterator it(this);
    for (it.First(); it.Valid(); it.Next()) {
      switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].type) {
      case typeHot:
      case typeCold:
      case typeInsul:
      case typeSlipH:
      case typeSolid:
        UpdateCellDirichlet_V(v, 0, it);
        break;
      case typeIn:
      case typeInV:
        UpdateCellDirichlet_V(v, _velocity[1], it);
        break;
      case typeSlipV:
      case typeOut:
        UpdateCellNeumann(v, it);
        break;
      case typeInH:
        UpdateCellDirichlet_V(
            v,
            _velocity[1] *
                _cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].factor,
            it);
        break;
      default:
        break;
      };
    }
  } else {
    BoundaryIterator it(this);
    if (_comm && _comm->isBottom()) {
      it.SetBoundary(0);
      for (it.First(); it.Valid(); it.Next())
        v->Cell(it) = 0.0;
    }
    if (_comm && _comm->isTop()) {
      it.SetBoundary(2);
      for (it.First(); it.Valid(); it.Next()) {
        v->Cell(it) = 0.0;
        v->Cell(it.Down()) = 0.0;
      }
    }
    if (_comm && _comm->isLeft()) {
      it.SetBoundary(1);
      for (it.First(); it.Valid(); it.Next())
        v->Cell(it) = -v->Cell(it.Right());
    }
    if (_comm && _comm->isRight()) {
      it.SetBoundary(3);
      for (it.First(); it.Valid(); it.Next())
        v->Cell(it) = -v->Cell(it.Left());
    }
  }
}
//------------------------------------------------------------------------------
void Geometry::Update_P(Grid *p) const {
  if (_cell) {
    Iterator it(this);
    for (it.First(); it.Valid(); it.Next()) {
      switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].type) {
      case typeHot:
      case typeCold:
      case typeInsul:
      case typeIn:
      case typeInH:
      case typeInV:
      case typeSolid:
        UpdateCellNeumann_P(p, it);
        break;
      case typeSlipH:
      case typeSlipV:
        p->Cell(it) = _pressure;
        break;
      case typeOut:
        p->Cell(it) = 0.0;
        break;
      default:
        break;
      };
    }
  } else {
	  printf("Hi \n");
    BoundaryIterator it(this);
    if (_comm && _comm->isBottom()) {
      it.SetBoundary(0);
      for (it.First(); it.Valid(); it.Next())
        p->Cell(it) = p->Cell(it.Top());
    }
    if (_comm && _comm->isTop()) {
      it.SetBoundary(2);
      for (it.First(); it.Valid(); it.Next())
        p->Cell(it) = p->Cell(it.Down());
    }
    if (_comm && _comm->isLeft()) {
      it.SetBoundary(1);
      for (it.First(); it.Valid(); it.Next())
        p->Cell(it) = p->Cell(it.Right());
    }
    if (_comm && _comm->isRight()) {
      it.SetBoundary(3);
      for (it.First(); it.Valid(); it.Next())
        p->Cell(it) = p->Cell(it.Left());
    }
  }
}

//------------------------------------------------------------------------------

void Geometry::Update_T(Grid *t, real_t hot, real_t cold) const {
  if (_cell) {
/*    Iterator it(this);
    for (it.First(); it.Valid(); it.Next()) {
      switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].type) {
      case typeInsul:
      case typeIn:
      case typeInH:
      case typeInV:
      case typeSolid:
      case typeSlipH:
      case typeSlipV:
      case typeOut:
    	  UpdateCellNeumann_P(t, it);
        break;
      case typeHot:
    	  UpdateCellDirichlet_T(t, hot, it);
    	break;
      case typeCold:
    	  UpdateCellDirichlet_T(t, cold, it);
      default:
        break;
      };
    }

  } else {*/
    BoundaryIterator it(this);
    if (_comm && _comm->isBottom()) {
      it.SetBoundary(0);
      for (it.First(); it.Valid(); it.Next())
          switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].type) {
          case typeHot:
          t->Cell(it) = 2.0 * hot - t->Cell(it.Top());
          break;
          case typeCold:
          t->Cell(it) = 2.0 * cold - t->Cell(it.Top());
          break;
          default:
          t->Cell(it) = t->Cell(it.Top());
          break;
          }
    }
    if (_comm && _comm->isTop()) {
      it.SetBoundary(2);
      for (it.First(); it.Valid(); it.Next())
          switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].type) {
          case typeHot:
          t->Cell(it) = 2.0 * hot - t->Cell(it.Down());
          break;
          case typeCold:
          t->Cell(it) = 2.0 * cold - t->Cell(it.Down());
          break;
          default:
          t->Cell(it) = t->Cell(it.Down());
          break;
          }
    }
    if (_comm && _comm->isLeft()) {
      it.SetBoundary(1);
      for (it.First(); it.Valid(); it.Next())
          switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].type) {
          case typeHot:
          t->Cell(it) = 2.0 * hot - t->Cell(it.Right());
          break;
          case typeCold:
          t->Cell(it) = 2.0 * cold - t->Cell(it.Right());
          break;
          default:
          t->Cell(it) = t->Cell(it.Right());
          break;
          }
    }
    if (_comm && _comm->isRight()) {
      it.SetBoundary(3);
      for (it.First(); it.Valid(); it.Next())
          switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1] * _size[0]].type) {
          case typeHot:
          t->Cell(it) = 2.0 * hot - t->Cell(it.Left());
          break;
          case typeCold:
          t->Cell(it) = 2.0 * cold - t->Cell(it.Left());
          break;
          default:
          t->Cell(it) = t->Cell(it.Left());
          break;
          }
    }
  }
}

//------------------------------------------------------------------------------
