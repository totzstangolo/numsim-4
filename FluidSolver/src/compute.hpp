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

#include "typedef.hpp"
#include "precice/SolverInterface.hpp"
//------------------------------------------------------------------------------
#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP
//------------------------------------------------------------------------------
class Compute {
public:
  /// Creates a compute instance with given geometry and parameter
  Compute(const Geometry *geom, const Parameter *param);
  /// Deletes all grids
  ~Compute();

  /// Execute one time step of the fluid simulation (with or without debug info)
  // @ param printInfo print information about current solver state (residual
  // etc.)
  double TimeStep(bool printInfo, precice::SolverInterface *interface,
                  int temperatureID, int heatfluxID, int N, int *vertexIDs,
                  double *vertices,double *temperature,double *heatflux,
                  double &precice_dt,const std::string& coric,
                  const std::string& cowic);

  /// Returns the simulated time in total
  const real_t &GetTime() const;

  /// Returns the pointer to U
  const Grid *GetU() const;
  /// Returns the pointer to V
  const Grid *GetV() const;
  /// Returns the pointer to P
  const Grid *GetP() const;
  /// Returns the pointer to RHS
  const Grid *GetT() const;
  /// Returns the pointer to T
  const Grid *GetRHS() const;

  /// Computes and returns the absolute velocity
  const Grid *GetVelocity();
  /// Computes and returns the vorticity
  const Grid *GetVorticity();
  /// Computes and returns the stream line values
  const Grid *GetStream();
  /// Timestep
  double getTimeStep(double dt);
  /// get temperature at coupling interface
  void GetCoupling_T(double *temperature,int N);
  /// set temperature at coupling interface
  void set_coupl_temp(double *heatflux, int N) const;
  /// fill vertices information
  void Vertices(double *vertices, double *temperature, double *heatflux, int N, int dim);

private:
  // testing
  index_t _iter_count;
  real_t _max_dt;

  // current timestep
  real_t _t;

  // donor-cell diffusion condition (p. 27)
  real_t _dtlimit;

  // limit for residual
  real_t _epslimit;

  // velocities
  Grid *_u;
  Grid *_v;

  // pressure and temperature
  Grid *_p;
  Grid *_T;

  // prel. vel
  Grid *_F;
  Grid *_G;

  // right-hand side
  Grid *_rhs;

  // container for interpolating whichever values
  Grid* _tmp;
  Solver *_solver;

  const Geometry *_geom;
  const Parameter *_param;

  ///Checkpoint data structures
  real_t _tmp_dt;
  real_t _tmp_precice_dt;
  double* _tmp_temperature;
  double* _tmp_heatflux;
  Grid* _tmp_p;
  Grid* _tmp_u;
  Grid* _tmp_v;
  Grid* _tmp_T;
  Grid* _tmp_G;
  Grid* _tmp_F;
  Grid* _tmp_rhs;

  /// Compute the new velocites u,v
  void NewVelocities(const real_t &dt);
  /// Compute the temporary velocites F,G
  void MomentumEqu(const real_t &dt);
  /// Compute the modified temporary velocites ~F,~G
  void ModMomentumEqu(const real_t &dt);
  /// Compute the temporary velocites F,G
  void TempEqu(const real_t &dt);
  /// Compute the RHS of the poisson equation
  void RHS(const real_t &dt);
};
//------------------------------------------------------------------------------
#endif // __COMPUTE_HPP
