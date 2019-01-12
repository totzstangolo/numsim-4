
#include "typedef.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"
#include "compute.hpp"
#include "iterator.hpp"
#include "solver.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <algorithm>

using namespace precice;
using namespace precice::constants;
using namespace std;

Compute::Compute(const Geometry *geom, const Parameter *param){
	_geom = geom;
	_param = param;
	//debug init
	_max_dt = 1.0;

	// Initializing time, dtlimit and epslimit_
	_t = 0.0;
	multi_real_t h = _geom->Mesh();

	multi_real_t h2;
	h2[0] = h[0]*h[0];
	h2[1] = h[1]*h[1];
	_dtlimit = h2[0]*h2[1]*_param->Re()/(2*(h2[0]+h2[1]));

	if (_param->Pr()>0)
		_dtlimit = std::min<real_t>(_param->Re()*_param->Pr()/(2*(1/h2[0]+1/h2[1])),_dtlimit);

	_epslimit = _param->Eps() * _param->Eps() * _geom->Size()[0] * _geom->Size()[1];

	// creating grids with offset
	multi_real_t compute_offset_x;
	compute_offset_x[0] = 0.0;
	compute_offset_x[1] = -0.5 * h[1];
	_u = new Grid(_geom, compute_offset_x);
	_F = new Grid(_geom, compute_offset_x);
	_u->Initialize(0);
	_geom->Update_U(_u);

	multi_real_t compute_offset_y;
	compute_offset_y[0] = -0.5 * h[0];
	compute_offset_y[1] = 0.0;
	_v = new Grid(_geom, compute_offset_y);
	_G = new Grid(_geom, compute_offset_y);
	_v->Initialize(0);
	_geom->Update_V(_v);

	multi_real_t compute_offset_p;
	compute_offset_p[0] = -0.5 * h[0];
	compute_offset_p[1] = -0.5 * h[1];
	_p = new Grid(_geom, compute_offset_p);
	_T = new Grid(_geom, compute_offset_p);
	_rhs = new Grid(_geom, compute_offset_p);
	_tmp = new Grid(_geom, compute_offset_p);
	_p->Initialize(0);
	_T->Initialize(0);

	//create solver (used script omega, not param omega)
	_solver = new SOR(_geom, _param->Omega());


}

Compute::~Compute(){
	delete _u;
	delete _F;
	delete _v;
	delete _G;
	delete _p;
	delete _T;
	delete _rhs;
	delete _tmp;
	delete _solver;
}

void Compute::set_coupl_temp(double *heatflux, int N) const{
	_geom->UpdateCoupling_T(_T,heatflux,N,_param->Expl());
}

void Compute::GetCoupling_T(double* temperature, int N){
	_geom->GetCoupling_T(_T,temperature,N,_param->Expl());
}

double Compute::getTimeStep(double dt){
	if(dt == 0){
		dt = abs(std::min<real_t>(_geom->Mesh()[0]/_u->AbsMax(),_geom->Mesh()[1]/_v->AbsMax()));
		dt = std::min<real_t>(dt,_dtlimit);
		dt *= _param->Tau();
	}

	return std::min<real_t>(_max_dt, dt);

}

void Compute::Vertices(double *vertices, double *temperature, double *heatflux, int N, int dim){
	if(_param->Expl()){
		for (int i = 0; i < N; i++) {
	          temperature[i] = 0.0;
	          heatflux[i]    = 1.0;
	          for (int j = 0; j < dim; j++){
	              if(j==0) {
	                  vertices[i * dim + j] = (i * (1 - j) + _geom->Origin()[0] + 0.5)*_geom->Mesh()[0];
	              }else{
	                  if(j==1){
	                      vertices[i * dim + j] = 0.25*_geom->Mesh()[1];
	                  }else{
	                      vertices[i * dim + j] = 0;
	                  }
	              }
	          }
	    }
	}else{
		for(int i=0;i<N;i++){
			temperature[i] = 303.0;
			heatflux[i]    = 0.0;
		}
		double *temp=new double[N*dim];
		int count = 0;
		BoundaryIterator it(_geom);
		for(int k = 1; k<=4; k++){
			it.SetBoundary(k-1,true);
			it.First();
			for (int i = 0; i < _geom->Coup() && it.Valid(); i++) {
		          for (int j = 0; j < dim; j++){
		              if(j==0) {
						if(k-1==0 || k-1==2) {
							// x-coordinate of upper or lower bound
							vertices[k*(i*dim+j)] = (it.Pos()[0]+0.5)*_geom->Mesh()[0];
						}else{
							// x-coordinate of sides
						  	vertices[k*(i*dim+j)] = it.Pos()[0]*_geom->Mesh()[0];
					    }
		              }else{
		                  if(j==1){
							  if(k-1==1 || k-1==3) {
								// y-coordinate of sides
	  							vertices[k*(i*dim+j)] = (it.Pos()[1]+0.5)*_geom->Mesh()[1];
	  						}else{
								if(k-1==0){
									// y-coordinate of lower bound
									vertices[k*(i*dim+j)] = _geom->Mesh()[1];
								}else{
									// y-coordinate of upper bound
									vertices[k*(i*dim+j)] = it.Pos()[1]*_geom->Mesh()[1];
								}
	  					    }
		                  }else{
		                      vertices[k*(i*dim+j)] = 0;
		                  }
		              }
					  temp[count] = vertices[k*(i*dim+j)];
					  count++;
		          }
				  it.Next();
		          // countThat++;
		    }
		}
		for(int i=0;i<3*N;i++){
			vertices[i] = temp[i];
			if((i+1)%3!=0) vertices[i]+=0.15;
		}
	}
}

double Compute::TimeStep(bool printInfo, SolverInterface *interface,int temperatureID, int heatfluxID, int N, int *vertexIDs,
double *vertices,double *temperature,double *heatflux,double &precice_dt){
	(*interface).readBlockScalarData(heatfluxID,N,vertexIDs,heatflux); // read new heatflux from preCICE buffers
	set_coupl_temp(heatflux,N);
	// Compute like in script page 23
	//compute dt
	real_t dt = _param->Dt();
	dt = abs(std::min<real_t>(_geom->Mesh()[0]/_u->AbsMax(),_geom->Mesh()[1]/_v->AbsMax()));
	dt = std::min<real_t>(dt,_dtlimit);
	dt *= _param->Tau();

	if(printInfo) {
		printf("dt: %f \n", dt);
	}

	_max_dt = std::min<real_t>(_max_dt, dt);

	/*
	///////////////////////////
	**/

	dt = std::min(precice_dt,_max_dt);
	/**
	///////////////////////////
	*/

	// compute temperature and update bound
	TempEqu(dt);
	_geom->Update_T(_T, _param->T_h(), _param->T_c(), false);

	// compute FG and update bound.
	MomentumEqu(dt);
	ModMomentumEqu(dt);
	_geom->Update_U(_F);
	_geom->Update_V(_G);

	// compute rhs and update bound.
	RHS(dt);
	_geom->Update_P(_rhs);

	// Solver, relative eps had bad performance
	index_t index = 0;
	real_t res = 1;
	while(index < _param->IterMax() && res > _epslimit) {
		index++;
		res = _solver->Cycle(_p, _rhs);
		_geom->Update_P(_p);
	}
	if(printInfo) {
		printf("iterations: %d / %d \n", index, _param->IterMax());
	printf("pmax: %f\n", _p->Max());
	}
	_iter_count += index;

	//compute uv and update bound.
	NewVelocities(dt);
	_geom->Update_U(_u);
	_geom->Update_V(_v);

	_t += dt;

	if(printInfo)
		printf("time: %f\n", _t);

	GetCoupling_T(temperature,N);
	(*interface).writeBlockScalarData(temperatureID,N,vertexIDs,temperature); // write new temperature to preCICE buffers
    precice_dt = (*interface).advance(dt); // advance coupling
	return dt;
}

  /// Returns the simulated time in total
const real_t& Compute::GetTime() const{
    return _t;
}

  /// Returns the pointer to U
const Grid* Compute::GetU() const{
    return _u;
}
  /// Returns the pointer to V
const Grid* Compute::GetV() const{
    return _v;
}
  /// Returns the pointer to P
const Grid* Compute::GetP() const{
    return _p;
}

/// Returns the pointer to T
const Grid* Compute::GetT() const{
  return _T;
}

  /// Returns the pointer to RHS
const Grid* Compute::GetRHS() const{
    return _rhs;
}

/// Computes and returns the absolute velocity
const Grid* Compute::GetVelocity(){
	multi_real_t mid;
	Iterator iterator(_geom);
	for (iterator.First(); iterator.Valid(); iterator.Next()){
		mid[0] = _u->Cell(iterator) + _u->Cell(iterator.Left());
		mid[1] = _v->Cell(iterator) + _v->Cell(iterator.Down());
		_tmp->Cell(iterator) = sqrt((mid[0] * mid[0]) + (mid[1] * mid[1]));
	}
	return _tmp;
}

/// Computes and returns the vorticity
const Grid* Compute::GetVorticity(){
	// TODO real implementation needed
    Grid *grid = new Grid(_geom);
	return grid;
}

/// Computes and returns the stream line values
const Grid* Compute::GetStream(){
	// TODO real implementation needed
    Grid *grid = new Grid(_geom);
    return grid;
}

void Compute::NewVelocities(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_u->Cell(intIterator) = _F->Cell(intIterator) - (dt * _p->dx_r(intIterator));
		_v->Cell(intIterator) = _G->Cell(intIterator) - (dt * _p->dy_r(intIterator));
	}
}

/// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_F->Cell(intIterator) = _u->Cell(intIterator) +
			(
				 (_u->dxx(intIterator) + _u->dyy(intIterator))/ _param->Re()
				- _u->DC_udu_x(intIterator, _param->Alpha())
				- _u->DC_vdu_y(intIterator, _param->Alpha(), _v
			)) * dt;
		_G->Cell(intIterator) = _v->Cell(intIterator) +
			(
			 (_v->dyy(intIterator) + _v->dxx(intIterator))/ _param->Re()
				- _v->DC_vdv_y(intIterator, _param->Alpha())
				- _v->DC_udv_x(intIterator, _param->Alpha(), _u)
				) * dt;
	}

}

/// Compute the temporary velocites F,G
void Compute::ModMomentumEqu(const real_t &dt){
	InteriorIterator intIterator(_geom);
	real_t tij = 0.0;
	real_t ti1j = 0.0;
	real_t tij1 = 0.0;
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		tij = _T->Cell(intIterator);
		ti1j = _T->Cell(intIterator.Right());
		tij1 = _T->Cell(intIterator.Top());
		_F->Cell(intIterator) = _F->Cell(intIterator) -
			dt*_param->Beta()*_param->Gx()*0.5*(tij+ti1j);
		_G->Cell(intIterator) = _G->Cell(intIterator) -
			dt*_param->Beta()*_param->Gy()*0.5*(tij+tij1);
	}
}

void Compute::TempEqu(const real_t &dt){
	InteriorIterator intIterator(_geom);
	real_t txx = 0.0;
	real_t tyy = 0.0;
	real_t utx = 0.0;
	real_t vty = 0.0;
	real_t pref = 1/(_param->Re()*_param->Pr());
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		txx = _T->dxx(intIterator);
		tyy = _T->dyy(intIterator);
		utx = _T->DC_duT_x(intIterator,_param->Alpha(),_u);
		vty = _T->DC_dvT_y(intIterator,_param->Alpha(),_v);
		_T->Cell(intIterator) = _T->Cell(intIterator) +
				(pref*(txx+tyy)-utx-vty)*dt;
	}
}

/// Compute the RHS of the poisson equation
void Compute::RHS(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_rhs->Cell(intIterator) = (_F->dx_l(intIterator) + _G->dy_l(intIterator)) / dt;
	}
}
