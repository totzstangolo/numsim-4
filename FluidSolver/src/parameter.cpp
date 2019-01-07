#include "parameter.hpp"
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>

Parameter::Parameter(){
// Load(file);

}
// Loads parameters from a file

void Parameter::Load(const char *file){
    std::cout << "Start Load" << std::endl;
FILE* handle = fopen(file,"r");
double inval;
char name[20000];
while (!feof(handle)) {
	if (!fscanf(handle, "%s = %lf\n", name, &inval)) continue;
	if (strcmp(name,"re") == 0) _re = inval;
	else if (strcmp(name,"omg") == 0) _omega = inval;
	else if (strcmp(name,"alpha") == 0) _alpha = inval;
    else if (strcmp(name,"beta") == 0) _beta = inval;
	else if (strcmp(name,"dt") == 0) _dt = inval;
	else if (strcmp(name,"tend") == 0) _tend = inval;
	else if (strcmp(name,"iter") == 0) _itermax = inval;
	else if (strcmp(name,"eps") == 0) _eps = inval;
	else if (strcmp(name,"tau") == 0) _tau = inval;
	else if (strcmp(name,"imax") == 0) _imax = inval;
	else if (strcmp(name,"jmax") == 0) _jmax = inval;
	else if (strcmp(name,"dt_value") == 0) _dt_value = inval;
	else if (strcmp(name,"gx") == 0) _gx = inval;
	else if (strcmp(name,"gy") == 0) _gy = inval;
    else if (strcmp(name,"pr") == 0) _pr = inval;
	else if (strcmp(name,"T_h") == 0) _T_h = inval;
	else if (strcmp(name,"T_c") == 0) _T_c = inval;
	else printf("Unknown parameter %s\n",name);
}
fclose(handle);
 }

 // Get  methods

const real_t &Parameter::Re()const{
return _re;
}

const real_t &Parameter::Omega()const{

return _omega;
}

const real_t &Parameter::Alpha()const{
return _alpha;
}
const real_t &Parameter::Dt()const{
return _dt;
}
const real_t &Parameter::Tend()const{
return _tend;
}
const index_t &Parameter::IterMax()const{
return _itermax;
}
const real_t &Parameter::Eps()const{
return _eps;
}
const real_t &Parameter::Tau()const{
return _tau;
}
const real_t &Parameter::Imax()const{
return _imax;
}
const real_t &Parameter::Jmax()const{
return _jmax;
}
const real_t &Parameter::Dt_value()const{
return _dt_value;
}
const real_t &Parameter::Gx()const{
return _gx;
}
const real_t &Parameter::Gy()const{
return _gy;
}
const real_t &Parameter::T_h()const{
return _T_h;
}
const real_t &Parameter::T_c()const{
return _T_c;
}
const real_t &Parameter::Beta()const{
return _beta;
}
const real_t &Parameter::Pr()const{
return _pr;
}
