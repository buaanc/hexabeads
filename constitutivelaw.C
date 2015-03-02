/*
 * constitutivelaw.cpp
 *
 *  Created on: Sep 29, 2014
 *      Author: miguel
 */

#define PI 3.14159265358979323846


#include "constitutivelaw.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <petscts.h>


#include <cassert>

ConstitutiveLaw::ConstitutiveLaw(){
	 _force_vector = new PetscScalar[4];
	 _partial_derivative_U_force_vector = new PetscScalar[4];
	 _unit_vector = new PetscScalar[2];

	 _partial_derivative_P = 0.0;

	 _internal_partial_derivative_P_bead = new PetscScalar[4];

	 _partial_derivative_alphaMax_coefficient = 0.0;
	 _delta = 0.0;
	 _RStar = 0.0;
	 _EStar = 0.0;
	 _Force = 0.0;
	 _distance = 0.0;
	 _Derivative_Force = 0.0;

	 _alphaMax = 0.0;
	 _alphaP = 0.0;
	 _AreaY = 0.0;
	 _alphaY = 0.0;

	 _problemdata = NULL;

	 _designdata = NULL;

	 _force_coefficients = NULL;
	 _partial_derivative_alphaMax = NULL;
	 _state_eq_partial_derivative_U = NULL;
	 _state_eq_partial_derivative_alphaMax = NULL;
	 _partial_derivative_P_alphaP_Initial = NULL;


	 dh_dalphaMax = 0.0;
	 dh_dU = arma::zeros<arma::colvec>(4);

	 residual_derivative_alphaMax_vec = arma::zeros<arma::colvec>(4);

	 _bead_to = NULL;
	 _bead_from = NULL;

	 _calculate_delta = NULL;

	 _is_delta_calculated = false;
	 _is_matprop_calculated = false;
	 _is_distance_calculated = false;
	 _plastic = false;
	 _is_alphaMax = false;

	 _is_force_coefficients_calculated = false;
	 _is_alphaP_initial_calculated = false;

	 _jacobian = arma::zeros<arma::mat>(4,4);
	 A = arma::zeros<arma::mat>(4,4);
	 A(0,0) = A(1,1) = -1.0;
	 A(2,2) = A(3,3) = -1.0;
	 A(2,0) = A(3,1) = 1.0;
	 A(0,2) = A(1,3) = 1.0;

	 unit_vector_attach = arma::zeros<arma::colvec>(4);

	 // 	alphaP = c1*delta - c2*alphaY + c3*alphaY*exp(-c4*(delta/alphaY-1.0));
	 c1 = 0.9481927511;
	 c2 = 25.934752932;
	 c3 = 24.9865601805;
	 c4 = 0.015;
	 // pAlpha = sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
	 d1 = 2.48;
	 d2 = 1.41;
	 d3 = 0.098;
	 // 	Anorm = pow(alphaNorm,1.137);
	 e1 = 1.137;
	 //  Anorm = (f1*alphaNorm - f2);
	 f1 = 2.37076627;
	 f2 = 59.955374401;
	 // alphaNorm < g1
	 g1 = 177.57;

	 // Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);
	 h1 = 1.35;

	 k2 = 4.0;

	 tol = 1e-10;

}

ConstitutiveLaw::ConstitutiveLaw(ProblemParameters * problemdata, DesignParameters * designdata){
	 _force_vector = new PetscScalar[4];
	 _partial_derivative_U_force_vector = new PetscScalar[4];
	 _unit_vector = new PetscScalar[2];
	 _delta = 0.0;
	 _distance = 0.0;

	 _partial_derivative_alphaMax_coefficient = 0.0;
	 _partial_derivative_P = 0.0;

	 _internal_partial_derivative_P_bead = new PetscScalar[4];

	 _RStar = 0.0;
	 _EStar = 0.0;
	 _Force = 0.0;
	 _Derivative_Force = 0.0;

	 _alphaMax = 0.0;
	 _alphaP = 0.0;
	 _AreaY = 0.0;
	 _alphaY = 0.0;

	 _is_delta_calculated = false;
	 _is_matprop_calculated = false;
	 _is_distance_calculated = false;
	 _plastic = false;
	 _is_alphaMax = false;

	 _is_force_coefficients_calculated = false;

	 _problemdata = problemdata;

	 _designdata = designdata;

	 _force_coefficients = NULL;

	 _partial_derivative_alphaMax = NULL;
	 _state_eq_partial_derivative_U = NULL;
	 _state_eq_partial_derivative_alphaMax = NULL;
	 _partial_derivative_P_alphaP_Initial = NULL;

	 _calculate_delta = NULL;

	 residual_derivative_alphaMax_vec = arma::zeros<arma::colvec>(4);

	 dh_dalphaMax = 0.0;
	 dh_dU = arma::zeros<arma::colvec>(4);

	 _bead_to = NULL;
	 _bead_from = NULL;


	 _jacobian = arma::zeros<arma::mat>(4,4);
	 unit_vector_attach = arma::zeros<arma::colvec>(4);
	 A = arma::zeros<arma::mat>(4,4);
	 A(0,0) = A(1,1) = -1.0;
	 A(2,2) = A(3,3) = -1.0;
	 A(2,0) = A(3,1) = 1.0;
	 A(0,2) = A(1,3) = 1.0;


	 // 	alphaP = c1*delta - c2*alphaY + c3*alphaY*exp(-c4*(delta/alphaY-1.0));
	 c1 = 0.9481927511;
	 c2 = 25.934752932;
	 c3 = 24.9865601805;
	 c4 = 0.015;
	 // pAlpha = sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
	 d1 = 2.48;
	 d2 = 1.41;
	 d3 = 0.098;
	 // 	Anorm = pow(alphaNorm,1.137);
	 e1 = 1.137;
	 //  Anorm = (f1*alphaNorm - f2);
	 f1 = 2.37076627;
	 f2 = 59.955374401;
	 // alphaNorm < g1
	 g1 = 177.57;

	 // Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);
	 h1 = 1.35;

	 k2 = 4.0;

	 tol = 1e-10;
}

ConstitutiveLaw::~ConstitutiveLaw(){
	delete [] _force_vector;
	delete [] _unit_vector;
	delete [] _partial_derivative_U_force_vector;
	delete [] _internal_partial_derivative_P_bead;
}

void ConstitutiveLaw::reset_beads(BeadsData & bead_from, BeadsData & bead_to, PetscInt & analysis,
									PetscScalar & U_i_from, PetscScalar & U_j_from, PetscScalar & U_i_to, PetscScalar & U_j_to){
	_bead_from = &bead_from;
	_bead_to = &bead_to;

	_bead_from->U_i = U_i_from;
	_bead_from->U_j = U_j_from;
	_bead_to->U_i = U_i_to;
	_bead_to->U_j = U_j_to;

	 _is_force_coefficients_calculated = false;

	get_material_properties();
	(*this.*_get_distance_AND_unit_vector)();

	assert((analysis == 2 || analysis == 1) && "Only plastic model, analysis has to be equal to 2 or 1");


	if (analysis == 1){
		_force_coefficients = &ConstitutiveLaw::calculate_force_dummy_plastic_coefficients;

		_partial_derivative_alphaMax= &ConstitutiveLaw::partial_derivative_alphaMax_force_vector_dummy_plastic;

		_state_eq_partial_derivative_alphaMax = &ConstitutiveLaw::state_eq_partial_derivative_alphaMax_dummy_plastic;

		_state_eq_partial_derivative_U = &ConstitutiveLaw::state_eq_partial_derivative_U_dummy_plastic;

		_partial_derivative_P_alphaP_Initial = &ConstitutiveLaw::_partial_derivative_P_alphaP_Initial_dummy_plastic;
	}
	else if (analysis == 2){
		_alphaY = pow(PI/2,2.0)*pow((_problemdata->sigmayielding*1.6/_EStar),2.0)*_RStar;

		_AreaY = PI*_RStar*_alphaY;

		_force_coefficients = &ConstitutiveLaw::calculate_force_plastic_coefficients;

		_partial_derivative_alphaMax= &ConstitutiveLaw::partial_derivative_alphaMax_force_vector_plastic;

		_state_eq_partial_derivative_alphaMax = &ConstitutiveLaw::state_eq_partial_derivative_alphaMax_plastic;

		_state_eq_partial_derivative_U = &ConstitutiveLaw::state_eq_partial_derivative_U_plastic;

		_partial_derivative_P_alphaP_Initial = &ConstitutiveLaw::_partial_derivative_P_alphaP_Initial_plastic;
	}
}

PetscScalar * ConstitutiveLaw::force_vector(){

	assert(_is_distance_calculated);
	if (!_is_force_coefficients_calculated){
		// This is the way to call a function pointer
		(*this.*_force_coefficients)();
		 _is_force_coefficients_calculated = true;
	}

	_force_vector[0] = _Force*_unit_vector[0];
	_force_vector[1] = _Force*_unit_vector[1];

	_force_vector[2] = -1.0*_Force*_unit_vector[0];
	_force_vector[3] = -1.0*_Force*_unit_vector[1];

	return _force_vector;
}

/* ----------------------------------------------------------
 *
 *
 *  Partial derivatives with respect to the implicit variable
 *
 *
 * ----------------------------------------------------------
 */
PetscScalar * ConstitutiveLaw::partial_derivative_U_force_vector(){

	if (!_is_force_coefficients_calculated){
		// This is the way to call a function pointer
		(*this.*_force_coefficients)();
		 _is_force_coefficients_calculated = true;
	}

	_partial_derivative_U_force_vector[0] = _Derivative_Force*_unit_vector[0];
	_partial_derivative_U_force_vector[1] = _Derivative_Force*_unit_vector[1];

	_partial_derivative_U_force_vector[2] = -1.0*_Derivative_Force*_unit_vector[0];
	_partial_derivative_U_force_vector[3] = -1.0*_Derivative_Force*_unit_vector[1];

	return _partial_derivative_U_force_vector;
}

/* ---------------------------------------------------
 *
 *
 *  Partial derivatives with respect to the parameter
 *
 *
 * ----------------------------------------------------
 */

PetscScalar * ConstitutiveLaw::internal_force_partial_derivative_P_bead(){

	// Only for plastic analysis, for now
	assert(_is_alphaMax);

	/* Calculate partial derivative of alphaP w.r.t. alphaMax */

	assert(_is_alphaP_initial_calculated);

	/* Call the force coefficients */

	if (!_is_force_coefficients_calculated){
		// This is the way to call a function pointer
		(*this.*_force_coefficients)();
		 _is_force_coefficients_calculated = true;
	}

	_partial_derivative_P = _Derivative_Force * _derivative_alphaP_Initial;

	assert(_is_distance_calculated);

	// We assemble the internal force vector for the "from" bead or the "to" bead
	// depending on the value of bead

	_internal_partial_derivative_P_bead[0] = _partial_derivative_P*_unit_vector[0];
	_internal_partial_derivative_P_bead[1] = _partial_derivative_P*_unit_vector[1];
	_internal_partial_derivative_P_bead[2] = -1.0*_partial_derivative_P*_unit_vector[0];
	_internal_partial_derivative_P_bead[3] = -1.0*_partial_derivative_P*_unit_vector[1];


	return _internal_partial_derivative_P_bead;
}

PetscScalar & ConstitutiveLaw::partial_derivative_P(){

	// Only for plastic analysis, for now
	assert(_is_alphaMax);

	/* Calculate partial derivative of alphaP w.r.t. alphaMax */

	assert(_is_alphaP_initial_calculated);

	/* Call the force coefficients */

	if (!_is_force_coefficients_calculated){
		// This is the way to call a function pointer
		(*this.*_force_coefficients)();
		 _is_force_coefficients_calculated = true;
	}

	_partial_derivative_P = _Derivative_Force * _derivative_alphaP_Initial;

	return _partial_derivative_P;
}



void ConstitutiveLaw::_partial_derivative_P_alphaP_Initial_plastic(PetscScalar & alphaMax_Initial) {

	_alphaP_Initial = c1*alphaMax_Initial - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaMax_Initial/_alphaY-1.0));

	if (_alphaP_Initial < 0.0)
		_alphaP_Initial = 0;

	_derivative_alphaP_Initial = c1 - c4*c3*exp(-c4*(alphaMax_Initial/_alphaY-1.0));
}

void ConstitutiveLaw::_partial_derivative_P_alphaP_Initial_dummy_plastic(PetscScalar & alphaMax_Initial) {

	_alphaP_Initial = alphaMax_Initial*(k2 - k1)/ k2;

	_derivative_alphaP_Initial = (k2 - k1)/ k2;;
}

/* ---------------------------------------------------
 *
 *
 *  Partial derivatives with respect to alphaMax
 *
 *
 * ----------------------------------------------------
 */

arma::colvec & ConstitutiveLaw::residual_derivative_alphaMax(){
	assert(_is_distance_calculated);

	// This is the way to call a function pointer
	(*this.*_partial_derivative_alphaMax)();



	residual_derivative_alphaMax_vec(0) = _partial_derivative_alphaMax_coefficient*_unit_vector[0];
	residual_derivative_alphaMax_vec(1) = _partial_derivative_alphaMax_coefficient*_unit_vector[1];

	residual_derivative_alphaMax_vec(2) = -1.0*_partial_derivative_alphaMax_coefficient*_unit_vector[0];
	residual_derivative_alphaMax_vec(3) = -1.0*_partial_derivative_alphaMax_coefficient*_unit_vector[1];

	return residual_derivative_alphaMax_vec;
}

PetscScalar & ConstitutiveLaw::partial_derivative_alphaMax_force_vector(){

	// This is the way to call a function pointer
	(*this.*_partial_derivative_alphaMax)();

	return _partial_derivative_alphaMax_coefficient;
}

void ConstitutiveLaw::partial_derivative_alphaMax_force_vector_dummy_plastic(){
	// Only for plastic analysis
	assert(_problemdata->analysis == 1);
	assert(_is_alphaMax);

	// Negative deltas produce negative forces --> compression
	// Positive deltas produce no forces --> tension

	k1 = (_bead_from->density_design + _bead_to->density_design)/2.0;

	_alphaP = _alphaMax*(k2 - k1)/ k2;


	PetscReal Felastic, FPlastic;


	Felastic = k2*(_delta - _alphaP);

	FPlastic = k1*_delta;

	_partial_derivative_alphaMax_coefficient = 0.0;

//	std::cout<<"Felastic = "<<std::setprecision(10)<<Felastic<<std::endl;
//
//	std::cout<<"FPlastic = "<<std::setprecision(10)<<FPlastic<<std::endl;
//
//	std::cout<<"_delta = "<<_delta<<std::endl;
//
//	std::cout<<"_alphaMax = "<<_alphaMax<<std::endl;

	if (_delta >= _alphaP){
		if (Felastic < FPlastic - tol){
			_partial_derivative_alphaMax_coefficient = k2*(- (k2 - k1)/ k2);

		}
		else {
			_partial_derivative_alphaMax_coefficient = 0.0;
		}
	}
	else{
		_partial_derivative_alphaMax_coefficient = 0.0;
	}

	//std::cout<<"_partial_derivative_alphaMax_coefficient = "<<_partial_derivative_alphaMax_coefficient<<std::endl;


}

void ConstitutiveLaw::partial_derivative_alphaMax_force_vector_plastic(){

	// Only for plastic analysis
	assert(_problemdata->analysis == 2);
	assert(_is_alphaMax);

	PetscScalar pAlpha, Anorm;
	PetscScalar alphaNorm = _alphaMax/_alphaY;

	PetscScalar DpAlphaDalphaM, DAnormDalphaM, DFmaxElasticDalphaM, DalphaPDalphaM;

	PetscScalar FMaxElastic, Felastic, FPlastic;

	// _alphaMax refers to the level reached in the previous timestep
	// if we are still in the elastic regime, alphaMax = 0
	_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
	if (_alphaP < 0.0)
		_alphaP = 0.0;

	_partial_derivative_alphaMax_coefficient = 0.0;

	if (_delta >= _alphaP)
	{
		if (_delta <= _alphaY && _alphaP == 0.0)
		{
			_partial_derivative_alphaMax_coefficient = 0.0;
		}
		else
		{
			if (_alphaP == 0)
			{
				_partial_derivative_alphaMax_coefficient = 0.0;
			}
			else{
				//Elastic force
				alphaNorm = _alphaMax/_alphaY;
				// We update alphaP because it changes when we change alphaY
				_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
				if (_alphaP<0.0)
					_alphaP = 0.0;


				pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
				if (alphaNorm<g1)
					Anorm = pow(alphaNorm,e1);
				else
					Anorm = (f1*alphaNorm - f2);
				FMaxElastic = pAlpha*Anorm;
				Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);

				//Plastic force
				alphaNorm = _delta/_alphaY;
				_alphaP = c1*_delta - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
				if (_alphaP<0.0)
					_alphaP = 0.0;

				pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
				if (alphaNorm<g1)
				{
					Anorm = pow(alphaNorm,e1);
				}
				else
				{
					Anorm = (f1*alphaNorm - f2);
				}
				FPlastic = pAlpha*Anorm;

				//Yield criteria
				if (Felastic < FPlastic - tol)
				{
					//std::cout<<"Elastic"<<std::endl;
					// We need to get back the _alphaP corresponding to _alphaMax instead of _delta
					alphaNorm = _alphaMax/_alphaY;
					_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
					DalphaPDalphaM = c1 - c4*c3*exp(-c4*(alphaNorm-1.0));

					if (_alphaP<0.0)
						_alphaP = 0.0;

					pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
					DpAlphaDalphaM = _problemdata->sigmayielding*_AreaY*(d2*d3*exp(-d3*(alphaNorm-1.0))/_alphaY);
					if (alphaNorm<g1)
					{
						Anorm = pow(alphaNorm,e1);
						DAnormDalphaM = 1.0/_alphaY*e1*pow(alphaNorm,e1-1.0);
					}
					else
					{
						Anorm = (f1*alphaNorm - f2);
						DAnormDalphaM = f1/_alphaY;
					}

					DFmaxElasticDalphaM = DpAlphaDalphaM*Anorm + DAnormDalphaM*pAlpha;

					_partial_derivative_alphaMax_coefficient = DFmaxElasticDalphaM*pow((_delta - _alphaP)/(_alphaMax - _alphaP),h1)
							+ h1*FMaxElastic*pow((_delta - _alphaP)/(_alphaMax - _alphaP),h1 - 1.0)*
							(DalphaPDalphaM*(_delta - _alphaMax) - (_delta - _alphaP))/(pow(_alphaMax - _alphaP,2.0));

				}
				else
				{
					//std::cout<<"Plastic"<<std::endl;
					_partial_derivative_alphaMax_coefficient = 0.0;
				}
			}
		}
	}
	else{
		_partial_derivative_alphaMax_coefficient = 0.0;
	}


}


/* ---------------------------------------------------
 *
 *
 *  Force and derivative force coefficients
 *
 *
 * ----------------------------------------------------
 */
PetscScalar & ConstitutiveLaw::force_coefficient(){


	if (!_is_force_coefficients_calculated){
		// This is the way to call a function pointer
		(*this.*_force_coefficients)();
		 _is_force_coefficients_calculated = true;
	}


	return _Force;
}

void ConstitutiveLaw::calculate_force_dummy_plastic_coefficients(){

	// Need to check that the material properties and delta have been calculated
	assert(_is_delta_calculated);
	assert(_is_matprop_calculated);

	// Negative deltas produce negative forces --> compression
	// Positive deltas produce no forces --> tension

	k1 = (_bead_from->density_design + _bead_to->density_design)/2.0;

	_alphaP = _alphaMax*(k2 - k1)/ k2;


	PetscReal Felastic, FPlastic;

#ifdef DEBUG
			std::cout<<"_delta = "<<_delta<<std::endl;
			std::cout<<"_alphaMax = "<<_alphaMax<<std::endl;
			std::cout<<"_alphaP = "<<_alphaP<<std::endl;
#endif


	Felastic = k2*(_delta - _alphaP);

	FPlastic = k1*_delta;

	if (_delta >= _alphaP){
		if (Felastic < FPlastic - tol){
			_Force = Felastic;
			_Derivative_Force = k2;

#ifdef DEBUG
			std::cout<<"AQUI ELASTIC DUMMY"<<std::endl;
#endif

		}
		else {
			_Force = FPlastic;
			_alphaMax = _delta;
			_Derivative_Force = k1;

#ifdef DEBUG
			std::cout<<"AQUI PLASTIC DUMMY"<<std::endl;
#endif


		}
	}
	else{
		_Force = 0.0;
		_Derivative_Force = 0.0;
	}
}


void ConstitutiveLaw::calculate_force_plastic_coefficients(){
	assert(_is_alphaMax);

	PetscScalar DP, pAlpha, Anorm, DA;
	PetscScalar alphaNorm = _alphaMax/_alphaY;

	PetscScalar FMaxElastic, Felastic, FPlastic;

	// _alphaMax refers to the level reached in the previous timestep
	// if we are still in the elastic regime, alphaMax = 0
	_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
	if (_alphaP < 0.0)
		_alphaP = 0.0;

	if (_delta >= _alphaP)
	{
		if (_delta <= _alphaY && _alphaP == 0.0)
		{

			PetscScalar ke = 4.0/3.0* _EStar * sqrt(_RStar);

			_Force = ke*pow(fabs(_delta),1.5);

			_Derivative_Force = 1.5*ke*pow(fabs(_delta),0.5);
		}
		else
		{
			if (_alphaP == 0)
			{
				// Update the alphaP
				alphaNorm = _delta/_alphaY;
				_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
				if (_alphaP<0.0)
					_alphaP = 0.0;

				pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
				DP = _problemdata->sigmayielding*_AreaY*d2*d3*exp(-d3*(alphaNorm-1.0))/_alphaY;
				if (alphaNorm<g1)
				{
					Anorm = pow(alphaNorm,e1);
					DA = 1.0/_alphaY*e1*pow(alphaNorm,e1 - 1.0);
				}
				else
				{
					Anorm = (f1*alphaNorm - f2);
					DA = f1/_alphaY;
				}
				_Force = pAlpha*Anorm;
				_Derivative_Force = DP*Anorm + DA*pAlpha;

				// Update _alphaMax
				_alphaMax = _delta;
			}
			else{
				//Elastic force
				alphaNorm = _alphaMax/_alphaY;
				// We update alphaP because it changes when we change alphaY
				_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
				if (_alphaP<0.0)
					_alphaP = 0.0;


				pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
				if (alphaNorm<g1)
					Anorm = pow(alphaNorm,e1);
				else
					Anorm = (f1*alphaNorm - f2);
				FMaxElastic = pAlpha*Anorm;
				Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);

				//Plastic force
				alphaNorm = _delta/_alphaY;
				_alphaP = c1*_delta - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
				if (_alphaP<0.0)
					_alphaP = 0.0;

				pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
				DP = _problemdata->sigmayielding*_AreaY*d2*d3*exp(-d3*(alphaNorm-1.0))/_alphaY;
				if (alphaNorm<g1)
				{
					Anorm = pow(alphaNorm,e1);
					DA = 1.0/_alphaY*e1*pow(alphaNorm,e1 - 1.0);
				}
				else
				{
					Anorm = (f1*alphaNorm - f2);
					DA = f1/_alphaY;
				}
				FPlastic = pAlpha*Anorm;

				//Yield criteria
				if (Felastic < FPlastic - tol)
				{
					// We need to get back the _alphaP corresponding to _alphaMax instead of _delta
					alphaNorm = _alphaMax/_alphaY;
					_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
					if (_alphaP<0.0)
						_alphaP = 0.0;
					_Force = Felastic;
					_Derivative_Force = FMaxElastic*h1*1.0/(_alphaMax - _alphaP)*pow((_delta - _alphaP)/(_alphaMax - _alphaP),h1 - 1.0);


				}
				else
				{
					//Update alphaP
					alphaNorm = _delta/_alphaY;
					_alphaP = c1*_delta - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
					if (_alphaP<0.0)
						_alphaP = 0.0;

					_Force = FPlastic;
					_Derivative_Force = DP*Anorm + DA*pAlpha;

					// Update _alphaMax
					_alphaMax = _delta;
				}
			}
		}
	}
	else{
		_Force = 0.0;
		_Derivative_Force = 0.0;
	}
}


/* ---------------------------------------------------
 *
 *
 *  				Printing information
 *
 *
 * ----------------------------------------------------
 */
void ConstitutiveLaw::print_force_vector(){

	std::cout<<"F(0) = "<<_force_vector[0]<<std::endl;
	std::cout<<"F(1) = "<<_force_vector[1]<<std::endl;
	std::cout<<"F(2) = "<<_force_vector[2]<<std::endl;
	std::cout<<"F(3) = "<<_force_vector[3]<<std::endl;

}

void ConstitutiveLaw::print_jacobian_matrix(){

	std::cout<<"Derivative Force = "<<_Derivative_Force<<std::endl;
	std::cout<<"Bead from = "<<_bead_from->node_id<<std::endl;
	std::cout<<"Bead to = "<<_bead_to->node_id<<std::endl;

	std::cout<<"Printing local jacobian"<<std::endl;
	_jacobian.print();
}

/* ---------------------------------------------------
 *
 *
 *  	Jacobian and Matrix-vector products
 *
 *
 * ----------------------------------------------------
 */

arma::mat & ConstitutiveLaw::jacobian_matrix(){

	assert(_is_delta_calculated);

	// This is the way to call a function pointer
	(*this.*_jacobian_function)();

	return _jacobian;
}

void ConstitutiveLaw::jacobian_matrix_small_def(){
	if (!_is_force_coefficients_calculated){
		// This is the way to call a function pointer
		(*this.*_force_coefficients)();
		 _is_force_coefficients_calculated = true;
	}


	_jacobian(0,0) = _Derivative_Force*_unit_vector[0]*_unit_vector[0];
	_jacobian(0,1) = _Derivative_Force*_unit_vector[0]*_unit_vector[1];
	_jacobian(0,2) = _Derivative_Force*_unit_vector[0]*(-1.0)*_unit_vector[0];
	_jacobian(0,3) = _Derivative_Force*_unit_vector[0]*(-1.0)*_unit_vector[1];

	_jacobian(1,0) = _Derivative_Force*_unit_vector[1]*_unit_vector[0];
	_jacobian(1,1) = _Derivative_Force*_unit_vector[1]*_unit_vector[1];
	_jacobian(1,2) = _Derivative_Force*_unit_vector[1]*(-1.0)*_unit_vector[0];
	_jacobian(1,3) = _Derivative_Force*_unit_vector[1]*(-1.0)*_unit_vector[1];

	_jacobian(2,0) = _Derivative_Force*(-1.0)*_unit_vector[0]*_unit_vector[0];
	_jacobian(2,1) = _Derivative_Force*(-1.0)*_unit_vector[0]*_unit_vector[1];
	_jacobian(2,2) = _Derivative_Force*(-1.0)*_unit_vector[0]*(-1.0)*_unit_vector[0];
	_jacobian(2,3) = _Derivative_Force*(-1.0)*_unit_vector[0]*(-1.0)*_unit_vector[1];

	_jacobian(3,0) = _Derivative_Force*(-1.0)*_unit_vector[1]*_unit_vector[0];
	_jacobian(3,1) = _Derivative_Force*(-1.0)*_unit_vector[1]*_unit_vector[1];
	_jacobian(3,2) = _Derivative_Force*(-1.0)*_unit_vector[1]*(-1.0)*_unit_vector[0];
	_jacobian(3,3) = _Derivative_Force*(-1.0)*_unit_vector[1]*(-1.0)*_unit_vector[1];
}

void ConstitutiveLaw::jacobian_matrix_large_def(){
	if (!_is_force_coefficients_calculated){
		// This is the way to call a function pointer
		(*this.*_force_coefficients)();
		 _is_force_coefficients_calculated = true;
	}

	 unit_vector_attach(0) = _unit_vector[0];
	 unit_vector_attach(1) = _unit_vector[1];
	 unit_vector_attach(2) = -(1.0)* _unit_vector[0];
	 unit_vector_attach(3) = -(1.0)*_unit_vector[1];

	/* First contribution of the stiffness matrix
	 *     dF 								  dF
	 *  ------- * t_x \crossbox D \alpha = ------- * t_x \crossbox t_x
	 *  d \alpha 							d \alpha
	 */

	 _jacobian = _Derivative_Force*unit_vector_attach*arma::trans(unit_vector_attach);
	/* Second contribution of the stiffness matrix
	 *         d t		 1
	 *  F * ------- = F --- (I - t \crossbox t)
	 *         d U		 d
	 */

	 _jacobian += _Force/_distance * (A + unit_vector_attach*arma::trans(unit_vector_attach));

}

void ConstitutiveLaw::jacobian_vector_product(arma::colvec & vector, arma::colvec & result){

	arma::mat & jacobian = jacobian_matrix();


	result = arma::trans(jacobian)*vector;

}

void ConstitutiveLaw::jacobian_alphaMax_vector_product(arma::colvec & vector, PetscScalar & result){

	arma::colvec & jacobian_alpha = residual_derivative_alphaMax();


	result = arma::dot(jacobian_alpha,vector);

}

void ConstitutiveLaw::check_jacobian(){

	arma::mat & jacobian_analytical = jacobian_matrix();

	PetscReal dt = 1e-8;

	arma::colvec jacobian_fd = arma::zeros(4);

	for (unsigned int i = 0; i < 4; i++){


		if (i==0)
			_bead_from->U_i += dt;
		else if (i==1)
			_bead_from->U_j += dt;
		else if (i==2)
			_bead_to->U_i += dt;
		else if (i==3)
			_bead_to->U_j += dt;


		/*
		 * Reinitialize delta, distance and unit vector
		 */
		_is_force_coefficients_calculated = false;
		(*this.*_get_distance_AND_unit_vector)();
		calculate_delta();

		/*
		 * Call force vector
		 */

		PetscScalar * force_vector_fd = force_vector();

		jacobian_fd(0) = force_vector_fd[0];
		jacobian_fd(1) = force_vector_fd[1];
		jacobian_fd(2) = force_vector_fd[2];
		jacobian_fd(3) = force_vector_fd[3];


		if (i==0)
			_bead_from->U_i -= 2.0*dt;
		else if (i==1)
			_bead_from->U_j -= 2.0*dt;
		else if (i==2)
			_bead_to->U_i -= 2.0*dt;
		else if (i==3)
			_bead_to->U_j -= 2.0*dt;

		_is_force_coefficients_calculated = false;
		(*this.*_get_distance_AND_unit_vector)();
		calculate_delta();


		PetscScalar * force_vector_fd_dos = force_vector();

		jacobian_fd(0) -= force_vector_fd_dos[0];
		jacobian_fd(1) -= force_vector_fd_dos[1];
		jacobian_fd(2) -= force_vector_fd_dos[2];
		jacobian_fd(3) -= force_vector_fd_dos[3];

		jacobian_fd *= 1/(2.0*dt);

		/*
		 * Comparing the column
		 */

		std::cout<<"F(0) FD = "<<jacobian_fd(0)<<" J(0,"<<i<<") = "<<jacobian_analytical(0,i)<<std::endl;
		std::cout<<"F(1) FD = "<<jacobian_fd(1)<<" J(1,"<<i<<") = "<<jacobian_analytical(1,i)<<std::endl;
		std::cout<<"F(2) FD = "<<jacobian_fd(2)<<" J(2,"<<i<<") = "<<jacobian_analytical(2,i)<<std::endl;
		std::cout<<"F(3) FD = "<<jacobian_fd(3)<<" J(3,"<<i<<") = "<<jacobian_analytical(3,i)<<std::endl;


		if (i==0)
			_bead_from->U_i -= dt;
		else if (i==1)
			_bead_from->U_j -= dt;
		else if (i==2)
			_bead_to->U_i -= dt;
		else if (i==3)
			_bead_to->U_j -= dt;



	}


}
/* ---------------------------------------------------
 *
 *
 *  	State equation partial derivatives
 *
 *		They return:
 *			- dh_dalphaMax: Scalar, it's either 1 or zero. Diagonal member for the correspondent row
 *			- dh_dU: Vector, dimensions 1 by n_dofs per element (4)
 * ----------------------------------------------------
 */

void ConstitutiveLaw::state_eq_partial_derivative_alphaMax_plastic(){

	assert(_is_alphaMax);

		PetscScalar pAlpha, Anorm;
		PetscScalar alphaNorm = _alphaMax/_alphaY;

		PetscScalar FMaxElastic, Felastic, FPlastic;

		// _alphaMax refers to the level reached in the previous timestep
		// if we are still in the elastic regime, alphaMax = 0
		_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
		if (_alphaP < 0.0)
			_alphaP = 0.0;

		if (_delta >= _alphaP)
		{
			if (_delta <= _alphaY && _alphaP == 0.0)
			{
				dh_dalphaMax = 0;
			}
			else
			{
				if (_alphaP == 0)
				{
					dh_dalphaMax = 0.0;
				}
				else{
					//Elastic force
					alphaNorm = _alphaMax/_alphaY;
					// We update alphaP because it changes when we change alphaY
					_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
					if (_alphaP<0.0)
						_alphaP = 0.0;


					pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
					if (alphaNorm<g1)
						Anorm = pow(alphaNorm,e1);
					else
						Anorm = (f1*alphaNorm - f2);
					FMaxElastic = pAlpha*Anorm;
					Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);

					//Plastic force
					alphaNorm = _delta/_alphaY;
					_alphaP = c1*_delta - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
					if (_alphaP<0.0)
						_alphaP = 0.0;

					pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
					if (alphaNorm<g1)
					{
						Anorm = pow(alphaNorm,e1);
					}
					else
					{
						Anorm = (f1*alphaNorm - f2);
					}
					FPlastic = pAlpha*Anorm;

					//Yield criteria
					if (Felastic < FPlastic - tol)
					{
						dh_dalphaMax = 1.0;
					}
					else
					{
						dh_dalphaMax = 0.0;
					}
				}
			}
		}
		else{
			dh_dalphaMax = 1.0;
		}
}

void ConstitutiveLaw::state_eq_partial_derivative_alphaMax_dummy_plastic(){

	// Need to check that the material properties and delta have been calculated
	assert(_is_delta_calculated);
	assert(_is_matprop_calculated);

	// Negative deltas produce negative forces --> compression
	// Positive deltas produce no forces --> tension

	k1 = (_bead_from->density_design + _bead_to->density_design)/2.0;

	_alphaP = _alphaMax*(k2 - k1)/ k2;


	PetscReal Felastic, FPlastic;


	Felastic = k2*(_delta - _alphaP);

	FPlastic = k1*_delta;


	if (_delta >= _alphaP){
		if (Felastic < FPlastic - tol){
			dh_dalphaMax = 1.0;

		}
		else {
			dh_dalphaMax = 0.0;
		}
	}
	else{
			dh_dalphaMax = 1.0;
	}
}


void ConstitutiveLaw::state_eq_partial_derivative_U_plastic(){

	assert(_is_alphaMax);
	assert(_is_distance_calculated);

		PetscScalar pAlpha, Anorm;
		PetscScalar alphaNorm = _alphaMax/_alphaY;

		PetscScalar FMaxElastic, Felastic, FPlastic;

		// _alphaMax refers to the level reached in the previous timestep
		// if we are still in the elastic regime, alphaMax = 0
		_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
		if (_alphaP < 0.0)
			_alphaP = 0.0;

		if (_delta >= _alphaP)
		{
			if (_delta <= _alphaY && _alphaP == 0.0)
			{
				dh_dU(0) = 0.0;
				dh_dU(1) = 0.0;

				dh_dU(2) = 0.0;
				dh_dU(3) = 0.0;
			}
			else
			{
				if (_alphaP == 0)
				{
					dh_dU(0) = _unit_vector[0];
					dh_dU(1) = _unit_vector[1];

					dh_dU(2) = -1.0*_unit_vector[0];
					dh_dU(3) = -1.0*_unit_vector[1];
				}
				else{
					//Elastic force
					alphaNorm = _alphaMax/_alphaY;
					// We update alphaP because it changes when we change alphaY
					_alphaP = c1*_alphaMax - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
					if (_alphaP<0.0)
						_alphaP = 0.0;


					pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
					if (alphaNorm<g1)
						Anorm = pow(alphaNorm,e1);
					else
						Anorm = (f1*alphaNorm - f2);
					FMaxElastic = pAlpha*Anorm;
					Felastic = FMaxElastic*pow(((_delta - _alphaP)/(_alphaMax - _alphaP)),h1);

					//Plastic force
					alphaNorm = _delta/_alphaY;
					_alphaP = c1*_delta - c2*_alphaY + c3*_alphaY*exp(-c4*(alphaNorm-1.0));
					if (_alphaP<0.0)
						_alphaP = 0.0;

					pAlpha = _problemdata->sigmayielding*_AreaY*(d1 - d2*exp(-d3*(alphaNorm-1.0)));
					if (alphaNorm<g1)
					{
						Anorm = pow(alphaNorm,e1);
					}
					else
					{
						Anorm = (f1*alphaNorm - f2);
					}
					FPlastic = pAlpha*Anorm;

					//Yield criteria
					if (Felastic < FPlastic - tol)
					{
						dh_dU(0) = 0.0;
						dh_dU(1) = 0.0;

						dh_dU(2) = 0.0;
						dh_dU(3) = 0.0;
					}
					else
					{
						dh_dU(0) = _unit_vector[0];
						dh_dU(1) = _unit_vector[1];

						dh_dU(2) = -1.0*_unit_vector[0];
						dh_dU(3) = -1.0*_unit_vector[1];
					}
				}
			}
		}
		else{
			dh_dU(0) = 0.0;
			dh_dU(1) = 0.0;

			dh_dU(2) = 0.0;
			dh_dU(3) = 0.0;
		}
}



void ConstitutiveLaw::state_eq_partial_derivative_U_dummy_plastic(){
	// Need to check that the material properties and delta have been calculated
	assert(_is_delta_calculated);
	assert(_is_matprop_calculated);

	// Negative deltas produce negative forces --> compression
	// Positive deltas produce no forces --> tension

	k1 = (_bead_from->density_design + _bead_to->density_design)/2.0;

	_alphaP = _alphaMax*(k2 - k1)/ k2;


	PetscReal Felastic, FPlastic;


	Felastic = k2*(_delta - _alphaP);

	FPlastic = k1*_delta;

	if (_delta >= _alphaP){
		if (Felastic < FPlastic - tol){
			dh_dU(0) = 0.0;
			dh_dU(1) = 0.0;

			dh_dU(2) = 0.0;
			dh_dU(3) = 0.0;
		}
		else {
			dh_dU(0) = _unit_vector[0];
			dh_dU(1) = _unit_vector[1];

			dh_dU(2) = -1.0*_unit_vector[0];
			dh_dU(3) = -1.0*_unit_vector[1];
		}
	}
	else{
		dh_dU(0) = 0.0;
		dh_dU(1) = 0.0;

		dh_dU(2) = 0.0;
		dh_dU(3) = 0.0;
	}


}


void ConstitutiveLaw::state_eq_product_alphaMax(const PetscScalar & input, PetscScalar & result){

	// This is the way to call a function pointer
	(*this.*_state_eq_partial_derivative_alphaMax)();

	result = dh_dalphaMax*input;

}

void ConstitutiveLaw::state_eq_vector_product_U(const PetscScalar & gamma, arma::colvec & result){

	// This is the way to call a function pointer
	(*this.*_state_eq_partial_derivative_U)();

	result = gamma*dh_dU;

}


/* ---------------------------------------------------
 *
 *
 *  				Set values
 *
 *
 * ----------------------------------------------------
 */

void ConstitutiveLaw::set_alphaMax(const PetscScalar & alphaMax){
	_alphaMax = alphaMax;
	_is_alphaMax = true;
}


void ConstitutiveLaw::calculate_delta(){

	assert(_is_distance_calculated);


	(*this.*_calculate_delta)();


}

void ConstitutiveLaw::calculate_delta_small_def(){

	// For small deformations, delta is calculated as the dot product of (t,-t) y el vector de desplazamientos
	// t is the unit vector

	_delta = 0.0;

	_delta += _unit_vector[0]*_bead_from->U_i;
	_delta += _unit_vector[1]*_bead_from->U_j;

	_delta += -1.0*_unit_vector[0]*_bead_to->U_i;
	_delta += -1.0*_unit_vector[1]*_bead_to->U_j;

	_is_delta_calculated = true;

	/* Add the alphaPInitial contribution */

	assert(_is_alphaP_initial_calculated);

	_delta += _alphaP_Initial;

}

void ConstitutiveLaw::calculate_delta_large_def(){

	_delta = 0.0;

	PetscScalar inc_X, inc_Y, coord_dist_X, coord_dist_Y;

	coord_dist_X = -_bead_from->coordx + _bead_to->coordx;
	coord_dist_Y = -_bead_from->coordy + _bead_to->coordy;


	inc_X = (-_bead_from->coordx + _bead_to->coordx) + (-_bead_from->U_i + _bead_to->U_i);
	inc_Y = (-_bead_from->coordy + _bead_to->coordy) + (-_bead_from->U_j + _bead_to->U_j);

	_delta = sqrt(pow(coord_dist_X,2.0) + pow(coord_dist_Y,2.0)) - sqrt(pow(inc_X,2.0) + pow(inc_Y,2.0));

	_is_delta_calculated = true;

	/* Add the alphaPInitial contribution */

	assert(_is_alphaP_initial_calculated);

	_delta += _alphaP_Initial;

}

void ConstitutiveLaw::vector_beads_small_def(){
	PetscScalar inc_X, inc_Y;

	inc_X = (-_bead_from->coordx + _bead_to->coordx);
	inc_Y = (-_bead_from->coordy + _bead_to->coordy);

	_distance = sqrt(pow(inc_X,2.0) + pow(inc_Y, 2.0));

	_unit_vector[0] = inc_X / _distance;
	_unit_vector[1] = inc_Y / _distance;

	_is_distance_calculated = true;
}

void ConstitutiveLaw::vector_beads_large_def(){
	PetscScalar inc_X, inc_Y;

	inc_X = (-_bead_from->coordx + _bead_to->coordx) + (-_bead_from->U_i + _bead_to->U_i);
	inc_Y = (-_bead_from->coordy + _bead_to->coordy) + (-_bead_from->U_j + _bead_to->U_j);

	_distance = sqrt(pow(inc_X,2.0) + pow(inc_Y, 2.0));

	_unit_vector[0] = inc_X / _distance;
	_unit_vector[1] = inc_Y / _distance;

	_is_distance_calculated = true;

}


void ConstitutiveLaw::get_material_properties(){

	PetscScalar E_to, E_from;

	if (_bead_to->type==2)	{ //Wall bead
		_RStar = _bead_from->R;
	}
	else if (_bead_from->type==2)  { //Wall bead
		_RStar = _bead_to->R;
	}
	else {
		_RStar = (_bead_to->R*_bead_from->R)/(_bead_to->R+_bead_from->R);
	}

	E_from = _bead_from->E;
	E_to = _bead_to->E;

	_EStar = (E_to*E_from)/(E_to*(1.0-pow(_bead_from->nu,2.0))+E_from*(1-pow(_bead_to->nu,2.0)));

	_is_matprop_calculated = true;
}



