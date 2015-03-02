/*
 * constitutivelaw.h
 *
 *  Created on: Sep 29, 2014
 *      Author: miguel
 */

#ifndef CONSTITUTIVELAW_H_
#define CONSTITUTIVELAW_H_

#include "beadsdata.h"
#include "problemparameters.h"
#include "designparameters.h"
#include "petscsys.h"
#include "armadillo"

#include <cassert>
class ConstitutiveLaw
{
	public:

		ConstitutiveLaw();
		ConstitutiveLaw(ProblemParameters * problemdata, DesignParameters * designdata);
		~ConstitutiveLaw();

		PetscReal * force_vector();

		PetscReal & force_coefficient();

		void calculate_delta();

		/* ---------------------------------------------------
		 *
		 *
		 *  	Setters
		 *
		 *
		 * ----------------------------------------------------
		 */

		void reset_beads(BeadsData & bead_to, BeadsData & bead_from, PetscInt & analysis,
				PetscScalar & U_i_from, PetscScalar & U_j_from, PetscScalar & U_i_to, PetscScalar & U_j_to);
		void set_alphaMax(const PetscReal & alphaMax);

		void set_data(ProblemParameters * problemdata, DesignParameters * designdata){
			_problemdata = problemdata;
			_designdata = designdata;
		}

		void set_large_deformations(bool large_deformations){
			if (large_deformations){
				_jacobian_function = &ConstitutiveLaw::jacobian_matrix_large_def;
				_get_distance_AND_unit_vector = &ConstitutiveLaw::vector_beads_large_def;
				_calculate_delta = &ConstitutiveLaw::calculate_delta_large_def;
			}
			else{
				_jacobian_function = &ConstitutiveLaw::jacobian_matrix_small_def;
				_get_distance_AND_unit_vector = &ConstitutiveLaw::vector_beads_small_def;
				_calculate_delta = &ConstitutiveLaw::calculate_delta_small_def;
			}
		}


		void set_alphaP_initial(PetscScalar & alphaMax_Initial){
			/* We need to have calculated alphaY*/
			assert(_is_matprop_calculated);

			(*this.*_partial_derivative_P_alphaP_Initial)(alphaMax_Initial);

			_is_alphaP_initial_calculated = true;

		}



		/* ---------------------------------------------------
		 *
		 *
		 *  	Getters
		 *
		 *
		 * ----------------------------------------------------
		 */

		PetscReal & get_alphaMax(){ return _alphaMax;}

		PetscReal & get_delta(){ return _delta;}

		PetscScalar & get_alphaP_initial(){

			/* We need to have calculated _alphaP_Initial*/
			assert(_is_alphaP_initial_calculated);

			return _alphaP_Initial;

		}

		PetscScalar & get_derivative_force_coefficient() {

			assert(_is_force_coefficients_calculated);

			return _Derivative_Force;
		}




		/* ---------------------------------------------------
		 *
		 *
		 *  	Jacobian and Matrix-vector products
		 *
		 *
		 * ----------------------------------------------------
		 */

		arma::mat & jacobian_matrix();

		void jacobian_vector_product(arma::colvec & vector, arma::colvec & result);

		void state_eq_product_alphaMax(const PetscScalar & input, PetscScalar & result);

		void state_eq_vector_product_U(const PetscReal & gamma, arma::colvec & result);

		void jacobian_alphaMax_vector_product(arma::colvec & vector, PetscScalar & result);

		void check_jacobian();



		/* ---------------------------------------------------
		 *
		 *
		 *  Partial derivatives
		 *
		 *
		 * ----------------------------------------------------
		 */


		PetscReal * partial_derivative_U_force_vector();

		PetscReal & partial_derivative_P();

		PetscReal & partial_derivative_alphaMax_force_vector();

		arma::colvec & residual_derivative_alphaMax();


		PetscReal * internal_force_partial_derivative_P_bead();





		/* ---------------------------------------------------
		 *
		 *
		 *  				Printing information
		 *
		 *
		 * ----------------------------------------------------
		 */
		void print_force_vector();

		void print_jacobian_matrix();


		ProblemParameters * _problemdata;
		DesignParameters * _designdata;





	private:


		void get_material_properties();


		/* ---------------------------------------------------
		 *
		 *
		 *  			Function pointers
		 *
		 *
		 * ----------------------------------------------------
		 */
		// Function pointer to calculate the force coefficients.
		// Assigned when resetting the beads
		void (ConstitutiveLaw::*_force_coefficients) ();

		void (ConstitutiveLaw::*_partial_derivative_alphaMax) ();

		void (ConstitutiveLaw::*_partial_derivative_P_alphaP_Initial) (PetscScalar & alphaMax_Initial);


		/* ---------------------------------------------------
		 *
		 *
		 *  Force and derivative force coefficients
		 *
		 *
		 * ----------------------------------------------------
		 */

		void calculate_force_dummy_plastic_coefficients();

		void calculate_force_linear_coefficients();

		void calculate_force_plastic_coefficients();

		/* ---------------------------------------------------
		 *
		 *
		 *  Jacobian and vector deformation functions and their pointers
		 *
		 *
		 * ----------------------------------------------------
		 */

		void (ConstitutiveLaw::*_jacobian_function) ();

		void jacobian_matrix_small_def();

		void jacobian_matrix_large_def();

		void (ConstitutiveLaw::*_get_distance_AND_unit_vector) ();

		void vector_beads_small_def();

		void vector_beads_large_def();

		void (ConstitutiveLaw::*_calculate_delta) ();

		void calculate_delta_large_def();

		void calculate_delta_small_def();

		/* ---------------------------------------------------
		 *
		 *
		 *  Partial derivatives with respect to the parameter
		 *
		 *
		 * ----------------------------------------------------
		 */


		void _partial_derivative_P_alphaP_Initial_plastic(PetscScalar & alphaMax_Initial);

		void _partial_derivative_P_alphaP_Initial_dummy_plastic(PetscScalar & alphaMax_Initial);


		/* ---------------------------------------------------
		 *
		 *
		 *  Partial derivatives with respect to alphaMax
		 *
		 *
		 * ----------------------------------------------------
		 */
		arma::colvec residual_derivative_alphaMax_vec;

		void partial_derivative_alphaMax_force_vector_plastic();

		void partial_derivative_alphaMax_force_vector_linear();

		void partial_derivative_alphaMax_force_vector_dummy_plastic();

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

		void (ConstitutiveLaw::*_state_eq_partial_derivative_alphaMax) ();

		void (ConstitutiveLaw::*_state_eq_partial_derivative_U) ();

		void state_eq_partial_derivative_alphaMax_plastic();

		void state_eq_partial_derivative_alphaMax_linear();

		void state_eq_partial_derivative_alphaMax_dummy_plastic();

		void state_eq_partial_derivative_U_plastic();

		void state_eq_partial_derivative_U_linear();

		void state_eq_partial_derivative_U_dummy_plastic();



		PetscReal dh_dalphaMax;
		arma::colvec dh_dU;




		BeadsData * _bead_to, * _bead_from;


		PetscReal * _force_vector, * _unit_vector, *_partial_derivative_U_force_vector , * _internal_partial_derivative_P_bead;

		PetscReal _partial_derivative_alphaMax_coefficient, _partial_derivative_P;

		PetscReal _distance, _RStar, _EStar, _Force,_Derivative_Force;

		PetscReal _delta;

		arma::mat  _jacobian, A;

		arma::colvec unit_vector_attach;

		bool _is_matprop_calculated, _is_distance_calculated, _is_delta_calculated;

		bool _is_force_coefficients_calculated;

		bool _is_alphaP_initial_calculated;

		// Parameters for the plastic evaluation
		bool _plastic;
		PetscReal _alphaY, _AreaY;

		// Plastic state variable. It needs to be set before evaluation
		// _derivative_alphaP_Initial is derivative w.r.t. initial alphaMax
		PetscReal _alphaMax, _alphaP, _alphaP_Initial, _derivative_alphaP_Initial;
		bool _is_alphaMax;

		// Constants for the plastic evaluation
		// 	alphaP = c1*delta - c2*alphaY + c3*alphaY*exp(-c4*(delta/alphaY-1.0));
		PetscReal c1, c2, c3, c4;
		PetscReal d1, d2, d3;
		PetscReal e1;
		PetscReal f1, f2;
		PetscReal g1;
		PetscReal h1;


		/*
		 * Coefficients for the dummy plastic law
		 */
		PetscReal k1, k2, tol;




};



#endif /* CONSTITUTUVELAW_H_ */
