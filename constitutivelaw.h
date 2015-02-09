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
class ConstitutiveLaw
{
	public:

		ConstitutiveLaw();
		ConstitutiveLaw(ProblemParameters * problemdata, DesignParameters * designdata);
		~ConstitutiveLaw();

		double * force_vector();

		double & force_coefficient();

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
		void set_alphaMax(const double & alphaMax);

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

		/* ---------------------------------------------------
		 *
		 *
		 *  	Getters
		 *
		 *
		 * ----------------------------------------------------
		 */

		double & get_alphaMax(){ return _alphaMax;}

		double & get_delta(){ return _delta;}




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

		void state_eq_vector_product_U(const double & gamma, arma::colvec & result);

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


		double * partial_derivative_U_force_vector();

		double * partial_derivative_P();

		double & partial_derivative_alphaMax_force_vector();

		arma::colvec & residual_derivative_alphaMax();


		double * internal_force_partial_derivative_P_bead(unsigned int & bead);





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

		void (ConstitutiveLaw::*_partial_derivative_P_function) ();

		void (ConstitutiveLaw::*_partial_derivative_alphaMax) ();


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

		void partial_derivative_P_plastic();

		void partial_derivative_P_linear();

		void partial_derivative_P_dummy_plastic();


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



		double dh_dalphaMax;
		arma::colvec dh_dU;




		BeadsData * _bead_to, * _bead_from;


		double * _force_vector, * _unit_vector, *_partial_derivative_U_force_vector, * _partial_derivative_P , * _internal_partial_derivative_P_bead;

		double _partial_derivative_alphaMax_coefficient;

		double _distance, _delta, _RStar, _EStar, _Force,_Derivative_Force;

		arma::mat  _jacobian, A;

		arma::colvec unit_vector_attach;

		bool _is_matprop_calculated, _is_distance_calculated, _is_delta_calculated;

		bool _is_force_coefficients_calculated;

		// Parameters for the plastic evaluation
		bool _plastic;
		double _alphaY, _AreaY;

		// Plastic state variable. It needs to be set before evaluation
		double _alphaMax, _alphaP;
		bool _is_alphaMax;

		// Constants for the plastic evaluation
		// 	alphaP = c1*delta - c2*alphaY + c3*alphaY*exp(-c4*(delta/alphaY-1.0));
		double c1, c2, c3, c4;
		double d1, d2, d3;
		double e1;
		double f1, f2;
		double g1;
		double h1;


		/*
		 * Coefficients for the dummy plastic law
		 */
		PetscScalar k1, k2;




};



#endif /* CONSTITUTUVELAW_H_ */
