/*
 * beadsdata.h
 *
 *  Created on: Sep 26, 2014
 *      Author: miguel
 */

#ifndef BEADSDATA_H_
#define BEADSDATA_H_

#include <vector>
#include <iostream>
#include <petscts.h>
class BeadsData
{
	public:

		BeadsData(){
			constrained_dof_x = false;
			constrained_dof_y = false;
		}
		// Material properties and coordinates
		double E, R, rho, nu, coordx, coordy, coordz, density_design;

		// type = 1, regular bead; type = 2, wall (completely constrained); type = 3, striker
		unsigned int type, node_id;

		bool constrained_dof_x, constrained_dof_y;
		double constraint_value_x, constraint_value_y;

		// Displacement values
		PetscReal U_i, U_j;
		void print_bead_information(){
			std::cout<<"Bead # "<<node_id<<std::endl;
			std::cout<<"E = "<<E<<std::endl;
			std::cout<<"R = "<<R<<std::endl;
			std::cout<<"rho = "<<rho<<std::endl;
			std::cout<<"nu = "<<nu<<std::endl;
			std::cout<<"coordx = "<<coordx<<std::endl;
			std::cout<<"coordy = "<<coordy<<std::endl;
			std::cout<<"type = "<<type<<std::endl;
			std::cout<<"constrained_dof_x = "<<constrained_dof_x<<std::endl;
			std::cout<<"constrained_dof_y = "<<constrained_dof_y<<std::endl;
		}

		void get_mass(double & mass){
			mass = rho * 4.0/3.0*PETSC_PI*PetscPowRealInt(R,3);
		}

		void get_mass_derivative(double & mass_derivative){

			mass_derivative = rho * 4.0/3.0*PETSC_PI*PetscPowRealInt(R,3);
		}

};

typedef class BeadsData * arrayBeads;

#endif /* BEADSDATA_H_ */
