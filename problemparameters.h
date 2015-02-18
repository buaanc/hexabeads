/*
 * problemparameters.h
 *
 *  Created on: Sep 26, 2014
 *      Author: miguel
 */
#include <vector>
#include "petscsys.h"
#include <utility>
#include <math.h>
#include <iostream>
#ifndef PROBLEMPARAMETERS_H_
#define PROBLEMPARAMETERS_H_

#define PI 3.14159265358979323846


class ProblemParameters
{
	public:
		ProblemParameters(){
			connectivity = NULL;
			N_Constraints = 0;
			N_Nodes = 0;
			N_Elements = 0;
			N_InitialConditions = 0;
			N_DesignVariables = 0;
			peak_load = 22;

			striker_DOF = 0;
			finiteDifference = false;
			matlabhistory = 0;
			loading = 1;
			timestep = 0.1;
			transient = 1;
			loadtime = 10;
			striker_ID = 1;
			analysis = 2;

			large_deformations = 0;
			TotalTime = 1;
			striker = true;

			sigmayielding = 0.55;
			InitialVelocity = 0.1;


			printstephistory = true;



		}
		~ProblemParameters() {
			if (connectivity != NULL)
				delete [] connectivity;
		}
		unsigned int N_Nodes, N_Elements, N_Constraints, N_InitialConditions, N_DesignVariables;

		std::vector<std::pair<unsigned int, unsigned int> > constraints_list;


		PetscScalar TotalTime;

		double timestep, sigmayielding;

		int * connectivity;


		bool finiteDifference;

		PetscInt analysis;

		bool transient;

		bool large_deformations;

		PetscInt loading;

		PetscInt matlabhistory;

		PetscInt printdisplacements;

		bool printstephistory;

		unsigned int striker_ID, striker_DOF;
		bool striker;

		PetscReal loadtime;
		PetscReal peak_load;


		double InitialVelocity;


		void print_connectivity(){
			if (connectivity != NULL)
				for (unsigned int i = 0; i< N_Elements; i++)
					std::cout<<"Element = "<<i<<" : "<<connectivity[2*i]<<" to "<<connectivity[2*i+1]<<std::endl;
		}

		void get_external_force(PetscReal & t, PetscReal & F_value){

			if (t < loadtime){
				F_value = peak_load/2.0 + (peak_load/2.0)*sin(PI/(loadtime/2.0)*t - PI/2.0);
				//F_value = peak_load*sin(PI/loadtime * t);
			}
			else
				F_value = 0.0;

//			loadtime = 40;
//			if (t < loadtime){
//				//F_value = peak_load/2.0 + (peak_load/2.0)*sin(PI/(loadtime/2.0)*t - PI/2.0);
//				F_value = 10*t;
//			}
//			else if (t >= loadtime && t < 2*loadtime)
//				F_value = 800 - 10*t;
//			else
//				F_value = 0;
		}

		void problem_info(){
			std::cout<<"##### PRINTING PROBLEM INFO ######"<<std::endl;
			std::cout<<"Number of elements = "<<N_Elements<<std::endl;
			std::cout<<"Number of nodes = "<<N_Nodes<<std::endl;
			std::cout<<"Number of constraints = "<<N_Constraints<<std::endl;
			std::cout<<"Number of initial conditions = "<<N_InitialConditions<<std::endl;
			std::cout<<"Total time = "<<TotalTime<<std::endl;
			std::cout<<"Striker ID = "<<striker_ID<<std::endl;
			std::cout<<"striker_DOF = "<<striker_DOF<<std::endl;
			std::cout<<"Transient ? "<<transient<<std::endl;
		}
};



#endif /* PROBLEMPARAMETERS_H_ */
