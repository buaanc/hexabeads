/*
 * designparameters.h
 *
 *  Created on: Sep 26, 2014
 *      Author: miguel
 */

#ifndef DESIGNPARAMETERS_H_
#define DESIGNPARAMETERS_H_


class DesignParameters
{
	public:
		unsigned int SIMP_parameter;

		double volfrac_min;

		double penalty_weight;

		unsigned int N_ObjFunc, N_ConstrFunc, N_DesignVariables;

		// MaxOrMin type: 1 for maximization, 0 for minimization
		std::vector<unsigned int> MaxOrMin;

		std::vector<std::vector<int> > TargetArea_list;

		// Constraint type: 1 for volume, 2 for SRV
		std::vector<unsigned int> ConstraintsType;

		// Constraint bund, number of intruders allow
		std::vector<double> ConstraintsBounds;

		// Equality constraint tolerance, 0 por inequality
		std::vector<double> ConstraintsTolerance;

		unsigned int SpaceNorm, TimeNorm;

		double InitialGuess;

		PetscScalar alphaMaxInitial;

		PetscInt objfunctionType;

		PetscReal initial_design;


		double Obj_reltol, KKT_tol, DesVar_reltol;
		void design_info(){
			std::cout<<"##### PRINTING Design INFO ######"<<std::endl;
			std::cout<<"Number of objective functions = "<<N_ObjFunc<<std::endl;
			for (unsigned int i = 0;i<N_ObjFunc; i++){
				std::cout<<"Number of target elements for Obj Funct # "<<i<<" = "<<TargetArea_list[i].size()<<std::endl;
				for (unsigned int j = 0;j<TargetArea_list[i].size(); j++){
					std::cout<<TargetArea_list[i][j]<<" ";
				}
				std::cout<<std::endl;
			}
			std::cout<<"Number of constraints = "<<N_ConstrFunc<<std::endl;
			std::cout<<"Number of design variables = "<<N_DesignVariables<<std::endl;
		}
};


#endif /* DESIGNPARAMETERS_H_ */
