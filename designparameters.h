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

		unsigned int  N_DesignVariables;

		PetscInt N_ObjFunc;
		// MaxOrMin type: 1 for maximization, 0 for minimization
		std::vector<unsigned int> MaxOrMin;

		std::vector<std::vector<int> > TargetArea_list;

		unsigned int SpaceNorm, TimeNorm;

		double InitialGuess;

		PetscInt objfunctionType;

		PetscReal initial_design;

		PetscReal xmax, xmin;

		// Print forces in obj function
		bool printforces;


		double Obj_reltol, KKT_tol, DesVar_reltol;
		void design_info(){
			std::cout<<"##### PRINTING Design INFO ######"<<std::endl;
			std::cout<<"Number of objective functions = "<<N_ObjFunc<<std::endl;
			for (PetscInt i = 0;i<N_ObjFunc; i++){
				std::cout<<"Number of target elements for Obj Funct # "<<i<<" = "<<TargetArea_list[i].size()<<std::endl;
				for (unsigned int j = 0;j<TargetArea_list[i].size(); j++){
					std::cout<<TargetArea_list[i][j]<<" ";
				}
				std::cout<<std::endl;
			}
			std::cout<<"Number of design variables = "<<N_DesignVariables<<std::endl;
		}
};


#endif /* DESIGNPARAMETERS_H_ */
