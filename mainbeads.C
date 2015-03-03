/*
 * mainbeads.C
 *
 *  Created on: Oct 7, 2014
 *      Author: miguel
 */
static char help[] = "This example demonstrates the use of DMNetwork interface for solving a nonlinear electric power grid problem.\n\
                      The available solver options are in the pfoptions file and the data files are in the datafiles directory.\n\
                      The data file format used is from the MatPower package (http://www.pserc.cornell.edu//matpower/).\n\
                      Run this program: mpiexec -n <n> ./PF\n\
                      mpiexec -n <n> ./PF -pfdata <filename>\n";


#include "system.h"
#include "quadbeads.h"
#include "IpIpoptApplication.hpp"
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char ** argv)
{
	PetscInitialize(&argc,&argv,"beadsoptions",help);
	PetscErrorCode       ierr;

	PetscInt sensitivities;
	ierr =PetscOptionsInt("-sensitivities", "1 for running just the sensitivities and the FD comparison, 0 for running the optimization","", 0, &sensitivities, NULL);

	if (sensitivities){
		System QuadBeadsDummy;



		/*
		 * Initialize and read data
		 */
		ierr = QuadBeadsDummy.InitializeData();

		/*
		 * Set up the DM Network
		 */
		ierr = QuadBeadsDummy.SetUpDMNetwork();

		/*
		 * Create PETSc Vectors and matrices
		 */
		ierr = QuadBeadsDummy.CreatePETScObjects();


		QuadBeadsDummy.problemData.printstephistory = true;

		/*
		 * Set up solver
		 */
		ierr = QuadBeadsDummy.SetUpSolver();

		QuadBeadsDummy.ProcessOptimalSolution();


		/*
		 * Solve
		 */
		ierr = QuadBeadsDummy.Solve();


		/*
		 * Calculate Obj Function
		 */
		ierr = QuadBeadsDummy.ProcessObjFunctionTime();

		/*
		 * Fix the value of FunctionValueGlobalOriginal to use it in the FD;
		 */
		QuadBeadsDummy.FunctionValueGlobalOriginal = QuadBeadsDummy.FunctionValueGlobal;

		/*
		 * Adjoint Analysis
		 */
		std::cout<<"Adjoint Analysis"<<std::endl;
		ierr = QuadBeadsDummy.AdjointSolve();

		ierr = QuadBeadsDummy.ProcessGradient();

		/*
		 * Check sensitivities
		 */
		ierr = QuadBeadsDummy.FiniteDifference();

		ierr = QuadBeadsDummy.CheckGradient();

	}
	else{

//		SmartPtr<TNLP> OptProblem = new QuadBeads;
//
//
//		// Create a new instance of IpoptApplication
//		//  (use a SmartPtr, not raw)
//		// We are using the factory, since this allows us to compile this
//		// example with an Ipopt Windows DLL
//		SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
//		// Change some options
//		// Note: The following choices are only examples, they might not be
//		//       suitable for your optimization problem.
//		app->Options()->SetNumericValue("tol", 1e-6);
//		app->Options()->SetStringValue("mu_strategy", "adaptive");
//		app->Options()->SetStringValue("hessian_approximation", "limited-memory");
//		// The following overwrites the default name (ipopt.opt) of the
//		// options file
//		// app->Options()->SetStringValue("option_file_name", "hs071.opt");
//
//		// Intialize the IpoptApplication and process the options
//		ApplicationReturnStatus status;
//		status = app->Initialize();
//		if (status != Solve_Succeeded) {
//		std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
//		return (int) status;
//		}
//
//		// Ask Ipopt to solve the problem
//		status = app->OptimizeTNLP(OptProblem);
//
//		if (status == Solve_Succeeded) {
//		std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
//		}
//		else {
//		std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
//		}

	}




	return 0;

}
