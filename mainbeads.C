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
#include "MMA.h"
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char ** argv)
{
	PetscInitialize(&argc,&argv,"beadsoptions",help);
	PetscErrorCode       ierr;

	PetscInt sensitivities;
	ierr =PetscOptionsInt("-sensitivities", "1 for running just the sensitivities and the FD comparison, 0 for running the optimization","", 0, &sensitivities, NULL);

	if (sensitivities){
		System hexabeadsDummy;



		/*
		 * Initialize and read data
		 */
		ierr = hexabeadsDummy.InitializeData();

		/*
		 * Set up the DM Network
		 */
		ierr = hexabeadsDummy.SetUpDMNetwork();

		/*
		 * Create PETSc Vectors and matrices
		 */
		ierr = hexabeadsDummy.CreatePETScObjects();


		hexabeadsDummy.problemData.printstephistory = true;

		/*
		 * Set up solver
		 */
		ierr = hexabeadsDummy.SetUpSolver();




		/*
		 * Solve
		 */
		ierr = hexabeadsDummy.Solve();


		/*
		 * Calculate Obj Function
		 */
		ierr = hexabeadsDummy.ProcessObjFunctionTime();

		/*
		 * Fix the value of FunctionValueGlobalOriginal to use it in the FD;
		 */
		hexabeadsDummy.FunctionValueGlobalOriginal = hexabeadsDummy.FunctionValueGlobal;

		/*
		 * Adjoint Analysis
		 */
		std::cout<<"Adjoint Analysis"<<std::endl;
		ierr = hexabeadsDummy.AdjointSolve();

		ierr = hexabeadsDummy.ProcessGradient();

		/*
		 * Check sensitivities
		 */
		ierr = hexabeadsDummy.FiniteDifference();

		ierr = hexabeadsDummy.CheckGradient();

		hexabeadsDummy.ProcessOptimalSolution();

	}
	else{

		System hexabeadsDummy;

		/*
		 * Initialize and read data
		 */
		ierr = hexabeadsDummy.InitializeData();

		/*
		 * Set up the DM Network
		 */
		ierr = hexabeadsDummy.SetUpDMNetwork();

		/*
		 * Create PETSc Vectors and matrices
		 */
		ierr = hexabeadsDummy.CreatePETScObjects();

		/*
		 * Set up solver
		 */
		ierr = hexabeadsDummy.SetUpSolver();

		MMA *mma;
		PetscInt itr=0;

		Vec design_variables = hexabeadsDummy.get_design_variables();
		Vec Xold = hexabeadsDummy.get_xold();
		Vec gradient = hexabeadsDummy.get_gradient();

		PetscScalar * ObjFunction = hexabeadsDummy.get_ObjFunction();
		PetscInt n_design_variables = hexabeadsDummy.get_number_of_design_variables();
		mma = new MMA(n_design_variables,1.0,design_variables);

		Vec Xmin = hexabeadsDummy.get_xmin_vec();
		Vec Xmax = hexabeadsDummy.get_xmax_vec();

		Vec * constraint_gradient = hexabeadsDummy.get_constraint_gradient();


		PetscInt maxItr = 400;

		PetscReal movlim = 0.2;

		// STEP 7: OPTIMIZATION LOOP
		  PetscScalar ch = 1.0;
		  double t1,t2;
		  while (itr < maxItr && ch > 1e-6){
			// Update iteration counter
			itr++;

			// start timer
			t1 = MPI_Wtime();

			ierr = hexabeadsDummy.OptimizerRoutine();

			PetscScalar * constraint = hexabeadsDummy.get_constraint();

			// Sets outer movelimits on design variables
			ierr = mma->SetOuterMovelimit(hexabeadsDummy.get_xmin(),hexabeadsDummy.get_xmax(),movlim,design_variables,Xmin,Xmax); CHKERRQ(ierr);

			// Update design by MMA
			ierr = mma->Update(design_variables,gradient,constraint,constraint_gradient,Xmin,Xmax); CHKERRQ(ierr);

			// Update the vector alphaMax with the new design
			hexabeadsDummy.InitialDesignSolution();


			// Inf norm on the design change
			ch = mma->DesignChange(design_variables,Xold);


			// stop timer
			t2 = MPI_Wtime();

			// Print to screen
			PetscPrintf(PETSC_COMM_WORLD,"It.: %i, obj.: %f, g[0]: %f, ch.: %f, time: %f\n",
						itr,*ObjFunction,constraint[0], ch,t2-t1);



		  }


			hexabeadsDummy.ProcessOptimalSolution();


			/*
			 * Print out the final design to MATLAB
			 */
			// File name
			PetscViewer    viewer;
			std::string filehistory("PostProcessing/OptimalDesign.mat");
			// Open and set Binary file
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,filehistory.c_str(),FILE_MODE_WRITE,&viewer);
			PetscViewerSetFormat(viewer, PETSC_VIEWER_BINARY_MATLAB);
			// Print
			VecView(design_variables,viewer);
			PetscViewerDestroy(&viewer);

			delete mma;

	}




	return 0;

}
