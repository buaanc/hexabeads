/*
 * system.h
 *
 *  Created on: Oct 7, 2014
 *      Author: miguel
 */

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "problemparameters.h"
#include "designparameters.h"
#include "beadsdata.h"
#include "constitutivelaw.h"
#include <petscdmnetwork.h>

#include <iostream>
#include <fstream>
#include <utility>
#include <iomanip>
#include <sstream>
#include "vechistory.h"
#include "vechistoryRK.h"

#include <mpi.h>

#include <time.h>

	/*
	 *
	 * Equilibrium equations
	 *
	 * M \ddot{u} + N(u) = F^{ext}
	 *
	 * Initial conditions, we might need to modify them
	 * u 	 	= 0 	everywhere
	 * \dot{u} 	= 0		everywhere
	 *
	 * Boundary conditions
	 * u 	 	= 0 	wall beads
	 *
	 * We transform it into a first order differential equation because the
	 * matrix M is diagonal
	 *
	 * \dot{v_i} = 1.0 / m_i * (F^ext - N(u))
	 * \dot{u_i} = v_i
	 *
	 * u 	 	= 0 	for	wall beads
	 * v 		= 0		for	wall beads
	 *
	 * Boundary conditions
	 * u 	 	= 0 	wall beads
	 * v 	 	= 0 	wall beads
	 *
	 */

struct _p_EDGEDATA{
  PetscInt      VariableGroup;
};

typedef struct _p_EDGEDATA *EDGEDATA;

class System
{
	public:

	System();
	/*
	 * Initialize data
	 */
	PetscErrorCode InitializeData();
	/*
	 * Reading data
	 */
	PetscErrorCode ReadBeadsData(arrayBeads & arrayBeads);

	/*
	 * Set up network
	 */
	PetscErrorCode SetUpDMNetwork();

	/*
	 * Create PETSc Vectors and matrices
	 */
	PetscErrorCode CreatePETScObjects();

	/*
	 * Set up solver
	 */
	PetscErrorCode SetUpSolver();
	/*
	 * Set up the solver for the FD so we can
	 * use the same time steps than the original
	 * solve
	 */
	PetscErrorCode FDSetUpSolver();

	/*
	 * Solve
	 */
	PetscErrorCode Solve();

	/*
	 * Adjoint Solve, it has to be called after: Solve(), CalculateObjFunction() and PartialDerivatives()
	 */
	PetscErrorCode AdjointSolve();

	/*
	 * Calculating the RHS parameter derivative for each design variable
	 * - DF is a vector corresponding to the design variable we're in
	 */
	PetscErrorCode RHSParameterDerivatives(PetscReal & t, Vec & U, Vec & DF,  std::vector<PetscInt> & DesignElements, Vec & alphaMax);

	PetscErrorCode RHSParameterDerivatives_FD(PetscReal & t, Vec & U, Vec & DF, PetscInt Design_Index, Vec & alphaMax);

	/*
	 * Destructor
	 */
	~System();

//	PetscErrorCode FormRHSFunction(TS ts,PetscReal t, Vec U, Vec F,void *appctx);
//
//
//	/*
//	 * Function for the monitor. Routine that is called once per time step
//	 */
//	PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal ptime,Vec U,void *ctx);
	/*
	 * Initial solution
	 */
	PetscErrorCode InitialSolution(DM networkdm,Vec X);

	PetscErrorCode InitialDesignSolution();


	/*
	 * Objective function calculation
	 */

	PetscErrorCode ObjFunctionTime(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode (System::*CalculateObjFunctionTime) (Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode CalculateObjFunctionTime_P_Norm(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode CalculateObjFunctionTime_Simple(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);

	PetscErrorCode ProcessObjFunctionTime();


	PetscErrorCode ProcessGradient();

	/*
	 * Process optimal solution
	 */

	void ProcessOptimalSolution(){
		problemData.printdisplacements = 1;
		designData.printforces = true;
		problemData.printstephistory = true;
	}

	PetscErrorCode PrintForcesTargetArea(Vec & U, Vec & alphaMax, PetscInt & ptime);



	/*
	 * Calculate partial derivatives and print them
	 */
	PetscErrorCode Partial_Derivatives();
	PetscErrorCode (System::*PartialDerivatives) ();
	PetscErrorCode PartialDerivatives_P_Norm();
	PetscErrorCode PartialDerivatives_Simple();

	PetscErrorCode EvaluatePartialImplicitDerivatives(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode EvaluatePartialExplicitDerivatives(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode (System::*PartialImplicitDerivatives) (Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode (System::*PartialExplicitDerivatives) (Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode EvaluatePartialImplicitDerivatives_P_Norm(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode EvaluatePartialExplicitDerivatives_P_Norm(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode EvaluatePartialImplicitDerivatives_Simple(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);
	PetscErrorCode EvaluatePartialExplicitDerivatives_Simple(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);

	//PetscErrorCode PrintPartialDerivatives();

	/*
	 * Finite difference check
	 */
	PetscErrorCode FD_parameter(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);

	PetscErrorCode FD_variables(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);

	PetscErrorCode FD_statevariables(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct);

	PetscErrorCode RestartSolver();

	PetscErrorCode FD_SRV();


	/*
	 * Global Finite Differente
	 */
	PetscErrorCode FiniteDifference();

	PetscErrorCode CheckGradient();
	bool _is_fd_calculated;
	bool _is_gradient_calculated;

	/*
	 * Public members
	 */
	PetscMPIInt rank;

	arrayBeads 		     allbeadsData;

	EDGEDATA			edgeData;

	PetscInt             numEdges, numVertices;

	ProblemParameters 	problemData;
	DesignParameters 	designData;
	ConstitutiveLaw 	elementFunction;
	std::ofstream 		steps_history;


	PetscInt             componentkey[4];

	/*
	 * PETSc Objects
	 */
	Vec X,F, alphaMax, Felement;
	Mat J;
	DM                   networkdm;

	/*
	 * TS object
	 */
	TS ts;

	/*
	 * PETSc Adjoint Objects
	 * Lambda is the adjoint variable associated with the implicit variables
	 * Gamma is the adjoint variable associated with the state variables
	 * Phi is the adjoint variable associated with the RK stage evaluations of the implicit variables
	 * K_RK is the variable associated with the RK stage function evaluations of the implicit variables
	 * U_RK is the variable associated with the RK stage of the implicit variables (Capital Y in the notes)
	 * Omega is the vector to which will apply the linear operator to obtain Phi for each stage
	 * Omega = b_i * lambda_{n-1} + sum^i_j a_ji * Phi_j
	 *
	 * JACP holds the RHS derivatives for each design variable.
	 *
	 * Mu is the adjoint variable for the parameters
	 * V_i is the RK stage variable for the variable Mu
	 *
	 * W_i is the RK stage variable for the variable Gamma
	 *
	 * DN is an auxiliary vector to calculate the RHS derivative w.r.t the para meter
	 */
	Vec Lambda, Gamma, Omega, * Phi, * K_RK, * U_RK, *JACP,  DN, Mu, * V_i, * W_i;

	Vec TMP, TMPalpha;

	/*
	 * Local number of parameters in processor
	 */
	PetscInt local_N_parameters, * parameters_indices;
	/*
	 * Bogacki–Shampine method parameters
	 */
	PetscScalar *c, *b, *a, *ones;

	/*
	 * Auxiliar variable for the VecMDot function
	 */
	PetscScalar * Mdot_product_result;

	/*
	 * Vector where we keep the coefficients for the derivatives
	 * of each objective function
	 */
	std::vector<PetscReal> CapitalOmega, FunctionObjetivoPartial;
	/*
	 * Implicit variables history
	 */
	VecHistory Uhistory, alphaMaxhistory;
	VecHistoryRungeKutta U_RungeKuttahistory;
	bool SaveStages;
	/*
	 * Time history
	 */
	std::vector<PetscReal> timeHistory;
	std::vector<PetscReal> timeStepHistory;

	std::vector<PetscReal> matlab_timestephistory;
	std::vector<PetscReal>::iterator matlab_timestephistory_iterator;
	PetscInt total_timesteps_matlab;

	bool solve_in_fd;
	std::vector<PetscReal>::iterator timeStepHistory_iterator , timeStepHistory_iterator_end;
	PetscScalar previous_time;
	/*
	 * bool to check that partial derivatives have been calculated, otherwise
	 * it doesn't make sense to do the FD check
	 */
	bool _is_partial_derivatives_calculated;
	/*
	 * Partial derivatives w.r.t U
	 */
	Vec dFdU, dFdU_partial;
	/*
	 * Partial derivatives w.r.t alphaMax
	 */
	Vec dFdalphaMax, dFdalphaMax_partial;

	/*
	 * Partial derivatives w.r.t alphaMax
	 */
	Vec dFdP, dFdP_partial;

	/*
	 * Function Gradient
	 */
	std::vector<PetscScalar> Gradient;
	std::vector<PetscScalar> Gradient_FD;

	/*
	 * Objective function value and derivative coefficients
	 */
	PetscReal FunctionValue, FunctionValueBeforePnorm, FunctionValueGlobal, FunctionValueGlobalBeforePenalty;\

	/*
	 * Use this FunctionValueGlobalOriginal for the FD;
	 */
	PetscReal FunctionValueGlobalOriginal;
	PetscReal LocalVolConstraintValue, VolConstraintValue, LocalSRVConstraintValue, SRVConstraintValue;
	PetscReal LocalPenaltyTermValue, PenaltyTermValue;
	PetscReal theta;
	std::vector<PetscReal> beta;
	std::vector<std::vector<PetscReal> > gamma;


	arma::colvec VolGradient;
	Vec SRVGradient, PenaltyGradient;

	/*
	 * Container with the element variables. It is a container
	 * of containers, one for each variable, which contains
	 * all the elements that depend on that variable
	 *
	 * Example: 3 variables, 7 elements
	 *
	 * 	Variable 1:		Elements: 2,3,4
	 * 	Variable 2:		Elements: 1,5
	 * 	Variable 1:		Elements: 6,7
	 */
	std::vector<std::vector<PetscInt> > design_elements;

	/*
	 * Design variables
	 */
	std::vector<PetscScalar> alphaMaxInit;


	/*
	 * Number of stages for the RK Algorithm
	 */
	unsigned int Stages;


	/*
	 * PETSc Stages
	 */
	PetscLogStage solvestage, adjointstage;

	/*
	 * Set up stages
	 */
	PetscErrorCode SetUpStages();


	clock_t start, end;
	double  totalcputime;




	private:

		PetscErrorCode PerturbPETScVector(Vec U, PetscInt & Index, PetscReal & value);

		void trim1(std::string& str)
		{
			std::string::size_type pos=str.find_last_not_of(' ');
			if (pos != std::string::npos) {
				str.erase(pos+1);
				pos = str.find_first_not_of(' ');
				if (pos != std::string::npos) str.erase(0,pos);
			}
			else str.erase(str.begin(),str.end());
		}


};



#endif /* SYSTEM_H_ */
