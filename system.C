/*
 * system.C
 *
 *  Created on: Oct 7, 2014
 *      Author: miguel
 */
#include "system.h"
#include <dirent.h>
#include <cassert>
#include <petsctime.h>

void SetRKDefault(System & sistema)
{
	sistema.c = new PetscScalar[sistema.Stages];
	sistema.b = new PetscScalar[sistema.Stages];
	sistema.a = new PetscScalar[sistema.Stages*sistema.Stages];
	sistema.ones = new PetscScalar[sistema.Stages];

	sistema.c[0] = 0.0;
	sistema.c[1] = 1.0/2.0;
	sistema.c[2] = 3.0/4.0;
	sistema.c[3] = 1.0;

	sistema.b[0] = 2.0/9.0;
	sistema.b[1] = 1.0/3.0;
	sistema.b[2] = 4.0/9.0;
	sistema.b[3] = 0.0;

	sistema.a[0] = 0.0;
	sistema.a[1] = 0.0;
	sistema.a[2] = 0.0;
	sistema.a[3] = 0.0;
	sistema.a[4] = 1.0/2.0;
	sistema.a[5] = 0.0;
	sistema.a[6] = 0.0;
	sistema.a[7] = 0.0;
	sistema.a[8] = 0.0;
	sistema.a[9] = 3.0/4.0;
	sistema.a[10] = 0.0;
	sistema.a[11] = 0.0;
	sistema.a[12] = 2.0/9.0;
	sistema.a[13] = 1.0/3.0;
	sistema.a[14] = 4.0/9.0;
	sistema.a[15] = 0.0;

	sistema.ones[0] = 1.0;
	sistema.ones[1] = 1.0;
	sistema.ones[2] = 1.0;
	sistema.ones[3] = 1.0;

}


PetscErrorCode FormRHSFunction(TS ts,PetscReal t, Vec U, Vec F,void *appctx)
{
  PetscErrorCode ierr;
  DM             networkdm;
  Vec           localU,localF;
  PetscInt      e;
  PetscInt      eStart,eEnd,vfrom,vto;
  const PetscScalar *uarr;
  PetscScalar U_i_from, U_j_from, U_i_to, U_j_to;
  PetscScalar   *farr, *alphaMaxarr;
  PetscInt      offsetfrom,offsetto;
  DMNetworkComponentGenericDataType *arr;

  arrayBeads beadsdataFrom, beadsdataTo ;

  System * sistema = (System*) appctx;

  PetscFunctionBegin;
  ierr = TSGetDM(ts,&networkdm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localF);CHKERRQ(ierr);
  ierr = VecSet(F,0.0);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);


  ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);

  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);



#ifdef DEBUG
  	  std::cout<<"             "<<std::endl;
  	  std::cout<<"  RHS EVAL  "<<std::endl;
	  std::cout<<"time = "<<t<<std::endl;

	  std::cout<<"alphaMax en RHS"<<std::endl;
	  VecView(sistema->alphaMax,PETSC_VIEWER_STDOUT_WORLD);


#endif


  // Get the components of the alphaMax Vec
  VecGetArray(sistema->alphaMax,&alphaMaxarr);
  // Index for alphaMaxarr

  /* Edge Data necessary to know which variable group it belongs to */
  EDGEDATA elementData;
  PetscInt key, offsetd_eStart;



  unsigned int k = 0;

  for (e=eStart; e < eEnd; e++) {
    PetscInt    offsetd,key;

	const PetscInt *cone;
	DMNetworkGetConnectedNodes(networkdm,e,&cone);
	vfrom = cone[0];
	vto   = cone[1];


	ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
	// Cast component data
	beadsdataFrom = (arrayBeads)(arr + offsetd);

	ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
	// Cast component data
	beadsdataTo = (arrayBeads)(arr + offsetd);

	ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
	ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

	// Grab displacement components
	U_i_from = uarr[offsetfrom + 2];
	U_j_from = uarr[offsetfrom + 3];

	U_i_to = uarr[offsetto + 2];
	U_j_to = uarr[offsetto + 3];

	// Reset the beads for the constitutive law object
	sistema->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, sistema->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);



	// Set initial design value, first, we need to know which variable group it belongs to
	DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd_eStart);
	elementData = (EDGEDATA)(arr + offsetd_eStart);
	sistema->elementFunction.set_alphaP_initial(sistema->alphaMaxInit[elementData->VariableGroup]);

#ifdef DEBUG
			std::cout<<"variable value = "<<std::setprecision(10)<<sistema->alphaMaxInit[elementData->VariableGroup]<<std::endl;

			double & alphaPInit = sistema->elementFunction.get_alphaP_initial();
			std::cout<<"alphaP Init = "<<std::setprecision(10)<<alphaPInit<<std::endl;


#endif

	// Set the value of alphaMax
	sistema->elementFunction.set_alphaMax(alphaMaxarr[k]);

	sistema->elementFunction.calculate_delta();
	PetscReal * force_vector = sistema->elementFunction.force_vector();



#ifdef DEBUG

	std::cout<<"Derivative Coef = "<<std::setprecision(10)<<sistema->elementFunction.get_derivative_force_coefficient()<<std::endl;
	//sistema->elementFunction.check_jacobian();
#endif

	k++;

#ifdef DEBUG
	std::cout<<"alphaMax RHS = "<<alphaMaxarr[k-1]<<" element = "<<k-1<<std::endl;
	std::cout<<"delta RHS = "<<sistema->elementFunction.get_delta()<<" element = "<<k-1<<std::endl;

	std::cout<<"force coefficient = "<<std::setprecision(20)<<sistema->elementFunction.force_coefficient()<<std::endl;
#endif

	farr[offsetfrom] -= force_vector[0];
	farr[offsetfrom + 1] -= force_vector[1];

	farr[offsetfrom + 2] = uarr[offsetfrom];
	farr[offsetfrom + 3] = uarr[offsetfrom + 1];

	farr[offsetto] -= force_vector[2];
	farr[offsetto + 1] -= force_vector[3];

	farr[offsetto + 2] = uarr[offsetto];
	farr[offsetto + 3] = uarr[offsetto + 1];

	if (beadsdataFrom->constrained_dof_x){
		farr[offsetfrom] = 0.0;
		farr[offsetfrom + 2] = 0.0;
	}
	if (beadsdataFrom->constrained_dof_y){
		farr[offsetfrom + 1] = 0.0;
		farr[offsetfrom + 3] = 0.0;
	}
	if (beadsdataTo->constrained_dof_x){
		farr[offsetto] = 0.0;
		farr[offsetto + 2] = 0.0;
	}
	if (beadsdataTo->constrained_dof_y){
		farr[offsetto + 1] = 0.0;
		farr[offsetto + 3] = 0.0;
	}

  }

   /*
    * Adding the value of the objective function
    */

   if (eStart == 0){ // Conditional so we only add the new DOF in the processor where the first element is
  	 ierr = DMNetworkGetVariableOffset(networkdm,eStart,&offsetfrom);CHKERRQ(ierr);

  	 for (PetscInt j = 0; j < sistema->designData.N_ObjFunc; j++){
  		 sistema->ObjFunctionTime(U,sistema->alphaMax,j);
  		 farr[offsetfrom + j] = sistema->FunctionValue;
  	 }

   }

#ifdef DEBUG
   	  std::cout<<"Obj Function RHS = "<<sistema->FunctionValue<<std::endl;
#endif

  ierr = VecRestoreArray(sistema->alphaMax,&alphaMaxarr);CHKERRQ(ierr);

  // We need to loop over the beads to multiply the farr vector (only the velocity components) with the inverse of the mass
  PetscInt vStart, vEnd, v;
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  for (v=vStart; v < vEnd; v++) {
    PetscInt    offsetd,key;
    PetscReal mass;
	ierr = DMNetworkGetComponentTypeOffset(networkdm,v,0,&key,&offsetd);CHKERRQ(ierr);
	// Cast component data
	beadsdataFrom = (arrayBeads)(arr + offsetd);

	ierr = DMNetworkGetVariableOffset(networkdm,v,&offsetfrom);CHKERRQ(ierr);

	beadsdataFrom->get_mass(mass);



	if (beadsdataFrom->type == 3){
		PetscReal Fext;
		sistema->problemData.get_external_force(t,Fext);
		farr[offsetfrom + sistema->problemData.striker_DOF - 1] += Fext;
	}


	farr[offsetfrom] *= 1.0/mass;
	farr[offsetfrom + 1] *= 1.0/mass;

#ifdef DEBUG
  std::cout<<"mass = "<<mass<<std::endl;
#endif
  }





  ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);


  ierr = DMLocalToGlobalBegin(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localF);CHKERRQ(ierr);


#ifdef DEBUG
  std::cout<<"F en RHS"<<std::endl;
  VecView(F,PETSC_VIEWER_STDOUT_WORLD);
#endif

  if (sistema->SaveStages)
	  sistema->U_RungeKuttahistory.insert(U);


  PetscFunctionReturn(0);

  return ierr;
}

PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal ptime,Vec U,void *ctx)
{
	  PetscErrorCode ierr;
	  DM             networkdm;
	  Vec           localU;
	  PetscInt      e;
	  PetscInt      eStart,eEnd,vfrom,vto;
	  const PetscScalar *uarr;
	  PetscScalar U_i_from, U_j_from, U_i_to, U_j_to;
	  PetscScalar    *alphaMaxarr;
	  PetscInt      offsetfrom,offsetto;
	  DMNetworkComponentGenericDataType *arr;

	  arrayBeads beadsdataFrom, beadsdataTo;

	  System * sistema = (System*) ctx;

	  PetscFunctionBegin;
	  ierr = TSGetDM(ts,&networkdm);CHKERRQ(ierr);
	  ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);

#ifdef DEBUG
//	  std::cout<<"primal time = "<<std::setprecision(20)<<ptime<<std::endl;
//	  std::cout<<"Solution"<<std::endl;
//	  VecView(U,PETSC_VIEWER_STDOUT_WORLD);
#endif


	  ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	  ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);


	  ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);


	  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
	  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	  // Get the components of the alphaMax Vec
	  VecGetArray(sistema->alphaMax,&alphaMaxarr);

	  // Allocate memory for the next RK stages
	  sistema->U_RungeKuttahistory.new_timestep(U);


#ifdef DEBUG
	  if (sistema->problemData.matlabhistory == 1){
		  ierr = TSSetTimeStep(ts,*sistema->matlab_timestephistory_iterator);
		  std::cout<<"time step = "<<*sistema->matlab_timestephistory_iterator<<std::endl;
		  std::cout<<"time = "<<ptime<<std::endl;
		  sistema->matlab_timestephistory_iterator++;
	  }
#endif





	  if (sistema->solve_in_fd && sistema->problemData.matlabhistory != 1 ){

		  if (sistema->timeStepHistory_iterator != sistema->timeStepHistory_iterator_end){
			  ierr = TSSetTimeStep(ts,*sistema->timeStepHistory_iterator);
#ifdef DEBUG
		  std::cout<<"time step que metimos en el FD = "<<*sistema->timeStepHistory_iterator<<std::endl;
		  std::cout<<"time = "<<ptime<<std::endl;
#endif
		  }
		  sistema->timeStepHistory_iterator++;
	  }

#ifdef DEBUG
		  std::cout<<"time = "<<ptime<<std::endl;
#endif

	  /* Edge Data necessary to know which variable group it belongs to */
	  EDGEDATA elementData;
	  PetscInt key, offsetd_eStart;


	  // Index for alphaMaxarr
	  unsigned int k = 0;

	  for (e=eStart; e < eEnd; e++) {
	    PetscInt    offsetd,key;

		const PetscInt *cone;
		DMNetworkGetConnectedNodes(networkdm,e,&cone);
		vfrom = cone[0];
		vto   = cone[1];


		ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
		// Cast component data
		beadsdataFrom = (arrayBeads)(arr + offsetd);

		ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
		// Cast component data
		beadsdataTo = (arrayBeads)(arr + offsetd);


		ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
		ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

		// Grab displacement components
		U_i_from = uarr[offsetfrom + 2];
		U_j_from = uarr[offsetfrom + 3];

		U_i_to = uarr[offsetto + 2];
		U_j_to = uarr[offsetto + 3];

		// Reset the beads for the constitutive law object
		sistema->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, sistema->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

		// Set initial design value, first, we need to know which variable group it belongs to
		DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd_eStart);
		elementData = (EDGEDATA)(arr + offsetd_eStart);
		sistema->elementFunction.set_alphaP_initial(sistema->alphaMaxInit[elementData->VariableGroup]);

		// Set the value of alphaMax. We need to do this to update
		sistema->elementFunction.set_alphaMax(alphaMaxarr[k]);
		sistema->elementFunction.calculate_delta();
		PetscReal * force_vector = sistema->elementFunction.force_vector();
		alphaMaxarr[k] = sistema->elementFunction.get_alphaMax();
		k++;

#ifdef DEBUG
//		std::cout<<"delta = "<<sistema->elementFunction.get_delta()<<" element = "<<k-1<<std::endl;
//		std::cout<<"_alphaMax = "<<alphaMaxarr[k-1]<<" element = "<<k-1<<std::endl;
#endif


#ifdef VERIF

		PetscScalar force_norm = sqrt(pow(force_vector[0],2.0) + pow(force_vector[1],2.0));
		VecSetValue(sistema->Felement,e,force_norm,INSERT_VALUES);
#endif
	  }

	  // Restore the modified vectors
	  ierr = VecRestoreArray(sistema->alphaMax,&alphaMaxarr);
	  ierr = VecAssemblyBegin(sistema->alphaMax);
	  ierr = VecAssemblyEnd(sistema->alphaMax);

#ifdef VERIF
	  ierr = VecAssemblyBegin(sistema->Felement);
	  ierr = VecAssemblyEnd(sistema->Felement);
#endif

	  ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
	  ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);

	  /*
	   * Save current time
	   */
	  sistema->timeHistory.push_back(ptime);


	  /*
	   * Print the vector to a matlab file
	   */
#ifdef VERIF
	  // File name
	  PetscViewer    viewer;
	  std::string filehistory("ForceHistory/");
	  std::ostringstream temp;
	  temp << step;
	  filehistory.append(temp.str());
	  filehistory.append(".mat");
	  // Open and set Binary file
	  PetscViewerBinaryOpen(PETSC_COMM_WORLD,filehistory.c_str(),FILE_MODE_WRITE,&viewer);
	  PetscViewerSetFormat(viewer, PETSC_VIEWER_BINARY_MATLAB);
	  // Print
	  VecView(sistema->Felement,viewer);
	  PetscViewerDestroy(&viewer);
#endif


	  /*
	   * Saving displacements
	   */
	  if (!sistema->solve_in_fd && sistema->problemData.printdisplacements == 1){

		  // File name
		  PetscViewer    viewer;
		  std::string filehistory("PostProcessing/Displacements/");
		  std::ostringstream temp;
		  temp << step;
		  filehistory.append(temp.str());
		  filehistory.append(".mat");
		  // Open and set Binary file
		  PetscViewerBinaryOpen(PETSC_COMM_WORLD,filehistory.c_str(),FILE_MODE_WRITE,&viewer);
		  PetscViewerSetFormat(viewer, PETSC_VIEWER_BINARY_MATLAB);
		  // Print
		  VecView(U,viewer);
		  PetscViewerDestroy(&viewer);
	  }

	  /*
	   * Make a copy of the current alphaMax Vec for the history array.
	   * This Vec belongs to the current time step. In this step, we will
	   * modify it, and it will be saved as the alphaMax of the next time step.
	   */
#ifdef PETSC_USE_LOG
	  PetscLogStagePush(sistema->solvestage);
#endif


	/*
	 * Save variable history
	 */
	  sistema->alphaMaxhistory.insert(sistema->alphaMax);
	  sistema->Uhistory.insert(U);



#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif
	  /*
	   * Print time step history
	   */

	  if (!sistema->solve_in_fd && ptime > 0){
		  PetscReal dt;
		  dt = ptime - sistema->previous_time;
		  sistema->timeStepHistory.push_back(dt);
#ifdef DEBUG
		  std::cout<<"time step saved = "<<dt<<" time = "<<ptime<<std::endl;
#endif

		  /*
		   * Print step to file only if final evaluation
		   */
		  if (sistema->problemData.printstephistory){
			  sistema->steps_history<<std::setprecision(9)<<dt;
			  sistema->steps_history<<"\n";
		  }


#ifdef DEBUG
		  std::cout<<"time step = "<<dt<<std::endl;
#endif
	  }


	  sistema->previous_time = ptime;


	  /*
	   * Print forces at the target area if necessary
	   */

	  if (sistema->designData.printforces){
		  sistema->PrintForcesTargetArea(U,sistema->alphaMax,step);
	  }


	return(0);
}

System::System(){

	elementFunction.set_data(&problemData,&designData);
	_is_partial_derivatives_calculated = false;


	_is_fd_calculated = false;
	_is_gradient_calculated = false;

	SaveStages = true;

	Mdot_product_result = NULL;

	FunctionValueGlobalOriginal = 0.0;


	designData.printforces = false;


}

System::~System()
{

#ifdef PETSC_USE_LOG
	  PetscLogStagePush(solvestage);
#endif
	Uhistory.destroyhistory();
	alphaMaxhistory.destroyhistory();
	U_RungeKuttahistory.destroyhistory();

	VecDestroy(&X);
	VecDestroy(&F);

	VecDestroy(&Felement);
	VecDestroy(&alphaMax);



	TSDestroy(&ts);
	DMDestroy(&networkdm);

#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif

#ifdef PETSC_USE_LOG
	  PetscLogStagePush(adjointstage);
#endif
	/*
	 * Destroying vectors used in the Adjoint Analysis
	 */
	VecDestroy(&Lambda);
	VecDestroy(&Gamma);

	VecDestroy(&Omega);
	VecDestroy(&Mu);

	VecDestroy(&TMP);
	VecDestroy(&TMPalpha);

	VecDestroy(&dFdalphaMax);
	VecDestroy(&dFdU);
	VecDestroy(&dFdP);

	VecDestroy(&dFdalphaMax_partial);
	VecDestroy(&dFdU_partial);
	VecDestroy(&dFdP_partial);


	VecDestroyVecs(Stages,&Phi);
	VecDestroyVecs(Stages,&K_RK);
	VecDestroyVecs(Stages,&U_RK);
	VecDestroyVecs(Stages,&V_i);
	VecDestroyVecs(Stages,&W_i);

	delete [] parameters_indices;

	VecDestroyVecs(designData.N_DesignVariables,&JACP);

	VecDestroy(&DN);

#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif

	delete [] c;
	delete [] b;
	delete [] a;
	delete [] ones;
	delete [] Mdot_product_result;
	PetscFinalize();
}

PetscErrorCode System::InitializeData()
{

	 PetscErrorCode       ierr;


	  PetscLogStage        stage1;

	  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

	  /*
	   * First, delete the contents of the folder ForceHistory
	   */

#ifdef VERIF
	  if (!rank){
		  system("exec rm -r ForceHistory/*");
	  }
#endif

	  if (!rank){
		  system("exec rm -r PostProcessing/ForcesTargetArea/*");
		  system("exec rm -r PostProcessing/Displacements/*");
	  }
	  /*
	   * Initialize data
	   */
	  elementFunction.set_data(&problemData,&designData);


	  problemData.peak_load = 22;
	  problemData.loadtime =10;
	  ierr =PetscOptionsReal("-peak_load", "Introduce peak load value in KN","", 22, &problemData.peak_load, NULL);
	  ierr =PetscOptionsReal("-loadtime", "Introduce load time in ms","", 10, &problemData.loadtime, NULL);



	  /* Create an empty network object */
	  ierr = DMNetworkCreate(PETSC_COMM_WORLD,&networkdm);CHKERRQ(ierr);
	  /* Register the components in the network */
	  ierr = DMNetworkRegisterComponent(networkdm,"beadsdata",sizeof(class BeadsData),&componentkey[0]);CHKERRQ(ierr);
	  ierr = DMNetworkRegisterComponent(networkdm,"globalindexedge",sizeof(struct _p_EDGEDATA),&componentkey[1]);CHKERRQ(ierr);

	  ierr = PetscLogStageRegister("Read Data",&stage1);CHKERRQ(ierr);
	  PetscLogStagePush(stage1);


	  /*
	   * READ THE DATA
	   */

	  ierr =PetscOptionsInt("-analysis", "0 for linear model, 1 for non-linear elastic, 2 for plastic","", 0, &problemData.analysis, NULL);
	  char                 dynamic[PETSC_MAX_PATH_LEN];
	  ierr = PetscOptionsGetString(PETSC_NULL,"-dynamic",dynamic,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);
	  std::string transient(dynamic);
	  if (transient.compare("1") == 0){
		problemData.transient = true;
	  }
	  else{
		problemData.transient = false;
	  }

	  if (rank == 0) {
	    /*    READ DATA */
	    /* Only rank 0 reads the data */
	    //ierr = PetscOptionsGetString(PETSC_NULL,"-pfdata",pfdata_file,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);
	    ierr = ReadBeadsData(allbeadsData);
	    numEdges = problemData.N_Elements;
	    numVertices = problemData.N_Nodes;

	  }



	  /*
	   * Set up stages
	   */
	  SetUpStages();

	  /*
	   *
	   */
	  ierr =PetscOptionsInt("-objfunction", "0 for simple obj function (sum of the displacements), 1 for the p-norm of the force","", 0, &designData.objfunctionType, NULL);
	  ierr =PetscOptionsInt("-loading", "0 for prescribed displacement, 1 for loading","", 0, &problemData.loading, NULL);
	  ierr =PetscOptionsInt("-matlabhistory", "0 for ignoring the time step history from matlab, 1 for using it","", 0, &problemData.matlabhistory, NULL);


	  if (designData.objfunctionType == 1){

		  CalculateObjFunctionTime = &System::CalculateObjFunctionTime_P_Norm;
		  PartialImplicitDerivatives = &System::EvaluatePartialImplicitDerivatives_P_Norm;
		  PartialExplicitDerivatives = &System::EvaluatePartialExplicitDerivatives_P_Norm;
	  }
	  else{

		  CalculateObjFunctionTime = &System::CalculateObjFunctionTime_Simple;
		  PartialImplicitDerivatives = &System::EvaluatePartialImplicitDerivatives_Simple;
		  PartialExplicitDerivatives = &System::EvaluatePartialExplicitDerivatives_Simple;
	  }



	return ierr;
}

PetscErrorCode System::SetUpStages(){
	  PetscErrorCode       ierr;

#ifdef PETSC_USE_LOG
	  ierr = PetscLogStageRegister("Primal problem",&solvestage);CHKERRQ(ierr);

	  ierr = PetscLogStageRegister("Adjoint problem",&adjointstage);CHKERRQ(ierr);

#endif
	  return ierr;

}

PetscErrorCode System::SetUpDMNetwork()
{
	 /*
	  * Set up the network and add components
	  */

	  PetscLogStage stage2;
	  PetscErrorCode       ierr;
	  PetscInt             i;
	  PetscInt             size;

	  PetscLogStagePop();
	  MPI_Barrier(PETSC_COMM_WORLD);
#ifdef PETSC_USE_LOG
	  ierr = PetscLogStageRegister("Create network",&stage2);CHKERRQ(ierr);
	  PetscLogStagePush(stage2);
#endif
	  /* Set number of nodes/edges */
	  ierr = DMNetworkSetSizes(networkdm,numVertices,numEdges,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
	  /* Add edge connectivity */
	  ierr = DMNetworkSetEdgeList(networkdm,problemData.connectivity);CHKERRQ(ierr);
	  /* Set up the network layout */
	  ierr = DMNetworkLayoutSetUp(networkdm);CHKERRQ(ierr);

	  /* Add network components */
	  PetscInt vStart, vEnd;

	  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);

	  for (i = vStart; i < vEnd; i++) {
	    ierr = DMNetworkAddComponent(networkdm,i,componentkey[0],&allbeadsData[i-vStart]);CHKERRQ(ierr);
	    /* Add number of variables */
	    if (problemData.transient == false){
	    	ierr = DMNetworkAddNumVariables(networkdm,i,2);CHKERRQ(ierr);
	    }
	    else{
	    	ierr = DMNetworkAddNumVariables(networkdm,i,4);CHKERRQ(ierr);
	    }
	  }

	  /*
	   * Add degrees of freedom corresponding to the objective functions
	   * It will be added to the first N edges, where N is the number of objective function
	   */
	  PetscInt eStart, eEnd;
	  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
	  if (eStart == 0){ // Conditional so we only add the new DOF in the processor where the first element is
		  PetscInt N = designData.N_ObjFunc;
		  ierr = DMNetworkAddNumVariables(networkdm,eStart,N);CHKERRQ(ierr);
	  }

	  /*
	   * Add a component to each edge that will be the global index;
	   */
	  for (i = eStart; i < eEnd; i++) {
	    ierr = DMNetworkAddComponent(networkdm,i,componentkey[1],&edgeData[i]);CHKERRQ(ierr);
	  }

	  /* Set up DM for use */
	  ierr = DMSetUp(networkdm);CHKERRQ(ierr);

	  /* We can delete the content of allbeadsData once it is set inside the DMNetwork*/
	  if (!rank){
	  	  PetscFree(allbeadsData);
	  	  PetscFree(edgeData);
	  }


	  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	  if (size > 1) {
	    DM distnetworkdm;
	    /* Network partitioning and distribution of data */
	    ierr = DMNetworkDistribute(networkdm,0,&distnetworkdm);CHKERRQ(ierr);
	    ierr = DMDestroy(&networkdm);CHKERRQ(ierr);
	    networkdm = distnetworkdm;
	  }

	  /*************************************************************************/
	  /*            Assigning initial values for the variables                 */
	  /*************************************************************************/
	  PetscInt    offsetd;
	  DMNetworkComponentGenericDataType *arr;
	  arrayBeads beadsdata;
	  // We need to loop over the local vertices
	  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
	  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	  std::ifstream finitialValue("input_InitialVariableValues.txt");

	  for (unsigned int i = 0; i<designData.N_DesignVariables; i++){

		  PetscScalar design_value = 0;

			if (designData.initial_design != 0)
				design_value = designData.initial_design;
			else{
				if (finitialValue.good())
					finitialValue>>design_value;
			}
		  alphaMaxInit.push_back(design_value);


	  }

//	  std::cout<<"size alphaMax = "<<alphaMaxInit.size()<<std::endl;
//	  std::cout<<"N_DesignVariables = "<<designData.N_DesignVariables<<std::endl;
//
//	  std::cout<<"elemtns = "<<problemData.N_Elements<<std::endl;




#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif


	  return ierr;

}

PetscErrorCode System::CreatePETScObjects()
{
	PetscErrorCode ierr;
#ifdef PETSC_USE_LOG
	  PetscLogStagePush(solvestage);
#endif
	ierr = DMCreateGlobalVector(networkdm,&X);CHKERRQ(ierr);
	ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

	PetscBool  flag_algorithm = PETSC_FALSE, flag_RK = PETSC_FALSE;

	char                 algorithm[PETSC_MAX_PATH_LEN] = "rk";
	char                 RK_method[PETSC_MAX_PATH_LEN] = "5dp";
	ierr = PetscOptionsGetString(PETSC_NULL,"-ts_type",algorithm,PETSC_MAX_PATH_LEN-1,&flag_algorithm);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(PETSC_NULL,"-ts_rk_type",RK_method,PETSC_MAX_PATH_LEN-1,&flag_RK);CHKERRQ(ierr);


	std::string algorithm_type(algorithm);

	std::string method_RK(RK_method);


	/*
	 * If the algorithm was not specified, by default we pick RK BS. Same case if we don't specify the RK method
	 */
	if (flag_algorithm){
		if (algorithm_type.compare("euler") == 0){
			Stages = 2;
		}
		else{
			if (flag_RK){
				if (method_RK.compare("5dp") == 0){
					Stages = 7;
				}
				else{
					Stages = 4;
				}
			}
			else
				Stages = 4;
		}
	}
	else{
		Stages = 4;
	}

	U_RungeKuttahistory.SetStages(Stages);




	// Create vectors for the spring variables
	VecCreate(PETSC_COMM_WORLD,&Felement);
	VecSetType(Felement,VECMPI);
	VecSetSizes(Felement,PETSC_DECIDE,problemData.N_Elements);
	VecDuplicate(Felement,&alphaMax);

	InitialDesignSolution();
	// We need to save the very first plastic deformation state
	//this->alphaMaxhistory.insert(alphaMax);

#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif

#ifdef PETSC_USE_LOG
	  PetscLogStagePush(adjointstage);
#endif
	/*
	 * Adjoint objects
	 */
	VecDuplicate(X,&Lambda);
	VecDuplicate(X,&Omega);
	VecDuplicateVecs(X,Stages,&Phi);
	VecDuplicateVecs(X,Stages,&K_RK);
	VecDuplicateVecs(X,Stages,&U_RK);
	VecDuplicate(alphaMax,&Gamma);

	VecDuplicate(Omega,&TMP);
	VecDuplicate(Gamma,&TMPalpha);


	VecDuplicateVecs(alphaMax, Stages, &W_i);

	VecDuplicate(X,&dFdU);
	VecDuplicate(X,&dFdU_partial);

	VecDuplicate(alphaMax,&dFdalphaMax);
	VecDuplicate(alphaMax,&dFdalphaMax_partial);

	VecCreate(PETSC_COMM_WORLD,&Mu);
	VecSetType(Mu,VECMPI);
	VecSetSizes(Mu,PETSC_DECIDE,designData.N_DesignVariables);
	VecDuplicateVecs(Mu,Stages,&V_i);

	VecDuplicate(Mu,&dFdP);
	VecDuplicate(Mu,&dFdP_partial);


	/*
	 * Number of parameters in this processor
	 */
	VecGetLocalSize(Mu,&local_N_parameters);
	parameters_indices = new PetscInt[local_N_parameters];


	/*
	 * RHS derivatives
	 */
	VecDuplicateVecs(X,local_N_parameters,&JACP);
	VecDuplicate(X,&DN);

#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif


	return ierr;

}

PetscErrorCode System::SetUpSolver()
{
	PetscErrorCode ierr;
	PetscInt       time_steps_max = 5000000;

	previous_time = 0;

#ifdef PETSC_USE_LOG
	  PetscLogStagePush(solvestage);
#endif

	ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
	ierr = TSSetType(ts,TSRK); CHKERRQ(ierr);
	ierr = TSSetDM(ts,networkdm);CHKERRQ(ierr);
	ierr = TSSetRHSFunction(ts,NULL,FormRHSFunction,this);CHKERRQ(ierr);
	ierr = TSSetDuration(ts,time_steps_max,problemData.TotalTime);
	ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
	ierr = TSMonitorSet(ts, MyTSMonitor, this, NULL);CHKERRQ(ierr);
	ierr = TSSetInitialTimeStep(ts,0.0,problemData.timestep);

	//timeStepHistory.push_back(problemData.timestep);

	ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);

	/*
	 * Set Initial solution
	 */
	InitialSolution(networkdm,X);
	ierr = TSSetSolution(ts,X);CHKERRQ(ierr);

	/*
	 * Set initial plastic deformation to initial deformation
	 */
	InitialDesignSolution();

	/*
	 * Error tolerances
	 */
	PetscReal abstol = 1e-5, reltol = 1e-5;

	PetscOptionsGetReal(NULL,"-abstol",&abstol,NULL);
	PetscOptionsGetReal(NULL,"-reltol",&reltol,NULL);
	ierr = TSSetTolerances(ts,abstol,NULL,reltol,NULL);

	solve_in_fd = false;

#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif

		if (problemData.matlabhistory == 1){
			std::fstream myfile ("matlab_timestep.txt", std::ios_base::in);
			PetscReal temp,totaltime = 0.0;
			total_timesteps_matlab = 0;
			while ( myfile >> temp)
			{
			  matlab_timestephistory.push_back(temp);
			  total_timesteps_matlab++;
			  totaltime += temp;
			  std::cout<<"h = "<<temp<<std::endl;
			}
			std::cout<<"Total time = "<<totaltime<<std::endl;
			std::cout<<"N of time steps = "<<total_timesteps_matlab<<std::endl;
			matlab_timestephistory_iterator = matlab_timestephistory.begin();

			ierr = TSSetDuration(ts,total_timesteps_matlab,problemData.TotalTime);
		}


		PetscBool  flag_algorithm = PETSC_FALSE, flag_RK = PETSC_FALSE;

		char                 algorithm[PETSC_MAX_PATH_LEN] = "rk";
		char                 RK_method[PETSC_MAX_PATH_LEN] = "5dp";
		ierr = PetscOptionsGetString(PETSC_NULL,"-ts_type",algorithm,PETSC_MAX_PATH_LEN-1,&flag_algorithm);CHKERRQ(ierr);
		ierr = PetscOptionsGetString(PETSC_NULL,"-ts_rk_type",RK_method,PETSC_MAX_PATH_LEN-1,&flag_RK);CHKERRQ(ierr);


		std::string algorithm_type(algorithm);

		std::string method_RK(RK_method);


		if (flag_algorithm){
			if (algorithm_type.compare("euler") == 0){
				c = new PetscScalar[Stages];
				b = new PetscScalar[Stages];
				a = new PetscScalar[(Stages)*(Stages)];
				ones = new PetscScalar[Stages];

				c[0] = 0.0;
				c[1] = 0.0;

				b[0] = 1.0;
				b[1] = 0.0;

				a[0] = 0.0;
				a[1] = 0.0;
				a[2] = 0.0;
				a[3] = 0.0;

				ones[0] = 1.0;
				ones[1] = 0.0;
			}
			else{
				if (flag_RK){
					if (method_RK.compare("5dp") == 0){
						b = new PetscScalar[Stages];
						c = new PetscScalar[Stages];
						a = new PetscScalar[Stages*Stages];
						ones = new PetscScalar[Stages];


						a[0] = 0;
						a[1] = 0;
						a[2] = 0;
						a[3] = 0;
						a[4] = 0;
						a[5] = 0;
						a[6] = 0;

						a[7] = 1.0/5.0;
						a[8] = 0;
						a[9] = 0;
						a[10] = 0;
						a[11] = 0;
						a[12] = 0;
						a[13] = 0;

						a[14] = 3.0/40;
						a[15] = 9.0/40;
						a[16] = 0;
						a[17] = 0;
						a[18] = 0;
						a[19] = 0;
						a[20] = 0;

						a[21] = 44.0/45.0;
						a[22] = -56.0/15.0;
						a[23] = 32.0/9.0;
						a[24] = 0;
						a[25] = 0;
						a[26] = 0;
						a[27] = 0;

						a[28] = 19372.0/6561.0;
						a[29] = -25360.0/2187.0;
						a[30] = 64448.0/6561.0;
						a[31] = -212.0/729.0;
						a[32] = 0;
						a[33] = 0;
						a[34] = 0;

						a[35] = 9017.0/3168.0;
						a[36] = -355.0/33.0;
						a[37] = 46732.0/5247.0;
						a[38] = 49.0/176.0;
						a[39] = -5103.0/18656.0;
						a[40] = 0;
						a[41] = 0;

						a[42] = 35.0/384.0;
						a[43] = 0;
						a[44] = 500.0/1113.0;
						a[45] = 125.0/192.0;
						a[46] = -2187.0/6784.0;
						a[47] = 11.0/84.0;
						a[48] = 0;


						b[0] = 35.0/384.0;
						b[1] = 0;
						b[2] = 500.0/1113.0;
						b[3] = 125.0/192.0;
						b[4] = -2187.0/6784.0;
						b[5] = 11.0/84.0;
						b[6] = 0;

						for (unsigned int i = 0; i<Stages; i++){
							c[i] = 0.0;
							for (unsigned int j = 0; j<Stages; j++)
								c[i] += a[i*Stages+j];
						}

						ones[0] = 1.0;
						ones[1] = 1.0;
						ones[2] = 1.0;
						ones[3] = 1.0;
						ones[4] = 1.0;
						ones[5] = 1.0;
						ones[6] = 1.0;
					}
					else{
						SetRKDefault(*this);
					}
				}
				else
					SetRKDefault(*this);
			}
		}
		else{
			SetRKDefault(*this);
		}


	return ierr;


}

PetscErrorCode System::FDSetUpSolver()
{
	PetscErrorCode ierr;

#ifdef PETSC_USE_LOG
	  PetscLogStagePush(solvestage);
#endif

	ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);

	/*
	 * Set Initial solution
	 */
	InitialSolution(networkdm,X);
	ierr = TSSetSolution(ts,X);CHKERRQ(ierr);

	/*
	 * Set initial plastic deformation to the design values
	 */
	InitialDesignSolution();

	/*
	 * Error tolerances
	 */
	PetscReal abstol = 1e-5, reltol = 1e-5;

	PetscOptionsGetReal(NULL,"-abstol",&abstol,NULL);
	PetscOptionsGetReal(NULL,"-reltol",&reltol,NULL);
	ierr = TSSetTolerances(ts,abstol,NULL,reltol,NULL);


	if (problemData.matlabhistory != 1){
		TSAdapt adapt;
		TSGetAdapt(ts,&adapt);
		if (!adapt)
			TSAdaptSetType(adapt,"none");

		solve_in_fd = true;

		timeStepHistory_iterator = timeStepHistory.begin();

		timeStepHistory_iterator_end = timeStepHistory.end();
	}
	ierr = TSSetDuration(ts,(PetscInt)timeStepHistory.size(),problemData.TotalTime);




#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif

//		if (problemData.matlabhistory == 1){
//			std::fstream myfile ("matlab_timestep.txt", std::ios_base::in);
//			PetscReal temp;
//			total_timesteps_matlab = 0;
//			while ( myfile >> temp)
//			{
//			  matlab_timestephistory.push_back(temp);
//			  total_timesteps_matlab++;
//			  std::cout<<"h = "<<temp<<std::endl;
//			}
//			matlab_timestephistory_iterator = matlab_timestephistory.begin();
//
//			ierr = TSSetDuration(ts,total_timesteps_matlab,problemData.TotalTime);
//		}

	return ierr;


}

PetscErrorCode System::Solve()
{
	PetscErrorCode ierr;
	PetscInt       steps;
	PetscReal ftime;
#ifdef PETSC_USE_LOG
	  PetscLogStagePush(solvestage);
#endif
	if (!solve_in_fd && problemData.printstephistory){
		// Open file to save the steps history
		steps_history.open("PostProcessing/Steps_history.txt",std::ofstream::out | std::ofstream::trunc);
	}

	start = clock();

	ierr = TSSolve(ts,X);
	ierr = TSGetTimeStepNumber(ts,&steps);

	end = clock();

	totalcputime = ((float)(end - start)) / CLOCKS_PER_SEC ;

	std::cout<<"Primal Problem time = "<<totalcputime<<std::endl;

	TSGetSolveTime(ts,&ftime);
	//PetscPrintf(PETSC_COMM_WORLD,"steps %D, ftime %g\n",steps,ftime);

	// Print total number of steps
	std::ofstream N_steps;
	N_steps.open("PostProcessing/N_steps.txt");
	N_steps<<steps<<" "<<ftime;
	N_steps.close();

	if (!solve_in_fd && problemData.printstephistory){
		steps_history.close();
	}

#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif

	return ierr;

}

PetscErrorCode System::RHSParameterDerivatives(PetscReal & t, Vec &  U, Vec & DF,  std::vector<PetscInt> & DesignElements, Vec & alphaMax){

	  PetscErrorCode ierr;
	  Vec           localU,localDF, localDN;
	  PetscInt      e;
	  PetscInt      eStart,eEnd,vfrom,vto;
	  const PetscScalar *uarr;
	  PetscScalar U_i_from, U_j_from, U_i_to, U_j_to;
	  PetscScalar   *farr, *alphaMaxarr, *dnarr;
	  PetscInt      offsetfrom,offsetto;
	  DMNetworkComponentGenericDataType *arr;

	  arrayBeads beadsdataFrom, beadsdataTo ;

	  PetscReal mass, mass_derivative;

	  /*
	   * We need to create a copy of the variable vector for the second contribution
	   * to the RHS derivative
	   */





	  PetscFunctionBeginUser;


	  ierr = VecSet(DF,0.0);CHKERRQ(ierr);
	  ierr = VecSet(DN,0.0);CHKERRQ(ierr);

	  ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
	  ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	  ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	  ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);


	  ierr = DMGetLocalVector(networkdm,&localDN);CHKERRQ(ierr);
	  ierr = DMGlobalToLocalBegin(networkdm,DN,INSERT_VALUES,localDN);CHKERRQ(ierr);
	  ierr = DMGlobalToLocalEnd(networkdm,DN,INSERT_VALUES,localDN);CHKERRQ(ierr);
	  ierr = VecGetArray(localDN,&dnarr);CHKERRQ(ierr);




	  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
	  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	  /*
	   * Iterate through the elements
	   */

	  std::vector<PetscInt>::iterator first_element = DesignElements.begin();
	  const std::vector<PetscInt>::iterator last_element = DesignElements.end();

	  /* Edge Data necessary to know which variable group it belongs to */
	  EDGEDATA elementData;
	  PetscInt key, offsetd_eStart;

	  VecGetArray(alphaMax,&alphaMaxarr);


	  unsigned int k = 0;
	  for (; first_element != last_element; ++first_element) {
		  PetscInt    offsetd,key;

		  e = *first_element - 1;

		  const PetscInt *cone;
		  DMNetworkGetConnectedNodes(networkdm,e,&cone);
		  vfrom = cone[0];
		  vto   = cone[1];


		  ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
		  // Cast component data
		  beadsdataFrom = (arrayBeads)(arr + offsetd);

		  ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
		  // Cast component data
		  beadsdataTo = (arrayBeads)(arr + offsetd);

		  ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
		  ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);



		  // Grab displacement components
		  U_i_from = uarr[offsetfrom + 2];
		  U_j_from = uarr[offsetfrom + 3];

		  U_i_to = uarr[offsetto + 2];
		  U_j_to = uarr[offsetto + 3];

		  // Reset the beads for the constitutive law object
		  this->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, this->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

		  // Set initial design value, first, we need to know which variable group it belongs to
		  DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd_eStart);
		  elementData = (EDGEDATA)(arr + offsetd_eStart);
		  this->elementFunction.set_alphaP_initial(alphaMaxInit[elementData->VariableGroup]);

		  // Set the value of alphaMax


		  this->elementFunction.set_alphaMax(alphaMaxarr[e]);

		  this->elementFunction.calculate_delta();
		  k++;


		  PetscReal * internal_force_derivative = this->elementFunction.internal_force_partial_derivative_P_bead();
		  dnarr[offsetfrom] += -1.0*internal_force_derivative[0];
		  dnarr[offsetfrom + 1] += -1.0*internal_force_derivative[1];

		  dnarr[offsetto] += -1.0*internal_force_derivative[2];
		  dnarr[offsetto + 1] += -1.0*internal_force_derivative[3];

		  /*
		   * Add the constrained degrees of freedom contribution
		   */
			if (beadsdataFrom->constrained_dof_x){
				dnarr[offsetfrom] = 0.0;
				dnarr[offsetfrom + 2] = 0.0;
			}
			if (beadsdataFrom->constrained_dof_y){
				dnarr[offsetfrom + 1] = 0.0;
				dnarr[offsetfrom + 3] = 0.0;
			}
			if (beadsdataTo->constrained_dof_x){
				dnarr[offsetto] = 0.0;
				dnarr[offsetto + 2] = 0.0;
			}
			if (beadsdataTo->constrained_dof_y){
				dnarr[offsetto + 1] = 0.0;
				dnarr[offsetto + 3] = 0.0;
			}

	  }

	  /*
	   * We need to loop over the beads to multiply the farr vector (only the velocity components)
	   * with the inverse of the mass. We need to have the contributions from all elements before
	   * multiplying by this factor
	   */
	  PetscInt vStart, vEnd, v;
	  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
	  for (v=vStart; v < vEnd; v++) {
	    PetscInt    offsetd,key;
		ierr = DMNetworkGetComponentTypeOffset(networkdm,v,0,&key,&offsetd);CHKERRQ(ierr);
		// Cast component data
		beadsdataFrom = (arrayBeads)(arr + offsetd);

		ierr = DMNetworkGetVariableOffset(networkdm,v,&offsetfrom);CHKERRQ(ierr);

		beadsdataFrom->get_mass(mass);
		dnarr[offsetfrom] *= 1.0/mass;
		dnarr[offsetfrom + 1] *= 1.0/mass;
	  }


	  VecRestoreArray(alphaMax,&alphaMaxarr);



	  ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);

	  ierr = VecRestoreArray(localDN,&dnarr);CHKERRQ(ierr);
	  ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);



	  ierr = DMLocalToGlobalBegin(networkdm,localDN,ADD_VALUES,DN);CHKERRQ(ierr);
	  ierr = DMLocalToGlobalEnd(networkdm,localDN,ADD_VALUES,DN);CHKERRQ(ierr);
	  ierr = DMRestoreLocalVector(networkdm,&localDN);CHKERRQ(ierr);

	  /*
	   * Sum both contributions
	   *
	   * VecAXPBY(Vec y,PetscScalar alpha,PetscScalar beta,Vec x)
	   * Computes y = alpha x + beta y.
	   */

	  VecAXPY(DF,1.0,DN);


	  PetscFunctionReturn(0);



	return ierr;


}

PetscErrorCode System::RHSParameterDerivatives_FD(PetscReal & t, Vec & U, Vec & DF, PetscInt Design_Index, Vec & alphaMax){
	  PetscErrorCode ierr;

	  assert(_is_partial_derivatives_calculated && "Need to calculate the derivatives first");

	  PetscReal original_design;
	  PetscReal dt = 1e-8;

	  // Vector for the FD partial derivatives for this time step
	  Vec Residual_partialP_FD, Residual_partialP_FD_TWO;
	  VecDuplicate(U,&Residual_partialP_FD);

	  VecDuplicate(U,&Residual_partialP_FD_TWO);



	  // Zero the vectors for a new calculation
	  VecSet(Residual_partialP_FD,0.0);
	  VecSet(Residual_partialP_FD_TWO,0.0);


	  // Change parameter data, first, keep original value
	  original_design = alphaMaxInit[Design_Index];


	  alphaMaxInit[Design_Index] += dt;

	  // Evaluate function
	  FormRHSFunction(ts,t,U,Residual_partialP_FD,this);

	  // Change the parameter again, now backwards
	  alphaMaxInit[Design_Index] -= 2.0*dt;


	  // Evaluate function
	  FormRHSFunction(ts,t,U,Residual_partialP_FD_TWO,this);

	  // Final value
	  VecAXPY(Residual_partialP_FD,-1.0,Residual_partialP_FD_TWO);
	  PetscScalar alpha = 1.0/(2.0*dt);
	  VecScale(Residual_partialP_FD,alpha);

	  // Recover density value
	  alphaMaxInit[Design_Index]= original_design;


	  std::cout<<"/*-----------Comparing residual derivatives -------------*/"<<std::endl;
	  std::cout<<"/*-- The first DOF of the FD is the objective function --*/"<<std::endl;
	  std::cout<<"/*--------------------Analytical -----------------------*/"<<std::endl;
	  VecView(DF,PETSC_VIEWER_STDOUT_WORLD);
	  std::cout<<"/*--------------------------FD -------------------------*/"<<std::endl;
	  VecView(Residual_partialP_FD,PETSC_VIEWER_STDOUT_WORLD);



	  VecDestroy(&Residual_partialP_FD);
	  VecDestroy(&Residual_partialP_FD_TWO);



	  return ierr;
}

PetscErrorCode System::AdjointSolve(){
	PetscErrorCode ierr;
	PetscInt      e,eStart,eEnd,vfrom,vto;
	DMNetworkComponentGenericDataType *arr;
	PetscInt      offsetfrom,offsetto;
	arrayBeads beadsdataFrom, beadsdataTo;

	unsigned int k;

	start = clock();


	Vec           localU, localOmega, localLambda;
	const PetscScalar *uarr;
	PetscScalar  * omegaarr, * tmparr, * tmpalphamaxarr, * gammaPrevarr, * lambdaarr;
	const PetscScalar *alphaMaxarr;
	PetscScalar U_i_from, U_j_from, U_i_to, U_j_to;

	/*
	 * Colvecs for the operations with the linear operator
	 * omega_local_U grabs the values of the original vector
	 * Only 4 components because it only takes care of the variable U
	 */
	arma::colvec omega_local_U(4), phi_local(4), omega_local_V(4); // Van al reves de hecho omega_local_V guarda los valores de los dofs de U

#ifdef PETSC_USE_LOG
	  PetscLogStagePush(adjointstage);
#endif

	/*
	 * Retrieve the last values of the partial derivatives, time step history and implicit variables
	 */

	std::vector<PetscReal>::reverse_iterator CurrentTime =  timeHistory.rbegin();
	std::vector<PetscReal>::reverse_iterator TimeStep = timeStepHistory.rbegin();
	std::vector<Vec>::reverse_iterator Implicit_Variable =  Uhistory.rbegin();
	std::vector<Vec>::reverse_iterator Implicit_Variable_previous =  Uhistory.rbegin();
	std::vector<Vec>::reverse_iterator Implicit_Variable_post =  Uhistory.rbegin();
	std::vector<Vec>::reverse_iterator Implicit_Variable_end =  Uhistory.rend();
	std::vector<Vec>::reverse_iterator alphaMaxVec = alphaMaxhistory.rbegin();
	std::vector<Vec>::reverse_iterator alphaMaxVec_previous = alphaMaxhistory.rbegin();
	std::vector<Vec>::reverse_iterator alphaMaxVec_end = alphaMaxhistory.rend();

	std::vector<std::vector<Vec> >::reverse_iterator Implicit_Variable_RkStage = U_RungeKuttahistory.rbegin();

#ifdef DEBUG
	std::cout<<"Size of Timesteps = "<<timeStepHistory.size()<<std::endl;
	std::cout<<"Size of Time history = "<<timeHistory.size()<<std::endl;
	std::cout<<"Size of Implicit_Variable history = "<<Uhistory.size()<<std::endl;
	std::cout<<"Size of alphaMax history = "<<alphaMaxhistory.size()<<std::endl;
#endif


	/*
	 * This last time step is actually a pointer to a time step before the very first time step
	 */
	std::vector<PetscReal>::reverse_iterator TimeStepEnd =  timeStepHistory.rend();
	CurrentTime++;
	Implicit_Variable++;
	alphaMaxVec++;
	Implicit_Variable_RkStage++;
	SaveStages = false;

	Implicit_Variable_previous++;Implicit_Variable_previous++;
	alphaMaxVec_previous++; alphaMaxVec_previous++;
	/*
	 * Start iterating backwards in time
	 */
	ierr = VecSet(Lambda,0.0);
	ierr = VecSet(Gamma,0.0);
	ierr = VecSet(Mu,0.0);

	/*
	 * Parameters indices that belong to this processor
	 */
	PetscInt firstParameter, lastParameter, index_dummy = 0;
	VecGetOwnershipRange(Mu,&firstParameter,&lastParameter);

	for (PetscInt design_number = firstParameter; design_number<lastParameter; design_number++){
		parameters_indices[index_dummy] = design_number;
#ifdef DEBUG
		std::cout<<"parameters_indices[index_dummy] = "<<parameters_indices[index_dummy]<<std::endl;
#endif
		index_dummy++;
	}
//	std::cout<<"Last Solution"<<std::endl;
//	VecView(*Uhistory.rbegin(),PETSC_VIEWER_STDOUT_WORLD);

	Vec localTMP;


	PetscScalar w = 0.0;
	PetscScalar h;

	PetscScalar result = 0;

	// Get the bounds [start, end) for the edges.
	ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

	PetscInt vStart, vEnd, v;

	ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
	ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	Vec U_RK_Stage;


	/* Edge Data necessary to know which variable group it belongs to */
	EDGEDATA elementData;
	PetscInt key, offsetd_eStart;



	unsigned int pasos = 0;

	//TimeStepEnd--;

	for(; TimeStep!=TimeStepEnd;++TimeStep, ++Implicit_Variable,++CurrentTime,++alphaMaxVec, ++Implicit_Variable_previous,
												++alphaMaxVec_previous, ++Implicit_Variable_RkStage, ++Implicit_Variable_post ){
		pasos++;
#ifdef DEBUG
		std::cout<<"Current time = "<<*CurrentTime<<std::endl;
		std::cout<<"Lambda"<<std::endl;
		VecView(Lambda,PETSC_VIEWER_STDOUT_WORLD);
		std::cout<<"Mu"<<std::endl;
		VecView(Mu,PETSC_VIEWER_STDOUT_WORLD);
		std::cout<<"Gamma-"<<std::endl;
		VecView(Gamma,PETSC_VIEWER_STDOUT_WORLD);
#endif
//		std::cout<<"Time step = "<<*TimeStep<<std::endl;
		h = *TimeStep;

//		std::cout<<"Current time = "<<*CurrentTime<<std::endl;
//		  std::cout<<"Solution"<<std::endl;
//		  VecView(*Implicit_Variable,PETSC_VIEWER_STDOUT_WORLD);


		/*
		 * Get alphaMax info
		 */
		VecGetArrayRead(*alphaMaxVec,&alphaMaxarr);


		  if (problemData.finiteDifference){
			  /*
			   * We need to copy the current alphaMaxVec to the sistema->alphaMax
			   * because it will be used in FormRHSFunction
			   */
			  VecCopy(*alphaMaxVec,alphaMax);
		  }



		/*
		 * Reinitialize the contributions of the RHS derivatives w.r.t. the parameter
		 */
		for (PetscInt i = 0; i<local_N_parameters; i++){
			VecSet(JACP[i],0.0);
		}

//		for (unsigned int i = 0; i<Stages; i++){
//			VecSet(K_RK[i],0.0);
//		}
//
//
//		for (unsigned int i = 0; i<Stages; i++){
//
//			/*
//			 * Copy the value of the current displacement into the variable U_RK
//			 */
//			VecCopy(*Implicit_Variable,U_RK[i]);
//			/*
//			 * Get the weight for the evaluation of the U_RK
//			 */
//			if (i > 0){
//				for (unsigned int j=0; j<i; j++){
//					w = h*a[i*Stages+j];
//
//					std::cout<<"a = "<<a[i*Stages+j]<<std::endl;
//
//					/*
//					 * VecAXPY(Vec y,PetscScalar alpha,Vec x)
//					 * y = alpha x + y.
//					 */
//					VecAXPY(U_RK[i],w,K_RK[j]);
//
//				}
//			}
//
//			/*
//			 * We need to evaluate the RHS Function to obtain K_RK for the next stage
//			 */
//			PetscScalar time_RK = *CurrentTime + h*c[i];
//
//			std::cout<<"From primal"<<std::endl;
//			VecView((*Implicit_Variable_RkStage)[i], PETSC_VIEWER_STDOUT_WORLD);
//
//			std::cout<<"Calculated"<<std::endl;
//			VecView(U_RK[i],PETSC_VIEWER_STDOUT_WORLD);
//
//			ierr = FormRHSFunction(ts,time_RK,U_RK[i],K_RK[i],this);CHKERRQ(ierr);
//
//
//
//		}







		for (unsigned int i = 0; i< Stages - 1; i++){
			VecSet(V_i[i],0.0);
			VecSet(Phi[i],0.0);
			VecSet(W_i[i],0.0);
		}


		for (int i = Stages - 2; i>=0; i--){


//			std::cout<<"------------------------------"<<std::endl;
//			std::cout<<"--------Stage "<<i<<"---------"<<std::endl;
//			std::cout<<"------------------------------"<<std::endl;
			U_RK_Stage = (*Implicit_Variable_RkStage)[i];

			/*
			 * Call partial derivatives of the objective function
			 */

			VecSet(dFdU,0.0);
			VecSet(dFdalphaMax,0.0);
			VecSet(dFdP,0.0);

			for (PetscInt NumbObjFunction = 0; NumbObjFunction < designData.N_ObjFunc; NumbObjFunction++){
				EvaluatePartialImplicitDerivatives(U_RK_Stage,*alphaMaxVec,NumbObjFunction);
				EvaluatePartialExplicitDerivatives(U_RK_Stage,*alphaMaxVec,NumbObjFunction);

				/*
				 * Check FD for the Partial derivatives
				 */
				  if (problemData.finiteDifference){
					  std::cout<<"Function # "<<NumbObjFunction<<std::endl;
					  FD_variables(U_RK_Stage,*alphaMaxVec,NumbObjFunction);
					  FD_statevariables(U_RK_Stage,*alphaMaxVec,NumbObjFunction);
					  FD_parameter(U_RK_Stage,*alphaMaxVec,NumbObjFunction);
				  }

				  /*
				   * Add the contribution of this objective function
				   */

				  VecAXPY(dFdU, CapitalOmega[NumbObjFunction], dFdU_partial);
				  VecAXPY(dFdalphaMax, CapitalOmega[NumbObjFunction], dFdalphaMax_partial);
				  VecAXPY(dFdP, CapitalOmega[NumbObjFunction], dFdP_partial);
			}

			/*
			 * For each stage, calculate the contribution of the RHS derivative
			 */
			PetscScalar time_RK = *CurrentTime + h*c[i];


			for (PetscInt design_number = 0; design_number<local_N_parameters; design_number++){
				RHSParameterDerivatives(time_RK,U_RK_Stage,JACP[design_number],design_elements[design_number],*alphaMaxVec);
				if (problemData.finiteDifference){
					RHSParameterDerivatives_FD(time_RK,U_RK_Stage,JACP[design_number],design_number,*alphaMaxVec);
				}
			}


			VecSet(Omega,0.0);

			if (i < (int)Stages - 2){
				for (unsigned int j=i + 1; j<Stages - 1; j++){
					w = h*a[j*Stages+i];
					/*
					 * VecAXPY(Vec y,PetscScalar alpha,Vec x)
					 * y = alpha x + y.
					 */
					VecAXPY(Omega,w,Phi[j]);

				}
			}


			/*
			 * Add the contribution from Lambda
			 */

			PetscScalar alpha = h*b[i];

			VecAXPY(Omega,alpha,Lambda);

			/*
			 * Calculate V_i [i] = JACP'*Omega + h*B(istage)*dRdP_S;
			 * We do not need to loop over the elements because JACP
			 * has already been calculated (JACP). Calculate the dot
			 * product JACP * Omega only for the indices owned by
			 * the current processor
			 */

			VecMDot(Omega,local_N_parameters,JACP,Mdot_product_result);

#ifdef DEBUG
//			for (int qq = 0; qq<local_N_parameters; qq++)
//				std::cout<<"Mdot_product_result = "<<Mdot_product_result[qq]<<std::endl;
#endif
			VecSetValues(V_i[i],local_N_parameters,parameters_indices,Mdot_product_result,ADD_VALUES);



			VecAssemblyBegin(V_i[i]);
			VecAssemblyEnd(V_i[i]);

#ifdef DEBUG
//			std::cout<<"V[i]"<<std::endl;
//			VecView(V_i[i],PETSC_VIEWER_STDOUT_WORLD);
#endif

			VecSet(TMP,0.0);
			VecSet(TMPalpha,0.0);



			/*
			 * Get local copies of the vectors
			 */
			ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
			ierr = DMGetLocalVector(networkdm,&localOmega);CHKERRQ(ierr);
			ierr = DMGetLocalVector(networkdm,&localTMP);CHKERRQ(ierr);

			ierr = DMGlobalToLocalBegin(networkdm,U_RK_Stage,INSERT_VALUES,localU);CHKERRQ(ierr);
			ierr = DMGlobalToLocalEnd(networkdm,U_RK_Stage,INSERT_VALUES,localU);CHKERRQ(ierr);

			ierr = DMGlobalToLocalBegin(networkdm,Omega,INSERT_VALUES,localOmega);CHKERRQ(ierr);
			ierr = DMGlobalToLocalEnd(networkdm,Omega,INSERT_VALUES,localOmega);CHKERRQ(ierr);

			ierr = DMGlobalToLocalBegin(networkdm,TMP,INSERT_VALUES,localTMP);CHKERRQ(ierr);
			ierr = DMGlobalToLocalEnd(networkdm,TMP,INSERT_VALUES,localTMP);CHKERRQ(ierr);


			ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);
			ierr = VecGetArray(localOmega,&omegaarr);CHKERRQ(ierr);
			ierr = VecGetArray(localTMP,&tmparr);CHKERRQ(ierr);

			ierr = VecGetArray(TMPalpha,&tmpalphamaxarr);CHKERRQ(ierr);


			/*
			 * Before calculating Phi, we're going to add the contribution
			 * of the mass to the equations. This factor comes from the
			 * jacobians:
			 *
			 * DF / DU          &&            DF / DalphaMax
			 */


			for (v=vStart; v < vEnd; v++) {
				PetscInt    offsetd,key;
			    PetscReal mass;
				ierr = DMNetworkGetComponentTypeOffset(networkdm,v,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataFrom = (arrayBeads)(arr + offsetd);

				ierr = DMNetworkGetVariableOffset(networkdm,v,&offsetfrom);CHKERRQ(ierr);

				beadsdataFrom->get_mass(mass);
				omegaarr[offsetfrom] *= 1.0/mass;
				omegaarr[offsetfrom + 1] *= 1.0/mass;
			}


			// Index for alphaMaxarr
			k = 0;

			for (e = eStart; e< eEnd; e++ ){
			    PetscInt    offsetd,key;

				const PetscInt *cone;
				DMNetworkGetConnectedNodes(networkdm,e,&cone);
				vfrom = cone[0];
				vto   = cone[1];


				ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataFrom = (arrayBeads)(arr + offsetd);

				ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataTo = (arrayBeads)(arr + offsetd);

				ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
				ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

				// Grab displacement components
				U_i_from = uarr[offsetfrom + 2];
				U_j_from = uarr[offsetfrom + 3];

				U_i_to = uarr[offsetto + 2];
				U_j_to = uarr[offsetto + 3];

				// Reset the beads for the constitutive law object
				this->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, this->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

				// Set initial design value, first, we need to know which variable group it belongs to
				DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd_eStart);
				elementData = (EDGEDATA)(arr + offsetd_eStart);
				this->elementFunction.set_alphaP_initial(alphaMaxInit[elementData->VariableGroup]);

				// Set the value of alphaMax. We need to do this to update
				this->elementFunction.set_alphaMax(alphaMaxarr[k]);


				this->elementFunction.calculate_delta();



#ifdef DEBUG
				PetscReal delta_adjoint = elementFunction.get_delta();
				std::cout<<"alphaMax en adjoint = "<<alphaMaxarr[k]<<" element = "<<k<<std::endl;
				std::cout<<"delta_adjoint en adjoint = "<<delta_adjoint<<" element = "<<k<<std::endl;
#endif
				/*
				 * Now we need to build the vector that will apply to the linear operator
				 * We grab only the values related to the variable U, nothing about V
				 */
				omega_local_U(0) = omegaarr[offsetfrom];
				omega_local_U(1) = omegaarr[offsetfrom + 1];
				omega_local_U(2) = omegaarr[offsetto];
				omega_local_U(3) = omegaarr[offsetto + 1];

				omega_local_V(0) = omegaarr[offsetfrom + 2];
				omega_local_V(1) = omegaarr[offsetfrom + 3];
				omega_local_V(2) = omegaarr[offsetto + 2];
				omega_local_V(3) = omegaarr[offsetto + 3];

				this->elementFunction.jacobian_vector_product(omega_local_U,phi_local);

				this->elementFunction.jacobian_alphaMax_vector_product(omega_local_U,result);

#ifdef DEBUG
//				std::cout<<"result for tmpalpha element "<<k<<" = "<<result<<std::endl;
#endif

				tmpalphamaxarr[k] -= result;
#ifdef DEBUG
//				std::cout<<"producto del jacobiano"<<std::endl;
//				phi_local.print();
#endif
				/*
				 * U variables
				 */
				tmparr[offsetfrom + 2] 	-= phi_local(0);
				tmparr[offsetfrom + 3] 	-= phi_local(1);
				tmparr[offsetto+ 2] 	-= phi_local(2);
				tmparr[offsetto + 3] 	-= phi_local(3);

				tmparr[offsetfrom] 		= omega_local_V(0);
				tmparr[offsetfrom + 1] 	= omega_local_V(1);
				tmparr[offsetto] 			= omega_local_V(2);
				tmparr[offsetto + 1] 		= omega_local_V(3);

				if (beadsdataFrom->constrained_dof_x){
					tmparr[offsetfrom] = 0.0;
					tmparr[offsetfrom + 2] = 0.0;
				}
				if (beadsdataFrom->constrained_dof_y){
					tmparr[offsetfrom + 1] = 0.0;
					tmparr[offsetfrom + 3] = 0.0;
				}
				if (beadsdataTo->constrained_dof_x){
					tmparr[offsetto] = 0.0;
					tmparr[offsetto + 2] = 0.0;
				}
				if (beadsdataTo->constrained_dof_y){
					tmparr[offsetto + 1] = 0.0;
					tmparr[offsetto + 3] = 0.0;
				}



				// Update element counter
				k++;
			}
			/*
			 * Vectors assembling
			 */


			ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
			ierr = VecRestoreArray(localOmega,&omegaarr);CHKERRQ(ierr);
			ierr = VecRestoreArray(localTMP,&tmparr);CHKERRQ(ierr);


			ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);


//			ierr = DMLocalToGlobalBegin(networkdm,localOmega,INSERT_VALUES,Omega);CHKERRQ(ierr);
//			ierr = DMLocalToGlobalEnd(networkdm,localOmega,INSERT_VALUES,Omega);CHKERRQ(ierr);
			ierr = DMRestoreLocalVector(networkdm,&localOmega);CHKERRQ(ierr);

			ierr = DMLocalToGlobalBegin(networkdm,localTMP,ADD_VALUES,TMP);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(networkdm,localTMP,ADD_VALUES,TMP);CHKERRQ(ierr);
			ierr = DMRestoreLocalVector(networkdm,&localTMP);CHKERRQ(ierr);

			ierr = VecRestoreArray(TMPalpha,&tmpalphamaxarr);




			/*
			 * Add contribution from the objective function
			 *
			 * VecWAXPY(Vec w,PetscScalar alpha,Vec x,Vec y)
			 * w = alpha x + y.
			 *
			 * VecAXPY(Vec y,PetscScalar alpha,Vec x)
			 *  y = alpha x + y.
			 *
			 *  VecAXPBY(Vec y,PetscScalar alpha,PetscScalar beta,Vec x)
			 *   y = alpha x + beta y.
			 */

			VecWAXPY(Phi[i],alpha,dFdU,TMP);

			VecAXPBY(V_i[i],alpha,1.0,dFdP);

			VecWAXPY(W_i[i],alpha,dFdalphaMax,TMPalpha);

//			std::cout<<"W_i["<<i<<"]"<<std::endl;
//			VecView(W_i[i],PETSC_VIEWER_STDOUT_WORLD);
		}

#ifdef DEBUG
//		std::cout<<"Phi[i]"<<std::endl;
//		VecView(Phi[0],PETSC_VIEWER_STDOUT_WORLD);
//
//		std::cout<<"W_i[i]"<<std::endl;
//		VecView(W_i[0],PETSC_VIEWER_STDOUT_WORLD);
#endif

		ierr =  VecRestoreArrayRead(*alphaMaxVec,&alphaMaxarr);CHKERRQ(ierr);


#ifdef DEBUG
		std::cout<<"Gamma before alpha"<<std::endl;
		VecView(Gamma,PETSC_VIEWER_STDOUT_WORLD);
#endif

		/*
		 * Now we need to calculate dHdalphaMax
		 */

		ierr = DMGetLocalVector(networkdm,&localLambda);CHKERRQ(ierr);
		ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);

		ierr = DMGlobalToLocalBegin(networkdm,*Implicit_Variable_post,INSERT_VALUES,localU);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(networkdm,*Implicit_Variable_post,INSERT_VALUES,localU);CHKERRQ(ierr);

		ierr = VecGetArrayRead(localU,&uarr);
		ierr = VecGetArray(Gamma,&gammaPrevarr);

		ierr =  VecGetArrayRead(*alphaMaxVec,&alphaMaxarr);CHKERRQ(ierr);

		// Index for alphaMaxarr
		k = 0;
		// Get the bounds [start, end) for the edges.
		ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
		for (e = eStart; e< eEnd; e++ ){

			PetscInt    offsetd,key;

			const PetscInt *cone;
			DMNetworkGetConnectedNodes(networkdm,e,&cone);
			vfrom = cone[0];
			vto   = cone[1];


			ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
			// Cast component data
			beadsdataFrom = (arrayBeads)(arr + offsetd);

			ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
			// Cast component data
			beadsdataTo = (arrayBeads)(arr + offsetd);

			ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
			ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

			// Grab displacement components
			U_i_from = uarr[offsetfrom + 2];
			U_j_from = uarr[offsetfrom + 3];

			U_i_to = uarr[offsetto + 2];
			U_j_to = uarr[offsetto + 3];

			// Reset the beads for the constitutive law object
			this->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, this->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

			// Set initial design value, first, we need to know which variable group it belongs to
			DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd_eStart);
			elementData = (EDGEDATA)(arr + offsetd_eStart);
			this->elementFunction.set_alphaP_initial(alphaMaxInit[elementData->VariableGroup]);

			// Set the value of alphaMax. We need to do this to update
			this->elementFunction.set_alphaMax(alphaMaxarr[k]);

			this->elementFunction.calculate_delta();

#ifdef DEBUG
			std::cout<<"alphaMax en dHdalphaMax = "<<alphaMaxarr[k]<<" element = "<<k<<std::endl;
			std::cout<<"delta_adjoint en dHdalphaMax = "<<this->elementFunction.get_delta()<<" element = "<<k<<std::endl;
#endif
			this->elementFunction.state_eq_product_alphaMax(gammaPrevarr[k],result);
#ifdef DEBUG
//			std::cout<<"result for gamma = "<<result<<std::endl;
#endif
			gammaPrevarr[k] = result;


			// Update element counter
			k++;
		}

		ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
		ierr = VecRestoreArray(Gamma,&gammaPrevarr);CHKERRQ(ierr);

#ifdef DEBUG
		std::cout<<"Gamma after alpha"<<std::endl;
		VecView(Gamma,PETSC_VIEWER_STDOUT_WORLD);
#endif

		ierr = DMRestoreLocalVector(networkdm,&localLambda);CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);

		ierr =  VecRestoreArrayRead(*alphaMaxVec,&alphaMaxarr);CHKERRQ(ierr);


		/*
		 * Update variables
		 * VecMAXPY(Vec y,PetscInt nv,const PetscScalar alpha[],Vec x[])
		 * y = y + sum alpha[j] x[j]
		 */

		VecMAXPY(Lambda,Stages,ones,Phi);CHKERRQ(ierr);
		VecMAXPY(Mu,Stages,ones,V_i);CHKERRQ(ierr);

#ifdef DEBUG
		std::cout<<"Mu"<<std::endl;
		VecView(Mu,PETSC_VIEWER_STDOUT_WORLD);
#endif

//		std::cout<<"GammaPrev"<<std::endl;
//		VecView(GammaPrev,PETSC_VIEWER_STDOUT_WORLD);

		VecMAXPY(Gamma,Stages,ones,W_i);CHKERRQ(ierr);
		//VecCopy(GammaPrev,Gamma);


//
//		std::cout<<"Lambda antes de gamma"<<std::endl;
//		VecView(Lambda,PETSC_VIEWER_STDOUT_WORLD);

		/*
		 * Add the sum to Lambda
		 */

//#ifdef DEBUG
//		std::cout<<"Gamma before lambda"<<std::endl;
//		VecView(Gamma,PETSC_VIEWER_STDOUT_WORLD);

//		std::cout<<"Lambda before gamma"<<std::endl;
//		VecView(Lambda,PETSC_VIEWER_STDOUT_WORLD);
//#endif


		if (alphaMaxVec_previous != alphaMaxVec_end) {

			ierr = DMGetLocalVector(networkdm,&localLambda);CHKERRQ(ierr);
			ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);

			ierr = DMGlobalToLocalBegin(networkdm,Lambda,INSERT_VALUES,localLambda);CHKERRQ(ierr);
			ierr = DMGlobalToLocalEnd(networkdm,Lambda,INSERT_VALUES,localLambda);CHKERRQ(ierr);

			ierr = DMGlobalToLocalBegin(networkdm,*Implicit_Variable,INSERT_VALUES,localU);CHKERRQ(ierr);
			ierr = DMGlobalToLocalEnd(networkdm,*Implicit_Variable,INSERT_VALUES,localU);CHKERRQ(ierr);

			ierr = VecGetArray(localLambda,&lambdaarr);
			ierr = VecGetArrayRead(localU,&uarr);
			ierr = VecGetArray(Gamma,&gammaPrevarr);

			ierr =  VecGetArrayRead(*alphaMaxVec_previous,&alphaMaxarr);CHKERRQ(ierr);

			// Index for alphaMaxarr
			k = 0;
			// Get the bounds [start, end) for the edges.
			ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
			for (e = eStart; e< eEnd; e++ ){

				PetscInt    offsetd,key;

				const PetscInt *cone;
				DMNetworkGetConnectedNodes(networkdm,e,&cone);
				vfrom = cone[0];
				vto   = cone[1];


				ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataFrom = (arrayBeads)(arr + offsetd);

				ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataTo = (arrayBeads)(arr + offsetd);

				ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
				ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

				// Grab displacement components
				U_i_from = uarr[offsetfrom + 2];
				U_j_from = uarr[offsetfrom + 3];

				U_i_to = uarr[offsetto + 2];
				U_j_to = uarr[offsetto + 3];


				// Reset the beads for the constitutive law object
				this->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, this->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

				// Set initial design value, first, we need to know which variable group it belongs to
				DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd_eStart);
				elementData = (EDGEDATA)(arr + offsetd_eStart);
				this->elementFunction.set_alphaP_initial(alphaMaxInit[elementData->VariableGroup]);

				// Set the value of alphaMax. We need to do this to update
				this->elementFunction.set_alphaMax(alphaMaxarr[k]);

				this->elementFunction.calculate_delta();



#ifdef DEBUG
				PetscReal delta_adjoint = elementFunction.get_delta();
				std::cout<<"alphaMax en adjoint PREVIOUS = "<<alphaMaxarr[k]<<" element = "<<k<<std::endl;
				std::cout<<"delta_adjoint en adjoint PREVIOUS = "<<delta_adjoint<<" element = "<<k<<std::endl;
#endif
				this->elementFunction.state_eq_vector_product_U(gammaPrevarr[k],phi_local);

				lambdaarr[offsetfrom + 2] 	+= phi_local(0);
				lambdaarr[offsetfrom + 3] 	+= phi_local(1);
				lambdaarr[offsetto+ 2] 		+= phi_local(2);
				lambdaarr[offsetto + 3] 	+= phi_local(3);

#ifdef DEBUG
//				std::cout<<"phi_local"<<std::endl;
//				phi_local.print();
//
//				if (phi_local(0) != 0)
//					std::cout<<"Transition reached element = "<<k<<std::endl;
#endif



				if (beadsdataFrom->constrained_dof_x){
					lambdaarr[offsetfrom] = 0.0;
					lambdaarr[offsetfrom + 2] = 0.0;
				}
				if (beadsdataFrom->constrained_dof_y){
					lambdaarr[offsetfrom + 1] = 0.0;
					lambdaarr[offsetfrom + 3] = 0.0;
				}
				if (beadsdataTo->constrained_dof_x){
					lambdaarr[offsetto] = 0.0;
					lambdaarr[offsetto + 2] = 0.0;
				}
				if (beadsdataTo->constrained_dof_y){
					lambdaarr[offsetto + 1] = 0.0;
					lambdaarr[offsetto + 3] = 0.0;
				}


				// Update element counter
				k++;
			}

			ierr = VecRestoreArray(localLambda,&lambdaarr);CHKERRQ(ierr);
			ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
			ierr = VecRestoreArray(Gamma,&gammaPrevarr);CHKERRQ(ierr);

//			ierr = DMLocalToGlobalBegin(networkdm,Lambda,ADD_VALUES,localLambda);CHKERRQ(ierr);
//			ierr = DMLocalToGlobalEnd(networkdm,Lambda,ADD_VALUES,localLambda);CHKERRQ(ierr);

			ierr = DMLocalToGlobalBegin(networkdm,localLambda,INSERT_VALUES,Lambda);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(networkdm,localLambda,INSERT_VALUES,Lambda);CHKERRQ(ierr);

			ierr = DMRestoreLocalVector(networkdm,&localLambda);CHKERRQ(ierr);
			ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);

			ierr =  VecRestoreArrayRead(*alphaMaxVec_previous,&alphaMaxarr);CHKERRQ(ierr);

		}

#ifdef DEBUG
		std::cout<<"Lambda after gamma"<<std::endl;
		VecView(Lambda,PETSC_VIEWER_STDOUT_WORLD);
#endif




	}




	end = clock();

	std::cout<<"NUmero de pasos dados en el adjoint = "<<pasos<<std::endl;

	totalcputime = ((float)(end - start)) / CLOCKS_PER_SEC ;

	std::cout<<"Adjoint Problem time = "<<totalcputime<<std::endl;

	/*
	 * Grab gradient, it is equal to the last value of Mu
	 */

#ifdef PETSC_USE_LOG
	  PetscLogStagePop();
#endif
	return ierr;
}

PetscErrorCode System::ProcessGradient() {

	PetscErrorCode ierr;

	_is_gradient_calculated = true;


	VecView(Mu,PETSC_VIEWER_STDOUT_WORLD);
	VecView(Gamma,PETSC_VIEWER_STDOUT_WORLD);
	VecAXPY(Mu,1.0,Gamma);

	return ierr;

}

PetscErrorCode System::PerturbPETScVector(Vec U, PetscInt & Index, PetscReal  & perturb){
	PetscErrorCode ierr;
	PetscScalar * uarr;
	/*
	 * Steps to modify a PETSc vector within DMNetwork
	 *
  	 * 	ierr = VecGetArray(U,&uarr);CHKERRQ(ierr);
  	 *
  	 * 	Modify uarr
  	 *
  	 * 	ierr = VecRestoreArray(U,&uarr);CHKERRQ(ierr);
  	 *
  	 * 	ierr = VecAssemblyBegin(U);
  	 *	ierr = VecAssemblyEnd(U);
	 */
	ierr = VecGetArray(U,&uarr);

	uarr[Index] += perturb;

	ierr = VecRestoreArray(U,&uarr);


	return ierr;

}


PetscErrorCode System::ObjFunctionTime(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){
	  PetscErrorCode ierr;
	  ierr = (*this.*CalculateObjFunctionTime)(U,alphaMax,NumObjFunct);
	  return ierr;
}

PetscErrorCode System::CalculateObjFunctionTime_P_Norm(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){

	  PetscErrorCode ierr;

	  PetscInt      eStart,eEnd,vfrom,vto;
	  const PetscScalar *uarr;
	  PetscScalar alphaMaxarr[1];
	  PetscScalar U_i_from, U_j_from, U_i_to, U_j_to;
	  PetscInt      offsetfrom,offsetto;
	  DMNetworkComponentGenericDataType *arr;
	  arrayBeads beadsdataFrom, beadsdataTo;


	ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	// Get the bounds [start, end) for the edges.
	ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

	FunctionValue = 0.0;

	PetscReal LocalFunctionValue = 0.0;
	PetscReal FValueTotalElement = 0.0;

	// Initialize derivative coefficients
	theta = 0.0;

	// To make the code cleaner
	unsigned int TimeNorm = designData.TimeNorm;
	unsigned int SpaceNorm = designData.SpaceNorm;

	Vec           localU;
	ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
	// Get displacement vector for this time step
	ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	// Copy local vector to an array
	ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);

	VecAssemblyBegin(alphaMax);
	VecAssemblyEnd(alphaMax);

	std::vector<int>::iterator 		target 	= designData.TargetArea_list[NumObjFunct].begin();
	const std::vector<int>::iterator 	end 	= designData.TargetArea_list[NumObjFunct].end();

	  /* Edge Data necessary to know which variable group it belongs to */
	  EDGEDATA elementData;
	  PetscInt key, offsetd_eStart;

	for (; target != end; ++target){
		// Only perform calculations if the element is within this processor


		if (*target - 1 >= eStart && *target - 1 <eEnd){




//						std::cout<<"Element obj function dentro = "<<*target<<std::endl;
					// Recover information for alphaMax
					// We need to call the values of alphaMax for this element, because we have passed
					// the conditional, we know the element is within this processor
					PetscInt elementIndex = *target - 1;

					VecGetValues(alphaMax,1,&elementIndex,alphaMaxarr);
					PetscInt    offsetd,key;


					// Beads that belong to this element/contact
					const PetscInt *cone;
					DMNetworkGetConnectedNodes(networkdm,elementIndex,&cone);
					vfrom = cone[0];
					vto   = cone[1];

					/*
					 * Get the data from these beads
					 */
					ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
					// Cast component data
					beadsdataFrom = (arrayBeads)(arr + offsetd);

					ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
					// Cast component data
					beadsdataTo = (arrayBeads)(arr + offsetd);

//						std::cout<<"coordenadas de los nodos del elemento"<<std::endl;

					// We need to correct the id to access the information in the DMNetwork

//						PetscInt bead_id_from = vfrom + 1 - problemData.N_Elements;
//						PetscInt bead_id_to = vto + 1 - problemData.N_Elements;

//						std::cout<<"bead_id_from = "<<bead_id_from<<std::endl;
//						std::cout<<"bead_id_to = "<<bead_id_to<<std::endl;
//						std::cout<<"x from = "<<beadsdataFrom->coordx<<std::endl;
//						std::cout<<"y from = "<<beadsdataFrom->coordy<<std::endl;
//
//						std::cout<<"x to = "<<beadsdataTo->coordx<<std::endl;
//						std::cout<<"y to = "<<beadsdataTo->coordy<<std::endl;

					ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
					ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

					// Grab displacement components
					U_i_from = uarr[offsetfrom + 2];
					U_j_from = uarr[offsetfrom + 3];

					U_i_to = uarr[offsetto + 2];
					U_j_to = uarr[offsetto + 3];

					// Reset the beads for the constitutive law object
					this->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, this->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

					// Set initial design value, first, we need to know which variable group it belongs to
					DMNetworkGetComponentTypeOffset(networkdm,elementIndex,0,&key,&offsetd_eStart);
					elementData = (EDGEDATA)(arr + offsetd_eStart);
					this->elementFunction.set_alphaP_initial(alphaMaxInit[elementData->VariableGroup]);


					this->elementFunction.set_alphaMax(*alphaMaxarr);
					this->elementFunction.calculate_delta();
					PetscReal F_element_element = this->elementFunction.force_coefficient();

					// Add contribution for this time step
					FValueTotalElement += PetscPowScalar(F_element_element,SpaceNorm);

//						std::cout<<"F_element_element for element "<<*target<<" = "<<F_element_element<<std::endl;
//						std::cout<<"FValueTotalElement  "<<*target<<" = "<<FValueTotalElement<<std::endl;


				}
	}
	LocalFunctionValue += PetscPowScalar(FValueTotalElement, TimeNorm/SpaceNorm);

//		std::cout<<"LocalFunctionValue = "<<LocalFunctionValue<<std::endl;

	// Gather values for the objective functions
	MPI_Allreduce(&LocalFunctionValue, &FunctionValue, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);

//	std::cout<<"FunctionValue = "<<FunctionValue<<std::endl;

	ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);

	return ierr;

}

PetscErrorCode System::CalculateObjFunctionTime_Simple(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){

	  PetscErrorCode ierr;

	  PetscInt      eStart,eEnd,vfrom;
	  const PetscScalar *uarr;
	  PetscScalar alphaMaxarr[1];
	  PetscScalar U_i_from;
	  PetscInt      offsetfrom;


	// Get the bounds [start, end) for the edges.
	ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

	FunctionValue = 0.0;

	PetscReal LocalFunctionValue = 0.0;
	PetscReal FValueTotalElement = 0.0;

	// Initialize derivative coefficients
	theta = 0.0;


	Vec           localU;
	ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
	// Get displacement vector for this time step
	ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	// Copy local vector to an array
	ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);

	VecAssemblyBegin(alphaMax);
	VecAssemblyEnd(alphaMax);


	std::vector<int>::iterator 		target 	= designData.TargetArea_list[NumObjFunct].begin();
	const std::vector<int>::iterator 	end 	= designData.TargetArea_list[NumObjFunct].end();

	for (; target != end; ++target){
		// Only perform calculations if the element is within this processor
		if (*target - 1 >= eStart && *target - 1 <eEnd){





					// Recover information for alphaMax
					// We need to call the values of alphaMax for this element, because we have passed
					// the conditional, we know the element is within this processor
					PetscInt elementIndex = *target - 1;

					VecGetValues(alphaMax,1,&elementIndex,alphaMaxarr);

					// Beads that belong to this element/contact
					const PetscInt *cone;
					DMNetworkGetConnectedNodes(networkdm,elementIndex,&cone);
					vfrom = cone[0];


					ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);


					// Grab displacement components
					U_i_from = uarr[offsetfrom + 2];


					// Add contribution for this time step
					FValueTotalElement += U_i_from;

				}
	}
	LocalFunctionValue += FValueTotalElement;

	// Gather values for the objective functions
	MPI_Allreduce(&LocalFunctionValue, &FunctionValue, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);

	ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);

	return ierr;
}
PetscErrorCode System::ProcessObjFunctionTime(){

	PetscErrorCode ierr;
	/*
	 * Grab the objective function from the solution vector in the last time step
	 */
	PetscReal FunctionValueFinalStep;
	Vec localU;
	const PetscReal * uarr;
	PetscInt eStart, eEnd, offsetfrom;

	ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localU);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localU);CHKERRQ(ierr);

	ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);

	/*
	* Initial value of the objective function
	*/
	ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
	if (eStart == 0){
	  ierr = DMNetworkGetVariableOffset(networkdm,eStart,&offsetfrom);CHKERRQ(ierr);

	  for (PetscInt j = 0; j < designData.N_ObjFunc; j++)
		  CapitalOmega[j] = uarr[offsetfrom + j];
	}

	ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);

	MPI_Bcast(&CapitalOmega.front(),1, MPI_REAL, 0,PETSC_COMM_WORLD);

	//FunctionValueBeforePnorm = FunctionValueFinalStep;

	if (designData.objfunctionType == 1){

		for (PetscInt j = 0; j < designData.N_ObjFunc; j++)
			FunctionObjetivoPartial[j] = PetscPowScalar(CapitalOmega[j],1.0/designData.TimeNorm);

	}
	else
		FunctionValueGlobal = FunctionValueFinalStep;

	if (rank == 0)
		//std::cout<<"FunctionValue = "<<FunctionValueGlobal<<std::endl;



	FunctionValueGlobal = 0.0;
	for (PetscInt j = 0; j < designData.N_ObjFunc; j++){
		if (designData.MaxOrMin[j])
			FunctionValueGlobal += -1.0*FunctionObjetivoPartial[j];
		else
			FunctionValueGlobal += FunctionObjetivoPartial[j];
	}


	/*
	 * Calculate the CapitalOmega coefficient that will be used
	 * in the adjoint analysis
	 */
	for (PetscInt j = 0; j < designData.N_ObjFunc; j++){
		if (designData.MaxOrMin[j]){
			if (CapitalOmega[j] != 0)
				CapitalOmega[j] = -1.0/designData.TimeNorm*PetscPowScalar(CapitalOmega[j],1.0/designData.TimeNorm - 1.0);
			else
				CapitalOmega[j] = 0;
		}
		else{
			if (CapitalOmega[j] != 0)
				CapitalOmega[j] = 1.0/designData.TimeNorm*PetscPowScalar(CapitalOmega[j],1.0/designData.TimeNorm - 1.0);
			else
				CapitalOmega[j] = 0;
		}
	}


	FunctionValueGlobalBeforePenalty = FunctionValueGlobal;


	std::cout<<"FunctionValueGlobal = "<<std::setprecision(20)<<FunctionValueGlobal<<std::endl;

	return ierr;
}


PetscErrorCode System::CheckGradient(){

	assert(_is_fd_calculated && "Need to calculate the FD gradient first");

	assert(_is_gradient_calculated && "Need to calculate the gradient first");

	PetscErrorCode ierr;
	/*
	 * Range of the local beads
	 */
	PetscInt vStart, vEnd;
	ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);



	const PetscScalar * muarr;

	ierr = VecGetArrayRead(Mu,&muarr);CHKERRQ(ierr);

	std::cout<<"DO NOT RUN THIS IN PARALLEL"<<std::endl;
	  for (unsigned int i = 0; i<designData.N_DesignVariables; i++){
		  std::cout<<"Gradient FD ["<<i<<"] = "<<std::setprecision(10)<<Gradient_FD[i]<<" Gradient ["<<i<<"] = "<<std::setprecision(10)<<muarr[i]<<std::endl;
	  }

	ierr =  VecRestoreArrayRead(Mu,&muarr);CHKERRQ(ierr);


	return ierr;

}

PetscErrorCode System::EvaluatePartialImplicitDerivatives(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){
	PetscErrorCode ierr;

	ierr = (*this.*PartialImplicitDerivatives)(U,alphaMax,NumObjFunct);

	return ierr;
}
PetscErrorCode System::EvaluatePartialExplicitDerivatives(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){
	PetscErrorCode ierr;

	ierr = (*this.*PartialExplicitDerivatives)(U,alphaMax,NumObjFunct);

	return ierr;
}

PetscErrorCode System::EvaluatePartialImplicitDerivatives_P_Norm(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){



	_is_partial_derivatives_calculated = true;
	PetscErrorCode ierr;

	PetscInt      eStart,eEnd,vfrom,vto;
	const PetscScalar *uarr;
	PetscScalar   *partialF_partialU_arr;
	PetscScalar alphaMaxarr[1];
	PetscScalar U_i_from, U_j_from, U_i_to, U_j_to;
	PetscInt      offsetfrom,offsetto;
	DMNetworkComponentGenericDataType *arr;
	arrayBeads beadsdataFrom, beadsdataTo;



	ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	// Get the bounds [start, end) for the edges.
	ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

	/*
	 * Zero out the vector
	 */
	VecSet(dFdU_partial,0.0);
	VecSet(dFdalphaMax_partial,0.0);


	/*
	 * Coefficients beta and gamma
	 * because of the chain rule
	 */
	PetscReal beta_coefficient_local 	= 0.0;
	PetscReal beta_coefficient_partial  = 0.0;
	PetscReal beta_coefficient  		= 0.0;

	unsigned int TimeNorm = designData.TimeNorm;
	unsigned int SpaceNorm = designData.SpaceNorm;

	Vec           localU, localdFdU;
	ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
	ierr = DMGetLocalVector(networkdm,&localdFdU);CHKERRQ(ierr);
	// Get displacement vector for this time step
	ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(networkdm,dFdU_partial,INSERT_VALUES,localdFdU);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(networkdm,dFdU_partial,INSERT_VALUES,localdFdU);CHKERRQ(ierr);

	// Copy local vector to an array
	ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);

	ierr = VecGetArray(localdFdU,&partialF_partialU_arr);CHKERRQ(ierr);

	VecAssemblyBegin(alphaMax);
	VecAssemblyEnd(alphaMax);

	std::vector<int>::iterator 		target 	= designData.TargetArea_list[NumObjFunct].begin();
	const std::vector<int>::iterator 	end 	= designData.TargetArea_list[NumObjFunct].end();

	// Derivative coefficient
	PetscScalar deriv_coef = 0.0;

	  /* Edge Data necessary to know which variable group it belongs to */
	  EDGEDATA elementData;
	  PetscInt key, offsetd_eStart;





	for (; target != end; ++target){
		// Only perform calculations if the element is within this processor
		if (*target - 1 >= eStart && *target - 1 <eEnd){




				// Recover information for alphaMax
				// We need to call the values of alphaMax for this element, because we have passed
				// the conditional, we know the element is within this processor
				PetscInt elementIndex = *target - 1;

				VecGetValues(alphaMax,1,&elementIndex,alphaMaxarr);




				// Beads that belong to this element/contact
				const PetscInt *cone;
				DMNetworkGetConnectedNodes(networkdm,elementIndex,&cone);
				vfrom = cone[0];
				vto   = cone[1];

				/*
				 * Get the data from these beads
				 */
				PetscInt    offsetd,key;
				ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataFrom = (arrayBeads)(arr + offsetd);

				ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataTo = (arrayBeads)(arr + offsetd);


				ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
				ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

				// Grab displacement components
				U_i_from = uarr[offsetfrom + 2];
				U_j_from = uarr[offsetfrom + 3];

				U_i_to = uarr[offsetto + 2];
				U_j_to = uarr[offsetto + 3];

				// Reset the beads for the constitutive law object
				this->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, this->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

				// Set initial design value, first, we need to know which variable group it belongs to
				DMNetworkGetComponentTypeOffset(networkdm,elementIndex,0,&key,&offsetd_eStart);
				elementData = (EDGEDATA)(arr + offsetd_eStart);
				this->elementFunction.set_alphaP_initial(alphaMaxInit[elementData->VariableGroup]);


				this->elementFunction.set_alphaMax(*alphaMaxarr);
				this->elementFunction.calculate_delta();

				PetscReal F_element_time = this->elementFunction.force_coefficient();

				// Calculating the partial derivatives coefficients


				/*
				 * Partial derivative w.r.t vector of variables
				 */
				PetscReal *  partial_F_element_time = this->elementFunction.partial_derivative_U_force_vector();

				// Derivative coefficient, from the chain rule in the objective function
				//deriv_coef = theta*beta[index_target_element]*gamma;
				deriv_coef = SpaceNorm * PetscPowScalar(F_element_time, SpaceNorm - 1.0);
				beta_coefficient_local += PetscPowScalar(F_element_time, SpaceNorm);

				partialF_partialU_arr[offsetfrom + 2] += deriv_coef*partial_F_element_time[0];
				partialF_partialU_arr[offsetfrom + 3] += deriv_coef*partial_F_element_time[1];

				partialF_partialU_arr[offsetfrom] = 0.0;
				partialF_partialU_arr[offsetfrom + 1] = 0.0;

				partialF_partialU_arr[offsetto + 2] += deriv_coef*partial_F_element_time[2];
				partialF_partialU_arr[offsetto + 3] += deriv_coef*partial_F_element_time[3];

				partialF_partialU_arr[offsetto] = 0.0;
				partialF_partialU_arr[offsetto + 1] = 0.0;

				if (beadsdataFrom->constrained_dof_x){
					partialF_partialU_arr[offsetfrom] = 0.0;
					partialF_partialU_arr[offsetfrom + 2] = 0.0;
				}
				if (beadsdataFrom->constrained_dof_y){
					partialF_partialU_arr[offsetfrom + 1] = 0.0;
					partialF_partialU_arr[offsetfrom + 3] = 0.0;
				}
				if (beadsdataTo->constrained_dof_x){
					partialF_partialU_arr[offsetto] = 0.0;
					partialF_partialU_arr[offsetto + 2] = 0.0;
				}
				if (beadsdataTo->constrained_dof_y){
					partialF_partialU_arr[offsetto + 1] = 0.0;
					partialF_partialU_arr[offsetto + 3] = 0.0;
				}





				/*
				 * Partial derivative w.r.t alphaMax
				 */

				PetscReal partialF_partialAlphaMax = this->elementFunction.partial_derivative_alphaMax_force_vector();

				partialF_partialAlphaMax = partialF_partialAlphaMax*deriv_coef;

				// Assemble into the corresponding time for partialF_partialU
				VecSetValue(dFdalphaMax_partial,elementIndex,(PetscScalar)(partialF_partialAlphaMax),ADD_VALUES);

		}
	}
	MPI_Allreduce(&beta_coefficient_local, &beta_coefficient_partial, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
	beta_coefficient =  TimeNorm/SpaceNorm * PetscPowScalar(beta_coefficient_partial, TimeNorm/SpaceNorm - 1.0);

	// Restore the modified vectors
//	ierr = VecRestoreArray(dFdU,&partialF_partialU_arr);
//	ierr = VecAssemblyBegin(dFdU);
//	ierr = VecAssemblyEnd(dFdU);

	ierr = VecAssemblyBegin(dFdalphaMax_partial);
	ierr = VecAssemblyEnd(dFdalphaMax_partial);

	ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);


	ierr = VecRestoreArray(localdFdU,&partialF_partialU_arr);
	ierr = DMLocalToGlobalBegin(networkdm,localdFdU,ADD_VALUES,dFdU_partial);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(networkdm,localdFdU,ADD_VALUES,dFdU_partial);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(networkdm,&localdFdU);CHKERRQ(ierr);

	VecScale(dFdU_partial,beta_coefficient);
//	std::cout<<"Result of the interpolated dFdU"<<std::endl;
//	VecView(dFdU_interp,PETSC_VIEWER_STDOUT_WORLD);

	VecScale(dFdalphaMax_partial,beta_coefficient);
//	std::cout<<"Result of the interpolated dFdalphaMax"<<std::endl;
//	VecView(dFdalphaMax_interp,PETSC_VIEWER_STDOUT_WORLD);


	return ierr;

}

PetscErrorCode System::EvaluatePartialExplicitDerivatives_P_Norm(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){


	_is_partial_derivatives_calculated = true;
	PetscErrorCode ierr;

	PetscInt      eStart,eEnd,vfrom,vto;
	const PetscScalar *uarr;
	PetscScalar alphaMaxarr[1];
	PetscScalar U_i_from, U_j_from, U_i_to, U_j_to;
	PetscInt      offsetfrom,offsetto;
	DMNetworkComponentGenericDataType *arr;
	arrayBeads beadsdataFrom, beadsdataTo;


	/*
	 * Zero out the vector
	 */
	VecSet(dFdP_partial,0.0);


	ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	// Get the bounds [start, end) for the edges.
	ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

	/*
	 * Coefficients beta and gamma
	 * because of the chain rule
	 */
	PetscReal beta_coefficient_local 	= 0.0;
	PetscReal beta_coefficient_partial  = 0.0;
	PetscReal beta_coefficient  		= 0.0;

	unsigned int TimeNorm = designData.TimeNorm;
	unsigned int SpaceNorm = designData.SpaceNorm;


	std::vector<int>::iterator 		target 	= designData.TargetArea_list[NumObjFunct].begin();
	const std::vector<int>::iterator 	end 	= designData.TargetArea_list[NumObjFunct].end();

	// Derivative coefficient
	PetscScalar deriv_coef = 0.0;

	Vec           localU;
	ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
	// Get displacement vector for this time step
	ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	// Copy local vector to an array
	ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);

	  /* Edge Data necessary to know which variable group it belongs to */
	  EDGEDATA elementData;
	  PetscInt key, offsetd_eStart;

	for (; target != end; ++target){
		// Only perform calculations if the element is within this processor
		if (*target - 1 >= eStart && *target - 1 <eEnd){

				// Recover information for alphaMax
				// We need to call the values of alphaMax for this element, because we have passed
				// the conditional, we know the element is within this processor
				PetscInt elementIndex = *target - 1;
				VecAssemblyBegin(alphaMax);
				VecAssemblyEnd(alphaMax);
				VecGetValues(alphaMax,1,&elementIndex,alphaMaxarr);




				// Beads that belong to this element/contact
				const PetscInt *cone;
				DMNetworkGetConnectedNodes(networkdm,elementIndex,&cone);
				vfrom = cone[0];
				vto   = cone[1];

				/*
				 * Get the data from these beads
				 */
				PetscInt    offsetd,key;
				ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataFrom = (arrayBeads)(arr + offsetd);

				ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataTo = (arrayBeads)(arr + offsetd);


				ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
				ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

				// Grab displacement components
				U_i_from = uarr[offsetfrom + 2];
				U_j_from = uarr[offsetfrom + 3];

				U_i_to = uarr[offsetto + 2];
				U_j_to = uarr[offsetto + 3];

				// Reset the beads for the constitutive law object
				this->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, this->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

				// Set initial design value, first, we need to know which variable group it belongs to
				DMNetworkGetComponentTypeOffset(networkdm,elementIndex,0,&key,&offsetd_eStart);
				elementData = (EDGEDATA)(arr + offsetd_eStart);
				this->elementFunction.set_alphaP_initial(alphaMaxInit[elementData->VariableGroup]);

				this->elementFunction.set_alphaMax(*alphaMaxarr);
				this->elementFunction.calculate_delta();

				PetscReal F_element_time = this->elementFunction.force_coefficient();

				// Derivative coefficient, from the chain rule in the objective function
				deriv_coef = SpaceNorm * PetscPowScalar(F_element_time, SpaceNorm - 1.0);
				beta_coefficient_local += PetscPowScalar(F_element_time, SpaceNorm);

				/*
				 * Partial derivative w.r.t parameter
				 */

				PetscReal & partialF_partialP_element = this->elementFunction.partial_derivative_P();
				PetscScalar value = deriv_coef*partialF_partialP_element;
				VecSetValue(dFdP_partial,elementData->VariableGroup,value,ADD_VALUES);

		}

	}
	VecAssemblyBegin(dFdP_partial);
	VecAssemblyEnd(dFdP_partial);

	ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);
	MPI_Allreduce(&beta_coefficient_local, &beta_coefficient_partial, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
	beta_coefficient =  TimeNorm/SpaceNorm * PetscPowScalar(beta_coefficient_partial, TimeNorm/SpaceNorm - 1.0);


	VecScale(dFdP_partial,beta_coefficient);
	return ierr;


}

PetscErrorCode System::EvaluatePartialImplicitDerivatives_Simple(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){


	_is_partial_derivatives_calculated = true;
	PetscErrorCode ierr;

	PetscInt      eStart,eEnd,vfrom;
	const PetscScalar *uarr;
	PetscScalar   *partialF_partialU_arr;
	PetscScalar alphaMaxarr[1];
	PetscInt      offsetfrom;
	DMNetworkComponentGenericDataType *arr;
	arrayBeads beadsdataFrom;



	ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	// Get the bounds [start, end) for the edges.
	ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

	/*
	 * Zero out the vector
	 */
	VecSet(dFdU_partial,0.0);
	VecSet(dFdalphaMax_partial,0.0);


	Vec           localU;
	ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
	// Get displacement vector for this time step
	ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	// Copy local vector to an array
	ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);


	std::vector<int>::iterator 		target 	= designData.TargetArea_list[NumObjFunct].begin();
	const std::vector<int>::iterator 	end 	= designData.TargetArea_list[NumObjFunct].end();




	for (; target != end; ++target){
		// Only perform calculations if the element is within this processor
		if (*target - 1 >= eStart && *target - 1 <eEnd){




				// Recover information for alphaMax
				// We need to call the values of alphaMax for this element, because we have passed
				// the conditional, we know the element is within this processor
				PetscInt elementIndex = *target - 1;
				VecAssemblyBegin(alphaMax);
				VecAssemblyEnd(alphaMax);
				VecGetValues(alphaMax,1,&elementIndex,alphaMaxarr);




				// Beads that belong to this element/contact
				const PetscInt *cone;
				DMNetworkGetConnectedNodes(networkdm,elementIndex,&cone);
				vfrom = cone[0];

				/*
				 * Get the data from these beads
				 */
				PetscInt    offsetd,key;
				ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
				// Cast component data
				beadsdataFrom = (arrayBeads)(arr + offsetd);




				ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);



				// Assemble into the corresponding time for partialF_partialU
				VecGetArray(dFdU_partial,&partialF_partialU_arr);

				// Derivative coefficient, from the chain rule in the objective function
				//deriv_coef = theta*beta[index_target_element]*gamma;

				partialF_partialU_arr[offsetfrom + 2] = 1.0;
				partialF_partialU_arr[offsetfrom + 3] = 0.0;

				partialF_partialU_arr[offsetfrom] = 0.0;
				partialF_partialU_arr[offsetfrom + 1] = 0.0;



				if (beadsdataFrom->constrained_dof_x){
					partialF_partialU_arr[offsetfrom] = 0.0;
					partialF_partialU_arr[offsetfrom + 2] = 0.0;
				}
				if (beadsdataFrom->constrained_dof_y){
					partialF_partialU_arr[offsetfrom + 1] = 0.0;
					partialF_partialU_arr[offsetfrom + 3] = 0.0;
				}

				// Restore the modified vectors
				ierr = VecRestoreArray(dFdU_partial,&partialF_partialU_arr);
				ierr = VecAssemblyBegin(dFdU_partial);
				ierr = VecAssemblyEnd(dFdU_partial);




		}
	}

	ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);


	return ierr;
}
PetscErrorCode System::EvaluatePartialExplicitDerivatives_Simple(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){

	VecSet(dFdP_partial,0.0);
	return(0);
}

PetscErrorCode System::PrintForcesTargetArea(Vec & U, Vec & alphaMax, PetscInt & step){
	  PetscErrorCode ierr;

		  PetscInt      eStart,eEnd,vfrom,vto;
		  const PetscScalar *uarr;
		  PetscScalar alphaMaxarr[1];
		  PetscScalar U_i_from, U_j_from, U_i_to, U_j_to;
		  PetscInt      offsetfrom,offsetto;
		  DMNetworkComponentGenericDataType *arr;
		  arrayBeads beadsdataFrom, beadsdataTo;


		ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

		// Get the bounds [start, end) for the edges.
		ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

		FunctionValue = 0.0;

		PetscReal LocalFunctionValue = 0.0;
		PetscReal FValueTotalElement = 0.0;

		// Initialize derivative coefficients
		theta = 0.0;

		// To make the code cleaner
		unsigned int TimeNorm = designData.TimeNorm;
		unsigned int SpaceNorm = designData.SpaceNorm;

		Vec           localU;
		ierr = DMGetLocalVector(networkdm,&localU);CHKERRQ(ierr);
		// Get displacement vector for this time step
		ierr = DMGlobalToLocalBegin(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(networkdm,U,INSERT_VALUES,localU);CHKERRQ(ierr);
		// Copy local vector to an array
		ierr = VecGetArrayRead(localU,&uarr);CHKERRQ(ierr);

		VecAssemblyBegin(alphaMax);
		VecAssemblyEnd(alphaMax);

		VecSet(Felement,0.0);

		  /* Edge Data necessary to know which variable group it belongs to */
		  EDGEDATA elementData;
		  PetscInt key, offsetd_eStart;



		for (PetscInt i = 0; i<designData.N_ObjFunc; i++){

			std::vector<int>::iterator 		target 	= designData.TargetArea_list[i].begin();
			const std::vector<int>::iterator 	end 	= designData.TargetArea_list[i].end();

			for (; target != end; ++target){
				// Only perform calculations if the element is within this processor


				if (*target - 1 >= eStart && *target - 1 <eEnd){




	//						std::cout<<"Element obj function dentro = "<<*target<<std::endl;
							// Recover information for alphaMax
							// We need to call the values of alphaMax for this element, because we have passed
							// the conditional, we know the element is within this processor
							PetscInt elementIndex = *target - 1;

							VecGetValues(alphaMax,1,&elementIndex,alphaMaxarr);
							PetscInt    offsetd;


							// Beads that belong to this element/contact
							const PetscInt *cone;
							DMNetworkGetConnectedNodes(networkdm,elementIndex,&cone);
							vfrom = cone[0];
							vto   = cone[1];

							/*
							 * Get the data from these beads
							 */
							ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&key,&offsetd);CHKERRQ(ierr);
							// Cast component data
							beadsdataFrom = (arrayBeads)(arr + offsetd);

							ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&key,&offsetd);CHKERRQ(ierr);
							// Cast component data
							beadsdataTo = (arrayBeads)(arr + offsetd);

	//						std::cout<<"coordenadas de los nodos del elemento"<<std::endl;

							// We need to correct the id to access the information in the DMNetwork

	//						PetscInt bead_id_from = vfrom + 1 - problemData.N_Elements;
	//						PetscInt bead_id_to = vto + 1 - problemData.N_Elements;

	//						std::cout<<"bead_id_from = "<<bead_id_from<<std::endl;
	//						std::cout<<"bead_id_to = "<<bead_id_to<<std::endl;
	//						std::cout<<"x from = "<<beadsdataFrom->coordx<<std::endl;
	//						std::cout<<"y from = "<<beadsdataFrom->coordy<<std::endl;
	//
	//						std::cout<<"x to = "<<beadsdataTo->coordx<<std::endl;
	//						std::cout<<"y to = "<<beadsdataTo->coordy<<std::endl;

							ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
							ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

							// Grab displacement components
							U_i_from = uarr[offsetfrom + 2];
							U_j_from = uarr[offsetfrom + 3];

							U_i_to = uarr[offsetto + 2];
							U_j_to = uarr[offsetto + 3];

							// Reset the beads for the constitutive law object
							this->elementFunction.reset_beads(*beadsdataFrom,*beadsdataTo, this->problemData.analysis,U_i_from, U_j_from, U_i_to, U_j_to);

							// Set initial design value, first, we need to know which variable group it belongs to
							DMNetworkGetComponentTypeOffset(networkdm,elementIndex,0,&key,&offsetd_eStart);
							elementData = (EDGEDATA)(arr + offsetd_eStart);
							this->elementFunction.set_alphaP_initial(alphaMaxInit[elementData->VariableGroup]);

							this->elementFunction.set_alphaMax(*alphaMaxarr);
							this->elementFunction.calculate_delta();
							PetscReal F_element_element = this->elementFunction.force_coefficient();


							VecSetValue(Felement,elementIndex,F_element_element,INSERT_VALUES);


							// Add contribution for this time step
							FValueTotalElement += PetscPowScalar(F_element_element,SpaceNorm);

	//						std::cout<<"F_element_element for element "<<*target<<" = "<<F_element_element<<std::endl;
	//						std::cout<<"FValueTotalElement  "<<*target<<" = "<<FValueTotalElement<<std::endl;


						}
			}
			LocalFunctionValue += PetscPowScalar(FValueTotalElement, TimeNorm/SpaceNorm);

	//		std::cout<<"LocalFunctionValue = "<<LocalFunctionValue<<std::endl;

			// Gather values for the objective functions
			MPI_Allreduce(&LocalFunctionValue, &FunctionValue, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
		}

	//	std::cout<<"FunctionValue = "<<FunctionValue<<std::endl;

		ierr = VecRestoreArrayRead(localU,&uarr);CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(networkdm,&localU);CHKERRQ(ierr);



		ierr = VecAssemblyBegin(Felement);
		ierr = VecAssemblyEnd(Felement);

		// File name
		PetscViewer    viewer;
		std::string filehistory("PostProcessing/ForcesTargetArea/");
		std::ostringstream temp;
		temp << step;
		filehistory.append(temp.str());
		filehistory.append(".mat");
		// Open and set Binary file
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,filehistory.c_str(),FILE_MODE_WRITE,&viewer);
		PetscViewerSetFormat(viewer, PETSC_VIEWER_BINARY_MATLAB);
		// Print
		VecView(Felement,viewer);
		PetscViewerDestroy(&viewer);

		return ierr;
}



PetscErrorCode System::FD_parameter(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){
	  PetscErrorCode ierr;

	  assert(_is_partial_derivatives_calculated && "Need to calculate the derivatives first");

	  // Bead pointer
	  arrayBeads DesignBead;
	  DMNetworkComponentGenericDataType *arr;

	  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	  PetscReal original_design;
	  PetscReal dt = 1e-8;

	  // Vector for the FD partial derivatives for this time step
	  Vec partialF_partialP_FD;
	  VecDuplicate(dFdP,&partialF_partialP_FD);
	  PetscScalar * dFdp_FD_arr;

	  // Zero the partialF_partialU_FD, we recycle it on each iteration
	  VecSet(partialF_partialP_FD,0.0);
	  ierr = VecGetArray(partialF_partialP_FD,&dFdp_FD_arr);

	  std::cout<<" /*----------------------------------------------------*/ "<<std::endl;
	  std::cout<<" Comparing FD of Partial Derivatives w.r.t the parameters "<<std::endl;
	  std::cout<<" /*----------------------------------------------------*/ "<<std::endl;
	  // We need to loop over the local vertices
	  PetscInt vStart, vEnd;
	  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
	  for (unsigned int i = 0; i<designData.N_DesignVariables; i++){

			  PetscReal FD_partial_deriv = 0.0;
			  PetscInt    offsetd,key;


			  // Change parameter data, first, keep original value
			  original_design = alphaMaxInit[i];

			  alphaMaxInit[i] += dt;

			  // Evaluate function
			  ObjFunctionTime(U,alphaMax,NumObjFunct);
			  FD_partial_deriv = FunctionValue;

			  // Change the parameter again, now backwards
			  alphaMaxInit[i] -= 2.0*dt;


			  // Evaluate function
			  ObjFunctionTime(U,alphaMax,NumObjFunct);
			  FD_partial_deriv -= FunctionValue;

			  // Final value
			  FD_partial_deriv *= 1.0/(2.0*dt);

			  // Recover density value
			  alphaMaxInit[i] = original_design;

			  // Insert value
			  dFdp_FD_arr[i] = FD_partial_deriv;
	  }
	  ierr = VecRestoreArray(partialF_partialP_FD,&dFdp_FD_arr);
	  ierr = VecAssemblyBegin(partialF_partialP_FD);
	  ierr = VecAssemblyEnd(partialF_partialP_FD);
	  std::cout<<"/*-----------Comparing partial derivatives -------------*/"<<std::endl;
	  std::cout<<"/*--------------------Analytical -----------------------*/"<<std::endl;
	  VecView(dFdP_partial,PETSC_VIEWER_STDOUT_WORLD);
	  std::cout<<"/*--------------------------FD -------------------------*/"<<std::endl;
	  VecView(partialF_partialP_FD,PETSC_VIEWER_STDOUT_WORLD);

	  VecDestroy(&partialF_partialP_FD);



	  return ierr;

}

PetscErrorCode System::FD_variables(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){


	  PetscErrorCode ierr;

	  assert(_is_partial_derivatives_calculated && "Need to calculate the derivatives first");


	  DMNetworkComponentGenericDataType *arr;

	  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	  PetscReal dt = 1e-10;

	  std::cout<<" /*---------------------------------------------------------------*/ "<<std::endl;
	  std::cout<<" Comparing FD of Partial Derivatives w.r.t the implicit variables "<<std::endl;
	  std::cout<<" /*---------------------------------------------------------------*/ "<<std::endl;

	  // Vector for the FD partial derivatives for this time step
	  Vec partialF_partialU_FD, Ucopy;
	  VecDuplicate(U,&partialF_partialU_FD);
	  VecDuplicate(U,&Ucopy);
	  VecCopy(U,Ucopy);
	  PetscScalar * dFdp_FD_arr;

	  // Zero the partialF_partialU_FD, we recycle it on each iteration
	  VecSet(partialF_partialU_FD,0.0);
	  ierr = VecGetArray(partialF_partialU_FD,&dFdp_FD_arr);

	  // We need to loop over the local vertices
	  PetscInt vStart, vEnd, v;
	  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
	  for (v=vStart; v < vEnd; v++) {
		  // For each time step, for each node, we perturb the four degrees of freedom
		  PetscInt    offset,key, offsetd;
			  /*
		   * Get the indices for this beads
			   */
		  ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr);
		  PetscInt index;
		  for (unsigned int j = 0; j<4 ; j++){
			  PetscReal FD_partial_deriv = 0.0;
			  PetscReal perturb;
			  index = offset + j;
			  // Perturb each index
			  perturb = dt;
			  PerturbPETScVector(Ucopy, index, perturb );

			  ObjFunctionTime(Ucopy,alphaMax,NumObjFunct);

			  FD_partial_deriv = FunctionValue;

			  // Go backwards
			  perturb = -2.0*dt;
			  PerturbPETScVector(Ucopy, index, perturb);

			  ObjFunctionTime(Ucopy,alphaMax,NumObjFunct);

			  FD_partial_deriv -= FunctionValue;

			  // Final value
			  FD_partial_deriv *= 1.0/(2.0*dt);

			  // Return vector to original value
			  perturb = dt;
			  PerturbPETScVector(Ucopy, index, perturb);

			  // Insert value
			  dFdp_FD_arr[index] = FD_partial_deriv;
		  }

		  // Apply the boundary conditions, first, get component
		  ierr = DMNetworkGetComponentTypeOffset(networkdm,v,0,&key,&offsetd);CHKERRQ(ierr);
		  arrayBeads Bead = (arrayBeads)(arr + offsetd);
		  if (Bead->constrained_dof_x){
			  dFdp_FD_arr[offset] = 0.0;
			  dFdp_FD_arr[offset + 2] = 0.0;
		  }
		  if (Bead->constrained_dof_y){
			  dFdp_FD_arr[offset + 1] = 0.0;
			  dFdp_FD_arr[offset + 3] = 0.0;
		  }

	  }
	  ierr = VecRestoreArray(partialF_partialU_FD,&dFdp_FD_arr);
	  ierr = VecAssemblyBegin(partialF_partialU_FD);
	  ierr = VecAssemblyEnd(partialF_partialU_FD);

	  // Compare the vectors
	  std::cout<<"/*-----------Comparing partial derivatives -------------*/"<<std::endl;
	  std::cout<<"/*--------------------Analytical -----------------------*/"<<std::endl;
	  VecView(dFdU_partial,PETSC_VIEWER_STDOUT_WORLD);
	  std::cout<<"/*--------------------------FD -------------------------*/"<<std::endl;
	  VecView(partialF_partialU_FD,PETSC_VIEWER_STDOUT_WORLD);

	  VecDestroy(&partialF_partialU_FD);


	  return ierr;
}

PetscErrorCode System::FD_statevariables(Vec & U, Vec & alphaMax, PetscInt & NumObjFunct){


	  PetscErrorCode ierr;

	  assert(_is_partial_derivatives_calculated && "Need to calculate the derivatives first");

	  std::vector<PetscScalar> partialF_partialP_FD;
	  partialF_partialP_FD.resize(designData.N_DesignVariables);


	  DMNetworkComponentGenericDataType *arr;

	  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	  PetscReal dt = 1e-9;

	  std::cout<<" /*---------------------------------------------------------------*/ "<<std::endl;
	  std::cout<<" Comparing FD of Partial Derivatives w.r.t the implicit variables "<<std::endl;
	  std::cout<<" /*---------------------------------------------------------------*/ "<<std::endl;


	  // Vector for the FD partial derivatives for this time step
	  Vec partialF_partialalphaMax_FD, alphaMaxCopy;
	  VecDuplicate(dFdalphaMax,&partialF_partialalphaMax_FD);
	  VecDuplicate(dFdalphaMax,&alphaMaxCopy);
	  VecCopy(alphaMax,alphaMaxCopy);

	  /*
	   * Iterate through time, i.e., displacement history and alphaMax history
	   */

		  // Zero the partialF_partialU_FD, we recycle it on each iteration
		  VecSet(partialF_partialalphaMax_FD,0.0);

		  // We need to loop over the local vertices
		  PetscInt eStart, eEnd,e;
		  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
		  for (e=eStart; e < eEnd; e++) {
			  // For each time step, for each element, we perturb the alphaMax
			  PetscReal perturb;
			  PetscReal FD_partial_deriv = 0.0;
			  perturb = dt;
			  PerturbPETScVector(alphaMaxCopy, e, perturb );

			  ObjFunctionTime(U,alphaMaxCopy,NumObjFunct);

			  FD_partial_deriv = FunctionValue;

			  // Go backwards
			  perturb = -2.0*dt;
			  PerturbPETScVector(alphaMaxCopy, e, perturb);

			  ObjFunctionTime(U,alphaMaxCopy,NumObjFunct);

			  FD_partial_deriv -= FunctionValue;

			  // Final value
			  FD_partial_deriv *= 1.0/(2.0*dt);

			  // Return vector to original value
			  perturb = dt;
			  PerturbPETScVector(alphaMaxCopy, e, perturb);

			  // Insert value
			  VecSetValue(partialF_partialalphaMax_FD,e,FD_partial_deriv,ADD_VALUES);

		  }
		  ierr = VecAssemblyBegin(partialF_partialalphaMax_FD);
		  ierr = VecAssemblyEnd(partialF_partialalphaMax_FD);

		  // Compare the vectors
		  std::cout<<"/*----Comparing partial derivatives w.r.t. alphaMax ----*/"<<std::endl;
		  std::cout<<"/*--------------------Analytical -----------------------*/"<<std::endl;
		  VecView(dFdalphaMax_partial,PETSC_VIEWER_STDOUT_WORLD);
		  std::cout<<"/*--------------------------FD -------------------------*/"<<std::endl;
		  VecView(partialF_partialalphaMax_FD,PETSC_VIEWER_STDOUT_WORLD);


	  VecDestroy(&partialF_partialalphaMax_FD);

	  return ierr;
}
PetscErrorCode System::RestartSolver()
{
	PetscErrorCode ierr;
	ierr = TSSetInitialTimeStep(ts,0.0,problemData.timestep);
	/*
	 * Set Initial solution
	 */
	this->InitialSolution(networkdm,X);
	ierr = TSSetSolution(ts,X);CHKERRQ(ierr);

	/*
	 * Set initial plastic deformation to zero
	 */
	InitialDesignSolution();

	/*
	 * Clear history
	 */
	timeHistory.clear();
	//timeStepHistory.clear();
	Uhistory.restart();
	alphaMaxhistory.restart();
	U_RungeKuttahistory.restart();

	previous_time = 0.0;

#ifdef DEBUG
	if(problemData.matlabhistory == 1){
		matlab_timestephistory_iterator = matlab_timestephistory.begin();
	}
#endif
	if (solve_in_fd && problemData.matlabhistory != 1){
		timeStepHistory_iterator = timeStepHistory.begin();
	}
	else{
		timeStepHistory.clear();
	}


	return ierr;
}


PetscErrorCode System::FiniteDifference(){
	  PetscErrorCode ierr;


#ifdef PETSC_USE_LOG
	PetscLogStage stage3;
	ierr = PetscLogStageRegister("Finite Diff",&stage3);CHKERRQ(ierr);
	PetscLogStagePush(stage3);
#endif

	  _is_fd_calculated = true;

	  Gradient_FD.resize(designData.N_DesignVariables);

	  // Bead pointer
	  arrayBeads DesignBead;
	  DMNetworkComponentGenericDataType *arr;

	  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

	  PetscReal original_design;
	  PetscReal dt = 1e-8;

	  /*
	   * Set up the solver to use the same time step history than the original
	   * solve
	   */
	  FDSetUpSolver();

	  if (rank == 0){
		  std::cout<<" /*----------------------------------------------------*/ "<<std::endl;
		  std::cout<<" 				Comparing Sensitivities		 				"<<std::endl;
		  std::cout<<" /*----------------------------------------------------*/ "<<std::endl;
	  }
	  // We need to loop over the local vertices
	  PetscInt vStart, vEnd;
	  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);

	  for (unsigned int i = 0; i<designData.N_DesignVariables; i++){
			  PetscInt    offsetd,key;


			  // Change parameter data, first, keep original value
			  original_design = alphaMaxInit[i];

			  alphaMaxInit[i] += dt;

			  //Gradient_FD[i] = -1.0*FunctionValueGlobalOriginal;
			  RestartSolver();
			  /*
			   * Solve
			   */
			  Solve();
			  /*
			   * Calculate Obj Function
			   */
			  ProcessObjFunctionTime();

			  //Gradient_FD[i] = FunctionValueGlobal;
			  Gradient_FD[i] = FunctionValueGlobal;

			  // Change the parameter again, now backwards
			  alphaMaxInit[i] -= 2.0*dt;

			  RestartSolver();
			  /*
			   * Solve
			   */
			  Solve();
			  /*
			   * Calculate Obj Function
			   */
			  ProcessObjFunctionTime();


			  Gradient_FD[i] -= FunctionValueGlobal;

			  // Final value
			  Gradient_FD[i] *= 1.0/(2.0*dt);

			  // Recover density value
			  alphaMaxInit[i] = original_design;

			  if (rank == 0)
				  std::cout<<std::setprecision(10)<<" dFdP_FD["<<i<<"] = "<<Gradient_FD[i]<<std::endl;
	  }

	  RestartSolver();


	  /*
	   * We are finished
	   */
	  solve_in_fd = false;

#ifdef PETSC_USE_LOG
	PetscLogStagePop();
#endif



	  return ierr;




}


PetscErrorCode System::ReadBeadsData(arrayBeads & _arrayBeads)
{


	PetscFunctionBegin;


	PetscErrorCode ierr;
	char                 data_dir[PETSC_MAX_PATH_LEN];
	ierr = PetscOptionsGetString(PETSC_NULL,"-size_data",data_dir,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);

	ierr =PetscOptionsReal("-initial_design", " !=0 for same value for all variables, 0 for reading the initial design from a file","", 0, &designData.initial_design, NULL); CHKERRQ(ierr);



	std::string data_directory(data_dir);

	std::string directory;



	if (data_directory.compare("hex") == 0){
		directory = "Hexagonal";
		std::cout<<"Reading data for large problem"<<std::endl;
	}
	else if (data_directory.compare("chain") == 0){
		directory = "Chain";
		std::cout<<"Reading data for chain problem"<<std::endl;
	}

	/*************************************************************************/
	/*             Get data info from input file                             */
	/*************************************************************************/

	std::string data(directory), elements(directory), nodes(directory), constraints(directory), initialcondition(directory), targetArea(directory);
	data.append("/input_Data.txt");
	elements.append("/input_Elements.txt");
	nodes.append("/input_Nodes.txt");
	constraints.append("/input_Constraints.txt");
	initialcondition.append("/input_InitialCondition.txt");
	targetArea.append("/input_TargetArea.txt");

	std::ifstream finp(data.c_str());
	std::string   keyword;
	std::size_t   found;

	if (finp.is_open()) {

		while(!finp.eof()) {

			std::getline(finp, keyword);
			trim1(keyword);

			found=keyword.find("##");
			if (found!=std::string::npos)  continue;
			else if (keyword == "#totaltime")  {
				finp >> problemData.TotalTime;
				continue;
			}
			else if (keyword == "#timestep")  {
				finp >> problemData.timestep;
				continue;
			}
			else if (keyword == "#finitedifference")  {
				finp >> problemData.finiteDifference;
				continue;
			}
			else if (keyword == "#yieldstress")  {
				finp >> problemData.sigmayielding;
				continue;
			}
			else if (keyword == "#striker")  {
				finp >> problemData.striker;
				continue;
			}
			else if (keyword == "#NumberOfObjectiveFunctions")  {
				finp >> designData.N_ObjFunc;
				continue;
			}
			else if (keyword == "#Objective")  {
				unsigned int mixmax;
				for (unsigned int i=0;i<designData.N_ObjFunc;i++) {
					finp >> mixmax;
					designData.MaxOrMin.push_back(mixmax);
				}
				continue;
			}
			else if (keyword == "#Weight")  {
				finp >> designData.penalty_weight;
				continue;
			}
			else if (keyword == "#Norm")  {
				finp >> designData.SpaceNorm;
				designData.TimeNorm = designData.SpaceNorm;
				continue;
			}
			else if (keyword == "#NumberOfOptimizationConstraints")  {
				finp >> designData.N_ConstrFunc;
				continue;
			}
			else if (keyword == "#Constraints")  {
				PetscReal tolerance;
				unsigned int consbound, constype;
				for (unsigned int i=0;i<designData.N_ConstrFunc;i++) {
					finp >> constype >> consbound >> tolerance;
					designData.ConstraintsBounds.push_back(consbound);
					designData.ConstraintsType.push_back(constype);
					designData.ConstraintsTolerance.push_back(tolerance);
				}
				continue;
			}
			else if (keyword == "#vfmin")  {
				finp >> designData.volfrac_min;
				continue;
			}
			else if (keyword == "#simp")  {
				finp >> designData.SIMP_parameter;
				continue;
			}
			else if (keyword == "#kkttol")  {
				finp >> designData.KKT_tol;
				continue;
			}
			else if (keyword == "#costtol")  {
				finp >> designData.Obj_reltol;
				continue;
			}
			else if (keyword == "#designtol")  {
				finp >> designData.DesVar_reltol;
				continue;
			}
		}
	} else {
		std::cout<<"Cannot open input_Data.txt. Abort"<<std::endl;
		exit(-1);
	}

	/*************************************************************************/
	/*                         Nodal Information                             */
	/*************************************************************************/
	std::ifstream finp_nodes(nodes.c_str());

	if (finp_nodes.is_open()) {
		finp_nodes>>problemData.N_Nodes;
		ierr = PetscMalloc(problemData.N_Nodes*sizeof(class BeadsData),&_arrayBeads);
		for (unsigned int i=0; i<problemData.N_Nodes; i++) {
			finp_nodes>>_arrayBeads[i].node_id;

			finp_nodes>>_arrayBeads[i].coordx;
			finp_nodes>>_arrayBeads[i].coordy;
			finp_nodes>>_arrayBeads[i].coordz;

			// Scaling
			_arrayBeads[i].coordx *= 1e3;
			_arrayBeads[i].coordy *= 1e3;
			_arrayBeads[i].coordz *= 1e3;

			finp_nodes>>_arrayBeads[i].type>>_arrayBeads[i].E>>_arrayBeads[i].rho>>_arrayBeads[i].nu>>_arrayBeads[i].R;

			//Scaling
			_arrayBeads[i].E   = _arrayBeads[i].E*1e-9;
			_arrayBeads[i].rho = _arrayBeads[i].rho*1e-3;
			_arrayBeads[i].R   = _arrayBeads[i].R*1e3;

			// Initialize density value
			_arrayBeads[i].density_design = 1.0;

			_arrayBeads[i].constrained_dof_y = false;
			_arrayBeads[i].constrained_dof_x = false;


		}
		finp_nodes.close();
	} else {
		printf("Cannot open input_Nodes.txt \n");
		exit(-1);
	}




	/*************************************************************************/
	/*                         Element Information                           */
	/*************************************************************************/
	std::ifstream finp_elem(elements.c_str());
	if (finp_elem.is_open()) {
		finp_elem>>problemData.N_Elements;
		ierr = PetscMalloc(problemData.N_Elements*sizeof(struct _p_EDGEDATA),&edgeData);
		problemData.connectivity = new int[2*problemData.N_Elements];
		unsigned int element_id;
		PetscReal alpha;
		for (unsigned int i=0; i<problemData.N_Elements; i++) {
			finp_elem>>element_id;
			finp_elem>>problemData.connectivity[2*i];
			problemData.connectivity[2*i] -= 1;
			finp_elem>>problemData.connectivity[2*i+1];
			problemData.connectivity[2*i+1] -= 1;
			finp_elem>>alpha;

			edgeData[i].VariableGroup = alpha - 1;

			if (design_elements.size() < alpha){
				std::vector<PetscInt> dummy;
				design_elements.push_back(dummy);
				design_elements[alpha - 1].push_back(element_id);
			}
			else{
				design_elements[alpha - 1].push_back(element_id);
			}
		}
		finp_elem.close();
	} else {
		printf("Cannot open input_Elements.txt \n");
		exit(-1);
	}


	designData.N_DesignVariables = design_elements.size();

	/*
	 * Set the number of variables for the arrays
	 */
	Mdot_product_result = new PetscScalar[designData.N_DesignVariables];

	for (unsigned int i = 0; i<designData.N_DesignVariables; i++)
		Mdot_product_result[i] = 0.0;

	/*************************************************************************/
	/*                         Constraint Information                        */
	/* ConstraintsId = [nodeid dof]                                           */
	/* ConstraintsValue = up  (u prescribed)                                  */
	/*************************************************************************/
	std::ifstream finp_bc(constraints.c_str());
	if (finp_bc.is_open()) {
		finp_bc>>problemData.N_Constraints;
		problemData.constraints_list.resize(problemData.N_Constraints);
		PetscReal constraint_value;
		unsigned int bead_number, dof;
		for (unsigned int i=0; i<problemData.N_Constraints; i++) {
			finp_bc>>bead_number;
			finp_bc>>dof;
			finp_bc>>constraint_value;
			problemData.constraints_list[i] = std::make_pair(bead_number,dof);
			if (dof == 1){
				_arrayBeads[bead_number-1].constraint_value_x = constraint_value;
				_arrayBeads[bead_number-1].constrained_dof_x = true;
			}
			else if (dof == 2){
				_arrayBeads[bead_number-1].constraint_value_y = constraint_value;
				_arrayBeads[bead_number-1].constrained_dof_y = true;
			}
		}
		finp_bc.close();
	} else {
		printf("Cannot open input_Constraints.txt \n");
		exit(-1);
	}


	/*************************************************************************/
	/*                  Initial Condition Information                        */
	/* InitialConditionId = [nodeid dof]                                     */
	/* InitialConditionValue = Initial Velocity                              */
	/*************************************************************************/
	std::ifstream finp_ic(initialcondition.c_str());
	if (finp_ic.is_open()) {
		finp_ic>>problemData.N_InitialConditions;
		if (problemData.transient){
			if (problemData.striker == 0) {
					finp_ic>>problemData.striker_ID;
					finp_ic>>problemData.striker_DOF;
					finp_ic>>problemData.InitialVelocity;
					problemData.InitialVelocity = problemData.InitialVelocity*1e-3;
			}
			else {
				finp_ic>>problemData.striker_ID;
				finp_ic>>problemData.striker_DOF;
			}
		}
		finp_ic.close();
	} else {
		printf("Cannot open input_InitialCondition.txt \n");
		exit(-1);
	}


	/*************************************************************************/
	/*                  Target Area Information				                 */
	/*************************************************************************/

	// Resize array that will hold the target areas
	designData.TargetArea_list.resize(designData.N_ObjFunc);
	std::ifstream finp_ta(targetArea.c_str());
	if (finp_ta.is_open()) {
		for (unsigned int i = 0;i<designData.N_ObjFunc; i++){
			int N_target_elem;
			finp_ta>>N_target_elem;
			designData.TargetArea_list[i].resize(N_target_elem);
			for (int j = 0;j<N_target_elem; j++)
				finp_ta>>designData.TargetArea_list[i][j];
		}
	} else {
		printf("Cannot open input_TargetArea.txt \n");
		exit(-1);
	}


	/*
	 * Fill the array with ones.
	 */
	CapitalOmega.resize(designData.N_ObjFunc);
	FunctionObjetivoPartial.resize(designData.N_ObjFunc);
	std::fill(CapitalOmega.begin(),CapitalOmega.end(),1.0);
	std::fill(FunctionObjetivoPartial.begin(),FunctionObjetivoPartial.end(),0.0);


	/*
	 * What kind of problem is it? Large or small deformations?
	 */
	PetscInt large = 0;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-large_deformations",&large,NULL);CHKERRQ(ierr);

	if (large == 1)
		problemData.large_deformations = true;
	else
		problemData.large_deformations = false;

	/*
	 * Set the constitutive law accordingly
	 */

	elementFunction.set_large_deformations(problemData.large_deformations);

	problemData.problem_info();
	designData.design_info();


	ierr =PetscOptionsInt("-printdisplacements", "1 to print displacements","", 0, &problemData.printdisplacements, NULL);






	PetscFunctionReturn(0);

	return ierr;
}

PetscErrorCode System::InitialSolution(DM networkdm,Vec X)
{
  PetscErrorCode ierr;
  Vec            localX;
  PetscScalar    *xarr;
  DMNetworkComponentGenericDataType *arr;

  PetscInt      v;
  PetscInt      vStart,vEnd, eStart, eEnd;
  PetscInt      offsetfrom;

  arrayBeads beadsdata;
  PetscFunctionBegin;

  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = VecSet(X,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  for (v=vStart; v < vEnd; v++) {
	PetscInt    offsetd,key;
	ierr = DMNetworkGetVariableOffset(networkdm,v,&offsetfrom);CHKERRQ(ierr);
	ierr = DMNetworkGetComponentTypeOffset(networkdm,v,0,&key,&offsetd);CHKERRQ(ierr);
	// Cast component data
	beadsdata = (arrayBeads)(arr + offsetd);
	// Velocity components equal to zero everywhere
	xarr[offsetfrom] = 0.0;
	xarr[offsetfrom + 1] = 0.0;

	if (problemData.loading == 0 && beadsdata->type == 3)
		xarr[offsetfrom + 2] = 0.5;

	// Displacement components
	if (beadsdata->constrained_dof_x)
		xarr[offsetfrom + 2] = 0.0;
	if (beadsdata->constrained_dof_y)
		xarr[offsetfrom + 3] = 0.0;

  }

  /*
   * Initial value of the objective function
   */
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  if (eStart == 0){
	  ierr = DMNetworkGetVariableOffset(networkdm,eStart,&offsetfrom);CHKERRQ(ierr);
	  xarr[offsetfrom] = 0.0;
  }

  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(networkdm,localX,INSERT_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localX,INSERT_VALUES,X);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);

  return ierr;
}


PetscErrorCode System::InitialDesignSolution(){

	PetscErrorCode ierr;


	PetscScalar * alphaMaxarr;


	VecGetArray(alphaMax, &alphaMaxarr);

	/* Edge Data necessary to know which variable group it belongs to */
	EDGEDATA elementData;
	PetscInt key, offsetd_eStart;
	DMNetworkComponentGenericDataType * arr;
	ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);



	for (unsigned int i = 0; i < problemData.N_Elements; i++){

		// Set initial design value, first, we need to know which variable group it belongs to
		DMNetworkGetComponentTypeOffset(networkdm,i,0,&key,&offsetd_eStart);
		elementData = (EDGEDATA)(arr + offsetd_eStart);

		alphaMaxarr[i] = alphaMaxInit[elementData->VariableGroup];

	}

	VecRestoreArray(alphaMax, &alphaMaxarr);


}



