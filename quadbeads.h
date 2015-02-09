/*
 * quadbeads.h
 *
 *  Created on: Jan 2, 2015
 *      Author: miguel
 */

#ifndef QUADBEADS_H_
#define QUADBEADS_H_

#include "IpTNLP.hpp"


#include "IpIteratesVector.hpp"
#include "IpIpoptData.hpp"
#include "IpJournalist.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpRestoIpoptNLP.hpp"

#include "assert.h"
#include <vector>

#include "system.h"
#include <petscdmnetwork.h>

using namespace Ipopt;

class QuadBeads : public TNLP
{
	public:

		QuadBeads();


		/** default destructor */
		virtual ~QuadBeads();

		/**@name Overloaded from TNLP */
		//@{
		/** Method to return some info about the nlp */
		virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
						  Index& nnz_h_lag, IndexStyleEnum& index_style);

		/** Method to return the bounds for my problem */
		virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
							 Index m, Number* g_l, Number* g_u);

		/** Method to return the starting point for the algorithm */
		virtual bool get_starting_point(Index n, bool init_x, Number* x,
								bool init_z, Number* z_L, Number* z_U,
								Index m, bool init_lambda,
								Number* lambda);

		/** Method to return the objective value */
		virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

		/** Method to return the gradient of the objective */
		virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

		/** Method to return the constraint residuals */
		virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

		/** Method to return:
		*   1) The structure of the jacobian (if "values" is NULL)
		*   2) The values of the jacobian (if "values" is not NULL)
		*/
		virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
						Index m, Index nele_jac, Index* iRow, Index *jCol,
						Number* values);

		//@}

		/** @name Solution Methods */
		//@{
		/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
		virtual void finalize_solution(SolverReturn status,
							   Index n, const Number* x, const Number* z_L, const Number* z_U,
							   Index m, const Number* g, const Number* lambda,
							   Number obj_value,
							   const IpoptData* ip_data,
							   IpoptCalculatedQuantities* ip_cq);

		virtual bool intermediate_callback(AlgorithmMode mode,
										  Index iter, Number obj_value,
										  Number inf_pr, Number inf_du,
										  Number mu, Number d_norm,
										  Number regularization_size,
										  Number alpha_du, Number alpha_pr,
										  Index ls_trials,
										  const IpoptData* ip_data,
										  IpoptCalculatedQuantities* ip_cq)
		 {
			OrigIpoptNLP* orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
			if (orignlp == NULL){
				RestoIpoptNLP* restonlp = dynamic_cast<RestoIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
				orignlp = dynamic_cast<OrigIpoptNLP*>(&restonlp->OrigIpNLP());
			}
			assert(ip_cq);
			TNLPAdapter* tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
			assert(tnlp_adapter);

			tnlp_adapter->ResortX(*(ip_data->curr()->x()), my_x);

			for (unsigned int k=0;k<Beads.designData.N_DesignVariables;k++)
				fout_volfrac<<my_x[k]<<"\t";

			fout_volfrac<<std::endl;
			iterations += 1;

			//delete orignlp;
			//delete tnlp_adapter;

			return true;
		 }


	private:

		QuadBeads(const QuadBeads&);
		QuadBeads& operator=(const QuadBeads&);

		PetscErrorCode transfer_densities(const Number *x);

		System Beads;
		double			* my_x;
		std::ofstream fout_volfrac;
		int iterations;


		VecScatter scatterctx;
		Vec seq_vec;

		// Jacobian information for the constraint gradient
		int     *jacIndexVars, *jacIndexCons;

};


#endif /* QUADBEADS_H_ */
