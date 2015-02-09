/*
 * vechistoryRK.h
 *
 *  Created on: Jan 12, 2015
 *      Author: miguel
 */

#ifndef VECHISTORYRK_H_
#define VECHISTORYRK_H_


#include <petscvec.h>
#include <cassert>
#include <iterator>
class VecHistoryRungeKutta{

	public:
		VecHistoryRungeKutta(){
			sizehistory = 0;
			current_index = 0;
			Runge_Kutta_current_index= 0;
			Stages = 4;

		}


		void SetStages(unsigned int & _Stages){
			Stages = _Stages;
		}

		void destroyhistory(){

			std::vector<std::vector<Vec> >::iterator begin = history.begin();
			std::vector<std::vector<Vec> >::iterator end = history.end();
			for (;begin!=end;begin++){
				std::vector<Vec>::iterator beginRK = (*begin).begin();
				std::vector<Vec>::iterator endRK = (*begin).end();
				for (;beginRK!=endRK;beginRK++)
					VecDestroy(&(*beginRK));
			}

			sizehistory = 0;
			current_index = 0;
			Runge_Kutta_current_index= 0;
		}

		/*
		 * We need to destroy the Vecs before calling the destructor because
		 * this is called after we have called PetscFinalize()
		 */

		unsigned int sizehistory;

		unsigned int current_index;

		unsigned int Runge_Kutta_current_index;

		unsigned int Stages;

		std::vector<std::vector<Vec> > history;



		/*
		 *
		 * Member functions
		 *
		 */

		std::vector<Vec> & operator[](unsigned int index){
			assert(index < sizehistory  && "index exceeds the size of the history");
			return history[index];
		}

		/*
		 *  ----------- ITERATORS ------------------
		 */
		std::vector<std::vector<Vec> >::iterator begin(){
			return history.begin();
		}

		/*
		 * end() and rbegin() cannot return the same thing by definition of the
		 * iterator. end() should give a pointer to one element beyond the current end
		 * rbegin() should point to the last element
		 */
		std::vector<std::vector<Vec> >::iterator end(){
			std::vector<std::vector<Vec> >::iterator  current_end = history.begin();
			std::advance(current_end,sizehistory);
			return current_end;
		}

		std::vector<std::vector<Vec> >::reverse_iterator rbegin(){
			std::vector<std::vector<Vec> >::reverse_iterator  current_end = history.rbegin();
			std::advance(current_end,history.size() - sizehistory);
			return current_end;
		}

		std::vector<std::vector<Vec> >::reverse_iterator rend(){
			return history.rend();
		}

		unsigned int & size(){

			return sizehistory;
		}

		void new_timestep(Vec vector){
			if (current_index >= sizehistory){
				std::vector<Vec> temp;
				for (unsigned int i=0 ; i<Stages; i++){
					Vec NewVec;
					VecDuplicate(vector,&NewVec);
					temp.push_back(NewVec);
				}
				history.push_back(temp);
				sizehistory++;
			}
			Runge_Kutta_current_index = 0;
			current_index++;

		}

		void insert(Vec vector){
			VecCopy(vector,(history[current_index-1][Runge_Kutta_current_index%Stages]));
			Runge_Kutta_current_index++;
		}
		/*
		 * We keep the vectors so we can recycle them
		 */
		void restart(){
			current_index = 0;
			Runge_Kutta_current_index = 0;
		}
};





#endif /* VECHISTORYRK_H_ */
