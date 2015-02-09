/*
 * vechistory.h
 *
 *  Created on: Dec 10, 2014
 *      Author: miguel
 */

#ifndef VECHISTORY_H_
#define VECHISTORY_H_

#include <petscvec.h>
#include <cassert>
#include <iterator>
#include <vector>
class VecHistory{

	public:
		VecHistory();

		VecHistory(VecHistory & History);

		~VecHistory();

		VecHistory & operator=(VecHistory & History);

		void destroyhistory();

		/*
		 * We need to destroy the Vecs before calling the destructor because
		 * this is called after we have called PetscFinalize()
		 */

		unsigned int sizehistory;

		unsigned int current_index;

		std::vector<Vec> history;

		Vec & operator[](unsigned int index);

		/*
		 *  ----------- ITERATORS ------------------
		 */
		std::vector<Vec>::iterator begin(){
			return history.begin();
		}


		/*
		 * end() and rbegin() cannot return the same thing by definition of the
		 * iterator. end() should give a pointer to one element beyond the current end
		 * rbegin() should point to the last element
		 */
		std::vector<Vec>::iterator end(){
			std::vector<Vec>::iterator  current_end = history.begin();
			std::advance(current_end,sizehistory);
			return current_end;
		}

		std::vector<Vec>::reverse_iterator rbegin(){
			std::vector<Vec>::reverse_iterator  current_end = history.rbegin();
			std::advance(current_end,history.size() - sizehistory);
			return current_end;
		}

		std::vector<Vec>::reverse_iterator rend(){
			return history.rend();
		}

		unsigned int & size(){ return sizehistory;}

		void insert(Vec vector);

		/*
		 * We keep the vectors so we can recycle them
		 */
		void restart(){ current_index = 0;}

	private:
		void copy(VecHistory & History);
};



#endif /* VECHISTORY_H_ */
