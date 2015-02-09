/*
 * vechistory.C
 *
 *  Created on: Feb 4, 2015
 *      Author: miguel
 */

#include "vechistory.h"


VecHistory::VecHistory(){
	sizehistory = 0;
	current_index = 0;
}

VecHistory::VecHistory(VecHistory & History){

	copy(History);

}

VecHistory::~VecHistory(){

	destroyhistory();

}

VecHistory & VecHistory::operator=(VecHistory & History){

	if (this != &History){

		destroyhistory();
		copy(History);

	}

	return *this;

}

void VecHistory::copy(VecHistory & History){
	sizehistory = History.sizehistory;
	current_index = History.current_index;


	std::vector<Vec>::iterator begin_copy = History.begin();
	std::vector<Vec>::iterator end_copy = History.end();
	for (;begin_copy!=end_copy;begin_copy++){
		Vec NewVec;
		VecDuplicate(*begin_copy,&NewVec);
		VecCopy(*begin_copy,NewVec);
		history.push_back(NewVec);
	}



}



void VecHistory::destroyhistory(){

	std::vector<Vec>::iterator begin = history.begin();
	std::vector<Vec>::iterator end = history.end();
	for (;begin!=end;begin++){
		VecDestroy(&(*begin));
	}

	sizehistory = 0;
	current_index = 0;
}

Vec & VecHistory::operator[](unsigned int index){
	assert(index < sizehistory  && "index exceeds the size of the history");
	return history[index];
}

void VecHistory::insert(Vec vector){
	if (current_index < sizehistory){
		VecCopy(vector,(history[current_index]));
		current_index++;
	}
	else{
		Vec NewVec;
		VecDuplicate(vector,&NewVec);
		VecCopy(vector,NewVec);
		history.push_back(NewVec);
		sizehistory++;
		current_index++;
	}
}



