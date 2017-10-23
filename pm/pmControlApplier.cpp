/*************************************************************************
*  Copyright (C) 2017 by Jakub Lengiewicz <jleng@ippt.pan.pl>            *
*  Copyright (C) 2017 by Pawel Chodkiewicz <pawel@chodkiewicz.com.pl>    * 
*  Copyright (C) 2017 by Pawel Holobut <pholob@ippt.pan.pl>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "pmControlApplier.hpp"
#include<core/Scene.hpp>
#include<pkg/pm/pmState.hpp>

void PmControlApplier::action(){
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){             
		PmState* st = static_cast<PmState*>(b->state.get());		
		if(st->setMinIndex) st->minIndex = 0; 		          
	}
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(I->isReal()){
			const shared_ptr<Body>& b1 = Body::byId(I->getId1(),scene);
			const shared_ptr<Body>& b2 = Body::byId(I->getId2(),scene);
			PmState* st1 = static_cast<PmState*>(b1->state.get());
			PmState* st2 = static_cast<PmState*>(b2->state.get());
			if(st1->setMinIndex){
				if(st1->minIndex == 0) 
					st1->minIndex = st2->index;
				else{
					int l = st2->index - st1->minIndex;
					if(l > 0){
						if(l <= indexRange || l >= maxLoopLength) st1->minIndex = st2->index;
					}
					else{
						if(l < -indexRange && l > -maxLoopLength) st1->minIndex = st2->index;
					}
				}
			}
			if(st2->setMinIndex){
				if(st2->minIndex == 0)
					st2->minIndex = st1->index;
				else{
					int l = st1->index - st2->minIndex;
					if(l > 0){
						if(l <= indexRange || l >= maxLoopLength) st2->minIndex = st1->index;
					}
					else{
						if(l < -indexRange && l > -maxLoopLength) st2->minIndex = st1->index;				
					}
				}
			}
		}		
	}
}
YADE_PLUGIN((PmControlApplier));
CREATE_LOGGER(PmControlApplier);

