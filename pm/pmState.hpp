/*************************************************************************
*  Copyright (C) 2017 by Jakub Lengiewicz <jleng@ippt.pan.pl>            *
*  Copyright (C) 2017 by Pawel Chodkiewicz <pawel@chodkiewicz.com.pl>    * 
*  Copyright (C) 2017 by Pawel Holobut <pholob@ippt.pan.pl>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#include<lib/serialization/Serializable.hpp>
#include<lib/multimethods/Indexable.hpp>
#include<core/Dispatcher.hpp>
#include<core/State.hpp>


class PmState: public State{
	public:
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(PmState,State,"PmState of a body (spatial configuration, internal variables).",
		((bool,setMinIndex,false,,"Should minIndex of a given body be set to the maximum value of indexes of its neighboring bodies at every step?"))
		((int,minIndex,0,,"If the value of index of a given neighboring body != minIndex, then the bodies do not cohere"))
		((int,index,0,,"Index of this body")) 
		((int,indexn,0,,"Second index of this body")),
 		createIndex();
        );
	REGISTER_CLASS_INDEX(PmState,State);
};
REGISTER_SERIALIZABLE(PmState);

