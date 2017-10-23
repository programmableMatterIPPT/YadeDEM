/*************************************************************************
*  Copyright (C) 2017 by Jakub Lengiewicz <jleng@ippt.pan.pl>            *
*  Copyright (C) 2017 by Pawel Chodkiewicz <pawel@chodkiewicz.com.pl>    * 
*  Copyright (C) 2017 by Pawel Holobut <pholob@ippt.pan.pl>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once
#include <core/GlobalEngine.hpp>
#include <stdexcept>


class PmControlApplier: public GlobalEngine{
	virtual void action();
	YADE_CLASS_BASE_DOC_ATTRS(PmControlApplier,GlobalEngine,"Class for controlling which neighboring bodies adhere to each other - based on their materials and indexes.",
		((int,maxLoopLength,100000,,"Convention: index1 is considered lower than index2 iff index2-indexRange<=index1<index2 or index1<=index2-maxLoopLength+1."))
		((int,indexRange,1,,"Convention: index1 is considered lower than index2 iff index2-indexRange<=index1<index2 or index1<=index2-maxLoopLength+1."))
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PmControlApplier);

