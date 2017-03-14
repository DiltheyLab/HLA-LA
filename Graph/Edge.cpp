/*
 * Edge.cpp
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#include "Edge.h"
#include <assert.h>
#include <string>

Edge::Edge() {
	count = 0.0;
	label = "";
	pgf_protect = false;
	seedChain_isGraphGap = false;
}

std::string Edge::getEmission()
{
	assert(emission != 0);
	std::string forReturn;
	forReturn.push_back(emission);
	return forReturn;
}
