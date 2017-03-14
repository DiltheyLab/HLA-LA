/*
 * Edge.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <vector>
#include <string>
#include <map>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Graph.h"
#include "Node.h"

using namespace std;

namespace boost {
namespace serialization {
class access;
}
}

class Node;

class Edge {
private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
 				ar & From;
				ar & To; 
				ar & count;
				ar & emission;
				ar & locus_id;
				ar & label;
				ar & pgf_protect;
				ar & seedChain_isGraphGap;
	}

public:
	Edge();

	Node* From;
	Node* To;

	double count;
	unsigned char emission;

	string locus_id;
	string label;

	bool pgf_protect;

	std::string getEmission();

	bool seedChain_isGraphGap;
};

#endif /* EDGE_H_ */
