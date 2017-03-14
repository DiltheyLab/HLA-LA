/*
 * Node.h
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef NODE_H_
#define NODE_H_


#include <vector>
#include <string>
#include <set>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>


#include "Graph.h"
#include "Edge.h"

namespace boost {
namespace serialization {
class access;
}
}

class nodeHoldsHaplo {
private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
				ar & haploID;
				ar & P;
				ar & indexInHaplotypeVector;
	}	
public:
	int haploID;
	double P;
	int indexInHaplotypeVector;
};



struct HaploAttachedToNode {
	Node* n;
	int indexInNodeHaploInfo;
	double P;
};

using namespace std;




class Node {
private:
    friend class boost::serialization::access;

    template<typename Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
 				// ar & g;
				// ar & Incoming_Edges;
				// ar & Outgoing_Edges;
				ar & haplotypes;
				ar & level;
				ar & terminal;
	}
public:
	Node();

	Graph* g;

	set<Edge*> Incoming_Edges;
	set<Edge*> Outgoing_Edges;

	vector<nodeHoldsHaplo> haplotypes;

	double Sum_Incoming();
	double Sum_Outgoing();

	unsigned int level;
	bool terminal;
};

#endif /* NODE_H_ */

