/*
 * Node.cpp
 *
 *  Created on: 25.05.2011
 *      Author: Alexander Dilthey
 */

#include "Node.h"

Node::Node() {
	terminal = false;
	level = 0;
}

double Node::Sum_Incoming()
{
	double sum = 0.0;
	for(set<Edge*>::iterator E = Incoming_Edges.begin(); E != Incoming_Edges.end(); E++)
	{
		sum += (*E)->count;
	}
	return sum;
}

double Node::Sum_Outgoing()
{
	double sum = 0.0;
	for(set<Edge*>::iterator E = Outgoing_Edges.begin(); E != Outgoing_Edges.end(); E++)
	{
		sum += (*E)->count;
	}
	return sum;
}
