/*
 * simpleGraphSimulator.h
 *
 *  Created on: 22.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef GRAPH_GRAPHSIMULATOR_SIMPLEGRAPHSIMULATOR_H_
#define GRAPH_GRAPHSIMULATOR_SIMPLEGRAPHSIMULATOR_H_

#include <map>
#include <vector>
#include <string>
#include <utility>

#include "../Graph.h"
#include "../../mapper/reads/PRGContigBAMAlignment.h"
#include "../../Graph/HaplotypePanel.h"
#include "../../pathFinder.h"

class simpleGraphSimulator {
protected:
	const pathFinder& pF;
	Graph* g;
	std::map<std::string, std::vector<std::vector<Edge*>>> contigs_by_category;
	void init();
	bool noGaps_switchContigs;
	unsigned int graphLength;
	unsigned int forPRG_mutations;
	unsigned int forPRG_largeGaps;
	unsigned int additionalForReads_mutations;

	double mutationDensity;
	double gapStartFrequency;
	double gapExpectedLength;

	std::set<Edge*> generated_edges;
	std::vector<Edge*> string_2_edgePath(std::string in);

	HaplotypePanel hP;

	std::map<std::string, std::pair<std::string, unsigned int>> PRGID_2_internalID;
	std::map<std::string, unsigned int> PRGID_2_INTid;

public:
	simpleGraphSimulator(bool noGaps_switchContigs_, const pathFinder& pF);
	virtual ~simpleGraphSimulator();
	std::map<std::string, std::vector<mapper::reads::PRGContigBAMAlignment>> simulateBAMAlignments(double contigCoverage, std::string qualityMatrixFile_, int read_length_, double insertSize_mean_, double insertSize_sd_, bool error);

	Graph* getGraph();

	void storeLikeRealPRG(std::string directory);

};


#endif /* GRAPH_GRAPHSIMULATOR_SIMPLEGRAPHSIMULATOR_H_ */
