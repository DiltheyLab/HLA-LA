/*
 * simulator.h
 *
 *  Created on: 16.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef SIMULATOR_SIMULATOR_H_
#define SIMULATOR_SIMULATOR_H_

#include <string>
#include <vector>
#include "../Graph/Graph.h"
#include "readSimulator.h"

namespace simulator {

class simulator {
public:
	simulator(std::string qualityMatrixFile_, int read_length_ = 101, double insertSize_mean_ = 100, double insertSize_sd_ = 10);
	void simulateFromGraph(Graph* g, int simulatedGraphGenomes, std::string outputDirectory, double haploidCoverage, bool withError) const;
	void simulateNormalGenome(std::string graphDir, std::string referenceGenomePath, std::string outputDirectory, double haploidCoverage, bool withError, unsigned int threads);
	void simulateFromDiploidHaplotypes(const std::string& h1, const std::string& h2, const std::string& outputDirectory, double haploidCoverage, bool withError) const;
	void simulateFromHaploidAlignedSequence(const std::string& S, const std::string& outputDirectory, double haploidCoverage, bool withError) const;

	virtual ~simulator();

protected:
	int read_length;
	double insertSize_mean;
	double insertSize_sd;
	unsigned int threads;
	std::string qualityMatrixFile;
	readSimulator* rS;
};

} /* namespace simulator */

#endif /* SIMULATOR_SIMULATOR_H_ */
