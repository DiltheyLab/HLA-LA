/*
 * linearALTs.h
 *
 *  Created on: 14.01.2016
 *      Author: AlexanderDilthey
 */

#ifndef LINEARALTS_LINEARALTS_H_
#define LINEARALTS_LINEARALTS_H_

#include <string>
#include <map>
#include <utility>
#include "../simulator/simulator.h"
#include "api/BamAlignment.h"
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include "../pathFinder.h"

namespace linearALTs {

class linearALTs {
public:
	linearALTs(std::string directory, const pathFinder& pF);
	virtual ~linearALTs();
	std::pair<std::string, std::string> simulateEqualLengthDiploidHaplotypes(std::string outputDirectory, const simulator::simulator& S, std::map<std::string, std::set<std::string>>& read2gene_R1, std::map<std::string, std::set<std::string>>& read2gene_R2, double haploidCoverage = 1, bool withError = true);
	void simulateGeneSequences(std::string outputDirectory, const simulator::simulator& S, bool simulateFromGraph, double haploidCoverage = 1, bool withError = true);

	std::pair<std::string, std::string> haplotypeLikelihoods(std::string outputDirectory, double insertSize_mean, double insertSize_sd, int fallBackReadLength, std::pair<std::string, std::string> printPair);
	void reads2Genes(std::string outputDirectory, double insertSize_mean, double insertSize_sd, int fallBackReadLength, std::map<std::string, std::set<std::string>>& true_read2gene_R1, std::map<std::string, std::set<std::string>>& true_read2gene_R2, std::map<std::string, std::map<std::string, int>>& ret_reads_trueGene_2_inferredGene, std::map<std::string, std::map<std::string, int>>& ret_reads_inferredGene_2_trueGene, const std::set<std::string>& restrictToHaplotypesForGenes, double T);

	bool alignedReadPair_strandsValid(const BamTools::BamAlignment& alignment_read1, const BamTools::BamAlignment& alignment_read2);
	int alignedReadPair_distance(const BamTools::BamAlignment& alignment_read1, const BamTools::BamAlignment& alignment_read2);

	void extractReads_equalLengthHaplotypes(std::string outputDirectory, std::string BAM);
	void extractReads_geneGraph(std::string outputDirectory, std::string BAM);
	void extractReads_extendedReferenceGenome(std::string outputDirectory, std::string BAM);

	std::string getGeneGraphPRGDirectory();

	static void extractReadsFromBAM(std::string outputDirectory, std::string BAM, const std::map<std::string, std::pair<int, int>>& regions, bool removePairNumbers = false);

protected:

	const pathFinder& pF;

	std::map<std::string, std::pair<int, int>> extendedReferenceGenome_coveredRegions;
	std::map<std::string, std::string> haplotypeID_2_extendedReferenceGenomeID;

	std::string fastaFile_equalLengthHaplotypes;
	std::map<std::string, std::string> equalLengthHaplotypes;
	std::vector<std::string> equalLengthHaplotypeIDs;
	std::map<std::string, std::pair<int, int>> equalLengthHaplotypes_inExtendedRef_0based;

	std::string geneSequencesFile;
	std::string simulateGeneSequences_alignedSequences;
	std::string geneGraphDirectory;
	std::string geneGraphFile;
	std::string fastaFile_explicitGenes;
	std::map<std::string, std::string> explicitGenes_haplotypeID_2_fastaID;
	std::set<std::string> explicitGenes_fastaIDs_genes;

	std::map<std::string, std::map<int, std::set<std::string>>> equalLengthHaplotypes_annotations_perPosition;

	static double scoreAlignment(const BamTools::BamAlignment& alignment);

	static std::string CIGARasString(const BamTools::BamAlignment& al);

	std::vector<double> processCollectedAlignments(const std::vector<BamTools::BamAlignment>& Alignments_R1, const std::vector<BamTools::BamAlignment>& Alignments_R2, std::pair<std::map<std::string, double>, std::map<std::string, double>>& hitGenes, size_t& processed_r1, size_t& processed_r2, const std::vector<std::string>& references_IDs, const std::set<std::string>& restrictToHaplotypes, const boost::math::normal& rnd_InsertSize, int fallBackReadLength, bool verbose);

};

} /* namespace linearALTs */


#endif /* LINEARALTS_LINEARALTS_H_ */

