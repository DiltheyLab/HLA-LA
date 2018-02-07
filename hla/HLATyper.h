/*
 * HLATyper.h
 *
 *  Created on: 28.10.2015
 *      Author: AlexanderDilthey
 */

#ifndef HLA_HLATYPER_H_
#define HLA_HLATYPER_H_

#include "../Graph/Graph.h"
#include "../mapper/reads/oneReadPair.h"
#include "../mapper/reads/verboseSeedChain.h"

#include "../simulator/readSimulator.h"

#include "oneExonPosition.h"
#include "../intervalTree/IntervalTree.h"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>

namespace hla {

class HLATyper {
protected:
	Graph* g;
	std::string graphDir;
	simulator::readSimulator* rS;

	IntervalTree<std::string>* interestingLevels;

	std::vector<std::string> loci_for_HLAtyping;
	std::map<std::string, std::vector<std::string> > loci_2_exons;

	std::vector<std::string> files_in_graphDir;
	std::vector<std::string> graphLoci;
	std::map<std::string, unsigned int> graphLocus_2_levels;


	int simulations_read_length;
	int simulations_haploidCoverage;
	std::string simulations_qualityMatrixFile;

	int threshold_reportColumn_forPresenceOfUnaccountedAlleles_minCoverage;
	double threshold_reportColumn_forPresenceOfUnaccountedAlleles_minAlleleFraction;
	bool highCoverage_filter_alleles;
	int highCoverage_minCoverage;
	double highCoverage_minAlleleFreq;

	bool filterFirst20;
	int filterFirst20N;
	double filterFirst20MinProp;
	int filterFirst20MinProp_limitKickOutPerRead = 2;

	std::set<std::string> graphGenes;
	std::map<std::string, std::pair<std::vector<std::string>, std::vector<std::string>>> segments_per_graphGene;
	std::map<std::string, std::vector<std::string>> segments_per_graphGene_all;
	std::map<std::string, std::pair<int, int>> graphgene_levelBoundaries;

	std::map<std::string, std::map<std::string, std::map<std::string, std::vector<std::string>>>> sequences_per_graphGeneSegment;
	std::map<std::string, std::map<std::string, std::vector<std::string>>> levelNames_per_graphGeneSegment;


	std::set<std::string> get_complete_genomic_types_per_gene(std::string gene);
	std::set<std::string> get_complete_exonic_types_per_gene(std::string gene);

	std::string find_file_for_exon(std::string locus, std::string exon);

	void oneReadAlignment_2_exonPositions_paired(const mapper::reads::verboseSeedChain& alignment, const mapper::reads::oneRead& read, std::vector<oneExonPosition>& ret_exonPositions, const mapper::reads::verboseSeedChain& paired_alignment, const mapper::reads::oneRead& paired_read, int read_1_or_2, int combined_exon_sequences_graphLevels_min, int combined_exon_sequences_graphLevels_max, const std::map<int, unsigned int>& graphLevel_2_exonPosition);

	void oneReadAlignment_2_exonPositions_unpaired(const mapper::reads::verboseSeedChain& alignment, const mapper::reads::oneRead& read, std::vector<oneExonPosition>& ret_exonPositions, int combined_exon_sequences_graphLevels_min, int combined_exon_sequences_graphLevels_max, const std::map<int, unsigned int>& graphLevel_2_exonPosition);

	static double alignmentFractionOK (const mapper::reads::verboseSeedChain& r);
	static std::pair<double, double> meanMedian (std::vector<double> L);
	static std::string printPerc (double v1, double v2);

	static double alignmentWeightedOKFraction(const mapper::reads::oneRead& underlyingRead, const mapper::reads::verboseSeedChain& alignment);

	static std::vector<oneExonPosition> removeDoublePositionsFromRead(const std::vector<oneExonPosition>& positions);

	static std::string removeLocusAndStar(std::string HLAType);
	
	static bool alleles_compatible(std::string inferredAllele, std::string trueAllele);

	bool can_translateToG_locus(std::string locus);
	std::string translate_allele_list_to_G_allele(const std::vector<std::string>& alleles, bool& ret_perfectly);
	void read_G_alleles();

	std::map<std::string, std::string> alleles_to_G;
	std::set<std::string> G_loci;

	static int kMerSeq_is_kMerKey(const std::string& kMer);
	static std::string kMer_canonical_representation(const std::string& kMer);
	static std::string kMer_canonical_representation(const std::string& kMer, bool& was_inverted);
	static double simpleChiSq(std::vector<double> observed, std::vector<double> expected);

public:
	HLATyper(Graph *g, std::string graphDir, std::string simulations_qualityMatrixFile_);
	virtual ~HLATyper();

	void HLATypeInference(const std::vector<mapper::reads::oneReadPair>& rawPairedReads, const std::vector<mapper::reads::verboseSeedChainPair>& alignedPairedReads, const std::vector<mapper::reads::oneRead>& rawUnpairedReads, const std::vector<mapper::reads::verboseSeedChain>& alignedUnpairedReads, double insertSize_mean, double insertSize_sd, std::string outputDirectory, std::string longReadsMode);
	std::set<std::string> getCompletelyDefinedHLAAlleles(std::string locus);

	void fill_loci_2_exons();

	void simulateOneIndividual(std::string outputDirectory, double insertSize_mean, double insertSize_sd, bool novelIntronExonRecombinats, bool withError = true);
	
	static void read_true_types(std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>>& forReturn, std::string file);
	static void read_inferred_types(std::string sampleID, std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>>& forReturn, std::string file);
	
	static std::map<std::string, std::pair<int, int>> evaluate_HLA_types(std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> truth, std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferred);

	bool intervalOverlapsWithGenes(int first, int second) const;
};

} /* namespace hla */

#endif /* HLA_HLATYPER_H_ */
