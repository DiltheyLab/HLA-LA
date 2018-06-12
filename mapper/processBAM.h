/*
 * processBAM.h
 *
 *  Created on: 21.09.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_PROCESSBAM_H_
#define MAPPER_PROCESSBAM_H_

#include <string>
#include <vector>
#include <set>
#include <map>
#include <utility>

#include "reads/verboseSeedChain.h"
#include "reads/seedChain.h"
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "reads/PRGContigBAMAlignment.h"
#include "reads/protoSeeds.h"
#include "aligner/extensionAligner.h"
#include "aligner/statistics.h"

#include "../simulator/trueReadLevels.h"

#include "../hla/HLATyper.h"

#include <boost/math/distributions/normal.hpp>

namespace mapper {

class oneInterestingInterval
{
public:
	int PRGid;
	std::string refID;
	int start_1based;
	int stop_1based;

	oneInterestingInterval(int PRGid_, std::string refID_, int start_, int stop_) : PRGid(PRGid_), refID(refID_), start_1based(start_), stop_1based(stop_)
	{

	}
};

class processBAM {
protected:
	std::string graphDir;
	BamTools::BamReader R;

	BamTools::BamReader R1;
	BamTools::BamReader R2;

	Graph* g;
	std::vector<std::string> g_level_names;
	
	aligner::extensionAligner* eA;
	std::vector<bool> inGraphGapStretch;

	std::set<std::string> regionIDs;
	std::map<std::string, std::string> extendedReferenceGenomeSequences;
	std::map<std::string, std::string> PRGonlyReferenceGenomeSequences;

	std::map<int, std::vector<int>> extendedReferenceGenome_levelTranslation;

	std::vector<std::map<int, int>> graphLevel_2_underlyingSequencePositions;

	std::map<int, std::string> PRGid_2_BAMid;
	std::map<std::string, int> BAMid_2_PRGid;

	std::map<std::string, std::vector<oneInterestingInterval>> interestingIntervals;

	std::map<std::string, std::vector<oneInterestingInterval>>::iterator interestingIntervals_iterator;
	unsigned int interestingIntervals_i;

	BamTools::BamAlignment currentAlignment;

	void updateBAMSelector();
	void updateBAMSelectors();

	void _loadMapping(int ID);

	static bool transformBAMreadToInternalAlignment(const std::string& referenceSequence, const std::vector<int>& reference2level, int reference2level_offset_0based, const BamTools::BamAlignment& al, const std::string& queryBases, const std::string& qualitiesString, reads::PRGContigBAMAlignment& forReturn);

	reads::verboseSeedChainPair alignOneReadPair(const reads::protoSeeds& protoSeed, const boost::math::normal& rnd_InsertSize, double max_insertsize_penalty_log, simulator::trueReadLevels* trueReadLevels, aligner::statistics* statisticsStore = 0, std::vector<reads::verboseSeedChainPair>* allFoundAlignments_forReturn = 0) const;
	reads::verboseSeedChain alignOneLongRead(const reads::protoSeeds& protoSeed, const simulator::trueReadLevels* trueReadLevels, aligner::statistics* statisticsStore, std::string longReadMode) const;

	bool paranoid;

	std::set<Edge*> edges_for_deletion;

	reads::verboseSeedChain alignment2Chain(const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al, const std::string& complete_read_sequence, const std::string& complete_read_qualities, bool verbose = false) const;
	std::pair<int, int> alignment_get_startstop_PRGcoordinates(const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al) const;
	bool alignment_touches_gene(const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al, const std::vector<std::string>& levels) const;
	std::set<std::string> alignment_get_segments(const std::tuple<std::string, int,  BamTools::BamAlignment, int>& al, const std::vector<std::string>& levels) const;

	static void cleanInitialAlignment(std::vector<int>& graph_aligned_levels, std::string& graph_aligned, std::string& sequence_aligned);
	static void restrictInitialAlignmentToNoGapAreas(std::vector<int>& graph_aligned_levels, std::string& graph_aligned, std::string& sequence_aligned, int& graphSpace_sequence_aligned_startInRaw, int& graphSpace_sequence_aligned_stopInRaw, const std::vector<bool>& inGraphGapStretch);

	void initBAM(std::string BAM, bool extendedReferenceGenome, bool multiSampleMapping);
	void initBAMs(std::string BAM1, std::string BAM2, bool extendedReferenceGenome);

	std::string _currentBAM;
	std::string _currentBAM1;
	std::string _currentBAM2;
	bool _currentBAM_isExtendedReferenceGenome;

	void protoSeedStatistics(const std::map<std::string, reads::protoSeeds>& seeds, aligner::statistics* statisticsStore = 0, std::string longReadsMode = "") const;
	void reduceNonGeneSeeds(std::map<std::string, reads::protoSeeds>& seeds) const;
	void sortChainsInSeeds(std::map<std::string, reads::protoSeeds>& seeds) const;

	unsigned int threads_for_HLAtyping;
	
	std::pair<double, double> calculateInsertSizeFromHistogram(const std::map<int, double>& IS_combined_counts, bool verbose = true);

	void alignReads_postSeedExtraction(std::map<std::string, reads::protoSeeds>& seeds, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, hla::HLATyper* HLATyper, int threads);

	static int getAlignmentScore(const BamTools::BamAlignment& al);
	std::vector<std::string> getReadIDs();

public:
	processBAM(std::string graphDir_, unsigned int maxThreads = 1);
	virtual ~processBAM();

	std::pair<double, double> estimateInsertSize(std::string BAM, bool extendedReferenceGenome);
	std::pair<double, double> estimateInsertSize_noGraph(std::string BAM);


	void alignReads(std::string BAM, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, bool extendedReferenceGenome, hla::HLATyper* HLATyper = 0, int threads = 1);
	void alignReads_and_inferHLA(std::string BAM, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, bool extendedReferenceGenome, hla::HLATyper* HLATyper, int threads, std::string longReadsMode);
	void alignReads2(std::string BAM1, std::string BAM2, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, bool extendedReferenceGenome, hla::HLATyper* HLATyper = 0, int threads = 1);
	void alignReadsMulti(std::vector<std::string> BAMs, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, bool extendedReferenceGenome, hla::HLATyper* HLATyper = 0, int threads = 1);

	std::map<std::string, reads::protoSeeds> extractSeeds(long long maximumIncludedReads = -1, std::set<std::string> limitToReadIDs = std::set<std::string>(), std::string longReadMode = "");
	std::map<std::string, reads::protoSeeds> extractSeeds2(std::set<std::string> limitToReadIDs = std::set<std::string>(), std::string longReadMode = "");

	size_t alignReads_postSeedExtraction_andStoreInto(std::map<std::string, reads::protoSeeds>& seeds, simulator::trueReadLevels* trueReadLevels, double insertSize_mean, double insertSize_sd, std::string outputDirectory, const hla::HLATyper* HLATyper, int threads, std::vector<std::vector<int>>& bases_per_level_perThread, std::vector<std::vector<mapper::reads::oneReadPair>>& HLA_raw_reads_perThread, std::vector<std::vector<mapper::reads::verboseSeedChainPair>>& HLA_alignments_perThread, aligner::statistics* statisticsStore);
	size_t alignReadsUnpaired_postSeedExtraction_andStoreInto(std::map<std::string, reads::protoSeeds>& seeds, simulator::trueReadLevels* trueReadLevels, std::string outputDirectory, const hla::HLATyper* HLATyper, int threads, std::vector<std::vector<int>>& bases_per_level_perThread, std::vector<std::vector<mapper::reads::oneRead>>& HLA_raw_reads_perThread, std::vector<std::vector<mapper::reads::verboseSeedChain>>& HLA_alignments_perThread, aligner::statistics* statisticsStore, std::string longReadMode);

	std::map<std::string, reads::protoSeeds> extractSeeds_2BAMs(long long maximumIncludedReads = -1);
	void extractSeedsInto(std::map<std::string, reads::protoSeeds>& seeds, int whichReader);

	static bool PRGContigAlignment2Seed(Graph *g, const reads::PRGContigBAMAlignment& PRGcontigAlignment, bool paranoid, reads::verboseSeedChain& graphSeed, reads::verboseSeedChain& sequenceSeed, const std::vector<bool>& inGraphGapStretch);

	static void assignMappingQualities(reads::verboseSeedChainPair& forReturn, const std::vector<std::pair<unsigned int, unsigned int>> combinations_indices, const std::vector<double>& combinations_LL, const std::pair<double, unsigned int>& combinations_max, const std::vector<reads::verboseSeedChain>& read1_extendedChains, const std::vector<reads::verboseSeedChain>& read2_extendedChains);
	static void assignMappingQualities_unpaired(reads::verboseSeedChain& forReturn, const std::vector<double>& LLs, const std::pair<double, unsigned int>& maxLL, const std::vector<reads::verboseSeedChain>& read1_extendedChains);

	Graph* getGraph();
};

} /* namespace mapper */

#endif /* MAPPER_PROCESSBAM_H_ */
