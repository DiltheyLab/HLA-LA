/*
 * HLATyper.cpp
 *
 *  Created on: 28.10.2015
 *      Author: AlexanderDilthey
 */

#include "assert.h"
#include <unordered_map>


#include "HLATyper.h"
#include "oneExonPosition.h"

#include "../mapper/aligner/alignerBase.h"
#include "../Utilities.h"
#include <cmath>
#include <utility>

#include <omp.h>
#include <boost/math/distributions/chi_squared.hpp>

// int max_mismatches_perRead = 2;
//double min_alignmentFraction_OK = 0.96; // measures all alignment positions but graph AND sequence gaps, separately for both reads
//double min_oneRead_weightedCharactersOK = 0.995; // one read, mismatches downweighted by quality
// double min_bothReads_weightedCharactersOK = 0.985; // both reads, mismatches downweighted by quality
// double min_bothReads_weightedCharactersOK = 0.95; // todo reinstate
double min_bothReads_weightedCharactersOK = 0.0;

double minimumMappingQuality = 0.0;
double minimumPerPositionMappingQuality = 0.7;
bool veryConservativeReadLikelihoods = true;

bool globalVerbose = false;

namespace hla {

HLATyper::HLATyper(Graph* g_, std::string graphDir_, std::string simulations_qualityMatrixFile_) : g(g_), graphDir(graphDir_), simulations_qualityMatrixFile(simulations_qualityMatrixFile_) {
	bool allLoci = true;
	if(allLoci)
	{
		loci_for_HLAtyping = {"A", "B", "C", "DQA1", "DQB1", "DRB1", "DPA1", "DPB1", "DRA", "DRB3", "DRB4", "E", "F", "G", "H", "K", "V"};
	}
	else
	{
		loci_for_HLAtyping = {"A", "B", "C", "DQA1", "DQB1", "DRB1"};
	}
	
	fill_loci_2_exons();	
	files_in_graphDir = Utilities::filesInDirectory(graphDir+"/PRG");

	for(unsigned int locusI = 0; locusI < loci_for_HLAtyping.size(); locusI++)
	{
		std::string locus = loci_for_HLAtyping.at(locusI);
		assert(loci_2_exons.at(locus).size());
		for(unsigned int exonI = 0; exonI < loci_2_exons.at(locus).size(); exonI++)
		{
			std::string exonID = loci_2_exons.at(locus).at(exonI);			
			std::string exonFile = find_file_for_exon(locus, exonID);
			if(! Utilities::fileReadable(exonFile))
			{
				std::cerr << "HLATypeInference(..): Locus " << locus << ", exon " << exonID << ": Can't read file " << exonFile << "\n";
			}
		}
	}
	
	threshold_reportColumn_forPresenceOfUnaccountedAlleles_minCoverage = 30;
	threshold_reportColumn_forPresenceOfUnaccountedAlleles_minAlleleFraction = 0.2;
	highCoverage_filter_alleles = false;
	highCoverage_minCoverage = 100; 
	highCoverage_minAlleleFreq = 0.2;

	filterFirst20 = true;
	filterFirst20N = 20;
	filterFirst20MinProp = 0.1;

	longReads_filterStrand = true;
	longReads_filterStrand_minAlleleCoverage = 100;
	longReads_filterStrand_minStrandFreq = 0.1;




	// translate location IDs to graph levels
	graphLoci = Graph::readGraphLoci(graphDir);

	for(unsigned int i = 0; i < graphLoci.size(); i++)
	{
		std::string locusID = graphLoci.at(i);
		assert(graphLocus_2_levels.count(locusID) == 0);
		graphLocus_2_levels[locusID] = i;
	}

	simulations_read_length = 101;
	simulations_haploidCoverage = 15;
	if(simulations_qualityMatrixFile != "")
	{
		rS = new simulator::readSimulator(simulations_qualityMatrixFile, simulations_read_length);
	}
	else
	{
		rS = 0;
	}

	std::ifstream segmentsStream;
	std::string segmentsFileName = graphDir + "/PRG/segments.txt";
	segmentsStream.open(segmentsFileName.c_str());
	assert(segmentsStream.is_open());
	std::string line;
	while(segmentsStream.good())
	{
		std::getline(segmentsStream, line);
		Utilities::eraseNL(line);
		if(line.length() > 0)
		{
			std::vector<std::string> split_by_underscore = Utilities::split(line, "_");
			if(split_by_underscore.at(1) != "gene")
			{
				continue;
			}

			std::string thisLocus = split_by_underscore.at(2);
			graphGenes.insert(thisLocus);

			std::string filePath = graphDir + "/PRG/" + line;

			std::string segmentID = split_by_underscore.at(4);
			if(split_by_underscore.size() > 5)
			{
				segmentID = segmentID + split_by_underscore.at(5);
			}
			
			std::string type;
			if(segmentID.find("intron") != std::string::npos)
			{
				type = "intron";
			}
			else
			{
				type = "exon";
				assert(segmentID.find("exon") != std::string::npos);
			}

			if(type == "intron")
			{
				segments_per_graphGene[thisLocus].first.push_back(segmentID);
			}
			else
			{
				segments_per_graphGene[thisLocus].second.push_back(segmentID);
			}
			segments_per_graphGene_all[thisLocus].push_back(segmentID);

			if(sequences_per_graphGeneSegment[thisLocus].count(segmentID) != 0)
			{
				std::cerr << "! (sequences_per_graphGeneSegment[thisLocus].count(segmentID) == 0)" << "\n";
				std::cerr << thisLocus << "\n";
				std::cerr << segmentID << "\n";
				std::cerr << line << "\n";
				std::cerr << std::flush;
			}
			assert(sequences_per_graphGeneSegment[thisLocus].count(segmentID) == 0);

			std::ifstream fileInputStream;
			fileInputStream.open(filePath.c_str());
			assert(fileInputStream.is_open());
			std::vector<std::string> file_lines;
			while(fileInputStream.good())
			{
				std::string line;
				std::getline(fileInputStream, line);
				Utilities::eraseNL(line);
				file_lines.push_back(line);
			}
			fileInputStream.close();

			std::string firstLine = file_lines.at(0);
			std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
			assert(firstLine_fields.at(0) == "IndividualID");

			std::vector<std::string> graph_level_names(firstLine_fields.begin() + 1, firstLine_fields.end());

			if(graphgene_levelBoundaries.count(thisLocus) == 0)
			{
				graphgene_levelBoundaries[thisLocus] = make_pair(-1, -1);
			}
			for(unsigned int lI = 0; lI < graph_level_names.size(); lI++)
			{
				std::string levelName = graph_level_names.at(lI);
				assert(graphLocus_2_levels.count(levelName));
				int levelI = graphLocus_2_levels.at(levelName);
				if((graphgene_levelBoundaries.at(thisLocus).first == -1) || (levelI < graphgene_levelBoundaries.at(thisLocus).first))
				{
					graphgene_levelBoundaries.at(thisLocus).first = levelI;
				}
				if((graphgene_levelBoundaries.at(thisLocus).second == -1) || (levelI > graphgene_levelBoundaries.at(thisLocus).second))
				{
					graphgene_levelBoundaries.at(thisLocus).second = levelI;
				}
			}
			levelNames_per_graphGeneSegment[thisLocus][segmentID] = graph_level_names;

			for(unsigned int lI = 1; lI < file_lines.size(); lI++)
			{
				if(file_lines.at(lI).length())
				{
					std::vector<std::string> line_fields = Utilities::split(file_lines.at(lI), " ");
					assert(line_fields.size() == firstLine_fields.size());
					std::string HLA_type = line_fields.at(0);
					std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

					sequences_per_graphGeneSegment[thisLocus][segmentID][HLA_type] = line_alleles;
				}
			}
		}
	}

	std::cout << "HLATyper: Information on " << segments_per_graphGene.size() << " total graph genes\n";
	for(std::set<std::string>::iterator graphGeneIt = graphGenes.begin(); graphGeneIt != graphGenes.end(); graphGeneIt++)
	{
		std::string gene = *graphGeneIt;
		std::cout << " - " << gene << "\n";
		std::cout << "\t\t" << segments_per_graphGene.at(gene).first.size() << " introns." << "\n";
		std::cout << "\t\t" << segments_per_graphGene.at(gene).second.size() << " exons." << "\n";

		std::set<std::string> complete_genomic_types = get_complete_genomic_types_per_gene(gene);
		std::cout << "\t\t" << complete_genomic_types.size() << " complete genomic types." << "\n";

		std::set<std::string> complete_exonic_types = get_complete_exonic_types_per_gene(gene);
		std::cout << "\t\t" << complete_exonic_types.size() << " complete exonic types." << "\n";

		assert(segments_per_graphGene_all.at(gene).size() ==
				(
						segments_per_graphGene.at(gene).first.size() +
						segments_per_graphGene.at(gene).second.size()
				)
			);
	}


	std::vector<Interval<std::string>> gene_intervals;
	for(std::map<std::string, std::pair<int, int>>::iterator geneBoundaryIt = graphgene_levelBoundaries.begin(); geneBoundaryIt != graphgene_levelBoundaries.end(); geneBoundaryIt++)
	{
		std::string gene = geneBoundaryIt->first;
		std::pair<int, int> B = geneBoundaryIt->second;
		assert(B.first < B.second);

		Interval<std::string> intervalForTree(B.first, B.second, gene);
		gene_intervals.push_back(intervalForTree);
	}

	if(graphDir.find("PRG_MHC_GRCh38_withIMGT") != std::string::npos)
	{
		Interval<std::string> intervalForTree(6261147, 6274591, "DRB5_notInGraph");
		gene_intervals.push_back(intervalForTree);
	}

	interestingLevels = new IntervalTree<std::string>(gene_intervals);
	// assert(graphLoci.size() == (g->NodesPerLevel.size() - 1));
	
	std::cout << find_file_for_exon("A", "exon_2") << "\n";

}

bool HLATyper::intervalOverlapsWithGenes(int first, int second) const
{
	assert(first >= 0);
	assert(second >= 0);
	assert(second >= first);

	std::vector<Interval<std::string>> found_intervals;
	interestingLevels->findOverlapping(first, second, found_intervals);
	return (found_intervals.size() != 0);
}

std::string HLATyper::intervalOverlapsWithWhichGenes(int first, int second) const
{
	assert(first >= 0);
	assert(second >= 0);
	assert(second >= first);

	std::vector<Interval<std::string>> found_intervals;
	interestingLevels->findOverlapping(first, second, found_intervals);

	assert(found_intervals.size() != 0);
	if(!(found_intervals.size() == 1))
	{
		std::cerr << "Multiple overlapping genes\n";
		std::cerr << first << "\n";
		std::cerr << second << "\n" << std::flush;
		for(auto oneI : found_intervals)
		{
			std::cerr << oneI.value << "\n" << std::flush;
		}
	}
	assert(found_intervals.size() == 1);

	return found_intervals.at(0).value;
}

std::set<std::string> HLATyper::get_complete_exonic_types_per_gene(std::string gene)
{
	assert(segments_per_graphGene.count(gene));

	if(segments_per_graphGene.at(gene).second.size() > 0)
	{
		std::set<std::string> forReturn;

		const std::vector<std::string>& segmentIDs = segments_per_graphGene.at(gene).second;
		std::map<std::string, int> counter;
		
		for(unsigned int segmentI = 0; segmentI < segmentIDs.size(); segmentI++)
		{
			std::string segmentID = segmentIDs.at(segmentI);
			for(std::map<std::string, std::vector<std::string>>::iterator sequenceIterator = sequences_per_graphGeneSegment.at(gene).at(segmentID).begin(); sequenceIterator != sequences_per_graphGeneSegment.at(gene).at(segmentID).end(); sequenceIterator++)
			{
				std::string sequenceID = sequenceIterator->first;
				if(counter.count(sequenceID) == 0)
				{
					counter[sequenceID] = 0;
				}
				counter.at(sequenceID)++;
			}
		}
		
		for(std::map<std::string, int>::iterator sequenceIt = counter.begin(); sequenceIt != counter.end(); sequenceIt++)
		{
			std::string sequenceID = sequenceIt->first;
			int count = sequenceIt->second;
			if(count == (int)segmentIDs.size())
			{
				if(sequenceID.find("*") != std::string::npos)
				{
					forReturn.insert(sequenceID);
				}
			}
		}
		
		for(std::set<std::string>::iterator typeIt = forReturn.begin(); typeIt != forReturn.end(); typeIt++)
		{
			std::string type = *typeIt;

			for(std::map<std::string, std::map<std::string, std::vector<std::string>>>::iterator segmentIt = sequences_per_graphGeneSegment.at(gene).begin(); segmentIt != sequences_per_graphGeneSegment.at(gene).end(); segmentIt++)
			{
				std::string segmentID = segmentIt->first;

				if(segmentID.find("exon") != std::string::npos)
				{
					if(!(sequences_per_graphGeneSegment.at(gene).at(segmentID).count(type)))
					{
						std::cerr << "! (sequences_per_graphGeneSegment.at(gene).at(segmentID).count(type))" << "\n";
						std::cerr << gene << "\n";
						std::cerr << segmentID << "\n";
						std::cerr << type << "\n";
					}
					assert(sequences_per_graphGeneSegment.at(gene).at(segmentID).count(type));
				}
			}
		}

		return forReturn;
	}
	else
	{
		return get_complete_genomic_types_per_gene(gene);
	}
}

std::string HLATyper::removeLocusAndStar(std::string HLAType)
{
	size_t pos = HLAType.find("*");
	if(pos == std::string::npos)
	{
		return HLAType;
	}
	else
	{
		std::string forReturn;
		forReturn = HLAType.substr(pos + 1);
		// std::cerr << HLAType << " " << forReturn << "\n" << std::flush;
		return forReturn;
	}
}

std::set<std::string> HLATyper::get_complete_genomic_types_per_gene(std::string gene)
{
	assert(segments_per_graphGene.count(gene));
	std::set<std::string> forReturn;
	bool firstAddToReturn = true;
	for(unsigned int isExon = 0; isExon <= 1; isExon++)
	{
		const std::vector<std::string>& segmentIDs = (isExon == 0) ? segments_per_graphGene.at(gene).first : segments_per_graphGene.at(gene).second;
		for(unsigned int segmentI = 0; segmentI < segmentIDs.size(); segmentI++)
		{
			std::string segmentID = segmentIDs.at(segmentI);

			for(std::map<std::string, std::vector<std::string>>::iterator sequenceIterator = sequences_per_graphGeneSegment.at(gene).at(segmentID).begin(); sequenceIterator != sequences_per_graphGeneSegment.at(gene).at(segmentID).end(); sequenceIterator++)
			{
				std::string sequenceID = sequenceIterator->first;
				if(firstAddToReturn)
				{
					assert(forReturn.count(sequenceID) == 0);
					if(sequenceID.find("*") != std::string::npos)
					{					
						forReturn.insert(sequenceID);
					}
				}
				else
				{
					if(forReturn.count(sequenceID) == 0)
					{
						forReturn.erase(sequenceID);
					}
				}
			}

			if(sequences_per_graphGeneSegment.at(gene).at(segmentID).size() > 0)
			{
				firstAddToReturn = false;
			}
		}
	}

	for(std::set<std::string>::iterator typeIt = forReturn.begin(); typeIt != forReturn.end(); typeIt++)
	{
		std::string type = *typeIt;

		for(std::map<std::string, std::map<std::string, std::vector<std::string>>>::iterator segmentIt = sequences_per_graphGeneSegment.at(gene).begin(); segmentIt != sequences_per_graphGeneSegment.at(gene).end(); segmentIt++)
		{
			std::string segmentID = segmentIt->first;

			assert(sequences_per_graphGeneSegment.at(gene).at(segmentID).count(type));
		}
	}

	return forReturn;
}

std::map<std::string, std::pair<int, int>> HLATyper::evaluate_HLA_types(std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>> truth, std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>> inferred)
{
	std::map<std::string, std::pair<int, int>> forReturn;
	
	for(std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>>::iterator individualIt = inferred.begin(); individualIt != inferred.end(); individualIt++)
	{
		std::string individualID = individualIt->first;
		for(std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>::iterator locusIt = inferred.at(individualID).begin(); locusIt != inferred.at(individualID).end(); locusIt++)
		{
			std::string locus = locusIt->first;
			
			std::set<std::string> inferredAllele_1 = inferred.at(individualID).at(locus).first;
			std::set<std::string> inferredAllele_2 = inferred.at(individualID).at(locus).second;
			
			if(truth.count(individualID) && truth.at(individualID).count(locus))
			{
				std::string trueAllele_1 = truth.at(individualID).at(locus).first;
				std::string trueAllele_2 = truth.at(individualID).at(locus).second;
				trueAllele_1 = removeLocusAndStar(trueAllele_1);						
				trueAllele_2 = removeLocusAndStar(trueAllele_2);						
				
				std::set<std::string> trueAlleles_1; trueAlleles_1.insert(trueAllele_1);
				std::set<std::string> trueAlleles_2; trueAlleles_2.insert(trueAllele_2);
				
				int max_correctAlleles;
				int max_correctAlleles_comparedAlleles;
				for(int invertOrder = 0; invertOrder <= 1; invertOrder++)
				{
					const std::set<std::string>& inferredAllele_compare_1 = (invertOrder == 0) ? inferredAllele_1 : inferredAllele_2;
					const std::set<std::string>& inferredAllele_compare_2 = (invertOrder == 0) ? inferredAllele_2 : inferredAllele_1;
					
					int this_correctAlleles = 0;
					int this_comparedAlleles = 0;
					
					if(inferredAllele_compare_1.size() && trueAlleles_1.size())
					{
						for(std::set<std::string>::const_iterator inferredAlleleIt = inferredAllele_compare_1.begin(); inferredAlleleIt != inferredAllele_compare_1.end(); inferredAlleleIt++)
						{
							std::string inferredAllele = removeLocusAndStar(*inferredAlleleIt);						
							bool secondBreak = false;
							
							for(std::set<std::string>::const_iterator trueAlleleIt = trueAlleles_1.begin(); trueAlleleIt != trueAlleles_1.end(); trueAlleleIt++)
							{
								std::string trueAllele = *trueAlleleIt;
								if(alleles_compatible(inferredAllele, trueAllele))
								{
									this_correctAlleles++;
									secondBreak = true;
									break;
								}
							}
							
							if(secondBreak)
							{
								break;
							}
						}
						this_comparedAlleles++;
					}

					if(inferredAllele_compare_2.size() && trueAlleles_2.size())
					{
						for(std::set<std::string>::const_iterator inferredAlleleIt = inferredAllele_compare_2.begin(); inferredAlleleIt != inferredAllele_compare_2.end(); inferredAlleleIt++)
						{
							std::string inferredAllele = removeLocusAndStar(*inferredAlleleIt);
							bool secondBreak = false;
							
							for(std::set<std::string>::const_iterator trueAlleleIt = trueAlleles_2.begin(); trueAlleleIt != trueAlleles_2.end(); trueAlleleIt++)
							{
								std::string trueAllele = *trueAlleleIt;
								if(alleles_compatible(inferredAllele, trueAllele))
								{
									this_correctAlleles++;
									secondBreak = true;
									break;
								}
							}
							
							if(secondBreak)
							{
								break;
							}							
						}
						this_comparedAlleles++;
					}

					if((invertOrder == 0) || (this_correctAlleles > max_correctAlleles))
					{
						max_correctAlleles = this_correctAlleles;
						max_correctAlleles_comparedAlleles = this_comparedAlleles;
					}
				}
				
				assert(max_correctAlleles <= max_correctAlleles_comparedAlleles);
				
				if(forReturn.count(locus) == 0)
				{
					forReturn[locus].first = 0;
					forReturn[locus].second = 0;
				}
				
				forReturn.at(locus).first += max_correctAlleles_comparedAlleles;
				forReturn.at(locus).second += max_correctAlleles;
			}
		}
	}
	
	std::cout << "HLATyper::evaluate_HLA_types(..) summary:\n";
	for(std::map<std::string, std::pair<int, int>>::iterator locusIt = forReturn.begin(); locusIt != forReturn.end(); locusIt++)
	{
		std::string locus = locusIt->first;
		std::cout << "\t" << locus << "\n";
		assert(locusIt->second.first > 0);
		double rateOK = locusIt->second.second / locusIt->second.first;
		std::cout << "\t\t" << locusIt->second.first << " evaluated." << "\n";
		std::cout << "\t\t" << locusIt->second.second << " correct." << "\n";
		std::cout << "\t\t"  << (rateOK*100) << "%" << "\n";
	}
	return forReturn;
}
HLATyper::~HLATyper() {
	// TODO Auto-generated destructor stub
}

bool HLATyper::alleles_compatible(std::string inferredAllele, std::string trueAllele)
{
	assert(inferredAllele.length());
	assert(trueAllele.length());
	
	if((trueAllele.substr(trueAllele.length() - 1, 1) == "G") || (trueAllele.substr(trueAllele.length() - 1, 1) == "g"))
	{
		trueAllele = trueAllele.substr(0, trueAllele.length() - 1);
	}

	if(!(inferredAllele.find(":") != std::string::npos))
	{
		std::cerr << "!(inferredAllele.find() != std::string::npos)" << "\n";
		std::cerr << "\t" << "inferredAllele" << ": " << inferredAllele << "\n" << std::flush;
	}
	assert(inferredAllele.find(":") != std::string::npos);
	if(trueAllele.find(":") == std::string::npos)
	{
		if(!(trueAllele.length() == 4))
		{
			std::cerr << "! (trueAllele.length() == 4)" << "\n";
			std::cerr << trueAllele << "\n";
			std::cerr << std::flush;
		}
		assert(trueAllele.length() == 4);
		trueAllele = trueAllele.substr(0, 2) + ":" + trueAllele.substr(2);
	}
	
	std::vector<std::string> inferred_groups = Utilities::split(inferredAllele, ":");
	std::vector<std::string> true_groups = Utilities::split(trueAllele, ":");
	
	assert(inferred_groups.size() >= 2);
	assert(true_groups.size() >= 2);
	
	if(inferred_groups.size() < true_groups.size())
	{
		return false;
	}
	else
	{
		for(unsigned int groupI = 0; groupI < true_groups.size(); groupI++)
		{
			if(inferred_groups.at(groupI) != true_groups.at(groupI))
			{
				return false;
			}
		}
		
		return true;
	}
}

void HLATyper::read_inferred_types(std::string sampleID, std::map<std::string, std::map<std::string, std::pair<std::set<std::string>, std::set<std::string>>>>& forReturn, std::string file)
{
	assert(forReturn.count(sampleID) == 0);
	
	std::ifstream fileInputStream;
	fileInputStream.open(file.c_str());
	assert(fileInputStream.is_open());
	assert(fileInputStream.good());
	std::string line;
	std::getline(fileInputStream, line);
	Utilities::eraseNL(line);
	std::vector<std::string> header_line_fields = Utilities::split(line, "\t");
	assert(header_line_fields.at(0) == "Locus");
	assert(header_line_fields.at(1) == "Chromosome");
	assert(header_line_fields.at(2) == "Allele");
	assert(header_line_fields.size() > 1);
	while(fileInputStream.good())
	{
		std::getline(fileInputStream, line);
		Utilities::eraseNL(line);
		if(line.length() == 0)
			continue;
		
		std::vector<std::string> line_fields = Utilities::split(line, "\t");
		assert(line_fields.size() == header_line_fields.size());
				
		std::string locus = line_fields.at(0);
		int chromosome = Utilities::StrtoI(line_fields.at(1));
		assert((chromosome == 1) || (chromosome == 2));
		std::vector<std::string> alleles = Utilities::split(line_fields.at(2), ";");
		assert(alleles.size() > 0);
		
		if(chromosome == 1)
		{
			assert(forReturn[sampleID][locus].first.size() == 0);
			forReturn[sampleID][locus].first.insert(alleles.begin(), alleles.end());
		}
		else
		{
			assert(forReturn[sampleID][locus].second.size() == 0);
			forReturn[sampleID][locus].second.insert(alleles.begin(), alleles.end());			
		}
	}	
}

void HLATyper::read_true_types(std::map<std::string, std::map<std::string, std::pair<std::string, std::string>>>& forReturn, std::string file)
{
	bool removeHLA = true;
	bool useSpace = false;
	
	std::ifstream fileInputStream;
	fileInputStream.open(file.c_str());
	assert(fileInputStream.is_open());
	assert(fileInputStream.good());
	std::string line;
	std::getline(fileInputStream, line);
	Utilities::eraseNL(line);
	std::vector<std::string> header_line_fields = Utilities::split(line, "\t");
	if(header_line_fields.size() == 1)
	{
		useSpace = true;
		header_line_fields = Utilities::split(line, " ");
		assert(header_line_fields.size() > 1);
	}
	assert(header_line_fields.at(0) == "IndividualID");
	assert(header_line_fields.size() > 1);
	while(fileInputStream.good())
	{
		std::getline(fileInputStream, line);
		Utilities::eraseNL(line);
		if(line.length() == 0)
			continue;
		
		std::vector<std::string> line_fields = (useSpace) ? Utilities::split(line, " ") : Utilities::split(line, "\t");
		assert(line_fields.size() == header_line_fields.size());
		
		std::string individualID = line_fields.at(0);
		assert(forReturn.count(individualID) == 0);
		for(unsigned int i = 1; i < header_line_fields.size(); i++)
		{
			std::string locus = header_line_fields.at(i);
			if(removeHLA)
			{
				if(locus.length() > 4)
				{
					if(locus.substr(0, 4) == "HLA-")
					{
						locus = locus.substr(4);
					}
				}
				if(locus.length() > 3)
				{
					if(locus.substr(0, 3) == "HLA")
					{
						locus = locus.substr(3);
					}
				}				
			}
			std::vector<std::string> alleles = Utilities::split(line_fields.at(i), "/");
			assert(alleles.size() == 2);
			forReturn[individualID][locus].first = alleles.at(0);
			forReturn[individualID][locus].second = alleles.at(1);
		}
	}
	fileInputStream.close();
}

void HLATyper::simulateOneIndividual(std::string outputDirectory, double insertSize_mean, double insertSize_sd, bool novelIntronExonRecombinats, bool withError)
{
	assert(rS != 0);

	std::string FASTQ1 = outputDirectory + "/R_1.fq";
	std::string FASTQ2 = outputDirectory + "/R_2.fq";
	std::string outputTypes = outputDirectory + "/HLAtypes.txt";
	std::string outputHaplotypes = outputDirectory + "/haplotypes.txt";

	std::string parametersFile = outputDirectory + "/parameters.txt";
	std::ofstream parametersStream;
	parametersStream.open(parametersFile.c_str());
	parametersStream << "graphDir: " << graphDir << "\n";
	parametersStream << "simulations_read_length: " << simulations_read_length << "\n";
	parametersStream << "insertSize_mean: " << insertSize_mean << "\n";
	parametersStream << "insertSize_sd: " << insertSize_sd << "\n";
	parametersStream << "simulations_haploidCoverage: " << simulations_haploidCoverage << "\n";
	parametersStream << "simulations_qualityMatrixFile: " << simulations_qualityMatrixFile << "\n";
	parametersStream << "withError: " << withError << "\n";
	parametersStream << "rS average error rates: " << rS->averageErrorRate_R1_R2().first << "\t" <<  rS->averageErrorRate_R1_R2().second <<  "\n";
	parametersStream.close();

	std::string fn_output_FASTQ_r1 = outputDirectory + "/R_1.fq";
	std::string fn_output_FASTQ_r2 = outputDirectory + "/R_2.fq";

	std::string fn_output_levels_r1 = outputDirectory + "/R_1.levels";
	std::string fn_output_levels_r2 = outputDirectory + "/R_2.levels";

	std::ofstream stream_FASTQ_r1;
	stream_FASTQ_r1.open(fn_output_FASTQ_r1.c_str());
	assert(stream_FASTQ_r1.is_open());

	std::ofstream stream_FASTQ_r2;
	stream_FASTQ_r2.open(fn_output_FASTQ_r2.c_str());
	assert(stream_FASTQ_r2.is_open());

	std::ofstream stream_levels_r1;
	stream_levels_r1.open(fn_output_levels_r1.c_str());
	assert(stream_levels_r1.is_open());

	std::ofstream stream_levels_r2;
	stream_levels_r2.open(fn_output_levels_r2.c_str());
	assert(stream_levels_r2.is_open());

	auto print_one_readPair = [&] (const simulator::oneReadPair& rP) -> void
	{
		stream_FASTQ_r1 << "@" << rP.reads.first.name << "\n";
		stream_FASTQ_r1 << rP.reads.first.sequence << "\n";
		stream_FASTQ_r1 << "+" << "\n";
		stream_FASTQ_r1 << rP.reads.first.quality << "\n";

		stream_FASTQ_r2 << "@" << rP.reads.second.name << "\n";
		stream_FASTQ_r2 << rP.reads.second.sequence << "\n";
		stream_FASTQ_r2 << "+" << "\n";
		stream_FASTQ_r2 << rP.reads.second.quality << "\n";

		stream_levels_r1 << "@" << rP.reads.first.name << "\n";
		stream_levels_r1 << Utilities::join(Utilities::ItoStr(rP.reads.first.coordinates_edgePath), " ") << "\n";
		stream_levels_r1 << rP.reads.first.underlyingEdgeLabels << "\n";
		stream_levels_r1 << Utilities::join(Utilities::ItoStr(rP.reads.first.fullAlignment_coordinates_edgePath), " ") << "\n";
		stream_levels_r1 << rP.reads.first.fullAlignment_underlyingEdgeLabels << "\n";
		stream_levels_r1 << rP.reads.first.fullAlignment_sequence << "\n";

		stream_levels_r2 << "@" << rP.reads.second.name << "\n";
		stream_levels_r2 << Utilities::join(Utilities::ItoStr(rP.reads.second.coordinates_edgePath), " ") << "\n";
		stream_levels_r2 << rP.reads.second.underlyingEdgeLabels << "\n";
		stream_levels_r2 << Utilities::join(Utilities::ItoStr(rP.reads.second.fullAlignment_coordinates_edgePath), " ") << "\n";
		stream_levels_r2 << rP.reads.second.fullAlignment_underlyingEdgeLabels << "\n";
		stream_levels_r2 << rP.reads.second.fullAlignment_sequence << "\n";
	};

	std::vector<std::string> processedGenes;
	std::map<std::string, std::vector<std::string>> underlyingTypes;
	std::map<std::string, std::vector<std::vector<std::string>>> underlyingHaplotypes;
	std::map<std::string, std::vector<std::string>> underlyingHaplotypes_levelIDs;

	for(std::set<std::string>::iterator geneIt = graphGenes.begin(); geneIt != graphGenes.end(); geneIt++)
	{
		std::string gene = *geneIt;

		std::vector<std::string> chooseTypeFrom;
		if(novelIntronExonRecombinats)
		{
			std::set<std::string> exonTypes = get_complete_exonic_types_per_gene(gene);
			chooseTypeFrom.insert(chooseTypeFrom.end(), exonTypes.begin(), exonTypes.end());
		}
		else
		{
			std::set<std::string> genomicTypes = get_complete_genomic_types_per_gene(gene);
			chooseTypeFrom.insert(chooseTypeFrom.end(), genomicTypes.begin(), genomicTypes.end());
		}
		assert(chooseTypeFrom.size() > 0);
		std::vector<std::string> selectedTypes;
		selectedTypes.push_back(Utilities::choose_uniformly_from_vector(chooseTypeFrom));
		selectedTypes.push_back(Utilities::choose_uniformly_from_vector(chooseTypeFrom));

		std::vector<std::vector<std::string>> haplotypes;
		std::vector<std::string> levelIDs;
		haplotypes.resize(selectedTypes.size());
		assert(haplotypes.size());

		const std::vector<std::string>& gene_segments = segments_per_graphGene_all.at(gene);
		assert(gene_segments.size());
		for(unsigned int segmentI = 0; segmentI < gene_segments.size(); segmentI++)
		{
			std::string segmentID = gene_segments.at(segmentI);
			bool isExon;
			if(segmentID.find("intron") != std::string::npos)
			{
				isExon = false;
			}
			else
			{
				if(!(segmentID.find("exon") != std::string::npos))
				{
					std::cerr << "! t(segmentID.find(exon) != std::string::npos)" << "\n";
					std::cerr << gene << "\n";
					std::cerr << segmentID << "\n";
					std::cerr << std::flush;
				}
				assert(segmentID.find("exon") != std::string::npos);
				isExon = true;
			}

			for(unsigned int haplotypeI = 0; haplotypeI < selectedTypes.size(); haplotypeI++)
			{
				std::string chosenType;
				if((! isExon) && (novelIntronExonRecombinats))
				{
					std::vector<std::string> choose_intron_type_from;
					for(std::map<std::string, std::vector<std::string>>::iterator typeIt = sequences_per_graphGeneSegment.at(gene).at(segmentID).begin(); typeIt != sequences_per_graphGeneSegment.at(gene).at(segmentID).end(); typeIt++)
					{
						std::string type = typeIt->first;
						choose_intron_type_from.push_back(type);
					}
					chosenType = Utilities::choose_uniformly_from_vector(choose_intron_type_from);
				}
				else
				{
					chosenType = selectedTypes.at(haplotypeI);
				}

				haplotypes.at(haplotypeI).insert(haplotypes.at(haplotypeI).end(), sequences_per_graphGeneSegment.at(gene).at(segmentID).at(chosenType).begin(), sequences_per_graphGeneSegment.at(gene).at(segmentID).at(chosenType).end());

				if(haplotypeI == 0)
				{
					levelIDs.insert(levelIDs.end(), levelNames_per_graphGeneSegment.at(gene).at(segmentID).begin(), levelNames_per_graphGeneSegment.at(gene).at(segmentID).end());
				}
			}
		}

		for(unsigned int haplotypeI = 0; haplotypeI < haplotypes.size(); haplotypeI++)
		{
			const std::vector<std::string>& haplotype = haplotypes.at(haplotypeI);
			assert(haplotype.size() == levelIDs.size());

			for(unsigned int characterI = 0; characterI < haplotype.size(); characterI++)
			{
				assert(haplotype.at(characterI).size() == 1);
			}
			std::string haplotypeString = Utilities::join(haplotype, "");

			std::cout << "Simulate for " << gene << ", string length " << haplotypeString.length() << "\n" << std::flush;
			
			std::vector<simulator::oneReadPair> reads = rS->simulate_paired_reads_from_string(
				haplotypeString,
				simulations_haploidCoverage,
				insertSize_mean,
				insertSize_sd,
				! withError,
				"PRG_" + gene + "_HAPLO_" + Utilities::ItoStr(haplotypeI)
			);

			for(size_t i = 0; i < reads.size(); i++)
			{
				for(int secondRead = 0; secondRead <= 1; secondRead++)
				{
					simulator::oneRead& R = (secondRead == 0) ? reads.at(i).reads.first : reads.at(i).reads.second;
					assert(R.coordinates_string.size() == R.sequence.length());
					for(unsigned int coordinateI = 0; coordinateI < R.coordinates_string.size(); coordinateI++)
					{
						int coordinate =  R.coordinates_string.at(coordinateI);
						if(coordinate != -1)
						{
							std::string levelID = levelIDs.at(coordinate);
							int correct_coordinate = graphLocus_2_levels.at(levelID);
							R.coordinates_string.at(coordinateI) = correct_coordinate;
						}
					}
				}
				print_one_readPair(reads.at(i));
			}
		}

		underlyingTypes[gene] = selectedTypes;
		underlyingHaplotypes[gene] = haplotypes;
		underlyingHaplotypes_levelIDs[gene] = levelIDs;
		processedGenes.push_back(gene);
	}

	std::ofstream trueTypesStream;
	trueTypesStream.open(outputTypes.c_str());
	assert(trueTypesStream.is_open());
	trueTypesStream << "IndividualID" << "\t" << Utilities::join(processedGenes, "\t") << "\n";
	std::vector<std::string> output_fields_trueTypes;
	output_fields_trueTypes.push_back(outputDirectory);
	for(unsigned int geneI = 0; geneI < processedGenes.size(); geneI++)
	{
		output_fields_trueTypes.push_back(Utilities::join(underlyingTypes.at(processedGenes.at(geneI)), "/"));
	}
	trueTypesStream << Utilities::join(output_fields_trueTypes, "\t") << "\n";
	trueTypesStream.close();

	std::ofstream trueHaplotypesStream;
	trueHaplotypesStream.open(outputHaplotypes.c_str());
	assert(trueHaplotypesStream.is_open());
	std::vector<std::string> haplotypes_header_fields;
	for(unsigned int geneI = 0; geneI < processedGenes.size(); geneI++)
	{
		std::string gene = processedGenes.at(geneI);
		haplotypes_header_fields.push_back("LevelIDs_" + gene);
		haplotypes_header_fields.push_back("Haplotype_" + gene);
	}
	trueHaplotypesStream << "IndividualID" << "\t" << Utilities::join(haplotypes_header_fields, "\t") << "\n";
	std::vector<std::string> output_fields_trueHaplotypes;
	output_fields_trueHaplotypes.push_back(outputDirectory);
	for(unsigned int geneI = 0; geneI < processedGenes.size(); geneI++)
	{
		output_fields_trueHaplotypes.push_back(Utilities::join(underlyingHaplotypes_levelIDs.at(processedGenes.at(geneI)), ";"));

		std::vector<std::string> output_per_haplotype;
		for(unsigned int haplotypeI = 0; haplotypeI < underlyingHaplotypes.at(processedGenes.at(geneI)).size(); haplotypeI++)
		{
			output_per_haplotype.push_back( Utilities::join(underlyingHaplotypes.at(processedGenes.at(geneI)).at(haplotypeI), ";") );
		}
		output_fields_trueHaplotypes.push_back(Utilities::join(output_per_haplotype, "/"));
	}
	trueHaplotypesStream << Utilities::join(output_fields_trueHaplotypes, "\t") << "\n";
	trueHaplotypesStream.close();

}

// td::string alignedReads_file, std::string graphDir, std::string sampleName, bool restrictToFullHaplotypes, std::string& forReturn_lociString, std::string& forReturn_starting_haplotype_1, std::string& forReturn_starting_haplotype_2,
void HLATyper::HLATypeInference(const std::vector<mapper::reads::oneReadPair>& rawPairedReads, const std::vector<mapper::reads::verboseSeedChainPair>& alignedPairedReads, const std::vector<mapper::reads::oneRead>& rawUnpairedReads, const std::vector<mapper::reads::verboseSeedChain>& alignedUnpairedReads, double insertSize_mean, double insertSize_sd, std::string outputDirectory, std::string longReadsMode)
{
	double insertionP = 0.001;
	double deletionP = 0.001;

	if(longReadsMode.length())
	{
		assert(rawPairedReads.size() == 0);
		insertionP = 0.075;
		deletionP = 0.075;
		
		highCoverage_filter_alleles = true;
		highCoverage_minCoverage = 1; 
		highCoverage_minAlleleFreq = 0.15;
	}
	else
	{
		assert(rawUnpairedReads.size() == 0);
	} 

	double log_likelihood_insertion = log(insertionP);
	double log_likelihood_insertion_actualAllele = log_likelihood_insertion + log(1.0/4.0);

	double log_likelihood_deletion = log(deletionP);
	double log_likelihood_nonInsertion = log(1 - insertionP);
	double log_likelihood_nonDeletion = log(1 - deletionP);
	double log_likelihood_match_mismatch = log(1 - insertionP - deletionP);

	assert(
			((rawPairedReads.size() == 0) && (alignedPairedReads.size() == 0)) ||
			((rawUnpairedReads.size() == 0) && (alignedUnpairedReads.size() == 0))
	);

	assert(
			(rawPairedReads.size() == alignedPairedReads.size()) &&
			(rawUnpairedReads.size() == alignedUnpairedReads.size())
	);

	assert(
			(rawPairedReads.size() > 0) ||
			(rawUnpairedReads.size() > 0)
	);

	bool longUnpairedReads = false;
	if(rawUnpairedReads.size() > 0)
	{
		longUnpairedReads = true;
	}

	std::string graph = graphDir + "/PRG/graph.txt";
	assert(Utilities::fileReadable(graph));

	// define locus -> exon

	int HLATypeInference_totalBases_used = 0;
	int HLAtypeInference_totalColumns = 0;

	// load reads
	std::cout << Utilities::timestamp() << "HLATyper::HLATypeInference(..): Load reads.\n" << std::flush;
	const std::vector<mapper::reads::verboseSeedChainPair>& alignments_paired = alignedPairedReads;
	const std::vector<mapper::reads::verboseSeedChain>& alignments_unpaired = alignedUnpairedReads;
	const std::vector<mapper::reads::oneReadPair>& alignments_originalReads_paired = rawPairedReads;
	const std::vector<mapper::reads::oneRead>& alignments_originalReads_unpaired = rawUnpairedReads;

	std::cout << Utilities::timestamp() << "HLATypeInference(..): Load reads -- done. Have " << alignments_paired.size() << " paired reads and " << alignments_unpaired.size() << " unpaired reads, long unpaired reads: " << longUnpairedReads << "\n" << std::flush;

	int k_for_kMer_index = 31;
	std::unordered_map<std::string, int> kMer_counts;
	// build kMer index of reads
	{
		auto add_sequence_to_index = [&](const std::string& S) -> void {
			std::vector<std::string> kMers = Utilities::partitionStringIntokMers(S, k_for_kMer_index);
			for(const auto& kMer : kMers)
			{
				std::string kMerKey = HLATyper::kMer_canonical_representation(kMer);
				if(kMer_counts.count(kMerKey) == 0)
				{
					kMer_counts[kMerKey] = 0;
				}
				kMer_counts.at(kMerKey)++;
			}
		};

		for(auto rP : alignments_originalReads_paired)
		{
			add_sequence_to_index(rP.reads.first.sequence);
			add_sequence_to_index(rP.reads.second.sequence);
		}

		for(auto uP : alignments_originalReads_unpaired)
		{
			add_sequence_to_index(uP.sequence);
		}

		std::cout << "Have " << kMer_counts.size() << " " << k_for_kMer_index << "-mers in index.\n" << std::flush;
	}

	// read alignment statistics

	int minAlignmentLength_unpaired = 1000;
	int alignmentStats_strandsValid = 0;
	int alignments_paired_perfect = 0;
	int alignments_paired_oneReadPerfect = 0;
	int alignments_unpaired_perfect = 0;
	int alignmentStats_strandsValid_and_distanceOK = 0;
	int alignments_unpaired_longEnough = 0;
	std::vector<double> alignmentStats_strandsValid_distances;
	double alignmentStats_paired_fractionOK_sum = 0;
	double alignmentStats_unpaired_fractionOK_sum = 0;

	for(unsigned int alignmentI = 0; alignmentI < alignments_paired.size(); alignmentI++)
	{
		const mapper::reads::verboseSeedChainPair& alignedReadPair = alignments_paired.at(alignmentI);

		if(mapper::aligner::alignerBase::alignedReadPair_strandsValid(alignedReadPair.chains.first, alignedReadPair.chains.second))
		{
			alignmentStats_strandsValid++;
			double pairsDistance = mapper::aligner::alignerBase::alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair);
			alignmentStats_strandsValid_distances.push_back(pairsDistance);
			if(abs(pairsDistance - insertSize_mean) <= (5 * insertSize_sd))
			{
				alignmentStats_strandsValid_and_distanceOK++;
			}
		}

		double fractionOK_1 = alignmentFractionOK(alignedReadPair.chains.first);
		double fractionOK_2 = alignmentFractionOK(alignedReadPair.chains.second);

		if(fractionOK_1 == 1)
		{
			alignments_paired_perfect++;
		}
		if(fractionOK_2 == 1)
		{
			alignments_paired_perfect++;
		}
		if((fractionOK_1 == 1) || (fractionOK_2 == 1))
		{
			alignments_paired_oneReadPerfect++;
		}

		alignmentStats_paired_fractionOK_sum += fractionOK_1;
		alignmentStats_paired_fractionOK_sum += fractionOK_2;
	}


	for(unsigned int alignmentI = 0; alignmentI < alignments_unpaired.size(); alignmentI++)
	{

		const mapper::reads::verboseSeedChain& alignedRead = alignments_unpaired.at(alignmentI);

		double fractionOK_1 = alignmentFractionOK(alignedRead);

		if(alignedRead.graph_aligned.size() >= minAlignmentLength_unpaired)
		{
			alignments_unpaired_longEnough++;
		}	
		
		if(fractionOK_1 == 1)
		{
			alignments_unpaired_perfect++;
		}

		alignmentStats_unpaired_fractionOK_sum += fractionOK_1;
	}

	std::pair<double, double> alignmentStats_distance_meanMedian = meanMedian(alignmentStats_strandsValid_distances);
	double alignmentStats_paired_fractionOK_avg = (alignments_paired.size() > 0) ? (alignmentStats_paired_fractionOK_sum / (2.0* (double)alignments_paired.size())) : 0;
	double alignmentStats_unpaired_fractionOK_avg = (alignments_unpaired.size() > 0) ? (alignmentStats_unpaired_fractionOK_sum / ((double)alignments_unpaired.size())) : 0;

	if(! Utilities::directoryExists(outputDirectory))
	{
		Utilities::makeDir(outputDirectory);
	}

	std::ofstream summaryStatisticsStream;
	std::string summaryStatisticsFilename = outputDirectory + "/" + "summaryStatistics.txt";
	summaryStatisticsStream.open(summaryStatisticsFilename.c_str());
	assert(summaryStatisticsStream.is_open());
	summaryStatisticsStream << "\nRead alignment statistics:\n";
	summaryStatisticsStream << "\t - Total number (paired) alignments:                 " << alignments_paired.size() << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs with strands OK:                  " << alignmentStats_strandsValid << " (" << printPerc(alignmentStats_strandsValid, alignments_paired.size()) << "%)\n";
	summaryStatisticsStream << "\t\t - Alignment pairs with strands OK && distance OK:   " << alignmentStats_strandsValid_and_distanceOK << " (" << printPerc(alignmentStats_strandsValid_and_distanceOK, alignments_paired.size()) << "%)\n";
	summaryStatisticsStream << "\t\t - Alignment pairs with strands OK, mean distance:   " << alignmentStats_distance_meanMedian.first << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs with strands OK, median distance: " << alignmentStats_distance_meanMedian.second << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs, average fraction alignment OK:   " << alignmentStats_paired_fractionOK_avg << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs, at least one alignment perfect:   " << alignments_paired_oneReadPerfect << "\n";
	summaryStatisticsStream << "\t\t - Single alignments, perfect (total):   " << alignments_paired_perfect << " (" << alignments_paired.size()*2 << ")\n";
	summaryStatisticsStream << "\t - Total number (unpaired) alignments:                 " << alignments_unpaired.size() << "\n";
	summaryStatisticsStream << "\t\t - Alignment pairs, average fraction alignment OK:   " << alignmentStats_unpaired_fractionOK_avg << "\n";
	summaryStatisticsStream << "\t\t - Single alignments, perfect (total):   " << alignments_unpaired_perfect << " (" << alignments_unpaired.size()*2 << ")\n";
	summaryStatisticsStream << "\t\t - Alignments with length >= " << minAlignmentLength_unpaired << ":   " << alignments_unpaired_longEnough << "\n";
	summaryStatisticsStream.close();

	std::string outputFN_bestGuess = outputDirectory + "/R1_bestguess.txt";
	std::ofstream bestGuess_outputStream;
	bestGuess_outputStream.open(outputFN_bestGuess.c_str());
	assert(bestGuess_outputStream.is_open());
	std::string fieldName_unaccountedAlleles = "NColumns_UnaccountedAllele_fGT" + Utilities::DtoStr(threshold_reportColumn_forPresenceOfUnaccountedAlleles_minAlleleFraction);

	bestGuess_outputStream << "Locus" << "\t" << "Chromosome" << "\t" << "Allele" << "\t" << "Q1" << "\t" << "Q2" << "\t" << "AverageCoverage" << "\t" << "CoverageFirstDecile" << "\t" << "MinimumCoverage" << "\t" << "proportionkMersCovered" << "\t" << "LocusAvgColumnError" << "\t" << fieldName_unaccountedAlleles << "\n";

	std::string outputFN_bestGuess_G = outputDirectory + "/R1_bestguess_G.txt";
	std::ofstream bestGuess_G_outputStream;
	bestGuess_G_outputStream.open(outputFN_bestGuess_G.c_str());
	assert(bestGuess_G_outputStream.is_open());
	bestGuess_G_outputStream << "Locus" << "\t" << "Chromosome" << "\t" << "Allele" << "\t" << "Q1" << "\t" << "Q2" << "\t" << "AverageCoverage" << "\t" << "CoverageFirstDecile" << "\t" << "MinimumCoverage" << "\t" << "proportionkMersCovered" << "\t" << "LocusAvgColumnError" << "\t" << fieldName_unaccountedAlleles << "\t" << "perfectG" << "\n";

	std::string outputFN_matchesPerReadHistogram = outputDirectory + "/histogram_matchesPerRead.txt";
	std::ofstream matchesPerReadHistogram_outputStream;
	matchesPerReadHistogram_outputStream.open(outputFN_matchesPerReadHistogram.c_str());
	assert(matchesPerReadHistogram_outputStream.is_open());
	matchesPerReadHistogram_outputStream << "Locus" << "\t" << "Level" << "Value" << "\n";

	std::vector<std::string> forReturn_starting_haplotype_1_vec;
	std::vector<std::string> forReturn_starting_haplotype_2_vec;
 
	for(unsigned int locusI = 0; locusI < loci_for_HLAtyping.size(); locusI++)
	{
		std::string locus = loci_for_HLAtyping.at(locusI);
		std::set<std::string> utilized_reads;
		
		//globalVerbose = (locus == "C");

		int HLATypeInference_thisLocus_bases_used = 0;

		std::string outputFN_allPairs = outputDirectory + "/R1_PP_"+locus+"_pairs.txt";
		std::string outputFN_columnError = outputDirectory + "/R1_columnIncompatibilities_"+locus+".txt";

		std::cout << Utilities::timestamp() << "HLATypeInference(..): Making inference for " << locus << "\n" << std::flush;

		std::vector<int> combined_exon_sequences_graphLevels;
		std::vector<int> combined_exon_sequences_graphLevels_individualExon;
		std::vector<int> combined_exon_sequences_graphLevels_individualExonPosition;
		int combined_exon_sequences_graphLevels_min = -1;
		int combined_exon_sequences_graphLevels_max = -1;
		
		std::vector<std::string> combined_exon_sequences_locusIDs;
		std::map<std::string, std::string> combined_exon_sequences;

		std::set<std::string> completeDefinedTypes;
		bool restrictToFullHaplotypes = false;
		if(restrictToFullHaplotypes)
		{
			completeDefinedTypes = getCompletelyDefinedHLAAlleles(locus);
			std::cout << Utilities::timestamp() << "HLATypeInference(..): restrictToFullHaplotypes in force, restrict to " << completeDefinedTypes.size() << " types.\n" << std::flush;
		}

		int thisLocus_totalColumns = 0;

		std::map<int, int> exon_lengths;

		assert(loci_2_exons.at(locus).size());
		for(unsigned int exonI = 0; exonI < loci_2_exons.at(locus).size(); exonI++)
		{
			std::string exonID = loci_2_exons.at(locus).at(exonI);
			std::cout << Utilities::timestamp() << "\tLocus" << locus << ", exon " << exonID << "\n" << std::flush;

			std::string exonFile = find_file_for_exon(locus, exonID);
			if(! Utilities::fileReadable(exonFile))
			{
				std::cerr << "HLATypeInference(..): Locus " << locus << ", exon " << exonID << ": Can't read file " << exonFile << "\n";
			}
			assert(Utilities::fileReadable(exonFile));

			std::ifstream exonInputStream;
			exonInputStream.open(exonFile.c_str());
			assert(exonInputStream.is_open());
			std::vector<std::string> exon_lines;
			while(exonInputStream.good())
			{
				std::string line;
				std::getline(exonInputStream, line);
				Utilities::eraseNL(line);
				exon_lines.push_back(line);
			}
			exonInputStream.close();

			std::string firstLine = exon_lines.at(0);
			std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
			assert(firstLine_fields.at(0) == "IndividualID");

			std::vector<std::string> exon_level_names(firstLine_fields.begin() + 1, firstLine_fields.end());
			std::string first_graph_locusID = exon_level_names.front();
			std::string last_graph_locusID = exon_level_names.back();

			assert(graphLocus_2_levels.count(first_graph_locusID));
			assert(graphLocus_2_levels.count(last_graph_locusID));

			unsigned int first_graph_level = graphLocus_2_levels.at(first_graph_locusID);
			unsigned int last_graph_level = graphLocus_2_levels.at(last_graph_locusID);

			std::cout << Utilities::timestamp() << "\tLocus" << locus << ", exon " << exonID << ": from " << first_graph_locusID << " (" << first_graph_level << ") to " << last_graph_locusID << " (" << last_graph_level << ").\n" << std::flush;

			assert(last_graph_level > first_graph_level);
			unsigned int expected_allele_length = last_graph_level - first_graph_level + 1;
			if(!(exon_level_names.size() == expected_allele_length))
			{
				std::cerr << "For locus " << locus << " exon " << exonID << " (" << exonFile << "), we have a problem with expected graph length.\n";
				std::cerr << "\t" << "expected_allele_length" << ": " << expected_allele_length << "\n";
				std::cerr << "\t" << "expected_allele_length" << ": " << expected_allele_length << "\n";
				std::cerr << std::flush;
			}
			assert(exon_level_names.size() == expected_allele_length);

			thisLocus_totalColumns += exon_level_names.size();

			combined_exon_sequences_locusIDs.insert(combined_exon_sequences_locusIDs.end(), exon_level_names.begin(), exon_level_names.end());
			assert(expected_allele_length);
			for(unsigned int lI = 0; lI < expected_allele_length; lI++)
			{
				unsigned int graphLevel = first_graph_level + lI;
				assert(graphLocus_2_levels.at(exon_level_names.at(lI)) == graphLevel);
				combined_exon_sequences_graphLevels.push_back(graphLevel);
				combined_exon_sequences_graphLevels_individualExon.push_back(exonI);
				combined_exon_sequences_graphLevels_individualExonPosition.push_back(lI);

				if((combined_exon_sequences_graphLevels_min == -1) || (combined_exon_sequences_graphLevels_min > (int)graphLevel))
				{
					combined_exon_sequences_graphLevels_min = graphLevel;
				}
				if((combined_exon_sequences_graphLevels_max == -1) || (combined_exon_sequences_graphLevels_max < (int)graphLevel))
				{
					combined_exon_sequences_graphLevels_max = graphLevel;
				}
			}

			exon_lengths[exonI] = expected_allele_length;

			for(unsigned int lI = 1; lI < exon_lines.size(); lI++)
			{
				if(exon_lines.at(lI).length())
				{
					std::vector<std::string> line_fields = Utilities::split(exon_lines.at(lI), " ");
					assert(line_fields.size() == firstLine_fields.size());
					std::string HLA_type = line_fields.at(0);
					
					if(HLA_type.find(":") == std::string::npos)
					{
						continue;
					}
					
					if(restrictToFullHaplotypes)
					{
						if(completeDefinedTypes.count(HLA_type) == 0)
						{
							continue;
						}
					}
					std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

					std::string HLA_type_sequence = Utilities::join(line_alleles, "");

					if(exonI == 0)
					{
						assert(combined_exon_sequences.count(HLA_type) == 0);
						combined_exon_sequences[HLA_type] = HLA_type_sequence;
					}
					else
					{
						assert(combined_exon_sequences.count(HLA_type));
						combined_exon_sequences.at(HLA_type) += HLA_type_sequence;
					}
				}
			}
			assert(combined_exon_sequences.size() > 0);
		}

		std::map<int, unsigned int> graphLevel_2_exonPosition;
		std::map<int, unsigned int> graphLevel_2_exonPosition_individualExon;
		std::map<int, unsigned int> graphLevel_2_exonPosition_individualExonPosition;
		for(unsigned int pI = 0; pI < combined_exon_sequences_graphLevels.size(); pI++)
		{
			int graphLevel = combined_exon_sequences_graphLevels.at(pI);
			assert(graphLevel >= 0);
			graphLevel_2_exonPosition[graphLevel] = pI;
			graphLevel_2_exonPosition_individualExon[graphLevel] = combined_exon_sequences_graphLevels_individualExon.at(pI);
			graphLevel_2_exonPosition_individualExonPosition[graphLevel] = combined_exon_sequences_graphLevels_individualExonPosition.at(pI);
		}

		std::cout << Utilities::timestamp() << "Have collected " << combined_exon_sequences.size() << " sequences -- first level " << combined_exon_sequences_graphLevels.front() << ", last level " << combined_exon_sequences_graphLevels.back() << ".\n" << std::flush;

		HLAtypeInference_totalColumns += thisLocus_totalColumns;

		std::map<std::string, unsigned int> HLAtype_2_clusterID;
		std::vector<std::set<std::string>> HLAtype_clusters;
		std::map<std::string, unsigned int> sequence_2_cluster;
		std::vector<std::string> cluster_2_sequence;

		size_t sequenceL;
		for(std::map<std::string, std::string>::iterator HLAtypeIt = combined_exon_sequences.begin(); HLAtypeIt != combined_exon_sequences.end(); HLAtypeIt++)
		{
			std::string HLAtypeID = HLAtypeIt->first;
			std::string sequence = HLAtypeIt->second;

			if(HLAtypeIt == combined_exon_sequences.begin())
			{
				sequenceL = sequence.length();
			}
			else
			{
				assert(sequenceL == sequence.length());
			}

			if(sequence_2_cluster.count(sequence))
			{
				unsigned int cluster = sequence_2_cluster.at(sequence);
				HLAtype_clusters.at(cluster).insert(HLAtypeID);
				HLAtype_2_clusterID[HLAtypeID] = cluster;
			}
			else
			{
				std::set<std::string> newCluster;
				newCluster.insert(HLAtypeID);
				HLAtype_clusters.push_back(newCluster);
				unsigned int newClusterID = HLAtype_clusters.size() - 1;

				assert(sequence_2_cluster.count(sequence) == 0);
				sequence_2_cluster[sequence] = newClusterID;

				HLAtype_2_clusterID[HLAtypeID] = newClusterID;
			}
		}

		std::cout << Utilities::timestamp() << "Clustered into " << HLAtype_clusters.size() << " identical (over exons considered) clusters." << "\n" << std::flush;

		for(unsigned int clusterI = 0; clusterI < HLAtype_clusters.size(); clusterI++)
		{
			for(std::set<std::string>::iterator typeIt = HLAtype_clusters.at(clusterI).begin(); typeIt != HLAtype_clusters.at(clusterI).end(); typeIt++)
			{
				std::string HLAtype = *typeIt;
				assert(HLAtype_2_clusterID.at(HLAtype) == clusterI);

				if(typeIt == HLAtype_clusters.at(clusterI).begin())
				{
					std::string sequence = combined_exon_sequences.at(HLAtype);
					cluster_2_sequence.push_back(sequence);
				}
			}
		}

		// now transform reads into sequences specifying exon genotype value

		std::cout << Utilities::timestamp() << "Compute exon positions and specified genotypes from " << alignments_paired.size() << " paired and " << alignments_unpaired.size() << " unpaired reads \n" << std::flush;
		std::vector< std::vector<oneExonPosition> > exonPositions_fromReads;


		unsigned int readPairs_OK = 0;
		unsigned int readPairs_broken = 0;

		assert(combined_exon_sequences_graphLevels_min != -1);
		assert(combined_exon_sequences_graphLevels_max != -1);

		for(unsigned int readPairI = 0; readPairI < alignments_paired.size(); readPairI++)
		{
			const mapper::reads::oneReadPair& originalReadPair = alignments_originalReads_paired.at(readPairI);
			const mapper::reads::verboseSeedChainPair& alignedReadPair = alignments_paired.at(readPairI);

			std::vector<oneExonPosition> read1_exonPositions;
			std::vector<oneExonPosition> read2_exonPositions;

			oneReadAlignment_2_exonPositions_paired(alignedReadPair.chains.first, originalReadPair.reads.first, read1_exonPositions, alignedReadPair.chains.second, originalReadPair.reads.second, 1, combined_exon_sequences_graphLevels_min, combined_exon_sequences_graphLevels_max, graphLevel_2_exonPosition);
			oneReadAlignment_2_exonPositions_paired(alignedReadPair.chains.second, originalReadPair.reads.second, read2_exonPositions, alignedReadPair.chains.first, originalReadPair.reads.first, 2, combined_exon_sequences_graphLevels_min, combined_exon_sequences_graphLevels_max, graphLevel_2_exonPosition);

			if(!((alignedReadPair.chains.first.mapQ >= 0) && (alignedReadPair.chains.first.mapQ <= 1)))
			{
					std::cerr << "alignedReadPair.chains.first.mapQ" << ": " << alignedReadPair.chains.first.mapQ << "\n";
					std::cerr << std::flush; 
			}
			
			assert((alignedReadPair.chains.first.mapQ >= 0) && (alignedReadPair.chains.first.mapQ <= 1));
			double mapQ_thisAlignment = alignedReadPair.chains.first.mapQ;
			if(
					mapper::aligner::alignerBase::alignedReadPair_strandsValid(alignedReadPair) &&
					(abs(mapper::aligner::alignerBase::alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair) - insertSize_mean) <= (5 * insertSize_sd)) &&
					(mapQ_thisAlignment >= minimumMappingQuality) &&
					((alignmentWeightedOKFraction(originalReadPair.reads.first, alignedReadPair.chains.first) >= min_bothReads_weightedCharactersOK) && (alignmentWeightedOKFraction(originalReadPair.reads.second, alignedReadPair.chains.second) >= min_bothReads_weightedCharactersOK))
					)
			{
				// good
				// std::cout << "\t\t" << "readPair " << readPairI << ", pairing OK.\n" << std::flush;

				std::vector<oneExonPosition> thisRead_exonPositions = read1_exonPositions;
				thisRead_exonPositions.insert(thisRead_exonPositions.end(), read2_exonPositions.begin(), read2_exonPositions.end());

				if(thisRead_exonPositions.size() > 0)
				{
					thisRead_exonPositions = removeDoublePositionsFromRead(thisRead_exonPositions);
					exonPositions_fromReads.push_back(thisRead_exonPositions);
				}

				readPairs_OK++;

				matchesPerReadHistogram_outputStream << locus << "\t" << "read" << alignmentWeightedOKFraction(originalReadPair.reads.first, alignedReadPair.chains.first) << "\n";
				matchesPerReadHistogram_outputStream << locus << "\t" << "read" << alignmentWeightedOKFraction(originalReadPair.reads.second, alignedReadPair.chains.second) << "\n";

				matchesPerReadHistogram_outputStream << locus << "\t" << "readPair" << (alignmentWeightedOKFraction(originalReadPair.reads.first, alignedReadPair.chains.first)+ alignmentWeightedOKFraction(originalReadPair.reads.second, alignedReadPair.chains.second))/2.0 << "\n";
			}
			else
			{
				// bad

				// std::cout << "\t\t" << "readPair " << readPairI << "/" < < alignments.size() << ", pairing FAILED.\n" << std::flush;

				if(mapper::aligner::alignerBase::alignedReadPair_strandsValid(alignedReadPair) &&
				(abs(mapper::aligner::alignerBase::alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair) - insertSize_mean) <= (5 * insertSize_sd)))
				{
					std::cout << "REJECTED MAPQ " << alignedReadPair.chains.first.mapQ << " GENOMIC " << alignedReadPair.chains.first.mapQ << "\n" << std::flush;
				}

				if(read1_exonPositions.size() || read2_exonPositions.size())
				{
					if(globalVerbose)
					{
						if(alignedReadPair.readID.find("HLA-C") != std::string::npos)
						{
							std::cout << "Report " << alignedReadPair.readID << "\n";
							std::cout << "\t" << "mapper::aligner::alignerBase::alignedReadPair_strandsValid(alignedReadPair)" << ": " << mapper::aligner::alignerBase::alignedReadPair_strandsValid(alignedReadPair) << "\n";
							std::cout << "\t" << "mapper::aligner::alignerBase::alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair)" << ": " << mapper::aligner::alignerBase::alignedReadPair_pairsDistanceInGraphLevels(alignedReadPair) << "\n";
							std::cout << "\t\t" << "insertSize_mean" << ": " << insertSize_mean << "\n";
							std::cout << "\t\t" << "insertSize_sd" << ": " << insertSize_sd << "\n";
							std::cout << "\t" << "mapQ_thisAlignment" << ": " << mapQ_thisAlignment << "\n";
							std::cout << "\t\t" << "minimumMappingQuality" << ": " << minimumMappingQuality << "\n";
							std::cout << "\t" << "alignmentWeightedOKFraction(originalReadPair.reads.first, alignedReadPair.chains.first)" << ": " << alignmentWeightedOKFraction(originalReadPair.reads.first, alignedReadPair.chains.first) << "\n";
							std::cout << "\t" << "alignmentWeightedOKFraction(originalReadPair.reads.second, alignedReadPair.chains.second)" << ": " << alignmentWeightedOKFraction(originalReadPair.reads.second, alignedReadPair.chains.second) << "\n";
							std::cout << "\n" << std::flush;
						}
					}
				}
				readPairs_broken++;
			}
		}

		size_t reads_with_alignments_taken = 0;
		for(unsigned int readI = 0; readI < alignments_unpaired.size(); readI++)
		{
			const mapper::reads::oneRead& originalRead = alignments_originalReads_unpaired.at(readI);
			const mapper::reads::verboseSeedChain& alignedRead = alignments_unpaired.at(readI);

			std::vector<oneExonPosition> read_exonPositions;
			oneReadAlignment_2_exonPositions_unpaired(alignedRead, originalRead, read_exonPositions, combined_exon_sequences_graphLevels_min, combined_exon_sequences_graphLevels_max, graphLevel_2_exonPosition);

			double mapQ_thisAlignment = alignedRead.mapQ;
			if((mapQ_thisAlignment >= minimumMappingQuality) && (alignedRead.graph_aligned.size() >= minAlignmentLength_unpaired))
			{
				// good
				std::vector<oneExonPosition> thisRead_exonPositions = read_exonPositions;

				if(thisRead_exonPositions.size() > 0)
				{
					exonPositions_fromReads.push_back(read_exonPositions);
					reads_with_alignments_taken++;
				}
				readPairs_OK++;
			}
			else
			{
				// bad
				std::cout << "Rejected mapQ " << alignedRead.mapQ << " / alignment length " << alignedRead.graph_aligned.size() << "\n" << std::flush;

				readPairs_broken++;
			}
		}

		std::cout << Utilities::timestamp() << "Mapped reads to exons. " << readPairs_OK << " pairs OK (and long enough when in long-read mode), " << readPairs_broken << " pairs broken. Took alignments from " << reads_with_alignments_taken << "unpaired reads.\n" << std::flush;

		// filtering

		// Pileup of mapped reads
		std::set<size_t> ignore_exonPositions_fromReads;

		// filtering of top quality distribution
		std::set<std::string> ignore_readIDs;
		std::map<unsigned int, std::map<std::string, int>> perPosition_allele_counts_postFiltering;
		std::map<unsigned int, std::set<std::string>> perPosition_ignore_alleles;

		if(filterFirst20 && (longReadsMode.length() == 0))
		{
			std::map<unsigned int, int> perRead_kickedOut;
			std::map<unsigned int, int> perRead_kickedOut_robust;

			unsigned int n_relevant_positions = graphLevel_2_exonPosition_individualExon.size();

			std::map<unsigned int, std::vector<std::string>> perPosition_alleles;
			std::map<unsigned int, std::vector<double>> perPosition_weightedOK;
			std::map<unsigned int, std::vector<unsigned int>> perPosition_readI;

			for(unsigned int readI = 0; readI < exonPositions_fromReads.size(); readI++)
			{
				std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(readI);
				for(unsigned int positionI = 0; positionI < individualPositions.size(); positionI++)
				{
					oneExonPosition& onePositionSpecifier = individualPositions.at(positionI);
					assert((onePositionSpecifier.mapQ_position >= 0) && (onePositionSpecifier.mapQ_position <= 1));
					if(onePositionSpecifier.mapQ_position < minimumPerPositionMappingQuality)
					{
						continue;
					}

					int position = onePositionSpecifier.positionInExon;
					std::string allele = onePositionSpecifier.genotype;

					double completeRead_weightedCharactersOK = (onePositionSpecifier.thisRead_WeightedCharactersOK+onePositionSpecifier.pairedRead_WeightedCharactersOK)/2.0;

					perPosition_alleles[position].push_back(allele);
					perPosition_weightedOK[position].push_back(completeRead_weightedCharactersOK);
					perPosition_readI[position].push_back(readI);

				}
			}

			size_t total_considered_positions = 0;
			size_t withKickedOut_considered_positions = 0;
			size_t total_considered_alleles = 0;
			size_t kickedOut_considered_alleles = 0;
			for(auto position : perPosition_alleles)
			{
				int n_alleles = perPosition_alleles.at(position.first).size();
				if(n_alleles < filterFirst20N)
					continue;
				std::vector<unsigned int> allele_indices;
				allele_indices.reserve(n_alleles);
				for(unsigned int i = 0; i < n_alleles; i++)
				{
					allele_indices.push_back(i);
				}
				std::sort(allele_indices.begin(), allele_indices.end(),
						[&](unsigned int a, unsigned int b)
						{
							return (perPosition_weightedOK.at(position.first).at(a) < perPosition_weightedOK.at(position.first).at(b));
						}
				);
				std::reverse(allele_indices.begin(), allele_indices.end()); assert(perPosition_weightedOK.at(position.first).at(allele_indices.at(0)) >= perPosition_weightedOK.at(position.first).at(allele_indices.at(1)));

				std::map<std::string, int> alleles_first20;
				std::set<std::string> kickedOutAlleles;
				std::set<std::string> passedAlleles;
				for(unsigned int i = 0; i < filterFirst20N; i++)
				{
					unsigned int iI = allele_indices.at(i);
					std::string allele = perPosition_alleles.at(position.first).at(iI);
					if(alleles_first20.count(allele) == 0)
					{
						alleles_first20[allele] = 0;
					}
					alleles_first20.at(allele)++;
				}
				

				bool kickedOneOut = false;
				for(unsigned int i = 0; i < n_alleles; i++)
				{
					std::string allele = perPosition_alleles.at(position.first).at(i);
					unsigned int readI = perPosition_readI.at(position.first).at(i);
					int first20_alleleCount = 0;
					if(alleles_first20.count(allele))
					{
						first20_alleleCount = alleles_first20.at(allele);
					}
					double first20_prop = (double) first20_alleleCount / (double) filterFirst20;

					total_considered_alleles++;
					if(first20_prop < filterFirst20MinProp)
					{
						kickedOutAlleles.insert(allele);
						perPosition_ignore_alleles[position.first].insert(allele);
						if(perRead_kickedOut.count(readI) == 0)
						{
							perRead_kickedOut[readI] = 0;
						}
						perRead_kickedOut.at(readI)++;
						kickedOneOut = true;
						kickedOut_considered_alleles++;
					}
					else
					{
						passedAlleles.insert(allele);
					}
				}
				
				// find out how many alleles of a specific type have been kicked out
				std::map<std::string, int> allele_kickedOut_howMany;
				for(unsigned int i = 0; i < n_alleles; i++)
				{
					std::string allele = perPosition_alleles.at(position.first).at(i);
					if(kickedOutAlleles.count(allele))
					{
						if(allele_kickedOut_howMany.count(allele) == 0)
						{
							allele_kickedOut_howMany[allele] = 0;
						}
						allele_kickedOut_howMany.at(allele)++;
					}
				}				
				
				for(unsigned int i = 0; i < n_alleles; i++)
				{
					std::string allele = perPosition_alleles.at(position.first).at(i);
					unsigned int readI = perPosition_readI.at(position.first).at(i);
					if(allele_kickedOut_howMany.count(allele) && (allele_kickedOut_howMany.at(allele) >= 2))
					{
						if(perRead_kickedOut_robust.count(readI) == 0)
						{
							perRead_kickedOut_robust[readI] = 0;
						}
						perRead_kickedOut_robust.at(readI)++;						
					}
				}
				
				std::vector<std::string> kickedOutAlleles_vec(kickedOutAlleles.begin(), kickedOutAlleles.end());
				std::vector<std::string> passedAlleles_vec(passedAlleles.begin(), passedAlleles.end());
				
				if(kickedOutAlleles_vec.size() || 1)
				{
					int thisPosition_individualExon = combined_exon_sequences_graphLevels_individualExon.at(position.first);
					int thisPosition_individualExonPosition = combined_exon_sequences_graphLevels_individualExonPosition.at(position.first);
					// std::cout << "Position " << position.first << " [e" << thisPosition_individualExon << " " << thisPosition_individualExonPosition << "] allele QC; passed " << Utilities::join(passedAlleles_vec, ", ") << "; kicked out " << Utilities::join(kickedOutAlleles_vec, ", ") << "\n";
				}
				
				total_considered_positions++;
				if(kickedOneOut)
				{
					withKickedOut_considered_positions++;
				}
			}

			std::map<int, int> reads_kickedOutPositions_histogram;
			int n_reads_kickedOut = 0;
			for(auto readKickedOut : perRead_kickedOut)
			{
				if(readKickedOut.second > filterFirst20MinProp_limitKickOutPerRead)
				{
					unsigned int readI = readKickedOut.first;
					const oneExonPosition& onePositionSpecifier = exonPositions_fromReads.at(readI).at(0);
					n_reads_kickedOut++;
					//ignore_readIDs.insert(onePositionSpecifier.thisRead_ID);
					//ignore_readIDs.insert(onePositionSpecifier.pairedRead_ID);
				}
				
				if(reads_kickedOutPositions_histogram.count(readKickedOut.second) == 0)
				{
					reads_kickedOutPositions_histogram[readKickedOut.second] = 0;
				}
				reads_kickedOutPositions_histogram.at(readKickedOut.second)++;
			}
			
			
			std::map<int, int> reads_kickedOutPositions_histogram_robust;
			int n_reads_kickedOut_robust = 0;
			for(auto readKickedOut : perRead_kickedOut_robust)
			{
				if(readKickedOut.second > filterFirst20MinProp_limitKickOutPerRead)
				{
					unsigned int readI = readKickedOut.first;
					const oneExonPosition& onePositionSpecifier = exonPositions_fromReads.at(readI).at(0);
					n_reads_kickedOut_robust++;
					ignore_readIDs.insert(onePositionSpecifier.thisRead_ID);
					ignore_readIDs.insert(onePositionSpecifier.pairedRead_ID);
				}
				
				if(reads_kickedOutPositions_histogram_robust.count(readKickedOut.second) == 0)
				{
					reads_kickedOutPositions_histogram_robust[readKickedOut.second] = 0;
				}
				reads_kickedOutPositions_histogram_robust.at(readKickedOut.second)++;
			}
			
			
			

			std::cout << "Summary of first-" << filterFirst20N << " filtering; threshold " << filterFirst20MinProp << "\n";
			std::cout << "\t" << "Considered positions (with enough coverage): " << total_considered_positions << "\n";
			std::cout << "\t" << "Positions with removed alleles: " << withKickedOut_considered_positions << "\n";
			std::cout << "\t" << "Total considered alleles: " << total_considered_alleles << "\n";
			std::cout << "\t" << "Removed alleles: " << kickedOut_considered_alleles << "\n";
			std::cout << "\t" << "Total reads: " <<  exonPositions_fromReads.size() << "\n" << std::flush;
			std::cout << "\t" << "Reads with more than " << filterFirst20MinProp_limitKickOutPerRead << " positions: " << n_reads_kickedOut << "\n" << std::flush;
			for(auto histCat : reads_kickedOutPositions_histogram)
			{
				std::cout << "\t\t" << histCat.first << ": " << histCat.second << "\n" << std::flush;
			}
			std::cout << "\t" << "Reads with more than " << filterFirst20MinProp_limitKickOutPerRead << " robust positions: " << n_reads_kickedOut_robust << "\n" << std::flush;
			for(auto histCat : reads_kickedOutPositions_histogram_robust)
			{
				std::cout << "\t\t" << histCat.first << ": " << histCat.second << "\n" << std::flush;
			}			
		}

		// filtering of low-frequency allelles and store strand bias

		std::map<unsigned int, std::map<std::string, double>> perPosition_allele_minStrandFreqs;
		std::map<unsigned int, std::map<std::string, double>> perPosition_allele_read1Freq;
		{ 
			unsigned int n_relevant_positions = graphLevel_2_exonPosition_individualExon.size();
			size_t coverage_over_relevant_positions = 0;
			std::set<size_t> read_with_good_positions;

			std::map<unsigned int, std::map<std::string, int>> perPosition_allele_counts;
			std::map<unsigned int, std::map<std::string, std::pair<int, int>>> perPosition_allele_counts_byStrand;
			std::map<unsigned int, std::map<std::string, int>> perPosition_allele_counts_strand1;
			std::map<unsigned int, std::map<std::string, std::vector<double>>> perPosition_allele_count_readWeightedOK;
			
			std::vector<double> allUsedBases_weightedOK;
			for(unsigned int readI = 0; readI < exonPositions_fromReads.size(); readI++)
			{
				std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(readI);
				for(unsigned int positionI = 0; positionI < individualPositions.size(); positionI++)
				{
					oneExonPosition& onePositionSpecifier = individualPositions.at(positionI);
					if(ignore_readIDs.count(onePositionSpecifier.thisRead_ID))
						continue;

					assert((onePositionSpecifier.mapQ_position >= 0) && (onePositionSpecifier.mapQ_position <= 1));
					if(onePositionSpecifier.mapQ_position < minimumPerPositionMappingQuality)
					{
						continue;
					}

					int position = onePositionSpecifier.positionInExon;
					std::string allele = onePositionSpecifier.genotype;

					if(perPosition_ignore_alleles.count(position) && (perPosition_ignore_alleles.at(position).count(allele)))
					{
						continue;
					}

						

					if(perPosition_allele_counts[position].count(allele) == 0)
					{
						perPosition_allele_counts[position][allele] = 0;
						perPosition_allele_counts_byStrand[position][allele] = make_pair(0, 0);
						perPosition_allele_counts_strand1[position][allele] = 0;
					}

					perPosition_allele_counts.at(position).at(allele)++;
					if(onePositionSpecifier.reverse)
					{
						perPosition_allele_counts_byStrand.at(position).at(allele).second++;
					}
					else
					{
						perPosition_allele_counts_byStrand.at(position).at(allele).first++;
					}
					
					if(onePositionSpecifier.fromFirstRead)
					{
						perPosition_allele_counts_strand1.at(position).at(allele)++;
					}
					
					double completeRead_weightedCharactersOK = (onePositionSpecifier.thisRead_WeightedCharactersOK+onePositionSpecifier.pairedRead_WeightedCharactersOK)/2.0;
					perPosition_allele_count_readWeightedOK[position][allele].push_back(completeRead_weightedCharactersOK);
					allUsedBases_weightedOK.push_back(completeRead_weightedCharactersOK);
				}
			}

			int positions_sufficientCoverage = 0;
			int positions_with_kickedOutAlleles = 0;
			int positions_with_kickedOutAlleles_strandFilter = 0;

			int kickedOutAlleles = 0;
			int strandFilter_allelesEnoughCoverage = 0;
			int kickedOutAlleles_strandFilter = 0;

			if(highCoverage_filter_alleles)
				std::cout << "Locus " << locus << ", highCoverage_filter_alleles active, evaluate " << perPosition_allele_counts.size() << " positions.\n" << std::flush;
			for(auto position : perPosition_allele_counts)
			{
				int count_position = 0;
				for(auto allele : perPosition_allele_counts.at(position.first))
				{
					count_position += allele.second;
				}
				if(count_position >= highCoverage_minCoverage)
				{
					positions_sufficientCoverage++;
					bool kickedOutAtLeastOneAllele = false;
					for(auto allele : perPosition_allele_counts.at(position.first))
					{
						double aF = (double)allele.second / (double)count_position;
						if((aF < highCoverage_minAlleleFreq) && highCoverage_filter_alleles)
						{
							perPosition_ignore_alleles[position.first].insert(allele.first);
							kickedOutAlleles += allele.second;
							kickedOutAtLeastOneAllele = true;
						}
						else
						{
							perPosition_allele_counts_postFiltering[position.first][allele.first] = allele.second;
						}
					}
					if(kickedOutAtLeastOneAllele)
					{
						positions_with_kickedOutAlleles++;
					}
				}

				bool kickedOutAtLeastOneAllele_sF = false;

				for(auto allele : perPosition_allele_counts_byStrand.at(position.first))
				{
					int totalCount = perPosition_allele_counts_byStrand.at(position.first).at(allele.first).first + perPosition_allele_counts_byStrand.at(position.first).at(allele.first).second;
					assert(totalCount == perPosition_allele_counts.at(position.first).at(allele.first));
					int minStrandCount =
							(perPosition_allele_counts_byStrand.at(position.first).at(allele.first).first < perPosition_allele_counts_byStrand.at(position.first).at(allele.first).second) ?
									perPosition_allele_counts_byStrand.at(position.first).at(allele.first).first : perPosition_allele_counts_byStrand.at(position.first).at(allele.first).second
					;
					int read1Count = perPosition_allele_counts_strand1.at(position.first).count(allele.first) ? perPosition_allele_counts_strand1.at(position.first).at(allele.first) : 0;

					double minStrandFreq = (double)minStrandCount/(double)totalCount;
					double read1Freq = (double)read1Count/(double)totalCount;
					perPosition_allele_minStrandFreqs[position.first][allele.first] = minStrandFreq;
					perPosition_allele_read1Freq[position.first][allele.first] = read1Freq;
					assert(minStrandFreq <= 0.5);
					if((longReadsMode.length() > 0) && longReads_filterStrand && (totalCount >= longReads_filterStrand_minAlleleCoverage))
					{
						strandFilter_allelesEnoughCoverage++;
						if(minStrandFreq < longReads_filterStrand_minStrandFreq)
						{
							perPosition_ignore_alleles[position.first].insert(allele.first);
							kickedOutAlleles_strandFilter++;
							kickedOutAtLeastOneAllele_sF = true;
						}
					}
					if(kickedOutAtLeastOneAllele_sF)
					{
						positions_with_kickedOutAlleles_strandFilter++;
					}
				}
			}

			if(highCoverage_filter_alleles)
			{
				std::cout << "\tPositions with coverage high enough for filtering: " << positions_sufficientCoverage << "\n";
				std::cout << "\tPositions with at least one allele removed with aF <= " << highCoverage_minAlleleFreq << ": " << positions_with_kickedOutAlleles << "\n";
				std::cout << "\tTotal number of (individual-read) alleles removed: " << kickedOutAlleles << "\n";
				std::cout << "\n" << std::flush;
			}

			if((longReadsMode.length() > 0) && longReads_filterStrand)
			{ 
				std::cout << "Locus " << locus << ", long-reads strand-filter (by allele) active.\n" << std::flush;
				std::cout << "\tAlleles with enough coverage (" << longReads_filterStrand_minAlleleCoverage << "): " << strandFilter_allelesEnoughCoverage << "\n";
				std::cout << "\tAlleles removed: " << kickedOutAlleles_strandFilter << "\n";
				std::cout << "\tPositions with alleles removed: " << positions_with_kickedOutAlleles_strandFilter << "\n";
				std::cout << "\n" << std::flush;
			}
		}

		// Pileup of mapped reads

		std::map<int, std::map<int, std::vector<oneExonPosition> > > pileUpPerPosition;
		for(unsigned int positionSpecifierI = 0; positionSpecifierI < exonPositions_fromReads.size(); positionSpecifierI++)
		{
			std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(positionSpecifierI);
			for(unsigned int positionI = 0; positionI < individualPositions.size(); positionI++)
			{
				oneExonPosition& onePositionSpecifier = individualPositions.at(positionI);

				assert((onePositionSpecifier.mapQ_position >= 0) && (onePositionSpecifier.mapQ_position <= 1));
				if(onePositionSpecifier.mapQ_position < minimumPerPositionMappingQuality)
				{
					if(globalVerbose)
					{
						if( onePositionSpecifier.thisRead_ID.find("HLA-C") != std::string::npos)
						{
							std::cout << "Report " << onePositionSpecifier.thisRead_ID << "\n";
							std::cout << "\t" << "onePositionSpecifier.mapQ_position" << ": " << onePositionSpecifier.mapQ_position << "\n";
							std::cout << "\t" << "minimumPerPositionMappingQuality" << ": " << minimumPerPositionMappingQuality << "\n";
							std::cout << "\n" << std::flush;
						}
					}
					continue;
				}

				if(perPosition_ignore_alleles.count(onePositionSpecifier.positionInExon) && (perPosition_ignore_alleles.at(onePositionSpecifier.positionInExon).count(onePositionSpecifier.genotype)))
				{
					continue;
				}

				if(ignore_readIDs.count(onePositionSpecifier.thisRead_ID))
				{
					continue;
				}
				
				if((longReadsMode.length() > 0) && (onePositionSpecifier.runningNovelGapEitherDirection >= 2))
				{
						continue;
				}				

				int individualExon = graphLevel_2_exonPosition_individualExon.at(onePositionSpecifier.graphLevel);
				int individualExonPosition = graphLevel_2_exonPosition_individualExonPosition.at(onePositionSpecifier.graphLevel);

				pileUpPerPosition[individualExon][individualExonPosition].push_back(onePositionSpecifier);

				matchesPerReadHistogram_outputStream << locus << "\t" << "base" << onePositionSpecifier.thisRead_WeightedCharactersOK << "\n";

			}
		}

		std::string fileName_pileUp = outputDirectory + "/R1_pileup_"+locus+".txt";
		std::ofstream pileUpStream;
		pileUpStream.open(fileName_pileUp.c_str());
		assert(pileUpStream.is_open());

		for(std::map<int, std::map<int, std::vector<oneExonPosition> > >::iterator exonIt = pileUpPerPosition.begin(); exonIt != pileUpPerPosition.end(); exonIt++)
		{
			int exon = exonIt->first;
			int exonL = exon_lengths.at(exon);
			for(int exonPos = 0; exonPos < exonL; exonPos++)
			{
				if(pileUpPerPosition.at(exon).count(exonPos))
				{
					// int exonPos = exonPosIt->first;
					std::vector<oneExonPosition> pileUp = pileUpPerPosition.at(exon).at(exonPos);

					std::vector<std::string> fieldsPerLine;
					fieldsPerLine.push_back(Utilities::ItoStr(exon));
					fieldsPerLine.push_back(Utilities::ItoStr(exonPos));
					fieldsPerLine.push_back(Utilities::ItoStr(pileUp.size()));

					std::vector<std::string> piledUpGenotypes;
					
					unsigned int thisPos_graphLevel;
					std::map<std::string, std::vector<int>> alleleCounts;
					for(unsigned int pI = 0; pI < pileUp.size(); pI++)
					{
						oneExonPosition piledPosition = pileUp.at(pI);

						std::vector<std::string> qualities_as_strings;
						for(unsigned int qI = 0; qI < piledPosition.qualities.size(); qI++)
						{
							char qC = piledPosition.qualities.at(qI);
							int qC_i = qC;
							qualities_as_strings.push_back(Utilities::ItoStr(qC_i));
						}

						std::string pileUpString = piledPosition.genotype
							+ " (" + Utilities::join(qualities_as_strings, ", ") + ")"
							+ " ["
							// + Utilities::DtoStr(piledPosition.thisRead_WeightedCharactersOK) + " "
							// + Utilities::DtoStr(piledPosition.pairedRead_WeightedCharactersOK) + " | "
							// + Utilities::DtoStr(piledPosition.thisRead_fractionOK) + " "
							// + Utilities::DtoStr(piledPosition.pairedRead_fractionOK) + " | "
							// + Utilities::ItoStr(piledPosition.pairs_strands_OK) + " "
							+ "pairsDistance " + Utilities::DtoStr(piledPosition.pairs_strands_distance) + " | "
							+ "alignmentLength " + Utilities::ItoStr(piledPosition.alignmentColumnsWithAtLeastOneNonGap) + " | "
							+ Utilities::DtoStr(piledPosition.mapQ_position) + " | "
							+ Utilities::DtoStr(piledPosition.mapQ) + " "
							+ Utilities::DtoStr(piledPosition.mapQ_genomic) + " | "
							+ Utilities::DtoStr(piledPosition.thisRead_WeightedCharactersOK) + " "
							+ Utilities::DtoStr(piledPosition.pairedRead_WeightedCharactersOK) + " | "
							+ piledPosition.thisRead_ID + " "
							+ piledPosition.pairedRead_ID
							+ "]";

						utilized_reads.insert(piledPosition.thisRead_ID);
						piledUpGenotypes.push_back(pileUpString);
						
						alleleCounts[piledPosition.genotype].push_back(piledPosition.alignmentColumnsWithAtLeastOneNonGap);
						
						if(pI == 0)
						{
							thisPos_graphLevel = piledPosition.positionInExon;
						}
						else
						{
							assert(thisPos_graphLevel == piledPosition.positionInExon);
						}
					}

					fieldsPerLine.push_back(Utilities::join(piledUpGenotypes, ", ")); 

					std::string components_gtSummary;
					for(auto a : alleleCounts)
					{
						long long combinedLengthts = 0;
						for(auto l : a.second)
						{
							combinedLengthts += l;
						}
						double avgL = (double)combinedLengthts / (double) a.second.size();
						
						std::string allele = a.first;   
						
						
						components_gtSummary += (a.first +"x" + Utilities::ItoStr(a.second.size()) + "[" + Utilities::DtoStr(avgL) + ";" + Utilities::DtoStr(perPosition_allele_minStrandFreqs.at(thisPos_graphLevel).at(allele)) + ";" + Utilities::DtoStr(perPosition_allele_read1Freq.at(thisPos_graphLevel).at(allele)) + "]");
					}
					
					pileUpStream << Utilities::join(fieldsPerLine, "\t") << "\t" << components_gtSummary << "\n";
				}
				else
				{
					std::vector<std::string> fieldsPerLine;
					fieldsPerLine.push_back(Utilities::ItoStr(exon));
					fieldsPerLine.push_back(Utilities::ItoStr(exonPos));
					fieldsPerLine.push_back(Utilities::ItoStr(0));
					pileUpStream << Utilities::join(fieldsPerLine, "\t") << "\n";
				}
			}
		}
		pileUpStream.close();

		std::string fileName_readIDs = outputDirectory + "/R1_readIDs_"+locus+".txt";
		std::ofstream readIDsStream;
		readIDsStream.open(fileName_readIDs.c_str());
		assert(readIDsStream.is_open());
		for(std::set<std::string>::iterator readIDit = utilized_reads.begin(); readIDit != utilized_reads.end(); readIDit++)
		{
			readIDsStream << *readIDit << "\n";
		}
		readIDsStream.close();
		// likelihoods for reads

		std::cout << Utilities::timestamp() << "Compute likelihoods for all exon-overlapping reads (" << exonPositions_fromReads.size() << "), conditional on underlying exons.\n" << std::flush;
		std::vector<std::vector<double> > likelihoods_perCluster_perRead;
		std::vector<std::vector<double> > likelihoods_perCluster_perObservedBase;
		std::vector<std::vector<int> > mismatches_perCluster_perRead;

		likelihoods_perCluster_perRead.resize(HLAtype_clusters.size());
		mismatches_perCluster_perRead.resize(HLAtype_clusters.size());
		likelihoods_perCluster_perObservedBase.resize(HLAtype_clusters.size());


		// std::cout << "\n\n" << HLAtype_2_clusterID.at("A*30:73N") << " " << HLAtype_2_clusterID.at("A*29:02:08") << "\n" << std::flush;
		// assert( 1 == 0 );

		std::set<int> printClusters;
		// printClusters.insert(1789);
		// printClusters.insert(1646);


		// can be removed later, will lead to problems during final normalization
//		assert(exonPositions_fromReads.size() > 0);
		for(unsigned int clusterI = 0; clusterI < HLAtype_clusters.size(); clusterI++)
		{
			std::vector<std::string> typesInCluster(HLAtype_clusters.at(clusterI).begin(), HLAtype_clusters.at(clusterI).end());
			std::string clusterName = Utilities::join(typesInCluster, "|");

			bool verbose = printClusters.count(clusterI);
			// verbose = (clusterI == 0);

			if(verbose)
			{
				std::cout << "CLUSTER " << clusterI << " " << clusterName << "\n";
			}

			std::string& clusterSequence = cluster_2_sequence.at(clusterI);

			// this is really readI
			for(unsigned int positionSpecifierI = 0; positionSpecifierI < exonPositions_fromReads.size(); positionSpecifierI++)
			{
				std::vector<oneExonPosition>& individualPositions = exonPositions_fromReads.at(positionSpecifierI);
				double log_likelihood_read = 0;
				int mismatches = 0;

				std::string readID;
				std::string read1_ID;

				// bool verbose = ((clusterI == 1127) && (readID == "@@B81EP5ABXX:8:2208:11879:23374#GATCAGAT/2") && 0);
				if(verbose)
					std::cout << "Likelihood calculation for read " << readID << " / cluster " << clusterI << "\n";


				for(unsigned int positionI = 0; positionI < individualPositions.size(); positionI++)
				{
					oneExonPosition& onePositionSpecifier = individualPositions.at(positionI);

					assert((onePositionSpecifier.mapQ_position >= 0) && (onePositionSpecifier.mapQ_position <= 1));
					if(onePositionSpecifier.mapQ_position < minimumPerPositionMappingQuality)
					{
						if(verbose)
						{
							std::cout << "\t" << positionI << " mapQ_position " << onePositionSpecifier.mapQ_position << " too low (below " << minimumPerPositionMappingQuality << ").\n" << std::flush;
						}
						continue;
					}

					if(perPosition_ignore_alleles.count(onePositionSpecifier.positionInExon) && (perPosition_ignore_alleles.at(onePositionSpecifier.positionInExon).count(onePositionSpecifier.genotype)))
					{
						continue;
					}

					if(ignore_readIDs.count(onePositionSpecifier.thisRead_ID))
					{
						continue;
					}

					if(clusterI == 0)
					{
						HLATypeInference_thisLocus_bases_used++;
						HLATypeInference_totalBases_used++;
					}

					double log_likelihood_position = 0;

					std::string exonGenotype = clusterSequence.substr(onePositionSpecifier.positionInExon, 1);
					std::string readGenotype = onePositionSpecifier.genotype;
					std::string readQualities = onePositionSpecifier.qualities;

					if(verbose)
					{
						std::cout << "\t" << positionI << " exon pos " << onePositionSpecifier.positionInExon << ": " << exonGenotype << " " << readGenotype << "\n" << std::flush;
					}

					assert(exonGenotype.length() == 1);
					assert(readGenotype.length() >= 1);
					unsigned int l_diff = readGenotype.length() - exonGenotype.length();
					if(exonGenotype == "_")
					{
						// assert(l_diff == 0);
						if(readGenotype == "_")
						{
							assert(onePositionSpecifier.graphLevel != -1);
							// likelihood 1 - intrinsic graph gap

							if(verbose)
							{
								std::cout << "\t\t" << "Intrinsic graph gap" << "\n";
							}

						}
						else
						{
							if(verbose)
							{
								std::cout << "\t\t" << "Insertion " << (1 + l_diff) << "\n";
							}

							assert(readGenotype.find("_") == std::string::npos);
							log_likelihood_position += (log_likelihood_insertion_actualAllele * (1 + l_diff));
						}
					}
					else
					{

						if(readGenotype.length() > 1)
						{
							std::string readGenotype_after1 = readGenotype.substr(1);
							assert(readGenotype_after1.find("_") == std::string::npos);
						}
						// score from first position match
						if(readGenotype.substr(0, 1) == "_")
						{
							log_likelihood_position += log_likelihood_deletion;

							if(verbose)
							{
								std::cout << "\t\t" << "Deletion" << "\n";
							}
						}
						else
						{
							log_likelihood_position += log_likelihood_match_mismatch;

							assert(readQualities.length());
							double pCorrect = Utilities::PhredToPCorrect(readQualities.at(0));
							if(veryConservativeReadLikelihoods)
							{
								if(pCorrect > 0.999)
									pCorrect = 0.999;
							}
							assert((pCorrect >= 0) && (pCorrect <= 1));

							if(pCorrect == 0)
							{
								pCorrect = 0.001;
							}

							if(exonGenotype == readGenotype.substr(0, 1))
							{
								if(verbose)
								{
									std::cout << "\t\t" << "Match " << pCorrect << "\n";
								}

								log_likelihood_position += log(pCorrect);
							}
							else
							{
								double pIncorrect = (1 - pCorrect)*(1.0/3.0);
								assert(pIncorrect <= 0.75);
								assert((pIncorrect > 0) && (pIncorrect < 1));
								log_likelihood_position += log(pIncorrect);

								if(verbose)
								{
									std::cout << "\t\t" << "Mismatch " << pIncorrect << "\n";
								}
							}
						}
						// if read allele is longer
						log_likelihood_position += (log_likelihood_insertion_actualAllele * l_diff);

						if(l_diff > 0)
						{
							if(verbose)
							{
								std::cout << "\t\t" << "Insertion " << l_diff << "\n";
							}
						}
					}

					if(readGenotype != "_")
					{
						if(readGenotype != exonGenotype)
						{
							mismatches++;
						}
					}


					log_likelihood_read += log_likelihood_position;

					if(verbose)
					{
						std::cout << "\t\t" << "Running log likelihood: " << log_likelihood_read << "\n";
					}


					likelihoods_perCluster_perObservedBase.at(clusterI).push_back(log_likelihood_position);
				}

				if(individualPositions.size() > 0)
				{
					if((clusterI == 1160) || (clusterI == 1127))
					{
						// std::cout << locus << " " << readID << " cluster " << clusterI << ": " << log_likelihood << "\n" << std::flush;
					}
				}

				// std::cout << "cluster " << clusterI << ", position sequence " << positionSpecifierI << ": " << log_likelihood << "\n";

				assert(exp(log_likelihood_read) >= 0);
				assert(exp(log_likelihood_read) <= 1);

				likelihoods_perCluster_perRead.at(clusterI).push_back(log_likelihood_read);

				mismatches_perCluster_perRead.at(clusterI).push_back(mismatches);

			}

			assert(likelihoods_perCluster_perRead.at(clusterI).size() == exonPositions_fromReads.size());
			assert(mismatches_perCluster_perRead.at(clusterI).size() == exonPositions_fromReads.size());
		}


		std::cout << Utilities::timestamp() << "Compute likelihoods for all exon cluster pairs (" << HLAtype_clusters.size() << "**2/2)\n" << std::flush;

		std::vector<double> LLs_completeReads;

		std::vector<double> Mismatches_avg;
		std::vector<double> Mismatches_min;

		std::vector<std::pair<unsigned int, unsigned int> > LLs_clusterIs;

		std::cout << "Threads: " << omp_get_num_threads() << "\n";

		size_t HLAtype_cluster_SIZE = HLAtype_clusters.size();

		#pragma omp parallel for schedule(dynamic)
		for(unsigned int clusterI1 = 0; clusterI1 < HLAtype_cluster_SIZE; clusterI1++)
		{
			if(omp_get_thread_num() == 0)
			{
				// std::cout << "\rClusterI1 = " << clusterI1 << std::flush;
			}

			std::vector<double> LLs_completeReads_perThread;

			std::vector<double> Mismatches_avg_perThread;
			std::vector<double> Mismatches_min_perThread;

			std::vector<std::pair<unsigned int, unsigned int> > LLs_clusterIs_perThread;

			for(unsigned int clusterI2 = clusterI1; clusterI2 < HLAtype_clusters.size(); clusterI2++)
			{

				double mismatches_sum_averages = 0;
				double mismatches_sum_min = 0;

				double pair_log_likelihood = 0;

				for(unsigned int readI = 0; readI < exonPositions_fromReads.size(); readI++)
				{
					double LL_thisRead_cluster1 = likelihoods_perCluster_perRead.at(clusterI1).at(readI);
					double LL_thisRead_cluster2 = likelihoods_perCluster_perRead.at(clusterI2).at(readI);

					double LL_average_2_2 = Utilities::logAvg(LL_thisRead_cluster1, LL_thisRead_cluster2);

					if(!((exp(LL_average_2_2) >= 0) && (exp(LL_average_2_2) <= 1)))
					{
						#pragma omp critical
						{
							std::cerr << "Problem with readI " << readI << " conditional on clusters " << clusterI1 << "/" << clusterI2 << "\n";
							std::cerr << "LL_average_2_2" << ": " << LL_average_2_2 << " exp: " << exp(LL_average_2_2) << "\n";
							std::cerr << "LL_thisRead_cluster1" << ": " << LL_thisRead_cluster1 << " exp: " << exp(LL_thisRead_cluster1) << "\n";
							std::cerr << "LL_thisRead_cluster2" << ": " << LL_thisRead_cluster2 << " exp: " << exp(LL_thisRead_cluster2) << "\n";
							std::cerr << std::flush;
						}
					}

					assert((exp(LL_average_2_2) >= 0) && (exp(LL_average_2_2) <= 1));

					int mismatches_cluster1 = mismatches_perCluster_perRead.at(clusterI1).at(readI);
					int mismatches_cluster2 = mismatches_perCluster_perRead.at(clusterI2).at(readI);

					mismatches_sum_averages += ((double)(mismatches_cluster1 + mismatches_cluster2) / 2.0);
					mismatches_sum_min += ((mismatches_cluster1 < mismatches_cluster2) ? mismatches_cluster1 : mismatches_cluster2);

					 pair_log_likelihood += LL_average_2_2;
				}

				LLs_completeReads_perThread.push_back(pair_log_likelihood);

				Mismatches_avg_perThread.push_back(mismatches_sum_averages);
				Mismatches_min_perThread.push_back(mismatches_sum_min);

				LLs_clusterIs_perThread.push_back(make_pair(clusterI1, clusterI2));
			}

			#pragma omp critical
			{
				LLs_completeReads.insert(LLs_completeReads.end(), LLs_completeReads_perThread.begin(), LLs_completeReads_perThread.end());

				Mismatches_avg.insert(Mismatches_avg.end(), Mismatches_avg_perThread.begin(), Mismatches_avg_perThread.end());
				Mismatches_min.insert(Mismatches_min.end(), Mismatches_min_perThread.begin(), Mismatches_min_perThread.end());

				LLs_clusterIs.insert(LLs_clusterIs.end(), LLs_clusterIs_perThread.begin(), LLs_clusterIs_perThread.end());
				// LLs_completeReads_indices.insert(LLs_completeReads_indices.end(), LLs_completeReads_indices_perThread.begin(), LLs_completeReads_indices_perThread.end());
			}
		}

		std::vector<size_t> LLs_completeReads_indices;
		for(size_t i = 0; i < LLs_completeReads.size(); i++)
		{
			LLs_completeReads_indices.push_back(i);
		}

		assert(Mismatches_avg.size() == LLs_completeReads.size());
		assert(Mismatches_min.size() == LLs_completeReads.size());

		// bool verbose = false;

		std::cout << "\n\n" << Utilities::timestamp() << "Sorting...\n\n" << std::flush;

		assert(LLs_completeReads.size() == Mismatches_avg.size());
		assert(LLs_completeReads.size() == LLs_completeReads_indices.size());

		std::sort(LLs_completeReads_indices.begin(), LLs_completeReads_indices.end(), [&](unsigned int a, unsigned int b) {
			if(!(a < LLs_completeReads.size()))
			{
				std::cerr << "a" << ": " << a << "\n";
				std::cerr << "LLs_completeReads.size()" << ": " << LLs_completeReads.size() << "\n" << std::flush;
			}
			assert(a < LLs_completeReads.size());
			assert(b < LLs_completeReads.size());
			assert(a < Mismatches_avg.size());
			assert(b < Mismatches_avg.size());


			if(LLs_completeReads.at(a) == LLs_completeReads.at(b))
			{
				return (Mismatches_avg.at(b) < Mismatches_avg.at(a));
			}
			else
			{
				return (LLs_completeReads.at(a) < LLs_completeReads.at(b));
			}
		});

		std::reverse(LLs_completeReads_indices.begin(), LLs_completeReads_indices.end());

		std::cout << "\n\n" << Utilities::timestamp() << "Sorting done!...\n\n" << std::flush;


		std::pair<double, unsigned int> maxPairI = Utilities::findVectorMax(LLs_completeReads);
		std::cout << Utilities::timestamp() << "Done. (One) maximum pair is " << maxPairI.second << " with LL = " << maxPairI.first << "\n" << std::flush;

		std::vector<double> LLs_normalized;
		double LL_max = maxPairI.first;
		double P_sum = 0;
		for(unsigned int cI = 0; cI < LLs_clusterIs.size(); cI++)
		{
			double LL = LLs_completeReads.at(cI);
			double P = exp(LL - LL_max);
			P_sum += P;
		}
		if(P_sum > 0)
		{
			for(unsigned int cI = 0; cI < LLs_clusterIs.size(); cI++)
			{
				double LL = LLs_completeReads.at(cI);
				double P = exp(LL - LL_max);
				double P_normalized = P / P_sum;
				if(!((P_normalized >= 0) && (P_normalized <= 1)))
				{
					std::cerr << "P_normalized: " << P_normalized << "\n";
					std::cerr << "P: " << P << "\n";
					std::cerr << "LL: " << LL << "\n";
					std::cerr << "LL_max: " << LL_max << "\n";

					std::cerr << std::flush;
				}
				assert(P_normalized >= 0);
				assert(P_normalized <= 1);
				LLs_normalized.push_back(P_normalized);
			}
		}
		else
		{
			for(unsigned int cI = 0; cI < LLs_clusterIs.size(); cI++)
			{
				LLs_normalized.push_back(1.0/(double)LLs_clusterIs.size());
			}

		}

		std::ofstream allPairsStream;
		allPairsStream.open(outputFN_allPairs.c_str());
		assert(allPairsStream.is_open());
		allPairsStream << "ClusterID" << "\t" << "P" << "\t" << "LL" << "\t" << "Mismatches_avg" << "\n";

		std::vector<std::string> LLs_identifiers;
		std::map<int, double> clusterI_overAllPairs;
		for(unsigned int cII = 0; cII < LLs_completeReads_indices.size(); cII++)
		{
			unsigned int cI = LLs_completeReads_indices.at(cII);
			std::pair<unsigned int, unsigned int>& clusters = LLs_clusterIs.at(cI);

			std::vector<std::string> cluster1_members(HLAtype_clusters.at(clusters.first).begin(), HLAtype_clusters.at(clusters.first).end());
			std::vector<std::string> cluster2_members(HLAtype_clusters.at(clusters.second).begin(), HLAtype_clusters.at(clusters.second).end());

			std::string id = Utilities::join(cluster1_members, ";") + "/" + Utilities::join(cluster2_members, ";");

			LLs_identifiers.push_back(id);

			allPairsStream << id << "\t" << LLs_normalized.at(cI) << "\t" << LLs_completeReads.at(cI) << "\t" << Mismatches_avg.at(cI) << "\n";

			if(clusterI_overAllPairs.count(clusters.first) == 0)
			{
				clusterI_overAllPairs[clusters.first] = 0;
			}
			clusterI_overAllPairs[clusters.first] += LLs_normalized.at(cI);

			if(clusters.second != clusters.first)
			{
				if(clusterI_overAllPairs.count(clusters.second) == 0)
				{
					clusterI_overAllPairs[clusters.second] = 0;
				}
				clusterI_overAllPairs[clusters.second] += LLs_normalized.at(cI);
			}
		}

		allPairsStream.close();

		std::pair<double, int> bestGuess_firstAllele = Utilities::findIntMapMax(clusterI_overAllPairs);
		std::string bestGuess_firstAllele_ID = Utilities::join(std::vector<std::string>(HLAtype_clusters.at(bestGuess_firstAllele.second).begin(), HLAtype_clusters.at(bestGuess_firstAllele.second).end()), ";");
		std::string bestGuess_firstAllele_oneType = *(HLAtype_clusters.at(bestGuess_firstAllele.second).begin());

		assert(bestGuess_firstAllele.first >= 0);
		assert(bestGuess_firstAllele.second >= 0);

		std::map<int, double> bestGuess_secondAllele_alternatives;
		std::map<int, double> bestGuess_secondAllele_alternatives_mismatches;

		for(unsigned int cI = 0; cI < LLs_clusterIs.size(); cI++)
		{
			std::pair<unsigned int, unsigned int>& clusters = LLs_clusterIs.at(cI);

			if((int)clusters.first == bestGuess_firstAllele.second)
			{
				assert(bestGuess_secondAllele_alternatives.count(clusters.second) == 0);
				bestGuess_secondAllele_alternatives[clusters.second] = LLs_normalized.at(cI);
				bestGuess_secondAllele_alternatives_mismatches[clusters.second] = Mismatches_min.at(cI);
			}
			else
			{
				if((int)clusters.second == bestGuess_firstAllele.second)
				{
					assert(bestGuess_secondAllele_alternatives.count(clusters.first) == 0);
					bestGuess_secondAllele_alternatives[clusters.first] = LLs_normalized.at(cI);
					bestGuess_secondAllele_alternatives_mismatches[clusters.first] = Mismatches_min.at(cI);
				}
			}
		}

		std::pair<double, int> oneBestGuess_secondAllele = Utilities::findIntMapMax(bestGuess_secondAllele_alternatives);
		assert(oneBestGuess_secondAllele.first >= 0);
		assert(oneBestGuess_secondAllele.first <= 1);

		std::map<int, double> mismatches_allBestGuessPairs;
		for(std::map<int, double>::iterator secondAlleleIt = bestGuess_secondAllele_alternatives.begin(); secondAlleleIt != bestGuess_secondAllele_alternatives.end(); secondAlleleIt++)
		{
			int cluster = secondAlleleIt->first;
			double LL_normalized = secondAlleleIt->second;
			if(LL_normalized == oneBestGuess_secondAllele.first)
			{
				mismatches_allBestGuessPairs[cluster] =  -1 * bestGuess_secondAllele_alternatives_mismatches.at(cluster);
			}
		}

		std::pair<double, int> bestGuess_secondAllele = Utilities::findIntMapMax(mismatches_allBestGuessPairs);
		std::string bestGuess_secondAllele_oneType = *(HLAtype_clusters.at(bestGuess_secondAllele.second).begin());
		std::string bestGuess_secondAllele_ID = Utilities::join(std::vector<std::string>(HLAtype_clusters.at(bestGuess_secondAllele.second).begin(), HLAtype_clusters.at(bestGuess_secondAllele.second).end()), ";");

		forReturn_starting_haplotype_1_vec.push_back(bestGuess_firstAllele_ID);
		forReturn_starting_haplotype_2_vec.push_back(bestGuess_secondAllele_ID);

		double locus_coverage = (double)HLATypeInference_thisLocus_bases_used / (double)thisLocus_totalColumns;;
		std::cout << "Locus " << locus << " " << HLATypeInference_thisLocus_bases_used << " bases across " << thisLocus_totalColumns << " columns utilized." << "\n";
		std::cout << "\tCoverage " << locus_coverage << "\n";

		std::vector<double> positionalCoverages;
		for(unsigned int pI = 0; pI < combined_exon_sequences_graphLevels.size(); pI++)
		{
			int graphLevel = combined_exon_sequences_graphLevels.at(pI);
			int individualExon = graphLevel_2_exonPosition_individualExon.at(graphLevel);
			int exonPosition = graphLevel_2_exonPosition_individualExonPosition.at(graphLevel);
			int coverage = pileUpPerPosition[individualExon][exonPosition].size();
			positionalCoverages.push_back(coverage);
		}
		assert(positionalCoverages.size() > 0);
		std::sort(positionalCoverages.begin(), positionalCoverages.end(), std::less<int>());
		if(positionalCoverages.size() > 1)
		{
			assert(positionalCoverages.at(0) <= positionalCoverages.at(1));
		}


		size_t allColumns_totalAlleles = 0;
		size_t allColumns_incompatibleAlleles = 0;
		std::vector<int> perColumn_totalAlleles;
		std::vector<int> perColumn_incomptibeAlleles;
		double average_perColumn_error_rate;
		double locus_minimumColumnP = -1;
		int n_columns_unaccounted_alleles = 0;

		double proportionkMersCovered_A1;
		double proportionkMersCovered_A2;
		// compute kMer coverages for both called alleles and statistics of deviant alleles
		{
			int clusterID_allele1 = HLAtype_2_clusterID.at(bestGuess_firstAllele_oneType);
			int clusterID_allele2 = HLAtype_2_clusterID.at(bestGuess_secondAllele_oneType);
			assert(clusterID_allele1 == bestGuess_firstAllele.second);
			assert(clusterID_allele2 == bestGuess_secondAllele.second);

			std::string combinedExonSequence_allele1 = combined_exon_sequences.at(bestGuess_firstAllele_oneType);
			std::string combinedExonSequence_allele2 = combined_exon_sequences.at(bestGuess_secondAllele_oneType);

			assert(combinedExonSequence_allele1.length() == combinedExonSequence_allele2.length());
			assert(combined_exon_sequences_graphLevels_individualExon.size() == combinedExonSequence_allele1.size());
			assert(combined_exon_sequences_graphLevels_individualExon.size() == combinedExonSequence_allele2.size());
			assert(combined_exon_sequences_graphLevels_individualExon.size() == combined_exon_sequences_graphLevels_individualExonPosition.size());

			std::vector<std::string> combinedExonSequence_allele1_byExon;
			std::vector<std::string> combinedExonSequence_allele2_byExon;

			for(unsigned int pI = 0; pI < combinedExonSequence_allele1.size(); pI++)
			{
				int thisPosition_individualExon = combined_exon_sequences_graphLevels_individualExon.at(pI);
				int thisPosition_individualExonPosition = combined_exon_sequences_graphLevels_individualExonPosition.at(pI);

				bool needNewExon = ((pI == 0) || (combined_exon_sequences_graphLevels_individualExon.at(pI) != combined_exon_sequences_graphLevels_individualExon.at(pI-1)));
				if(needNewExon)
				{
					combinedExonSequence_allele1_byExon.push_back("");
					combinedExonSequence_allele2_byExon.push_back("");
				}

				std::string underlyingAllele1 = combinedExonSequence_allele1.substr(pI, 1);
				std::string underlyingAllele2 = combinedExonSequence_allele2.substr(pI, 1);

				combinedExonSequence_allele1_byExon.back().push_back(underlyingAllele1.at(0));
				combinedExonSequence_allele2_byExon.back().push_back(underlyingAllele2.at(0));

				int totalAlleles = 0;
				int incompatibeAlleles = 0;
				for(auto piledUpAllele : pileUpPerPosition.at(thisPosition_individualExon).at(thisPosition_individualExonPosition))
				{
					totalAlleles++;
					if((piledUpAllele.genotype != underlyingAllele1) && (piledUpAllele.genotype != underlyingAllele2))
					{
						incompatibeAlleles++;
					}
				}

				allColumns_totalAlleles += totalAlleles;
				allColumns_incompatibleAlleles += incompatibeAlleles;

				perColumn_totalAlleles.push_back(totalAlleles);
				perColumn_incomptibeAlleles.push_back(incompatibeAlleles);

				if(perPosition_allele_counts_postFiltering.count(pI))
				{
					int totalCoverage_position = 0;
					for(auto alleleCount : perPosition_allele_counts_postFiltering.at(pI))
					{
						totalCoverage_position += alleleCount.second;
					}
					if(totalCoverage_position >= threshold_reportColumn_forPresenceOfUnaccountedAlleles_minCoverage)
					{
						for(auto alleleCount : perPosition_allele_counts_postFiltering.at(pI))
						{
							if(alleleCount.first == underlyingAllele1)
								continue;
							if(alleleCount.first == underlyingAllele2)
								continue;

							double aF = (double)alleleCount.second / (double)totalCoverage_position;
							if(aF >= threshold_reportColumn_forPresenceOfUnaccountedAlleles_minAlleleFraction)
							{
								n_columns_unaccounted_alleles++;
							}
						}
					}
				}
			}

			auto calculcatekMerPresence = [&](std::vector<std::string> exons) -> double {
				int kMers_total = 0;
				int kMers_present = 0;
				for(auto exonSeq : exons)
				{
					std::string exonSeqNoGaps = Utilities::removeGaps(exonSeq);
					std::vector<std::string> kMers = Utilities::partitionStringIntokMers(exonSeqNoGaps, k_for_kMer_index);
					for(auto kMer : kMers)
					{
						kMers_total++;
						if(kMer.find("*") == std::string::npos)
						{
							std::string kMerKey = kMer_canonical_representation(kMer);
							if(kMer_counts.count(kMerKey) && (kMer_counts.at(kMerKey) > 0))
							{
								kMers_present++;
							}
						}
					}
				}
				if(kMers_total == 0)
				{
					return -1;
					std::cerr << "kMers_total: " << kMers_total << "\n" << std::flush;
					for(unsigned int exonI = 0; exonI < exons.size(); exonI++)
					{
						std::cerr << "\texon " << exonI << ": " << exons.at(exonI) << "\n" << std::flush;
					}

				}
				assert(kMers_total > 0);
				return (double) kMers_present / (double) kMers_total;
			};

			proportionkMersCovered_A1 = calculcatekMerPresence(combinedExonSequence_allele1_byExon);
			proportionkMersCovered_A2 = calculcatekMerPresence(combinedExonSequence_allele2_byExon);


			std::ofstream columnErrorRateStream;
			columnErrorRateStream.open(outputFN_columnError.c_str());
			assert(columnErrorRateStream.is_open());
			columnErrorRateStream << Utilities::join({"Column", "Coverage", "ExpectedIncompatible", "ObservedIncompatible", "p"}, "\t") << "\n";

			if(allColumns_totalAlleles > 0)
			{
				average_perColumn_error_rate = (double)allColumns_incompatibleAlleles / (double)allColumns_totalAlleles;
			}
			else
			{
				average_perColumn_error_rate = 0;
			}

			for(unsigned int columnI = 0; columnI < combinedExonSequence_allele1.size(); columnI++)
			{
				std::vector<std::string> outputFields;
				outputFields.push_back(Utilities::ItoStr(columnI));

				int columnCoverage = perColumn_totalAlleles.at(columnI);
				outputFields.push_back(Utilities::ItoStr(columnCoverage));

				double expectedError = average_perColumn_error_rate * columnCoverage;
				outputFields.push_back(Utilities::DtoStr(expectedError));

				int observedError = perColumn_incomptibeAlleles.at(columnI);
				outputFields.push_back(Utilities::ItoStr(observedError));

				double p = 1;
				if(observedError > expectedError)
				{
					std::vector<double> observed;
						observed.push_back(columnCoverage - observedError);
						observed.push_back(observedError);

					std::vector<double> expected;
						expected.push_back(columnCoverage - expectedError);
						expected.push_back(expectedError);

					p = simpleChiSq(observed, expected);

				}

				outputFields.push_back(Utilities::DtoStr(p));
				columnErrorRateStream << Utilities::join(outputFields, "\t") << "\n";

				if((locus_minimumColumnP < 0) || (p < locus_minimumColumnP))
				{
					locus_minimumColumnP = p;
				}
			}
		}

		int index_for_decile = (int)((double)positionalCoverages.size() / 10.0);
		double firstDecileCoverage = positionalCoverages.at(index_for_decile);
		double minimumCoverage = positionalCoverages.at(0);
		bestGuess_outputStream << locus << "\t" << 1 << "\t" << bestGuess_firstAllele_ID << "\t" << bestGuess_firstAllele.first << "\t" << bestGuess_secondAllele.first << "\t" << locus_coverage << "\t" << firstDecileCoverage << "\t" << minimumCoverage << "\t" << proportionkMersCovered_A1 << "\t" << average_perColumn_error_rate << "\t" << n_columns_unaccounted_alleles << "\n";
		bestGuess_outputStream << locus << "\t" << 2 << "\t" << bestGuess_secondAllele_ID << "\t" << oneBestGuess_secondAllele.first << "\t" << bestGuess_secondAllele.first << "\t" << locus_coverage << "\t" << firstDecileCoverage << "\t" << minimumCoverage << "\t" << proportionkMersCovered_A2 << "\t" << average_perColumn_error_rate << "\t" << n_columns_unaccounted_alleles << "\n" << std::flush;

		if(can_translateToG_locus(locus))
		{
			bool a1_perfectly;
			bool a2_perfectly;
			std::string bestGuess_firstAllele_ID_G = translate_allele_list_to_G_allele({HLAtype_clusters.at(bestGuess_firstAllele.second).begin(), HLAtype_clusters.at(bestGuess_firstAllele.second).end()}, a1_perfectly);
			std::string bestGuess_secondAllele_ID_G = translate_allele_list_to_G_allele({HLAtype_clusters.at(bestGuess_secondAllele.second).begin(), HLAtype_clusters.at(bestGuess_secondAllele.second).end()}, a2_perfectly);

			bestGuess_G_outputStream << locus << "\t" << 1 << "\t" << bestGuess_firstAllele_ID_G << "\t" << bestGuess_firstAllele.first << "\t" << bestGuess_secondAllele.first << "\t" << locus_coverage << "\t" << firstDecileCoverage << "\t" << minimumCoverage << "\t" << proportionkMersCovered_A1 << "\t" << average_perColumn_error_rate << "\t" << n_columns_unaccounted_alleles <<  "\t" << a1_perfectly << "\n";
			bestGuess_G_outputStream << locus << "\t" << 2 << "\t" << bestGuess_secondAllele_ID_G << "\t" << oneBestGuess_secondAllele.first << "\t" << bestGuess_secondAllele.first << "\t" << locus_coverage << "\t" << firstDecileCoverage << "\t" << minimumCoverage << "\t" << proportionkMersCovered_A2 << "\t" << average_perColumn_error_rate << "\t" << n_columns_unaccounted_alleles << "\t" << a2_perfectly << "\n" << std::flush;
		}


		unsigned int maxPairPrint = (LLs_completeReads_indices.size() > 10) ? 10 : LLs_completeReads_indices.size();
		for(unsigned int LLi = 0; LLi < maxPairPrint; LLi++)
		{
			unsigned int pairIndex = LLs_completeReads_indices.at(LLi);
			std::pair<unsigned int, unsigned int>& clusters = LLs_clusterIs.at(pairIndex);
			// std::cout << "#" << (LLi+1) << ": " << clusters.first << " / " << clusters.second << ": " << LLs.at(pairIndex) << " absolute and " << exp(LLs.at(pairIndex) - maxPairI.first) << " relative." << "\n" << std::flush;

			std::vector<std::string> cluster1_members(HLAtype_clusters.at(clusters.first).begin(), HLAtype_clusters.at(clusters.first).end());
			std::vector<std::string> cluster2_members(HLAtype_clusters.at(clusters.second).begin(), HLAtype_clusters.at(clusters.second).end());

			// std::cout << "\tcluster " << clusters.first << ": " << Utilities::join(cluster1_members, ", ") << "\n";
			// std::cout << "\tcluster " << clusters.second << ": " << Utilities::join(cluster2_members, ", ") << "\n" << std::flush;

			// std::cout << "\tMismatches " << Mismatches_avg.at(pairIndex) << " avg / " << Mismatches_min.at(pairIndex) << " min\n" << std::flush;

			// bool equal = (LLs.at(pairIndex) == LLs.at(LLs_indices.at(0)));
			// std::cout << "\tEqual: " << equal << "\n" << std::flush;


		}

		assert(LLs_completeReads.at(LLs_completeReads_indices.at(0)) == maxPairI.first);
		if(LLs_completeReads_indices.size() > 1)
		{
			assert(LLs_completeReads.at(LLs_completeReads_indices.at(0)) >= LLs_completeReads.at(LLs_completeReads_indices.at(1)));
		}


	}

	bestGuess_outputStream.close();

//	forReturn_starting_haplotype_1 = Utilities::join(forReturn_starting_haplotype_1_vec, ",");
//	forReturn_starting_haplotype_2 = Utilities::join(forReturn_starting_haplotype_2_vec, ",");

	std::string outputFN_parameters = outputDirectory + "/R1_parameters.txt";
	std::ofstream outputFN_parameters_outputStream;
	outputFN_parameters_outputStream.open(outputFN_parameters.c_str());
	assert(outputFN_parameters_outputStream.is_open());

	outputFN_parameters_outputStream << "Loci" << " = " <<  Utilities::join(loci_for_HLAtyping, ",") << "\n";
//	outputFN_parameters_outputStream << "restrictToFullHaplotypes" << " = " << restrictToFullHaplotypes << "\n";
	outputFN_parameters_outputStream << "veryConservativeReadLikelihoods" << " = " << veryConservativeReadLikelihoods << "\n";
	outputFN_parameters_outputStream.close();

	double coverage = (double)HLATypeInference_totalBases_used / (double) HLAtypeInference_totalColumns;
	std::cout << "Sample " << " " << HLATypeInference_totalBases_used << " bases utilized (total) across" << HLAtypeInference_totalColumns <<  "columns\n" << std::flush;
	std::cout << "\tCoverage " << coverage << "\n" << std::flush;
}

void HLATyper::fill_loci_2_exons()
{
	std::vector<std::string> exons_A = {"exon_2", "exon_3"};
	loci_2_exons["A"] = exons_A;

	std::vector<std::string> exons_B = {"exon_2", "exon_3"};
	loci_2_exons["B"] = exons_B;

	std::vector<std::string> exons_C = {"exon_2", "exon_3"};
	loci_2_exons["C"] = exons_C;

	std::vector<std::string> exons_DQA1 = {"exon_2"};
	loci_2_exons["DQA1"] = exons_DQA1;

	std::vector<std::string> exons_DQB1 = {"exon_2"};
	loci_2_exons["DQB1"] = exons_DQB1;

	std::vector<std::string> exons_DRB1 = {"exon_2"};
	loci_2_exons["DRB1"] = exons_DRB1;

	loci_2_exons["DPA1"] = {"exon_2"};
	loci_2_exons["DPB1"] = {"exon_2"};
	loci_2_exons["DRA"] = {"exon_2"};
	loci_2_exons["DRB3"] = {"exon_2"};
	loci_2_exons["DRB4"] = {"exon_2"};

	loci_2_exons["E"] = {"exon_2", "exon_3"};
	loci_2_exons["F"] = {"exon_2", "exon_3"};
	loci_2_exons["G"] = {"exon_2", "exon_3"};
	loci_2_exons["H"] = {"exon_2", "exon_3"};
	loci_2_exons["J"] = {"exon_2", "exon_3"};
	loci_2_exons["K"] = {"exon_2", "exon_3"};
	loci_2_exons["L"] = {"exon_2", "exon_3"};
	loci_2_exons["V"] = {"exon_2", "exon_3"};
}

std::set<std::string> HLATyper::getCompletelyDefinedHLAAlleles(std::string locus)
{
	std::vector<std::string> files_in_order;
	std::vector<std::string> files_in_order_type;
	std::vector<int> files_in_order_number;

	std::set<std::string> loci;
	std::map<std::string, std::vector<std::string> > files_per_locus;
	std::map<std::string, std::vector<std::string> > files_per_locus_type;
	std::map<std::string, std::vector<std::string> > files_per_locus_number;

	std::ifstream segmentsStream;
	std::string segmentsFileName = graphDir + "/PRG/segments.txt";
	segmentsStream.open(segmentsFileName.c_str());
	assert(segmentsStream.is_open());
	std::string line;
	while(segmentsStream.good())
	{
		std::getline(segmentsStream, line);
		Utilities::eraseNL(line);
		if(line.length() > 0)
		{
			std::vector<std::string> split_by_underscore = Utilities::split(line, "_");
			if(split_by_underscore.at(1) != "gene")
			{
				continue;
			}

			std::string file_locus = split_by_underscore.at(2);
			loci.insert(file_locus);

			std::string filePath = graphDir + "/PRG/" + line;

			files_per_locus[file_locus].push_back(filePath);
			std::string type;
			if(split_by_underscore.at(4).find("intron") != std::string::npos)
			{
				type = "intron";
			}
			else
			{
				type = "exon";
				assert(split_by_underscore.at(4).find("exon") != std::string::npos);
			}
			files_per_locus_type[file_locus].push_back(type);

//			if(split_by_underscore.size() >= 4)
//			{
//				std::vector<std::string> split_by_colon = Utilities::split(split_by_underscore.at(3), ".");
//				assert(split_by_colon.size() == 2);
//				files_per_locus_number[file_locus].push_back(split_by_colon.at(0));
//			}
//			else
//			{
//				files_per_locus_number[file_locus].push_back("");
//			}

		}
	}
	assert(loci.size() > 0);
	assert(loci.count(locus));


	std::string arbitraryIntronFile;
	for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
	{
		if(files_per_locus_type.at(locus).at(fI) == "intron")
		{
			arbitraryIntronFile = files_per_locus.at(locus).at(fI);
			break;
		}
	}

	if(arbitraryIntronFile.length() == 0)
	{
		for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
		{
			if(files_per_locus_type.at(locus).at(fI) == "exon")
			{
				arbitraryIntronFile = files_per_locus.at(locus).at(fI);
				break;
			}
		}
	}

	if(!(arbitraryIntronFile.length()))
	{
		std::cerr << "Cannot find an arbitrary intron/exon file for locus " << locus << "\n" << std::flush;
	}
	assert(arbitraryIntronFile.length());

	std::set<std::string> availableTypes;

	std::ifstream fileInputStream;
	fileInputStream.open(arbitraryIntronFile.c_str());
	if(!fileInputStream.is_open())
	{
		std::cerr << "Cannot open file " << arbitraryIntronFile << "\n" << std::flush;
	}
	assert(fileInputStream.is_open());
	unsigned int lI = 0;
	while(fileInputStream.good())
	{
		std::string line;
		std::getline(fileInputStream, line);
		Utilities::eraseNL(line);
		std::vector<std::string> line_fields = Utilities::split(line, " ");
		if(lI == 0)
		{
			if(!(line_fields.size() > 0))
			{
				std::cerr << "First line of file weird: " << arbitraryIntronFile << "\n" << std::flush;
			}
			assert(line_fields.size() > 0);
			assert(line_fields.at(0) == "IndividualID");
		}
		else
		{
			if(line.length())
			{
				std::string type = line_fields.at(0);
				availableTypes.insert(type);
			}
		}
		lI++;
	}
	fileInputStream.close();

	// find type haplotypes

	std::map<std::string, std::set<std::string> > removeTypes;

	std::map<std::string, std::map<std::string, std::vector<std::string>> > loci_types_haplotypes;
	std::map<std::string, std::vector<std::string> > loci_graphLevelIDs;

	for(unsigned int fI = 0; fI < files_per_locus.at(locus).size(); fI++)
	{
		std::ifstream fileInputStream;
		fileInputStream.open(files_per_locus.at(locus).at(fI).c_str());
		assert(fileInputStream.is_open());
		std::vector<std::string> file_lines;
		while(fileInputStream.good())
		{
			std::string line;
			std::getline(fileInputStream, line);
			Utilities::eraseNL(line);
			file_lines.push_back(line);
		}
		fileInputStream.close();

		std::string firstLine = file_lines.at(0);
		std::vector<std::string> firstLine_fields = Utilities::split(firstLine, " ");
		assert(firstLine_fields.at(0) == "IndividualID");

		std::vector<std::string> graph_level_names(firstLine_fields.begin() + 1, firstLine_fields.end());
		loci_graphLevelIDs[locus].insert(loci_graphLevelIDs[locus].end(), graph_level_names.begin(), graph_level_names.end());

		for(unsigned int lI = 1; lI < file_lines.size(); lI++)
		{
			if(file_lines.at(lI).length())
			{
				std::vector<std::string> line_fields = Utilities::split(file_lines.at(lI), " ");
				assert(line_fields.size() == firstLine_fields.size());
				std::string HLA_type = line_fields.at(0);
				std::vector<std::string> line_alleles(line_fields.begin()+1, line_fields.end());

				assert((files_per_locus_type.at(locus).at(fI) == "intron") || (files_per_locus_type.at(locus).at(fI) == "exon"));
				{
					if(availableTypes.count(HLA_type))
					{
						if(loci_types_haplotypes[locus].count(HLA_type) == 0)
						{
							loci_types_haplotypes[locus][HLA_type].resize(0);
						}

						loci_types_haplotypes.at(locus).at(HLA_type).insert(loci_types_haplotypes.at(locus).at(HLA_type).end(), line_alleles.begin(), line_alleles.end());
					}
				}
			}
		}
	}

	// check that all haplotypes make sense
	for(std::set<std::string>::iterator HLA_typeIT = availableTypes.begin(); HLA_typeIT != availableTypes.end(); HLA_typeIT++)
	{
		std::string HLA_type = *HLA_typeIT;

		for(unsigned int gI = 0; gI < loci_types_haplotypes.at(locus).at(HLA_type).size(); gI++)
		{
			std::string& S = loci_types_haplotypes.at(locus).at(HLA_type).at(gI);
			if(S.find("*") != std::string::npos)
			{
				if(removeTypes[locus].count(HLA_type) == 0)
				{
					removeTypes[locus].insert(HLA_type);
					std::cerr << "Locus " << locus << " allele " << HLA_type << " out because of stars.\n\n" << std::flush;
					std::cerr << gI << "/" << loci_types_haplotypes.at(locus).at(HLA_type).size() << ": " <<  S << "\n\n" << std::flush;
				}
			}
		}

		if(!(loci_types_haplotypes.at(locus).at(HLA_type).size() == loci_graphLevelIDs.at(locus).size()))
		{
			std::cerr << "Length problem for locus " << locus << " allele " << HLA_type << ": " << loci_types_haplotypes.at(locus).at(HLA_type).size() << " vs " << loci_graphLevelIDs.at(locus).size() << "\n" << std::flush;
			removeTypes[locus].insert(HLA_type);
			continue;
		}
		assert(loci_types_haplotypes.at(locus).at(HLA_type).size() == loci_graphLevelIDs.at(locus).size());
	}

	if(removeTypes.count(locus) && (removeTypes.at(locus).size() > 0))
	{
		std::cerr << "\t\t" << locus << " before removal: " << availableTypes.size() << " types.\n";

		std::set<std::string> availableTypes_new;
		for(std::set<std::string>::iterator HLA_typeIT = availableTypes.begin(); HLA_typeIT != availableTypes.end(); HLA_typeIT++)
		{
			std::string HLA_type = *HLA_typeIT;
			if(removeTypes.at(locus).count(HLA_type) == 0)
			{
				availableTypes_new.insert(HLA_type);
			}
		}
		availableTypes = availableTypes_new;

		std::cerr << "\t\t\t" << locus << " after removal: " << availableTypes.size() << " types.\n" << std::flush;

		assert(availableTypes.size() > 0);
	}

	return availableTypes;
}


double HLATyper::alignmentFractionOK (const mapper::reads::verboseSeedChain& r)
{
	int positions_OK = 0;
	int positions_checked = 0;

	for(unsigned int i = 0; i < r.graph_aligned.size(); i++)
	{
		if((r.graph_aligned.at(i) == '_') && (r.sequence_aligned.at(i) == '_'))
		{
			continue;
		}

		positions_checked++;
		if(r.graph_aligned.at(i) == r.sequence_aligned.at(i))
		{
			positions_OK++;
		}
	}
	return double(positions_OK)/double(positions_checked);
}

std::pair<double, double> HLATyper::meanMedian (std::vector<double> L)
{
	std::sort(L.begin(), L.end());
	double S = 0;
	for(unsigned int i = 0; i < L.size(); i++)
	{
		S += L.at(i);
		// std::cout << L.at(i) << " ";
	}
	double mean = 0;
	double median = 0;
	if(L.size() > 0)
	{
		mean = S / (double)L.size();
		unsigned int medium_index = L.size() / 2;
		median = L.at(medium_index);
	}
	return make_pair(mean, median);
}


std::string HLATyper::printPerc (double v1, double v2)
{
	double perc = (v1/v2) * 100;
	return Utilities::DtoStr(perc);
}

std::string HLATyper::find_file_for_exon(std::string locus, std::string exon)
{
	std::vector<std::string> exon_parts = Utilities::split(exon, "_");
	assert(exon_parts.size() == 2);
	assert(exon_parts.at(0) == "exon");
	int exonN = Utilities::StrtoI(exon_parts.at(1));
	assert(exonN > 0);
	std::string forReturn;
	for(unsigned int fI = 0; fI < files_in_graphDir.size(); fI++)
	{
		assert(files_in_graphDir.at(fI).find("\\") == std::string::npos);
		std::vector<std::string> split_by_slash = Utilities::split(files_in_graphDir.at(fI), "/");
		if(split_by_slash.size() > 0)
		{
			std::string filename = split_by_slash.at(split_by_slash.size() - 1); // 146_gene_HLA-C_15_exon_8.txt
			std::vector<std::string> split_by_underscore = Utilities::split(filename, "_");
			// if(!(split_by_underscore.size() >= 3))
			// {
				// std::cerr << "Can't split according to underscores: " << files_in_graphDir.at(fI) << "\n" << std::flush;
			// }
			// assert(split_by_underscore.size() >= 3);

			
			if(split_by_underscore.size() >= 6)
			{
				if(split_by_underscore.at(1) == "gene")
				{
					std::string hlaLocus = "HLA-"+locus;
					if((split_by_underscore.at(2) == hlaLocus) || (split_by_underscore.at(2) == locus))
					{
						if(split_by_underscore.at(4) == "exon")
						{
							if(split_by_underscore.at(5) == (Utilities::ItoStr(exonN)+".txt"))
							{
								forReturn = files_in_graphDir.at(fI);
							}
						}
					}
				}
				
				/*
				if(split_by_underscore.at(split_by_underscore.size()-4).substr(split_by_underscore.at(split_by_underscore.size()-4).length() - (locus.length()+1)) == ("/"+locus))
				{
					if((split_by_underscore.at(split_by_underscore.size()-2)+"_"+split_by_underscore.at(split_by_underscore.size()-1)) == (exon + ".txt"))
					{
						forReturn = files_in_graphDir.at(fI);
					}
				}
				*/
			}
		}
	}
	if(! forReturn.length())
	{
		std::cerr << "find_file_for_exon -- problem -- " << locus << " -- " << exon << "\n";
		//std::cerr << Utilities::join(files_in_graphDir, " ") << "\n" << std::flush; // todo 
		//assert(forReturn.length()); // todo
	}
	return forReturn;
}


void HLATyper::oneReadAlignment_2_exonPositions_paired(const mapper::reads::verboseSeedChain& alignment, const mapper::reads::oneRead& read, std::vector<oneExonPosition>& ret_exonPositions, const mapper::reads::verboseSeedChain& paired_alignment, const mapper::reads::oneRead& paired_read, int read_1_or_2, int combined_exon_sequences_graphLevels_min, int combined_exon_sequences_graphLevels_max, const std::map<int, unsigned int>& graphLevel_2_exonPosition)
{
	int alignment_firstLevel = alignment.alignment_firstLevel();
	int alignment_lastLevel = alignment.alignment_lastLevel();

	assert(alignment_firstLevel <= alignment_lastLevel);
	double thisRead_fractionOK = alignmentFractionOK(alignment);
	double pairedRead_fractionOK = alignmentFractionOK(paired_alignment);

	double thisRead_WeightedCharactersOK = alignmentWeightedOKFraction(read, alignment);
	double pairedRead_WeightedCharactersOK = alignmentWeightedOKFraction(paired_read, paired_alignment);

	double pairs_strands_OK = mapper::aligner::alignerBase::alignedReadPair_strandsValid(alignment, paired_alignment);
	double pairs_strands_distance = mapper::aligner::alignerBase::alignedReadPair_pairsDistanceInGraphLevels(alignment, paired_alignment);

	// std::cout << "This alignment " << alignment_firstLevel << " - " << alignment_lastLevel << "\n";
	// std::cout << "\tvs combined exon " << combined_exon_sequences_graphLevels.front() << " - " << combined_exon_sequences_graphLevels.back() << "\n\n" << std::flush;


	if(globalVerbose && (read.name.find("HLA-C") != std::string::npos))
	{
		std::cout << "Read" << ": " << read.name << "\n";
		std::cout << "\t" << "alignment_firstLevel" << ": " << alignment_firstLevel << "\n";
		std::cout << "\t" << "alignment_lastLevel" << ": " << alignment_lastLevel << "\n";
		std::cout << "\t" << "combined_exon_sequences_graphLevels_min" << ": " << combined_exon_sequences_graphLevels_min << "\n";
		std::cout << "\t" << "combined_exon_sequences_graphLevels_max" << ": " << combined_exon_sequences_graphLevels_max << "\n";
		std::cout << "\n" << std::flush;
	}

	assert((combined_exon_sequences_graphLevels_min <= combined_exon_sequences_graphLevels_max));

	// if( ((alignment_firstLevel >= combined_exon_sequences_graphLevels_front) && (alignment_firstLevel <= combined_exon_sequences_graphLevels_back)) ||
	//	((alignment_lastLevel >= combined_exon_sequences_graphLevels_front) && (alignment_lastLevel <= combined_exon_sequences_graphLevels_back)) )
	if(Utilities::intervalsOverlap(alignment_firstLevel, alignment_lastLevel, combined_exon_sequences_graphLevels_min, combined_exon_sequences_graphLevels_max))
	{
		std::vector<oneExonPosition> readAlignment_exonPositions;

		int alignmentColumns_oneNonGap = 0;
		for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
		{
			std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
			std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);		
			
			if((sequenceCharacter != "_") || (sequenceCharacter != "_"))
			{
				alignmentColumns_oneNonGap++;
			}
		}
		
		
		std::vector<int> runningNovelGaps;
		runningNovelGaps.resize(alignment.sequence_aligned.length(), 0);
		{
			int runningNovelGap = 0;
			for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
			{
				std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
				std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
				
				if((graphCharacter != "_") && (sequenceCharacter != "_"))
				{
					runningNovelGap = 0;
				}
				else
				{
					if(!((graphCharacter == "_") && (sequenceCharacter == "_")))
					{
						runningNovelGap++;
					}					
				}
			
				if(runningNovelGap > runningNovelGaps.at(cI))
					runningNovelGaps.at(cI) = runningNovelGap;
			}
			
			runningNovelGap = 0;
			for(int cI = alignment.sequence_aligned.length()-1; cI >= 0; cI--)
			{
				std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
				std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
				
				if((graphCharacter != "_") && (sequenceCharacter != "_"))
				{
					runningNovelGap = 0;
				}
				else
				{
					if(!((graphCharacter == "_") && (sequenceCharacter == "_")))
					{
						runningNovelGap++;
					}					
				}
			
				if(runningNovelGap > runningNovelGaps.at(cI))
					runningNovelGaps.at(cI) = runningNovelGap;
			}		
		}
		
		int indexIntoOriginalReadData = -1;
		for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
		{
			std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
			std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
			int graphLevel = alignment.graph_aligned_levels.at(cI);
			unsigned char alignmentQualityCharacter_thisPosition = alignment.mapQ_perPosition.at(cI);
			double alignmentQuality_thisPosition = Utilities::PhredToPCorrect(alignmentQualityCharacter_thisPosition);

			if(graphLevel == -1)
			{
				// insertion relative to the graph - we need to extend last character

				assert(graphCharacter == "_");
				assert(sequenceCharacter != "_");

				indexIntoOriginalReadData++;
				int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
				if(alignment.reverse)
				{
					indexIntoOriginalReadData_correctlyAligned = read.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
				}
				assert(indexIntoOriginalReadData_correctlyAligned >= 0);
				assert(indexIntoOriginalReadData_correctlyAligned <(int) read.sequence.length());

				std::string underlyingReadCharacter = read.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
				if(alignment.reverse)
				{
					underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
				}
				assert(underlyingReadCharacter == sequenceCharacter);
				char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

				if(readAlignment_exonPositions.size() > 0)
				{
					readAlignment_exonPositions.back().genotype.append(sequenceCharacter);
					readAlignment_exonPositions.back().alignment_edgelabels.append(graphCharacter);
					readAlignment_exonPositions.back().qualities.push_back(qualityCharacter);
					if(!(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length()))
					{
						assert(readAlignment_exonPositions.back().genotype.length() == (readAlignment_exonPositions.back().qualities.length()+1));
						assert(readAlignment_exonPositions.back().genotype.at(0) == '_');
						readAlignment_exonPositions.back().genotype = readAlignment_exonPositions.back().genotype.substr(1);
						readAlignment_exonPositions.back().alignment_edgelabels = readAlignment_exonPositions.back().alignment_edgelabels.substr(1);
						assert(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length());

						// std::cerr << "readAlignment_exonPositions.back().genotype.length()" << ": " << readAlignment_exonPositions.back().genotype.length() << "\n";
						// std::cerr << "readAlignment_exonPositions.back().qualities.length()" << ": " << readAlignment_exonPositions.back().qualities.length()<< "\n";

						// std::cerr << std::flush;
					}
					assert(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length());
					assert(readAlignment_exonPositions.back().alignment_edgelabels.length() == readAlignment_exonPositions.back().genotype.length());
				}
			}
			else
			{
				if(sequenceCharacter != "_")
				{
					indexIntoOriginalReadData++;
					int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
					if(alignment.reverse)
					{
						indexIntoOriginalReadData_correctlyAligned = read.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
					}
					assert(indexIntoOriginalReadData_correctlyAligned >= 0);
					assert(indexIntoOriginalReadData_correctlyAligned < (int)read.sequence.length());

					std::string underlyingReadCharacter = read.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
					if(alignment.reverse)
					{
						underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
					}
					assert(underlyingReadCharacter == sequenceCharacter);

					if(graphCharacter == "_")
					{

						assert(graphLevel != -1);

						char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = sequenceCharacter;
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities.push_back(qualityCharacter);

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = paired_read.name;
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = (read_1_or_2 == 1) ? read.name : paired_read.name;
						thisPosition.alignmentColumnsWithAtLeastOneNonGap = alignmentColumns_oneNonGap;
						thisPosition.runningNovelGapEitherDirection = runningNovelGaps.at(cI);
						thisPosition.reverse = alignment.reverse;
						thisPosition.fromFirstRead = alignment.fromFirstRead;
						
						readAlignment_exonPositions.push_back(thisPosition);

					}
					else
					{
						// two well-defined characters
						char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = sequenceCharacter;
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities.push_back(qualityCharacter);

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = paired_read.name;
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;
						thisPosition.alignmentColumnsWithAtLeastOneNonGap = alignmentColumns_oneNonGap;
						thisPosition.runningNovelGapEitherDirection = runningNovelGaps.at(cI);
						thisPosition.reverse = alignment.reverse;
						thisPosition.fromFirstRead = alignment.fromFirstRead;
						
						thisPosition.read1_ID = (read_1_or_2 == 1) ? read.name : paired_read.name;

						readAlignment_exonPositions.push_back(thisPosition);
					}

				}
				else
				{
					assert(sequenceCharacter == "_");
					if(graphCharacter == "_")
					{
						assert(graphLevel != -1);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = "_";
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities = "";

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = paired_read.name;
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = (read_1_or_2 == 1) ? read.name : paired_read.name;
						thisPosition.alignmentColumnsWithAtLeastOneNonGap = alignmentColumns_oneNonGap;
						thisPosition.runningNovelGapEitherDirection = runningNovelGaps.at(cI);
						thisPosition.reverse = alignment.reverse;
						thisPosition.fromFirstRead = alignment.fromFirstRead;
						
						readAlignment_exonPositions.push_back(thisPosition);
					}
					else
					{
						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = "_";
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities = "";

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = paired_read.name;
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = (read_1_or_2 == 1) ? read.name : paired_read.name;
						thisPosition.alignmentColumnsWithAtLeastOneNonGap = alignmentColumns_oneNonGap;
						thisPosition.runningNovelGapEitherDirection = runningNovelGaps.at(cI);
						thisPosition.reverse = alignment.reverse;
						thisPosition.fromFirstRead = alignment.fromFirstRead;
						
						readAlignment_exonPositions.push_back(thisPosition);
					}
				}
			}
		}

		int alongReadMode = 0;
		int lastPositionInExon = -1;
		for(unsigned int posInAlignment = 0; posInAlignment < readAlignment_exonPositions.size(); posInAlignment++)
		{
			oneExonPosition& thisPosition = readAlignment_exonPositions.at(posInAlignment);
			assert(thisPosition.graphLevel != -1);

			// std::cout << alongReadMode << " " << thisPosition.graphLevel << " " << (int)graphLevel_2_exonPosition.count(thisPosition.graphLevel) << "\n" << std::flush;

			if(graphLevel_2_exonPosition.count(thisPosition.graphLevel))
			{
				// if(!((alongReadMode == 0) || (alongReadMode == 1)))
				// {
					// std::cerr << "alongReadMode" << ": " << alongReadMode << "\n";
					// std::cerr << "alignment.sequence_aligned.substr(cI, 1)" << ": " << alignment.sequence_aligned << "\n";
					// std::cerr << "alignment.graph_aligned.substr(cI, 1)" << ": " << alignment.graph_aligned << "\n";
					// std::cerr << "alignment.graph_aligned_levels)" << ": " << Utilities::join(Utilities::ItoStr(alignment.graph_aligned_levels), ", ") << "\n";
					// std::cerr << "alignment_firstLevel" << ": " << alignment_firstLevel << "\n";
					// std::cerr << "alignment_lastLevel" << ": " << alignment_lastLevel << "\n";
					// std::cerr << "combined_exon_sequences_graphLevels.front()" << ": " << combined_exon_sequences_graphLevels.front() << "\n";
					// std::cerr << "combined_exon_sequences_graphLevels.back()" << ": " << combined_exon_sequences_graphLevels.back() << "\n";
					// for(unsigned int i = 0; i < alignment.graph_aligned_levels.size(); i++)
					// {
						// std::cout << "\t" << i << " " << alignment.graph_aligned_levels.at(i) << " " << (int)graphLevel_2_exonPosition.count(thisPosition.graphLevel) << "\n";
					// }
					// std::cerr << std::flush;
				// }
				if(alongReadMode == 2)
				{
					lastPositionInExon = -1;
				}
				// assert((alongReadMode == 0) || (alongReadMode == 1));


				thisPosition.positionInExon = graphLevel_2_exonPosition.at(thisPosition.graphLevel);
				assert((lastPositionInExon == -1) || ((int)thisPosition.positionInExon == ((int)lastPositionInExon + 1)));
				lastPositionInExon = thisPosition.positionInExon;
				alongReadMode = 1;

				if(!((thisPosition.graphLevel >= combined_exon_sequences_graphLevels_min) && (thisPosition.graphLevel <= combined_exon_sequences_graphLevels_max)))
				{
					std::cout << "! ((thisPosition.graphLevel >= combined_exon_sequences_graphLevels_front) && (thisPosition.graphLevel <= combined_exon_sequences_graphLevels_back))" << "\n";
					std::cout << "\t" << "thisPosition.graphLevel" << ": " << thisPosition.graphLevel << "\n";
					std::cout << "\t" << "combined_exon_sequences_graphLevels_min" << ": " << combined_exon_sequences_graphLevels_min << "\n";
					std::cout << "\t" << "combined_exon_sequences_graphLevels_max" << ": " << combined_exon_sequences_graphLevels_max << "\n";
					std::cout << "\t" << "graphLevel_2_exonPosition.at(thisPosition.graphLevel)" << ": " << graphLevel_2_exonPosition.at(thisPosition.graphLevel) << "\n";
					// std::cout << "\t" << "combined_exon_sequences_graphLevels" << ": " << Utilities::join(Utilities::ItoStr(combined_exon_sequences_graphLevels), " ") << "\n";
					
					std::cout << "\n" << std::flush;
				}
				assert((thisPosition.graphLevel >= combined_exon_sequences_graphLevels_min) && (thisPosition.graphLevel <= combined_exon_sequences_graphLevels_max));
				ret_exonPositions.push_back(thisPosition);
			}
			else
			{
				if(alongReadMode == 1)
				{
					alongReadMode = 2;
				}
			}
		}

		// std::cout << "\n";
	}
}


void HLATyper::oneReadAlignment_2_exonPositions_unpaired(const mapper::reads::verboseSeedChain& alignment, const mapper::reads::oneRead& read, std::vector<oneExonPosition>& ret_exonPositions, int combined_exon_sequences_graphLevels_min, int combined_exon_sequences_graphLevels_max, const std::map<int, unsigned int>& graphLevel_2_exonPosition)
{
	int alignment_firstLevel = alignment.alignment_firstLevel();
	int alignment_lastLevel = alignment.alignment_lastLevel();

	double thisRead_fractionOK = alignmentFractionOK(alignment);
	double pairedRead_fractionOK = -1;

	double thisRead_WeightedCharactersOK = alignmentWeightedOKFraction(read, alignment);
	double pairedRead_WeightedCharactersOK = -1;

	double pairs_strands_OK = 1;
	double pairs_strands_distance = -1;

	// std::cout << "This alignment " << alignment_firstLevel << " - " << alignment_lastLevel << "\n";
	// std::cout << "\tvs combined exon " << combined_exon_sequences_graphLevels.front() << " - " << combined_exon_sequences_graphLevels.back() << "\n\n" << std::flush;


	if(globalVerbose && (read.name.find("HLA-C") != std::string::npos))
	{
		std::cout << "Read" << ": " << read.name << "\n";
		std::cout << "\t" << "alignment_firstLevel" << ": " << alignment_firstLevel << "\n";
		std::cout << "\t" << "alignment_lastLevel" << ": " << alignment_lastLevel << "\n";
		std::cout << "\t" << "combined_exon_sequences_graphLevels_min" << ": " << combined_exon_sequences_graphLevels_min << "\n";
		std::cout << "\t" << "combined_exon_sequences_graphLevels_max" << ": " << combined_exon_sequences_graphLevels_max << "\n";
		std::cout << "\n" << std::flush;
	}


	// if( ((alignment_firstLevel >= combined_exon_sequences_graphLevels_front) && (alignment_firstLevel <= combined_exon_sequences_graphLevels_back)) ||
	//	((alignment_lastLevel >= combined_exon_sequences_graphLevels_front) && (alignment_lastLevel <= combined_exon_sequences_graphLevels_back)) )
	assert((combined_exon_sequences_graphLevels_min <= combined_exon_sequences_graphLevels_max));
	if(Utilities::intervalsOverlap(alignment_firstLevel, alignment_lastLevel, combined_exon_sequences_graphLevels_min, combined_exon_sequences_graphLevels_max))
	{
		std::vector<oneExonPosition> readAlignment_exonPositions;

		int alignmentColumns_oneNonGap = 0;
		for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
		{
			std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
			std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);		
			
			if((sequenceCharacter != "_") || (sequenceCharacter != "_"))
			{
				alignmentColumns_oneNonGap++;
			}
		}
		
		std::vector<int> runningNovelGaps;
		runningNovelGaps.resize(alignment.sequence_aligned.length(), 0);
		{
			int runningNovelGap = 0;
			for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
			{
				std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
				std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
				
				if((graphCharacter != "_") && (sequenceCharacter != "_"))
				{
					runningNovelGap = 0;
				}
				else
				{
					if(!((graphCharacter == "_") && (sequenceCharacter == "_")))
					{
						runningNovelGap++;
					}					
				}
			
				if(runningNovelGap > runningNovelGaps.at(cI))
					runningNovelGaps.at(cI) = runningNovelGap;
			}
			
			runningNovelGap = 0;
			for(int cI = alignment.sequence_aligned.length()-1; cI >= 0; cI--)
			{
				std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
				std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
				
				if((graphCharacter != "_") && (sequenceCharacter != "_"))
				{
					runningNovelGap = 0;
				}
				else
				{
					if(!((graphCharacter == "_") && (sequenceCharacter == "_")))
					{
						runningNovelGap++;
					}					
				}
			
				if(runningNovelGap > runningNovelGaps.at(cI))
					runningNovelGaps.at(cI) = runningNovelGap;
			}		
		}
		
		int indexIntoOriginalReadData = alignment.sequence_begin - 1;; 
		for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
		{
			std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
			std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
			int graphLevel = alignment.graph_aligned_levels.at(cI);
			unsigned char alignmentQualityCharacter_thisPosition = alignment.mapQ_perPosition.at(cI);
			double alignmentQuality_thisPosition = Utilities::PhredToPCorrect(alignmentQualityCharacter_thisPosition);

			if(graphLevel == -1)
			{
				// insertion relative to the graph - we need to extend last character

				assert(graphCharacter == "_");
				assert(sequenceCharacter != "_");

				indexIntoOriginalReadData++;
				int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
				if(alignment.reverse)
				{
					indexIntoOriginalReadData_correctlyAligned = read.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
				}
				assert(indexIntoOriginalReadData_correctlyAligned >= 0);
				assert(indexIntoOriginalReadData_correctlyAligned <(int) read.sequence.length());

				std::string underlyingReadCharacter = read.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
				if(alignment.reverse)
				{
					underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
				}
				assert(underlyingReadCharacter == sequenceCharacter);
				char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

				if(readAlignment_exonPositions.size() > 0)
				{
					readAlignment_exonPositions.back().genotype.append(sequenceCharacter);
					readAlignment_exonPositions.back().alignment_edgelabels.append(graphCharacter);
					readAlignment_exonPositions.back().qualities.push_back(qualityCharacter);
					if(!(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length()))
					{
						assert(readAlignment_exonPositions.back().genotype.length() == (readAlignment_exonPositions.back().qualities.length()+1));
						assert(readAlignment_exonPositions.back().genotype.at(0) == '_');
						readAlignment_exonPositions.back().genotype = readAlignment_exonPositions.back().genotype.substr(1);
						readAlignment_exonPositions.back().alignment_edgelabels = readAlignment_exonPositions.back().alignment_edgelabels.substr(1);
						assert(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length());

						// std::cerr << "readAlignment_exonPositions.back().genotype.length()" << ": " << readAlignment_exonPositions.back().genotype.length() << "\n";
						// std::cerr << "readAlignment_exonPositions.back().qualities.length()" << ": " << readAlignment_exonPositions.back().qualities.length()<< "\n";

						// std::cerr << std::flush;
					}
					assert(readAlignment_exonPositions.back().genotype.length() == readAlignment_exonPositions.back().qualities.length());
					assert(readAlignment_exonPositions.back().alignment_edgelabels.length() == readAlignment_exonPositions.back().genotype.length());
				}
			}
			else
			{
				if(sequenceCharacter != "_")
				{
					indexIntoOriginalReadData++;
					int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
					if(alignment.reverse)
					{
						indexIntoOriginalReadData_correctlyAligned = read.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
					}
					assert(indexIntoOriginalReadData_correctlyAligned >= 0);
					assert(indexIntoOriginalReadData_correctlyAligned < (int)read.sequence.length());

					std::string underlyingReadCharacter = read.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
					if(alignment.reverse)
					{
						underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
					}
					assert(underlyingReadCharacter == sequenceCharacter);

					if(graphCharacter == "_")
					{

						assert(graphLevel != -1);

						char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = sequenceCharacter;
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities.push_back(qualityCharacter);

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = "";
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = read.name;
						
						thisPosition.alignmentColumnsWithAtLeastOneNonGap = alignmentColumns_oneNonGap;
						thisPosition.runningNovelGapEitherDirection = runningNovelGaps.at(cI);
						thisPosition.reverse = alignment.reverse;
						thisPosition.fromFirstRead = alignment.fromFirstRead;
						
						readAlignment_exonPositions.push_back(thisPosition);

					}
					else
					{
						// two well-defined characters
						char qualityCharacter = read.quality.at(indexIntoOriginalReadData_correctlyAligned);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = sequenceCharacter;
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities.push_back(qualityCharacter);

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = "";
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = read.name;

						thisPosition.alignmentColumnsWithAtLeastOneNonGap = alignmentColumns_oneNonGap;
						thisPosition.runningNovelGapEitherDirection = runningNovelGaps.at(cI);
						thisPosition.reverse = alignment.reverse;
						thisPosition.fromFirstRead = alignment.fromFirstRead;
						
						readAlignment_exonPositions.push_back(thisPosition);
					}

				}
				else
				{
					assert(sequenceCharacter == "_");
					if(graphCharacter == "_")
					{
						assert(graphLevel != -1);

						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = "_";
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities = "";

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = "";
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = read.name;

						thisPosition.alignmentColumnsWithAtLeastOneNonGap = alignmentColumns_oneNonGap;
						thisPosition.runningNovelGapEitherDirection = runningNovelGaps.at(cI);
						thisPosition.reverse = alignment.reverse;
						thisPosition.fromFirstRead = alignment.fromFirstRead;
						
						readAlignment_exonPositions.push_back(thisPosition);
					}
					else
					{
						oneExonPosition thisPosition;
						thisPosition.graphLevel = graphLevel;
						thisPosition.genotype = "_";
						thisPosition.alignment_edgelabels = graphCharacter;
						thisPosition.qualities = "";

						thisPosition.thisRead_ID = read.name;
						thisPosition.pairedRead_ID = "";
						thisPosition.thisRead_fractionOK = thisRead_fractionOK;
						thisPosition.pairedRead_fractionOK = pairedRead_fractionOK;
						thisPosition.pairs_strands_OK = pairs_strands_OK;
						thisPosition.pairs_strands_distance = pairs_strands_distance;
						thisPosition.thisRead_WeightedCharactersOK = thisRead_WeightedCharactersOK;
						thisPosition.pairedRead_WeightedCharactersOK = pairedRead_WeightedCharactersOK;

						thisPosition.mapQ = alignment.mapQ;
						thisPosition.mapQ_genomic = alignment.mapQ;
						thisPosition.mapQ_position = alignmentQuality_thisPosition;

						thisPosition.read1_ID = read.name;

						thisPosition.alignmentColumnsWithAtLeastOneNonGap = alignmentColumns_oneNonGap;
						thisPosition.runningNovelGapEitherDirection = runningNovelGaps.at(cI);
						thisPosition.reverse = alignment.reverse;
						thisPosition.fromFirstRead = alignment.fromFirstRead;
						
						readAlignment_exonPositions.push_back(thisPosition);
					}
				}
			}
		}

		int alongReadMode = 0;
		int lastPositionInExon = -1;
		for(unsigned int posInAlignment = 0; posInAlignment < readAlignment_exonPositions.size(); posInAlignment++)
		{
			oneExonPosition& thisPosition = readAlignment_exonPositions.at(posInAlignment);
			assert(thisPosition.graphLevel != -1);

			// std::cout << alongReadMode << " " << thisPosition.graphLevel << " " << (int)graphLevel_2_exonPosition.count(thisPosition.graphLevel) << "\n" << std::flush;

			if(graphLevel_2_exonPosition.count(thisPosition.graphLevel))
			{
				// if(!((alongReadMode == 0) || (alongReadMode == 1)))
				// {
					// std::cerr << "alongReadMode" << ": " << alongReadMode << "\n";
					// std::cerr << "alignment.sequence_aligned.substr(cI, 1)" << ": " << alignment.sequence_aligned << "\n";
					// std::cerr << "alignment.graph_aligned.substr(cI, 1)" << ": " << alignment.graph_aligned << "\n";
					// std::cerr << "alignment.graph_aligned_levels)" << ": " << Utilities::join(Utilities::ItoStr(alignment.graph_aligned_levels), ", ") << "\n";
					// std::cerr << "alignment_firstLevel" << ": " << alignment_firstLevel << "\n";
					// std::cerr << "alignment_lastLevel" << ": " << alignment_lastLevel << "\n";
					// std::cerr << "combined_exon_sequences_graphLevels.front()" << ": " << combined_exon_sequences_graphLevels.front() << "\n";
					// std::cerr << "combined_exon_sequences_graphLevels.back()" << ": " << combined_exon_sequences_graphLevels.back() << "\n";
					// for(unsigned int i = 0; i < alignment.graph_aligned_levels.size(); i++)
					// {
						// std::cout << "\t" << i << " " << alignment.graph_aligned_levels.at(i) << " " << (int)graphLevel_2_exonPosition.count(thisPosition.graphLevel) << "\n";
					// }
					// std::cerr << std::flush;
				// }
				if(alongReadMode == 2)
				{
					lastPositionInExon = -1;
				}
				// assert((alongReadMode == 0) || (alongReadMode == 1));


				thisPosition.positionInExon = graphLevel_2_exonPosition.at(thisPosition.graphLevel);
				assert((lastPositionInExon == -1) || ((int)thisPosition.positionInExon == ((int)lastPositionInExon + 1)));
				lastPositionInExon = thisPosition.positionInExon;
				alongReadMode = 1;

				ret_exonPositions.push_back(thisPosition);
			}
			else
			{
				if(alongReadMode == 1)
				{
					alongReadMode = 2;
				}
			}
		}

		// std::cout << "\n";
	}
}
 
double HLATyper::alignmentWeightedOKFraction(const mapper::reads::oneRead& underlyingRead, const mapper::reads::verboseSeedChain& alignment)
{
	int indexIntoOriginalReadData = alignment.sequence_begin - 1;

	int totalMismatches = 0;
	double weightedMismatches = 0;

	for(unsigned int cI = 0; cI < alignment.sequence_aligned.length(); cI++)
	{
		std::string sequenceCharacter = alignment.sequence_aligned.substr(cI, 1);
		std::string graphCharacter = alignment.graph_aligned.substr(cI, 1);
		int graphLevel = alignment.graph_aligned_levels.at(cI);

		if(sequenceCharacter != "_")
		{
			indexIntoOriginalReadData++;
			int indexIntoOriginalReadData_correctlyAligned = indexIntoOriginalReadData;
			if(alignment.reverse)
			{
				indexIntoOriginalReadData_correctlyAligned = underlyingRead.sequence.length() - indexIntoOriginalReadData_correctlyAligned - 1;
			}
			assert(indexIntoOriginalReadData_correctlyAligned >= 0);
			assert(indexIntoOriginalReadData_correctlyAligned < (int)underlyingRead.sequence.length());;

			std::string underlyingReadCharacter = underlyingRead.sequence.substr(indexIntoOriginalReadData_correctlyAligned, 1);
			if(alignment.reverse)
			{
				underlyingReadCharacter = Utilities::seq_reverse_complement(underlyingReadCharacter);
			}
			
			if(!(underlyingReadCharacter == sequenceCharacter))
			{
				std::cerr << "! (underlyingReadCharacter == sequenceCharacter)" << "\n";
				std::cerr << underlyingReadCharacter << "\n";
				std::cerr << sequenceCharacter << "\n";
				std::cerr << alignment.sequence_aligned << "\n";
				std::cerr << underlyingRead.sequence << "\n";
				std::cerr << std::flush;
			}
			assert(underlyingReadCharacter == sequenceCharacter);

			if(graphCharacter == "_")
			{
				totalMismatches++;
				weightedMismatches++;
			}
			else
			{
				// two well-defined characters
				char qualityCharacter = underlyingRead.quality.at(indexIntoOriginalReadData_correctlyAligned);
				double pCorrect = Utilities::PhredToPCorrect(qualityCharacter);
				assert((pCorrect >= 0) && (pCorrect <= 1));
				if(sequenceCharacter == graphCharacter)
				{
					// match!
				}
				else
				{
					weightedMismatches += pCorrect;
					totalMismatches++;
				}
			}

		}
		else
		{
			assert(sequenceCharacter == "_");
			if(graphCharacter == "_")
			{
				// sequence gap, graph gap - likelihood 1
				assert(graphLevel != -1);
			}
			else
			{
				// sequence gap, graph non gap - deletion in sequence
				totalMismatches++;
				weightedMismatches++;
			}
		}
	}

	double readLength = underlyingRead.sequence.length();
	assert(totalMismatches >= weightedMismatches);

	return  (1.0 - (weightedMismatches / readLength));
};

std::vector<oneExonPosition> HLATyper::removeDoublePositionsFromRead(const std::vector<oneExonPosition>& positions)
{
	std::map<int, std::vector<oneExonPosition>> positions_per_graphLevel;
	std::vector<oneExonPosition> forReturn;

	for(unsigned int i = 0; i < positions.size(); i++)
	{
		const oneExonPosition& p = positions.at(i);
		int graphLevel = p.graphLevel;
		positions_per_graphLevel[graphLevel].push_back(p);
	}

	auto getWorstQuality = [](const std::string& qualities) -> unsigned char {
		assert(qualities.size() > 0);
		unsigned char worstQ;
		for(unsigned int i = 0; i < qualities.size(); i++)
		{
			if((i == 0) || (qualities.at(i) < worstQ))
			{
				worstQ = qualities.at(i);
			}
		}
		return worstQ;
	};

	auto getBestExonPosition = [&](const std::vector<oneExonPosition>& alternatives) -> oneExonPosition {
		assert(alternatives.size() > 0);

		unsigned int bestI;
		unsigned char bestI_quality;

		for(unsigned int i = 0; i < alternatives.size(); i++)
		{
			const oneExonPosition& thisAlternative = alternatives.at(i);
			if(!(((thisAlternative.genotype == "_") || (thisAlternative.qualities.size() > 0))))
			{
				std::cerr << "thisAlternative.qualities.size() == 0" << "\n";
				std::cerr << "thisAlternative.graphLevel: " << thisAlternative.graphLevel << "\n";
				std::cerr << "thisAlternative.genotype: " << thisAlternative.genotype << "\n";
				std::cerr << std::flush;

			}
			assert((thisAlternative.genotype == "_") || (thisAlternative.qualities.size() > 0));
			unsigned char Q = (thisAlternative.genotype == "_") ? 0 : getWorstQuality(thisAlternative.qualities);
			if((i == 0) || (Q > bestI_quality))
			{
				bestI = i;
				bestI_quality = Q;
			}
		}

		return alternatives.at(bestI);
	};

	for(std::map<int, std::vector<oneExonPosition> >::iterator graphLevelIt = positions_per_graphLevel.begin(); graphLevelIt != positions_per_graphLevel.end(); graphLevelIt++)
	{
		oneExonPosition bestPositionAtLevel = getBestExonPosition(graphLevelIt->second);
		forReturn.push_back(bestPositionAtLevel);
	}

	// std::cout << "removeDoublePositionsFromRead(..): Removed " << (positions.size() - forReturn.size()) << " exon positions (from " << positions.size() << " to " << forReturn.size() << ").\n" << std::flush;

	return forReturn;
}


bool HLATyper::can_translateToG_locus(std::string locus)
{
	read_G_alleles();
	assert(alleles_to_G.size() > 0);
	assert(locus.find("*") == std::string::npos);
	return (G_loci.count(locus));
}


std::string HLATyper::translate_allele_list_to_G_allele(const std::vector<std::string>& alleles, bool& ret_perfectly)
{
	read_G_alleles();
	assert(alleles_to_G.size() > 0);

	std::map<std::string, int> g_groups;
	for(auto a : alleles)
	{
		std::vector<std::string> locus_and_allele = Utilities::split(a, "*");
		if(locus_and_allele.size() != 2)
		{
			std::string oneAllele = alleles_to_G.begin()->first;
			throw std::runtime_error("Weird allele: "+a+"; have only things like "+oneAllele);
		}
		std::string locus_without_star = locus_and_allele.at(0);
		assert(can_translateToG_locus(locus_without_star));
		if(! alleles_to_G.count(a))
		{
			std::string oneAllele = alleles_to_G.begin()->first;
			std::cerr << "Warning: Can't G-translate: " << a << "; have only things like " << oneAllele;
			continue;
		}
		assert(alleles_to_G.count(a));
		std::string g_group = alleles_to_G.at(a);
		if(g_groups.count(g_group) == 0)
		{
			g_groups[g_group] = 0;
		}
		g_groups.at(g_group)++;
	}

	if(g_groups.size() == 0)
	{
		std::string forReturn = Utilities::join(alleles, ";");
		ret_perfectly = false;
		return forReturn;
	}
	assert(g_groups.size());

	if(g_groups.size() == 1)
	{
		ret_perfectly = true;
		return g_groups.begin()->first;
	}
	else
	{
		ret_perfectly = false;

		std::vector<std::string> possible_groups_sorted = Utilities::get_map_keys_sorted_by_value(g_groups);
		if(possible_groups_sorted.size() > 1)
		{
			assert(g_groups.at(possible_groups_sorted.at(0)) >= g_groups.at(possible_groups_sorted.at(1)));
		}

		return possible_groups_sorted.at(0);
	}
}

void HLATyper::read_G_alleles()
{
	if(alleles_to_G.size() == 0)
	{
		std::vector<std::string> _test_split = Utilities::split("A*;01:01:22;", ";");
		assert(_test_split.size() == 3);

		std::string g_filename = "hla_nom_g.txt";
		std::ifstream gStream;
		gStream.open(g_filename.c_str());
		if(! gStream.is_open())
		{
			throw std::runtime_error("Can't open file " + g_filename + " - are you executing me from the right directory?");
		}

		std::string line;
		while(gStream.good())
		{
			std::getline(gStream, line);
			Utilities::eraseNL(line);
			if(line.length() == 0)
				continue;

			if(line.substr(0, 1) == "#")
				continue;

			std::vector<std::string> components = Utilities::split(line, ";");
			std::string locus_with_star = components.front();
			assert(locus_with_star.back() == '*');

			std::string locus_without_star = locus_with_star.substr(0, locus_with_star.size() - 1);
			G_loci.insert(locus_without_star);

			std::string g_code;
			if(components.back() != "")
			{
				g_code = components.back();
			}
			else
			{
				assert(components.size() == 3);
				g_code = components.at(1);
			}

			g_code = locus_with_star + g_code;
			std::vector<std::string> differentAlleles = Utilities::split(components.at(1), "/");

			for(auto a : differentAlleles)
			{
				std::string allele_with_locus = locus_with_star + a;
				alleles_to_G[allele_with_locus] = g_code;
			}
		}
	}
}


int HLATyper::kMerSeq_is_kMerKey(const std::string& kMer)
{
	char firstC = kMer.at(0);
	char lastC = kMer.at(kMer.size() - 1);
	char lastC_reverse = Utilities::reverse_char_nucleotide(lastC);

	if(firstC < lastC_reverse)
	{
		return 1;
	}
	else if(firstC == lastC_reverse)
	{
		return 0;
	}
	else
	{
		assert(firstC > lastC_reverse);
		return - 1;
	}
}

std::string HLATyper::kMer_canonical_representation(const std::string& kMer)
{
	bool ignore;
	return kMer_canonical_representation(kMer, ignore);
}

std::string HLATyper::kMer_canonical_representation(const std::string& kMer, bool& was_inverted)
{
	if(kMerSeq_is_kMerKey(kMer) > 0)
	{
		was_inverted = false;
		return kMer;
	}

	std::string kMer_revCmp = Utilities::seq_reverse_complement(kMer);
	if(kMer_revCmp < kMer)
	{
		was_inverted = true;
		return kMer_revCmp;
	}
	else
	{
		was_inverted = false;
		return kMer;
	}
}

double HLATyper::simpleChiSq(std::vector<double> observed, std::vector<double> expected)
{
        double statistic = 0;
        assert(observed.size() == expected.size());

        double observed_sum = 0;
        for(unsigned int i = 0; i < observed.size(); i++)
        {
                if(observed.at(i) < 0)
                {
                        std::cerr << "Error: the " << i << "-th element of observed is < 0\n" << std::flush;
                }
                assert(observed.at(i) >= 0);
                observed_sum += observed.at(i);
        }
        assert(observed_sum > 0);

        for(unsigned int i = 0; i < observed.size(); i++)
        {
                // observed.at(i) = observed.at(i)/observed_sum;
        }

        double expected_sum = 0;
        for(unsigned int i = 0; i < expected.size(); i++)
        {
                assert(expected.at(i) > 0);
                expected_sum += expected.at(i);
        }
        assert(expected_sum > 0);

        for(unsigned int i = 0; i < observed.size(); i++)
        {
                // expected.at(i) = expected.at(i)/expected_sum;
        }

        for(unsigned int i = 0; i < observed.size(); i++)
        {
                assert(expected.at(i) >= 0);
                double thisSummand = pow((observed.at(i) - expected.at(i)), 2)/expected.at(i);
                if(!(thisSummand >= 0))
                {
                        std::cerr << "thisSummand" << ": " << thisSummand << "\n";
                        std::cerr << "observed.at(i)" << ": " << observed.at(i) << "\n";
                        std::cerr << "expected.at(i)" << ": " << expected.at(i) << "\n";
                        std::cerr << "pow((observed.at(i) - expected.at(i)), 2)" << ": " << pow((observed.at(i) - expected.at(i)), 2) << "\n";
                        std::cerr << std::flush;
                }
                assert(thisSummand >= 0);
                statistic += thisSummand;
        }

        assert(statistic >= 0);
        // std::cout << "ChiSq statistic: " << statistic << "\n";

        int df = observed.size() - 1;
        assert(df > 0);


        boost::math::chi_squared chiSqdist(1);
        double pValue = 1 -  boost::math::cdf(chiSqdist, statistic);
        assert((pValue >= 0) && (pValue <= 1));
        return pValue;
}


} /* namespace hla */
