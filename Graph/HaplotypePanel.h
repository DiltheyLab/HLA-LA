/*
 * HaplotypePanel.h
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef HAPLOTYPEPANEL_H_
#define HAPLOTYPEPANEL_H_
#include <map>
#include <string>
#include <vector>
#include <string>

extern std::string haplotype_field_separator;

using namespace std;

class HaplotypePanel {
protected:
	map<string, int> haplotypeID_cache;

public:
	HaplotypePanel();

	void readFromFile(string filename, string positions);
	void writeToFile(string filename);

	vector<string> getUnorderedLoci();

	vector<string> HaplotypeIDs;
	map<string, basic_string<unsigned char> > HaplotypesByLoci;
	map<string, int> LocusPositions;
	map<string, string> LocusStrands;

	std::map<std::string, std::string> getAllHaplotypes();

	basic_string<unsigned char> getIndividualHaplotype(const vector<string>& loci, string haplotypeID, bool pseudoTwoHaplotypes);
	std::string getIndividualStringHaplotype(const vector<string>& loci, string haplotypeID, bool pseudoTwoHaplotypes);

	std::vector<std::string> getOrderedLoci();

	void addString(const std::vector<std::string>& loci, std::string ID, std::string haplotype);
};

#endif /* HAPLOTYPEPANEL_H_ */
