/*
 * LocusCodeAllocation.h
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef LOCUSCODEALLOCATION_H_
#define LOCUSCODEALLOCATION_H_

#include <string>
#include <map>
#include <vector>
#include <set>

using namespace std;

struct MutateAwayTo {
	unsigned char targetAllele;
	double p;
};

extern string separatorForSerialization;

class LocusCodeAllocation
{
protected:
	map<string, map<string, unsigned char> > coded_values;
	map<string, map<unsigned char, string> > coded_values_rev;
	map<string, set<unsigned char> > restrictedHLAcache;
	
public:
	LocusCodeAllocation();
	virtual ~LocusCodeAllocation();
	unsigned char doCode(string locus, string value);
	unsigned char invert(string locus, unsigned char emission);

	void removeLocus(string locus);

	string deCode(string locus, unsigned char value);
	vector<unsigned char> getAlleles(string locus);

	bool knowCode(string locus, unsigned char code);
	bool knowAllele(string locus, string allele);
	vector<string> getLoci();
	vector<string> getSNPs();
	
	set<unsigned char> getAlleles_2DrestrictOnSingletons (string locus);
	vector<unsigned char> getAlleles4DOnly(string locus);
	string alleleString2D(string allele);
	bool allele4D(string locus, unsigned char allele);
	bool locusIsHLA(string locus);

	vector<string> serializeIntoVector();
	void readFromVector(vector<string> lines);

	void takeLocusData(LocusCodeAllocation& otherAllocation);

};

#endif /* LOCUSCODEALLOCATION_H_ */
