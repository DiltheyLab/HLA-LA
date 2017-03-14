/*
 * LocusCodeAllocation.cpp
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#include "LocusCodeAllocation.h"
#include "../Utilities.h"
#include <boost/algorithm/string.hpp>
#include <assert.h>
#include <set>
#include <algorithm>
#include <iostream>
#include <exception>
#include <stdexcept>

LocusCodeAllocation::LocusCodeAllocation() {

}

LocusCodeAllocation::~LocusCodeAllocation() {
}

void LocusCodeAllocation::removeLocus(string locus)
{
	coded_values.erase(locus);
	coded_values_rev.erase(locus);
	restrictedHLAcache.erase(locus);
}

std::string LocusCodeAllocation::deCode(std::string locus, unsigned char value)
{
	if(coded_values_rev.count(locus) == 0)
	{
		throw std::runtime_error("Non-assigned locus "+locus);
		return 0;
	}
	if(coded_values_rev[locus].count(value) > 0)
	{
		return coded_values_rev[locus][value];
	}
	else
	{
		throw std::runtime_error("Non-assigned allele for "+locus+ " value: "+Utilities::ItoStr((unsigned int)value));
		return 0;
	}
}

void LocusCodeAllocation::takeLocusData(LocusCodeAllocation& otherAllocation)
{
	for(map<string, map<string, unsigned char> >::iterator locusIt = otherAllocation.coded_values.begin(); locusIt != otherAllocation.coded_values.end(); locusIt++)
	{
		string locus = locusIt->first;
		if(coded_values.count(locus) != 0)
		{
			cerr << "Want to take over other CODE, which has data for locus which is already present in this object. Locus: " << locus  << "\n";
		}
		assert(coded_values.count(locus) == 0);
		coded_values[locus] = otherAllocation.coded_values[locus];
		coded_values_rev[locus] = otherAllocation.coded_values_rev[locus];

	}
}

unsigned char LocusCodeAllocation::invert(string locus, unsigned char emission)
{
	string nucleotide = deCode(locus, emission);
	string invertedNucleotide;
	if(nucleotide == "A")
	{
		invertedNucleotide = "T";
	}
	else if(nucleotide == "C")
	{
		invertedNucleotide = "G";
	}
	else if(nucleotide == "G")
	{
		invertedNucleotide = "C";
	}
	else if(nucleotide == "T")
	{
		invertedNucleotide = "A";
	}
	else
	{
		throw std::runtime_error("Cannot invert this nucleotide: "+nucleotide);
		assert(1==0);
	}
	return doCode(locus, invertedNucleotide);
}
unsigned char LocusCodeAllocation::doCode(std::string locus, std::string value)
{
	if((value == "?") || (value == "????"))
	{
		return '0';
	}

	if(coded_values[locus].count(value) > 0)
	{
		return coded_values[locus][value];
	}
	else
	{
		unsigned char start = '0';
		int elements_in = coded_values[locus].size();
		unsigned char new_index = (unsigned char) start+elements_in+1;
		if(new_index > 250)
		{
			throw std::runtime_error("Too many alleles");
		}

		coded_values[locus][value] = new_index;
		coded_values_rev[locus][new_index] = value;

		return new_index;
	}
}

vector<string> LocusCodeAllocation::getLoci()
{
	vector<string> loci;
	for(map<string, map<string, unsigned char> >::iterator StringCodes = coded_values.begin(); StringCodes != coded_values.end(); StringCodes++)
	{
		string locus = StringCodes->first;
		loci.push_back(locus);
	}
	return loci;	
}

vector<string> LocusCodeAllocation::getSNPs()
{
	vector<string> allLoci = getLoci();
	vector<string> forReturn;
	for(vector<string>::iterator locusIt = allLoci.begin(); locusIt != allLoci.end(); locusIt++)
	{
		string locus = *locusIt;
		if(locus.substr(0, 2) == "rs")
		{
			forReturn.push_back(locus);
		}
	}
	return forReturn;
}

bool LocusCodeAllocation::knowCode(string locus, unsigned char code)
{
	if(coded_values_rev.count(locus) == 0)
	{
		return false;
	}

	return (coded_values_rev[locus].count(code) > 0);
}

bool LocusCodeAllocation::knowAllele(string locus, string allele)
{
	if(coded_values.count(locus) == 0)
	{
		return false;
	}

	return (coded_values[locus].count(allele) > 0);
}


vector<unsigned char> LocusCodeAllocation::getAlleles(string locus)
{
	vector<unsigned char> alleles;
	for(map<unsigned char, string>::iterator CodeAllele = coded_values_rev[locus].begin(); CodeAllele != coded_values_rev[locus].end(); CodeAllele++)
	{
		unsigned char allele = CodeAllele->first;
		assert(allele != '0');
		alleles.push_back(allele);
	}
	return alleles;
}



bool LocusCodeAllocation::locusIsHLA(string locus)
{
	if(locus.length() < 3)
		return false;

	return ((locus.substr(0, 3) == "HLA") || (locus.substr(0, 3) == "KIR"));
}

bool LocusCodeAllocation::allele4D(string locus, unsigned char allele)
{
	assert(locusIsHLA(locus) == true);
	string actualType = deCode(locus, allele);
	assert(actualType.length() == 4);
	return (actualType.substr(2,2) != "00");
}

string LocusCodeAllocation::alleleString2D(string allele)
{
	assert(allele.length() == 4);
	return allele.substr(0,2)+"00";
}

vector<unsigned char> LocusCodeAllocation::getAlleles4DOnly(string locus)
{
	assert(locusIsHLA(locus) == true);
	vector<unsigned char> all_alleles = getAlleles(locus);
	vector<unsigned char> alleles_4D;
	for(unsigned int i = 0; i < all_alleles.size(); i++)
	{
		if(allele4D(locus, all_alleles.at(i)))
		{
			alleles_4D.push_back(all_alleles.at(i));
		}
	}
	return alleles_4D;
}

set<unsigned char> LocusCodeAllocation::getAlleles_2DrestrictOnSingletons (string locus)
{
	if(restrictedHLAcache.count(locus) > 0)
	{
		return restrictedHLAcache[locus];
	}
	else
	{
		assert(locusIsHLA(locus) == true);
	
		set<string> alleles_2D_implied_by_4D;
		vector<unsigned char> alleles_4D = getAlleles4DOnly(locus);
		for(unsigned int i = 0; i < alleles_4D.size(); i++)
		{
			string allele_4D_string = deCode(locus, alleles_4D.at(i));
			string allele_2D_string = alleleString2D(allele_4D_string);
			alleles_2D_implied_by_4D.insert(allele_2D_string);
		}
	
		set<unsigned char> forReturn;
		vector<unsigned char> alleles_all = getAlleles(locus);
		for(unsigned int i = 0; i < alleles_all.size(); i++)
		{
			unsigned char allele = alleles_all.at(i);
			if(allele4D(locus, allele))
			{
				forReturn.insert(allele);
			}
			else
			{
				string alleleString = deCode(locus, allele);
				assert(alleleString.length() == 4);
				if(alleles_2D_implied_by_4D.count(alleleString) == 0)
				{
					forReturn.insert(allele);
				}
			}
		}
		
		restrictedHLAcache[locus] = forReturn;
		
		return forReturn;
	}
}

vector<string> LocusCodeAllocation::serializeIntoVector()
{
	vector<string> forReturn;
	for(map<string, map<string, unsigned char> >::iterator locusIt = coded_values.begin(); locusIt != coded_values.end(); locusIt++)
	{
		string locus = locusIt->first;
		for(map<string, unsigned char>::iterator codeIt = coded_values[locus].begin(); codeIt != coded_values[locus].end(); codeIt++)
		{
			string originalAllele = codeIt->first;
			unsigned char codedAllele = codeIt->second;

			vector<string> lineVector;
			lineVector.push_back(locus);
			lineVector.push_back(originalAllele);
			lineVector.push_back(Utilities::ItoStr((unsigned int)codedAllele));

			forReturn.push_back(boost::algorithm::join(lineVector, separatorForSerialization));
		}
	}

	return forReturn;
}

void LocusCodeAllocation::readFromVector(vector<string> lines)
{
	for(unsigned int i = 0; i < lines.size(); i++)
	{
		string line = lines.at(i);
		vector<string> lineFields;
		boost::iter_split(lineFields, line, boost::first_finder(separatorForSerialization, boost::is_iequal()));

		if(lineFields.size() != 3)
		{
			throw std::runtime_error("Cannot read CODE from line, expect 3 fields, got" + Utilities::ItoStr(lineFields.size())+ "! Line: "+line);
		}
		string locus = lineFields.at(0);
		string originalAllele = lineFields.at(1);
		string integerString = lineFields.at(2);
		int codedChar = Utilities::StrtoI(integerString);
		if((codedChar < 0) or (codedChar > 250))
		{
			throw std::runtime_error("Weird codedChar value: cannot convert back! "+integerString);
		}
		unsigned char codedAllele = (unsigned char) codedChar;

		coded_values[locus][originalAllele] = codedAllele;
		coded_values_rev[locus][codedAllele] = originalAllele;
	}
}

