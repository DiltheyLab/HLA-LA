/*
 * HaplotypePanel.cpp
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "assert.h"

#include "HaplotypePanel.h"
#include "../Utilities.h"

using namespace std;

string haplotype_field_separator = "_____";

HaplotypePanel::HaplotypePanel() {

}

void HaplotypePanel::addString(const std::vector<std::string>& loci, std::string ID, std::string haplotype)
{

	if(LocusPositions.size() == 0)
	{
		for(unsigned int i = 0; i < loci.size(); i++)
		{
			LocusPositions[loci.at(i)] = i;
		}
	}
	else
	{
		assert(haplotype.length() == LocusPositions.size());
		for(unsigned int i = 0; i < loci.size(); i++)
		{
			assert(LocusPositions.at(loci.at(i)) == i);
		}
	}

	assert(loci.size() == haplotype.length());
	for(unsigned int locusI = 0; locusI < haplotype.length(); locusI++)
	{
		std::string locusID = loci.at(locusI);
		assert(LocusPositions.count(locusID));
		unsigned char locusChar = haplotype.at(locusI);
		HaplotypesByLoci[locusID].push_back(locusChar);
	}
	assert(std::find(HaplotypeIDs.begin(), HaplotypeIDs.end(), ID) == HaplotypeIDs.end());
	HaplotypeIDs.push_back(ID);
}

std::map<std::string, std::string> HaplotypePanel::getAllHaplotypes()
{
	std::map<std::string, std::string> forReturn;
	std::vector<std::string> orderedLoci = getOrderedLoci();
	for(unsigned int haplotypeI = 0; haplotypeI < HaplotypeIDs.size(); haplotypeI++)
	{
		std::string haplotypeID = HaplotypeIDs.at(haplotypeI);
		std::string haplotype = getIndividualStringHaplotype(orderedLoci, haplotypeID, false);
		forReturn[haplotypeID] = haplotype;
	}
	return forReturn;
}

std::vector<std::string> HaplotypePanel::getOrderedLoci()
{
	std::vector<std::string> forReturn;
	forReturn = getUnorderedLoci();
	std::sort(forReturn.begin(), forReturn.end(), [&](const std::string& a, const std::string& b){
		return (LocusPositions.at(a) < LocusPositions.at(b));
	});
	if(forReturn.size() > 1)
	{
		for(unsigned int i = 0; i < (forReturn.size() - 1); i++)
		{
			std::string l1 = forReturn.at(i);
			std::string l2 = forReturn.at(i+1);
			assert(LocusPositions.at(l1) < LocusPositions.at(l2));
		}
	}
	return forReturn;
}
void HaplotypePanel::readFromFile(string filename, string positions)
{
	/*
	ifstream positionsStr;
	positionsStr.open (positions.c_str(), ios::in);
	if(positionsStr.is_open())
	{
		string line;
		while(positionsStr.good())
		{
			getline (positionsStr, line);
			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			vector<string> fields = Utilities::split(line, ' ');
			if((fields.size() != 2) and (fields.size() != 3))
				errEx("Strange format in positions file. Expect two/three fields: field name and position (integer), separated by a whitespace (strand optional). Problem occured in line: "+line);

			string locusID = fields.at(0);
			string posStr = fields.at(1);
			int pos = Utilities::StrtoI(posStr);

			LocusPositions[locusID] = pos;
			if(fields.size() == 3)
			{
				LocusStrands[locusID] = fields.at(2);
			}
		}
		positionsStr.close();
	}
	else
	{
		errEx("Cannot open positions file: "+positions);
	}

	ifstream haplotypes;
	haplotypes.open (filename.c_str(), ios::in);
	if (haplotypes.is_open())
	{
		string line;
		if(! haplotypes.good())
		{
			errEx("Read error: first line!");
		}
		getline (haplotypes, line);
		Utilities::eraseNL(line);

		vector<string> header_fields = Utilities::split(line, ' ');

		if(header_fields.at(0) != "IndividualID")
			errEx("File error: no IndividualID!");

		bool useChromosome = true;
		if(header_fields.at(1) != "Chromosome")
			useChromosome = false;

		int eraseFields = 2;
		if(! useChromosome)
			eraseFields = 1;

		vector<string> locus_IDs = header_fields;
		locus_IDs.erase(locus_IDs.begin(), locus_IDs.begin()+eraseFields);

		int lineI = 0;
		while(haplotypes.good())
		{
			getline (haplotypes, line);
			lineI++;
			Utilities::eraseNL(line);
			vector<string> line_fields = Utilities::split(line, ' ');

			if(line_fields.size() == 0)
				continue;

			if(line_fields.size() != header_fields.size())
				errEx("Wrong field size: " + Utilities::ItoStr(line_fields.size()) + " vs expected " + Utilities::ItoStr(header_fields.size()) + " line " + Utilities::ItoStr(lineI));

			string hID;
			if(useChromosome)
			{
				hID = line_fields.at(0)+haplotype_field_separator+line_fields.at(1);
			}
			else
			{
				hID = line_fields.at(0)+haplotype_field_separator+"1";
			}
			HaplotypeIDs.push_back(hID);

			for(unsigned int i = eraseFields; i < line_fields.size(); i++)
			{
				string locusID = locus_IDs.at(i-eraseFields);
				unsigned char thisV = CODE.doCode(locusID, line_fields[i]);
				HaplotypesByLoci[locusID].push_back(thisV);
			}

		}
		haplotypes.close();
	}
	else
	{
		errEx("Could not open file: "+filename);
	}

	*/
	assert( 1 == 0 );
}

void HaplotypePanel::writeToFile(string filename)
{
	/*
	ofstream haplotypes;
	haplotypes.open (filename.c_str(), ios::out | ios::trunc);
	string sep = " ";
	if (haplotypes.is_open())
	{
		vector<string> fieldNames = getLoci();

		haplotypes << "IndividualID" << sep << "Chromosome" << sep << Utilities::join(fieldNames, " ") << "\n";

		int haplotype_counter = -1;
		for(vector<string>::iterator haplotypeID = HaplotypeIDs.begin(); haplotypeID != HaplotypeIDs.end(); haplotypeID++)
		{
			haplotype_counter++;
			size_t found_separator = haplotypeID->find(haplotype_field_separator);
			if(found_separator == string::npos)
			{
				errEx("Cannot find separator!");
			}
			string ID = haplotypeID->substr(0, (int) found_separator);
			size_t start_chromosome = found_separator + haplotype_field_separator.length();
			string chromosome = haplotypeID->substr(start_chromosome, haplotypeID->length() - start_chromosome);
			vector<string> fields_to_print;
			fields_to_print.push_back(ID);
			fields_to_print.push_back(chromosome);
			for(vector<string>::iterator fieldName = fieldNames.begin(); fieldName != fieldNames.end(); fieldName++)
			{
				unsigned char value = HaplotypesByLoci[*fieldName][haplotype_counter];
				string decoded = CODE.deCode(*fieldName, value);
				fields_to_print.push_back(decoded);
			}

			haplotypes << Utilities::join(fields_to_print, " ") << "\n";
		}

		haplotypes.close();
	}
	else
	{
		errEx("Could not open file: "+filename);
	}
	*/
	assert( 1 == 0 );
}

std::string HaplotypePanel::getIndividualStringHaplotype(const vector<string>& loci, string haplotypeID, bool pseudoTwoHaplotypes)
{
	basic_string<unsigned char> S = getIndividualHaplotype(loci, haplotypeID, pseudoTwoHaplotypes);
	std::string forReturn;
	forReturn.reserve(S.length());
	for(unsigned int i = 0; i < S.length(); i++)
	{
		forReturn.push_back(S.at(i));
	}
	return forReturn;
}

basic_string<unsigned char> HaplotypePanel::getIndividualHaplotype(const vector<string>& loci, string haplotypeID, bool pseudoTwoHaplotypes)
{
	basic_string<unsigned char> forReturn;

	if(haplotypeID_cache.count(haplotypeID) == 0)
	{
		int found_id = -1;
		for(unsigned int i = 0; i < HaplotypeIDs.size(); i++)
		{
			if(HaplotypeIDs.at(i) == haplotypeID)
			{
				found_id = i;
			}
		}
		assert(found_id != -1);
		haplotypeID_cache[haplotypeID] = found_id;
	}

	assert(HaplotypeIDs.at(haplotypeID_cache[haplotypeID]) == haplotypeID);

	unsigned char missingDataChar = '0';

	for(unsigned int l = 0; l < loci.size(); l++)
	{
		string locus = loci.at(l);
		if(HaplotypesByLoci.count(locus) == 0)
		{
			cout << "Problem: requested locus " << locus << ", but not present in haplotype panel -- add as missing data!\n";
		}
		assert(HaplotypesByLoci.count(locus) > 0);
		basic_string<unsigned char> toAppend = HaplotypesByLoci[locus].substr(haplotypeID_cache[haplotypeID], 1);
		forReturn.append(toAppend);
		if(pseudoTwoHaplotypes == true)
		{
			forReturn.push_back(missingDataChar);
		}
	}

	if(pseudoTwoHaplotypes == false)
	{
		assert(forReturn.size() == loci.size());
	}
	else
	{
		assert(forReturn.size() == (loci.size()*2));
	}

	return forReturn;
}

vector<string> HaplotypePanel::getUnorderedLoci()
{
	vector<string> fieldNames;
	for(map<string, basic_string<unsigned char> >::iterator iter = HaplotypesByLoci.begin(); iter != HaplotypesByLoci.end(); ++iter)
	{
		fieldNames.push_back(iter->first);
	}
	return fieldNames;
}

