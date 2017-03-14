/*
 * protoSeeds.cpp
 *
 *  Created on: 21.09.2015
 *      Author: AlexanderDilthey
 */

#include "protoSeeds.h"

#include <iostream>
#include <exception>
#include <stdexcept>
#include "../../Utilities.h"

namespace mapper {

namespace reads {
protoSeeds::protoSeeds() {
	read1_havePrimary = false;
	read2_havePrimary = false;
}

void protoSeeds::takeAlignment(const BamTools::BamAlignment& a, int which, std::string referenceID, int reference2level_offset_0based, int whichReader)
{
	if(which == 1)
	{
		read1_alignments.push_back(make_tuple(referenceID, reference2level_offset_0based, a, whichReader));
		read1_havePrimary = (read1_havePrimary || a.IsPrimaryAlignment());
	}
	else
	{
		assert(which == 2);
		read2_alignments.push_back(make_tuple(referenceID, reference2level_offset_0based, a, whichReader));
		read2_havePrimary = (read2_havePrimary || a.IsPrimaryAlignment());
	}
}

void protoSeeds::refreshPrimaryStatus()
{
	read1_havePrimary = false;
	for(auto alignment : read1_alignments)
	{
		read1_havePrimary = (read1_havePrimary || (std::get<2>(alignment).IsPrimaryAlignment()));
	}
	read2_havePrimary = false;
	for(auto alignment : read2_alignments)
	{
		read2_havePrimary = (read2_havePrimary || (std::get<2>(alignment).IsPrimaryAlignment()));
	}
}
void protoSeeds::printAlignmentInfo(size_t i, const BamTools::BamAlignment& al)
{
	std::string CIGARstring;
	std::vector< BamTools::CigarOp > CIGAR_Compressed = al.CigarData;
	std::vector<char> CIGAR;
	CIGAR.reserve(al.AlignedBases.length()*1.1);
	for(unsigned int cigarI = 0; cigarI < CIGAR_Compressed.size(); cigarI++)
	{
		BamTools::CigarOp thisOp = CIGAR_Compressed.at(cigarI);
		if(thisOp.Type != 'P')
		{
			for(unsigned int thisOpChar = 0; thisOpChar < thisOp.Length; thisOpChar++)
			{
				CIGAR.push_back(thisOp.Type);
			}
		}

		std::string thisCigarOP = Utilities::ItoStr(thisOp.Length);
		thisCigarOP.push_back(thisOp.Type);
		CIGARstring.append(thisCigarOP);
	}

	std::cerr << "\t" << i << "\n";
	std::cerr << "\t\t" << "Name: " << al.Name << "\n";
	std::cerr << "\t\t" << "Primary: " << al.IsPrimaryAlignment() << "\n";
	std::cerr << "\t\t" << "First Mate: " << al.IsFirstMate() << "\n";
	std::cerr << "\t\t" << "IsReverseStrand: " << al.IsReverseStrand() << "\n";
	std::cerr << "\t\t" << "CIGAR: " << CIGARstring << "\n";
	std::cerr << "\t\t" << "Ref ID: " << al.RefID << "\n";
	std::cerr << "\t\t\t" << "Start position: " << al.Position << "\n";
	std::cerr << "\t\t\t" << "End position: " << al.GetEndPosition() << "\n";
	std::cerr << "\t\t\t" << "L: " << (al.GetEndPosition() - al.Position) << "\n";
	std::cerr << "\t\t" << "QueryBases.length(): " << al.QueryBases.length() << "\n";
	std::cerr << "\t\t" << "AlignedBases.length(): " << al.AlignedBases.length() << "\n";
	std::cerr << "\t\t" << "QueryBases: " << al.QueryBases << "\n";
	std::cerr << "\t\t" << "AlignedBases: " << al.AlignedBases << "\n";	
	//std::cerr << "\t\t" << "ErrorString.length(): " << al.ErrorString.length() << "\n";

	double mapQ = Utilities::PhredToPCorrect(al.MapQuality+33);
	assert(mapQ >= 0);
	assert(mapQ <= 1);
	std::cerr << "\t\t" << "Mapping quality: " << mapQ << "\n";

	std::vector<std::string> flags = {
			"template having multiple segments in sequencing",
			"each segment properly aligned according to the aligner",
			"segment unmapped",
			"next segment in the template unmapped",
			"SEQ being reverse complemented",
			"SEQ of the next segment in the template being reverse complemented",
			"the first segment in the template",
			"the last segment in the template",
			"secondary alignment",
			"not passing fiters, such as platform/vendor quality controls",
			"PCR or optical duplicate",
			"supplementary alignment"
	};

	std::vector<uint32_t> bits_flags = {
			1,
			2,
			4,
			8,
			16,
			32,
			64,
			128,
			256,
			512,
			1024,
			2014
	};

	assert(flags.size() == bits_flags.size());

	std::cerr << "\t\t" << "Flags:\n";

	for(unsigned int flagI = 0; flagI < flags.size(); flagI++)
	{
		uint32_t thisFlag_set = (al.AlignmentFlag & bits_flags.at(flagI));
		std::string isSet = (thisFlag_set == 0) ? "0" : "1";
		std::cerr << "\t\t\t" << flags.at(flagI) << "\t" << isSet << "\n";
	}

	std::vector<char> atoZ;
	for(char C = 'a'; C <= 'z'; C++)
	{
		atoZ.push_back(C);
	}

	for(char C = 'A'; C <= 'Z'; C++)
	{
		atoZ.push_back(C);
	}

	std::vector<char> atoZ_0to9 = atoZ;
	for(char C = '0'; C <= '9'; C++)
	{
		atoZ_0to9.push_back(C);
	}

	std::vector<std::string> tagNames;
	for(unsigned int i1 = 0; i1 < atoZ.size(); i1++)
	{
		for(unsigned int i2 = 0; i2 < atoZ_0to9.size(); i2++)
		{
			std::string tagName;
			tagName.push_back(atoZ.at(i1));
			tagName.push_back(atoZ_0to9.at(i2));
			tagNames.push_back(tagName);
		}
	}


	std::cerr << "\t\tTags:\n" << std::flush;
	for(unsigned int tagI = 0; tagI < tagNames.size(); tagI++)
	{
		std::string tagName = tagNames.at(tagI);
		if(al.HasTag(tagName))
		{
			char type;
			assert(al.GetTagType(tagName, type));

			std::cerr << "\t\t\t" << tagName << ": " << std::flush;
			if(type == BamTools::Constants::BAM_TAG_TYPE_ASCII)
			{
				char v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_ASCII)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_INT8)
			{
				int32_t v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_INT8)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_UINT8)
			{
				u_int32_t v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_UINT8)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_INT16)
			{
				int32_t v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_INT16)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_UINT16)
			{
				u_int32_t v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_UINT16)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_INT32)
			{
				int32_t v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_INT32)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_UINT32)
			{
				u_int32_t v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_UINT32)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_FLOAT)
			{
				double v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_FLOAT)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_STRING)
			{
				std::string v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_STRING)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_HEX)
			{
				char v;
				assert(al.GetTag(tagName, v));
				std::cerr << v << " (BAM_TAG_TYPE_HEX)";
			}
			else if(type == BamTools::Constants::BAM_TAG_TYPE_ARRAY)
			{
				throw std::runtime_error("BAM_TAG_TYPE_ARRAY not implemented!");
			}
			else
			{
				std::string typeAsString;
				typeAsString.push_back(type);
				throw std::runtime_error("BAM alignment printer: unknown TAG type: "+typeAsString);
			}

			std::cerr << "\n" << std::flush;
		}
	}


	std::cerr << "\n" << std::flush;

};

size_t protoSeeds::read1_getPrimaryAlignmentI() const
{
	bool found = false;
	size_t forReturn;


	for(size_t i = 0; i < read1_alignments.size(); i++)
	{
		if(std::get<2>(read1_alignments.at(i)).IsPrimaryAlignment())
		{
			/*
			if(found)
			{
				std::cerr << "It seems that there is more than one primary-read alignment for read 1" << "\n";
				std::cerr << "Read 1:" << "\n";
				for(size_t i2 = 0; i2 < read1_alignments.size(); i2++)
				{
					printAlignmentInfo(i2, read1_alignments.at(i2).second);
				}
				std::cerr << "Read 2:" << "\n";
				for(size_t i2 = 0; i2 < read2_alignments.size(); i2++)
				{
					printAlignmentInfo(i2, read2_alignments.at(i2).second);
				}
				std::cerr << std::flush;
			}
			assert(! found);
			*/
			found = true;
			forReturn = i;
			break;
		}
	}
	assert(found);
	return forReturn;
}

size_t protoSeeds::read2_getPrimaryAlignmentI() const
{
	bool found = false;
	size_t forReturn;
	for(size_t i = 0; i < read2_alignments.size(); i++)
	{
		if(std::get<2>(read2_alignments.at(i)).IsPrimaryAlignment())
		{
			/*
			if(found)
			{
				std::cerr << "It seems that there is more than one primary-read alignment for read 2" << "\n";
				std::cerr << "Read 1:" << "\n";
				for(size_t i2 = 0; i2 < read1_alignments.size(); i2++)
				{
					printAlignmentInfo(i2, read1_alignments.at(i2).second);
				}
				std::cerr << "Read 2:" << "\n";
				for(size_t i2 = 0; i2 < read2_alignments.size(); i2++)
				{
					printAlignmentInfo(i2, read2_alignments.at(i2).second);
				}
				std::cerr << std::flush;
			}

			assert(! found);
			*/
			found = true;
			forReturn = i;
			break;
		}
	}
	assert(found);
	return forReturn;
}

void protoSeeds::printDebug(const BamTools::BamReader& R) const
{
	std::cout << "Debug protoSeeds\n";
	std::cout << "Complete: " << isComplete() << "\n";
	std::cout << "read1_havePrimary: " << read1_havePrimary << "\n";
	std::cout << "read2_havePrimary: " << read2_havePrimary << "\n";
	std::cout << "\tRead 1, alignments: " << read1_alignments.size() <<  "\n";
	for(unsigned int i = 0; i < read1_alignments.size(); i++)
	{
		std::cout << "\t\t" << i << "\n";
		//std::get<2>(read1_alignments.at(i)).BuildCharData();
		// std::cout << "\t\t\tName: " << std::get<2>(read1_alignments.at(i)).Name << "\n";
		std::cout << "\t\t\tRefID: " << R.GetReferenceData().at(std::get<2>(read1_alignments.at(i)).RefID).RefName << " [" << std::get<2>(read1_alignments.at(i)).RefID << "]" << "\n";
		std::cout << "\t\t\tIsPrimaryAlignment: " << std::get<2>(read1_alignments.at(i)).IsPrimaryAlignment() << "\n";
		std::cout << "\t\t\tIsPaired: " << std::get<2>(read1_alignments.at(i)).IsPaired() << "\n";
		std::cout << "\t\t\tIsProperPair: " << std::get<2>(read1_alignments.at(i)).IsProperPair() << "\n";
		std::cout << "\t\t\tIsMateMapped: " << std::get<2>(read1_alignments.at(i)).IsMateMapped() << "\n";
		if(std::get<2>(read1_alignments.at(i)).IsPaired()  && std::get<2>(read1_alignments.at(i)).IsMateMapped())
		{
			std::cout << "\t\t\tMateRefID: " << R.GetReferenceData().at(std::get<2>(read1_alignments.at(i)).MateRefID).RefName << " [" << std::get<2>(read1_alignments.at(i)).MateRefID << "]" << "\n";

		}
		std::cout << std::flush;
	}
	std::cout << "\tRead 2, alignments: " << read2_alignments.size() <<  "\n";
	for(unsigned int i = 0; i < read2_alignments.size(); i++)
	{
		std::cout << "\t\t" << i << "\n";
		//std::get<2>(read2_alignments.at(i)).BuildCharData();
		std::cout << "\t\t\tName: " << std::get<2>(read2_alignments.at(i)).Name << "\n";
		std::cout << "\t\t\tRefID: " << R.GetReferenceData().at(std::get<2>(read2_alignments.at(i)).RefID).RefName << " [" << std::get<2>(read2_alignments.at(i)).RefID << "]" << "\n";
		std::cout << "\t\t\tIsPrimaryAlignment: " << std::get<2>(read2_alignments.at(i)).IsPrimaryAlignment() << "\n";
		std::cout << "\t\t\tIsPaired: " << std::get<2>(read2_alignments.at(i)).IsPaired() << "\n";
		std::cout << "\t\t\tIsProperPair: " << std::get<2>(read2_alignments.at(i)).IsProperPair() << "\n";
		std::cout << "\t\t\tIsMateMapped: " << std::get<2>(read2_alignments.at(i)).IsMateMapped() << "\n";
		if(std::get<2>(read2_alignments.at(i)).IsPaired() && std::get<2>(read2_alignments.at(i)).IsMateMapped())
		{
			std::cout << "\t\t\tMateRefID: " << R.GetReferenceData().at(std::get<2>(read2_alignments.at(i)).MateRefID).RefName << " [" << std::get<2>(read2_alignments.at(i)).MateRefID << "]" << "\n";
		}
		std::cout << std::flush;
	}
}

bool protoSeeds::isComplete() const
{
	return (read1_havePrimary && read2_havePrimary);
}

protoSeeds::~protoSeeds() {
	// TODO Auto-generated destructor stub
}

}

} /* namespace mapper */
