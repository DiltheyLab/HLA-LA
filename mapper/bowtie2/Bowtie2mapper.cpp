/*
 * Bowtie2mapper.cpp
 *
 *  Created on: 10.03.2016
 *      Author: AlexanderDilthey
 */

#include "Bowtie2mapper.h"

#include <assert.h>
#include "../../Utilities.h"

namespace mapper {
namespace bowtie2 {

Bowtie2mapper::Bowtie2mapper(const pathFinder& pF) {
	bowtie2_dir = pF.find_bowtie();
	samtools_bin = pF.find_samtools();

}

void Bowtie2mapper::mapUnpaired(std::string indexName, std::string FASTQ, std::string outputBAM)
{
	std::string bowtie2_bin = bowtie2_dir + "/bowtie2";
	assert(Utilities::fileExists(bowtie2_bin));

	assert(Utilities::fileExists(samtools_bin));
	assert(Utilities::fileExists(FASTQ));
	assert(ref_is_indexed(indexName));

	if(Utilities::fileExists(outputBAM))
	{
		Utilities::deleteFile(outputBAM);
	}

	std::string outputSAM = outputBAM + ".sam";
	if(Utilities::fileExists(outputSAM))
	{
		Utilities::deleteFile(outputSAM);
	}

	std::string outputUnsorted = outputBAM + ".unsorted";
	if(Utilities::fileExists(outputUnsorted))
	{
		Utilities::deleteFile(outputUnsorted);
	}

	//std::string bwa_mapping_cmd = bwa_bin + " mem " + indexedReferenceGenome + " " + FASTQ1 + " " + FASTQ2 + " > " + outputSAM;
	std::string bowtie_mapping_cmd = bowtie2_bin + " -a -L 25 --omit-sec-seq -x " + indexName + " -U " + FASTQ + " -S " + outputSAM;
	std::cerr << bowtie_mapping_cmd << "\n" << std::flush;
	int retCode = std::system(bowtie_mapping_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bowtie_mapping_cmd + " returned code "+Utilities::ItoStr(retCode));
	}
	assert(Utilities::fileExists(outputSAM));

	std::string bam_conversion_cmd = samtools_bin + " view -bS " + outputSAM + " > " + outputUnsorted;
	std::cerr << bam_conversion_cmd << "\n" << std::flush;
	retCode = std::system(bam_conversion_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bam_conversion_cmd + " returned code "+Utilities::ItoStr(retCode));
	}
	assert(Utilities::fileExists(outputUnsorted));

	assert(outputBAM.substr(outputBAM.length() - 4) == ".bam");
	std::string sortedBAM_noNAM = outputBAM.substr(0, outputBAM.length() - 4);

	std::string samtools_sort_cmd = samtools_bin + " sort " + outputUnsorted + " " + sortedBAM_noNAM;
	std::cerr << samtools_sort_cmd << "\n" << std::flush;
	retCode = std::system(samtools_sort_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + samtools_sort_cmd + " returned code "+Utilities::ItoStr(retCode));
	}
	assert(Utilities::fileExists(outputBAM));

	std::string samtools_index_cmd = samtools_bin + " index " + outputBAM;
	std::cerr << samtools_index_cmd << "\n" << std::flush;
	retCode = std::system(samtools_index_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + samtools_index_cmd + " returned code "+Utilities::ItoStr(retCode));
	}
	assert(Utilities::fileExists(outputBAM + ".bai"));

	Utilities::deleteFile(outputUnsorted);
}


void Bowtie2mapper::createIndex(std::string referenceGenomeFastaFile, std::string indexName)
{
	std::string fn_pos = indexName + ".1.bt2";
	std::string fn_rev = indexName + ".rev.1.bt2";

	if(Utilities::fileExists(fn_pos))
		Utilities::deleteFile(fn_pos);

	if(Utilities::fileExists(fn_rev))
			Utilities::deleteFile(fn_rev);

	//std::string dir_before = Utilities::getCWD();
	//Utilities::setCWD(bowtie2_dir);
	//std::string indexname_absolute = dir_before + "/" + indexName;

	std::string bowtie2_build_bin = bowtie2_dir + "/bowtie2-build";

	if(bowtie2_dir.find("C:/") != std::string::npos)
	{
		bowtie2_build_bin += "-s-debug";
	}

	assert(Utilities::fileExists(bowtie2_build_bin));
	assert(Utilities::fileExists(referenceGenomeFastaFile));

	std::string bowtie2_index_cmd = "" + bowtie2_build_bin + " --offrate 0 --ftabchars 14 " + referenceGenomeFastaFile + " " + indexName;

	std::cerr << bowtie2_index_cmd << "\n" << std::flush;

	int retCode = std::system(bowtie2_index_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bowtie2_index_cmd + " returned code "+Utilities::ItoStr(retCode));
	}

	// Utilities::setCWD(dir_before);


	std::cout << "Exists " << fn_pos << ": " << Utilities::fileExists(fn_pos) << "\n";
	std::cout << "Exists " << fn_rev << ": " << Utilities::fileExists(fn_rev) << "\n" << std::flush;;

	assert(ref_is_indexed(indexName));

	/*
	 *

	std::string dir_before = Utilities::getCWD();
	assert(Utilities::fileExists(referenceGenomeFastaFile));

	Utilities::setCWD(bowtie2_dir);

	std::string indexname_absolute = dir_before + "/" + indexName;
	std::string referenceGenomeFastaFile_absolute = dir_before + "/" + referenceGenomeFastaFile;

	std::string bowtie2_build_bin = "./bowtie2-build";

	assert(Utilities::fileExists(bowtie2_build_bin));

	std::string bowtie2_index_cmd = bowtie2_build_bin + " " + referenceGenomeFastaFile_absolute + " " + indexname_absolute;

	std::cerr << bowtie2_index_cmd << "\n" << std::flush;

	int retCode = std::system(bowtie2_index_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bowtie2_index_cmd + " returned code "+Utilities::ItoStr(retCode));
	}

	Utilities::setCWD(dir_before);
	*/


}

bool Bowtie2mapper::ref_is_indexed(std::string indexName)
{
	std::string fn_pos = indexName + ".1.bt2";
	std::string fn_rev = indexName + ".rev.1.bt2";

	return (Utilities::fileExists(fn_pos) && Utilities::fileExists(fn_rev));
}

void Bowtie2mapper::make_sure_ref_is_indexed(std::string referenceGenomeFastaFile, std::string indexName)
{
	if(! ref_is_indexed(indexName))
	{
		createIndex(referenceGenomeFastaFile, indexName);
	}
	assert(ref_is_indexed(indexName));
}


} /* namespace bowtie2 */
} /* namespace mapper */
