/*
 * BWAmapper.cpp
 *
 *  Created on: 28.10.2015
 *      Author: AlexanderDilthey
 */

#include "BWAmapper.h"

#include <assert.h>
#include "../../Utilities.h"

namespace mapper {
namespace bwa {

BWAmapper::BWAmapper(const pathFinder& pF,  unsigned int threads_) {
	bwa_bin = pF.find_bwa();
	samtools_bin = pF.find_samtools();
	threads = threads_;
}

BWAmapper::~BWAmapper() {
	// TODO Auto-generated destructor stub
}

void BWAmapper::index(std::string referenceGenomeFastaFile)
{
	assert(Utilities::fileExists(bwa_bin));
	assert(Utilities::fileExists(referenceGenomeFastaFile));

	std::string bwa_index_cmd = bwa_bin + " index " + referenceGenomeFastaFile;
	int retCode = std::system(bwa_index_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bwa_index_cmd + " returned code "+Utilities::ItoStr(retCode));
	}
}

void BWAmapper::make_sure_ref_is_indexed(std::string referenceGenomeFastaFile)
{
	if(! ref_is_indexed(referenceGenomeFastaFile))
	{
		index(referenceGenomeFastaFile);
	}
	bool isNowIndexed = ref_is_indexed(referenceGenomeFastaFile);
	if(! isNowIndexed)
	{
		std::cerr << "Reference " << referenceGenomeFastaFile << " is not indexe\n" << std::flush;
	}
	assert(isNowIndexed);
}

std::string BWAmapper::getPRGonlyReferenceGenomePath_forGraph(std::string PRG_graph_dir)
{
	return PRG_graph_dir + "/mapping_PRGonly/referenceGenome.fa";

}

std::string BWAmapper::getReferenceGenomeForMapping_forGraph(std::string PRG_graph_dir, bool mapAgainstCompleteGenome)
{
	std::string PRGonlyReferenceGenomePath = getPRGonlyReferenceGenomePath_forGraph(PRG_graph_dir);
	std::string extendedReferenceGenomePath;
	if(Utilities::fileExists(PRG_graph_dir + "/extendedReferenceGenomePath.txt"))
	{
		extendedReferenceGenomePath  = Utilities::getFirstLine(PRG_graph_dir + "/extendedReferenceGenomePath.txt");
	}
	else
	{
		extendedReferenceGenomePath  = PRG_graph_dir + "/extendedReferenceGenome/extendedReferenceGenome.fa";
		assert(Utilities::fileExists(extendedReferenceGenomePath));
	}

	std::string referenceGenomeForMapping = mapAgainstCompleteGenome ? extendedReferenceGenomePath : PRGonlyReferenceGenomePath;
	return referenceGenomeForMapping;
}

void BWAmapper::createRemappedBAM_forGraph(std::string PRG_graph_dir, std::string FASTQ1, std::string FASTQ2, std::string FASTQU, bool mapAgainstCompleteGenome, std::string longReads, std::string outputBAM, bool remap_with_a)
{
	assert((longReads == "") || (longReads == "ont2d") || (longReads == "pacbio"));

	if(longReads.length())
	{
		assert(FASTQU.length());
		assert(FASTQ1.length() == 0);
		assert(FASTQ2.length() == 0);
	}
	else
	{
		assert(FASTQU.length() == 0);
	}

	std::string referenceGenomeForMapping = getReferenceGenomeForMapping_forGraph(PRG_graph_dir, mapAgainstCompleteGenome);
	if(longReads.length() != 0)
	{
		mapLong(referenceGenomeForMapping, FASTQU, outputBAM, remap_with_a, longReads);
	}
	else
	{
		map(referenceGenomeForMapping, FASTQ1, FASTQ2, outputBAM, remap_with_a);
	}

	std::cout << Utilities::timestamp() << "Remapping done.\n" << std::flush;
	assert(Utilities::fileExists(outputBAM));
	assert(Utilities::fileExists(outputBAM+".bai"));
}

bool BWAmapper::ref_is_indexed(std::string referenceGenomeFastaFile)
{
	std::vector<std::string> suffixes = {".sa", ".ann", ".bwt"};
	for(unsigned int i = 0; i < suffixes.size(); i++)
	{
		if(! Utilities::fileExists(referenceGenomeFastaFile + suffixes.at(i)))
		{
			return false;
		}
	}

	return true;
}

void BWAmapper::mapUnpaired(std::string indexedReferenceGenome, std::string FASTQ, std::string outputBAM, bool withA)
{
	assert(Utilities::fileExists(bwa_bin));
	assert(Utilities::fileExists(samtools_bin));

	assert(Utilities::fileExists(FASTQ));

	make_sure_ref_is_indexed(indexedReferenceGenome);

	if(Utilities::fileExists(outputBAM))
	{
		Utilities::deleteFile(outputBAM);
	}

	std::string outputUnsorted = outputBAM + ".unsorted";
	if(Utilities::fileExists(outputUnsorted))
	{
		Utilities::deleteFile(outputUnsorted);
	}

	//std::string with_a = (withA) ? "-a -L30 -k15 " : "";
	std::string with_a = (withA) ? "-a " : "";

	//std::string bwa_mapping_cmd = bwa_bin + " mem " + indexedReferenceGenome + " " + FASTQ1 + " " + FASTQ2 + " > " + outputSAM;
	std::string bwa_mapping_cmd = bwa_bin + " mem -M -k 15 " + with_a + indexedReferenceGenome + " " + FASTQ  + " | " + samtools_bin + " view -Sb - > " + outputUnsorted;
	std::cout << "BWA mapping command: " << bwa_mapping_cmd << "\n" << std::flush;
	std::cerr << bwa_mapping_cmd << "\n" << std::flush;
	int retCode = std::system(bwa_mapping_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bwa_mapping_cmd + " returned code "+Utilities::ItoStr(retCode));
	}

	assert(outputBAM.substr(outputBAM.length() - 4) == ".bam");
	std::string sortedBAM_noNAM = outputBAM.substr(0, outputBAM.length() - 4);

	std::string samtools_sort_cmd = samtools_bin + " sort -o " + outputBAM + " " + outputUnsorted;
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

void BWAmapper::mapLong(std::string indexedReferenceGenome, std::string FASTQ, std::string outputBAM, bool withA, std::string longMode)
{
	assert((longMode == "pacbio") || (longMode == "ont2d"));

	assert(Utilities::fileExists(bwa_bin));
	assert(Utilities::fileExists(samtools_bin));

	assert(Utilities::fileExists(FASTQ));

	make_sure_ref_is_indexed(indexedReferenceGenome);

	if(Utilities::fileExists(outputBAM))
	{
		Utilities::deleteFile(outputBAM);
	}

	std::string outputUnsorted = outputBAM + ".unsorted";
	if(Utilities::fileExists(outputUnsorted))
	{
		Utilities::deleteFile(outputUnsorted);
	}

	//std::string with_a = (withA) ? "-a -L30 -k15 " : "";
	std::string with_a = (withA) ? "-a " : "";
	
	int threads_minus_1 = threads - 1;
	assert(threads_minus_1 >= 0);
	
	//std::string bwa_mapping_cmd = bwa_bin + " mem " + indexedReferenceGenome + " " + FASTQ1 + " " + FASTQ2 + " > " + outputSAM;
	std::string bwa_mapping_cmd = bwa_bin + " mem -t" + Utilities::ItoStr(threads) + " -x " + longMode + " -M " + with_a + indexedReferenceGenome + " " + FASTQ + " | " + samtools_bin + " view -@ " + Utilities::ItoStr(threads_minus_1) + " -Sb - > " + outputUnsorted;
	std::cerr << bwa_mapping_cmd << "\n" << std::flush;
	int retCode = std::system(bwa_mapping_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bwa_mapping_cmd + " returned code "+Utilities::ItoStr(retCode));
	}

	assert(outputBAM.substr(outputBAM.length() - 4) == ".bam");
	std::string sortedBAM_noNAM = outputBAM.substr(0, outputBAM.length() - 4);

	std::string samtools_sort_cmd = samtools_bin + " sort -@ " + Utilities::ItoStr(threads) + " -o " + outputBAM + " " + outputUnsorted;
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

void BWAmapper::map(std::string indexedReferenceGenome, std::string FASTQ1, std::string FASTQ2, std::string outputBAM, bool withA)
{
	assert(Utilities::fileExists(bwa_bin));
	assert(Utilities::fileExists(samtools_bin));

	assert(Utilities::fileExists(FASTQ1));
	assert(Utilities::fileExists(FASTQ2));

	make_sure_ref_is_indexed(indexedReferenceGenome);

	if(Utilities::fileExists(outputBAM))
	{
		Utilities::deleteFile(outputBAM);
	}

	std::string outputUnsorted = outputBAM + ".unsorted";
	if(Utilities::fileExists(outputUnsorted))
	{
		Utilities::deleteFile(outputUnsorted);
	}

	//std::string with_a = (withA) ? "-a -L30 -k15 " : ""; 
	std::string with_a = (withA) ? "-a " : ""; 

	//std::string bwa_mapping_cmd = bwa_bin + " mem " + indexedReferenceGenome + " " + FASTQ1 + " " + FASTQ2 + " > " + outputSAM;
	int threads_minus_1 = threads - 1;
	assert(threads_minus_1 >= 0);
	
	std::string bwa_mapping_cmd = bwa_bin + " mem -t" + Utilities::ItoStr(threads) + " -M " + with_a + indexedReferenceGenome + " " + FASTQ1 + " " + FASTQ2 + " | " + samtools_bin + " view -@ " + Utilities::ItoStr(threads_minus_1) + " -Sb - > " + outputUnsorted;
	std::cerr << bwa_mapping_cmd << "\n" << std::flush;
	int retCode = std::system(bwa_mapping_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bwa_mapping_cmd + " returned code "+Utilities::ItoStr(retCode));
	}

	assert(outputBAM.substr(outputBAM.length() - 4) == ".bam");
	std::string sortedBAM_noNAM = outputBAM.substr(0, outputBAM.length() - 4);

	std::string samtools_sort_cmd = samtools_bin + " sort -@ " + Utilities::ItoStr(threads) + " -o " + outputBAM + " " + outputUnsorted;
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

void BWAmapper::map_all_unpaired_unsorted(std::string indexedReferenceGenome, std::string FASTQ, std::string outputBAM)
{
	assert(Utilities::fileExists(bwa_bin));
	assert(Utilities::fileExists(samtools_bin));

	assert(Utilities::fileExists(FASTQ));

	make_sure_ref_is_indexed(indexedReferenceGenome);

	if(Utilities::fileExists(outputBAM))
	{
		Utilities::deleteFile(outputBAM);
	}

	std::string outputUnsorted = outputBAM + ".unsorted";
	if(Utilities::fileExists(outputUnsorted))
	{
		Utilities::deleteFile(outputUnsorted);
	}


	std::string bwa_mapping_cmd = bwa_bin + " mem -t" + Utilities::ItoStr(threads) + " -M -a " + indexedReferenceGenome + " " + FASTQ + " | " + samtools_bin + " view -Sb - > " + outputBAM;
	std::cerr << bwa_mapping_cmd << "\n" << std::flush;
	int retCode = std::system(bwa_mapping_cmd.c_str());
	if(retCode != 0)
	{
		throw std::runtime_error("Command " + bwa_mapping_cmd + " returned code "+Utilities::ItoStr(retCode));
	}

	assert(outputBAM.substr(outputBAM.length() - 4) == ".bam");
	std::string sortedBAM_noNAM = outputBAM.substr(0, outputBAM.length() - 4);

//	std::string samtools_sort_cmd = samtools_bin + " sort " + outputUnsorted + " " + sortedBAM_noNAM;
//	std::cerr << samtools_sort_cmd << "\n" << std::flush;
//	retCode = std::system(samtools_sort_cmd.c_str());
//	if(retCode != 0)
//	{
//		throw std::runtime_error("Command " + samtools_sort_cmd + " returned code "+Utilities::ItoStr(retCode));
//	}
//	assert(Utilities::fileExists(outputBAM));
//
//	std::string samtools_index_cmd = samtools_bin + " index " + outputBAM;
//	std::cerr << samtools_index_cmd << "\n" << std::flush;
//	retCode = std::system(samtools_index_cmd.c_str());
//	if(retCode != 0)
//	{
//		throw std::runtime_error("Command " + samtools_index_cmd + " returned code "+Utilities::ItoStr(retCode));
//	}
//	assert(Utilities::fileExists(outputBAM + ".bai"));

//	Utilities::deleteFile(outputUnsorted);
}

} /* namespace bwa */
} /* namespace mapper */
