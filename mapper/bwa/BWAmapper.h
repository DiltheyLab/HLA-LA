/*
 * BWAmapper.h
 *
 *  Created on: 28.10.2015
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_BWA_BWAMAPPER_H_
#define MAPPER_BWA_BWAMAPPER_H_

#include "../../pathFinder.h"

#include <string>

namespace mapper {
namespace bwa {

class BWAmapper {
protected:
	std::string bwa_bin;
	std::string samtools_bin;
	unsigned int threads;
public:
	BWAmapper(const pathFinder& pF, unsigned int threads_ = 1);
	virtual ~BWAmapper();
	void index(std::string referenceGenomeFastaFile);
	void map(std::string indexedReferenceGenome, std::string FASTQ1, std::string FASTQ2, std::string outputBAM, bool withA = false);
	void mapLong(std::string indexedReferenceGenome, std::string FASTQ, std::string outputBAM, bool withA, std::string longMode);

	void mapUnpaired(std::string indexedReferenceGenome, std::string FASTQ, std::string outputBAM, bool withA = false);
	void map_all_unpaired_unsorted(std::string indexedReferenceGenome, std::string FASTQ, std::string outputBAM);
	void make_sure_ref_is_indexed(std::string referenceGenomeFastaFile);
	bool ref_is_indexed(std::string referenceGenomeFastaFile);
};

} /* namespace bwa */
} /* namespace mapper */

#endif /* MAPPER_BWA_BWAMAPPER_H_ */
