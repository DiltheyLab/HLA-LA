/*
 * Bowtie2mapper.h
 *
 *  Created on: 10.03.2016
 *      Author: AlexanderDilthey
 */

#ifndef MAPPER_BOWTIE2_BOWTIE2MAPPER_H_
#define MAPPER_BOWTIE2_BOWTIE2MAPPER_H_

#include <string>
#include "../../pathFinder.h"


namespace mapper {
namespace bowtie2 {


class Bowtie2mapper {
protected:
	std::string bowtie2_dir;
	std::string samtools_bin;

public:
	Bowtie2mapper(const pathFinder& pF);


	void createIndex(std::string referenceGenomeFastaFile, std::string indexName);
	void mapUnpaired(std::string indexName, std::string FASTQ, std::string outputBAM);
	void make_sure_ref_is_indexed(std::string referenceGenomeFastaFile, std::string indexName);
	bool ref_is_indexed(std::string indexName);
};

} /* namespace bowtie2 */
} /* namespace mapper */

#endif /* MAPPER_BOWTIE2_BOWTIE2MAPPER_H_ */
