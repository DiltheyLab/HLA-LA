/*
 * pathFinder.h
 *
 *  Created on: Dec 20, 2016
 *      Author: diltheyat
 */

#ifndef PATHFINDER_H_
#define PATHFINDER_H_

#include <string>
#include "Utilities.h"
#include <exception>
#include <stdexcept>

class pathFinder {
public:
	pathFinder(const std::map<std::string, std::string>& arguments) : arguments(arguments)
	{

	}

	virtual ~pathFinder();

	std::string find_bwa() const
	{
		std::string forReturn;
		forReturn = arguments.count("bwa_bin") ? arguments.at("bwa_bin") : Utilities::findFileFromAlternatives({"C:/Users/AlexanderDilthey/bwa/bwa.exe", "/home/dilthey/bwa/bwa/bwa", "/apps/well/bwa/0.7.12/bwa"});

		if(! Utilities::fileExists(forReturn))
		{
			throw std::runtime_error("File " + forReturn + " not existing.");
		}
		assert(Utilities::fileExists(forReturn));
		return forReturn;
		// Utilities::findFileFromAlternatives({"C:/Users/AlexanderDilthey/bwa/bwa.exe", "/home/dilthey/bwa/bwa/bwa", "/apps/well/bwa/0.7.12/bwa"});
	}

	std::string find_bowtie() const
	{
		return Utilities::findFileFromAlternatives({"C:/Users/AlexanderDilthey/bowtie2-2.2.8", "/home/dilthey/bowtie2-2.2.8"});
	}

	std::string find_samtools() const
	{
		std::string forReturn;
		forReturn = arguments.count("samtools_bin") ? arguments.at("samtools_bin") : Utilities::findFileFromAlternatives({"C:/Users/AlexanderDilthey/samtools/samtools-0.1.19/samtools.exe", "/home/dilthey/samtools-0.1.18/samtools"});

		if(! Utilities::fileExists(forReturn))
		{
			throw std::runtime_error("File " + forReturn + " not existing.");
		}
		assert(Utilities::fileExists(forReturn));
		return forReturn;
		
		return forReturn;
	}
	
protected:
	const std::map<std::string, std::string>& arguments;
};

#endif /* PATHFINDER_H_ */
