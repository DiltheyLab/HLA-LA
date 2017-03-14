/*
 * Utilities.h
 *
 *  Created on: 23.05.2011
 *      Author: Alexander Dilthey
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <boost/random.hpp>
#include <algorithm>
#include <functional>

using namespace std;

class Utilities {
public:
	Utilities();
	virtual ~Utilities();

	static map<string, string> readFASTA(std::string file, bool fullIdentifier = false);
	static void writeFASTA(std::string file, const std::map<std::string, std::string>& sequences);

	static vector<string> split(const string &s, char delim, vector<string> &elems);
	static vector<string> split(const string &s, char delim);
	static string ItoStr(int i);
	static std::string BtoStr(bool b);
	static std::vector<std::string> BtoStr(std::vector<bool> b);

	static string DtoStr(double d);
	static string PtoStr(void* p);

	static int StrtoI(string s);
	static std::vector<int> StrtoI(std::vector<std::string> s);

	static vector<string> ItoStr(vector<int> i);
	static double StrtoD(string s);
	static bool StrtoB(string s);
	static long long StrtoLongLong(string s);

	static string join(vector<string> parts, string delim);
	static void eraseNL(string& s);
	static int chooseFromVector(vector<double>& v);

	static vector<string> split(string input, string delimiter);


	static std::string timestamp();

	static std::pair<double, unsigned int> findVectorMax(const std::vector<double>& v);
	static std::pair<double, int> findIntMapMax(std::map<int, double>& m);
	static std::pair<double, std::string> findStringMapMax(const std::map<std::string, double>& m);
	static std::pair<double, int> findIntMapMaxP_nonCritical(std::map<int, double>& m, unsigned int* thisSeed);
	static std::pair<double, unsigned int> findVectorMaxP(std::vector<double>& v);
	static std::pair<double, unsigned int> findVectorMaxP_nonCritical(std::vector<double>& v, unsigned int* thisSeed);

	static std::string generateRandomSequence(int length);
	static std::string generateRandomSequenceWithGaps(int length, double gapFrequency = 0.1);
	static std::pair<std::string,std::vector<int>> modifySequence(std::string sequence, std::vector<int> positionOrigin, double mutationFrequence = 0.1, double insertionFrequence = 0.05, int insertionMaxLength = 4, double deletionFrequence = 0.1, int deletionMaxLength = 5);

	static char randomNucleotide();
	static char randomNucleotide(boost::mt19937& rng);


	static int randomNumber(int max);
	static int randomNumber_nonCritical(int max, unsigned int* thisSeed);

	static double randomDouble();

	static std::string repeatString(std::string s, int repeatNumber);

	static bool extractBit(unsigned int number, unsigned int bit);

	static int readStatus(std::string statusFile);
	static void writeStatus(std::string statusFile, int status);

	static std::string JoinMapUInt2Str(std::map<std::string, unsigned int> M);

	static void check_map_is_normalized(std::map<char, double> m);
	static std::map<char, double> normalize_map(std::map<char, double> m);
	static int choose_from_nonnormalized_map(std::map<int, int> m);
	static char choose_from_normalized_map(std::map<char, double> m);
	static char choose_from_normalized_map(std::map<char, double> m, boost::mt19937& rng);
	static std::string choose_uniformly_from_vector(std::vector<std::string> v);

	static void normalize_vector(std::vector<double>& v);
	static std::vector<double> normalize_log_vector(const std::vector<double>& v);

	static std::string seq_reverse_complement(std::string sequence);
	static char reverse_char_nucleotide(char c);

	static bool oneBernoulliTrial(double p);
	static bool oneBernoulliTrial(double p, boost::mt19937& rng);

	static double PhredToPCorrect(unsigned char nucleotideQuality);
	static unsigned char PCorrectToPhred(double PCorrect);


	static bool fileReadable(std::string file);

	static bool directoryExists(std::string dir);
	static void makeDir(std::string dir);
	static bool fileExists(std::string filepath);
	static size_t fileLastWrite(std::string filepath);
	static void clearDirectory(std::string dir);
	static void deleteFile(std::string file);
	static std::vector<std::string> filesInDirectory(std::string path);

	static void make_or_clearDirectory(std::string dir);

	static std::string removeFROM(std::string readID);
		
	static double LogSumLogPs (const std::vector<double>& v);
	static double logAvg(double a, double b);

	static std::string getFirstLine(std::string file);
	static std::vector<std::string> getAllLines(std::string file);
	static void printToFile(std::string file, std::string what);

	static std::string removeGaps(std::string in);

	static int randomPoisson(double mean);
	
	static std::string findFileFromAlternatives(std::vector<std::string> alternatives);

	static bool intervalsOverlap(int x1, int x2, int y1, int y2);

	static double proportionN(std::string S);

	static std::vector<std::string> partitionStringIntokMers(std::string str, int k);

	static std::string getCWD();
	static void setCWD(std::string path);

	static bool sequence_all_Ns(const std::string S);
	static std::vector<std::string> get_map_keys_sorted_by_value(const std::map<std::string, int>& m);
};

template <typename Collection,typename Predicate>
Collection filterNot(Collection col,Predicate predicate ) {
    auto returnIterator = std::remove_if(col.begin(),col.end(),predicate);
    col.erase(returnIterator,std::end(col));
    return col;
}

template <typename Collection,typename Predicate>
Collection filter(Collection col,Predicate predicate) {
 //capture the predicate in order to be used inside function
 auto fnCol = filterNot(col,[predicate](typename Collection::value_type i) { return !predicate(i);});
 return fnCol;
}

template <typename Collection,typename unop>
  Collection mymap(Collection col,unop op) {
  std::transform(col.begin(),col.end(),col.begin(),op);
  return col;
}

template <typename Collection,typename Collection2, typename unop>
  Collection2 mymap2(Collection col,unop op) {
  std::transform(col.begin(),col.end(),col.begin(),op);
  return col;
}



extern unsigned int globalRandRSeed;
#endif /* UTILITIES_H_ */


