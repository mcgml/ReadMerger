#include <string>
#include <iostream>
#include <unordered_map>

using namespace std;

string ReverseComplement(const string& dna);
bool ReadMerger(const string& read1, const string& qual1, const string& read2, const string& qual2, pair<string, string>& mergedRead); //merge paired reads using gap-less alignment (reads with phasing errors will be discarded)