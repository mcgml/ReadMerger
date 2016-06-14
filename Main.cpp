/*ReadMerger.cpp : Defines the entry point for the console application.
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Status: Release
*/

/*NB
1-FASTQs must be adapter trimmed
2-Parameters used for merging expect overlap --sensitivity may need adjusting
3-Non-overlapping reads are discarded
*/

#include <iostream>
#include <fstream>
#include <string>
#include <Header.h>

using namespace std;

const string version = "v0.1";

bool isHeaderMatched(const string& header1, const string& header2){

	unsigned mismatches = 0, n, header1Len = header1.length(), header2Len = header2.length();

	if (header1Len != header2Len){
		return false;
	}

	for (n = 0; n < header1Len; n++){
		if (header1[n] != header2[n]){
			mismatches++;
		}

		if (mismatches > 1){
			return false;
		}
	}

	if (mismatches == 1){
		return true;
	} else {
		return false;
	}

}

int main(int argc, char* argv[]) {

	//check argument number is correct; print usage
	if (argc != 3) { //program r1 r2
		std::cerr << "\nProgram: ReadMerger " << version << endl;
		std::cerr << "Contact: Matthew Lyon, WRGL/UoS (matthew.lyon@salisbury.nhs.uk)\n" << endl;
		std::cerr << "Usage: ReadMerger <Read1.fastq> <Read2.fastq>\n" << endl;
		return -1;
	}

	string rLine1, rLine2, header1, header2, seq1, seq2, qual1, qual2, filename = argv[1];
	pair<string, string> mergedRead;
	unsigned lNo = 0;
	unsigned long totalReads = 0, mergedReads = 0;

	ifstream R1_in(argv[1]);
	ifstream R2_in(argv[2]);
	ofstream out((filename.substr(0, filename.find_first_of('_')) + "_merged.fastq").c_str());
	
	//parse FASTQs
	if (R1_in.is_open() && R2_in.is_open()) {
		while (R1_in.good() && R2_in.good()) {
			getline(R1_in, rLine1);
			getline(R2_in, rLine2);
			
			//skip empty lines
			if (rLine1 == "" || rLine2 == "") {
				continue;
			}

			lNo++;

			if (lNo == 1) {
				header1 = rLine1;
				header2 = rLine2;

				totalReads++;
				
				if (isHeaderMatched(header1, header2) == false){
					cerr << "ERROR: Read headers not paired. Check file input" << endl;
					return -1;
				}

			} else if (lNo == 2) {
				seq1 = rLine1;
				seq2 = rLine2;
			} else if (lNo == 4) {
				qual1 = rLine1;
				qual2 = rLine2;

				lNo = 0;
				
				//merge reads
				if (ReadMerger(seq1, qual1, seq2, qual2, mergedRead) == 1) { //merged successful
					mergedReads++;
					out << header1 << "\012" << mergedRead.first << "\012+\012" << mergedRead.second << "\012";
				}

			}
		}

	} else {
		std::cerr << "ERROR: Unable to open FASTQ file(s)" << endl;
		return -1;
	}

	cout << "TotalReads\tMergedReads\tFraction" << endl;
	cout << totalReads << "\t" << mergedReads << "\t" << (float)mergedReads / totalReads << endl;

	return 0;
}