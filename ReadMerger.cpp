#include "Header.h"
#include <unordered_map>
#include <string>

using namespace std;

bool ReadMerger(const string& read1, const string& qual1, const string& read2, const string& qual2, pair<string, string>& mergedRead) //merge paired reads using gap-less alignment (reads with phasing errors will be discarded)
{
	/*									Method
	R1 ---->	R1 ---->		B1 R1 ----> B2 R1 ---->  B3 R1 ---->   B4 R1 ---->
	R2   <----	R2   ----> (RC) B1 R2 ----> B2 R2  ----> B3 R2   ----> B4 R2    ----> etc
	*/

	//parameters
	unsigned minscore = 15, mismatch_penalty = 4, match_award = 1, phredOffset = 33, Qcap = 40; //use positive values
	double MisMatchDenominator = 20; //overlap length / MisMatchDenominatorless; than 5% mismatches

	unsigned readpos = 0, read1len = read1.length(), read2len = read2.length(), n, bestpos, mismatches;
	double bestscore = 0, secondBestScore = 0;
	int score, Q1, Q2;
	string read2Rev = ReverseComplement(read2);

	//match base by base reads and score
	while (readpos < read1len) { //iterate over read1
		score = 0;
		mismatches = 0;

		//Fix R1 in place, start R1 base 1 at R2 base 1, move R2 left to right one base at a time and check for matches/mismatches against R1

		for (n = 0; n + readpos < read1len && n < read2len; n++) { //stop loop when read2 (readpos) extends beyond the length of the read1 OR get to the end of R2

			if (read1[n + readpos] == read2Rev[n]) {
				score += match_award;
			} else {
				score -= mismatch_penalty;
				mismatches++;
			}

			if (mismatches > (read1len - readpos) / MisMatchDenominator) {  //length of potential overlap over maxmismatchdenominator
				break; //stop checking if read exceeds maximum mismatches for the whole overlap; improves preformance and accuracy
			}

		}

		if (score > bestscore) {
			secondBestScore = bestscore;
			bestscore = score;
			bestpos = readpos;
		}

		readpos++; //start read on next base
	}

	//check best alignment & merge
	if (bestscore > minscore && secondBestScore / bestscore < 0.9 && read2len + bestpos >= read1len) { //score is adequate, score is sufficently higher than the second best & R1 does not have adapter

		//attach start of read1
		mergedRead.first = read1.substr(0, bestpos);
		mergedRead.second = qual1.substr(0, bestpos);

		//take consensus across overlap
		for (n = 0; n + bestpos < read1len; n++) {

			Q1 = qual1[n + bestpos] - phredOffset;
			Q2 = qual2[n] - phredOffset;

			if (read1[n + bestpos] == read2Rev[n]) { //base is the same; match
				mergedRead.first += read1[n + bestpos];

				if (Q1 + Q2 > Qcap) { //add quality scores together; cap at 40
					mergedRead.second += toascii(Qcap + phredOffset);
				} else {
					mergedRead.second += toascii(Q1 + Q2 + phredOffset);
				}

			} else { //bases not the same; mismatch

				/*CLC BIO:If the two scores of the input reads are approximately equal, the resulting score will be very low which will reflect the fact that it is a very unreliable base.
				On the other hand, if one score is very low and the other is high, it is likely that the base with the high quality score is indeed correct,
				and this will be reflected in a relatively high quality score.*/

				if (Q1 >= Q2){ //R1 is more likely to be correct even if Qscores are the same
					mergedRead.first += read1[n + bestpos]; //use highest scoring base
					mergedRead.second += toascii((Q1 - Q2) + phredOffset); //
				} else {
					mergedRead.first += read2Rev[n]; //use highest scoring base
					mergedRead.second += toascii((Q2 - Q1) + phredOffset);
				}

			}

		}

		//attach end of read2
		mergedRead.first += read2Rev.substr(read1len - bestpos, std::string::npos);
		mergedRead.second += qual2.substr(read1len - bestpos, std::string::npos);

		return true;

	} else {
		return false; //no merge occured
	}

}