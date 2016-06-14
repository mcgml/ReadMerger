#include <string>
#include <Header.h>

using namespace std;

string ReverseComplement(const string& dna) {

	string revcomp;

	for (short n = (dna.length() - 1); n > -1; n--) {

		if (dna[n] == 'A') {
			revcomp += 'T';
		} else if (dna[n] == 'T') {
			revcomp += 'A';
		} else if (dna[n] == 'G') {
			revcomp += 'C';
		} else if (dna[n] == 'C') {
			revcomp += 'G';
		} else if (dna[n] == 'a') {
			revcomp += 't';
		} else if (dna[n] == 't') {
			revcomp += 'a';
		} else if (dna[n] == 'g') {
			revcomp += 'c';
		} else if (dna[n] == 'c') {
			revcomp += 'g';
		} else {
			revcomp += dna[n];
		}

	}

	return revcomp;
}