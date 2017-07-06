#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "kseq.h"

using namespace std;

void usage() {
	cout << "Program: distmat - compute distance matrix from a FASTA alignment file" << endl;
	cout << "Usage:   distmat <inp.aln> <out.pref>" << endl;
	return;
}

int main (int argc, const char **argv) {
	if (argc!=2){
		usage();
		return EXIT_FAILURE;
	}
}
