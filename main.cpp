#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <cmath>
#include <zlib.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

void usage() {
	cout << "Program: distmat - compute distance matrix from a FASTA alignment file" << endl;
	cout << "Usage:   distmat <inp.aln>" << endl;
	return;
}

/*
 * Compute number of differences between two sequences (N's ignored) and the number of non-N positions.
 */
void compare_seqs(const string &s1, const string &s2, int &non_ns, int &non_ns_diff) {
	assert (s1.size()==s2.size());

	non_ns=0;
	non_ns_diff=0;

	for (int i=0;i<s1.size();i++){
		if (s1[i]!='N' && s2[i]!='N') {
			non_ns++;

			if (s1[i]!=s2[i]) {
				non_ns_diff++;
			}
		}
	}
}

template <typename T, typename U>
void dist_mat(const T &seqs, U &non_ns, U &non_ns_diff){
	for (int i=0;i<seqs.size();i++){
		for (int j=0;j<=i;j++){
			compare_seqs(seqs[i], seqs[j], non_ns[i][j], non_ns_diff[i][j]);
			non_ns[j][i]=non_ns[i][j];
			non_ns_diff[j][i]=non_ns_diff[i][j];
		}
	}
}

int main (int argc, const char **argv) {
	if (argc!=2){
		usage();
		return EXIT_FAILURE;
	}

	vector<string> names, seqs;

	gzFile fp;
	kseq_t *seq;
	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);

	int count=0;
	while ((l = kseq_read(seq)) >= 0) {
		names.push_back(seq->name.s);
		string s(seq->seq.s);
		for (auto & c: s) {
			c = toupper(c);
		}
		seqs.push_back(s);
		count++;
	}
	kseq_destroy(seq);
	gzclose(fp);

	int len=seqs[0].size();

	cerr << "Constructing empty matrices" << std::endl;

	vector<vector<int>> non_ns(count, vector<int>(count, 0));
	vector<vector<int>> non_ns_diff(count, vector<int>(count, 0));

	cerr << "Computing distance matrices" << std::endl;

	dist_mat(seqs, non_ns, non_ns_diff);

	/*
	 * print the output
	 */

	cout << "#taxid";
	for (int j=0;j<count;j++){
		cout << "\t" << names[j];
	}
	cout << endl;

	for (int i=0;i<count;i++){
		cout << names[i];
		for (int j=0;j<count;j++){
			cout << "\t" << (int) round( (1.0*non_ns_diff[i][j]/non_ns[i][j])*len );
		}
		cout << endl;
	}

	return 0;
}
