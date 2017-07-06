#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>
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
 * Compute number of differences between two sequences (N's ignored).
 */
int dist(const string &s1, const string &s2) {
	assert (s1.size()==s2.size());
	int dist=0;
	for (int i=0;i<s1.size();i++){
		if (s1[i]!=s2[i] && s1[i]!='N' && s2[i] != 'N'){
			dist++;
		}
	}
	return dist;
}

/*
 * Compute divergence between two sequences (N's ommited from stats)
 */
float diver(const string &s1, const string &s2) {
	assert (s1.size()==s2.size());
	int poss=0;
	int diff=0;
	for (int i=0;i<s1.size();i++){
		if (s1[i]=='N' || s2[j]=='N') {
			continue;
		}

		poss++;

		if (s1[i]!=s2[i]){
			diff++;
		}
	}
	return 1.0*diff/poss;
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

	int d[count][count];

	for (int i=0;i<count;i++){
		for (int j=0;j<=i;j++){
			d[i][j]=d[j][i]=dist(seqs[i],seqs[j]);
		}
	}

	cout << "#taxid";
	for (int j=0;j<count;j++){
		cout << "\t" << names[j];
	}
	cout << endl;

	for (int i=0;i<count;i++){
		cout << names[i];
		for (int j=0;j<count;j++){
			cout << "\t" << d[i][j];
		}
		cout << endl;
	}

	return 0;
}
