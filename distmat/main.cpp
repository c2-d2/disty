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


static const uint8_t nt256_nt4[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};


KSEQ_INIT(gzFile, gzread)

using namespace std;

void usage() {
	cout << "Program: distmat - compute distance matrix from a FASTA alignment file" << endl;
    cout << "Usage:   distmat <inp.aln>" << endl << endl;
    cout << "Options:" << endl;
    cout << "  -s   skip columns with N's" << endl;
    cout << "  -n   don't normalize (to account for skipped columns)" << endl;
	return;
}

/*
 * Compare sequences
 */
template<typename T>
void compare_seqs(const string &s1, const string &s2, T (&counts)[5][5]) {
    // const string &skipped,
	auto l1=s1.size();
    auto l2=s2.size();
	assert (l1==l2);

    for (int i=0;i<5;i++){
        for (int j=0;j<5;j++){
            counts[i][j]=0;
        }
    }

	for (int i=0;i<l1;i++){
        ++counts[nt256_nt4[s1[i]]][nt256_nt4[s2[i]]];
	}
}

template <typename T, typename U>
void dist_mat(const T &seqs, U &matrix){
    int counts[5][5];
	auto no_seqs=seqs.size();
	//int len=seqs[0].size();
	for (int i=0;i<no_seqs;i++){
		for (int j=0;j<=i;j++){
			compare_seqs(seqs[i], seqs[j], counts);
            matrix[i][j]=matrix[j][i]=counts[0][0]+counts[1][1]+counts[2][2]+counts[3][3]+counts[4][4];
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
		cerr << "Usage: " << argv[0] << " <in.seq>" << endl;
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


	cerr << "Constructing empty matrices" << std::endl;

	vector<vector<int>> distance_matrix(count, vector<int>(count, 0));

	cerr << "Computing distance matrices" << std::endl;

	dist_mat(seqs, distance_matrix);

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
            cout << "\t" << distance_matrix[i][j];

            //cout << "\t" << (int) round( (1.0*non_ns_diff[i][j]/non_ns[i][j])*len );
		}
		cout << endl;
	}

	return 0;
}
