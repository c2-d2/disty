/*
	The MIT License

	Copyright (c) 2017 Karel Brinda <kbrinda@hsph.harvard.edu>

	Permission is hereby granted, free of charge, to any person obtaining
	a copy of this software and associated documentation files (the
	"Software"), to deal in the Software without restriction, including
	without limitation the rights to use, copy, modify, merge, publish,
	distribute, sublicense, and/or sell copies of the Software, and to
	permit persons to whom the Software is furnished to do so, subject to
	the following conditions:

	The above copyright notice and this permission notice shall be
	included in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
	BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
	ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
	CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
 */



#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <cmath>
#include <zlib.h>
#include <getopt.h>

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

/*
 * Print usage.
 */
void usage() {
    cerr<<
    "\n"
    "Program: distmat - compute distance matrix from a FASTA alignment file\n"
    "\n"
    "Usage:   distmat <inp.aln>\n"
    "\n"
    "Options:\n"
    "  -s   skip columns with N's\n"
    "  -n   don't normalize (to account for skipped columns)\n"
    << endl;
    return;
}

/*
 * Parse arguments.
 */
void parse_arguments(const int argc, const char **argv, string &fasta_fn) {
    if (argc==1){
        usage();
        exit(1);
    }
    
    int c;
    while ((c = getopt(argc, (char *const *)argv, "hsn")) >= 0) {
        switch (c) {
            case 'h': {
                usage();
                exit(0);
                break;
            }
            case 's': {
                //skip_N_cols=true;
                break;
            }
            case 'n': {
                //comp_abs=true;
                break;
            }
            case '?': {
                cerr << "Unknown error" << endl;
                exit(1);
            }
            default: {
                cerr << "Unknown option " << c << endl;
                exit(1);
            }
        }
    }
}


/*
 * Load sequences and convert nucleotides to upper case.
 */
template <typename T>
void load_sequences(const string &fasta_fn, T &names, T &seqs) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);

    int len=0; // length of sequencing (for checking)

    while ((l = kseq_read(seq)) >= 0) {
        names.push_back(seq->name.s);
        string s(seq->seq.s);
        for (auto & c: s) {
            c = toupper(c);
        }
        if(len!=0){
            assert(len==s.size());
        }
        else{
            len=(int)s.size();
        }

        for(char &a: s){
            assert ((int)a<128);
        }

        seqs.push_back(s);
    }
    kseq_destroy(seq);
    gzclose(fp);

    assert(seqs.size()>0);
}


/*
 * Compute pileup (len x 128).
 */
template <typename T, typename U>
void compute_pileup(const T &seqs, U &pileup) {
    auto len=seqs[0].size();
    for(int i=0; i<len; i++){
        for(int j=0; j<128; i++){
            pileup[i][j]=0;
        }
    }

    for (const auto &seq: seqs){
        for(int i=0; i<len; i++){
            int c=seq[i];
            ++pileup[i][c];
        }
    }
}


/*
 * Compute consensus string (the most common char at each position, except N).
 */
template <typename T>
void compute_consensus(const T &pileup, string &consensus) {
    assert(pileup.size()==consensus.size());
    
    consensus = string('A', pileup[0].size());
    for(int i=0; i<pileup.size(); i++){
        char c='N';
        int max_freq=-1;
        const auto &column=pileup[i];
        
        for(int d=0;d<128;d++){
            if(d!='N'){
                if(column[d]>max_freq){
                    max_freq=column[d];
                    c=(char)d;
                }
            }
        }
        
        consensus[i]=c;
    }
}


/*
 * Compute mask.
 *
 * 0 - position ignored
 * N - position non-ignored, containg Ns
 * 1 - position non-ignored
 */
template <typename T>
void compute_mask(const T &pileup, string &mask) {
    //todo
    mask = string('1', pileup[0].size());
}


/*
 * Compute pair matrix (128 x 128).
 */
template <typename T>
void compute_pair_matrix(const string &seq1, const string &seq2, T &pair_matrix){
    for(int i=0; i<128; i++){
        for(int j=0; j<128; j++){
            pair_matrix[i][j]=0;
        }
    }
    int len=(int)seq1.size();
    for (int i=0;i<len;i++){
        ++pair_matrix[seq1[i]][seq2[i]];
    }
}


/*
 * Compute distances.
 */
template <typename T, typename U>
T compute_jaccard_distance(U &pair_matrix) {
}


/*
 * Print distance matrix.
 */
template <typename T, typename U>
void compute_jaccard_distance(const T &names, const U &distance_matrix, int count) {
    cout << "#taxid";
    for (const string& name : names){
        cout << "\t" << name;
    }
    cout << endl;

    for (int i=0;i<count;i++){
        cout << names[i];
        for (int j=0;j<count;j++){
            cout << "\t" << distance_matrix[i][j];
        }
        cout << endl;
    }
}



/* Compare sequences */
template<typename T>
void compare_seqs(const string &s1, const string &s2, T (&counts)[5][5], const string &ncols, bool comp_abs) {
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
        if(ncols[i]!='N'){
            ++counts[nt256_nt4[s1[i]]][nt256_nt4[s2[i]]];
        }
    }
}

template <typename T, typename U>
void dist_mat(const T &seqs, U &matrix, const string &ncols, bool comp_abs){
    int counts[5][5];
    auto no_seqs=seqs.size();
    int len=(int)seqs[0].size();
    for (int i=0;i<no_seqs;i++){
        for (int j=0;j<=i;j++){
            compare_seqs(seqs[i], seqs[j], counts, ncols, comp_abs);
            matrix[i][j]=len-(matrix[j][i]=counts[0][0]+counts[1][1]+counts[2][2]+counts[3][3]+counts[4][4]);
        }
    }
}


int main (int argc, const char **argv) {
    bool skip_N_cols=false;
    bool comp_abs=false;
    string fasta_fn;
    parse_arguments(argc, argv, fasta_fn);

    vector<string> names, seqs;
    load_sequences(argv[1], names, seqs);

    int count=(int)seqs.size();
    int len=(int)seqs[0].size();


    cerr << "Constructing empty matrices" << endl;
    vector<vector<int>> distance_matrix(count, vector<int>(count, 0));

    cerr << "Computing distance matrices" << endl;

    string ncols=string(len, 'A');
    if (skip_N_cols){
        for (auto& seq : seqs) {
            for(int i=0;i<seq.size();i++){
                if(seq[i]=='N'){
                    ncols[i]='N';
                }
            }

        }
    }
    dist_mat(seqs, distance_matrix, ncols, comp_abs);

    /*
     * print the output
     */


    return 0;
}
