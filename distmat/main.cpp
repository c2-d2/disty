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
#include <numeric>
#include <cmath>
#include <zlib.h>
#include <getopt.h>

#include "kseq.h"

using namespace std;


enum class input_t {
    ACGT,
    BINARY,
    _MAX
};

enum class n_strategy_t {
    IGNORE_PAIRWISE,
    IGNORE_PAIRWISE_NORM,
    IGNORE_GLOBALLY,
    REPLACE_MAJOR,
    REPLACE_CLOSEST,
    _MAX
};

string USAGE=
"\n"
"Program: distmat - compute distance matrix from a FASTA alignment file\n"
"\n"
"Usage:   distmat <inp.aln>\n"
"\n"
"Options:\n"
"  -n  FLOAT  skip columns having frequency of N > FLOAT [1.00]\n"
"  -i  INT    input format [0]\n"
"                 0: ACGT\n"
"                 1: 01\n"
"  -s  INT    strategy to deal with N's [0]\n"
"                 0: ignore pairwisely\n"
"                 1: ignore pairwisely and normalize\n"
"                 2: ignore globally\n"
"                 3: replace by the major allele (not implemented yet)\n"
"                 4: replace by the closest individual (not implemented yet)\n"
"  -h         print help message and exit\n";

struct params_t {
    string fasta_fn;
    input_t input;
    n_strategy_t n_strategy;
    float skip_n;

    params_t()
    :fasta_fn(""), input(input_t::ACGT), n_strategy(n_strategy_t::IGNORE_PAIRWISE), skip_n(1.0)
    {}
};

struct pair_char_t {
    int matches;
    int mismatches;
    int unknown;
};

static const uint8_t acgt_nt256_nt4[] = {
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


/*
static const uint8_t acgt_nt256_nt16[] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    1 , 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};*/

//static const uint8_t acgt_nt16_nt4[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };


static const uint8_t binary_nt256_nt4[] = {
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    0,1,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,

    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
    2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2
};


KSEQ_INIT(gzFile, gzread)


/*
 * Parse arguments.
 */
void parse_arguments(int argc, const char **argv, params_t &params) {
    if (argc==1){
        cerr << USAGE << endl;
        exit(1);
    }

    int c;
    while ((c = getopt(argc, (char *const *)argv, "hi:s:n:")) >= 0) {
        switch (c) {
            case 'h': {
                cout << USAGE << endl;
                exit(0);
            }
            case 'i': {
                int val=atoi(optarg);
                assert(val>=0);
                assert(val<(int)input_t::_MAX);
                params.input=static_cast<input_t>(val);
                break;
            }
            case 's': {
                int val=atoi(optarg);
                assert(val>=0);
                assert(val<(int)n_strategy_t::_MAX);
                params.n_strategy=static_cast<n_strategy_t>(val);
                break;
            }
            case 'n': {
                float val=atof(optarg);
                assert(val>=0.0);
                assert(val<=1.0);
                params.skip_n=val;
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

    argc -= optind;
    argv += optind;

    if(argc != 1){
        cerr << USAGE << endl;
        exit(1);
    }
    else {
        params.fasta_fn=string(argv[0]);
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
    assert (fp != nullptr);
    seq = kseq_init(fp);

    int len=0; // length of sequences (for checking)

    while ((l = kseq_read(seq)) >= 0) {
        names.push_back(seq->name.s);
        string s(seq->seq.s);
        for (auto & c: s) {
            c = toupper(c); 
        }
        if(len!=0){
            assert(len==static_cast<int>(s.size()));
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

template <typename T>
void print_sequences(T &seqs) {
    cerr << endl;
    for (auto const &s: seqs){
        cerr << s << endl;
    }
    cerr << endl;
}

/*
 * Compute pileup (len x 128).
 */
template <typename T, typename U>
void compute_pileup(const T &seqs, U &pileup) {
    assert(seqs[0].size()==pileup.size());
    assert(pileup[0].size()==128);
    auto len=seqs[0].size();
    for(int i=0; i<static_cast<int>(len); i++){
        for(int c=0; c<128; c++){
            pileup[i][c]=0;
        }
    }

    for (const auto &seq: seqs){
        for(int i=0; i<static_cast<int>(len); i++){
            unsigned char c=seq[i];
            ++pileup[i][c];
        }
    }
}

template <typename T>
void print_pileup(const T &pileup){
    assert(pileup[0].size()==128);
    for (int i=0;i<pileup.size();i++){
        cout << i;
        for(int c=0;c<128;c++){
            cout << "\t" << pileup[i][c];
        }
        cout << endl;
    }
}



/*
 * Compute consensus string (the most common char at each position, except N).
 */
template <typename T>
void compute_consensus(const T &pileup, string &consensus) {
    assert(pileup.size()==consensus.size());
    assert(pileup[0].size()==128);

    for(int i=0; i<static_cast<int>(pileup.size()); i++){
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

void print_consensus(const string &consensus){
    cout << consensus << endl;
}

/*
 * Compute mask.
 *
 * if N >= skip_n, then mask the column
 *
 * 0 - position ignored
 * N - position non-ignored, containg Ns
 * 1 - position non-ignored
 */
template <typename T>
void compute_mask(string &mask, const T &pileup, int n_thres) {
    assert(pileup.size()==mask.size());
    assert(pileup[0].size()==128);

    int column_sum=accumulate(pileup[0].begin(), pileup[0].end(), 0);

    int masked_columns=0;

    for(int i=0; i<static_cast<int>(pileup.size()); i++){
        int ns=pileup[i]['n']+pileup[i]['N'];
        if (ns >= n_thres){
            mask[i]='0';
            masked_columns++;
        }
        else{
            if(ns>0)
            {
                mask[i]='N';
            }
            else{
                mask[i]='1';
            }
        }
    }

    cerr << "Number of masked columns: " << masked_columns << " (out of " << pileup.size() << " positions, threshold: " << n_thres << " Ns, number of samples: " << column_sum << ")" << endl;
}

void print_mask(const string &mask){
    cerr << mask << endl;
}


/*
 * Compute pair matrix (128 x 128).
 */
template <typename T>
void compute_pair_matrix(const string &seq1, const string &seq2, const string &mask, T &pair_matrix){
    assert(seq1.size()==seq2.size());
    assert(seq1.size()==mask.size());
    assert(pair_matrix.size()==128);
    assert(pair_matrix[0].size()==128);

    for(int i=0; i<128; i++){
        for(int j=0; j<128; j++){
            pair_matrix[i][j]=0;
        }
    }
    int len=(int)seq1.size();
    for (int i=0;i<len;i++){
        if(mask[i]!='0'){
            ++pair_matrix[seq1[i]][seq2[i]];
        }
    }
}


/*
 * Compute characteristics of a pair matrix.
 *
 * (matches, mismatches, unknown)
 */

template <typename T>
void pair_matrix_char_acgt(T &pair_matrix, pair_char_t &pair_char) {
    assert(pair_matrix.size()==128);
    assert(pair_matrix[0].size()==128);

    pair_char.matches=0;
    pair_char.mismatches=0;
    pair_char.unknown=0;

    for(unsigned char i=0;i<128;i++){
        char n1_nt4=acgt_nt256_nt4[i];
        for(unsigned char j=0;j<128;j++){
            char n2_nt4=acgt_nt256_nt4[j];

            if (n1_nt4==4 || n2_nt4==4){
                pair_char.unknown+=pair_matrix[i][j];
            }
            else{
                if(n1_nt4==n2_nt4){
                    pair_char.matches+=pair_matrix[i][j];
                }
                else {
                    pair_char.mismatches+=pair_matrix[i][j];
                }
            }
        }
    }
    //cout << pair_char.matches << " " << pair_char.mismatches << " " << pair_char.unknown << endl;
}

template <typename T>
void pair_matrix_char_binary(T &pair_matrix, pair_char_t &pair_char) {
    assert(pair_matrix.size()==128);
    assert(pair_matrix[0].size()==128);

    pair_char.matches=0;
    pair_char.mismatches=0;
    pair_char.unknown=0;

    for(unsigned char i=0;i<128;i++){
        char n1_nt4=binary_nt256_nt4[i];
        for(unsigned char j=0;j<128;j++){
            char n2_nt4=binary_nt256_nt4[j];

            if (n1_nt4==2 || n2_nt4==2){
                pair_char.unknown++;
            }
            else{
                if(n1_nt4+n2_nt4==1)
                {
                    pair_char.matches++;
                }
                else {
                    pair_char.mismatches++;
                }
            }
        }
    }
}

void print_pair_matrix_char(pair_char_t &pair_char){
    cerr << pair_char.matches << "\t" << pair_char.mismatches << "\t" << pair_char.unknown << endl;
}

/*
 * Compute distance.
 */

int distance(const pair_char_t &pair_char) {
    return pair_char.mismatches;
}

int distance_norm(const pair_char_t &pair_char) {
    float multiplicator=1.0*(pair_char.matches+pair_char.matches+pair_char.unknown)/(pair_char.matches+pair_char.mismatches);
    return round(multiplicator * pair_char.mismatches);
}


template <typename T, typename U>
void print_distance_matrix(const T &distance_matrix, const U &names){
    assert(distance_matrix.size() == distance_matrix[0].size());
    assert(distance_matrix.size() == names.size());
    int count=static_cast<int>(distance_matrix.size());

    //cout << "";
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


int main (int argc, const char **argv) {
    params_t params;
    parse_arguments(argc, argv, params);

    cerr << "Loading sequences from " << params.fasta_fn << endl;
    vector<string> names, seqs;
    load_sequences(params.fasta_fn, names, seqs);
    //print_sequences(seqs);

    int count=(int)seqs.size();
    int len=(int)seqs[0].size();

    cerr << "Computing pileup" << endl;
    vector<vector<int>> pileup(len, vector<int>(128));
    compute_pileup(seqs, pileup);
    //print_pileup(pileup);

    cerr << "Computing consensus" << endl;
    string consensus(len, '?');
    compute_consensus(pileup, consensus);
    //print_consensus(consensus);


    cerr << "Computing mask" << endl;
    string mask(len, '?');
    if(params.n_strategy==n_strategy_t::IGNORE_GLOBALLY){
        compute_mask(mask, pileup, 1);
    }
    else{
        int min_n=ceil(count*params.skip_n);
        compute_mask(mask, pileup, min_n);
    }
    //print_mask(mask);


    cerr << "Computing distance matrix" << endl;
    vector<vector<int>> pair_matrix(128, vector<int>(128, 0));
    vector<vector<int>> distance_matrix(count, vector<int>(count, 0));
    pair_char_t pair_char;

    /*
     * For each pair:
     */

    for(int i=0;i<count;i++){
        for(int j=0;j<=i;j++){
            //cerr << "\n(" << i << "," << j << ")" << endl;

            /*
             * 1) Pair matrix
             */
            compute_pair_matrix(seqs[i], seqs[j], mask, pair_matrix);

            /*
             * 2) Characteristics
             */

            if (params.input==input_t::ACGT){
                pair_matrix_char_acgt(pair_matrix, pair_char);
            }
            else{
                if (params.input==input_t::BINARY){
                    pair_matrix_char_binary(pair_matrix, pair_char);
                }
                else{
                    assert(false);
                }
            }

            //print_pair_matrix_char(pair_char);

            /*
             * 3) Distance
             */

            int dist=0;
            if(params.n_strategy==n_strategy_t::IGNORE_PAIRWISE_NORM){
                dist=distance_norm(pair_char);
            }
            else{
                dist=distance(pair_char);
            }
            distance_matrix[i][j]=distance_matrix[j][i]=dist;
        }
    }


    print_distance_matrix(distance_matrix, names);
    
    
    return 0;
}
