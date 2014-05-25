/*++

 Ovdje direktno korisim LandauVishkin te isprobavam njegovu brzinu.
 Ispadne da je brzi od mojeg Myersa, sto god pokrenem ovaj uvijek radi za 0.000s!
 Probao sam samo za jako slicne proteine, trebam probati i druge stvari.
 S druge strane, ne radi za MAX_K >= 100000 jer ne moze toliko memorije alocirati kaze on.
 On cini se unaprijed zauzme memorije s obzirom na MAX_K.
 Takodjer ako je taj MAX_K velik (npr 10000) onda on postane sporiji (al to je vjerojatno zbog alokacije memorije).

--*/

#include "stdafx.h"
#include "options.h"
#include "FASTA.h"
#include "GenomeIndex.h"
#include "SingleAligner.h"
#include "PairedAligner.h"
#include "exit.h"
#include "SeedSequencer.h"

#include "LandauVishkin.h"

#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <climits>
#include <queue>


using namespace std;


int readFastaSequences(const char* path, vector< vector<char> >* seqs,
                           unsigned char* letterIdx, char* idxToLetter, bool* inAlphabet, int &alphabetLength);


int main(int argc, char * const argv[]) {
    int kArg = -1;
    bool silent = false;
    int option;
    while ((option = getopt(argc, argv, "k:s")) >= 0) { // : is not a delimiter but indicator of parameter
        switch (option) {
        case 'k': kArg = atoi(optarg); break;
        case 's': silent = true; break;
        }
    }
    if (optind + 2 != argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: snap [options...] <queries.fasta> <target.fasta>\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "\t-s  If specified, there will be no score or alignment output (silent mode).\n");
        fprintf(stderr, "\t-k K  Sequences with score > K will be discarded."
                " Smaller k, faster calculation.\n");
        return 1;
    }


    // Alphabet information, will be constructed on fly while reading sequences
    unsigned char letterIdx[128]; //!< letterIdx[c] is index of letter c in alphabet
    char idxToLetter[128]; //!< numToLetter[i] is letter that has index i in alphabet
    bool inAlphabet[128]; // inAlphabet[c] is true if c is in alphabet
    for (int i = 0; i < 128; i++) {
        inAlphabet[i] = false;
    }
    int alphabetLength = 0;

    int readResult;
    // Read queries
    char* queriesFilepath = argv[optind];
    vector< vector<char> >* querySequences = new vector< vector<char> >();
    printf("Reading queries fasta file...\n");
    readResult = readFastaSequences(queriesFilepath, querySequences, letterIdx, idxToLetter,
                                    inAlphabet, alphabetLength);
    if (readResult) {
        printf("Error: There is no file with name %s\n", queriesFilepath);
        delete querySequences;
        return 1;
    }
    int numQueries = querySequences->size();

    // Read target
    char* targetFilepath = argv[optind+1];    
    vector< vector<char> >* targetSequences = new vector< vector<char> >();
    printf("Reading target fasta file...\n");
    readResult = readFastaSequences(targetFilepath, targetSequences, letterIdx, idxToLetter,
                                    inAlphabet, alphabetLength);
    if (readResult) {
        printf("Error: There is no file with name %s\n", targetFilepath);
        delete querySequences;
        delete targetSequences;
        return 1;
    }
    char* target = (*targetSequences)[0].data();
    int targetLength = (*targetSequences)[0].size();

    printf("Alphabet: ");
    for (int c = 0; c < 128; c++)
        if (inAlphabet[c])
            printf("%c ", c);
    printf("\n");



    // ----------------------------- MAIN CALCULATION ----------------------------- //
    printf("\nSearching...\n");
    int* scores = new int[numQueries];
    int k = kArg;
    clock_t start = clock();

    LandauVishkin<1>* lv = new LandauVishkin<>;

    for (int i = 0; i < numQueries; i++) {
        char* query = (*querySequences)[i].data();
        int queryLength = (*querySequences)[i].size();
        // Calculate score
        /**
         * They are calculating SHW here!
         * Compute the edit distance between two strings, if it is <= k, or return -1 otherwise.
         * int computeEditDistance(
         *       const char* text,
         *       int textLen,
         *       const char* pattern,
         *       int patternLen,
         *       int k)
         **/
        scores[i] = lv->computeEditDistance(target, targetLength, query, queryLength, k);        
    }

    if (!silent) {
        printf("\n");

        printf("Scores (score, position) \n");
        for (int i = 0; i < numQueries; i++)
            if (scores[i] > -1) {
                printf("%d: %d\n", i, scores[i]);
            }
        
    }

    clock_t finish = clock();
    double cpuTime = ((double)(finish-start))/CLOCKS_PER_SEC;
    printf("\nCpu time of searching: %lf\n", cpuTime);
    // ---------------------------------------------------------------------------- //

    delete querySequences;
    delete targetSequences;
    delete[] scores;

    return 0;
}



/** Reads sequences from fasta file.
 * Function is passed current alphabet information and will update it if needed.
 * @param [in] path Path to fasta file containing sequences.
 * @param [out] seqs Sequences will be stored here, each sequence as vector of indexes from alphabet.
 * @param [inout] letterIdx  Array of length 128. letterIdx[c] is index of letter c in alphabet.
 * @param [inout] inAlphabet  Array of length 128. inAlphabet[c] is true if c is in alphabet.
 * @param [inout] alphabetLength
 * @return 0 if all ok, positive number otherwise.
 */
int readFastaSequences(const char* path, vector< vector<char> >* seqs,
                       unsigned char* letterIdx, char* idxToLetter, bool* inAlphabet, int &alphabetLength) {
    seqs->clear();
    
    FILE* file = fopen(path, "r");
    if (file == 0)
        return 1;

    bool inHeader = false;
    bool inSequence = false;
    int buffSize = 4096;
    char buffer[buffSize];
    while (!feof(file)) {
        int read = fread(buffer, sizeof(char), buffSize, file);
        for (int i = 0; i < read; ++i) {
            char c = buffer[i];
            if (inHeader) { // I do nothing if in header
                if (c == '\n')
                    inHeader = false;
            } else {
                if (c == '>') {
                    inHeader = true;
                    inSequence = false;
                } else {
                    if (c == '\r' || c == '\n')
                        continue;
                    // If starting new sequence, initialize it.
                    if (inSequence == false) {
                        inSequence = true;
                        seqs->push_back(vector<char>());
                    }

                    if (!inAlphabet[c]) { // I construct alphabet on fly
                        inAlphabet[c] = true;
                        letterIdx[c] = alphabetLength;
                        idxToLetter[alphabetLength] = c;
                        alphabetLength++;
                    }
                    seqs->back().push_back(letterIdx[c]);
                }
            }
        }
    }

    fclose(file);
    return 0;
}
