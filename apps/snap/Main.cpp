/*++

Module Name:

    snap.cpp

Abstract:

    Main entry point for the snap binary. Calls into the indexer, single-end aligner or
    paired-end aligner functions defined in other source files in this directory.

Authors:

    Matei Zaharia & Bill Bolosky, February, 2012

Environment:

    User mode service.

Revision History:

    Adapted from cSNAP, which was in turn adapted from the scala prototype

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


using namespace std;

const char *SNAP_VERSION = "1.0beta.10";

static void usage()
{
    fprintf(stderr,
            "Usage: snap <command> [<options>]\n"
            "Commands:\n"
            "   index    build a genome index\n"
            "   single   align single-end reads\n"
            "   paired   align paired-end reads\n"
            "Type a command without arguments to see its help.\n");
    soft_exit(1);
}

int main(int argc, const char **argv)
{
    //InitializeSeedSequencers();  // TODO: Do I need this

    LandauVishkin<1>* lv = new LandauVishkin<>;

    const char* text = "ACTGACACACGGGCCC";
    int textLen = 16;
    const char* pattern = "TGACACAGGGC";
    int patternLen = 11;

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
    int ed = lv->computeEditDistance(text, textLen, pattern, patternLen, 10);
    printf("%d\n", ed);

    int ed2 = lv->computeEditDistance(text, textLen, text, textLen, 5);
    printf("%d\n", ed2);
}
