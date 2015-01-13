// Switch error rate calculation
// Author: Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <set>


void parseCmdLine(int argc, char **argv);
void printUsage(char **argv);
FILE * openFile(char *filename);
void readTrioParentsPairs(char *filename, int *otherParentIdx, int numSamples);
int  getLocalAnc(FILE *hapmixFile);

// When run with -s option:
// Number of individuals in the estimated phase file to skip before the
// individuals for the comparison are reached
int skipNumInEst = 0;
// When run with -t option:
// Do trio parents appear in succession in the estimated and true phgeno files
// with the transmitted haplotype first for each parent?  If so, when the truth
// set has both parents het and the child is also het (indicated by the
// transmitted haplotypes for the two parents having different genotypes), skip
// the site and do not count switch errors for it: its pedigree-based phase is
// ambiguous.
bool trioParentsInSuccession = false;
// Specified with -p option:
// Filename that specifies which pairs of samples in the input estimated and
// true phgeno files are parents of trios.
char *trioParentsFilename = NULL;
// Specified with -v option:
bool verbose = false;
// When run with -o option:
// File listing the individual numbers (starting from 0) to omit from
// <estimated phgeno> and thus not compare them with <true phgeno>.  Numbers
// begin after any individuals have been skipped.
char *omitIndFile = NULL;
// Set of individuals to omit from <estimated phgeno>, read from omitIndFile
std::set<int> omitIndSet;
// When run with -l option:
// Stratifies switch errors by local ancestry status:
// homozygous POP1 (typically European), heterozygous, homozygous POP2
// We only count sites where the HAPMIX has > .9 posterior probability of
// having the indicated local ancestry status.  We also require both sides of
// a switch error to occur in the confident block of local ancestry for the
// switch error to count -- switch errors that straddle local ancestry classes
// or that occur in regions where HAPMIX is not confident of the local ancestry
// do not count.
char *hapmixLocalAncFilesPrefix = NULL;
// Specified with -c option and necessary for local ancestry:
int chrom;


int main(int argc, char **argv) {
  parseCmdLine(argc, argv);

  // should only have 3 more arguments besides any options that were already
  // parsed with parseCmdLine(), which updates optind to be the first index
  // in argv for a non-option argument.
  if (argc - optind != 3) {
    printUsage(argv);
  }


  int numSamples = atoi(argv[optind]);

  char *estGenoFile = argv[optind+1];
  char *trueGenoFile = argv[optind+2];

  FILE *estG  = openFile(estGenoFile);
  FILE *trueG = openFile(trueGenoFile);

  // When trioParentsFilename != NULL, otherParentIdx[idx] == the index of the
  // spouse of idx
  int *otherParentIdx = NULL;

  if (trioParentsFilename != NULL) {
    otherParentIdx = new int[numSamples];
    for(int i = 0; i < numSamples; i++) { otherParentIdx[i] = -1; } // init
    readTrioParentsPairs(trioParentsFilename, otherParentIdx, numSamples);
  }
  
  if (omitIndFile != NULL) {
    FILE *omitIn = openFile(omitIndFile);
    int id;
    while (fscanf(omitIn, "%d", &id) == 1) {
      assert(id >= 0);
      omitIndSet.insert(id);
    }
    fclose(omitIn);
  }

  bool useLocalAnc = false;
  FILE **hapmixFiles = NULL;
  // when using local ancestry, specifies previous ancestry class:
  // -1 => ambiguous/unknown; 0 => homozy POP1, 1 => het, 2 => homozy POP2
  int *prevLocalAnc = NULL;
  if (hapmixLocalAncFilesPrefix != NULL) {
    useLocalAnc = true;

    if (numSamples > 1000) {
      fprintf(stderr, "Warning: limitations on the number of open files may prevent\n");
      fprintf(stderr, "the program from opening all HAPMIX output files and cause a crash\n");
      fprintf(stderr, "try running ulimit -n if this occurs\n\n");
    }

    hapmixFiles = new FILE*[numSamples];
    prevLocalAnc = new int[numSamples];
    char *filename = new char[strlen(hapmixLocalAncFilesPrefix)+10];
    for(int i = 0; i < numSamples; i++) {
      sprintf(filename, "%s.%d.%d", hapmixLocalAncFilesPrefix, i, chrom);
      hapmixFiles[i] = openFile(filename);
      prevLocalAnc[i] = -1;
    }
    delete [] filename;
  }


  char *allEstAlleles = new char[2*numSamples];
  char *allTrueAlleles = new char[2*numSamples];
  int *homologsInverted = new int[numSamples];
  int *prevSwitchError = new int[numSamples];
  int *indivNumSwitches = new int[numSamples];
  for(int i = 0; i < numSamples; i++) {
    homologsInverted[i] = -1;
    prevSwitchError[i] = 0;
    indivNumSwitches[i] = 0;
  }

  int numMissing = 0;
  int numSwitchErrors = 0;
  int totalHetSites = 0;
  int numMarkers = 0;

  // when doing switch errors stratified by ancestry; last class is ambiguous:
  int numAncClassSwitchErrors[4], totalAncClassHetSites[4];
  for(int i = 0; i < 4; i++) {
    numAncClassSwitchErrors[i] = totalAncClassHetSites[i] = 0;
  }

  bool oneHapTruthWarn = false;

  // read in each sample
  int c;
  while ((c = fgetc(estG)) != EOF) {
    ungetc(c, estG);

    numMarkers++;

    //////////////////////////////////////////////////////////////////////////
    // skip the initial samples in the file that the user specified to skip
    // (if not set, <skipNumInEst> defaults to 0)
    int i = 0;
    while (i < 2 * skipNumInEst && (c = fgetc(estG)) != '\n') {
      i++;
    }
    assert(i == 2 * skipNumInEst);

    //////////////////////////////////////////////////////////////////////////
    // read all the estimated alleles for the current SNP:
    i = 0;
    int curHap = 0; // the real haplotype number (including omitted ones)
    while ((c = fgetc(estG)) != '\n' && i < 2*numSamples) {
      int curSamp = curHap / 2;
      curHap++;
      if (omitIndSet.find(curSamp) != omitIndSet.end()) {
	continue; // omit this sample / haplotype
      }

      allEstAlleles[i] = c;
      i++;
    }
    assert(i == 2*numSamples);
    if (c == '\n')
      // put endline back so that we don't read the next line later
      ungetc(c, estG);

    // read all the true alleles for the current SNP:
    i = 0;
    while ((c = fgetc(trueG)) != '\n' && i < 2*numSamples) {
      allTrueAlleles[i] = c;
      i++;
    }
    assert(i == 2*numSamples);
    if (c == '\n')
      // put endline back so that we don't read the next line later
      ungetc(c, trueG);


    //////////////////////////////////////////////////////////////////////////
    // compare phase for this SNP:
    for(int samp = 0; samp < numSamples; samp++) {
      int prevAndCurAncClass = -1;
      if (useLocalAnc) {
	int curLocalAncClass = getLocalAnc(hapmixFiles[samp]);
	int prevLocalAncClass = prevLocalAnc[samp];
	// update for next SNP:
	prevLocalAnc[samp] = curLocalAncClass;

	if (curLocalAncClass == prevLocalAncClass && curLocalAncClass >= 0)
	  // same class across sites
	  prevAndCurAncClass = curLocalAncClass;
	else
	  // differing ancestry class: ambiguous, won't count anc class switches
	  prevAndCurAncClass = -1;
      }

      char estAlleles[2];
      char trueAlleles[2];

      // read alleles for these samples
      for(int h = 0; h < 2; h++) {
//	estAlleles[h] = fgetc(estG);
//	trueAlleles[h] = fgetc(trueG);
	estAlleles[h] = allEstAlleles[samp*2+h];
	trueAlleles[h] = allTrueAlleles[samp*2+h];
	if (!(trueAlleles[h] == '0' || trueAlleles[h] == '1' || trueAlleles[h] == '9')) {
	  printf("At locus %d!\n", numMarkers-1);
	}
	assert(trueAlleles[h] == '0' || trueAlleles[h] == '1' || trueAlleles[h] == '9');
      }

      if (trueAlleles[0] == '9' || trueAlleles[1] == '9') {//missing data? skip!
	if (!oneHapTruthWarn &&
	    (trueAlleles[0] != '9' || trueAlleles[1] != '9')) {
	  fprintf(stderr, "Warning: missing data for only one haplotype in truth set\n");
	  oneHapTruthWarn = true;
	}
	continue;
      }

      // output from phaser should not have missing data:
      assert(estAlleles[0] != '9' && estAlleles[1] != '9');


      // check the trio is triple het at this site; if so, skip it:
      if ((trioParentsInSuccession && samp % 2 == 0) || // first of two parents?
	  (otherParentIdx != NULL)) {
	// Note: this code assumes the parents appear side by side in the file
	// but they need not be so; will modify this later
	int otherParent;
	if (trioParentsInSuccession)
	  otherParent = samp+1;
	else
	  otherParent = otherParentIdx[samp];

	if (trueAlleles[0] != trueAlleles[1] &&
	    allTrueAlleles[otherParent*2] != allTrueAlleles[otherParent*2+1] &&
	    trueAlleles[0] != allTrueAlleles[otherParent*2]) {
	  // Note: trueAlleles[0] != allTrueAlleles[otherParent*2] implies the
	  // child is heterozygous since the first haplotype in each case is
	  // the one that was transmitted to the child

	  // Triple het site!  Ambiguous.
	  if (trioParentsInSuccession) {
	    // if next person is the other parent: skip this sample and the
	    // next one -- incrementing here and continuing suffices since
	    // continuing will increment samp once more
	    samp++;
	  }

	  continue;
	}
      }


      // Missing estimated haplotype?
      if (estAlleles[0] == '?' || estAlleles[1] == '?') {
	assert(estAlleles[0] == estAlleles[1]);
	numMissing++;
	continue;
      }

      if (homologsInverted[samp] == -1) {
	// haven't yet encountered heterozygous locus, so identify which
	// homolog we're on if this is heterozygous:
	if (trueAlleles[0] == trueAlleles[1]) {
	  // homozygous, so no info yet
//	  printf("est: %c/%c, true: %c/%c, samp: %d, locus: %d\n",
//		 estAlleles[0], estAlleles[1], trueAlleles[0], trueAlleles[1],
//		 samp, numMarkers-1);
	  if (estAlleles[0] != estAlleles[1] ||
					      estAlleles[0] != trueAlleles[0]) {
	    printf("At locus %d, samp %d: true: %d/%d est: %d/%d.\n",
		   numMarkers-1, samp, trueAlleles[0], trueAlleles[1],
		   estAlleles[0], estAlleles[1]);
	  }
	  assert(estAlleles[0] == estAlleles[1] &&
		 estAlleles[0] == trueAlleles[0]);
	}
	else {
	  // don't count the first het site as a het site: it cannot be switched
	  // relative to the previous locus since there isn't a previous locus
	  //totalHetSites++;

	  // can now determine the correspondence of the homologs:
	  if (estAlleles[0] == trueAlleles[0]) {
	    homologsInverted[samp] = 0;
	    assert(estAlleles[1] == trueAlleles[1]);
	  }
	  else {
	    assert(estAlleles[0] == trueAlleles[1] &&
		   estAlleles[1] == trueAlleles[0]);
	    homologsInverted[samp] = 1;
	  }
	}
      }
      else {
	// check for switch error:
	int h0 = homologsInverted[samp];
	int h1 = 1 - homologsInverted[samp];

	if (trueAlleles[0] != trueAlleles[1]) {
	  totalHetSites++;
	  if (prevAndCurAncClass >= 0)
	    totalAncClassHetSites[ prevAndCurAncClass ]++;
	  else
	    totalAncClassHetSites[ 3 ]++;
	}


	if (estAlleles[h0] == trueAlleles[0]) {
	  // match!
	  assert(estAlleles[h1] == trueAlleles[1]);
	}
	else {
	  // mismatch!
//	  printf("mismatch: est: %c/%c, true: %c/%c, samp: %d, locus: %d\n",
//		 estAlleles[h0], estAlleles[h1], trueAlleles[0], trueAlleles[1],
//		 samp, numMarkers-1);
	  if (estAlleles[h0] != trueAlleles[1] ||
					    estAlleles[h1] != trueAlleles[0]) {
	    fprintf(stderr, "samp = %d locus = %d; est = %c/%c true = %c/%c\n",
		    samp, numMarkers-1, estAlleles[h0], estAlleles[h1],
		    trueAlleles[0], trueAlleles[1]);
	  }

	  assert(estAlleles[h0] == trueAlleles[1] &&
		 estAlleles[h1] == trueAlleles[0]);
	  numSwitchErrors++;
	  if (prevAndCurAncClass >= 0)
	    numAncClassSwitchErrors[ prevAndCurAncClass ]++;
	  else
	    numAncClassSwitchErrors[ 3 ]++;
	  // invert for next locus:
	  homologsInverted[samp] = 1 - homologsInverted[samp];
	  // calculate length of chunk:
	  int locus = numMarkers - 1;
	  if (verbose) {
	    int blockLength = locus - prevSwitchError[samp];
	    fprintf(stderr, "%d %d %d %d\n", samp, indivNumSwitches[samp],
		    locus, blockLength);
	  }
	  prevSwitchError[samp] = locus; // set new switch err location
	  indivNumSwitches[samp]++;
	}
      }
    }

    // read to endline:
    while ((c = getc(estG)) != '\n');
    while ((c = getc(trueG)) != '\n');
  }

  if (verbose) { // print length of last chunk
    for(int samp = 0; samp < numSamples; samp++) {
      // calculate length of chunk:
      int locus = numMarkers - 1;
      int blockLength = locus - prevSwitchError[samp];
      fprintf(stderr, "%d %d %d %d\n", samp, indivNumSwitches[samp], locus,
	      blockLength);
    }
  }


  // subtract 1 here since the first het locus is arbitrary
  int denom1 = totalHetSites;
  int denom2 = numSamples * numMarkers;
  printf("switch %d / %d = %lf\n",
         numSwitchErrors, denom1, (double) numSwitchErrors / denom1);
  if (numMissing > 0)
    printf("missing %d / %d = %lf\n",
	   numMissing, denom2, (double) numMissing / denom2);

  if (useLocalAnc) {
    printf("Homozy_POP1:  %d / %d = %lf\n",
	   numAncClassSwitchErrors[ 0 ], totalAncClassHetSites[ 0 ],
	   (double) numAncClassSwitchErrors[ 0 ]  / totalAncClassHetSites[ 0 ]);
    printf("Heterozygous: %d / %d = %lf\n",
	   numAncClassSwitchErrors[ 1 ], totalAncClassHetSites[ 1 ],
	   (double) numAncClassSwitchErrors[ 1 ]  / totalAncClassHetSites[ 1 ]);
    printf("Homozy_POP2:  %d / %d = %lf\n",
	   numAncClassSwitchErrors[ 2 ], totalAncClassHetSites[ 2 ],
	   (double) numAncClassSwitchErrors[ 2 ]  / totalAncClassHetSites[ 2 ]);
    printf("Ambiguous:    %d / %d = %lf\n",
	   numAncClassSwitchErrors[ 3 ], totalAncClassHetSites[ 3 ],
	   (double) numAncClassSwitchErrors[ 3 ]  / totalAncClassHetSites[ 3 ]);
  }

  return 0;
}

void parseCmdLine(int argc, char **argv) {
  int c;

  while ((c = getopt(argc, argv, "s:tp:vo:l:c:")) != -1) {
    switch (c) {
      case 's':
	skipNumInEst = atoi(optarg);
	break;
      case 't':
	trioParentsInSuccession = true;
	break;
      case 'p':
	trioParentsFilename = optarg;
	break;
      case 'v': // verbose: print info about where switches occur
	verbose = true;
	break;
      case 'o':
	omitIndFile = optarg;
	break;
      case 'l':
	hapmixLocalAncFilesPrefix = optarg;
	break;
      case 'c':
	chrom = atoi(optarg);
	break;

      default:
	abort();
    }
  }
}

void printUsage(char **argv) {
  fprintf(stderr, "Usage: %s [OPTIONS] <num estimated> <estimated phgeno> <true phgeno>\n\n",
	  argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -s <#>       Skip specified number of samples in estimated file\n");
  fprintf(stderr, "  -t           Trio aware, trio parents in succession; omits triple hets\n");
  fprintf(stderr, "  -p <file>    Trio aware, <file> gives parent relationships; omits triple hets\n");
  fprintf(stderr, "  -v           Verbose: prints switch point information to stderr\n");
  fprintf(stderr, "  -o <file>    Skips/omits given ind numbers <estimated phgeno> from comparison\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -l <prefix>  Local ancestry aware, <prefix> specifies HAPMIX local ancestry\n");
  fprintf(stderr, "  -c <#>       For local ancestry: need suffix of chromosome number\n");
  exit(1);
}

void readTrioParentsPairs(char *filename, int *otherParentIdx, int numSamples) {
  FILE *in = openFile(filename);
  int id1, id2;

  int numPairs = 0;
  while (fscanf(in, "%d %d", &id1, &id2) == 2) {
    numPairs++;

    assert(otherParentIdx[id1] == -1 && otherParentIdx[id2] == -1);
    assert(id1 != id2 && id1 < numSamples && id2 < numSamples);
    otherParentIdx[id1] = id2;
    otherParentIdx[id2] = id1;
  }

  assert(numPairs * 2 == numSamples);

  fclose(in);
}

FILE * openFile(char *filename) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    fprintf(stderr, "Error: Couldn't open %s\n", filename);
    perror(filename);
    exit(1);
  }

  return file;
}

int getLocalAnc(FILE *hapmixFile) {
  int pos;
  float homozyPop1, het, homozyPop2;
  int ret = fscanf(hapmixFile, "%d %f %f %f", &pos, &homozyPop1, &het,
                   &homozyPop2);
  if (ret != 4) {
    fprintf(stderr, "Error, malformed line in local ancestry file\n");
    exit(1);
  }
  float sum = homozyPop1 + het + homozyPop2;
  assert(sum >= 0.997f && sum <= 1.003f);

  if (homozyPop1 > 0.9f)
    return 0;
  else if (het > 0.9f)
    return 1;
  else if (homozyPop2 > 0.9f)
    return 2;
  else
    return -1;
}
