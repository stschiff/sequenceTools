import std.stdio;
import std.getopt;
import std.range;
import std.algorithm;
import std.exception;
import std.string;
import std.random;
import std.conv;
import gzip;

string fwdTag, bwdTag;
string fwdFilename, bwdFilename;
string fakePrefix;
string fakePrefixQual;
string outPrefix;
long startNr, endNr;
int maxMismatch = 1;
int overlapSize = 20;
int minLength = 30;
int phredQualASCIIOffset = 33;
long[char[2]] mismatchDict;
long totalMismatches;
char[char] revNuc;
long[int] lengthDist;


void main(string[] args) {
    try {
        readArgs(args);
    }
    catch(Exception e) {
        stderr.writeln(e.msg);
        printHelp();
        return;
    }
    run();
}

void readArgs(string[] args) {
    getopt(args, std.getopt.config.caseSensitive,
           "fwdTag|f", &fwdTag,
           "bwdTag|b", &bwdTag,
           "startNr", &startNr,
           "endNr", &endNr,
           "outPrefix|o", &outPrefix,
           "fakePrefix|p", &fakePrefix,
           "minLength", &minLength,
           "maxMismatch", &maxMismatch,
           "minOverlapSize", &overlapSize);
    enforce(args.length == 3, "need two input files");
    fwdFilename = args[1];
    bwdFilename = args[2];
    enforce(endNr >= startNr);
    if(fakePrefix.length > 0) {
        auto tmp = new char[fakePrefix.length];
        tmp[] = makeQualChar(40); // highest Quality on Illumina machines
        fakePrefixQual = tmp.idup;
    }
}

void printHelp() {
    stderr.writeln("./filterTrim [options] <ForwardReads.fastq> <BackwardReads.fastq>
Options:
--fwdTag, -f
--bwdTag, -b
--fakePrefix, -p
--outPrefix, -o
--startNr
--endNr
--maxMismatch [1]
--minLength [30]
--minOverlapSize [20]");
}

void run() {
    revNuc['A'] = 'T';
    revNuc['T'] = 'A';
    revNuc['C'] = 'G';
    revNuc['G'] = 'C';
    revNuc['N'] = 'N';
    
    
    auto fwdFr = new GzipByLine(fwdFilename);
    auto bwdFr = new GzipByLine(bwdFilename);
    
    auto fwdR = new ChunkRange!(typeof(fwdFr), string, 4)(fwdFr);
    auto bwdR = new ChunkRange!(typeof(bwdFr), string, 4)(bwdFr);
    
    long fwdTagFails, bwdTagFails, nrPassed, nrMerged, nrTotal;
    long fwdTagImperfects, bwdTagImperfects;
    auto fwdOutF = new GzipOut(outPrefix ~ "_1.fastq.gz");
    auto bwdOutF = new GzipOut(outPrefix ~ "_2.fastq.gz");
    auto mergedOutF = new GzipOut(outPrefix ~ "_merged.fastq.gz");
    auto statsF = File(outPrefix ~ "_stats.txt", "w");
    
    auto cnt = 0;
    foreach(pair; zip(fwdR, bwdR)) {
        nrTotal += 1;
        cnt += 1;
        if(cnt % 10000 == 0)
            stderr.writeln("processing read nr. ", cnt);
        if(endNr > 0) {
            if(cnt < startNr)
                continue;
            if(cnt > endNr)
                break;
        }
        enforce(pair[0][0][0 .. $ - 1] == pair[1][0][0 .. $ - 1], "read pairs not in equal order");
        auto fwdSeq = fakePrefix ~ pair[0][1];
        auto fwdQual = fakePrefixQual ~ pair[0][3];
        auto bwdSeq = pair[1][1];
        auto bwdQual = pair[1][3];
        auto fwdTagSeq = fwdSeq[0 .. fwdTag.length];
        auto bwdTagSeq = bwdSeq[0 .. bwdTag.length];
        auto fwdTagFail = !seqMatch(fwdTagSeq, fwdTag);
        auto bwdTagFail = !seqMatch(bwdTagSeq, bwdTag);
        fwdTagFails += fwdTagFail;
        bwdTagFails += bwdTagFail;
        fwdTagImperfects += canFind(fwdTagSeq, 'N');
        bwdTagImperfects += canFind(bwdTagSeq, 'N');
        if(fwdTagFail || bwdTagFail)
            continue;
        nrPassed += 1;
        auto readEndDist = searchReadEndDist(fwdSeq, bwdSeq, overlapSize);
        if(readEndDist >= 0) {
            if(readEndDist <= fwdTag.length + bwdTag.length)
                continue;
            auto merged = mergeReads([fwdSeq, fwdQual], [bwdSeq, bwdQual], readEndDist);
            auto mergedSeq = merged[0];
            auto mergedQual = merged[1];
            auto l = to!int(mergedSeq.length);
            if(l !in lengthDist)
                lengthDist[l] = 0;
            lengthDist[l] += 1;
            if(l >= minLength) {
                nrMerged += 1;
                mergedOutF.compress(pair[0][0][0 .. $ - 2] ~ "\n");
                mergedOutF.compress(mergedSeq ~ "\n");
                mergedOutF.compress("+\n");
                mergedOutF.compress(mergedQual ~ "\n");
            }
        }
        else {
            fwdOutF.compress(pair[0][0] ~ "\n");
            fwdOutF.compress(fwdSeq[fwdTag.length .. $] ~ "\n");
            fwdOutF.compress(pair[0][2] ~ "\n");
            fwdOutF.compress(fwdQual[fwdTag.length .. $] ~ "\n");
            bwdOutF.compress(pair[1][0] ~ "\n");
            bwdOutF.compress(bwdSeq[bwdTag.length .. $] ~ "\n");
            bwdOutF.compress(pair[1][2] ~ "\n");
            bwdOutF.compress(bwdQual[bwdTag.length .. $] ~ "\n");
        }
    }
    mergedOutF.finish();
    fwdOutF.finish();
    bwdOutF.finish();
    logPrint(statsF, "Forward strand tag fails:\t%s", fwdTagFails);
    logPrint(statsF, "Backward strand tag fails:\t%s", bwdTagFails);
    logPrint(statsF, "Forward strand tag contains N:\t%s", fwdTagImperfects);
    logPrint(statsF, "Backward strand tag contains N:\t%s", bwdTagImperfects);
    logPrint(statsF, "Total reads:\t%s", nrTotal);
    logPrint(statsF, "Passed reads:\t%s", nrPassed);
    logPrint(statsF, "Merged reads:\t%s", nrMerged);
    logPrint(statsF, "Nr of mismatches:\t%s", totalMismatches);
    foreach(key, val; mismatchDict)
        logPrint(statsF, "Mismatch\t%s,%s\t%s", key[0], key[1], val);
    foreach(key; sort(lengthDist.keys()))
        logPrint(statsF, "Length\t%s\t%s", key, lengthDist[key]);
  
}

void logPrint(S...)(File f, S s) {
    f.writefln(s);
    stderr.writefln(s);
}

class ChunkRange(R, T, size_t L) {
  R input;
  T[L] frontElem;
  bool empty_;
  
  this(R range) {
    input = range;
    loadFront();
  }
  
  void loadFront() {
    if(!input.empty) {
      foreach(i; 0 .. L) {
        frontElem[i] = input.front;
        input.popFront();
      }
    }
    else {
      empty_ = true;
    }
  }
  
  @property T[L] front() {
    return frontElem;
  }
  
  @property bool empty() {
    return empty_;
  }
  
  void popFront() {
    loadFront();
  }
}

bool seqMatch(in char[] seq1, in char[] seq2, int nrMismatches=0) {
  enforce(seq1.length == seq2.length);
  auto sum = 0UL;
  foreach(pair; zip(seq1, seq2)) {
    if(pair[0] != 'N' && pair[1] != 'N' && pair[0] != pair[1])
      sum += 1;
    if(sum > nrMismatches)
      return false;
  }
  return true;
}

unittest {
  assert(seqMatch("AACCT", "ACCCT") == false);
  assert(seqMatch("AACCT", "ACCCT", 1) == true);
}

// the readEndDistance is defined as the distance from the start of the forward read to the end of backward read 
int searchReadEndDist(in char[] fwdSeq, in char[] bwdSeq, int overlapSize) {
  auto minDist = overlapSize;
  auto maxDist = cast(int)(fwdSeq.length + bwdSeq.length - overlapSize);
  int[] ret;
  foreach(dist; minDist .. maxDist + 1) {
    auto end = min(dist, cast(int)bwdSeq.length);
    auto offset = max(0, dist - cast(int)fwdSeq.length);
    auto cmpSeq = bwdSeq[offset .. end].reverseComplement();
    auto fwdOffset = max(0, dist - cast(int)bwdSeq.length);
    // stderr.writeln(dist, " ", cmpSeq.length, " ", fwdOffset, " ", cmpSeq.length, " ", fwdSeq.length);
    if(seqMatch(fwdSeq[fwdOffset .. fwdOffset + cast(int)cmpSeq.length], cmpSeq, maxMismatch))
      ret ~= dist;
  }
  if(ret.length == 0 || ret.length > 1)
      return -1;
  else
      return ret[0];
}

unittest {
  assert(searchReadEndDist("ACCTGCCTGC", "GCTGGTAATC", 4) == 6);
  assert(searchReadEndDist("ACCTGCCTGC", "AGCCAGGTCC", 4) == 8);
  assert(searchReadEndDist("ACCTGCCTGC", "GCAGGCTGGT", 8) == 10);
  assert(searchReadEndDist("ACCTGCCTGC", "TTTGCAGGCT", 4) == 13);
  assert(searchReadEndDist("ACCTGCCTGC", "TTGGTTGGCC", 4) == -1);
}

string[2] mergeReads(string[2] fwd, string[2] bwd, int readEndDist) {
  char[] mergedSeq;
  char[] mergedQual;
  
  auto revcmp = bwd[0].reverseComplement();
  auto revQual = bwd[1].retro.map!"a.to!char()"().array;
  foreach(pos; cast(int)fwdTag.length .. readEndDist - cast(int)bwdTag.length) {
    if(pos < readEndDist - cast(int)revcmp.length) {
      mergedSeq ~= fwd[0][pos];
      mergedQual ~= fwd[1][pos];
    }
    else {
      auto revOffset = pos - (readEndDist - cast(int)revcmp.length);
      if(pos >= cast(int)fwd[0].length) {
        mergedSeq ~= revcmp[revOffset];
        mergedQual ~= revQual[revOffset];
      }
      else {
        auto qual1 = getQualPhred(fwd[1][pos]);
        auto qual2 = getQualPhred(revQual[revOffset]);
        if(fwd[0][pos] != revcmp[revOffset]) {
          if(fwd[0][pos] != 'N' && revcmp[revOffset] != 'N') {
            totalMismatches += 1;
            char[2] pair = [fwd[0][pos], revcmp[revOffset]];
            if(pair !in mismatchDict)
              mismatchDict[pair] = 1;
            else
              mismatchDict[pair] += 1;
          }
          if(qual1 > qual2)
            mergedSeq ~= fwd[0][pos];
          if(qual1 < qual2)
            mergedSeq ~= revcmp[revOffset];
          if(qual1 == qual2) {
            auto rn = uniform(0.0, 1.0);
            mergedSeq ~= rn < 0.5 ? fwd[0][pos] : revcmp[revOffset];
          }
          mergedQual ~= makeQualChar(min(qual1, qual2));
        }
        else {
          mergedSeq ~= fwd[0][pos];
          mergedQual ~= makeQualChar(min(40,qual1 + qual2));
        }
      }
    }
  }

  // printout for debug
  // if(readEndDist >= bwd[0].length) {
  //   stderr.writeln("forward:  ", fwd[0]);
  //   auto pad = new char[max(0, readEndDist - cast(int)bwd[0].length)];
  //   pad[] = ' ';
  //   stderr.writeln("backward: ", pad[] ~ revcmp);
  //   pad.length = fwdTag.length;
  //   pad[] = ' ';
  //   stderr.writeln("merged:   ", pad ~ mergedSeq);
  //   stderr.writeln("");
  // }
  // else {
  //   auto pad = new char[bwd[0].length - readEndDist];
  //   pad[] = ' ';
  //   stderr.writeln("forward:  ", pad ~ fwd[0]);
  //   stderr.writeln("backward: ", revcmp);
  //   pad.length = bwd[0].length - readEndDist + fwdTag.length;
  //   pad[] = ' ';
  //   stderr.writeln("merged:   ", pad ~ mergedSeq);
  //   stderr.writeln("");
  // }
  
  return [mergedSeq.idup, mergedQual.idup];
}

char[] reverseComplement(in char[] seq) {
  auto ret = seq.dup;
  auto n = ret.length;
  foreach(i; 0 .. n) {
    ret[i] = revNuc[seq[n - 1 - i]];
  }
  return ret;
}

unittest {
  assert(reverseComplement("ACCTT") == "AAGGT");
  assert(reverseComplement("ACNTT") == "AANGT");
}

int getQualPhred(char qualChar) {
  return cast(int)qualChar - phredQualASCIIOffset;
}

char makeQualChar(int qual) {
  return cast(char)(qual + phredQualASCIIOffset);
}