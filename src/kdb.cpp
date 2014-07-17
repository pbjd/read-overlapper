#include <cstdint>
#include <cstdlib>
#include <cstddef>
#include <cstdio>
#include <cmath>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using std::string;
using std::vector;

typedef uint16_t iid_t;

struct KmerRecord {
    iid_t iid;
    uint16_t offs;
};

typedef vector<KmerRecord> KmerEntries;

struct TargetRecord {
    string name;
    uint16_t length;
    uint32_t vecOffs;
};

class KmerPair {
public:
    KmerPair(uint16_t qoffs, uint16_t toffs);
    int offsDiff();
    uint16_t getQoffs();
    uint16_t getToffs();
    string repr();

private:
    uint16_t qoffs_;
    uint16_t toffs_;
}; 

KmerPair::KmerPair(uint16_t qoffs, uint16_t toffs) :
        qoffs_(qoffs), toffs_(toffs) {}

int KmerPair::offsDiff() {
    return qoffs_ - toffs_;
}

uint16_t KmerPair::getQoffs() {
    return qoffs_;
}

uint16_t KmerPair::getToffs() {
    return toffs_;
}

string KmerPair::repr() {
    char buf[32];
    string fmt;
    sprintf(buf, "qoff: %d toff: %d", qoffs_, toffs_);
    fmt = buf;
    return fmt;
}

class ScorePair {
public:
    ScorePair() : score_(-1) {};

    void setLengths(unsigned int qlen, unsigned int tlen);

    void addFwd(KmerPair kpair);

    void addRev(KmerPair kpair); 

    int getScore();

    void tally();
private:
    vector<KmerPair> fwd_;
    vector<KmerPair> rev_;
    uint32_t qlen_;
    uint32_t tlen_;
    int score_;
};

void ScorePair::setLengths(unsigned int qlen, unsigned int tlen) {
    qlen_ = qlen;
    tlen_ = tlen;
}

void ScorePair::addFwd(KmerPair kpair) {
    fwd_.push_back(kpair);
}

void ScorePair::addRev(KmerPair kpair) {
    rev_.insert(rev_.begin(),kpair);
}

int ScorePair::getScore() {
    return score_; 
}

void ScorePair::tally() {
    if (fwd_.size() < 5 && rev_.size() < 5)
        return; 

    uint32_t bins = std::ceil((float)(qlen_ + tlen_) / 200.0);
    uint32_t mid = tlen_ / 200;
    uint32_t maxcount = 0;
    vector<uint32_t> hist(bins, 0);
    vector<KmerPair>::iterator it;

    if (rev_.size() > 4) {
        int prevQoffs = -1;
        for (it = rev_.begin(); it != rev_.end(); ++it) {
            if (it->getQoffs() == prevQoffs)
                continue;
            uint32_t slot = (it->offsDiff())/200 + mid;
            if (slot >= bins) {
                string kpr = it->repr();
                fprintf(stderr, "Reverse slot: %d, bins: %d, mid: %d, qlen: %d, tlen: %d, %s\n", slot, bins, mid, qlen_, tlen_, kpr.data());
                exit(-1);
            }
            hist[slot]++;
            prevQoffs = it->getQoffs();
        }
        maxcount = *std::max_element(hist.begin(), hist.end());
        std::fill(hist.begin(), hist.end(), 0);
    }
    
    if (fwd_.size() > 4) {
        int prevQoffs = -1;
        for (it = fwd_.begin(); it != fwd_.end(); ++it) {
            if (it->getQoffs() == prevQoffs)
                continue;
            uint32_t slot = (it->offsDiff())/200 + mid;
            if (slot >= bins) {
                string kpr = it->repr();
                fprintf(stderr, "Forward slot: %d, bins: %d, mid: %d, qlen: %d, tlen: %d, %s\n", slot, bins, mid, qlen_, tlen_, kpr.data());
                exit(-1);
            }
            hist[slot]++;
            prevQoffs = it->getQoffs();
        }
        maxcount = std::max(maxcount, *std::max_element(hist.begin(), hist.end()));
    }

    score_ = maxcount;
}

string revComp(string seq) {
    const string bases = "ACTG";
    string::iterator curr = seq.begin();
    for (; curr != seq.end(); ++curr) {
        char& c = *curr;
        c = c == 'T' ? bases[0] :
            c == 'G' ? bases[1] :
            c == 'A' ? bases[2] :
            c == 'C' ? bases[3] : c;
    }
    return string(seq.rbegin(), seq.rend());
}

uint16_t lenFromName(const string name) {
    std::istringstream start(name.substr(0,name.find('_')));
    std::istringstream end(name.substr(name.find('_')+1));
    uint16_t i, j;
    start >> i;
    end >> j;
    return j - i;
}

uint32_t bitv(const string kmer) {
    uint32_t slot = 0;
    for (int i = 0; i < kmer.length(); i++) {
        slot <<= 2;
        switch (kmer[i]) {
            case 'A':
                //slot |= 0;
                break;
            case 'C':
                slot |= 1;
                break;
            case 'G':
                slot |= 2;
                break;
            default:  // T
                slot |= 3;
        }     
    }
    return slot;
}

uint8_t baseToShort(const char base) {
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        default:  // T
            return 3;
    }
}

void addEntry(vector< KmerEntries* > *kdb, const string &seq, const iid_t iid) {
    std::unordered_map<uint32_t,bool> seen;
    for (uint16_t offs = 0; offs < seq.length() - 14; offs++) {
        uint32_t slot = bitv(seq.substr(offs, 14));

        if (seen.count(slot))
            continue;

        if ((*kdb)[slot] == NULL)
            (*kdb)[slot] = new KmerEntries();

        KmerRecord value;
        value.iid = iid;
        value.offs = offs;
        (*kdb)[slot]->push_back(value);
        seen[slot] = 1;
    }
}

void addEntry2(vector< KmerEntries* > *kdb, const string &seq, const iid_t iid) {
    uint32_t slot = bitv(seq.substr(0, 14));
    string::const_iterator curr(seq.begin() + 14);
    string::const_iterator end(seq.end() - 14);
    uint16_t offs = 1;
    uint32_t prevSlot = slot;
    for (; curr != end; ++curr) {
        slot = (slot & 0x3FFFFFF) << 2;
        slot |= baseToShort(*curr);

        // handles simple tandem repeats, e.g., AAAAAAAAAAAAAAAAAAAAAAAAAAAA
        if (slot == prevSlot) 
            continue;

        if ((*kdb)[slot] == NULL)
            (*kdb)[slot] = new KmerEntries();

        KmerRecord value;
        value.iid = iid;
        value.offs = offs++;
        (*kdb)[slot]->push_back(value);
        prevSlot = slot;
    }
}

void findPairs(vector< KmerEntries* > *kdb, const string &seq, std::unordered_map<iid_t,ScorePair> &unsortedPairs) {
    std::unordered_map<uint32_t,bool> seen;
    uint16_t qlen = seq.length();
    for (uint16_t qoffs = 0; qoffs < qlen - 14; qoffs += 7) {
        string fseq = seq.substr(qoffs, 14);
        uint32_t slot = bitv(fseq);
        KmerEntries* targets = (*kdb)[slot];
        if (targets != NULL && !seen.count(slot)) {
            KmerEntries::iterator it;
            for (it = (*targets).begin(); it != (*targets).end(); ++it) {
                KmerRecord target = *it;
                unsortedPairs[target.iid].addFwd(KmerPair(qoffs, target.offs));
            }
            seen[slot] = 1;
        }

        string rseq = revComp(fseq);
        slot = bitv(rseq);
        targets = (*kdb)[slot];
        if (targets != NULL && !seen.count(slot)) {
            KmerEntries::iterator it;
            for (it = (*targets).begin(); it != (*targets).end(); ++it) {
                KmerRecord target = *it;
                unsortedPairs[target.iid].addRev(KmerPair(qlen-(qoffs+14), target.offs));
            }
            seen[slot] = 1;
        }
    }
}

void findPairs2(vector< KmerEntries* > *kdb, const string &seq, std::unordered_map<iid_t,ScorePair> &unsortedPairs) {
    std::unordered_map<uint32_t,bool> seen;
    string rseq = revComp(seq); 
    uint32_t fslot = bitv(seq.substr(0, 14));
    uint32_t rslot = bitv(rseq.substr(0, 14));
    string::const_iterator fcurr(seq.begin() + 7);
    string::const_iterator rcurr(rseq.begin() + 7);
    uint16_t qoffs = 7;
    uint16_t qlen = seq.length();
    for (; qoffs < qlen - 14; fcurr += 7, rcurr += 7) {
        fslot = (fslot & 0x3FFF) << 14;
        fslot |= baseToShort(*(fcurr+7)) << 12;
        fslot |= baseToShort(*(fcurr+8)) << 10;
        fslot |= baseToShort(*(fcurr+9)) << 8;
        fslot |= baseToShort(*(fcurr+10))<< 6;
        fslot |= baseToShort(*(fcurr+11))<< 4;
        fslot |= baseToShort(*(fcurr+12))<< 2;
        fslot |= baseToShort(*(fcurr+13));
        KmerEntries* targets = (*kdb)[fslot];
        if (targets != NULL && !seen.count(fslot)) {
            KmerEntries::iterator curr((*targets).begin());
            KmerEntries::iterator end((*targets).end());
            for (; curr != end; ++curr) {
                KmerRecord target = *curr;
                unsortedPairs[target.iid].addFwd(KmerPair(qoffs, target.offs));
            }
            seen[fslot] = 1;
        }

        rslot = (rslot & 0x3FFF) << 14;
        rslot |= baseToShort(*(rcurr+7)) << 12;
        rslot |= baseToShort(*(rcurr+8)) << 10;
        rslot |= baseToShort(*(rcurr+9)) << 8;
        rslot |= baseToShort(*(rcurr+10))<< 6;
        rslot |= baseToShort(*(rcurr+11))<< 4;
        rslot |= baseToShort(*(rcurr+12))<< 2;
        rslot |= baseToShort(*(rcurr+13));
        targets = (*kdb)[rslot];
        if (targets != NULL && !seen.count(rslot)) {
            KmerEntries::iterator curr((*targets).begin());
            KmerEntries::iterator end((*targets).end());
            for (; curr != end; ++curr) {
                KmerRecord target = *curr;
                unsortedPairs[target.iid].addRev(KmerPair(qlen-(qoffs+14), target.offs));
            }
            seen[rslot] = 1;
        }

        qoffs += 7;
    }
}

void findPairs3(vector< KmerEntries* > *kdb, const string &seq, std::unordered_map<iid_t,ScorePair> &unsortedPairs) {
    std::unordered_map<uint32_t,bool> seen;
    uint32_t slot = bitv(seq.substr(0, 14));
    string::const_iterator fcurr(seq.begin() + 7);
    uint16_t qoffs = 7;
    uint16_t qlen = seq.length();
    for (; qoffs < qlen - 14; fcurr += 7) {
        slot = (slot & 0x3FFF) << 14;
        slot |= baseToShort(*(fcurr+7)) << 12;
        slot |= baseToShort(*(fcurr+8)) << 10;
        slot |= baseToShort(*(fcurr+9)) << 8;
        slot |= baseToShort(*(fcurr+10))<< 6;
        slot |= baseToShort(*(fcurr+11))<< 4;
        slot |= baseToShort(*(fcurr+12))<< 2;
        slot |= baseToShort(*(fcurr+13));
        KmerEntries* targets = (*kdb)[slot];
        if (targets != NULL && !seen.count(slot)) {
            KmerEntries::iterator curr((*targets).begin());
            KmerEntries::iterator end((*targets).end());
            for (; curr != end; ++curr) {
                KmerRecord target = *curr;
                unsortedPairs[target.iid].addFwd(KmerPair(qoffs, target.offs));
            }
            seen[slot] = 1;
        }

        // this little gem reverse complements the binary k-mer slot
        uint32_t rslot = 0;
        for (int i = 0; i < 14; i++) {
            rslot <<= 2;
            rslot |= slot >> (i * 2) & 0x3 ^ 0x3;
        }

        targets = (*kdb)[rslot];
        if (targets != NULL && !seen.count(rslot)) {
            KmerEntries::iterator curr((*targets).begin());
            KmerEntries::iterator end((*targets).end());
            for (; curr != end; ++curr) {
                KmerRecord target = *curr;
                unsortedPairs[target.iid].addRev(KmerPair(qlen-(qoffs+14), target.offs));
            }
            seen[rslot] = 1;
        }

        qoffs += 7;
    }
}

void printPairs(vector< KmerEntries* > *kdb, const string qname, const string &seq, const vector<TargetRecord> &targets) {
    std::unordered_map<iid_t, ScorePair> pairs;

    //findPairs(kdb, seq, pairs);
    //findPairs2(kdb, seq, pairs);
    findPairs3(kdb, seq, pairs);

    std::multimap<int, TargetRecord> ranked;
    uint16_t qlen = seq.length();
    std::unordered_map<iid_t, ScorePair>::iterator pcurr(pairs.begin());
    std::unordered_map<iid_t, ScorePair>::iterator pend(pairs.end());
    for (; pcurr != pend; ++pcurr) {
        uint32_t tid = pcurr->first;
        ScorePair score = pcurr->second;
        TargetRecord tRec = targets[tid];
        score.setLengths(qlen, tRec.length);
        score.tally();
        ranked.insert(std::pair<int, TargetRecord>(score.getScore(), tRec));
    }

    std::multimap<int, TargetRecord>::reverse_iterator rankCur(ranked.rbegin()); 
    std::multimap<int, TargetRecord>::reverse_iterator rankEnd(ranked.rend()); 
    int count = 0;
    for (; rankCur != rankEnd; ++rankCur) {

        if (rankCur->first < 5)
            continue;

        TargetRecord tRec = rankCur->second;
        string tname = tRec.name.substr(tRec.name.find('/')+1);
        std::cout << qname << ":" << tname << " " << rankCur->first << std::endl;

        if (++count == 15)
            break;
    }
}

void storeSeq(vector<bool> &seqs, const string &seq) {
    for (int i = 0; i < seq.length(); i++) {
        switch (seq[i]) {
            case 'A':
                seqs.push_back(0);
                seqs.push_back(0);
                break;
            case 'C':
                seqs.push_back(0);
                seqs.push_back(1);
                break;
            case 'G':
                seqs.push_back(1);
                seqs.push_back(0);
                break;
            default:  // T
                seqs.push_back(1);
                seqs.push_back(1);
        }     
    }
}

void getSeq(const vector<bool> &seqs, const TargetRecord tRec, string &seq) {
    vector<bool>::const_iterator start = seqs.begin() + tRec.vecOffs;
    vector<bool>::const_iterator end = seqs.begin() + tRec.vecOffs + tRec.length;
    vector<bool>::const_iterator curr;
    seq.clear();
    for (curr = start; curr != end; curr+=2) {
        uint16_t base = *curr << 1 | *(curr+1);
        switch (base) {
            case 0:
                seq.append("A");
                break;
            case 1:
                seq.append("C");
                break;
            case 2:
                seq.append("G");
                break;
            default:
                seq.append("T");
                
        }
    }    
}

void matchKmers(string qfa, string tfa) {
    vector<TargetRecord> targets;
    vector<bool> targetSeqs;
    vector< KmerEntries* > *kdb = new vector< KmerEntries* >();
    kdb->reserve(1<<(2*14));

    for (int i = 0; i < kdb->capacity(); i++)
        (*kdb)[i] = NULL;
    
    std::ifstream ifs;
    ifs.open(tfa, std::ifstream::in);
    iid_t iid = 0;
    string seq;
    string line;
    while (ifs >> line) {
        if (line[0] == '>') {
            TargetRecord tRec;
            tRec.name = line.substr(1);
            tRec.length = lenFromName(tRec.name.substr(tRec.name.rfind("/")+1));
            tRec.vecOffs = targetSeqs.size();
            targets.push_back(tRec);
            if (seq.length() > 0) {
                //storeSeq(targetSeqs, seq);
                //addEntry(kdb, seq, iid);
                addEntry2(kdb, seq, iid);
                seq.clear();
                iid++;
            }
        } else {
            seq += line;
        }
    }

    if (seq.length() > 0) {
        //storeSeq(targetSeqs, seq);
        //addEntry(kdb, seq, iid);
        addEntry2(kdb, seq, iid);
        seq.clear();
    }

    ifs.close();

    ifs.open(qfa, std::ifstream::in);
    string qname;
    while (ifs >> line) {
        if (line[0] == '>') {
            if (seq.length() > 0) {
                printPairs(kdb, qname, seq, targets);
                seq.clear();
            }
            qname = line.substr(line.find('/')+1);
        } else {
            seq += line;
        }
    } 
    if (seq.length() > 0) 
        printPairs(kdb, qname, seq, targets);
    
    ifs.close();
}

int main(int argc, char* argv[]) {
    string qfa(argv[1]);
    string tfa(argv[2]);

    matchKmers(qfa, tfa);
}
