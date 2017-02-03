/*
    EDSM: Elastic Degenerate String Matching

    Copyright (C) 2017 Chang Liu, Solon P. Pissis, Ahmad Retha and Fatima Vayani.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <ctime>
#include <divsufsort64.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_trees.hpp>
#include "edsm.hpp"

using namespace sdsl;
using namespace std;

/**
* @constructor
*/
EDSM::EDSM()
{
    this->m = 0;
    this->B = 0;
    this->f = 0;
    this->F = 0;
    this->d = 0;
    this->D = 0;
    this->Np = 0;
    this->Nm = 0;
    this->pos = 0;
    this->duration = 0;
    this->kmpBT = NULL;
    this->primed = false;
    this->chr2idx[(int)'N'] = 0;
    this->chr2idx[(int)'A'] = 1;
    this->chr2idx[(int)'C'] = 2;
    this->chr2idx[(int)'G'] = 3;
    this->chr2idx[(int)'T'] = 4;
}

/**
* The constructor can be given a pattern
*
* @constructor
* @param P The determinate pattern string to search the segments for
*/
EDSM::EDSM(const string & P) : EDSM()
{
    this->setPattern(P);
}

/**
* @destructor
*/
EDSM::~EDSM()
{
    if (this->kmpBT != NULL) {
        delete [] this->kmpBT;
    }
}

/**
* Set the pattern to search for
*
* @param P A determinate pattern consisting of A, C, G or T characters
*/
void EDSM::setPattern(const string & P)
{
    clock_t start = clock();

    if (P.length() == 0)
    {
        cerr << "Error: Pattern length 0! Aborting!" << endl;
        return;
    }
    else if (P.length() > WORDSIZE)
    {
        cerr << "Error: Pattern length greater than word size (" << WORDSIZE << ")! Aborting!" << endl;
        return;
    }
    unsigned int j;
    for (j = 0; j < P.length(); j++) {
        if (!(P[j] == 'A' || P[j] == 'C' || P[j] == 'G' || P[j] == 'T')) {
            cerr << "Error: Invalid character in pattern: '" << P[j] << "'!" << endl;
            return;
        }
    }
    this->P = P;
    this->m = j;

    //reset search state
    this->primed = false;
    this->d = 0;
    this->D = 0;

    //construct the border table for the KMP search
    this->constructKMPBT();

    //initialize I bitvector by first clearing it then filling it with positions of P characters
    for (j = 0; j < 5ul; j++) {
        I[j] = 0;
    }
    for (j = 0; j < this->m; j++)
    {
        if (this->chr2idx[(int)this->P[j]] > 0) {
            this->I[this->chr2idx[(int)this->P[j]]] = this->I[this->chr2idx[(int)this->P[j]]] | (1ul << (this->m - j));
        }
    }

    //construct the suffix tree of P
    construct_im(this->STp, this->P.c_str(), sizeof(char));

    this->duration += clock() - start;
}

/**
* Construct the border table of P for the KMP search
*/
void EDSM::constructKMPBT()
{
    if (this->kmpBT != NULL)
    {
        delete [] this->kmpBT;
    }

    this->kmpBT = new int[this->m];
    this->kmpBT[0] = -1;
    int i, j;
    for (i = 1; i < (int)this->m; i++)
    {
        j = this->kmpBT[i - 1];
        while (j >= 0)
        {
            if (this->P[j] == this->P[i - 1])
            {
                break;
            }
            else
            {
                j = this->kmpBT[j];
            }
        }
        this->kmpBT[i] = j + 1;
    }
}

/**
* Get a list of the positions of matches found.
*
* @return A vector of ints of the position numbers where matches were found
*/
vector<int> EDSM::getMatches() const
{
    return this->matches;
}

/**
* Clear the matches list
*/
void EDSM::clearMatches()
{
    this->matches.clear();
}

/**
* Get the total number of determinate segments searched so far
*/
unsigned int EDSM::getd() const
{
    return this->d;
}

/**
* Get the total number of degenerate segments searched so far
*/
unsigned int EDSM::getD() const
{
    return this->D;
}

/**
* Get the total length of determinate segments searched so far
*/
unsigned int EDSM::getf() const
{
    return this->f;
}

/**
* Get the total length of degenerate segments searched so far
*/
unsigned int EDSM::getF() const
{
    return this->F;
}

/**
* Get the total number of strings analyzed that are shorter than m
*/
unsigned int EDSM::getNp() const
{
    return this->Np;
}

/**
* Get the total length of the strings analyzed that are shorter than m
*/
unsigned int EDSM::getNm() const
{
    return this->Nm;
}

/**
* Returns the execution duration of EDSM-BV in seconds
*/
double EDSM::getDuration() const
{
    return this->duration / (double) CLOCKS_PER_SEC;
}

/**
* Finds a node u (explicitNode) in the STp where (a) is a substring of (P), then
* proceeds to encode the children of u into a bit-vector (M).
*
* @param a The substring to find in P
* @return A bit-vector
*/
WORD EDSM::occVector(const string & a)
{
    node_t explicitNode = this->STp.root();
    string::const_iterator it;
    uint64_t char_pos = 0;
    unsigned int j = 0;
    for (it = a.begin(); it != a.end(); ++it)
    {
        if (forward_search(this->STp, explicitNode, it - a.begin(), *it, char_pos) > 0) {
            j++;
        } else {
            break;
        }
    }

    //if a is present in p
    if (j == a.length())
    {
        WORD M = 0;
        this->recFindAllChildNodes(explicitNode, M);
        return M;
    }

    return 0;
}

/**
* Recursively finds leaves in the tree from node u and updates M
*
* @param u The starting node
* @param M The bitvector where leaves string lengths are encoded to
*/
void EDSM::recFindAllChildNodes(const node_t & u, WORD & M)
{
    if (this->STp.is_leaf(u))
    {
        int l = (int)this->m - this->STp.sn(u); //sn(u) gets the suffix index from the leaf node u.
        M = M | (1ul << l);
    }
    else
    {
        for (const auto & child : this->STp.children(u)) {
            this->recFindAllChildNodes(child, M);
        }
    }
}

/**
* Computes the border table of every suffix of s in S that is also a prefix of P
*
* TODO: there's likely a better more memory efficient way to code this
*
* @param S A segment (contains one or more strings)
* @return Border table bit-vector representation
*/
WORD EDSM::computePrefixBorderTable(const Segment & S)
{
    unsigned int i, km;
    int j, h;
    Segment::const_iterator s;

    string k = this->P + "_";
    for (s = S.begin(); s != S.end(); ++s) {
        if ((*s).length() > 0) {
            if ((*s).length() >= this->m) {
                k += (*s).substr((*s).length() - this->m + 1);
            } else {
                k += *s;
            }
            k += '#';
        }
    }

    km = k.length();
    int * BT = new int[km];

    BT[0] = 0;
    for (j = 1; j < (int)this->m; j++)
    {
        BT[j - 1] = this->kmpBT[j];
    }
    BT[j - 1] = 0;

    j = 0;
    for (i = this->m + 1; i < km; i++)
    {
        BT[i - 1] = j;
        while (j >= 0 && k[j] != k[i]) {
            j = (j == 0) ? -1 : BT[j - 1];
        }
        j++;
    }
    BT[km - 1] = j;

    //See Example 12 in the paper. This is based on KMP algorithm.
    WORD B = 0;
    j = (int)this->m;
    for (s = S.begin(); s != S.end(); ++s)
    {
        if ((*s).length() >= this->m) {
            j += 1 + (*s).length() - ((*s).length() - (int)this->m + 1);
        } else {
            j += 1 + (*s).length();
        }
        if (BT[j - 1] > 0) {
            B = B | (1ul << -(BT[j - 1] - (int)this->m));
            for (h = BT[j - 1] - 1; BT[h] > 0; h = BT[h - 1]) {
                B = B | (1ul << -(BT[h] - (int)this->m));
            }
        }
    }

    delete [] BT;

    return B;
}

/**
* @deprecated Because EDSM::computePrefixBorderTable() has better Big O complexity
* (linear). This method is sometimes faster in practice but it has quadratic time
* complexity.
*
* Goes through the strings in a segment and returns a bitvector representing the
* matching prefixes of the pattern in the suffix of the strings.
*
* @param S A segment
* @return A bitvector representing matches of prefix of P
*/
WORD EDSM::computeSegmentPrefixMatches(const Segment & S)
{
    int i, j, m;
    WORD B = 0;

    Segment::const_iterator s;
    for (s = S.begin(); s != S.end(); ++s) {
        if (*s == EPSILON) {
            continue;
        }
        m = (*s).length();
        if (m >= (int)this->m) {
            i = m - (int)this->m + 1;
        } else {
            i = 0;
        }
        while (i < m)
        {
            j = 0;
            while (i + j < m && this->P[j] == (*s)[i + j])
            {
                if (i + j == m - 1) {
                    B = B | (1ul << -(j + 1 - (int)this->m));
                }
                j++;
            }
            i++;
        }
    }

    return B;
}

/**
* KMP search
*
* @param needle The pattern being searched for
* @param haystack The sequence that we are searching in
* @param B The border table of the needle
* @param i The starting position in the haystack from which to start searching from
* @return -1 if the needle is not found or the index of the last character of the needle found in the haystack
*/
int EDSM::KMP(const string & needle, const string & haystack, int * B, int i)
{
    int m, n;
    m = needle.length();
    n = haystack.length();

    int j = 0;
    while (i < n)
    {
        if (j == -1)
        {
            j = 0;
            i++;
        }
        else if (haystack[i] == needle[j])
        {
            j++;
            if (j == m) {
                return i;
            }
            i++;
        }
        else
        {
            j = B[j];
        }
    }

    return -1;
}

/**
* Report a match
*
* NOTE: We do not report the index of the match currently... no reason to.
*
* @param s The segment it was discovered in
* @param i The index of the match
*/
void EDSM::report(const int s, const int i)
{
    this->matches.push_back(s);
}

/**
* Search for P in S
*
* @param S The first segment and subsequent segments
* @return Match Found or not
*/
bool EDSM::searchNextSegment(const Segment & S)
{
    //start timer
    clock_t start = clock();

    if (this->m == 0) {
        cerr << "Please set a pattern before searching!" << endl;
        return false;
    }

    // set initial match found values
    unsigned int j;
    int kmpStartPos, matchIdx = -1;
    bool reportOnce = true;
    bool matchFound = false;
    bool isDeterminateSegment = false;

    // bit-vectors to temporarily hold the state of a search mid-processing
    WORD B1, B2;

    // define iterator for a segment. 'StringI' iterator is used to access all the strings in a segment
    Segment::const_iterator stringI;

    // Do the first segment, then search the next segments afterwards
    if (!this->primed)
    {
        this->B = 0;

        for (stringI = S.begin(); stringI != S.end(); ++stringI)
        {
            //keep track of f/F counts
            if (S.size() == 1)
            {
                this->f += (*stringI).length();
                isDeterminateSegment = true;
            }
            else if (*stringI != EPSILON)
            {
                this->F += (*stringI).length();
                isDeterminateSegment = false;
            }
        }

        for (stringI = S.begin(); stringI != S.end(); ++stringI)
        {
            if ((*stringI).length() >= this->m)
            {
                kmpStartPos = 0;
                while ((matchIdx = this->KMP(this->P, *stringI, this->kmpBT, kmpStartPos)) != -1)
                {
                    if (reportOnce && matchFound) {
                        break;
                    }
                    if (isDeterminateSegment) {
                        this->report(this->pos + matchIdx, matchIdx + this->m);
                    } else {
                        this->report(this->pos, matchIdx + this->m);
                    }
                    matchFound = true;
                    kmpStartPos = 2 + matchIdx - (int)this->m;
                    if (reportOnce) {
                        break;
                    }
                }
            }
        }

        //this->B = this->computeSegmentPrefixMatches(S);
        this->B = this->computePrefixBorderTable(S);
        this->primed = true;
    }
    else
    {
        //B1 = this->computeSegmentPrefixMatches(S);
        B1 = this->computePrefixBorderTable(S);

        for (stringI = S.begin(); stringI != S.end(); ++stringI)
        {
            if (*stringI == EPSILON) {
                B1 = B1 | this->B;
                continue;
            }

            //keep track of f/F counts
            if (S.size() == 1)
            {
                this->f += (*stringI).length();
                isDeterminateSegment = true;
            }
            else
            {
                this->F += (*stringI).length();
                isDeterminateSegment = false;
            }
        }

        for (stringI = S.begin(); stringI != S.end(); ++stringI)
        {
            if (*stringI != EPSILON)
            {
                if ((*stringI).length() >= this->m)
                {
                    kmpStartPos = 0;
                    matchIdx = -1;
                    while ((matchIdx = this->KMP(this->P, *stringI, this->kmpBT, kmpStartPos)) != -1)
                    {
                        if (reportOnce && matchFound) {
                            break;
                        }
                        if (isDeterminateSegment) {
                            this->report(this->pos + matchIdx, matchIdx);
                        } else {
                            this->report(this->pos, matchIdx);
                        }
                        matchFound = true;
                        kmpStartPos = 2 + matchIdx - (int)this->m;
                        if (reportOnce) {
                            break;
                        }
                    }
                }
                if (this->B != 0)
                {
                    B2 = this->B;
                    for (j = 0; j < min((unsigned int)(*stringI).length(), this->m - 1); j++)
                    {
                        B2 = B2 & this->I[this->chr2idx[(int)(*stringI)[j]]];
                        B2 = B2 >> 1;
                        if (B2 & 1ul) {
                            if (!(reportOnce && matchFound))
                            {
                                if (isDeterminateSegment) {
                                    this->report((int)(this->pos + j), 0);
                                } else {
                                    this->report((int)this->pos, (int)j);
                                }
                                matchFound = true;
                            }
                        }
                    }

                    if ((*stringI).length() < this->m)
                    {
                        this->Np++;
                        this->Nm += (*stringI).length();
                        B2 = this->B & this->occVector(*stringI);
                        B1 = B1 | (B2 >> (*stringI).length());
                    }
                }
            }
        }

        this->B = B1;
    }

    //increment the segment counter and position counter
    if (isDeterminateSegment) {
        this->d++;
        this->pos += S[0].length();
    } else {
        this->D++;
        this->pos++;
    }

    this->duration += clock() - start;

    return matchFound;
}
