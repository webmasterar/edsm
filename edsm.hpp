
#ifndef __EDSM__
#define __EDSM__

#include <cstdlib>
#include <string>
#include <vector>
#include <sdsl/util.hpp>
#include <sdsl/suffix_trees.hpp>

typedef sdsl::cst_sct3<> cst_t;
typedef sdsl::cst_sada<> cst_d;
typedef cst_t::node_type node_t;
typedef cst_t::size_type size_type;
typedef cst_t::char_type char_t;
typedef std::vector<std::string> Segment;
typedef std::vector<Segment> GenIndSeq;
#define WORD unsigned long int
#define WORDSIZE sizeof(WORD) * 8
#define EPSILON "E"
#define BUFFERSIZE 1000000

class EDSM
{
private:

    void recFindAllChildNodes(const node_t & n, WORD & x);

protected:

    /**
     * @var matches All matches found as tuples of segment number and ending-position
     */
    //std::vector<std::tuple<const int, const int>>
    std::vector<int> matches;

    /**
     * @var P The determinate pattern to find in T
     */
    std::string P;

    /**
     * @var m The length of P
     */
    unsigned int m;

    /**
     * @var f The total length of determinate segments searched
     */
    unsigned int f;

    /**
     * @var F the total length of degenerate segments searched
     */
    unsigned int F;

    /**
     * @var d The number of determinate segments searched so far
     */
    unsigned int d;

    /**
     * @var D The number of degenerate segments searched so far
     */
    unsigned int D;

    /**
     * @var pos The current position in the input file being read. Degenerate
     * positions/segments count as 1 position, whereas determinate segments are
     * counted as N positions (as many characters as there are in the segment).
     */
    unsigned int pos;

    /**
     * @var duration The amount of time spent by EDSM-BV
     */
    double duration;

    /**
     * @var STp The suffix tree (array) of P
     */
    cst_t STp;

    /**
     * @var kmpBT The border table for P used in the KMP search
     */
    int * kmpBT;

    /**
     * @var chr2idx A lookup array used in combination with I to convert characters to indexes of I
     */
    int chr2idx[117] = {0};

    /**
     * @var primed Has the algorithm been primed with an initial segment to search?
     */
    bool primed;

    /**
     * @var B Bitvector maintaining the current state of the search
     */
    WORD B;

    /**
     * @var I A bitvector representing the positions of the letters in P
     *
     * Imagine P="ATAACTG", then I would hold:
     * I[*]: 00000000
     * I[A]: 10110000
     * I[C]: 00001000
     * I[G]: 00000010
     * I[T]: 01000100
     *
     * Note that index 'A' is actually 1 to 'T' 4. Also, the bits are shifted one position further to the left.
     */
    WORD I[5] = {0};

    WORD computeSegmentPrefixMatches(const Segment & S);

    WORD computePrefixBorderTable(const Segment & S);

    void constructKMPBT();

    int KMP(const std::string & needle, const std::string & haystack, int * B, int i);

    void report(const int s, const int i);

    WORD occVector(const std::string & a);

    void setPattern(const std::string & P);

public:

    EDSM();

    EDSM(const std::string & P);

    ~EDSM();

    bool searchNextSegment(const Segment & S);

    std::vector<int> getMatches() const;

    void clearMatches();

    double getDuration() const;

    unsigned int getd() const;

    unsigned int getD() const;

    unsigned int getf() const;

    unsigned int getF() const;

};

#endif

