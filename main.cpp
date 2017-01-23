
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include "edsm.hpp"
#include <Variant.h>

using namespace std;
using namespace vcflib;

char * BUFF = new char[BUFFERSIZE];
int BUFFLIMIT = 0;
int POS = 0;

/*
* Buffered file reading char by char
*
* @param f opened file handle
*/
char getNextChar(ifstream & f)
{
    if (BUFFLIMIT == 0) {
        f.read(BUFF, BUFFERSIZE);
        BUFFLIMIT = f.gcount();
        POS = 0;
        if (BUFFLIMIT == 0) {
            return '\0';
        }
    }
    char c = BUFF[POS];
    if (++POS == BUFFLIMIT) {
        BUFFLIMIT = 0;
    }
    return c;
}

int main(int argc, char * argv[])
{
    string help = "There are two ways to run Elastic Degenerate String Matching (EDSM) ---\n\
    \tUsage: ./edsm seq.txt pattern\n\
    \tUsage: ./edsm reference.fasta variants.vcf pattern";

    if (argc == 1 || (argc == 2 && (strcmp("--help", argv[1]) == 0 || strcmp("-h", argv[1]) == 0))) {
        cout << help << endl;
        return 0;
    }

    if (!(argc == 4 || argc == 3)) {
        cerr << "Invalid number of arguments!" << endl;
        cout << help << endl;
        return 1;
    }

	//pattern p
    string p;
    if (argc == 3) {
        p = argv[2];
    } else {
        p = argv[3];
    }
    //if user is passing pattern file instead of literal pattern try to read the pattern from the file
    if (p.find('.') != string::npos) {
        ifstream pf(p.c_str(), ios::in);
        if (!pf.good()) {
            cerr << "Error: Failed to open pattern file!" << endl;
            return 1;
        }
        char c;
        p = "";
        while(pf.get(c))
        {
            if (!(c == '\0' || c == '\r' || c == '\n')) {
                p += c;
            }
        }
        pf.close();
    }

    EDSM edsm(p);

    cout << "EDSM-BV searching..." << endl << endl;

    if (argc == 4)
    {
        string refName = argv[1];
        ifstream rf(refName.c_str(), ios::in);
        if (!rf.good()) {
            cerr << "Error: Failed to open reference file!" << endl;
            return 1;
        }

        string vcfName = argv[2];
        VariantCallFile vf;
        vf.open(vcfName);
        if (!vf.is_open()) {
            cerr << "Error: Failed to open variants file!" << endl;
            return 1;
        }

        //initialize fasta file reading and segment creation helper variables
        string tBuff = "";
        tBuff.reserve(BUFFERSIZE + p.length());
        char c;
        unsigned int rfIdx = 1, vfIdx = 0, i = 0;
        Segment segment;

        //skip first line of fasta file
        getline(rf, tBuff);

        //create variables for reading through vcf records and looking for duplicates
        Variant var(vf), vBuffer(vf), var2(vf);
        bool hasMoreVariants = true;
        Segment vAlleles;

        //read first variant and possibly successive duplicates for the same position, removing duplicate alleles
        hasMoreVariants = vf.getNextVariant(var);
        if (hasMoreVariants)
        {
            vfIdx = (unsigned int) var.position;
            for (const auto & a : var.alleles) {
                if (a[0] != '<') {
                    vAlleles.push_back(a);
                }
            }
            while (true) {
                hasMoreVariants = vf.getNextVariant(var2);
                if (hasMoreVariants)
                {
                    if (var2.position == var.position) {
                        for (const auto & a : var2.alt) {
                            if (a[0] != '<') {
                                vAlleles.push_back(a);
                            }
                        }
                    } else {
                        vBuffer = var2;
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
        }

        //go through the reference sequence
        while ((c = getNextChar(rf)) != '\0')
        {
            if (!(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')) {
                continue;
            }

            if (rfIdx != vfIdx)
            {
                tBuff += c;
                i++;
                if (i >= BUFFERSIZE) {
                    segment.clear();
                    segment.push_back(tBuff);
                    edsm.searchNextSegment(segment);
                    tBuff = "";
                    i = 0;
                }
            }
            else
            {
                segment.clear();
                if (tBuff.length() > 0) {
                    segment.push_back(tBuff);
                    edsm.searchNextSegment(segment);
                    tBuff = "";
                }

                //then search current variant
                if (vAlleles.size() > 0) {
                    edsm.searchNextSegment(vAlleles);
                    vAlleles.clear();
                }

                //fetch the next variant to be searched for when its position comes up
                if (vBuffer.alleles.size() > 0)
                {
                    vfIdx = (unsigned int) vBuffer.position;
                    for (const auto & a : vBuffer.alleles) {
                        if (a[0] != '<') {
                            vAlleles.push_back(a);
                        }
                    }

                    vBuffer.alleles.clear();

                    while (true) {
                        hasMoreVariants = vf.getNextVariant(var2);
                        if (hasMoreVariants)
                        {
                            if (vfIdx == (unsigned int) var2.position) {
                                for (const auto & a : var2.alt) {
                                    if (a[0] != '<') {
                                        vAlleles.push_back(a);
                                    }
                                }
                            } else {
                                vBuffer = var2;
                                break;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }

            rfIdx++;
        }
        if (tBuff.length() > 0)
        {
            segment.push_back(tBuff);
            edsm.searchNextSegment(segment);
            tBuff = "";
        }

        rf.close();
    }
    else
    {
        //open sequence file
        ifstream eds(argv[1], ios::in);
        if (!eds.good()) {
            cerr << "Error. Unable to open sequence file!" << endl;
            return 1;
        }

        //initialize variables required for searching
        Segment tempSeg;
        string x = "";
        x.reserve(BUFFERSIZE + p.length());
        char c = 0;
        bool inDegSeg = false;
        unsigned int i = 0, j = 0;

		//go through the sequence file
        while ((c = getNextChar(eds)) != '\0')
        {
            if (i == 0 && c == '{')
            {
                inDegSeg = true;
            }
            else if (c == ',')
            {
                tempSeg.push_back(x);
                x = "";
                j = 0;
            }
            else if (c == '}' || (c == '{' && i > 0))
            {
                if (x.length() > 0) {
                    tempSeg.push_back(x);
                    x = "";
                    j = 0;
                    edsm.searchNextSegment(tempSeg);
                    tempSeg.clear();
                }
                inDegSeg = (c == '{');
            }
            else if (c != '{')
            {
                switch (c) {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case EPSILON[0]:
                    x += c;
                    j++;
                    if (!inDegSeg && j == BUFFERSIZE)
                    {
                        tempSeg.push_back(x);
                        x = "";
                        j = 0;
                        edsm.searchNextSegment(tempSeg);
                        tempSeg.clear();
                    }
                    else if (inDegSeg && j == (BUFFERSIZE + p.length()))
                    {
                        tempSeg.push_back(x.substr(0, BUFFERSIZE));
                        x = x.substr(BUFFERSIZE - 1, p.length());
                        j = p.length();
                    }
                    break;
                }
            }
            i++;
        }
        if (x != "") {
            tempSeg.push_back(x);
            edsm.searchNextSegment(tempSeg);
            x = "";
            tempSeg.clear();
        }

        eds.close();
    }


    //output results

    cout << "No. determinate bases (f): " << edsm.getf() << endl;
    cout << "No. degenerate bases (F): " << edsm.getF() << endl;
    cout << "No. determinate segments (d): " << edsm.getd() << endl;
    cout << "No. degenerate segments (D): " << edsm.getD() << endl;
    cout << "EDSM-BV searching duration: " << edsm.getDuration() << "s." << endl << endl;

    if (edsm.getMatches().size() >= 1)
    {
        cout << "Matches found: " << edsm.getMatches().size() << endl << endl;
        cout << "Positions" << endl << "---------" << endl;
        for (const auto & a : edsm.getMatches()) {
            cout << a << endl;
        }
    }
    else
    {
        cout << "No matches found." << endl;
    }

    delete [] BUFF;

    return 0;
}

