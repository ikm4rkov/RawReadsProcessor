#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <limits.h>
#include <cstdlib>
#include <cstring>
#include <map>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

// defining read for reading, 3 for wrinting, one empty for last buffer drop and vector for cleaning IDs
std::vector<std::string> one_read = { "", "", "", "" };
std::vector<std::string> first_read = { "", "", "", "" };
std::vector<std::string> second_read = { "", "", "", "" };
std::vector<std::string> dna_read = { "", "", "", "" };
std::vector<std::string> dna_read_clear = { "", "", "", "" };
std::vector<std::string> rna_read = { "", "", "", "" };
std::vector<std::string> rna_read_clear = { "", "", "", "" };
std::vector<std::string> read_nb = { "", "", "", "" };
std::vector<std::string> empty_vec = {};
std::vector<std::string> split_id = {};
std::vector<std::string> name_split = {};
std::vector<std::string> cond_codes = {};
std::set<string> condi_codes;
map<string, int>count_SE_codes;
map<string, int>count_PE_codes;
// that way with a buffer was unsuccessful
/*int buffersize = 10*1024*1024;
char *buffer1 = (char*) malloc(buffersize*sizeof(char));
int curpoi = 0;*/

// global variables for bridge coordinates
int brstart = 0;
int brend = 0;
int last_start;
int last_end;

// search function, returns nothing but can change values
void bitap_approximate_search(string pattern, string text, int mismatch_limit, std::vector<int> &match_vector) {
    int mask_length = pattern.length();
    unsigned long *buffer;
    unsigned long pattern_mask[255];
    int i, mismatch_iter;

    if (pattern[0] == '\0') {
        cout << "Empty pattern";
        return;
    }
    if (mask_length > 63) {
        cout << "Too long pattern";
        return;
    }

    /* Allocating memory for the match matrix and filling it with zeroes*/
    buffer = (unsigned long *) malloc((mismatch_limit + 1) * sizeof ( *buffer ));
    for (i=0; i <= mismatch_limit; ++i)
        buffer[i] = ~1;

    /* Filling patterns with ones */
    for (i=0; i <= 255; ++i)
        pattern_mask[i] = ~0;
    /* Building patterns */
    for (i=0; i < mask_length; ++i) {
        if (pattern[i] == 'N') {
            pattern_mask['A']  &= ~(1UL << i);
            pattern_mask['G']  &= ~(1UL << i);
            pattern_mask['C']  &= ~(1UL << i);
            pattern_mask['T']  &= ~(1UL << i);
        }
        if (pattern[i] == 'W') {
            pattern_mask['A']  &= ~(1UL << i);
            pattern_mask['T']  &= ~(1UL << i);
        }
        if (pattern[i] == 'S') {
            pattern_mask['C']  &= ~(1UL << i);
            pattern_mask['G']  &= ~(1UL << i);
        }
        if (pattern[i] == 'M') {
            pattern_mask['A']  &= ~(1UL << i);
            pattern_mask['C']  &= ~(1UL << i);
        }
        if (pattern[i] == 'K') {
            pattern_mask['G']  &= ~(1UL << i);
            pattern_mask['T']  &= ~(1UL << i);
        }
        if (pattern[i] == 'R') {
            pattern_mask['A']  &= ~(1UL << i);
            pattern_mask['G']  &= ~(1UL << i);
        }
        if (pattern[i] == 'Y') {
            pattern_mask['C']  &= ~(1UL << i);
            pattern_mask['T']  &= ~(1UL << i);
        }
        if (pattern[i] == 'B') {
            pattern_mask['G']  &= ~(1UL << i);
            pattern_mask['C']  &= ~(1UL << i);
            pattern_mask['T']  &= ~(1UL << i);
        }
        if (pattern[i] == 'D') {
            pattern_mask['A']  &= ~(1UL << i);
            pattern_mask['G']  &= ~(1UL << i);
            pattern_mask['T']  &= ~(1UL << i);
        }
        if (pattern[i] == 'H') {
            pattern_mask['A']  &= ~(1UL << i);
            pattern_mask['C']  &= ~(1UL << i);
            pattern_mask['T']  &= ~(1UL << i);
        }
        if (pattern[i] == 'V') {
            pattern_mask['A']  &= ~(1UL << i);
            pattern_mask['G']  &= ~(1UL << i);
            pattern_mask['C']  &= ~(1UL << i);
        }
        pattern_mask[pattern[i]] &= ~(1UL << i);
    }
    unsigned long old_buffer;
    for (i=0; text[i] != '\0'; ++i) {
        old_buffer = buffer[0];
        /* Making a step */
        buffer[0] |= pattern_mask[text[i]];
        buffer[0] <<= 1;
        /* Checking mask length + 1 column, counting mismatch number */
        for (mismatch_iter=1; mismatch_iter <= mismatch_limit; ++mismatch_iter) {
            unsigned long tmp = buffer[mismatch_iter];
            buffer[mismatch_iter] = (old_buffer & (buffer[mismatch_iter] | pattern_mask[text[i]])) << 1;
            old_buffer = tmp;
        }

        /* If there's a zero we have a match */
        if (0 == (buffer[mismatch_limit] & (1UL << mask_length))) {
            match_vector.push_back(i - mask_length + 1);
        }
    }

    /* Free memory for the buffer */
    free(buffer);
}

/* Unsuccessful print with buffer
void print1f_simple(FILE *filei, std::vector<string> read){
        std::string catstr = read[0]+read[1]+read[2]+read[3];
        if (curpoi+catstr.length() > buffersize){
                fwrite(buffer1, curpoi, 1, filei);
                curpoi = 0;
        }
                memcpy(buffer1 + curpoi, catstr.c_str(), catstr.length());
                curpoi+= catstr.length();
}
*/

// defining buffers for effiecnt writing
std::string dnabuffer = "";
std::string rnabuffer = "";
std::string nbbuffer = "";
std::string typesbuffer = "";

// function that writes reads to a file
void print1f_simple(int filei, std::vector<string> read, std::string &wbuffer) {
    if (!read.empty())
        wbuffer += read[0] + "\n" + read[1] + "\n" + read[2] + "\n" + read[3] + "\n" ;
    if (wbuffer.length() > 10000000 || wbuffer.length() != 0 && read.empty()) {
        //if (read.empty()) cout << "writing " << wbuffer.length() << " bytes..." << endl;
        //cout <<  wbuffer.length() << " <- length, isempty -> " << read.empty() << endl;
        //cout << filei << endl;
        write(filei, wbuffer.c_str(), wbuffer.length());
        wbuffer="";
    }
}

// function that writes statistics (string) to a file
void print1f_types(int filei, string line, std::string &wbuffer) {
    if (line.length() != 0)
        wbuffer += line;
    if (wbuffer.length() > 10000000 || wbuffer.length() != 0 && line.length() == 0) {
        //if (read.empty()) cout << "writing " << wbuffer.length() << " bytes..." << endl;
        write(filei, wbuffer.c_str(), wbuffer.length());
        wbuffer="";
    }
}



// adds reverse complement to an empty string
void rcomplement(string a, string &out)
{
    // i = 'ATGCUMRWSYKVHDBN'
    //o = 'TACGAKYWSRMBDHVN'
    std::map<char, char> conv = { { 'A', 'T'}, { 'G', 'C' }, { 'C', 'G' }, { 'T', 'A' }, { 'U', 'A' }, { 'M', 'K' }, { 'R', 'Y' }, { 'W', 'W' }, { 'S', 'S' }, { 'Y', 'R' }, { 'K', 'M' }, { 'V', 'B' }, { 'H', 'D' }, { 'D', 'H' }, { 'B', 'V' }, { 'N', 'N' } };
    for (int i = a.length()-1; i >= 0; i--) {
        out.push_back(conv[a[i]]);
    }
}


// parser for a description sequence, changes the value of success
void parser(int &success,  string describe_seq, string restofline, string restoflineq, std::vector<int> &match_vector) {
    // defining iterator for description sequnce, logical variables for states, and auxillary stings for operators
    int i;
    int cnt = 0;
    bool DNAu = false;
    bool RNAu = false;
    bool MMNu = false;
    bool bridgeu = false;
    bool bridgeb = true;
    string brseq = "";
    string Cut = "";
    string CheckDNA = "";
    string CheckRNA = "";
    string MMN = "";
    //  parsing description sequence
    for (i=0; i<=describe_seq.length(); i++) {
        // cout << "DSI: " << describe_seq[i] << endl;
        string allseq = restofline;
        string allq = restoflineq;
        // states
        if (describe_seq[i] == '*') {
            DNAu = true;
            RNAu = false;
            MMNu = false;
            bridgeu = false;
        }
        if (describe_seq[i] == '.') {
            DNAu = false;
            RNAu = true;
            MMNu = false;
            bridgeu = false;
        }
        if (describe_seq[i] == 'b') {
            bridgeu = true;
            MMNu = false;
        }
        if (describe_seq[i] == '(') {
            MMNu = true;
        }
        // operators
        // cout << "DSI: " << describe_seq[i] << " MMNu: " << MMNu << endl;
        if (bridgeu && describe_seq[i] != 'b' && !MMNu) {
            brseq += describe_seq[i];
        }
        if (MMNu && describe_seq[i] != '(') {
            MMN += describe_seq[i];
            //cout << describe_seq[i] << " dsi"  << endl;
        }
        // cout << describe_seq[i] << ": ";
        // cout << "bridgeu, DNAu, RNAu, desc_seq[i]: " << bridgeu << " " << DNAu << " " << RNAu << " " << describe_seq[i] << endl;
        // when ) is met, time to search
        if (describe_seq[i] == ')') {
            // cout << i << " MMN " << MMN << endl;
            // MMN = MMN.substr(0, MMN.length() - 1);
            // Cut branch
            // bridge cut branch
            if (bridgeu) {
                // cout << i << endl;
                // cnt += 1;
                // cout << Cut << " " << MMN << endl;
                // cout << "<<<<<<<<<<<<<<<<<<<<" << endl;
                //cout << Cut << " " << restofline << " " << stoi(MMN) << " " << endl;
                bitap_approximate_search(brseq, restofline, stoi(MMN), match_vector);
                //cout << match_vector.size() << endl;
                if (match_vector.size() == 0) {
                    cnt += 1;
                    // cout << "NF<" << endl;
                    success = 0;
                    brstart = -1;
                    brend = -1;
                    break;
                }
                brstart = match_vector[0];
                brend =  match_vector[0] + brseq.length();
                // cout << DNAu << endl;
                if (DNAu) {
                    // cout << dna_read[1] << " bf ";
                    dna_read[1] += restofline.substr(0, match_vector[0]);
                    dna_read[3] += restoflineq.substr(0, match_vector[0]);
                    // cout << dna_read[1] << " af " << endl;
                }
                if (RNAu) {
                    rna_read[1] += restofline.substr(0, match_vector[0]);
                    rna_read[3] += restoflineq.substr(0, match_vector[0]);
                }
                restofline = restofline.substr(match_vector[0] + brseq.length(), strlen(restofline.c_str())+1);
                restoflineq = restoflineq.substr(match_vector[0] + brseq.length(), strlen(restoflineq.c_str())+1);
                // cout << "restofline " << restofline << " i: " << i << endl;
                match_vector.clear();
                bridgeu = "";
            }
            // changing to a zero state
            DNAu = false;
            RNAu = false;
            MMNu = false;
            bridgeu = false;
            MMN = "";
        }
        // if no more operators then append parts
        if (describe_seq[i] == '\0' && success == 1) {
            if (DNAu) {
                dna_read[1] += restofline.substr(0, restofline.length());
                dna_read[3] += restoflineq.substr(0, restofline.length());
                DNAu = false;
            }
            if (RNAu) {
                rna_read[1] += restofline.substr(0, restofline.length());
                rna_read[3] += restoflineq.substr(0, restofline.length());
                RNAu = false;
            }
        }
    }
}

// function for spliting a line on a character
size_t split(const std::string &txt, std::vector<std::string> &strs, char ch)
{
    size_t pos = txt.find( ch );
    size_t initialPos = 0;
    strs.clear();

    // Decompose statement
    while( pos != std::string::npos ) {
        strs.push_back( txt.substr( initialPos, pos - initialPos ) );
        initialPos = pos + 1;

        pos = txt.find( ch, initialPos );
    }

    // Add the last one
    strs.push_back( txt.substr( initialPos, std::min( pos, txt.size() ) - initialPos + 1 ) );

    return strs.size();
}

// main
int main(int argc, char** argv)
{
    // options parser
    int SEu = 0;
    int PEu = 0;
    int TSVu = 0;
    int sepu = 0;
    int min_length_rna = 0;
    int min_length_dna = 0;
    string inputfilepath = "";
    string inputfilepath1 = "";
    string inputfilepath2 = "";
    string describe_seq = "";
    string bridge_codes = "10,01,20,02";
    int index;
    int c;
    string types0 = "";
    opterr = 0;
    while ((c = getopt (argc, argv, "hvseptu:i:k:j:d:l:m:")) != -1)
        switch (c) {
        case 'h':
            cout << "Tool fot processing reads in .FASTQ fomat" << endl;
            cout << "Parameters: " << endl;
            cout << "-h: help and description option" << endl;
            cout << "-i input file path, must be in .FASTQ format" << endl;
            cout << "-s single end reads mode, the ouput will be in 3 .FASTQ files - DNA- parts, RNA- parts and NB (no bridge) reads. !can be run with -p, combines both outputs" << endl;
            cout << "-p paired end mode, the output will be in .tsv file: read_id\t0|1|2(1 for forward bridge, 2 for reverse, 0 for no bridge)\tstart_position\tend_position. !can be run with -s, combines both outputs" << endl;
            cout << "-t .TSV statistic format" << endl;
            cout << "-u conditional codes for bridge strand; comma separated without whitespaces; default is 10,01,20,02 meaning one bridge of any strand, 10 would mean only main strand file 1. When SE mode 10=1 mode, 20=2 mode.";
            cout << "-d description sequence, must be silenced: * - start of DNA- part, . - start of RNA- part, < - start for cut part, ! - start for checked DNA- part. ? - start for checked RNA- part, ([0-9]+) - max mismatch limit, b - is the start for bridge (-p mode only). Example *<AGTC(1). will find a bridge AGTC and separate DNA- and RNA- parts." << endl;
            cout << "-l optional min length filter for a total DNA- and RNA- parts, default is 0." << endl;
            cout << "-u bridge codes for linker orientation; default 10,01,20,02 10 for paired end mode, 1,2 for single end; example forward linker in R1, 02 reverse complement linker in R2, for single end 1 - forward linker, 2 - reverse complement" << endl;
            cout << "Usage:" << endl;
            cout << "Single-end: ./alpha1 -s -i <path to the .fastq> -d <description sequence> -e -t -l <int; default 0> -m <int; default 0> -u <string; default 1,2>" << endl;
            cout << "Paired-ned: ./alpha1 -p -j <path to the R1 .fastq> -k <path the to R2 .fastq> -d <description sequence> -e -t -l <int; default 0> -m <int; default 0> -u <string; default 1,2>" << endl;
            break;
        case 'v':
            cout << "Version 1.3; to get more details on options and example code use -h" << endl;
            break;
        case 's':
            SEu = 1;
            break;
        case 'e':
            sepu = 1;
            break;
        case 'p':
            PEu = 1;
            break;
        case 't':
            TSVu = 1;
            break;
        case 'u':
            bridge_codes = optarg;
            break;
        case 'i':
            inputfilepath = optarg;
            break;
        case 'k':
            inputfilepath1 = optarg;
            break;
        case 'j':
            inputfilepath2 = optarg;
            break;
        case 'd':
            describe_seq = optarg;
            break;
        case 'l':
            min_length_rna = stoi(optarg);
            break;
        case 'm':
            min_length_dna = stoi(optarg);
            break;
            if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr,"Unknown option character `\\x%x'.\n", optopt);
            return 1;
        default:
            abort ();
        }
    //cout << "k: " << inputfilepath1 << "; j: " << inputfilepath2 << "; d: " << describe_seq << endl;
    //printf ("SEu = %d, PEu = %d, file2= %s, inputfilepath =%s, describe_seq= %s\n", SEu, PEu, inputfilepath2, inputfilepath1, describe_seq);
    //cout << inputfilepath << " <-in_file desc_seq-> " << describe_seq << endl;
    for (index = optind; index < argc; index++)
        printf ("Non-option argument %s\n", argv[index]);




    // bridge sequence, auxiliary line of reading, counter of string in one read, array of strings containting one read, vector of matches' positions
    // string describe_seq = argv[2];
    cout << "Description sequence: " << describe_seq << endl;
    string line;
    vector<int> match_vector_forward;
    vector<int> match_vector_reverse;
    vector<int> match_vector_forward1;
    vector<int> match_vector_reverse1;
    vector<int> match_vector_forward2;
    vector<int> match_vector_reverse2;
    short line_cnt_mod = 0;
    //int min_length = stoi(argv[3]);
    // int q = INT_MAX;
    int i;
    unsigned long long progress = 0;
    //cout << inputfilepath << endl;
    string s1 = "", s2 = "", s3 = "", s4 = "", s5 = "";
    // reading file and opening 4 output files
    split(bridge_codes, cond_codes, ',');
    if (SEu == 1) {
        //cout << "opening" << endl;
        s1 = inputfilepath;
        split(s1, name_split, '.');
        s1 = name_split[0] + ".DNA.fastq";
        s2 = name_split[0] + ".RNA.fastq";
        s3 = name_split[0] + ".garbage.fastq";
        s4 = name_split[0] + ".types.tsv";
        s5 = name_split[0] + ".codes.tsv";
    }
    if (PEu == 1) {
        s1 = inputfilepath1;
        split(s1, name_split, '.');
        s1 = name_split[0] + ".DNA.fastq";
        s2 = name_split[0] + ".RNA.fastq";
        s3 = name_split[0] + ".garbage.fastqlike";
        s4 = name_split[0] + ".types.tsv";
        s5 = name_split[0] + ".codes.tsv";
    }
    // FR nomenclature branch
    for(i=0; i<cond_codes.size(); i++) {
        if (SEu == 1) {
            if (cond_codes[i] == "F") condi_codes.insert("1");
            if (cond_codes[i] == "R") condi_codes.insert("2");
        }
        else if (PEu == 1) {
            if (cond_codes[i] == "F0") condi_codes.insert("10");
            if (cond_codes[i] == "R0") condi_codes.insert("20");
            if (cond_codes[i] == "0F") condi_codes.insert("01");
            if (cond_codes[i] == "0R") condi_codes.insert("02");
        }
    }
    // 10 nomenclature branch and default
    for(i=0; i<cond_codes.size(); i++) {
        if (SEu == 1) {
            if (cond_codes[i] == "10") condi_codes.insert("1");
            if (cond_codes[i] == "20") condi_codes.insert("2");
            if (cond_codes[i] == "1") condi_codes.insert("1");
            if (cond_codes[i] == "2") condi_codes.insert("2");
        }
        else if (PEu == 1 ||  cond_codes[i] == "10" ||  cond_codes[i] == "20" ||  cond_codes[i] == "01" ||  cond_codes[i] == "02" ) condi_codes.insert(cond_codes[i]);
    }
    /*cout << "condi_codes" << endl;
    for(auto it = condi_codes.begin(); it != condi_codes.end(); it++)
    {
    cout << *it << endl;
    }
    for(const auto& elem : condi_codes)
    {
    std::cout << elem  << "\n";
    }*/
    fstream infile;
    fstream infile1;
    fstream infile2;
    infile.open(inputfilepath, std::ifstream::in);
    //cout << infile.is_open() << endl;
    //cout << inputfilepath << endl;
    infile1.open(inputfilepath1, std::ios::in);
    infile2.open(inputfilepath2, std::ios::in);
    //cout << infile1.is_open() << endl;
    //cout << infile2.is_open() << endl;
    /*ofstream dnafile(s1);
    ofstream rnafile(s2);
    ofstream types(s4);

                        ofstream file_nb(s3);*/
    //ofstream types(s4);
    //FILE *types;

//      FILE *dnafile = fopen( s1.c_str() , "w" );
//      FILE *rnafile = fopen( s2.c_str() , "w" );
//      FILE *file_nb = fopen( s3.c_str() , "w" );
    //dnafile.open(s1);
    //rnafile.open(s2);



    // opening output files to !append
    const int dnafile = open(s1.c_str(), O_CREAT | O_WRONLY, 0644);
    const int rnafile = open(s2.c_str(), O_CREAT | O_WRONLY, 0644);
    const int file_nb = open(s3.c_str(), O_CREAT | O_WRONLY, 0644);
    const int types = open(s4.c_str(), O_CREAT | O_WRONLY, 0644);
    const int codes = open(s5.c_str(), O_CREAT | O_WRONLY, 0644);



    // checking if the program has rights
    if (access (s1.c_str(), W_OK) ||  access (s2.c_str(), W_OK) || access (s3.c_str(), W_OK) || access (s4.c_str(), W_OK)) {
        cerr << "Cannot write to output file" << endl;
        exit(1);
    }
    // reading input file
    if (SEu == 1) {
        //cout << "SEu" << endl;
        //cout << infile.is_open() << endl;
        if (infile.is_open())
        {
            //cout << "infile is open" << endl;
            while ( getline (infile,line) )
            {
                one_read[line_cnt_mod] = line;
                // checking whether read is fully read and if so searching for the bridge pattern
                // cout << line_cnt_mod << endl;
                if (line_cnt_mod == 3)  {
                    // progress counter
                    if (progress % 1000000 == 0 && progress != 0)
                        cout << "Processed " << progress / 1000000 << " million reads" << endl;
                    progress += 1;
                    string read_seq = one_read[1];
                    //cout << one_read[0] << endl;
                    string restofline = read_seq;
                    string restoflineq = one_read[3];
                    string restofline_rc = "";
                    string restoflineq_rc = "";
                    rcomplement(restofline, restofline_rc);
                    bool DNAu = false;
                    bool RNAu = false;
                    bool Cutu = false;
                    bool CheckDNAu = false;
                    bool CheckRNAu = false;
                    bool MMNu = false;
                    bool RNA_Length_Success = true;
                    bool DNA_Length_Success = true;
                    int successf = 1;
                    int successr = 1;
                    string Cut = "";
                    string CheckDNA = "";
                    string CheckRNA = "";
                    string MMN = "";
                    for(i=0; i<4; i++) {
                        dna_read[i] = "";
                    }
                    for(i=0; i<4; i++) {
                        rna_read[i] = "";
                    }
                    for(i=0; i<4; i++) {
                        read_nb[i] = "";
                    }
                    // spliting ID to SRR%d.%d
                    split( one_read[0], split_id, ' ' );
                    dna_read[0] = split_id[0];
                    dna_read[2] = "+\0";
                    rna_read[0] = split_id[0];
                    rna_read[2] = "+\0";
                    split_id.clear();
                    //cout << "succ, desseq, readseq, readseqq, matchvsize " << success << " " << describe_seq << " " << restofline << " " << restoflineq << " " << match_vector_forward.size() << endl;
                    // trying forward orientation of read
                    parser(successf, describe_seq,  restofline, restoflineq,  match_vector_forward);
                    RNA_Length_Success = rna_read[1].length() >= min_length_rna;
                    DNA_Length_Success = dna_read[1].length() >= min_length_dna;
                    if (!successf) {
                        dna_read[1] = "";
                        dna_read[3] = "";
                        rna_read[1] = "";
                        rna_read[3] = "";
                        split( one_read[0], split_id, ' ' );
                        dna_read[0] = split_id[0];
                        dna_read[2] = "+\0";
                        rna_read[0] = split_id[0];
                        rna_read[2] = "+\0";
                    }
                    else {
                        last_start = brstart;
                        last_end = brend;
                    }
                    for(i = restofline.length() - 1; i >= 0; i--)
                    {
                        restoflineq_rc += restoflineq[i];
                    }
                    parser(successr, describe_seq,  restofline_rc, restoflineq_rc,  match_vector_reverse);
                    if (successr) {
                        RNA_Length_Success = rna_read[1].length() >= min_length_rna;
                        DNA_Length_Success = dna_read[1].length() >= min_length_dna;
                    }
                    if (successf + successr == 1) {
                        if (successf == 1) {
                            if (DNA_Length_Success && RNA_Length_Success) {
                                int is_in = condi_codes.count("1");
                                if (sepu == 1 && is_in > 0) {
                                    print1f_simple(dnafile, dna_read, dnabuffer);
                                    print1f_simple(rnafile, rna_read, rnabuffer);
                                }
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tF-DR\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                    count_SE_codes["F-DR"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                            }
                            if (DNA_Length_Success && !RNA_Length_Success) {
                                split( one_read[0], split_id, ' ' );
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tF-D0\t" +  std::to_string(last_start + 1) + "\t" + std::to_string(last_end) +  "\n";
                                    count_SE_codes["F-D0"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                for(i=0; i<4; i++) {
                                    read_nb[i] = "";
                                }
                            }
                            if (!DNA_Length_Success && RNA_Length_Success) {
                                split( one_read[0], split_id, ' ' );
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tF-0R\t" +  std::to_string(last_start + 1) + "\t" + std::to_string(last_end) +  "\n";
                                    count_SE_codes["F-0R"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                for(i=0; i<4; i++) {
                                    read_nb[i] = "";
                                }
                            }
                            if (!DNA_Length_Success && !RNA_Length_Success) {
                                split( one_read[0], split_id, ' ' );
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tF-00\t" +  std::to_string(last_start + 1) + "\t" + std::to_string(last_end) +  "\n";
                                    count_SE_codes["F-00"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                for(i=0; i<4; i++) {
                                    read_nb[i] = "";
                                }
                            }
                        }
                        else {
                            if (DNA_Length_Success && RNA_Length_Success) {
                                int is_in = condi_codes.count("1");
                                if (sepu == 1 && is_in > 0) {
                                    print1f_simple(dnafile, dna_read, dnabuffer);
                                    print1f_simple(rnafile, rna_read, rnabuffer);
                                }
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tR-DR\t" + std::to_string(one_read[1].size() - brend + 1)  + "\t" +  std::to_string(one_read[1].size() - brstart) + "\n";
                                    count_SE_codes["R-DR"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                            }
                            if (DNA_Length_Success && !RNA_Length_Success) {
                                split( one_read[0], split_id, ' ' );
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tR-D0\t" + std::to_string(one_read[1].size() - brend + 1)  + "\t" +  std::to_string(one_read[1].size() - brstart) +  "\n";
                                    count_SE_codes["R-D0"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                for(i=0; i<4; i++) {
                                    read_nb[i] = "";
                                }
                            }
                            if (!DNA_Length_Success && RNA_Length_Success) {
                                split( one_read[0], split_id, ' ' );
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tR-0R\t" + std::to_string(one_read[1].size() - brend + 1)  + "\t" +  std::to_string(one_read[1].size() - brstart) +  "\n";
                                    count_SE_codes["R-0R"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                for(i=0; i<4; i++) {
                                    read_nb[i] = "";
                                }
                            }
                            if (!DNA_Length_Success && !RNA_Length_Success) {
                                split( one_read[0], split_id, ' ' );
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                                if (TSVu == 1) {
                                    types0 = dna_read[0] + "\tR-00\t" + std::to_string(one_read[1].size() - brend + 1)  + "\t" + std::to_string(one_read[1].size() - brstart) + "\n";
                                    count_SE_codes["R-00"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                for(i=0; i<4; i++) {
                                    read_nb[i] = "";
                                }
                            }
                        }
                    }
                    if (successf + successr == 2) {
                        split( one_read[0], split_id, ' ' );
                        read_nb[0] = split_id[0];
                        read_nb[1] = one_read[1];
                        read_nb[3] = one_read[3];
                        read_nb[2] = "+";
                        split_id.clear();
                        if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                        if (TSVu == 1) {
                            types0 = dna_read[0] +  "\tFR\t-1\t-1" +  "\n";
                            count_SE_codes["FR"] += 1;
                            brstart = 0;
                            brend = 0;
                            print1f_types(types, types0, typesbuffer);
                            types0 = "";
                        }
                        for(i=0; i<4; i++) {
                            read_nb[i] = "";
                        }
                    }
                    if (successf + successr == 0) {
                        split( one_read[0], split_id, ' ' );
                        read_nb[0] = split_id[0];
                        read_nb[1] = one_read[1];
                        read_nb[3] = one_read[3];
                        read_nb[2] = "+";
                        split_id.clear();
                        if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                        if (TSVu == 1) {
                            types0 = dna_read[0] +  "\t0\t-1\t-1" +  "\n";
                            count_SE_codes["0"] += 1;
                            brstart = 0;
                            brend = 0;
                            print1f_types(types, types0, typesbuffer);
                            types0 = "";
                        }
                        for(i=0; i<4; i++) {
                            read_nb[i] = "";
                        }
                    }

                    one_read = { "", "", "", "" };
                    line_cnt_mod = -1;
                }
                line_cnt_mod += 1;

            }
            match_vector_forward.clear();
            match_vector_reverse.clear();
            line_cnt_mod = 0;
            infile.close();
            // cout << cnt << endl;
        }
    }
    if (PEu==1) {
        string line1, line2;
        //cout << "PEu" << endl;
        //cout << infile1.is_open() << endl;
        //cout << infile2.is_open() << endl;
        while(getline(infile1, line1) && getline(infile2, line2)) {
            //cout << "lines1,2 " << line1 << " " << line2 << endl;
            first_read[line_cnt_mod] = line1;
            second_read[line_cnt_mod] = line2;
            //cout << line_cnt_mod << endl;
            if (line_cnt_mod == 3) {
                if (progress % 1000000 == 0 && progress != 0) cout << "Processed " << progress / 1000000 << " million reads" << endl;
                //cout << "PEu reads" << endl;
                split( first_read[0], split_id, ' ' );
                first_read[0] = split_id[0];
                split_id.clear();
                split( second_read[0], split_id, ' ' );
                second_read[0] = split_id[0];
                split_id.clear();
                if (first_read[0] != second_read[0]) {
                    cerr << first_read[0] << " " << second_read[0] << endl;
                    cerr << "input files are not sorted, different read IDs in the lines" << endl;
                    exit(1);
                }
                bool RNA_Length_Success = true;
                bool DNA_Length_Success = true;
                progress += 1;
                //cout << "first_read length" << first_read.size() << endl;
                string read_seq1 = first_read[1];
                string restofline1 = read_seq1;
                string restoflineq1 = first_read[3];
                string restofline_rc1 = "";
                string restoflineq_rc1 = "";
                //cout << "prerc" << endl;
                rcomplement(restofline1, restofline_rc1);
                string read_seq2 = second_read[1];
                string restofline2 = read_seq2;
                string restoflineq2 = second_read[3];
                string restofline_rc2 = "";
                string restoflineq_rc2 = "";
                rcomplement(restofline2, restofline_rc2);
                //cout << "postrc" << endl;
                for(i=0; i<4; i++) dna_read[i] = "";
                for(i=0; i<4; i++) rna_read[i] = "";
                for(i=0; i<4; i++) read_nb[i] = "";
                //cout << "presplit" << endl;
                // spliting ID to SRR%d.%d
                split( first_read[0], split_id, ' ' );
                dna_read[0] = split_id[0];
                dna_read[2] = "+\0";
                rna_read[0] = split_id[0];
                rna_read[2] = "+\0";
                //cout << "postsplit" << endl;
                int success1F = 1;
                int success1R = 1;
                int success2F = 1;
                int success2R = 1;
                for(i = restofline1.length() - 1; i >= 0; i--)
                    restoflineq_rc1 += restoflineq1[i];
                for(i = restofline2.length() - 1; i >= 0; i--)
                    restoflineq_rc2 += restoflineq2[i];
                //cout << "preparser" << endl;
                //cout << success1F << " " << describe_seq << " " << restofline1 << " " << restoflineq1 << " " << match_vector_forward1.size() << endl;
                parser(success1F, describe_seq,  restofline1, restoflineq1,  match_vector_forward1);
                match_vector_forward1.clear();
                if (success1F != 1) {
                    dna_read[1] = "";
                    dna_read[3] = "";
                    rna_read[1] = "";
                    rna_read[3] = "";
                }
                else {
                    for(i = 3; i >= 0; i--) {
                        dna_read_clear[i] = dna_read[i];
                        rna_read_clear[i] = rna_read[i];
                    }
                    last_start = brstart;
                    last_end = brend;
                }
                //cout << "parser1" << endl;
                //cout << success1R << " " << describe_seq << " " << restofline_rc1 << " " << restoflineq_rc1 << " " << match_vector_reverse1.size() << endl;
                parser(success1R, describe_seq,  restofline_rc1, restoflineq_rc1,  match_vector_reverse1);
                //cout << brstart << endl;
                match_vector_reverse1.clear();
                if (success1R != 1) {
                    dna_read[1] = "";
                    dna_read[3] = "";
                    rna_read[1] = "";
                    rna_read[3] = "";
                }
                else {
                    for(i = 3; i >= 0; i--) {
                        dna_read_clear[i] = dna_read[i];
                        rna_read_clear[i] = rna_read[i];
                    }
                    last_start = brstart;
                    last_end = brend;
                }
                //cout << "parser2" << endl;
                parser(success2F, describe_seq,  restofline2, restoflineq2,  match_vector_forward2);
                //cout << brstart << endl;
                match_vector_forward2.clear();
                if (success2F != 1) {
                    dna_read[1] = "";
                    dna_read[3] = "";
                    rna_read[1] = "";
                    rna_read[3] = "";
                }
                else {
                    for(i = 3; i >= 0; i--) {
                        dna_read_clear[i] = dna_read[i];
                        rna_read_clear[i] = rna_read[i];
                    }
                    last_start = brstart;
                    last_end = brend;
                }
                //cout << "parser3" << endl;
                parser(success2R, describe_seq,  restofline_rc2, restoflineq_rc2,  match_vector_reverse2);
                //cout << brstart << endl;
                match_vector_reverse2.clear();
                if (success2R != 1) {
                    dna_read[1] = "";
                    dna_read[3] = "";
                    rna_read[1] = "";
                    rna_read[3] = "";
                }
                else
                    for(i = 3; i >= 0; i--) {
                        dna_read_clear[i] = dna_read[i];
                        rna_read_clear[i] = rna_read[i];
                    }
                //cout << "postparser" << endl;
                //printf("%d, %d, %d, %d\n", success1F , success1R , success2F , success2R);
                if (success1F + success1R + success2F + success2R == 1) {
                    RNA_Length_Success = rna_read_clear[1].length() >= min_length_rna;
                    DNA_Length_Success = dna_read_clear[1].length() >= min_length_dna;
                    if (RNA_Length_Success && DNA_Length_Success) {
                        if (success1F == 1) {
                            int is_in = condi_codes.count("10");
                            if (sepu == 1 && is_in) {
                                print1f_simple(dnafile, dna_read_clear, dnabuffer);
                                print1f_simple(rnafile, rna_read_clear, rnabuffer);
                            }
                            if (TSVu == 1) {
                                types0 = dna_read[0] +  "\tF-0-DR\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                count_PE_codes["F-0-DR"] += 1;
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                        if (success1R == 1) {
                            int is_in = condi_codes.count("20");
                            if (sepu == 1 && is_in) {
                                print1f_simple(dnafile, dna_read_clear, dnabuffer);
                                print1f_simple(rnafile, rna_read_clear, rnabuffer);
                            }
                            if (TSVu == 1) {
                                types0 = dna_read[0] +  "\tR-0-DR\t" + std::to_string(first_read[1].size() - last_end + 1) + "\t" + std::to_string(first_read[1].size() - last_start) + "\n";
                                count_PE_codes["R-0-DR"] += 1;
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                        if (success2F == 1) {
                            int is_in = condi_codes.count("01");
                            if (sepu == 1 && is_in) {
                                print1f_simple(dnafile, dna_read_clear, dnabuffer);
                                print1f_simple(rnafile, rna_read_clear, rnabuffer);
                            }
                            if (TSVu == 1) {
                                types0 = dna_read[0] +  "\t0-F-DR\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                count_PE_codes["0-F-DR"] += 1;
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                        // if (success2R == 1)
                        //      cout << rna_read[1].length() << " " << dna_read[1].length() << endl;
                        if (success2R == 1) {
                            int is_in = condi_codes.count("02");
                            //printf("success2R 1");
                            if (sepu == 1 && is_in) {
                                print1f_simple(dnafile, dna_read_clear, dnabuffer);
                                print1f_simple(rnafile, rna_read_clear, rnabuffer);
                            }
                            if (TSVu == 1) {
                                types0 = dna_read[0] +  "\t0-R-DR\t" + std::to_string(first_read[1].size() - brend + 1) + "\t" + std::to_string(first_read[1].size() - brstart) + "\n";
                                count_PE_codes["0-R-DR"] += 1;
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                    }
                    else {
                        if (success1F == 1) {
                            if (sepu == 1) {
                                read_nb[0] = split_id[0];
                                read_nb[1] = read_seq1 + "|" + read_seq2;
                                read_nb[3] = first_read[3] + "|" + second_read[3];
                                read_nb[2] = "+";
                                print1f_simple(file_nb, read_nb, nbbuffer);
                                split_id.clear();
                            }
                            if (TSVu == 1) {
                                if (DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tF-0-D0"  + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                    count_PE_codes["F-0-D0"] += 1;
                                }
                                if (!DNA_Length_Success && RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tF-0-0R" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                    count_PE_codes["F-0-0R"] += 1;
                                }
                                if (!DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tF-0-00"  + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                    count_PE_codes["F-0-00"] += 1;
                                }
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                        if (success1R == 1) {
                            if (sepu == 1) {
                                read_nb[0] = split_id[0];
                                read_nb[1] = read_seq1 + "|" + read_seq2;
                                read_nb[3] = first_read[3] + "|" + second_read[3];
                                read_nb[2] = "+";
                                print1f_simple(file_nb, read_nb, nbbuffer);
                                split_id.clear();
                            }
                            if (TSVu == 1) {
                                if (DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tR-0-D0\t" + std::to_string(first_read[1].size() - last_end + 1) + "\t" + std::to_string(first_read[1].size() - last_start) + "\n";
                                    count_PE_codes["R-0-D0"] += 1;
                                }
                                if (!DNA_Length_Success && RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tR-0-0R\t" + std::to_string(first_read[1].size() - last_end + 1) + "\t" + std::to_string(first_read[1].size() - last_start) + "\n";
                                    count_PE_codes["R-0-0R"] += 1;
                                }
                                if (!DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tR-0-00\t" + std::to_string(first_read[1].size() - last_end + 1) + "\t" + std::to_string(first_read[1].size() - last_start) + "\n";
                                    count_PE_codes["R-0-00"] += 1;
                                }
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                        if (success2F == 1) {
                            if (sepu == 1) {
                                read_nb[0] = split_id[0];
                                read_nb[1] = read_seq1 + "|" + read_seq2;
                                read_nb[3] = first_read[3] + "|" + second_read[3];
                                read_nb[2] = "+";
                                print1f_simple(file_nb, read_nb, nbbuffer);
                                split_id.clear();

                            }
                            if (TSVu == 1) {
                                if (DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-F-D0\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                    count_PE_codes["0-F-D0"] += 1;
                                }
                                if (!DNA_Length_Success && RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-F-0R\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                    count_PE_codes["0-F-0R"] += 1;
                                }
                                if (!DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-F-00\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\n";
                                    count_PE_codes["0-F-00"] += 1;
                                }
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                        // if (success2R == 1)
                        //  cout << rna_read[1].length() << " " << dna_read[1].length() << endl;
                        if (success2R == 1) {
                            //printf("success2R 1");
                            if (sepu == 1) {
                                read_nb[0] = split_id[0];
                                read_nb[1] = read_seq1 + "|" + read_seq2;
                                read_nb[3] = first_read[3] + "|" + second_read[3];
                                read_nb[2] = "+";
                                print1f_simple(file_nb, read_nb, nbbuffer);
                                split_id.clear();
                            }
                            if (TSVu == 1) {
                                if (DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-R-D0\t" + std::to_string(first_read[1].size() - brend + 1) + "\t" + std::to_string(first_read[1].size() - brstart) + "\n";
                                    count_PE_codes["0-R-D0"] += 1;
                                }
                                if (!DNA_Length_Success && RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-R-0R\t" + std::to_string(first_read[1].size() - brend + 1) + "\t" + std::to_string(first_read[1].size() - brstart) + "\n";
                                    count_PE_codes["0-R-0R"] += 1;
                                }
                                if (!DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-R-00\t" + std::to_string(first_read[1].size() - brend + 1) + "\t" + std::to_string(first_read[1].size() - brstart) + "\n";
                                    count_PE_codes["0-R-00"] += 1;
                                }
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                    }
                }
                else {
                    if (sepu == 1) {
                        read_nb[0] = split_id[0];
                        read_nb[1] = read_seq1 + "|" + read_seq2;
                        read_nb[3] = first_read[3] + "|" + second_read[3];
                        read_nb[2] = "+";
                        print1f_simple(file_nb, read_nb, nbbuffer);
                        split_id.clear();
                    }
                    if (TSVu == 1) {
                        if (success1F + success1R + success2F + success2R == 0) {
                            types0 = read_nb[0] +  "\t0-0\t-1\t-1" + "\n";
                            count_PE_codes["0-0"] += 1;
                        }
                        if (success1F && !success1R && success2F && !success2R) {
                            types0 = read_nb[0] +  "\tF-F\t-1\t-1" + "\n";
                            count_PE_codes["F-F"] += 1;
                        }
                        if (!success1F && success1R && !success2F && success2R) {
                            types0 = read_nb[0] +  "\tR-R\t-1\t-1" + "\n";
                            count_PE_codes["R-R"] += 1;
                        }
                        if (success1F && !success1R && !success2F && success2R) {
                            types0 = read_nb[0] +  "\tF-R\t-1\t-1" + "\n";
                            count_PE_codes["F-R"] += 1;
                        }
                        if (!success1F && success1R && success2F && !success2R) {
                            types0 = read_nb[0] +  "\tR-F\t-1\t-1" + "\n";
                            count_PE_codes["R-F"] += 1;
                        }
                        if (!success1F && !success1R && success2F && success2R) {
                            types0 = read_nb[0] +  "\t0-FR\t-1\t-1" + "\n";
                            count_PE_codes["0-FR"] += 1;
                        }
                        if (success1F && success1R && !success2F && !success2R) {
                            types0 = read_nb[0] +  "\tFR-0\t-1\t-1" + "\n";
                            count_PE_codes["FR-0"] += 1;
                        }
                        if (!success1F && success1R && success2F && success2R) {
                            types0 = read_nb[0] +  "\tR-FR\t-1\t-1" + "\n";
                            count_PE_codes["R-FR"] += 1;
                        }
                        if (success1F && !success1R && success2F && success2R) {
                            types0 = read_nb[0] +  "\tF-FR\t-1\t-1" + "\n";
                            count_PE_codes["F-FR"] += 1;
                        }
                        if (success1F && success1R && !success2F && success2R) {
                            types0 = read_nb[0] +  "\tFR-R\t-1\t-1" + "\n";
                            count_PE_codes["FR-R"] += 1;
                        }
                        if (success1F && success1R && success2F && !success2R) {
                            types0 = read_nb[0] +  "\tFR-F\t-1\t-1" + "\n";
                            count_PE_codes["FR-F"] += 1;
                        }
                        if (success1F && success1R && success2F && success2R) {
                            types0 = read_nb[0] +  "\tFR-FR\t-1\t-1" + "\n";
                            count_PE_codes["FR-FR"] += 1;
                        }

                        brstart = 0;
                        brend = 0;
                        last_start = 0;
                        last_end = 0;
                        print1f_types(types, types0, typesbuffer);
                        types0 = "";
                    }
                }
                first_read = { "", "", "", "" };
                second_read = { "", "", "", "" };
                line_cnt_mod = -1;

            }
            line_cnt_mod += 1;
        }
        match_vector_forward1.clear();
        match_vector_reverse1.clear();
        match_vector_forward2.clear();
        match_vector_reverse2.clear();
        line_cnt_mod = 0;
        infile1.close();
        infile2.close();
    }

    if (SEu == 0 && PEu == 0) cout << "No mode specified";
    /*dnafile.close();
                rnafile.close();
    types.close();
    file_nb();*/
    //fclose(dnafile);
    //fclose(rnafile);
    //fclose(file_nb);
    // printing if buffer is not empty in the end
    print1f_simple(dnafile, empty_vec, dnabuffer);
    print1f_simple(rnafile, empty_vec, rnabuffer);
    print1f_simple(file_nb, empty_vec, nbbuffer);
    print1f_types(types, types0, typesbuffer);
    close(dnafile);
    close(rnafile);
    close(file_nb);
    close(types);
    string code_out;
    if (SEu == 1)
        for (auto element :count_SE_codes)  {
            code_out = element.first + "\t" + to_string(element.second) + "\n";
            write(codes, code_out.c_str(), code_out.length());
        }
    if (PEu == 1)
        for (auto element :count_PE_codes)  {
            code_out = element.first + "\t" + to_string(element.second) + "\n";
            write(codes, code_out.c_str(), code_out.length());
        }
    close(codes);
    //fclose(types);
    //types.close();

    return 0;
}
