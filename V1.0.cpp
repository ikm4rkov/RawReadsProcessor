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
std::vector<std::string> dna_read = { "", "", "", "" };
std::vector<std::string> rna_read = { "", "", "", "" };
std::vector<std::string> read_nb = { "", "", "", "" };
std::vector<std::string> empty_vec = {};
std::vector<std::string> split_id = {};

// that way with a buffer was unsuccessful
/*int buffersize = 10*1024*1024;
char *buffer1 = (char*) malloc(buffersize*sizeof(char));
int curpoi = 0;*/

// global vairables fro bridge coordinates
int brstart = 0;
int brend = 0;

// search function, returns nothing but can change values
void bitap_approximate_search(string pattern, string text, int mismatch_limit, std::vector<int> &match_vector){
    int mask_length = pattern.length();
    unsigned long *buffer;
    unsigned long pattern_mask[255];
    int i, mismatch_iter;

    if (pattern[0] == '\0'){
        cout << "Empty pattern";
        return;
    }
    if (mask_length > 63){
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
    for (i=0; i < mask_length; ++i)
         pattern_mask[pattern[i]] &= ~(1UL << i);
    unsigned long old_buffer;
    for (i=0; text[i] != '\0'; ++i){
        old_buffer = buffer[0];
    /* Making a step */
                                buffer[0] |= pattern_mask[text[i]];
                                buffer[0] <<= 1;
        /* Checking mask length + 1 column, counting mismatch number */
        for (mismatch_iter=1; mismatch_iter <= mismatch_limit; ++mismatch_iter){
            unsigned long tmp = buffer[mismatch_iter];
            buffer[mismatch_iter] = (old_buffer & (buffer[mismatch_iter] | pattern_mask[text[i]])) << 1;
            old_buffer = tmp;
        }

        /* If there's a zero we have a match */
        if (0 == (buffer[mismatch_limit] & (1UL << mask_length))){
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
                std::map<char, char> conv = { { 'A' , 'T'}, { 'G' , 'C' }, { 'C' , 'G' }, { 'T' , 'A' }, { 'U' , 'A' }, { 'M' , 'K' }, { 'R' , 'Y' }, { 'W' , 'W' }, { 'S' , 'S' }, { 'Y' , 'R' }, { 'K' , 'M' }, { 'V' , 'B' }, { 'H' , 'D' }, { 'D' , 'H' }, { 'B' , 'V' }, { 'N' , 'N' } };
                for (int i = a.length()-1; i >= 0; i--) {
                    out.push_back(conv[a[i]]);
                }
}


// parser for a description sequence, changes the value of success
void parser(bool &success, string describe_seq, string restofline, string restoflineq, std::vector<int> &match_vector){
                // defining iterator for description sequnce, logical variables for states, and auxillary stings for operators
                int i;
                int cnt = 0;
                bool DNAu = false;
                bool RNAu = false;
                bool Cutu = false;
                bool CheckDNAu = false;
                bool CheckRNAu = false;
                bool MMNu = false;
                bool bridgeu = false;
                string brseq = "";
                string Cut = "";
                string CheckDNA = "";
                string CheckRNA = "";
                string MMN = "";
                //  parsing description sequence
                for (i=0; i<=describe_seq.length(); i++){
                                // cout << "DSI: " << describe_seq[i] << endl;
                                string allseq = restofline;
                                string allq = restoflineq;
                                // states
                                if (describe_seq[i] == '*'){
                                                DNAu = true;
                                                RNAu = false;
                                                Cutu = false;
                                                CheckDNAu = false;
                                                CheckRNAu = false;
                                                MMNu = false;
                                                bridgeu = false;
                                }
                                if (describe_seq[i] == '.'){
                                                DNAu = false;
                                                RNAu = true;
                                                Cutu = false;
                                                MMNu = false;
                                                bridgeu = false;
                                }
                                if (describe_seq[i] == 'b'){
                                                bridgeu = true;
                                                Cutu = false;
                                                CheckDNAu = false;
                                                CheckRNAu = false;
                                                MMNu = false;
                                }
                                if (describe_seq[i] == '<'){
                                                Cutu = true;
                                                bridgeu = false;
                                                CheckDNAu = false;
                                                CheckRNAu = false;
                                                MMNu = false;
                                }
                                if (describe_seq[i] == '?'){
                                                Cutu = false;
                                                bridgeu = false;
                                                CheckDNAu = true;
                                                CheckRNAu = false;
                                                MMNu = false;
                                }
                                if (describe_seq[i] == '!'){
                                                Cutu = false;
                                                bridgeu = false;
                                                CheckDNAu = false;
                                                CheckRNAu = true;
                                                MMNu = false;
                                }
                                if (describe_seq[i] == '('){
                                                MMNu = true;
                                }
                                // operators
                                // cout << "DSI: " << describe_seq[i] << " MMNu: " << MMNu << endl;
                                if (Cutu && describe_seq[i] != '<' && !MMNu){
                                                Cut += describe_seq[i];
                                }
                                if (bridgeu && describe_seq[i] != 'b' && !MMNu){
                                                brseq += describe_seq[i];
                                }
                                if (CheckDNAu && describe_seq[i] != '?' && !MMNu){
                                                CheckDNA += describe_seq[i];
                                }
                                if (CheckRNAu && describe_seq[i] != '!' && !MMNu){
                                                CheckRNA += describe_seq[i];
                                }
                                if (MMNu && describe_seq[i] != '('){
                                                MMN += describe_seq[i];
                                                //cout << describe_seq[i] << " dsi"  << endl;
                                }
                                // cout << describe_seq[i] << ": ";
                                // cout << "bridgeu, DNAu, RNAu, desc_seq[i]: " << bridgeu << " " << DNAu << " " << RNAu << " " << describe_seq[i] << endl;
                                // when ) is met, time to search
                                if (describe_seq[i] == ')'){
                                                // cout << i << " MMN " << MMN << endl;
                                                // MMN = MMN.substr(0, MMN.length() - 1);
                                                // Cut branch
                                                if (Cutu){
                                                                // cout << i << endl;
                                                                // cnt += 1;
                                                                // cout << Cut << " " << MMN << endl;
                                                                // cout << "<<<<<<<<<<<<<<<<<<<<" << endl;
                                                                //cout << Cut << " " << restofline << " " << stoi(MMN) << " " << endl;
                                                                bitap_approximate_search(Cut, restofline, stoi(MMN), match_vector);
                                                                // if not found break
                                                                if (match_vector.size() == 0){
                                                                                cnt += 1;
                                                                                // cout << "NF<" << endl;
                                                                                success = false;
                                                                                brstart = -1;
                                                                                brend = -1;
                                                                                break;
                                                                }
                                                                // filling DNA- and RNA- parts, cutting read sequence, chaning state and clearing match vector
                                                                if (DNAu){
                                                                                // cout << dna_read[1] << " bf ";
                                                                                dna_read[1] += restofline.substr(0, match_vector[0]);
                                                                                dna_read[3] += restoflineq.substr(0, match_vector[0]);
                                                                                // cout << dna_read[1] << " af " << endl;
                                                                }
                                                                if (RNAu){
                                                                                rna_read[1] += restofline.substr(0, match_vector[0]);
                                                                                rna_read[3] += restoflineq.substr(0, match_vector[0]);
                                                                }
                                                                                restofline = restofline.substr(match_vector[0] + Cut.length(), strlen(restofline.c_str())+1);
                                                                                restoflineq = restoflineq.substr(match_vector[0] + Cut.length(), strlen(restoflineq.c_str())+1);
                                                                                // cout << "restofline " << restofline << " i: " << i << endl;
                                                                                match_vector.clear();
                                                                Cut = "";
                                                }
                                                // bridge cut branch
                                                 if (bridgeu){
                                                                // cout << i << endl;
                                                                // cnt += 1;
                                                                // cout << Cut << " " << MMN << endl;
                                                                // cout << "<<<<<<<<<<<<<<<<<<<<" << endl;
                                                                //cout << Cut << " " << restofline << " " << stoi(MMN) << " " << endl;
                                                                bitap_approximate_search(brseq, restofline, stoi(MMN), match_vector);
                                                                if (match_vector.size() == 0){
                                                                                cnt += 1;
                                                                                // cout << "NF<" << endl;
                                                                                success = false;
                                                                                brstart = -1;
                                                                                brend = -1;
                                                                                break;
                                                                }
                                                                brstart = match_vector[0];
                                                                brend =  match_vector[0] + brseq.length();
                                                                // cout << DNAu << endl;
                                                                if (DNAu){
                                                                                // cout << dna_read[1] << " bf ";
                                                                                dna_read[1] += restofline.substr(0, match_vector[0]);
                                                                                dna_read[3] += restoflineq.substr(0, match_vector[0]);
                                                                                // cout << dna_read[1] << " af " << endl;
                                                                }
                                                                if (RNAu){
                                                                                rna_read[1] += restofline.substr(0, match_vector[0]);
                                                                                rna_read[3] += restoflineq.substr(0, match_vector[0]);
                                                                }
                                                                                restofline = restofline.substr(match_vector[0] + brseq.length(), strlen(restofline.c_str())+1);
                                                                                restoflineq = restoflineq.substr(match_vector[0] + brseq.length(), strlen(restoflineq.c_str())+1);
                                                                                // cout << "restofline " << restofline << " i: " << i << endl;
                                                                                match_vector.clear();
                                                                bridgeu = "";
                                                }
                                                // checking DNA
                                                if (CheckDNAu){
                                                                // cout << CheckDNA << endl;
                                                                // cout << "?????????????????" << endl;
                                                                // cout << CheckDNA << " " << restofline << " " << stoi(MMN) << " " << endl;
                                                                bitap_approximate_search(CheckDNA, restofline, stoi(MMN), match_vector);
                                                                if (match_vector.size() == 0){
                                                                                cnt += 1;
                                                                                success = false;
                                                                                // cout << "Not found ?" << endl;
                                                                                break;
                                                                }
                                                                // cout << restofline.substr(0, match_vector_forward[0]) << endl;
                                                                if (DNAu){
                                                                                dna_read[1] += restofline.substr(0, match_vector[0]);
                                                                                dna_read[3] += restoflineq.substr(0, match_vector[0]);
                                                                }
                                                                if (RNAu){
                                                                                rna_read[1] += restofline.substr(0, match_vector[0]);
                                                                                rna_read[3] += restoflineq.substr(0, match_vector[0]);
                                                                }
                                                                dna_read[1] += restofline.substr(match_vector[0], CheckDNA.length());
                                                                dna_read[3] += restoflineq.substr(match_vector[0], CheckDNA.length());
                                                                restofline = restofline.substr(match_vector[0] + CheckDNA.length(), strlen(restofline.c_str())+1);
                                                                restoflineq = restoflineq.substr(match_vector[0] + CheckDNA.length(), strlen(restoflineq.c_str())+1);
                                                                match_vector.clear();
                                                                CheckDNA = "";
                                                }
                                                // checking RNA
                                                if (CheckRNAu){
                                                                // cout << CheckRNA << endl;
                                                                // cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                                                                // cout << CheckRNA << " " << restofline << " " << stoi(MMN) << " " << endl;
                                                                bitap_approximate_search(CheckRNA, restofline, stoi(MMN), match_vector);
                                                                if (match_vector.size() == 0){
                                                                                cnt += 1;
                                                                                success = false;
                                                                                // cout << "Not found !" << endl;
                                                                                break;
                                                                }
                                                                if (DNAu){
                                                                                dna_read[1] += restofline.substr(0, match_vector[0]);
                                                                                dna_read[3] += restoflineq.substr(0, match_vector[0]);
                                                                }
                                                                if (RNAu){
                                                                                rna_read[1] += restofline.substr(0, match_vector[0]);
                                                                                rna_read[3] += restoflineq.substr(0, match_vector[0]);
                                                                }
                                                                // cout << "rread_1: " << rna_read[1] << endl;
                                                                // cout << restofline.substr(match_vector_forward[0], match_vector_forward[0] + CheckRNA.length()) << " " << match_vector_forward[0] << " " << CheckRNA.length() << endl;
                                                                rna_read[1] += restofline.substr(match_vector[0], CheckRNA.length());
                                                                rna_read[3] += restoflineq.substr(match_vector[0], CheckRNA.length());
                                                                // cout << "rread_2: " << rna_read[1] << endl;
                                                                restofline = restofline.substr(match_vector[0] + CheckRNA.length(), strlen(restofline.c_str())+1);
                                                                restoflineq = restoflineq.substr(match_vector[0] + CheckRNA.length(), strlen(restoflineq.c_str())+1);
                                                                match_vector.clear();
                                                                CheckRNA = "";
                                                }
                                // changing to a zero state
                                DNAu = false;
                                RNAu = false;
                                Cutu = false;
                                CheckDNAu = false;
                                CheckRNAu = false;
                                MMNu = false;
                                bridgeu = false;
                                MMN = "";
                                }
                                // if no more operators then append parts
                                if (describe_seq[i] == '\0' && success){
                                                if (DNAu){
                                                                dna_read[1] += restofline.substr(0, restofline.length());
                                                                dna_read[3] += restoflineq.substr(0, restofline.length());
                                                                DNAu = false;
                                                }
                                                if (RNAu){
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
        int min_length = 0;
        string inputfilepath = "";
        string describe_seq = "";
        int index;
        int c;
        string types0 = "";
        opterr = 0;
        while ((c = getopt (argc, argv, "hspi:d:l:")) != -1)
                switch (c){
                        case 'h':
                                cout << "Tool fot processing reads in .FASTQ fomat" << endl;
                                cout << "Parameters: " << endl;
                                cout << "-h: help and description option" << endl;
                                cout << "-i input file path, must be in .FASTQ format" << endl;
                                cout << "-s single end reads mode, the ouput will be in 3 .FASTQ files - DNA- parts, RNA- parts and NB (no bridge) reads. !can be run with -p, combines both outputs" << endl;
                                cout << "-p paired end mode, the output will be in .tsv file: read_id\t0|1|2(1 for forward bridge, 2 for reverse, 0 for no bridge)\tstart_position\tend_position. !can be run with -s, combines both outputs" << endl;
                                cout << "-d description sequence, must be silenced: * - start of DNA- part, . - start of RNA- part, < - start for cut part, ! - start for checked DNA- part. ? - start for checked RNA- part, ([0-9]+) - max mismatch limit, b - is the start for bridge (-p mode only). Example *<AGTC(1). will find a bridge AGTC and separate DNA- and RNA- parts." << endl;
                                cout << "-l optional min length filter for a total DNA- and RNA- parts, default is 0." << endl;
                                break;
                        case 's':
                                SEu = 1;
                                break;
                        case 'p':
                                PEu = 1;
                                break;
                        case 'i':
                                inputfilepath = optarg;
                                break;
                        case 'd':
                                describe_seq = optarg;
                                break;
                        case 'l':
                                min_length = stoi(optarg);
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
        //printf ("SEu = %d, PEu = %d, min_length= %d, inputfilepath =%s, describe_seq= %s\n", SEu, PEu, min_length, inputfilepath, describe_seq);
        //cout << inputfilepath << " <-in_file desc_seq-> " << describe_seq << endl;
        for (index = optind; index < argc; index++)
                printf ("Non-option argument %s\n", argv[index]);




        // bridge sequence, auxiliary line of reading, counter of string in one read, array of strings containting one read, vector of matches' positions
       // string describe_seq = argv[2];
        cout << "Description sequence: " << describe_seq << endl;
        string line;
        short line_cnt_mod = 0;
        //int min_length = stoi(argv[3]);
        // int q = INT_MAX;
                                int i;
                                unsigned long long progress = 0;
        vector<int> match_vector_forward;
                                vector<int> match_vector_reverse;
        // reading file and opening 4 output files
        ifstream infile(inputfilepath);
        string s1 = inputfilepath; s1.append("_DNA");
        string s2 = inputfilepath; s2.append("_RNA");
                                string s3 = inputfilepath; s3.append("_NB");
        string s4 = inputfilepath; s4.append("_types.tsv");
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
        const int dnafile = open(s1.c_str(), O_CREAT | O_WRONLY | O_DIRECT, 0644);
        const int rnafile = open(s2.c_str(), O_CREAT | O_WRONLY | O_DIRECT, 0644);
        const int file_nb = open(s3.c_str(), O_CREAT | O_WRONLY | O_DIRECT, 0644);
        const int types = open(s4.c_str(), O_CREAT | O_WRONLY | O_DIRECT, 0644);



        // checking if the program has rights
        if (access (s1.c_str(), W_OK) ||  access (s2.c_str(), W_OK)){
                                cerr << "Cannot write to output file" << endl;
                                exit(1);
                                }
        // reading input file
        if (infile.is_open())
        {
              while ( getline (infile,line) )
              {
              // checking whether read is fully read and if so searching for the bridge pattern
                                                        // cout << line << endl;
                                                                        if ((line[0] == '@') && (line_cnt_mod != 0) && (one_read[line_cnt_mod - 1] != "+"))
                                                        {
                                                                // progress counter
                                                                if (progress % 1000000 == 0 && progress != 0)
                                                                       cout << "Processed " << progress / 1000000 << " million reads" << endl;
                                                                                progress += 1;
                                                                                string read_seq = one_read[1];
                                                                                // cout << one_read[0] << endl;
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
                                                                                bool success = true;
                                                                                string Cut = "";
                                                                                string CheckDNA = "";
                                                                                string CheckRNA = "";
                                                                                string MMN = "";
                                                                                for(i=0;i<4;i++){
                        dna_read[i] = "";
                    }
                                                                                for(i=0;i<4;i++){
                        rna_read[i] = "";
                    }
                                                                                for(i=0;i<4;i++){
                        read_nb[i] = "";
                    }
                                                                                // spliting ID to SRR%d.%d
                                                                                split( one_read[0], split_id, ' ' );
                                                                                dna_read[0] = split_id[0];
                                                                                dna_read[2] = "+\0";
                                                                                rna_read[0] = split_id[0];
                                                                                rna_read[2] = "+\0";
                                                                                split_id.clear();
                                                                                // cout << "DNARBF " << dna_read[1] << endl;
                                                                                // cout << "succ, desseq, readseq, readseqq, matchvsize " << success << " " << describe_seq << " " << restofline << " " << restoflineq << " " << match_vector_forward.size() << endl;
                                                                                // trying forward orientation of read
                                                                                parser(success, describe_seq,  restofline, restoflineq,  match_vector_forward);
                                                                                // cout << "DNARAF " << dna_read[1] << endl;
                                                                                // cout << "success1 " << success << endl;
                                                                                if (success && dna_read[1].length() >= min_length && rna_read[1].length() >= min_length){
                                                                                        //      cout << "f1" << line << endl;
                                                                                                // cout << "seq " << rna_read[1] << endl;
                                                                                                // writing depending on mode
                                                                                                if (SEu == 1){
                                                                                                        print1f_simple(dnafile, dna_read, dnabuffer);
                                                                                                        print1f_simple(rnafile, rna_read, rnabuffer);
                                                                                                }
                                                                                                if (PEu == 1){
                                                                                                        types0 = dna_read[0] +  "\t1\t" + std::to_string(brstart) + "\t" + std::to_string(brend) + "\n";
                                                                                                        brstart = 0;
                                                                                                        brend = 0;
                                                                                                        print1f_types(types, types0, typesbuffer);
                                                                                                        types0 = "";
                                                                                                }
                                                                                }
                                                                                else{
                                                                                                // same for the reverse orientaion of read
                                                                                                for(i = restofline.length() - 1; i >= 0; i--)
                                                                                                {
                                                                                                                restoflineq_rc += restoflineq[i];
                                                                                                }
                                                                                                success = true;
                                                                                                dna_read[1] = "";
                                                                                                dna_read[3] = "";
                                                                                                rna_read[1] = "";
                                                                                                rna_read[3] = "";
                                                                                                split( one_read[0], split_id, ' ' );
                                                                                                dna_read[0] = split_id[0];
                                                                                                dna_read[2] = "+\0";
                                                                                                rna_read[0] = split_id[0];
                                                                                                rna_read[2] = "+\0";
                                                                                                split_id.clear();
                                                                                                parser(success, describe_seq,  restofline_rc, restoflineq_rc,  match_vector_reverse);
                                                                                                // cout << "success2 " << success << endl;
                                                                                                if (success && dna_read[1].length() >= min_length && rna_read[1].length() >= min_length){
                                                                                                                // cout << "f2" << endl;
                                                                                                                if (SEu == 1){
                                                                                                                        print1f_simple(dnafile, dna_read, dnabuffer);
                                                                                                                        print1f_simple(rnafile, rna_read, rnabuffer);
                                                                                                                }
                                                                                                                if (PEu == 1){
                                                                                                                        types0 = dna_read[0] +  "\t2\t" + std::to_string(brstart) + "\t" + std::to_string(brend) + "\n";
                                                                                                                        brstart = 0;
                                                                                                                        brend = 0;
                                                                                                                        print1f_types(types, types0, typesbuffer);
                                                                                                                        types0 = "";
                                                                                                                }



                                                                                                }
                                                                                                else{
                                                                                                                // writing if any condition is not true
                                                                                                                split( one_read[0], split_id, ' ' );
                                                                                                                read_nb[0] = split_id[0];
                                                                                                                read_nb[1] = read_seq;
                                                                                                                read_nb[3] = one_read[3];
                                                                                                                read_nb[2] = "+";
                                                                                                                split_id.clear();
                                                                                                                if (SEu == 1) print1f_simple(file_nb, read_nb, nbbuffer);
                                                                                                                if (PEu == 1){
                                                                                                                types0 = dna_read[0] +  "\t0\t" + std::to_string(brstart) + "\t" + std::to_string(brend) + "\n";
                                                                                                                brstart = 0;
                                                                                                                brend = 0;
                                                                                                                print1f_types(types, types0, typesbuffer);

                                                                                                                types0 = "";
                                                                                                                }
                                                                                                                for(i=0;i<4;i++){
                                                                                                                                read_nb[i] = "";
                                                                                                                }
                                                                                        }
                                                                        }
                                                                /*for(i=0;i<4;i++){
                                                                        one_read[i] = "";
                                                                }*/
                                                                one_read = { "", "", "", "" };
                                                                line_cnt_mod = 0;
            }
                                                one_read[line_cnt_mod] = line;
                                                line_cnt_mod += 1;

                                }
        match_vector_forward.clear();
        match_vector_reverse.clear();
        line_cnt_mod = 0;
                                infile.close();
                                // cout << cnt << endl;
        }
        else cout << "Unable to open file";
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
        close(dnafile); close(rnafile); close(file_nb); close(types);

        //fclose(types);
        //types.close();

                                return 0;
}
