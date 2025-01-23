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



// Global vectors to store read data (single-end and paired-end)
// Each vector contains 4 strings (identifier, sequence, plus, quality)
// Used to store different types of reads
std::vector<std::string> one_read = { "", "", "", "" };  // Single-end read
std::vector<std::string> first_read = { "", "", "", "" };  // First of paired-end reads
std::vector<std::string> second_read = { "", "", "", "" };  // Second of paired-end reads
std::vector<std::string> dna_read = { "", "", "", "" };  // DNA part of the read
std::vector<std::string> dna_read_clear = { "", "", "", "" };  // Cleaned DNA part
std::vector<std::string> rna_read = { "", "", "", "" };  // RNA part of the read
std::vector<std::string> rna_read_clear = { "", "", "", "" };  // Cleaned RNA part
std::vector<std::string> read_nb = { "", "", "", "" };  // Read number
std::vector<std::string> empty_vec = {};  // Empty vector
std::vector<std::string> split_id = {};  // Split identifiers
std::vector<std::string> name_split = {};  // Split names
std::vector<std::string> cond_codes = {};  // Conditional sequence codes
std::set<string> condi_codes;  // Set of conditional sequence codes
map<string, int> count_SE_codes;  // Counter for single-end reads
map<string, int> count_PE_codes;  // Counter for paired-end reads

// Global variables for linker coordinates
int brstart = 0;  // Start of the linker
int brend = 0;  // End of the linker
int last_start;  // Last start position
int last_end;  // Last end position

// Approximate search function using the Shift-Or algorithm
// Returns a vector of positions where matches are found with errors (less than or equal to mismatch_limit)
void bitap_approximate_search(string pattern, string text, int mismatch_limit, std::vector<int> &match_vector) {
    int mask_length = pattern.length();  // Length of the pattern
    unsigned long *buffer;  // Buffer to store state
    unsigned long pattern_mask[255];  // Masks for characters
    int i, mismatch_iter;

    if (pattern.empty()) {
        cout << "Empty pattern provided" << endl;
        return;
    }
    if (mask_length > 63) {
        cout << "Pattern length exceeds the limit (63)" << endl;
        return;
    }

    // Allocate memory for the buffer and initialize it with ones (inverted bits)
    buffer = (unsigned long *) malloc((mismatch_limit + 1) * sizeof(*buffer));
    for (i = 0; i <= mismatch_limit; ++i) {
        buffer[i] = ~1;
    }

    // Initialize masks for characters
    for (i = 0; i <= 255; ++i) {
        pattern_mask[i] = ~0;
    }

    // Build masks for each character of the pattern (extended nucleotide codes)
    for (i = 0; i < mask_length; ++i) {
        if (pattern[i] == 'N') {
            pattern_mask['A'] &= ~(1UL << i);
            pattern_mask['G'] &= ~(1UL << i);
            pattern_mask['C'] &= ~(1UL << i);
            pattern_mask['T'] &= ~(1UL << i);
        }
        if (pattern[i] == 'W') {
            pattern_mask['A'] &= ~(1UL << i);
            pattern_mask['T'] &= ~(1UL << i);
        }
        if (pattern[i] == 'S') {
            pattern_mask['C'] &= ~(1UL << i);
            pattern_mask['G'] &= ~(1UL << i);
        }
        if (pattern[i] == 'M') {
            pattern_mask['A'] &= ~(1UL << i);
            pattern_mask['C'] &= ~(1UL << i);
        }
        if (pattern[i] == 'K') {
            pattern_mask['G'] &= ~(1UL << i);
            pattern_mask['T'] &= ~(1UL << i);
        }
        if (pattern[i] == 'R') {
            pattern_mask['A'] &= ~(1UL << i);
            pattern_mask['G'] &= ~(1UL << i);
        }
        if (pattern[i] == 'Y') {
            pattern_mask['C'] &= ~(1UL << i);
            pattern_mask['T'] &= ~(1UL << i);
        }
        if (pattern[i] == 'B') {
            pattern_mask['G'] &= ~(1UL << i);
            pattern_mask['C'] &= ~(1UL << i);
            pattern_mask['T'] &= ~(1UL << i);
        }
        if (pattern[i] == 'D') {
            pattern_mask['A'] &= ~(1UL << i);
            pattern_mask['G'] &= ~(1UL << i);
            pattern_mask['T'] &= ~(1UL << i);
        }
        if (pattern[i] == 'H') {
            pattern_mask['A'] &= ~(1UL << i);
            pattern_mask['C'] &= ~(1UL << i);
            pattern_mask['T'] &= ~(1UL << i);
        }
        if (pattern[i] == 'V') {
            pattern_mask['A'] &= ~(1UL << i);
            pattern_mask['G'] &= ~(1UL << i);
            pattern_mask['C'] &= ~(1UL << i);
        }
        pattern_mask[pattern[i]] &= ~(1UL << i);
    }

    // Main loop for the search
    unsigned long old_buffer;
    for (i = 0; text[i] != '\0'; ++i) {
        old_buffer = buffer[0];
        buffer[0] |= pattern_mask[text[i]];
        buffer[0] <<= 1;

        // Algorithm step with error counting
        for (mismatch_iter = 1; mismatch_iter <= mismatch_limit; ++mismatch_iter) {
            unsigned long tmp = buffer[mismatch_iter];
            buffer[mismatch_iter] = (old_buffer & (buffer[mismatch_iter] | pattern_mask[text[i]])) << 1;
            old_buffer = tmp;
        }

        // If a match is found, add the position to the vector
        if (0 == (buffer[mismatch_limit] & (1UL << mask_length))) {
            match_vector.push_back(i - mask_length + 1);
        }
    }

    // Free memory
    free(buffer);
}

// Global buffers for writing
std::string dnabuffer = "";  // Buffer for DNA part
std::string rnabuffer = "";  // Buffer for RNA part
std::string nbbuffer = "";  // Buffer for read numbers
std::string typesbuffer = "";  // Buffer for read types

// Function that writes reads to a file with buffering
void print1f_simple(int filei, std::vector<string> read, std::string &wbuffer) {
    // Add read lines to the buffer if read is not empty
    if (!read.empty()) {
        wbuffer += read[0] + "\n" + read[1] + "\n" + read[2] + "\n" + read[3] + "\n";
    }

    // Check if the buffer is large enough to write to file
    if (wbuffer.length() > 10000000 || (wbuffer.length() != 0 && read.empty())) {
        // Write the buffer contents to the file
        write(filei, wbuffer.c_str(), wbuffer.length());
        // Clear the buffer after writing
        wbuffer.clear();
    }
}

// Function that writes statistics (string) to a file with buffering
void print1f_types(int filei, string line, std::string &wbuffer) {
    // Add line to the buffer if not empty
    if (!line.empty()) {
        wbuffer += line;
    }

    // Check if the buffer is large enough or if the line is empty (flush buffer)
    if (wbuffer.length() > 10000000 || (wbuffer.length() != 0 && line.empty())) {
        // Write the buffer contents to the file
        write(filei, wbuffer.c_str(), wbuffer.length());
        // Clear the buffer after writing
        wbuffer.clear();
    }
}

// Function to add reverse complement to a string
void rcomplement(const string &a, string &out) {
    // Mapping of nucleotides to their complements (including extended codes)
    std::map<char, char> conv = {
        { 'A', 'T'}, { 'G', 'C' }, { 'C', 'G' }, { 'T', 'A' },
        { 'U', 'A' }, { 'M', 'K' }, { 'R', 'Y' }, { 'W', 'W' },
        { 'S', 'S' }, { 'Y', 'R' }, { 'K', 'M' }, { 'V', 'B' },
        { 'H', 'D' }, { 'D', 'H' }, { 'B', 'V' }, { 'N', 'N' }
    };

    // Traverse the input string in reverse and convert each character
    for (int i = a.length() - 1; i >= 0; i--) {
        out.push_back(conv[a[i]]);
    }
}

// Parser for description sequence (handles pattern matching)
void parser(int &success, const string &describe_seq, string &restofline, string &restoflineq, std::vector<int> &match_vector) {
    int cnt = 0;  // Counter for operations (e.g., successful match attempts)
    bool DNAu = false, RNAu = false, MMNu = false, bridgeu = false;
    string brseq = "", MMN = "";

    // Iterate over the description sequence to identify states and extract components
    for (int i = 0; i <= describe_seq.length(); i++) {
        string allseq = restofline;  // Copy of the rest of the sequence (for comparison)
        string allq = restoflineq;   // Copy of the quality string

        // Update states based on the description sequence characters
        if (describe_seq[i] == '*') {
            DNAu = true;
            RNAu = false;
            MMNu = false;
            bridgeu = false;
        } else if (describe_seq[i] == '.') {
            DNAu = false;
            RNAu = true;
            MMNu = false;
            bridgeu = false;
        } else if (describe_seq[i] == 'b') {
            bridgeu = true;
            MMNu = false;
        } else if (describe_seq[i] == '(') {
            MMNu = true;
        }

        // Collect bridge sequence if in bridge mode (not yet searching)
        if (bridgeu && describe_seq[i] != 'b' && !MMNu) {
            brseq += describe_seq[i];
        }

        // Collect MMN value for Hamming distance if in MMN mode
        if (MMNu && describe_seq[i] != '(') {
            MMN += describe_seq[i];
        }

        // If the closing parenthesis is reached, perform search
        if (describe_seq[i] == ')') {
            // Perform the bitap approximate search with the bridge sequence
            bitap_approximate_search(brseq, restofline, stoi(MMN), match_vector);

            if (match_vector.empty()) {
                success = 0;  // No match found
                brstart = -1;
                brend = -1;
                break;
            }

            // Update match positions and handle the sequence based on type (DNA/RNA)
            brstart = match_vector[0];
            brend = match_vector[0] + brseq.length();

            if (DNAu) {
                dna_read[1] += restofline.substr(0, match_vector[0]);
                dna_read[3] += restoflineq.substr(0, match_vector[0]);
            }

            if (RNAu) {
                rna_read[1] += restofline.substr(0, match_vector[0]);
                rna_read[3] += restoflineq.substr(0, match_vector[0]);
            }

            // Update the rest of the line after the match
            restofline = restofline.substr(match_vector[0] + brseq.length());
            restoflineq = restoflineq.substr(match_vector[0] + brseq.length());
            match_vector.clear();
            bridgeu = false;  // Reset bridge state
        }

        // If no more operators, append the remaining part of the sequence
        if (describe_seq[i] == '\0' && success == 1) {
            if (DNAu) {
                dna_read[1] += restofline;
                dna_read[3] += restoflineq;
                DNAu = false;
            } else if (RNAu) {
                rna_read[1] += restofline;
                rna_read[3] += restoflineq;
                RNAu = false;
            }
        }
    }
}

// Function for splitting a line by a character
size_t split(const std::string &txt, std::vector<std::string> &strs, char ch)
{
    size_t pos = txt.find(ch);
    size_t initialPos = 0;
    strs.clear();  // Очищаем вектор перед заполнением

    // Разбиваем строку
    while (pos != std::string::npos) {
        strs.push_back(txt.substr(initialPos, pos - initialPos));  // Добавляем подстроку в вектор
        initialPos = pos + 1;  // Обновляем начальную позицию для следующего сегмента

        pos = txt.find(ch, initialPos);  // Находим следующий символ разделителя
    }

    // Добавляем последнюю часть строки
    strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos));

    return strs.size();  // Возвращаем количество элементов в векторе
}

// General function to process and print read data
void processRead(std::vector<std::string>& one_read,
                 std::vector<std::string>& split_id,
                 const std::string& read_type,
                 const std::string& types0,
                 bool sepu, bool TSVu,
                 int& brstart, int& brend,
                 std::map<std::string, int>& count_SE_codes,
                 std::vector<std::string>& read_nb,
                 const int types, std::string& typesbuffer,
                 const int file_nb) {  // Добавили file_nb как параметр
    // Split the read ID and prepare the read_nb array
    split(one_read[0], split_id, ' ');
    read_nb[0] = split_id[0];
    read_nb[1] = one_read[1];
    read_nb[3] = one_read[3];
    read_nb[2] = "+";
    split_id.clear();

    // If output to file is enabled, print the read
    if (sepu == 1) {
        print1f_simple(file_nb, read_nb, nbbuffer);
    }

    // If TSV output is enabled, generate the tabular output
    if (TSVu == 1) {
        count_SE_codes[read_type] += 1;  // Update the read type counter
        brstart = 0;  // Reset the start position
        brend = 0;    // Reset the end position
        print1f_types(types, types0, typesbuffer);  // Print the formatted types output
    }

    // Clear the read_nb array for the next iteration
    for (int i = 0; i < 4; i++) {
        read_nb[i] = "";
    }
}

std::string get_filename_without_extension(const std::string &filepath) {
    size_t pos = filepath.find_last_of("/\\");
    std::string filename = (pos == std::string::npos) ? filepath : filepath.substr(pos + 1);
    size_t ext_pos = filename.find_last_of(".");
    return (ext_pos == std::string::npos) ? filename : filename.substr(0, ext_pos);
}


// Main function
int main(int argc, char** argv)
{
    // Options parser
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
    string output_dir = ".";  // Variable to store the output directory path
    int index;
    int c;
    string types0 = "";
    opterr = 0;

    // Parsing options
    while ((c = getopt(argc, argv, "hVsepvti:k:j:d:l:m:o:u:")) != -1) {
        switch (c) {
        case 'h':
            cout << "Tool for processing reads in .FASTQ format" << endl;
            cout << "Parameters: " << endl;
            cout << "-h: help and description option" << endl;
            cout << "-i input file path, must be in .FASTQ format" << endl;
            cout << "-s single-end reads mode, the output will be in 3 .FASTQ files - DNA- parts, RNA- parts, and NB (no bridge) reads. Can be used with -p to combine outputs" << endl;
            cout << "-p paired-end mode, output in .tsv file: read_id\t0|1|2 (1 for forward bridge, 2 for reverse, 0 for no bridge)" << endl;
            cout << "-t .TSV statistic format" << endl;
            cout << "-u conditional codes for bridge strand, comma-separated without whitespaces; default is 10,01,20,02" << endl;
            cout << "-d description sequence with special symbols for DNA/RNA parts" << endl;
            cout << "-l optional min length filter for total DNA- and RNA- parts, default is 0." << endl;
            cout << "-o output directory for saving the result files" << endl;
            cout << "-V version information" << endl;
            cout << "Usage: ./alpha1 -s -i <path to .fastq> -d <description sequence> -e -t -l <int; default 0> -m <int; default 0> -u <string; default 1,2>" << endl;
            break;
        case 'V':
            cout << "Version 1.3.3 beautified code, added -o option for output file, -V version for version control; for more details on options and examples use -h" << endl;
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
        case 'o':
            output_dir = optarg;  // Store the output directory path
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
        case '?':
            if (optopt == 'c')
                fprintf(stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint(optopt))
                fprintf(stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            return 1;
        default:
            abort();
        }
    }

    // Add logic for handling the output directory if needed
    if (!output_dir.empty()) {
        cout << "Output will be saved to directory: " << output_dir << endl;
    } else {
        cout << "No output directory specified. Results will be saved in the current directory." << endl;
    }
    // Output non-option arguments, if any
    for (int index = optind; index < argc; index++) {
        printf("Non-option argument %s\n", argv[index]);
    }

    // Print the description sequence
    cout << "Description sequence: " << describe_seq << endl;

    // Auxiliary variables
    string line; // Line for reading
    int i; // Iterator

    // Vectors to store match positions for different sequences
    vector<int> match_vector_forward;
    vector<int> match_vector_reverse;
    vector<int> match_vector_forward1;
    vector<int> match_vector_reverse1;
    vector<int> match_vector_forward2;
    vector<int> match_vector_reverse2;

    // Counter for lines in one read
    short line_cnt_mod = 0;

    // Progress tracking variable
    unsigned long long progress = 0;

    // Split the bridge codes (linker orientation codes)
    vector<string> cond_codes;
    split(bridge_codes, cond_codes, ',');

    // Variables for output file paths
    string s1 = "", s2 = "", s3 = "", s4 = "", s5 = "";
		// Handling output file paths based on SE (Single End) or PE (Paired End) mode
    if (SEu == 1) {  // Single End mode
        std::string filename = get_filename_without_extension(inputfilepath);  // Получаем имя файла без расширения
        s1 = output_dir + "/" + filename + ".DNA.fastq";
        s2 = output_dir + "/" + filename + ".RNA.fastq";
        s3 = output_dir + "/" + filename + ".garbage.fastq";
        s4 = output_dir + "/" + filename + ".types.tsv";
        s5 = output_dir + "/" + filename + ".codes.tsv";
    }
    if (PEu == 1) {  // Paired End mode
        std::string filename1 = get_filename_without_extension(inputfilepath1);  // Получаем имя для первого файла
        std::string filename2 = get_filename_without_extension(inputfilepath2);  // Получаем имя для второго файла
        s1 = output_dir + "/" + filename1 + ".DNA.fastq";
        s2 = output_dir + "/" + filename2 + ".RNA.fastq";
        s3 = output_dir + "/" + filename1 + ".garbage.fastqlike";  // Например, для второго файла
        s4 = output_dir + "/" + filename1 + ".types.tsv";
        s5 = output_dir + "/" + filename1 + ".codes.tsv";
		}


    // Process linker orientation nomenclature (F/R -> 1/2 or 10/20)
    unordered_set<string> condi_codes;
    for (int i = 0; i < cond_codes.size(); i++) {
        if (SEu == 1) {
            if (cond_codes[i] == "F") condi_codes.insert("1");  // Forward orientation
            if (cond_codes[i] == "R") condi_codes.insert("2");  // Reverse orientation
        }
        else if (PEu == 1) {
            if (cond_codes[i] == "F0") condi_codes.insert("10");  // Forward in first read
            if (cond_codes[i] == "R0") condi_codes.insert("20");  // Reverse in first read
            if (cond_codes[i] == "0F") condi_codes.insert("01");  // Forward in second read
            if (cond_codes[i] == "0R") condi_codes.insert("02");  // Reverse in second read
        }
    }

    // 10 nomenclature handling and default cases
    for (int i = 0; i < cond_codes.size(); i++) {
        if (SEu == 1) {
            if (cond_codes[i] == "10") condi_codes.insert("1");  // Forward in single end
            if (cond_codes[i] == "20") condi_codes.insert("2");  // Reverse in single end
            if (cond_codes[i] == "1") condi_codes.insert("1");   // Default forward
            if (cond_codes[i] == "2") condi_codes.insert("2");   // Default reverse
        }
        else if (PEu == 1 || cond_codes[i] == "10" || cond_codes[i] == "20" || cond_codes[i] == "01" || cond_codes[i] == "02") {
            condi_codes.insert(cond_codes[i]);  // Insert other valid codes
        }
    }

    // Check read permissions for input files (to ensure they are accessible)
    if (SEu == 1 && access(inputfilepath.c_str(), R_OK) != 0) {
        cerr << "Cannot read input file: " << inputfilepath << endl;
        exit(1);
    }

    if (PEu == 1) {
        if (access(inputfilepath1.c_str(), R_OK) != 0) {
            cerr << "Cannot read input file 1: " << inputfilepath1 << endl;
            exit(1);
        }
        if (access(inputfilepath2.c_str(), R_OK) != 0) {
            cerr << "Cannot read input file 2: " << inputfilepath2 << endl;
            exit(1);
        }
    }

    // Create and open input files
    fstream infile;
    fstream infile1;
    fstream infile2;
    infile.open(inputfilepath, std::ifstream::in);
    infile1.open(inputfilepath1, std::ifstream::in);
    infile2.open(inputfilepath2, std::ifstream::in);

    // Open output files for writing (using low-level file opening)
    const int dnafile = open(s1.c_str(), O_CREAT | O_WRONLY, 0644);
    const int rnafile = open(s2.c_str(), O_CREAT | O_WRONLY, 0644);
    const int file_nb = open(s3.c_str(), O_CREAT | O_WRONLY, 0644);
    const int types = open(s4.c_str(), O_CREAT | O_WRONLY, 0644);
    const int codes = open(s5.c_str(), O_CREAT | O_WRONLY, 0644);

    // Check write permissions for output files
		// cout << s1 << " " << s2 << " " << s3 << " " << s4 << " " << s5 << endl;
    if (access(s1.c_str(), W_OK) != 0 || access(s2.c_str(), W_OK) != 0 || access(s3.c_str(), W_OK) != 0 || access(s4.c_str(), W_OK) != 0 || access(s5.c_str(), W_OK) != 0) {
        cerr << "Cannot write to output file" << endl;
        exit(1);
    }

    // Reading the input file in Single End mode (SEu == 1)
    if (SEu == 1) {
        // Check if the input file is open
        if (infile.is_open()) {
            // Loop through each line in the input file
            while (getline(infile, line)) {
                // Store the current line in the one_read array (line_cnt_mod keeps track of the line number)
                one_read[line_cnt_mod] = line;

                // When we reach the 4th line of a read (line_cnt_mod == 3)
                if (line_cnt_mod == 3) {
                    // Display progress every million reads processed
                    if (progress % 1000000 == 0 && progress != 0)
                        cout << "Processed " << progress / 1000000 << " million reads" << endl;
                    progress += 1; // Increment the progress counter

                    // Retrieve the read sequence (second line of the read) and its quality scores (fourth line)
                    string read_seq = one_read[1];
                    string restofline = read_seq;
                    string restoflineq = one_read[3];
                    string restofline_rc = "";
                    string restoflineq_rc = "";
										// Reverse the sequence and the quality strings for reverse orientation
										rcomplement(restofline, restofline_rc);
										for (i = restofline_rc.length() - 1; i >= 0; i--) {
										    restoflineq_rc += restoflineq[i];
										}

                    // Flags to track the validity of different conditions
                    bool DNAu = false;
                    bool RNAu = false;
                    bool Cutu = false;
                    bool CheckDNAu = false;
                    bool CheckRNAu = false;
                    bool MMNu = false;
                    bool RNA_Length_Success = true;
                    bool DNA_Length_Success = true;

                    // Counters for success (forward and reverse orientation)
                    int successf = 1;
                    int successr = 1;

                    // Temporary strings for handling sequences and checks
                    string Cut = "";
                    string CheckDNA = "";
                    string CheckRNA = "";
                    string MMN = "";

                    // Clear the read arrays (for DNA, RNA, and general reads)
                    for (i = 0; i < 4; i++) {
                        dna_read[i] = "";
                        rna_read[i] = "";
                        read_nb[i] = "";
                    }

                    // Split the read ID into parts (e.g., SRR%d.%d format)
                    split(one_read[0], split_id, ' ');
                    dna_read[0] = split_id[0];
                    dna_read[2] = "+\0";  // Mark the DNA read as forward
                    rna_read[0] = split_id[0];
                    rna_read[2] = "+\0";  // Mark the RNA read as forward
                    split_id.clear();  // Clear the split ID vector

                    // Try parsing the read sequence in forward orientation (+ strand)
                    parser(successf, describe_seq, restofline, restoflineq, match_vector_forward);

                    // Check the length of the RNA and DNA sequences for validity
                    RNA_Length_Success = rna_read[1].length() >= min_length_rna;
                    DNA_Length_Success = dna_read[1].length() >= min_length_dna;

                    // If the forward orientation parsing failed, reset the read arrays
                    if (!successf) {
                        dna_read[1] = "";
                        dna_read[3] = "";
                        rna_read[1] = "";
                        rna_read[3] = "";

                        // Split the read ID again and reset the read arrays
                        split(one_read[0], split_id, ' ');
                        dna_read[0] = split_id[0];
                        dna_read[2] = "+\0";
                        rna_read[0] = split_id[0];
                        rna_read[2] = "+\0";
                    }
                    else {
                        // Update the last successful bridge start and end positions (if any)
                        last_start = brstart;
                        last_end = brend;
                    }

										//if (dna_read[0] == "@SRR10010328.15282064") cout << successr << " " << describe_seq << " " << restofline_rc << " " << restoflineq_rc << " " << match_vector_reverse.size() <<  endl;

                    // Try parsing the reverse orientation of the read (- strand)
                    parser(successr, describe_seq, restofline_rc, restoflineq_rc, match_vector_reverse);

                    // Check the length of RNA and DNA sequences for reverse orientation
                    if (successr) {
                        RNA_Length_Success = rna_read[1].length() >= min_length_rna;
                        DNA_Length_Success = dna_read[1].length() >= min_length_dna;
                    }
										//if (dna_read[0] == "@SRR10010328.15282064") cout << successf << " " << successr << endl;

// Check if the read matches exactly once, either in forward (successf == 1) or reverse (successr == 1) orientation
                    if (successf + successr == 1) {
                        // Handle the case where the read matches in the forward orientation
                        if (successf == 1) {
                            // Only process the read if both DNA and RNA lengths are valid
                            if (DNA_Length_Success && RNA_Length_Success) {
                                // Check if the condition for splitting (via "1" code) is met
                                int is_in = condi_codes.count("1");
                                // Debug output for checking the read
                                // cout << "Processing forward: " << one_read[0] << " Length: " << dna_read[1].length() << ", " << rna_read[1].length() << endl;

                                // If SEPU is enabled, print the sequences to the respective files
                                if (sepu == 1 && is_in > 0) {
                                    print1f_simple(dnafile, dna_read, dnabuffer);
                                    print1f_simple(rnafile, rna_read, rnabuffer);
                                }

                                // If TSVu is enabled, create a tabular output for the read type and increment the appropriate count
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tF-DR\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" +
                                             std::to_string(dna_read[1].length()) + "\t" + std::to_string(rna_read[1].length()) + "\n";
                                    count_SE_codes["F-DR"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                            }

                            // If only the DNA part is valid (RNA part failed)
                            if (DNA_Length_Success && !RNA_Length_Success) {
                                // Split the ID and prepare the read for output
                                split(one_read[0], split_id, ' ');
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                // Debug output for checking the read
                                // cout << "Processing F-D0: " << one_read[0] << endl;

                                // Print to output file if SEPU is enabled
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);

                                // If TSVu is enabled, generate a table output
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tF-D0\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" +
                                             std::to_string(dna_read[1].length()) + "\t" + std::to_string(rna_read[1].length()) + "\n";
                                    count_SE_codes["F-D0"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                // Clear the temporary read array
                                for (i = 0; i < 4; i++) {
                                    read_nb[i] = "";
                                }
                            }

                            // If only the RNA part is valid (DNA part failed)
                            if (!DNA_Length_Success && RNA_Length_Success) {
                                // Split the ID and prepare the read for output
                                split(one_read[0], split_id, ' ');
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                // Debug output for checking the read
                                // cout << "Processing F-0R: " << one_read[0] << endl;

                                // Print to output file if SEPU is enabled
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);

                                // If TSVu is enabled, generate a table output
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tF-0R\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" +
                                             std::to_string(dna_read[1].length()) + "\t" + std::to_string(rna_read[1].length()) + "\n";
                                    count_SE_codes["F-0R"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                // Clear the temporary read array
                                for (i = 0; i < 4; i++) {
                                    read_nb[i] = "";
                                }
                            }

                            // If neither DNA nor RNA part is valid
                            if (!DNA_Length_Success && !RNA_Length_Success) {
                                // Split the ID and prepare the read for output
                                split(one_read[0], split_id, ' ');
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                // Debug output for checking the read
                                // cout << "Processing F-00: " << one_read[0] << endl;

                                // Print to output file if SEPU is enabled
                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);

                                // If TSVu is enabled, generate a table output
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tF-00\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" +
                                             std::to_string(dna_read[1].length()) + "\t" + std::to_string(rna_read[1].length()) + "\n";
                                    count_SE_codes["F-00"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                                // Clear the temporary read array
                                for (i = 0; i < 4; i++) {
                                    read_nb[i] = "";
                                }
                            }
                        }

                        // Handle the case where the read matches in the reverse orientation
                        else {
                            // Only process the read if both DNA and RNA lengths are valid
                            if (DNA_Length_Success && RNA_Length_Success) {
                                int is_in = condi_codes.count("2");
                                // Debug output for checking the reverse read
                                // cout << "Processing reverse: " << one_read[0] << " Length: " << dna_read[1].length() << ", " << rna_read[1].length() << endl;

                                // If SEPU is enabled, print the sequences to the respective files
                                if (sepu == 1 && is_in > 0) {
                                    print1f_simple(dnafile, dna_read, dnabuffer);
                                    print1f_simple(rnafile, rna_read, rnabuffer);
                                }

                                // If TSVu is enabled, create a tabular output for the reverse read
                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tR-DR\t" + std::to_string(one_read[1].size() - brend + 1)  + "\t" +
                                             std::to_string(one_read[1].size() - brstart) + "\t" + std::to_string(dna_read[1].length()) +
                                             "\t" + std::to_string(rna_read[1].length()) + "\n";
                                    count_SE_codes["R-DR"] += 1;
                                    brstart = 0;
                                    brend = 0;
                                    print1f_types(types, types0, typesbuffer);
                                    types0 = "";
                                }
                            }

                            // If only the DNA part is valid (RNA part failed)
                            if (DNA_Length_Success && !RNA_Length_Success) {
                                split(one_read[0], split_id, ' ');
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                // Debug output for checking the reverse F-D0
                                // cout << "Processing R-D0: " << one_read[0] << endl;

                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);

                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tR-D0\t" + std::to_string(one_read[1].size() - brend + 1)  + "\t" +
                                             std::to_string(one_read[1].size() - brstart) + "\t" + std::to_string(dna_read[1].length()) +
                                             "\t" + std::to_string(rna_read[1].length()) +  "\n";
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

                            // If only the RNA part is valid (DNA part failed)
                            if (!DNA_Length_Success && RNA_Length_Success) {
                                split(one_read[0], split_id, ' ');
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                // Debug output for checking the reverse F-0R
                                // cout << "Processing R-0R: " << one_read[0] << endl;

                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);

                                if (TSVu == 1) {
                                    types0 = dna_read[0] +  "\tR-0R\t" + std::to_string(one_read[1].size() - brend + 1)  + "\t" +
                                             std::to_string(one_read[1].size() - brstart) + "\t" + std::to_string(dna_read[1].length()) +
                                             "\t" + std::to_string(rna_read[1].length()) +  "\n";
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

                            // If neither DNA nor RNA part is valid
                            if (!DNA_Length_Success && !RNA_Length_Success) {
                                split(one_read[0], split_id, ' ');
                                read_nb[0] = split_id[0];
                                read_nb[1] = one_read[1];
                                read_nb[3] = one_read[3];
                                read_nb[2] = "+";
                                split_id.clear();
                                // Debug output for checking the reverse F-00
                                // cout << "Processing R-00: " << one_read[0] << endl;

                                if (sepu == 1) print1f_simple(file_nb, read_nb, nbbuffer);

                                if (TSVu == 1) {
                                    types0 = dna_read[0] + "\tR-00\t" + std::to_string(one_read[1].size() - brend + 1)  + "\t" +
                                             std::to_string(one_read[1].size() - brstart) + "\t" + std::to_string(dna_read[1].length()) +
                                             "\t" + std::to_string(rna_read[1].length()) + "\n";
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
										else {  
												// Handle cases where more than one linker is found
												if (successf + successr > 1) {
                            // General case for multiple linkers found (FR)
                            processRead(one_read, split_id, "FR", dna_read[0] + "\tFR\t-1\t-1\t-1\t-1", sepu, TSVu, brstart, brend, count_SE_codes, read_nb, types, typesbuffer, file_nb);
                        }

                        // Handle cases where no linker is found
                        if (successf + successr == 0) {
                            // General case for no linkers found (0)
                            processRead(one_read, split_id, "0", dna_read[0] + "\t0\t-1\t-1\t-1\t-1", sepu, TSVu, brstart, brend, count_SE_codes, read_nb, types, typesbuffer, file_nb);
                        }

                    }
										// Reset read data after processing
                    one_read = { "", "", "", "" };  // Clear the read data
										line_cnt_mod = -1;
                }
                line_cnt_mod += 1;  // Increment line counter
            }
            match_vector_forward.clear();  // Clear forward match vector
            match_vector_reverse.clear();  // Clear reverse match vector
            line_cnt_mod = 0;  // Reset the line counter
            infile.close();  // Close the input file
// cout << cnt << endl;  // Debug output (if needed)

        }

    }
    if (PEu == 1) {  // Check if paired-end mode is enabled
        string line1, line2;  // Strings to hold the lines from both files
        while (getline(infile1, line1) && getline(infile2, line2)) {  // Process each line from both files
            // Store the current lines from both files as paired reads
						first_read[line_cnt_mod] = line1;
            second_read[line_cnt_mod] = line2;

            if (line_cnt_mod == 3) {  // Process the 3rd line (after reading 3 lines)
                if (progress % 1000000 == 0 && progress != 0) {
                    cout << "Processed " << progress / 1000000 << " million reads" << endl;
                }

                // Split the read IDs from both lines and check if they match
                split(first_read[0], split_id, ' ');
                first_read[0] = split_id[0];  // Get the ID for the first read
                split_id.clear();
                split(second_read[0], split_id, ' ');
                second_read[0] = split_id[0];  // Get the ID for the second read
                split_id.clear();

                if (first_read[0] != second_read[0]) {  // Ensure the read IDs match between the two files
                    cerr << first_read[0] << " " << second_read[0] << endl;
                    cerr << "Input files are not sorted, different read IDs in the lines" << endl;
                    exit(1);  // Exit if IDs do not match
                }

                bool RNA_Length_Success = true;  // Flags for RNA sequence success
                bool DNA_Length_Success = true;  // Flags for DNA sequence success

                progress++;  // Increment progress counter

                // Extract the DNA sequences and their qualities for both reads
                string read_seq1 = first_read[1];
                string restofline1 = read_seq1;
                string restoflineq1 = first_read[3];
                string restofline_rc1, restoflineq_rc1;  // Reverse complement of the first read
                rcomplement(restofline1, restofline_rc1);  // Compute reverse complement for the first read

                string read_seq2 = second_read[1];
                string restofline2 = read_seq2;
                string restoflineq2 = second_read[3];
                string restofline_rc2, restoflineq_rc2;  // Reverse complement of the second read
                rcomplement(restofline2, restofline_rc2);  // Compute reverse complement for the second read

                // Reset read arrays for both DNA and RNA
                fill(begin(dna_read), end(dna_read), "");
                fill(begin(rna_read), end(rna_read), "");
                fill(begin(read_nb), end(read_nb), "");

                // Split the ID for both reads (SRR%d.%d format)
                split(first_read[0], split_id, ' ');
                dna_read[0] = split_id[0];
                dna_read[2] = "+";  // Set the strand information for DNA
                rna_read[0] = split_id[0];
                rna_read[2] = "+";  // Set the strand information for RNA

                // Flags for successful parsing of both forward and reverse strands
                int success1F = 1, success1R = 1, success2F = 1, success2R = 1;

                // Reverse complement quality strings for both reads
                restoflineq_rc1 = string(restoflineq1.rbegin(), restoflineq1.rend());
                restoflineq_rc2 = string(restoflineq2.rbegin(), restoflineq2.rend());

                // Parse the first read (forward strand)
                parser(success1F, describe_seq, restofline1, restoflineq1, match_vector_forward1);
                match_vector_forward1.clear();  // Clear the match vector

                // If the forward strand parsing failed, clear the corresponding arrays
                if (success1F != 1) {
                    dna_read[1] = "";
                    dna_read[3] = "";
                    rna_read[1] = "";
                    rna_read[3] = "";
                } else {
                    // Store successful DNA and RNA reads
                    copy(begin(dna_read), end(dna_read), begin(dna_read_clear));
                    copy(begin(rna_read), end(rna_read), begin(rna_read_clear));
                    last_start = brstart;
                    last_end = brend;
                }

                // Parse the first read (reverse strand)
                parser(success1R, describe_seq, restofline_rc1, restoflineq_rc1, match_vector_reverse1);
                match_vector_reverse1.clear();  // Clear the match vector

                // If reverse strand parsing failed, clear the corresponding arrays
                if (success1R != 1) {
                    dna_read[1] = "";
                    dna_read[3] = "";
                    rna_read[1] = "";
                    rna_read[3] = "";
                } else {
                    // Store successful DNA and RNA reads
                    copy(begin(dna_read), end(dna_read), begin(dna_read_clear));
                    copy(begin(rna_read), end(rna_read), begin(rna_read_clear));
                    last_start = brstart;
                    last_end = brend;
                }

                // Parse the second read (forward strand)
                parser(success2F, describe_seq, restofline2, restoflineq2, match_vector_forward2);
                match_vector_forward2.clear();  // Clear the match vector

                // If the forward strand parsing failed, clear the corresponding arrays
                if (success2F != 1) {
                    dna_read[1] = "";
                    dna_read[3] = "";
                    rna_read[1] = "";
                    rna_read[3] = "";
                } else {
                    // Store successful DNA and RNA reads
                    copy(begin(dna_read), end(dna_read), begin(dna_read_clear));
                    copy(begin(rna_read), end(rna_read), begin(rna_read_clear));
                    last_start = brstart;
                    last_end = brend;
                }

                // Parse the second read (reverse strand)
                parser(success2R, describe_seq, restofline_rc2, restoflineq_rc2, match_vector_reverse2);
                match_vector_reverse2.clear();  // Clear the match vector

                // If reverse strand parsing failed, clear the corresponding arrays
                if (success2R != 1) {
                    dna_read[1] = "";
                    dna_read[3] = "";
                    rna_read[1] = "";
                    rna_read[3] = "";
                } else {
                    // Store successful DNA and RNA reads
                    copy(begin(dna_read), end(dna_read), begin(dna_read_clear));
                    copy(begin(rna_read), end(rna_read), begin(rna_read_clear));
                }

                if (success1F + success1R + success2F + success2R == 1) {  // Check if exactly one of the four conditions is true
                    RNA_Length_Success = rna_read_clear[1].length() >= min_length_rna;  // Check if RNA length meets the minimum requirement
                    DNA_Length_Success = dna_read_clear[1].length() >= min_length_dna;  // Check if DNA length meets the minimum requirement

                    if (RNA_Length_Success && DNA_Length_Success) {  // Proceed only if both RNA and DNA lengths are sufficient
                        if (success1F == 1) {  // If the first forward read was successful
                            int is_in = condi_codes.count("10");  // Check if condition "10" exists in condi_codes
                            if (sepu == 1 && is_in) {  // If SEPU flag is set and condition "10" exists
                                print1f_simple(dnafile, dna_read_clear, dnabuffer);  // Print the DNA sequence to file
                                print1f_simple(rnafile, rna_read_clear, rnabuffer);  // Print the RNA sequence to file
                            }
                            if (TSVu == 1) {  // If TSVu flag is set
                                types0 = dna_read[0] +  "\tF-0-DR\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";  // Format types0 string
                                count_PE_codes["F-0-DR"] += 1;  // Increment count for "F-0-DR"
                                brstart = 0;  // Reset start position
                                brend = 0;    // Reset end position
                                print1f_types(types, types0, typesbuffer);  // Print types0 to types buffer
                                types0 = "";  // Clear types0
                            }
                        }

                        // The following blocks are similar for success1R, success2F, and success2R
                        if (success1R == 1) {  // If the first reverse read was successful
                            int is_in = condi_codes.count("20");
                            if (sepu == 1 && is_in) {
                                print1f_simple(dnafile, dna_read_clear, dnabuffer);
                                print1f_simple(rnafile, rna_read_clear, rnabuffer);
                            }
                            if (TSVu == 1) {
                                types0 = dna_read[0] +  "\tR-0-DR\t" + std::to_string(first_read[1].size() - last_end + 1) + "\t" + std::to_string(first_read[1].size() - last_start) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                count_PE_codes["R-0-DR"] += 1;
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }

                        if (success2F == 1) {  // If the second forward read was successful
                            int is_in = condi_codes.count("01");
                            if (sepu == 1 && is_in) {
                                print1f_simple(dnafile, dna_read_clear, dnabuffer);
                                print1f_simple(rnafile, rna_read_clear, rnabuffer);
                            }
                            if (TSVu == 1) {
                                types0 = dna_read[0] +  "\t0-F-DR\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                count_PE_codes["0-F-DR"] += 1;
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }

                        if (success2R == 1) {  // If the second reverse read was successful
                            int is_in = condi_codes.count("02");
                            if (sepu == 1 && is_in) {
                                print1f_simple(dnafile, dna_read_clear, dnabuffer);
                                print1f_simple(rnafile, rna_read_clear, rnabuffer);
                            }
                            if (TSVu == 1) {
                                types0 = dna_read[0] +  "\t0-R-DR\t" + std::to_string(first_read[1].size() - brend + 1) + "\t" + std::to_string(first_read[1].size() - brstart) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                count_PE_codes["0-R-DR"] += 1;
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }
                    } else {  // If either RNA or DNA length is insufficient
                        if (success1F == 1) {  // If the first forward read was successful
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
                                    types0 = dna_read[0] +  "\tF-0-D0"  + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                    count_PE_codes["F-0-D0"] += 1;
                                }
                                if (!DNA_Length_Success && RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tF-0-0R" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                    count_PE_codes["F-0-0R"] += 1;
                                }
                                if (!DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tF-0-00"  + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                    count_PE_codes["F-0-00"] += 1;
                                }
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }

                        // Repeating structure for success1R, success2F, and success2R
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
                                    types0 = dna_read[0] +  "\tR-0-D0\t" + std::to_string(first_read[1].size() - last_end + 1) + "\t" + std::to_string(first_read[1].size() - last_start) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                    count_PE_codes["R-0-D0"] += 1;
                                }
                                if (!DNA_Length_Success && RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tR-0-0R\t" + std::to_string(first_read[1].size() - last_end + 1) + "\t" + std::to_string(first_read[1].size() - last_start) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                    count_PE_codes["R-0-0R"] += 1;
                                }
                                if (!DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\tR-0-00\t" + std::to_string(first_read[1].size() - last_end + 1) + "\t" + std::to_string(first_read[1].size() - last_start) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
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
                                    types0 = dna_read[0] +  "\t0-F-D0\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length())  + "\n";
                                    count_PE_codes["0-F-D0"] += 1;
                                }
                                if (!DNA_Length_Success && RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-F-0R\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length())  + "\n";
                                    count_PE_codes["0-F-0R"] += 1;
                                }
                                if (!DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-F-00\t" + std::to_string(last_start + 1) + "\t" + std::to_string(last_end) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length())  + "\n";
                                    count_PE_codes["0-F-00"] += 1;
                                }
                                brstart = 0;
                                brend = 0;
                                print1f_types(types, types0, typesbuffer);
                                types0 = "";
                            }
                        }

                        if (success2R == 1) {  // Process second reverse read with length checks
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
                                    types0 = dna_read[0] +  "\t0-R-D0\t" + std::to_string(first_read[1].size() - brend + 1) + "\t" + std::to_string(first_read[1].size() - brstart) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length())  + "\n";
                                    count_PE_codes["0-R-D0"] += 1;
                                }
                                if (!DNA_Length_Success && RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-R-0R\t" + std::to_string(first_read[1].size() - brend + 1) + "\t" + std::to_string(first_read[1].size() - brstart) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
                                    count_PE_codes["0-R-0R"] += 1;
                                }
                                if (!DNA_Length_Success && !RNA_Length_Success) {
                                    types0 = dna_read[0] +  "\t0-R-00\t" + std::to_string(first_read[1].size() - brend + 1) + "\t" + std::to_string(first_read[1].size() - brstart) + "\t" + std::to_string(dna_read_clear[1].length()) + "\t" + std::to_string(rna_read_clear[1].length()) + "\n";
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
                    if (sepu == 1) {  // If SEPU flag is set
                        read_nb[0] = split_id[0];  // Set the first element of the read
                        read_nb[1] = read_seq1 + "|" + read_seq2;  // Concatenate read sequences
                        read_nb[3] = first_read[3] + "|" + second_read[3];  // Combine additional info for first and second read
                        read_nb[2] = "+";  // Set the orientation of the reads
                        print1f_simple(file_nb, read_nb, nbbuffer);  // Print the data to the output file
                        split_id.clear();  // Clear split_id vector
                    }

                    if (TSVu == 1) {  // If TSVu flag is set
                        // Case 1: No successful read pairs
                        if (success1F + success1R + success2F + success2R == 0) {
                            types0 = read_nb[0] + "\t0-0\t-1\t-1\t-1\t-1" + "\n";  // Assign type 0-0 (no successful reads)
                            count_PE_codes["0-0"] += 1;  // Increment the count for this type
                        }

                        // Case 2: First forward and second forward reads are successful, others are not
                        if (success1F && !success1R && success2F && !success2R) {
                            types0 = read_nb[0] + "\tF-F\t-1\t-1\t-1\t-1" + "\n";  // Assign type F-F
                            count_PE_codes["F-F"] += 1;  // Increment the count for this type
                        }

                        // Case 3: First reverse and second reverse reads are successful, others are not
                        if (!success1F && success1R && !success2F && success2R) {
                            types0 = read_nb[0] + "\tR-R\t-1\t-1\t-1\t-1" + "\n";  // Assign type R-R
                            count_PE_codes["R-R"] += 1;  // Increment the count for this type
                        }

                        // Case 4: First forward and second reverse reads are successful
                        if (success1F && !success1R && !success2F && success2R) {
                            types0 = read_nb[0] + "\tF-R\t-1\t-1\t-1\t-1" + "\n";  // Assign type F-R
                            count_PE_codes["F-R"] += 1;  // Increment the count for this type
                        }

                        // Case 5: First reverse and second forward reads are successful
                        if (!success1F && success1R && success2F && !success2R) {
                            types0 = read_nb[0] + "\tR-F\t-1\t-1\t-1\t-1" + "\n";  // Assign type R-F
                            count_PE_codes["R-F"] += 1;  // Increment the count for this type
                        }

                        // Case 6: Only second forward and second reverse reads are successful
                        if (!success1F && !success1R && success2F && success2R) {
                            types0 = read_nb[0] + "\t0-FR\t-1\t-1\t-1\t-1" + "\n";  // Assign type 0-FR (second pair successful)
                            count_PE_codes["0-FR"] += 1;  // Increment the count for this type
                        }

                        // Case 7: Both first and second forward reads are successful, others are not
                        if (success1F && success1R && !success2F && !success2R) {
                            types0 = read_nb[0] + "\tFR-0\t-1\t-1\t-1\t-1" + "\n";  // Assign type FR-0
                            count_PE_codes["FR-0"] += 1;  // Increment the count for this type
                        }

                        // Case 8: First reverse and both second forward and second reverse reads are successful
                        if (!success1F && success1R && success2F && success2R) {
                            types0 = read_nb[0] + "\tR-FR\t-1\t-1\t-1\t-1" + "\n";  // Assign type R-FR
                            count_PE_codes["R-FR"] += 1;  // Increment the count for this type
                        }

                        // Case 9: First forward and both second forward and second reverse reads are successful
                        if (success1F && success1R && !success2F && success2R) {
                            types0 = read_nb[0] + "\tF-FR\t-1\t-1\t-1\t-1" + "\n";  // Assign type F-FR
                            count_PE_codes["F-FR"] += 1;  // Increment the count for this type
                        }

                        // Case 10: First forward and both first reverse and second reverse reads are successful
                        if (success1F && success1R && success2F && !success2R) {
                            types0 = read_nb[0] + "\tFR-R\t-1\t-1\t-1\t-1" + "\n";  // Assign type FR-R
                            count_PE_codes["FR-R"] += 1;  // Increment the count for this type
                        }

                        // Case 11: First forward and both first reverse and second forward reads are successful
                        if (success1F && success1R && success2F && success2R) {
                            types0 = read_nb[0] + "\tFR-F\t-1\t-1\t-1\t-1" + "\n";  // Assign type FR-F
                            count_PE_codes["FR-F"] += 1;  // Increment the count for this type
                        }

                        // Case 12: All pairs of reads are successful
                        if (success1F && success1R && success2F && success2R) {
                            types0 = read_nb[0] + "\tFR-FR\t-1\t-1\t-1\t-1" + "\n";  // Assign type FR-FR
                            count_PE_codes["FR-FR"] += 1;  // Increment the count for this type
                        }

                        // Reset variables to prepare for the next iteration
                        brstart = 0;
                        brend = 0;
                        last_start = 0;
                        last_end = 0;
                        print1f_types(types, types0, typesbuffer);  // Print the types information
                        types0 = "";  // Clear the types0 string
                    }

                }
								// Reset read variables for the next iteration
                first_read = { "", "", "", "" };
                second_read = { "", "", "", "" };
                line_cnt_mod = -1;  // Reset line count modifier

            }
						line_cnt_mod += 1;  // Increment line count modifier for the next iteration
        }
// Clear the match vectors to release memory and reset their content
        match_vector_forward1.clear();  // Clear forward match vector for the first read
        match_vector_reverse1.clear();  // Clear reverse match vector for the first read
        match_vector_forward2.clear();  // Clear forward match vector for the second read
        match_vector_reverse2.clear();  // Clear reverse match vector for the second read

        line_cnt_mod = 0;  // Reset the line count modifier to start fresh for the next run

// Close the input files
        infile1.close();  // Close the first input file
        infile2.close();  // Close the second input file
    }

// Check if no processing mode was specified (both SEu and PEu are 0)
    if (SEu == 0 && PEu == 0)
        cout << "No mode specified";  // Print a message if no mode is set

// If there is any remaining data in the buffers, print it out to the corresponding files
    print1f_simple(dnafile, empty_vec, dnabuffer);  // Print remaining DNA data to the DNA file
    print1f_simple(rnafile, empty_vec, rnabuffer);  // Print remaining RNA data to the RNA file
    print1f_simple(file_nb, empty_vec, nbbuffer);  // Print remaining non-base data to the non-base file
    print1f_types(types, types0, typesbuffer);  // Print remaining type data to the types file

// Close the output files after all data has been written
    close(dnafile);  // Close the DNA file
    close(rnafile);  // Close the RNA file
    close(file_nb);  // Close the non-base file
    close(types);    // Close the types file

// Prepare and write the code counts to the 'codes' file
    string code_out;

// If SE mode is active, write SE code counts to the 'codes' file
    if (SEu == 1) {
        for (auto element : count_SE_codes) {  // Loop through all SE code counts
            code_out = element.first + "\t" + to_string(element.second) + "\n";  // Format the code and count
            write(codes, code_out.c_str(), code_out.length());  // Write the formatted string to the 'codes' file
        }
    }

// If PE mode is active, write PE code counts to the 'codes' file
    if (PEu == 1) {
        for (auto element : count_PE_codes) {  // Loop through all PE code counts
            code_out = element.first + "\t" + to_string(element.second) + "\n";  // Format the code and count
            write(codes, code_out.c_str(), code_out.length());  // Write the formatted string to the 'codes' file
        }
    }

// Close the 'codes' file after writing the code counts
    close(codes);


// Return 0 indicating successful execution
    return 0;
}
