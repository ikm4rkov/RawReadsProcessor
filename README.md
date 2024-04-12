Version 1.0 (CPP) - only single-end file input, old verison (old read forming method, will skip the last read unless you add '\n@'.
Version 1.1 (CPP) - single-end and paired-end input. The algorithm is non-greedy, that's why EndsProcessor was added. ConvertDesqSeq splits the description sequence into two for EndsProcessor.
Allows grammar:
  * DNA-part start
  . RNA-part start
  < cut part start
  ? DNA check start
  ! RNA check part
  b bridge check start (difference from cut is the coordinates are saved, implies there is only one bridge)
  ( start for max mismatch limit
  ) end for max mismatch limit
  0-9 digits inside brackets, are max mismatch limit
Options:
  -s single-end mode (one input file; incompatible with paired-end mode -p)
  -p paired end modd (two input files; incompatible with single-end mode -s)
  -e separate; creates _RNA, _DNA and _gargbage files, first two are separated by bridge and processed with all grammar commands, last one is where either condition is not met.
  -t tabular; creates  types .TSV file in form of <ID>\t<BT (bridge true) or BF (bridge false)\t<bridge orientation, e.g. 10 means forward, and no reverse complemnet)\t<start coordinate of the bridge>\t<end coordinate of the bridge>
  -i input file for single-end mode (.FASTQ)
  -j first input file for paired-end mode (.FASTQ)
  -k second input file for paired-end mode (.FASTQ)
  -d description sequence (most symbols are special characters in Bash, so it is better to set this in apostrophes '; follows grammer) e.g.
  -l min length for DNA part (default 0)
  -k min length for RNA part (default 0)
If both -s and -p are abscent, program couldn't write (no rights) or reads in R1 and R2 are not sorted (ID from R1 doesn't correspond to ID from R2 at the same string) the program will exit and write the reason.

Endsprocessor (CPP) - position-based editor for restriction sites and other sequences.
Allows grammar:
  * DNA-part start
  . RNA-part start
  + add given characters
  - remove given number of characters
  s substitute given set of symbols
  | separates subsrquence which will be removed and the one which will be added, for example AG|AGCT means AG will be subsituted into AGCT.
  [ contain operator value, like +[TC] will add TC, -[2] remove 2 characters, s[AG|AGTC] substitute AG into AGCT
  ] closing operator value
Unlike the first program has a contants arguments order: DNA.FASTQ RNA.FASTQ DNA_desq_seq RNA_desc_seq.
This program creates _RS files which are results of its grammar operations and in addition will create oligos .TSV file. Oligos file is in format <Oligonucleotide (2-3 length)>\t<Description e.g. DNA_3'_Last2>\t<Count>. It can be processed and visualised to check what do the reads have at their ends and how does it correlate with an experimental design.

ConvertDesqSeq - is a simple sed command to remove unrelated to EndsProcessor operators. Implies your description sequence is in file Input, but you can change in to stdin.

