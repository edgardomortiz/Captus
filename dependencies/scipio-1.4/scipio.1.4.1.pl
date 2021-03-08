#!/usr/bin/perl -w
################################################################################
#
#   Name:    SCIPIO - Eukaryotic Gene Identification
#   Project: Gene Prediction with Protein Family Patterns 
#   Author:  Oliver Keller, Klas Hatje
#   Date:    2013-05-24
#   Version: 1.4.1
#
#
#   scipio.pl: [<options>] <target> <query> 
#
#   <target> is a collection of DNA sequences (a genome or part of it)
#   <query> is a collection of protein query sequences
#   both in FASTA format
#
#   This script runs BLAT and, for each queried protein sequence, returns
#   the correct, i.e. best-scoring location, if one was found.
#   If BLAT only returns partial matches on multiple targets, it tries to 
#   assemble them by matching to the ends of each target.
#   
#   Also, codons that are split by an intron (and therefore not recognized
#   by BLAT) are added to the reported location, giving an exact gene structure
#   of the query protein.
#
#   The script produces output in YAML format (see documentation)
#   
#   For usage information, run scipio.pl without parameters
#
#   Counting of bps/aas: interval counting (see documentation)
#


use strict;
use File::stat;
use Getopt::Long;
use List::Util('sum', 'max', 'min', 'reduce');
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use YAML;
if ($YAML::VERSION gt "0.60") {
    use YAML::Dumper::Base;
}
 

my $QPFX="prot";
my $TPFX="dna";
my $TRPFX="trans";
my $QNPFX="nucl";


################################## user options ################################

my $BLAT_BIN = "blat";      # name of blat executable, dont run if empty
my $BLAT_OUTPUT = "";       # name of blat  output file
my $BLAT_TILESIZE;
my $ADDED_MATCHING_MAXLEN;

my $BLAT_MIN_IDENTITY;      # minimal hit accuracy in %
my $BLAT_MIN_SCORE = 15;
my $BLAT_PARAMS = "";       # parameters passed to BLAT directly
my $FORCE_NEW=0;            # overwrite previous output file, always run blat; automatically
                            # chosen if no BLAT output is specified;

my $KEEP_BLATOUT;           # keep the BLAT output, automatically chosen 
                            # when BLAT output is specified 
my $TR_TABLE = 1;
my $VERBOSE = 0;
# my $JOIN_FS = 1;
my $SPLIT_ON_FS = 0;
# my $SPLIT_STRING = "_split:";  # split mode disabled
my $MIN_IDENTITY;
my $DEFAULT_MIN_IDENTITY = 90;
my $MIN_COVERAGE = 60;
my $MIN_BEST_SCORE = 0.3;
my $MAX_MISMATCH = 0;
# my $ADD_MARGINS = 1;
my $MAX_CHECK_OVERLAP = 300;
my $DEFAULT_REGION_SIZE = 1000;
my $REGION_SIZE = $DEFAULT_REGION_SIZE;     # this much is shown in log
# my $SHOWALL_SCORE = 1.0;
my $MULTIPLE_RESULTS = 0;
my $ALLOW_MISSING_TARGETS = 0;
my $PRED_OVERLAP=2;    # how many codons of prediction are checked for intron start/end

my $HIT_BLESS_LIST =
    "ID,status,reason,${QPFX}_len,${QPFX}_start,${QPFX}_end,${QPFX}_seq,".
    "target,target_len,strand,${TPFX}_start,${TPFX}_end,matches,mismatches,".
    "undetermined,unmatched,additional,score,upstream,upstream_gap,matchings,stopcodon,downstream,downstream_gap";
my %MATCH_BLESS_LISTS = ( 
    exon =>   
    "type,${QNPFX}_start,${QNPFX}_end,${TPFX}_start,${TPFX}_end,seq,${QPFX}_start,${QPFX}_end,".
    "seqshifts,mismatchlist,undeterminedlist,inframe_stopcodons,translation,overlap",
    intron => "type,${QNPFX}_start,${TPFX}_start,${TPFX}_end,seq",
    gap =>    "type,${QNPFX}_start,${QNPFX}_end,${TPFX}_start,${TPFX}_end,seq",
    any =>    "type,${QNPFX}_start,${QNPFX}_end,${TPFX}_start,${TPFX}_end,seq");

my $HIDE_UNDEF = 0;         # show also undefined keys and empty lists
my $HIDE_DEFAULTS = 0;      # show also keys with default values
my $DONT_DIE = 1;
my $SHOW_BLATLINE = 0;
my $BLAT_ONLY = 0;
my $ADDED_INTRON_MAXLEN = 500;
my $ONE_TARGET_ONLY = 0;
my %VALUE_FORMAT = ( "strand" => { '+' => '+', '-' => '-' });
my $VERBAL_STRANDS = 0;
my $SHOW_FRAME_IN_TRANS = 0;
my $MIN_INTRON_LEN = 22; # maximal length (in bases) of non GT..AG gap (-> seqshift)

my $MIN_SEQ_COVERAGE = 0;
my $HITBUFFER_SIZE = 75000;  # this much is loaded around a blat hit; 
                             # also the maximum length of an intron across different targets
my $ASSEMBLE_SIZE = $HITBUFFER_SIZE;

# cost for non-matches in NW-alignment
my $NW_MISMATCH = 1.0; 
my $NW_PROT_INS = 1.5;  # one extra aminoacid
my $NW_PROT_GAP = 1.1;  # one aminoacid missing
my $NW_FS       = 2.5;  # one or two extra nucleotides / one or two nucleotides missing
my $NW_STP      = 2.5;

my $NW_INTRON    = 2.0;
my $NW_ANY_SPL   = 1.6;
# my $NW_DSS = { "gt" => 0.5, "gc" => 0.6, "ga" => 0.7, "gg" => 0.7,
# 	       "gn" => 0.6, "nt" => 0.8, "nn" => 0.9 };
# my $NW_ASS = { "ag" => 0.5, "ng" => 0.8, "an" => 0.9, "nn" => 0.9 };

my $DSS_PENALTY = { "gt" => 0.0, "gc" => 0.2, "ga" => 1.2, "gg" => 1.2,
		    "gn" => 1.05, "nt" => 1.05, "nn" => 1.1 };
my $ASS_PENALTY = { "ag" => 0.0, "ng" => 1.05, "an" => 1.05, "nn" => 1.1 };
my $GOOD_INTRON = 1.0;

sub get_splice_penalty() {
    my ($dss,$ass) = @_;
    return 1.1 if "$dss $ass" eq "at ac";
    my $dssres = $DSS_PENALTY->{$dss};
    my $assres = $ASS_PENALTY->{$ass};
    $dssres = 1.5 unless defined $dssres;
    $assres = 1.7 unless defined $assres;
    $dssres += $assres;
    return 2.0 if $dssres > 2.0;
    return $dssres;
}

my $MAX_GAPLEN = 6;    # maximal size (in amino acids) of gap in query that we try to close

my %PARAMETER = ("blat_output=s" => \$BLAT_OUTPUT,
		 "overwrite|force_new" => \$FORCE_NEW,
		 "blat_only" => \$BLAT_ONLY,
		 "verbose" => \$VERBOSE,
		 #
		 "partial_target" => \$ALLOW_MISSING_TARGETS,
	   	 "keep_blat_output" => \$KEEP_BLATOUT,
		 "show_blatline" => \$SHOW_BLATLINE,              # undocumented

		 "best_size|min_score=f" => \$MIN_BEST_SCORE,
		 "min_identity=f" => \$MIN_IDENTITY,
		 "min_coverage=f" => \$MIN_COVERAGE,
		 "multiple_results" => \$MULTIPLE_RESULTS,
		 "single_target_hits|chromosome" => \$ONE_TARGET_ONLY,
		 #
	   	 "max_mismatch=i" => \$MAX_MISMATCH,

		 "plusstrand" => \$VALUE_FORMAT{strand}{'+'},
		 "minusstrand" => \$VALUE_FORMAT{strand}{'-'},
		 "verbal_strands" => \$VERBAL_STRANDS,
		 "show=s" => \$HIT_BLESS_LIST,
		 (map { ("show_${_}=s" => \$MATCH_BLESS_LISTS{$_}) } keys %MATCH_BLESS_LISTS),
		 "hide_undef" => \$HIDE_UNDEF,
		 "hide_defaults" => \$HIDE_DEFAULTS,
		 #
	   	 "region_size=i" => \$REGION_SIZE,

		 "blat_bin=s" => \$BLAT_BIN,
		 "blat_params=s" => \$BLAT_PARAMS,
		 "blat_tilesize=i" => \$BLAT_TILESIZE,
		 "blat_score=i" => \$BLAT_MIN_SCORE,
		 "blat_identity=f" => \$BLAT_MIN_IDENTITY,

		 "nw_insert_penalty=f" => \$NW_PROT_INS,
		 "nw_gap_penalty=f" => \$NW_PROT_GAP,
		 "nw_frameshift_penalty=f" => \$NW_FS,
		 "nw_intron_penalty=f" => \$NW_INTRON,
		 "nw_stop_penalty=f" => \$NW_STP,
		 "exhaust_align_size=i" => \$ADDED_INTRON_MAXLEN, 
		 "exhaust_gap_size=i" => \$ADDED_MATCHING_MAXLEN,


	   	 "split_on_fs" => \$SPLIT_ON_FS,
		 "max_assemble_size=i" => \$ASSEMBLE_SIZE,
		 "min_intron_len=i" => \$MIN_INTRON_LEN,
		 "min_dna_coverage=f" => \$MIN_SEQ_COVERAGE,
		 "transtable=i" => \$TR_TABLE,   # undocumented until BLAT supports this
		 "max_move_exon=i" => \$PRED_OVERLAP,
		 "gap_to_close=i" => \$MAX_GAPLEN,
		 "accepted_intron_penalty=f" => \$GOOD_INTRON,

#		 "margins!" => \$ADD_MARGINS,    # removed from 2010-06-21
		 "trans_shift12" => \$SHOW_FRAME_IN_TRANS # undocumented
);


my $GFF_REVERSE=1;
my $MAX_OVERLAP=5;     # how much overlap we allow for hits on different targets
my $MIN_REMAIN=4;      # how many bases have to be left after intron offsets

my $INTERNAL_MIN_MATCH = 15;
my $IFS_MALUS = 1.5;

local $YAML::InlineSeries = 10;

####################### function headers #######################################

sub usage
{
    my $s = shift;
    print STDERR "\n $s" if (defined $s);
    print STDERR "
Scipio (v1.4) - mapping protein sequences to genomes

    usage:
    scipio.pl [<options>] <target> <query> 
   
    <target> is a DNA file
    <query> is a protein file
    both in FASTA format

    General options:
                --blat_output=<filename>   name of BLAT output file (*.psl) to read from (without
                                           this option, BLAT would be started to create a temporary
                                           output file). With this option, if the given file exists, 
                                           BLAT will not be run, instead the given file will be used;
                                           if the file does not exist, BLAT will be started to create
                                           it.
                --overwrite                use this option together with --blat_output to force
                                           overwriting the specified file with a new BLAT run
                --blat_only                only call BLAT, do not actually run Scipio
                --verbose                  show verbose information (including progress)
 
    Options controlling the choice of hits
                --min_score=<value>        minimal score of a hit; in case of partial hits, of the
                                           best-scoring partial hit (the score is the number of 
                                           matches minus mismatches, divided by total query length;
                                           default for min_score is 0.3)
                --min_identity=<value>     minimal identity in any hit in % (default is 90)
                --min_coverage=<value>     minimal portion of the query that a hit must be aligned
                                           to (by either matches or mismatches, after internal gaps
                                           have been subtracted). In partial hits, only the mapped
                                           subcorresponding partial query sequence is considered 
                                           (default is 60)
                --multiple_results         if this option is given, multiple hits for the same query
                                           will be shown when their scores exceed the minimal score;
                                           multiple hits are named <queryname>_(1) etc.
                --single_target_hits       do not compose hits from multiple targets
                 
    Options controlling the output format:
                --plusstrand=<value>       value shown for plus strand (default: '+')
                --minusstrand=<value>      value shown for minus strand (default: '-'; some YAML
                                           parsers don't like this)
                --verbal_strands           short for --plusstrand=forward --minusstrand=backward
                --show=<list>              which keys are to be shown in YAML output
                --show_intron=<list>       (for details see documentation)
                --show_exon=<list>
                --show_gap=<list>
                --hide_undef               hide keys with undefined or empty list values
                --hide_defaults            hide some keys that have default values
 
    Options passed to BLAT (ignored when BLAT is not run):
                --blat_bin=<name>          name of BLAT executable, defaults to \"blat\"
                --blat_params=<params>     parameters passed to BLAT; see BLAT documentation
                --blat_tilesize=<..>       value for BLAT parameter -tileSize
                --blat_score=<..>          value for BLAT parameter -minScore
                --blat_identity=<..>       value for BLAT parameter -minIdentity; (default is 
                                           0.9*min_identity=81)

    For more options see documentation.

 ";
    exit 1;
}

sub run_bug;                # die on bug; continue if --continue_on_bug is set
sub revcom;                 # returns reverse complement of $1
# sub extract_offset;         # if $1 is <target>:<offset>, strips and returns offset, otherwise 0
sub with_number;            # returns ($1 $2), and adds "s" if $1 is a number not equal to 1
sub diff_str;               # compares strings $1, $2, in $3 positions, starting with position $4 aligned to $5
                            # returns "|" for matches, "X" for mismatches, and " " if "X" or "-" is found
sub str_diff_list;          # returns a list of all positions where diff_str finds mismatches
sub str_diff_count_all;        # counts the "X"s and spaces in diff_str's result
sub str_diff_count_undeterm;   # counts the spaces in diff_str's result
sub str_diff_count_mismatches; # counts the "X"s in diff_str's result
sub str_pos;                # return the occurences of pattern $1 in str

sub extend_to_full_name;    # returns full name of a target (the %full_name value given the cut targetname as key)
sub save_target_region;     # marks a region around a BLAT hit. ($1,$2,$3) = targetname, from, to
sub cut_target_regions;     # cuts regions to proposed length 
sub fill_target_region;     # fills a region with the sequence. ($1,$2,$3) = target, ref to seq, offset
sub read_targets;           # read target sequences from FASTA file, fills marked regions
sub validate_targets;       # removes empty regions, check if all sequences have been found
sub nu_to_aa;               # division by 3; rounded

# revtranslate_to_pattern ($seq)
# return a regular expression based on reverse translation of the amino acid sequence $seq.
# example:  revtranslate_to_pattern("AS")
# returns   "(gct|gcc|gca|gcg)(tct|tcc|tca|tcg|agt|agc)"
sub revtranslate_to_pattern; 

# split_codons($aa, $dna, $backward)
# return a list of patterns modelling a splice site (dss if $backward is true, ass otherwise)
# $aa contains one amino acid, $dna the other splice site (which is known), together with 
# one half (possibly empty) of a codon. The returned list contains all patterns that complete
# the known splice site to a sequence consisting of a codon coding for $aa,
# with an intron inserted into it
# if $backward is true, return a dss (always ending "g[ct]") for an ass given in $dna
# if $backward is false, return an ass (always starting "ag") for a dss given in $dna
# example:  split_codons("A", "gcgt",0)
# returns   ('agt','aggct','agc','aggcc','aga','aggca','agg','aggcg')
sub split_codons;

# get_splicesite_patterns_forward($suffix, $dna)
# take the first two amino acids of $suffix; each of them is split with split_codons
# $dna is representing the dna sequence starting with the first possible splice point
# return a pattern that matches all possible ass's that complete $dna to a valid intron
sub get_splicesite_patterns_forward;
sub get_splicesite_patterns_backward;

# subseq ($targetname, $from, $count, $complement)
# returns subsequence of a target; 
# if $complement is true, $from must be negative to be inside the target, and -$from 
# refers to the subsequence end (so $count characters are taken beginning with -$from-$count);
# in this case the reverse complement is returned; unknown bases are given as ".", whereas
# positions outside the target are returned as spaces. It is guaranteed that the length
# of the result equals $count. If $count is negative, undef is returned.
sub subseq;    
  
# get_intron_offset ($gap_aa, $dss_seq, $ass_seq, $frameshift_mode, preferred_offset) 
# this routine concatenates a prefix of $dss_seq to a suffix of $ass_seq, where
# $dss_seq and $ass_seq are nucleotide sequences with length of 3*length($gap_aa)+2,
# and tries to align the resulting translated sequence to $gap_aa.
# Returned as offset is the length of the left part (the prefix). In frameshift mode,
# only multiples of three can be taken as offsets.
# If multiple offsets yield the minimum number of errors, then the one closest to
# preferred_offset is taken
#
# it returns a list with four values: 
#  - best offset, 
#  - best number of mismatches
#  - translated sequence
#  - did we find a splice site pattern
#
# if frameshift_mode is false, then splice side patterns are traded off against more
# errors (INTRON mode); in frameshift mode, splice side patterns are ignored
# 
sub get_intron_offset;

# calculate_diff_strings ($targetname, $queryname, 
#                         \@targetlocations, \@querylocations, $complement);
# returns a list of mismatch strings, one for each exon
sub calculate_diff_strings;

# check if successive BLAT hits are compatible to each other
sub incompatible;

# replace_in_list($from, $to, $list, $newelem)
# replace elements x with (from < x <= to) from @$list
# by elements found in @$newelems
# returns the change in number of elements
#
# used in replace_in_last/next
sub replace_in_list;

# replace_in_last( $prevhit, $dnaseq, $protseq, $dna_add, $nsub)
# replaces the last $nsub nucleotides of the last exon found in $prevhit
# by $dnaseq, and changes the coordinates of the hit accordingly
#
# NOTE: - if $dnaseq is nonempty, the hit MUST end on a codon boundary after
#       subtracting $nsub codons
#       - the hit MUST not have a seqshift in the last $nsub nucleotides
#
# the exon matching is changed the following way:
# - the translation is set to the translation of $dnaseq.$dna_add 
#   if the translation is longer then $protseq, only the prefix of the corresponding
#   length will be compared to $protseq;
# - the mismatch positions between $protseq and the translation are added to
#   the last matching
# - translation is concatenated to the last matching
# - prot_end, nucl_end and dna_end are increased by the length of $protseq
sub replace_in_last;

# replace_in_next( $nexthit, $dnaseq, $protseq, $dna_add, $nsub)
# replaces the first $nsub nucleotides of the first exon found in $nexthit 
# by $dnaseq,  and changes the coordinates of the hit accordingly
#
# NOTE: - $dnaseq MUST to start at the reading frame specified by the
#       length of $dnaseq
#
# the exon matching is changed the following way:
# - the translation is set to the translation of $dna_add.$dnaseq;
#   if the translation is longer then $protseq, only the suffix of the corresponding length 
#   will be compared to $protseq; here, if dnaseq doesnt start in frame 0, 
# - the mismatch positions between $protseq and the translation are added to
#   the first matching
# - translation is concatenated with first matching
# - prot_start, nucl_start and dna_start are decreased by the length of $protseq
sub replace_in_next;

# add initial ($hit, $type, $qlen)
# inserts an initial matching of type $type to the hit
sub add_initial;

# add_final ($hit, $type, $qlen)
# inserts a terminal matching of type $type to the hit
sub add_final;

# set_problem_field ($hit)
# sets the problem and reason keys according to the condition of the hit
sub set_problem_field;

#
# nw_align ($dnaseq, $protseq)
# find the optimal alignment of the protein sequence versus the dna
# this is used to fill small gaps
sub nw_align;

##################################### main #####################################

### Command Line
&usage unless(@ARGV);
my $commandline = join(" ",grep {/^--/ && !/^--blat_o/ && !/^--verbose/ && !/^--ov/} @ARGV);

&GetOptions (%PARAMETER) or &usage;

my %dbase;
my %queries;
my %full_name;
my $ttable = Bio::Tools::CodonTable->new( -id => $TR_TABLE );

# prepare verbose output
my $comments;
if ($VERBOSE) {
    open COMMENT, ">&", \*STDERR;
    select COMMENT;
    $| = 1;
    select STDOUT; 
} else {
    open COMMENT, '>', \$comments;
}

$HITBUFFER_SIZE=$ASSEMBLE_SIZE if ($ASSEMBLE_SIZE > $HITBUFFER_SIZE);
$HITBUFFER_SIZE+=$REGION_SIZE;

unless (defined $KEEP_BLATOUT) {
    $KEEP_BLATOUT = ($BLAT_OUTPUT ne "");
}

unless (defined $MIN_IDENTITY) {
    $MIN_IDENTITY = $DEFAULT_MIN_IDENTITY / 100;
} elsif ($MIN_IDENTITY <= 1) {
    print STDERR "Warning: multiplied min_identity by 100 (should be specified in % now)\n";
} else {
    $MIN_IDENTITY /= 100;
}
if ($MIN_COVERAGE > 1) {
    $MIN_COVERAGE /= 100;
}


if ($MIN_SEQ_COVERAGE > 0.2) {
    print STDERR "Warning: min_seq_coverage>0.2 will not give results. Ignored.\n";
    $MIN_SEQ_COVERAGE = 0;
} else {
    $MIN_SEQ_COVERAGE /= 100;
}

# if ($SHOWALL_SCORE < $MIN_BEST_SCORE) {
#     $MIN_BEST_SCORE = $SHOWALL_SCORE;
# }
    
if ($VERBAL_STRANDS) {
    $VALUE_FORMAT{strand} = { '+' => "forward", '-' => "backward" };
}

unless (defined $ADDED_MATCHING_MAXLEN) {
    $ADDED_MATCHING_MAXLEN = defined $BLAT_TILESIZE ? 3 * $BLAT_TILESIZE : 15;
} 

# show parameters in verbose mode
while (my ($key, $ref) = each %PARAMETER) {
    $key =~ s/[!:=].*$//;
    print COMMENT "$key=$$ref\n" if (defined $$ref);
}
print COMMENT "\n";
# parameters passed to BLAT
$BLAT_MIN_IDENTITY = 90 * $MIN_IDENTITY unless (defined $BLAT_MIN_IDENTITY);
$BLAT_PARAMS.=" -tileSize=$BLAT_TILESIZE" if (defined $BLAT_TILESIZE && $BLAT_PARAMS !~ /tileSize/);
$BLAT_PARAMS.=" -minIdentity=$BLAT_MIN_IDENTITY" unless ($BLAT_PARAMS =~ /minIdentity/);
$BLAT_PARAMS.=" -minScore=$BLAT_MIN_SCORE" unless ($BLAT_PARAMS =~ /minScore/);

my ($dbasefile, $queryfile) = @ARGV;

&usage ("Please specify target file in FASTA format!\n") unless ($dbasefile && -f $dbasefile);
&usage ("Please specify query file in FASTA format!\n") unless ($queryfile && -f $queryfile);


### Open query file and store sequences in hash
my $querySeq = Bio::SeqIO->new(-file=>$queryfile, -format=>"fasta") 
    or usage("Queryfile did not contain sequences.\n");
my $unique_queries=1;
my $queries_with_problems=0;
while (my $seq = $querySeq->next_seq()) {
    my ($newid, $newseq) = ($seq->id(), $seq->seq());
    unless (defined $newseq && defined $newid) {
	print STDERR "ERROR: Could not read query sequence".(defined $newid ? " '$newid'" : "").".\n";
	next;
    }
    $newid =~ s/_\(\d+\)$//;
    $newseq =~ s/\*$//;               # strip final stop codon from sequence
    $newseq =~ tr/a-z/A-Z/;           # make seqence upper case
    if ($newseq =~ s/[-]//g || $newseq =~ /\*/) {
	$queries_with_problems=1;
    }
    my $oldseq = $queries{$newid};
    if (defined $oldseq && $oldseq ne $newseq) {
	$unique_queries=0;
	my $len = length($newseq);
	$newid = "${newid}_[$len]";
	if (length($oldseq) == $len || 
	    defined $queries{$newid} && $queries{$newid} ne $newseq) {
	    print STDERR "Neither query sequence names nor sequence lengths unique. Giving up.\n";
	    exit 1;
	}
    } 
    $queries{$newid} = $newseq;
}

print STDERR "Warning: Query sequence names not unique. Trying to identify by sequence length...\n"
    unless ($unique_queries);
    

### Determine BLAT output file and run BLAT
unless ($FORCE_NEW) {
    if ($BLAT_OUTPUT && -f $BLAT_OUTPUT) {
	print COMMENT "Not running BLAT. Using previously created output file \"$BLAT_OUTPUT\".\n";
    } else {
	$FORCE_NEW = 1;
    }
}

if ($queries_with_problems) {
    print STDERR "Warning: Query sequences contain gaps and/or premature stopcodons.\n"
	."It will be necessary to remove them if blat is run separately!\n";
    if ($FORCE_NEW) {
	$queryfile = "$queryfile.temp";
	open TEMPFILE, ">$queryfile";
	while (my ($qname, $qseq) = each %queries) {
	    $qseq =~ s/\*/X/g; # Replace '*' by 'X' for BLAT
	    print TEMPFILE ">$qname\n$qseq\n";
	}
    }
}
	
if ($FORCE_NEW) {
    # do run BLAT (output file not specified, or not existing)
    unless ($BLAT_OUTPUT) {
	my ($n,$next)=("","a"); 
	do {  
	      $BLAT_OUTPUT="Scipio${$}${n}_blat.psl";
	      $n=$next; $next++;
	} while (-f $BLAT_OUTPUT);
    }
    my $showtime = time(); 
    print COMMENT 
	"Running BLAT to produce file \"$BLAT_OUTPUT\": (started at ".
	join (":", map { sprintf "%02d",$_ } (localtime($showtime))[2,1,0]).")...\n";
    my $BLAT_EXESTRING = "$BLAT_BIN -t=dnax -q=prot -noHead $BLAT_PARAMS \"$dbasefile\" \"$queryfile\" \"$BLAT_OUTPUT\"";
    print STDERR "$BLAT_EXESTRING\n" if ($SHOW_BLATLINE);
    open (BLRUN, "$BLAT_EXESTRING |");

    print COMMENT while (<BLRUN>);
    my $blat_exitstatus = $? >> 8;
    printf STDERR "Warning: blat returned exit status %d\n", $blat_exitstatus 
	if ($blat_exitstatus);
    $showtime = time() - $showtime;
    print COMMENT "...done in $showtime seconds\n";
    if ($queries_with_problems) {
	# delete temporary query file
	unlink $queryfile;


	$queryfile =~ s/\.temp$//;
    }
}
if ($BLAT_ONLY) {
    exit(-1);
}
open BLOUT, $BLAT_OUTPUT or die ("Could not access output file \"$BLAT_OUTPUT\"\n");


my $hitcount=0;
### Open target file and load sequence parts around hits
while (<BLOUT>)
{
    next if (/^\#/ || /^$/);
    my @result = split "\t";
    next unless @result>=21;
    my ($targetname,$targetsize,$targetstart,$targetend) = @result[13..16];
    next if (($targetstart.$targetsize.$targetend) =~ /\D/);
#    my $offset = extract_offset($targetname);
    $targetstart -= $HITBUFFER_SIZE; 
    $targetstart = 0 if ($targetstart<0);
    $targetend += $HITBUFFER_SIZE;
    &save_target_region($targetname, $targetstart, $targetend, $targetsize);
    $hitcount = $.;
}
&cut_target_regions();

print COMMENT 
    "We have to check ".with_number($hitcount," BLAT hit")." in ".
    with_number(scalar keys %dbase," target sequence").".\n";
print COMMENT "(Rerun BLAT with other options to reduce / increase this number).\n\n";
print COMMENT "Reading sequences from \"$dbasefile\":\n";

seek BLOUT, 0,0;
$.=0;
&read_targets($dbasefile);
&validate_targets;

my %hits = map { $_ => [] } (keys %queries);

###################### postprocessing of BLAT hits ######################
my $dot_stands_for = max (1,$hitcount / 100);

print COMMENT "Processing BLAT hits:\n";
BLATLINES: while(<BLOUT>)
{
    next if (/^\#/ || /^$/);  # skip comments and empty lines
    print COMMENT "." unless ($. % $dot_stands_for);
    printf COMMENT "%2.0f%%", (($.) / $hitcount)*100  unless ($. % ($dot_stands_for * 10));
    printf COMMENT "\n" unless ($. % ($dot_stands_for * 50));

    # @result[..]=
    # 0:  number of matches
    # 1:  number of mismatches
    # 2:  number of repeated matches
    # 3:  number of unspecified nucleotides
    # 4:  number of gaps in query (not including zero size gaps)
    # 5:  total size of gaps in query
    # 6:  number of gaps in target
    # 7:  total size of gaps in target
    # 8:  strand
    # 9:  name of query
    # 10: total length of query
    # 11: start of aligned part of query 
    # 12: end of aligned part of query  
    # 13: name of target
    # 14: size of complete target
    # 15: start of aligned part of target
    # 16: end of aligned part of target
    # 17: number of blocks
    # 18: size of blocks in codons/aas, comma separated
    # 19: start positions of blocks in query
    # 20: start positions of blocks in target

    my @result = split "\t";
    next unless @result>=21;  # skip lines with less than 21 entries

    # read line of BLAT output file
    # the prefix 'q' is used for amino acid positions in the query
    my ($BLAT_matchcount, $BLAT_mismatchcount, $gapcount, $strand, $queryname, $querysize, 
	$targetname, $targetsize, $block_sizes, $qFroms, $tFroms) 
	= @result[0,1,5,8,9,10,
		  13,14,18..20];


    my ($undetermined, $insertions) = (0,0);
    &extend_to_full_name($targetname);
    next unless (exists  $dbase{$targetname});

    my $queryseq = $queries{$queryname};
    if (defined $queryseq && $querysize != length($queryseq))
    {
	if ($unique_queries) {
	    print STDERR "Warning: query length mismatch. This will produce unpredictable results!\n";
	    $unique_queries = 0;
	} else {
	    $queryname = "${queryname}_[$querysize]";
	    $queryseq = $queries{$queryname};
	}
    }	
    unless (defined $queryseq)
    {
	print STDERR "No query sequence '$queryname' of length $querysize found.\nSkipped.\n";
	next;
    }
    $strand =~ s/^.//;
    my $complement = ($strand =~ /-/);
 
    chomp $tFroms;
    s/\,$// foreach ($block_sizes, $qFroms, $tFroms);
    my @block_sizes = split ",",$block_sizes;
    next if (($BLAT_mismatchcount >= ($MAX_MISMATCH + @block_sizes) && $MAX_MISMATCH) ||
	     ($BLAT_mismatchcount-1) / ($BLAT_matchcount + $BLAT_mismatchcount) > 1-$MIN_IDENTITY);
    my @qFroms = split ",", $qFroms;
    my @tFroms = split ",", $tFroms;
    my @qTos=(); 
    my @tTos=();
    foreach (0..$#block_sizes) {
	$tFroms[$_] += $complement ? -$targetsize : 0;
	push @qTos, $qFroms[$_]+($block_sizes[$_]);
	push @tTos, $tFroms[$_]+($block_sizes[$_]*3);
    }
 

    # ignore if internal partial match is too small
    next BLATLINES if ($qFroms[0] > 5 && $qTos[-1] < $querysize-5 && $BLAT_matchcount < $INTERNAL_MIN_MATCH);


    ### Use Needleman/Wunsch to add small missing exons
    ### (this is SLOW and should not be tried for long introns)
    foreach (my $inno=0; $inno<$#block_sizes; $inno++) {
	my $qgaplen = $qFroms[$inno+1] - $qTos[$inno];
	my $intronlen = $tFroms[$inno+1] - $tTos[$inno];
	my (@tFromsAdded, @tTosAdded, 
	    @qFromsAdded, @qTosAdded, @added_sizes);

        # Short dna parts used not to be examined here. Why?? 
	# Skip N/W if number of residues exceeds the gap size
	next if ($qgaplen*3 > $intronlen || $qgaplen > $ADDED_MATCHING_MAXLEN); # || $intronlen < $MIN_INTRON_LEN);
	if ($intronlen < $ADDED_INTRON_MAXLEN) {
	    # try to align missing part with Needleman/Wunsch

	    my $left_overlap = max(0,min($PRED_OVERLAP,$block_sizes[$inno]-3));
	    my $right_overlap = max(0,min($PRED_OVERLAP,$block_sizes[$inno+1]-3));
	    my $both_overlaps = $left_overlap + $right_overlap;
	    my $gapstr=substr($queryseq, $qTos[$inno]-$left_overlap, $qgaplen+$both_overlaps);
	    my $intronstr=&subseq($targetname, $tTos[$inno]-3*$left_overlap, $intronlen+3*$both_overlaps, $complement); 

	    my @added_exons = &nw_align($intronstr, $gapstr);
	    my ($cutleft, $cutright) = (0,0);
	    if (@added_exons > 1) {
		$cutleft = $left_overlap - (shift @added_exons)->[2];
		$cutright = $right_overlap - (pop @added_exons)->[2];
		foreach (@added_exons) {
		    my ($tpos, $qpos, $nextsize, $nextscore) = @$_;
		    $tpos += ($tTos[$inno] - 3 * $left_overlap);
		    $qpos += ($qTos[$inno] - $left_overlap);
		    push @tFromsAdded, $tpos;
		    push @qFromsAdded, $qpos;
		    push @tTosAdded, $tpos + 3*$nextsize;
		    push @qTosAdded, $qpos + $nextsize;
		    push @added_sizes, $nextsize;
		}
		$tTos[$inno] -= 3*$cutleft;
		$_->[$inno] -= $cutleft foreach (\@qTos, \@block_sizes);
		$tFroms[$inno + 1] += 3*$cutright;
		$qFroms[$inno + 1] += $cutright;
		$block_sizes[$inno + 1] -= $cutright;
		$gapcount += ($cutleft+$cutright);
	    } 	    
	} elsif ($qgaplen >= 4) {
	    # try to find exact hits in missing part that are close to neighbouring exons
	    foreach ( [0,0], [1,0], [0,1], [1,1], [2,1], [1,2] ) {
		my $patlen = $qgaplen - sum(@$_);
		next if $patlen < 3 || $patlen < $qgaplen/2
		    || $patlen*3 > $intronlen;
		my ($left_cutoff, $right_cutoff) = @$_;
		my $gapstr=substr($queryseq, $qTos[$inno] + $left_cutoff, $patlen);
		my $pattern = &revtranslate_to_pattern($gapstr);
		if ($pattern) {
		    my $checkrange = 3*$patlen + $MIN_INTRON_LEN;
		    my $intFrom = $tTos[$inno];
		    my $location;
		    if ($intronlen < $checkrange) {
			$checkrange = $intronlen;
		    }
		    if (&subseq($targetname, $intFrom, $checkrange, $complement) =~ /$pattern/) {
			$location = $-[0];
		    } elsif ($checkrange < $intronlen) {
			$intFrom = $intFrom + $intronlen - $checkrange;
			if (&subseq($targetname, $intFrom, $checkrange, $complement) =~ /$pattern/) {
			    $location = $-[0];
			}
		    }
		    next unless defined $location;

#		    $gapcount -= $patlen;  now taken care of below!
		    @tFromsAdded = ($intFrom + $location);
		    @tTosAdded = ($intFrom + $location + 3*$patlen);
		    @qFromsAdded = ($qTos[$inno] + $left_cutoff);
		    @qTosAdded = ($qFroms[$inno+1] - $right_cutoff);
		    @added_sizes = ($patlen);
		    last;
		}
	    }
	}
	if (@tFromsAdded) {
	    splice(@tFroms,      $inno+1, 0, @tFromsAdded);
	    splice(@tTos,        $inno+1, 0, @tTosAdded);
	    splice(@qFroms,      $inno+1, 0, @qFromsAdded);
	    splice(@qTos,        $inno+1, 0, @qTosAdded);
	    splice(@block_sizes, $inno+1, 0, @added_sizes);
	    $gapcount -= sum (@added_sizes);
	    $inno += scalar @tFromsAdded;
	}
    }

    ### Process Introns
    my ($querystart, $targetstart) = ($qFroms[0], $tFroms[0]);
    my @targetlocations = ($targetstart);
    my @querylocations = ($querystart*3);
    my @types = ("exon");
    
    for (my $inno=0; $inno<$#block_sizes; $inno++) {
	my ($tfrom, $tto, $qfrom, $qto) = ($tTos[$inno], $tFroms[$inno+1], 
					   $qTos[$inno], $qFroms[$inno+1]);

	# enlarge overlap if borders fit
	if ($tto < $tfrom) {
	    while ($ttable->translate(&subseq($targetname, $tfrom, 3, $complement)) eq substr($queryseq, $qfrom, 1)) {
		$tfrom += 3;
		$qfrom++;
		$gapcount--;
	    }
	    while ($ttable->translate(&subseq($targetname, $tto-3, 3, $complement)) eq substr($queryseq, $qto-1, 1)) {
		$tto -= 3;
		$qto--;
		$gapcount--;
	    }
	}

	# the following two rare cases may happen if we are very unlucky;
	# so deal with it to avoid problems
	if ($tfrom +3 > $tTos[$inno+1]) {
	    # remove the next match if it is contained completely in the current one
	    if ($inno == $#block_sizes-1) {
		$gapcount -= ($qto-$qfrom);
		splice (@$_, -1) foreach (\@tTos, \@qTos); # for query-/targetend
		last;
	    }
	    $gapcount += $qTos[$inno+1] - $qto;
	    $inno++;
	    ($tto, $qto) = map { $_->[$inno+1] } (\@tFroms, \@qFroms);
	} elsif ($targetlocations[-1] > $tto-3) {
	    # remove current match if it is contained completely in the next
	    if (@targetlocations == 1) { # first match is deleted
		$gapcount -= ($qto-$qfrom);
		@targetlocations = ($targetstart = $tto);
		@querylocations = (($querystart = $qto)*3);
		next;
	    } else {
		$gapcount += $qfrom - int($querylocations[-1]/3);
		my $mintto = $targetlocations[-3]+3;
		if ($mintto > $tto) {
		    # doubtful if that ever happens...
		    my $cutoff = $mintto -$tto;
		    $tto = $targetlocations[-3] + (-$cutoff) % 3;
		    $_ += ($cutoff + 2)/3 foreach ($qto, $gapcount);
		}
		($tfrom, $qfrom) = map { $_->[-2] } (\@targetlocations, \@querylocations);
		splice (@$_, -2, 2) foreach (\@targetlocations, \@querylocations, \@types);
		my $frame = $qfrom%3;
		$_ -= $frame foreach ($tfrom, $qfrom);
		$qfrom/=3;
	    } 
	}

	my $qgaplen=($qto-$qfrom);
	my $offset = 0;
	my $type = "gap";
	my $additional = $tto-$tfrom-$qgaplen*3;
	if ($qgaplen <= $MAX_GAPLEN || $tfrom > $tto) {
	    ### DEBUG PART
#	    print COMMENT "WARNING. Have query overlaps.\n" if ($qgaplen < 0);
            ### END OF DEBUG PART
	    if ($additional == 0) {
		$gapcount -= $qgaplen; next; # gap was closed completely: no location needed
	    } elsif ($additional < 0) {
		### case 1: no additional bases for an intron
		### here we just fill up optimally as much of target sequence
		### as we can; this is always the case when qgaplen > MAX_GAPLEN
		$type = "seqshift" if ($qgaplen <= $MAX_GAPLEN);

		# discard overlapping parts
		my $tframe = ($tto-$tfrom) % 3;
		if ($tfrom > $tto) { 
		    my $tgap = $tfrom - $tto + $tframe;
		    $tto += $tgap;
		    $tfrom -= $tgap;
		    $qto += $tgap/3;
		    $qfrom -= $tgap/3;
		    $qgaplen = $qto-$qfrom;
		    $gapcount += $tgap/3*2;
		}
		# find optimal location for inserting frameshift
		my $tseq = &subseq($targetname, $tfrom, $tto-$tfrom, $complement);
		my $closed_gaplen = $tto-$tfrom-$tframe;
		my ($max_mismatches, $best_gappos) = ($qgaplen, 0);
#		for (my $tgappos=0; $tgappos <= $closed_gaplen; $tgappos+=3) {
		for (my $tgappos=$closed_gaplen; $tgappos >= 0; $tgappos-=3) {
		    my $gap_trans = $ttable->translate(
			substr($tseq,0,$tgappos).substr($tseq,$tgappos+$tframe));
		    my $rightlen = ($closed_gaplen-$tgappos)/3;
		    my $querypart = substr($queryseq,$qfrom,$tgappos/3).substr($queryseq,$qto-$rightlen,$rightlen);
		    my $new_mismatches = &str_diff_count_mismatches($gap_trans, $querypart);
		    if ($new_mismatches < $max_mismatches) { 
			($max_mismatches, $best_gappos) = ($new_mismatches, $tgappos);
			last if ($new_mismatches == 0);
		    }
		}
		# fill up to optimal location
		$tfrom += $best_gappos;
		$qfrom += $best_gappos/3;
		$tto = $tfrom + $tframe;
		$qto -= ($closed_gaplen-$best_gappos)/3;
		$gapcount -= $closed_gaplen/3;
	    } else { 
                ### case 2: we have a gap of nucleotides 
                ### this is an intron or additional nucleotides (e.g. frameshift)

		# calculate marginal parts to be checked; always check matching regions
		my $max_overlap = int($qto - $querylocations[-1]/3)-1;
		my $left_overlap = max(0, $max_overlap - $qgaplen);
		if ($PRED_OVERLAP <= $left_overlap) { 
		    for (my $i = 1; $i < $max_overlap; $i++) {
			if ($ttable->translate(&subseq($targetname, $tto - 3*$i, 3, $complement)) 
			    ne substr($queryseq, $qto-$i,1)) {
			    $left_overlap = max($i-$qgaplen, $PRED_OVERLAP);
			    last;
			}
		    }
		}
		my $right_overlap = $block_sizes[$inno+1]-1;
		$max_overlap = $right_overlap + $qgaplen;
		if ($PRED_OVERLAP <= $right_overlap) { 
		    for (my $i = 0; $i < $max_overlap-1; $i++) {
			if ($ttable->translate(&subseq($targetname, $tfrom + 3*$i, 3, $complement)) 
			    ne substr($queryseq, $qfrom+$i,1)) {
			    $right_overlap = max($i-$qgaplen+1, $PRED_OVERLAP);
			    last;
			}
		    }
		}
		my $gapstart=($qfrom-$left_overlap);
		my $qgaplenplus=$qgaplen + $left_overlap + $right_overlap;
		### DEBUG PART
		if ($qgaplenplus < 0) {
		    print COMMENT "ERROR. There should not be query overlaps. Skipped this hit (line $.)\n";
		    next BLATLINES;
		}
		### END OF DEBUG PART

		# find optimal location for inserting intron
		my ($dss_seq, $ass_seq, $gap_aa) = 
		    (&subseq($targetname, $tfrom-3*$left_overlap, $qgaplenplus*3+2, $complement),
		     &subseq($targetname, $tto-($qgaplenplus-$right_overlap)*3-2, 
			     $qgaplenplus*3+2, $complement),
		     substr($queryseq."*", $gapstart, $qgaplenplus));
		### DEBUG PART
		if ($gap_aa eq "*") {
		    &run_bug("Internal error. This shouldn't happen!");
		}
		### END OF DEBUG PART
		my ($get_offset, $mismatchadd, $gap_trans, $splicesite_score) = 
		    &get_intron_offset($gap_aa, $dss_seq, $ass_seq,
				       $additional < $MIN_INTRON_LEN,
				       3*($qgaplenplus-$right_overlap)   );
		$get_offset -= 3*$left_overlap;
		$mismatchadd=$qgaplen if ($qgaplen < $mismatchadd);
		
		# close the gap in query, if possible with few mismatches; otherwise keep the gap
		# (short gaps less than MIN_INTRON_LEN are always closed)
		if (($mismatchadd <= 4  || $tto-$tfrom < $MIN_INTRON_LEN) && $tfrom+$get_offset >= $targetlocations[-1]+2)
		{   
		    $qto=$qfrom;
		    $tto -= $qgaplen*3;
		    $gapcount -= $qgaplen; 
		    $offset = $get_offset;
		    $type = $additional < $MIN_INTRON_LEN ? "seqshift" : $splicesite_score>2.0-$GOOD_INTRON ? "intron" : "intron?";
		}
	    }
	    
	}
	push @types, ($type, "exon"); 
	push @targetlocations, ($tfrom+$offset, $tto+$offset);
	push @querylocations, ($qfrom*3+$offset, $qto*3+$offset);
    } # end process introns

    my ($queryend, $targetend) = ($qTos[-1], $tTos[-1]);

    ### Try to add suffix that wasn't found by BLAT
    # adding more matches to last exon
    if ($queryend < $querysize && $queryend >= $querysize - $MAX_GAPLEN
	&& ($queryend-$PRED_OVERLAP)*3 > $querylocations[-1]) {
	$queryend -= $PRED_OVERLAP;
	$targetend -= 3*$PRED_OVERLAP;
    }
    my $targetbound = ($complement ? 0 : $targetsize)-3*$PRED_OVERLAP;
    while ($queryend < $querysize && $targetend <= $targetbound-3) {
	my $aa = substr($queryseq, $queryend, 1);
	my $codon = &subseq($targetname, $targetend, 3, $complement);
	last unless ($ttable->translate($codon) eq $aa);
	$queryend++; $targetend+=3;
    }
    # separate new exon
    my $qgaplen = $querysize - $queryend;
    if (0 < $qgaplen && $qgaplen <= $MAX_GAPLEN) # && $ADD_MARGINS)
    {
	my $suffix = substr($queryseq, $queryend-1)."*";
	my $dss =  &subseq($targetname, $targetend-3, 7, $complement);
	my $pattern = &get_splicesite_patterns_forward($suffix,$dss);
	my $target = &subseq($targetname, $targetend, $ASSEMBLE_SIZE, $complement);
	my ($tfrom, $tto, $prev_qnto, $qnfrom, $type);
	if ($pattern && $target =~ /$pattern/)
	{
	    ($tfrom, $tto) = map { $_ += $targetend } ($-[0]+2, $+[0]-3);
	    $qnfrom = $prev_qnto = $querysize*3 - ($tto-$tfrom);
	    $targetend += $qgaplen*3 - ($tto-$tfrom);
	    $type = "intron";
	} elsif ($ttable->translate(substr($target, 0, 3*$MAX_GAPLEN+3)) =~ /(\*)/) {
	    my $added = $-[0];
	    if ($added == $qgaplen) {
		$queryend = $querysize;
		$targetend += $qgaplen*3;
	    } else {
		$prev_qnto = $queryend*3;
		$type = "seqshift";
		$tto = $targetend + 3*$added;
		$tfrom = max ($targetend, $tto - 3*$qgaplen);
		$qnfrom = max (3*($querysize-$added), $prev_qnto);
		$gapcount += max(0, $qgaplen-$added);
	    }
	}		
	if (defined $type) {
	    push @querylocations, ($prev_qnto, $qnfrom);
	    push @targetlocations, ($targetend, $tfrom);
	    push @types, ($type, "exon");
	    $targetend = $tto;
	    $queryend = $querysize;
	} 
    }

    ### Try to add prefix that wasn't found by BLAT
    # adding more matches to first exon
    if ($querystart > 0 && $querystart <= $MAX_GAPLEN
	&& ($querystart + $PRED_OVERLAP) < $queryend
	&& (@querylocations < 2 || ($querystart + $PRED_OVERLAP)*3 < $querylocations[1])) {
	$querystart += $PRED_OVERLAP;
	$_ += 3*$PRED_OVERLAP foreach ($targetstart, $targetlocations[0], $querylocations[0]);
    }
    $targetbound -= ($targetsize-6*$PRED_OVERLAP);
    while ($querystart > 0 && $targetstart >= $targetbound) {
	my $aa = substr($queryseq, $querystart-1, 1);
	my $codon = &subseq($targetname, $targetstart-3, 3, $complement);
	last unless ($ttable->translate($codon) eq $aa);
	$querystart--;
	$_-=3 foreach ($targetstart, $targetlocations[0], $querylocations[0]);
    }
    # separate new exon
    if (0 < $querystart && $querystart <= $MAX_GAPLEN) { # && $ADD_MARGINS)
	my $ass = &subseq($targetname, $targetstart-4, 7, $complement);
	my $prefix = substr($queryseq, 0, $querystart+1);
	my $searchstart = $targetstart - $ASSEMBLE_SIZE;
	$searchstart = max( ($complement ? -$targetsize : 0), $searchstart);
	my $target = &subseq($targetname, $searchstart, $targetstart-$searchstart, $complement);
	my $pattern = &get_splicesite_patterns_backward($prefix, $ass);
	my ($tfrom, $tto, $qnto, $type);
	if ($pattern && $target =~ /^.*($pattern)/) {
	    ($targetstart, $tto) = map {$_ + $searchstart } ($-[1], $+[1]-2);
	    $querylocations[0] = $qnto = $tto - $targetstart;
	    $targetlocations[0] += ($qnto - 3*$querystart);
	    $type = "intron";
	} else {
	    # translate target up to MAX_GAPLEN before matchstart
	    # revert it for technical reasons
	    $target = reverse $ttable->translate(substr($target, -3*$MAX_GAPLEN));
	    $target =~ s/\*.*$//g; # remove everything from first stop codon
	    # qfroms contain candidates for new matchstarts (having M)
	    # given in distance from query start (0 meaning an M fitting without gaps)
	    my @qfroms = ();
	    push @qfroms, $+[0]-$querystart while ($target =~ /M/g);
	    my $qfrom = reduce { abs($a) < abs($b) ? $a : $b } @qfroms;
	    if (defined $qfrom) {
		if ($qfrom == 0) {
		    $_ -= $querylocations[0] foreach ($targetstart, $targetlocations[0]);
		    $querystart = $querylocations[0] = 0;
		} elsif ($qfrom > 0) {
		    $tto = $targetstart - $qfrom*3;
		    $targetstart -= ($qfrom + $querystart)*3;
		    $qnto = $querystart*3;
		    $type = "seqshift";
		} else {
		    $tto = $targetstart;
		    $qnto = ($querystart+$qfrom)*3;
		    $targetstart -= $qnto;
		    $type = "seqshift";
		    $gapcount -= $qfrom;
		}
	    }
	}
	if (defined $type) {
	    unshift @querylocations, (0, $qnto);
	    unshift @targetlocations, ($targetstart, $tto);
	    unshift @types, ("exon", $type);
	    $querystart = 0;
	}
    }
    push @targetlocations, $targetend;
    push @querylocations, $queryend*3;

    # if last match reaches end of query, check if there is a stop codon
    my $stopcodon;
    if ($queryend == $querysize)
    {
	$stopcodon = &subseq($targetname, $targetend, 3, $complement);
	$stopcodon = undef unless($ttable->is_ter_codon($stopcodon));
    } 

    # matchsize is the length of the aligned part of the query
    # matchsize differs from querysize only when partial hits are found
    my $matchsize = ($queryend-$querystart);
	
    my @diff_strings = 
	&calculate_diff_strings($targetname, $queryname, \@targetlocations, \@querylocations, $complement);
    
    ### debug part
    if (length(join "", map { /:(.*)$/; $1 } @diff_strings) + $gapcount != $matchsize)
    {
	print STDERR "query:$queryname\ntarget:$targetname\n";
	&run_bug("Incorrect calculation of unmatched aa's in line $.!\n");
    }
    ### end of debug part

    my @matchings = ();
    my @mismatchpos = ();
    my $mismatchcount = 0;

    ### Construct matchings hash
    foreach (@types)
    {
	my $previous = $matchings[-1];
	my ($tfrom, $qnfrom) = (shift @targetlocations, shift @querylocations);
	my ($tto, $qnto) =  ($targetlocations[0], $querylocations[0]);
	my ($qfrom, $qto) = map ( &nu_to_aa($_), ($qnfrom, $qnto) );
	my $current = { "${TPFX}_start" => $tfrom,
			"${TPFX}_end" => $tto,
			"${QNPFX}_start" =>  $qnfrom,
			"${QPFX}_start" => $qfrom,
			"${QNPFX}_end" => $qnto,
			"${QPFX}_end" => $qto };
	if (/seqshift/) {
	    my $translation = $ttable->translate(&subseq($targetname, $tfrom, $tto-$tfrom, $complement)."--");
	    my $frameshift = ($tto-$tfrom) % 3;
	    if ($frameshift && $SHOW_FRAME_IN_TRANS) {
		&run_bug("Frameshift not translated to 'X'!\n") unless s/X$/$frameshift/;
	    }
	    if ($SPLIT_ON_FS) {
		$current->{translation} = $translation;
	    } else {
		my $translen = length($translation);
		if ($qnto == $qnfrom) {
		    $insertions += $translen;
		} else {
		    $gapcount -= $translen;
		    $undetermined += $translen;
		    push @{$previous->{undeterminedlist}}, ($qfrom+1 .. $qfrom + $translen);
		    push @{$previous->{gaplist}}, ($qfrom + $translen +1 .. $qto);
		}
		$previous->{translation} .= $translation;
		push @{$previous->{seqshifts}}, $current;
		while ($translation =~ /[*]/g) {
		    push @{$previous->{inframe_stopcodons}}, $qfrom;
		}
		next;
	    }
	}
	if (/exon/) 
	{
	    my $diffstr;
	    (shift @diff_strings) =~ /^(.*):(.*)$/;
	    ($current->{translation}, $diffstr) = ($1, $2);
	    $current->{$_} = [] 
		foreach("seqshifts","mismatchlist","gaplist",
			"undeterminedlist","inframe_stopcodons");
	    while ($diffstr =~ /([X ])/g) {
		my $qpos =  $-[0] + $qfrom + 1;
		if ($1 eq 'X') {
		    push @{$current->{mismatchlist}}, $qpos;
		    ### DEBUG PART ###
		    push @mismatchpos, $qpos;
		    ### END OF DEBUG PART ###
		    $mismatchcount++;
		} elsif ($1 eq ' ') {
		    push @{$current->{undeterminedlist}}, $qpos;
		    $undetermined++;
		}
	    }
	    while ($current->{translation} =~ /[*]/g) {
		push @{$current->{inframe_stopcodons}}, $-[0] + $qfrom + 1;
	    }
	    if ($previous && $previous->{type} eq "exon")
	    {   # add to previous exon...
		@$previous{"${TPFX}_end","${QNPFX}_end","${QPFX}_end"} = ($tto, $qnto, $qto);
		push (@{$previous->{$_}}, @{$current->{$_}} )
		    foreach ("mismatchlist","undeterminedlist","inframe_stopcodons");
		$previous->{translation} .= $current->{translation};
		# and don't put on list
		next;
	    }
	} elsif (/intron/) {
	    @$current{"${QNPFX}_pos", "${QPFX}_pos"} = ($qnfrom, $qfrom) ;
	}
	$current->{type} = $_;
	push @matchings, $current;
    } # end foreach (@types)

    ### debug part
    my @mismatchcounts = map { scalar @{$_->{mismatchlist}} } grep { $_->{type} eq "exon" }  @matchings;
    if ($mismatchcount != sum (@mismatchcounts,0))
    {
	print STDERR "ID: $.\nMismatchcount: $mismatchcount\n".
	    scalar(@mismatchpos),"\naufgeschlüsselt: (".join("+",@mismatchcounts).")\n";
	print STDERR "Exons: ".join("/", @targetlocations)."\n";
	&run_bug("Incorrect calculation of mismatches!");
    }
    &run_bug("Schlimm!") unless defined ($qTos[-1]);
    ### end of debug part

    # matchcount is the number of matched positions of the query sequence
    my $matchcount = $matchsize - $mismatchcount - $gapcount - $undetermined; 

    # the score is the fraction of matches minus mismatches; 
    # any query should have at least one hit with a score of $MIN_BEST_SCORE (=0.3)
    my $score = ($matchcount-$mismatchcount)/$querysize;
    
    ### Collect all data for this BLAT hit
    my $hitref = { "${QPFX}_len" => $querysize, "${QPFX}_start" =>  $querystart, "${QPFX}_end" => $queryend,  
		   "target" => $targetname, "target_len" => $dbase{$targetname}{length},
		   "${TPFX}_start" => $targetstart, "${TPFX}_end" => $targetend, 
		   "strand" => $VALUE_FORMAT{strand}{$complement?'-':'+'}, "complement" => $complement, 
		   "stopcodon" => $stopcodon,
		   "score" => $score, "unmatched" => $gapcount,
		   "undetermined" => $undetermined, "additional" => $insertions,
		   "matches" => $matchcount, "mismatches" => $mismatchcount, "mismatchcounts" => \@mismatchcounts,
		   "ID" => $., "matchings" => \@matchings };
    next if (($mismatchcount >= $MAX_MISMATCH && $MAX_MISMATCH) || 
	     $matchcount/($mismatchcount+$matchcount) < $MIN_IDENTITY ||
	     ($mismatchcount+$matchcount)/$matchsize < $MIN_COVERAGE);
    push @{$hits{$queryname}}, $hitref;
    $dbase{$targetname}{users}{$.} = 1;
} # end BLATLINES
print COMMENT "done\n";
print COMMENT sum (map {scalar @{$hits{$_}}} keys %hits)." of $hitcount hits saved.\n";
unlink $BLAT_OUTPUT unless ($KEEP_BLATOUT);

# free memory of unused sequences
my $delcount = 0;
foreach (keys %dbase) {
    unless (keys %{$dbase{$_}{users}}) { 
#	print COMMENT "Deleting sequence $_ - not needed anymore.\n";
	$delcount++;
	delete $dbase{$_};
    } 
} 


print COMMENT "Discarded $delcount unused sequences.\n" if ($delcount);
print COMMENT "Now assembling BLAT hits...\n";
#print COMMENT "[";
$delcount = 0;

### Sort BLAT hits
# before postprocessing is completed, all results are stored in %hits which is a hash 
# indexed by query names whose values are references to lists, each list entry
# referring to the information about one hit for the query

# after postprocessing, only one hitref (or sequence of consecutive hits) will remain
# for each query, and sequences are added

my @missing = ();
my @hitkeyqueries = keys %hits;

foreach my $queryname (@hitkeyqueries) {
    my @usercounts = map { scalar keys %{$_->{users}} } values %dbase;
#    print COMMENT "$queryname: (".join(",",@usercounts).") => ".sum(@usercounts);
    my $hitlist = $hits{$queryname};

    my $is_copy = $queryname =~ s/_\((\d+)\)$//;
    my $altnames = $is_copy ? $1 : 0;

    if (@$hitlist) {
	@$hitlist = sort { $b->{"score"} <=> $a->{"score"} } @$hitlist;
	if (grep { ! defined $_->{"${QPFX}_len"} } @$hitlist) {
	    print STDERR "${QPFX}_len undefined!\n";
	}
#	unless (grep { ($_->{"${QPFX}_end"} - $_->{"${QPFX}_start"}) / $_->{"${QPFX}_len"} >= $MIN_BEST_SCORE } @$hitlist) {
	if ($hitlist->[0]{score} < $MIN_BEST_SCORE) {
	    printf COMMENT
 		"Hits for '$queryname' discarded since best piece has low score. Set --min_score=x, x < %5.3f to make it visible.\n", 
		$hitlist->[0]{score};
 	    @$hitlist = ();
	}
    }


    unless (@$hitlist)
    {
	push @missing, $queryname; 
	delete $hits{$queryname};
	next;
    }

    my $incomp_hits = [];
    for (my $i = 1; $i<@$hitlist; $i++) {
	foreach (@$hitlist[0..($i-1)])	{
	    if ($ONE_TARGET_ONLY || &incompatible($_, $hitlist->[$i]) || &high_margins($hitlist->[$i], $_) ) {
		my $removed_hit = splice(@$hitlist, $i, 1);
		if ($removed_hit->{score}>=1 || 
		    ($MULTIPLE_RESULTS && (@$incomp_hits || $removed_hit->{score} >= $MIN_BEST_SCORE))) {
		    push @$incomp_hits, $removed_hit;
		} else {
		    &remove_user(@{$removed_hit}{"target","ID"});
		}
		$i--; last;  # do not increment $i since list was spliced
	    }
	}
    } #
    if (@$incomp_hits) {
	my $altqueryname = "${queryname}_(".($altnames + 1).")";
	push @hitkeyqueries, $altqueryname;
	$hits{$altqueryname} = $incomp_hits;
    } elsif ($altnames) {
	print COMMENT 
	    "Added $altnames additional high scoring hit".
	    ($altnames==1 ? " as query" : "s as queries")." ${queryname}_(1)".
	    ($altnames==1 ? "" : $altnames==2 ? ",_(2)" : ",...,_($altnames)").".\n";
    }

    @$hitlist = sort { $a->{"${QPFX}_start"} <=> $b->{"${QPFX}_start"} } @$hitlist;
    for (my $i = 0; $i < @$hitlist-1; ) {
	my ($qfrom, $qto, $nextqfrom, $nextqto) = map { @$_{"{QPFX}_start","{QPFX}_end"} } @$hitlist[$i, $i+1];
	if ($qfrom = $nextqfrom) {
	    splice(@$hitlist, $i + ($qto < $nextqto ? 0:1), 1);
	} else {
	    $i++;
	}
    }

    my $queryseq = $queries{$queryname};
    

    # adjusting with neighbouring hits
    foreach (0..$#$hitlist) {
	my $current = $hitlist->[$_];
	my ($targetname, $complement, $tfrom, $tto, $tlen, $qto, $qlen) 
	    = @{$current}{"target", "complement", "${TPFX}_start", "${TPFX}_end", 
			  "target_len", "${QPFX}_end", "${QPFX}_len"};
	if ($_==0) {
	    my $upstream_size = min ($REGION_SIZE, $tfrom + ($complement ? $tlen : 0));
	    $current->{upstream} = &subseq($targetname, $tfrom - $upstream_size, $upstream_size, $complement);
	    my $qfrom = $current->{"${QPFX}_start"};
	    push @missing, "${queryname}[1..$qfrom]" if ($qfrom > 0 && !$is_copy);
	} 
	# no. of nucleotides on previous/next target sequence not covered by the hits
	my ($prevnucs, $postnucs); 
	$prevnucs = ($complement ? 0 : $tlen) - $tto;   

	if ($_< $#$hitlist) {
	    my $nexthit = $hitlist->[$_+1];
	    my ($nexttfrom, $nextqfrom, $nexttlen) = @{$nexthit}{"${TPFX}_start", "${QPFX}_start", "target_len"};
	    while (defined($current->{matchings}[-1]{"${QPFX}_start"}) && defined($nextqfrom) && $current->{matchings}[-1]{"${QPFX}_start"} >= $nextqfrom) {
		last if scalar(@{$current->{matchings}}) <= 2;
		splice(@{$current->{matchings}},-2,2);
		my @ch_keys = ("${TPFX}_end", "${QPFX}_end");
		($tto, $qto) = @$current{@ch_keys} = @{$current->{matchings}[-1]}{@ch_keys};
		$prevnucs = ($complement ? 0 : $tlen) - $tto;   
	    }
	    while (defined($nexthit->{matchings}[-1]{"${QPFX}_end"}) && defined($qto) && $nexthit->{matchings}[-1]{"${QPFX}_end"} <= $qto) {
		last if scalar(@{$nexthit->{matchings}}) <= 2;
		splice(@{$nexthit->{matchings}},0,2);
		my @ch_keys =  ("${TPFX}_start", "${QPFX}_start");
		($nexttfrom, $nextqfrom) = @$nexthit{@ch_keys} = @{$nexthit->{matchings}[0]}{@ch_keys};
	    }

	    $postnucs = ($nexthit->{complement}? $nexttlen : 0) + $nexttfrom; 

	    my $t_o = $current->{target_overlaps};                   # copy to variable to prevent autovivification of key
	    my $has_overlap = grep { $_ == $nexthit->{ID} } @$t_o;   # true if ID of next hit is among overlaps


	    my $qgaplen = $nextqfrom-$qto;                           # no. of query positions not covered by hits
	    my $additional = $prevnucs+$postnucs-$qgaplen*3;         # no. of nucleotides left after filling the qgap
#	    if ($qgaplen <= $MAX_GAPLEN || $has_overlap) {
		if ($has_overlap && $qgaplen < 0 && $prevnucs<3 && $postnucs<3) {   
                    # case 1: nucleotides have been counted twice, no intron here, 
		    #         but additional (= length of overlap) is >0
		    #         (relevant only when there are query overlaps and both
                    #          targets end inside exon)

		    #  eeeeeee..]
                    #      [..eeeeeee
		    #   ||||||||||||
                    #   AAAAAAAAAAAA             

		    # uncovered nucleotides
		    my $postrem =
			&subseq($nexthit->{target}, $nexttfrom-$postnucs, $postnucs, $nexthit->{complement});

		    # this is to be passed to &replace_in_next later
		    my $dna_add = "n" x (-$postnucs%3);
		    $current->{matchings}[-1]{overlap} = $additional;
		    $additional-=$prevnucs;
		    

		    # remove the overlap from the preceding exon
		    &replace_in_last($current, "", "", $postrem, $additional);
		    $qto-=&nu_to_aa($additional);
		    $prevnucs = 0;
		    if ($postnucs % 3 == 2) { 
			$dna_add = &subseq($targetname, $tto-$additional-1,1,$complement);
		    }
		    my $postaas = &nu_to_aa($postnucs);
		    &replace_in_next($nexthit, $postrem, substr($queryseq, $nextqfrom-$postaas, $postaas),
				     $dna_add);
		    $nextqfrom -= $postaas;
		    $postnucs = 0;
		} elsif ($additional <= 0) { 
                    # case 2: $additional <= 0: we don't have introns here, rather too few nucleotides to fill the gap    
                    #       
		    #  eee..]-----[..eee      additional is minus length of the gap --; prevnucs
                    #  |||           |||      and postnucs are the lengths of the .. regions
		    #  AAAAAAAAAAAAAAAAA 

		    # uncovered nucleotides
		    my $postrem = 
			&subseq($nexthit->{target}, $nexttfrom-$postnucs, $postnucs, $nexthit->{complement});

		    # this is to be passed to &replace_in_next later
		    my $dna_add = "n" x (-$postnucs%3);

		    # close small gaps despite possible mismatches
		    # in preceding contig
		    if ($prevnucs <= 4) {
			my $prevrem = &subseq($targetname, $tto, $prevnucs, $complement);
			&replace_in_last($current, $prevrem, substr($queryseq, $qto, &nu_to_aa($prevnucs)),
					 ($postnucs <= 4 && $additional==0) ? $postrem : "nn");
			$qto += &nu_to_aa($prevnucs);
			$prevnucs = 0;
			if ($additional == 0) { $dna_add = $prevrem; }
		    } 
		    # in following contig
		    if ($postnucs <= 4) {
			my $postaas = &nu_to_aa($postnucs);
			&replace_in_next($nexthit, $postrem, substr($queryseq, $nextqfrom-$postaas, $postaas),
					 $dna_add);
			$nextqfrom -= $postaas;
			$postnucs = 0;
		    }
		} elsif ($qgaplen <= $MAX_GAPLEN && $prevnucs < $ASSEMBLE_SIZE && $postnucs < $ASSEMBLE_SIZE) { 
		    # case 3: $additional > 0, qgaplen small (so no exon is missing in between, and prevnucs/postnucs not too big
		    # 
		    # eee.....][.....eee
		    # |||            |||
		    # AAA--(intron)--AAA

		    my $left_border = $current->{matchings}[-1]{"${QPFX}_start"} + $MIN_REMAIN;
		    my $right_border = $nexthit->{matchings}[0]{"${QPFX}_end"} - $MIN_REMAIN;
		    my $left_overlap = min ( max ( $PRED_OVERLAP, $PRED_OVERLAP-$qgaplen ), max ($qto - $left_border, 0) );
		    my $right_overlap = min ( max ( $PRED_OVERLAP, $PRED_OVERLAP-$qgaplen ), max ($right_border - $nextqfrom, 0));
		    my $qgaplenplus = $qgaplen + $left_overlap + $right_overlap;
		    if ($left_border > $qto || $right_border < $nextqfrom || $qgaplenplus < 0) {
			&run_bug("This is a bug and shouldn't happen. Please report.");
		    }
		    my $asslen = min ( $postnucs+3*$right_overlap, 3*$qgaplenplus +2);
		    my ($dss_seq, $ass_seq, $gap_aa) =
			(&subseq($targetname, $tto-3*$left_overlap, 
				 min ($prevnucs+3*$left_overlap, 3*$qgaplenplus +2),
				 $complement),
			 &subseq($nexthit->{target}, $nexttfrom+3*$right_overlap-$asslen, $asslen, $nexthit->{complement}),
			 substr($queryseq, $qto-$left_overlap, $qgaplenplus));
		    my $preferred_offset = 3*int($qgaplenplus/2);
		    my ($get_offset, $mismatchadd, $gap_trans, $splicesite_score) =
			#always look for splice sites
			&get_intron_offset($gap_aa, $dss_seq."xx", "xx".$ass_seq, 0 , $preferred_offset); # todo: $additional < $MIN_INTRON_LEN
		    if (($mismatchadd <= 2 && $splicesite_score > 0.5) || $qgaplen<=0 ) { # only add sure introns in the end
			my $asslen = 3*$qgaplenplus - $get_offset;
			my ($dss, $ass) = (substr($dss_seq."xx",0,$get_offset+2),substr("xx".$ass_seq,-$asslen-2));
			my $dss_inpart = substr($dss,-2,2,"");
			my $ass_inpart = substr($ass,0,2,"");
			
			# for replace_in_last we just set $dna_add to length 2
			&replace_in_last($current, $dss,
					 substr($gap_aa, 0, &nu_to_aa($get_offset)),
					 substr($ass,0,2), 3*$left_overlap );
			$prevnucs += ($get_offset - 3*$left_overlap);
			# for replace_in_next we set $dna_add such that its length equals the reading frame
			&replace_in_next($nexthit, $ass,
					 substr($gap_aa, &nu_to_aa($get_offset)),
					 substr($dss,-($get_offset%3), $get_offset%3),
 					 3*$right_overlap);
			$postnucs += (3*$qgaplen+ 3*$left_overlap-$get_offset);
			$qto += &nu_to_aa($get_offset) - $left_overlap;
			$nextqfrom -= ($qgaplenplus - $right_overlap - &nu_to_aa($get_offset));
			## DEBUG PART
			if ($qto != $nextqfrom) {
			    &run_bug("Internal error. Serious bug with intron offset calculation. Please report!");
			}
			## END OF DEBUG PART
			my ($dsstype, $asstype) = ("intron", "intron");
			if ($splicesite_score <= 0.5) {
			    $dsstype .= "?" unless ($dss_inpart =~ /^g[ct]$/);
			    $asstype .= "?" unless ($ass_inpart eq "ag");
			}
			&add_final($current,$dsstype) if $prevnucs > 0; $prevnucs=0;
			&add_initial($nexthit,$asstype) if $postnucs > 0; $postnucs=0;
		    } # mismatchadd <= 2
		    ## DEBUG PART
		    elsif ($qto == $nextqfrom && $mismatchadd > 2) {
			print STDERR "current: ".$current->{ID}." next hit: ".$nexthit->{ID}."\n";
			print STDERR "prevnucs=$prevnucs postnucs=$postnucs mismatchadd=$mismatchadd\n";
 			substr($dss_seq,$left_overlap*3,0)=".";
 			substr($ass_seq,-$right_overlap*3,0) = "." if ($right_overlap);
			print STDERR "to=nextfrom=$qto additional=$additional length=".length($queryseq)."\n";
			print STDERR "dss=$dss_seq ass=$ass_seq match=".substr($queryseq, $qto-$left_overlap, $left_overlap)."()".substr($queryseq,$nextqfrom,$right_overlap)."\n";
			print STDERR "offset=$get_offset\n";
			&run_bug("Internal error. Adding no aas should not add errors. Please report!");
		    }
		    ## END OF DEBUG PART
		}  # additional > 0
#	    } # gaplen <= MAX_GAPLEN
	    if ($qto < $nextqfrom) {
		# add gap sequences in all cases in which we have a gap; 
		# empty gap sequences are needed later to determine that there is a gap to neighbouring hits
		$current->{downstream_gap} = &subseq($targetname, $tto, $prevnucs, $complement);
#		my $next_upstream_size = min ($REGION_SIZE, $postnucs);
		$nexthit->{upstream_gap} = &subseq($nexthit->{target}, $nexttfrom-$postnucs, $postnucs, $nexthit->{complement});
		push @missing, "${queryname}[".($qto+1)."..$nextqfrom]";
	    } elsif ($prevnucs > 0 || $postnucs > 0)  {
		### DEBUG PART
		&run_bug ("Internal error: gap is closed - we should not have regions to show here. Please report!");
		### END OF DEBUG PART
	    }
	} else { # $_==$#$hitlist  or no assembling
	    $current->{downstream} = &subseq($targetname, $tto, min($REGION_SIZE, $prevnucs), $complement);
	    push @missing, "${queryname}[".($qto+1)."..$qlen]" if ($qto < $qlen && !$is_copy)
	}  # end if ($_ < $#$hitlist

	foreach (@{$current->{matchings}}) {
	    my ($reftfrom, $reftto, $refqnfrom, $refqnto) = \@$_{"${TPFX}_start","${TPFX}_end", "${QNPFX}_start", "${QNPFX}_end"};
	    my ($qfrom, $qto) = @$_{"${QPFX}_start","${QPFX}_end"} = map ( &nu_to_aa($_), ($$refqnfrom, $$refqnto) );
	    $_->{seq} = &subseq($targetname, $$reftfrom, $$reftto-$$reftfrom, $complement);
	    if ($_->{type} eq "intron") {
		if ($_->{seq} =~ s/\.{3}(\.+)/.../) { $$reftto -= length($1); }
		elsif ($_->{seq} =~ s/^(\.+)\.{3}/.../) { $$reftfrom += length($1); }
	    } elsif ($_->{type} eq "exon" || $_->{type} eq "gap" ) {
		$_->{"${QPFX}_seq"} = substr($queryseq, $qfrom, $qto-$qfrom);
	    }
	}
	&remove_user($targetname,$current->{ID});
	my $qfrom = $current->{ "${QPFX}_start" };
	## DEBUG PART
	if ($qto != $current->{"${QPFX}_end"}) {
	    print STDERR "";
	    &run_bug("BUG: coordinates mismatch. id:".
		$current->{ID}." qto=$qto prot_end=".$current->{"${QPFX}_end"}.
		". Please report!");
	}
	## END OF DEBUG PART
	$current->{"${QPFX}_seq"} = substr($queryseq, $qfrom, $qto-$qfrom);
	$current->{score} = sprintf("%5.3f",($current->{matches} - $current->{mismatches})/$qlen);
	&set_problem_field($current);

	# choose keys for YAML output
	my @bless_list = split ",",$HIT_BLESS_LIST;
	my %default_values = ( 
	    matches => $current->{"${QPFX}_len"}, 
	    mismatches => 0,
	    unmatched => 0,
	    undetermined => 0,
	    additional => 0,
	    upstream_gap => "",
	    downstream_gap => "",
	    "${QPFX}_start" => 0,
	    "${QPFX}_end" => $current->{"${QPFX}_len"} );
	if ($HIDE_UNDEF) {
	    @bless_list = grep { defined $current->{$_} } @bless_list;
	}
	if ($HIDE_DEFAULTS) {
	    foreach (@bless_list) {
		print STDERR "\$current->$_\n" unless defined $current->{$_};
		print STDERR "\$default_values{$_}\n" if (exists $default_values{$_} && !defined $default_values{$_});
	    }
	    @bless_list = grep { ! (exists $default_values{$_} && $current->{$_} eq $default_values{$_}) } @bless_list;
	}
	YAML::Bless($current)->keys(\@bless_list);
	foreach my $match (@{$current->{matchings}}) {
	    my $type = $match->{type};
	    $type = "intron" if ($type =~ /intron/);
	    $type = "any" unless (exists $MATCH_BLESS_LISTS{$type});
	    my @bless_list = split ",",$MATCH_BLESS_LISTS{$type};
	    if ($type eq "exon") {
		foreach my $fs (@{$match->{seqshifts}}) {
		    YAML::Bless($fs)->keys( [ grep { exists $fs->{$_} } @bless_list ] );
		}
	    }
	    if ($HIDE_UNDEF) {
		foreach (values %$match) {
		    $_ = undef if ( ref($_) eq 'ARRAY' && !@$_);
		}
		@bless_list = grep { defined $match->{$_} } @bless_list;
	    }
	    YAML::Bless($match)->keys(\@bless_list);
	}
    } # end foreach (@hitlist)
    print STDERR "";
} # end foreach (@hitkeyqueries)


###################### output of results ################################

print "### Scipio v1.4 output\n"; 
print "# query file    $queryfile\n";
print "# target file   $dbasefile\n";
print "# BLAT output   $BLAT_OUTPUT\n" if $KEEP_BLATOUT;
my @time = localtime(time); $time[5]+=1900;
printf "# timestamp     %4d-%02d-%02d %02d:%02d:%02d\n", reverse @time[0..5];
print "# parameters    ".($commandline ? "$commandline\n" : "(default)\n");

my @querybless = grep { defined $hits{$_} } sort keys %hits;
YAML::Bless(\%hits)->keys(\@querybless);
print YAML::Dump(\%hits);


if (@missing) {
     print "### not found: ".(join ",",@missing)."\n\n";
}
print "### end of Scipio output\n";
exit scalar @missing;


####################### function definitions ####################################

sub run_bug {
    my $errmess = shift;
    if ($DONT_DIE) { print STDERR "$errmess\n" }
    else { die $errmess }
}

sub revcom {
    my $result = reverse shift;
    $result =~ tr/tcyawmhgksbrdvnTCYAWMHGKSBRDVN/agrtwkdcmsvyhbnAGRTWKDCMSVYHBN/;
    return $result;
}

sub with_number {
    my ($n, $s, $alt) = @_;

    $alt = "${s}s" unless (defined $alt);
    return ($n==1) ? "$n$s" : "$n$alt";
}

sub diff_str {
    my ($s1, $s2, $count, $s1from, $s2from) = @_;
    $s1from = 0 unless (defined $s1from);
    $count = length($s1)-$s1from unless (defined $count);

    $s2from = $s1from unless (defined $s2from);
    if (length($s1)<$s1from+$count || length($s2)<$s2from+$count)  {
	print STDERR "diff_str returns undef: input is\n" .
	    "s1: $s1\ns2: $s2\ncount: $count\nfrom1 :$s1from\nfrom2: $s2from\n";
    }
    return undef if (length($s1)<$s1from+$count || length($s2)<$s2from+$count);
    my @s = split //, substr($s1, $s1from, $count);
    my @t = split //, substr($s2, $s2from, $count);
    
    my $result=""; 
    foreach (@s)  {  # $_: translation <-> $comp: query
	my $comp = shift(@t);
	### DEBUG part
	&run_bug ("'-' in query or translation!!!") if ($comp eq '-' || $_ eq '-');
	### end of DEBUG part
	# count 'X' as undetermined
	# others as match/mismatch
	$result .= ($comp =~ /^[X12-]$/ || /^[X12-]$/) ?  " " : ( $comp eq $_ ? "|" : "X" );
	
    }
    return $result;
}

sub str_diff_list {
    my @result=(); my $pattern=shift;
    local $_ = &diff_str(@_);
    push @result, $-[0] while(/$pattern/g);
    return @result;
}

sub str_diff_count_all {
    local $_ = &diff_str(@_);
    s/\|//g;
    return length;
}

sub str_diff_count_undeterm {
    local $_ = &diff_str(@_);
    s/[^ ]//g;
    return length;
}

sub str_diff_count_mismatches {
    local $_ = &diff_str(@_);
    s/[^X]//g;
    return length;
}

sub str_pos { 
    my @result=(); my $pattern=shift;
    local $_ = shift;
    push @result, $-[0] while(/$pattern/g);
    return @result;
}

sub extend_to_full_name {
    my $result = $full_name{$_[0]};
    $_[0] = $result if (defined $result);
}

sub save_target_region {
    my ($targetname, $from, $to, $total_length) = @_;
    $from = 0 if ($from<0);
    return if $to<=$from;

    unless (exists $dbase{$targetname})  {
	$dbase{$targetname} = { "length" => $total_length, "users" => {} };
    }
    my $target = $dbase{$targetname};
    $target->{length} = $total_length if ($target->{length} < $total_length);

    foreach (keys %$target)  {
	my ($piecefrom,$pieceto) = /^(\d+),(\d+)$/;
	next unless (defined $pieceto);
	return if ($piecefrom <= $from && $to <= $pieceto);
	if ($from <= $piecefrom && $pieceto <= $to) {
	    delete $target->{$_};
	    next;
	}
	if ($piecefrom < $from && $from <= $pieceto) {
	    delete $target->{$_};
	    $from = $piecefrom;
	}
	elsif ($piecefrom <= $to && $to < $pieceto) {
	    delete $target->{$_};
	    $to = $pieceto;
	}
    }
    $target->{"$from,$to"} = "";
}

sub cut_target_regions {
    foreach my $target (values %dbase) {
	my $length = $target->{length};
	foreach (keys %$target) {
	    my ($piecefrom,$pieceto) = /^(\d+),(\d+)$/;
	    next unless (defined $pieceto && $pieceto>$length);
#	    my $piece = $target->{$_};
	    delete $target->{$_};
	    $target->{"$piecefrom,$length"}="" if ($piecefrom < $length);
	}
    }
}
	
sub fill_target_region {
    my ($targetname, $sref, $offset) = @_;
    my $target = $dbase{$targetname};
    my $end = $offset+length($$sref);
    
    foreach (keys %$target) {
	my ($piecefrom,$pieceto) = /^(\d+),(\d+)$/;
	next unless (defined $pieceto);
	next if ($offset >= $pieceto || $piecefrom >= $end);
	if ($target->{$_}) { # we already have this target; something is wrong
	    die "Error: multiple target sequences with the same name '$targetname' found! Aborting...\n"
	}
	if ($piecefrom < $offset) {
	    delete $target->{$_};
	    $target->{"$piecefrom,$offset"} = "";
	    $piecefrom=$offset;
	}
	if ($pieceto > $end) {
	    delete $target->{$_};
	    $target->{"$end,$pieceto"} = "";
	    $pieceto=$end;
	}
	$target->{"$piecefrom,$pieceto"}=substr($$sref,$piecefrom-$offset,$pieceto-$piecefrom);
    }
}	

sub read_targets {
    my ($targetname, $offset);
    my $seq = "";
    open FASTAFILE, shift;
    my $cutoff_warning_no=0;
    my $seqcount=0;

    while (<FASTAFILE>) {
	if (s/^>\s*//) { # new target sequence
	    if ($seq) {  # save old sequence
		&fill_target_region($targetname, \$seq, $offset);
		$seq = "";
	    }
	    chomp;
	    if (/^(\S*)(\s+.*)$/) { # new target name contains white spaces
		if (exists $dbase{$1}) {
		    # the keys for %dbase come from the .psl file in which they were cut off
		    # now reconstruct them if they were present
		    $cutoff_warning_no++;
		    $full_name{$1}=$_;
		    $dbase{$_}=$dbase{$1};
		    delete $dbase{$1};
		}
	    }
	    $offset = 0;
	    $targetname = $_;
	    if (exists $dbase{$targetname} && $seqcount <= 80) {
		print COMMENT ($seqcount == 80) ? "[...] " : ".";
		$seqcount++;
	    }
	} elsif (exists $dbase{$targetname}) { # ignore sequences not among BLAT results
	    chomp;
	    s/\s//g;
	    $seq.=lc;
	    ### don't let $seq become larger than 20MBases
  	    if (length($seq) > 20000000) {
  		&fill_target_region($targetname, \$seq, $offset);
  		$offset += length($seq);
  		$seq = "";
 	    }
	}
    }
    if ($seq) {
	&fill_target_region($targetname, \$seq, $offset);
    }
    print COMMENT "done\n";	
    #if ($cutoff_warning_no) {
	#print STDERR "Warning: Some target sequence names stripped at white spaces.\n";
    #}
}

sub validate_targets {
    while (my ($targetname, $target) = each %dbase) {
	foreach (keys %$target) {
	    if (/^(\d+),(\d+)$/ && ($2 > $1) && $target->{$_} eq "")
	    {
		die "Missing or incomplete target sequence '$targetname'. (You can use --partial_target\n".
		    " to enforce a scipio run with missing targets ignored.)\n".
		    "Nonexistent targets in psl file" unless ($ALLOW_MISSING_TARGETS);
		delete $dbase{$targetname};
		last;
	    }
	}
    }
}

sub nu_to_aa {
    # if changing the assignment of nucleotides, check in particular: postnucs(l. 863)
    my $arg = shift;
    my $residue = ($arg+1) % 3 -1;   # -1, 0, 1
    return ($arg-$residue)/3;
}

sub revtranslate_to_pattern {
    my $seq = uc shift;
    return $seq =~ /[^ABCDEFGHIKLMNPQRSTVWYZ*]/? undef :  
	join ("", map { "(".join("|",$ttable->revtranslate($_)).")" } split(//, $seq));
}

sub split_codons {
    my ($aa, $dna, $backward) = @_;
    my @result=();
    my @codons = $ttable->revtranslate($aa);
    foreach (@codons) {
	/^(.)(.)(.)$/;
	my ($c1,$c2,$c3) = ($1,$2,$3);
	if ($backward) {
	    if ($dna =~ /ag$c3$/) { 
		push @result, "$c1${c2}g[ct]"; 
	    } else {
		if ($dna =~ /ag$c2$c3$/) { push @result, "${c1}g[ct]"; }
		if ($dna =~ /ag$/) { push @result, "${_}g[ct]"; }
	    }
	} else {
	    if ($dna =~ /^$1g[ct]/) {
		push @result, "ag$c2$c3"; 
	    } else {
		if ($dna =~ /^${c1}${c2}g[ct]/) { push @result, "ag$c3"; }
		if ($dna =~ /^g[ct]/) { push @result, "ag$_"; }
	    }
	}
    }
    return @result;
}

sub get_splicesite_patterns_forward {
    my ($suffix, $dna) = @_;
    $suffix =~ s/^(.)(.)//;
    my $longpatt = &revtranslate_to_pattern($suffix);
    return undef unless defined $longpatt;
    my $pattern2 = &revtranslate_to_pattern($2);
    return undef unless (defined $pattern2);
     my @result = &split_codons($1, $dna, 0);
    @result = ( "(".join("|", @result).")$pattern2") if @result;
    push @result, &split_codons($2, substr($dna,3),0);
    return undef unless (@result);
    return "(".join("|", @result).")$longpatt";
}
   
sub get_splicesite_patterns_backward {
    my ($prefix, $dna) = @_;
    $prefix =~ s/(.)(.)$//;
    my $longpatt = &revtranslate_to_pattern($prefix);
    return undef unless (defined $longpatt);
    my $pattern1 = &revtranslate_to_pattern($1);
    return undef unless (defined $pattern1);
    my @result = &split_codons($2, $dna, 1);
    @result = ( "$pattern1(".join("|",@result).")" ) if @result;
    substr($dna,-3)="";
    push @result, &split_codons($1, $dna, 1);
    return undef unless (@result);
    return "$longpatt(".join("|", @result).")";
}

sub subseq {
    my ($targetname, $from, $count, $complement) = @_;
 
    return undef if ($count<0);
    return "" if ($count==0);
    
    my $target = $dbase{$targetname};
    return undef unless (defined $target);
    $from=-$from-$count if ($complement);
    
    my $prelap = max ( 0, -$from);
    my $postlap = max ( 0, $from+$count-$target->{length} );
    if ($count < $prelap+$postlap) {
	return " " x $count;
    }
    my $result = (" " x $prelap).("." x ($count-$prelap-$postlap)).(" " x $postlap);
    foreach (keys %$target) {
	next unless (/^(\d+),(\d+)$/);
	my ($piecefrom,$pieceto) = ($1, $2);
	next unless ($from<$pieceto && $from+$count>$piecefrom);
	    
	my ($pieceoffset,$offset,$piececount) = ($from-$piecefrom,0,$count);
	if ($from<$piecefrom) {
	    $offset = $piecefrom-$from;
	    $pieceoffset = 0;
	    $piececount -= $offset;
	}
	if ($from+$count>$pieceto) {
	    $piececount = $pieceto-$piecefrom-$pieceoffset;
	}
 	substr($result, $offset, $piececount) = substr($$target{$_},$pieceoffset,$piececount);
    }
    return $complement ? 
	&revcom($result) : $result;
}

sub remove_user {
    my ($targetname, $ID) = @_;
#    my $oldkeycount = keys %{$dbase{$targetname}{users}};
    delete $dbase{$targetname}{users}{$ID};
    my $newkeycount = keys %{$dbase{$targetname}{users}};
    unless ($newkeycount) {
	delete $dbase{$targetname};
#	print COMMENT "Sequence '$targetname' discarded - not needed anymore.\n";
    } # elsif ($oldkeycount == $newkeycount) {
# 	print COMMENT "No user removed for '$targetname'\n";
#     } else {
# 	print COMMENT "$newkeycount users remain for '$targetname'\n";
#     }
}
    
sub get_intron_offset {   
    my ($gap_aa, $dss_seq, $ass_seq, $frameshift_mode, $preferred_offset) = @_;
    my $gaplen=length($gap_aa)*3;
    my ($bestoffset,$bestmismatchcount,$gap_trans,$splicesite_score) = (0,length($gap_aa),"-" x $gaplen,0); 
    
    my $ssb = $frameshift_mode ? 0 : 2; # in intron mode, account for splice sites
    
    my ($leftend, $rightend) = ($gaplen+$ssb-length($ass_seq), length($dss_seq)-$ssb);
    $leftend=0 if ($leftend < 0);
    $rightend=$gaplen if ($rightend > $gaplen);
    foreach my $offset  (reverse $leftend..$rightend) {
	next if ($frameshift_mode && $offset % 3);
	my $dna_comp = substr($dss_seq, 0, $offset);
	my $iscore_cand;
	$dna_comp .= substr($ass_seq, $offset-$gaplen) if ($offset < $gaplen);
	$dna_comp =  $ttable->translate($dna_comp);


	unless ($frameshift_mode) {
	    $iscore_cand = 2.0 - &get_splice_penalty(substr($dss_seq, $offset, 2),substr($ass_seq, $offset-$gaplen-2, 2));
	    $iscore_cand -= 0.01 unless $offset % 3 == 0;
	}

	my $mismatchcount = &str_diff_count_mismatches($dna_comp, $gap_aa); #no of mismatches between dna_comp and gap
	$mismatchcount += 0.5 * &str_diff_count_undeterm($dna_comp, $gap_aa); # punishment for undetermined
	$mismatchcount += ($dna_comp =~ /[*]/) * $IFS_MALUS;                  # punishment for stop codons
	my $optimiser = $mismatchcount - $bestmismatchcount;
	$optimiser -= ($iscore_cand - $splicesite_score) unless ($frameshift_mode);
	if ($optimiser < 0 || ($optimiser == 0 && $offset >= 2*$preferred_offset - $bestoffset)) {
	    ($bestoffset,$bestmismatchcount,$gap_trans,$splicesite_score)
		=($offset,$mismatchcount,$dna_comp,$iscore_cand);
	    last if ($mismatchcount==0 && ($frameshift_mode || $iscore_cand >= 2.0) && $offset <= $preferred_offset);
	}
    }
   
    return ($bestoffset, $bestmismatchcount, $gap_trans, $splicesite_score);
}  ### end sub &get_intron_offset

sub calculate_diff_strings {
    my ($targetname, $queryname, $targetlocations, $querylocations, $complement) = @_;
    my $transcript = ""; 
    for (my $i=0; $i<@$targetlocations; $i+=2) {
	my ($from, $to) = @$targetlocations[$i, $i+1]; 
	$transcript .= &subseq ($targetname, $from, $to-$from, $complement);
    }
    my $product = $ttable->translate($transcript);
    my @result = ();
    for (my $i=0; $i<@$querylocations; $i+=2) {
	my ($from, $to) = map { &nu_to_aa($_) } @$querylocations[$i, $i+1];
	my $partial_product = substr($product, 0, $to-$from, "");
	my $diff_str = &diff_str($partial_product, $queries{$queryname}, $to-$from, 0, $from);
	### DEBUG part ###
	unless (defined $diff_str) {
	    print STDERR "query:$queryname\n";
	    print STDERR "$targetname:$transcript\n";
	    &run_bug("diff_str returned undef!");
	}
	### end of DEBUG part ###
	push @result, "$partial_product:$diff_str";
    }
    return @result;
}
    

# Two hits on different targets are incompatible 
# - if they overlap more than MAX_OVERLAP
# Two hits on the same target are incompatible
# - if they overlap, or are in reversed order or on different strands
# 
sub incompatible {
    my ($hitref1, $hitref2) = @_;
    my $q1from = $hitref1->{"${QPFX}_start"};
    my $q2from = $hitref2->{"${QPFX}_start"};
    ($q2from, $hitref1, $hitref2) = ($q1from, $hitref2, $hitref1)
	if ($q1from > $q2from);
    my ($t1from, $sq1);
    ($q1from, $t1from, $sq1) = @{$hitref1->{matchings}[-1]}{"${QPFX}_start", "${TPFX}_start", "seqshifts"};
    my ($q2to,$t2to, $sq2) = @{$hitref2->{matchings}[0]}{"${QPFX}_end","${TPFX}_end", "seqshifts"};
    my ($q1to, $t1to, $target1, $comp1) = @$hitref1{"${QPFX}_end", "${TPFX}_end", "target", "complement"};
    my ($t2from, $target2, $comp2) = @$hitref2{"${TPFX}_start","target", "complement"};

    my ($length1, $length2) = map { $_->{length}} @dbase{$target1, $target2};
 
    # hits on the same target must not overlap and must be in correct order
    return 1 if (($target1 eq $target2) && ($q1to > $q2from  || $comp1 ne $comp2 || $t1to > $t2from));

    # hits with a huge distance to the sequence end are incompatible to any other
    my ($t1min, $t1end) = $comp1 ? (-$length1,0) : (0,$length1);
    my ($t2start, $t2max) = $comp2 ? (-$length2,0) : (0,$length2);
    return 1 if ($t1end - $t1to >= $ASSEMBLE_SIZE || $t2from - $t2start >= $ASSEMBLE_SIZE);

 
    if ($q1to >= $q2from -1) {
	# there could be overlapping targets; check if they are identical

	# do not allow frameshifts here
	foreach (@$sq1) { return 1 if ($_->{"${QPFX}_end"} >= $q2from) }
	foreach (@$sq2) { return 1 if ($_->{"${QPFX}_start"} <= $q1to) }


	my $offset =  3*($q1to-$q2from) - $t1to + $t2from;
	my $t1start = $t2start - $offset;
    
	# do not allow large overlaps (e.g., duplicate targets)
	my $count = $t1end - $t1start;

	if ($t1start >= $t1min && $t2start + $count <= $t2max && 0 < $count && $count <= $MAX_CHECK_OVERLAP) {
	    my $acceptable_errors = min ($count / 15, $count / 100 +1);
	    my ($subseq1, $subseq2) = (&subseq($target1, $t1start, $count, $comp1), &subseq($target2, $t2start, $count, $comp2));
	    if ( &str_diff_count_all($subseq1, $subseq2) <= $acceptable_errors  ) {
		# yes we found an overlap
		push @{$hitref1->{target_overlaps}}, $hitref2->{ID};
		push @{$hitref2->{target_overlaps}}, $hitref1->{ID};
		$t1start++ unless ($comp1);
		$t2start++ unless ($comp2);
		print COMMENT 
		    "We found an overlap of target sequences:\n A: ${target1}\n B:$target2\n".
		     " A[".abs($t1start)."..".abs($t1start+$count-1)."] ".($comp1?"(comp)":"")." equals\n".
		     " B[".abs($t2start)."..".abs($t2start+$count-1)."] ".($comp2?"(comp)":"")."\n".
		return 0;
	    }
	}    
    }

    # no target overlap: ok if query overlap is no more than MAX_OVERLAP, 
    # and non-overlapping parts have at least the minimum length
    return 
	( $q1to  >  $q2from+$MAX_OVERLAP ||
	  $q1to > $q2to - $MIN_REMAIN ||
	  $q1from > $q2from - $MIN_REMAIN );

}

sub high_margins {
    my ($hit, $neighbor) = @_;
    
    my ($querystart, $queryend, $targetsize, 
	$tailsize, $headsize, $complement) = 
	@{$hit}{"${QPFX}_start", "${QPFX}_end", "target_len",
		"${TPFX}_start", "${TPFX}_end", "complement" };
    if ($complement) {
	$headsize += $targetsize;
    } else {
	$tailsize -= $targetsize;
    }
    my $marginsize = $querystart > $neighbor->{"${QPFX}_start"} ?
	$headsize : -$tailsize;
    return ($queryend-$querystart) < $MIN_SEQ_COVERAGE * $marginsize;
} 

sub replace_in_list {
    my ($from, $to, $listref, $newelems) = @_;
    my $oldsize = @$listref;
    @$listref = ( defined $from ? grep { $_ <= $from } @$listref : (), 
		  @$newelems,
		  defined $to ? grep { $_ > $to } @$listref : ());
    return @$listref - $oldsize;
}
    


sub replace_in_last {
    my ($prevhit, $dnaseq, $protseq, $dna_add, $nsub) = @_;
    $dna_add="" unless (defined $dna_add);
    $nsub=0 unless (defined $nsub);

    my $prevmatch = $prevhit->{matchings}[-1];
    my $oldend = $prevmatch->{"${QNPFX}_end"};
    my $newend = $oldend - $nsub;

    my $psub = &nu_to_aa($oldend) - &nu_to_aa($newend);
    
    my ($ncount, $pcount) = map  length ,($dnaseq, $protseq);
 
    my $translation = substr($ttable->translate($dnaseq.$dna_add),0, $pcount);
    $prevhit->{"${QPFX}_end"} += ($pcount - $psub);
    map { $_ += $ncount-$nsub } (@{$prevmatch}{"${QNPFX}_end","${TPFX}_end"}, 
				 $prevhit->{"${TPFX}_end"});
    substr(substr($prevmatch->{translation},-$psub-1),1) = $translation;
    my %newlists = ( mismatchlist => [ &str_diff_list('X',$translation, $protseq) ],
		     undeterminedlist => [ &str_diff_list(' ',$translation, $protseq) ],
		     inframe_stopcodons => [ &str_pos("[*]", $translation) ] );
    map { $_ += &nu_to_aa($newend) + 1 } @$_ foreach values %newlists;
    my (undef, $mismatchadd, $undetermadd) = 	map {
	&replace_in_list(&nu_to_aa($newend), undef, $prevmatch->{$_}, $newlists{$_}); 
    } sort keys %newlists;

    $prevhit->{mismatches} += $mismatchadd; 
    $prevhit->{undetermined} += $undetermadd;
    $prevhit->{matches} += ($pcount-$psub-$mismatchadd-$undetermadd);
}
    
sub replace_in_next {
    my ($nexthit, $dnaseq, $protseq, $dna_add, $nsub) = @_;
    $dna_add="" unless (defined $dna_add);
    $nsub=0 unless (defined $nsub);

    my $nextmatch = $nexthit->{matchings}[0];
    my $oldstart = $nextmatch->{"${QNPFX}_start"};
    my $newstart = $oldstart + $nsub;
    my $psub = &nu_to_aa($newstart) - &nu_to_aa($oldstart);
    
    my ($ncount, $pcount) = map  length ,($dnaseq, $protseq);
    my $translation = substr($ttable->translate($dna_add.$dnaseq), -$pcount, $pcount);

    $nexthit->{"${QPFX}_start"} -= ($pcount -$psub);
    map { $_ -= $ncount-$nsub } (@{$nextmatch}{"${QNPFX}_start","${TPFX}_start"}, 
				 $nexthit->{"${TPFX}_start"});
    substr($nextmatch->{translation},0,$psub)=$translation;
    my %newlists = ( mismatchlist => [ &str_diff_list('X',$translation, $protseq) ],
		     undeterminedlist => [ &str_diff_list(' ',$translation, $protseq) ],
		     inframe_stopcodons => [ &str_pos("[*]", $translation) ] );
    map { $_ += &nu_to_aa($newstart) - $pcount + 1 } @$_ foreach values %newlists;
    my (undef, $mismatchadd, $undetermadd) = 	map {
	&replace_in_list(undef, &nu_to_aa($newstart),  $nextmatch->{$_}, $newlists{$_}); 
    } sort keys %newlists;

    $nexthit->{mismatches} += $mismatchadd;
    $nexthit->{undetermined} += $undetermadd;
    $nexthit->{matches} += ($pcount-$psub-$mismatchadd-$undetermadd);
}

# this might be needed later when the replace_in_xxx methods are rewritten
#
# sub update_translation {
#     my ($match, $count, $protseq, $dna, $left) = @_;
#     my $translation = $ttable->translate($dna);
#     my %newlists = ( mismatchlist => [ &str_diff_list('X',$translation, $protseq) ],
# 		     undeterminedlist => [ &str_diff_list(' ',$translation, $protseq) ],
# 		     inframe_stopcodons => [ &str_pos("[*]", $translation) ] );
#     my $trref = \$match->{translation};
#     substr($$trref, $left ? 0 : length($$trref)-$count, $count) = $translation;
#     if ($left) {
# 	my $newstart = $match->{"{${QPFX}_start"} + length($protseq);
# 	map { $_ += $newstart - $count +1 } @$_
# 	    foreach values %newlists;
# 	foreach keys %newlists {
# 	    replace_in_list(undef, $newstart, $match->{$_}, $newlists{$_});
# 	}
#     } else {
# 	my $from = $match->{"${QPFX}_end"}-$count;
# 	map { $_ += $from +1 } @$_ foreach values %newlists;
# 	foreach keys %newlists {
# 	    replace_in_list($from, undef, $match->{$_}, $newlists{$_});
# 	}
#     }
# }


sub add_initial {
    my ($hit, $type) = @_;
    my $first_match = $hit->{matchings}[0];
    my $targetname = $hit->{target};
    my $complement = $hit->{complement};
    $hit->{"${TPFX}_start"} = $complement ? -$dbase{$targetname}{length} : 0;
    unshift @{$hit->{matchings}}, { "${TPFX}_start" => $hit->{"${TPFX}_start"},
				    "${TPFX}_end" => $first_match->{"${TPFX}_start"},
				    type => $type,
				    "${QNPFX}_start" => $first_match->{"${QNPFX}_start"},
				    "${QNPFX}_end" => $first_match->{"${QNPFX}_start"} }
}

sub add_final {
    my ($hit, $type) = @_;
    my $last_match = $hit->{matchings}[-1];
    my $targetname = $hit->{target};
    my $complement = $hit->{complement};
    my $newend = $complement ? 0 : $dbase{$targetname}{length};
    $hit->{"${TPFX}_end"} = $newend;
    push @{$hit->{matchings}}, { "${TPFX}_start" => $last_match->{"${TPFX}_end"},
				 "${TPFX}_end" => $newend,
				 type => $type,
				 "${QNPFX}_start" => $last_match->{"${QNPFX}_end"},
				 "${QNPFX}_end" => $last_match->{"${QNPFX}_end"} }
}

sub set_problem_field {
    my $current = shift;
    my @problems;
    my ($qfrom, $qto, $qlen) = @$current{"${QPFX}_start","${QPFX}_end","${QPFX}_len"};
    foreach my $bad_type  (grep { $_ ne "exon" && $_ ne "intron" } 
			    map { $_->{type} } @{$current->{matchings}}) {
	$bad_type = $bad_type eq "intron?" ? "bad intron" : $bad_type;
	push @problems, $bad_type unless (grep { $_ eq $bad_type } @problems);
    }
    push @problems, "sequence shift" 
	if (grep  { (defined $_->{seqshifts} && @{$_->{seqshifts}}) } 
	    @{$current->{matchings}}); 
    push @problems, "mismatches" if ($current->{mismatches});
#   unnecessary because this implies a frameshift or a gap
#    push @problems, "unmatched" if ($current->{unmatched});
    push @problems, "missing stopcodon" if ($qto == $qlen && !$current->{stopcodon});
    push @problems, "in-frame stopcodon" 
	if (grep  { (defined $_->{inframe_stopcodons} && @{$_->{inframe_stopcodons}}) } 
	    @{$current->{matchings}}); 
    push @problems, "gap to querystart" if (defined $current->{upstream} && $qfrom>0);
    push @problems, "gap to previous hit" if (defined $current->{upstream_gap});
    push @problems, "gap to next hit" if (defined $current->{downstream_gap});
    push @problems, "gap to queryend" if (defined $current->{downstream} && $qto<$qlen);
#   add unmatched/undetermined only unless there are other problems
#   in most cases unmatched implies a frameshift or a gap
    unless (@problems) {
	if ($current->{undetermined}) {
	    push @problems, "undetermined residues";
	} elsif ($current->{unmatched}) {
	    push @problems, "unmatched residues";
	}
    }
    if (@problems) {
	$current->{status} = "incomplete";
	$current->{reason} = join "/",@problems;
    } elsif ($qfrom>0 || $qto < $qlen) {
	$current->{status} = "partial";
    } else {
	$current->{status} = "auto";
    }
}

sub nw_splice {
    my ($type, $patt) = @_;
    return (exists $type->{$patt} ? $type->{$patt}: $NW_ANY_SPL);
}
 
sub nw_minimize {
    my ($target, $tkey, $matrix, $score, $i, $j, $k, $m) = @_;
    my $targetref = $target->{$tkey};
    my $sourceref = $matrix->[$i][$j]{$k};
    return unless (defined $sourceref);
    $score += $sourceref->[0];
    return if (defined $targetref && $targetref->[0] <= $score);
    if ($k ne "") {
	($i, $j)  = @{$sourceref}[1,2];
    }
    @{$target->{$tkey}} = ($score, $i, $j );
}

sub nw_align {
    my ($dnaseq, $prot) = @_;
    my ($n, $m) = map { length } ($dnaseq, $prot);
    my @matrix = ( [] );


    # the dp algorithm
    foreach (0..$m) {
	my $currscore = [ $_ * $NW_PROT_INS, 0, 0 ];
	my $currentry = { "" => $currscore };
	push @{$matrix[0]}, $currentry;
    }
	
    foreach my $i (1..$n) {
	my ($fs_minprev, $ttrans, $stopmalus) = (0, "", 0);
	if ($i>2) {
	    $fs_minprev = $i-2;
	    $ttrans = $ttable->translate(substr($dnaseq, $i-3, 3));
	    $stopmalus = $NW_STP if $ttrans eq  "*";
	} 
	for (my $j=0; $j <= $m; $j++) {
	    my $currentry = {};
	    foreach ($fs_minprev .. $i-1) {
		&nw_minimize($currentry, "", \@matrix, $NW_FS, $_, $j, "");
	    }
	    if ($i>=3) {
		&nw_minimize($currentry, "", \@matrix, $NW_PROT_GAP + $stopmalus, $i-3, $j, "");
	    }
	    foreach ("I0",
		     "I1a","I1c","I1g","I1t",
		     "I2a","I2c","I2g","I2t") {
		&nw_minimize($currentry, $_, \@matrix, 0, $i-1, $j, $_);
	    }
	    if ($i>=2) {
#		&nw_minimize($currentry, "", \@matrix, $NW_FS, $i-2, $j, "");
		&nw_minimize($currentry, "", \@matrix, &nw_splice($ASS_PENALTY, substr($dnaseq,$i-2,2)), $i-2, $j, "I0");
		my $dsspos = $i - ($MIN_INTRON_LEN -2);
		if ($dsspos >= 0) {
		    my $dss = &nw_splice($DSS_PENALTY, substr($dnaseq, $dsspos, 2)) + $NW_INTRON;
		    &nw_minimize($currentry, "I0", \@matrix, $dss, $dsspos, $j, "");
		    if ($dsspos>0) {
			&nw_minimize($currentry, "I1".(lc substr($dnaseq, $dsspos-1, 1)), \@matrix, $dss, $dsspos-1, $j, "");
		    }
		    if ($j<$m && $dsspos>1) {
			my $codon = substr($dnaseq, $dsspos-2, 2);
			foreach my $nuc ("a","c","g","t") {
			    die "dss not defined" unless defined($dss);
			    die "NW_MISMATCH not defined" unless defined($NW_MISMATCH);
			    my $codontrans = $ttable->translate("$codon$nuc");
			    my $codonscore = ($codontrans eq uc substr($prot, $j, 1)) ? 0 :
				($NW_MISMATCH + ($codontrans eq "*" ? $NW_STP : 0));
			    &nw_minimize($currentry, "I2$nuc", \@matrix,
					 $dss + $codonscore, $dsspos-2, $j, "");
			}
		    }
		}
	    }
	    if ($j>0) {
		&nw_minimize($currentry, "", \@matrix, $NW_PROT_INS, $i, $j-1, "");
		if ($i>=3) {
		    my $tscore = 0;
		    if ($ttrans ne uc substr($prot, $j-1,1)) {
			$tscore = $NW_MISMATCH + $stopmalus;
		    }
		    &nw_minimize($currentry, "", \@matrix, $tscore, $i-3, $j-1, "");
		}
		foreach ($fs_minprev .. $i-1) {
		    &nw_minimize($currentry, "", \@matrix, $NW_FS, $_, $j-1, "");
		}
		if ($i>=3)  {
		    &nw_minimize($currentry, "", \@matrix, &nw_splice($ASS_PENALTY, substr($dnaseq, $i-3, 2)), 
				   $i-3, $j-1, "I2".(lc substr($dnaseq, $i-1, 1)));
		} 
		if ($i>=4) {
		    my $codon = substr($dnaseq, $i-2, 2);
		    foreach my $nuc ("a", "c", "g", "t") {
			my $codontrans = $ttable->translate("$nuc$codon");
			my $codonscore = ($codontrans eq uc substr($prot, $j-1, 1)) ? 0 :
			    ($NW_MISMATCH + ($codontrans eq "*" ? $NW_STP : 0));
			&nw_minimize($currentry, "", \@matrix,
				       &nw_splice($ASS_PENALTY, substr($dnaseq, $i-4,2)) + $codonscore,
				       $i-4, $j-1, "I1$nuc");
		    }
		}
	    }
	    $matrix[$i][$j] = $currentry;
	}
    }

    # backtracking
    my @result = ();
    my ($tpos, $qpos, $size) = ($n, $m, 0);
    my $endscore = $matrix[$n][$m]{""}[0];
    while ($tpos >  0 || $qpos > 0) {
	my $next = $matrix[$tpos][$qpos]{""};
	my ($score,$nexttpos, $nextqpos) = @$next;
	my $newsize = $qpos - $nextqpos;
	if ($newsize == 1 && $nexttpos == $tpos - 3) {
	    # match or mismatch
	    # we enlarge the matching by 1
	    $size++;
	} else {
	    # we reached the start of a matching
	    # at tpos/qpos
	    if ($size > 0) {
		unshift (@result, [$tpos, $qpos, $size, $endscore - $score]);
	    }
	    # if we have an intron splicing a codon, count this codon for the 
	    # preceding match
	    $size= ($newsize==1 && $nexttpos < $tpos - $MIN_INTRON_LEN) ? 1: 0;
	    $endscore = $score;
	}
	($tpos, $qpos) = ($nexttpos, $nextqpos);
    }
    unshift (@result, [0, 0, $size, $endscore]);
 
    ($tpos, $qpos, $size) = @{$result[-1]};
    if ($tpos + 3*$size == $n && $qpos + $size == $m) {
	return @result;
    }
    return (@result, [ $n, $m, 0, 0 ] );
}

	  
 

################################################################################
