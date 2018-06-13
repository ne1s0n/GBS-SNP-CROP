#!/usr/bin/perl

#########################################################################################################################
# GBS-SNP-CROP, Step 5. For description, please see Melo et al. (2016) BMC Bioinformatics DOI 10.1186/s12859-016-0879-y
##########################################################################################################################

##########################################################################################
# Requirement 1: BWA aligner (Li & Durbin, 2009)
# Requirement 2: SAMTools (Li et al., 2009)
##########################################################################################

use strict;
no warnings 'uninitialized';
use Getopt::Long qw(GetOptions);
use Parallel::ForkManager;

my $Usage = "Usage: perl GBS-SNP-CROP-5.pl -d <data type, PE = Paired-End or SR = Single-End> -b <barcode-ID file name>  -ref <reference FASTA file> -Q <Phred score> -q <mapping quality score>\n"
." -f <SAMTools -f flag> -F <SAMTools _F flag> -t <threads> -Opt <any additional desired SAMTools options>.\n"
."-ir <Index Reference genome, either 'TRUE' or 'FALSE' (true is default)>\n";

my $Manual = "Please see UserManual on GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($dataType,$barcodesID_file,$Reference,$phred_Q,$map_q,$f,$F,$threads,$sam_add,$index_reference);

GetOptions(
'd=s' => \$dataType,          	# string - "PE" or "SE"
'b=s' => \$barcodesID_file,     # file
'ref=s' => \$Reference,         # file
'Q=s' => \$phred_Q,             # numeric
'q=s' => \$map_q,               # numeric
'f=s' => \$f,               	  # numeric 
'F=s' => \$F,                	  # numeric
't=s' => \$threads,             # numeric
'Opt=s' => \$sam_add,           # string
'ir=s' => \$index_reference,    # string - "TRUE" or "FALSE"
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 5, v.3.0\n#################################\n";
my $pm = new Parallel::ForkManager($threads);
my $sttime = time;

my @files = ();
open my $BAR, "<", "$barcodesID_file" or die "Can't find barcode_ID file\n";
while(<$BAR>) {
	my $barcodesID = $_;
	chomp $barcodesID;
	my @barcode = split("\t", $barcodesID);
	my $barcode_list = $barcode[0];
	my $TaxaNames = $barcode[1];
	push @files, $TaxaNames;
}
close $BAR;
chomp (@files);

#ensuring the value of $index_reference flag
if (! defined $index_reference){
	$index_reference = "TRUE";
}
$index_reference = uc($index_reference);
if ($index_reference eq "TRUE"){
	$index_reference = 1;
}elsif ($index_reference eq "FALSE"){
	$index_reference = 0;
}else{
	die ("Parameter -ir can only accept values 'TRUE' and 'FALSE'\n");
}

#ensuring the value of $dataType flag
if (! defined $dataType){
	die ("Parameter -d is required, allowed values: 'PE' and 'SE'\n");
}
$dataType = uc($dataType);
if ($dataType ne 'SE' and $dataType ne 'PE'){
	die ("Parameter -d is required, allowed values: 'PE' and 'SE'\n");
}

sub main {
   	my $dir = "alignments";
  	unless(-e $dir, or mkdir $dir) {die "Directory $dir does not exist and cannot be created.\n";}
}   
main();

############################
# 1. Reference genome indexing (BWA and samtools)
############################
if ($index_reference){
	# 1.1 Index
	print "\nIndexing reference FASTA file ...";
	system ( "bwa index -a bwtsw $Reference" );
	print "\nDONE.\n\n";
		
	# 1.2 Index reference FASTA file
	print "\nIndexing the reference genome FASTA file ...";
	system ( "samtools faidx $Reference" );
	print "\nDONE.\n\n";
}else{
	print "\nSkipping indexing of reference genome\n";
}

foreach my $file (@files) {
	print "\nProcessing sample $file\n";
	#########################################
	# Mapping on reference genome via BWA-mem
	#########################################
	my $BWA_out = join(".","$file","sam");
	my $input_R1 = join (".", "$file","R1","fq","gz");
	my $input_R2 = '';
	if ($dataType eq "PE") {
		$input_R2 = join (".", "$file","R2","fq","gz");
	}
	print " - mapping paired $input_R1 $input_R2 file(s) to $Reference ...";
	system ( "bwa mem -t $threads -M $Reference $input_R1 $input_R2 > $BWA_out" );
	print "DONE.\n";

	#########################################
	# SAMTools procedures
	#########################################

	#--------------- 2.1 SAM to BAM
	print " - SAM to BAM...";
	my $input_sam = join (".", "$file","sam");
	my $view_out = join(".","$file","bam");
	#	print "\nProcessing $input_sam file ...";

	if ($F > 0 && $f > 0 && ($sam_add ne '0') ) {
		system ( "samtools view -@ $threads -b -q$phred_Q -f$f -F$F $sam_add $input_sam > $view_out" );
	} elsif ($F > 0 && $f == 0 && ($sam_add ne '0') ) {
		system ( "samtools view -@ $threads -b -q$phred_Q -F$F $sam_add $input_sam > $view_out" );
	} elsif ($f > 0 && $F == 0 && ($sam_add ne '0') ) {
		system ( "samtools view -@ $threads -b -q$phred_Q -f$f $sam_add $input_sam > $view_out" );
	} elsif ($f > 0 && $F > 0 && ($sam_add eq '0') ) {
		system ( "samtools view -@ $threads -b -q$phred_Q -f$f -F$F $input_sam > $view_out" );
	} elsif ($F > 0 && $f == 0 && ($sam_add eq '0') ) {
		system ( "samtools view -@ $threads -b -q$phred_Q -F$F $input_sam > $view_out" );
	} elsif ($f > 0 && $F == 0 && ($sam_add eq '0') ) {
		system ( "samtools view -@ $threads -b -q$phred_Q -f$f $input_sam > $view_out" );
	} else {
		print "Unable to proceeed; please re-check the syntax of all declared SAMTools flags and options...";
	}
		
	#removing SAM file
	unlink ($input_sam);
	
	print "DONE.\n";

	#--------------- 2.2 Sorting BAM files
	print " - sorting BAM file...";
	my $input_bam = join (".", "$file","bam");
	my $sort_out = join(".","$file","sorted.bam");
	#	print "\nSorting $input_bam file ...";
	system ( "samtools sort -@ $threads $input_bam -o $sort_out" );
		
	#removing unsorted BAM file
	unlink ($input_bam);

	print "DONE.\n";

	#--------------- 2.3 Index sorted BAM files
	print " - indexing the sorted BAM file...";
	my $input_sorted = join (".","$file","sorted","bam");
	#	print "\nIndexing $input_sorted file ...";
	system ( "samtools index -@ $threads -c $input_sorted" );
	print "DONE.\n";

	#--------------- 2.4 Mpileup SNPs discovery
	print " - producing the mpileup file...";
	my $input = join (".", "$file","sorted","bam");
	my $mpileup = join (".", "$file","mpileup");
	#	print "Producing mpileup file from $file ...\n";
	system ("samtools mpileup -Q$phred_Q -q$map_q -B -C 50 -f $Reference $input > $mpileup");
		
	#removing sorted BAM file
	unlink ($input);
	print "DONE.\n";
		
	#--------------- 2.5 Mpileup compression
	print " - compressing the final mpileup...";
	system("gzip $mpileup");
	print "DONE.\n";
}

system ( "mv *bam* ./alignments" );
system ( "rm *.sam" );
print "Elapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
