#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use 5.010;

# usage
my $thisScript = basename $0;
die "Usage: perl $thisScript in.gtf out.bed\n" if(not @ARGV);

# command
my $file_gtf = shift @ARGV;
my $file_bed = shift @ARGV;

# global options
my $TYPE = "exon";

#
my %transcript    = ();
my @transcript_id = ();

#
open GTF, $file_gtf or die "";
open BED, ">$file_bed" or die "";
while(<GTF>){
	chomp ;
	next if($_=~m/^#/);
	my ($chr, $type, $start, $end, $strand, $attr) = (split /\t/, $_)[0,2,3,4,6,8];
	next if($type ne $TYPE);
	#
	my $tid = $1 if($attr=~m/transcript_id "(.*?)";/);
	my $gid = $1 if($attr=~m/gene_id "(.*?)";/);
	$start = $start - 1;
	print BED "$chr\t$start\t$end\t$gid\t$tid\t$strand\n" if($tid and $gid);
}
close GTF;
close BED;

