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
while(<GTF>){
	chomp ;
	next if($_=~m/^#/);
	my ($chr, $type, $start, $end, $strand,$attr) = (split /\t/, $_)[0,2,3,4,6,8];
	next if($type ne $TYPE);
	#
	my $tid = $1 if($attr=~m/transcript_id "(.*?)";/);
	my $gid = $1 if($attr=~m/gene_id "(.*?)";/);
	if($tid and $gid and not exists $transcript{$tid}){
		push @transcript_id, $tid;
		$transcript{$tid}{'gene'}   = $gid;
		$transcript{$tid}{'chr'}    = $chr;
		$transcript{$tid}{'strand'} = $strand;
		$transcript{$tid}{'start'}  = $start;
		$transcript{$tid}{'end'}    = $end;
		$transcript{$tid}{$TYPE}    = "$start-$end";
	} elsif($transcript{$tid}){
		$transcript{$tid}{'start'} = $start if($start < $transcript{$tid}{'start'});
		$transcript{$tid}{'end'}   = $end   if($end   > $transcript{$tid}{'end'});
		$transcript{$tid}{$TYPE} .= "; $start-$end";
	}
}
close GTF;

#
open BED, ">$file_bed" or die "";
foreach my $tid (@transcript_id){
	my $chr           = $transcript{$tid}{'chr'};
	my $start         = $transcript{$tid}{'start'} - 1;
	my $end           = $transcript{$tid}{'end'};
	my $gene_id       = $transcript{$tid}{'gene'};
	my $transcript_id = $tid;
	my $strand        = $transcript{$tid}{'strand'};
	my $exons         = $transcript{$tid}{$TYPE};
	#
	my @array = split /; /, $exons;
	my $exon_length = 0;
	my @exon_length = ();
	#
	for(my $i=0;$i<scalar @array;$i++){
		my ($exon_start, $exon_end) = split /-/, $array[$i];
		my $len = abs($exon_end - $exon_start) + 1;
		$exon_length += $len;
		push @exon_length, $len;
	}
	my $exon_length_str = join "; ", @exon_length;
	my $exon_number     = scalar @exon_length;
	#
	print BED "$chr\t$start\t$end\t$gene_id\t$transcript_id\t$strand\t$exon_length\t$exon_number\t$exons\t$exon_length_str\n";
}
close BED;

