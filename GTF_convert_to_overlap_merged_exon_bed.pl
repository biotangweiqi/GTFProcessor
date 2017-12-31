#!/usr/bin/perl
# author: biotang
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
my %genes    = ();
my @gene_ids = ();
foreach my $tid (@transcript_id){
	my $gid = $transcript{$tid}{'gene'};
	if($gid and not exists $genes{$gid}){
		push @gene_ids, $gid;
		$genes{$gid}{'chr'}    = $transcript{$tid}{'chr'};
		$genes{$gid}{'start'}  = $transcript{$tid}{'start'};
		$genes{$gid}{'end'}    = $transcript{$tid}{'end'};
		$genes{$gid}{'strand'} = $transcript{$tid}{'strand'};
		$genes{$gid}{'transcripts'} = $tid;
		$genes{$gid}{$TYPE} = $transcript{$tid}{$TYPE};
	} elsif(exists $genes{$gid}){
		$genes{$gid}{'transcripts'} .= "; $tid";
		$genes{$gid}{$TYPE}         .= "; $transcript{$tid}{$TYPE}";
		$genes{$gid}{'start'} = $transcript{$tid}{'start'} if($genes{$gid}{'start'} > $transcript{$tid}{'start'});
		$genes{$gid}{'end'}   = $transcript{$tid}{'end'}   if($genes{$gid}{'end'}   < $transcript{$tid}{'end'});
	}
}

# output
open BED, ">$file_bed" or die "";
foreach my $gid (@gene_ids){
	my $chr            = $genes{$gid}{'chr'};
	my $gene_id        = $gid;
	my $strand         = $genes{$gid}{'strand'};
	#
	my @segments = split /; /, $genes{$gid}{$TYPE};
	my ($arr_start, $arr_end) = &merge_overlap(\@segments);
	#
	for(my $i=0;$i<scalar @$arr_start;$i++){
		my $start = $$arr_start[$i] - 1;
		my $end   = $$arr_end[$i];
		print BED "$chr\t$start\t$end\t$gene_id\t$strand\n";
	}
}
close BED;

# functions
sub merge_overlap{
	my $arr_segments = shift @_;
	# get start and end
	my @start = ();
	my @end   = ();
	foreach my $seg (@$arr_segments){
		my ($s, $e) = split /-/, $seg;
		push @start, $s;
		push @end, $e;
	}
	# order and sort
	my $size  = scalar @start;
	my @index = (0 .. $size-1); #print join ",", @index; print "\n";
	@index = sort {$start[$a] <=> $start[$b] or $end[$a] <=> $end[$b]} @index;
	#print join ",", @index; print "\n";
	#
	my @new_start = ();
	my @new_end   = ();
	foreach my $i (@index){
		push @new_start, $start[$i];
		push @new_end,   $end[$i];
		#print "$start[$i]-$end[$i]\n";
	}
	#print "\n";
	@start = @new_start;
	@end   = @new_end;

	# merge overlap
	my @mark_start = ();
	my @mark_end   = ();
	#
	my $mk_s = $start[0];
	my $mk_e = $end[0];
	for(my $i = 0; $i < $size; $i++){
		#print "$i novo: $mk_s-$mk_e\n";
		if($i+1 == $size){
			push @mark_start, $mk_s;
			push @mark_end,   $mk_e;
			#print "record the last one: $mk_s-$mk_e\n";
			last;		
		}
		#
		for(my $j = $i+1;$j<$size;$j++){
			#print "Error:$start[$j]-$end[$j] < $mk_s-$mk_e\n" if($start[$j] < $mk_s);
			if($start[$j] >= $mk_s and $start[$j] <= $mk_e){ #
				next if($end[$j] <= $mk_e); # contained
				$mk_e = $end[$j] if($end[$j] > $mk_e); # overlaped
				#print "$i $j expand: $mk_s-$mk_e\n";
			} 
			elsif($start[$j] > $mk_e){
				#print "$i $j record: $mk_s-$mk_e\n";
				push @mark_start, $mk_s;
				push @mark_end,   $mk_e;
				# new segment
				$mk_s = $start[$j];
				$mk_e = $end[$j];
				$i    = $j-1;
				last;
			}
		}
	}
	#
	return(\@mark_start, \@mark_end);
}

