#!/usr/bin/perl
use strict;
use warnings;

# usage: perl get_intron_size.pl gff3_file >output

my $input=$ARGV[0];
my ($eachline,@exons);
my $first=0;
open (IN, "<$input") or die ("no such file!");
while(defined($eachline=<IN>)){
  if($eachline=~/\tmRNA\t/){
    $first++;
    if($first != 1){
      print_intron(@exons);
      @exons=();
    }
  }elsif($eachline=~/\tCDS\t/){
    my @eachline=split(/\t/,$eachline);
    push (@exons, $eachline[3],$eachline[4]);
  }
}
print_intron(@exons);

sub print_intron{
  my (@exons)=@_;
  if(scalar(@exons)>2){
    my @ordered_exons=sort {$a<=>$b} @exons;
    for (my $i=1;$i<=scalar(@ordered_exons)-3;$i=$i+2){
      my $each_intron_size=$ordered_exons[$i+1]-$ordered_exons[$i]-1;
      print "$each_intron_size\t";
    }
  }else{print "0";}
  print "\n";
}