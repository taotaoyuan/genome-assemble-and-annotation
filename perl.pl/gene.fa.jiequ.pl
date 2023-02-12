#!/usr/bin perl
die "perl $0 <genome.fa> <weizhi.txt> <OUT> " unless(@ARGV==3 );
use Math::BigFloat;
use Bio::SeqIO;
use Bio::Seq;
$in = Bio::SeqIO -> new(-file => "$ARGV[0]",
                                  -format => 'Fasta');
$out = Bio::SeqIO -> new(-file => ">$ARGV[2]",
                                  -format => 'Fasta');
my %keep=() ;
open IN,"$ARGV[0]" or die "$!";
my%ref=();
while ( my $seq = $in->next_seq() ) {
     my($id,$sequence,$desc)=($seq->id,$seq->seq,$seq->desc);

         $ref{$id}=$seq;

}

$in->close();
open IN,"$ARGV[1]" or die "$!";
while (<IN>) {
     chomp;
     next if /^#/;
     my @a= split /\t/;
     my$seq=$ref{$a[0]};

     my$seq_string=$seq->subseq($a[1],$a[2]);
     my$newseqobj1=Bio::Seq -> new(-seq => $seq_string,
                     -id => "$a[4]"
     ) ;
     if( $a[3]  eq "-" ){
               my$reseq = $newseqobj1 ->revcom();
               $out->write_seq($reseq);
     }elsif ( $a[3]  eq "+" ){

              $out->write_seq($newseqobj1);
     }

}
