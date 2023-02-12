#!/usr/bin/perl
#author   :Gou Xiang Jian (862137261@qq.com)
#date     :2018/06/18
#function :Translation of DNA sequence into protein sequence by using a six-frame translation method
#use      :perl six_frame_tran_method.pl <fasta file>
#test     :perl six_frame_tran_method.pl gffread_transcripts-info-cpc.fa
#result   :run time is about 20 seconds (fasta file is 16.4 Mb, 47104 lines, 23552 sequences)
use strict;
use warnings;
my @id_seq;
{
    local $/;
    my $file_info = <>; #read file from the command line
    foreach(split />/,$file_info){
        if($_){
            my ($id,$seq) = /(.*)\n([\d\D]*)\Z/;
            $seq =~ s/\n//g;
            push @id_seq,{$id => $seq};
        }
    }
}
&produce_six_frame($_) foreach @id_seq;
my %codon_code = &codon_list;
$ARGV =~ s/(.*)\.(fa(?:sta)?)\z/$1-pro.$2/ or die qq(the suffix of the input file isn't fa/fasta.\n);
open OUT,'>',$ARGV or die qq(can't generate $ARGV:$!);
foreach my $ref (@id_seq){
    print OUT ">$_\n",&trans_to_protein($ref->{$_}),"\n" foreach sort keys %$ref;
}
close OUT;
sub produce_six_frame{
    my $ref = shift;
    my ($key) = keys %$ref;
    my $value = $ref->{$key};
    delete $ref->{$key};
    my $an_seq = reverse $value;
    $an_seq =~ tr/AGTCagtc/TCAGtcag/;
    foreach(1 .. 6){
        my $id  = "${key}_$_";
        my $seq = $_ <= 3 ? substr $value,$_-1 : substr $an_seq,$_-4;
        $ref->{$id} = $seq;
    }
}
sub trans_to_protein{
    my $seq = uc $_[0];
    my $num = length($seq) % 3 ? int( length($seq)/3 )+1 : length($seq)/3;
    my ($sub,$protein_seq) = (0,'');
    foreach(1 .. $num){
        my $codon = $_ != $num ? substr $seq,$sub,3 : substr $seq,$sub;
        $sub += 3;
        if(length $codon == 3){
            if($codon !~ /N/){
                $protein_seq .= $codon_code{$codon};
            }
            else{
                my $n_num = $codon =~ s/N/N/g;
                if($n_num >= 2){
                    $protein_seq .= 'X'; # X indicates any amino acid
                }
                else{
                    my @codon = ();
                    foreach(qw/A T G C/){
                        my $replace = $codon =~ s/N/$_/r;
                        push @codon,$replace;
                    }
                    $protein_seq .= ($codon_code{$codon[0]} eq $codon_code{$codon[1]} and $codon_code{$codon[1]} eq $codon_code{$codon[2]} and  $codon_code{$codon[2]} eq $codon_code{$codon[3]}) ? $codon_code{$codon[0]} : 'X';
                }
            }
        }
        elsif(length $codon == 2){
            if($codon =~ /N/){
                $protein_seq .= 'X';
            }
            else{
                my @codon = ();
                push @codon,$codon.$_ foreach qw/A T G C/;
                $protein_seq .= ($codon_code{$codon[0]} eq $codon_code{$codon[1]} and $codon_code{$codon[1]} eq $codon_code{$codon[2]} and  $codon_code{$codon[2]} eq $codon_code{$codon[3]}) ? $codon_code{$codon[0]} : 'X';
            }
        }
        else{
            $protein_seq .= 'X';
        }
    }
    return $protein_seq;
}
sub codon_list{
    my %temp = (
        GCA => 'A', GCC => 'A', GCG => 'A', GCT => 'A',                             # Ala   4
        TGC => 'C', TGT => 'C',                                                     # Cys   2
        GAC => 'D', GAT => 'D',                                                     # Asp   2
        GAA => 'E', GAG => 'E',                                                     # Glu   2
        TTT => 'F', TTC => 'F',                                                     # Phe   2
        GGA => 'G', GGC => 'G', GGG => 'G', GGT => 'G',                             # Gly   4
        CAC => 'H', CAT => 'H',                                                     # His   2
        ATA => 'I', ATC => 'I', ATT => 'I',                                         # Ile   3
        AAA => 'K', AAG => 'K',                                                     # Lys   2
        TTA => 'L', TTG => 'L', CTA => 'L', CTC => 'L', CTG => 'L', CTT => 'L',     # Leu   6
        ATG => 'M',                                                                 # Met   1
        AAC => 'N', AAT => 'N',                                                     # Asp   2
        CCA => 'P', CCC => 'P', CCG => 'P', CCT => 'P',                             # Pro   4
        CAA => 'Q', CAG => 'Q',                                                     # Glu   2
        AGA => 'R', AGG => 'R', CGA => 'R', CGC => 'R', CGG => 'R', CGT => 'R',     # Arg   6
        AGC => 'S', AGT => 'S', TCA => 'S', TCC => 'S', TCG => 'S', TCT => 'S',     # Ser   6
        ACA => 'T', ACC => 'T', ACG => 'T', ACT => 'T',                             # Thr   4
        GTA => 'V', GTC => 'V', GTG => 'V', GTT => 'V',                             # Val   4
        TGG => 'W',                                                                 # Trp   1
        TAC => 'Y', TAT => 'Y',                                                     # Tyr   2
        TAA => '*', TAG => '*', TGA => '*',                                         # Stop  3
    );
    return %temp;
}
