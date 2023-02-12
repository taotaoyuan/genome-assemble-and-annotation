use Data::Dumper;
use Getopt::Long;
use strict;
use Cwd qw(abs_path getcwd);


my %opts;

GetOptions (\%opts,"id=s","tandem=s","od=s","name=s"); 

my $od=$opts{od};
$od||=getcwd;
$od=abs_path($od);
unless(-d $od){    mkdir $od;}

####get target gene id

my $gene;
my @info;
my %hashG;
open (IN,"$opts{id}") || die "open $opts{id} failed\n";
while(<IN>){
    chomp;
    @info=split(/\s+/,$_);
    $gene=$info[0];
    $hashG{$gene}=$gene;
}
close(IN);


#######select tandem


my $Agene;
my $Bgene;
open(OUT,">$od/$opts{name}.tandem")||die "open $od/$opts{name}.tandem failed\n";
open (IN,"$opts{tandem}") || die "open $opts{tandem} failed\n";
while(<IN>){
    chomp;
    @info=split(/,/,$_);
    $Agene=$info[0];
    $Bgene=$info[1];
    if(exists $hashG{$Agene} && exists $hashG{$Bgene}){
        print OUT $Agene."\t".$Bgene."\n";
    }

}
close(IN);
close(OUT);