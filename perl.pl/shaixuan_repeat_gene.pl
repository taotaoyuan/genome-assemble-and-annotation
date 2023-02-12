use Getopt::Long;
my %opts;
use Data::Dumper;
GetOptions (\%opts,"in1=s","in2=s","out=s","h");
if (! defined($opts{in1}) ||! defined($opts{in2})|! defined($opts{out}) || defined($opts{h})){
	&USAGE;
}
open (IN1,"$opts{in1}") || die "open $opts{in} failed\n";
open (IN2,"$opts{in2}") || die "open $opts{ina} failed\n";	
open (OUT,">$opts{out}") || die "open $opts{out} failed\n";
 my %count;
 while(<IN1>){
	chomp;
	my @line = split("\t",$_);
	$count{$line[0]}= $line[1];
	#print "$count{$line[0]\n}";
}

while( <IN2>){

		chomp($_);
		my @line1 =split ("\t",$_);
		#print @line1;
		#print "\n"；
		my $a = $count{$line[0]} > $count{$line1[1]} ? $count{$line1[0]}:$count{$line1[1]};
		if(($line1[0] ne $line1[1]) && ($line1[2] > 75 )&& ($line1[3]*1 > 0.75*$a)){
			print OUT $_."\t$a\n";

		}
		#print $count1{$line1[0]};
		
}


close(IN1);
close(IN2);
close(OUT);
sub USAGE {
		print "usage: perl test1.pl -in1 cds_length -in2 result.txt -out shaixuan_result.txt";
	exit;
}
