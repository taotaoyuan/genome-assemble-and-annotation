use Getopt::Long;
use strict;
use Cwd qw(abs_path getcwd);




my %opts;

GetOptions (\%opts,"list=s","od=s","colline=s","gff=s","name=s"); 



my $od=$opts{od};
$od||=getcwd;
$od=abs_path($od);
unless(-d $od){	mkdir $od;}



#############gff for cir text###########3




open (IN,"$opts{gff}") || die "open $opts{gff} failed\n";
my %gff;
my @info;
my $chr;
my $start;
my $end;
my $gene;
while(<IN>){
	chomp;
	next if /^#/;
	
	@info=split(/\t/,$_);
	
	#next unless($info[2]=~/gene/);
	#($gene)=($info[8]=~/ID=([^;]+)/);
	
	$gene=$info[1];
	$chr=$info[0];
	$start=$info[2];
	$end=$info[3];
	$gff{$gene}=$chr."\t".$start."\t".$end;
}

close(IN);



####################### list ##############


my %list;
my $pair;
my $Len;
my $Agene;
my $Bgene;
my %text;

open (IN,"$opts{list}") || die "open $opts{list} failed\n";
open (OUT,">$od/$opts{name}.txt") || die "open $od/$opts{name}.txt failed\n";
open (OUTL,">$od/$opts{name}.link.txt") || die "open $od/$opts{name}.link.txt failed\n";

while(<IN>){
	chomp;
	@info=split(/\t/,$_);
	$Len = @info;
	print $Len;
	
	my $len=1;
	while($len<$Len){
			$pair=$info[$len];   ####�ų�0��Ҳ���ǵ�һ��λ�ã��˴�Ϊ����
			($Agene,$Bgene)=split(/:/,$pair,2);
			print OUT $Agene."\t".$Bgene."\n";
			
	        if(exists $gff{$Agene} && exists $gff{$Bgene}){ #######ʵ���Ͽ϶����ڣ�
	        	
	        	print OUTL $gff{$Agene}."\t".$gff{$Bgene}."\n";
	        	
	        	$text{$Agene}=$gff{$Agene}."\t".$Agene;   ##��ȥ�ؿ���ֱ�ӿ�ʼ���
	        	$text{$Bgene}=$gff{$Bgene}."\t".$Bgene;
	        	
	        }
	        $len=$len+1;
	}

}
close(IN);
close(OUT);


##ȥ���ظ�ID��text
open (OUT,">$od/$opts{name}.text.txt") || die "open $od/$opts{name}.text.txt failed\n";
my $loc;
while(($gene,$loc)=each %text){
	print OUT $loc."\n";
}
close(OUT);

######### collinearity for genome block colline #####

open (IN,"$opts{colline}") || die "open $opts{colline} failed\n";
open (OUT,">$od/genome.blocklink.txt") || die "open $od/genome.blocklink.txt failed\n";
open (OUTA,">$od/genome.align.blocklink.txt") || die "open $od/genome.align.blocklink.txt failed\n";
my $n;
my $align;
my $colline;
my %block;
my $Agene1S;
my $AgeneNE;
my $Bgene1S;
my $BgeneNE;
my $Achr;
my $Bchr;
while(<IN>){
	chomp;
	if(/^#/){
		if(/Alignment/){
			$n=1;
			$_=~/Alignment ([^:]*)/;
			$align="Alignment".$1;
		}
		next;
	}
	
	$colline=$_;
	@info=split("\t",$colline);
	$Agene=$info[1];
	$Bgene=$info[2];
	
	if(exists $gff{$Agene} && exists $gff{$Bgene} ){
	
	    if($n ==1 ){
		
		    ($chr,$start,$end)=split(/\t/,$gff{$Agene});
	    	$Agene1S=$start;
		    $Achr=$chr;
		
		    ($chr,$start,$end)=split(/\t/,$gff{$Bgene});
		    $Bgene1S=$start;
		    $Bchr=$chr;
		
	    }else{
		
				
		    ($chr,$start,$end)=split(/\t/,$gff{$Agene});
		    $AgeneNE=$end;
		
		    ($chr,$start,$end)=split(/\t/,$gff{$Bgene});
		    $BgeneNE=$end;
		
	    }
	}
	$n=$n+1;
	$block{$align}=$Achr."\t".$Agene1S."\t".$AgeneNE."\t".$Bchr."\t".$Bgene1S."\t".$BgeneNE;	
		
	
}

close(IN);

my $block_info;

while(($align,$block_info)=each %block){
	print OUT $block_info."\n";
	print OUTA $align."\t".$block_info."\n";
}
close(OUT);
close(OUTA);






