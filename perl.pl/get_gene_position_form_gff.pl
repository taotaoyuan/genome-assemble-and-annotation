open IN,"$ARGV[0]" or die "$!";
open OUT,">$ARGV[1]" or die "$!";
while(<IN>){
		chomp;
		
		next if /^#/;
		@tmp=split(/\t/);
		if($tmp[2]=~/gene/ && $tmp[0]=~/^\d+/ && $tmp[-1]=~/protein_coding/){
				my($id)=($tmp[-1]=~/ID=([^;]+)/);
				print OUT "$tmp[0]\t$id\t$tmp[3]\t$tmp[4]\n";
		}
}
close(IN);
close(OUT);
