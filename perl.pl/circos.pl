#!/usr/bin/perl


use Config::General;
use Data::Dumper;
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
#use PerlIO::gzip;
use feature 'state' ;

# here's where the perceptual Brewer palettes come in handy
# a 9 color diverging spectral palette is used
#color  = spectral-9-div-1,spectral-9-div-2,spectral-9-div-3,spectral-9-div-4,spectral-9-div-5,spectral-9-div-6,spectral-9-div-7,spectral-9-div-8,spectral-9-div-9
my@color_sequential=qw(
blues-6-seq-1
bugn-6-seq-1
bupu-6-seq-1
gnbu-6-seq-1
greens-6-seq-1
greys-6-seq-1
oranges-6-seq-1
orrd-6-seq-1
pubu-6-seq-1
pubugn-6-seq-1
purd-6-seq-1
purples-6-seq-1
rdpu-6-seq-1
reds-6-seq-1
ylgn-6-seq-1
ylgnbu-6-seq-1
ylorbr-6-seq-1
ylorrd-6-seq-1
);

my@color_diverging=qw(
brbg-6-div-1
piyg-6-div-1
prgn-6-div-1
puor-6-div-1
rdbu-6-div-1
rdgy-6-div-1
rdylbu-6-div-1
rdylgn-6-div-1
spectral-6-div-1
);

my @color_qualitative=qw(
accent-6-qual-1
dark2-6-qual-1
paired-6-qual-1
pastel1-6-qual-1
pastel2-6-qual-1
set1-6-qual-1
set2-6-qual-1
set3-6-qual-1
);
my @color_set1=qw(set1-9-qual-1 set1-9-qual-2 set1-9-qual-3 set1-9-qual-4 set1-9-qual-5 set1-9-qual-6 set1-9-qual-7 set1-9-qual-8 set1-9-qual-9);
my @color_set2=qw(set2-8-qual-1 set2-8-qual-2 set2-8-qual-3 set2-8-qual-4 set2-8-qual-5 set2-8-qual-6 set2-8-qual-7 set2-8-qual-8);
my @color_set3=qw(set3-12-qual-1 set3-12-qual-2 set3-12-qual-3 set3-12-qual-4 set3-12-qual-5 set3-12-qual-6 set3-12-qual-7 set3-12-qual-8 set3-12-qual-9 set3-12-qual-10 set3-12-qual-11  set3-12-qual-12);
my @color_paired=qw(paired-12-qual-1 paired-12-qual-2 paired-12-qual-3 paired-12-qual-4 paired-12-qual-5 paired-12-qual-6 paired-12-qual-7 paired-12-qual-8 paired-12-qual-9 paired-12-qual-10 paired-12-qual-11  paired-12-qual-12);

my($chr,$outfile,@circles,$od,@type,$space,@thickness,@link);
#-------------chr initial var-------------------------------------------------#
my($chromosomes_units,$chr_spacing,$chr_radius,$chr_thickness,$chr_fill,$chr_stroke_color,$chr_stroke_thickness,$chr_show_label,
	$chr_label_radius,$chr_label_size,$chr_label_parallel,$gap_thickness,$max_circle_out_pisition,$circle_width)=(1000000,'10u','0.90r','30p','yes','dgrey','2p','yes','60p','45','yes','0.01','0.99','0.1');

#------------ticks of chr initial var-------------------------------------------#
my($show_ticks,$show_tick_labels,$tick_label_size,$tick_label_offset,$tick_format)=('yes','yes','25p','10p','%d');
#type:histogram ,scatter,line,heatmap

sub USAGE {
my$msg=shift;
my $usage=<<"USAGE";
    Name:  
        $0
    Usage:
        draw whole genome info by circos, support line,scatter,heatmap,histogram,link;
    Options:
         --h for help

        --chr <character>  input chr ideogram file , required
        --outfile <character>  output file name,default of "circos"
        --od  <character>  outpath,required
        --circle <character>	input files of each Tracks;
        --type <character>	the type of  plot data Tracks to draw, corresponding to --circle; [line,scatter,heatmap,histogram,highlight,text,link];
        	
        --chromosomes_units	<character>	default	1000000
        --chr_spacing	<character>	default	'10u'
        --chr_radius	<character>	default	'0.90r'
        --chr_thickness	<character>	default	'20p'
        --chr_fill	<character>	default	'yes'
        --chr_stroke_color	<character>	default	'dgrey'
        --chr_stroke_thickness	<character>	default	'2p'
        --chr_show_label	<character>	default	'yes'
        --chr_label_radius	<character>	default	'60p'
        --chr_label_size	<character>	default	'30'
        --chr_label_parallel	<character>	default	'yes'
        --gap_thickness	<character>	default	'0.015'
        --circle_width <character> default 0.1 
        
	--thick <character> default '0.1' for each track,
        --show_ticks	<character>	default	'yes'
        --show_tick_labels	<character>	default	'yes'
        --tick_label_size	<character>	default	'20p'
        --tick_label_offset	<character>	default	'10p'
        --tick_format	<character>	default	'%d'
    
    Usage example:
        (1)perl circos.pl --chr Target.gene.scaffold.length --od ./
        (2)perl circos.pl --chr Target.gene.scaffold.length --od ./  --chromosomes_units 1000 
            --circle R3_R2.Target.gene.scaffold.fst_scatter --type line --circle R3_R2.Target.gene.scaffold.fst_scatter --type histogram 
            --circle R3_R2.Target.gene.scaffold.fst_scatter --type scatter
        (3)perl circos.pl --chr Target.gene.scaffold.length -od ./  --chromosomes_units 1000 --circle R3_R2.Target.gene.scaffold.fst_scatter --type line 
            --circle R3_R2.Target.gene.scaffold.fst_scatter --type histogram --circle R3_R2.Target.gene.scaffold.fst_scatter 
            --type scatter --circle target.gene.pos_hightlight --type highlight --thick 0.1 --thick 0.2 --thick 0.2
    Data format:
        Ref chr:
            format: chr    length;
                hs1 46975000
                hs2 36975000
            
        Heatmaps: 
            The data format for the heatmap is the same as for other plots: chr, start, end, value, and a list of optional parameters.
                ...
                hs7 36975000 36999999 33
                hs7 37000000 37024999 50
                hs7 37025000 37049999 60 color=blue
                hs7 37050000 37074999 44
                ...
        Histograms,scatter,line:
            In line and scatter plots, the data point is placed at the midpoint of the point's span. Thus, if you define a data point
                hs1 1000 2000 0.5
                hs1 5000 5500 0.25
                hs1 9000 9250 0.75
        Highlight:
            Sample highlight input format:chr start end ;
                hs1 1000 2000 
                hs1 5000 5500 
                hs1 9000 9250
        Text format:
                hs1 1000 2000 label1 
                hs1 5000 5500  label2
                hs1 9000 9250 label3
        Link format:
                hs1 1000 2000 hs7 36975000 36999999 
                hs1 5000 5500  hs7 37000000 37024999
                hs1 9000 9250 hs7 37025000 37049999
                
        To get more data format or other information,please go to website: http://circos.ca/tutorials/lessons/;
        
    Author:
        huangls
    Version:
        2.0;2015-4-22
USAGE
    print $usage,"\n";
    print $msg,"\n";
    exit;
}

my $man = 0;
my $help = 0;

GetOptions(
            "help|h|?"=>\$help,
            "man"=>\$man,
            "chr:s"=>\$chr,
            "outfile:s"=>\$outfile,
            "circle:s"=>\@circles,
            "type:s"=>\@type,           
            "space"=>\$space,
            "thick"=>\@thickness,
            "od:s"=>\$od,
           
            
            'chromosomes_units:s'=>\$chromosomes_units,
			'chr_spacing:s'=>\$chr_spacing,
			'chr_radius:s'=>\$chr_radius,
			'chr_thickness:s'=>\$chr_thickness,
			'chr_fill:s'=>\$chr_fill,
			'chr_stroke_color:s'=>\$chr_stroke_color,
			'chr_stroke_thickness:s'=>\$chr_stroke_thickness,
			'chr_show_label:s'=>\$chr_show_label,
			'chr_label_radius:s'=>\$chr_label_radius,
			'chr_label_size:s'=>\$chr_label_size,
			'chr_label_parallel:s'=>\$chr_label_parallel,
			'gap_thickness:s'=>\$gap_thickness,
			'circle_width:s'=>\$circle_width,
			#'max_circle_out_pisition:s'=>\$max_circle_out_pisition,
            
            'show_ticks:s'=>\$show_ticks,
			'show_tick_labels:s'=>\$show_tick_labels,
			'tick_label_size:s'=>\$tick_label_size,
			'tick_label_offset:s'=>\$tick_label_offset,
			'tick_format:s'=>\$tick_format,
            
            ) || pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

&USAGE unless ($chr);
$od||=getcwd();
unless (-d $od) {
	mkdir $od;
}
$od = abs_path($od);
&USAGE unless (scalar(@circles) == scalar(@type));
#--------------------------initial draw type--------------------------#
$outfile||="circos";

chomp $chromosomes_units;
my$len_char=split(//,$chromosomes_units);
my$tick_multiplier='1e-'.($len_char-1);
#print $tick_multiplier;
my $len_circle = @circles;
if($len_circle){
	if($len_circle != @type ){
		&USAGE("please check,--type and --circle must be equal!!!");
	}else{
		my$type_str="line,scatter,heatmap,histogram,highlight,text,link";
		for my$i (@type){
			if ($type_str!~/$i/){
				&USAGE ("Exit: don't support this type of \"$i\" to draw circos \n type of text,link,line,scatter,heatmap,histogram,highlight were supported so far;");
				
			}
		}
	}
	if($len_circle == @thickness){
		
	}else{
		#print "there is too many thickness of circle,extra element will be omited\n";
		@thickness=();
		my$sum;
		for (my $i=0;$i<$len_circle;$i++){
			
			if($type[$i] eq "link"){
				next;
			}else{
				$sum=$sum+$circle_width+$gap_thickness;
				
			}
					
		}
		my$link_r=1-$sum-$gap_thickness;
		for (my $i=0;$i<$len_circle;$i++){
			
			if($type[$i] eq "link"){
				push @thickness,$link_r;
			}else{
				push @thickness,$circle_width;
				
			}
					
		}
		
		
		
		if(1-$sum<= 0){
			&USAGE("too many tracks to draw");
		}

	}
	
}

$chr=&format_chr($chr);




#my $conf = Config::General->new("$infile");
#my %config = $conf->getall;
#print Dumper(\%config);
#---------------------------------------------set initial config hash---------------------------------#
 my %config = (

			karyotype =>$chr,
			chromosomes_units => $chromosomes_units,
			spacing => $chr_spacing,
			ideogram=>{
				spacing=>{
					default => '0.005r'
				},
				# Ideogram position, fill and outline
				radius           => $chr_radius,
				thickness        => $chr_thickness,
				fill             => $chr_fill,
				stroke_color     => $chr_stroke_color,
				stroke_thickness => $chr_stroke_thickness,

				# Minimum definition for ideogram labels.

				show_label       => $chr_show_label,
				# see etc/fonts.conf for list of font names
				label_font       => 'default' ,
				label_radius     => 'dims(image,radius)-'.$chr_label_radius,
				label_size       => $chr_label_size,
				label_parallel   => $chr_label_parallel,	

			},
			show_ticks          => $show_ticks,
			show_tick_labels    =>$show_tick_labels,
			ticks=>{
				radius           => '1r',
				color            => 'black',
				thickness        => '2p',

				# the tick label is derived by multiplying the tick position
				# by 'multiplier' and casting it in 'format':
				#
				# sprintf(format,position*multiplier)
				#

				multiplier       => $tick_multiplier,

				# %d   - integer
				# %f   - float
				# %.1f - float with one decimal
				# %.2f - float with two decimals
				#
				# for other formats, see http://perldoc.perl.org/functions/sprintf.html

				format           => '%d',

				tick=>[{
					spacing        => '5u',
					size           => '10p',
				},
				{
					spacing        => '25u',
					size           => '15p',
					show_label     =>$show_tick_labels,
					label_size     =>$tick_label_size,
					label_offset   => $tick_label_offset,
					format         => $tick_format,
					thickness      => '4p',
					color          =>'black',
				}]																
			},
);
#---------------------------plots flag intial----------------------#           
 my %plots=(
 	
 		type =>'histogram',
 		#fill_color => 'vdgrey',
 		#color =>'vdgrey',
 	);
#--------------------------highlights plot ---------------------------#
#<highlights>

#z          = 0
#fill_color = green

#<highlight>
#file       = data/3/random.highlights.txt
#r0         = 0.8r
#r1         = 0.8r + 200p
#</highlight>

#<highlight>
#file       = data/3/random.highlights.z_by_size.txt
#r0         = 0.6r
#r1         = 0.6r + 200p
#stroke_thickness = 4
#stroke_color     = white
#</highlight>

#<highlight>
#file             = data/3/random.highlights.zr_by_size.txt
#stroke_thickness = 2
#stroke_color     = black
#</highlight>

#</highlights>


my %highlights=(
	z          => 0,
	fill_color => 'green',
);
my%links=(
#set globl default value
	ribbon           => "no",
	color            => "black",
	thickness        => "1",
	radius           => "0.40r",
	bezier_radius    => "0r",
	crest                => 0.5,
	bezier_radius_purity => 0.75
);
#--------------------------set axes-------------------------#
           
my %axes=(
# Show axes only on ideograms that have data for this track

axis=>{
	spacing   => '0.2r',
	color     => 'grey',
	thickness => '1',
	},	
);
my %backgrounds=(
 				background=>{
 				color => 'vvlgrey'
 				},
 			 );
my %one_link_rules=(
				rule=>[{
						condition  => "var(intrachr)",
   						color       => "red"
				},{
						condition  => "var(interchr)",
   						color       => "green"
				}],

);
#--------------------------------------set plots--------------------------------------
my $color='black';

if (@circles){
	
	my @R=&get_circle_position($gap_thickness,\@thickness);
	my($a)=($R[-1]=~/^(0\.\d+)r$/);
	while($a< 0.3){
		for (my $i=0;$i<@thickness;$i++){
			$thickness[$i]=$thickness[$i]-0.001;
		}
		@R=&get_circle_position($gap_thickness,\@thickness);
		($a)=($R[-1]=~/^(0\.\d+)r$/);
	}
	
	for (my $circle=0;$circle<@circles;$circle++){
		
		my $plot_id="plot";
		my $highlight_id="highlight";
		my $r0=shift @R;
		my $r1=shift @R;
		$circles[$circle]=abs_path($circles[$circle]);
		if($type[$circle] eq 'scatter'){
			my($min,$max)=&get_min_max($circles[$circle]);
			my %plot=&get_plot($type[$circle],$color_paired[$circle],$circles[$circle],2,$r0,$r1,$min,$max);
			
			$plot{'axes'}={%axes};
			$plots{'backgrounds'}={%backgrounds};
			push @{$plots{$plot_id}},{%plot};
			
		}elsif($type[$circle] eq 'line'){
			my($min,$max)=&get_min_max($circles[$circle]);
            my %plot=&get_plot($type[$circle],$color_set2[$circle],$circles[$circle],2,$r0,$r1,$min,$max);
            $plot{'axes'}={%axes};
            $plot{'backgrounds'}={%backgrounds};
            push @{$plots{$plot_id}},{%plot};
			
		}elsif($type[$circle] eq 'highlight'){
			 my%highlight= &get_highlight($circles[$circle],$r0,$r1,4,$color_set2[$circle],$color_set2[$circle]);
			 push @{$highlights{$highlight_id}} ,{%highlight} ;
			
		}elsif($type[$circle] eq 'heatmap'){
			my %plot=&get_plot($type[$circle],$color_set1[$circle],$circles[$circle],2,$r0,$r1);
			push @{$plots{$plot_id}},{%plot};
			#$plots{$plot_id}{'axes'}={%axes};
			#$plots{$plot_id}{'backgrounds'}={%backgrounds};
		}elsif($type[$circle] eq 'histogram'){
			my($min,$max)=&get_min_max($circles[$circle]);
			my %plot=&get_plot($type[$circle],$color_set1[$circle],$circles[$circle],2,$r0,$r1,$min,$max);
			
			$plot{'axes'}={%axes};
			#$plots{'backgrounds'}={%backgrounds};	
			push @{$plots{$plot_id}},{%plot};		
		
		}elsif($type[$circle] eq 'text'){
			my %plot=&get_plot($type[$circle],$color_set2[$circle],$circles[$circle],2,$r0,$r1);
			#$plot{'axes'}={%axes};
			#$plot{'backgrounds'}={%backgrounds};
			push @{$plots{$plot_id}},{%plot};
			
		}elsif($type[$circle] eq 'link'){
			#warn "haven't add link function";
			#delete $type[$circle];
			#delete $circles[$circle];
			#$file,$r0,$fill_color
			my%link= &get_link($circles[$circle],$r0,$color_set2[$circle]);
			push @{$links{"link"}},{%link} ;
			
		}else{
			die "Exit: don't support this type of \"$type[$circle]\" to draw circos \n type of text,link,line,scatter,heatmap,histogram,highlight were supported so far;";
			delete $type[$circle];
			delete $circles[$circle];
		}			
	}
}






#----------------------------------plot  config----------------------------------------------#
if("@type"=~/link/){
	if(@{$links{"link"}}==1){
		push @{$links{"link"}->[0]->{"rules"}},{%one_link_rules};
	}
	$config{'links'} = {%links};
}
if("@type"=~/line|scatter|heatmap|histogram|highlight|text/){
	$config{'plots'} = {%plots};
}
if("@type"=~/highlight/){
	$config{'highlights'} = {%highlights};
}
my$conf = Config::General->new(-SplitPolicy=> 'custom' , -SplitDelimiter => '\s*=\s*',-StoreDelimiter=>'=',-SaveSorted => 1);
$conf->save_file("$od/$outfile.conf", \%config);

&add_include();
#print Dumper(\%config);




#-----------------------sub ---------------------------------------------------------#


sub get_plot{
	my($type,$color,$file,$thickness,$r0,$r1,$min,$max)=@_;
	$file=abs_path($file);
	my %plot_demo;
	if ($type eq 'scatter'){
		 %plot_demo=(	
 			type    => $type,
			#color => 'vdgrey',
			color => $color,
			min => $min,
			max => $max,
			#thickness     => 2,
			thickness     => $thickness,
			file             => $file,
			#r0               => '0.90r',
			#r1               => '0.98r',
 			r0               => $r0,
			r1               => $r1,
		
			);			
	}elsif($type eq 'line' ){
		%plot_demo=(
			type      => $type,
			thickness => $thickness,
			max_gap => '1u',
			file    => $file,
			color   => $color,
			min     => $min,
			max     => $max,
			r0      => $r0,
			r1      => $r1,
			#fill_color = vdgrey_a3
			#orientation = in
		);
	}elsif($type eq 'histogram'){
		%plot_demo=(
			type      => $type,
			file      => $file,
			r1        => $r1,
			r0        => $r0,
			max       => $max,
			min       => $min,
			stroke_type => 'both',# outline | bin | both
			fill_color=>$color,
			thickness   => $thickness,
			color       => $color,
			extend_bin  => 'no',
			#fill_color       = red
			#background_color = lred
		);
	}elsif($type eq 'heatmap'){
		%plot_demo=(	
			
			type    => $type,
			file    => $file,
			# color list
			color   => 'spectral-5-div',
			r1      => $r1,
			r0      => $r0,);
	}elsif($type eq 'text'){
			%plot_demo=(	
			
			type    => $type,
			file    => $file,
			# color list
			color   => $color,
			r1      => $r1,
			r0      => $r0,
			label_font => "light",
			# padding  - text margin in angular direction
			# rpadding - text margin in radial direction
			
			rpadding   => "5p",
			
			# Short lines can be placed before the label to connect them to the
			# label's position. This is most useful when the labels are
			# rearranged.
			
			show_links     => "no",
			link_dims      => "0p,2p,5p,2p,2p",
			link_thickness => "2p",
			link_color     => "black",
			);
		
	}
#-----------------------scatter plot -------------------------------
# show - as with highlights or links, determines whether to draw the plot or not
# type - determines the type of plot and can be one of scatter, line, histogram, heatmap, etc.
# file - file that contains plot data
 # min/max - the range of the plot axis, data outside this range are clipped
   # r0/r1 - the inner and outer radius of the plot track, which may be formatted in absolute or relative (or mix) form, just like for highlights
   #glyph - shape of glyph to use for the scatter plot and can be one of circle, rectangle, or triangle
   #glyph_size - size of the glyph, in pixels
   #color - if used, the the glyph will be solid and of this color
   #stroke_color - if used, the glyph will have an outline of this color
   #stroke_thickness - if used, the glyph's outline will have this thickness, in pixels 

#show  = yes
#type  = scatter

#file  = data/6/snp.density.txt
#r1    = 0.75r
#r0    = 0.90r
#max   = 1.0
#min   = 0.0

#glyph            = rectangle
#glyph_size       = 8
#color            = red
#stroke_color     = dred
#stroke_thickness = 1  
    return %plot_demo;
}

sub get_highlight{
	my ($file,$r0,$r1,$stroke_thickness,$stroke_color,$fill_color)=@_;
	$file=abs_path($file);
	my %highlight=(
		file       => $file, 
		r0         => $r0, 
		r1         => $r1, 
		stroke_thickness => $stroke_thickness ,
		stroke_color     => $stroke_color ,
		fill_color => $fill_color,
	);
	return %highlight;
}
sub get_link{
	my ($file,$r0,$fill_color)=@_;
	$file=abs_path($file);
	state $a = 10;
	$a=$a+10;
	my%link=(
		color=>$fill_color,
		z      => $a,
		file=>$file,
		radius => $r0,
		crest  => 0.5,
		bezier_radius        => "0r",
		bezier_radius_purity => 0.75,
		thickness     => 2
	
	);
	
	return %link;
	
	
}
sub get_circle_position{
	
	my($gap_between_circle,$thickness_of_circle)=@_;
	print Dumper($thickness_of_circle);
	my @circle_position;
	my $s= 1-$gap_between_circle;
	print $s,"\n";
	for (my $i=0;$i<@$thickness_of_circle;$i++){
		print $i,"\n";
		if($type[$i] eq "link"){
			push @circle_position,$thickness_of_circle->[$i];
			push @circle_position,$thickness_of_circle->[$i];
			next;
		}
		my $e=$s-$$thickness_of_circle[$i];
		push @circle_position,$e;
		push @circle_position,$s;
		$s=$e-$gap_between_circle;			
	}	
	print Dumper(\@circle_position);
	@circle_position=&add_scale(@circle_position,'r');
	print "@circle_position\n";
	return @circle_position;
}
sub add_scale{
	my $scale=pop @_;
	my @res;
	for (@_){
		
		push @res,$_.$scale;
	}
	if(@res==1){
		return $res[0];
	}else{
		return @res;
	}
}

sub format_chr{
	my$file=shift @_;
	my$i=0;
	my$col_len;
	open INC,"$file"  or die "$!";
	open OCHR,">$od/chr.info" or die "$!";
	while (my $line = <INC> ){
		chomp $line;
		my @tmp=split(/\s+/,$line);
		$col_len=@tmp if $.==1;
		last if ($col_len==7);
		next if($line=~/^#/ or $line=~/^$/ or  $line=~/^\s*$/);
		$i++;
		$i=1 if ($i>24);
		print OCHR "chr\t-\t$tmp[0]\t$tmp[0]\t0\t$tmp[1]\tchr$i\n";
	
	}
	close INC;
	close OCHR;
	if($col_len==7){
		return "$file";
		`rm $od/chr.info`;
	}else{
		return "$od/chr.info";
	}
	
}




sub get_min_max(){
	my @array;
	my $file=shift;
	#my $index=shift;
	
	open (IN,"$file") or die "$!";
	while(<IN>){
		chomp;
		next if (/^#/);
		my @tmp=split(/\s+/);
		push @array,$tmp[3];
	}
	my ($min)=(sort {$a<=>$b} @array);
	my ($max)=(sort {$b<=>$a} @array);
	close IN;
	return ($min,$max);
}



sub add_include(){
my $include=<<END;
<colors>
<<include etc/colors.conf>>
<<include etc/brewer.conf>>
#<<include etc/colors_fonts_patterns.conf>>
#<<include colors.ucsc.conf>>
#<<include colors.hsv.conf>>
</colors>

<fonts>
<<include etc/fonts.conf>>
</fonts>

<image>
<<include etc/image.conf>>
</image>
<<include etc/housekeeping.conf>>
END
#return $include;
open OUT ,">>$od/$outfile.conf" or die "$!";
print OUT $include;
close OUT;
#`sed  -i 's/_[0-9]//g'  $od/$outfile.conf`;
print  "/share/work/biosoft/circos/circos-0.69/bin/circos -conf $od/$outfile.conf -outputdir $od -outputfile $outfile\n";
system ("/share/work/biosoft/circos/circos-0.69/bin/circos -conf $od/$outfile.conf -outputdir $od -outputfile $outfile");
}




