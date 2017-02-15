#!/usr/bin/perl -w
#print " 0 - mips sources\n 1 - AGN\n 2 - A901\n 3 - P. Barmby's X-ray sources\n 4 - 2006 data \n 5 - 2007 data\n 6 - AGN 2008\n";
#
# Originally written by CNAW.  Minor changes by BJW.
#
# Path where ds9 regions files must be written
# **EDIT**
$regions_path = '/home/bjw/field2mips/data2008b/';
#
# Path to raw data directory. Input catalogues to plugcat will be 
# written there
# BJW - also going to use this to find the obsdir data subdirectories
# **EDIT**
$raw_data_path ='/home/bjw/field2mips/data2008b/';
#
#
#
$agn = 0 ;
$setup = 0 ;
#
# CNAW Mips 2007 sources
#
#if($#ARGV == -1) {
#    $agn = 0 ;
#    $source='mips_2008';
#}
##
## AGN
##
#if($#ARGV > -1 ) {
#    $field =  $ARGV[0] ;
#    $field =~ s/\s//g  ;
#    if($#ARGV == 1) {
#	$setup =   $ARGV[1];
#	$setup =~ s/\s//g  ;
#	$source='agn_2007';
#    }
#    if($field eq '1' || $field eq 'agn' || 
#       $field eq 'AGN' || $field eq 'Agn') {
#	$agn = 1 ;}
##
## A901 
##
#    if($field eq '2' ||  $field eq 'A901' ||  $field eq 'a901') {
#	$agn = 2 ;
#	$setup = $ARGV[1] ;
#    }
## 
## P. Barmby's X-ray sources
##
#    if($field eq '3') {
#	$agn = 3 ;
#    }
##
## 2006 EGS 
##
#    if($field eq '4') {
#	$agn = 4 ;
#    }
##
## 2007 EGS
##
#    if($field eq '5') {
#	$agn = 5 ;
#	$source='mips_2007';
#    }
##
## 2008 AGN
##
#    if($field eq '6') {
#	$agn = 6 ;
#	$source='agn_2008';
#    }
#}
#print "agn is $agn\n";
#
# Generate DS9 regions files for objects configured by the xfitfibs 
# program for hectospec observations.
#  **EDIT**
$pi = 4.0* atan2(1.0, 1.0) ;
$iraf_geomap    = '/home/bjw/field2mips/data2008b/geomap.csh';
#
# Set up these variables to make statistics later
#
foreach ($k = 0 ; $k <= 10 ; $k++) {
    $priority[$k] = 0 ;
    $configured[$k] = 0 ;
}
# 
# Original catalog provided to xfitfibs, valid for all fields
#
#if($agn == 0) {
#    $cat_root_name = 'hectospec_2008_egs';
#    $obsdir[0] = '2008.0310';
#    $config[0] = 1;
#    $obsdir[1] = '2008.0315';
#    $config[1] = 2;
#    $cat = '/home/cnaw/xfitfibs/'.$cat_root_name.'.cat';
#}
#if($agn == 5) {
#    $cat_root_name = 'hectospec_2007_egs';
#    $obsdir[0] = '2007.0508';
#    $obsdir[1] = '2007.0510';
#    $obsdir[2] = '2007.0612';
#    $obsdir[3] = '2007.0614';
#    $obsdir[4] = '2007.0615';
#    $config[0] = 1;
#    $config[1] = 2;
#    $config[2] = 3;
#    $config[3] = 4;
#    $config[4] = 5;
#    $cat = '/home/cnaw/xfitfibs/'.$cat_root_name.'.cat';
#}
#if($agn == 1) {
#    $cat_root_name = 'agn_2007';
#    $obsdir[0] = '2007.0509';
#    $obsdir[1] = '2007.0514';
#    $obsdir[2] = '2007.0717';
#    $config[0] = 1;
#    $config[1] = 2;
#    $config[2] = 5;
#    $cat = '/home/cnaw/xfitfibs/'.$cat_root_name.'.cat';
#} 
##
## A901 data were organised differently
##
#if($agn == 2) {
#    if($setup == 1) {
#	$obsdir[0] = '2006.0322';
#	$config[0] = 1;
#	$cat_root_name = 'a901_faint_xfitfibs';
#    }
#    if($setup == 2) {
#	$cat_root_name = 'a901_bright_xfitfibs';
#	$obsdir[0] = '2006.0402';
#	$config[0] = 2;
#    }
#    $cat = '/data/cnaw/mmt/'.$obsdir[0].'/'.$cat_root_name.'.cat';
#} 
##
## Pauline Barmby's setups
##
#if($agn == 3) {
#    $cat_root_name = 'master4_sb';
#    $obsdir[0] = '2007.0314';
#    $config[0] = 1;
#    $obsdir[1] = '2007.0316';
#    $config[1] = 2;
#    $obsdir[2] = '2007.0419';
#    $config[2] = 3;
#    $cat = '/data/cnaw/mmt/pbarmby/obsplan/master4_sb.cat';
#}
##
## 2006 setups
##
#if($agn == 4) {
#    $cat_root_name = 'new_egs_hecto_spec';
#    $obsdir[0] = '2006.0502';
#    $config[0] = 2;
#    $obsdir[1] = '2006.0504';
#    $config[1] = 15;
#    $obsdir[2] = '2006.0505';
#    $config[2] = 16;
#    $cat = '/home/cnaw/mmt/newest/new_egs_hecto_spec.cat';
##
## But we don't want to deal with this:
##
#    die ;
#}
##
## 2008 AGN
##
#if($agn == 6) {
#    $cat_root_name = 'agn_2008';
#    $obsdir[0] = '2008.0311';
#    $config[0] = 1;
#    $obsdir[1] = '2008.0314';
#    $config[1] = 2;
#    $cat = '/home/cnaw/xfitfibs/'.$cat_root_name.'.cat';
#}
#
#  BJW - just set the cat, date and config w/o trying to switch on the dataset
#  This needs to find the the input <something>.cat file that you
#  fed to xfitfibs to generate the Hectospec configs.
#  **EDIT**
#$cat_root_name = 'finalhecto2.all.edit.hms';
$cat_root_name = 'field2_2008b';
$obsdir[0] = '2008.0505';
$config[0] = 1;
$obsdir[1] = '2008.0531';
$config[1] = 2;
$obsdir[2] = '2008.0602';
$config[2] = 3;
$obsdir[3] = '2008.0603';
$config[3] = 4;
$obsdir[4] = '2008.0604';
$config[4] = 5;
$cat = '/home/bjw/field2mips/hectospec3/'.$cat_root_name.'.cat';


#
# Read map file, which tells us where each fiber was effectively placed.
# Note that the sky (rank=-1) positions are bogus.
#
# 
# Find the map file
#
foreach ($k = 0 ; $k <= $#config ; $k++) {
    $dir =$raw_data_path.$obsdir[$k].'/'.$cat_root_name.'_'.$config[$k].'*_map';
    @list = `ls $dir`;
    @clist = `ls $dir`;
    print "First map file for this configuration is\n$list[0]\n";
    if($#list >= 0 ) {
	chop $list[0] ;
	$map[$k] = $list[0] ;
	$cat[$k] = $clist[0];
    } else {
	print "No map files !\n";
	die;
    }
#
# Since 2007 there is only one map file. Copy each half to two subfiles so
# hsred works
#
    @mapname = split('\/',$map[$k]) ;
    $map1 = $map[$k].'_1';
    $map2 = $map[$k].'_2';
    open(MAP,"<$map[$k]") || die "cannot open $map[$k]";
    open(MAP1,">$map1")   || die "cannot open $map1";
    open(MAP2,">$map2")   || die "cannot open $map2";
#
# Copy header
#
    $line = <MAP> ;
    print MAP1 $line;
    print MAP2 $line;
    $line = <MAP> ;
    print MAP1 $line;
    print MAP2 $line;
#
# read lines
#
    $j = 0 ;
    while(<MAP>) {
	$j++ ;
	if($j <= 150) {
	    print MAP1 $_;
	} else {
	    print MAP2 $_;
	}
    }
    close (MAP) ;
    close (MAP1) ;
    close (MAP2) ;
#
# Create plugdir if it does not exist
    $plugdir = $raw_data_path.$obsdir[$k].'/plugdir' ;
    if ( ! -e $plugdir) {
	$command  ='mkdir '.$plugdir;
	system($command);
    } else {
	print "file $plugdir exists\n";
    }
#
# Write idl reduction script. This assumes that the "uberextraction" will
# be made (i.e., response correction, remove telluric absorption)
#
    $script = $raw_data_path.$obsdir[$k].'/reduce.idl';
    open(SCRIPT, ">$script") || die "cannot open $script";

    print SCRIPT "retall\n.r /home/bjw/idl/hsred/new_hs_readcat\n";
    print SCRIPT "hs_preproc\nhs_calibproc, /doall\n";
    $plugcat = "cat='./plugdir/".$cat_root_name.'_'.$config[$k]."_cat'";
    print SCRIPT $plugcat,"\n";
    print SCRIPT "cols = ['beam','ra','dec','rmag','fiber','rank','id']\n";
    print SCRIPT "format= ['L','D','D','D','L','L','A']\n";
    $mapfile =  "mapfile='".$mapname[$#mapname]."'";
    print SCRIPT $mapfile,"\n";
    print SCRIPT "hs_readcat,cat,mapfile,cols,format\n";
    $line = "files = findfile('". $cat_root_name."*.fits')" ;
    print SCRIPT $line,"\n";
    $plugdir   = "plugfile='./plugdir/".$cat_root_name.'_'.$config[$k]."'";
#    $line   = "hs_extract, files,".$plugdir.",/plugcat, /docosmic,\$";
    $line   = "hs_extract, files,".$plugdir.",/plugcat,\$";
    print SCRIPT $line,"\n";
    $line = " /uberextract, outname='".$obsdir[$k].".fits'";
    print SCRIPT $line,"\n";
    $line = "hs_reduce1d,'".$raw_data_path.$obsdir[$k]."/reduction/0000/spHect-".$obsdir[$k].".fits'";
    print SCRIPT $line,"\n";
    close(SCRIPT);
}
#
# Read source catalogue since the CFG file uses the _order_ in 
# this catalogue rather than object ID to identify a source
#
print "$cat\n";
open(CAT,"<$cat") || die "cannot open $cat";
<CAT>;
<CAT>;
$n=0 ;
while(<CAT>) {
    chop $_ ;
    $n++ ;
    $id = sprintf("%05d",$n) ;
    
    $catalog{$id} = $_ ;
    @junk = split('\t',$_);
    $rank = $junk[3];
# Count number of objects in each rank level ;
    if($rank ne '    ' && $rank ne ' ' && $rank ne '   ' ) {
	$priority[$rank] = $priority[$rank] + 1;
    }
}
close(CAT);
#
# Find how many objects are in each priority class. Note that this does not 
# identify objects that appear multiple times in the list
#
foreach ($k = 0 ; $k <= 10 ; $k++) {
    print "$k, $priority[$k]\n";
    $nset[$k] = 0 ;
}
#
# Read the map files, separate objects and sky positions, write
# input files to geomap to calculate sky positions, output object to
# chartable file
#
foreach ($obs=0 ; $obs <= $#map ; $obs++) {
    open(CAT,"<$map[$obs]") || die "cannot open $map[$obs]";
    print "<CAT> is $map[$obs]\n";
    <CAT> ;
    <CAT> ;
#
# File output in chart format
#
#    if ($agn == 0) {
	$name = sprintf("config_%02d",$config[$obs]);
	$name = sprintf("config_%02d",$config[$obs]);
	$chart  =$name.'.chart';
#    } 
#    if ($agn == 1) {
#	$name = sprintf("agn_%02d",$config[$obs]);
#	$name = sprintf("agn_%02d",$config[$obs]);
#	$chart  =$name.'.chart';
#    }
#    if ($agn == 2) {
#	$name = sprintf("a901_%02d",$config[$obs]);
#	$name = sprintf("a901_%02d",$config[$obs]);
#	$chart  =$name.'.chart';
#    }
#    if( $agn == 3) {
#	$name = sprintf("pbarmby_%02d",$config[$obs]);
#	$name = sprintf("pbarmby_%02d",$config[$obs]);
#	$chart  =$name.'.chart';
#    }
#    if ($agn == 6) {
#	$name = sprintf("agn_2008_%02d",$config[$obs]);
#	$name = sprintf("agn_2008_%02d",$config[$obs]);
#	$chart  =$name.'.chart';
#    }
    open(CHART,">$chart") || die "cannot open $chart";
#
# Output file in ds9 regions format
#
    $ds9  =$regions_path.$name.'.reg';
    open(DS9,">$ds9") || die "cannot open $ds9";
#
# Output file to calculate Ra, Dec for random sky positions placed
# by xfitfibs
#
    $geomap = 'geomap_'.$name.'.dat';
    open(GEO,">$geomap") || die "cannot open $geomap";
#
# Output file to containing the platex, platey for sky positions placed
# by xfitfibs
#
    $coords = 'coords_'.$name.'.dat';
    open(COORDS,">$coords") || die "cannot open $coords";
    foreach ($k = 0 ; $k <= 10 ; $k++) {
	$nset[$k] = 0 ;
    }
#
    $ns     = -1 ;
    $k      = -1 ;
    while(<CAT>) {
	$line = $_;
	chop $line ;
#	print "$line\n";
	($aperture, $beam, $id, $rahms, $decdms, $target, $fiber, 
	 $platex, $platey) 
	    = split(' ',$line);
	($rah, $ram,$ras) = split(':',$rahms) ;
	$ra = $rah + $ram/60.+$ras/3600.;
	$ra = sprintf("%12.8f",$ra * 15.0);
	($rah, $ram,$ras) = split(':',$decdms) ;
	$dec = sprintf("%12.8f",$rah + $ram/60.+$ras/3600.);
	$rmag = 0.0 ;
	$color = 'color=black' ;
	$symbol= '# point=cross';
	$k++ ;
# Sky position
	if($target == -1 ) {
	    $color = 'color=blue' ;
	    $symbol= '# point=cross';
	    $minus1++;
	    print COORDS join(' ',$platex, $platey),"\n";
	    $ns++;
	    $sky[$ns] = join(' ',$aperture, $platex, $platey, 
			     $rahms, $decdms,$id, $target);
	    $nset[$#nset] = $nset[$#nset]+1 ;
	    $rank  = -1 ;
	    $id    ='sky';
	    $rmag   = 0.0;
	}
# parked fibers:
	if($target == 0  ) {
	    $color = 'color=black' ;
	    $symbol= '# point=cross';
	    print COORDS join(' ',$platex, $platey),"\n";
	    $ns++;
	    $sky[$ns] = join(' ',$aperture, $platex, $platey, 
			     $rahms, $decdms,$id, $target);
	    $nset[0] = $nset[0] + 1;
	    $sample  ='parked';
	    $id    ='parked';
	    $rmag   = 0.0;
	}

# Find rank for galaxy
	if($target > 0) {
	    $galaxy = sprintf("%05d",$target);
	    if( ! exists($catalog{$galaxy})) { 
		print "\n\n$galaxy, skipping $catalog{$galaxy}\n";
		next ;
	    }
# This is where we read the input catalog for ID, rank, sample, and mag.
# If your input catalog differs in format (extra fields or spaces) 
# you may have to change some of this code.
	    @data = split(' ',$catalog{$galaxy});
	    $id   = $data[2] ;
# For dealing with a catalog where the guide stars had no rank printed
#	    $ncols = @data ;
	    $foo = $data[3] ;
#	    if($ncols = 6) {
	    if($foo eq "   ") {
		$rank   = 0;
	    } else {
		$rank   = $data[3];
	    }
	    $sample = $data[4];
	    $rmag   = $data[6];
	    if($rank eq " ") {
		    print "$line\n$catalog{$galaxy}\n";
		next;
	    }
#
# This will allow deriving precise coordinate for sky fibers
# (only necessary for 2006 data since 2007 map files have good 
#  sky coordinates)
#
	    print GEO join(' ',$platex, $platey, $rahms,$decdms),"\n";
#
# Set colours for objects according to rank
# Fstars (for EGS data):
	    if($rank == 1) {
		$color = 'color=yellow' ;
		$symbol= '# point=boxcircle';
	    }
# High priority galaxies
	    if($rank == 2) {
		$symbol= '# point=diamond';	
		$color = 'color=red' ;
	    }
	    if($rank == 3) {
		$color = 'color=yellow' ;
		$symbol= '# point=diamond';
	    }
	    if($rank == 4) {
		$color = 'color=red' ;
		$symbol= '# point=circle';
	    }
	    if($rank == 5) {
		$color = 'color=white' ;
		$symbol= '# point=circle';
	    }
	    if($rank == 6) {
		$symbol= '# point=box';
		$color = 'color=black' ;
	    }
	    $nset[$rank] = $nset[$rank] + 1;
	    $configured[$rank] = $configured[$rank]+1;

	    $id =~ s/\s//g;
#
# Output stars and galaxies to ds9 regions file
#	
	    $fiber_ds9   = join('_',$aperture,$id);
	    $line=join(' ','fk5;point(',$ra, $dec,$symbol,
		       $color,'text={',$fiber_ds9,'}');
	    print DS9 $line,"\n";
	    $chart = sprintf("%3d %3d %12s %12s %8.3f %8.3f %5d %5.1f %10.5f %10.5f %5.2f %s",
			     $aperture, $fiber, $rahms, $decdms, $platex,$platey,$target, $rank, $ra/15.,$dec, $rmag, $id) ;
#	    print "$chart\n";
	    print CHART  $chart,"\n";
	}
#
# Store data that will be written in the plugcat
#
	$plug_data[$k] = join(' ', $aperture, $beam, $ra/15., $dec,
			      $rmag, $rank,$fiber, $platex, $platey,$id) ;
#	print "$id $plug_data[$k]\n";
    }
#
# The catalogue used in the reduction should not contain sky positions
# which are placed in another catalogue
#
    close(GEO) ;
    close(COORDS) ;
    $symbol = '# point=boxcircle';
    $color  =  $color = 'color=cyan' ;
#
# Run the geomap IRAF scripts
#
    $skies = $name.'.sky';
    $log   = 'log';
    system("\\rm hectospec.map $skies results $log") ;
    open( IRAF, "| $iraf_geomap > $log") ;
    print IRAF <<END;
geomap
$geomap
hectospec.map
-320.
320.
-320.
320.
config
results
general
legendre
6
6
half
6
6
half
1
3
real
no
no
stdgraph
bye
END
;
      close (IRAF) ;
      open( IRAF, "| $iraf_geomap >> $log") ;
      print IRAF <<END;
geoxytran
$coords
$skies
hectospec.map
config
geometric
forward
real
1
2
%12.8f
%12.8f
12
INDEF
INDEF
INDEF
INDEF
INDEF
INDEF
INDEF
INDEF
INDEF
INDEF
bye
END
 ;
    close (IRAF);
#    
# This file contains the transformed RA, Dec of sky positions.
#
    open(SKY,"<$skies") || die "cannot open $skies";
    $n = -1 ;
    while(<SKY>) {
	chop $_ ;
	($ra, $dec) = split(' ',$_) ;
	$rahms  = &Dec2Sex($ra,12,3) ;
	$decdms = &Dec2Sex($dec,12,3) ;
	$rahms  =~ s/\s//g;
	$decdms =~ s/\s//g;
	$ra = $ra * 15. ;
	$n++ ;
	($aperture,$platex,$platey, $raorg,$decorg,$id, $target)  = split(' ',$sky[$n]);
#
# Recover original data 
#
	$k = $aperture -1 ;
	($iap, $beam, $bogus_ra, $bogus_dec, $rmag, $rank, $fiber,$id) =
	    split(' ',$plug_data[$k] ) ;
	if($target == -1) {
	    $id ='sky_'.sprintf("%03d",$k);
	    $plug_data[$k] = join(' ', $aperture, $beam, $ra/15., $dec,
				  $rmag, $rank,$fiber, $platex, $platey,$id) ;
	    $sky_id = join('_',$config[$obs],$aperture);
	}
	if($target ==  0) {
	    $sky_id ='parked_'.sprintf("%03d",$aperture);
	}

	$line=join(' ','fk5;point(',$ra, $dec,$symbol,$color,'text={',$sky_id,'}');
#
# Write sky (and parked ?) to ds9 regions file
#
	print DS9 $line,"\n";
	$chart = sprintf("%3d %3d %12s %12s %8.3f %8.3f %5d %5.1f %10.5f %10.5f %5.2f %s",
			 $aperture, $fiber, $rahms, $decdms, $platex,$platey,$target, $rank, $ra/15.,$dec, $rmag, $sky_id) ;
	print CHART  $chart,"\n";
	
    }
    close(SKY) ;
    close(DS9) ;
    close(CHART) ;
    
# Names for ranks?
#    @names = ('parked','Spec-ph','high','scuba+irac','deep2',
#	      'outriggers', 'sdss','sky','-1');
    @names = ('parked','Spec-ph','70um','24um_hiz','24um_lowz',
	      'other5', 'other6','sky','-1');
    $line = join(' ',@names);
#    print "$line\n";

    $line = join(' ',@nset);
#    print "$line\n";
#
# Output plugcat files, creating  the 5 required files
#
    $plugall  =$raw_data_path.$obsdir[$obs].'/plugdir/'.$cat_root_name.'_'.$config[$obs].'_all';
    open(PLUG_ALL,">$plugall") || die "cannot open $plugall";

    $plugtarg  =$raw_data_path.$obsdir[$obs].'/plugdir/'.$cat_root_name.'_'.$config[$obs].'_targets';
    open(TARGETS,">$plugtarg") || die "cannot open $plugtarg";

    $plugsky  =$raw_data_path.$obsdir[$obs].'/plugdir/'.$cat_root_name.'_'.$config[$obs].'_sky';
    open(PLUG_SKY,">$plugsky") || die "cannot open $plugsky";

    $plugu  =$raw_data_path.$obsdir[$obs].'/plugdir/'.$cat_root_name.'_'.$config[$obs].'_unused';
    open(PLUG_UNUSED,">$plugu") || die "cannot open $plugu";

    $plugcat  =$raw_data_path.$obsdir[$obs].'/plugdir/'.$cat_root_name.'_'.$config[$obs].'_cat';
    open(PLUGCAT,">$plugcat") || die "cannot open $plugcat";
#
# Loop through objects in this setup and write the catalogues
#

    foreach ($k=0 ; $k <= $#plug_data ; $k++) {
	($aperture, $beam, $ra, $dec, $rmag, $rank, $fiber, $plx, $ply,$id) 
	    = split(' ',$plug_data[$k] ) ;
#	print "$k $id, $plug_data[$k]\n";
	$dmin = 0.0 ;
	$junk = 0   ;
	$line = sprintf("%5d %7.3f %12.6f %12.6f %7.3f %7.3f %4d %4d %4d %4d",
			$aperture, $dmin, $ra, $dec, $rmag, $rmag, $junk, $junk, 
			$junk, $junk);
	print PLUG_ALL $line,"\n";
	if($beam == 0 ) {
	    $line = join(' ',$aperture,'SKY') ;
	    print PLUG_SKY $line,"\n";
	}
	if($beam == 1 ) {
	    print TARGETS $line,"\n";
	}
	if($beam == 2 ) {
	    $line = join(' ',$aperture,'unused') ;
	    print PLUG_UNUSED $line,"\n";
	}
	$line =  sprintf("%5d %12.6f %12.6f %7.3f %3d %2d %s",
			$aperture, $ra, $dec, $rmag, $fiber,$beam, $id);
	print PLUGCAT $line,"\n";
    }
    close(PLUG_ALL)   ;
    close(TARGETS)   ;
    close(PLUG_SKY)   ;
    close(PLUG_UNUSED)   ;
    close(PLUGCAT)   ;
    close(CAT);
}
#
# list total number of objects observed in a given rank
#
foreach ($k = 0 ; $k <= 10 ; $k++) {
    if($priority[$k] == 0) { next;}
    $percent = $configured[$k]/$priority[$k] ;
    $line = sprintf("%2d %5d %5d %5.2f",
		    $k, $priority[$k], $configured[$k], $percent);
    print "$line\n";
}
#
sub Dec2Sex{
#-----------------------------------------------------------------------
# dec2sex / GDW / 9 Mar 93
#
# Converts from decimal to sexagesimal (e.g. +HH:MM:SS.SS)
# Arguments: value, string width, number of decimal places
#-----------------------------------------------------------------------
    my( $a, $w, $d) ;
    my( $pm);
    my( $str);
    $d = 2;
    ( $a, $w, $d) = @_;
    if( $a<0){ $pm="-" } else { $pm=" " }
    $a = abs($a);
    my( $dd) = int($a);
    my( $mm) = int(($a - $dd)*60);
    my( $ss) = ($a - $dd - $mm/60)*3600;
    $ss = sprintf("%.${d}f",$ss);
    if($ss < 10.00) {
#       $str = sprintf("%1s%02d:%02d:%0.${d}f",$pm,$dd,$mm,$ss);
        $str = sprintf("%1s%02d:%02d:%1s%.${d}f",$pm,$dd,$mm,'0',$ss);
    } else {
#       $str = sprintf("%1s%02d:%02d:%0.${d}f",$pm,$dd,$mm,$ss);
        $str = sprintf("%1s%02d:%02d:%.${d}f",$pm,$dd,$mm,$ss);
    }
    return sprintf("%${w}s",$str);
}

#-----------------------------------------------------------------------------
# Add offsets between object coordinates and image coordinates.
# This is done because the currently reconstructed DEEP2 images
# have an WCS offset.
#
sub Deep2_coords {
    my( $ra, $dec) = @_ ;
    my( $pi) = 4.0* atan2(1.0, 1.0) ;
    my($ra_corr, $dec_corr)= (0.0, 0.0);
#    $ra_corr = $ra   + 0.18/3600. ;
#    $dec_corr = $dec + 3.16/(3600.*cos($dec*$pi/180.)) ;
#    $dec_corr = $dec - 2.42/3600. ;
#    $ra_corr =  $ra + 2.75/(3600.*cos($dec*$pi/180.));
    return $ra_corr, $dec_corr;
}
