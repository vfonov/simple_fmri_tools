#! /usr/bin/env perl
#
use strict;
use warnings "all";
use File::Basename;
use File::Temp qw/ tempdir /;
use Getopt::Long;
use POSIX qw(floor);

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my ($hi,$lo,$slice);

GetOptions (
  "hi=f"      =>  \$hi,
  "lo=f"      =>  \$lo,
  "slice=n"   =>  \$slice,
  'verbose' => \$verbose,
  'clobber' => \$clobber
    );
    
die "Usage: $me <input.mnc> <sample.txt> <output.mnc> --hi <hz> --lo <hz> --slice <n> (slice number of sample) --verbose --clobber\n" if $#ARGV<2;
my ($in,$sample,$output)=@ARGV;


die "${output} exists!\n" if (-e $output) && !$clobber; 

my $tmpdir = tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

do_cmd('mincreshape',$in,"$tmpdir/in.mnc",'-dimorder','zspace,yspace,xspace,time');

my @args=('/data/ipl/proj01/nihpd/analysis_vladimir/data/quarantine/bin2/time_correlate',"$tmpdir/in.mnc",$sample,"$tmpdir/out.mnc");
push @args,'--hi',$hi if $hi;
push @args,'--low',$lo if $lo;
push @args,'--shift',$slice if $slice;

do_cmd(@args);
do_cmd('mincreshape',"$tmpdir/out.mnc",$output,'-dimorder','time,zspace,yspace,xspace','-clobber');

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake){
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}