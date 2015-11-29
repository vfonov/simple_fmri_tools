#! /usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use POSIX qw(floor);

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $time_fft=dirname($0).'/time_fft';
my $threshold=10;
my $low_freq=0;
my $hi_freq=0;

GetOptions (    
	      "verbose" => \$verbose,
        "clobber" => \$clobber,
        "low=f"   => \$low_freq,
        "hi=f"    => \$hi_freq
); 

die "Usage: $me <input> <output> [--verbose --clobber --low <freq Hz> --hi <freq Hz>] \n" if $#ARGV<1;

my ($input,$output)=@ARGV;

check_file($output) unless $clobber;

my @dims=split(/\s+/,`mincinfo -dimnames $input`);
my $input_type=`mincinfo -vartype image $input`;
chomp($input_type);
my $time_len=`mincinfo -dimlength time $input`+0;
die "This programm will work only on 4D minc file, one dimension should be time!\n" if $#dims<3 || $time_len<1;
#sanitity check
die "This programm expects regular spacing in time dimension!\n" if `mincinfo -attvalue time:spacing $input` !~ /^regular.*/;
my $time_step=`mincinfo -attvalue time:step $input`+0.0;
my ($hi,$low)=(0,0);

#converting from hertz to acquisition frequency
$hi= POSIX::floor($time_step*$time_len*$hi_freq);
$low=POSIX::floor($time_step*$time_len*$low_freq);

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

my $file=$input;
my $out="$tmpdir/output.mnc";
#print join(',',@dims),"\n";
if($dims[4]!~'time') #simplest case
{
  $file="$tmpdir/input.mnc";
  do_cmd('mincreshape','-dimorder','zspace,yspace,xspace,time',$input,$file);
}

do_cmd($time_fft,$file,$out,'--low',$low,'--hi',$hi);

do_cmd('mincreshape','-clobber','-'.$input_type,'-dimorder',join(',',@dims),$out,$output);
  
sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}
sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
