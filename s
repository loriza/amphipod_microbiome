#!/usr/bin/env perl

use strict;
use warnings;
use utf8;
use File::Basename;
use Getopt::Std;
use Cwd qw( abs_path );
use File::Path qw( mkpath );
use FindBin qw($RealBin);

#===============================================================================
our ( $opt_i, $opt_l, $opt_o, $opt_G, $opt_D, $opt_h );

#===============================================================================
&usage() if ( 0 == @ARGV );
&usage() unless ( getopts('i:l:o:GDh') );
&usage() if ( defined $opt_h );
unless ( 0 == @ARGV ) {
    &usage("with undefined options: @ARGV");
}
&usage("lack input with: -i")  unless ( defined $opt_i );
&usage("lack input with: -l")  unless ( defined $opt_l );
&usage("lack output with: -o") unless ( defined $opt_o );

#===============================================================================
my ( $input, $genomelength, $prefix, $Catalog, $debug );
$input        = $opt_i;
$genomelength = $opt_l;
$prefix       = $opt_o;
$Catalog      = ( defined $opt_G ) ? 1 : 0;
$debug        = ( defined $opt_D ) ? 1 : 0;

#===============================================================================
my ( %gomlen, %genome, %gomReads );
$prefix =~ s/.sam$/.abundance/ if ( abs_path($prefix) eq abs_path($input) );
if ($Catalog) {
    &readcataloglength($genomelength);
    &loadsamfileDiamond( $input, $prefix );
}
else {
    &readgenomelength($genomelength);
    &loadsamfile( $input, $prefix );
}

#===============================================================================
sub loadsamfile {
    my ( $samfile, $prefix ) = @_;
    my ( @info, $IN, $sumReads );
    &openfile( $samfile, \$IN );
    while (<$IN>) {
        next if (/^@/);
        @info = split /\t/;
        unless ( $info[1] & 0x4 or $info[2] eq '*' ) {
            ++$sumReads;
            ++$gomReads{ $genome{ $info[2] } };
        }
        last;
    }
    while (<$IN>) {
        chomp;
        @info = split /\t/;
        unless ( $info[1] & 0x4 or $info[2] eq '*' ) {
            ++$sumReads;
            ++$gomReads{ $genome{ $info[2] } };
        }
    }
    close $IN;
    &printAbundance( $sumReads, $prefix );
}

#===============================================================================
sub loadsamfileDiamond {
    my ( $samfile, $prefix ) = @_;
    my ( @info, $IN, $reads, $score, $gene, $max, $sumReads );
    &openfile( $samfile, \$IN );
    while (<$IN>) {
        next if (/^@/);
        @info = split /\t/;
        next if ( $info[1] & 0x4 or $info[2] eq '*' );
        $reads = $info[0];
        $gene  = $info[2];
        $max   = ( split /:/, $info[11] )[2];
        last;
    }
    while (<$IN>) {
        chomp;
        @info = split /\t/;
        next if ( $info[1] & 0x4 or $info[2] eq '*' );
        $score = ( split /:/, $info[11] )[2];
        if ( $reads ne $info[0] ) {
            ++$sumReads;
            ++$gomReads{$gene};
            $reads = $info[0];
            $gene  = $info[2];
            $max   = $score;
        }
        elsif ( $score > $max ) {
            $max  = $score;
            $gene = $info[2];
        }
    }
    if ( $info[0] eq $reads ) {
        ++$sumReads;
        ++$gomReads{$gene};
    }
    close $IN;
    &printAbundanceCatalog( $sumReads, $prefix );
}

#===============================================================================
sub printAbundance {
    my ( $sumReads, $prefix ) = @_;
    my ($i);
    open OT, ">$prefix" or die "write $prefix $!\n";
    print OT "#Genome\tLength\tReads\tAbundance\n";
    print OT "#Total Aligned Reads\t$sumReads\n";
    foreach $i ( keys %gomlen ) {
        print OT $i, "\t", $gomlen{$i};
        print OT ( defined $gomReads{$i} )
          ? sprintf( "\t%d\t%g\n", $gomReads{$i}, $gomReads{$i} / $sumReads )
          : "\t0\t0\n";
    }
    close OT;
}

#===============================================================================
sub printAbundanceCatalog {
    my ( $sumReads, $prefix ) = @_;
    my ($i);
    open OT, ">$prefix" or die "write $prefix $!\n";
    print OT "#GeneID\tGeneLength\tReads\tAbundance\n";
    print OT "#Total Aligned Reads\t$sumReads\n";
    foreach $i ( keys %gomlen ) {
        print OT $i, "\t", $gomlen{$i};
        print OT ( defined $gomReads{$i} )
          ? sprintf( "\t%d\t%g\n", $gomReads{$i}, $gomReads{$i} / $sumReads )
          : "\t0\t0\n";
    }
    close OT;
}

#===============================================================================
sub readgenomelength {
    my ($genomelength) = @_;
    my ( @info, $IN );
    &openfile( $genomelength, \$IN );
    while (<$IN>) {
        chomp;
        next if (/^#/);
        @info = split /\t/;
        last if ( scalar @info < 3 );
        $genome{ $info[0] } = $info[2];
        unless ( defined $gomlen{ $info[2] } ) {
            $gomlen{ $info[2] } = $info[1];
        }
        else {
            $gomlen{ $info[2] } += $info[1];
        }
    }
    close $IN;
}

#===============================================================================
sub readcataloglength {
    my ($cataloglength) = @_;
    my ( @info, $IN );
    &openfile( $cataloglength, \$IN );
    while (<$IN>) {
        chomp;
        next if (/^#/);
        @info = split /\t/;
        $gomlen{ $info[0] } = $info[1];
    }
    close $IN;
}

#===============================================================================
sub openfile {
    my ( $infile, $IN ) = @_;
    if ( $infile =~ /\.gz$/ ) {
        open $$IN, "gzip -dc $infile |" or die "gzip -dc $infile $!\n";
    }
    else {
        open $$IN, "<$infile" or die "read $infile $!\n";
    }
}

#===============================================================================
sub systemCommand {
    my ( $cmd, $cmdpath ) = @_;

    unless ($cmdpath) {
        open FD,
          "find -L $RealBin -name $cmd -type f -print -quit 2>/dev/null |"
          or die "find -L $RealBin -name $cmd -type f -print -quit $!\n";
        $cmdpath = <FD>;
        close FD;
        chomp($cmdpath) if ($cmdpath);
    }

    unless ($cmdpath) {
        open SH, "which $cmd 2>/dev/null|" or die "which $cmd $!\n";
        while (<SH>) {
            chomp;
            $cmdpath = $_;
        }
        close SH;
    }

    &usage("program:\t$cmd not found with current bin/\$PATH")
      unless ($cmdpath);
    &usage("program:\t$cmd not executable of $cmdpath\n")
      unless ( -x $cmdpath );

    return $cmdpath;
}

#===============================================================================
sub localScript {
    my ( $src, $srcpath ) = @_;

    unless ($srcpath) {
        open FD,
          "find -L $RealBin -name $src -type f -print -quit 2>/dev/null |"
          or die "find -L $RealBin -name $src -type f -print -quit $!\n";
        $srcpath = <FD>;
        close FD;
        chomp($srcpath) if ($srcpath);
    }

    &usage("file:\t$src not found with current bin/")
      unless ($srcpath);
    &usage("file:\t$src not readable of $srcpath\n")
      unless ( -e $srcpath );

    if ( $srcpath =~ /\.pl$/ ) {
        return "perl $srcpath";
    }
    elsif ( $srcpath =~ /\.py$/ ) {
        return "python $srcpath";
    }
    else {
        return $srcpath;
    }
}

#===============================================================================
sub subsystem {
    my ( $cmd, $log, $war ) = @_;
    my ($time);
    if ($debug) {
        $log = ( defined $log ) ? "1>$log" : "";
        $war = ( defined $war ) ? "2>$war" : "";
        print STDOUT "start by: ", $time = localtime, "\n";
        print STDOUT "$cmd $log $war\n";
    }
    else {
        $log = ( defined $log ) ? "1>$log" : "1>/dev/null";
        $war = ( defined $war ) ? "2>$war" : "2>/dev/null";
    }
    if ( system("$cmd $log $war") ) {
        print STDERR "sh $cmd $!\n";
        exit;
    }
    print STDOUT "finish by: ", $time = localtime, "\n" if ($debug);
}

#===============================================================================
sub versionByDate {
    my ($time) = @_;
    my (@info);
    @info = split /\s+/, $time;
    $info[3] =~ s/^(\d+):(\d+):\d+/$1.$2/;
    return "$info[-1].$info[1]$info[2]";
}

#===============================================================================
sub submkdir {
    my ($subdir) = @_;
    mkpath( $subdir, 0, 0751 )
      or die "mkdir -p $subdir $!\n"
      unless ( -d $subdir );
}

#===============================================================================
sub usage {
    my ($reason) = @_;

    if ( defined $reason ) {
        print STDERR <<STE;

================================================================================
    $reason
================================================================================
STE
    }

    print STDERR <<STE;

    Last modify: $Version
    Contact: $Contact

    Usage:
    \$$Script [options]
    -i  <str>   input samfile for abundance, e.g. required
    -l  <str>   input index length for abundance, required
    -o  <str>   output file for result, e.g. required
    -G  gene    samfile is gene catalog, default genome
    -D  debug   do show information by stdout, default not show
    -h  help
STE
    exit;
}
