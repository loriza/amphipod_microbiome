#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Getopt::Std;
use FindBin qw($RealBin);

#===============================================================================
our ( $opt_i, $opt_o, $opt_p, $opt_c, $opt_r, $opt_w, $opt_e, $opt_h );

#===============================================================================
&usage if ( 0 == @ARGV );
&usage unless ( getopts('i:o:p:c:r:w:eh') );
&usage if ( defined $opt_h );
unless ( 0 == @ARGV ) {
    &usage("with undefined options: @ARGV");
}
&usage("lack input profiling table with: -i") unless ( defined $opt_i );
&usage("lack input phenotype table with: -p") unless ( defined $opt_p );
&usage("lack output prefix with: -o")         unless ( defined $opt_o );

#===============================================================================
my ( $profile, $phenotype, $prefix, $column, $arrange, $width, $execute,
    $LEFSe );
$LEFSe     = $RealBin;
$profile   = $opt_i;
$phenotype = $opt_p;
$prefix = ( $opt_o ne $opt_i && $opt_o ne $opt_p ) ? $opt_o : $opt_o . ".LEFSe";
$column = ( defined $opt_c ) ? $opt_c : 0;
$arrange = ( defined $opt_r ) ? $opt_r : 0;
$width   = ( defined $opt_w ) ? $opt_w : 11;
$execute = ( defined $opt_e ) ? 1 : 0;

#===============================================================================
my $condapath = "/work/install/miniconda/bin";
my $condaenvs = "LEfSe";

#===============================================================================
my %nextlevel = (
    "k" => "p",
    "p" => "c",
    "c" => "o",
    "o" => "f",
    "f" => "g",
    "g" => "s"
);

#===============================================================================
my ( %subsample, %grouporder );
&open_phenotype( $phenotype, \%subsample, $column, $arrange );
&profile2LEFSe( $profile, \%subsample, $prefix );
&LEFSe2plot( $prefix, $execute );

#===============================================================================
sub LEFSe2plot {
    my ( $prefix, $execute ) = @_;
    open SH, ">$prefix.sh" or die "write $prefix.sh $!\n";
    print SH <<shell;
export PATH="$condapath:\$PATH"
source activate $condaenvs
python $LEFSe/format_input.py $prefix $prefix.format -c 1 -u 2 -o 1000000
python $LEFSe/run_lefse.py $prefix.format $prefix.res #-l 2.5
python $LEFSe/plot_cladogram.py $prefix.res $prefix.cladogram.png --format png \\
--max_lev 7 --labeled_start_lev 2 --labeled_stop_lev 7 --abrv_start_lev 3 \\
--abrv_stop_lev 7 --right_space_prop 0.3 --left_space_prop 0 \\
--title "log10 (LDA) >2" --width $width
python $LEFSe/plot_cladogram.py $prefix.res $prefix.cladogram.pdf --format pdf \\
--max_lev 7 --labeled_start_lev 2 --labeled_stop_lev 7 --abrv_start_lev 3 \\
--abrv_stop_lev 7 --right_space_prop 0.3 --left_space_prop 0 \\
--title "log10 (LDA) >2" --width $width
#$LEFSe/plot_res.py --format svg $prefix.res $prefix.svg
python $LEFSe/plot_res.py --format png $prefix.res $prefix.png
python $LEFSe/plot_res.py --format pdf $prefix.res $prefix.pdf
mkdir -p $prefix.biomarkers_images
python $LEFSe/plot_features.py $prefix.format $prefix.res $prefix.biomarkers_images/
shell
    close SH;
    system("sh $prefix.sh 1>$prefix.sh.log ") if ($execute);
}

#===============================================================================
sub profile2LEFSe {
    my ( $profile, $subsample, $prefix ) = @_;
    my (
        @info, @title, @column, $hcolumn, $rcolumn,
        $row,  $i,     @tax,    $species, $lastlevel
    );
    if ( $profile =~ /\.gz/ ) {
        open PR, "gzip -dc $profile |" or die "gzip -dc $profile $!\n";
    }
    else {
        open PR, "<$profile" or die "read $profile $!\n";
    }
    open OT, ">$prefix" or die "write $prefix $!\n";
    print OT "State";
    chomp( $_ = <PR> );
    chomp( $_ = <PR> ) if ( $_ =~ /^# Constructed from biom file/ );
    @title   = split /\t/;
    $hcolumn = scalar @title;
    &usage("table is not '\\t' as separate") if ( 1 eq $hcolumn );

    for ( $i = 1 ; $i < @title ; ++$i ) {
        push( @column, $i ) if ( defined $$subsample{ $title[$i] } );
    }
    &usage("unmatched sample name between profile and phenotype")
      if ( 2 > scalar @column );
    foreach $i (@column) {
        print OT "\t", $$subsample{ $title[$i] };
    }
    print OT "\nSample";
    foreach $i (@column) {
        print OT "\t", $title[$i];
    }
    while (<PR>) {
        chomp;
        ++$row;
        next if (/^#/);
        @info    = split /\t/;
        $rcolumn = scalar @info;
        &usage(
"data row $row do not match head's columns: $rcolumn\t$hcolumn\n$_\n"
        ) if ( $hcolumn ne $rcolumn );
        $info[0] =~ s/\s+//g;
        @tax = split /;/, $info[0];
        $info[0] = $tax[0];

        for ( $i = 1 ; $i < @tax ; ++$i ) {
            if ( $tax[$i] =~ /[kpcofg]__\S+/ ) {
                $info[0] .= "|" . $tax[$i];
            }
            elsif ( $tax[$i] =~ /s__(\S+)/ ) {
                $species = $tax[ $i - 1 ] . "_" . $1;
                $species =~ s/g__/s__/;
                $info[0] .= "|" . $species;
            }
            elsif ( $tax[$i] =~ /Other/ ) {
                $species = $tax[ $i - 1 ] . "_Other";
                $species =~ s/^([kpcofg])__/$nextlevel{$1}__/;
                $info[0] .= "|" . $species;
                last;
            }
            elsif ( $tax[$i] =~ /([kpcofgs])__$/ ) {
                $species = $tax[ $i - 1 ] . "_" . $1 . "__";
                $species =~ s/^([kpcofg])__/$nextlevel{$1}__/;
                $info[0] .= "|" . $species;
                last;
            }
        }
        print OT "\n", $info[0];

        foreach $i (@column) {
            print OT "\t", $info[$i];
        }
    }
    print OT "\n";
    close OT;
    close PR;
}

#===============================================================================
sub open_phenotype {
    my ( $phenotype, $subsample, $column, $arrange ) = @_;
    my ( @info, $head, $col, $i, %grouporder );
    open LS, "<$phenotype" or die "read $phenotype $!\n";
    chomp( $head = <LS> );
    @info = split /\t/, $head;
    unless ($column) {
        $col = 1;
    }
    else {
        for ( $i = 1 ; $i < @info ; ++$i ) {
            if ( $info[$i] eq $column ) {
                $col = $i;
                last;
            }
        }
        unless ($col) {
            &usage("column title for LEFSe not matched: $column \n$head");
        }
    }
    if ($arrange) {
        &arrangement( $arrange, \%grouporder );
    }
    while (<LS>) {
        chomp;
        next if (/^#/);
        @info = split /\t/;
        $subsample->{ $info[0] } =
          ($arrange) ? $grouporder{ $info[$col] } : $info[$col];
    }
    close LS;
}

#===============================================================================
sub arrangement {
    my ( $arrange, $grouporder ) = @_;
    my ( @info, $i );
    $arrange =~ s/,$//g;
    @info = split /,/, $arrange;
    for ( $i = 0 ; $i < @info ; ++$i ) {
        $grouporder->{ $info[$i] } = $i + 1 . "." . $info[$i];
    }
}

#===============================================================================
sub usage {
    my ($reason) = @_;
    print STDERR "
  ==============================================================================
  $reason
  ==============================================================================
  " if ( defined $reason );
    print STDERR "
  Last modify: $Version
  Contact: $Contact

  Usage:
  \$perl $0 [options]
  -i  --input     input profile table, required
  -p  --phenotype input phenotype file, required
  -o  --prefix    output result prefix, required
  -c  --column    column title for phenotype file, default 2nd column
  -r  --arrange   group arrangement, connected by ',', default ASCII
  -w  --width     width of image for cladogram, default 11
  -e  --execute   execute [prefix.sh], default not
  -h  --help
  \n";
    exit;
}

#===============================================================================
