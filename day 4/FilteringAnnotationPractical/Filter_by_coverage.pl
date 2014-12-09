use strict;
use warnings;

my $infile = "GATK_filtered.vcf";
my $outfile = "GATK_coverage_filtered.vcf";
open(OUT, ">$outfile");

my $count_pass=0;
my $count_filt=0;

open(FILE, $infile)||die "Cannot open $infile";
while(<FILE>)
{
	chomp;
	if ($_ =~ /^#/) ### ignore and print header lines to outfile
	{
		print OUT "$_\n";
	}
	else
	{
		my @line = split /\t/, $_;
		my $info = $line[9]; ### select last column of file with genotype and allelic depth information
		my ($gt, $ad, $gq, $pl) = split /:/, $info; ### split column by ':'
		my @depths = split /,/, $ad; ### split allelic depth into depths for each allele
		my $total_cov = 0;
		foreach my $depth (@depths)
		{
			$total_cov = $total_cov + $depth; ### sum allelic depths to get total coverage
		}
		if ($total_cov < 50){} ### remove if coverage is <50
		else
		{
			print OUT "$_\n"; ### print high coverage variants to outfile
		}
	}
}
