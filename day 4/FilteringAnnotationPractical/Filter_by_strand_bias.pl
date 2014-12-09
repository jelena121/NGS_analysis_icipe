use strict;
use warnings;

my $infile = "GATK_coverage_filtered.vcf";
my $outfile = "FINAL_filtered.vcf";
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
		my $ref_cov = $depths[0];
		my $alt_cov = $depths[1];
		my $balance;
		if ($ref_cov == 0)
		{
			$balance = 1;
		}
		else
		{
			$balance = $alt_cov/($alt_cov+$ref_cov); ### calculate allelic balance
		}
		if ((($gt =~ /0\/1/)||($gt =~ /1\/0/))&&(($balance < 0.3)||($balance > 0.7))){} ### ignore if heterozygote with allelic balance < 0.3 or > 0.7
		elsif ((($gt =~ /1\/1/)||($gt =~ /0\/0/))&&($balance < 0.8)){} ### ignore if homozygote with allelic balance < 0.8
		else
		{
			print OUT "$_\n"; ### print OK variants to outfile
		}
	}
}
