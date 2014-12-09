use strict;
use warnings;

my $infile = "variants.vcf";
my $outfile = "GATK_filtered.vcf";
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
		if ($line[6] =~ /PASS/) ### look for PASS in FILTER column
		{
			print OUT "$_\n"; ### print PASS variants to outfile
			$count_pass++; ### count variants that PASS filters
		}
		else
		{
			$count_filt++; ### count variants that are filtered out
		}
	}
}

print "$count_pass varaints PASSED, $count_filt varaints Filtered\n"; ### print to screen numbers of filtered and retained variants
