use strict;
use warnings;
use File::Basename;

if (@ARGV != 6)
{
	die "Usage: perl $0 hb_universe mm_universe sb_universe DEG_dir DEG_list Outdir\n";
}

my $argc = 0;
my $hb_universe = $ARGV[$argc++];
my $mm_universe = $ARGV[$argc++];
my $sb_universe = $ARGV[$argc++];
my $DEG_dir = $ARGV[$argc++];
my $DEG_list = $ARGV[$argc++];
my $Outdir = $ARGV[$argc++];

my @hb_universe_gene = ();
open IN,$hb_universe;
while (my $line = <IN>)
{
	chomp $line;
	push(@hb_universe_gene, $line);
}
close IN;
my $hb_universe_size = scalar @hb_universe_gene;

my @mm_universe_gene = ();
open IN,$mm_universe;
while (my $line = <IN>)
{
	chomp $line;
	push(@mm_universe_gene, $line);
}
close IN;
my $mm_universe_size = scalar @mm_universe_gene;

my @sb_universe_gene = ();
open IN,$sb_universe;
while (my $line = <IN>)
{
	chomp $line;
	push(@sb_universe_gene, $line);
}
close IN;
my $sb_universe_size = scalar @sb_universe_gene;

open IN,$DEG_list;
while (my $line = <IN>)
{
	my %seen = ();
	chomp $line;
	my $deg_file = "$DEG_dir/$line";
	open OUT,">$Outdir/$line";
	open DEG,$deg_file;
	while (my $deg_gene = <DEG>)
	{
		if ($line =~ /hb/)
		{
			my $rand_ind = int(rand($hb_universe_size));
			while ($seen{$rand_ind})
			{
				$rand_ind = int(rand($hb_universe_size));
			}
			$seen{$rand_ind} = 1;
			print OUT "$hb_universe_gene[$rand_ind]\n";
		}
		if ($line =~ /mm/)
		{
			my $rand_ind = int(rand($mm_universe_size));
			while ($seen{$rand_ind})
			{
				$rand_ind = int(rand($mm_universe_size));
			}
			$seen{$rand_ind} = 1;
			print OUT "$mm_universe_gene[$rand_ind]\n";
		}
		if ($line =~ /sb/)
		{
			my $rand_ind = int(rand($sb_universe_size));
			while ($seen{$rand_ind})
			{
				$rand_ind = int(rand($sb_universe_size));
			}
			$seen{$rand_ind} = 1;
			print OUT "$sb_universe_gene[$rand_ind]\n";
		}
	}
	close OUT;
	close DEG;
}
close IN;