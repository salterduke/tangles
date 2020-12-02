#!/usr/bin/perl -w
# Transfer the probilities of leaf node sets to corresponding protein pairs

open(IN, "lns_count.txt") || die "Can't open input file!\n";
while (<IN>) {
	chomp;
	($count, $pos, $neg, @lns) = split;
	$lns = join "<->", @lns;
	$count{$lns} = $count;
}
close(IN);

open(IN, "all_lns.txt") || die "Can't open input file!\n";
while (<IN>) {
	chomp;
	($orf_1, $orf_2, $int, @lns) = split;
	$orf_2 =~ s/://;
	$lns = join "<->", @lns;
	print "$orf_1\t$orf_2\t$count{$lns}\n";
}
close(IN);

