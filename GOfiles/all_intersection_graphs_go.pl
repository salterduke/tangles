#!/usr/bin/perl -w
# Compute the leaf node sets (and intersection graphs) for all gene pairs
# Count how many times leaf node sets occur among all gene pairs

open(IN, "process_orf_graph.txt");
while (<IN>) {
	chomp;
	($orf, @edge) = split;
	for $edge (@edge) {
		$graph{$orf}{$edge} = 1;
	}
}
close(IN);


# All ORFs in data
@orf = sort(keys %graph);
$N = $#orf;


# Loop over all ORF pairs

open(INT, ">all_intersection_graph.txt");
open(LEA, ">all_lns.txt");
open(COU, ">lns_count.txt");
for $i (0 .. $N) {
	$fraction_done = 100.0*$i/$N;
	printf("%.2f ", $fraction_done); print "%\n";
	
	# Get graph for $orf[$i]
	@edge = sort(keys %{$graph{$orf[$i]}});
	
	for $j ($i + 1 .. $N) {
		
		# Find common edges -> intersection
		
		@intersection = ();
		%parent = ();
		%child = ();
		for $edge (@edge) {
			if (defined $graph{$orf[$j]}{$edge}) {
				push @intersection, $edge;
				
				($parent, $child) = split /->/, $edge;
				$parent{$parent} = 1;
				$child{$child} = 1;
			}
		}
		unless (defined $intersection[0]) {
			$intersection[0] = "NA";
		} 
		
		
		# Identify leaf nodes (must be a subset of the children nodes)
		
		@lns = ();
		@child = sort(keys %child);
		for $child (@child) {
			unless (defined $parent{$child}) {	# If the child is NOT a parent
				push @lns, $child;
			}
		}
		unless (defined $lns[0]) {
			$lns[0] = "NA";
		} 		
		
		
		# Are ORFs a MIPS complexes pair or a negative pair? -- deleted by HY
#		($orf[$i], $orf[$j]) = sort($orf[$i], $orf[$j]);
#		$pair = join ":", ($orf[$i], $orf[$j]); 
		$int = 0;
		
		print INT "$orf[$i]\t$orf[$j]\t$int\t@intersection\n";
		print LEA "$orf[$i]\t$orf[$j]\t$int\t@lns\n";
		
		
		# Count the frequency of specific leaf node sets
		
		$leaf_node_set = join "<->", @lns;
		$count{$leaf_node_set}++;
		
		if ($int == 1) {
			$count_int{$leaf_node_set}++;
		} elsif ($int == -1) {
			$count_neg{$leaf_node_set}++;
		}
	}
}
close(INT);
close(LEA);


# Output of leaf node set frequencies

@leaf_node_set = sort(keys %count);
for $leaf_node_set (@leaf_node_set) {
	if (defined $count_int{$leaf_node_set}) {
		$count_int = $count_int{$leaf_node_set};
	} else {
		$count_int = 0;
	}

	if (defined $count_neg{$leaf_node_set}) {
		$count_neg = $count_neg{$leaf_node_set};
	} else {
		$count_neg = 0;
	}

	print COU "$count{$leaf_node_set}\t$count_int\t$count_neg\t$leaf_node_set\n";
}
close(COU);
