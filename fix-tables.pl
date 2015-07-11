#!/usr/bin/env perl


my $str = do { local $/; <STDIN> };

$str =~ s/\\usepackage{longtable}//g;
$str =~ s/\\endhead//g;
$str =~ s/\\caption{([^}]*)}\n\\end{longtable}/\\end{longtable}\n\\bottomcaption{$1}/gm;
$str =~ s/\\begin{longtable}/\\begin{center}\\begin{supertabular}/g;
$str =~ s/\\end{longtable}/\\end{supertabular}\\end{center}/g;

print($str);
