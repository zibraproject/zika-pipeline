#!/usr/bin/env Rscript

args=commandArgs(TRUE)
stats = read.table(args[1], header=T, sep="\t")

perc = mean(stats$matches/(stats$matches + stats$insertions + stats$deletions + stats$mismatches))
med = median(stats$matches/(stats$matches + stats$insertions + stats$deletions + stats$mismatches))
mod = mode(stats$matches/(stats$matches + stats$insertions + stats$deletions + stats$mismatches))
ma = max(stats$matches/(stats$matches + stats$insertions + stats$deletions + stats$mismatches))
print(perc)
print(med)
print(mod)
print(ma)

