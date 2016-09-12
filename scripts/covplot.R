#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
args=commandArgs(TRUE)
a=read.table(args[1], header=F)
b=summarise(group_by(a, V4, V5, V11, V12), n=n())
p=ggplot(subset(b, n>5), aes(x=V11, xend=V12, y=n, yend=n, label=V4)) + geom_segment() + theme_bw(base_size=14) + xlab("Postion in genome") + ylab("Number of reads") + ggtitle(args[1])
ggsave(paste(args[1],'.png',sep=''), width=10, height=6)
