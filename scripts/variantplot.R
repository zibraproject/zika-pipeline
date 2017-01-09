library(ggplot2)
variants=read.table('Control.variants.tab', sep='\t', header=T)
p=ggplot(variants, aes(x=pos, y=supportfraction, color=vartype)) + geom_point() + facet_wrap(~set, ncol=1)
ggsave('Control.variants.png', p)
