#!/usr/bin/env Rscript

## clustertree.R
##
## implement Prosperi's clustering algorithm on a tree.
## source article:
## A novel methodology for large-scale phylogeny partition.
## http://www.ncbi.nlm.nih.gov/pubmed/21610724
##
## code by glenn lawyer October 2012
######################################################
require(ape, quietly=TRUE)     ## ape and gieger to read the tree
require(geiger, quietly=TRUE)
require(igraph, quietly=TRUE)  # igraph for the search

## Prosperi method:
## Perform DFS. Label a subtree a cluster if both
##   -its median pairwise patristic ditance (MPPD) is below a percentile
##       of the whole-tree p-distanwww.scipirate.com  ce and
##   -it is not a leaf. 

## Given
##  a tree and a threshold, and (optionally)
##  a vector of the MPPD of the subtree at each internal node, 
##  the reliability threshold and vector,
##  and specify if you want the membership for all nodes or
##  only the tips
## Return
##  a vector indicating, for each internal node, which
##  cluster it belongs to
prosperi.cluster <- function(tree,thresh,distvec=NULL,
                             rthresh=NULL,reliabilityvec=NULL,
                             retval=c("tips","all")){
  ## take care of optional arguments
  retval <- match.arg(retval)
  if(is.null(distvec)){
    #cat("Calculating MPPD for each node ...\n")
    distvec <- pdist.clusttree(tree,mode='all')
  }
  if(is.null(rthresh) || is.null(reliabilityvec)){
    reliabilityvec=NULL
    #cat("No reliability thresholding.\n")
  }
  ## set up clustering
  ntips<-Ntip(tree)
  cnum <- 0 ## cluster number
  assign <- rep(0,ntips+tree$Nnode) ## cluster assignment
  igraph.tree <- graph.edgelist(tree$edge) ## tree in igraph form
  dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out',
                   order=TRUE,dist=TRUE)
  ## travese the tree in depth first order
  for(i in 1:length(dfs$order)){
    node <- dfs$order[i]
    ## skip leaves
    if(node < ntips+1){ next }
    ## skip unreliable nodes (if reliability measure is available)
    if(! is.null(reliabilityvec) &&
       reliabilityvec[node-ntips] >= rthresh){ next } 
    ## If the node's subtree is below the threshold, mark it and
    ## its subtree as members of a new cluster
    if(distvec[node-ntips]<=thresh && assign[node]<=0){
      cnum <- cnum+1
      subtree <- graph.dfs(igraph.tree,node,
                           neimode='out',unreachable=FALSE)$order
      subtree <- subtree[! is.na(subtree)]
      assign[subtree] <- cnum
    }}
  ans <- list(membership=assign,allcsize=table(assign),
              leafclustsize=table(assign[1:ntips]),
              ntips=ntips,threshold=thresh)
  if(retval=="tips"){ans$membership <- ans$membership[1:ntips]}
  class(ans) <- c(class(ans),'p.cluster')
  return(ans)
}



testplot <- function(testtree,thresh,plotfile=NULL){
  graphics.off()
  clustering <- prosperi.cluster(testtree,thresh)$membership
  if(class(plotfile)=="character"){ png(plotfile) }    
  plot(testtree,show.tip.label=FALSE)
  tiplabels(rep("  ",Ntip(testtree)),
            bg=clustering)
  if(class(plotfile)=="character"){ dev.off() }

}


##############################################################
##############################################################
## Three helper functions

## Given
##   a node, tree, and distance matrix
## Return
##   median pairwise patristic distance (MPPD) of its leaves
get.node.leaf.MPPD <- function(node,tree,distmat){
  nlist <- tips(tree,node)
  foo <- distmat[nlist,nlist] 
  return(median(foo[upper.tri(foo,diag=FALSE)]))
}


## Given
##   a node, tree, and distance matrix
## Return
##   median pairwise patristic distance (MPPD) of all of its decendants
get.node.full.MPPD <- function(node,tree,distmat){
  nlist <- tips(tree,node)
  elist <- tree$edge[which.edge(tree,nlist),2]
  foo <- distmat[elist,elist] 
  return(median(foo[upper.tri(foo,diag=FALSE)]))
}


## Given
##   a tree and a distance matrix
## Return
##   a vector giving the median pairwise
##   patristic distance of the subtree under
##   each internal node
## SLOW!! May be a good idea to save/cache results
pdist.clusttree <- function(tree,distmat=NULL,mode=c('leaf','all')){
  mode <- match.arg(mode)
  if(is.null(distmat)){
    if(mode=='leaf'){ distmat <-  cophenetic.phylo(tree) }
    else{ distmat <-  dist.nodes(tree) }
  }
  ntips<- Ntip(tree)
  nint <- tree$Nnode ## number of internal nodes
  if(mode=='leaf'){
    return(sapply((ntips+1):(ntips+nint),get.node.leaf.MPPD,tree,distmat))
  }
  else{
    return(sapply((ntips+1):(ntips+nint),get.node.full.MPPD,tree,distmat))
  }
}

args=commandArgs(TRUE)
tree = read.tree(args[1])
clust=prosperi.cluster(tree, as.numeric(args[2]), retval='tips')

df <- data.frame(tree$tip.label, clust$membership)
names(df) <- c("node", "group")
write.table(df, "", sep="\t", row.names = F)

