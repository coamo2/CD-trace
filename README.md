# CD-trace
 Compositional Data Network Analysis via Lasso Penalized D-Trace Loss


a)	The R code named “ mycd.R” is an implementation of Algorithm 1 in our manuscript.

b)	The folder named “SpiecEasi-master” is an implementation of S-E(glasso) and S-E(mb) (Kurtz et al, 2015), which is available at https://github.com/zdk123/SpiecEasi by Kurtz et al (2015). We also use the R function “make_graph” in this package to generate the hub, block and scale-free graph in our simulations. The R function “graph2prec” is used to convert graph topologies into the precision matrix.

c)	The R code named “gcoda.R” is an implementation of gCoda (Fang et al, 2017), which is available from https://github.com/huayingfang/gCoda by Fang et al (2017).

d)	The R code named “sim.R” is for simulations under different scenarios. The graph topologies are generated and converted into precision matrix with function “make_graph” and “graph2prec” in SpiecEasi-master. The networks are estimated according to CD-trace, gCoda, S-E(glasso) and S-E(mb) with the help of above-mentioned codes. The codes for AUC calculation, ROC figures and tables are also included in this file. 
