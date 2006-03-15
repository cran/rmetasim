#
#
#
#


#He
exp.het.landscape <- function(Rland)
  {
    tot <- 0
    rl <- matrix(0,nrow=Rland$intparam$habitats,ncol=length(Rland$loci))
    for (j in 1:length(unique(populations(Rland))))
      {
        for (loc in 1:length(Rland$loci))
          {
            tab <- table(landscape.locus(loc,Rland)[populations(Rland)==j,c(-1:-(landscape.democol()))])
            sctab <- tab/sum(tab)
            rl[j,loc] <- 1 - sum(sctab^2)
          }
      }
    rl
  }

#Ho
obs.het.landscape <- function(Rland)
  {
    tot <- 0
    rl <- matrix(0,nrow=Rland$intparam$habitats,ncol=length(Rland$loci))
    for (j in unique(populations(Rland)))
      {
        for (loc in 1:length(Rland$loci))
          {
            if (ploidy(Rland)[loc]==1) #obs het doesn't make sense for a haploid locus
              {
                rl[j,loc] <- NA
              }  else {
                freq.df <- data.frame(table(landscape.locus(loc,Rland)[populations(Rland)==j,c(-1,-2,-3,-5)],landscape.locus(loc,Rland)[populations(Rland)==j,c(-1,-2,-3,-4)]))
                rl[j,loc] <- (1-sum(freq.df[as.character(freq.df[,1])==as.character(freq.df[,2]),3])/sum(freq.df[,3]))
              }
          }
      }
    rl
  }


FWright.landscape <- function (Rland)
  {
    1-obs.het.landscape(Rland)/exp.het.landscape(Rland)
  }


allelefreq.landscape <- function(Rland,tbl.out=FALSE)
  {
    rv <- NULL
    for (i in 1:length(Rland$loci))
      {
        pops <- vector("list",length(Rland$loci))
        for (j in unique(populations(Rland)))
          {
            alleles     <- landscape.locus(i,Rland)[populations(Rland)==j,c(-1:-(landscape.democol()))]
            freqtbl     <- table(alleles)
            scframe     <- data.frame(freqtbl/sum(freqtbl))
            scframe$pop <- rep(j,dim(scframe)[1])
            scframe$loc <- rep(i,dim(scframe)[1])
            rv <- rbind(rv,scframe)
          }
      }
    rownames(rv) <- 1:dim(rv)[1]
    if (tbl.out==TRUE)
      {
        xtabs(Freq~pop+alleles+loc,rv)
      } else {
      rv
    }
  }

allelecount.landscape <- function(Rland,tbl.out=FALSE)
  {
    rv <- NULL
    for (i in 1:length(Rland$loci))
      {
        pops <- vector("list",length(Rland$loci))
        for (j in unique(populations(Rland)))
          {
            alleles     <- landscape.locus(i,Rland)[populations(Rland)==j,c(-1:-(landscape.democol()))]
            freqtbl     <- table(alleles)
            scframe     <- data.frame(freqtbl)
            scframe$pop <- rep(j,dim(scframe)[1])
            scframe$loc <- rep(i,dim(scframe)[1])
            rv <- rbind(rv,scframe)
          }
      }
    rownames(rv) <- 1:dim(rv)[1]
    if (tbl.out==TRUE)
      {
        xtabs(Freq~pop+alleles+loc,rv)
      } else {
      rv
    }
  }
