#
#
#
#routines that operate on loci, both the actual states and the indices stored in
#individual genotypes

#returns a individual x ploidy matrix of aindices
landscape.locus <- function(lnum=1,Rland)
  {
    if(is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          Rland$individuals[,c(TRUE,TRUE,TRUE,locusvec(Rland)==lnum)]
        }
  }

#returns a individual x ploidy matrix of states
landscape.states <- function(lnum=1,Rland)
  {
    if(is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          lmat <- as.data.frame(Rland$individuals[,c(TRUE,TRUE,TRUE,locusvec(Rland)==lnum)])
          st <- states(lnum,Rland)
          lmat[,4] <- st$state[c(lmat[,4]+1)]

          if (ploidy(Rland)[lnum]==2)
            {
              lmat[,5] <- st$state[c(lmat[,4]+1)]
            }
          lmat
        }
  }

landscape.states.old <- function(lnum=1,Rland)
  {
    if(is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          lmat <- as.data.frame(Rland$individuals[,c(TRUE,TRUE,TRUE,locusvec(Rland)==lnum)])
          st <- states(lnum,Rland)
          stvec <- rep(NA,max(st$aindex+1))
          for (i in 1:length(stvec))
            {
              stvec[st$aindex[i]+1] <- st$state[i]
            }

          lmat[,4] <- stvec[lmat[,4]+1]

          if (ploidy(Rland)[lnum]==2)
            {
              lmat[,5] <- stvec[lmat[,5]+1]
            }
          lmat
        }
  }


#returns a vector of ploidys for all loci
ploidy<- function(Rland)
  {
    ploidy<-c();
    for (i in 1:Rland$intparam$locusnum)
      {
        ploidy<-c(ploidy,Rland$loci[[i]]$ploidy);
      }
    ploidy
  }

#returns a vector of locus ids
locusvec<- function(Rland)
  {
    p<-ploidy(Rland);
    lv<-c();
    for (i in  1:Rland$intparam$locusnum)
      {
        lv<-c(lv,rep(i,p[i]));
      }
    lv
  }
#
#takes a locus and returns the states and their indices
#
states<-function(lnum=1,Rland)
  {
    if (is.landscape(Rland))
      if (lnum<=Rland$intparam$locusnum)
        {
          ain<-c();
          sta<-c();
          locin <- landscape.locus(lnum,Rland)[,c(-1,-2,-3)]
#          print(locin)
          ainds <- unique(c(locin))
#          print(ainds)
          for (i in 1:length(Rland$loci[[lnum]]$alleles))
            {
              if (Rland$loci[[lnum]]$alleles[[i]]$aindex %in% ainds)
                {
                  ain<-c(ain,Rland$loci[[lnum]]$alleles[[i]]$aindex);
                  sta<-c(sta,Rland$loci[[lnum]]$alleles[[i]]$state);
                }
            }
          list(aindex=ain,state=sta);
          
        }
  }

indxfreq <-function(lnum=1,Rland)
  {
    lv<-landscape.locus(lnum,Rland)[,c(FALSE,FALSE,FALSE,rep(TRUE,(ncol(locus(lnum,Rland))-3)))];
    if (ploidy(Rland)[lnum]==1)
      {
        table(populations(Rland),lv)
      }
    else
      {
        lv2<-c(lv[,1],lv[,2]);
        table(rep(populations(Rland),2),lv2)
      }
    
  }





