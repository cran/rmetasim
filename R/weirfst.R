#
# Implements an ANOVA approach to estimating pop structure statistics.  Implemented
# straight from Weir 1996 Genetic Data Analysis
#



# These functions depend upon a data frame in the following format
# which is produced by metatrans with the 'r' input options
#
# The class field can be effectively ignored for the current functions
# but it could be required in the future.
#
#pop class individual locus aindex allele
#000 001 000 000 000 064
#000 001 000 001 000 002
#000 001 000 002 000 001
#000 001 000 002 001 007
#000 001 000 003 000 002
#000 001 000 003 001 002
#000 001 000 004 000 047
#000 001 000 004 001 563
#000 001 000 005 000 000
#000 001 000 005 001 000


                                        #
#
#this function returns the ploidy of the loci in the genotypic matrix in an ordered vector
weir.ploidy <- function(genmat)
  {
    attach(genmat)
    pl <- apply((table(locus,aindex)>0),1,sum)
    detach(genmat)
    pl
  }

#constructs an indicator variable for each allele,individual,locus,population
#
indicator <- function(genmat,loc)
{
    genmat<-genmat[genmat$locus==loc,]
    attach(genmat)
    indic<- table(aindex,allele,individual,pop)
    detach(genmat)
    indic
}

#sums indicator value across alleles within individuals
#expects an indicator matrix created by indicator

allelesum <- function(indicatormat)
{
  if (dim(indicatormat)[1]<2)
    {
      print ("must send data from a diploid locus to this function")
    }
  im<-indicatormat[1,,,] + indicatormat[2,,,]
  if (length(dim(im))<3)
    {
      im<-array(im,c(1,dim(im)))
      dimnames(im)<-dimnames(indicatormat)[2:4]
    }
  im
}
#sums indicator value across alleles within individuals
#expects an indicator matrix created by alleleindicator


#
#
#
allelevec <- function(allele,allelevec)
  {
    1 * (allelevec==allele)
  }

#
#makes an alleleindicator dataframe that is suitable for anova, etc
#
alleleindicator <-function(genmat)
{
  ual <- unique(genmat$allele)
  m <- matrix(genmat$allele,ncol=1)
  d <- as.data.frame(t(apply(m,1,allelevec,ual)))
  names(d)<-make.names(ual)
  cbind(genmat,d)
}

#
# calculate the allele frequencies for each population for each locus
# expects an alleleindicator matrix (output of alleleindicator)
#

alfreqp <- function (alind)
  {
    pl<-unique(alind$pop)
    ll<-unique(alind$locus)
    retmat<-matrix(0,nrow=(length(pl)*length(ll)),ncol=(ncol(alind)-4))
    t<-1
    for (i in ll)
      {
        thislocus<-alind[alind$locus==i,]
        for (j in pl)
          {
            retmat[t,]<-array(c(i,j,apply(thislocus[thislocus$pop==j,][,7:ncol(thislocus)],2,mean)))
            t<-t+1
          }
        rm <-as.data.frame(retmat)
        names(rm)<-make.names(c("locus","pop",names(alind)[7:ncol(alind)]))
      }
    rm
  }

#calculates the average allele freq for each locus across populations
#expects alleleindicator output
#
alfreql <- function (alind)
  {
    ll<-unique(alind$locus)
    retmat<-matrix(0,nrow=(length(ll)),ncol=(ncol(alind)-4))
    t<-1
    for (i in ll)
      {
        thislocus<-alind[alind$locus==i,]
        retmat[t,]<-array(c(i,1,apply(thislocus[,7:ncol(thislocus)],2,mean)))
        t<-t+1
        rm <-as.data.frame(retmat)
        names(rm)<-make.names(c("locus","dummy",names(alind)[7:ncol(alind)]))
      }
    rm
  }

#this function takes a group of rows representing all individuals sampled for a locus
#from an allelefreq mat and eliminates the alleles not involved in this locus
cleanallelefrq <- function (afreqmat) # this is the result of afreq
  {
    al<-apply(afreqmat,2,sum)
    ab<- (al>0)
    ab[1]<-TRUE
    ab[2]<-TRUE
    afreqmat[,ab]
  }

cleanalleleind <- function(alind)
  {
    al<-apply(alind,2,sum)
    ab<- (al>0)
    ab[1:6]<-TRUE
    alind[,ab]
  }

#determines the expected proportion of the homozygotes for each allele in a all pops
Hexp <- function(genin)
	{
		allelefreq(genin) ^ 2
	}



#determines the observed proportion of the homozygotes for each allele in a all pops

Homobs <- function(alind)
  {
    asum<-aggregate(alind[,7:22],list(alind$locus,alind$pop,alind$individual),sum)
    nm<-c("locus","pop","individual",names(asum)[4:ncol(asum)])
    names(asum)<-nm
    asum<-cbind(asum[,1:3],(asum[,4:ncol(asum)]==2)*1)
    homobs<-aggregate(asum[,4:ncol(asum)],list(asum$pop,asum$locus),mean)
    nm<-c("pop","locus",names(homobs)[3:ncol(homobs)])
    names(homobs)<-nm
    homobs
  }

                                        #
#functions that calculate population structure statistics
#
#
#


#
# uses the approach on pg 173 of Weir (1996) to calculate theta for
# haploid data
#
# takes a genmat (from metatrans)  that has been trimmed to include only haploid loci
# and returns an estimate for each locus in an array
#
hapstruct <- function (hapmat)
{
  alind <- alleleindicator(hapmat)
  alfrqp <- alfreqp(alind)
  alfrql <- alfreql(alind)

  npops <- length(unique(alind$pop))
  nsize <- array(0,npops)
  T1array <- array(0,unique(alind$locus))
  T2array <- array(0,unique(alind$locus))
  t<-1
  nss <-0
  for (i in unique(alind$pop))
    {
      nsize[t]<-length(unique(alind[alind$pop==i,]$individual))
      nss <- nss + nsize[t]^2
      t<-t+1
    }
  nmean <- mean(nsize)

  nc <- (1/(npops - 1)) * (sum(nsize) - nss/sum(nsize))
  
#this is the loop across all of the loci

  li <- 1 #index of T1 and T2 arrays
  for (i in unique(alfrql$locus))
    {
      tmp1<- alfrql[alfrql$locus==i,]
      T1<-0
      T2<-0
      allelenames <- names(cleanallelefrq(alfrql[alfrql$locus==i,]))
      allelenames <- allelenames[3:length(allelenames)]
      #loop over alleles within a locus
      for (al in allelenames)
        {
          #now loop across populations
          sinter <- 0
          t<-1
          for (j in unique(alfrqp$pop))
            {
              tmp2<- tmp1[tmp1$pop==j,]
              pAj <- tmp2[,names(tmp2)==al]
              pAdot <- alfrql[alfrql$locus==i,names(alfrql)==al]
              sinter<-sinter + (nsize[t] * (pAj - pAdot)^2)
              t<-t+1
            }
          s2A <- (1/((npops-1)*nmean)) * sinter
          
          T1 <- T1 + s2A - (1/(nmean-1)) * (pAdot*(1-pAdot) - ((npops-1)/npops)*s2A)
          T2 <- T2 + ((nc - 1)/(nmean-1))*(pAdot*(1-pAdot))+(1+((npops-1)*(nmean-nc)/(nmean - 1)))*(s2A/npops)
        }
      
      T1array[li]<-T1
      T2array[li]<-T2
      li<-li+1
    }

#  print(T1array)
#  print(T2array)

  theta <- T1array/T2array
  for (i in 1:length(theta))
    {
      if (is.na(theta[i]))
      {
        theta[i]<-0
      }
    }
  theta
}




#
# uses the approach on pg 173 of Weir (1996) to calculate f, theta, Fit for
# diploid data
#
# takes a genmat (from metatrans)  that has been trimmed to include only diploid loci
# and returns an estimate for each locus in an array
#
dipstruct <- function (dipmat)
{
  alind <- alleleindicator(dipmat)

#  print("calculated indicator matrix")
  
  npops <- length(unique(alind$pop))
  nsize <- array(0,npops)
  t<-1
  nss <-0
  for (i in unique(alind$pop))
    {
      nsize[t]<-length(unique(alind[alind$pop==i,]$individual))
      nss <- nss + nsize[t]^2
      t<-t+1
    }
  
#  print("calculated nss")

  nmean <- mean(nsize)
  nc <- (1/(npops - 1)) * (sum(nsize) - nss/sum(nsize))
#this is the loop across all of the loci
  li <- 1 #index of S arrays

  S1array<-c(rep(0,length(unique(alind$locus))))
  S2array<-S1array
  S3array<-S1array

#  print("Looping across loci")

  for (i in unique(alind$locus))
    {
      tmp1<- cleanalleleind(alind[alind$locus==i,])
#      print(tmp1)
      S1<-0
      S2<-0
      S3<-0

#      print("Looping across alleles")
      
      for (al in names(tmp1)[7:ncol(tmp1)])
        {
          selvec<-names(tmp1)==al
          selvec[1:6]<-TRUE
#          print(selvec)
          tmp2<-tmp1[,selvec]
#          print(tmp2)
          if (length(unique(tmp2$aindex))==1)
            {
              anmat<-anova(lm(tmp2[,7]~as.factor(tmp2$pop)))
#              print(anmat)
              S1<- S1 + anmat[1,3] - anmat[2,3]
              S2<- S2 + anmat[1,3] + (nc-1)*anmat[2,3]
              S3<- NA
            }
          else
            {
              anmat<-anova(lm(tmp2[,7]~as.factor(tmp2$pop)/as.factor(tmp2$individual)))
#              print(anmat)
              S1<- S1 + anmat[1,3] - anmat[2,3]
              S2<- S2 + anmat[1,3] + (nc-1)*anmat[2,3]+nc*anmat[3,3]
              S3<- S3 + 2*nc*anmat[3,3]
            }
        }
      S1array[li]<-S1
      S2array[li]<-S2
      S3array[li]<-S3
      li<-li+1
    }


  theta <- S1array/S2array
  Fhat <- S3array/S2array
  for (i in 1:length(theta))
    {
      if (is.na(theta[i]))
      {
        theta[i]<-0
      }
    }
  cbind(theta,Fhat)
}

popstruct<-function(genin)
  {
    dipstruct(genin)
  }

autopopstruct<-function()
  {
    genin<- read.table(file="tmp.rin",header=TRUE)
    dipstruct(genin)
  }

summarize <-function()
  {
    tmp<-read.table(file="theta.summary",header=TRUE)
    thetas <- tmp[order(tmp$locus,tmp$gen,tmp$iteration),]
    postscript(file="theta.ps",horizontal=F,paper="letter")
    attach(thetas)
    coplot(theta~gen|as.factor(locus)*as.factor(iteration),show.given=F,xlab=c("Generation","Columns: Loci"),ylab=(c("\Theta","Rows: Independent iteration")),panel=lines)
    dev.off()

    #calculate some summary stats and then graphically report them
    thetameans <- aggregate(thetas,list(gen,locus),mean)[,c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE)]
    theta25<- aggregate(thetas$theta,list(gen,locus),quantile,1/4)$x
    theta75<- aggregate(thetas$theta,list(gen,locus),quantile,3/4)$x
    thetasd<- aggregate(thetas$theta,list(gen,locus),sd)
    n <- length(unique(thetas$iteration))
    df <- n-1
    t <- qt(0.975,df)
    stderr <-  thetasd$x / n
    attach(thetameans)
    conf.low <- theta-(t*stderr)
    conf.high<- theta+(t*stderr)
    data.frame(gen,locus,theta,stderr,conf.low,conf.high,df,theta25,theta75)
  }

plottmeans <- function()
  {
    thetasum<-summarize()
    postscript(file="thetameans.ps",paper="letter")
    layout(matrix(unique(thetasum$locus),ceiling(length(unique(thetasum$locus))/2),2,T))
    cvec <- c("A","B","C","D","E","F","G","H","I","J","K","L");
    for (i in unique(thetasum$locus))
      {
        plotsum<-thetasum[(locus==i),]
        attach(plotsum)
        ymax <- max(theta75,theta25,theta)
        ymin <- min(theta75,theta25,theta)
        plot(gen,theta,type="l",lwd=2,xlim=c(10000,30000),ylim=c(ymin,ymax),ylab="theta",xlab="Time")
        text(10050,(0.9*ymax),cvec[i])
#        lines(gen,conf.low,lty=3)
#        lines(gen,conf.high,lty=3)
        lines(gen,theta25,lty=3)
        lines(gen,theta75,lty=3)
        detach(plotsum)
      }
    dev.off()
  }

####################
# this function just converts a landscape object into a data frame that can be analyzed by R
#
#

landscape.to.rtheta.old <- function(Rland,loci=NULL,pops=NULL,samp=20)
  {
    if (is.null(loci))
      {
        loci <- 1:length(Rland$loci)
      }
    gn <- NULL
    p <- landscape.populations(Rland)

    for (loc in loci)
      {
        if (landscape.ploidy(Rland)[loc]==1)
          {
            thislocus <- cbind(p,rep(loc,length(p)),rep(0,length(p)),c(1:length(p)),landscape.locus(loc,Rland)[,c(1,4)])
          }
        else
          {
            thislocus <- rbind(cbind(p,rep(loc,length(p)),rep(0,length(p)),c(1:length(p)),landscape.locus(loc,Rland)[,c(1,4)]),cbind(p,rep(loc,length(p)),rep(1,length(p)),c(1:length(p)),landscape.locus(loc,Rland)[,c(1,5)]))
          }
        gn <- rbind(gn,thislocus)
      }
    gn <- as.data.frame(gn)
    names(gn) <- c("pop","locus","aindex","individual","class","allele")
    gn <- gn[,c(1,5,4,2,3,6)] #keep the columns in a familiar order
    if (is.null(pops))
      {
        pops <- unique(p)
      }

    np <- NULL
    for (pop in pops)
      {
        pm <- gn[gn$pop%in%pop,]
        if (length(unique(pm$individual))>samp)
          {
            pm <-pm[pm$individual%in%sample(unique(pm$individual),samp),] #sample inds
          }
        np <- rbind(np,pm)
      }
    np
  }


#this version is implemented mostly in C++ and is faster than the R version.
landscape.to.rtheta <- function(Rland, numi=0)
  {
    if (is.landscape(Rland))
      {
        rv <- data.frame(matrix(.Call("l2w",Rland,numi,PACKAGE = "rmetasim"),ncol=6,byrow=TRUE))
        names(rv) <- c("pop","class","individual","locus","aindex","allele")
        rv$pop <- rv$pop+1
        rv$locus <- rv$locus+1
        rv
      }
  }
