#  Allan Strand 9/27/01
#
#

###
### coerce landscape tries to make sure that the landscape object is in a form that is palatable to 
### C++ code in rmetasim 
###

coerce.landscape <- function(rland)
  {
    rland$intparam <- lapply(rland$intparam,as.integer)
    rland$switchparam <- lapply(rland$switchparam,as.integer)
    rland$loci <- lapply(rland$loci,function(x)
                         {
                           x$type <- as.integer(x$type)
                           x$ploidy <- as.integer(x$ploidy)
                           x$trans <- as.integer(x$trans)
                           x$alleles <- lapply(x$alleles,function(y,typ)
                                               {
                                                 y$aindex <- as.integer(y$aindex)
                                                 y$birth <- as.integer(y$birth)
                                                 if (typ!=253)
                                                   {
                                                     y$state <- as.integer(y$state)
                                                   }
                                                 y
                                               },
                                               typ=x$type)
                           x
                         })
    rland$individuals <- matrix(as.integer(rland$individuals),nrow=dim(rland$individuals)[1])
    rland
  }
####
#default seed is same as calling environment.  Zero or any positive integer
#is assigned to R's random number seed.  The type of RNG is inherited from the
#calling environment
#
sim.landscape <- function(Rland, numit, seed=-1, compress=FALSE, adj.lambda=0)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- coerce.landscape(Rland)
        .Call("iterate_landscape",as.integer(numit),Rland,as.integer(compress),as.integer(adj.lambda),PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }


survive.landscape <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- coerce.landscape(Rland)
        .Call("survive_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

reproduce.landscape <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- coerce.landscape(Rland)
        .Call("reproduce_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }


carry.landscape <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- coerce.landscape(Rland)
        .Call("carry_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

extinct.landscape <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- coerce.landscape(Rland)
        .Call("extinct_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

advance.landscape <- function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland <- coerce.landscape(Rland)
        .Call("advance_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

compress.landscape <- function(Rland)
    {
    if (is.landscape(Rland))
      {
        Rland <- coerce.landscape(Rland)
        .Call("compress_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
        Rland
      }
  }
