#  Allan Strand 9/27/01
#
#

#
#default seed is same as calling environment.  Zero or any positive integer
#is assigned to R's random number seed.  The type of RNG is inherited from the
#calling environment
#
simulate.landscape <- function(Rland, numit, seed=-1, compress=FALSE, adj.lambda=0)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
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
        .Call("compress_landscape",Rland,PACKAGE = "rmetasim")
      }
    else
      {
        print("Rland not a landscape object...exiting")
        Rland
      }
  }
