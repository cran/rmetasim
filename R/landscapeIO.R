read.landscape <- function(fn = "filename")
  {
    if (file.exists(fn))
      {
        .Call("read_landscape",fn,PACKAGE = "rmetasim")
      }
    else
      {
        print (paste("Filename: ",fn,"does not appear to exist"))
        NULL
      }
  }


write.landscape <- function(rland, fn = "filename")
  {
    if (is.landscape(rland))
      .Call("write_landscape",fn,rland,PACKAGE = "rmetasim")
  }


