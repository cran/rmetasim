seedMigrationRate <- 1
pollenMigrationRate <- 1
habitats <- 10
stages <- 2

rland <- NULL
rland <- new.landscape.empty()
rland <- new.intparam.land(rland, h=habitats, s=stages, totgen=5000)
rland <- new.switchparam.land(rland, mp=0)
rland <- new.floatparam.land(rland)

S <- matrix(c(0,0,1,0), nrow = 2, byrow = TRUE)
R <- matrix(c(0,1.1,0,0), nrow = 2, byrow = TRUE)
M <- matrix(c(0,0,0,1), nrow = 2, byrow = TRUE)
rland <- new.local.demo(rland,S,R,M)

rland <- new.epoch.island(rland,0,c(0,0),c(0,0),
                          seedMigrationRate, c(1,0), c(1,0),
                          pollenMigrationRate, c(0,1), c(0,1))

rland <- new.locus(rland,type=2,ploidy=1,transmission=1,numalleles=4,allelesize=100)
for (x in 1:9) {
  rland <- new.locus(rland,type=2,ploidy=2,transmission=0,numalleles=4,allelesize=100)
}
for (x in 1:10) {
  rland <- new.locus(rland,type=1,ploidy=2,transmission=0,numalleles=4)
}

# add individuals
rland <- new.individuals(rland,
                         c(5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0))
