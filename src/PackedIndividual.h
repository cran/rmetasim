/*

$Modified: astrand$

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the site object type.

This is the object that defines an individual.

there is a base object: Individual



*/

#ifndef  PACKED_INDIVIDUAL_H
#define PACKED_INDIVIDUAL_H

/*
includes
*/

#include <metasim.h>
#include <TransMat.h>
#include <FastAllele.h>
#include <FastSeqAllele.h>
#include <AlleleTbl.h>

using namespace std;

/**

The indivdual class is essentially abstract.  It mainly contains a
demographic class variable to pass to its descendents

 */
class PackedIndividual {

private:

  /// time click in which the individual is born
  int gen;
  /// the time click the ind was last modified by survival.
  int changed;
  /// the time the individual last reproduced
  int lastrep;
  /// the number of offspring produced last gen
  int noff;
  ///demographic age or stage of the individual
  short cl;
  ///sexuality 0=herm, 1=female, 2=male
  short sex;
  ///The number of loci
  int nloc;
  ///array of the ploidy for each locus
  short PL[MAXLOCI];
  ///linearized matrix of diploid genotypes of the ind.
  short G[MAXLOCI * MAXPLOIDY];

public:
  PackedIndividual(int c=0, int sx=0, int g=0, int nl=0);
  ~PackedIndividual();

  ///randomly sets the class of an individual (for certain types of initialization)
  int RandomizeClass(int numclass=1);

  void SetClass(int newclass=0);
  int GetClass();

  void SetSex(int newsex=0);
  int GetSex();

  void SetGen(int newgen=0);
  int GetGen();

  void SetLoci(AlleleLookTbl &Atbls);

  void resetLoci(AlleleLookTbl &Atbls);

  inline  int GetLastRep()
    {
      return lastrep;
    }

  inline void SetLastRep(int lr)
    {
      lastrep = lr;
    }
  inline  int GetNumOff()
    {
      return noff;
    }

  inline void SetNumOff(int nf)
    {
      noff = nf;
    }

  inline  int GetAllele(int l, int a)
    {
      return G[ ((l * MAXPLOIDY) + a) ];
    }

  inline void SetAllele(int l, int a, int al)
    {
      assert (al<MAXALLELES);
      G[ ((l * MAXPLOIDY) + a) ] = al;
    }


  inline void WriteState(int l, int a, AlleleLookTbl &Atbls, ostream & stream = cout)
    {
      Atbls[l]->WriteAlleleState(G [ ((l * MAXPLOIDY) + a) ], stream );
    }


  inline void swap_allele(int l)
    {
      int tmpi;
      tmpi = G[((l * MAXPLOIDY) + 0)];
      G[((l * MAXPLOIDY) + 0)] = G[((l * MAXPLOIDY) + 1)];
      G[((l * MAXPLOIDY) + 1)] = tmpi;
    }
  ///sets the survival change flag to true
  inline void Change(int t)
    {
      changed = t;
    }
  inline int GetChanged()
    {
      return changed;
    }

  PackedIndividual MakeGamete(AlleleLookTbl &Atbls);


  int GetRandAlleleIndex(int l);
  void SetRandGenotype(AlleleLookTbl &Atbls);
  int IsGenotypeSet();

  /**
     Updates the allele tables.  If t>0 then also runs the mutation algorithms on each allele.
   */
  void Birth(int t,AlleleLookTbl &Atbls);
  /**
     Used for transferring individuals among classes
   */
  void Growth(AlleleLookTbl &Atbls);
  /**
     Updates the allele tables as an individual is removed from the population
   */
  void Death(int t, AlleleLookTbl &Atbls);

  PackedIndividual repro_sex(PackedIndividual & SO1, PackedIndividual & SO2, int t, AlleleLookTbl &Atbls);
  PackedIndividual repro_asex(PackedIndividual & SO, int t=0);


  friend ostream & operator<<(ostream & stream, PackedIndividual &ind);
  friend istream & operator>>(istream & stream, PackedIndividual &ind);

  size_t Sizeof();


}; // end PackedIndividual


#endif /*PACKED_INDIVIDUAL*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
