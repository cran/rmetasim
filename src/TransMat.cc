/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the site object type.


In this file the transition matrix class is defined
*/


/*
includes
*/

#include <TransMat.h>

TransMat::TransMat (size_t s )
{
  size_t lt;
  size = s;
  tm.resize(s);
  for (lt=0; lt<s; lt++)
    {
      tm[lt].resize(s);
    }
}

TransMat::~TransMat ()
{
#ifdef RDEBUG
   cout << "destructing TransMat" << "\n";
#endif
   size_t i;
   for (i=0;i<tm.size();i++)
     {
       tm[i].resize(0);
     }
   tm.resize(0);
#ifdef RDEBUG
   cout << "finished destructing TransMat" << "\n";
#endif
}



void TransMat::SetSize(size_t sz) 
{
  size_t i;

  size = sz; 

#ifdef RDEBUG
  cerr << "Resizing a first dimension of a transmatrix of size"<<sz <<endl;
#endif

  tm.resize(sz);

#ifdef RDEBUG
  cerr << "Finished resizing the first dimension" <<endl;
#endif
  for (i=0;i<sz;i++)
    {
      tm[i].resize(sz);
    }
}

///Sets an entire TransMat of size s from a 2d array pointed to by a.  The colummns should represent from and rows, to
void TransMat::SetMat(TransMat a)
{
  size_t i,j;
  for (i=0;i<a.Size();i++)
    {
      for (j=0;j<a.Size();j++)
	{
	  SetElement(j,i,a.GetElement(j,i));
	}
    }
  /// don't have an error function here probably need to check for reasonable values.
}

//

/*
TransMat & TransMat::operator= (TransMat &T)
{
  int i, j;
  int s = T.Size();
  for (i=0; i<s; i++)
    for (j=0; j<s; j++)
      {
	this->SetElement(i,j,T.GetElement(i,j));
      }
  return *this;
}
*/


///Implementation of the random state algorithm
void TransMat::SetRandomToStateVec (double eigenratio)
{
  size_t sz = Size();
  double p[sz + 1];

  size_t i;

  for (i=0;i<sz;i++)
    {
      SetToState(i);
      p[i]=Value()*eigenratio;
      assert(p[i]>=0);
    }
  RandLibObj.SetDiscreteLookup(p,sz+1);
}

///Implementation of the random state algorithm
void TransMat::SetRandomFromStateVec ()
{
  size_t sz = Size();
  double p[sz + 1];

  size_t i;

  for (i=0;i<sz;i++)
    {
      SetFromState(i);
      p[i]=Value();
      assert(p[i]>=0);
    }
  RandLibObj.SetDiscreteLookup(p,sz+1);
}

size_t TransMat::RandomState()
{
  SetToState(RandLibObj.PickMultinomial());
  if (Size() == GetToState()) 
    { 
      SetToState(-1); 
    }
  return GetToState();
}



size_t TransMat::PoissonOffspring(double eigenratio)
{
  float val;
  size_t num = 0;
  val = Value()*eigenratio;
  if (val>0.0)
    {
      num = RandLibObj.poisson(val);
    }
  return num;
}

int TransMat::AnyFrom(size_t fs)
{
  size_t i;
  double tot = 0.0;
  SetFromState(fs);
  for (i=0;i<size;i++)
    {
      SetToState(i);
      tot = tot + Value();
    }
  return (tot>0);
}

ostream &operator<<(ostream &stream, TransMat & TM)
{
  size_t i,j;
  stream.precision(3);
  stream << TM.Size() << endl;
  for (i=0;i<TM.Size();i++)
    {
      TM.SetToState(i);
      for (j=0;j<TM.Size();j++)
	{
	  TM.SetFromState(j);
	  stream << TM.Value() << " ";
	}
      stream << endl;
    }
  return stream;
}

istream &operator>>(istream &stream, TransMat & TM)
{
  size_t i,j, n=199;
  stream >> n;
  TM.SetSize(n);
  for (j=0;j<n;j++)
    {
      for (i=0;i<n;i++)
	{
	  stream >> TM.tm[j][i] ;
	}
    }
  return stream;
}


ostream &operator<<(ostream & stream, DemoVec &d)
{
  int i,sz;
  sz = d.v.size();
  stream << sz << endl ;
  for (i=0;i<sz;i++)
    {
      stream << sz << " " ;
    }
  stream << endl ;
  return stream  ; 
}

istream &operator>>(istream & stream, DemoVec &d)
{
  
  int i,sz;
  int tmp;

  sz = d.size()   ;
  stream >> sz    ;
  d.resize(sz)    ;
  for (i=0;i<sz;i++)
    {
      stream >> tmp;
      d.Set(tmp,i) ;
    }

  return stream ; 
}



/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
