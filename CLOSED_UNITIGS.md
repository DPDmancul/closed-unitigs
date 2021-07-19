# Closed Unitigs

## Support
We define the support of a _k_-mer as its count:
> given _m_ ∈ _K_ (set of _k_-mers)  
> supp(_m_) = # occourencies of _m_ in _K_

and the support of a unitig as the minimum support of its _k_-mers, which is equal to the minimum support of a substring of it:
> supp(_u_) = min<sub>_m_ ⊆ _u_</sub> supp(_m_)  
> <small>(⊆ denotes a substring of at least _k_ characters)</small>

### Antimonotonicity of the support
> _m_ ⊆ _u_ ⇒ supp(_u_) ⩽ supp(_m_)

From the definition of support.

## Closed unitig
A unitig is closed if adding a _k_-mer to it its support will be lower:
> _u_ is closed ⇔ supp(_v_) < supp(_u_) ∀ _v_ ⊃ _u_

### Property
For each _k_-mer _m_ there exists a closed unitig _u_ with the same support.

Dim.  
  * _m_ closed, then _u_=_m_;
  * _m_ not closed, then ∃ _v_ ⊃ _m_ s.t. supp(_v_) = supp(_m_): reiterate for _v_.

This always terminates since the unitig containg all _k_-mers is closed.

## Lossless compression via closed unitigs
From the set of closed untigs with their supports we can extract all the _k_-mers and their counts:
  * All the substrings of lenght _k_ of the closed unitigs are the _k_-mers;
  * Each _k_-mer has a closed untig with the same support, and since the support of a unitig is the same of minimum support of its _k_-mers, then the support of a _k_-mer is the maximum support of a closed unitig containing it. 