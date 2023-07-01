#ifndef __HEADER__
#define __HEADER__

#include <Rcpp.h>

#define CGAL_EIGEN3_ENABLED 1

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/polynomial_utils.h>
#include <CGAL/Polynomial/Monomial_representation.h>

typedef CGAL::Polynomial_type_generator<int, 5>::Type Poly_5;
typedef CGAL::Polynomial_traits_d<Poly_5>             PT_5;
//typedef PT_5::Coefficient_type                        PT_4;  
//typedef PT_5::Innermost_coefficient_type              Integer;

typedef CGAL::Polynomial_type_generator<int, 2>::Type Poly_2;
typedef CGAL::Polynomial_traits_d<Poly_2>             PT_2;

typedef std::tuple<int, int, int> Expo3;
typedef std::map<Expo3, Poly_2>   XYZ;

#endif
