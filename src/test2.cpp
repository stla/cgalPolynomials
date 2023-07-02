#include "cgalPolynomials.h"

// [[Rcpp::export]]
void test2(Rcpp::IntegerMatrix Powers, Rcpp::IntegerVector Coeffs) {

  CGAL::IO::set_pretty_mode(std::cout);

  int nterms = Coeffs.size();
  std::list<Monomial_5> terms;

  for(int i = 0; i < nterms; i++) {
    Rcpp::IntegerVector powers = Powers(Rcpp::_, i);
    terms.push_back(
      std::make_pair(
        CGAL::Exponent_vector(powers.begin(), powers.end()), 
        Coeffs(i)
      )
    );
  }
  Poly_5 P = construct_polynomial(terms.begin(), terms.end());

  std::list<Monomial_5> monoms;
  PT_5::Monomial_representation mrepr;
  mrepr(P, std::back_inserter(monoms));

  XYZ Result;
  Poly_2 A = PT_2::Shift()(Poly_2(1), 1, 0); // A
  Poly_2 B = PT_2::Shift()(Poly_2(1), 1, 1); // B

  std::list<Monomial_5>::iterator it_monoms;
  for(it_monoms = monoms.begin(); it_monoms != monoms.end(); it_monoms++) {
    CGAL::Exponent_vector allPowers = (*it_monoms).first;
    int coef = (*it_monoms).second;
    Expo3 powersXYZ = {allPowers[0], allPowers[1], allPowers[2]};
    Poly_2 PAB = 
      coef * CGAL::ipower(A, allPowers[3]) * CGAL::ipower(B, allPowers[4]);
    if(Result.count(powersXYZ) == 0) {
      Result.emplace(powersXYZ, PAB);
    } else {
      Result.at(powersXYZ) += PAB;
    }
  }

  for(const auto& [key, value] : Result) {
    std::cout << "x^" << std::get<0>(key) 
              << "y^" << std::get<1>(key) 
              << "z^" << std::get<2>(key) 
              << ": " << value << "\n";
  }

}