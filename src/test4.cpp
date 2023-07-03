#include "cgalPolynomials.h"

// [[Rcpp::export]]
Rcpp::StringMatrix test4(
  Rcpp::IntegerMatrix Powers, Rcpp::NumericVector Coeffs
) {

  PT_9::Construct_polynomial construct_polynomial;

  int nterms = Coeffs.size();
  std::list<Monomial_9> terms;

  for(int i = 0; i < nterms; i++) {
    Rcpp::IntegerVector powers = Powers(Rcpp::_, i);
    terms.push_back(
      std::make_pair(
        CGAL::Exponent_vector(powers.begin(), powers.end()), 
        Coeffs(i)
      )
    );
  }
  Poly_9 P = construct_polynomial(terms.begin(), terms.end());

  std::list<Monomial_9> monoms;
  PT_9::Monomial_representation mrepr;
  mrepr(P, std::back_inserter(monoms));

  XYZ6 Result;
  Poly_6 w0    = PT_6::Shift()(Poly_6(1), 1, 0); // 
  Poly_6 sqrt3 = PT_6::Shift()(Poly_6(1), 1, 1); // 
  Poly_6 A     = PT_6::Shift()(Poly_6(1), 1, 2); // 
  Poly_6 B     = PT_6::Shift()(Poly_6(1), 1, 3); // 
  Poly_6 C     = PT_6::Shift()(Poly_6(1), 1, 4); // 
  Poly_6 D     = PT_6::Shift()(Poly_6(1), 1, 5); // 

  std::list<Monomial_9>::iterator it_monoms;
  for(it_monoms = monoms.begin(); it_monoms != monoms.end(); it_monoms++) {
    CGAL::Exponent_vector allPowers = (*it_monoms).first;
    double coef = (*it_monoms).second;
    Expo3 powersXYZ = {allPowers[0], allPowers[1], allPowers[2]};
    Poly_6 PAB = 
      coef * CGAL::ipower(w0, allPowers[3]) * CGAL::ipower(sqrt3, allPowers[4]) * 
      CGAL::ipower(A, allPowers[5]) * CGAL::ipower(B, allPowers[6]) * 
      CGAL::ipower(C, allPowers[7]) * CGAL::ipower(D, allPowers[8]);
    if(Result.count(powersXYZ) == 0) {
      Result.emplace(powersXYZ, PAB);
    } else {
      Result.at(powersXYZ) += PAB;
    }
  }

  int nmonoms = Result.size();
  Rcpp::StringMatrix POVRay(2, nmonoms);
  int i = 0;
  for(const auto& [key, value] : Result) {
    std::string powers = "xyz(" + 
      std::to_string(std::get<0>(key)) + ", " +
      std::to_string(std::get<1>(key)) + ", " +
      std::to_string(std::get<2>(key)) + "): ";
    std::stringstream buffer;
    CGAL::IO::set_pretty_mode(buffer);
    buffer << value << "," << std::endl;
    std::string coefab = buffer.str();
    Rcpp::StringVector monom = Rcpp::StringVector::create(powers, coefab);
    POVRay(Rcpp::_, i++) = monom;
  }

  return Rcpp::transpose(POVRay);
}