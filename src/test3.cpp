#include "cgalPolynomials.h"

// [[Rcpp::export]]
Rcpp::StringMatrix test3(
  Rcpp::IntegerMatrix Powers, Rcpp::IntegerVector Coeffs
) {

  PT_5::Construct_polynomial construct_polynomial;

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