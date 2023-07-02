#include "cgalPolynomials.h"

// [[Rcpp::export]]
void test() {

  CGAL::IO::set_pretty_mode(std::cout);

  Poly_5 x = PT_5::Shift()(Poly_5(1), 1, 0); // x_0^1
  Poly_5 y = PT_5::Shift()(Poly_5(1), 1, 1); // x_1^1
  Poly_5 z = PT_5::Shift()(Poly_5(1), 1, 2); // x_2^1
  Poly_5 a = PT_5::Shift()(Poly_5(1), 1, 3); // x_3^1
  Poly_5 b = PT_5::Shift()(Poly_5(1), 1, 4); // x_4^1

  Poly_5 P = CGAL::ipower(x, 8)*CGAL::ipower(a, 2)+2*CGAL::ipower(x, 6)*CGAL::ipower(y, 2)*CGAL::ipower(a, 2)+2*CGAL::ipower(x, 6)*CGAL::ipower(y, 2)*a*b+(-2)*CGAL::ipower(x, 6)*CGAL::ipower(a, 2)*b-2*CGAL::ipower(x, 6)*CGAL::ipower(a, 2)+2*CGAL::ipower(x, 6)*a*b*CGAL::ipower(z, 2)+4*CGAL::ipower(x, 5)*y*CGAL::ipower(a, 2)*z-4*CGAL::ipower(x, 5)*y*a*b*z+CGAL::ipower(x, 4)*CGAL::ipower(y, 4)*CGAL::ipower(a, 2)+4*CGAL::ipower(x, 4)*CGAL::ipower(y, 4)*a*b+CGAL::ipower(x, 4)*CGAL::ipower(y, 4)*CGAL::ipower(b, 2)+(-4)*CGAL::ipower(x, 4)*CGAL::ipower(y, 2)*CGAL::ipower(a, 2)*b+(-2)*CGAL::ipower(x, 4)*CGAL::ipower(y, 2)*CGAL::ipower(a, 2)*CGAL::ipower(z, 2)-2*CGAL::ipower(x, 4)*CGAL::ipower(y, 2)*CGAL::ipower(a, 2)+(-2)*CGAL::ipower(x, 4)*CGAL::ipower(y, 2)*a*CGAL::ipower(b, 2)+10*CGAL::ipower(x, 4)*CGAL::ipower(y, 2)*a*b*CGAL::ipower(z, 2)-4*CGAL::ipower(x, 4)*CGAL::ipower(y, 2)*a*b-2*CGAL::ipower(x, 4)*CGAL::ipower(y, 2)*CGAL::ipower(b, 2)*CGAL::ipower(z, 2)+CGAL::ipower(x, 4)*CGAL::ipower(a, 2)*CGAL::ipower(b, 2)-2*CGAL::ipower(x, 4)*CGAL::ipower(a, 2)*b+CGAL::ipower(x, 4)*CGAL::ipower(a, 2)+(-2)*CGAL::ipower(x, 4)*a*CGAL::ipower(b, 2)*CGAL::ipower(z, 2)+2*CGAL::ipower(x, 4)*a*b*CGAL::ipower(z, 2)+CGAL::ipower(x, 4)*CGAL::ipower(b, 2)*CGAL::ipower(z, 4)+4*CGAL::ipower(x, 3)*CGAL::ipower(y, 3)*CGAL::ipower(a, 2)*z-4*CGAL::ipower(x, 3)*CGAL::ipower(y, 3)*CGAL::ipower(b, 2)*z+4*CGAL::ipower(x, 3)*y*CGAL::ipower(a, 2)*b*z-4*CGAL::ipower(x, 3)*y*CGAL::ipower(a, 2)*z+(-4)*CGAL::ipower(x, 3)*y*a*CGAL::ipower(b, 2)*z+(-4)*CGAL::ipower(x, 3)*y*a*b*CGAL::ipower(z, 3)+4*CGAL::ipower(x, 3)*y*a*b*z+4*CGAL::ipower(x, 3)*y*CGAL::ipower(b, 2)*CGAL::ipower(z, 3)+2*CGAL::ipower(x, 2)*CGAL::ipower(y, 6)*a*b+2*CGAL::ipower(x, 2)*CGAL::ipower(y, 6)*CGAL::ipower(b, 2)+(-2)*CGAL::ipower(x, 2)*CGAL::ipower(y, 4)*CGAL::ipower(a, 2)*b-2*CGAL::ipower(x, 2)*CGAL::ipower(y, 4)*CGAL::ipower(a, 2)*CGAL::ipower(z, 2)+(-4)*CGAL::ipower(x, 2)*CGAL::ipower(y, 4)*a*CGAL::ipower(b, 2)+10*CGAL::ipower(x, 2)*CGAL::ipower(y, 4)*a*b*CGAL::ipower(z, 2)-4*CGAL::ipower(x, 2)*CGAL::ipower(y, 4)*a*b+(-2)*CGAL::ipower(x, 2)*CGAL::ipower(y, 4)*CGAL::ipower(b, 2)*CGAL::ipower(z, 2)-2*CGAL::ipower(x, 2)*CGAL::ipower(y, 4)*CGAL::ipower(b, 2)+2*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*CGAL::ipower(a, 2)*CGAL::ipower(b, 2)+(-2)*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*CGAL::ipower(a, 2)*b*CGAL::ipower(z, 2)-2*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*CGAL::ipower(a, 2)*b+6*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*CGAL::ipower(a, 2)*CGAL::ipower(z, 2)+(-2)*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*a*CGAL::ipower(b, 2)*CGAL::ipower(z, 2)-2*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*a*CGAL::ipower(b, 2)+2*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*a*b*CGAL::ipower(z, 4)-8*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*a*b*CGAL::ipower(z, 2)+2*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*a*b+6*CGAL::ipower(x, 2)*CGAL::ipower(y, 2)*CGAL::ipower(b, 2)*CGAL::ipower(z, 2)+4*x*CGAL::ipower(y, 5)*a*b*z-4*x*CGAL::ipower(y, 5)*CGAL::ipower(b, 2)*z+4*x*CGAL::ipower(y, 3)*CGAL::ipower(a, 2)*b*z-4*x*CGAL::ipower(y, 3)*CGAL::ipower(a, 2)*CGAL::ipower(z, 3)+(-4)*x*CGAL::ipower(y, 3)*a*CGAL::ipower(b, 2)*z+4*x*CGAL::ipower(y, 3)*a*b*CGAL::ipower(z, 3)-4*x*CGAL::ipower(y, 3)*a*b*z+4*x*CGAL::ipower(y, 3)*CGAL::ipower(b, 2)*z+CGAL::ipower(y, 8)*CGAL::ipower(b, 2)+(-2)*CGAL::ipower(y, 6)*a*CGAL::ipower(b, 2)+2*CGAL::ipower(y, 6)*a*b*CGAL::ipower(z, 2)-2*CGAL::ipower(y, 6)*CGAL::ipower(b, 2)+CGAL::ipower(y, 4)*CGAL::ipower(a, 2)*CGAL::ipower(b, 2)-2*CGAL::ipower(y, 4)*CGAL::ipower(a, 2)*b*CGAL::ipower(z, 2)+CGAL::ipower(y, 4)*CGAL::ipower(a, 2)*CGAL::ipower(z, 4)+(-2)*CGAL::ipower(y, 4)*a*CGAL::ipower(b, 2)+2*CGAL::ipower(y, 4)*a*b*CGAL::ipower(z, 2)+CGAL::ipower(y, 4)*CGAL::ipower(b, 2);

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
    //std::cout << "exponent: "<< (*it_monoms).first << std::endl;
    //std::cout << "coef: "<< (*it_monoms).second << std::endl;
    //std::cout << "exponents of a and b: "<< (*it_monoms).first[3] << " and " << (*it_monoms).first[4] << std::endl;
  }

  for(const auto& [key, value] : Result) {
    std::cout << "x^" << std::get<0>(key) 
              << "y^" << std::get<1>(key) 
              << "z^" << std::get<2>(key) 
              << ": " << value << "\n";
  }

  PT_5::Degree degree;
  std::cout << "degree of P with respect to x: "<< degree(P, 0) << std::endl;
  std::cout << "degree of P with respect to y: "<< degree(P, 1) << std::endl;
  std::cout << "degree of P with respect to z: "<< degree(P, 2) << std::endl;
  std::cout << "degree of P with respect to a: "<< degree(P, 3) << std::endl;
  std::cout << "degree of P with respect to b: "<< degree(P, 4) << std::endl;

  /*PT_5::Leading_coefficient lcoeff;
  std::cout << "Leading coefficient with respect to x:           "
            << lcoeff(P)
            << std::endl;

  PT_5::Get_coefficient get_coefficient;
  std::cout << "Coefficient of x: "<< get_coefficient(P, 0) << std::endl;
  std::cout << "Coefficient of y: "<< get_coefficient(P, 1) << std::endl;
  */

  PT_5::Get_innermost_coefficient get_icoeff;
  std::vector<int> exp = {4, 0, 2, 1, 2};
  std::cout << "Innermost coefficient of monomial x^4 y^0 z^2 a^1 b^2:        "
            << get_icoeff(P, CGAL::Exponent_vector(exp.begin(), exp.end())) 
            << std::endl;  
}