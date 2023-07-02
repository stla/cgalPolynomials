library(spray)

x <- lone(1, 5)
y <- lone(2, 5)
z <- lone(3, 5)
a <- lone(4, 5)
b <- lone(5, 5)

P <- ((x*x+y*y+1)*(a*x*x+b*y*y)+z*z*(b*x*x+a*y*y)-2*(a-b)*x*y*z-a*b*(x*x+y*y))^2-4*(x*x+y*y)*(a*x*x+b*y*y-x*y*z*(a-b))^2

Coeffs <- as.integer(P[["value"]])
Powers <- t(P[["index"]])

cgalPolynomials:::test2(Powers, Coeffs)

