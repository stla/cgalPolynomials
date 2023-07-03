library(spray)

x <- lone(1, 5)
y <- lone(2, 5)
z <- lone(3, 5)
a <- lone(4, 5)
b <- lone(5, 5)

P <- ((x*x+y*y+1)*(a*x*x+b*y*y)+z*z*(b*x*x+a*y*y)-2*(a-b)*x*y*z-a*b*(x*x+y*y))^2-4*(x*x+y*y)*(a*x*x+b*y*y-x*y*z*(a-b))^2

Coeffs <- as.integer(P[["value"]])
Powers <- t(P[["index"]])

pov <- cgalPolynomials:::test3(Powers, Coeffs)

ab <- pov[, 2L]
ab <- gsub("\n$", "", gsub("y", "b", gsub("x", "a", ab)))
ab <- gsub("([ab])\\^(\\d+)", "pow(\\1,\\2)", x = ab)

prcode <- paste0(pov[, 1L], ab)

writeLines(prcode, "SMS2.txt")
