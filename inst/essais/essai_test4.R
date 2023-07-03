library(spray)

f <- function(x, y, z, w, sqrt3) {
  ((x^2+y^2+z^2+w^2+145/3)^2-4*(9*z^2+16*w^2))^2*((x^2+y^2+z^2+w^2+145/3)^2+296*(x^2+y^2)-4*(9*z^2+16*w^2))-16*(x^2+y^2)*(x^2+y^2+z^2+w^2+145/3)^2*(37*(x^2+y^2+z^2+w^2+145/3)^2-1369*(x^2+y^2)-7*(225*z^2+448*w^2))-sqrt3*16/9*(x^3-3*x*y^2)*(110*(x^2+y^2+z^2+w^2+145/3)^3-148*(x^2+y^2+z^2+w^2+145/3)*(110*x^2+110*y^2-297*z^2+480*w^2))-64*(x^2+y^2)*(3*(729*z^4+4096*w^4)+168*(x^2+y^2)*(15*z^2-22*w^2))+64*(12100/27*(x^3-3*x*y^2)^2-7056*(3*x^2*y-y^3)^2)-592240896*z^2*w^2
}

x <- lone(1, 9)
y <- lone(2, 9)
z <- lone(3, 9)
w0 <- lone(4, 9)
sqrt3 <- lone(5, 9)
a <- lone(6, 9)
b <- lone(7, 9)
c <- lone(8, 9)
d <- lone(9, 9)

# substitutions to be done
X = a*x - b*y - c*z - d*w0
Y = a*y + b*x + c*w0 - d*z
Z = a*z - b*w0 + c*x + d*y
W = a*w0 + b*z - c*y + d*x

# perform the substitutions
P = f(X, Y, Z, W, sqrt3)


#
Coeffs <- P[["value"]]
Powers <- t(P[["index"]])

pov <- cgalPolynomials:::test4(Powers, Coeffs)

ab <- pov[, 2L]
ab <- gsub("z", "d",
           gsub("y", "e",
                gsub("x", "f",
                     gsub("w3", "a",
                          gsub("\n$", "",
                               gsub("w2", "b",
                                    gsub("w1", "c", ab)))))))
ab <- gsub("([abcdef])\\^(\\d+)", "pow(\\1,\\2)", x = ab)

prcode <- paste0(pov[, 1L], ab)

writeLines(prcode, "Julia.txt")
