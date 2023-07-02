library(Ryacas)


# solid Möbius strip
library(Ryacas)
library(partitions)
library(spray)
x <- lone(1, 5)
y <- lone(2, 5)
z <- lone(3, 5)
a <- lone(4, 5)
b <- lone(5, 5)

expr <- ((x*x+y*y+1)*(a*x*x+b*y*y)+z*z*(b*x*x+a*y*y)-2*(a-b)*x*y*z-a*b*(x*x+y*y))^2-4*(x*x+y*y)*(a*x*x+b*y*y-x*y*z*(a-b))^2
powers <- expr$index
coeffs <- as.integer(expr$value)

yac_assign(expr, "POLY")
yac_str("Degree(POLY, x)") # 8
yac_str("Degree(POLY, y)") # 8
yac_str("Degree(POLY, z)") # 4




gsub("(\\(.+?\\))\\^(\\d+)", "CGAL::ipower(\\1, \\2)", x = expr)

expanded <- yac_str("ExpandBrackets(POLY)")
gsub("([xyzab])\\^(\\d+)", "CGAL::ipower(\\1, \\2)", x = expanded)

f <- function(comp){
  if(comp[3L] > 4L) return(NULL)
  xyz <- sprintf("xyz(%s): ", toString(comp))
  coef <- yac_str(sprintf(
    "ExpandBrackets(Coef(Coef(Coef(POLY, x, %d), y, %d), z, %d))",
    comp[1L], comp[2L], comp[3L]
  ))
  if(coef == "0") return(NULL)
  coef <- gsub("([ab])\\^(\\d+)", "pow(\\1,\\2)", x = coef)
  paste0(xyz, coef, ",")
}

for(deg in 0L:8L){
  comps <- compositions(deg, 3L)
  povray <- apply(comps, 2L, f, simplify = FALSE)
  cat(
    unlist(povray), sep = "\n", file = "SMS.txt", append = deg > 0L
  )
}





yac_str("aaa := Coef(POLY, a, {0,1,2})")
yac_str("Coef(aaa[1], b, {0,1,2})")
yac_str("Coef(aaa[2], b, {0,1,2})")
yac_str("Coef(aaa[3], b, {0,1,2})")
