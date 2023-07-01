library(Ryacas)

p <- 
  "x^2 + a*x^2 + 2/3*y + b*y + x*z + a*x*z"

yac_str("p := x^2 + a*x^2 + 2/3*y + b*y + 5*x*z + a*x*z")
yac_str("Coef(p, x, 2)")
yac_str("Expand(p, {x,y,z})")


library(Ryacas)
library(partitions)

yac_str("poly := x^2 + a*x^2 + 2/3*y + b*y + 5*x*z + a*x*z")

f <- function(comp){
  xyz <- sprintf("xyz(%s): ", toString(comp))
  coef <- yac_str(sprintf(
    "Expand(Coef(Coef(Coef(poly, x, %d), y, %d), z, %d), {a, b})",
    comp[1L], comp[2L], comp[3L]
  ))
  if(coef == "0") return(NULL)
  paste0(xyz, coef, ",")
}
for(deg in 0L:2L){
  comps <- compositions(deg, 3L)
  povray <- apply(comps, 2L, f, simplify = FALSE)
  cat(
    unlist(povray), sep = "\n", file = "povray.txt", append = deg > 0L
  )
}

# yac_str("Set(X, a*x - b*y - c*z - d*w)")

library(Ryacas)
library(partitions)

yac_str("XX := a*x - b*y - c*z - d*w")
yac_str("Y  := a*y + b*x + c*w - d*z")
yac_str("Z  := a*z - b*w + c*x + d*y")
yac_str("W  := a*w + b*z - c*y + d*x")
expr <- "((XX^2+Y^2+Z^2+W^2+145/3)^2-4*(9*Z^2+16*W^2))^2*((XX^2+Y^2+Z^2+W^2+145/3)^2+296*(XX^2+Y^2)-4*(9*Z^2+16*W^2)) -16*(XX^2+Y^2)*(XX^2+Y^2+Z^2+W^2+145/3)^2*(37*(XX^2+Y^2+Z^2+W^2+145/3)^2-1369*(XX^2+Y^2)-7*(225*Z^2+448*W^2)) -16*sqrt(3)/9*(XX^3-3*XX*Y^2)*(110*(XX^2+Y^2+Z^2+W^2+145/3)^3 -148*(XX^2+Y^2+Z^2+W^2+145/3)*(110*XX^2+110*Y^2-297*Z^2+480*W^2)) -64*(XX^2+Y^2)*(3*(729*Z^4+4096*W^4)+168*(XX^2+Y^2)*(15*Z^2-22*W^2)) +64*(12100/27*(XX^3-3*XX*Y^2)^2 -7056*(3*XX^2*Y-Y^3)^2) -592240896*Z^2*W^2"
yac_assign(expr, "POLY")
yac_str("Degree(POLY, x)") # 12
yac_str("Degree(POLY, y)") # 12
yac_str("Degree(POLY, z)") # 12
yac_str("Degree(POLY, w)") # 12

f <- function(comp){
  xyzw <- sprintf("xyzw(%s): ", toString(comp))
  coef <- yac_str(sprintf(
    "ExpandBrackets(Coef(Coef(Coef(Coef(POLY, x, %d), y, %d), z, %d), w, %d))",
    comp[1L], comp[2L], comp[3L], comp[4L]
  ))
  if(coef == "0") return(NULL)
  paste0(xyzw, coef, ",")
}

for(deg in 0L:12L){
  comps <- compositions(deg, 4L)
  povray <- apply(comps, 2L, f, simplify = FALSE)
  cat(
    unlist(povray), sep = "\n", file = "povray.txt", append = deg > 0L
  )
}

# solid Möbius strip
library(Ryacas)
library(partitions)

expr <- "((x*x+y*y+1)*(a*x*x+b*y*y)+z*z*(b*x*x+a*y*y)-2*(a-b)*x*y*z-a*b*(x*x+y*y))^2-4*(x*x+y*y)*(a*x*x+b*y*y-x*y*z*(a-b))^2"
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