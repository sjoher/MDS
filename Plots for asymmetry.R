# Probeer paper van Okada (1990) te hanteren voor het berekenen van de afstanden tussen de ellipsen

library(smacof)
library(asymmetry)
library(pracma)
library(plotrix)
options(digits = 16)
data("Englishtowns")

final_config <- rbind(c(-146.005302992741662, 98.991143016455567),
                      c(-90.222656454876500, 49.835161949828340),
                      c(-83.807165965806163, -138.977142208711342),
                      c(18.247421886387102, -38.851441893646516),
                      c(254.055437076186678, 54.100947272837423),
                      c(79.387267439036705, -62.261202558365177),
                      c(112.493005205769634, 17.792023590857671),
                      c(-144.148006193955695, 19.370510830744031))

ind_drift <- function(data, fc){
  symres <- symdecomp(data)
  M <- symres$M
  N <- symres$N
  n <- nrow(M)
  ind <- expand.grid(1:n, 1:n)[, 2:1]
  indmat <- ind[-which(ind[, 1] == ind[, 2]), ]
  Amat <- t(apply(indmat, 1, function(ij) fc[ij[2], ] - fc[ij[1],]))
  Bmat <- t(apply(Amat, 1, function(ab) ab/c(sqrt((t(ab) %*% ab)))))
  diag(N) <- NA
  nij <- as.vector(na.omit(c(N)))
  Cmat <- Bmat * nij
  vectors <- final_config[rep(1:nrow(final_config), times = c(7,7,7,7,7,7,7,7)), ] + Cmat
  return(list(fc, vectors, Cmat))
}

Min_enclosing_ellips <- function(P, tolerance) {
  d = size(P, 1)
  N = size(P, 2)
  
  Q = zeros(d + 1, N)
  Q[1:d, ] = P[1:d, 1:N]
  Q[d + 1, ] = ones(1, N)
  
  count = 1
  err = 1
  u = (1 / N) * rep(1, N)
  
  while (err > tolerance) {
    X = Q %*% diag(u) %*% t(Q)
    M = diag(t(Q) %*% solve(X) %*% Q)
    j = which.max(M)
    maximum = M[j]
    step_size = (maximum - d - 1) / ((d + 1) * (maximum - 1))
    new_u = (1 - step_size) * u
    new_u[j] = new_u[j] + step_size
    count = count + 1
    err = norm(matrix(new_u - u))
    u = new_u
  }
  
  U = diag(u)
  c = P %*% matrix(u)
  A = (1 / d) * solve(P %*% U %*% t(P) - c %*% t(c))
  
  svd_A <- svd(A)
  U <- svd_A$u
  Q <- svd_A$d
  # V = rotation matrix
  V <- svd_A$v
  r1 = 1/sqrt(Q[1])
  r2 = 1/sqrt(Q[2])
  theta <- atan2(V[2,1], V[1,1])
  theta <- rad2deg(theta)
  returnList <- list("center" = c, "a" = r1, "b" = r2, "theta" = theta)
}

names <- expand.grid(rownames(Englishtowns), rownames(Englishtowns))
names <- names[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
result <- ind_drift(Englishtowns, final_config)
rownames(result[[3]]) <- paste(as.character(names[,1]), as.character(names[,2]))
vectors <- result[[2]]


result_kendal <- Min_enclosing_ellips(t(vectors[1:7,]), 0.001)
result_manchester <- Min_enclosing_ellips(t(vectors[8:14,]), 0.001)
result_norwich <- Min_enclosing_ellips(t(vectors[15:21,]), 0.001)
result_oxford <- Min_enclosing_ellips(t(vectors[22:28,]), 0.001)
result_penzance <- Min_enclosing_ellips(t(vectors[29:35,]), 0.001)
result_southampton <- Min_enclosing_ellips(t(vectors[36:42,]), 0.001)
result_taunton <- Min_enclosing_ellips(t(vectors[43:49,]), 0.001)
result_york <- Min_enclosing_ellips(t(vectors[50:56,]), 0.001)

centers <- as.data.frame(matrix(rbind(result_kendal$center, result_manchester$center, result_norwich$center, result_oxford$center, result_penzance$center, result_southampton$center, result_taunton$center, result_york$center),ncol = 2, byrow = T))

rotation <- as.data.frame(matrix(rbind(result_kendal$theta, result_manchester$theta, result_norwich$theta, result_oxford$theta, result_penzance$theta, result_southampton$theta, result_taunton$theta, result_york$theta),ncol = 1, byrow = T))
colnames(rotation) <- c("Angle")
rownames(rotation) <- rownames(Englishtowns)

ellipse_a_b <- as.data.frame(rbind(
  c(result_kendal$a, result_kendal$b),
  c(result_manchester$a, result_manchester$b), 
  c(result_norwich$a, result_norwich$b), 
  c(result_oxford$a, result_oxford$b),
  c(result_penzance$a, result_penzance$b), 
  c(result_southampton$a, result_southampton$b), 
  c(result_taunton$a, result_taunton$b), 
  c(result_york$a, result_york$b))) 
colnames(ellipse_a_b) <- c("a", "b")
rownames(ellipse_a_b) <- rownames(Englishtowns)
#a_b <- ellipse_a_b/2


##### Math stackexchange kendal -> Manchester
a_b <- ellipse_a_b

# Generaized
combinations <- expand.grid(1:8, 1:8)
combinations <- combinations[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
colnames(combinations) <- c("From", "To")
combinations <- combinations[c("To", "From")]
colnames(combinations) <- c("From", "To")
rownames(combinations) <- 1:nrow(combinations)
combinations <- combinations[-c(8,15:16,22:24,29:32,36:40,43:48,50:56),]
rownames(combinations) <- 1:nrow(combinations)

ellips_dist_1 <- numeric(nrow(combinations))
ellips_dist_2 <- numeric(nrow(combinations))
R <- numeric(nrow(combinations))
u <- numeric(nrow(combinations))
U <- numeric(nrow(combinations))
d_x <- numeric(nrow(combinations))
d_y <- numeric(nrow(combinations))
e_x <- numeric(nrow(combinations))
e_y <- numeric(nrow(combinations))
D_x <- numeric(nrow(combinations))
D_y <- numeric(nrow(combinations))
E_x <- numeric(nrow(combinations))
E_y <- numeric(nrow(combinations))

for(i in 1:nrow(combinations)){
  # combinations[i,1] = From
  # combinations[i,2] = To
  R[i] <- ((centers[combinations[i,1],2] - centers[combinations[i,2],2])/(centers[combinations[i,1],1] - centers[combinations[i,2],1]))
  u[i] <- (((cos(deg2rad(rotation[combinations[i,1],1])) + R[i]*sin(deg2rad(rotation[combinations[i,1],1])))^2)/((a_b[combinations[i,1],1])^2)) + (((cos(deg2rad(rotation[combinations[i,1],1])) - R[i]*sin(deg2rad(rotation[combinations[i,1],1])))^2)/((a_b[combinations[i,1],2])^2))
  U[i] <- (((cos(deg2rad(rotation[combinations[i,2],1])) + R[i]*sin(deg2rad(rotation[combinations[i,2],1])))^2)/((a_b[combinations[i,2],1])^2)) + (((cos(deg2rad(rotation[combinations[i,2],1])) - R[i]*sin(deg2rad(rotation[combinations[i,2],1])))^2)/((a_b[combinations[i,2],2])^2))
  d_x[i] <- centers[combinations[i,1],1] - (1/sqrt(u[i]))
  d_y[i] <- centers[combinations[i,1],2] - (R[i]/sqrt(u[i]))
  e_x[i] <- centers[combinations[i,1],1] + (1/sqrt(u[i]))
  e_y[i] <- centers[combinations[i,1],2] + (R[i]/sqrt(u[i]))
  D_x[i] <- centers[combinations[i,2],1] - (1/sqrt(U[i]))
  D_y[i] <- centers[combinations[i,2],2] - (R[i]/sqrt(U[i]))
  E_x[i] <- centers[combinations[i,2],1] + (1/sqrt(U[i]))
  E_y[i] <- centers[combinations[i,2],2] + (R[i]/sqrt(U[i]))
  #ellips_dist_1[i] <- sqrt((E_y[i] - e_y[i])^2 + (E_x[i] - e_x[i])^2)
  #ellips_dist_2[i] <- sqrt((d_y[i] - D_y[i])^2 + (d_x[i] - D_x[i])^2)
  ellips_dist_1[i] <- sqrt((E_y[i] - e_y[i])^2 + (E_x[i] - e_x[i])^2)
  ellips_dist_2[i] <- sqrt((d_y[i] - D_y[i])^2 + (d_x[i] - D_x[i])^2)
}  
new_distance_radius_df <- as.data.frame(matrix(NA, nrow = nrow(combinations), ncol = 2))
for(i in 1:nrow(combinations)){
  for(j in 1:2){
    new_distance_radius_df[i,j] <- rownames(Englishtowns)[combinations[i,j]]
  }
}
new_distance_radius_df <- cbind(new_distance_radius_df, ellips_dist_1, ellips_dist_2)
colnames(new_distance_radius_df) <- c("From", "To", "Row distance", "Column distance")

dist_ellips_mat <- matrix(0, nrow = 8, ncol = 8)
rownames(dist_ellips_mat) <- colnames(dist_ellips_mat) <- rownames(Englishtowns)
dist_ellips_mat[1,2] <- new_distance_radius_df[1,3]
dist_ellips_mat[2,1] <- new_distance_radius_df[1,4]
dist_ellips_mat[1,3] <- new_distance_radius_df[2,3]
dist_ellips_mat[3,1] <- new_distance_radius_df[2,4]
dist_ellips_mat[1,4] <- new_distance_radius_df[3,3]
dist_ellips_mat[4,1] <- new_distance_radius_df[3,4]
dist_ellips_mat[1,5] <- new_distance_radius_df[4,3]
dist_ellips_mat[5,1] <- new_distance_radius_df[4,4]
dist_ellips_mat[1,6] <- new_distance_radius_df[5,3]
dist_ellips_mat[6,1] <- new_distance_radius_df[5,4]
dist_ellips_mat[1,7] <- new_distance_radius_df[6,3]
dist_ellips_mat[7,1] <- new_distance_radius_df[6,4]
dist_ellips_mat[1,8] <- new_distance_radius_df[7,3]
dist_ellips_mat[8,1] <- new_distance_radius_df[7,4]
dist_ellips_mat[2,3] <- new_distance_radius_df[8,3]
dist_ellips_mat[3,2] <- new_distance_radius_df[8,4]
dist_ellips_mat[2,4] <- new_distance_radius_df[9,3]
dist_ellips_mat[4,2] <- new_distance_radius_df[9,4]
dist_ellips_mat[2,5] <- new_distance_radius_df[10,3]
dist_ellips_mat[5,2] <- new_distance_radius_df[10,4]
dist_ellips_mat[2,6] <- new_distance_radius_df[11,3]
dist_ellips_mat[6,2] <- new_distance_radius_df[11,4]
dist_ellips_mat[2,7] <- new_distance_radius_df[12,3]
dist_ellips_mat[7,2] <- new_distance_radius_df[12,4]
dist_ellips_mat[2,8] <- new_distance_radius_df[13,3]
dist_ellips_mat[8,2] <- new_distance_radius_df[13,4]
dist_ellips_mat[3,4] <- new_distance_radius_df[14,3]
dist_ellips_mat[4,3] <- new_distance_radius_df[14,4]
dist_ellips_mat[3,5] <- new_distance_radius_df[15,3]
dist_ellips_mat[5,3] <- new_distance_radius_df[15,4]
dist_ellips_mat[3,6] <- new_distance_radius_df[16,3]
dist_ellips_mat[6,3] <- new_distance_radius_df[16,4]
dist_ellips_mat[3,7] <- new_distance_radius_df[17,3]
dist_ellips_mat[7,3] <- new_distance_radius_df[17,4]
dist_ellips_mat[3,8] <- new_distance_radius_df[18,3]
dist_ellips_mat[8,3] <- new_distance_radius_df[18,4]
dist_ellips_mat[4,5] <- new_distance_radius_df[19,3]
dist_ellips_mat[5,4] <- new_distance_radius_df[19,4]
dist_ellips_mat[4,6] <- new_distance_radius_df[20,3]
dist_ellips_mat[6,4] <- new_distance_radius_df[20,4]
dist_ellips_mat[4,7] <- new_distance_radius_df[21,3]
dist_ellips_mat[7,4] <- new_distance_radius_df[21,4]
dist_ellips_mat[4,8] <- new_distance_radius_df[22,3]
dist_ellips_mat[8,4] <- new_distance_radius_df[22,4]
dist_ellips_mat[5,6] <- new_distance_radius_df[23,3]
dist_ellips_mat[6,5] <- new_distance_radius_df[23,4]
dist_ellips_mat[5,7] <- new_distance_radius_df[24,3]
dist_ellips_mat[7,5] <- new_distance_radius_df[24,4]
dist_ellips_mat[5,8] <- new_distance_radius_df[25,3]
dist_ellips_mat[8,5] <- new_distance_radius_df[25,4]
dist_ellips_mat[6,7] <- new_distance_radius_df[26,3]
dist_ellips_mat[7,6] <- new_distance_radius_df[26,4]
dist_ellips_mat[6,8] <- new_distance_radius_df[27,3]
dist_ellips_mat[8,6] <- new_distance_radius_df[27,4]
dist_ellips_mat[7,8] <- new_distance_radius_df[28,3]
dist_ellips_mat[8,7] <- new_distance_radius_df[28,4]

### Stress
Englishtown_distances <- matrix(Englishtowns, byrow = T, ncol = 1)
Englishtown_distances <- Englishtown_distances[-which(Englishtown_distances == 0)]
Ellips_distances <- matrix(dist_ellips_mat, byrow = T, ncol = 1)
Ellips_distances <- Ellips_distances[-which(Ellips_distances == 0)]

sqrt((sum((Englishtown_distances - Ellips_distances)^2)/sum((Ellips_distances)^2)))

### Plot zonder RDM
# Bij kendal -> york is de afstand d-> D
plot(c(-160,-80), c(0,110), type="n", main="Ellipses", xlab = "Dim 1", ylab = "Dim 2", asp = 1)
text(final_config[1,1],final_config[1,2],"A", cex=1.5, pos=3,col="red")
text(final_config[8,1],final_config[8,2],"B", cex=1.5, pos=3,col="red")
draw.ellipse(x = result_kendal$center[1], y =result_kendal$center[2], a = result_kendal$a, b = result_kendal$b, angle = result_kendal$theta, border = "blue")
draw.ellipse(x = result_york$center[1], y =result_york$center[2], a = result_york$a, b = result_york$b, angle = result_york$theta, border = "blue")
points(x = result_kendal$center[1], y =result_kendal$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_york$center[1], y =result_york$center[2], pch=1, col = "black", cex = 0.5)
arrows(d_x[7], d_y[7], D_x[7], D_y[7], col = "black")
segments(d_x[7], d_y[7], e_x[7], e_y[7], lty = "dashed", col = "black")
text(e_x[7], e_y[7],"e", cex=1, pos=2,col="red")
text(d_x[7], d_y[7],"d", cex=1, pos=2,col="red")
text(E_x[7], E_y[7],"E", cex=1, pos=2,col="red")
text(D_x[7], D_y[7],"D", cex=1, pos=2,col="red")

# Bij kendal -> manchester is het van e naar E
plot(c(-160,-80), c(0,110), type="n", main="Ellipses", xlab = "Dim 1", ylab = "Dim 2", asp = 1)
text(final_config[1,1],final_config[1,2],"A", cex=1.5, pos=3,col="red")
text(final_config[2,1],final_config[2,2],"B", cex=1.5, pos=3,col="red")
draw.ellipse(x = result_kendal$center[1], y =result_kendal$center[2], a = result_kendal$a, b = result_kendal$b, angle = result_kendal$theta, border = "blue")
draw.ellipse(x = result_manchester$center[1], y =result_manchester$center[2], a = result_manchester$a, b = result_manchester$b, angle = result_manchester$theta, border = "blue")
points(x = result_kendal$center[1], y =result_kendal$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_manchester$center[1], y =result_manchester$center[2], pch=1, col = "black", cex = 0.5)
arrows(e_x[1], e_y[1], E_x[1], E_y[1], col = "black")


plot(c(-200,000), c(-200,100), type="n", main="Ellipses", xlab = "Dim 1", ylab = "Dim 2", asp = 1)
points(x = d_x, y =d_y, pch=1, col = "black", cex = 0.5)
points(x = e_x, y = e_y, pch=1, col = "black", cex = 0.5)
points(x = D_x, y =D_y, pch=1, col = "black", cex = 0.5)
points(x = E_x, y = E_y, pch=1, col = "black", cex = 0.5)
text(final_config[2,1],final_config[2,2],"Manchester", cex=1, pos=3,col="red")
text(final_config[3,1],final_config[3,2],"Norwich", cex=1, pos=3,col="red")
text(d_x, d_y,"d", cex=1, pos=3,col="red")
text(e_x, e_y,"e", cex=1, pos=3,col="red")
text(D_x, D_y, "D", cex=1, pos=3,col="red")
text(E_x, E_y, "E", cex=1, pos=3,col="red")
draw.ellipse(x = result_manchester$center[1], y =result_manchester$center[2], a = result_manchester$a, b = result_manchester$b, angle = result_manchester$theta, border = "blue")
draw.ellipse(x = result_norwich$center[1], y =result_norwich$center[2], a = result_norwich$a, b = result_norwich$b, angle = result_norwich$theta, border = "blue")
points(x = result_norwich$center[1], y =result_norwich$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_norwich$center[1], y =result_norwich$center[2], pch=1, col = "black", cex = 0.5)
arrows(e_x, e_y, E_x, E_y)
segments(result_manchester$center[1], result_manchester$center[2], e_x, e_y, lty = "dashed")


radii <- c(13.875000000000000, 12.125000000000000, 0.000000000000000, 8.375000000000000, 22.875000000000000, 8.625000000000000, 14.125000000000000, 8.000000000000000) 

combinations <- expand.grid(1:8, 1:8)
combinations <- combinations[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
colnames(combinations) <- c("From", "To")
dist_radius <- numeric(nrow(combinations))
for(i in 1:nrow(combinations)){
  dist_radius[i] <- abs((sqrt((final_config[combinations[i,2],1] - final_config[combinations[i,1],1])^2 + (final_config[combinations[i,2],2] - final_config[combinations[i,1],2])^2)) + radii[combinations[i,2]] - radii[combinations[i,1]])
}
new_distance_radius_df <- as.data.frame(matrix(NA, nrow = nrow(combinations), ncol = 2))
for(i in 1:nrow(combinations)){
  for(j in 1:2){
    new_distance_radius_df[i,j] <- rownames(Englishtowns)[combinations[i,j]]
  }
}
new_distance_radius_df <- cbind(new_distance_radius_df, dist_radius)
colnames(new_distance_radius_df) <- c("From", "To", "Radius distance")

Englishtown_distances <- matrix(Englishtowns, byrow = T, ncol = 1)
Englishtown_distances <- Englishtown_distances[-which(Englishtown_distances == 0)]
sqrt((sum((Englishtown_distances - new_distance_radius_df[,3])^2)/sum((new_distance_radius_df[,3])^2)))



R <- ((centers[1,2] - centers[2,2])/(centers[1,1] - centers[2,1]))
u <- (((cos(deg2rad(rotation[1,1])) + R*sin(deg2rad(rotation[1,1])))^2)/((a_b[1,1])^2)) + (((cos(deg2rad(rotation[1,1])) - R*sin(deg2rad(rotation[1,1])))^2)/((a_b[1,2])^2))
U <- (((cos(deg2rad(rotation[2,1])) + R*sin(deg2rad(rotation[2,1])))^2)/((a_b[2,1])^2)) + (((cos(deg2rad(rotation[2,1])) - R*sin(deg2rad(rotation[2,1])))^2)/((a_b[2,2])^2))
d_x <- centers[1,1] - (1/sqrt(u))
d_y <- centers[1,2] - (R/sqrt(u))
e_x <- centers[1,1] + (1/sqrt(u))
e_y <- centers[1,2] + (R/sqrt(u))
D_x <- centers[2,1] - (1/sqrt(U))
D_y <- centers[2,2] - (R/sqrt(U))
E_x <- centers[2,1] + (1/sqrt(U))
E_y <- centers[2,2] + (R/sqrt(U))

plot(c(-160,-80), c(30,110), type="n", main="Ellipses", xlab = "Dim 1", ylab = "Dim 2", asp = 1)
points(x = d_x, y =d_y, pch=1, col = "black", cex = 0.5)
points(x = e_x, y = e_y, pch=1, col = "black", cex = 0.5)
points(x = D_x, y =D_y, pch=1, col = "black", cex = 0.5)
points(x = E_x, y = E_y, pch=1, col = "black", cex = 0.5)
text(final_config[1,1],final_config[1,2],"Kendal", cex=1, pos=3,col="red")
text(final_config[2,1],final_config[2,2],"Manchester", cex=1, pos=3,col="red")
text(d_x, d_y,"d", cex=1, pos=3,col="red")
text(e_x, e_y,"e", cex=1, pos=3,col="red")
text(D_x, D_y, "D", cex=1, pos=3,col="red")
text(E_x, E_y, "E", cex=1, pos=3,col="red")
draw.ellipse(x = result_kendal$center[1], y =result_kendal$center[2], a = result_kendal$a, b = result_kendal$b, angle = result_kendal$theta, border = "blue")
draw.ellipse(x = result_manchester$center[1], y =result_manchester$center[2], a = result_manchester$a, b = result_manchester$b, angle = result_manchester$theta, border = "blue")
points(x = result_kendal$center[1], y =result_kendal$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_manchester$center[1], y =result_manchester$center[2], pch=1, col = "black", cex = 0.5)
arrows(e_x, e_y, E_x, E_y)
segments(result_kendal$center[1], result_kendal$center[2], e_x, e_y, lty = "dashed")

# Kendal manchester
sqrt((E_y - e_y)^2 + (E_x - e_x)^2)
# Manchester kendal
sqrt((d_y - D_y)^2 + (d_x - D_x)^2)


# Manchester - Norwich en Norwhich manchester
R <- ((centers[2,2] - centers[3,2])/(centers[2,1] - centers[3,1]))
u <- (((cos(deg2rad(rotation[2,1])) + R*sin(deg2rad(rotation[2,1])))^2)/((a_b[2,1])^2)) + (((cos(deg2rad(rotation[2,1])) - R*sin(deg2rad(rotation[2,1])))^2)/((a_b[2,2])^2))
U <- (((cos(deg2rad(rotation[3,1])) + R*sin(deg2rad(rotation[3,1])))^2)/((a_b[3,1])^2)) + (((cos(deg2rad(rotation[3,1])) - R*sin(deg2rad(rotation[3,1])))^2)/((a_b[3,2])^2))
d_x <- centers[2,1] - (1/sqrt(u))
d_y <- centers[2,2] - (R/sqrt(u))
e_x <- centers[2,1] + (1/sqrt(u))
e_y <- centers[2,2] + (R/sqrt(u))
D_x <- centers[3,1] - (1/sqrt(U))
D_y <- centers[3,2] - (R/sqrt(U))
E_x <- centers[3,1] + (1/sqrt(U))
E_y <- centers[3,2] + (R/sqrt(U))

plot(c(-200,000), c(-200,100), type="n", main="Ellipses", xlab = "Dim 1", ylab = "Dim 2", asp = 1)
points(x = d_x, y =d_y, pch=1, col = "black", cex = 0.5)
points(x = e_x, y = e_y, pch=1, col = "black", cex = 0.5)
points(x = D_x, y =D_y, pch=1, col = "black", cex = 0.5)
points(x = E_x, y = E_y, pch=1, col = "black", cex = 0.5)
text(final_config[2,1],final_config[2,2],"Manchester", cex=1, pos=3,col="red")
text(final_config[3,1],final_config[3,2],"Norwich", cex=1, pos=3,col="red")
text(d_x, d_y,"d", cex=1, pos=3,col="red")
text(e_x, e_y,"e", cex=1, pos=3,col="red")
text(D_x, D_y, "D", cex=1, pos=3,col="red")
text(E_x, E_y, "E", cex=1, pos=3,col="red")
draw.ellipse(x = result_manchester$center[1], y =result_manchester$center[2], a = result_manchester$a, b = result_manchester$b, angle = result_manchester$theta, border = "blue")
draw.ellipse(x = result_norwich$center[1], y =result_norwich$center[2], a = result_norwich$a, b = result_norwich$b, angle = result_norwich$theta, border = "blue")
points(x = result_norwich$center[1], y =result_norwich$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_norwich$center[1], y =result_norwich$center[2], pch=1, col = "black", cex = 0.5)
arrows(e_x, e_y, E_x, E_y)
segments(result_manchester$center[1], result_manchester$center[2], e_x, e_y, lty = "dashed")

#  Manchester Norwich
sqrt((E_y - e_y)^2 + (E_x - e_x)^2)
# Norwich Manchester
sqrt((d_y - D_y)^2 + (d_x - D_x)^2)

# Slide vector
slide_vec <- slidevector(Englishtowns)
plot(slide_vec)
skewsymmetry(Englishtowns)
summary(slide_vec)

-8.1229 -4.9596


### Plot zonder RDM
plot(c(-160,-80), c(30,110), type="n", main="Ellipses", xlab = "Dim 1", ylab = "Dim 2", asp = 1)
text(final_config[1,1],final_config[1,2],"A", cex=1.5, pos=3,col="red")
text(final_config[2,1],final_config[2,2],"B", cex=1.5, pos=3,col="red")
draw.ellipse(x = result_kendal$center[1], y =result_kendal$center[2], a = result_kendal$a, b = result_kendal$b, angle = result_kendal$theta, border = "blue")
draw.ellipse(x = result_manchester$center[1], y =result_manchester$center[2], a = result_manchester$a, b = result_manchester$b, angle = result_manchester$theta, border = "blue")
points(x = result_kendal$center[1], y =result_kendal$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_manchester$center[1], y =result_manchester$center[2], pch=1, col = "black", cex = 0.5)
arrows(e_x[1], e_y[1], E_x[1], E_y[1], col = "black")
segments(d_x[1], d_y[1], e_x[1], e_y[1], lty = "dashed", col = "black")
text(e_x[1], e_y[1],"e", cex=1, pos=3,col="red")
text(d_x[1], d_y[1],"d", cex=1, pos=3,col="red")
text(E_x[1], E_y[1],"E", cex=1, pos=3,col="red")
text(D_x[1], D_y[1],"D", cex=1, pos=3,col="red")

?segments
### Plot zonder RDM
plot(c(-300,300), c(-300,300), type="n", main="Ellipses", xlab = "Dim 1", ylab = "Dim 2", asp = 1)
text(final_config[1,1],final_config[1,2],row.names(Englishtowns)[1], cex=0.5, pos=3,col="red")
text(final_config[2,1],final_config[2,2],row.names(Englishtowns)[2], cex=0.5, pos=3,col="red")
text(final_config[3,1],final_config[3,2],row.names(Englishtowns)[3], cex=0.5, pos=3,col="red")
text(final_config[4,1],final_config[4,2],row.names(Englishtowns)[4], cex=0.5, pos=3,col="red")
text(final_config[5,1],final_config[5,2],row.names(Englishtowns)[5], cex=0.5, pos=3,col="red")
text(final_config[6,1],final_config[6,2],row.names(Englishtowns)[6], cex=0.5, pos=3,col="red")
text(final_config[7,1],final_config[7,2],row.names(Englishtowns)[7], cex=0.5, pos=3,col="red")
text(final_config[8,1],final_config[8,2],row.names(Englishtowns)[8], cex=0.5, pos=3,col="red")
draw.ellipse(x = result_kendal$center[1], y =result_kendal$center[2], a = result_kendal$a, b = result_kendal$b, angle = result_kendal$theta, border = "blue")
draw.ellipse(x = result_manchester$center[1], y =result_manchester$center[2], a = result_manchester$a, b = result_manchester$b, angle = result_manchester$theta, border = "blue")
draw.ellipse(x = result_norwich$center[1], y =result_norwich$center[2], a = result_norwich$a, b = result_norwich$b, angle = result_norwich$theta, border = "blue")
draw.ellipse(x = result_oxford$center[1], y =result_oxford$center[2], a = result_oxford$a, b = result_oxford$b, angle = result_oxford$theta, border = "blue")
draw.ellipse(x = result_penzance$center[1], y =result_penzance$center[2], a = result_penzance$a, b = result_penzance$b, angle = result_penzance$theta, border = "blue")
draw.ellipse(x = result_southampton$center[1], y =result_southampton$center[2], a = result_southampton$a, b = result_southampton$b, angle = result_southampton$theta, border = "blue")
draw.ellipse(x = result_taunton$center[1], y =result_taunton$center[2], a = result_taunton$a, b = result_taunton$b, angle = result_taunton$theta, border = "blue")
draw.ellipse(x = result_york$center[1], y =result_york$center[2], a = result_york$a, b = result_york$b, angle = result_york$theta, border = "blue")
segments(final_config[1,1], final_config[1,2], vectors[1,1], vectors[1,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[2,1], vectors[2,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[3,1], vectors[3,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[4,1], vectors[4,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[5,1], vectors[5,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[6,1], vectors[6,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[7,1], vectors[7,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[8,1], vectors[8,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[9,1], vectors[9,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[10,1], vectors[10,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[11,1], vectors[11,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[12,1], vectors[12,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[13,1], vectors[13,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[14,1], vectors[14,2], col = "green")
segments(final_config[3,1], final_config[3,2], vectors[15,1], vectors[15,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[16,1], vectors[16,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[17,1], vectors[17,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[18,1], vectors[18,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[19,1], vectors[19,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[20,1], vectors[20,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[21,1], vectors[21,2], col = "green")
segments(final_config[4,1], final_config[4,2],  vectors[22,1], vectors[22,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[23,1], vectors[23,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[24,1], vectors[24,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[25,1], vectors[25,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[26,1], vectors[26,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[27,1], vectors[27,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[28,1], vectors[28,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[29,1], vectors[29,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[30,1], vectors[30,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[31,1], vectors[31,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[32,1], vectors[32,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[33,1], vectors[33,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[34,1], vectors[34,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[35,1], vectors[35,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[36,1], vectors[36,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[37,1], vectors[37,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[38,1], vectors[38,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[39,1], vectors[39,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[40,1], vectors[40,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[41,1], vectors[41,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[42,1], vectors[42,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[43,1], vectors[43,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[44,1], vectors[44,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[45,1], vectors[45,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[46,1], vectors[46,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[47,1], vectors[47,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[48,1], vectors[48,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[49,1], vectors[49,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[50,1], vectors[50,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[51,1], vectors[51,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[52,1], vectors[52,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[53,1], vectors[53,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[54,1], vectors[54,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[55,1], vectors[55,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[56,1], vectors[56,2], col = "green")
points(x = result_kendal$center[1], y =result_kendal$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_manchester$center[1], y =result_manchester$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_norwich$center[1], y =result_norwich$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_oxford$center[1], y =result_oxford$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_penzance$center[1], y =result_penzance$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_southampton$center[1], y =result_southampton$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_taunton$center[1], y =result_taunton$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_york$center[1], y =result_york$center[2], pch=1, col = "black", cex = 0.5)

### Plot
plot(c(-300,300), c(-300,300), type="n", main="Ellipses", xlab = "Dim 1", ylab = "Dim 2", asp = 1)
text(final_config[1,1],final_config[1,2],row.names(Englishtowns)[1], cex=0.5, pos=3,col="red")
text(final_config[2,1],final_config[2,2],row.names(Englishtowns)[2], cex=0.5, pos=3,col="red")
text(final_config[3,1],final_config[3,2],row.names(Englishtowns)[3], cex=0.5, pos=3,col="red")
text(final_config[4,1],final_config[4,2],row.names(Englishtowns)[4], cex=0.5, pos=3,col="red")
text(final_config[5,1],final_config[5,2],row.names(Englishtowns)[5], cex=0.5, pos=3,col="red")
text(final_config[6,1],final_config[6,2],row.names(Englishtowns)[6], cex=0.5, pos=3,col="red")
text(final_config[7,1],final_config[7,2],row.names(Englishtowns)[7], cex=0.5, pos=3,col="red")
text(final_config[8,1],final_config[8,2],row.names(Englishtowns)[8], cex=0.5, pos=3,col="red")
draw.ellipse(x = result_kendal$center[1], y =result_kendal$center[2], a = result_kendal$a, b = result_kendal$b, angle = result_kendal$theta, border = "blue")
draw.ellipse(x = result_manchester$center[1], y =result_manchester$center[2], a = result_manchester$a, b = result_manchester$b, angle = result_manchester$theta, border = "blue")
draw.ellipse(x = result_norwich$center[1], y =result_norwich$center[2], a = result_norwich$a, b = result_norwich$b, angle = result_norwich$theta, border = "blue")
draw.ellipse(x = result_oxford$center[1], y =result_oxford$center[2], a = result_oxford$a, b = result_oxford$b, angle = result_oxford$theta, border = "blue")
draw.ellipse(x = result_penzance$center[1], y =result_penzance$center[2], a = result_penzance$a, b = result_penzance$b, angle = result_penzance$theta, border = "blue")
draw.ellipse(x = result_southampton$center[1], y =result_southampton$center[2], a = result_southampton$a, b = result_southampton$b, angle = result_southampton$theta, border = "blue")
draw.ellipse(x = result_taunton$center[1], y =result_taunton$center[2], a = result_taunton$a, b = result_taunton$b, angle = result_taunton$theta, border = "blue")
draw.ellipse(x = result_york$center[1], y =result_york$center[2], a = result_york$a, b = result_york$b, angle = result_york$theta, border = "blue")
for(i in 1:8){
  draw.circle(x = final_config[i,1], y = final_config[i,2], radii[i], border = "black")
}
segments(final_config[1,1], final_config[1,2], vectors[1,1], vectors[1,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[2,1], vectors[2,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[3,1], vectors[3,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[4,1], vectors[4,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[5,1], vectors[5,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[6,1], vectors[6,2], col = "green")
segments(final_config[1,1], final_config[1,2], vectors[7,1], vectors[7,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[8,1], vectors[8,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[9,1], vectors[9,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[10,1], vectors[10,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[11,1], vectors[11,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[12,1], vectors[12,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[13,1], vectors[13,2], col = "green")
segments(final_config[2,1], final_config[2,2],  vectors[14,1], vectors[14,2], col = "green")
segments(final_config[3,1], final_config[3,2], vectors[15,1], vectors[15,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[16,1], vectors[16,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[17,1], vectors[17,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[18,1], vectors[18,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[19,1], vectors[19,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[20,1], vectors[20,2], col = "green")
segments(final_config[3,1], final_config[3,2],  vectors[21,1], vectors[21,2], col = "green")
segments(final_config[4,1], final_config[4,2],  vectors[22,1], vectors[22,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[23,1], vectors[23,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[24,1], vectors[24,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[25,1], vectors[25,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[26,1], vectors[26,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[27,1], vectors[27,2], col = "green")
segments(final_config[4,1], final_config[4,2], vectors[28,1], vectors[28,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[29,1], vectors[29,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[30,1], vectors[30,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[31,1], vectors[31,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[32,1], vectors[32,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[33,1], vectors[33,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[34,1], vectors[34,2], col = "green")
segments(final_config[5,1], final_config[5,2], vectors[35,1], vectors[35,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[36,1], vectors[36,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[37,1], vectors[37,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[38,1], vectors[38,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[39,1], vectors[39,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[40,1], vectors[40,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[41,1], vectors[41,2], col = "green")
segments(final_config[6,1], final_config[6,2], vectors[42,1], vectors[42,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[43,1], vectors[43,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[44,1], vectors[44,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[45,1], vectors[45,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[46,1], vectors[46,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[47,1], vectors[47,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[48,1], vectors[48,2], col = "green")
segments(final_config[7,1], final_config[7,2], vectors[49,1], vectors[49,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[50,1], vectors[50,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[51,1], vectors[51,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[52,1], vectors[52,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[53,1], vectors[53,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[54,1], vectors[54,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[55,1], vectors[55,2], col = "green")
segments(final_config[8,1], final_config[8,2], vectors[56,1], vectors[56,2], col = "green")
points(x = result_kendal$center[1], y =result_kendal$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_manchester$center[1], y =result_manchester$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_norwich$center[1], y =result_norwich$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_oxford$center[1], y =result_oxford$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_penzance$center[1], y =result_penzance$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_southampton$center[1], y =result_southampton$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_taunton$center[1], y =result_taunton$center[2], pch=1, col = "black", cex = 0.5)
points(x = result_york$center[1], y =result_york$center[2], pch=1, col = "black", cex = 0.5)

