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
a_b <- ellipse_a_b/2

drift_end <- as.data.frame(rbind(c(-144.495919939308, 91.436733536520),
                                 c(-91.973704113374, 43.078391284222),
                                 c(-99.234374053059, -159.809724844431),
                                 c(11.840497370732, -46.249609811980),
                                 c(226.348375478697,44.396132505969),
                                 c(73.836615365397,-70.670991656576),
                                 c(101.506191259054,10.498897775017),
                                 c(-152.605532124492,13.185188015017)))

# Iterative algorithm results
combinations <- expand.grid(1:8, 1:8)
combinations <- combinations[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
colnames(combinations) <- c("From", "To")
ptp_dist <- numeric(nrow(combinations))
ctc_dist <- numeric(nrow(combinations))
ctp_dist <- numeric(nrow(combinations))
ptc_dist <- numeric(nrow(combinations))
etc_dist <- numeric(nrow(combinations))
cte_dist <- numeric(nrow(combinations))
etp_dist <- numeric(nrow(combinations))
pte_dist <- numeric(nrow(combinations))
for(i in 1:nrow(combinations)){
  ptp_dist[i] <- (sqrt((final_config[combinations[i,2],1] - final_config[combinations[i,1],1])^2 + (final_config[combinations[i,2],2] - final_config[combinations[i,1],2])^2))
  ctc_dist[i] <- (sqrt((centers[combinations[i,2],1] - centers[combinations[i,1],1])^2 + (centers[combinations[i,2],2] - centers[combinations[i,1],2])^2))
  ctp_dist[i] <- (sqrt((final_config[combinations[i,2],1] - centers[combinations[i,1],1])^2 + (final_config[combinations[i,2],2] - centers[combinations[i,1],2])^2))
  ptc_dist[i] <- (sqrt((centers[combinations[i,2],1] - final_config[combinations[i,1],1])^2 + (centers[combinations[i,2],2] - final_config[combinations[i,1],2])^2))
  etc_dist[i] <- (sqrt((centers[combinations[i,2],1] - drift_end[combinations[i,1],1])^2 + (centers[combinations[i,2],2] - drift_end[combinations[i,1],2])^2))
  cte_dist[i] <- (sqrt((drift_end[combinations[i,2],1] - centers[combinations[i,1],1])^2 + (drift_end[combinations[i,2],2] - centers[combinations[i,1],2])^2))
  etp_dist[i] <- (sqrt((final_config[combinations[i,2],1] - drift_end[combinations[i,1],1])^2 + (final_config[combinations[i,2],2] - drift_end[combinations[i,1],2])^2))
  pte_dist[i] <- (sqrt((drift_end[combinations[i,2],1] - final_config[combinations[i,1],1])^2 + (drift_end[combinations[i,2],2] - final_config[combinations[i,1],2])^2))
  }
all_dist_df <- as.data.frame(matrix(NA, nrow = nrow(combinations), ncol = 2))
for(i in 1:nrow(combinations)){
  for(j in 1:2){
    all_dist_df[i,j] <- rownames(Englishtowns)[combinations[i,j]]
  }
}
all_dist_df <- cbind(all_dist_df, ptp_dist, ctc_dist, ctp_dist, ptc_dist, etc_dist, cte_dist, etp_dist, pte_dist)
colnames(all_dist_df) <- c("From", "To", "Object to object distance", "Centroid to centroid distance", "Centroid to object distance", "Object to centroid distance", "End point to center", "Center to end point", "End point to object", "Object to end point")

# Centroids obtainer as a mean from drift vectors
vec_kendal <- colMeans(vectors[1:7,])
vec_manchester <- colMeans(vectors[8:14,])
vec_norwich <- colMeans(vectors[15:21,])
vec_oxford <- colMeans(vectors[22:28,])
vec_penzance <- colMeans(vectors[29:35,])
vec_southampton <- colMeans(vectors[36:42,])
vec_taunton <- colMeans(vectors[43:49,])
vec_york <- colMeans(vectors[50:56,])

vectorcentroids <- as.data.frame(matrix(rbind(vec_kendal, vec_manchester, vec_norwich, vec_oxford, vec_penzance, vec_southampton, vec_taunton, vec_york),ncol = 2, byrow = F))

combinations <- expand.grid(1:8, 1:8)
combinations <- combinations[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
colnames(combinations) <- c("From", "To")
ptp_dist <- numeric(nrow(combinations))
ctc_dist <- numeric(nrow(combinations))
ctp_dist <- numeric(nrow(combinations))
ptc_dist <- numeric(nrow(combinations))
etc_dist <- numeric(nrow(combinations))
cte_dist <- numeric(nrow(combinations))
etp_dist <- numeric(nrow(combinations))
pte_dist <- numeric(nrow(combinations))
for(i in 1:nrow(combinations)){
  ptp_dist[i] <- (sqrt((final_config[combinations[i,2],1] - final_config[combinations[i,1],1])^2 + (final_config[combinations[i,2],2] - final_config[combinations[i,1],2])^2))
  ctc_dist[i] <- (sqrt((vectorcentroids[combinations[i,2],1] - vectorcentroids[combinations[i,1],1])^2 + (vectorcentroids[combinations[i,2],2] - vectorcentroids[combinations[i,1],2])^2))
  # maybe scale by half 
  ctp_dist[i] <- (sqrt((final_config[combinations[i,2],1] - vectorcentroids[combinations[i,1],1])^2 + (final_config[combinations[i,2],2] - vectorcentroids[combinations[i,1],2])^2))
  ptc_dist[i] <- (sqrt((vectorcentroids[combinations[i,2],1] - final_config[combinations[i,1],1])^2 + (vectorcentroids[combinations[i,2],2] - final_config[combinations[i,1],2])^2))
  etc_dist[i] <- (sqrt((vectorcentroids[combinations[i,2],1] - drift_end[combinations[i,1],1])^2 + (vectorcentroids[combinations[i,2],2] - drift_end[combinations[i,1],2])^2))
  cte_dist[i] <- (sqrt((drift_end[combinations[i,2],1] - vectorcentroids[combinations[i,1],1])^2 + (drift_end[combinations[i,2],2] - vectorcentroids[combinations[i,1],2])^2))
  etp_dist[i] <- (sqrt((final_config[combinations[i,2],1] - drift_end[combinations[i,1],1])^2 + (final_config[combinations[i,2],2] - drift_end[combinations[i,1],2])^2))
  pte_dist[i] <- (sqrt((drift_end[combinations[i,2],1] - final_config[combinations[i,1],1])^2 + (drift_end[combinations[i,2],2] - final_config[combinations[i,1],2])^2))
}
vector_centroid_dist_df <- as.data.frame(matrix(NA, nrow = nrow(combinations), ncol = 2))
for(i in 1:nrow(combinations)){
  for(j in 1:2){
    vector_centroid_dist_df[i,j] <- rownames(Englishtowns)[combinations[i,j]]
  }
}
vector_centroid_dist_df <- cbind(vector_centroid_dist_df, ptp_dist, ctc_dist, ctp_dist, ptc_dist, etc_dist, cte_dist, etp_dist, pte_dist)
colnames(vector_centroid_dist_df) <- c("From", "To", "Point to point distance", "Centroid to centroid distance", "Centroid to point distance", "Point to centroid distance", "End point to center", "Center to end point", "End point to object", "Object to end point")


# Stress = distances(original data) - lower dimensional distances
Englishtown_distances <- matrix(Englishtowns, byrow = T, ncol = 1)
Englishtown_distances <- Englishtown_distances[-which(Englishtown_distances == 0)]

sqrt((sum((Englishtown_distances - all_dist_df[,3])^2)/sum((all_dist_df[,3])^2)))
sqrt((sum((Englishtown_distances - all_dist_df[,4])^2)/sum((all_dist_df[,4])^2)))
sqrt((sum((Englishtown_distances - all_dist_df[,5])^2)/sum((all_dist_df[,5])^2)))
sqrt((sum((Englishtown_distances - all_dist_df[,6])^2)/sum((all_dist_df[,6])^2)))
sqrt((sum((Englishtown_distances - all_dist_df[,7])^2)/sum((all_dist_df[,7])^2)))
sqrt((sum((Englishtown_distances - all_dist_df[,8])^2)/sum((all_dist_df[,7])^2)))
sqrt((sum((Englishtown_distances - all_dist_df[,9])^2)/sum((all_dist_df[,9])^2)))
sqrt((sum((Englishtown_distances - all_dist_df[,10])^2)/sum((all_dist_df[,10])^2)))

# Stress for center
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,3])^2)/sum((vector_centroid_dist_df[,3])^2)))
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,4])^2)/sum((vector_centroid_dist_df[,4])^2)))
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,5])^2)/sum((vector_centroid_dist_df[,5])^2)))
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,6])^2)/sum((vector_centroid_dist_df[,6])^2)))
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,7])^2)/sum((vector_centroid_dist_df[,7])^2)))
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,8])^2)/sum((vector_centroid_dist_df[,8])^2)))
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,9])^2)/sum((vector_centroid_dist_df[,9])^2)))
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,10])^2)/sum((vector_centroid_dist_df[,10])^2)))

# Methode met ellips
combinations <- expand.grid(1:8, 1:8)
combinations <- combinations[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
colnames(combinations) <- c("From", "To")
theta <- numeric(nrow(combinations))
dist <- numeric(nrow(combinations))
from_x_new <- numeric(nrow(combinations))
from_y_new <- numeric(nrow(combinations))
d_intersect <- numeric(nrow(combinations))
theta_new <- numeric(nrow(combinations))
max_dist <- numeric(nrow(combinations))
for(i in 1:nrow(combinations)){
  # Het deel om de afstand te berekenen van het punt waar je naartoe gaat tot het snijpunt
  theta[i] <- deg2rad(360-rotation[combinations[i,2],1]) 
  from_x_new[i] = final_config[combinations[i,2],1] + (centers[combinations[i,1],1] - final_config[combinations[i,2],1])*cos(theta[i]) - (centers[combinations[i,1],2] - final_config[combinations[i,2],2])*sin(theta[i])
  from_y_new[i] = final_config[combinations[i,2],2] + (final_config[combinations[i,1],1] - final_config[combinations[i,2],1])*sin(theta[i]) + (centers[combinations[i,1],2] - final_config[combinations[i,2],2])*cos(theta[i])
  theta_new[i] <- atan2(c(from_y_new[i] - final_config[combinations[i,2],2]), c(from_x_new[i] - final_config[combinations[i,2],1]))
  d_intersect[i] <- sqrt(((a_b[combinations[i,2],1]^2)*(cos(theta_new[i])^2))+((a_b[combinations[i,2],2]^2)*(sin(theta_new[i])^2)))
  max_dist[i] <- sqrt((final_config[combinations[i,2],1] - from_x_new[i])^2 + (final_config[combinations[i,2],2] - from_y_new[i])^2) + d_intersect[i]
}
methode_1_df <- as.data.frame(matrix(NA, nrow = nrow(combinations), ncol = 2))
for(i in 1:nrow(combinations)){
  for(j in 1:2){
    methode_1_df[i,j] <- rownames(Englishtowns)[combinations[i,j]]
  }
}
methode_1_df <- cbind(methode_1_df, max_dist)
colnames(methode_1_df) <- c("From", "To", "Distance")
sqrt((sum((Englishtown_distances - methode_1_df[,3])^2)/sum((methode_1_df[,3])^2)))

360-rotation

# Alt method
M <- function(x, y, psi, to_x, to_y, r_x, r_y){
  ((((x - to_x)*cos(psi) + (y - to_y)*sin(psi))^2)/r_x^2) + ((((x - to_x)*sin(psi) - (y - to_y)*cos(psi))^2)/r_y^2)
} 
combinations <- expand.grid(1:8, 1:8)
combinations <- combinations[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
colnames(combinations) <- c("From", "To")
psi <- numeric(nrow(combinations))
s_1 <- numeric(nrow(combinations))
s_2 <- numeric(nrow(combinations))
d_qs <- numeric(nrow(combinations))
d_pq <- numeric(nrow(combinations))
d_ps <- numeric(nrow(combinations))
for(i in 1:nrow(combinations)){
  psi[i] <- rotation[combinations[i,2],1]
  s_1[i] <- final_config[combinations[i,2],1] + (1/M(centers[combinations[i,1],1], centers[combinations[i,1],2], psi[i], final_config[combinations[i,2],1], final_config[combinations[i,2],2], a_b[combinations[i,2],1], a_b[combinations[i,2],2]))*(centers[combinations[i,1],1] - centers[combinations[i,2],1])
  s_2[i] <- final_config[combinations[i,2],2] + (1/M(centers[combinations[i,1],1], centers[combinations[i,1],2], psi[i], final_config[combinations[i,2],1], final_config[combinations[i,2],2], a_b[combinations[i,2],1], a_b[combinations[i,2],2]))*(centers[combinations[i,1],2] - centers[combinations[i,2],2])
  d_qs[i] <- sqrt((s_1[i] - final_config[combinations[i,2],1])^2 + (s_2[i] - final_config[combinations[i,2],2])^2)
  d_pq[i] <- sqrt((final_config[combinations[i,2],1] - centers[combinations[i,1],1])^2 + (final_config[combinations[i,2],2] - centers[combinations[i,1],2])^2) 
  d_ps[i] <- d_pq[i] + d_qs[i]
}
methode_2_df <- as.data.frame(matrix(NA, nrow = nrow(combinations), ncol = 2))
for(i in 1:nrow(combinations)){
  for(j in 1:2){
    methode_2_df[i,j] <- rownames(Englishtowns)[combinations[i,j]]
  }
}
methode_2_df <- cbind(methode_2_df, d_ps)
colnames(methode_2_df) <- c("From", "To", "Distance")
sqrt((sum((Englishtown_distances - methode_2_df[,3])^2)/sum((methode_2_df[,3])^2)))



# Laagste stress als je de afstand van centroid naar punt pakt, en dan als centroid de centroid van het iteratieve algoritme. Maar, zelfs de centroid naar punt manier van de centroid als gemiddelde van drift vectors is beter dan de andere


###
radii <- c(13.875000000000000, 12.125000000000000, 0.000000000000000, 8.375000000000000, 22.875000000000000, 8.625000000000000, 14.125000000000000, 8.000000000000000) 
green_red_radii <- c(2.875000000000, 1.125000000000, -11.000000000000, -2.625000000000, 11.875000000000, -2.375000000000, 3.125000000000, -3.000000000000)


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

#### Het nieuwe RDM met de negatieve en positieve radii
combinations <- expand.grid(1:8, 1:8)
combinations <- combinations[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
colnames(combinations) <- c("From", "To")
dist_radius <- numeric(nrow(combinations))
for(i in 1:nrow(combinations)){
  if(green_red_radii[combinations[i,1]] >= 0){
    dist_radius[i] <- abs((sqrt((final_config[combinations[i,1],1] - final_config[combinations[i,2],1])^2 + (final_config[combinations[i,1],2] - final_config[combinations[i,2],2])^2)) + green_red_radii[combinations[i,2]] - green_red_radii[combinations[i,1]])
  } else{
    dist_radius[i] <- abs((sqrt((final_config[combinations[i,1],1] - final_config[combinations[i,2],1])^2 + (final_config[combinations[i,1],2] - final_config[combinations[i,2],2])^2)) + green_red_radii[combinations[i,2]] + green_red_radii[combinations[i,1]])
  }
}
distance_green_red_radius_df <- as.data.frame(matrix(NA, nrow = nrow(combinations), ncol = 2))
for(i in 1:nrow(combinations)){
  for(j in 1:2){
    distance_green_red_radius_df[i,j] <- rownames(Englishtowns)[combinations[i,j]]
  }
}
distance_green_red_radius_df <- cbind(distance_green_red_radius_df, dist_radius)
colnames(distance_green_red_radius_df) <- c("From", "To", "Radius distance")

sqrt((sum((Englishtown_distances - new_distance_radius_df[,3])^2)/sum((new_distance_radius_df[,3])^2)))

sqrt((sum((Englishtown_distances - distance_green_red_radius_df[,3])^2)/sum((vector_centroid_dist_df[,3])^2)))
Englishtowns



# Nog eentje met herschaling
combinations <- expand.grid(1:8, 1:8)
combinations <- combinations[-c(1, 10, 19, 28, 37 ,46, 55, 64),]
colnames(combinations) <- c("From", "To")
ctc_dist <- numeric(nrow(combinations))
ctc_cor1_dist <- numeric(nrow(combinations))
ctc_cor2_dist <- numeric(nrow(combinations))
for(i in 1:nrow(combinations)){
  ctc_dist[i] <- (sqrt((vectorcentroids[combinations[i,2],1] - vectorcentroids[combinations[i,1],1])^2 + (vectorcentroids[combinations[i,2],2] - vectorcentroids[combinations[i,1],2])^2))
  ctc_cor1_dist[i] <- (sqrt((vectorcentroids[combinations[i,2],1] - vectorcentroids[combinations[i,1],1])^2 + (vectorcentroids[combinations[i,2],2] - vectorcentroids[combinations[i,1],2])^2)) + 0.5*(sqrt((vectorcentroids[combinations[i,2],1] - final_config[combinations[i,1],1])^2 + (vectorcentroids[combinations[i,2],2] - final_config[combinations[i,1],2])^2))
  ctc_cor2_dist[i] <- (sqrt((vectorcentroids[combinations[i,2],1] - vectorcentroids[combinations[i,1],1])^2 + (vectorcentroids[combinations[i,2],2] - vectorcentroids[combinations[i,1],2])^2)) + 0.5*(sqrt((vectorcentroids[combinations[i,1],1] - final_config[combinations[i,2],1])^2 + (vectorcentroids[combinations[i,1],2] - final_config[combinations[i,2],2])^2))
}
vector_centroid_dist_df <- as.data.frame(matrix(NA, nrow = nrow(combinations), ncol = 2))
for(i in 1:nrow(combinations)){
  for(j in 1:2){
    vector_centroid_dist_df[i,j] <- rownames(Englishtowns)[combinations[i,j]]
  }
}
vector_centroid_dist_df <- cbind(vector_centroid_dist_df, ctc_dist, ctc_cor1_dist, ctc_cor2_dist)
colnames(vector_centroid_dist_df) <- c("From", "To", "Centroid to centroid distance", "Correction + 0.5 distance point i centroid i", "Correction + 0.5 distance point j centroid j")
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,4])^2)/sum((all_dist_df[,4])^2)))
sqrt((sum((Englishtown_distances - vector_centroid_dist_df[,5])^2)/sum((all_dist_df[,5])^2)))

#distance_Etown <- matrix(dist(Englishtowns), byrow = T, ncol = 1)

#m <- matrix(NA,8,8)
#m[lower.tri(m,diag=F)] <- matrix(dist(Englishtowns), byrow = T, ncol = 1)

#makeSymm <- function(m) {
#  m[upper.tri(m)] <- t(m)[upper.tri(m)]
#  return(m)
#}
#m <- makeSymm(m)
#rownames(m) <- rownames(Englishtowns)
#colnames(m) <- rownames(Englishtowns)
#m <- matrix(m, byrow = T, ncol = 1)
#m <- na.omit(m)
#m <- matrix(m)

# Stress for iterative
#sqrt((sum((m - all_dist_df[,3])^2))/sum((all_dist_df[,3])^2))
#sqrt((sum((m - all_dist_df[,4])^2))/sum((all_dist_df[,4])^2))
#sqrt((sum((m - all_dist_df[,5])^2))/sum((all_dist_df[,5])^2))
#sqrt((sum((m - all_dist_df[,6])^2))/sum((all_dist_df[,6])^2))

# Stress for center
#sqrt((sum((m - vector_centroid_dist_df[,3])^2))/sum((vector_centroid_dist_df[,3])^2))
#sqrt((sum((m - vector_centroid_dist_df[,4])^2))/sum((vector_centroid_dist_df[,4])^2))
#sqrt((sum((m - vector_centroid_dist_df[,5])^2))/sum((vector_centroid_dist_df[,5])^2))
#sqrt((sum((m - vector_centroid_dist_df[,6])^2))/sum((vector_centroid_dist_df[,6])^2))



