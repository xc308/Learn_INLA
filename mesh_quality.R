#*******#
# Mesh
#*******#

# ref: https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/

# spatially varying mesh 
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)


# for plot with ggplot2
install.packages("inlabru")
library(inlabru)


## The algorithm first builds a basic mesh, 
## including any points the user specifies, 
## as well as any boundary curves 
## (if no boundary curve is given, a default boundary will be added), 
## connecting all the points into a Delauney triangulation.

## Then mesh quality is defined by two:
  # the minmum allowed angle of tri
  # the maximum allowed tri edge length

## as long as any tri does not fullfil this criteria, 
# a new mesh point is added, n a way that is guaranteed to locally fix the problem, 
# and a new Delaunay triangulation is obtained


## the convergence of this algo is guaranteed if the 
  # the minimum angle is no larger than 21 deg
  # the max length is strictly positive


## A basic mesh with regular interior triangles is as follows:

grd <- expand.grid(1:10, 1:10) # data.frame':	100 obs. of  2 variables:
str(grd)

loc <- as.matrix(expand.grid(1:10, 1:10))
bd1 <- inla.nonconvex.hull(points = loc, convex = 1, concave = 10)
bd2 <- inla.nonconvex.hull(points = loc, convex = -1, concave = 10)
bd3 <- inla.nonconvex.hull(points = loc, convex = -1, concave = -10)


# convex:
# The desired extension radius. 
# Also determines the smallest allowed convex curvature radius. 
# Negative values are interpreted as fractions of the approximate initial set diameter.


# concave:
# The desired minimal concave curvature radius.


par(mfrow = c(1, 3))
mesh1 <- inla.mesh.create(loc = loc, boundary = list(bd1),
                         refine = list(max.edge = Inf))

mesh2 <- inla.mesh.create(loc = loc, boundary = list(bd2),
                          refine = list(max.edge = Inf))

mesh3 <- inla.mesh.create(loc = loc, boundary = list(bd3),
                          refine = list(max.edge = Inf))


plt_fun <- function(mesh) {
  ggplot() + 
    gg(mesh) +
    coord_equal()
}

plt_fun(mesh = mesh1)
plt_fun(mesh = mesh2)
plt_fun(mesh = mesh3)


plt_mesh_bd(loc = loc, boundary = bd1)
plt_mesh_bd(loc = loc, boundary = bd1, max.edge = Inf)
plt_mesh_bd(loc = loc, boundary = bd2, max.edge = Inf)


# The refine = list(max.edge = Inf) setting makes fmesher enforce the 
# default minimum angle criterion (21 degrees) 
# but ignore the edge length criterion


mesh_a <-inla.mesh.create(loc = loc, 
                 boundary = list(bd1), 
                 refine = list(max.edge = 0.5))


mesh_b <-inla.mesh.create(loc = loc, 
                          boundary = list(bd1), 
                          refine = list(max.edge = 5))


mesh_c <-inla.mesh.create(loc = loc, 
                          boundary = list(bd1), 
                          refine = list(max.edge = 50))

# b,c similar



mesh_d <-inla.mesh.create(loc = loc, 
                          boundary = list(bd1), 
                          refine = list(max.edge = 0.05))



mesh_e <-inla.mesh.create(loc = loc, 
                          boundary = list(bd1), 
                          refine = list(max.edge = Inf))



ggplot() +
  gg(mesh_b) +
  coord_equal()



#==============================#
# Spatially Varying edge length
#==============================#

# first define a function that computes the 
# desired max.edge length as a function of loc
# and feed the output to that inla.mesh.create()
# using quality.spec parameter instead of edge

qual_loc <- function(loc) {
  pmax(0.05, (loc[, 1] * 2 + loc[, 2]) / 16)
}



qual_loc <- function(loc) {
  pmax(0.05, (loc[, 1] * 2 + loc[, 2]) / 16)
}


mesh_spv <- inla.mesh.create(loc = loc, boundary = list(bd1),
                 refine = list(max.edge = Inf),
                 quality.spec = list(loc = qual_loc(loc),
                                     segm = qual_loc(bd1$loc)))


ggplot() + gg(mesh_spv) + coord_equal()

## This gave a smooth transition between large and small triangles!


#=========================#
# Change boundary settings
#=========================#

# NA_real_ or Inf to make it not care about edge lengths near the boundary:

qual_bnd <- function(loc) {
  rep(Inf, nrow(loc))
}

mesh_infbd <- inla.mesh.create(loc = loc, 
                 boundary = list(bd1),
                 refine = list(max.edge = Inf),
                 quality.spec = list(loc = qual_loc(loc),
                                     segm = qual_bnd(bd1$loc)))
ggplot() + gg(mesh_infbd) + coord_equal()



#==========================#
# More complicated settings
#==========================#

# recommend setting the max.n.strict and max.n values in the refine parameter list, 
# that prohibits adding infinitely many triangles

qual_bnd <- function(loc) {
  pmax(0.1, 1 - abs(loc[, 2] / 10)^2)
}


mesh_complicate <- inla.mesh.create(loc = loc,
                 boundary = list(bd1),
                 refine = list(max.edge = Inf, 
                               max.n.strict = 5000),
                 quality.spec = list(loc = qual_bnd(loc),
                                     segm = qual_bnd(bd1$loc)))


ggplot() + gg(mesh_complicate) + coord_equal()



#=====================#
# inla.mesh.assessment
#=====================#

# assess the finite element approximation errors
out <- inla.mesh.assessment(mesh_complicate, 
                     spatial.range = 5,
                     alpha = 2, 
                     dims = c(200, 200))

print(names(out))

# [1] "sd"       "sd.dev"   "edge.len"

ggplot() + gg(out, aes(color = edge.len)) + coord_equal()



#============#
# Good mesh
#============#

# good mesh has std.dev close to 1
# this means nominal variance speicfied by
# continuous domain model and the variance 
# of discrete model are similar

ggplot() +
  gg(out, aes(color = sd.dev)) +
  scale_color_gradient(range(out$sd.dev, na.rm = T))


# standard deviation of the discretised model is smaller in the triangle interiors
# than at the vertices, 
# but the ratio is close to 1,
# more uniform where the triangles are small compared 
# with the spatial correlation length that we set to 5.


#The most problematic points seem to be in the 
# rapid transition between large and small triangles.
  # can improve things by increasing the minimum angle criterion:

mesh_minagl <- inla.mesh.create(loc = loc, 
                 boundary = list(bd1),
                 refine = list(min.angle = 30,
                               max.edge = Inf,
                               max.n.strict = 5000),
                 quality.spec = list(loc = qual_bnd(loc),
                                     segm = qual_bnd(bd1$loc)))

ggplot() + gg(mesh_minagl) + coord_equal()

out_minagl <- inla.mesh.assessment(mesh = mesh_minagl, 
                     spatial.range = 5,
                     alpha = 2, 
                     dims = c(200, 200))


ggplot() + gg(out_minagl, aes(color = out_minagl$sd.dev)) +
  scale_color_gradient(range(out_minagl$sd.dev, na.rm = T))

















