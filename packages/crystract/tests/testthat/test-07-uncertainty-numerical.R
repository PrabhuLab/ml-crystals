test_that("distance and angle gradients agree with independent finite differences", {
  for (angles in list(c(90, 90, 90), c(79, 83, 74))) {
    for (shared_parent in c(TRUE, FALSE)) {
    cell <- data.table::data.table(
      `_cell_length_a`=7, `_cell_length_b`=8, `_cell_length_c`=9,
      `_cell_angle_alpha`=angles[1], `_cell_angle_beta`=angles[2], `_cell_angle_gamma`=angles[3])
    atoms <- data.table::data.table(Label=c("A", "B"), x_a=c(.1,.2), y_b=c(.05,.12), z_c=c(.07,.3))
    # Both neighbours derive from B; each sample regenerates that dependency.
    images <- data.table::data.table(Label=c("B1","B2"), x_a=c(.2,.8), y_b=c(.12,.88), z_c=c(.3,.7),
       ParentJacobian=list(list(B=diag(3)),list(B=-diag(3))))
    bond <- data.table::data.table(Atom1="A", Atom2="B2", Distance=1)
    angle <- data.table::data.table(CentralAtom="A", Neighbor1="B1", Neighbor2="B2", Angle=1)
    keys <- c(paste0("cell_",c("a","b","c","alpha","beta","gamma")),"A:x","A:y","A:z","B:x","B:y","B:z")
    p <- c(7,8,9,angles*pi/180,.1,.05,.07,.2,.12,.3)
    if (!shared_parent) {
      atoms <- rbind(atoms, data.table::data.table(Label="C", x_a=.8, y_b=.88, z_c=.7))
      images$ParentJacobian[[2]] <- list(C=diag(3))
      keys <- c(keys, "C:x", "C:y", "C:z")
      p <- c(p, .8, .88, .7)
    }
    n_input <- length(p)
    evaluate <- function(p) {
      ca<-cos(p[4]); cb<-cos(p[5]); cg<-cos(p[6]); sg<-sin(p[6])
      D<-1-ca^2-cb^2-cg^2+2*ca*cb*cg
      M<-rbind(c(p[1],p[2]*cg,p[3]*cb), c(0,p[2]*sg,p[3]*(ca-cb*cg)/sg), c(0,0,p[3]*sqrt(D)/sg))
      u<-as.numeric(M%*%(p[10:12]-p[7:9]));v<-as.numeric(M%*%((if (length(p)==15) p[13:15] else 1-p[10:12])-p[7:9]))
      c(distance=sqrt(sum(v^2)),angle=acos(sum(u*v)/sqrt(sum(u^2)*sum(v^2)))*180/pi)
    }
    reference <- NULL
    for (h in c(1e-5, 3e-6)) {
      J <- sapply(seq_along(p),function(i) {
        plus<-minus<-p;plus[i]<-plus[i]+h;minus[i]<-minus[i]-h
        (evaluate(plus)-evaluate(minus))/(2*h)
      })
      if (!is.null(reference)) expect_equal(J,reference,tolerance=1e-7)
      reference<-J
      for (i in seq_along(p)) {
        Sigma<-matrix(0,n_input,n_input,dimnames=list(keys,keys));Sigma[i,i]<-1
        expect_equal(propagate_distance_error(bond,atoms,cell,images,covariance=Sigma)$DistanceError,unname(abs(J[1,i])),tolerance=1e-7)
        expect_equal(propagate_angle_error(angle,atoms,images,cell,covariance=Sigma)$AngleError,unname(abs(J[2,i])),tolerance=1e-7)
      }
    }
    # A correlated input matrix also checks cross-terms and relative signs.
    Sigma <- diag(n_input)*1e-6 + outer(seq_len(n_input),seq_len(n_input))*1e-8
    dimnames(Sigma)<-list(keys,keys)
    expected<-sqrt(diag(J%*%Sigma%*%t(J)))
    expect_equal(propagate_distance_error(bond,atoms,cell,images,covariance=Sigma)$DistanceError,unname(expected[1]),tolerance=1e-7)
    expect_equal(propagate_angle_error(angle,atoms,images,cell,covariance=Sigma)$AngleError,unname(expected[2]),tolerance=1e-7)
    }
  }
})

test_that("merged sites retain both parent derivatives through cell expansion", {
  cell <- data.table::data.table(`_cell_length_a`=10, `_cell_length_b`=10, `_cell_length_c`=10,
    `_cell_angle_alpha`=90, `_cell_angle_beta`=90, `_cell_angle_gamma`=90)
  sites <- data.table::data.table(Label=c("A1","B1"),x_a=c(.1,.100001),y_b=.2,z_c=.3,
    ParentJacobian=list(list(A=diag(3)),list(B=diag(3))))
  merged <- merge_sites_pbc(sites,cell,tol=1e-4)
  expect_equal(nrow(merged),1)
  expect_equal(merged$ParentJacobian[[1]]$A,diag(3)/2)
  expect_equal(merged$ParentJacobian[[1]]$B,diag(3)/2)
  images <- expand_transformed_coords(merged)
  expect_true(all(vapply(images$ParentJacobian,identical,logical(1),merged$ParentJacobian[[1]])))
})
