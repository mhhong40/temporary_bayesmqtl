#include "utils.h"
// Not used in practice as there's a bug here that causes RStudio to crash... fix/find out ASAP

// Parameter theta_st is the most computationally expensive to update, so its update is coded here.
// Hélène -- this is the code you wrote for locus's coreLoop(). I directly adapted it here to save time, sorry.
// [[Rcpp::export]]
void coreLoop(const MapMat z,
              const MapMat X,
              MapMat mat_x_m1,
              MapMat theta_vb,
              MapArr2D mu_theta_vb,
              const MapArr1D tau_vb,
              const MapArr2D sig2_theta_vb,
              const double c = 1) {

  for (int s = 0; s < X.cols(); ++s) {

   // For some reason, this code runs very slowly -- what's the performance bottleneck?
   mat_x_m1.noalias() -= X.col(s) * theta_vb.row(s); // take away the contribution of the sth predictor

   mu_theta_vb.row(s) = c * sig2_theta_vb.row(s) * tau_vb *
     ((z - mat_x_m1).transpose() * X.col(s)).array();

   theta_vb.row(s) = mu_theta_vb.row(s);

   mat_x_m1.noalias() += X.col(s) * theta_vb.row(s); // add back s's contribution for the next predictor.
  }

}
