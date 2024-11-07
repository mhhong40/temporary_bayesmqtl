## Running a basic check on bayesmqtl using a subset of the first replicate of the mixture model-based dataset.
devtools::load_all()
X_full <- readRDS("main_ls copy.rds")
X_full <- X_full[[1]]$snps

Y_mix_reps <- readRDS("Y_mix_paper copy.rds")
Y_mix <- Y_mix_reps[[1]]

indices <- sample(1:20000, 600)
Y_mix <- Y_mix[, indices]
indices <- sample(1:1000, 400)
X <- X_full[, indices]
rm(indices)
# Test bayesmqtl() performance and runtime, run serially because of the small problem size.

# User must set list_hyper, list_init... no default settings yet, sorry.
n <- nrow(X)
p <- ncol(X)
d <- ncol(Y_mix)
list_hyper <- set_hyper(d, p, X, kappa = 1, eta = 1, lambda = 1)

theta <- matrix(rnorm(p*d), nrow = p, ncol = d)
list_init <- set_init(d, p, theta, tau = 1)

bounds <- quantile(colMeans(Y_mix), probs = seq(0, 1, 0.01)) # a very simple partition
bounds <- bounds[-1] # remove the 0th percentile since that is a lower, not upper bound
bin_ind <- generate_bin_ind(Y_mix, bounds)

bin_mix_prop <- rep(0.5, length(bin_ind)) # ** mix_prop are VERY important!
list_shape <- estimate_shape_params(Y_mix, bin_ind, bin_mix_prop)

# Do externally because for some reason, these vectors don't update and
# get passed into the core function when within bayesmqtl()?

# TO DO: Turn this into its own function and pass into bayesmqtl as list_shape
alpha_0 <- beta_0 <- alpha_1 <- beta_1 <- mix_prop <- rep(0, d)
for (j in 1:length(bin_ind)) {

  alpha_0 <- replace(alpha_0, bin_ind[[j]], list_shape$alpha_0[j])
  beta_0 <- replace(beta_0, bin_ind[[j]], list_shape$beta_0[j])
  alpha_1 <- replace(alpha_1, bin_ind[[j]], list_shape$alpha_1[j])
  beta_1 <- replace(beta_1, bin_ind[[j]], list_shape$beta_1[j])
  mix_prop <- replace(mix_prop, bin_ind[[j]], bin_mix_prop[j])
}

# VI on the membership probability linear regression coefficients
z <- matrix(rnorm(n*d), nrow = n, ncol = d) # test to see
z <- estimate_mem_probs_(Y_mix, z, mix_prop, alpha_0, beta_0, alpha_1, beta_1)
z <- qnorm(z)
z <- scale(z, center = TRUE, scale = FALSE)

# Need to scale X in order to use n - 1 as ||X||^2
X <- scale(X)

anneal <- c(1, 5, 30)
vb <- bayesmqtl(Y_mix, X, z, bin_ind, mix_prop, alpha_0, beta_0, alpha_1, beta_1,
                list_hyper, list_init, tol = 0.1, anneal = anneal)
