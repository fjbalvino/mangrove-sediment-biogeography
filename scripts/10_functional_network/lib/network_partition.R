# Deterministic multistart partitioning; the caller's sampling RNG is preserved.
fit_network_partition <- function(g, method = "leiden", resolution = 1,
                                  seeds = 2001:2050, weights = igraph::E(g)$weight) {
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv)
           else if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))
             rm(".Random.seed", envir=.GlobalEnv)}, add=TRUE)
  v <- igraph::V(g)$name
  if (igraph::ecount(g) == 0L) return(list(membership=setNames(seq_along(v),v),
    best_seed=NA_integer_, Q=NA_real_, fits=data.frame(), community=NULL))
  if (!method %in% c("leiden","louvain")) seeds <- seeds[1L]
  Q <- numeric(length(seeds)); sizes <- integer(length(seeds))
  best <- NULL; bestQ <- -Inf; best_seed <- NA_integer_
  memberships <- matrix(NA_integer_,length(v),length(seeds),dimnames=list(v,NULL))
  for (j in seq_along(seeds)) {
    set.seed(seeds[j])
    cl <- switch(method,
      leiden=igraph::cluster_leiden(g,objective_function="modularity",weights=weights,
                                   resolution=resolution,n_iterations=-1),
      louvain=igraph::cluster_louvain(g,weights=weights,resolution=resolution),
      infomap=igraph::cluster_infomap(g,e.weights=weights),
      walktrap=igraph::cluster_walktrap(g,weights=weights),
      stop("Unknown community method: ",method))
    m <- as.integer(igraph::membership(cl)); memberships[,j] <- m
    Q[j] <- igraph::modularity(g,m,weights=weights,resolution=resolution)
    sizes[j] <- sum(table(m[igraph::degree(g)>0])>=10L)
    if(Q[j]>bestQ) { best <- cl; bestQ <- Q[j]; best_seed <- seeds[j] }
  }
  list(membership=setNames(as.integer(igraph::membership(best)),v),
       best_seed=best_seed,Q=bestQ,community=best,memberships=memberships,
       fits=data.frame(method=method,resolution=resolution,seed=seeds,
                       modularity=Q,n_major_modules=sizes,selected=seeds==best_seed))
}
