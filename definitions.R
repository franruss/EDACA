#function definitions
distToSim  <- function(x) {1-2*x^2}
simToDist <- function(x) {sqrt(1/2*(1-(x)))}
simFun <- function(x, y) {corCosine(x, y)} #Calculates the cosine correlation(s) between two vectors or conformable matrices.
#For centered data the cosine and Pearson's correlation coefficients are equivalent.  
#cosine similarity is defined as similarity: cos(??)=(??? XY) \div (???(??? X^2)\times ???(??? Y^2))

ntp_mine <- function (emat, templates, nPerm = 1000, distance = "cosine", 
                      nCores = 1, seed = NULL, verbose = getOption("verbose"), 
                      doPlot = FALSE) 
{
  if (class(emat)[1] == "ExpressionSet") 
    emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
  if (is.data.frame(emat) | is.vector(emat)) 
    emat <- as.matrix(emat)
  if (is.null(rownames(emat))) 
    stop("missing emat rownames, check input")
  if (is.null(templates$class) | is.null(templates$probe)) {
    stop("missing columns in templates, check input")
  }
  if (is.character(templates$class)) 
    templates$class <- factor(templates$class)
  if (is.factor(templates$probe) | is.integer(templates$probe)) {
    warning("templates$probe coerced to character", call. = FALSE)
    templates$probe <- as.character(templates$probe)
  }
  keepP <- stats::complete.cases(emat)
  if (sum(!keepP) > 0) {emat <- emat[keepP, , drop = FALSE]}
  keepT <- templates$probe %in% rownames(emat)
  if (sum(!keepT) > 0) {
    if (isTRUE(verbose)){
      message(paste0(sum(!keepT), "/", length(keepT), 
                     " templates features not in emat, discarded"))
    }
    templates <- templates[keepT, ]
  }
  if (min(table(templates$class)) < 2) {
    message("<2 matched features/class")
    stop("check templates$probe is matchable against rownames(emat)", 
         call. = FALSE)
  }
  if (min(table(templates$class)) < 5) {
    warning("<5 matched features/class - unstable predictions", 
            call. = FALSE)
  }
  N <- ncol(emat)
  K <- nlevels(templates$class)
  S <- nrow(templates)
  P <- nrow(emat)
  class.names <- levels(templates$class)
  templates$class <- as.numeric(templates$class) #qui trasformo tamplate_class in numeri per costruire dopo twat
  emat.mean <- round(mean(emat), 2)
  if (abs(emat.mean) > 1) {
    isnorm <- " <- check feature centering!"
    emat.sd <- round(stats::sd(emat), 2)
    warning(paste0("emat mean=", emat.mean, "; sd=", 
                   emat.sd, isnorm), call. = FALSE)
  }
  feat.class <- paste(range(table(templates$class)), collapse = "-")
  mm <- match(templates$probe, rownames(emat), nomatch = 0) # prendi gli indici dei predittori da qui
  if (!all(rownames(emat)[mm] == templates$probe)) {stop("error matching probes, check rownames(emat) and templates$probe")}
  pReplace <- length(templates$probe) > length(unique(templates$probe))
  tmat <- matrix(rep(templates$class, K), ncol = K)
  for (k in (1:K)){tmat[, k] <- as.numeric(tmat[, k] == k)}
  tmat # design matrix
  if (K == 2){tmat[tmat == 0] <- -1}
  res=NULL
  #ntpFUN <- function(n) { # n corre su tutti i campioni
  for (n in (1:N)){
    n.sim <- as.vector(simFun( emat[mm, n, drop = FALSE] ,tmat)) #Perch? calcolare la correlazione contro una design matrix? 
    #Per sapere in quali campioni i geni upregolati sono pi? vicini a 1 e quelli down regolati pi? vicini a zero.
    #plot(emat[mm, n, drop = FALSE] ,tmat[,1])
    n.sim.perm.max <- apply(simFun(matrix(emat[, n][sample.int(P, S * nPerm, replace = TRUE)], ncol = nPerm), tmat), 1, max) 
    # adesso fai delle permutazioni e crea una matrice di valori random partendo da quelli osservati
    #e prendi il massimo per ogni riga
    n.ntp <- which.max(n.sim) #seleziona la classe che ha ottenuto il massimo score come similitudine nei valori reali
    n.sim.ranks <- rank(-c(n.sim[n.ntp], (n.sim.perm.max)))
    n.pval <- n.sim.ranks[1]/length(n.sim.ranks)
    res[[n]] <- (c(n.ntp, simToDist(n.sim), n.pval))
  }
  #res <- lapply(seq_len(N), ntpFUN)
  res <- data.frame(do.call(rbind, res))
  colnames(res) <- c("prediction", paste0("d.", class.names), "p.value")
  res$prediction <- factor(class.names[res$prediction], levels = class.names)
  rownames(res) <- colnames(emat)
  res$p.value[res$p.value < 1/nPerm] <- 1/nPerm
  res$FDR <- stats::p.adjust(res$p.value, "fdr")
  return(res)
}
##### END definitions