
standardize.axis <- function(x, orig)
(x-mean(orig, na.rm=TRUE))/sd(orig, na.rm=TRUE)


plot.res <- function(mod, mod.name){
    ## function tp plot diagnostic figures for mcmc
    pdf(sprintf("figures/diagnostics/%s_Diag.pdf", mod.name),
        height=11, width=8.5)
    plot(mod,  N = 4, ask = FALSE)
    dev.off()
}

plot.res <- function(mod, mod.name){
    ## function tp plot diagnostic figures for mcmc
    pdf(sprintf("figures/diagnostics/%s_Diag.pdf", mod.name),
        height=11, width=8.5)
    plot(mod,  nvariables = 4, ask = FALSE)
    dev.off()
}
