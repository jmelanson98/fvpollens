library(glmmTMB)
library(lme4)
library(performance)

run_plot_freq_model_diagnostics <- function(this_formula, #brms model formula
                                            this_data, #data frame, subsetted to correct weights!
                                            num_chains=1, 
                                            num_iter=10^4, 
                                            launch.shiny=FALSE,
                                            examine.pairs=FALSE,
                                            this_family, #model family
                                            fig.path =
                                            "figures/diagnostics",
                                            site.lat="lat", ...
                                            ){

    
    ## function to run frequentist version of brms models and plot diagnostics

    ## currently works with all families included in ?family, as well as
    ## zero inflated neg binomial, negbinomial, Gamma, inverse gaussian,
    ## poisson, quasibinomial will need to add other families as needed
    ## uses other glm packages to model frequentist versions of the
    ## bayesian models

    if (this_family == 'gaussian'){
                                        #run model
        this_model_output <- brms::brm(this_formula,
                                       data = this_data, 
                                       chains = num_chains, 
                                       iter = num_iter, family=this_family,
                                       thin=1,
                                       init=0,
                                       open_progress = FALSE,
                                       control = list(adapt_delta = 0.99),
                                       save_pars = save_pars(all = TRUE))

        
        print(summary(this_model_output))
                                        # return a list of single plots
        diagnostic.plots <- plot(check_model(this_model_output, panel
                                             = TRUE))
    } else if (this_family=='negbinomial') {
        
        this_model_output <- glmer.nb(this_formula, data=this_data)
        diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
        
        
        

    } else if (this_family=='zero_inflated_binomial') {
        
        this_model_output <- glmmTMB(this_formula, data=this_data,
                                     ziformula = ~1,
                                     family = binomial)
        diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
        
        
        

    } else if (this_family=='hurdle_gamma') {
        
        this_model_output <- glmer(this_formula, data=this_data, family=Gamma)
        diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
        
        

    } else if (this_family=='hurdle_negbinomial') {
        
        this_model_output <- glmmTMB(this_formula, data=this_data,
                                     ziformula = ~1,
                                     family = truncated_nbinom2())
        diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
        

    } else if (this_family=='students') {
        
        this_model_output <- glmmTMB(this_formula, data=this_data, family = t_family())
        diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
        
        

    } else if (this_family=='negbinomial') {
        this_model_output <- glmer.nb(this_formula, data=this_data, ...)
        diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
        

    }else if (this_family=='bernoulli') {
        this_model_output <- glmmTMB(this_formula, data=this_data, family='binomial')
        diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
        
        

    } else if (this_family == 'Gamma'){
                                        #run model
        this_model_output <- glmer(this_formula, data=this_data, family=Gamma(link=log))
        
                                        # return a list of single plots
        diagnostic.plots <- plot(check_model(this_model_output, 
                                             panel = TRUE))


    } else if (this_family == 'inverse.gaussian'){
                                        #run model
        this_model_output <- brms::brm(this_formula,
                                       data = this_data, 
                                       chains = num_chains, 
                                       iter = num_iter, family=this_family,
                                       thin=1,
                                       init=0,
                                       open_progress = FALSE,
                                       control = list(adapt_delta = 0.99),
                                       save_pars = save_pars(all = TRUE))
                                        # return a list of single plots
        diagnostic.plots <- plot(check_model(this_model_output, 
                                             panel = TRUE))


    } else if (this_family == 'poisson'){
                                        #run model
        this_model_output <- brms::brm(this_formula,
                                       data = this_data, 
                                       chains = num_chains, 
                                       iter = num_iter, family=this_family,
                                       thin=1,
                                       init=0,
                                       open_progress = FALSE,
                                       control = list(adapt_delta = 0.99),
                                       save_pars = save_pars(all = TRUE))
                                        # return a list of single plots
        diagnostic.plots <- plot(check_model(this_model_output, 
                                             panel = TRUE), type = "discrete_dots")

    } else if (this_family == 'quasibinomial'){
                                        #run model
        this_model_output <- brms::brm(this_formula,
                                       data = this_data, 
                                       chains = num_chains, 
                                       iter = num_iter, family=this_family,
                                       thin=1,
                                       init=0,
                                       open_progress = FALSE,
                                       control = list(adapt_delta = 0.99),
                                       save_pars = save_pars(all = TRUE))
                                        # return a list of single plots
        diagnostic.plots <- plot(check_model(this_model_output, 
                                             panel = TRUE))
    } 
    if (launch.shiny == TRUE){shinystan::launch_shinystan(this_model_output)}
    if (examine.pairs == TRUE){pairs(this_model_output)}

    if(class(this_formula)[1] == "brmsformula"){
        file.name <- paste0(as.character(this_formula)[4], site.lat, ".pdf")
    }else{
        file.name <- paste0(as.character(this_formula)[2], site.lat, ".pdf")
    }
    ggsave(diagnostic.plots, file= file.path(fig.path, file.name))
    return(diagnostic.plots)
}
