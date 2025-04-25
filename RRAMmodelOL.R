
library(RTMB)

##Censored bounds
##We use the same code as TMB as a workaround and make it into an RTMB tape
TMB::compile("pnormtape.cpp")
dyn.load(TMB::dynlib("pnormtape"))
##make some fake data to make it 
parm = list(x=1,mu=0,sd=1,lower=0.5,upper=2)
data = list()
pnormobj = TMB::MakeADFun(data,parameters = parm,DLL="pnormtape")

pnormtape = RTMB::GetTape(pnormobj)

##nice function to use
censored_bounds <- function(x, mu, sd, lower, upper){
    pnormtape(c(x,mu,sd,lower,upper))
}


##Convert from optimizer space to -1 to 1 for correlation parameters
rho_trans <- function(x){
    2/(1+exp(-2*x))-1
}

## ##Projection F solver
solveF <- function(catch_option,Finit,S,M,N){
    niter = 7
    F = Finit
    for(i in 1:niter){
        Z = F*S+M
        eZ = exp(-1*Z)
        Mort = 1 - eZ
        CF = F * sum((N* Mort* S/Z))
        derCF = sum(N*S*(Mort*M/Z+eZ*F*S)/Z)
        F = F + ((catch_option-CF)/derCF)
        F
    }
    F
}
    
##aggregate catch function
aggregate_catch <- function(catch_vec,catch_key,n_size){
    aggr_catch = numeric(n_size)

    for(i in 1:n_size){
        aggr_catch[i] = 0
    }

    for(i in 1:length(catch_vec)){
        aggr_catch[catch_key[i]] = aggr_catch[catch_key[i]]+catch_vec[i]
    }

    aggr_catch
    
}

##Convert proportions to CRLs
make_CRL_vec <- function(proportions){
    S = length(proportions)
    crl = numeric(S-1)

    for(ss in 1:(S-1)){
        pi_b = proportions[ss]/sum(proportions[ss:S])
        crl[ss] = log(pi_b/(1-pi_b))
    }

    crl
    
}

##Growth model functions


mu_len <- function(gf,t){
    mu_t = with(gf,ell+(L-ell)*((1-k^(t-1)))/(1-k^(M-1)))
    mu_t
}

sd_len <- function(gf,t){
    sd_t = with(gf,s+(S-s)*((t-ell)/(L-ell)));
    sd_t
}

mean_len_at_age <- function(gf,t,y,gf_ext){
    if(t >= gf$M & gf_ext){
        weights = numeric(gf$pg_ext)
        mus = numeric(gf$pg_ext)
        for(i in 1:length(mus)){
            weights[i] = dexp(i+(t-gf$M)-1,1/(gf$as_pg[y]))
            mus[i] = mu_len(gf,t+i-1)
        }
        REPORT(weights)
        musclus = mus
        REPORT(musclus)
        mu_t = sum(mus*weights)/sum(weights)
        ##mu_t = mu_len(gf,t)
    }else{
        mu_t = mu_len(gf,t)
    }
    mu_t
}

sd_len_at_age <- function(gf,t,y,gf_ext){
    if(t >= gf$M & gf_ext){
        weights = numeric(gf$pg_ext)
        sds = numeric(gf$pg_ext)
        for(i in 1:length(sds)){
            weights[i] = dexp(i+(t-gf$M)-1,1/(gf$as_pg[y]))
            sds[i] = sd_len(gf,t+i-1)
        }
        sd_t = sum(sds*weights)/sum(weights)
        ##sd_t = sd_len(gf,t)
    }else{
        sd_t = sd_len(gf,t)
    }
    sd_t
}

prob_len_at_age <- function(gf_mus,gf_sds,length_bins,ages){
    A = length(ages)
    L = length(length_bins)

    pla = matrix(NA,L,A)

    for(a in 1:A){
        for(l in 1:L){
            if(l == 1){
                pla[l,a] = pnorm(length_bins[l],gf_mus[a],gf_sds[a])
            }else if(l == L){
                pla[l,a] = 1-pnorm(length_bins[l-1],gf_mus[a],gf_sds[a])
            }else{
                pla[l,a] = pnorm(length_bins[l],gf_mus[a],gf_sds[a])-pnorm(length_bins[l-1],gf_mus[a],gf_sds[a])
            }
            
        }
    }

    pla

}


prob_len_at_age_gamma <- function(gf_mus,gf_sds,length_bins,ages){
    A = length(ages)
    L = length(length_bins)

    pla = matrix(NA,L,A)

    for(a in 1:A){
        for(l in 1:L){
            scaly = gf_sds[a]^2/gf_mus[a]
            shapy = gf_mus[a]^2/gf_sds[a]^2
            
            if(l == 1){
                pla[l,a] = pgamma(length_bins[l],shape=shapy,scale=scaly)
            }else if(l == L){
                pla[l,a] = 1-pgamma(length_bins[l-1],shape=shapy,scale=scaly)
            }else{
                pla[l,a] = pgamma(length_bins[l],shape=shapy,scale=scaly)-pgamma(length_bins[l-1],shape=shapy,scale=scaly)
            }
            
        }
    }

    pla

}


##Upper bounded inverse transform
ub_inv <- function(y, b){
    x = b - exp(y)
    x
}

##Assume we have a vector where we want x_k < x_{k+1}
ordered_inv_transform <- function(y){
    x = y
    for(i in 2:length(x)){
        x[i] = x[i-1]+exp(y[i])
    }
    x
}


##Multivariate t distribution
dmvt <- function(x,mu,sigma,df,give_log=FALSE){
    p = length(x)
    top_d = gamma((df+p)/2)
    bottom_d = gamma(df/2)*df^(p/2)*pi^(p/2)*exp(determinant(sigma,logarithm = TRUE)$modulus)^(1/2)
    sig_inv = solve(sigma)
    big_b = 1+1/df*t(x-mu)%*%sig_inv%*%(x-mu)
    res = (top_d/bottom_d)*(big_b)^(-(df+p)/2)
    if(give_log){
        res = log(res)
    }
    res
}

##dat = exampleS$data

##Fix all survey sds on one common thing
## for(x in 1:length(dat$survey_list)){
##     dat$survey_list[[x]]$mmap = rep(0,length(dat$survey_list[[x]]$mmap))
## }
## exampleS$parameters$log_sd_survey = exampleS$parameters$log_sd_survey[1] 

rram_model <- function(parameters,dat){

    ##I don't like the way getAll works...
    ##parameters
    log_Q_max_sd = parameters$log_Q_max_sd
    log_ell = parameters$log_ell
    log_L = parameters$log_L
    log_k = parameters$log_k
    log_S = parameters$log_S
    log_s = parameters$log_s
    log_recruit_sd = parameters$log_recruit_sd
    log_N0_sd = parameters$log_N0_sd
    log_Fy = parameters$log_Fy
    log_Fy_sd = parameters$log_Fy_sd
    log_Qmax = parameters$log_Qmax
    log_QL50 = parameters$log_QL50
    log_QL95 = parameters$log_QL95
    log_landings_sd = parameters$log_landings_sd
    engel.log_Qmax = parameters$engel.log_Qmax
    engel.log_QL50 = parameters$engel.log_QL50
    engel.log_QL95 = parameters$engel.log_QL95
    log_N_a = parameters$log_N_a
    log_surv_sd = parameters$log_surv_sd
    log_b_beta1 = parameters$log_b_beta1
    log_sel_scale = parameters$log_sel_scale
    log_sel_shape = parameters$log_sel_shape
    log_sd_catch_prop_m = parameters$log_sd_catch_prop_m
    logit_rhoC = parameters$logit_rhoC
    log_init_a_pg = parameters$log_init_a_pg
    logit_rhoS = parameters$logit_rhoS
    log_sd_survey = parameters$log_sd_survey

    ##asymmetric logistic parameters
    log_delta_survey = parameters$log_delta_survey
    delta_survey = exp(log_delta_survey)
    log_b_beta2 = parameters$log_b_beta2

    

    ##Data
    ##Not all necessarily used
    start_year = dat$start_year
    end_year = dat$end_year
    Y = dat$Y
    A = dat$A
    L = dat$L
    L3 = dat$L3
    inf_length = dat$inf_length
    weightsF = dat$weightsF
    weightsM = dat$weightsM
    maturityF = dat$maturityF
    maturityM = dat$maturityM
    survey_index = dat$survey_index
    survey_year = dat$survey_year
    survey_length = dat$survey_length
    survey_type = dat$survey_type
    landing_year = dat$landing_year
    landing_nums = dat$landing_nums
    bin_adjust = dat$bin_adjust
    adj_landings = dat$adj_landings
    discards = dat$discards
    log_landingsL = dat$log_landingsL
    log_landingsU = dat$log_landingsU
    dum_indices = dat$dum_indices
    Qpriors = dat$Qpriors
    vagesindices = dat$vagesindices
    vyearsindices = dat$vyearsindices
    kages = dat$kages
    kyears = dat$kyears
    catch_list = dat$catch_list
    base_M = dat$base_M
    Q_prior = dat$Q_prior
    mora_year = dat$mora_year
    rounding_bit = dat$rounding_bit
    catch_type = dat$catch_type
    catch_sig_sw = dat$catch_sig_sw
    pg_ext = dat$pg_ext
    survey_list = dat$survey_list
    survey_size = dat$survey_size
    r_proj = dat$r_proj
    landing_proj_y = dat$landing_proj_y
    cfg = dat$cfg
    proj_type = dat$proj_type
    proj_years = dat$proj_years
    og_Y = dat$og_Y
    supplied_F = dat$supplied_F
    given_catch = dat$given_catch
    gf_ext = dat$gf_ext
    sel_type = dat$sel_type
    l_dist_type = dat$l_dist_type
    s_dist_type = dat$s_dist_type
    c_dist_type = dat$c_dist_type
    plus_s = dat$plus_s

    if(l_dist_type == "normal"){
        plaa = prob_len_at_age
    }else{
        plaa = prob_len_at_age_gamma
    }
    
    
    start_age = 1
    fall_adjust = 10.5/12
    mid_adjust = 0.5

    nll = 0.0

    N_a = exp(log_N_a)

    REPORT(log_N_a)

    ##Calculate the value of the mean age in the plus group
    if(gf_ext){
        as_pg = numeric(Y)
        as_pg[1] = ub_inv(log_init_a_pg,pg_ext)
        for(y in 2:Y){
            ##We can do this because log_N_a is a parameter (random effect)
            as_pg[y] = (exp(log_N_a[A-1,y-1])*1+exp(log_N_a[A,y-1])*(as_pg[y-1]+1))/(exp(log_N_a[A-1,y-1])+exp(log_N_a[A,y-1]));

        }
        REPORT(as_pg)
    }else{
        as_pg = 0
    }
    ##Create the vector to hold the GF data and parameters and do transformations

    ## We want ell ALWAYS less than L
    ell_L = exp(ordered_inv_transform(c(log_ell,log_L)))
    gf_par = list(ell = ell_L[1],
               L = ell_L[2],
               s = exp(log_s),
               S = exp(log_S),
               k = exp(log_k),
               M = A,
               pg_ext = pg_ext,
               as_pg = as_pg)

    
    len_bins = 1:L3+bin_adjust
    ages_jan = start_age:(A)
    ages_mid = ages_jan + mid_adjust
    ages_fall = ages_jan + fall_adjust

    len_sds = matrix(NA,A,Y)
    len_mus = matrix(NA,A,Y)

    len_sds_mid = matrix(NA,A,Y)
    len_mus_mid = matrix(NA,A,Y)

    len_sds_fall = matrix(NA,A,Y)
    len_mus_fall = matrix(NA,A,Y)
    

    for(a in 1:A){
        for(y in 1:Y){
            len_sds[a,y] = sd_len_at_age(gf_par,ages_jan[a],y,gf_ext)
            len_mus[a,y] = mean_len_at_age(gf_par,ages_jan[a],y,gf_ext)
            len_sds_mid[a,y] = sd_len_at_age(gf_par,ages_mid[a],y,gf_ext)
            len_mus_mid[a,y] = mean_len_at_age(gf_par,ages_mid[a],y,gf_ext)
            len_sds_fall[a,y] = sd_len_at_age(gf_par,ages_fall[a],y,gf_ext)
            len_mus_fall[a,y] = mean_len_at_age(gf_par,ages_fall[a],y,gf_ext)
        }
    }


    REPORT(len_mus)
    REPORT(len_sds)
    
    ##Lorenzen M
    M = matrix(NA,A,Y)
    
    mus = matrix(NA,A,Y)
    for(y in 1:Y){
        mu_f = mean_len_at_age(gf_par,A+0.5,y,gf_ext)
        for(a in 1:A){
            mu = mean_len_at_age(gf_par,a+0.5,y,gf_ext)
            mus[a,y] = mu
            M[a,y] = base_M*((mu^(-1))/(mu_f^(-1)))
        }
    }

    REPORT(mus)
    REPORT(M)
    REPORT(mu_f)        

    ##Probability of length at age
    plas_jan = list(rep(NA,Y))
    plas_mid = list(rep(NA,Y))
    plas_fall = list(rep(NA,Y))

    for(y in 1:Y){
        plas_jan[[y]] = plaa(len_mus[,y],len_sds[,y],len_bins,ages_jan)
        plas_mid[[y]] = plaa(len_mus_mid[,y],len_sds_mid[,y],len_bins,ages_mid)
        plas_fall[[y]] = plaa(len_mus_fall[,y],len_sds_fall[,y],len_bins,ages_fall)
    }

    ##Selectivity at length
    S_ly = matrix(NA,L3,Y)

    logi_key = dat$logi_key

    for(y in 1:(mora_year)){

        ##Logistic Pre-moratorium time period
        if(rounding_bit != FALSE){
            b_beta1 = exp(log_b_beta1)
            b_beta2 = b_beta1+rounding_bit*b_beta1
        }else{
            u_b_beta = exp(ordered_inv_transform(c(log_b_beta1[logi_key[y]],log_b_beta2[logi_key[y]])))
            b_beta1 = u_b_beta[1]
            b_beta2 = u_b_beta[2]
        }
        r_b_beta2 = -1*(b_beta2-b_beta1)/(log(0.05/0.95))


        
        for(l in 1:L3){
            l_bin = len_bins[l]
            interior = -(1/r_b_beta2)*(l_bin-b_beta1)
            bot = (1+exp(interior))
            S_ly[l,y] = (1/bot)
        }
    }
    ##Gamma shaped post-moratorium
    sel_shape = exp(log_sel_shape)
    sel_scale = exp(log_sel_scale)
    gamma_key = dat$gamma_key
    for(y in (mora_year+1):Y){
        for(l in 1:L3){
            l_bin = len_bins[l]
            S_ly[l,y] = dgamma(l_bin,shape=sel_shape[gamma_key[y]],scale=sel_scale[gamma_key[y]])
        }
    }
        
        ##Adjust gamma shaped to be between 0 and 1
        for(y in (mora_year:Y)){
            S_ly[,y] = S_ly[,y]/max(S_ly[,y])
        }
        

    
    ##Convert S_ly to S_ay
    S_ay = matrix(NA,A,Y)
    for(y in 1:Y){
        pla_y = plas_mid[[y]]
        for(a in 1:A){
            S_ay[a,y] = 0
            for(l in 1:L3){
                S_ay[a,y] = S_ay[a,y] + pla_y[l,a]*S_ly[l,y]
            }
        }
    }

    ##The overall F value
    F = matrix(NA,A,Y)

    Fy_nll = numeric(length(log_Fy))
    Fy_nll[1] = 0
    for(y in 2:length(log_Fy)){
        nll = nll - dnorm(log_Fy[y],log_Fy[y-1],exp(log_Fy_sd),TRUE)
        Fy_nll[y] = dnorm(log_Fy[y],log_Fy[y-1],exp(log_Fy_sd),TRUE)
    }

    REPORT(Fy_nll)
    Fy = exp(log_Fy)

    for(a in 1:A){
        for(y in 1:og_Y){
            F[a,y] = Fy[y]*S_ay[a,y]
        }
    }

    ##Peice for projections
    if(proj_type == 0){
        proj_Fy = exp(supplied_F)
    }

    if(Y > og_Y){
        for(a in 1:A){
            for(y in (og_Y+1):Y){
                if(proj_type == 0){
                    yc = y-og_Y
                    F[a,y] = proj_Fy[yc]*S_ay[a,y]
                }else{
                    F[a,y] = 0
                }
            }
        }
    }
    ##End peice for projections

    Z = F+M

    ##Cohort Effect N
    ## In theory we can do correlated N stuff here
    ## We don't, but we COULD
    ##surv_sd = exp(log_surv_sd)
    N0_sd = exp(log_N0_sd)

    ##predicted mean
    mpredN = matrix(NA,A,Y)
    ##dynamics one
    mpredN2 = matrix(NA,A,Y)
    ##observed one
    mcomp = matrix(NA,A,Y)

    mpredN2[1,1] = log_N_a[1,1]
    mpredN[1,1] = log_N_a[1,1]
    mcomp[1,1] = log_N_a[1,1]

    ##Recruits
    for(y in 2:Y){
        mpredN[1,y] = log_N_a[1,y-1]
        mpredN2[1,y] = log_N_a[1,y-1]
        mcomp[1,y] = log_N_a[1,y]
    }

    ##initial numbers at age
    for(a in 2:A){
        mpredN[a,1] = -Z[a-1,1]
        mpredN2[a,1] = log_N_a[a,1]
        mcomp[a,1] = log_N_a[a,1] - log_N_a[a-1,1]
        if(a == A){
            mpredN[a,1] = -log(1-exp(-Z[a-1,1]))
        }
    }

    ## Peice for projections
    len_gc = length(given_catch)  
    solved_Fs = numeric(len_gc)
    MMs = list()
    SSs = list()
    NNs = list()
    MM = numeric(A)
    SS = numeric(A)
    NN = numeric(A)
    ##End peice for projections

    
    

    ##Cohort Dynamics
    for(y in 2:Y){

        ##Peice for projections
        ##Are we in a projection year?
        if(y > og_Y){
            if(proj_type == 1){
                for(a in 1:A){
                    MM[a] = M[a,y-1]
                    SS[a] = S_ay[a,y-1]
                    NN[a] = exp(log_N_a[a,y-1])
                }
            
            
            cur_catch = given_catch[y-og_Y]
            nFy = solveF(cur_catch,0.1,SS,MM,NN)
            proj_F = nFy*SS
            solved_Fs[y-og_Y] = nFy
            for(a in 1:A){
                Z[a,y] = proj_F[a]+M[a,y]
                F[a,y] = proj_F[a]
            }
            
            SSs[[y-og_Y]] = SS
            MMs[[y-og_Y]] = MM
            NNs[[y-og_Y]] = NN
            }
        }

        for(a in 2:A){
            mcomp[a,y] = log_N_a[a,y]
            mpredN2[a,y] = log_N_a[a-1,y-1]-Z[a-1,y-1]
            if(a == A){
                mpredN2[a,y] = log(exp(mpredN2[a,y])+exp(log_N_a[a,y-1]-Z[a,y-1]))
            }
            mpredN[a,y] = mpredN2[a,y]
        }
        
    }

    resNraw = matrix(NA,A,Y)
    
    ##Build sigma for the recruits
    rec_sigma = diag(rep(exp(log_recruit_sd)^2,Y))

    ##Build sigma for every other age class
    surv_sd = 0.27 - exp(log_surv_sd)
    
    ev_sigma = diag(c(exp(log_N0_sd)^2,rep(surv_sd^2,Y-1)))
    plus_sigma = ev_sigma
    if(plus_s){
        p_surv_sd = exp(parameters$log_p_surv_sd)
        plus_sigma = diag(c(exp(log_N0_sd)^2,rep(p_surv_sd^2,Y-1)))
    }
    N_nll = numeric(A)

    N_nll[1] = dmvnorm(mcomp[1,],mpredN[1,],rec_sigma,TRUE)
    nll = nll - dmvnorm(mcomp[1,],mpredN[1,],rec_sigma,TRUE)
    resNraw[1,] = mcomp[1,]-mpredN[1,]

    if(plus_s){
        AAA = 18
        BBB = A
    }else{
        AAA = A
    }
    
    for(a in 2:AAA){
        nll = nll - dmvnorm(mcomp[a,],mpredN[a,],ev_sigma,TRUE)
        N_nll[a] = dmvnorm(mcomp[a,],mpredN[a,],ev_sigma,TRUE)
        resNraw[a,] = mcomp[a,]-mpredN[a,]
    }
    if(plus_s){
        for(a in (AAA+1):BBB){
            nll = nll - dmvnorm( mcomp[a,],mpredN[a,],plus_sigma,TRUE)
            N_nll[a] = dmvnorm( mcomp[a,],mpredN[a,],plus_sigma,TRUE)
            resNraw[a,] = mcomp[a,]-mpredN[a,]
        }
    }

    REPORT(N_nll)
    REPORT(rec_sigma)
    REPORT(ev_sigma)
    REPORT(plus_sigma)
    REPORT(mcomp)
    REPORT(mpredN)
    REPORT(resNraw)
    
    log_recruit = log_N_a[1,]
    N = exp(log_N_a) 
    
    ##Adjusting N for different times
    N_mid = matrix(NA,A,Y)
    N_fall = matrix(NA,A,Y)

    for(a in 1:A){
        for(y in 1:Y){
            N_mid[a,y] = N[a,y]*exp(-Z[a,y]*0.5)
            N_fall[a,y] = N[a,y]*exp(-Z[a,y]*fall_adjust)
        }
    }

    ##Convert to length

    NL = matrix(NA,L3,Y)
    NL_mid = matrix(NA,L3,Y)
    NL_fall = matrix(NA,L3,Y)
    
    for(y in 1:Y){
        NL[,y] = plas_jan[[y]]%*%N[,y]
        NL_mid[,y] = plas_mid[[y]]%*%N_mid[,y]
        NL_fall[,y] = plas_fall[[y]]%*%N_fall[,y]
    }

    ##Observation Equations ##

    ##Setting up Qs
    QLM = matrix(NA,L3,2)
    Qmax = exp(log_Qmax)
    engel.Qmax = exp(engel.log_Qmax)
    QL5095 = exp(ordered_inv_transform(c(log_QL50,log_QL95)))
    QL50 = QL5095[1]
    QL95 = QL5095[2]
    engel.QL5095 = exp(ordered_inv_transform(c(engel.log_QL50,engel.log_QL95)))
    engel.QL50 = engel.QL5095[1]
    engel.QL95 = engel.QL5095[2]
    
    ##QL50 = exp(log_QL50)
    ##engel.QL50 = exp(engel.log_QL50)
    ##QL95 = exp(log_QL95)
    ##engel.QL95 = exp(engel.log_QL95)
    QK = -1*log(0.05/0.95)/(QL95-QL50)
    engel.QK = -1*log(0.05/0.95)/(engel.QL95-engel.QL50)

    
    for(l in 1:L3){
        len = len_bins[l]
        QLM[l,1] = (engel.Qmax/(1+exp(-engel.QK*(len-engel.QL50)))^(delta_survey[1]))
        QLM[l,2] = Qmax/(1+exp(-QK*(len-QL50)))^(delta_survey[2])
    }
    Q_rho = log(QLM[,1]/QLM[,2])

    Q_nll = numeric(L3)
    if(Q_prior == TRUE){
        Q_max_sd = exp(log_Q_max_sd)
        for(l in 1:L3){
            Q_nll[l] = dnorm(Q_rho[l],log(Qpriors[l]),Q_max_sd,TRUE)
            nll = nll- dnorm(Q_rho[l],log(Qpriors[l]),Q_max_sd,TRUE)
        }

        REPORT(QLM)
        ADREPORT(Q_max_sd)
    }

    REPORT(Q_nll)
    REPORT(Q_rho)

    ##Survey
    n_years = length(survey_list)

    rhoS = rho_trans(logit_rhoS)
    sd_survey = exp(log_sd_survey)

    log_exp_index = numeric(survey_size)
    exp_index = numeric(survey_size)
    log_survey_resids = numeric(survey_size)
    std_log_survey_resids = numeric(survey_size)
    pearson_survey_resids = numeric(survey_size)
    total_survey_abundance = numeric(n_years)
    
    Ssigmas = list()

    survey_nll = numeric(n_years)
    for(yy in 1:n_years){
        curr_list = survey_list[[yy]]
        ##Since it's off by one
        year = curr_list$year+1
        ##Since it's off by one 
        stype = curr_list$type+1
        projy = curr_list$projy

        s_ind = curr_list$indices
        s_len = curr_list$lengths+1
        s_place = curr_list$place+1
        s_map = curr_list$mmap+1

        nlen = curr_list$nlen

        exp_ind = numeric(length(s_ind))
        c_sd = numeric(length(s_ind))

        y_min_len = min(s_len)
        y_max_len = max(s_len)
        

        for(i in 1:length(s_ind)){
            Q = QLM[s_len[i],stype]
            if(s_len[i] <= y_min_len){
                NLminus = 0
                for(l in 1:y_min_len){
                    Qt = QLM[l,stype]
                    NLminus = NLminus + Qt*NL_fall[l,year]
                }
                exp_ind[i] = NLminus
            }else if(s_len[i] != y_max_len){
                exp_ind[i] = Q*NL_fall[s_len[i],year]
            }else{
                NLplus = 0
                for(l in (y_max_len):L3){
                    QLt = QLM[l,stype]
                    NLplus = NLplus + QLt*NL_fall[l,year]
                }
                exp_ind[i] = NLplus
            }
        
            exp_index[s_place[i]] = exp_ind[i]
            log_exp_index[s_place[i]] = log(exp_ind[i])
            ##Not in a projection
            if(proj_years == 0){
                c_sd[i] = sd_survey[s_map[i]]
            }
            total_survey_abundance[yy] = total_survey_abundance[yy] + exp_ind[i] 
        }

        log_diff = log(s_ind)-log(exp_ind)
        Ssigma = matrix(NA,nlen,nlen)
 
        for(ii in 1:nlen){
            for(jj in 1:ii){
               Ssigma[ii,jj] = (rhoS[stype]^(ii-jj))*c_sd[ii]^2
               Ssigma[jj,ii] = Ssigma[ii,jj]
            }
        }

        if(s_dist_type == "normal"){
            survey_nll[yy] = dmvnorm(log_diff,Sigma=Ssigma,log=TRUE)
            nll = nll - dmvnorm(log_diff,Sigma=Ssigma,log=TRUE)
        }else{
            ##We want t_df > 3
            t_df = exp(parameters$log_t_df)+3
            muuu = rep(0,length(log_diff))
            survey_nll[yy] = dmvt(log_diff,muuu,Ssigma,t_df,TRUE)
            nll = nll - dmvt(log_diff,muuu,Ssigma,t_df,TRUE)
            ##Fix for cov, df/(df-2)*Sigma
            Ssigma = (t_df/(t_df-2))*Ssigma
            ##ADREPORT(t_df)
        }

        
       LinvD = solve(chol(Ssigma))
       std_resids = LinvD%*%log_diff
       
        for(i in 1:length(s_ind)){
            log_survey_resids[s_place[i]] = log_diff[i]
            std_log_survey_resids[s_place[i]] = std_resids[i]
            pearson_survey_resids[s_place[i]] = log_diff[i]/c_sd[i]
        }

        Ssigmas[[yy]] = Ssigma
        
    }

    ##I borke this here
    if(r_proj == 1){
        ttotal_survey_abundance = numeric(Y)
        for(y in 1:n_years){
            ttotal_survey_abundance[y] = total_survey_abundance[y]
        }
        for(y in (n_years+1):Y){
            ttotal_survey_abundance[y] = 0
        }
        
        diffy = Y-n_years
        addspace = nrow(QLM)*diffy
        texp_index = numeric(length(exp_index)+addspace)
        texp_index[1:length(exp_index)] = exp_index
        qq = length(exp_index)+1
        for(jj in (n_years+1):Y){
            for(zz in 1:nrow(QLM)){
                texp_index[qq] = QLM[zz,2]*NL_fall[zz,jj]
                ttotal_survey_abundance[jj] = ttotal_survey_abundance[jj] + texp_index[qq]
                qq = qq+1
            }
        }
        total_survey_abundance = ttotal_survey_abundance
        exp_index = texp_index
        log_exp_index = log(exp_index)
    }

    REPORT(log_survey_resids)
    REPORT(std_log_survey_resids)
    REPORT(pearson_survey_resids)
    REPORT(survey_nll)
    REPORT(log_exp_index)

    ## Catch at age and length
    catch_at_age = matrix(NA,A,Y)
    for(a in 1:A){
        for(y in 1:Y){
            catch_at_age[a,y] = N_mid[a,y]*(1-exp(-Z[a,y]))*(F[a,y]/Z[a,y])
        }
    }

    REPORT(catch_at_age)

    catch_at_length = matrix(NA,L3,Y)
    for(y in 1:Y){
        catch_at_length[,y] = plas_mid[[y]]%*%catch_at_age[,y]
    }

    REPORT(catch_at_length)

    catch_biomass_mat = weightsF*catch_at_length
    expected_landings = colSums(catch_biomass_mat)
    log_expected_landings = log(expected_landings)

    REPORT(catch_biomass_mat)
    REPORT(log_expected_landings)
    
    ##handling landings
    log_landing_resids = numeric(length(landing_nums))
    std_landing_resids = numeric(length(landing_nums))
    landings_sd = exp(log_landings_sd)

    landing_nll = numeric(length(landing_nums))
    for(i in 1:length(landing_nums)){
        log_landing_resids[i] = log(landing_nums[i])-log_expected_landings[i]
        std_landing_resids[i] = log_landing_resids[i]/landings_sd
        landing_nll[i] = censored_bounds(log(landing_nums[i]),log_expected_landings[i],landings_sd,-log_landingsL[i],log_landingsU[i])
        nll = nll - censored_bounds(log(landing_nums[i]),log_expected_landings[i],landings_sd,-log_landingsL[i],log_landingsU[i])
    }
    REPORT(landing_nll)
    REPORT(log_landing_resids)
    REPORT(std_landing_resids)

    ## handling the catch proportion 
    n_years2 = length(catch_list)

    obs_crls = list()
    exp_props = list()
    diffs = list()
    sigmas = list()
    agg_catches = list()
    years = numeric(n_years)
    keys = list()
    props = list()
    agg_props = list()
    std_diffs = list()
    logit_exp_props = list()

    rhoC = rho_trans(logit_rhoC)

    catch_nll = numeric(n_years2)
    for(i in 1:n_years2){
        curr_list = catch_list[[i]]
        year = curr_list$year+1
        ysize = curr_list$ysize
        key = curr_list$key+1
        prop = curr_list$prop
        c_map = curr_list$mmap+1

        years[[i]] = year
        keys[[i]] = key
        props[[i]] = prop
        
        curr_c_at_l = catch_at_length[,year]
        agg_prop = aggregate_catch(prop,key,ysize)
        obs_crl = make_CRL_vec(agg_prop)
        
        
        agg_c_at_l = aggregate_catch(curr_c_at_l,key,ysize)

        exp_prop = agg_c_at_l/sum(agg_c_at_l)
        ##logit_exp_prop = logit(exp_prop)
        exp_crl = make_CRL_vec(exp_prop)

        diff = obs_crl-exp_crl

        nlen = length(diff)
        
        log_sd_catch_prop = log_sd_catch_prop_m[c_map]
        sd_catch_prop = exp(log_sd_catch_prop)
        
        sigma = matrix(NA,nlen,nlen)
        for(ii in 1:nlen){
            for(jj in 1:ii){
               sigma[ii,jj] = (rhoC^(ii-jj))*sd_catch_prop[ii]^2
               sigma[jj,ii] = sigma[ii,jj]
            }
        }

        LinvD = solve(chol(sigma))
        std_diff = LinvD%*%diff


        
        obs_crls[[i]] = obs_crl
        exp_props[[i]] = exp_prop
        logit_exp_props[[i]] = qlogis(exp_prop) 
        diffs[[i]] = diff
        sigmas[[i]] = sigma
        agg_catches[[i]] = agg_c_at_l
        agg_props[[i]] = agg_prop
        std_diffs[[i]] = std_diff

        if(c_dist_type == "normal"){
            catch_nll[i] = dmvnorm(diff,Sigma=sigma,log=TRUE)
            nll = nll - dmvnorm(diff,Sigma=sigma,log=TRUE)
        }else{
            log_c_t_df = parameters$log_c_t_df
            c_t_df = exp(log_c_t_df)+3
            muuu = rep(0,length(diff))
            catch_nll[i] = dmvt(diff,muuu,sigma,c_t_df,TRUE)
            nll = nll - dmvt(diff,muuu,sigma,c_t_df,TRUE)
            ##Fix for cov, df/(df-2)*Sigma
            sigma = (c_t_df/(c_t_df-2))*sigma
            sigmas[[i]] = sigma
            ##ADREPORT(t_df)
   
        }
    }

    REPORT(catch_nll)
    biomass_mat = weightsF*NL
    ssb_mat = maturityF*biomass_mat
    
    tot_ssb = colSums(ssb_mat)
    log_tot_ssb = log(tot_ssb)
    

    REPORT(agg_props)
    REPORT(std_diffs)
    REPORT(exp_props)
    REPORT(logit_exp_props)
    REPORT(diffs)
    REPORT(sigmas)
    REPORT(rhoS)
    REPORT(Ssigmas)
    REPORT(NL)
    REPORT(NL_mid)
    REPORT(NL_fall)
    REPORT(N_mid)
    REPORT(N_fall)
    REPORT(N)
    REPORT(rec_sigma)
    REPORT(ev_sigma)
    REPORT(plus_sigma)
    REPORT(F)
    REPORT(Z)
    REPORT(S_ay)
    REPORT(S_ly)
    REPORT(ages_jan)
    REPORT(len_sds)
    REPORT(len_mus)
    REPORT(plas_jan)
    REPORT(plas_mid)
    REPORT(plas_fall)

    ##ADREPORTing
    ell = ell_L[1]
    L = ell_L[2]
    k = gf_par$k
    s = gf_par$s
    S = gf_par$S
    Fy_sd = exp(log_Fy_sd)
    recruit_sd = exp(log_recruit_sd)
    catch_sds = exp(log_sd_catch_prop_m)
    log_total_survey_abundance = log(total_survey_abundance)
    log_recruit = log_N_a[1,]
    tot_ssb = exp(log_tot_ssb)

    ##Calculate Survey < 15cm abundance
    survey15 = numeric(n_years)
    for(yy in 1:n_years){
        curr_list = survey_list[[yy]]
        ##Since it's off by one
        year = curr_list$year+1
        ##Since it's off by one 
        stype = curr_list$type+1
        for(l in 1:15){
            survey15[yy] = survey15[yy] + QLM[l,stype]*NL[l,yy]
        }
    }
    log_survey15 = log(survey15)
    

    ##Add on A for reporting since it's a deviation 
    as_pgA = as_pg+A
    init_a_pg = as_pgA[1]
    
    ##Old school VB parameters
    K = -log(k)
    L_inf_top = L-ell*k^(A-1)
    L_inf_bot = 1-k^(A-1)
    L_inf = L_inf_top/L_inf_bot

    t_z_p = (L-ell)/(L-ell*k^(A-1))
    t_zero = 1-(1/log(k))*log(t_z_p)

    ##Needed for eaiser report generation
    exp_prop_key = lapply(1:length(logit_exp_props),function(x){
        y = rep(x,length(logit_exp_props[[x]]))
        y
    })
    REPORT(exp_prop_key)

    ##This is really gross but I can't seem to use unlist, RTMB can't ADREPORT a list...
    prop_count = 0
    for(i in 1:length(logit_exp_props)){
        for(j in 1:length(logit_exp_props[[i]])){
            prop_count = prop_count+1
        }
    }
    l_exp_props = numeric(prop_count)
    jj = 0
    for(i in 1:length(logit_exp_props)){
        for(j in 1:length(logit_exp_props[[i]])){
            jj = jj+1
            l_exp_props[jj] = logit_exp_props[[i]][j]
        }
    }

    Fbar = colSums(F)/A
    log_Fbar = log(Fbar)
    
    REPORT(delta_survey)

    ##Parameters
    ADREPORT(ell)
    ADREPORT(L)
    ADREPORT(k)
    ADREPORT(s)
    ADREPORT(S)
    ADREPORT(b_beta1)
    ADREPORT(r_b_beta2)
    ADREPORT(sel_shape)
    ADREPORT(sel_scale)
    
    ADREPORT(Fy_sd)
    ADREPORT(surv_sd)
    ADREPORT(N0_sd)
    ADREPORT(recruit_sd)
    ADREPORT(rhoS)
    ADREPORT(sd_survey)
    ADREPORT(landings_sd)
    ADREPORT(rhoC)
    ADREPORT(catch_sds)
    ADREPORT(init_a_pg)
    ADREPORT(delta_survey)
    ADREPORT(log_survey15)
    if(s_dist_type != "normal"){
        ADREPORT(t_df)
    }

    if(c_dist_type != "normal"){
        ADREPORT(c_t_df)
    }

    if(plus_s){
        ADREPORT(p_surv_sd)
    }
    
    ##VB 'Parameters'
    ADREPORT(K)
    ADREPORT(L_inf)
    ADREPORT(t_zero)

    ##Quantities of interest
    ADREPORT(log_total_survey_abundance)
    ADREPORT(log_recruit)
    ADREPORT(log_exp_index)
    ADREPORT(log_expected_landings)
    ADREPORT(log_tot_ssb)
    ADREPORT(tot_ssb)
    ADREPORT(Q_rho)
    ADREPORT(QLM)
    ADREPORT(l_exp_props)
    ADREPORT(log_Fbar)
    ADREPORT(Fbar)
    
    
    
    nll
}


run_projections <- function(tmb.data,sdreport,map,proj_years,hessian,weights,maturity,randoms,old_report,old_opt,supplied_F=NULL,given_catch=NULL){
    tdat = tmb.data
    gparms = as.list(sdreport,"Estimate")
    ttmap = map
    tdat$gamma_key = c(tdat$gamma_key,rep(tdat$gamma_key[tdat$Y],proj_years))


    if(is.null(supplied_F) & is.null(given_catch)){

        oldNalen = length(gparms$log_N_a)
        lastCol = gparms$log_N_a[,ncol(gparms$log_N_a)]
        gparms$log_N_a = cbind(gparms$log_N_a,matrix(rep(lastCol,proj_years),nrow=tmb.data$A,proj_years))

        supplied_F = rep(gparms$log_Fy[tmb.data$Y],proj_years)
        goose = as.character(ttmap$log_Fy)
        goose[(tmb.data$Y+1):((tmb.data$Y+proj_years))] = (tmb.data$Y+1):((tmb.data$Y+proj_years))
        goose = as.factor(goose)
        ttmap$log_Fy = goose

        gparms$log_Fy = c(gparms$log_Fy,supplied_F)
        tdat$orig_end_year = tdat$end_year
        tdat$end_year = tdat$end_year+proj_years
        tdat$Y = tdat$Y+proj_years
        ##Yes this is on purpose
        tdat$og_Y = tdat$Y
        tdat$maturityF = cbind(tdat$maturityF,maturity)
        tdat$weightsF = cbind(tdat$weightsF,weights)
        tdat$r_proj = 1


        ##Make sure the model gets the right data
        rram_wrapper_wrap <- function(dat){
            rram_wrap <- function(parameters){
                rram_model(parameters,dat)
            }
        }

        rram_to_run <- rram_wrapper_wrap(tdat)    

        obj2 = RTMB::MakeADFun(rram_to_run,gparms,random=randoms,map=ttmap)
        obj2$fn(old_opt$par)
        rep2 = obj2$report()
        sdr2 = RTMB::sdreport(obj2,hessian.fixed=hessian)

        ret = list(obj=obj2,report=rep2,sdr=sdr2,ssdr=summary(sdr2),tmb.data=tdat)
        
        
    }else{

        oldNalen = length(gparms$log_N_a)
        lastCol = gparms$log_N_a[,ncol(gparms$log_N_a)]
        gparms$log_N_a = cbind(gparms$log_N_a,matrix(rep(lastCol,proj_years),nrow=tmb.data$A,proj_years))

        tdat$orig_end_year = tdat$end_year
        tdat$end_year = tdat$end_year+proj_years
        tdat$og_Y = tmb.data$Y
        ##Again on purpose
        tdat$Y = tdat$Y+proj_years
        tdat$maturityF = cbind(tdat$maturityF,maturity)
        tdat$weightsF = cbind(tdat$weightsF,weights)
        tdat$r_proj = 1
        if(!is.null(supplied_F)){
            tdat$proj_type = 0
            tdat$supplied_F = supplied_F
            randos = randoms
        }else{
            tdat$proj_type = 1
            tdat$given_catch = given_catch
            randos = randoms
        }


        rram_wrapper_wrap <- function(dat){
            rram_wrap <- function(parameters){
                rram_model(parameters,dat)
            }
        }

        rram_to_run <- rram_wrapper_wrap(tdat)    

        obj2 = RTMB::MakeADFun(rram_to_run,gparms,random=randoms,map=ttmap)
        obj2$fn(old_opt$par)
        rep2 = obj2$report()
        sdr2 = RTMB::sdreport(obj2,hessian.fixed=hessian)
        ssdr2 = summary(sdr2)
        
        ret = list(obj=obj2,report=rep2,sdr=sdr2,ssdr=ssdr2,tmb.data=tdat)
        
        

    }

    ##Check if random effects have varied greatly from the model fit and if so throw a waring
    proFy = ret$report$log_Fy[1:tmb.data$Y]
    modFy = old_report$log_Fy
    FyC = all.equal(proFy,modFy)
    if(FyC != TRUE){
        warning("F Random Effects have changed from model fit!")
    }

    proN = ret$report$N[,1:tmb.data$Y]
    modN = old_report$N
    NyC = all.equal(proN,modN)
    if(NyC != TRUE){
        warning("N Random Effects have changed from model fit!")
    }
    
    

ret 
}
