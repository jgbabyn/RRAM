#' Generate prior for conv. factor
#'
#' This uses the inverse of the conversion factor between
#' Engel and Campelen Trawl from Warren 1997
#' @param x the value of the length
#' @export
redfish_rat <- function(x){
    lny = 6.7580137+0.006839*x-1.927210*log(x)
    1/exp(lny)
}


#' Transform monotone increasing vector to optimizer scale
#'
#' @param x the vector to transform
#' @export
ordered_transform <- function(x){
    y = numeric(length(x))
    y[1] = x[1]
    for(i in 2:length(x)){
        y[i] = log(x[i]-x[i-1])
    }
    y
}


#' Build model data and parameters and map!
#' @param weight_array array of weights with dimnames for male/female, year, length, etc.
#' @param maturity_array array of maturity ogives to use, with appropriate dimnames
#' @param survey_df data.frame of survey data to use
#' @param landings_df data.frame of landings data as well as upper and lower bounds
#' @param base_M the base_M underlying value assumption
#' @param upper_landing_bounds column of landings upper multiplier for bounds
#' @param lower_landing_bounds column of landings lower multiplier for bounds
#' @param catch_prop matrix with named row/cols of catch proportions
#' @param agg_key specifies how to aggregate catch
#' @param years vector of years to include
#' @param ages vector of ages to use in the model
#' @param lengths vector of lengths to include
#' @param tmb.map optional named list, use to manually provide TMB map
#' @param random optional, override the default random effect parameters
#' @param start.parms optional named list of start parameters, used to override the starting parameters
#' @param data optional named list of data to give the TMB model, used to directly override the data
#' @param inf_length the max length used in the model
#' @param Q_prior_max the maximum of the Q prior length to use
#' @param pg_ext controls how fine the plus group growth extension goes
#' @param rounding_bit how much 95% selectivity for catch is ahead of 50%
#' @param survey_sd_map list by year and length specifying how to internally map the survey sds
#' @param catch_prop_map list by year and length specifying how to internally map the catch prop sds
#' @export
build_data_and_parameters <- function(weight_array,maturity_array,survey_df,landings_df,base_M,
                            catch_prop,
                            agg_key,years=1983:2021,ages=1:20,lengths=7:45,tmb.map=NULL,random=NULL,start.parms=NULL,
                            data=NULL,
                            inf_length=60
                           ,Q_prior_max=35,pg_ext=60,rounding_bit=0.01,survey_sd_map = NULL,catch_prop_map=NULL){

    ##orginal data for retros
    orig_data = list()
    orig_data$weight_array=weight_array
    orig_data$maturity_array = maturity_array
    orig_data$survey_df = survey_df
    orig_data$landings_df = landings_df
    orig_data$base_M = base_M
    orig_data$catch_prop = catch_prop
    orig_data$agg_key = agg_key
    orig_data$years = years
    orig_data$ages = ages
    orig_data$lengths = lengths
    if(!missing(tmb.map)){
        orig_data$tmb.map = tmb.map
    }else{
        orig_data$tmb.map = NULL
    }
    if(!missing(start.parms)){
        orig_data$start.parms = start.parms
    }else{
        orig_data$start.parms = NULL
    }
    if(!missing(data)){
        orig_data$data = data
    }else{
        orig_data$data = NULL
    }
    orig_data$inf_length = inf_length
    orig_data$Q_prior_max = Q_prior_max
    orig_data$pg_ext = pg_ext
    orig_data$rounding_bit = rounding_bit
    if(!missing(survey_sd_map)){
        orig_data$survey_sd_map = survey_sd_map
    }else{
        orig_data$survey_sd_map = NULL
    }
    if(!missing(catch_prop_map)){
        orig_data$catch_prop_map = catch_prop_map
    }else{
        orig_data$catch_prop_map = NULL
    }

    bin_adjust = 0.5
    Y = length(seq(years[1],years[length(years)]))
    A = length(ages)
    L = length(lengths)
    L3 =length(1:inf_length)

    start_length = lengths[1]
    end_length = lengths[length(lengths)]

    start_age = ages[1]
    end_age = ages[length(ages)]

    start_year = years[1]
    end_year = years[length(years)]

    weightsF = weight_array[as.character(1:inf_length),as.character(start_year:end_year),"female"]
    weightsM = weight_array[as.character(1:inf_length),as.character(start_year:end_year),"male"]

    maturityF = maturity_array[as.character(1:inf_length),as.character(start_year:end_year),"female"]
    maturityM = maturity_array[as.character(1:inf_length),as.character(start_year:end_year),"male"]

        survey2 = survey_df |>
        dplyr::filter(survey.year >= start_year) |>
        dplyr::filter(survey.year <= end_year) |>
        dplyr::mutate(length = dplyr::case_when(length <= start_length ~ start_length,
                                  length >= end_length ~ end_length,
                                  TRUE ~ length)) |>
        dplyr::group_by(survey.year,length) |>
        dplyr::mutate(total=sum(total),total.lcl=sum(total.lcl),total.ucl=sum(total.ucl),
               mean=sum(mean),mean.lcl=sum(mean.lcl),mean.ucl=sum(mean.ucl)) |>
        dplyr::distinct(survey.year,length,.keep_all=TRUE) |>
        dplyr::filter(total > 0)

    add.index <- function(data,scale=1){
        data$rsu = data$sample.units/mean(unique(data$sample.units))
        data$index = data$total/(scale*data$rsu)
        data
    }

    survey3 = add.index(survey2,1e6)
    
    survey_index = survey3$index
    survey_year = as.integer(survey3$survey.year-start_year)
    survey_length = as.integer(survey3$length-1)
    survey_type = ifelse(survey3$survey.year < 1995,0,1)

    landings_df = landings_df |>
        dplyr::filter(Years >= start_year) |>
        dplyr::filter(Years <= end_year)

    log_landingsL = log(landings_df$lower_bound)
    log_landingsU = log(landings_df$upper_bound)
    

    landing_year = as.integer(landings_df$Years-start_year)
    og_landings = landings_df$Total
    adj_landings = ifelse(is.na(landings_df$landadj),landings_df$Total,landings_df$landadj)
    discards = landings_df$discards

    landing_nums = og_landings
    


    agesindices = matrix(1:A,nrow=A,ncol=Y)
    vagesindices = as.vector(agesindices)
    yearsindices = matrix(1:Y,nrow=A,ncol=Y,byrow=TRUE)
    vyearsindices = as.vector(yearsindices)

    kages = matrix(agesindices,nrow=A*Y,ncol=A*Y,byrow=FALSE)
    kyears = matrix(yearsindices,nrow=A*Y,ncol=A*Y,byrow=FALSE)

    

    tmb.data = list(
        start_length = start_length,
        end_length = end_length,
        start_year = start_year,
        end_year = end_year,
        Y = Y,
        A = A,
        L = L,
        L2 = length(start_length:inf_length),
        L3 = L3,
        inf_length = inf_length,
        weightsF = weightsF,
        weightsM = weightsM,
        maturityF = maturityF,
        maturityM = maturityM,
        survey_index = survey_index,
        survey_year = survey_year,
        survey_length = survey_length,
        survey_type = survey_type,
        landing_year = landing_year,
        landing_nums = landing_nums,
        bin_adjust = bin_adjust,
        adj_landings = adj_landings,
        discards = discards,
        log_landingsL = log_landingsL,
        log_landingsU = log_landingsU,
        dum_indices = ifelse(is.na(discards),0,1),
        Qpriors = c(redfish_rat(1:(Q_prior_max-1)),rep(redfish_rat(Q_prior_max),length(Q_prior_max:inf_length))),
        vagesindices = vagesindices-1,
        vyearsindices = vyearsindices-1,
        kages = kages-1,
        kyears = kyears-1
    )

    modDat = list(
        survey=survey3,
        landings=landings_df,
        discards = landings_df,
        weightsF=weightsF,
        weightsM=weightsM,
        maturityF=maturityF,
        maturityM=maturityM
    )
    
    catch_list = list()
    if(!is.null(agg_key) & !is.null(catch_prop)){
        for(i in 1:ncol(agg_key)){
            temp_key = list()
            temp_key$year = as.numeric(colnames(agg_key)[i])-tmb.data$start_year
            temp_key$key = seq(1:tmb.data$L3)
            temp_key$key[temp_key$key < agg_key[1,i]] = agg_key[1,i]
            temp_key$key[temp_key$key > agg_key[2,i]] = agg_key[2,i]
            temp_key$key = temp_key$key-min(temp_key$key)
            temp_key$ysize = length(unique(temp_key$key))
            if(!is.null(catch_prop_map)){
                temp_key$mmap = catch_prop_map[[i]]
            }else{
                temp_key$mmap = rep(0,temp_key$ysize)
                names(temp_key$mmap) = agg_key[1,i]:agg_key[2,i]
            }
            propt = numeric(length(temp_key$key))
            propt[as.numeric(rownames(catch_prop))] = catch_prop[,colnames(agg_key)[i],drop=TRUE]
            propt2 = numeric(inf_length)
            propt2[1:inf_length] = propt[1:inf_length]
            propt2[inf_length] = sum(propt[inf_length:length(propt)])
            temp_key$prop = propt2
            catch_list[[i]] = temp_key
            names(catch_list)[i] = colnames(agg_key)[i]
        }
        ##Check if years in years
        CLyears = colnames(agg_key)[colnames(agg_key) %in% years]
        catch_list = catch_list[CLyears]
        
    }
    tmb.data$catch_list = catch_list

    mapp = list()
    parms = list()
    randos = c("log_N_a","log_F_y")

    tmb.data$base_M = base_M

    tmb.data$Q_prior = 1
    parms$log_Q_max_sd = log(0.5)

    mapp$log_Q_max_sd = as.factor(NA)

      ell_L = c(7,45)
    tran_ell_L = ordered_transform(log(ell_L))

    s_S = c(1,2)

    
    parms$log_ell = tran_ell_L[1]
    parms$log_L = tran_ell_L[2]
    parms$log_k = log(0.5)
    parms$log_S = s_S[2]
    parms$log_s = s_S[1]
    parms$log_recruit_sd = log(0.3)
    
    
    parms$log_N0_sd = log(0.3)
    parms$log_Fy = rep(log(0.2),tmb.data$Y)
    parms$log_Fy_sd = log(0.1)
    parms$log_Qmax = log(1)
    parms$log_QL50 = log(12)
    parms$log_QL95 = log(20)
 
    parms$log_landings_sd = log(0.02)

    mapp$log_QL95 = as.factor(NA)
    mapp$log_QL50 = as.factor(NA)
    mapp$log_N0_sd = as.factor(NA)

    mapp$log_landings_sd = as.factor(NA)

    
    parms$engel.log_Qmax = log(1)
    parms$engel.log_QL50 = log(12)
    parms$engel.log_QL95 = log(20)

    parms$log_N_a = matrix(log(5),nrow=tmb.data$A,ncol=tmb.data$Y)
    parms$log_N_a[tmb.data$A,] = log(4.9)
    parms$log_N_a[1,] = seq(log(5),log(5.1),length.out = tmb.data$Y)
    parms$log_surv_sd = log(0.20)
    mapp$log_surv_sd = as.factor(NA)

    tmb.data$mora_year = 1997-min(years)

    parms$log_b_beta1 = log(25)
    tmb.data$rounding_bit = rounding_bit

    parms$log_sel_scale = log(1)
    parms$log_sel_shape = log(15)

    parms$log_delta_survey = c(log(1),log(1))
    mapp$log_delta_survey = as.factor(c(1,1))
    parms$log_delta_catch = log(1)

    tmb.data$catch_type = 1

    catch_map = lapply(catch_list,function(x){
        x$mmap
    })
    catch_map = unlist(catch_map)
    
    parms$log_sd_catch_prop_m = rep(log(1),length(unique(catch_map)))
    tkeys = do.call(cbind,lapply(catch_list,function(x){x$key}))

    tmb.data$catch_sig_sw = 1
    parms$logit_rhoC = plogis(0.5)

    tmb.data$pg_ext = pg_ext
    parms$log_init_a_pg = log(pg_ext-pg_ext/2)

       ##Creating the survey_list
    ## This makes life easier for AR survey indices

    s_ind = split(tmb.data$survey_index,tmb.data$survey_year)
    s_len = split(tmb.data$survey_length,tmb.data$survey_year)
    s_type = split(tmb.data$survey_type,tmb.data$survey_year)
    s_place = split(1:length(tmb.data$survey_index),tmb.data$survey_year)
    survey_list = list()
    for(i in 1:length(s_ind)){
        tlist = list()
        tlist$year = as.numeric(names(s_ind)[i])
        tlist$lengths = s_len[[i]]
        tlist$indices = s_ind[[i]]
        tlist$type = unique(s_type[[i]])
        tlist$place = s_place[[i]]-1
        tlist$nlen = length(tlist$indices)
        tlist$projy = 0
        if(is.null(survey_sd_map)){
            tlist$mmap = rep(0,length(s_ind[[i]]))
        }else{
            tlist$mmap = survey_sd_map[[i]]
        }
        
        survey_list[[i]] = tlist
    }
    tmb.data$survey_list = survey_list
    ## This is just the overall size of the survey indices because we need to recreate some vectors to not mess up plotting
    tmb.data$survey_size = length(tmb.data$survey_index)
    sur_map = lapply(survey_list,function(x){
        x$mmap
    })
    sur_map = unlist(sur_map)
    parms$logit_rhoS = 0.1
    
    tmb.data$r_proj = 0
    tmb.data$landing_proj_y = rep(0,length(tmb.data$landing_nums))
    tmb.data$cfg = list(sparse=FALSE,trace=TRUE)
    tmb.data$proj_type = 0
    tmb.data$proj_years = 0
    tmb.data$og_Y = tmb.data$Y
    tmb.data$supplied_F = rep(0.1,tmb.data$proj_years)


    ##parms$log_sd_survey = log(0.5)
    parms$log_sd_survey = rep(log(0.5),length(unique(sur_map)))
    

    if(!is.null(start.parms)){
        parms[names(start.parms)] = start.parms
    }

    if(!is.null(tmb.map)){
        mapp[names(tmb.map)] = tmb.map
    }

    if(!is.null(random)){
        randos = random
    }

        if(!is.null(data)){
        tmb.data[names(data)] = data
    }

    
    ##Used possibly for projections
    tmb.data$given_catch = 1

    
    ret = list(tmb.data=tmb.data,mod.data=modDat,parameters=parms,map=mapp,orig_data=orig_data)
    ret

}
