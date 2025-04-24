library(tidyverse)

##Setup landing bounds, need to provide values all the way back to '59
b1L = 0.95
b1U = 1.25
landing_b1L = rep(b1L,62)
landing_b1U = rep(b1U,62)

##Load data and maps and things
survey = read.csv("unconvertedSurvey.csv")
landings = read.csv("landings.csv")
landings$lower_bound = landing_b1L
landings$upper_bound = landing_b1U

weight_array = readRDS("lweights.rds")
maturity_array = readRDS("lmaturity.rds")
catch_stuff = readRDS("catch_and_key.rds")

##Set things for TMB map
tmap = list(log_Qmax =as.factor(NA),log_S=as.factor(1),log_surv_sd=as.factor(NA),log_delta_survey=as.factor(c(1,NA)),log_t_df=as.factor(1))


##Controls the lengths used in the survey each year 
survey_l_key = data.frame(survey.year=1983:2020,min_length=7,max_length=37)
survey_l_key$min_length[survey_l_key$survey.year < 1995] = 7


##The example isn't using cmap though...
cmap = readRDS("cmap.rds")


##Source file to make data and parameters
source("makeData.R")



### Create the map for the survey SDs
d_and_p = build_data_and_parameters(weight_array,maturity_array,survey,
                          landings,0.05,catch_stuff$prop_catch,catch_stuff$agg_key,
                          years=1983:2020,ages=1:20,survey_l_key=survey_l_key,
                          tmb.map=tmap,survey_sd_map = NULL,catch_prop_map = NULL,rounding_bit = 0.05,gf_ext=TRUE,sel_type = "old",l_dist="normal",start.parms=start.parms,s_dist="goof")

##Get the lengths in each year
lens = lapply(1:length(d_and_p$tmb.data$survey_list),function(x){
    yy = d_and_p$tmb.data$survey_list[[x]]
    df = data.frame(year=x,length=yy$lengths)
    df})

##setup how you want sds for length
case_sdlen_type <- function(year,length){
    ##goose = rep(0,length(year))
     goose = case_when(length <= 10 ~ 0,
               length > 10 & length < 36 ~ 1,
               length >= 36 ~ 2)
    if(year[1] > 12){
     #   goose = rep(1,length(year))
        goose = case_when(
            length <= 10 ~ 3,
            length > 10 & length < 36 ~ 4,
                         length >= 36 ~ 5)
    }  
    goose
}



##make map 
neomap = lapply(lens,function(x){
    case_sdlen_type(x$year,x$length)})

##Build d_and_p with parameters
d_and_p = build_data_and_parameters(weight_array,maturity_array,survey,
                          landings,0.05,catch_stuff$prop_catch,catch_stuff$agg_key,
                          years=1983:2020,ages=1:20,survey_l_key=survey_l_key,
                          tmb.map=tmap,survey_sd_map = neomap,catch_prop_map = NULL,rounding_bit = 0.05,gf_ext=TRUE,sel_type = "old",l_dist="normal",start.parms=start.parms,s_dist="normal")



##Source the model file
A = d_and_p$tmb.data$A
source("RRAMmodelOL.R")

dat = d_and_p$tmb.data

##Make sure the model gets the right data
rram_wrapper_wrap <- function(dat){
    rram_wrap <- function(parameters){
        rram_model(parameters,dat)
    }
}

rram_to_run <- rram_wrapper_wrap(d_and_p$tmb.data)    


saveRDS(neomap,"neomap.rds")

obj = MakeADFun(rram_to_run,d_and_p$parameters,random=c("log_N_a","log_Fy"),map=d_and_p$map)
opt = nlminb(obj$par,obj$fn,obj$gr,control=list(iter.max=2000,eval.max=2000,trace=FALSE))
opt = nlminb(opt$par,obj$fn,obj$gr,control=list(rel.tol=1e-8))
##Get the hessian so we can use it for projections and sdreport and only call it once
hess = optimHess(opt$par,obj$fn,obj$gr)



repp = obj$report()

sdr = sdreport(obj,hessian.fixed = hess)
ssdr = summary(sdr)

## Put things in a way to be used by the reporting system
tmb.data = d_and_p$tmb.data
modDat = d_and_p$mod.data
orig_data = list(years = 1983:2020,
                 ages = 1:20,
                 length= 7:45,
                 agg_key=catch_stuff$agg_key)
outdat = list(opt=opt,ssdr=ssdr,report=repp,orig_data=orig_data)

##Source the file with plotting, report gen and some utilities code
##Uses plotly + ggplot
source("utilities.R")


##create the report!
create_report("test2",outdat,"./",tmb.data,modDat)

## Do projections


##Assuming projection weights and maturity are just the average of the last 5 years of data...
##M is just whatever the model was doing

proj_years = 5
weightsP = matrix(rowMeans(tmb.data$weightsF[,c("2016","2017","2018","2019","2020")]),ncol=proj_years,nrow=60)
maturityP = matrix(rowMeans(tmb.data$maturityF[,c("2016","2017","2018","2019","2020")]),ncol=proj_years,nrow=60)

##Add a given catch amount
proj1 = run_projections(tmb.data,sdr,d_and_p$map,proj_years,hess,weightsP,maturityP,c("log_N_a","log_Fy"),repp,opt,given_catch = c(200,50,200,75,65))

##Given specific log F_y
proj2 = run_projections(tmb.data,sdr,d_and_p$map,proj_years,hess,weightsP,maturityP,c("log_N_a","log_Fy"),repp,opt,supplied_F = c(log(0.1),log(0.2),log(0.1),log(0.1),log(0.3)))

##Project the terminal F
proj3 = run_projections(tmb.data,sdr,d_and_p$map,proj_years,hess,weightsP,maturityP,c("log_N_a","log_Fy"),repp,opt)



create_projection_report("test_pro",proj1,"./",proj1$tmb.data,modDat)
create_projection_report("test_pro2",proj2,"./",proj2$tmb.data,modDat)
create_projection_report("test_pro3",proj3,"./",proj2$tmb.data,modDat)





my_retros = run_peels(d_and_p$orig_data)



create_retro_report("test_retro","./",my_retros)
