library("simulation.utils")

parameter_list <- list()

parameter_list$envt_change_1 <- c(0.1,0.5,0.9)
parameter_list$envt_change_2 <- 0.5
parameter_list$maternal_resources_1 <- 20
#... to be completed
parameter_list$survival_strength_envt1 <- c(0.1,0.5,1,2)
parameter_list$survival_strength_envt2 <- 0.5
parameter_list$max_gen <- 10000
   

make.batch.file(
        parameter_list=parameter_list
        ,executable_path = "./coop_size.exe"
        ,n_replicates = 10)

