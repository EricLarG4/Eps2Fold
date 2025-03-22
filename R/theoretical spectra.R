#Dependencies----

#librarian::shelf(tidyverse, data.table, tictoc, readxl)

library('tidyverse')
library('data.table')
library('readxl')

#UV parameters database----

#parameters for the extinction coefficient of oligodeoxynucleotides at 260 nm
#source: https://doi.org/10.1016/j.bpc.2007.12.004
nn.260 <- fread('data/db/epsilon260.csv')

#parameter ratio for ssDNA
nn.param <- fread('data/db/epsilondb.csv') %>% 
  pivot_longer(cols = 2:ncol(.),
                values_to = "ratio",
                names_to = "nn") %>% 
  as.data.table()


#Contributions of pairs of nt at 260 nm----
dt.contributR <- function(input.seq){
  
  #initialisation with first and last nucleotides alone
  epsilon.calc <- data.table(
    position = c(1, str_length(input.seq)), #first and last position
    nn = c(substr(input.seq, 1, 1), #first and last nt name
           substr(input.seq, str_length(input.seq), str_length(input.seq)))
  )
  
  #extraction of all nt pairs position and name
  for (i in 1:(str_length(input.seq) - 1)) {
    buffer <- data.table(position = paste(i, i + 1, sep = '-'),
                         nn = substr(input.seq, i, i + 1)
    )
    epsilon.calc <- rbindlist(list(epsilon.calc, buffer)) 
  }
  
  #association of each nt pair with its 260 nm parameter from db
  epsilon.calc[nn.260, on = 'nn', contrib := epsij*1000]
  
  epsilon.calc
}


#Contributions of pairs of nt at all wavelength---- 

#Determined from contribution at 260 nm (dt.contributR()) and ratios from db (nn.param)

dt.spec.calcR <- function(input.wl, input.contrib){
  
  #Select wavelength in db
  nn.param <- nn.param[wl == input.wl]
  
  #Join ratio from db to each pair or nt
  input.contrib[nn.param, on = .(nn), ratio := ratio]
  
  #Calculate the contribution of each pair of nt from ratio
  input.contrib[,by = position, ratioed.contrib := ratio*contrib]
  
  #Sum across sequences to get epsilon
  eps <- input.contrib[, sum(ratioed.contrib)]
  
}


# # #Example----
# # 
# ##Example data----
# oligo.table <- data.frame(
#   oligo = rep(c("22AG", "24TTG"), 310-220+1),
#   seq = rep(c("AGGGTTAGGGTTAGGGTTAGGG", "TTGGGTTAGGGTTAGGGTTAGGGA"), 310-220+1),
#   wl = 220:310
# )


# ##Example calculation----
# tic()
# oligo.spec <- oligo.table %>%
#   group_by(seq, wl) %>% 
#   mutate(
#     eps = dt.spec.calcR(input.contrib = dt.contributR(input.seq = seq), input.wl = wl)
#   )
# toc()


# ##Example plot----
# ggplot(oligo.spec, aes(x = wl, y = eps, color = oligo)) +
#   geom_line()
