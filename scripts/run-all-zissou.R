library(tidyverse)
library(here)

# go into each source file and change as needed. will add in control file at some point.
# files were separated in order to allow for cloud computing of indidual components as needed

message("Prepare to wait a loooooong time")
source(here("scripts","run-zissou.R"))
message("finished estimating channel island effects")
source(here("scripts","sim-mpa-outcomes.R"))
message("finished simulating MPA outcomes")
source(here("scripts","sim-did.R"))
message("finished simulation testing difference in difference")
source(here("scripts","zissou-figs.R"))
message("finished making figures")
knitr::knit(here("documenrts","zissou-pnas","zissou-pnas-v3.Rmd"))