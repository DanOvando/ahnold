
p = rnorm(length(parm.names))
Rprof()

b <- proc.time()
for (i in 1:100)
{
  a = mlpa_delta_likelihood(parm = p , Data = Data)
}
proc.time() - b

Rprof(NULL)
RProfData<- readProfileData('Rprof.out')
flatProfile(RProfData,byTotal=TRUE)
