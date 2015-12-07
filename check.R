
p = rnorm(length(parm.names))
Rprof()

for (i in 1:100)
{
  a = mlpa_delta_likelihood(parm = p , Data = Data,reg_model = 'tobit')
}

Rprof(NULL)
RProfData<- readProfileData('Rprof.out')
flatProfile(RProfData,byTotal=TRUE)
