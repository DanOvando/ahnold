load('FUCKIT.Rdata')

p =Initial.Values
Rprof()


p <- Initial.Values

a <- proc.time()

for (i in 1:100)
{
  a = mlpa_delta_likelihood(parm = p , Data = Data)
}
proc.time() - b

Rprof(NULL)
RProfData<- readProfileData('Rprof.out')
flatProfile(RProfData,byTotal=TRUE)
