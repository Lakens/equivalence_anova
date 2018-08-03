Fstat=3.25
df1=4
df2=50
alpha = 0.05
tol = 1e-09

#Calculate lower ncp
ll <- optimize(function(q) (pf(Fstat, 
                               df1, 
                               df2, 
                               ncp = q) - (1-alpha/2))^2,
               interval = c(0,
                            Fstat * df1),
               tol = tol)$minimum
ll
#Calculate upper ncp
ul <- optimize(function(q) (pf(Fstat, 
                               df1, 
                               df2, 
                               ncp = q) - alpha/2)^2,
               interval = c(Fstat * df1,
                            Fstat * df1 * 10),
               tol = tol)$minimum
ul

#This does the same as the conf.limits.ncf function in MBESS
library(MBESS)
ul_ll <- conf.limits.ncf(F.value = Fstat, conf.level = .95, df.1 = df1, df.2 = df2)

#The values are acceptably close
ul_ll[[1]] - ll
ul_ll[[3]] - ul

#From example Harlan
Fstat=46.99195
df1=1
df2=147
N=150
conf.level=0.90

ci.pvaf(F.value=Fstat, df.1=df1, df.2=df2, N=N, conf.level=0.90)

#Now not using MBESS
alpha = 0.1
ll <- optimize(function(q) (pf(Fstat, 
                               df1, 
                               df2, 
                               ncp = q) - (1-alpha/2))^2,
               interval = c(0,
                            Fstat * df1),
               tol = tol)$minimum
ll

#Then, the only thing the ci.pvaf function does is:
#Lower.Limit.Proportion.of.Variance.Accounted.for
ll/(ll+N)

#Calculate upper ncp
ul <- optimize(function(q) (pf(Fstat, 
                               df1, 
                               df2, 
                               ncp = q) - alpha/2)^2,
               interval = c(Fstat * df1,
                            Fstat * df1 * 10),
               tol = tol)$minimum
ul
#Upper.Limit.Proportion.of.Variance.Accounted.for
ul/(ul+N)
