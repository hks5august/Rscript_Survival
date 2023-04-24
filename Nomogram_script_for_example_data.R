library(rms)
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
d <- data.frame(age = rnorm(n, 50, 10),
                blood.pressure = rnorm(n, 120, 15),
                cholesterol = rnorm(n, 200, 25),
                sex = factor(sample(c('female','male'), n,TRUE)))


L <- 0.01 + .1*(d$sex=='male') + -.02*d$age + 
  0.01*d$cholesterol + -0.01*d$blood.pressure
p <- plogis(L)
d$y <- rbinom(n = 1000, size = 1, prob = p)
table(d$y)

f <- lrm(y ~ age + sex + cholesterol + blood.pressure,
         data = d)
f

#Create a nomogram
ddist <- datadist(d)
options(datadist='ddist')

nom <- nomogram(f, fun=plogis, funlabel="Risk of Death")
plot(nom)