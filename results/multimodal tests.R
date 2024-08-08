# Hartigan's dip test
# J. A. Hartigan and P. M. Hartigan (1985) The Dip Test of Unimodality; Annals of Statistics 13, 70–84.
library (diptest)

x1 <- c(rnorm(40,0),
        rnorm (40, 4))
dip.test(x1)


x2 <- c(rnorm(20,0),
        rnorm (20, 0))
dip.test(x2)

library (flexmix)
data("NPreg")
m0 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 1)
m1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)
parameters(m1, component = 1)
summary (m0)
summary (m1)

m2 <- flexmix (x1 ~ 1, k =2)
parameters(m2, component = 2)



