from numpy import linspace, exp
from matplotlib.pyplot import plot, show, xlabel, ylabel

alpha = 2
r = linspace(-2, 5)
func = exp(-alpha*r)
plot(r, func)
xlabel("r")
ylabel(r"$\exp(-\alpha r)$")
show()
