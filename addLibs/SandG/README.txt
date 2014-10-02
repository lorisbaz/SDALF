
This is a set of MATLAB m-files implementing the mixture
fitting algorithm described in the paper
M. Figueiredo and A.K.Jain, "Unsupervised learning of
finite mixture models",  IEEE Transaction on Pattern Analysis
and Machine Intelligence, vol. 24, no. 3, pp. 381-396, March 2002.

It consists of a main MATLAB function called "mixtures4.m" and
three auxiliary functions: "uninorm.m", "multinorm.m", and
"elipsnorm.m", which are called by the main program.

For instructions type "help mixtures4" at the MATLAB prompt,
or read the first few lines of the "mixtures4.m" file.

Also included are two simple demos which exemplify how to use the
program. The first one, "demo1.m", uses the three component mixture
from the paper
N. Ueda and R. Nakano, "Deterministic annealing EM algorithm",
Neural Networks, vol. 11, pp. 271-282, 1998.

The second one, "demo2.m", uses the "Simulated Set 2" from the
book by McLachlan and Peel, 2000, page 218. These demos call a
function "genmix.m", which generates samples from a Gaussian
mixture; this function is also included.
