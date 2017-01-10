# QVE code

This contains some code that implements several algorithms to solve quadratic vector equations from my papers.

A quadratic vector equation is an equation of the form
$$
Mx = a + b(x \otimes x),
$$
where $a\in\mathbb{R}^n$, $M\in\mathbb{R}^{n\times n}$, $b\in\mathbb{R}^{n\times n^2}$.

The main case we are interested in is the one in which $a, b\geq 0$ (elementwise) and $M$ is an M-matrix. In this case, one can define a minimal solution, and most algorithms here converge to it. For more detail, see [QVE].

Algorithms `depth`, `order`, `thicknesses` and `qve_newton` are from [EP]. These are guaranteed to converge to the minimal nonnegative solution, under the hypotheses in [QVE].

Algorithm `perron_iteration` is from [P], and algorithm `perron_newton` is from [PN]. These two algorithms assume that $M=I$ (which is not restrictive: just set $a\leftarrow M^{-1}a$, $b\leftarrow M^{-1}b$), and that the equation has a solution $e=[1,1,1,\dots,1]^T$ that we are not interested with, and a second real solution that we are interested in computing. It is possible to adapt them to "simplify" other known solutions than $e$. (In lack of better ideas, you can always rescale $a,b$ so that $e$ is a solution.)

## Bibliography

[QVE]: [Quadratic vector equations](http://dx.doi.org/10.1016/j.laa.2011.05.036) (Federico Poloni), In Linear Algebra and its Applications, volume 438, 2013.

[P]: [A Perron iteration for the solution of a quadratic vector equation arising in Markovian binary trees](http://dx.doi.org/10.1137/100796765) (Beatrice Meini, Federico Poloni), In SIAM J. Matrix Anal. Appl., volume 32, 2011.

[PN]: [On the solution of a quadratic vector equation arising in Markovian binary trees](http://dx.doi.org/10.1002/nla.809) (Dario A. Bini, Beatrice Meini, Federico Poloni), In Numer. Linear Algebra Appl., volume 18, 2011.

[EP]: S. Hautphenne, G. Latouche, and M.-A. Remiche. [Algorithmic approach to the extinction probability of branching processes.](http://dx.doi.org/10.1007/s11009-009-9141-7) Methodology and Computing in Applied Probability, 2011, 13(1):171-192.

PDFs of my papers can be obtained free of charge at http://pages.di.unipi.it/fpoloni/publications/publications.php.

## License

Feel free to use this code in your research. Just remember cite me if it is appropriate. I mean in your paper, not "summon me in court". Which reminds me: there is no warranty for this program, it's provided "as is", we are not liable for damage, and all that kind of stuff you read in licenses.
