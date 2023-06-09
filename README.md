# Fourier Transform code for 2D radial functions

Uses Scipy's `fft` package to do fourier transforms of 2D (real) radial functions, according to the definition

$$ f(\mathbf{b}) = \int \frac{d^2 \mathbf{Q}}{(2\pi)^2} e^{-i\mathbf{Q}\cdot\mathbf{b}}F(|\mathbf{Q}|) $$

The function $F(|\mathbf{Q}|)$ is defined inside `func(Q)`, while the number of points is defined as `N`, and the domain over which $F(|\mathbf{Q}|)$ is defined $|\mathbf{Q}| \in [0, Q_\text{max}]$ is defined as `Qrange` $= Q_\text{max}$.

The fourier transformed function $f(\mathbf{b})$ is plotted, and data points are exported as `output.csv`. The $b$ values will be $b = \left(0, \frac{2\pi}{Q_\text{max}}, \frac{4\pi}{Q_\text{max}}, ..., \frac{(N-1)(2\pi)}{Q_\text{max}}\right)$.

In this particular case, $F(|\mathbf{Q}|)$ corresponds to some distribution in 2D momentum space $[\text{GeV}]$ while $f(\mathbf{b})$ corresponds to some distribution in 2D position space $[\text{fm}]$. The function in `func(Q)` can be changed to suit your specific case.

More detailed notes on how the code works are in `FT_notes.ipynb`.