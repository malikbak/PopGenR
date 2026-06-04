# variance.d: Calculate Variance of the Genetic Diversity Index (D)

This function calculates the variance of the genetic diversity index (D)
using the number of samples and the observed diversity.

## Usage

``` r
variance.d(n, S)
```

## Arguments

- n:

  An integer representing the number of samples.

- S:

  A numeric value representing the observed genetic diversity.

## Value

A numeric value representing the variance of the genetic diversity index
(D).

## Details

The function uses the following formulas to compute the variance:

\$\$a1 = \sum\_{i=1}^{n-1} \frac{1}{i}\$\$ \$\$a2 = \sum\_{i=1}^{n-1}
\frac{1}{i^2}\$\$

Then it calculates coefficients \\b1\\, \\b2\\, \\c1\\, and \\c2\\ based
on these sums, and finally computes the variance using:

\$\$var = e1 \cdot S + e2 \cdot S \cdot (S - 1)\$\$

where \\e1\\ and \\e2\\ are derived from \\c1\\ and \\c2\\.

## Note

This function assumes that the input values are valid and that the
calculations adhere to the statistical properties of genetic diversity.

## Examples

``` r
# Example usage:
n_samples <- 10
observed_diversity <- 0.5
diversity_variance <- variance.d(n_samples, observed_diversity)

# Output: A numeric value representing the variance of the genetic diversity index
```
