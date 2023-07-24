# Readme/Notes for DAO and Halo mass function code

## James' email
I'm attaching two abominations of code that seem to work, which you can use to calculate the DAO and diffusion scales, as well as the Sheth-Tormen halo mass function as described in dark recombination. Note the (non-ADM) cosmological parameters are currently set by global variables...Sorry it's a bit of a mess, but I basically ran out of time. You're of course free to ignore, or use as a reference point to redo it.

You could calculate the fraction of the total matter density in a given mass range of halos as (I think):
$f=( \int M \frac{dn}{d\log M} d\log M) /(\Omega_M \rho_c)$

*(Michael's Note: I have modified the original files to try to get them working. So far a failure. The original code can be found in `api_original.jl` and `daohmf_original.jl`)*

## Comments from meeting
(This is word for word from my meeting notes and is likely wrong. For example, I think it was probably Angulo 2013 that James was referencing)
Use Dark Recombination constraints for DAO -see Angulo (sp?) 2023

If find updated WDM constraints, see if they have scale

Look for state of the art WDM constraints - compare with Angulo 2023