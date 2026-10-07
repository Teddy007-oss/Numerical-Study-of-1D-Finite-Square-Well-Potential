# EXTREME CASE 3: narrow AND deep well with fixed strength g = 2 a V0 = 10 eV Å.
# As a -> 0 the ground state wavefunction must turn into the delta-well wavefunction
#   psi(x) = sqrt(kappa) exp(-kappa |x|),   kappa = m g / hbar^2
# usage:  gnuplot plots/delta_limit.gp
load "plots/style.gp"
set output "figures/delta_limit.png"

hbar_sq = 7.6199682
g = 10.0
kappa = g/hbar_sq
psi_delta(x) = sqrt(kappa)*exp(-kappa*abs(x))

tags = "V0_10_a_0.5 V0_100_a_0.05 V0_1000_a_0.005"
array names[3] = ["V_0 = 10 eV,  a = 0.5 Å", "V_0 = 100 eV,  a = 0.05 Å", "V_0 = 1000 eV,  a = 0.005 Å"]
array shade[3] = ["#a9cdf5", "#5d9ee8", "#174a8c"]

set title sprintf("Delta-function limit at fixed strength g = 2aV_0 = %g eV·Å", g)
set xlabel "x (Å)"
set ylabel "ψ_1(x)  (1/√Å)"
set xrange [-3:3]
set yrange [0:*]
set samples 2001
set key top right

plot for [i=1:3] "data/wf_".word(tags,i).".dat" using 1:3 with lines lw 2.5 lc rgb shade[i] \
         title "finite well  ".names[i], \
     psi_delta(x) with lines ls 3 dt 2 lw 2.5 title "delta well  √κ e^{−κ|x|}"
