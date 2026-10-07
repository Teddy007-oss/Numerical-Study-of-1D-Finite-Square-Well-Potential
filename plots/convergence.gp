# Error |E_k - E| after each iteration k, hybrid vs bisection, for one energy level.
# usage:  gnuplot -e "name='even_6.0_6.4'" plots/convergence.gp
load "plots/style.gp"
file = "data/convergence_".name.".dat"
set output "figures/convergence_".name.".png"

# the converged root = last hybrid iterate
stats file index 0 using 2 nooutput
stats file index 0 every ::STATS_records-1 using 2 nooutput
root = STATS_min

err(E) = abs(E - root) > 1e-16 ? abs(E - root) : 1e-16

set title sprintf("Convergence to E = %.10f eV", root)
set xlabel "iteration k"
set ylabel "|E_k − E|  (eV)"
set logscale y
set format y "10^{%L}"
set yrange [1e-16:1]
set xrange [0:*]
set key top right

plot file index 1 using 1:(err($2)) with linespoints ls 5 lc rgb ink pt 6 ps 1 title "bisection", \
     file index 0 using 1:(err($2)) with linespoints ls 3 lw 2.5 pt 7 ps 1.5 title "hybrid (bisection + quadratic)"
