# Every bound-state energy as a function of one well parameter (V0 or a), even and odd together.
# A new state appears every time z0 = a sqrt(2 m V0)/hbar passes a multiple of pi/2.
# usage:  gnuplot -e "sweep='V0'" plots/spectrum.gp     (data/sweep_V0.dat, a fixed)
#         gnuplot -e "sweep='a'"  plots/spectrum.gp     (data/sweep_a.dat,  V0 fixed)
load "plots/style.gp"
file = "data/sweep_".sweep.".dat"
set output "figures/spectrum_vs_".sweep.".png"

# columns of the sweep file: n parity E it_hybrid it_bisection V0 a
xcol = (sweep eq "V0") ? 6 : 7
stats file using 7 nooutput
afix = STATS_min
stats file using 6 nooutput
Vfix = STATS_min

if (sweep eq "V0") {
    set title sprintf("Bound-state energies vs well depth  (a = %g Å)", afix)
    set xlabel "well depth V_0 (eV)"
} else {
    set title sprintf("Bound-state energies vs well half-width  (V_0 = %g eV)", Vfix)
    set xlabel "well half-width a (Å)"
}
set ylabel "E_n (eV)"
set key tmargin center horizontal
set xrange [0:*]
set yrange [0:*]

top(x) = (sweep eq "V0") ? x : Vfix
plot top(x) with lines ls 5 title "E = V_0 (top of the well)", \
     file using xcol:(strcol(2) eq "even" ? $3 : NaN) with points pt 7 ps 0.5 lc rgb even_c title "even states", \
     file using xcol:(strcol(2) eq "odd"  ? $3 : NaN) with points pt 7 ps 0.5 lc rgb odd_c  title "odd states"
