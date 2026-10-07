# f_even(E) and f_odd(E) on the same plot; the roots found by the hybrid method are marked.
# usage:  gnuplot -e "V0=10; a=3; tag='V0_10_a_3'" plots/transcendental.gp
load "plots/style.gp"
set output "figures/transcendental_".tag.".png"

tfile = "data/transcendental_".tag.".dat"
efile = "data/E_".tag.".dat"

set title sprintf("Transcendental equations, V_0 = %g eV, a = %g Å", V0, a)
set xlabel "E (eV)"
set ylabel "f(E)  (1/Å)"
set xrange [0:V0]
set key tmargin center horizontal
set xzeroaxis lc rgb ink lw 1

plot tfile using 1:2 with lines ls 1 title "f_{even}(E) = α sin αa − β cos αa", \
     tfile using 1:3 with lines ls 2 title "f_{odd}(E) = α cos αa + β sin αa", \
     efile using 3:(strcol(2) eq "even" ? 0 : NaN) with points pt 7 ps 1.6 lc rgb even_c title "even roots", \
     efile using 3:(strcol(2) eq "odd"  ? 0 : NaN) with points pt 5 ps 1.4 lc rgb odd_c  title "odd roots", \
     efile using 3:(0):(sprintf("%.4f", $3)) with labels offset 0,-1.3 tc rgb ink font ",10" notitle
