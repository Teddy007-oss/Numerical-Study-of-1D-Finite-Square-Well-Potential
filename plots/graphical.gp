# The textbook graphical solution (Griffiths eq. 2.156 / thesis Fig. 2), with even AND odd states,
# and the hybrid-method roots marked on it.  In terms of z = alpha a  and  z0 = a sqrt(2 m V0)/hbar:
#   even:  tan z = sqrt((z0/z)^2 - 1)        odd:  -cot z = sqrt((z0/z)^2 - 1)
# usage:  gnuplot -e "V0=10; a=3; tag='V0_10_a_3'" plots/graphical.gp
load "plots/style.gp"
set output "figures/graphical_".tag.".png"

efile = "data/E_".tag.".dat"
hbar_sq = 7.6199682
z0 = a*sqrt(2*V0/hbar_sq)
zof(E) = a*sqrt(2*E/hbar_sq)
rhs(z) = (z < z0) ? sqrt((z0/z)**2 - 1) : NaN

set title sprintf("Graphical solution, V_0 = %g eV, a = %g Å  (z_0 = %.4f)", V0, a, z0)
set xlabel "z = αa"
set ylabel ""
set xrange [0:z0*1.08]
set yrange [0:8]
set samples 4000
set key tmargin center horizontal
set arrow 1 from z0, graph 0 to z0, graph 1 nohead ls 5
set label 1 "z_0" at z0, graph 0.96 offset 0.5,0 tc rgb ink

plot (tan(x) > 0 ? tan(x) : NaN) with lines ls 1 title "tan z  (even)", \
     (-1/tan(x) > 0 ? -1/tan(x) : NaN) with lines ls 2 title "−cot z  (odd)", \
     rhs(x) with lines ls 4 title "√((z_0/z)^2 − 1)", \
     efile using (zof($3)):(strcol(2) eq "even" ? rhs(zof($3)) : NaN) with points pt 7 ps 1.6 lc rgb even_c title "hybrid method, even", \
     efile using (zof($3)):(strcol(2) eq "odd"  ? rhs(zof($3)) : NaN) with points pt 5 ps 1.4 lc rgb odd_c  title "hybrid method, odd"
