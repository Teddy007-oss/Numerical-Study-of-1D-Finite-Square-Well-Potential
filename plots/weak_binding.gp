# EXTREME CASE 2: very shallow OR very narrow well (z0 -> 0).
# Only the even ground state survives, and it should approach the bound state of a delta-function well
# of the same strength g = 2 a V0, whose binding energy is  E_delta = m g^2 / (2 hbar^2).
# The finite-well binding energy is E_B = V0 - E_1.  The ratio E_B / E_delta must go to 1.
# Making the well shallower (a fixed) and making it narrower (V0 fixed) must give the SAME curve,
# because E_B / E_delta depends only on z0.
# usage:  gnuplot plots/weak_binding.gp      (reads data/shallow.dat and data/narrow.dat)
load "plots/style.gp"
set output "figures/weak_binding.png"

hbar_sq = 7.6199682
z0(V0, a)     = a*sqrt(2*V0/hbar_sq)
Edelta(V0, a) = (2*a*V0)**2/(2*hbar_sq)

set title "Weak-binding limit: ground state vs delta-function well"
set xlabel "well strength z_0 = a√(2mV_0)/ħ"
set ylabel "E_B / E_δ   (binding energy / delta-well binding energy)"
set logscale x
set xrange [0.01:10]
set yrange [0:1.05]
set key bottom left

plot 1 with lines ls 3 dt 2 title "delta-function well", \
     "data/shallow.dat" using (z0($6,$7)):($1 == 1 ? ($6 - $3)/Edelta($6,$7) : NaN) with lines ls 1 lw 3 title "shallower well (a = 3 Å, V_0 decreasing)", \
     "data/narrow.dat"  using (z0($6,$7)):($1 == 1 ? ($6 - $3)/Edelta($6,$7) : NaN) with points pt 6 ps 1.1 lc rgb odd_c title "narrower well (V_0 = 10 eV, a decreasing)"
