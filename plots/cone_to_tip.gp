# EXTREME CASE 3b: the ground state of a narrowing well (V0 = 10 eV), in the style of thesis Figs. 11-12.
# Inside the well psi = cos(alpha x) (rounded top), outside psi ~ exp(-beta|x|) (straight exponential flanks).
# As a -> 0 the rounded part shrinks and the curve turns from a smooth cone into the pointed
# delta-function cusp.  Each curve is scaled to psi(0) = 1 so the shapes can be compared.
# usage:  gnuplot plots/cone_to_tip.gp      (reads data/E_V0_10_a_<a>.dat for the widths below)
load "plots/style.gp"
set output "figures/cone_to_tip.png"

hbar_sq = 7.6199682
V0 = 10.0
widths = "1 0.5 0.25 0.1"
array shade[4] = ["#a9cdf5", "#5d9ee8", "#2a78d6", "#123f75"]   # light = wide well, dark = narrow well

# ground-state energy of each well (row 1, column 3 of its energies file)
array E[4]
do for [i=1:4] {
    stats "data/E_V0_10_a_".word(widths,i).".dat" every ::0::0 using 3 nooutput
    E[i] = STATS_min
}

alpha(E)  = sqrt(2*E/hbar_sq)
beta(E)   = sqrt(2*(V0 - E)/hbar_sq)
psi(x, a, E) = abs(x) <= a ? cos(alpha(E)*x) : cos(alpha(E)*a)*exp(-beta(E)*(abs(x) - a))

set title "Ground state of a narrowing well (V_0 = 10 eV): from rounded cone to delta-function tip"
set xlabel "x (Å)"
set ylabel "ψ_1(x) / ψ_1(0)"
set xrange [-10:10]
set yrange [0:1.05]
set samples 4001
set key top right

set multiplot
plot for [i=1:4] psi(x, real(word(widths,i)), E[i]) with lines lw 2.5 lc rgb shade[i] \
         title sprintf("a = %s Å   (E_1 = %.4f eV)", word(widths,i), E[i])

# inset: zoom on the top, where the rounded (cos) part meets the exponential flanks
set origin 0.08, 0.30
set size 0.30, 0.45
set title "zoom on the top" font ",10"
set xrange [-1.2:1.2]
set yrange [0.6:1.02]
set xlabel ""
set ylabel ""
set xtics 0.5 font ",9"
set ytics 0.1 font ",9"
set object 1 rectangle from graph 0,0 to graph 1,1 behind fc rgb "white" fs solid noborder
unset key
plot for [i=1:4] psi(x, real(word(widths,i)), E[i]) with lines lw 2.5 lc rgb shade[i]
unset multiplot
