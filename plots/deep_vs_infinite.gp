# EXTREME CASE 1: a very deep well should behave like the infinite square well of width 2a,
#   E_n(infinite) = n^2 pi^2 hbar^2 / (2 m (2a)^2)
# The ratio E_n / E_n(infinite) is plotted for wells of increasing depth; it must approach 1.
# usage:  gnuplot plots/deep_vs_infinite.gp      (reads data/deep.dat)
load "plots/style.gp"
file = "data/deep.dat"
set output "figures/deep_vs_infinite.png"

hbar_sq = 7.6199682
Einf(n, a) = n**2*pi**2*hbar_sq/(8*a**2)

depths = "10 100 1000 10000"
array shade[4] = ["#a9cdf5", "#5d9ee8", "#2a78d6", "#123f75"]   # one hue, light = shallow, dark = deep

stats file using 7 nooutput
afix = STATS_min

set title sprintf("Deep-well limit: finite well vs infinite well  (a = %g Å)", afix)
set xlabel "state number n"
set ylabel "E_n / E_n^{(∞)}"
set logscale x
set xrange [0.8:150]
set yrange [0:1.05]
set key bottom left

plot 1 with lines ls 3 dt 2 title "infinite square well", \
     for [i=1:4] file using 1:($6 == real(word(depths,i)) ? $3/Einf($1, $7) : NaN) \
         with linespoints lw 2 pt 7 ps 0.7 lc rgb shade[i] title sprintf("V_0 = %s eV", word(depths,i))
