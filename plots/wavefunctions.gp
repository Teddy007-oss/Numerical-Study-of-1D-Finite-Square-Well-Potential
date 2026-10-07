# All bound states (even AND odd) of one well on the same plot.
# Each wavefunction is drawn on top of its own energy level E_n, inside the potential V(x).
# usage:  gnuplot -e "V0=10; a=3; tag='V0_10_a_3'" plots/wavefunctions.gp
load "plots/style.gp"

efile  = "data/E_".tag.".dat"
wffile = "data/wf_".tag.".dat"
set output "figures/wf_".tag.".png"

# read the energies (column 3 of the energies file)
stats efile using 3 nooutput
N = STATS_records
array E[N]
do for [i=1:N] {
    stats efile every ::i-1::i-1 using 3 nooutput
    E[i] = STATS_min
}

# scale every wavefunction to the same height h so the levels do not overlap
array M[N]
do for [i=1:N] {
    stats wffile using (abs(column(i+2))) nooutput
    M[i] = STATS_max
}
h = 0.4*V0/(N+1)

# states alternate even, odd, even, ... starting with the even ground state
col(i) = (i%2 == 1) ? even_c : odd_c
par(i) = (i%2 == 1) ? "even" : "odd"

stats wffile using 1 nooutput
L = STATS_max

set title sprintf("Finite square well, V_0 = %g eV, a = %g Å:  %d bound state%s", V0, a, N, N>1 ? "s" : "")
set xlabel "x (Å)"
set ylabel "Energy (eV)   [ψ_n drawn on its level E_n]"
set xrange [-L:L]
set yrange [*:V0 + 2.2*h]
set key tmargin center horizontal

# label each level on the right (below the line if the level is close to the top of the well)
do for [i=1:N] {
    set label i sprintf("n=%d %s  E=%.7g eV", i, par(i), E[i]) at graph 0.985, first E[i] right \
        offset 0, ((V0 - E[i]) < 0.05*V0 ? -0.8 : 0.6) tc rgb ink font ",10" front
}

plot wffile using 1:2 with lines ls 4 title "V(x)", \
     for [i=1:N] E[i] with lines ls 5 notitle, \
     for [i=1:N] wffile using 1:(E[i] + h*column(i+2)/M[i]) with lines lw 2.5 lc rgb col(i) notitle, \
     NaN with lines ls 1 lw 2.5 title "even states", \
     NaN with lines ls 2 lw 2.5 title (N > 1 ? "odd states" : "")
