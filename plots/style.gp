# Common look for every figure (loaded by the other scripts)
set terminal pngcairo size 1000,650 enhanced font "Sans,12" linewidth 1
set encoding utf8

# colours: even states = blue, odd states = orange, reference/theory curves = aqua, text and axes = grey
even_c = "#2a78d6"
odd_c  = "#eb6834"
ref_c  = "#1baf7a"
ink    = "#52514e"
light  = "#d8d7d2"

set border 3 lc rgb ink lw 1
set tics nomirror tc rgb ink
set xlabel tc rgb ink
set ylabel tc rgb ink
set title tc rgb "#0b0b0b"
set grid lc rgb "#ecebe7" lw 1
set key tc rgb ink box lc rgb light opaque samplen 2.5 spacing 1.2

set style line 1 lc rgb even_c lw 2
set style line 2 lc rgb odd_c  lw 2
set style line 3 lc rgb ref_c  lw 2
set style line 4 lc rgb "#0b0b0b" lw 2.5
set style line 5 lc rgb ink lw 1 dt 2
