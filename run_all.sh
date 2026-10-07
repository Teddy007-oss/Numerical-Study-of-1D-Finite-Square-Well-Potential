#!/bin/bash
# Builds every program, runs every case of the study and draws every figure.
# usage:  ./run_all.sh            (needs gfortran and gnuplot)
#         GNUPLOT=/path/to/gnuplot ./run_all.sh
set -e
export LC_NUMERIC=C
cd "$(dirname "$0")"
GNUPLOT=${GNUPLOT:-gnuplot}
mkdir -p bin data figures

echo "=== 1. Root finders tested on simple equations ==="
for p in bisection quadratic_approximation hybrid; do
    gfortran root_finding_tests/$p.f90 -o bin/test_$p
done
echo "bisection,  f(x) = x^2 - 1 on [0, 3]:";          printf "0\n3\n" | bin/test_bisection | tail -1
echo "quadratic approximation,  f(x) = 2x - 6:";        bin/test_quadratic_approximation | tail -1
echo "hybrid,  f(x) = 5x^2 + 12x + 7 on [-1.2, -2.8]:"; bin/test_hybrid | tail -1

echo "=== 2. Building the finite square well programs ==="
for p in transcendental energies wavefunctions convergence; do
    gfortran -O2 -J bin finite_well/roots.f90 finite_well/well.f90 finite_well/$p.f90 -o bin/$p
done

# energies + wavefunctions (+ the all-states plot) for one well:  one_well V0 a
one_well () {
    tag="V0_$1_a_$2"
    echo "$1 $2" | bin/energies      > data/E_$tag.dat
    echo "$1 $2" | bin/wavefunctions > data/wf_$tag.dat
    $GNUPLOT -e "V0=$1; a=$2; tag='$tag'" plots/wavefunctions.gp
    grep "bound states" data/E_$tag.dat | sed "s/^#/  V0=$1 a=$2:/"
}

echo "=== 3. The reference well V0 = 10 eV, a = 3 A ==="
one_well 10 3
echo "10 3" | bin/transcendental > data/transcendental_V0_10_a_3.dat
$GNUPLOT -e "V0=10; a=3; tag='V0_10_a_3'" plots/transcendental.gp
$GNUPLOT -e "V0=10; a=3; tag='V0_10_a_3'" plots/graphical.gp
cat data/E_V0_10_a_3.dat

echo "=== 4. Convergence: hybrid vs bisection (thesis brackets) ==="
echo "10 3 even 6.0 6.4" | bin/convergence > data/convergence_even_6.0_6.4.dat
echo "10 3 odd 9.5 10.0" | bin/convergence > data/convergence_odd_9.5_10.0.dat
tail -1 data/convergence_even_6.0_6.4.dat
tail -1 data/convergence_odd_9.5_10.0.dat
$GNUPLOT -e "name='even_6.0_6.4'"  plots/convergence.gp
$GNUPLOT -e "name='odd_9.5_10.0'" plots/convergence.gp

echo "=== 5. Varying the depth V0 (a = 3 A) and the width a (V0 = 10 eV) ==="
for V0 in 6 2 1 0.5; do one_well $V0 3;   done
for a  in 2 1 0.5;    do one_well 10 $a;  done
rm -f data/sweep_V0.dat data/sweep_a.dat
for V0 in $(seq 0.1 0.1 20); do echo "$V0 3"  | bin/energies | grep -v "^#" >> data/sweep_V0.dat; done
for a  in $(seq 0.05 0.05 6); do echo "10 $a" | bin/energies | grep -v "^#" >> data/sweep_a.dat;  done
$GNUPLOT -e "sweep='V0'" plots/spectrum.gp
$GNUPLOT -e "sweep='a'"  plots/spectrum.gp

echo "=== 6. Extreme cases ==="
# deep well -> infinite square well
rm -f data/deep.dat
for V0 in 10 100 1000 10000; do echo "$V0 3" | bin/energies | grep -v "^#" >> data/deep.dat; done
$GNUPLOT plots/deep_vs_infinite.gp
# shallow well and narrow well -> delta-function well (z0 from about 5 down to 0.01)
rm -f data/shallow.dat data/narrow.dat
for k in $(seq 0 60); do
    s=$(awk -v k=$k 'BEGIN{printf "%.6e", 10^(-2.7*k/60)}')        # 1 down to 0.002
    echo "$(awk -v s=$s 'BEGIN{print 10*s*s}') 3" | bin/energies | awk '$1 == 1' >> data/shallow.dat   # ground state only
    echo "10 $(awk -v s=$s 'BEGIN{print 3*s}')"   | bin/energies | awk '$1 == 1' >> data/narrow.dat
done
$GNUPLOT plots/weak_binding.gp
# fixed strength g = 2 a V0 = 10 eV A, a -> 0
for w in "10 0.5" "100 0.05" "1000 0.005"; do
    set -- $w
    echo "$1 $2" | bin/wavefunctions > data/wf_V0_$1_a_$2.dat
done
$GNUPLOT plots/delta_limit.gp
# wavefunctions of the extreme wells: deep, very shallow, very narrow
one_well 50 3
one_well 0.01 3
one_well 10 0.1
# narrowing well at fixed V0 = 10 eV: rounded cone -> pointed delta-function tip (thesis Figs. 11-12 style)
echo "10 0.25" | bin/energies > data/E_V0_10_a_0.25.dat
$GNUPLOT plots/cone_to_tip.gp

echo "=== done: figures are in figures/, numbers in data/ ==="
