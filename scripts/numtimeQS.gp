#!/bin/bash
reset 
# DEFINING COMMAND TO READ FROM files
at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
# DONE

# READING SPECTRUM RANGE (omin omax) FROM "../data/input/nanosphere_eV.dat"
file="../data/input/nanosphere_eV.dat" ;
row=1 ; col=5
omin = at(file,row,col)
row=1 ; col=6
omax = at(file,row,col)
# DONE


# READING GAIN QUANTITY (G) FROM "../data/input/nanosphere_eV.dat"
row=1 ; col=4
G = abs(at(file,row,col))
# DONE

# READING GAIN THRESHOLD (Gth) FROM "../data/output/frohlich.dat"
file="../data/output/frohlich.dat" ;
row=1 ; col=2
Gth = abs(at(file,row,col))
# DONE

# READING FREQUENCY (ome) FROM "../data/output/omega.dat"
file="../data/output/omega.dat" ;
row=1 ; col=1
omega = abs(at(file,row,col))
# DONE

# READING TOTAL TIME AND TIME_PUMP_ON (ome T tpump) FROM "../data/input/time.dat"
file="../data/input/time.dat" ;
row=1 ; col=1
T = abs(at(file,row,col))
row=1 ; col=2
tpump = at(file,row,col)
# DONE

# READING STATIONARY QUASI STATIC ALPHA (Re(alph_QS), Im(alph_QS)) FROM "../data/output/alpha.dat"
file="../data/output/alpha.dat" ;
row=1 ; col=1
reaQS = at(file,row,col)
row=1 ; col=2
imaQS = at(file,row,col)

# DONE
set print "-"
stats "../data/output/stationary.dat" using 2 name "Y2" nooutput
stats "../data/output/stationary.dat" using 3 name "Y3" nooutput
set print

# Calcola il valore massimo tra le due colonne
ymax = (Y2_max > Y3_max) ? Y2_max : Y3_max
ymin = (Y2_min < Y3_min) ? Y2_min : Y3_min

# Imposta il terminale e la gamma y
set term unknown
if (ymax > 200) {
    set yrange [-200:200]
    set y2range [-200:200]
    } else {
    set yrange [1.1*ymin:1.1*ymax]
    set y2range [1.1*ymin:1.1*ymax]
    }

set xrange [0:T];
set term pdf color enhanced size 20cm, 8cm font ",16";
set output '../img/intime.pdf';
set style line 6 lt 1 lc 0 lw 4;
set style line 7 lt 1 lc 8 lw 4;
set size 2,1.;
set title sprintf("G = %.2f",G);
set multiplot
set size 0.56,1.;

set ylabel "Normalized Polarizability {/Symbol a}/(4{/Symbol p}a^3)"

# set grid ytics
set label 22 at graph 0.02, 0.08 "{/=24 (a)}";

if (G!=0) { 
    set label 33 at first 0.22, graph 0.95 "{/=14 PUMP ON}";
        set object 66 rect from first tpump, graph 0 to first T, graph 1  fc rgb "light-blue" fillstyle transparent solid;
    }

set xlabel "Time (ps)"

plot "../data/output/numtime.dat" u 1:2 w l ls 6 t "",\
     "../data/output/numtime.dat" u 1:3 w l ls 7 t "";
     
dome=(omax-omin)/8.
set xtics omin+dome, dome
set xrange [omin:omax];
set label 22 at graph 0.02, 0.08 "{/=24 (b)}";
unset label 33
set title sprintf("ℏ{/Symbol w} = %.2f",omega) tc rgb "blue";
unset ylabel

unset ytics
set y2label "Normalized Polarizability {/Symbol a}/(4{/Symbol p}a^3)"
set y2tics mirror
set xlabel "Energy (eV)"

set origin 0.52,0
set rmargin 17

set object 66 rect from graph 0, graph 0 to graph 1, graph 1
set arrow 22 nohead from first omega, graph 0 to first omega, graph 1 lc rgb "blue" lw 2 front
set arrow 23 nohead from graph 0, first reaQS to first omega, first reaQS lc 0 lw 2 dt 2
set arrow 24 nohead from graph 0, first imaQS to first omega, first imaQS ls 7 lw 2 dt 2

plot "../data/output/stationary.dat" u 1:2 axes x1y2 w l ls 6 t "{/Symbol a}' ",\
     "../data/output/stationary.dat" u 1:3 axes x1y2 w l ls 7 t "{/Symbol a}''"; 

unset multiplot
unset term
! convert -density 200 ../img/intime.pdf ../img/intime.png
