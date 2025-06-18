#!/bin/bash

infile=../data/input/nanosphere_eV.dat
echo I N T I M E
echo
echo "> Compiling sources..."
g++ -Wall -I/usr/local/include -L/usr/local/lib ../src/frohlich.cxx -o ../bin/fro -lgsl -lgslcblas -lm
g++ -Wall -I/usr/local/include -L/usr/local/lib ../src/single.cxx -o ../bin/sgl -lgsl -lgslcblas -lm -larmadillo
g++ -Wall -I/usr/local/include -L/usr/local/lib ../src/eV2cfc.cxx -o ../bin/eV2cfc -lgsl -lgslcblas -lm -lcomplex_bessel -larmadillo
echo "> ..done!"
read r1 Dome ome_0 G omemi omema mtl mdl active sol < $infile
echo ">"
info=`./../bin/fro`
om0=`echo $info | cut -d ' ' -f1`
Gth=`echo $info | cut -d ' ' -f2`
echo $om0 $Gth > ../data/output/frohlich.dat

echo "> Using parameters read form '$infile'" 
echo "> if we center the gain emission in ome0 = $om0"
echo "> The gain needed for the singular behavior is Gth = $Gth"

echo ">"
echo "> Positioning the gain emission in ome0 = $om0"
echo $r1 $Dome $om0 $G $omemi $omema $mtl $mdl $active $sol > $infile
echo "> Done!"
echo ">"



file=("nopump" "hlfGth" "dwnGth" "fllGth" "ovrGth")
# file=("hlfGth" "ovrGth")
mG=("0." "0.5" "0.9" "1." "1.1")
# mG=("0." "1.1")
for j in {0..4..1} 
# for j in {0..1..1} 
    do 
        echo echo "> Plasmon with G = ${mG[$j]}*Gth"
        wG=`echo ${mG[$j]}*$Gth|bc -l`
        echo $r1 $Dome $om0 $wG $omemi $omema $mtl $mdl $active $sol > $infile
        rm -fr "../img/${file[$j]}"
        mkdir  "../img/${file[$j]}"
        mkdir  "../img/${file[$j]}/pdf"
        mkdir  "../img/${file[$j]}/png"
        
        N=100;
        dig="${#N}";
        frm="%0${dig}d "
        dom=`echo "($omema-$omemi)/$N"|bc -l`
        for (( i=0; i<=$N; i++ ))
        do 
            ome=`echo "$omemi+$i*$dom"|bc -l`
            ./../bin/sgl $ome
            ./../bin/eV2cfc $ome
            echo $ome > ../data/output/omega.dat
            gnuplot timeQS.gp
            name=`printf $frm $i`
            echo $name
            mv ../img/intime.pdf ../img/${file[$j]}/pdf/"im-$name.pdf"
            mv ../img/intime.png ../img/${file[$j]}/png/"im-$name.png"
        done
        convert -dispose  previous -delay 20 ../img/${file[$j]}/png/*.png -loop 0 ../img/${file[$j]}.gif
    done
exit

rm -fr "../img/nopump"
mkdir  "../img/nopump"
mkdir  "../img/nopump/pdf"
mkdir  "../img/nopump/png"

N=100;
dig="${#N}";
frm="%0${dig}d "

dom=`echo "($omema-$omemi)/$N"|bc -l`
for (( i=0; i<=$N; i++ ))
  do 
     ome=`echo "$omemi+$i*$dom"|bc -l`
     ./../bin/anl $ome
     gnuplot time.gnp
     name=`printf $frm $i`
     echo $name
     mv ../img/intime.pdf ../img/nopump/pdf/"$name.pdf"
     mv ../img/intime.png ../img/nopump/png/"$name.png"
 done
 convert -dispose  previous -delay 20 ../img/nopump/png/*.png -loop 0 ../img/nopump.gif
echo $r1 $Dome $ome_0 $G $omemi $omema $mtl $mdl $active $sol > $infile
