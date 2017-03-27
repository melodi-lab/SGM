
ms2file="../python/data/plasm-1-ch3.ms2"
xcorr="../python/plasm-1-ch3/xcorrIdent-plasm-1-ch3.txt"

nn=$(grep -c "S" $ms2file)
if [ "$nn" = "0" ]
then
	echo asdf
	exit
fi
n=$(grep -c "." $xcorr)

if [ "$n" = "1" ]
then
	echo asdf
	exit
fi


./OUT_SGM_sup.out -i ../python/plasm-1-ch3/msms-fastSequestTransform-plasm-1-ch3.txt -m ../python/data/plasm-1-ch3.ms2 -o test.txt -n $n -N $nn
