

nn=$(grep -c "S" $ms2file)
charge=$2
index=$1


if [ "$nn" = "0" ]
then
	echo "not enough spectra"
	exit
fi
n=$(grep -c "." $xcorr)

if [ "$n" = "1" ]
then
	echo "not enough spectra"
	exit
fi


./OUT_SGM_sup.out -i ../data/encode/$1-ch$2/msms-fastSequestTransform-$1-ch$2.txt -m ../data/ms2file/$1-ch$2.ms2 -o ../data/result/$1-ch$2.txt
