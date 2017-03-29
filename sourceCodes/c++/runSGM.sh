
index=$1
charge=$2
msmsfile=../data/encode/$1-ch$2/msms-fastSequestTransform-$1-ch$2.txt 
ms2file=../data/ms2file/$1-ch$2.ms2
xcorr=../data/encode/$1-ch$2/xcorrIdent-$1-ch$2.txt

nn=$(grep -c "S" $ms2file)


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
n=$((($n-1)/2))
echo $nn
echo $n

./OUT_SGM_sup.out -i $msmsfile -m $ms2file -o ../data/result/$1-ch$2.txt -c $2 -n $n -N $nn
