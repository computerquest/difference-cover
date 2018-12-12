while true
do
	a=$(nproc --all)
	b=$(grep 'cpu ' /proc/stat | awk '{usage=($2+$4)/($2+$4+$5)*100} END {print usage}')
	b=${b%.*}
	c=$((100-$b))
	d=$(($a*$c))
	e=$(($d/100-4))
	if (( $e <= 0 )); then
		e=1
	fi
	echo $e
        mpirun -np $e blah solutions.txt 20 1000 1000000
done
