while true
do
	a=$(nproc --all)
	b=$(top -b -n2 -p 1 | fgrep "Cpu(s)" | tail -1 | awk -F'id,' -v prefix="$prefix" '{ split($1, vs, ","); v=vs[length(vs)]; sub("%", "", v); printf "%s%.1f\n", prefix, 100 - v }')
	b=${b%.*}
	c=$((100-$b))
	d=$(($a*$c))
	e=$(($d/100-4))
	if (( $e <= 0 )); then
		e=1
	fi
	echo $e
        mpirun -np $e blah solutions.txt 100 10000 235000000
done
