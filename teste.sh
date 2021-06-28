# Eduardo Vudala Senoski GRR20195689

END=10
testar () {
    SUM=0
    for i in $(seq 1 $END);
    do
        start=$(date +%s.%N)
        mpirun --hostfile hosts_file.txt -np $1 $2
        dur=$(echo "$(date +%s.%N) - $start" | bc)
        SUM=$(echo "$SUM + $dur" | bc)
    done
    echo -n "Average: "
    echo "scale=9; $SUM / $END" | bc
}

testar 8 ./dna80