# Eduardo Vudala Senoski GRR20195689

export OMP_PROC_BIND=true
export OMP_DYNAMIC=false

END=10
testar () {
    export OMP_NUM_THREADS=$1
    SUM=0
    for i in $(seq 1 $END);
    do
        start=$(date +%s.%N)
        $2
        dur=$(echo "$(date +%s.%N) - $start" | bc)
        SUM=$(echo "$SUM + $dur" | bc)
    done
    echo -n "Average: "
    echo "scale=9; $SUM / $END" | bc
}

testar 8 ./dna