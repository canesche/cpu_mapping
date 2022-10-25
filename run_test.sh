set -e

GRAPH=(
	atax
)

# 0: mesh
# 1: 1-hop
# 2: chess
# 3: hex

ARCH=(
0
)

# build and constructor the code
mkdir -p build && cd build && cmake .. && make && cd ..

NUM_THREADS=8

export OMP_NUM_THREADS=$NUM_THREADS

mkdir -p results

SIZE=3

for ((i=0; i < ${#GRAPH[@]}; i++)) do

    echo "GRAPH: ${GRAPH[i]}"
    DOT="bench/lisa/dac/"${GRAPH[i]}".dot"
    
    # bench ngrids program arch
    ./build/main $DOT 1000 0 0 

done
