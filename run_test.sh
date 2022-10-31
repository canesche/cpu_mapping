set -e

GRAPH=(
	atax  
    cholesky  
    gemm    
    gesummv  
    symm  
    trmm
    bicg  
    doitgen   
    gemver  
    mvt      
    syrk
)

ARCH=( # 0 = mesh, 1 = 1-hop, 2 = chess, 3 = hex
    0 
    #1
)

SIZE=(
    1 
    10 
    100 
    1000
)

PROG=( # 0 = yolt, 1 = yott, 2 = sa
    0
    1
    2
)

# build and constructor the code
mkdir -p build && cd build && cmake .. && make && cd ..

NUM_THREADS=16

export OMP_NUM_THREADS=$NUM_THREADS

mkdir -p results

for ((l=0; l < ${#PROG[@]}; l++)) do
    echo "Tool: "${PROG[l]}
    for ((k=0; k < ${#SIZE[@]}; k++)) do
        echo "SIZE: "${SIZE[k]}
        for ((j=0; j < ${#ARCH[@]}; j++)) do
            
            if [ ${ARCH[j]} == 0 ]; then
                mkdir -p results/yolt/mesh/${SIZE[k]}
                mkdir -p results/yott/mesh/${SIZE[k]}
                mkdir -p results/sa/mesh/${SIZE[k]}
                echo "ARCH: MESH"
            elif [ ${ARCH[j]} == 1 ]; then
                mkdir -p results/yolt/1hop/${SIZE[k]}
                mkdir -p results/yott/1hop/${SIZE[k]}
                mkdir -p results/sa/1hop/${SIZE[k]}
                echo "ARCH: 1HOP"
            fi
            
            for ((i=0; i < ${#GRAPH[@]}; i++)) do

                echo "GRAPH: ${GRAPH[i]}"
                DOT="bench/lisa/dac/"${GRAPH[i]}".dot"
        
                # bench ngrids program arch
                ./build/main $DOT ${SIZE[k]} ${PROG[l]} ${ARCH[j]} 
            done
        done
    done
done
