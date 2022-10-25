# CPU Mapping

CPU Mapping - This project presents an mapping (Placement and routing) for dataflow to three tools: YOLT, YOTT, and SA.

**YOLT** means **Y**ou **O**nly **L**ook **T**wice

**YOTT** means **Y**ou **O**nly **T**raversal **T**wice

**SA** means **S**imulated **A**nnealing

# How to compiler

    mkdir -p build
    cd build
    cmake ..
    make -j4

# How to execute a test

You just need execute the script.

    ./run_test.sh

# How to execute the tools

Directory build after you compile you'll have an executable main.

    ./build/main <DATAFLOW> <NTIMES> <TOOL <ARCH>

**DATAFLOW** is a dot, you can look example in bench directory

**NTIMES** is a number of times your project will try to mapping a code, how bigger is number better it will the solution, but more time the solution will take to finish.

**TOOL** There are three tools available in this project: 

    0: YOLT
    1: YOTT
    2: SA

**ARCH** There are two architectures available in this project:

    0: MESH
    1: 1-HOP

## Paper to cite

YOLT tool:

    @article{canesche2020traversal,
        title={Traversal: A fast and adaptive graph-based placement and routing for cgras}, 
        author={Canesche, Michael and Menezes, Marcelo and Carvalho, Westerley and Torres, Frank Sill and Jamieson, Peter and Nacif, Jos{\'e} Augusto and Ferreira, Ricardo}, 
        journal={IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems},
        volume={40},
        number={8},
        pages={1600--1612},
        year={2020},
        publisher={IEEE}
    }

YOTT tool:

    @article{canesche2021you,
        title={You Only Traverse Twice: A YOTT Placement, Routing, and Timing Approach for CGRAs},
        author={Canesche, Michael and Carvalho, Westerley and Reis, Lucas and Oliveira, Matheus and Magalhaes, Salles and Jamieson, Peter and Nacif, Jaugusto M and Ferreira, Ricardo},
        journal={ACM Transactions on Embedded Computing Systems (TECS)},
        volume={20},
        number={5s},
        pages={1--25},
        year={2021},
        publisher={ACM New York, NY}
    }

SA Tool:

    @inproceedings{carvalho2020design,
        title={A design exploration of scalable mesh-based fully pipelined accelerators},
        author={Carvalho, Westerley and Canesche, Michael and Reis, Lucas and Torres, Frank and Silva, Lucas and Jamieson, Peter and Nacif, Jos{\'e} and Ferreira, Ricardo},
        booktitle={2020 International Conference on Field-Programmable Technology (ICFPT)},
        pages={233--236},
        year={2020},
        organization={IEEE}
    }