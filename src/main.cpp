// #define ARCH 0  // 0: MESH, 1: 1-hop, 2: CHESS, 3: hex

#include <main.h>

int main(int argc, char* argv[]) {
    int value = 1656949570; //time(NULL);
    //printf("Seed %d\n", value);
    srand (value);
    
    string bench = "";
    unsigned int F_ID = 0, LIMIT = 0, ARCH = 0, NGRIDS = 1000, program = 0, TEMP = 100;
    
    if (argc > 4) {
        bench = argv[1]; // benchmark
        NGRIDS = stoi(argv[2]); // number of placement
        program = stoi(argv[3]); // algorithm: yoto, yott, sa
        if (program < 0 || program > 2) {
            cerr << "ERROR: Algorithm not supported!\n Yoto: 0, YOTT: 1, SA: 2" << endl;
            return 1;
        }
        ARCH = stoi(argv[4]); // architecture: mesh, 1-hop
        if (ARCH < 0 || ARCH > 1) {
            cerr << "ERROR: Architecture not supported!\n Mesh: 0, 1-hop: 1" << endl;
            return 1;
        }
    } else {
        cout << "ERROR: Number of parameters isn't enough!" << endl;
        return 1;
    }

    switch (program) {
        case 0:
            yolt_main(bench, NGRIDS, ARCH);
            break;
        case 1:
            yott_main(bench, NGRIDS, ARCH);
            break;
        case 2:
            sa_main(bench, NGRIDS, ARCH, LIMIT, F_ID, TEMP);
            break;
        default:
            break;
    }
    
    return 0;
}
