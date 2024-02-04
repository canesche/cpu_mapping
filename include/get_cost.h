#ifndef __GET_COST_H
#define __GET_COST_H

int cost_local(int pos_a_i, int pos_a_j, int pos_b_i, int pos_b_j, int GRID_SIZE, int arch) {
    if (pos_a_i == pos_b_i && pos_a_j == pos_b_j) return 1;

    int diff_i = abs(pos_a_i - pos_b_i);
    int diff_j = abs(pos_a_j - pos_b_j);

    if (arch == 0) {
        return max(1, diff_i + diff_j);
    } else if (arch == 1) {
        return max(1, diff_i/2 + diff_i%2 + diff_j/2 + diff_j%2);
    } 
    return 1000; // error 
}

#endif