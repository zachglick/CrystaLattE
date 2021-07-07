import os, sys
import numpy as np


if __name__ == '__main__':

    lines_zach = open('components_mixed.csv', 'r').readlines()
    lines_zach = [line.split(',') for line in lines_zach]
    e_elst = 0.0
    for line in lines_zach[1:]:
        e_elst += int(line[1]) * float(line[2])
    print(e_elst * 2625.5 / 2.0)

    e_elst2 = 0.0
    for line in lines_zach[1:]:
        key = line[0]
        lines_carlos = open(f'/theoryfs2/ds/glick/FromCarlos/libefp/Benzene/{key}.out', 'r').readlines()
        #print(lines_carlos)
        for line in lines_carlos:
            if line.strip().startswith('ELECTROSTATIC ENERGY'):
            #if line.strip().startswith('CHARGE PENETRATION ENERGY'):
                e_elst2 += float(line.strip().split()[-1])
                print(key, float(line.strip().split()[-1]) * 2625.5 / 2.0)
    print(e_elst2 * 2625.5 / 2.0)
