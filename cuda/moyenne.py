import sys
with open("exec_time.data") as f :
    sum = 0
    for line in f.readlines() :
        sum += float(line);
    print("Average execution time for nbody_brute_force with {} particles and TFINAL = {} was {}".format(sys.argv[1], sys.argv[2], sum))
