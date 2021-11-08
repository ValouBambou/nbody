nparticles=$1
tfinal = $2
for i in {1..100}
rm exec_time.data
do 
    echo "Running nbody_brute_force with $nparticles particles."
    exec_time = $(./nbody_brute_force $1 $2| grep 'Simulation' | cut -f3 -d ' ')
    touch exec_time.data
    printf ("$exec_time\n") >> exec_time.data
done
python3 ./moyenne.py $1 $2


