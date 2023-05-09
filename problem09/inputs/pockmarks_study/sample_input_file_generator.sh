# bash setup
set -x
set -e

id=$1

#rm -rf ${id}*.ini

#generate input files
generate_input_file(){
echo "
[grid.yasp]
LZ=25.
NZ=250

[before_storm]
dt_initial=60
output_time_interval=3600

[time]
dt_initial=2
time_end=79200.

[adaptive_time_control]
flag=false
dt_min=1.e-6
dt_max=1.
max_newton_steps=6
min_newton_steps=3

[output]
path_name=tests202007/set20200707
file_name=${id}${1}
time_interval=450

[newton]
max_iterations=12
abs_error=1.e-5

[initial]
nxch4=${7}

[free_gas_pocket]
z=${2}
dz=5.
Sg_max=${3}

[sediment]
number_of_layers=3

[sediment.layer0]
name=background
z=${2}
por=0.5
K=${4}
pentry=5000.
lambda=1.2
swr=0.
sgr=0.
beta=1.

[sediment.layer1]
name=barrier
z="$( bc <<<"${2} +0.5" )"
por=0.5
K=${4}
pentry=${5}
lambda=1.2
swr=0.
sgr=0.
beta=1.

[sediment.layer2]
name=background
z=25.
por=0.5
K=${4}
pentry=5000.
lambda=1.2
swr=0.
sgr=0.
beta=1.

[water_column]
wave_start_time=43200.
average_height=${6}
wave_amplitude=10.
wave_period=12.

[reference_state]
salinity=0.
temperature=10.

[gravity]
flag=true
magnitude=9.81
" > ${id}${1}.ini
}

##################
#   PARAMETERS   
##################
ns=0
#pe=30000.0
for z in 23.0 20.0 15.0; do
  for sg in 0.1 0.3 0.5; do
    for K in 1e-10 1e-12 1e-14; do
      for pe in 30000.0; do
        for H in 25.0 40.0 60.0; do
          for nxch4 in 0.30 0.60 0.90; do
	    ns=$((ns+1));
	    echo "Generating input file for Scenario ${ns}:"
	    generate_input_file "$ns" "$z" "$sg" "$K" "$pe" "$H" "$nxch4"
          done
        done
      done
    done
  done
done
