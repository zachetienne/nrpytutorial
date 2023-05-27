#!/usr/bin/env bash

config=$(readlink -f "$1")
bin=$(readlink -f "$2")
output=${bin%/bin/*}/output

row_counter=0
while IFS=, read -r nx0 nx1 nx2 cfl_factor; do
    row_counter=$((row_counter+1))
    output_target=$output/run-$row_counter
    [[ -e $output_target ]] && mv --no-target-directory --force --backup=t $output_target $output_target.bak
    mkdir $output_target
    echo -e "\n$row_counter is running with parm set: $nx0 $nx1 $nx2 $cfl_factor"
    echo "NX0=$nx0", "NX1=$nx1", "NX2=$nx2", "CFL_FACTOR=$cfl_factor" > $output_target/run.config
    cp $output/SFID.txt $output_target/.
    (cd $output_target && $bin $nx0 $nx1 $nx2 $cfl_factor)
    echo -e "Succeed\n"
done < <(tr -d "[:blank:]" < $config | grep -v '#')
