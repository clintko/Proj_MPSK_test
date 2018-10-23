N=$1
SEED=$2

treatments=(01 02 03 05 06 07 09 10 11)
patients=(F6901PRY_C G6904VJT_E K6902C85_A)

for trt in "${treatments[@]}"; do
    for patient in "${patients[@]}"; do
        out=mpsk_ep3_${trt}_${patient}_SEED${SEED}_N${N}
        echo $out
        Rscript mpsk_batch_ep3.R -o ${out} -t ${trt} -p ${patient} -n ${N} -s ${SEED} > ${out}.txt 2>&1
    done
done