#!/usr/bin/env bash
set -e

echo "=========================================================="
echo "RUNNING ALL DIALGA GROVER RESOURCE ESTIMATIONS"
echo "=========================================================="

mkdir -p results

python3 dialga128_ct_lw_full_grover.py --resource | tee results/dialga128_lw.txt
python3 dialga128_ct_ld_full_grover.py --resource | tee results/dialga128_ld.txt
python3 dialga128d_ct_lw_full_grover.py --resource | tee results/dialga128d_lw.txt
python3 dialga128d_ct_ld_full_grover.py --resource | tee results/dialga128d_ld.txt

python3 dialga256_ct_lw_full_grover.py --resource | tee results/dialga256_lw.txt
python3 dialga256_ct_ld_full_grover.py --resource | tee results/dialga256_ld.txt
python3 dialga256d_ct_lw_full_grover.py --resource | tee results/dialga256d_lw.txt
python3 dialga256d_ct_ld_full_grover.py --resource | tee results/dialga256d_ld.txt

echo "=========================================================="
echo "ALL RESOURCE ESTIMATIONS COMPLETED SUCCESSFULLY"
echo "Results saved in results/"
echo "=========================================================="
