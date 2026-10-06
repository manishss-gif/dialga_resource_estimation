#!/usr/bin/env bash
set -e

echo "=========================================================="
echo "RUNNING SHARED GROVER COMMON REGRESSION TESTS"
echo "=========================================================="
python3 test_dialaga_grover_common.py

echo "=========================================================="
echo "RUNNING REVERSIBLE ENCRYPTION KAT TESTS (ALL 8 PROFILES)"
echo "=========================================================="
python3 test_ct_profiles.py

echo "=========================================================="
echo "RUNNING COMPLETE GROVER FUNCTIONAL VALIDATION SUITES"
echo "=========================================================="
python3 dialaga128_ct_lw_full_grover.py --all-tests
python3 dialaga128_ct_ld_full_grover.py --all-tests
python3 dialaga128d_ct_lw_full_grover.py --all-tests
python3 dialaga128d_ct_ld_full_grover.py --all-tests

python3 dialaga256_ct_lw_full_grover.py --all-tests
python3 dialaga256_ct_ld_full_grover.py --all-tests
python3 dialaga256d_ct_lw_full_grover.py --all-tests
python3 dialaga256d_ct_ld_full_grover.py --all-tests

echo "=========================================================="
echo "ALL FUNCTIONAL VALIDATION TESTS PASSED SUCCESSFULLY"
echo "=========================================================="
