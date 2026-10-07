"""
Regression and correctness tests comparing common Grover utilities
against the AND-tree implementation in dialga128_ct_lw_grover.py.
"""

import mpmath as mp
from projectq import MainEngine
from projectq.backends import ClassicalSimulator
from projectq.ops import Measure, All, X

# Import base functions under aliases
from dialga128_ct_lw_grover import (
    build_and_tree as base_build_and_tree,
    uncompute_and_tree as base_uncompute_and_tree
)

# Import common functions under aliases
from dialga_grover_common import (
    build_and_tree as common_build_and_tree,
    uncompute_and_tree as common_uncompute_and_tree,
    grover_diffusion,
    get_grover_iterations,
    grover_success_probability,
    expected_spurious_keys,
    format_base_2,
    calculate_exhaustive_totals,
    test_small_diffusion,
    test_toy_complete_grover
)


# =====================================================================
# Regression Tests F1 through F9
# =====================================================================


from test_advanced_counter import AdvancedResourceCounter

mp.mp.dps = 100


def test_and_tree_functional():
    print("F1: AND Tree Functional Comparison...")
    # Test for small sizes
    sizes = [2, 3, 4, 5, 8, 16]
    sim = ClassicalSimulator()
    eng = MainEngine(sim, [])
    
    for n in sizes:
        # We test several input patterns (e.g. all 0s, all 1s, and a few intermediate ones)
        patterns = [0, (1 << n) - 1]
        if n > 2:
            patterns.append(1)
            patterns.append((1 << (n - 1)) - 1)
            patterns.append(1 << (n - 1))
            
        for pattern in patterns:
            # Test Base AND Tree
            inputs_base = eng.allocate_qureg(n)
            cmp_ancillas_base = eng.allocate_qureg(n - 1)
            cmp_helpers_base = eng.allocate_qureg(n // 2)
            
            # Load pattern
            for j in range(n):
                if (pattern >> j) & 1:
                    X | inputs_base[j]
            eng.flush()
            
            root_base, history_base = base_build_and_tree(eng, inputs_base, cmp_ancillas_base, cmp_helpers_base)
            eng.flush()
            
            expected_val = 1 if (pattern == (1 << n) - 1) else 0
            val_base = sim.read_bit(root_base)
            assert val_base == expected_val, f"Base AND tree error for n={n}, pattern={pattern}: expected {expected_val}, got {val_base}"
            
            base_uncompute_and_tree(history_base)
            eng.flush()
            
            # Assert all ancillas and helpers are cleaned up
            for q in cmp_ancillas_base:
                assert sim.read_bit(q) == 0, "Base AND tree ancilla dirty after uncomputation"
            for q in cmp_helpers_base:
                assert sim.read_bit(q) == 0, "Base AND tree helper dirty after uncomputation"
            for j in range(n):
                expected_input = (pattern >> j) & 1
                assert sim.read_bit(inputs_base[j]) == expected_input, "Base AND tree input corrupted"
                
            All(Measure) | inputs_base
            All(Measure) | cmp_ancillas_base
            All(Measure) | cmp_helpers_base
            eng.flush()
            
            # Test Common AND Tree
            inputs_common = eng.allocate_qureg(n)
            cmp_ancillas_common = eng.allocate_qureg(n - 1)
            cmp_helpers_common = eng.allocate_qureg(n // 2)
            
            for j in range(n):
                if (pattern >> j) & 1:
                    X | inputs_common[j]
            eng.flush()
            
            root_common, history_common = common_build_and_tree(eng, inputs_common, cmp_ancillas_common, cmp_helpers_common)
            eng.flush()
            
            val_common = sim.read_bit(root_common)
            assert val_common == expected_val, f"Common AND tree error for n={n}, pattern={pattern}: expected {expected_val}, got {val_common}"
            
            common_uncompute_and_tree(history_common)
            eng.flush()
            
            for q in cmp_ancillas_common:
                assert sim.read_bit(q) == 0, "Common AND tree ancilla dirty after uncomputation"
            for q in cmp_helpers_common:
                assert sim.read_bit(q) == 0, "Common AND tree helper dirty after uncomputation"
            for j in range(n):
                expected_input = (pattern >> j) & 1
                assert sim.read_bit(inputs_common[j]) == expected_input, "Common AND tree input corrupted"
                
            All(Measure) | inputs_common
            All(Measure) | cmp_ancillas_common
            All(Measure) | cmp_helpers_common
            eng.flush()
            
    print("  AND Tree Functional Comparison: SUCCESS\n")


def test_and_tree_resources():
    print("F2: AND Tree Resource Comparison...")
    
    for n in [256, 384]:
        # Measure Base AND Tree resources
        rc_base = AdvancedResourceCounter()
        rc_base.is_recording = True
        eng_base = MainEngine(rc_base, [])
        inputs_base = eng_base.allocate_qureg(n)
        cmp_ancillas_base = eng_base.allocate_qureg(n - 1)
        cmp_helpers_base = eng_base.allocate_qureg(n // 2)
        
        root_base, history_base = base_build_and_tree(eng_base, inputs_base, cmp_ancillas_base, cmp_helpers_base)
        base_uncompute_and_tree(history_base)
        eng_base.flush()
        rc_base.is_recording = False
        
        # Measure Common AND Tree resources
        rc_common = AdvancedResourceCounter()
        rc_common.is_recording = True
        eng_common = MainEngine(rc_common, [])
        inputs_common = eng_common.allocate_qureg(n)
        cmp_ancillas_common = eng_common.allocate_qureg(n - 1)
        cmp_helpers_common = eng_common.allocate_qureg(n // 2)
        
        root_common, history_common = common_build_and_tree(eng_common, inputs_common, cmp_ancillas_common, cmp_helpers_common)
        common_uncompute_and_tree(history_common)
        eng_common.flush()
        rc_common.is_recording = False
        
        # Compare gate counts and depths
        for metric in ["cnot", "t", "measure", "oneq", "t_depth", "full_depth"]:
            val_base = rc_base.get_cnot() if metric == "cnot" else \
                       rc_base.get_t() if metric == "t" else \
                       rc_base.get_measure() if metric == "measure" else \
                       rc_base.get_1q_clifford() if metric == "oneq" else \
                       rc_base.t_depth if metric == "t_depth" else \
                       rc_base.pure_circuit_depth
            
            val_common = rc_common.get_cnot() if metric == "cnot" else \
                         rc_common.get_t() if metric == "t" else \
                         rc_common.get_measure() if metric == "measure" else \
                         rc_common.get_1q_clifford() if metric == "oneq" else \
                         rc_common.t_depth if metric == "t_depth" else \
                         rc_common.pure_circuit_depth
                         
            assert val_base == val_common, f"AND tree resource mismatch for n={n}, metric={metric}: base={val_base}, common={val_common}"
            
        # Assert workspace sizes
        if n == 384:
            assert len(history_common) == 383, f"Expected 383 AND nodes, got {len(history_common)}"
            assert len(cmp_helpers_common) == 192, f"Expected 192 helpers, got {len(cmp_helpers_common)}"
        elif n == 256:
            assert len(history_common) == 255, f"Expected 255 AND nodes, got {len(history_common)}"
            assert len(cmp_helpers_common) == 128, f"Expected 128 helpers, got {len(cmp_helpers_common)}"
            
        # Cleanup
        All(Measure) | inputs_base
        All(Measure) | cmp_ancillas_base
        All(Measure) | cmp_helpers_base
        eng_base.flush()
        
        All(Measure) | inputs_common
        All(Measure) | cmp_ancillas_common
        All(Measure) | cmp_helpers_common
        eng_common.flush()
        
    print("  AND Tree Resource Comparison: SUCCESS\n")


def test_diffusion_resources():
    print("F3: Diffusion Comparison...")
    
    # Check 256-bit diffusion resources
    rc = AdvancedResourceCounter()
    rc.is_recording = True
    eng = MainEngine(rc, [])
    
    key_bits = eng.allocate_qureg(256)
    cmp_ancillas = eng.allocate_qureg(255)
    cmp_helpers = eng.allocate_qureg(128)
    
    grover_diffusion(eng, key_bits, cmp_ancillas, cmp_helpers)
    eng.flush()
    rc.is_recording = False
    
    cnot = rc.get_cnot()
    oneq = rc.get_1q_clifford()
    t = rc.get_t()
    measure = rc.get_measure()
    t_depth = rc.t_depth
    full_depth = rc.pure_circuit_depth
    
    # Assert against golden baseline
    assert cnot == 2550, f"Diffusion CNOT mismatch: got {cnot}, expected 2550"
    assert oneq == 3065, f"Diffusion 1qClifford mismatch: got {oneq}, expected 3065"
    assert t == 1020, f"Diffusion T mismatch: got {t}, expected 1020"
    assert measure == 255, f"Diffusion Measure mismatch: got {measure}, expected 255"
    assert t_depth == 8, f"Diffusion T-depth mismatch: got {t_depth}, expected 8"
    assert full_depth == 63, f"Diffusion Full-depth mismatch: got {full_depth}, expected 63"
    
    All(Measure) | key_bits
    All(Measure) | cmp_ancillas
    All(Measure) | cmp_helpers
    eng.flush()
    
    # Run functional correctness tests (Simulator-based)
    test_small_diffusion()
    print("  Diffusion Comparison: SUCCESS\n")


def test_toy_grover_regression():
    print("F4: Toy Complete Grover Regression...")
    test_toy_complete_grover()
    print("  Toy Complete Grover Regression: SUCCESS\n")


def test_grover_iteration_count():
    print("F5: Grover Iteration Count Validation...")
    N_G = get_grover_iterations()
    assert N_G == 267257146016241686964920093290467695825, f"Iteration count mismatch: {N_G}"
    print("  Grover Iteration Count Validation: SUCCESS\n")


def test_spurious_keys():
    print("F6: Spurious-Key Function Validation...")
    spurious = expected_spurious_keys(256, 3, 128)
    assert spurious != 0, "Spurious keys must not be integer zero"
    assert spurious == mp.mpf(2) ** (-128), f"Spurious keys value mismatch: {spurious}"
    print("  Spurious-Key Function Validation: SUCCESS\n")


def test_format_base_2_regression():
    print("F7: format_base_2 Regression...")
    
    # Positive integer
    s1 = format_base_2(2**146)
    assert s1.strip() == "1.0000 * 2^146", f"format_base_2(2**146) failed: {s1}"
    
    # Positive float/mpf less than 1
    s2 = format_base_2(mp.mpf(2)**-128)
    assert s2.strip() == "1.0000 * 2^-128", f"format_base_2(2**-128) failed: {s2}"
    
    # Zero
    s3 = format_base_2(0)
    assert s3.strip() == "0", f"format_base_2(0) failed: {s3}"
    
    print("  format_base_2 Regression: SUCCESS\n")


def test_exhaustive_scaling_regression():
    print("F8: Exhaustive Scaling Regression...")
    prep = {
        "cnot": 2,
        "oneq": 3,
        "t": 5,
        "measure": 7,
        "t_depth": 11,
        "full_depth": 13,
        "width": 17
    }
    iteration = {
        "cnot": 19,
        "oneq": 23,
        "t": 29,
        "measure": 31,
        "t_depth": 37,
        "full_depth": 41,
        "width": 43
    }
    N = 5
    readout_measure = 256
    readout_full_depth = 1

    totals = calculate_exhaustive_totals(
        prep, iteration, readout_measure, readout_full_depth,
        N_G=N, delta_t=37, delta_full=41
    )
    
    assert totals['cnot'] == 2 + 5 * 19
    assert totals['oneq'] == 3 + 5 * 23
    assert totals['t'] == 5 + 5 * 29
    assert totals['measure'] == 7 + 5 * 31 + 256
    assert totals['t_depth_cons'] == 11 + 5 * 37
    assert totals['full_depth_cons'] == 13 + 5 * 41 + 1
    assert totals['width'] == max(17, 43, 256)
    
    print("  Exhaustive Scaling Regression: SUCCESS\n")


def test_grover_success_probability_regression():
    print("F9: Grover Success Probability Regression...")
    N_G = get_grover_iterations()
    p = grover_success_probability(N_G)
    assert 0 < p <= 1
    assert (1 - p) < mp.mpf("1e-77")
    print("  Grover Success Probability Regression: SUCCESS\n")


def run_all_tests():
    print("==============================================================")
    print("TEST SUITE: DIALGA GROVER COMMON VS BASE REGRESSION")
    print("==============================================================")
    test_and_tree_functional()
    test_and_tree_resources()
    test_diffusion_resources()
    test_toy_grover_regression()
    test_grover_iteration_count()
    test_spurious_keys()
    test_format_base_2_regression()
    test_exhaustive_scaling_regression()
    test_grover_success_probability_regression()
    print("==============================================================")
    print("ALL TESTS PASSED SUCCESSFULLY!")
    print("==============================================================")


if __name__ == "__main__":
    run_all_tests()
