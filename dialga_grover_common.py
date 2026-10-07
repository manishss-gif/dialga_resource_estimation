"""
Common Grover framework utilities for the Dialga family.
This module contains functionality shared across all Dialga variants (128, 128†, 256, 256†),
such as diffusion, base-2 formatting, N_G calculation, correctness tests, and resource scaling.
"""

import sys
import mpmath as mp
from projectq import MainEngine
from projectq.backends import ClassicalSimulator, Simulator
from projectq.ops import X, Measure, All, Z, CNOT, H
from corefunctions import and_gate, and_gate_dagger

mp.mp.dps = 100


# =====================================================================
# AND-Tree primitives
# =====================================================================

def build_and_tree(eng, current_level, cmp_ancillas, cmp_helpers):
    num_inputs = len(current_level)
    assert num_inputs >= 2, "build_and_tree requires at least 2 input qubits"
    assert len(cmp_ancillas) >= num_inputs - 1, f"build_and_tree requires at least {num_inputs - 1} ancillas, got {len(cmp_ancillas)}"
    assert len(cmp_helpers) >= num_inputs // 2, f"build_and_tree requires at least {num_inputs // 2} helpers, got {len(cmp_helpers)}"
    
    history = []
    anc_idx = 0
    level_nodes = current_level[:]
    
    while len(level_nodes) > 1:
        next_level = []
        hlp_idx = 0
        for idx in range(0, len(level_nodes) - 1, 2):
            a = level_nodes[idx]
            b = level_nodes[idx + 1]
            c = cmp_ancillas[anc_idx]
            anc_idx += 1
            d = cmp_helpers[hlp_idx]
            hlp_idx += 1
            
            and_gate(a, b, c, d)
            history.append((a, b, c))
            next_level.append(c)
            
        if len(level_nodes) % 2 == 1:
            next_level.append(level_nodes[-1])
            
        level_nodes = next_level
        
    assert len(history) == num_inputs - 1, f"Expected {num_inputs - 1} AND nodes, got {len(history)}"
    return level_nodes[0], history


def uncompute_and_tree(history):
    for a, b, c in reversed(history):
        and_gate_dagger(a, b, c)


# =====================================================================
# Grover Diffusion Operator
# =====================================================================

def grover_diffusion(eng, key_bits, cmp_ancillas, cmp_helpers):
    """
    Implements the Grover diffusion transform on key_bits:
    H^⊗k -> X^⊗k -> phase flip on |1...1> -> X^⊗k -> H^⊗k
    This represents inversion-about-the-mean up to a global phase.
    """
    k = len(key_bits)
    assert k >= 2, "grover_diffusion requires at least 2 key bits"
    assert len(cmp_ancillas) >= k - 1, f"grover_diffusion requires at least {k - 1} ancillas, got {len(cmp_ancillas)}"
    assert len(cmp_helpers) >= k // 2, f"grover_diffusion requires at least {k // 2} helpers, got {len(cmp_helpers)}"
    
    is_classical = eng.backend.__class__.__name__ == 'ClassicalSimulator'
    
    # 1. Hadamards and X gates (skipped on ClassicalSimulator)
    if not is_classical:
        for q in key_bits:
            H | q
        for q in key_bits:
            X | q
        
    # 3. Balanced AND tree on the key bits to check if they are all 1.
    root, history = build_and_tree(eng, key_bits, cmp_ancillas[:k-1], cmp_helpers[:k//2])
    
    # 4. Phase flip on root
    if is_classical:
        pass
    else:
        Z | root
    
    # 5. Uncompute AND tree
    uncompute_and_tree(history)
    
    # 6. X gates and Hadamards (skipped on ClassicalSimulator)
    if not is_classical:
        for q in key_bits:
            X | q
        for q in key_bits:
            H | q


# =====================================================================
# Mathematical Calculations & Scaling
# =====================================================================

def get_grover_iterations():
    """Returns the exact N_G = floor((pi/4) * 2^128)."""
    N_G = int(mp.floor((mp.pi / 4) * (mp.mpf(2) ** 128)))
    assert N_G == 267257146016241686964920093290467695825
    return N_G


def grover_success_probability(N_G):
    """Computes the success probability of Grover search for k=256."""
    theta = mp.asin(mp.mpf(1) / mp.sqrt(mp.mpf(2) ** 256))
    p_success = mp.sin((2 * N_G + 1) * theta) ** 2
    return p_success


def expected_spurious_keys(k, r, n):
    """Computes the expected number of spurious matching keys."""
    return mp.mpf(2) ** (k - r * n)


def format_base_2(val):
    """Helper to format large integers and fractions in base-2 scientific notation (mantissa * 2^exponent)."""
    if val == 0:
        return "0"
    
    # Convert input to arbitrary precision mpmath float
    val_mp = mp.mpf(val)
    
    # Calculate exponent: exp = floor(log2(val))
    exp = int(mp.floor(mp.log(val_mp, 2)))
    
    # Mantissa: mantissa = val / 2^exp
    mantissa = val_mp / (mp.mpf(2) ** exp)
    
    return f"{mantissa:.4f} * 2^{exp}"


def calculate_exhaustive_totals(prep, iter_res, readout_measure, readout_full_depth, N_G, delta_t=None, delta_full=None):
    """Scales one-iteration and setup resources to complete key-search totals."""
    assert N_G >= 0, "N_G must be non-negative"
    assert readout_measure >= 0, "readout_measure must be non-negative"
    assert readout_full_depth >= 0, "readout_full_depth must be non-negative"
    
    # Setup added once, iteration scaled by N_G
    cnot_total = prep["cnot"] + N_G * iter_res["cnot"]
    oneq_total = prep["oneq"] + N_G * iter_res["oneq"]
    t_total = prep["t"] + N_G * iter_res["t"]
    
    # Readout added once
    measure_total = prep["measure"] + N_G * iter_res["measure"] + readout_measure
    
    conservative_t_depth = prep["t_depth"] + N_G * iter_res["t_depth"]
    conservative_full_depth = prep["full_depth"] + N_G * iter_res["full_depth"] + readout_full_depth
    
    # Widths treated as maximum (non-additive)
    width_total = max(prep["width"], iter_res["width"], 256)
    
    affine_t_depth = None
    if delta_t is not None:
        affine_t_depth = prep["t_depth"] + iter_res["t_depth"] + (N_G - 1) * delta_t
        
    affine_full_depth = None
    if delta_full is not None:
        affine_full_depth = prep["full_depth"] + iter_res["full_depth"] + (N_G - 1) * delta_full + readout_full_depth
        
    return {
        "cnot": cnot_total,
        "oneq": oneq_total,
        "t": t_total,
        "measure": measure_total,
        "t_depth_cons": conservative_t_depth,
        "full_depth_cons": conservative_full_depth,
        "t_depth_affine": affine_t_depth,
        "full_depth_affine": affine_full_depth,
        "width": width_total
    }


# =====================================================================
# Functional Correctness Validation Tests
# =====================================================================

def get_ideal_diffusion_matrix(n):
    """Generates the ideal Grover diffusion matrix D = 2|s><s| - I."""
    size = 1 << n
    D = [[0.0 for _ in range(size)] for _ in range(size)]
    for r in range(size):
        for c in range(size):
            val = 2.0 / size
            if r == c:
                val -= 1.0
            D[r][c] = val
    return D


def get_key_state(sim, key_bits, ancillas, helpers):
    """Helper to extract the state vector of the key bits from the simulator."""
    qubit_map, state_vector = sim.cheat()
    n = len(key_bits)
    key_amplitudes = [0.0] * (1 << n)
    
    key_indices = [qubit_map[q.id] for q in key_bits]
    ancilla_indices = [qubit_map[q.id] for q in ancillas]
    helper_indices = [qubit_map[q.id] for q in helpers]
    
    for k_state in range(len(state_vector)):
        amp = state_vector[k_state]
        if abs(amp) < 1e-9:
            continue
        # Verify all ancilla and helper qubits are in state |0>
        for idx in ancilla_indices + helper_indices:
            if (k_state >> idx) & 1:
                raise ValueError(f"Ancilla/Helper qubit at index {idx} in state vector is not 0 (amplitude {amp})")
        # Extract key value
        key_val = 0
        for j, idx in enumerate(key_indices):
            bit = (k_state >> idx) & 1
            key_val |= (bit << j)
        key_amplitudes[key_val] = amp
    return key_amplitudes


def test_small_diffusion():
    """Validates diffusion on small systems (2, 3, 4 qubits) against the ideal matrix."""
    print("Running small system diffusion tests...")
    import projectq.backends._sim._pysim as pysim
    pysim._USE_REFCHECK = False

    for n in [2, 3, 4]:
        size = 1 << n
        ideal_D = get_ideal_diffusion_matrix(n)
        measured_M = [[0.0 for _ in range(size)] for _ in range(size)]
        
        for c in range(size):
            sim = Simulator()
            eng = MainEngine(sim, [])
            
            key_bits = eng.allocate_qureg(n)
            cmp_ancillas = eng.allocate_qureg(n - 1)
            cmp_helpers = eng.allocate_qureg(n // 2)
            
            # Load basis state c
            for j in range(n):
                if (c >> j) & 1:
                    X | key_bits[j]
            eng.flush()
            
            # Run diffusion
            grover_diffusion(eng, key_bits, cmp_ancillas, cmp_helpers)
            eng.flush()
            
            key_state = get_key_state(sim, key_bits, cmp_ancillas, cmp_helpers)
            for r in range(size):
                measured_M[r][c] = key_state[r]
                
            # Cleanup for engine destruction
            All(Measure) | key_bits
            All(Measure) | cmp_ancillas
            All(Measure) | cmp_helpers
            eng.flush()
            
        # Find global phase factor by comparing entry (0, 0)
        phase_factor = measured_M[0][0] / ideal_D[0][0]
        assert abs(abs(phase_factor) - 1.0) < 1e-5, f"Phase factor magnitude is not 1 for n={n}: {abs(phase_factor)}"
        
        # Verify measured matrix matches ideal matrix up to the global phase factor
        for r in range(size):
            for c in range(size):
                diff = abs(measured_M[r][c] - phase_factor * ideal_D[r][c])
                assert diff < 1e-5, f"Matrix mismatch at ({r}, {c}) for n={n}: measured={measured_M[r][c]}, expected={phase_factor * ideal_D[r][c]}, diff={diff}"
                
        print(f"  Diffusion verification for n={n} qubits: SUCCESS")
    print("All diffusion correctness tests passed!\n")


def toy_phase_oracle(eng, key_bits, target, cmp_ancillas, cmp_helpers, test_flag=None):
    """Toy oracle that marks the target state in a 3-qubit key space."""
    for j in range(3):
        if not ((target >> j) & 1):
            X | key_bits[j]
            
    root, history = build_and_tree(eng, key_bits, cmp_ancillas[:2], cmp_helpers[:1])
    
    if test_flag is not None:
        CNOT | (root, test_flag)
    else:
        Z | root
        
    uncompute_and_tree(history)
    for j in range(3):
        if not ((target >> j) & 1):
            X | key_bits[j]


def test_toy_complete_grover():
    """Runs a complete 3-qubit Grover search to verify amplitude amplification success probability."""
    print("Running Toy Complete Grover test (3 qubits)...")
    import projectq.backends._sim._pysim as pysim
    pysim._USE_REFCHECK = False

    target = 5  # target state |101>
    
    sim = Simulator()
    eng = MainEngine(sim, [])
    
    key_bits = eng.allocate_qureg(3)
    cmp_ancillas = eng.allocate_qureg(2)
    cmp_helpers = eng.allocate_qureg(1)
    
    # 1. Setup uniform superposition
    for q in key_bits:
        H | q
    eng.flush()
    
    # 2. Run 2 Grover iterations (optimal for 3 qubits, 1 target)
    for _ in range(2):
        toy_phase_oracle(eng, key_bits, target, cmp_ancillas, cmp_helpers)
        grover_diffusion(eng, key_bits, cmp_ancillas, cmp_helpers)
        eng.flush()
        
    # Read amplitudes
    key_state = get_key_state(sim, key_bits, cmp_ancillas, cmp_helpers)
    
    # Cleanup for engine destruction
    All(Measure) | key_bits
    All(Measure) | cmp_ancillas
    All(Measure) | cmp_helpers
    eng.flush()
    
    target_prob = abs(key_state[target]) ** 2
    expected_prob = 121.0 / 128.0
    print(f"  Target state: |{target:03b}>")
    print(f"  Measured target probability: {target_prob:.6f}")
    print(f"  Expected target probability: {expected_prob:.6f}")
    
    assert abs(target_prob - expected_prob) < 1e-5, f"Toy Grover probability mismatch: got {target_prob}, expected {expected_prob}"
    print("  Toy Grover complete search validation: SUCCESS\n")


# =====================================================================
# Driver shared by the eight *_full_grover.py scripts
# =====================================================================

K_BITS = 256
N_BITS = 128
R_TUPLES = 3
COMPARED_BITS = R_TUPLES * N_BITS
KEY_HEX = "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00"


class GroverVariant:
    """Cipher variant, circuit profile and reference values for one script.

    tuples           three (plaintext, tweak blocks, ciphertext) under KEY_HEX
    allocate_work    returns the S-box work qubits of one encryption lane
    encrypt          encrypt(eng, state, tweaks, key, work)
    encrypt_inverse  same arguments, undoes encrypt
    oracle           oracle(eng, keys, states, tweaks, works,
                            cmp_ancillas, cmp_helpers, ciphertexts, test_flag)
    diffusion_workspace  (product ancillas, helpers) given to D_256 when it is
                         measured on its own
    reference        expected (CNOT, 1qClifford, T, Measure, T-depth, depth, width)
                     for "encryption", "oracle", "diffusion" and "iteration"
    """

    def __init__(self, name, profile, tuples, allocate_work, encrypt,
                 encrypt_inverse, oracle, diffusion_workspace, reference):
        self.name = name
        self.profile = profile
        self.title = f"{name} {profile}"
        self.tuples = tuples
        self.num_tweaks = len(tuples[0][1])
        self.ciphertexts = [c for _, _, c in tuples]
        self.allocate_work = allocate_work
        self.encrypt = encrypt
        self.encrypt_inverse = encrypt_inverse
        self.oracle = oracle
        self.diffusion_workspace = diffusion_workspace
        self.reference = reference


def columns(per_lane):
    """[[a1, b1], [a2, b2], [a3, b3]] -> [[a1, a2, a3], [b1, b2, b3]]"""
    return [list(c) for c in zip(*per_lane)]


def _allocate_lanes(eng, v):
    from corefunctions import allocate_cells
    return {
        "keys": [allocate_cells(eng, 32) for _ in range(R_TUPLES)],
        "states": [allocate_cells(eng, 16) for _ in range(R_TUPLES)],
        "tweaks": [[allocate_cells(eng, 16) for _ in range(v.num_tweaks)] for _ in range(R_TUPLES)],
        "works": [v.allocate_work(eng) for _ in range(R_TUPLES)],
        "cmp_ancillas": eng.allocate_qureg(COMPARED_BITS - 1),
        "cmp_helpers": eng.allocate_qureg(COMPARED_BITS // 2),
    }


def _load_public(v, lanes):
    from corefunctions import load_hex
    for lane, (plaintext, tweaks, _) in enumerate(v.tuples):
        load_hex(plaintext, lanes["states"][lane])
        for reg, value in zip(lanes["tweaks"][lane], tweaks):
            load_hex(value, reg)


def _run_oracle(eng, v, lanes, test_flag=None):
    v.oracle(eng, lanes["keys"], lanes["states"], lanes["tweaks"], lanes["works"],
             lanes["cmp_ancillas"], lanes["cmp_helpers"], v.ciphertexts, test_flag)


def _run_iteration(eng, v, lanes, test_flag=None):
    _run_oracle(eng, v, lanes, test_flag)
    key0 = [q for cell in lanes["keys"][0] for q in cell]
    grover_diffusion(eng, key0, lanes["cmp_ancillas"], lanes["cmp_helpers"])


def _check_restored(sim, v, lanes, key0_hex):
    """Asserts that the oracle left every register as it found it."""
    from corefunctions import flatten, read_hex
    assert read_hex(sim, lanes["keys"][0]) == key0_hex, "candidate key changed"
    for lane, (plaintext, tweaks, _) in enumerate(v.tuples):
        assert read_hex(sim, lanes["states"][lane]) == plaintext, f"plaintext {lane + 1} changed"
        for reg, value in zip(lanes["tweaks"][lane], tweaks):
            assert read_hex(sim, reg) == value, f"tweak {lane + 1} changed"
    workspace = [lanes["keys"][1:], lanes["works"], lanes["cmp_ancillas"], lanes["cmp_helpers"]]
    assert all(sim.read_bit(q) == 0 for q in flatten(workspace)), "workspace not returned to 0"


def _release(eng, regs):
    from corefunctions import flatten
    for q in flatten(regs):
        Measure | q
    eng.flush()


# ---------------------------------------------------------------------
# Functional tests (computational-basis states)
# ---------------------------------------------------------------------

def test_encryption(v):
    """Encrypts the three oracle tuples and compares with the expected ciphertexts."""
    from corefunctions import allocate_cells, flatten, load_hex, read_hex
    print("Encryption of the three oracle tuples (tuple 1 is the official test vector)")
    for i, (plaintext, tweaks, ciphertext) in enumerate(v.tuples, 1):
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        state = allocate_cells(eng, 16)
        tweak_regs = [allocate_cells(eng, 16) for _ in tweaks]
        key = allocate_cells(eng, 32)
        work = v.allocate_work(eng)
        load_hex(plaintext, state)
        for reg, value in zip(tweak_regs, tweaks):
            load_hex(value, reg)
        load_hex(KEY_HEX, key)
        eng.flush()

        v.encrypt(eng, state, tweak_regs, key, work)
        eng.flush()
        result = read_hex(sim, state)
        assert result == ciphertext, f"tuple {i}: expected {ciphertext}, got {result}"
        assert all(sim.read_bit(q) == 0 for q in flatten(work)), f"tuple {i}: work qubits not returned to 0"
        print(f"    tuple {i}  {result}  PASS")
        _release(eng, [state, tweak_regs, key, work])


def test_round_trip(v):
    """Checks that ENC followed by ENC^dagger restores all registers."""
    from corefunctions import allocate_cells, flatten, load_hex, read_hex
    print("ENC followed by ENC^dagger")
    for i, (plaintext, tweaks, _) in enumerate(v.tuples, 1):
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        state = allocate_cells(eng, 16)
        tweak_regs = [allocate_cells(eng, 16) for _ in tweaks]
        key = allocate_cells(eng, 32)
        work = v.allocate_work(eng)
        load_hex(plaintext, state)
        for reg, value in zip(tweak_regs, tweaks):
            load_hex(value, reg)
        load_hex(KEY_HEX, key)
        eng.flush()

        v.encrypt(eng, state, tweak_regs, key, work)
        v.encrypt_inverse(eng, state, tweak_regs, key, work)
        eng.flush()
        assert read_hex(sim, state) == plaintext, f"tuple {i}: state not restored"
        assert [read_hex(sim, r) for r in tweak_regs] == list(tweaks), f"tuple {i}: tweak not restored"
        assert read_hex(sim, key) == KEY_HEX, f"tuple {i}: key changed"
        assert all(sim.read_bit(q) == 0 for q in flatten(work)), f"tuple {i}: work qubits not returned to 0"
        print(f"    tuple {i}  PASS")
        _release(eng, [state, tweak_regs, key, work])


def test_phase_oracle(v):
    """Correct key sets the predicate, a one-bit-off key does not; all registers restored."""
    from corefunctions import load_hex, read_hex
    print("Phase oracle predicate and cleanup")
    sim = ClassicalSimulator()
    eng = MainEngine(sim, [])
    lanes = _allocate_lanes(eng, v)
    flag = eng.allocate_qubit()[0]
    load_hex(KEY_HEX, lanes["keys"][0])
    _load_public(v, lanes)
    eng.flush()

    for label, flip, expected in (("correct key", False, 1), ("wrong key", True, 0)):
        if flip:
            X | lanes["keys"][0][0][0]
            eng.flush()
        key0_hex = read_hex(sim, lanes["keys"][0])
        _run_oracle(eng, v, lanes, flag)
        eng.flush()
        assert sim.read_bit(flag) == expected, f"{label}: predicate {sim.read_bit(flag)}, expected {expected}"
        _check_restored(sim, v, lanes, key0_hex)
        print(f"    {label:<12} predicate = {expected}, registers restored  PASS")
        if expected:
            X | flag
            eng.flush()
    _release(eng, [list(lanes.values()), flag])


def test_iteration(v):
    """One iteration (oracle + diffusion) leaves the workspace clean.

    On ClassicalSimulator the diffusion runs only its AND tree, so this checks
    register bookkeeping; the diffusion itself is tested by test_small_diffusion.
    """
    from corefunctions import load_hex
    print("One Grover iteration: workspace cleanup")
    sim = ClassicalSimulator()
    eng = MainEngine(sim, [])
    lanes = _allocate_lanes(eng, v)
    flag = eng.allocate_qubit()[0]
    load_hex(KEY_HEX, lanes["keys"][0])
    _load_public(v, lanes)
    eng.flush()

    _run_iteration(eng, v, lanes, flag)
    eng.flush()
    assert sim.read_bit(flag) == 1, "correct key not marked"
    _check_restored(sim, v, lanes, KEY_HEX)
    print("    PASS")
    _release(eng, [list(lanes.values()), flag])


def run_functional_tests(v):
    from test_dialga_grover_common import test_and_tree_functional, test_and_tree_resources
    print("=" * 70)
    print(f"{v.title}: functional tests")
    print("=" * 70)
    test_encryption(v)
    test_round_trip(v)
    test_phase_oracle(v)
    test_iteration(v)
    test_and_tree_functional()
    test_and_tree_resources()
    test_small_diffusion()
    test_toy_complete_grover()


# ---------------------------------------------------------------------
# Resource report
# ---------------------------------------------------------------------

def _check_reference(label, measured, expected):
    from test_advanced_counter import METRICS, as_tuple
    got = as_tuple(measured)
    if got != tuple(expected):
        rows = [f"{name}: {g} (expected {e})" for (_, name), g, e in zip(METRICS, got, expected) if g != e]
        raise AssertionError(f"{label} differs from the reference values: " + ", ".join(rows))


def _print_totals(totals, oneq_with_loading):
    rows = [
        ("#CNOT", totals["cnot"]),
        ("#1qClifford", totals["oneq"]),
        ("#T", totals["t"]),
        ("#Measure", totals["measure"]),
        ("T-depth", totals["t_depth_cons"]),
        ("Full depth", totals["full_depth_cons"]),
    ]
    if totals["t_depth_affine"] is not None:
        rows.append(("T-depth (affine)", totals["t_depth_affine"]))
    if totals["full_depth_affine"] is not None:
        rows.append(("Full depth (affine)", totals["full_depth_affine"]))
    for label, value in rows:
        print(f"    {label:<20} {format_base_2(value):<16} {value}")
    print(f"    {'Width':<20} {totals['width']}")
    print(f"    #1qClifford including the X gates that load the public plaintexts and tweaks:")
    print(f"    {'':<20} {format_base_2(oneq_with_loading):<16} {oneq_with_loading}")


def run_resource_report(v):
    from corefunctions import allocate_cells, load_hex
    from test_advanced_counter import count_resources, print_resources

    def build_encryption(eng):
        state = allocate_cells(eng, 16)
        tweaks = [allocate_cells(eng, 16) for _ in range(v.num_tweaks)]
        key = allocate_cells(eng, 32)
        v.encrypt(eng, state, tweaks, key, v.allocate_work(eng))

    def build_diffusion(eng):
        key = eng.allocate_qureg(K_BITS)
        ancillas, helpers = v.diffusion_workspace
        grover_diffusion(eng, key, eng.allocate_qureg(ancillas), eng.allocate_qureg(helpers))

    def build_iterations(m):
        def build(eng):
            lanes = _allocate_lanes(eng, v)
            for _ in range(m):
                _run_iteration(eng, v, lanes)
        return build

    def build_setup_key_only(eng):
        for q in eng.allocate_qureg(K_BITS):
            H | q

    def build_setup_with_loading(eng):
        key = eng.allocate_qureg(K_BITS)
        lanes = {
            "states": [allocate_cells(eng, 16) for _ in range(R_TUPLES)],
            "tweaks": [[allocate_cells(eng, 16) for _ in range(v.num_tweaks)] for _ in range(R_TUPLES)],
        }
        for q in key:
            H | q
        _load_public(v, lanes)

    enc = count_resources(build_encryption)
    oracle = count_resources(lambda eng: _run_oracle(eng, v, _allocate_lanes(eng, v)))
    diffusion = count_resources(build_diffusion)
    composed = [count_resources(build_iterations(m)) for m in (1, 2, 3)]
    iteration = composed[0]
    setup = count_resources(build_setup_key_only)
    setup_loading = count_resources(build_setup_with_loading)

    print("=" * 70)
    print(f"{v.title}: Grover key-search resources")
    print("=" * 70)
    print_resources("[1] Reversible encryption (paper Table 4)", enc)
    print_resources("[2] Three-tuple phase oracle U_f (paper Table 2)", oracle)
    print_resources("[3] Diffusion D_256 (paper Table 3)", diffusion)
    print_resources("[4] One Grover iteration G = D_256 U_f (paper Table 5)", iteration)
    _check_reference("encryption", enc, v.reference["encryption"])
    _check_reference("phase oracle", oracle, v.reference["oracle"])
    _check_reference("diffusion", diffusion, v.reference["diffusion"])
    _check_reference("iteration", iteration, v.reference["iteration"])
    print("    [1]-[4] agree with the reference values in the script")

    print("[5] Composition of m consecutive iterations")
    print(f"    {'m':<3} {'#CNOT':>9} {'#1qClifford':>12} {'#T':>9} {'#Measure':>9} {'T-depth':>8} {'Full depth':>11}")
    for m, r in enumerate(composed, 1):
        print(f"    {m:<3} {r['cnot']:>9} {r['oneq']:>12} {r['t']:>9} {r['measure']:>9} {r['t_depth']:>8} {r['full_depth']:>11}")
    for key in ("cnot", "oneq", "t", "measure"):
        assert composed[1][key] == 2 * iteration[key] and composed[2][key] == 3 * iteration[key], f"{key} not additive"
    steps = {}
    for key, label in (("t_depth", "T-depth"), ("full_depth", "full depth")):
        d12 = composed[1][key] - composed[0][key]
        d23 = composed[2][key] - composed[1][key]
        steps[key] = d12 if d12 == d23 else None
        print(f"    {label} step: {d12}, {d23}" + ("" if d12 == d23 else " (not constant, affine estimate omitted)"))
    print("    gate counts are additive in m")

    N_G = get_grover_iterations()
    p_fail = 1 - grover_success_probability(N_G)
    print("[6] Search parameters")
    print(f"    k = {K_BITS}, n = {N_BITS}, r = {R_TUPLES}, compared bits = {COMPARED_BITS}")
    print(f"    N_G = floor(pi/4 * 2^128) = {N_G}")
    print(f"    expected spurious keys = 2^(k - rn) = {format_base_2(expected_spurious_keys(K_BITS, R_TUPLES, N_BITS))}")
    print(f"    Grover failure probability = {mp.nstr(p_fail, 3)}")

    totals = calculate_exhaustive_totals(setup, iteration, readout_measure=K_BITS, readout_full_depth=1,
                                         N_G=N_G, delta_t=steps["t_depth"], delta_full=steps["full_depth"])
    totals_loading = calculate_exhaustive_totals(setup_loading, iteration, readout_measure=K_BITS,
                                                 readout_full_depth=1, N_G=N_G)
    print("[7] Complete key search (paper Table 7)")
    print("    setup H^256 on the key, N_G iterations, readout of the 256 key qubits;")
    print("    depths are N_G times one iteration, affine values use the step from [5]")
    _print_totals(totals, totals_loading["oneq"])


def grover_main(v, argv):
    """Command line of the *_full_grover.py scripts.

    --resource    resource report (default)
    --all-tests   functional tests followed by the resource report
    """
    mode = argv[0] if argv else "--resource"
    if mode == "--resource":
        run_resource_report(v)
    elif mode == "--all-tests":
        run_functional_tests(v)
        run_resource_report(v)
        print(f"{v.title}: all tests passed")
    else:
        print("usage: [--resource | --all-tests]")


# =====================================================================
# CLI Entry Point for Common Module Correctness Testing
# =====================================================================

if __name__ == "__main__":
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == "--test-diffusion":
            test_small_diffusion()
        elif arg == "--test-toy-grover":
            test_toy_complete_grover()
        else:
            print("Usage: python3 dialga_grover_common.py [--test-diffusion | --test-toy-grover]")
    else:
        test_small_diffusion()
        test_toy_complete_grover()
