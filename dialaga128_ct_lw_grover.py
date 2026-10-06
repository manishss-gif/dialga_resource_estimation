from projectq import MainEngine
from projectq.backends import ClassicalSimulator
from projectq.ops import X, Measure, All, Z, CNOT

from qdialga128_ct_lw import dialga128_encrypt_ct_low_width
from corefunctions import (
    cf_constants,
    cm_constants,
    cb_constants,
    add_tweak_128,
    add_key_128,
    add_round_constant,
    permute_cells_middle,
    sbox_8bitCT,
    round_functionCT,
    round_function_invCT,
    and_gate,
    and_gate_dagger
)

# === 1. Explicit Inverse of Quantum Encryption ===
def dialga128_encrypt_ct_low_width_inverse(eng, state, tweak, key, ancilla, helper):
    """Performs the exact inverse of dialga128_encrypt_ct_low_width in-place."""
    K0 = key[:16]
    K1 = key[16:]
    
    # Reverse of === 5. Post-whitening ===
    add_key_128(eng, state, K1)
    add_key_128(eng, state, K0)
    for cell_idx in range(16):
        sbox_8bitCT(eng, state[cell_idx], cell_idx, ancilla, helper)
        
    # Reverse of === 4. Backward Rounds (Rb: 4 iterations) ===
    alpha = 5
    beta = 4
    for i in range(beta, 0, -1):
        # Reverse second state round update inside the iteration
        add_round_constant(eng, state, cb_constants[2 * i - 1])
        tweak_perm = permute_cells_middle(tweak)
        add_tweak_128(eng, state, tweak_perm)
        state[:] = round_function_invCT(eng, state, 2 * (alpha + i) + 1, ancilla, helper)
        
        # Reverse tweak state update
        if i < beta:
            tweak[:] = round_functionCT(eng, tweak, (alpha - i - 1) % 4, ancilla, helper)
        K_tweak = K0 if ((alpha - i - 1) % 2 == 0) else K1
        add_key_128(eng, tweak, K_tweak)
        
        # Reverse first state round update inside the iteration
        add_round_constant(eng, state, cb_constants[2 * (i - 1)])
        K_round = K0 if ((alpha - i - 1) % 2 == 0) else K1
        add_key_128(eng, state, K_round)
        state[:] = round_function_invCT(eng, state, 2 * (alpha + i), ancilla, helper)
        
    # Reverse of === 3. Middle Rounds (Rm: 1 iteration) ===
    add_round_constant(eng, state, cm_constants[1])
    tweak_perm = permute_cells_middle(tweak)
    add_tweak_128(eng, state, tweak_perm)
    state[:] = round_function_invCT(eng, state, 2 * alpha + 1, ancilla, helper)
    
    tweak[:] = round_functionCT(eng, tweak, (alpha - 1) % 4, ancilla, helper)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak, K_tweak)
    
    for cell_idx in range(16):
        sbox_8bitCT(eng, tweak[cell_idx], cell_idx, ancilla, helper)
    add_tweak_128(eng, state, tweak)
    for cell_idx in range(16):
        sbox_8bitCT(eng, tweak[cell_idx], cell_idx, ancilla, helper)
        
    add_round_constant(eng, state, cm_constants[0])
    K_rm = K1 if (alpha % 2 == 0) else K0
    add_key_128(eng, state, K_rm)
    state[:] = round_function_invCT(eng, state, 2 * alpha, ancilla, helper)
    
    # Reverse of === 2. Forward Rounds (Rf: 5 iterations) ===
    for i in range(alpha, 0, -1):
        # Reverse second state round update
        add_round_constant(eng, state, cf_constants[2 * i - 1])
        add_tweak_128(eng, state, tweak)
        state[:] = round_function_invCT(eng, state, 2 * i - 1, ancilla, helper)
        
        # Reverse tweak state update
        if i == 1:
            add_key_128(eng, tweak, K0)
        else:
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak, K_tweak)
            tweak[:] = round_function_invCT(eng, tweak, (i - 1) % 4, ancilla, helper)
            
        # Reverse first state round update
        add_round_constant(eng, state, cf_constants[2 * (i - 1)])
        K_round = K1 if (i % 2 == 1) else K0
        add_key_128(eng, state, K_round)
        state[:] = round_function_invCT(eng, state, 2 * i - 2, ancilla, helper)
        
    # Reverse of === 1. Pre-whitening ===
    add_key_128(eng, state, K1)
    add_key_128(eng, state, K0)
    add_tweak_128(eng, state, tweak)

# Helper to load hex string into qubits
def load_hex_to_state(hex_str, state_reg):
    bin_str = bin(int(hex_str, 16))[2:].zfill(len(hex_str) * 4)
    for cell_idx in range(len(state_reg)):
        for bit_idx in range(8):
            val = int(bin_str[8 * cell_idx + bit_idx])
            if val == 1:
                X | state_reg[cell_idx][bit_idx]

# === 2. Comparator Equality Mask ===
def apply_equality_mask(eng, candidate_bits, expected_bits):
    assert len(candidate_bits) == len(expected_bits), "Comparator input length mismatch"
    for q, val in zip(candidate_bits, expected_bits):
        if val == 0:
            X | q

# === 3. Balanced AND Tree and Uncomputation ===
def build_and_tree(eng, current_level, cmp_ancillas, cmp_helpers):
    num_inputs = len(current_level)
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

# === 4. Actual Grover Phase Oracle ===
def dialga128_grover_oracle_ct_lw(
    eng,
    key0,
    key1,
    key2,
    states,
    tweaks,
    enc_ancillas,
    enc_helpers,
    cmp_ancillas,
    cmp_helpers,
    expected_ciphertexts,
    test_flag=None
):
    assert len(expected_ciphertexts) == 3, "Expected exactly 3 target ciphertexts"
    for c_hex in expected_ciphertexts:
        assert len(c_hex) == 32, "Each expected ciphertext must be exactly 32 hex chars"

    # Coherent Key Fanout
    flat_key0 = []
    flat_key1 = []
    flat_key2 = []
    for i in range(32):
        for j in range(8):
            flat_key0.append(key0[i][j])
            flat_key1.append(key1[i][j])
            flat_key2.append(key2[i][j])
            
    for j in range(256):
        CNOT | (flat_key0[j], flat_key1[j])
        CNOT | (flat_key1[j], flat_key2[j])
        
    # Parallel Encryptions on 3 branches
    dialga128_encrypt_ct_low_width(
        eng, states[0], tweaks[0], key0, enc_ancillas[0], enc_helpers[0]
    )
    dialga128_encrypt_ct_low_width(
        eng, states[1], tweaks[1], key1, enc_ancillas[1], enc_helpers[1]
    )
    dialga128_encrypt_ct_low_width(
        eng, states[2], tweaks[2], key2, enc_ancillas[2], enc_helpers[2]
    )
    
    # Flatten outputs for comparator (3 * 128 = 384 bits)
    candidate_bits = []
    for state_reg in states:
        for cell_idx in range(16):
            for bit_idx in range(8):
                candidate_bits.append(state_reg[cell_idx][bit_idx])
                
    expected_bits = []
    for c_hex in expected_ciphertexts:
        bin_str = bin(int(c_hex, 16))[2:].zfill(128)
        for cell_idx in range(16):
            for bit_idx in range(8):
                expected_bits.append(int(bin_str[8 * cell_idx + bit_idx]))
                
    assert len(candidate_bits) == 384, f"Expected 384 candidate ciphertext bits, got {len(candidate_bits)}"
    assert len(expected_bits) == 384, f"Expected 384 reference ciphertext bits, got {len(expected_bits)}"
    assert len(candidate_bits) == len(expected_bits)

    # Workspace assertions for comparator
    assert len(cmp_ancillas) >= len(candidate_bits) - 1, f"Expected at least {len(candidate_bits) - 1} cmp_ancillas"
    assert len(cmp_helpers) >= len(candidate_bits) // 2, f"Expected at least {len(candidate_bits) // 2} cmp_helpers"

    # Transform equality to all-ones
    apply_equality_mask(eng, candidate_bits, expected_bits)
    
    # Compute Balanced AND tree
    match, history = build_and_tree(eng, candidate_bits, cmp_ancillas, cmp_helpers)
    assert len(history) == 383, f"Expected 383 AND nodes, got {len(history)}"
    
    # Phase mark Z (or CNOT in test mode)
    if test_flag is not None:
        CNOT | (match, test_flag)
    else:
        Z | match
        
    # Uncompute comparator
    uncompute_and_tree(history)
    apply_equality_mask(eng, candidate_bits, expected_bits)
    
    # Parallel Inverse Encryptions
    dialga128_encrypt_ct_low_width_inverse(
        eng, states[0], tweaks[0], key0, enc_ancillas[0], enc_helpers[0]
    )
    dialga128_encrypt_ct_low_width_inverse(
        eng, states[1], tweaks[1], key1, enc_ancillas[1], enc_helpers[1]
    )
    dialga128_encrypt_ct_low_width_inverse(
        eng, states[2], tweaks[2], key2, enc_ancillas[2], enc_helpers[2]
    )
    
    # Reverse Key Fanout
    for j in reversed(range(256)):
        CNOT | (flat_key1[j], flat_key2[j])
        CNOT | (flat_key0[j], flat_key1[j])

# === 5. Unit and Correctness Tests ===
def test_encrypt_inverse_round_trip():
    print("Running ENC -> ENC_INV Round-Trip Tests...")
    
    # Helper to read state as hex
    def read_state_as_hex(sim, state_reg):
        bin_str = ""
        for cell_idx in range(len(state_reg)):
            for bit_idx in range(8):
                bin_str += str(sim.read_bit(state_reg[cell_idx][bit_idx]))
        return hex(int(bin_str, 2))[2:].zfill(len(state_reg) * 2)

    # 1. Test official KAT roundtrip
    sim = ClassicalSimulator()
    eng = MainEngine(sim, [])
    
    p_hex = "00112233445566778899aabbccddeeff"
    t_hex = "2233445566778899aabbccddee00ff11"
    k_hex = "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00"
    
    key = [eng.allocate_qureg(8) for _ in range(32)]
    tweak = [eng.allocate_qureg(8) for _ in range(16)]
    state = [eng.allocate_qureg(8) for _ in range(16)]
    ancilla = eng.allocate_qubit()[0]
    helper = eng.allocate_qubit()[0]
    
    load_hex_to_state(p_hex, state)
    load_hex_to_state(t_hex, tweak)
    load_hex_to_state(k_hex, key)
    eng.flush()
    
    dialga128_encrypt_ct_low_width(eng, state, tweak, key, ancilla, helper)
    eng.flush()
    
    dialga128_encrypt_ct_low_width_inverse(eng, state, tweak, key, ancilla, helper)
    eng.flush()
    
    assert read_state_as_hex(sim, state) == p_hex, "KAT state not restored"
    assert read_state_as_hex(sim, tweak) == t_hex, "KAT tweak not restored"
    assert read_state_as_hex(sim, key) == k_hex, "KAT key not restored"
    assert sim.read_bit(ancilla) == 0, "KAT ancilla not clean"
    assert sim.read_bit(helper) == 0, "KAT helper not clean"
    
    for reg in state + tweak + key:
        All(Measure) | reg
    Measure | ancilla
    Measure | helper
    eng.flush()
    print("  Official KAT validation: SUCCESS")
    
    # 2. Test 20 deterministic random values
    import random
    rng = random.Random(20260807)
    sim = ClassicalSimulator()
    eng = MainEngine(sim, [])
    
    key = [eng.allocate_qureg(8) for _ in range(32)]
    tweak = [eng.allocate_qureg(8) for _ in range(16)]
    state = [eng.allocate_qureg(8) for _ in range(16)]
    ancilla = eng.allocate_qubit()[0]
    helper = eng.allocate_qubit()[0]
    
    for idx in range(20):
        p = "".join(rng.choice("0123456789abcdef") for _ in range(32))
        t = "".join(rng.choice("0123456789abcdef") for _ in range(32))
        k = "".join(rng.choice("0123456789abcdef") for _ in range(64))
        
        load_hex_to_state(p, state)
        load_hex_to_state(t, tweak)
        load_hex_to_state(k, key)
        eng.flush()
        
        dialga128_encrypt_ct_low_width(eng, state, tweak, key, ancilla, helper)
        eng.flush()
        
        dialga128_encrypt_ct_low_width_inverse(eng, state, tweak, key, ancilla, helper)
        eng.flush()
        
        assert read_state_as_hex(sim, state) == p, f"Random {idx} state mismatch"
        assert read_state_as_hex(sim, tweak) == t, f"Random {idx} tweak mismatch"
        assert read_state_as_hex(sim, key) == k, f"Random {idx} key mismatch"
        assert sim.read_bit(ancilla) == 0, f"Random {idx} ancilla mismatch"
        assert sim.read_bit(helper) == 0, f"Random {idx} helper mismatch"
        
        # Reset registers to 0 by un-applying load_hex_to_state
        load_hex_to_state(p, state)
        load_hex_to_state(t, tweak)
        load_hex_to_state(k, key)
        eng.flush()
        
    for reg in state + tweak + key:
        All(Measure) | reg
    Measure | ancilla
    Measure | helper
    eng.flush()
    print("  20 Random round-trip validation: SUCCESS")
    print("All ENC -> ENC_INV tests passed!\n")


if __name__ == "__main__":
    test_encrypt_inverse_round_trip()
