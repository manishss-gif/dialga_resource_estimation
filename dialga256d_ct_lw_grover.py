"""
Phase oracle and inverse encryption for Dialga-256† (256d) CT-LW.
"""

from projectq import MainEngine
from projectq.backends import ClassicalSimulator
from projectq.ops import X, Measure, All, Z, CNOT
from qdialga256d_ct_lw import dialga256d_encrypt_ct_low_width
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
    round_function_invCT
)
from dialga_grover_common import build_and_tree, uncompute_and_tree


# === 1. Reversible Inverse Encryption ===
def dialga256d_encrypt_ct_low_width_inverse(eng, state, tweak0, tweak1, key, ancilla, helper):
    """Performs the exact inverse of dialga256d_encrypt_ct_low_width in-place."""
    K0 = key[:16]
    K1 = key[16:]
    
    # === 5. Post-whitening (reversed) ===
    add_key_128(eng, state, K1)
    add_key_128(eng, state, K0)
    tweak1_perm = permute_cells_middle(tweak1)
    add_tweak_128(eng, state, tweak1_perm)
    for cell_idx in range(16):
        sbox_8bitCT(eng, state[cell_idx], cell_idx, ancilla, helper)
        
    # === 4. Backward Rounds (Rb: 3 iterations, reversed i = beta down to 1) ===
    alpha = 4
    beta = 3
    for i in range(beta, 0, -1):
        if i == beta:
            # Reverse update tweak1
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
            
            # Reverse round 2(alpha+i)+1
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            state[:] = round_function_invCT(eng, state, 2 * (alpha + i) + 1, ancilla, helper)
            
            # Reverse round 2(alpha+i)
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            state[:] = round_function_invCT(eng, state, 2 * (alpha + i), ancilla, helper)
            
            # Reverse update tweak0
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
        else:
            # Reverse update tweak1
            tweak1[:] = round_functionCT(eng, tweak1, (alpha - i - 1) % 4, ancilla, helper)
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
            
            # Reverse round 2(alpha+i)+1
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            state[:] = round_function_invCT(eng, state, 2 * (alpha + i) + 1, ancilla, helper)
            
            # Reverse round 2(alpha+i)
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            state[:] = round_function_invCT(eng, state, 2 * (alpha + i), ancilla, helper)
            
            # Reverse update tweak0
            tweak0[:] = round_functionCT(eng, tweak0, (alpha - i - 1) % 4, ancilla, helper)
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
            
    # === 3. Middle Rounds (Rm: 1 iteration, reversed) ===
    # Reverse update tweak1
    tweak1[:] = round_functionCT(eng, tweak1, (alpha - 1) % 4, ancilla, helper)
    K_tweak_1 = K0 if (alpha % 2 == 0) else K1
    add_key_128(eng, tweak1, K_tweak_1)
    
    # Reverse round 2alpha+1
    add_round_constant(eng, state, cm_constants[1])
    tweak0_perm = permute_cells_middle(tweak0)
    add_tweak_128(eng, state, tweak0_perm)
    state[:] = round_function_invCT(eng, state, 2 * alpha + 1, ancilla, helper)
    
    # Reverse update tweak0
    tweak0[:] = round_functionCT(eng, tweak0, (alpha - 1) % 4, ancilla, helper)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak0, K_tweak)
    for cell_idx in range(16):
        sbox_8bitCT(eng, tweak0[cell_idx], cell_idx, ancilla, helper)
        
    # Reverse round 2alpha
    add_round_constant(eng, state, cm_constants[0])
    add_tweak_128(eng, state, tweak0)
    add_tweak_128(eng, state, tweak1)
    state[:] = round_function_invCT(eng, state, 2 * alpha, ancilla, helper)
    
    # Reverse tweak0 Sbox prep
    for cell_idx in range(16):
        sbox_8bitCT(eng, tweak0[cell_idx], cell_idx, ancilla, helper)
        
    # === 2. Forward Rounds (Rf: 4 iterations, reversed i = alpha down to 1) ===
    for i in range(alpha, 0, -1):
        if i == 1:
            # Reverse update tweak1
            add_key_128(eng, tweak1, K1)
            
            # Reverse round 2i-1 (round 1)
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            add_tweak_128(eng, state, tweak0)
            state[:] = round_function_invCT(eng, state, 2 * i - 1, ancilla, helper)
            
            # Reverse round 2i-2 (round 0)
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            add_key_128(eng, state, K1)
            add_tweak_128(eng, state, tweak1)
            state[:] = round_function_invCT(eng, state, 2 * i - 2, ancilla, helper)
            
            # Reverse update tweak0
            add_key_128(eng, tweak0, K0)
        else:
            # Reverse update tweak1
            K_tweak_1 = K1 if (i % 2 == 1) else K0
            add_key_128(eng, tweak1, K_tweak_1)
            tweak1[:] = round_function_invCT(eng, tweak1, (i - 1) % 4, ancilla, helper)
            
            # Reverse round 2i-1
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            add_tweak_128(eng, state, tweak0)
            state[:] = round_function_invCT(eng, state, 2 * i - 1, ancilla, helper)
            
            # Reverse round 2i-2
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            add_tweak_128(eng, state, tweak1)
            state[:] = round_function_invCT(eng, state, 2 * i - 2, ancilla, helper)
            
            # Reverse update tweak0
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak0, K_tweak)
            tweak0[:] = round_function_invCT(eng, tweak0, (i - 1) % 4, ancilla, helper)
            
    # === 1. Pre-whitening (reversed) ===
    add_key_128(eng, state, K1)
    add_key_128(eng, state, K0)
    add_tweak_128(eng, state, tweak0)


# === 2. Phase Oracle ===
def dialga256d_grover_oracle_ct_lw(
    eng,
    key0,
    key1,
    key2,
    states,
    tweaks0,
    tweaks1,
    enc_ancillas,
    enc_helpers,
    cmp_ancillas,
    cmp_helpers,
    expected_ciphertexts,
    test_flag=None
):
    assert len(expected_ciphertexts) == 3, "expected_ciphertexts must have exactly 3 elements"
    for c_hex in expected_ciphertexts:
        assert len(c_hex) == 32, f"Each expected ciphertext must be a 32-character hex string, got {len(c_hex)} chars"

    # Key copy fanout
    flat_key0 = [q for r in key0 for q in r]
    flat_key1 = [q for r in key1 for q in r]
    flat_key2 = [q for r in key2 for q in r]
    
    for j in range(256):
        CNOT | (flat_key0[j], flat_key1[j])
        CNOT | (flat_key1[j], flat_key2[j])
        
    # Parallel Encryptions
    dialga256d_encrypt_ct_low_width(
        eng, states[0], tweaks0[0], tweaks1[0], key0, enc_ancillas[0], enc_helpers[0]
    )
    dialga256d_encrypt_ct_low_width(
        eng, states[1], tweaks0[1], tweaks1[1], key1, enc_ancillas[1], enc_helpers[1]
    )
    dialga256d_encrypt_ct_low_width(
        eng, states[2], tweaks0[2], tweaks1[2], key2, enc_ancillas[2], enc_helpers[2]
    )
    
    # Flatten candidate bits
    candidate_bits = [q for state_reg in states for r in state_reg for q in r]
    
    # Expected bits parsing
    expected_bits = []
    for c_hex in expected_ciphertexts:
        bin_str = bin(int(c_hex, 16))[2:].zfill(128)
        for cell_idx in range(16):
            for bit_idx in range(8):
                expected_bits.append(int(bin_str[8 * cell_idx + bit_idx]))
                
    assert len(candidate_bits) == 384, f"Expected 384 candidate bits, got {len(candidate_bits)}"
    assert len(expected_bits) == 384, f"Expected 384 expected bits, got {len(expected_bits)}"
    assert len(candidate_bits) == len(expected_bits)
    
    assert len(cmp_ancillas) >= 383, f"Expected at least 383 comparator ancillas, got {len(cmp_ancillas)}"
    assert len(cmp_helpers) >= 192, f"Expected at least 192 comparator helpers, got {len(cmp_helpers)}"

    # Transform equality to all-ones
    for c_bit, ref_bit in zip(candidate_bits, expected_bits):
        if ref_bit == 0:
            X | c_bit
            
    # Compute AND tree
    match, history = build_and_tree(eng, candidate_bits, cmp_ancillas, cmp_helpers)
    assert len(history) == 383, f"Expected AND tree history size to be 383, got {len(history)}"
    
    # Phase flip
    if test_flag is not None:
        CNOT | (match, test_flag)
    else:
        Z | match
        
    # Uncompute AND tree and mask
    uncompute_and_tree(history)
    for c_bit, ref_bit in zip(candidate_bits, expected_bits):
        if ref_bit == 0:
            X | c_bit
            
    # Parallel Decryptions
    dialga256d_encrypt_ct_low_width_inverse(
        eng, states[0], tweaks0[0], tweaks1[0], key0, enc_ancillas[0], enc_helpers[0]
    )
    dialga256d_encrypt_ct_low_width_inverse(
        eng, states[1], tweaks0[1], tweaks1[1], key1, enc_ancillas[1], enc_helpers[1]
    )
    dialga256d_encrypt_ct_low_width_inverse(
        eng, states[2], tweaks0[2], tweaks1[2], key2, enc_ancillas[2], enc_helpers[2]
    )
    
    # Unfanout keys
    for j in reversed(range(256)):
        CNOT | (flat_key1[j], flat_key2[j])
        CNOT | (flat_key0[j], flat_key1[j])


# === Helper for loading hex strings ===
def load_hex_to_state(hex_str, state_reg):
    bin_str = bin(int(hex_str, 16))[2:].zfill(len(hex_str) * 4)
    for cell_idx in range(len(state_reg)):
        for bit_idx in range(8):
            val = int(bin_str[8 * cell_idx + bit_idx])
            if val == 1:
                X | state_reg[cell_idx][bit_idx]


def test_encrypt_inverse_round_trip():
    print("Testing Dialga-256† Encryption/Inverse Round-Trip...")
    sim = ClassicalSimulator()
    eng = MainEngine(sim, [])
    
    state = [eng.allocate_qureg(8) for _ in range(16)]
    tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
    tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
    key = [eng.allocate_qureg(8) for _ in range(32)]
    ancilla = eng.allocate_qubit()[0]
    helper = eng.allocate_qubit()[0]
    
    p_hex = "00112233445566778899aabbccddeeff"
    t0_hex = "2233445566778899aabbccddeeff1122"
    t1_hex = "33445566778899aabbccddeeff001100"
    k_hex = "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00"
    
    load_hex_to_state(p_hex, state)
    load_hex_to_state(t0_hex, tweak0)
    load_hex_to_state(t1_hex, tweak1)
    load_hex_to_state(k_hex, key)
    eng.flush()
    
    # Forward
    dialga256d_encrypt_ct_low_width(eng, state, tweak0, tweak1, key, ancilla, helper)
    eng.flush()
    
    # Inverse
    dialga256d_encrypt_ct_low_width_inverse(eng, state, tweak0, tweak1, key, ancilla, helper)
    eng.flush()
    
    def read_state_as_hex(state_reg):
        bin_str = ""
        for cell_idx in range(len(state_reg)):
            for bit_idx in range(8):
                bin_str += str(sim.read_bit(state_reg[cell_idx][bit_idx]))
        return hex(int(bin_str, 2))[2:].zfill(len(state_reg) * 2)
        
    actual_p = read_state_as_hex(state)
    assert actual_p == p_hex, f"Round-trip failed: expected {p_hex}, got {actual_p}"
    
    # Cleanup
    for r in state + tweak0 + tweak1 + key:
        All(Measure) | r
    Measure | ancilla
    Measure | helper
    eng.flush()
    print("  Round-Trip test: SUCCESS\n")


if __name__ == "__main__":
    test_encrypt_inverse_round_trip()
