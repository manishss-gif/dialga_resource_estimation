"""
Phase oracle and reversible encryption inverse for Dialga-256d CT-LD profile.
"""

from projectq import MainEngine
from projectq.backends import ClassicalSimulator
from projectq.ops import X, Measure, All, Z, CNOT

from qdialga256d_ct_ld import dialga256d_encrypt_ct_low_depth
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


# === 1. Reversible Inverse of CT-LD Quantum Encryption ===
def dialga256d_encrypt_ct_low_depth_inverse(eng, state, tweak0, tweak1, key, state_ancillae, state_helpers, tweak_ancillae, tweak_helpers):
    """Performs the exact inverse of dialga256d_encrypt_ct_low_depth in-place."""
    K0 = key[:16]
    K1 = key[16:]
    
    # Reverse of === 5. Post-whitening ===
    add_key_128(eng, state, K1)
    add_key_128(eng, state, K0)
    tweak1_perm = permute_cells_middle(tweak1)
    add_tweak_128(eng, state, tweak1_perm)
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, state[cell_idx], cell_idx,
            state_ancillae[2 * cell_idx], state_helpers[2 * cell_idx],
            state_ancillae[2 * cell_idx + 1], state_helpers[2 * cell_idx + 1]
        )
        
    # Reverse of === 4. Backward Rounds (Rb: 3 iterations) ===
    alpha = 4
    beta = 3
    for i in range(beta, 0, -1):
        if i == beta:
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
            
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            state[:] = round_function_invCT(eng, state, 2 * (alpha + i) + 1, state_ancillae, state_helpers)
            
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            state[:] = round_function_invCT(eng, state, 2 * (alpha + i), state_ancillae, state_helpers)
            
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
        else:
            tweak1[:] = round_functionCT(eng, tweak1, (alpha - i - 1) % 4, tweak_ancillae, tweak_helpers)
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
            
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            state[:] = round_function_invCT(eng, state, 2 * (alpha + i) + 1, state_ancillae, state_helpers)
            
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            state[:] = round_function_invCT(eng, state, 2 * (alpha + i), state_ancillae, state_helpers)
            
            tweak0[:] = round_functionCT(eng, tweak0, (alpha - i - 1) % 4, tweak_ancillae, tweak_helpers)
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
            
    # Reverse of === 3. Middle Rounds (Rm: 1 iteration) ===
    tweak1[:] = round_functionCT(eng, tweak1, (alpha - 1) % 4, tweak_ancillae, tweak_helpers)
    K_tweak_1 = K0 if (alpha % 2 == 0) else K1
    add_key_128(eng, tweak1, K_tweak_1)
    
    add_round_constant(eng, state, cm_constants[1])
    tweak0_perm = permute_cells_middle(tweak0)
    add_tweak_128(eng, state, tweak0_perm)
    state[:] = round_function_invCT(eng, state, 2 * alpha + 1, state_ancillae, state_helpers)
    
    tweak0[:] = round_functionCT(eng, tweak0, (alpha - 1) % 4, tweak_ancillae, tweak_helpers)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak0, K_tweak)
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, tweak0[cell_idx], cell_idx,
            tweak_ancillae[2 * cell_idx], tweak_helpers[2 * cell_idx],
            tweak_ancillae[2 * cell_idx + 1], tweak_helpers[2 * cell_idx + 1]
        )
        
    add_round_constant(eng, state, cm_constants[0])
    add_tweak_128(eng, state, tweak0)
    add_tweak_128(eng, state, tweak1)
    state[:] = round_function_invCT(eng, state, 2 * alpha, state_ancillae, state_helpers)
    
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, tweak0[cell_idx], cell_idx,
            tweak_ancillae[2 * cell_idx], tweak_helpers[2 * cell_idx],
            tweak_ancillae[2 * cell_idx + 1], tweak_helpers[2 * cell_idx + 1]
        )
        
    # Reverse of === 2. Forward Rounds (Rf: 4 iterations) ===
    for i in range(alpha, 0, -1):
        if i == 1:
            add_key_128(eng, tweak1, K1)
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            add_tweak_128(eng, state, tweak0)
            state[:] = round_function_invCT(eng, state, 2 * i - 1, state_ancillae, state_helpers)
            
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            add_key_128(eng, state, K1)
            add_tweak_128(eng, state, tweak1)
            state[:] = round_function_invCT(eng, state, 2 * i - 2, state_ancillae, state_helpers)
            add_key_128(eng, tweak0, K0)
        else:
            K_tweak_1 = K1 if (i % 2 == 1) else K0
            add_key_128(eng, tweak1, K_tweak_1)
            tweak1[:] = round_function_invCT(eng, tweak1, (i - 1) % 4, tweak_ancillae, tweak_helpers)
            
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            add_tweak_128(eng, state, tweak0)
            state[:] = round_function_invCT(eng, state, 2 * i - 1, state_ancillae, state_helpers)
            
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            add_tweak_128(eng, state, tweak1)
            state[:] = round_function_invCT(eng, state, 2 * i - 2, state_ancillae, state_helpers)
            
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak0, K_tweak)
            tweak0[:] = round_function_invCT(eng, tweak0, (i - 1) % 4, tweak_ancillae, tweak_helpers)
            
    # Reverse of === 1. Pre-whitening ===
    add_key_128(eng, state, K1)
    add_key_128(eng, state, K0)
    add_tweak_128(eng, state, tweak0)


def load_hex_to_state(hex_str, state_reg):
    bin_str = bin(int(hex_str, 16))[2:].zfill(len(hex_str) * 4)
    for cell_idx in range(len(state_reg)):
        for bit_idx in range(8):
            val = int(bin_str[8 * cell_idx + bit_idx])
            if val == 1:
                X | state_reg[cell_idx][bit_idx]


def apply_equality_mask(eng, candidate_bits, expected_bits):
    assert len(candidate_bits) == len(expected_bits), "Comparator input length mismatch"
    for q, val in zip(candidate_bits, expected_bits):
        if val == 0:
            X | q


# === 2. Dialga-256d CT-LD Phase Oracle ===
def dialga256d_grover_oracle_ct_ld(
    eng,
    key0, key1, key2,
    states, tweaks0, tweaks1,
    state_ancillae_list, state_helpers_list,
    tweak_ancillae_list, tweak_helpers_list,
    cmp_ancillas, cmp_helpers,
    expected_ciphertexts,
    test_flag=None
):
    """Executes the complete Dialga-256d CT-LD Grover phase oracle."""
    assert len(expected_ciphertexts) == 3, "expected_ciphertexts must have 3 hex strings"
    for c_hex in expected_ciphertexts:
        assert len(c_hex) == 32, "Each ciphertext hex string must be 32 chars"
    assert len(states) == 3 and len(tweaks0) == 3 and len(tweaks1) == 3, "Phase oracle requires 3 parallel tuples"
    assert len(state_ancillae_list) == 3 and len(state_helpers_list) == 3
    assert len(tweak_ancillae_list) == 3 and len(tweak_helpers_list) == 3
    assert len(cmp_ancillas) >= 383, "Comparator requires at least 383 ancillas"
    assert len(cmp_helpers) >= 192, "Comparator requires at least 192 helpers"
    
    flat_k0 = [q for cell in key0 for q in cell]
    flat_k1 = [q for cell in key1 for q in cell]
    flat_k2 = [q for cell in key2 for q in cell]
    
    # 1. Coherent Key Fanout (512 CNOTs)
    for j in range(256):
        CNOT | (flat_k0[j], flat_k1[j])
        CNOT | (flat_k1[j], flat_k2[j])
        
    # 2. Parallel CT-LD Forward Encryptions (3 disjoint branches)
    dialga256d_encrypt_ct_low_depth(
        eng, states[0], tweaks0[0], tweaks1[0], key0,
        state_ancillae_list[0], state_helpers_list[0],
        tweak_ancillae_list[0], tweak_helpers_list[0]
    )
    dialga256d_encrypt_ct_low_depth(
        eng, states[1], tweaks0[1], tweaks1[1], key1,
        state_ancillae_list[1], state_helpers_list[1],
        tweak_ancillae_list[1], tweak_helpers_list[1]
    )
    dialga256d_encrypt_ct_low_depth(
        eng, states[2], tweaks0[2], tweaks1[2], key2,
        state_ancillae_list[2], state_helpers_list[2],
        tweak_ancillae_list[2], tweak_helpers_list[2]
    )
    
    # 3. Concatenate 384 Candidate Ciphertext Qubits
    candidate_bits = []
    for s in states:
        for cell in s:
            candidate_bits.extend(cell)
    assert len(candidate_bits) == 384
    
    expected_bits = []
    for c_hex in expected_ciphertexts:
        bin_str = bin(int(c_hex, 16))[2:].zfill(128)
        expected_bits.extend([int(b) for b in bin_str])
    assert len(expected_bits) == 384
    assert len(candidate_bits) == len(expected_bits)
    
    # 4. Equality Masking
    apply_equality_mask(eng, candidate_bits, expected_bits)
    
    # 5. Balanced 384-Input AND Tree
    match_qubit, history = build_and_tree(eng, candidate_bits, cmp_ancillas, cmp_helpers)
    assert len(history) == 383, f"Expected 383 AND nodes, got {len(history)}"
    
    # 6. Phase Mark
    if test_flag is not None:
        CNOT | (match_qubit, test_flag)
    else:
        Z | match_qubit
        
    # 7. Uncompute AND Tree & Equality Mask
    uncompute_and_tree(history)
    apply_equality_mask(eng, candidate_bits, expected_bits)
    
    # 8. Parallel CT-LD ENC† Uncomputation (3 disjoint branches)
    dialga256d_encrypt_ct_low_depth_inverse(
        eng, states[0], tweaks0[0], tweaks1[0], key0,
        state_ancillae_list[0], state_helpers_list[0],
        tweak_ancillae_list[0], tweak_helpers_list[0]
    )
    dialga256d_encrypt_ct_low_depth_inverse(
        eng, states[1], tweaks0[1], tweaks1[1], key1,
        state_ancillae_list[1], state_helpers_list[1],
        tweak_ancillae_list[1], tweak_helpers_list[1]
    )
    dialga256d_encrypt_ct_low_depth_inverse(
        eng, states[2], tweaks0[2], tweaks1[2], key2,
        state_ancillae_list[2], state_helpers_list[2],
        tweak_ancillae_list[2], tweak_helpers_list[2]
    )
    
    # 9. Reverse Key Fanout (512 CNOTs)
    for j in reversed(range(256)):
        CNOT | (flat_k1[j], flat_k2[j])
        CNOT | (flat_k0[j], flat_k1[j])


def test_encrypt_inverse_round_trip():
    """Validates that CT-LD inverse encryption exactly uncomputes forward encryption."""
    print("Testing Dialga-256d CT-LD Encryption/Inverse Round-Trip...")
    sim = ClassicalSimulator()
    eng = MainEngine(sim, [])
    
    state = [eng.allocate_qureg(8) for _ in range(16)]
    tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
    tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
    key = [eng.allocate_qureg(8) for _ in range(32)]
    state_ancillae = eng.allocate_qureg(32)
    state_helpers = eng.allocate_qureg(32)
    tweak_ancillae = eng.allocate_qureg(32)
    tweak_helpers = eng.allocate_qureg(32)
    
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
    dialga256d_encrypt_ct_low_depth(eng, state, tweak0, tweak1, key, state_ancillae, state_helpers, tweak_ancillae, tweak_helpers)
    eng.flush()
    
    # Inverse
    dialga256d_encrypt_ct_low_depth_inverse(eng, state, tweak0, tweak1, key, state_ancillae, state_helpers, tweak_ancillae, tweak_helpers)
    eng.flush()
    
    def read_state_as_hex(state_reg):
        bin_str = ""
        for cell_idx in range(len(state_reg)):
            for bit_idx in range(8):
                bin_str += str(sim.read_bit(state_reg[cell_idx][bit_idx]))
        return hex(int(bin_str, 2))[2:].zfill(len(state_reg) * 2)
        
    actual_p = read_state_as_hex(state)
    actual_t0 = read_state_as_hex(tweak0)
    actual_t1 = read_state_as_hex(tweak1)
    actual_k = read_state_as_hex(key)
    assert actual_p == p_hex, f"Round-trip failed for state: expected {p_hex}, got {actual_p}"
    assert actual_t0 == t0_hex, f"Round-trip failed for tweak0: expected {t0_hex}, got {actual_t0}"
    assert actual_t1 == t1_hex, f"Round-trip failed for tweak1: expected {t1_hex}, got {actual_t1}"
    assert actual_k == k_hex, f"Round-trip failed for key: expected {k_hex}, got {actual_k}"
    
    for q in state_ancillae + state_helpers + tweak_ancillae + tweak_helpers:
        assert sim.read_bit(q) == 0, "CT-LD work qubit dirty after round trip"
        
    for r in state + tweak0 + tweak1 + key:
        All(Measure) | r
    All(Measure) | state_ancillae
    All(Measure) | state_helpers
    All(Measure) | tweak_ancillae
    All(Measure) | tweak_helpers
    eng.flush()
    print("  Round-Trip test: SUCCESS\n")


if __name__ == "__main__":
    test_encrypt_inverse_round_trip()
