"""
Dialga-256 encryption (20 rounds), low-depth Clifford+T profile (CT-LD).

The state and the tweak each get 32 product ancillas and 32 helpers,
so the 32 S-boxes of a SubCell layer are evaluated in parallel.
"""

import sys

from corefunctions import (
    add_key_128,
    add_round_constant,
    add_tweak_128,
    cb_constants,
    cf_constants,
    cm_constants,
    permute_cells_middle,
    round_functionCT,
    round_function_invCT,
    sbox_8bitCT,
    encryption_main,
)


def dialga256_encrypt_ct_low_depth(eng, state, tweak0, tweak1, key, state_ancillae, state_helpers, tweak_ancillae, tweak_helpers):
    """Performs the complete Dialga-256 encryption in-place on the state using Clifford+T low-depth parallel Sbox and separate tweak helper registers."""
    K0 = key[:16]
    K1 = key[16:]
    
    # === 1. Pre-whitening ===
    add_tweak_128(eng, state, tweak0)
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    # === 2. Forward Rounds (Rf: 5 iterations) ===
    alpha = 5
    for i in range(1, alpha + 1):
        if i == 1:
            # Update tweak0
            add_key_128(eng, tweak0, K0)
            # Round 2i-2 (round 0) on state
            state[:] = round_functionCT(eng, state, 2 * i - 2, state_ancillae, state_helpers)
            add_tweak_128(eng, state, tweak1)
            add_key_128(eng, state, K1)
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            
            # Round 2i-1 (round 1) on state
            state[:] = round_functionCT(eng, state, 2 * i - 1, state_ancillae, state_helpers)
            add_tweak_128(eng, state, tweak0)
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            # Update tweak1
            add_key_128(eng, tweak1, K1)
        else:
            # Update tweak0 (runs in parallel with state Round 2i-2 since registers are independent)
            tweak0[:] = round_functionCT(eng, tweak0, (i - 1) % 4, tweak_ancillae, tweak_helpers)
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak0, K_tweak)
            
            # Round 2i-2 on state
            state[:] = round_functionCT(eng, state, 2 * i - 2, state_ancillae, state_helpers)
            add_tweak_128(eng, state, tweak1)
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            
            # Round 2i-1 on state
            state[:] = round_functionCT(eng, state, 2 * i - 1, state_ancillae, state_helpers)
            add_tweak_128(eng, state, tweak0)
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            
            # Update tweak1
            tweak1[:] = round_functionCT(eng, tweak1, (i - 1) % 4, tweak_ancillae, tweak_helpers)
            K_tweak_1 = K1 if (i % 2 == 1) else K0
            add_key_128(eng, tweak1, K_tweak_1)
            
    # === 3. Middle Rounds (Rm: 1 iteration) ===
    # Tweak0 Sbox prep
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, tweak0[cell_idx], cell_idx,
            tweak_ancillae[2 * cell_idx], tweak_helpers[2 * cell_idx],
            tweak_ancillae[2 * cell_idx + 1], tweak_helpers[2 * cell_idx + 1]
        )
        
    # Round 2alpha (round 10) on state
    state[:] = round_functionCT(eng, state, 2 * alpha, state_ancillae, state_helpers)
    add_tweak_128(eng, state, tweak1)
    add_tweak_128(eng, state, tweak0)
    add_round_constant(eng, state, cm_constants[0])
    
    # Update tweak0 (undo S-box prep & add key & invert round)
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, tweak0[cell_idx], cell_idx,
            tweak_ancillae[2 * cell_idx], tweak_helpers[2 * cell_idx],
            tweak_ancillae[2 * cell_idx + 1], tweak_helpers[2 * cell_idx + 1]
        )  # Restore tweak0 (SubCell^-1)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak0, K_tweak)
    tweak0[:] = round_function_invCT(eng, tweak0, (alpha - 1) % 4, tweak_ancillae, tweak_helpers)
    
    # Round 2alpha+1 (round 11) on state
    state[:] = round_functionCT(eng, state, 2 * alpha + 1, state_ancillae, state_helpers)
    tweak0_perm = permute_cells_middle(tweak0)
    add_tweak_128(eng, state, tweak0_perm)
    add_round_constant(eng, state, cm_constants[1])
    
    # Update tweak1
    K_tweak_1 = K0 if (alpha % 2 == 0) else K1
    add_key_128(eng, tweak1, K_tweak_1)
    tweak1[:] = round_function_invCT(eng, tweak1, (alpha - 1) % 4, tweak_ancillae, tweak_helpers)
    
    # === 4. Backward Rounds (Rb: 4 iterations) ===
    beta = 4
    for i in range(1, beta + 1):
        if i == beta:
            # Update tweak0
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
            
            # Round 2(alpha+i) on state
            state[:] = round_functionCT(eng, state, 2 * (alpha + i), state_ancillae, state_helpers)
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            
            # Round 2(alpha+i)+1 on state
            state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, state_ancillae, state_helpers)
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            
            # Update tweak1
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
        else:
            # Update tweak0 (runs in parallel with state Round 2(alpha+i) since registers are independent)
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
            tweak0[:] = round_function_invCT(eng, tweak0, (alpha - i - 1) % 4, tweak_ancillae, tweak_helpers)
            
            # Round 2(alpha+i) on state
            state[:] = round_functionCT(eng, state, 2 * (alpha + i), state_ancillae, state_helpers)
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            
            # Round 2(alpha+i)+1 on state
            state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, state_ancillae, state_helpers)
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            
            # Update tweak1
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
            tweak1[:] = round_function_invCT(eng, tweak1, (alpha - i - 1) % 4, tweak_ancillae, tweak_helpers)
            
    # === 5. Post-whitening ===
    # Final SubCell on state
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, state[cell_idx], cell_idx,
            state_ancillae[2 * cell_idx], state_helpers[2 * cell_idx],
            state_ancillae[2 * cell_idx + 1], state_helpers[2 * cell_idx + 1]
        )
    # XOR MS(St1) ^ K0 ^ K1
    tweak1_perm = permute_cells_middle(tweak1)
    add_tweak_128(eng, state, tweak1_perm)
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    return state


def allocate_work(eng):
    """Product ancillas and helpers for the state and the tweak SubCell layers."""
    state_ancillae = eng.allocate_qureg(32)
    state_helpers = eng.allocate_qureg(32)
    tweak_ancillae = eng.allocate_qureg(32)
    tweak_helpers = eng.allocate_qureg(32)
    return [state_ancillae, state_helpers, tweak_ancillae, tweak_helpers]


def encrypt(eng, state, tweaks, key, work):
    dialga256_encrypt_ct_low_depth(eng, state, *tweaks, key, *work)


# Official test vector from the Dialga specification.
KAT = {
    "plaintext": "00112233445566778899aabbccddeeff",
    "tweaks": (
        "2233445566778899aabbccddeeff1122",
        "33445566778899aabbccddeeff001100",
    ),
    "key": "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00",
    "ciphertext": "e129ff0920371753db65532540d06881",
}


if __name__ == "__main__":
    encryption_main("Dialga-256 CT-LD", allocate_work, encrypt, KAT, sys.argv[1:])
