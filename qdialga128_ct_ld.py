"""
Dialga-128 encryption (20 rounds), low-depth Clifford+T profile (CT-LD).

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


def dialga128_encrypt_ct_low_depth(eng, state, tweak, key, state_ancillae, state_helpers, tweak_ancillae, tweak_helpers):
    """Performs one complete Dialga-128 forward encryption in-place on the state using Clifford+T Sbox with parallel Sbox layers and separate tweak helper registers."""
    K0 = key[:16]
    K1 = key[16:]
    
    # === 1. Pre-whitening ===
    add_tweak_128(eng, state, tweak)
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    # === 2. Forward Rounds (Rf: 5 iterations) ===
    alpha = 5
    for i in range(1, alpha + 1):
        # Round 2i-2 on state
        state[:] = round_functionCT(eng, state, 2 * i - 2, state_ancillae, state_helpers)
        # XOR key and constant
        K_round = K1 if (i % 2 == 1) else K0
        add_key_128(eng, state, K_round)
        add_round_constant(eng, state, cf_constants[2 * (i - 1)])
        
        # Update tweak state St (runs in parallel with state Round 2i-2 since registers are independent)
        if i == 1:
            add_key_128(eng, tweak, K0)
        else:
            tweak[:] = round_functionCT(eng, tweak, (i - 1) % 4, tweak_ancillae, tweak_helpers)
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak, K_tweak)
            
        # Round 2i-1 on state
        state[:] = round_functionCT(eng, state, 2 * i - 1, state_ancillae, state_helpers)
        # XOR tweak and constant (NO separate key addition here!)
        add_tweak_128(eng, state, tweak)
        add_round_constant(eng, state, cf_constants[2 * i - 1])
        
    # === 3. Middle Rounds (Rm: 1 iteration) ===
    # Round 2alpha (round 10) on state
    state[:] = round_functionCT(eng, state, 2 * alpha, state_ancillae, state_helpers)
    # XOR key, constant, and SubCell(St)
    K_rm = K1 if (alpha % 2 == 0) else K0
    add_key_128(eng, state, K_rm)
    add_round_constant(eng, state, cm_constants[0])
    
    # XOR SubCell(St) to state in-place:
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, tweak[cell_idx], cell_idx,
            tweak_ancillae[2 * cell_idx], tweak_helpers[2 * cell_idx],
            tweak_ancillae[2 * cell_idx + 1], tweak_helpers[2 * cell_idx + 1]
        )
    add_tweak_128(eng, state, tweak)
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, tweak[cell_idx], cell_idx,
            tweak_ancillae[2 * cell_idx], tweak_helpers[2 * cell_idx],
            tweak_ancillae[2 * cell_idx + 1], tweak_helpers[2 * cell_idx + 1]
        )  # Restore tweak
        
    # Update tweak state St
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak, K_tweak)
    tweak[:] = round_function_invCT(eng, tweak, (alpha - 1) % 4, tweak_ancillae, tweak_helpers)
    
    # Round 2alpha+1 (round 11) on state
    state[:] = round_functionCT(eng, state, 2 * alpha + 1, state_ancillae, state_helpers)
    # XOR MS(St) and constant
    tweak_perm = permute_cells_middle(tweak)
    add_tweak_128(eng, state, tweak_perm)
    add_round_constant(eng, state, cm_constants[1])
    
    # === 4. Backward Rounds (Rb: 4 iterations) ===
    beta = 4
    for i in range(1, beta + 1):
        # Round 2(alpha+i) on state (rounds 12, 14, 16, 18)
        state[:] = round_functionCT(eng, state, 2 * (alpha + i), state_ancillae, state_helpers)
        # XOR key and constant
        K_round = K0 if ((alpha - i - 1) % 2 == 0) else K1
        add_key_128(eng, state, K_round)
        add_round_constant(eng, state, cb_constants[2 * (i - 1)])
        
        # Update tweak state St (runs in parallel with state Round 2(alpha+i) since registers are independent)
        K_tweak = K0 if ((alpha - i - 1) % 2 == 0) else K1
        add_key_128(eng, tweak, K_tweak)
        if i < beta:
            tweak[:] = round_function_invCT(eng, tweak, (alpha - i - 1) % 4, tweak_ancillae, tweak_helpers)
            
        # Round 2(alpha+i)+1 on state (rounds 13, 15, 17, 19)
        state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, state_ancillae, state_helpers)
        # XOR MS(St) and constant
        tweak_perm = permute_cells_middle(tweak)
        add_tweak_128(eng, state, tweak_perm)
        add_round_constant(eng, state, cb_constants[2 * i - 1])
        
    # === 5. Post-whitening ===
    # Final SubCell on state
    for cell_idx in range(16):
        sbox_8bitCT(
            eng, state[cell_idx], cell_idx,
            state_ancillae[2 * cell_idx], state_helpers[2 * cell_idx],
            state_ancillae[2 * cell_idx + 1], state_helpers[2 * cell_idx + 1]
        )
    # XOR K0 ^ K1
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
    dialga128_encrypt_ct_low_depth(eng, state, *tweaks, key, *work)


# Official test vector from the Dialga specification.
KAT = {
    "plaintext": "00112233445566778899aabbccddeeff",
    "tweaks": (
        "2233445566778899aabbccddee00ff11",
    ),
    "key": "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00",
    "ciphertext": "a4a1ea948919d8996e13b1b365bb0ce6",
}


if __name__ == "__main__":
    encryption_main("Dialga-128 CT-LD", allocate_work, encrypt, KAT, sys.argv[1:])
