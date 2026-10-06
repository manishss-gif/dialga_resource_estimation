"""
Dialga-128† encryption (16 rounds), low-width Clifford+T profile (CT-LW).

One product ancilla and one helper qubit are shared by every 4-bit
S-box, so the S-boxes of a SubCell layer are evaluated one after another.
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


def dialga128d_encrypt_ct_low_width(eng, state, tweak, key, ancilla, helper):
    """Performs one complete Dialga-128d (16-round) forward encryption in-place on the state using Clifford+T low-width sequential Sbox."""
    K0 = key[:16]
    K1 = key[16:]
    
    # === 1. Pre-whitening ===
    add_tweak_128(eng, state, tweak)
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    # === 2. Forward Rounds (Rf: 4 iterations) ===
    alpha = 4
    for i in range(1, alpha + 1):
        # Round 2i-2 on state
        state[:] = round_functionCT(eng, state, 2 * i - 2, ancilla, helper)
        # XOR key and constant
        K_round = K1 if (i % 2 == 1) else K0
        add_key_128(eng, state, K_round)
        add_round_constant(eng, state, cf_constants[2 * (i - 1)])
        
        # Update tweak state St
        if i == 1:
            add_key_128(eng, tweak, K0)
        else:
            tweak[:] = round_functionCT(eng, tweak, (i - 1) % 4, ancilla, helper)
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak, K_tweak)
            
        # Round 2i-1 on state
        state[:] = round_functionCT(eng, state, 2 * i - 1, ancilla, helper)
        # XOR tweak and constant (NO separate key addition here!)
        add_tweak_128(eng, state, tweak)
        add_round_constant(eng, state, cf_constants[2 * i - 1])
        
    # === 3. Middle Rounds (Rm: 1 iteration) ===
    # Round 2alpha (round 8) on state
    state[:] = round_functionCT(eng, state, 2 * alpha, ancilla, helper)
    # XOR key, constant, and SubCell(St)
    K_rm = K1 if (alpha % 2 == 0) else K0
    add_key_128(eng, state, K_rm)
    add_round_constant(eng, state, cm_constants[0])
    
    # XOR SubCell(St) to state in-place:
    for cell_idx in range(16):
        sbox_8bitCT(eng, tweak[cell_idx], cell_idx, ancilla, helper)
    add_tweak_128(eng, state, tweak)
    for cell_idx in range(16):
        sbox_8bitCT(eng, tweak[cell_idx], cell_idx, ancilla, helper)  # Restore tweak
        
    # Update tweak state St
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak, K_tweak)
    tweak[:] = round_function_invCT(eng, tweak, (alpha - 1) % 4, ancilla, helper)
    
    # Round 2alpha+1 (round 9) on state
    state[:] = round_functionCT(eng, state, 2 * alpha + 1, ancilla, helper)
    # XOR MS(St) and constant
    tweak_perm = permute_cells_middle(tweak)
    add_tweak_128(eng, state, tweak_perm)
    add_round_constant(eng, state, cm_constants[1])
    
    # === 4. Backward Rounds (Rb: 3 iterations) ===
    beta = 3
    for i in range(1, beta + 1):
        # Round 2(alpha+i) on state (rounds 10, 12, 14)
        state[:] = round_functionCT(eng, state, 2 * (alpha + i), ancilla, helper)
        # XOR key and constant
        K_round = K0 if ((alpha - i - 1) % 2 == 0) else K1
        add_key_128(eng, state, K_round)
        add_round_constant(eng, state, cb_constants[2 * (i - 1)])
        
        # Update tweak state St
        K_tweak = K0 if ((alpha - i - 1) % 2 == 0) else K1
        add_key_128(eng, tweak, K_tweak)
        if i < beta:
            tweak[:] = round_function_invCT(eng, tweak, (alpha - i - 1) % 4, ancilla, helper)
            
        # Round 2(alpha+i)+1 on state (rounds 11, 13, 15)
        state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, ancilla, helper)
        # XOR MS(St) and constant
        tweak_perm = permute_cells_middle(tweak)
        add_tweak_128(eng, state, tweak_perm)
        add_round_constant(eng, state, cb_constants[2 * i - 1])
        
    # === 5. Post-whitening ===
    # Final SubCell on state
    for cell_idx in range(16):
        sbox_8bitCT(eng, state[cell_idx], cell_idx, ancilla, helper)
    # XOR K0 ^ K1
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    return state


def allocate_work(eng):
    """Product ancilla and helper qubit shared by all S-boxes."""
    ancilla = eng.allocate_qubit()[0]
    helper = eng.allocate_qubit()[0]
    return [ancilla, helper]


def encrypt(eng, state, tweaks, key, work):
    dialga128d_encrypt_ct_low_width(eng, state, *tweaks, key, *work)


# Official test vector from the Dialga specification.
KAT = {
    "plaintext": "00112233445566778899aabbccddeeff",
    "tweaks": (
        "2233445566778899aabbccddee00ff11",
    ),
    "key": "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00",
    "ciphertext": "838407143af9a876fbdc6be378e9045b",
}


if __name__ == "__main__":
    encryption_main("Dialga-128† CT-LW", allocate_work, encrypt, KAT, sys.argv[1:])
