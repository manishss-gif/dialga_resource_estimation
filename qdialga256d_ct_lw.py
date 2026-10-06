"""
Dialga-256† encryption (16 rounds), low-width Clifford+T profile (CT-LW).

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


def dialga256d_encrypt_ct_low_width(eng, state, tweak0, tweak1, key, ancilla, helper):
    """Performs the complete Dialga-256d (16-round Dialga-256†) encryption in-place on the state using Clifford+T low-width sequential Sbox."""
    K0 = key[:16]
    K1 = key[16:]
    
    # === 1. Pre-whitening ===
    add_tweak_128(eng, state, tweak0)
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    # === 2. Forward Rounds (Rf: 4 iterations) ===
    alpha = 4
    for i in range(1, alpha + 1):
        if i == 1:
            # Update tweak0
            add_key_128(eng, tweak0, K0)
            # Round 2i-2 (round 0) on state
            state[:] = round_functionCT(eng, state, 2 * i - 2, ancilla, helper)
            add_tweak_128(eng, state, tweak1)
            add_key_128(eng, state, K1)
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            
            # Round 2i-1 (round 1) on state
            state[:] = round_functionCT(eng, state, 2 * i - 1, ancilla, helper)
            add_tweak_128(eng, state, tweak0)
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            # Update tweak1
            add_key_128(eng, tweak1, K1)
        else:
            # Update tweak0
            tweak0[:] = round_functionCT(eng, tweak0, (i - 1) % 4, ancilla, helper)
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak0, K_tweak)
            
            # Round 2i-2 on state
            state[:] = round_functionCT(eng, state, 2 * i - 2, ancilla, helper)
            add_tweak_128(eng, state, tweak1)
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            
            # Round 2i-1 on state
            state[:] = round_functionCT(eng, state, 2 * i - 1, ancilla, helper)
            add_tweak_128(eng, state, tweak0)
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            
            # Update tweak1
            tweak1[:] = round_functionCT(eng, tweak1, (i - 1) % 4, ancilla, helper)
            K_tweak_1 = K1 if (i % 2 == 1) else K0
            add_key_128(eng, tweak1, K_tweak_1)
            
    # === 3. Middle Rounds (Rm: 1 iteration) ===
    # Tweak0 Sbox prep
    for cell_idx in range(16):
        sbox_8bitCT(eng, tweak0[cell_idx], cell_idx, ancilla, helper)
        
    # Round 2alpha (round 8) on state
    state[:] = round_functionCT(eng, state, 2 * alpha, ancilla, helper)
    add_tweak_128(eng, state, tweak1)
    add_tweak_128(eng, state, tweak0)
    add_round_constant(eng, state, cm_constants[0])
    
    # Update tweak0 (undo S-box prep & add key & invert round)
    for cell_idx in range(16):
        sbox_8bitCT(eng, tweak0[cell_idx], cell_idx, ancilla, helper)  # Restore tweak0 (SubCell^-1)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak0, K_tweak)
    tweak0[:] = round_function_invCT(eng, tweak0, (alpha - 1) % 4, ancilla, helper)
    
    # Round 2alpha+1 (round 9) on state
    state[:] = round_functionCT(eng, state, 2 * alpha + 1, ancilla, helper)
    tweak0_perm = permute_cells_middle(tweak0)
    add_tweak_128(eng, state, tweak0_perm)
    add_round_constant(eng, state, cm_constants[1])
    
    # Update tweak1
    K_tweak_1 = K0 if (alpha % 2 == 0) else K1
    add_key_128(eng, tweak1, K_tweak_1)
    tweak1[:] = round_function_invCT(eng, tweak1, (alpha - 1) % 4, ancilla, helper)
    
    # === 4. Backward Rounds (Rb: 3 iterations) ===
    beta = 3
    for i in range(1, beta + 1):
        if i == beta:
            # Update tweak0
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
            
            # Round 2(alpha+i) on state (round 10, 12, 14)
            state[:] = round_functionCT(eng, state, 2 * (alpha + i), ancilla, helper)
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            
            # Round 2(alpha+i)+1 on state (round 11, 13, 15)
            state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, ancilla, helper)
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            
            # Update tweak1
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
        else:
            # Update tweak0
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
            tweak0[:] = round_function_invCT(eng, tweak0, (alpha - i - 1) % 4, ancilla, helper)
            
            # Round 2(alpha+i) on state
            state[:] = round_functionCT(eng, state, 2 * (alpha + i), ancilla, helper)
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            
            # Round 2(alpha+i)+1 on state
            state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, ancilla, helper)
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            
            # Update tweak1
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
            tweak1[:] = round_function_invCT(eng, tweak1, (alpha - i - 1) % 4, ancilla, helper)
            
    # === 5. Post-whitening ===
    # Final SubCell on state
    for cell_idx in range(16):
        sbox_8bitCT(eng, state[cell_idx], cell_idx, ancilla, helper)
    # XOR MS(St1) ^ K0 ^ K1
    tweak1_perm = permute_cells_middle(tweak1)
    add_tweak_128(eng, state, tweak1_perm)
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    return state


def allocate_work(eng):
    """Product ancilla and helper qubit shared by all S-boxes."""
    ancilla = eng.allocate_qubit()[0]
    helper = eng.allocate_qubit()[0]
    return [ancilla, helper]


def encrypt(eng, state, tweaks, key, work):
    dialga256d_encrypt_ct_low_width(eng, state, *tweaks, key, *work)


# Official test vector from the Dialga specification.
KAT = {
    "plaintext": "00112233445566778899aabbccddeeff",
    "tweaks": (
        "2233445566778899aabbccddeeff1122",
        "33445566778899aabbccddeeff001100",
    ),
    "key": "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00",
    "ciphertext": "833367ce5fd49a7b1f0dd9ec62720018",
}


if __name__ == "__main__":
    encryption_main("Dialga-256† CT-LW", allocate_work, encrypt, KAT, sys.argv[1:])
