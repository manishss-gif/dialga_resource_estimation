# Independent Classical Reference Implementation of Dialga-256

# Bit permutation tables Pb,i
pb_tables = [
    [4, 1, 6, 3, 0, 5, 2, 7],  # Pb,0
    [1, 6, 7, 0, 5, 2, 3, 4],  # Pb,1
    [2, 3, 4, 1, 6, 7, 0, 5],  # Pb,2
    [7, 4, 1, 2, 3, 0, 5, 6]   # Pb,3
]

# Cell permutation tables pi_i
pi_tables = [
    [7, 0, 13, 10, 5, 2, 15, 8, 4, 3, 14, 9, 6, 1, 12, 11],  # pi_0
    [13, 0, 10, 7, 11, 6, 12, 1, 2, 15, 5, 8, 4, 9, 3, 14],  # pi_1
    [7, 13, 10, 0, 6, 12, 11, 1, 5, 15, 8, 2, 4, 14, 9, 3],  # pi_2
    [13, 8, 6, 3, 14, 11, 5, 0, 12, 9, 7, 2, 15, 10, 4, 1]   # pi_3
]

# Cell permutation table for the middle round pi_m
pi_m = [0, 10, 5, 15, 14, 4, 11, 1, 9, 3, 12, 6, 7, 13, 2, 8]

# Round constants (hexadecimal form)
cf_constants = [
    "243f6a8885a308d313198a2e03707344",  # cf0
    "a4093822299f31d0082efa98ec4e6c89",  # cf1
    "452821e638d01377be5466cf34e90c6c",  # cf2
    "c0ac29b7c97c50dd3f84d5b5b5470917",  # cf3
    "9216d5d98979fb1bd1310ba698dfb5ac",  # cf4
    "2ffd72dbd01adfb7b8e1afed6a267e96",  # cf5
    "ba7c9045f12c7f9924a19947b3916cf7",  # cf6
    "0801f2e2858efc16636920d871574e69",  # cf7
    "a458fea3f4933d7e0d95748f728eb658",  # cf8
    "718bcd5882154aee7b54a41dc25a59b5"   # cf9
]

cm_constants = [
    "9c30d5392af26013c5d1b023286085f0",  # cm0
    "ca417918b8db38ef8e79dcb0603a180e"   # cm1
]

cb_constants = [
    "6c9e0e8bb01e8a3ed71577c1bd314b27",  # cb0
    "78af2fda55605c60e65525f3aa55ab94",  # cb1
    "5748986263e8144055ca396a2aab10b6",  # cb2
    "b4cc5c341141e8cea15486af7c72e993",  # cb3
    "b3ee1411636fbc2a2ba9c55d741831f6",  # cb4
    "ce5c3e169b87931eafd6ba336c24cf5c",  # cb5
    "7a325381289586773b8f48986b4bb9af",  # cb6
    "c4bfe81b6628219361d809ccfb21a991"   # cb7
]

# Official 4-bit S-box Sb0
Sb0 = [0xc, 0xa, 0xd, 0x3, 0xe, 0xb, 0xf, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6]

def hex_to_bytes(hex_str):
    return [int(hex_str[2*i : 2*i+2], 16) for i in range(len(hex_str)//2)]

def bytes_to_hex(byte_list):
    return "".join(f"{x:02x}" for x in byte_list)

def validate_hex_inputs_256(plaintext_hex, tweak_hex, key_hex):
    """Audits input types, sizes, and hexadecimal characters for Dialga-256."""
    if not isinstance(plaintext_hex, str) or not isinstance(tweak_hex, str) or not isinstance(key_hex, str):
        raise ValueError("Inputs must be hexadecimal strings.")
        
    for name, hex_str in [("plaintext", plaintext_hex), ("tweak", tweak_hex), ("key", key_hex)]:
        if not all(c in "0123456789abcdefABCDEF" for c in hex_str):
            raise ValueError(f"Input '{name}' must contain only valid hexadecimal characters.")
            
    if len(plaintext_hex) != 32:
        raise ValueError(f"plaintext must be exactly 128 bits (32 hex characters), got {len(plaintext_hex)}")
    if len(tweak_hex) != 64:
        raise ValueError(f"tweak must be exactly 256 bits (64 hex characters), got {len(tweak_hex)}")
    if len(key_hex) != 64:
        raise ValueError(f"key must be exactly 256 bits (64 hex characters), got {len(key_hex)}")

def classical_sbox_8bit(val, cell_idx):
    """Classical mathematical reference of the 8-bit S-box SSb_j."""
    p = pb_tables[cell_idx % 4]
    in_bits = [(val >> (7 - i)) & 1 for i in range(8)]
    perm_bits = [in_bits[p[i]] for i in range(8)]
    hi = (perm_bits[0] << 3) | (perm_bits[1] << 2) | (perm_bits[2] << 1) | perm_bits[3]
    lo = (perm_bits[4] << 3) | (perm_bits[5] << 2) | (perm_bits[6] << 1) | perm_bits[7]
    hi_out = Sb0[hi]
    lo_out = Sb0[lo]
    out_bits_unperm = [
        (hi_out >> 3) & 1, (hi_out >> 2) & 1, (hi_out >> 1) & 1, hi_out & 1,
        (lo_out >> 3) & 1, (lo_out >> 2) & 1, (lo_out >> 1) & 1, lo_out & 1
    ]
    out_bits = [0] * 8
    for i in range(8):
        out_bits[p[i]] = out_bits_unperm[i]
    out_val = 0
    for i in range(8):
        out_val |= (out_bits[i] << (7 - i))
    return out_val

def invert_permutation(p):
    inv = [0] * len(p)
    for idx, val in enumerate(p):
        inv[val] = idx
    return inv

pi_tables_inv = [invert_permutation(p) for p in pi_tables]
pi_m_inv = invert_permutation(pi_m)

def permute_cells(state, round_idx):
    p = pi_tables[round_idx % 4]
    return [state[idx] for idx in p]

def permute_cells_inv(state, round_idx):
    return [state[idx] for idx in pi_tables_inv[round_idx % 4]]

def classical_matrix_mul(state):
    out = [0] * 16
    for col in range(4):
        c0 = state[4 * col]
        c1 = state[4 * col + 1]
        c2 = state[4 * col + 2]
        c3 = state[4 * col + 3]
        out[4 * col]     = c1 ^ c2 ^ c3
        out[4 * col + 1] = c0 ^ c2 ^ c3
        out[4 * col + 2] = c0 ^ c1 ^ c3
        out[4 * col + 3] = c0 ^ c1 ^ c2
    return out

def classical_round_function(state, round_idx):
    state = list(state)
    for cell_idx in range(16):
        state[cell_idx] = classical_sbox_8bit(state[cell_idx], cell_idx)
    state = permute_cells(state, round_idx)
    state = classical_matrix_mul(state)
    return state

def classical_round_function_inv(state, round_idx):
    state = classical_matrix_mul(state)
    state = permute_cells_inv(state, round_idx)
    for cell_idx in range(16):
        state[cell_idx] = classical_sbox_8bit(state[cell_idx], cell_idx)
    return state

# === Modular helper functions for Dialga-256 ===

def classical_rf_iteration_256(state, tweak0, tweak1, K0, K1, i):
    """Applies iteration i of Dialga-256 forward rounds Rf."""
    state = list(state)
    tweak0 = list(tweak0)
    tweak1 = list(tweak1)
    
    if i == 1:
        # Update tweak0
        for j in range(16):
            tweak0[j] ^= K0[j]
        # Round 2i-2 (round 0) on state
        state = classical_round_function(state, 2 * i - 2)
        cf_bytes = hex_to_bytes(cf_constants[2 * (i - 1)])
        for j in range(16):
            state[j] ^= tweak1[j] ^ K1[j] ^ cf_bytes[j]
            
        # Round 2i-1 (round 1) on state
        state = classical_round_function(state, 2 * i - 1)
        cf_bytes_1 = hex_to_bytes(cf_constants[2 * i - 1])
        for j in range(16):
            state[j] ^= tweak0[j] ^ cf_bytes_1[j]
        # Update tweak1
        for j in range(16):
            tweak1[j] ^= K1[j]
    else:
        # Update tweak0
        tweak0 = classical_round_function(tweak0, (i - 1) % 4)
        K_tweak = K0 if (i % 2 == 1) else K1
        for j in range(16):
            tweak0[j] ^= K_tweak[j]
            
        # Round 2i-2 on state
        state = classical_round_function(state, 2 * i - 2)
        cf_bytes = hex_to_bytes(cf_constants[2 * (i - 1)])
        for j in range(16):
            state[j] ^= tweak1[j] ^ cf_bytes[j]
            
        # Round 2i-1 on state
        state = classical_round_function(state, 2 * i - 1)
        cf_bytes_1 = hex_to_bytes(cf_constants[2 * i - 1])
        for j in range(16):
            state[j] ^= tweak0[j] ^ cf_bytes_1[j]
            
        # Update tweak1
        tweak1 = classical_round_function(tweak1, (i - 1) % 4)
        K_tweak_1 = K1 if (i % 2 == 1) else K0
        for j in range(16):
            tweak1[j] ^= K_tweak_1[j]
            
    return state, tweak0, tweak1

def classical_rm_256(state, tweak0, tweak1, K0, K1, alpha):
    """Applies the middle rounds Rm for Dialga-256."""
    state = list(state)
    tweak0 = list(tweak0)
    tweak1 = list(tweak1)
    
    # Tweak0 Sbox prep
    for j in range(16):
        tweak0[j] = classical_sbox_8bit(tweak0[j], j)
        
    # Round 2alpha on state
    state = classical_round_function(state, 2 * alpha)
    cm0_bytes = hex_to_bytes(cm_constants[0])
    for j in range(16):
        state[j] ^= tweak1[j] ^ tweak0[j] ^ cm0_bytes[j]
        
    # Restore tweak0 (undo S-box prep & add key & invert round)
    for j in range(16):
        tweak0[j] = classical_sbox_8bit(tweak0[j], j)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    for j in range(16):
        tweak0[j] ^= K_tweak[j]
    tweak0 = classical_round_function_inv(tweak0, (alpha - 1) % 4)
    
    # Round 2alpha+1 on state
    state = classical_round_function(state, 2 * alpha + 1)
    tweak0_perm = [tweak0[idx] for idx in pi_m]
    cm1_bytes = hex_to_bytes(cm_constants[1])
    for j in range(16):
        state[j] ^= tweak0_perm[j] ^ cm1_bytes[j]
        
    # Update tweak1
    K_tweak_1 = K0 if (alpha % 2 == 0) else K1
    for j in range(16):
        tweak1[j] ^= K_tweak_1[j]
    tweak1 = classical_round_function_inv(tweak1, (alpha - 1) % 4)
    
    return state, tweak0, tweak1

def classical_rb_iteration_256(state, tweak0, tweak1, K0, K1, alpha, i, beta):
    """Applies iteration i of the backward rounds Rb for Dialga-256."""
    state = list(state)
    tweak0 = list(tweak0)
    tweak1 = list(tweak1)
    
    if i == beta:
        # Update tweak0
        K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
        for j in range(16):
            tweak0[j] ^= K_tweak_0[j]
            
        # Round 2(alpha+i) on state
        state = classical_round_function(state, 2 * (alpha + i))
        tweak1_perm = [tweak1[idx] for idx in pi_m]
        cb_bytes = hex_to_bytes(cb_constants[2 * (i - 1)])
        for j in range(16):
            state[j] ^= tweak1_perm[j] ^ cb_bytes[j]
            
        # Round 2(alpha+i)+1 on state
        state = classical_round_function(state, 2 * (alpha + i) + 1)
        tweak0_perm = [tweak0[idx] for idx in pi_m]
        cb_bytes_1 = hex_to_bytes(cb_constants[2 * i - 1])
        for j in range(16):
            state[j] ^= tweak0_perm[j] ^ cb_bytes_1[j]
            
        # Update tweak1
        K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
        for j in range(16):
            tweak1[j] ^= K_tweak_1[j]
    else:
        # Update tweak0
        K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
        for j in range(16):
            tweak0[j] ^= K_tweak_0[j]
        tweak0 = classical_round_function_inv(tweak0, (alpha - i - 1) % 4)
        
        # Round 2(alpha+i) on state
        state = classical_round_function(state, 2 * (alpha + i))
        tweak1_perm = [tweak1[idx] for idx in pi_m]
        cb_bytes = hex_to_bytes(cb_constants[2 * (i - 1)])
        for j in range(16):
            state[j] ^= tweak1_perm[j] ^ cb_bytes[j]
            
        # Round 2(alpha+i)+1 on state
        state = classical_round_function(state, 2 * (alpha + i) + 1)
        tweak0_perm = [tweak0[idx] for idx in pi_m]
        cb_bytes_1 = hex_to_bytes(cb_constants[2 * i - 1])
        for j in range(16):
            state[j] ^= tweak0_perm[j] ^ cb_bytes_1[j]
            
        # Update tweak1
        K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
        for j in range(16):
            tweak1[j] ^= K_tweak_1[j]
        tweak1 = classical_round_function_inv(tweak1, (alpha - i - 1) % 4)
        
    return state, tweak0, tweak1

# === Complete Dialga-256 encryption/decryption ===

def classical_dialga256_encrypt(plaintext_hex, tweak_hex, key_hex, trace=False):
    """Complete classical reference encryption of Dialga-256."""
    validate_hex_inputs_256(plaintext_hex, tweak_hex, key_hex)
    
    state = hex_to_bytes(plaintext_hex)
    tweak = hex_to_bytes(tweak_hex)
    tweak0 = tweak[:16]
    tweak1 = tweak[16:]
    key = hex_to_bytes(key_hex)
    K0 = key[:16]
    K1 = key[16:]
    
    trace_dict = {}
    if trace:
        trace_dict["initial_state"] = bytes_to_hex(state)
        trace_dict["initial_tweak0"] = bytes_to_hex(tweak0)
        trace_dict["initial_tweak1"] = bytes_to_hex(tweak1)
        
    # === 1. Pre-whitening ===
    for j in range(16):
        state[j] ^= tweak0[j] ^ K0[j] ^ K1[j]
        
    # === 2. Forward Rounds (Rf) ===
    alpha = 5
    for i in range(1, alpha + 1):
        state, tweak0, tweak1 = classical_rf_iteration_256(state, tweak0, tweak1, K0, K1, i)
        if trace:
            trace_dict[f"rf_iteration_{i}_state"] = bytes_to_hex(state)
            trace_dict[f"rf_iteration_{i}_tweak0"] = bytes_to_hex(tweak0)
            trace_dict[f"rf_iteration_{i}_tweak1"] = bytes_to_hex(tweak1)
            
    # === 3. Middle Rounds (Rm) ===
    if trace:
        trace_dict["rm_state_before"] = bytes_to_hex(state)
        trace_dict["rm_tweak0_before"] = bytes_to_hex(tweak0)
        trace_dict["rm_tweak1_before"] = bytes_to_hex(tweak1)
        
    state, tweak0, tweak1 = classical_rm_256(state, tweak0, tweak1, K0, K1, alpha)
    if trace:
        trace_dict["rm_state_after"] = bytes_to_hex(state)
        trace_dict["rm_tweak0_after"] = bytes_to_hex(tweak0)
        trace_dict["rm_tweak1_after"] = bytes_to_hex(tweak1)
        
    # === 4. Backward Rounds (Rb) ===
    beta = 4
    for i in range(1, beta + 1):
        state, tweak0, tweak1 = classical_rb_iteration_256(state, tweak0, tweak1, K0, K1, alpha, i, beta)
        if trace:
            trace_dict[f"rb_iteration_{i}_state"] = bytes_to_hex(state)
            trace_dict[f"rb_iteration_{i}_tweak0"] = bytes_to_hex(tweak0)
            trace_dict[f"rb_iteration_{i}_tweak1"] = bytes_to_hex(tweak1)
            
    # === 5. Post-whitening ===
    for j in range(16):
        state[j] = classical_sbox_8bit(state[j], j)
    tweak1_perm = [tweak1[idx] for idx in pi_m]
    for j in range(16):
        state[j] ^= tweak1_perm[j] ^ K0[j] ^ K1[j]
        
    ciphertext = bytes_to_hex(state)
    if trace:
        trace_dict["final_ciphertext"] = ciphertext
        return ciphertext, trace_dict
    return ciphertext

def classical_dialga256_decrypt(ciphertext_hex, tweak_hex, key_hex):
    """Complete classical reference decryption of Dialga-256."""
    validate_hex_inputs_256(ciphertext_hex, tweak_hex, key_hex)
    
    state = hex_to_bytes(ciphertext_hex)
    tweak = hex_to_bytes(tweak_hex)
    tweak0 = tweak[:16]
    tweak1 = tweak[16:]
    key = hex_to_bytes(key_hex)
    K0 = key[:16]
    K1 = key[16:]
    
    # === 5. Inverse Post-whitening ===
    tweak1_perm = [tweak1[idx] for idx in pi_m]
    for j in range(16):
        state[j] ^= tweak1_perm[j] ^ K0[j] ^ K1[j]
    for j in range(16):
        state[j] = classical_sbox_8bit(state[j], j)
        
    # === 4. Inverse Backward Rounds (Rb) ===
    alpha = 5
    beta = 4
    for i in range(beta, 0, -1):
        if i == beta:
            # Undo tweak1 update
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            for j in range(16):
                tweak1[j] ^= K_tweak_1[j]
                
            # Undo Round 2(alpha+i)+1
            cb_bytes_1 = hex_to_bytes(cb_constants[2 * i - 1])
            tweak0_perm = [tweak0[idx] for idx in pi_m]
            for j in range(16):
                state[j] ^= tweak0_perm[j] ^ cb_bytes_1[j]
            state = classical_round_function_inv(state, 2 * (alpha + i) + 1)
            
            # Undo Round 2(alpha+i)
            cb_bytes = hex_to_bytes(cb_constants[2 * (i - 1)])
            tweak1_perm = [tweak1[idx] for idx in pi_m]
            for j in range(16):
                state[j] ^= tweak1_perm[j] ^ cb_bytes[j]
            state = classical_round_function_inv(state, 2 * (alpha + i))
            
            # Undo tweak0 update
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            for j in range(16):
                tweak0[j] ^= K_tweak_0[j]
        else:
            # Undo tweak1 update
            tweak1 = classical_round_function(tweak1, (alpha - i - 1) % 4)
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            for j in range(16):
                tweak1[j] ^= K_tweak_1[j]
                
            # Undo Round 2(alpha+i)+1
            cb_bytes_1 = hex_to_bytes(cb_constants[2 * i - 1])
            tweak0_perm = [tweak0[idx] for idx in pi_m]
            for j in range(16):
                state[j] ^= tweak0_perm[j] ^ cb_bytes_1[j]
            state = classical_round_function_inv(state, 2 * (alpha + i) + 1)
            
            # Undo Round 2(alpha+i)
            cb_bytes = hex_to_bytes(cb_constants[2 * (i - 1)])
            tweak1_perm = [tweak1[idx] for idx in pi_m]
            for j in range(16):
                state[j] ^= tweak1_perm[j] ^ cb_bytes[j]
            state = classical_round_function_inv(state, 2 * (alpha + i))
            
            # Undo tweak0 update
            tweak0 = classical_round_function(tweak0, (alpha - i - 1) % 4)
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            for j in range(16):
                tweak0[j] ^= K_tweak_0[j]
                
    # === 3. Inverse Middle Rounds (Rm) ===
    # Undo tweak1 update
    tweak1 = classical_round_function(tweak1, (alpha - 1) % 4)
    K_tweak_1 = K0 if (alpha % 2 == 0) else K1
    for j in range(16):
        tweak1[j] ^= K_tweak_1[j]
        
    # Undo Round 2alpha+1
    cm1_bytes = hex_to_bytes(cm_constants[1])
    tweak0_perm = [tweak0[idx] for idx in pi_m]
    for j in range(16):
        state[j] ^= tweak0_perm[j] ^ cm1_bytes[j]
    state = classical_round_function_inv(state, 2 * alpha + 1)
    
    # Undo tweak0 update
    tweak0 = classical_round_function(tweak0, (alpha - 1) % 4)
    K_tweak_0 = K0 if ((alpha - 1) % 2 == 0) else K1
    for j in range(16):
        tweak0[j] ^= K_tweak_0[j]
    for j in range(16):
        tweak0[j] = classical_sbox_8bit(tweak0[j], j)
        
    # Undo Round 2alpha
    cm0_bytes = hex_to_bytes(cm_constants[0])
    for j in range(16):
        state[j] ^= tweak1[j] ^ tweak0[j] ^ cm0_bytes[j]
    state = classical_round_function_inv(state, 2 * alpha)
    
    # Undo tweak0 Sbox prep
    for j in range(16):
        tweak0[j] = classical_sbox_8bit(tweak0[j], j)
        
    # === 2. Inverse Forward Rounds (Rf) ===
    for i in range(alpha, 0, -1):
        if i == 1:
            # Undo tweak1 update
            for j in range(16):
                tweak1[j] ^= K1[j]
                
            # Undo Round 2i-1
            cf_bytes_1 = hex_to_bytes(cf_constants[2 * i - 1])
            for j in range(16):
                state[j] ^= tweak0[j] ^ cf_bytes_1[j]
            state = classical_round_function_inv(state, 2 * i - 1)
            
            # Undo Round 2i-2
            cf_bytes = hex_to_bytes(cf_constants[2 * (i - 1)])
            for j in range(16):
                state[j] ^= tweak1[j] ^ K1[j] ^ cf_bytes[j]
            state = classical_round_function_inv(state, 2 * i - 2)
            
            # Undo tweak0 update
            for j in range(16):
                tweak0[j] ^= K0[j]
        else:
            # Undo tweak1 update
            K_tweak_1 = K1 if (i % 2 == 1) else K0
            for j in range(16):
                tweak1[j] ^= K_tweak_1[j]
            tweak1 = classical_round_function_inv(tweak1, (i - 1) % 4)
            
            # Undo Round 2i-1
            cf_bytes_1 = hex_to_bytes(cf_constants[2 * i - 1])
            for j in range(16):
                state[j] ^= tweak0[j] ^ cf_bytes_1[j]
            state = classical_round_function_inv(state, 2 * i - 1)
            
            # Undo Round 2i-2
            cf_bytes = hex_to_bytes(cf_constants[2 * (i - 1)])
            for j in range(16):
                state[j] ^= tweak1[j] ^ cf_bytes[j]
            state = classical_round_function_inv(state, 2 * i - 2)
            
            # Undo tweak0 update
            K_tweak_0 = K0 if (i % 2 == 1) else K1
            for j in range(16):
                tweak0[j] ^= K_tweak_0[j]
            tweak0 = classical_round_function_inv(tweak0, (i - 1) % 4)
            
    # === 1. Inverse Pre-whitening ===
    for j in range(16):
        state[j] ^= tweak0[j] ^ K0[j] ^ K1[j]
        
    return bytes_to_hex(state)

# === Dagger Variant (Dialga-256†) ===

def classical_dialga256d_encrypt(plaintext_hex, tweak_hex, key_hex, trace=False):
    """Complete classical reference encryption of Dialga-256 dagger variant."""
    validate_hex_inputs_256(plaintext_hex, tweak_hex, key_hex)
    
    state = hex_to_bytes(plaintext_hex)
    tweak = hex_to_bytes(tweak_hex)
    tweak0 = tweak[:16]
    tweak1 = tweak[16:]
    key = hex_to_bytes(key_hex)
    K0 = key[:16]
    K1 = key[16:]
    
    trace_dict = {}
    if trace:
        trace_dict["initial_state"] = bytes_to_hex(state)
        trace_dict["initial_tweak0"] = bytes_to_hex(tweak0)
        trace_dict["initial_tweak1"] = bytes_to_hex(tweak1)
        
    # === 1. Pre-whitening ===
    for j in range(16):
        state[j] ^= tweak0[j] ^ K0[j] ^ K1[j]
        
    # === 2. Forward Rounds (Rf) ===
    alpha = 4
    for i in range(1, alpha + 1):
        state, tweak0, tweak1 = classical_rf_iteration_256(state, tweak0, tweak1, K0, K1, i)
        if trace:
            trace_dict[f"rf_iteration_{i}_state"] = bytes_to_hex(state)
            trace_dict[f"rf_iteration_{i}_tweak0"] = bytes_to_hex(tweak0)
            trace_dict[f"rf_iteration_{i}_tweak1"] = bytes_to_hex(tweak1)
            
    # === 3. Middle Rounds (Rm) ===
    if trace:
        trace_dict["rm_state_before"] = bytes_to_hex(state)
        trace_dict["rm_tweak0_before"] = bytes_to_hex(tweak0)
        trace_dict["rm_tweak1_before"] = bytes_to_hex(tweak1)
        
    state, tweak0, tweak1 = classical_rm_256(state, tweak0, tweak1, K0, K1, alpha)
    if trace:
        trace_dict["rm_state_after"] = bytes_to_hex(state)
        trace_dict["rm_tweak0_after"] = bytes_to_hex(tweak0)
        trace_dict["rm_tweak1_after"] = bytes_to_hex(tweak1)
        
    # === 4. Backward Rounds (Rb) ===
    beta = 3
    for i in range(1, beta + 1):
        state, tweak0, tweak1 = classical_rb_iteration_256(state, tweak0, tweak1, K0, K1, alpha, i, beta)
        if trace:
            trace_dict[f"rb_iteration_{i}_state"] = bytes_to_hex(state)
            trace_dict[f"rb_iteration_{i}_tweak0"] = bytes_to_hex(tweak0)
            trace_dict[f"rb_iteration_{i}_tweak1"] = bytes_to_hex(tweak1)
            
    # === 5. Post-whitening ===
    for j in range(16):
        state[j] = classical_sbox_8bit(state[j], j)
    tweak1_perm = [tweak1[idx] for idx in pi_m]
    for j in range(16):
        state[j] ^= tweak1_perm[j] ^ K0[j] ^ K1[j]
        
    ciphertext = bytes_to_hex(state)
    if trace:
        trace_dict["final_ciphertext"] = ciphertext
        return ciphertext, trace_dict
    return ciphertext

def classical_dialga256d_decrypt(ciphertext_hex, tweak_hex, key_hex):
    """Complete classical reference decryption of Dialga-256 dagger variant."""
    validate_hex_inputs_256(ciphertext_hex, tweak_hex, key_hex)
    
    state = hex_to_bytes(ciphertext_hex)
    tweak = hex_to_bytes(tweak_hex)
    tweak0 = tweak[:16]
    tweak1 = tweak[16:]
    key = hex_to_bytes(key_hex)
    K0 = key[:16]
    K1 = key[16:]
    
    # === 5. Inverse Post-whitening ===
    tweak1_perm = [tweak1[idx] for idx in pi_m]
    for j in range(16):
        state[j] ^= tweak1_perm[j] ^ K0[j] ^ K1[j]
    for j in range(16):
        state[j] = classical_sbox_8bit(state[j], j)
        
    # === 4. Inverse Backward Rounds (Rb) ===
    alpha = 4
    beta = 3
    for i in range(beta, 0, -1):
        if i == beta:
            # Undo tweak1 update
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            for j in range(16):
                tweak1[j] ^= K_tweak_1[j]
                
            # Undo Round 2(alpha+i)+1
            cb_bytes_1 = hex_to_bytes(cb_constants[2 * i - 1])
            tweak0_perm = [tweak0[idx] for idx in pi_m]
            for j in range(16):
                state[j] ^= tweak0_perm[j] ^ cb_bytes_1[j]
            state = classical_round_function_inv(state, 2 * (alpha + i) + 1)
            
            # Undo Round 2(alpha+i)
            cb_bytes = hex_to_bytes(cb_constants[2 * (i - 1)])
            tweak1_perm = [tweak1[idx] for idx in pi_m]
            for j in range(16):
                state[j] ^= tweak1_perm[j] ^ cb_bytes[j]
            state = classical_round_function_inv(state, 2 * (alpha + i))
            
            # Undo tweak0 update
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            for j in range(16):
                tweak0[j] ^= K_tweak_0[j]
        else:
            # Undo tweak1 update
            tweak1 = classical_round_function(tweak1, (alpha - i - 1) % 4)
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            for j in range(16):
                tweak1[j] ^= K_tweak_1[j]
                
            # Undo Round 2(alpha+i)+1
            cb_bytes_1 = hex_to_bytes(cb_constants[2 * i - 1])
            tweak0_perm = [tweak0[idx] for idx in pi_m]
            for j in range(16):
                state[j] ^= tweak0_perm[j] ^ cb_bytes_1[j]
            state = classical_round_function_inv(state, 2 * (alpha + i) + 1)
            
            # Undo Round 2(alpha+i)
            cb_bytes = hex_to_bytes(cb_constants[2 * (i - 1)])
            tweak1_perm = [tweak1[idx] for idx in pi_m]
            for j in range(16):
                state[j] ^= tweak1_perm[j] ^ cb_bytes[j]
            state = classical_round_function_inv(state, 2 * (alpha + i))
            
            # Undo tweak0 update
            tweak0 = classical_round_function(tweak0, (alpha - i - 1) % 4)
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            for j in range(16):
                tweak0[j] ^= K_tweak_0[j]
                
    # === 3. Inverse Middle Rounds (Rm) ===
    # Undo tweak1 update
    tweak1 = classical_round_function(tweak1, (alpha - 1) % 4)
    K_tweak_1 = K0 if (alpha % 2 == 0) else K1
    for j in range(16):
        tweak1[j] ^= K_tweak_1[j]
        
    # Undo Round 2alpha+1
    cm1_bytes = hex_to_bytes(cm_constants[1])
    tweak0_perm = [tweak0[idx] for idx in pi_m]
    for j in range(16):
        state[j] ^= tweak0_perm[j] ^ cm1_bytes[j]
    state = classical_round_function_inv(state, 2 * alpha + 1)
    
    # Undo tweak0 update
    tweak0 = classical_round_function(tweak0, (alpha - 1) % 4)
    K_tweak_0 = K0 if ((alpha - 1) % 2 == 0) else K1
    for j in range(16):
        tweak0[j] ^= K_tweak_0[j]
    for j in range(16):
        tweak0[j] = classical_sbox_8bit(tweak0[j], j)
        
    # Undo Round 2alpha
    cm0_bytes = hex_to_bytes(cm_constants[0])
    for j in range(16):
        state[j] ^= tweak1[j] ^ tweak0[j] ^ cm0_bytes[j]
    state = classical_round_function_inv(state, 2 * alpha)
    
    # Undo tweak0 Sbox prep
    for j in range(16):
        tweak0[j] = classical_sbox_8bit(tweak0[j], j)
        
    # === 2. Inverse Forward Rounds (Rf) ===
    for i in range(alpha, 0, -1):
        if i == 1:
            # Undo tweak1 update
            for j in range(16):
                tweak1[j] ^= K1[j]
                
            # Undo Round 2i-1
            cf_bytes_1 = hex_to_bytes(cf_constants[2 * i - 1])
            for j in range(16):
                state[j] ^= tweak0[j] ^ cf_bytes_1[j]
            state = classical_round_function_inv(state, 2 * i - 1)
            
            # Undo Round 2i-2
            cf_bytes = hex_to_bytes(cf_constants[2 * (i - 1)])
            for j in range(16):
                state[j] ^= tweak1[j] ^ K1[j] ^ cf_bytes[j]
            state = classical_round_function_inv(state, 2 * i - 2)
            
            # Undo tweak0 update
            for j in range(16):
                tweak0[j] ^= K0[j]
        else:
            # Undo tweak1 update
            K_tweak_1 = K1 if (i % 2 == 1) else K0
            for j in range(16):
                tweak1[j] ^= K_tweak_1[j]
            tweak1 = classical_round_function_inv(tweak1, (i - 1) % 4)
            
            # Undo Round 2i-1
            cf_bytes_1 = hex_to_bytes(cf_constants[2 * i - 1])
            for j in range(16):
                state[j] ^= tweak0[j] ^ cf_bytes_1[j]
            state = classical_round_function_inv(state, 2 * i - 1)
            
            # Undo Round 2i-2
            cf_bytes = hex_to_bytes(cf_constants[2 * (i - 1)])
            for j in range(16):
                state[j] ^= tweak1[j] ^ cf_bytes[j]
            state = classical_round_function_inv(state, 2 * i - 2)
            
            # Undo tweak0 update
            K_tweak_0 = K0 if (i % 2 == 1) else K1
            for j in range(16):
                tweak0[j] ^= K_tweak_0[j]
            tweak0 = classical_round_function_inv(tweak0, (i - 1) % 4)
            
    # === 1. Inverse Pre-whitening ===
    for j in range(16):
        state[j] ^= tweak0[j] ^ K0[j] ^ K1[j]
        
    return bytes_to_hex(state)
