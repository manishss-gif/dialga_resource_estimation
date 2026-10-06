# Independent Classical Reference Implementation of Dialga-128
# Based directly on the Dialga specification paper.

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

def validate_hex_inputs(plaintext_hex, tweak_hex, key_hex):
    """Audits input types, sizes, and hexadecimal characters."""
    if not isinstance(plaintext_hex, str) or not isinstance(tweak_hex, str) or not isinstance(key_hex, str):
        raise ValueError("Inputs must be hexadecimal strings.")
        
    for name, hex_str in [("plaintext", plaintext_hex), ("tweak", tweak_hex), ("key", key_hex)]:
        if not all(c in "0123456789abcdefABCDEF" for c in hex_str):
            raise ValueError(f"Input '{name}' must contain only valid hexadecimal characters.")
            
    if len(plaintext_hex) != 32:
        raise ValueError(f"plaintext must be exactly 128 bits (32 hex characters), got {len(plaintext_hex)}")
    if len(tweak_hex) != 32:
        raise ValueError(f"tweak must be exactly 128 bits (32 hex characters), got {len(tweak_hex)}")
    if len(key_hex) != 64:
        raise ValueError(f"key must be exactly 256 bits (64 hex characters), got {len(key_hex)}")

def classical_sbox_8bit(val, cell_idx):
    """Classical mathematical reference of the 8-bit S-box SSb_j."""
    p = pb_tables[cell_idx % 4]
    
    # 1. Convert byte to bits
    in_bits = [(val >> (7 - i)) & 1 for i in range(8)]
    # 2. Permute input bits
    perm_bits = [in_bits[p[i]] for i in range(8)]
    # 3. Split into high and low nibbles
    hi = (perm_bits[0] << 3) | (perm_bits[1] << 2) | (perm_bits[2] << 1) | perm_bits[3]
    lo = (perm_bits[4] << 3) | (perm_bits[5] << 2) | (perm_bits[6] << 1) | perm_bits[7]
    # 4. Apply 4-bit S-box
    hi_out = Sb0[hi]
    lo_out = Sb0[lo]
    # 5. Reassemble bits
    out_bits_unperm = [
        (hi_out >> 3) & 1, (hi_out >> 2) & 1, (hi_out >> 1) & 1, hi_out & 1,
        (lo_out >> 3) & 1, (lo_out >> 2) & 1, (lo_out >> 1) & 1, lo_out & 1
    ]
    # 6. Apply inverse bit permutation
    out_bits = [0] * 8
    for i in range(8):
        out_bits[p[i]] = out_bits_unperm[i]
    # 7. Convert bits back to byte
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
    """Permutes the 16 cells of the state according to pi_table."""
    p = pi_tables[round_idx % 4]
    return [state[idx] for idx in p]

def permute_cells_inv(state, round_idx):
    """Applies the inverse cell permutation."""
    return [state[idx] for idx in pi_tables_inv[round_idx % 4]]

def classical_matrix_mul(state):
    """Classical column-wise matrix multiplication."""
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
    """Applies a single round classically: MatrixMul o Permute o SubCell."""
    state = list(state)
    for cell_idx in range(16):
        state[cell_idx] = classical_sbox_8bit(state[cell_idx], cell_idx)
    state = permute_cells(state, round_idx)
    state = classical_matrix_mul(state)
    return state

def classical_round_function_inv(state, round_idx):
    """Applies single inverse round classically."""
    state = classical_matrix_mul(state)
    state = permute_cells_inv(state, round_idx)
    for cell_idx in range(16):
        state[cell_idx] = classical_sbox_8bit(state[cell_idx], cell_idx)
    return state

# === Modular helper functions for intermediate rounds ===

def classical_rf_iteration(state, tweak, K0, K1, i):
    """Applies iteration i of the forward rounds Rf."""
    # Round 2i-2 on state
    state = classical_round_function(state, 2 * i - 2)
    # XOR key and constant
    K_round = K1 if (i % 2 == 1) else K0
    cf_bytes = hex_to_bytes(cf_constants[2 * (i - 1)])
    for j in range(16):
        state[j] ^= K_round[j] ^ cf_bytes[j]
        
    # Update tweak St
    if i == 1:
        for j in range(16):
            tweak[j] ^= K0[j]
    else:
        tweak = classical_round_function(tweak, (i - 1) % 4)
        K_tweak = K0 if (i % 2 == 1) else K1
        for j in range(16):
            tweak[j] ^= K_tweak[j]
            
    # Round 2i-1 on state
    state = classical_round_function(state, 2 * i - 1)
    cf_bytes_1 = hex_to_bytes(cf_constants[2 * i - 1])
    for j in range(16):
        state[j] ^= tweak[j] ^ cf_bytes_1[j]
        
    return state, tweak

def classical_rm(state, tweak, K0, K1, alpha):
    """Applies the middle rounds Rm."""
    # Tweak Sbox prep
    for j in range(16):
        tweak[j] = classical_sbox_8bit(tweak[j], j)
        
    # Round 2alpha (round 10) on state
    state = classical_round_function(state, 2 * alpha)
    K_rm = K1 if (alpha % 2 == 0) else K0
    cm0_bytes = hex_to_bytes(cm_constants[0])
    for j in range(16):
        state[j] ^= K_rm[j] ^ cm0_bytes[j] ^ tweak[j]
        
    # Tweak Sbox restore
    for j in range(16):
        tweak[j] = classical_sbox_8bit(tweak[j], j)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    for j in range(16):
        tweak[j] ^= K_tweak[j]
    tweak = classical_round_function_inv(tweak, (alpha - 1) % 4)
    
    # Round 2alpha+1 (round 11) on state
    state = classical_round_function(state, 2 * alpha + 1)
    tweak_perm = [tweak[idx] for idx in pi_m]
    cm1_bytes = hex_to_bytes(cm_constants[1])
    for j in range(16):
        state[j] ^= tweak_perm[j] ^ cm1_bytes[j]
        
    return state, tweak

def classical_rb_iteration(state, tweak, K0, K1, alpha, i, beta):
    """Applies iteration i of the backward rounds Rb."""
    # Round 2(alpha+i) on state
    state = classical_round_function(state, 2 * (alpha + i))
    K_round = K0 if ((alpha - i - 1) % 2 == 0) else K1
    cb_bytes = hex_to_bytes(cb_constants[2 * (i - 1)])
    for j in range(16):
        state[j] ^= K_round[j] ^ cb_bytes[j]
        
    # Update tweak state St
    K_tweak = K0 if ((alpha - i - 1) % 2 == 0) else K1
    for j in range(16):
        tweak[j] ^= K_tweak[j]
    if i < beta:
        tweak = classical_round_function_inv(tweak, (alpha - i - 1) % 4)
        
    # Round 2(alpha+i)+1 on state
    state = classical_round_function(state, 2 * (alpha + i) + 1)
    tweak_perm = [tweak[idx] for idx in pi_m]
    cb_bytes_1 = hex_to_bytes(cb_constants[2 * i - 1])
    for j in range(16):
        state[j] ^= tweak_perm[j] ^ cb_bytes_1[j]
        
    return state, tweak

# === Complete Dialga-128 encryption/decryption ===

def classical_dialga128_encrypt(plaintext_hex, tweak_hex, key_hex, trace=False):
    """Complete classical reference encryption of Dialga-128."""
    validate_hex_inputs(plaintext_hex, tweak_hex, key_hex)
    
    state = hex_to_bytes(plaintext_hex)
    tweak = hex_to_bytes(tweak_hex)
    key = hex_to_bytes(key_hex)
    K0 = key[:16]
    K1 = key[16:]
    
    trace_dict = {}
    if trace:
        trace_dict["initial_state"] = bytes_to_hex(state)
        trace_dict["initial_tweak"] = bytes_to_hex(tweak)
        
    # === 1. Pre-whitening ===
    for i in range(16):
        state[i] ^= tweak[i] ^ K0[i] ^ K1[i]
        
    # === 2. Forward Rounds (Rf) ===
    alpha = 5
    for i in range(1, alpha + 1):
        state, tweak = classical_rf_iteration(state, tweak, K0, K1, i)
        if trace:
            trace_dict[f"rf_iteration_{i}_state"] = bytes_to_hex(state)
            trace_dict[f"rf_iteration_{i}_tweak"] = bytes_to_hex(tweak)
            
    # === 3. Middle Rounds (Rm) ===
    state, tweak = classical_rm(state, tweak, K0, K1, alpha)
    if trace:
        trace_dict["rm_state"] = bytes_to_hex(state)
        trace_dict["rm_tweak"] = bytes_to_hex(tweak)
        
    # === 4. Backward Rounds (Rb) ===
    beta = 4
    for i in range(1, beta + 1):
        state, tweak = classical_rb_iteration(state, tweak, K0, K1, alpha, i, beta)
        if trace:
            trace_dict[f"rb_iteration_{i}_state"] = bytes_to_hex(state)
            trace_dict[f"rb_iteration_{i}_tweak"] = bytes_to_hex(tweak)
            
    # === 5. Post-whitening ===
    for j in range(16):
        state[j] = classical_sbox_8bit(state[j], j)
    for j in range(16):
        state[j] ^= K0[j] ^ K1[j]
        
    ciphertext = bytes_to_hex(state)
    if trace:
        trace_dict["final_ciphertext"] = ciphertext
        return ciphertext, trace_dict
    return ciphertext

def classical_dialga128_decrypt(ciphertext_hex, tweak_hex, key_hex):
    """Complete classical reference decryption of Dialga-128."""
    validate_hex_inputs(ciphertext_hex, tweak_hex, key_hex)
    
    state = hex_to_bytes(ciphertext_hex)
    tweak = hex_to_bytes(tweak_hex)
    key = hex_to_bytes(key_hex)
    K0 = key[:16]
    K1 = key[16:]
    
    # === 5. Inverse Post-whitening ===
    for j in range(16):
        state[j] ^= K0[j] ^ K1[j]
    for j in range(16):
        state[j] = classical_sbox_8bit(state[j], j)
        
    # === 4. Inverse Backward Rounds (Rb) ===
    alpha = 5
    beta = 4
    for i in range(beta, 0, -1):
        # Round 2(alpha+i)+1 on state
        cb_bytes_1 = hex_to_bytes(cb_constants[2 * i - 1])
        tweak_perm = [tweak[idx] for idx in pi_m]
        for j in range(16):
            state[j] ^= tweak_perm[j] ^ cb_bytes_1[j]
        state = classical_round_function_inv(state, 2 * (alpha + i) + 1)
        
        # Undo tweak state St update
        if i < beta:
            tweak = classical_round_function(tweak, (alpha - i - 1) % 4)
        K_tweak = K0 if ((alpha - i - 1) % 2 == 0) else K1
        for j in range(16):
            tweak[j] ^= K_tweak[j]
            
        # Round 2(alpha+i) on state
        cb_bytes = hex_to_bytes(cb_constants[2 * (i - 1)])
        K_round = K0 if ((alpha - i - 1) % 2 == 0) else K1
        for j in range(16):
            state[j] ^= K_round[j] ^ cb_bytes[j]
        state = classical_round_function_inv(state, 2 * (alpha + i))
        
    # === 3. Inverse Middle Rounds (Rm) ===
    # Round 2alpha+1 on state
    cm1_bytes = hex_to_bytes(cm_constants[1])
    tweak_perm = [tweak[idx] for idx in pi_m]
    for j in range(16):
        state[j] ^= tweak_perm[j] ^ cm1_bytes[j]
    state = classical_round_function_inv(state, 2 * alpha + 1)
    
    # Undo tweak state St update
    tweak = classical_round_function(tweak, (alpha - 1) % 4)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    for j in range(16):
        tweak[j] ^= K_tweak[j]
    for j in range(16):
        tweak[j] = classical_sbox_8bit(tweak[j], j)
        
    # Round 2alpha on state
    cm0_bytes = hex_to_bytes(cm_constants[0])
    K_rm = K1 if (alpha % 2 == 0) else K0
    for j in range(16):
        state[j] ^= K_rm[j] ^ cm0_bytes[j] ^ tweak[j]
    state = classical_round_function_inv(state, 2 * alpha)
    
    # Undo tweak Sbox prep
    for j in range(16):
        tweak[j] = classical_sbox_8bit(tweak[j], j)
        
    # === 2. Inverse Forward Rounds (Rf) ===
    for i in range(alpha, 0, -1):
        # Round 2i-1 on state
        cf_bytes_1 = hex_to_bytes(cf_constants[2 * i - 1])
        for j in range(16):
            state[j] ^= tweak[j] ^ cf_bytes_1[j]
        state = classical_round_function_inv(state, 2 * i - 1)
        
        # Undo tweak state St update
        if i == 1:
            for j in range(16):
                tweak[j] ^= K0[j]
        else:
            K_tweak = K0 if (i % 2 == 1) else K1
            for j in range(16):
                tweak[j] ^= K_tweak[j]
            tweak = classical_round_function_inv(tweak, (i - 1) % 4)
            
        # Round 2i-2 on state
        cf_bytes = hex_to_bytes(cf_constants[2 * (i - 1)])
        K_round = K1 if (i % 2 == 1) else K0
        for j in range(16):
            state[j] ^= K_round[j] ^ cf_bytes[j]
        state = classical_round_function_inv(state, 2 * i - 2)
        
    # === 1. Inverse Pre-whitening ===
    for j in range(16):
        state[j] ^= tweak[j] ^ K0[j] ^ K1[j]
        
    return bytes_to_hex(state)
