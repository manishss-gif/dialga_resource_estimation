from projectq.ops import CNOT, Toffoli, X, Measure, H, T, Tdag, S, Sdag

# =====================================================================
# SECTION A: SPECIFICATION CONSTANTS AND PERMUTATIONS
# =====================================================================

# Bit permutation tables Pb,j
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

# Pre-computed inverse permutation tables
def invert_permutation(p):
    inv = [0] * len(p)
    for idx, val in enumerate(p):
        inv[val] = idx
    return inv

pi_tables_inv = [invert_permutation(p) for p in pi_tables]

def permute_cells(state, round_idx):
    """Permutes the 16 cells of the state list according to pi_table."""
    p = pi_tables[round_idx % 4]
    return [state[idx] for idx in p]

def permute_cells_inv(state, round_idx):
    """Applies the inverse cell permutation using pre-computed tables."""
    return [state[idx] for idx in pi_tables_inv[round_idx % 4]]

def permute_cells_middle(state):
    """Permutes the 16 cells according to pi_m."""
    return [state[idx] for idx in pi_m]

def matrix_mul(eng, state):
    """Applies the linear layer column-wise multiplication in-place using exactly 8 CNOTs per bit."""
    for col in range(4):
        c0 = state[4 * col]
        c1 = state[4 * col + 1]
        c2 = state[4 * col + 2]
        c3 = state[4 * col + 3]
        
        for bit in range(8):
            CNOT | (c0[bit], c1[bit])
            CNOT | (c1[bit], c0[bit])
            CNOT | (c2[bit], c3[bit])
            CNOT | (c0[bit], c2[bit])
            CNOT | (c3[bit], c0[bit])
            CNOT | (c0[bit], c1[bit])
            CNOT | (c1[bit], c2[bit])
            CNOT | (c2[bit], c3[bit])

def add_key_128(eng, state, key):
    """XORs a 128-bit key register (16 cells of 8 qubits) into the state in-place."""
    for cell_idx in range(16):
        for bit_idx in range(8):
            CNOT | (key[cell_idx][bit_idx], state[cell_idx][bit_idx])

def add_tweak_128(eng, state, tweak):
    """XORs a 128-bit tweak register (16 cells of 8 qubits) into the state in-place."""
    for cell_idx in range(16):
        for bit_idx in range(8):
            CNOT | (tweak[cell_idx][bit_idx], state[cell_idx][bit_idx])

def add_round_constant(eng, state, constant_hex):
    """XORs a 128-bit hexadecimal constant into the 16-cell state."""
    bin_str = bin(int(constant_hex, 16))[2:].zfill(128)
    for cell_idx in range(16):
        cell_bin = bin_str[8 * cell_idx : 8 * cell_idx + 8]
        for bit_idx in range(8):
            if cell_bin[bit_idx] == '1':
                X | state[cell_idx][bit_idx]

def get_engine_list():
    """Returns the list of compilation engines for ProjectQ."""
    from projectq.setups import default
    return default.get_engine_list()


# =====================================================================
# SECTION B: ORDINARY TOFFOLI-LEVEL IMPLEMENTATION
# =====================================================================

def SubCell(eng, x0, x1, x2, x3):
    """Applies the optimized 4-bit S-box in-place using 4 Toffolis, 4 CNOTs, and 2 NOTs. (Ancilla-free)"""
    f = [x2, x1, x3, x0]
    
    # Circuit gates
    CNOT | (f[2], f[1])
    X | f[3]
    Toffoli | (f[2], f[1], f[3])
    Toffoli | (f[3], f[1], f[2])
    CNOT | (f[0], f[3])
    CNOT | (f[2], f[1])
    Toffoli | (f[3], f[1], f[0])
    Toffoli | (f[3], f[0], f[1])
    X | f[1]
    CNOT | (f[0], f[3])
    
    return f[1], f[3], f[2], f[0]

def sbox_8bit(eng, cell, cell_idx):
    """Applies the Toffoli-level 8-bit S-box layer in-place (Ancilla-free)."""
    p = pb_tables[cell_idx % 4]
    p0, p1, p2, p3, p4, p5, p6, p7 = p
    
    # First 4-bit S-box
    c_out1 = SubCell(eng, cell[p0], cell[p1], cell[p2], cell[p3])
    cell[p0], cell[p1], cell[p2], cell[p3] = c_out1
    
    # Second 4-bit S-box
    c_out2 = SubCell(eng, cell[p4], cell[p5], cell[p6], cell[p7])
    cell[p4], cell[p5], cell[p6], cell[p7] = c_out2

def round_function(eng, state, round_idx):
    """Applies a single round of Dialga in-place using ordinary Toffoli S-box."""
    for cell_idx in range(16):
        sbox_8bit(eng, state[cell_idx], cell_idx)
    state[:] = permute_cells(state, round_idx)
    matrix_mul(eng, state)
    return state

def round_function_inv(eng, state, round_idx):
    """Applies a single inverse round of Dialga using ordinary Toffoli S-box."""
    matrix_mul(eng, state)
    state[:] = permute_cells_inv(state, round_idx)
    for cell_idx in range(16):
        sbox_8bit(eng, state[cell_idx], cell_idx)
    return state


# =====================================================================
# SECTION C: CLIFFORD+T LOW-WIDTH IMPLEMENTATION
# =====================================================================

def and_gate(a, b, c, d):
    """Measurement-assisted logical-AND computation: |a,b,0,0> -> |a,b,ab,0>."""
    is_rec = getattr(c.engine.main_engine.backend, 'is_recording', False)
    is_sim = c.engine.main_engine.backend.__class__.__name__ == 'Simulator'
    if not (is_rec or is_sim):
        Toffoli | (a, b, c)
    else:
        H    | c
        CNOT | (b,   d)
        CNOT | (c,   a)
        CNOT | (c,   b)
        CNOT | (a,   d)
        Tdag | a
        Tdag | b
        T    | c        
        T    | d      
        CNOT | (a,   d)
        CNOT | (c,   b)
        CNOT | (c,   a)
        CNOT | (b,   d)
        H    | c
        S    | c


def and_gate_dagger(a, b, c):
    """Measurement-assisted logical-AND uncomputation (cleanup) with feed-forward corrections."""
    is_rec = getattr(c.engine.main_engine.backend, 'is_recording', False)
    is_sim = c.engine.main_engine.backend.__class__.__name__ == 'Simulator'
    if not (is_rec or is_sim):
        Toffoli | (a, b, c)
    else:
        H | c
        Measure | c
        if is_rec or int(c):
            X | c
            S | a
            S | b
            CNOT | (a, b)
            Sdag | b
            CNOT | (a, b)

def SubCellWithAncillaCT_low_width(eng, x0, x1, x2, x3, ancilla, helper):
    """Applies the optimized 4-bit S-box in-place using low-width Clifford+T gates."""
    f = [x2, x1, x3, x0]
    
    # Circuit gates
    CNOT | (f[2], f[1])
    X | f[3]
    
    # 1. Decompose Toffoli | (f[2], f[1], f[3])
    and_gate(f[2], f[1], ancilla, helper)
    CNOT | (ancilla, f[3])
    and_gate_dagger(f[2], f[1], ancilla)
 
    # 2. Decompose Toffoli | (f[3], f[1], f[2])
    and_gate(f[3], f[1], ancilla, helper)
    CNOT | (ancilla, f[2])
    and_gate_dagger(f[3], f[1], ancilla)
 
    CNOT | (f[0], f[3])
    CNOT | (f[2], f[1])
 
    # 3. Decompose Toffoli | (f[3], f[1], f[0])
    and_gate(f[3], f[1], ancilla, helper)
    CNOT | (ancilla, f[0])  
    and_gate_dagger(f[3], f[1], ancilla)
 
    # 4. Decompose Toffoli | (f[3], f[0], f[1])
    and_gate(f[3], f[0], ancilla, helper)
    CNOT | (ancilla, f[1])
    and_gate_dagger(f[3], f[0], ancilla)
 
    X | f[1]
    CNOT | (f[0], f[3])
    
    return f[1], f[3], f[2], f[0]

def sbox_8bitCT(eng, cell, cell_idx, ancilla0, helper0, ancilla1=None, helper1=None):
    """Applies the 8-bit S-box layer in-place. If ancilla1 and helper1 are provided, runs in parallel (low-depth); otherwise sequentially (low-width)."""
    p = pb_tables[cell_idx % 4]
    p0, p1, p2, p3, p4, p5, p6, p7 = p
    
    # First 4-bit S-box
    c_out1 = SubCellWithAncillaCT_low_width(eng, cell[p0], cell[p1], cell[p2], cell[p3], ancilla0, helper0)
    cell[p0], cell[p1], cell[p2], cell[p3] = c_out1
    
    # Second 4-bit S-box
    anc1 = ancilla1 if (ancilla1 is not None) else ancilla0
    hlp1 = helper1 if (helper1 is not None) else helper0
    c_out2 = SubCellWithAncillaCT_low_width(eng, cell[p4], cell[p5], cell[p6], cell[p7], anc1, hlp1)
    cell[p4], cell[p5], cell[p6], cell[p7] = c_out2


def round_functionCT(eng, state, round_idx, ancilla, helper):
    """Applies a single round of Dialga in-place. If ancilla/helper are indexable lists/quregs, runs in parallel (low-depth); otherwise sequentially."""
    is_parallel = hasattr(ancilla, '__len__') or isinstance(ancilla, list)
    
    for cell_idx in range(len(state)):
        if is_parallel:
            sbox_8bitCT(
                eng, state[cell_idx], cell_idx,
                ancilla[2 * cell_idx], helper[2 * cell_idx],
                ancilla[2 * cell_idx + 1], helper[2 * cell_idx + 1]
            )
        else:
            sbox_8bitCT(eng, state[cell_idx], cell_idx, ancilla, helper)
            
    state[:] = permute_cells(state, round_idx)
    matrix_mul(eng, state)
    return state


def round_function_invCT(eng, state, round_idx, ancilla, helper):
    """Applies a single inverse round of Dialga in-place. If ancilla/helper are indexable lists/quregs, runs in parallel (low-depth); otherwise sequentially."""
    is_parallel = hasattr(ancilla, '__len__') or isinstance(ancilla, list)
    
    matrix_mul(eng, state)
    state[:] = permute_cells_inv(state, round_idx)
    
    for cell_idx in range(len(state)):
        if is_parallel:
            sbox_8bitCT(
                eng, state[cell_idx], cell_idx,
                ancilla[2 * cell_idx], helper[2 * cell_idx],
                ancilla[2 * cell_idx + 1], helper[2 * cell_idx + 1]
            )
        else:
            sbox_8bitCT(eng, state[cell_idx], cell_idx, ancilla, helper)
            
    return state


# =====================================================================
# SECTION D: REGISTER HELPERS AND ENCRYPTION SCRIPT ENTRY POINT
# =====================================================================

def allocate_cells(eng, num_cells):
    """Allocates num_cells byte cells of 8 qubits each."""
    return [eng.allocate_qureg(8) for _ in range(num_cells)]


def flatten(regs):
    """Yields every qubit in a nested list of qubits and registers."""
    if hasattr(regs, "id"):
        yield regs
    else:
        for r in regs:
            yield from flatten(r)


def load_hex(hex_str, cells):
    """Sets the computational-basis value of a cell register from a hex string."""
    bits = bin(int(hex_str, 16))[2:].zfill(8 * len(cells))
    for cell_idx, cell in enumerate(cells):
        for bit_idx in range(8):
            if bits[8 * cell_idx + bit_idx] == "1":
                X | cell[bit_idx]


def read_hex(sim, cells):
    """Reads a cell register from a ClassicalSimulator as a hex string."""
    bits = "".join(str(sim.read_bit(q)) for cell in cells for q in cell)
    return hex(int(bits, 2))[2:].zfill(2 * len(cells))


def encryption_main(title, allocate_work, encrypt, kat, argv):
    """Command line of the qdialga*_ct_*.py scripts.

    --test       checks the circuit against the official test vector
    --resource   prints the encryption resources (default)
    """
    from projectq import MainEngine
    from projectq.backends import ClassicalSimulator
    from test_advanced_counter import count_resources, print_resources

    mode = argv[0] if argv else "--resource"
    num_tweaks = len(kat["tweaks"])

    if mode == "--test":
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        state = allocate_cells(eng, 16)
        tweaks = [allocate_cells(eng, 16) for _ in range(num_tweaks)]
        key = allocate_cells(eng, 32)
        work = allocate_work(eng)
        load_hex(kat["plaintext"], state)
        for reg, value in zip(tweaks, kat["tweaks"]):
            load_hex(value, reg)
        load_hex(kat["key"], key)
        eng.flush()

        encrypt(eng, state, tweaks, key, work)
        eng.flush()

        ciphertext = read_hex(sim, state)
        print(f"{title} encryption test vector")
        print(f"    computed   {ciphertext}")
        print(f"    expected   {kat['ciphertext']}")
        assert ciphertext == kat["ciphertext"], "ciphertext mismatch"
        assert all(sim.read_bit(q) == 0 for q in flatten(work)), "work qubits not returned to 0"
        print("    PASS")
        for q in flatten([state, tweaks, key, work]):
            Measure | q
        eng.flush()

    elif mode == "--resource":
        res = count_resources(lambda eng: encrypt(
            eng,
            allocate_cells(eng, 16),
            [allocate_cells(eng, 16) for _ in range(num_tweaks)],
            allocate_cells(eng, 32),
            allocate_work(eng),
        ))
        print_resources(f"{title} encryption resources", res)

    else:
        print("usage: [--test | --resource]")