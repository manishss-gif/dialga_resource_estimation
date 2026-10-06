import unittest
import random
import projectq.backends._sim._pysim
projectq.backends._sim._pysim._USE_REFCHECK = False

from projectq import MainEngine
from projectq.backends import ClassicalSimulator
from projectq.ops import X, Measure, All

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
)

from qdialga128_ct_lw import dialga128_encrypt_ct_low_width
from qdialga128_ct_ld import dialga128_encrypt_ct_low_depth
from qdialga128d_ct_lw import dialga128d_encrypt_ct_low_width
from qdialga128d_ct_ld import dialga128d_encrypt_ct_low_depth
from qdialga256_ct_lw import dialga256_encrypt_ct_low_width
from qdialga256_ct_ld import dialga256_encrypt_ct_low_depth
from qdialga256d_ct_lw import dialga256d_encrypt_ct_low_width
from qdialga256d_ct_ld import dialga256d_encrypt_ct_low_depth

from classical_dialga_reference import classical_dialga128_encrypt, hex_to_bytes, bytes_to_hex, classical_rf_iteration, classical_rm, classical_rb_iteration, classical_sbox_8bit
from classical_dialga256_reference import classical_dialga256_encrypt, classical_dialga256d_encrypt

# Expected intermediate round values from specification
expected_rounds_128d = {
    "r1": "9f95f3ff7a092a2c465dfdf31225ea00",
    "r2": "98871f6b568e38c69b2df2b8fecc46c4",
    "r3": "1d02f8badc70a734a02946d7584c6699",
    "r4": "0924e3460275180560387c9d89ccaef7",
    "r5": "7df1476adedc82fe587839dc0e43f31c",
    "r6": "dd68245888f9cbed8c05065ddcffc56b",
    "r7": "9da77dade1660eaf5a2877126d7bdeeb",
    "r8": "4c31411fc08dd78da9c94db8175ca087",
    "r9": "52ffe58da193d2b26faaa208042c1dfe",
    "r10": "d7beb80383bb057a8a470b61e6f1cd19",
    "r11": "d0f79a0ce3b4f519c8b06141af26a4ca",
    "r12": "ea549d58f095e1159158679953b8be8a",
    "r13": "9ef0091d1fd39ce3f836b8d14d5040eb",
    "r14": "efef7cf1a11eae91f1ed0ac3d8773621",
    "r15": "87e87173805a0eac97de2f230432bcd0",
    "r16": "92a33e3c3115979441131a892119bed7"
}

expected_rounds_256d = {
    "r1": "532ea66649c5b3c29bf7207bed25fb00",
    "r2": "68f551577d90d741e47ef1c492b309b1",
    "r3": "d2f4cd4495ca249118faa7ef77b29488",
    "r4": "019584fc6aadc3dc7159deaf8b9c598d",
    "r5": "a0f832285c988302a78436058aa1b5ea",
    "r6": "28bdf5cffe3cf26b495faf0486e5f6a6",
    "r7": "300aa5780bd7749acecd0bbc5e62300f",
    "r8": "dc6ddc2f5b2bf7778d678b646e6be7be",
    "r9": "c2d77bf88079b4e40658ad332701b32e",
    "r10": "76946ec3a4e948cb34cb8200bf8aef5c",
    "r11": "3ee9efed8d0d12fe76509075814d86ad",
    "r12": "2f4f205788c4e36cd129048945f9abee",
    "r13": "8432386053953d7a5b440167f1b550ef",
    "r14": "17c2edb3182d494b264a9b154d3952f7",
    "r15": "a15eaf8e0aa0d125c9714c300d5b8436",
    "r16": "4daa1d40e36de6bdda58801d83a4aa6b"
}


# Classical reference for Dialga-128d (16 rounds)
def classical_dialga128d_encrypt(plaintext_hex, tweak_hex, key_hex):
    state = hex_to_bytes(plaintext_hex)
    tweak = hex_to_bytes(tweak_hex)
    key = hex_to_bytes(key_hex)
    K0 = key[:16]
    K1 = key[16:]
    
    # 1. Pre-whitening
    for i in range(16):
        state[i] ^= tweak[i] ^ K0[i] ^ K1[i]
        
    # 2. Forward Rounds (Rf)
    alpha = 4
    for i in range(1, alpha + 1):
        state, tweak = classical_rf_iteration(state, tweak, K0, K1, i)
        
    # 3. Middle Rounds (Rm)
    state, tweak = classical_rm(state, tweak, K0, K1, alpha)
    
    # 4. Backward Rounds (Rb)
    beta = 3
    for i in range(1, beta + 1):
        state, tweak = classical_rb_iteration(state, tweak, K0, K1, alpha, i, beta)
        
    # 5. Post-whitening
    for j in range(16):
        state[j] = classical_sbox_8bit(state[j], j)
    for j in range(16):
        state[j] ^= K0[j] ^ K1[j]
        
    return bytes_to_hex(state)


# Helper functions for quantum testing
def load_hex_to_state(hex_str, state_reg):
    bin_str = bin(int(hex_str, 16))[2:].zfill(len(hex_str) * 4)
    for cell_idx in range(len(state_reg)):
        for bit_idx in range(8):
            val = int(bin_str[8 * cell_idx + bit_idx])
            if val == 1:
                X | state_reg[cell_idx][bit_idx]


def read_state_as_hex(sim, state_reg):
    bin_str = ""
    for cell_idx in range(len(state_reg)):
        for bit_idx in range(8):
            bin_str += str(sim.read_bit(state_reg[cell_idx][bit_idx]))
    return hex(int(bin_str, 2))[2:].zfill(len(state_reg) * 2)


def measure_all(*regs):
    from projectq.ops import Measure
    for r in regs:
        if isinstance(r, list):
            for item in r:
                if isinstance(item, list) or hasattr(item, '__len__'):
                    for q in item:
                        Measure | q
                else:
                    Measure | item
        else:
            Measure | r


# Quantum intermediate-state tracing wrappers
def dialga128d_encrypt_traced_ct(eng, sim, state, tweak, key, state_anc, state_hlp, tweak_anc, tweak_hlp):
    K0 = key[:16]
    K1 = key[16:]
    trace = {}
    
    # 1. Pre-whitening
    add_tweak_128(eng, state, tweak)
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    # 2. Forward Rounds (Rf)
    alpha = 4
    for i in range(1, alpha + 1):
        state[:] = round_functionCT(eng, state, 2 * i - 2, state_anc, state_hlp)
        K_round = K1 if (i % 2 == 1) else K0
        add_key_128(eng, state, K_round)
        add_round_constant(eng, state, cf_constants[2 * (i - 1)])
        eng.flush()
        trace[f"r{2*i - 1}"] = read_state_as_hex(sim, state)
        
        if i == 1:
            add_key_128(eng, tweak, K0)
        else:
            tweak[:] = round_functionCT(eng, tweak, (i - 1) % 4, tweak_anc, tweak_hlp)
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak, K_tweak)
            
        state[:] = round_functionCT(eng, state, 2 * i - 1, state_anc, state_hlp)
        add_tweak_128(eng, state, tweak)
        add_round_constant(eng, state, cf_constants[2 * i - 1])
        eng.flush()
        trace[f"r{2*i}"] = read_state_as_hex(sim, state)
        
    # 3. Middle Rounds (Rm)
    state[:] = round_functionCT(eng, state, 2 * alpha, state_anc, state_hlp)
    K_rm = K1 if (alpha % 2 == 0) else K0
    add_key_128(eng, state, K_rm)
    add_round_constant(eng, state, cm_constants[0])
    
    # Tweak SubCell XOR
    for cell_idx in range(16):
        if isinstance(tweak_anc, list) or hasattr(tweak_anc, '__len__'):
            sbox_8bitCT(eng, tweak[cell_idx], cell_idx, tweak_anc[2*cell_idx], tweak_hlp[2*cell_idx], tweak_anc[2*cell_idx+1], tweak_hlp[2*cell_idx+1])
        else:
            sbox_8bitCT(eng, tweak[cell_idx], cell_idx, tweak_anc, tweak_hlp)
    add_tweak_128(eng, state, tweak)
    for cell_idx in range(16):
        if isinstance(tweak_anc, list) or hasattr(tweak_anc, '__len__'):
            sbox_8bitCT(eng, tweak[cell_idx], cell_idx, tweak_anc[2*cell_idx], tweak_hlp[2*cell_idx], tweak_anc[2*cell_idx+1], tweak_hlp[2*cell_idx+1])
        else:
            sbox_8bitCT(eng, tweak[cell_idx], cell_idx, tweak_anc, tweak_hlp)
    eng.flush()
    trace["r9"] = read_state_as_hex(sim, state)
    
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak, K_tweak)
    tweak[:] = round_function_invCT(eng, tweak, (alpha - 1) % 4, tweak_anc, tweak_hlp)
    
    state[:] = round_functionCT(eng, state, 2 * alpha + 1, state_anc, state_hlp)
    tweak_perm = permute_cells_middle(tweak)
    add_tweak_128(eng, state, tweak_perm)
    add_round_constant(eng, state, cm_constants[1])
    eng.flush()
    trace["r10"] = read_state_as_hex(sim, state)
    
    # 4. Backward Rounds (Rb)
    beta = 3
    for i in range(1, beta + 1):
        state[:] = round_functionCT(eng, state, 2 * (alpha + i), state_anc, state_hlp)
        K_round = K0 if ((alpha - i - 1) % 2 == 0) else K1
        add_key_128(eng, state, K_round)
        add_round_constant(eng, state, cb_constants[2 * (i - 1)])
        eng.flush()
        trace[f"r{2*alpha + 2*i + 1}"] = read_state_as_hex(sim, state)
        
        K_tweak = K0 if ((alpha - i - 1) % 2 == 0) else K1
        add_key_128(eng, tweak, K_tweak)
        if i < beta:
            tweak[:] = round_function_invCT(eng, tweak, (alpha - i - 1) % 4, tweak_anc, tweak_hlp)
            
        state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, state_anc, state_hlp)
        tweak_perm = permute_cells_middle(tweak)
        add_tweak_128(eng, state, tweak_perm)
        add_round_constant(eng, state, cb_constants[2 * i - 1])
        eng.flush()
        trace[f"r{2*alpha + 2*i + 2}"] = read_state_as_hex(sim, state)
        
    return trace


def dialga256d_encrypt_traced_ct(eng, sim, state, tweak0, tweak1, key, state_anc, state_hlp, tweak_anc, tweak_hlp):
    K0 = key[:16]
    K1 = key[16:]
    trace = {}
    
    # 1. Pre-whitening
    add_tweak_128(eng, state, tweak0)
    add_key_128(eng, state, K0)
    add_key_128(eng, state, K1)
    
    # 2. Forward Rounds (Rf)
    alpha = 4
    for i in range(1, alpha + 1):
        if i == 1:
            add_key_128(eng, tweak0, K0)
            state[:] = round_functionCT(eng, state, 2 * i - 2, state_anc, state_hlp)
            add_tweak_128(eng, state, tweak1)
            add_key_128(eng, state, K1)
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            eng.flush()
            trace[f"r{2*i - 1}"] = read_state_as_hex(sim, state)
            
            state[:] = round_functionCT(eng, state, 2 * i - 1, state_anc, state_hlp)
            add_tweak_128(eng, state, tweak0)
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            add_key_128(eng, tweak1, K1)
            eng.flush()
            trace[f"r{2*i}"] = read_state_as_hex(sim, state)
        else:
            tweak0[:] = round_functionCT(eng, tweak0, (i - 1) % 4, tweak_anc, tweak_hlp)
            K_tweak = K0 if (i % 2 == 1) else K1
            add_key_128(eng, tweak0, K_tweak)
            
            state[:] = round_functionCT(eng, state, 2 * i - 2, state_anc, state_hlp)
            add_tweak_128(eng, state, tweak1)
            add_round_constant(eng, state, cf_constants[2 * (i - 1)])
            eng.flush()
            trace[f"r{2*i - 1}"] = read_state_as_hex(sim, state)
            
            state[:] = round_functionCT(eng, state, 2 * i - 1, state_anc, state_hlp)
            add_tweak_128(eng, state, tweak0)
            add_round_constant(eng, state, cf_constants[2 * i - 1])
            
            tweak1[:] = round_functionCT(eng, tweak1, (i - 1) % 4, tweak_anc, tweak_hlp)
            K_tweak_1 = K1 if (i % 2 == 1) else K0
            add_key_128(eng, tweak1, K_tweak_1)
            eng.flush()
            trace[f"r{2*i}"] = read_state_as_hex(sim, state)
            
    # 3. Middle Rounds (Rm)
    for cell_idx in range(16):
        if isinstance(tweak_anc, list) or hasattr(tweak_anc, '__len__'):
            sbox_8bitCT(eng, tweak0[cell_idx], cell_idx, tweak_anc[2*cell_idx], tweak_hlp[2*cell_idx], tweak_anc[2*cell_idx+1], tweak_hlp[2*cell_idx+1])
        else:
            sbox_8bitCT(eng, tweak0[cell_idx], cell_idx, tweak_anc, tweak_hlp)
            
    state[:] = round_functionCT(eng, state, 2 * alpha, state_anc, state_hlp)
    add_tweak_128(eng, state, tweak1)
    add_tweak_128(eng, state, tweak0)
    add_round_constant(eng, state, cm_constants[0])
    eng.flush()
    trace["r9"] = read_state_as_hex(sim, state)
    
    for cell_idx in range(16):
        if isinstance(tweak_anc, list) or hasattr(tweak_anc, '__len__'):
            sbox_8bitCT(eng, tweak0[cell_idx], cell_idx, tweak_anc[2*cell_idx], tweak_hlp[2*cell_idx], tweak_anc[2*cell_idx+1], tweak_hlp[2*cell_idx+1])
        else:
            sbox_8bitCT(eng, tweak0[cell_idx], cell_idx, tweak_anc, tweak_hlp)
    K_tweak = K0 if ((alpha - 1) % 2 == 0) else K1
    add_key_128(eng, tweak0, K_tweak)
    tweak0[:] = round_function_invCT(eng, tweak0, (alpha - 1) % 4, tweak_anc, tweak_hlp)
    
    state[:] = round_functionCT(eng, state, 2 * alpha + 1, state_anc, state_hlp)
    tweak0_perm = permute_cells_middle(tweak0)
    add_tweak_128(eng, state, tweak0_perm)
    add_round_constant(eng, state, cm_constants[1])
    
    K_tweak_1 = K0 if (alpha % 2 == 0) else K1
    add_key_128(eng, tweak1, K_tweak_1)
    tweak1[:] = round_function_invCT(eng, tweak1, (alpha - 1) % 4, tweak_anc, tweak_hlp)
    eng.flush()
    trace["r10"] = read_state_as_hex(sim, state)
    
    # 4. Backward Rounds (Rb)
    beta = 3
    for i in range(1, beta + 1):
        if i == beta:
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
            
            state[:] = round_functionCT(eng, state, 2 * (alpha + i), state_anc, state_hlp)
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            eng.flush()
            trace[f"r{2*alpha + 2*i + 1}"] = read_state_as_hex(sim, state)
            
            state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, state_anc, state_hlp)
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
            eng.flush()
            trace[f"r{2*alpha + 2*i + 2}"] = read_state_as_hex(sim, state)
        else:
            K_tweak_0 = K0 if ((alpha - i - 1) % 2 == 0) else K1
            add_key_128(eng, tweak0, K_tweak_0)
            tweak0[:] = round_function_invCT(eng, tweak0, (alpha - i - 1) % 4, tweak_anc, tweak_hlp)
            
            state[:] = round_functionCT(eng, state, 2 * (alpha + i), state_anc, state_hlp)
            tweak1_perm = permute_cells_middle(tweak1)
            add_tweak_128(eng, state, tweak1_perm)
            add_round_constant(eng, state, cb_constants[2 * (i - 1)])
            eng.flush()
            trace[f"r{2*alpha + 2*i + 1}"] = read_state_as_hex(sim, state)
            
            state[:] = round_functionCT(eng, state, 2 * (alpha + i) + 1, state_anc, state_hlp)
            tweak0_perm = permute_cells_middle(tweak0)
            add_tweak_128(eng, state, tweak0_perm)
            add_round_constant(eng, state, cb_constants[2 * i - 1])
            
            K_tweak_1 = K0 if ((alpha - i) % 2 == 0) else K1
            add_key_128(eng, tweak1, K_tweak_1)
            tweak1[:] = round_function_invCT(eng, tweak1, (alpha - i - 1) % 4, tweak_anc, tweak_hlp)
            eng.flush()
            trace[f"r{2*alpha + 2*i + 2}"] = read_state_as_hex(sim, state)
            
    return trace


class TestCTProfiles(unittest.TestCase):
    # ========================================================
    # PART E: OFFICIAL KATs VALIDATION FOR ALL VARIANTS
    # ========================================================
    def test_dialga128_kats_lw(self):
        """Dialga-128 official KATs LW"""
        kats = [
            ("00112233445566778899aabbccddeeff", "2233445566778899aabbccddee00ff11", "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00", "a4a1ea948919d8996e13b1b365bb0ce6"),
            ("0123456789abcdef0123456789abcdef", "00001111222233334444555566668888", "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210", "dc355a6376d9617723efb9a98c1b4864")
        ]
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        for p, t, k, expected_c in kats:
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            ancilla = eng.allocate_qubit()[0]
            helper = eng.allocate_qubit()[0]
            
            load_hex_to_state(p, state)
            load_hex_to_state(t, tweak)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga128_encrypt_ct_low_width(eng, state, tweak, key, ancilla, helper)
            eng.flush()
            
            # Assert cleanup before deallocate
            self.assertEqual(sim.read_bit(ancilla), 0)
            self.assertEqual(sim.read_bit(helper), 0)
            
            c = read_state_as_hex(sim, state)
            self.assertEqual(c, expected_c)
            
            measure_all(state, tweak, key, ancilla, helper)
            eng.flush()

    def test_dialga128_kats_ld(self):
        """Dialga-128 official KATs LD"""
        kats = [
            ("00112233445566778899aabbccddeeff", "2233445566778899aabbccddee00ff11", "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00", "a4a1ea948919d8996e13b1b365bb0ce6"),
            ("0123456789abcdef0123456789abcdef", "00001111222233334444555566668888", "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210", "dc355a6376d9617723efb9a98c1b4864")
        ]
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        for p, t, k, expected_c in kats:
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            state_anc = eng.allocate_qureg(32)
            state_hlp = eng.allocate_qureg(32)
            tweak_anc = eng.allocate_qureg(32)
            tweak_hlp = eng.allocate_qureg(32)
            
            load_hex_to_state(p, state)
            load_hex_to_state(t, tweak)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga128_encrypt_ct_low_depth(eng, state, tweak, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            
            # Assert cleanup before deallocate
            for q in list(state_anc) + list(state_hlp) + list(tweak_anc) + list(tweak_hlp):
                self.assertEqual(sim.read_bit(q), 0)
                
            c = read_state_as_hex(sim, state)
            self.assertEqual(c, expected_c)
            
            measure_all(state, tweak, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()

    def test_dialga128d_kats_lw(self):
        """Dialga-128† official KATs LW"""
        kats = [
            ("00112233445566778899aabbccddeeff", "2233445566778899aabbccddee00ff11", "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00", "838407143af9a876fbdc6be378e9045b"),
            ("0123456789abcdef0123456789abcdef", "00001111222233334444555566668888", "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210", "a16812d9738333d238c23e67ac20ef16")
        ]
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        for p, t, k, expected_c in kats:
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            ancilla = eng.allocate_qubit()[0]
            helper = eng.allocate_qubit()[0]
            
            load_hex_to_state(p, state)
            load_hex_to_state(t, tweak)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga128d_encrypt_ct_low_width(eng, state, tweak, key, ancilla, helper)
            eng.flush()
            
            self.assertEqual(sim.read_bit(ancilla), 0)
            self.assertEqual(sim.read_bit(helper), 0)
            
            c = read_state_as_hex(sim, state)
            self.assertEqual(c, expected_c)
            
            measure_all(state, tweak, key, ancilla, helper)
            eng.flush()

    def test_dialga128d_kats_ld(self):
        """Dialga-128† official KATs LD"""
        kats = [
            ("00112233445566778899aabbccddeeff", "2233445566778899aabbccddee00ff11", "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00", "838407143af9a876fbdc6be378e9045b"),
            ("0123456789abcdef0123456789abcdef", "00001111222233334444555566668888", "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210", "a16812d9738333d238c23e67ac20ef16")
        ]
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        for p, t, k, expected_c in kats:
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            state_anc = eng.allocate_qureg(32)
            state_hlp = eng.allocate_qureg(32)
            tweak_anc = eng.allocate_qureg(32)
            tweak_hlp = eng.allocate_qureg(32)
            
            load_hex_to_state(p, state)
            load_hex_to_state(t, tweak)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga128d_encrypt_ct_low_depth(eng, state, tweak, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            
            for q in list(state_anc) + list(state_hlp) + list(tweak_anc) + list(tweak_hlp):
                self.assertEqual(sim.read_bit(q), 0)
                
            c = read_state_as_hex(sim, state)
            self.assertEqual(c, expected_c)
            
            measure_all(state, tweak, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()

    def test_dialga256_kats_lw(self):
        """Dialga-256 official KATs LW"""
        kats = [
            ("00112233445566778899aabbccddeeff", "2233445566778899aabbccddeeff1122", "33445566778899aabbccddeeff001100", "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00", "e129ff0920371753db65532540d06881"),
            ("0123456789abcdef0123456789abcdef", "00001111222233334444555566667777", "88889999aaaabbbbccccddddeeeeffff", "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210", "b30b9d8a8d6dd98d235a4aaeda5b06ca")
        ]
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        for p, t0, t1, k, expected_c in kats:
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
            tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            ancilla = eng.allocate_qubit()[0]
            helper = eng.allocate_qubit()[0]
            
            load_hex_to_state(p, state)
            load_hex_to_state(t0, tweak0)
            load_hex_to_state(t1, tweak1)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga256_encrypt_ct_low_width(eng, state, tweak0, tweak1, key, ancilla, helper)
            eng.flush()
            
            self.assertEqual(sim.read_bit(ancilla), 0)
            self.assertEqual(sim.read_bit(helper), 0)
            
            c = read_state_as_hex(sim, state)
            self.assertEqual(c, expected_c)
            
            measure_all(state, tweak0, tweak1, key, ancilla, helper)
            eng.flush()

    def test_dialga256_kats_ld(self):
        """Dialga-256 official KATs LD"""
        kats = [
            ("00112233445566778899aabbccddeeff", "2233445566778899aabbccddeeff1122", "33445566778899aabbccddeeff001100", "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00", "e129ff0920371753db65532540d06881"),
            ("0123456789abcdef0123456789abcdef", "00001111222233334444555566667777", "88889999aaaabbbbccccddddeeeeffff", "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210", "b30b9d8a8d6dd98d235a4aaeda5b06ca")
        ]
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        for p, t0, t1, k, expected_c in kats:
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
            tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            state_anc = eng.allocate_qureg(32)
            state_hlp = eng.allocate_qureg(32)
            tweak_anc = eng.allocate_qureg(32)
            tweak_hlp = eng.allocate_qureg(32)
            
            load_hex_to_state(p, state)
            load_hex_to_state(t0, tweak0)
            load_hex_to_state(t1, tweak1)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga256_encrypt_ct_low_depth(eng, state, tweak0, tweak1, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            
            for q in list(state_anc) + list(state_hlp) + list(tweak_anc) + list(tweak_hlp):
                self.assertEqual(sim.read_bit(q), 0)
                
            c = read_state_as_hex(sim, state)
            self.assertEqual(c, expected_c)
            
            measure_all(state, tweak0, tweak1, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()

    def test_dialga256d_kats_lw(self):
        """Dialga-256† official KATs LW"""
        kats = [
            ("00112233445566778899aabbccddeeff", "2233445566778899aabbccddeeff1122", "33445566778899aabbccddeeff001100", "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00", "833367ce5fd49a7b1f0dd9ec62720018"),
            ("0123456789abcdef0123456789abcdef", "00001111222233334444555566667777", "88889999aaaabbbbccccddddeeeeffff", "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210", "0f9991a2b1f311c8c786ad0434edc15a")
        ]
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        for p, t0, t1, k, expected_c in kats:
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
            tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            ancilla = eng.allocate_qubit()[0]
            helper = eng.allocate_qubit()[0]
            
            load_hex_to_state(p, state)
            load_hex_to_state(t0, tweak0)
            load_hex_to_state(t1, tweak1)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga256d_encrypt_ct_low_width(eng, state, tweak0, tweak1, key, ancilla, helper)
            eng.flush()
            
            self.assertEqual(sim.read_bit(ancilla), 0)
            self.assertEqual(sim.read_bit(helper), 0)
            
            c = read_state_as_hex(sim, state)
            self.assertEqual(c, expected_c)
            
            measure_all(state, tweak0, tweak1, key, ancilla, helper)
            eng.flush()

    def test_dialga256d_kats_ld(self):
        """Dialga-256† official KATs LD"""
        kats = [
            ("00112233445566778899aabbccddeeff", "2233445566778899aabbccddeeff1122", "33445566778899aabbccddeeff001100", "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00", "833367ce5fd49a7b1f0dd9ec62720018"),
            ("0123456789abcdef0123456789abcdef", "00001111222233334444555566667777", "88889999aaaabbbbccccddddeeeeffff", "fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210", "0f9991a2b1f311c8c786ad0434edc15a")
        ]
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        for p, t0, t1, k, expected_c in kats:
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
            tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            state_anc = eng.allocate_qureg(32)
            state_hlp = eng.allocate_qureg(32)
            tweak_anc = eng.allocate_qureg(32)
            tweak_hlp = eng.allocate_qureg(32)
            
            load_hex_to_state(p, state)
            load_hex_to_state(t0, tweak0)
            load_hex_to_state(t1, tweak1)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga256d_encrypt_ct_low_depth(eng, state, tweak0, tweak1, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            
            for q in list(state_anc) + list(state_hlp) + list(tweak_anc) + list(tweak_hlp):
                self.assertEqual(sim.read_bit(q), 0)
                
            c = read_state_as_hex(sim, state)
            self.assertEqual(c, expected_c)
            
            measure_all(state, tweak0, tweak1, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()

    # ========================================================
    # PART F: RANDOM DIFFERENTIAL TESTING
    # ========================================================
    def test_dialga128_random(self):
        """Dialga-128 random equivalence"""
        rng = random.Random(20260807)
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        
        for idx in range(100):
            p = f"{rng.getrandbits(128):032x}"
            t = f"{rng.getrandbits(128):032x}"
            k = f"{rng.getrandbits(256):064x}"
            
            expected_c = classical_dialga128_encrypt(p, t, k)
            
            # LW
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            ancilla = eng.allocate_qubit()[0]
            helper = eng.allocate_qubit()[0]
            
            load_hex_to_state(p, state)
            load_hex_to_state(t, tweak)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga128_encrypt_ct_low_width(eng, state, tweak, key, ancilla, helper)
            eng.flush()
            c_lw = read_state_as_hex(sim, state)
            
            # Assert cleanup before deallocate
            anc_val = sim.read_bit(ancilla)
            hlp_val = sim.read_bit(helper)
            if anc_val != 0:
                self.fail(f"Cleanup check failed: variant=dialga128, profile=LW, test-case index={idx}, workspace register=ancilla, observed value={anc_val}")
            if hlp_val != 0:
                self.fail(f"Cleanup check failed: variant=dialga128, profile=LW, test-case index={idx}, workspace register=helper, observed value={hlp_val}")
            
            measure_all(state, tweak, key, ancilla, helper)
            eng.flush()
            
            # LD
            key_ld = [eng.allocate_qureg(8) for _ in range(32)]
            tweak_ld = [eng.allocate_qureg(8) for _ in range(16)]
            state_ld = [eng.allocate_qureg(8) for _ in range(16)]
            state_anc = eng.allocate_qureg(32)
            state_hlp = eng.allocate_qureg(32)
            tweak_anc = eng.allocate_qureg(32)
            tweak_hlp = eng.allocate_qureg(32)
            
            load_hex_to_state(p, state_ld)
            load_hex_to_state(t, tweak_ld)
            load_hex_to_state(k, key_ld)
            eng.flush()
            
            dialga128_encrypt_ct_low_depth(eng, state_ld, tweak_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            c_ld = read_state_as_hex(sim, state_ld)
            
            # Assert cleanup before deallocate
            for q_idx, q in enumerate(state_anc):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga128, profile=LD, test-case index={idx}, workspace register=state_anc, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(state_hlp):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga128, profile=LD, test-case index={idx}, workspace register=state_hlp, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(tweak_anc):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga128, profile=LD, test-case index={idx}, workspace register=tweak_anc, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(tweak_hlp):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga128, profile=LD, test-case index={idx}, workspace register=tweak_hlp, workspace index={q_idx}, observed value={val}")
            
            measure_all(state_ld, tweak_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            
            self.assertEqual(c_lw, expected_c)
            self.assertEqual(c_ld, expected_c)

    def test_dialga128d_random(self):
        """Dialga-128† random equivalence"""
        rng = random.Random(20260807)
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        
        for idx in range(100):
            p = f"{rng.getrandbits(128):032x}"
            t = f"{rng.getrandbits(128):032x}"
            k = f"{rng.getrandbits(256):064x}"
            
            expected_c = classical_dialga128d_encrypt(p, t, k)
            
            # LW
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            ancilla = eng.allocate_qubit()[0]
            helper = eng.allocate_qubit()[0]
            
            load_hex_to_state(p, state)
            load_hex_to_state(t, tweak)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga128d_encrypt_ct_low_width(eng, state, tweak, key, ancilla, helper)
            eng.flush()
            c_lw = read_state_as_hex(sim, state)
            
            # Assert cleanup
            anc_val = sim.read_bit(ancilla)
            hlp_val = sim.read_bit(helper)
            if anc_val != 0:
                self.fail(f"Cleanup check failed: variant=dialga128d, profile=LW, test-case index={idx}, workspace register=ancilla, observed value={anc_val}")
            if hlp_val != 0:
                self.fail(f"Cleanup check failed: variant=dialga128d, profile=LW, test-case index={idx}, workspace register=helper, observed value={hlp_val}")
            
            measure_all(state, tweak, key, ancilla, helper)
            eng.flush()
            
            # LD
            key_ld = [eng.allocate_qureg(8) for _ in range(32)]
            tweak_ld = [eng.allocate_qureg(8) for _ in range(16)]
            state_ld = [eng.allocate_qureg(8) for _ in range(16)]
            state_anc = eng.allocate_qureg(32)
            state_hlp = eng.allocate_qureg(32)
            tweak_anc = eng.allocate_qureg(32)
            tweak_hlp = eng.allocate_qureg(32)
            
            load_hex_to_state(p, state_ld)
            load_hex_to_state(t, tweak_ld)
            load_hex_to_state(k, key_ld)
            eng.flush()
            
            dialga128d_encrypt_ct_low_depth(eng, state_ld, tweak_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            c_ld = read_state_as_hex(sim, state_ld)
            
            # Assert cleanup
            for q_idx, q in enumerate(state_anc):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga128d, profile=LD, test-case index={idx}, workspace register=state_anc, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(state_hlp):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga128d, profile=LD, test-case index={idx}, workspace register=state_hlp, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(tweak_anc):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga128d, profile=LD, test-case index={idx}, workspace register=tweak_anc, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(tweak_hlp):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga128d, profile=LD, test-case index={idx}, workspace register=tweak_hlp, workspace index={q_idx}, observed value={val}")
            
            measure_all(state_ld, tweak_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            
            self.assertEqual(c_lw, expected_c)
            self.assertEqual(c_ld, expected_c)

    def test_dialga256_random(self):
        """Dialga-256 random equivalence"""
        rng = random.Random(20260807)
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        
        for idx in range(100):
            p = f"{rng.getrandbits(128):032x}"
            t0 = f"{rng.getrandbits(128):032x}"
            t1 = f"{rng.getrandbits(128):032x}"
            k = f"{rng.getrandbits(256):064x}"
            
            tweak_full = t0 + t1
            expected_c = classical_dialga256_encrypt(p, tweak_full, k)
            
            # LW
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
            tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            ancilla = eng.allocate_qubit()[0]
            helper = eng.allocate_qubit()[0]
            
            load_hex_to_state(p, state)
            load_hex_to_state(t0, tweak0)
            load_hex_to_state(t1, tweak1)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga256_encrypt_ct_low_width(eng, state, tweak0, tweak1, key, ancilla, helper)
            eng.flush()
            c_lw = read_state_as_hex(sim, state)
            
            # Assert cleanup
            anc_val = sim.read_bit(ancilla)
            hlp_val = sim.read_bit(helper)
            if anc_val != 0:
                self.fail(f"Cleanup check failed: variant=dialga256, profile=LW, test-case index={idx}, workspace register=ancilla, observed value={anc_val}")
            if hlp_val != 0:
                self.fail(f"Cleanup check failed: variant=dialga256, profile=LW, test-case index={idx}, workspace register=helper, observed value={hlp_val}")
            
            measure_all(state, tweak0, tweak1, key, ancilla, helper)
            eng.flush()
            
            # LD
            key_ld = [eng.allocate_qureg(8) for _ in range(32)]
            tweak0_ld = [eng.allocate_qureg(8) for _ in range(16)]
            tweak1_ld = [eng.allocate_qureg(8) for _ in range(16)]
            state_ld = [eng.allocate_qureg(8) for _ in range(16)]
            state_anc = eng.allocate_qureg(32)
            state_hlp = eng.allocate_qureg(32)
            tweak_anc = eng.allocate_qureg(32)
            tweak_hlp = eng.allocate_qureg(32)
            
            load_hex_to_state(p, state_ld)
            load_hex_to_state(t0, tweak0_ld)
            load_hex_to_state(t1, tweak1_ld)
            load_hex_to_state(k, key_ld)
            eng.flush()
            
            dialga256_encrypt_ct_low_depth(eng, state_ld, tweak0_ld, tweak1_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            c_ld = read_state_as_hex(sim, state_ld)
            
            # Assert cleanup
            for q_idx, q in enumerate(state_anc):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga256, profile=LD, test-case index={idx}, workspace register=state_anc, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(state_hlp):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga256, profile=LD, test-case index={idx}, workspace register=state_hlp, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(tweak_anc):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga256, profile=LD, test-case index={idx}, workspace register=tweak_anc, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(tweak_hlp):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga256, profile=LD, test-case index={idx}, workspace register=tweak_hlp, workspace index={q_idx}, observed value={val}")
            
            measure_all(state_ld, tweak0_ld, tweak1_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            
            self.assertEqual(c_lw, expected_c)
            self.assertEqual(c_ld, expected_c)

    def test_dialga256d_random(self):
        """Dialga-256† random equivalence"""
        rng = random.Random(20260807)
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        
        for idx in range(100):
            p = f"{rng.getrandbits(128):032x}"
            t0 = f"{rng.getrandbits(128):032x}"
            t1 = f"{rng.getrandbits(128):032x}"
            k = f"{rng.getrandbits(256):064x}"
            
            tweak_full = t0 + t1
            expected_c = classical_dialga256d_encrypt(p, tweak_full, k)
            
            # LW
            key = [eng.allocate_qureg(8) for _ in range(32)]
            tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
            tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
            state = [eng.allocate_qureg(8) for _ in range(16)]
            ancilla = eng.allocate_qubit()[0]
            helper = eng.allocate_qubit()[0]
            
            load_hex_to_state(p, state)
            load_hex_to_state(t0, tweak0)
            load_hex_to_state(t1, tweak1)
            load_hex_to_state(k, key)
            eng.flush()
            
            dialga256d_encrypt_ct_low_width(eng, state, tweak0, tweak1, key, ancilla, helper)
            eng.flush()
            c_lw = read_state_as_hex(sim, state)
            
            # Assert cleanup
            anc_val = sim.read_bit(ancilla)
            hlp_val = sim.read_bit(helper)
            if anc_val != 0:
                self.fail(f"Cleanup check failed: variant=dialga256d, profile=LW, test-case index={idx}, workspace register=ancilla, observed value={anc_val}")
            if hlp_val != 0:
                self.fail(f"Cleanup check failed: variant=dialga256d, profile=LW, test-case index={idx}, workspace register=helper, observed value={hlp_val}")
            
            measure_all(state, tweak0, tweak1, key, ancilla, helper)
            eng.flush()
            
            # LD
            key_ld = [eng.allocate_qureg(8) for _ in range(32)]
            tweak0_ld = [eng.allocate_qureg(8) for _ in range(16)]
            tweak1_ld = [eng.allocate_qureg(8) for _ in range(16)]
            state_ld = [eng.allocate_qureg(8) for _ in range(16)]
            state_anc = eng.allocate_qureg(32)
            state_hlp = eng.allocate_qureg(32)
            tweak_anc = eng.allocate_qureg(32)
            tweak_hlp = eng.allocate_qureg(32)
            
            load_hex_to_state(p, state_ld)
            load_hex_to_state(t0, tweak0_ld)
            load_hex_to_state(t1, tweak1_ld)
            load_hex_to_state(k, key_ld)
            eng.flush()
            
            dialga256d_encrypt_ct_low_depth(eng, state_ld, tweak0_ld, tweak1_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            c_ld = read_state_as_hex(sim, state_ld)
            
            # Assert cleanup
            for q_idx, q in enumerate(state_anc):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga256d, profile=LD, test-case index={idx}, workspace register=state_anc, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(state_hlp):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga256d, profile=LD, test-case index={idx}, workspace register=state_hlp, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(tweak_anc):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga256d, profile=LD, test-case index={idx}, workspace register=tweak_anc, workspace index={q_idx}, observed value={val}")
            for q_idx, q in enumerate(tweak_hlp):
                val = sim.read_bit(q)
                if val != 0:
                    self.fail(f"Cleanup check failed: variant=dialga256d, profile=LD, test-case index={idx}, workspace register=tweak_hlp, workspace index={q_idx}, observed value={val}")
            
            measure_all(state_ld, tweak0_ld, tweak1_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            eng.flush()
            
            self.assertEqual(c_lw, expected_c)
            self.assertEqual(c_ld, expected_c)

    # ========================================================
    # PART H: DAGGER INTERMEDIATE-STATE TRACING TESTS
    # ========================================================
    def test_dialga128d_round_trace(self):
        """Dialga-128† round trace"""
        p_hex = "00112233445566778899aabbccddeeff"
        t_hex = "2233445566778899aabbccddee00ff11"
        k_hex = "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00"
        
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        
        # LW Trace
        key = [eng.allocate_qureg(8) for _ in range(32)]
        tweak = [eng.allocate_qureg(8) for _ in range(16)]
        state = [eng.allocate_qureg(8) for _ in range(16)]
        ancilla = eng.allocate_qubit()[0]
        helper = eng.allocate_qubit()[0]
        
        load_hex_to_state(p_hex, state)
        load_hex_to_state(t_hex, tweak)
        load_hex_to_state(k_hex, key)
        eng.flush()
        
        trace_lw = dialga128d_encrypt_traced_ct(eng, sim, state, tweak, key, ancilla, helper, ancilla, helper)
        
        # Verify workspace cleanup for LW
        self.assertEqual(sim.read_bit(ancilla), 0)
        self.assertEqual(sim.read_bit(helper), 0)
        
        measure_all(state, tweak, key, ancilla, helper)
        eng.flush()
        
        # LD Trace
        key_ld = [eng.allocate_qureg(8) for _ in range(32)]
        tweak_ld = [eng.allocate_qureg(8) for _ in range(16)]
        state_ld = [eng.allocate_qureg(8) for _ in range(16)]
        state_anc = eng.allocate_qureg(32)
        state_hlp = eng.allocate_qureg(32)
        tweak_anc = eng.allocate_qureg(32)
        tweak_hlp = eng.allocate_qureg(32)
        
        load_hex_to_state(p_hex, state_ld)
        load_hex_to_state(t_hex, tweak_ld)
        load_hex_to_state(k_hex, key_ld)
        eng.flush()
        
        trace_ld = dialga128d_encrypt_traced_ct(eng, sim, state_ld, tweak_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
        
        # Verify workspace cleanup for LD
        for q in list(state_anc) + list(state_hlp) + list(tweak_anc) + list(tweak_hlp):
            self.assertEqual(sim.read_bit(q), 0)
            
        measure_all(state_ld, tweak_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
        eng.flush()
        
        for r in range(1, 17):
            r_key = f"r{r}"
            self.assertEqual(trace_lw[r_key], expected_rounds_128d[r_key])
            self.assertEqual(trace_ld[r_key], expected_rounds_128d[r_key])

    def test_dialga256d_round_trace(self):
        """Dialga-256† round trace"""
        p_hex = "00112233445566778899aabbccddeeff"
        t_hex = "2233445566778899aabbccddeeff112233445566778899aabbccddeeff001100"
        k_hex = "00112233445566778899aabbccddeeff112233445566778899aabbccddeeff00"
        
        sim = ClassicalSimulator()
        eng = MainEngine(sim, [])
        
        # LW Trace
        key = [eng.allocate_qureg(8) for _ in range(32)]
        tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
        tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
        state = [eng.allocate_qureg(8) for _ in range(16)]
        ancilla = eng.allocate_qubit()[0]
        helper = eng.allocate_qubit()[0]
        
        load_hex_to_state(p_hex, state)
        load_hex_to_state(t_hex[:32], tweak0)
        load_hex_to_state(t_hex[32:], tweak1)
        load_hex_to_state(k_hex, key)
        eng.flush()
        
        trace_lw = dialga256d_encrypt_traced_ct(eng, sim, state, tweak0, tweak1, key, ancilla, helper, ancilla, helper)
        
        # Verify workspace cleanup for LW
        self.assertEqual(sim.read_bit(ancilla), 0)
        self.assertEqual(sim.read_bit(helper), 0)
        
        measure_all(state, tweak0, tweak1, key, ancilla, helper)
        eng.flush()
        
        # LD Trace
        key_ld = [eng.allocate_qureg(8) for _ in range(32)]
        tweak0_ld = [eng.allocate_qureg(8) for _ in range(16)]
        tweak1_ld = [eng.allocate_qureg(8) for _ in range(16)]
        state_ld = [eng.allocate_qureg(8) for _ in range(16)]
        state_anc = eng.allocate_qureg(32)
        state_hlp = eng.allocate_qureg(32)
        tweak_anc = eng.allocate_qureg(32)
        tweak_hlp = eng.allocate_qureg(32)
        
        load_hex_to_state(p_hex, state_ld)
        load_hex_to_state(t_hex[:32], tweak0_ld)
        load_hex_to_state(t_hex[32:], tweak1_ld)
        load_hex_to_state(k_hex, key_ld)
        eng.flush()
        
        trace_ld = dialga256d_encrypt_traced_ct(eng, sim, state_ld, tweak0_ld, tweak1_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
        
        # Verify workspace cleanup for LD
        for q in list(state_anc) + list(state_hlp) + list(tweak_anc) + list(tweak_hlp):
            self.assertEqual(sim.read_bit(q), 0)
            
        measure_all(state_ld, tweak0_ld, tweak1_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
        eng.flush()
        
        for r in range(1, 17):
            r_key = f"r{r}"
            self.assertEqual(trace_lw[r_key], expected_rounds_256d[r_key])
            self.assertEqual(trace_ld[r_key], expected_rounds_256d[r_key])

    # ========================================================
    # PART J & K: INVARIANTS & SANITY CHECKS (RESOURCE INVARIANTS)
    # ========================================================
    def test_gate_count_invariants(self):
        """LD/LW gate-count invariants for all variants"""
        from test_advanced_counter import AdvancedResourceCounter
        
        variants = ["dialga128", "dialga128d", "dialga256", "dialga256d"]
        expected_widths = {
            "dialga128_lw": 514, "dialga128_ld": 640,
            "dialga128d_lw": 514, "dialga128d_ld": 640,
            "dialga256_lw": 642, "dialga256_ld": 768,
            "dialga256d_lw": 642, "dialga256d_ld": 768
        }
        
        for var in variants:
            rc_lw = AdvancedResourceCounter()
            rc_lw.is_recording = True
            eng_lw = MainEngine(rc_lw, [])
            
            if var == "dialga128":
                key = [eng_lw.allocate_qureg(8) for _ in range(32)]
                tweak = [eng_lw.allocate_qureg(8) for _ in range(16)]
                state = [eng_lw.allocate_qureg(8) for _ in range(16)]
                anc = eng_lw.allocate_qubit()[0]
                hlp = eng_lw.allocate_qubit()[0]
                dialga128_encrypt_ct_low_width(eng_lw, state, tweak, key, anc, hlp)
            elif var == "dialga128d":
                key = [eng_lw.allocate_qureg(8) for _ in range(32)]
                tweak = [eng_lw.allocate_qureg(8) for _ in range(16)]
                state = [eng_lw.allocate_qureg(8) for _ in range(16)]
                anc = eng_lw.allocate_qubit()[0]
                hlp = eng_lw.allocate_qubit()[0]
                dialga128d_encrypt_ct_low_width(eng_lw, state, tweak, key, anc, hlp)
            elif var == "dialga256":
                key = [eng_lw.allocate_qureg(8) for _ in range(32)]
                tweak0 = [eng_lw.allocate_qureg(8) for _ in range(16)]
                tweak1 = [eng_lw.allocate_qureg(8) for _ in range(16)]
                state = [eng_lw.allocate_qureg(8) for _ in range(16)]
                anc = eng_lw.allocate_qubit()[0]
                hlp = eng_lw.allocate_qubit()[0]
                dialga256_encrypt_ct_low_width(eng_lw, state, tweak0, tweak1, key, anc, hlp)
            elif var == "dialga256d":
                key = [eng_lw.allocate_qureg(8) for _ in range(32)]
                tweak0 = [eng_lw.allocate_qureg(8) for _ in range(16)]
                tweak1 = [eng_lw.allocate_qureg(8) for _ in range(16)]
                state = [eng_lw.allocate_qureg(8) for _ in range(16)]
                anc = eng_lw.allocate_qubit()[0]
                hlp = eng_lw.allocate_qubit()[0]
                dialga256d_encrypt_ct_low_width(eng_lw, state, tweak0, tweak1, key, anc, hlp)
                
            eng_lw.flush()
            
            rc_ld = AdvancedResourceCounter()
            rc_ld.is_recording = True
            eng_ld = MainEngine(rc_ld, [])
            
            if var == "dialga128":
                key_ld = [eng_ld.allocate_qureg(8) for _ in range(32)]
                tweak_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                state_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                state_anc = eng_ld.allocate_qureg(32)
                state_hlp = eng_ld.allocate_qureg(32)
                tweak_anc = eng_ld.allocate_qureg(32)
                tweak_hlp = eng_ld.allocate_qureg(32)
                dialga128_encrypt_ct_low_depth(eng_ld, state_ld, tweak_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            elif var == "dialga128d":
                key_ld = [eng_ld.allocate_qureg(8) for _ in range(32)]
                tweak_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                state_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                state_anc = eng_ld.allocate_qureg(32)
                state_hlp = eng_ld.allocate_qureg(32)
                tweak_anc = eng_ld.allocate_qureg(32)
                tweak_hlp = eng_ld.allocate_qureg(32)
                dialga128d_encrypt_ct_low_depth(eng_ld, state_ld, tweak_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            elif var == "dialga256":
                key_ld = [eng_ld.allocate_qureg(8) for _ in range(32)]
                tweak0_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                tweak1_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                state_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                state_anc = eng_ld.allocate_qureg(32)
                state_hlp = eng_ld.allocate_qureg(32)
                tweak_anc = eng_ld.allocate_qureg(32)
                tweak_hlp = eng_ld.allocate_qureg(32)
                dialga256_encrypt_ct_low_depth(eng_ld, state_ld, tweak0_ld, tweak1_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
            elif var == "dialga256d":
                key_ld = [eng_ld.allocate_qureg(8) for _ in range(32)]
                tweak0_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                tweak1_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                state_ld = [eng_ld.allocate_qureg(8) for _ in range(16)]
                state_anc = eng_ld.allocate_qureg(32)
                state_hlp = eng_ld.allocate_qureg(32)
                tweak_anc = eng_ld.allocate_qureg(32)
                tweak_hlp = eng_ld.allocate_qureg(32)
                dialga256d_encrypt_ct_low_depth(eng_ld, state_ld, tweak0_ld, tweak1_ld, key_ld, state_anc, state_hlp, tweak_anc, tweak_hlp)
                
            eng_ld.flush()
            
            # Verify identical logical gate operations
            self.assertEqual(rc_lw.get_cnot(), rc_ld.get_cnot(), f"{var} CNOT count invariant failed")
            self.assertEqual(rc_lw.get_1q_clifford(), rc_ld.get_1q_clifford(), f"{var} 1qClifford count invariant failed")
            self.assertEqual(rc_lw.get_t(), rc_ld.get_t(), f"{var} T count invariant failed")
            self.assertEqual(rc_lw.get_measure(), rc_ld.get_measure(), f"{var} Measure count invariant failed")
            
            # Verify depth and width trade-offs
            self.assertLess(rc_ld.t_depth, rc_lw.t_depth, f"{var} T-depth LD < LW invariant failed")
            self.assertLess(rc_ld.pure_circuit_depth, rc_lw.pure_circuit_depth, f"{var} FullDepth LD < LW invariant failed")
            self.assertGreater(rc_ld.max_width, rc_lw.max_width, f"{var} Width LD > LW invariant failed")
            
            # Check expected widths
            self.assertEqual(rc_lw.max_width, expected_widths[f"{var}_lw"], f"{var}_lw width mismatch")
            self.assertEqual(rc_ld.max_width, expected_widths[f"{var}_ld"], f"{var}_ld width mismatch")

    def test_resource_regression(self):
        """Exact resource regression tests for all 8 combinations"""
        from test_advanced_counter import AdvancedResourceCounter
        
        expected_resources = {
            "dialga128_lw": {"cnot": 59392, "oneq": 34979, "t": 15872, "measure": 3968, "t_depth": 3968, "full_depth": 47632, "width": 514},
            "dialga128_ld": {"cnot": 59392, "oneq": 34979, "t": 15872, "measure": 3968, "t_depth": 88, "full_depth": 1212, "width": 640},
            "dialga128d_lw": {"cnot": 47872, "oneq": 28197, "t": 12800, "measure": 3200, "t_depth": 3200, "full_depth": 38404, "width": 514},
            "dialga128d_ld": {"cnot": 47872, "oneq": 28197, "t": 12800, "measure": 3200, "t_depth": 72, "full_depth": 991, "width": 640},
            "dialga256_lw": {"cnot": 75264, "oneq": 43683, "t": 19968, "measure": 4992, "t_depth": 4992, "full_depth": 59909, "width": 642},
            "dialga256_ld": {"cnot": 75264, "oneq": 43683, "t": 19968, "measure": 4992, "t_depth": 88, "full_depth": 1213, "width": 768},
            "dialga256d_lw": {"cnot": 59904, "oneq": 34725, "t": 15872, "measure": 3968, "t_depth": 3968, "full_depth": 47621, "width": 642},
            "dialga256d_ld": {"cnot": 59904, "oneq": 34725, "t": 15872, "measure": 3968, "t_depth": 72, "full_depth": 992, "width": 768}
        }
        
        for key_name, expected in expected_resources.items():
            rc = AdvancedResourceCounter()
            rc.is_recording = True
            eng = MainEngine(rc, [])
            
            if "dialga128_" in key_name:
                key = [eng.allocate_qureg(8) for _ in range(32)]
                tweak = [eng.allocate_qureg(8) for _ in range(16)]
                state = [eng.allocate_qureg(8) for _ in range(16)]
                if "_lw" in key_name:
                    anc = eng.allocate_qubit()[0]
                    hlp = eng.allocate_qubit()[0]
                    dialga128_encrypt_ct_low_width(eng, state, tweak, key, anc, hlp)
                else:
                    state_anc = eng.allocate_qureg(32)
                    state_hlp = eng.allocate_qureg(32)
                    tweak_anc = eng.allocate_qureg(32)
                    tweak_hlp = eng.allocate_qureg(32)
                    dialga128_encrypt_ct_low_depth(eng, state, tweak, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            elif "dialga128d_" in key_name:
                key = [eng.allocate_qureg(8) for _ in range(32)]
                tweak = [eng.allocate_qureg(8) for _ in range(16)]
                state = [eng.allocate_qureg(8) for _ in range(16)]
                if "_lw" in key_name:
                    anc = eng.allocate_qubit()[0]
                    hlp = eng.allocate_qubit()[0]
                    dialga128d_encrypt_ct_low_width(eng, state, tweak, key, anc, hlp)
                else:
                    state_anc = eng.allocate_qureg(32)
                    state_hlp = eng.allocate_qureg(32)
                    tweak_anc = eng.allocate_qureg(32)
                    tweak_hlp = eng.allocate_qureg(32)
                    dialga128d_encrypt_ct_low_depth(eng, state, tweak, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            elif "dialga256_" in key_name:
                key = [eng.allocate_qureg(8) for _ in range(32)]
                tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
                tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
                state = [eng.allocate_qureg(8) for _ in range(16)]
                if "_lw" in key_name:
                    anc = eng.allocate_qubit()[0]
                    hlp = eng.allocate_qubit()[0]
                    dialga256_encrypt_ct_low_width(eng, state, tweak0, tweak1, key, anc, hlp)
                else:
                    state_anc = eng.allocate_qureg(32)
                    state_hlp = eng.allocate_qureg(32)
                    tweak_anc = eng.allocate_qureg(32)
                    tweak_hlp = eng.allocate_qureg(32)
                    dialga256_encrypt_ct_low_depth(eng, state, tweak0, tweak1, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
            elif "dialga256d_" in key_name:
                key = [eng.allocate_qureg(8) for _ in range(32)]
                tweak0 = [eng.allocate_qureg(8) for _ in range(16)]
                tweak1 = [eng.allocate_qureg(8) for _ in range(16)]
                state = [eng.allocate_qureg(8) for _ in range(16)]
                if "_lw" in key_name:
                    anc = eng.allocate_qubit()[0]
                    hlp = eng.allocate_qubit()[0]
                    dialga256d_encrypt_ct_low_width(eng, state, tweak0, tweak1, key, anc, hlp)
                else:
                    state_anc = eng.allocate_qureg(32)
                    state_hlp = eng.allocate_qureg(32)
                    tweak_anc = eng.allocate_qureg(32)
                    tweak_hlp = eng.allocate_qureg(32)
                    dialga256d_encrypt_ct_low_depth(eng, state, tweak0, tweak1, key, state_anc, state_hlp, tweak_anc, tweak_hlp)
                    
            eng.flush()
            
            self.assertEqual(rc.get_cnot(), expected["cnot"], f"{key_name} CNOT mismatch")
            self.assertEqual(rc.get_1q_clifford(), expected["oneq"], f"{key_name} 1qClifford mismatch")
            self.assertEqual(rc.get_t(), expected["t"], f"{key_name} T count mismatch")
            self.assertEqual(rc.get_measure(), expected["measure"], f"{key_name} Measure mismatch")
            self.assertEqual(rc.t_depth, expected["t_depth"], f"{key_name} T-depth mismatch")
            self.assertEqual(rc.pure_circuit_depth, expected["full_depth"], f"{key_name} FullDepth mismatch")
            self.assertEqual(rc.max_width, expected["width"], f"{key_name} Width mismatch")

    # ========================================================
    # PART L: ADVANCED DEPTH COUNTER UNIT TESTS
    # ========================================================
    def test_counter_semantics(self):
        """Resource-counter validation"""
        from test_advanced_counter import AdvancedResourceCounter
        from projectq.ops import T, Tdag, CNOT, H, Measure
        from corefunctions import and_gate, and_gate_dagger, sbox_8bitCT
        
        # Scenario 1 (Test A): two T gates on disjoint qubits -> T-depth 1
        rc = AdvancedResourceCounter()
        rc.is_recording = True
        eng = MainEngine(rc, [])
        q0 = eng.allocate_qubit()[0]
        q1 = eng.allocate_qubit()[0]
        T | q0
        T | q1
        eng.flush()
        self.assertEqual(rc.t_depth, 1)
        
        # Scenario 2 (Test B): two sequential T gates on same qubit -> T-depth 2
        rc = AdvancedResourceCounter()
        rc.is_recording = True
        eng = MainEngine(rc, [])
        q0 = eng.allocate_qubit()[0]
        T | q0
        Tdag | q0
        eng.flush()
        self.assertEqual(rc.t_depth, 2)
        
        # Scenario 3 (Test C): T -- CNOT -- T dependency -> T-depth 2
        rc = AdvancedResourceCounter()
        rc.is_recording = True
        eng = MainEngine(rc, [])
        q0 = eng.allocate_qubit()[0]
        q1 = eng.allocate_qubit()[0]
        T | q0
        CNOT | (q0, q1)
        T | q1
        eng.flush()
        self.assertEqual(rc.t_depth, 2)

        # Counter test D: independent CNOTs -> same FullDepth layer (1)
        rc = AdvancedResourceCounter()
        rc.is_recording = True
        eng = MainEngine(rc, [])
        q0 = eng.allocate_qubit()[0]
        q1 = eng.allocate_qubit()[0]
        q2 = eng.allocate_qubit()[0]
        q3 = eng.allocate_qubit()[0]
        CNOT | (q0, q1)
        CNOT | (q2, q3)
        eng.flush()
        self.assertEqual(rc.pure_circuit_depth, 1)

        # Counter test E: dependent CNOTs -> FullDepth increases (2)
        rc = AdvancedResourceCounter()
        rc.is_recording = True
        eng = MainEngine(rc, [])
        q0 = eng.allocate_qubit()[0]
        q1 = eng.allocate_qubit()[0]
        q2 = eng.allocate_qubit()[0]
        CNOT | (q0, q1)
        CNOT | (q1, q2)
        eng.flush()
        self.assertEqual(rc.pure_circuit_depth, 2)

        # Counter test F: measurement increases pure_circuit_depth
        rc = AdvancedResourceCounter()
        rc.is_recording = True
        eng = MainEngine(rc, [])
        q0 = eng.allocate_qubit()[0]
        H | q0
        Measure | q0
        eng.flush()
        self.assertEqual(rc.pure_circuit_depth, 2)

        # Counter test G: independent ANDs (disjoint qubits) -> T-depth equal to one AND
        rc1 = AdvancedResourceCounter()
        rc1.is_recording = True
        eng1 = MainEngine(rc1, [])
        a1 = eng1.allocate_qubit()[0]
        b1 = eng1.allocate_qubit()[0]
        anc1 = eng1.allocate_qubit()[0]
        hlp1 = eng1.allocate_qubit()[0]
        and_gate(a1, b1, anc1, hlp1)
        eng1.flush()
        t_depth_1and = rc1.t_depth
        All(Measure) | [a1, b1, anc1, hlp1]
        eng1.flush()
        
        rc2 = AdvancedResourceCounter()
        rc2.is_recording = True
        eng2 = MainEngine(rc2, [])
        a2 = eng2.allocate_qubit()[0]
        b2 = eng2.allocate_qubit()[0]
        anc2 = eng2.allocate_qubit()[0]
        hlp2 = eng2.allocate_qubit()[0]
        a3 = eng2.allocate_qubit()[0]
        b3 = eng2.allocate_qubit()[0]
        anc3 = eng2.allocate_qubit()[0]
        hlp3 = eng2.allocate_qubit()[0]
        and_gate(a2, b2, anc2, hlp2)
        and_gate(a3, b3, anc3, hlp3)
        eng2.flush()
        t_depth_2and = rc2.t_depth
        All(Measure) | [a2, b2, anc2, hlp2, a3, b3, anc3, hlp3]
        eng2.flush()
        self.assertEqual(t_depth_1and, t_depth_2and)

        # Counter test H: shared-workspace ANDs (reused workspace) -> dependency depth increases
        rc3 = AdvancedResourceCounter()
        rc3.is_recording = True
        eng3 = MainEngine(rc3, [])
        a4 = eng3.allocate_qubit()[0]
        b4 = eng3.allocate_qubit()[0]
        a5 = eng3.allocate_qubit()[0]
        b5 = eng3.allocate_qubit()[0]
        anc_shared = eng3.allocate_qubit()[0]
        hlp_shared = eng3.allocate_qubit()[0]
        and_gate(a4, b4, anc_shared, hlp_shared)
        and_gate_dagger(a4, b4, anc_shared)
        t_depth_roundtrip = rc3.t_depth
        and_gate(a5, b5, anc_shared, hlp_shared)
        eng3.flush()
        t_depth_combined = rc3.t_depth
        self.assertGreater(t_depth_combined, t_depth_roundtrip)
        All(Measure) | [a4, b4, a5, b5, anc_shared, hlp_shared]
        eng3.flush()

        # Counter test I: LW versus LD 8-bit S-box resources
        rc_lw = AdvancedResourceCounter()
        rc_lw.is_recording = True
        eng_lw = MainEngine(rc_lw, [])
        cell_lw = eng_lw.allocate_qureg(8)
        anc_lw = eng_lw.allocate_qubit()[0]
        hlp_lw = eng_lw.allocate_qubit()[0]
        sbox_8bitCT(eng_lw, cell_lw, 0, anc_lw, hlp_lw)
        eng_lw.flush()
        
        rc_ld = AdvancedResourceCounter()
        rc_ld.is_recording = True
        eng_ld = MainEngine(rc_ld, [])
        cell_ld = eng_ld.allocate_qureg(8)
        anc0_ld = eng_ld.allocate_qubit()[0]
        hlp0_ld = eng_ld.allocate_qubit()[0]
        anc1_ld = eng_ld.allocate_qubit()[0]
        hlp1_ld = eng_ld.allocate_qubit()[0]
        sbox_8bitCT(eng_ld, cell_ld, 0, anc0_ld, hlp0_ld, anc1_ld, hlp1_ld)
        eng_ld.flush()
        
        self.assertEqual(rc_lw.get_cnot(), rc_ld.get_cnot())
        self.assertEqual(rc_lw.get_1q_clifford(), rc_ld.get_1q_clifford())
        self.assertEqual(rc_lw.get_t(), rc_ld.get_t())
        self.assertEqual(rc_lw.get_measure(), rc_ld.get_measure())
        
        self.assertLess(rc_ld.t_depth, rc_lw.t_depth)
        self.assertLess(rc_ld.pure_circuit_depth, rc_lw.pure_circuit_depth)
        self.assertGreater(rc_ld.max_width, rc_lw.max_width)
        
        All(Measure) | list(cell_lw) + [anc_lw, hlp_lw]
        eng_lw.flush()
        All(Measure) | list(cell_ld) + [anc0_ld, hlp0_ld, anc1_ld, hlp1_ld]
        eng_ld.flush()


if __name__ == "__main__":
    unittest.main()
