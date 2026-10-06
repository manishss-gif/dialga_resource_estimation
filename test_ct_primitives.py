import unittest
import projectq.backends._sim._pysim
projectq.backends._sim._pysim._USE_REFCHECK = False

from projectq import MainEngine
from projectq.backends import Simulator
from projectq.ops import X, Measure, All, H
from corefunctions import (
    and_gate,
    and_gate_dagger,
    SubCellWithAncillaCT_low_width,
    sbox_8bitCT
)
from classical_dialga_reference import classical_sbox_8bit

# Expected Sb0 table values
Sb0 = [0xC, 0xA, 0xD, 0x3, 0xE, 0xB, 0xF, 0x7,
       0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6]


class TestCTPrimitives(unittest.TestCase):
    def test_ct_logical_and_forward(self):
        """CT logical AND forward"""
        sim = Simulator()
        eng = MainEngine(backend=sim, engine_list=[])
        
        a = eng.allocate_qubit()[0]
        b = eng.allocate_qubit()[0]
        ancilla = eng.allocate_qubit()[0]
        helper = eng.allocate_qubit()[0]
        
        for a_val in [0, 1]:
            for b_val in [0, 1]:
                if a_val:
                    X | a
                if b_val:
                    X | b
                eng.flush()
                
                # Apply logical AND
                and_gate(a, b, ancilla, helper)
                eng.flush()
                
                # Measure all to check correctness of forward gate
                All(Measure) | [a, b, ancilla, helper]
                eng.flush()
                
                self.assertEqual(int(a), a_val)
                self.assertEqual(int(b), b_val)
                self.assertEqual(int(ancilla), a_val & b_val)
                self.assertEqual(int(helper), 0)
                
                # Reset
                for q in [a, b, ancilla, helper]:
                    if int(q) == 1:
                        X | q
                eng.flush()
                
        # Cleanup final deallocate
        All(Measure) | [a, b, ancilla, helper]
        eng.flush()

    def test_ct_logical_and_roundtrip(self):
        """CT logical AND roundtrip"""
        sim = Simulator()
        eng = MainEngine(backend=sim, engine_list=[])
        
        a = eng.allocate_qubit()[0]
        b = eng.allocate_qubit()[0]
        ancilla = eng.allocate_qubit()[0]
        helper = eng.allocate_qubit()[0]
        
        for a_val in [0, 1]:
            for b_val in [0, 1]:
                # Repeat to exercise stochastic measurement-assisted uncomputation.
                for trial in range(10):
                    if a_val:
                        X | a
                    if b_val:
                        X | b
                    eng.flush()
                    
                    # Apply logical AND and uncomputation
                    and_gate(a, b, ancilla, helper)
                    and_gate_dagger(a, b, ancilla)
                    
                    # Measure all to check cleanup
                    All(Measure) | [a, b, ancilla, helper]
                    eng.flush()
                    
                    self.assertEqual(int(ancilla), 0)
                    self.assertEqual(int(helper), 0)
                    self.assertEqual(int(a), a_val)
                    self.assertEqual(int(b), b_val)
                    
                    # Reset
                    for q in [a, b, ancilla, helper]:
                        if int(q) == 1:
                            X | q
                    eng.flush()
                    
        # Cleanup final deallocate
        All(Measure) | [a, b, ancilla, helper]
        eng.flush()

    def test_ct_and_roundtrip_superposition(self):
        """Phase-sensitive superposition AND/ANDdag roundtrip test"""
        sim = Simulator()
        eng = MainEngine(backend=sim, engine_list=[])
        
        a = eng.allocate_qubit()[0]
        b = eng.allocate_qubit()[0]
        ancilla = eng.allocate_qubit()[0]
        helper = eng.allocate_qubit()[0]
        
        # Prepare controls in superposition: |+> |+> |0> |0>
        H | a
        H | b
        eng.flush()
        
        # Get mapping and state vector after preparing superposition
        mapping_init, state_init = sim.cheat()
        
        # Apply AND and AND-dagger
        and_gate(a, b, ancilla, helper)
        and_gate_dagger(a, b, ancilla)
        eng.flush()
        
        # Get mapping and state vector after roundtrip
        mapping_final, state_final = sim.cheat()
        
        self.assertEqual(mapping_init, mapping_final)
        
        # Find global phase factor by comparing a non-zero amplitude
        phase_factor = None
        for i in range(len(state_init)):
            if abs(state_init[i]) > 1e-5:
                phase_factor = state_final[i] / state_init[i]
                break
                
        self.assertIsNotNone(phase_factor)
        self.assertAlmostEqual(abs(phase_factor), 1.0, places=5)
        
        # Verify state_final == phase_factor * state_init
        for i in range(len(state_init)):
            diff = state_final[i] - phase_factor * state_init[i]
            self.assertAlmostEqual(abs(diff), 0.0, places=5)
            
        # Clean/reset
        All(Measure) | [a, b, ancilla, helper]
        eng.flush()

    def test_ct_4bit_sbox(self):
        """CT 4-bit S-box"""
        sim = Simulator()
        eng = MainEngine(backend=sim, engine_list=[])
        
        cell = eng.allocate_qureg(4)
        ancilla = eng.allocate_qubit()[0]
        helper = eng.allocate_qubit()[0]
        
        for x in range(16):
            bin_str = bin(x)[2:].zfill(4)
            for b_idx in range(4):
                if bin_str[b_idx] == '1':
                    X | cell[b_idx]
            eng.flush()
            
            # Run 4-bit Sbox
            c_out = SubCellWithAncillaCT_low_width(eng, cell[0], cell[1], cell[2], cell[3], ancilla, helper)
            cell[0], cell[1], cell[2], cell[3] = c_out
            
            # Measure all at the end
            All(Measure) | list(cell) + [ancilla, helper]
            eng.flush()
            
            self.assertEqual(int(ancilla), 0)
            self.assertEqual(int(helper), 0)
            
            # Read outputs
            val_out = 0
            for i in range(4):
                val_out = (val_out << 1) | int(cell[i])
            self.assertEqual(val_out, Sb0[x])
            
            # Reset qubits to 0 state
            for q in list(cell) + [ancilla, helper]:
                if int(q) == 1:
                    X | q
            eng.flush()
            
        All(Measure) | list(cell) + [ancilla, helper]
        eng.flush()

    def test_ct_4bit_sbox_superposition(self):
        """Phase-sensitive superposition 4-bit S-box test"""
        sim = Simulator()
        eng = MainEngine(backend=sim, engine_list=[])
        
        cell = eng.allocate_qureg(4)
        ancilla = eng.allocate_qubit()[0]
        helper = eng.allocate_qubit()[0]
        
        # Prepare uniform superposition of 4 input qubits: (1/4) sum |x> |0> |0>
        for q in cell:
            H | q
        eng.flush()
        
        # Get mapping and state vector
        mapping_init, state_init = sim.cheat()
        
        # Apply S-box
        c_out = SubCellWithAncillaCT_low_width(eng, cell[0], cell[1], cell[2], cell[3], ancilla, helper)
        cell[0], cell[1], cell[2], cell[3] = c_out
        eng.flush()
        
        mapping_final, state_final = sim.cheat()
        self.assertEqual(mapping_init, mapping_final)
        
        pos_cell = [mapping_init[q.id] for q in cell]
        pos_anc = mapping_init[ancilla.id]
        pos_hlp = mapping_init[helper.id]
        
        # Reconstruct state values from state vectors
        init_amps = {}
        for idx, amp in enumerate(state_init):
            if abs(amp) > 1e-5:
                val_cell = 0
                for b_idx, pos in enumerate(pos_cell):
                    bit = (idx >> pos) & 1
                    val_cell = (val_cell << 1) | bit
                val_anc = (idx >> pos_anc) & 1
                val_hlp = (idx >> pos_hlp) & 1
                self.assertEqual(val_anc, 0)
                self.assertEqual(val_hlp, 0)
                init_amps[val_cell] = amp
                
        final_amps = {}
        for idx, amp in enumerate(state_final):
            if abs(amp) > 1e-5:
                val_cell = 0
                for b_idx, pos in enumerate(pos_cell):
                    bit = (idx >> pos) & 1
                    val_cell = (val_cell << 1) | bit
                val_anc = (idx >> pos_anc) & 1
                val_hlp = (idx >> pos_hlp) & 1
                self.assertEqual(val_anc, 0, f"Ancilla not clean at output state {val_cell}")
                self.assertEqual(val_hlp, 0, f"Helper not clean at output state {val_cell}")
                final_amps[val_cell] = amp
                
        self.assertEqual(len(init_amps), 16)
        self.assertEqual(len(final_amps), 16)
        
        # Find global phase factor using x = 0
        phase_factor = final_amps[Sb0[0]] / init_amps[0]
        self.assertAlmostEqual(abs(phase_factor), 1.0, places=5)
        
        for x in range(16):
            diff = final_amps[Sb0[x]] - phase_factor * init_amps[x]
            self.assertAlmostEqual(abs(diff), 0.0, places=5, msg=f"Phase mismatch at x={x}")
            
        # Clean/reset
        All(Measure) | list(cell) + [ancilla, helper]
        eng.flush()

    def test_ct_8bit_sbox_lw(self):
        """CT 8-bit S-box LW"""
        sim = Simulator()
        eng = MainEngine(backend=sim, engine_list=[])
        
        cell = eng.allocate_qureg(8)
        ancilla = eng.allocate_qubit()[0]
        helper = eng.allocate_qubit()[0]
        
        for cell_idx in range(4):
            for x in range(256):
                bin_str = bin(x)[2:].zfill(8)
                for b_idx in range(8):
                    if bin_str[b_idx] == '1':
                        X | cell[b_idx]
                eng.flush()
                
                # Apply 8-bit Sbox in sequential LW mode
                sbox_8bitCT(eng, cell, cell_idx, ancilla, helper)
                
                # Measure all at the end
                All(Measure) | list(cell) + [ancilla, helper]
                eng.flush()
                
                self.assertEqual(int(ancilla), 0)
                self.assertEqual(int(helper), 0)
                
                # Read outputs
                val_out = 0
                for i in range(8):
                    val_out = (val_out << 1) | int(cell[i])
                expected = classical_sbox_8bit(x, cell_idx)
                self.assertEqual(val_out, expected)
                
                # Reset qubits to 0 state
                for q in list(cell) + [ancilla, helper]:
                    if int(q) == 1:
                        X | q
                eng.flush()
                
        All(Measure) | list(cell) + [ancilla, helper]
        eng.flush()

    def test_ct_8bit_sbox_ld(self):
        """CT 8-bit S-box LD"""
        sim = Simulator()
        eng = MainEngine(backend=sim, engine_list=[])
        
        cell = eng.allocate_qureg(8)
        anc0 = eng.allocate_qubit()[0]
        hlp0 = eng.allocate_qubit()[0]
        anc1 = eng.allocate_qubit()[0]
        hlp1 = eng.allocate_qubit()[0]
        
        for cell_idx in range(4):
            for x in range(256):
                bin_str = bin(x)[2:].zfill(8)
                for b_idx in range(8):
                    if bin_str[b_idx] == '1':
                        X | cell[b_idx]
                eng.flush()
                
                # Apply 8-bit Sbox in parallel LD mode
                sbox_8bitCT(eng, cell, cell_idx, anc0, hlp0, anc1, hlp1)
                
                # Measure all at the end
                All(Measure) | list(cell) + [anc0, hlp0, anc1, hlp1]
                eng.flush()
                
                self.assertEqual(int(anc0), 0)
                self.assertEqual(int(hlp0), 0)
                self.assertEqual(int(anc1), 0)
                self.assertEqual(int(hlp1), 0)
                
                # Read outputs
                val_out = 0
                for i in range(8):
                    val_out = (val_out << 1) | int(cell[i])
                expected = classical_sbox_8bit(x, cell_idx)
                self.assertEqual(val_out, expected)
                
                # Reset qubits to 0 state
                for q in list(cell) + [anc0, hlp0, anc1, hlp1]:
                    if int(q) == 1:
                        X | q
                eng.flush()
                
        All(Measure) | list(cell) + [anc0, hlp0, anc1, hlp1]
        eng.flush()


if __name__ == "__main__":
    unittest.main()

