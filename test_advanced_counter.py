# test_advanced_counter.py
# Generic AdvancedResourceCounter class to track T-depth and pure circuit depth.

from projectq.backends import ResourceCounter
from projectq.ops import TGate, Allocate, Deallocate, DaggeredGate
from projectq.meta import get_control_count

class AdvancedResourceCounter(ResourceCounter):
    def __init__(self):
        super().__init__()
        # key: qubit id, T-depth of this qubit
        self._t_depth_of_qubit = {}
        self._previous_max_t_depth = 0
        # key: qubit id, Toffoli-depth of this qubit
        self._toffoli_depth_of_qubit = {}
        self._previous_max_toffoli_depth = 0
        # key: qubit id, pure circuit depth of this qubit (excluding Allocate/Deallocate/Measure)
        self._circuit_depth_of_qubit = {}
        self._previous_max_circuit_depth = 0
        # key: qubit id, DAG depth of this qubit (including Allocate/Deallocate/Measure)
        self._dag_depth_of_qubit = {}
        self._previous_max_dag_depth = 0
        
        # Flag to control recording
        self.is_recording = False
        # Custom gate counts for the recorded portion
        self.recorded_gate_counts = {}

    @property
    def t_depth(self):
        if self._t_depth_of_qubit:
            current_max = max(self._t_depth_of_qubit.values())
            return max(current_max, self._previous_max_t_depth)
        return self._previous_max_t_depth

    @property
    def toffoli_depth(self):
        if self._toffoli_depth_of_qubit:
            current_max = max(self._toffoli_depth_of_qubit.values())
            return max(current_max, self._previous_max_toffoli_depth)
        return self._previous_max_toffoli_depth

    @property
    def pure_circuit_depth(self):
        if self._circuit_depth_of_qubit:
            current_max = max(self._circuit_depth_of_qubit.values())
            return max(current_max, self._previous_max_circuit_depth)
        return self._previous_max_circuit_depth

    @property
    def recorded_dag_depth(self):
        if self._dag_depth_of_qubit:
            current_max = max(self._dag_depth_of_qubit.values())
            return max(current_max, self._previous_max_dag_depth)
        return self._previous_max_dag_depth

    def _add_cmd(self, cmd):
        # Base count (always updates self.gate_counts)
        super()._add_cmd(cmd)

        gate = cmd.gate
        
        # If we are not recording, do not count/process
        if not self.is_recording:
            return

        # Record in our custom gate counts
        ctrl_cnt = get_control_count(cmd)
        gate_key = (gate, ctrl_cnt)
        self.recorded_gate_counts[gate_key] = self.recorded_gate_counts.get(gate_key, 0) + 1

        is_t_or_tdag = False
        if isinstance(gate, DaggeredGate):
            actual_gate = gate._gate
        else:
            actual_gate = gate
            
        if isinstance(actual_gate, TGate):
            is_t_or_tdag = True

        if cmd.gate == Allocate:
            qubit_id = cmd.qubits[0][0].id
            self._t_depth_of_qubit[qubit_id] = 0
            self._toffoli_depth_of_qubit[qubit_id] = 0
            self._circuit_depth_of_qubit[qubit_id] = 0
            self._dag_depth_of_qubit[qubit_id] = 1
        elif cmd.gate == Deallocate:
            qubit_id = cmd.qubits[0][0].id
            t_depth = self._t_depth_of_qubit.pop(qubit_id, 0)
            self._previous_max_t_depth = max(self._previous_max_t_depth, t_depth)
            
            toff_depth = self._toffoli_depth_of_qubit.pop(qubit_id, 0)
            self._previous_max_toffoli_depth = max(self._previous_max_toffoli_depth, toff_depth)
            
            c_depth = self._circuit_depth_of_qubit.pop(qubit_id, 0)
            self._previous_max_circuit_depth = max(self._previous_max_circuit_depth, c_depth)
            
            d_depth = self._dag_depth_of_qubit.pop(qubit_id, 0) + 1
            self._previous_max_dag_depth = max(self._previous_max_dag_depth, d_depth)
        else:
            qubit_ids = set()
            for qureg in cmd.all_qubits:
                for qubit in qureg:
                    qubit_ids.add(qubit.id)
            if qubit_ids:
                # 1. Update T-depth
                max_t = max(self._t_depth_of_qubit.get(qubit_id, 0) for qubit_id in qubit_ids)
                if is_t_or_tdag:
                    for qubit_id in qubit_ids:
                        self._t_depth_of_qubit[qubit_id] = max_t + 1
                else:
                    for qubit_id in qubit_ids:
                        self._t_depth_of_qubit[qubit_id] = max_t

                # 1.5 Update Toffoli-depth
                max_toff = max(self._toffoli_depth_of_qubit.get(qubit_id, 0) for qubit_id in qubit_ids)
                gate_name = self.format_gate_name(gate_key)
                if gate_name == 'Toffoli':
                    for qubit_id in qubit_ids:
                        self._toffoli_depth_of_qubit[qubit_id] = max_toff + 1
                else:
                    for qubit_id in qubit_ids:
                        self._toffoli_depth_of_qubit[qubit_id] = max_toff

                # 2. Update Pure Circuit Depth (includes Measure since is_recording is True)
                max_c = max(self._circuit_depth_of_qubit.get(qubit_id, 0) for qubit_id in qubit_ids)
                for qubit_id in qubit_ids:
                    self._circuit_depth_of_qubit[qubit_id] = max_c + 1
                    
                # 3. Update DAG Depth
                max_d = max(self._dag_depth_of_qubit.get(qubit_id, 0) for qubit_id in qubit_ids)
                for qubit_id in qubit_ids:
                    self._dag_depth_of_qubit[qubit_id] = max_d + 1

    def format_gate_name(self, gate_key):
        gate, ctrl_cnt = gate_key
        
        from projectq.ops import DaggeredGate
        is_dagger = False
        if isinstance(gate, DaggeredGate):
            is_dagger = True
            gate = gate._gate
            
        class_name = gate.__class__.__name__
        name = class_name
        
        if class_name == 'AllocateQubitGate':
            name = 'Allocate'
        elif class_name == 'DeallocateQubitGate':
            name = 'Deallocate'
        elif class_name == 'MeasureGate':
            name = 'Measure'
        elif class_name == 'XGate':
            if ctrl_cnt == 0:
                name = 'X'
            elif ctrl_cnt == 1:
                name = 'CNOT'
            elif ctrl_cnt == 2:
                name = 'Toffoli'
            else:
                name = f'MCX ({ctrl_cnt} controls)'
        elif class_name == 'HGate':
            name = 'H'
        elif class_name == 'TGate':
            name = 'T'
        elif class_name == 'SGate':
            name = 'S'
        elif class_name == 'SwapGate':
            name = 'Swap'
        elif class_name == 'ZGate':
            name = 'Z'
            
        if is_dagger:
            name += r'^\dagger'
            
        return name

    def get_cnot(self):
        cnot_count = 0
        for gate_key, count in self.recorded_gate_counts.items():
            name = self.format_gate_name(gate_key)
            if name == 'CNOT' or 'MCX' in name or 'C^' in name:
                cnot_count += count
        return cnot_count

    def get_1q_clifford(self):
        cliff_count = 0
        for gate_key, count in self.recorded_gate_counts.items():
            name = self.format_gate_name(gate_key)
            if name in ('H', 'S', 'S^\\dagger', 'X', 'Z'):
                cliff_count += count
        return cliff_count

    def get_t(self):
        t_count = 0
        for gate_key, count in self.recorded_gate_counts.items():
            name = self.format_gate_name(gate_key)
            if name in ('T', 'T^\\dagger'):
                t_count += count
        return t_count

    def get_toffoli(self):
        toffoli_count = 0
        for gate_key, count in self.recorded_gate_counts.items():
            name = self.format_gate_name(gate_key)
            if name == 'Toffoli':
                toffoli_count += count
        return toffoli_count
        
    def get_measure(self):
        measure_count = 0
        for gate_key, count in self.recorded_gate_counts.items():
            name = self.format_gate_name(gate_key)
            if 'Measure' in name:
                measure_count += count
        return measure_count

    def get_t_depth(self):
        return self.t_depth

    def get_toffoli_depth(self):
        return self.toffoli_depth

    def get_full_depth(self):
        return self.pure_circuit_depth

    def get_width(self):
        return self.max_width


METRICS = (
    ("cnot", "#CNOT"),
    ("oneq", "#1qClifford"),
    ("t", "#T"),
    ("measure", "#Measure"),
    ("t_depth", "T-depth"),
    ("full_depth", "Full depth"),
    ("width", "Width"),
)


def count_resources(build):
    """Runs build(eng) on a fresh counter and returns the seven reported metrics."""
    from projectq import MainEngine

    rc = AdvancedResourceCounter()
    rc.is_recording = True
    eng = MainEngine(backend=rc, engine_list=[])
    build(eng)
    eng.flush()
    rc.is_recording = False
    return {
        "cnot": rc.get_cnot(),
        "oneq": rc.get_1q_clifford(),
        "t": rc.get_t(),
        "measure": rc.get_measure(),
        "t_depth": rc.t_depth,
        "full_depth": rc.pure_circuit_depth,
        "width": rc.max_width,
    }


def as_tuple(res):
    return tuple(res[key] for key, _ in METRICS)


def print_resources(title, res):
    print(title)
    for key, label in METRICS:
        print(f"    {label:<12} {res[key]}")
