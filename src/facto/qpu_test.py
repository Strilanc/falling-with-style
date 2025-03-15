import collections
import math
import random
from cmath import phase

import pytest
from matplotlib import pyplot as plt

from facto.qpu import QPU
from facto.main_make_plot import generate_qasm_for_shor_circuit, simulate_sampling_frequency_of_periodic_signal


@pytest.mark.parametrize('n', [1, 2, 5, 10, 20])
def test_iadd_no_carry(n: int):
    for _ in range(20):
        qpu = QPU()
        m = 20
        v1s = []
        v2s = []
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
            v2s.append(qpu.alloc_qubit())
        for q in v1s:
            q.debug_val = random.randrange(2 ** m)
        for q in v2s:
            q.debug_val = random.randrange(2 ** m)
        expected_1s = qpu.debug_read_vals(v1s, m)
        expected_2s = [(b + a) % (2 ** (n)) for a, b in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals(v2s, m))]
        qpu.iadd(offset=v1s, target=v2s)
        assert qpu.debug_read_vals(v1s, m) == expected_1s
        assert qpu.debug_read_vals(v2s, m) == expected_2s


def test_iadd_carry():
    qpu = QPU()
    for _ in range(20):
        n = 10
        m = 20
        v1s = []
        v2s = []
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
            v2s.append(qpu.alloc_qubit())
        v2s.append(qpu.alloc_qubit())
        for q in v1s:
            q.debug_val = random.randrange(2 ** m)
        for q in v2s:
            q.debug_val = random.randrange(2 ** m)
        expected_1s = qpu.debug_read_vals(v1s, m)
        expected_2s = [(b + a) % (2 ** (n + 1)) for a, b in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals(v2s, m))]
        qpu.iadd(offset=v1s, target=v2s)
        assert qpu.debug_read_vals(v1s, m) == expected_1s
        assert qpu.debug_read_vals(v2s, m) == expected_2s


@pytest.mark.parametrize('n', [1, 2, 5, 10, 20])
def test_isub_no_carry(n: int):
    qpu = QPU()
    for _ in range(20):
        m = 20
        v1s = []
        v2s = []
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
            v2s.append(qpu.alloc_qubit())
        for q in v1s:
            q.debug_val = random.randrange(2 ** m)
        for q in v2s:
            q.debug_val = random.randrange(2 ** m)
        expected_1s = qpu.debug_read_vals(v1s, m)
        expected_2s = [(b - a) % (2 ** (n)) for a, b in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals(v2s, m))]
        qpu.isub(offset=v1s, target=v2s)
        assert qpu.debug_read_vals(v1s, m) == expected_1s
        assert qpu.debug_read_vals(v2s, m) == expected_2s


@pytest.mark.parametrize('n', [1, 2, 5, 10, 20])
def test_isub_carry(n: int):
    qpu = QPU()
    for _ in range(20):
        m = 20
        v1s = []
        v2s = []
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
            v2s.append(qpu.alloc_qubit())
        v2s.append(qpu.alloc_qubit())
        for q in v1s:
            q.debug_val = random.randrange(2 ** m)
        for q in v2s:
            q.debug_val = random.randrange(2 ** m)
        expected_1s = qpu.debug_read_vals(v1s, m)
        expected_2s = [(b - a) % (2 ** (n + 1)) for a, b in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals(v2s, m))]
        qpu.isub(offset=v1s, target=v2s)
        assert qpu.debug_read_vals(v1s, m) == expected_1s
        assert qpu.debug_read_vals(v2s, m) == expected_2s


@pytest.mark.parametrize('n', [1, 2, 5, 10, 20])
def test_isub_const(n: int):
    qpu = QPU()
    for _ in range(20):
        m = 20
        v1s = []
        control = qpu.alloc_qubit()
        control.debug_val = random.randrange(2 ** m)
        const = random.randrange(2**n)
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
        for q in v1s:
            q.debug_val = random.randrange(2 ** m)
        expected_1s = [(v - const * c) % (2 ** n) for v, c in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals([control], m))]
        qpu.isub_const(target=v1s, control=control, offset_const=const)
        assert qpu.debug_read_vals(v1s, m) == expected_1s


@pytest.mark.parametrize('n', [1, 2, 5, 10, 20])
def test_flip_if_lt(n: int):
    qpu = QPU()
    for _ in range(200):
        m = 20
        v1s = []
        v2s = []
        gs = []
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
            v2s.append(qpu.alloc_qubit())
        for q in v1s:
            q.debug_val = random.randrange(2 ** m)
        for q in v2s:
            q.debug_val = random.randrange(2 ** m)
        gs.append(qpu.alloc_qubit())
        gs[0].debug_val = random.randrange(2 ** m)
        expected_1s = qpu.debug_read_vals(v1s, m)
        expected_2s = qpu.debug_read_vals(v2s, m)
        expected_gs = [g ^ (a < b) for a, b, g in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals(v2s, m), qpu.debug_read_vals(gs, m))]
        qpu.flip_if_lt(lhs=v1s, rhs=v2s, out=gs[0])

        assert qpu.debug_read_vals(v1s, m) == expected_1s
        assert qpu.debug_read_vals(v2s, m) == expected_2s
        assert qpu.debug_read_vals(gs, m) == expected_gs


@pytest.mark.parametrize('n', [1, 2, 5, 10, 20])
def test_flip_if_lt_const(n: int):
    qpu = QPU()
    for _ in range(200):
        m = 20
        v1s = []
        const = random.randrange(2**n)
        gs = []
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
        gs.append(qpu.alloc_qubit())
        control = qpu.alloc_qubit()
        for q in v1s:
            q.debug_val = random.randrange(2 ** m)
        for q in gs:
            q.debug_val = random.randrange(2 ** m)
        control.debug_val = random.randrange(2 ** m)
        expected_1s = qpu.debug_read_vals(v1s, m)
        expected_cs = qpu.debug_read_vals([control], m)
        expected_gs = [g ^ (c and (a < const)) for a, g, c in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals(gs, m), qpu.debug_read_vals([control], m))]
        qpu.flip_if_lt_const(lhs=v1s, rhs_const=const, out=gs[0], control=control)
        assert qpu.debug_read_vals(v1s, m) == expected_1s
        assert qpu.debug_read_vals(gs, m) == expected_gs
        assert qpu.debug_read_vals([control], m) == expected_cs


@pytest.mark.parametrize('n', [1, 2, 5, 10, 20])
def test_flip_if_ge_const(n: int):
    qpu = QPU()
    for _ in range(200):
        m = 20
        v1s = []
        const = random.randrange(2**n)
        gs = []
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
        gs.append(qpu.alloc_qubit())
        for q in v1s:
            q.debug_val = random.randrange(2 ** m)
        for q in gs:
            q.debug_val = random.randrange(2 ** m)
        expected_1s = qpu.debug_read_vals(v1s, m)
        expected_gs = [g ^ (a >= const) for a, g in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals(gs, m))]
        qpu.flip_if_ge_const(lhs=v1s, rhs_const=const, out=gs[0])
        assert qpu.debug_read_vals(v1s, m) == expected_1s
        assert qpu.debug_read_vals(gs, m) == expected_gs


@pytest.mark.parametrize('n', [1, 2, 5, 10, 20])
def test_iadd_const_mod_if(n: int):
    qpu = QPU()
    for _ in range(20):
        m = 20
        v1s = []
        modulus = random.randrange(2**(n - 1), 2**n)
        control = qpu.alloc_qubit()
        const = random.randrange(modulus)
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
        control.debug_val = random.randrange(2 ** m)
        for k in range(m):
            qpu.debug_set_val(v1s, mi=k, val=random.randrange(modulus))
        expected_1s = [(v + const * c) % modulus for v, c in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals([control], m))]
        qpu.iadd_const_mod_if(target=v1s, control=control, offset_const=const, modulus=modulus)
        assert qpu.debug_read_vals(v1s, m) == expected_1s


@pytest.mark.parametrize('n', [3, 4, 5, 10, 20])
def test_imul_const_mod_if(n: int):
    for _ in range(20):
        qpu = QPU()
        m = 20
        v1s = []
        v2s = []
        cs = []
        while True:
            modulus = random.randrange(2 ** (n - 1), 2 ** n)
            const = random.randrange(2, modulus)
            if math.gcd(const, modulus) == 1:
                break
        for k in range(n):
            v1s.append(qpu.alloc_qubit())
            v2s.append(qpu.alloc_qubit())
        cs.append(qpu.alloc_qubit())
        cs[-1].debug_val = random.randrange(2 ** m)
        for k in range(m):
            qpu.debug_set_val(v1s, mi=k, val=random.randrange(1, modulus))
        expected_vs1 = [q * const % modulus if c else q for q, c in zip(qpu.debug_read_vals(v1s, m), qpu.debug_read_vals(cs, m))]
        expected_cs = [q for q in qpu.debug_read_vals(cs, m)]
        qpu.imul_const_mod_if(target=v1s, const_factor=const, modulus=modulus, control=cs[0])
        assert qpu.debug_read_vals(cs, m) == expected_cs
        assert qpu.debug_read_vals(v1s, m) == expected_vs1


def test_inv_qft_measure():
    qpu = QPU()
    qpu.recorded = []
    qs = [qpu.alloc_qubit() for _ in range(4)]
    qpu.inv_qft_measure(qs)
    assert qpu.recorded == [
        ('H', 3),
        ('M', 3),
        ('CPHASE', -0.5, 2, 3),
        ('H', 2),
        ('M', 2),
        ('CPHASE', -0.5, 1, 2),
        ('CPHASE', -0.25, 1, 3),
        ('H', 1),
        ('M', 1),
        ('CPHASE', -0.5, 0, 1),
        ('CPHASE', -0.25, 0, 2),
        ('CPHASE', -0.125, 0, 3),
        ('H', 0),
        ('M', 0),
    ]
    assert qpu.to_qasm() == """
OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
h q[3];
measure q[3] -> c[0];
cp(pi * -0.5) q[2], q[3];
h q[2];
measure q[2] -> c[1];
cp(pi * -0.5) q[1], q[2];
cp(pi * -0.25) q[1], q[3];
h q[1];
measure q[1] -> c[2];
cp(pi * -0.5) q[0], q[1];
cp(pi * -0.25) q[0], q[2];
cp(pi * -0.125) q[0], q[3];
h q[0];
measure q[0] -> c[3];
    """.strip()


def test_generate_qasm_for_shor_circuit():
    qasm = generate_qasm_for_shor_circuit(n=15, g=2)
    assert qasm.endswith("""
ccx q[11], q[7], q[12];
h q[7];
measure q[7] -> c[0];
cp(pi * -0.5) q[6], q[7];
h q[6];
measure q[6] -> c[1];
cp(pi * -0.5) q[5], q[6];
cp(pi * -0.25) q[5], q[7];
h q[5];
measure q[5] -> c[2];
cp(pi * -0.5) q[4], q[5];
cp(pi * -0.25) q[4], q[6];
cp(pi * -0.125) q[4], q[7];
h q[4];
measure q[4] -> c[3];
cp(pi * -0.5) q[3], q[4];
cp(pi * -0.25) q[3], q[5];
cp(pi * -0.125) q[3], q[6];
cp(pi * -0.0625) q[3], q[7];
h q[3];
measure q[3] -> c[4];
cp(pi * -0.5) q[2], q[3];
cp(pi * -0.25) q[2], q[4];
cp(pi * -0.125) q[2], q[5];
cp(pi * -0.0625) q[2], q[6];
cp(pi * -0.03125) q[2], q[7];
h q[2];
measure q[2] -> c[5];
cp(pi * -0.5) q[1], q[2];
cp(pi * -0.25) q[1], q[3];
cp(pi * -0.125) q[1], q[4];
cp(pi * -0.0625) q[1], q[5];
cp(pi * -0.03125) q[1], q[6];
cp(pi * -0.015625) q[1], q[7];
h q[1];
measure q[1] -> c[6];
cp(pi * -0.5) q[0], q[1];
cp(pi * -0.25) q[0], q[2];
cp(pi * -0.125) q[0], q[3];
cp(pi * -0.0625) q[0], q[4];
cp(pi * -0.03125) q[0], q[5];
cp(pi * -0.015625) q[0], q[6];
cp(pi * -0.0078125) q[0], q[7];
h q[0];
measure q[0] -> c[7];
    """.strip())

    assert qasm.startswith("""
OPENQASM 2.0;
include "qelib1.inc";
qreg q[24];
creg c[8];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
x q[8];
ccx q[0], q[8], q[16];
    """.strip())

    qasm = generate_qasm_for_shor_circuit(n=15, g=2, cphase_angle_cutoff=2**-3)

    assert qasm.endswith("""
ccx q[11], q[7], q[12];
h q[7];
measure q[7] -> c[0];
cp(pi * -0.5) q[6], q[7];
h q[6];
measure q[6] -> c[1];
cp(pi * -0.5) q[5], q[6];
cp(pi * -0.25) q[5], q[7];
h q[5];
measure q[5] -> c[2];
cp(pi * -0.5) q[4], q[5];
cp(pi * -0.25) q[4], q[6];
cp(pi * -0.125) q[4], q[7];
h q[4];
measure q[4] -> c[3];
cp(pi * -0.5) q[3], q[4];
cp(pi * -0.25) q[3], q[5];
cp(pi * -0.125) q[3], q[6];
h q[3];
measure q[3] -> c[4];
cp(pi * -0.5) q[2], q[3];
cp(pi * -0.25) q[2], q[4];
cp(pi * -0.125) q[2], q[5];
h q[2];
measure q[2] -> c[5];
cp(pi * -0.5) q[1], q[2];
cp(pi * -0.25) q[1], q[3];
cp(pi * -0.125) q[1], q[4];
h q[1];
measure q[1] -> c[6];
cp(pi * -0.5) q[0], q[1];
cp(pi * -0.25) q[0], q[2];
cp(pi * -0.125) q[0], q[3];
h q[0];
measure q[0] -> c[7];
    """.strip())


def test_cost():
    import qiskit
    import qiskit_ibm_runtime

    service = qiskit_ibm_runtime.QiskitRuntimeService(
        channel="local",
    )

    backend = service.backend("fake_sherbrooke")

    b = 4
    g = 2
    n = 2**b - 1
    qc = qiskit.QuantumCircuit.from_qasm_str(generate_qasm_for_shor_circuit(n=n, g=g, cphase_angle_cutoff=0))
    qc_transpiled, = qiskit.compiler.transpile([qc], backend=backend)
    counter = collections.Counter()
    for e in qc_transpiled:
        counter[e.name] += 1
    assert 43000 < counter['ecr'] < 45000
