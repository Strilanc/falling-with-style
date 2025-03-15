import collections
import fractions
import math
import os
import random
from typing import Literal

import numpy as np
import sympy
from matplotlib import pyplot as plt

from facto.qpu import QPU


def brute_force_period(*, g: int, n: int):
    x = g
    p = 1
    while x != 1:
        x *= g
        x %= n
        p += 1
    return p


def simulate_sampling_frequency_of_periodic_signal(*, period: int, qubits: int):
    # Make the periodic signal.
    state_vector = np.zeros(1 << qubits, dtype=np.complex128)
    state_vector[::period] = 1
    state_vector /= np.sqrt(np.sum(state_vector * state_vector))
    # Fourier transform to get the frequencies.
    state_vector = np.fft.fft(state_vector, norm='ortho')
    # Sample with probability equal to squared magnitude.
    probs = np.real(state_vector * np.conj(state_vector))
    return int(np.random.choice(range(len(probs)), p=probs))


QPU_RECORDED = {}
QPU_INDEX = collections.Counter()

def sample_shor_quantum_circuit(
        n: int,
        g: int,
        *,
        strategy: Literal['recorded_from_qpu', 'rng', 'simulate_ideal'],
) -> int | None:
    if strategy == 'rng':
        return random.randrange(1 << (2 * n.bit_length()))
    elif strategy == 'simulate_ideal':
        period = brute_force_period(g=g, n=n)
        return simulate_sampling_frequency_of_periodic_signal(
            period=period,
            qubits=2 * n.bit_length(),
        )
    elif strategy == 'recorded_from_qpu':
        key = (n, g)
        recorded = QPU_RECORDED.get(key, [])
        k = QPU_INDEX[key]
        if k < len(recorded):
            return recorded[k]
        return None
    else:
        raise NotImplementedError(f'{strategy=}')


def attempt_quantum_factoring(
        n: int,
        *,
        strategy: Literal['recorded_from_qpu', 'rng', 'simulate_ideal'],
) -> tuple[Literal['is_prime', 'factored_by_even', 'factored_by_gcd_hit', 'factored_by_period', 'collect_qpu_sample_for_g', 'retry'], int]:
    if sympy.isprime(n):
        return 'is_prime', 0

    # Shor's algorithm requires n to be odd, so even integers are solved by explicit checking.
    if n % 2 == 0:
        return 'factored_by_even', 2

    # Pick a random generator. If it happens to share factors with n, exit early.
    g = random.randrange(2, n - 2)
    gcd = math.gcd(g, n)
    if gcd != 1:
        return 'factored_by_gcd_hit', gcd

    # Perform quantum period finding.
    sample = sample_shor_quantum_circuit(n, g, strategy=strategy)
    if sample is None:
        return 'collect_qpu_sample_for_g', g  # Failed to get a quantum sample.
    fraction = fractions.Fraction(sample, 1 << (2 * n.bit_length()))
    potential_period = fraction.limit_denominator(n - 1).denominator

    # Try to extract a quadratic residue by looking halfway through the cycle
    h = pow(g, potential_period // 2, n)
    gcd = math.gcd(h + 1, n)
    if 1 < gcd < n:
        return 'factored_by_period', gcd

    # This attempt failed. The caller should try again.
    return 'retry', 0


def generate_qasm_for_shor_circuit(*, n: int, g: int, cphase_angle_cutoff: float = 0) -> str:
    qpu = QPU()
    qpu.recorded = []
    qpu.perform_shor_quantum_part(
        g=g,
        modulus=n,
    )
    return qpu.to_qasm(cphase_angle_cutoff=cphase_angle_cutoff)


def run_problems():
    max_value = 256
    # primes = list(sympy.primerange(3, max_value))
    # problems = sorted({
    #     p1 * p2
    #     for p1, p2 in itertools.combinations(primes, 2)
    #     if p1 * p2 < max_value
    # })
    problems = range(15, max_value)
    collects = []
    for strategy in ['recorded_from_qpu']:  #, 'rng', 'simulate_ideal']:
        successes = 0
        attempts = 0
        for n in problems:
            # Force the same random seed for every problem, to ensure recorded shots are relevant.
            random.seed(2025_04_01)
            np.random.seed(2025_04_01)

            while True:
                attempts += 1
                kind, v = attempt_quantum_factoring(n, strategy=strategy)
                if kind.startswith('factored'):
                    assert 1 < v < n
                    assert v * (n // v) == n
                    successes += 1
                    break
                elif kind == 'retry':
                    pass
                elif kind == 'is_prime':
                    break
                elif kind == 'collect_qpu_sample_for_g':
                    collects.append((n, v))
                    break
                else:
                    raise NotImplementedError(f'{kind=}')
        # print(strategy, successes / attempts)

    return collects


def plot_shots():
    fig, ax = plt.subplots()
    fig: plt.Figure
    ax: plt.Axes
    for strategy in ['recorded_from_qpu', 'rng', 'simulate_ideal']:
        xs = []
        ys = []
        samples = 0
        problems = range(15, 256)
        for n in problems:
            # Force the same random seed for every problem, to ensure recorded shots are relevant.
            random.seed(2025_04_01)
            np.random.seed(2025_04_01)

            while True:
                kind, v = attempt_quantum_factoring(n, strategy=strategy)
                if kind.startswith('factored'):
                    assert 1 < v < n
                    assert v * (n // v) == n
                    if kind == 'factored_by_period':
                        samples += 1
                    break
                elif kind == 'retry':
                    samples += 1
                elif kind == 'is_prime':
                    break
                else:
                    raise NotImplementedError(f'{kind=}')
            xs.append(n - 1)
            ys.append(samples)
            xs.append(n)
            ys.append(samples)
        print(strategy, ys[-1])
        ax.plot(xs, ys, linewidth=3, label={
            'rng': 'Random Number Generator',
            'recorded_from_qpu': 'Real Quantum Computer',
            'simulate_ideal': 'Simulated Noiseless Quantum Computer',
        }.get(strategy, strategy))
    ax.set_xlabel("Factored All Numbers Up To", fontsize=24)
    ax.set_ylabel("Number of Samples Collected", fontsize=24)
    ax.grid()
    ax.grid(which='minor')
    ax.set_xlim(15, 256)
    ax.set_ylim(0, 130)
    ax.legend(fontsize=24)
    ax.set_xticks(range(0, 256, 50))
    ax.set_xticklabels([str(e) for e in range(0, 256, 50)], fontsize=18)
    ax.set_yticks(range(0, 130, 10))
    ax.set_yticklabels([str(e) for e in range(0, 130, 10)], fontsize=18)
    fig.set_dpi(200)
    fig.set_size_inches(12, 8)
    fig.tight_layout()
    fig.savefig('plot.png')
    plt.show()


def collect_job_results():
    from qiskit_ibm_runtime import QiskitRuntimeService

    token = os.environ['QISKIT_TOKEN']
    assert token is not None
    service = QiskitRuntimeService(
        # channel="local",
        channel="ibm_quantum",
        token=token,
    )

    for job_id, n, g, v, _ in ACTIVE_JOBS:
        job = service.job(job_id)
        if str(job.status()) == 'DONE':
            result = job.result()
            for sub_result in result:
                # print(sub_result.metadata['circuit_metadata'])
                # print(sub_result.data['c'])
                for shot in sub_result.data['c'].array:
                    bits = []
                    for k, b in enumerate(shot[::-1]):
                        for k2 in range(8):
                            bits.append((b >> k2) & 1)
                    bits = bits[:n.bit_length() * 2]
                    bit_str = '0b' + ''.join(str(e) for e in bits)
                    print(f"""    ({job_id!r}, {n}, {g}, {v}, {bit_str}),""")
        else:
            print(job_id, "not done", job.status())


def send_jobs_to_ibm(
        to_collect: list[tuple[int, int]],
        *,
        cphase_angle_cutoff: float,
        circuit_version: int,
):
    import qiskit
    import qiskit_ibm_runtime

    token = os.environ['QISKIT_TOKEN']
    assert token is not None
    service = qiskit_ibm_runtime.QiskitRuntimeService(
        # channel="local",
        channel="ibm_quantum",
        token=token,
    )

    # backend = service.backend("fake_sherbrooke")
    backend = service.backend("ibm_sherbrooke")

    sampler = qiskit_ibm_runtime.SamplerV2(mode=backend)
    for n, g in to_collect:
        v = circuit_version
        qc = qiskit.QuantumCircuit.from_qasm_str(generate_qasm_for_shor_circuit(n=n, g=g, cphase_angle_cutoff=cphase_angle_cutoff))
        qc_transpiled, = qiskit.compiler.transpile([qc], backend=backend)
        qc_transpiled.metadata['n'] = n
        qc_transpiled.metadata['g'] = g
        qc_transpiled.metadata['v'] = v
        job = sampler.run([qc_transpiled], shots=1)
        with open('tmp.txt', 'a') as f:
            print(f"""        ({job.job_id()!r}, {n}, {g}, {v}, None)""", file=f)
        print(f"""        ({job.job_id()!r}, {n}, {g}, {v}, None),""")


def store_collected_shots_from_jobs():
    for job_id, n, g, v, result in [
        ('cza169stp60g008h13r0', 15, 7, 1, 0b01000111),
        ('cza16m2b7tt0008g48ag', 35, 27, 1, 0b111010011111),
        ('cza16yvhfwp00087wzz0', 39, 23, 1, 0b111110101010),
        ('cza1ekjhfwp00087x0cg', 35, 24, 1, 0b000100100000),
        ('cza1ey3hfwp00087x0dg', 39, 32, 1, 0b001011111100),
        ('cza1f9d4spc00086tpwg', 51, 47, 1, 0b010111000100),
        ('cza1jc94spc00086tq5g', 35, 12, 1, 0b101011110110),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b101000001010),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b001011010101),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b001011101110),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b111100101110),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b110011111010),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b101001100001),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b010111101011),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b010110001111),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b100000000110),
        ('cza1kcdb7tt0008g49c0', 51, 23, 1, 0b001111001111),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b000000001100),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b110010100010),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b011000000010),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b000100100100),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b101001011011),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b000101011011),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b110111011010),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b101111100111),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b001001001110),
        ('cza1kqetp60g008h14s0', 55, 52, 1, 0b000011011001),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b010111110011),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b100011110010),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b100100111111),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b010010001100),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b111110011000),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b101010001100),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b001000101000),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b111011110110),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b000000101011),
        ('cza1m2g4spc00086tqc0', 57, 52, 1, 0b110001010110),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b10100000010001),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b01010110001000),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b11000111010000),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b00110000100000),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b01110001000101),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b11110011000100),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b00001011010101),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b00110110000100),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b10111000011110),
        ('cza1mmjb7tt0008g49h0', 77, 45, 1, 0b00100001011111),
        ('cza1rfh4spc00086tqzg', 35, 29, 2, 0b011000110111),
        ('cza1rtkhfwp00087x1eg', 51, 32, 2, 0b110000010110),
        ('cza1s5m4spc00086tr1g', 55, 47, 2, 0b101000110110),
        ('cza1sgpqadq0008c7830', 57, 47, 2, 0b000110010001),
        ('cza1t30qadq0008c7840', 91, 45, 2, 0b00100000010001),
        ('cza1yd94spc00086trag', 35, 17, 3, 0b111100100100),
        ('cza1yrkb7tt0008g4af0', 55, 23, 3, 0b110000100110),
        ('cza1z3wb7tt0008g4agg', 57, 23, 3, 0b111111000110),
        ('cza1zpp4spc00086trf0', 95, 92, 3, 0b11100011010110),

        ('cza2a7rb7tt0008g4bz0', 55, 32, 3, 0b101100111000),
        ('cza2arv4spc00086tsmg', 111, 103, 3, 0b10000110111110),
        ('cza2bhe4spc00086tsp0', 115, 103, 3, 0b10010011110111),
        ('cza2c1rqadq0008c79z0', 119, 103, 3, 0b11001000111100),
        ('cza2ckatp60g008h17r0', 123, 103, 3, 0b00001110111011),
        ('cza2dbnr3jrg008nwbw0', 129, 103, 3, 0b1010101101100011),
        ('cza2e50r3jrg008nwby0', 133, 89, 3, 0b1110011110110100),
        ('cza2eyk4spc00086tsvg', 141, 89, 3, 0b0001110111110010),
        ('cza2frf4spc00086tsxg', 143, 89, 3, 0b0010001101101100),
        ('cza2gj2r3jrg008nwc60', 145, 89, 3, 0b0010000001001101),
        ('cza2hbxqadq0008c7akg', 155, 89, 3, 0b1000110000000010),
        ('cza2j50qadq0008c7ap0', 159, 89, 3, 0b0111101100111000),
        ('cza2jz3hfwp00087x4ag', 161, 89, 3, 0b1001000001110011),
        ('cza2krf4spc00086tt8g', 177, 89, 3, 0b0010111000001000),
        ('cza2mjthfwp00087x4f0', 183, 89, 3, 0b0000001110100111),
        ('cza2ncnr3jrg008nwch0', 185, 89, 3, 0b0101110111110011),
        ('cza2p6r4spc00086ttk0', 187, 183, 3, 0b0100011111001001),
        ('cza2q0m4spc00086ttmg', 203, 183, 3, 0b0101000010001001),
        ('cza2qv7b7tt0008g4ddg', 205, 183, 3, 0b0110111101001001),
        ('cza2rn2tp60g008h18zg', 209, 205, 3, 0b1010001100001101),
        ('cza2sf54spc00086ttw0', 213, 205, 3, 0b0000011000101101),
        ('cza2t914spc00086tv2g', 217, 205, 3, 0b0110000110101000),
        ('cza2v44tp60g008h19cg', 219, 205, 3, 0b1101101000000110),
        ('cza2vy74spc00086tvf0', 221, 205, 3, 0b0100000010001110),
        ('cza2wr3hfwp00087x5hg', 237, 205, 3, 0b0101110001001010),
        ('cza2xhptp60g008h19x0', 247, 205, 3, 0b0101011110101001),
        ('cza2yb14spc00086tvvg', 249, 205, 3, 0b0010100111011001),
        ('cza2z544spc00086tvy0', 253, 205, 3, 0b0111100011100111),
        ('czaagsktp60g008h2wtg', 25, 12, 4, 0b1010011000),
        ('czaah4wb7tt0008g6170', 45, 23, 4, 0b001010101100),
        ('czaahg64spc00086wedg', 49, 23, 4, 0b111001110111),
        ('czaahvq4spc00086weeg', 63, 52, 4, 0b001101010111),
        ('czaajdhqadq0008c8zg0', 99, 92, 4, 0b11010110100010),
        ('czaak04hfwp00087ysy0', 105, 92, 4, 0b01000001110111),
        ('czaakjyqadq0008c8zp0', 111, 92, 4, 0b00110111000011),
        ('czaam5rhfwp00087yt70', 117, 103, 4, 0b01100000100000),
        ('czaamrb4spc00086wf2g', 119, 92, 4, 0b01101111001111),
        ('czaan9xb7tt0008g623g', 121, 103, 4, 0b10010011100011),
        ('czaanx7hfwp00087ytt0', 123, 92, 4, 0b01110101010011),
        ('czaapga4spc00086wfjg', 125, 103, 4, 0b10101010101110),
        ('czaaqbdtp60g008h2ybg', 129, 92, 4, 0b1110110001010110),
        ('czaar6r4spc00086wg60', 133, 124, 4, 0b0100111100010000),
        ('czaas24qadq0008c90zg', 135, 89, 4, 0b0001101100011000),
        ('czaasxzqadq0008c916g', 143, 124, 4, 0b0010000001110011),
        ('czaatsbr3jrg008ny4a0', 147, 89, 4, 0b1111010011110110),
        ('czaavn6b7tt0008g6400', 153, 89, 4, 0b1001010110110010),
        ('czaawj2r3jrg008ny4z0', 159, 124, 4, 0b1111000101110100),
        ('czaaxcxtp60g008h30bg', 161, 124, 4, 0b0100110011100110),
        ('czaay81tp60g008h30r0', 165, 89, 4, 0b1100001111111001),
        ('czaaz2w4spc00086wjkg', 169, 89, 4, 0b1101110110010101),
        ('czaazxqqadq0008c93f0', 171, 89, 4, 0b1100110111111001),
        ('czab0rvr3jrg008ny6ng', 175, 89, 4, 0b1010011101011111),
        ('czab1k6qadq0008c94a0', 177, 124, 4, 0b0111100111110011),
        ('czab2e1r3jrg008ny79g', 187, 89, 4, 0b1100101011001010),
        ('czab38n4spc00086wm80', 203, 89, 4, 0b0110110000111110),
        ('czab4384spc00086wmeg', 205, 89, 4, 0b0110100100110100),
        ('czab4y34spc00086wmkg', 209, 183, 4, 0b0000000010110001),
        ('czab5sf4spc00086wmt0', 217, 183, 4, 0b1100000001011010),
        ('czab6matp60g008h33f0', 231, 205, 4, 0b0111010110010100),
        ('czab7fxr3jrg008ny8cg', 243, 205, 4, 0b0101001000100100),
        ('czab8ahr3jrg008ny8kg', 247, 183, 4, 0b1100001110010100),
        ('czab964r3jrg008ny8r0', 253, 183, 4, 0b0111011011101100),

        ('czabfyztp60g008h3580', 25, 17, 4, 0b1011101111),
        ('czabga1qadq0008c97rg', 49, 32, 4, 0b111001101010),
        ('czabgn2tp60g008h35ag', 63, 47, 4, 0b001000110001),
        ('czabh7cqadq0008c97y0', 117, 92, 4, 0b11110000011100),
        ('czabhszqadq0008c982g', 119, 45, 4, 0b11000100001000),
        ('czabjc1hfwp00087z1fg', 121, 92, 4, 0b10110110110110),
        ('czabjzkb7tt0008g6a00', 125, 92, 4, 0b01100010010010),
        ('czabkt74spc00086wqw0', 133, 18, 4, 0b0100010101011101),
        ('czabmnthfwp00087z1t0', 135, 124, 4, 0b0001011100101111),
        ('czabnhphfwp00087z200', 143, 18, 4, 0b1000000110010111),
        ('czabpc14spc00086wreg', 147, 124, 4, 0b0111110000111011),
        ('czabq6wb7tt0008g6ang', 161, 18, 4, 0b0011110101100111),
        ('czabr18r3jrg008nyba0', 165, 124, 4, 0b1010111010010001),
        ('czabrtvb7tt0008g6b2g', 169, 124, 4, 0b1011101001100101),
        ('czabsp64spc00086wrzg', 187, 124, 4, 0b1010100111101101),
        ('czabthttp60g008h375g', 203, 124, 4, 0b0100100010000110),
        ('czabvd5r3jrg008nyc50', 205, 124, 4, 0b0001110101010100),
        ('czabw9htp60g008h37e0', 209, 89, 4, 0b0101011111000010),
        ('czabx5wqadq0008c99r0', 217, 89, 4, 0b0110001101010110),
        ('czaby0rr3jrg008nycc0', 247, 89, 4, 0b1010011011110001),

        ('czac2xbr3jrg008nyd00', 25, 4, 4, 0b0111001010),
        ('czac37wr3jrg008nyd3g', 49, 6, 4, 0b110111011011),
        ('czac3s7b7tt0008g6d00', 121, 45, 4, 0b10101001001110),
        ('czac4kjb7tt0008g6db0', 143, 28, 4, 0b0101111010101111),
        ('czac5edr3jrg008nydj0', 169, 18, 4, 0b0000001101110000),
        ('czac69sqadq0008c9bg0', 187, 18, 4, 0b0001000011110101),
        ('czac75cb7tt0008g6dw0', 203, 18, 4, 0b0001011011000010),
        ('czac810b7tt0008g6e4g', 205, 18, 4, 0b0111011011011110),
        ('czac8wvr3jrg008nyef0', 209, 124, 4, 0b1101101110110000),
        ('czac9rq4spc00086ww90', 247, 223, 4, 0b1001011010001001),

        ('czaced1qadq0008c9ec0', 49, 8, 4, 0b011111010111),
        ('czacezvb7tt0008g6gag', 121, 112, 4, 0b01011101011001),
        ('czacfxzr3jrg008nygk0', 143, 129, 4, 0b0010000000101000),
        ('czacgsk4spc00086wy60', 209, 18, 4, 0b1000100110000111),
        ('czachn6hfwp00087z870', 247, 124, 4, 0b1101100110000101),

        ('czacn2wr3jrg008nyhsg', 49, 33, 4, 0b110011011011),
        ('czacnm6hfwp00087z910', 121, 63, 4, 0b10110000011110),
        ('czacpehtp60g008h3e1g', 209, 28, 4, 0b0000010001011011),
        ('czacqfnb7tt0008g6j1g', 247, 18, 4, 0b0001100100111011),

        ('czacwwvr3jrg008nykv0', 49, 34, 4, 0b011010011111),
        ('czacxdxtp60g008h3g0g', 121, 10, 4, 0b00101100000011),
        ('czacy7gqadq0008c9jh0', 209, 129, 4, 0b0100011100000100),

        ('czad060hfwp00087zc20', 121, 15, 4, 0b01000001011010),
        ('czad114hfwp00087zc90', 209, 131, 4, 0b1110000111111000),

        ('czad2qatp60g008h3hfg', 121, 65, 4, 0b00010101001101),
    ]:
        QPU_RECORDED.setdefault((n, g), []).append(result)


ACTIVE_JOBS: list[tuple[str, int, int, int, None]] = [

]


def main():
    # Plot collected data.
    store_collected_shots_from_jobs()
    plot_shots()

    # # Collect results from job ids listed in ACTIVE_JOBS and print them.
    # # They should be manually put into the array in `store_collected_shots_from_jobs`.
    # collect_job_results()

    # # Determine jobs to submit, submit them, and print them.
    # # They should be manually put into the ACTIVE_JOBS array.
    # store_collected_shots_from_jobs()
    # to_collect = run_problems()
    # print("to collect:", len(to_collect), to_collect)
    # send_jobs_to_ibm(to_collect, cphase_angle_cutoff=0, circuit_version=4)


if __name__ == '__main__':
    main()
