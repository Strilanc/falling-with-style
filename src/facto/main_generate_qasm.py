import argparse

from facto.qpu import QPU


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', required=True, type=int, help="The number to factor.")
    parser.add_argument('--g', required=True, type=int, help="The randomly chosen g, coprime to n.")
    args = parser.parse_args()
    qpu = QPU()
    qpu.recorded = []
    qpu.perform_shor_quantum_part(
        g=args.g,
        modulus=args.n,
    )
    print(qpu.to_qasm(cphase_angle_cutoff=0))



if __name__ == '__main__':
    main()