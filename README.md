Code for 'Falling with Style: Factoring up to 255 "with" a Quantum Computer'
============================================================================

This repository contains the code used to collect the data shown in the
april fools paper 'Falling with Style: Factoring up to 255 "with" a
Quantum Computer'.

Usage instructions
==================

Install python dependencies:

```bash
pip install -r requirements.txt
```

Run the tests:

```bash
PYTHONPATH=src pytest src
```

Generate a factoring circuit:

```bash
PYTHONPATH=src python src/facto/main_generate_qasm.py --n 15 --g 2
```

Make the plot shown in the paper:

```bash
PYTHONPATH=src python src/facto/main_make_plot.py
```

Recollect the data from the quantum computer (very manual):

1. Consider changing the hardcoded random seed (currently `2025_04_01`) used in `src/facto/main_make_plot.py`.
2. Set the environment variable `QISKIT_TOKEN` to your ibm quantum token.
3. Delete the contents of the array in `store_collected_shots_from_jobs` in `src/facto/main_make_plot.py`
4. Modify the main method of `src/facto/main_make_plot.py`; comment the plot call and uncomment the call to `run_problems` and `send_jobs_to_ibm`.
5. Once the jobs finish submitting a few minutes later, add the printed contents to the `ACTIVE_JOBS` array, then remodify main to only have the call to `collect_job_results` uncommented.
6. Add the printed job results to the array in `store_collected_shots_from_job`, and clear `ACTIVE_JOBS`.
7. Goto step 4, unless no jobs were submitted (indicating no more samples requested for any number)
8. Remodify main to only call the plot command, run it, and look at your data
