"""
Grover key search for Dialga-128 (20 rounds) with the CT-LD circuit.

    python3 dialga128_ct_ld_full_grover.py --resource     resource report (default)
    python3 dialga128_ct_ld_full_grover.py --all-tests    functional tests, then the report
"""

import sys

from dialga_grover_common import GroverVariant, columns, grover_main
from qdialga128_ct_ld import allocate_work, encrypt
from dialga128_ct_ld_grover import dialga128_encrypt_ct_low_depth_inverse, dialga128_grover_oracle_ct_ld


def encrypt_inverse(eng, state, tweaks, key, work):
    dialga128_encrypt_ct_low_depth_inverse(eng, state, *tweaks, key, *work)


def oracle(eng, keys, states, tweaks, works, cmp_ancillas, cmp_helpers, ciphertexts, test_flag=None):
    dialga128_grover_oracle_ct_ld(
        eng, *keys, states, *columns(tweaks), *columns(works),
        cmp_ancillas, cmp_helpers, ciphertexts, test_flag,
    )


VARIANT = GroverVariant(
    name="Dialga-128",
    profile="CT-LD",
    # (plaintext, tweak blocks, ciphertext) under the common key; tuple 1 is the official test vector
    tuples=[
        ("00112233445566778899aabbccddeeff",
         ("2233445566778899aabbccddee00ff11",),
         "a4a1ea948919d8996e13b1b365bb0ce6"),
        ("0123456789abcdef0123456789abcdef",
         ("00001111222233334444555566668888",),
         "f939d243557ddf38b1caef37ba506d90"),
        ("fedcba9876543210fedcba9876543210",
         ("ffffeeeeddddccccbbbbaaaa99998888",),
         "e3245f25f0d3e155641bff798a38127a"),
    ],
    allocate_work=allocate_work,
    encrypt=encrypt,
    encrypt_inverse=encrypt_inverse,
    oracle=oracle,
    diffusion_workspace=(255, 128),
    # (#CNOT, #1qClifford, #T, #Measure, T-depth, full depth, width)
    reference={
        "encryption": (59392, 34979, 15872, 3968, 88, 1212, 640),
        "oracle": (361206, 213309, 96764, 24191, 177, 2434, 2495),
        "diffusion": (2550, 3065, 1020, 255, 8, 63, 639),
        "iteration": (363756, 216374, 97784, 24446, 185, 2497, 2495),
    },
)


if __name__ == "__main__":
    grover_main(VARIANT, sys.argv[1:])
