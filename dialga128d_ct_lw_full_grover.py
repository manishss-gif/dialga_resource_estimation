"""
Grover key search for Dialga-128† (16 rounds) with the CT-LW circuit.

    python3 dialga128d_ct_lw_full_grover.py --resource     resource report (default)
    python3 dialga128d_ct_lw_full_grover.py --all-tests    functional tests, then the report
"""

import sys

from dialga_grover_common import GroverVariant, columns, grover_main
from qdialga128d_ct_lw import allocate_work, encrypt
from dialga128d_ct_lw_grover import dialga128d_encrypt_ct_low_width_inverse, dialga128d_grover_oracle_ct_lw


def encrypt_inverse(eng, state, tweaks, key, work):
    dialga128d_encrypt_ct_low_width_inverse(eng, state, *tweaks, key, *work)


def oracle(eng, keys, states, tweaks, works, cmp_ancillas, cmp_helpers, ciphertexts, test_flag=None):
    dialga128d_grover_oracle_ct_lw(
        eng, *keys, states, *columns(tweaks), *columns(works),
        cmp_ancillas, cmp_helpers, ciphertexts, test_flag,
    )


VARIANT = GroverVariant(
    name="Dialga-128†",
    profile="CT-LW",
    # (plaintext, tweak blocks, ciphertext) under the common key; tuple 1 is the official test vector
    tuples=[
        ("00112233445566778899aabbccddeeff",
         ("2233445566778899aabbccddee00ff11",),
         "838407143af9a876fbdc6be378e9045b"),
        ("0123456789abcdef0123456789abcdef",
         ("00001111222233334444555566668888",),
         "d260ca248048ff53640bc3032558ef53"),
        ("fedcba9876543210fedcba9876543210",
         ("ffffeeeeddddccccbbbbaaaa99998888",),
         "89ae99b018f90ed820303a8bc207963e"),
    ],
    allocate_work=allocate_work,
    encrypt=encrypt,
    encrypt_inverse=encrypt_inverse,
    oracle=oracle,
    diffusion_workspace=(255, 128),
    # (#CNOT, #1qClifford, #T, #Measure, T-depth, full depth, width)
    reference={
        "encryption": (47872, 28197, 12800, 3200, 3200, 38404, 514),
        "oracle": (292086, 172661, 78332, 19583, 6400, 76807, 2117),
        "diffusion": (2550, 3065, 1020, 255, 8, 63, 639),
        "iteration": (294636, 175726, 79352, 19838, 6408, 76870, 2117),
    },
)


if __name__ == "__main__":
    grover_main(VARIANT, sys.argv[1:])
