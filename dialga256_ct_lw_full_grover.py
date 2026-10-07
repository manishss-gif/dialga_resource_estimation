"""
Grover key search for Dialga-256 (20 rounds) with the CT-LW circuit.

    python3 dialga256_ct_lw_full_grover.py --resource     resource report (default)
    python3 dialga256_ct_lw_full_grover.py --all-tests    functional tests, then the report
"""

import sys

from dialga_grover_common import GroverVariant, columns, grover_main
from qdialga256_ct_lw import allocate_work, encrypt
from dialga256_ct_lw_grover import dialga256_encrypt_ct_low_width_inverse, dialga256_grover_oracle_ct_lw


def encrypt_inverse(eng, state, tweaks, key, work):
    dialga256_encrypt_ct_low_width_inverse(eng, state, *tweaks, key, *work)


def oracle(eng, keys, states, tweaks, works, cmp_ancillas, cmp_helpers, ciphertexts, test_flag=None):
    dialga256_grover_oracle_ct_lw(
        eng, *keys, states, *columns(tweaks), *columns(works),
        cmp_ancillas, cmp_helpers, ciphertexts, test_flag,
    )


VARIANT = GroverVariant(
    name="Dialga-256",
    profile="CT-LW",
    # (plaintext, tweak blocks, ciphertext) under the common key; tuple 1 is the official test vector
    tuples=[
        ("00112233445566778899aabbccddeeff",
         ("2233445566778899aabbccddeeff1122", "33445566778899aabbccddeeff001100"),
         "e129ff0920371753db65532540d06881"),
        ("0123456789abcdef0123456789abcdef",
         ("00001111222233334444555566668888", "9999aaaabbbbccccddddeeeeffff0000"),
         "c1bc3fa252cc52c8561fa825efd93a06"),
        ("fedcba9876543210fedcba9876543210",
         ("ffffeeeeddddccccbbbbaaaa99998888", "77776666555544443333222211110000"),
         "ad077bdcc5ae51550b3a91df7ef491ef"),
    ],
    allocate_work=allocate_work,
    encrypt=encrypt,
    encrypt_inverse=encrypt_inverse,
    oracle=oracle,
    diffusion_workspace=(255, 128),
    # (#CNOT, #1qClifford, #T, #Measure, T-depth, full depth, width)
    reference={
        "encryption": (75264, 43683, 19968, 4992, 4992, 59909, 642),
        "oracle": (456438, 265545, 121340, 30335, 9984, 119815, 2501),
        "diffusion": (2550, 3065, 1020, 255, 8, 63, 639),
        "iteration": (458988, 268610, 122360, 30590, 9992, 119878, 2501),
    },
)


if __name__ == "__main__":
    grover_main(VARIANT, sys.argv[1:])
