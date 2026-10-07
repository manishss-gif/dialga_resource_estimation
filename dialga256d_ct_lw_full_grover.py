"""
Grover key search for Dialga-256† (16 rounds) with the CT-LW circuit.

    python3 dialga256d_ct_lw_full_grover.py --resource     resource report (default)
    python3 dialga256d_ct_lw_full_grover.py --all-tests    functional tests, then the report
"""

import sys

from dialga_grover_common import GroverVariant, columns, grover_main
from qdialga256d_ct_lw import allocate_work, encrypt
from dialga256d_ct_lw_grover import dialga256d_encrypt_ct_low_width_inverse, dialga256d_grover_oracle_ct_lw


def encrypt_inverse(eng, state, tweaks, key, work):
    dialga256d_encrypt_ct_low_width_inverse(eng, state, *tweaks, key, *work)


def oracle(eng, keys, states, tweaks, works, cmp_ancillas, cmp_helpers, ciphertexts, test_flag=None):
    dialga256d_grover_oracle_ct_lw(
        eng, *keys, states, *columns(tweaks), *columns(works),
        cmp_ancillas, cmp_helpers, ciphertexts, test_flag,
    )


VARIANT = GroverVariant(
    name="Dialga-256†",
    profile="CT-LW",
    # (plaintext, tweak blocks, ciphertext) under the common key; tuple 1 is the official test vector
    tuples=[
        ("00112233445566778899aabbccddeeff",
         ("2233445566778899aabbccddeeff1122", "33445566778899aabbccddeeff001100"),
         "833367ce5fd49a7b1f0dd9ec62720018"),
        ("0123456789abcdef0123456789abcdef",
         ("00001111222233334444555566668888", "9999aaaabbbbccccddddeeeeffff0000"),
         "346db22e68eda5a1a02edded78cb2904"),
        ("fedcba9876543210fedcba9876543210",
         ("ffffeeeeddddccccbbbbaaaa99998888", "77776666555544443333222211110000"),
         "402546b6ceb2545408b815501250360c"),
    ],
    allocate_work=allocate_work,
    encrypt=encrypt,
    encrypt_inverse=encrypt_inverse,
    oracle=oracle,
    diffusion_workspace=(255, 128),
    # (#CNOT, #1qClifford, #T, #Measure, T-depth, full depth, width)
    reference={
        "encryption": (59904, 34725, 15872, 3968, 3968, 47621, 642),
        "oracle": (364278, 211835, 96764, 24191, 7936, 95239, 2501),
        "diffusion": (2550, 3065, 1020, 255, 8, 63, 639),
        "iteration": (366828, 214900, 97784, 24446, 7944, 95302, 2501),
    },
)


if __name__ == "__main__":
    grover_main(VARIANT, sys.argv[1:])
