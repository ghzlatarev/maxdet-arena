from __future__ import annotations

import itertools
import random
import unittest

from maxdet.exact import bareiss_determinant, modular_determinant


def determinant_by_permutations(matrix: list[list[int]]) -> int:
    total = 0
    for permutation in itertools.permutations(range(len(matrix))):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(len(matrix))
            for right in range(left + 1, len(matrix))
        )
        product = 1
        for row, column in enumerate(permutation):
            product *= matrix[row][column]
        total += -product if inversions % 2 else product
    return total


class ExactArithmeticTests(unittest.TestCase):
    def test_bareiss_matches_definition_on_seeded_small_matrices(self) -> None:
        randomizer = random.Random(2301)
        for order in range(1, 7):
            for _ in range(8):
                matrix = [
                    [randomizer.randint(-4, 4) for _ in range(order)]
                    for _ in range(order)
                ]
                expected = determinant_by_permutations(matrix)
                self.assertEqual(bareiss_determinant(matrix), expected)

    def test_modular_algorithm_matches_exact_residue(self) -> None:
        randomizer = random.Random(2302)
        for order in range(1, 8):
            for _ in range(8):
                matrix = [
                    [randomizer.randint(-20, 20) for _ in range(order)]
                    for _ in range(order)
                ]
                exact = bareiss_determinant(matrix)
                for prime in (101, 1009, 998244353):
                    self.assertEqual(
                        modular_determinant(matrix, prime),
                        exact % prime,
                    )

    def test_pivoting_and_singular_cases(self) -> None:
        self.assertEqual(bareiss_determinant([[0, 1], [1, 0]]), -1)
        self.assertEqual(bareiss_determinant([[1, 2], [2, 4]]), 0)


if __name__ == "__main__":
    unittest.main()
