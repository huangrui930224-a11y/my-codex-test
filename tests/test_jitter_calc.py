import math
import unittest

from src.utils.jitter_calc import calculate_rms_jitter


class TestCalculateRmsJitter(unittest.TestCase):
    def test_calculate_rms_jitter_with_mixed_values(self) -> None:
        samples = [1.0, -1.0, 1.0, -1.0]
        self.assertAlmostEqual(calculate_rms_jitter(samples), 1.0)

    def test_calculate_rms_jitter_with_decimal_values(self) -> None:
        samples = [0.1, 0.2, 0.3]
        expected = math.sqrt((0.01 + 0.04 + 0.09) / 3)
        self.assertAlmostEqual(calculate_rms_jitter(samples), expected)

    def test_calculate_rms_jitter_raises_on_empty_samples(self) -> None:
        with self.assertRaises(ValueError):
            calculate_rms_jitter([])


if __name__ == "__main__":
    unittest.main()
