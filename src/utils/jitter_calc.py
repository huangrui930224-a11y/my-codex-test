"""Jitter calculation utilities."""

from __future__ import annotations

import math
from typing import Iterable


def calculate_rms_jitter(jitter_samples: Iterable[float]) -> float:
    """Calculate RMS (root mean square) jitter.

    The RMS jitter is defined as:

        sqrt((x1^2 + x2^2 + ... + xn^2) / n)

    where each ``x`` is a jitter sample (typically a time deviation from an
    ideal edge position).

    Args:
        jitter_samples: Iterable of jitter values.

    Returns:
        The RMS jitter as a float.

    Raises:
        ValueError: If ``jitter_samples`` is empty.
    """

    samples = [float(sample) for sample in jitter_samples]
    if not samples:
        raise ValueError("jitter_samples must not be empty")

    mean_square = sum(sample * sample for sample in samples) / len(samples)
    return math.sqrt(mean_square)
