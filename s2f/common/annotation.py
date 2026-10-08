"""What counts as an uncharacterized protein, defined once for every module.

M1's `--seq-cap` exempts these proteins and M2's `annotation_gap` pays for them. If the two
used different tests, the cap could still remove proteins the score rewards (issue #85).
"""

from __future__ import annotations

import re

UNCHARACTERIZED = re.compile(
    r"hypothetical|uncharacteri[sz]ed|putative|unknown function|unnamed|"
    r"\bDUF\d+|protein of unknown|predicted protein",
    re.IGNORECASE,
)


def is_uncharacterized(product: str | None) -> bool:
    """True for a hypothetical or otherwise uncharacterized product.

    An empty product is not matched, as M2 has always scored it.
    """
    return bool(UNCHARACTERIZED.search(product or ""))
