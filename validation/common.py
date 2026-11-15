"""Shared helpers for the verification / validation scripts."""
from __future__ import annotations

import json
import os
import sys
import copy

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from aharrsim import case_from_dict, Simulator, DAY, YEAR, INCH  # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)


def save(name, obj):
    def conv(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        raise TypeError(type(o))
    with open(os.path.join(RESULTS, name + ".json"), "w") as fh:
        json.dump(obj, fh, indent=1, default=conv)


def load(name):
    with open(os.path.join(RESULTS, name + ".json")) as fh:
        return json.load(fh)


def run(d, progress=False):
    case = case_from_dict(copy.deepcopy(d))
    sim = Simulator(case)
    res = sim.run(progress=progress)
    return sim, res


def C(TK):
    return np.asarray(TK) - 273.15
