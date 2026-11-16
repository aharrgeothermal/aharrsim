"""aharrsim - transient wellbore-reservoir-plant simulator for geothermal wells.

The core is a coaxial closed loop: cold water is pumped down the annulus and returns
hot up an (insulated) centre tubing; the rock supplies heat by conduction only.  The
same solver runs production, injection and U-tube wells (aharrsim.simulator), chains of
wells and reservoir elements (aharrsim.systems, aharrsim.reservoir), flow through a porous
or fractured block between two screens on one well (aharrsim.porous) and surface plants
(aharrsim.plant).  Water/steam properties are IAPWS-95 (CoolProp).  See README.md for the
model description, the verification and the benchmark results.
"""
from .case import Case, case_from_dict, DAY, YEAR, INCH
from .simulator import Simulator, Results
from .props import WaterIAPWS, ConstantFluid

__version__ = "2.4.1"


def load_case(path):
    """Read a YAML or JSON case file."""
    import json
    with open(path) as fh:
        text = fh.read()
    if str(path).endswith((".yaml", ".yml")):
        import yaml
        d = yaml.safe_load(text)
    else:
        d = json.loads(text)
    return case_from_dict(d)


def run_case(case_or_dict, progress=False):
    """Build and run a case; returns (Simulator, Results)."""
    case = case_from_dict(case_or_dict) if isinstance(case_or_dict, dict) else case_or_dict
    sim = Simulator(case)
    res = sim.run(progress=progress)
    return sim, res
