"""Development-only oracle for fs-gas chemistry/tests/hydrogen.rs.

Run with Python and Cantera 3.2.0 in an isolated development environment:
    python3 data/reference/thermochem/hydrogen_cantera32.py [h2o2.yaml]
The optional input is the exact pinned Cantera source; otherwise use the
installed release's data file. Writes JSON to stdout, never changes Rust data.
No Cantera/Python dependency enters the FrankenSim production runtime.
"""
import hashlib
import json
import math
from pathlib import Path
import sys

import cantera as ct

assert ct.__version__ == "3.2.0"
SOURCE = "4a8358eb80cfeb50474386b5f9ec0b3a83519889"
BLOB = "65e6321a60bafb1feb4f3b1be5347e6d1e89c937"
NAMES = ["AR", "H", "H2", "H2O", "H2O2", "HO2", "N2", "O", "O2", "OH"]
path = Path(sys.argv[1]) if len(sys.argv) > 1 else next(
    Path(d) / "h2o2.yaml" for d in ct.get_data_directories()
    if (Path(d) / "h2o2.yaml").is_file()
)
raw = path.read_bytes()
assert hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest() == BLOB


def mechanism():
    gas = ct.Solution(str(path), "ohmech")
    # Explicitly reproduce ContinuousAt300K: preserve low-region coefficients,
    # Cp and all kinetic laws; reconcile only upper h/s integration constants.
    for i, species in enumerate(gas.species()):
        a = species.thermo.coeffs.copy()
        t, high, low = a[0], a[1:8], a[8:15]
        def h_r(c):
            return sum(c[j] * t ** (j + 1) / (j + 1) for j in range(5)) + c[5]
        def s_r(c):
            return c[0] * math.log(t) + sum(c[j] * t ** j / j for j in range(1, 5)) + c[6]
        a[6] += h_r(low) - h_r(high)
        a[7] += s_r(low) - s_r(high)
        species.thermo = ct.NasaPoly2(species.thermo.min_temp,
            species.thermo.max_temp, species.thermo.reference_pressure, a)
        gas.modify_species(i, species)
    return gas


cases = []
for temperature, duration in [(1000.0, 1e-3), (1200.0, 1e-4)]:
    gas = mechanism()
    gas.TPX = temperature, 101325.0, {"H2": 1.0, "O2": 1.0, "N2": 3.76}
    reactor = ct.IdealGasReactor(gas, energy="on", volume=0.001, clone=True)
    target = reactor.mass * reactor.phase.int_energy_mass
    network = ct.ReactorNet([reactor])
    network.rtol, network.atol = 1e-11, 1e-24
    network.max_steps, network.max_time_step = 200000, duration / 50.0
    network.advance(duration)
    phase = reactor.phase
    cases.append({
        "initial_temperature_k": temperature, "duration_s": duration,
        "temperature_k": phase.T, "pressure_pa": phase.P,
        "target_energy_j": target,
        "energy_residual_j": reactor.mass * phase.int_energy_mass - target,
        "amounts_mol": [float(phase.concentrations[phase.species_index(n)]
            * reactor.volume * 1000.0) for n in NAMES],
    })
print(json.dumps({"cantera": ct.__version__, "source_commit": SOURCE,
    "source_blob": BLOB, "species": NAMES, "cases": cases}, indent=2))
