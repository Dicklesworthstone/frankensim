"""Development-only oracle for the complete sourced GRI-Mech 3.0 mechanism.

Run in an isolated Python environment with Cantera 3.2.0:
    python3 data/reference/thermochem/gri30_cantera32.py [gri30.yaml]
The optional path must contain the exact pinned release source. Otherwise use
the installed release's data file. Writes JSON to stdout; never rewrites Rust
data. Neither Cantera nor Python enters FrankenSim's production runtime.

The original-data source checks exercise every species and all 325 channels.
The constant-volume reactors independently reproduce the explicitly selected
ContinuousAt300K caloric mode, including the three non-1000 K NASA seams.
"""
import hashlib
import json
import math
from pathlib import Path
import sys

import cantera as ct

assert ct.__version__ == "3.2.0"
SOURCE = "4a8358eb80cfeb50474386b5f9ec0b3a83519889"
BLOB = "9a62004ec17bbec8ac11e97785932d426dd58610"
path = Path(sys.argv[1]) if len(sys.argv) > 1 else next(
    Path(d) / "gri30.yaml" for d in ct.get_data_directories()
    if (Path(d) / "gri30.yaml").is_file()
)
raw = path.read_bytes()
assert hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest() == BLOB


def mechanism(continuous):
    gas = ct.Solution(str(path), "gri30")
    if continuous:
        # Preserve each low region and Cp exactly; change only the upper h/s
        # integration constants. The same elemental reference remains in use.
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


gas = mechanism(False)
names = sorted(gas.species_names)
indices = [gas.species_index(n) for n in names]
elements = sorted(gas.element_names)
concentrations = [0.01 * (i + 1) for i in range(len(names))]
source_checks = []
for temperature in [700.0, 1800.0]:
    gas.TPX = temperature, sum(concentrations) * ct.gas_constant * temperature / 1000.0, dict(zip(names, concentrations))
    source_checks.append({
        "temperature_k": temperature,
        "concentrations_mol_m3": concentrations,
        "thermodynamics_cp_h_s": [[
            float(gas.species(i).thermo.cp(temperature) / 1000.0),
            float(gas.species(i).thermo.h(temperature) / 1000.0),
            float(gas.species(i).thermo.s(temperature) / 1000.0),
        ] for i in indices],
        "net_rates_mol_m3_s": (gas.net_rates_of_progress * 1000.0).tolist(),
        "forward_rates_mol_m3_s": (gas.forward_rates_of_progress * 1000.0).tolist(),
        "reverse_rates_mol_m3_s": (gas.reverse_rates_of_progress * 1000.0).tolist(),
    })

cases = []
# Lean methane/air (phi=0.5), without radical seeds: induction, ignition and
# postignition nitrogen/carbon chemistry. All trajectories stay in 300--3000 K.
for temperature, duration in [(1300.0, 2e-3), (1300.0, 1e-2), (1500.0, 2e-3)]:
    gas = mechanism(True)
    gas.TPX = temperature, 101325.0, {"CH4": 1.0, "O2": 4.0, "N2": 15.04}
    reactor = ct.IdealGasReactor(gas, energy="on", volume=0.001, clone=True)
    target = reactor.mass * reactor.phase.int_energy_mass
    initial_entropy = reactor.mass * reactor.phase.entropy_mass
    initial_mass = reactor.mass
    network = ct.ReactorNet([reactor])
    network.rtol, network.atol = 1e-11, 1e-24
    network.max_steps, network.max_time_step = 200000, duration / 100.0
    network.advance(duration)
    phase = reactor.phase
    cases.append({
        "initial_temperature_k": temperature, "duration_s": duration,
        "temperature_k": phase.T, "pressure_pa": phase.P,
        "target_energy_j": target,
        "energy_residual_j": reactor.mass * phase.int_energy_mass - target,
        "entropy_change_j_k": reactor.mass * phase.entropy_mass - initial_entropy,
        "initial_mass_kg": initial_mass,
        "amounts_mol": [float(phase.concentrations[i] * reactor.volume * 1000.0) for i in indices],
    })
print(json.dumps({"cantera": ct.__version__, "source_commit": SOURCE,
    "source_blob": BLOB, "species": names, "elements": elements,
    "molar_masses_kg_mol": [float(gas.molecular_weights[i] / 1000.0) for i in indices],
    "element_counts": [[gas.n_atoms(i, element) for i in indices] for element in elements],
    "source_checks": source_checks, "cases": cases}, indent=2))
