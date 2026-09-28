//! Seed rows for the workspace registry (bead
//! frankensim-ext-benchmark-vv-registry-f1gq).
//!
//! G1 analytic entries pin an authored canonical spec as the deck: the
//! exact parameterization, the closed form where it is beyond doubt, and a
//! derivation-required discipline where quoting a formula from memory
//! would itself be the mnemonic trap the registry exists to close. G2
//! benchmark entries are seeded as known targets in a deliberately
//! UNPINNED state: they refuse citation until exact edition, license, deck
//! hash, QoIs, and acceptance envelopes are pinned.

use crate::{
    AcceptanceEnvelope, DeckPin, Edition, LicenseState, OracleBinding, PrimaryReference, Qoi,
    RegistryEntry, RegistryTier,
};

/// Edition pin shared by all authored G1 specs in this seed.
const AUTHORED_V1: Edition = Edition::Exact {
    version: "fs-vvreg authored spec v1",
};

/// License for authored spec text (part of this repository).
const REPO_LICENSE: LicenseState = LicenseState::Spdx {
    id: "LicenseRef-FrankenSim-MIT-with-OpenAI-Anthropic-Rider",
};

const fn tol(atol: f64, rtol: f64) -> AcceptanceEnvelope {
    AcceptanceEnvelope::Tolerance { atol, rtol }
}

pub(crate) fn entries() -> Vec<RegistryEntry> {
    let mut rows = Vec::new();
    rows.extend_from_slice(G1_CONTACT_IMPACT);
    rows.extend_from_slice(G1_MECHANISMS);
    rows.extend_from_slice(G1_ELECTROMAGNETICS);
    rows.extend_from_slice(G1_THERMO_FLOW);
    rows.extend_from_slice(G2_UNPINNED);
    rows
}

pub(crate) fn references() -> Vec<PrimaryReference> {
    PRIMARY_REFERENCES.to_vec()
}

const G1_CONTACT_IMPACT: &[RegistryEntry] = &[
    RegistryEntry {
        id: "g1-hertz-sphere-plane",
        tier: RegistryTier::G1Analytic,
        family: "Hertz contact",
        title: "Hertz normal contact: elastic sphere on a rigid-backed elastic half-space",
        edition: AUTHORED_V1,
        source: "Hertz (1882), J. reine angew. Math. 92; standard contact-mechanics closed form",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: sphere radius R, normal load F, effective modulus \
E* with 1/E* = (1-nu1^2)/E1 + (1-nu2^2)/E2.\n\
CLOSED FORM: contact radius a = (3 F R / (4 E*))^(1/3); peak pressure \
p0 = 3 F / (2 pi a^2); normal approach delta = a^2 / R.\n\
ASSUMPTIONS: frictionless, linear elastic, small strain, half-space \
geometry, a << R.\n\
ORACLE DISCIPLINE: evaluate the closed form directly from the pinned \
inputs; no fitted constants.",
        },
        qois: &[
            Qoi {
                name: "contact_radius",
                unit: "m",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "peak_pressure",
                unit: "Pa",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "normal_approach",
                unit: "m",
                envelope: tol(0.0, 1e-10),
            },
        ],
        notes: "Valid only inside the Hertz assumption set; no friction, adhesion, \
plasticity, or finite-thickness claims.",
    },
    RegistryEntry {
        id: "g1-hertz-cylinder-plane",
        tier: RegistryTier::G1Analytic,
        family: "Hertz contact",
        title: "Hertz line contact: elastic cylinder on an elastic half-space",
        edition: AUTHORED_V1,
        source: "Hertz (1882); standard line-contact closed form",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: cylinder radius R, length L, normal load F, \
effective modulus E* as in the sphere spec.\n\
CLOSED FORM: contact half-width b = sqrt(4 F R / (pi L E*)); peak \
pressure p0 = 2 F / (pi b L).\n\
ASSUMPTIONS: frictionless, linear elastic, plane strain, half-space, \
b << R, end effects excluded.\n\
ORACLE DISCIPLINE: the rigid-body approach of a 2-D line contact \
depends on a declared far-field datum and is NOT a registry QoI.",
        },
        qois: &[
            Qoi {
                name: "contact_half_width",
                unit: "m",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "peak_pressure",
                unit: "Pa",
                envelope: tol(0.0, 1e-10),
            },
        ],
        notes: "Approach/compliance is excluded on purpose: 2-D line-contact approach \
has no datum-free closed form.",
    },
    RegistryEntry {
        id: "g1-bouncing-ball-impact-map",
        tier: RegistryTier::G1Analytic,
        family: "impact maps",
        title: "Bouncing ball with constant restitution: impact-time map and Zeno accumulation",
        edition: AUTHORED_V1,
        source: "classical impact-map derivation (Newtonian restitution)",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: drop height h0, gravity g, restitution \
coefficient e in (0,1); release from rest.\n\
CLOSED FORM: first impact t1 = sqrt(2 h0 / g); k-th post-release impact \
time t_k = t1 * (1 + 2 e (1 - e^(k-1)) / (1 - e)); accumulation (Zeno) \
time t_inf = t1 * (1 + e) / (1 - e).\n\
ASSUMPTIONS: instantaneous impacts, constant e, no drag, planar \
vertical motion.\n\
ORACLE DISCIPLINE: event-driven simulators must hit t_k within the \
envelope for k = 3 and converge to t_inf.",
        },
        qois: &[
            Qoi {
                name: "impact_time_3",
                unit: "s",
                envelope: tol(0.0, 1e-12),
            },
            Qoi {
                name: "accumulation_time",
                unit: "s",
                envelope: tol(0.0, 1e-12),
            },
        ],
        notes: "A Zeno-accumulation fixture: integrators that step past the \
accumulation without event handling fail the map, not the ODE.",
    },
    RegistryEntry {
        id: "g1-block-incline-stick-slip",
        tier: RegistryTier::G1Analytic,
        family: "Coulomb friction",
        title: "Rigid block on an incline: stick/slip onset and sliding acceleration",
        edition: AUTHORED_V1,
        source: "classical Coulomb-friction statics",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: incline angle theta, static coefficient mu_s, \
kinetic coefficient mu_k <= mu_s, gravity g.\n\
CLOSED FORM: slip onset at tan(theta) = mu_s for the ideal rigid block; \
sliding acceleration a = g (sin(theta) - mu_k cos(theta)) for \
theta above onset.\n\
ASSUMPTIONS: ideal rigid block, uniform Coulomb friction, no tipping \
(geometry pinned so sliding precedes tipping).\n\
ORACLE DISCIPLINE: the onset angle is an exact algebraic statement; \
solvers must reproduce the stick set and the slip acceleration.",
        },
        qois: &[
            Qoi {
                name: "critical_angle",
                unit: "rad",
                envelope: tol(1e-12, 0.0),
            },
            Qoi {
                name: "sliding_acceleration",
                unit: "m/s^2",
                envelope: tol(0.0, 1e-12),
            },
        ],
        notes: "The ideal-rigid-block boundary is load-bearing: compliant or tipping \
variants are different entries, not tolerance relaxations.",
    },
];

const G1_MECHANISMS: &[RegistryEntry] = &[
    RegistryEntry {
        id: "g1-bennett-linkage-mobility",
        tier: RegistryTier::G1Analytic,
        family: "overconstrained linkages",
        title: "Bennett linkage: exactly parameterized family with mobility one",
        edition: AUTHORED_V1,
        source: "Bennett (1903), Engineering 76; overconstrained 4R spatial linkage family",
        license: REPO_LICENSE,
        oracle: OracleBinding::DerivationRequired {
            obligation: "the symbolic closure/nullity oracle for the loop-closure \
Jacobian must be pinned; numeric rank thresholding alone is refused",
        },
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: skew 4R loop with link lengths a, b and twists \
alpha, beta satisfying the Bennett conditions: opposite links equal \
(lengths and twists) and a / sin(alpha) = b / sin(beta).\n\
ORACLE: mobility DOF = 1 for the exact family even though the \
Chebychev-Grubler-Kutzbach count gives -2; the dependent-constraint \
nullity must be established by a symbolic closure/nullity oracle on \
the loop-closure Jacobian, never by numeric rank thresholding alone.\n\
ASSUMPTIONS: exact parameter satisfaction; perturbed parameters leave \
the family and lock the loop.\n\
ORACLE DISCIPLINE: closure residual of the assembled configuration \
must vanish to round-off.",
        },
        qois: &[
            Qoi {
                name: "mobility_dof",
                unit: "1",
                envelope: AcceptanceEnvelope::Interval { lo: 1.0, hi: 1.0 },
            },
            Qoi {
                name: "closure_residual",
                unit: "1",
                envelope: tol(1e-12, 0.0),
            },
        ],
        notes: "Mobility is a property of the exact family: no claim for perturbed \
parameters, where the correct answer is a locked structure.",
    },
    RegistryEntry {
        id: "g1-involute-constant-ratio",
        tier: RegistryTier::G1Analytic,
        family: "gear kinematics",
        title: "Involute gear pair: constant transmission ratio under declared assumptions",
        edition: AUTHORED_V1,
        source: "classical involute conjugate-action derivation",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: base radii rb1, rb2, declared operating center \
distance, rigid ideal involute profiles in continuous external mesh.\n\
CLOSED FORM: the transmission-ratio MAGNITUDE |omega1 / omega2| = \
rb2 / rb1, constant and center-distance-insensitive for ideal \
involutes; the SIGNED external-mesh ratio is -rb2 / rb1 (opposite \
rotation senses); kinematic transmission error (angular deviation \
from ideal conjugate motion) is identically zero.\n\
ASSUMPTIONS: rigid profiles, external mesh, no manufacturing \
deviation, declared center-distance and single/double-tooth contact \
assumptions, no backlash reversal events.\n\
ORACLE DISCIPLINE: the ratio QoI is the UNSIGNED magnitude; sign \
convention is fixed by the external-mesh declaration above. \
Transmission error is an angular-position QoI, not an interval width \
(see the gear-TE primary reference).",
        },
        qois: &[
            Qoi {
                name: "transmission_ratio_magnitude",
                unit: "1",
                envelope: tol(0.0, 1e-12),
            },
            Qoi {
                name: "transmission_error",
                unit: "rad",
                envelope: tol(1e-12, 0.0),
            },
        ],
        notes: "Ideal-profile kinematics only: loaded, compliant, or modified-profile \
transmission error belongs to G2 decks with pinned data.",
    },
    RegistryEntry {
        id: "g1-geneva-closed-form",
        tier: RegistryTier::G1Analytic,
        family: "intermittent mechanisms",
        title: "External Geneva mechanism: closed-form wheel kinematics",
        edition: AUTHORED_V1,
        source: "classical external-Geneva kinematic derivation",
        license: REPO_LICENSE,
        oracle: OracleBinding::DerivationRequired {
            obligation: "closed-form wheel kinematics must be derived from the pinned \
tangency geometry and symbolically checked (zero wheel velocity at slot \
entry/exit)",
        },
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: n slots (n >= 3), tangency-condition geometry: \
drive crank radius r = C sin(pi / n) with center distance C; drive \
angle measured from slot-entry tangency.\n\
ORACLE: the wheel angle, angular velocity ratio, and angular \
acceleration ratio versus drive angle have closed forms that MUST be \
derived from the pinned geometry and symbolically checked (entry and \
exit tangency give zero wheel velocity); they are not to be quoted \
from memory into test code.\n\
ASSUMPTIONS: rigid links, ideal pin/slot contact, purely kinematic \
(no dynamics, no clearances).\n\
ORACLE DISCIPLINE: peak ratios are evaluated on the derived closed \
form; the envelope binds candidate kinematics against that derivation.",
        },
        qois: &[
            Qoi {
                name: "max_wheel_velocity_ratio",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "max_wheel_acceleration_ratio",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
        ],
        notes: "Derivation-required entry: the deck pins geometry and boundary \
conditions, not a memorized formula.",
    },
    RegistryEntry {
        id: "g1-epitrochoid-closed-form",
        tier: RegistryTier::G1Analytic,
        family: "trochoid geometry",
        title: "Epitrochoid curve: exact parameterization for rotary-housing geometry",
        edition: AUTHORED_V1,
        source: "classical roulette derivation; see the Wankel-geometry primary reference",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: fixed circle radius R, rolling circle radius r \
(external rolling), generating point at distance d from the rolling \
center.\n\
CLOSED FORM: x(t) = (R + r) cos(t) - d cos(((R + r) / r) t); \
y(t) = (R + r) sin(t) - d sin(((R + r) / r) t).\n\
ASSUMPTIONS: pure rolling without slip; housing-form specializations \
(e.g. two-lobe Wankel bore) fix R / r and d and MUST be re-derived and \
symbolically checked against this parameterization, never quoted \
mnemonically.\n\
ORACLE DISCIPLINE: candidate curve generators are checked pointwise \
against the closed form over a pinned parameter grid.",
        },
        qois: &[Qoi {
            name: "curve_point_max_deviation",
            unit: "m",
            envelope: tol(1e-12, 0.0),
        }],
        notes: "Geometry only: no claim about seal kinematics, chamber volume, or \
combustion behavior.",
    },
];

const G1_ELECTROMAGNETICS: &[RegistryEntry] = &[
    RegistryEntry {
        id: "g1-coaxial-cable-fields",
        tier: RegistryTier::G1Analytic,
        family: "canonical electrostatics",
        title: "Coaxial cable: per-length capacitance and inductance",
        edition: AUTHORED_V1,
        source: "classical TEM transmission-line closed forms",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: inner radius a, outer radius b > a, uniform \
permittivity eps and permeability mu between conductors.\n\
CLOSED FORM: C' = 2 pi eps / ln(b / a); L' = mu ln(b / a) / (2 pi).\n\
ASSUMPTIONS: perfect conductors, static/TEM fields, infinite length \
(no end effects), homogeneous dielectric.\n\
ORACLE DISCIPLINE: field solvers must recover both per-length \
constants from the pinned geometry.",
        },
        qois: &[
            Qoi {
                name: "capacitance_per_length",
                unit: "F/m",
                envelope: tol(0.0, 1e-12),
            },
            Qoi {
                name: "inductance_per_length",
                unit: "H/m",
                envelope: tol(0.0, 1e-12),
            },
        ],
        notes: "No skin-effect, loss, or finite-conductivity claims.",
    },
    RegistryEntry {
        id: "g1-sphere-uniform-field",
        tier: RegistryTier::G1Analytic,
        family: "canonical electrostatics",
        title: "Conducting sphere in a uniform field: enhancement and induced dipole",
        edition: AUTHORED_V1,
        source: "classical separation-of-variables electrostatics",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: grounded conducting sphere radius a in a uniform \
applied field E0 in vacuum (permittivity eps0).\n\
CLOSED FORM: surface field E(theta) = 3 E0 cos(theta) (peak \
enhancement factor exactly 3 at the poles); induced dipole moment \
p = 4 pi eps0 a^3 E0.\n\
ASSUMPTIONS: electrostatics, isolated sphere, uniform far field.\n\
ORACLE DISCIPLINE: the enhancement factor is dimensionless and exact; \
mesh-converged solvers are checked against it and the dipole moment.",
        },
        qois: &[
            Qoi {
                name: "peak_field_enhancement",
                unit: "1",
                envelope: tol(0.0, 1e-12),
            },
            Qoi {
                name: "induced_dipole_moment",
                unit: "C*m",
                envelope: tol(0.0, 1e-12),
            },
        ],
        notes: "Conductor case only; the dielectric-sphere interior-field variant is \
a separate future entry.",
    },
    RegistryEntry {
        id: "g1-helmholtz-coil",
        tier: RegistryTier::G1Analytic,
        family: "canonical magnetostatics",
        title: "Helmholtz coil pair: center field and axial uniformity",
        edition: AUTHORED_V1,
        source: "classical Biot-Savart loop superposition",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: two coaxial filamentary loops of radius R, \
separation R, equal currents I in the same sense, vacuum \
permeability mu0.\n\
CLOSED FORM: center field B = (4/5)^(3/2) mu0 I / R; the first and \
second axial derivatives of B vanish at the midpoint (the Helmholtz \
condition).\n\
ASSUMPTIONS: filamentary loops, magnetostatics, exact spacing.\n\
ORACLE DISCIPLINE: solvers must recover the center field and the \
vanishing second derivative.",
        },
        qois: &[
            Qoi {
                name: "center_field",
                unit: "T",
                envelope: tol(0.0, 1e-12),
            },
            Qoi {
                name: "axial_second_derivative_at_center",
                unit: "T/m^2",
                envelope: tol(1e-12, 0.0),
            },
        ],
        notes: "Filamentary idealization: finite winding cross-sections are outside \
this entry.",
    },
];

const G1_THERMO_FLOW: &[RegistryEntry] = &[
    RegistryEntry {
        id: "g1-stefan-problem",
        tier: RegistryTier::G1Analytic,
        family: "moving boundaries",
        title: "One-phase Stefan problem: similarity interface motion",
        edition: AUTHORED_V1,
        source: "classical Neumann similarity solution",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: semi-infinite solid at the melt temperature, \
wall step temperature above melting, thermal diffusivity alpha, \
Stefan number St = c dT / L.\n\
CLOSED FORM: interface position s(t) = 2 lambda sqrt(alpha t) where \
lambda solves the transcendental relation \
sqrt(pi) lambda exp(lambda^2) erf(lambda) = St.\n\
ASSUMPTIONS: one-phase, constant properties, planar front, no \
convection or density change.\n\
ORACLE DISCIPLINE: lambda is obtained by a bracketed root solve of \
the pinned transcendental relation; front-tracking or enthalpy \
solvers are checked against s(t).",
        },
        qois: &[
            Qoi {
                name: "interface_coefficient_lambda",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "interface_position",
                unit: "m",
                envelope: tol(0.0, 1e-10),
            },
        ],
        notes: "One-phase idealization; two-phase and density-change variants are \
separate entries.",
    },
    RegistryEntry {
        id: "g1-riemann-sod",
        tier: RegistryTier::G1Analytic,
        family: "exact Riemann problems",
        title: "Sod shock tube: exact Riemann solution states",
        edition: AUTHORED_V1,
        source: "Sod (1978), J. Comput. Phys. 27; exact Riemann solver oracle",
        license: REPO_LICENSE,
        oracle: OracleBinding::DerivationRequired {
            obligation: "the exact-solver construction (pressure-function root \
bracketing, wave-pattern branch relations) and the promised (x, t) \
comparison sample set must be pinned as executable deck data",
        },
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: ideal gas gamma = 1.4; left state \
(rho, u, p) = (1, 0, 1); right state (0.125, 0, 0.1); dimensionless \
units.\n\
ORACLE: the exact self-similar solution (left rarefaction, contact, \
right shock) computed by a convergent exact Riemann solver: star \
pressure and velocity from the pressure-function root, shock speed \
from Rankine-Hugoniot.\n\
ASSUMPTIONS: calorically perfect gas, 1-D, no viscosity.\n\
ORACLE DISCIPLINE: the oracle is the exact solver, not tabulated \
digits; candidate codes are compared at pinned (x, t) sample points \
and on the star-region constants.",
        },
        qois: &[
            Qoi {
                name: "star_pressure",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "star_velocity",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "shock_speed",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
        ],
        notes: "Exact-solution comparison only; scheme convergence-order claims are a \
different gate.",
    },
    RegistryEntry {
        id: "g1-riemann-lax",
        tier: RegistryTier::G1Analytic,
        family: "exact Riemann problems",
        title: "Lax problem: exact Riemann solution states",
        edition: AUTHORED_V1,
        source: "Lax (1954) initial data as standardized in shock-tube test suites",
        license: REPO_LICENSE,
        oracle: OracleBinding::DerivationRequired {
            obligation: "the exact-solver construction and the pinned (x, t) \
comparison sample set must be pinned as executable deck data, as in \
the Sod entry",
        },
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: ideal gas gamma = 1.4; left state \
(rho, u, p) = (0.445, 0.698, 3.528); right state (0.5, 0, 0.571); \
dimensionless units.\n\
ORACLE: exact Riemann solution via the same convergent exact solver \
discipline as the Sod entry.\n\
ASSUMPTIONS: calorically perfect gas, 1-D, no viscosity.\n\
ORACLE DISCIPLINE: star-region constants and pinned (x, t) samples.",
        },
        qois: &[
            Qoi {
                name: "star_pressure",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "star_velocity",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
            Qoi {
                name: "shock_speed",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
        ],
        notes: "Stronger-wave companion to Sod; same exact-solver oracle discipline.",
    },
    RegistryEntry {
        id: "g1-isentropic-nozzle",
        tier: RegistryTier::G1Analytic,
        family: "compressible flow",
        title: "Quasi-1D isentropic nozzle: area-Mach relation and choking",
        edition: AUTHORED_V1,
        source: "classical gas-dynamics area-Mach relation",
        license: REPO_LICENSE,
        oracle: OracleBinding::DerivationRequired {
            obligation: "the area law A(x), exit/throat area ratio, and subsonic-vs-\
supersonic branch selection are delegated to the consuming deck and \
must be pinned there before citation",
        },
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: ideal gas gamma = 1.4, converging-diverging area \
law A(x) pinned by the consuming deck, exit/throat area ratio given, \
supersonic branch declared downstream of the throat.\n\
CLOSED FORM: A / A* = (1 / M) * ((2 / (gamma + 1)) * \
(1 + ((gamma - 1) / 2) M^2))^((gamma + 1) / (2 (gamma - 1))); choked \
flow has M = 1 exactly at the throat.\n\
ASSUMPTIONS: quasi-1D, isentropic, calorically perfect gas, declared \
branch selection (subsonic vs supersonic root).\n\
ORACLE DISCIPLINE: exit Mach is the declared-branch root of the pinned \
relation, obtained by a bracketed solve.",
        },
        qois: &[
            Qoi {
                name: "throat_mach",
                unit: "1",
                envelope: tol(1e-12, 0.0),
            },
            Qoi {
                name: "supersonic_exit_mach",
                unit: "1",
                envelope: tol(0.0, 1e-10),
            },
        ],
        notes: "Branch selection is part of the deck: quoting an exit Mach without \
the declared branch is not a citation.",
    },
    RegistryEntry {
        id: "g1-otto-cycle",
        tier: RegistryTier::G1Analytic,
        family: "air-standard cycles",
        title: "Air-standard Otto cycle: thermal efficiency closed form",
        edition: AUTHORED_V1,
        source: "classical air-standard cycle analysis",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: compression ratio r > 1, ratio of specific heats \
gamma; legs: isentropic compression, isochoric heat addition, \
isentropic expansion, isochoric rejection.\n\
CLOSED FORM: eta = 1 - r^(1 - gamma).\n\
ASSUMPTIONS: air-standard (ideal gas, constant properties, reversible \
legs, no combustion chemistry).\n\
ORACLE DISCIPLINE: efficiency evaluated from the closed form at pinned \
(r, gamma).",
        },
        qois: &[Qoi {
            name: "thermal_efficiency",
            unit: "1",
            envelope: tol(0.0, 1e-12),
        }],
        notes: "Air-standard idealization only; no claim about real engine cycles \
(see the CFR G2 entry for measured traces).",
    },
    RegistryEntry {
        id: "g1-diesel-cycle",
        tier: RegistryTier::G1Analytic,
        family: "air-standard cycles",
        title: "Air-standard Diesel cycle: thermal efficiency closed form",
        edition: AUTHORED_V1,
        source: "classical air-standard cycle analysis",
        license: REPO_LICENSE,
        oracle: OracleBinding::SelfContained,
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: compression ratio r > 1, cutoff ratio rho > 1, \
ratio of specific heats gamma; legs: isentropic compression, isobaric \
heat addition to the cutoff, isentropic expansion, isochoric \
rejection.\n\
CLOSED FORM: eta = 1 - r^(1 - gamma) * (rho^gamma - 1) / \
(gamma (rho - 1)).\n\
ASSUMPTIONS: air-standard as in the Otto entry.\n\
ORACLE DISCIPLINE: efficiency evaluated from the closed form at pinned \
(r, rho, gamma); the Otto limit rho -> 1 must be recovered.",
        },
        qois: &[Qoi {
            name: "thermal_efficiency",
            unit: "1",
            envelope: tol(0.0, 1e-12),
        }],
        notes: "The rho -> 1 Otto limit is a mandatory internal cross-check of the \
consuming test.",
    },
    RegistryEntry {
        id: "g1-atkinson-cycle",
        tier: RegistryTier::G1Analytic,
        family: "air-standard cycles",
        title: "Air-standard Atkinson cycle: derivation-required efficiency",
        edition: AUTHORED_V1,
        source: "classical air-standard cycle analysis (full-expansion Atkinson)",
        license: REPO_LICENSE,
        oracle: OracleBinding::DerivationRequired {
            obligation: "closed-form efficiency must be derived from the pinned cycle \
legs and checked by first-law leg balances internal to the declared \
topology",
        },
        deck: DeckPin::AuthoredSpec {
            spec: "PARAMETERIZATION: compression ratio r, expansion ratio e > r, \
ratio of specific heats gamma; legs: isentropic compression (ratio r), \
isochoric heat addition, isentropic expansion (ratio e) to the initial \
pressure, isobaric heat rejection closing the cycle.\n\
ORACLE: the closed-form efficiency MUST be derived from the pinned \
legs and symbolically checked (first-law leg balances summing to the \
net work; state continuity around the closed cycle), never quoted \
from memory; published Atkinson formulas vary with parameterization \
and are a known mnemonic trap.\n\
ASSUMPTIONS: air-standard as in the Otto entry; e strictly greater \
than r (at e = r the pinned isobaric-rejection topology degenerates \
and does NOT reduce to the constant-volume-rejection Otto cycle).\n\
ORACLE DISCIPLINE: checks are internal to the declared topology: \
per-leg first-law balances, cycle closure of the state point, and \
efficiency strictly inside (0, 1) for valid (r, e, gamma).",
        },
        qois: &[Qoi {
            name: "thermal_efficiency",
            unit: "1",
            envelope: tol(0.0, 1e-12),
        }],
        notes: "Derivation-required entry: the deck pins the cycle legs and the \
internal first-law/closure checks, not a formula. No Otto-limit claim: \
the pinned topology degenerates at e = r.",
    },
];

/// G2 rows seeded as known targets, deliberately unpinned: exact edition,
/// license, deck hash, and acceptance data must be added before citation.
const G2_UNPINNED: &[RegistryEntry] = &[
    RegistryEntry {
        id: "g2-team-10",
        tier: RegistryTier::G2Benchmark,
        family: "TEAM",
        title: "TEAM problem 10: steel plates around a coil, nonlinear transient eddy current",
        edition: Edition::Unpinned,
        source: "COMPUMAG TEAM benchmark suite, official problem definition 10",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "average_flux_density_probe",
            unit: "T",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Family name only: exact geometry revision, excitation, material law, circuit, QoI set, and acceptance data must be pinned before citation.",
    },
    RegistryEntry {
        id: "g2-team-13",
        tier: RegistryTier::G2Benchmark,
        family: "TEAM",
        title: "TEAM problem 13: nonlinear magnetostatic model with thin steel channels",
        edition: Edition::Unpinned,
        source: "COMPUMAG TEAM benchmark suite, official problem definition 13",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "flux_density_probe_set",
            unit: "T",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "The problem number alone is not a deck: pin geometry revision, B-H \
curve, excitation, QoIs, and acceptance data.",
    },
    RegistryEntry {
        id: "g2-team-20",
        tier: RegistryTier::G2Benchmark,
        family: "TEAM",
        title: "TEAM problem 20: 3-D static force on a pole piece",
        edition: Edition::Unpinned,
        source: "COMPUMAG TEAM benchmark suite, official problem definition 20",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "pole_force",
            unit: "N",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Pin geometry revision, coil excitation levels, material law, force \
QoI definition, and acceptance data before citation.",
    },
    RegistryEntry {
        id: "g2-team-24",
        tier: RegistryTier::G2Benchmark,
        family: "TEAM",
        title: "TEAM problem 24: nonlinear time-transient rotational device",
        edition: Edition::Unpinned,
        source: "COMPUMAG TEAM benchmark suite, official problem definition 24",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "torque_transient",
            unit: "N*m",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Pin geometry revision, drive circuit, material law, measured \
transient QoIs, and acceptance data before citation.",
    },
    RegistryEntry {
        id: "g2-team-30a",
        tier: RegistryTier::G2Benchmark,
        family: "TEAM",
        title: "TEAM problem 30a: induction motor analysis (single-phase variant)",
        edition: Edition::Unpinned,
        source: "COMPUMAG TEAM benchmark suite, official problem definition 30a",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "torque_versus_slip",
            unit: "N*m",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Pin geometry, winding excitation, conductivity set, slip sweep, and \
acceptance data before citation.",
    },
    RegistryEntry {
        id: "g2-team-30b",
        tier: RegistryTier::G2Benchmark,
        family: "TEAM",
        title: "TEAM problem 30b: induction motor analysis (three-phase variant)",
        edition: Edition::Unpinned,
        source: "COMPUMAG TEAM benchmark suite, official problem definition 30b",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "torque_versus_slip",
            unit: "N*m",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Pin geometry, winding excitation, conductivity set, slip sweep, and \
acceptance data before citation.",
    },
    RegistryEntry {
        id: "g2-iftomm-rectangular-bricard",
        tier: RegistryTier::G2Benchmark,
        family: "IFToMM",
        title: "IFToMM rectangular Bricard linkage: redundantly constrained multibody benchmark",
        edition: Edition::Unpinned,
        source: "IFToMM Library of Computational Benchmark Problems",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "joint_trajectory_reference",
            unit: "m",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Redundantly constrained with a configuration-dependent \
dependent-constraint subset: the pinned deck must include the exact \
input/result artifacts, not just the mechanism name.",
    },
    RegistryEntry {
        id: "g2-iftomm-andrews-squeezer",
        tier: RegistryTier::G2Benchmark,
        family: "IFToMM",
        title: "IFToMM Andrews squeezer mechanism: stiff multibody dynamics benchmark",
        edition: Edition::Unpinned,
        source: "IFToMM Library of Computational Benchmark Problems",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "angle_trajectory_reference",
            unit: "rad",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Pin the exact benchmark input/result artifact version and tolerance \
bands before citation.",
    },
    RegistryEntry {
        id: "g2-nafems-thermal-set",
        tier: RegistryTier::G2Benchmark,
        family: "NAFEMS",
        title: "NAFEMS thermal benchmark set",
        edition: Edition::Unpinned,
        source: "NAFEMS benchmark index (exact case IDs and reports to be pinned)",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "case_temperature_probe",
            unit: "K",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "NAFEMS artifacts are licensed: pin exact case ID, report edition, \
license terms, and storage location before citation.",
    },
    RegistryEntry {
        id: "g2-nafems-r0083-acoustic",
        tier: RegistryTier::G2Benchmark,
        family: "NAFEMS",
        title: "NAFEMS R0083 acoustic benchmark",
        edition: Edition::Unpinned,
        source: "NAFEMS report R0083 (exact edition and case selection to be pinned)",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "acoustic_pressure_probe",
            unit: "Pa",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Licensed report: pin edition, case, license, deck hash, QoIs, and \
acceptance data before citation.",
    },
    RegistryEntry {
        id: "g2-cfr-engine-pressure-traces",
        tier: RegistryTier::G2Benchmark,
        family: "CFR",
        title: "CFR engine pressure traces and motored p-V loops",
        edition: Edition::Unpinned,
        source: "CFR cooperative fuel research engine datasets (versioned/licensed source to be pinned)",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[
            Qoi {
                name: "cylinder_pressure_trace",
                unit: "Pa",
                envelope: AcceptanceEnvelope::Unpinned,
            },
            Qoi {
                name: "motored_pv_loop_work",
                unit: "J",
                envelope: AcceptanceEnvelope::Unpinned,
            },
        ],
        notes: "Measured-data entry: pin dataset version, license, measurement \
uncertainty, and acceptance envelopes before citation.",
    },
    RegistryEntry {
        id: "g2-nasa-caa-benchmarks",
        tier: RegistryTier::G2Benchmark,
        family: "NASA CAA",
        title: "NASA computational aeroacoustics benchmark decks",
        edition: Edition::Unpinned,
        source: "NASA CAA workshop benchmark proceedings (exact deck selection to be pinned)",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "acoustic_waveform_probe",
            unit: "Pa",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Pin the exact workshop volume, case deck, QoIs, and acceptance data \
before citation.",
    },
    RegistryEntry {
        id: "g2-sandia-ecn-spray",
        tier: RegistryTier::G2Benchmark,
        family: "ECN",
        title: "Sandia Engine Combustion Network spray benchmark",
        edition: Edition::Unpinned,
        source: "Engine Combustion Network data archive (exact configuration to be pinned)",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[
            Qoi {
                name: "liquid_penetration_length",
                unit: "m",
                envelope: AcceptanceEnvelope::Unpinned,
            },
            Qoi {
                name: "ignition_delay",
                unit: "s",
                envelope: AcceptanceEnvelope::Unpinned,
            },
        ],
        notes: "Pin exact configuration (e.g. Spray A revision), diagnostic \
uncertainty, and blind/held-out QoI partitions before citation.",
    },
    RegistryEntry {
        id: "g2-oscillating-cylinder-flow",
        tier: RegistryTier::G2Benchmark,
        family: "moving-boundary flow",
        title: "Oscillating-cylinder flow: moving-boundary benchmark deck",
        edition: Edition::Unpinned,
        source: "canonical oscillating-cylinder studies (exact dataset to be pinned)",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[
            Qoi {
                name: "force_coefficient_history",
                unit: "1",
                envelope: AcceptanceEnvelope::Unpinned,
            },
            Qoi {
                name: "shedding_frequency",
                unit: "Hz",
                envelope: AcceptanceEnvelope::Unpinned,
            },
        ],
        notes: "Pin the exact study, Reynolds/Keulegan-Carpenter parameters, data \
license, and acceptance envelopes before citation.",
    },
    RegistryEntry {
        id: "g2-gerotor-moving-boundary",
        tier: RegistryTier::G2Benchmark,
        family: "moving-boundary flow",
        title: "Gerotor pump: moving-boundary internal flow benchmark deck",
        edition: Edition::Unpinned,
        source: "gerotor flow benchmark literature (exact dataset to be pinned)",
        license: LicenseState::Unpinned,
        oracle: OracleBinding::Unpinned,
        deck: DeckPin::Unpinned,
        qois: &[Qoi {
            name: "flow_ripple",
            unit: "m^3/s",
            envelope: AcceptanceEnvelope::Unpinned,
        }],
        notes: "Pin the exact geometry (tied to the epitrochoid G1 entry's derived \
profile discipline), operating point, dataset license, and acceptance \
data before citation.",
    },
];

/// The complete primary-reference seed from the bead: anchors definitions
/// and benchmark provenance; never authority-by-citation.
const PRIMARY_REFERENCES: &[PrimaryReference] = &[
    PrimaryReference {
        index: 1,
        key: "feec-stability-afw",
        citation: "Arnold, Falk & Winther, 'Finite element exterior calculus: from Hodge theory to numerical stability'",
        locator: "arXiv:0906.4325",
        anchors: "dd=0 is insufficient; subcomplex property, bounded cochain projection, and formulation stability are explicit gates",
        boundary: "does not certify any particular discretization in this workspace without those gates being tested",
    },
    PrimaryReference {
        index: 2,
        key: "sheaf-spectra-hansen-ghrist",
        citation: "Hansen & Ghrist, 'Toward a Spectral Theory of Cellular Sheaves'",
        locator: "arXiv:1808.01513",
        anchors: "sheaf Laplacians are meaningful only after explicit stalks/restrictions are constructed",
        boundary: "L1 spectral machinery stays generic; no claim transfers without constructed restriction maps",
    },
    PrimaryReference {
        index: 3,
        key: "sheaf-cosheaf-curry",
        citation: "Curry, 'Sheaves, Cosheaves and Applications'",
        locator: "arXiv:1303.3255",
        anchors: "restriction and aggregation directions are distinct structures",
        boundary: "the power pairing and commuting diagram are NEW theorem obligations, not borrowed authority",
    },
    PrimaryReference {
        index: 4,
        key: "nasa9-thermo-mcbride",
        citation: "McBride, Zehe & Gordon, NASA Glenn coefficients for calculating thermodynamic properties of individual species",
        locator: "NASA/TP-2002-211556",
        anchors: "NASA-9 polynomials are standard-state MOLAR thermodynamics; amount-of-substance and reference metadata are load-bearing",
        boundary: "phase/EOS derivations beyond standard state must be derived and tested separately",
    },
    PrimaryReference {
        index: 5,
        key: "entropy-stable-tadmor",
        citation: "Tadmor, 'Entropy Stable Schemes'",
        locator: "doi:10.1016/bs.hna.2016.09.006",
        anchors: "entropy-conservative spatial flux construction and its sign conventions",
        boundary: "the spatial flux theorem is not the fully discrete/source/boundary theorem; the exact scheme must be pinned",
    },
    PrimaryReference {
        index: 6,
        key: "port-hamiltonian-dirac-cervera",
        citation: "Cervera, van der Schaft & Banos, 'Interconnection of port-Hamiltonian systems and composition of Dirac structures'",
        locator: "doi:10.1016/j.automatica.2006.08.014",
        anchors: "Dirac interconnection relations are power-conserving",
        boundary: "dissipation is a separate resistive structure, not part of the Dirac composition claim",
    },
    PrimaryReference {
        index: 7,
        key: "rigidity-maxwell-calladine-rocks",
        citation: "Rocks et al., 'Integrating local energetics into Maxwell-Calladine constraint counting'",
        locator: "arXiv:2208.07419",
        anchors: "mechanisms = ker J; self-stresses = ker J^T; index and energy structure are distinct",
        boundary: "constraint counting alone does not decide mobility of exact overconstrained families (see the Bennett entry)",
    },
    PrimaryReference {
        index: 8,
        key: "contact-ipc-li",
        citation: "Li et al., 'Incremental Potential Contact' and 'Convergent IPC'",
        locator: "doi:10.1145/3386569.3392425; arXiv:2307.15908",
        anchors: "nonintersection/convergence claims enumerate candidate-set, CCD, solve, and refinement assumptions",
        boundary: "pin the exact preprint version/hash; smoothness claims are regime-specific",
    },
    PrimaryReference {
        index: 9,
        key: "codim-ipc-li",
        citation: "Li et al., 'Codimensional Incremental Potential Contact'",
        locator: "arXiv:2012.04457",
        anchors: "shell/rod/seal thickness behavior under a declared codimensional model",
        boundary: "claims are conditional on the declared candidates, CCD, and accepted solve",
    },
    PrimaryReference {
        index: 10,
        key: "nonsmooth-contact-acary",
        citation: "Acary, 'Energy conservation and dissipation properties of time-integration methods for nonsmooth elastodynamics with contact'",
        locator: "arXiv:1410.2499",
        anchors: "Moreau-Jean energy behavior is law- and parameter-dependent",
        boundary: "pairwise restitution is not a global energy theorem; the exact global impact-law source must be pinned per bead",
    },
    PrimaryReference {
        index: 11,
        key: "validated-flowpipes-walawska-wilczak",
        citation: "Walawska & Wilczak, implicit validated variational-equation enclosures",
        locator: "arXiv:1509.07388",
        anchors: "true-flow enclosure for validated integration",
        boundary: "DAE regularity and finite-guard isolability are separate obligations with their own pinned sources",
    },
    PrimaryReference {
        index: 12,
        key: "gear-te-athavale",
        citation: "Athavale, Krishnaswami & Kuo, gear transmission error analysis",
        locator: "SAE 2001-01-1507",
        anchors: "transmission error is angular-position deviation from ideal conjugate motion, not an interval width",
        boundary: "an executable G2 case needs a complete reproducible deck, not the paper alone",
    },
    PrimaryReference {
        index: 13,
        key: "wankel-seals-handschuh-owen",
        citation: "Handschuh & Owen, rotary engine apex seal analysis",
        locator: "NASA/TM-2010-216353",
        anchors: "reduced seal-load/friction reference with explicit assumptions",
        boundary: "does not prove universal finite-radius housing geometry claims",
    },
    PrimaryReference {
        index: 14,
        key: "iftomm-benchmark-library",
        citation: "IFToMM Library of Computational Benchmark Problems",
        locator: "IFToMM benchmark library (exact artifact per entry)",
        anchors: "exact input/result artifacts for Andrews, Bricard, and slider-crank benchmarks",
        boundary: "a mechanism name without the exact artifact version is not citable",
    },
    PrimaryReference {
        index: 15,
        key: "team-em-benchmarks",
        citation: "COMPUMAG TEAM benchmark official problem definitions 10/13/20/24/30a/30b",
        locator: "COMPUMAG TEAM suite (exact revision per entry)",
        anchors: "official definitions of geometry, excitation, material law, circuit, QoI, and acceptance data",
        boundary: "the problem number alone is not a deck; pin the exact revision",
    },
    PrimaryReference {
        index: 16,
        key: "nafems-thermal",
        citation: "NAFEMS thermal benchmark index and guide",
        locator: "NAFEMS benchmark index (exact case ID/report/license per entry)",
        anchors: "case identity for thermal benchmarks",
        boundary: "licensed artifacts: exact case ID, report edition, and license are required before citation",
    },
    PrimaryReference {
        index: 17,
        key: "nasa9-cantera-oracle",
        citation: "Cantera 3.2 NASA9 documentation (development-only comparison oracle)",
        locator: "Cantera 3.2 release documentation",
        anchors: "oracle receipts pin release, species/mechanism hash, temperature region, and reference pressure",
        boundary: "never a production dependency; oracle receipts only",
    },
    PrimaryReference {
        index: 18,
        key: "nonholonomic-modin-verdier",
        citation: "Modin & Verdier, 'What makes nonholonomic integrators work?'",
        locator: "doi:10.1007/s00211-020-01126-y",
        anchors: "rolling constraints are Pfaffian and generally nonintegrable",
        boundary: "no ordinary Hamilton/RATTLE claims are inherited; the exact discrete Lagrange-d'Alembert source is pinned per bead",
    },
    PrimaryReference {
        index: 19,
        key: "hcurl-ams-hiptmair-xu",
        citation: "Hiptmair & Xu, auxiliary space preconditioning in H(curl), plus hypre AMS documentation",
        locator: "doi:10.1137/060660588",
        anchors: "the auxiliary-space component list for H(curl) preconditioning",
        boundary: "boundary, topology, and coefficient assumptions must be explicit per use",
    },
    PrimaryReference {
        index: 20,
        key: "switched-descriptor-yildiz",
        citation: "Yildiz, MNA-based unified ideal switch model",
        locator: "doi:10.1142/S0218126613500461",
        anchors: "inconsistent switching causes impulses and discontinuous device states",
        boundary: "impulse-free continuity is conditional, not a default assumption",
    },
    PrimaryReference {
        index: 21,
        key: "nonlinear-magnetic-energy-mandlmayr-egger",
        citation: "Mandlmayr & Egger, nonlinear magnetic field energy analysis",
        locator: "arXiv:2311.02380",
        anchors: "convex energy/coenergy structure versus merely strongly monotone nonintegrable laws",
        boundary: "energy-based claims require the convexity hypothesis to be checked, not assumed",
    },
    PrimaryReference {
        index: 22,
        key: "reacting-entropy-ching",
        citation: "Ching, Johnson & Kercher, positivity-preserving entropy-stable reacting flow, Parts I/II",
        locator: "arXiv:2211.16254; arXiv:2211.16297",
        anchors: "constructions tied to exact thermodynamics, chemistry, source terms, and discretization",
        boundary: "no transfer to a different thermodynamic or chemistry model without re-derivation",
    },
    PrimaryReference {
        index: 23,
        key: "acoustics-burton-miller-fwh",
        citation: "Burton & Miller combined-field exterior formulation; Ffowcs Williams & Hawkings moving-surface aeroacoustics; NASA CAA proceedings; NAFEMS R0083",
        locator: "doi:10.1016/0022-247X(84)90146-X; doi:10.1098/rsta.1969.0031",
        anchors: "combined-field uniqueness for exterior acoustics and moving-surface aeroacoustic analogy definitions",
        boundary: "benchmark citations require the pinned G2 decks, not the formulation papers alone",
    },
    PrimaryReference {
        index: 24,
        key: "tribology-ehl-hamrock-dowson",
        citation: "Hamrock & Dowson EHL film-thickness correlations, plus independent optical-film/traction datasets",
        locator: "NASA-TP-1342",
        anchors: "film correlations are regime-bounded engineering fits",
        boundary: "full EHL must verify balances and validate held-out data; correlations are not verification oracles",
    },
    PrimaryReference {
        index: 25,
        key: "fatigue-astm-nasgro",
        citation: "ASTM E466-21 and E647-24 test methods, plus a pinned NASA NASGRO reference artifact as a development oracle",
        locator: "ASTM E466-21; ASTM E647-24",
        anchors: "method conformance for fatigue and crack-growth testing",
        boundary: "method conformance never validates component life",
    },
    PrimaryReference {
        index: 26,
        key: "assurance-iso12100-iec60034",
        citation: "ISO 12100:2010, IEC 60034-1:2026, IEC 60034-2-1:2024",
        locator: "ISO 12100:2010; IEC 60034-1:2026; IEC 60034-2-1:2024",
        anchors: "scope, rating, and loss-test methods as conformance artifacts",
        boundary: "simulation conformance is not regulatory approval",
    },
    PrimaryReference {
        index: 27,
        key: "vvuq-asme-gum-nasa7009",
        citation: "ASME VVUQ standards registry, BIPM/JCGM GUM, NASA-STD-7009",
        locator: "ASME VVUQ; JCGM GUM; NASA-STD-7009",
        anchors: "context of use, solution verification, experimental uncertainty, calibration/validation split, and prediction assessment as separate artifacts",
        boundary: "naming the standard does not discharge any of those artifacts",
    },
    PrimaryReference {
        index: 28,
        key: "combustion-sandia-ecn",
        citation: "Sandia Engine Combustion Network data archive",
        locator: "ECN data archive (exact configuration per entry)",
        anchors: "versioned spray/combustion configurations with diagnostic uncertainty",
        boundary: "blind/held-out QoI partitions must be pinned and versioned before validation claims",
    },
    PrimaryReference {
        index: 29,
        key: "interop-fmi-ssp",
        citation: "FMI 3.0.2 and SSP 2.0 interoperability standards",
        locator: "FMI 3.0.2; SSP 2.0",
        anchors: "conformance and quarantined-adapter behavior for co-simulation interchange",
        boundary: "foreign runtimes are never admitted to the trusted graph by conformance alone",
    },
    PrimaryReference {
        index: 30,
        key: "gear-iso6336",
        citation: "ISO 6336-1:2019 with exact applicable parts, editions, and worked examples",
        locator: "ISO 6336-1:2019",
        anchors: "formula conformance for gear load capacity inside the standard's scope",
        boundary: "drive-life prediction needs independent validation; conformance is scope-bounded",
    },
];
