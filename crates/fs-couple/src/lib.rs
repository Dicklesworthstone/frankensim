//! fs-couple — multiphysics composition through port-Hamiltonian Dirac
//! structures. Layer: L3.
//!
//! Ad-hoc FSI staggering suffers added-mass instabilities and energy drift.
//! The implemented DIRAC interconnection is LOSSLESS BY CONSTRUCTION:
//! power-conjugate [`Port`]s use equal effort and opposite flow, so their net
//! interface power is exactly zero. [`EnergyAudit`] records caller-supplied
//! interface balances as a G0 bug alarm. Neither invariant alone proves that
//! the coupled components, discretizations, transfers, iterations, time
//! integrators, sources, or a finite accounting window are passive.
//!
//! [`PortSchema`] v2 is the dependency-light, versioned description carried by
//! new coupling relations. It makes identity, dimensions, shape, coordinates,
//! clock, power pairing, and conservation roles explicit. The four primitive
//! relation descriptors distinguish conservative junctions, storage,
//! dissipation, and sources/reservoirs instead of smuggling all four claims
//! into a lossless topology.
//!
//! [`StreamPort`] is the distinct multi-conserved-flux bundle for transported
//! mass, constituent amount, momentum, energy, and entropy. Its single
//! [`StreamEnergyChart`] is structural, and alternate thermodynamic charts
//! require context-bound exact pressure-work or Euler/Legendre crosswalks.
//! [`AccountingWindowSpec`] and [`WindowAuditReport`] separately audit
//! pre-integrated, window-bound storage and port contributions without
//! pretending that a logical clock tick is a physical duration.
//!
//! For the hard, strongly-coupled cases, [`AitkenRelaxation`] gives dynamic
//! interface relaxation: on the classic ADDED-MASS-INSTABILITY fixture (a light
//! structure in a dense fluid) naive staggering diverges, while Aitken-relaxed
//! coupling converges — demonstrated by [`iterate_fixed_relaxation`] vs
//! [`iterate_aitken`]. Deterministic; depends only on the neutral `fs-iface`
//! vocabulary and `fs-qty`'s six-base dimension vector.

pub mod acoustic_realize;
pub mod air_path;
pub mod bakeoff;
pub mod bernoulli_aperture;
pub mod bowed_string;
pub mod brass_loop;
pub mod broadband_radiation;
pub mod cavity_phs;
/// Shared Helmholtz storage and neck flow in the reciprocal modal network.
pub mod cavity_modal;
pub mod cinematic_stem;
/// Finite electrical storage and per-resistor heating on the shared thermal network.
pub mod circuit_enthalpy;
/// Original hereditary material history coupled to inertia and finite self-heating.
pub mod thermovisco_dynamics;
pub mod description;
/// Shared transactional inlet load on the original material DG propagator.
pub mod dg_port;
pub mod driving_point;
/// Material-derived axial and transverse continuum of a tensioned fiber ribbon.
pub mod fiber_ribbon;
/// Finite source gas on the shared uniform-body thermal callback.
pub mod finite_gas_thermal;
pub mod flue_loop;
pub mod gesture_path;
pub mod glottis;
pub mod iqn_ils;
pub mod material_assembly;
/// Material-resolved branch heat and finite charge/flux storage on Kirchhoff circuits.
pub mod material_circuit;
/// Pack-backed material records and fixture labels for file-driven render inputs.
pub mod material_source;
pub mod mm_line;
/// Finite material moisture transport and its native mechanical consumers.
pub mod moisture;
/// Time-domain modal realization. The exact zero-order-hold core is owned by `fs_modal::time` (IRP M1.1-W1a, owner decision OD10) and
/// every former name of this path still resolves. The opt-in broadband observer stays here: it needs `fs-vfit`, `fs-phs` and the
/// broadband radiation owner, none of which `fs-modal` may depend on.
pub mod modal_acoustic_time {
    pub use fs_modal::time::*;
    /// Opt-in broadband pressure observation using the original exact modal step.
    pub mod broadband;
}
pub mod pcm_wav;
pub mod piano_vertical;
pub mod pickup;
/// Source-bound reciprocal longitudinal piezoelectric sensors and actuators.
pub mod piezoelectric;
/// Material-driven biased magnetostrictive coil/core motion on the shared stepper.
pub mod magnetostrictive;
/// Full source-derived viscous tensors on the shared DKT and modal runtime.
pub mod plate_viscosity;
/// Accepted-step port schemas and physical accounting adapters.
pub mod port_adoption;
pub mod reed_bore;
pub mod render;
pub mod selector;
/// Independent characterized sheet bending and inertia on the original DKT owner.
pub mod sheet_plate;
pub mod speaker;
/// Shared enthalpy and material fields driving native three-dimensional solid equilibrium.
pub mod solid_enthalpy;
pub mod stribeck_friction;
pub mod string_specimen;
pub mod strip_enthalpy;
/// Full tensor material memory coupled to structural inertia and force ports.
pub mod tensor_material;
pub mod thin_plate;
/// Material-bound curved shear shells with shared mechanical and pressure ports.
pub mod shear_shell;
pub mod tract;
/// Source-bound dry sliding work and finite-body temperature feedback.
pub mod tribothermal;
/// Data-driven glottal source plus tract area function, observed at a baffled microphone.
pub mod voice_performance;
// TOMBSTONE-3ez8g.13.5: dead module, deny new consumers (see its header).
pub mod traveling_wave_line;
pub mod unilateral_contact;
pub mod vibroacoustic;
pub mod wind_line;

/// Canonical domain-neutral ports, relations, and accounting.
pub use fs_ports::*;
