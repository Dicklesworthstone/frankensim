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
pub mod cinematic_stem;
pub mod driving_point;
pub mod flue_loop;
pub mod glottis;
pub mod iqn_ils;
pub mod material_assembly;
pub mod mm_line;
pub mod modal_acoustic_time;
pub mod pcm_wav;
pub mod piano_vertical;
pub mod pickup;
pub mod reed_bore;
pub mod render;
pub mod selector;
pub mod speaker;
pub mod stribeck_friction;
pub mod string_specimen;
pub mod thin_plate;
pub mod tract;
// TOMBSTONE-3ez8g.13.5: dead module, deny new consumers (see its header).
pub mod traveling_wave_line;
pub mod unilateral_contact;
pub mod vibroacoustic;
pub mod wind_line;

/// Canonical domain-neutral ports, relations, and accounting.
pub use fs_ports::*;
