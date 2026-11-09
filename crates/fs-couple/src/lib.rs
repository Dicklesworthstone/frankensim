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

/// Nuclear irradiation, point kinetics and thermal coupling with their original owners.
#[cfg(feature = "nuclear-thermal")]
pub mod fission {
    pub use fs_nuclear::{
        DecayHeatGroup, EnergyChannels, FissionBudget, FissionError, FissionModel,
        FissionState, FissionStep,
    };
    pub use fs_nuclear::{burnup, cross_sections, diffusion, endf, irradiation, kinetics, radiation};
    pub use fs_nuclear::thermal::*;
}

/// Radioactive daughter inventories and deposited heat on the original thermal owner.
#[cfg(feature = "nuclear-thermal")]
pub mod radioactive {
    pub use fs_nuclear::EnergyChannels;
    pub use fs_nuclear::radiation;
    pub use fs_nuclear::endf::decay::{
        DecayEvaluation, EndfDecayReadLimits, EndfDecaySelection,
    };
    pub use fs_nuclear::endf::decay::model::{
        DecayEnergyFractions, DecayEnergyRouting, DecayModelEntry, EvaluatedDecayAssembly,
    };
    pub use fs_nuclear::endf::decay::model::emissions::{
        BurnupDecayEmission, DecayEmissionCoverage, DecayEmissionLimits, DecayEmissionStep,
        PreparedDecayEmission, UnreportedDecayEmission,
    };
    pub use fs_nuclear::endf::decay::spectra::DecayRadiation;
    pub use fs_nuclear::radiation::decay::{
        BurnupDecayRadiation, DecayRadiationClosure, DecayRadiationStep, PreparedDecayRadiation,
    };
    pub use fs_nuclear::radiation::grouped::{
        GroupedRadiationResponse, GroupedRadiationResponseLimits,
        GroupedRadiationSourceResponse, RadiationGroupResponse,
    };
    pub use fs_nuclear::radiation::transport::{
        PhotonCollision, PhotonOrdinate, PhotonSlabCell, PhotonSlabLimits,
        PhotonSlabResponse, gauss_legendre_4, prepare_slab_response,
    };
    pub use fs_nuclear::radiation::transport::material::{PhotonMassCollision, PhotonMassMaterial};
    pub use fs_nuclear::radiation::transport::refinement::{
        PhotonSlabRefinement, RefinedPhotonSlabResponse, refine_slab_response,
    };
    pub use fs_nuclear::decay::{
        DecayBranch, DecayBudget, DecayError, DecayInventory, DecayModel,
        DecayNuclide, DecayStep, PreparedDecay,
    };
    pub use fs_nuclear::decay::thermal::*;
    pub use fs_nuclear::decay::production::{
        DecayProductionChange, DecayProductionSchedule, DecayProductionStep,
        PreparedDecayProduction,
    };
    pub use fs_nuclear::decay::activation::{
        CaptureActivationBudget, CaptureActivationError, CaptureActivationState,
        CaptureActivationStep, CaptureProduct,
    };
    pub use fs_nuclear::thermal::{
        CaptureActivationCellSource, CaptureActivationHeating,
        DecayProductionCellSource, DecayProductionHeating,
        DecayRadiationEmitter, DecayRadiationHeatSource,
        PhotonSlabGeometry, PhotonSlabLayer,
        FissionHeatingAudit, FissionThermalError, NuclearHeating, RadiationReceiver,
    };
}

pub mod acoustic_realize;
/// Original local support mechanics coupled to coherent elastic-wave reflection.
pub mod elastic_reflection;
/// Compact outlet radiation reaction through the original impedance-wall field.
pub mod compact_wall_radiation;
pub mod air_path;
/// Original admitted frequency-dependent loads on accepted shell assemblies.
pub mod assembly_harmonic;
/// Joint native shear-shell equilibrium, signed membrane prestress and modes.
pub mod body_assembly;
/// Original-coordinate tensor-Prony shell loading and coherent acoustic receivers.
pub mod tensor_shell_radiation;
pub mod bakeoff;
/// Source-owned shell resonances retaining the original complex BEM pressure equations.
pub mod bem_loaded_modes;
/// Prepared force participation and analytic gas/material pole derivatives.
pub mod bem_pole_response;
pub mod bernoulli_aperture;
pub mod bowed_string;
pub mod brass_loop;
/// Sample-indexed pressure, pre-stress and valve controls on the original voice.
pub mod brass_sequence;
pub mod broadband_radiation;
/// Bulk acoustic viscous heat on original fixed-volume internal-energy cells.
pub mod bulk_viscous_heating;
pub mod cavity_phs;
/// Shared material-bound cavitating films and opt-in bubble acoustic histories.
pub mod cavitation;
/// Shared Helmholtz storage and neck flow in the reciprocal modal network.
pub mod cavity_modal;
pub mod cinematic_stem;
/// Finite electrical storage and per-resistor heating on the shared thermal network.
pub mod circuit_enthalpy;
/// Typed material-law residuals and exact graph tangents on the shared nonlinear solver.
pub mod constitutive_graph;
/// Full opening/shear path work coupled to the actual mechanical compliance.
pub mod cohesive_contact;
/// Original hereditary material history coupled to inertia and finite self-heating.
pub mod thermovisco_dynamics;
pub mod description;
/// Shared transactional inlet load on the original material DG propagator.
pub mod dg_port;
pub mod driving_point;
/// Solved electrostatic free-volume-charge loads on the native solid mesh.
pub mod electrostatic;
/// Material-derived axial and transverse continuum of a tensioned fiber ribbon.
pub mod fiber_ribbon;
/// Finite source gas on the shared uniform-body thermal callback.
pub mod finite_gas_thermal;
/// Closed-channel superfluid dissipation on original finite thermal cells.
#[cfg(feature = "superfluid")]
pub mod superfluid;
pub mod finite_solid;
/// Native rigid impacts on finite plastic or viscoelastic material histories.
pub mod inelastic_impact;
/// Original signed frame modes and power-conjugate physical vibration ports.
pub mod frame_vibration;
/// Retained nonlinear frame fibers with source-bound finite caloric feedback.
pub mod frame_enthalpy;
/// Finite-inertia guided walls coupled to native free-surface liquid momentum.
pub mod free_surface_wall;
pub mod flue_loop;
pub mod gesture_path;
mod gas_field_input;
pub mod glottis;
pub mod iqn_ils;
/// Characteristic lines whose taps move between frozen shapes in one render, with the outgoing history carried.
pub mod keyframe_line;
/// Source-bound incandescent heating on the original finite enthalpy owner.
pub mod incandescent;
/// Shared native impulse and work stage for finite-inertia induction rotors.
pub mod induction_rotor;
/// Reciprocal three-dimensional eddy-current and finite rotor-inertia feedback.
pub mod spatial_induction;
/// Source-bound mean and twice-frequency induction forcing of native linear solids.
pub mod spatial_induction_vibration;
/// Implicit guided-body motion driven by original material Stokes reactions.
#[cfg(any(feature = "experimental-stokes-motion", feature = "laminar-enthalpy"))]
pub mod stokes_motion;
/// Original BDM laminar momentum and source rheology on native finite enthalpy cells.
#[cfg(feature = "laminar-enthalpy")]
pub mod laminar_enthalpy;
/// Declared gas chemistry, finite atmospheric heat and native moving liquid.
pub mod liquid_gas;
/// Exact reciprocal drag sources on native liquid and finite rigid momenta.
pub mod liquid_momentum;
/// Conservative matter and momentum transfer between original finite liquid owners.
pub mod liquid_transfer;
/// Reciprocal molecular streams on original finite gas inventories and bulk momenta.
pub mod molecular_exchange;
/// Molecular acoustic volume absorption on original finite gas internal energy.
pub mod molecular_heating;
/// Finite gas and rigid-wall molecular shear with original caloric feedback.
pub mod molecular_shear;
pub mod material_assembly;
/// Nonlinear membrane material loss deposited in original source-bound caloric cells.
pub mod membrane_thermal;
/// Moving rigid damper networks and their actual heat on the finite thermal timeline.
pub mod rigid_damping_heat;
/// Material-resolved branch heat and finite charge/flux storage on Kirchhoff circuits.
pub mod material_circuit;
/// One fluid source for native immersed motion, sourced additional inertia and gap films.
pub mod material_immersed;
/// Pack-backed material records and fixture labels for file-driven render inputs.
pub mod material_source;
/// Full geometry-derived BEM radiation on circular duct modal flow coordinates.
pub mod modal_mouth;
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
/// Original resolved compressible flow coupled to the original causal DG bore.
pub mod nsf_dg;
/// Resolved fluid wall heat coupled to original finite spatial solid storage.
pub mod nsf_thermal;
pub mod pcm_wav;
pub mod piano_vertical;
pub mod pickup;
/// Source-bound reciprocal longitudinal piezoelectric sensors and actuators.
pub mod piezoelectric;
/// Native field/finite-inertia feedback for invariant circular planar conductors.
pub mod planar_induction;
/// Finite-amplitude periodic regimes of complete nonlinear pHS ODEs via fs-orbit.
pub mod phs_orbit;
/// Actual heat, port-work and wave-flux resonance of complete pHS models.
pub mod phs_power_resonance;
/// Material-driven biased magnetostrictive coil/core motion on the shared stepper.
pub mod magnetostrictive;
/// Full source-derived viscous tensors on the shared DKT and modal runtime.
pub mod plate_viscosity;
/// Reciprocal point impacts between a spatial plate and a finite rigid body.
pub mod plate_impact;
/// Sustained normal forces on two original linear plate trajectories.
pub mod plate_sustained;
/// Accepted-step port schemas and physical accounting adapters.
pub mod port_adoption;
pub mod reed_bore;
pub mod render;
/// Finite propellant storage and nozzle impulse on the original variable-mass vehicle.
pub mod rocket_blowdown;
/// Typed scenario potentials and normal currents on native material fields.
pub mod scenario_electric;
/// Canonical magnetic histories on the original material field owners.
pub mod scenario_magnetic;
/// Typed SI pressure histories on the original open-boundary owners.
pub mod scenario_pressure;
/// Atomic mixed velocity-inlet and pressure-outlet control updates.
pub mod scenario_open_boundary;
/// Scenario instrument response on the original accepted finite-cell heat path.
pub mod scenario_sensor;
pub mod scenario_velocity;
pub mod scenario_wall_velocity;
pub mod selector;
/// Independent characterized sheet bending and inertia on the original DKT owner.
pub mod sheet_plate;
/// Retained plastic shell histories with reciprocal source-bound layer heat.
pub mod shell_plastic_enthalpy;
/// Accepted Euler wall impulses on the original structural pressure surfaces.
#[cfg(feature = "experimental-shock-loading")]
pub mod shock_loading;
pub mod speaker;
/// Shared enthalpy and material fields driving native three-dimensional solid equilibrium.
pub mod solid_enthalpy;
/// Reciprocal native CutFEM elasticity and the original finite caloric timeline.
pub mod cut_thermoelastic;
/// Finite sublimation, deposition and reciprocal molecular recoil on original caloric stores.
pub mod solid_gas;
/// Shared conservative gravity and its analytic derivatives on deforming solids.
pub mod solid_gravity;
/// Original admitted rational radiation on the complete hereditary solid equations.
pub mod solid_rational;
/// Original stationary operating points on the shared Newton--Krylov owner.
pub mod stationary_equilibrium;
pub mod stribeck_friction;
pub mod string_specimen;
pub mod strip_enthalpy;
/// Full tensor material memory coupled to structural inertia and force ports.
pub mod tensor_material;
pub mod thin_plate;
/// Material-bound curved shear shells with shared mechanical and pressure ports.
pub mod shear_shell;
/// Source-resolved layer temperatures driving compatible native shell deformation.
pub mod laminate_thermal;
/// Coherent finite-film absorption deposited once into original thermal cells.
#[cfg(feature = "layered-wave-heating")]
pub mod layered_wave_heating;
/// Original circular-wave absorption evolving source-bound electron/heavy gas energy.
#[cfg(feature = "ionized-wave-heating")]
pub mod ionized_wave_heating;
/// Analytic porous acoustic absorption deposited once into existing thermal cells.
pub mod porous_heating;
/// Finite dry-skeleton recoil under native periodic Darcy--Forchheimer resistance.
pub mod porous_momentum;
pub mod tract;
/// Seeded turbulence noise: shared generator and filter, and the voice's opt-in glottal aspiration source.
pub mod turbulence_noise;
/// Source-bound dry sliding work and finite-body temperature feedback.
pub mod tribothermal;
/// Physical friction heat through a common interface in the existing enthalpy network.
pub mod friction_heat;
/// Data-driven glottal source plus tract area function, observed at a baffled microphone.
pub mod voice_performance;
// TOMBSTONE-3ez8g.13.5: dead module, deny new consumers (see its header).
pub mod traveling_wave_line;
pub mod unilateral_contact;
/// The airway upstream of a valve (the trachea and bronchial tree below the glottis): reflectance and source FIRs.
pub mod upstream_airway;
pub mod vibroacoustic;
pub mod wind_line;

/// Canonical domain-neutral ports, relations, and accounting.
pub use fs_ports::*;

/// Driven finite gas chambers with prescribed volume and explicit pressure work.
pub mod prescribed_gas;

/// Source-bound large-rotation rod dynamics on the existing mechanical/time owners.
#[cfg(feature = "rod-dynamics")]
pub mod rod_motion;
