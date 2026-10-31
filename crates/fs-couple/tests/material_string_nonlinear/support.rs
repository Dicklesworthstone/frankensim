use fs_couple::bowed_string::runtime::reciprocal::bank::MaterialStringSource;
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::ImpactConfig;
use fs_couple::string_specimen::{KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA, with_uniform_circular_material_state};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::PrestressedString;

pub fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena, fs_exec::StreamKey { seed: 6, kernel_id: 5, tile: 0, iteration: 0 },
        fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic), &gate))
}
pub fn material(rho: f64, young: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [("density", Density::DIMS, rho), ("young_modulus", Pressure::DIMS, young),
        ("poisson_ratio", Dims::NONE, 0.3), (ETA, DynViscosity::DIMS, 10.0),
        ("G", Pressure::DIMS, 60e9), ("eta_s", DynViscosity::DIMS, 1000.0)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, 200_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source: "synthetic nonlinear string".into(), license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
pub fn template(n_modes: usize) -> PrestressedString {
    PrestressedString { length_m: 0.5, tension_n: 20.0, lin_density_kg_m: 1.0, axial_stiffness_n: 1.0, width_m: 1.0,
        n_modes, damping_ratio: 0.0, rayleigh: None, bending_stiffness_n_m2: 0.0, kelvin_voigt_bending: None,
        relaxation_bending: None, polarization_detune: 0.0, moving_end: false,
        tension_loss_factor: None, bending_loss_factor: None }
}
pub fn source(rho: f64, young: f64, modes: usize) -> MaterialStringSource {
    MaterialStringSource::Homogeneous(with_uniform_circular_material_state(template(modes), 0.0003, &material(rho, young))
        .unwrap().with_kelvin_voigt_bending_loss().unwrap())
}
pub fn gas() -> GasState { GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap() }
pub fn initial(n: usize, amplitude_m: f64) -> Vec<ModalAcousticState> {
    let mut state = vec![ModalAcousticState::default(); n];
    let mu = 7800.0 * core::f64::consts::PI * 0.0003_f64.powi(2);
    state[0].displacement_m_sqrt_kg = amplitude_m * (mu * 0.5 * 0.5).sqrt();
    state
}
pub fn config() -> ImpactConfig {
    ImpactConfig { dt_s: 1.0/48_000.0, max_steps: 10_000, maximum_energy_j: 2.0,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8, maximum_generalized_force: 1e4 }
}
/// Clock for bodies with twist modes. ImpactSystem's midpoint block admits only
/// the one-cent band, and the layered torsion pole (about 26000 rad/s) needs
/// dt <= 1.5e-6 s.
#[allow(dead_code)]
pub const TORSION_RATE_HZ: u32 = 768_000;
#[allow(dead_code)]
pub fn config_at(rate_hz: u32) -> ImpactConfig {
    ImpactConfig { dt_s: 1.0/f64::from(rate_hz), ..config() }
}
