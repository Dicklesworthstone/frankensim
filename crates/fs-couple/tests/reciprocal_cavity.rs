//! Source -> reciprocal string/body/cavity -> scheduled pressure integration.
//! Synthetic solids/gases are wiring fixtures, not measured named instruments.
use fs_couple::bowed_string::runtime::reciprocal::{
    ReciprocalMaterialString, ReciprocalStringOptions,
    plate::{MaterialBodyCavity, MaterialPlateBody},
};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig, ModalAttachment};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive, ModalFrictionConfig};
use fs_couple::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use fs_couple::string_specimen::{ResolvedStringSpecimen, with_uniform_circular_material_state,
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA};
use fs_couple::thin_plate::{PlateMaterialModel, PlateThicknessConstraint, with_uniform_plate_material_state};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_scenario::{AmbientGas, HelmholtzCavity, NeckEnd, PrestressedString, ThinPlate,
    gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue}};
use fs_tribo::FrictionLaw;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 0, kernel_id: 83, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate))
}
fn material() -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new(); let mut requirements = Vec::new();
    for (name, value, dims) in [("density", 7800.0, Density::DIMS),
        ("young_modulus", 200e9, Pressure::DIMS), ("poisson_ratio", 0.3, Dims::NONE),
        (ETA, 10.0, DynViscosity::DIMS)] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0).with("omega", 0.0, 100_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "synthetic reciprocal-cavity regression".into(),
                license: "CC0-1.0".into(), artifact: None }, observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic".into(), phase: "solid".into(), process: "test".into(), revision: 0,
    }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn specimen() -> ResolvedStringSpecimen {
    with_uniform_circular_material_state(PrestressedString {
        length_m: 0.5, tension_n: 20.0, lin_density_kg_m: 1.0, axial_stiffness_n: 0.0,
        width_m: 0.01, n_modes: 4, damping_ratio: 0.0, rayleigh: None,
        bending_stiffness_n_m2: 0.0, kelvin_voigt_bending: None, relaxation_bending: None,
        polarization_detune: 0.0, moving_end: false, tension_loss_factor: None, bending_loss_factor: None,
    }, 0.0003, &material()).unwrap().with_kelvin_voigt_bending_loss().unwrap()
}
fn body() -> MaterialPlateBody {
    let specimen = with_uniform_plate_material_state(ThinPlate {
        length_m: 0.4, width_m: 0.3, thickness_m: 1.0, density_kg_m3: 1.0,
        e1_pa: 1.0, e2_pa: 1.0, nu12: 0.0, g12_pa: 1.0, material_angle_rad: 0.0,
        damping_ratio: 0.0, thermoelastic: None, kelvin_voigt_bending: None,
        n_modes: 2, geometric_nonlinearity: false, pretension_n_m: 0.0, clamped: false,
    }, &material(), PlateMaterialModel::Isotropic, PlateThicknessConstraint::FixedThickness(0.001))
        .unwrap().with_kelvin_voigt_bending_loss().unwrap();
    MaterialPlateBody { specimen, bridge_station_fraction: 0.8, stiffness_n_m: 3000.0,
        damping_n_s_m: 0.1, rest_extension_m: 0.0,
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0 }
}
fn gas() -> GasState { GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap() }
fn cavity(gas: GasState, volume: f64) -> MaterialBodyCavity {
    MaterialBodyCavity {
        cavity: HelmholtzModalCavity::new(HelmholtzCavity {
            volume_m3: volume, neck_radius_m: 0.01, neck_length_m: 0.02,
            inner_end: NeckEnd::Flanged, outer_end: NeckEnd::Flanged,
        }, gas, CavityNeckLoss::RadiationAtResonance).unwrap(),
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0,
    }
}
fn options() -> ReciprocalStringOptions {
    ReciprocalStringOptions {
        sample_rate_hz: 48_000, bow_station_fraction: 0.2, connections: vec![],
        string_budget: ModalAcousticTimeBudget::audible_reference(),
        network: ModalCouplingConfig {
            max_modes: 8, max_connections: 4, max_setup_terms: 1024,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1.0e6, maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1.0e-10, energy_absolute_tolerance_j: 1.0e-12,
            energy_relative_tolerance: 1.0e-8,
        },
        friction: ModalFrictionConfig {
            max_iterations: 96, maximum_force_n: 100.0, maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1.0e-10, force_relative_tolerance: 1.0e-9, sticking_tolerance_m: 1.0e-12,
        },
    }
}
fn law() -> FrictionLaw { FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 } }
fn build(cx: &Cx<'_>, gate: &CancelGate, gas: GasState, volume: f64) -> ReciprocalMaterialString {
    ReciprocalMaterialString::from_material_plate_with_cavity(cx, specimen(), gas, body(),
        cavity(gas, volume), law(), options(), gate).unwrap()
}
fn drive() -> FrictionDrive { FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.1 } }

#[test]
fn legacy_air_and_explicit_gas_keep_the_same_reciprocal_pressure_bits() {
    with_cx(|cx, gate| {
        let ambient = AmbientGas { temperature_k: 293.15, pressure_pa: 101_325.0, relative_humidity: 0.0 };
        let mut a = ReciprocalMaterialString::from_material_plate(cx, specimen(), ambient,
            body(), law(), options(), gate).unwrap();
        let mut b = ReciprocalMaterialString::from_material_plate_with_gas(cx, specimen(), gas(),
            body(), law(), options(), gate).unwrap();
        assert_eq!(a.source().ambient(), Some(&ambient));
        assert_eq!(b.source().ambient(), None);
        assert_eq!(a.source().gas(), b.source().gas());
        let zeros = vec![0.0; a.system().mode_count()];
        for _ in 0..128 {
            assert_eq!(a.system_mut().step(&zeros, drive()).unwrap(), b.system_mut().step(&zeros, drive()).unwrap());
        }
    });
}

#[test]
fn cavity_geometry_changes_string_motion_and_pressure_and_survives_bow_release() {
    with_cx(|cx, gate| {
        let mut a = build(cx, gate, gas(), 0.006);
        let mut b = build(cx, gate, gas(), 0.02);
        assert_eq!(a.system().components().len(), 3);
        assert_eq!(a.source().body_connections().len(), 2);
        let zeros = vec![0.0; a.system().mode_count()];
        let (mut motion_difference, mut pressure_difference, mut neck_motion) = (0.0_f64, 0.0_f64, 0.0_f64);
        let mut work_minus_loss = [0.0; 2];
        for index in 0..1280 {
            let input = if index < 1024 { drive() } else { FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 } };
            let fa = *a.system_mut().step(&zeros, input).unwrap();
            let fb = *b.system_mut().step(&zeros, input).unwrap();
            for (i, f) in [fa, fb].iter().enumerate() {
                assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
                work_minus_loss[i] += f.drive_work_j - f.network_dissipation_j - f.friction_dissipation_j;
            }
            pressure_difference += (fa.observer_pressure_pa - fb.observer_pressure_pa).abs();
            motion_difference += (a.system().components()[0].states()[0].displacement_m_sqrt_kg
                - b.system().components()[0].states()[0].displacement_m_sqrt_kg).abs();
            neck_motion = neck_motion.max(a.system().components()[2].states()[0].velocity_m_sqrt_kg_per_s.abs());
        }
        assert!(motion_difference > 1.0e-10);
        assert!(pressure_difference > 1.0e-9);
        assert!(neck_motion > 1.0e-10);
        for (i, model) in [&a, &b].iter().enumerate() {
            let energy = model.system().total_energy_j().unwrap();
            assert!(energy > 0.0, "bow lift must not erase stored body/cavity vibration");
            assert!((energy - work_minus_loss[i]).abs() < 1.0e-10);
        }
    });
}

#[test]
fn one_gas_substitution_reaches_drag_cavity_loading_neck_radiation_and_sound() {
    with_cx(|cx, gate| {
        let replacement = GasState::try_new(&GasSpec {
            molar_mass: 0.018, gamma: 1.3, sutherland_beta: 1.1e-6,
            sutherland_s: 80.0, conductivity: ConductivityModel::Eucken,
        }, 293.15, 101_325.0).unwrap();
        let mut a = build(cx, gate, gas(), 0.01);
        let mut b = build(cx, gate, replacement, 0.01);
        assert_eq!(b.source().gas(), &replacement);
        let ca = &a.source().material_cavity().unwrap().cavity;
        let cb = &b.source().material_cavity().unwrap().cavity;
        assert_ne!(ca.neck_mass_kg(), cb.neck_mass_kg());
        assert_ne!(ca.compression_stiffness_n_m(), cb.compression_stiffness_n_m());
        assert_ne!(ca.neck_resistance_n_s_m(), cb.neck_resistance_n_s_m());
        assert_eq!(a.system().components()[0].modes()[0].angular_frequency_rad_s,
            b.system().components()[0].modes()[0].angular_frequency_rad_s);
        assert_ne!(a.system().components()[0].modes()[0].damping_ratio,
            b.system().components()[0].modes()[0].damping_ratio);
        assert_eq!(a.system().components()[1].modes(), b.system().components()[1].modes());
        assert_ne!(a.system().acceleration_pressure_weights(), b.system().acceleration_pressure_weights());
        let zeros = vec![0.0; a.system().mode_count()];
        let mut difference = 0.0;
        for _ in 0..256 {
            difference += (a.system_mut().step(&zeros, drive()).unwrap().observer_pressure_pa
                - b.system_mut().step(&zeros, drive()).unwrap().observer_pressure_pa).abs();
        }
        assert!(difference > 1.0e-9);
    });
}

#[test]
fn scheduled_cavity_sound_retains_sources_and_matches_every_callback_partition() {
    with_cx(|cx, gate| {
        let schedule = GestureSchedule::try_new(700, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.01, normal_force_n: 0.1, station: 0.2 }, events: vec![],
        }]).unwrap();
        let cfg = FrictionGestureConfig { sample_rate_hz: 48_000, samples: 257, max_block: 128,
            max_work: 1000, max_events: 16, station_fraction: 0.2 };
        let mut direct = build(cx, gate, gas(), 0.01);
        let zeros = vec![0.0; direct.system().mode_count()];
        let expected: Vec<u64> = (0..257).map(|_| direct.system_mut().step(&zeros, drive())
            .unwrap().observer_pressure_pa.to_bits()).collect();
        for size in [1, 37, 128] {
            let mut model = build(cx, gate, gas(), 0.01).into_scheduled(&schedule, "bow", cfg).unwrap();
            assert_eq!(model.source().material_cavity().unwrap().cavity.gas(), model.source().gas());
            assert_eq!(model.source().specimen().material().identity(), material().identity());
            let mut output = vec![0.0_f64; 257];
            for chunk in output.chunks_mut(size) { model.renderer_mut().block(chunk).unwrap(); }
            assert_eq!(output.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), expected);
        }
    });
}

#[test]
fn stale_cavity_gas_missing_capacity_and_conflicting_bridges_refuse() {
    with_cx(|cx, gate| {
        let other = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 80_000.0).unwrap();
        assert!(ReciprocalMaterialString::from_material_plate_with_cavity(cx, specimen(), gas(), body(),
            cavity(other, 0.01), law(), options(), gate).is_err());
        for which in 0..3 {
            let mut options = options();
            match which {
                0 => options.network.max_modes = 6,
                1 => options.network.max_connections = 2,
                _ => options.connections.push(fs_couple::bowed_string::runtime::reciprocal::InteriorStringConnection {
                    station_fraction: 0.8, body: ModalAttachment { component: 0, shapes: vec![1.0, 0.0] },
                    stiffness_n_m: 3000.0, damping_n_s_m: 0.0, rest_extension_m: 0.0,
                }),
            }
            assert!(ReciprocalMaterialString::from_material_plate_with_cavity(cx, specimen(), gas(), body(),
                cavity(gas(), 0.01), law(), options, gate).is_err());
        }
    });
}
