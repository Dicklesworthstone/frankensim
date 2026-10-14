use super::*;
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::schedule::force::coupled::contact::friction::FrictionDrive;
use crate::string_specimen::{KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, with_uniform_circular_material_state};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, DynViscosity, Pressure, QuantitySpec};
use fs_scenario::{PrestressedString, gesture::{GestureTarget, GestureTrack, GestureValue}};

// Synthetic coefficients and body transfers test wiring, not real wood/metal data.
fn specimen(rho: f64, young: f64) -> ResolvedStringSpecimen {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, rho), ("young_modulus", Pressure::DIMS, young),
        (KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS, 10.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 0.0, 100_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance {
                source: "synthetic reciprocal string".into(), license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_with_quantity(key, QuantitySpec::dimensional(dims), ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic wire".into(), phase: "solid".into(), process: "synthetic".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let material = resolve_material_state_point(&card,
        &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap();
    with_uniform_circular_material_state(PrestressedString {
        length_m: 0.5, tension_n: 20.0, lin_density_kg_m: 1.0, axial_stiffness_n: 1.0,
        width_m: 1.0, n_modes: 4, damping_ratio: 0.0, rayleigh: None,
        bending_stiffness_n_m2: 1.0, kelvin_voigt_bending: None, relaxation_bending: None,
        polarization_detune: 0.0, moving_end: false,
    }, 0.0003, &material).unwrap().with_kelvin_voigt_bending_loss().unwrap()
}

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, fs_exec::StreamKey {
            seed: 5, kernel_id: 5, tile: 0, iteration: 0,
        }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        f(&cx, &gate)
    })
}
fn ambient() -> AmbientGas {
    AmbientGas { temperature_k: 293.15, pressure_pa: 101_325.0, relative_humidity: 0.0 }
}
fn options(mass: f64) -> ReciprocalStringOptions {
    ReciprocalStringOptions {
        sample_rate_hz: 48_000, bow_station_fraction: 0.2,
        connections: vec![InteriorStringConnection {
            station_fraction: 0.8,
            body: ModalAttachment { component: 0, shapes: vec![1.0 / mass.sqrt()] },
            stiffness_n_m: 3000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0,
        }],
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
            force_absolute_tolerance_n: 1.0e-10, force_relative_tolerance: 1.0e-9,
            sticking_tolerance_m: 1.0e-12,
        },
    }
}
fn body(mass: f64, initial_velocity: f64) -> ModalAcousticTimeModel {
    let mut body = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: (40_000.0 / mass).sqrt(), damping_ratio: 0.01,
        pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    body.restore_states(&[ModalAcousticState {
        displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: initial_velocity,
    }]).unwrap();
    body
}
fn law() -> FrictionLaw { FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 } }
fn build(cx: &Cx<'_>, gate: &CancelGate, mass: f64, initial_velocity: f64) -> ReciprocalMaterialString {
    ReciprocalMaterialString::from_specimen(cx, specimen(7800.0, 200e9), ambient(),
        vec![body(mass, initial_velocity)], law(), options(mass), gate).unwrap()
}

#[test]
fn receiver_motion_drives_the_string_back_with_the_bow_released() {
    with_cx(|cx, gate| {
        let mut model = build(cx, gate, 0.2, 0.01);
        let initial = model.system().total_energy_j().unwrap();
        assert_eq!(model.system().components()[1].states()[0].velocity_m_sqrt_kg_per_s, 0.01);
        let external = vec![0.0; model.system().mode_count()];
        for _ in 0..128 {
            let frame = model.system_mut().step(&external, FrictionDrive {
                speed_m_s: 0.0, normal_force_n: 0.0,
            }).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert_eq!(frame.drive_work_j, 0.0);
        }
        assert!(model.system().components()[0].states().iter()
            .any(|s| s.velocity_m_sqrt_kg_per_s.abs() > 1.0e-8));
        assert!(model.system().total_energy_j().unwrap() <= initial + 1.0e-12);
    });
}

#[test]
fn receiver_mass_changes_bowed_string_motion_not_only_its_pressure_filter() {
    with_cx(|cx, gate| {
        let mut a = build(cx, gate, 0.05, 0.0);
        let mut b = build(cx, gate, 0.4, 0.0);
        let zero = vec![0.0; a.system().mode_count()];
        let mut difference = 0.0;
        for _ in 0..512 {
            let drive = FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.1 };
            for model in [&mut a, &mut b] {
                let frame = model.system_mut().step(&zero, drive).unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
                assert!(frame.friction_dissipation_j >= 0.0);
            }
            difference += (a.system().components()[0].states()[0].displacement_m_sqrt_kg
                - b.system().components()[0].states()[0].displacement_m_sqrt_kg).abs();
        }
        assert!(difference > 1.0e-10, "body did not react back: {difference}");
    });
}

#[test]
fn material_swap_rebuilds_modes_losses_and_physical_port_shapes() {
    with_cx(|cx, gate| {
        let a = build(cx, gate, 0.2, 0.0);
        let replacement = specimen(15_600.0, 100e9);
        let id = replacement.material().identity();
        let b = ReciprocalMaterialString::from_specimen(cx, replacement, ambient(),
            vec![body(0.2, 0.0)], law(), options(0.2), gate).unwrap();
        assert_eq!(b.source().specimen().material().identity(), id);
        let ma = a.system().components()[0].modes();
        let mb = b.system().components()[0].modes();
        assert!(mb[0].angular_frequency_rad_s < ma[0].angular_frequency_rad_s);
        assert_ne!(ma[0].damping_ratio, mb[0].damping_ratio);
        for (&pa, &pb) in a.system().port().left.shapes.iter().zip(&b.system().port().left.shapes) {
            assert!((pb / pa - 0.5_f64.sqrt()).abs() < 1.0e-12);
        }
        assert!(mb.iter().all(|m| m.pressure_per_modal_velocity == fs_math::c64::C64::ZERO));
    });
}

#[test]
fn endpoint_zero_and_unknown_body_connections_refuse_before_assembly() {
    with_cx(|cx, gate| {
        for station in [0.0, 1.0, f64::NAN] {
            let mut spec = options(0.2);
            spec.connections[0].station_fraction = station;
            assert!(ReciprocalMaterialString::from_specimen(cx, specimen(7800.0, 200e9), ambient(),
                vec![body(0.2, 0.0)], law(), spec, gate).is_err());
        }
        for which in 0..3 {
            let mut spec = options(0.2);
            match which {
                0 => spec.connections[0].body.shapes[0] = 0.0,
                1 => spec.connections[0].body.component = 1,
                _ => spec.connections[0].stiffness_n_m = 0.0,
            }
            if which == 2 { spec.connections[0].damping_n_s_m = 0.0; }
            assert!(ReciprocalMaterialString::from_specimen(cx, specimen(7800.0, 200e9), ambient(),
                vec![body(0.2, 0.0)], law(), spec, gate).is_err());
        }
    });
}

#[test]
fn existing_scheduler_retains_source_and_checks_the_actual_bow_station() {
    with_cx(|cx, gate| {
        let schedule = GestureSchedule::try_new(700, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.01, normal_force_n: 0.1, station: 0.2 },
            events: vec![],
        }]).unwrap();
        let config = FrictionGestureConfig {
            sample_rate_hz: 48_000, samples: 129, max_block: 128, max_work: 1000,
            max_events: 16, station_fraction: 0.2,
        };
        let mut direct = build(cx, gate, 0.2, 0.0);
        let identity = direct.source().specimen().material().identity();
        let zero = vec![0.0; direct.system().mode_count()];
        let expected: Vec<_> = (0..129).map(|_| direct.system_mut().step(&zero,
            FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.1 }).unwrap().observer_pressure_pa.to_bits()).collect();
        for partition in [1, 37, 128] {
            let mut performance = build(cx, gate, 0.2, 0.0).into_scheduled(&schedule, "bow", config).unwrap();
            assert_eq!(performance.source().specimen().material().identity(), identity);
            let mut output = vec![0.0; 129];
            for block in output.chunks_mut(partition) { performance.renderer_mut().block(block).unwrap(); }
            assert_eq!(output.iter().map(|p| p.to_bits()).collect::<Vec<_>>(), expected);
            assert_eq!(performance.renderer().samples_rendered(), 129);
        }
        let wrong = FrictionGestureConfig {
            sample_rate_hz: 48_000, samples: 129, max_block: 128, max_work: 1000,
            max_events: 16, station_fraction: 0.3,
        };
        assert!(build(cx, gate, 0.2, 0.0).into_scheduled(&schedule, "bow", wrong).is_err());
    });
}
