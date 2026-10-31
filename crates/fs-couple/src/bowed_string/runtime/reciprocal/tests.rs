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
        ("poisson_ratio", fs_qty::Dims::NONE, 0.3),
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
        polarization_detune: 0.0, moving_end: false, tension_loss_factor: None, bending_loss_factor: None,
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


fn material_plate(rho: f64) -> super::plate::MaterialPlateBody {
    use crate::thin_plate::{PlateMaterialModel, PlateThicknessConstraint, with_uniform_plate_material_state};
    let material = specimen(rho, 200e9);
    let plate = with_uniform_plate_material_state(fs_scenario::ThinPlate {
        length_m: 0.4, width_m: 0.3, thickness_m: 1.0, density_kg_m3: 1.0,
        e1_pa: 1.0, e2_pa: 1.0, nu12: 0.0, g12_pa: 1.0, material_angle_rad: 0.0,
        damping_ratio: 0.0, thermoelastic: None, kelvin_voigt_bending: None,
        n_modes: 2, geometric_nonlinearity: false, pretension_n_m: 0.0, clamped: false,
    }, material.material(), PlateMaterialModel::Isotropic,
        PlateThicknessConstraint::FixedThickness(0.001)).unwrap()
        .with_kelvin_voigt_bending_loss().unwrap();
    super::plate::MaterialPlateBody {
        specimen: plate, bridge_station_fraction: 0.8, stiffness_n_m: 3000.0,
        damping_n_s_m: 0.1, rest_extension_m: 0.0,
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0,
    }
}
fn plate_options() -> ReciprocalStringOptions {
    let mut config = options(0.2);
    config.connections.clear();
    config
}
fn with_plate(cx: &Cx<'_>, gate: &CancelGate, rho: f64) -> ReciprocalMaterialString {
    ReciprocalMaterialString::from_material_plate(cx, specimen(7800.0, 200e9), ambient(),
        material_plate(rho), law(), plate_options(), gate).unwrap()
}

#[test]
fn source_plate_generates_actual_body_modes_conjugate_bridge_and_acceleration_pressure() {
    with_cx(|cx, gate| {
        let body = material_plate(7800.0);
        let expected = crate::thin_plate::certified_radiators(body.specimen.plate()).unwrap();
        let material_id = body.specimen.material().identity();
        let gas = crate::acoustic_realize::gas_state(ambient()).unwrap();
        let mut model = ReciprocalMaterialString::from_material_plate(cx, specimen(7800.0, 200e9),
            ambient(), body, law(), plate_options(), gate).unwrap();
        assert_eq!(model.source().material_plate().unwrap().specimen.material().identity(), material_id);
        let weights = model.system().acceleration_pressure_weights().unwrap();
        let modes = model.system().components()[1].modes();
        let bridge = &model.source().options().connections[0].body.shapes;
        assert_eq!(modes.len(), expected.len());
        for (k, radiator) in expected.iter().enumerate() {
            assert_eq!(modes[k].angular_frequency_rad_s, radiator.omega);
            assert_eq!(modes[k].damping_ratio, radiator.zeta);
            assert!((bridge[k] - radiator.drive_participation / radiator.mass_kg.sqrt()).abs() < 1.0e-12);
            assert_eq!(weights[4 + k], radiator.radiate_mass_normalized(1.0, gas.density, 1.0));
        }
        assert!(weights[..4].iter().all(|w| *w == 0.0));
        let zero = vec![0.0; model.system().mode_count()];
        let mut peak = 0.0_f64;
        for _ in 0..256 {
            let frame = model.system_mut().step(&zero,
                FrictionDrive { speed_m_s: 0.01, normal_force_n: 0.1 }).unwrap();
            peak = peak.max(frame.observer_pressure_pa.abs());
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(peak > 1.0e-12, "material body never emitted physical pressure");
    });
}

#[test]
fn body_material_swap_changes_realized_mechanics_and_scheduled_pressure_without_code_changes() {
    with_cx(|cx, gate| {
        let a = with_plate(cx, gate, 7800.0);
        let b = with_plate(cx, gate, 15_600.0);
        let wa = a.system().components()[1].modes()[0].angular_frequency_rad_s;
        let wb = b.system().components()[1].modes()[0].angular_frequency_rad_s;
        assert!((wb / wa - 0.5_f64.sqrt()).abs() < 1.0e-6);
        let source_a = a.source().material_plate().unwrap().specimen.material().identity();
        let source_b = b.source().material_plate().unwrap().specimen.material().identity();
        assert_ne!(source_a, source_b);
        let schedule = GestureSchedule::try_new(700, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.01, normal_force_n: 0.1, station: 0.2 },
            events: vec![],
        }]).unwrap();
        let config = FrictionGestureConfig {
            sample_rate_hz: 48_000, samples: 256, max_block: 256, max_work: 1000,
            max_events: 16, station_fraction: 0.2,
        };
        let mut a = a.into_scheduled(&schedule, "bow", config).unwrap();
        let mut b = b.into_scheduled(&schedule, "bow", config).unwrap();
        assert_eq!(a.source().material_plate().unwrap().specimen.material().identity(), source_a);
        assert_eq!(b.source().material_plate().unwrap().specimen.material().identity(), source_b);
        let mut pa = vec![0.0; 256];
        let mut pb = vec![0.0; 256];
        a.renderer_mut().block(&mut pa).unwrap();
        for block in pb.chunks_mut(37) { b.renderer_mut().block(block).unwrap(); }
        assert!(pa.iter().zip(&pb).any(|(x, y)| (x - y).abs() > 1.0e-12));
        let mut partitioned = with_plate(cx, gate, 7800.0).into_scheduled(&schedule, "bow", config).unwrap();
        let mut check = vec![0.0; 256];
        for block in check.chunks_mut(37) { partitioned.renderer_mut().block(block).unwrap(); }
        assert_eq!(pa.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            check.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    });
}

#[test]
fn material_plate_constructor_never_ignores_authored_connections_or_invalid_observers() {
    with_cx(|cx, gate| {
        assert!(ReciprocalMaterialString::from_material_plate(cx, specimen(7800.0, 200e9), ambient(),
            material_plate(7800.0), law(), options(0.2), gate).is_err());
        for distance in [0.0, -1.0, f64::NAN] {
            let mut body = material_plate(7800.0);
            body.listener_m = distance;
            assert!(ReciprocalMaterialString::from_material_plate(cx, specimen(7800.0, 200e9), ambient(),
                body, law(), plate_options(), gate).is_err());
        }
    });
}
