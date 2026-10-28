//! Original six-coordinate native bodies exchange full cohesive path work.
use fs_couple::cohesive_contact::CohesiveContactConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::cohesive::{
    CohesiveModalConfig, CohesiveModalError, CohesiveModalSystem,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_material::cohesive::{BilinearCohesive, CohesiveViscosity, mixed_mode::MixedModeParameters};
use fs_math::{c64::C64, det};

const MA: f64 = 2.0;
const MB: f64 = 3.0;
const AREA: f64 = 1e-6;
const RATE: u32 = 5_000;

fn law() -> BilinearCohesive {
    BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 4e8,
            normal_strength_pa: 2e5,
            tangential_strength_pa: 1.6e5,
            mode_i_fracture_energy_j_m2: 100.0,
            mode_ii_fracture_energy_j_m2: 400.0,
            bk_exponent: 2.0,
        },
        2e9,
        0.02,
    )
    .unwrap()
}
fn config() -> CohesiveModalConfig {
    CohesiveModalConfig {
        contact: CohesiveContactConfig {
            max_iterations: 24,
            max_material_evaluations: 400,
            max_backtracks: 12,
            jump_absolute_tolerance_m: 1e-14,
            jump_relative_tolerance: 1e-11,
            maximum_jump_increment_m: 0.02,
            maximum_force_n: 100.0,
            material_energy_tolerance_j_m2: 1e-8,
        },
        force_absolute_tolerance_n: 1e-9,
        force_relative_tolerance: 1e-9,
    }
}
fn columns() -> [Vec<f64>; 3] {
    core::array::from_fn(|i| {
        let mut column = vec![0.0; 6];
        column[i] = -1.0 / det::sqrt(MA);
        column[3 + i] = 1.0 / det::sqrt(MB);
        column
    })
}
fn network(relative: [f64; 3], energy_ceiling: f64, cross_axis_spring: bool) -> CoupledModalSystem {
    let center = [0.012, -0.006, 0.003];
    let models = [MA, MB]
        .into_iter()
        .enumerate()
        .map(|(side, mass)| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                RATE,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO
                    };
                    3
                ],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            let states: Vec<_> = (0..3)
                .map(|i| ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: det::sqrt(mass)
                        * (center[i]
                            + if side == 0 {
                                -MB / (MA + MB)
                            } else {
                                MA / (MA + MB)
                            } * relative[i]),
                })
                .collect();
            model.restore_states(&states).unwrap();
            model
        })
        .collect();
    let links = if cross_axis_spring {
        vec![ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![
                    0.7 / det::sqrt(MA),
                    0.3 / det::sqrt(MA),
                    -0.2 / det::sqrt(MA),
                ],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![
                    0.7 / det::sqrt(MB),
                    0.3 / det::sqrt(MB),
                    -0.2 / det::sqrt(MB),
                ],
            },
            stiffness_n_m: 1e8,
            damping_n_s_m: 100.0,
            rest_extension_m: 0.0,
        }]
    } else {
        vec![]
    };
    CoupledModalSystem::new(
        models,
        links,
        ModalCouplingConfig {
            max_modes: 6,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: energy_ceiling,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 1e5,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-10,
        },
        &CancelGate::new(),
    )
    .unwrap()
}
fn system(law: &BilinearCohesive, relative: [f64; 3]) -> CohesiveModalSystem<'_> {
    CohesiveModalSystem::new(
        network(relative, 1.0, false),
        law,
        AREA,
        columns(),
        config(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn states(system: &CohesiveModalSystem<'_>) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}
fn momentum(system: &CohesiveModalSystem<'_>) -> [f64; 3] {
    core::array::from_fn(|i| {
        det::sqrt(MA) * system.components()[0].states()[i].velocity_m_sqrt_kg_per_s
            + det::sqrt(MB) * system.components()[1].states()[i].velocity_m_sqrt_kg_per_s
    })
}
fn relative_velocity(system: &CohesiveModalSystem<'_>) -> [f64; 3] {
    core::array::from_fn(|i| {
        system.components()[1].states()[i].velocity_m_sqrt_kg_per_s / det::sqrt(MB)
            - system.components()[0].states()[i].velocity_m_sqrt_kg_per_s / det::sqrt(MA)
    })
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn native_two_body_motion_preserves_three_momenta_and_complete_debonding_energy() {
    for law in [
        law(),
        law().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 100.0).unwrap()),
    ] {
        let mut system = system(&law, [0.05, 0.04, -0.025]);
        let initial = system.total_energy_j().unwrap();
        let initial_momentum = momentum(&system);
        let mut losses = 0.0;
        let mut fracture = 0.0;
        let mut viscous = 0.0;
        let mut partial = false;
        for _ in 0..500 {
            let frame = *system.step(&[0.0; 6], [0.0; 3]).unwrap();
            partial |= frame.damage > 0.0 && frame.damage < 1.0;
            losses += frame.network_dissipation_j
                + frame.fracture_dissipation_j
                + frame.viscous_dissipation_j;
            fracture += frame.fracture_dissipation_j;
            viscous += frame.viscous_dissipation_j;
            assert_eq!(frame.external_work_j, 0.0);
            assert_eq!(frame.prescribed_work_j, 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert!(frame.material_evaluations <= config().contact.max_material_evaluations);
            for i in 0..3 {
                close(momentum(&system)[i], initial_momentum[i], 3e-14);
            }
            close(system.total_energy_j().unwrap() + losses, initial, 3e-11);
        }
        assert!(partial && system.point().failed());
        assert_eq!(system.samples_rendered(), 500);
        close(fracture, AREA * system.point().dissipated_j_m2(), 1e-13);
        assert_eq!(viscous > 0.0, law.viscosity().is_some());
    }
}

#[test]
fn native_cross_axis_compliance_and_prescribed_port_recover_exact_elastic_work() {
    let law = law();
    let mut system = CohesiveModalSystem::new(
        network([0.0; 3], 1.0, true),
        &law,
        AREA,
        columns(),
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let s = system.compliance_m_per_n();
    assert!(s[0][1].abs() > 1e-10 && s[0][2].abs() > 1e-10);
    let target = [8e-5, 6e-5, -4e-5];
    let stiffness = [1e9, 4e8, 4e8];
    let expected_force: [f64; 3] = core::array::from_fn(|i| 0.5 * AREA * stiffness[i] * target[i]);
    let drive = core::array::from_fn(|i| {
        (target[i] + (0..3).map(|j| s[i][j] * expected_force[j]).sum::<f64>())
            / system.sample_period_s()
    });
    let before = system.total_energy_j().unwrap();
    let frame = *system.step(&[0.0; 6], drive).unwrap();
    for i in 0..3 {
        close(frame.jump_m[i], target[i], 3e-13);
        close(frame.mean_force_n[i], expected_force[i], 3e-9);
    }
    close(
        system.total_energy_j().unwrap() - before + frame.network_dissipation_j,
        frame.prescribed_work_j,
        2e-12,
    );
    assert_eq!(frame.fracture_dissipation_j, 0.0);
}

#[test]
fn native_driven_mixed_mode_rupture_and_reclosure_keep_fracture_separate_from_viscosity() {
    let law = law().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 100.0).unwrap());
    let mut system = system(&law, [0.0; 3]);
    let mut previous_damage = 0.0;
    let mut history = None;
    for target in [
        [4e-4, 3e-4, -1e-4],
        [5e-5, -4e-5, 2e-5],
        [1e-3, 0.012, -0.003],
        [-5e-4, 0.006, 0.0],
    ] {
        let old = system.point();
        let dt = system.sample_period_s();
        let (_, mean) = old
            .trial_over_time_with_mean_traction(target, dt, 1e-8)
            .unwrap();
        let velocity = relative_velocity(&system);
        let s = system.compliance_m_per_n();
        let drive = core::array::from_fn(|i| {
            (target[i] - old.jump_m()[i] - dt * velocity[i]
                + (0..3).map(|j| s[i][j] * AREA * mean[j]).sum::<f64>())
                / dt
        });
        let frame = *system.step(&[0.0; 6], drive).unwrap();
        assert!(frame.damage >= previous_damage);
        previous_damage = frame.damage;
        if history.is_none() {
            history = system.point().mixed_mode_history();
        }
        assert_eq!(system.point().mixed_mode_history(), history);
        for i in 0..3 {
            close(frame.jump_m[i], target[i], 5e-12);
        }
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
    let frame = system.last_frame().unwrap();
    assert!(system.point().failed());
    assert_eq!(frame.fracture_dissipation_j, 0.0);
    assert_eq!(frame.viscous_dissipation_j, 0.0);
    assert_eq!([frame.mean_force_n[1], frame.mean_force_n[2]], [0.0; 2]);
    assert!(frame.mean_force_n[0] < 0.0 && frame.endpoint_traction_pa[0] < 0.0);
}

#[test]
fn late_total_energy_and_cancellation_refusal_preserve_native_motion_damage_and_replay() {
    let law = law();
    // Center-of-mass motion is retained; the ceiling leaves enough for it but
    // rejects the later cohesive elastic energy after native mechanics stages.
    let ceiling =
        0.5 * (MA + MB) * (0.012_f64.powi(2) + 0.006_f64.powi(2) + 0.003_f64.powi(2)) + 1e-7;
    let build = || {
        CohesiveModalSystem::new(
            network([0.0; 3], ceiling, false),
            &law,
            AREA,
            columns(),
            config(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let mut system = build();
    let mut replay = build();
    for target in [&mut system, &mut replay] {
        target.step(&[0.0; 6], [0.02, 0.015, -0.005]).unwrap();
    }
    let before = states(&system);
    let point = system.point();
    let last = system.last_frame().copied();
    let gate = CancelGate::new();
    gate.request();
    assert!(system.step_under_gate(&[0.0; 6], [0.01; 3], &gate).is_err());
    let refused = system.step(&[0.0; 6], [1.0, 0.0, 0.0]);
    assert!(matches!(
        refused,
        Err(CohesiveModalError::Network(ModalCouplingError::Budget {
            what: "total energy including cohesive interface",
            ..
        }))
    ));
    assert_eq!(states(&system), before);
    assert_eq!(system.point(), point);
    assert_eq!(system.last_frame().copied(), last);
    let drive = [0.01, -0.01, 0.002];
    assert_eq!(
        *system.step(&[0.0; 6], drive).unwrap(),
        *replay.step(&[0.0; 6], drive).unwrap()
    );
    assert_eq!(states(&system), states(&replay));
}

#[test]
fn native_source_binding_retains_original_ordered_receipts_and_rejects_reversed_materials() {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
        UncertaintyModel,
    };
    use fs_material::cohesive::source::*;
    use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
    use fs_qty::Pressure;
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e9,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured native cohesive fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured fixture".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit bond".into(),
            third_body: None,
            environment: "fixed source query".into(),
            history: "initial bond".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let source = resolve_bilinear_cohesive(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap();
    assert!(
        CohesiveModalSystem::from_source(
            network([0.0; 3], 1.0, false),
            &source,
            [&card.surface_b().material, &card.surface_a().material],
            AREA,
            columns(),
            config(),
            &CancelGate::new()
        )
        .is_err()
    );
    let mut system = CohesiveModalSystem::from_source(
        network([0.01, 0.02, 0.0], 1.0, false),
        &source,
        [&card.surface_a().material, &card.surface_b().material],
        AREA,
        columns(),
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let receipts = source.resolved().properties().as_ptr();
    system.step(&[0.0; 6], [0.0; 3]).unwrap();
    assert!(std::ptr::eq(system.source().unwrap(), &source));
    assert!(std::ptr::eq(system.law(), source.law()));
    assert_eq!(
        system.source().unwrap().resolved().properties().as_ptr(),
        receipts
    );
}
