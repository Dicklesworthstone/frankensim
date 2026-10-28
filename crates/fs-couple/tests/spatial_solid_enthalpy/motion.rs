//! G1/G3/G4: nodal heat and continuing pressure motion use the same owners.
use super::*;
use fs_couple::solid_enthalpy::spatial::motion::{
    SpatialHeatedValveMotion, SpatialValveMotionError,
};
use fs_couple::solid_enthalpy::valve::HeatedValveReduction;
use fs_math::c64::C64;
use fs_modal::time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_modal::transfer::ModalTransferControls;
use fs_solid::reduce::OrificeSpec;
use fs_solid::reduce::thermal::ThermalPressureLoadedValve;
use std::convert::Infallible;

const RATE: u32 = 100;
const DT: f64 = 1.0 / RATE as f64;
const DEAD_FORCE: f64 = -0.1;

fn observation(valve: &ThermalPressureLoadedValve) -> Result<Vec<C64>, Infallible> {
    // Diagnostic 2 Pa s/m times actual free-node velocity, in the new basis.
    // This supplies no acoustic radiation reaction or material damping.
    Ok(valve
        .valve()
        .source_modes()
        .iter()
        .map(|mode| C64::new(2.0 * mode.phi[0], 0.0))
        .collect())
}

fn controls() -> ModalTransferControls {
    ModalTransferControls {
        maximum_dofs: 12,
        maximum_modes: 12,
        maximum_factor_bytes: 1 << 20,
        maximum_mass_orthogonality_error: 1e-9,
        maximum_relative_eigen_residual: 1e-9,
        maximum_displacement_residual_m_sqrt_kg: 1e-9,
        maximum_kinetic_residual_j: 1e-12,
        maximum_elastic_residual_j: 1e-12,
        maximum_abs_parameter_work_j: 1.0,
        energy_tolerance_j: 1e-9,
    }
}

fn with_run<R>(
    cp: f64,
    f: impl FnOnce(&mut SpatialHeatedValveMotion<'_, '_, '_, '_, '_>, &ThermalBoundary, &Cx<'_>) -> R,
) -> R {
    with_cx(|cx| {
        let (card, caloric) = fixture(2.0, 1000.0, 1e-5, cp, 500.0, 750.0);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6.0, cx);
        let context = point(300.0);
        let source = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let fixed: Vec<_> = (0..11).collect();
        let pressure = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 1.]];
        let reduction = HeatedValveReduction {
            fixed_dofs: &fixed,
            orifice: OrificeSpec {
                face_nodes: &[1, 2, 3],
                opening_axis: [0., 0., 1.],
                width_axis: [1., 0., 0.],
                opposing_plane_offset_m: -0.5,
            },
            equilibrium_nodal_forces_n: &[(11, DEAD_FORCE)],
            nodal_load_per_pa_m2: &pressure,
            window_hz: (0.01, 20.0),
            retained_compliance_floor: 0.999,
            loss_factor: 0.0,
            source_id: "manufactured spatial-heat tetrahedron",
            static_solve: TetStaticSolveConfig::default(),
            maximum_nodal_displacement_m: 0.1,
        };
        let mut run = SpatialHeatedValveMotion::try_new(
            &source,
            vec![0.0; 4],
            reduction,
            RATE,
            &[ModalAcousticState {
                displacement_m_sqrt_kg: 3e-4,
                velocity_m_sqrt_kg_per_s: 4e-3,
            }],
            observation,
            ModalAcousticTimeBudget::audible_reference(),
            cx,
        )
        .unwrap();
        f(&mut run, &boundary(&mesh), cx)
    })
}

// Independent one-free-DOF affine tetrahedron integrals. V=1/6,
// M=rho V/10=0.1; C_zz=1.2E and thermal stress=-2E epsilon,
// K_e=E/5; u_eq=5 epsilon/3 + 5 F/E; K_g=V sigma_zz=F.
fn pencil(temperature: f64) -> (f64, f64) {
    let rise = temperature - 300.0;
    let e = 1000.0 * (1.0 - 0.0025 * rise);
    let strain = 1e-5 * (rise + 0.005 * rise * rise);
    (
        5.0 * strain / 3.0 + 5.0 * DEAD_FORCE / e,
        e / 5.0 + DEAD_FORCE,
    )
}

#[test]
fn g1_nodal_heat_changes_native_motion_and_observation_without_resetting_state() {
    with_run(10.0, |run, boundary, cx| {
        let input = ScalarField::uniform("uniform finite heating", 1200.0).unwrap();
        let initial = run.motion().nodal_motion(cx).unwrap();
        let mut u = initial.displacement_m[3][2];
        let mut v = initial.velocity_m_s[3][2];
        let mut temperature = 300.0;
        for step in 1..=8 {
            let pressure = if step <= 4 { 0.04 } else { 0.0 };
            let (equilibrium, k) = pencil(temperature);
            let omega = (k / 0.1).sqrt();
            let (sin, cos) = (omega * DT).sin_cos();
            let driven_equilibrium = equilibrium + pressure / k;
            let next_u = driven_equilibrium + (u - driven_equilibrium) * cos + v * sin / omega;
            let next_v = -(u - driven_equilibrium) * omega * sin + v * cos;
            let result = run
                .step(
                    pressure,
                    boundary,
                    None,
                    &input,
                    thermal_config(),
                    controls(),
                    observation,
                    cx,
                )
                .unwrap();
            // rho cp dT/dt = volumetric generation, with sealed faces.
            temperature += 1200.0 * DT / (6.0 * 10.0);
            let (next_equilibrium, next_k) = pencil(temperature);
            let old_energy = 0.05 * next_v.powi(2) + 0.5 * k * (next_u - equilibrium).powi(2);
            let new_energy =
                0.05 * next_v.powi(2) + 0.5 * next_k * (next_u - next_equilibrium).powi(2);
            let nodes = run.motion().nodal_motion(cx).unwrap();
            near(nodes.displacement_m[3][2], next_u, 2e-12);
            near(nodes.velocity_m_s[3][2], next_v, 2e-11);
            near(result.endpoint_observer_pressure_pa, 2.0 * next_v, 4e-11);
            near(
                result.pressure_step.modal.input_work_j,
                pressure * (next_u - u),
                1e-12,
            );
            near(
                result.pressure_step.swept_volume_change_m3,
                next_u - u,
                2e-12,
            );
            near(
                result.transfer.parameter_work_j,
                new_energy - old_energy,
                1e-11,
            );
            near(result.transfer.retained_modal_energy_j, new_energy, 1e-11);
            near(
                result.thermal.stored_energy_change_j,
                1200.0 * DT / 6.0,
                1e-8,
            );
            assert!(result.thermal.energy_residual_j.abs() < 1e-8);
            for &t in &result.thermal.temperature {
                near(t, temperature, 1e-8);
            }
            assert_eq!(
                run.specific_enthalpy_j_kg(),
                result.thermal.specific_enthalpy_j_kg
            );
            assert_eq!(run.motion().accepted_steps(), step);
            assert_eq!(result.time_s, step as f64 * DT);
            let state = run.state(cx).unwrap();
            near(state.samples()[0].phase.temperature_k(), temperature, 1e-8);
            near(
                run.motion().valve().equilibrium().pencil().total_mass_kg,
                1.0,
                1e-12,
            );
            (u, v) = (next_u, next_v);
        }
    });
}

#[test]
fn g3_zero_heat_retains_original_trajectory_and_heat_capacity_changes_real_motion() {
    with_run(10.0, |run, boundary, cx| {
        let input = ScalarField::uniform("no heating", 0.0).unwrap();
        let mut unchanged = run.motion().clone();
        for _ in 0..8 {
            let expected = unchanged.step_pressure(0.04, cx).unwrap();
            let result = run
                .step(
                    0.04,
                    boundary,
                    None,
                    &input,
                    thermal_config(),
                    controls(),
                    observation,
                    cx,
                )
                .unwrap();
            for (a, b) in run
                .motion()
                .model()
                .states()
                .iter()
                .zip(unchanged.model().states())
            {
                near(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg, 1e-13);
                near(
                    a.velocity_m_sqrt_kg_per_s,
                    b.velocity_m_sqrt_kg_per_s,
                    1e-12,
                );
            }
            near(
                result.endpoint_observer_pressure_pa,
                expected.modal.observer_pressure_pa,
                1e-11,
            );
            near(result.transfer.parameter_work_j, 0.0, 1e-14);
            assert_eq!(run.specific_enthalpy_j_kg(), &[0.0; 4]);
        }
    });
    let endpoint = |cp| {
        with_run(cp, |run, boundary, cx| {
            let input = ScalarField::uniform("same heating", 1200.0).unwrap();
            let mut pressure = 0.0;
            for _ in 0..8 {
                pressure = run
                    .step(
                        0.04,
                        boundary,
                        None,
                        &input,
                        thermal_config(),
                        controls(),
                        observation,
                        cx,
                    )
                    .unwrap()
                    .endpoint_observer_pressure_pa;
            }
            (
                run.motion().nodal_motion(cx).unwrap(),
                pressure,
                run.state(cx).unwrap().samples()[0].phase.temperature_k(),
            )
        })
    };
    let quick = endpoint(10.0);
    let slow = endpoint(20.0);
    assert!(quick.2 > slow.2);
    assert!((quick.0.displacement_m[3][2] - slow.0.displacement_m[3][2]).abs() > 1e-8);
    assert!((quick.1 - slow.1).abs() > 1e-8);
}

#[test]
fn g4_late_refusals_keep_heat_motion_observation_and_retry_identical() {
    with_run(10.0, |run, boundary, cx| {
        let input = ScalarField::uniform("finite heating", 1200.0).unwrap();
        let before_h = run.specific_enthalpy_j_kg().to_vec();
        let before = run.motion().clone();
        let assert_unchanged = |run: &SpatialHeatedValveMotion<'_, '_, '_, '_, '_>| {
            assert_eq!(run.specific_enthalpy_j_kg(), before_h);
            assert_eq!(run.motion().model().states(), before.model().states());
            assert_eq!(run.motion().valve().identity(), before.valve().identity());
            assert_eq!(run.motion().accepted_steps(), before.accepted_steps());
        };
        assert!(matches!(
            run.step(
                0.04,
                boundary,
                None,
                &input,
                thermal_config(),
                controls(),
                |_| Err::<Vec<C64>, _>("observation refused"),
                cx
            ),
            Err(SpatialValveMotionError::Observation("observation refused"))
        ));
        assert_unchanged(run);
        let mut tight = controls();
        tight.maximum_abs_parameter_work_j = 0.0;
        assert!(
            run.step(
                0.04,
                boundary,
                None,
                &input,
                thermal_config(),
                tight,
                observation,
                cx
            )
            .is_err()
        );
        assert_unchanged(run);
        let outside = ScalarField::uniform("source leaves supported state", 1e9).unwrap();
        assert!(
            run.step(
                0.04,
                boundary,
                None,
                &outside,
                thermal_config(),
                controls(),
                observation,
                cx
            )
            .is_err()
        );
        assert_unchanged(run);
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let failed = run.step(
                0.04,
                boundary,
                None,
                &input,
                thermal_config(),
                controls(),
                |valve| {
                    gate.request();
                    observation(valve)
                },
                &cancelled,
            );
            assert!(failed.is_err());
        });
        assert_unchanged(run);
        let actual = run
            .step(
                0.04,
                boundary,
                None,
                &input,
                thermal_config(),
                controls(),
                observation,
                cx,
            )
            .unwrap();
        let actual_motion = run.motion().model().states().to_vec();
        with_run(10.0, |fresh, boundary, cx| {
            let expected = fresh
                .step(
                    0.04,
                    boundary,
                    None,
                    &input,
                    thermal_config(),
                    controls(),
                    observation,
                    cx,
                )
                .unwrap();
            assert_eq!(actual_motion, fresh.motion().model().states());
            assert_eq!(
                actual.thermal.specific_enthalpy_j_kg,
                expected.thermal.specific_enthalpy_j_kg
            );
            assert_eq!(
                actual.endpoint_observer_pressure_pa,
                expected.endpoint_observer_pressure_pa
            );
            assert_eq!(
                actual.transfer.parameter_work_j,
                expected.transfer.parameter_work_j
            );
        });
    });
}
