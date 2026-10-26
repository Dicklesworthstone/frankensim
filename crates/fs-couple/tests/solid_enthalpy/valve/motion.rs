//! G1/G3/G4: finite source heat changes the same continuing native trajectory.

use super::*;
use fs_couple::solid_enthalpy::valve::motion::{HeatedValveMotion, HeatedValveMotionError};
use fs_math::c64::C64;
use fs_modal::time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_modal::transfer::ModalTransferControls;
use fs_solid::reduce::thermal::ThermalPressureLoadedValve;
use std::cell::Cell;
use std::convert::Infallible;

const RATE: u32 = 20;
const DT: f64 = 1.0 / RATE as f64;

fn observations(valve: &ThermalPressureLoadedValve) -> Result<Vec<C64>, Infallible> {
    // Explicit diagnostic p = 2 Pa s/m * physical free-node velocity. This is
    // a source-basis observation, not a fitted acoustic radiation model.
    Ok(valve
        .valve()
        .source_modes()
        .iter()
        .map(|mode| C64::new(2.0 * mode.phi[0], 0.0))
        .collect())
}

fn transfer_controls() -> ModalTransferControls {
    ModalTransferControls {
        maximum_dofs: 24,
        maximum_modes: 24,
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

fn with_motion<R>(
    cp: f64,
    power_w: f64,
    origin_s: f64,
    f: impl FnOnce(&mut HeatedValveMotion<'_, '_, '_, '_>, &Cx<'_>) -> R,
) -> R {
    with_cx(|cx| {
        let fixed: Vec<_> = (0..24).filter(|&d| d != 21).collect();
        let mut pressure = [[0.0; 3]; 8];
        pressure[7][0] = 1.0;
        let mut input = reduction(&fixed, &pressure);
        input.loss_factor = 0.0;
        let context = point(300.0);
        let (card, caloric) = material(1000.0, 1e-5, cp, 0.0, false);
        let cells: Vec<_> = (0..6).map(|id| (id, power_w)).collect();
        let thermal = initial(caloric.curve(), &cells, 1.0);
        let thermal = if origin_s == 0.0 {
            thermal
        } else {
            // Birth at the requested epoch, so the empty prefix contains no
            // heat or material aging and does not alter the reference mass.
            let events = thermal
                .cells()
                .values()
                .cloned()
                .map(|cell| TimedEvent {
                    time_s: origin_s,
                    event: Event::Birth(cell),
                })
                .collect();
            DepositionSimulation::new(
                events,
                &[(MaterialId(0), caloric.curve() as &dyn ThermalLaw)],
                DepositionBudget {
                    max_cells: 8,
                    max_links: 8,
                    max_events: 16,
                    max_substeps: 1000,
                },
            )
            .unwrap()
            .advance_to(
                origin_s,
                StepConfig {
                    max_dt_s: origin_s,
                    energy_tolerance_j: 1e-8,
                },
                || Ok(()),
            )
            .unwrap()
            .0
        };
        let source = SolidEnthalpy::try_new(
            &thermal,
            &CUBE,
            &CUBE_TETS,
            (0..6)
                .map(|id| super::super::source(id, &card, caloric.curve(), &context))
                .collect(),
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let mut run = HeatedValveMotion::try_new(
            &source,
            thermal,
            input,
            RATE,
            &[ModalAcousticState {
                displacement_m_sqrt_kg: 3e-4,
                velocity_m_sqrt_kg_per_s: 4e-3,
            }],
            observations,
            ModalAcousticTimeBudget::audible_reference(),
            cx,
        )
        .unwrap();
        f(&mut run, cx)
    })
}

// Independent one-free-DOF cube integrals from the manufactured specimen:
// M=rho/10=0.6 kg, K_e=2E/3, f_thermal=E*epsilon/3 and
// K_g=E*(u_equilibrium-3*epsilon)/3. No reducer or transfer helper is reused.
fn pencil_at(temperature_k: f64) -> (f64, f64) {
    let rise = temperature_k - 300.0;
    let strain = 1e-5 * (rise + 0.005 * rise * rise);
    let e = 1000.0 * (1.0 - 0.002 * rise);
    let equilibrium = 0.5 * strain + FORCE[0].1 / (2.0 * e / 3.0);
    let stiffness = 2.0 * e / 3.0 + e * (equilibrium - 3.0 * strain) / 3.0;
    (equilibrium, stiffness)
}

#[test]
fn g1_finite_heating_and_continuing_motion_match_an_independent_piecewise_oscillator() {
    with_motion(10.0, 200.0, 0.0, |run, cx| {
        let initial = run.motion().nodal_motion(cx).unwrap();
        let mut u = initial.displacement_m[7][0];
        let mut v = initial.velocity_m_s[7][0];
        let mut temperature = 300.0;
        for step in 1..=8 {
            let pressure = if step <= 4 { 0.04 } else { 0.0 };
            let (equilibrium, k) = pencil_at(temperature);
            let omega = (k / 0.6).sqrt();
            let (sin, cos) = (omega * DT).sin_cos();
            let driven_equilibrium = equilibrium + pressure / k;
            let next_u = driven_equilibrium + (u - driven_equilibrium) * cos + v * sin / omega;
            let next_v = -(u - driven_equilibrium) * omega * sin + v * cos;
            let frame = run
                .step(
                    pressure,
                    thermal_config(),
                    transfer_controls(),
                    observations,
                    cx,
                )
                .unwrap();
            temperature += 200.0 * DT / 10.0;
            let (next_equilibrium, next_k) = pencil_at(temperature);
            let old_energy = 0.3 * next_v.powi(2) + 0.5 * k * (next_u - equilibrium).powi(2);
            let new_energy =
                0.3 * next_v.powi(2) + 0.5 * next_k * (next_u - next_equilibrium).powi(2);
            let nodes = run.motion().nodal_motion(cx).unwrap();
            near(nodes.displacement_m[7][0], next_u, 2e-12);
            near(nodes.velocity_m_s[7][0], next_v, 2e-11);
            near(frame.endpoint_observer_pressure_pa, 2.0 * next_v, 4e-11);
            near(
                frame.pressure_step.modal.input_work_j,
                pressure * (next_u - u),
                1e-12,
            );
            near(
                frame.pressure_step.swept_volume_change_m3,
                next_u - u,
                2e-12,
            );
            near(
                frame.pressure_step.modal.total_modal_energy_j,
                old_energy,
                1e-11,
            );
            near(frame.transfer.retained_modal_energy_j, new_energy, 1e-11);
            near(
                frame.transfer.parameter_work_j,
                new_energy - old_energy,
                1e-11,
            );
            assert!(frame.transfer.kinetic_residual_j < 1e-20);
            assert!(frame.transfer.elastic_residual_j < 1e-20);
            // The explicit thermal boundary supplies this heat. Signed modal
            // parameter work remains a separate account, never a caloric debit.
            near(frame.thermal.external_heat_j, 6.0 * 200.0 * DT, 1e-10);
            near(
                frame.thermal.stored_energy_change_j,
                frame.thermal.external_heat_j,
                1e-9,
            );
            near(
                run.thermal().temperature_k(CellId(0)).unwrap(),
                temperature,
                1e-11,
            );
            assert_eq!(run.motion().accepted_steps(), step);
            assert_eq!(frame.time_s, step as f64 * DT);
            assert_eq!(run.time_s(), frame.time_s);
            assert_eq!(
                run.initial_time_s() + frame.pressure_step.time_s,
                frame.time_s
            );
            (u, v) = (next_u, next_v);
        }
    });
}

#[test]
fn g3_heat_capacity_only_material_swap_changes_actual_retained_motion_and_observation() {
    let run = |cp| {
        with_motion(cp, 200.0, 0.0, |run, cx| {
            let mut pressure = 0.0;
            for _ in 0..8 {
                pressure = run
                    .step(0.0, thermal_config(), transfer_controls(), observations, cx)
                    .unwrap()
                    .endpoint_observer_pressure_pa;
            }
            (
                run.motion().nodal_motion(cx).unwrap(),
                pressure,
                run.thermal().temperature_k(CellId(0)).unwrap(),
                run.motion().valve().equilibrium().pencil().total_mass_kg,
            )
        })
    };
    let quick = run(10.0);
    let slow = run(20.0);
    assert!(quick.2 > slow.2);
    assert!((quick.0.displacement_m[7][0] - slow.0.displacement_m[7][0]).abs() > 1e-8);
    assert!((quick.0.velocity_m_s[7][0] - slow.0.velocity_m_s[7][0]).abs() > 1e-8);
    assert!((quick.1 - slow.1).abs() > 1e-8);
    near(quick.3, 6.0, 1e-12);
    assert_eq!(quick.3, slow.3);
}

#[test]
fn g4_late_source_observation_transfer_and_cancellation_refusals_keep_both_owners() {
    let clean = with_motion(10.0, 200.0, 0.0, |run, cx| {
        let frame = run
            .step(
                0.04,
                thermal_config(),
                transfer_controls(),
                observations,
                cx,
            )
            .unwrap();
        (
            run.motion().model().states().to_vec(),
            run.motion().valve().identity(),
            run.thermal().cells().clone(),
            frame.endpoint_observer_pressure_pa,
        )
    });
    with_motion(10.0, 200.0, 0.0, |run, cx| {
        let states = run.motion().model().states().to_vec();
        let identity = run.motion().valve().identity();
        let cells = run.thermal().cells().clone();
        let modes = run.motion().model().modes().to_vec();
        let unchanged = |run: &HeatedValveMotion<'_, '_, '_, '_>| {
            assert_eq!(run.motion().model().states(), states);
            assert_eq!(run.motion().model().modes(), modes);
            assert_eq!(run.motion().valve().identity(), identity);
            assert_eq!(run.thermal().cells(), &cells);
            assert_eq!(run.motion().accepted_steps(), 0);
            assert_eq!(run.time_s(), 0.0);
        };
        let visited = Cell::new(false);
        let failed = run.step(
            0.04,
            thermal_config(),
            transfer_controls(),
            |valve| {
                visited.set(true);
                assert_ne!(valve.identity(), identity);
                Err::<Vec<C64>, _>("late observer refusal")
            },
            cx,
        );
        assert!(matches!(
            failed,
            Err(HeatedValveMotionError::Observation("late observer refusal"))
        ));
        assert!(visited.get());
        unchanged(run);
        let no_work = ModalTransferControls {
            maximum_abs_parameter_work_j: 0.0,
            ..transfer_controls()
        };
        assert!(matches!(
            run.step(0.04, thermal_config(), no_work, observations, cx),
            Err(HeatedValveMotionError::Motion(_))
        ));
        unchanged(run);

        // Cancellation occurs after the complete thermal advance/reduction,
        // from the observation callback, before source transfer can publish.
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled_cx = Cx::new(
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
                thermal_config(),
                transfer_controls(),
                |valve| {
                    gate.request();
                    observations(valve)
                },
                &cancelled_cx,
            );
            assert!(failed.is_err());
        });
        unchanged(run);

        let frame = run
            .step(
                0.04,
                thermal_config(),
                transfer_controls(),
                observations,
                cx,
            )
            .unwrap();
        assert_eq!(run.motion().model().states(), clean.0);
        assert_eq!(run.motion().valve().identity(), clean.1);
        assert_eq!(run.thermal().cells(), &clean.2);
        assert_eq!(frame.endpoint_observer_pressure_pa, clean.3);
    });
}

#[test]
fn g4_nonzero_origin_keeps_one_interval_and_coarse_absolute_time_refuses_atomically() {
    with_motion(10.0, 200.0, 1.0, |run, cx| {
        let frame = run
            .step(0.0, thermal_config(), transfer_controls(), observations, cx)
            .unwrap();
        assert_eq!(run.initial_time_s(), 1.0);
        assert_eq!(frame.time_s, 1.0 + DT);
        assert_eq!(frame.pressure_step.time_s, DT);
        near(frame.thermal.external_heat_j, 6.0 * 200.0 * DT, 1e-10);
        near(
            run.thermal().temperature_k(CellId(0)).unwrap(),
            301.0,
            1e-11,
        );
    });
    with_motion(10.0, 200.0, 1e13, |run, cx| {
        let states = run.motion().model().states().to_vec();
        let cells = run.thermal().cells().clone();
        let identity = run.motion().valve().identity();
        let called = Cell::new(false);
        let result = run.step(
            0.0,
            thermal_config(),
            transfer_controls(),
            |valve| {
                called.set(true);
                observations(valve)
            },
            cx,
        );
        assert!(matches!(result, Err(HeatedValveMotionError::Clock)));
        assert!(!called.get());
        assert_eq!(run.motion().model().states(), states);
        assert_eq!(run.thermal().cells(), &cells);
        assert_eq!(run.motion().valve().identity(), identity);
        assert_eq!(run.motion().accepted_steps(), 0);
        assert_eq!(run.time_s(), 1e13);
    });
}
