//! Independent finite-volume momentum and caloric references on synthetic laws.

use super::*;
use crate::graph::rheology::{CarreauPipeBranch, CarreauPipeGraph};
use crate::graph::tests::with_cx;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::fluid_transport::rheology::pipe::LaminarCarreauPipe;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_qty::VolumetricFlowRate;

const RADIUS: f64 = 0.01;
const LENGTH: f64 = 0.5;
const DENSITY: f64 = 1000.0;
const CAPACITY: f64 = 0.02;

fn near(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        a.is_finite()
            && b.is_finite()
            && (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

fn model(compliance: Option<f64>) -> TransientCarreauGraph {
    rheological_model(compliance, false)
}

fn rheological_model(compliance: Option<f64>, shear_thinning: bool) -> TransientCarreauGraph {
    let source = ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            0.2,
            if shear_thinning { 0.01 } else { 0.2 },
            0.8,
            0.5,
            1.2,
            [0.0, 1e7],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 9000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        if shear_thinning {
            ViscosityShiftMode::VerticalAndRate
        } else {
            ViscosityShiftMode::Vertical
        },
    );
    let graph = CarreauPipeGraph::new(
        2,
        vec![CarreauPipeBranch {
            name: "original".into(),
            from: 0,
            to: 1,
            pipe: LaminarCarreauPipe::new(source, 300.0, RADIUS, LENGTH, 1e5).unwrap(),
        }],
    )
    .unwrap();
    TransientCarreauGraph::new(
        graph,
        DENSITY,
        &[2],
        vec![
            None,
            compliance.map(|compliance_m3_pa| {
                crate::graph::rheology::transient::LinearHydraulicStorage {
                    compliance_m3_pa,
                    minimum_pressure: Pressure::new(-1e5),
                    maximum_pressure: Pressure::new(1e5),
                }
            }),
        ],
    )
    .unwrap()
}

fn storage() -> TabulatedHeatStorage {
    TabulatedHeatStorage::try_new(&[[280.0, CAPACITY], [360.0, CAPACITY]]).unwrap()
}
fn boundary(node: usize, pressure: f64) -> FixedPressure {
    FixedPressure {
        node,
        pressure: Pressure::new(pressure),
    }
}
fn quiet() -> ThermalCarreauBoundary {
    ThermalCarreauBoundary {
        prescribed_power_w: 0.0,
        bath_conductance_w_k: 0.0,
        bath_temperature_k: 300.0,
    }
}
fn controls() -> TransientThermalCarreauControls {
    let mut controls = TransientThermalCarreauControls::default();
    controls.hydraulics.absolute_flow_tolerance = VolumetricFlowRate::new(2e-15);
    controls.hydraulics.relative_flow_tolerance = 1e-11;
    controls.hydraulics.radial.absolute_velocity_tolerance_m_s = 1e-14;
    controls.hydraulics.radial.relative_velocity_tolerance = 1e-12;
    controls.temperature_tolerance_k = 2e-11;
    controls.absolute_energy_tolerance_j = 1e-12;
    controls.relative_energy_tolerance = 1e-10;
    controls
}

fn history(
    model: &TransientCarreauGraph,
    velocities: [f64; 2],
    pressure: Option<f64>,
) -> TransientCarreauState {
    with_cx(|cx| {
        model.state_from_velocities(
            vec![velocities.to_vec()],
            &[None, pressure.map(Pressure::new)],
            2.0,
            cx,
        )
    })
    .unwrap()
}
fn snapshot(state: &TransientCarreauState) -> (u64, Vec<u64>, Vec<Option<u64>>) {
    (
        state.accepted_time_s().to_bits(),
        state.branches()[0]
            .velocities_m_s()
            .iter()
            .map(|value| value.to_bits())
            .collect(),
        state
            .storage_pressures()
            .iter()
            .map(|pressure| pressure.map(|value| value.value().to_bits()))
            .collect(),
    )
}

// Two equal-width annuli give masses rho*A_i*L. Their interior shear stiffness
// is 2*pi*mu*L; the outer no-slip half-cell gives 8*pi*mu*L. The optional third
// unknown is the actual midpoint pressure of a compliant end reservoir.
fn reference_step(
    old: [f64; 2],
    pressure: f64,
    old_storage: Option<(f64, f64)>,
    dt: f64,
    temperature: f64,
) -> ([f64; 2], Option<f64>, f64, f64, f64) {
    let area = core::f64::consts::PI * (0.5 * RADIUS).powi(2);
    let areas = [area, 3.0 * area];
    let masses = areas.map(|value| DENSITY * LENGTH * value);
    let viscosity =
        0.2 * ((9000.0 / 8.314_462_618_153_24) * (1.0 / temperature - 1.0 / 300.0)).exp();
    let face = 2.0 * core::f64::consts::PI * viscosity * LENGTH;
    let wall = 4.0 * face;
    let stiffness = [[face, -face], [-face, face + wall]];
    let n = if old_storage.is_some() { 3 } else { 2 };
    let mut matrix = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for i in 0..2 {
        matrix[i][i] = masses[i];
        rhs[i] = masses[i] * old[i] + dt * areas[i] * pressure;
        for j in 0..2 {
            matrix[i][j] += 0.5 * dt * stiffness[i][j];
            rhs[i] -= 0.5 * dt * stiffness[i][j] * old[j];
        }
        if old_storage.is_some() {
            matrix[i][2] = dt * areas[i];
            matrix[2][i] = -0.5 * areas[i];
        }
    }
    if let Some((compliance, old_pressure)) = old_storage {
        matrix[2][2] = 2.0 * compliance / dt;
        rhs[2] =
            2.0 * compliance * old_pressure / dt + 0.5 * (areas[0] * old[0] + areas[1] * old[1]);
    }
    for column in 0..n {
        for row in column + 1..n {
            let ratio = matrix[row][column] / matrix[column][column];
            for j in column..n {
                matrix[row][j] -= ratio * matrix[column][j];
            }
            rhs[row] -= ratio * rhs[column];
        }
    }
    let mut solved = [0.0; 3];
    for i in (0..n).rev() {
        let known: f64 = (i + 1..n).map(|j| matrix[i][j] * solved[j]).sum();
        solved[i] = (rhs[i] - known) / matrix[i][i];
    }
    let next = [solved[0], solved[1]];
    let mid = [0.5 * (old[0] + next[0]), 0.5 * (old[1] + next[1])];
    let heat = dt * (face * (mid[1] - mid[0]).powi(2) + wall * mid[1].powi(2));
    let kinetic: f64 = (0..2)
        .map(|i| 0.5 * masses[i] * (next[i] * next[i] - old[i] * old[i]))
        .sum();
    let work = dt * pressure * (areas[0] * mid[0] + areas[1] * mid[1]);
    (
        next,
        old_storage.map(|(_, old_pressure)| 2.0 * solved[2] - old_pressure),
        heat,
        kinetic,
        work,
    )
}

fn reference_temperature(
    old: [f64; 2],
    pressure: f64,
    old_storage: Option<(f64, f64)>,
    dt: f64,
) -> f64 {
    let mut lo = 300.0;
    let mut hi = 305.0;
    let residual = |temperature: f64| {
        CAPACITY * (temperature - 300.0)
            - reference_step(old, pressure, old_storage, dt, 0.5 * (300.0 + temperature)).2
    };
    assert!(residual(lo) <= 0.0 && residual(hi) > 0.0);
    for _ in 0..64 {
        let midpoint = 0.5 * (lo + hi);
        if residual(midpoint) < 0.0 {
            lo = midpoint;
        } else {
            hi = midpoint;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_startup_coastdown_and_recovering_head_follow_original_heat_and_radial_history() {
    let caloric = storage();
    for (old, pressure) in [([0.0, 0.0], 600.0), ([0.3, 0.1], 0.0), ([0.3, 0.1], -100.0)] {
        let model = model(None);
        let state = history(&model, old, None);
        let before = snapshot(&state);
        let mut workspace = model.workspace().unwrap();
        let owner = TransientThermalCarreauGraph::new(model, &caloric).unwrap();
        let dt = 0.05;
        let expected_temperature = reference_temperature(old, pressure, None, dt);
        let expected = reference_step(
            old,
            pressure,
            None,
            dt,
            0.5 * (300.0 + expected_temperature),
        );
        let result = with_cx(|cx| {
            owner.advance_with_heads(
                &state,
                &mut workspace,
                dt,
                &[boundary(0, 0.0), boundary(1, 0.0)],
                &[],
                &[Pressure::new(pressure)],
                quiet(),
                controls(),
                cx,
            )
        })
        .unwrap();
        assert_eq!(snapshot(&state), before);
        assert_eq!(owner.temperature_k(), 300.0);
        assert!(core::ptr::eq(result.next.storage(), &caloric));
        near(
            result.next.temperature_k(),
            expected_temperature,
            0.0,
            3e-10,
        );
        for (&actual, expected) in result.state.branches()[0]
            .velocities_m_s()
            .iter()
            .zip(expected.0)
        {
            near(actual, expected, 3e-10, 2e-13);
        }
        near(result.report.viscous_heat_j, expected.2, 3e-10, 1e-14);
        near(
            result.report.hydraulics.kinetic_energy_change_j,
            expected.3,
            3e-10,
            1e-14,
        );
        near(
            result.report.hydraulics.branch_head_work_j,
            expected.4,
            3e-10,
            1e-14,
        );
        near(
            result.report.stored_energy_change_j,
            expected.2,
            3e-9,
            1e-12,
        );
        assert!(result.report.viscous_heat_j > 0.0);
        assert!(
            result.report.thermal_residual_j.abs() <= result.report.thermal_residual_allowance_j
        );
        assert!(
            result.report.combined_energy_residual_j.abs()
                <= result.report.combined_residual_allowance_j
        );
        if pressure > 0.0 {
            assert!(result.report.hydraulics.kinetic_energy_change_j > 0.0);
            assert!(result.report.hydraulics.branch_head_work_j > result.report.viscous_heat_j);
        } else if pressure == 0.0 {
            assert_eq!(result.report.hydraulics.branch_head_work_j, 0.0);
            near(
                result.report.viscous_heat_j,
                -result.report.hydraulics.kinetic_energy_change_j,
                1e-10,
                1e-13,
            );
        } else {
            assert!(result.report.hydraulics.branch_head_work_j < 0.0);
        }
        // The returned temperature owner accepts the retained state and SAME
        // workspace on the next interval; no steady-profile reset intervenes.
        with_cx(|cx| {
            result.next.advance(
                &result.state,
                &mut workspace,
                0.01,
                &[boundary(0, 0.0), boundary(1, 0.0)],
                &[],
                quiet(),
                controls(),
                cx,
            )
        })
        .unwrap();
    }
}

#[test]
fn g1_compliance_storage_is_not_deposited_as_heat() {
    let caloric = storage();
    let compliance = 2e-10;
    let old_pressure = 12.0;
    let old = [0.0, 0.0];
    let model = model(Some(compliance));
    let state = history(&model, old, Some(old_pressure));
    let mut workspace = model.workspace().unwrap();
    let owner = TransientThermalCarreauGraph::new(model, &caloric).unwrap();
    let dt = 0.01;
    let pressure = 600.0;
    let expected_temperature =
        reference_temperature(old, pressure, Some((compliance, old_pressure)), dt);
    let expected = reference_step(
        old,
        pressure,
        Some((compliance, old_pressure)),
        dt,
        0.5 * (300.0 + expected_temperature),
    );
    let result = with_cx(|cx| {
        owner.advance(
            &state,
            &mut workspace,
            dt,
            &[boundary(0, pressure)],
            &[],
            quiet(),
            controls(),
            cx,
        )
    })
    .unwrap();
    let next_pressure = expected.1.unwrap();
    let stored = 0.5 * compliance * (next_pressure * next_pressure - old_pressure * old_pressure);
    near(
        result.next.temperature_k(),
        expected_temperature,
        0.0,
        3e-10,
    );
    near(
        result.state.storage_pressures()[1].unwrap().value(),
        next_pressure,
        3e-10,
        1e-9,
    );
    near(
        result.report.hydraulics.storage_energy_change_j,
        stored,
        3e-10,
        1e-14,
    );
    near(
        result.report.hydraulics.reservoir_work_j,
        expected.4,
        3e-10,
        1e-14,
    );
    near(result.report.viscous_heat_j, expected.2, 3e-10, 1e-14);
    near(
        result.report.stored_energy_change_j,
        expected.2,
        3e-8,
        1e-12,
    );
    assert!(stored > result.report.viscous_heat_j);
    assert!(result.report.hydraulics.reservoir_work_j > result.report.viscous_heat_j);
}

#[test]
fn g1_zero_flow_cooling_and_balanced_steady_limit_retain_native_behavior() {
    let caloric = storage();
    let model = model(None);
    let state = history(&model, [0.0, 0.0], None);
    let mut workspace = model.workspace().unwrap();
    let owner = TransientThermalCarreauGraph::new(model, &caloric).unwrap();
    let ports = [boundary(0, 0.0), boundary(1, 0.0)];
    let idle = with_cx(|cx| {
        owner.advance(
            &state,
            &mut workspace,
            1.0,
            &ports,
            &[],
            quiet(),
            controls(),
            cx,
        )
    })
    .unwrap();
    assert_eq!(idle.next.temperature_k(), 300.0);
    assert_eq!(idle.report.viscous_heat_j, 0.0);
    assert_eq!(idle.report.hydraulic_solves, 1);
    let bath = ThermalCarreauBoundary {
        bath_conductance_w_k: 0.002,
        bath_temperature_k: 280.0,
        ..quiet()
    };
    let cooled = with_cx(|cx| {
        owner.advance(
            &state,
            &mut workspace,
            1.0,
            &ports,
            &[],
            bath,
            controls(),
            cx,
        )
    })
    .unwrap();
    near(
        cooled.next.temperature_k(),
        ((CAPACITY - 0.001) * 300.0 + 0.002 * 280.0) / (CAPACITY + 0.001),
        0.0,
        3e-10,
    );
    assert_eq!(cooled.report.viscous_heat_j, 0.0);
    assert!(
        cooled.state.branches()[0]
            .velocities_m_s()
            .iter()
            .all(|&u| u == 0.0)
    );

    // Exact two-annulus steady solution K*u=A*dp. A compensating declared
    // cooling power holds its temperature, isolating the steady mechanical limit.
    let pressure = 600.0;
    let area = core::f64::consts::PI * (0.5 * RADIUS).powi(2);
    let face = 2.0 * core::f64::consts::PI * 0.2 * LENGTH;
    let outer = area * pressure / face;
    let velocities = [outer + area * pressure / face, outer];
    let state = history(owner.graph(), velocities, None);
    let mut native_state = state.clone();
    let mut native_workspace = owner.graph().workspace().unwrap();
    let dt = 0.0625;
    let native = with_cx(|cx| {
        owner.graph().step(
            &mut native_state,
            &mut native_workspace,
            &[boundary(0, pressure), boundary(1, 0.0)],
            &[],
            dt,
            controls().hydraulics,
            cx,
        )
    })
    .unwrap();
    let balanced = ThermalCarreauBoundary {
        prescribed_power_w: -native.viscous_heat_j / dt,
        ..quiet()
    };
    let result = with_cx(|cx| {
        owner.advance(
            &state,
            &mut workspace,
            dt,
            &[boundary(0, pressure), boundary(1, 0.0)],
            &[],
            balanced,
            controls(),
            cx,
        )
    })
    .unwrap();
    assert_eq!(result.next.temperature_k(), 300.0);
    assert_eq!(snapshot(&result.state), snapshot(&native_state));
    near(
        result.report.hydraulics.reservoir_work_j,
        result.report.viscous_heat_j,
        1e-10,
        1e-13,
    );
    for (&actual, expected) in result.state.branches()[0]
        .velocities_m_s()
        .iter()
        .zip(velocities)
    {
        near(actual, expected, 1e-10, 1e-13);
    }
}

#[test]
fn g3_nonlinear_heat_feedback_uses_the_actual_original_midpoint_trajectory() {
    let caloric = storage();
    let model = rheological_model(None, true);
    let state = history(&model, [0.3, 0.1], None);
    let mut workspace = model.workspace().unwrap();
    let owner = TransientThermalCarreauGraph::new(model, &caloric).unwrap();
    let ports = [boundary(0, 200.0), boundary(1, 0.0)];
    let dt = 0.05;
    let result = with_cx(|cx| {
        owner.advance(
            &state,
            &mut workspace,
            dt,
            &ports,
            &[],
            quiet(),
            controls(),
            cx,
        )
    })
    .unwrap();
    assert!(result.next.temperature_k() > owner.temperature_k());
    let midpoint = with_cx(|cx| {
        owner
            .graph()
            .at_temperature(result.report.midpoint_temperature_k, cx)
    })
    .unwrap();
    let mut original_state = state.clone();
    let mut original_workspace = owner.graph().workspace().unwrap();
    let original = with_cx(|cx| {
        midpoint.step(
            &mut original_state,
            &mut original_workspace,
            &ports,
            &[],
            dt,
            controls().hydraulics,
            cx,
        )
    })
    .unwrap();
    assert_eq!(snapshot(&result.state), snapshot(&original_state));
    assert_eq!(result.report.hydraulics.branches, original.branches);
    assert_eq!(result.report.viscous_heat_j, original.viscous_heat_j);
    near(
        caloric
            .energy_change_j(300.0, result.next.temperature_k())
            .unwrap(),
        original.viscous_heat_j,
        0.0,
        result.report.thermal_residual_allowance_j,
    );
    let mut frozen_state = state.clone();
    let frozen = with_cx(|cx| {
        owner.graph().step(
            &mut frozen_state,
            &mut original_workspace,
            &ports,
            &[],
            dt,
            controls().hydraulics,
            cx,
        )
    })
    .unwrap();
    assert_ne!(snapshot(&result.state), snapshot(&frozen_state));
    assert_ne!(original.viscous_heat_j, frozen.viscous_heat_j);
}

#[test]
fn g4_refused_thermal_trials_preserve_temperature_history_and_retry() {
    let caloric = storage();
    let model = model(None);
    let state = history(&model, [0.3, 0.1], None);
    let before = snapshot(&state);
    let mut workspace = model.workspace().unwrap();
    let owner = TransientThermalCarreauGraph::new(model, &caloric).unwrap();
    let ports = [boundary(0, 0.0), boundary(1, 0.0)];
    for bad in [
        TransientThermalCarreauControls {
            maximum_hydraulic_solves: 1,
            ..controls()
        },
        TransientThermalCarreauControls {
            maximum_total_branch_solves: 1,
            ..controls()
        },
        TransientThermalCarreauControls {
            maximum_temperature_change_k: 1e-9,
            ..controls()
        },
    ] {
        assert!(
            with_cx(|cx| owner.advance(
                &state,
                &mut workspace,
                0.05,
                &ports,
                &[],
                quiet(),
                bad,
                cx
            ))
            .is_err()
        );
        assert_eq!(snapshot(&state), before);
        assert_eq!(owner.temperature_k(), 300.0);
    }
    let gate = CancelGate::new();
    gate.request();
    let cancelled = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 719,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        owner.advance(
            &state,
            &mut workspace,
            0.05,
            &ports,
            &[],
            quiet(),
            controls(),
            &cx,
        )
    });
    assert!(matches!(
        cancelled,
        Err(ThermalCarreauError::Hydraulic(
            crate::graph::rheology::CarreauGraphError::Graph(GraphError::Interrupted)
        ))
    ));
    assert_eq!(snapshot(&state), before);
    let retry = with_cx(|cx| {
        owner.advance(
            &state,
            &mut workspace,
            0.05,
            &ports,
            &[],
            quiet(),
            controls(),
            cx,
        )
    })
    .unwrap();
    let mut fresh = owner.graph().workspace().unwrap();
    let reference = with_cx(|cx| {
        owner.advance(
            &state,
            &mut fresh,
            0.05,
            &ports,
            &[],
            quiet(),
            controls(),
            cx,
        )
    })
    .unwrap();
    assert_eq!(retry.next.temperature_k(), reference.next.temperature_k());
    assert_eq!(snapshot(&retry.state), snapshot(&reference.state));
    assert_eq!(retry.report.viscous_heat_j, reference.report.viscous_heat_j);
}
