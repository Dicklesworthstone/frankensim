//! Synthetic source laws: independent caloric, Poiseuille and nonlinear inverses.
use fs_airflow::graph::rheology::thermal::{
    ThermalCarreauBoundary, ThermalCarreauControls, ThermalCarreauError, ThermalCarreauGraph,
    ThermalCarreauReport, ThermalCarreauStep,
};
use fs_airflow::graph::rheology::{
    CarreauGraphError, CarreauPipeBranch, CarreauPipeGraph, FlowInjection,
};
use fs_airflow::graph::{FixedPressure, GraphError};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::fluid_transport::rheology::pipe::LaminarCarreauPipe;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::thermal::transport::storage::TabulatedHeatStorage;
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_qty::{Pressure, VolumetricFlowRate};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 719,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn law(constant: bool, mode: ViscosityShiftMode) -> ThermalCarreauYasuda {
    let (zero, infinity) = if constant { (3.0, 3.0) } else { (12.0, 0.0) };
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(zero, infinity, 0.6, 0.5, 1.0, [0.0, 1e8]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        mode,
    )
}

fn graph(
    constant: bool,
    mode: ViscosityShiftMode,
    temperature: f64,
    ceiling: f64,
    lengths: &[f64],
) -> CarreauPipeGraph {
    let original = law(constant, mode);
    CarreauPipeGraph::new(
        lengths.len() + 1,
        lengths
            .iter()
            .enumerate()
            .map(|(i, &length)| CarreauPipeBranch {
                name: format!("pipe-{i}"),
                from: i,
                to: i + 1,
                pipe: LaminarCarreauPipe::new(original, temperature, 0.01, length, ceiling)
                    .unwrap(),
            })
            .collect(),
    )
    .unwrap()
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
fn controls() -> ThermalCarreauControls {
    let mut c = ThermalCarreauControls::default();
    c.hydraulics.absolute_flow_tolerance = VolumetricFlowRate::new(2e-15);
    c.hydraulics.relative_flow_tolerance = 1e-11;
    c.hydraulics.pipe.relative_constitutive_tolerance = 1e-13;
    c.hydraulics.pipe.relative_quadrature_tolerance = 1e-11;
    c.temperature_tolerance_k = 1e-10;
    c.absolute_energy_tolerance_j = 1e-12;
    c.relative_energy_tolerance = 1e-10;
    c
}
fn advance<'a>(
    owner: &ThermalCarreauGraph<'a>,
    dt: f64,
    p: &[FixedPressure],
    q: &[FlowInjection],
    thermal: ThermalCarreauBoundary,
    c: ThermalCarreauControls,
) -> Result<ThermalCarreauStep<'a>, ThermalCarreauError> {
    with_gate(&CancelGate::new(), |cx| {
        owner.advance(dt, p, q, thermal, c, cx)
    })
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        a.is_finite()
            && b.is_finite()
            && (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn audit(r: &ThermalCarreauReport) {
    assert!(r.viscous_heat_j >= 0.0);
    assert!(r.thermal_residual_j.abs() <= r.thermal_residual_allowance_j);
    close(
        r.stored_energy_change_j,
        r.viscous_heat_j + r.prescribed_heat_j + r.bath_heat_j,
        0.0,
        r.thermal_residual_allowance_j + 1e-12,
    );
    close(
        r.combined_energy_residual_j,
        r.thermal_residual_j - r.hydraulic_work_residual_j,
        1e-10,
        1e-12,
    );
    close(
        r.hydraulic_work_residual_j,
        r.duration_s * r.midpoint_hydraulics.power_balance_residual_w
            - r.constitutive_residual_work_j
            - r.continuity_residual_work_j,
        1e-10,
        1e-12,
    );
    close(
        r.viscous_heat_j,
        r.duration_s * r.midpoint_hydraulics.dissipated_power_w,
        0.0,
        0.0,
    );
}
fn shift(temperature: f64) -> f64 {
    ((12_000.0 / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K) * (1.0 / temperature - 1.0 / 300.0)).exp()
}
fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn g0_one_inventory_thermal_owner_refuses_distinct_regional_temperatures() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 20.0], [360.0, 20.0]]).unwrap();
    let uniform = graph(true, ViscosityShiftMode::Vertical, 300.0, 1e5, &[0.5, 0.7]);
    let regional = with_gate(&CancelGate::new(), |cx| {
        uniform.at_temperatures(&[290.0, 330.0], cx)
    })
    .unwrap();
    // Regional hydraulic viscosity is admissible, but this caloric owner has
    // only one inventory temperature and cannot discard the second region.
    let hydraulics = with_gate(&CancelGate::new(), |cx| {
        regional.solve(
            &[boundary(0, 1000.0), boundary(2, 0.0)],
            controls().hydraulics,
            cx,
        )
    })
    .unwrap();
    assert!(hydraulics.dissipated_power_w > 0.0);
    assert!(matches!(
        ThermalCarreauGraph::new(regional, &storage, 0.0),
        Err(ThermalCarreauError::Hydraulic(CarreauGraphError::Graph(
            GraphError::InvalidInput(_)
        )))
    ));
    assert!(ThermalCarreauGraph::new(uniform, &storage, 0.0).is_ok());
}

#[test]
fn g1_cooling_has_exact_caloric_midpoint_and_second_order_time_refinement() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 20.0], [360.0, 20.0]]).unwrap();
    let original = ThermalCarreauGraph::new(
        graph(
            false,
            ViscosityShiftMode::VerticalAndRate,
            320.0,
            1e5,
            &[0.5],
        ),
        &storage,
        2.0,
    )
    .unwrap();
    let ports = [boundary(0, 0.0), boundary(1, 0.0)];
    let bath = ThermalCarreauBoundary {
        bath_conductance_w_k: 4.0,
        ..quiet()
    };
    let step = advance(&original, 1.0, &ports, &[], bath, controls()).unwrap();
    let expected = (18.0 * 320.0 + 4.0 * 300.0) / 22.0;
    close(step.next.temperature_k(), expected, 0.0, 2e-10);
    assert_eq!(step.next.time_s(), 3.0);
    assert_eq!(original.time_s(), 2.0);
    assert_eq!(original.temperature_k(), 320.0);
    assert_eq!(step.report.viscous_heat_j, 0.0);
    assert_eq!(step.report.reservoir_work_j, 0.0);
    audit(&step.report);
    let exact = 300.0 + 20.0 * (-0.4_f64).exp();
    let mut errors = Vec::new();
    for n in [1, 2, 4] {
        let mut current = original.clone();
        let c = ThermalCarreauControls {
            maximum_temperature_change_k: 10.0,
            ..controls()
        };
        for _ in 0..n {
            current = advance(&current, 2.0 / f64::from(n), &ports, &[], bath, c)
                .unwrap()
                .next;
        }
        errors.push((current.temperature_k() - exact).abs());
    }
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.9 && pair[0] / pair[1] < 4.2);
    }
}

#[test]
fn g1_arrhenius_poiseuille_self_heating_crosses_an_original_caloric_knot() {
    let storage = TabulatedHeatStorage::try_new(&[
        [280.0, 0.02],
        [300.0, 0.02],
        [300.5, 0.04],
        [360.0, 0.04],
    ])
    .unwrap();
    let owner = ThermalCarreauGraph::new(
        graph(true, ViscosityShiftMode::Vertical, 300.0, 1e5, &[0.5]),
        &storage,
        0.0,
    )
    .unwrap();
    let pressure = 10_000.0;
    let dt = 0.1;
    let power = |t: f64| {
        core::f64::consts::PI * 0.01_f64.powi(4) * pressure * pressure
            / (8.0 * 3.0 * shift(t) * 0.5)
    };
    let energy = |t: f64| {
        let d = t - 300.0;
        if d <= 0.5 {
            0.02 * d + 0.02 * d * d
        } else {
            0.015 + 0.04 * (d - 0.5)
        }
    };
    let expected = root(300.0, 305.0, |t| energy(t) - dt * power(0.5 * (300.0 + t)));
    let step = advance(
        &owner,
        dt,
        &[boundary(0, pressure), boundary(1, 0.0)],
        &[],
        quiet(),
        controls(),
    )
    .unwrap();
    assert!(step.next.temperature_k() > 300.5);
    close(step.next.temperature_k(), expected, 0.0, 3e-9);
    close(
        step.report.stored_energy_change_j,
        energy(expected),
        2e-8,
        1e-11,
    );
    close(
        step.report.viscous_heat_j,
        dt * power(step.report.midpoint_temperature_k),
        2e-10,
        1e-12,
    );
    assert_eq!(
        step.next.graph().branches()[0].pipe.temperature_k(),
        step.next.temperature_k()
    );
    assert_eq!(
        step.report.midpoint_hydraulics.branches[0]
            .profile
            .pipe()
            .temperature_k(),
        step.report.midpoint_temperature_k
    );
    assert_eq!(
        step.next.graph().branches()[0].pipe.law(),
        owner.graph().branches()[0].pipe.law()
    );
    audit(&step.report);
}

// Exact inverse of tau=12*gamma/sqrt(1+0.6*gamma), independently of native queries.
fn inverse(stress: f64) -> f64 {
    let z = 0.6 * stress / 12.0;
    (z * z + z * (z * z + 4.0).sqrt()) / 1.2
}
fn nonlinear_flow(pressure: f64, length: f64, temperature: f64, mode: ViscosityShiftMode) -> f64 {
    let factor = shift(temperature);
    let wall = pressure * 0.01 / (2.0 * length);
    let n = 2048;
    let mut integral = 0.0;
    for i in 0..=n {
        let x = f64::from(i) / f64::from(n);
        let rate = match mode {
            ViscosityShiftMode::Vertical => inverse(wall * x / factor),
            ViscosityShiftMode::VerticalAndRate => inverse(wall * x) / factor,
        };
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * rate * x * x;
    }
    core::f64::consts::PI * 0.01_f64.powi(3) * integral / (3.0 * f64::from(n))
}

#[test]
fn g1_original_nonlinear_series_network_retains_both_thermal_shift_calibrations() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 0.02], [360.0, 0.02]]).unwrap();
    let pressure = 6000.0;
    let mut temperatures = Vec::new();
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        let owner =
            ThermalCarreauGraph::new(graph(false, mode, 300.0, 1e5, &[0.3, 0.7]), &storage, 0.0)
                .unwrap();
        let expected = root(300.0, 305.0, |t| {
            0.02 * (t - 300.0) - pressure * nonlinear_flow(pressure, 1.0, 0.5 * (300.0 + t), mode)
        });
        let step = advance(
            &owner,
            1.0,
            &[boundary(0, pressure), boundary(2, 0.0)],
            &[],
            quiet(),
            controls(),
        )
        .unwrap();
        close(step.next.temperature_k(), expected, 0.0, 2e-8);
        let expected_flow = nonlinear_flow(pressure, 1.0, step.report.midpoint_temperature_k, mode);
        for branch in &step.report.midpoint_hydraulics.branches {
            close(
                branch.profile.flow().volumetric_flow_m3_s,
                expected_flow,
                2e-8,
                1e-14,
            );
        }
        close(
            step.report.midpoint_hydraulics.pressures[1].value(),
            0.7 * pressure,
            2e-10,
            1e-9,
        );
        audit(&step.report);
        temperatures.push(step.next.temperature_k());
    }
    // Vertical-and-rate shifting must not collapse to vertical scaling alone.
    assert!((temperatures[0] - temperatures[1]).abs() > 1e-5);
}

#[test]
fn g1_flow_controlled_heating_preserves_signed_supply_work_without_double_deposition() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 0.01], [360.0, 0.01]]).unwrap();
    let owner = ThermalCarreauGraph::new(
        graph(true, ViscosityShiftMode::Vertical, 300.0, 1e5, &[0.5]),
        &storage,
        0.0,
    )
    .unwrap();
    let q = 2e-6;
    let injection = FlowInjection {
        node: 0,
        flow: VolumetricFlowRate::new(q),
    };
    let at_reservoir = FlowInjection {
        node: 1,
        flow: VolumetricFlowRate::new(5e-6),
    };
    let ports = [boundary(1, 2000.0)];
    let native = advance(&owner, 1.0, &ports, &[injection], quiet(), controls()).unwrap();
    let step = advance(
        &owner,
        1.0,
        &ports,
        &[injection, at_reservoir],
        quiet(),
        controls(),
    )
    .unwrap();
    let power =
        |t: f64| 8.0 * 3.0 * shift(t) * 0.5 * q * q / (core::f64::consts::PI * 0.01_f64.powi(4));
    let expected = root(300.0, 301.0, |t| {
        0.01 * (t - 300.0) - power(0.5 * (300.0 + t))
    });
    close(step.next.temperature_k(), expected, 0.0, 2e-9);
    assert!(step.report.reservoir_work_j < 0.0 && step.report.injection_work_j > 0.0);
    assert_eq!(step.next.temperature_k(), native.next.temperature_k());
    assert_eq!(step.report.viscous_heat_j, native.report.viscous_heat_j);
    close(
        step.report.reservoir_work_j + step.report.injection_work_j,
        step.report.viscous_heat_j,
        2e-10,
        1e-12,
    );
    audit(&step.report);
}

#[test]
fn g0_final_source_admission_and_caloric_limits_preserve_the_initial_owner() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 10.0], [360.0, 10.0]]).unwrap();
    let owner = ThermalCarreauGraph::new(
        graph(true, ViscosityShiftMode::Vertical, 300.0, 5.14, &[0.5]),
        &storage,
        4.0,
    )
    .unwrap();
    let ports = [boundary(0, 1500.0), boundary(1, 0.0)];
    let forcing = ThermalCarreauBoundary {
        prescribed_power_w: 20.0,
        ..quiet()
    };
    let c = ThermalCarreauControls {
        maximum_temperature_change_k: 3.0,
        ..controls()
    };
    // Midpoint trials at <=301.5 K fit the rate ceiling; the ~302 K physical
    // endpoint does not. Admitting only the midpoint would publish bad history.
    with_gate(&CancelGate::new(), |cx| {
        assert!(
            graph(true, ViscosityShiftMode::Vertical, 301.5, 5.14, &[0.5])
                .solve(&ports, c.hydraulics, cx)
                .is_ok()
        );
        assert!(
            graph(true, ViscosityShiftMode::Vertical, 302.0, 5.14, &[0.5])
                .solve(&ports, c.hydraulics, cx)
                .is_err()
        );
    });
    assert!(matches!(
        advance(&owner, 1.0, &ports, &[], forcing, c),
        Err(ThermalCarreauError::Hydraulic(_))
    ));
    assert_eq!(owner.temperature_k(), 300.0);
    assert_eq!(owner.time_s(), 4.0);
    let small = advance(&owner, 0.01, &ports, &[], forcing, c).unwrap();
    let retry = advance(&owner, 0.01, &ports, &[], forcing, c).unwrap();
    assert_eq!(small.next.temperature_k(), retry.next.temperature_k());
    assert_eq!(small.report.viscous_heat_j, retry.report.viscous_heat_j);

    let limited = TabulatedHeatStorage::try_new(&[[299.0, 10.0], [300.2, 10.0]]).unwrap();
    let owner = ThermalCarreauGraph::new(
        graph(true, ViscosityShiftMode::Vertical, 300.0, 1e5, &[0.5]),
        &limited,
        0.0,
    )
    .unwrap();
    assert!(matches!(
        advance(&owner, 1.0, &ports, &[], forcing, c),
        Err(ThermalCarreauError::UnbracketedTemperature)
    ));
    assert_eq!(owner.temperature_k(), 300.0);
}

#[test]
fn g0_stiff_equilibrium_overshoot_budget_and_cancellation_publish_no_next_state() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 20.0], [360.0, 20.0]]).unwrap();
    let owner = ThermalCarreauGraph::new(
        graph(true, ViscosityShiftMode::Vertical, 320.0, 1e5, &[0.5]),
        &storage,
        3.0,
    )
    .unwrap();
    let ports = [boundary(0, 0.0), boundary(1, 0.0)];
    let bath = ThermalCarreauBoundary {
        bath_conductance_w_k: 4.0,
        ..quiet()
    };
    let c = ThermalCarreauControls {
        maximum_temperature_change_k: 40.0,
        ..controls()
    };
    assert!(matches!(
        advance(&owner, 20.0, &ports, &[], bath, c),
        Err(ThermalCarreauError::EquilibriumCrossing)
    ));
    let short = ThermalCarreauControls {
        maximum_hydraulic_solves: 2,
        ..c
    };
    assert!(matches!(
        advance(&owner, 1.0, &ports, &[], bath, short),
        Err(ThermalCarreauError::HydraulicSolveBudget)
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        with_gate(&gate, |cx| owner.advance(1.0, &ports, &[], bath, c, cx)),
        Err(ThermalCarreauError::Hydraulic(CarreauGraphError::Graph(
            GraphError::Interrupted
        )))
    ));
    assert_eq!(owner.temperature_k(), 320.0);
    assert_eq!(owner.time_s(), 3.0);
    let first = advance(&owner, 1.0, &ports, &[], bath, c).unwrap();
    let retry = advance(&owner, 1.0, &ports, &[], bath, c).unwrap();
    assert_eq!(first.next.temperature_k(), retry.next.temperature_k());
    assert_eq!(first.report.hydraulic_solves, retry.report.hydraulic_solves);
}
