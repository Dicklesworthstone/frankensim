use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{NsfSpec, WallThermal};
use fs_material::gas::GasSpec;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

fn control(tolerance: f64, maximum_endpoint_evaluations: usize) -> NsfCouplingControl {
    NsfCouplingControl {
        maximum_endpoint_evaluations,
        force_relative_tolerance: tolerance,
        velocity_relative_tolerance: tolerance,
        caloric_relative_tolerance: tolerance,
        energy_relative_tolerance: tolerance,
    }
}

fn spec(dims: [usize; 3], workers: usize) -> NsfSpec {
    NsfSpec {
        dims,
        periodic: [true, false, true],
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 300.0,
        pressure_pa: 100_000.0,
        dx_m: 0.001,
        temperature_window_k: [250.0, 400.0],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.3,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012 N2 ratio; endpoint control regression",
        maximum_cells: dims.into_iter().product(),
        wall: WallThermal::Adiabatic,
        workers,
        transport_points: 17,
    }
}

fn convecting_gas(workers: usize) -> NsfGrid3<D3q27> {
    let dims = [16, 12, 4];
    NsfGrid3::new(
        spec(dims, workers),
        |_| false,
        |x| {
            let px = core::f64::consts::TAU * (x[0] as f64 + 0.5) / dims[0] as f64;
            let py = core::f64::consts::PI * (x[1] as f64 + 0.5) / dims[1] as f64;
            let temperature = 300.0 + 25.0 * px.sin() * py.sin();
            (
                300.0 / temperature,
                [
                    0.045 * py.sin(),
                    0.01 * px.sin() * py.sin(),
                    0.006 * px.cos(),
                ],
                temperature,
            )
        },
    )
    .unwrap()
    .with_gravity_si([0.0, -9.80665, 0.0])
    .unwrap()
}

#[test]
fn accepted_report_measures_the_published_thermal_flow_endpoint() {
    let controls = control(1e-9, 64);
    let mut serial = convecting_gas(1).with_coupling_control(controls).unwrap();
    let mut pooled = convecting_gas(3).with_coupling_control(controls).unwrap();
    for _ in 0..4 {
        let old_e = serial.e.clone();
        let old_div = serial.div.clone();
        let (old_rho, old_momentum) = serial.macro_fields(&serial.f);
        let old_u = serial.velocity_field(&old_rho, &old_momentum, &serial.force);
        let ledger = serial.step().unwrap();
        assert_eq!(ledger, pooled.step().unwrap());
        assert_eq!(serial.state_hash(), pooled.state_hash());
        assert_eq!(serial.last_coupling_report(), pooled.last_coupling_report());
        let report = serial.last_coupling_report().unwrap();
        assert!(controls.accepts(report), "{report:?}");
        assert!(report.maximum_tolerance_ratio <= 1.0);
        assert_eq!(report.endpoint_evaluations, report.corrector_passes + 1);
        assert_eq!(
            report.residual_probes,
            report.endpoint_evaluations - COUPLING_PASSES
        );
        assert!(report.endpoint_evaluations <= controls.maximum_endpoint_evaluations);
        assert!(ledger.residual_rel < 2e-12);

        // Re-evaluate the PUBLISHED state, independent of the iteration history.
        // In particular the energy flux uses its held physical velocity and
        // Guo force, not the velocity of a further hypothetical corrector.
        let (rho, momentum) = serial.macro_fields(&serial.f);
        let u = serial.velocity_field(&rho, &momentum, &serial.force);
        let force = serial.force_field(&rho, &serial.theta, None);
        let required_u = serial.velocity_field(&rho, &momentum, &force);
        let material = serial.material_field(&rho, &serial.theta).unwrap();
        let stress = serial.stress_field(
            &serial.f,
            &serial.force,
            &rho,
            &serial.theta,
            &material,
            &[],
        );
        let (div, _, _) = serial
            .flux_divergence(&serial.e, &rho, &u, &serial.theta, &stress, &[], None)
            .unwrap();
        let cv = 1.0 / (serial.law.gamma - 1.0);
        let mut measured = [0.0_f64; 4];
        for i in 0..serial.n {
            for a in 0..3 {
                measured[0] =
                    measured[0].max((serial.force[i][a] - force[i][a]).abs() / (rho[i] * THETA0));
                measured[1] = measured[1].max((u[i][a] - required_u[i][a]).abs() / THETA0.sqrt());
            }
            let thermal_scale = rho[i] * THETA0 * cv;
            let ke = 0.5 * (u[i][0] * u[i][0] + u[i][1] * u[i][1] + u[i][2] * u[i][2]);
            let caloric = rho[i] * (cv * serial.theta[i] + ke);
            measured[2] = measured[2].max((serial.e[i] - caloric).abs() / thermal_scale);
            let expected = old_e[i] - 0.5 * (div[i] + old_div[i])
                + serial
                    .gravity
                    .unwrap()
                    .work(old_rho[i], old_u[i], rho[i], u[i]);
            measured[3] = measured[3].max((serial.e[i] - expected).abs() / thermal_scale);
        }
        for (actual, recorded) in measured.into_iter().zip(report.residuals()) {
            assert!(
                (actual - recorded).abs() <= 16.0 * f64::EPSILON,
                "published residual {actual:e}, report {recorded:e}"
            );
        }
    }
}

#[test]
fn loose_control_accepts_the_original_three_pass_candidate_without_rebinding() {
    let mut fixed = convecting_gas(1);
    let mut checked = convecting_gas(1)
        .with_coupling_control(control(1.0, 4))
        .unwrap();
    for _ in 0..3 {
        assert_eq!(fixed.step().unwrap(), checked.step().unwrap());
        assert_eq!(fixed.f, checked.f);
        assert_eq!(fixed.e, checked.e);
        assert_eq!(fixed.theta, checked.theta);
        assert_eq!(fixed.force, checked.force);
        assert_eq!(fixed.div, checked.div);
        assert_eq!(fixed.cumulative_ledger(), checked.cumulative_ledger());
        let report = checked.last_coupling_report().unwrap();
        assert_eq!(
            (
                report.corrector_passes,
                report.endpoint_evaluations,
                report.residual_probes
            ),
            (3, 4, 1)
        );
        assert_eq!(fixed.last_coupling_report(), None);
    }
    checked.set_coupling_control(None).unwrap();
    assert_eq!(fixed.state_hash(), checked.state_hash());
    assert_eq!(fixed.step().unwrap(), checked.step().unwrap());
    assert_eq!(fixed.state_hash(), checked.state_hash());
    assert_eq!(checked.last_coupling_report(), None);
}

#[test]
fn exhausted_accuracy_budget_preserves_state_drive_and_last_successful_report() {
    let mut grid = convecting_gas(1)
        .with_coupling_control(control(1.0, 4))
        .unwrap();
    grid.step().unwrap();
    grid.set_coupling_control(Some(control(1e-30, 4))).unwrap();
    let hash = grid.state_hash();
    let before = (
        grid.theta.clone(),
        grid.force.clone(),
        grid.div.clone(),
        grid.steps,
        grid.cumulative_ledger(),
        grid.last_coupling_report(),
    );
    let field = vec![[30.0, -20.0, 10.0]; grid.n];
    let error = grid
        .step_with_body_acceleration_si(&field, || false)
        .unwrap_err();
    let NsfError::CouplingNotConverged { control, report } = error else {
        panic!("expected exhausted endpoint budget, got {error:?}");
    };
    assert!(!control.accepts(report));
    assert!(report.maximum_tolerance_ratio > 1.0);
    assert_eq!(
        (
            report.corrector_passes,
            report.endpoint_evaluations,
            report.residual_probes
        ),
        (3, 4, 1)
    );
    assert_eq!(grid.state_hash(), hash);
    assert_eq!(
        (
            grid.theta.clone(),
            grid.force.clone(),
            grid.div.clone(),
            grid.steps,
            grid.cumulative_ledger(),
            grid.last_coupling_report()
        ),
        before
    );
    assert!(grid.body_acceleration.is_empty());
    let selected = grid.coupling_control();
    let invalid = NsfCouplingControl {
        maximum_endpoint_evaluations: 3,
        ..control
    };
    assert!(grid.set_coupling_control(Some(invalid)).is_err());
    assert_eq!(grid.coupling_control(), selected);
    assert_eq!(grid.last_coupling_report(), before.5);
}

#[test]
fn cancellation_after_endpoint_checks_keeps_the_original_owner() {
    let controls = control(1e-9, 64);
    let mut counted = convecting_gas(1).with_coupling_control(controls).unwrap();
    let mut polls = 0;
    counted
        .step_body_acceleration(None, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(counted.last_coupling_report().is_some());
    // Cancel at the actual last pre-publication checkpoint, after all probes.
    let mut grid = convecting_gas(1).with_coupling_control(controls).unwrap();
    let before = (
        grid.state_hash(),
        grid.theta.clone(),
        grid.force.clone(),
        grid.div.clone(),
        grid.cumulative_ledger(),
    );
    let mut calls = 0;
    assert!(
        grid.step_body_acceleration(None, &mut || {
            calls += 1;
            calls == polls
        })
        .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(
        (
            grid.state_hash(),
            grid.theta.clone(),
            grid.force.clone(),
            grid.div.clone(),
            grid.cumulative_ledger()
        ),
        before
    );
    assert_eq!(grid.steps(), 0);
    assert_eq!(grid.last_coupling_report(), None);
    grid.step().unwrap();
    assert_eq!(grid.state_hash(), counted.state_hash());
    assert_eq!(grid.last_coupling_report(), counted.last_coupling_report());
}

#[test]
fn zero_velocity_heat_source_uses_fixed_reference_norms_and_no_extra_source_samples() {
    let calls = Arc::new(AtomicUsize::new(0));
    let sampled = calls.clone();
    let mut input = spec([4, 3, 2], 1);
    input.periodic = [true; 3];
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        input,
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(move |_, _| {
            sampled.fetch_add(1, Ordering::Relaxed);
            (0.0, [0.0; 3], 0.001)
        })),
    )
    .unwrap()
    .with_coupling_control(control(1e-11, 16))
    .unwrap();
    let initial_calls = calls.load(Ordering::Relaxed);
    let ledger = grid.step().unwrap();
    assert_eq!(calls.load(Ordering::Relaxed) - initial_calls, 2 * grid.n);
    let report = grid.last_coupling_report().unwrap();
    assert!(report.residuals().into_iter().all(f64::is_finite));
    assert!(grid.coupling_control().unwrap().accepts(report));
    assert!(ledger.residual_rel < 2e-12);
    let heat = 0.001 * grid.n as f64 * grid.law.energy_unit_j();
    assert!((ledger.source_energy_j - heat).abs() < 2e-12 * heat);
    assert!(grid.cell([0; 3]).2 > 300.0);
}
