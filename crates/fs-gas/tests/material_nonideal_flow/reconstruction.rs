//! Native original-EOS MUSCL momentum transport and complete-step atomicity.
use super::*;

fn muscl() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::Muscl,
        ..controls()
    }
}

fn advance(flow: &mut RealGasEuler1d, finish: f64, control: StepControl) {
    let mut steps = 0;
    while flow.time().value() < finish {
        let report = flow
            .step(Time::new(finish - flow.time().value()), control, || false)
            .unwrap();
        check_balance(report);
        assert!(
            !report.first_order_retry,
            "smooth original-source field lost reconstruction order"
        );
        steps += 1;
        assert!(steps < 2048);
    }
}

fn acoustic_momentum_error(n: usize, control: StepControl) -> f64 {
    let source = gas();
    let reference = source.reference_state();
    let temperature_slope = reference.temperature_k() * reference.root().pressure_derivatives().0
        / (100.0 * reference.specific_heat_cv_j_kg_k());
    let phase_half_width = core::f64::consts::PI / n as f64;
    let average = phase_half_width.sin() / phase_half_width;
    let amplitude = 1e-4;
    let states = (0..n)
        .map(|i| {
            let phase = core::f64::consts::TAU * (i as f64 + 0.5) / n as f64;
            let drho = amplitude * average * phase.sin();
            state(
                &source,
                10.0 + drho,
                DENSE_SOUND_SPEED * drho / 10.0,
                400.0 + temperature_slope * drho,
            )
        })
        .collect();
    let mut flow = uniform_mesh_field(source, states);
    let finish = 0.2 / DENSE_SOUND_SPEED;
    let before = flow.totals(|| false).unwrap();
    advance(&mut flow, finish, control);
    let after = flow.totals(|| false).unwrap();
    near(after.momentum, before.momentum, 1e-10);
    flow.states()
        .iter()
        .enumerate()
        .map(|(i, q)| {
            let phase =
                core::f64::consts::TAU * ((i as f64 + 0.5) / n as f64 - DENSE_SOUND_SPEED * finish);
            let drho = amplitude * average * phase.sin();
            let expected = (10.0 + drho) * DENSE_SOUND_SPEED * drho / 10.0;
            (q.momentum - expected).abs()
        })
        .sum::<f64>()
        / (n as f64 * amplitude * DENSE_SOUND_SPEED)
}

#[test]
fn muscl_resolves_dense_acoustic_momentum_and_refines() {
    let errors = [24, 48, 96].map(|n| acoustic_momentum_error(n, muscl()));
    assert!(
        errors[1] < 0.55 * errors[0],
        "MUSCL momentum errors {errors:?}"
    );
    assert!(
        errors[2] < 0.55 * errors[1],
        "MUSCL momentum errors {errors:?}"
    );
    let first_order = acoustic_momentum_error(48, controls());
    assert!(
        errors[1] < 0.3 * first_order,
        "MUSCL {} must materially reduce momentum diffusion relative to {first_order}",
        errors[1]
    );
}

#[test]
fn muscl_preserves_uniform_dense_flow_on_original_nonuniform_mesh() {
    let source = gas();
    let q = state(&source, 20.0, 65.0, 430.0);
    let positions = vec![0.0, 0.07, 0.20, 0.39, 0.62, 1.0];
    let volumes = positions.windows(2).map(|x| 0.07 * (x[1] - x[0])).collect();
    let mesh = Mesh1d::new(positions, vec![0.07; 6], volumes).unwrap();
    let make = || {
        RealGasEuler1d::new_on_mesh(
            source.clone(),
            mesh.clone(),
            vec![q; 5],
            RealGasBoundaries::Periodic,
            || false,
        )
        .unwrap()
    };
    let mut first = make();
    let mut second = make();
    for _ in 0..4 {
        let a = first.step(Time::new(1e-5), controls(), || false).unwrap();
        let b = second.step(Time::new(1e-5), muscl(), || false).unwrap();
        assert_eq!(a, b);
        assert_eq!(second.states(), &[q; 5]);
        assert_eq!(first.primitives(), second.primitives());
        assert_eq!(b.limited_cells, 0);
    }
}

#[test]
fn muscl_transport_retains_original_formation_energy_gauge() {
    let original = gas();
    let shifted = bind(
        common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]),
        domain(),
    );
    let make = |source: RealGas| {
        let states = (0..12)
            .map(|i| {
                let phase = core::f64::consts::TAU * (i as f64 + 0.5) / 12.0;
                state(
                    &source,
                    20.0 + 0.5 * phase.cos(),
                    35.0 + 4.0 * phase.sin(),
                    450.0 + 2.0 * phase.cos(),
                )
            })
            .collect();
        uniform_mesh_field(source, states)
    };
    let mut a = make(original);
    let mut b = make(shifted);
    for _ in 0..4 {
        let ra = a.step(Time::new(1e-5), muscl(), || false).unwrap();
        let rb = b.step(Time::new(1e-5), muscl(), || false).unwrap();
        check_balance(ra);
        check_balance(rb);
        for (qa, qb) in a.states().iter().zip(b.states()) {
            for (x, y) in [
                (qa.density, qb.density),
                (qa.momentum, qb.momentum),
                (qa.energy, qb.energy),
            ] {
                near(x, y, 2e-11 * (1.0 + x.abs() + y.abs()));
            }
        }
        for (pa, pb) in a.primitives().iter().zip(b.primitives()) {
            near(pa.temperature, pb.temperature, 3e-8);
            near(pa.pressure, pb.pressure, 2e-5);
        }
    }
}

#[test]
fn cancelled_muscl_traces_and_late_publication_preserve_accepted_field() {
    let source = gas();
    let initial = (0..16)
        .map(|i| {
            let phase = core::f64::consts::TAU * (i as f64 + 0.5) / 16.0;
            state(
                &source,
                12.0 + 0.2 * phase.sin(),
                10.0 + phase.cos(),
                430.0 + phase.sin(),
            )
        })
        .collect::<Vec<_>>();
    let mut baseline = uniform_mesh_field(source.clone(), initial.clone());
    let mut polls = 0;
    let expected = baseline
        .step(Time::new(1e-5), muscl(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 32);
    for stop in [1, polls / 4, polls / 2, polls] {
        let mut flow = uniform_mesh_field(source.clone(), initial.clone());
        let original_primitive = flow.primitives().to_vec();
        let mut calls = 0;
        let result = flow.step(Time::new(1e-5), muscl(), || {
            calls += 1;
            calls == stop
        });
        assert_eq!(result, Err(RealGasFlowError::Cancelled));
        assert_eq!(flow.states(), initial);
        assert_eq!(flow.primitives(), original_primitive);
        assert_eq!(flow.time().value(), 0.0);
        let retry = flow.step(Time::new(1e-5), muscl(), || false).unwrap();
        assert_eq!(retry, expected);
        assert_eq!(flow.states(), baseline.states());
        assert_eq!(flow.primitives(), baseline.primitives());
    }
}

#[test]
fn muscl_closed_compression_keeps_pressure_impulse_and_total_energy() {
    let source = gas();
    let states = (0..20)
        .map(|i| state(&source, 20.0, if i < 10 { 200.0 } else { -100.0 }, 430.0))
        .collect();
    let mut flow = RealGasEuler1d::new(
        source,
        Length::new(0.05),
        0.07,
        states,
        RealGasBoundaries::Ends {
            left: RealGasBoundary::Reflecting,
            right: RealGasBoundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let before = flow.totals(|| false).unwrap();
    let mut boundary_impulse = 0.0;
    for _ in 0..8 {
        let report = flow.step(Time::new(2e-6), muscl(), || false).unwrap();
        check_balance(report);
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
        boundary_impulse += report.boundary_change.momentum;
    }
    let after = flow.totals(|| false).unwrap();
    near(after.mass, before.mass, 1e-13);
    near(after.energy, before.energy, 2e-8);
    near(after.momentum - before.momentum, boundary_impulse, 2e-11);
}
