use super::*;

fn controls(step: f64) -> StabilizedDissipationControl {
    StabilizedDissipationControl {
        transport: DissipationControl {
            maximum_step_s: step,
            ..DissipationControl::default()
        },
        ..StabilizedDissipationControl::default()
    }
}

fn heat(count: usize) -> MixtureEuler1d {
    field(gas(1, 0.0, false), count, Boundaries::Periodic, |x| {
        (
            vec![1.0],
            0.0,
            2.0 + 0.2 * (core::f64::consts::TAU * x).sin(),
        )
    })
}

// Independent eigenmode oracle: evaluate the Legendre stability polynomial,
// not the implementation's state/flux Runge--Kutta recurrence.
fn amplification(stages: usize, z: f64) -> f64 {
    let s = stages as f64;
    let x = 1.0 + 4.0 * z / (s * s + s - 2.0);
    let (mut older, mut previous) = (1.0, x);
    for j in 2..=stages {
        let j = j as f64;
        let next = ((2.0 * j - 1.0) * x * previous - (j - 1.0) * older) / j;
        older = previous;
        previous = next;
    }
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    1.0 - b + b * previous
}

#[test]
fn g1_rkl2_resolves_the_extended_heat_spectrum_with_fewer_rhs_calls() {
    let n = 128;
    let setup = || {
        field(gas(1, 0.0, false), n, Boundaries::Periodic, |x| {
            let checker = if (x * n as f64) as usize % 2 == 0 {
                1.0
            } else {
                -1.0
            };
            (
                vec![1.0],
                0.0,
                2.0 + 0.2 * (core::f64::consts::TAU * x).sin() + 0.01 * checker,
            )
        })
    };
    let mut flow = setup();
    let mut ordinary = setup();
    let mut op = operator(&flow, 0.0, 0.0, 0.05, 0.0);
    let mut explicit = operator(&ordinary, 0.0, 0.0, 0.05, 0.0);
    let dt = 0.05;
    let report = op
        .advance_stabilized_for(&mut flow, Time::new(dt), controls(dt), || false)
        .unwrap();
    let explicit = explicit
        .advance_for(
            &mut ordinary,
            Time::new(dt),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 2);
    assert_eq!(report.operator_evaluations, report.maximum_stages_used + 1);
    assert!(report.operator_evaluations * 2 < 2 * explicit.accepted_steps);
    assert!(report.transport.maximum_parabolic_number > 1.0);
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
    let scale = 4.0 * 0.02 * (n * n) as f64;
    let low = amplification(
        report.maximum_stages_used,
        -dt * scale * (core::f64::consts::PI / n as f64).sin().powi(2),
    );
    let high = amplification(report.maximum_stages_used, -dt * scale);
    assert!(low.abs() <= 1.0 && high.abs() <= 1.0);
    for i in 0..n {
        let x = (i as f64 + 0.5) / n as f64;
        let checker = if i % 2 == 0 { 1.0 } else { -1.0 };
        let expected = 2.0 + 0.2 * low * (core::f64::consts::TAU * x).sin() + 0.01 * high * checker;
        let p = flow
            .gas()
            .primitive(flow.cell(i).unwrap(), || false)
            .unwrap();
        assert!(
            (p.temperature - expected).abs() < 3e-11,
            "cell {i}: {p:?}, expected {expected}"
        );
    }
    assert_eq!(flow.time().value(), 0.0);
    assert!(report.transport.absolute_energy_residual_j.abs() < 3e-12);
    eprintln!(
        "RKL2 heat: {} RHS calls versus {} SSP RK2 RHS calls",
        report.operator_evaluations,
        2 * explicit.accepted_steps
    );
}

fn temporal_error(step: f64) -> f64 {
    let n = 96;
    let mut flow = heat(n);
    let mut op = operator(&flow, 0.0, 0.0, 0.05, 0.0);
    op.advance_stabilized_for(&mut flow, Time::new(0.16), controls(step), || false)
        .unwrap();
    let eigenvalue =
        -4.0 * 0.02 * (n * n) as f64 * (core::f64::consts::PI / n as f64).sin().powi(2);
    let amplitude = 0.2 * (eigenvalue * 0.16).exp();
    (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) / n as f64;
            let actual = flow
                .gas()
                .primitive(flow.cell(i).unwrap(), || false)
                .unwrap()
                .temperature;
            (actual - (2.0 + amplitude * (core::f64::consts::TAU * x).sin())).abs() / n as f64
        })
        .sum()
}

#[test]
fn g1_rkl2_heat_is_second_order_in_physical_macrostep() {
    let coarse = temporal_error(0.04);
    let medium = temporal_error(0.02);
    let fine = temporal_error(0.01);
    eprintln!("RKL2 temporal L1 errors: {coarse:e}, {medium:e}, {fine:e}");
    assert!(coarse < 1e-5);
    assert!(medium < 0.31 * coarse);
    assert!(fine < 0.31 * medium);
}

#[test]
fn g1_rkl2_binary_maxwell_stefan_mode_keeps_mass_frame_and_temperature() {
    let n = 96;
    let mut flow = field(gas(2, 0.0, false), n, Boundaries::Periodic, |x| {
        let a = 0.5 + 0.2 * (core::f64::consts::TAU * x).sin();
        (vec![a, 1.0 - a], 0.0, 2.0)
    });
    let mut op = operator(&flow, 0.0, 0.0, 0.0, 0.03);
    let report = op
        .advance_stabilized_for(&mut flow, Time::new(0.1), controls(0.01), || false)
        .unwrap();
    let eigenvalue =
        -4.0 * 0.03 * (n * n) as f64 * (core::f64::consts::PI / n as f64).sin().powi(2);
    let amplitude = 0.2 * (eigenvalue * 0.1).exp();
    let mut error = 0.0;
    for i in 0..n {
        let x = (i as f64 + 0.5) / n as f64;
        let q = flow.cell(i).unwrap();
        let p = flow.gas().primitive(q, || false).unwrap();
        error += (q.partial_densities[0] - (0.5 + amplitude * (core::f64::consts::TAU * x).sin()))
            .abs()
            / n as f64;
        assert!((p.density - 1.0).abs() < 3e-12);
        assert!((p.temperature - 2.0).abs() < 5e-11);
    }
    assert!(error < 1e-6, "binary temporal error: {error:e}");
    assert!(report.maximum_stages_used > 2);
    assert!(
        report.transport.conservation_residual.species_masses[..2]
            .iter()
            .all(|v| v.abs() < 2e-13)
    );
    assert!(report.transport.maximum_mass_frame_residual_kg_m2_s < 3e-14);
}

#[test]
fn g0_rkl2_unequal_species_preserve_enthalpy_cross_effects_and_energy_references() {
    use fs_thermochem::MassSpecificEnergyQuantityV1;
    let initial = |x: f64| {
        let a = 0.5 + 0.2 * (core::f64::consts::TAU * x).sin();
        (
            vec![a, 2.0 * (1.0 - a)],
            0.0,
            2.0 + 0.1 * (core::f64::consts::TAU * x).cos(),
        )
    };
    let mut ordinary = field(gas(2, 0.0, true), 64, Boundaries::Periodic, initial);
    let mut shifted = field(gas(2, 1e18, true), 64, Boundaries::Periodic, initial);
    let mut explicit = field(gas(2, 0.0, true), 64, Boundaries::Periodic, initial);
    let make = |flow: &MixtureEuler1d| {
        let original = law(flow, 0.0, 0.0, 0.02, 0.03);
        let domain = original.domain();
        let heats = flow
            .gas()
            .source()
            .models()
            .iter()
            .enumerate()
            .map(|(i, model)| {
                (
                    model.species().clone(),
                    MassSpecificEnergyQuantityV1::new(if i == 0 { 0.4 } else { -0.2 }),
                )
            })
            .collect();
        let coupled =
            SoretDufourTransport::new(original, "manufactured-cross-transport", 1, domain, heats)
                .unwrap();
        DissipationOperator::new(flow, coupled, || false).unwrap()
    };
    let mut first = make(&ordinary);
    let mut second = make(&shifted);
    let mut reference = make(&explicit);
    let report = first
        .advance_stabilized_for(&mut ordinary, Time::new(0.03), controls(0.003), || false)
        .unwrap();
    second
        .advance_stabilized_for(&mut shifted, Time::new(0.03), controls(0.003), || false)
        .unwrap();
    reference
        .advance_for(
            &mut explicit,
            Time::new(0.03),
            DissipationControl {
                maximum_step_s: 0.0001,
                ..DissipationControl::default()
            },
            || false,
        )
        .unwrap();
    assert_eq!(ordinary.states, shifted.states);
    assert!(
        ordinary
            .states
            .iter()
            .zip(&explicit.states)
            .all(|(a, b)| (a - b).abs() < 2e-5)
    );
    assert!(report.maximum_stages_used > 2);
    assert!(
        report.transport.conservation_residual.species_masses[..2]
            .iter()
            .all(|v| v.abs() < 2e-13)
    );
    assert!(report.transport.absolute_energy_residual_j.abs() < 3e-12);
    assert!(report.transport.maximum_mass_frame_residual_kg_m2_s < 3e-14);
}

#[test]
fn g0_rkl2_stationary_wall_traction_uses_the_stage_boundary_ledger() {
    let mut flow = field(
        gas(1, 0.0, false),
        40,
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        |_| (vec![1.0], 0.3, 2.0),
    );
    let initial_kinetic = kinetic(&flow);
    let mut op = operator(&flow, 0.03, 0.0, 0.01, 0.0);
    let report = op
        .advance_stabilized_for(&mut flow, Time::new(0.01), controls(0.002), || false)
        .unwrap()
        .transport;
    assert!(report.left_boundary_change.momentum < 0.0);
    assert!(report.right_boundary_change.momentum < 0.0);
    assert_eq!(report.boundary_change.energy, 0.0);
    assert_eq!(report.boundary_change.mass(), 0.0);
    assert!(
        (report.boundary_change.momentum
            - report.left_boundary_change.momentum
            - report.right_boundary_change.momentum)
            .abs()
            < 2e-15
    );
    assert!(report.conservation_residual.momentum.abs() < 3e-13);
    assert!(report.absolute_energy_residual_j.abs() < 3e-12);
    assert!(kinetic(&flow) < initial_kinetic);
}

#[test]
fn g4_rkl2_late_cancellation_and_both_work_budgets_are_atomic_and_reusable() {
    let mut reference = heat(48);
    let mut reference_op = operator(&reference, 0.0, 0.0, 0.05, 0.0);
    let mut polls = 0;
    let control = controls(0.01);
    reference_op
        .advance_stabilized_for(&mut reference, Time::new(0.03), control, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut flow = heat(48);
    let original = flow.states.clone();
    let mut op = operator(&flow, 0.0, 0.0, 0.05, 0.0);
    let mut calls = 0;
    let error = op
        .advance_stabilized_for(&mut flow, Time::new(0.03), control, || {
            calls += 1;
            calls >= polls - 2
        })
        .unwrap_err();
    assert_eq!(error, DissipationError::Cancelled);
    assert_eq!(flow.states, original);
    let mut limited = control;
    limited.transport.maximum_steps = 1;
    let error = op
        .advance_stabilized_for(&mut flow, Time::new(0.03), limited, || false)
        .unwrap_err();
    assert!(matches!(
        error,
        DissipationError::WorkLimit {
            accepted_steps: 1,
            ..
        }
    ));
    assert_eq!(flow.states, original);
    limited = control;
    limited.maximum_operator_evaluations = 2;
    let error = op
        .advance_stabilized_for(&mut flow, Time::new(0.03), limited, || false)
        .unwrap_err();
    assert!(matches!(
        error,
        DissipationError::OperatorWorkLimit {
            operator_evaluations: 2,
            ..
        }
    ));
    assert_eq!(flow.states, original);
    assert_eq!(flow.time().value(), 0.0);
    op.advance_stabilized_for(&mut flow, Time::new(0.03), control, || false)
        .unwrap();
    assert!(
        flow.states
            .iter()
            .zip(&reference.states)
            .all(|(a, b)| (a - b).abs() < 2e-12)
    );
    let accepted = flow.states.clone();
    let noop = op
        .advance_stabilized_for(&mut flow, Time::new(0.0), control, || false)
        .unwrap();
    assert_eq!(noop.operator_evaluations, 0);
    assert_eq!(noop.maximum_stages_used, 0);
    assert_eq!(flow.states, accepted);
    for invalid in [
        StabilizedDissipationControl {
            maximum_stages: 1,
            ..control
        },
        StabilizedDissipationControl {
            maximum_stages: 257,
            ..control
        },
        StabilizedDissipationControl {
            maximum_operator_evaluations: 0,
            ..control
        },
    ] {
        assert!(matches!(
            op.advance_stabilized_for(&mut flow, Time::new(0.03), invalid, || false),
            Err(DissipationError::Invalid(_))
        ));
        assert_eq!(flow.states, accepted);
    }
}

#[test]
fn g0_rkl2_nonuniform_heat_uses_original_physical_face_distances_and_volumes() {
    let g = gas(1, 0.0, false);
    let widths = [0.2, 0.3, 0.5];
    let temperature = [2.1, 2.35, 2.75];
    let states = temperature
        .iter()
        .map(|t| g.conserved(&[1.0], 0.0, *t).unwrap())
        .collect();
    let mesh = Mesh1d::new(
        vec![0.0, 0.2, 0.5, 1.0],
        vec![0.7; 4],
        widths.iter().map(|w| 0.7 * w).collect(),
    )
    .unwrap();
    let mut flow = MixtureEuler1d::new_on_mesh(
        g,
        mesh,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    let dt = 1e-3;
    let alpha = 0.05 / 2.5;
    let first = [alpha / widths[0], 0.0, -alpha / widths[2]];
    let left = alpha * (first[1] - first[0]) / 0.25;
    let right = alpha * (first[2] - first[1]) / 0.4;
    let second = [
        left / widths[0],
        (right - left) / widths[1],
        -right / widths[2],
    ];
    let mut op = operator(&flow, 0.0, 0.0, 0.05, 0.0);
    let report = op
        .advance_stabilized_for(&mut flow, Time::new(dt), controls(dt), || false)
        .unwrap();
    assert_eq!(report.maximum_stages_used, 2);
    for i in 0..3 {
        let actual = flow
            .gas()
            .primitive(flow.cell(i).unwrap(), || false)
            .unwrap()
            .temperature;
        let expected = temperature[i] + dt * first[i] + 0.5 * dt * dt * second[i];
        assert!((actual - expected).abs() < 2e-12);
    }
    assert!(report.transport.absolute_energy_residual_j.abs() < 2e-13);
}

#[test]
fn g4_rkl2_rechecks_actual_state_coefficients_and_halves_growing_parabolic_trials() {
    struct HeatingSensitive(ConstantTransport);
    impl TransportLaw for HeatingSensitive {
        fn source(&self) -> &fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1 {
            self.0.source()
        }
        fn evaluate(
            &mut self,
            temperature: Temperature,
            pressure: fs_qty::Pressure,
            mole_fractions: &[f64],
            binary_diffusivities: &mut [f64],
            cancelled: &mut dyn FnMut() -> bool,
        ) -> Result<TransportCoefficients, DissipationError> {
            let mut coefficients = self.0.evaluate(
                temperature,
                pressure,
                mole_fractions,
                binary_diffusivities,
                cancelled,
            )?;
            // An authored manufactured law: wall dissipation raises the actual
            // temperature and thus conductivity after the initial rate bound.
            coefficients.thermal_conductivity = ThermalConductivityQuantityV1::new(
                0.001 * (1.0 + 100_000.0 * (temperature.value() - 2.0).max(0.0)),
            );
            Ok(coefficients)
        }
    }
    let mut flow = field(
        gas(1, 0.0, false),
        12,
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        |_| (vec![1.0], 0.5, 2.0),
    );
    let sensitive = HeatingSensitive(law(&flow, 0.03, 0.0, 0.001, 0.0));
    let mut op = DissipationOperator::new(&flow, sensitive, || false).unwrap();
    let report = op
        .advance_stabilized_for(&mut flow, Time::new(0.02), controls(0.02), || false)
        .unwrap();
    assert!(report.transport.rejected_steps > 0);
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert!(report.transport.absolute_energy_residual_j.abs() < 3e-12);
    assert!(report.transport.conservation_residual.momentum.abs() < 3e-13);
    assert_eq!(report.transport.elapsed.value(), 0.02);
    assert_eq!(flow.time().value(), 0.0);
}
