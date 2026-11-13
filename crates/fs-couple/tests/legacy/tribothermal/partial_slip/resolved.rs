//! Finite-body and reversible-contact references for conduction-derived routing.
use super::*;
use fs_couple::tribothermal::partial_slip::resolved::*;

fn solve_controls() -> ResolvedCaloricContactSolve {
    let mut contact = controls();
    contact.newton.linear_restart = 3;
    ResolvedCaloricContactSolve {
        contact,
        temperature_scale_j_per_k: 1.0,
        interface_energy_tolerance_j: 2e-8,
    }
}
fn transport(
    state: &DepositionSimulation<'_>,
    _: f64,
) -> Result<ResolvedCaloricContactSource<()>, Infallible> {
    Ok(ResolvedCaloricContactSource {
        conductivity_w_m_k: [3.0, 1.0],
        transport: other(state)?,
    })
}
fn resolved<'a>(
    thermal: &DepositionSimulation<'_>,
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    model: AffineContactStiffness,
    strain: [f64; 3],
) -> CaloricPartialSlip<'a> {
    let prior = contact(thermal, curves, model.clone(), strain);
    // A square metre face, matching the original thermal interface. The native
    // patch accepts a declared second moment; no ellipse is inferred.
    let checkpoint = PartialSlipCheckpoint::new(
        NormalPatchView::new(
            "patch",
            "held square",
            "synthetic patch",
            NormalPatchAuthority::SyntheticFixture,
            10.0,
            0.5,
            0.5,
            1.0 / 6.0,
        )
        .unwrap(),
        prior.checkpoint().interface().clone(),
        prior.checkpoint().law().clone(),
        prior.checkpoint().state().clone(),
    )
    .unwrap();
    CaloricPartialSlip::new_with_resolved_interface(
        thermal,
        model,
        checkpoint,
        [curves[0].curve(), curves[1].curve()],
        7,
        [3.0, 1.0],
    )
    .unwrap()
}

#[test]
fn physical_slip_heats_a_common_node_and_finite_bodies_follow_their_conductivities() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 1.5);
    let old = resolved(&thermal, &curves, model(0.4, zero_slopes()), [0.0; 3]);
    let step = old
        .step_resolved_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("resolved slide"),
            solve_controls(),
            &[0.0; 2],
            300.0,
            transport,
            || false,
        )
        .unwrap();
    let r = step.report();
    // Independent Jenkins work, then an independent two-capacity/common-node
    // solve. The half-cell conductances are (6,2) W/K, not heat fractions.
    let heat = 0.6 * (0.1 - 0.6 / 750.0) + 3.4 * (0.1 - 3.4 / 250.0);
    let numerical = 0.5 * (0.6f64.powi(2) / 750.0 + 3.4f64.powi(2) / 250.0);
    let g = [6.0, 2.0];
    let capacity = [10.0, 20.0];
    let admittance = [
        g[0] * capacity[0] / (capacity[0] + g[0]),
        g[1] * capacity[1] / (capacity[1] + g[1]),
    ];
    let rise = heat / (admittance[0] + admittance[1]);
    near(r.interface_temperature_k, 300.0 + rise);
    near(r.caloric.contact.dissipation.heat_j, heat);
    near(
        r.caloric.contact.generalized_work.numerical_dissipation_j,
        numerical,
    );
    near(r.caloric.thermal.stored_energy_change_j, heat);
    for i in 0..2 {
        let body_rise = g[i] / (capacity[i] + g[i]) * rise;
        near(step.contact().temperature_k()[i], 300.0 + body_rise);
        near(r.body_heat_j[i], capacity[i] * body_rise);
        near(r.caloric.physical_heat_j[i], heat * g[i] / (g[0] + g[1]));
        near(
            r.conduction_heat_j[i],
            r.body_heat_j[i] - r.caloric.physical_heat_j[i],
        );
    }
    near(r.conduction_heat_j[0] + r.conduction_heat_j[1], 0.0);
    near(r.interface_energy_residual_j, 0.0);
    near(r.caloric.total_energy_residual_j, 0.0);
    assert!(numerical > 0.0);
    assert!(
        r.body_heat_j[0] > 0.6 * heat,
        "old authored 25% split must not survive"
    );
    assert_eq!(old.checkpoint().state().accepted_steps(), 1);
    assert_eq!(step.contact().checkpoint().state().accepted_steps(), 2);
}

#[test]
fn reversible_entropy_exchange_stays_with_each_body_while_the_node_conducts() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 1.5);
    let old = resolved(&thermal, &curves, model(1000.0, slopes()), [0.1, 0.0, 0.0]);
    let step = old
        .step_resolved_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("resolved elasticity"),
            solve_controls(),
            &[0.0; 2],
            300.0,
            transport,
            || false,
        )
        .unwrap();
    let r = step.report();
    let ds = [
        0.5 * (0.2f64.powi(2) - 0.1f64.powi(2)),
        0.25 * (0.2f64.powi(2) - 0.1f64.powi(2)),
    ];
    // C_i*(T_i-300) = 1.5*(T_j-T_i) - T_i*Delta(s_i).
    let a = 10.0 + 1.5 + ds[0];
    let b = 20.0 + 1.5 + ds[1];
    let determinant = a * b - 1.5 * 1.5;
    let expected = [
        (3000.0 * b + 1.5 * 6000.0) / determinant,
        (6000.0 * a + 1.5 * 3000.0) / determinant,
    ];
    for i in 0..2 {
        near(step.contact().temperature_k()[i], expected[i]);
        near(r.caloric.reversible_heat_j[i], -expected[i] * ds[i]);
        near(r.caloric.physical_heat_j[i], 0.0);
        near(r.body_heat_j[i], 1.5 * (expected[1 - i] - expected[i]));
    }
    near(
        r.interface_temperature_k,
        (6.0 * expected[0] + 2.0 * expected[1]) / 8.0,
    );
    near(
        r.caloric.contact_internal_energy_change_j,
        0.5 * 1450.0 * (0.2f64.powi(2) - 0.1f64.powi(2)),
    );
    near(r.caloric.total_energy_residual_j, 0.0);
    assert!(r.caloric.reversible_heat_j.iter().all(|q| *q < 0.0));
    assert!(r.caloric.contact.generalized_work.numerical_dissipation_j > 0.0);
}

#[test]
fn zero_slip_keeps_signed_hot_to_cold_conduction_and_creates_no_heat() {
    let curves = curves();
    let thermal = initial(&curves, [400.0, 300.0], 1.5);
    let old = resolved(&thermal, &curves, model(1000.0, slopes()), [0.1, 0.05, 0.2]);
    let step = old
        .step_resolved_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.0),
            &work("resolved conduction"),
            solve_controls(),
            &[0.0; 2],
            350.0,
            transport,
            || false,
        )
        .unwrap();
    let r = step.report();
    let transfer = 150.0 / (1.0 + 1.5 / 10.0 + 1.5 / 20.0);
    near(step.contact().temperature_k()[0], 400.0 - transfer / 10.0);
    near(step.contact().temperature_k()[1], 300.0 + transfer / 20.0);
    near(r.body_heat_j[0], -transfer);
    near(r.body_heat_j[1], transfer);
    near(r.caloric.contact.dissipation.heat_j, 0.0);
    near(r.caloric.thermal.stored_energy_change_j, 0.0);
    near(r.caloric.contact_internal_energy_change_j, 0.0);
    near(r.caloric.total_energy_residual_j, 0.0);
    assert_eq!(
        step.contact()
            .checkpoint()
            .state()
            .region_elastic_displacement_m(),
        old.checkpoint().state().region_elastic_displacement_m()
    );
}

#[test]
fn routing_admission_refusal_and_retry_preserve_both_original_histories() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 1.5);
    let old = resolved(&thermal, &curves, model(0.4, zero_slopes()), [0.0; 3]);
    let checkpoint = old.checkpoint().clone();
    let cells = thermal.cells().clone();
    let mut budget = solve_controls();
    budget.contact.maximum_contact_trials = 1;
    assert!(
        old.step_resolved_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("budget"),
            budget,
            &[0.0; 2],
            300.0,
            transport,
            || false
        )
        .is_err()
    );
    assert!(matches!(
        old.step_resolved_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("cancel"),
            solve_controls(),
            &[0.0; 2],
            300.0,
            transport,
            || true
        ),
        Err(CaloricContactError::Cancelled)
    ));
    // Reject competing direct-link transport instead of silently doubling heat.
    let competing = old.step_resolved_candidate(
        &thermal,
        1.0,
        frame(),
        motion(0.1),
        &work("competing"),
        solve_controls(),
        &[0.0; 2],
        300.0,
        |state, t| {
            let mut source = transport(state, t)?;
            source.transport.conductances_w_k.insert(7, 1.5);
            Ok::<_, Infallible>(source)
        },
        || false,
    );
    assert!(competing.is_err());
    assert!(
        old.step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("wrong owner"),
            controls(),
            &[0.0; 2],
            other,
            || false
        )
        .is_err()
    );
    assert_eq!(old.checkpoint(), &checkpoint);
    assert_eq!(thermal.cells(), &cells);
    let run = || {
        old.step_resolved_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("retry"),
            solve_controls(),
            &[0.0; 2],
            300.0,
            transport,
            || false,
        )
        .unwrap()
    };
    let a = run();
    let b = run();
    assert_eq!(a.contact().checkpoint(), b.contact().checkpoint());
    assert_eq!(a.thermal().cells(), b.thermal().cells());
    assert_eq!(
        a.report().interface_temperature_k.to_bits(),
        b.report().interface_temperature_k.to_bits()
    );
    let prior = contact(&thermal, &curves, model(0.4, zero_slopes()), [0.0; 3]);
    assert!(
        prior
            .clone()
            .with_resolved_interface(&thermal, 7, [1.0, 1.0])
            .is_err()
    );
    let upgraded = prior
        .clone()
        .with_resolved_interface(&thermal, 7, [3.0, 1.0])
        .unwrap();
    assert_eq!(prior.checkpoint(), upgraded.checkpoint());
}
