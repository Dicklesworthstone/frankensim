//! Complete open D-C-A-C-D entropy with original stationary bath efforts.
//! Native explicit SSP and RKL2 keep their own conservative boundary quadrature.

use super::*;
use crate::multicomponent::Conserved;
use crate::multicomponent::dissipation::{DissipationError, StabilizedDissipationControl};

const CALORICS: [[f64; 3]; 2] = [[3.5, 0.2, 0.0], [4.0, 0.1, 1.5]];
const WITHOUT_DATUM: [[f64; 3]; 2] = [[3.5, 0.2, 0.0], [4.0, 0.1, 0.0]];
const T0: f64 = 0.001;

fn fixture() -> (ReactingEuler1d, DissipationOperator<ConstantTransport>, [Conserved; 2]) {
    let (simulation, _) = variable_caloric();
    let mechanism = simulation.mechanism().clone();
    let chemistry = simulation.chemistry_controls();
    let mut flow = simulation.into_flow();
    flow.set_boundaries(
        Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        || false,
    ).unwrap();
    // M=R for both original source species, so p=(rho_A+rho_B)*T.
    let baths = [
        flow.gas().conserved(&[0.8, 0.2], 0.0, 1.3).unwrap(),
        flow.gas().conserved(&[0.3, 0.7], 0.0, 0.9).unwrap(),
    ];
    let operator = molecular(&flow, 0.5, 0.02, 0.4);
    (
        ReactingEuler1d::new(flow, mechanism, chemistry, || false).unwrap(),
        operator, baths,
    )
}

fn budget_fixture(hot: bool) -> (ReactingEuler1d, DissipationOperator<ConstantTransport>, Conserved) {
    let mechanism = if hot { mechanism(40.0, 500.0, 0.0) } else { mechanism(0.0, 2.0, 0.0) };
    let rho = if hot { [0.99, 0.01] } else { [0.5, 0.5] };
    let mut flow = field(&mechanism, 4, 0.0, |_| rho);
    flow.set_boundaries(
        Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        || false,
    ).unwrap();
    let bath = flow.gas().conserved(&rho, 0.0, 1.0).unwrap();
    let operator = molecular(&flow, 0.05, 0.01, 0.1);
    (
        ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap(),
        operator, bath,
    )
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-10 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}",
    );
}

// At stationary exterior states, -v_b.dQ + psi_b*A*dt reduces to
// E_ref/T_b - sum[(e_i,ref/T_b+1-s_i,b)*dm_i], since u_b=psi_b=0.
// The SAME effort contracts actual inviscid and molecular transfers, including
// partial enthalpy and the bath's viscous energy share. No Fourier-only heat
// surrogate or production entropy routine enters this independent oracle.
fn supply(
    rho: [f64; 2],
    temperature: f64,
    calorics: &[[f64; 3]; 2],
    transfer: &Totals,
) -> f64 {
    let mut input = transfer.energy / temperature;
    for (k, &[a, b, c]) in calorics.iter().enumerate() {
        let internal = (a - 1.0) * (temperature - T0)
            + 0.5 * b * (temperature * temperature - T0 * T0);
        let entropy = (a - 1.0) * temperature.ln() + b * temperature + c - rho[k].ln();
        input -= (internal / temperature + 1.0 - entropy) * transfer.species_masses[k];
    }
    input
}

fn inputs(step: &ReactiveStepReport, calorics: &[[f64; 3]; 2]) -> [f64; 2] {
    let halves = step.dissipation.unwrap();
    [
        [
            step.transport.left_boundary_change,
            halves.first.left_boundary_change,
            halves.second.left_boundary_change,
        ].iter().map(|q| supply([0.8, 0.2], 1.3, calorics, q)).sum(),
        [
            step.transport.right_boundary_change,
            halves.first.right_boundary_change,
            halves.second.right_boundary_change,
        ].iter().map(|q| supply([0.3, 0.7], 0.9, calorics, q)).sum(),
    ]
}

#[test]
fn g1_open_ssp_and_rkl_compositions_include_both_molecular_baths_and_real_reaction_entropy() {
    for stabilized in [false, true] {
        let (mut simulation, mut operator, baths) = fixture();
        let before = source_entropy(simulation.flow(), &CALORICS);
        let stripped_before = source_entropy(simulation.flow(), &WITHOUT_DATUM);
        let report = if stabilized {
            let result = simulation.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
                &mut operator, Time::new(0.002), controls(),
                StabilizedDissipationControl::default(), MixtureEntropyControl::default(),
                [Some(baths[0].as_cell()), Some(baths[1].as_cell())], || false,
            ).unwrap();
            assert_eq!(result.coupled.coupled.coupling_retries, 0);
            assert_eq!(result.operator_evaluations,
                result.dissipation.first.operator_evaluations + result.dissipation.second.operator_evaluations);
            assert_eq!(result.dissipation.first.transport.elapsed.value(), 0.001);
            assert_eq!(result.dissipation.second.transport.elapsed.value(), 0.001);
            result.coupled
        } else {
            simulation.step_entropy_consistent_with_dissipation_and_exterior_states(
                &mut operator, Time::new(0.002), controls(),
                DissipationControl::default(), MixtureEntropyControl::default(),
                [Some(baths[0].as_cell()), Some(baths[1].as_cell())], || false,
            ).unwrap()
        };
        let step = report.coupled;
        let halves = step.dissipation.unwrap();
        let measured_halves = report.dissipation_entropy.unwrap();
        assert_eq!(step.elapsed.value(), 0.002);
        assert_eq!(simulation.flow().time(), step.elapsed);
        assert_eq!(halves.first.elapsed.value(), 0.5 * step.elapsed.value());
        assert_eq!(halves.second.elapsed.value(), 0.5 * step.elapsed.value());
        assert_eq!(step.chemistry.cell_solves, 2 * simulation.flow().cell_count());
        assert_eq!(step.transport.hll_fallback_faces, 0);
        assert!(step.chemistry.species_change_kg[1] > 1e-4);
        for k in 0..2 {
            near(
                step.after.species_masses[k] - step.before.species_masses[k],
                step.transport.boundary_change.species_masses[k]
                    + halves.first.boundary_change.species_masses[k]
                    + halves.second.boundary_change.species_masses[k]
                    + step.chemistry.species_change_kg[k],
            );
        }
        near(step.absolute_energy_residual_j, 0.0);
        near(step.mass_residual_kg, 0.0);
        near(step.conservation_residual.momentum, 0.0);
        assert!(report.entropy.maximum_relative_species_residual < 1e-12);

        let boundary = inputs(&step, &CALORICS);
        let measured = source_entropy(simulation.flow(), &CALORICS) - before;
        near(report.entropy.source_entropy_change_j_k, measured);
        near(report.entropy.source_entropy_generation_j_k, measured - boundary.iter().sum::<f64>());
        near(report.entropy.entropy_generation_j_k,
            report.entropy.source_entropy_generation_j_k - report.entropy.reference_entropy_residual_j_k);
        for end in 0..2 {
            near(report.entropy.boundary_entropy_input_j_k[end], boundary[end]);
            near(boundary[end], report.transport_entropy.boundary_entropy_input_j_k[end]
                + measured_halves.first.boundary_entropy_input_j_k[end]
                + measured_halves.second.boundary_entropy_input_j_k[end]);
        }
        for (physical, entropy) in [
            (halves.first, measured_halves.first), (halves.second, measured_halves.second),
        ] {
            near(entropy.boundary_entropy_input_j_k[0],
                supply([0.8, 0.2], 1.3, &CALORICS, &physical.left_boundary_change));
            near(entropy.boundary_entropy_input_j_k[1],
                supply([0.3, 0.7], 0.9, &CALORICS, &physical.right_boundary_change));
            assert!(entropy.entropy_generation_j_k >= -1e-10);
        }
        assert!(report.transport_entropy.entropy_generation_j_k >= -1e-10);
        assert!(report.entropy.entropy_generation_j_k > 0.0);
        assert!((measured - step.chemistry.entropy_change_j_k
            - report.transport_entropy.source_entropy_change_j_k
            - measured_halves.first.source_entropy_change_j_k
            - measured_halves.second.source_entropy_change_j_k).abs() < 2e-8);

        // Stripping B's a9 only in the independent oracle removes genuine
        // chemical entropy; its total molecular+advective boundary contribution
        // cancels against the physical species transfer, not chemical creation.
        let stripped = source_entropy(simulation.flow(), &WITHOUT_DATUM) - stripped_before
            - inputs(&step, &WITHOUT_DATUM).iter().sum::<f64>();
        let reaction_datum = 1.5 * step.chemistry.species_change_kg[1];
        assert!(reaction_datum > 1e-4);
        near(report.entropy.source_entropy_generation_j_k - stripped, reaction_datum);
    }
}

#[test]
fn g4_open_rkl_one_work_budget_covers_both_halves_and_abandoned_coupled_attempts() {
    for hot in [false, true] {
        let (mut simulation, mut operator, bath) = budget_fixture(hot);
        let initial = simulation.accepted.states.clone();
        let duration = Time::new(if hot { 0.2 } else { 0.002 });
        let error = simulation.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
            &mut operator, duration, controls(),
            StabilizedDissipationControl {
                maximum_stages: 2, maximum_operator_evaluations: 3,
                ..StabilizedDissipationControl::default()
            },
            MixtureEntropyControl::default(),
            [Some(bath.as_cell()), Some(bath.as_cell())], || false,
        ).unwrap_err();
        // The first spatially uniform D half consumes exactly three calls.
        // Cold: the second D half must share the spent counter.
        // Hot: acoustic shortening discards D-C-A, but that work remains spent.
        assert!(matches!(error, ReactiveError::Dissipation(
            DissipationError::OperatorWorkLimit { operator_evaluations: 3, .. }
        )));
        assert_eq!(simulation.accepted.states, initial);
        assert_eq!(simulation.flow().time().value(), 0.0);
        let recovered = simulation.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
            &mut operator, duration, controls(),
            StabilizedDissipationControl::default(), MixtureEntropyControl::default(),
            [Some(bath.as_cell()), Some(bath.as_cell())], || false,
        ).unwrap();
        let accepted_work = recovered.dissipation.first.operator_evaluations
            + recovered.dissipation.second.operator_evaluations;
        if hot {
            assert!(recovered.coupled.coupled.coupling_retries > 0);
            assert!(recovered.operator_evaluations > accepted_work);
        } else {
            assert_eq!(recovered.coupled.coupled.coupling_retries, 0);
            assert_eq!(recovered.operator_evaluations, 6);
            assert_eq!(accepted_work, 6);
        }
        assert_eq!(simulation.flow().time(), recovered.coupled.coupled.elapsed);
        near(recovered.coupled.coupled.absolute_energy_residual_j, 0.0);
        assert!(recovered.coupled.entropy.entropy_generation_j_k >= -1e-10);
    }
}

#[test]
fn g4_open_rkl_late_cancellation_restores_the_whole_composition_and_reuses_both_owners() {
    let (mut reference, mut reference_operator, baths) = fixture();
    let mut polls = 0;
    let expected = reference.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
        &mut reference_operator, Time::new(0.002), controls(),
        StabilizedDissipationControl::default(), MixtureEntropyControl::default(),
        [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
        || { polls += 1; false },
    ).unwrap();
    assert!(polls > 100);
    let (mut simulation, mut operator, baths) = fixture();
    let initial = simulation.accepted.states.clone();
    let mut visits = 0;
    let error = simulation.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
        &mut operator, Time::new(0.002), controls(),
        StabilizedDissipationControl::default(), MixtureEntropyControl::default(),
        [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
        || { visits += 1; visits >= polls },
    ).unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(simulation.accepted.states, initial);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let actual = simulation.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
        &mut operator, Time::new(0.002), controls(),
        StabilizedDissipationControl::default(), MixtureEntropyControl::default(),
        [Some(baths[0].as_cell()), Some(baths[1].as_cell())], || false,
    ).unwrap();
    assert_eq!(actual.coupled.coupled.elapsed, expected.coupled.coupled.elapsed);
    assert_eq!(actual.operator_evaluations, expected.operator_evaluations);
    for (a, b) in simulation.accepted.states.iter().zip(&reference.accepted.states) { near(*a, *b); }
    near(actual.coupled.entropy.entropy_generation_j_k, expected.coupled.entropy.entropy_generation_j_k);
    for end in 0..2 {
        near(actual.coupled.entropy.boundary_entropy_input_j_k[end],
            expected.coupled.entropy.boundary_entropy_input_j_k[end]);
    }
}

#[test]
fn g4_moving_molecular_bath_refuses_before_any_coupled_retry_or_publication() {
    for stabilized in [false, true] {
        let (mut simulation, mut operator, baths) = fixture();
        let initial = simulation.accepted.states.clone();
        let moving = simulation.flow().gas().conserved(&[0.8, 0.2], 0.1, 1.3).unwrap();
        let error = if stabilized {
            simulation.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
                &mut operator, Time::new(0.002), controls(),
                StabilizedDissipationControl::default(), MixtureEntropyControl::default(),
                [Some(moving.as_cell()), Some(baths[1].as_cell())], || false,
            ).unwrap_err()
        } else {
            simulation.step_entropy_consistent_with_dissipation_and_exterior_states(
                &mut operator, Time::new(0.002), controls(),
                DissipationControl::default(), MixtureEntropyControl::default(),
                [Some(moving.as_cell()), Some(baths[1].as_cell())], || false,
            ).unwrap_err()
        };
        assert!(!matches!(error, ReactiveError::RetryLimit { .. }));
        assert_eq!(simulation.accepted.states, initial);
        assert_eq!(simulation.flow().time().value(), 0.0);
    }
}

#[test]
fn g0_absent_baths_keep_closed_ssp_and_rkl_compositions_bit_for_bit() {
    for stabilized in [false, true] {
        let (mut closed, mut closed_operator) = variable_caloric();
        let (mut through_open_api, mut open_operator) = variable_caloric();
        let (expected, actual) = if stabilized {
            let expected = closed.step_entropy_consistent_with_stabilized_dissipation(
                &mut closed_operator, Time::new(0.002), controls(),
                StabilizedDissipationControl::default(), MixtureEntropyControl::default(), || false,
            ).unwrap();
            let actual = through_open_api.step_entropy_consistent_with_stabilized_dissipation_and_exterior_states(
                &mut open_operator, Time::new(0.002), controls(),
                StabilizedDissipationControl::default(), MixtureEntropyControl::default(), [None; 2], || false,
            ).unwrap();
            assert_eq!(actual.dissipation, expected.dissipation);
            assert_eq!(actual.operator_evaluations, expected.operator_evaluations);
            (expected.coupled, actual.coupled)
        } else {
            let expected = closed.step_entropy_consistent_with_dissipation(
                &mut closed_operator, Time::new(0.002), controls(),
                DissipationControl::default(), MixtureEntropyControl::default(), || false,
            ).unwrap();
            let actual = through_open_api.step_entropy_consistent_with_dissipation_and_exterior_states(
                &mut open_operator, Time::new(0.002), controls(),
                DissipationControl::default(), MixtureEntropyControl::default(), [None; 2], || false,
            ).unwrap();
            (expected, actual)
        };
        assert_eq!(actual.coupled, expected.coupled);
        assert_eq!(through_open_api.accepted.states, closed.accepted.states);
        assert_eq!(through_open_api.flow().time(), closed.flow().time());
        assert_eq!(actual.entropy.boundary_entropy_input_j_k, [0.0; 2]);
        assert_eq!(actual.entropy.entropy_generation_j_k, expected.entropy.entropy_change_j_k);
        assert_eq!(actual.entropy.source_entropy_change_j_k, expected.entropy.source_entropy_change_j_k);
        let open = actual.dissipation_entropy.unwrap();
        let closed = expected.dissipation_entropy.unwrap();
        for (a, e) in [(open.first, closed.first), (open.second, closed.second)] {
            assert_eq!(a.boundary_entropy_input_j_k, [0.0; 2]);
            assert_eq!(a.entropy_generation_j_k, e.entropy_change_j_k);
            assert_eq!(a.source_entropy_change_j_k, e.source_entropy_change_j_k);
            assert_eq!(a.minimum_step_entropy_generation_j_k, e.minimum_step_entropy_change_j_k);
        }
    }
}
