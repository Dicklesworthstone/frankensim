//! Open RKL2 thermodynamics from actual native boundary transfers and stages.

use super::super::{EntropyTracker, ReservoirEntropyAccount};
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationError, DissipationOperator,
    MolecularEntropyControl, StabilizedDissipationControl, TransportCoefficients, TransportDomain,
};
use crate::multicomponent::tests::{gas, source};
use crate::multicomponent::{Boundaries, Boundary, Conserved, MixtureEuler1d, MixtureGas};
use fs_qty::{DynViscosity, Length, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn controls(step: f64) -> StabilizedDissipationControl {
    StabilizedDissipationControl {
        transport: DissipationControl {
            maximum_step_s: step,
            ..DissipationControl::default()
        },
        ..StabilizedDissipationControl::default()
    }
}

// Every fixture source has M=R and Cv=2.5; density is one.
fn field(
    gas: MixtureGas,
    count: usize,
    initial: impl Fn(usize) -> (Vec<f64>, f64),
) -> MixtureEuler1d {
    let states = (0..count)
        .map(|i| {
            let (rho, t) = initial(i);
            gas.conserved(&rho, 0.0, t).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        0.7,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap()
}

fn law(flow: &MixtureEuler1d, k: f64, d: f64) -> ConstantTransport {
    let n = flow.gas().species_count();
    let mut binary = vec![d; n * n];
    for i in 0..n {
        binary[i * n + i] = 0.0;
    }
    ConstantTransport::new(
        flow.gas().source().clone(),
        "open-rkl2-entropy-fixture",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.0),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(k),
        },
        binary,
    )
    .unwrap()
}

fn temperature(flow: &MixtureEuler1d, cell: usize) -> f64 {
    flow.gas()
        .primitive(flow.cell(cell).unwrap(), || false)
        .unwrap()
        .temperature
}

fn binary(formation: f64) -> MixtureEuler1d {
    field(
        gas(vec![
            source("A", R, 3.5, 0.0, formation),
            source("B", R, 3.5, 0.0, -2.0 * formation),
        ]),
        2,
        |_| (vec![0.5, 0.5], 2.0),
    )
}

fn composition_baths(flow: &MixtureEuler1d) -> [Conserved; 2] {
    [
        flow.gas().conserved(&[0.8, 0.2], 0.0, 2.0).unwrap(),
        flow.gas().conserved(&[0.2, 0.8], 0.0, 2.0).unwrap(),
    ]
}

fn mixing(x: f64) -> f64 {
    -x * x.ln() - (1.0 - x) * (1.0 - x).ln()
}

// Independent Legendre polynomial, rather than the native state/flux recurrence.
fn amplification(stages: usize, z: f64) -> f64 {
    let s = stages as f64;
    let x = 1.0 + 4.0 * z / (s * s + s - 2.0);
    let (mut previous, mut value) = (1.0, x);
    for j in 2..=stages {
        let j = j as f64;
        let next = ((2.0 * j - 1.0) * x * value - (j - 1.0) * previous) / j;
        previous = value;
        value = next;
    }
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    1.0 - b + b * value
}

#[test]
fn g1_multistage_dirichlet_heat_keeps_actual_bath_entropy_and_transient_storage() {
    let n = 64;
    let dt = 0.05;
    let mut flow = field(gas(vec![source("A", R, 3.5, 0.0, 0.0)]), n, |i| {
        let x = (i as f64 + 0.5) / n as f64;
        (vec![1.0], 2.0 + x + 0.2 * (core::f64::consts::PI * x).sin())
    });
    let baths = [
        flow.gas().conserved(&[1.0], 0.0, 2.0).unwrap(),
        flow.gas().conserved(&[1.0], 0.0, 3.0).unwrap(),
    ];
    let initial: Vec<_> = (0..n).map(|i| temperature(&flow, i)).collect();
    let mut op = DissipationOperator::new(&flow, law(&flow, 0.05, 0.0), || false).unwrap();
    let (report, account) = op
        .advance_stabilized_entropy_consistent_with_reservoirs_for(
            &mut flow,
            Time::new(dt),
            controls(dt),
            MolecularEntropyControl::default(),
            [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
            || false,
        )
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 2);
    assert_eq!(report.operator_evaluations, report.maximum_stages_used + 1);
    let eigenvalue =
        -4.0 * 0.02 * (n * n) as f64 * (core::f64::consts::PI / (2 * n) as f64).sin().powi(2);
    let amplitude = amplification(report.maximum_stages_used, dt * eigenvalue);
    let mut energy = 0.0;
    let mut state_entropy = 0.0;
    for (i, &old) in initial.iter().enumerate() {
        let x = (i as f64 + 0.5) / n as f64;
        let expected = 2.0 + x + 0.2 * amplitude * (core::f64::consts::PI * x).sin();
        assert!((temperature(&flow, i) - expected).abs() < 4e-11);
        energy += 0.7 * 2.5 / n as f64 * (expected - old);
        state_entropy += 0.7 * 2.5 / n as f64 * (expected / old).ln();
    }
    // The linear steady profile contributes opposite constant heat transfers.
    // The symmetric Dirichlet mode loses half of its represented energy at
    // each end, including the signed native multistage quadrature.
    let ends = [
        -0.7 * 0.05 * dt + 0.5 * energy,
        0.7 * 0.05 * dt + 0.5 * energy,
    ];
    let bath_entropy = [ends[0] / 2.0, ends[1] / 3.0];
    assert!((report.transport.left_boundary_change.energy - ends[0]).abs() < 3e-11);
    assert!((report.transport.right_boundary_change.energy - ends[1]).abs() < 3e-11);
    for (actual, expected) in account.boundary_entropy_input_j_k.iter().zip(bath_entropy) {
        assert!((actual - expected).abs() < 2e-11);
    }
    assert!(state_entropy < 0.0);
    assert!((account.source_entropy_change_j_k - state_entropy).abs() < 2e-11);
    let generation = state_entropy - bath_entropy.iter().sum::<f64>();
    assert!(generation > 0.0);
    assert!((account.entropy_generation_j_k - generation).abs() < 3e-11);
    assert!(report.transport.absolute_energy_residual_j.abs() < 3e-11);
    assert_eq!(flow.time().value(), 0.0);
}

// Original equal-M binary Maxwell--Stefan face at uniform temperature/density.
fn species_flux(left: f64, right: f64, distance: f64) -> f64 {
    let x = 0.5 * (left + right);
    -0.03 * x * (1.0 - x) * ((right / (1.0 - right)).ln() - (left / (1.0 - left)).ln()) / distance
}

#[test]
fn g0_signed_two_stage_baths_preserve_source_entropy_and_absolute_species_throughput() {
    let dt = 0.4;
    let first = species_flux(0.8, 0.5, 0.25);
    let middle = 0.5 + (dt / 3.0) * first / 0.5;
    let second = species_flux(0.8, middle, 0.25);
    let interior = species_flux(middle, 1.0 - middle, 0.5);
    // Independent two-stage Butcher weights are [-1/2,3/2]. The endpoint
    // evaluation admits the final source state but is not a third quadrature term.
    let final_a = 0.5 + (dt / 0.5) * (-0.5 * first + 1.5 * (second - interior));
    let mass = 0.7 * dt * (-0.5 * first + 1.5 * second);
    let throughput = 2.0 * 0.7 * dt * (0.5 * first.abs() + 1.5 * second.abs());
    let state_entropy = 0.7 * (mixing(final_a) - mixing(0.5));
    let bath_entropy = -mass * 4.0_f64.ln();
    let mut reference = None;
    for formation in [0.0, 1e18] {
        let mut flow = binary(formation);
        let baths = composition_baths(&flow);
        let mut op = DissipationOperator::new(&flow, law(&flow, 0.0, 0.03), || false).unwrap();

        // Use the real native evolution with its ordinary account initialization
        // so the accepted absolute throughput can be inspected directly. No stage
        // or flux is injected, and both finite-state admissions still run.
        let admitted = crate::multicomponent::dissipation::reservoir::admit(
            flow.gas(),
            flow.boundaries,
            [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
            &mut || false,
        )
        .unwrap();
        let mut tracker =
            EntropyTracker::new(&flow, MolecularEntropyControl::default(), &mut || false).unwrap();
        tracker.boundary = Some(
            ReservoirEntropyAccount::new(
                flow.gas(),
                tracker.reference_temperature,
                &tracker.reference_density,
                &tracker.reference_entropy,
                admitted,
                &mut || false,
            )
            .unwrap(),
        );
        let report = crate::multicomponent::dissipation::stabilized::advance_entropy_consistent(
            &mut op,
            &mut flow,
            Time::new(dt),
            controls(dt),
            admitted,
            &mut tracker,
            &mut || false,
        )
        .unwrap();
        let account = tracker.completed_reservoir_report.unwrap();
        assert_eq!(report.maximum_stages_used, 2);
        assert_eq!(report.operator_evaluations, 3);
        assert_eq!(report.transport.rejected_steps, 0);
        assert!((flow.cell(0).unwrap().partial_densities[0] - final_a).abs() < 2e-12);
        assert!((flow.cell(1).unwrap().partial_densities[0] - (1.0 - final_a)).abs() < 2e-12);
        assert!((report.transport.left_boundary_change.species_masses[0] - mass).abs() < 2e-13);
        assert!((report.transport.right_boundary_change.species_masses[0] + mass).abs() < 2e-13);
        let boundary = tracker.boundary.as_ref().unwrap();
        for k in 0..2 {
            assert!((boundary.total.throughput[k].to_f64() - throughput).abs() < 3e-13);
            assert!(boundary.total.species[k].to_f64().abs() < 3e-13);
        }
        // This would fail if opposite ends or the negative anchor were canceled
        // before taking absolute stage contributions.
        assert!(throughput > 2.0 * (2.0 * mass.abs()));
        for actual in account.boundary_entropy_input_j_k {
            assert!((actual - bath_entropy).abs() < 3e-12);
        }
        assert!(state_entropy < 0.0);
        assert!((account.source_entropy_change_j_k - state_entropy).abs() < 3e-12);
        assert!(
            (account.entropy_generation_j_k - (state_entropy - 2.0 * bath_entropy)).abs() < 3e-12
        );
        assert!(account.entropy_generation_j_k > 0.0);
        for i in 0..2 {
            assert!((temperature(&flow, i) - 2.0).abs() < 2e-12);
        }
        if let Some((states, previous)) = &reference {
            assert_eq!(&flow.states, states);
            assert_eq!(&account, previous);
        } else {
            assert!(report.transport.absolute_energy_residual_j.abs() < 3e-12);
            reference = Some((flow.states.clone(), account));
        }
    }
}

#[test]
fn g4_open_cancellation_and_incompatible_bath_leave_no_partial_transfer_before_reuse() {
    let control = controls(0.2);
    let mut reference = binary(0.0);
    let baths = composition_baths(&reference);
    let mut op =
        DissipationOperator::new(&reference, law(&reference, 0.0, 0.03), || false).unwrap();
    let mut polls = 0;
    let (_, expected) = op
        .advance_stabilized_entropy_consistent_with_reservoirs_for(
            &mut reference,
            Time::new(1.2),
            control,
            MolecularEntropyControl::default(),
            [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 2);

    let mut flow = binary(0.0);
    let original = flow.states.clone();
    let mut op = DissipationOperator::new(&flow, law(&flow, 0.0, 0.03), || false).unwrap();
    let mut count = 0;
    let error = op
        .advance_stabilized_entropy_consistent_with_reservoirs_for(
            &mut flow,
            Time::new(1.2),
            control,
            MolecularEntropyControl::default(),
            [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
            || {
                count += 1;
                count >= polls - 2
            },
        )
        .unwrap_err();
    assert_eq!(error, DissipationError::Cancelled);
    assert_eq!(flow.states, original);
    assert_eq!(flow.time().value(), 0.0);

    let absent = flow.gas().conserved(&[1.0, 0.0], 0.0, 2.0).unwrap();
    let error = op
        .advance_stabilized_entropy_consistent_with_reservoirs_for(
            &mut flow,
            Time::new(1.2),
            control,
            MolecularEntropyControl::default(),
            [Some(absent.as_cell()), Some(baths[1].as_cell())],
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, DissipationError::Invalid(_)));
    assert_eq!(flow.states, original);

    let (_, actual) = op
        .advance_stabilized_entropy_consistent_with_reservoirs_for(
            &mut flow,
            Time::new(1.2),
            control,
            MolecularEntropyControl::default(),
            [Some(baths[0].as_cell()), Some(baths[1].as_cell())],
            || false,
        )
        .unwrap();
    assert!(
        flow.states
            .iter()
            .zip(&reference.states)
            .all(|(a, b)| (a - b).abs() < 2e-12)
    );
    assert!((actual.entropy_generation_j_k - expected.entropy_generation_j_k).abs() < 3e-12);
    for (actual, expected) in actual
        .boundary_entropy_input_j_k
        .iter()
        .zip(expected.boundary_entropy_input_j_k)
    {
        assert!((actual - expected).abs() < 3e-12);
    }
}
