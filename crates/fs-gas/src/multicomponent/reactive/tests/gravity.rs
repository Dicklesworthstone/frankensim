//! G0/G1/G3/G4: chemical evolution in a conservative gravitational potential.

use super::*;
use crate::dynamics::IsothermalAtmosphere;
use crate::multicomponent::{IsothermalMixtureGravity, Mesh1d, Reconstruction};
use fs_qty::{Acceleration, Pressure};

fn atmosphere(mechanism: &IdealGasKineticsV1, acceleration: f64) -> IsothermalMixtureGravity {
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let positions: Vec<_> = (0..=8).map(|i| (i as f64 / 8.0).powf(1.3)).collect();
    let volumes = positions.windows(2).map(|x| 0.7 * (x[1] - x[0])).collect();
    IsothermalMixtureGravity::try_new(
        gas,
        Mesh1d::new(positions, vec![0.7; 9], volumes).unwrap(),
        IsothermalAtmosphere {
            acceleration: Acceleration::new(acceleration),
            temperature: Temperature::new(2.0),
            reference_pressure: Pressure::new(3.0),
            reference_position: Length::new(0.0),
        },
        &[0.8, 0.2],
        || false,
    )
    .unwrap()
}

fn simulation(
    reference: &IsothermalMixtureGravity,
    mechanism: &IdealGasKineticsV1,
) -> ReactingEuler1d {
    ReactingEuler1d::new(
        reference.initial_flow(|| false).unwrap(),
        mechanism.clone(),
        chemistry(),
        || false,
    )
    .unwrap()
}

fn controls() -> ReactiveStepControl {
    ReactiveStepControl {
        transport: StepControl {
            reconstruction: Reconstruction::FirstOrder,
            ..StepControl::default()
        },
        ..ReactiveStepControl::default()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn g1_isomer_relaxation_retains_hydrostatic_rest_and_actual_gravity_mass() {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    for g in [-1.5, 1.5] {
        let reference = atmosphere(&mechanism, g);
        let mut flow = simulation(&reference, &mechanism);
        let mass = flow.flow().totals().species_masses[..2].iter().sum::<f64>();
        for _ in 0..10 {
            let step = flow
                .step_with_isothermal_gravity(Time::new(0.001), controls(), &reference, || false)
                .unwrap();
            near(
                step.physical_gravity_impulse_kg_m_s,
                g * mass * step.reactive.elapsed.value(),
                2e-12,
            );
            near(
                step.reactive.transport.source_change.impulse,
                step.physical_gravity_impulse_kg_m_s
                    + step.reference_traction_impulse_defect_kg_m_s,
                2e-12,
            );
            near(step.total_energy_residual_j, 0.0, 3e-10);
            near(step.reactive.mass_residual_kg, 0.0, 3e-12);
            near(step.potential_energy_change_j, 0.0, 3e-11);
            assert!(step.reactive.chemistry.species_change_kg[1] > 0.0);
        }
        let expected_fraction = 0.5 + 0.3 * (-4.0 * flow.flow().time().value()).exp();
        for i in 0..flow.flow().cell_count() {
            let q = flow.flow().cell(i).unwrap();
            let primitive = flow.flow().gas().primitive(q, || false).unwrap();
            let (_, original) = reference.reference_cell(i).unwrap();
            near(primitive.density, original.density, 3e-10);
            near(primitive.temperature, 2.0, 3e-9);
            near(primitive.velocity, 0.0, 3e-9);
            near(
                q.partial_densities[0] / primitive.density,
                expected_fraction,
                2e-6,
            );
        }
    }
}

#[test]
fn g1_heat_release_replays_native_gravity_and_closes_chemical_plus_potential_energy() {
    let mechanism = mechanism(12.0, 80.0, 0.0);
    let reference = atmosphere(&mechanism, -3.0);
    let mut flow = simulation(&reference, &mechanism);
    let mut manual = simulation(&reference, &mechanism);
    let step = flow
        .step_with_isothermal_gravity(Time::new(0.1), controls(), &reference, || false)
        .unwrap();
    let h = step.reactive.elapsed.value();
    assert!(h < 0.1);
    assert_eq!(flow.flow().time(), step.reactive.elapsed);

    // Independently sequence the original chemistry and gravity owners at the
    // finally accepted duration, rather than reusing the coupled entry point.
    let mut chemistry_report = ChemistryReport::default();
    react_cells(
        &mut manual.accepted,
        &mut manual.reactor,
        &mut manual.amounts,
        0.5 * h,
        &mut chemistry_report,
        &mut || false,
    )
    .unwrap();
    let transport = manual
        .accepted
        .step_with_isothermal_gravity(Time::new(h), controls().transport, &reference, || false)
        .unwrap();
    assert_eq!(transport.transport.elapsed.value(), h);
    react_cells(
        &mut manual.accepted,
        &mut manual.reactor,
        &mut manual.amounts,
        0.5 * h,
        &mut chemistry_report,
        &mut || false,
    )
    .unwrap();
    for (actual, expected) in flow.accepted.states.iter().zip(&manual.accepted.states) {
        near(*actual, *expected, 3e-10);
    }
    assert!(step.reactive.chemistry.referenced_energy_change_j > 0.0);
    assert!(
        flow.accepted
            .states
            .chunks_exact(4)
            .any(|q| q[2].abs() > 1e-7)
    );
    near(step.total_energy_residual_j, 0.0, 3e-9);
    near(step.reactive.absolute_energy_residual_j, 0.0, 3e-9);
    near(
        step.potential_energy_change_j + step.reactive.transport.source_change.work,
        0.0,
        3e-9,
    );
    near(
        step.potential_energy_after_j,
        reference.potential_energy(flow.flow(), || false).unwrap(),
        1e-14,
    );
}

#[test]
fn g3_formation_reference_shift_and_zero_gravity_retain_physical_motion() {
    let original = mechanism(8.0, 5.0, 0.0);
    let shifted = mechanism(8.0, 5.0, 1e4);
    let a = atmosphere(&original, -2.0);
    let b = atmosphere(&shifted, -2.0);
    let mut left = simulation(&a, &original);
    let mut right = simulation(&b, &shifted);
    let h = Time::new(0.001);
    let ra = left
        .step_with_isothermal_gravity(h, controls(), &a, || false)
        .unwrap();
    let rb = right
        .step_with_isothermal_gravity(h, controls(), &b, || false)
        .unwrap();
    assert_eq!(ra.reactive.elapsed, rb.reactive.elapsed);
    for i in 0..left.flow().cell_count() {
        let pa = left
            .flow()
            .gas()
            .primitive(left.flow().cell(i).unwrap(), || false)
            .unwrap();
        let pb = right
            .flow()
            .gas()
            .primitive(right.flow().cell(i).unwrap(), || false)
            .unwrap();
        near(pa.velocity, pb.velocity, 3e-8);
        near(pa.temperature, pb.temperature, 3e-7);
        near(pa.density, pb.density, 3e-8);
    }
    near(ra.total_energy_residual_j, 0.0, 3e-9);
    near(rb.total_energy_residual_j, 0.0, 3e-8);

    let zero = atmosphere(&original, 0.0);
    let mut ordinary = simulation(&zero, &original);
    let mut gravity = simulation(&zero, &original);
    let expected = ordinary.step(h, controls(), || false).unwrap();
    let actual = gravity
        .step_with_isothermal_gravity(h, controls(), &zero, || false)
        .unwrap();
    assert_eq!(actual.reactive, expected);
    assert_eq!(gravity.accepted.states, ordinary.accepted.states);
    assert_eq!(actual.potential_energy_change_j, 0.0);
    assert_eq!(actual.physical_gravity_impulse_kg_m_s, 0.0);
}

#[test]
fn g4_late_gravity_account_cancellation_and_wrong_source_are_atomic() {
    let mechanism = mechanism(8.0, 5.0, 0.0);
    let reference = atmosphere(&mechanism, -2.0);
    let mut expected = simulation(&reference, &mechanism);
    let mut polls = 0;
    let h = Time::new(0.001);
    let accepted = expected
        .step_with_isothermal_gravity(h, controls(), &reference, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut flow = simulation(&reference, &mechanism);
    let original = flow.accepted.states.clone();
    let mut attempted = 0;
    assert_eq!(
        flow.step_with_isothermal_gravity(h, controls(), &reference, || {
            attempted += 1;
            attempted >= polls - 1
        })
        .unwrap_err(),
        ReactiveError::Cancelled,
    );
    assert_eq!(flow.accepted.states, original);
    assert_eq!(flow.flow().time().value(), 0.0);
    let foreign = atmosphere(&super::mechanism(8.0, 5.0, 1.0), -2.0);
    assert!(
        flow.step_with_isothermal_gravity(h, controls(), &foreign, || false)
            .is_err()
    );
    assert_eq!(flow.accepted.states, original);
    let replay = flow
        .step_with_isothermal_gravity(h, controls(), &reference, || false)
        .unwrap();
    assert_eq!(replay.reactive.elapsed, accepted.reactive.elapsed);
    for (actual, expected) in flow.accepted.states.iter().zip(&expected.accepted.states) {
        near(*actual, *expected, 3e-10);
    }
}
