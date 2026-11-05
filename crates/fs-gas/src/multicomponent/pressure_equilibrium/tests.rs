use super::*;
use crate::multicomponent::tests::{gas as mixture, source};

fn gas(variable: bool, references: [f64; 2]) -> MixtureGas {
    mixture(vec![
        source(
            "A",
            R,
            3.5,
            if variable { 0.25 } else { 0.0 },
            references[0],
        ),
        source(
            "B",
            2.0 * R,
            5.0,
            if variable { 0.7 } else { 0.0 },
            references[1],
        ),
        source("absent-C", 0.5 * R, 4.5, 0.2, 7.0),
    ])
}
fn control(budget: f64) -> ContactEquilibriumControl {
    ContactEquilibriumControl {
        relative_tolerance: 1e-8,
        maximum_energy_defect_j: budget,
    }
}
fn workspace(cells: usize, budget: f64, velocity: [f64; 3]) -> ContactWorkspace {
    let bytes = ContactWorkspace::required_workspace_bytes(cells).unwrap();
    let mut work = ContactWorkspace::new(cells, bytes, || false).unwrap();
    work.admit(cells, control(budget)).unwrap();
    assert_eq!(work.workspace_bytes(), bytes);
    work.begin_attempt(control(budget)).unwrap();
    for i in 0..cells {
        work.set_reference(i, 1.0, velocity, 1.25_f64.sqrt(), true)
            .unwrap();
    }
    work
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}
fn mixed_energy(gas: &MixtureGas, velocity: &[f64]) -> f64 {
    let a = gas.conserved(&[1.0, 0.0, 0.0], velocity[0], 1.0).unwrap();
    let b = gas.conserved(&[0.0, 1.0, 0.0], velocity[0], 1.0).unwrap();
    let mut energy = 0.25 * a.energy + 0.75 * b.energy;
    for v in &velocity[1..] {
        energy += 0.5 * v * v;
    }
    energy
}
fn pressure(gas: &MixtureGas, densities: &[f64], velocity: &[f64], energy: f64) -> f64 {
    let mut reduced = energy;
    for v in &velocity[1..] {
        reduced -= 0.5 * v * v;
    }
    gas.primitive(
        Cell {
            partial_densities: densities,
            momentum: velocity[0],
            energy: reduced,
        },
        || false,
    )
    .unwrap()
    .pressure
}

#[test]
fn g0_exact_contact_mixing_and_final_rk_average_retain_the_explicit_energy_defect() {
    let gas = gas(false, [0.0, 0.0]);
    let velocity = [0.4, 0.0, 0.0];
    let densities = [0.25, 0.75, 0.0];
    let original = mixed_energy(&gas, &velocity[..1]);
    close(
        pressure(&gas, &densities, &velocity[..1], original),
        145.0 / 136.0,
    );
    let mut work = workspace(1, 1.0, velocity);
    let first = work
        .project_energy(
            &gas,
            0,
            &densities,
            &velocity[..1],
            original,
            1.0,
            0.5,
            || false,
        )
        .unwrap();
    close(pressure(&gas, &densities, &velocity[..1], first), 1.0);
    close(first - original, -9.0 / 40.0);
    close(work.receipt().signed_energy_j, -9.0 / 80.0);
    close(work.receipt().absolute_energy_j, 9.0 / 80.0);

    // A donor A cell advects into B with Courant number1/4. The second FE
    // stage and final SSP average need their own original-EOS projections.
    let donor = gas
        .conserved(&[1.0, 0.0, 0.0], velocity[0], 1.0)
        .unwrap()
        .energy;
    let initial = gas
        .conserved(&[0.0, 1.0, 0.0], velocity[0], 1.0)
        .unwrap()
        .energy;
    work.set_reference(0, 1.0, velocity, 1.2, false).unwrap();
    work.observe_pair(0, 1.0, velocity, 1.0, velocity, false)
        .unwrap();
    let second_raw = 0.75 * first + 0.25 * donor;
    let second = work
        .project_energy(
            &gas,
            0,
            &[0.4375, 0.5625, 0.0],
            &velocity[..1],
            second_raw,
            1.0,
            0.5,
            || false,
        )
        .unwrap();
    let averaged = 0.5 * initial + 0.5 * second;
    let final_densities = [0.21875, 0.78125, 0.0];
    assert!(pressure(&gas, &final_densities, &velocity[..1], averaged) > 1.01);
    let final_energy = work
        .project_energy(
            &gas,
            0,
            &final_densities,
            &velocity[..1],
            averaged,
            1.0,
            1.0,
            || false,
        )
        .unwrap();
    close(
        pressure(&gas, &final_densities, &velocity[..1], final_energy),
        1.0,
    );
    let boundary_change = 0.125 * ((donor - initial) + (donor - first));
    close(
        final_energy - initial - boundary_change,
        work.receipt().signed_energy_j,
    );
    close(work.receipt().signed_energy_j, -31.0 / 130.0);
    close(
        work.receipt().absolute_energy_j,
        0.5 * (first - original).abs()
            + 0.5 * (second - second_raw).abs()
            + (final_energy - averaged).abs(),
    );
    assert_eq!(work.receipt().corrected_cells, 1);
    let receipt = work.receipt();
    assert_eq!(
        work.project_energy(
            &gas,
            0,
            &final_densities,
            &velocity[..1],
            final_energy,
            1.0,
            1.0,
            || false,
        )
        .unwrap(),
        final_energy
    );
    assert_eq!(work.receipt(), receipt);
}

#[test]
fn g3_original_variable_calorics_preserve_species_datums_and_all_kinetic_components() {
    let densities = [0.25, 0.75, 0.0];
    for variable in [false, true] {
        for velocity in [[0.4, 0.0, 0.0], [3.4, -1.5, 2.0], [-2.6, 1.5, -2.0]] {
            for dimensions in 1..=3 {
                let mut admitted_velocity = velocity;
                admitted_velocity[dimensions..].fill(0.0);
                let mut baseline = None;
                for references in [[0.0, 0.0], [1e18, -1e18], [-1e18, 1e18]] {
                    let gas = gas(variable, references);
                    let momenta = &admitted_velocity[..dimensions];
                    let original = mixed_energy(&gas, momenta);
                    let mut work = workspace(1, 2.0, admitted_velocity);
                    let projected = work
                        .project_energy(&gas, 0, &densities, momenta, original, 0.7, 1.0, || false)
                        .unwrap();
                    let slope = if variable { 0.325 } else { 0.0 };
                    let thermal =
                        2.125 * (1.6 - 0.001) + 0.5 * slope * (1.6_f64.powi(2) - 0.001_f64.powi(2));
                    let kinetic: f64 = momenta.iter().map(|v| 0.5 * v * v).sum();
                    close(projected, thermal + kinetic);
                    close(pressure(&gas, &densities, momenta, projected), 1.0);
                    close(
                        work.receipt().signed_energy_j,
                        0.7 * if variable { -0.36525 } else { -0.225 },
                    );
                    assert_eq!(work.receipt().corrected_cells, 1);
                    if let Some(expected) = baseline {
                        assert_eq!((projected, work.receipt()), expected);
                    } else {
                        baseline = Some((projected, work.receipt()));
                    }
                }
            }
        }
    }
}

#[test]
fn g4_recipient_reference_hll_and_acoustic_velocity_guards_remain_excluded() {
    let gas = gas(false, [0.0, 0.0]);
    let velocity = [100.0, -2.0, 0.0];
    let mut work = workspace(5, 10.0, velocity);
    // Equal left/right face pressures do not make a nonconstant cell stencil
    // an equilibrium contact: both must match the recipient cell itself.
    work.observe_pair(0, 1.01, velocity, 1.01, velocity, false)
        .unwrap();
    work.observe_pair(1, 1.0, velocity, 1.0, velocity, true)
        .unwrap();
    let mut accelerated = velocity;
    accelerated[0] += 1e-7;
    // This difference exceeds tolerance*c but is smaller than tolerance*|u|.
    work.observe_pair(2, 1.0, velocity, 1.0, accelerated, false)
        .unwrap();
    work.set_reference(3, 1.01, velocity, 1.0, false).unwrap();
    work.exclude(4).unwrap();
    let energy = mixed_energy(&gas, &velocity[..2]);
    for i in 0..5 {
        work.set_reference(i, 1.0, velocity, 1.25_f64.sqrt(), false)
            .unwrap();
        work.observe_pair(i, 1.0, velocity, 1.0, velocity, false)
            .unwrap();
        assert_eq!(
            work.project_energy(
                &gas,
                i,
                &[0.25, 0.75, 0.0],
                &velocity[..2],
                energy,
                1.0,
                1.0,
                || false,
            )
            .unwrap(),
            energy
        );
    }
    assert_eq!(work.receipt(), ContactCorrection::default());
    work.begin_attempt(control(10.0)).unwrap();
    work.set_reference(0, 1.0, velocity, 1.25_f64.sqrt(), true)
        .unwrap();
    let accelerated_energy = mixed_energy(&gas, &accelerated[..2]);
    assert_eq!(
        work.project_energy(
            &gas,
            0,
            &[0.25, 0.75, 0.0],
            &accelerated[..2],
            accelerated_energy,
            1.0,
            1.0,
            || false,
        )
        .unwrap(),
        accelerated_energy
    );
    assert_eq!(work.receipt(), ContactCorrection::default());
}

#[test]
fn g4_energy_budget_source_window_precision_and_cancellation_are_atomic() {
    let gas = gas(false, [0.0, 0.0]);
    let velocity = [0.4, 0.0, 0.0];
    let densities = [0.25, 0.75, 0.0];
    let original = mixed_energy(&gas, &velocity[..1]);
    let target = gas.conserved(&densities, velocity[0], 1.0).unwrap().energy;
    let mut work = workspace(2, 0.3, velocity);
    work.project_energy(
        &gas,
        0,
        &densities,
        &velocity[..1],
        original,
        1.0,
        1.0,
        || false,
    )
    .unwrap();
    let previous = work.receipt();
    assert!(matches!(work.project_energy(
        &gas, 1, &densities, &velocity[..1], target - 0.225, 1.0, 1.0, || false,
    ), Err(Error::ContactEnergyBudget { observed_j, maximum_j })
        if observed_j > 0.44 && maximum_j == 0.3));
    assert_eq!(work.receipt(), previous);

    let mut domain = gas.domain();
    domain.maximum_temperature = 1.5;
    let narrow = MixtureGas::new(gas.source().clone(), domain, || false).unwrap();
    let mut work = workspace(1, 1.0, velocity);
    assert!(matches!(
        work.project_energy(
            &narrow,
            0,
            &densities,
            &velocity[..1],
            original,
            1.0,
            1.0,
            || false,
        ),
        Err(Error::EnergyDomain)
    ));
    assert_eq!(work.receipt(), ContactCorrection::default());
    work.project_energy(
        &gas,
        0,
        &densities,
        &velocity[..1],
        original,
        1.0,
        1.0,
        || false,
    )
    .unwrap();
    assert_eq!(work.receipt().corrected_cells, 1);

    let mut work = workspace(1, 1.0, [1e10, 0.0, 0.0]);
    assert!(matches!(
        work.project_energy(&gas, 0, &densities, &[1e10], 5e19, 1.0, 1.0, || false,),
        Err(Error::Invalid(_))
    ));
    assert_eq!(work.receipt(), ContactCorrection::default());

    let mut probe = workspace(1, 1.0, velocity);
    let mut polls = 0;
    let expected = probe
        .project_energy(
            &gas,
            0,
            &densities,
            &velocity[..1],
            original,
            1.0,
            1.0,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in 1..=polls {
        let mut work = workspace(1, 1.0, velocity);
        let mut count = 0;
        assert!(matches!(
            work.project_energy(
                &gas,
                0,
                &densities,
                &velocity[..1],
                original,
                1.0,
                1.0,
                || {
                    count += 1;
                    count >= cutoff
                },
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(work.receipt(), ContactCorrection::default());
        assert_eq!(
            work.project_energy(
                &gas,
                0,
                &densities,
                &velocity[..1],
                original,
                1.0,
                1.0,
                || false,
            )
            .unwrap(),
            expected
        );
        assert_eq!(work.receipt(), probe.receipt());
    }

    let bytes = ContactWorkspace::required_workspace_bytes(1).unwrap();
    assert!(ContactWorkspace::new(1, bytes - 1, || false).is_err());
    assert!(ContactWorkspace::required_workspace_bytes(usize::MAX).is_err());
    assert!(matches!(
        ContactWorkspace::new(1, bytes, || true),
        Err(Error::Cancelled)
    ));
    assert!(probe.admit(2, control(1.0)).is_err());
    let previous = probe.receipt();
    assert!(
        probe
            .begin_attempt(ContactEquilibriumControl {
                relative_tolerance: 1e-5,
                maximum_energy_defect_j: 1.0,
            })
            .is_err()
    );
    assert_eq!(probe.receipt(), previous);
}
