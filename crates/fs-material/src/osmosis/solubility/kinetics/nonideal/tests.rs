//! Pure-pair calibration, multisolute phase evolution and common chemical work.

use super::super::KineticSolubilityError;
use super::*;
use fs_qty::VolumetricFlowRate;
use fs_thermochem::{
    SpeciesId,
    solution::solubility::{
        IdealSolubilityV1,
        kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1},
    },
};

const T: f64 = 300.0;
const V: f64 = 1e-3;
const VBAR: f64 = 18e-6;

fn ideal(
    dissolved: &[f64],
    solid: &[f64],
    saturation: &[f64],
    particles: bool,
) -> KineticSolubilitySolution {
    let components = dissolved
        .iter()
        .zip(solid)
        .zip(saturation)
        .enumerate()
        .map(|(i, ((&aq, &crystal), &sat))| {
            let law = IdealSolubilityV1::new(
                SpeciesId::new("water").unwrap(),
                SpeciesId::new(format!("component-{i}")).unwrap(),
                Temperature::new(T),
                AmountConcentration::new(sat),
                0.0,
                [Temperature::new(280.0), Temperature::new(340.0)],
            )
            .unwrap();
            let model = SolubilityKineticsV1::new(
                law,
                VolumetricFlowRate::new(V * 0.1 / (i + 1) as f64),
                PrecipitationPolicyV1::ExistingSolidOnly,
            )
            .unwrap();
            let model = if particles {
                model
                    .with_particle_surface(Amount::new(crystal * V))
                    .unwrap()
            } else {
                model
            };
            (
                model,
                SolubilityKineticStateV1::from_phases(
                    Amount::new(aq * V),
                    Amount::new(crystal * V),
                )
                .unwrap(),
            )
        })
        .collect();
    KineticSolubilitySolution::new(Temperature::new(T), VBAR, Volume::new(V), components).unwrap()
}

fn source(
    reference: f64,
    cap: f64,
    interaction: Vec<f64>,
    count: usize,
) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, VBAR, reference, vec![cap; count], interaction).unwrap())
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn original_binary_saturation_and_exact_particle_trajectory_survive_activity_binding() {
    for diagonal in [-0.002, 0.01] {
        let mut original = ideal(&[10.0], &[200.0], &[50.0], true);
        let mut bound = original
            .clone()
            .with_nonideal_thermodynamics(source(1.0, 300.0, vec![diagonal], 1))
            .unwrap();
        assert_eq!(
            bound.conditional_saturation(0, &[200.0]).unwrap().value(),
            50.0
        );
        let amounts = bound.solution().solute_moles().as_ptr();
        let phases = bound.last_steps().as_ptr();
        original.advance(Time::new(0.7), || false).unwrap();
        bound.advance(Time::new(0.7), || false).unwrap();
        assert_eq!(bound.last_steps(), original.last_steps());
        assert_eq!(bound.models(), original.models());
        assert_eq!(bound.solution().solute_moles().as_ptr(), amounts);
        assert_eq!(bound.last_steps().as_ptr(), phases);
        let concentration = bound.solution().solute_moles()[0] / V;
        close(
            bound.solution().osmotic_pressure_pa().unwrap(),
            R * T * concentration * (1.0 + 0.5 * diagonal * concentration),
            2e-15,
        );
    }
}

#[test]
fn cross_activities_shift_conditional_saturation_in_both_directions() {
    for cross in [-0.001, 0.001] {
        let bound = ideal(&[10.0, 100.0], &[200.0, 100.0], &[50.0, 80.0], false)
            .with_nonideal_thermodynamics(source(1.0, 300.0, vec![0.0, cross, cross, 0.01], 2))
            .unwrap();
        let saturation = bound
            .conditional_saturation(0, &[10.0, 100.0])
            .unwrap()
            .value();
        close(saturation, 50.0 * (-cross * 100.0).exp(), 2e-15);
        assert_eq!(saturation < 50.0, cross > 0.0);
        let mut chemical = [None; 2];
        bound
            .solid_relative_chemical_potentials_into(&[saturation, 100.0], &mut chemical)
            .unwrap();
        assert!(chemical[0].unwrap().abs() < 2e-12);
        bound
            .solid_relative_chemical_potentials_into(&[0.0, 100.0], &mut chemical)
            .unwrap();
        assert_eq!(chemical[0], None);
    }
    let bound = ideal(&[10.0, 80.0], &[200.0, 100.0], &[50.0, 80.0], false)
        .with_nonideal_thermodynamics(source(1.0, 300.0, vec![0.01, 0.005, 0.005, 0.02], 2))
        .unwrap();
    let mut bracket = [0.0, 50.0];
    for _ in 0..100 {
        let c: f64 = 0.5 * (bracket[0] + bracket[1]);
        if (c / 50.0).ln() + 0.01 * (c - 50.0) + 0.005 * 80.0 > 0.0 {
            bracket[1] = c;
        } else {
            bracket[0] = c;
        }
    }
    close(
        bound
            .conditional_saturation(0, &[10.0, 80.0])
            .unwrap()
            .value(),
        0.5 * (bracket[0] + bracket[1]),
        8e-15,
    );
}

// Independent continuous ODE. With zero self interactions the conditional
// concentration is csat_i*exp(-B_ij*c_j), so this oracle needs no material root.
fn continuous(duration: f64, particles: bool) -> [f64; 4] {
    let rhs = |x: [f64; 4]| {
        let area = if particles {
            [
                (x[2] / 150.0).powf(2.0 / 3.0),
                (x[3] / 100.0).powf(2.0 / 3.0),
            ]
        } else {
            [1.0; 2]
        };
        let a = 0.1 * area[0] * (50.0 * (-0.001 * x[1]).exp() - x[0]);
        let b = 0.05 * area[1] * (80.0 * (-0.001 * x[0]).exp() - x[1]);
        [a, b, -a, -b]
    };
    let add = |x: [f64; 4], a: [f64; 4], h: f64| core::array::from_fn(|i| x[i] + h * a[i]);
    let h = duration / 32768.0;
    let mut x = [0.0, 10.0, 150.0, 100.0];
    for _ in 0..32768 {
        let a = rhs(x);
        let b = rhs(add(x, a, 0.5 * h));
        let c = rhs(add(x, b, 0.5 * h));
        let d = rhs(add(x, c, h));
        for i in 0..4 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn sequential_phase_steps_refine_to_independent_coupled_dynamics_and_descend_in_energy() {
    for particles in [false, true] {
        let expected = continuous(4.0, particles);
        let mut errors = Vec::new();
        for count in [8, 32] {
            let mut owner = ideal(&[0.0, 10.0], &[150.0, 100.0], &[50.0, 80.0], particles)
                .with_nonideal_thermodynamics(source(1.0, 300.0, vec![0.0, 0.001, 0.001, 0.0], 2))
                .unwrap();
            for _ in 0..count {
                let before = owner.solution().solute_moles().to_vec();
                owner
                    .advance(Time::new(4.0 / f64::from(count)), || false)
                    .unwrap();
                assert!(
                    owner
                        .solid_relative_inventory_free_energy_change_j(
                            V,
                            &before,
                            V,
                            owner.solution().solute_moles()
                        )
                        .unwrap()
                        < 0.0
                );
            }
            let actual = [
                owner.solution().solute_moles()[0] / V,
                owner.solution().solute_moles()[1] / V,
                owner.last_steps()[0].state().solid().value() / V,
                owner.last_steps()[1].state().solid().value() / V,
            ];
            errors.push(
                actual
                    .into_iter()
                    .zip(expected)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max),
            );
            close(actual[0] + actual[2], 150.0, 5e-15);
            close(actual[1] + actual[3], 110.0, 5e-15);
            // The first stage receipt retains its actual target/transfer even
            // though the second component subsequently shifts that target.
            assert_ne!(
                owner.last_steps()[0]
                    .equilibrium()
                    .saturation_concentration
                    .value(),
                owner
                    .conditional_saturation(0, &[actual[0], actual[1]])
                    .unwrap()
                    .value()
            );
            assert!(owner.last_steps()[0].transferred_to_solution().value() > 0.0);
        }
        assert!(errors[0] < 0.1 && errors[1] < 0.3 * errors[0], "{errors:?}");
    }
}

#[test]
fn one_calibrated_potential_supplies_phase_forces_work_and_osmotic_pressure() {
    let owner = ideal(&[30.0, 40.0], &[200.0, 100.0], &[50.0, 80.0], false)
        .with_nonideal_thermodynamics(source(1.0, 300.0, vec![0.01, 0.003, 0.003, 0.02], 2))
        .unwrap();
    let energy = |volume: f64, n: [f64; 2]| {
        let c = [n[0] / volume, n[1] / volume];
        R * T
            * (n[0] * ((c[0] / 50.0).ln() - 1.0 - 0.01 * 50.0)
                + n[1] * ((c[1] / 80.0).ln() - 1.0 - 0.02 * 80.0)
                + 0.5 * volume * (0.01 * c[0] * c[0] + 0.006 * c[0] * c[1] + 0.02 * c[1] * c[1]))
    };
    let before = [0.03, 0.04];
    let after = [0.035, 0.025];
    for volume in [V, 0.0009, 0.002] {
        close(
            owner
                .solid_relative_inventory_free_energy_change_j(V, &before, volume, &after)
                .unwrap(),
            energy(volume, after) - energy(V, before),
            2e-14,
        );
    }
    let h = V * 1e-6;
    let plus = owner
        .solid_relative_inventory_free_energy_change_j(V, &before, V + h, &before)
        .unwrap();
    let minus = owner
        .solid_relative_inventory_free_energy_change_j(V, &before, V - h, &before)
        .unwrap();
    close(
        -(plus - minus) / (2.0 * h),
        owner.solution().osmotic_pressure_pa().unwrap(),
        3e-10,
    );
    let mut potential = [None; 2];
    owner
        .solid_relative_chemical_potentials_into(&[30.0, 40.0], &mut potential)
        .unwrap();
    for i in 0..2 {
        let mut after = before;
        let mut prior = before;
        after[i] += 1e-8;
        prior[i] -= 1e-8;
        let work = owner
            .solid_relative_inventory_free_energy_change_j(V, &prior, V, &after)
            .unwrap();
        close(work / (2e-8), potential[i].unwrap(), 3e-10);
    }
}

#[test]
fn original_sources_traces_and_support_are_preserved_by_transport_and_rollback() {
    let matrix = vec![0.01, 0.003, 0.003, 0.02];
    let mut owner = ideal(&[10.0, 20.0], &[200.0, 100.0], &[50.0, 80.0], true)
        .with_nonideal_thermodynamics(source(1.0, 300.0, matrix.clone(), 2))
        .unwrap();
    let mut shifted = ideal(&[10.0, 20.0], &[200.0, 100.0], &[50.0, 80.0], true)
        .with_nonideal_thermodynamics(source(1e200, 300.0, matrix, 2))
        .unwrap();
    let before = owner.clone();
    let mut calls = 0;
    assert_eq!(
        owner.advance(Time::new(1.0), || {
            calls += 1;
            calls == 3
        }),
        Err(KineticSolubilityError::Cancelled)
    );
    assert_eq!(owner, before);
    owner.advance(Time::new(0.5), || false).unwrap();
    shifted.advance(Time::new(0.5), || false).unwrap();
    assert_eq!(owner.last_steps(), shifted.last_steps());
    assert_eq!(
        owner
            .solid_relative_inventory_free_energy_change_j(
                V,
                &[0.01, 0.02],
                V,
                owner.solution().solute_moles()
            )
            .unwrap(),
        shifted
            .solid_relative_inventory_free_energy_change_j(
                V,
                &[0.01, 0.02],
                V,
                shifted.solution().solute_moles()
            )
            .unwrap()
    );
    let diluted = owner.with_solvent_volume_m3(2.0 * V).unwrap();
    assert_eq!(
        diluted.solution().nonideal_model(),
        owner.solution().nonideal_model()
    );
    for (a, b) in diluted.last_steps().iter().zip(owner.last_steps()) {
        assert_eq!(a.state(), b.state());
    }
    assert!(owner.with_temperature(Temperature::new(310.0)).is_err());
    assert_eq!(owner.with_temperature(Temperature::new(T)).unwrap(), owner);
    assert!(owner.copy_state_from(&shifted).is_err());
    let old_solid = owner.last_steps()[0].state().solid();
    owner
        .replace_dissolved_inventories(V, &[0.015, 0.025])
        .unwrap();
    assert_eq!(owner.last_steps()[0].state().solid(), old_solid);
    assert_eq!(
        owner.models()[0].particle_reference_solid(),
        Some(Amount::new(0.2))
    );

    let mut limited = ideal(&[10.0, 20.0], &[200.0, 200.0], &[80.0, 80.0], false)
        .with_nonideal_thermodynamics(source(1.0, 100.0, vec![0.0, -0.004, -0.004, 0.0], 2))
        .unwrap();
    let original = limited.clone();
    assert!(limited.advance(Time::new(100.0), || false).is_err());
    assert_eq!(limited, original);
    assert!(
        limited
            .replace_dissolved_inventories(V, &[0.01, 0.1])
            .is_err()
    );
    assert_eq!(limited, original);
    assert!(
        ideal(&[10.0], &[100.0], &[101.0], false)
            .with_nonideal_thermodynamics(source(1.0, 100.0, vec![0.0], 1))
            .is_err()
    );
}

#[test]
fn conditional_root_keeps_positive_traces_and_refuses_unsupported_saturation() {
    close(
        conditional_root(1e-250, 0.0, 100.0, 1.0).unwrap(),
        1e-250 * (-100.0_f64).exp(),
        2e-14,
    );
    let trace = conditional_root(1e-250, 1e250, 10.0, 2e-250).unwrap();
    assert!(trace > 0.0 && trace < 1e-250);
    assert!(((trace / 1e-250).ln() + 1e250 * (trace - 1e-250) + 10.0).abs() < 1e-12);
    assert!(conditional_root(1e-300, 0.0, 100.0, 1.0).is_err());
    assert!(conditional_root(80.0, 0.0, -0.4, 100.0).is_err());
}
