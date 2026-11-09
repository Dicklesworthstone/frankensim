//! Conserved crystals, original MS transport and solid-relative free energy.

use super::*;
use crate::osmosis::coupled::CoupledSolutionControls;
use fs_material::osmosis::solubility::kinetics::KineticSolubilitySolution;
use fs_qty::{Amount, AmountConcentration, Temperature, Volume, VolumetricFlowRate};
use fs_thermochem::{
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    solution::solubility::{
        IdealSolubilityV1,
        kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
    },
};

fn kinetic_cell(
    dissolved: &[f64],
    solid: &[f64],
    saturation: &[f64],
    rates: &[f64],
    particles: bool,
) -> OsmoticCompartment {
    let components = dissolved
        .iter()
        .zip(solid)
        .zip(saturation)
        .zip(rates)
        .enumerate()
        .map(|(index, (((&dissolved, &solid), &saturation), &rate))| {
            let law = IdealSolubilityV1::new(
                SpeciesId::new("water").unwrap(),
                SpeciesId::new(format!("solute-{index}")).unwrap(),
                Temperature::new(T),
                AmountConcentration::new(saturation),
                0.0,
                [Temperature::new(280.0), Temperature::new(350.0)],
            )
            .unwrap();
            let model = SolubilityKineticsV1::new(
                law,
                VolumetricFlowRate::new(rate * VOLUME),
                PrecipitationPolicyV1::ExistingSolidOnly,
            )
            .unwrap();
            let model = if particles {
                model
                    .with_particle_surface(Amount::new(200.0 * VOLUME))
                    .unwrap()
            } else {
                model
            };
            (
                model,
                SolubilityKineticStateV1::from_phases(
                    Amount::new(dissolved * VOLUME),
                    Amount::new(solid * VOLUME),
                )
                .unwrap(),
            )
        })
        .collect();
    let original =
        KineticSolubilitySolution::new(Temperature::new(T), VBAR, Volume::new(VOLUME), components)
            .unwrap();
    OsmoticCompartment::with_kinetic_solubility(
        original,
        OsmoticPressureVolume::new(VOLUME, 100000.0, 1e8).unwrap(),
        0.1 * VOLUME,
    )
    .unwrap()
}

fn prepared(reference: f64, particles: bool) -> OsmoticNetwork {
    let model = source(reference, vec![0.0], vec![300.0]);
    OsmoticNetwork::with_kinetic_solution_diffusion(
        vec!["solute-0".into()],
        vec![
            kinetic_cell(&[0.0], &[200.0], &[50.0], &[0.1], particles),
            cell(&model, &[0.0]),
        ],
        vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()],
        model,
        physical(),
        1 << 22,
    )
    .unwrap()
}

fn complete_bits(system: &OsmoticNetwork) -> Vec<u64> {
    let mut result = bits(system);
    for cell in system.compartments() {
        if let Some(phases) = cell.kinetic_solubility_snapshot() {
            for phase in phases.last_steps() {
                result.push(phase.state().solid().value().to_bits());
                result.push(phase.state().total_solute().value().to_bits());
                result.push(phase.transferred_to_solution().value().to_bits());
            }
        }
    }
    result
}

fn phases(system: &OsmoticNetwork) -> KineticSolubilitySolution {
    system.compartments()[0]
        .kinetic_solubility_snapshot()
        .unwrap()
}

// Independent continuous equations in concentration units. The original
// binary MS law in stationary solvent gives N = D*dc*(1+c_log*Vbar)/L;
// the particle source separately obeys the changing-area finite-bath law.
fn oracle(duration: f64, particles: bool) -> [f64; 3] {
    let rhs = |x: [f64; 3]| {
        let [left, right, solid] = x;
        let mean = if left == 0.0 || right == 0.0 {
            0.0
        } else if left == right {
            left
        } else {
            (left - right) / (left / right).ln()
        };
        let molecular = 0.01 * (left - right) * (1.0 + mean * VBAR);
        let area = if particles {
            (solid / 200.0).powf(2.0 / 3.0)
        } else {
            1.0
        };
        let dissolution = 0.1 * area * (50.0 - left);
        [dissolution - molecular, molecular, -dissolution]
    };
    let add = |x: [f64; 3], a: [f64; 3], scale: f64| {
        core::array::from_fn(|index| x[index] + scale * a[index])
    };
    let mut x = [0.0, 0.0, 200.0];
    let h = duration / 32768.0;
    for _ in 0..32768 {
        let a = rhs(x);
        let b = rhs(add(x, a, 0.5 * h));
        let c = rhs(add(x, b, 0.5 * h));
        let d = rhs(add(x, c, h));
        for index in 0..3 {
            x[index] += h * (a[index] + 2.0 * b[index] + 2.0 * c[index] + d[index]) / 6.0;
        }
    }
    x
}

#[test]
fn particle_dissolution_and_native_molecular_diffusion_refine_to_independent_ode() {
    for particles in [false, true] {
        let exact = oracle(5.0, particles);
        let mut errors = Vec::new();
        for relative in [1e-3, 6.25e-5] {
            let mut system = prepared(1.0, particles);
            let report = system.advance(5.0, numerical(relative), || false).unwrap();
            let retained = phases(&system);
            let actual = [
                concentration(&system, 0, 0),
                concentration(&system, 1, 0),
                retained.last_steps()[0].state().solid().value() / VOLUME,
            ];
            errors.push(
                actual
                    .iter()
                    .zip(exact)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max),
            );
            assert!(report.phase_evaluations >= 6 * report.accepted_steps);
            assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
            assert_eq!(report.water_balance_residual_m3, 0.0);
            assert!(
                report.solute_balance_residual_mol[0].abs()
                    <= physical().amount_balance_tolerance_mol
            );
            assert!((actual.iter().sum::<f64>() - 200.0).abs() < 1e-9);
            assert_eq!(
                retained.solution().solute_moles(),
                system.compartments()[0].solution().solute_moles()
            );
            if particles {
                assert_eq!(
                    retained.models()[0].particle_reference_solid(),
                    Some(Amount::new(200.0 * VOLUME))
                );
                assert!(
                    retained.models()[0]
                        .conductance_at(retained.last_steps()[0].state().solid())
                        .unwrap()
                        .value()
                        < 0.1 * VOLUME
                );
            }
        }
        assert!(errors[0] < 0.1 && errors[1] < 0.6 * errors[0], "{errors:?}");
    }
}

#[test]
fn solid_chemical_reference_allows_dissolution_while_bare_mixing_energy_increases() {
    let model = source(1.0, vec![0.0], vec![300.0]);
    let mut system = OsmoticNetwork::with_kinetic_solution_diffusion(
        vec!["solute-0".into()],
        vec![kinetic_cell(&[10.0], &[100.0], &[50.0], &[0.1], false)],
        Vec::new(),
        model,
        physical(),
        1 << 22,
    )
    .unwrap();
    let report = system.advance(1.0, numerical(1e-5), || false).unwrap();
    let final_c = concentration(&system, 0, 0);
    let mixing = |c: f64| R * T * VOLUME * c * (c.ln() - 1.0);
    assert!(mixing(final_c) > mixing(10.0));
    let complete = |c: f64| R * T * VOLUME * c * ((c / 50.0).ln() - 1.0);
    let expected = complete(final_c) - complete(10.0);
    assert!(expected < 0.0);
    assert!((report.diffusion_free_energy_change_j.unwrap() - expected).abs() < 1e-12);
    assert!(report.solute_balance_residual_mol[0].abs() <= physical().amount_balance_tolerance_mol);
}

#[test]
fn original_reference_and_existing_implicit_binary_stages_preserve_complete_phases() {
    let mut ordinary = prepared(1.0, true);
    let mut shifted = prepared(1e200, true);
    let a = ordinary.advance(2.0, numerical(1e-4), || false).unwrap();
    let b = shifted.advance(2.0, numerical(1e-4), || false).unwrap();
    assert_eq!(complete_bits(&ordinary), complete_bits(&shifted));
    assert_eq!(
        a.diffusion_free_energy_change_j,
        b.diffusion_free_energy_change_j
    );
    let mut implicit = prepared(1.0, true)
        .with_implicit_binary_diffusion(ImplicitBinaryDiffusionControls::default(), 1 << 22)
        .unwrap();
    let report = implicit.advance(2.0, numerical(1e-5), || false).unwrap();
    let exact = oracle(2.0, true);
    assert!((concentration(&implicit, 0, 0) - exact[0]).abs() < 0.01);
    assert!((concentration(&implicit, 1, 0) - exact[1]).abs() < 0.01);
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
    assert!(report.solute_balance_residual_mol[0].abs() <= physical().amount_balance_tolerance_mol);
}

fn two_component() -> OsmoticNetwork {
    let model = source(1.0, vec![0.0; 4], vec![300.0; 2]);
    OsmoticNetwork::with_kinetic_solution_diffusion(
        vec!["solute-0".into(), "solute-1".into()],
        vec![
            kinetic_cell(
                &[2.0, 1.0],
                &[200.0, 60.0],
                &[50.0, 20.0],
                &[0.1, 0.03],
                true,
            ),
            cell(&model, &[1.0, 2.0]),
        ],
        vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()],
        model,
        physical(),
        1 << 24,
    )
    .unwrap()
}

#[test]
fn shared_multicomponent_friction_transports_only_dissolved_stocks() {
    for implicit in [false, true] {
        let system = two_component();
        let mut system = if implicit {
            system
                .with_implicit_multicomponent_diffusion(
                    ImplicitMulticomponentDiffusionControls::default(),
                    1 << 24,
                )
                .unwrap()
        } else {
            system
        };
        let report = system.advance(1.0, numerical(1e-4), || false).unwrap();
        let retained = phases(&system);
        for (index, total) in [203.0, 63.0].into_iter().enumerate() {
            let phase = retained.last_steps()[index].state();
            let actual = phase.solid().value()
                + phase.dissolved().value()
                + system.compartments()[1].solution().solute_moles()[index];
            assert!((actual / VOLUME - total).abs() < 1e-8);
            assert!(
                report.solute_balance_residual_mol[index].abs()
                    <= physical().amount_balance_tolerance_mol
            );
        }
        assert!(
            report.diffusion_free_energy_change_j.unwrap() <= physical().free_energy_tolerance_j
        );
    }
}

#[test]
fn cancelled_or_exhausted_work_keeps_both_phase_stocks_and_clock_atomic() {
    let mut system = prepared(1.0, true);
    let before = complete_bits(&system);
    let mut polls = 0;
    let error = system.advance(5.0, numerical(1e-5), || {
        polls += 1;
        polls == 64
    });
    assert!(matches!(error, Err(OsmoticNetworkError::Cancelled)));
    assert_eq!(before, complete_bits(&system));
    let mut limited = numerical(1e-5);
    limited.maximum_step_s = 0.01;
    limited.maximum_steps = 1;
    assert!(matches!(
        system.advance(1.0, limited, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(before, complete_bits(&system));
    let mut fresh = prepared(1.0, true);
    system.advance(1.0, numerical(1e-5), || false).unwrap();
    fresh.advance(1.0, numerical(1e-5), || false).unwrap();
    assert_eq!(complete_bits(&system), complete_bits(&fresh));
}

#[test]
fn unbound_nonideal_saturation_incompatible_references_and_ionic_or_water_modes_refuse() {
    let interacting = source(1.0, vec![0.001], vec![300.0]);
    assert!(
        OsmoticNetwork::with_kinetic_solution_diffusion(
            vec!["solute-0".into()],
            vec![kinetic_cell(&[1.0], &[100.0], &[50.0], &[0.1], false)],
            Vec::new(),
            interacting,
            physical(),
            1 << 22,
        )
        .is_err()
    );
    let model = source(1.0, vec![0.0], vec![300.0]);
    assert!(
        OsmoticNetwork::with_kinetic_solution_diffusion(
            vec!["solute-0".into()],
            vec![
                kinetic_cell(&[1.0], &[100.0], &[50.0], &[0.1], false),
                kinetic_cell(&[1.0], &[100.0], &[60.0], &[0.1], false),
            ],
            Vec::new(),
            model.clone(),
            physical(),
            1 << 22,
        )
        .is_err()
    );
    let other_solvent = SolutionMaxwellStefan::new(
        model.clone(),
        (
            SpeciesId::new("other-solvent").unwrap(),
            MolarMass::new(0.018),
        ),
        vec![(SpeciesId::new("solute-0").unwrap(), MolarMass::new(0.04))],
        vec![0.0, 1e-9, 1e-9, 0.0],
    )
    .unwrap();
    assert!(
        OsmoticNetwork::with_kinetic_solution_diffusion(
            vec!["solute-0".into()],
            vec![
                kinetic_cell(&[1.0], &[100.0], &[50.0], &[0.1], false),
                cell(&model, &[1.0])
            ],
            vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, other_solvent).unwrap()],
            model,
            physical(),
            1 << 22,
        )
        .is_err()
    );
    let ionic = || ZeroCurrentDiffusion::new(vec![1, -1], 1e-15, 1e-16).unwrap();
    assert!(
        two_component()
            .with_zero_current_diffusion(ionic(), 1 << 24)
            .is_err()
    );
    assert!(
        two_component()
            .with_electroneutral_diffusion(
                ionic(),
                ElectroneutralDiffusionControls::default(),
                1 << 24
            )
            .is_err()
    );
    let original = prepared(1.0, false);
    assert!(
        OsmoticNetwork::with_coupled_solution_transport(
            vec!["solute-0".into()],
            original.cells.clone(),
            Vec::new(),
            original.diffusion.clone(),
            CoupledSolutionControls {
                diffusion: physical(),
                water_balance_tolerance_m3: 1e-18
            },
            1 << 22,
        )
        .is_err()
    );
}

#[test]
fn explicitly_bound_nonideal_phases_share_molecular_forces_and_solid_reference() {
    let model = source(1.0, vec![0.002, 0.003, 0.003, 0.001], vec![300.0; 2]);
    let bound_phases = kinetic_cell(
        &[10.0, 40.0],
        &[200.0, 0.0],
        &[50.0, 20.0],
        &[0.1, 0.0],
        true,
    )
    .kinetic_solubility_snapshot()
    .unwrap()
    .with_nonideal_thermodynamics(model.clone())
    .unwrap();
    let left = OsmoticCompartment::with_kinetic_solubility(
        bound_phases,
        OsmoticPressureVolume::new(VOLUME, 1e5, 1e8).unwrap(),
        0.1 * VOLUME,
    )
    .unwrap();
    let initial = [vec![10.0, 40.0], vec![1.0, 10.0]];
    let energy = |c: &[f64]| {
        let ideal = c
            .iter()
            .zip([50.0, 20.0])
            .zip([0.002, 0.001])
            .map(|((&c, saturation), diagonal)| {
                c * ((c / saturation).ln() - 1.0 - diagonal * saturation)
            })
            .sum::<f64>();
        R * T
            * VOLUME
            * (ideal
                + 0.5 * (0.002 * c[0] * c[0] + 2.0 * 0.003 * c[0] * c[1] + 0.001 * c[1] * c[1]))
    };
    for implicit in [false, true] {
        let system = OsmoticNetwork::with_kinetic_solution_diffusion(
            vec!["solute-0".into(), "solute-1".into()],
            vec![left.clone(), cell(&model, &initial[1])],
            vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()],
            model.clone(),
            physical(),
            1 << 24,
        )
        .unwrap();
        let mut system = if implicit {
            system
                .with_implicit_multicomponent_diffusion(
                    ImplicitMulticomponentDiffusionControls::default(),
                    1 << 24,
                )
                .unwrap()
        } else {
            system
        };
        let report = system.advance(1.0, numerical(1e-4), || false).unwrap();
        let retained = phases(&system);
        assert_eq!(retained.solution().nonideal_model(), Some(model.as_ref()));
        assert_eq!(retained.last_steps()[1].state().solid().value(), 0.0);
        assert!(retained.last_steps()[0].state().solid().value() < 200.0 * VOLUME);
        for (i, total) in [211.0, 50.0].into_iter().enumerate() {
            let actual = retained.last_steps()[i].state().solid().value()
                + retained.last_steps()[i].state().dissolved().value()
                + system.compartments()[1].solution().solute_moles()[i];
            assert!((actual / VOLUME - total).abs() < 1e-8);
            assert!(
                report.solute_balance_residual_mol[i].abs()
                    <= physical().amount_balance_tolerance_mol
            );
        }
        let actual: [Vec<f64>; 2] = core::array::from_fn(|cell| {
            (0..2)
                .map(|species| concentration(&system, cell, species))
                .collect()
        });
        let expected =
            energy(&actual[0]) + energy(&actual[1]) - energy(&initial[0]) - energy(&initial[1]);
        assert!(expected < 0.0);
        assert!((report.diffusion_free_energy_change_j.unwrap() - expected).abs() < 2e-12);
    }
}
