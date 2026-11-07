use super::*;
use fs_material::osmosis::{IdealDiluteSolution, nonideal::VirialSolutionModel};
use std::sync::Arc;

const TEMPERATURE: f64 = 300.0;
const RT: f64 = R * TEMPERATURE;
const SATURATION: f64 = 1e-5;

fn law() -> SurfactantLaw {
    SurfactantLaw::Langmuir {
        clean_tension_n_m: 0.2,
        thermal_energy_j_mol: RT,
        saturation_concentration_mol_m2: SATURATION,
    }
}
fn bulk(reference: f64, amounts: Vec<f64>, interactions: Vec<f64>) -> AqueousSolution {
    let model = VirialSolutionModel::new(
        TEMPERATURE,
        18e-6,
        reference,
        vec![20.0; amounts.len()],
        interactions,
    )
    .unwrap();
    AqueousSolution::nonideal(Arc::new(model), 1e-6, amounts).unwrap()
}
fn controls(step: f64, faces: usize, edges: usize) -> SolubleSurfaceControls {
    SolubleSurfaceControls {
        maximum_step_s: step,
        maximum_substeps: 2_000,
        maximum_root_evaluations: 200_000,
        maximum_surface_pair_exchanges: 20_000,
        maximum_faces: faces,
        maximum_edges: edges,
        root_tolerance_mol: 1e-20,
        amount_tolerance_mol: 1e-15,
        energy_tolerance_j: 1e-11,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (allowance {tolerance:e})"
    );
}
fn energy(solution: &AqueousSolution, areas: &[f64], amounts: &[f64], offset: f64) -> f64 {
    let mut concentration = vec![0.0; solution.solute_moles().len()];
    solution.concentrations_into(&mut concentration).unwrap();
    solution.solvent_volume_m3()
        * solution
            .nonideal_model()
            .unwrap()
            .mixing_free_energy_density_j_m3(&concentration)
            .unwrap()
        + areas
            .iter()
            .zip(amounts)
            .map(|(&area, &amount)| {
                SurfactantFace {
                    law: law(),
                    amount_mol: amount,
                }
                .at_area(area)
                .unwrap()
                .free_energy_j
                    + offset * amount
            })
            .sum::<f64>()
}

#[test]
fn finite_bath_adsorption_and_desorption_match_independent_implicit_quadratic() {
    let area = 0.2;
    let capacity = area * SATURATION;
    let reference = 2.0;
    let offset = RT * 2.0_f64.ln();
    let kd = 3.0;
    let dt = 0.4;
    for (mobile, occupied) in [(4e-6, 0.0), (0.0, 1e-6)] {
        let mut solution = bulk(reference, vec![mobile, 1e-6], vec![0.0; 4]);
        let mut amounts = [occupied];
        let initial_energy = energy(&solution, &[area], &amounts, offset);
        let total = mobile + occupied;
        let a = kd * (-offset / RT).exp() / (solution.solvent_volume_m3() * reference);
        // h*a*n² - (1+h*kd+h*a*(M+C))*n + n0+h*a*M*C = 0.
        let b = 1.0 + dt * kd + dt * a * (total + capacity);
        let c = occupied + dt * a * total * capacity;
        let expected = 2.0 * c / (b + (b * b - 4.0 * dt * a * c).sqrt());
        let report = advance_soluble_surface(
            &mut solution,
            &[law()],
            &[area],
            &mut amounts,
            SolubleSurfaceKinetics {
                bulk_species: 0,
                standard_potential_offset_j_mol: offset,
                desorption_rates_s_inverse: &[kd],
            },
            &[],
            dt,
            controls(dt, 1, 0),
            || false,
        )
        .unwrap();
        close(amounts[0], expected, 3e-20);
        close(solution.solute_moles()[0], total - expected, 3e-20);
        assert_eq!(solution.solute_moles()[1], 1e-6);
        assert!(amounts[0] > 0.0 && amounts[0] < capacity);
        assert!(solution.solute_moles()[0] > 0.0);
        close(
            report.gross_adsorption_mol - report.gross_desorption_mol,
            amounts[0] - occupied,
            3e-20,
        );
        close(
            report.free_energy_change_j,
            energy(&solution, &[area], &amounts, offset) - initial_energy,
            2e-16,
        );
        assert!(report.exchange_free_energy_dissipation_j > 0.0);
        assert_eq!(report.surface_pair_exchanges, 0);
        assert_eq!(report.surface_diffusion_free_energy_dissipation_j, 0.0);
    }
}

#[test]
fn finite_bath_kinetics_converges_at_first_order_to_the_exact_riccati_flow() {
    let area = 0.2;
    let capacity = area * SATURATION;
    let mobile = 3e-6;
    let initial_surface = 0.1e-6;
    let total = mobile + initial_surface;
    let kd = 2.0;
    let duration = 0.7;
    let a: f64 = kd / 1e-6;
    let b: f64 = kd + a * (total + capacity);
    let c: f64 = a * total * capacity;
    let discriminant = (b * b - 4.0 * a * c).sqrt();
    let low = 2.0 * c / (b + discriminant);
    let high = (b + discriminant) / (2.0 * a);
    let ratio =
        (initial_surface - low) / (initial_surface - high) * (-discriminant * duration).exp();
    let exact = (low - ratio * high) / (1.0 - ratio);
    let errors = [16, 32, 64].map(|steps| {
        let mut solution = bulk(1.0, vec![mobile], vec![0.0]);
        let mut amounts = [initial_surface];
        advance_soluble_surface(
            &mut solution,
            &[law()],
            &[area],
            &mut amounts,
            SolubleSurfaceKinetics {
                bulk_species: 0,
                standard_potential_offset_j_mol: 0.0,
                desorption_rates_s_inverse: &[kd],
            },
            &[],
            duration,
            controls(duration / steps as f64, 1, 0),
            || false,
        )
        .unwrap();
        (amounts[0] - exact).abs()
    });
    assert!(
        errors[0] / errors[1] > 1.8 && errors[0] / errors[1] < 2.2,
        "{errors:?}"
    );
    assert!(
        errors[1] / errors[2] > 1.8 && errors[1] / errors[2] < 2.2,
        "{errors:?}"
    );
}

#[test]
fn finite_equilibrium_retains_gross_turnover_and_zero_time_retains_both_inventories() {
    let original = bulk(1.0, vec![1e-6], vec![0.0]);
    let initial_surface = [0.5 * (0.2 * SATURATION)];
    let kinetics = SolubleSurfaceKinetics {
        bulk_species: 0,
        standard_potential_offset_j_mol: 0.0,
        desorption_rates_s_inverse: &[1.0],
    };
    for duration in [0.0, 1e6] {
        let mut solution = original.clone();
        let mut amounts = initial_surface;
        let mut limits = controls(1e6, 1, 0);
        limits.maximum_root_evaluations = 1;
        limits.maximum_surface_pair_exchanges = 0;
        let report = advance_soluble_surface(
            &mut solution,
            &[law()],
            &[0.2],
            &mut amounts,
            kinetics,
            &[],
            duration,
            limits,
            || false,
        )
        .unwrap();
        assert_eq!(solution, original);
        assert_eq!(amounts, initial_surface);
        assert_eq!(report.free_energy_change_j, 0.0);
        assert_eq!(report.gross_adsorption_mol, report.gross_desorption_mol);
        if duration > 0.0 {
            close(report.gross_adsorption_mol, 1.0, 2e-15);
            assert_eq!(report.root_evaluations, 1);
        } else {
            assert_eq!(report.gross_adsorption_mol, 0.0);
            assert_eq!(report.root_evaluations, 0);
        }
    }
}

#[test]
fn nonideal_activity_co_solutes_and_chemical_reference_are_retained() {
    let interactions = vec![0.1, 0.03, 0.03, 0.02];
    let areas = [0.2, 0.3];
    let old = [0.2e-6, 0.9e-6];
    let kd = [0.4, 1.1];
    let offset = RT * 0.5_f64.ln();
    let dt = 0.3;
    let mut results = Vec::new();
    for reference in [1.0, 4.0] {
        let mut solution = bulk(reference, vec![3e-6, 2e-6], interactions.clone());
        let original_model = solution.nonideal_model().unwrap().clone();
        let shifted_offset = offset - RT * reference.ln();
        let mut amounts = old;
        let before = energy(&solution, &areas, &amounts, shifted_offset);
        let report = advance_soluble_surface(
            &mut solution,
            &[law(); 2],
            &areas,
            &mut amounts,
            SolubleSurfaceKinetics {
                bulk_species: 0,
                standard_potential_offset_j_mol: shifted_offset,
                desorption_rates_s_inverse: &kd,
            },
            &[],
            dt,
            controls(dt, 2, 0),
            || false,
        )
        .unwrap();
        assert_eq!(solution.nonideal_model().unwrap(), &original_model);
        assert_eq!(solution.solute_moles()[1], 2e-6);
        let concentration = solution.solute_moles()[0] / solution.solvent_volume_m3();
        let activity = concentration / reference * (0.1 * concentration + 0.03 * 2.0).exp();
        let factor = activity * (-shifted_offset / RT).exp();
        for i in 0..2 {
            let expected = (old[i] + dt * kd[i] * factor * areas[i] * SATURATION)
                / (1.0 + dt * kd[i] * (1.0 + factor));
            close(amounts[i], expected, 2e-20);
        }
        close(
            solution.solute_moles()[0] + amounts.iter().sum::<f64>(),
            3e-6 + old.iter().sum::<f64>(),
            2e-20,
        );
        close(
            report.free_energy_change_j,
            energy(&solution, &areas, &amounts, shifted_offset) - before,
            3e-16,
        );
        results.push((
            solution.solute_moles()[0],
            amounts,
            report.free_energy_change_j,
        ));
    }
    close(results[0].0, results[1].0, 3e-20);
    for i in 0..2 {
        close(results[0].1[i], results[1].1[i], 3e-20);
    }
    close(results[0].2, results[1].2, 4e-16);

    let mut ideal = bulk(1.0, vec![3e-6, 2e-6], vec![0.0; 4]);
    let mut amounts = old;
    advance_soluble_surface(
        &mut ideal,
        &[law(); 2],
        &areas,
        &mut amounts,
        SolubleSurfaceKinetics {
            bulk_species: 0,
            standard_potential_offset_j_mol: offset,
            desorption_rates_s_inverse: &kd,
        },
        &[],
        dt,
        controls(dt, 2, 0),
        || false,
    )
    .unwrap();
    assert!(amounts.iter().sum::<f64>() < results[0].1.iter().sum::<f64>());
}

#[test]
fn lateral_diffusion_supplies_sealed_faces_from_one_finite_adsorbing_bath() {
    let areas = [0.1, 0.2, 0.3];
    let mut amounts = [0.05e-6, 1e-6, 2.7e-6];
    let mut solution = bulk(1.0, vec![1e-6], vec![0.0]);
    let total = solution.solute_moles()[0] + amounts.iter().sum::<f64>();
    let capacity = areas.iter().sum::<f64>() * SATURATION;
    let bath_scale = solution.solvent_volume_m3();
    // At equilibrium N_bulk=B*theta/(1-theta) and N_surface=C*theta.
    let b = total + capacity + bath_scale;
    let theta = 2.0 * total / (b + (b * b - 4.0 * capacity * total).sqrt());
    let before = energy(&solution, &areas, &amounts, 0.0);
    let edges = [
        SurfacePartitionDiffusionEdge {
            faces: [0, 1],
            conductance_mol_s: 4e-6,
        },
        SurfacePartitionDiffusionEdge {
            faces: [1, 2],
            conductance_mol_s: 2e-6,
        },
    ];
    let report = advance_soluble_surface(
        &mut solution,
        &[law(); 3],
        &areas,
        &mut amounts,
        SolubleSurfaceKinetics {
            bulk_species: 0,
            standard_potential_offset_j_mol: 0.0,
            desorption_rates_s_inverse: &[2.0, 0.0, 0.0],
        },
        &edges,
        100.0,
        controls(0.1, 3, 2),
        || false,
    )
    .unwrap();
    for i in 0..3 {
        close(amounts[i] / (areas[i] * SATURATION), theta, 2e-7);
        assert!(amounts[i] > 0.0 && amounts[i] < areas[i] * SATURATION);
    }
    close(
        solution.solute_moles()[0] + amounts.iter().sum::<f64>(),
        total,
        1e-15,
    );
    close(
        report.free_energy_change_j,
        energy(&solution, &areas, &amounts, 0.0) - before,
        1e-13,
    );
    assert!(report.exchange_free_energy_dissipation_j > 0.0);
    assert!(report.surface_diffusion_free_energy_dissipation_j > 0.0);
    assert!(report.gross_adsorption_mol > 0.0 && report.gross_desorption_mol > 0.0);
    assert_eq!(report.surface_pair_exchanges, 6_000);
    close(
        report.free_energy_change_j
            + report.exchange_free_energy_dissipation_j
            + report.surface_diffusion_free_energy_dissipation_j,
        0.0,
        1e-11,
    );
}

#[test]
fn source_support_tension_work_and_late_cancellation_refuse_atomically() {
    let original = bulk(1.0, vec![3e-6], vec![0.0]);
    let original_amounts = [0.1e-6, 0.7e-6];
    let areas = [0.2, 0.3];
    let edges = [SurfacePartitionDiffusionEdge {
        faces: [0, 1],
        conductance_mol_s: 1e-6,
    }];
    let kinetics = SolubleSurfaceKinetics {
        bulk_species: 0,
        standard_potential_offset_j_mol: 0.0,
        desorption_rates_s_inverse: &[1.0, 2.0],
    };
    let limits = controls(0.1, 2, 1);
    let mut accepted = original.clone();
    let mut accepted_amounts = original_amounts;
    let mut polls = 0;
    let report = advance_soluble_surface(
        &mut accepted,
        &[law(); 2],
        &areas,
        &mut accepted_amounts,
        kinetics,
        &edges,
        0.3,
        limits,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(report.root_evaluations > 3);
    let mut retry = original.clone();
    let mut amounts = original_amounts;
    let mut restricted = limits;
    restricted.maximum_root_evaluations = report.root_evaluations - 1;
    assert!(matches!(
        advance_soluble_surface(
            &mut retry,
            &[law(); 2],
            &areas,
            &mut amounts,
            kinetics,
            &edges,
            0.3,
            restricted,
            || false
        ),
        Err(SolubleSurfaceError::RootBudget { .. })
    ));
    assert_eq!(retry, original);
    assert_eq!(amounts, original_amounts);
    let mut current = 0;
    assert!(matches!(
        advance_soluble_surface(
            &mut retry,
            &[law(); 2],
            &areas,
            &mut amounts,
            kinetics,
            &edges,
            0.3,
            limits,
            || {
                current += 1;
                current == polls
            }
        ),
        Err(SolubleSurfaceError::Cancelled)
    ));
    assert_eq!(retry, original);
    assert_eq!(amounts, original_amounts);

    let datum_free: AqueousSolution =
        IdealDiluteSolution::new(TEMPERATURE, 18e-6, 1e-6, vec![3e-6])
            .unwrap()
            .into();
    let mut refused = datum_free.clone();
    assert!(matches!(
        advance_soluble_surface(
            &mut refused,
            &[law(); 2],
            &areas,
            &mut amounts,
            kinetics,
            &[],
            0.1,
            limits,
            || false
        ),
        Err(SolubleSurfaceError::Invalid(_))
    ));
    assert_eq!(refused, datum_free);
    let model = VirialSolutionModel::new(TEMPERATURE, 18e-6, 1.0, vec![0.01], vec![0.0]).unwrap();
    let mut limited = AqueousSolution::nonideal(Arc::new(model), 1e-6, vec![0.0]).unwrap();
    let original_limited = limited.clone();
    let mut loaded = [1e-6];
    assert!(matches!(
        advance_soluble_surface(
            &mut limited,
            &[law()],
            &[0.2],
            &mut loaded,
            SolubleSurfaceKinetics {
                bulk_species: 0,
                standard_potential_offset_j_mol: RT * 100.0_f64.ln(),
                desorption_rates_s_inverse: &[2.0]
            },
            &[],
            1.0,
            controls(1.0, 1, 0),
            || false
        ),
        Err(SolubleSurfaceError::BulkSupport { .. })
    ));
    assert_eq!(limited, original_limited);
    assert_eq!(loaded, [1e-6]);
    let weak_tension = SurfactantLaw::Langmuir {
        clean_tension_n_m: 0.001,
        thermal_energy_j_mol: RT,
        saturation_concentration_mol_m2: SATURATION,
    };
    let mut clean = [0.0];
    assert!(matches!(
        advance_soluble_surface(
            &mut retry,
            &[weak_tension],
            &[0.2],
            &mut clean,
            SolubleSurfaceKinetics {
                bulk_species: 0,
                standard_potential_offset_j_mol: 0.0,
                desorption_rates_s_inverse: &[2.0]
            },
            &[],
            1.0,
            controls(1.0, 1, 0),
            || false
        ),
        Err(SolubleSurfaceError::Surface {
            error: SurfactantError::NegativeTension,
            ..
        })
    ));
    assert_eq!(retry, original);
    assert_eq!(clean, [0.0]);
}
