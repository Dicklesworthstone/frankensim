use super::*;

fn ideal() -> SurfactantLaw {
    SurfactantLaw::IdealDilute {
        clean_tension_n_m: 10.0,
        thermal_energy_j_mol: 0.6,
        reference_concentration_mol_m2: 1.0,
    }
}

fn controls(step: f64) -> SurfaceDiffusionControls {
    SurfaceDiffusionControls {
        maximum_step_s: step,
        maximum_substeps: 256,
        maximum_faces: 16,
        maximum_edges: 32,
        amount_tolerance_mol: 1e-12,
        energy_increase_tolerance_j: 1e-12,
    }
}

fn edge(i: usize, j: usize, conductance: f64) -> SurfaceDiffusionEdge {
    SurfaceDiffusionEdge {
        faces: [i, j],
        conductance_m2_s: conductance,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_unequal_area_pair_matches_the_exact_diffusion_mode_and_energy_change() {
    let areas = [2.0, 5.0];
    let before = [1.6, 0.5];
    let mut amounts = before;
    let duration = 0.3;
    let conductance = 0.7;
    let report = diffuse_surface_surfactant(
        ideal(),
        &areas,
        &mut amounts,
        &[edge(0, 1, conductance)],
        duration,
        controls(duration),
        || false,
    )
    .unwrap();
    let decay = (-conductance * (1.0 / areas[0] + 1.0 / areas[1]) * duration).exp();
    close(amounts[0], 2.0 * (0.3 + 0.5 * decay), 8e-16);
    close(amounts[1], 5.0 * (0.3 - 0.2 * decay), 8e-16);
    close(report.amount_before_mol, 2.1, 5e-16);
    close(report.amount_after_mol, 2.1, 5e-16);
    close(report.amount_defect_mol, 0.0, 5e-16);
    assert_eq!(report.substeps, 1);
    assert_eq!(report.pair_exchanges, 1);
    assert!(report.free_energy_dissipation_j > 0.0);
    close(
        report.free_energy_change_j,
        report.free_energy_after_j - report.free_energy_before_j,
        8e-15,
    );
    for i in 0..2 {
        assert!((0.1..=0.8).contains(&(amounts[i] / areas[i])));
    }
}

#[test]
fn g0_stiff_langmuir_exchange_preserves_packing_mass_and_free_energy() {
    let law = SurfactantLaw::Langmuir {
        clean_tension_n_m: 4.0,
        thermal_energy_j_mol: 0.7,
        saturation_concentration_mol_m2: 1.0,
    };
    let areas = [1.0, 2.0, 0.5, 3.0, 1.0];
    let mut amounts = [0.05, 1.5, 0.06, 2.4, 0.33];
    let initial = amounts;
    let edges = [edge(0, 1, 1e6), edge(1, 2, 2e6), edge(2, 3, 5e5)];
    let mut losses = 0.0;
    for _ in 0..8 {
        let report = diffuse_surface_surfactant(
            law,
            &areas,
            &mut amounts,
            &edges,
            0.4,
            controls(0.2),
            || false,
        )
        .unwrap();
        assert!(report.free_energy_change_j <= 1e-13);
        losses += report.free_energy_dissipation_j;
        close(report.amount_defect_mol, 0.0, 1e-14);
        assert_eq!(amounts[4].to_bits(), initial[4].to_bits());
        for i in 0..4 {
            let concentration = amounts[i] / areas[i];
            assert!((0.05 - 1e-14..=0.8 + 1e-14).contains(&concentration));
            assert!(
                SurfactantFace {
                    law,
                    amount_mol: amounts[i]
                }
                .at_area(areas[i])
                .is_ok()
            );
        }
    }
    close(amounts[..4].iter().sum(), initial[..4].iter().sum(), 3e-14);
    assert!(losses > 0.1);
}

#[test]
fn g1_palindromic_sweeps_converge_at_second_order_to_a_three_face_spectrum() {
    let duration: f64 = 1.0;
    // Equal-area chain eigenvectors (1,1,1), (1,0,-1), (1,-2,1)
    // have decay rates 0, 1 and 3. This oracle is independent of pair sweeps.
    let slow = 0.25 * (-duration).exp();
    let fast = 0.25 * (-3.0 * duration).exp();
    let exact = [1.0 + slow + fast, 1.0 - 2.0 * fast, 1.0 - slow + fast];
    let mut previous = None;
    for step in [0.2, 0.1, 0.05] {
        let mut amounts = [1.5, 0.5, 1.0];
        diffuse_surface_surfactant(
            ideal(),
            &[1.0; 3],
            &mut amounts,
            &[edge(0, 1, 1.0), edge(1, 2, 1.0)],
            duration,
            controls(step),
            || false,
        )
        .unwrap();
        let error = amounts
            .iter()
            .zip(exact)
            .fold(0.0_f64, |sum, (a, b)| sum.hypot(a - b));
        if let Some(old) = previous {
            assert!((3.7..4.3).contains(&(old / error)));
        }
        assert!(error > 0.0);
        previous = Some(error);
    }
    assert!(previous.unwrap() < 1e-4);
}

#[test]
fn g0_uniform_concentration_and_empty_inventory_remain_exactly_stationary() {
    let areas = [1.0, 2.0, 4.0];
    let edges = [edge(0, 1, 0.3), edge(1, 2, 0.8)];
    for initial in [[0.125, 0.25, 0.5], [0.0; 3]] {
        let mut amounts = initial;
        let report = diffuse_surface_surfactant(
            ideal(),
            &areas,
            &mut amounts,
            &edges,
            2.0,
            controls(0.01),
            || false,
        )
        .unwrap();
        assert_eq!(amounts.map(f64::to_bits), initial.map(f64::to_bits));
        assert_eq!(report.amount_defect_mol, 0.0);
        assert_eq!(report.free_energy_change_j, 0.0);
    }
    let mut amounts = [0.0, 0.5, 1.0];
    let report = diffuse_surface_surfactant(
        ideal(),
        &areas,
        &mut amounts,
        &edges,
        0.0,
        controls(0.1),
        || false,
    )
    .unwrap();
    assert_eq!(amounts, [0.0, 0.5, 1.0]);
    assert_eq!(report.substeps, 0);
}

#[test]
fn g3_connection_orientation_and_area_amount_scaling_preserve_the_diffusion() {
    let areas = [1.0, 2.0, 4.0];
    let initial = [0.7, 0.4, 0.4];
    let mut first = initial;
    let report = diffuse_surface_surfactant(
        ideal(),
        &areas,
        &mut first,
        &[edge(0, 1, 0.3), edge(1, 2, 0.8)],
        0.4,
        controls(0.1),
        || false,
    )
    .unwrap();
    let mut reversed = initial;
    diffuse_surface_surfactant(
        ideal(),
        &areas,
        &mut reversed,
        &[edge(1, 0, 0.3), edge(2, 1, 0.8)],
        0.4,
        controls(0.1),
        || false,
    )
    .unwrap();
    for i in 0..3 {
        close(first[i], reversed[i], 4e-16);
    }
    let scale = 8.0;
    let mut scaled = initial.map(|n| scale * n);
    let scaled_report = diffuse_surface_surfactant(
        ideal(),
        &areas.map(|a| scale * a),
        &mut scaled,
        &[edge(0, 1, scale * 0.3), edge(1, 2, scale * 0.8)],
        0.4,
        controls(0.1),
        || false,
    )
    .unwrap();
    for i in 0..3 {
        close(first[i], scaled[i] / scale, 4e-16);
    }
    close(
        report.free_energy_change_j,
        scaled_report.free_energy_change_j / scale,
        4e-15,
    );
}

#[test]
fn g4_work_source_and_late_cancellation_refusals_preserve_every_amount() {
    let initial = [1.5, 0.5, 1.0];
    let areas = [1.0; 3];
    let edges = [edge(0, 1, 1.0), edge(1, 2, 1.0)];
    let mut amounts = initial;
    let mut budget = controls(0.1);
    budget.maximum_substeps = 1;
    assert!(matches!(
        diffuse_surface_surfactant(ideal(), &areas, &mut amounts, &edges, 1.0, budget, || false,),
        Err(SurfaceDiffusionError::WorkBudget)
    ));
    assert_eq!(amounts, initial);
    for invalid in [
        edge(0, 0, 1.0),
        edge(0, 3, 1.0),
        edge(0, 1, -1.0),
        edge(0, 1, f64::NAN),
    ] {
        assert!(matches!(
            diffuse_surface_surfactant(
                ideal(),
                &areas,
                &mut amounts,
                &[invalid],
                0.1,
                controls(0.1),
                || false,
            ),
            Err(SurfaceDiffusionError::InvalidEdge { .. })
        ));
        assert_eq!(amounts, initial);
    }
    let mut invalid_areas = areas;
    invalid_areas[2] = 0.0;
    assert!(matches!(
        diffuse_surface_surfactant(
            ideal(),
            &invalid_areas,
            &mut amounts,
            &edges,
            0.1,
            controls(0.1),
            || false,
        ),
        Err(SurfaceDiffusionError::Constitutive { face: 2, .. })
    ));
    assert_eq!(amounts, initial);
    let mut polls = 0;
    let result = diffuse_surface_surfactant(
        ideal(),
        &areas,
        &mut amounts,
        &edges,
        0.1,
        controls(0.1),
        || {
            polls += 1;
            polls == 8 // Final publication checkpoint, after candidate evaluation.
        },
    );
    assert!(matches!(result, Err(SurfaceDiffusionError::Cancelled)));
    assert_eq!(polls, 8);
    assert_eq!(amounts.map(f64::to_bits), initial.map(f64::to_bits));
}
