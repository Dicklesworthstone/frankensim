//! G0/G1/G3 binary molecular diffusion: physical analytic and refinement checks.

use super::*;
use crate::conjugate::{
    DensityProjectionControls, FlowFace, FlowField, MaterialAdvectionControls,
    MaterialAdvectionInput, MaterialVolumeSource, Voxel,
};

fn controls(substeps: usize) -> BinarySpeciesDiffusionControls {
    BinarySpeciesDiffusionControls {
        substeps,
        maximum_pair_updates: 1 << 22,
        maximum_storage_bytes: 1 << 22,
        maximum_mass_defect_kg: [1e-11; 2],
        maximum_partition_defect_kg: 1e-11,
    }
}

fn carrier(domain: &VoxelDomain, mass: &[f64]) -> MaterialAdvectionStep {
    let gate = CancelGate::new_clock_free();
    let faces = [FlowFace::Wall; 6];
    let flow = FlowField::from_face_velocity(domain, faces, |_, _| 0.0).unwrap();
    flow.project_and_advect_material(
        domain,
        faces,
        MaterialAdvectionInput {
            mass_kg: mass,
            specific_enthalpy_j_kg: &vec![0.0; mass.len()],
            volume_sources: &vec![MaterialVolumeSource::None; mass.len()],
            boundary_inflow: [None; 6],
        },
        MaterialAdvectionControls {
            projection: DensityProjectionControls {
                duration_s: 1.0,
                relative_tolerance: 1e-12,
                maximum_volume_defect_m3_s: 1e-10,
                maximum_energy_defect_j: 1e-9,
            },
            max_cells: mass.len(),
            max_connections: 6 * mass.len(),
            max_iterations: 100,
            relative_tolerance: 1e-12,
            maximum_mass_defect_kg: 1e-11,
            maximum_enthalpy_defect_j: 1e-9,
        },
        &gate,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}

#[test]
fn graded_unequal_density_pair_matches_exact_fick_relaxation_on_every_axis() {
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![2.0], vec![2.0]];
        widths[axis] = vec![0.3, 1.4];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mass = [1.5, 7.25];
        let state = carrier(&domain, &mass);
        assert_eq!(state.report.mass_routes, 0);
        let initial = [[1.2, 0.3], [1.45, 5.8]];
        let d = [0.03, 0.4];
        let rho = [mass[0] / domain.volume(0), mass[1] / domain.volume(1)];
        let g = domain.face_area(0, axis) / (0.15 / (rho[0] * d[0]) + 0.7 / (rho[1] * d[1]));
        let h = 1.0 / (1.0 / mass[0] + 1.0 / mass[1]);
        for duration in [1e-4, 0.4, 1e6] {
            let q = h
                * (initial[0][0] / mass[0] - initial[1][0] / mass[1])
                * (-(-g * duration / h).exp_m1());
            for substeps in [1, 17] {
                let result = state
                    .diffuse_binary_species(
                        &domain,
                        &initial,
                        &d,
                        duration,
                        controls(substeps),
                        &CancelGate::new_clock_free(),
                    )
                    .unwrap();
                for (c, sign) in [(0, -1.0), (1, 1.0)] {
                    near(
                        result.species_mass_kg[c][0],
                        initial[c][0] + sign * q,
                        3e-14,
                    );
                    near(
                        result.species_mass_kg[c][1],
                        initial[c][1] - sign * q,
                        3e-14,
                    );
                }
                assert_eq!(result.report.diffusive_faces, 1);
                assert_eq!(result.report.pair_updates, 2 * substeps);
                for k in 0..2 {
                    near(
                        result.report.mass_before_kg[k],
                        result.report.mass_after_kg[k],
                        3e-14,
                    );
                }
            }
        }
    }
}

#[test]
fn uniform_mass_fraction_is_stationary_at_nonuniform_density() {
    let domain = VoxelDomain::graded([vec![0.5, 1.0, 2.0], vec![2.0], vec![3.0]]).unwrap();
    let mass = [1.0, 16.0, 0.25];
    let state = carrier(&domain, &mass);
    let initial = mass.map(|m| [0.25 * m, 0.75 * m]);
    let result = state
        .diffuse_binary_species(
            &domain,
            &initial,
            &[0.01, 0.5, 1.0],
            1e6,
            controls(5),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // Diffusing partial density instead of mass fraction would move this field.
    assert_eq!(result.species_mass_kg, initial);
}

#[test]
fn zero_barriers_solids_and_exterior_are_impermeable() {
    let mut domain = VoxelDomain::new(5, 1, 1, 1.0).unwrap();
    domain.set(4, 0, 0, Voxel::Solid(0));
    let mass = [1.0, 1.0, 1.0, 1.0, 0.0];
    let state = carrier(&domain, &mass);
    let initial = [[1.0, 0.0], [0.0, 1.0], [0.2, 0.8], [0.7, 0.3], [0.0; 2]];
    let d = [0.1, 0.1, 0.0, 1.0, 0.0];
    let result = state
        .diffuse_binary_species(
            &domain,
            &initial,
            &d,
            1e6,
            controls(3),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.report.diffusive_faces, 1);
    near(result.species_mass_kg[0][0], 0.5, 1e-15);
    near(result.species_mass_kg[1][0], 0.5, 1e-15);
    assert_eq!(&result.species_mass_kg[2..], &initial[2..]);
    let closed = state
        .diffuse_binary_species(
            &domain,
            &initial,
            &[0.0; 5],
            1.0,
            controls(1_000_000),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(closed.species_mass_kg, initial);
    assert_eq!(closed.report.pair_updates, 0);
}

#[test]
fn trace_species_diffuses_when_majority_fractions_round_to_one() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0, 1.0]);
    let initial = [[1.0, 1e-200], [1.0, 0.0]];
    let mut ctl = controls(8);
    ctl.maximum_mass_defect_kg[1] = 1e-213;
    let result = state
        .diffuse_binary_species(
            &domain,
            &initial,
            &[0.2; 2],
            0.7,
            ctl,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let transfer = 0.5e-200 * (-(-0.28_f64).exp_m1());
    near(result.species_mass_kg[0][1], 1e-200 - transfer, 2e-214);
    near(result.species_mass_kg[1][1], transfer, 2e-214);
    near(result.report.mass_after_kg[1], 1e-200, 2e-214);
    let swapped = initial.map(|[a, b]| [b, a]);
    ctl.maximum_mass_defect_kg = [1e-213, 1e-11];
    let other = state
        .diffuse_binary_species(
            &domain,
            &swapped,
            &[0.2; 2],
            0.7,
            ctl,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    for c in 0..2 {
        assert_eq!(
            other.species_mass_kg[c],
            [result.species_mass_kg[c][1], result.species_mass_kg[c][0],]
        );
    }
}

fn mode_error(cells: usize, substeps: usize, continuum: bool) -> f64 {
    let h = 1.0 / cells as f64;
    let domain = VoxelDomain::graded([vec![h; cells], vec![1.0], vec![1.0]]).unwrap();
    let mass = vec![h; cells];
    let state = carrier(&domain, &mass);
    let mode = |c: usize| (core::f64::consts::PI * (c as f64 + 0.5) * h).cos();
    let initial: Vec<_> = (0..cells)
        .map(|c| {
            let a = 0.5 + 0.2 * mode(c);
            [h * a, h * (1.0 - a)]
        })
        .collect();
    let result = state
        .diffuse_binary_species(
            &domain,
            &initial,
            &vec![0.1; cells],
            0.2,
            controls(substeps),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let lambda = if continuum {
        0.1 * core::f64::consts::PI.powi(2)
    } else {
        0.4 / h.powi(2) * (0.5 * core::f64::consts::PI * h).sin().powi(2)
    };
    (0..cells)
        .map(|c| {
            (result.species_mass_kg[c][0] / h - (0.5 + 0.2 * (-lambda * 0.2).exp() * mode(c))).abs()
        })
        .fold(0.0, f64::max)
}

#[test]
fn symmetric_sweeps_converge_quadratically_to_the_discrete_neumann_mode() {
    let errors = [8, 16, 32].map(|steps| mode_error(8, steps, false));
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
            "{errors:?}"
        );
    }
}

#[test]
fn spatial_refinement_recovers_continuum_diffusion_length_and_decay() {
    let errors = [8, 16, 32].map(|cells| mode_error(cells, 2 * cells * cells, true));
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
            "{errors:?}"
        );
    }
}

#[test]
fn strong_multiface_steps_preserve_positivity_partition_and_each_species() {
    let domain = VoxelDomain::new(8, 1, 1, 1.0).unwrap();
    let mass = [1.0, 2.0, 3.0, 4.0, 0.5, 0.25, 8.0, 16.0];
    let state = carrier(&domain, &mass);
    let initial: Vec<_> = mass
        .iter()
        .enumerate()
        .map(|(c, &m)| if c % 2 == 0 { [m, 0.0] } else { [0.0, m] })
        .collect();
    let result = state
        .diffuse_binary_species(
            &domain,
            &initial,
            &[100.0; 8],
            1e6,
            controls(3),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    for (m, species) in mass.iter().zip(&result.species_mass_kg) {
        assert!(species.iter().all(|&x| x >= 0.0 && x <= *m));
        near(species[0] + species[1], *m, 5e-15);
    }
    for k in 0..2 {
        near(
            result.report.mass_after_kg[k],
            result.report.mass_before_kg[k],
            1e-14,
        );
    }
    assert!(result.report.maximum_pair_relaxation > 0.999);
}

#[test]
fn invalid_data_budgets_and_changed_carrier_refuse_without_publishing() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let initial = [[1.0, 0.0], [0.0, 1.0]];
    let gate = CancelGate::new_clock_free();
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            state
                .diffuse_binary_species(&domain, &initial, &[bad, 0.1], 1.0, controls(1), &gate,)
                .is_err()
        );
    }
    for ctl in [
        BinarySpeciesDiffusionControls {
            substeps: 0,
            ..controls(1)
        },
        BinarySpeciesDiffusionControls {
            maximum_storage_bytes: 1,
            ..controls(1)
        },
        BinarySpeciesDiffusionControls {
            maximum_pair_updates: 1,
            ..controls(1)
        },
    ] {
        assert!(
            state
                .diffuse_binary_species(&domain, &initial, &[0.1; 2], 1.0, ctl, &gate,)
                .is_err()
        );
    }
    assert!(
        state
            .diffuse_binary_species(
                &domain,
                &[[1.0, 0.1], [0.0, 1.0]],
                &[0.1; 2],
                1.0,
                controls(1),
                &gate,
            )
            .is_err()
    );
    let mut changed = state.clone();
    changed.mass_kg[0] = 2.0;
    assert!(
        changed
            .diffuse_binary_species(&domain, &initial, &[0.1; 2], 1.0, controls(1), &gate,)
            .is_err()
    );
    let mut changed_domain = domain.clone();
    changed_domain.set(1, 0, 0, Voxel::Solid(0));
    assert!(
        state
            .diffuse_binary_species(
                &changed_domain,
                &initial,
                &[0.1; 2],
                1.0,
                controls(1),
                &gate,
            )
            .is_err()
    );
    let original = state.clone();
    let error = state
        .diffuse_species_checked(
            &domain,
            &initial,
            &[0.1; 2],
            1.0,
            controls(4),
            &gate,
            |substep| {
                if substep == 2 {
                    Err(ChtError::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .unwrap_err();
    assert_eq!(error, ChtError::Cancelled);
    assert_eq!(state, original);
    assert_eq!(initial, [[1.0, 0.0], [0.0, 1.0]]);
}

fn common_controls<const N: usize>(substeps: usize) -> SpeciesDiffusionControls<N> {
    SpeciesDiffusionControls {
        substeps,
        maximum_pair_updates: 1 << 22,
        maximum_storage_bytes: 1 << 22,
        maximum_mass_defect_kg: [1e-11; N],
        maximum_partition_defect_kg: 1e-11,
    }
}

#[test]
fn four_species_share_the_exact_graded_variable_density_pair_relaxation() {
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![2.0], vec![2.0]];
        widths[axis] = vec![0.3, 1.4];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mass = [1.5, 7.25];
        let state = carrier(&domain, &mass);
        let initial = [[0.75, 0.45, 0.3, 0.0], [0.725, 1.45, 2.175, 2.9]];
        let d = [0.03, 0.4];
        let rho = [mass[0] / domain.volume(0), mass[1] / domain.volume(1)];
        let conductance =
            domain.face_area(0, axis) / (0.15 / (rho[0] * d[0]) + 0.7 / (rho[1] * d[1]));
        let harmonic_mass = 1.0 / (1.0 / mass[0] + 1.0 / mass[1]);
        for duration in [1e-4, 0.4, 1e6] {
            let relaxation = -(-conductance * duration / harmonic_mass).exp_m1();
            for substeps in [1, 17] {
                let result = state
                    .diffuse_species(
                        &domain,
                        &initial,
                        &d,
                        duration,
                        common_controls(substeps),
                        &CancelGate::new_clock_free(),
                    )
                    .unwrap();
                for k in 0..4 {
                    let q = harmonic_mass
                        * relaxation
                        * (initial[0][k] / mass[0] - initial[1][k] / mass[1]);
                    for (c, sign) in [(0, -1.0), (1, 1.0)] {
                        near(
                            result.species_mass_kg[c][k],
                            initial[c][k] + sign * q,
                            3e-14,
                        );
                        assert!(result.species_mass_kg[c][k] >= 0.0);
                    }
                    near(
                        result.report.mass_before_kg[k],
                        result.report.mass_after_kg[k],
                        3e-14,
                    );
                }
                assert_eq!(result.report.diffusive_faces, 1);
                assert_eq!(result.report.pair_updates, 2 * substeps);
                assert_eq!(result.report.component_updates, 8 * substeps);
                for (c, &m) in mass.iter().enumerate() {
                    near(result.species_mass_kg[c].iter().sum(), m, 3e-14);
                }
            }
        }
    }
}

#[test]
fn common_mixing_preserves_trace_inventory_and_species_permutation() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0, 2.0, 4.0]);
    let initial = [
        [1.0, 1e-200, 0.0, 0.0],
        [0.4, 0.0, 1.6, 0.0],
        [0.0, 0.0, 0.0, 4.0],
    ];
    let d = [0.3, 0.4, 0.2];
    let gate = CancelGate::new_clock_free();
    let mut ctl = common_controls(7);
    ctl.maximum_mass_defect_kg[1] = 1e-213;
    let result = state
        .diffuse_species(&domain, &initial, &d, 0.7, ctl, &gate)
        .unwrap();
    near(result.report.mass_after_kg[1], 1e-200, 1e-213);
    assert!(result.species_mass_kg.iter().all(|m| m[1] > 0.0));
    let permutation = [2, 0, 3, 1];
    let permuted = initial.map(|m| permutation.map(|k| m[k]));
    ctl.maximum_mass_defect_kg = permutation.map(|k| ctl.maximum_mass_defect_kg[k]);
    let other = state
        .diffuse_species(&domain, &permuted, &d, 0.7, ctl, &gate)
        .unwrap();
    for c in 0..3 {
        assert_eq!(
            other.species_mass_kg[c],
            permutation.map(|k| result.species_mass_kg[c][k])
        );
    }
    // Refining the model's viscosity or Schmidt changes only the shared rate.
    let twice_d = d.map(|value| 2.0 * value);
    let same_rate_time = state
        .diffuse_species(&domain, &initial, &twice_d, 0.35, common_controls(7), &gate)
        .unwrap();
    assert_eq!(same_rate_time.species_mass_kg, result.species_mass_kg);
}

#[test]
fn common_mixing_has_exact_uniform_zero_and_single_species_controls() {
    let domain = VoxelDomain::graded([vec![0.5, 1.0, 2.0], vec![2.0], vec![3.0]]).unwrap();
    let mass = [1.0, 16.0, 0.25];
    let state = carrier(&domain, &mass);
    let initial = mass.map(|m| [0.125 * m, 0.25 * m, 0.5 * m, 0.125 * m]);
    let gate = CancelGate::new_clock_free();
    let uniform = state
        .diffuse_species(
            &domain,
            &initial,
            &[0.01, 0.5, 1.0],
            1e6,
            common_controls(5),
            &gate,
        )
        .unwrap();
    assert_eq!(uniform.species_mass_kg, initial);
    let distinct = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 8.0, 8.0, 0.0],
        [0.0, 0.0, 0.0, 0.25],
    ];
    let disabled = state
        .diffuse_species(
            &domain,
            &distinct,
            &[0.0; 3],
            1e6,
            common_controls(1_000_000),
            &gate,
        )
        .unwrap();
    assert_eq!(disabled.species_mass_kg, distinct);
    assert_eq!(disabled.report.pair_updates, 0);
    assert_eq!(disabled.report.component_updates, 0);
    let single = mass.map(|m| [m]);
    let pure = state
        .diffuse_species(&domain, &single, &[0.2; 3], 1e6, common_controls(3), &gate)
        .unwrap();
    assert_eq!(pure.species_mass_kg, single);
}

#[test]
fn all_common_species_recover_the_neumann_diffusion_mode_with_second_order_sweeps() {
    let cells = 8;
    let h = 1.0 / cells as f64;
    let domain = VoxelDomain::graded([vec![h; cells], vec![1.0], vec![1.0]]).unwrap();
    let state = carrier(&domain, &vec![h; cells]);
    let mode = |c: usize| (core::f64::consts::PI * (c as f64 + 0.5) * h).cos();
    let mean = [0.2, 0.3, 0.1, 0.4];
    let amplitude = [0.05, -0.1, 0.025, 0.025];
    let initial: Vec<[f64; 4]> = (0..cells)
        .map(|c| core::array::from_fn(|k| h * (mean[k] + amplitude[k] * mode(c))))
        .collect();
    let decay = (-0.4 / h.powi(2) * (0.5 * core::f64::consts::PI * h).sin().powi(2) * 0.2).exp();
    let errors = [8, 16, 32].map(|substeps| {
        let result = state
            .diffuse_species(
                &domain,
                &initial,
                &vec![0.1; cells],
                0.2,
                common_controls(substeps),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let mut error = [0.0_f64; 4];
        for c in 0..cells {
            for k in 0..4 {
                let exact = mean[k] + amplitude[k] * decay * mode(c);
                error[k] = error[k].max((result.species_mass_kg[c][k] / h - exact).abs());
            }
        }
        error
    });
    for k in 0..4 {
        for pair in errors.windows(2) {
            let order_ratio = pair[0][k] / pair[1][k];
            assert!(
                (3.8..4.2).contains(&order_ratio),
                "component {k}: {errors:?}"
            );
        }
    }
}

#[test]
fn common_mixing_admits_component_bounds_and_charges_component_storage() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let gate = CancelGate::new_clock_free();
    assert!(
        state
            .diffuse_species(&domain, &[[]; 2], &[0.1; 2], 1.0, common_controls(1), &gate,)
            .is_err()
    );
    assert!(
        state
            .diffuse_species(
                &domain,
                &[[1.0 / 129.0; 129]; 2],
                &[0.1; 2],
                1.0,
                common_controls(1),
                &gate,
            )
            .is_err()
    );
    let initial = [[1.0 / MAX_MATERIAL_SPECIES as f64; MAX_MATERIAL_SPECIES]; 2];
    let result = state
        .diffuse_species(&domain, &initial, &[0.1; 2], 1.0, common_controls(1), &gate)
        .unwrap();
    assert_eq!(result.species_mass_kg, initial);
    assert_eq!(result.report.component_updates, 2 * MAX_MATERIAL_SPECIES);
    let mut small = common_controls(1);
    small.maximum_storage_bytes = 4096;
    assert!(
        state
            .diffuse_species(&domain, &initial, &[0.1; 2], 1.0, small, &gate)
            .is_err()
    );
    let mut binary = controls(1);
    binary.maximum_storage_bytes = 4096;
    assert!(
        state
            .diffuse_binary_species(
                &domain,
                &[[1.0, 0.0], [0.0, 1.0]],
                &[0.1; 2],
                1.0,
                binary,
                &gate,
            )
            .is_ok()
    );

    let initial = [[0.25, 0.0, 0.25, 0.5], [0.0, 1.0, 0.0, 0.0]];
    let transposed = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    assert!(
        state
            .diffuse_species(
                &transposed,
                &initial,
                &[0.1; 2],
                1.0,
                common_controls(1),
                &gate,
            )
            .is_err()
    );
    let before = state.clone();
    let error = state
        .diffuse_species_checked(
            &domain,
            &initial,
            &[0.1; 2],
            1.0,
            common_controls(4),
            &gate,
            |substep| {
                if substep == 2 {
                    Err(ChtError::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .unwrap_err();
    assert_eq!(error, ChtError::Cancelled);
    assert_eq!(state, before);
}

fn adaptive_controls<const N: usize>(tolerance: [f64; N]) -> AdaptiveSpeciesDiffusionControls<N> {
    AdaptiveSpeciesDiffusionControls {
        solver: common_controls(1),
        absolute_mass_fraction_tolerance: tolerance,
        relative_tolerance: 0.0,
        maximum_refinements: 20,
    }
}

#[test]
fn adaptive_diffusion_resolves_a_stable_but_inaccurate_multiface_step() {
    let cells = 8;
    let h = 1.0 / cells as f64;
    let duration = 2.0;
    let d = 0.1;
    let domain = VoxelDomain::graded([vec![h; cells], vec![1.0], vec![1.0]]).unwrap();
    let state = carrier(&domain, &vec![h; cells]);
    let initial: Vec<_> = (0..cells).map(|c| {
        let y = 0.5 + 0.4 * (core::f64::consts::PI * (c as f64 + 0.5) * h).cos();
        [h * y, h * (1.0 - y)]
    }).collect();
    // Independent eigenmode of the cell-centred Neumann diffusion matrix.
    let decay = (-4.0 * d / h.powi(2)
        * (0.5 * core::f64::consts::PI * h).sin().powi(2) * duration).exp();
    let error = |values: &[[f64; 2]]| -> f64 {
        values.iter().zip(&initial).map(|(actual, original)| {
            (actual[0] / h - (0.5 + (original[0] / h - 0.5) * decay)).abs()
        }).fold(0.0, f64::max)
    };
    let gate = CancelGate::new_clock_free();
    let unresolved = state.diffuse_binary_species(
        &domain, &initial, &vec![d; cells], duration, controls(1), &gate,
    ).unwrap();
    assert!(error(&unresolved.species_mass_kg) > 0.2);
    assert!(unresolved.species_mass_kg.iter().flatten().all(|&m| m >= 0.0));
    let result = state.diffuse_species_adaptive(
        &domain, &initial, &vec![d; cells], duration, adaptive_controls([2e-7; 2]), &gate,
    ).unwrap();
    assert!(result.refinements > 1);
    assert!(result.maximum_error_ratio <= 1.0);
    assert!(error(&result.step.species_mass_kg) < 3e-7);
    let fine_steps = result.step.report.substeps;
    // Reuse the previous fine solution as coarse: no repeated coarse solve.
    assert_eq!(result.total_pair_updates, 2 * (cells - 1) * (2 * fine_steps - 1));
    assert_eq!(result.total_component_updates, 2 * result.total_pair_updates);
    for m in &result.step.species_mass_kg {
        assert!(m.iter().all(|&value| value >= 0.0 && value <= h));
        near(m[0] + m[1], h, 2e-15);
    }
    for k in 0..2 {
        near(result.step.report.mass_before_kg[k], result.step.report.mass_after_kg[k], 1e-14);
    }
}

#[test]
fn adaptive_diffusion_resolves_trace_species_with_its_own_tolerance() {
    let cells = 6;
    let h = 1.0 / cells as f64;
    let domain = VoxelDomain::graded([vec![h; cells], vec![1.0], vec![1.0]]).unwrap();
    let state = carrier(&domain, &vec![h; cells]);
    let mode = |c: usize| (core::f64::consts::PI * (c as f64 + 0.5) * h).cos();
    let initial: Vec<_> = (0..cells).map(|c| {
        [0.6 * h, 0.4 * h, h * 1e-200 * (1.0 + 0.4 * mode(c))]
    }).collect();
    let mut ctl = adaptive_controls([1e-3, 1e-3, 1e-206]);
    ctl.solver.maximum_mass_defect_kg[2] = 1e-212;
    let gate = CancelGate::new_clock_free();
    let result = state.diffuse_species_adaptive(
        &domain, &initial, &vec![0.1; cells], 0.5, ctl, &gate,
    ).unwrap();
    assert!(result.refinements > 1);
    let decay = (-0.4 / h.powi(2)
        * (0.5 * core::f64::consts::PI * h).sin().powi(2) * 0.5).exp();
    for c in 0..cells {
        near(result.step.species_mass_kg[c][2] / h,
            1e-200 * (1.0 + 0.4 * decay * mode(c)), 2e-206);
    }
    near(result.step.report.mass_after_kg[2], 1e-200, 1e-212);
    let permutation = [2, 0, 1];
    let permuted: Vec<_> = initial.iter().map(|m| permutation.map(|k| m[k])).collect();
    ctl.absolute_mass_fraction_tolerance = permutation.map(|k| ctl.absolute_mass_fraction_tolerance[k]);
    ctl.solver.maximum_mass_defect_kg = permutation.map(|k| ctl.solver.maximum_mass_defect_kg[k]);
    let other = state.diffuse_species_adaptive(
        &domain, &permuted, &vec![0.1; cells], 0.5, ctl, &gate,
    ).unwrap();
    assert_eq!(other.refinements, result.refinements);
    for c in 0..cells {
        assert_eq!(other.step.species_mass_kg[c], permutation.map(|k| result.step.species_mass_kg[c][k]));
    }
}

#[test]
fn adaptive_diffusion_budgets_include_discarded_solves_and_live_comparisons() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let initial = [[1.0, 0.0], [0.0, 1.0]];
    let gate = CancelGate::new_clock_free();
    let mut ctl = adaptive_controls([1e-6; 2]);
    ctl.solver.maximum_pair_updates = 5;
    // The fine solve needs only four pair updates, but its coarse comparison
    // has already consumed two. A five-update aggregate budget must refuse.
    assert!(state.diffuse_species_adaptive(
        &domain, &initial, &[0.1; 2], 1.0, ctl, &gate,
    ).is_err());
    ctl.solver.maximum_pair_updates = 6;
    let accepted = state.diffuse_species_adaptive(
        &domain, &initial, &[0.1; 2], 1.0, ctl, &gate,
    ).unwrap();
    assert_eq!(accepted.refinements, 1);
    assert_eq!(accepted.total_pair_updates, 6);
    assert_eq!(accepted.total_component_updates, 12);
    ctl.solver.maximum_storage_bytes = 3000;
    assert!(state.diffuse_binary_species(
        &domain, &initial, &[0.1; 2], 1.0, ctl.solver, &gate,
    ).is_ok());
    assert!(state.diffuse_species_adaptive(
        &domain, &initial, &[0.1; 2], 1.0, ctl, &gate,
    ).is_err());
}

#[test]
fn adaptive_diffusion_refinement_refusal_and_cancellation_publish_nothing() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 3]);
    let before = state.clone();
    let initial = [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]];
    let gate = CancelGate::new_clock_free();
    let mut ctl = adaptive_controls([1e-12; 2]);
    ctl.maximum_refinements = 1;
    assert!(state.diffuse_species_adaptive(
        &domain, &initial, &[0.3; 3], 1.0, ctl, &gate,
    ).is_err());
    ctl.maximum_refinements = 8;
    let error = state.diffuse_species_adaptive_checked(
        &domain, &initial, &[0.3; 3], 1.0, ctl, &gate,
        |refinement| if refinement == 2 { Err(ChtError::Cancelled) } else { Ok(()) },
    ).unwrap_err();
    assert_eq!(error, ChtError::Cancelled);
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        ctl.absolute_mass_fraction_tolerance[1] = invalid;
        assert!(state.diffuse_species_adaptive(
            &domain, &initial, &[0.3; 3], 1.0, ctl, &gate,
        ).is_err());
    }
    assert_eq!(state, before);
    assert_eq!(initial, [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]]);
}

#[test]
fn adaptive_diffusion_preserves_the_zero_coefficient_identity_without_refinement() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 3]);
    let initial = [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]];
    let mut ctl = adaptive_controls([1e-12; 2]);
    ctl.solver.substeps = 1_000_000;
    ctl.solver.maximum_pair_updates = 0;
    let result = state.diffuse_species_adaptive(
        &domain, &initial, &[0.0; 3], 1.0, ctl, &CancelGate::new_clock_free(),
    ).unwrap();
    assert_eq!(result.step.species_mass_kg, initial);
    assert_eq!(result.refinements, 0);
    assert_eq!(result.total_pair_updates, 0);
    assert_eq!(result.estimated_mass_fraction_error, [0.0; 2]);
}

#[test]
fn partial_enthalpy_diffuses_with_the_exact_binary_mass_transfer() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let initial = [[0.8, 0.2], [0.2, 0.8]];
    let partial_h = [-200_000.0, 500_000.0];
    let initial_h = initial.map(|m| m[0] * partial_h[0] + m[1] * partial_h[1]);
    let calls = std::cell::Cell::new(0);
    let result = state.diffuse_species_with_enthalpy(
        &domain, &initial, &initial_h, &[0.15; 2], 2.0,
        |cells, axis| {
            assert_eq!(cells, [0, 1]);
            assert_eq!(axis, 0);
            calls.set(calls.get() + 1);
            Ok(partial_h)
        }, controls(7), 1e-8, &gate,
    ).unwrap();
    assert_eq!(calls.get(), 1);
    let mass_transfer = 0.5 * 0.6 * (-(-0.6_f64).exp_m1());
    let energy_transfer = mass_transfer * (partial_h[0] - partial_h[1]);
    near(result.diffusive_enthalpy_change_j[0], -energy_transfer, 1e-9);
    near(result.diffusive_enthalpy_change_j[1], energy_transfer, 1e-9);
    // Isothermal ideal mixing must retain its original caloric temperature:
    // endpoint mixture enthalpy equals the transported caloric inventory.
    for c in 0..2 {
        let m = result.species.species_mass_kg[c];
        near(result.specific_enthalpy_j_kg[c],
            m[0] * partial_h[0] + m[1] * partial_h[1], 2e-9);
    }
    let mass_only = state.diffuse_binary_species(
        &domain, &initial, &[0.15; 2], 2.0, controls(7), &gate,
    ).unwrap();
    assert_eq!(result.species, mass_only);
    assert!(result.report.absolute_enthalpy_defect_j < 1e-8);
}

#[test]
fn multicomponent_enthalpy_mixing_preserves_isothermal_caloric_inventory() {
    let domain = VoxelDomain::graded([vec![0.5, 1.0, 2.0], vec![1.0], vec![1.0]]).unwrap();
    let mass = [0.5, 2.0, 1.0];
    let state = carrier(&domain, &mass);
    let initial = [[0.35, 0.1, 0.05], [0.2, 1.4, 0.4], [0.3, 0.2, 0.5]];
    let h = [-1e6, 3e5, 8e5];
    let initial_h: Vec<_> = initial.iter().zip(mass).map(|(m, carrier)|
        m.iter().zip(h).map(|(mk, hk)| mk * hk).sum::<f64>() / carrier
    ).collect();
    let ctl = SpeciesDiffusionControls {
        substeps: 16, maximum_pair_updates: 1000, maximum_storage_bytes: 1 << 20,
        maximum_mass_defect_kg: [1e-12; 3], maximum_partition_defect_kg: 1e-12,
    };
    let gate = CancelGate::new_clock_free();
    let result = state.diffuse_species_with_enthalpy(
        &domain, &initial, &initial_h, &[0.1; 3], 4.0, |_, _| Ok(h),
        ctl, 1e-7, &gate,
    ).unwrap();
    for c in 0..3 {
        let expected = result.species.species_mass_kg[c].iter().zip(h)
            .map(|(mk, hk)| mk * hk).sum::<f64>();
        near(mass[c] * result.specific_enthalpy_j_kg[c], expected, 1e-8);
    }
    // A common source datum must not turn rounding of sum(j_k) into heat.
    let uniform = state.diffuse_species_with_enthalpy(
        &domain, &initial, &[1e200; 3], &[0.1; 3], 4.0,
        |_, _| Ok([1e200; 3]), ctl, 0.0, &gate,
    ).unwrap();
    assert_eq!(uniform.diffusive_enthalpy_change_j, [0.0; 3]);
    assert_eq!(uniform.specific_enthalpy_j_kg, [1e200; 3]);
    assert_eq!(uniform.report.absolute_enthalpy_defect_j, 0.0);
}

#[test]
fn diffusive_enthalpy_is_invariant_to_a_common_resolvable_datum() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let initial = [[0.75, 0.25], [0.25, 0.75]];
    let run = |datum: f64| state.diffuse_species_with_enthalpy(
        &domain, &initial, &[datum + 125.0, datum + 175.0], &[0.2; 2], 0.5,
        |_, _| Ok([datum + 100.0, datum + 200.0]), controls(3), 0.001,
        &CancelGate::new_clock_free(),
    ).unwrap();
    let native = run(0.0);
    let shifted = run(1e12);
    assert_eq!(native.species, shifted.species);
    assert_eq!(native.diffusive_enthalpy_change_j, shifted.diffusive_enthalpy_change_j);
    for c in 0..2 {
        near(shifted.specific_enthalpy_j_kg[c] - 1e12, native.specific_enthalpy_j_kg[c], 0.0002);
    }
}

#[test]
fn enthalpy_diffusion_refuses_unrepresented_heat_and_failed_face_sources_atomically() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let initial = [[0.75, 0.25], [0.25, 0.75]];
    let gate = CancelGate::new_clock_free();
    let run = |initial_h, limit, ctl, partial_h| state.diffuse_species_with_enthalpy(
        &domain, &initial, initial_h, &[0.2; 2], 0.5, |_, _| Ok(partial_h),
        ctl, limit, &gate,
    );
    assert!(run(&[1e30; 2], 1e-12, controls(2), [100.0, 200.0]).is_err());
    assert!(run(&[125.0, 175.0], 1e-9, controls(2), [f64::NAN, 200.0]).is_err());
    let mut constrained = controls(2);
    constrained.maximum_storage_bytes = 1;
    assert!(run(&[125.0, 175.0], 1e-9, constrained, [100.0, 200.0]).is_err());
    let source_error = state.diffuse_species_with_enthalpy(
        &domain, &initial, &[125.0, 175.0], &[0.2; 2], 0.5,
        |_, _| Err(refuse("original source is out of its thermal domain")),
        controls(2), 1e-9, &gate,
    );
    assert!(source_error.is_err());
    let cancelled = state.diffuse_species_with_enthalpy(
        &domain, &initial, &[125.0, 175.0], &[0.2; 2], 0.5,
        |_, _| { gate.request(); Ok([100.0, 200.0]) },
        controls(2), 1e-9, &gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(state.mass_kg, [1.0; 2]);
    assert_eq!(initial, [[0.75, 0.25], [0.25, 0.75]]);
}

#[test]
fn rounded_carrier_flux_cannot_create_unbudgeted_chemical_heat() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    // The majority cannot represent the trace correction at this scale. Its
    // independent mass budget is permissive, but the chemical energy budget
    // must still see the amplified mass-frame defect.
    let initial = [[0.0, 0.5, 1e-20, 0.5], [0.0, 0.25, 0.0, 0.75]];
    let ctl = SpeciesDiffusionControls {
        substeps: 1, maximum_pair_updates: 10, maximum_storage_bytes: 1 << 20,
        maximum_mass_defect_kg: [1e-12, 1e-12, 1e-30, 1e-12],
        maximum_partition_defect_kg: 1e-12,
    };
    let result = state.diffuse_species_with_enthalpy(
        &domain, &initial, &[0.0; 2], &[0.2; 2], 0.5,
        |_, _| Ok([0.0, 1e20, 1e20, 1e20]), ctl, 1e-10,
        &CancelGate::new_clock_free(),
    );
    assert!(result.is_err());
}

#[test]
fn opposing_species_formation_terms_retain_small_diffusive_heat() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let state = carrier(&domain, &[1.0; 2]);
    let initial = [[0.5, 0.375, 0.125], [0.5, 0.125, 0.375]];
    // In both charts h_mix = 4*Y_2: Y_0 stays at one half, so the
    // negative and positive formation terms cancel in each original cell.
    let initial_h = [0.5, 1.5];
    let gate = CancelGate::new_clock_free();
    let run = |h| state.diffuse_species_with_enthalpy(
        &domain, &initial, &initial_h, &[0.2; 2], 0.5, |_, _| Ok(h),
        common_controls(1), 1e-13, &gate,
    ).unwrap();
    let baseline = run([0.0, 0.0, 4.0]);
    let transfer = 0.5 * 0.25 * (-(-0.2_f64).exp_m1());
    near(baseline.diffusive_enthalpy_change_j[0], 4.0 * transfer, 1e-14);
    for formation in [2.0_f64.powi(53), 2.0_f64.powi(54)] {
        // The first chart loses product low parts; the second also loses
        // the +4 while centering h_2-h_0. Either loss changes actual heat,
        // even though a mass or paired-energy conservation audit can pass.
        let result = run([-formation, formation, formation + 4.0]);
        assert_eq!(result.species, baseline.species);
        for c in 0..2 {
            near(result.diffusive_enthalpy_change_j[c],
                baseline.diffusive_enthalpy_change_j[c], 1e-14);
            near(result.specific_enthalpy_j_kg[c],
                4.0 * result.species.species_mass_kg[c][2], 1e-14);
        }
        assert!(result.report.absolute_enthalpy_defect_j < 1e-13);
    }
}

#[test]
fn mass_product_roundoff_is_charged_to_the_enthalpy_budget() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let epsilon = 2.0_f64.powi(-27);
    let mass = 1.0 + epsilon;
    let state = carrier(&domain, &[mass; 2]);
    let initial = [[1.0, epsilon], [0.0, mass]];
    let gate = CancelGate::new_clock_free();
    let run = |limit| state.diffuse_species_with_enthalpy(
        &domain, &initial, &[0.0; 2], &[1e30; 2], 1.0,
        |_, _| Ok([0.0, 2.0_f64.powi(55)]), controls(1), limit, &gate,
    );
    // Full pair relaxation transfers exactly half a kilogram, hence 2^54 J.
    // Each returned |h| is 2^54-2^27; its EXACT represented mass product is
    // 2^54-1, although ordinary f64 multiplication rounds it back to 2^54.
    // The two one-joule local errors must not hide behind zero global change.
    assert!(run(1.0).is_err());
    let accepted = run(2.0).unwrap();
    assert_eq!(accepted.diffusive_enthalpy_change_j,
        [2.0_f64.powi(54), -2.0_f64.powi(54)]);
    assert_eq!(accepted.report.stored_enthalpy_change_j, 0.0);
    assert_eq!(accepted.report.transported_enthalpy_change_j, 0.0);
    assert_eq!(accepted.report.absolute_enthalpy_defect_j, 2.0);
    assert_eq!(state.mass_kg, [mass; 2]);
}
