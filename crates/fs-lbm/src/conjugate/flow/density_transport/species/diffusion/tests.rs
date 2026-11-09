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
        let g = domain.face_area(0, axis)
            / (0.15 / (rho[0] * d[0]) + 0.7 / (rho[1] * d[1]));
        let h = 1.0 / (1.0 / mass[0] + 1.0 / mass[1]);
        for duration in [1e-4, 0.4, 1e6] {
            let q = h * (initial[0][0] / mass[0] - initial[1][0] / mass[1])
                * (-(-g * duration / h).exp_m1());
            for substeps in [1, 17] {
                let result = state.diffuse_binary_species(
                    &domain, &initial, &d, duration, controls(substeps),
                    &CancelGate::new_clock_free(),
                ).unwrap();
                for (c, sign) in [(0, -1.0), (1, 1.0)] {
                    near(result.species_mass_kg[c][0], initial[c][0] + sign * q, 3e-14);
                    near(result.species_mass_kg[c][1], initial[c][1] - sign * q, 3e-14);
                }
                assert_eq!(result.report.diffusive_faces, 1);
                assert_eq!(result.report.pair_updates, 2 * substeps);
                for k in 0..2 {
                    near(result.report.mass_before_kg[k], result.report.mass_after_kg[k], 3e-14);
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
    let result = state.diffuse_binary_species(
        &domain, &initial, &[0.01, 0.5, 1.0], 1e6, controls(5),
        &CancelGate::new_clock_free(),
    ).unwrap();
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
    let result = state.diffuse_binary_species(
        &domain, &initial, &d, 1e6, controls(3), &CancelGate::new_clock_free(),
    ).unwrap();
    assert_eq!(result.report.diffusive_faces, 1);
    near(result.species_mass_kg[0][0], 0.5, 1e-15);
    near(result.species_mass_kg[1][0], 0.5, 1e-15);
    assert_eq!(&result.species_mass_kg[2..], &initial[2..]);
    let closed = state.diffuse_binary_species(
        &domain, &initial, &[0.0; 5], 1.0, controls(1_000_000),
        &CancelGate::new_clock_free(),
    ).unwrap();
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
    let result = state.diffuse_binary_species(
        &domain, &initial, &[0.2; 2], 0.7, ctl, &CancelGate::new_clock_free(),
    ).unwrap();
    let transfer = 0.5e-200 * (-(-0.28_f64).exp_m1());
    near(result.species_mass_kg[0][1], 1e-200 - transfer, 2e-214);
    near(result.species_mass_kg[1][1], transfer, 2e-214);
    near(result.report.mass_after_kg[1], 1e-200, 2e-214);
    let swapped = initial.map(|[a, b]| [b, a]);
    ctl.maximum_mass_defect_kg = [1e-213, 1e-11];
    let other = state.diffuse_binary_species(
        &domain, &swapped, &[0.2; 2], 0.7, ctl, &CancelGate::new_clock_free(),
    ).unwrap();
    for c in 0..2 {
        assert_eq!(other.species_mass_kg[c], [
            result.species_mass_kg[c][1], result.species_mass_kg[c][0],
        ]);
    }
}

fn mode_error(cells: usize, substeps: usize, continuum: bool) -> f64 {
    let h = 1.0 / cells as f64;
    let domain = VoxelDomain::graded([vec![h; cells], vec![1.0], vec![1.0]]).unwrap();
    let mass = vec![h; cells];
    let state = carrier(&domain, &mass);
    let mode = |c: usize| (core::f64::consts::PI * (c as f64 + 0.5) * h).cos();
    let initial: Vec<_> = (0..cells).map(|c| {
        let a = 0.5 + 0.2 * mode(c);
        [h * a, h * (1.0 - a)]
    }).collect();
    let result = state.diffuse_binary_species(
        &domain, &initial, &vec![0.1; cells], 0.2, controls(substeps),
        &CancelGate::new_clock_free(),
    ).unwrap();
    let lambda = if continuum {
        0.1 * core::f64::consts::PI.powi(2)
    } else {
        0.4 / h.powi(2) * (0.5 * core::f64::consts::PI * h).sin().powi(2)
    };
    (0..cells).map(|c| {
        (result.species_mass_kg[c][0] / h - (0.5 + 0.2 * (-lambda * 0.2).exp() * mode(c))).abs()
    }).fold(0.0, f64::max)
}

#[test]
fn symmetric_sweeps_converge_quadratically_to_the_discrete_neumann_mode() {
    let errors = [8, 16, 32].map(|steps| mode_error(8, steps, false));
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2, "{errors:?}");
    }
}

#[test]
fn spatial_refinement_recovers_continuum_diffusion_length_and_decay() {
    let errors = [8, 16, 32].map(|cells| mode_error(cells, 2 * cells * cells, true));
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2, "{errors:?}");
    }
}

#[test]
fn strong_multiface_steps_preserve_positivity_partition_and_each_species() {
    let domain = VoxelDomain::new(8, 1, 1, 1.0).unwrap();
    let mass = [1.0, 2.0, 3.0, 4.0, 0.5, 0.25, 8.0, 16.0];
    let state = carrier(&domain, &mass);
    let initial: Vec<_> = mass.iter().enumerate().map(|(c, &m)| {
        if c % 2 == 0 { [m, 0.0] } else { [0.0, m] }
    }).collect();
    let result = state.diffuse_binary_species(
        &domain, &initial, &[100.0; 8], 1e6, controls(3), &CancelGate::new_clock_free(),
    ).unwrap();
    for (m, species) in mass.iter().zip(&result.species_mass_kg) {
        assert!(species.iter().all(|&x| x >= 0.0 && x <= *m));
        near(species[0] + species[1], *m, 5e-15);
    }
    for k in 0..2 {
        near(result.report.mass_after_kg[k], result.report.mass_before_kg[k], 1e-14);
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
        assert!(state.diffuse_binary_species(
            &domain, &initial, &[bad, 0.1], 1.0, controls(1), &gate,
        ).is_err());
    }
    for ctl in [
        BinarySpeciesDiffusionControls { substeps: 0, ..controls(1) },
        BinarySpeciesDiffusionControls { maximum_storage_bytes: 1, ..controls(1) },
        BinarySpeciesDiffusionControls { maximum_pair_updates: 1, ..controls(1) },
    ] {
        assert!(state.diffuse_binary_species(
            &domain, &initial, &[0.1; 2], 1.0, ctl, &gate,
        ).is_err());
    }
    assert!(state.diffuse_binary_species(
        &domain, &[[1.0, 0.1], [0.0, 1.0]], &[0.1; 2], 1.0, controls(1), &gate,
    ).is_err());
    let mut changed = state.clone();
    changed.mass_kg[0] = 2.0;
    assert!(changed.diffuse_binary_species(
        &domain, &initial, &[0.1; 2], 1.0, controls(1), &gate,
    ).is_err());
    let mut changed_domain = domain.clone();
    changed_domain.set(1, 0, 0, Voxel::Solid(0));
    assert!(state.diffuse_binary_species(
        &changed_domain, &initial, &[0.1; 2], 1.0, controls(1), &gate,
    ).is_err());
    let original = state.clone();
    let error = state.diffuse_binary_species_checked(
        &domain, &initial, &[0.1; 2], 1.0, controls(4), &gate,
        |substep| if substep == 2 { Err(ChtError::Cancelled) } else { Ok(()) },
    ).unwrap_err();
    assert_eq!(error, ChtError::Cancelled);
    assert_eq!(state, original);
    assert_eq!(initial, [[1.0, 0.0], [0.0, 1.0]]);
}
