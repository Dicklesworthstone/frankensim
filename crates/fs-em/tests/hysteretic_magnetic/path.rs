//! G1/G3/G4: continuous-path references, independent gauges and atomic refusal.

use super::*;
use fs_em::planar::nonlinear::hysteretic::path::{MagneticPathControls, MagneticPathError};

fn path_controls() -> MagneticPathControls {
    MagneticPathControls {
        field: controls(),
        relative_tolerance: 0.0,
        absolute_polarization_tolerance_t: 1e-4,
        absolute_heat_tolerance_j: 1e-3,
        absolute_work_tolerance_j: 0.02,
        energy_tolerance_j: 2e-7,
        maximum_field_solves: 30000,
        maximum_newton_iterations: 60000,
        maximum_refinement_depth: 24,
    }
}

fn linear_problem() -> PlanarHystereticMagnetostatics {
    PlanarHystereticMagnetostatics::new(
        square(0.7),
        vec![source(0.01, false); 4],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn linear_path_resolves_work_without_inventing_physical_heat() {
    let owner = linear_problem();
    let virgin = owner.demagnetized(0.0, &mut || false).unwrap();
    let solve = |tolerance| {
        let mut c = path_controls();
        c.absolute_work_tolerance_j = tolerance;
        owner
            .affine_path_candidate(
                &virgin,
                &[0.0; 4],
                &[36.0; 4],
                &boundary(0.0),
                c,
                &mut || false,
            )
            .unwrap()
    };
    let coarse = solve(0.02);
    let fine = solve(0.002);
    let exact_storage = 0.7 * 0.005 * 6.0_f64.powi(2);
    for candidate in [&coarse, &fine] {
        let report = candidate.report();
        near(candidate.state().stored_energy_j(), exact_storage);
        near(report.account.stored_energy_change_j, exact_storage);
        assert_eq!(report.account.irreversible_heat_j, 0.0);
        assert!(report.account.numerical_dissipation_j > 0.0);
        near(
            report.account.supplied_work_j,
            exact_storage + report.account.numerical_dissipation_j,
        );
        near(report.account.boundary_work_j, 0.0);
        assert!(report.summary.maximum_accepted_defect_ratio <= 1.0);
        assert!(report.summary.field_solves > report.summary.accepted_substeps);
        assert!(report.summary.newton_iterations > 0);
    }
    assert!(
        fine.report().account.numerical_dissipation_j
            < 0.2 * coarse.report().account.numerical_dissipation_j
    );
    assert!(fine.report().summary.accepted_substeps > coarse.report().summary.accepted_substeps);
    assert!(fine.report().summary.estimated_supplied_work_error_j <= 0.002 * (1.0 + 1e-12));
    let saved = fine.report().summary;
    let (native, summary) = fine.into_magnetic_candidate();
    assert_eq!(summary, saved);
    assert_eq!(native.step().irreversible_heat_j, 0.0);
    near(native.state().stored_energy_j(), exact_storage);
}

#[test]
fn saturated_major_loop_retains_exact_heat_and_native_minor_loop_memory() {
    let owner = problem(0.7);
    let virgin = owner.demagnetized(0.0, &mut || false).unwrap();
    let (mut state, _) = owner
        .candidate(&virgin, &[72.0; 4], &boundary(0.0), controls(), &mut || {
            false
        })
        .unwrap()
        .into_parts();
    let start = state.clone();
    let (mut heat, mut numerical, mut supplied) = (0.0, 0.0, 0.0);
    let mut c = path_controls();
    c.absolute_work_tolerance_j = 0.15;
    for (from, to) in [(72.0, -72.0), (-72.0, 72.0)] {
        let path = owner
            .affine_path_candidate(&state, &[from; 4], &[to; 4], &boundary(0.0), c, &mut || {
                false
            })
            .unwrap();
        heat += path.report().account.irreversible_heat_j;
        numerical += path.report().account.numerical_dissipation_j;
        supplied += path.report().account.supplied_work_j;
        near(path.report().account.energy_residual_j, 0.0);
        (state, _) = path.into_parts();
    }
    // Scalar play gives total variation 4*P_sat and Q=H_c*variation*volume.
    near(heat, 0.7 * 4.0 * 2.0);
    near(state.stored_energy_j(), start.stored_energy_j());
    near(supplied, heat + numerical);
    assert!(numerical > 0.0 && numerical < 0.3);
    for (a, b) in state.cell_states().iter().zip(start.cell_states()) {
        for axis in 0..3 {
            near(a.polarizations_t()[0][axis], b.polarizations_t()[0][axis]);
        }
    }
    // Nontrivial nested reversals must be separate authored path segments.
    let mut current = 72.0;
    let mut remanence = None;
    for (i, h) in [0.0, 6.0, 1.0, 5.5, 1.0, 6.0, 0.0].into_iter().enumerate() {
        let path = owner
            .affine_path_candidate(
                &state,
                &[current; 4],
                &[6.0 * h; 4],
                &boundary(0.0),
                c,
                &mut || false,
            )
            .unwrap();
        (state, _) = path.into_parts();
        if i == 0 {
            remanence = Some(state.cell_states().to_vec());
        }
        if i == 6 {
            for (a, b) in state.cell_states().iter().zip(remanence.as_ref().unwrap()) {
                for axis in 0..3 {
                    near(a.polarizations_t()[0][axis], b.polarizations_t()[0][axis]);
                }
            }
        }
        current = 6.0 * h;
    }
}

fn trace(flux: [f64; 2]) -> Vec<(usize, f64)> {
    vertices()
        .iter()
        .take(4)
        .enumerate()
        .map(|(i, &[x, y])| (i, flux[0] * y - flux[1] * x))
        .collect()
}

#[test]
fn vector_loading_converges_toward_independently_uniform_native_history() {
    let depth = 0.4;
    // In B-control the vector-play radius is mu*chi*hc/(mu+chi).
    // This fit retains directional memory over the four-radius excursion;
    // mu=0.01 would exhaust it over 81 radii and compare endpoint roundoff.
    let law = source(0.2, true);
    let mesh = TriComplex2::from_indexed_triangles(
        "magnetic-path-rotation",
        vertices()[..4].to_vec(),
        vec![[0, 1, 2], [0, 2, 3]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    let owner =
        PlanarHystereticMagnetostatics::new(mesh, vec![law.clone(); 2], limits(), &mut || false)
            .unwrap();
    let virgin = owner.demagnetized(0.0, &mut || false).unwrap();
    let a = [0.8, 0.0];
    let b = [-0.6, 0.8];
    let (initial, _) = owner
        .candidate(&virgin, &[0.0; 2], &trace(a), controls(), &mut || false)
        .unwrap()
        .into_parts();
    // This reference advances the original constitutive B-driven law uniformly,
    // without the adaptive algorithm, field solver or its substep selection.
    let mut reference = initial.cell_states()[0].clone();
    let (mut reference_heat, mut reference_work) = (0.0, 0.0);
    for i in 1..=4096 {
        let t = i as f64 / 4096.0;
        let flux = [
            a[0] * (1.0 - t) + b[0] * t,
            a[1] * (1.0 - t) + b[1] * t,
            0.0,
        ];
        let step = law
            .law()
            .candidate_flux_density(&reference, flux, controls().material, &mut || false)
            .unwrap();
        reference_heat += depth * step.irreversible_heat_j_m3;
        reference_work += depth * step.endpoint_work_j_m3;
        reference = step.state;
    }
    let run = |scale| {
        let mut c = path_controls();
        c.absolute_polarization_tolerance_t = scale * 0.002;
        c.absolute_heat_tolerance_j = scale * 0.01;
        c.absolute_work_tolerance_j = scale * 0.05;
        owner
            .affine_path_candidate(&initial, &[0.0; 2], &[0.0; 2], &trace(b), c, &mut || false)
            .unwrap()
    };
    let coarse = run(1.0);
    let fine = run(0.1);
    let error = |state: &fs_em::planar::nonlinear::hysteretic::PlanarHystereticState<'_>| {
        let p = state.cell_states()[0].polarizations_t()[0];
        let q = reference.polarizations_t()[0];
        (p[0] - q[0]).hypot(p[1] - q[1])
    };
    assert!(error(fine.state()) < error(coarse.state()));
    assert!(error(fine.state()) < 0.0003);
    assert!(
        (fine.report().account.irreversible_heat_j - reference_heat).abs()
            < (coarse.report().account.irreversible_heat_j - reference_heat).abs()
    );
    assert!(
        (fine.report().account.supplied_work_j - reference_work).abs()
            < (coarse.report().account.supplied_work_j - reference_work).abs()
    );
    assert!(fine.report().summary.accepted_substeps > coarse.report().summary.accepted_substeps);
    assert_eq!(fine.report().summary.newton_iterations, 0);
    near(fine.report().account.impressed_current_work_j, 0.0);
    near(
        fine.report().account.boundary_work_j,
        fine.report().account.supplied_work_j,
    );
    near(fine.report().account.energy_residual_j, 0.0);
}

#[test]
fn continuation_reaches_admissible_boundary_loading_after_coarse_domain_refusal() {
    let owner = problem(0.6);
    let initial = owner.demagnetized(0.0, &mut || false).unwrap();
    let target = trace([1.25, 0.0]);
    // The existing endpoint guess retains interior A=0, creating a spurious
    // initial B=2.5 T; the true uniform endpoint B=1.25 T is inside the fit.
    assert!(matches!(
        owner.candidate(&initial, &[0.0; 4], &target, controls(), &mut || false,),
        Err(fs_em::planar::nonlinear::NonlinearMagneticError::History {
            error:
                fs_material::field_storage::hysteresis::vector::VectorHysteresisError::FluxRange { .. },
            ..
        })
    ));
    let mut c = path_controls();
    c.absolute_work_tolerance_j = 0.1;
    let path = owner
        .affine_path_candidate(&initial, &[0.0; 4], &[0.0; 4], &target, c, &mut || false)
        .unwrap();
    assert!(path.report().summary.failed_field_solves > 0);
    assert!(path.report().summary.rejected_intervals > 0);
    for (i, potential) in path.state().potential_wb_m().iter().enumerate() {
        near(*potential, 1.25 * vertices()[i][1]);
    }
    for cell in path.state().cell_states() {
        near(cell.flux_density_t()[0], 1.25);
        near(cell.field_a_m()[0], 25.0);
        near(cell.polarizations_t()[0][0], 1.0);
    }
    near(path.report().account.irreversible_heat_j, 0.6 * 2.0);
    near(path.report().account.energy_residual_j, 0.0);
    let mut refused = c;
    refused.maximum_refinement_depth = 0;
    assert!(matches!(
        owner.affine_path_candidate(
            &initial,
            &[0.0; 4],
            &[0.0; 4],
            &target,
            refused,
            &mut || false,
        ),
        Err(MagneticPathError::FieldRefinementLimit {
            error: fs_em::planar::nonlinear::NonlinearMagneticError::History { .. },
            ..
        })
    ));
    refused = c;
    refused.maximum_field_solves = 1;
    assert!(matches!(
        owner.affine_path_candidate(
            &initial,
            &[0.0; 4],
            &[0.0; 4],
            &target,
            refused,
            &mut || false,
        ),
        Err(MagneticPathError::FieldSolveBudget)
    ));
    refused = c;
    refused.maximum_newton_iterations = 0;
    assert!(matches!(
        owner.affine_path_candidate(
            &initial,
            &[0.0; 4],
            &[0.0; 4],
            &target,
            refused,
            &mut || false,
        ),
        Err(MagneticPathError::NewtonBudget)
    ));
    assert_eq!(initial.stored_energy_j(), 0.0);
}

#[test]
fn source_admission_work_caps_and_late_cancellation_preserve_accepted_history() {
    let owner = linear_problem();
    let initial = owner.demagnetized(0.0, &mut || false).unwrap();
    let c = path_controls();
    let run = |c, cancelled: &mut dyn FnMut() -> bool| {
        owner.affine_path_candidate(
            &initial,
            &[0.0; 4],
            &[12.0; 4],
            &boundary(0.0),
            c,
            &mut || cancelled(),
        )
    };
    assert!(matches!(
        owner.affine_path_candidate(
            &initial,
            &[12.0; 4],
            &[24.0; 4],
            &boundary(0.0),
            c,
            &mut || false,
        ),
        Err(MagneticPathError::InitialEquilibrium { .. })
    ));
    let mut exhausted = c;
    exhausted.maximum_field_solves = 2;
    assert!(matches!(
        run(exhausted, &mut || false),
        Err(MagneticPathError::FieldSolveBudget)
    ));
    exhausted = c;
    exhausted.maximum_newton_iterations = 0;
    assert!(matches!(
        run(exhausted, &mut || false),
        Err(MagneticPathError::NewtonBudget)
    ));
    exhausted = c;
    exhausted.maximum_refinement_depth = 0;
    exhausted.absolute_work_tolerance_j = 1e-8;
    assert!(matches!(
        run(exhausted, &mut || false),
        Err(MagneticPathError::RefinementLimit { .. })
    ));
    let mut calls = 0;
    let accepted = run(c, &mut || {
        calls += 1;
        false
    })
    .unwrap();
    let mut interrupted = 0;
    assert!(matches!(
        run(c, &mut || {
            interrupted += 1;
            interrupted == calls - 1
        }),
        Err(MagneticPathError::Cancelled)
    ));
    assert_eq!(initial.stored_energy_j(), 0.0);
    assert!(initial.potential_wb_m().iter().all(|&v| v == 0.0));
    assert!(
        initial
            .cell_states()
            .iter()
            .all(|s| s.field_a_m() == [0.0; 3])
    );
    let replay = run(c, &mut || false).unwrap();
    assert_eq!(accepted.state().cell_states(), replay.state().cell_states());
    assert_eq!(
        accepted.state().potential_wb_m(),
        replay.state().potential_wb_m()
    );
    assert_eq!(accepted.report().summary, replay.report().summary);
    assert_eq!(
        accepted.report().account.supplied_work_j.to_bits(),
        replay.report().account.supplied_work_j.to_bits()
    );
}

#[test]
fn independently_shifted_component_gauges_preserve_the_loading_path() {
    let mut coordinates = vertices();
    coordinates.extend(vertices().iter().map(|&[x, y]| [x + 2.0, y]));
    let mesh = TriComplex2::from_indexed_triangles(
        "magnetic-path-independent-gauges",
        coordinates,
        vec![
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
            [5, 6, 9],
            [6, 7, 9],
            [7, 8, 9],
            [8, 5, 9],
        ],
        Metric2::planar(0.7).unwrap(),
    )
    .unwrap();
    let owner = PlanarHystereticMagnetostatics::new(
        mesh,
        vec![source(0.01, false); 8],
        limits(),
        &mut || false,
    )
    .unwrap();
    let initial = owner.demagnetized(0.0, &mut || false).unwrap();
    let currents = [12.0, 12.0, 12.0, 12.0, -24.0, -24.0, -24.0, -24.0];
    let run = |shift| {
        let target = [0, 1, 2, 3, 5, 6, 7, 8].map(|v| {
            (
                v,
                if shift {
                    if v < 5 { 7.0 } else { -11.0 }
                } else {
                    0.0
                },
            )
        });
        owner
            .affine_path_candidate(
                &initial,
                &[0.0; 8],
                &currents,
                &target,
                path_controls(),
                &mut || false,
            )
            .unwrap()
    };
    let reference = run(false);
    let shifted = run(true);
    for (i, (&a, &b)) in shifted
        .state()
        .potential_wb_m()
        .iter()
        .zip(reference.state().potential_wb_m())
        .enumerate()
    {
        near(a - b, if i < 5 { 7.0 } else { -11.0 });
    }
    for (a, b) in shifted
        .state()
        .magnetic_flux_density_t()
        .iter()
        .zip(reference.state().magnetic_flux_density_t())
    {
        near(a[0], b[0]);
        near(a[1], b[1]);
    }
    near(
        shifted.report().account.supplied_work_j,
        reference.report().account.supplied_work_j,
    );
    near(
        shifted.report().account.numerical_dissipation_j,
        reference.report().account.numerical_dissipation_j,
    );
    near(shifted.report().account.boundary_work_j, 0.0);
    assert_eq!(shifted.report().account.irreversible_heat_j, 0.0);
}
