//! G1/G3: actual source H(B), isolated currents, spatial Joule heat and the
//! separate backward-Euler loss. Manufactured laws are not metal qualification.

#[path = "nonlinear_eddy/common.rs"]
mod common;
use common::{boundary, conductor, conductors, controls, mesh, near, problem, source};
use fs_em::planar::{
    eddy::{
        EddyConductor,
        nonlinear::{NonlinearEddyError, PlanarNonlinearEddy},
    },
    nonlinear::NonlinearMagneticError,
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_rep_mesh::{Metric2, TriComplex2};

#[test]
fn linear_limit_matches_independent_backward_euler_and_separates_heat_from_numerical_loss() {
    for sigma in [18.0, 36.0] {
        for (offset, sign) in [(0.0, 1.0), (100.0, -1.0)] {
            let card = source(&[(0.0, 0.0), (4.0, 4.0)], sigma);
            let problem = problem(&card, mesh(0.3));
            let initial = problem
                .solve(
                    &[sign * 12.0; 4],
                    &boundary(offset),
                    &[offset; 5],
                    controls().nonlinear,
                    &mut || false,
                )
                .unwrap();
            let state = PlanarNonlinearEddy::new(
                &initial,
                conductors(&card, 4, 300.0, ElectricalPropertyBasis::Resistivity),
                &[],
                &boundary(offset),
                0.0,
                controls(),
                &mut || false,
            )
            .unwrap();
            let dt = 0.1;
            let candidate = state
                .step_to_candidate(
                    dt,
                    state.conductors().to_vec(),
                    &[0.0; 4],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            let c = sigma * 0.3 / 18.0;
            let k = 4.0 * 0.3;
            let a = c / (c + dt * k);
            let delta = a - 1.0;
            near(candidate.state().potential_wb_m()[4] - offset, sign * a);
            let step = candidate.step();
            near(step.stored_before_j, 0.5 * k);
            near(step.stored_after_j, 0.5 * k * a * a);
            near(step.heat_j, c * delta * delta / dt);
            near(step.numerical_dissipation_j, 0.5 * k * delta * delta);
            assert!(step.numerical_dissipation_j > 0.0);
            assert!(step.stored_before_j - step.stored_after_j > step.heat_j);
            near(step.energy_residual_j, 0.0);
            near(step.supplied_work_j, 0.0);
            near(step.region_currents_a[0].1, 0.0);
            for cell in &step.cells {
                near(cell.heat_j, step.heat_j / 4.0);
                for (local, &vertex) in problem.mesh().faces()[cell.source.triangle]
                    .iter()
                    .enumerate()
                {
                    let e = sign * delta / dt * if vertex == 4 { -2.0 / 3.0 } else { 1.0 / 3.0 };
                    near(cell.electric_field_v_m[local], e);
                    near(cell.current_density_a_m2[local], sigma * e);
                }
            }
            let field = candidate.state().field(&mut || false).unwrap();
            near(field.energy_j(), step.stored_after_j);
            for (before, after) in initial
                .magnetic_flux_density_t()
                .iter()
                .zip(field.magnetic_flux_density_t())
            {
                for axis in 0..2 {
                    near(after[axis], a * before[axis]);
                }
            }
            assert!(step.report.residual_norm_a_m <= step.report.residual_target_a_m);
        }
    }
}

#[test]
fn changing_only_h_b_data_changes_actual_transient_field_and_joule_heat() {
    let mut heat = Vec::new();
    for saturating in [false, true] {
        let knots = if saturating {
            vec![(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]
        } else {
            vec![(0.0, 0.0), (4.0, 4.0)]
        };
        let card = source(&knots, 18.0);
        let problem = problem(&card, mesh(0.3));
        let initial = problem
            .solve(
                &[12.0; 4],
                &boundary(0.0),
                &[0.0; 5],
                controls().nonlinear,
                &mut || false,
            )
            .unwrap();
        let state = PlanarNonlinearEddy::new(
            &initial,
            conductors(&card, 4, 300.0, ElectricalPropertyBasis::Conductivity),
            &[],
            &boundary(0.0),
            0.0,
            controls(),
            &mut || false,
        )
        .unwrap();
        let dt = 0.02;
        let next = state
            .step_to_candidate(
                dt,
                state.conductors().to_vec(),
                &[0.0; 4],
                controls(),
                &mut || false,
            )
            .unwrap();
        let a0 = if saturating { 0.625 } else { 1.0 };
        // Independent scalar P1 equation: 2*depth*H(2a) + C44*(a-a0)/dt=0.
        let a1 = if saturating {
            (15.0 * a0 + 1.8) / 19.8
        } else {
            15.0 * a0 / 16.2
        };
        near(next.state().potential_wb_m()[4], a1);
        let delta = a1 - a0;
        near(next.step().heat_j, 15.0 * delta * delta);
        near(
            next.step().numerical_dissipation_j,
            if saturating {
                2.4 * delta * delta
            } else {
                0.6 * delta * delta
            },
        );
        let b = 2.0 * a1;
        let psi = if saturating {
            0.5 + (b - 1.0) + 2.0 * (b - 1.0).powi(2)
        } else {
            0.5 * b * b
        };
        near(next.step().stored_after_j, 0.3 * psi);
        near(next.step().energy_residual_j, 0.0);
        heat.push(next.step().heat_j);
        assert_eq!(
            next.step().cells[0]
                .source
                .conductor
                .material()
                .card_identity(),
            card.content_hash()
        );
    }
    assert!((heat[0] - heat[1]).abs() > 1e-3);
}

#[test]
fn time_refinement_reduces_the_reported_numerical_loss_and_recovers_continuum_decay() {
    let card = source(&[(0.0, 0.0), (4.0, 4.0)], 18.0);
    let problem = problem(&card, mesh(0.3));
    let initial = problem
        .solve(
            &[12.0; 4],
            &boundary(0.0),
            &[0.0; 5],
            controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let original = PlanarNonlinearEddy::new(
        &initial,
        conductors(&card, 4, 300.0, ElectricalPropertyBasis::Resistivity),
        &[],
        &boundary(0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let mut previous_error = f64::INFINITY;
    let mut previous_loss = f64::INFINITY;
    for steps in [1, 2, 4, 8] {
        let mut state = original.clone();
        let mut loss = 0.0;
        let mut heat = 0.0;
        for i in 0..steps {
            let candidate = state
                .step_to_candidate(
                    (i + 1) as f64 * 0.1 / steps as f64,
                    state.conductors().to_vec(),
                    &[0.0; 4],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            loss += candidate.step().numerical_dissipation_j;
            heat += candidate.step().heat_j;
            state = candidate.into_parts().0;
        }
        let a = (1.0 + 0.4 / steps as f64).powi(-steps);
        near(state.potential_wb_m()[4], a);
        near(state.energy_j() + heat + loss, original.energy_j());
        let error = (a - (-0.4_f64).exp()).abs();
        assert!(error < previous_error && loss < previous_loss);
        previous_error = error;
        previous_loss = loss;
    }
}

#[test]
fn conductivity_refresh_changes_dynamics_without_changing_the_magnetic_source() {
    let card = source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)], 18.0);
    let problem = problem(&card, mesh(0.3));
    let initial = problem
        .solve(
            &[12.0; 4],
            &boundary(0.0),
            &[0.0; 5],
            controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let state = PlanarNonlinearEddy::new(
        &initial,
        conductors(&card, 4, 300.0, ElectricalPropertyBasis::Resistivity),
        &[],
        &boundary(0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let cold = state
        .step_to_candidate(
            0.02,
            state.conductors().to_vec(),
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    let hot = state
        .step_to_candidate(
            0.02,
            conductors(&card, 4, 500.0, ElectricalPropertyBasis::Resistivity),
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(hot.state().potential_wb_m()[4] < cold.state().potential_wb_m()[4]);
    near(hot.step().energy_residual_j, 0.0);
    assert!(std::ptr::eq(hot.state().problem(), state.problem()));
    let switched = conductors(&card, 4, 300.0, ElectricalPropertyBasis::Conductivity);
    assert!(
        state
            .step_to_candidate(0.02, switched, &[0.0; 4], controls(), &mut || false)
            .is_err()
    );
    let changed_card = source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)], 36.0);
    assert!(
        state
            .step_to_candidate(
                0.02,
                conductors(
                    &changed_card,
                    4,
                    300.0,
                    ElectricalPropertyBasis::Resistivity
                ),
                &[0.0; 4],
                controls(),
                &mut || false
            )
            .is_err()
    );
}

#[test]
fn balanced_excitation_reaches_conducting_and_insulating_cells_at_the_actual_endpoint() {
    let n = 4;
    let mut vertices = Vec::new();
    let mut boundary = Vec::new();
    for y in 0..=n {
        for x in 0..=n {
            if x == 0 || y == 0 || x == n || y == n {
                boundary.push((vertices.len(), 0.0));
            }
            vertices.push([x as f64 / n as f64, y as f64 / n as f64]);
        }
    }
    let mut faces = Vec::new();
    let mut conducting = Vec::new();
    let mut insulating = Vec::new();
    let mut current = Vec::new();
    for y in 0..n {
        for x in 0..n {
            let a = (y * (n + 1) + x) as u32;
            let b = a + 1;
            let d = a + (n + 1) as u32;
            let c = d + 1;
            for face in [[a, b, c], [a, c, d]] {
                if x < n / 2 {
                    conducting.push(faces.len());
                } else {
                    insulating.push(faces.len());
                }
                current.push(if x < n / 2 { 12.0 } else { -12.0 });
                faces.push(face);
            }
        }
    }
    let card = source(&[(0.0, 0.0), (0.02, 0.02), (4.0, 19.92)], 18.0);
    let mesh = TriComplex2::from_indexed_triangles(
        "conducting-and-insulating-endpoints",
        vertices,
        faces,
        Metric2::planar(0.3).unwrap(),
    )
    .unwrap();
    let problem = problem(&card, mesh);
    let zero_current = vec![0.0; current.len()];
    let initial = problem
        .solve(
            &zero_current,
            &boundary,
            &vec![0.0; problem.mesh().vertices().len()],
            controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let material = conductor(&card, 300.0, ElectricalPropertyBasis::Resistivity);
    let conductors = conducting
        .into_iter()
        .map(|triangle| EddyConductor {
            triangle,
            region: 7,
            conductor: material.clone(),
        })
        .collect();
    let state = PlanarNonlinearEddy::new(
        &initial,
        conductors,
        &insulating,
        &boundary,
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let candidate = state
        .step_to_candidate(
            0.05,
            state.conductors().to_vec(),
            &current,
            controls(),
            &mut || false,
        )
        .unwrap();
    let step = candidate.step();
    assert!(step.stored_after_j > 0.0 && step.heat_j > 0.0 && step.numerical_dissipation_j > 0.0);
    near(
        step.supplied_work_j,
        step.stored_after_j + step.heat_j + step.numerical_dissipation_j,
    );
    let independent_load_work = problem
        .mesh()
        .faces()
        .iter()
        .enumerate()
        .map(|(cell, face)| {
            current[cell] * problem.mesh().face_measure(cell).unwrap() / 3.0
                * face
                    .iter()
                    .map(|&v| candidate.state().potential_wb_m()[v as usize])
                    .sum::<f64>()
        })
        .sum::<f64>();
    near(step.supplied_work_j, independent_load_work);
    assert!(
        step.cells
            .iter()
            .all(|cell| !insulating.contains(&cell.source.triangle))
    );
    // Recover the endpoint weak magnetic equation independently in pure-air
    // free hats, where no conduction capacity can hide an algebraic error.
    let field = candidate.state().field(&mut || false).unwrap();
    for y in 1..n {
        let vertex = y * (n + 1) + 3;
        let mut residual = 0.0;
        for (cell, face) in problem.mesh().faces().iter().enumerate() {
            let Some(local) = face.iter().position(|&v| v as usize == vertex) else {
                continue;
            };
            let p = face.map(|v| problem.mesh().vertices()[v as usize]);
            let determinant = (p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
                - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
            let j = (local + 1) % 3;
            let k = (local + 2) % 3;
            let grad = [
                (p[j][1] - p[k][1]) / determinant,
                (p[k][0] - p[j][0]) / determinant,
            ];
            let h = field.magnetic_field_a_m()[cell];
            residual += problem.mesh().face_measure(cell).unwrap()
                * (h[0] * grad[1] - h[1] * grad[0] - current[cell] / 3.0);
        }
        near(residual, 0.0);
    }
    let next = candidate.into_parts().0;
    let reversed: Vec<_> = current.iter().map(|j| -j).collect();
    let reversed = next
        .step_to_candidate(
            0.1,
            next.conductors().to_vec(),
            &reversed,
            controls(),
            &mut || false,
        )
        .unwrap();
    near(reversed.step().energy_residual_j, 0.0);
    near(reversed.step().region_currents_a[0].1, 0.0);
}

#[test]
fn source_balance_nonlinear_budget_and_late_cancellation_preserve_accepted_history() {
    let card = source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)], 18.0);
    let problem = problem(&card, mesh(0.3));
    let initial = problem
        .solve(
            &[12.0; 4],
            &boundary(0.0),
            &[0.0; 5],
            controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let state = PlanarNonlinearEddy::new(
        &initial,
        conductors(&card, 4, 300.0, ElectricalPropertyBasis::Resistivity),
        &[],
        &boundary(0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let original = state.potential_wb_m().to_vec();
    assert!(
        state
            .step_to_candidate(
                0.02,
                state.conductors().to_vec(),
                &[1.0; 4],
                controls(),
                &mut || false
            )
            .is_err()
    );
    let mut limited = controls();
    limited.nonlinear.max_newton_iterations = 0;
    assert!(matches!(
        state.step_to_candidate(
            0.02,
            state.conductors().to_vec(),
            &[0.0; 4],
            limited,
            &mut || false
        ),
        Err(NonlinearEddyError::Magnetic(
            NonlinearMagneticError::Nonconvergence(_)
        ))
    ));
    let mut calls = 0;
    let successful = state
        .step_to_candidate(
            0.02,
            state.conductors().to_vec(),
            &[0.0; 4],
            controls(),
            &mut || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut visited = 0;
    let cancelled = state.step_to_candidate(
        0.02,
        state.conductors().to_vec(),
        &[0.0; 4],
        controls(),
        &mut || {
            visited += 1;
            visited == calls
        },
    );
    assert!(matches!(cancelled, Err(NonlinearEddyError::Cancelled)));
    assert_eq!(state.potential_wb_m(), original);
    assert_eq!(state.time_s(), 0.0);
    let retry = state
        .step_to_candidate(
            0.02,
            state.conductors().to_vec(),
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        retry.state().potential_wb_m(),
        successful.state().potential_wb_m()
    );
    assert_eq!(
        retry.step().heat_j.to_bits(),
        successful.step().heat_j.to_bits()
    );
}
