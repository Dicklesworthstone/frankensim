//! G0/G1 prescribed-field induction: independent Faraday and tetrahedron oracles.

use super::*;

#[test]
fn uniform_insulating_field_ramp_has_signed_boundary_work_and_no_fictitious_heat() {
    let (positions, tets) = grid(2);
    let count = tets.len();
    let (material, _) = source(2.0, None, None);
    let problem = SpatialVectorMagnetostatics::new(
        positions,
        tets,
        vec![material; count],
        assembly(),
        &mut || false,
    )
    .unwrap();
    let current = vec![[0.0; 3]; count];
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let before = [0.12, -0.16, 0.2];
    let after = [0.3, -0.08, 0.45];
    let initial = problem
        .solve(
            &current,
            &trace(&problem, before, false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(
        &initial,
        vec![None; count],
        preparation(),
        1e-9,
        &mut || false,
    )
    .unwrap();
    let forward = state
        .step_with_boundary_candidate(
            &current,
            &velocity,
            &trace(&problem, after, false),
            0.05,
            controls(),
            &mut || false,
        )
        .unwrap();
    let energy = |b: [f64; 3]| b.iter().map(|v| v * v).sum::<f64>() / 4.0;
    for b in &forward.step().endpoint_flux_density_t {
        for j in 0..3 {
            near(b[j], after[j], 2e-9);
        }
    }
    near(forward.step().stored_after_j, energy(after), 2e-10);
    near(
        forward.step().magnetic_boundary_work_j,
        energy(after) - energy(before),
        2e-10,
    );
    assert_eq!(forward.step().heat_j, 0.0);
    assert_eq!(forward.step().impressed_work_j, 0.0);
    assert_eq!(forward.step().terminal_work_j, 0.0);
    let reverse = forward
        .state()
        .step_with_boundary_candidate(
            &current,
            &velocity,
            &trace(&problem, before, false),
            0.05,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(reverse.step().magnetic_boundary_work_j < 0.0);
    near(
        reverse.step().magnetic_boundary_work_j,
        energy(before) - energy(after),
        2e-10,
    );
    near(reverse.step().energy_residual_j, 0.0, 2e-10);
    assert_eq!(state.time_s(), 0.0);
}

#[test]
fn boundary_faraday_field_drives_the_exact_one_tet_charge_projection_and_ohmic_heat() {
    let (material, conductor) = source(2.0, Some(5.0), None);
    let problem = SpatialVectorMagnetostatics::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        vec![[0, 1, 2, 3]],
        vec![material],
        assembly(),
        &mut || false,
    )
    .unwrap();
    let current = [[0.0; 3]];
    let velocity = [[0.0; 3]; 4];
    let initial = problem
        .solve(
            &current,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(&initial, vec![conductor], preparation(), 1e-9, &mut || {
        false
    })
    .unwrap();
    let result = state
        .step_with_boundary_candidate(
            &current,
            &velocity,
            &trace(&problem, [0.0, 0.0, 0.2], false),
            0.1,
            controls(),
            &mut || false,
        )
        .unwrap();
    let step = result.step();
    let cell = &step.cells[0];
    // For Bdot=(0,0,2), E_F=(y,-x,0). The P1 floating-potential
    // equations remove its volume mean, giving E=(y-1/4,1/4-x,0).
    // The simplex moments E[x²]=E[y²]=1/10 and E[x]=E[y]=1/4
    // give heat = sigma*dt*V*(var(x)+var(y)) = 1/160 joule.
    for (local, &vertex) in problem.mesh().tets[0].iter().enumerate() {
        let p = problem.positions_m()[vertex as usize];
        let expected = [p[1] - 0.25, 0.25 - p[0], 0.0];
        for (j, &e) in expected.iter().enumerate() {
            near(cell.laboratory_electric_field_v_m[local][j], e, 2e-10);
            near(cell.current_density_a_m2[local][j], 5.0 * e, 2e-10);
        }
    }
    near(step.heat_j, 1.0 / 160.0, 2e-11);
    near(step.stored_after_j, 0.2_f64.powi(2) / 24.0, 2e-11);
    near(
        step.magnetic_boundary_work_j,
        step.stored_after_j + step.heat_j,
        2e-11,
    );
    near(step.maximum_charge_residual_a, 0.0, 2e-10);
    assert_eq!(step.impressed_work_j, 0.0);
    assert_eq!(step.mechanical_work_j, 0.0);
    // Independently integrate the returned laboratory E along each oriented
    // face. This checks the SIGN of Faraday's law, not just positive heating.
    for face in &problem.mesh().faces {
        let positions = face.map(|v| problem.positions_m()[v as usize]);
        let electric = face.map(|v| {
            let local = problem.mesh().tets[0].iter().position(|&w| w == v).unwrap();
            cell.laboratory_electric_field_v_m[local]
        });
        let mut circulation = 0.0;
        for k in 0..3 {
            let next = (k + 1) % 3;
            for j in 0..3 {
                circulation += 0.5
                    * (electric[k][j] + electric[next][j])
                    * (positions[next][j] - positions[k][j]);
            }
        }
        let a = std::array::from_fn(|j| positions[1][j] - positions[0][j]);
        let b = std::array::from_fn(|j| positions[2][j] - positions[0][j]);
        near(circulation, -cross(a, b)[2], 2e-10);
    }
}

#[test]
fn applied_ramp_then_hold_retains_diffusion_and_matches_the_original_fixed_trace_path() {
    let (problem, conductors) = cube(2, None);
    let current = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let initial = problem
        .solve(
            &current,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, conductors, preparation(), 1e-9, &mut || false).unwrap();
    let endpoint = trace(&problem, [0.0, 0.0, 0.1], false);
    let ramp = state
        .step_with_boundary_candidate(
            &current,
            &velocity,
            &endpoint,
            0.05,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(ramp.step().heat_j > 1e-6);
    assert!(ramp.step().magnetic_boundary_work_j > ramp.step().heat_j);
    let hold = ramp
        .state()
        .step_candidate(&current, &velocity, 0.05, controls(), &mut || false)
        .unwrap();
    let explicit_hold = ramp
        .state()
        .step_with_boundary_candidate(
            &current,
            &velocity,
            &endpoint,
            0.05,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        hold.state().potential_wb(),
        explicit_hold.state().potential_wb()
    );
    assert_eq!(hold.step().heat_j, explicit_hold.step().heat_j);
    assert_eq!(hold.step().magnetic_boundary_work_j, 0.0);
    assert_eq!(explicit_hold.step().magnetic_boundary_work_j, 0.0);
    assert!(hold.step().heat_j > 1e-9);
    assert!(hold.step().stored_after_j < hold.step().stored_before_j);
    near(
        hold.step().stored_before_j - hold.step().stored_after_j,
        hold.step().heat_j,
        2e-9,
    );
    assert_eq!(ramp.state().tree_edges(), state.tree_edges());
    assert_eq!(state.time_s(), 0.0);
    assert!(state.potential_wb().iter().all(|&v| v == 0.0));
}

#[test]
fn a_time_dependent_boundary_gradient_is_cancelled_by_charge_and_does_not_heat() {
    let (problem, conductors) = cube(2, None);
    let b = [0.2, -0.1, 0.3];
    let current = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let initial = problem
        .solve(
            &current,
            &trace(&problem, b, false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, conductors, preparation(), 1e-9, &mut || false).unwrap();
    let gauge = state
        .step_with_boundary_candidate(
            &current,
            &velocity,
            &trace(&problem, b, true),
            0.1,
            controls(),
            &mut || false,
        )
        .unwrap();
    for field in &gauge.step().endpoint_flux_density_t {
        for j in 0..3 {
            near(field[j], b[j], 2e-8);
        }
    }
    assert!(gauge.step().heat_j < 1e-18);
    near(gauge.step().magnetic_boundary_work_j, 0.0, 2e-9);
    for c in &gauge.step().cells {
        assert!(
            c.current_density_a_m2
                .iter()
                .flatten()
                .all(|j| j.abs() < 2e-8)
        );
    }
}

#[test]
fn boundary_support_conflicts_and_cancellation_leave_the_accepted_history_untouched() {
    let (problem, conductors) = cube(1, None);
    let current = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let initial = problem
        .solve(
            &current,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, conductors, preparation(), 1e-9, &mut || false).unwrap();
    let valid = trace(&problem, [0.0, 0.0, 0.1], false);
    let mut inconsistent = valid.clone();
    inconsistent.potential[0].edge_integrals_wb[0] += 0.1;
    let mut nonfinite = valid.clone();
    nonfinite.potential[0].edge_integrals_wb[0] = f64::NAN;
    for boundary in [VectorMagneticBoundary::default(), inconsistent, nonfinite] {
        assert!(matches!(
            state.step_with_boundary_candidate(
                &current,
                &velocity,
                &boundary,
                0.05,
                controls(),
                &mut || false,
            ),
            Err(SpatialMagneticError::Invalid(_))
        ));
    }
    let mut polls = 0;
    assert!(matches!(
        state.step_with_boundary_candidate(
            &current,
            &velocity,
            &valid,
            0.05,
            controls(),
            &mut || {
                polls += 1;
                polls == 50
            },
        ),
        Err(SpatialMagneticError::Cancelled)
    ));
    assert_eq!(state.time_s(), 0.0);
    assert!(state.potential_wb().iter().all(|&v| v == 0.0));
}
