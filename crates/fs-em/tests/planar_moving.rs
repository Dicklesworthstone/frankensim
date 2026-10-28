//! G1/G3/G4: independently assembled moving-conductor MQS and physical work.
//! Synthetic material coefficients test the model, not a measured device.
#[path = "winding_ohmic/common.rs"]
mod material;

use fs_em::planar::{
    PlanarMagnetostatics,
    eddy::{EddyConductor, EddyControls, EddyError, PlanarEddy, moving::MovingEddyControls},
};
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_rep_mesh::{Metric2, TriComplex2};
use material::{bundle, card, limits, linear, near};

fn controls() -> MovingEddyControls {
    MovingEddyControls {
        relative_tolerance: 1e-10,
        absolute_residual_a_m: 1e-12,
        restart: 8,
        max_cycles: 50,
        max_krylov_scalars: 10000,
        max_cell_magnetic_peclet: 20.0,
        energy_tolerance_j: 1e-8,
        current_balance_tolerance_a: 1e-8,
    }
}
fn stationary_controls() -> EddyControls {
    EddyControls {
        linear: linear(),
        energy_tolerance_j: 1e-8,
        current_balance_tolerance_a: 1e-8,
    }
}

fn problem(mu: f64, sigma: f64, depth: f64) -> PlanarMagnetostatics {
    // Two free vertices (5,6); all other vertices lie on the rectangle boundary.
    let vertices: Vec<_> = (0..3)
        .flat_map(|j| (0..4).map(move |i| [f64::from(i), f64::from(j) / 2.0]))
        .collect();
    let mut faces = Vec::new();
    for j in 0..2 {
        for i in 0..3 {
            let a = 4 * j + i;
            faces.extend([[a, a + 1, a + 5], [a, a + 5, a + 4]]);
        }
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "moving-conductor-two-interior-modes",
        vertices,
        faces,
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    let basis = ElectricalPropertyBasis::Conductivity;
    let source = card(sigma.recip(), 0.0, mu, basis);
    PlanarMagnetostatics::new(
        mesh,
        vec![bundle(&source, 300.0, basis); 12],
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn conductors(problem: &PlanarMagnetostatics) -> Vec<EddyConductor> {
    problem
        .materials()
        .iter()
        .enumerate()
        .map(|(triangle, source)| EddyConductor {
            triangle,
            region: 9,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                source.clone(),
                ElectricalPropertyBasis::Conductivity,
            )
            .unwrap(),
        })
        .collect()
}
fn boundary(problem: &PlanarMagnetostatics, affine_x: f64, offset: f64) -> Vec<(usize, f64)> {
    problem
        .mesh()
        .vertices()
        .iter()
        .enumerate()
        .filter(|(v, _)| *v != 5 && *v != 6)
        .map(|(v, p)| (v, affine_x * p[0] + offset))
        .collect()
}

// Independent elementary triangle formulas; do not call the production field
// gradient, conductivity projection, Hodge or Krylov implementation.
fn geometry(points: [[f64; 2]; 3], depth: f64) -> (f64, [[f64; 2]; 3]) {
    let det = (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
        - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1]);
    let gradients = std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        [
            (points[j][1] - points[k][1]) / det,
            (points[k][0] - points[j][0]) / det,
        ]
    });
    (det.abs() * depth / 2.0, gradients)
}

fn dense_reference(
    problem: &PlanarMagnetostatics,
    a0: &[f64],
    velocity: &[[f64; 2]],
    mu: f64,
    sigma: f64,
    depth: f64,
    dt: f64,
) -> Vec<f64> {
    let n = a0.len();
    let mut k = vec![vec![0.0; n]; n];
    let mut c = k.clone();
    let mut v = k.clone();
    let mut moments = vec![0.0; n];
    let mut advected_moments = vec![0.0; n];
    let mut weight = 0.0;
    for face in problem.mesh().faces() {
        let f = face.map(|v| v as usize);
        let (volume, g) = geometry(f.map(|i| problem.mesh().vertices()[i]), depth);
        weight += sigma * volume;
        let velocity = f.map(|v| velocity[v]);
        for i in 0..3 {
            moments[f[i]] += sigma * volume / 3.0;
            advected_moments[f[i]] += sigma * volume / 3.0
                * velocity
                    .iter()
                    .map(|vel| vel[0] * g[i][0] + vel[1] * g[i][1])
                    .sum::<f64>();
            for j in 0..3 {
                k[f[i]][f[j]] += volume / mu * (g[i][0] * g[j][0] + g[i][1] * g[j][1]);
                c[f[i]][f[j]] += sigma * volume / 12.0 * if i == j { 2.0 } else { 1.0 };
                for (q, vel) in velocity.iter().enumerate() {
                    v[f[i]][f[j]] += sigma * volume / 12.0
                        * if i == q { 2.0 } else { 1.0 }
                        * (vel[0] * g[j][0] + vel[1] * g[j][1]);
                }
            }
        }
    }
    for i in 0..n {
        for j in 0..n {
            c[i][j] -= moments[i] * moments[j] / weight;
            v[i][j] -= moments[i] * advected_moments[j] / weight;
        }
    }
    let free = [5, 6];
    let matrix: [[f64; 2]; 2] =
        free.map(|i| free.map(|j| 0.5 * k[i][j] + c[i][j] / dt + 0.5 * v[i][j]));
    let rhs = free.map(|i| -(0..n).map(|j| (k[i][j] + v[i][j]) * a0[j]).sum::<f64>());
    let det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    let delta = [
        (rhs[0] * matrix[1][1] - matrix[0][1] * rhs[1]) / det,
        (matrix[0][0] * rhs[1] - rhs[0] * matrix[1][0]) / det,
    ];
    let mut a1 = a0.to_vec();
    for (i, d) in free.into_iter().zip(delta) {
        a1[i] += d;
    }
    a1
}

#[test]
fn moving_field_matches_independent_nonsymmetric_system_and_three_way_work() {
    for (mu, sigma, depth) in [(1.0, 18.0, 0.3), (2.0, 7.0, 0.7)] {
        let problem = problem(mu, sigma, depth);
        let boundary = boundary(&problem, 0.0, 0.0);
        let initial = problem
            .solve(&[12.0; 12], &boundary, linear(), &mut || false)
            .unwrap();
        let state = PlanarEddy::new(
            &initial,
            conductors(&problem),
            &[],
            &boundary,
            0.0,
            stationary_controls(),
            &mut || false,
        )
        .unwrap();
        // Includes a shear component so motional EMF is P1 within each cell.
        let velocity: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|p| [0.2 + 0.1 * p[1], 0.0])
            .collect();
        let dt = 0.08;
        let next = state
            .step_with_eulerian_velocity_candidate(
                dt,
                &velocity,
                &[0.0; 12],
                controls(),
                &mut || false,
            )
            .unwrap();
        let reference = dense_reference(
            &problem,
            initial.potential_wb_m(),
            &velocity,
            mu,
            sigma,
            depth,
            dt,
        );
        for (&actual, expected) in next.state().potential_wb_m().iter().zip(reference) {
            near(actual, expected);
        }
        let stationary = state
            .step_to_candidate(
                dt,
                state.conductors().to_vec(),
                &[0.0; 12],
                stationary_controls(),
                &mut || false,
            )
            .unwrap();
        assert!(
            (next.state().potential_wb_m()[5] - stationary.state().potential_wb_m()[5]).abs()
                > 1e-5
        );
        assert!(
            (next.state().potential_wb_m()[6] - stationary.state().potential_wb_m()[6]).abs()
                > 1e-5
        );

        let stage: Vec<_> = initial
            .potential_wb_m()
            .iter()
            .zip(next.state().potential_wb_m())
            .map(|(a, b)| 0.5 * (a + b))
            .collect();
        let rate: Vec<_> = initial
            .potential_wb_m()
            .iter()
            .zip(next.state().potential_wb_m())
            .map(|(a, b)| (b - a) / dt)
            .collect();
        let mut mean = 0.0;
        let mut volume_sum = 0.0;
        let mut fields = Vec::new();
        for face in problem.mesh().faces() {
            let f = face.map(|v| v as usize);
            let (volume, g) = geometry(f.map(|i| problem.mesh().vertices()[i]), depth);
            let grad: [f64; 2] =
                std::array::from_fn(|axis| (0..3).map(|i| stage[f[i]] * g[i][axis]).sum());
            let u = f.map(|i| velocity[i][0] * grad[0] + velocity[i][1] * grad[1]);
            let r: [f64; 3] = std::array::from_fn(|i| rate[f[i]] + u[i]);
            mean += volume * r.iter().sum::<f64>() / 3.0;
            volume_sum += volume;
            fields.push((volume, grad, u, r));
        }
        mean /= volume_sum;
        let mut heat = 0.0;
        let mut work = 0.0;
        for (actual, (volume, grad, u, r)) in next.step().cells.iter().zip(fields) {
            let mut cell_heat = 0.0;
            let mut cell_work = 0.0;
            // Independent degree-two triangle quadrature, not the mass formula.
            for q in 0..3 {
                let bary: [f64; 3] =
                    std::array::from_fn(|i| if i == q { 2.0 / 3.0 } else { 1.0 / 6.0 });
                let rate = (0..3).map(|i| bary[i] * r[i]).sum::<f64>() - mean;
                let motion = (0..3).map(|i| bary[i] * u[i]).sum::<f64>();
                cell_heat += dt * volume / 3.0 * sigma * rate * rate;
                cell_work -= dt * volume / 3.0 * sigma * rate * motion;
            }
            near(actual.joule.heat_j, cell_heat);
            near(actual.mechanical_work_j, cell_work);
            for i in 0..3 {
                near(actual.joule.electric_field_v_m[i], mean - r[i]);
                near(
                    actual.laboratory_electric_field_v_m[i] + actual.motional_electric_field_v_m[i],
                    actual.joule.electric_field_v_m[i],
                );
                for (axis, &grad) in grad.iter().enumerate() {
                    near(
                        actual.lorentz_force_density_n_m3[i][axis],
                        sigma * (mean - r[i]) * grad,
                    );
                }
            }
            heat += cell_heat;
            work += cell_work;
        }
        near(next.step().heat_j, heat);
        near(next.step().mechanical_work_j, work);
        assert!(work < -1e-5); // The imposed flow does work against magnetic braking.
        near(
            next.step().stored_after_j - next.step().stored_before_j + heat + work,
            0.0,
        );
        near(next.step().energy_residual_j, 0.0);
        near(next.step().region_currents_a[0].1, 0.0);
        assert!(next.step().report.as_ref().unwrap().converged_euclidean());
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(state.potential_wb_m(), initial.potential_wb_m());
    }
}

#[test]
fn uniform_motional_emf_is_cancelled_by_floating_charge_without_false_heating() {
    let problem = problem(1.0, 18.0, 0.3);
    let boundary = boundary(&problem, 1.0, 0.0);
    let initial = problem
        .solve(&[0.0; 12], &boundary, linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem),
        &[],
        &boundary,
        0.0,
        stationary_controls(),
        &mut || false,
    )
    .unwrap();
    let mut controls = controls();
    controls.max_cell_magnetic_peclet = 100.0;
    let next = state
        .step_with_eulerian_velocity_candidate(
            0.1,
            &[[0.5, 0.0]; 12],
            &[0.0; 12],
            controls,
            &mut || false,
        )
        .unwrap();
    for (&a, &b) in next
        .state()
        .potential_wb_m()
        .iter()
        .zip(state.potential_wb_m())
    {
        near(a, b);
    }
    near(next.step().heat_j, 0.0);
    near(next.step().mechanical_work_j, 0.0);
    for cell in &next.step().cells {
        for i in 0..3 {
            near(cell.motional_electric_field_v_m[i], -0.5);
            near(cell.laboratory_electric_field_v_m[i], 0.5);
            near(cell.joule.current_density_a_m2[i], 0.0);
        }
    }
}

#[test]
fn zero_velocity_retains_stationary_decay_and_continuation() {
    let problem = problem(1.0, 18.0, 0.3);
    let boundary = boundary(&problem, 0.0, 0.0);
    let initial = problem
        .solve(&[12.0; 12], &boundary, linear(), &mut || false)
        .unwrap();
    let mut moving = PlanarEddy::new(
        &initial,
        conductors(&problem),
        &[],
        &boundary,
        0.0,
        stationary_controls(),
        &mut || false,
    )
    .unwrap();
    let mut stationary = moving.clone();
    for step in 1..=6 {
        let next = moving
            .step_with_eulerian_velocity_candidate(
                f64::from(step) * 0.1,
                &[[0.0; 2]; 12],
                &[0.0; 12],
                controls(),
                &mut || false,
            )
            .unwrap();
        let reference = stationary
            .step_to_candidate(
                f64::from(step) * 0.1,
                stationary.conductors().to_vec(),
                &[0.0; 12],
                stationary_controls(),
                &mut || false,
            )
            .unwrap();
        for (&a, &b) in next
            .state()
            .potential_wb_m()
            .iter()
            .zip(reference.state().potential_wb_m())
        {
            near(a, b);
        }
        near(next.step().heat_j, reference.step().heat_j);
        near(next.step().mechanical_work_j, 0.0);
        moving = next.into_parts().0;
        stationary = reference.into_parts().0;
    }
}

#[test]
fn velocity_resolution_work_and_cancellation_refusals_keep_accepted_history() {
    let problem = problem(1.0, 18.0, 0.3);
    let boundary = boundary(&problem, 0.0, 0.0);
    let initial = problem
        .solve(&[12.0; 12], &boundary, linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem),
        &[],
        &boundary,
        0.0,
        stationary_controls(),
        &mut || false,
    )
    .unwrap();
    for velocity in [
        vec![[0.0; 2]; 11],
        vec![[f64::NAN, 0.0]; 12],
        vec![[100.0, 0.0]; 12],
    ] {
        assert!(
            state
                .step_with_eulerian_velocity_candidate(
                    0.1,
                    &velocity,
                    &[0.0; 12],
                    controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    let mut bounded = controls();
    bounded.max_krylov_scalars = 1;
    assert!(
        state
            .step_with_eulerian_velocity_candidate(
                0.1,
                &[[0.2, 0.0]; 12],
                &[0.0; 12],
                bounded,
                &mut || false
            )
            .is_err()
    );
    bounded = controls();
    bounded.restart = 1;
    bounded.max_cycles = 1;
    assert!(
        state
            .step_with_eulerian_velocity_candidate(
                0.1,
                &[[0.2, 0.0]; 12],
                &[0.0; 12],
                bounded,
                &mut || false
            )
            .is_err()
    );
    assert!(
        state
            .step_with_eulerian_velocity_candidate(
                0.1,
                &[[0.2, 0.0]; 12],
                &[1.0; 12],
                controls(),
                &mut || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        state.step_with_eulerian_velocity_candidate(
            0.1,
            &[[0.2, 0.0]; 12],
            &[0.0; 12],
            controls(),
            &mut || {
                polls += 1;
                polls == 110
            }
        ),
        Err(EddyError::Cancelled)
    ));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb_m(), initial.potential_wb_m());
}

#[test]
fn changing_balanced_impressed_source_closes_actual_field_heat_and_motion_work() {
    let problem = problem(1.0, 18.0, 0.3);
    let boundary = boundary(&problem, 0.0, 0.0);
    let initial = problem
        .solve(&[0.0; 12], &boundary, linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem),
        &[],
        &boundary,
        0.0,
        stationary_controls(),
        &mut || false,
    )
    .unwrap();
    let current: Vec<_> = problem
        .mesh()
        .faces()
        .iter()
        .map(|face| {
            let x = face
                .iter()
                .map(|&v| problem.mesh().vertices()[v as usize][0])
                .sum::<f64>()
                / 3.0;
            if x < 1.5 { 6.0 } else { -6.0 }
        })
        .collect();
    let next = state
        .step_with_eulerian_velocity_candidate(
            0.15,
            &[[0.25, 0.0]; 12],
            &current,
            controls(),
            &mut || false,
        )
        .unwrap();
    let mut source_work = 0.0;
    for (&current, face) in current.iter().zip(problem.mesh().faces()) {
        let f = face.map(|v| v as usize);
        let (volume, _) = geometry(f.map(|i| problem.mesh().vertices()[i]), 0.3);
        source_work += 0.5 * current * volume / 3.0
            * f.iter()
                .map(|&i| next.state().potential_wb_m()[i] - state.potential_wb_m()[i])
                .sum::<f64>();
    }
    assert!(source_work > 0.0);
    near(next.step().supplied_work_j, source_work);
    near(
        next.step().stored_after_j + next.step().heat_j + next.step().mechanical_work_j,
        source_work,
    );
    assert_eq!(next.state().impressed_current_density_a_m2(), current);
}
