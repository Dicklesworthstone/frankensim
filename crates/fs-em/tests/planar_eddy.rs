//! G1/G3/G4: independent P1 induction, source work and material-only changes.
//! These are manufactured coefficients, not a measured metal qualification.
#[path = "winding_ohmic/common.rs"]
mod material;
#[path = "planar_eddy/thermal.rs"]
mod thermal;

use fs_em::planar::{
    PlanarMagnetostatics,
    eddy::{EddyConductor, EddyControls, PlanarEddy},
};
use fs_matdb::MaterialCard;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_rep_mesh::{Metric2, TriComplex2};
use material::{bundle, card, limits, linear, near, point};

fn controls() -> EddyControls {
    EddyControls {
        linear: linear(),
        energy_tolerance_j: 1e-8,
        current_balance_tolerance_a: 1e-9,
    }
}

fn square(
    source: &MaterialCard,
    basis: ElectricalPropertyBasis,
    depth: f64,
) -> PlanarMagnetostatics {
    let mesh = TriComplex2::from_indexed_triangles(
        "eddy-current-four-triangle-mode",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    PlanarMagnetostatics::new(
        mesh,
        vec![bundle(source, 300.0, basis); 4],
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn conductors(
    problem: &PlanarMagnetostatics,
    basis: ElectricalPropertyBasis,
) -> Vec<EddyConductor> {
    problem
        .materials()
        .iter()
        .enumerate()
        .map(|(triangle, material)| EddyConductor {
            triangle,
            region: 7,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                material.clone(),
                basis,
            )
            .unwrap(),
        })
        .collect()
}

fn boundary(offset: f64) -> Vec<(usize, f64)> {
    (0..4).map(|vertex| (vertex, offset)).collect()
}

#[test]
fn isolated_conductor_decay_uses_consistent_mass_and_the_floating_electric_potential() {
    let mut outcomes = Vec::new();
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        for sigma in [18.0, 36.0] {
            let source = card(1.0 / sigma, 0.0, 1.0, basis);
            let depth = 0.3;
            let problem = square(&source, basis, depth);
            // The original magnetostatic solve constructs a nonzero initial field.
            let initial = problem
                .solve(&[12.0; 4], &boundary(0.0), linear(), &mut || false)
                .unwrap();
            near(initial.potential_wb_m()[4], 1.0);
            let state = PlanarEddy::new(
                &initial,
                conductors(&problem, basis),
                &[],
                &boundary(0.0),
                0.0,
                controls(),
                &mut || false,
            )
            .unwrap();
            let dt = 0.1;
            let next = state
                .step_to_candidate(
                    dt,
                    state.conductors().to_vec(),
                    &[0.0; 4],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            // Analytic P1 moments on a unit square: integral(N4)=1/3,
            // integral(N4^2)=1/6. Eliminating the region's floating potential
            // gives C44=sigma*depth*(1/6-1/9), not sigma*depth/6.
            let c = sigma * depth / 18.0;
            let k = 4.0 * depth;
            let amplitude = (c - dt * k / 2.0) / (c + dt * k / 2.0);
            let derivative = (amplitude - 1.0) / dt;
            near(next.state().potential_wb_m()[4], amplitude);
            let step = next.step();
            near(step.stored_before_j, 0.5 * k);
            near(step.stored_after_j, 0.5 * k * amplitude * amplitude);
            near(step.heat_j, dt * c * derivative * derivative);
            near(step.stored_before_j - step.stored_after_j, step.heat_j);
            near(step.supplied_work_j, 0.0);
            near(step.energy_residual_j, 0.0);
            assert_eq!(step.region_currents_a.len(), 1);
            assert_eq!(step.region_currents_a[0].0, 7);
            near(step.region_currents_a[0].1, 0.0);
            for cell in &step.cells {
                let face = problem.mesh().faces()[cell.source.triangle];
                for (local, &vertex) in face.iter().enumerate() {
                    let electric = derivative * if vertex == 4 { -2.0 / 3.0 } else { 1.0 / 3.0 };
                    near(cell.electric_field_v_m[local], electric);
                    near(cell.current_density_a_m2[local], sigma * electric);
                }
                near(cell.volume_m3, depth / 4.0);
                near(cell.heat_j, step.heat_j / 4.0);
                near(cell.power_w, cell.heat_j / dt);
                assert_eq!(cell.source.conductor.source_basis(), basis);
                assert_eq!(
                    cell.source.conductor.material().card_identity(),
                    source.content_hash()
                );
            }
            let field = next.state().field(&mut || false).unwrap();
            near(field.energy_j(), step.stored_after_j);
            for (a, b) in field
                .magnetic_flux_density_t()
                .iter()
                .zip(initial.magnetic_flux_density_t())
            {
                for axis in 0..2 {
                    near(a[axis], amplitude * b[axis]);
                }
            }
            near(state.potential_wb_m()[4], 1.0);
            assert_eq!(state.time_s(), 0.0);
            outcomes.push(amplitude);
        }
    }
    assert!(outcomes[1] > outcomes[0]); // Higher conductivity delays diffusion.
    near(outcomes[0], outcomes[2]);
    near(outcomes[1], outcomes[3]);
}

#[test]
fn retained_induction_has_second_order_midpoint_decay_and_pinned_temperature_refresh() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let source = card(1.0 / 18.0, 0.001, 1.0, basis);
    let problem = square(&source, basis, 0.5);
    let initial = problem
        .solve(&[12.0; 4], &boundary(0.0), linear(), &mut || false)
        .unwrap();
    let mut errors = Vec::new();
    for n in [10, 20, 40] {
        let mut state = PlanarEddy::new(
            &initial,
            conductors(&problem, basis),
            &[],
            &boundary(0.0),
            0.0,
            controls(),
            &mut || false,
        )
        .unwrap();
        let mut heat = 0.0;
        for step in 1..=n {
            let candidate = state
                .step_to_candidate(
                    0.5 * step as f64 / n as f64,
                    state.conductors().to_vec(),
                    &[0.0; 4],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            heat += candidate.step().heat_j;
            state = candidate.into_parts().0;
        }
        near(state.energy_j() + heat, initial.energy_j());
        errors.push((state.potential_wb_m()[4] - (-2.0_f64).exp()).abs());
        assert_eq!(state.time_s(), 0.5);
    }
    assert!(errors[0] > 3.8 * errors[1]);
    assert!(errors[1] > 3.8 * errors[2]);

    let state = PlanarEddy::new(
        &initial,
        conductors(&problem, basis),
        &[],
        &boundary(0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let refreshed: Vec<_> = state
        .conductors()
        .iter()
        .map(|old| EddyConductor {
            triangle: old.triangle,
            region: old.region,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                old.conductor
                    .material()
                    .at_query_point(&source, &point(400.0))
                    .unwrap(),
                basis,
            )
            .unwrap(),
        })
        .collect();
    let hot = state
        .step_to_candidate(0.1, refreshed, &[0.0; 4], controls(), &mut || false)
        .unwrap();
    let c = (18.0 / 1.1) * 0.5 / 18.0;
    let k = 4.0 * 0.5;
    near(
        hot.state().potential_wb_m()[4],
        (c - 0.1 * k / 2.0) / (c + 0.1 * k / 2.0),
    );
    near(hot.step().stored_before_j, initial.energy_j());
    for cell in hot.state().conductors() {
        near(cell.conductor.resistivity_ohm_m(), 1.1 / 18.0);
        assert_eq!(
            cell.conductor.material().card_identity(),
            source.content_hash()
        );
    }
}

#[test]
fn balanced_impressed_current_builds_the_same_retained_field_and_then_decays_as_heat() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [
        card(1.0 / 18.0, 0.0, 1.0, basis),
        card(1.0 / 18.0, 0.0, 1.0, basis),
    ];
    let (problem, profiles, boundary_nodes) = material::problem(&cards, basis, 0.3);
    let bc: Vec<_> = boundary_nodes.iter().map(|&node| (node, 0.0)).collect();
    let zero = vec![0.0; problem.mesh().faces().len()];
    let initial = problem.solve(&zero, &bc, linear(), &mut || false).unwrap();
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem, basis),
        &[],
        &bc,
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let candidate = state
        .step_to_candidate(
            0.001,
            state.conductors().to_vec(),
            &profiles[0],
            controls(),
            &mut || false,
        )
        .unwrap();
    let driven = candidate.step();
    assert!(driven.supplied_work_j > 0.0);
    assert!(driven.stored_after_j > 0.0);
    assert!(driven.heat_j > 0.0);
    near(
        driven.supplied_work_j,
        driven.stored_after_j + driven.heat_j,
    );
    let independent_work: f64 = problem
        .mesh()
        .faces()
        .iter()
        .enumerate()
        .map(|(i, face)| {
            // This fixture's uniform 8x8 mesh has triangle volume depth/128.
            // The accepted source starts at zero and ramps to this endpoint.
            0.5 * profiles[0][i] * 0.3 / 128.0
                * face
                    .iter()
                    .map(|&v| candidate.state().potential_wb_m()[v as usize] / 3.0)
                    .sum::<f64>()
        })
        .sum();
    near(driven.supplied_work_j, independent_work);
    let (running, _) = candidate.into_parts();
    let ramp_down = running
        .step_to_candidate(
            0.002,
            running.conductors().to_vec(),
            &zero,
            controls(),
            &mut || false,
        )
        .unwrap();
    near(
        ramp_down.step().stored_after_j - ramp_down.step().stored_before_j
            + ramp_down.step().heat_j,
        ramp_down.step().supplied_work_j,
    );
    assert_eq!(ramp_down.state().impressed_current_density_a_m2(), zero);
    let (unforced, _) = ramp_down.into_parts();
    let coast = unforced
        .step_to_candidate(
            0.003,
            unforced.conductors().to_vec(),
            &zero,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(coast.step().stored_after_j < coast.step().stored_before_j);
    near(
        coast.step().stored_before_j - coast.step().stored_after_j,
        coast.step().heat_j,
    );
    near(coast.step().supplied_work_j, 0.0);
    near(state.energy_j(), 0.0);
}

#[test]
fn insulating_algebraic_field_follows_both_source_endpoints_without_alternating_history() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [
        card(1.0 / 18.0, 0.0, 1.0, basis),
        card(1.0 / 18.0, 0.0, 1.0, basis),
    ];
    let (problem, profiles, nodes) = material::problem(&cards, basis, 0.3);
    let bc: Vec<_> = nodes.into_iter().map(|v| (v, 0.0)).collect();
    let insulators: Vec<_> = (0..problem.mesh().faces().len()).collect();
    let zero = vec![0.0; insulators.len()];
    let initial = problem.solve(&zero, &bc, linear(), &mut || false).unwrap();
    let excited = problem
        .solve(&profiles[0], &bc, linear(), &mut || false)
        .unwrap();
    let mut limits = controls();
    limits.current_balance_tolerance_a = 1e-8;
    // Electrical roles are explicit: this section has no modeled conductor.
    // Every free magnetic coordinate is therefore an algebraic air coordinate.
    assert!(
        PlanarEddy::new(&excited, vec![], &insulators, &bc, 0.0, limits, &mut || {
            false
        })
        .is_err()
    );
    let consistent = PlanarEddy::new_with_impressed_current(
        &excited,
        vec![],
        &insulators,
        &bc,
        &profiles[0],
        0.0,
        limits,
        &mut || false,
    )
    .unwrap();
    assert_eq!(consistent.potential_wb_m(), excited.potential_wb_m());
    let mut state = PlanarEddy::new(&initial, vec![], &insulators, &bc, 0.0, limits, &mut || {
        false
    })
    .unwrap();
    for (time, input, expected) in [
        (0.001, profiles[0].as_slice(), excited.potential_wb_m()),
        (0.002, zero.as_slice(), initial.potential_wb_m()),
    ] {
        let candidate = state
            .step_to_candidate(time, vec![], input, limits, &mut || false)
            .unwrap();
        for (&actual, &reference) in candidate.state().potential_wb_m().iter().zip(expected) {
            near(actual, reference);
        }
        assert_eq!(candidate.state().impressed_current_density_a_m2(), input);
        assert_eq!(candidate.step().heat_j, 0.0);
        near(
            candidate.step().stored_after_j - candidate.step().stored_before_j,
            candidate.step().supplied_work_j,
        );
        state = candidate.into_parts().0;
    }
    near(state.energy_j(), 0.0);
}

#[test]
fn field_and_material_history_survive_source_budget_cancel_and_topology_refusals() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let source = card(1.0 / 18.0, 0.001, 1.0, basis);
    let problem = square(&source, basis, 0.5);
    let initial = problem
        .solve(&[12.0; 4], &boundary(0.0), linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem, basis),
        &[],
        &boundary(0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let old = state.potential_wb_m().to_vec();
    let energy = state.energy_j();
    assert!(
        state
            .step_to_candidate(
                0.1,
                state.conductors().to_vec(),
                &[1.0; 4],
                controls(),
                &mut || false
            )
            .is_err()
    ); // No hidden impressed-current return.
    let mut starved = controls();
    starved.linear.max_iterations = 0;
    assert!(
        state
            .step_to_candidate(
                0.1,
                state.conductors().to_vec(),
                &[0.0; 4],
                starved,
                &mut || false
            )
            .is_err()
    );
    let mut checks = 0;
    assert!(
        state
            .step_to_candidate(
                0.1,
                state.conductors().to_vec(),
                &[0.0; 4],
                controls(),
                &mut || {
                    checks += 1;
                    checks > 8
                }
            )
            .is_err()
    );
    assert!(checks > 8);
    let other = card(1.0 / 36.0, 0.001, 1.0, basis);
    let mut replaced = state.conductors().to_vec();
    replaced[0].conductor =
        ResolvedConductorState::try_from_material_with_basis(bundle(&other, 300.0, basis), basis)
            .unwrap();
    assert!(
        state
            .step_to_candidate(0.1, replaced, &[0.0; 4], controls(), &mut || false)
            .is_err()
    );
    let mut disconnected = state.conductors().to_vec();
    disconnected[1].region = 8;
    disconnected[3].region = 8;
    assert!(
        PlanarEddy::new(
            &initial,
            disconnected,
            &[],
            &boundary(0.0),
            0.0,
            controls(),
            &mut || false
        )
        .is_err()
    ); // Opposite triangles share only a point.
    assert_eq!(state.potential_wb_m(), old);
    near(state.energy_j(), energy);
    assert_eq!(state.time_s(), 0.0);
    let a = state
        .step_to_candidate(
            0.1,
            state.conductors().to_vec(),
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    let b = state
        .step_to_candidate(
            0.1,
            state.conductors().to_vec(),
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(a.state().potential_wb_m(), b.state().potential_wb_m());
    near(a.step().heat_j, b.step().heat_j);
}

#[test]
fn conductivity_contrast_uses_fewer_iterations_without_changing_field_or_heat() {
    use fs_solver::{LinearOp, checked_cg};
    use fs_sparse::precond::IdentityPrecond;

    struct Reference(Vec<Vec<f64>>);
    impl LinearOp for Reference {
        fn n(&self) -> usize {
            self.0.len()
        }
        fn apply(&self, x: &[f64], y: &mut [f64]) {
            for (row, value) in self.0.iter().zip(y) {
                *value = row.iter().zip(x).map(|(a, b)| a * b).sum();
            }
        }
        // Its one instance below is a symmetric matrix (stiffness, mass and the moment outer product are each assembled symmetrically), so its transpose is itself.
        fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
            self.apply(x, y);
        }
    }

    // A connected conductor with five decades of source-resolved conductivity
    // contrast. These are manufactured coefficients, not measured metal data.
    // The reference assembles elementary P1 integrals independently of the
    // production FEEC Hodge, projected-mass operator and Jacobi preconditioner.
    let n = 12;
    let depth = 0.3;
    let mu = 1e-6;
    let dt = 1e-5;
    let basis = ElectricalPropertyBasis::Conductivity;
    let mut vertices = Vec::new();
    let mut boundary = Vec::new();
    let mut free = Vec::new();
    for y in 0..=n {
        for x in 0..=n {
            let v = vertices.len();
            if x == 0 || y == 0 || x == n || y == n {
                boundary.push((v, 0.0));
            } else {
                free.push(v);
            }
            vertices.push([x as f64 / n as f64, y as f64 / n as f64]);
        }
    }
    let sources: Vec<_> = (0..n)
        .map(|x| {
            let sigma = 100.0 * 3.0_f64.powi(x as i32);
            bundle(&card(sigma.recip(), 0.0, mu, basis), 300.0, basis)
        })
        .collect();
    let mut faces = Vec::new();
    let mut materials = Vec::new();
    let mut sigma = Vec::new();
    let mut current = Vec::new();
    for y in 0..n {
        for (x, source) in sources.iter().enumerate() {
            let a = (y * (n + 1) + x) as u32;
            let d = a + (n + 1) as u32;
            faces.extend([[a, a + 1, d + 1], [a, d + 1, d]]);
            materials.extend([source.clone(), source.clone()]);
            sigma.extend([100.0 * 3.0_f64.powi(x as i32); 2]);
            current.extend([if y < n / 2 { 1e4 } else { -1e4 }; 2]);
        }
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "induction-conductivity-contrast",
        vertices,
        faces,
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    let problem = PlanarMagnetostatics::new(mesh, materials, limits(), &mut || false).unwrap();
    let count = problem.mesh().vertices().len();
    let mut stiffness = vec![vec![0.0; count]; count];
    let mut mass = stiffness.clone();
    let mut moments = vec![0.0; count];
    let mut load = vec![0.0; count];
    let mut weight = 0.0;
    let mut geometry = Vec::new();
    for (cell, face) in problem.mesh().faces().iter().enumerate() {
        let f = face.map(|v| v as usize);
        let p = f.map(|v| problem.mesh().vertices()[v]);
        let det =
            (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
        let volume = det.abs() * depth / 2.0;
        let g: [[f64; 2]; 3] = std::array::from_fn(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            [(p[j][1] - p[k][1]) / det, (p[k][0] - p[j][0]) / det]
        });
        geometry.push((volume, g));
        weight += sigma[cell] * volume;
        for i in 0..3 {
            moments[f[i]] += sigma[cell] * volume / 3.0;
            load[f[i]] += 0.5 * current[cell] * volume / 3.0;
            for j in 0..3 {
                stiffness[f[i]][f[j]] += volume / mu * (g[i][0] * g[j][0] + g[i][1] * g[j][1]);
                mass[f[i]][f[j]] += sigma[cell] * volume / 12.0 * if i == j { 2.0 } else { 1.0 };
            }
        }
    }
    let reference = Reference(
        free.iter()
            .map(|&i| {
                free.iter()
                    .map(|&j| {
                        0.5 * stiffness[i][j] + (mass[i][j] - moments[i] * moments[j] / weight) / dt
                    })
                    .collect()
            })
            .collect(),
    );
    let rhs: Vec<_> = free.iter().map(|&v| load[v]).collect();
    let unpreconditioned = checked_cg(&reference, &IdentityPrecond, &rhs, linear(), |_| {
        Ok::<(), ()>(())
    })
    .unwrap();
    assert!(unpreconditioned.report.converged_euclidean());
    let mut expected = vec![0.0; count];
    for (&v, &value) in free.iter().zip(&unpreconditioned.x) {
        expected[v] = value;
    }

    let zero = vec![0.0; current.len()];
    let initial = problem
        .solve(&zero, &boundary, linear(), &mut || false)
        .unwrap();
    let control = EddyControls {
        current_balance_tolerance_a: 1e-7,
        ..controls()
    };
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem, basis),
        &[],
        &boundary,
        0.0,
        control,
        &mut || false,
    )
    .unwrap();
    let next = state
        .step_to_candidate(
            dt,
            state.conductors().to_vec(),
            &current,
            control,
            &mut || false,
        )
        .unwrap();
    let report = next.step().report.as_ref().unwrap();
    assert!(report.converged_euclidean());
    eprintln!(
        "conductivity contrast: Jacobi {} vs identity {} CG iterations",
        report.iters, unpreconditioned.report.iters,
    );
    assert!(
        3 * report.iters < unpreconditioned.report.iters,
        "Jacobi {} vs identity {} iterations",
        report.iters,
        unpreconditioned.report.iters,
    );
    let field_scale = expected.iter().map(|v| v.abs()).fold(0.0, f64::max);
    for (&actual, &expected) in next.state().potential_wb_m().iter().zip(&expected) {
        assert!((actual - expected).abs() < 2e-8 * field_scale);
    }

    let mean = moments
        .iter()
        .zip(&expected)
        .map(|(m, a)| m * a)
        .sum::<f64>()
        / weight;
    let field = next.state().field(&mut || false).unwrap();
    let mut heat = 0.0;
    let mut stored = 0.0;
    for (cell, face) in problem.mesh().faces().iter().enumerate() {
        let f = face.map(|v| v as usize);
        let (volume, g) = geometry[cell];
        let grad: [f64; 2] =
            std::array::from_fn(|axis| (0..3).map(|i| expected[f[i]] * g[i][axis]).sum());
        let b = [grad[1], -grad[0]];
        for (&actual, expected) in field.magnetic_flux_density_t()[cell].iter().zip(b) {
            assert!((actual - expected).abs() < 2e-7 * field_scale);
        }
        let electric = f.map(|v| (mean - expected[v]) / dt);
        let q = dt * sigma[cell] * volume / 12.0
            * (electric.iter().sum::<f64>().powi(2) + electric.iter().map(|e| e * e).sum::<f64>());
        heat += q;
        stored += volume / (2.0 * mu) * (b[0] * b[0] + b[1] * b[1]);
    }
    let work = load.iter().zip(&expected).map(|(f, a)| f * a).sum::<f64>();
    assert!((next.step().heat_j - heat).abs() < 2e-8 * heat);
    assert!((next.step().stored_after_j - stored).abs() < 2e-8 * stored);
    assert!((next.step().supplied_work_j - work).abs() < 2e-8 * work);
    assert!((stored + heat - work).abs() < 1e-9 * work);
    near(next.step().energy_residual_j, 0.0);
    near(next.step().region_currents_a[0].1, 0.0);
}
