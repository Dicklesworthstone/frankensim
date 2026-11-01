//! G0/G1/G3: original enriched director inertia and thermal midpoint motion.
//! Synthetic coefficients are numerical oracles, not material validation.

use fs_material::sheet::SheetMatrix;
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::transient::EnrichedMitcShell;
use fs_plate::shell::mitc::thermal::{FacetThermalStrain, ShellThermalStrain};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_sparse::Csr;
use std::sync::Arc;

const I0: f64 = 2.0;
const I1: f64 = 0.1;
const I2: f64 = 0.04;
const AREA: f64 = 0.5;

fn make_model(interior_only: bool, mass: ShellMass) -> Arc<MitcShellModel> {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    let a = SheetMatrix::orthotropic(2.0e4, 1.0e4, 0.2, 3.0e3).unwrap();
    let d = SheetMatrix::orthotropic(20.0, 10.0, 0.2, 3.0).unwrap();
    let section =
        ShearShellSection::new_coupled(a, [0.0; 9], d, [120.0, 0.0, 80.0], 1.0, I0, I1, I2, 0.0)
            .unwrap();
    let constraints = if interior_only {
        (0..3)
            .map(|node| ShellConstraint {
                node,
                // Leave one numerical drilling coordinate free. Physical nodal
                // director motion is zero, so the two interior modes are isolated.
                mask: [true, true, true, true, true, node != 0],
                frame: ConstraintFrame::Global,
            })
            .collect()
    } else {
        vec![]
    };
    Arc::new(
        assemble_mitc3p(
            &mesh,
            &directors,
            &[section],
            &constraints,
            mass,
            MitcBudget {
                max_nodes: 3,
                max_triangles: 1,
            },
            1e-3,
            || false,
        )
        .unwrap(),
    )
}

fn fields(scale: f64) -> Vec<FacetThermalStrain> {
    vec![FacetThermalStrain::AtNodes([
        ShellThermalStrain::try_new([0.0; 3], [0.0; 3]).unwrap(),
        ShellThermalStrain::try_new([0.0; 3], [scale, 0.0, 0.0]).unwrap(),
        ShellThermalStrain::try_new([0.0; 3], [0.0, 0.4 * scale, 0.0]).unwrap(),
    ])]
}
fn entry(matrix: &Csr, row: usize, column: usize) -> f64 {
    let (cols, vals) = matrix.row(row);
    cols.iter()
        .zip(vals)
        .find_map(|(&c, &v)| (c == column).then_some(v))
        .unwrap_or(0.0)
}
fn near(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_exact_enriched_inertia_retains_cross_moment_and_original_nodal_stiffness() {
    let model = make_model(false, ShellMass::Consistent);
    let enriched = EnrichedMitcShell::try_new(Arc::clone(&model), || false).unwrap();
    let n = model.pencil().free;
    assert_eq!(enriched.dimension(), n + 2);
    assert!(Arc::ptr_eq(&model, &enriched.shared_model()));
    // On the reference triangle b=27 lambda0 lambda1 lambda2, so
    // mean(b)=9/20 and mean(b²)=81/280, independently of implementation quadrature.
    for i in 0..2 {
        for j in 0..2 {
            near(
                entry(enriched.mass(), n + i, n + j),
                if i == j {
                    AREA * I2 * 81.0 / 280.0
                } else {
                    0.0
                },
                2e-14,
                1e-16,
            );
        }
    }
    for i in 0..n {
        for j in 0..n {
            assert_eq!(
                entry(enriched.stiffness(), i, j),
                entry(&model.pencil().k, i, j)
            );
        }
    }
    let (translation, rotation) = (0.7, -0.2);
    let mut velocity = vec![0.0; 18];
    for node in 0..3 {
        velocity[6 * node + 1] = translation;
    }
    let state = enriched
        .state(
            0.0,
            &[0.0; 18],
            &velocity,
            &[[0.0; 2]],
            &[[rotation, 0.0]],
            || false,
        )
        .unwrap();
    // theta_x cross e_z = -e_y. The I1 cross term is essential.
    let expected = 0.5
        * AREA
        * (I0 * translation * translation - 2.0 * I1 * translation * rotation * 9.0 / 20.0
            + I2 * rotation * rotation * 81.0 / 280.0);
    near(state.kinetic_energy_j(), expected, 4e-14, 1e-15);
    let mut drilling = [0.0; 18];
    for node in 0..3 {
        drilling[6 * node + 5] = 0.3;
    }
    // A director-axis rotation has no physical material velocity. The native
    // transient gauge retains potential stiffness but must not manufacture
    // inertia. The general modal pencil separately regularizes relative-spin
    // inertia; the physical trajectory deliberately does not inherit it.
    assert!(entry(&model.pencil().m, 5, 5) > 0.0);
    near(entry(enriched.mass(), 5, 5), 0.0, 0.0, 1e-16);
    assert!(
        enriched
            .state(0.0, &[0.0; 18], &drilling, &[[0.0; 2]], &[[0.0; 2]], || {
                false
            })
            .is_err()
    );
}

#[test]
fn g1_heated_internal_rotations_obey_momentum_instead_of_jumping_to_static_lift() {
    let model = make_model(true, ShellMass::Consistent);
    let enriched = EnrichedMitcShell::try_new(Arc::clone(&model), || false).unwrap();
    let cold = model.prepare_thermal(&fields(0.0), || false).unwrap();
    let hot = model.prepare_thermal(&fields(0.002), || false).unwrap();
    let initial = enriched
        .equilibrium_state(0.0, &cold, &[0.0; 18], || false)
        .unwrap();
    let static_hot = enriched
        .equilibrium_state(0.0, &hot, &[0.0; 18], || false)
        .unwrap();
    let load0 = enriched
        .generalized_load(&cold, &[0.0; 18], || false)
        .unwrap();
    let load1 = enriched
        .generalized_load(&hot, &[0.0; 18], || false)
        .unwrap();
    let dt = 1e-4;
    let prepared = enriched.prepare_midpoint(dt, || false).unwrap();
    let step = prepared
        .advance(&initial, &load0, &load1, || false)
        .unwrap();
    let n = model.pencil().free;
    assert_eq!(n, 1);
    let mb = AREA * I2 * 81.0 / 280.0;
    let h = 0.5 * dt;
    let aa = mb + h * h * entry(enriched.stiffness(), n, n);
    let bb = h * h * entry(enriched.stiffness(), n, n + 1);
    let cc = mb + h * h * entry(enriched.stiffness(), n + 1, n + 1);
    let rhs = [0.5 * h * load1[n], 0.5 * h * load1[n + 1]];
    let determinant = aa * cc - bb * bb;
    let vmid = [
        (cc * rhs[0] - bb * rhs[1]) / determinant,
        (aa * rhs[1] - bb * rhs[0]) / determinant,
    ];
    for i in 0..2 {
        near(
            step.state.bubble_coordinates()[0][i],
            dt * vmid[i],
            2e-12,
            1e-16,
        );
        near(
            step.state.bubble_velocity()[0][i],
            2.0 * vmid[i],
            2e-12,
            1e-16,
        );
    }
    let actual = step.state.bubble_coordinates()[0]
        .iter()
        .map(|x| x * x)
        .sum::<f64>();
    let quasistatic = static_hot.bubble_coordinates()[0]
        .iter()
        .map(|x| x * x)
        .sum::<f64>();
    assert!(quasistatic > 1e-16);
    assert!(actual < 1e-5 * quasistatic);
    assert!(step.state.kinetic_energy_j() > 0.0);
    near(step.report.energy_residual_j, 0.0, 0.0, 1e-17);

    let response = enriched.recover(&hot, &step.state, || false).unwrap();
    assert!(response.maximum_director_motion_rad > 0.0);
    let point = response.strain_fields[0]
        .at_barycentric([0.2, 0.3, 0.5])
        .unwrap();
    assert!(point.gradient_per_m.iter().any(|x| x.abs() > 0.0));
    // Exact generalized material residual, including nonzero bubble force.
    for row in 0..enriched.dimension() {
        let (cols, vals) = enriched.stiffness().row(row);
        let kq: f64 = cols
            .iter()
            .zip(vals)
            .map(|(&c, &v)| v * step.state.coordinates()[c])
            .sum();
        near(
            response.material_generalized_force[row] + response.numerical_gauge_force[row],
            kq - load1[row],
            5e-12,
            1e-12,
        );
    }
    assert!(
        response.bubble_material_force[0]
            .iter()
            .any(|x| x.abs() > 1e-8)
    );
}

#[test]
fn g1_native_midpoint_conserves_unforced_full_energy_and_retains_static_limit() {
    let model = make_model(true, ShellMass::Consistent);
    let enriched = EnrichedMitcShell::try_new(Arc::clone(&model), || false).unwrap();
    let cold = model.prepare_thermal(&fields(0.0), || false).unwrap();
    let hot = model.prepare_thermal(&fields(0.001), || false).unwrap();
    let hot_state = enriched
        .equilibrium_state(0.0, &hot, &[0.0; 18], || false)
        .unwrap();
    let hot_load = enriched
        .generalized_load(&hot, &[0.0; 18], || false)
        .unwrap();
    let prepared = enriched.prepare_midpoint(0.002, || false).unwrap();
    let static_step = prepared
        .advance(&hot_state, &hot_load, &hot_load, || false)
        .unwrap();
    for (a, b) in static_step
        .state
        .coordinates()
        .iter()
        .zip(hot_state.coordinates())
    {
        near(*a, *b, 2e-11, 1e-14);
    }
    let condensed = hot.recover(&[0.0; 18], || false).unwrap();
    let full = enriched.recover(&hot, &hot_state, || false).unwrap();
    near(
        full.response.strain_energy_j,
        condensed.strain_energy_j,
        3e-13,
        1e-15,
    );
    for f in full.bubble_material_force.iter().flatten() {
        near(*f, 0.0, 0.0, 1e-12);
    }

    let initial = enriched
        .state(
            0.0,
            &[0.0; 18],
            &[0.0; 18],
            &[[0.002, -0.001]],
            &[[0.0; 2]],
            || false,
        )
        .unwrap();
    let reference = enriched.recover(&cold, &initial, || false).unwrap();
    let energy = initial.kinetic_energy_j()
        + reference.response.strain_energy_j
        + reference.drilling_gauge_energy_j;
    let load = vec![0.0; enriched.dimension()];
    let mut state = initial;
    for _ in 0..80 {
        let step = prepared.advance(&state, &load, &load, || false).unwrap();
        near(step.report.energy_residual_j, 0.0, 0.0, 1e-14);
        state = step.state;
    }
    let response = enriched.recover(&cold, &state, || false).unwrap();
    near(
        state.kinetic_energy_j()
            + response.response.strain_energy_j
            + response.drilling_gauge_energy_j,
        energy,
        2e-10,
        1e-14,
    );
    assert!(state.kinetic_energy_j() > 0.0);
}

#[test]
fn g0_refused_or_cancelled_proposals_leave_original_motion_and_factor_usable() {
    let model = make_model(true, ShellMass::Consistent);
    let enriched = EnrichedMitcShell::try_new(Arc::clone(&model), || false).unwrap();
    let thermal = model.prepare_thermal(&fields(0.0), || false).unwrap();
    let state = enriched
        .equilibrium_state(0.0, &thermal, &[0.0; 18], || false)
        .unwrap();
    let prepared = enriched.prepare_midpoint(0.01, || false).unwrap();
    let mut force = vec![0.0; enriched.dimension()];
    force[1] = 0.01;
    let expected = prepared.advance(&state, &force, &force, || false).unwrap();
    assert!(matches!(
        prepared.advance(&state, &force, &force, || true),
        Err(PlateError::Cancelled)
    ));
    let mut count = 0;
    prepared
        .advance(&state, &force, &force, || {
            count += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        prepared.advance(&state, &force, &force, || {
            calls += 1;
            calls == count
        }),
        Err(PlateError::Cancelled)
    ));
    assert_eq!(state.time_s(), 0.0);
    assert!(state.coordinates().iter().all(|x| *x == 0.0));
    let repeated = prepared.advance(&state, &force, &force, || false).unwrap();
    assert_eq!(repeated.state.coordinates(), expected.state.coordinates());
    assert_eq!(repeated.state.velocity(), expected.state.velocity());
    let mut forbidden = [0.0; 18];
    forbidden[0] = 0.001;
    assert!(
        enriched
            .state(
                0.0,
                &[0.0; 18],
                &forbidden,
                &[[0.0; 2]],
                &[[0.0; 2]],
                || false
            )
            .is_err()
    );
    assert!(enriched.prepare_midpoint(0.0, || false).is_err());
    assert!(
        prepared
            .advance(&state, &force[..2], &force, || false)
            .is_err()
    );
    let mut gauge_force = force.clone();
    gauge_force[0] = 0.01;
    assert!(
        prepared
            .advance(&state, &gauge_force, &gauge_force, || false)
            .is_err()
    );
    let other = make_model(true, ShellMass::Consistent);
    let foreign = other.prepare_thermal(&fields(0.0), || false).unwrap();
    assert!(
        enriched
            .generalized_load(&foreign, &[0.0; 18], || false)
            .is_err()
    );
    assert!(EnrichedMitcShell::try_new(make_model(false, ShellMass::Lumped), || false).is_err());
}

#[test]
fn g3_rotating_geometry_directors_and_oblique_supports_preserves_enriched_motion() {
    let build = |frame: [[f64; 3]; 3]| {
        let rotate = |point: [f64; 3]| {
            core::array::from_fn(|axis| (0..3).map(|i| point[i] * frame[i][axis]).sum())
        };
        let mesh = ShellMesh::new(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
                .map(rotate)
                .to_vec(),
            vec![[0, 1, 2]],
        )
        .unwrap();
        let directors = NodalDirectors::new(&mesh, vec![frame[2]; 3], 0.9).unwrap();
        let section = ShearShellSection::new_coupled(
            SheetMatrix::orthotropic(2.0e4, 1.0e4, 0.2, 3.0e3).unwrap(),
            [0.0; 9],
            SheetMatrix::orthotropic(20.0, 10.0, 0.2, 3.0).unwrap(),
            [120.0, 0.0, 80.0],
            1.0,
            I0,
            I1,
            I2,
            0.0,
        )
        .unwrap();
        let supports: Vec<_> = [
            [true, true, true, false, false, true],
            [false, true, false, false, false, false],
            [false, false, false, false, false, false],
        ]
        .into_iter()
        .enumerate()
        .map(|(node, mask)| ShellConstraint {
            node,
            mask,
            frame: ConstraintFrame::Local(frame),
        })
        .collect();
        Arc::new(
            assemble_mitc3p(
                &mesh,
                &directors,
                &[section],
                &supports,
                ShellMass::Consistent,
                MitcBudget {
                    max_nodes: 3,
                    max_triangles: 1,
                },
                1e-3,
                || false,
            )
            .unwrap(),
        )
    };
    let original = build([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    let rotated = build([[0.6, 0.8, 0.0], [-0.64, 0.48, 0.6], [0.48, -0.36, 0.8]]);
    let a = EnrichedMitcShell::try_new(Arc::clone(&original), || false).unwrap();
    let b = EnrichedMitcShell::try_new(Arc::clone(&rotated), || false).unwrap();
    assert_eq!(a.dimension(), b.dimension());
    let free = original.pencil().free;
    let q: Vec<_> = (0..free).map(|i| 1e-5 * (i as f64 - 5.0)).collect();
    let v: Vec<_> = (0..free).map(|i| 2e-4 * (3.0 - i as f64)).collect();
    let compatible = |reduced: &[f64]| {
        let mut world = original.expand(reduced).unwrap();
        let spin = 0.5 * (world[7] - world[1] - world[12] + world[0]);
        for node in [1, 2] {
            world[6 * node + 5] = spin;
        }
        original.project_load(&world).unwrap()
    };
    let q = compatible(&q);
    let v = compatible(&v);
    let first = a
        .state(
            0.0,
            &original.expand(&q).unwrap(),
            &original.expand(&v).unwrap(),
            &[[0.0002, -0.0001]],
            &[[0.003, -0.002]],
            || false,
        )
        .unwrap();
    let second = b
        .state(
            0.0,
            &rotated.expand(&q).unwrap(),
            &rotated.expand(&v).unwrap(),
            &[[0.0002, -0.0001]],
            &[[0.003, -0.002]],
            || false,
        )
        .unwrap();
    near(
        first.kinetic_energy_j(),
        second.kinetic_energy_j(),
        4e-12,
        1e-16,
    );
    let heat_a = original.prepare_thermal(&fields(0.0007), || false).unwrap();
    let heat_b = rotated.prepare_thermal(&fields(0.0007), || false).unwrap();
    let load_a = a.generalized_load(&heat_a, &[0.0; 18], || false).unwrap();
    let load_b = b.generalized_load(&heat_b, &[0.0; 18], || false).unwrap();
    let next_a = a
        .prepare_midpoint(0.001, || false)
        .unwrap()
        .advance(&first, &load_a, &load_a, || false)
        .unwrap();
    let next_b = b
        .prepare_midpoint(0.001, || false)
        .unwrap()
        .advance(&second, &load_b, &load_b, || false)
        .unwrap();
    for (a, b) in next_a
        .state
        .coordinates()
        .iter()
        .zip(next_b.state.coordinates())
    {
        near(*a, *b, 2e-9, 1e-14);
    }
    for (a, b) in next_a.state.velocity().iter().zip(next_b.state.velocity()) {
        near(*a, *b, 2e-9, 1e-12);
    }
    let response_a = a.recover(&heat_a, &next_a.state, || false).unwrap();
    let response_b = b.recover(&heat_b, &next_b.state, || false).unwrap();
    near(
        response_a.response.strain_energy_j,
        response_b.response.strain_energy_j,
        2e-9,
        1e-14,
    );
}
