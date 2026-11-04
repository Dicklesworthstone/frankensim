//! Independent physical checks for the TET10/P1 finite-strain mixed owner.
//!
//! The affine cases check exact constitutive response and pressure, not locking.
//! Non-affine Cook-panel refinement is a separate, explicitly bounded diagnostic.

use std::collections::BTreeMap;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::mixed::{
    MixedLawField, MixedUpError, MixedUpEvaluation, MixedUpProblem, MixedUpSettings,
};
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solid::linear3::{TetAssemblyBudget, TetElasticError};
use fs_solver::nonlinear::NewtonKrylovConfig;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    with_gate(&gate, f);
}

fn with_gate(gate: &CancelGate, f: impl FnOnce(&Cx<'_>)) {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 61,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn law(mu: f64, bulk: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 0.5 * mu,
            c01: 0.0,
            kappa: bulk,
        },
        5.0,
    )
    .unwrap()
}

fn settings(length: f64, stress: f64) -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 4,
        maximum_newton_iterations: 30,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: length,
        stress_scale_pa: stress,
        maximum_assembly_bytes: 64 * 1024 * 1024,
        maximum_factor_bytes: 32 * 1024 * 1024,
    }
}

fn problem<'a, L: MixedFiniteStrainPotential>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: &'a L,
    supports: &'a [(usize, f64)],
    forces: &'a [(usize, f64)],
    volumetric: VolumetricTreatment,
) -> MixedUpProblem<'a, L> {
    MixedUpProblem {
        nodes_m: nodes,
        elements,
        material: MixedLawField::Uniform(material),
        volumetric,
        prescribed_m: supports,
        nodal_forces_n: forces,
        body_force_n_m3: [0.0; 3],
        pressure_gauge: None,
        budget: TetAssemblyBudget::standard(),
        quadrature_points_per_dir: 4,
        maximum_quadrature_points: 20_000,
        settings: settings(1.0, 1.0),
    }
}

// Conforming Kuhn topology only; the production promotion owns every P2 node.
fn block(n: [usize; 3], size: [f64; 3]) -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let id = |i, j, k| (k * (n[1] + 1) + j) * (n[0] + 1) + i;
    let mut nodes = Vec::new();
    for k in 0..=n[2] {
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                nodes.push([
                    size[0] * i as f64 / n[0] as f64,
                    size[1] * j as f64 / n[1] as f64,
                    size[2] * k as f64 / n[2] as f64,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                for perm in [
                    [0, 1, 2],
                    [0, 2, 1],
                    [1, 0, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [2, 1, 0],
                ] {
                    let mut cell = [i, j, k];
                    let mut tet = [id(i, j, k); 4];
                    for (step, axis) in perm.into_iter().enumerate() {
                        cell[axis] += 1;
                        tet[step + 1] = id(cell[0], cell[1], cell[2]);
                    }
                    tets.push(tet);
                }
            }
        }
    }
    (nodes, tets)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Exact constant traction on straight quadratic triangle faces: corner
/// integrals vanish; each of the three midside integrals is A/3.
fn face_load(
    nodes: &[[f64; 3]],
    elements: &[[usize; 10]],
    selected: impl Fn([f64; 3]) -> bool,
    traction: [f64; 3],
) -> Vec<(usize, f64)> {
    let mut faces = BTreeMap::<[usize; 3], (usize, [usize; 3])>::new();
    for tet in elements {
        for omit in 0..4 {
            let local: Vec<_> = (0..4).filter(|&a| a != omit).collect();
            let mut key = [tet[local[0]], tet[local[1]], tet[local[2]]];
            key.sort_unstable();
            let midsides: Vec<_> = TET10_EDGES
                .iter()
                .enumerate()
                .filter(|(_, edge)| !edge.contains(&omit))
                .map(|(edge, _)| tet[4 + edge])
                .collect();
            let entry = faces
                .entry(key)
                .or_insert((0, [midsides[0], midsides[1], midsides[2]]));
            entry.0 += 1;
        }
    }
    let mut force = BTreeMap::<usize, f64>::new();
    for (corners, (count, midsides)) in faces {
        if count != 1 || !corners.into_iter().all(|node| selected(nodes[node])) {
            continue;
        }
        let [a, b, c] = corners.map(|node| nodes[node]);
        let normal = cross(sub(b, a), sub(c, a));
        let area = 0.5 * normal[0].hypot(normal[1]).hypot(normal[2]);
        for node in midsides {
            for axis in 0..3 {
                *force.entry(3 * node + axis).or_default() += area * traction[axis] / 3.0;
            }
        }
    }
    force
        .into_iter()
        .filter(|(_, value)| *value != 0.0)
        .collect()
}

fn rollers(nodes: &[[f64; 3]]) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, point)| {
            (0..3).filter_map(move |axis| {
                (point[axis].abs() < 1e-12).then_some((3 * node + axis, 0.0))
            })
        })
        .collect()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn exact_incompressible_traction_recovers_uniaxial_stretch_pressure_and_reaction() {
    with_cx(|cx| {
        let mu = 2.5;
        let material = law(mu, 40.0 * mu);
        let (p1, tets) = block([2, 1, 1], [2.0, 1.0, 1.0]);
        let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let supports = rollers(&nodes);
        for stretch in [1.2_f64, 1.5] {
            let lateral = stretch.sqrt().recip();
            let nominal = mu * (stretch - stretch.powi(-2));
            let expected_p = mu * (stretch * stretch - stretch.recip()) / 3.0;
            let forces = face_load(
                &nodes,
                &elements,
                |p| (p[0] - 2.0).abs() < 1e-12,
                [nominal, 0.0, 0.0],
            );
            near(forces.iter().map(|(_, value)| value).sum(), nominal, 2e-14);
            let mut p = problem(
                &nodes,
                &elements,
                &material,
                &supports,
                &forces,
                VolumetricTreatment::ExactIncompressible,
            );
            p.settings = settings(2.0, mu);
            let result = p.solve(cx).unwrap();
            for (point, u) in nodes.iter().zip(&result.displacement_m) {
                for axis in 0..3 {
                    let scale = if axis == 0 { stretch } else { lateral };
                    near(u[axis], (scale - 1.0) * point[axis], 3e-9);
                }
            }
            for pressure in &result.pressure_pa {
                near(*pressure, expected_p, 2e-8 * mu);
            }
            near(result.min_det_f, 1.0, 3e-10);
            let reaction: f64 = result
                .reactions_n
                .iter()
                .filter(|(d, _)| d % 3 == 0)
                .map(|(_, v)| v)
                .sum();
            near(reaction, -nominal, 2e-9 * mu);
            assert_eq!(
                result.kkt_inertia.positive,
                3 * nodes.len() - supports.len()
            );
            assert_eq!(result.kkt_inertia.negative, result.pressure_pa.len());
        }
    });
}

// Scalar lateral traction equation from the declared decoupled NH energy,
// independent of the material stress and of the mixed implementation.
fn finite_bulk_uniaxial(mu: f64, bulk: f64, stretch: f64) -> (f64, f64, f64) {
    let (mut lo, mut hi) = (stretch.sqrt().recip(), 1.0);
    for _ in 0..80 {
        let lateral = f64::midpoint(lo, hi);
        let j = stretch * lateral * lateral;
        let residual = bulk * (j - 1.0) * j
            - mu / 3.0 * j.powf(-2.0 / 3.0) * (stretch * stretch - lateral * lateral);
        if residual > 0.0 {
            hi = lateral;
        } else {
            lo = lateral;
        }
    }
    let lateral = f64::midpoint(lo, hi);
    let j = stretch * lateral * lateral;
    let pressure = mu / (3.0 * j.powf(5.0 / 3.0)) * (stretch * stretch - lateral * lateral);
    let nominal = mu * j.powf(-2.0 / 3.0) * (stretch - lateral * lateral / stretch);
    (lateral, pressure, nominal)
}

#[test]
fn finite_bulk_traction_preserves_original_law_over_eight_bulk_decades() {
    with_cx(|cx| {
        let (p1, tets) = block([1, 1, 1], [1.0; 3]);
        let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let supports = rollers(&nodes);
        let mu = 2.0;
        let stretch = 1.2;
        for ratio in [1.0, 1e2, 1e4, 1e6, 1e8] {
            let material = law(mu, mu * ratio);
            let (lateral, pressure, nominal) = finite_bulk_uniaxial(mu, mu * ratio, stretch);
            let forces = face_load(
                &nodes,
                &elements,
                |p| (p[0] - 1.0).abs() < 1e-12,
                [nominal, 0.0, 0.0],
            );
            let mut p = problem(
                &nodes,
                &elements,
                &material,
                &supports,
                &forces,
                VolumetricTreatment::OriginalMaterial,
            );
            p.settings = settings(1.0, mu);
            let result = p.solve(cx).unwrap();
            let mut error = 0.0_f64;
            for (point, u) in nodes.iter().zip(&result.displacement_m) {
                for axis in 0..3 {
                    let scale = if axis == 0 { stretch } else { lateral };
                    error = error.max((u[axis] - (scale - 1.0) * point[axis]).abs());
                }
            }
            assert!(
                error < 4e-9,
                "bulk/shear={ratio:e}, displacement error={error:e}"
            );
            for value in &result.pressure_pa {
                near(*value, pressure, 3e-8 * mu);
            }
            near(result.min_det_f, stretch * lateral * lateral, 1e-9);
            println!("mixed bulk/shear={ratio:e} uniaxial displacement_error={error:e}");
        }
    });
}

fn residual(evaluation: &MixedUpEvaluation) -> Vec<f64> {
    evaluation
        .internal_forces_n
        .iter()
        .chain(&evaluation.pressure_residual_m3)
        .copied()
        .collect()
}

#[test]
fn full_mixed_tangent_matches_finite_differences_and_superposed_rotation() {
    with_cx(|cx| {
        let p1 = [
            [0.0, 0.0, 0.0],
            [1.1, 0.1, 0.0],
            [0.2, 1.0, 0.1],
            [0.0, 0.2, 0.9],
        ];
        let (nodes, elements) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let material = law(2.5, 80.0);
        let p = problem(
            &nodes,
            &elements,
            &material,
            &[],
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        let prepared = p.prepare(cx).unwrap();
        let u: Vec<_> = nodes
            .iter()
            .map(|x| {
                [
                    0.08 * x[0] * x[0] + 0.03 * x[1],
                    -0.04 * x[1] + 0.02 * x[0] * x[2],
                    0.05 * x[2],
                ]
            })
            .collect();
        let pressure: Vec<_> = prepared
            .pressure_nodes()
            .iter()
            .map(|&node| 0.3 + 0.2 * nodes[node][0])
            .collect();
        let evaluated = prepared.evaluate(&u, &pressure, true, cx).unwrap();
        let tangent = evaluated.tangent.as_ref().unwrap();
        let count = 3 * nodes.len() + pressure.len();
        let direction: Vec<_> = (0..count)
            .map(|i| ((i * 7 % 13) as f64 - 6.0) / 17.0)
            .collect();
        let h = 1e-6;
        let shifted = |sign: f64| {
            let moved: Vec<[f64; 3]> = u
                .iter()
                .enumerate()
                .map(|(n, v)| {
                    core::array::from_fn(|axis| v[axis] + sign * h * direction[3 * n + axis])
                })
                .collect();
            let shifted_p: Vec<_> = pressure
                .iter()
                .enumerate()
                .map(|(i, p)| p + sign * h * direction[3 * nodes.len() + i])
                .collect();
            residual(&prepared.evaluate(&moved, &shifted_p, false, cx).unwrap())
        };
        let plus = shifted(1.0);
        let minus = shifted(-1.0);
        let mut scale = 0.0_f64;
        let mut error = 0.0_f64;
        for row in 0..count {
            let (columns, values) = tangent.row(row);
            let action: f64 = columns
                .iter()
                .zip(values)
                .map(|(&col, &value)| value * direction[col])
                .sum();
            let fd = (plus[row] - minus[row]) / (2.0 * h);
            scale = scale.max(fd.abs());
            error = error.max((action - fd).abs());
            for col in 0..count {
                near(tangent.get(row, col), tangent.get(col, row), 2e-12);
            }
        }
        assert!(
            error < 3e-7 * scale,
            "mixed tangent relative error={:e}",
            error / scale
        );
        let rotate = |v: [f64; 3]| [-v[1], v[0], v[2]];
        let rotated: Vec<_> = nodes
            .iter()
            .zip(&u)
            .map(|(x, u)| {
                let moved = core::array::from_fn(|i| x[i] + u[i]);
                let q = rotate(moved);
                core::array::from_fn(|i| q[i] - x[i])
            })
            .collect();
        let transformed = prepared.evaluate(&rotated, &pressure, false, cx).unwrap();
        near(transformed.potential_j, evaluated.potential_j, 2e-12);
        near(transformed.min_det_f, evaluated.min_det_f, 2e-13);
        for node in 0..nodes.len() {
            let original = core::array::from_fn(|i| evaluated.internal_forces_n[3 * node + i]);
            let expected = rotate(original);
            for axis in 0..3 {
                near(
                    transformed.internal_forces_n[3 * node + axis],
                    expected[axis],
                    3e-12,
                );
            }
        }
        for (a, b) in transformed
            .pressure_residual_m3
            .iter()
            .zip(&evaluated.pressure_residual_m3)
        {
            near(*a, *b, 2e-14);
        }
    });
}

#[test]
fn pressure_constraint_integrates_quartic_p2_kinematics_without_underintegration() {
    with_cx(|cx| {
        let corners = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (nodes, elements) = promote_p1_to_p2(&corners, &[[0, 1, 2, 3]], None, cx).unwrap();
        let material = law(1.0, 100.0);
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &[],
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        let [a, b, c] = [0.2, 0.3, -0.1];
        let u: Vec<_> = nodes
            .iter()
            .map(|x| {
                [
                    0.5 * a * x[0] * x[0],
                    0.5 * b * x[1] * x[1],
                    0.5 * c * x[2] * x[2],
                ]
            })
            .collect();
        // ∫ λ0 [(1+ax)(1+by)(1+cz)−1] dV, by exact simplex monomial integrals.
        let exact = (a + b + c) / 120.0 + (a * b + a * c + b * c) / 720.0 + a * b * c / 5040.0;
        for order in [4, 5] {
            p.quadrature_points_per_dir = order;
            let prepared = p.prepare(cx).unwrap();
            near(prepared.reference_volume_m3(), 1.0 / 6.0, 2e-15);
            let evaluated = prepared.evaluate(&u, &[0.0; 4], false, cx).unwrap();
            near(evaluated.pressure_residual_m3[0], exact, 3e-16);
        }
        p.quadrature_points_per_dir = 3;
        assert!(
            p.prepare(cx).is_err(),
            "cubic J times P1 pressure needs a full quartic rule"
        );
    });
}

#[test]
fn nonaffine_cook_panel_refines_and_stays_finite_at_nearly_incompressible_bulk() {
    with_cx(|cx| {
        let mut tips = Vec::new();
        for n in [2, 4, 6] {
            let (mut p1, tets) = block([n, n, 1], [1.0; 3]);
            for point in &mut p1 {
                let [s, t, z] = *point;
                *point = [4.8 * s, 4.4 * s + (4.4 - 2.8 * s) * t, 0.1 * z];
            }
            let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
            let material = law(1.0, 1e6);
            let supports: Vec<_> = nodes
                .iter()
                .enumerate()
                .flat_map(|(node, point)| {
                    (0..3).filter_map(move |axis| {
                        (axis == 2 || point[0].abs() < 1e-12).then_some((3 * node + axis, 0.0))
                    })
                })
                .collect();
            let forces = face_load(
                &nodes,
                &elements,
                |p| (p[0] - 4.8).abs() < 1e-12,
                [0.0, 1e-4, 0.0],
            );
            near(forces.iter().map(|(_, value)| value).sum(), 1.6e-5, 2e-18);
            let tip = nodes
                .iter()
                .position(|p| {
                    (p[0] - 4.8).abs() < 1e-12 && (p[1] - 5.2).abs() < 1e-12 && p[2].abs() < 1e-12
                })
                .unwrap();
            let mut p = problem(
                &nodes,
                &elements,
                &material,
                &supports,
                &forces,
                VolumetricTreatment::ExactIncompressible,
            );
            p.settings = settings(4.8, 1.0);
            p.settings.load_steps = 1;
            let exact = p.solve(cx).unwrap();
            let value = exact.displacement_m[tip][1];
            assert!(value > 0.0 && value < 0.02 && exact.min_det_f > 0.99);
            tips.push(value);
            p.volumetric = VolumetricTreatment::OriginalMaterial;
            let finite = p.solve(cx).unwrap();
            let bulk_difference = (finite.displacement_m[tip][1] - value).abs() / value;
            assert!(
                bulk_difference < 2e-4,
                "bulk penalty changes mixed bending by {bulk_difference:e}"
            );
            println!(
                "mixed Cook n={n} elements={} tip_m={value:e} bulk1e6_vs_exact={bulk_difference:e}",
                elements.len()
            );
        }
        let coarse_change = (tips[1] - tips[0]).abs();
        let fine_change = (tips[2] - tips[1]).abs();
        assert!(
            fine_change < coarse_change,
            "non-affine refinement is not resolving: {tips:?}"
        );
        assert!(
            fine_change / tips[2] < 0.1,
            "finest-pair diagnostic is still coarse: {tips:?}"
        );
        // This bounded ladder is not the bead's full MMS/inf-sup/Cook benchmark.
    });
}

#[test]
fn already_equilibrated_compressed_bar_cannot_bypass_kkt_inertia_refusal() {
    with_cx(|cx| {
        let (p1, tets) = block([4, 1, 1], [10.0, 1.0, 1.0]);
        let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let mu = 2.0;
        let material = law(mu, 100.0);
        let stretch = 0.7_f64;
        let lateral = stretch.sqrt().recip();
        let u: Vec<[f64; 3]> = nodes
            .iter()
            .map(|p| {
                [
                    (stretch - 1.0) * p[0],
                    (lateral - 1.0) * p[1],
                    (lateral - 1.0) * p[2],
                ]
            })
            .collect();
        let displacement = &u;
        let supports: Vec<_> = nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p[0].abs() < 1e-12 || (p[0] - 10.0).abs() < 1e-12)
            .flat_map(|(node, _)| {
                (0..3).map(move |axis| (3 * node + axis, displacement[node][axis]))
            })
            .collect();
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        p.settings = settings(10.0, mu);
        p.settings.load_steps = 1;
        let prepared = p.prepare(cx).unwrap();
        let pressure =
            vec![mu * (stretch * stretch - stretch.recip()) / 3.0; prepared.pressure_nodes().len()];
        let result = p.solve_from(&u, &pressure, cx);
        assert!(
            matches!(result, Err(MixedUpError::KktInertia { .. })),
            "a zero-residual compressed equilibrium must still pass stability admission: {result:?}"
        );
    });
}

#[test]
fn pressure_gauge_requires_a_true_constant_null_mode_on_free_displacements() {
    with_cx(|cx| {
        let material = law(2.0, 100.0);
        let (p1, tets) = block([2, 2, 2], [1.0; 3]);
        let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let boundary: Vec<_> = nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p.iter().any(|&x| x == 0.0 || x == 1.0))
            .flat_map(|(node, _)| (0..3).map(move |axis| (3 * node + axis, 0.0)))
            .collect();
        let free = 3 * nodes.len() - boundary.len();
        assert!(
            free > 0,
            "the gauge must be checked against interior P2 motion"
        );
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &boundary,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        p.settings.load_steps = 1;
        p.pressure_gauge = Some((0, 0.0));
        let result = p.solve(cx).unwrap();
        assert_eq!(result.kkt_inertia.positive, free);
        assert_eq!(result.kkt_inertia.negative, result.pressure_pa.len() - 1);
        assert!(
            result
                .displacement_m
                .iter()
                .flatten()
                .all(|u| u.abs() < 1e-12)
        );
        assert!(result.pressure_pa.iter().all(|p| p.abs() < 1e-12));
        near(result.min_det_f, 1.0, 1e-12);

        // Natural traction determines hydrostatic pressure; this same pin
        // would remove a physical pressure equation and must be refused.
        let supports = rollers(&nodes);
        let forces = face_load(&nodes, &elements, |p| p[0] == 1.0, [0.1, 0.0, 0.0]);
        p.prescribed_m = &supports;
        p.nodal_forces_n = &forces;
        assert!(matches!(
            p.solve(cx),
            Err(MixedUpError::IncompatiblePressureGauge { node: 0 })
        ));
    });
}

#[test]
fn quadrature_budget_constraints_and_precancelled_work_are_refused() {
    with_cx(|cx| {
        let material = law(2.0, 100.0);
        let (p1, tets) = block([1, 1, 1], [1.0; 3]);
        let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &[],
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        let required = 64 * elements.len();
        p.maximum_quadrature_points = required;
        assert!(
            p.prepare(cx).is_ok(),
            "an exactly full quadrature budget is valid"
        );
        p.maximum_quadrature_points = required - 1;
        assert!(matches!(p.prepare(cx),
            Err(MixedUpError::Geometry(TetElasticError::BudgetExceeded {
                what: "mixed_quadrature_points", requested, maximum
            })) if requested == required && maximum + 1 == required));
        p.maximum_quadrature_points = required;
        for constraints in [
            vec![(0, 0.0), (0, 0.0)],
            vec![(3 * nodes.len(), 0.0)],
            vec![(0, f64::NAN)],
        ] {
            let mut bad = p;
            bad.prescribed_m = &constraints;
            assert!(matches!(
                bad.prepare(cx),
                Err(MixedUpError::InvalidInput { .. })
            ));
        }
        let gate = CancelGate::new();
        gate.request();
        with_gate(&gate, |cancelled| {
            assert!(matches!(p.prepare(cancelled), Err(MixedUpError::Cancelled)));
            assert!(matches!(p.solve(cancelled), Err(MixedUpError::Cancelled)));
        });
    });
}

fn modal_settings() -> fs_solid::hyper3::mixed::vibration::MixedVibrationSettings {
    fs_solid::hyper3::mixed::vibration::MixedVibrationSettings {
        spectrum: fs_solid::hyper3::vibration::VibrationSettings {
            maximum_free_dofs: 128,
            relative_residual_tolerance: 1e-8,
        },
        maximum_mass_assembly_bytes: 8 * 1024 * 1024,
        maximum_condensation_bytes: 32 * 1024 * 1024,
        pressure_residual_tolerance: 1e-10,
    }
}

#[test]
fn finite_bulk_corner_vibration_has_exact_p2_mass_and_algebraic_pressure() {
    with_cx(|cx| {
        let (length, rho, mu) = (0.4_f64, 7.0, 2.0);
        let p1 = [
            [0.0; 3],
            [length, 0.0, 0.0],
            [0.0, length, 0.0],
            [0.0, 0.0, length],
        ];
        let (nodes, elements) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != 3)
            .map(|d| (d, 0.0))
            .collect();
        let volume = length.powi(3) / 6.0;
        // On an affine Tet10, N1=lambda1(2lambda1-1), so
        // integral N1²=V/70 and integral (dN1/dX1)²=3V/(5L²).
        // div(delta u) is P1 and the native pressure space represents it exactly.
        // These are homogeneous patch identities, not a general locking ladder.
        for ratio in [1.0, 1e4, 1e8] {
            let bulk = ratio * mu;
            let material = law(mu, bulk);
            let p = problem(
                &nodes,
                &elements,
                &material,
                &fixed,
                &[],
                VolumetricTreatment::OriginalMaterial,
            );
            let spectrum = p.solve_vibration(&[rho], modal_settings(), cx).unwrap();
            assert_eq!(spectrum.free_dofs(), &[3]);
            assert_eq!(spectrum.mass_kg().len(), 1, "pressure must have no inertia");
            assert_eq!(spectrum.equilibrium().pressure_nodes, [0, 1, 2, 3]);
            let stiffness = (bulk + 4.0 * mu / 3.0) * 3.0 * volume / (5.0 * length * length);
            let mass = rho * volume / 70.0;
            near(spectrum.stiffness_n_m()[0], stiffness, 2e-11 * stiffness);
            near(spectrum.mass_kg()[0], mass, 2e-13 * mass);
            let mode = &spectrum.modes()[0];
            near(
                mode.displacement.lambda,
                42.0 * (bulk + 4.0 * mu / 3.0) / (rho * length * length),
                2e-11 * stiffness / mass,
            );
            near(mode.displacement.phi[0].powi(2) * mass, 1.0, 2e-12);
            for (point, coefficient) in [-1.0, 3.0, -1.0, -1.0].into_iter().enumerate() {
                let expected = coefficient * bulk / length;
                near(
                    spectrum.pressure_response_pa_per_m()[point],
                    expected,
                    2e-11 * expected.abs(),
                );
                near(
                    mode.pressure_shape[point] / mode.displacement.phi[0],
                    expected,
                    2e-11 * expected.abs(),
                );
            }
            assert!(
                mode.pressure_equation_relative_residual
                    <= modal_settings().pressure_residual_tolerance
            );
            assert!(mode.displacement.interval.0 > 0.0);
        }
    });
}

fn modal_source_card(mu: f64, bulk: f64) -> fs_matdb::ConstitutiveModelCard {
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
    use fs_qty::Dims;
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let artifact = fs_blake3::hash_bytes(b"synthetic mixed modal material; not measured data");
    ConstitutiveModelCard {
        law: LawId("mooney-rivlin-decoupled".into()),
        law_version: 1,
        parameters: [
            ("stretch_limit", 5.0, Dims::NONE),
            ("c10", 0.5 * mu, pa),
            ("c01", 0.0, pa),
            ("kappa", bulk, pa),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: fs_evidence::ValidityDomain::unconstrained(),
        sources: vec![artifact],
        provenance: Provenance {
            source: "manufactured finite-compliance fixture".into(),
            license: "test-only".into(),
            artifact: Some(artifact),
        },
    }
}

#[test]
fn sourced_prestretched_mixed_modes_match_full_material_tangent_and_pressure_equations() {
    use fs_material::hyper::cards::ResolvedHyperelastic;
    with_cx(|cx| {
        let (mu, bulk, rho, length, stretch, shear) = (2.0, 30.0, 7.0, 0.4_f64, 1.2, 0.25);
        let p1 = [
            [0.0; 3],
            [length, 0.0, 0.0],
            [0.0, length, 0.0],
            [0.0, 0.0, length],
        ];
        let (nodes, elements) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let initial: Vec<_> = nodes
            .iter()
            .map(|x| [(stretch - 1.0) * x[0] + shear * x[1], 0.0, 0.0])
            .collect();
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|d| !(3..6).contains(d))
            .map(|d| (d, initial[d / 3][d % 3]))
            .collect();
        let card = modal_source_card(mu, bulk);
        let material = ResolvedHyperelastic::try_from_model_card(&card).unwrap();
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &fixed,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        p.settings.load_steps = 1;
        let pressure = [bulk * (stretch - 1.0); 4];
        // Constant stress integrates to zero force at this P2 corner. The
        // other nodes hold the homogeneous finite stretch/shear equilibrium.
        let spectrum = p
            .solve_vibration_from(&initial, &pressure, &[rho], modal_settings(), cx)
            .unwrap();
        assert_eq!(spectrum.free_dofs(), &[3, 4, 5]);
        let MixedLawField::Uniform(retained) = spectrum.problem().material else {
            panic!("original sourced material was replaced")
        };
        assert!(core::ptr::eq(retained, &material));
        assert_eq!(retained.model_card(), &card);
        assert_eq!(spectrum.reference_density_kg_m3(), &[rho]);
        let f = [stretch, shear, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let tangent = material.law().tangent(&f).unwrap();
        let volume = length.powi(3) / 6.0;
        let factor = 3.0 * volume / (5.0 * length * length);
        let mut expected_k = vec![0.0; 9];
        let mut expected_m = vec![0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                expected_k[3 * i + j] = factor * tangent[3 * i][3 * j];
                near(
                    spectrum.stiffness_n_m()[3 * i + j],
                    expected_k[3 * i + j],
                    2e-8 * factor * bulk,
                );
                expected_m[3 * i + j] = if i == j { rho * volume / 70.0 } else { 0.0 };
                near(
                    spectrum.mass_kg()[3 * i + j],
                    expected_m[3 * i + j],
                    2e-13 * rho * volume,
                );
            }
        }
        assert!(
            expected_k[1].abs() > factor * bulk * shear * 0.5,
            "finite shear must exercise pressure cross-coupling"
        );
        let expected_modes = fs_modal::eigh_gen_dense(&expected_k, &expected_m, 3).unwrap();
        let full = spectrum.mixed_tangent();
        let nu = 3 * nodes.len();
        for (mode, expected) in spectrum.modes().iter().zip(expected_modes) {
            near(
                mode.displacement.lambda,
                expected.lambda,
                3e-8 * expected.lambda.abs(),
            );
            // Verify the uncondensed physical pressure equation independently
            // of the adapter's reported residual and without any pressure mass.
            for row in 0..4 {
                let mut residual = 0.0_f64;
                let mut scale = 0.0_f64;
                for j in 0..3 {
                    let term = full.get(nu + row, 3 + j) * mode.displacement.phi[j];
                    residual += term;
                    scale += term.abs();
                }
                for j in 0..4 {
                    let term = full.get(nu + row, nu + j) * mode.pressure_shape[j];
                    residual += term;
                    scale += term.abs();
                }
                assert!(residual.abs() <= 1e-10 * scale.max(f64::MIN_POSITIVE));
            }
        }
        for (point, coefficient) in [-1.0, 3.0, -1.0, -1.0].into_iter().enumerate() {
            for (axis, cofactor) in [1.0, -shear, 0.0].into_iter().enumerate() {
                near(
                    spectrum.pressure_response_pa_per_m()[3 * point + axis],
                    coefficient * bulk * cofactor / length,
                    2e-8 * bulk / length,
                );
            }
        }
        let doubled = p
            .solve_vibration_from(&initial, &pressure, &[2.0 * rho], modal_settings(), cx)
            .unwrap();
        for (mode, heavier) in spectrum.modes().iter().zip(doubled.modes()) {
            near(
                heavier.displacement.lambda,
                0.5 * mode.displacement.lambda,
                2e-11 * mode.displacement.lambda.abs(),
            );
        }
    });
}

#[test]
fn rotated_isochoric_modes_keep_uncancelled_pressure_residual_scale() {
    with_cx(|cx| {
        let p1 = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (nodes, elements) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let initial: Vec<_> = nodes
            .iter()
            .map(|x| [-0.2 * x[0] - 0.6 * x[1], 0.6 * x[0] - 0.2 * x[1], 0.0])
            .collect();
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|d| !(3..6).contains(d))
            .map(|d| (d, initial[d / 3][d % 3]))
            .collect();
        let (mu, bulk, rho) = (2.0, 200.0, 7.0);
        let material = law(mu, bulk);
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &fixed,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        p.settings.load_steps = 1;
        let spectrum = p
            .solve_vibration_from(&initial, &[0.0; 4], &[rho], modal_settings(), cx)
            .unwrap();
        // F is a rigid rotation with first column (0.8,0.6,0). Transverse
        // eigenvectors cancel the nonzero pressure couplings, so |B*phi|
        // alone cannot scale their pressure-equation backward error.
        for mode in &spectrum.modes()[..2] {
            near(
                mode.displacement.lambda,
                42.0 * mu / rho,
                1e-8 * 42.0 * mu / rho,
            );
            let phi_scale = mode
                .displacement
                .phi
                .iter()
                .map(|v| v.abs())
                .fold(0.0, f64::max);
            assert!(
                mode.pressure_shape
                    .iter()
                    .all(|p| p.abs() <= 1e-8 * bulk * phi_scale)
            );
            assert!(
                mode.pressure_equation_relative_residual
                    <= modal_settings().pressure_residual_tolerance
            );
        }
        near(
            spectrum.modes()[2].displacement.lambda,
            42.0 * (bulk + 4.0 * mu / 3.0) / rho,
            1e-8 * 42.0 * bulk / rho,
        );
    });
}

#[test]
fn mixed_vibration_refuses_exact_constraints_storage_caps_and_cancellation() {
    use fs_solid::hyper3::mixed::vibration::MixedVibrationError;
    with_cx(|cx| {
        let p1 = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (nodes, elements) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != 3)
            .map(|d| (d, 0.0))
            .collect();
        let material = law(2.0, 200.0);
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &fixed,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        assert!(matches!(
            p.solve_vibration(&[1.0], modal_settings(), cx),
            Err(MixedVibrationError::Mixed(
                MixedUpError::InvalidInput { .. }
            ))
        ));
        p.volumetric = VolumetricTreatment::OriginalMaterial;
        for kind in 0..3 {
            let mut settings = modal_settings();
            match kind {
                0 => settings.spectrum.maximum_free_dofs = 0,
                1 => settings.maximum_condensation_bytes = 0,
                _ => settings.maximum_mass_assembly_bytes = 0,
            }
            assert!(matches!(
                p.solve_vibration(&[1.0], settings, cx),
                Err(MixedVibrationError::Mixed(
                    MixedUpError::InvalidInput { .. }
                ))
            ));
        }
        let gate = CancelGate::new();
        gate.request();
        with_gate(&gate, |cancelled| {
            assert!(matches!(
                p.solve_vibration(&[1.0], modal_settings(), cancelled),
                Err(MixedVibrationError::Mixed(MixedUpError::Cancelled))
            ));
        });
    });
}


// Independent unit-tetra polynomial model. Every P2 vector field that vanishes
// on z=0 and has zero divergence is spanned by
// (z,0,0), (xz,0,-z²/2), (yz,0,0), (z²,0,0),
// (0,z,0), (0,xz,0), (0,yz,-z²/2), (0,z²,0).
// The integrals below use exact simplex monomials, not native shapes, pressure
// matrices, quadrature, QR or the production projected stiffness/mass.
fn clamped_face_solenoidal_pencil(mu: f64, rho: f64) -> (Vec<f64>, Vec<f64>) {
    let mut fields = [[[0.0; 4]; 3]; 8]; // coefficients of z, xz, yz, z²
    for axis in 0..2 {
        for monomial in 0..4 {
            fields[4 * axis + monomial][axis][monomial] = 1.0;
        }
    }
    fields[1][2][3] = -0.5;
    fields[6][2][3] = -0.5;
    // Integral products of (1,x,y,z), then (z,xz,yz,z²), over the unit tetra.
    let linear = [
        [1.0 / 6.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0],
        [1.0 / 24.0, 1.0 / 60.0, 1.0 / 120.0, 1.0 / 120.0],
        [1.0 / 24.0, 1.0 / 120.0, 1.0 / 60.0, 1.0 / 120.0],
        [1.0 / 24.0, 1.0 / 120.0, 1.0 / 120.0, 1.0 / 60.0],
    ];
    let quadratic = [
        [1.0 / 60.0, 1.0 / 360.0, 1.0 / 360.0, 1.0 / 120.0],
        [1.0 / 360.0, 1.0 / 1260.0, 1.0 / 2520.0, 1.0 / 840.0],
        [1.0 / 360.0, 1.0 / 2520.0, 1.0 / 1260.0, 1.0 / 840.0],
        [1.0 / 120.0, 1.0 / 840.0, 1.0 / 840.0, 1.0 / 210.0],
    ];
    let mut gradients = [[[[0.0; 4]; 3]; 3]; 8];
    for i in 0..8 {
        for a in 0..3 {
            let p = fields[i][a];
            gradients[i][a][0] = [0.0, 0.0, 0.0, p[1]];
            gradients[i][a][1] = [0.0, 0.0, 0.0, p[2]];
            gradients[i][a][2] = [p[0], p[1], p[2], 2.0 * p[3]];
        }
    }
    let (mut stiffness, mut mass) = (vec![0.0; 64], vec![0.0; 64]);
    for i in 0..8 {
        for j in 0..8 {
            for a in 0..3 {
                for p in 0..4 {
                    for q in 0..4 {
                        mass[8 * i + j] += rho * fields[i][a][p]
                            * fields[j][a][q] * quadratic[p][q];
                        for b in 0..3 {
                            let ei = 0.5 * (gradients[i][a][b][p] + gradients[i][b][a][p]);
                            let ej = 0.5 * (gradients[j][a][b][q] + gradients[j][b][a][q]);
                            stiffness[8 * i + j] += 2.0 * mu * ei * ej * linear[p][q];
                        }
                    }
                }
            }
        }
    }
    (stiffness, mass)
}

#[test]
fn exact_incompressible_modes_match_independent_polynomials_and_objectivity() {
    use fs_material::hyper::cards::ResolvedHyperelastic;
    with_cx(|cx| {
        let p1 = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (nodes, elements) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let (mu, rho) = (2.0, 7.0);
        let (k, m) = clamped_face_solenoidal_pencil(mu, rho);
        let expected = fs_modal::eigh_gen_dense(&k, &m, 8).unwrap();
        for (cosine, sine, bulk) in [(1.0, 0.0, 20.0), (0.8, 0.6, 2e8)] {
            let initial: Vec<_> = nodes
                .iter()
                .map(|x| [
                    (cosine - 1.0) * x[0] - sine * x[1],
                    sine * x[0] + (cosine - 1.0) * x[1],
                    0.0,
                ])
                .collect();
            let fixed: Vec<_> = nodes
                .iter()
                .enumerate()
                .filter(|(_, x)| x[2] == 0.0)
                .flat_map(|(node, _)| {
                    let initial = &initial;
                    (0..3).map(move |a| (3 * node + a, initial[node][a]))
                })
                .collect();
            let card = modal_source_card(mu, bulk);
            let material = ResolvedHyperelastic::try_from_model_card(&card).unwrap();
            let mut p = problem(
                &nodes, &elements, &material, &fixed, &[],
                VolumetricTreatment::ExactIncompressible,
            );
            p.settings.load_steps = 1;
            let spectrum = if sine == 0.0 {
                p.solve_incompressible_vibration(&[rho], modal_settings(), cx)
            } else {
                p.solve_incompressible_vibration_from(
                    &initial, &[0.0; 4], &[rho], modal_settings(), cx,
                )
            }.unwrap();
            assert_eq!(spectrum.free_dofs().len(), 12);
            assert_eq!(spectrum.dimension(), 8);
            assert_eq!(spectrum.independent_pressure_indices(), &[0, 1, 2, 3]);
            assert_eq!(spectrum.basis().len(), 12 * 8);
            assert_eq!(spectrum.modes().len(), 8);
            assert_eq!(spectrum.mass_kg().len(), 12 * 12, "no pressure inertia");
            assert_eq!(spectrum.reference_density_kg_m3(), &[rho]);
            let MixedLawField::Uniform(retained) = spectrum.problem().material else {
                panic!("original source material was replaced")
            };
            assert!(core::ptr::eq(retained, &material));
            assert_eq!(retained.model_card(), &card);
            for (mode, oracle) in spectrum.modes().iter().zip(&expected) {
                near(mode.coordinates.lambda, oracle.lambda, 2e-8 * oracle.lambda);
                assert!(mode.coordinates.interval.0 > 0.0);
                assert!(mode.constraint_relative_residual
                    <= modal_settings().pressure_residual_tolerance);
                assert!(mode.momentum_relative_residual
                    <= modal_settings().spectrum.relative_residual_tolerance);
                let mut physical_norm = 0.0;
                for i in 0..12 {
                    for j in 0..12 {
                        physical_norm += mode.displacement_shape[i]
                            * spectrum.mass_kg()[12 * i + j] * mode.displacement_shape[j];
                    }
                }
                near(physical_norm, 1.0, 2e-11);
                assert_eq!(mode.pressure_shape.len(), 4);
            }
            // Neither the unused declared finite bulk modulus nor a rigid
            // superposed rotation changes this exact-constraint shear spectrum.
            // This is a one-element polynomial identity, not a locking study.
        }
    });
}

#[test]
fn exact_incompressible_modes_preserve_the_original_pressure_gauge() {
    with_cx(|cx| {
        let p1 = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (nodes, elements) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|d| ![0, 3, 7, 11].contains(d))
            .map(|d| (d, 0.0))
            .collect();
        let (mu, rho) = (2.0, 7.0);
        let material = law(mu, 20.0);
        let mut p = problem(
            &nodes, &elements, &material, &fixed, &[],
            VolumetricTreatment::ExactIncompressible,
        );
        p.settings.load_steps = 1;
        // Every P2 corner shape has zero integrated gradient. Constant pressure
        // is therefore a true gauge for these four free DOFs. The remaining
        // three independent constraints leave phi proportional to [1,-1,-1,-1].
        let mut normalized_pressures = Vec::new();
        for pin in [0, 1] {
            p.pressure_gauge = Some((pin, 0.0));
            let spectrum = p
                .solve_incompressible_vibration(&[rho], modal_settings(), cx)
                .unwrap();
            assert_eq!(spectrum.problem().pressure_gauge, Some((pin, 0.0)));
            assert_eq!(spectrum.free_dofs(), &[0, 3, 7, 11]);
            assert_eq!(spectrum.dimension(), 1);
            assert_eq!(spectrum.equilibrium().kkt_inertia.positive, 4);
            assert_eq!(spectrum.equilibrium().kkt_inertia.negative, 3);
            assert!(!spectrum.independent_pressure_indices().contains(&pin));
            let mode = &spectrum.modes()[0];
            // Independent exact integrals for this solenoidal shape:
            // stiffness=26 mu V/5; mass=11 rho V/210.
            near(mode.coordinates.lambda, 1092.0 * mu / (11.0 * rho), 2e-8);
            let a = mode.displacement_shape[0];
            for &value in &mode.displacement_shape[1..] {
                near(value / a, -1.0, 2e-11);
            }
            near(a * a * 11.0 * rho / 1260.0, 1.0, 2e-11);
            assert_eq!(mode.pressure_shape[pin], 0.0);
            let full = spectrum.mixed_tangent();
            let nu = 3 * nodes.len();
            for row in 0..4 {
                let mut residual = 0.0;
                let mut magnitude = 0.0;
                for (&dof, &value) in spectrum.free_dofs().iter().zip(&mode.displacement_shape) {
                    let term = full.get(nu + row, dof) * value;
                    residual += term;
                    magnitude += term.abs();
                }
                assert!(residual.abs() <= 1e-10 * magnitude);
                for col in 0..4 {
                    assert_eq!(full.get(nu + row, nu + col), 0.0);
                }
            }
            normalized_pressures.push(mode.pressure_shape.iter().map(|p| p / a).collect::<Vec<_>>());
        }
        let shift = normalized_pressures[1][0] - normalized_pressures[0][0];
        for i in 1..4 {
            near(normalized_pressures[1][i] - normalized_pressures[0][i], shift, 2e-8);
        }
    });
}

#[test]
fn exact_incompressible_vibration_retains_domain_capacity_and_gauge_refusals() {
    use fs_solid::hyper3::mixed::vibration::MixedVibrationError;
    with_cx(|cx| {
        let p1 = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (nodes, elements) = promote_p1_to_p2(&p1, &[[0, 1, 2, 3]], None, cx).unwrap();
        let fixed: Vec<_> = nodes.iter().enumerate()
            .filter(|(_, x)| x[2] == 0.0)
            .flat_map(|(node, _)| (0..3).map(move |a| (3 * node + a, 0.0)))
            .collect();
        let material = law(2.0, 20.0);
        let mut p = problem(
            &nodes, &elements, &material, &fixed, &[],
            VolumetricTreatment::OriginalMaterial,
        );
        assert!(matches!(
            p.solve_incompressible_vibration(&[1.0], modal_settings(), cx),
            Err(MixedVibrationError::Mixed(MixedUpError::InvalidInput { .. }))
        ));
        p.volumetric = VolumetricTreatment::ExactIncompressible;
        for kind in 0..3 {
            let mut settings = modal_settings();
            match kind {
                0 => settings.spectrum.maximum_free_dofs = 8, // 12 before reduction
                1 => settings.maximum_condensation_bytes = 0,
                _ => settings.maximum_mass_assembly_bytes = 0,
            }
            assert!(matches!(
                p.solve_incompressible_vibration(&[1.0], settings, cx),
                Err(MixedVibrationError::Mixed(MixedUpError::InvalidInput { .. }))
            ));
        }
        // With free face traction the mean pressure is physical, not a gauge.
        p.pressure_gauge = Some((0, 0.0));
        assert!(matches!(
            p.solve_incompressible_vibration(&[1.0], modal_settings(), cx),
            Err(MixedVibrationError::Mixed(MixedUpError::IncompatiblePressureGauge { node: 0 }))
        ));
        p.pressure_gauge = None;
        let gate = CancelGate::new();
        gate.request();
        with_gate(&gate, |cancelled| {
            assert!(matches!(
                p.solve_incompressible_vibration(&[1.0], modal_settings(), cancelled),
                Err(MixedVibrationError::Mixed(MixedUpError::Cancelled))
            ));
        });
    });
}
