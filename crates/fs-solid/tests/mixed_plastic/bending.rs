//! Small non-affine plastic bending diagnostic. The two axial mesh rungs and
//! bulk comparison do not establish an inf-sup theorem or a continuum bound.

use super::*;
use fs_solid::linear3::tet10::TET10_EDGES;
use std::collections::BTreeMap;

const LENGTH: f64 = 2.0;
const HEIGHT: f64 = 1.0;
const WIDTH: f64 = 0.5;

fn cantilever(n: usize, cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    let id = |i, j, k| (k * 2 + j) * (n + 1) + i;
    let mut points = Vec::new();
    for k in 0..=1 {
        for j in 0..=1 {
            for i in 0..=n {
                points.push([
                    LENGTH * i as f64 / n as f64,
                    HEIGHT * j as f64,
                    WIDTH * k as f64,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for i in 0..n {
        for permutation in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            let mut cell = [i, 0, 0];
            let mut tet = [id(i, 0, 0); 4];
            for (step, axis) in permutation.into_iter().enumerate() {
                cell[axis] += 1;
                tet[step + 1] = id(cell[0], cell[1], cell[2]);
            }
            tets.push(tet);
        }
    }
    promote_p1_to_p2(&points, &tets, None, cx).unwrap()
}

// Every selected face lies on the outer x=L plane. For a straight quadratic
// triangle the constant-traction corner weights are zero and each midside
// receives A/3; lumping this load onto corners would change the bending test.
fn end_traction(nodes: &[[f64; 3]], elements: &[[usize; 10]], traction: f64) -> Vec<(usize, f64)> {
    let mut loads = BTreeMap::<usize, f64>::new();
    for tet in elements {
        for omit in 0..4 {
            let corners: Vec<_> = (0..4).filter(|&a| a != omit).map(|a| tet[a]).collect();
            if !corners.iter().all(|&node| nodes[node][0] == LENGTH) {
                continue;
            }
            let [a, b, c] = [nodes[corners[0]], nodes[corners[1]], nodes[corners[2]]];
            let area = 0.5 * ((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])).abs();
            for (edge, endpoints) in TET10_EDGES.iter().enumerate() {
                if !endpoints.contains(&omit) {
                    *loads.entry(3 * tet[4 + edge] + 1).or_default() += traction * area / 3.0;
                }
            }
        }
    }
    let loads: Vec<_> = loads.into_iter().collect();
    near(
        loads.iter().map(|(_, force)| force).sum(),
        traction * HEIGHT * WIDTH,
        2e-15,
    );
    loads
}

fn response(n: usize, bulk: f64, treatment: VolumetricTreatment, cx: &Cx<'_>) -> f64 {
    let (nodes, elements) = cantilever(n, cx);
    assert_eq!(elements.len(), 6 * n);
    let fixed: Vec<_> = nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0)
        .flat_map(|(node, _)| (0..3).map(move |axis| (3 * node + axis, 0.0)))
        .collect();
    let material = law(bulk);
    let mut solid = MixedPlasticTetSolid::try_new(
        &nodes,
        &elements,
        MixedPlasticLawField::Uniform(&material),
        treatment,
        TetAssemblyBudget::standard(),
        4,
        elements.len() * 64,
        cx,
    )
    .unwrap();
    let mut solver = settings();
    solver.length_scale_m = LENGTH;
    let mut loaded_mean = 0.0;
    let mut accumulated_heat = 0.0;
    // These are separate physical endpoints; each Newton solve sees its own
    // fixed previous history. The second load exceeds the bending yield scale.
    for traction in [0.15, 0.3] {
        let forces = end_traction(&nodes, &elements, traction);
        let total_force = traction * HEIGHT * WIDTH;
        let step = solid.advance(&fixed, &forces, [0.0; 3], None,
            solver, controls(), cx).unwrap_or_else(|error| {
                panic!("plastic bending n={n}, bulk={bulk:e}, {treatment:?}, traction={traction}: {error}")
            });
        assert!(step.equilibrium.min_det_f > 0.9);
        assert!(step.solid.pressure_pa().iter().all(|p| p.is_finite()));
        assert!(
            step.equilibrium
                .reactions_n
                .iter()
                .all(|(_, f)| f.is_finite())
        );
        for axis in 0..3 {
            let reaction: f64 = step
                .equilibrium
                .reactions_n
                .iter()
                .filter(|(d, _)| d % 3 == axis)
                .map(|(_, value)| value)
                .sum();
            near(reaction, if axis == 1 { -total_force } else { 0.0 }, 2e-7);
        }
        loaded_mean = forces
            .iter()
            .map(|(d, force)| force * step.solid.displacement_m()[d / 3][d % 3])
            .sum::<f64>()
            / total_force;
        assert!(
            loaded_mean > 0.0 && loaded_mean < 0.3,
            "bounded upward bending response required, got {loaded_mean:e}"
        );
        near(
            step.energy.equilibrium_work_residual_j,
            0.0,
            controls().equilibrium_work_tolerance_j,
        );
        near(
            step.energy.pressure_constraint_work_j,
            0.0,
            controls().equilibrium_work_tolerance_j,
        );
        assert!(
            step.solid
                .quadrature_states()
                .iter()
                .zip(solid.quadrature_states())
                .all(|(new, old)| new.alpha >= old.alpha)
        );
        accumulated_heat += step.energy.plastic_dissipation_j;
        solid = step.solid;
    }
    let minimum_alpha = solid
        .quadrature_states()
        .iter()
        .map(|s| s.alpha)
        .fold(f64::INFINITY, f64::min);
    let maximum_alpha = solid
        .quadrature_states()
        .iter()
        .map(|s| s.alpha)
        .fold(0.0, f64::max);
    assert!(
        maximum_alpha > 1e-3 && maximum_alpha - minimum_alpha > 1e-3,
        "nonuniform yielded history required, min={minimum_alpha:e}, max={maximum_alpha:e}"
    );
    assert!(accumulated_heat > 0.0 && accumulated_heat.is_finite());
    assert!(solid.pressure_pa().iter().any(|p| p.abs() > 1e-3));
    if treatment == VolumetricTreatment::ExactIncompressible {
        assert_eq!(solid.bulk_energy_j(), 0.0);
    }
    println!(
        "mixed plastic bending elements={} bulk={bulk:e} treatment={treatment:?} mean_tip_m={loaded_mean:e} alpha_min={minimum_alpha:e} alpha_max={maximum_alpha:e} plastic_heat_j={accumulated_heat:e}",
        elements.len()
    );
    loaded_mean
}

#[test]
fn nonaffine_plastic_bending_retains_deflection_at_the_incompressible_limit() {
    with_cx(|cx| {
        let mut tips = Vec::new();
        for n in [1, 2] {
            let finite = response(n, 8e8, VolumetricTreatment::OriginalMaterial, cx);
            let exact = response(n, 1e300, VolumetricTreatment::ExactIncompressible, cx);
            let relative = (finite - exact).abs() / exact;
            assert!(
                relative < 3e-5,
                "high-bulk mixed response must approach exact incompressibility: n={n}, relative={relative:e}"
            );
            tips.push(exact);
        }
        let relative_change = (tips[1] - tips[0]).abs() / tips[1];
        assert!(
            relative_change < 0.35,
            "bounded two-rung bending response changed excessively: {tips:?}"
        );
        // Axial refinement only: this is a small regression against collapse
        // of the plastic bending response, not a converged 3-D Cook benchmark.
    });
}
