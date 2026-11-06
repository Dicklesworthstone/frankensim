//! Continuous fluid density and exact pressure moments through native solid owners.
use super::*;

const BETA: f64 = 2.4;

// Independent divided differences for the volume mean of positive affine
// depth and squared depth. Fully wet values use their elementary simplex
// integrals, avoiding cancellation when the surface is far above the body.
fn positive_mean(q: [f64; 4], squared: bool) -> f64 {
    if q.iter().all(|&x| x >= 0.0) {
        if squared {
            let sum: f64 = q.iter().sum();
            return (sum * sum + q.iter().map(|x| x * x).sum::<f64>()) / 20.0;
        }
        return q.iter().sum::<f64>() / 4.0;
    }
    (0..4)
        .map(|i| {
            let denominator = (0..4)
                .filter(|&j| i != j)
                .map(|j| q[i] - q[j])
                .product::<f64>();
            if squared {
                q[i].max(0.0).powi(5) / (10.0 * denominator)
            } else {
                q[i].max(0.0).powi(4) / (4.0 * denominator)
            }
        })
        .sum()
}

fn density_potential(u: &[f64; 12], beta: f64) -> f64 {
    let x: [[f64; 3]; 4] = std::array::from_fn(|node| {
        std::array::from_fn(|axis| NODES[node][axis] + u[3 * node + axis])
    });
    let volume = dot(sub(x[1], x[0]), cross(sub(x[2], x[0]), sub(x[3], x[0]))) / 6.0;
    let magnitude = dot(GRAVITY, GRAVITY).sqrt();
    let depth = x.map(|point| dot(GRAVITY, sub(point, ORIGIN)) / magnitude);
    let mut mean = 0.5 * beta * positive_mean(depth, true);
    let mut previous_density = 0.0;
    for layer in LAYERS {
        mean += (layer.density_kg_m3 - previous_density)
            * positive_mean(depth.map(|d| d - layer.top_depth_m), false);
        previous_density = layer.density_kg_m3;
    }
    volume * magnitude * mean
}

#[test]
fn linear_density_changes_archimedes_force_torque_and_heave_stiffness() {
    with_cx(|_, cx| {
        let profile = StratifiedHydrostaticPressure::new(
            [0.0, 0.0, 2.0],
            [0.0, 0.0, -10.0],
            &[HydrostaticDensityLayer {
                top_depth_m: 0.0,
                density_kg_m3: 1000.0,
            }],
            1,
            cx,
        )
        .unwrap()
        .with_density_gradient(100.0)
        .unwrap();
        let surface = surface(cx);
        let force = surface.force_n(&[0.0; 12], profile.loading(), cx).unwrap();
        let tangent = surface
            .force_tangent(&[0.0; 12], profile.loading(), cx)
            .unwrap();
        let mut resultant = [0.0; 3];
        let mut torque = [0.0; 3];
        for node in 0..4 {
            let f = std::array::from_fn(|axis| force[3 * node + axis]);
            let moment = cross(NODES[node], f);
            for axis in 0..3 {
                resultant[axis] += f[axis];
                torque[axis] += moment[axis];
            }
        }
        // V=1/6, mean(z)=1/4, mean(x*z)=1/20:
        // mean(rho)=1175 and mean(rho*x)=mean(rho*y)=295.
        near(resultant[0], 0.0, 5e-11);
        near(resultant[1], 0.0, 5e-11);
        near(resultant[2], 10.0 * 1175.0 / 6.0, 5e-11);
        near(torque[0], 10.0 * 295.0 / 6.0, 5e-11);
        near(torque[1], -10.0 * 295.0 / 6.0, 5e-11);
        near(torque[2], 0.0, 5e-11);
        let mut heave = 0.0;
        for row in [2, 5, 8, 11] {
            for column in [2, 5, 8, 11] {
                heave += tangent.get(row, column);
            }
        }
        near(heave, -10.0 * 100.0 / 6.0, 5e-11);
        near(profile.pressure_pa([0.0; 3], cx).unwrap(), 22000.0, 1e-10);
        near(profile.density_at_depth_kg_m3(2.0).unwrap(), 1200.0, 1e-13);
    });
}

#[test]
fn deforming_partial_immersion_matches_independent_potential_and_tangent() {
    with_cx(|_, cx| {
        let profile = profile(cx).with_density_gradient(BETA).unwrap();
        let surface = surface(cx);
        let (before, after) = endpoints();
        let force = surface.force_n(&after, profile.loading(), cx).unwrap();
        let tangent = surface
            .force_tangent(&after, profile.loading(), cx)
            .unwrap();
        let epsilon = 2e-6;
        for column in 0..12 {
            let mut plus = after;
            let mut minus = after;
            plus[column] += epsilon;
            minus[column] -= epsilon;
            near(
                force[column],
                -(density_potential(&plus, BETA) - density_potential(&minus, BETA))
                    / (2.0 * epsilon),
                3e-8,
            );
            let forward = surface.force_n(&plus, profile.loading(), cx).unwrap();
            let backward = surface.force_n(&minus, profile.loading(), cx).unwrap();
            for row in 0..12 {
                near(
                    tangent.get(row, column),
                    (forward[row] - backward[row]) / (2.0 * epsilon),
                    3e-8,
                );
                near(tangent.get(row, column), tangent.get(column, row), 3e-12);
            }
        }
        let work = surface
            .work_j(&before, &after, profile.loading(), cx)
            .unwrap();
        near(
            work,
            density_potential(&before, BETA) - density_potential(&after, BETA),
            3e-12,
        );
        near(
            surface
                .work_j(&after, &before, profile.loading(), cx)
                .unwrap(),
            -work,
            3e-12,
        );
        let middle: [f64; 12] = std::array::from_fn(|i| 0.4 * before[i] + 0.6 * after[i]);
        near(
            surface
                .work_j(&before, &middle, profile.loading(), cx)
                .unwrap()
                + surface
                    .work_j(&middle, &after, profile.loading(), cx)
                    .unwrap(),
            work,
            3e-12,
        );
        for factor in [0.0, 0.25, 2.0] {
            let loading = TetPressure::Stratified {
                profile: &profile,
                load_factor: factor,
            };
            let scaled = surface.force_n(&after, loading, cx).unwrap();
            for i in 0..12 {
                near(scaled[i], factor * force[i], 3e-12);
            }
            near(
                surface.work_j(&before, &after, loading, cx).unwrap(),
                factor * work,
                3e-12,
            );
        }
    });
}

#[test]
fn native_discrete_gradient_closes_density_gradient_work_and_endpoint_tangent() {
    with_cx(|_, cx| {
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 3.0,
            },
            5.0,
        )
        .unwrap();
        let problem = HyperTetProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            material: HyperelasticLawField::Uniform(&material),
            prescribed_m: &[],
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings::default(),
        };
        let prepared = problem.prepare(cx).unwrap();
        let surface = surface(cx);
        let profile = profile(cx).with_density_gradient(BETA).unwrap();
        for submerged_shift in [0.0, -4.0] {
            let (mut before, mut after) = endpoints();
            for node in 0..4 {
                before[3 * node + 2] += submerged_shift;
                after[3 * node + 2] += submerged_shift;
            }
            let force = prepared
                .energy_momentum_force_with_pressure(
                    &before,
                    &after,
                    &surface,
                    profile.loading(),
                    cx,
                )
                .unwrap();
            let reverse = prepared
                .energy_momentum_force_with_pressure(
                    &after,
                    &before,
                    &surface,
                    profile.loading(),
                    cx,
                )
                .unwrap();
            let work = (0..12)
                .map(|i| force[i] * (after[i] - before[i]))
                .sum::<f64>();
            let stored = problem.stored_energy_j(&after, cx).unwrap()
                - problem.stored_energy_j(&before, cx).unwrap();
            near(
                work,
                stored + density_potential(&after, BETA) - density_potential(&before, BETA),
                2e-11,
            );
            for i in 0..12 {
                near(force[i], reverse[i], 2e-11);
            }
            let tangent = prepared
                .energy_momentum_tangent_with_pressure(
                    &before,
                    &after,
                    &surface,
                    profile.loading(),
                    cx,
                )
                .unwrap();
            let epsilon = 2e-6;
            for column in 0..12 {
                let mut plus = after;
                let mut minus = after;
                plus[column] += epsilon;
                minus[column] -= epsilon;
                let forward = prepared
                    .energy_momentum_force_with_pressure(
                        &before,
                        &plus,
                        &surface,
                        profile.loading(),
                        cx,
                    )
                    .unwrap();
                let backward = prepared
                    .energy_momentum_force_with_pressure(
                        &before,
                        &minus,
                        &surface,
                        profile.loading(),
                        cx,
                    )
                    .unwrap();
                for row in 0..12 {
                    near(
                        tangent.get(row, column),
                        (forward[row] - backward[row]) / (2.0 * epsilon),
                        2e-6,
                    );
                }
            }
        }
    });
}

#[test]
fn coplanar_emergence_keeps_the_exact_density_gradient_cap_work() {
    with_cx(|_, cx| {
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 3.0,
            },
            5.0,
        )
        .unwrap();
        let problem = HyperTetProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            material: HyperelasticLawField::Uniform(&material),
            prescribed_m: &[],
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings::default(),
        };
        let prepared = problem.prepare(cx).unwrap();
        let surface = surface(cx);
        let profile = StratifiedHydrostaticPressure::new(
            [0.0; 3],
            [0.0, 0.0, -1.0],
            &[HydrostaticDensityLayer {
                top_depth_m: 0.0,
                density_kg_m3: 1.3,
            }],
            1,
            cx,
        )
        .unwrap()
        .with_density_gradient(BETA)
        .unwrap();
        let a: f64 = 0.1;
        let b = 1.0 - a;
        let before: [f64; 12] = std::array::from_fn(|d| if d % 3 == 2 { -a } else { 0.0 });
        let after = [0.0; 12];
        // Integrate p(d)*(b+d)^2/2 over d in [0,a]; the horizontal
        // tetrahedron slices provide this area exactly.
        let cap_potential = 0.5
            * (1.3 * (b * b * a * a / 2.0 + 2.0 * b * a.powi(3) / 3.0 + a.powi(4) / 4.0)
                + 0.5
                    * BETA
                    * (b * b * a.powi(3) / 3.0 + 2.0 * b * a.powi(4) / 4.0 + a.powi(5) / 5.0));
        let force = prepared
            .energy_momentum_force_with_pressure(&before, &after, &surface, profile.loading(), cx)
            .unwrap();
        let tangent = prepared
            .energy_momentum_tangent_with_pressure(&before, &after, &surface, profile.loading(), cx)
            .unwrap();
        let work = (0..12)
            .map(|d| force[d] * (after[d] - before[d]))
            .sum::<f64>();
        near(work, -cap_potential, 2e-13);
        near(
            surface
                .work_j(&before, &after, profile.loading(), cx)
                .unwrap(),
            cap_potential,
            2e-13,
        );
        assert_eq!(
            surface.force_n(&after, profile.loading(), cx).unwrap(),
            vec![0.0; 12]
        );
        for row in 0..12 {
            assert!(force[row].is_finite());
            for column in 0..12 {
                assert!(tangent.get(row, column).is_finite());
            }
        }
    });
}

#[test]
fn density_gradient_zero_limits_physical_overflow_and_cancellation_are_explicit() {
    with_cx(|gate, cx| {
        let base = profile(cx);
        for beta in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(base.clone().with_density_gradient(beta).is_err());
        }
        let surface = surface(cx);
        let (before, after) = endpoints();
        let baseline_force = surface.force_n(&after, base.loading(), cx).unwrap();
        let baseline_work = surface.work_j(&before, &after, base.loading(), cx).unwrap();
        for identical in [
            base.clone().with_density_gradient(0.0).unwrap(),
            base.clone()
                .with_density_gradient(BETA)
                .unwrap()
                .with_density_gradient(0.0)
                .unwrap(),
        ] {
            assert_eq!(
                surface.force_n(&after, identical.loading(), cx).unwrap(),
                baseline_force
            );
            assert_eq!(
                surface
                    .work_j(&before, &after, identical.loading(), cx)
                    .unwrap(),
                baseline_work
            );
        }
        let field = base.with_density_gradient(BETA).unwrap();
        assert_eq!(field.density_at_depth_kg_m3(-0.1).unwrap(), 0.0);
        near(
            field.density_at_depth_kg_m3(0.2).unwrap(),
            2.1 + BETA * 0.2,
            1e-15,
        );
        near(
            field.density_at_depth_kg_m3(0.55).unwrap(),
            0.9 + BETA * 0.55,
            1e-15,
        );
        let zero = StratifiedHydrostaticPressure::new(ORIGIN, [0.0; 3], &LAYERS, 3, cx)
            .unwrap()
            .with_density_gradient(BETA)
            .unwrap();
        assert_eq!(
            surface.force_n(&after, zero.loading(), cx).unwrap(),
            vec![0.0; 12]
        );
        assert_eq!(
            surface.work_j(&before, &after, zero.loading(), cx).unwrap(),
            0.0
        );
        // Small gravity does not excuse an unrepresentable physical density.
        let overflowing = StratifiedHydrostaticPressure::new(
            [0.0, 0.0, 2.0],
            [0.0, 0.0, -1e-100],
            &LAYERS[..1],
            1,
            cx,
        )
        .unwrap()
        .with_density_gradient(1e308)
        .unwrap();
        assert!(overflowing.density_at_depth_kg_m3(2.0).is_err());
        assert!(overflowing.pressure_pa([0.0; 3], cx).is_err());
        assert!(
            surface
                .force_n(&[0.0; 12], overflowing.loading(), cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            surface.force_n(&after, field.loading(), cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            surface.work_j(&before, &after, field.loading(), cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}
