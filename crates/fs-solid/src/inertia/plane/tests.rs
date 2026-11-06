use super::*;
use crate::inertia::BackwardEulerInertia;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 72,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn close(a: f64, b: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= 4e-13 * a.abs().max(b.abs()).max(1.0),
        "{a:.17e} versus {b:.17e}"
    );
}

fn triangle() -> Mesh2 {
    Mesh2 {
        nodes: vec![[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]],
        elems: vec![vec![0, 1, 2]],
        patches: vec![],
    }
}

fn warped_quad() -> Mesh2 {
    Mesh2 {
        nodes: vec![[0.0, 0.0], [2.0, 0.0], [2.5, 1.75], [0.0, 1.5]],
        elems: vec![vec![0, 1, 2, 3]],
        patches: vec![],
    }
}

#[test]
fn p1_mass_matches_exact_degree_two_coefficients_and_retains_nonrigid_inertia() {
    with_cx(false, |cx| {
        let mesh = triangle();
        let mass = ReferenceMass::from_plane(&mesh, 0.2, &[7.0], 1 << 20, cx).unwrap();
        let inventory = 3.0 * 0.2 * 7.0;
        close(mass.reference_mass_kg(), inventory);
        assert_eq!(mass.ndofs(), 6);
        for a in 0..3 {
            for b in 0..3 {
                for i in 0..2 {
                    for j in 0..2 {
                        let expected = if i == j {
                            inventory / 12.0 * if a == b { 2.0 } else { 1.0 }
                        } else {
                            0.0
                        };
                        close(mass.matrix().get(2 * a + i, 2 * b + j), expected);
                    }
                }
                assert_eq!(
                    mass.matrix().get(2 * a, 2 * b).to_bits(),
                    mass.matrix().get(2 * b, 2 * a).to_bits()
                );
            }
        }
        // This velocity has zero centroid value. Centroid-only material
        // quadrature would wrongly assign it zero kinetic energy.
        let v = [1.0, 0.0, -1.0, 0.0, 0.0, 0.0];
        let p = mass.action(&v, cx).unwrap();
        let kinetic: f64 = 0.5 * v.iter().zip(p).map(|(x, y)| x * y).sum::<f64>();
        close(kinetic, inventory / 12.0);
        assert!(kinetic > 0.0);
        // Independent simplex moment: mean x²=2²/6 and mean y²=3²/6.
        let v = mesh.nodes.as_flattened();
        let p = mass.action(v, cx).unwrap();
        close(
            0.5 * v.iter().zip(p).map(|(x, y)| x * y).sum::<f64>(),
            0.5 * inventory * (4.0 + 9.0) / 6.0,
        );
    });
}

#[test]
fn warped_q1_mass_matches_independent_polynomial_moments_and_original_density_inventory() {
    with_cx(false, |cx| {
        let mesh = warped_quad();
        let rho = 3.7;
        let thickness = 0.11;
        let mass = ReferenceMass::from_plane(&mesh, thickness, &[rho], 1 << 20, cx).unwrap();
        // Independent map on [0,1]²:
        // x=2s+0.5st, y=1.5t+0.25st, detJ=3+0.5s+0.75t.
        // Integrate shape pairs by elementary polynomial antiderivatives.
        let ends = [(0, 0), (1, 0), (1, 1), (0, 1)];
        let moment = |a: usize, b: usize, degree: usize| match (a == b, a, degree) {
            (true, _, 0) => 1.0 / 3.0,
            (false, _, 0) => 1.0 / 6.0,
            (true, 1, 1) => 1.0 / 4.0,
            (_, _, 1) => 1.0 / 12.0,
            _ => unreachable!(),
        };
        close(
            mass.reference_mass_kg(),
            rho * thickness * (3.0 + 0.25 + 0.375),
        );
        for a in 0..4 {
            for b in 0..4 {
                let s0 = moment(ends[a].0, ends[b].0, 0);
                let s1 = moment(ends[a].0, ends[b].0, 1);
                let t0 = moment(ends[a].1, ends[b].1, 0);
                let t1 = moment(ends[a].1, ends[b].1, 1);
                let expected = rho * thickness * (3.0 * s0 * t0 + 0.5 * s1 * t0 + 0.75 * s0 * t1);
                for i in 0..2 {
                    close(mass.matrix().get(2 * a + i, 2 * b + i), expected);
                    close(mass.matrix().get(2 * a + i, 2 * b + (1 - i)), 0.0);
                }
            }
        }
        let mut shifted = mesh.clone();
        for point in &mut shifted.nodes {
            point[0] += 4294967296.0;
            point[1] -= 2147483648.0;
        }
        let shifted = ReferenceMass::from_plane(&shifted, thickness, &[rho], 1 << 20, cx).unwrap();
        assert_eq!(mass.matrix(), shifted.matrix());
        assert_eq!(
            mass.reference_mass_kg().to_bits(),
            shifted.reference_mass_kg().to_bits()
        );
        let mut rotated = mesh.clone();
        for p in &mut rotated.nodes {
            *p = [-p[1], p[0]];
        }
        let rotated = ReferenceMass::from_plane(&rotated, thickness, &[rho], 1 << 20, cx).unwrap();
        let thicker =
            ReferenceMass::from_plane(&mesh, 5.0 * thickness, &[rho], 1 << 20, cx).unwrap();
        for i in 0..8 {
            for j in 0..8 {
                close(rotated.matrix().get(i, j), mass.matrix().get(i, j));
                close(thicker.matrix().get(i, j), 5.0 * mass.matrix().get(i, j));
            }
        }
        // Independent densities remain attached to their original cells.
        let strip = Mesh2::quads(3.0, 1.0, 2, 1);
        let heterogeneous =
            ReferenceMass::from_plane(&strip, 0.1, &[2.0, 5.0], 1 << 20, cx).unwrap();
        close(heterogeneous.reference_mass_kg(), 1.5 * 0.1 * (2.0 + 5.0));
        let mut translation = vec![0.0; heterogeneous.ndofs()];
        for node in 0..strip.nodes.len() {
            translation[2 * node] = 1.0;
        }
        let momentum = heterogeneous.action(&translation, cx).unwrap();
        close(
            momentum.iter().step_by(2).sum(),
            heterogeneous.reference_mass_kg(),
        );
        for (node, p) in strip.nodes.iter().enumerate() {
            let expected = if p[0] == 0.0 {
                1.5 * 0.1 * 2.0 / 4.0
            } else if p[0] == 3.0 {
                1.5 * 0.1 * 5.0 / 4.0
            } else {
                1.5 * 0.1 * (2.0 + 5.0) / 4.0
            };
            close(momentum[2 * node], expected);
        }
    });
}

#[test]
fn plane_motion_uses_all_components_and_the_original_backward_euler_work_identity() {
    with_cx(false, |cx| {
        for mesh in [triangle(), warped_quad()] {
            let mass = ReferenceMass::from_plane(&mesh, 0.2, &[7.0], 1 << 20, cx).unwrap();
            let n = mass.ndofs();
            let dt = 0.125;
            let v0 = [0.25, -0.5];
            let acceleration = [0.75, 0.25];
            let next = [v0[0] + dt * acceleration[0], v0[1] + dt * acceleration[1]];
            let old_u = vec![0.0; n];
            let old_v: Vec<_> = (0..n).map(|i| v0[i % 2]).collect();
            let u: Vec<_> = (0..n).map(|i| dt * next[i % 2]).collect();
            let inertia = BackwardEulerInertia::new(&mass, &old_u, &old_v, dt, cx).unwrap();
            let (velocity, energy) = inertia.motion_flat(&u, cx).unwrap();
            assert_eq!(velocity.len(), n);
            for i in 0..n {
                close(velocity[i], next[i % 2]);
            }
            let m = mass.reference_mass_kg();
            let square = |x: [f64; 2]| x[0] * x[0] + x[1] * x[1];
            close(energy.kinetic_before_j, 0.5 * m * square(v0));
            close(energy.kinetic_after_j, 0.5 * m * square(next));
            close(
                energy.kinetic_change_j,
                0.5 * m * (square(next) - square(v0)),
            );
            close(
                energy.numerical_loss_j,
                0.5 * m * dt * dt * square(acceleration),
            );
            close(
                energy.inertial_work_j,
                m * dt * (next[0] * acceleration[0] + next[1] * acceleration[1]),
            );
            close(energy.balance_error_j, 0.0);
            // The legacy carrier refuses an incomplete triple rather than
            // silently truncating this four-node plane's last two components.
            if !n.is_multiple_of(3) {
                assert!(inertia.motion(&u, cx).is_err());
            }
            let mut tangent = Coo::new(n, n);
            inertia.add_tangent(&mut tangent, cx).unwrap();
            let tangent = tangent.assemble();
            let mut filtered = Coo::new(n, n);
            inertia
                .add_tangent_where(&mut filtered, |i| i >= 2, cx)
                .unwrap();
            let filtered = filtered.assemble();
            for i in 0..n {
                for j in 0..n {
                    close(tangent.get(i, j), mass.matrix().get(i, j) / (dt * dt));
                    close(
                        filtered.get(i, j),
                        if i >= 2 && j >= 2 {
                            tangent.get(i, j)
                        } else {
                            0.0
                        },
                    );
                }
            }
            let direction: Vec<_> = (0..n).map(|i| 0.01 * (i as f64 - 3.0)).collect();
            let h = 1e-6;
            let plus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u + h * d).collect();
            let minus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u - h * d).collect();
            let (mut fp, mut fm) = (vec![0.0; n], vec![0.0; n]);
            inertia.add_force(&plus, &mut fp, cx).unwrap();
            inertia.add_force(&minus, &mut fm, cx).unwrap();
            for i in 0..n {
                let expected: f64 = (0..n).map(|j| tangent.get(i, j) * direction[j]).sum();
                assert!(((fp[i] - fm[i]) / (2.0 * h) - expected).abs() < 2e-8);
            }
            let flight: Vec<_> = old_v.iter().map(|v| dt * v).collect();
            let (v, energy) = inertia.motion_flat(&flight, cx).unwrap();
            assert_eq!(v, old_v);
            assert_eq!(energy.numerical_loss_j, 0.0);
            close(energy.kinetic_before_j, energy.kinetic_after_j);
            close(energy.inertial_work_j, 0.0);
        }
    });
}

#[test]
fn plane_mass_refuses_coverage_budget_massless_nodes_and_hidden_corner_folds() {
    with_cx(false, |cx| {
        let mesh = triangle();
        for density in [&[][..], &[0.0][..], &[f64::NAN][..], &[7.0, 7.0][..]] {
            assert!(ReferenceMass::from_plane(&mesh, 0.2, density, 1 << 20, cx).is_err());
        }
        for thickness in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(ReferenceMass::from_plane(&mesh, thickness, &[7.0], 1 << 20, cx).is_err());
        }
        assert!(ReferenceMass::from_plane(&mesh, 0.2, &[7.0], 0, cx).is_err());
        let mut unused = mesh.clone();
        unused.nodes.push([2.0, 3.0]);
        assert!(ReferenceMass::from_plane(&unused, 0.2, &[7.0], 1 << 20, cx).is_err());
        for conn in [vec![0, 2, 1], vec![0, 1, 1], vec![0, 1, 3], vec![0, 1]] {
            let mut bad = mesh.clone();
            bad.elems[0] = conn;
            assert!(ReferenceMass::from_plane(&bad, 0.2, &[7.0], 1 << 20, cx).is_err());
        }
        let mut bad = mesh.clone();
        bad.nodes[0][0] = f64::NAN;
        assert!(ReferenceMass::from_plane(&bad, 0.2, &[7.0], 1 << 20, cx).is_err());
        let fold = Mesh2 {
            nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.4, 0.4], [0.0, 1.0]],
            elems: vec![vec![0, 1, 2, 3]],
            patches: vec![],
        };
        assert!(
            quad_points(4).into_iter().all(|(xi, eta, _)| shapes_at(
                &fold.nodes,
                &fold.elems[0],
                xi,
                eta
            )
            .2 > 0.0)
        );
        assert!(checked_reference_element(&fold.nodes, &fold.elems[0]).is_none());
        assert!(ReferenceMass::from_plane(&fold, 0.2, &[7.0], 1 << 20, cx).is_err());
    });
}

#[test]
fn plane_mass_and_flat_motion_observe_cancellation_before_publishing() {
    with_cx(false, |cx| {
        let mesh = warped_quad();
        let mass = ReferenceMass::from_plane(&mesh, 0.2, &[7.0], 1 << 20, cx).unwrap();
        let zero = [0.0; 8];
        let inertia = BackwardEulerInertia::new(&mass, &zero, &zero, 0.1, cx).unwrap();
        with_cx(true, |cancelled| {
            assert!(matches!(
                ReferenceMass::from_plane(&mesh, 0.2, &[7.0], 1 << 20, cancelled),
                Err(InertiaError::Cancelled)
            ));
            assert!(matches!(
                inertia.motion_flat(&zero, cancelled),
                Err(InertiaError::Cancelled)
            ));
            assert!(matches!(
                inertia.add_tangent_where(&mut Coo::new(8, 8), |_| true, cancelled),
                Err(InertiaError::Cancelled)
            ));
        });
    });
}
