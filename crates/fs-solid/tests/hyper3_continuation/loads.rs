//! Native follower loads: independent pressure equilibria, current derivatives,
//! original gravity/static owners, physical reactions and transactional refusal.
use super::*;
use fs_solid::hyper3::continuation::HyperTetPathLoads;
use fs_solid::hyper3::gravity::{GravityField3, GravitySample3, ReferenceGravityDensity};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn affine_gravity(x: [f64; 3]) -> Result<GravitySample3, HyperTetError> {
    let b = [0.03, -0.04, 0.08];
    let h = [[-0.3, 0.08, -0.04], [0.08, -0.2, 0.06], [-0.04, 0.06, -0.4]];
    let hx: [f64; 3] = core::array::from_fn(|i| (0..3).map(|j| h[i][j] * x[j]).sum::<f64>());
    Ok(GravitySample3 {
        potential_j_kg: -(0..3).map(|i| x[i] * (b[i] + 0.5 * hx[i])).sum::<f64>(),
        acceleration_m_s2: core::array::from_fn(|i| b[i] + hx[i]),
        acceleration_gradient_s2: h,
    })
}

#[test]
fn pure_pressure_crosses_the_native_hydrostatic_fold_and_resumes_identically() {
    with_cx(|_, cx| {
        let x = nodes(1.0);
        let (mu, lame, pressure) = (2.0, 5.0, 8.0);
        let law =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda: lame }, 6.0).unwrap();
        let p = problem(&x, &law, &DIAGONAL_SUPPORTS, &[]);
        let surface = ClosedTetPressureSurface::new(&x, &TETS, p.budget, 4, cx).unwrap();
        let path = p
            .prepare_continuation_with_loads(
                HyperTetPathScales {
                    displacement_m: 1.0,
                    force_n: pressure / 6.0,
                },
                HyperTetPathLoads {
                    pressure: Some((&surface, TetPressure::Uniform(-pressure))),
                    gravity: None,
                },
                cx,
            )
            .unwrap();
        assert_eq!(path.free_dofs(), &[3, 7, 11]);
        assert!(path.reference_force_n().iter().all(|&f| f == 0.0));
        assert!(core::ptr::eq(
            path.held_loads().pressure.unwrap().0,
            &surface
        ));
        let controls = ArcSettings {
            ds: 0.06,
            ds_max: 0.08,
            ..settings()
        };
        let initial = path.start(&controls, cx).unwrap();
        let result = path.advance(&initial, &controls, 36, cx).unwrap();
        let middle = path.advance(&initial, &controls, 13, cx).unwrap();
        same_state(&result, &path.advance(&middle, &controls, 23, cx).unwrap());
        // F=sI: P=mu(s-1/s)+3*lame*ln(s)/s; each free pressure force
        // is pressure*s^2/6. This is the original NH law's analytic branch,
        // not a calibrated inflation limit or a stability certificate.
        let load = |s: f64| (mu * (s - 1.0 / s) + 3.0 * lame * s.ln() / s) / (pressure * s * s);
        let (mut lo, mut hi) = (1.0, 3.0);
        for _ in 0..60 {
            let s: f64 = 0.5 * (lo + hi);
            if 3.0 * mu + 3.0 * lame - mu * s * s - 9.0 * lame * s.ln() > 0.0 {
                lo = s;
            } else {
                hi = s;
            }
        }
        let fold = 0.5 * (lo + hi);
        let s = 1.0 + result.u[0];
        assert!(s > fold + 0.2 && result.tangent.1 < 0.0);
        assert!(result.lambda < load(fold));
        assert!(
            result
                .events
                .iter()
                .any(|e| matches!(e, PathEvent::LimitPoint { .. }))
        );
        for &(lambda, q) in &result.trace {
            near(lambda, load(1.0 + q), 2e-9);
        }
        for &q in &result.u {
            near(1.0 + q, s, 2e-9);
        }
        for derivative in path.load_vector_at(&result.u, result.lambda, cx).unwrap() {
            near(derivative, s * s, 2e-9);
        }
        let response = path.response(&result, cx).unwrap();
        let energy =
            (1.5 * mu * (s * s - 1.0) - 3.0 * mu * s.ln() + 4.5 * lame * s.ln().powi(2)) / 6.0;
        near(response.strain_energy_j, energy, 2e-9);
        assert!(
            response
                .free_force_residual_n
                .iter()
                .all(|r| r.abs() < 1e-9)
        );
    });
}

#[test]
fn pressure_variants_and_graded_gravity_have_current_q_and_parameter_derivatives() {
    with_cx(|_, cx| {
        let x = nodes(1.0);
        let law = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 8.0,
                lambda: 14.0,
            },
            5.0,
        )
        .unwrap();
        // Six free components retain shear and the nonsymmetric patch tangent.
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let forces = [(3, 0.07), (0, 0.11)];
        let mut p = problem(&x, &law, &fixed, &forces);
        p.body_force_n_m3 = [0.03, -0.02, 0.05];
        let surface = ClosedTetPressureSurface::new(&x, &TETS, p.budget, 4, cx).unwrap();
        let face_values = [0.1, 0.5, -0.2, 0.0];
        let density = [2.0, 3.0, 5.0, 7.0];
        let fields = [
            TetPressure::Uniform(0.22),
            TetPressure::PerFace(&face_values),
            TetPressure::Affine {
                origin_m: [0.1, -0.2, 0.3],
                pressure_pa: 0.4,
                gradient_pa_m: [0.05, -0.08, 0.1],
            },
            TetPressure::FreeSurface {
                origin_m: [0.0; 3],
                pressure_pa: 0.13,
                gradient_pa_m: [0.0, 0.0, -0.4],
            },
        ];
        let scales = HyperTetPathScales {
            displacement_m: 0.4,
            force_n: 2.3,
        };
        let q = [0.27, 0.06, -0.14, -0.04, 0.05, 0.18];
        let direction = [0.7, -0.2, 0.3, 0.15, -0.6, 0.4];
        let h = 2e-6;
        for pressure in fields {
            let path = p
                .prepare_continuation_with_loads(
                    scales,
                    HyperTetPathLoads {
                        pressure: Some((&surface, pressure)),
                        gravity: Some((ReferenceGravityDensity::Nodal(&density), &affine_gravity)),
                    },
                    cx,
                )
                .unwrap();
            assert_eq!(path.free_dofs(), &[3, 6, 7, 9, 10, 11]);
            let mut u = vec![0.0; 12];
            for (&dof, &value) in path.free_dofs().iter().zip(&q) {
                u[dof] = scales.displacement_m * value;
            }
            let pressure_force = surface.force_n(&u, pressure, cx).unwrap();
            let gravity = p
                .gravity_load_with_nodal_density(&density, &u, &affine_gravity, cx)
                .unwrap();
            let mut dead: Vec<f64> = (0..12).map(|d| p.body_force_n_m3[d % 3] / 24.0).collect();
            for &(dof, force) in &forces {
                dead[dof] += force;
            }
            for (&actual, &expected) in path.reference_force_n().iter().zip(&dead) {
                near(actual, expected, 2e-15);
            }
            for lambda in [-0.6, 0.0, 0.7] {
                let derivative = path.load_vector_at(&q, lambda, cx).unwrap();
                let rp = path.residual(&q, lambda + h, cx).unwrap();
                let rm = path.residual(&q, lambda - h, cx).unwrap();
                for (i, &dof) in path.free_dofs().iter().enumerate() {
                    let expected =
                        (dead[dof] + pressure_force[dof] + gravity.force_n[dof]) / scales.force_n;
                    near(derivative[i], expected, 2e-12);
                    near((rm[i] - rp[i]) / (2.0 * h), expected, 2e-8);
                }
                let plus: Vec<_> = q.iter().zip(direction).map(|(q, d)| q + h * d).collect();
                let minus: Vec<_> = q.iter().zip(direction).map(|(q, d)| q - h * d).collect();
                let rp = path.residual(&plus, lambda, cx).unwrap();
                let rm = path.residual(&minus, lambda, cx).unwrap();
                let tangent = path.tangent(&q, lambda, cx).unwrap();
                let mut action = vec![0.0; q.len()];
                tangent.spmv(&direction, &mut action);
                for i in 0..q.len() {
                    near((rp[i] - rm[i]) / (2.0 * h), action[i], 3e-8);
                }
                if matches!(pressure, TetPressure::PerFace(_)) && lambda != 0.0 {
                    assert!((0..q.len()).any(|i| {
                        (0..i).any(|j| (tangent.get(i, j) - tangent.get(j, i)).abs() > 1e-8)
                    }));
                }
            }
        }
    });
}

#[test]
fn accepted_loaded_path_matches_static_owner_and_complete_physical_reactions() {
    with_cx(|_, cx| {
        let x = nodes(1.0);
        let (mu, lame) = (2.4, 6.0);
        let law =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda: lame }, 5.0).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 11).map(|d| (d, 0.0)).collect();
        let forces = [(11, 0.05), (0, 0.17)];
        let cell_density = [6.0];
        let nodal_density = [3.0, 3.0, 3.0, 15.0];
        let surface = ClosedTetPressureSurface::new(&x, &TETS, Default::default(), 4, cx).unwrap();
        let (p0, dp, gradient, initial_lambda) = (0.18, 0.12, -0.8, 0.7);
        let initial_stretch: f64 = 1.2;
        let internal = |s: f64| (mu * (s - 1.0 / s) + lame * s.ln() / s) / 6.0;
        for graded in [false, true] {
            // Independent simplex moments: total mass=1, share=integral rho*N3,
            // mass=integral rho*N3^2. Grading changes both retained moments.
            let (share, mass) = if graded { (0.325, 0.15) } else { (0.25, 0.1) };
            let constant = (internal(initial_stretch) / initial_lambda - (0.05 + 0.12 / 24.0)
                + p0 / 6.0
                + dp * initial_stretch / 12.0
                - gradient * mass * initial_stretch)
                / share;
            let field = |x: [f64; 3]| {
                Ok(GravitySample3 {
                    potential_j_kg: -constant * x[2] - 0.5 * gradient * x[2] * x[2],
                    acceleration_m_s2: [0.0, 0.0, constant + gradient * x[2]],
                    acceleration_gradient_s2: [[0.0; 3], [0.0; 3], [0.0, 0.0, gradient]],
                })
            };
            let field: &GravityField3<'_> = &field;
            let mut p = problem(&x, &law, &fixed, &forces);
            p.body_force_n_m3 = [0.0, 0.0, 0.12];
            let density = if graded {
                ReferenceGravityDensity::Nodal(&nodal_density)
            } else {
                ReferenceGravityDensity::PerElement(&cell_density)
            };
            let pressure = TetPressure::Affine {
                origin_m: [0.0; 3],
                pressure_pa: p0,
                gradient_pa_m: [0.0, 0.0, dp],
            };
            let path = p
                .prepare_continuation_with_loads(
                    HyperTetPathScales {
                        displacement_m: 1.0,
                        force_n: 1.0,
                    },
                    HyperTetPathLoads {
                        pressure: Some((&surface, pressure)),
                        gravity: Some((density, field)),
                    },
                    cx,
                )
                .unwrap();
            let (held_density, held_field) = path.held_loads().gravity.unwrap();
            assert!(core::ptr::eq(held_field, field));
            match held_density {
                ReferenceGravityDensity::Nodal(rho) => {
                    assert!(graded);
                    assert_eq!(rho, &nodal_density);
                }
                ReferenceGravityDensity::PerElement(rho) => {
                    assert!(!graded);
                    assert_eq!(rho, &cell_density);
                }
            }
            let controls = settings();
            let mut initial_u = [0.0; 12];
            initial_u[11] = initial_stretch - 1.0;
            let initial = path
                .start_from(&initial_u, initial_lambda, &controls, cx)
                .unwrap();
            let state = path.advance(&initial, &controls, 4, cx).unwrap();
            let response = path.response(&state, cx).unwrap();
            let s = 1.0 + response.displacement_m[3][2];
            assert!(state.lambda > initial_lambda);
            let applied_free = 0.05 + 0.12 / 24.0 - p0 / 6.0 - dp * s / 12.0
                + constant * share
                + gradient * mass * s;
            near(internal(s), state.lambda * applied_free, 2e-10);
            // Complete force balance includes supported dead force, all reference
            // gravity mass and the closed-surface affine-pressure resultant.
            let total = [
                0.17,
                0.0,
                0.05 + 0.12 / 6.0 + constant + gradient * s * share - dp * s / 6.0,
            ];
            for (axis, force) in total.into_iter().enumerate() {
                let reaction: f64 = response
                    .reactions_n
                    .iter()
                    .filter(|(dof, _)| dof % 3 == axis)
                    .map(|(_, f)| f)
                    .sum();
                near(reaction + state.lambda * force, 0.0, 2e-9);
            }
            let scaled_forces = forces.map(|(dof, f)| (dof, state.lambda * f));
            let scaled_field = |x| {
                let sample = field(x)?;
                Ok(GravitySample3 {
                    potential_j_kg: state.lambda * sample.potential_j_kg,
                    acceleration_m_s2: sample.acceleration_m_s2.map(|g| state.lambda * g),
                    acceleration_gradient_s2: sample
                        .acceleration_gradient_s2
                        .map(|row| row.map(|g| state.lambda * g)),
                })
            };
            let scaled_pressure = TetPressure::Affine {
                origin_m: [0.0; 3],
                pressure_pa: state.lambda * p0,
                gradient_pa_m: [0.0, 0.0, state.lambda * dp],
            };
            let loaded = HyperTetProblem {
                nodal_forces_n: &scaled_forces,
                body_force_n_m3: p.body_force_n_m3.map(|f| state.lambda * f),
                ..p
            };
            let native = if graded {
                loaded.solve_with_nodal_density_gravity_and_pressure(
                    &nodal_density,
                    &scaled_field,
                    &surface,
                    scaled_pressure,
                    cx,
                )
            } else {
                loaded.solve_with_gravity_and_pressure(
                    &cell_density,
                    &scaled_field,
                    &surface,
                    scaled_pressure,
                    cx,
                )
            }
            .unwrap();
            for (actual, expected) in response
                .displacement_m
                .iter()
                .flatten()
                .zip(native.displacement_m.iter().flatten())
            {
                near(*actual, *expected, 2e-9);
            }
            for ((dof, actual), (other_dof, expected)) in
                response.reactions_n.iter().zip(&native.reactions_n)
            {
                assert_eq!(dof, other_dof);
                near(*actual, *expected, 2e-9);
            }
            near(response.strain_energy_j, native.strain_energy_j, 2e-10);
        }
    });
}

#[test]
fn load_binding_and_trial_domain_refusals_leave_the_input_checkpoint_unchanged() {
    with_cx(|_, cx| {
        let x = nodes(1.0);
        let law = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 3.2,
                lambda: 8.0,
            },
            5.0,
        )
        .unwrap();
        let p = problem(&x, &law, &DIAGONAL_SUPPORTS, &[]);
        let scales = HyperTetPathScales {
            displacement_m: 1.0,
            force_n: 1.0,
        };
        let wrong_x = nodes(1.1);
        let wrong_surface =
            ClosedTetPressureSurface::new(&wrong_x, &TETS, p.budget, 4, cx).unwrap();
        assert!(
            p.prepare_continuation_with_loads(
                scales,
                HyperTetPathLoads {
                    pressure: Some((&wrong_surface, TetPressure::Uniform(0.0))),
                    gravity: None,
                },
                cx
            )
            .is_err()
        );
        let surface = ClosedTetPressureSurface::new(&x, &TETS, p.budget, 4, cx).unwrap();
        assert!(
            p.prepare_continuation_with_loads(
                scales,
                HyperTetPathLoads {
                    pressure: Some((&surface, TetPressure::PerFace(&[0.2; 3]))),
                    gravity: None,
                },
                cx
            )
            .is_err()
        );
        for density in [
            ReferenceGravityDensity::PerElement(&[]),
            ReferenceGravityDensity::PerElement(&[0.0]),
            ReferenceGravityDensity::Nodal(&[1.0; 3]),
            ReferenceGravityDensity::Nodal(&[1.0, 1.0, f64::NAN, 1.0]),
        ] {
            assert!(
                p.prepare_continuation_with_loads(
                    scales,
                    HyperTetPathLoads {
                        pressure: None,
                        gravity: Some((density, &affine_gravity))
                    },
                    cx
                )
                .is_err()
            );
        }
        let density = [3.0, 3.0, 3.0, 15.0];
        let bounded = |x: [f64; 3]| {
            if x[2] > 1.01 {
                return Err(HyperTetError::InvalidInput {
                    what: "outside original gravity domain".into(),
                });
            }
            Ok(GravitySample3 {
                potential_j_kg: -x[2],
                acceleration_m_s2: [0.0, 0.0, 1.0],
                acceleration_gradient_s2: [[0.0; 3]; 3],
            })
        };
        let path = p
            .prepare_continuation_with_loads(
                scales,
                HyperTetPathLoads {
                    pressure: None,
                    gravity: Some((ReferenceGravityDensity::Nodal(&density), &bounded)),
                },
                cx,
            )
            .unwrap();
        let controls = ArcSettings {
            ds: 0.5,
            ds_min: 0.5,
            ds_max: 0.5,
            ..settings()
        };
        let initial = path.start(&controls, cx).unwrap();
        let saved = initial.clone();
        assert!(matches!(path.advance(&initial, &controls, 1, cx),
            Err(HyperTetContinuationError::Solid(HyperTetError::InvalidInput { what }))
                if what == "outside original gravity domain"));
        same_state(&initial, &saved);
        assert!(path.load_vector_at(&[0.0; 2], 0.0, cx).is_err());
        assert!(path.load_vector_at(&[0.0, 0.0, 0.1], 0.0, cx).is_err());
        let tight = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 3.2,
                lambda: 8.0,
            },
            1.01,
        )
        .unwrap();
        let tight_path = problem(&x, &tight, &DIAGONAL_SUPPORTS, &[])
            .prepare_continuation_with_loads(
                scales,
                HyperTetPathLoads {
                    pressure: Some((&surface, TetPressure::Uniform(-8.0))),
                    gravity: None,
                },
                cx,
            )
            .unwrap();
        let initial = tight_path.start(&controls, cx).unwrap();
        let saved = initial.clone();
        assert!(matches!(
            tight_path.advance(&initial, &controls, 1, cx),
            Err(HyperTetContinuationError::Solid(
                HyperTetError::MaterialRefused { element: 0, .. }
            ))
        ));
        same_state(&initial, &saved);
    });
}

#[test]
fn cancellation_inside_a_retained_gravity_callback_publishes_no_path_step() {
    with_cx(|gate, cx| {
        let x = nodes(1.0);
        let law = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 3.2,
                lambda: 8.0,
            },
            5.0,
        )
        .unwrap();
        let armed = Cell::new(false);
        let calls = Cell::new(0_usize);
        let field = |x| {
            calls.set(calls.get() + 1);
            if armed.get() {
                gate.request();
            }
            affine_gravity(x)
        };
        let density = [2.0];
        let path = problem(&x, &law, &DIAGONAL_SUPPORTS, &[])
            .prepare_continuation_with_loads(
                HyperTetPathScales {
                    displacement_m: 1.0,
                    force_n: 1.0,
                },
                HyperTetPathLoads {
                    pressure: None,
                    gravity: Some((ReferenceGravityDensity::PerElement(&density), &field)),
                },
                cx,
            )
            .unwrap();
        let controls = settings();
        let initial = path.start(&controls, cx).unwrap();
        let saved = initial.clone();
        let before = calls.get();
        armed.set(true);
        assert!(matches!(
            path.advance(&initial, &controls, 4, cx),
            Err(HyperTetContinuationError::Solid(HyperTetError::Cancelled))
        ));
        assert!(calls.get() > before);
        same_state(&initial, &saved);
    });
}
