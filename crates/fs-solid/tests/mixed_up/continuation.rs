//! Native mixed path regressions. These remain spatial discretization checks,
//! not a general stability, inf-sup or material-validation certificate.
use super::*;
use fs_material::hyper::{OgdenTerm, Volumetric};
use fs_solid::continuation::{ArcSettings, PathEvent};
use fs_solid::hyper3::gravity::GravitySample3;
use fs_solid::hyper3::mixed::continuation::MixedContinuationError;
use fs_solid::hyper3::mixed::pressure::UniformBoundaryPressure;

fn arc(ds: f64) -> ArcSettings {
    ArcSettings {
        ds,
        ds_min: 1e-7,
        ds_max: ds,
        max_corrector_iters: 16,
        tol: 5e-9,
    }
}
fn maximum(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |m, value| m.max(value.abs()))
}

#[test]
fn g1_incompressible_ogden_traction_path_crosses_analytic_fold_and_resumes() {
    with_cx(|cx| {
        // One positive Ogden term with alpha<1 has a tensile load maximum.
        // P(s)=s^(-1/2)-s^(-5/4), so P'(s)=0 at s=(5/2)^(4/3).
        // This is a declared constitutive fixture, not a fitted material card.
        let material = Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![OgdenTerm {
                    mu_pa: 1.0,
                    alpha: 0.5,
                }],
                volumetric: Volumetric::Quadratic { kappa: 80.0 },
            },
            8.0,
        )
        .unwrap();
        let (p1, tets) = block([1, 1, 1], [1.0; 3]);
        let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let supports = rollers(&nodes);
        let forces = face_load(&nodes, &elements, |x| x[0] == 1.0, [1.0, 0.0, 0.0]);
        let native = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &forces,
            VolumetricTreatment::ExactIncompressible,
        );
        let path = native.prepare_continuation(cx).unwrap();
        let initial_stretch = 2.9_f64;
        let lateral = initial_stretch.sqrt().recip();
        let u: Vec<_> = nodes
            .iter()
            .map(|x| {
                [
                    (initial_stretch - 1.0) * x[0],
                    (lateral - 1.0) * x[1],
                    (lateral - 1.0) * x[2],
                ]
            })
            .collect();
        let p = (initial_stretch.sqrt() - initial_stretch.powf(-0.25)) / 3.0;
        let load = initial_stretch.powf(-0.5) - initial_stretch.powf(-1.25);
        let controls = arc(0.04);
        let start = path
            .start_from(
                &u,
                &vec![p; path.pressure_nodes().len()],
                load,
                &controls,
                cx,
            )
            .unwrap();
        let straight = path.advance(&start, &controls, 80, cx).unwrap();
        let half = path.advance(&start, &controls, 40, cx).unwrap();
        let resumed = path.advance(&half, &controls, 40, cx).unwrap();
        assert_eq!(straight.u, resumed.u);
        assert_eq!(straight.lambda, resumed.lambda);
        assert_eq!(straight.tangent, resumed.tangent);
        assert_eq!(straight.trace, resumed.trace);
        assert_eq!(straight.events, resumed.events);
        assert!(
            straight
                .events
                .iter()
                .any(|e| matches!(e, PathEvent::LimitPoint { .. }))
        );
        let response = path.response(&straight, cx).unwrap();
        let right = nodes.iter().position(|x| *x == [1.0, 0.0, 0.0]).unwrap();
        let stretch = 1.0 + response.displacement_m[right][0];
        let fold = 2.5_f64.powf(4.0 / 3.0);
        assert!(
            stretch > fold,
            "path must reach the descending tensile branch"
        );
        near(
            straight.lambda,
            stretch.powf(-0.5) - stretch.powf(-1.25),
            3e-7,
        );
        assert!(straight.lambda < straight.trace.iter().map(|x| x.0).fold(0.0_f64, f64::max));
        let lateral = stretch.sqrt().recip();
        for (x, actual) in nodes.iter().zip(&response.displacement_m) {
            for i in 0..3 {
                let expected = x[i] * (if i == 0 { stretch } else { lateral } - 1.0);
                near(actual[i], expected, 3e-6);
            }
        }
        let expected_pressure = (stretch.sqrt() - stretch.powf(-0.25)) / 3.0;
        for &pressure in &response.pressure_pa {
            near(pressure, expected_pressure, 3e-6);
        }
        near(response.min_det_f, 1.0, 3e-7);
        assert!(maximum(&response.free_force_residual_n) < 1e-7);
        assert!(maximum(&response.pressure_residual_m3) < 1e-8);
        let reaction: f64 = response
            .reactions_n
            .iter()
            .filter(|(d, _)| d % 3 == 0)
            .map(|(_, r)| r)
            .sum();
        near(reaction, -straight.lambda, 3e-7);
    });
}

#[test]
fn g0_full_mixed_tangent_and_parameter_column_retain_moving_supports_and_native_loads() {
    with_cx(|cx| {
        let material = law(2.0, 80.0);
        let (p1, tets) = block([1, 1, 1], [1.0; 3]);
        let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let supports: Vec<_> = nodes
            .iter()
            .enumerate()
            .filter(|(_, x)| x[0] == 0.0)
            .flat_map(|(node, x)| {
                [
                    (3 * node, 0.03 * x[2]),
                    (3 * node + 1, -0.02 * x[1]),
                    (3 * node + 2, 0.01 * x[2]),
                ]
            })
            .collect();
        let forces = face_load(&nodes, &elements, |x| x[0] == 1.0, [0.12, -0.03, 0.02]);
        let densities = vec![1.7; elements.len()];
        let gravity = |x: [f64; 3]| {
            Ok(GravitySample3 {
                potential_j_kg: 0.15 * x[0] * x[0] - 0.1 * x[1] * x[1]
                    + 0.1 * x[2]
                    + 0.2 * x[2] * x[2],
                acceleration_m_s2: [-0.3 * x[0], 0.2 * x[1], -0.1 - 0.4 * x[2]],
                acceleration_gradient_s2: [[-0.3, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, -0.4]],
            })
        };
        for volumetric in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            let native = problem(&nodes, &elements, &material, &supports, &forces, volumetric);
            let held = native
                .prepare_with_gravity(
                    &densities,
                    &gravity,
                    Some(UniformBoundaryPressure::new(0.09).unwrap()),
                    cx,
                )
                .unwrap();
            let path = held.into_continuation(cx).unwrap();
            let nu = 3 * nodes.len();
            let q: Vec<_> = path
                .free_dofs()
                .iter()
                .map(|&d| {
                    if d < nu {
                        0.015 * (d % 3 + 1) as f64 * nodes[d / 3][(d % 3 + 1) % 3]
                    } else {
                        0.2 + 0.01 * (d - nu) as f64
                    }
                })
                .collect();
            let lambda = 0.43;
            let k = path.tangent(&q, lambda, cx).unwrap();
            let h = 1e-6;
            for seed in [0.2, 0.7, 1.3] {
                let direction: Vec<_> = (0..q.len())
                    .map(|i| ((i + 1) as f64 * seed).sin())
                    .collect();
                let plus: Vec<_> = q.iter().zip(&direction).map(|(q, d)| q + h * d).collect();
                let minus: Vec<_> = q.iter().zip(&direction).map(|(q, d)| q - h * d).collect();
                let rp = path.residual(&plus, lambda, cx).unwrap();
                let rm = path.residual(&minus, lambda, cx).unwrap();
                let mut action = vec![0.0; q.len()];
                k.spmv(&direction, &mut action);
                let scale = maximum(&action).max(1.0);
                for i in 0..q.len() {
                    near((rp[i] - rm[i]) / (2.0 * h), action[i], 5e-7 * scale);
                }
            }
            let parameter = path.load_vector_at(&q, lambda, cx).unwrap();
            let rp = path.residual(&q, lambda + h, cx).unwrap();
            let rm = path.residual(&q, lambda - h, cx).unwrap();
            let scale = maximum(&parameter).max(1.0);
            for i in 0..q.len() {
                near(-(rp[i] - rm[i]) / (2.0 * h), parameter[i], 5e-7 * scale);
            }
            assert!(
                parameter
                    .iter()
                    .enumerate()
                    .any(|(i, v)| path.free_dofs()[i] >= nu && v.abs() > 1e-5),
                "moving supports must reach the pressure-constraint rows"
            );
            // Cancellation is checked even when q,lambda already has a cached tangent.
            path.tangent(&q, lambda, cx).unwrap();
            let gate = CancelGate::new();
            gate.request();
            with_gate(&gate, |cancelled| {
                assert!(matches!(
                    path.tangent(&q, lambda, cancelled),
                    Err(MixedContinuationError::Mixed(MixedUpError::Cancelled))
                ));
            });
        }
        let native = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &forces,
            VolumetricTreatment::OriginalMaterial,
        );
        let direct = native
            .prepare_with_uniform_boundary_pressure(0.09, cx)
            .unwrap()
            .into_continuation(cx)
            .unwrap();
        assert_eq!(direct.pressure_nodes().len(), 8);
        let mut tiny = native;
        tiny.settings.maximum_factor_bytes = 1;
        assert!(tiny.prepare_continuation(cx).is_err());
    });
}

#[test]
fn g0_moving_pressure_gauge_is_a_datum_and_omitted_volume_equation_cannot_disappear() {
    with_cx(|cx| {
        let material = law(2.0, 100.0);
        let (p1, tets) = block([2, 2, 2], [1.0; 3]);
        let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let boundary: Vec<_> = nodes
            .iter()
            .enumerate()
            .filter(|(_, x)| x.iter().any(|&v| v == 0.0 || v == 1.0))
            .flat_map(|(node, x)| {
                [
                    (3 * node, 0.04 * x[1]),
                    (3 * node + 1, 0.0),
                    (3 * node + 2, 0.0),
                ]
            })
            .collect();
        let mut native = problem(
            &nodes,
            &elements,
            &material,
            &boundary,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        native.pressure_gauge = Some((0, 0.3));
        let path = native.prepare_continuation(cx).unwrap();
        let controls = arc(0.1);
        let start = path.start(&controls, cx).unwrap();
        let next = path.advance(&start, &controls, 3, cx).unwrap();
        let response = path.response(&next, cx).unwrap();
        for (x, u) in nodes.iter().zip(&response.displacement_m) {
            near(u[0], next.lambda * 0.04 * x[1], 2e-7);
            near(u[1], 0.0, 2e-7);
            near(u[2], 0.0, 2e-7);
        }
        for p in &response.pressure_pa {
            near(*p, next.lambda * 0.3, 2e-6);
        }
        assert!(maximum(&response.pressure_residual_m3) < 1e-9);
        // The physical pressure datum must remain in the exact load column.
        let h = 1e-6;
        let load = path.load_vector_at(&next.u, next.lambda, cx).unwrap();
        let rp = path.residual(&next.u, next.lambda + h, cx).unwrap();
        let rm = path.residual(&next.u, next.lambda - h, cx).unwrap();
        for i in 0..load.len() {
            near(-(rp[i] - rm[i]) / (2.0 * h), load[i], 5e-7);
        }

        let expansion: Vec<_> = boundary
            .iter()
            .map(|&(d, _)| (d, 0.03 * nodes[d / 3][d % 3]))
            .collect();
        native.prescribed_m = &expansion;
        let incompatible = native.prepare_continuation(cx).unwrap();
        let previous = incompatible.start(&controls, cx).unwrap();
        assert!(matches!(
            incompatible.advance(&previous, &controls, 1, cx),
            Err(MixedContinuationError::Mixed(
                MixedUpError::IncompatiblePressureGauge { node: 0 }
            ))
        ));
        let mut malformed = previous.clone();
        malformed.u.push(0.0);
        assert!(incompatible.advance(&malformed, &controls, 0, cx).is_err());
        malformed = previous.clone();
        malformed.tangent.0.push(0.0);
        assert!(incompatible.advance(&malformed, &controls, 0, cx).is_err());
        assert_eq!(previous.step, 0);
        assert_eq!(previous.lambda, 0.0);
        assert!(previous.u.iter().all(|v| *v == 0.0));

        // A traction boundary fixes pressure physically and cannot use this gauge.
        let supports = rollers(&nodes);
        native.prescribed_m = &supports;
        let bad_gauge = native.prepare_continuation(cx).unwrap();
        assert!(matches!(
            bad_gauge.start(&controls, cx),
            Err(MixedContinuationError::Mixed(
                MixedUpError::IncompatiblePressureGauge { node: 0 }
            ))
        ));
    });
}
