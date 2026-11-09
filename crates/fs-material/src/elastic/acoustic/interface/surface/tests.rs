use super::*;

fn material(scale: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        2500.0 * scale,
        2500.0 * (4000.0_f64.powi(2) - (4.0 / 3.0) * 2000.0_f64.powi(2)) * scale,
        2500.0 * 2000.0_f64.powi(2) * scale,
        ContentHash([19; 32]),
        0.002,
    )
    .unwrap()
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

const MODES: [ElasticMode; 3] = [ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh];
const CONDITIONS: [ElasticSurfaceCondition; 2] = [
    ElasticSurfaceCondition::TractionFree,
    ElasticSurfaceCondition::Clamped,
];

// G1: phase is meaningful in the shared modal basis. A free normal P wave
// reverses its modal sign while doubling its Cartesian normal velocity.
#[test]
fn g1_normal_free_and_clamped_reflection_keep_phase_and_reaction() {
    let state = material(1.0);
    for condition in CONDITIONS {
        let owner = IsotropicElasticSurface::new(&state, condition).unwrap();
        assert!(core::ptr::eq(owner.state(), &state));
        let prepared = owner.prepare_incidence(ElasticMode::P, 1.0).unwrap();
        assert_eq!(prepared.condition(), condition);
        for (i, mode) in MODES.into_iter().enumerate() {
            let sign = if i == 0 { -1.0 } else { 1.0 }
                * if condition == ElasticSurfaceCondition::TractionFree {
                    1.0
                } else {
                    -1.0
                };
            let coefficients = prepared.velocity_coefficients(mode).unwrap();
            for (j, coefficient) in coefficients.into_iter().enumerate() {
                close(
                    coefficient,
                    C64::from_re(if i == j { sign } else { 0.0 }),
                    2e-14,
                );
            }
            let mut incoming = [C64::ZERO; 3];
            incoming[i] = C64::new(0.003, -0.004);
            let response = prepared.drive(incoming).unwrap();
            assert_eq!(
                response.material_state_identity,
                state.material_state_identity()
            );
            let wave = prepared.incoming_mode(mode).unwrap();
            for j in 0..3 {
                let (velocity, traction) = if condition == ElasticSurfaceCondition::TractionFree {
                    (wave.velocity[j] * incoming[i].scale(2.0), C64::ZERO)
                } else {
                    (C64::ZERO, wave.traction_pa_s_m[j] * incoming[i].scale(2.0))
                };
                close(response.boundary_velocity_m_s[j], velocity, 1e-15);
                close(response.boundary_traction_pa[j], traction, 2e-9);
            }
            assert!(
                (response.total_incoming_power_w_m2 - response.total_reflected_power_w_m2).abs()
                    < 1e-10
            );
        }
    }
}

// Independent closed-form half-space coefficients, not the production
// traction solve. Both formulas use the documented Cartesian polarizations.
fn analytic(condition: ElasticSurfaceCondition, cp: f64, cs: f64, p: f64) -> [[C64; 2]; 2] {
    let sp = cp * p;
    let ss = cs * p;
    let qp = C64::from_re(1.0 - sp * sp).sqrt();
    let qs = C64::from_re(1.0 - ss * ss).sqrt();
    match condition {
        ElasticSurfaceCondition::TractionFree => {
            let d = 1.0 - 2.0 * ss * ss;
            let conversion = (qp * qs).scale(4.0 * cs.powi(3) / cp * p * p);
            let denominator = C64::from_re(d * d) + conversion;
            [
                [
                    (conversion - C64::from_re(d * d)) / denominator,
                    qs.scale(4.0 * cs * cs / cp * p * d) / denominator,
                ],
                [
                    qp.scale(4.0 * cs * p * d) / denominator,
                    (C64::from_re(d * d) - conversion) / denominator,
                ],
            ]
        }
        ElasticSurfaceCondition::Clamped => {
            let denominator = qp * qs + C64::from_re(sp * ss);
            [
                [
                    (qp * qs - C64::from_re(sp * ss)) / denominator,
                    qs.scale(-2.0 * ss) / denominator,
                ],
                [
                    qp.scale(-2.0 * sp) / denominator,
                    (C64::from_re(sp * ss) - qp * qs) / denominator,
                ],
            ]
        }
        ElasticSurfaceCondition::Impedance(_) | ElasticSurfaceCondition::TensorImpedance(_) => {
            unreachable!("ideal boundary fixture")
        }
    }
}

#[test]
fn g1_oblique_conversion_matches_closed_forms_across_critical_angle() {
    let state = material(1.0);
    for condition in CONDITIONS {
        let owner = IsotropicElasticSurface::new(&state, condition).unwrap();
        let [cp, cs] = owner.phase_speeds_m_s();
        for incidence in [ElasticMode::P, ElasticMode::Sv] {
            for cosine in [0.999, 0.95, 0.87, 0.86, 0.7, 0.5, 0.01] {
                let prepared = owner.prepare_incidence(incidence, cosine).unwrap();
                let expected = analytic(condition, cp, cs, prepared.tangential_slowness_s_m());
                for (column, mode) in [ElasticMode::P, ElasticMode::Sv].into_iter().enumerate() {
                    if let Some(coefficients) = prepared.velocity_coefficients(mode) {
                        for row in 0..2 {
                            close(coefficients[row], expected[row][column], 3e-12);
                        }
                        assert_eq!(coefficients[2], C64::ZERO);
                    }
                }
            }
        }
    }
}

// G0/G1: all propagating input channels share one conserved slowness. Their
// coherent combination closes work without an incoherent energy replacement.
#[test]
fn g0_coherent_reflection_is_power_unitary_and_satisfies_the_boundary() {
    let state = material(1.0);
    for condition in CONDITIONS {
        let owner = IsotropicElasticSurface::new(&state, condition).unwrap();
        for mode in [ElasticMode::P, ElasticMode::Sv] {
            for cosine in [1.0, 0.95, 0.8, 0.3, 1e-6] {
                let prepared = owner.prepare_incidence(mode, cosine).unwrap();
                let incoming = core::array::from_fn(|i| {
                    if prepared.incoming_mode(MODES[i]).is_some() {
                        C64::new((i + 1) as f64 * 0.001, (3 - i) as f64 * -0.0007)
                    } else {
                        C64::ZERO
                    }
                });
                let response = prepared.drive(incoming).unwrap();
                let power = response.total_incoming_power_w_m2;
                assert!((response.total_reflected_power_w_m2 / power - 1.0).abs() < 2e-11);
                assert!(response.boundary_normal_flux_w_m2.abs() / power < 2e-11);
                let constrained = match condition {
                    ElasticSurfaceCondition::TractionFree => response.boundary_traction_pa,
                    ElasticSurfaceCondition::Clamped => response.boundary_velocity_m_s,
                    ElasticSurfaceCondition::Impedance(_)
                    | ElasticSurfaceCondition::TensorImpedance(_) => {
                        unreachable!("ideal boundary fixture")
                    }
                };
                let scale = match condition {
                    ElasticSurfaceCondition::TractionFree => {
                        state.density_kg_m3() * owner.speeds[0]
                    }
                    ElasticSurfaceCondition::Clamped => 1.0,
                    ElasticSurfaceCondition::Impedance(_)
                    | ElasticSurfaceCondition::TensorImpedance(_) => {
                        unreachable!("ideal boundary fixture")
                    }
                };
                assert!(constrained.iter().all(|x| x.abs() < 2e-14 * scale));
                // Independent weighted inner products of the retained columns.
                for a in 0..3 {
                    for b in 0..3 {
                        if let (Some(ca), Some(cb)) = (
                            prepared.velocity_coefficients(MODES[a]),
                            prepared.velocity_coefficients(MODES[b]),
                        ) {
                            let product = (0..3).fold(C64::ZERO, |sum, i| {
                                sum + (ca[i].conj() * cb[i])
                                    .scale(prepared.reflected[i].normal_power_factor_pa_s_m)
                            });
                            let expected = if a == b {
                                prepared.incoming[a].normal_power_factor_pa_s_m
                            } else {
                                0.0
                            };
                            close(
                                product,
                                C64::from_re(expected),
                                3e-12 * state.density_kg_m3() * owner.speeds[0],
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn g1_supercritical_shear_keeps_decaying_p_field_and_reflection_phase() {
    let state = material(1.0);
    for condition in CONDITIONS {
        let prepared = IsotropicElasticSurface::new(&state, condition)
            .unwrap()
            .prepare_incidence(ElasticMode::Sv, 0.5)
            .unwrap();
        let p = prepared.reflected_mode(ElasticMode::P);
        assert!(p.normal_slowness_s_m.im < 0.0); // decays toward z -> -infinity
        assert_eq!(p.normal_power_factor_pa_s_m, 0.0);
        assert!(prepared.incoming_mode(ElasticMode::P).is_none());
        let coefficients = prepared.velocity_coefficients(ElasticMode::Sv).unwrap();
        assert!(coefficients[0].abs() > 0.1);
        assert!(coefficients[1].im.abs() > 0.1);
        assert!((coefficients[1].norm_sq() - 1.0).abs() < 1e-13);
        let response = prepared
            .drive([C64::ZERO, C64::from_re(0.002), C64::ZERO])
            .unwrap();
        assert_eq!(response.reflected_normal_power_w_m2[0], 0.0);
        assert!(
            prepared
                .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
                .is_err()
        );
    }
}

#[test]
fn g3_common_material_scaling_preserves_reflection_and_scales_power() {
    let first = material(1.0);
    let second = material(100.0);
    for condition in CONDITIONS {
        let a = IsotropicElasticSurface::new(&first, condition)
            .unwrap()
            .prepare_incidence(ElasticMode::P, 0.7)
            .unwrap();
        let b = IsotropicElasticSurface::new(&second, condition)
            .unwrap()
            .prepare_incidence(ElasticMode::P, 0.7)
            .unwrap();
        let input = [C64::new(0.002, -0.003), C64::new(-0.001, 0.002), C64::ZERO];
        let ra = a.drive(input).unwrap();
        let rb = b.drive(input).unwrap();
        for i in 0..3 {
            close(
                ra.reflected_velocity_m_s[i],
                rb.reflected_velocity_m_s[i],
                1e-14,
            );
        }
        assert!(
            (rb.total_reflected_power_w_m2 / ra.total_reflected_power_w_m2 - 100.0).abs() < 2e-10
        );
    }
}

#[test]
fn g0_g4_source_domains_and_cancelled_preparation_refuse() {
    let state = material(1.0);
    let anisotropic = ElasticAcousticState::try_new_mandel(
        *state.stiffness_mandel_pa(),
        state.density_kg_m3(),
        ContentHash([20; 32]),
        0.002,
    )
    .unwrap();
    for condition in CONDITIONS {
        assert!(IsotropicElasticSurface::new(&anisotropic, condition).is_err());
        let owner = IsotropicElasticSurface::new(&state, condition).unwrap();
        for cosine in [0.0, -0.1, 1.1, f64::NAN, f64::INFINITY] {
            assert!(owner.prepare_incidence(ElasticMode::P, cosine).is_err());
        }
        for stop in [1, 3, 5] {
            let mut polls = 0;
            assert!(
                owner
                    .prepare_incidence_cancellable(ElasticMode::P, 0.8, &mut || {
                        polls += 1;
                        polls == stop
                    })
                    .is_err()
            );
            assert_eq!(polls, stop);
        }
        let prepared = owner.prepare_incidence(ElasticMode::P, 0.8).unwrap();
        for bad in [C64::new(f64::NAN, 0.0), C64::from_re(1e6)] {
            assert!(prepared.drive([bad, C64::ZERO, C64::ZERO]).is_err());
        }
        assert_eq!(
            prepared
                .drive([C64::ZERO; 3])
                .unwrap()
                .total_reflected_power_w_m2,
            0.0
        );
    }
}

fn support(normal: C64, tangent: C64) -> ElasticSurfaceCondition {
    ElasticSurfaceCondition::Impedance(ElasticSurfaceImpedance::new(normal, tangent).unwrap())
}

// G1: all three channels decouple at normal incidence. Each reflection follows
// its scalar traction/velocity equation, including the P polarization reversal.
#[test]
fn g1_passive_support_matches_normal_scalar_reflection_and_loss() {
    let state = material(1.0);
    let zp = 2500.0 * 4000.0;
    let zs = 2500.0 * 2000.0;
    for (normal, tangent) in [
        (C64::from_re(zp), C64::from_re(zs)),
        (C64::new(0.4 * zp, 2.0 * zp), C64::new(3.0 * zs, -0.5 * zs)),
        (C64::new(0.0, -zp), C64::new(0.0, 2.0 * zs)),
        (C64::ZERO, C64::ZERO),
    ] {
        let prepared = IsotropicElasticSurface::new(&state, support(normal, tangent))
            .unwrap()
            .prepare_incidence(ElasticMode::P, 1.0)
            .unwrap();
        for (index, mode) in MODES.into_iter().enumerate() {
            let z = if index == 0 { normal } else { tangent };
            let reference = if index == 0 { zp } else { zs };
            let cartesian = (C64::from_re(reference) - z) / (C64::from_re(reference) + z);
            let modal = if index == 0 { -cartesian } else { cartesian };
            let coefficients = prepared.velocity_coefficients(mode).unwrap();
            for (output, coefficient) in coefficients.into_iter().enumerate() {
                close(
                    coefficient,
                    if output == index { modal } else { C64::ZERO },
                    3e-14,
                );
            }
            let mut incoming = [C64::ZERO; 3];
            incoming[index] = C64::new(0.002, -0.001);
            let response = prepared.drive(incoming).unwrap();
            let expected_velocity = incoming[index] * (C64::ONE + cartesian);
            let expected_loss = 0.5 * z.re * expected_velocity.norm_sq();
            assert!((response.absorbed_normal_power_w_m2 - expected_loss).abs() < 2e-13);
            if z.re == 0.0 {
                assert_eq!(response.absorbed_normal_power_w_m2, 0.0);
            }
        }
    }
}

// G1: the SH boundary sees rho*c_s*cos(theta), not normal-incidence rho*c_s.
#[test]
fn g1_oblique_support_sh_reflection_uses_the_normal_wave_impedance() {
    let state = material(1.0);
    let z = C64::new(1.2e6, -4.0e6);
    for cosine in [1.0, 0.9, 0.5, 0.01] {
        let prepared = IsotropicElasticSurface::new(&state, support(C64::from_re(1e7), z))
            .unwrap()
            .prepare_incidence(ElasticMode::Sh, cosine)
            .unwrap();
        let reference = C64::from_re(2500.0 * 2000.0 * cosine);
        let expected = (reference - z) / (reference + z);
        close(
            prepared.velocity_coefficients(ElasticMode::Sh).unwrap()[2],
            expected,
            3e-14,
        );
        let response = prepared
            .drive([C64::ZERO, C64::ZERO, C64::from_re(0.003)])
            .unwrap();
        assert!(response.absorbed_normal_power_w_m2 > 0.0);
        assert!(
            (response.total_incoming_power_w_m2
                - response.total_reflected_power_w_m2
                - response.absorbed_normal_power_w_m2)
                .abs()
                < 2e-13
        );
    }
}

// G0/G3: coherent mode conversion, supercritical P evanescence and reactive
// traction all use the original elastic tensor and close physical wall work.
#[test]
fn g0_coherent_support_reflection_closes_work_across_critical_angles() {
    let state = material(1.0);
    for resistance in [0.0, 1.0e-7, 1.0e6] {
        let normal = C64::new(resistance, 1.7e7);
        let tangent = C64::new(0.5 * resistance, -2.3e6);
        for mode in [ElasticMode::P, ElasticMode::Sv] {
            for cosine in [1.0, 0.95, 0.87, 0.86, 0.5, 0.05] {
                let prepared = IsotropicElasticSurface::new(&state, support(normal, tangent))
                    .unwrap()
                    .prepare_incidence(mode, cosine)
                    .unwrap();
                let incoming = core::array::from_fn(|i| {
                    if prepared.incoming_mode(MODES[i]).is_some() {
                        C64::new((i + 1) as f64 * 0.001, (3 - i) as f64 * -0.0007)
                    } else {
                        C64::ZERO
                    }
                });
                let response = prepared.drive(incoming).unwrap();
                let power = response.total_incoming_power_w_m2;
                assert!(
                    (power
                        - response.total_reflected_power_w_m2
                        - response.absorbed_normal_power_w_m2)
                        .abs()
                        < 3e-12 * power
                );
                assert!(
                    (response.boundary_normal_flux_w_m2 - response.absorbed_normal_power_w_m2)
                        .abs()
                        < 3e-12 * power
                );
                let mut loss = 0.0;
                for axis in 0..3 {
                    let z = if axis == 2 { normal } else { tangent };
                    let velocity = response.boundary_velocity_m_s[axis];
                    close(response.boundary_traction_pa[axis], -z * velocity, 1e-9);
                    loss += 0.5 * z.re * velocity.norm_sq();
                }
                assert!((loss - response.absorbed_normal_power_w_m2).abs() <= 2e-14 * loss);
                if resistance == 0.0 {
                    assert_eq!(response.absorbed_normal_power_w_m2, 0.0);
                } else {
                    assert!(response.absorbed_normal_power_w_m2 > 0.0);
                }
            }
        }
    }
}

#[test]
fn g1_stiff_support_retains_wall_motion_and_absorption_below_reflection_roundoff() {
    let state = material(1.0);
    let impedance = 1e30;
    let reference = 1e7;
    let input = 0.001;
    let prepared = IsotropicElasticSurface::new(
        &state,
        support(C64::from_re(impedance), C64::from_re(impedance)),
    )
    .unwrap()
    .prepare_incidence(ElasticMode::P, 1.0)
    .unwrap();
    let response = prepared
        .drive([C64::from_re(input), C64::ZERO, C64::ZERO])
        .unwrap();
    let velocity = 2.0 * input * (reference / impedance) / (1.0 + reference / impedance);
    let loss = 0.5 * impedance * velocity * velocity;
    assert!((response.boundary_velocity_m_s[2].re / velocity - 1.0).abs() < 2e-14);
    assert!((response.absorbed_normal_power_w_m2 / loss - 1.0).abs() < 2e-14);
    assert!(response.absorbed_normal_power_w_m2 > 0.0);
    assert_eq!(response.reflected_velocity_m_s[0], C64::from_re(input));
}

#[test]
fn g0_support_admission_rejects_active_and_nonfinite_values() {
    for bad in [
        C64::from_re(-1.0),
        C64::new(f64::NAN, 0.0),
        C64::new(0.0, f64::INFINITY),
    ] {
        assert!(ElasticSurfaceImpedance::new(bad, C64::ZERO).is_err());
        assert!(ElasticSurfaceImpedance::new(C64::ZERO, bad).is_err());
    }
    let normal = C64::new(1.0, 2.0);
    let tangent = C64::new(0.0, -3.0);
    let impedance = ElasticSurfaceImpedance::new(normal, tangent).unwrap();
    assert_eq!(impedance.normal_pa_s_m(), normal);
    assert_eq!(impedance.tangential_pa_s_m(), tangent);
}

mod tensor_support;
