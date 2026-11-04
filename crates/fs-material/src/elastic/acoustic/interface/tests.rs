use super::*;

fn material(density: f64, cp: f64, cs: f64, identity: u8) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        density,
        density * (cp * cp - (4.0 / 3.0) * cs * cs),
        density * cs * cs,
        ContentHash([identity; 32]),
        0.02,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:e}"
    );
}

fn complex_close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}

fn side(i: usize) -> ElasticSide {
    if i == 0 {
        ElasticSide::First
    } else {
        ElasticSide::Second
    }
}
fn mode(i: usize) -> ElasticMode {
    [ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh][i]
}

// Independent Cartesian displacement-gradient -> original Mandel tensor ->
// traction reconstruction. No production traction or flux coefficient is used.
fn tensor_traction(state: &ElasticAcousticState, wave: ElasticWaveTrace, p: f64) -> [C64; 3] {
    let slow = [C64::from_re(p), C64::ZERO, wave.normal_slowness_s_m];
    let strain: [[C64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (slow[i] * wave.velocity[j] + slow[j] * wave.velocity[i]).scale(-0.5)
        })
    });
    let e = [
        strain[0][0],
        strain[1][1],
        strain[2][2],
        strain[0][1].scale(core::f64::consts::SQRT_2),
        strain[1][2].scale(core::f64::consts::SQRT_2),
        strain[2][0].scale(core::f64::consts::SQRT_2),
    ];
    let stress: [C64; 6] = core::array::from_fn(|i| {
        (0..6).fold(C64::ZERO, |a, j| {
            a + e[j].scale(state.stiffness_mandel_pa()[i][j])
        })
    });
    [
        stress[5].scale(core::f64::consts::FRAC_1_SQRT_2),
        stress[4].scale(core::f64::consts::FRAC_1_SQRT_2),
        stress[2],
    ]
}

fn check_tensor_fields(
    prepared: &PreparedElasticInterface<'_>,
    response: &ElasticInterfaceResponse,
) {
    let p = prepared.tangential_slowness_s_m();
    for s in 0..2 {
        let mut traction = [C64::ZERO; 3];
        let mut velocity = [C64::ZERO; 3];
        for m in 0..3 {
            let i = 3 * s + m;
            for (wave, amplitude) in [
                (prepared.incoming[i], response.incoming_velocity_m_s[s][m]),
                (prepared.outgoing[i], response.outgoing_velocity_m_s[s][m]),
            ] {
                let stress = tensor_traction(prepared.states[s], wave, p);
                for c in 0..3 {
                    let scale = wave.traction_pa_s_m[c].abs().max(1.0);
                    complex_close(stress[c], wave.traction_pa_s_m[c], 3e-12 * scale);
                    traction[c] = traction[c] + stress[c] * amplitude;
                    velocity[c] = velocity[c] + wave.velocity[c] * amplitude;
                }
            }
        }
        let flux = -0.5
            * (0..3)
                .map(|c| (traction[c] * velocity[c].conj()).re)
                .sum::<f64>();
        close(
            flux,
            response.boundary_normal_flux_w_m2[s],
            3e-11 * response.total_incoming_power_w_m2.max(1e-30),
        );
        for c in 0..3 {
            complex_close(
                traction[c],
                response.boundary_traction_pa[s][c],
                1e-10 * traction[c].abs().max(1.0),
            );
            complex_close(velocity[c], response.boundary_velocity_m_s[s][c], 2e-12);
            complex_close(
                response.boundary_velocity_m_s[0][c],
                response.boundary_velocity_m_s[1][c],
                2e-12,
            );
            complex_close(
                response.boundary_traction_pa[0][c],
                response.boundary_traction_pa[1][c],
                3e-11 * response.boundary_traction_pa[0][c].abs().max(1.0),
            );
        }
    }
}

#[test]
fn g1_normal_incidence_retains_original_p_and_shear_impedances() {
    let first = material(2400.0, 4200.0, 2100.0, 1);
    let second = material(3100.0, 6100.0, 3300.0, 2);
    let owner = IsotropicElasticInterface::new(&first, &second).unwrap();
    let prepared = owner
        .prepare_incidence(ElasticSide::First, ElasticMode::P, 1.0)
        .unwrap();
    assert_eq!(
        prepared.material_state_identities(),
        [ContentHash([1; 32]), ContentHash([2; 32])]
    );
    assert!(core::ptr::eq(owner.state(ElasticSide::Second), &second));
    for s in 0..2 {
        for m in 0..3 {
            let zi = owner.impedances[s][m.min(1)];
            let zt = owner.impedances[1 - s][m.min(1)];
            let r = (zt - zi) / (zt + zi) * if m == 0 { 1.0 } else { -1.0 };
            let t = 2.0 * zi / (zt + zi);
            let coefficients = prepared.velocity_coefficients(side(s), mode(m)).unwrap();
            for output in 0..3 {
                complex_close(
                    coefficients[s][output],
                    C64::from_re(if output == m { r } else { 0.0 }),
                    1e-14,
                );
                complex_close(
                    coefficients[1 - s][output],
                    C64::from_re(if output == m { t } else { 0.0 }),
                    1e-14,
                );
            }
            let mut source = [[C64::ZERO; 3]; 2];
            source[s][m] = C64::new(0.003, -0.004);
            let response = prepared.drive(source).unwrap();
            let input = 0.5 * zi * 25e-6;
            close(response.total_incoming_power_w_m2, input, 1e-12 * input);
            close(
                response.outgoing_normal_power_w_m2[s][m] / input,
                r * r,
                2e-14,
            );
            close(
                response.outgoing_normal_power_w_m2[1 - s][m] / input,
                4.0 * zi * zt / (zi + zt).powi(2),
                2e-14,
            );
            check_tensor_fields(&prepared, &response);
        }
    }
}

#[test]
fn g0_identical_medium_is_exactly_transparent_including_grazing_shear() {
    let first = material(2600.0, 6000.0, 3000.0, 3);
    let second = material(2600.0, 6000.0, 3000.0, 4);
    let owner = IsotropicElasticInterface::new(&first, &second).unwrap();
    for cosine in [1.0, 0.9, 0.5, 1e-8, 1e-14] {
        for kind in [ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh] {
            let prepared = owner
                .prepare_incidence(ElasticSide::First, kind, cosine)
                .unwrap();
            for s in 0..2 {
                for m in 0..3 {
                    if let Some(column) = prepared.velocity_coefficients(side(s), mode(m)) {
                        assert_eq!(column[s], [C64::ZERO; 3]);
                        for j in 0..3 {
                            assert_eq!(column[1 - s][j], if m == j { C64::ONE } else { C64::ZERO });
                        }
                    }
                }
            }
            let mut source = [[C64::ZERO; 3]; 2];
            source[0][kind.index()] = C64::new(0.001, 0.002);
            let response = prepared.drive(source).unwrap();
            assert_eq!(response.outgoing_velocity_m_s[1], source[0]);
            assert_eq!(
                response.total_incoming_power_w_m2,
                response.total_outgoing_power_w_m2
            );
            if cosine > 1e-8 {
                check_tensor_fields(&prepared, &response);
            }
        }
    }
}

#[test]
fn g1_full_tensor_continuity_and_coherent_two_sided_work() {
    let first = material(2500.0, 4000.0, 2000.0, 5);
    let second = material(2800.0, 6200.0, 3400.0, 6);
    let owner = IsotropicElasticInterface::new(&first, &second).unwrap();
    for cosine in [0.999, 0.9, 0.8, 0.6, 0.3, 0.01] {
        for kind in [ElasticMode::P, ElasticMode::Sv] {
            let prepared = owner
                .prepare_incidence(ElasticSide::First, kind, cosine)
                .unwrap();
            let input = core::array::from_fn(|s| {
                core::array::from_fn(|m| {
                    if prepared.incoming_mode(side(s), mode(m)).is_some() {
                        C64::new((m + 1) as f64 * 0.001, (1 + s + m) as f64 * -0.0007)
                    } else {
                        C64::ZERO
                    }
                })
            });
            let result = prepared.drive(input).unwrap();
            close(
                result.total_incoming_power_w_m2,
                result.total_outgoing_power_w_m2,
                2e-12 * result.total_incoming_power_w_m2,
            );
            check_tensor_fields(&prepared, &result);
        }
    }
}

#[test]
fn g1_shear_total_reflection_keeps_evanescent_psv_phase_and_zero_normal_power() {
    let first = material(2500.0, 4000.0, 2000.0, 7);
    let second = material(2700.0, 6000.0, 3200.0, 8);
    let owner = IsotropicElasticInterface::new(&first, &second).unwrap();
    let prepared = owner
        .prepare_incidence(ElasticSide::First, ElasticMode::Sv, 0.5)
        .unwrap();
    let coefficients = prepared
        .velocity_coefficients(ElasticSide::First, ElasticMode::Sv)
        .unwrap();
    assert!(coefficients[0][0].abs() > 0.1);
    assert!(coefficients[0][1].im.abs() > 0.1);
    close(coefficients[0][1].norm_sq(), 1.0, 2e-14);
    assert!(coefficients[1][0].abs() > 0.1 && coefficients[1][1].abs() > 0.1);
    assert!(
        prepared
            .outgoing_mode(ElasticSide::First, ElasticMode::P)
            .normal_slowness_s_m
            .im
            < 0.0
    );
    assert!(
        prepared
            .outgoing_mode(ElasticSide::Second, ElasticMode::P)
            .normal_slowness_s_m
            .im
            > 0.0
    );
    assert!(
        prepared
            .incoming_mode(ElasticSide::First, ElasticMode::P)
            .is_none()
    );
    assert!(
        prepared
            .incoming_mode(ElasticSide::Second, ElasticMode::Sv)
            .is_none()
    );
    let result = prepared
        .drive([[C64::ZERO, C64::from_re(0.01), C64::ZERO], [C64::ZERO; 3]])
        .unwrap();
    assert_eq!(result.outgoing_normal_power_w_m2[0][0], 0.0);
    assert_eq!(result.outgoing_normal_power_w_m2[1], [0.0; 3]);
    close(
        result.total_outgoing_power_w_m2,
        result.total_incoming_power_w_m2,
        2e-12 * result.total_incoming_power_w_m2,
    );
    check_tensor_fields(&prepared, &result);
    let sh = prepared
        .velocity_coefficients(ElasticSide::First, ElasticMode::Sh)
        .unwrap();
    close(sh[0][2].norm_sq(), 1.0, 2e-14);
    assert!(sh[1][2].abs() > 0.0);
    assert_eq!(sh[0][0], C64::ZERO);
    assert_eq!(sh[0][1], C64::ZERO);
}

#[test]
fn g1_critical_branch_retains_the_original_speed_pair_without_angle_snapping() {
    // Exact raw speed values decide the sign, including an almost matched
    // grazing channel and the two adjacent floats around a critical angle.
    assert_eq!(normal_ratio(1.0, 1.0, 1e-15).unwrap(), C64::from_re(1e-15));
    assert!(normal_ratio(1.0, 1.0_f64.next_up(), 1e-9).unwrap().im > 0.0);
    assert!(normal_ratio(1.0, 1.0_f64.next_down(), 1e-9).unwrap().re > 0.0);
    let below = normal_ratio(3.0, 5.0, 0.8_f64.next_down()).unwrap();
    let above = normal_ratio(3.0, 5.0, 0.8_f64.next_up()).unwrap();
    assert!(below.im > 0.0 && above.re > 0.0);
    let first = material(3.0, 3000.0, 1400.0, 9);
    let second = material(5.0, 5000.0, 2300.0, 10);
    let owner = IsotropicElasticInterface::new(&first, &second).unwrap();
    let critical = (1.0 - (owner.speeds[0][0] / owner.speeds[1][0]).powi(2)).sqrt();
    let mut columns = Vec::new();
    for cosine in [critical - 1e-10, critical + 1e-10] {
        let p = owner
            .prepare_incidence(ElasticSide::First, ElasticMode::P, cosine)
            .unwrap();
        columns.push(
            p.velocity_coefficients(ElasticSide::First, ElasticMode::P)
                .unwrap(),
        );
        let response = p
            .drive([[C64::from_re(0.01), C64::ZERO, C64::ZERO], [C64::ZERO; 3]])
            .unwrap();
        check_tensor_fields(&p, &response);
        if cosine < critical {
            assert_eq!(response.outgoing_normal_power_w_m2[1][0], 0.0);
            assert!(response.outgoing_normal_power_w_m2[1][1] > 0.0);
        }
    }
    for s in 0..2 {
        for m in 0..3 {
            complex_close(columns[0][s][m], columns[1][s][m], 5e-4);
        }
    }
}

#[test]
fn g1_grazing_zero_lame_source_retains_its_small_normal_traction() {
    let first =
        ElasticAcousticState::from_bulk_and_shear(6000.0, 2e9, 3e9, ContentHash([22; 32]), 0.01)
            .unwrap();
    assert_eq!(first.stiffness_mandel_pa()[0][1], 0.0);
    let owner = IsotropicElasticInterface::new(&first, &first).unwrap();
    let cp = owner.phase_speeds_m_s(ElasticSide::First)[0];
    for cosine in [1e-6, 1e-10, 1e-14] {
        let p = owner
            .prepare_incidence(ElasticSide::First, ElasticMode::P, cosine)
            .unwrap();
        let wave = p.outgoing_mode(ElasticSide::Second, ElasticMode::P);
        // With lambda=0, sigma_zz = 2*mu*epsilon_zz. No subtraction of
        // nearly equal speed-ratio terms may set the normal-traction scale.
        let expected = -(6e9 / cp) * cosine * cosine;
        close(wave.traction_pa_s_m[2].re, expected, expected.abs() * 2e-14);
        assert_eq!(wave.traction_pa_s_m[2].im, 0.0);
    }
}

#[test]
fn g0_flux_normalized_operator_is_unitary_and_reciprocal() {
    let first = material(2300.0, 4500.0, 2500.0, 11);
    let second = material(2900.0, 6100.0, 3100.0, 12);
    let owner = IsotropicElasticInterface::new(&first, &second).unwrap();
    for cosine in [1.0, 0.95, 0.85, 0.7, 0.4] {
        let p = owner
            .prepare_incidence(ElasticSide::First, ElasticMode::P, cosine)
            .unwrap();
        let weight: [f64; 6] = core::array::from_fn(|i| p.incoming[i].normal_power_factor_pa_s_m);
        for i in 0..6 {
            if weight[i] == 0.0 {
                continue;
            }
            for j in 0..6 {
                if weight[j] == 0.0 {
                    continue;
                }
                let inner = (0..6).fold(C64::ZERO, |a, k| {
                    a + (p.scattering[k][i].conj() * p.scattering[k][j]).scale(weight[k])
                });
                complex_close(
                    inner,
                    C64::from_re(if i == j { weight[i] } else { 0.0 }),
                    3e-12 * (weight[i] * weight[j]).sqrt(),
                );
                // Polarization phase conventions can change a reciprocal
                // coefficient's sign; reciprocal modal NORMAL power is equal.
                close(
                    p.scattering[j][i].norm_sq() * weight[j] / weight[i],
                    p.scattering[i][j].norm_sq() * weight[i] / weight[j],
                    3e-12,
                );
            }
        }
    }
}

#[test]
fn g3_common_material_scaling_preserves_conversion_and_scales_actual_work() {
    let first = material(1700.0, 4100.0, 2300.0, 13);
    let second = material(3800.0, 5500.0, 2900.0, 14);
    let scaled_first = material(17000.0, 4100.0, 2300.0, 15);
    let scaled_second = material(38000.0, 5500.0, 2900.0, 16);
    let p = IsotropicElasticInterface::new(&first, &second)
        .unwrap()
        .prepare_incidence(ElasticSide::First, ElasticMode::P, 0.9)
        .unwrap();
    let q = IsotropicElasticInterface::new(&scaled_first, &scaled_second)
        .unwrap()
        .prepare_incidence(ElasticSide::First, ElasticMode::P, 0.9)
        .unwrap();
    let sources = [
        [
            C64::new(0.001, -0.002),
            C64::new(-0.003, 0.001),
            C64::from_re(0.001),
        ],
        [
            C64::from_re(-0.002),
            C64::new(0.0003, 0.002),
            C64::from_re(0.002),
        ],
    ];
    let a = p.drive(sources).unwrap();
    let b = q.drive(sources).unwrap();
    for s in 0..2 {
        for m in 0..3 {
            complex_close(
                a.outgoing_velocity_m_s[s][m],
                b.outgoing_velocity_m_s[s][m],
                1e-16,
            );
        }
    }
    close(
        b.total_incoming_power_w_m2,
        10.0 * a.total_incoming_power_w_m2,
        1e-11 * a.total_incoming_power_w_m2,
    );
    close(
        b.total_outgoing_power_w_m2,
        10.0 * a.total_outgoing_power_w_m2,
        1e-11 * a.total_incoming_power_w_m2,
    );
    assert_eq!(
        b.material_state_identities,
        [ContentHash([15; 32]), ContentHash([16; 32])]
    );
}

#[test]
fn g1_auxetic_and_high_contrast_solids_retain_coupled_shear_response() {
    let first =
        ElasticAcousticState::from_bulk_and_shear(1000.0, 1e6, 5e9, ContentHash([17; 32]), 0.01)
            .unwrap();
    let second = material(9000.0, 8000.0, 3900.0, 18);
    assert!(first.stiffness_mandel_pa()[0][1] < 0.0);
    let interface = IsotropicElasticInterface::new(&first, &second).unwrap();
    for cosine in [1.0, 0.9, 0.6, 0.02] {
        for kind in [ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh] {
            let p = interface
                .prepare_incidence(ElasticSide::First, kind, cosine)
                .unwrap();
            let mut input = [[C64::ZERO; 3]; 2];
            input[0][kind.index()] = C64::new(0.001, 0.001);
            let r = p.drive(input).unwrap();
            check_tensor_fields(&p, &r);
            close(
                r.total_incoming_power_w_m2,
                r.total_outgoing_power_w_m2,
                2e-11 * r.total_incoming_power_w_m2,
            );
        }
    }
}

#[test]
fn g0_g4_invalid_excitation_original_domain_and_cancellation_refuse() {
    let first = material(2500.0, 4000.0, 2000.0, 19);
    let second = material(2700.0, 6000.0, 3200.0, 20);
    let arbitrary = ElasticAcousticState::try_new_mandel(
        *first.stiffness_mandel_pa(),
        first.density_kg_m3(),
        ContentHash([21; 32]),
        0.02,
    )
    .unwrap();
    assert!(IsotropicElasticInterface::new(&arbitrary, &second).is_err());
    let owner = IsotropicElasticInterface::new(&first, &second).unwrap();
    for cosine in [0.0, -0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(
            owner
                .prepare_incidence(ElasticSide::First, ElasticMode::P, cosine)
                .is_err()
        );
    }
    for stop in [1, 5, 12, 25] {
        let mut count = 0;
        let result = owner.prepare_incidence_cancellable(
            ElasticSide::First,
            ElasticMode::P,
            0.9,
            &mut || {
                count += 1;
                count == stop
            },
        );
        assert!(result.is_err(), "cancel poll {stop}");
    }
    let prepared = owner
        .prepare_incidence(ElasticSide::First, ElasticMode::Sv, 0.5)
        .unwrap();
    let mut input = [[C64::ZERO; 3]; 2];
    input[1][0] = C64::ONE;
    assert!(prepared.drive(input).is_err());
    input = [[C64::ZERO; 3]; 2];
    input[0][1] = C64::new(f64::NAN, 0.0);
    assert!(prepared.drive(input).is_err());
    input[0][1] = C64::from_re(10000.0);
    assert!(prepared.drive(input).is_err());
    let silent = prepared.drive([[C64::ZERO; 3]; 2]).unwrap();
    assert_eq!(silent.total_incoming_power_w_m2, 0.0);
    assert_eq!(silent.total_outgoing_power_w_m2, 0.0);
    assert_eq!(silent.boundary_normal_flux_w_m2, [0.0; 2]);
}
