use super::*;

fn material(rho: f64, cp: f64, cs: f64, id: u8) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        rho,
        rho * (cp * cp - 4.0 / 3.0 * cs * cs),
        rho * cs * cs,
        ContentHash([id; 32]),
        0.02,
    )
    .unwrap()
}

fn mode(i: usize) -> ElasticMode {
    [ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh][i]
}
fn side(i: usize) -> ElasticSide {
    if i == 0 {
        ElasticSide::First
    } else {
        ElasticSide::Second
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:e}"
    );
}
fn cc(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:?} != {b:?}; tolerance={tolerance:e}"
    );
}
fn phase(x: f64) -> C64 {
    C64::new(x.cos(), x.sin())
}
fn single(side: usize, kind: usize, amplitude: C64) -> [[C64; 3]; 2] {
    let mut a = [[C64::ZERO; 3]; 2];
    a[side][kind] = amplitude;
    a
}

fn check_fields(prepared: &PreparedElasticLayer<'_>, response: &ElasticLayerResponse<'_, '_>) {
    let power = response.exterior.total_incoming_power_w_m2;
    let h = prepared.thickness_m();
    for fraction in [0.0, 0.17, 0.4, 0.5, 0.77, 1.0] {
        let field = response.field_at(h * fraction).unwrap();
        // Independently contract the ORIGINAL Mandel tensor with the returned
        // interior strain, including all coherent mode-conversion cross terms.
        let stress: [C64; 6] = core::array::from_fn(|i| {
            (0..6).fold(C64::ZERO, |a, j| {
                a + field.strain_mandel[j].scale(prepared.layer.stiffness_mandel_pa()[i][j])
            })
        });
        let traction = [
            stress[5].scale(core::f64::consts::FRAC_1_SQRT_2),
            stress[4].scale(core::f64::consts::FRAC_1_SQRT_2),
            stress[2],
        ];
        for i in 0..3 {
            cc(
                field.traction_pa[i],
                traction[i],
                3e-10 * traction[i].abs().max(1.0),
            );
        }
        let flux = -0.5
            * (0..3)
                .map(|i| (field.velocity_m_s[i].conj() * traction[i]).re)
                .sum::<f64>();
        close(field.normal_flux_w_m2, flux, 1e-10 * power.max(1e-30));
        close(
            flux,
            response.exterior.boundary_normal_flux_w_m2[0],
            2e-10 * power.max(1e-30),
        );
        let norm = field
            .strain_mandel
            .iter()
            .map(|x| x.norm_sq())
            .sum::<f64>()
            .sqrt();
        assert!(norm <= response.layer_strain_norm_bound * (1.0 + 8e-15));
    }
    for s in 0..2 {
        let field = response.layer_faces.unwrap()[s];
        for c in 0..3 {
            cc(
                field.velocity_m_s[c],
                response.exterior.boundary_velocity_m_s[s][c],
                2e-12,
            );
            cc(
                field.traction_pa[c],
                response.exterior.boundary_traction_pa[s][c],
                2e-10 * field.traction_pa[c].abs().max(1.0),
            );
        }
    }
}

#[test]
fn normal_slab_matches_independent_fabry_perot_for_every_mode_and_side() {
    let first = material(2400.0, 4200.0, 2100.0, 1);
    let layer = material(1900.0, 4800.0, 2600.0, 2);
    let second = material(3100.0, 6100.0, 3300.0, 3);
    let f = 1300.0;
    for h in [1e-9, 0.013, 0.3, 1.7, 4.8] {
        let owner = IsotropicElasticLayer::new(&first, &layer, &second, h).unwrap();
        let p = owner
            .prepare_incidence(ElasticSide::First, ElasticMode::P, 1.0, f)
            .unwrap();
        for s in 0..2 {
            for m in 0..3 {
                let z0 = owner.exterior.impedances[s][m.min(1)];
                let z1 = owner.impedances[m.min(1)];
                let z2 = owner.exterior.impedances[1 - s][m.min(1)];
                let sign = if m == 0 { 1.0 } else { -1.0 };
                let r01 = sign * (z1 - z0) / (z1 + z0);
                let r12 = sign * (z2 - z1) / (z2 + z1);
                let e = phase(core::f64::consts::TAU * f * h / owner.speeds[m.min(1)]);
                let denominator = C64::ONE + (e * e).scale(r01 * r12);
                let r = (C64::from_re(r01) + (e * e).scale(r12)) / denominator;
                let t = e.scale((2.0 * z0 / (z0 + z1)) * (2.0 * z1 / (z1 + z2))) / denominator;
                let a = p.velocity_coefficients(side(s), mode(m)).unwrap();
                cc(a[s][m], r, 3e-13);
                cc(a[1 - s][m], t, 3e-13);
                for j in 0..3 {
                    if j != m {
                        assert_eq!(a[0][j], C64::ZERO);
                        assert_eq!(a[1][j], C64::ZERO);
                    }
                }
                let response = p.drive(single(s, m, C64::new(0.001, -0.002))).unwrap();
                check_fields(&p, &response);
            }
        }
    }
}

#[test]
fn zero_thickness_and_original_homogeneous_material_limits_are_exact() {
    let a = material(2600.0, 6000.0, 3000.0, 4);
    let b = material(1900.0, 4000.0, 2200.0, 5);
    let c = material(2600.0, 6000.0, 3000.0, 6);
    let interface = IsotropicElasticInterface::new(&a, &b).unwrap();
    let absent = IsotropicElasticLayer::new(&a, &c, &b, 0.0).unwrap();
    for cosine in [1.0, 0.8, 0.1] {
        let original = interface
            .prepare_incidence(ElasticSide::First, ElasticMode::P, cosine)
            .unwrap();
        let p = absent
            .prepare_incidence(ElasticSide::First, ElasticMode::P, cosine, 1000.0)
            .unwrap();
        for s in 0..2 {
            for m in 0..3 {
                assert_eq!(
                    p.velocity_coefficients(side(s), mode(m)),
                    original.velocity_coefficients(side(s), mode(m))
                );
            }
        }
        let r = p.drive(single(0, 0, C64::from_re(0.001))).unwrap();
        assert!(r.layer_faces.is_none());
        assert!(r.field_at(0.0).is_err());
    }
    for h in [0.001, 2.0] {
        let owner = IsotropicElasticLayer::new(&a, &c, &a, h).unwrap();
        for cosine in [1.0, 0.8, 1e-8] {
            let p = owner
                .prepare_incidence(ElasticSide::First, ElasticMode::Sv, cosine, 800.0)
                .unwrap();
            for s in 0..2 {
                for m in 0..3 {
                    if let Some(out) = p.velocity_coefficients(side(s), mode(m)) {
                        assert_eq!(out[s], [C64::ZERO; 3]);
                        let q = p
                            .incoming_mode(side(s), mode(m))
                            .unwrap()
                            .normal_slowness_s_m
                            .re
                            .abs();
                        cc(
                            out[1 - s][m],
                            phase(core::f64::consts::TAU * 800.0 * q * h),
                            2e-14,
                        );
                    }
                }
            }
            let r = p.drive(single(0, 1, C64::new(0.001, 0.002))).unwrap();
            check_fields(&p, &r);
        }
    }
}

#[test]
fn oblique_coherent_conversion_closes_work_and_interior_displacement_gradients() {
    let a = material(2300.0, 4300.0, 2300.0, 7);
    let b = material(3500.0, 6200.0, 3300.0, 8);
    let c = material(2700.0, 4900.0, 2600.0, 9);
    for h in [0.005, 0.8] {
        let p = IsotropicElasticLayer::new(&a, &b, &c, h)
            .unwrap()
            .prepare_incidence(ElasticSide::First, ElasticMode::P, 0.89, 1900.0)
            .unwrap();
        let input = core::array::from_fn(|s| {
            core::array::from_fn(|m| C64::new((m + 1) as f64 * 0.001, (s + m + 1) as f64 * -0.0004))
        });
        let result = p.drive(input).unwrap();
        check_fields(&p, &result);
        assert!(
            p.velocity_coefficients(ElasticSide::First, ElasticMode::P)
                .unwrap()[0][1]
                .abs()
                > 0.01
        );
        // Four-point differentiation of actual interior velocities reconstructs
        // strain independently of the production modal strain traces.
        let z = 0.43 * h;
        let dz = h * 1e-4;
        let v = [-2.0, -1.0, 1.0, 2.0].map(|d| result.field_at(z + d * dz).unwrap().velocity_m_s);
        let derivative: [C64; 3] = core::array::from_fn(|i| {
            (v[0][i] - v[1][i].scale(8.0) + v[2][i].scale(8.0) - v[3][i]).scale(1.0 / (12.0 * dz))
        });
        let f = result.field_at(z).unwrap();
        let q = p.tangential_slowness_s_m();
        let iw = C64::new(0.0, 1.0 / (core::f64::consts::TAU * 1900.0));
        let e = [
            f.velocity_m_s[0].scale(-q),
            C64::ZERO,
            iw * derivative[2],
            f.velocity_m_s[1].scale(-q * core::f64::consts::FRAC_1_SQRT_2),
            (iw * derivative[1]).scale(core::f64::consts::FRAC_1_SQRT_2),
            (iw * derivative[0] - f.velocity_m_s[2].scale(q))
                .scale(core::f64::consts::FRAC_1_SQRT_2),
        ];
        for i in 0..6 {
            cc(
                f.strain_mandel[i],
                e[i],
                3e-9 * result.layer_strain_norm_bound,
            );
        }
        let zero = p.drive([[C64::ZERO; 3]; 2]).unwrap();
        assert_eq!(zero.field_at(0.3 * h).unwrap().velocity_m_s, [C64::ZERO; 3]);
    }
}

#[test]
fn finite_layer_scattering_is_flux_unitary_and_reciprocal() {
    let a = material(2300.0, 4500.0, 2500.0, 10);
    let b = material(3800.0, 7000.0, 4000.0, 11);
    let c = material(2900.0, 6100.0, 3100.0, 12);
    for cosine in [1.0, 0.95, 0.8, 0.6, 0.2] {
        let p = IsotropicElasticLayer::new(&a, &b, &c, 0.71)
            .unwrap()
            .prepare_incidence(ElasticSide::First, ElasticMode::P, cosine, 2400.0)
            .unwrap();
        let e = &p.exterior;
        let w = e.incoming.map(|trace| trace.normal_power_factor_pa_s_m);
        for i in 0..6 {
            if w[i] > 0.0 {
                for j in 0..6 {
                    if w[j] > 0.0 {
                        let inner = (0..6).fold(C64::ZERO, |s, k| {
                            s + (e.scattering[k][i].conj() * e.scattering[k][j]).scale(w[k])
                        });
                        cc(
                            inner,
                            C64::from_re(if i == j { w[i] } else { 0.0 }),
                            3e-11 * (w[i] * w[j]).sqrt(),
                        );
                        close(
                            e.scattering[j][i].norm_sq() * w[j] / w[i],
                            e.scattering[i][j].norm_sq() * w[i] / w[j],
                            3e-11,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn evanescent_layer_transmits_interference_power_and_reaches_opaque_limit() {
    let outer = material(2500.0, 4000.0, 2000.0, 13);
    let layer = material(3000.0, 8000.0, 4500.0, 14);
    let mut transmitted = Vec::new();
    for h in [0.001, 0.08, 0.8, 40.0] {
        let p = IsotropicElasticLayer::new(&outer, &layer, &outer, h)
            .unwrap()
            .prepare_incidence(ElasticSide::First, ElasticMode::Sv, 0.7, 1000.0)
            .unwrap();
        assert!(
            p.layer_normal_slownesses_s_m()
                .unwrap()
                .iter()
                .all(|q| q.re == 0.0 && q.im > 0.0)
        );
        let r = p.drive(single(0, 2, C64::from_re(0.002))).unwrap();
        check_fields(&p, &r);
        let t = r.exterior.outgoing_normal_power_w_m2[1][2] / r.exterior.total_incoming_power_w_m2;
        let z0 = p
            .incoming_mode(ElasticSide::First, ElasticMode::Sh)
            .unwrap()
            .normal_slowness_s_m
            .scale(0.5 * outer.stiffness_mandel_pa()[4][4]);
        let z1 = p.layer_normal_slownesses_s_m().unwrap()[2]
            .scale(0.5 * layer.stiffness_mandel_pa()[4][4]);
        let reflection = (z0 - z1) / (z0 + z1);
        let decay =
            (-core::f64::consts::TAU * 1000.0 * p.layer_normal_slownesses_s_m().unwrap()[2].im * h)
                .exp();
        let analytic = ((C64::ONE - reflection * reflection).scale(decay))
            / (C64::ONE - (reflection * reflection).scale(decay * decay));
        cc(
            p.velocity_coefficients(ElasticSide::First, ElasticMode::Sh)
                .unwrap()[1][2],
            analytic,
            2e-12 * analytic.abs(),
        );
        assert!(t > 0.0 && t < 1.0);
        assert!(r.field_at(0.5 * h).unwrap().normal_flux_w_m2 > 0.0);
        transmitted.push(t);
    }
    assert!(transmitted.windows(2).all(|w| w[0] > w[1]));
    assert!(transmitted[0] > 0.99);
    assert!(transmitted[3] < 1e-50, "transmitted={transmitted:?}");
}

#[test]
fn critical_internal_shear_has_the_finite_linear_field_limit() {
    let outer = material(2500.0, 5000.0, 3000.0, 15);
    let layer = material(3200.0, 9000.0, 5000.0, 16);
    let h = 0.1;
    let f = 1700.0;
    let owner = IsotropicElasticLayer::new(&outer, &layer, &outer, h).unwrap();
    let critical = (1.0 - (owner.exterior.speeds[0][1] / owner.speeds[1]).powi(2)).sqrt();
    let mut prior = None;
    for cosine in [critical.next_down(), critical, critical.next_up()] {
        let p = owner
            .prepare_incidence(ElasticSide::First, ElasticMode::Sh, cosine, f)
            .unwrap();
        assert!(p.basis.unwrap()[2].regularized);
        let r = p.drive(single(0, 2, C64::from_re(0.003))).unwrap();
        check_fields(&p, &r);
        let q = p
            .incoming_mode(ElasticSide::First, ElasticMode::Sh)
            .unwrap()
            .normal_slowness_s_m
            .re;
        let outer_mu = outer.stiffness_mandel_pa()[4][4] * 0.5;
        let inner_mu = layer.stiffness_mandel_pa()[4][4] * 0.5;
        let a = core::f64::consts::TAU * f * h * outer_mu * q / inner_mu;
        let t_limit = C64::from_re(2.0) / C64::new(2.0, -a);
        let t = p
            .velocity_coefficients(ElasticSide::First, ElasticMode::Sh)
            .unwrap()[1][2];
        cc(t, t_limit, 2e-13);
        if let Some(old) = prior {
            cc(t, old, 2e-13);
        }
        prior = Some(t);
        // The derivative basis at exact q=0 has a constant term and a linear
        // term, not two copies of a tangential plane wave.
        let mut exact = p.basis.unwrap()[2];
        exact.cosine = C64::ZERO;
        let left = exact.at(0.0, h);
        let mid = exact.at(0.5 * h, h);
        let right = exact.at(h, h);
        for parity in 0..2 {
            for component in 0..12 {
                cc(
                    mid[parity].0[component],
                    (left[parity].0[component] + right[parity].0[component]).scale(0.5),
                    1e-9,
                );
            }
        }
        assert!((right[1].0[1] - left[1].0[1]).abs() > 0.0);
    }
}

#[test]
fn interior_antinode_cannot_escape_the_original_layer_strain_ceiling() {
    let rho = 4096.0;
    let cp = 4096.0;
    let cs = 2048.0;
    let outer = material(rho, cp, cs, 17);
    let loose = material(rho, cp, cs, 18);
    let limited = ElasticAcousticState::from_bulk_and_shear(
        rho,
        rho * (cp * cp - 4.0 / 3.0 * cs * cs),
        rho * cs * cs,
        ContentHash([19; 32]),
        1e-6,
    )
    .unwrap();
    let f = 1024.0;
    let h = cp / (2.0 * f);
    let input = [[C64::from_re(0.003), C64::ZERO, C64::ZERO]; 2];
    let p = IsotropicElasticLayer::new(&outer, &loose, &outer, h)
        .unwrap()
        .prepare_incidence(ElasticSide::First, ElasticMode::P, 1.0, f)
        .unwrap();
    let r = p.drive(input).unwrap();
    assert!(
        r.layer_faces
            .unwrap()
            .iter()
            .all(|face| face.strain_mandel.iter().all(|e| e.abs() < 1e-18))
    );
    assert!(r.field_at(h / 2.0).unwrap().strain_mandel[2].abs() > 1e-6);
    let q = IsotropicElasticLayer::new(&outer, &limited, &outer, h)
        .unwrap()
        .prepare_incidence(ElasticSide::First, ElasticMode::P, 1.0, f)
        .unwrap();
    assert!(q.drive(input).is_err());
    assert!(q.drive([[C64::ZERO; 3]; 2]).is_ok());
}

#[test]
fn layer_refuses_invalid_phases_sources_depths_and_cancellation() {
    let a = material(2500.0, 4000.0, 2000.0, 20);
    let b = material(3000.0, 7000.0, 3900.0, 21);
    for h in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(IsotropicElasticLayer::new(&a, &b, &a, h).is_err());
    }
    let owner = IsotropicElasticLayer::new(&a, &b, &a, 0.2).unwrap();
    for f in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            owner
                .prepare_incidence(ElasticSide::First, ElasticMode::P, 0.9, f)
                .is_err()
        );
    }
    for cutoff in [1, 4, 9, 20] {
        let mut polls = 0;
        assert!(
            owner
                .prepare_incidence_cancellable(
                    ElasticSide::First,
                    ElasticMode::P,
                    0.9,
                    1000.0,
                    &mut || {
                        polls += 1;
                        polls >= cutoff
                    }
                )
                .is_err()
        );
    }
    let p = owner
        .prepare_incidence(ElasticSide::First, ElasticMode::Sv, 0.5, 1000.0)
        .unwrap();
    assert!(
        p.incoming_mode(ElasticSide::First, ElasticMode::P)
            .is_none()
    );
    assert!(p.drive(single(0, 0, C64::ONE)).is_err());
    assert!(p.drive(single(0, 1, C64::from_re(f64::NAN))).is_err());
    let r = p.drive(single(0, 1, C64::from_re(0.001))).unwrap();
    for z in [-1e-30, 0.2_f64.next_up(), f64::NAN] {
        assert!(r.field_at(z).is_err());
    }
    assert!(r.field_at(0.1).is_ok());
}
