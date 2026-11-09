use super::*;
use crate::elastic::acoustic::interface::surface::ElasticSurfaceImpedance;

fn material(rho: f64, cp: f64, cs: f64, id: u8, limit: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        rho,
        rho * (cp * cp - 4.0 / 3.0 * cs * cs),
        rho * cs * cs,
        ContentHash([id; 32]),
        limit,
    )
    .unwrap()
}
fn solid() -> ElasticAcousticState {
    material(2500.0, 4000.0, 2000.0, 31, 0.02)
}
fn scalar(normal: C64, shear: C64) -> ElasticSurfaceCondition {
    ElasticSurfaceCondition::Impedance(ElasticSurfaceImpedance::new(normal, shear).unwrap())
}
fn cc(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:?} != {b:?}");
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn phase(q: C64, omega_h: f64) -> C64 {
    C64::new((q.re * omega_h).cos(), (q.re * omega_h).sin()).scale((-q.im * omega_h).exp())
}
fn tensor() -> ElasticSurfaceCondition {
    let matrix = [
        [
            C64::new(1.5e6, -2e7),
            C64::new(-1e5, -4e6),
            C64::new(-1.5e5, 4e6),
        ],
        [
            C64::new(-1e5, -4e6),
            C64::new(2e6, -3e7),
            C64::new(-4e5, 3.5e7),
        ],
        [
            C64::new(-1.5e5, 4e6),
            C64::new(-4e5, 3.5e7),
            C64::new(3e6, -9e7),
        ],
    ];
    ElasticSurfaceCondition::TensorImpedance(ElasticSurfaceImpedanceTensor::new(matrix).unwrap())
}

fn check_fields(response: &ElasticCoatingResponse<'_, '_>) {
    let h = response.prepared.thickness_m;
    let power = response.reflection.total_incoming_power_w_m2;
    for fraction in [0.0, 0.23, 0.5, 0.81, 1.0] {
        let field = response.field_at(h * fraction).unwrap();
        let stress: [C64; 6] = core::array::from_fn(|i| {
            (0..6).fold(C64::ZERO, |s, j| {
                s + field.strain_mandel[j]
                    .scale(response.prepared.coating.stiffness_mandel_pa()[i][j])
            })
        });
        let traction = [
            stress[5].scale(core::f64::consts::FRAC_1_SQRT_2),
            stress[4].scale(core::f64::consts::FRAC_1_SQRT_2),
            stress[2],
        ];
        for i in 0..3 {
            cc(field.traction_pa[i], traction[i], 2e-8);
        }
        close(
            field.normal_flux_w_m2,
            response.reflection.absorbed_normal_power_w_m2,
            3e-10 * power,
        );
        let strain = field
            .strain_mandel
            .iter()
            .map(|s| s.norm_sq())
            .sum::<f64>()
            .sqrt();
        assert!(strain <= response.coating_strain_norm_bound * (1.0 + 2e-14));
    }
    let faces = response.coating_faces.unwrap();
    for i in 0..3 {
        cc(
            faces[0].velocity_m_s[i],
            response.reflection.boundary_velocity_m_s[i],
            2e-15,
        );
        cc(
            faces[0].traction_pa[i],
            response.reflection.boundary_traction_pa[i],
            3e-8,
        );
        cc(
            faces[1].velocity_m_s[i],
            response.support_velocity_m_s[i],
            2e-15,
        );
        cc(
            faces[1].traction_pa[i],
            response.support_traction_pa[i],
            3e-8,
        );
    }
}

#[test]
fn absent_coating_is_the_exact_original_surface_for_every_backing() {
    let solid = solid();
    let coating = material(1800.0, 6000.0, 3000.0, 32, 0.02);
    for backing in [
        ElasticSurfaceCondition::TractionFree,
        ElasticSurfaceCondition::Clamped,
        scalar(C64::new(2e6, 3e6), C64::new(4e6, -1e6)),
        tensor(),
    ] {
        let original = IsotropicElasticSurface::new(&solid, backing)
            .unwrap()
            .prepare_incidence(ElasticMode::Sv, 0.5)
            .unwrap();
        let prepared = IsotropicElasticCoating::new(&solid, &coating, 0.0, backing)
            .unwrap()
            .prepare_incidence(ElasticMode::Sv, 0.5, 1000.0)
            .unwrap();
        let input = [C64::ZERO, C64::new(0.001, -0.0004), C64::from_re(0.0006)];
        let a = original.drive(input).unwrap();
        let b = prepared.drive(input).unwrap();
        assert_eq!(
            a.reflected_velocity_m_s,
            b.reflection.reflected_velocity_m_s
        );
        assert_eq!(
            a.absorbed_normal_power_w_m2,
            b.reflection.absorbed_normal_power_w_m2
        );
        assert_eq!(a.boundary_velocity_m_s, b.support_velocity_m_s);
        assert_eq!(a.boundary_traction_pa, b.support_traction_pa);
        assert!(b.coating_faces.is_none());
        assert!(b.field_at(0.0).is_err());
        assert_eq!(
            prepared.material_state_identities(),
            [ContentHash([31; 32]), ContentHash([32; 32])]
        );
    }
}

#[test]
fn normal_coating_matches_independent_impedance_transform_and_quarter_wave_matching() {
    let solid = solid();
    let coating = material(1800.0, 3000.0, 1500.0, 32, 0.02);
    let omega = 1700.0;
    let hz = omega / core::f64::consts::TAU;
    let z = [C64::new(2e6, 3e6), C64::new(4e6, -1e6), C64::new(4e6, -1e6)];
    for h in [1e-9, 0.013, 0.7, 2.0] {
        let prepared = IsotropicElasticCoating::new(&solid, &coating, h, scalar(z[0], z[1]))
            .unwrap()
            .prepare_incidence(ElasticMode::P, 1.0, hz)
            .unwrap();
        for (m, mode) in [ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh]
            .into_iter()
            .enumerate()
        {
            let (d0, d1, speed) = if m == 0 {
                (1e7, 5.4e6, 3000.0)
            } else {
                (5e6, 2.7e6, 1500.0)
            };
            let theta: f64 = omega * h / speed;
            let a = C64::from_re(theta.cos()) - C64::new(0.0, theta.sin() / d1) * z[m];
            let b = C64::new(0.0, d1 * theta.sin()) - z[m].scale(theta.cos());
            let rear = C64::from_re(2.0 * d0 * 0.001) / (a.scale(d0) - b);
            let expected = (a * rear - C64::from_re(0.001)).scale(if m == 0 { -1.0 } else { 1.0 });
            let mut input = [C64::ZERO; 3];
            input[m] = C64::from_re(0.001);
            let response = prepared.drive(input).unwrap();
            cc(
                response.reflection.reflected_velocity_m_s[m],
                expected,
                3e-15,
            );
            cc(
                response.support_velocity_m_s[if m == 0 { 2 } else { m - 1 }],
                rear,
                3e-15,
            );
            close(
                response.reflection.absorbed_normal_power_w_m2,
                0.5 * z[m].re * rear.norm_sq(),
                2e-11,
            );
            assert!(prepared.velocity_coefficients(mode).is_some());
            check_fields(&response);
        }
    }
    // Independent quarter-wave transformer: Zin=Zcoating²/Zback=Zincident.
    let h = core::f64::consts::PI * 3000.0 / (2.0 * omega);
    let prepared = IsotropicElasticCoating::new(
        &solid,
        &coating,
        h,
        scalar(C64::from_re(2_916_000.0), C64::from_re(2.7e6)),
    )
    .unwrap()
    .prepare_incidence(ElasticMode::P, 1.0, hz)
    .unwrap();
    let response = prepared
        .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
        .unwrap();
    cc(
        response.reflection.reflected_velocity_m_s[0],
        C64::ZERO,
        3e-15,
    );
    cc(
        response.support_velocity_m_s[2],
        C64::new(0.0, 0.00185185185185185185),
        3e-15,
    );
    close(response.reflection.absorbed_normal_power_w_m2, 5.0, 2e-11);
    check_fields(&response);
}

#[test]
fn homogeneous_tensor_coating_retains_each_roundtrip_phase_and_support_motion() {
    let state = solid();
    let omega = 1700.0;
    let h = 0.7;
    let surface = IsotropicElasticSurface::new(&state, tensor())
        .unwrap()
        .prepare_incidence(ElasticMode::Sv, 0.5)
        .unwrap();
    let prepared = IsotropicElasticCoating::new(&state, &state, h, tensor())
        .unwrap()
        .prepare_incidence(ElasticMode::Sv, 0.5, omega / core::f64::consts::TAU)
        .unwrap();
    let modes = [ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh];
    let input = [
        C64::ZERO,
        C64::new(0.001, 0.0002),
        C64::new(-0.0003, 0.0006),
    ];
    let rear_input = core::array::from_fn(|i| {
        surface
            .incoming_mode(modes[i])
            .map(|w| phase(w.normal_slowness_s_m, omega * h) * input[i])
            .unwrap_or(C64::ZERO)
    });
    let direct = surface.drive(rear_input).unwrap();
    let response = prepared.drive(input).unwrap();
    for i in 0..3 {
        let expected = phase(
            -surface.reflected_mode(modes[i]).normal_slowness_s_m,
            omega * h,
        ) * direct.reflected_velocity_m_s[i];
        cc(
            response.reflection.reflected_velocity_m_s[i],
            expected,
            3e-15,
        );
        cc(
            response.support_velocity_m_s[i],
            direct.boundary_velocity_m_s[i],
            3e-15,
        );
    }
    close(
        response.reflection.absorbed_normal_power_w_m2,
        direct.absorbed_normal_power_w_m2,
        2e-11,
    );
    check_fields(&response);
}

#[test]
fn oblique_tensor_at_internal_p_criticality_matches_independent_continuum_exponential() {
    let state = solid();
    let coating = material(1800.0, 6000.0, 3000.0, 32, 0.02);
    // Independent 70-digit state ODE exponential y'=-i omega A y, not the
    // production mode basis. Cartesian components are vxyz,txyz at both faces.
    let xx = C64::new(1565245.37122186327, -21770499.5486048268);
    let xy = C64::new(-110136.843637422769, -4200573.72624464334);
    let xz = C64::new(-146927.215553463142, 3933299.11822961398);
    let yy = C64::new(2139620.11087995800, -33587504.03581203597);
    let yz = C64::new(-413412.792404436721, 36346344.46456713227);
    let zz = C64::new(2795134.51789817873, -93759643.47440666665);
    let backing = ElasticSurfaceCondition::TensorImpedance(
        ElasticSurfaceImpedanceTensor::new([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]).unwrap(),
    );
    let prepared = IsotropicElasticCoating::new(&state, &coating, 0.4, backing)
        .unwrap()
        .prepare_incidence(
            ElasticMode::P,
            5.0_f64.sqrt() / 3.0,
            1700.0 / core::f64::consts::TAU,
        )
        .unwrap();
    let response = prepared
        .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
        .unwrap();
    let reflected = [
        C64::new(0.000477320409248758921, -0.000088341450532083812),
        C64::new(-0.001036714836951965164, 0.000088050392703874206),
        C64::new(0.000059365158627390571, 0.000233698939208309649),
    ];
    let rear = [
        C64::new(0.0000494309540591284862, 0.000239899380385387949),
        C64::new(0.000104171475547877968, 0.000217539282860588838),
        C64::new(0.0000644739463404441506, 0.000247418044554392257),
    ];
    for i in 0..3 {
        cc(
            response.reflection.reflected_velocity_m_s[i],
            reflected[i],
            3e-14,
        );
        cc(response.support_velocity_m_s[i], rear[i], 3e-14);
    }
    close(
        response.reflection.absorbed_normal_power_w_m2,
        0.1600204802140179672,
        3e-11,
    );
    check_fields(&response);
}

#[test]
fn critical_sh_linear_field_and_deep_evanescent_transport_keep_positive_backing_work() {
    let state = solid();
    let coating = material(1800.0, 6000.0, 3000.0, 32, 0.02);
    let z = 2e6;
    let omega = 1700.0;
    let backing = scalar(C64::from_re(z), C64::from_re(z));
    let h = 0.4;
    let cosine = 5.0_f64.sqrt() / 3.0;
    let prepared = IsotropicElasticCoating::new(&state, &coating, h, backing)
        .unwrap()
        .prepare_incidence(ElasticMode::Sh, cosine, omega / core::f64::consts::TAU)
        .unwrap();
    let response = prepared
        .drive([C64::ZERO, C64::ZERO, C64::from_re(0.001)])
        .unwrap();
    // At qSH=0 the shear traction is constant and velocity is linear.
    let a = C64::new(1.0, -omega * h * z / (1800.0 * 3000.0_f64.powi(2)));
    let d = 5e6 * cosine;
    let rear = C64::from_re(2.0 * d * 0.001) / (a.scale(d) + C64::from_re(z));
    cc(response.support_velocity_m_s[1], rear, 3e-15);
    cc(
        response.reflection.reflected_velocity_m_s[2],
        a * rear - C64::from_re(0.001),
        3e-15,
    );
    check_fields(&response);
    for fraction in [0.0, 0.25, 0.75, 1.0] {
        let expected = C64::new(
            1.0,
            -omega * (h - fraction * h) * z / (1800.0 * 3000.0_f64.powi(2)),
        ) * rear;
        cc(
            response.field_at(fraction * h).unwrap().velocity_m_s[1],
            expected,
            3e-15,
        );
    }
    for h in [0.4, 80.0] {
        let cosine: f64 = 0.5;
        let p = (1.0 - cosine * cosine).sqrt() / 2000.0;
        let q = (p * p - 1.0 / 3000.0_f64.powi(2)).sqrt();
        let d0 = 5e6 * cosine;
        let d1 = C64::new(0.0, 1800.0 * 3000.0_f64.powi(2) * q);
        let reflection = (d1 - C64::from_re(z)) / (d1 + C64::from_re(z));
        let e = (-omega * q * h).exp();
        let forward = C64::from_re(2.0 * d0 * 0.001)
            / (C64::from_re(d0) + d1 + (C64::from_re(d0) - d1) * reflection.scale(e * e));
        let expected = forward.scale(e) * (C64::ONE + reflection);
        let prepared = IsotropicElasticCoating::new(&state, &coating, h, backing)
            .unwrap()
            .prepare_incidence(ElasticMode::Sh, cosine, omega / core::f64::consts::TAU)
            .unwrap();
        let response = prepared
            .drive([C64::ZERO, C64::ZERO, C64::from_re(0.001)])
            .unwrap();
        cc(response.support_velocity_m_s[1] / expected, C64::ONE, 3e-12);
        close(
            response.reflection.absorbed_normal_power_w_m2 / (0.5 * z * expected.norm_sq()),
            1.0,
            3e-12,
        );
        assert!(response.reflection.absorbed_normal_power_w_m2 > 0.0);
    }
}

#[test]
fn stiff_motion_original_domains_and_bounded_cancellation_remain_explicit() {
    let state = solid();
    let omega = 1700.0;
    let h = 4000.0 * core::f64::consts::PI / omega;
    let owner = IsotropicElasticCoating::new(
        &state,
        &state,
        h,
        scalar(C64::from_re(1e30), C64::from_re(1e30)),
    )
    .unwrap();
    let prepared = owner
        .prepare_incidence(ElasticMode::P, 1.0, omega / core::f64::consts::TAU)
        .unwrap();
    let response = prepared
        .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
        .unwrap();
    let rear = -0.002 * (1e7 / 1e30) / (1.0 + 1e7 / 1e30);
    close(response.support_velocity_m_s[2].re / rear, 1.0, 3e-13);
    close(
        response.reflection.absorbed_normal_power_w_m2 / (0.5 * 1e30 * rear * rear),
        1.0,
        3e-13,
    );
    assert!(response.field_at(-1.0).is_err());
    assert!(response.field_at(h.next_up()).is_err());
    assert!(
        prepared
            .drive([C64::from_re(1e9), C64::ZERO, C64::ZERO])
            .is_err()
    );
    assert_eq!(
        prepared
            .drive([C64::ZERO; 3])
            .unwrap()
            .reflection
            .absorbed_normal_power_w_m2,
        0.0
    );
    let limited = material(1800.0, 3000.0, 1500.0, 32, 1e-12);
    let limited =
        IsotropicElasticCoating::new(&state, &limited, 0.4, ElasticSurfaceCondition::TractionFree)
            .unwrap()
            .prepare_incidence(ElasticMode::P, 1.0, 1000.0)
            .unwrap();
    assert!(
        limited
            .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
            .is_err()
    );
    let mut polls = 0;
    assert!(
        owner
            .prepare_incidence_cancellable(ElasticMode::P, 0.8, 1000.0, &mut || {
                polls += 1;
                polls == 20
            })
            .is_err()
    );
    assert_eq!(polls, 20);
    for h in [-1.0, f64::NAN] {
        assert!(
            IsotropicElasticCoating::new(&state, &state, h, ElasticSurfaceCondition::Clamped)
                .is_err()
        );
    }
    for hz in [0.0, -1.0, f64::INFINITY] {
        assert!(owner.prepare_incidence(ElasticMode::P, 1.0, hz).is_err());
    }
}
