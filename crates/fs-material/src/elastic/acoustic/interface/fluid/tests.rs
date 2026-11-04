use super::*;
use std::cell::Cell;

fn solid(rho: f64, cp: f64, cs: f64, ceiling: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        rho,
        rho * (cp * cp - 4.0 * cs * cs / 3.0),
        rho * cs * cs,
        ContentHash([61; 32]),
        ceiling,
    )
    .unwrap()
}
fn limits() -> FluidDriveLimits {
    FluidDriveLimits::try_new(0.05, 0.05).unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.17e} != {b:.17e}, tolerance {tol:e}"
    );
}
fn cc(a: C64, b: C64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:?} != {b:?}, tolerance {tol:e}");
}
fn channel(i: usize) -> FluidElasticChannel {
    if i == 0 {
        FluidElasticChannel::Fluid
    } else {
        FluidElasticChannel::Solid([ElasticMode::P, ElasticMode::Sv, ElasticMode::Sh][i - 1])
    }
}

#[test]
fn normal_incidence_matches_acoustic_impedances_and_shear_free_reflection() {
    let fluid = FluidAcousticState::try_new(998.0, 1480.0).unwrap();
    let body = solid(2700.0, 6100.0, 3100.0, 0.02);
    let owner = FluidElasticInterface::new(&fluid, &body).unwrap();
    let p = owner
        .prepare_incidence(FluidElasticChannel::Fluid, 1.0)
        .unwrap();
    let zf = fluid.characteristic_impedance_pa_s_m();
    let zs = owner.impedances[0];
    let r = (zs - zf) / (zs + zf);
    let from_fluid = p.velocity_coefficients(channel(0)).unwrap();
    let from_solid = p.velocity_coefficients(channel(1)).unwrap();
    cc(from_fluid[0], C64::from_re(r), 2e-15);
    cc(from_fluid[1], C64::from_re(2.0 * zf / (zs + zf)), 2e-15);
    cc(from_solid[0], C64::from_re(2.0 * zs / (zs + zf)), 2e-15);
    cc(from_solid[1], C64::from_re(-r), 2e-15);
    for i in 2..4 {
        let coefficients = p.velocity_coefficients(channel(i)).unwrap();
        for j in 0..4 {
            cc(
                coefficients[j],
                if i == j { C64::ONE } else { C64::ZERO },
                2e-15,
            );
        }
    }
    for i in 0..4 {
        let mut drive = [C64::ZERO; 4];
        drive[i] = C64::new(0.003, -0.004);
        let response = p.drive(drive, limits()).unwrap();
        close(
            response.total_incoming_power_w_m2,
            response.total_outgoing_power_w_m2,
            2e-12 * response.total_incoming_power_w_m2,
        );
        cc(
            response.boundary_velocity_m_s[0][2],
            response.boundary_velocity_m_s[1][2],
            1e-16,
        );
        cc(
            response.fluid_pressure_pa,
            -response.solid_traction_pa[2],
            1e-9,
        );
        cc(response.solid_traction_pa[0], C64::ZERO, 1e-9);
        cc(response.solid_traction_pa[1], C64::ZERO, 1e-9);
    }
}

fn sqrt_branch(x: f64) -> C64 {
    if x >= 0.0 {
        C64::from_re(x.sqrt())
    } else {
        C64::new(0.0, (-x).sqrt())
    }
}

// Independent isotropic Lamé stress equations and the shear-free null vector.
// This oracle neither calls native wave traces nor solves the production LU.
fn fluid_incidence_oracle(rho_f: f64, cf: f64, rho: f64, cp: f64, cs: f64, cos: f64) -> [C64; 3] {
    let p = (1.0 - cos * cos).sqrt() / cf;
    let sp = p * cp;
    let ss = p * cs;
    let qp = sqrt_branch(1.0 - sp * sp);
    let qs = sqrt_branch(1.0 - ss * ss);
    let mu = rho * cs * cs;
    let lambda = rho * (cp * cp - 2.0 * cs * cs);
    let tx_p = qp.scale(-2.0 * mu * p);
    let tx_s = (qs * qs - C64::from_re(ss * ss)).scale(-mu / cs);
    let tz_p = C64::from_re(-lambda / cp) - (qp * qp).scale(2.0 * mu / cp);
    let tz_s = qs.scale(2.0 * mu * ss / cs);
    let vz = qp * tx_s + tx_p.scale(ss);
    let tz = tz_p * tx_s - tz_s * tx_p;
    let z = rho_f * cf;
    let scale = C64::from_re(2.0 * z * cos) / (vz.scale(z) - tz.scale(cos));
    [
        (tz.scale(cos) + vz.scale(z)) / (tz.scale(cos) - vz.scale(z)),
        tx_s * scale,
        -tx_p * scale,
    ]
}

#[test]
fn oblique_fluid_excitation_matches_independent_stress_mode_conversion_and_evanescence() {
    let fluid = FluidAcousticState::try_new(1000.0, 1200.0).unwrap();
    let body = solid(2600.0, 3400.0, 1800.0, 0.02);
    let owner = FluidElasticInterface::new(&fluid, &body).unwrap();
    for cosine in [1.0, 0.99, 0.97, 0.9, 0.8, 0.7, 0.4, 0.05] {
        let p = owner.prepare_incidence(channel(0), cosine).unwrap();
        let c = p.velocity_coefficients(channel(0)).unwrap();
        let expected = fluid_incidence_oracle(1000.0, 1200.0, 2600.0, 3400.0, 1800.0, cosine);
        for i in 0..3 {
            cc(c[i], expected[i], 2e-11 * expected[i].abs().max(1.0));
        }
        let response = p
            .drive(
                [C64::new(0.004, 0.003), C64::ZERO, C64::ZERO, C64::ZERO],
                limits(),
            )
            .unwrap();
        let input = response.total_incoming_power_w_m2;
        close(response.total_outgoing_power_w_m2, input, 2e-12 * input);
        if cosine < 0.7 {
            assert!(p.outgoing[1].normal_slowness_s_m.im > 0.0);
            assert!(p.outgoing[2].normal_slowness_s_m.im > 0.0);
            assert_eq!(response.outgoing_normal_power_w_m2[1], 0.0);
            assert_eq!(response.outgoing_normal_power_w_m2[2], 0.0);
            close(c[0].abs(), 1.0, 2e-14);
        }
    }
}

fn tensor_traction(state: &ElasticAcousticState, w: ElasticWaveTrace, p: f64) -> [C64; 3] {
    let s = [C64::from_re(p), C64::ZERO, w.normal_slowness_s_m];
    let e: [[C64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| (s[i] * w.velocity[j] + s[j] * w.velocity[i]).scale(-0.5))
    });
    let mandel = [
        e[0][0],
        e[1][1],
        e[2][2],
        e[0][1].scale(2.0_f64.sqrt()),
        e[1][2].scale(2.0_f64.sqrt()),
        e[2][0].scale(2.0_f64.sqrt()),
    ];
    let stress: [C64; 6] = core::array::from_fn(|i| {
        (0..6).fold(C64::ZERO, |v, j| {
            v + mandel[j].scale(state.stiffness_mandel_pa()[i][j])
        })
    });
    [
        stress[5].scale(0.5_f64.sqrt()),
        stress[4].scale(0.5_f64.sqrt()),
        stress[2],
    ]
}

#[test]
fn coherent_two_sided_sources_preserve_tensor_continuity_unitarity_and_reciprocity() {
    let fluid = FluidAcousticState::try_new(1200.0, 1800.0).unwrap();
    let body = solid(3100.0, 3500.0, 1700.0, 0.02);
    let p = FluidElasticInterface::new(&fluid, &body)
        .unwrap()
        .prepare_incidence(channel(0), 0.97)
        .unwrap();
    let weights = p.incoming.map(|w| w.normal_power_factor_pa_s_m);
    for a in 0..4 {
        for b in 0..4 {
            let dot = (0..4).fold(C64::ZERO, |sum, i| {
                sum + (p.scattering[i][a].conj() * p.scattering[i][b])
                    .scale(weights[i] / (weights[a] * weights[b]).sqrt())
            });
            cc(dot, if a == b { C64::ONE } else { C64::ZERO }, 2e-13);
            close(
                p.scattering[a][b].abs() * (weights[a] / weights[b]).sqrt(),
                p.scattering[b][a].abs() * (weights[b] / weights[a]).sqrt(),
                2e-13,
            );
        }
    }
    let incoming = [
        C64::new(0.003, -0.002),
        C64::new(-0.004, 0.001),
        C64::new(0.001, 0.002),
        C64::new(0.002, 0.003),
    ];
    let response = p.drive(incoming, limits()).unwrap();
    let mut traction = [C64::ZERO; 3];
    let mut velocity = [C64::ZERO; 3];
    for mode in 1..4 {
        for (w, a) in [
            (p.incoming[mode], incoming[mode]),
            (p.outgoing[mode], response.outgoing_velocity_m_s[mode]),
        ] {
            let t = tensor_traction(&body, w, p.tangential_slowness_s_m());
            for j in 0..3 {
                traction[j] = traction[j] + t[j] * a;
                velocity[j] = velocity[j] + w.velocity[j] * a;
            }
        }
    }
    cc(traction[0], C64::ZERO, 3e-10);
    cc(traction[1], C64::ZERO, 3e-10);
    cc(traction[2], -response.fluid_pressure_pa, 3e-10);
    cc(velocity[2], response.boundary_velocity_m_s[0][2], 3e-17);
    let flux = -0.5
        * (0..3)
            .map(|j| (traction[j] * velocity[j].conj()).re)
            .sum::<f64>();
    close(
        flux,
        response.boundary_normal_flux_w_m2[0],
        3e-12 * response.total_incoming_power_w_m2,
    );
    // Slip is a physical degree of freedom, not a welded fluid constraint.
    assert!((velocity[0] - response.boundary_velocity_m_s[0][0]).abs() > 1e-5);
    // Coherent cross-source cancellation changes each outgoing channel.
    let fluid_alone = p
        .drive([incoming[0], C64::ZERO, C64::ZERO, C64::ZERO], limits())
        .unwrap();
    let solids_alone = p
        .drive([C64::ZERO, incoming[1], incoming[2], incoming[3]], limits())
        .unwrap();
    cc(
        response.outgoing_velocity_m_s[0],
        fluid_alone.outgoing_velocity_m_s[0] + solids_alone.outgoing_velocity_m_s[0],
        2e-17,
    );
    assert!(
        (response.outgoing_normal_power_w_m2[0]
            - fluid_alone.outgoing_normal_power_w_m2[0]
            - solids_alone.outgoing_normal_power_w_m2[0])
            .abs()
            > 1e-4
    );
}

#[test]
fn solid_incidence_keeps_evanescent_fluid_pressure_without_false_radiated_power() {
    let fluid = FluidAcousticState::try_new(1600.0, 6000.0).unwrap();
    let body = solid(2700.0, 4000.0, 2000.0, 0.02);
    let p = FluidElasticInterface::new(&fluid, &body)
        .unwrap()
        .prepare_incidence(channel(2), 0.8)
        .unwrap();
    assert!(p.incoming_mode(channel(0)).is_none());
    assert!(p.incoming_mode(channel(1)).is_none());
    assert!(p.outgoing_mode(channel(0)).normal_slowness_s_m.im < 0.0);
    let response = p
        .drive(
            [
                C64::ZERO,
                C64::ZERO,
                C64::new(0.003, 0.004),
                C64::from_re(0.001),
            ],
            limits(),
        )
        .unwrap();
    assert_eq!(response.outgoing_normal_power_w_m2[0], 0.0);
    assert_eq!(response.outgoing_normal_power_w_m2[1], 0.0);
    assert!(response.fluid_pressure_pa.abs() > 1.0);
    close(
        response.boundary_normal_flux_w_m2[0],
        0.0,
        1e-13 * response.total_incoming_power_w_m2,
    );
    assert!(
        p.drive(
            [C64::from_re(1e-6), C64::ZERO, C64::ZERO, C64::ZERO],
            limits()
        )
        .is_err()
    );
}

#[test]
fn critical_neighborhoods_and_common_density_scaling_keep_the_physical_branches() {
    let fluid = FluidAcousticState::try_new(1000.0, 1100.0).unwrap();
    let body = solid(2500.0, 2200.0, 1100.0, 0.02);
    let owner = FluidElasticInterface::new(&fluid, &body).unwrap();
    let critical = (0.75_f64).sqrt();
    for cosine in [
        1.0,
        0.99,
        critical,
        critical.next_up(),
        critical.next_down(),
        0.5,
        1e-8,
    ] {
        let p = owner.prepare_incidence(channel(0), cosine).unwrap();
        let response = p
            .drive(
                [C64::from_re(0.001), C64::ZERO, C64::ZERO, C64::ZERO],
                limits(),
            )
            .unwrap();
        close(
            response.total_incoming_power_w_m2,
            response.total_outgoing_power_w_m2,
            1e-11 * response.total_incoming_power_w_m2,
        );
        if cosine == critical.next_up() {
            assert!(p.outgoing[1].normal_slowness_s_m.re > 0.0);
        }
        if cosine == critical.next_down() {
            assert!(p.outgoing[1].normal_slowness_s_m.im > 0.0);
        }
    }
    let original = owner
        .prepare_incidence(channel(0), 0.96)
        .unwrap()
        .velocity_coefficients(channel(0))
        .unwrap();
    for scale in [1e-6, 1e6] {
        let f = FluidAcousticState::try_new(1000.0 * scale, 1100.0).unwrap();
        let s = solid(2500.0 * scale, 2200.0, 1100.0, 0.02);
        let scaled = FluidElasticInterface::new(&f, &s)
            .unwrap()
            .prepare_incidence(channel(0), 0.96)
            .unwrap()
            .velocity_coefficients(channel(0))
            .unwrap();
        for i in 0..4 {
            cc(scaled[i], original[i], 3e-13);
        }
    }
}

#[derive(Debug)]
struct OriginalFluid {
    state: FluidAcousticState,
    evaluations: Cell<usize>,
}
impl FluidAcoustics for OriginalFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.evaluations.set(self.evaluations.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn original_owners_limits_and_cancellation_remain_explicit_and_reusable() {
    let fluid = OriginalFluid {
        state: FluidAcousticState::try_new(1100.0, 1500.0).unwrap(),
        evaluations: Cell::new(0),
    };
    let body = solid(3000.0, 5000.0, 2500.0, 0.02);
    let owner = FluidElasticInterface::new(&fluid, &body).unwrap();
    let p = owner.prepare_incidence(channel(0), 0.98).unwrap();
    assert!(core::ptr::eq(p.fluid_source(), &fluid));
    assert!(core::ptr::eq(p.solid_state(), &body));
    assert_eq!(p.fluid_state(), fluid.state);
    let incoming = [C64::from_re(0.01), C64::ZERO, C64::ZERO, C64::ZERO];
    let expected = p.drive(incoming, limits()).unwrap();
    for _ in 0..8 {
        assert_eq!(
            p.drive(incoming, limits()).unwrap().outgoing_velocity_m_s,
            expected.outgoing_velocity_m_s
        );
    }
    assert_eq!(fluid.evaluations.get(), 1);
    assert!(
        p.drive(incoming, FluidDriveLimits::try_new(1e-10, 0.1).unwrap())
            .is_err()
    );
    assert!(
        p.drive(incoming, FluidDriveLimits::try_new(0.1, 1e-10).unwrap())
            .is_err()
    );
    let fragile = solid(3000.0, 5000.0, 2500.0, 1e-10);
    let fragile_owner = FluidElasticInterface::new(&fluid, &fragile).unwrap();
    assert!(
        fragile_owner
            .prepare_incidence(channel(0), 1.0)
            .unwrap()
            .drive(incoming, limits())
            .is_err()
    );
    for bad in [f64::NAN, 0.0, -1.0, 1.0001] {
        assert!(owner.prepare_incidence(channel(0), bad).is_err());
    }
    for bad in [f64::NAN, 0.0, -1.0, 1.0] {
        assert!(FluidDriveLimits::try_new(bad, 0.1).is_err());
    }
    for stop in 1..8 {
        let mut polls = 0;
        assert!(
            owner
                .prepare_incidence_cancellable(channel(0), 0.98, &mut || {
                    polls += 1;
                    polls == stop
                })
                .is_err()
        );
        assert_eq!(
            owner
                .prepare_incidence(channel(0), 0.98)
                .unwrap()
                .velocity_coefficients(channel(0)),
            p.velocity_coefficients(channel(0))
        );
    }
}
