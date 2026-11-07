use super::*;

fn weakening(cutoff: f64) -> FrictionLaw {
    FrictionLaw::VelocityDependent {
        static_mu: 0.8,
        mu_zero: 0.6,
        slope_per_speed: -0.6 / cutoff,
    }
}

#[test]
fn intersecting_active_disks_match_independent_whole_surface_traction() {
    // Both positions of the stationary point and both exterior ray topologies.
    for (d, k) in [
        (0.25, 0.9),
        (0.8, 0.5),
        (1.0, 0.7),
        (1.4, 0.6),
        (1.4, 1.2),
        (10.0, 9.5),
    ] {
        let law = weakening(k);
        let model = HertzFrictionPatch::new(2.0, 0.4, law.clone()).unwrap();
        let slip = [0.6 * d, -0.8 * d];
        let actual = model.evaluate(slip, 2.5, rate_grid(24, 24)).unwrap();
        let reference = spatial_rate_law(&law, [2.0, 0.4], slip, 2.5, [512, 1024]);
        close(actual.force_n[0], reference[0], 3e-6);
        close(actual.force_n[1], reference[1], 3e-6);
        close(actual.spin_moment_nm, reference[2], 2e-6);
        close(actual.dissipated_power_w, reference[3], 3e-6);
        close(actual.power_residual_w, 0.0, 3e-14);
        assert!(actual.dissipated_power_w > 0.0);
        let reversed = model
            .evaluate(slip.map(|v| -v), -2.5, rate_grid(24, 24))
            .unwrap();
        assert_eq!(reversed.force_n, actual.force_n.map(|v| -v));
        assert_eq!(reversed.spin_moment_nm, -actual.spin_moment_nm);
        assert_eq!(reversed.dissipated_power_w, actual.dissipated_power_w);
        assert_eq!(
            actual,
            model.evaluate(slip, 2.5, rate_grid(24, 24)).unwrap()
        );
    }
}

#[test]
fn shrinking_exterior_crescents_retain_the_independent_cubic_load_limit() {
    // At fixed d>1, k=d-1+eps. Local Cartesian edge coordinates give
    // p~sqrt(2z), mu/mu0~(eps-z-d*y^2/[2(d-1)])/(d-1).
    // Integrating that parabola gives F=-N*mu0*eps^3/[4sqrt(d(d-1))],
    // M=-a*F and D=-a*|omega|*(d-1)*F to leading order.
    let d = 1.5_f64;
    for requested in [1e-4, 1e-6, 1e-8] {
        let k = d - 1.0 + requested;
        let eps = k - (d - 1.0);
        let model = HertzFrictionPatch::new(2.0, 0.4, weakening(k)).unwrap();
        let actual = model.evaluate([d, 0.0], 2.5, rate_grid(12, 12)).unwrap();
        let force = -2.0 * 0.6 * eps.powi(3) / (4.0 * (d * (d - 1.0)).sqrt());
        let relative = 20.0 * eps + 2e-6;
        close(actual.force_n[0], force, force.abs() * relative);
        close(
            actual.spin_moment_nm,
            -0.4 * force,
            0.4 * force.abs() * relative,
        );
        close(
            actual.dissipated_power_w,
            -(d - 1.0) * force,
            force.abs() * relative,
        );
        close(actual.power_residual_w, 0.0, force.abs() * 2e-14);
    }
    let closed = HertzFrictionPatch::new(2.0, 0.4, weakening(d - 1.0)).unwrap();
    let zero = closed.evaluate([d, 0.0], 2.5, rate_grid(1, 1)).unwrap();
    assert_eq!(zero.force_n, [0.0; 2]);
    assert_eq!(zero.spin_moment_nm, 0.0);
    assert_eq!(zero.dissipated_power_w, 0.0);
}

#[test]
fn microscopic_active_disks_on_the_rim_keep_positive_force_and_heat() {
    // d=1, k<<1: p~sqrt(2*q*cos(phi)). Independent beta integrals give
    // integral cos(phi)^(3/2)=sqrt(pi)*Gamma(5/4)/(2*Gamma(7/4)),
    // integral cos(phi)^(1/2)=sqrt(pi)*Gamma(3/4)/(2*Gamma(5/4)).
    let a_force = 0.874_019_184_764_04;
    let a_heat = 1.198_140_234_735_592;
    for k in [1e-8_f64, 1e-12, 1e-20] {
        let edge = k.recip();
        let model = HertzFrictionPatch::new(2.0, 1.0, weakening(1.0)).unwrap();
        let actual = model
            .evaluate([edge, 0.0], edge, rate_grid(24, 24))
            .unwrap();
        let scale = (3.0 / PI) * 2.0 * 0.6 * 2.0_f64.sqrt() * k.powf(2.5);
        let force = -scale * (4.0 / 35.0) * a_force;
        let heat = scale * (4.0 / 63.0) * a_heat;
        close(actual.force_n[0], force, force.abs() * 2e-6);
        close(actual.spin_moment_nm, -force, force.abs() * 2e-6);
        close(actual.dissipated_power_w, heat, heat * 2e-6);
        assert!(actual.dissipated_power_w > 0.0);
        // The large opposing translation/spin powers can cancel beyond their
        // representable relative resolution; directly integrated heat is kept.
        assert!(actual.power_residual_w.is_finite());
    }
}

#[test]
fn active_intersection_refines_and_crosses_geometric_branches_continuously() {
    let model = HertzFrictionPatch::new(2.0, 0.4, weakening(1.2)).unwrap();
    let reference = model.evaluate([1.4, 0.0], 2.5, rate_grid(64, 64)).unwrap();
    let mut errors = Vec::new();
    for panels in [1, 2, 4] {
        let result = model
            .evaluate([1.4, 0.0], 2.5, rate_grid(panels, panels))
            .unwrap();
        errors.push(
            (result.force_n[0] - reference.force_n[0])
                .abs()
                .max((result.spin_moment_nm - reference.spin_moment_nm).abs())
                .max((result.dissipated_power_w - reference.dissipated_power_w).abs()),
        );
    }
    assert!(errors[1] < errors[0] / 3.0, "{errors:?}");
    assert!(errors[2] < errors[1] / 3.0, "{errors:?}");
    for (cutoff, d) in [(0.4, 0.6), (0.25, 1.0)] {
        let model = HertzFrictionPatch::new(2.0, 0.4, weakening(cutoff)).unwrap();
        let a = model
            .evaluate([d - 1e-8, 0.0], 2.5, rate_grid(32, 32))
            .unwrap();
        let b = model
            .evaluate([d + 1e-8, 0.0], 2.5, rate_grid(32, 32))
            .unwrap();
        close(a.force_n[0], b.force_n[0], 2e-7);
        close(a.spin_moment_nm, b.spin_moment_nm, 2e-7);
        close(a.dissipated_power_w, b.dissipated_power_w, 2e-7);
    }
}

#[test]
fn narrow_exterior_angles_retain_the_independent_half_patch_limit() {
    // k=d>>1 leaves mu~mu0*y/d on the leading half of the unit patch.
    // The Hertz marginal is 3(1-y^2)/4: <y_+>=3/16 and <y_+^2>=1/10.
    // atan2 geometry must retain this region when its crossing cosine rounds
    // to one; a direct acos would erase the angular interval.
    for d in [1e4_f64, 1e8] {
        let model = HertzFrictionPatch::new(2.0, 0.4, weakening(d)).unwrap();
        let actual = model.evaluate([d, 0.0], 2.5, rate_grid(24, 24)).unwrap();
        let force = -3.0 * 2.0 * 0.6 / (16.0 * d);
        let moment = 2.0 * 0.6 * 0.4 / (10.0 * d);
        let heat = 3.0 * 2.0 * 0.6 / 16.0;
        let relative = 5.0 / d + 3e-6;
        close(actual.force_n[0], force, force.abs() * relative);
        close(actual.spin_moment_nm, moment, moment * relative);
        close(actual.dissipated_power_w, heat, heat * relative);
        close(actual.power_residual_w, 0.0, 3e-14);
    }
}
