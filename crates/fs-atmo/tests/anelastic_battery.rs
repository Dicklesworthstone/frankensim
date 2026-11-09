//! G0/G1/G3/G5: fixed-background mass continuity for complete modal atmospheres.

use fs_atmo::fetch_mass_consistent::ExponentialDensityProfile;
use fs_atmo::mann::MANN_TARGET_V1;
use fs_atmo::mann_field::fixed_step::FixedStepMannField;
use fs_atmo::mann_field::wall::MannPlaneWall;
use fs_atmo::mann_field::{HomogeneousMannField, MannQuadrature, MannSpectralModel};
use fs_atmo::{
    AnelasticFlowSample, Atmosphere, DEC17_AIR, FlatSiteLogLaw, FlowSample, TurbulenceField,
};

fn atmosphere(amplitude: f64) -> Atmosphere {
    Atmosphere {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.005,
            displacement_height_m: 0.02,
            reference_height_m: 10.0,
            reference_speed_mps: 12.0,
        },
        turbulence: TurbulenceField::build(42, 24, amplitude, 30.0, 12.0).unwrap(),
        air: DEC17_AIR,
    }
}

fn profile(rho: f64, k: f64) -> ExponentialDensityProfile {
    ExponentialDensityProfile::new(rho, 10.0, k, [0.0, 200.0]).unwrap()
}

fn same_bits(a: FlowSample, b: FlowSample) {
    for (actual, expected) in
        a.u.iter()
            .chain(a.grad.iter().flatten())
            .zip(b.u.iter().chain(b.grad.iter().flatten()))
    {
        assert_eq!(actual.to_bits(), expected.to_bits());
    }
}

fn check_continuity(k: f64, sample: impl Fn([f64; 3]) -> AnelasticFlowSample) {
    check_continuity_at_step(k, 1e-3, sample);
}

fn check_continuity_at_step(k: f64, step: f64, sample: impl Fn([f64; 3]) -> AnelasticFlowSample) {
    for p in [[3.7, -8.1, 4.9], [-23.0, 18.0, 33.0], [48.0, -9.0, 170.0]] {
        let s = sample(p);
        let diagonal = [s.flow.grad[0][0], s.flow.grad[1][1], s.flow.grad[2][2]];
        let expected_divergence = k * s.flow.u[2];
        let scale = diagonal.iter().map(|x| x.abs()).sum::<f64>() + expected_divergence.abs();
        assert!(scale > 1e-6, "nontrivial divergence cancellation required");
        assert!(
            (diagonal.iter().sum::<f64>() - expected_divergence).abs() < 2e-12 * scale,
            "weighted analytic continuity failed at {p:?}"
        );

        // Differentiate the reported physical flux, independently of the
        // analytic derivative and its density correction.
        let mut mass_divergence = 0.0;
        let mut mass_scale = 0.0;
        for j in 0..3 {
            let mut lo = p;
            let mut hi = p;
            lo[j] -= step;
            hi[j] += step;
            let low = sample(lo);
            let high = sample(hi);
            let derivative = (high.mass_flux_kg_m2_s[j] - low.mass_flux_kg_m2_s[j]) / (2.0 * step);
            mass_divergence += derivative;
            mass_scale += derivative.abs();
            for i in 0..3 {
                let velocity_derivative = (high.flow.u[i] - low.flow.u[i]) / (2.0 * step);
                assert!(
                    (velocity_derivative - s.flow.grad[i][j]).abs()
                        < 1e-7 * (1.0 + s.flow.grad[i][j].abs()),
                    "derivative ({i},{j}) at {p:?}"
                );
            }
        }
        assert!(
            mass_divergence.abs() < 2e-7 * mass_scale,
            "physical mass divergence {mass_divergence:e} at {p:?}, scale {mass_scale:e}"
        );
    }
}

#[test]
fn complete_frozen_and_evolving_fields_conserve_background_mass() {
    let wind = atmosphere(1.4);
    let k = 0.004;
    let density = profile(1.2, k);
    for tick in [0, 103, 2007] {
        check_continuity(k, |[x, y, h]| {
            wind.sample_anelastic_flow(&density, x, y, h, tick).unwrap()
        });
    }
    // The original solenoidal field with this density is a live falsifier:
    // div(rho*u) = -k*rho*w is nonzero wherever it moves vertically.
    let original = wind.sample_flow(3.7, -8.1, 4.9, 103).unwrap();
    let missing_mass_term = k * density.density_kg_m3(4.9).unwrap() * original.u[2];
    assert!(missing_mass_term.abs() > 1e-5);

    let mut evolving = wind.into_evolving(7, vec![0.8; 24]).unwrap();
    for tick in [-3840, -3800, -2900] {
        evolving.turbulence.advance_to(tick).unwrap();
        check_continuity(k, |[x, y, h]| {
            evolving.sample_anelastic_flow(&density, x, y, h).unwrap()
        });
    }
}

#[test]
fn uniform_density_is_bitwise_original_flow_and_mean_is_not_rescaled() {
    let wind = atmosphere(1.4);
    let uniform = profile(1.7, 0.0);
    let evolving = wind.clone().into_evolving(7, vec![0.8; 24]).unwrap();
    for h in [0.0, 0.025, 4.9, 10.0, 170.0] {
        for tick in [0, 2007] {
            let sample = wind
                .sample_anelastic_flow(&uniform, 3.7, -8.1, h, tick)
                .unwrap();
            same_bits(sample.flow, wind.sample_flow(3.7, -8.1, h, tick).unwrap());
            assert_eq!(sample.density_kg_m3, 1.7);
            assert_eq!(sample.mass_flux_kg_m2_s, sample.flow.u.map(|v| 1.7 * v));
        }
        let sample = evolving
            .sample_anelastic_flow(&uniform, 3.7, -8.1, h)
            .unwrap();
        same_bits(sample.flow, evolving.sample_flow(3.7, -8.1, h).unwrap());

        let mean_only = atmosphere(0.0);
        let stratified = profile(1.2, 0.004);
        let sample = mean_only
            .sample_anelastic_flow(&stratified, 3.7, -8.1, h, 103)
            .unwrap();
        assert_eq!(
            sample.flow.u[0].to_bits(),
            mean_only.mean.speed(h).to_bits()
        );
        assert_eq!(
            sample.flow.grad[0][2].to_bits(),
            mean_only.mean.dspeed_dh(h).to_bits()
        );
        assert_eq!(sample.flow.u[2], 0.0);
    }
}

#[test]
fn reference_amplitude_wall_condition_and_density_scale_are_preserved() {
    let wind = atmosphere(1.4);
    let density = profile(1.2, 0.004);
    let heavier = profile(3.6, 0.004);
    for tick in [0, 103, 2007] {
        let reference = wind
            .sample_anelastic_flow(&density, 3.7, -8.1, 10.0, tick)
            .unwrap();
        assert_eq!(
            reference.flow.u,
            wind.sample_flow(3.7, -8.1, 10.0, tick).unwrap().u
        );
        let wall = wind
            .sample_anelastic_flow(&density, 3.7, -8.1, 0.0, tick)
            .unwrap();
        assert_eq!(wall.flow.u[2], 0.0);
        assert_eq!(wall.mass_flux_kg_m2_s[2], 0.0);
        for h in [4.9, 33.0, 170.0] {
            let a = wind
                .sample_anelastic_flow(&density, 3.7, -8.1, h, tick)
                .unwrap();
            let b = wind
                .sample_anelastic_flow(&heavier, 3.7, -8.1, h, tick)
                .unwrap();
            for i in 0..3 {
                assert!((a.flow.u[i] - b.flow.u[i]).abs() < 2e-14 * (1.0 + a.flow.u[i].abs()));
                assert!(
                    (b.mass_flux_kg_m2_s[i] - 3.0 * a.mass_flux_kg_m2_s[i]).abs()
                        < 2e-14 * (1.0 + b.mass_flux_kg_m2_s[i].abs())
                );
            }
        }
    }
}

#[test]
fn stratified_queries_preserve_ou_state_clock_and_replay() {
    let mut observed = atmosphere(1.4).into_evolving(7, vec![0.8; 24]).unwrap();
    let mut unobserved = observed.clone();
    let density = profile(1.2, 0.004);
    for tick in [-3840, -3800, -2900] {
        observed.turbulence.advance_to(tick).unwrap();
        unobserved.turbulence.advance_to(tick).unwrap();
        for h in [0.0, 4.9, 33.0, 170.0] {
            observed
                .sample_anelastic_flow(&density, 3.7, -8.1, h)
                .unwrap();
        }
        assert_eq!(
            observed.turbulence.checkpoint(),
            unobserved.turbulence.checkpoint()
        );
        let checkpoint = observed.turbulence.checkpoint();
        let mut resumed = observed.clone();
        resumed.turbulence =
            fs_atmo::evolving::EvolvingTurbulenceField::resume(checkpoint).unwrap();
        assert_eq!(
            observed
                .sample_anelastic_flow(&density, 3.7, -8.1, 33.0)
                .unwrap(),
            resumed
                .sample_anelastic_flow(&density, 3.7, -8.1, 33.0)
                .unwrap()
        );
    }
}

#[test]
fn invalid_density_queries_refuse_without_changing_the_evolving_field() {
    let wind = atmosphere(1.4);
    let evolving = wind.clone().into_evolving(7, vec![0.8; 24]).unwrap();
    let before = evolving.turbulence.checkpoint();
    let density = profile(1.2, 0.004);
    assert_eq!(
        wind.sample_anelastic_flow(&density, 0.0, 0.0, 201.0, 0)
            .unwrap_err()
            .code,
        "density-profile-domain"
    );
    assert_eq!(
        evolving
            .sample_anelastic_flow(&density, 0.0, 0.0, 201.0)
            .unwrap_err()
            .code,
        "density-profile-domain"
    );
    assert_eq!(
        wind.sample_anelastic_flow(&density, f64::NAN, 0.0, 10.0, 0)
            .unwrap_err()
            .code,
        "non-finite-input"
    );
    assert_eq!(
        evolving
            .sample_anelastic_flow(&density, 0.0, f64::INFINITY, 10.0)
            .unwrap_err()
            .code,
        "non-finite-input"
    );
    let excessive_mass_flux = profile(1e308, 0.0);
    assert_eq!(
        wind.sample_anelastic_flow(&excessive_mass_flux, 0.0, 0.0, 10.0, 0)
            .unwrap_err()
            .code,
        "anelastic-flow-range"
    );
    assert_eq!(
        evolving
            .sample_anelastic_flow(&excessive_mass_flux, 0.0, 0.0, 10.0)
            .unwrap_err()
            .code,
        "anelastic-flow-range"
    );
    assert_eq!(evolving.turbulence.checkpoint(), before);
}

fn mann_field() -> HomogeneousMannField {
    let model = MannSpectralModel::build(
        41,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 24,
            min_kl: 0.1,
            max_kl: 12.0,
        },
        [12.0, -0.3, 0.0],
    )
    .unwrap();
    HomogeneousMannField::at_anchor(model, 57, vec![0.8; 24]).unwrap()
}

#[test]
fn full_mann_tensor_fields_conserve_mass_with_both_clocks_and_wall_observations() {
    let mut original = mann_field();
    original.advance_to(-3800).unwrap();
    let mut fixed =
        FixedStepMannField::from_snapshot(original.clone(), 123, 0.0079, 88 * 24).unwrap();
    let k = 0.004;
    let density = profile(1.2, k);
    let wall = MannPlaneWall::new(3.0).unwrap();
    for step in [0, 3] {
        fixed.advance_to(step, 4 * 24 * 3, || false).unwrap();
        // The canonical Mann spectrum reaches shorter wavelengths than the
        // legacy 30 m modal basis, so resolve its gradients with a finer stencil.
        check_continuity_at_step(k, 1e-5, |p| original.sample_anelastic(p, &density).unwrap());
        check_continuity_at_step(k, 1e-5, |p| fixed.sample_anelastic(p, &density).unwrap());
        check_continuity_at_step(k, 1e-5, |p| {
            original
                .sample_wall_blocked_anelastic(p, wall, &density)
                .unwrap()
        });
        check_continuity_at_step(k, 1e-5, |p| {
            fixed
                .sample_wall_blocked_anelastic(p, wall, &density)
                .unwrap()
        });
    }
    // The old homogeneous velocity combined directly with this density fails
    // the mass constraint; the product-rule term is physically necessary.
    let p = [3.7, -8.1, 4.9];
    let old = original.sample(p).unwrap();
    assert!((k * density.density_kg_m3(p[2]).unwrap() * old.u[2]).abs() > 1e-5);
}

#[test]
fn mann_density_observations_retain_uniform_bits_reference_amplitude_and_wall_flux() {
    let original = mann_field();
    let fixed = FixedStepMannField::from_snapshot(original.clone(), 123, 0.0079, 88 * 24).unwrap();
    let uniform = profile(1.7, 0.0);
    let density = profile(1.2, 0.004);
    let heavier = profile(3.6, 0.004);
    let wall = MannPlaneWall::new(3.0).unwrap();
    for height in [3.0, 10.0, 33.0, 170.0] {
        let p = [3.7, -8.1, height];
        same_bits(
            original.sample_anelastic(p, &uniform).unwrap().flow,
            original.sample(p).unwrap(),
        );
        same_bits(
            fixed.sample_anelastic(p, &uniform).unwrap().flow,
            fixed.sample(p).unwrap(),
        );
        same_bits(
            original
                .sample_wall_blocked_anelastic(p, wall, &uniform)
                .unwrap()
                .flow,
            original.sample_wall_blocked(p, wall).unwrap(),
        );
        same_bits(
            fixed
                .sample_wall_blocked_anelastic(p, wall, &uniform)
                .unwrap()
                .flow,
            fixed.sample_wall_blocked(p, wall).unwrap(),
        );
        for blocked in [false, true] {
            let observe = |rho: &ExponentialDensityProfile| {
                if blocked {
                    original.sample_wall_blocked_anelastic(p, wall, rho)
                } else {
                    original.sample_anelastic(p, rho)
                }
                .unwrap()
            };
            let a = observe(&density);
            let b = observe(&heavier);
            let unscaled = observe(&uniform);
            let scale = (0.004 * (height - 10.0)).exp();
            for i in 0..3 {
                assert!((a.flow.u[i] - scale * unscaled.flow.u[i]).abs() < 2e-14);
                assert!((a.flow.u[i] - b.flow.u[i]).abs() < 2e-14);
                assert!((b.mass_flux_kg_m2_s[i] - 3.0 * a.mass_flux_kg_m2_s[i]).abs() < 4e-14);
            }
            if height == 10.0 {
                assert_eq!(a.flow.u, unscaled.flow.u);
            }
            if blocked && height == wall.height_m() {
                assert_eq!(a.flow.u[2], 0.0);
                assert_eq!(a.mass_flux_kg_m2_s[2], 0.0);
            }
            let fixed_sample = if blocked {
                fixed.sample_wall_blocked_anelastic(p, wall, &density)
            } else {
                fixed.sample_anelastic(p, &density)
            }
            .unwrap();
            assert_eq!(a, fixed_sample);
        }
    }
}

#[test]
fn mann_density_queries_are_read_only_replayable_and_refuse_unsupported_states() {
    let original = mann_field();
    let mut fixed =
        FixedStepMannField::from_snapshot(original.clone(), 123, 0.0079, 88 * 24).unwrap();
    fixed.advance_to(3, 4 * 24 * 3, || false).unwrap();
    let before_original = original.checkpoint();
    let before_fixed = fixed.checkpoint();
    let density = profile(1.2, 0.004);
    let wall = MannPlaneWall::new(3.0).unwrap();
    let resumed_original = HomogeneousMannField::resume(before_original.clone()).unwrap();
    let resumed_fixed = FixedStepMannField::resume(before_fixed.clone(), 88 * 24).unwrap();
    let p = [3.7, -8.1, 33.0];
    assert_eq!(
        original.sample_anelastic(p, &density),
        resumed_original.sample_anelastic(p, &density)
    );
    assert_eq!(
        fixed.sample_wall_blocked_anelastic(p, wall, &density),
        resumed_fixed.sample_wall_blocked_anelastic(p, wall, &density)
    );
    for bad in [[0.0, 0.0, 201.0], [0.0, 0.0, -1.0], [f64::NAN, 0.0, 10.0]] {
        assert!(original.sample_anelastic(bad, &density).is_err());
        assert!(fixed.sample_anelastic(bad, &density).is_err());
    }
    let below_wall = [0.0, 0.0, 2.0];
    assert!(
        original
            .sample_wall_blocked_anelastic(below_wall, wall, &density)
            .is_err()
    );
    assert!(
        fixed
            .sample_wall_blocked_anelastic(below_wall, wall, &density)
            .is_err()
    );
    let extreme = ExponentialDensityProfile::new(1.0, 0.0, 3.6, [0.0, 200.0]).unwrap();
    for failure in [
        original.sample_anelastic([0.0, 0.0, 200.0], &extreme),
        fixed.sample_anelastic([0.0, 0.0, 200.0], &extreme),
    ] {
        assert_eq!(failure.unwrap_err().code, "anelastic-flow-range");
    }
    assert_eq!(original.checkpoint(), before_original);
    assert_eq!(fixed.checkpoint(), before_fixed);
}
