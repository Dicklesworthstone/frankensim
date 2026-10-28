use super::*;

fn support() -> WetPoreSupport {
    WetPoreSupport {
        liquid_saturation: (0.0, 1.0),
        total_porosity: (0.01, 1.0),
        mean_gas_pore_radius_m: (1e-9, 1e-3),
    }
}
fn law(diffusivity: GasDiffusivityLaw) -> WetPoreGasLaw {
    WetPoreGasLaw::try_new(
        0.1,
        0.2,
        0.5,
        diffusivity,
        support(),
        "authored analytic fixture; no measured material authority",
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

#[test]
fn vg_mualem_matches_independent_m_half_polynomial_and_residual_normalization() {
    let law = law(GasDiffusivityLaw::ConstantTransportTortuosity { factor: 2.0 });
    for (saturation, effective) in [(0.275, 0.25), (0.45, 0.5), (0.625, 0.75)] {
        let state = law.sample(saturation, 0.4, 1e-6).unwrap();
        // With m=1/2 the capillary term reduces exactly to 1-Se^2.
        let expected = f64::sqrt(1.0 - effective) * (1.0 - effective * effective);
        close(state.gas_relative_permeability(), expected, 2e-14);
        close(state.effective_liquid_saturation(), effective, 2e-15);
        close(state.gas_filled_porosity(), 0.4 * (1.0 - saturation), 2e-15);
        close(
            state.effective_diffusivity_factor(),
            0.2 * (1.0 - saturation),
            2e-15,
        );
        assert_eq!(state.transport_tortuosity_factor(), Some(2.0));
        assert_eq!(state.law(), &law);
    }
}

#[test]
fn exact_dry_residual_and_fully_wet_limits_keep_mobility_and_diffusion_distinct() {
    let law = law(GasDiffusivityLaw::ConstantTransportTortuosity { factor: 4.0 });
    for saturation in [0.0, 0.05, 0.1] {
        let state = law.sample(saturation, 0.4, 1e-6).unwrap();
        assert_eq!(state.gas_relative_permeability(), 1.0);
        assert_eq!(state.effective_liquid_saturation(), 0.0);
    }
    for saturation in [0.8, 0.9] {
        let state = law.sample(saturation, 0.4, 1e-6).unwrap();
        assert_eq!(state.gas_relative_permeability(), 0.0);
        assert_eq!(state.effective_liquid_saturation(), 1.0);
        assert!(
            state.effective_diffusivity_factor() > 0.0,
            "residual immobility is not an inferred diffusion threshold"
        );
    }
    let flooded = law.sample(1.0, 0.4, 1e-6).unwrap();
    assert_eq!(flooded.gas_relative_permeability(), 0.0);
    assert_eq!(flooded.gas_filled_porosity(), 0.0);
    assert_eq!(flooded.effective_diffusivity_factor(), 0.0);
    assert_eq!(flooded.transport_tortuosity_factor(), None);
}

#[test]
fn millington_quirk_has_independent_exact_cube_root_fixtures() {
    let law = law(GasDiffusivityLaw::MillingtonQuirk1961);
    let dry = law.sample(0.0, 0.125, 1e-6).unwrap();
    // phi=1/8 gives phi^(4/3)=1/16 and transport tau=2.
    close(dry.effective_diffusivity_factor(), 1.0 / 16.0, 2e-15);
    close(dry.transport_tortuosity_factor().unwrap(), 2.0, 2e-15);
    let wet = law.sample(0.875, 0.125, 1e-6).unwrap();
    // theta=1/64, theta^(10/3)/phi^2 = 1/16384.
    close(wet.effective_diffusivity_factor(), 1.0 / 16384.0, 2e-14);
    close(wet.transport_tortuosity_factor().unwrap(), 256.0, 2e-14);
    assert_eq!(wet.gas_relative_permeability(), 0.0);
    assert!(wet.effective_diffusivity_factor() > 0.0);
    let flooded = law.sample(1.0, 0.125, 1e-6).unwrap();
    assert_eq!(flooded.effective_diffusivity_factor(), 0.0);
    assert_eq!(flooded.transport_tortuosity_factor(), None);
}

#[test]
fn prescribed_geometry_evolution_preserves_independent_permeability_and_pore_radius() {
    let law = law(GasDiffusivityLaw::MillingtonQuirk1961);
    let a = law.sample(0.3, 0.125, 1e-6).unwrap();
    let b = law.sample(0.3, 1.0, 2e-6).unwrap();
    assert_eq!(a.gas_relative_permeability(), b.gas_relative_permeability());
    close(
        b.effective_diffusivity_factor() / a.effective_diffusivity_factor(),
        16.0,
        5e-15,
    );
    assert_eq!(b.mean_gas_pore_radius_m() / a.mean_gas_pore_radius_m(), 2.0);
    // The DGM Knudsen length factor uses current r and absolute D_eff/D_bulk
    // once. It scales 32-fold here, while the independently supplied skeleton
    // permeability is not inferred from either changed geometry input.
    close(
        (b.mean_gas_pore_radius_m() * b.effective_diffusivity_factor())
            / (a.mean_gas_pore_radius_m() * a.effective_diffusivity_factor()),
        32.0,
        5e-15,
    );
    assert_eq!(
        a,
        law.sample(0.3, 0.125, 1e-6).unwrap(),
        "sampling another state must not mutate calibration or history"
    );
}

#[test]
fn transport_is_bounded_monotone_and_repeatable_across_supported_saturation() {
    for diffusivity in [
        GasDiffusivityLaw::MillingtonQuirk1961,
        GasDiffusivityLaw::ConstantTransportTortuosity { factor: 1.0 },
    ] {
        let law = law(diffusivity);
        let mut previous_kr = 1.0;
        let mut previous_d = 1.0;
        for i in 0..=100 {
            let state = law.sample(f64::from(i) / 100.0, 0.42, 3e-6).unwrap();
            assert!((0.0..=previous_kr).contains(&state.gas_relative_permeability()));
            assert!((0.0..=previous_d).contains(&state.effective_diffusivity_factor()));
            assert!(state.effective_diffusivity_factor() <= state.gas_filled_porosity());
            assert_eq!(state, law.sample(f64::from(i) / 100.0, 0.42, 3e-6).unwrap());
            previous_kr = state.gas_relative_permeability();
            previous_d = state.effective_diffusivity_factor();
        }
    }
}

#[test]
fn adjacent_wet_endpoint_retains_positive_flow_without_cancellation() {
    let law = law(GasDiffusivityLaw::MillingtonQuirk1961);
    let upper = 1.0 - law.residual_gas_saturation();
    let below = f64::from_bits(upper.to_bits() - 1);
    let state = law.sample(below, 0.4, 1e-6).unwrap();
    let gas = (upper - below) / (upper - law.residual_liquid_saturation());
    // m=1/2: 1-Se^2 = gas*(2-gas), evaluated without 1-Se cancellation.
    let expected = gas.sqrt() * gas * (2.0 - gas);
    close(state.gas_relative_permeability(), expected, 5e-14);
    assert!(state.gas_relative_permeability() > 0.0);
    assert_eq!(
        law.sample(upper, 0.4, 1e-6)
            .unwrap()
            .gas_relative_permeability(),
        0.0
    );
    let before_flooding = f64::from_bits(1.0_f64.to_bits() - 1);
    assert!(
        law.sample(before_flooding, 0.4, 1e-6)
            .unwrap()
            .effective_diffusivity_factor()
            > 0.0
    );
}

#[test]
fn calibration_domains_are_enforced_before_endpoint_extension() {
    let mut bounds = support();
    bounds.liquid_saturation = (0.2, 0.7);
    bounds.total_porosity = (0.3, 0.4);
    bounds.mean_gas_pore_radius_m = (1e-6, 2e-6);
    let law = WetPoreGasLaw::try_new(
        0.1,
        0.2,
        0.5,
        GasDiffusivityLaw::MillingtonQuirk1961,
        bounds,
        "bounded specimen fit",
    )
    .unwrap();
    assert!(law.sample(0.2, 0.3, 1e-6).is_ok());
    assert!(law.sample(0.7, 0.4, 2e-6).is_ok());
    for (s, phi, r) in [
        (0.0, 0.35, 1.5e-6),
        (1.0, 0.35, 1.5e-6),
        (0.4, 0.2, 1.5e-6),
        (0.4, 0.5, 1.5e-6),
        (0.4, 0.35, 3e-6),
        (f64::NAN, 0.35, 1.5e-6),
        (0.4, f64::INFINITY, 1.5e-6),
        (0.4, 0.35, f64::NAN),
    ] {
        assert!(
            law.sample(s, phi, r).is_err(),
            "invalid sample {s}, {phi}, {r}"
        );
    }
}

#[test]
fn invalid_calibrations_do_not_become_default_materials() {
    let d = GasDiffusivityLaw::MillingtonQuirk1961;
    for (sl, sg, m) in [
        (-0.1, 0.2, 0.5),
        (0.1, -0.2, 0.5),
        (0.5, 0.5, 0.5),
        (0.8, 0.3, 0.5),
        (0.1, 0.2, 0.0),
        (0.1, 0.2, 1.0),
        (0.1, 0.2, f64::NAN),
    ] {
        assert!(WetPoreGasLaw::try_new(sl, sg, m, d, support(), "fixture").is_err());
    }
    for factor in [0.0, 0.9, f64::NAN, f64::INFINITY] {
        assert!(
            WetPoreGasLaw::try_new(
                0.1,
                0.2,
                0.5,
                GasDiffusivityLaw::ConstantTransportTortuosity { factor },
                support(),
                "fixture"
            )
            .is_err()
        );
    }
    assert!(WetPoreGasLaw::try_new(0.1, 0.2, 0.5, d, support(), " ").is_err());
    let mut invalid_support = support();
    invalid_support.total_porosity = (0.0, 1.0);
    assert!(WetPoreGasLaw::try_new(0.1, 0.2, 0.5, d, invalid_support, "fixture").is_err());
    invalid_support = support();
    invalid_support.mean_gas_pore_radius_m = (2e-6, 1e-6);
    assert!(WetPoreGasLaw::try_new(0.1, 0.2, 0.5, d, invalid_support, "fixture").is_err());
}

#[test]
fn stable_porosity_scaling_avoids_premature_underflow_and_refuses_true_loss() {
    let mut bounds = support();
    bounds.total_porosity = (1e-300, 1.0);
    let mq = WetPoreGasLaw::try_new(
        0.1,
        0.2,
        0.5,
        GasDiffusivityLaw::MillingtonQuirk1961,
        bounds,
        "numerical range fixture",
    )
    .unwrap();
    let tiny = mq.sample(0.0, 1e-180, 1e-6).unwrap();
    close(tiny.effective_diffusivity_factor(), 1e-240, 2e-13);
    // Directly forming theta^(10/3) would already be zero at phi=1e-180.
    assert_eq!(1e-180_f64.powf(10.0 / 3.0), 0.0);
    assert!(
        mq.sample(0.0, 1e-250, 1e-6).is_err(),
        "a positive but unrepresentable transport factor must not become a sealed material"
    );
    let high_tau = WetPoreGasLaw::try_new(
        0.1,
        0.2,
        0.5,
        GasDiffusivityLaw::ConstantTransportTortuosity { factor: 1e308 },
        bounds,
        "numerical range fixture",
    )
    .unwrap();
    assert!(high_tau.sample(0.0, 1e-300, 1e-6).is_err());
    assert_eq!(
        mq.sample(1.0, 1e-300, 1e-6)
            .unwrap()
            .effective_diffusivity_factor(),
        0.0
    );
}
