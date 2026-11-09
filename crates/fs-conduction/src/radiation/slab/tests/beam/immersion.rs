use super::*;

fn immersed(depths: &[f64], indices: &[f64], exterior: f64, mu: f64, s: f64) -> SlabBeamResponse {
    SlabBeamResponse::new_in_medium(
        layers(depths, indices),
        exterior,
        "declared surrounding dielectric",
        mu,
        s,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_immersed_beams_match_polarized_fresnel_and_absorbing_slab_series() {
    // Independent single-layer geometric series, including illumination from
    // a higher-index medium below the critical angle. No native slab helper
    // supplies the reference Fresnel factors or internal attenuation.
    for (exterior, index, mu) in [(1.33_f64, 1.5_f64, 0.4_f64), (1.5, 1.2, 0.8)] {
        let transmitted_cosine = (1.0 - (exterior / index).powi(2) * (1.0 - mu * mu)).sqrt();
        let rs = ((exterior * mu - index * transmitted_cosine)
            / (exterior * mu + index * transmitted_cosine))
            .powi(2);
        let rp = ((index * mu - exterior * transmitted_cosine)
            / (index * mu + exterior * transmitted_cosine))
            .powi(2);
        let t = (-0.3 / transmitted_cosine).exp();
        for s in [0.0, 0.37, 1.0] {
            let response = immersed(&[0.3], &[index], exterior, mu, s);
            let mut expected = [0.0; 3];
            for (weight, r) in [(s, rs), (1.0 - s, rp)] {
                let denominator = 1.0 - r * r * t * t;
                expected[0] += weight * (r + (1.0 - r).powi(2) * t * t * r / denominator);
                expected[1] += weight * (1.0 - r).powi(2) * t / denominator;
                expected[2] += weight * (1.0 - r) * (1.0 - t) * (1.0 + r * t) / denominator;
            }
            for side in 0..2 {
                close(response.reflectances()[side], expected[0], 3e-14);
                close(response.layer_absorptance()[0][side], expected[2], 3e-14);
            }
            close(response.transmittance(), expected[1], 3e-14);
            close(response.layer_direction_cosine(0).unwrap(), transmitted_cosine, 2e-14);
            let mut heat = [0.0];
            let report = response.exchange([7.0, 11.0], &mut heat).unwrap();
            close(heat[0], 18.0 * expected[2], 3e-14);
            close(
                response.layer_absorbed_flux(0, [0.0, 1.0], [7.0, 11.0]).unwrap(),
                heat[0],
                3e-14,
            );
            assert!(report.balance_residual_w_m2.abs() < 1e-13);
        }
    }
}

#[test]
fn g1_total_reflection_has_no_invented_gap_transmission_or_absorption() {
    for absorption in [0.0, 0.3, 1e100] {
        let response = immersed(&[absorption], &[1.0], 1.5, 0.5, 0.3);
        assert_eq!(response.reflectances(), [1.0; 2]);
        assert_eq!(response.transmittance(), 0.0);
        assert_eq!(response.layer_absorptance(), &[[0.0; 2]]);
        assert_eq!(response.layer_direction_cosine(0), Some(0.0));
        let mut heat = [99.0];
        let result = response.exchange([7.0, 11.0], &mut heat).unwrap();
        assert_eq!(result.reflected_w_m2, [7.0, 11.0]);
        assert_eq!(result.balance_residual_w_m2, 0.0);
        assert_eq!(heat, [0.0]);
        for depth in [0.0, 0.5, 1.0] {
            assert_eq!(
                response.layer_fluxes_at_depth(0, depth, [7.0, 11.0]).unwrap(),
                [0.0; 2],
            );
        }
        assert_eq!(
            response.layer_absorbed_flux(0, [0.0, 1.0], [7.0, 11.0]).unwrap(),
            0.0,
        );
    }

    // A low-index interior gap separates two illuminated absorbing regions.
    // Each accessible layer has a matched exterior and one perfect return:
    // R=exp(-2 alpha d/mu), A=1-R, with no coupling across the blocked gap.
    let response = immersed(&[0.2, 0.0, 0.4], &[1.5, 1.0, 1.5], 1.5, 0.5, 0.5);
    assert_eq!(response.transmittance(), 0.0);
    let mut heat = [0.0; 3];
    response.exchange([7.0, 11.0], &mut heat).unwrap();
    close(heat[0], 7.0 * -(-0.8_f64).exp_m1(), 3e-14);
    assert_eq!(heat[1], 0.0);
    close(heat[2], 11.0 * -(-1.6_f64).exp_m1(), 3e-14);
    assert_eq!(response.layer_absorptance()[0][1], 0.0);
    assert_eq!(response.layer_absorptance()[2][0], 0.0);
    for (side, depth) in [0.2, 0.4].into_iter().enumerate() {
        close(response.reflectances()[side], (-4.0_f64 * depth).exp(), 3e-14);
    }
}

#[test]
fn g3_immersed_power_obeys_index_scaling_reversal_and_layer_subdivision() {
    let original = immersed(&[0.2, 0.9], &[1.4, 2.1], 1.33, 0.4, 0.7);
    let scaled = immersed(&[0.2, 0.9], &[2.8, 4.2], 2.66, 0.4, 0.7);
    let reversed = immersed(&[0.9, 0.2], &[2.1, 1.4], 1.33, 0.4, 0.7);
    let split = immersed(&[0.1, 0.1, 0.9], &[1.4, 1.4, 2.1], 1.33, 0.4, 0.7);
    for response in [&scaled, &reversed, &split] {
        close(response.transmittance(), original.transmittance(), 5e-14);
    }
    for side in 0..2 {
        close(original.reflectances()[side], scaled.reflectances()[side], 3e-14);
        close(original.reflectances()[side], reversed.reflectances()[1 - side], 3e-14);
        close(original.reflectances()[side], split.reflectances()[side], 3e-14);
        close(
            original.layer_absorptance()[0][side],
            split.layer_absorptance()[0][side] + split.layer_absorptance()[1][side],
            3e-14,
        );
    }
    // A fully matched immersed slab has Beer attenuation with no n² multiplier
    // on the already projected physical beam flux.
    let matched = immersed(&[0.2], &[1.5], 1.5, 0.5, 0.37);
    assert_eq!(matched.reflectances(), [0.0; 2]);
    close(matched.transmittance(), (-0.4_f64).exp(), 2e-14);
    assert_eq!(original.exterior_refractive_index(), 1.33);
    assert_eq!(original.exterior_index_source(), Some("declared surrounding dielectric"));
    assert_eq!(beam(&[0.2], &[1.5], 0.5, 0.5).exterior_index_source(), None);
}

#[test]
fn g0_immersed_critical_angle_and_environment_admission_are_explicit() {
    let critical_cosine = (1.0_f64 - (1.0_f64 / 1.5).powi(2)).sqrt();
    let below = immersed(&[0.0], &[1.0], 1.5, critical_cosine + 1e-7, 0.5);
    let above = immersed(&[0.0], &[1.0], 1.5, critical_cosine - 1e-7, 0.5);
    assert!(below.transmittance() > 0.0);
    assert!(below.layer_direction_cosine(0).unwrap() > 0.0);
    assert_eq!(above.transmittance(), 0.0);
    assert_eq!(above.reflectances(), [1.0; 2]);
    for index in [0.0, 0.9, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(SlabBeamResponse::new_in_medium(
            layers(&[0.2], &[1.5]), index, "source", 0.5, 0.5, || false,
        ).is_err());
    }
    assert!(SlabBeamResponse::new_in_medium(
        layers(&[0.2], &[1.5]), 1.33, " ", 0.5, 0.5, || false,
    ).is_err());
    assert!(matches!(
        SlabBeamResponse::new_in_medium(
            layers(&[0.2], &[1.5]), 1.33, "source", 0.5, 0.5, || true,
        ),
        Err(ConductionError::Cancelled { .. }),
    ));

    // The original diffuse thermal owner assumes index-one exterior radiance.
    // An immersed beam must not silently be combined with that other medium.
    let slab = AbsorbingSlab::new(
        layers(&[0.2], &[1.5]), SlabQuadrature::default(), || false,
    ).unwrap();
    let response = immersed(&[0.2], &[1.5], 1.33, 0.5, 0.5);
    let mut heat = [99.0];
    assert!(slab.layer_net_fluxes_with_beam(
        [0.0; 2], &[500.0], &response, [7.0, 0.0], &mut heat,
    ).is_err());
    assert_eq!(heat, [99.0]);
}
