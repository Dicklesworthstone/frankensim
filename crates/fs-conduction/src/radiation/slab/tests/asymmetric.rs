use super::*;

pub(super) fn between(depths: &[f64], indices: &[f64], exterior: [f64; 2]) -> AbsorbingSlab {
    AbsorbingSlab::new_between_media(
        depths
            .iter()
            .zip(indices)
            .map(|(&depth, &index)| {
                AbsorbingLayer::new(1.0, depth, [100.0, 3000.0], "bulk absorption")
                    .unwrap()
                    .with_refractive_index(index, "absolute layer index")
                    .unwrap()
            })
            .collect(),
        exterior,
        ["front optical source", "back optical source"],
        SlabQuadrature::default(),
        || false,
    )
    .unwrap()
}

// Independent one-layer intensity series, with ordinary Fresnel amplitudes.
// This does not call the production reflection elimination or its helpers.
fn one_layer(
    index: f64,
    depth: f64,
    exterior: [f64; 2],
    side: usize,
    mu: f64,
    s: f64,
) -> ([f64; 2], f64, [f64; 2], [f64; 2]) {
    let p = exterior[side] * (1.0 - mu * mu).sqrt();
    let cosine = |n: f64| {
        if p < n {
            (1.0 - (p / n).powi(2)).sqrt()
        } else {
            0.0
        }
    };
    let outside = exterior.map(cosine);
    let inside = cosine(index);
    assert!(inside > 0.0);
    let t = (-depth / inside).exp();
    let mut reflected = [0.0; 2];
    let mut transmitted = 0.0;
    let mut absorbed = [0.0; 2];
    for (pol, weight) in [s, 1.0 - s].into_iter().enumerate() {
        let r: [f64; 2] = core::array::from_fn(|face| {
            if outside[face] == 0.0 {
                return 1.0;
            }
            let (a, b) = if pol == 0 {
                (exterior[face] * outside[face], index * inside)
            } else {
                (exterior[face] * inside, index * outside[face])
            };
            ((a - b) / (a + b)).powi(2)
        });
        let e = r.map(|r| 1.0 - r);
        let denominator = 1.0 - r[0] * r[1] * t * t;
        transmitted += weight * e[0] * e[1] * t / denominator;
        for face in 0..2 {
            if outside[face] > 0.0 {
                reflected[face] +=
                    weight * (r[face] + e[face] * e[face] * t * t * r[1 - face] / denominator);
                absorbed[face] +=
                    weight * e[face] * (1.0 - t) * (1.0 + r[1 - face] * t) / denominator;
            }
        }
    }
    (reflected, transmitted, absorbed, outside)
}

#[test]
fn g1_two_exterior_beams_match_independent_polarized_return_series() {
    for exterior in [[1.0, 1.33], [1.7, 1.1], [1.33, 2.0]] {
        let slab = between(&[0.37], &[2.1], exterior);
        for side in 0..2 {
            for mu in [0.2, 0.65, 1.0] {
                for s in [0.0, 0.4, 1.0] {
                    let beam = slab.beam_response_from_side(side, mu, s, || false).unwrap();
                    let (r, t, a, cosines) = one_layer(2.1, 0.37, exterior, side, mu, s);
                    assert_eq!(beam.incident_side(), side);
                    assert_eq!(beam.direction_cosine(), mu);
                    assert_eq!(beam.exterior_refractive_indices(), exterior);
                    close(beam.transmittance(), t, 5e-14);
                    for face in 0..2 {
                        close(beam.reflectances()[face], r[face], 5e-14);
                        close(beam.layer_absorptance()[0][face], a[face], 5e-14);
                        close(
                            beam.exterior_direction_cosines()[face],
                            cosines[face],
                            5e-14,
                        );
                    }
                    if cosines[1 - side] > 0.0 {
                        let reverse = slab
                            .beam_response_from_side(1 - side, cosines[1 - side], s, || false)
                            .unwrap();
                        close(reverse.transmittance(), beam.transmittance(), 2e-13);
                    }
                }
            }
        }
    }
}

#[test]
fn g1_diffuse_transmission_uses_each_exterior_mode_measure() {
    // A transparent layer matched to the incident medium leaves ONE interface.
    let slab = between(&[0.0], &[1.0], [1.0, 1.5]);
    let panels = 32_768;
    let mut expected_r = 0.0;
    for i in 0..panels {
        let mu = (i as f64 + 0.5) / panels as f64;
        let transmitted_mu = (1.0 - (1.0 - mu * mu) / 2.25).sqrt();
        let rs = ((mu - 1.5 * transmitted_mu) / (mu + 1.5 * transmitted_mu)).powi(2);
        let rp = ((1.5 * mu - transmitted_mu) / (1.5 * mu + transmitted_mu)).powi(2);
        expected_r += mu * (rs + rp) / panels as f64;
    }
    close(slab.reflectances()[0], expected_r, 2e-9);
    close(slab.transmittances()[0], 1.0 - expected_r, 2e-9);
    close(slab.transmittances()[1], (1.0 - expected_r) / 2.25, 2e-9);
    close(
        slab.reflectances()[1],
        1.0 - slab.transmittances()[1],
        3e-14,
    );
    assert_eq!(slab.face_coupling(), &[[0.0; 2]]);
    assert!(slab.reflectances()[1] > slab.reflectances()[0] + 0.4);
    assert_eq!(
        slab.exterior_index_sources(),
        [Some("front optical source"), Some("back optical source")]
    );
}

#[test]
fn g0_asymmetric_thermal_equilibrium_and_complete_energy_balance() {
    for exterior in [[1.0, 1.7], [2.5, 1.33], [1.1, 1.8]] {
        let slab = between(&[0.2, 0.9, 0.4], &[1.4, 2.1, 1.2], exterior);
        let mut heat = [0.0; 3];
        let bath = exterior.map(|index| index * index * 7.0);
        slab.layer_net_fluxes(bath, &[7.0; 3], &mut heat).unwrap();
        assert_eq!(heat, [0.0; 3]);
        let outgoing = slab.outgoing_fluxes(bath, &[7.0; 3]).unwrap();
        for side in 0..2 {
            close(outgoing[side], bath[side], 3e-12);
        }
        let h = [11.0, 5.0];
        let b = [3.0, 7.0, 2.0];
        slab.layer_net_fluxes(h, &b, &mut heat).unwrap();
        let outgoing = slab.outgoing_fluxes(h, &b).unwrap();
        close(
            heat.iter().sum::<f64>() + outgoing.iter().sum::<f64>(),
            16.0,
            3e-12,
        );
        let t = slab.transmittances();
        close(
            exterior[0].powi(2) * t[0],
            exterior[1].powi(2) * t[1],
            3e-12,
        );
        assert!(slab.quadrature().maximum_partition_residual < 2e-13);
    }
}

#[test]
fn g3_asymmetric_reversal_subdivision_and_index_scaling_preserve_heat() {
    let slab = between(&[0.2, 0.9, 0.4], &[1.4, 2.1, 1.2], [1.0, 1.7]);
    let reverse = between(&[0.4, 0.9, 0.2], &[1.2, 2.1, 1.4], [1.7, 1.0]);
    let split = between(&[0.2, 0.3, 0.6, 0.4], &[1.4, 2.1, 2.1, 1.2], [1.0, 1.7]);
    let scaled = between(&[0.2, 0.9, 0.4], &[2.8, 4.2, 2.4], [2.0, 3.4]);
    let mut heat = [0.0; 3];
    let mut back = [0.0; 3];
    let mut divided = [0.0; 4];
    let mut larger = [0.0; 3];
    slab.layer_net_fluxes([11.0, 5.0], &[3.0, 7.0, 2.0], &mut heat)
        .unwrap();
    reverse
        .layer_net_fluxes([5.0, 11.0], &[2.0, 7.0, 3.0], &mut back)
        .unwrap();
    split
        .layer_net_fluxes([11.0, 5.0], &[3.0, 7.0, 7.0, 2.0], &mut divided)
        .unwrap();
    scaled
        .layer_net_fluxes([44.0, 20.0], &[3.0, 7.0, 2.0], &mut larger)
        .unwrap();
    for i in 0..3 {
        close(heat[i], back[2 - i], 3e-11);
        close(4.0 * heat[i], larger[i], 3e-11);
    }
    close(divided[0], heat[0], 3e-11);
    close(divided[1] + divided[2], heat[1], 3e-11);
    close(divided[3], heat[2], 3e-11);
    for side in 0..2 {
        close(
            slab.transmittances()[side],
            reverse.transmittances()[1 - side],
            3e-12,
        );
        close(
            slab.transmittances()[side],
            scaled.transmittances()[side],
            3e-12,
        );
        close(
            slab.reflectances()[side],
            reverse.reflectances()[1 - side],
            3e-12,
        );
    }
}

#[test]
fn g3_asymmetric_depth_queries_retain_local_and_whole_layer_exchange() {
    let slab = between(&[0.2, 0.9, 0.4], &[1.4, 2.1, 1.2], [1.0, 1.7]);
    let h = [11.0, 5.0];
    let b = [3.0, 7.0, 2.0];
    let mut heat = [0.0; 3];
    slab.layer_net_fluxes(h, &b, &mut heat).unwrap();
    for layer in 0..3 {
        let mut net = 0.0;
        let mut absorbed = [0.0; 2];
        for interval in [[0.0, 0.25], [0.25, 0.8], [0.8, 1.0]] {
            let depth = slab
                .depth_response(layer, interval, SlabQuadrature::default(), || false)
                .unwrap();
            let diffuse = slab
                .diffuse_depth_response(layer, interval, SlabQuadrature::default(), || false)
                .unwrap();
            assert_eq!(depth.exterior_refractive_indices(), [1.0, 1.7]);
            assert_eq!(diffuse.exterior_refractive_indices(), [1.0, 1.7]);
            net += depth.net_flux(h, &b).unwrap();
            for side in 0..2 {
                absorbed[side] += diffuse.absorptance()[side];
            }
            assert_eq!(
                depth.net_flux([7.0, 1.7 * 1.7 * 7.0], &[7.0; 3]).unwrap(),
                0.0
            );
        }
        close(net, heat[layer], 4e-11);
        for side in 0..2 {
            close(absorbed[side], slab.face_coupling()[layer][side], 4e-11);
        }
    }
}

#[test]
fn g0_closed_beam_ports_and_asymmetric_source_admission_are_explicit() {
    let slab = between(&[0.2], &[1.5], [1.5, 1.0]);
    let beam = slab.beam_response_from_side(0, 0.5, 0.5, || false).unwrap();
    assert_eq!(beam.exterior_direction_cosines(), [0.5, 0.0]);
    assert_eq!(beam.transmittance(), 0.0);
    close(beam.reflectances()[0], (-0.8_f64).exp(), 3e-14);
    let mut heat = [7.0];
    let exchange = beam.exchange([2.0, 0.0], &mut heat).unwrap();
    close(exchange.balance_residual_w_m2, 0.0, 3e-14);
    heat[0] = 7.0;
    assert!(beam.exchange([2.0, 1.0], &mut heat).is_err());
    assert_eq!(heat, [7.0]);
    assert!(beam.layer_absorbed_flux(0, [0.0, 1.0], [0.0, 1.0]).is_err());
    assert!(slab.beam_response_from_side(2, 0.5, 0.5, || false).is_err());
    for exterior in [[0.9, 1.5], [1.0, f64::INFINITY], [f64::NAN, 1.5]] {
        assert!(
            AbsorbingSlab::new_between_media(
                slab.layers().to_vec(),
                exterior,
                ["front", "back"],
                SlabQuadrature::default(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        AbsorbingSlab::new_between_media(
            slab.layers().to_vec(),
            [1.5, 1.0],
            ["front", " "],
            SlabQuadrature::default(),
            || false
        )
        .is_err()
    );
    assert!(
        AbsorbingSlab::new_between_media(
            slab.layers().to_vec(),
            [1.5, 1.0],
            ["front", "back"],
            SlabQuadrature::default(),
            || true
        )
        .is_err()
    );
    let exhausted = SlabQuadrature {
        maximum_evaluations: 0,
        ..SlabQuadrature::default()
    };
    assert!(
        AbsorbingSlab::new_between_media(
            slab.layers().to_vec(),
            [1.5, 1.0],
            ["front", "back"],
            exhausted,
            || false
        )
        .is_err()
    );
    let wrong = SlabBeamResponse::new_in_medium(
        slab.layers().to_vec(),
        1.5,
        "different exterior",
        0.5,
        0.5,
        || false,
    )
    .unwrap();
    assert!(
        slab.layer_net_fluxes_with_beam([0.0; 2], &[0.0], &wrong, [1.0, 0.0], &mut heat)
            .is_err()
    );
    assert_eq!(heat, [7.0]);
}
