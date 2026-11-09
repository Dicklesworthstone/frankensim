use super::*;

fn immersed(depths: &[f64], indices: &[f64], exterior: f64) -> AbsorbingSlab {
    AbsorbingSlab::new_in_medium(
        depths
            .iter()
            .zip(indices)
            .map(|(&depth, &index)| {
                AbsorbingLayer::new(1.0, depth, [100.0, 3000.0], "absorbing material")
                    .unwrap()
                    .with_refractive_index(index, "absolute layer index")
                    .unwrap()
            })
            .collect(),
        exterior,
        "declared exterior dielectric",
        SlabQuadrature::default(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_matched_immersed_slab_retains_beer_exchange_and_absolute_mode_density() {
    let vacuum = slab(&[0.2, 0.7]);
    let medium = immersed(&[0.2, 0.7], &[1.5, 1.5], 1.5);
    close(medium.transmittance(), vacuum.transmittance(), 3e-10);
    assert_eq!(medium.reflectances(), [0.0; 2]);
    for (a, b) in medium
        .face_coupling()
        .iter()
        .flatten()
        .zip(vacuum.face_coupling().iter().flatten())
    {
        close(*a, *b, 3e-10);
    }
    for (a, b) in medium.layer_coupling().iter().zip(vacuum.layer_coupling()) {
        close(*a, 2.25 * b, 3e-10);
    }
    let b = [700.0, 500.0];
    let emission = vacuum.emission_to_faces(&b).unwrap();
    let immersed_emission = medium.emission_to_faces(&b).unwrap();
    for side in 0..2 {
        close(immersed_emission[side], 2.25 * emission[side], 3e-10);
    }
    let mut vacuum_heat = [0.0; 2];
    let mut medium_heat = [0.0; 2];
    vacuum
        .layer_net_fluxes([100.0, 900.0], &b, &mut vacuum_heat)
        .unwrap();
    medium
        .layer_net_fluxes([225.0, 2025.0], &b, &mut medium_heat)
        .unwrap();
    for i in 0..2 {
        close(medium_heat[i], 2.25 * vacuum_heat[i], 3e-10);
    }
}

#[test]
fn g1_escape_cone_matches_independent_fresnel_reciprocity_integral() {
    // An opaque low-index layer viewed through a higher-index exterior.
    // Integrate from the opposite, low-index side, where there is no critical
    // singularity. Etendue reciprocity gives A_ext = (n_layer/n_ext)^2 A_layer.
    let exterior = 1.5_f64;
    let index = 1.0_f64;
    let panels = 16_384;
    let mut escape = 0.0;
    for k in 1..=panels {
        let mu = k as f64 / panels as f64;
        let transmitted = (1.0 - (index / exterior).powi(2) * (1.0 - mu * mu)).sqrt();
        let rs =
            ((index * mu - exterior * transmitted) / (index * mu + exterior * transmitted)).powi(2);
        let rp =
            ((exterior * mu - index * transmitted) / (exterior * mu + index * transmitted)).powi(2);
        let weight = if k == panels {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        escape += weight * 2.0 * mu * (1.0 - 0.5 * (rs + rp)) / (3.0 * panels as f64);
    }
    let material = immersed(&[1000.0], &[index], exterior);
    let expected_absorption = escape * (index / exterior).powi(2);
    assert_eq!(material.transmittance(), 0.0);
    for side in 0..2 {
        close(
            material.face_coupling()[0][side],
            expected_absorption,
            2e-10,
        );
        close(
            material.reflectances()[side],
            1.0 - expected_absorption,
            2e-10,
        );
        close(
            material.emission_to_faces(&[700.0]).unwrap()[side],
            700.0 * escape,
            2e-10,
        );
    }
    assert!(material.reflectance() > 0.5); // The exterior TIR cone is retained.
}

#[test]
fn g0_immersed_equilibrium_energy_balance_and_reciprocity_close() {
    for exterior in [1.0_f64, 1.33, 1.7, 2.5] {
        let material = immersed(&[0.2, 0.0, 0.8], &[1.4, 1.0, 2.1], exterior);
        let equilibrium = exterior.powi(2) * 500.0;
        let mut heat = [99.0; 3];
        material
            .layer_net_fluxes([equilibrium; 2], &[500.0; 3], &mut heat)
            .unwrap();
        assert_eq!(heat, [0.0; 3]);
        let emission = material.emission_to_faces(&[500.0; 3]).unwrap();
        for side in 0..2 {
            close(
                emission[side]
                    + equilibrium * (material.reflectances()[side] + material.transmittance()),
                equilibrium,
                2e-12,
            );
        }
        let b = [750.0, 500.0, 350.0];
        let h = [150.0, 800.0];
        material.layer_net_fluxes(h, &b, &mut heat).unwrap();
        let emitted = material.emission_to_faces(&b).unwrap();
        let outgoing: f64 = (0..2)
            .map(|side| {
                emitted[side]
                    + material.reflectances()[side] * h[side]
                    + material.transmittance() * h[1 - side]
            })
            .sum();
        close(heat.iter().sum(), h.iter().sum::<f64>() - outgoing, 2e-11);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(
                    material.layer_coupling()[3 * i + j],
                    material.layer_coupling()[3 * j + i]
                );
            }
        }
    }
}

#[test]
fn g3_immersed_scaling_reversal_and_depth_queries_share_one_boundary() {
    let material = immersed(&[0.2, 0.7], &[1.2, 2.1], 1.7);
    let scaled = immersed(&[0.2, 0.7], &[2.4, 4.2], 3.4);
    let reversed = immersed(&[0.7, 0.2], &[2.1, 1.2], 1.7);
    close(material.transmittance(), scaled.transmittance(), 2e-10);
    close(material.transmittance(), reversed.transmittance(), 2e-10);
    for side in 0..2 {
        close(
            material.reflectances()[side],
            scaled.reflectances()[side],
            2e-10,
        );
        close(
            material.reflectances()[side],
            reversed.reflectances()[1 - side],
            2e-10,
        );
        for i in 0..2 {
            close(
                material.face_coupling()[i][side],
                scaled.face_coupling()[i][side],
                2e-10,
            );
            close(
                material.face_coupling()[i][side],
                reversed.face_coupling()[1 - i][1 - side],
                2e-10,
            );
        }
    }
    for (a, b) in scaled
        .layer_coupling()
        .iter()
        .zip(material.layer_coupling())
    {
        close(*a, 4.0 * b, 3e-10);
    }
    for i in 0..2 {
        let full = material
            .diffuse_depth_response(i, [0.0, 1.0], SlabQuadrature::default(), || false)
            .unwrap();
        let front = material
            .diffuse_depth_response(i, [0.0, 0.3], SlabQuadrature::default(), || false)
            .unwrap();
        let back = material
            .diffuse_depth_response(i, [0.3, 1.0], SlabQuadrature::default(), || false)
            .unwrap();
        for side in 0..2 {
            close(
                full.absorptance()[side],
                material.face_coupling()[i][side],
                2e-8,
            );
            close(
                front.absorptance()[side] + back.absorptance()[side],
                full.absorptance()[side],
                2e-8,
            );
        }
        assert_eq!(full.exterior_refractive_index(), 1.7);
        assert_eq!(
            full.exterior_index_source(),
            Some("declared exterior dielectric")
        );
    }
    let beam = material.beam_response(0.8, 0.3, || false).unwrap();
    assert_eq!(
        beam.exterior_refractive_index(),
        material.exterior_refractive_index()
    );
    assert_eq!(
        beam.exterior_index_source(),
        material.exterior_index_source()
    );
    let mut beam_heat = [0.0; 2];
    beam.exchange([7.0, 11.0], &mut beam_heat).unwrap();
    let mut combined_heat = [0.0; 2];
    material
        .layer_net_fluxes_with_beam(
            [1.7_f64.powi(2) * 500.0; 2],
            &[500.0; 2],
            &beam,
            [7.0, 11.0],
            &mut combined_heat,
        )
        .unwrap();
    assert_eq!(combined_heat, beam_heat);
    let wrong_boundary =
        SlabBeamResponse::new(material.layers().to_vec(), 0.8, 0.3, || false).unwrap();
    let mut unchanged = [99.0; 2];
    assert!(
        material
            .layer_net_fluxes_with_beam(
                [0.0; 2],
                &[500.0; 2],
                &wrong_boundary,
                [7.0, 11.0],
                &mut unchanged
            )
            .is_err()
    );
    assert_eq!(unchanged, [99.0; 2]);
}

#[test]
fn g0_immersed_admission_and_exact_transparent_limit() {
    let layer = AbsorbingLayer::new(1.0, 0.0, [100.0, 3000.0], "clear")
        .unwrap()
        .with_refractive_index(1.5, "index")
        .unwrap();
    let controls = SlabQuadrature {
        maximum_evaluations: 0,
        ..SlabQuadrature::default()
    };
    let clear =
        AbsorbingSlab::new_in_medium(vec![layer.clone()], 1.5, "medium", controls, || false)
            .unwrap();
    assert_eq!(clear.reflectances(), [0.0; 2]);
    assert_eq!(clear.transmittance(), 1.0);
    assert_eq!(clear.quadrature().angular_evaluations, 0);
    assert_eq!(clear.exterior_index_source(), Some("medium"));
    for exterior in [0.0, 0.9, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            AbsorbingSlab::new_in_medium(vec![layer.clone()], exterior, "medium", controls, || {
                false
            })
            .is_err()
        );
    }
    assert!(
        AbsorbingSlab::new_in_medium(vec![layer.clone()], 1.5, " ", controls, || false).is_err()
    );
    assert!(matches!(
        AbsorbingSlab::new_in_medium(vec![layer], 1.5, "medium", controls, || true),
        Err(ConductionError::Cancelled { .. })
    ));
    // Finite escape can remain representable when hypothetical unweighted
    // exterior exitance is too large; exact zero coupling stays exact.
    close(
        exterior_difference_product(1e-300, 0.0, 1e300, 1e100).unwrap(),
        -1e100,
        2e-15,
    );
    assert_eq!(
        exterior_difference_product(0.0, 0.0, 1e300, 1e100).unwrap(),
        0.0
    );
}

#[test]
fn g3_immersed_net_depth_heat_preserves_whole_layer_and_cell_balance() {
    for exterior in [1.33_f64, 1.7, 2.5] {
        let material = immersed(&[0.2, 0.0, 0.8], &[1.4, 1.0, 2.1], exterior);
        let h = [150.0, 800.0];
        let b = [750.0, 500.0, 350.0];
        let mut whole = [0.0; 3];
        material.layer_net_fluxes(h, &b, &mut whole).unwrap();
        for layer in 0..3 {
            let cell = |interval| {
                material
                    .depth_response(layer, interval, SlabQuadrature::default(), || false)
                    .unwrap()
            };
            let full = cell([0.0, 1.0]);
            close(full.net_flux(h, &b).unwrap(), whole[layer], 3e-8);
            close(
                cell([0.0, 0.3]).net_flux(h, &b).unwrap()
                    + cell([0.3, 1.0]).net_flux(h, &b).unwrap(),
                whole[layer],
                3e-8,
            );
            assert_eq!(
                full.net_flux([exterior.powi(2) * 500.0; 2], &[500.0; 3])
                    .unwrap(),
                0.0
            );
            assert_eq!(full.exterior_refractive_index(), exterior);
            assert_eq!(
                full.exterior_index_source(),
                material.exterior_index_source()
            );
            for side in 0..2 {
                close(
                    full.face_coupling()[side],
                    material.face_coupling()[layer][side],
                    2e-10,
                );
            }
            for source in 0..3 {
                close(
                    full.layer_coupling()[source],
                    material.layer_coupling()[3 * layer + source],
                    2e-10,
                );
            }
        }
    }
}
