use super::*;

fn cell(slab: &AbsorbingSlab, layer: usize, interval: [f64; 2]) -> SlabDepthResponse {
    slab.depth_response(layer, interval, SlabQuadrature::default(), || false)
        .unwrap()
}

fn indexed(parts: &[(f64, f64, f64)]) -> AbsorbingSlab {
    AbsorbingSlab::new(
        parts
            .iter()
            .map(|&(d, a, n)| {
                AbsorbingLayer::new(d, a, [100.0, 3000.0], "original absorption")
                    .unwrap()
                    .with_refractive_index(n, "original index")
                    .unwrap()
            })
            .collect(),
        SlabQuadrature::default(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_net_depth_thermal_emission_matches_e3_and_exact_equilibrium() {
    let slab = slab(&[1.0, 1.0]);
    let target = cell(&slab, 1, [0.0, 0.25]);
    // Independent 70-digit E1 series/E3 recurrence:
    // 2[E3(0)-E3(1)-E3(.25)+E3(1.25)]. The source is the adjacent
    // whole layer; the receiver is only the first quarter of the target.
    close(target.layer_coupling()[0], 0.288_392_509_004_879_56, 2e-11);
    assert_eq!(target.layer_coupling()[1], 0.0);
    let external = slab
        .diffuse_depth_response(1, [0.0, 0.25], SlabQuadrature::default(), || false)
        .unwrap();
    for side in 0..2 {
        close(
            target.face_coupling()[side],
            external.absorptance()[side],
            2e-11,
        );
    }
    close(
        target.net_flux([100.0; 2], &[200.0, 100.0]).unwrap(),
        100.0 * 0.288_392_509_004_879_56,
        2e-9,
    );
    for background in [0.0, 700.0, f64::MAX, -f64::MAX] {
        assert_eq!(
            target.net_flux([background; 2], &[background; 2]).unwrap(),
            0.0
        );
    }
    // Signed deviations have the same result as their common-background form.
    close(
        target.net_flux([30.0, -20.0], &[15.0, -5.0]).unwrap(),
        target.net_flux([1030.0, 980.0], &[1015.0, 995.0]).unwrap(),
        1e-13,
    );
    assert_eq!(target.layers().as_ptr(), slab.layers().as_ptr());
    assert_eq!(target.layer_index(), 1);
    assert_eq!(target.interval_m(), [0.0, 0.25]);
    assert!(target.quadrature().accepted_panels > 0);
}

#[test]
fn g3_net_depth_full_layers_and_additive_cells_match_original_heat() {
    for slab in [
        slab(&[0.2, 0.9, 0.4]),
        indexed(&[(1.0, 0.2, 1.4), (1.0, 0.9, 2.1), (1.0, 0.4, 1.2)]),
    ] {
        let h = [1100.0, 230.0];
        let b = [550.0, 870.0, 310.0];
        let mut original = [0.0; 3];
        slab.layer_net_fluxes(h, &b, &mut original).unwrap();
        for target in 0..3 {
            let whole = cell(&slab, target, [0.0, 1.0]);
            for side in 0..2 {
                close(
                    whole.face_coupling()[side],
                    slab.face_coupling()[target][side],
                    2e-10,
                );
            }
            for source in 0..3 {
                close(
                    whole.layer_coupling()[source],
                    slab.layer_coupling()[target * 3 + source],
                    2e-10,
                );
            }
            let mut sum = 0.0;
            for interval in [[0.0, 0.04], [0.04, 0.35], [0.35, 0.8], [0.8, 1.0]] {
                sum += cell(&slab, target, interval).net_flux(h, &b).unwrap();
            }
            close(whole.net_flux(h, &b).unwrap(), original[target], 3e-7);
            close(sum, original[target], 3e-7);
        }
    }
}

#[test]
fn g3_net_depth_reciprocal_cells_close_heat_without_new_interfaces() {
    let slab = indexed(&[(1.0, 0.7, 2.1), (1.0, 1.2, 1.5)]);
    let intervals = [
        [0.0, 0.2],
        [0.2, 0.7],
        [0.7, 1.0],
        [0.0, 0.3],
        [0.3, 0.8],
        [0.8, 1.0],
    ];
    // Spatial refinement is an independent oracle only. Equal adjacent
    // indices introduce no Fresnel interface; the production responses below
    // all retain the original two-layer stack and source identities.
    let refined = indexed(&[
        (0.2, 0.7, 2.1),
        (0.5, 0.7, 2.1),
        (0.3, 0.7, 2.1),
        (0.3, 1.2, 1.5),
        (0.5, 1.2, 1.5),
        (0.2, 1.2, 1.5),
    ]);
    let h = [190.0, 30.0];
    let b = [120.0, 450.0];
    let mut oracle = [0.0; 6];
    refined
        .layer_net_fluxes(h, &[120.0, 120.0, 120.0, 450.0, 450.0, 450.0], &mut oracle)
        .unwrap();
    let mut absorbed = 0.0;
    for (i, interval) in intervals.into_iter().enumerate() {
        let response = cell(&slab, i / 3, interval);
        let heat = response.net_flux(h, &b).unwrap();
        close(heat, oracle[i], 2e-7);
        absorbed += heat;
    }
    // Independently prepared directions, not a mirrored matrix assignment.
    let a = cell(&refined, 1, [0.0, 0.5]);
    let b_cell = cell(&refined, 4, [0.0, 0.5]);
    close(a.layer_coupling()[4], b_cell.layer_coupling()[1], 2e-11);
    assert!(a.layer_coupling()[4] > 0.0);
    let emitted = slab.emission_to_faces(&b).unwrap();
    let outgoing: f64 = (0..2)
        .map(|side| {
            emitted[side] + slab.reflectances()[side] * h[side] + slab.transmittance() * h[1 - side]
        })
        .sum();
    close(absorbed, h[0] + h[1] - outgoing, 5e-7);
}

#[test]
fn g1_net_depth_includes_trapped_n_squared_and_weak_interval_exchange() {
    let slab = refractive(&[1000.0, 1000.0], 1.5);
    let boundary = cell(&slab, 0, [0.9, 1.0]);
    close(boundary.layer_coupling()[1], 2.25, 2e-10);
    close(
        boundary.net_flux([100.0; 2], &[100.0, 200.0]).unwrap(),
        225.0,
        2e-8,
    );
    let weak = refractive(&[1e-20, 3e-20], 1.5);
    let narrow = cell(&weak, 0, [0.4, 0.400_000_001]);
    assert!(narrow.layer_coupling()[1] > 0.0);
    close(
        narrow.layer_coupling()[1] / ((0.400_000_001 - 0.4) * weak.layer_coupling()[1]),
        1.0,
        2e-8,
    );
    assert!(narrow.net_flux([0.0; 2], &[0.0, 1e29]).unwrap() > 0.0);
    assert_eq!(narrow.net_flux([300.0; 2], &[300.0; 2]).unwrap(), 0.0);
}

#[test]
fn g0_net_depth_validates_fields_and_bounded_preparation() {
    let slab = slab(&[0.0, 1.0]);
    let no_work = SlabQuadrature {
        maximum_evaluations: 0,
        ..SlabQuadrature::default()
    };
    for (layer, interval) in [(0, [0.0, 1.0]), (1, [0.4, 0.4])] {
        let empty = slab
            .depth_response(layer, interval, no_work, || false)
            .unwrap();
        assert_eq!(empty.net_flux([100.0, -100.0], &[20.0, 30.0]).unwrap(), 0.0);
        assert_eq!(empty.quadrature().angular_evaluations, 0);
    }
    assert!(
        slab.depth_response(1, [0.0, 1.0], no_work, || false)
            .is_err()
    );
    assert!(
        slab.depth_response(2, [0.0, 1.0], SlabQuadrature::default(), || false)
            .is_err()
    );
    assert!(
        slab.depth_response(1, [0.8, 0.2], SlabQuadrature::default(), || false)
            .is_err()
    );
    assert!(matches!(
        slab.depth_response(1, [0.0, 1.0], SlabQuadrature::default(), || true),
        Err(ConductionError::Cancelled { .. })
    ));
    let response = cell(&slab, 1, [0.0, 1.0]);
    assert!(response.net_flux([0.0; 2], &[1.0]).is_err());
    assert!(response.net_flux([f64::NAN, 0.0], &[0.0; 2]).is_err());
    assert!(response.net_flux([0.0; 2], &[0.0, f64::INFINITY]).is_err());
    assert!(response.net_flux([f64::MAX; 2], &[0.0; 2]).is_err());
}
