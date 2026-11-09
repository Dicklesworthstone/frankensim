use super::*;

fn response(slab: &AbsorbingSlab, layer: usize, interval: [f64; 2]) -> SlabDiffuseDepthResponse {
    slab.diffuse_depth_response(layer, interval, SlabQuadrature::default(), || false)
        .unwrap()
}

#[test]
fn g1_diffuse_depth_absorption_matches_exponential_integral_oracles() {
    let slab = slab(&[1.0]);
    // Independent 70-digit E1 power series and E3 recurrence, not angular
    // quadrature: A_front = 2[E3(k lo)-E3(k hi)], and reverse depths for back.
    for (interval, expected) in [
        (
            [0.0, 0.25],
            [0.350_631_748_043_712_7, 0.090_149_411_052_300_28],
        ),
        (
            [0.25, 0.5],
            [0.206_159_523_405_930_37, 0.133_675_383_102_536_35],
        ),
        ([0.0, 1.0], [0.780_616_065_604_479_7; 2]),
    ] {
        let prepared = response(&slab, 0, interval);
        for side in 0..2 {
            close(prepared.absorptance()[side], expected[side], 2e-11);
        }
        close(
            prepared.absorbed_flux([120.0, 35.0]).unwrap(),
            120.0 * expected[0] + 35.0 * expected[1],
            3e-9,
        );
        assert_eq!(prepared.layers(), slab.layers());
        assert_eq!(prepared.layers().as_ptr(), slab.layers().as_ptr());
        assert_eq!(prepared.layer_index(), 0);
        assert_eq!(prepared.interval_m(), interval);
        assert!(prepared.quadrature().angular_evaluations >= 24);
        assert_eq!(prepared.quadrature().angular_evaluations % 24, 0);
        assert!(prepared.quadrature().accepted_panels > 0);
        for side in 0..2 {
            let controls = prepared.quadrature().controls;
            assert!(
                prepared.quadrature().estimated_absorptance_error[side]
                    <= controls.absolute_coefficient_error
                        + controls.relative_coefficient_error * prepared.absorptance()[side]
            );
        }
    }
    let thick = AbsorbingSlab::new(
        vec![AbsorbingLayer::new(10.0, 1.0, [100.0, 3000.0], "thick material").unwrap()],
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    let surface = response(&thick, 0, [0.0, 0.1]);
    close(surface.absorptance()[0], 0.167_417_084_183_442_48, 2e-11);
    close(surface.absorptance()[1], 8.091_736_706_588_389e-7, 1e-14);
    // One percent of the physical depth absorbs over sixteen percent of front
    // irradiation. Uniformly distributing whole-layer heat misses this profile.
    assert!(surface.absorptance()[0] > 0.16);
}

#[test]
fn g3_diffuse_cells_sum_to_original_matched_and_fresnel_layer_heat() {
    for indices in [[1.0; 3], [1.4, 2.1, 1.2]] {
        let slab = AbsorbingSlab::new(
            [0.2, 0.9, 0.4]
                .into_iter()
                .zip(indices)
                .enumerate()
                .map(|(i, (depth, index))| {
                    AbsorbingLayer::new(
                        1.0,
                        depth,
                        [100.0, 3000.0],
                        format!("cell material {i}"),
                    )
                    .unwrap()
                    .with_refractive_index(index, format!("cell index {i}"))
                    .unwrap()
                })
                .collect(),
            SlabQuadrature::default(),
            || false,
        )
        .unwrap();
        for layer in 0..3 {
            let mut summed = [0.0; 2];
            for interval in [[0.0, 0.03], [0.03, 0.3], [0.3, 0.85], [0.85, 1.0]] {
                let cell = response(&slab, layer, interval);
                for side in 0..2 {
                    summed[side] += cell.absorptance()[side];
                }
            }
            let whole = response(&slab, layer, [0.0, 1.0]);
            for side in 0..2 {
                close(summed[side], whole.absorptance()[side], 2e-10);
                close(summed[side], slab.face_coupling()[layer][side], 2e-10);
            }
        }
    }
}

#[test]
fn g3_diffuse_cell_reversal_retains_original_refractive_interfaces() {
    let layers: Vec<_> = [(0.7, 0.8, 1.5), (1.3, 0.1, 1.1), (0.4, 2.0, 1.8)]
        .into_iter()
        .enumerate()
        .map(|(i, (d, k, n))| {
            AbsorbingLayer::new(d, k, [100.0, 3000.0], format!("original material {i}"))
                .unwrap()
                .with_refractive_index(n, format!("original index {i}"))
                .unwrap()
        })
        .collect();
    let reversed = AbsorbingSlab::new(
        layers.iter().rev().cloned().collect(),
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    let original = AbsorbingSlab::new(layers, SlabQuadrature::default(), || false).unwrap();
    let forward = response(&original, 0, [0.12, 0.4]);
    let backward = response(&reversed, 2, [0.7 - 0.4, 0.7 - 0.12]);
    for side in 0..2 {
        close(
            forward.absorptance()[side],
            backward.absorptance()[1 - side],
            2e-10,
        );
    }
    close(
        forward.absorbed_flux([27.0, 13.0]).unwrap(),
        backward.absorbed_flux([13.0, 27.0]).unwrap(),
        6e-9,
    );
}

#[test]
fn g1_diffuse_narrow_weak_and_shallow_opaque_cells_keep_positive_heat() {
    let weak = slab(&[1e-20]);
    let interval = [0.4, 0.400_000_001];
    let cell = response(&weak, 0, interval);
    let thin_limit = 2e-20 * (interval[1] - interval[0]);
    for fraction in cell.absorptance() {
        assert!(fraction > 0.0);
        close(fraction / thin_limit, 1.0, 2e-13);
    }
    assert!(cell.absorbed_flux([1e29, 0.0]).unwrap() > 1.0);
    // Neither the original infinite optical depth nor an interval fraction of
    // it is used for this finite, shallow penetration calculation.
    let opaque = AbsorbingSlab::new(
        vec![AbsorbingLayer::new(1e308, 1e308, [100.0, 3000.0], "opaque material").unwrap()],
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    let shallow = response(&opaque, 0, [0.0, 1e-309]);
    close(shallow.absorptance()[0], 0.167_417_084_183_442_48, 2e-11);
    assert_eq!(shallow.absorptance()[1], 0.0);
}

#[test]
fn g0_diffuse_zero_width_and_clear_cells_need_no_angular_work() {
    let no_work = SlabQuadrature {
        maximum_evaluations: 0,
        ..SlabQuadrature::default()
    };
    for (slab, interval) in [
        (slab(&[1.0]), [0.3, 0.3]),
        (refractive(&[0.0], 1.5), [0.0, 1.0]),
    ] {
        let cell = slab
            .diffuse_depth_response(0, interval, no_work, || false)
            .unwrap();
        assert_eq!(cell.absorptance(), [0.0; 2]);
        assert_eq!(cell.absorbed_flux([f64::MAX; 2]).unwrap(), 0.0);
        assert_eq!(cell.quadrature().angular_evaluations, 0);
        assert_eq!(cell.quadrature().accepted_panels, 0);
        assert_eq!(cell.quadrature().estimated_absorptance_error, [0.0; 2]);
    }
}

#[test]
fn g0_diffuse_depth_admission_cancellation_and_work_limits_refuse() {
    let slab = slab(&[1.0]);
    for (layer, interval) in [
        (1, [0.0, 1.0]),
        (0, [-0.1, 0.2]),
        (0, [0.2, 0.1]),
        (0, [0.0, 1.1]),
        (0, [0.0, f64::NAN]),
        (0, [f64::INFINITY; 2]),
    ] {
        assert!(
            slab.diffuse_depth_response(layer, interval, SlabQuadrature::default(), || false)
                .is_err()
        );
    }
    for controls in [
        SlabQuadrature {
            maximum_evaluations: 23,
            ..SlabQuadrature::default()
        },
        SlabQuadrature {
            maximum_evaluations: 65_537,
            ..SlabQuadrature::default()
        },
        SlabQuadrature {
            maximum_depth: 0,
            ..SlabQuadrature::default()
        },
        SlabQuadrature {
            absolute_coefficient_error: 0.0,
            relative_coefficient_error: 0.0,
            ..SlabQuadrature::default()
        },
        SlabQuadrature {
            absolute_coefficient_error: 1e-15,
            relative_coefficient_error: 1e-14,
            maximum_depth: 1,
            ..SlabQuadrature::default()
        },
    ] {
        assert!(
            slab.diffuse_depth_response(0, [0.0, 0.01], controls, || false)
                .is_err()
        );
    }
    let mut checkpoints = 0;
    assert!(matches!(
        slab.diffuse_depth_response(0, [0.0, 1.0], SlabQuadrature::default(), || {
            checkpoints += 1;
            checkpoints == 7
        }),
        Err(ConductionError::Cancelled { .. })
    ));
    let cell = response(&slab, 0, [0.0, 1.0]);
    for irradiation in [
        [-1.0, 0.0],
        [f64::NAN, 0.0],
        [0.0, f64::INFINITY],
        [f64::MAX; 2],
    ] {
        assert!(cell.absorbed_flux(irradiation).is_err());
    }
}
