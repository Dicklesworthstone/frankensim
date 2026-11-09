use super::*;

pub(super) fn partition(slab: &AbsorbingSlab, cut: f64) -> SlabCellResponse {
    slab.cell_response(
        (0..slab.layer_count())
            .flat_map(|layer| {
                [
                    SlabCell {
                        layer,
                        interval_m: [0.0, cut],
                    },
                    SlabCell {
                        layer,
                        interval_m: [cut, slab.layers()[layer].thickness_m()],
                    },
                ]
            })
            .collect(),
        SlabQuadrature::default(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_same_material_cell_exchange_matches_independent_exponential_integrals() {
    let original = slab(&[1.0]);
    let response = partition(&original, 0.5);
    // Independent 70-digit E1 series/E3 recurrence:
    // G01=2[E3(0)-2E3(.5)+E3(1)], C0f=1-2E3(.5).
    close(response.cell_coupling()[1], 0.332_966_477_294_806_44, 2e-11);
    close(
        response.face_coupling()[0][0],
        0.556_791_271_449_643_1,
        2e-11,
    );
    close(
        response.face_coupling()[0][1],
        0.223_824_794_154_836_64,
        2e-11,
    );
    assert_eq!(response.cell_coupling()[1], response.cell_coupling()[2]);
    assert_eq!(response.cell_coupling()[0], 0.0);
    assert_eq!(response.cell_coupling()[3], 0.0);
    assert_eq!(response.slab().layer_count(), 1);
    assert_eq!(response.cell_count(), 2);
    assert_eq!(
        response.slab().layers().as_ptr(),
        original.layers().as_ptr()
    );
    let mut heat = [0.0; 2];
    response
        .net_fluxes([100.0; 2], &[200.0, 100.0], &mut heat)
        .unwrap();
    close(heat[1], 100.0 * 0.332_966_477_294_806_44, 2e-9);
    assert!(heat[0] < 0.0 && heat[1] > 0.0);
    for background in [0.0, 600.0, f64::MAX, -f64::MAX] {
        response
            .net_fluxes([background; 2], &[background; 2], &mut heat)
            .unwrap();
        assert_eq!(heat, [0.0; 2]);
    }
}

#[test]
fn g3_cell_face_and_cross_couplings_sum_to_original_layers_and_depth_heating() {
    for original in [
        slab(&[0.2, 0.8, 0.4]),
        asymmetric::between(&[0.3, 0.8, 0.2], &[1.5, 1.0, 2.0], [1.8, 1.2]),
    ] {
        let response = partition(&original, 0.3);
        let n = response.cell_count();
        for layer in 0..original.layer_count() {
            for side in 0..2 {
                let sum: f64 = response
                    .cells()
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| c.layer == layer)
                    .map(|(i, _)| response.face_coupling()[i][side])
                    .sum();
                close(sum, original.face_coupling()[layer][side], 2e-9);
            }
            for other in layer + 1..original.layer_count() {
                let sum: f64 = response
                    .cells()
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| c.layer == layer)
                    .flat_map(|(i, _)| {
                        response
                            .cells()
                            .iter()
                            .enumerate()
                            .filter(move |(_, c)| c.layer == other)
                            .map(move |(j, _)| (i, j))
                    })
                    .map(|(i, j)| response.cell_coupling()[i * n + j])
                    .sum();
                close(
                    sum,
                    original.layer_coupling()[layer * original.layer_count() + other],
                    2e-9,
                );
            }
        }
        let exitances = [300.0, 150.0, 450.0];
        let cell_exitances: Vec<_> = response
            .cells()
            .iter()
            .map(|c| exitances[c.layer])
            .collect();
        let mut heat = vec![0.0; n];
        response
            .net_fluxes([73.0, 81.0], &cell_exitances, &mut heat)
            .unwrap();
        let mut whole = [0.0; 3];
        original
            .layer_net_fluxes([73.0, 81.0], &exitances, &mut whole)
            .unwrap();
        for layer in 0..3 {
            let total: f64 = response
                .cells()
                .iter()
                .enumerate()
                .filter(|(_, c)| c.layer == layer)
                .map(|(i, _)| heat[i])
                .sum();
            close(total, whole[layer], 2e-7);
        }
        for (i, cell) in response.cells().iter().enumerate() {
            let old = original
                .depth_response(
                    cell.layer,
                    cell.interval_m,
                    SlabQuadrature::default(),
                    || false,
                )
                .unwrap();
            for side in 0..2 {
                close(
                    response.face_coupling()[i][side],
                    old.face_coupling()[side],
                    2e-9,
                );
            }
            close(
                heat[i],
                old.net_flux([73.0, 81.0], &exitances).unwrap(),
                2e-7,
            );
        }
    }
}

#[test]
fn g1_cell_returns_retain_asymmetric_reversal_and_all_trapped_modes() {
    let original = asymmetric::between(&[1000.0], &[2.0], [1.5, 1.1]);
    let response = partition(&original, 0.4);
    close(response.cell_coupling()[1], 4.0, 3e-13);
    let reverse = partition(&asymmetric::between(&[1000.0], &[2.0], [1.1, 1.5]), 0.6);
    for i in 0..2 {
        for side in 0..2 {
            close(
                response.face_coupling()[i][side],
                reverse.face_coupling()[1 - i][1 - side],
                1e-12,
            );
        }
    }
    let mut heat = [0.0; 2];
    // Use the actual authored index product for an exactly representable bath.
    let equilibrium = original
        .exterior_refractive_indices()
        .map(|x| x * x * 500.0);
    response
        .net_fluxes(equilibrium, &[500.0; 2], &mut heat)
        .unwrap();
    assert_eq!(heat, [0.0; 2]);
    response
        .net_fluxes([300.0, 200.0], &[700.0, 450.0], &mut heat)
        .unwrap();
    let mut reversed = [0.0; 2];
    reverse
        .net_fluxes([200.0, 300.0], &[450.0, 700.0], &mut reversed)
        .unwrap();
    close(heat[0], reversed[1], 2e-9);
    close(heat[1], reversed[0], 2e-9);
}

#[test]
fn g1_weak_narrow_cells_resolve_positive_exchange_through_original_reflections() {
    let controls = SlabQuadrature {
        absolute_coefficient_error: 0.0,
        relative_coefficient_error: 1e-8,
        maximum_evaluations: 65_536,
        maximum_depth: 48,
    };
    let make = |widths: &[f64]| {
        AbsorbingSlab::new_refractive(
            widths
                .iter()
                .map(|&width| {
                    AbsorbingLayer::new(
                        width,
                        1e-20,
                        [100.0, 3000.0],
                        "same original weak absorption",
                    )
                    .unwrap()
                })
                .collect(),
            1.5,
            "same original real index",
            controls,
            || false,
        )
        .unwrap()
    };
    let original = make(&[1.0]);
    let narrow = 1e-9;
    let response = original
        .cell_response(
            vec![
                SlabCell {
                    layer: 0,
                    interval_m: [0.0, narrow],
                },
                SlabCell {
                    layer: 0,
                    interval_m: [narrow, 1.0],
                },
            ],
            controls,
            || false,
        )
        .unwrap();
    // Independent existing layer-exchange equations at a zero-reflection
    // equal-index subdivision are an oracle, not the production cell method.
    let reference = make(&[narrow, 1.0 - narrow]);
    let exchange = response.cell_coupling()[1];
    assert!(exchange > 0.0 && exchange < 1e-27);
    assert!((exchange / reference.layer_coupling()[1] - 1.0).abs() < 2e-8);
    for i in 0..2 {
        for side in 0..2 {
            let expected = reference.face_coupling()[i][side];
            assert!(response.face_coupling()[i][side] > 0.0);
            assert!((response.face_coupling()[i][side] / expected - 1.0).abs() < 2e-8);
        }
    }
}

#[test]
fn g0_cell_partition_work_limits_and_heat_publication_are_explicit() {
    let original = slab(&[0.3]);
    let interval = |start, end| SlabCell {
        layer: 0,
        interval_m: [start, end],
    };
    for cells in [
        vec![],
        vec![interval(0.0, 0.5)],
        vec![interval(0.1, 1.0)],
        vec![interval(0.0, 0.6), interval(0.5, 1.0)],
        vec![interval(0.0, 0.4), interval(0.5, 1.0)],
        vec![interval(0.0, 0.0), interval(0.0, 1.0)],
        vec![interval(0.0, f64::NAN)],
        vec![SlabCell {
            layer: 1,
            interval_m: [0.0, 1.0],
        }],
        vec![interval(0.0, 1.0); MAX_SLAB_CELLS + 1],
    ] {
        assert!(
            original
                .cell_response(cells, SlabQuadrature::default(), || false)
                .is_err()
        );
    }
    let zero_work = SlabQuadrature {
        maximum_evaluations: 0,
        ..SlabQuadrature::default()
    };
    assert!(
        original
            .cell_response(vec![interval(0.0, 1.0)], zero_work, || false)
            .is_err()
    );
    assert!(
        original
            .cell_response(vec![interval(0.0, 1.0)], SlabQuadrature::default(), || true)
            .is_err()
    );
    let clear = slab(&[0.0])
        .cell_response(
            vec![interval(0.0, 0.2), interval(0.2, 1.0)],
            zero_work,
            || false,
        )
        .unwrap();
    assert_eq!(clear.quadrature().angular_evaluations, 0);
    assert_eq!(clear.cell_coupling(), &[0.0; 4]);
    let response = partition(&original, 0.5);
    let mut out = [37.0; 2];
    assert!(
        response
            .net_fluxes([10.0; 2], &[100.0, f64::NAN], &mut out)
            .is_err()
    );
    assert_eq!(out, [37.0; 2]);
}
