//! G1 ballistic calorimetry and an independent two-stream transport solution;
//! G3 spatial refinement; G4 refusal and cancellation. Coefficients are
//! manufactured physical laws, not evaluated shielding data.

use super::*;
use crate::radiation::RadiationDeposition;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn normal_pair() -> Vec<PhotonOrdinate> {
    vec![PhotonOrdinate {
        mu: 1.0,
        weight: 1.0,
    }]
}

fn mono(cell_id: u64, thickness_m: f64, heat: f64, scatter: f64) -> PhotonSlabCell {
    PhotonSlabCell {
        cell_id,
        thickness_m,
        groups: vec![PhotonCollision {
            total_per_m: heat + scatter,
            heat_per_m: heat,
            redistribution_per_m: vec![scatter],
        }],
    }
}

fn prepare(
    cells: Vec<PhotonSlabCell>,
    source_ids: Vec<u64>,
    ordinates: Vec<PhotonOrdinate>,
) -> PhotonSlabResponse {
    let edges = (0..=cells[0].groups.len())
        .map(|group| group as f64 * 1e-13)
        .collect();
    prepare_slab_response(
        "manufactured slab transport".into(),
        DecayRadiation::Gamma,
        edges,
        cells,
        source_ids,
        ordinates,
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn apply(response: &PhotonSlabResponse, emitted: &[(u64, &[f64])]) -> RadiationDeposition {
    response.response.apply(emitted, 1e-10, || Ok(())).unwrap()
}

fn heat(result: &RadiationDeposition, cell_id: u64) -> f64 {
    result
        .receiver_energy_j
        .get(&cell_id)
        .copied()
        .unwrap_or(0.0)
}

fn escape_from_uniform_source(optical_thickness: f64) -> f64 {
    if optical_thickness == 0.0 {
        1.0
    } else {
        -(-optical_thickness).exp_m1() / optical_thickness
    }
}

#[test]
fn g1_vacuum_layers_export_every_joule_without_creating_heat() {
    let response = prepare(
        vec![mono(40, 0.2, 0.0, 0.0), mono(7, 3.0, 0.0, 0.0)],
        vec![40, 7],
        gauss_legendre_4(),
    );
    let result = apply(&response, &[(7, &[2.0]), (40, &[5.0])]);
    near(result.escaped_j, 7.0, 2e-14);
    assert_eq!(result.deposited_j, 0.0);
    assert!(result.receiver_energy_j.values().all(|&value| value == 0.0));
    assert_eq!(response.response.source_ids().collect::<Vec<_>>(), [7, 40]);
}

#[test]
fn g1_uniform_absorbing_slab_matches_exact_thin_and_thick_escape() {
    // Average exp(-kappa*x) over a uniform source and its equal +/- directions.
    // This independent Beer-Lambert integral includes within-source absorption.
    for tau in [1e-8, 0.35, 3.0, 100.0] {
        let response = prepare(vec![mono(5, 1.0, tau, 0.0)], vec![5], normal_pair());
        let result = apply(&response, &[(5, &[7.0])]);
        let escape = escape_from_uniform_source(tau);
        near(result.escaped_j, 7.0 * escape, 5e-14);
        near(heat(&result, 5), 7.0 * (1.0 - escape), 5e-14);
        assert!(result.deposited_j > 0.0);
        near(result.deposited_j + result.escaped_j, 7.0, 5e-14);
    }
}

#[test]
fn g1_layered_shielding_keeps_geometric_order_and_source_superposition() {
    let response = prepare(
        vec![
            mono(900, 0.5, 0.5, 0.0),
            mono(8, 1.0, 0.6, 0.0),
            mono(77, 2.0, 0.7, 0.0),
        ],
        vec![77, 900],
        normal_pair(),
    );
    let a = apply(&response, &[(900, &[3.0]), (77, &[0.0])]);
    let b = apply(&response, &[(77, &[5.0]), (900, &[0.0])]);
    let both = apply(&response, &[(900, &[3.0]), (77, &[5.0])]);
    let source_left = escape_from_uniform_source(0.25);
    let source_right = escape_from_uniform_source(1.4);
    let middle_transmission = (-0.6_f64).exp();
    near(heat(&a, 900), 3.0 * (1.0 - source_left), 3e-14);
    near(
        heat(&a, 8),
        1.5 * source_left * (1.0 - middle_transmission),
        3e-14,
    );
    near(
        heat(&a, 77),
        1.5 * source_left * middle_transmission * (1.0 - (-1.4_f64).exp()),
        3e-14,
    );
    near(
        a.escaped_j,
        1.5 * source_left * (1.0 + (-2.0_f64).exp()),
        3e-14,
    );
    near(heat(&b, 77), 5.0 * (1.0 - source_right), 3e-14);
    near(
        heat(&b, 8),
        2.5 * source_right * (1.0 - middle_transmission),
        3e-14,
    );
    near(
        heat(&b, 900),
        2.5 * source_right * middle_transmission * (1.0 - (-0.25_f64).exp()),
        3e-14,
    );
    near(
        b.escaped_j,
        2.5 * source_right * (1.0 + (-0.85_f64).exp()),
        3e-14,
    );
    for id in [900, 8, 77] {
        near(heat(&both, id), heat(&a, id) + heat(&b, id), 3e-14);
    }
    near(both.escaped_j, a.escaped_j + b.escaped_j, 3e-14);
    near(both.deposited_j + both.escaped_j, 8.0, 5e-14);
}

#[test]
fn g0_elastic_multiple_scattering_never_becomes_absorbed_heat() {
    let response = prepare(
        vec![mono(11, 0.5, 0.0, 2.0), mono(22, 0.5, 0.0, 2.0)],
        vec![11, 22],
        gauss_legendre_4(),
    );
    let result = apply(&response, &[(11, &[2.0]), (22, &[5.0])]);
    assert_eq!(result.deposited_j, 0.0);
    near(result.escaped_j, 7.0, 1e-12);
    assert!(result.roundoff_j > 0.0);
    // Preparation corrections remain chargeable even when the final escape
    // coefficient and all subsequent unit-energy products are exactly one.
    assert!(matches!(
        response.response.apply(&[(11, &[1.0]), (22, &[0.0])], 0.0, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
}

#[test]
fn g1_inelastic_downscatter_separates_local_energy_loss_and_lower_group_escape() {
    let response = prepare(
        vec![PhotonSlabCell {
            cell_id: 3,
            thickness_m: 0.7,
            groups: vec![
                PhotonCollision {
                    total_per_m: 0.0,
                    heat_per_m: 0.0,
                    redistribution_per_m: vec![0.0, 0.0],
                },
                PhotonCollision {
                    total_per_m: 2.0,
                    heat_per_m: 1.2,
                    redistribution_per_m: vec![0.8, 0.0],
                },
            ],
        }],
        vec![3],
        normal_pair(),
    );
    // Each high-group collision deposits 60% and releases the remaining 40%
    // into a transparent lower group. No second collision can hide that energy.
    let result = apply(&response, &[(3, &[0.0, 4.0])]);
    let uncollided = escape_from_uniform_source(1.4);
    near(heat(&result, 3), 4.0 * 0.6 * (1.0 - uncollided), 4e-14);
    near(
        result.escaped_j,
        4.0 * (uncollided + 0.4 * (1.0 - uncollided)),
        4e-14,
    );
    near(result.deposited_j + result.escaped_j, 4.0, 4e-14);
    let transparent = apply(&response, &[(3, &[1.0, 0.0])]);
    assert_eq!(transparent.deposited_j, 0.0);
    near(transparent.escaped_j, 1.0, 1e-14);
}

fn uniformly_emitting_slab(cells: usize, absorption: f64, scattering: f64) -> f64 {
    let response = prepare(
        (0..cells)
            .map(|i| mono(i as u64, 2.0 / cells as f64, absorption, scattering))
            .collect(),
        (0..cells as u64).collect(),
        normal_pair(),
    );
    let source_energy = [1.0 / cells as f64];
    let inputs: Vec<_> = (0..cells as u64)
        .map(|id| (id, &source_energy[..]))
        .collect();
    let result = apply(&response, &inputs);
    near(result.deposited_j + result.escaped_j, 1.0, 2e-13);
    result.escaped_j
}

#[test]
fn g3_ballistic_subdivision_is_invariant_and_scattering_converges_to_two_stream_ode() {
    let ballistic = escape_from_uniform_source(1.6);
    for cells in [1, 8, 32] {
        near(uniformly_emitting_slab(cells, 0.8, 0.0), ballistic, 2e-13);
    }

    // For constant a, t and isotropic uniform source q, the independent
    // two-stream equations give phi'' = a*t*phi - t*q (mu=1). Vacuum incoming
    // intensities at x=0,L give this closed-form total escape fraction.
    let absorption = 0.4_f64;
    let total = 4.0_f64;
    let z = 2.0 * (absorption * total).sqrt();
    let ratio = (absorption / total).sqrt();
    let tangent = (z * 0.5).tanh();
    let exact = 2.0 * tangent / (z * (1.0 + ratio * tangent));
    let errors = [1, 8, 64].map(|cells| {
        (uniformly_emitting_slab(cells, absorption, total - absorption) - exact).abs()
    });
    assert!(errors[1] < errors[0] * 0.5, "{errors:?}");
    assert!(errors[2] < errors[1] * 0.5, "{errors:?}");
    assert!(errors[2] < 5e-4, "{errors:?}");
}

#[test]
fn g4_malformed_transport_resource_exhaustion_and_late_cancellation_refuse() {
    let attempt = |cells, edges, sources, ordinates, limits| {
        prepare_slab_response(
            "manufactured refusal".into(),
            DecayRadiation::Gamma,
            edges,
            cells,
            sources,
            ordinates,
            limits,
            || Ok(()),
        )
    };
    let edges = vec![0.0, 1e-13];
    for cells in [
        vec![mono(1, 0.0, 1.0, 0.0)],
        vec![mono(1, f64::NAN, 1.0, 0.0)],
        vec![mono(1, 1.0, 1.0, 0.0), mono(1, 1.0, 1.0, 0.0)],
        vec![mono(1, 1.0, -1.0, 2.0)],
        vec![PhotonSlabCell {
            cell_id: 1,
            thickness_m: 1.0,
            groups: vec![PhotonCollision {
                total_per_m: 2.0,
                heat_per_m: 1.0,
                redistribution_per_m: vec![0.0],
            }],
        }],
    ] {
        assert!(
            attempt(
                cells,
                edges.clone(),
                vec![1],
                normal_pair(),
                PhotonSlabLimits::default()
            )
            .is_err()
        );
    }
    for sources in [vec![2], vec![1, 1]] {
        assert!(
            attempt(
                vec![mono(1, 1.0, 1.0, 0.0)],
                edges.clone(),
                sources,
                normal_pair(),
                PhotonSlabLimits::default()
            )
            .is_err()
        );
    }
    assert!(
        attempt(
            vec![mono(1, 1.0, 1.0, 0.0)],
            vec![0.0, 1e-13, 2e-13],
            vec![1],
            normal_pair(),
            PhotonSlabLimits::default()
        )
        .is_err()
    );
    // An upward energy transfer needs an additional material-energy owner;
    // it cannot disappear from a descending-group cascade solve.
    assert!(
        attempt(
            vec![PhotonSlabCell {
                cell_id: 1,
                thickness_m: 1.0,
                groups: vec![
                    PhotonCollision {
                        total_per_m: 1.0,
                        heat_per_m: 0.0,
                        redistribution_per_m: vec![0.0, 1.0],
                    },
                    PhotonCollision {
                        total_per_m: 0.0,
                        heat_per_m: 0.0,
                        redistribution_per_m: vec![0.0, 0.0],
                    },
                ],
            }],
            vec![0.0, 1e-13, 2e-13],
            vec![1],
            normal_pair(),
            PhotonSlabLimits::default(),
        )
        .is_err()
    );
    for ordinate in [
        PhotonOrdinate {
            mu: 0.0,
            weight: 1.0,
        },
        PhotonOrdinate {
            mu: 1.1,
            weight: 1.0,
        },
        PhotonOrdinate {
            mu: 0.5,
            weight: 0.5,
        },
    ] {
        assert!(
            attempt(
                vec![mono(1, 1.0, 1.0, 0.0)],
                edges.clone(),
                vec![1],
                vec![ordinate],
                PhotonSlabLimits::default()
            )
            .is_err()
        );
    }
    let mut limits = PhotonSlabLimits::default();
    limits.max_work = 0;
    assert!(matches!(
        attempt(
            vec![mono(1, 1.0, 1.0, 0.0)],
            edges.clone(),
            vec![1],
            normal_pair(),
            limits
        ),
        Err(RadiationError::Budget(_))
    ));
    let mut limits = PhotonSlabLimits::default();
    limits.max_values = 0;
    assert!(matches!(
        attempt(
            vec![mono(1, 1.0, 1.0, 0.0)],
            edges,
            vec![1],
            normal_pair(),
            limits
        ),
        Err(RadiationError::Budget(_))
    ));

    let run = |checkpoint: &mut dyn FnMut() -> Result<(), RadiationError>| {
        prepare_slab_response(
            "manufactured cancellation".into(),
            DecayRadiation::Gamma,
            vec![0.0, 1e-13],
            (0..8).map(|i| mono(i, 0.25, 0.4, 3.6)).collect(),
            (0..8).collect(),
            gauss_legendre_4(),
            PhotonSlabLimits::default(),
            checkpoint,
        )
    };
    let mut total_checkpoints = 0;
    run(&mut || {
        total_checkpoints += 1;
        Ok(())
    })
    .unwrap();
    assert!(total_checkpoints > 1);
    let mut observed = 0;
    let cancelled = run(&mut || {
        observed += 1;
        if observed == total_checkpoints {
            Err(RadiationError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(cancelled, Err(RadiationError::Cancelled)));
    assert_eq!(observed, total_checkpoints);
}
