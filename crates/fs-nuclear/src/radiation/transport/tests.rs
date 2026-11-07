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

fn beam(source_id: u64, boundary: PhotonSlabBoundary, mu: f64) -> PhotonSlabSource {
    PhotonSlabSource {
        source_id,
        geometry: PhotonSourceGeometry::Collimated { boundary, mu },
    }
}

fn prepare_sources(
    cells: Vec<PhotonSlabCell>,
    sources: Vec<PhotonSlabSource>,
) -> PhotonSlabResponse {
    prepare_angular_sources(cells, sources, normal_pair())
}

fn prepare_angular_sources(
    cells: Vec<PhotonSlabCell>,
    sources: Vec<PhotonSlabSource>,
    ordinates: Vec<PhotonOrdinate>,
) -> PhotonSlabResponse {
    let edges = (0..=cells[0].groups.len())
        .map(|group| group as f64 * 1e-13)
        .collect();
    prepare_slab_response_with_sources(
        "manufactured internal and external sources".into(),
        DecayRadiation::Gamma,
        edges,
        cells,
        sources,
        ordinates,
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn diffuse(source_id: u64, boundary: PhotonSlabBoundary) -> PhotonSlabSource {
    PhotonSlabSource {
        source_id,
        geometry: PhotonSourceGeometry::Diffuse { boundary },
    }
}

#[test]
fn g1_g3_diffuse_absorption_refines_to_the_exponential_integral_solution() {
    // Constant incident radiance gives T(tau) = 2*integral_0^1 mu*exp(-tau/mu)dmu
    // = 2*E_3(tau), DLMF https://dlmf.nist.gov/8.19.E3. At tau=1, recurrence
    // gives 2*E_3(1)=E_1(1); this independently evaluated value is not a beam sum.
    let exact = 0.219_383_934_395_520_27;
    let errors = [8, 32, 128].map(|count| {
        let ordinates = (0..count)
            .map(|i| PhotonOrdinate {
                mu: (i as f64 + 0.5) / count as f64,
                weight: 1.0 / count as f64,
            })
            .collect();
        let response = prepare_angular_sources(
            vec![mono(5, 1.0, 1.0, 0.0)],
            vec![diffuse(99, PhotonSlabBoundary::Left)],
            ordinates,
        );
        let result = apply(&response, &[(99, &[1.0])]);
        near(result.deposited_j + result.escaped_j, 1.0, 2e-14);
        (result.escaped_j - exact).abs()
    });
    assert!(errors[1] < errors[0] / 10.0, "{errors:?}");
    assert!(errors[2] < errors[1] / 10.0, "{errors:?}");
    assert!(errors[2] < 4e-6, "{errors:?}");
}

#[test]
fn g1_g3_diffuse_illumination_preserves_vacuum_energy_and_layer_orientation() {
    // The original mirrored full-range rule does not have first moment 1/2.
    // A constant-radiance source still supplies exactly the declared joules.
    for opacity in [0.0, 1e-12] {
        let result = apply(
            &prepare_angular_sources(
                vec![mono(5, 1.0, opacity, 0.0)],
                vec![diffuse(99, PhotonSlabBoundary::Left)],
                gauss_legendre_4(),
            ),
            &[(99, &[7.0])],
        );
        near(result.deposited_j + result.escaped_j, 7.0, 3e-14);
        if opacity == 0.0 {
            assert_eq!(result.deposited_j, 0.0);
        } else {
            assert!(result.deposited_j > 0.0);
        }
    }
    let cells = vec![mono(50, 0.5, 0.6, 0.0), mono(3, 2.0, 0.7, 0.0)];
    let left = prepare_angular_sources(
        cells.clone(),
        vec![diffuse(99, PhotonSlabBoundary::Left)],
        gauss_legendre_4(),
    );
    let reversed = prepare_angular_sources(
        cells.into_iter().rev().collect(),
        vec![diffuse(99, PhotonSlabBoundary::Right)],
        gauss_legendre_4(),
    );
    let a = apply(&left, &[(99, &[3.0])]);
    let b = apply(&reversed, &[(99, &[3.0])]);
    assert_eq!(a.receiver_energy_j, b.receiver_energy_j);
    assert_eq!(a.escaped_j, b.escaped_j);
}

#[test]
fn g1_diffuse_multigroup_scattering_matches_cosine_weighted_illumination() {
    let cells = vec![PhotonSlabCell {
        cell_id: 5,
        thickness_m: 1.0,
        groups: vec![
            PhotonCollision {
                total_per_m: 1.0,
                heat_per_m: 0.3,
                redistribution_per_m: vec![0.7, 0.0],
            },
            PhotonCollision {
                total_per_m: 2.0,
                heat_per_m: 0.4,
                redistribution_per_m: vec![0.8, 0.8],
            },
        ],
    }];
    let ordinates = gauss_legendre_4();
    let response = prepare_angular_sources(
        cells.clone(),
        vec![diffuse(99, PhotonSlabBoundary::Left)],
        ordinates.clone(),
    );
    let result = apply(&response, &[(99, &[2.0, 7.0])]);
    let normalization: f64 = ordinates.iter().map(|o| o.mu * o.weight).sum();
    let mut expected_heat = 0.0;
    for ordinate in &ordinates {
        let response = prepare_angular_sources(
            cells.clone(),
            vec![beam(99, PhotonSlabBoundary::Left, ordinate.mu)],
            ordinates.clone(),
        );
        let energy = apply(&response, &[(99, &[2.0, 7.0])]);
        expected_heat += ordinate.mu * ordinate.weight / normalization * energy.deposited_j;
    }
    near(result.deposited_j, expected_heat, 5e-14);
    near(result.deposited_j + result.escaped_j, 9.0, 5e-14);
    let elastic = prepare_angular_sources(
        vec![mono(5, 1.0, 0.0, 2.0)],
        vec![diffuse(99, PhotonSlabBoundary::Left)],
        ordinates,
    );
    let elastic = apply(&elastic, &[(99, &[7.0])]);
    assert_eq!(elastic.deposited_j, 0.0);
    near(elastic.escaped_j, 7.0, 5e-14);
}

#[test]
fn g1_incident_beam_matches_thin_and_thick_oblique_beer_lambert_absorption() {
    for mu in [1.0, 0.37] {
        for tau in [1e-12, 0.7, 3.0, 100.0] {
            let response = prepare_sources(
                vec![mono(5, 1.0, tau, 0.0)],
                vec![beam(99, PhotonSlabBoundary::Left, mu)],
            );
            let result = apply(&response, &[(99, &[7.0])]);
            let absorbed = 7.0 * (-(-tau / mu).exp_m1());
            near(heat(&result, 5), absorbed, 5e-14 * absorbed);
            near(result.escaped_j, 7.0 * (-tau / mu).exp(), 5e-14);
            near(result.deposited_j + result.escaped_j, 7.0, 5e-14);
            assert!(heat(&result, 5) > 0.0);
            assert!(!result.receiver_energy_j.contains_key(&99));
        }
    }
}

#[test]
fn g1_beam_layer_heat_and_reversed_illumination_follow_ordered_path_lengths() {
    let cells = vec![
        mono(50, 0.5, 0.6, 0.0),
        mono(7, 0.8, 0.0, 0.0),
        mono(3, 2.0, 0.7, 0.0),
    ];
    let mu = 0.4;
    let left = prepare_sources(cells.clone(), vec![beam(99, PhotonSlabBoundary::Left, mu)]);
    let result = apply(&left, &[(99, &[3.0])]);
    let front_pass = (-0.3_f64 / mu).exp();
    let rear_pass = (-1.4_f64 / mu).exp();
    near(heat(&result, 50), 3.0 * (1.0 - front_pass), 3e-14);
    assert_eq!(heat(&result, 7), 0.0);
    near(
        heat(&result, 3),
        3.0 * front_pass * (1.0 - rear_pass),
        3e-14,
    );
    near(result.escaped_j, 3.0 * front_pass * rear_pass, 3e-14);

    let reversed = prepare_sources(
        cells.into_iter().rev().collect(),
        vec![beam(99, PhotonSlabBoundary::Right, mu)],
    );
    let reversed = apply(&reversed, &[(99, &[3.0])]);
    assert_eq!(result.receiver_energy_j, reversed.receiver_energy_j);
    assert_eq!(result.escaped_j, reversed.escaped_j);
}

#[test]
fn g1_beam_repeated_scattering_matches_single_cell_collision_series() {
    // A first collision creates local heat h and an isotropic outgoing fraction
    // a. Each ensuing flight collides inside the cell with probability c, so
    // total heat is q*h * sum_k (a*c)^k. This tests the infinite scattering
    // series, rather than treating total beam attenuation as absorbed energy.
    let thickness = 0.8;
    let absorption = 0.4;
    let scattering = 1.6;
    let mu = 0.6;
    let tau = thickness * (absorption + scattering);
    let response = prepare_sources(
        vec![mono(5, thickness, absorption, scattering)],
        vec![beam(99, PhotonSlabBoundary::Right, mu)],
    );
    let result = apply(&response, &[(99, &[7.0])]);
    let first_collision = -(-tau / mu).exp_m1();
    let uniform_collision = 1.0 - escape_from_uniform_source(tau);
    let expected = 7.0 * first_collision * 0.2 / (1.0 - 0.8 * uniform_collision);
    near(heat(&result, 5), expected, 5e-14);
    near(result.escaped_j, 7.0 - expected, 5e-14);
    assert!(heat(&result, 5) < 7.0 * first_collision);

    let elastic = prepare_sources(
        vec![mono(5, thickness, 0.0, 2.0), mono(8, 0.3, 0.0, 0.8)],
        vec![beam(99, PhotonSlabBoundary::Left, mu)],
    );
    let result = apply(&elastic, &[(99, &[7.0])]);
    assert_eq!(result.deposited_j, 0.0);
    near(result.escaped_j, 7.0, 8e-14);
}

#[test]
fn g1_beam_downscatter_reuses_lower_group_absorption_without_losing_energy() {
    let response = prepare_sources(
        vec![PhotonSlabCell {
            cell_id: 5,
            thickness_m: 1.0,
            groups: vec![
                PhotonCollision {
                    total_per_m: 1.0,
                    heat_per_m: 1.0,
                    redistribution_per_m: vec![0.0, 0.0],
                },
                PhotonCollision {
                    total_per_m: 2.0,
                    heat_per_m: 1.2,
                    redistribution_per_m: vec![0.8, 0.0],
                },
            ],
        }],
        vec![beam(99, PhotonSlabBoundary::Left, 0.5)],
    );
    let result = apply(&response, &[(99, &[2.0, 7.0])]);
    let high_collision = -(-4.0_f64).exp_m1();
    let low_escape = escape_from_uniform_source(1.0);
    let expected_heat =
        2.0 * (-(-2.0_f64).exp_m1()) + 7.0 * high_collision * (0.6 + 0.4 * (1.0 - low_escape));
    near(heat(&result, 5), expected_heat, 5e-14);
    near(result.escaped_j, 9.0 - expected_heat, 5e-14);
    near(result.deposited_j + result.escaped_j, 9.0, 5e-14);
}

#[test]
fn g3_internal_and_external_sources_share_one_response_and_superpose() {
    let cells = vec![mono(5, 0.6, 0.4, 0.8), mono(9, 0.3, 0.7, 0.2)];
    let response = prepare_sources(
        cells.clone(),
        vec![
            beam(99, PhotonSlabBoundary::Left, 0.8),
            PhotonSlabSource {
                source_id: 42,
                geometry: PhotonSourceGeometry::UniformLayer { cell_id: 5 },
            },
        ],
    );
    let old = prepare(cells, vec![5], normal_pair());
    let internal = apply(&response, &[(42, &[3.0]), (99, &[0.0])]);
    let old = apply(&old, &[(5, &[3.0])]);
    assert_eq!(internal.receiver_energy_j, old.receiver_energy_j);
    assert_eq!(internal.escaped_j, old.escaped_j);
    let external = apply(&response, &[(42, &[0.0]), (99, &[2.0])]);
    let both = apply(&response, &[(42, &[3.0]), (99, &[2.0])]);
    for cell in [5, 9] {
        near(
            heat(&both, cell),
            heat(&internal, cell) + heat(&external, cell),
            3e-14,
        );
    }
    near(
        both.escaped_j,
        internal.escaped_j + external.escaped_j,
        3e-14,
    );
    near(both.deposited_j + both.escaped_j, 5.0, 5e-14);
}

#[test]
fn g4_beam_sources_refuse_invalid_incidence_identity_budgets_and_cancellation() {
    let run = |sources, limits, checkpoint: &mut dyn FnMut() -> Result<(), RadiationError>| {
        prepare_slab_response_with_sources(
            "manufactured beam admission".into(),
            DecayRadiation::Gamma,
            vec![0.0, 1e-13],
            vec![mono(5, 1.0, 0.4, 0.6)],
            sources,
            normal_pair(),
            limits,
            checkpoint,
        )
    };
    let limits = PhotonSlabLimits::default();
    for mu in [0.0, -0.5, 1.1, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            run(
                vec![beam(99, PhotonSlabBoundary::Left, mu)],
                limits,
                &mut || Ok(())
            ),
            Err(RadiationError::InvalidInput(_))
        ));
    }
    let source = beam(99, PhotonSlabBoundary::Left, 0.5);
    assert!(matches!(
        run(vec![source, source], limits, &mut || Ok(())),
        Err(RadiationError::InvalidInput(_))
    ));
    assert!(matches!(
        run(
            vec![PhotonSlabSource {
                source_id: 99,
                geometry: PhotonSourceGeometry::UniformLayer { cell_id: 77 },
            }],
            limits,
            &mut || Ok(())
        ),
        Err(RadiationError::InvalidInput(_))
    ));
    let mut count = 0;
    let completed = run(vec![source], limits, &mut || {
        count += 1;
        Ok(())
    })
    .unwrap();
    let mut short = limits;
    short.max_work = completed.work - 1;
    assert!(matches!(
        run(vec![source], short, &mut || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    short = limits;
    short.max_values = 0;
    assert!(matches!(
        run(vec![source], short, &mut || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    short = limits;
    short.response.max_sources = 0;
    assert!(matches!(
        run(vec![source], short, &mut || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    // Cancel during preparation, including after all physical work is finished;
    // no partial response may be published at either boundary.
    for stop in [1, count / 2, count] {
        let mut observed = 0;
        assert!(matches!(
            run(vec![source], limits, &mut || {
                observed += 1;
                if observed == stop {
                    Err(RadiationError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(RadiationError::Cancelled)
        ));
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
        response
            .response
            .apply(&[(11, &[1.0]), (22, &[0.0])], 0.0, || Ok(())),
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
fn g3_uncoupled_groups_keep_the_same_physical_residual_allowance() {
    let run = |count: usize| {
        let groups = (0..count)
            .map(|incident| {
                let mut outgoing = vec![0.0; count];
                outgoing[incident] = 0.9;
                PhotonCollision {
                    total_per_m: 1.0,
                    heat_per_m: 0.1,
                    redistribution_per_m: outgoing,
                }
            })
            .collect();
        let prepared = prepare(
            vec![PhotonSlabCell {
                cell_id: 7,
                thickness_m: 1.0,
                groups,
            }],
            vec![7],
            normal_pair(),
        );
        let mut emission = vec![0.0; count];
        emission[count - 1] = 1.0;
        prepared
            .response
            .apply(&[(7, &emission)], 1e-12, || Ok(()))
            .unwrap()
    };
    let single = run(1);
    let grouped = run(64);
    // Only one group is active in each row. The correctly rounded scalar solve
    // has a nonzero FMA residual that cannot be improved by correction; adding
    // independent energy groups must not divide its physical energy allowance.
    assert_eq!(grouped.receiver_energy_j, single.receiver_energy_j);
    assert_eq!(grouped.escaped_j, single.escaped_j);
    near(grouped.deposited_j + grouped.escaped_j, 1.0, 1e-14);
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
    // Mirroring a positive ordinate cannot silently erase either hemisphere.
    assert!(matches!(
        attempt(
            vec![mono(1, 1.0, 1.0, 0.0)],
            edges.clone(),
            vec![1],
            vec![
                PhotonOrdinate {
                    mu: 1.0,
                    weight: 1.0
                },
                PhotonOrdinate {
                    mu: 1.0,
                    weight: f64::from_bits(1)
                },
            ],
            PhotonSlabLimits::default()
        ),
        Err(RadiationError::Arithmetic)
    ));
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
