//! G1/G3: independent slab calorimetry and two-stream convergence. G4: one
//! shared refinement budget and all-or-nothing cancellation of prepared maps.

use super::*;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn mono(id: u64, width: f64, absorption: f64, scattering: f64) -> PhotonSlabCell {
    PhotonSlabCell {
        cell_id: id,
        thickness_m: width,
        groups: vec![PhotonCollision {
            total_per_m: absorption + scattering,
            heat_per_m: absorption,
            redistribution_per_m: vec![scattering],
        }],
    }
}

fn accuracy(tolerance: f64, max_levels: usize) -> PhotonSlabRefinement {
    PhotonSlabRefinement {
        absolute_tolerance: tolerance,
        relative_tolerance: 0.0,
        max_levels,
    }
}

fn run(
    cells: Vec<PhotonSlabCell>,
    sources: Vec<u64>,
    accuracy: PhotonSlabRefinement,
    limits: PhotonSlabLimits,
    checkpoint: impl FnMut() -> Result<(), RadiationError>,
) -> Result<RefinedPhotonSlabResponse, RadiationError> {
    let edges = (0..=cells[0].groups.len())
        .map(|i| i as f64 * 1e-13)
        .collect();
    refine_slab_response(
        "manufactured refined slab".into(),
        DecayRadiation::Gamma,
        edges,
        cells,
        sources,
        vec![PhotonOrdinate {
            mu: 1.0,
            weight: 1.0,
        }],
        limits,
        accuracy,
        checkpoint,
    )
}

#[test]
fn g1_ballistic_response_converges_after_one_split_and_retains_uniform_source_heat() {
    let refined = run(
        vec![mono(91, 0.8, 0.7, 0.0)],
        vec![91],
        accuracy(1e-12, 1),
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(refined.levels, 1);
    assert_eq!(refined.subdivisions_per_layer, 2);
    assert!(refined.max_absolute_change <= 1e-12);
    let output = refined
        .transport
        .response
        .apply(&[(91, &[8.0])], 1e-10, || Ok(()))
        .unwrap();
    let escape = -(-0.56_f64).exp_m1() / 0.56;
    near(output.escaped_j, 8.0 * escape, 5e-14);
    near(output.receiver_energy_j[&91], 8.0 * (1.0 - escape), 5e-14);
    near(output.escaped_j + output.deposited_j, 8.0, 5e-14);
}

#[test]
fn g3_tighter_spatial_tolerance_improves_the_independent_two_stream_solution() {
    let solve = |tolerance| {
        run(
            vec![mono(6, 2.0, 0.4, 3.6)],
            vec![6],
            accuracy(tolerance, 7),
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap()
    };
    let loose = solve(1e-2);
    let tight = solve(8e-4);
    // One layer and two equal half-layers have the same mean collision source
    // by symmetry. Their agreement must never admit scattering at level one.
    assert!(
        loose.levels >= 2,
        "coarse scattering symmetry plateau was accepted"
    );
    assert!(
        tight.levels >= 2,
        "coarse scattering symmetry plateau was accepted"
    );
    // The uniform-source two-stream ODE is phi''=a*t*phi-t*q for mu=1.
    // Vacuum incoming boundaries give this closed form independently of the
    // collision-probability discretization used by the production solver.
    let z = 2.0 * (0.4_f64 * 4.0).sqrt();
    let r = (0.4_f64 / 4.0).sqrt();
    let tangent = (z * 0.5).tanh();
    let exact = 2.0 * tangent / (z * (1.0 + r * tangent));
    let escape = |result: &RefinedPhotonSlabResponse| {
        result
            .transport
            .response
            .apply(&[(6, &[1.0])], 1e-10, || Ok(()))
            .unwrap()
            .escaped_j
    };
    assert!(tight.levels > loose.levels);
    assert!(tight.transport.work > loose.transport.work);
    assert!(tight.max_absolute_change <= 8e-4);
    assert!((escape(&tight) - exact).abs() < (escape(&loose) - exact).abs());
    near(escape(&tight), exact, 5e-4);
    // Preparation errors survive reduction and remain charged in actual joules.
    let preparation = tight.transport.response.preparation_roundoff(6, 0).unwrap();
    assert!(preparation > 0.0);
    assert!(matches!(
        tight
            .transport
            .response
            .apply(&[(6, &[1.0])], preparation * 0.5, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
}

#[test]
fn g4_original_ids_zero_receivers_and_shared_refinement_limits_are_preserved() {
    let layer = |cell_id, thickness_m, heats: [f64; 2]| PhotonSlabCell {
        cell_id,
        thickness_m,
        groups: heats
            .into_iter()
            .map(|heat| PhotonCollision {
                total_per_m: heat,
                heat_per_m: heat,
                redistribution_per_m: vec![0.0; 2],
            })
            .collect(),
    };
    let cells = vec![
        layer(900, 0.7, [0.3, 0.8]),
        layer(u64::MAX, 1.4, [0.5, 1.2]),
        layer(5, 0.9, [0.0, 0.0]),
    ];
    let sources = vec![u64::MAX, 900];
    let config = accuracy(1e-12, 1);
    let mut checkpoints = 0;
    let refined = run(
        cells.clone(),
        sources.clone(),
        config,
        PhotonSlabLimits::default(),
        || {
            checkpoints += 1;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(
        refined.transport.response.source_ids().collect::<Vec<_>>(),
        [900, u64::MAX]
    );
    assert_eq!(
        refined
            .transport
            .response
            .receiver_ids()
            .collect::<Vec<_>>(),
        [5, 900, u64::MAX]
    );
    let coarse = prepare_slab_response(
        "same original ballistic geometry".into(),
        DecayRadiation::Gamma,
        vec![0.0, 1e-13, 2e-13],
        cells.clone(),
        sources.clone(),
        vec![PhotonOrdinate {
            mu: 1.0,
            weight: 1.0,
        }],
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let emission: &[(u64, &[f64])] = &[(900, &[2.0, 1.0]), (u64::MAX, &[3.0, 4.0])];
    let expected = coarse.response.apply(emission, 1e-10, || Ok(())).unwrap();
    let actual = refined
        .transport
        .response
        .apply(emission, 1e-10, || Ok(()))
        .unwrap();
    for id in [900, u64::MAX, 5] {
        near(
            actual.receiver_energy_j[&id],
            expected.receiver_energy_j[&id],
            1e-13,
        );
    }
    assert_eq!(actual.receiver_energy_j[&5], 0.0);
    near(actual.escaped_j, expected.escaped_j, 1e-13);

    let mut limited = PhotonSlabLimits::default();
    limited.max_work = refined.transport.work - 1;
    assert!(matches!(
        run(cells.clone(), sources.clone(), config, limited, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    let mut limited = PhotonSlabLimits::default();
    limited.max_cells = cells.len();
    assert!(matches!(
        run(cells.clone(), sources.clone(), config, limited, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    let mut limited = PhotonSlabLimits::default();
    limited.max_values = 0;
    assert!(matches!(
        run(cells.clone(), sources.clone(), config, limited, || Ok(())),
        Err(RadiationError::Budget(_))
    ));
    assert!(matches!(
        run(
            vec![mono(6, 2.0, 0.4, 3.6)],
            vec![6],
            accuracy(1e-12, 1),
            PhotonSlabLimits::default(),
            || Ok(()),
        ),
        Err(RadiationError::Budget(_))
    ));

    let mut observed = 0;
    let cancelled = run(cells, sources, config, PhotonSlabLimits::default(), || {
        observed += 1;
        if observed == checkpoints {
            Err(RadiationError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(cancelled, Err(RadiationError::Cancelled)));
    assert_eq!(observed, checkpoints);
}

fn run_sources(
    cells: Vec<PhotonSlabCell>,
    sources: Vec<PhotonSlabSource>,
    accuracy: PhotonSlabRefinement,
    limits: PhotonSlabLimits,
    checkpoint: impl FnMut() -> Result<(), RadiationError>,
) -> Result<RefinedPhotonSlabResponse, RadiationError> {
    refine_slab_response_with_sources(
        "manufactured refined illumination".into(),
        DecayRadiation::Gamma,
        vec![0.0, 1e-13],
        cells,
        sources,
        vec![PhotonOrdinate {
            mu: 1.0,
            weight: 1.0,
        }],
        limits,
        accuracy,
        checkpoint,
    )
}

fn beam(id: u64, boundary: PhotonSlabBoundary, mu: f64) -> PhotonSlabSource {
    PhotonSlabSource {
        source_id: id,
        geometry: PhotonSourceGeometry::Collimated { boundary, mu },
    }
}

#[test]
fn g1_refined_oblique_beams_keep_full_incident_energy_and_original_receivers() {
    let sources = vec![
        beam(u64::MAX, PhotonSlabBoundary::Left, 0.5),
        beam(8, PhotonSlabBoundary::Right, 0.25),
    ];
    let mut limits = PhotonSlabLimits::default();
    // A spatial split must not invent another boundary source or charge its
    // source budget as if it were an internal emitter in every sublayer.
    limits.response.max_sources = 2;
    let result = run_sources(
        vec![mono(900, 0.8, 0.7, 0.0), mono(3, 0.4, 1.1, 0.0)],
        sources,
        accuracy(1e-12, 1),
        limits,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(result.subdivisions_per_layer, 2);
    assert_eq!(
        result.transport.response.source_ids().collect::<Vec<_>>(),
        [8, u64::MAX]
    );
    assert_eq!(
        result.transport.response.receiver_ids().collect::<Vec<_>>(),
        [3, 900]
    );
    let energy = result
        .transport
        .response
        .apply(&[(u64::MAX, &[7.0]), (8, &[4.0])], 1e-10, || Ok(()))
        .unwrap();
    let left_a = 7.0 * -(-1.12_f64).exp_m1();
    let left_b = 7.0 * (-1.12_f64).exp() * -(-0.88_f64).exp_m1();
    let right_b = 4.0 * -(-1.76_f64).exp_m1();
    let right_a = 4.0 * (-1.76_f64).exp() * -(-2.24_f64).exp_m1();
    near(energy.receiver_energy_j[&900], left_a + right_a, 5e-14);
    near(energy.receiver_energy_j[&3], left_b + right_b, 5e-14);
    near(
        energy.escaped_j,
        7.0 * (-2.0_f64).exp() + 4.0 * (-4.0_f64).exp(),
        5e-14,
    );
    near(energy.deposited_j + energy.escaped_j, 11.0, 5e-14);
}

#[test]
fn g3_incident_scattering_refines_to_independent_two_stream_boundary_solution() {
    let solve = |tolerance| {
        run_sources(
            vec![mono(21, 2.0, 0.4, 3.6)],
            vec![beam(100, PhotonSlabBoundary::Left, 1.0)],
            accuracy(tolerance, 7),
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap()
    };
    let loose = solve(1e-2);
    let tight = solve(8e-4);
    // Independent two-stream ODE: I+'=-c I+ + d I-,
    // I-'=-d I+ + c I-, with I+(0)=1 and I-(L)=0.
    // exp(M L) gives R and T directly, without collision probabilities.
    let k = (0.4_f64 * 4.0).sqrt();
    let z = 2.0 * k;
    let c = 0.4 + 0.5 * 3.6;
    let d = 0.5 * 3.6;
    let reflected = d * z.tanh() / (k + c * z.tanh());
    let transmitted = 1.0 / (z.cosh() + c / k * z.sinh());
    let exact = reflected + transmitted;
    let escape = |result: &RefinedPhotonSlabResponse| {
        result
            .transport
            .response
            .apply(&[(100, &[1.0])], 1e-10, || Ok(()))
            .unwrap()
            .escaped_j
    };
    assert!(loose.levels >= 2);
    assert!(tight.levels > loose.levels);
    assert!((escape(&tight) - exact).abs() < (escape(&loose) - exact).abs());
    near(escape(&tight), exact, 6e-4);
}

#[test]
fn g3_mixed_internal_and_boundary_sources_preserve_superposition_and_cancellation() {
    let cells = vec![mono(81, 0.4, 0.8, 0.0), mono(17, 0.7, 0.3, 0.0)];
    let sources = vec![
        beam(600, PhotonSlabBoundary::Right, 0.6),
        PhotonSlabSource {
            source_id: 20,
            geometry: PhotonSourceGeometry::UniformLayer { cell_id: 81 },
        },
        // Same physical emitting layer with a distinct source identity.
        PhotonSlabSource {
            source_id: 7,
            geometry: PhotonSourceGeometry::UniformLayer { cell_id: 81 },
        },
    ];
    let config = accuracy(1e-12, 1);
    let mut checkpoints = 0;
    let refined = run_sources(
        cells.clone(),
        sources.clone(),
        config,
        PhotonSlabLimits::default(),
        || {
            checkpoints += 1;
            Ok(())
        },
    )
    .unwrap();
    let coarse = prepare_slab_response_with_sources(
        "same mixed source slab".into(),
        DecayRadiation::Gamma,
        vec![0.0, 1e-13],
        cells.clone(),
        sources.clone(),
        vec![PhotonOrdinate {
            mu: 1.0,
            weight: 1.0,
        }],
        PhotonSlabLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let input: &[(u64, &[f64])] = &[(600, &[2.0]), (20, &[3.0]), (7, &[5.0])];
    let expected = coarse.response.apply(input, 1e-10, || Ok(())).unwrap();
    let actual = refined
        .transport
        .response
        .apply(input, 1e-10, || Ok(()))
        .unwrap();
    for id in [81, 17] {
        near(
            actual.receiver_energy_j[&id],
            expected.receiver_energy_j[&id],
            5e-14,
        );
    }
    near(actual.escaped_j, expected.escaped_j, 5e-14);
    let mut observed = 0;
    assert!(matches!(
        run_sources(cells, sources, config, PhotonSlabLimits::default(), || {
            observed += 1;
            if observed == checkpoints {
                Err(RadiationError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(RadiationError::Cancelled)
    ));
}
