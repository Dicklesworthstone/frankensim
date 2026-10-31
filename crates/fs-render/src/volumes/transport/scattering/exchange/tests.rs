use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{
    CellPhaseFunction, MassScatteringSpectrum, PhaseFunction,
};
use crate::volumes::transport::{AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_exec::{Budget, CancelGate, ExecMode};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 17,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}
fn geometry<'a>(cx: &Cx<'_>, density: &'a [f64], cell: [f64; 3]) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([density.len(), 1, 1], density, [0.0; 3], cell),
        MassAbsorptionSpectrum::constant(0.0).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: density.len(),
            max_steps_per_segment: 32,
        },
    )
    .unwrap()
}
fn config(histories: u32, collisions: u32) -> RadiativeExchangeConfig {
    RadiativeExchangeConfig {
        histories_per_emitter_lane: histories,
        max_collisions_per_history: collisions,
        max_total_segments: u64::from(histories) * u64::from(collisions) * 8,
        max_material_cells: 8,
        max_emitter_cells: 8,
        max_pairs: 28,
        seed: 1821,
        source_id: 43,
    }
}

// Independent deterministic integration over uniformly emitting volume and
// solid angle. A direct finite-box chord and exp attenuation are used, rather
// than the production DDA, exponential collision draw or implicit capture.
fn halton(mut index: usize, base: usize) -> f64 {
    let mut value = 0.0;
    let mut weight = 1.0;
    while index > 0 {
        weight /= base as f64;
        value += weight * (index % base) as f64;
        index /= base;
    }
    value
}
fn sample_cube(index: usize, offset: f64) -> ([f64; 3], [f64; 3]) {
    let p = [
        offset + halton(index, 2),
        halton(index, 3),
        halton(index, 5),
    ];
    let z = 1.0 - 2.0 * halton(index, 7);
    let r = (1.0 - z * z).sqrt();
    let phi = 2.0 * std::f64::consts::PI * halton(index, 11);
    (p, [r * phi.cos(), r * phi.sin(), z])
}
fn exit_distance(p: [f64; 3], d: [f64; 3], width: f64) -> f64 {
    (0..3)
        .map(|axis| {
            let upper = if axis == 0 { width } else { 1.0 };
            if d[axis] > 0.0 {
                (upper - p[axis]) / d[axis]
            } else if d[axis] < 0.0 {
                -p[axis] / d[axis]
            } else {
                f64::INFINITY
            }
        })
        .fold(f64::INFINITY, f64::min)
}
fn cube_escape(absorption: f64, samples: usize) -> f64 {
    (1..=samples)
        .map(|i| {
            let (p, d) = sample_cube(i, 0.0);
            (-absorption * exit_distance(p, d, 1.0)).exp()
        })
        .sum::<f64>()
        / samples as f64
}
fn layer_capture(source: usize, a: [f64; 2], samples: usize) -> f64 {
    (1..=samples)
        .map(|i| {
            let (p, d) = sample_cube(i, source as f64);
            if (source == 0 && d[0] <= 0.0) || (source == 1 && d[0] >= 0.0) {
                return 0.0;
            }
            let interface = (1.0 - p[0]) / d[0];
            let exit = exit_distance(p, d, 2.0);
            if interface >= exit {
                return 0.0;
            }
            (-a[source] * interface).exp() * -(-a[1 - source] * (exit - interface)).exp_m1()
        })
        .sum::<f64>()
        / samples as f64
}

#[test]
fn g1_pure_absorption_emission_escape_matches_independent_finite_cube_integral() {
    with_cx(|cx, _| {
        let base = geometry(cx, &[1.0], [1.0; 3]);
        for absorption in [0.01, 0.7, 6.0] {
            let a = [[absorption]];
            let field = CellScatteringGrid::try_new(
                cx,
                &base,
                [1500.0],
                &a,
                &[[0.0]],
                CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
                (1000.0, 2000.0),
                ContentHash([7; 32]),
            )
            .unwrap();
            let r = field
                .estimate_thermal_exchange(cx, config(50_000, 1))
                .unwrap();
            r.validate_with_checkpoint(&mut || Ok(()), 2e-11).unwrap();
            assert!(r.pairs.is_empty());
            let cell = &r.cells[0];
            let capacity = 4.0 * std::f64::consts::PI * absorption;
            assert_eq!(cell.emission_capacity_m2_sr[0], capacity);
            assert_eq!(cell.sampled_histories, [50_000]);
            assert_eq!(cell.unresolved.mean_m2_sr, [0.0]);
            // Total moments observe every zero history individually; sparse
            // cell moments fold missing histories into batches. Both estimate
            // this same Bernoulli capture distribution, with different final
            // binary64 reduction roundoff rather than bitwise equality.
            for (total, local) in cell
                .absorbed
                .mean_m2_sr
                .iter()
                .chain(cell.absorbed.standard_error_m2_sr.iter())
                .zip(
                    cell.self_return
                        .mean_m2_sr
                        .iter()
                        .chain(cell.self_return.standard_error_m2_sr.iter()),
                )
            {
                let roundoff = 64.0 * f64::EPSILON * total.abs().max(local.abs());
                assert!(
                    (total - local).abs() <= roundoff,
                    "total {total:e}, self return {local:e}, coefficient roundoff {roundoff:e}"
                );
            }
            let reference = cube_escape(absorption, 131_072);
            let coarse = cube_escape(absorption, 65_536);
            assert!((reference - coarse).abs() < 3e-4);
            let observed = cell.escape.mean_m2_sr[0] / capacity;
            let se = cell.escape.standard_error_m2_sr[0] / capacity;
            assert!(
                (observed - reference).abs() < 6.0 * se + 3e-4,
                "alpha {absorption}: native {observed}, independent {reference}, SE {se}"
            );
            if absorption == 0.01 {
                assert!(observed > 0.98);
            }
        }
    });
}

#[test]
fn g1_two_material_pair_is_reciprocal_positive_and_matches_volume_integrals() {
    with_cx(|cx, _| {
        let base = geometry(cx, &[1.0, 1.0], [1.0; 3]);
        let a = [[0.2], [2.0]];
        let field = CellScatteringGrid::try_new(
            cx,
            &base,
            [1500.0],
            &a,
            &[[0.0]; 2],
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (1000.0, 2000.0),
            ContentHash([8; 32]),
        )
        .unwrap();
        let r = field
            .estimate_thermal_exchange(cx, config(90_000, 1))
            .unwrap();
        r.validate_with_checkpoint(&mut || Ok(()), 3e-11).unwrap();
        assert_eq!(r.pairs.len(), 1);
        let p = &r.pairs[0];
        assert_eq!((p.cell_a, p.cell_b), (0, 1));
        let forward = 4.0 * std::f64::consts::PI * 0.2 * layer_capture(0, [0.2, 2.0], 131_072);
        let reverse = 4.0 * std::f64::consts::PI * 2.0 * layer_capture(1, [0.2, 2.0], 131_072);
        assert!((forward - reverse).abs() < 2e-3);
        let reference = 0.5 * (forward + reverse);
        assert!(
            (p.exchange.mean_m2_sr[0] - reference).abs()
                < 6.0 * p.exchange.standard_error_m2_sr[0] + 2e-3
        );
        assert!(
            p.reciprocity_residual_m2_sr[0].abs()
                < 12.0 * p.exchange.standard_error_m2_sr[0] + 2e-3
        );
        let cold_to_hot = p.exchange.mean_m2_sr[0] * (2.0 - 5.0);
        let hot_to_cold = -cold_to_hot;
        assert_eq!(cold_to_hot + hot_to_cold, 0.0);
        assert!(cold_to_hot < 0.0);
        assert_eq!(p.exchange.mean_m2_sr[0] * (5.0 - 5.0), 0.0);
    });
}

#[test]
fn g3_optical_depth_preserving_size_change_scales_exchange_by_area_and_replays() {
    with_cx(|cx, _| {
        let make = |length, a| {
            let base = geometry(cx, &[1.0], [length; 3]);
            let medium = ScatteringGrid::try_new(
                cx,
                base,
                MassScatteringSpectrum::constant(0.0).unwrap(),
                PhaseFunction::Isotropic,
            )
            .unwrap();
            let field_a = [[a]];
            let field = CellScatteringGrid::try_new(
                cx,
                medium.absorption_grid(),
                [550.0],
                &field_a,
                &[[0.0]],
                CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
                (380.0, 780.0),
                ContentHash([9; 32]),
            )
            .unwrap();
            field
                .estimate_thermal_exchange(cx, config(4096, 1))
                .unwrap()
        };
        let small = make(1.0, 0.8);
        let large = make(2.0, 0.4);
        assert_eq!(small, make(1.0, 0.8));
        assert_eq!(large.cell_volume_m3, 8.0 * small.cell_volume_m3);
        assert_eq!(
            large.cells[0].emission_capacity_m2_sr[0],
            4.0 * small.cells[0].emission_capacity_m2_sr[0]
        );
        assert_eq!(
            large.cells[0].escape.mean_m2_sr[0],
            4.0 * small.cells[0].escape.mean_m2_sr[0]
        );
        assert_eq!(
            large.cells[0].self_return.mean_m2_sr[0],
            4.0 * small.cells[0].self_return.mean_m2_sr[0]
        );
    });
}

#[test]
fn g3_scattering_tail_stays_unresolved_and_conservative_scatter_does_not_emit() {
    with_cx(|cx, _| {
        let base = geometry(cx, &[1.0], [1.0; 3]);
        let phase = PhaseFunction::henyey_greenstein(0.7).unwrap();
        let a = [[0.2]];
        let s = [[3.0]];
        let field = CellScatteringGrid::try_new(
            cx,
            &base,
            [550.0],
            &a,
            &s,
            CellPhaseFunction::Uniform(&phase),
            (380.0, 780.0),
            ContentHash([10; 32]),
        )
        .unwrap();
        let shallow = field
            .estimate_thermal_exchange(cx, config(8000, 1))
            .unwrap();
        let deep = field
            .estimate_thermal_exchange(cx, config(8000, 64))
            .unwrap();
        shallow
            .validate_with_checkpoint(&mut || Ok(()), 1e-11)
            .unwrap();
        deep.validate_with_checkpoint(&mut || Ok(()), 1e-11)
            .unwrap();
        assert!(shallow.cells[0].unresolved.mean_m2_sr[0] > 0.1);
        assert!(deep.cells[0].unresolved.mean_m2_sr[0] < 1e-8);
        assert!(
            deep.cells[0].self_return.mean_m2_sr[0] > shallow.cells[0].self_return.mean_m2_sr[0]
        );
        let conservative = CellScatteringGrid::try_new(
            cx,
            &base,
            [550.0],
            &[[0.0]],
            &s,
            CellPhaseFunction::Uniform(&phase),
            (380.0, 780.0),
            ContentHash([11; 32]),
        )
        .unwrap();
        let empty = conservative
            .estimate_thermal_exchange(
                cx,
                RadiativeExchangeConfig {
                    max_total_segments: 0,
                    max_emitter_cells: 0,
                    max_pairs: 0,
                    ..config(2, 1)
                },
            )
            .unwrap();
        assert!(empty.cells.is_empty() && empty.pairs.is_empty());
        assert_eq!(empty.traced_segments, 0);
        empty.validate_with_checkpoint(&mut || Ok(()), 0.0).unwrap();
    });
}

#[test]
fn g0_work_storage_corruption_and_callback_errors_refuse_without_partial_operator() {
    with_cx(|cx, _| {
        let base = geometry(cx, &[1.0, 1.0], [1.0; 3]);
        let a = [[0.8]; 2];
        let s = [[0.2]; 2];
        let field = CellScatteringGrid::try_new(
            cx,
            &base,
            [550.0],
            &a,
            &s,
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (380.0, 780.0),
            ContentHash([12; 32]),
        )
        .unwrap();
        for invalid in [
            RadiativeExchangeConfig {
                max_emitter_cells: 1,
                ..config(128, 16)
            },
            RadiativeExchangeConfig {
                max_total_segments: 255,
                ..config(128, 16)
            },
            RadiativeExchangeConfig {
                max_pairs: 0,
                ..config(128, 16)
            },
            RadiativeExchangeConfig {
                histories_per_emitter_lane: 1,
                ..config(128, 16)
            },
            RadiativeExchangeConfig {
                max_collisions_per_history: 0,
                ..config(128, 16)
            },
        ] {
            assert!(matches!(
                field.estimate_thermal_exchange(cx, invalid),
                Err(VolumeTransportError::Budget)
            ));
        }
        let mut checks = 0;
        let cancelled = field.estimate_thermal_exchange_with_checkpoint(
            &mut || {
                checks += 1;
                if checks == 35 {
                    Err(VolumeTransportError::Cancelled)
                } else {
                    Ok(())
                }
            },
            config(128, 16),
        );
        assert!(matches!(cancelled, Err(VolumeTransportError::Cancelled)));
        assert!(matches!(
            field.estimate_thermal_exchange_with_checkpoint(
                &mut || Err(VolumeTransportError::Temperature),
                config(128, 16)
            ),
            Err(VolumeTransportError::Temperature)
        ));
        let valid = field
            .estimate_thermal_exchange(cx, config(128, 16))
            .unwrap();
        valid
            .validate_with_checkpoint(&mut || Ok(()), 1e-11)
            .unwrap();
        let mut corrupt = valid.clone();
        corrupt.cells[0].escape.mean_m2_sr[0] += 0.1;
        assert!(
            corrupt
                .validate_with_checkpoint(&mut || Ok(()), 1e-11)
                .is_err()
        );
        let mut corrupt = valid.clone();
        corrupt.pairs[0].exchange.mean_m2_sr[0] = -1.0;
        assert!(
            corrupt
                .validate_with_checkpoint(&mut || Ok(()), 1e-11)
                .is_err()
        );
        let mut corrupt = valid.clone();
        corrupt.cells.swap(0, 1);
        assert!(
            corrupt
                .validate_with_checkpoint(&mut || Ok(()), 1e-11)
                .is_err()
        );
        assert_eq!(
            valid,
            field
                .estimate_thermal_exchange(cx, config(128, 16))
                .unwrap()
        );
    });
}

#[test]
fn g0_reciprocal_averaging_preserves_representable_subnormals_and_refuses_lost_uncertainty() {
    let minimum = f64::from_bits(1);
    let estimate = ExchangeCoefficientEstimate {
        mean_m2_sr: [minimum],
        standard_error_m2_sr: [minimum],
    };
    let mut pairs = PairTallies::new(1);
    pairs.insert(0, 1, estimate).unwrap();
    pairs.insert(1, 0, estimate).unwrap();
    let result = pairs.finish(&mut || Ok(())).unwrap();
    assert_eq!(result[0].exchange.mean_m2_sr, [minimum]);
    assert_eq!(result[0].exchange.standard_error_m2_sr, [minimum]);
    let mut pairs = PairTallies::new(1);
    pairs.insert(0, 1, estimate).unwrap();
    assert!(matches!(
        pairs.finish(&mut || Ok(())),
        Err(VolumeTransportError::Unrepresentable)
    ));
}

#[test]
fn g0_material_domains_include_pure_scatterers_and_transparent_cells_without_photon_work() {
    with_cx(|cx, _| {
        let density = [0.0, 1.0, 1.0, 1.0];
        let base = geometry(cx, &density, [1.0; 3]);
        let a = [[10.0], [0.0], [0.0], [0.7]];
        let s = [[8.0], [1.0], [0.0], [0.0]];
        let field = CellScatteringGrid::try_new(
            cx,
            &base,
            [550.0],
            &a,
            &s,
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (380.0, 780.0),
            ContentHash([23; 32]),
        )
        .unwrap();
        let r = field
            .estimate_thermal_exchange(cx, config(256, 32))
            .unwrap();
        assert_eq!(r.material_cells, vec![1, 2, 3]);
        assert_eq!(r.cells.len(), 1);
        assert_eq!(r.cells[0].cell_index, 3);
        assert_eq!(r.cells[0].sampled_histories, [256]);
        r.validate_with_checkpoint(&mut || Ok(()), 1e-11).unwrap();
        assert!(matches!(
            field.estimate_thermal_exchange(
                cx,
                RadiativeExchangeConfig {
                    max_material_cells: 2,
                    ..config(256, 32)
                }
            ),
            Err(VolumeTransportError::Budget)
        ));
        let mut corrupt = r.clone();
        corrupt.material_cells.swap(0, 1);
        assert!(
            corrupt
                .validate_with_checkpoint(&mut || Ok(()), 1e-11)
                .is_err()
        );
        let mut corrupt = r;
        corrupt.material_cells.pop();
        assert!(
            corrupt
                .validate_with_checkpoint(&mut || Ok(()), 1e-11)
                .is_err()
        );
        let vacuum = geometry(cx, &[0.0], [1.0; 3]);
        let field = CellScatteringGrid::try_new(
            cx,
            &vacuum,
            [550.0],
            &[[1.0]],
            &[[1.0]],
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (380.0, 780.0),
            ContentHash([24; 32]),
        )
        .unwrap();
        let empty = field
            .estimate_thermal_exchange(
                cx,
                RadiativeExchangeConfig {
                    max_material_cells: 0,
                    max_emitter_cells: 0,
                    max_pairs: 0,
                    max_total_segments: 0,
                    ..config(2, 1)
                },
            )
            .unwrap();
        assert!(empty.material_cells.is_empty() && empty.cells.is_empty());
        assert_eq!(empty.traced_segments, 0);
        empty.validate_with_checkpoint(&mut || Ok(()), 0.0).unwrap();
    });
}
