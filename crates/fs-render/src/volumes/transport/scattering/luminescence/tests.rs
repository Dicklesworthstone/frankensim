use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{
    CellPhaseFunction, MassScatteringSpectrum, PhaseFunction,
};
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
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

fn geometry<'a>(cx: &Cx<'_>, density: &'a [f64], absorption: f64) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([density.len(), 1, 1], density, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(absorption).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: density.len(),
            max_steps_per_segment: 32,
        },
    )
    .unwrap()
}

fn medium<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    absorption: f64,
    scattering: f64,
) -> ScatteringGrid<'a> {
    ScatteringGrid::try_new(
        cx,
        geometry(cx, density, absorption),
        MassScatteringSpectrum::constant(scattering).unwrap(),
        PhaseFunction::Isotropic,
    )
    .unwrap()
}

fn config(histories: u32, collisions: u32) -> LuminescenceTransportConfig {
    LuminescenceTransportConfig {
        histories_per_emitter_lane: histories,
        max_collisions_per_history: collisions,
        max_total_segments: u64::from(histories) * u64::from(collisions) * 32,
        max_material_cells: 8,
        max_emitter_cells: 8,
        max_absorption_links: 64,
        seed: 1821,
        source_id: 43,
    }
}

#[test]
fn transparent_material_and_zero_density_cells_emit_with_exact_escape() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.0, 0.0, 2.0], 0.0, 0.0);
        let result = medium
            .estimate_luminescence_transfer(cx, &[500.0, 650.0], &[2, 0, 1], config(32, 1))
            .unwrap();
        result
            .validate_with_checkpoint(&mut || Ok(()), 0.0)
            .unwrap();
        assert_eq!(result.medium_identity, medium.content_identity());
        assert_eq!(
            result.absorption_identity,
            medium.absorption_grid().content_identity()
        );
        assert_eq!(result.material_cells, [0, 2]);
        assert_eq!(result.grid_cell_count, 3);
        assert_eq!(result.wavelengths_nm(), &[500.0, 650.0]);
        assert_eq!(result.traced_segments, 3 * 2 * 32);
        for (index, row) in result.rows.iter().enumerate() {
            assert_eq!(row.emitter_cell, index);
            assert!(row.absorption.is_empty());
            assert_eq!(row.absorbed, FractionEstimate::zero());
            assert_eq!(row.unresolved, FractionEstimate::zero());
            assert_eq!(row.escape.mean, [1.0; 2]);
            assert_eq!(row.escape.standard_error, [0.0; 2]);
        }
        assert_eq!(
            result,
            medium
                .estimate_luminescence_transfer(cx, &[500.0, 650.0], &[0, 1, 2], config(32, 1))
                .unwrap()
        );
        // The native transport packet is not limited by the separate emission
        // source's table size. Duplicate wavelengths remain independent lanes.
        let wide = medium
            .estimate_luminescence_transfer(cx, &[550.0; 17], &[1], config(2, 1))
            .unwrap();
        wide.validate_with_checkpoint(&mut || Ok(()), 0.0).unwrap();
        assert_eq!(wide.rows[0].escape.mean, [1.0; 17]);
    });
}

#[test]
fn opaque_source_returns_energy_to_itself_without_lte_capacity_or_symmetrization() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.0], 1e12, 0.0);
        let result = medium
            .estimate_luminescence_transfer(cx, &[550.0], &[0], config(64, 1))
            .unwrap();
        result
            .validate_with_checkpoint(&mut || Ok(()), 0.0)
            .unwrap();
        let row = &result.rows[0];
        assert_eq!(row.absorption.len(), 1);
        assert_eq!(row.absorption[0].cell_index, 0);
        assert_eq!(row.absorption[0].fraction.mean, [1.0]);
        assert_eq!(row.absorbed.mean, [1.0]);
        assert_eq!(row.escape.mean, [0.0]);
        assert_eq!(row.unresolved.mean, [0.0]);
        // The row acts on radiant joules directly, without multiplying by
        // volume, absorptivity, solid angle, or a spectral-bin width again.
        let emitted_j = 7.25;
        assert_eq!(emitted_j * row.absorption[0].fraction.mean[0], emitted_j);
    });
}

// Independent volume/solid-angle quadrature using finite-box chords and the
// analytic Beer law, not the production DDA or sampled extinction histories.
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

fn emitted_point_direction(index: usize) -> ([f64; 3], [f64; 3]) {
    let p = [halton(index, 2), halton(index, 3), halton(index, 5)];
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

fn box_escape(absorption: f64, samples: usize) -> f64 {
    (1..=samples)
        .map(|i| {
            let (p, d) = emitted_point_direction(i);
            (-absorption * exit_distance(p, d, 1.0)).exp()
        })
        .sum::<f64>()
        / samples as f64
}

fn neighboring_capture(absorption: f64, samples: usize) -> f64 {
    (1..=samples)
        .map(|i| {
            let (p, d) = emitted_point_direction(i);
            if d[0] <= 0.0 {
                return 0.0;
            }
            let interface = (1.0 - p[0]) / d[0];
            let exit = exit_distance(p, d, 2.0);
            if interface >= exit {
                0.0
            } else {
                -(-absorption * (exit - interface)).exp_m1()
            }
        })
        .sum::<f64>()
        / samples as f64
}

#[test]
fn homogeneous_source_matches_independent_uniform_box_integral() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.0], 1.3, 0.0);
        let result = medium
            .estimate_luminescence_transfer(cx, &[550.0], &[0], config(12_000, 1))
            .unwrap();
        result
            .validate_with_checkpoint(&mut || Ok(()), 2e-12)
            .unwrap();
        let reference = box_escape(1.3, 65_536);
        assert!((reference - box_escape(1.3, 32_768)).abs() < 1e-3);
        let estimate = &result.rows[0].escape;
        assert!(
            (estimate.mean[0] - reference).abs() < 6.0 * estimate.standard_error[0] + 1e-3,
            "observed {:?}, independent {reference}",
            estimate
        );
    });
}

#[test]
fn cell_local_spectra_retain_directed_neighbor_absorption_and_replay() {
    with_cx(|cx, _| {
        let geometry = geometry(cx, &[1.0, 1.0], 0.0);
        let field = CellScatteringGrid::try_new(
            cx,
            &geometry,
            [500.0, 650.0],
            &[[0.0, 4.0], [4.0, 0.0]],
            &[[0.0; 2]; 2],
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (380.0, 780.0),
            ContentHash([7; 32]),
        )
        .unwrap();
        let options = config(8192, 1);
        let result = field
            .estimate_luminescence_transfer(cx, &[1, 0], options)
            .unwrap();
        result
            .validate_with_checkpoint(&mut || Ok(()), 2e-12)
            .unwrap();
        assert_eq!(result.absorption_identity, geometry.content_identity());
        assert_eq!(result.medium_identity, field.content_identity());
        assert_eq!(result.material_cells, [0, 1]);
        let row = &result.rows[0];
        assert_eq!(row.absorption.len(), 2);
        assert_eq!(row.absorption[0].cell_index, 0);
        assert_eq!(row.absorption[1].cell_index, 1);
        assert_eq!(row.absorption[0].fraction.mean[0], 0.0);
        assert_eq!(row.absorption[1].fraction.mean[1], 0.0);
        let forward = &row.absorption[1].fraction;
        let reverse = &result.rows[1].absorption[0].fraction;
        assert!(forward.mean[0] > 0.05);
        assert_eq!(reverse.mean[0], 0.0);
        assert!(reverse.mean[1] > 0.05);
        let reference = neighboring_capture(4.0, 65_536);
        assert!((reference - neighboring_capture(4.0, 32_768)).abs() < 1e-3);
        assert!((forward.mean[0] - reference).abs() < 6.0 * forward.standard_error[0] + 1e-3);
        assert_eq!(
            result,
            field
                .estimate_luminescence_transfer(cx, &[0, 1], options)
                .unwrap()
        );
        let unordered_reports: [fn(&mut LuminescenceTransferOperator<2>); 4] = [
            |r| r.rows.swap(0, 1),
            |r| r.material_cells.swap(0, 1),
            |r| r.rows[0].absorption.swap(0, 1),
            |r| r.rows[0].absorption[1].cell_index = 0,
        ];
        for reorder in unordered_reports {
            let mut changed = result.clone();
            reorder(&mut changed);
            assert_eq!(
                changed.validate_with_checkpoint(&mut || Ok(()), 2e-12),
                Err(VolumeTransportError::Grid)
            );
        }
        let changed_source = field
            .estimate_luminescence_transfer(
                cx,
                &[0, 1],
                LuminescenceTransportConfig {
                    source_id: options.source_id + 1,
                    ..options
                },
            )
            .unwrap();
        assert_ne!(
            result.rows[0].escape.mean,
            changed_source.rows[0].escape.mean
        );
    });
}

#[test]
fn elastic_scattering_tail_is_radiation_and_never_absorption() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.0], 0.0, 8.0);
        let result = medium
            .estimate_luminescence_transfer(cx, &[550.0], &[0], config(2048, 2))
            .unwrap();
        result
            .validate_with_checkpoint(&mut || Ok(()), 2e-12)
            .unwrap();
        let row = &result.rows[0];
        assert!(row.absorption.is_empty());
        assert_eq!(row.absorbed.mean, [0.0]);
        assert!(row.escape.mean[0] > 0.0);
        assert!(row.unresolved.mean[0] > 0.5);
        assert!((row.escape.mean[0] + row.unresolved.mean[0] - 1.0).abs() < 2e-12);
        assert!(result.traced_segments > 2048);
        assert!(result.traced_segments <= 4096);
    });
}

#[test]
fn public_report_admission_preserves_explicit_fraction_allowance() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.0], 1e12, 0.0);
        let report = medium
            .estimate_luminescence_transfer(cx, &[550.0], &[0], config(16, 1))
            .unwrap();
        let corruptions: [fn(&mut LuminescenceTransferOperator<1>); 16] = [
            |r| r.rows[0].absorbed.mean[0] = f64::NAN,
            |r| r.rows[0].escape.mean[0] = -1e-20,
            |r| r.rows[0].unresolved.standard_error[0] = f64::INFINITY,
            |r| r.rows[0].absorption[0].fraction.standard_error[0] = -1.0,
            |r| r.rows[0].absorption[0].fraction.mean[0] = 0.5,
            |r| r.rows[0].escape.mean[0] = 0.25,
            |r| r.rows[0].absorption[0].cell_index = 1,
            |r| r.rows[0].emitter_cell = 1,
            |r| r.rows[0].precision_limited_histories[0] = 17,
            |r| r.material_cells.clear(),
            |r| r.material_cells.push(0),
            |r| r.medium_identity = ContentHash([0; 32]),
            |r| r.absorption_identity = ContentHash([0; 32]),
            |r| r.wavelengths_nm[0] = 200.0,
            |r| r.traced_segments = 15,
            |r| r.config.max_absorption_links = 0,
        ];
        for corrupt in corruptions {
            let mut changed = report.clone();
            corrupt(&mut changed);
            assert!(
                changed
                    .validate_with_checkpoint(&mut || Ok(()), 1e-12)
                    .is_err()
            );
        }
        for allowance in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(
                report
                    .validate_with_checkpoint(&mut || Ok(()), allowance)
                    .is_err()
            );
        }
        let mut rounded = report.clone();
        rounded.rows[0].absorbed.mean[0] += 4e-10;
        rounded.rows[0].absorption[0].fraction.mean[0] += 4e-10;
        rounded
            .validate_with_checkpoint(&mut || Ok(()), 1e-9)
            .unwrap();
        assert!(
            rounded
                .validate_with_checkpoint(&mut || Ok(()), 1e-12)
                .is_err()
        );
        assert_eq!(
            report.validate_with_checkpoint(&mut || Err(VolumeTransportError::Cancelled), 0.0),
            Err(VolumeTransportError::Cancelled)
        );
    });
}

#[test]
fn invalid_sources_storage_work_and_cancellation_refuse_without_partial_output() {
    with_cx(|cx, gate| {
        let medium = medium(cx, &[1.0], 1e12, 0.0);
        let options = config(16, 1);
        for cells in [&[0, 0][..], &[1][..]] {
            assert_eq!(
                medium.estimate_luminescence_transfer(cx, &[550.0], cells, options),
                Err(VolumeTransportError::Grid)
            );
        }
        for invalid in [
            LuminescenceTransportConfig {
                histories_per_emitter_lane: 1,
                ..options
            },
            LuminescenceTransportConfig {
                max_collisions_per_history: 0,
                ..options
            },
            LuminescenceTransportConfig {
                max_total_segments: 15,
                ..options
            },
            LuminescenceTransportConfig {
                max_material_cells: 0,
                ..options
            },
            LuminescenceTransportConfig {
                max_emitter_cells: 0,
                ..options
            },
            LuminescenceTransportConfig {
                max_absorption_links: 0,
                ..options
            },
        ] {
            assert_eq!(
                medium.estimate_luminescence_transfer(cx, &[550.0], &[0], invalid),
                Err(VolumeTransportError::Budget)
            );
        }
        let two_cells = ScatteringGrid::try_new(
            cx,
            geometry(cx, &[1.0, 1.0], 1e12),
            MassScatteringSpectrum::constant(0.0).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        // The directed-link allowance is global, not renewed for each emitter.
        assert_eq!(
            two_cells.estimate_luminescence_transfer(
                cx,
                &[550.0],
                &[0, 1],
                LuminescenceTransportConfig {
                    max_absorption_links: 1,
                    ..options
                },
            ),
            Err(VolumeTransportError::Budget)
        );
        let empty = medium
            .estimate_luminescence_transfer(cx, &[550.0], &[], options)
            .unwrap();
        assert!(empty.rows.is_empty());
        assert_eq!(empty.material_cells, [0]);
        assert_eq!(empty.traced_segments, 0);
        empty.validate_with_checkpoint(&mut || Ok(()), 0.0).unwrap();
        let field = CellScatteringGrid::try_new(
            cx,
            medium.absorption_grid(),
            [550.0],
            &[[1.0]],
            &[[0.0]],
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (380.0, 780.0),
            ContentHash([4; 32]),
        )
        .unwrap();
        let scattering = CellScatteringGrid::try_new(
            cx,
            medium.absorption_grid(),
            [550.0],
            &[[0.0]],
            &[[8.0]],
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (380.0, 780.0),
            ContentHash([5; 32]),
        )
        .unwrap();
        // The initial work lower bound fits, but actual continuations exhaust
        // the global free-flight allowance. No incomplete operator is returned.
        assert_eq!(
            scattering.estimate_luminescence_transfer(
                cx,
                &[0],
                LuminescenceTransportConfig {
                    max_collisions_per_history: 2,
                    max_total_segments: 16,
                    ..options
                },
            ),
            Err(VolumeTransportError::Budget)
        );
        let mut checks = 0;
        let cancelled = field.estimate_luminescence_transfer_with_checkpoint(
            &mut || {
                checks += 1;
                if checks == 20 {
                    Err(VolumeTransportError::Cancelled)
                } else {
                    Ok(())
                }
            },
            &[0],
            options,
        );
        assert_eq!(cancelled, Err(VolumeTransportError::Cancelled));
        gate.request();
        assert_eq!(
            medium.estimate_luminescence_transfer(cx, &[550.0], &[0], options),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
