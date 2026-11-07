//! Continued wake shedding without truncating aged physical vorticity.

use fs_vpm::farfield3d::{FarField, FarRetirementLimits, WakeCoreEvolutionMode, hybrid_kinematics};
use fs_vpm::filament3d::{FilamentWake, SelfInducedShed, WakeRateCertificate, WakeSegment};

const DT: f64 = 0.05;
const LIMITS: FarRetirementLimits = FarRetirementLimits {
    max_cells: 128,
    max_segments: 4096,
};

fn sheet(rows: usize, cap: usize) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 1.0 / DT,
            n_stations: 3,
            max_rows: cap,
        },
        [-0.9, -0.3, 0.3, 0.9]
            .map(|y| [0.0, y, 1.0 + 0.03 * y * y])
            .to_vec(),
    )
    .unwrap();
    for r in 0..rows {
        let t = r as f64;
        wake.shed(&[0.1 + 0.02 * t, -0.15 + 0.01 * t, 0.08], [-0.2, 0.0, 0.0])
            .unwrap();
        wake.rows.last_mut().unwrap().core2_m2 = 0.04 + 0.003 * t;
    }
    wake
}

fn empty(mode: WakeCoreEvolutionMode) -> FarField {
    FarField {
        mode,
        shed_dt_s: DT,
        cells: Vec::new(),
    }
}

fn sources(near: &FilamentWake, far: &FarField) -> Vec<WakeSegment> {
    let mut all = Vec::new();
    for cell in &far.cells {
        for (&(start, end, circulation), &core) in cell.segs.iter().zip(&cell.source_core2_m2) {
            all.push(WakeSegment {
                start,
                end,
                circulation,
                core2_m2: core + cell.diffusion_offset_m2,
            });
        }
    }
    near.visit_segments_with_core_growth(far.near_core_growth_m2().unwrap(), |s| all.push(s));
    all
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 4e-12 * a.abs().max(b.abs()).max(1.0),
        "{a:?} != {b:?}"
    );
}

fn same_sources(a: &[WakeSegment], b: &[WakeSegment]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.circulation, b.circulation);
        close(a.core2_m2, b.core2_m2);
        for axis in 0..3 {
            close(a.start[axis], b.start[axis]);
            close(a.end[axis], b.end[axis]);
        }
    }
}

fn same_fields(a: (&FilamentWake, &FarField), b: (&FilamentWake, &FarField)) {
    for p in [[-0.5, 0.15, 1.1], [-2.0, -0.2, 0.4], [1.0, 2.0, 3.0]] {
        let x = hybrid_kinematics(a.0, a.1, p).unwrap();
        let y = hybrid_kinematics(b.0, b.1, p).unwrap();
        for axis in 0..3 {
            close(x.velocity[axis], y.velocity[axis]);
            close(x.vorticity()[axis], y.vorticity()[axis]);
            for derivative in 0..3 {
                close(x.gradient[axis][derivative], y.gradient[axis][derivative]);
            }
        }
    }
}

#[test]
fn repeated_transfer_preserves_exact_sources_cores_impulses_and_existing_buffers() {
    let mode = WakeCoreEvolutionMode::CoreSpreading {
        nu_eff_m2ps: 0.0125,
    };
    let mut near = sheet(8, 32);
    let mut far = FarField::aggregate(&mut near, 2, 1, mode, DT).unwrap();
    // Existing cells already have elapsed diffusion and extra material vertices.
    far.shed_prescribed(&mut near, &[0.2, -0.1, 0.05], [-0.03, 0.01, 0.0])
        .unwrap();
    near.rows.last_mut().unwrap().core2_m2 = 0.07;
    far.refine_segments(0.15, LIMITS.max_segments, || false)
        .unwrap();

    for groups in [2, 1] {
        let before_near = near.clone();
        let before_far = far.clone();
        let before_sources = sources(&near, &far);
        let before_impulse = far
            .hydrodynamic_impulse_per_density(&near, [0.3, -0.2, 0.1], 4096, || false)
            .unwrap();
        let pointers: Vec<_> = far
            .cells
            .iter()
            .map(|c| (c.segs.as_ptr(), c.source_core2_m2.as_ptr()))
            .collect();
        let report = far
            .retire_near_rows(&mut near, 2, groups, LIMITS, || false)
            .unwrap();
        assert_eq!(report.rows_transferred, 2);
        assert_eq!(report.near_rows_remaining, before_near.rows.len() - 2);
        assert_eq!(report.cells_added, 2 / groups);
        assert_eq!(&far.cells[..before_far.cells.len()], &before_far.cells);
        for (cell, &(segments, cores)) in far.cells.iter().zip(&pointers) {
            assert_eq!(cell.segs.as_ptr(), segments);
            assert_eq!(cell.source_core2_m2.as_ptr(), cores);
        }
        assert_eq!(sources(&near, &far), before_sources);
        assert_eq!(near.ticks, before_near.ticks);
        assert_eq!(near.rows, before_near.rows[2..]);
        assert_eq!(
            near.downstream_gamma.as_ref().unwrap(),
            &before_near.rows[1].gamma
        );
        same_fields((&near, &far), (&before_near, &before_far));
        assert_eq!(
            far.hydrodynamic_impulse_per_density(&near, [0.3, -0.2, 0.1], 4096, || false)
                .unwrap(),
            before_impulse
        );
        for cell in &far.cells[before_far.cells.len()..] {
            assert_eq!(cell.diffusion_offset_m2, 0.0);
        }
    }

    // Future growth must start from the captured core, not reapply its old age.
    let before = sources(&near, &far);
    far.advance_core(DT).unwrap();
    for row in &mut near.rows {
        row.age_ticks += 1.0;
    }
    for (old, new) in before.iter().zip(sources(&near, &far)) {
        close(new.core2_m2, old.core2_m2 + 4.0 * 0.0125 * DT);
    }
}

#[test]
fn bounded_near_wake_continues_self_induced_motion_beyond_its_row_cap() {
    for mode in [
        WakeCoreEvolutionMode::Frozen,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.01 },
    ] {
        let mut near = sheet(4, 4);
        let mut far = empty(mode);
        let mut full = near.clone();
        full.cert.max_rows = 64;
        let reference_far = empty(mode);
        let config = SelfInducedShed {
            ambient_velocity_mps: [-0.3, 0.01, 0.02],
            emitted_core_radius_m: 0.2,
            core_mode: mode,
            substeps: 4,
            max_displacement_fraction: 1.0,
            max_work: 10_000_000,
        };
        for interval in 0..12 {
            if near.rows.len() == near.cert.max_rows {
                far.retire_near_rows(&mut near, 2, 1, LIMITS, || false)
                    .unwrap();
            }
            let gamma = [0.1 + 0.002 * interval as f64, -0.12, 0.07];
            full.shed_self_induced(&gamma, config, || false).unwrap();
            far.shed_self_induced(&mut near, &gamma, config, || false)
                .unwrap();
            assert!(near.rows.len() <= near.cert.max_rows);
            assert_eq!(near.ticks, full.ticks);
            same_sources(&sources(&near, &far), &sources(&full, &reference_far));
            same_fields((&near, &far), (&full, &reference_far));
        }
        assert_eq!(near.cert.max_rows, 4);
        assert_eq!(full.rows.len(), 16);
        assert!(far.cells.len() >= 12);
    }
}

#[test]
fn total_capacity_and_every_cancellation_checkpoint_are_transactional() {
    let mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.01 };
    let mut original_near = sheet(6, 16);
    let original_far = FarField::aggregate(&mut original_near, 2, 1, mode, DT).unwrap();
    let mut near = original_near.clone();
    let mut far = original_far.clone();
    let mut polls = 0;
    let report = far
        .retire_near_rows(&mut near, 2, 1, LIMITS, || {
            polls += 1;
            false
        })
        .unwrap();
    let exact = FarRetirementLimits {
        max_cells: report.far_cells,
        max_segments: report.far_segments,
    };
    for limits in [
        FarRetirementLimits {
            max_cells: exact.max_cells - 1,
            ..exact
        },
        FarRetirementLimits {
            max_segments: exact.max_segments - 1,
            ..exact
        },
    ] {
        let mut near = original_near.clone();
        let mut far = original_far.clone();
        assert_eq!(
            far.retire_near_rows(&mut near, 2, 1, limits, || false)
                .unwrap_err()
                .code,
            "far-retirement-capacity"
        );
        assert_eq!(near, original_near);
        assert_eq!(far, original_far);
    }
    for stop in 1..=polls {
        let mut near = original_near.clone();
        let mut far = original_far.clone();
        let pointer = far.cells.as_ptr();
        let mut seen = 0;
        assert_eq!(
            far.retire_near_rows(&mut near, 2, 1, exact, || {
                seen += 1;
                seen == stop
            })
            .unwrap_err()
            .code,
            "far-retirement-cancelled"
        );
        assert_eq!(near, original_near);
        assert_eq!(far, original_far);
        assert_eq!(far.cells.as_ptr(), pointer);
        let retried = far
            .retire_near_rows(&mut near, 2, 1, exact, || false)
            .unwrap();
        assert_eq!(retried, report);
    }
}

#[test]
fn empty_groups_identity_and_invalid_clock_or_sources_are_explicit() {
    let mode = WakeCoreEvolutionMode::Frozen;
    let mut near = sheet(4, 8);
    let mut far = empty(mode);
    let before = near.clone();
    let identity = far
        .retire_near_rows(
            &mut near,
            0,
            1,
            FarRetirementLimits {
                max_cells: 0,
                max_segments: 0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(identity.rows_transferred, 0);
    assert_eq!(near, before);
    assert!(far.cells.is_empty());
    for (rows, group) in [(1, 0), (3, 2), (4, 1), (usize::MAX, 1)] {
        assert!(
            far.retire_near_rows(&mut near, rows, group, LIMITS, || false)
                .is_err()
        );
        assert_eq!(near, before);
        assert!(far.cells.is_empty());
    }
    far.shed_dt_s = DT * 2.0;
    assert!(
        far.retire_near_rows(&mut near, 1, 1, LIMITS, || false)
            .is_err()
    );
    assert_eq!(near, before);
    assert!(far.cells.is_empty());
    far.shed_dt_s = DT;
    near.rows[0].age_ticks = f64::NAN;
    assert!(
        far.retire_near_rows(&mut near, 1, 1, LIMITS, || false)
            .is_err()
    );
    assert!(near.rows[0].age_ticks.is_nan());
    assert!(far.cells.is_empty());
    near = before;
    for row in &mut near.rows {
        row.gamma.fill(0.0);
    }
    let zero = far
        .retire_near_rows(
            &mut near,
            2,
            1,
            FarRetirementLimits {
                max_cells: 2,
                max_segments: 0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(zero.far_cells, 2);
    assert_eq!(zero.far_segments, 0);
}
