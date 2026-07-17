//! D3Q19 performance-model and tropical-attribution conformance (bead 712t).

use fs_lbm::d3q19::{E3, TILE};
use fs_lbm::perf::{
    BGK_EQUILIBRIUM_FLOPS, BGK_FLOPS_PER_CELL, BGK_FORCE_RELAX_FLOPS, BGK_MACRO_VELOCITY_FLOPS,
    D3Q19_DISTRIBUTIONS, D3Q19_HALO_LINKS_PER_TILE, D3Q19_LINKS_PER_TILE,
    D3Q19_LOCAL_LINKS_PER_TILE, D3Q19_PERF_MODEL_VERSION, D3q19PerfRow, D3q19TrafficModel,
    DISTRIBUTION_BYTES, EvidenceClass, KernelClass, LaneShape, OccupancyClass, PerfGateVerdict,
    PerfModelError, RATIO_PPM, ReferenceIsa, SPARSE_TILE_CELLS, SPARSE_TILE_EDGE, TaskSample,
    ThreadingClass, attribute_critical_path, critical_path_is_stable,
};

fn task(id: u32, class: KernelClass, wall_ns: u64, predecessors: &[u32]) -> TaskSample {
    TaskSample {
        id,
        class,
        wall_ns,
        predecessors: predecessors.to_vec(),
    }
}

fn representative_tasks() -> Vec<TaskSample> {
    vec![
        task(10, KernelClass::Activation, 10, &[]),
        task(20, KernelClass::Collide, 50, &[10]),
        task(30, KernelClass::Halo, 20, &[20]),
        task(40, KernelClass::Stream, 40, &[20]),
        task(50, KernelClass::Stream, 15, &[30, 40]),
    ]
}

fn critical() -> fs_lbm::perf::CriticalPathAttribution {
    attribute_critical_path(&representative_tasks()).expect("representative DAG admits")
}

fn report_only_row() -> D3q19PerfRow {
    D3q19PerfRow {
        reference_isa: ReferenceIsa::AppleMClass,
        shape: LaneShape::memory_resident(
            OccupancyClass::SparseTenPercent,
            ThreadingClass::AllCore,
            12,
        )
        .expect("shape admits"),
        glups: 0.72,
        dispersion_ppm: 41_000,
        floor_glups: None,
        evidence: EvidenceClass::ReportOnly {
            reason: "first quiet-window calibration has not authorized a floor".to_owned(),
        },
        placement_identity: "tile-edge=4;workers=12;placement=host".to_owned(),
        critical_paths: vec![critical(); 5],
    }
}

#[test]
fn memory_resident_shapes_round_sparse_occupancy_up_to_whole_tiles() {
    assert_eq!(SPARSE_TILE_EDGE, 4);
    assert_eq!(SPARSE_TILE_CELLS, 64);
    let dense =
        LaneShape::memory_resident(OccupancyClass::DenseActive, ThreadingClass::SingleThread, 1)
            .expect("dense shape");
    assert_eq!(dense.total_tiles().unwrap(), 32 * 32 * 32);
    assert_eq!(dense.active_tiles().unwrap(), 32 * 32 * 32);
    assert_eq!(dense.active_cells().unwrap(), 128 * 128 * 128);
    assert_eq!(dense.allocated_population_bytes().unwrap(), 956_301_312);

    let sparse = LaneShape::memory_resident(
        OccupancyClass::SparseTenPercent,
        ThreadingClass::AllCore,
        32,
    )
    .expect("sparse shape");
    assert_eq!(sparse.active_tiles().unwrap(), 3_277);
    assert_eq!(sparse.active_cells().unwrap(), 209_728);
    assert!(
        sparse.active_cells().unwrap() * RATIO_PPM as usize
            >= dense.active_cells().unwrap() * 100_000
    );
    assert!(
        (sparse.active_cells().unwrap() - SPARSE_TILE_CELLS) * RATIO_PPM as usize
            < dense.active_cells().unwrap() * 100_000,
        "one fewer whole tile would undershoot ten percent"
    );
}

#[test]
fn shape_validation_refuses_zero_misalignment_and_fake_single_thread_rows() {
    let zero = LaneShape {
        dims: [128, 0, 128],
        occupancy: OccupancyClass::DenseActive,
        threading: ThreadingClass::SingleThread,
        workers: 1,
    };
    assert!(matches!(
        zero.validate(),
        Err(PerfModelError::InvalidDimensions(_))
    ));
    let misaligned = LaneShape {
        dims: [127, 128, 128],
        ..zero
    };
    assert!(matches!(
        misaligned.validate(),
        Err(PerfModelError::InvalidDimensions(_))
    ));
    let fake_serial = LaneShape {
        dims: [128; 3],
        workers: 2,
        ..zero
    };
    assert!(matches!(
        fake_serial.validate(),
        Err(PerfModelError::InvalidWorkers { .. })
    ));
}

#[test]
fn arithmetic_intensity_header_counts_population_and_sparse_metadata() {
    let model = D3q19TrafficModel::default();
    assert_eq!(D3Q19_DISTRIBUTIONS, 19);
    assert_eq!(DISTRIBUTION_BYTES, 8);
    assert_eq!(BGK_MACRO_VELOCITY_FLOPS, 142);
    assert_eq!(BGK_EQUILIBRIUM_FLOPS, 271);
    assert_eq!(BGK_FORCE_RELAX_FLOPS, 611);
    assert_eq!(BGK_FLOPS_PER_CELL, 1_024);
    assert_eq!(D3Q19_LINKS_PER_TILE, 1_216);
    assert_eq!(D3Q19_LOCAL_LINKS_PER_TILE, 784);
    assert_eq!(D3Q19_HALO_LINKS_PER_TILE, 432);
    assert_eq!(
        D3Q19_LOCAL_LINKS_PER_TILE + D3Q19_HALO_LINKS_PER_TILE,
        D3Q19_LINKS_PER_TILE
    );
    let mut independently_counted_halo_links = 0usize;
    for lz in 0..TILE {
        for ly in 0..TILE {
            for lx in 0..TILE {
                for &(ex, ey, ez) in &E3 {
                    let sx = lx as i64 - i64::from(ex);
                    let sy = ly as i64 - i64::from(ey);
                    let sz = lz as i64 - i64::from(ez);
                    let tile = TILE as i64;
                    if !(0..tile).contains(&sx)
                        || !(0..tile).contains(&sy)
                        || !(0..tile).contains(&sz)
                    {
                        independently_counted_halo_links += 1;
                    }
                }
            }
        }
    }
    assert_eq!(independently_counted_halo_links, D3Q19_HALO_LINKS_PER_TILE);
    assert_eq!(
        model.population_bytes_per_cell().to_bits(),
        608.0f64.to_bits()
    );
    assert_eq!(
        model.sparse_overhead_bytes_per_cell().to_bits(),
        108.25f64.to_bits()
    );
    assert_eq!(model.bytes_per_cell().to_bits(), 716.25f64.to_bits());
    assert_eq!(
        model.arithmetic_intensity().to_bits(),
        (1_024.0f64 / 716.25).to_bits()
    );
    let receipt = model.receipt_json();
    assert!(receipt.contains(D3Q19_PERF_MODEL_VERSION));
    assert!(receipt.contains("\"links_per_tile\":1216"));
    assert!(receipt.contains("\"local_links_per_tile\":784"));
    assert!(receipt.contains("\"halo_links_per_tile\":432"));
    assert!(receipt.contains("\"bytes_per_cell\":716.250000"));
    assert!(receipt.contains("\"sparse_overhead_bytes_per_cell\":108.250000"));
}

#[test]
fn plan_targets_are_reported_by_reference_family_not_used_as_floor() {
    assert_eq!(ReferenceIsa::AppleMClass.plan_target_glups(), Some(1.0));
    assert_eq!(
        ReferenceIsa::ThreadripperClass.plan_target_glups(),
        Some(0.6)
    );
    assert_eq!(ReferenceIsa::Other.plan_target_glups(), None);
    let row = report_only_row();
    assert_eq!(row.plan_target_met(), Some(false));
    assert_eq!(row.floor_glups, None);
    assert_eq!(row.gate_verdict().unwrap(), PerfGateVerdict::ReportOnly);
}

#[test]
fn max_plus_path_names_the_true_dominant_kernel_class() {
    let result = critical();
    assert_eq!(result.path, [10, 20, 40, 50]);
    assert_eq!(result.makespan_ns, 115);
    assert_eq!(result.class_wall_ns, [10, 50, 0, 55]);
    assert_eq!(result.dominant_class, KernelClass::Stream);
    assert_eq!(result.dominant_share_ppm, 478_261);
    result.validate().expect("derived attribution validates");
    let json = result.receipt_json();
    assert!(json.contains("\"dominant_class\":\"stream\""));
    assert!(json.contains("\"path\":[10,20,40,50]"));
}

#[test]
fn max_plus_attribution_is_invariant_to_task_and_predecessor_enumeration() {
    let expected = critical();
    let mut permuted = representative_tasks();
    permuted.reverse();
    permuted
        .iter_mut()
        .find(|sample| sample.id == 50)
        .expect("join task")
        .predecessors
        .reverse();
    assert_eq!(
        attribute_critical_path(&permuted).expect("permutation admits"),
        expected
    );
}

#[test]
fn exact_path_ties_break_to_the_smaller_task_identity() {
    let tied = [
        task(1, KernelClass::Activation, 10, &[]),
        task(2, KernelClass::Collide, 20, &[1]),
        task(3, KernelClass::Halo, 20, &[1]),
        task(4, KernelClass::Stream, 5, &[3, 2]),
    ];
    let result = attribute_critical_path(&tied).expect("tied DAG admits");
    assert_eq!(result.path, [1, 2, 4]);
    assert_eq!(result.dominant_class, KernelClass::Collide);
}

#[test]
fn malformed_timing_graphs_fail_closed() {
    let duplicate_task = [
        task(1, KernelClass::Collide, 1, &[]),
        task(1, KernelClass::Stream, 1, &[]),
    ];
    assert_eq!(
        attribute_critical_path(&duplicate_task),
        Err(PerfModelError::DuplicateTask(1))
    );

    let duplicate_edge = [
        task(1, KernelClass::Collide, 1, &[]),
        task(2, KernelClass::Stream, 1, &[1, 1]),
    ];
    assert!(matches!(
        attribute_critical_path(&duplicate_edge),
        Err(PerfModelError::DuplicatePredecessor { .. })
    ));

    let missing = [task(2, KernelClass::Stream, 1, &[99])];
    assert!(matches!(
        attribute_critical_path(&missing),
        Err(PerfModelError::MissingPredecessor { .. })
    ));

    let cycle = [
        task(1, KernelClass::Collide, 1, &[2]),
        task(2, KernelClass::Stream, 1, &[1]),
    ];
    assert_eq!(
        attribute_critical_path(&cycle),
        Err(PerfModelError::CyclicTaskGraph)
    );

    let zero = [task(1, KernelClass::Collide, 0, &[])];
    assert_eq!(
        attribute_critical_path(&zero),
        Err(PerfModelError::ZeroWallSample(1))
    );
}

#[test]
fn deterministic_repetitions_require_the_same_path_and_dominant_class() {
    let first = critical();
    let mut faster = representative_tasks();
    for sample in &mut faster {
        sample.wall_ns *= 2;
    }
    let second = attribute_critical_path(&faster).expect("scaled repetition");
    assert!(critical_path_is_stable(&[first.clone(), second]));

    let mut changed = representative_tasks();
    changed
        .iter_mut()
        .find(|sample| sample.id == 30)
        .unwrap()
        .wall_ns = 80;
    let changed = attribute_critical_path(&changed).expect("changed repetition");
    assert!(!critical_path_is_stable(&[first, changed]));
    assert!(!critical_path_is_stable(&[]));
}

#[test]
fn citable_rows_require_receipt_floor_repetitions_and_stable_attribution() {
    let mut row = report_only_row();
    row.evidence = EvidenceClass::Citable {
        admission_receipt: "ab".repeat(32),
    };
    assert!(row.validate().is_err(), "citable row needs a floor");
    row.floor_glups = Some(0.5);
    row.critical_paths.truncate(2);
    assert!(row.validate().is_err(), "citable row needs at least 3 reps");
    row.critical_paths = vec![critical(); 5];
    let mut unstable_tasks = representative_tasks();
    unstable_tasks
        .iter_mut()
        .find(|sample| sample.id == 30)
        .unwrap()
        .wall_ns = 80;
    row.critical_paths[4] =
        attribute_critical_path(&unstable_tasks).expect("alternate attribution");
    assert!(row.validate().is_err(), "citable row needs stable path");
    row.critical_paths = vec![critical(); 5];
    row.validate().expect("fully admitted row");
    assert_eq!(row.gate_verdict().unwrap(), PerfGateVerdict::FloorMet);
    row.glups = 0.49;
    assert_eq!(row.gate_verdict().unwrap(), PerfGateVerdict::FloorMiss);

    row.evidence = EvidenceClass::Citable {
        admission_receipt: "not-a-hash".to_owned(),
    };
    assert!(
        row.validate().is_err(),
        "citable receipt is exact lowercase hex"
    );
}

#[test]
fn environment_invalid_is_neither_floor_pass_nor_floor_miss() {
    let mut row = report_only_row();
    row.floor_glups = Some(0.5);
    row.evidence = EvidenceClass::EnvironmentInvalid {
        reason: "post-run axes drift exceeded the admitted band".to_owned(),
    };
    assert_eq!(
        row.gate_verdict().unwrap(),
        PerfGateVerdict::EnvironmentInvalid
    );
}

#[test]
fn receipt_json_retains_model_shape_target_floor_admission_and_critical_path() {
    let row = report_only_row();
    let json = row.receipt_json().expect("report-only receipt");
    for fragment in [
        "\"metric\":\"lbm-d3q19-sweep\"",
        "\"model_version\":\"d3q19-sparse-sweep-v2-local-halo\"",
        "\"dims\":[128,128,128]",
        "\"occupancy\":\"sparse-ten-percent\"",
        "\"active_tiles\":3277",
        "\"threading\":\"all-core\"",
        "\"plan_target_glups\":1.000000",
        "\"plan_target_met\":false",
        "\"floor_glups\":null",
        "\"evidence_class\":\"report_only\"",
        "\"gate_verdict\":\"report_only\"",
        "\"critical_path_stable\":true",
        "\"dominant_class\":\"stream\"",
    ] {
        assert!(json.contains(fragment), "missing {fragment} from {json}");
    }
}

#[test]
fn forged_favorable_critical_path_totals_refuse() {
    let mut forged = critical();
    forged.class_wall_ns[KernelClass::Collide as usize] += 1;
    assert!(forged.validate().is_err());

    let mut row = report_only_row();
    row.critical_paths[0] = forged;
    assert!(row.validate().is_err());
}
