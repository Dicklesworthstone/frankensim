//! Public-surface seal for the production-run protocol (bead fz2.5).
//!
//! `ProductionProbe`/`ProductionRun` expose no public constructor taking
//! axes, kernels, or a post-probe, so a forged run cannot be built from
//! outside the crate at all — that half of the seal is the type system.
//! What CAN still happen from out here is the old attack: a custom kernel
//! wearing a production name (`gemm-f64`), recorded through the public
//! `record_run` path with caller-supplied axes (including the pre-probe
//! passed twice as the post-probe). This suite proves such evidence is
//! permanently stamped `"protocol":"custom-registry"` and never wears the
//! production stamp, keeping it non-citable no matter how faithfully the
//! kernel imitates the production registry.

use std::sync::atomic::{AtomicU32, Ordering};

use fs_roofline::{
    AxisBaselinePolicy, BaselineAxes, BaselineCandidate, BaselineIdentity, KernelSpec, MachineAxes,
    RooflineKernel, Staleness, TargetAxis, Threading, finalize_registry_tuning, promote_baseline,
    record_run, run_registry, staleness_at,
};

static NEXT_DB: AtomicU32 = AtomicU32::new(0);

fn temp_db(tag: &str) -> String {
    let n = NEXT_DB.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!(
            "fs-roofline-seal-{tag}-{}-{n}.db",
            std::process::id()
        ))
        .display()
        .to_string()
}

fn cleanup_db(path: &str) {
    for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
        let _ = std::fs::remove_file(format!("{path}{suffix}"));
    }
}

fn synthetic_axes(fingerprint: u64) -> MachineAxes {
    // Roofs far above any real machine (bead xjhz).
    MachineAxes {
        fingerprint,
        cpu_brand: "synthetic".to_string(),
        logical_cpus: 8,
        bandwidth_single_gbs: 100_000.0,
        bandwidth_all_core_gbs: 400_000.0,
        peak_single_gflops: 50_000.0,
        peak_all_core_gflops: 300_000.0,
    }
}

fn trusted_baseline(axes: &MachineAxes) -> (BaselineAxes, BaselineIdentity) {
    let identity =
        BaselineIdentity::current(axes, "test-firmware").expect("valid synthetic identity");
    let candidates: Vec<_> = (0_u64..3)
        .map(|ordinal| {
            BaselineCandidate::from_receipt(
                axes.clone(),
                identity.clone(),
                fs_blake3::hash_domain(
                    "fs-roofline.seal-baseline-source.v1",
                    &ordinal.to_le_bytes(),
                ),
            )
            .expect("valid synthetic candidate")
        })
        .collect();
    let baseline = promote_baseline(
        &candidates,
        "test-operator",
        "deterministic seal fixture",
        20_000,
        90,
    )
    .expect("valid synthetic baseline");
    (baseline, identity)
}

/// A caller-controlled kernel that can optionally wear a production identity.
struct ForgedKernel {
    name: &'static str,
    version: &'static str,
    value: u64,
}

impl RooflineKernel for ForgedKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: self.name,
            version: self.version,
            bytes_per_elem: 8.0,
            flops_per_elem: 2.0,
            threading: Threading::AllCore,
            target_axis: TargetAxis::ComputePeak,
            target_fraction: None,
        }
    }

    fn elements(&self) -> usize {
        4096
    }

    fn run_once(&mut self) -> Result<(), String> {
        for _ in 0..1024 {
            self.value = std::hint::black_box(
                self.value
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1),
            );
        }
        Ok(())
    }
}

#[test]
fn forged_name_through_the_public_path_is_stamped_custom_registry() {
    let db = temp_db("forged");
    let ledger = fs_ledger::Ledger::open(&db).expect("open ledger");
    let axes = synthetic_axes(0xF0F0);
    let (baseline, identity) = trusted_baseline(&axes);
    let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);

    // The old attack, end to end: forged production name, caller-supplied
    // axes, and the PRE-probe passed twice (a drifted post-probe silently
    // discarded). The public path still records it as evidence...
    let mut registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ForgedKernel {
        name: "gemm-f64",
        version: "3",
        value: 1,
    })];
    let mut results =
        run_registry(&mut registry, 0, 1, &axes).expect("bounded forged-registry run");
    let version = results[0].version.clone();
    let mut finalized =
        finalize_registry_tuning(&mut registry, &axes, &axes, policy, &results).expect("finalize");
    assert!(
        !finalized.admitted(),
        "a forged GEMM lacks the required execution-decision binding"
    );
    let op = record_run(&ledger, &axes, &axes, policy, &mut finalized, &mut results)
        .expect("public path records");

    // ...but the operation is permanently stamped custom-registry, never
    // production-v2, and carries no run nonce: non-citable by construction.
    let op_row = ledger.op(op).unwrap().expect("op row");
    let recorded_at = op_row.t_end.expect("finished custom operation");
    let ir = op_row.ir;
    assert!(
        ir.contains("\"protocol\":\"custom-registry\""),
        "public-path evidence must be stamped custom-registry: {ir}"
    );
    assert!(!ir.contains("production-v2"));
    assert!(!ir.contains("run_nonce"));
    assert_eq!(
        staleness_at(
            &ledger,
            "gemm-f64",
            &version,
            axes.fingerprint,
            policy.baseline_hash(),
            recorded_at + 1,
        )
        .expect("custom-row staleness"),
        Staleness::NeverMeasured,
        "the independently rejected forged GEMM must not publish a tune row"
    );
    cleanup_db(&db);
}

#[test]
fn admitted_custom_registry_row_cannot_classify_as_fresh() {
    let db = temp_db("admitted-custom");
    let ledger = fs_ledger::Ledger::open(&db).expect("open ledger");
    let axes = synthetic_axes(0xC057);
    let (baseline, identity) = trusted_baseline(&axes);
    let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);

    let mut registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ForgedKernel {
        name: "caller-controlled-f64",
        version: "1",
        value: 1,
    })];
    let mut results =
        run_registry(&mut registry, 0, 1, &axes).expect("bounded cloned-registry run");
    let version = results[0].version.clone();
    let mut finalized =
        finalize_registry_tuning(&mut registry, &axes, &axes, policy, &results).expect("finalize");
    assert!(
        finalized.admitted(),
        "fixture must cross measurement admission to exercise provenance refusal"
    );
    let op = record_run(&ledger, &axes, &axes, policy, &mut finalized, &mut results)
        .expect("public path records admitted custom row");
    let op_row = ledger.op(op).unwrap().expect("op row");
    let recorded_at = op_row.t_end.expect("finished custom operation");
    assert!(op_row.ir.contains("\"protocol\":\"custom-registry\""));
    assert_eq!(
        ledger
            .tune_rows("caller-controlled-f64")
            .expect("custom tune rows")
            .len(),
        1,
        "fixture must publish exactly one row before staleness rejects its provenance"
    );
    assert_eq!(
        staleness_at(
            &ledger,
            "caller-controlled-f64",
            &version,
            axes.fingerprint,
            policy.baseline_hash(),
            recorded_at + 1,
        )
        .expect("custom-row staleness"),
        Staleness::NeverMeasured,
        "custom-registry rows live outside the production staleness namespace"
    );
    cleanup_db(&db);
}
