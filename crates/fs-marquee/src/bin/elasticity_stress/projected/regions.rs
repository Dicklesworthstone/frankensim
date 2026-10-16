//! Reuse the authored-region parser, rasterizer and exporter of the multi-load CLI.
use super::{field as write_field, writer};
use fs_topols::GridSdf;
use fs_topols::design_regions::{DesignRegionStage, PreparedDesignRegions};
use std::error::Error;
use std::ops::ControlFlow;
use std::path::Path;

// The same source serves both binaries; no parallel CSV/rasterization protocol.
#[path = "../../elasticity_robust/projected/design_regions.rs"]
mod authoring;

fn numbers(values: &[f64]) -> String {
    format!("[{}]", values.iter().map(|value| format!("{value:.17e}")).collect::<Vec<_>>().join(","))
}

pub(super) struct Regions(authoring::Authoring);
impl Regions {
    pub(super) fn prepared(&self) -> &PreparedDesignRegions { &self.0.prepared }
    pub(super) fn export(&self, output: &Path) -> Result<String, Box<dyn Error>> { self.0.export(output) }
}

pub(super) fn options(args: &[String]) -> Result<(Vec<String>, Option<String>), Box<dyn Error>> {
    authoring::options(args)
}

pub(super) fn load_controlled<B>(
    path: &Path, field: &GridSdf, fixed: &[(usize, f64)],
    control: impl FnMut(DesignRegionStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, Regions>, Box<dyn Error>> {
    match authoring::load_controlled(path, field, fixed, control)? {
        ControlFlow::Continue(authored) => Ok(ControlFlow::Continue(Regions(authored))),
        ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{checkpoint, run_at};
    use fs_marquee::level_set_csv::read_field;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::time::Instant;

    static NEXT: AtomicUsize = AtomicUsize::new(0);
    struct Scratch(PathBuf);
    impl Scratch {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!("fs-stress-regions-{}-{}",
                std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }
    impl Drop for Scratch {
        fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); }
    }
    const REGIONS: &str = "material,0.25,0.375,0.375,0.5,0.02\nvoid,0.625,0.625,0.75,0.75,0.02\n";
    fn args(output: &Path, regions: &Path, pause: bool) -> Vec<String> {
        let mut args = vec![output.to_string_lossy().into_owned()];
        args.extend(["1e12", "3", "2", "0.6", "8", "0", "300", "--checkpoint", "--design-regions"].map(str::to_string));
        args.push(regions.to_string_lossy().into_owned());
        if pause { args.extend(["--pause-after", "0"].map(str::to_string)); }
        args
    }
    fn fixed(field: &GridSdf) -> Vec<(usize, f64)> {
        field.nodes().iter().copied().enumerate().filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect()
    }
    fn matches_policy(field: &GridSdf, policy: &checkpoint::Policy) {
        for &(index, value) in &policy.fixed {
            assert_eq!(field.nodes()[index].to_bits(), value.to_bits(), "fixed node {index}");
        }
        for j in 3..=4 { for i in 2..=3 { assert!(field.node(i, j) <= -0.02); } }
        for j in 5..=6 { for i in 5..=6 { assert!(field.node(i, j) >= 0.02); } }
    }

    #[test]
    fn regions_are_authored_and_preserved_through_real_stress_search_and_disk_resume() {
        let scratch = Scratch::new();
        let regions = scratch.0.join("regions.csv");
        std::fs::write(&regions, REGIONS).unwrap();
        let start = scratch.0.join("start");
        assert_eq!(run_at(&args(&start, &regions, true), Instant::now()).unwrap(), 6);
        let checkpoint_path = start.join("checkpoint-0000.fscp");
        let bytes = std::fs::read(&checkpoint_path).unwrap();
        assert!(bytes[65..].starts_with(b"fs-marquee-projected-stress-checkpoint-v2\n"));
        let executable = checkpoint::executable().unwrap();
        let ControlFlow::Continue((_, policy)) = checkpoint::load_controlled(&checkpoint_path,
            &executable, |_| ControlFlow::<()>::Continue(())).unwrap()
        else { panic!("uninterrupted checkpoint recovery stopped") };
        assert!(policy.fixed.len() > 18);
        let baseline = read_field(&start.join("baseline-level-set.csv"), 8).unwrap();
        matches_policy(&baseline, &policy);
        let source = read_field(&start.join("input-level-set.csv"), 8).unwrap();
        assert!(source.node(5, 5) < 0.0 && baseline.node(5, 5) > 0.0,
            "the void must actually remove authored input material");
        let summary = std::fs::read_to_string(start.join("summary.json")).unwrap();
        assert!(summary.contains("whole_intersected_cells"));
        assert!(start.join("design-regions.csv").is_file());
        assert!(start.join("design-region-level-set.csv").is_file());

        let full = scratch.0.join("full");
        let full_exit = run_at(&args(&full, &regions, false), Instant::now()).unwrap();
        assert!(matches!(full_exit, 0 | 11));
        assert!(!std::fs::read(full.join("attempts.jsonl")).unwrap().is_empty(), "must attempt real candidates");
        // Resume is self-contained and must not re-author from an external file.
        std::fs::rename(&regions, regions.with_extension("retained-input")).unwrap();
        let resumed = scratch.0.join("resumed");
        let resume_args = vec!["--resume".into(), checkpoint_path.to_string_lossy().into_owned(),
            resumed.to_string_lossy().into_owned(), "--wall-seconds".into(), "300".into()];
        assert_eq!(run_at(&resume_args, Instant::now()).unwrap(), full_exit);
        for name in ["level-set.csv", "trajectory.jsonl", "attempts.jsonl", "stress-checks.jsonl"] {
            assert_eq!(std::fs::read(full.join(name)).unwrap(), std::fs::read(resumed.join(name)).unwrap(), "{name}");
        }
        for output in [&full, &resumed] {
            matches_policy(&read_field(&output.join("level-set.csv"), 8).unwrap(), &policy);
            for item in std::fs::read_dir(output).unwrap() {
                let path = item.unwrap().path();
                if path.file_name().unwrap().to_string_lossy().starts_with("accepted-") {
                    matches_policy(&read_field(&path, 8).unwrap(), &policy);
                }
            }
        }
        let summary = std::fs::read_to_string(resumed.join("summary.json")).unwrap();
        assert!(summary.contains(&format!("\"fixed_node_count\":{}", policy.fixed.len())));
        assert!(summary.contains("\"design_regions\":null"), "do not fabricate historical source records");
        assert_eq!(std::fs::read(checkpoint_path).unwrap(), bytes);
    }

    #[test]
    fn conflicting_regions_cannot_modify_the_load_or_publish_a_study() {
        let scratch = Scratch::new();
        for (i, text) in [
            "void,0.875,0.375,1,0.625,0.02\n",
            "material,0.25,0.375,0.375,0.5,0.02\nvoid,0.25,0.375,0.375,0.5,0.02\n",
        ].iter().enumerate() {
            let regions = scratch.0.join(format!("conflict-{i}.csv"));
            std::fs::write(&regions, text).unwrap();
            let output = scratch.0.join(format!("refused-{i}"));
            let error = run_at(&args(&output, &regions, true), Instant::now()).unwrap_err();
            assert!(error.to_string().contains("conflict"), "{error}");
            assert!(!output.exists());
            assert_eq!(std::fs::read_to_string(regions).unwrap(), *text);
        }
    }

    #[test]
    fn authoring_can_cancel_between_regions_without_changing_the_source() {
        let scratch = Scratch::new();
        let path = scratch.0.join("regions.csv");
        std::fs::write(&path, REGIONS).unwrap();
        let source = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
        let before: Vec<_> = source.nodes().iter().map(|v| v.to_bits()).collect();
        let result = load_controlled(&path, &source, &fixed(&source), |stage| {
            if matches!(stage, DesignRegionStage::Rasterize { region: 1, .. }) {
                ControlFlow::Break("stop authoring")
            } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break("stop authoring")));
        assert_eq!(source.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>(), before);
        let sync = authoring::load(&path, &source, &fixed(&source)).unwrap();
        let ControlFlow::Continue(retry) = load_controlled(&path, &source, &fixed(&source),
            |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("retry stopped") };
        assert_eq!(sync.prepared.geometry.nodes(), retry.prepared().geometry.nodes());
        assert_eq!(sync.prepared.fixed_nodes, retry.prepared().fixed_nodes);
    }

    #[test]
    fn expanded_checkpoints_refuse_downgrade_and_changed_prescribed_nodes_before_physics() {
        let scratch = Scratch::new();
        let regions = scratch.0.join("regions.csv");
        std::fs::write(&regions, REGIONS).unwrap();
        let output = scratch.0.join("study");
        assert_eq!(run_at(&args(&output, &regions, true), Instant::now()).unwrap(), 6);
        let original = std::fs::read(output.join("checkpoint-0000.fscp")).unwrap();
        let executable = checkpoint::executable().unwrap();
        let magic_length = b"fs-marquee-projected-stress-checkpoint-v2\n".len();
        let mut downgraded = original.clone();
        downgraded[65 + magic_length - 2] = b'1';
        let mut changed = original.clone();
        let fixed_start = 65 + magic_length + 64 + 35 * 8;
        let count_start = 65 + magic_length + 64 + 7 * 8;
        let count = u64::from_le_bytes(original[count_start..count_start + 8].try_into().unwrap()) as usize;
        let record = (0..count).find(|index| {
            let start = fixed_start + index * 16;
            let node = u64::from_le_bytes(original[start..start + 8].try_into().unwrap());
            node % 9 != 0 && node % 9 != 8
        }).unwrap();
        let start = fixed_start + record * 16 + 8;
        changed[start..start + 8].copy_from_slice(&0.0f64.to_bits().to_le_bytes());
        for (i, mut bytes) in [downgraded, changed].into_iter().enumerate() {
            let hash = fs_ledger::hash_bytes(&bytes[65..]).to_string();
            bytes[..64].copy_from_slice(hash.as_bytes());
            let path = scratch.0.join(format!("invalid-{i}.fscp"));
            std::fs::write(&path, bytes).unwrap();
            assert!(checkpoint::load_controlled(&path, &executable,
                |_| -> ControlFlow<()> { panic!("invalid fixed masks must refuse before any PDE work") }).is_err());
        }
    }
}
