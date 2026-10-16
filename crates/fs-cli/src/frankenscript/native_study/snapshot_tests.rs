use super::*;
use std::sync::atomic::{AtomicU64, Ordering};

const SOURCE: &str = r#"(fsim-uncertainty-study
 :version 1 :project "cooling-reference.fsim" :samples 4 :seed 29
 :wall-time 120s :method monte-carlo :correlation independent :qoi "temperature-max"
 :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
 :materials ("aa6061.fsmcdpk") :interfaces ()
 :parameters (
  (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
  (uniform :name "ambient" :target convection-temperature :entity "air" :low 294K :high 300K)))"#;

struct Fixture { root: PathBuf, source: PathBuf, pins: StudyPins }
impl Fixture {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let root = loop {
            let path = std::env::temp_dir().join(format!("fs-script-snapshot-{}-{}",
                std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
            match std::fs::create_dir(&path) {
                Ok(()) => break path,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("fixture: {error}"),
            }
        };
        let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
        for name in ["cooling-reference.fsim", "plate.stl", "aa6061.fsmcdpk"] {
            std::fs::copy(reference.join(name), root.join(name)).unwrap();
        }
        let source = root.join("study.fsim");
        std::fs::write(&source, SOURCE).unwrap();
        let base = crate::read_project_for_solve(&root.join("cooling-reference.fsim"), OutputMode::Json)
            .unwrap_or_else(|output| panic!("{}", output.stderr));
        let study = fs_project::uncertainty::UncertaintyStudy::parse(SOURCE).unwrap();
        let pins = StudyPins { project: base.hash(),
            source: Some(fs_blake3::hash_bytes(study.canonical().as_bytes())), wall_seconds: Some(120.0) };
        Self { root, source, pins }
    }
}

#[test]
fn repeated_source_uses_one_snapshot_and_still_enforces_each_pin() {
    let f = Fixture::new();
    let mut cache = SnapshotCache::new(Some(64 * 1024 * 1024));
    let first = cache.prepare(f.source.clone(), f.pins).unwrap();
    let charged = cache.used;
    assert_eq!(charged, first.input_bytes());
    // A subsequent clause refers to the same admitted model, not newly read
    // source or assets. Repeated paths cannot consume duplicate storage.
    std::fs::rename(&f.source, f.root.join("moved-study.fsim")).unwrap();
    std::fs::rename(f.root.join("plate.stl"), f.root.join("moved-plate.stl")).unwrap();
    let second = cache.prepare(f.source.clone(), f.pins).unwrap();
    assert!(Rc::ptr_eq(&first, &second));
    assert_eq!(cache.used, charged);
    assert_eq!(cache.models.len(), 1);
    for pins in [
        StudyPins { source: Some(fs_blake3::hash_bytes(b"changed")), ..f.pins },
        StudyPins { project: fs_blake3::hash_bytes(b"another project"), ..f.pins },
        StudyPins { wall_seconds: Some(60.0), ..f.pins },
    ] {
        assert!(cache.prepare(f.source.clone(), pins).is_err());
        assert_eq!(cache.used, charged);
        assert_eq!(cache.models.len(), 1);
    }
}

#[test]
fn distinct_studies_share_the_aggregate_cap_and_failed_loads_do_not_poison_it() {
    let f = Fixture::new();
    let mut cache = SnapshotCache::new(Some(64 * 1024 * 1024));
    let first = cache.prepare(f.source.clone(), f.pins).unwrap();
    let used = cache.used;
    let second_path = f.root.join("second.fsim");
    std::fs::write(&second_path, SOURCE).unwrap();
    cache.limit = used; // no space for a second model, even with identical bytes
    assert!(cache.prepare(second_path.clone(), f.pins).is_err());
    assert_eq!(cache.used, used);
    assert_eq!(cache.models.len(), 1);
    assert!(Rc::ptr_eq(&first, &cache.prepare(f.source.clone(), f.pins).unwrap()));
    cache.limit = 16 * 1024 * 1024;
    let wrong = StudyPins { source: Some(fs_blake3::hash_bytes(b"bad pin")), ..f.pins };
    assert!(cache.prepare(second_path.clone(), wrong).is_err());
    assert_eq!(cache.used, used);
    assert!(!cache.models.contains_key(&second_path));
    let second = cache.prepare(second_path, f.pins).unwrap();
    assert!(!Rc::ptr_eq(&first, &second));
    assert_eq!(cache.used, used + second.input_bytes());
    assert_eq!(cache.models.len(), 2);
}

#[test]
fn no_memory_grant_does_not_become_an_unlimited_snapshot_budget() {
    let f = Fixture::new();
    for memory in [None, Some(0), Some(3)] {
        let mut cache = SnapshotCache::new(memory);
        let error = cache.prepare(f.source.clone(), f.pins).unwrap_err();
        assert!(error.contains("input-storage allowance"), "{error}");
        assert_eq!(cache.used, 0);
        assert!(cache.models.is_empty());
    }
}
