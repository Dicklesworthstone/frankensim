use super::*;

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-adaptive.fsim"));
const DECLARATION: &str = "(design-regions :solid (((0.0 0.0 0.0) (0.5 0.5 0.5))) :void (((0.5 0.5 0.5) (1.0 1.0 1.0))))";
fn source(declaration: &str) -> String {
    let mut text = BASE.replace(":updates-per-stage 3", ":updates-per-stage 1")
        .replace(":schedule ((1.0 1.0) (2.0 2.0) (3.0 8.0))", ":schedule ((1.0 1.0) (2.0 2.0))");
    let last = text.rfind(')').unwrap();
    text.insert_str(last, &format!("\n  {declaration}\n"));
    text
}

#[test]
fn declared_si_regions_are_explicit_bounded_and_identity_bearing() {
    let spec = spec::parse(&source(DECLARATION)).unwrap();
    assert_eq!(spec.regions.len(), 2);
    assert_eq!(spec.regions[0], RegionBox { lower: [0; 3], upper: [1; 3], kind: PhysicalRegion3::Solid });
    assert_eq!(spec.regions[1], RegionBox { lower: [1; 3], upper: [2; 3], kind: PhysicalRegion3::Void });
    assert_eq!(spec::parse(&spec.canonical).unwrap().id, spec.id);
    assert_eq!(spec::parse(&spec.canonical).unwrap().regions, spec.regions);
    assert!(spec::parse(BASE).unwrap().regions.is_empty());
    let other = DECLARATION.replace("(0.5 0.5 0.5))) :void", "(1.0 0.5 0.5))) :void");
    assert_ne!(spec::parse(&source(&other)).unwrap().id, spec.id);
    for invalid in [
        DECLARATION.replace("design-regions", "ignored-regions"),
        DECLARATION.replace(":solid", ":unknown"),
        DECLARATION.replace("(0.5 0.5 0.5))) :void", "(0.25 0.5 0.5))) :void"),
        DECLARATION.replace("(0.0 0.0 0.0)", "(-0.5 0.0 0.0)"),
        DECLARATION.replace("(1.0 1.0 1.0)", "(0.5 1.0 1.0)"),
        "(design-regions :solid (((0 0 0) (1 1 1))) :void (((0 0 0) (0.5 0.5 0.5))))".into(),
        "(design-regions :solid () :void () :extra ())".into(),
    ] {
        assert!(spec::parse(&source(&invalid)).is_err(), "must refuse invalid/partial region: {invalid}");
    }
    let boxes = "((0 0 0) (0.5 0.5 0.5)) ".repeat(33);
    assert!(spec::parse(&source(&format!("(design-regions :solid ({boxes}) :void ())"))).is_err());
}

#[test]
fn a_translated_physical_box_maps_planes_to_the_same_initial_cell_indices() {
    let mut spec = spec::parse(BASE).unwrap();
    spec.bounds = ([-1.0, 2.0, 4.0], [1.0, 4.0, 6.0]);
    let node = fs_ir::sexpr::parse("(design-regions :solid (((-1 2 4) (0 3 5))) :void ())").unwrap();
    assert_eq!(parse(&node, &spec).unwrap(),
        [RegionBox { lower: [0; 3], upper: [1; 3], kind: PhysicalRegion3::Solid }]);
}

fn assert_retained_regions(spec: &Spec, result: &Computation) {
    assert_eq!(result.status, "completed", "{:?}", result.report);
    assert_eq!(result.completed_stages, 2);
    assert!(result.report.refinements.iter().any(|r| r.installed));
    let labels = result.study.physical_regions().expect("physical labels survive adaptive replacement");
    let leaves = result.study.operator().elasticity().leaves();
    let accepted = result.report.continuation.last.as_ref().unwrap();
    assert_eq!(labels.len(), leaves.len());
    assert_eq!(accepted.projected_rho.len(), leaves.len());
    let mut solid = 0; let mut void = 0; let mut design = 0;
    for ((leaf, label), density) in leaves.iter().zip(labels).zip(&accepted.projected_rho) {
        let ancestor = leaf.index().map(|i| i >> (u32::from(leaf.level()) - spec.level));
        let expected = spec.regions.iter().find(|region| (0..3).all(|a|
            region.lower[a] <= ancestor[a] && ancestor[a] < region.upper[a]))
            .map_or(PhysicalRegion3::Design, |region| region.kind);
        assert_eq!(*label, expected);
        match label {
            PhysicalRegion3::Solid => { solid += 1; assert_eq!(*density, 1.0); }
            PhysicalRegion3::Void => { void += 1; assert_eq!(*density, 0.0); }
            PhysicalRegion3::Design => { design += 1; assert!((0.0..=1.0).contains(density)); }
        }
    }
    assert!(solid > 0 && void > 0 && design > 0);
    assert_eq!(accepted.displacements.len(), 2);
    assert_ne!(accepted.displacements[0], accepted.displacements[1]);
    for stage in &result.report.continuation.stages {
        assert!(stage.history.len() > 1, "each stage must make an actual design update");
        assert!(stage.history.iter().all(|row| row.volume_fraction <= spec.volume + 1e-8));
        assert!(stage.history.windows(2).all(|rows| rows[1].compliance <= rows[0].compliance));
    }
}

#[test]
fn actual_adaptive_study_preserves_material_regions_and_replays_accepted_fields() {
    let spec = spec::parse(&source(DECLARATION)).unwrap();
    let first = compute(&spec, &CancelGate::new()).unwrap();
    assert_retained_regions(&spec, &first);
    let second = compute(&spec, &CancelGate::new()).unwrap();
    assert_retained_regions(&spec, &second);
    let a = first.report.continuation.last.as_ref().unwrap();
    let b = second.report.continuation.last.as_ref().unwrap();
    assert_eq!(a.rho, b.rho);
    assert_eq!(a.projected_rho, b.projected_rho);
    assert_eq!(a.displacements, b.displacements);
    assert_eq!(first.study.physical_regions(), second.study.physical_regions());
}

#[test]
fn pre_cancelled_region_study_never_publishes_a_baseline() {
    let spec = spec::parse(&source(DECLARATION)).unwrap();
    let gate = CancelGate::new(); gate.request();
    let error = match compute(&spec, &gate) {
        Err(error) => error,
        Ok(_) => panic!("cancelled region study cannot be published"),
    };
    assert_eq!(error.code, "cli-study-sdf3-cancelled");
}
