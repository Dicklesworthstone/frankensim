use super::*;
use fs_topopt::pipeline::LoadCase;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-adaptive.fsim"));
fn physical() -> Spec {
    spec::parse(&FIXTURE.replace("curved-height-sdf", "physical-curved-height-sdf")).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * b.abs().max(1e-20), "{a:e} != {b:e}");
}

#[test]
fn physical_dimensions_and_support_are_explicit_and_legacy_envelope_is_unchanged() {
    let old = spec::parse(FIXTURE).unwrap();
    assert!(!old.physical);
    assert!(spec::parse(&FIXTURE.replace("(1.0 1.0 1.0)", "(0.1 0.05 0.04)")).is_err());
    let source = FIXTURE.replace("curved-height-sdf", "physical-curved-height-sdf")
        .replace("(1.0 1.0 1.0)", "(0.1 0.05 0.04)")
        .replace(":height-m 0.7", ":height-m 0.028")
        .replace(":curvature-per-m 0.1", ":curvature-per-m 1.0")
        .replace(":filter-radius-m 0.15", ":filter-radius-m 0.01")
        .replace(":fixed-boundary left", ":fixed-boundary right");
    let admitted = spec::parse(&source).unwrap();
    assert_eq!(admitted.bounds, ([0.0; 3], [0.1, 0.05, 0.04]));
    assert_eq!(admitted.fixed, FixedFace::Right);
    assert_ne!(admitted.id, old.id);
    assert_eq!(spec::parse(&admitted.canonical).unwrap().id, admitted.id);
    for bad in [
        source.replace("(0.1 0.05 0.04)", "(-0.1 0.05 0.04)"),
        source.replace("(0.1 0.05 0.04)", "(0.1 0.05 0.000001)"),
        source.replace(":height-m 0.028", ":height-m 0.7"),
        source.replace(":curvature-per-m 1.0", ":curvature-per-m 100.0"),
        source.replace(":filter-radius-m 0.01", ":filter-radius-m 0.0"),
        source.replace(":fixed-boundary right", ":fixed-boundary top"),
    ] {
        assert!(spec::parse(&bad).is_err(), "unsupported geometry/support must refuse");
    }
}

#[test]
fn physical_field_enclosures_and_tangents_follow_translated_si_coordinates() {
    let domain = PhysicalDomain {
        bounds: ([0.25, -0.5, 1.0], [0.375, -0.375, 1.125]),
        height: 0.0875, curvature: 0.8,
    };
    for i in 0..32 {
        let x = 0.25 + (i as f64 + 0.5) / 256.0;
        let lo = [x - 1.0 / 1024.0, -0.45, 1.07];
        let hi = [x + 1.0 / 1024.0, -0.40, 1.10];
        let enclosed = domain.enclose(lo, hi);
        let derivative = domain.derivative_enclose(lo, hi, HeightAxis::X);
        for px in [lo[0], x, hi[0]] {
            for pz in [lo[2], hi[2]] {
                let value = domain.value([px, -0.42, pz]);
                assert!(enclosed.lo() <= value && value <= enclosed.hi());
            }
            let dx = -domain.curvature * ((domain.bounds.1[0] - px) - (px - domain.bounds.0[0]));
            assert!(derivative.lo() <= dx && dx <= derivative.hi());
        }
    }
}

fn response(spec: &Spec) -> (f64, Vec<f64>, Vec<[f64; 3]>, f64, Vec<bool>) {
    validate(spec).unwrap();
    let tree = Octree3::uniform(1, 4, 1000).unwrap();
    let domain = PhysicalDomain { bounds: spec.bounds, height: spec.height, curvature: spec.curvature };
    let mut poll = |_| ControlFlow::Continue(());
    let mut quadrature = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let op = AdaptiveElasticity3::build(
        HexCell::try_new(spec.bounds.0, spec.bounds.1).unwrap(), &tree, &domain,
        &IsotropicElastic::new(spec.youngs, spec.poisson, 1.0).unwrap(),
        &|p| spec.fixed.contains(p, spec.bounds), ElasticityOptions3::default(), &mut quadrature,
    ).unwrap();
    let nodes = op.nodes().to_vec();
    let fixed = op.fixed().to_vec();
    let volume: f64 = op.volumes().iter().sum();
    let force = op.body_load(&|_| spec.loads[0].0, || ControlFlow::Continue(())).unwrap();
    let mut study = CutDensityStudy3::new(op, spec.radius, spec.schedule[0]);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = study.evaluate(&vec![0.5; study.cells()],
        &[LoadCase { force: &force, weight: 1.0 }], &mut control).unwrap();
    (result.objective.compliance, result.objective.displacements[0].clone(), nodes, volume, fixed)
}

#[test]
fn real_cut_equilibrium_obeys_body_load_length_and_material_scaling() {
    let base = physical();
    let a = response(&base);
    let scale = 0.125_f64;
    let mut scaled = base.clone();
    scaled.bounds = ([0.25, -0.5, 1.0], [0.375, -0.375, 1.125]);
    scaled.height *= scale;
    scaled.curvature /= scale;
    scaled.radius *= scale;
    let b = response(&scaled);
    near(b.0, a.0 * scale.powi(5), 1e-6);
    near(b.3, a.3 * scale.powi(3), 1e-8);
    assert_eq!(a.1.len(), b.1.len());
    let norm = a.1.iter().copied().map(f64::abs).fold(0.0_f64, f64::max);
    for (u, v) in a.1.iter().zip(&b.1) {
        assert!((v / scale.powi(2) - u).abs() <= norm * 1e-6);
    }
    let mut stiff = scaled.clone();
    stiff.youngs *= 16.0;
    near(response(&stiff).0, b.0 / 16.0, 1e-7);
    let mut loaded = scaled;
    loaded.loads[0].0 = loaded.loads[0].0.map(|v| 2.0 * v);
    near(response(&loaded).0, 4.0 * b.0, 1e-7);
}

#[test]
fn translated_right_and_bottom_clamps_pin_the_actual_physical_nodes() {
    for face in [FixedFace::Right, FixedFace::Bottom] {
        let mut spec = physical();
        spec.bounds = ([0.25, -0.5, 1.0], [1.25, 0.5, 2.0]);
        spec.fixed = face;
        let (_, field, nodes, _, fixed) = response(&spec);
        assert!(fixed.iter().any(|v| *v));
        assert!(fixed.iter().any(|v| !v));
        for (i, (&p, &pin)) in nodes.iter().zip(&fixed).enumerate() {
            assert_eq!(pin, face.contains(p, spec.bounds));
            if pin { assert_eq!(&field[3 * i..3 * i + 3], &[0.0; 3]); }
        }
    }
}
