//! Independent area antiderivatives, shadow projections and original-owner checks.
use super::*;
use crate::instances::GeometryInstance;
use crate::photometry::source::{SourceSpectrum, SPECTRAL_BINS};
use crate::spectral::SpectralReflectance;
use fs_evidence::NumericalKind;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use std::f64::consts::PI;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 173, kernel_id: 0x4152_5348, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate))
}
fn material() -> Material {
    Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0).unwrap() }
}
fn panel(x0: f64, x1: f64, y0: f64, y1: f64, z: f64) -> Primitive {
    Primitive { shape: Shape::Mesh(TriMesh::new(
        vec![[x0,y0,z],[x1,y0,z],[x1,y1,z],[x0,y1,z]], vec![[0,2,1],[0,3,2]],
    )), material: material(), emission: None }
}
fn emitter() -> PhotometricAreaEmitter {
    PhotometricAreaEmitter::try_new(0, 1000.0, 1.0, SourceSpectrum::equal_energy(), EmissionSidedness::Front).unwrap()
}
fn probe(x: f64, y: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(x,y,0.0),Vec3::new(0.0,0.0,1.0)).unwrap()
}
fn mass(source: &PhotometricAreaEmitter) -> f64 {
    source.spectrum().bins_w_per_lumen_nm().iter().enumerate()
        .map(|(i,&p)| p*observer::xyz_at_bin(i)[1]*BIN_WIDTH_NM*LUMINOUS_EFFICACY_LM_W).sum()
}
fn rectangle(x0: f64, x1: f64, y0: f64, y1: f64, h: f64) -> f64 {
    // Independent double antiderivative of h²/(x²+y²+h²)².
    let f = |x: f64, y: f64| {
        let ax = (x*x+h*h).sqrt(); let by = (y*y+h*h).sqrt();
        0.5*(x/ax*(y/ax).atan()+y/by*(x/by).atan())
    };
    f(x1,y1)-f(x0,y1)-f(x1,y0)+f(x0,y0)
}
fn check(report: OccludedAreaIlluminance, expected: f64) {
    assert_eq!(report.numerical.kind, NumericalKind::Enclosure);
    assert!(report.numerical.lo <= expected && expected <= report.numerical.hi, "{report:?}, reference {expected}");
    assert!((report.lux-expected).abs() <= report.numerical_error_bound_lux);
    assert!(report.total_error_bound_lux.is_none());
    if report.termination == AreaOcclusionTermination::ToleranceMet {
        assert!(report.numerical_error_bound_lux <= report.requested_error_lux);
    }
}
fn strip_scene() -> Vec<Primitive> {
    vec![panel(-1.0,1.0,-1.0,1.0,2.0), panel(-0.25,0.25,-2.0,2.0,1.0), panel(-4.0,4.0,-4.0,4.0,0.0)]
}

#[test]
fn clear_panel_and_receiver_floor_keep_the_analytical_reference() {
    with_cx(|cx,_| {
        let primitives = [panel(-1.0,1.0,-1.0,1.0,2.0),panel(-4.0,4.0,-4.0,4.0,0.0)];
        let sources = [emitter()];
        let reference = OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let r = reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget::default()).unwrap();
        check(r,1000.0/(4.0*PI)*mass(&sources[0])*rectangle(-1.0,1.0,-1.0,1.0,2.0));
        assert_eq!(r.termination,AreaOcclusionTermination::ToleranceMet);
        assert_eq!(r.evaluated_patches,2); assert_eq!(r.uncertain_visibility_patches,0);
        assert!(r.numerical_error_bound_lux<1e-8);
    });
}

#[test]
fn definite_blocker_is_zero_without_sampling_a_shadow_ray() {
    with_cx(|cx,_| {
        // One triangle contains the entire projected aperture strictly inside.
        let blocker = Primitive { shape:Shape::Mesh(TriMesh::new(
            vec![[-10.0,-10.0,1.0],[10.0,-10.0,1.0],[0.0,10.0,1.0]],vec![[0,1,2]])),material:material(),emission:None };
        let primitives = [panel(-1.0,1.0,-1.0,1.0,2.0),blocker]; let sources=[emitter()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let r=reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget {absolute_error_lux:0.0,..AreaOcclusionBudget::default()}).unwrap();
        check(r,0.0);assert_eq!(r.numerical_error_bound_lux,0.0);
        assert_eq!(r.termination,AreaOcclusionTermination::ToleranceMet);
        assert_eq!(r.evaluated_patches,2);
    });
}

#[test]
fn partial_strip_shadow_encloses_independent_visible_rectangle_integrals() {
    with_cx(|cx,_| {
        let primitives=strip_scene();let sources=[emitter()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let expected=1000.0/(4.0*PI)*mass(&sources[0])
            *(rectangle(-1.0,1.0,-1.0,1.0,2.0)-rectangle(-0.5,0.5,-1.0,1.0,2.0));
        let p=probe(0.0,0.0);
        let coarse=reference.illuminance(cx,&p,AreaOcclusionBudget {max_patch_evaluations:2,..AreaOcclusionBudget::default()}).unwrap();
        check(coarse,expected);assert_eq!(coarse.termination,AreaOcclusionTermination::PatchBudget);
        let fine=reference.illuminance(cx,&p,AreaOcclusionBudget {absolute_error_lux:1.0,..AreaOcclusionBudget::default()}).unwrap();
        check(fine,expected);assert_eq!(fine.termination,AreaOcclusionTermination::ToleranceMet);
        assert!(fine.numerical.lo>=coarse.numerical.lo && fine.numerical.hi<=coarse.numerical.hi);
        assert!(fine.uncertain_visibility_patches>0); // Small unresolved penumbrae remain included.
        assert!(fine.evaluated_patches>2 && fine.evaluated_patches<=16_384);
        assert_eq!(coarse.numerical_model_identity,fine.numerical_model_identity);
        let depth=reference.illuminance(cx,&p,AreaOcclusionBudget {max_depth:0,..AreaOcclusionBudget::default()}).unwrap();
        check(depth,expected);assert_eq!(depth.termination,AreaOcclusionTermination::ResolutionLimit);
    });
}

#[test]
fn different_planes_of_one_emitting_primitive_still_occlude_each_other() {
    with_cx(|cx,_| {
        // Two equal aligned panels in ONE mesh. The nearer one covers every
        // direction to the farther one. Excluding the whole source is unsound.
        let primitives=[Primitive {shape:Shape::Mesh(TriMesh::new(
            vec![[-1.0,-1.0,2.0],[1.0,-1.0,2.0],[1.0,1.0,2.0],[-1.0,1.0,2.0],
                 [-1.0,-1.0,1.0],[1.0,-1.0,1.0],[1.0,1.0,1.0],[-1.0,1.0,1.0]],
            vec![[0,2,1],[0,3,2],[4,6,5],[4,7,6]])),material:material(),emission:None}];
        let sources=[emitter()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let r=reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget {absolute_error_lux:2.0,..AreaOcclusionBudget::default()}).unwrap();
        let expected=1000.0/(8.0*PI)*mass(&sources[0])*rectangle(-1.0,1.0,-1.0,1.0,1.0);
        check(r,expected);assert!(r.numerical_error_bound_lux<3.0);
        let unoccluded=reference.sources.illuminance(cx,&probe(0.0,0.0)).unwrap();
        assert!(r.numerical.hi<unoccluded.numerical.lo);
    });
}

#[test]
fn rigid_instances_and_passive_source_filters_retain_original_owners() {
    with_cx(|cx,_| {
        let transform=RigidTransform::try_new([0.0,0.0,1.0,0.0],[8.0,-4.0,2.0]).unwrap();
        let original=strip_scene();
        let primitives:Vec<_>=original.into_iter().enumerate().map(|(i,p)| {
            let Shape::Mesh(mesh)=p.shape else {unreachable!()};
            Primitive {shape:Shape::Instance(GeometryInstance::try_new(i as u64+1,ContentHash([7;32]),SharedGeometry::mesh(mesh),transform).unwrap()),material:material(),emission:None}
        }).collect();
        let sources=[emitter().filtered([0.25;SPECTRAL_BINS]).unwrap()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let p=IlluminanceProbe::try_new(Point3::new(8.0,-4.0,2.0),Vec3::new(0.0,0.0,1.0)).unwrap();
        let r=reference.illuminance(cx,&p,AreaOcclusionBudget {absolute_error_lux:1.0,..AreaOcclusionBudget::default()}).unwrap();
        let expected=250.0/(4.0*PI)*mass(&sources[0])*(rectangle(-1.0,1.0,-1.0,1.0,2.0)-rectangle(-0.5,0.5,-1.0,1.0,2.0));
        check(r,expected);
        let changed=[panel(-1.0,1.0,-1.0,1.0,2.0),panel(-0.4,0.4,-2.0,2.0,1.0)];
        let changed=OccludedAreaReference::new(&changed,&sources,AreaReferenceBudget::default(),cx).unwrap()
            .illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget {max_patch_evaluations:2,..AreaOcclusionBudget::default()}).unwrap();
        assert_ne!(r.numerical_model_identity,changed.numerical_model_identity);
    });
}

#[test]
fn cancellation_and_chunking_preserve_the_published_prefix() {
    with_cx(|cx,gate| {
        let primitives=strip_scene();let sources=[emitter()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let grid=CalculationGrid::try_new(Point3::new(-0.25,0.0,0.0),Vec3::new(0.25,0.0,0.0),Vec3::new(0.0,1.0,0.0),3,1).unwrap();
        let budget=AreaOcclusionBudget {max_patch_evaluations:128,..AreaOcclusionBudget::default()};
        let mut full=OccludedAreaJob::for_grid(&reference,grid,budget,cx).unwrap();full.advance(cx,usize::MAX).unwrap();
        let mut split=OccludedAreaJob::for_grid(&reference,grid,budget,cx).unwrap();split.advance(cx,1).unwrap();
        gate.request();assert!(matches!(split.advance(cx,1),Err(AreaOcclusionError::Cancelled)));
        assert_eq!(split.results(),&full.results()[..1]);
        with_cx(|resume,_| {
            while !split.is_complete(){split.advance(resume,1).unwrap();}
            assert_eq!(split.results(),full.results());
            let probes=[probe(0.0,0.0),IlluminanceProbe::try_new(Point3::new(0.0,0.0,2.0),Vec3::new(0.0,0.0,1.0)).unwrap()];
            let mut failing=OccludedAreaJob::new(&reference,&probes,budget,resume).unwrap();
            assert!(matches!(failing.advance(resume,2),Err(AreaOcclusionError::Reference(AreaReferenceError::UnresolvedSourcePlane{..}))));
            assert_eq!(failing.results().len(),1);
        });
    });
}

#[test]
fn numerical_goals_work_limits_and_output_errors_stay_explicit() {
    with_cx(|cx,_| {
        let primitives=strip_scene();let sources=[emitter()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        for budget in [AreaOcclusionBudget {absolute_error_lux:f64::NAN,..AreaOcclusionBudget::default()},
            AreaOcclusionBudget {absolute_error_lux:-1.0,..AreaOcclusionBudget::default()},
            AreaOcclusionBudget {ray_epsilon_m:0.0,..AreaOcclusionBudget::default()},
            AreaOcclusionBudget {max_depth:33,..AreaOcclusionBudget::default()},
            AreaOcclusionBudget {max_patch_evaluations:MAX_AREA_SHADOW_PATCHES+1,..AreaOcclusionBudget::default()}] {
            assert!(matches!(reference.illuminance(cx,&probe(0.0,0.0),budget),Err(AreaOcclusionError::InvalidOptions)));
        }
        assert!(matches!(reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget {max_patch_evaluations:1,..AreaOcclusionBudget::default()}),Err(AreaOcclusionError::Budget)));
        let r=reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget {max_patch_evaluations:2,..AreaOcclusionBudget::default()}).unwrap();
        let mut csv=Vec::new();OccludedAreaIlluminance::write_csv_header(&mut csv).unwrap();r.write_csv_row(0,&mut csv).unwrap();
        let csv=String::from_utf8(csv).unwrap();let rows:Vec<_>=csv.lines().collect();
        assert_eq!(rows[0].split(',').count(),17);assert_eq!(rows[1].split(',').count(),17);
        assert!(rows[1].contains(",patch-budget,"));assert!(rows[1].ends_with(",unknown,2"));
        struct Refuse;impl Write for Refuse {fn write(&mut self,_:&[u8])->io::Result<usize>{Err(io::Error::other("test sink"))}fn flush(&mut self)->io::Result<()>{Ok(())}}
        assert!(r.write_csv_row(0,&mut Refuse).is_err());
    });
}

#[test]
fn source_blackout_is_exact_without_promoting_physical_accuracy() {
    with_cx(|cx,_| {
        let primitives=strip_scene();let sources=[emitter().filtered([0.0;SPECTRAL_BINS]).unwrap()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let r=reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget {absolute_error_lux:0.0,..AreaOcclusionBudget::default()}).unwrap();
        check(r,0.0);assert_eq!(r.termination,AreaOcclusionTermination::ToleranceMet);
    });
}

#[test]
fn native_shadowed_area_estimator_remains_a_separate_sampling_claim() {
    use crate::photometry::area::PhotometricEmitters;
    use crate::photometry::path::{IlluminanceSampling, PathIlluminanceJob, PathMedium};
    use crate::photometry::TransportOptions;
    use crate::tracer::{Camera, Sampler, Scene};
    with_cx(|cx,_| {
        let scene=Scene {primitives:strip_scene(),lights:Vec::new(),environment:None,
            camera:Camera {eye:Point3::new(0.0,0.0,0.0),forward:Vec3::new(0.0,0.0,1.0),up:Vec3::new(0.0,1.0,0.0),half_tan:0.0}};
        let sources=[emitter()];let probes=[probe(0.0,0.0)];
        let reference=OccludedAreaReference::new(&scene.primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let r=reference.illuminance(cx,&probes[0],AreaOcclusionBudget {absolute_error_lux:1.0,..AreaOcclusionBudget::default()}).unwrap();
        let mut paths=PathIlluminanceJob::new_with_emitters(&scene,PhotometricEmitters {points:&[],surfaces:&sources},PathMedium::Vacuum,
            &probes,TransportOptions::default(),IlluminanceSampling {samples_per_probe:32_768,max_depth:1,sampler:Sampler::OwenSobol,seed:173},cx).unwrap();
        paths.advance(cx,usize::MAX).unwrap();
        let estimate=paths.probe_estimate(0).unwrap();
        // Authored regression allowance, NOT a probability or sampling bound.
        assert!(estimate.illuminance_lux>=r.numerical.lo-1.0 && estimate.illuminance_lux<=r.numerical.hi+1.0);
        assert_eq!(estimate.numerical.kind,NumericalKind::NoClaim);
        assert!(estimate.total_error_bound_lux.is_none());
    });
}

mod haze;
