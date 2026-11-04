//! Independent area-integral references and original-owner integration checks.
use super::*;
use crate::instances::GeometryInstance;
use crate::photometry::area::PhotometricEmitters;
use crate::photometry::path::{IlluminanceSampling, PathIlluminanceJob, PathMedium};
use crate::photometry::source::{SourceSpectrum, SPECTRAL_BINS};
use crate::photometry::TransportOptions;
use crate::spectral::SpectralReflectance;
use crate::tracer::{Camera, Sampler, Scene};
use fs_evidence::NumericalKind;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use std::f64::consts::PI;

pub(super) fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 173, kernel_id: 0x4152_4546, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate))
}
pub(super) fn material() -> Material {
    Material::SpectralDiffuse { reflectance: SpectralReflectance::constant(0.0).unwrap() }
}
pub(super) fn panel(a: f64, b: f64, h: f64, front: bool) -> Primitive {
    Primitive { shape: Shape::Mesh(TriMesh::new(
        vec![[-a,-b,h],[a,-b,h],[a,b,h],[-a,b,h]],
        if front { vec![[0,2,1],[0,3,2]] } else { vec![[0,1,2],[0,2,3]] },
    )), material: material(), emission: None }
}
fn source(sides: EmissionSidedness) -> PhotometricAreaEmitter {
    PhotometricAreaEmitter::try_new(0,1000.0,1.0,SourceSpectrum::equal_energy(),sides).unwrap()
}
pub(super) fn probe(p: [f64;3], n: [f64;3]) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(p[0],p[1],p[2]),Vec3::new(n[0],n[1],n[2])).unwrap()
}
pub(super) fn mass(source: &PhotometricAreaEmitter) -> f64 {
    source.spectrum().bins_w_per_lumen_nm().iter().enumerate()
        .map(|(i,&p)| p*observer::xyz_at_bin(i)[1]*BIN_WIDTH_NM*LUMINOUS_EFFICACY_LM_W).sum()
}
pub(super) fn check(r: AreaReferenceIlluminance, expected: f64) {
    assert_eq!(r.numerical.kind,NumericalKind::Enclosure);
    assert!(r.numerical.lo <= expected && r.numerical.hi >= expected,"{r:?}; expected={expected}");
    assert!((r.lux-expected).abs() <= r.numerical_error_bound_lux);
    assert!(r.total_error_bound_lux.is_none());
}
// Independent antiderivative of h²/(x²+y²+h²)² over the source rectangle.
pub(super) fn rectangle_factor(a: f64,b: f64,h: f64,x: f64,y: f64) -> f64 {
    let f=|x:f64,y:f64| {
        let ax=(x*x+h*h).sqrt(); let by=(y*y+h*h).sqrt();
        0.5*(x/ax*(y/ax).atan()+y/by*(x/by).atan())
    };
    f(a-x,b-y)-f(-a-x,b-y)-f(a-x,-b-y)+f(-a-x,-b-y)
}

#[test]
fn rectangle_near_off_axis_and_far_field_enclose_an_independent_area_integral() {
    with_cx(|cx,_| {
        for h in [0.125,1.0,3.0,64.0] {
            let primitives=[panel(0.5,0.75,h,true)]; let sources=[source(EmissionSidedness::Front)];
            let reference=UnoccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
            for [x,y] in [[0.0,0.0],[0.25,-0.125],[1.0,0.5]] {
                let r=reference.illuminance(cx,&probe([x,y,0.0],[0.0,0.0,1.0])).unwrap();
                check(r,1000.0/(1.5*PI)*mass(&sources[0])*rectangle_factor(0.5,0.75,h,x,y));
                assert_eq!(r.coarse_triangles,0);
                assert!(r.numerical_error_bound_lux < 1e-7);
            }
        }
    });
}

#[test]
fn horizon_clipping_keeps_the_positive_half_of_a_panel() {
    with_cx(|cx,_| {
        let primitives=[panel(1.0,1.0,2.0,true)]; let sources=[source(EmissionSidedness::Front)];
        let reference=UnoccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let r=reference.illuminance(cx,&probe([0.0;3],[1.0,0.0,0.0])).unwrap();
        // Integrate 2*x/(x²+y²+4)² over x in [0,1], y in [-1,1].
        // Integrating the unclipped polygon first would incorrectly yield zero.
        let s=5.0_f64.sqrt();
        let factor=2.0*((0.5_f64).atan()/2.0-(1.0/s).atan()/s);
        check(r,1000.0/(4.0*PI)*mass(&sources[0])*factor);
        assert!(r.lux>1.0);
        assert_eq!(r.coarse_triangles,0);
        assert!(r.numerical_error_bound_lux<1e-8);
    });
}

#[test]
fn original_sidedness_black_sources_and_passive_spectra_are_preserved() {
    with_cx(|cx,_| {
        let primitives=[panel(1.0,1.0,2.0,false)];
        let p=probe([0.0;3],[0.0,0.0,1.0]);
        let front=[source(EmissionSidedness::Front)];
        let r=UnoccludedAreaReference::new(&primitives,&front,AreaReferenceBudget::default(),cx).unwrap().illuminance(cx,&p).unwrap();
        check(r,0.0); assert_eq!(r.numerical_error_bound_lux,0.0);
        let both=[source(EmissionSidedness::Both)];
        let r=UnoccludedAreaReference::new(&primitives,&both,AreaReferenceBudget::default(),cx).unwrap().illuminance(cx,&p).unwrap();
        check(r,1000.0/(8.0*PI)*mass(&both[0])*rectangle_factor(1.0,1.0,2.0,0.0,0.0));
        let filtered=[both[0].filtered([0.25;SPECTRAL_BINS]).unwrap()];
        let f=UnoccludedAreaReference::new(&primitives,&filtered,AreaReferenceBudget::default(),cx).unwrap().illuminance(cx,&p).unwrap();
        assert!((f.lux/r.lux-0.25).abs()<1e-12);
        assert_ne!(r.numerical_model_identity,f.numerical_model_identity);
        let black=[both[0].filtered([0.0;SPECTRAL_BINS]).unwrap()];
        // An exactly off source does not require a defined near-field direction.
        let r=UnoccludedAreaReference::new(&primitives,&black,AreaReferenceBudget::default(),cx).unwrap().illuminance(cx,&probe([1.0,1.0,2.0],[0.0,0.0,1.0])).unwrap();
        check(r,0.0);
    });
}

#[test]
fn rigid_instance_reference_uses_actual_geometry_not_its_caller_supplied_hash() {
    with_cx(|cx,_| {
        let transform=RigidTransform::try_new([0.0,0.0,1.0,0.0],[8.0,-4.0,2.0]).unwrap();
        let sources=[source(EmissionSidedness::Front)];
        let make=|a| {
            let Shape::Mesh(mesh)=panel(a,0.75,3.0,true).shape else { unreachable!() };
            [Primitive { shape: Shape::Instance(GeometryInstance::try_new(1,ContentHash([7;32]),
                SharedGeometry::mesh(mesh),transform).unwrap()), material:material(),emission:None }]
        };
        let p=probe([8.0,-4.0,2.0],[0.0,0.0,1.0]);
        let shapes=make(0.5);
        let r=UnoccludedAreaReference::new(&shapes,&sources,AreaReferenceBudget::default(),cx).unwrap().illuminance(cx,&p).unwrap();
        check(r,1000.0/(1.5*PI)*mass(&sources[0])*rectangle_factor(0.5,0.75,3.0,0.0,0.0));
        let changed=make(1.0);
        let changed=UnoccludedAreaReference::new(&changed,&sources,AreaReferenceBudget::default(),cx).unwrap().illuminance(cx,&p).unwrap();
        assert_ne!(r.numerical_model_identity,changed.numerical_model_identity);
    });
}

#[test]
fn native_path_estimator_can_be_compared_without_promoting_its_claim() {
    with_cx(|cx,_| {
        let scene=Scene { primitives:vec![panel(1.0,1.0,2.0,true)],lights:Vec::new(),environment:None,
            camera:Camera {eye:Point3::new(0.0,0.0,0.0),forward:Vec3::new(0.0,0.0,1.0),up:Vec3::new(0.0,1.0,0.0),half_tan:0.0} };
        let sources=[source(EmissionSidedness::Front)]; let probes=[probe([0.2,0.1,0.0],[0.0,0.0,1.0])];
        let reference=UnoccludedAreaReference::new(&scene.primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let r=reference.illuminance(cx,&probes[0]).unwrap();
        let mut job=PathIlluminanceJob::new_with_emitters(&scene,PhotometricEmitters {points:&[],surfaces:&sources},PathMedium::Vacuum,
            &probes,TransportOptions::default(),IlluminanceSampling {samples_per_probe:16_384,max_depth:1,sampler:Sampler::OwenSobol,seed:173},cx).unwrap();
        job.advance(cx,usize::MAX).unwrap();
        let estimate=job.probe_estimate(0).unwrap();
        // This is a regression tolerance, NOT a sampling confidence statement.
        assert!((estimate.illuminance_lux-r.lux).abs()<0.02*r.lux);
        assert_eq!(estimate.numerical.kind,NumericalKind::NoClaim);
        assert!(estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn cancellation_grid_partitioning_and_failed_receivers_keep_the_accepted_prefix() {
    with_cx(|cx,gate| {
        let shapes=[panel(0.5,0.5,2.0,true)]; let sources=[source(EmissionSidedness::Front)];
        let reference=UnoccludedAreaReference::new(&shapes,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let grid=CalculationGrid::try_new(Point3::new(-0.5,-0.5,0.0),Vec3::new(0.5,0.0,0.0),Vec3::new(0.0,0.5,0.0),3,3).unwrap();
        let mut full=UnoccludedAreaJob::for_grid(&reference,grid,cx).unwrap();full.advance(cx,usize::MAX).unwrap();
        let mut split=UnoccludedAreaJob::for_grid(&reference,grid,cx).unwrap();split.advance(cx,2).unwrap();
        gate.request();assert!(matches!(split.advance(cx,5),Err(AreaReferenceError::Cancelled)));
        assert_eq!(split.results(),&full.results()[..2]);
        with_cx(|resume,_| {
            while !split.is_complete(){split.advance(resume,2).unwrap();}
            assert_eq!(split.results(),full.results());
            let probes=[probe([0.0;3],[0.0,0.0,1.0]),probe([0.0,0.0,2.0],[0.0,0.0,1.0])];
            let mut job=UnoccludedAreaJob::new(&reference,&probes,resume).unwrap();
            assert!(matches!(job.advance(resume,2),Err(AreaReferenceError::UnresolvedSourcePlane {..})));
            assert_eq!(job.results().len(),1);assert!(!job.is_complete());
        });
    });
}

#[test]
fn storage_limits_duplicate_faces_and_duplicate_bindings_refuse() {
    with_cx(|cx,_| {
        let shapes=[panel(1.0,1.0,2.0,true)];let sources=[source(EmissionSidedness::Front)];
        for budget in [AreaReferenceBudget {max_emitters:0,..AreaReferenceBudget::default()},AreaReferenceBudget {max_triangles:1,..AreaReferenceBudget::default()}] {
            assert!(matches!(UnoccludedAreaReference::new(&shapes,&sources,budget,cx),Err(AreaReferenceError::Budget)));
        }
        let duplicate=[sources[0].clone(),sources[0].clone()];
        assert!(matches!(UnoccludedAreaReference::new(&shapes,&duplicate,AreaReferenceBudget::default(),cx),Err(AreaReferenceError::InvalidGeometry)));
        let duplicate=[Primitive {shape:Shape::Mesh(TriMesh::new(vec![[0.0,0.0,1.0],[1.0,0.0,1.0],[0.0,1.0,1.0]],vec![[0,1,2],[2,1,0]])),material:material(),emission:None}];
        assert!(matches!(UnoccludedAreaReference::new(&duplicate,&sources,AreaReferenceBudget::default(),cx),Err(AreaReferenceError::InvalidGeometry)));
        let reference=UnoccludedAreaReference::new(&shapes,&sources,AreaReferenceBudget {max_outputs:0,..AreaReferenceBudget::default()},cx).unwrap();
        assert!(matches!(UnoccludedAreaJob::new(&reference,&[probe([0.0;3],[0.0,0.0,1.0])],cx),Err(AreaReferenceError::Budget)));
    });
}

#[test]
fn a_numerically_ambiguous_horizon_keeps_a_conservative_range() {
    with_cx(|cx,_| {
        let shapes=[panel(1.0,1.0,2.0,true)];let sources=[source(EmissionSidedness::Front)];
        let reference=UnoccludedAreaReference::new(&shapes,&sources,AreaReferenceBudget::default(),cx).unwrap();
        // Two corners lie on x+y=0; their rounded dot intervals cannot decide
        // a sign. No sample at a guessed side substitutes for the full range.
        let r=reference.illuminance(cx,&probe([0.0;3],[1.0,1.0,0.0])).unwrap();
        assert!(r.coarse_triangles>0);assert_eq!(r.numerical.lo,0.0);assert!(r.numerical.hi>0.0);
        assert!(r.total_error_bound_lux.is_none());
    });
}

#[test]
fn unselected_blockers_are_explicitly_outside_this_reference_scope() {
    with_cx(|cx,_| {
        let sources=[source(EmissionSidedness::Front)];let p=probe([0.0;3],[0.0,0.0,1.0]);
        let shapes=[panel(1.0,1.0,2.0,true)];
        let r=UnoccludedAreaReference::new(&shapes,&sources,AreaReferenceBudget::default(),cx).unwrap().illuminance(cx,&p).unwrap();
        let shapes=[panel(1.0,1.0,2.0,true),panel(10.0,10.0,1.0,true)];
        let unchanged=UnoccludedAreaReference::new(&shapes,&sources,AreaReferenceBudget::default(),cx).unwrap().illuminance(cx,&p).unwrap();
        assert_eq!(r,unchanged);assert!(r.lux>0.0); // Not a full-scene shadow certificate.
        let mut csv=Vec::new();AreaReferenceIlluminance::write_csv_header(&mut csv).unwrap();r.write_csv_row(0,&mut csv).unwrap();
        let csv=String::from_utf8(csv).unwrap();let rows:Vec<_>=csv.lines().collect();
        assert_eq!(rows[0].split(',').count(),13);assert_eq!(rows[1].split(',').count(),13);
        assert!(rows[1].contains(",unoccluded-vacuum-direct-area,"));
        assert!(rows[1].ends_with(",unknown,unknown,unknown,unknown,1"));
        struct Refuse;impl Write for Refuse {fn write(&mut self,_:&[u8])->io::Result<usize>{Err(io::Error::other("test sink"))}fn flush(&mut self)->io::Result<()>{Ok(())}}
        assert!(r.write_csv_row(0,&mut Refuse).is_err());
    });
}

#[test]
fn black_emitting_cavity_closes_the_independent_hemisphere_energy_identity() {
    with_cx(|cx,_| {
        // Six 4 m² inward faces. Every interior direction terminates on one
        // uniformly bright face, hence E=pi*L=flux/24 before any reflection.
        let primitives=[Primitive { shape:Shape::Mesh(TriMesh::new(
            vec![[-1.0,-1.0,-1.0],[1.0,-1.0,-1.0],[1.0,1.0,-1.0],[-1.0,1.0,-1.0],
                 [-1.0,-1.0,1.0],[1.0,-1.0,1.0],[1.0,1.0,1.0],[-1.0,1.0,1.0]],
            vec![[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,3,7],[0,7,4],
                 [1,5,6],[1,6,2],[0,4,5],[0,5,1],[3,2,6],[3,6,7]],
        )),material:material(),emission:None }];
        let sources=[source(EmissionSidedness::Front)];
        let reference=UnoccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        for (p,n) in [([0.0;3],[0.0,0.0,1.0]),([0.25,-0.125,0.5],[0.2,0.35,1.0])] {
            let r=reference.illuminance(cx,&probe(p,n)).unwrap();
            check(r,sources[0].flux_lm()/24.0*mass(&sources[0]));
            assert_eq!(r.coarse_triangles,0);
            assert!(r.numerical_error_bound_lux<1e-7);
        }
    });
}
