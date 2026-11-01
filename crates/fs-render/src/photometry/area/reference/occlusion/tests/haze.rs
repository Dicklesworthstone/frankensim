//! Independent two-dimensional surface integrals, not a second edge-integral oracle.
use super::*;
use crate::photometry::medium::{HomogeneousMedium, PhaseFunction};

// 70-digit Gauss-Legendre surface integration, independently split at the exact
// planar shadow/horizon boundaries. Orders 32 and 48 agreed within 2e-37 lx.
// Source area 4 m², nominal flux 1000 lm, height 2 m, extinction 1/4 m^-1:
// E = 1000/(4*pi) integral_visible [4 exp(-sqrt(x²+y²+4)/4)/(x²+y²+4)²] dx dy.
const CLEAR_LUX: f64 = 35.058_514_038_683_9;
const SHADOW_LUX: f64 = 15.435_136_245_708_191;
const OFF_AXIS_SHADOW_LUX: f64 = 19.623_377_792_975_71;
const VERTICAL_LUX: f64 = 4.038_488_724_789_262;

fn medium() -> HomogeneousMedium {
    HomogeneousMedium::try_new(0.125,0.125,PhaseFunction::HenyeyGreenstein {g:0.7}).unwrap()
}

#[test]
fn clear_panel_refines_extinction_without_any_shadow_ambiguity() {
    with_cx(|cx,_| {
        let primitives=[panel(-1.0,1.0,-1.0,1.0,2.0),panel(-4.0,4.0,-4.0,4.0,0.0)];
        let sources=[emitter()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap()
            .with_homogeneous_medium(&medium());
        let coarse=reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget {
            max_patch_evaluations:2,..AreaOcclusionBudget::default()
        }).unwrap();
        check(coarse,CLEAR_LUX*mass(&sources[0]));
        assert_eq!(coarse.termination,AreaOcclusionTermination::PatchBudget);
        assert_eq!(coarse.uncertain_visibility_patches,0);
        assert_eq!(coarse.varying_transmittance_patches,2);
        let fine=reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget::default()).unwrap();
        check(fine,CLEAR_LUX*mass(&sources[0]));
        assert_eq!(fine.termination,AreaOcclusionTermination::ToleranceMet);
        assert!(fine.evaluated_patches>2 && fine.evaluated_patches<=16_384);
        assert!(fine.numerical.lo>=coarse.numerical.lo && fine.numerical.hi<=coarse.numerical.hi);
        assert_eq!(fine.numerical_model_identity,coarse.numerical_model_identity);
        assert_eq!(fine.absorption_per_m,Some(0.125)); assert_eq!(fine.scattering_per_m,Some(0.125));
        let depth=reference.illuminance(cx,&probe(0.0,0.0),AreaOcclusionBudget {
            max_depth:0,..AreaOcclusionBudget::default()
        }).unwrap();
        assert_eq!(depth.termination,AreaOcclusionTermination::ResolutionLimit);
        check(depth,CLEAR_LUX*mass(&sources[0]));
    });
}

#[test]
fn penumbra_and_haze_share_one_error_account_not_two_tolerance_passes() {
    with_cx(|cx,_| {
        let primitives=strip_scene();let sources=[emitter()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap()
            .with_homogeneous_medium(&medium());
        for (x,truth) in [(0.0,SHADOW_LUX),(0.5,OFF_AXIS_SHADOW_LUX)] {
            let budget=AreaOcclusionBudget {absolute_error_lux:1.0,..AreaOcclusionBudget::default()};
            let value=reference.illuminance(cx,&probe(x,0.0),budget).unwrap();
            check(value,truth*mass(&sources[0]));
            assert_eq!(value.termination,AreaOcclusionTermination::ToleranceMet);
            assert!(value.uncertain_visibility_patches>0 && value.varying_transmittance_patches>0);
        }
    });
}

#[test]
fn vertical_receiver_keeps_horizon_clipping_and_complete_path_attenuation() {
    with_cx(|cx,_| {
        let primitives=[panel(-1.0,1.0,-1.0,1.0,2.0),panel(-4.0,4.0,-4.0,4.0,0.0)];
        let sources=[emitter()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap()
            .with_homogeneous_medium(&medium());
        let p=IlluminanceProbe::try_new(Point3::new(0.0,0.0,0.0),Vec3::new(1.0,0.0,0.0)).unwrap();
        let budget=AreaOcclusionBudget {max_patch_evaluations:1024,..AreaOcclusionBudget::default()};
        let value=reference.illuminance(cx,&p,budget).unwrap();
        check(value,VERTICAL_LUX*mass(&sources[0]));
        // The floor is not a blocker and the epsilon changes visibility only.
        let longer=reference.illuminance(cx,&p,AreaOcclusionBudget {ray_epsilon_m:0.25,..budget}).unwrap();
        assert_eq!(value.numerical,longer.numerical);
        assert_eq!(value.evaluated_patches,longer.evaluated_patches);
        assert_ne!(value.numerical_model_identity,longer.numerical_model_identity);
    });
}

#[test]
fn vacuum_and_phase_changes_preserve_the_intended_arithmetic_and_scope() {
    with_cx(|cx,_| {
        let primitives=[panel(-1.0,1.0,-1.0,1.0,2.0)];let sources=[emitter()];
        let make=||OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap();
        let p=probe(0.0,0.0);let budget=AreaOcclusionBudget {max_patch_evaluations:128,..AreaOcclusionBudget::default()};
        let original=make().illuminance(cx,&p,budget).unwrap();
        let vacuum=make().with_homogeneous_medium(&HomogeneousMedium::vacuum()).illuminance(cx,&p,budget).unwrap();
        assert_eq!(original,vacuum);assert_eq!(vacuum.varying_transmittance_patches,0);
        let isotropic=HomogeneousMedium::try_new(0.125,0.125,PhaseFunction::Isotropic).unwrap();
        let a=make().with_homogeneous_medium(&isotropic).illuminance(cx,&p,budget).unwrap();
        let b=make().with_homogeneous_medium(&medium()).illuminance(cx,&p,budget).unwrap();
        assert_eq!(a,b); // Phase affects the explicitly omitted in-scattered term only.
        let absorption=HomogeneousMedium::try_new(0.25,0.0,PhaseFunction::Isotropic).unwrap();
        let c=make().with_homogeneous_medium(&absorption).illuminance(cx,&p,budget).unwrap();
        assert_eq!(a.numerical,c.numerical);assert_eq!(a.evaluated_patches,c.evaluated_patches);
        assert_ne!(a.numerical_model_identity,c.numerical_model_identity);
        assert_ne!(a.numerical_model_identity,vacuum.numerical_model_identity);
        let reset=make().with_homogeneous_medium(&medium()).with_homogeneous_medium(&HomogeneousMedium::vacuum())
            .illuminance(cx,&p,budget).unwrap();
        assert_eq!(reset,vacuum); // Replacing a medium does not retain a stale bound or identity.
    });
}

#[test]
fn optically_thick_positive_light_never_becomes_a_false_exact_blackout() {
    with_cx(|cx,_| {
        let primitives=[panel(-1.0,1.0,-1.0,1.0,2.0)];let sources=[emitter()];
        let thick=HomogeneousMedium::try_new(512.0,0.0,PhaseFunction::Isotropic).unwrap();
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap()
            .with_homogeneous_medium(&thick);
        let budget=AreaOcclusionBudget {absolute_error_lux:0.0,max_patch_evaluations:2,..AreaOcclusionBudget::default()};
        let report=reference.illuminance(cx,&probe(0.0,0.0),budget).unwrap();
        assert_eq!(report.numerical.lo,0.0);
        assert!(report.numerical.hi>0.0 && report.numerical.hi<f64::MIN_POSITIVE);
        assert!(report.numerical_error_bound_lux>0.0);
        assert_eq!(report.termination,AreaOcclusionTermination::PatchBudget);
        assert!(report.total_error_bound_lux.is_none());
        let black=[emitter().filtered([0.0;SPECTRAL_BINS]).unwrap()];
        let report=OccludedAreaReference::new(&primitives,&black,AreaReferenceBudget::default(),cx).unwrap()
            .with_homogeneous_medium(&thick).illuminance(cx,&probe(0.0,0.0),budget).unwrap();
        assert_eq!(report.numerical.hi,0.0);assert_eq!(report.numerical_error_bound_lux,0.0);
        assert_eq!(report.varying_transmittance_patches,0);
    });
}

#[test]
fn rigid_placement_filters_and_cancel_resume_keep_haze_bounds_replayable() {
    with_cx(|cx,gate| {
        let placement=RigidTransform::try_new([0.0,0.0,1.0,0.0],[8.0,-4.0,2.0]).unwrap();
        let primitives:Vec<_>=strip_scene().into_iter().enumerate().map(|(i,p)| {
            let Shape::Mesh(mesh)=p.shape else {unreachable!()};
            Primitive {shape:Shape::Instance(GeometryInstance::try_new(i as u64+1,ContentHash([7;32]),
                SharedGeometry::mesh(mesh),placement).unwrap()),material:material(),emission:None}
        }).collect();
        let sources=[emitter().filtered([0.25;SPECTRAL_BINS]).unwrap()];
        let reference=OccludedAreaReference::new(&primitives,&sources,AreaReferenceBudget::default(),cx).unwrap()
            .with_homogeneous_medium(&medium());
        let probes=[IlluminanceProbe::try_new(Point3::new(8.0,-4.0,2.0),Vec3::new(0.0,0.0,1.0)).unwrap(),
            IlluminanceProbe::try_new(Point3::new(7.5,-4.0,2.0),Vec3::new(0.0,0.0,1.0)).unwrap()];
        let budget=AreaOcclusionBudget {max_patch_evaluations:128,..AreaOcclusionBudget::default()};
        let mut full=OccludedAreaJob::new(&reference,&probes,budget,cx).unwrap();full.advance(cx,usize::MAX).unwrap();
        check(full.results()[0],0.25*SHADOW_LUX*mass(&sources[0]));
        check(full.results()[1],0.25*OFF_AXIS_SHADOW_LUX*mass(&sources[0]));
        let mut split=OccludedAreaJob::new(&reference,&probes,budget,cx).unwrap();split.advance(cx,1).unwrap();
        gate.request();assert!(matches!(split.advance(cx,1),Err(AreaOcclusionError::Cancelled)));
        assert_eq!(split.results(),&full.results()[..1]);
        with_cx(|resume,_| {split.advance(resume,usize::MAX).unwrap();});
        assert_eq!(split.results(),full.results());
        let mut bytes=Vec::new();OccludedAreaIlluminance::write_csv_header(&mut bytes).unwrap();
        full.results()[0].write_csv_row(0,&mut bytes).unwrap();
        let csv=String::from_utf8(bytes).unwrap();let rows:Vec<_>=csv.lines().collect();
        let fields:Vec<_>=rows[1].split(',').collect();
        assert_eq!(rows[0].split(',').count(),18);assert_eq!(fields.len(),18);
        assert_eq!(fields[1],"opaque-mesh-homogeneous-gray-direct-area");
        assert_eq!(&fields[12..14],&["0.125","0.125"]);
        assert_eq!(fields[15],"unknown");assert_eq!(fields[16],"3");
    });
}

#[test]
fn native_absorbing_grid_path_estimate_is_checked_without_promoting_its_certificate() {
    use crate::photometry::area::PhotometricEmitters;
    use crate::photometry::path::{IlluminanceSampling,PathIlluminanceJob,PathMedium};
    use crate::photometry::TransportOptions;
    use crate::tracer::{Camera,Sampler,Scene};
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::{AbsorbingGrid,AbsorbingGridBudget,MassAbsorptionSpectrum};
    with_cx(|cx,_| {
        let scene=Scene {primitives:strip_scene(),lights:Vec::new(),environment:None,
            camera:Camera {eye:Point3::new(0.0,0.0,0.0),forward:Vec3::new(0.0,0.0,1.0),
                up:Vec3::new(0.0,1.0,0.0),half_tan:0.0}};
        let density=[1.0];
        let grid=AbsorbingGrid::try_new(cx,
            VolumeGrid::new([1,1,1],&density,[-8.0,-8.0,-1.0],[16.0,16.0,4.0]),
            MassAbsorptionSpectrum::table(&[(380.0,0.25),(780.0,0.25)]).unwrap(),None,
            AbsorbingGridBudget {max_grid_cells:1,max_steps_per_segment:128}).unwrap();
        let sources=[emitter()];let probes=[probe(0.0,0.0)];
        let gray=HomogeneousMedium::try_new(0.25,0.0,PhaseFunction::Isotropic).unwrap();
        let reference=OccludedAreaReference::new(&scene.primitives,&sources,AreaReferenceBudget::default(),cx).unwrap()
            .with_homogeneous_medium(&gray);
        let report=reference.illuminance(cx,&probes[0],AreaOcclusionBudget {
            absolute_error_lux:1.0,..AreaOcclusionBudget::default()
        }).unwrap();
        check(report,SHADOW_LUX*mass(&sources[0]));
        let mut paths=PathIlluminanceJob::new_with_emitters(&scene,
            PhotometricEmitters {points:&[],surfaces:&sources},PathMedium::Absorbing(&grid),&probes,
            TransportOptions::default(),IlluminanceSampling {samples_per_probe:32_768,max_depth:1,
                sampler:Sampler::OwenSobol,seed:173},cx).unwrap();
        paths.advance(cx,usize::MAX).unwrap();
        let estimate=paths.probe_estimate(0).unwrap();
        // Fixed authored regression allowance, NOT a confidence or sampling bound.
        assert!(estimate.illuminance_lux>=report.numerical.lo-1.0 && estimate.illuminance_lux<=report.numerical.hi+1.0);
        assert_eq!(estimate.numerical.kind,NumericalKind::NoClaim);
        assert!(estimate.total_error_bound_lux.is_none());
    });
}
