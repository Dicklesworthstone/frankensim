//! Independent line integrals and original-owner integration checks.
use super::*;
use crate::charts::{Ray, TriMesh};
use crate::instances::RigidTransform;
use crate::photometry::{TransportOptions, medium::{HomogeneousMedium, ScatterSettings}};
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::spectral::SpectralReflectance;
use crate::tracer::{Material, Primitive, Shape};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Vec3;
use core::f64::consts::PI;

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate)->T)->T {
    let gate=CancelGate::new_clock_free();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey{seed:173,kernel_id:0x5343_454e,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        f(&cx,&gate)
    })
}
fn fixture(position:[f64;3])->PhotometricSource {
    PhotometricSource::try_new(RigidTransform::try_new([0.0,0.0,0.0,1.0],position).unwrap(),
        400.0*PI,1.0,AngularDistribution::isotropic(),SourceSpectrum::equal_energy()).unwrap()
}
fn segment()->ScatterSegment {
    ScatterSegment::try_new(Point3::new(0.0,0.0,0.0),Point3::new(0.0,0.0,2.0)).unwrap()
}
fn budget()->ScatterEnclosureBudget {
    ScatterEnclosureBudget { absolute_error_cd_m2:1e-4,max_evaluations:16_384,max_depth:40 }
}
fn contains(value:EnclosedSingleScatter,truth:f64) {
    assert!(value.numerical.lo<=truth && value.numerical.hi>=truth,"{value:?} misses {truth}");
    assert!((value.luminance_cd_m2-truth).abs()<=value.numerical_error_bound_cd_m2);
    assert_eq!(value.numerical.kind,fs_evidence::NumericalKind::Enclosure);
    assert!(value.total_error_bound_cd_m2.is_none());
    assert!(value.omitted_multiple_scattering_bound_cd_m2.is_none());
}
fn opaque(vertices:Vec<[f64;3]>,triangles:Vec<[u32;3]>)->Primitive {
    Primitive { shape:Shape::Mesh(TriMesh::new(vertices,triangles)),
        material:Material::SpectralDiffuse{reflectance:SpectralReflectance::constant(0.0).unwrap()},emission:None }
}

#[test]
fn collinear_airlight_has_an_independent_closed_form_and_native_comparison() {
    with_cx(|cx,_| {
        let sources=[fixture([0.0,0.0,5.0])];
        let scene=PhotometricScene::try_new(&[],&sources,
            HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        let reference=HomogeneousScatterReference::new(&scene,cx).unwrap();
        let value=reference.enclose(cx,segment(),budget()).unwrap();
        // r=5-s, so BOTH-leg optical distance is constant 5; integral ds/r².
        let truth=20.0/(4.0*PI)*(-1.5_f64).exp()*(1.0/3.0-1.0/5.0);
        contains(value,truth);
        assert_eq!(value.termination,ScatterEnclosureTermination::ToleranceMet);
        assert_eq!(value.uncertain_visibility_panels,0);
        let native=scene.luminance_along_ray(cx,&Ray{origin:segment().camera_m(),dir:Vec3::new(0.0,0.0,1.0)},2.0,
            ScatterSettings { initial_panels:32,max_panels:16_384,absolute_tolerance_cd_m2:1e-10,relative_tolerance:1e-9 }).unwrap();
        contains(value,native.single_scatter.luminance_cd_m2);
        assert_eq!(native.numerical.kind,fs_evidence::NumericalKind::NoClaim);
    });
}

#[test]
fn original_forward_backward_rayleigh_and_tabulated_phases_keep_their_meaning() {
    with_cx(|cx,_| {
        let sources=[fixture([0.0,0.0,5.0])];
        let tabulated=PhaseFunction::tabulated(&[(-1.0,1.0),(0.0,1.0),(1.0,1.0)]).unwrap();
        let cases=[(PhaseFunction::Isotropic,1.0/(4.0*PI)),
            (PhaseFunction::HenyeyGreenstein{g:0.5},6.0/(4.0*PI)),
            (PhaseFunction::HenyeyGreenstein{g:-0.5},(2.0/9.0)/(4.0*PI)),
            (PhaseFunction::Rayleigh,1.5/(4.0*PI)),
            (PhaseFunction::RayleighDepolarized{depolarization:0.2},(3.0/2.2)/(4.0*PI)),
            (tabulated,1.0/(4.0*PI))];
        for (phase,p) in cases {
            let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.1,0.2,phase).unwrap(),TransportOptions::default(),cx).unwrap();
            let result=HomogeneousScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
            contains(result,20.0*p*(-1.5_f64).exp()*(1.0/3.0-1.0/5.0));
        }
        // A narrow interior knot peak must survive even when both endpoints are dark.
        let phase=PhaseFunction::tabulated(&[(-1.0,0.0),(-0.01,0.0),(0.0,1.0),(0.01,0.0),(1.0,0.0)]).unwrap();
        let bound=phase::enclose(&phase,Interval::new(-0.1,0.1),cx).unwrap();
        assert_eq!(bound.lo(),0.0);
        assert!(bound.hi()>=phase.eval_cos(0.0).unwrap());
        for g in [fs_math::next_down(1.0),-fs_math::next_down(1.0)] {
            let phase=PhaseFunction::HenyeyGreenstein{g};
            let peak=phase::enclose(&phase,Interval::point(g.signum()),cx).unwrap();
            assert!(finite(peak) && peak.lo()>0.0);
        }
    });
}

#[test]
fn narrow_fixture_support_is_not_lost_between_samples() {
    with_cx(|cx,_| {
        let q=0.5_f64.sqrt();
        let law=AngularDistribution::uniform_cone(0.01).unwrap();
        let sources=[PhotometricSource::try_new(
            RigidTransform::try_new([0.0,-q,0.0,q],[1.0,0.0,0.731]).unwrap(),
            100.0/law.peak_per_steradian(),1.0,law,SourceSpectrum::equal_energy()).unwrap()];
        let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        let reference=HomogeneousScatterReference::new(&scene,cx).unwrap();
        let coarse=reference.enclose(cx,segment(),ScatterEnclosureBudget{max_evaluations:1,..budget()}).unwrap();
        // Independent 80-digit quadrature between the geometric cone intersections.
        let truth=0.018_937_362_075_047_183;
        contains(coarse,truth);
        assert_eq!(coarse.termination,ScatterEnclosureTermination::PanelBudget);
        let refined=reference.enclose(cx,segment(),budget()).unwrap();
        contains(refined,truth);
        assert!(refined.numerical.lo>0.0);
        assert_eq!(refined.termination,ScatterEnclosureTermination::ToleranceMet);
    });
}

#[test]
fn light_shadows_and_camera_blockage_bound_the_complete_line_integral() {
    with_cx(|cx,_| {
        let sources=[fixture([2.0,0.0,0.0])];
        let light_screen=[opaque(vec![[1.0,-1.0,0.2],[1.0,1.0,0.2],[1.0,1.0,0.6],[1.0,-1.0,0.6]],vec![[0,1,2],[0,2,3]])];
        let camera_screen=[opaque(vec![[-1.0,-1.0,1.0],[1.0,-1.0,1.0],[0.0,1.0,1.0]],vec![[0,1,2]])];
        // Independent point-light airlight quadrature over explicitly projected
        // visible intervals: [0,.4] U [1.2,2], and [0,1], respectively.
        for (geometry,truth) in [(&light_screen[..],0.137_851_370_774_233_26),(&camera_screen[..],0.172_288_032_886_401_63)] {
            let scene=PhotometricScene::try_new(geometry,&sources,HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
            let result=HomogeneousScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
            contains(result,truth);
            assert_eq!(result.termination,ScatterEnclosureTermination::ToleranceMet);
        }
    });
}

#[test]
fn work_limits_return_bounds_and_singular_segments_are_not_skipped() {
    with_cx(|cx,_| {
        let sources=[fixture([0.0,0.0,5.0])];
        let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        let reference=HomogeneousScatterReference::new(&scene,cx).unwrap();
        let coarse=reference.enclose(cx,segment(),ScatterEnclosureBudget{max_depth:0,..budget()}).unwrap();
        assert_eq!(coarse.termination,ScatterEnclosureTermination::ResolutionLimit);
        assert!(coarse.numerical_error_bound_cd_m2>coarse.requested_error_cd_m2);
        for settings in [ScatterEnclosureBudget{max_evaluations:0,..budget()},ScatterEnclosureBudget{max_depth:49,..budget()},ScatterEnclosureBudget{absolute_error_cd_m2:f64::NAN,..budget()}] {
            assert!(matches!(reference.enclose(cx,segment(),settings),Err(ScatterEnclosureError::Budget)));
        }
        let singular=[fixture([0.0,0.0,1.0])];
        let scene=PhotometricScene::try_new(&[],&singular,HomogeneousMedium::try_new(0.0,0.1,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        assert!(matches!(HomogeneousScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()),Err(ScatterEnclosureError::UnresolvedSourceSeparation{source:0})));
    });
}

#[test]
fn exact_no_scattering_and_positive_underflow_remain_distinct() {
    with_cx(|cx,_| {
        let sources=[fixture([0.0,0.0,5.0])];
        for (absorption,scattering) in [(1.0,0.0),(200.0,0.2)] {
            let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(absorption,scattering,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
            let result=HomogeneousScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
            assert_eq!(result.numerical.lo,0.0);
            if scattering==0.0 { assert_eq!(result.numerical.hi,0.0); }
            else { assert!(result.numerical.hi>0.0 && result.numerical.hi<f64::MIN_POSITIVE); }
            assert!(result.total_error_bound_cd_m2.is_none());
        }
    });
}

#[test]
fn progressive_replay_provenance_and_failed_receiver_are_transactional() {
    with_cx(|cx,gate| {
        let sources=[fixture([0.0,0.0,5.0])];
        let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        let reference=HomogeneousScatterReference::new(&scene,cx).unwrap();
        let segments=[segment(),ScatterSegment::try_new(Point3::new(0.1,0.0,0.0),Point3::new(0.1,0.0,2.0)).unwrap()];
        let mut whole=ScatterEnclosureJob::new(&reference,&segments,budget(),cx).unwrap();
        whole.advance(cx,usize::MAX).unwrap();
        let mut chunked=ScatterEnclosureJob::new(&reference,&segments,budget(),cx).unwrap();
        chunked.advance(cx,1).unwrap();gate.request();
        assert!(matches!(chunked.advance(cx,1),Err(ScatterEnclosureError::Transport(PhotometryError::Cancelled))));
        assert_eq!(chunked.results(),&whole.results()[..1]);
        with_cx(|fresh,_| {chunked.advance(fresh,1).unwrap();});
        assert_eq!(chunked.results(),whole.results());
        assert_ne!(whole.results()[0].numerical_model_identity,whole.results()[1].numerical_model_identity);
    });
    with_cx(|cx,_| {
        let sources=[fixture([0.0,0.0,5.0])];
        let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        let reference=HomogeneousScatterReference::new(&scene,cx).unwrap();
        let segments=[segment(),ScatterSegment::try_new(Point3::new(0.0,0.0,4.0),Point3::new(0.0,0.0,6.0)).unwrap()];
        let mut job=ScatterEnclosureJob::new(&reference,&segments,budget(),cx).unwrap();
        assert!(job.advance(cx,2).is_err());assert_eq!(job.results().len(),1);
        assert!(!job.is_complete());
        let changed=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Rayleigh).unwrap(),TransportOptions::default(),cx).unwrap();
        let other=HomogeneousScatterReference::new(&changed,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
        assert_ne!(other.numerical_model_identity,job.results()[0].numerical_model_identity);
    });
}

#[test]
fn spectral_grid_is_not_silently_replaced_by_gray_haze() {
    use crate::photometry::transport::PhotometricMedium;
    use crate::volumes::VolumeGrid;
    use crate::volumes::transport::{AbsorbingGrid,AbsorbingGridBudget,MassAbsorptionSpectrum};
    with_cx(|cx,_| {
        let density=[1.0];
        let grid=AbsorbingGrid::try_new(cx,VolumeGrid::new([1,1,1],&density,[0.0;3],[1.0;3]),
            MassAbsorptionSpectrum::table(&[(380.0,0.0),(780.0,1.0)]).unwrap(),None,
            AbsorbingGridBudget{max_grid_cells:1,max_steps_per_segment:32}).unwrap();
        let sources=[fixture([0.0,0.0,5.0])];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Absorbing(&grid),TransportOptions::default(),cx).unwrap();
        assert!(matches!(HomogeneousScatterReference::new(&scene,cx),Err(ScatterEnclosureError::HomogeneousOnly)));
    });
}

#[test]
fn csv_keeps_component_bounds_separate_from_unbounded_multiple_scattering() {
    with_cx(|cx,_| {
        let sources=[fixture([0.0,0.0,5.0])];
        let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        let result=HomogeneousScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
        let mut bytes=Vec::new();EnclosedSingleScatter::write_csv_header(&mut bytes).unwrap();result.write_csv_row(0,&mut bytes).unwrap();
        let csv=String::from_utf8(bytes).unwrap();let rows:Vec<_>=csv.lines().collect();
        assert_eq!(rows[0].split(',').count(),15);let fields:Vec<_>=rows[1].split(',').collect();assert_eq!(fields.len(),15);
        assert_eq!(&fields[11..13],&["unknown","unknown"]);assert_eq!(fields[7],"tolerance-met");
        assert_eq!(fields[2].parse::<f64>().unwrap().to_bits(),result.luminance_cd_m2.to_bits());
        struct Refuse;impl Write for Refuse {
            fn write(&mut self,_:&[u8])->io::Result<usize>{Err(io::Error::other("test sink"))}
            fn flush(&mut self)->io::Result<()>{Ok(())}
        }
        assert!(result.write_csv_row(0,&mut Refuse).is_err());
    });
}

#[test]
fn a_direction_box_containing_zero_does_not_erase_a_separated_beam() {
    with_cx(|cx,_| {
        let law=AngularDistribution::cosine_power(0.5,core::f64::consts::FRAC_PI_2).unwrap();
        let sources=[PhotometricSource::try_new(
            RigidTransform::try_new([0.0,0.0,0.0,1.0],[0.0,0.5,0.0]).unwrap(),
            100.0/law.peak_per_steradian(),1.0,law,SourceSpectrum::equal_energy()).unwrap()];
        let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.1,0.2,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        let ray=ScatterSegment::try_new(Point3::new(-1.0,-1.0,-1.0),Point3::new(1.0,1.0,1.0)).unwrap();
        let result=HomogeneousScatterReference::new(&scene,cx).unwrap().enclose(cx,ray,
            ScatterEnclosureBudget{max_evaluations:1,..budget()}).unwrap();
        // Independent 80-digit integral: the real line misses the source even
        // though its axis-aligned box contains it and a normalization is unresolved.
        contains(result,2.123_493_525_137_603_6);
        assert_eq!(result.termination,ScatterEnclosureTermination::PanelBudget);
    });
}
