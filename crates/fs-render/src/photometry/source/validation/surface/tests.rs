//! Independent physical references and end-to-end use of original scene owners.
use super::*;
use crate::charts::{Ray, TriMesh};
use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use crate::photometry::{TransportOptions, medium::{HomogeneousMedium, PhaseFunction}, transport::PhotometricMedium};
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::tracer::{Primitive, Shape};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum, CellPalette, CellPaletteBudget};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, ScatteringGrid, ScatteringOptics};
use fs_evidence::NumericalKind;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Vec3;
use core::f64::consts::PI;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 173, kernel_id: 0x5355_5246, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate))
}
fn floor(rho: SpectralReflectance) -> Primitive {
    Primitive { shape: Shape::Mesh(TriMesh::new(
        vec![[-2.0,-2.0,0.0], [2.0,-2.0,0.0], [0.0,2.0,0.0]], vec![[0,1,2]])),
        material: Material::SpectralDiffuse { reflectance: rho }, emission: None }
}
fn source(position: [f64;3], cd: f64, spectrum: SourceSpectrum) -> PhotometricSource {
    PhotometricSource::try_new(RigidTransform::try_new([0.0,0.0,0.0,1.0],position).unwrap(),
        cd*4.0*PI, 1.0, AngularDistribution::isotropic(), spectrum).unwrap()
}
fn center(z: f64) -> DiffuseSurfaceSample {
    DiffuseSurfaceSample::try_new(0,[1.0,1.0,2.0],Point3::new(0.0,0.0,z)).unwrap()
}
fn check(r: EnclosedDiffuseLuminance, expected: f64) {
    assert_eq!(r.numerical.kind,NumericalKind::Enclosure);
    assert!(r.numerical.lo<=expected && expected<=r.numerical.hi, "{r:?}, expected={expected}");
    assert!((r.luminance_cd_m2-expected).abs()<=r.numerical_error_bound_cd_m2);
    assert!(r.total_error_bound_cd_m2.is_none());
}
fn screen(x0: f64,x1: f64,z: f64) -> Primitive {
    Primitive {shape:Shape::Mesh(TriMesh::new(vec![
        [x0,-1.0,z],[x1,-1.0,z],[x1,2.0,z],[x0,2.0,z]],vec![[0,1,2],[0,2,3]])),
        material:Material::SpectralDiffuse{reflectance:SpectralReflectance::constant(0.0).unwrap()},emission:None}
}

#[test]
fn lambertian_reference_is_rho_e_over_pi_not_observer_inverse_square() {
    with_cx(|cx,_| {
        let geometry=[floor(SpectralReflectance::constant(0.25).unwrap())];
        let sources=[source([0.0,0.0,3.0],100.0,SourceSpectrum::equal_energy())];
        let scene=PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let reference=DiffuseSurfaceReference::new(&scene,0,cx).unwrap();
        for height in [1.0,2.0,8.0] {
            let result=reference.enclose(cx,center(height)).unwrap();
            check(result,25.0/(9.0*PI));
            assert!(result.numerical_error_bound_cd_m2<1e-10);
            assert_eq!(result.uncertain_light_connections,0);
            assert!(!result.uncertain_observer_visibility);
        }
        // The original camera path remains untouched and supplies its nominal
        // value through the SAME material and geometry, not this reference.
        let ray=Ray{origin:Point3::new(0.0,0.0,2.0),dir:Vec3::new(0.0,0.0,-1.0)};
        let native=scene.direct_surface_luminance_along_ray(cx,&ray,2.5).unwrap();
        check(reference.enclose(cx,center(2.0)).unwrap(),native.surface_cd_m2);
        assert_eq!(native.numerical.kind,NumericalKind::NoClaim);
        // Exact relative weights, not a rounded sum-to-one admission rule.
        let extreme=DiffuseSurfaceSample::try_new(0,[f64::MAX/2.0,f64::MAX/2.0,f64::MAX],Point3::new(0.0,0.0,2.0)).unwrap();
        check(reference.enclose(cx,extreme).unwrap(),25.0/(9.0*PI));
    });
}

#[test]
fn both_homogeneous_legs_attenuate_without_shortening_receiver_distance() {
    with_cx(|cx,_| {
        let geometry=[floor(SpectralReflectance::constant(0.5).unwrap())];
        let sources=[source([0.0,0.0,3.0],100.0,SourceSpectrum::equal_energy())];
        let medium=HomogeneousMedium::try_new(0.125,0.125,PhaseFunction::HenyeyGreenstein{g:0.7}).unwrap();
        let mut last=None;
        for epsilon in [1e-6,0.01] {
            let scene=PhotometricScene::try_new(&geometry,&sources,medium.clone(),TransportOptions{ray_epsilon_m:epsilon,..TransportOptions::default()},cx).unwrap();
            let r=DiffuseSurfaceReference::new(&scene,0,cx).unwrap().enclose(cx,center(2.0)).unwrap();
            check(r,50.0/(9.0*PI)*(-1.25_f64).exp());
            assert!(r.numerical_error_bound_cd_m2<1e-10);
            if let Some((lo,hi))=last {assert_eq!((r.numerical.lo,r.numerical.hi),(lo,hi));}
            last=Some((r.numerical.lo,r.numerical.hi));
        }
    });
}

#[test]
fn observer_and_light_blockers_are_independent_and_edge_ambiguity_is_not_black() {
    with_cx(|cx,_| {
        let sources=[source([2.0,0.0,3.0],100.0,SourceSpectrum::equal_energy())];
        let expected=150.0/(13.0_f64.powf(1.5)*PI);
        for (blocker,blocked,uncertain) in [
            (screen(0.8,1.2,1.5),true,false), // light at x=1
            (screen(-0.2,0.2,1.0),true,false), // observer at x=0
            (screen(0.0,0.2,1.0),false,true), // observer exactly on edge
            (screen(1.5,1.8,1.0),false,false),
        ] {
            let geometry=[floor(SpectralReflectance::constant(0.5).unwrap()),blocker];
            let scene=PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
            let result=DiffuseSurfaceReference::new(&scene,0,cx).unwrap().enclose(cx,center(2.0)).unwrap();
            if blocked {check(result,0.0);assert_eq!(result.numerical.hi,0.0);}
            else {
                check(result,expected);
                assert_eq!(result.uncertain_observer_visibility,uncertain);
                if uncertain {assert_eq!(result.numerical.lo,0.0);}
            }
        }
    });
}

#[test]
fn indexed_spectral_haze_and_colored_material_use_both_real_column_integrals() {
    with_cx(|cx,_| {
        let density=[1.0,2.0,4.0];
        let cap=CellPaletteBudget{max_models:2,max_cells:3};
        let a=CellPalette::try_new(cx,vec![
            MassAbsorptionSpectrum::constant(0.125).unwrap(),
            MassAbsorptionSpectrum::table(&[(380.0,0.0),(780.0,0.25)]).unwrap()],Arc::from([0,1,0]),cap).unwrap();
        let absorption=AbsorbingGrid::try_new_indexed(cx,VolumeGrid::new([1,1,3],&density,[-0.5,-0.5,1.0],[1.0;3]),
            a,None,AbsorbingGridBudget{max_grid_cells:3,max_steps_per_segment:64}).unwrap();
        let s=CellPalette::try_new(cx,vec![
            ScatteringOptics{scattering:MassScatteringSpectrum::rayleigh(400.0,0.25).unwrap(),phase:PhaseFunction::Rayleigh},
            ScatteringOptics{scattering:MassScatteringSpectrum::constant(0.125).unwrap(),phase:PhaseFunction::Isotropic}],
            Arc::from([1,0,1]),cap).unwrap();
        let grid=ScatteringGrid::try_new_indexed(cx,absorption,s).unwrap();
        let geometry=[floor(SpectralReflectance::table(&[(380.0,0.125),(780.0,0.875)]).unwrap())];
        let mut bins=[0.0;SPECTRAL_BINS];bins[4]=1.0;bins[44]=2.0;
        let sources=[source([0.0,0.0,5.0],1000.0,SourceSpectrum::try_from_bins(bins).unwrap())];
        let scene=PhotometricScene::try_with_medium(&geometry,&sources,PhotometricMedium::Scattering(&grid),TransportOptions::default(),cx).unwrap();
        let reference=DiffuseSurfaceReference::new(&scene,0,cx).unwrap();
        let result=reference.enclose(cx,center(2.5)).unwrap();
        // Light column is [1,2,4]; observer column is [1,1,0]. Thus
        // tau = (2+4)*(1/8+1/8) + 3*((lambda-380)/1600+(1/4)*(400/lambda)^4).
        let mut spectral=0.0;
        for bin in [4,44] {
            let lambda=wavelength_nm(bin);let ratio=400.0/lambda;
            let tau=1.5+3.0*((lambda-380.0)/1600.0+0.25*ratio.powi(4));
            let rho=0.125+0.75*(lambda-380.0)/400.0;
            spectral+=sources[0].spectrum().bins_w_per_lumen_nm()[bin]*photopic_weight_nm(bin)*rho*(-tau).exp();
        }
        check(result,40.0/PI*683.0*spectral);
        assert!(result.numerical_error_bound_cd_m2<1e-8);
        let limited=reference.with_grid_budget(GridEnclosureBudget{max_nodes_per_segment:1}).unwrap();
        assert!(matches!(limited.enclose(cx,center(2.5)),Err(SurfaceEnclosureError::PhysicalMedium(VolumeTransportError::Budget))));
    });
}

#[test]
fn rigid_instance_and_back_face_keep_the_same_native_surface_response() {
    with_cx(|cx,_| {
        let placement=RigidTransform::try_new([0.0,1.0,0.0,0.0],[8.0,-4.0,2.0]).unwrap();
        let base=floor(SpectralReflectance::constant(0.5).unwrap());
        let Shape::Mesh(mesh)=base.shape else {unreachable!()};
        let geometry=[Primitive{shape:Shape::Instance(GeometryInstance::try_new(7,ContentHash([7;32]),
            SharedGeometry::mesh(mesh),placement).unwrap()),material:base.material,emission:None}];
        // The placed triangle faces -Z. Source/observer are on its back (+Z).
        let sources=[source([8.0,-4.0,5.0],100.0,SourceSpectrum::equal_energy()).filtered([0.25;SPECTRAL_BINS]).unwrap()];
        let scene=PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let r=DiffuseSurfaceReference::new(&scene,0,cx).unwrap().enclose(cx,
            DiffuseSurfaceSample::try_new(0,[1.0,1.0,2.0],Point3::new(8.0,-4.0,4.0)).unwrap()).unwrap();
        check(r,12.5/(9.0*PI));
        assert!(r.numerical_error_bound_cd_m2<1e-8);
    });
}

#[test]
fn cancellation_invalid_samples_and_storage_limits_preserve_completed_results() {
    with_cx(|cx,gate| {
        let geometry=[floor(SpectralReflectance::constant(0.5).unwrap())];
        let sources=[source([0.0,0.0,3.0],100.0,SourceSpectrum::equal_energy())];
        let scene=PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let reference=DiffuseSurfaceReference::new(&scene,0,cx).unwrap();
        let samples=[center(2.0),center(4.0),DiffuseSurfaceSample::try_new(0,[1.0,2.0,1.0],Point3::new(0.0,0.0,2.0)).unwrap()];
        let mut all=DiffuseSurfaceJob::new(&reference,&samples,cx).unwrap();all.advance(cx,usize::MAX).unwrap();
        let mut split=DiffuseSurfaceJob::new(&reference,&samples,cx).unwrap();split.advance(cx,1).unwrap();
        gate.request();assert!(matches!(split.advance(cx,1),Err(SurfaceEnclosureError::Transport(PhotometryError::Cancelled))));
        assert_eq!(split.results(),&all.results()[..1]);
        with_cx(|resume,_| {
            split.advance(resume,usize::MAX).unwrap();assert_eq!(split.results(),all.results());
            let invalid=DiffuseSurfaceSample::try_new(999,[1.0;3],Point3::new(0.0,0.0,2.0)).unwrap();
            let failure_samples=[center(2.0),invalid];
            let mut failure=DiffuseSurfaceJob::new(&reference,&failure_samples,resume).unwrap();
            assert!(matches!(failure.advance(resume,2),Err(SurfaceEnclosureError::InvalidSample)));
            assert_eq!(failure.results().len(),1);
            assert!(!failure.is_complete());
            let coplanar=DiffuseSurfaceSample::try_new(0,[1.0,1.0,2.0],Point3::new(1.0,0.0,0.0)).unwrap();
            assert!(matches!(reference.enclose(resume,coplanar),Err(SurfaceEnclosureError::UnresolvedObserverSide)));
            let short=PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions{max_outputs:1,..TransportOptions::default()},resume).unwrap();
            let short=DiffuseSurfaceReference::new(&short,0,resume).unwrap();
            assert!(matches!(DiffuseSurfaceJob::new(&short,&samples,resume),Err(SurfaceEnclosureError::Budget)));
        });
    });
    for weights in [[0.0;3],[-1.0,1.0,1.0],[f64::NAN,1.0,1.0],[f64::INFINITY,1.0,1.0]] {
        assert!(matches!(DiffuseSurfaceSample::try_new(0,weights,Point3::new(0.0,0.0,2.0)),Err(SurfaceEnclosureError::InvalidSample)));
    }
}

#[test]
fn exact_black_and_observer_variants_keep_model_identity_and_shared_preparation() {
    with_cx(|cx,_| {
        // Coincident declarations share the same opaque boundary, not extra light.
        let geometry=[floor(SpectralReflectance::constant(0.0).unwrap()),floor(SpectralReflectance::constant(0.5).unwrap())];
        let sources=[source([0.0,0.0,3.0],100.0,SourceSpectrum::equal_energy())];
        let scene=PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let black=DiffuseSurfaceReference::new(&scene,0,cx).unwrap();
        let gray=black.for_primitive(1,cx).unwrap();
        assert!(Arc::ptr_eq(&black.prepared,&gray.prepared));
        let a=black.enclose(cx,center(2.0)).unwrap();let b=gray.enclose(cx,center(2.0)).unwrap();
        check(a,0.0);assert_eq!(a.numerical.hi,0.0);assert_eq!(a.numerical_error_bound_cd_m2,0.0);
        check(b,50.0/(9.0*PI));
        assert_ne!(a.numerical_model_identity,b.numerical_model_identity);
        let c=gray.enclose(cx,center(4.0)).unwrap();
        assert_ne!(b.numerical_model_identity,c.numerical_model_identity);
        assert_eq!(b,gray.enclose(cx,center(2.0)).unwrap());
        assert!(matches!(black.for_primitive(2,cx),Err(SurfaceEnclosureError::InvalidSample)));
    });
}

#[test]
fn retained_reflectance_interpolant_and_positive_underflow_remain_enclosed() {
    let rho=SpectralReflectance::table(&[(380.0,0.0),(580.0,1.0),(780.0,0.0)]).unwrap();
    for (lambda,expected) in [(380.0,0.0),(430.0,0.25),(580.0,1.0),(680.0,0.5),(780.0,0.0)] {
        let r=reflectance_enclosure(&rho,lambda);assert!(r.lo()<=expected&&r.hi()>=expected);
    }
    with_cx(|cx,_| {
        let geometry=[floor(SpectralReflectance::constant(0.5).unwrap())];
        let sources=[source([0.0,0.0,3.0],100.0,SourceSpectrum::equal_energy())];
        let medium=HomogeneousMedium::try_new(200.0,0.0,PhaseFunction::Isotropic).unwrap();
        let scene=PhotometricScene::try_new(&geometry,&sources,medium,TransportOptions::default(),cx).unwrap();
        let r=DiffuseSurfaceReference::new(&scene,0,cx).unwrap().enclose(cx,center(2.0)).unwrap();
        assert!(r.numerical.hi>0.0&&r.numerical.hi<f64::MIN_POSITIVE);
        assert!(r.numerical_error_bound_cd_m2>0.0);
        assert!(r.total_error_bound_cd_m2.is_none());
    });
}

#[test]
fn csv_retains_component_scope_units_unknowns_and_writer_errors() {
    with_cx(|cx,_| {
        let geometry=[floor(SpectralReflectance::constant(0.5).unwrap())];
        let sources=[source([0.0,0.0,3.0],100.0,SourceSpectrum::equal_energy())];
        let scene=PhotometricScene::try_new(&geometry,&sources,HomogeneousMedium::vacuum(),TransportOptions::default(),cx).unwrap();
        let r=DiffuseSurfaceReference::new(&scene,0,cx).unwrap().enclose(cx,center(2.0)).unwrap();
        let mut bytes=Vec::new();EnclosedDiffuseLuminance::write_csv_header(&mut bytes).unwrap();r.write_csv_row(3,&mut bytes).unwrap();
        let csv=String::from_utf8(bytes).unwrap();let rows:Vec<_>=csv.lines().collect();
        let fields:Vec<_>=rows[1].split(',').collect();assert_eq!(fields.len(),13);assert_eq!(rows[0].split(',').count(),13);
        assert_eq!(fields[1],"direct-diffuse-surface");assert_eq!(fields[11],"unknown");
        assert_eq!(fields[4].parse::<f64>().unwrap().to_bits(),r.luminance_cd_m2.to_bits());
        struct Refuse;impl Write for Refuse {
            fn write(&mut self,_:&[u8])->io::Result<usize>{Err(io::Error::other("test sink"))}
            fn flush(&mut self)->io::Result<()>{Ok(())}
        }
        assert!(r.write_csv_row(0,&mut Refuse).is_err());
    });
}
