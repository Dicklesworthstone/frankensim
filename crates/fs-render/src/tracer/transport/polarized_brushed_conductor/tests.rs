use super::*;
use crate::conductor::{CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorSource};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::surface_detail::{BrushingPattern, RadialMark};
use crate::tracer::{Material, opaque_bsdf_eval};

fn medium() -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(1.33).unwrap(), BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
}
fn optics() -> ConductorOptics {
    let source = ConductorSource::try_new(fs_blake3::ContentHash([92;32]), ConductorDataStatus::Representative).unwrap();
    ConductorOptics::try_new(core::array::from_fn(|i|
        ConductorIorSample::try_new(380.0+50.0*i as f64, 0.2, 3.0).unwrap()), source).unwrap()
        .with_coating(CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.45).unwrap(),125e-9,source).unwrap())
}
fn beam(wo: Vec3, state: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(Point3::new(0.4,0.0,0.0),
        PolarizationFrame::try_new(wo.scale(-1.0),Vec3::new(1.0,0.0,0.0)).unwrap(),550.0,medium(),state).unwrap()
}
fn detail(alpha: f64, ratio: f64, tilt: f64) -> SurfaceDetail {
    SurfaceDetail::try_new(alpha,ratio,BrushingPattern::Radial,tilt,None).unwrap()
}
fn frame(detail: &SurfaceDetail) -> ShadingFrame {
    shading_frame(detail,Point3::new(0.4,0.0,0.0),Vec3::new(0.0,0.0,1.0)).unwrap()
}
fn near(a: f64,b: f64) { assert!((a-b).abs() <= 3e-9*a.abs().max(b.abs()).max(1.0), "{a:.17e} != {b:.17e}"); }
fn same(a: Stokes,b: Stokes) {
    for (a,b) in [a.i,a.q,a.u,a.v].into_iter().zip([b.i,b.q,b.u,b.v]) { near(a,b); }
}

#[test]
fn g1_brushed_forward_power_matches_reversed_native_radiance_bsdf() {
    let law=optics(); let mut accepted=0; let mut rejected=0; let mut above_one=false;
    for alpha in [0.05,0.3] { for ratio in [1.0,3.0] { for tilt in [0.0,0.2,0.35] {
        let detail=detail(alpha,ratio,tilt); let frame=frame(&detail);
        let material=Material::BrushedConductor { optics:law,detail };
        for cosine in [0.005_f64,0.4,1.0] {
            let source=beam(Vec3::new(0.0,(1.0-cosine*cosine).sqrt(),cosine),Stokes::unpolarized(1.0));
            let wo=source.ray().dir.scale(-1.0);
            for draw in [[0.0,0.0],[0.13,0.27],[0.6,0.9],[0.99,0.4]] {
                let native=sample_ggx_vndf_reflection_aniso(brushed_local(&frame,wo),alpha,alpha*ratio,draw[0],draw[1])
                    .map(|(local,pdf)| (brushed_world(&frame,local),pdf))
                    .filter(|(wi,_)| frame.geometric_normal.dot(*wi)>0.0);
                let actual=source.sample_brushed_conductor(&law,&detail,&frame,draw).unwrap();
                match (native,actual) {
                    (None,None)=>rejected+=1,
                    (Some((wi,pdf)),Some(out))=> {
                        accepted+=1; assert_eq!(out.pdf,pdf); assert!(!out.delta);
                        // Reverse the effective radiance law, not the light ray.
                        let expected=opaque_bsdf_eval(&material,frame.geometric_normal,wi,wo,550.0,
                            Some(medium()),Some(&frame)).unwrap()*frame.geometric_normal.dot(wi)/pdf;
                        near(out.beam.power().i,expected);
                        above_one |= out.beam.power().i>1.0;
                        assert_eq!(out.beam.medium(),source.medium());
                        assert_eq!(out.event,DielectricEvent::Reflection);
                        for (a,b) in [wi.x,wi.y,wi.z].into_iter().zip([
                            out.beam.ray().dir.x,out.beam.ray().dir.y,out.beam.ray().dir.z]) { near(a,b); }
                    }
                    other=>panic!("native support differs: {other:?}"),
                }
            }
        }
    } } }
    assert!(accepted>0 && rejected>0 && above_one);
}

#[test]
fn g3_polarization_reference_changes_do_not_rotate_the_brushing() {
    let d=detail(0.1,4.0,0.2); let f=frame(&d); let law=optics();
    for state in [Stokes::try_new(2.0,0.4,-0.7,0.3).unwrap(),Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()] {
        let source=beam(Vec3::new(0.0,0.6,0.8),state);
        let changed=source.with_reference_axis(Vec3::new(1.0,2.0,3.0)).unwrap();
        let a=source.sample_brushed_conductor(&law,&d,&f,[0.13,0.27]).unwrap().unwrap();
        let b=changed.sample_brushed_conductor(&law,&d,&f,[0.13,0.27]).unwrap().unwrap();
        same(a.beam.power(),b.beam.power()); assert_eq!(a.pdf,b.pdf); assert_eq!(a.beam.ray(),b.beam.ray());
        if source.power().degree() == 1.0 { near(a.beam.power().degree(),1.0); }
    }
}

#[test]
fn g3_instanced_brushing_moves_with_geometry_and_not_world_axes() {
    let d=detail(0.1,4.0,0.15); let law=optics();
    let source=beam(Vec3::new(0.0,0.6,0.8),Stokes::try_new(1.0,0.0,0.0,1.0).unwrap());
    let n=Vec3::new(0.0,0.0,1.0); let local=source.ray().origin;
    let identity=RigidTransform::identity();
    let angle=0.7_f64;
    let placed=RigidTransform::try_new([0.0,(angle/2.0).sin(),0.0,(angle/2.0).cos()],[2.0,-3.0,4.0]).unwrap();
    let moved=PolarizedBeam::try_new(placed.transform_point(local),
        PolarizationFrame::try_new(placed.transform_vector(source.ray().dir),placed.transform_vector(source.frame().horizontal())).unwrap(),
        550.0,medium(),source.power()).unwrap();
    let a=source.sample_instanced_brushed_conductor(&law,&d,identity,local,n,[0.13,0.27]).unwrap().unwrap();
    let b=moved.sample_instanced_brushed_conductor(&law,&d,placed,local,placed.transform_vector(n),[0.13,0.27]).unwrap().unwrap();
    same(a.beam.power(),b.beam.power()); near(a.pdf,b.pdf);
    let expected=placed.transform_vector(a.beam.ray().dir);
    near(expected.x,b.beam.ray().dir.x); near(expected.y,b.beam.ray().dir.y); near(expected.z,b.beam.ray().dir.z);
    assert_eq!(b.beam.ray().origin,moved.ray().origin);
}

#[test]
fn g3_zero_detail_and_undefined_radial_poles_reuse_isotropic_samples_exactly() {
    let law=optics(); let n=Vec3::new(0.0,0.0,1.0);
    let source=beam(Vec3::new(0.0,0.6,0.8),Stokes::unpolarized(1.0));
    for raw_normal in [n,Vec3::new(0.2,0.3,1.0),Vec3::new(0.02,-0.03,-1.0)] {
        assert_eq!(source.sample_instanced_brushed_conductor(&law,&SurfaceDetail::none(0.2),
            RigidTransform::identity(),source.ray().origin,raw_normal,[0.13,0.27]).unwrap(),
            source.sample_rough_conductor(raw_normal,&law,ConductorSurface::try_rough(0.2).unwrap(),[0.13,0.27]).unwrap());
    }
    let expected=source.sample_rough_conductor(n,&law,ConductorSurface::try_rough(0.2).unwrap(),[0.13,0.27]).unwrap();
    for (d,point) in [(SurfaceDetail::none(0.2),Point3::new(0.4,0.2,0.0)),
        (detail(0.2,3.0,0.2),Point3::new(0.0,0.0,2.0))] {
        assert_eq!(source.sample_instanced_brushed_conductor(&law,&d,RigidTransform::identity(),point,n,[0.13,0.27]).unwrap(),expected);
    }
}

#[test]
fn g1_radial_mark_changes_the_native_anisotropic_proposal() {
    let law=optics();
    let d=SurfaceDetail::try_new(0.12,3.0,BrushingPattern::Radial,0.1,
        Some(RadialMark { azimuth_rad:0.0,half_width_rad:0.1 })).unwrap();
    let marked=frame(&d); assert!(marked.mark_covered);
    let mut plain=marked; plain.mark_covered=false;
    let source=beam(Vec3::new(0.0,0.6,0.8),Stokes::unpolarized(1.0));
    let draw=[0.13,0.27];
    let a=source.sample_brushed_conductor(&law,&d,&marked,draw).unwrap().unwrap();
    let b=source.sample_brushed_conductor(&law,&d,&plain,draw).unwrap().unwrap();
    assert_ne!(a.pdf,b.pdf); assert_ne!(a.beam.ray().dir,b.beam.ray().dir);
    let (ax,ay)=d.alpha_pair_for_mark(true);
    let native=sample_ggx_vndf_reflection_aniso(brushed_local(&marked,source.ray().dir.scale(-1.0)),ax,ay,draw[0],draw[1]).unwrap();
    assert_eq!(a.pdf,native.1);
}

#[test]
fn g0_invalid_frames_and_uniforms_refuse_without_mutating_the_beam() {
    let source=beam(Vec3::new(0.0,0.6,0.8),Stokes::unpolarized(1.0)); let original=source;
    let law=optics(); let d=detail(0.1,3.0,0.2); let f=frame(&d);
    let mut broken=[f;4];
    broken[0].tangent=broken[0].bitangent;
    broken[1].bitangent=broken[1].bitangent.scale(-1.0);
    broken[2].shading_normal.x=f64::NAN;
    broken[3].geometric_normal=Vec3::new(0.0,0.0,2.0);
    for frame in broken { assert!(source.sample_brushed_conductor(&law,&d,&frame,[0.1,0.2]).is_err()); }
    for draw in [[-0.1,0.0],[0.0,1.0],[f64::NAN,0.0]] {
        assert_eq!(source.sample_brushed_conductor(&law,&d,&f,draw),Err(ConductorError::InvalidSample.into()));
    }
    assert!(source.sample_brushed_conductor(&law,&SurfaceDetail::none(f64::NAN),&f,[0.1,0.2]).is_err());
    assert_eq!(source,original);
}

#[test]
fn g0_shading_horizon_is_a_zero_sample_not_a_reoriented_shading_frame() {
    let source=beam(Vec3::new(0.0,-(1.0-0.01_f64*0.01).sqrt(),0.01),Stokes::unpolarized(1.0));
    let d=detail(0.1,3.0,0.2); let f=frame(&d);
    assert!(f.geometric_normal.dot(source.ray().dir.scale(-1.0))>0.0);
    assert!(f.shading_normal.dot(source.ray().dir.scale(-1.0))<0.0);
    assert_eq!(source.sample_brushed_conductor(&optics(),&d,&f,[0.1,0.2]).unwrap(),None);
}
