use super::*;
use super::super::{Interaction, PathTermination, SurfaceBinding};
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::{Stokes, dielectric::transport::PolarizationFrame};
use crate::volumes::{VolumeGrid, transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum}};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction, ScatteringGrid};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| f(&Cx::new(
        &gate, arena, StreamKey { seed:31,kernel_id:7,tile:0,iteration:0 },
        Budget::INFINITE, ExecMode::Deterministic,
    )))
}
fn glass(n:f64,a:f64) -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(a).unwrap(), GlassProvenance::Custom)
}
fn launch(material:DielectricGlass,wavelength:f64,power:Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(Point3::new(0.0,0.0,0.0),
        PolarizationFrame::try_new(Vec3::new(0.0,0.0,1.0),Vec3::new(1.0,0.0,0.0)).unwrap(),
        wavelength,material,power).unwrap()
}
fn plane(id:u64,axis:usize,position:f64) -> GeometryInstance {
    let vertices = if axis==0 {
        vec![[position,-10.0,-10.0],[position,10.0,-10.0],[position,0.0,10.0]]
    } else {
        vec![[-10.0,-10.0,position],[10.0,-10.0,position],[0.0,10.0,position]]
    };
    GeometryInstance::try_new(id,fs_blake3::hash_domain("org.frankensim.test.rayleigh-receiver",&id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vertices,vec![[0,1,2]])),RigidTransform::identity()).unwrap()
}
fn detector(id:u64) -> SurfaceBinding {
    SurfaceBinding { object_id:id,interaction:Interaction::Detector { reference_axis:Vec3::new(0.0,1.0,0.0) } }
}
fn settings() -> RayleighTraceSettings {
    RayleighTraceSettings { surfaces:TraceSettings {
        max_surface_events:8,max_distance_m:4.0,surface_epsilon_m:1e-9,query_offset_m:1e-8,
    },max_volume_events:8 }
}
fn grid<'a>(cx:&Cx<'_>,density:&'a [f64],absorption:f64,scattering:f64) -> ScatteringGrid<'a> {
    let absorption=AbsorbingGrid::try_new(cx,
        VolumeGrid::new([1,1,1],density,[-2.0;3],[4.0;3]),
        MassAbsorptionSpectrum::constant(absorption).unwrap(),None,
        AbsorbingGridBudget { max_grid_cells:1,max_steps_per_segment:16 }).unwrap();
    ScatteringGrid::try_new(cx,absorption,MassScatteringSpectrum::constant(scattering).unwrap(),PhaseFunction::Rayleigh).unwrap()
}
fn near(a:f64,b:f64) { assert!((a-b).abs()<3e-11*a.abs().max(b.abs()).max(1.0),"{a:.17e} != {b:.17e}"); }

#[test]
fn volume_collision_preempts_geometry_and_reaches_original_receiver() {
    with_cx(|cx| {
        let density=[1.0];let medium=grid(cx,&density,0.0,1.0);let media=[glass(1.0,0.0)];
        let geometry=InstanceScene::try_new(vec![plane(7,0,1.0)]).unwrap();
        let scene=PolarizedScene::try_new(&geometry,&media,vec![detector(7)]).unwrap();
        let bindings=[Some(RayleighGrid::try_new(cx,&medium).unwrap())];
        let participating=RayleighScene::try_new(cx,&scene,&bindings).unwrap();
        let input=launch(media[0],550.0,Stokes::unpolarized(1.0));
        let mut calls=0;
        let out=participating.trace_with_samples(cx,input,0,settings(),|_| panic!("receiver needs no draw"),|flight| {
            assert_eq!(flight,calls); calls+=1;
            if flight==0 { [(-0.5_f64).exp(),0.5,0.0] } else { [0.0,0.5,0.0] }
        }).unwrap();
        assert_eq!(calls,2);assert_eq!(out.free_flight_draws(),2);assert_eq!(out.volume_events(),1);
        assert_eq!(out.path().termination(),PathTermination::Detected);
        assert_eq!(out.path().surface_events(),1);assert_eq!(out.path().last_hit().unwrap().object_id,7);
        near(out.path().traveled_m(),1.5);
        near(out.path().beam().ray().origin.x,1.0);near(out.path().beam().ray().origin.z,0.5);
        near(out.path().beam().power().i,1.0);near(out.path().beam().power().q,1.0);
        near(out.path().beam().power().u,0.0);near(out.path().beam().power().v,0.0);
        assert_eq!(participating.trace(cx,input,0,settings(),31,42).unwrap(),
            participating.trace(cx,input,0,settings(),31,42).unwrap());
    });
}

#[test]
fn successive_scattering_planes_retain_dipole_extinction() {
    with_cx(|cx| {
        let density=[1.0];let medium=grid(cx,&density,0.0,1.0);let media=[glass(1.0,0.0)];
        let geometry=InstanceScene::try_new(vec![plane(7,0,1.0)]).unwrap();
        let scene=PolarizedScene::try_new(&geometry,&media,vec![detector(7)]).unwrap();
        let bindings=[Some(RayleighGrid::try_new(cx,&medium).unwrap())];
        let participating=RayleighScene::try_new(cx,&scene,&bindings).unwrap();
        let out=participating.trace_with_samples(cx,launch(media[0],550.0,Stokes::unpolarized(1.0)),0,settings(),
            |_| panic!("no surface before extinction"),|flight| match flight {
                0=>[(-0.5_f64).exp(),0.5,0.0],
                1=>[(-0.25_f64).exp(),0.5,0.25],
                _=>panic!("a y-polarized dipole cannot scatter along y"),
            }).unwrap();
        assert_eq!(out.path().termination(),PathTermination::Extinguished);
        assert_eq!(out.volume_events(),2);assert_eq!(out.path().surface_events(),0);
        near(out.path().traveled_m(),0.75);near(out.path().beam().power().i,0.0);
    });
}

#[test]
fn explicit_medium_binding_excludes_glass_and_preserves_surface_only_paths() {
    with_cx(|cx| {
        let density=[1.0];let medium=grid(cx,&density,0.1,0.0);
        let media=[glass(1.0,0.0),glass(1.5,0.2)];
        let geometry=InstanceScene::try_new(vec![plane(1,2,0.5),plane(2,2,1.5),plane(7,2,2.0)]).unwrap();
        let scene=PolarizedScene::try_new(&geometry,&media,vec![
            SurfaceBinding { object_id:1, interaction:Interaction::Dielectric { inside:0,outside:1,surface:DielectricSurface::SMOOTH } },
            SurfaceBinding { object_id:2, interaction:Interaction::Dielectric { inside:1,outside:0,surface:DielectricSurface::SMOOTH } },
            detector(7),
        ]).unwrap();
        let input=launch(media[0],550.0,Stokes::try_new(1.0,0.0,0.0,1.0).unwrap());
        let empty=[None,None];
        let unbound=RayleighScene::try_new(cx,&scene,&empty).unwrap();
        for seed in [1,31,97] {
            assert_eq!(*unbound.trace(cx,input,0,settings(),seed,42).unwrap().path(),
                scene.trace(cx,input,0,settings().surfaces,seed,42).unwrap());
        }
        let bindings=[Some(RayleighGrid::try_new(cx,&medium).unwrap()),None];
        let participating=RayleighScene::try_new(cx,&scene,&bindings).unwrap();
        let mut surfaces=0;
        let out=participating.trace_with_samples(cx,input,0,settings(),|event| {
            assert_eq!(event,surfaces);surfaces+=1;[0.3,0.7,0.999999]
        },|_| panic!("pure absorption consumes no random draw")).unwrap();
        assert_eq!(surfaces,2);assert_eq!(out.volume_events(),0);assert_eq!(out.free_flight_draws(),0);
        assert_eq!(out.path().termination(),PathTermination::Detected);assert_eq!(out.path().medium_index(),0);
        near(out.path().traveled_m(),2.0);
        // Conditional transmission weights cancel the scalar branch probability.
        // Only 1 m of grid absorption and 1 m of native glass absorption remain.
        near(out.path().beam().power().i,(-0.3_f64).exp());
        near(out.path().beam().power().v,(-0.3_f64).exp());
    });
}

#[test]
fn absorbing_cell_in_molecular_palette_terminates_without_invented_phase() {
    use crate::volumes::transport::{CellPalette,CellPaletteBudget};
    use crate::volumes::transport::scattering::ScatteringOptics;
    with_cx(|cx| {
        let density=[1.0,1.0];
        let absorption=AbsorbingGrid::try_new(cx,
            VolumeGrid::new([1,1,2],&density,[-1.0,-1.0,0.0],[2.0,2.0,1.0]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),None,
            AbsorbingGridBudget { max_grid_cells:2,max_steps_per_segment:16 }).unwrap();
        let palette=CellPalette::try_new(cx,vec![
            ScatteringOptics { scattering:MassScatteringSpectrum::constant(0.0).unwrap(),phase:PhaseFunction::Isotropic },
            ScatteringOptics { scattering:MassScatteringSpectrum::constant(1.0).unwrap(),phase:PhaseFunction::Rayleigh },
        ],vec![0_u32,1].into(),CellPaletteBudget { max_models:2,max_cells:2 }).unwrap();
        let medium=ScatteringGrid::try_new_indexed(cx,absorption,palette).unwrap();
        let rayleigh=RayleighGrid::try_new(cx,&medium).unwrap();
        let input=launch(glass(1.0,0.0),550.0,Stokes::unpolarized(1.0));
        let out=rayleigh.sample_beam_segment(cx,input,3.0,[(-0.25_f64).exp(),0.5,0.0]).unwrap();
        assert_eq!(out.event.unwrap().cell_index,0);assert_eq!(out.phase_pdf,None);
        assert_eq!(out.beam.power(),Stokes::unpolarized(0.0));
    });
}

#[test]
fn incompatible_carrier_refuses_before_consuming_flight_stream() {
    with_cx(|cx| {
        let density=[1.0];let medium=grid(cx,&density,0.0,1.0);
        let geometry=InstanceScene::try_new(vec![plane(7,2,1.0)]).unwrap();
        for material in [glass(1.33,0.0),glass(1.0,0.2)] {
            let media=[material];
            let scene=PolarizedScene::try_new(&geometry,&media,vec![detector(7)]).unwrap();
            let bindings=[Some(RayleighGrid::try_new(cx,&medium).unwrap())];
            let participating=RayleighScene::try_new(cx,&scene,&bindings).unwrap();
            assert!(matches!(RayleighScene::try_new(cx,&scene,&[]),Err(SceneError::InvalidBindings)));
            let input=launch(material,550.0,Stokes::unpolarized(1.0));
            assert_eq!(participating.trace_with_samples(cx,input,0,settings(),
                |_| panic!("no surface reached"),|_| panic!("invalid carrier consumes no draw")),
                Err(SceneError::Scattering(PolarizedScatteringError::CarrierMedium)));
        }
    });
}

#[cfg(feature="photometry")]
#[test]
fn calibrated_emission_scattering_and_meter_keep_power_and_truncation() {
    use crate::photometry::source::{AngularDistribution,PhotometricSource,SourceSpectrum,SPECTRAL_BINS};
    use super::super::measurement::MeterError;
    with_cx(|cx| {
        let density=[1.0];let medium=grid(cx,&density,0.0,1.0);let media=[glass(1.0,0.0)];
        let geometry=InstanceScene::try_new(vec![plane(7,0,1.0)]).unwrap();
        let scene=PolarizedScene::try_new(&geometry,&media,vec![detector(7)]).unwrap();
        let bindings=[Some(RayleighGrid::try_new(cx,&medium).unwrap())];
        let participating=RayleighScene::try_new(cx,&scene,&bindings).unwrap();
        let source=PhotometricSource::try_new(RigidTransform::identity(),1000.0,1.0,
            AngularDistribution::isotropic(),SourceSpectrum::try_from_bins([1.0;SPECTRAL_BINS]).unwrap()).unwrap()
            .with_unpolarized_optics();
        let mut meter=scene.power_meter(7).unwrap();let mut limited=scene.power_meter(7).unwrap();
        for bin in 0..SPECTRAL_BINS {
            let beam=source.emit_beam(Vec3::new(0.0,0.0,1.0),bin,4.0*core::f64::consts::PI,media[0]).unwrap();
            let draws=|flight| if flight==0 { [(-0.5_f64).exp(),0.5,0.0] } else { [0.0,0.5,0.0] };
            let out=participating.trace_with_samples(cx,beam,0,settings(),|_| [0.5;3],draws).unwrap();
            meter.record(bin,out.path()).unwrap();
            let mut cap=settings();cap.max_volume_events=1;
            let truncated=participating.trace_with_samples(cx,beam,0,cap,|_| [0.5;3],draws).unwrap();
            assert_eq!(truncated.path().termination(),PathTermination::VolumeLimit);
            assert!(truncated.path().beam().power().i>0.0);
            limited.record(bin,truncated.path()).unwrap();
        }
        // This is one explicitly weighted path per bin, not angular convergence
        // or a claim that all fixture light physically reaches this receiver.
        let reading=meter.reading().unwrap().require_complete().unwrap();
        near(reading.photopic_stokes_lm().i,1000.0);near(reading.photopic_stokes_lm().q,1000.0);
        near(reading.linear_analyzer(0.0).unwrap().luminous_lm,1000.0);
        near(reading.linear_analyzer(core::f64::consts::FRAC_PI_2).unwrap().luminous_lm,0.0);
        let incomplete=limited.reading().unwrap();
        assert_eq!(incomplete.launched_paths(),SPECTRAL_BINS as u64);
        assert_eq!(incomplete.detected_paths(),0);
        assert_eq!(incomplete.require_complete(),Err(MeterError::UnresolvedPaths { count:SPECTRAL_BINS as u64 }));
    });
}
