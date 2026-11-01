//! Physical cell integrals, independent spectra and complete-ray regressions.
use super::*;
use crate::charts::TriMesh;
use crate::instances::RigidTransform;
use crate::photometry::{TransportOptions, medium::HomogeneousMedium, transport::PhotometricMedium};
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum,
    SPECTRAL_BINS, LUMINOUS_EFFICACY_LM_W, photopic_weight_nm, wavelength_nm};
use crate::spectral::SpectralReflectance;
use crate::tracer::{Material, Primitive, Shape};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, CellPalette,
    CellPaletteBudget, MassAbsorptionSpectrum, VolumeTransportError};
use crate::volumes::transport::enclosure::GridEnclosureBudget;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, ScatteringGrid, ScatteringOptics};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::sync::Arc;
use core::f64::consts::PI;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena|f(&Cx::new(&gate,arena,
        StreamKey{seed:173,kernel_id:0x4753_4341,tile:0,iteration:0},
        Budget::INFINITE,ExecMode::Deterministic),&gate))
}
fn fixture(blue:f64, x:f64)->PhotometricSource {
    let mut bins=[0.0;SPECTRAL_BINS];
    bins[4]=blue/(LUMINOUS_EFFICACY_LM_W*photopic_weight_nm(4));
    bins[44]=(1.0-blue)/(LUMINOUS_EFFICACY_LM_W*photopic_weight_nm(44));
    PhotometricSource::try_new(RigidTransform::try_new([0.0,0.0,0.0,1.0],[x,0.0,4.0]).unwrap(),
        400.0*PI,1.0,AngularDistribution::isotropic(),SourceSpectrum::try_from_bins(bins).unwrap()).unwrap()
}
fn segment()->ScatterSegment {
    ScatterSegment::try_new(Point3::new(0.0,0.0,-1.0),Point3::new(0.0,0.0,3.5)).unwrap()
}
fn budget()->ScatterEnclosureBudget {
    ScatterEnclosureBudget{absolute_error_cd_m2:1e-3,max_evaluations:16_384,max_depth:40}
}
fn field<'a>(cx:&Cx<'_>,rho:&'a [f64])->ScatteringGrid<'a> {
    let cap=CellPaletteBudget{max_models:2,max_cells:3};
    let a=CellPalette::try_new(cx,vec![MassAbsorptionSpectrum::constant(0.05).unwrap(),
        MassAbsorptionSpectrum::table(&[(380.0,0.01),(780.0,0.2)]).unwrap()],Arc::from([0,1,0]),cap).unwrap();
    let a=AbsorbingGrid::try_new_indexed(cx,
        VolumeGrid::new([1,1,3],rho,[-2.0,-0.5,0.0],[4.0,1.0,1.0]),a,None,
        AbsorbingGridBudget{max_grid_cells:3,max_steps_per_segment:64}).unwrap();
    let s=CellPalette::try_new(cx,vec![
        ScatteringOptics{scattering:MassScatteringSpectrum::constant(0.1).unwrap(),phase:PhaseFunction::Isotropic},
        ScatteringOptics{scattering:MassScatteringSpectrum::rayleigh(400.0,0.2).unwrap(),
            phase:PhaseFunction::spectral_table(&[(380.0,PhaseFunction::Isotropic),(780.0,PhaseFunction::Rayleigh)]).unwrap()},
    ],Arc::from([0,1,0]),cap).unwrap();
    ScatteringGrid::try_new_indexed(cx,a,s).unwrap()
}
fn contains(result:EnclosedSingleScatter, truth:f64) {
    assert!(result.numerical.lo<=truth && truth<=result.numerical.hi,"{result:?} excludes {truth}");
    assert!((result.luminance_cd_m2-truth).abs()<=result.numerical_error_bound_cd_m2);
    assert!(result.total_error_bound_cd_m2.is_none());
    assert!(result.omitted_multiple_scattering_bound_cd_m2.is_none());
}
fn collinear_reference(source:&PhotometricSource,rho:[f64;3])->f64 {
    // Physical z antiderivative, independent of grid traversal/phase helpers.
    // Both optical legs together cross every layer exactly once.
    let mut total=0.0;
    for bin in [4,44] {
        let lambda=wavelength_nm(bin);let r=400.0/lambda;let t=(lambda-380.0)/400.0;
        let ks=[0.1,0.2*r*r*r*r,0.1];let ka=[0.05,0.01+0.19*t,0.05];
        let phase=[1.0,1.0+0.5*t,1.0];
        let tau=(0..3).map(|i|rho[i]*(ka[i]+ks[i])).sum::<f64>();
        let integral=(0..3).map(|i|rho[i]*ks[i]*phase[i]
            *(1.0/(3.0-i as f64)-1.0/(4.0-i as f64))).sum::<f64>()/(4.0*PI);
        let weight=source.spectrum().bins_w_per_lumen_nm()[bin]*photopic_weight_nm(bin)*683.0;
        total+=weight*100.0*(-tau).exp()*integral;
    }
    total
}
fn opaque(vertices:Vec<[f64;3]>, triangles:Vec<[u32;3]>)->Primitive {
    Primitive{shape:Shape::Mesh(TriMesh::new(vertices,triangles)),
        material:Material::SpectralDiffuse{reflectance:SpectralReflectance::constant(0.0).unwrap()},emission:None}
}

#[test]
fn layered_spectral_glow_encloses_the_independent_two_leg_antiderivative() {
    with_cx(|cx,_| {
        let density=[1.0,2.0,0.5];let grid=field(cx,&density);let sources=[fixture(0.75,0.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&grid),TransportOptions::default(),cx).unwrap();
        let reference=SingleScatterReference::new(&scene,cx).unwrap();
        let truth=collinear_reference(&sources[0],density);
        let coarse=reference.enclose(cx,segment(),ScatterEnclosureBudget{max_evaluations:1,..budget()}).unwrap();
        contains(coarse,truth);assert_eq!(coarse.termination,ScatterEnclosureTermination::PanelBudget);
        let fine=reference.enclose(cx,segment(),budget()).unwrap();
        contains(fine,truth);assert_eq!(fine.termination,ScatterEnclosureTermination::ToleranceMet);
        assert!(fine.numerical_error_bound_cd_m2<=1e-3);
        assert_eq!(fine.medium_identity,Some(grid.content_identity()));
        assert_eq!(coarse.numerical_model_identity,fine.numerical_model_identity);
    });
}

#[test]
fn blue_and_red_sources_retain_different_extinction_and_phase_weights() {
    with_cx(|cx,_| {
        let density=[1.0,2.0,0.5];let grid=field(cx,&density);let mut results=Vec::new();
        for blue in [1.0,0.0] {
            let sources=[fixture(blue,0.0)];
            let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&grid),TransportOptions::default(),cx).unwrap();
            let value=SingleScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
            contains(value,collinear_reference(&sources[0],density));results.push(value);
        }
        assert!(results[0].numerical.lo>results[1].numerical.hi);
        assert_ne!(results[0].numerical_model_identity,results[1].numerical_model_identity);
        assert_eq!(results[0].medium_identity,results[1].medium_identity);
    });
}

#[test]
fn a_thin_scattering_layer_cannot_disappear_between_camera_samples() {
    with_cx(|cx,_| {
        let density=[0.0,0.0,8.0,0.0,0.0,0.0,0.0,0.0];
        let a=AbsorbingGrid::try_new(cx,VolumeGrid::new([1,1,8],&density,[-0.5,-0.5,0.25],[1.0,1.0,0.125]),
            MassAbsorptionSpectrum::constant(0.0).unwrap(),None,AbsorbingGridBudget{max_grid_cells:8,max_steps_per_segment:64}).unwrap();
        let grid=ScatteringGrid::try_new(cx,a,MassScatteringSpectrum::constant(0.5).unwrap(),PhaseFunction::Isotropic).unwrap();
        let sources=[fixture(0.75,0.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&grid),TransportOptions::default(),cx).unwrap();
        let reference=SingleScatterReference::new(&scene,cx).unwrap();
        let truth=100.0*(-0.5_f64).exp()*4.0/(4.0*PI)*(1.0/3.375-1.0/3.5);
        let coarse=reference.enclose(cx,segment(),ScatterEnclosureBudget{max_evaluations:1,..budget()}).unwrap();
        contains(coarse,truth);
        let fine=reference.enclose(cx,segment(),budget()).unwrap();contains(fine,truth);
        assert!(fine.numerical.lo>0.0);assert_eq!(fine.termination,ScatterEnclosureTermination::ToleranceMet);
    });
}

#[test]
fn source_and_camera_blockers_bound_the_same_original_cell_integral() {
    with_cx(|cx,_| {
        let density=[1.0,2.0,0.5];let grid=field(cx,&density);
        // Light projection onto x in [1/4,3/8] at z=2 excludes z in [4/5,4/3].
        // The camera x=0 is clear. Independent 90-digit physical-z quadrature.
        let screen=[opaque(vec![[0.25,-0.5,2.0],[0.375,-0.5,2.0],[0.375,0.5,2.0],[0.25,0.5,2.0]],vec![[0,1,2],[0,2,3]])];
        let sources=[fixture(0.75,1.0)];
        let scene=PhotometricScene::try_with_medium(&screen,&sources,PhotometricMedium::Scattering(&grid),TransportOptions::default(),cx).unwrap();
        let value=SingleScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
        contains(value,0.230_728_879_459_855_68);
        assert_eq!(value.termination,ScatterEnclosureTermination::ToleranceMet);
        // For a collinear source, a broad opaque triangle blocks the light leg
        // below it and the camera leg above it: the full direct airlight is zero.
        let screen=[opaque(vec![[-2.0,-2.0,1.5],[2.0,-2.0,1.5],[0.0,3.0,1.5]],vec![[0,1,2]])];
        let sources=[fixture(0.75,0.0)];
        let scene=PhotometricScene::try_with_medium(&screen,&sources,PhotometricMedium::Scattering(&grid),TransportOptions::default(),cx).unwrap();
        let value=SingleScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
        contains(value,0.0);assert_eq!(value.termination,ScatterEnclosureTermination::ToleranceMet);
    });
}

#[test]
fn grid_work_refusals_and_cancellation_never_publish_a_partial_receiver() {
    with_cx(|cx,gate| {
        let density=[1.0,2.0,0.5];let grid=field(cx,&density);let sources=[fixture(0.75,0.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&grid),TransportOptions::default(),cx).unwrap();
        let reference=SingleScatterReference::new(&scene,cx).unwrap();
        let segments=[segment(),ScatterSegment::try_new(Point3::new(0.0,0.0,-1.0),Point3::new(0.0,0.0,2.5)).unwrap()];
        let coarse=ScatterEnclosureBudget{max_evaluations:1,..budget()};
        let mut full=ScatterEnclosureJob::new(&reference,&segments,coarse,cx).unwrap();
        full.advance(cx,usize::MAX).unwrap();
        let mut split=ScatterEnclosureJob::new(&reference,&segments,coarse,cx).unwrap();
        split.advance(cx,1).unwrap();gate.request();
        assert!(matches!(split.advance(cx,1),Err(ScatterEnclosureError::Transport(PhotometryError::Cancelled))));
        assert_eq!(split.results(),&full.results()[..1]);
        with_cx(|resume,_| {split.advance(resume,usize::MAX).unwrap();});
        assert_eq!(split.results(),full.results());
        with_cx(|resume,_| {
            let limited=SingleScatterReference::with_grid_budget(&scene,GridEnclosureBudget{max_nodes_per_segment:1},resume).unwrap();
            let mut job=ScatterEnclosureJob::new(&limited,&segments,coarse,resume).unwrap();
            assert!(matches!(job.advance(resume,1),Err(ScatterEnclosureError::PhysicalMedium(VolumeTransportError::Budget))));
            assert!(job.results().is_empty());
        });
    });
}

#[test]
fn absorbing_only_and_empty_fields_are_zero_without_promoting_total_accuracy() {
    with_cx(|cx,_| {
        let density=[0.0,0.0,0.0];let grid=field(cx,&density);let sources=[fixture(0.75,0.0)];
        for medium in [PhotometricMedium::Scattering(&grid),PhotometricMedium::Absorbing(grid.absorption_grid())] {
            let scene=PhotometricScene::try_with_medium(&[],&sources,medium,TransportOptions::default(),cx).unwrap();
            assert!(matches!(HomogeneousScatterReference::new(&scene,cx),Err(ScatterEnclosureError::HomogeneousOnly)));
            let value=SingleScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),budget()).unwrap();
            contains(value,0.0);assert_eq!(value.numerical.lo,0.0);assert_eq!(value.numerical.hi,0.0);
            assert_eq!(value.numerical_error_bound_cd_m2,0.0);
        }
        let scene=PhotometricScene::try_new(&[],&sources,HomogeneousMedium::try_new(0.05,0.1,PhaseFunction::Rayleigh).unwrap(),TransportOptions::default(),cx).unwrap();
        let old=HomogeneousScatterReference::new(&scene,cx).unwrap();
        let general=SingleScatterReference::new(&scene,cx).unwrap();
        assert_eq!(old.enclose(cx,segment(),budget()).unwrap(),general.enclose(cx,segment(),budget()).unwrap());
    });
}

#[test]
fn grid_csv_binds_original_medium_and_keeps_higher_order_error_unknown() {
    with_cx(|cx,_| {
        let density=[1.0,2.0,0.5];let grid=field(cx,&density);let sources=[fixture(0.75,0.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&grid),TransportOptions::default(),cx).unwrap();
        let result=SingleScatterReference::new(&scene,cx).unwrap().enclose(cx,segment(),ScatterEnclosureBudget{max_evaluations:1,..budget()}).unwrap();
        let mut bytes=Vec::new();EnclosedSingleScatter::write_csv_header(&mut bytes).unwrap();result.write_csv_row(0,&mut bytes).unwrap();
        let text=String::from_utf8(bytes).unwrap();let rows:Vec<_>=text.lines().collect();
        let row:Vec<_>=rows[1].split(',').collect();assert_eq!(row.len(),15);assert_eq!(rows[0].split(',').count(),15);
        assert_eq!(row[1],"physical-grid-spectral-single-scatter");assert_eq!(&row[11..13],&["unknown","unknown"]);
        assert_eq!(row[13],grid.content_identity().to_string());assert_eq!(row[14],"2");
        let changed=[1.0,2.0,0.25];let changed_grid=field(cx,&changed);
        let changed_scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&changed_grid),TransportOptions::default(),cx).unwrap();
        let changed=SingleScatterReference::new(&changed_scene,cx).unwrap().enclose(cx,segment(),ScatterEnclosureBudget{max_evaluations:1,..budget()}).unwrap();
        assert_ne!(result.numerical_model_identity,changed.numerical_model_identity);
    });
}
