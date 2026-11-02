use super::*;
use crate::instances::RigidTransform;
use crate::photometry::{TransportOptions, medium::HomogeneousMedium, transport::PhotometricMedium};
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, ScatteringGrid};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use core::f64::consts::PI;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate,arena,
        StreamKey {seed:173,kernel_id:0x5654_4149,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic),&gate))
}
fn field<'a>(cx:&Cx<'_>,rho:&'a [f64],a:f64,s:f64)->ScatteringGrid<'a> {
    let absorbing=AbsorbingGrid::try_new(cx,
        VolumeGrid::new([1,1,rho.len()],rho,[0.0;3],[1.0,2.0,2.0/rho.len() as f64]),
        MassAbsorptionSpectrum::constant(a).unwrap(),None,
        AbsorbingGridBudget {max_grid_cells:rho.len(),max_steps_per_segment:32}).unwrap();
    ScatteringGrid::try_new(cx,absorbing,MassScatteringSpectrum::constant(s).unwrap(),PhaseFunction::Isotropic).unwrap()
}
fn source(z:f64)->PhotometricSource {
    PhotometricSource::try_new(RigidTransform::try_new([0.0,0.0,0.0,1.0],[0.5,1.0,z]).unwrap(),
        400.0*PI,1.0,AngularDistribution::isotropic(),SourceSpectrum::equal_energy()).unwrap()
}
fn segment()->ScatterSegment {
    ScatterSegment::try_new(Point3::new(0.5,1.0,-1.0),Point3::new(0.5,1.0,2.5)).unwrap()
}
fn limits()->ScatterEnclosureBudget {
    ScatterEnclosureBudget {absolute_error_cd_m2:1e-4,max_evaluations:4096,max_depth:32}
}

#[test]
fn external_fixture_tail_matches_an_independent_geometric_envelope() {
    with_cx(|cx,_| {
        let data=[0.125];let m=field(cx,&data,0.0,0.5);let sources=[source(3.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&m),TransportOptions::default(),cx).unwrap();
        let plain=SingleScatterReference::new(&scene,cx).unwrap();
        assert!(plain.volume_remainder().is_none());
        let reference=SingleScatterReference::new(&scene,cx).unwrap().with_volume_remainder(1,cx).unwrap();
        let bound=reference.volume_remainder().unwrap();
        // D=3, source/grid separation=1, I_peak=100, sigma_s=1/16.
        // h=q=3/16, B1=100*(3/16)/(4*pi), tail=B1*(3/13).
        let first=100.0*(3.0/16.0)/(4.0*PI);
        let tail=first*3.0/13.0;
        assert!(bound.upper_cd_m2()>=tail);
        assert!(bound.upper_cd_m2()<tail*(1.0+1e-11));
        assert!(bound.first_order_upper_cd_m2()>=first);
        assert!(bound.min_source_distance_m().unwrap()<=1.0);
        assert!(bound.max_contraction()<0.188);
        let expected=plain.enclose(cx,segment(),limits()).unwrap();
        let result=reference.enclose(cx,segment(),limits()).unwrap();
        assert_eq!(result.numerical,expected.numerical);
        assert_eq!(result.evaluations,expected.evaluations);
        assert_eq!(result.numerical_model_identity,expected.numerical_model_identity);
        assert_eq!(result.omitted_multiple_scattering_bound_cd_m2,Some(bound.upper_cd_m2()));
        let all=result.all_volume_orders_enclosure().unwrap();
        assert_eq!(all.lo,result.numerical.lo);
        assert!(all.hi>=result.numerical.hi+bound.upper_cd_m2());
        assert!(result.total_error_bound_cd_m2.is_none());
        // Independent first-order z integral, not the tail used as its oracle.
        let single=100.0*0.0625*(-0.125_f64).exp()*(1.0-1.0/3.0)/(4.0*PI);
        assert!(result.numerical.lo<=single && single<=result.numerical.hi);
    });
}

#[test]
fn unknown_singular_or_noncontractive_models_stay_unknown() {
    with_cx(|cx,_| {
        for (rho,a,s,z) in [(10.0,0.0,1.0,3.0),(0.125,0.0,0.5,1.0)] {
            let data=[rho];let m=field(cx,&data,a,s);let sources=[source(z)];
            let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&m),TransportOptions::default(),cx).unwrap();
            let reference=SingleScatterReference::new(&scene,cx).unwrap().with_volume_remainder(1,cx).unwrap();
            assert!(reference.volume_remainder().is_none());
        }
        let sources=[source(3.0)];
        let scene=PhotometricScene::try_new(&[],&sources,
            HomogeneousMedium::try_new(1.0,0.1,PhaseFunction::Isotropic).unwrap(),TransportOptions::default(),cx).unwrap();
        let reference=SingleScatterReference::new(&scene,cx).unwrap().with_volume_remainder(1,cx).unwrap();
        assert!(reference.volume_remainder().is_none());
        assert!(reference.enclose(cx,segment(),limits()).unwrap().all_volume_orders_enclosure().is_none());
    });
}

#[test]
fn exact_zero_scattering_has_zero_tail_without_claiming_physical_accuracy() {
    with_cx(|cx,_| {
        let data=[1.0];let m=field(cx,&data,0.5,0.0);let sources=[source(1.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&m),TransportOptions::default(),cx).unwrap();
        let reference=SingleScatterReference::new(&scene,cx).unwrap().with_volume_remainder(1,cx).unwrap();
        let bound=reference.volume_remainder().unwrap();
        assert_eq!(bound.upper_cd_m2(),0.0);assert!(bound.min_source_distance_m().is_none());
        let value=reference.enclose(cx,segment(),limits()).unwrap();
        assert_eq!(value.all_volume_orders_enclosure(),Some(value.numerical));
        assert!(value.total_error_bound_cd_m2.is_none());
    });
}

#[test]
fn tail_is_global_and_progressive_replay_does_not_refit_it() {
    with_cx(|cx,gate| {
        let data=[0.125,0.0625];let m=field(cx,&data,0.25,0.5);let sources=[source(3.0),source(4.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&m),TransportOptions::default(),cx).unwrap();
        let reference=SingleScatterReference::new(&scene,cx).unwrap().with_volume_remainder(2,cx).unwrap();
        let segments=[segment(),ScatterSegment::try_new(Point3::new(0.4,1.0,-1.0),Point3::new(0.4,1.0,2.5)).unwrap()];
        let coarse=ScatterEnclosureBudget {max_evaluations:2,..limits()};
        let mut full=ScatterEnclosureJob::new(&reference,&segments,coarse,cx).unwrap();full.advance(cx,2).unwrap();
        let mut split=ScatterEnclosureJob::new(&reference,&segments,coarse,cx).unwrap();split.advance(cx,1).unwrap();
        gate.request();assert!(split.advance(cx,1).is_err());assert_eq!(split.results().len(),1);
        with_cx(|resume,_| {split.advance(resume,1).unwrap();});
        assert_eq!(split.results(),full.results());
        assert_eq!(full.results()[0].volume_remainder,full.results()[1].volume_remainder);
        assert_eq!(full.results()[0].termination,ScatterEnclosureTermination::PanelBudget);
        // Obtaining a tail never upgrades the requested single-order tolerance.
        assert!(full.results()[0].omitted_multiple_scattering_bound_cd_m2.is_some());
    });
}

#[test]
fn preparation_limits_and_cancellation_do_not_publish_partial_bounds() {
    with_cx(|cx,gate| {
        let data=[0.125,0.125];let m=field(cx,&data,0.0,0.5);let sources=[source(3.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&m),TransportOptions::default(),cx).unwrap();
        let a=SingleScatterReference::new(&scene,cx).unwrap();
        assert!(matches!(a.with_volume_remainder(1,cx),Err(ScatterEnclosureError::Budget)));
        let b=SingleScatterReference::new(&scene,cx).unwrap();gate.request();
        assert!(matches!(b.with_volume_remainder(2,cx),Err(ScatterEnclosureError::Transport(PhotometryError::Cancelled))));
    });
}

#[test]
fn csv_populates_only_the_existing_volume_tail_column() {
    with_cx(|cx,_| {
        let data=[0.125];let m=field(cx,&data,0.0,0.5);let sources=[source(3.0)];
        let scene=PhotometricScene::try_with_medium(&[],&sources,PhotometricMedium::Scattering(&m),TransportOptions::default(),cx).unwrap();
        let reference=SingleScatterReference::new(&scene,cx).unwrap().with_volume_remainder(1,cx).unwrap();
        let result=reference.enclose(cx,segment(),ScatterEnclosureBudget {max_evaluations:1,..limits()}).unwrap();
        let mut csv=Vec::new();EnclosedSingleScatter::write_csv_header(&mut csv).unwrap();result.write_csv_row(0,&mut csv).unwrap();
        let text=String::from_utf8(csv).unwrap();let rows:Vec<_>=text.lines().collect();let row:Vec<_>=rows[1].split(',').collect();
        assert_eq!(rows[0].split(',').count(),15);assert_eq!(row.len(),15);
        assert_eq!(row[11].parse::<f64>().unwrap(),result.omitted_multiple_scattering_bound_cd_m2.unwrap());
        assert_eq!(row[12],"unknown");assert_eq!(row[13],m.content_identity().to_string());
        assert_eq!(row[14],"2");assert_eq!(result.volume_remainder.unwrap().semantics_version(),1);
    });
}
