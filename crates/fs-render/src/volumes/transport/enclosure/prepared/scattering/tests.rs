use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGridBudget,CellPalette,CellPaletteBudget};
use crate::volumes::transport::scattering::{MassScatteringSpectrum,ScatteringOptics};
use fs_exec::{Budget,CancelGate,ExecMode,StreamKey};
use std::sync::Arc;
use core::f64::consts::PI;
fn with_cx<R>(f:impl FnOnce(&Cx<'_>,&CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena|f(&Cx::new(&gate,arena,StreamKey{seed:173,kernel_id:0x5343_4f4c,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic),&gate))
}
fn ends(z:f64)->[Interval;3] { [0.5,0.5,z].map(Interval::point) }
#[test]
fn independent_cells_keep_scattering_and_phase_paired_at_each_wavelength() {
    with_cx(|cx,_| {
        let data=[1.0,2.0,4.0];let cap=CellPaletteBudget{max_models:2,max_cells:3};
        let a=AbsorbingGrid::try_new(cx,VolumeGrid::new([1,1,3],&data,[0.0;3],[1.0;3]),
            MassAbsorptionSpectrum::constant(0.1).unwrap(),None,AbsorbingGridBudget{max_grid_cells:3,max_steps_per_segment:32}).unwrap();
        let phase=PhaseFunction::spectral_table(&[(380.0,PhaseFunction::Isotropic),(780.0,PhaseFunction::Rayleigh)]).unwrap();
        let optics=CellPalette::try_new(cx,vec![
            ScatteringOptics{scattering:MassScatteringSpectrum::constant(0.25).unwrap(),phase:PhaseFunction::Isotropic},
            ScatteringOptics{scattering:MassScatteringSpectrum::rayleigh(400.0,0.5).unwrap(),phase},
        ],Arc::from([0,1,0]),cap).unwrap();
        let medium=ScatteringGrid::try_new_indexed(cx,a,optics).unwrap();
        let prepared=medium.prepare_extinction_enclosure(cx,&[400.0,600.0],6).unwrap()
            .into_single_scattering(cx,6).unwrap();
        assert!(prepared.has_scattering());
        let result=prepared.scattering_integral_box(cx,ends(-1.0),ends(4.0),Interval::point(1.0),GridEnclosureBudget::default()).unwrap();
        for (lane,lambda) in [400.0,600.0].into_iter().enumerate() {
            let r=400.0/lambda;let t=(lambda-380.0)/400.0;
            let truth=(1.25+r*r*r*r*(1.0+0.5*t))/(4.0*PI);
            assert!(result[lane].contains(truth),"{result:?} {truth}");
            assert!(result[lane].width()<1e-12);
        }
        let reverse=prepared.scattering_integral_box(cx,ends(4.0),ends(-1.0),Interval::point(1.0),GridEnclosureBudget::default()).unwrap();
        for lane in 0..2 {assert!(result[lane].intersect(reverse[lane]).is_some());}
    });
}
#[test]
fn a_thin_occupied_cell_is_integrated_even_without_a_midpoint_sample() {
    with_cx(|cx,_| {
        let data=[0.0,0.0,8.0,0.0,0.0,0.0,0.0,0.0];
        let a=AbsorbingGrid::try_new(cx,VolumeGrid::new([1,1,8],&data,[0.0,0.0,0.25],[1.0,1.0,0.125]),
            MassAbsorptionSpectrum::constant(0.0).unwrap(),None,AbsorbingGridBudget{max_grid_cells:8,max_steps_per_segment:32}).unwrap();
        let medium=ScatteringGrid::try_new(cx,a,MassScatteringSpectrum::constant(0.5).unwrap(),PhaseFunction::Isotropic).unwrap();
        let prepared=medium.prepare_extinction_enclosure(cx,&[500.0],1).unwrap().into_single_scattering(cx,2).unwrap();
        let result=prepared.scattering_integral_box(cx,ends(0.0),ends(2.0),Interval::point(0.0),GridEnclosureBudget::default()).unwrap();
        assert!(result[0].contains(0.5/(4.0*PI)));assert!(result[0].lo()>0.0);
    });
}
#[test]
fn response_preserves_extinction_bits_and_declared_work_limits() {
    with_cx(|cx,gate| {
        let data=[1.0;8];let a=AbsorbingGrid::try_new(cx,VolumeGrid::new([1,1,8],&data,[0.0;3],[1.0;3]),
            MassAbsorptionSpectrum::constant(0.125).unwrap(),None,AbsorbingGridBudget{max_grid_cells:8,max_steps_per_segment:32}).unwrap();
        let medium=ScatteringGrid::try_new(cx,a,MassScatteringSpectrum::constant(0.25).unwrap(),PhaseFunction::Rayleigh).unwrap();
        let old=medium.prepare_extinction_enclosure(cx,&[400.0,600.0],2).unwrap();
        let expected=old.transmittance_box(cx,ends(-1.0),ends(9.0),GridEnclosureBudget::default()).unwrap();
        assert!(matches!(old.into_single_scattering(cx,3),Err(VolumeTransportError::Budget)));
        let new=medium.prepare_extinction_enclosure(cx,&[400.0,600.0],2).unwrap().into_single_scattering(cx,4).unwrap();
        assert_eq!(expected,new.extinction().transmittance_box(cx,ends(-1.0),ends(9.0),GridEnclosureBudget::default()).unwrap());
        assert!(matches!(new.scattering_integral_box(cx,ends(-1.0),ends(9.0),Interval::point(0.0),GridEnclosureBudget{max_nodes_per_segment:1}),Err(VolumeTransportError::Budget)));
        gate.request();
        assert!(matches!(new.scattering_integral_box(cx,ends(-1.0),ends(9.0),Interval::point(0.0),GridEnclosureBudget::default()),Err(VolumeTransportError::Cancelled)));
    });
}
