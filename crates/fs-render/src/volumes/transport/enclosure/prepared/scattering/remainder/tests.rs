use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGridBudget, CellPalette, CellPaletteBudget};
use crate::volumes::transport::scattering::{ScatteringOptics, MassScatteringSpectrum};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::sync::Arc;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 173, kernel_id: 0x5441_494c, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic), &gate));
}
fn medium<'a>(cx: &Cx<'_>, data: &'a [f64], a: f64, s: f64) -> ScatteringGrid<'a> {
    let grid = AbsorbingGrid::try_new(cx,
        VolumeGrid::new([1,1,data.len()], data, [0.0;3], [1.0,2.0,2.0/data.len() as f64]),
        MassAbsorptionSpectrum::constant(a).unwrap(), None,
        AbsorbingGridBudget { max_grid_cells: data.len(), max_steps_per_segment: 32 }).unwrap();
    ScatteringGrid::try_new(cx, grid, MassScatteringSpectrum::constant(s).unwrap(), PhaseFunction::Isotropic).unwrap()
}
fn tight(upper: f64, exact: f64) {
    assert!(upper >= exact, "{upper} < {exact}");
    assert!(upper <= exact*(1.0+2e-12), "{upper} too wide for {exact}");
}

#[test]
fn finite_support_contracts_even_without_absorption() {
    with_cx(|cx,_| {
        // D=3, sigma_s=1/16: the independent thin-chord bound is 3/16.
        let m=medium(cx,&[0.125],0.0,0.5);
        let p=m.prepare_extinction_enclosure(cx,&[450.0,650.0],2).unwrap()
            .into_single_scattering(cx,4).unwrap();
        let e=p.volume_contraction_envelope(cx,1).unwrap().unwrap();
        tight(e.diameter_m(),3.0);
        for lane in 0..2 {
            tight(e.scattering_integral_upper()[lane],3.0/16.0);
            tight(e.contraction_upper()[lane],3.0/16.0);
            // Independent sum B1*(3/16)/(1-3/16), with B1=13.
            tight(e.tail_upper(lane,13.0).unwrap(),3.0);
        }
        assert_eq!(e.bounds_m()[0].lo(),0.0);
        assert!(e.bounds_m()[2].contains(2.0));
    });
}

#[test]
fn original_absorption_gain_beats_the_diameter_bound_in_dense_media() {
    with_cx(|cx,_| {
        let m=medium(cx,&[100.0],3.0,1.0);
        let p=m.prepare_extinction_enclosure(cx,&[500.0],1).unwrap()
            .into_single_scattering(cx,2).unwrap();
        let e=p.volume_contraction_envelope(cx,1).unwrap().unwrap();
        tight(e.contraction_upper()[0],0.25);
        tight(e.tail_upper(0,3.0).unwrap(),1.0);
    });
}

#[test]
fn indexed_fields_ignore_empty_conservative_materials_and_keep_colored_peaks() {
    with_cx(|cx,_| {
        let data=[1.0,0.0]; let cap=CellPaletteBudget { max_models:2,max_cells:2 };
        let a=CellPalette::try_new(cx,vec![MassAbsorptionSpectrum::constant(3.0).unwrap(),
            MassAbsorptionSpectrum::constant(0.0).unwrap()],Arc::from([0,1]),cap).unwrap();
        let a=AbsorbingGrid::try_new_indexed(cx,VolumeGrid::new([1,1,2],&data,[0.0;3],[1.0;3]),a,None,
            AbsorbingGridBudget {max_grid_cells:2,max_steps_per_segment:32}).unwrap();
        let s=CellPalette::try_new(cx,vec![
            ScatteringOptics { scattering:MassScatteringSpectrum::rayleigh(400.0,1.0).unwrap(),phase:PhaseFunction::Rayleigh },
            ScatteringOptics { scattering:MassScatteringSpectrum::constant(100.0).unwrap(),phase:PhaseFunction::henyey_greenstein(0.9).unwrap() },
        ],Arc::from([0,1]),cap).unwrap();
        let m=ScatteringGrid::try_new_indexed(cx,a,s).unwrap();
        let p=m.prepare_extinction_enclosure(cx,&[400.0,700.0],8).unwrap().into_single_scattering(cx,8).unwrap();
        let e=p.volume_contraction_envelope(cx,2).unwrap().unwrap();
        for peak in e.phase_peak_per_sr() {tight(*peak,3.0/(8.0*core::f64::consts::PI));}
        assert!(e.contraction_upper().iter().all(|q| *q<0.3));
        assert!(e.scattering_integral_upper()[1]<e.scattering_integral_upper()[0]);
    });
}

#[test]
fn conservative_dense_fields_do_not_invent_a_contraction() {
    with_cx(|cx,_| {
        let m=medium(cx,&[10.0],0.0,1.0);
        let p=m.prepare_extinction_enclosure(cx,&[500.0],1).unwrap().into_single_scattering(cx,2).unwrap();
        let e=p.volume_contraction_envelope(cx,1).unwrap().unwrap();
        assert!(e.contraction_upper()[0]>=1.0);
        assert!(e.tail_upper(0,1.0).is_none());
        assert_eq!(e.tail_upper(0,0.0),Some(0.0));
    });
}

#[test]
fn tail_arithmetic_retains_positive_underflow_and_rejects_bad_inputs() {
    let e=VolumeContractionEnvelope::<1> {
        bounds_m:[Interval::new(0.0,1.0);3],diameter_m:2.0,
        scattering_integral_upper:[0.5],phase_peak_per_sr:[1.0],contraction_upper:[0.5],
    };
    tight(e.tail_upper(0,8.0).unwrap(),8.0);
    assert!(e.tail_upper(0,f64::from_bits(1)).unwrap()>0.0);
    for value in [-1.0,f64::NAN,f64::INFINITY] {assert!(e.tail_upper(0,value).is_none());}
    assert!(e.tail_upper(1,1.0).is_none());
    let near=VolumeContractionEnvelope {contraction_upper:[1.0_f64.next_down()],..e};
    assert!(near.tail_upper(0,f64::MAX).is_none());
}

#[test]
fn empty_fields_and_budget_cancellation_keep_their_distinct_meanings() {
    with_cx(|cx,gate| {
        let m=medium(cx,&[0.0,0.0],0.0,1.0);
        let p=m.prepare_extinction_enclosure(cx,&[500.0],1).unwrap().into_single_scattering(cx,2).unwrap();
        let e=p.volume_contraction_envelope(cx,2).unwrap().unwrap();
        assert_eq!(e.contraction_upper(),&[0.0]);
        assert_eq!(e.tail_upper(0,f64::MAX),Some(0.0));
        assert!(matches!(p.volume_contraction_envelope(cx,1),Err(VolumeTransportError::Budget)));
        gate.request();
        assert!(matches!(p.volume_contraction_envelope(cx,2),Err(VolumeTransportError::Cancelled)));
    });
}
