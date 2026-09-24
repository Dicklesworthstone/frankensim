use super::*;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use fs_couple::bernoulli_aperture::dynamic::relaxation::{PlateRelaxationSpec,PlateRelaxationRegion,InitialApertureMemory};
use fs_couple::bernoulli_aperture::plate::closure::PlateClosureSpec;
use fs_couple::bernoulli_aperture::performance::file::PlateValvePerformance;
use fs_couple::bernoulli_aperture::tube::{ApertureTube,TubeDrive,UniformTubeSpec};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_material::{gas::GasState,visco::GeneralizedMaxwell};
use fs_scenario::gesture::GestureSchedule;
const INPUT:&str=include_str!("../../examples/plate-valve.performance");
fn load(text:&str)->PlateValvePerformance {PlateValvePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new()).unwrap()}
fn render(text:&str)->Vec<f64> {
    let p=load(text);let mut values=vec![0.0;p.info().samples as usize];let mut r=p.into_renderer();
    for out in values.chunks_mut(37) {r.block(out).unwrap();}values
}
#[test]
fn complete_mesh_material_and_gesture_file_matches_independent_physical_construction() {
    let source=load(INPUT);let info=source.info();
    assert_eq!((info.nodes,info.triangles,info.sections,info.memory_branches),(21,24,1,2));
    let (chart,mut options)=fixture(4e9,900.0);options.damping_ratio=0.0;
    options.eigenvalue_window=((core::f64::consts::TAU*0.2).powi(2),(core::f64::consts::TAU*450.0).powi(2));
    let plate=PlateApertureReduction::from_chart(chart,options,&CancelGate::new()).unwrap();
    let closure=PlateClosureSpec{nodal_rest_gap_m:plate.chart().mesh.nodes.iter().map(|p|0.0002*(0.05+1.9*p.1/0.01)).collect(),
        lay_triangles:(0..24).collect(),stiffness_pa_per_m_alpha:1e12,alpha:2.0,internal_loss_s_per_m:0.5,
        provenance:"independent construction, identical supplied physical fields".into(),max_penetration_m:0.0002};
    let gas=GasState::try_new_moist_air(293.15,101325.0,0.0).unwrap();
    let spec=UniformTubeSpec{length_m:0.25,radius_m:0.007,sound_speed_m_s:gas.sound_speed,
        terminal_reflection:-0.8,max_length_error_m:0.002,max_wave_memory_bytes:1048576};
    let valve=DynamicAperture::from_plate_with_closure(plate,closure,gas.density,spec.characteristic_impedance(gas.density).unwrap(),
        1.0/96000.0,9602,ApertureState{opening_m:0.0002,opening_velocity_m_s:0.0}).unwrap()
        .with_plate_relaxation(PlateRelaxationSpec{regions:vec![PlateRelaxationRegion{triangles:(0..24).collect(),
            material:GeneralizedMaxwell::new(4e9,vec![(2e9,0.001),(1e9,0.01)]).unwrap(),poisson_ratio:0.3,band_hz:(0.0,10000.0),
            provenance:"independent supplied material".into()}],max_branches:64,max_dt_over_tau:0.1,max_angular_step:0.5},InitialApertureMemory::Relaxed).unwrap();
    let mut direct=ApertureTube::new(valve,spec).unwrap();
    assert_eq!(info.represented_tube_length_m,direct.represented_length_m());
    let gesture=GestureSchedule::from_canonical_bytes(INPUT.split_once("\nschedule\n").unwrap().1.as_bytes()).unwrap();
    let mut expected=Vec::new();
    for n in 0..9602 {
        expected.push(direct.step(TubeDrive{upstream_pressure_pa:gesture.sample("mouth",n*700/96000).unwrap(),body_flow_m3_s:0.0}).unwrap().aperture.bore_pressure_pa);
    }
    let actual=render(INPUT);
    assert_eq!(actual.iter().map(|p|p.to_bits()).collect::<Vec<_>>(),expected.iter().map(|p|p.to_bits()).collect::<Vec<_>>());
    let specimen=source.renderer().system().aperture().plate_reduction().unwrap();
    assert!(specimen.material_chart().is_none(),"numeric file must not mint material receipts");
    assert_eq!(specimen.chart().mesh.nodes.len(),21);
}

#[test]
fn actual_regional_thickness_and_supplied_history_change_pressure_not_only_metadata() {
    let a=render(INPUT);
    let heavy=INPUT.replace("0.0003 900 4000000000","0.0003 1800 4000000000");
    let b=render(&heavy);
    let initial=INPUT.replace("initial 0.0002 0","initial 0.00021 0");
    let c=render(&initial);let d=render(&initial.replace("memory_initial relaxed","memory_initial unrelaxed"));
    assert!(a.iter().zip(&b).any(|(x,y)|(x-y).abs()>1e-6));
    assert!(c.iter().zip(&d).any(|(x,y)|(x-y).abs()>1e-6));
    let mut regional=INPUT.replace("sections 1\nsection isotropic 0.0003 900 4000000000 0.3",
        "sections 2\nsection isotropic 0.0003 900 4000000000 0.3\nsection isotropic 0.0002 900 4000000000 0.3")
        .replace("relaxation 1","relaxation 2")
        .replace("compile_limits", "region 1 4000000000 0.3 0 10000 1\nbranch 1000000000 0.004\ncompile_limits");
    regional=regional.lines().map(|line| {
        if line.starts_with("triangle ") {
            let mut fields:Vec<_>=line.split_whitespace().map(str::to_string).collect();
            let node:usize=fields[1].parse().unwrap();if node%7>=3 {fields[4]="1".into();}fields.join(" ")
        } else {line.into()}
    }).collect::<Vec<_>>().join("\n")+"\n";
    let p=load(&regional);assert_eq!(p.info().sections,2);assert_eq!(p.info().memory_branches,3);
    let r=p.renderer().system().aperture().relaxation().unwrap();
    assert_eq!(r.spec().regions.len(),2);
    assert!(r.region_bending_stiffness_n_m().iter().all(|x|*x>0.0));
    assert!(render(&regional).iter().zip(a).any(|(x,y)|(x-y).abs()>1e-6));
}

#[test]
fn invalid_mesh_sections_history_counts_and_out_of_window_gestures_are_not_silently_repaired() {
    for text in [
        INPUT.replace("nodes 21","nodes 18446744073709551615"),
        INPUT.replace("triangles 24","triangles 2049"),
        INPUT.replace("triangle 0 1 8 0 1","triangle 0 1 9999 0 1"),
        INPUT.replace("triangle 0 1 8 0 1","triangle 0 8 1 0 1"),
        INPUT.replace("triangle 0 1 8 0 1","triangle 0 1 8 1 1"),
        INPUT.replace("triangle 0 1 8 0 1","triangle 0 1 8 0 2"),
        INPUT.replace("edge 6 13","edge 0 8"),
        INPUT.replace("region 0 4000000000","region 0 3000000000"),
        INPUT.replace("memory_initial relaxed","memory_initial explicit 3 0 0 0"),
        INPUT.replace("memory_initial relaxed","memory_initial guessed"),
        INPUT.replace("plate clamped 0 0 ","plate clamped 0 0.02 "),
        INPUT.replace("region 0 4000000000 0.3 0 10000 2","region 0 4000000000 0.3 0 10000 65"),
        INPUT.replace("events\t5","events\t18446744073709551615"),
        INPUT.replace("event\t6e-2","event\t2e0"),
        INPUT.replace("ambient 293.15","ambient NaN"),
        INPUT.replace("observation inlet","observation exterior"),
        format!("{INPUT}ignored\n"),
    ] {assert!(PlateValvePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new()).is_err(),"accepted {text}");}
    let gate=CancelGate::new();gate.request();
    assert!(PlateValvePerformance::from_bytes(INPUT.as_bytes(),37,&gate).is_err());
    assert!(PlateValvePerformance::from_bytes(INPUT.as_bytes(),0,&CancelGate::new()).is_err());
}

#[path = "file_radiation.rs"]
mod radiation;

#[path = "duct_graph.rs"]
mod duct_graph;
