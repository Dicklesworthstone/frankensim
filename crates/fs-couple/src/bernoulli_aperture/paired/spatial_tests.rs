//! Physical-path regressions; numerical plate sections are authored, not measured cane.
use super::*;
use crate::bernoulli_aperture::plate::{PlateApertureOptions, closure::PlateClosureSpec};
use crate::bernoulli_aperture::dynamic::relaxation::PlateRelaxationRegion;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::bernoulli_aperture::paired_network::PairedApertureNetwork;
use crate::bernoulli_aperture::tube::TubeDrive;
use fs_exec::CancelGate;
use fs_material::{gas::GasSpec, visco::GeneralizedMaxwell};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};

const RATE: u32=100_000;
fn gas() -> GasState { GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap() }
fn plate(e:f64,density:f64,length:f64) -> PlateApertureReduction {
    let mesh=PlateMesh::rectangle(length,0.01,6,2);
    let root=mesh.nodes.iter().enumerate().filter(|(_,p)|p.0==0.0).map(|(i,_)|i).collect();
    let slit_edges=mesh.boundary_edges().into_iter().filter(|&(a,b)|
        (mesh.nodes[a].0-length).abs()<1e-12 && (mesh.nodes[b].0-length).abs()<1e-12)
        .map(|(a,b)|[a,b]).collect();
    let chart=PlateChart::with_boundary_and_regions(mesh,
        PlateSection::isotropic(e,0.3,0.0003,density).unwrap(),root,vec![]).unwrap();
    PlateApertureReduction::from_chart(chart,PlateApertureOptions {
        assembly:AssemblyOptions {pretension:0.0,support:EdgeSupport::Clamped},
        eigenvalue_window:(1.0,(core::f64::consts::TAU*450.0).powi(2)),mode_index:0,
        eigensolver:SliceOptions::default(),max_nodes:100,max_triangles:100,slit_edges,
        rest_opening_m:0.0002,damping_ratio:0.0,max_slit_mode_variation:0.25,max_slope:0.1,
    },&CancelGate::new_clock_free()).unwrap()
}
fn profile(p:&PlateApertureReduction) -> PlateClosureSpec {
    PlateClosureSpec {
        nodal_rest_gap_m:p.chart().mesh.nodes.iter().map(|n|0.00002+0.00036*n.1/0.01).collect(),
        lay_triangles:(0..p.chart().mesh.tris.len()).collect(),
        stiffness_pa_per_m_alpha:1e12,alpha:2.0,internal_loss_s_per_m:0.5,
        provenance:"authored matching blade/blade contact for numerical regression".into(),
        max_penetration_m:0.0005,
    }
}
fn material(p:&PlateApertureReduction,e:f64,tau:f64) -> PlateRelaxationSpec {
    PlateRelaxationSpec {regions:vec![PlateRelaxationRegion {
        triangles:(0..p.chart().mesh.tris.len()).collect(),
        material:GeneralizedMaxwell::new(e,vec![(0.5*e,tau)]).unwrap(),poisson_ratio:0.3,
        band_hz:(0.0,1000.0),provenance:"synthetic spatial blade material".into(),
    }],max_branches:4,max_dt_over_tau:0.1,max_angular_step:0.1}
}
fn raw(plates:[PlateApertureReduction;2],z:f64) -> PairedAperture {
    let contact=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],0.0,2.0,
        "explicit uniform no-contact comparison, replaced by selected spatial law".into()).unwrap();
    // The spatial momentum residual has a measured evaluation floor near
    // 3.9e-15 m/s at |v| = 1.7e-3 m/s (ts1 probe, 2026-10-17): a 1e-13
    // relative request (1.5e-15 m/s) is below what the arithmetic resolves.
    // 1e-11 keeps tenfold headroom, as the other paired fixtures do.
    PairedAperture::from_plates(plates,gas(),z,1.0/f64::from(RATE),4096,
        PairedApertureLimits {maximum_abs_displacement_m:0.001,maximum_abs_velocity_m_s:10.0,
            maximum_energy_j:1.0,max_evaluations:256,velocity_tolerance_m_s:1e-15,relative_tolerance:1e-11},
        [BladeState::default();2],contact).unwrap()
}
fn model(z:f64) -> PairedAperture {
    let plates=[plate(4e9,900.0,0.025),plate(8e9,1300.0,0.025)];
    let closure=profile(&plates[0]);
    let laws=[material(&plates[0],4e9,0.001),material(&plates[1],8e9,0.003)];
    raw(plates,z).with_blade_relaxation(0,laws[0].clone(),InitialApertureMemory::Unrelaxed).unwrap()
        .with_blade_relaxation(1,laws[1].clone(),InitialApertureMemory::Relaxed).unwrap()
        .with_plate_closure(closure).unwrap()
}

#[test]
fn spatial_contact_flow_and_both_material_histories_share_one_accepted_step() {
    let mut m=model(1e6);
    let mut contacted=false; let mut partially_open=false;
    for n in 0..800 {
        let before=m.blades.map(|s|s.displacement_m);
        let f=m.step(ApertureDrive {upstream_pressure_pa:if n<600 {50.0}else{0.0},
            ..ApertureDrive::default()}).unwrap();
        let after=f.blades.map(|s|s.displacement_m);
        let mid=[f64::midpoint(before[0],after[0]),f64::midpoint(before[1],after[1])];
        let c=m.plate_closure().unwrap();
        let area=c.open_area_m2(mid).unwrap(); let probe=c.probe(after).unwrap();
        assert!((f.aperture.midpoint_opening_m*m.spec.width_m-area).abs()<1e-18);
        contacted|=f.contact_force_n>0.0;
        partially_open|=area>0.0 && probe.active_lay_points>0;
        assert!(f.momentum_residual_n.iter().all(|r|r.abs()<1e-9));
        let scale=f.aperture.stored_energy_j+f.aperture.storage_change_j.abs()
            +f.aperture.dissipated_energy_j+f.aperture.pressure_work_j.abs();
        assert!(f.aperture.balance_residual_j().abs()<3e-9*scale.max(1e-30));
        assert_eq!(m.accepted_steps(),n+1);
    }
    assert!(contacted && partially_open);
    assert!(m.blade_relaxation(0).unwrap().dissipated_energy_j()>0.0);
    assert!(m.blade_relaxation(1).unwrap().dissipated_energy_j()>0.0);
    assert!(m.stored_energy_j()>0.0,"pressure release must not reset vibration or memory");
}

#[test]
fn spatial_geometry_changes_response_without_changing_either_blade_material() {
    let a=plate(4e9,900.0,0.025); let b=plate(8e9,1300.0,0.025);
    let varied=profile(&a);
    let mut flat=varied.clone();flat.nodal_rest_gap_m.fill(0.0002);
    let mut x=raw([a.clone(),b.clone()],1e6).with_plate_closure(varied).unwrap();
    let mut y=raw([a,b],1e6).with_plate_closure(flat).unwrap();
    assert_eq!(x.spec(),y.spec());
    let mut changed=false;
    for _ in 0..600 {
        let drive=ApertureDrive {upstream_pressure_pa:50.0,..ApertureDrive::default()};
        let a=x.step(drive).unwrap(); let b=y.step(drive).unwrap();
        changed|=(a.aperture.bore_pressure_pa-b.aperture.bore_pressure_pa).abs()>1e-8;
    }
    assert!(changed,"clearance data must enter pressure/structure feedback, not only observations");
}

#[test]
fn refused_spatial_network_step_keeps_both_blades_memories_and_waves() {
    let g=gas();
    let graph=TubeNetworkSpec {nodes:vec![NetworkNode::Inlet,NetworkNode::Termination {reflection:-0.7}],
        sections:vec![TubeSection {nodes:[0,1],length_m:0.25,radius_m:0.007,max_length_error_m:0.005}],
        sound_speed_m_s:g.sound_speed,max_wave_memory_bytes:1<<20};
    let mut aperture=model(graph.inlet_impedance_with_gases(&[g]).unwrap());
    // No speculative root evaluation may be published as a physical step.
    aperture.spec.limits.max_evaluations=1;
    let mut n=PairedApertureNetwork::new(aperture,graph,vec![g],&CancelGate::new_clock_free()).unwrap();
    let states=*n.aperture().blade_states();let energy=n.stored_energy_j();
    let history:Vec<_>=(0..2).map(|i|n.aperture().blade_relaxation(i).unwrap().memory_sqrt_j().to_vec()).collect();
    let node=n.node_frame(1).copied();
    assert!(n.step(TubeDrive {upstream_pressure_pa:50.0,body_flow_m3_s:0.0}).is_err());
    assert_eq!(*n.aperture().blade_states(),states);assert_eq!(n.stored_energy_j(),energy);
    assert_eq!(n.aperture().accepted_steps(),0);assert_eq!(n.node_frame(1).copied(),node);
    for i in 0..2 { assert_eq!(n.aperture().blade_relaxation(i).unwrap().memory_sqrt_j(),history[i]); }
}

#[test]
fn mismatched_surfaces_late_selection_and_initial_penetration_refuse() {
    let a=plate(4e9,900.0,0.025);let b=plate(4e9,900.0,0.030);let p=profile(&a);
    assert!(a.compile_paired_closure(&b,p.clone()).is_err());
    let mut bad=p.clone();bad.max_penetration_m=1e-6;bad.nodal_rest_gap_m.fill(-1e-5);
    assert!(raw([a.clone(),a.clone()],1e6).with_plate_closure(bad).is_err());
    let mut late=raw([a.clone(),a.clone()],1e6);late.step(ApertureDrive::default()).unwrap();
    assert!(late.with_plate_closure(p.clone()).is_err());
    let once=raw([a.clone(),a],1e6).with_plate_closure(p.clone()).unwrap();
    assert!(once.with_plate_closure(p).is_err());
}
