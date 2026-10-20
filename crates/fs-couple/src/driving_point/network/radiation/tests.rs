use super::*;
use crate::driving_point::network::tests::{graph, schema, gas};

fn budget()->ExteriorNetworkBudget {
    ExteriorNetworkBudget{max_panels:128,max_field_entries:4096,max_solve_terms:2_000_000}
}
fn exterior()->(SpherePanels,Vec<OpeningPatch>) {
    let surface=SpherePanels::icosphere(0.03,1).unwrap();
    let a=surface.normals()[0];
    let opposite=surface.normals().iter().enumerate().min_by(|(_,b),(_,c)| {
        let dot=|v:&[f64;3]|a[0]*v[0]+a[1]*v[1]+a[2]*v[2];
        dot(b).total_cmp(&dot(c))
    }).unwrap().0;
    let patches=vec![OpeningPatch{port:"right".into(),panels:vec![0]},
        OpeningPatch{port:"left".into(),panels:vec![opposite]}];
    (surface,patches)
}
fn model(gas:GasState,points:Vec<[f64;3]>)->RadiatingNetwork {
    let (surface,patches)=exterior();
    let driver=NetworkDrivingPoint::new_with_coupled_ports(graph(),"inlet".into(),vec![],
        patches.iter().map(|p|p.port.clone()).collect(),schema()).unwrap();
    RadiatingNetwork::new(driver,surface,patches,gas,points,0.001,ContentHash([7;32]),budget(),
        &CancelGate::new_clock_free()).unwrap()
}
fn receivers()->Vec<[f64;3]> {vec![[0.13,0.05,0.11],[-0.1,0.12,0.04]]}
fn close(a:C64,b:C64) {assert!((a-b).abs()<=2e-9*b.abs().max(1e-12),"{a:?} != {b:?}");}

#[test]
fn same_exterior_solution_supplies_mutual_load_and_coherent_microphone() {
    let m=model(gas(),receivers()); let gate=CancelGate::new_clock_free();
    let omega=core::f64::consts::TAU*1200.0; let flow=C64::new(1e-6,-2e-6);
    let sample=m.solve(omega,LossModel::Lossless,flow,&gate).unwrap();
    assert_eq!(sample.boundary.labels,m.driving_point().coupled_ports());
    assert!(sample.boundary.z[1].abs()>0.0);
    assert_eq!(sample.inlet.outward_volume_flow,-flow);
    assert!(sample.receiver_fields.iter().any(|f|f.pressure_pa.abs()>1e-12));
    let mut velocity=vec![C64::ZERO;m.surface().areas().len()];
    for (patch,&q) in m.patches().iter().zip(&sample.outward_flows_m3_s) {
        let area=patch.panels.iter().map(|&i|m.surface().areas()[i]).sum::<f64>();
        for &i in &patch.panels {velocity[i]=q.scale(area.recip());}
    }
    // Independent assembled boundary solve, not summing the unit-source results.
    let medium=Medium{density:m.gas().density,sound_speed:m.gas().sound_speed};
    let direct=helmholtz::solve_radiation(m.surface(),omega/medium.sound_speed,medium,&velocity,Formulation::BurtonMiller).unwrap();
    let fields=laplace::observe_frequency(m.surface(),&direct,m.receivers(),0.001).unwrap();
    for (a,b) in sample.receiver_fields.iter().zip(fields) {
        close(a.pressure_pa,b.pressure_pa);
        for i in 0..3 {close(a.gradient_pa_per_m[i],b.gradient_pa_per_m[i]);}
    }
    let surface_power=0.5*direct.pressure.iter().zip(&velocity).zip(m.surface().areas())
        .map(|((p,v),a)|(*p*v.conj()).re*a).sum::<f64>();
    assert!((sample.exterior_surface_power_w-surface_power).abs()<2e-9*surface_power.abs().max(1e-20));
    let inlet_power=0.5*(sample.inlet.pressure*flow.conj()).re;
    assert!((inlet_power-sample.exterior_surface_power_w).abs()
        < 2e-9*inlet_power.abs().max(1e-20), "lossless bore must supply the same exterior work");
    for (row,patch) in m.patches().iter().enumerate() {
        let p=(0..m.patches().len()).fold(C64::ZERO,|a,j|a+sample.boundary.z[row*m.patches().len()+j]*sample.outward_flows_m3_s[j]);
        let index=if patch.port=="left" {1} else {2};
        close(sample.response.port_pressures[index]*flow,p);
    }
    // Dropping mutual entries changes the source load, not merely its observer.
    let mut diagonal=sample.boundary.clone(); diagonal.z[1]=C64::ZERO;diagonal.z[2]=C64::ZERO;
    let independent=m.driving_point().response_coupled(&m.gas(),omega,LossModel::Lossless,&diagonal).unwrap();
    assert!((sample.response.impedance-independent.impedance).abs()>1e-6*sample.response.impedance.abs());
}

#[test]
fn shared_gas_scales_both_mechanical_loading_and_actual_sound() {
    let gate=CancelGate::new_clock_free(); let w=core::f64::consts::TAU*1200.0;
    let a=model(gas(),receivers());
    let lower=GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(),293.15,50_662.5).unwrap();
    let b=model(lower,receivers());
    let fa=a.solve(w,LossModel::Lossless,C64::new(1e-6,0.0),&gate).unwrap();
    let fb=b.solve(w,LossModel::Lossless,C64::new(1e-6,0.0),&gate).unwrap();
    let scale=b.gas().density/a.gas().density;
    assert_ne!(a.boundary_hash(),b.boundary_hash());
    for (&x,&y) in fa.boundary.z.iter().zip(&fb.boundary.z) {close(y,x.scale(scale));}
    close(fb.inlet.pressure,fa.inlet.pressure.scale(scale));
    for (&x,&y) in fa.outward_flows_m3_s.iter().zip(&fb.outward_flows_m3_s) {close(y,x);}
    for (x,y) in fa.receiver_fields.iter().zip(&fb.receiver_fields) {close(y.pressure_pa,x.pressure_pa.scale(scale));}
    let hot=GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(),350.0,101_325.0).unwrap();
    let c=model(hot,receivers()).solve(w,LossModel::Bessel,C64::new(1e-6,0.0),&gate).unwrap();
    assert!((c.receiver_fields[0].pressure_pa-fa.receiver_fields[0].pressure_pa).abs()>1e-8);
}

#[test]
fn observer_moves_without_changing_the_loading_and_cancelled_solve_is_retryable() {
    let gate=CancelGate::new_clock_free(); let w=core::f64::consts::TAU*1200.0;
    let a=model(gas(),receivers()); let b=model(gas(),vec![[0.2,0.11,0.17]]);
    assert_eq!(a.boundary_hash(),b.boundary_hash());
    let stopped=CancelGate::new_clock_free();stopped.request();
    assert!(matches!(a.solve(w,LossModel::Lossless,C64::ONE,&stopped),Err(ExteriorNetworkError::Cancelled)));
    let x=a.solve(w,LossModel::Lossless,C64::new(1e-6,0.0),&gate).unwrap();
    let y=b.solve(w,LossModel::Lossless,C64::new(1e-6,0.0),&gate).unwrap();
    assert_eq!(x.boundary.z,y.boundary.z);assert_eq!(x.response.impedance,y.response.impedance);
    assert_eq!(x.outward_flows_m3_s,y.outward_flows_m3_s);
    assert!((x.receiver_fields[0].pressure_pa-y.receiver_fields[0].pressure_pa).abs()>1e-12);
}

#[test]
fn ambiguous_geometry_incoherent_gas_and_budget_overflow_refuse() {
    let gate=CancelGate::new_clock_free();
    for case in 0..7 {
        let (surface,mut patches)=exterior();let mut gas=gas();let mut limit=budget();
        let mut points=receivers();let mut source=ContentHash([7;32]);
        match case {
            0=>patches[1].panels=patches[0].panels.clone(),
            1=>{let duplicate=patches[0].panels[0];patches[0].panels.push(duplicate);},
            2=>patches[0].panels=vec![surface.areas().len()],
            3=>gas.density*=2.0,
            4=>limit.max_field_entries=1,
            5=>points=vec![[0.0;3]],
            _=>source=ContentHash([0;32]),
        }
        let driver=NetworkDrivingPoint::new_with_coupled_ports(graph(),"inlet".into(),vec![],
            vec!["right".into(),"left".into()],schema()).unwrap();
        assert!(RadiatingNetwork::new(driver,surface,patches,gas,points,0.001,source,limit,&gate).is_err(),"accepted case {case}");
    }
}
