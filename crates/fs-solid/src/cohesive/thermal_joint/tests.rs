use super::*;
use fs_material::thermovisco::{Limits,ThermalKnot};
use fs_solver::NewtonKrylovConfig;

fn material(branches: Vec<(f64,f64)>) -> ThermoMaxwell {
    ThermoMaxwell::new(1e9,branches,[280.0,520.0].into_iter().map(|t| ThermalKnot {
        temperature_k: t,modulus_scale: 1.0,relaxation_scale: 1.0,free_strain: (t-300.0)*1e-4,
    }).collect(),Limits { max_strain: 0.04,max_strain_rate_s: 1.0 }).unwrap()
}
fn bond() -> BilinearCohesive { BilinearCohesive::new(1e10,2e5,10.0,2e10,0.01).unwrap() }
fn budget() -> Budget { Budget { max_members: 8,max_branch_states: 32,max_krylov_basis_entries: 64 } }
fn config() -> Config {
    Config { equilibrium: super::super::Config {
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-10,relative_tolerance: 1e-12,
            linear_restart: 4,forcing_maximum: 0.01,..NewtonKrylovConfig::default() },
        max_newton_iterations: 80,max_operator_evaluations: 2000,force_tolerance_n: 1e-7,
        local_energy_tolerance_j_m2: 1e-7,total_cohesive_energy_tolerance_j: 1e-9,
        maximum_displacement_m: 0.005 },extension_tolerance_m: 1e-11,energy_tolerance_j: 1e-7 }
}
fn joint<'a>(material: &'a ThermoMaxwell,bond: &'a BilinearCohesive) -> Joint<'a> {
    Joint::new(0.0,vec![Member::new(1,0.01,1e-4,material,500.0).unwrap()],7,1e-4,bond,budget()).unwrap()
}
fn near(a: f64,b: f64,tol: f64) { assert!((a-b).abs()<tol,"{a:e} != {b:e}"); }

// G1: exact series compliance and thermal eigenstrain, before cohesive softening.
#[test]
fn constrained_cooling_solves_opening_and_reaction_from_temperature() {
    let material=material(vec![]); let law=bond(); let initial=joint(&material,&law);
    let (next,r)=initial.advance_controlled(1.0,0.0,&BTreeMap::from([(1,480.0)]),config(),&mut || false).unwrap();
    let gap=2e-5*1e7/(1e7+1e6);
    near(r.opening_m,gap,1e-13); near(r.force_n,1e6*gap,1e-7);
    near(next.members()[&1].point().total_strain(),-gap/0.01,1e-12);
    near(r.endpoint_work_j,0.0,1e-15);
    assert!(r.bulk_energy_j.thermal_strain_work>0.0);
    assert_eq!(next.bond().damage(),0.0);
    assert!(r.energy_residual_j.abs()<config().energy_tolerance_j);
    near(next.extension_m(),0.0,1e-15);
}

#[test]
fn cooling_can_break_the_joint_and_reheating_does_not_heal_it() {
    let material=material(vec![]); let law=bond(); let initial=joint(&material,&law);
    let (broken,r)=initial.advance_controlled(1.0,0.0,&BTreeMap::from([(1,300.0)]),config(),&mut || false).unwrap();
    assert!(broken.bond().failed()); near(r.opening_m,0.0002,1e-12); near(r.force_n,0.0,1e-7);
    near(r.cohesive_energy_j.dissipated,0.001,1e-10);
    let (hot,r)=broken.advance_controlled(2.0,0.0,&BTreeMap::from([(1,500.0)]),config(),&mut || false).unwrap();
    assert!(hot.bond().failed()); near(r.opening_m,0.0,1e-12);
    let (compressed,r)=hot.advance_controlled(3.0,-1e-5,&BTreeMap::from([(1,500.0)]),config(),&mut || false).unwrap();
    near(r.opening_m,-1e-5*1e7/(1e7+2e6),1e-12);
    assert!(r.force_n<0.0 && compressed.bond().failed());
    assert_eq!(r.cohesive_energy_j.dissipated,0.0);
}

// G1: a Maxwell branch in series with an elastic cohesive spring has an exact
// effective relaxation time. Compare the discrete BE recurrence, not raw tau.
#[test]
fn stress_relaxation_uses_member_history_and_joint_compliance() {
    let material=material(vec![(2e9,2.0)]); let law=bond(); let initial=joint(&material,&law);
    let temperatures=BTreeMap::from([(1,500.0)]);
    let (mut state,loaded)=initial.advance_controlled(0.1,1e-5,&temperatures,config(),&mut || false).unwrap();
    let q0=state.members()[&1].point().branch_strains()[0];
    let gamma: f64=1e-4/(1e6*0.01);
    let tau_eff=2.0*(1.0+gamma*3e9)/(1.0+gamma*1e9);
    let mut final_force=loaded.force_n;
    for step in 1..=20 {
        let (next,r)=state.advance_controlled(0.1+f64::from(step)*0.1,1e-5,&temperatures,config(),&mut || false).unwrap();
        assert!(r.bulk_energy_j.dissipated>0.0);
        state=next; final_force=r.force_n;
    }
    near(state.members()[&1].point().branch_strains()[0],q0*(1.0+0.1/tau_eff).powi(-20),1e-13);
    assert!(final_force<loaded.force_n); assert_eq!(state.bond().damage(),0.0);
}

// G3: axial partition changes neither represented compliance nor extensive work.
#[test]
fn member_partition_retains_force_opening_and_energy() {
    let material=material(vec![(2e9,2.0)]); let law=bond();
    let whole=joint(&material,&law);
    let split=Joint::new(0.0,vec![Member::new(2,0.004,1e-4,&material,500.0).unwrap(),
        Member::new(1,0.006,1e-4,&material,500.0).unwrap()],7,1e-4,&law,budget()).unwrap();
    let (a,ar)=whole.advance_controlled(1.0,0.0,&BTreeMap::from([(1,480.0)]),config(),&mut || false).unwrap();
    let (b,br)=split.advance_controlled(1.0,0.0,&BTreeMap::from([(1,480.0),(2,480.0)]),config(),&mut || false).unwrap();
    near(ar.force_n,br.force_n,1e-7); near(ar.opening_m,br.opening_m,1e-12);
    near(a.stored_energy_j(),b.stored_energy_j(),1e-10);
    near(ar.bulk_energy_j.dissipated,br.bulk_energy_j.dissipated,1e-10);
}

#[test]
fn failed_or_cancelled_proposals_preserve_both_material_histories() {
    let material=material(vec![(2e9,2.0)]); let law=bond(); let initial=joint(&material,&law);
    let saved=initial.clone(); let temperatures=BTreeMap::from([(1,480.0)]);
    for bad in [BTreeMap::new(),BTreeMap::from([(2,480.0)]),BTreeMap::from([(1,600.0)])] {
        assert!(initial.advance_controlled(1.0,0.0,&bad,config(),&mut || false).is_err());
    }
    assert!(initial.advance_controlled(0.0,0.0,&temperatures,config(),&mut || false).is_err());
    let mut limited=config(); limited.equilibrium.max_operator_evaluations=1;
    assert!(initial.advance_controlled(1.0,0.0,&temperatures,limited,&mut || false).is_err());
    for allowed in [0,2,10] {
        let mut calls=0;
        assert!(initial.advance_controlled(1.0,0.0,&temperatures,config(),&mut || { calls+=1; calls>allowed }).is_err());
    }
    assert_eq!(initial,saved);
}

#[test]
fn checkpoint_replay_and_nonuniform_member_temperatures_preserve_compatibility() {
    let material=material(vec![(2e9,2.0)]); let law=bond();
    let state=Joint::new(0.0,vec![Member::new(2,0.005,1e-4,&material,450.0).unwrap(),
        Member::new(1,0.005,2e-4,&material,500.0).unwrap()],7,1e-4,&law,budget()).unwrap();
    let t=BTreeMap::from([(1,480.0),(2,430.0)]);
    let (prefix,_)=state.advance_controlled(1.0,0.0,&t,config(),&mut || false).unwrap();
    let (a,ar)=prefix.advance_controlled(2.0,0.0,&t,config(),&mut || false).unwrap();
    let (b,br)=prefix.clone().advance_controlled(2.0,0.0,&t,config(),&mut || false).unwrap();
    assert_eq!(a,b); assert_eq!(ar.force_n,br.force_n);
    let extension=a.members().values().map(|m| m.reference_length_m()*m.point().total_strain()).sum::<f64>()
        +a.bond().jump_m()[0];
    near(extension,0.0,config().extension_tolerance_m);
    for member in a.members().values() { near(member.reference_area_m2()*member.point().stress_pa(),ar.force_n,1e-7); }
}

#[test]
fn thermal_joint_uses_real_interval_for_bond_viscosity_and_energy() {
    use fs_material::cohesive::CohesiveViscosity;
    let material=material(vec![]);let law=bond().with_viscosity(CohesiveViscosity::new(1e9,1e9,0.1).unwrap());
    let initial=joint(&material,&law);let temperatures=BTreeMap::from([(1,480.0)]);
    for dt in [0.1,1.0,10.0] {
        let (next,r)=initial.advance_controlled(dt,0.0,&temperatures,config(),&mut ||false).unwrap();
        let expected_gap=2e-5*1e7/(1e7+1e-4*(1e10+1e9/dt));
        near(r.opening_m,expected_gap,1e-13);
        near(r.force_n,1e-4*(1e10+1e9/dt)*expected_gap,1e-7);
        near(r.bond_viscous_dissipated_j,1e-4*1e9*expected_gap*expected_gap/dt,1e-12);
        near(r.cohesive_energy_j.dissipated,0.0,1e-14);
        assert_eq!(next.bond().damage(),0.0);
        assert!(r.energy_residual_j.abs()<config().energy_tolerance_j);
    }
}

#[test]
fn thermal_joint_fracture_accounts_viscous_path_loss_before_failure() {
    use fs_material::cohesive::CohesiveViscosity;
    let material=material(vec![(2e8,0.4)]);let law=bond().with_viscosity(CohesiveViscosity::new(1e8,1e8,0.1).unwrap());
    let initial=joint(&material,&law);
    let (failed,r)=initial.advance_controlled(1.0,0.0,&BTreeMap::from([(1,300.0)]),config(),&mut ||false).unwrap();
    assert!(failed.bond().failed());assert!(r.bond_viscous_dissipated_j>0.0);
    near(r.force_n,0.0,1e-7);near(r.cohesive_energy_j.dissipated,0.001,1e-10);
    assert!(r.energy_residual_j.abs()<config().energy_tolerance_j);
    let (reclosed,r)=failed.advance_controlled(2.0,-1e-5,&BTreeMap::from([(1,500.0)]),config(),&mut ||false).unwrap();
    assert!(reclosed.bond().failed());assert_eq!(r.bond_viscous_dissipated_j,0.0);
    assert_eq!(r.cohesive_energy_j.dissipated,0.0);
}
