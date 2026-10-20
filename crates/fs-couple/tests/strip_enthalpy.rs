//! Original material, compatible section, thermal network and nonlinear owners.
//! Synthetic fixtures check equations, not calibrated material fidelity.
#[path = "strip_enthalpy/common.rs"]
mod common;
use common::*;
use std::collections::BTreeMap;
use fs_couple::strip_enthalpy::*;
use fs_solid::deposition::{AxialControl, BendingControl, Control};
use fs_thermal::enthalpy::*;

fn stretched(strain: f64) -> Control {
    Control {axial:AxialControl::Strain(strain),bending:BendingControl::CurvaturePerM(0.0)}
}

#[test]
fn layer_accounts_preserve_the_identical_section_solve_and_actual_heat_partition() {
    let materials = [material(5.0,3e8,1.0,0.0),material(10.0,1e8,1.0,0.0)];
    let (_,strip) = initial(&materials,0.0,vec![]);
    let temps = BTreeMap::from([(LAYERS[0],301.0),(LAYERS[1],303.0)]);
    let (a,r) = strip.advance_to(0.2,&temps,stretched(0.01),config().mechanics,||Ok(())).unwrap();
    let (b,s,layers) = strip.advance_with_layer_energy(0.2,&temps,stretched(0.01),config().mechanics,||Ok(())).unwrap();
    assert_eq!(a,b);
    assert_eq!(r.axial_force_n.to_bits(),s.axial_force_n.to_bits());
    assert_eq!(r.energy_j.dissipated.to_bits(),s.energy_j.dissipated.to_bits());
    assert_eq!(s.material_point_updates,4);
    let q = 0.01/1.2;
    for (i,modulus) in [3e8,1e8].into_iter().enumerate() {
        near(layers[&LAYERS[i]].dissipated,VOLUME*modulus*0.2*q*q,1e-14);
    }
    near(layers.values().map(|e|e.dissipated).sum(),r.energy_j.dissipated,1e-14);
}

#[test]
fn actual_layer_heat_enters_a_conservative_two_cell_transport_solve() {
    let materials = [material(5.0,3e8,1.0,0.0),material(10.0,1e8,1.0,0.0)];
    let (heat,strip) = initial(&materials,0.002,vec![]);
    let result = step_candidate(&heat,&strip,&bindings(&materials),0.2,stretched(0.01),config(),
        &[0.0;2],coefficients,||false).unwrap();
    // Independent constant-shift Maxwell heat and 2x2 backward-Euler network.
    let q = 0.01/1.2;
    let qa = VOLUME*3e8*0.2*q*q; let qb = VOLUME*1e8*0.2*q*q;
    let ca = MASS*5.0; let cb = MASS*10.0; let gdt = 0.002*0.2;
    let determinant = (ca+gdt)*(cb+gdt)-gdt*gdt;
    let ta = ((cb+gdt)*qa+gdt*qb)/determinant;
    let tb = (gdt*qa+(ca+gdt)*qb)/determinant;
    near(result.thermal().temperature_k(CELLS[0]).unwrap(),300.0+ta,1e-6);
    near(result.thermal().temperature_k(CELLS[1]).unwrap(),300.0+tb,1e-6);
    near(result.report().deposited_heat_j,qa+qb,1e-12);
    near(result.report().non_mechanical_heat_j,0.0,1e-14);
    assert!(result.report().mechanical.energy_j.numerical_dissipation>0.0);
    assert!(result.report().total_energy_residual_j.abs()<config().total_energy_tolerance_j);
    assert_eq!(heat.time_s(),0.0); assert_eq!(strip.axial_strain(),0.0);
}

#[test]
fn a_caloric_only_substitution_changes_compatible_curvature_and_both_heat_histories() {
    let run = |upper_cp| {
        let materials = [material(5.0,3e8,0.05,0.0),material(upper_cp,3e8,0.05,0.0)];
        let (heat,strip) = initial(&materials,0.0002,vec![]);
        let control = Control {axial:AxialControl::Strain(0.01),bending:BendingControl::MomentNm(0.0)};
        let next = step_candidate(&heat,&strip,&bindings(&materials),0.2,control,config(),
            &[0.0;2],coefficients,||false).unwrap();
        assert!(next.report().mechanical.moment_residual_nm.abs()<config().mechanics.moment_tolerance_nm);
        let endpoint = BTreeMap::from([
            (LAYERS[0],next.thermal().temperature_k(CELLS[0]).unwrap()),
            (LAYERS[1],next.thermal().temperature_k(CELLS[1]).unwrap())]);
        let independently_replayed = strip.advance_to(0.2,&endpoint,control,config().mechanics,||Ok(())).unwrap().0;
        assert_eq!(next.strip(),&independently_replayed);
        (endpoint[&LAYERS[0]],endpoint[&LAYERS[1]],next.strip().curvature_per_m())
    };
    let same = run(5.0); let changed = run(50.0);
    near(same.2,0.0,1e-7);
    assert!(changed.0>changed.1+0.1 && changed.2.abs()>1e-4);
    assert!((changed.0-same.0).abs()>1e-5,"compatibility feeds back into the unchanged lower material too");
}

#[test]
fn lossless_material_does_not_thermalize_numerical_loss() {
    let materials = [material(5.0,0.0,1.0,0.0),material(10.0,0.0,1.0,0.0)];
    let (heat,strip) = initial(&materials,0.002,vec![]);
    let next = step_candidate(&heat,&strip,&bindings(&materials),0.2,stretched(0.01),config(),
        &[0.0;2],coefficients,||false).unwrap();
    assert_eq!(next.report().deposited_heat_j,0.0);
    assert_eq!(next.thermal().energy_j(),heat.energy_j());
    assert!(next.report().mechanical.energy_j.numerical_dissipation>0.0);
    assert!(next.report().total_energy_residual_j.abs()<config().total_energy_tolerance_j);
}

#[test]
fn unloading_and_external_cooling_retain_layer_memory_and_close_imposed_work() {
    let materials = [material(5.0,3e8,0.05,0.0),material(10.0,3e8,0.05,0.0)];
    let (mut heat,mut strip) = initial(&materials,0.002,vec![]);
    let (mut work,mut other_heat,mut numerical) = (0.0,0.0,0.0);
    for i in 1..=8 {
        let control = Control {axial:AxialControl::ForceN(if i<=4 {40.0}else{0.0}),
            bending:BendingControl::MomentNm(0.0)};
        let solve = || step_candidate(&heat,&strip,&bindings(&materials),0.2*i as f64,control,config(),
            &[0.0;2],|s| {
                let mut c = coefficients(s).unwrap();
                c.sources_w.insert(CELLS[0],-0.0001);
                Ok::<_,core::convert::Infallible>(c)
            },||false).unwrap();
        let next = solve(); let replay = solve();
        assert_eq!(next.strip(),replay.strip()); assert_eq!(next.thermal().cells(),replay.thermal().cells());
        let r = next.report();
        near(r.non_mechanical_heat_j,-0.00002,1e-12);
        if i>4 { assert_eq!(r.external_mechanical_work_j,0.0); }
        work += r.external_mechanical_work_j; other_heat += r.non_mechanical_heat_j;
        numerical += r.mechanical.energy_j.numerical_dissipation;
        near(next.strip().stored_energy_j()+next.thermal().energy_j()+numerical-work-other_heat,0.0,
            i as f64*config().total_energy_tolerance_j);
        let (t,m,_) = next.into_parts(); heat=t; strip=m;
    }
    assert!(strip.layers()[&LAYERS[0]].material_points()[0].branch_strains()[0].abs()>1e-8);
    assert!(strip.axial_strain()>0.0);
}

#[test]
fn mappings_caloric_physics_and_the_final_trial_budget_are_not_silently_bypassed() {
    let materials = [material(5.0,3e8,0.05,0.0),material(10.0,3e8,0.05,0.0)];
    let (heat,strip) = initial(&materials,0.002,vec![]);
    let before = heat.cells().clone(); let old = strip.clone(); let map = bindings(&materials);
    let bad = [map[0],LayerThermalBinding {cell:map[0].cell,..map[1]}];
    for mapping in [&map[..1],&bad[..]] {
        assert!(step_candidate(&heat,&strip,mapping,0.2,stretched(0.01),config(),&[0.0;2],coefficients,||false).is_err());
    }
    let foreign = [map[0],LayerThermalBinding {caloric:map[0].caloric,..map[1]}];
    assert!(step_candidate(&heat,&strip,&foreign,0.2,stretched(0.01),config(),&[0.0;2],coefficients,||false).is_err());
    let normal = step_candidate(&heat,&strip,&map,0.2,stretched(0.01),config(),&[0.0;2],coefficients,||false).unwrap();
    assert!(matches!(step_candidate(&heat,&strip,&map,0.2,stretched(0.01),
        StripEnthalpyConfig {maximum_structural_trials:normal.report().structural_trials-1,..config()},
        &[0.0;2],coefficients,||false),Err(StripEnthalpyError::StructuralBudget)));
    let flag = std::cell::Cell::new(false); let calls = std::cell::Cell::new(0usize);
    let late = step_candidate(&heat,&strip,&map,0.2,stretched(0.01),config(),&[0.0;2],|s| {
        calls.set(calls.get()+1);
        if calls.get()==normal.report().structural_trials {flag.set(true);}
        coefficients(s)
    },||flag.get());
    assert!(matches!(late,Err(StripEnthalpyError::Cancelled)));
    assert_eq!(heat.cells(),&before); assert_eq!(strip,old);
    let unsupported = [material(5.0,3e8,0.05,0.001),material(10.0,3e8,0.05,0.0)];
    let (h,s) = initial(&unsupported,0.0,vec![]);
    assert!(matches!(step_candidate(&h,&s,&bindings(&unsupported),0.2,stretched(0.01),config(),
        &[0.0;2],coefficients,||false),Err(StripEnthalpyError::CaloricAdmission(_))));
    assert!(matches!(step_candidate(&heat,&strip,&map,0.2,stretched(0.01),config(),&[0.0;2],
        |_|Err::<ThermalCoefficients<()>,_>("missing thermal source"),||false),
        Err(StripEnthalpyError::Source("missing thermal source"))));
}

#[test]
fn a_later_layer_is_born_once_on_the_solved_shape_after_endpoint_acceptance() {
    let materials = [material(5.0,3e8,0.05,0.0),material(10.0,3e8,0.05,0.0)];
    let mut born = cell(CellId(99),2); born.material_id=MaterialId(0);
    let event = TimedEvent {time_s:0.4,event:Event::Birth(born)};
    let (heat,strip) = initial(&materials,0.002,vec![event]);
    let map = bindings(&materials);
    assert!(step_candidate(&heat,&strip,&map,0.5,stretched(0.01),config(),&[0.0;2],coefficients,||false).is_err());
    let left = step_candidate(&heat,&strip,&map,0.4,stretched(0.01),config(),&[0.0;2],coefficients,||false).unwrap();
    assert_eq!(left.thermal().cells().len(),2);
    let (right, birth_audit) = left.thermal().advance_to(0.4,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap();
    near(birth_audit.incoming_mass_kg,MASS,1e-16);
    let next_strip = left.strip().deposit(99,geometry(2),materials[0].mechanical.law(),
        right.temperature_k(CellId(99)).unwrap()).unwrap();
    for id in LAYERS { assert_eq!(next_strip.layers()[&id],left.strip().layers()[&id]); }
    for p in next_strip.layers()[&99].material_points() { assert_eq!(p.stress_pa(),0.0); }
    let complete = [map[0],map[1],LayerThermalBinding {layer:99,cell:CellId(99),..map[0]}];
    let next = step_candidate(&right,&next_strip,&complete,0.6,stretched(0.015),config(),
        &[0.0;3],coefficients,||false).unwrap();
    assert_eq!(next.strip().layers().len(),3);
    assert_eq!(next.report().mechanical.material_point_updates,6);
    assert_eq!(next.strip().layers()[&99].birth_time_s,0.4);
    assert!(next.report().layer_energy_j[&99].dissipated>0.0);
    assert_eq!(heat.cells().len(),2); assert_eq!(strip.layers().len(),2);
}