//! Spatial heat and sag use the original Maxwell, beam and enthalpy owners.
#[path = "strip_enthalpy/common.rs"]
mod common;
use common::*;
use std::collections::BTreeMap;
use fs_couple::strip_enthalpy::{LayerThermalBinding, StripEnthalpyError};
use fs_couple::strip_enthalpy::bridge::*;
use fs_couple::circuit_enthalpy::material::{ConductionLink,MaterialCell,MaterialNetwork};
use fs_material::state_point::{MaterialPropertySelection,THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{PropertyKey,QueryPoint};
use fs_solid::bridge::{self as solid,Bridge,Support};
use fs_thermal::enthalpy::*;

const AXIAL:[u64;2]=[20,3];
const HEAT:[[CellId;2];2]=[[CellId(90),CellId(4)],[CellId(50),CellId(12)]];
const LINKS:[ConductionLink;4]=[
    ConductionLink {link:5,distance_a_m:0.0005,distance_b_m:0.0005},
    ConductionLink {link:6,distance_a_m:0.0005,distance_b_m:0.0005},
    ConductionLink {link:7,distance_a_m:0.025,distance_b_m:0.025},
    ConductionLink {link:8,distance_a_m:0.025,distance_b_m:0.025},
];
fn solve_config()->BridgeEnthalpyConfig {
    BridgeEnthalpyConfig {feedback:common::config(),maximum_deflection_m:0.009}
}
fn initial_bridge(m:&[Material;2],mass_scale:f64)->(DepositionSimulation<'_>,Bridge<'_>) {
    let mut events=Vec::new();
    for (i,ids) in HEAT.iter().enumerate() {
        for (j,&id) in ids.iter().enumerate() {
            let z=geometry(j).centroid_z_m;
            events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
                id,segment_id:AXIAL[i],body_id:1,material_id:MaterialId(j as u64),
                mass_kg:MASS*0.5*mass_scale,specific_h_j_kg:0.0,
                geometry:CellGeometry {start_m:[0.05*i as f64,0.0,z],
                    end_m:[0.05*(i+1) as f64,0.0,z],volume_m3:VOLUME*0.5,
                    surface_area_m2:0.00112},boundary:boundary()})});
        }
    }
    for (id,a,b,area) in [(5,HEAT[0][0],HEAT[0][1],0.0005),
        (6,HEAT[1][0],HEAT[1][1],0.0005),(7,HEAT[0][0],HEAT[1][0],AREA),
        (8,HEAT[0][1],HEAT[1][1],AREA)] {
        events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
            id,a,b,area_m2:area,conductance_w_k:9.0})});
    }
    let thermal=DepositionSimulation::new(events,&[
        (MaterialId(0),m[0].caloric.curve() as &dyn ThermalLaw),
        (MaterialId(1),m[1].caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:8,max_links:8,max_events:16,max_substeps:1000}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
    let mut bridge=Bridge::new(&[(AXIAL[0],0.05),(AXIAL[1],0.05)],Support::CantileverLeft,
        solid::Budget {max_cells:2,max_layers_per_cell:2,max_branch_states:8}).unwrap();
    for j in 0..2 {
        bridge=bridge.deposit_bonded_material_layer(LAYERS[j],geometry(j),&m[j].mechanical,&m[j].caloric,
            &BTreeMap::from([(AXIAL[0],300.0),(AXIAL[1],300.0)]),&mut ||false).unwrap();
    }
    (thermal,bridge)
}
fn map(m:&[Material;2])->Vec<BridgeThermalBinding<'_>> {
    let mut map=Vec::new();
    for i in 0..2 { for j in 0..2 {
        map.push(BridgeThermalBinding {axial_cell:AXIAL[i],section:LayerThermalBinding {
            layer:LAYERS[j],cell:HEAT[i][j],material:&m[j].mechanical,caloric:&m[j].caloric}});
    }}
    map
}
fn transport_cells<'a>(m:&'a [Material;2],context:&'a QueryPoint,key:&'a PropertyKey)->Vec<MaterialCell<'a>> {
    let mut result=Vec::new();
    for row in HEAT { for j in 0..2 {
        result.push(MaterialCell {cell:row[j],material_id:MaterialId(j as u64),curve:m[j].caloric.curve(),
            point:ThermalMaterialPoint {material:&m[j].card,context,temperature_axis:"T",liquid_fraction_axis:None},
            conductivity:Some((key,MaterialPropertySelection::SingleClaimOnly))});
    }}
    result
}
fn materials(k:f64,branch:f64)->[Material;2] {
    [0.05,0.5].map(|cp|material_with_conductivity(cp,branch,0.05,0.0,vec![(280.0,k),(320.0,k)]))
}
fn key()->PropertyKey {PropertyKey::new("conductivity",THERMAL_CONDUCTIVITY_DIMS)}

#[test]
fn spatial_heating_uses_actual_cell_work_and_replays_the_original_bridge() {
    let m=materials(0.0002,3e8); let (heat,bridge)=initial_bridge(&m,1.0);
    let context=point(300.0); let key=key(); let cells=transport_cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&cells,resistors:&[],links:&LINKS};
    let prepared=network.prepare(&heat,||false).unwrap(); let bindings=map(&m);
    let run=||step_transport_candidate(&heat,&bridge,&bindings,&prepared,0.2,9.81,&[1.0;2],solve_config(),&[0.0;4],||false).unwrap();
    let next=run(); let replay=run();
    assert_eq!(next.bridge(),replay.bridge()); assert_eq!(next.thermal().cells(),replay.thermal().cells());
    let temperatures=bindings.iter().map(|b|((b.axial_cell,b.section.layer),
        next.thermal().temperature_k(b.section.cell).unwrap())).collect();
    let manual=bridge.advance_controlled(0.2,&temperatures,9.81,&[1.0;2],
        solid::Config {section:solve_config().feedback.mechanics,maximum_deflection_m:0.009},&mut ||false).unwrap();
    assert_eq!(next.bridge(),&manual.0);
    let r=next.report();
    near(r.deposited_heat_j,r.layer_energy_j.values().map(|e|e.dissipated).sum(),1e-12);
    assert!(r.layer_energy_j[&(AXIAL[0],LAYERS[0])].dissipated
        >10.0*r.layer_energy_j[&(AXIAL[1],LAYERS[0])].dissipated);
    assert!(next.thermal().temperature_k(HEAT[0][0]).unwrap()>next.thermal().temperature_k(HEAT[1][0]).unwrap());
    assert_eq!(r.mechanical.material_point_updates,8);
    assert_eq!(r.source.cells.len(),4); assert_eq!(r.source.links.len(),4);
    near(r.mechanical.vertical_reactions_n[0],bridge.mass_kg()*9.81+0.1,1e-12);
    assert!(r.total_energy_residual_j.abs()<solve_config().feedback.total_energy_tolerance_j);
    assert_eq!(bridge.time_s(),0.0); assert_eq!(heat.time_s(),0.0);
}

#[test]
fn conductivity_only_substitution_changes_solved_temperatures_and_sag() {
    let run=|k| {
        let m=materials(k,3e8); let (heat,bridge)=initial_bridge(&m,1.0);
        let context=point(300.0); let key=key(); let cells=transport_cells(&m,&context,&key);
        let network=MaterialNetwork {cells:&cells,resistors:&[],links:&LINKS};
        let prepared=network.prepare(&heat,||false).unwrap();
        let next=step_transport_candidate(&heat,&bridge,&map(&m),&prepared,0.2,9.81,&[1.0;2],solve_config(),&[0.0;4],||false).unwrap();
        (next.thermal().temperature_k(HEAT[0][0]).unwrap()-next.thermal().temperature_k(HEAT[0][1]).unwrap(),
            next.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1],next.bridge().mass_kg())
    };
    let low=run(0.0002); let high=run(0.02);
    assert!(low.0>high.0+0.01);
    assert!((low.1-high.1).abs()>1e-9);
    assert_eq!(low.2,high.2,"a conductivity substitution does not change gravity's reference mass");
}

#[test]
fn numerical_loss_does_not_heat_a_lossless_beam() {
    let m=materials(0.002,0.0); let (heat,bridge)=initial_bridge(&m,1.0);
    // Reduce load for the softer equilibrium-only specimen, retaining gravity work.
    let next=step_candidate(&heat,&bridge,&map(&m),0.2,1.0,&[0.1;2],solve_config(),&[0.0;4],coefficients,||false).unwrap();
    assert_eq!(next.report().deposited_heat_j,0.0);
    assert_eq!(next.thermal().cells(),heat.cells());
    assert!(next.report().mechanical.energy_j.numerical_dissipation>0.0);
    assert!(next.report().total_energy_residual_j.abs()<solve_config().feedback.total_energy_tolerance_j);
}

#[test]
fn affine_spring_material_cools_original_spatial_cells_and_changes_sag_without_dashpots() {
    let slope=-0.005;
    let m=[0.05,0.5].map(|cp|material_with_affine_scale(cp,0.0,0.05,0.0,
        vec![(280.0,0.002),(320.0,0.002)],None,vec![(280.0,0.1),(320.0,0.1)],slope));
    let (heat,bridge)=initial_bridge(&m,1.0);
    let next=step_candidate(&heat,&bridge,&map(&m),0.2,1.0,&[0.1;2],solve_config(),&[0.0;4],coefficients,||false).unwrap();
    let r=next.report();
    assert_eq!(r.mechanical.energy_j.dissipated,0.0);
    assert_eq!(r.mechanical.energy_j.modulus_work,0.0,"initial springs are stress-free");
    assert!(r.mechanical.energy_j.numerical_dissipation>0.0);
    assert!(r.reversible_heat_j<0.0 && r.deposited_heat_j<0.0);
    near(r.deposited_heat_j,r.reversible_heat_j,1e-12);
    near(next.thermal().energy_j()-heat.energy_j(),r.deposited_heat_j,1e-9);
    let mut expected_u=0.0;let mut expected_ts=0.0;
    for (i,cell) in next.bridge().cells().iter().enumerate() {
        for (j,&layer) in LAYERS.iter().enumerate() {
            let t=next.thermal().temperature_k(HEAT[i][j]).unwrap();
            assert!(t<300.0);
            // Independent oracle for the authored S(T), evaluated at the actual
            // original Gauss histories. This fixture starts with U=s=0.
            for point in cell.section().layers()[&layer].material_points() {
                let xi=point.stored_energy_j_m3()/(1.0+slope*(t-300.0));
                let weight=0.5*AREA*cell.length_m();
                expected_u+=weight*(1.0-300.0*slope)*xi;
                expected_ts+=weight*(-slope)*t*xi;
            }
        }
    }
    near(r.mechanical_internal_energy_change_j,expected_u,1e-12);
    near(r.mechanical_entropy_energy_change_j,expected_ts,1e-12);
    near(r.reversible_heat_j+expected_ts,0.0,1e-12);
    near(expected_u+next.thermal().energy_j()+r.mechanical.energy_j.numerical_dissipation
        -r.mechanical.load_work_j,0.0,solve_config().feedback.total_energy_tolerance_j);
    let constant=materials(0.002,0.0);let (old_heat,old_beam)=initial_bridge(&constant,1.0);
    let baseline=step_candidate(&old_heat,&old_beam,&map(&constant),0.2,1.0,&[0.1;2],solve_config(),&[0.0;4],coefficients,||false).unwrap();
    assert!((next.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1]
        -baseline.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1]).abs()>1e-9);
    assert_eq!(next.bridge().mass_kg(),bridge.mass_kg());
    assert_eq!(next.bridge().mass_kg(),baseline.bridge().mass_kg());
    assert_eq!(heat.time_s(),0.0);assert_eq!(bridge.time_s(),0.0);
}

#[test]
fn source_expansion_bends_the_actual_beam_and_retains_thermal_maxwell_history() {
    let alpha: f64=1e-4;
    let cp=5.0;
    let dt=0.2;
    let a=1.0/(1.0+dt);
    let mut final_tip=[0.0;2];
    for (variant,upper_sign) in [1.0,-1.0].into_iter().enumerate() {
        let alphas=[alpha,upper_sign*alpha];
        let m=alphas.map(|value|material(cp,3e8,1.0,40.0*value));
        let (mut heat,mut bridge)=initial_bridge(&m,1.0);
        let mut old_rise=0.0;
        let mut old_q=0.0;
        let mut old_internal=0.0;
        let mut old_ts=0.0;
        let mut accumulated_heat=0.0;
        let mut numerical=0.0;
        for (step,heat_per_layer) in [0.01,-0.004].into_iter().enumerate() {
            let next=step_candidate(&heat,&bridge,&map(&m),(step+1) as f64*dt,0.0,&[0.0;2],
                solve_config(),&[0.0;4],|s|coefficients(s).map(|mut c| {
                    for row in HEAT { for cell in row {
                        c.sources_w.insert(cell,heat_per_layer/(2.0*dt));
                    }}
                    c
                }),||false).unwrap();
            // Equal heating and mirror-related layers leave all four thermal
            // cells isothermal, but opposite source expansion bends the bonded
            // beam. Its finite layer thickness leaves real residual stress.
            // With E_inf=1e8, E_branch=3e8 and tau=1, the native Gauss strains
            // factor as e=c*(T-300), q=c*q_scalar. The caloric equation is a
            // scalar quadratic; use its small root without cancellation.
            let old_temperature=300.0+old_rise;
            let factor=if upper_sign<0.0 {0.25*VOLUME*alpha*alpha}else{0.0};
            let quadratic=factor*3e8*dt*a*a;
            let linear=MASS*cp-factor*old_temperature*(1e8+3e8*a)-2.0*quadratic*old_q;
            let constant=heat_per_layer+quadratic*old_q*old_q
                -factor*old_temperature*3e8*(1.0-a)*old_q;
            let increment=2.0*constant/(linear+(linear*linear-4.0*quadratic*constant).sqrt());
            let rise=old_rise+increment;
            let temperature=300.0+rise;
            let q_scalar=a*(old_q+increment);
            let strain=0.5*(alphas[0]+alphas[1])*rise;
            // Integral of z² over both native rectangular layers, including
            // their thickness moments, gives kappa=750*(alpha0-alpha1)*ΔT.
            let curvature=750.0*(alphas[0]-alphas[1])*rise;
            let mut internal=0.0;
            let mut ts=0.0;
            for (i,cell) in next.bridge().cells().iter().enumerate() {
                near(cell.section().axial_strain(),strain,1e-11);
                near(cell.section().curvature_per_m(),curvature,1e-8);
                for (j,&layer) in LAYERS.iter().enumerate() {
                    near(next.thermal().temperature_k(HEAT[i][j]).unwrap(),temperature,1e-6);
                    let shape=geometry(j);
                    let weight=0.5*AREA*cell.length_m();
                    for (point,sign) in cell.section().layers()[&layer].material_points().iter().zip([-1.0,1.0]) {
                        let z=shape.centroid_z_m+sign*shape.thickness_m/(2.0*3.0_f64.sqrt());
                        let coefficient=0.5*(alphas[0]+alphas[1])
                            -z*750.0*(alphas[0]-alphas[1])-alphas[j];
                        let e=coefficient*rise;
                        let q=coefficient*q_scalar;
                        let stress=1e8*e+3e8*q;
                        near(point.total_strain(),strain-z*curvature,1e-11);
                        near(point.branch_strains()[0],q,1e-11);
                        near(point.stress_pa(),stress,1e-3);
                        let entropy_energy=weight*temperature*alphas[j]*stress;
                        ts+=entropy_energy;
                        internal+=weight*0.5*(1e8*e*e+3e8*q*q)+entropy_energy;
                    }
                }
            }
            let r=next.report();
            near(next.bridge().deformed_position_m(LENGTH,0.0).unwrap()[0],LENGTH*(1.0+strain),1e-11);
            near(next.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1],0.5*curvature*LENGTH*LENGTH,1e-10);
            near(r.mechanical_internal_energy_change_j,internal-old_internal,1e-9);
            near(r.mechanical_entropy_energy_change_j,ts-old_ts,1e-9);
            near(r.reversible_heat_j,-r.mechanical.energy_j.thermal_strain_work-(ts-old_ts),1e-9);
            near(r.deposited_heat_j,r.mechanical.energy_j.dissipated+r.reversible_heat_j,1e-10);
            accumulated_heat+=2.0*heat_per_layer;
            numerical+=r.mechanical.energy_j.numerical_dissipation;
            near(internal+next.thermal().energy_j()+numerical-accumulated_heat,0.0,4e-9);
            assert_eq!(r.mechanical.load_work_j,0.0);
            assert_eq!(next.bridge().mass_kg(),bridge.mass_kg());
            if upper_sign<0.0 {
                assert!(r.mechanical.energy_j.dissipated>0.0);
                assert!(r.mechanical.energy_j.thermal_strain_work.abs()>1e-8);
                assert!(if step==0 {r.reversible_heat_j>0.0}else{r.reversible_heat_j<0.0});
            }
            old_rise=rise;old_q=q_scalar;old_internal=internal;old_ts=ts;
            let (h,b,_)=next.into_parts();heat=h;bridge=b;
        }
        final_tip[variant]=bridge.deformed_position_m(LENGTH,0.0).unwrap()[1];
    }
    near(final_tip[0],0.0,1e-10);
    assert!(final_tip[1]>5e-4,"changing only upper-layer expansion changes the actual heated shape");
}

#[test]
fn unloading_retains_spatial_memory_and_closes_the_whole_history() {
    let m=materials(0.002,3e8); let (heat,bridge)=initial_bridge(&m,1.0);
    let context=point(300.0); let key=key(); let cells=transport_cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&cells,resistors:&[],links:&LINKS};
    let prepared=network.prepare(&heat,||false).unwrap(); let bindings=map(&m);
    let first=step_transport_candidate(&heat,&bridge,&bindings,&prepared,0.2,9.81,&[1.0;2],solve_config(),&[0.0;4],||false).unwrap();
    let next=step_transport_candidate(first.thermal(),first.bridge(),&bindings,&prepared,0.4,0.0,&[0.0;2],solve_config(),&[0.0;4],||false).unwrap();
    assert_eq!(next.report().mechanical.load_work_j,0.0);
    assert!(next.report().deposited_heat_j>0.0);
    near(next.bridge().stored_energy_j()+next.thermal().energy_j()
        +first.report().mechanical.energy_j.numerical_dissipation
        +next.report().mechanical.energy_j.numerical_dissipation
        -first.report().mechanical.load_work_j,0.0,2.0*solve_config().feedback.total_energy_tolerance_j);
    assert!(next.bridge().cells()[0].section().layers()[&LAYERS[0]].material_points()[0]
        .branch_strains()[0].abs()>1e-10);
}

#[test]
fn mass_mapping_source_and_final_cancellation_fail_without_partial_publication() {
    let m=materials(0.002,3e8); let (heat,bridge)=initial_bridge(&m,1.0);
    let before=heat.cells().clone(); let old=bridge.clone(); let bindings=map(&m);
    let mut bad=bindings.clone(); bad[1].section.cell=bad[0].section.cell;
    assert!(step_candidate(&heat,&bridge,&bad,0.2,9.81,&[1.0;2],solve_config(),&[0.0;4],coefficients,||false).is_err());
    let (wrong_mass,_)=initial_bridge(&m,2.0);
    assert!(step_candidate(&wrong_mass,&bridge,&bindings,0.2,9.81,&[1.0;2],solve_config(),&[0.0;4],coefficients,||false).is_err());
    let context=point(300.0); let key=key(); let cells=transport_cells(&m,&context,&key);
    let network=MaterialNetwork {cells:&cells,resistors:&[],links:&LINKS};
    let prepared=network.prepare(&heat,||false).unwrap();
    let mut polls=0;
    step_transport_candidate(&heat,&bridge,&bindings,&prepared,0.2,9.81,&[1.0;2],solve_config(),&[0.0;4],||{polls+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(step_transport_candidate(&heat,&bridge,&bindings,&prepared,0.2,9.81,&[1.0;2],solve_config(),&[0.0;4],
        ||{actual+=1;actual==polls}),Err(BridgeEnthalpyError::Coupled(StripEnthalpyError::Cancelled))));
    let altered=point(300.0).with("lot",17.0).unwrap(); let cells=transport_cells(&m,&altered,&key);
    let other=MaterialNetwork {cells:&cells,resistors:&[],links:&LINKS};
    let prepared=other.prepare(&heat,||false).unwrap();
    assert!(step_transport_candidate(&heat,&bridge,&bindings,&prepared,0.2,9.81,&[1.0;2],solve_config(),&[0.0;4],||false).is_err());
    assert_eq!(heat.cells(),&before); assert_eq!(bridge,old);
}
