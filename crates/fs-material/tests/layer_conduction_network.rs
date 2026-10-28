//! Original source resolution and actual implicit enthalpy history, not a second integrator.
use std::sync::Arc;
use fs_evidence::ValidityDomain;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::thermal::{source::ThermalMaterialPoint, conduction::{
    LayerConductionConfig,
    source::{LayerSourceError, ResolvedConductionLayer},
    network::{LayerLinkBinding, PreparedLayerConduction},
}};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_matdb::{ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{*, backward_euler::BackwardEulerStep};

fn config() -> LayerConductionConfig {
    LayerConductionConfig { max_iterations:96, absolute_heat_tolerance_w:1e-10, relative_heat_tolerance:1e-12 }
}
fn source() -> MaterialCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS),
        value:PropertyValue::Curve { abscissa:"T".into(), abscissa_dims:Temperature::DIMS,
            knots:vec![(200.0,4.0),(450.0,6.5),(1000.0,12.0)], dims:THERMAL_CONDUCTIVITY_DIMS },
        validity:ValidityDomain::unconstrained().with("T",200.0,1000.0),
        uncertainty:UncertaintyModel::Unstated, interpolation:InterpolationPolicy::LinearInside,
        observations:vec![], provenance:Provenance { source:"manufactured linear conductivity".into(),
            license:"CC0-1.0".into(), artifact:None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry:"synthetic".into(), phase:"solid".into(),
        process:"numerical fixture".into(),revision:0 },claims,vec![]).unwrap()
}
fn curve(card:&MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(card.content_hash(),SolidLiquidPhase::Solid,
        [(0.0,200.0),(80_000.0,1000.0)].into_iter().map(|(h,t)|EnthalpyPhaseKnot {
            specific_enthalpy_j_kg:h, temperature_k:t, liquid_mass_fraction:0.0, bulk_density_kg_m3:1000.0,
        }).collect()).unwrap()
}
fn layer<'a>(card:&'a MaterialCard, curve:&'a EquilibriumEnthalpyPhaseCurve,
    context:&'a QueryPoint) -> ResolvedConductionLayer<'a> {
    ResolvedConductionLayer::new(curve,ThermalMaterialPoint { material:card,context,
        temperature_axis:"T",liquid_fraction_axis:None },
        &PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),SolidLiquidPhase::Solid,
        [200.0,1000.0],||Ok(())).unwrap()
}
fn initial(curve:&EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events:Vec<_>=[(1,600.0),(2,300.0)].into_iter().map(|(id,t)|TimedEvent {
        time_s:0.0,event:Event::Birth(DepositedCell {
            id:CellId(id),segment_id:id,body_id:1,material_id:MaterialId(1),mass_kg:0.01,
            specific_h_j_kg:100.0*(t-200.0),
            geometry:CellGeometry { start_m:[id as f64,0.0,0.0],end_m:[id as f64,0.1,0.0],
                volume_m3:1e-5,surface_area_m2:0.1 },
            boundary:Boundary { ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0 },
        }),
    }).collect();
    events.push(TimedEvent { time_s:0.0,event:Event::Contact(ThermalLink {
        id:7,a:CellId(1),b:CellId(2),area_m2:0.004,conductance_w_k:0.0,
    }) });
    DepositionSimulation::new(events,&[(MaterialId(1),curve)],DepositionBudget {
        max_cells:2,max_links:1,max_events:3,max_substeps:100,
    }).unwrap().advance_to(0.0,StepConfig {max_dt_s:1.0,energy_tolerance_j:1e-8},||Ok(())).unwrap().0
}
fn binding<'a,'m>(layer:&'a ResolvedConductionLayer<'m>) -> LayerLinkBinding<'a,'m> {
    LayerLinkBinding { link:7,layers:[layer;2],distances_m:[0.02,0.03],contact_area_resistance_m2_k_w:0.0 }
}

#[test]
fn prepared_source_retains_every_breakpoint_and_shares_receipts_across_flux_trials() {
    let card=source();let curve=curve(&card);let context=QueryPoint::new().with("T",300.0).unwrap();
    let layer=layer(&card,&curve,&context);
    assert_eq!(layer.conductivity_curve().knots(),&[[200.0,4.0],[450.0,6.5],[1000.0,12.0]]);
    assert_eq!(layer.source_line().knots.len(),3);
    assert!(layer.source_line().knots.iter().all(|(_,a)|a.receipt.selected==layer.source_line().knots[0].1.receipt.selected));
    let state=initial(&curve);
    let prepared=PreparedLayerConduction::new(&state,vec![binding(&layer)],1,||Ok(())).unwrap();
    let first=prepared.resolve(&state,config(),||Ok(())).unwrap();
    assert_eq!(first.sources_w.len(),0);
    assert!((first.payload[0].applied_power_a_to_b_w-156.0).abs()<1e-9);
    assert!(Arc::ptr_eq(&first.payload[0].source_lines[0],layer.source_line()));
    let second=prepared.resolve(&state,config(),||Ok(())).unwrap();
    assert_eq!(first.payload,second.payload);
    assert_eq!(state.links()[&7].conductance_w_k,0.0);
    assert_eq!(state.time_s(),0.0);
}

#[test]
fn real_backward_euler_step_uses_integrated_flux_and_closes_internal_energy_once() {
    let card=source();let curve=curve(&card);let context=QueryPoint::new().with("T",300.0).unwrap();
    let layer=layer(&card,&curve,&context);let old=initial(&curve);
    let prepared=PreparedLayerConduction::new(&old,vec![binding(&layer)],1,||Ok(())).unwrap();
    // Equal 1 J/K capacities conserve Ta+Tb=900 K. Integrated linear k gives
    // G=A/L*(2+0.005*(Ta+Tb))=0.52 W/K throughout this trajectory.
    let dt=0.25;let difference=300.0/(1.0+2.0*0.52*dt);
    let exchange=0.5*(300.0-difference);
    let step=BackwardEulerStep::new(&old,dt,||Ok(())).unwrap();
    let trial=step.evaluate(&[-exchange,exchange],|s|prepared.resolve(s,config(),||Ok(())),||Ok(())).unwrap();
    assert!(trial.maximum_cell_residual_j()<1e-8);
    let (next,used,audit)=trial.accept(1e-8,1e-8,||Ok(())).unwrap();
    assert!((next.temperature_k(CellId(1)).unwrap()-(450.0+0.5*difference)).abs()<1e-10);
    assert!((next.temperature_k(CellId(2)).unwrap()-(450.0-0.5*difference)).abs()<1e-10);
    assert!((used[0].applied_power_a_to_b_w*dt-exchange).abs()<1e-8);
    assert_eq!(audit.external_heat_j,0.0);
    assert!(audit.residual_j.abs()<1e-8);
    assert_eq!(next.links()[&7],old.links()[&7]);
    assert_eq!(old.time_s(),0.0);
    assert_eq!(next.time_s(),dt);
    // Reusing a genuinely new endpoint still re-solves k(T), not a cached flux.
    assert!(prepared.resolve(&next,config(),||Ok(())).unwrap().payload[0].applied_power_a_to_b_w<156.0);
    let forged=step.evaluate(&[0.0,0.0],|s|prepared.resolve(s,config(),||Ok(())),||Ok(())).unwrap();
    assert!(forged.accept(1e-8,1e-8,||Ok(())).is_err());
}

#[test]
fn source_spans_cannot_cross_a_melting_plateau_or_use_foreign_caloric_objects() {
    let card=source();let context=QueryPoint::new().with("T",300.0).unwrap();
    let mixed=EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),
        [(0.0,200.0,0.0),(40_000.0,600.0,0.0),(60_000.0,600.0,1.0),(100_000.0,1000.0,1.0)]
        .into_iter().map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:f,bulk_density_kg_m3:1000.0}).collect()).unwrap();
    let make=|range|ResolvedConductionLayer::new(&mixed,ThermalMaterialPoint {material:&card,
        context:&context,temperature_axis:"T",liquid_fraction_axis:None},
        &PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),SolidLiquidPhase::Solid,range,||Ok(()));
    assert!(make([200.0,1000.0]).is_err());
    let solid=make([200.0,600.0]).unwrap();
    assert_eq!(solid.temperature_at_enthalpy(40_000.0).unwrap(),600.0);
    assert!(solid.temperature_at_enthalpy(40_000.0_f64.next_up()).is_err());
    let curve=curve(&card);let cloned=curve.clone();let layer=layer(&card,&cloned,&context);
    let state=initial(&curve);
    assert!(PreparedLayerConduction::new(&state,vec![binding(&layer)],1,||Ok(())).is_err());
}

#[test]
fn changed_contact_and_late_cancellation_refuse_without_overwriting_heat_history() {
    let card=source();let curve=curve(&card);let context=QueryPoint::new().with("T",300.0).unwrap();
    let layer=layer(&card,&curve,&context);let state=initial(&curve);
    let prepared=PreparedLayerConduction::new(&state,vec![binding(&layer)],1,||Ok(())).unwrap();
    let old=state.contact_patch(7).unwrap();
    let mut new=old;new.covered_area_m2=[0.005,0.006];
    let changed=state.replace_contact_patches(&[ContactPatchUpdate {id:7,expected:Some(old),replacement:Some(new)}],||Ok(())).unwrap();
    assert!(prepared.resolve(&changed,config(),||Ok(())).is_err());
    let mut polls=0;
    let expected=prepared.resolve(&state,config(),||{polls+=1;Ok(())}).unwrap();
    let mut calls=0;
    assert!(matches!(prepared.resolve(&state,config(),||{calls+=1;
        if calls==polls {Err(DepositionError::Cancelled)} else {Ok(())}
    }),Err(LayerSourceError::Thermal(DepositionError::Cancelled))));
    assert_eq!(prepared.resolve(&state,config(),||Ok(())).unwrap().payload,expected.payload);
    assert_eq!(state.contact_patch(7),Some(old));
    assert_eq!(state.time_s(),0.0);
    assert!(PreparedLayerConduction::new(&state,vec![binding(&layer)],0,||Ok(())).is_err());
}
