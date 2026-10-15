//! Fully data-bound electrical/thermal material substitutions on shared history.
#[path = "circuit_enthalpy/common.rs"]
mod common;
use common::*;
use fs_couple::circuit_enthalpy::material::*;
use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_material::state_point::{MaterialPropertySelection,THERMAL_CONDUCTIVITY_DIMS};
use fs_matdb::{MaterialCard,PropertyKey,QueryPoint};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve,ResolvedHeatCapacityCurve};
use fs_thermal::enthalpy::*;

fn bound_cells<'a>(cards:&'a [MaterialCard;2],curves:&'a [ResolvedHeatCapacityCurve;2],
    context:&'a QueryPoint,key:&'a PropertyKey) -> [MaterialCell<'a>;2] {
    [0,1].map(|i|MaterialCell {cell:DESTINATIONS[i].cell,material_id:MaterialId(i as u64),curve:curves[i].curve(),
        point:ThermalMaterialPoint {material:&cards[i],context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((key,MaterialPropertySelection::SingleClaimOnly))})
}
fn branches() -> [CellResistor;2] {
    [0,1].map(|i|CellResistor {destination:DESTINATIONS[i],selection:MaterialPropertySelection::SingleClaimOnly,
        length_m:1.0,area_m2:1.0})
}
fn links() -> [ConductionLink;1] { [ConductionLink {link:7,distance_a_m:0.5,distance_b_m:0.5}] }

#[test]
fn conductivity_only_material_change_propagates_through_other_body_and_circuit() {
    let run=|k| {
        let (a,ca)=material_with_conductivity(1.0,1.0,0.005,k);
        let (b,cb)=material_with_conductivity(3.0,2.0,0.0,2.0);
        let cards=[a,b];let curves=[ca,cb];let context=point(300.0);
        let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cells=bound_cells(&cards,&curves,&context,&key);let resistors=branches();let links=links();
        let binding=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
        // Authored link coefficient is deliberately wrong: the source binding
        // must replace it at every trial rather than use it as a fallback.
        let mut thermal=initial(&curves,999.0,None);let mut electrical=circuit(&thermal,&cards,&curves);
        let total0=thermal.energy_j()+50.0;
        for n in 1..=40 {
            let candidate=binding.step_candidate(&thermal,n as f64*0.05,&electrical,&[],config(),&[0.0;2],||false).unwrap();
            let report=candidate.report();
            near(report.source.links[0].conductance_w_k,1.0/(0.5/k+0.5/2.0));
            for sample in &report.source.cells {
                near(sample.phase.temperature_k(),candidate.thermal().temperature_k(sample.cell).unwrap());
                assert!(sample.conductivity.is_some());
            }
            near(candidate.thermal().energy_j()+report.electrical.stored_after_j,total0);
            let (t,c,_)=candidate.into_parts();thermal=t;electrical=c;
        }
        (thermal.temperature_k(CellId(30)).unwrap(),thermal.temperature_k(CellId(10)).unwrap(),electrical.charges_c()[0])
    };
    let a=run(0.2);let b=run(2.0);
    assert!(b.0>a.0 && b.1<a.1);
    assert!((a.2-b.2).abs()>1e-6,"conductivity must alter temperature-dependent current history");
}

#[test]
fn wrong_caloric_instance_identity_geometry_or_missing_transport_refuses() {
    let (a,ca)=material(1.0,1.0,0.005);let (b,cb)=material(3.0,2.0,0.0);
    let cards=[a,b];let curves=[ca,cb];let context=point(300.0);
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let mut cells=bound_cells(&cards,&curves,&context,&key);let mut resistors=branches();let links=links();
    let thermal=initial(&curves,1.0,None);let electrical=circuit(&thermal,&cards,&curves);
    let copied=curves[0].curve().clone();cells[0].curve=&copied;
    assert!(MaterialNetwork {cells:&cells,resistors:&resistors,links:&links}.resolve(&thermal,||false).is_err());
    cells[0].curve=curves[0].curve();cells[0].point.material=&cards[1];
    assert!(MaterialNetwork {cells:&cells,resistors:&resistors,links:&links}.resolve(&thermal,||false).is_err());
    cells[0].point.material=&cards[0];cells[0].conductivity=None;
    assert!(MaterialNetwork {cells:&cells,resistors:&resistors,links:&links}.resolve(&thermal,||false).is_err());
    cells[0].conductivity=Some((&key,MaterialPropertySelection::SingleClaimOnly));resistors[0].area_m2=2.0;
    assert!(MaterialNetwork {cells:&cells,resistors:&resistors,links:&links}.resolve(&thermal,||false).is_err());
    resistors[0].area_m2=1.0;
    let binding=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let mut polls=0;
    binding.step_candidate(&thermal,0.1,&electrical,&[],config(),&[0.0;2],||{polls+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(binding.step_candidate(&thermal,0.1,&electrical,&[],config(),&[0.0;2],
        ||{actual+=1;actual==polls}),Err(fs_couple::circuit_enthalpy::CircuitEnthalpyError::Cancelled)));
    assert_eq!(thermal.time_s(),0.0);assert_eq!(electrical.charges_c(),&[10.0]);
}

fn phase_card(r:f64) -> MaterialCard {
    use fs_matdb::*;use fs_qty::Dims;use fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS;
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS),
        value:PropertyValue::Curve {abscissa:"melt".into(),abscissa_dims:Dims([0;6]),
            knots:vec![(0.0,r),(1.0,2.0*r)],dims:ELECTRICAL_RESISTIVITY_DIMS},
        validity:fs_evidence::ValidityDomain::unconstrained().with("T",200.0,800.0).with("melt",0.0,1.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,observations:vec![],
        provenance:Provenance {source:"authored phase-resistivity fixture".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"declared by chart".into(),
        process:"authored".into(),revision:0},claims,vec![]).unwrap()
}
#[test]
fn latent_energy_data_changes_resistance_and_other_body_heat_at_identical_kelvin() {
    use fs_material::phase::{EnthalpyPhaseKnot,SolidLiquidPhase};
    use fs_phs::circuit::{CircuitGraph,Branch};use fs_couple::material_circuit::MaterialCircuit;
    let run=|latent:f64| {
        let cards=[phase_card(1.0),phase_card(3.0)];
        let knot=|h,t,x|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:x,bulk_density_kg_m3:1.0};
        // Explicit supplied phase-chart data, not latent heat inferred from Cp.
        let curves=[EquilibriumEnthalpyPhaseCurve::try_new(cards[0].content_hash(),vec![
            knot(-100.0,250.0,0.0),knot(0.0,300.0,0.0),knot(latent,300.0,1.0),knot(latent+1000.0,800.0,1.0)]).unwrap(),
            EquilibriumEnthalpyPhaseCurve::try_single_phase(cards[1].content_hash(),SolidLiquidPhase::Solid,
                vec![knot(-100.0,250.0,0.0),knot(1000.0,800.0,0.0)]).unwrap()];
        let context=QueryPoint::new().with("T",300.0).unwrap().with("melt",0.0).unwrap();
        let cells=[0,1].map(|i|MaterialCell {cell:DESTINATIONS[i].cell,material_id:MaterialId(i as u64),curve:&curves[i],
            point:ThermalMaterialPoint {material:&cards[i],context:&context,temperature_axis:"T",liquid_fraction_axis:Some("melt")},conductivity:None});
        let resistors=branches();let binding=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
        let events=[0,1].map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id:DESTINATIONS[i].cell,segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:1.0,
            specific_h_j_kg:if i==0 {0.25*latent} else {0.0},geometry:CellGeometry {
                start_m:[i as f64,0.0,0.0],end_m:[i as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:6.0},boundary:boundary()})}).to_vec();
        let thermal=DepositionSimulation::new(events,&[(MaterialId(0),&curves[0] as &dyn ThermalLaw),
            (MaterialId(1),&curves[1] as &dyn ThermalLaw)],DepositionBudget {max_cells:2,max_links:0,max_events:2,max_substeps:10}).unwrap()
            .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-7},||Ok(())).unwrap().0;
        let circuit=MaterialCircuit::new(CircuitGraph {node_count:3,branches:vec![
            (1,0,Branch::Capacitor {farads:1.0}),(1,2,Branch::Resistor {ohms:0.0}),(2,0,Branch::Resistor {ohms:0.0})],transformers:vec![]},
            binding.resolve(&thermal,||false).unwrap().resistors,vec![],vec![10.0],&[]).unwrap();
        let candidate=binding.step_candidate(&thermal,0.25,&circuit,&[],config(),&[0.0;2],||false).unwrap();
        let first=candidate.report().source.cells.iter().find(|s|s.cell==CellId(30)).unwrap();
        assert_eq!(first.phase.temperature_k(),300.0);assert!(first.phase.liquid_mass_fraction()>0.25);
        near(candidate.report().electrical.resistors[0].material.conductor.resistance_ohm(),1.0+first.phase.liquid_mass_fraction());
        (first.phase.liquid_mass_fraction(),candidate.report().electrical.resistors[1].heat_j)
    };
    let a=run(20.0);let b=run(80.0);
    assert!(a.0>b.0 && a.1<b.1);
}
