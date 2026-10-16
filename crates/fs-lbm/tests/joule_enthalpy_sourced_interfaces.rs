#![cfg(feature = "enthalpy-coupling")]
//! Immutable interface source cards through the shared transient current owner.
//! Synthetic numerical fixtures, not measured contact or alloy properties.
use std::cell::Cell;
use fs_blake3::ContentHash;
use fs_conduction::transient::deposition::*;
use fs_conduction::observers::ObserverCouplingError;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ElectricalContact, ElectricSetup, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError, JouleEnthalpyInput};
use fs_lbm::conjugate::enthalpy::transport::{EnthalpyTransport, EnthalpyTransportSample};
use fs_lbm::conjugate::enthalpy::transport::interfaces::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec,
    SystemContext, UncertaintyModel};
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_DIMS,
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, resolve_electrical_interface_state};
use fs_material::phase::SolidLiquidPhase;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Temperature;
use fs_lbm::conjugate::ChtError;

fn config(dt: f64) -> StepConfig { StepConfig {max_dt_s:dt, energy_tolerance_j:1e-8} }
const IDS: [CellId; 3] = [CellId(11), CellId(23), CellId(41)];
fn initial(law: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events: Vec<_> = IDS.iter().enumerate().map(|(c,&id)| TimedEvent {
        time_s:0.0, event:Event::Birth(DepositedCell {
            id, segment_id:c as u64, body_id:1, material_id:MaterialId(1), mass_kg:1.0,
            specific_h_j_kg:if c==1 {40.0} else {0.0},
            geometry:CellGeometry {start_m:[c as f64,0.5,0.5],end_m:[c as f64+1.0,0.5,0.5],
                volume_m3:1.0,surface_area_m2:6.0},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0},
        }),
    }).collect();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:7,a:IDS[0],b:IDS[1],area_m2:1.0,conductance_w_k:0.0,
    })});
    DepositionSimulation::new(events,&[(MaterialId(1),law as &dyn ThermalLaw)],
        DepositionBudget {max_cells:3,max_links:1,max_events:4,max_substeps:10000}).unwrap()
        .advance_to(0.0,config(0.125),||Ok(())).unwrap().0
}
fn setup(drive: ElectrodeDrive) -> ElectricSetup {
    ElectricSetup {conductors:vec![],electrodes:vec![
        Electrode {lo:[0,0,0],hi:[1,1,1],drive},
        Electrode {lo:[2,0,0],hi:[3,1,1],drive:ElectrodeDrive::Voltage(0.0)},
    ]}
}
fn face(resistance: f64) -> ElectricalContact {
    ElectricalContact {cells:[0,1],resistance_ohm_m2:resistance,heat_fraction_to_first:0.25}
}
fn input<'a>(domain:&'a VoxelDomain, setup:&'a ElectricSetup, contacts:&'a [ElectricalContact]) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {domain,setup,contacts,cell_ids:&IDS,max_cells:3,max_current_solves:10001,
        power_tolerance_w:1e-8,energy_tolerance_j:1e-8}
}
fn near(a:f64,b:f64) {assert!((a-b).abs()<1e-7,"{a} != {b}");}

fn surface(c: usize) -> SurfaceSpec {
    SurfaceSpec {material:MaterialStateId {chemistry:if c==1 {"specimen"} else {"electrode"}.into(),
        phase:"solid".into(),process:"synthetic fixture".into(),revision:0},texture_frame:"declared finish".into()}
}
fn sensible() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([0xE3;32]),SolidLiquidPhase::Solid,
        [(0.0,300.0),(200.0,400.0)].map(|(h,t)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,
            temperature_k:t,liquid_mass_fraction:0.0,bulk_density_kg_m3:1.0}).to_vec()).unwrap()
}
fn source(scale: f64, upper: f64) -> InterfaceSystemCard {
    let dims=ELECTRICAL_CONTACT_RESISTANCE_DIMS;let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,dims),
        value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(300.0,3.0*scale),(upper,scale*(3.0+0.03*(upper-300.0)))],dims},
        validity:ValidityDomain::unconstrained().with("T",300.0,upper).with("pressure",1.0,2.0),
        uncertainty:UncertaintyModel::RelativeHalfWidth {fraction:0.01,confidence:0.95},
        interpolation:InterpolationPolicy::LinearInside,observations:vec![],
        provenance:Provenance {source:"synthetic interface, not measured properties".into(),license:"test-only".into(),artifact:None},
    }).unwrap();
    InterfaceSystemCard::assemble(surface(0),surface(1),SystemContext {medium:"dry".into(),third_body:None,
        environment:"synthetic".into(),history:"virgin".into()},claims,vec![]).unwrap()
}
fn response(card:&InterfaceSystemCard, template:&ElectricalContact, sides:[&EnthalpyTransportSample<SurfaceSpec>;2])
    ->Result<ElectricalInterfaceResponse<fs_lbm::conjugate::material_coupling::interfaces::MaterialContact>,JouleEnthalpyError> {
    // This fixture EXPLICITLY uses the specimen-side temperature, never an
    // implicit mean. Full source pressure support is supplied independently.
    let point=QueryPoint::new().with("T",sides[1].temperature_k).unwrap().with("pressure",1.5).unwrap();
    let source=resolve_electrical_interface_state(card,&point,MaterialPropertySelection::SingleClaimOnly)
        .map_err(|e|ChtError::InvalidInput {field:"fixture.interface",reason:e.to_string()})?;
    ElectricalInterfaceResponse::try_from_source(template,source,sides.map(|s|s.response.material()),template.heat_fraction_to_first)
}

#[test]
fn source_data_changes_current_and_heat_with_matching_final_receipts() {
    let law=sensible();let old=initial(&law);let domain=VoxelDomain::new(3,1,1,1.0).unwrap();
    let electric=setup(ElectrodeDrive::Voltage(6.0));let faces=[face(0.0)];let gate=CancelGate::new();
    let mut results=Vec::new();
    for scale in [1.0,2.0] {
        let card=source(scale,400.0);
        let (next,(),r)=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.5,config(0.015625),input(&domain,&electric,&faces),&[],&gate,&(),
            |c,_,_|EnthalpyTransport::try_new(surface(c),Some(2.0),1.0),|_,_|Ok((0.0,())),
            |_,template,sides|response(&card,template,sides)).unwrap();
        let endpoint=&r.endpoint;let resolved=endpoint.contacts[0].response.material().source();
        let temperature=next.temperature_k(IDS[1]).unwrap();
        assert_eq!(resolved.resolved().state_coordinate("T"),Some(temperature));
        assert_eq!(resolved.resolved().state_coordinate("pressure"),Some(1.5));
        assert_eq!(resolved.resolved().card_identity(),card.content_hash());
        assert_eq!(endpoint.contacts[0].applied().resistance_ohm_m2,resolved.resistance_ohm_m2());
        near(resolved.resistance_ohm_m2(),scale*(3.0+0.03*(temperature-300.0)));
        near(r.transport.joule.endpoint.electric.electrodes[0].current_a,6.0/(2.0+resolved.resistance_ohm_m2()));
        let last=r.last_applied.as_ref().unwrap();
        assert!(last.contacts[0].response.material().source().resistance_ohm_m2()<resolved.resistance_ohm_m2());
        // Fresh physical solve from the retained endpoint, not a stale receipt.
        let fresh=electric.solve_with_cell_contacts(&domain,&r.transport.joule.endpoint.temperature_k,&[endpoint.contacts[0].applied()],
            |_,_,_|Ok(Some(2.0)),&gate).unwrap();
        assert_eq!(fresh,r.transport.joule.endpoint.electric);
        near(next.energy_j()-old.energy_j(),r.transport.joule.deposited_joule_j);
        results.push((temperature,r.transport.joule.endpoint.electric.electrodes[0].current_a,r.transport.joule.deposited_joule_j));
    }
    assert!(results[0].0>results[1].0 && results[0].1>results[1].1 && results[0].2>results[1].2);
}

#[test]
fn source_domain_failure_after_heating_and_surface_mismatch_publish_no_candidate() {
    let law=sensible();let old=initial(&law);let domain=VoxelDomain::new(3,1,1,1.0).unwrap();
    let electric=setup(ElectrodeDrive::Voltage(6.0));let faces=[face(0.0)];let gate=CancelGate::new();
    let calls=Cell::new(0);let card=source(1.0,320.1); // starts at 320 K, fails only after heat
    let result=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config(0.125),input(&domain,&electric,&faces),&[],&gate,&(),
        |c,_,_|EnthalpyTransport::try_new(surface(c),Some(2.0),1.0),|_,_|Ok((0.0,())),|_,template,sides| {
            calls.set(calls.get()+1);response(&card,template,sides)
        });
    assert!(matches!(result,Err(ObserverCouplingError::First(JouleEnthalpyError::Electric(ChtError::InvalidInput {field:"fixture.interface",..})))));
    assert_eq!(calls.get(),2);
    let card=source(1.0,400.0);
    for swap in [false,true] {
        let result=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config(0.125),input(&domain,&electric,&faces),&[],&gate,&(),
            |c,_,_| {
                let mut actual=surface(if swap {1-c.min(1)} else {c});
                if !swap && c==1 {actual.material.process="different processing history".into();}
                EnthalpyTransport::try_new(actual,Some(2.0),1.0)
            },|_,_|Ok((0.0,())),|_,template,sides|response(&card,template,sides));
        assert!(result.is_err());
    }
    near(old.energy_j(),40.0);assert_eq!(old.time_s(),0.0);
}

#[test]
fn source_contact_feedback_refines_to_analytic_enthalpy_solution() {
    let law=sensible();let old=initial(&law);let domain=VoxelDomain::new(3,1,1,1.0).unwrap();
    let electric=setup(ElectrodeDrive::Current(1.0));let faces=[face(0.0)];let gate=CancelGate::new();
    let card=source(1.0,400.0);
    // T=300+h/2, R=3+.03*(T-300), I=1: h'=4.25+.01125*h.
    let exact=(40.0+4.25/0.01125)*f64::exp(0.01125)-4.25/0.01125;
    let mut errors=Vec::new();
    for dt in [0.0625,0.03125,0.015625] {
        let (next,(),r)=advance_joule_enthalpy_with_material_interfaces_observed(&old,1.0,config(dt),
            input(&domain,&electric,&faces),&[],&gate,&(),
            |c,_,_|EnthalpyTransport::try_new(surface(c),Some(2.0),1.0),|_,_|Ok((0.0,())),
            |_,template,sides|response(&card,template,sides)).unwrap();
        errors.push((next.cells()[&IDS[1]].specific_h_j_kg-exact).abs());
        near(next.energy_j()-old.energy_j(),r.transport.joule.deposited_joule_j);
    }
    for error in errors.windows(2) {assert!(error[0]/error[1]>1.9 && error[0]/error[1]<2.1);}
}
