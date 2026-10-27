//! Same-source coupled radiation, thermal contacts, electrical losses and motion.
//! Synthetic cards; no measured-material or renderer-fidelity claim.
#[path = "../../fs-conduction/tests/radiation_material/common.rs"] mod radiation_fixture;
#[path = "../../fs-conduction/tests/contact_material_state/common.rs"] mod contact_fixture;
#[path = "material_exchange/gas.rs"] mod gas_fixture;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::NetworkSolve;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::{MaterialContact, exchange::*};
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell};
use fs_couple::circuit_enthalpy::material::*;
use fs_couple::circuit_enthalpy::contacts::material::*;
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_conduction::transient::deposition::contact::{apply_material_contacts, ContactCell, MaterialContactLink};
use fs_conduction::transient::deposition::radiation::{MaterialRadiation, RadiatingCell};
use fs_conduction::SURFACE_EMISSIVITY_PROPERTY;
use fs_evidence::ValidityDomain;
use fs_exec::{Cx, CancelGate};
use fs_material::conductor::{ElectricalPropertyBasis, ELECTRICAL_RESISTIVITY_DIMS};
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, ELECTRICAL_CONTACT_RESISTANCE_DIMS};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_matdb::*;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;

fn with_cx<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, fs_exec::StreamKey { seed: 73, kernel_id: 59, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        run(&cx, &gate)
    })
}
fn card(emissivity: [f64;2], r: f64) -> MaterialCard {
    let base = radiation_fixture::card(emissivity[0], emissivity[1]);
    let mut claims = base.claims().clone();
    for (name, dims, a, b) in [
        ("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, r, 2.0*r),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.2*r, 0.4*r),
    ] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0,a),(500.0,b)], dims },
            validity: ValidityDomain::unconstrained().with("T",300.0,500.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "synthetic coupled exchange".into(), license: "CC0-1.0".into(), artifact: None } }).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn fixture<R>(emissivity: [f64;2], first_h: f64,
    run: impl FnOnce(&MaterialExchangeNetwork<'_>, DepositionSimulation<'_>, MaterialCircuit) -> R) -> R
{
    let cards = [card(emissivity,1.0),card([0.6,0.6],2.0)];
    let curves = cards.each_ref().map(|c| radiation_fixture::curve(c,100.0));
    let context = radiation_fixture::context();
    let ids = [CellId(30),CellId(10)];
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0,1].map(|i| MaterialCell { cell: ids[i], material_id: MaterialId(i as u64), curve: &curves[i],
        point: radiation_fixture::binding(&cards[i],&context), conductivity: Some((&key,MaterialPropertySelection::SingleClaimOnly)) });
    let resistors = [0,1].map(|i| CellResistor { destination: ResistorCell { branch:i+1,cell:ids[i] },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m:1.0, area_m2:1.0 });
    let links = [ConductionLink { link:7, distance_a_m:0.5, distance_b_m:0.5 }];
    let bulk = MaterialNetwork { cells:&cells, resistors:&resistors, links:&links };
    let (thermal_card,_) = contact_fixture::source(&cards,0.2);
    let mut claims = thermal_card.claims().clone();
    claims.insert_claim(PropertyClaim { key:PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,ELECTRICAL_CONTACT_RESISTANCE_DIMS),
        value:PropertyValue::Curve { abscissa:"T".into(), abscissa_dims:Temperature::DIMS,
            knots:vec![(300.0,0.04),(500.0,0.06)], dims:ELECTRICAL_CONTACT_RESISTANCE_DIMS },
        validity:ValidityDomain::unconstrained().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),300.0,500.0)
            .with("pressure",200000.0,200000.0),
        uncertainty:UncertaintyModel::Unstated, interpolation:InterpolationPolicy::LinearInside, observations:vec![],
        provenance:Provenance { source:"synthetic common thermal/electrical interface".into(),license:"CC0-1.0".into(),artifact:None } }).unwrap();
    let interface = contact_fixture::interface(&cards,claims);
    let interface_context = contact_fixture::point();
    let interface_point = contact_fixture::binding(&cards,&interface,&interface_context);
    let electrical_contacts = [MaterialElectricalContact { link:7, branch:1, point:interface_point,
        selection:MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first:0.25 }];
    let thermal_contacts = [MaterialContactLink { link:7,
        sides:[0,1].map(|i| ContactCell { cell:ids[i], material_id:MaterialId(i as u64),curve:&curves[i] }),
        point:interface_point, selection:ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) }];
    let bases = [(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)];
    let electrical = ElectricalContactNetwork { bulk:&bulk,contacts:&electrical_contacts,electrical_bases:&bases };
    let radiating = [0,1].map(|i| RadiatingCell { surface:i,cell:ids[i],material_id:MaterialId(i as u64),curve:&curves[i],
        material:cells[i].point,selection:MaterialPropertySelection::SingleClaimOnly });
    let enclosure = radiation_fixture::enclosure([
        radiation_fixture::emissivity(&cards[0],&curves[0],first_h),
        radiation_fixture::emissivity(&cards[1],&curves[1],-50.0),
    ]);
    let radiation = MaterialRadiation { enclosure:&enclosure,cells:&radiating,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8 };
    let assembly = MaterialExchangeNetwork { electrical:&electrical,thermal_contacts:&thermal_contacts,radiation:Some(&radiation) };
    let mut events = [0,1].map(|i| TimedEvent { time_s:0.0,event:Event::Birth(DepositedCell {
        id:ids[i],segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:1.0,
        specific_h_j_kg:if i==0 {first_h} else {-50.0},geometry:CellGeometry {
            start_m:[i as f64,0.0,0.0],end_m:[i as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:1.5 },
        boundary:Boundary { ambient_k:400.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0 },
    }) }).to_vec();
    events.push(TimedEvent { time_s:0.0,event:Event::Contact(ThermalLink { id:7,a:ids[0],b:ids[1],area_m2:0.1,conductance_w_k:999.0 }) });
    let thermal = DepositionSimulation::new(events,&[(MaterialId(0),&curves[0] as &dyn ThermalLaw),(MaterialId(1),&curves[1] as &dyn ThermalLaw)],
        DepositionBudget { max_cells:2,max_links:1,max_events:3,max_substeps:1000 }).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-7},||Ok(())).unwrap().0;
    let initial = electrical.resolve(&thermal,||false).unwrap();
    let circuit = MaterialCircuit::new_with_contacts(CircuitGraph { node_count:3,branches:vec![
        (1,0,Branch::Capacitor {farads:0.01}),(1,2,Branch::Resistor {ohms:999.0}),(2,0,Branch::Resistor {ohms:999.0})],transformers:vec![] },
        initial.bulk.resistors,initial.contacts,vec![],vec![1.0],&[]).unwrap();
    run(&assembly,thermal,circuit)
}
fn controls() -> NetworkSolve {
    NetworkSolve { network:CircuitEnthalpySolve { maximum_cells:2,maximum_circuit_trials:300,maximum_newton_iterations:32,
        difference_step_j:1e-4,cell_energy_tolerance_j:1e-7,thermal_energy_tolerance_j:1e-7,total_energy_tolerance_j:1e-7,
        circuit:MaterialCircuitSolve {energy_tolerance_j:1e-7,power_tolerance_w:1e-7},
        newton:fs_solver::NewtonKrylovConfig {absolute_tolerance:1e-9,relative_tolerance:1e-14,forcing_minimum:1e-4,
            linear_restart:4,..fs_solver::NewtonKrylovConfig::default()} },
        maximum_chambers:2,maximum_contacts:2,maximum_boundary_solves:600,total_energy_tolerance_j:1e-6 }
}
fn near(a:f64,b:f64) { assert!((a-b).abs()<2e-6,"{a} != {b}"); }

#[test]
fn data_only_emissivity_reaches_both_bodies_circuit_and_both_gas_motions() {
    let run = |e| fixture([e,e],150.0,|assembly,mut thermal,mut circuit| with_cx(|cx,_| {
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let mut gases = gas_fixture::pistons();
        let before = gases.clone(); let energy0 = thermal.energy_j()+50.0; let mut work = 0.0;
        for _ in 0..40 {
            let target = thermal.time_s()+0.00005; let dt = target-thermal.time_s();
            let chambers = [gas_fixture::chamber(10,&gases[0],dt),gas_fixture::chamber(20,&gases[1],dt)];
            let result = prepared.step_candidate(cx,&thermal,target,&circuit,&[],&chambers,&gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
            let source = &result.material().material;
            assert_eq!(source.contacts.len(),1);
            assert!(source.contacts[0].conductance_w_k<source.contacts[0].bulk_conductance_w_k);
            let radiation = source.radiation.as_ref().unwrap();
            assert!(radiation.net_heat_w().abs()<1e-8);
            for surface in radiation.surfaces() {
                assert_eq!(surface.phase,source.bulk.cells.iter().find(|c|c.cell==surface.cell).unwrap().phase);
            }
            for chamber in result.chambers() { for surface in chamber.boundary.source() {
                assert_eq!(source.bulk.cells.iter().find(|c|c.cell==surface.source.0.cell).unwrap().conductivity.as_ref(),Some(&surface.source.1));
            }}
            near(result.step().external_network_heat_j,0.0);
            assert!(result.step().total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
            let (t,c,report,_) = result.into_parts();
            assert!(report.electrical.contacts[0].heat_j>0.0);
            for (i,chamber) in report.source.chambers.into_iter().enumerate() {
                work+=chamber.boundary.step().piston.exterior_pressure_work_j;
                gases[i]=chamber.boundary.into_parts().0;
            }
            let gas_energy:f64=gases.iter().zip(&before).map(|(g,o)|g.gas().mass_kg()
                *(g.gas().state().specific_internal_energy_j_kg()-o.gas().state().specific_internal_energy_j_kg())
                +0.5*g.mass_kg()*g.velocity_m_s().powi(2)).sum();
            assert!((t.energy_j()+report.electrical.stored_after_j+gas_energy+work-energy0).abs()<4e-5);
            thermal=t; circuit=c;
        }
        (thermal.temperature_k(CellId(30)).unwrap(),thermal.temperature_k(CellId(10)).unwrap(),circuit.charges_c()[0],gases)
    }));
    let low=run(0.2); let high=run(0.8);
    assert!((low.0-high.0).abs()>1e-3 && (low.1-high.1).abs()>1e-3);
    assert!((low.2-high.2).abs()>1e-8);
    for i in 0..2 { assert!((low.3[i].position_m()-high.3[i].position_m()).abs()>1e-13); }
}

#[test]
fn omitted_extra_mechanisms_reduce_to_original_contact_thermal_solver() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit| with_cx(|cx,_| {
        let simple=MaterialExchangeNetwork {electrical:assembly.electrical,thermal_contacts:&[],radiation:None};
        let prepared=simple.prepare(cx,&thermal).unwrap();
        let combined=prepared.step_candidate(cx,&thermal,0.00005,&circuit,&[],&[],&[],controls(),&[0.0;2]).unwrap();
        let original=assembly.electrical.step_candidate(&thermal,0.00005,&circuit,&[],controls().network,&[0.0;2],||false).unwrap();
        assert_eq!(combined.thermal().cells(),original.thermal().cells());
        assert_eq!(combined.circuit(),original.circuit());
        assert!(combined.material().material.contacts.is_empty());
        assert!(combined.material().material.radiation.is_none());
    }));
}

#[test]
fn combined_coefficients_match_original_owners_and_do_not_accumulate_between_trials() {
    fixture([0.2,0.8],25.0,|assembly,thermal,circuit| with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let first=prepared.resolve(cx,&thermal).unwrap();
        let second=prepared.resolve(cx,&thermal).unwrap();
        assert_eq!(first.bulk.conductances_w_k,second.bulk.conductances_w_k);
        assert_eq!(first.bulk.sources_w,second.bulk.sources_w);
        let mut original=assembly.electrical.resolve(&thermal,||false).unwrap();
        let contacts=apply_material_contacts(&thermal,assembly.thermal_contacts,&mut original.bulk.conductances_w_k,||false).unwrap();
        let radiation=assembly.radiation.unwrap().resolve(cx,&thermal).unwrap();
        radiation.add_to(&mut original.bulk.sources_w).unwrap();
        assert_eq!(first.bulk.conductances_w_k,original.bulk.conductances_w_k);
        assert_eq!(first.bulk.sources_w,original.bulk.sources_w);
        near(first.bulk.source.contacts[0].conductance_w_k,contacts[0].conductance_w_k);
        let candidate=prepared.step_candidate(cx,&thermal,0.00005,&circuit,&[],&[],&[],controls(),&[0.0;2]).unwrap();
        let sample=&candidate.material().material;
        let surface=sample.radiation.as_ref().unwrap().surfaces().iter().find(|s|s.cell==CellId(30)).unwrap();
        assert_eq!(surface.phase.temperature_k(),400.0);
        let fraction=surface.phase.liquid_mass_fraction();
        assert!((fraction-0.25).abs()>1e-7);
        near(surface.material.property(SURFACE_EMISSIVITY_PROPERTY).unwrap().value_si(),0.2+0.6*fraction);
        near(sample.contacts[0].interface.value_m2_k_per_w(),0.2*(1.0+fraction));
    }));
}

#[test]
fn mismatched_surface_conditions_duplicates_area_and_cancellation_publish_nothing() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit| with_cx(|cx,gate| {
        let before_cells=thermal.cells().clone(); let before_circuit=circuit.clone();
        let duplicates=[assembly.thermal_contacts[0];2];
        assert!(MaterialExchangeNetwork {electrical:assembly.electrical,thermal_contacts:&duplicates,radiation:assembly.radiation}
            .prepare(cx,&thermal).is_err());
        let original=assembly.radiation.unwrap();
        let changed=original.cells[0].material.context.clone().with("lot",8.0).unwrap();
        let mut surfaces=original.cells.iter().map(|s|RadiatingCell {surface:s.surface,cell:s.cell,material_id:s.material_id,
            curve:s.curve,material:s.material,selection:s.selection.clone()}).collect::<Vec<_>>();
        surfaces[0].material.context=&changed;
        let inconsistent=MaterialRadiation {enclosure:original.enclosure,cells:&surfaces,
            power_tolerance_w:original.power_tolerance_w,radiosity_tolerance_w_m2:original.radiosity_tolerance_w_m2};
        assert!(MaterialExchangeNetwork {electrical:assembly.electrical,thermal_contacts:assembly.thermal_contacts,radiation:Some(&inconsistent)}
            .prepare(cx,&thermal).is_err());
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let gases=gas_fixture::pistons(); let chambers=[gas_fixture::chamber(10,&gases[0],0.00005),gas_fixture::chamber(20,&gases[1],0.00005)];
        let mut too_large=gas_fixture::patches();too_large[0].geometry.area_m2=2.0;
        assert!(prepared.step_candidate(cx,&thermal,0.00005,&circuit,&[],&chambers,&too_large,controls(),&[0.0;2]).is_err());
        gate.request();
        assert!(matches!(prepared.resolve(cx,&thermal),Err(MaterialExchangeError::Cancelled)));
        assert!(matches!(prepared.step_candidate(cx,&thermal,0.00005,&circuit,&[],&[],&[],controls(),&[0.0;2]),Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.cells(),&before_cells);assert_eq!(thermal.time_s(),0.0);assert_eq!(circuit,before_circuit);
    }));
}

#[path = "material_exchange/reactive.rs"]
mod reactive_tests;

#[path = "material_exchange/coupled_material.rs"]
mod coupled_material;

#[path = "material_exchange/coupled_core.rs"]
mod coupled_core;
