//! G1/G3: original field U/S, independently calculated RC/RL and finite-layer heat.
//! Synthetic source properties, not a measured-material qualification.
use super::*;
use crate::circuit_enthalpy::field_storage::FieldThermalBinding;
use crate::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_material::field_storage::{ResolvedCapacitor, ResolvedInductor,
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY,
    MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
    thermal::ThermalFieldStorage};

#[derive(Clone, Copy)]
struct Spec { magnetic: bool, slope: f64 }
fn spec() -> Spec { Spec { magnetic: false, slope: 0.005 } }
fn coefficient(s: Spec, t: f64) -> f64 {
    (if s.magnetic { 0.2 } else { CAPACITANCE }) * (1.0 + s.slope * (t-300.0))
}
fn initial_coordinate(s: Spec) -> f64 { if s.magnetic { 0.002 } else { INITIAL_CHARGE } }
fn physical(s: Spec, t: f64, z: f64) -> (f64, f64, f64) {
    let psi = z*z/(2.0*coefficient(s,t));
    let entropy = psi*s.slope/(1.0+s.slope*(t-300.0));
    (psi, entropy, psi+t*entropy)
}
fn source(layer: usize, s: Spec) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (card, _) = powered_source(layer,1.0,false);
    let mut claims = ClaimSet::new();
    for (_, claim) in card.claims().claims_ordered() { claims.insert_claim(claim.clone()).unwrap(); }
    let (name,dims) = if s.magnetic { (MAGNETIC_PERMEABILITY_PROPERTY,MAGNETIC_PERMEABILITY_DIMS) }
        else { (ELECTRIC_PERMITTIVITY_PROPERTY,ELECTRIC_PERMITTIVITY_DIMS) };
    let value = coefficient(s,300.0)/area(layer); // Original area and 1 m gap/path.
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
        value:if s.slope==0.0 {PropertyValue::Scalar {value,dims}} else {
            PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(290.0,value*(1.0-10.0*s.slope)),(340.0,value*(1.0+40.0*s.slope))],dims}},
        validity:ValidityDomain::unconstrained().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),290.0,340.0),
        interpolation:if s.slope==0.0 {InterpolationPolicy::ConstantWithinValidity} else {InterpolationPolicy::LinearInside},
        uncertainty:UncertaintyModel::Unstated,observations:vec![],provenance:Provenance {
            source:"manufactured field/relaxation specimen".into(),license:"test-only".into(),artifact:None},
    }).unwrap();
    let card=MaterialCard::assemble(card.id().clone(),claims,vec![]).unwrap();
    let curve=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",
        lower:&point(290.0),upper:&point(340.0),selection:SelectionPolicy::SingleClaimOnly},
        SolidLiquidPhase::Solid,0.0,1e-6).unwrap();
    (card,curve)
}
fn fixture<R>(s:Spec,run:impl FnOnce(&MaterialFluidFrame<'_,FluidTransportState>,
    &MaterialNetwork<'_>,&MaterialNetwork<'_>,DepositionSimulation<'_>,
    &[(usize,ElectricalPropertyBasis)],MaterialCircuit,Vec<FieldThermalBinding>)->R)->R {
    let (a,ac)=source(0,s);let (b,bc)=source(1,s);let cards=[a,b];let curves=[ac,bc];
    let f=frame(&cards);let m=model(&f);let query=point(300.0);let heat=thermal(&curves,1.0);
    let cells:[_;2]=core::array::from_fn(|i|MaterialCell {cell:CellId(i as u64),material_id:MaterialId(i as u64),
        curve:curves[i].curve(),point:ThermalMaterialPoint {material:&cards[i],context:&query,
            temperature_axis:"T",liquid_fraction_axis:None},conductivity:None});
    let resistors:[_;2]=core::array::from_fn(|i|CellResistor {destination:ResistorCell {branch:i+1,cell:CellId(i as u64)},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:area(i)});
    let plain=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    let powered=MaterialNetwork {cells:plain.cells,resistors:&resistors,links:plain.links};
    let bases=[(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)];
    let (branch,element,law)=if s.magnetic {
        let value=ResolvedInductor::resolve(&cards[0],&query,MaterialPropertySelection::SingleClaimOnly,area(0),1.0,1).unwrap();
        let law=ThermalFieldStorage::inductor(&cards[0],&value,"T",[290.0,340.0]).unwrap();
        (Branch::Inductor {henries:f64::NAN},MaterialStorageElement::Inductor(value),law)
    } else {
        let value=ResolvedCapacitor::resolve(&cards[0],&query,MaterialPropertySelection::SingleClaimOnly,area(0),1.0).unwrap();
        let law=ThermalFieldStorage::capacitor(&cards[0],&value,"T",[290.0,340.0]).unwrap();
        (Branch::Capacitor {farads:f64::NAN},MaterialStorageElement::Capacitor(value),law)
    };
    let resolved=powered.prepare_with_bases(&heat,&bases,||false).unwrap().resolve(&heat,||false).unwrap();
    let z=initial_coordinate(s);
    let circuit=MaterialCircuit::new_with_storage(CircuitGraph {node_count:3,branches:vec![
        (1,0,branch),(1,2,Branch::Resistor {ohms:f64::NAN}),(2,0,Branch::Resistor {ohms:f64::NAN})],transformers:vec![]},
        resolved.resistors,vec![MaterialStorage {branch:0,element}],
        if s.magnetic {vec![z]} else {vec![]},if s.magnetic {vec![]} else {vec![z]},&[]).unwrap();
    let fields=vec![FieldThermalBinding {branch:0,cell:CellId(0),law,maximum_abs_coordinate:1.01*z}];
    run(&m,&plain,&powered,heat,&bases,circuit,fields)
}
fn coordinate(circuit:&MaterialCircuit,s:Spec)->f64 {
    if s.magnetic {circuit.fluxes_wb()[0]} else {circuit.charges_c()[0]}
}
fn solve()->FrameThermalSolve {
    let mut c=controls();c.maximum_material_trials=256;c.maximum_newton_iterations=20;
    c.difference_step_j=1e-10;c.cell_energy_tolerance_j=1e-12;c.thermal_energy_tolerance_j=1e-12;
    c.newton.absolute_tolerance=1e-13;c
}
fn copy_fields(fields:&[FieldThermalBinding])->Vec<FieldThermalBinding> {
    fields.iter().map(|f|FieldThermalBinding {branch:f.branch,cell:f.cell,law:f.law.clone(),
        maximum_abs_coordinate:f.maximum_abs_coordinate}).collect()
}

#[test]
fn actual_source_fields_drive_rc_rl_caloric_exchange_and_continuing_material_motion() {
    for magnetic in [false,true] {
        let mut responses=Vec::new();
        for slope in [0.0,0.005,-0.005] {
            let s=Spec {magnetic,slope};
            fixture(s,|m,plain,powered,mut heat,bases,circuit,fields| {
                let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
                let initial=native.total_energy_j().unwrap()+caloric(&heat)+physical(s,300.0,initial_coordinate(s)).2;
                let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
                owner.bind_material_field_circuit(&heat,powered,circuit, bases,fields,&gate).unwrap();
                let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
                force.resize(owner.mechanical().mode_count(),0.0);
                let (mut work,mut outside,mut reversible)=(0.0,0.0,0.0);
                for _ in 0..12 {
                    let cells=heat.cells().clone();let z0=coordinate(owner.electrical().unwrap(),s);
                    let t0=heat.temperature_k(CellId(0)).unwrap();
                    let report=owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
                        solve(),&[0.0;2],&gate).unwrap();
                    let t=heat.temperature_k(CellId(0)).unwrap();let z=coordinate(owner.electrical().unwrap(),s);
                    let e=report.electrical.as_ref().unwrap();let r:f64=e.resistors.iter().map(|r|r.material.conductor.resistance_ohm()).sum();
                    let ratio=if magnetic {r/(512.0*coefficient(s,t))} else {1.0/(512.0*r*coefficient(s,t))};
                    assert!((z-z0*(1.0-ratio)/(1.0+ratio)).abs()<1e-12);
                    let before=physical(s,t0,z0);let after=physical(s,t,z);let frozen=physical(s,t,z0);
                    let qrev=before.0-frozen.0-(t*after.1-t0*before.1);
                    let field=report.field_storage.as_ref().unwrap();
                    assert!((field.reversible_heat_j-qrev).abs()<1e-12);
                    assert!((field.field_internal_energy_j[0]-before.2).abs()<1e-12);
                    assert!((field.field_internal_energy_j[1]-after.2).abs()<1e-12);
                    for i in 0..2 {
                        let id=CellId(i as u64);let c=&heat.cells()[&id];
                        let loss=report.layers.iter().find(|l|l.destination.cell==id).unwrap().heat_j;
                        let joule=e.resistors.iter().find(|r|r.material.branch==i+1).unwrap().heat_j;
                        assert!((c.mass_kg*(c.specific_h_j_kg-cells[&id].specific_h_j_kg)-loss-joule
                            -if i==0 {qrev} else {0.0}).abs()<3e-12);
                    }
                    let entropy=report.field_entropy.unwrap();
                    assert!((entropy.storage_change_j_k-(after.1-before.1)).abs()<1e-13);
                    assert!((entropy.reversible_transfer_j_k-qrev/t).abs()<1e-13);
                    assert!(!report.entropy.has_unclassified_sources);
                    assert!(report.entropy.accounting_residual_j_k.abs()<1e-12);
                    assert!(report.external_thermal_heat_j.abs()<1e-12);
                    work+=report.mechanical.external_work_j+e.native.supplied;
                    outside+=report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
                    assert!((owner.mechanical().total_energy_j().unwrap()+caloric(&heat)+after.2-initial+outside-work).abs()<2e-10);
                    reversible+=qrev;
                }
                if slope==0.0 {assert_eq!(reversible,0.0);} else {assert!(reversible*slope>0.0);}
                responses.push(placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0]);
            });
        }
        assert!((responses[0]-responses[1]).abs()>1e-12);
        assert!((responses[0]-responses[2]).abs()>1e-12);
    }
}

#[test]
fn complete_original_field_capacity_and_receiver_are_required_before_binding() {
    fixture(spec(),|m,plain,powered,heat,bases,circuit,fields| {
        let gate=CancelGate::new_clock_free();let (_,material,native)=native(m);
        let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
        assert!(owner.bind_material_circuit(&heat,powered,circuit.clone(),bases,&gate).is_err());
        for fault in 0..4 {
            let mut bad=copy_fields(&fields);
            match fault {0=>bad.clear(),1=>bad[0].cell=CellId(1),2=>bad[0].maximum_abs_coordinate=0.0,
                _=>bad[0].maximum_abs_coordinate=1.0}
            assert!(owner.bind_material_field_circuit(&heat,powered,circuit.clone(),bases,bad,&gate).is_err());
            assert!(owner.electrical().is_none());assert_eq!(owner.mechanical().samples_rendered(),0);
        }
        owner.bind_material_field_circuit(&heat,powered,circuit.clone(),bases,fields,&gate).unwrap();
        assert_eq!(owner.electrical().unwrap(),&circuit);
    });
}

#[test]
fn failed_field_thermal_trial_retains_charge_flux_sources_and_every_material_history() {
    for magnetic in [false,true] {
        fixture(Spec {magnetic,..spec()},|m,plain,powered,mut heat,bases,circuit,fields| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
            owner.bind_material_field_circuit(&heat,powered,circuit.clone(),bases,fields,&gate).unwrap();
            let cells=heat.cells().clone();let memory=owner.mechanical().material_memory().to_vec();
            let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
            let source=owner.material().members()[0].source.clone();let force=vec![0.0;owner.mechanical().mode_count()];
            let mut bad=solve();bad.maximum_material_trials=1;
            assert!(owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},bad,&[0.0;2],&gate).is_err());
            let cancelled=CancelGate::new_clock_free();cancelled.request();
            assert!(owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},solve(),&[0.0;2],&cancelled).is_err());
            assert_eq!(owner.electrical().unwrap(),&circuit);assert_eq!(heat.cells(),&cells);
            assert_eq!(owner.mechanical().material_memory(),memory);assert_eq!(owner.material().members()[0].source,source);
            assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
            assert_eq!(owner.mechanical().samples_rendered(),0);assert_eq!(heat.time_s(),0.0);
            owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},solve(),&[0.0;2],&gate).unwrap();
            assert_ne!(coordinate(owner.electrical().unwrap(),Spec {magnetic,..spec()}),initial_coordinate(Spec {magnetic,..spec()}));
        });
    }
}
