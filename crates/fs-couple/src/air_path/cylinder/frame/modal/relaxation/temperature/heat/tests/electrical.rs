//! Existing native electrical and Maxwell histories share their ORIGINAL cells.
use super::*;
use super::super::electrical::RelaxationElectricalInput;
use crate::circuit_enthalpy::{ResistorCell, material::CellResistor};
use crate::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_material::conductor::{ElectricalPropertyBasis, ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_DIMS};
use fs_phs::circuit::{Branch, CircuitGraph};

const CAPACITANCE: f64 = 0.001;
const INITIAL_CHARGE: f64 = 0.0002;
fn area(layer: usize) -> f64 {
    core::f64::consts::PI * (RADII[layer].powi(2) - if layer == 0 { 0.0 } else { RADII[layer-1].powi(2) })
}
fn powered_source(layer: usize, scale: f64, conductivity: bool) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = source(layer, true);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() { claims.insert_claim(claim.clone()).unwrap(); }
    let (name, dims) = if conductivity {
        ("electrical_conductivity", ELECTRICAL_CONDUCTIVITY_DIMS)
    } else { ("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS) };
    let rho = [0.2, 0.6][layer] * scale;
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: [290.0, 300.0, 340.0].into_iter().map(|t| {
                let resistivity = rho * (1.0 + 0.02 * (t-300.0));
                (t, if conductivity { 1.0/resistivity } else { resistivity })
            }).collect(), dims },
        validity: ValidityDomain::unconstrained().with("T",290.0,340.0),
        interpolation: InterpolationPolicy::LinearInside, uncertainty: UncertaintyModel::Unstated,
        observations: vec![], provenance: Provenance { source: "manufactured electrical/relaxation fixture".into(),
            license: "test-only".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(290.0), upper: &point(340.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, curve)
}
fn fixture_powered<R>(scale: f64, conductivity: bool,
    run: impl FnOnce(&MaterialFluidFrame<'_,FluidTransportState>, &MaterialNetwork<'_>,
        &MaterialNetwork<'_>, DepositionSimulation<'_>, &[(usize,ElectricalPropertyBasis)]) -> R) -> R {
    let (a,ac) = powered_source(0,scale,conductivity);
    let (b,bc) = powered_source(1,1.0,false);
    let cards = [a,b]; let curves = [ac,bc]; let f = frame(&cards); let m = model(&f);
    let query = point(300.0);
    let cells: [_;2] = core::array::from_fn(|i| MaterialCell {
        cell: CellId(i as u64), material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context: &query, temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: None,
    });
    let resistors: [_;2] = core::array::from_fn(|i| CellResistor {
        destination: ResistorCell { branch: i+1, cell: CellId(i as u64) },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: area(i),
    });
    let plain = MaterialNetwork { cells: &cells, links: &[], resistors: &[] };
    let powered = MaterialNetwork { cells: plain.cells, links: plain.links, resistors: &resistors };
    let bases = [(1, if conductivity { ElectricalPropertyBasis::Conductivity } else { ElectricalPropertyBasis::Resistivity }),
        (2, ElectricalPropertyBasis::Resistivity)];
    run(&m, &plain, &powered, thermal(&curves,1.0), &bases)
}
fn circuit(source: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>,
    bases: &[(usize,ElectricalPropertyBasis)], current: Option<f64>) -> MaterialCircuit {
    let resistors = source.prepare_with_bases(thermal,bases,||false).unwrap()
        .resolve(thermal,||false).unwrap().resistors;
    let mut branches = vec![(1,0,Branch::Capacitor { farads: CAPACITANCE }),
        (1,2,Branch::Resistor { ohms: 999.0 }), (2,0,Branch::Resistor { ohms: 999.0 })];
    if current.is_some() { branches.push((0,1,Branch::CurrentSource { port: 0 })); }
    MaterialCircuit::new(CircuitGraph { node_count: 3, branches, transformers: vec![] },
        resistors, vec![], vec![INITIAL_CHARGE], &current.into_iter().collect::<Vec<_>>()).unwrap()
}
fn electrical_controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve { energy_tolerance_j: 1e-11, power_tolerance_w: 1e-10 }
}

#[test]
fn finite_charge_drives_same_layer_temperature_rates_and_motion_without_duplicate_heat() {
    let mut outcomes = Vec::new();
    for scale in [0.5,2.0] {
        fixture_powered(scale,false,|m,plain,powered,mut heat,bases| {
            let gate = CancelGate::new_clock_free(); let (placement,material,native) = native(m);
            let initial = native.total_energy_j().unwrap() + caloric(&heat)
                + 0.5*INITIAL_CHARGE.powi(2)/CAPACITANCE;
            let mut owner = FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
            owner.bind_material_circuit(&heat,powered,circuit(powered,&heat,bases,None),bases,&gate).unwrap();
            let mut force = placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
            force.resize(owner.mechanical().mode_count(),0.0);
            let (mut outside,mut work,mut joule) = (0.0,0.0,0.0);
            for _ in 0..12 {
                let before = heat.cells().clone(); let q0 = owner.electrical().unwrap().charges_c()[0];
                let report = owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {
                    ports: &[], solve: electrical_controls() },controls(),&[0.0;2],&gate).unwrap();
                let electrical = report.electrical.as_ref().unwrap();
                let resistance: f64 = electrical.resistors.iter().map(|r| r.material.conductor.resistance_ohm()).sum();
                let dt = owner.mechanical().sample_period_s();
                let a = dt/(2.0*resistance*CAPACITANCE);
                let expected = q0*(1.0-a)/(1.0+a);
                let q1 = owner.electrical().unwrap().charges_c()[0];
                assert!((q1-expected).abs()<1e-12);
                let current = (q0+q1)/(2.0*CAPACITANCE*resistance);
                for loss in &electrical.resistors {
                    let cell = CellId((loss.material.branch-1) as u64);
                    let relaxation = report.layers.iter().find(|r|r.destination.cell==cell).unwrap().heat_j;
                    let delta = heat.cells()[&cell].mass_kg*(heat.cells()[&cell].specific_h_j_kg-before[&cell].specific_h_j_kg);
                    assert!((loss.heat_j-dt*current*current*loss.material.conductor.resistance_ohm()).abs()<1e-12);
                    assert!((delta-relaxation-loss.heat_j).abs()<1e-10);
                    near(loss.material.conductor.material().query().axes()["T"],heat.temperature_k(cell).unwrap());
                }
                near(report.layers.iter().map(|l|l.heat_j).sum(),report.deposited_material_heat_j);
                near(report.deposited_joule_heat_j,electrical.heat_j);
                joule += report.deposited_joule_heat_j;
                work += report.mechanical.external_work_j + electrical.native.supplied;
                outside += report.external_mechanical_loss_j + report.mechanical.numerical_dissipation_j;
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)+electrical.stored_after_j-initial+outside-work,0.0);
                near(report.total_energy_residual_j,0.0);
                assert!(!report.rates.changes.is_empty());
            }
            assert!(joule>0.0 && owner.electrical().unwrap().charges_c()[0]<INITIAL_CHARGE);
            near(owner.mechanical().material_memory()[0],0.007*(-12.0_f64/256.0/0.7).exp());
            outcomes.push((placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0],
                heat.temperature_k(CellId(0)).unwrap(),owner.electrical().unwrap().charges_c()[0]));
        });
    }
    assert!((outcomes[0].0-outcomes[1].0).abs()>1e-12);
    assert!((outcomes[0].1-outcomes[1].1).abs()>1e-5);
    assert!((outcomes[0].2-outcomes[1].2).abs()>1e-8);
}

#[test]
fn conductivity_only_source_and_held_current_close_native_supply_and_combined_storage() {
    fixture_powered(1.0,true,|m,plain,powered,mut heat,bases| {
        let gate = CancelGate::new_clock_free(); let (_,material,native) = native(m);
        let initial = native.total_energy_j().unwrap()+caloric(&heat)+0.5*INITIAL_CHARGE.powi(2)/CAPACITANCE;
        let mut owner = FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
        let electrical = circuit(powered,&heat,bases,Some(0.001));
        owner.bind_material_circuit(&heat,powered,electrical,bases,&gate).unwrap();
        let force = vec![0.0;owner.mechanical().mode_count()];
        let (mut supplied,mut outside) = (0.0,0.0);
        for _ in 0..8 {
            let q0 = owner.electrical().unwrap().charges_c()[0];
            let r = owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {
                ports: &[0.001], solve: electrical_controls() },controls(),&[0.0;2],&gate).unwrap();
            let e = r.electrical.as_ref().unwrap(); let q1 = owner.electrical().unwrap().charges_c()[0];
            assert!((e.native.supplied-0.001*(q0+q1)/(2.0*CAPACITANCE)/256.0).abs()<1e-12);
            let conductor = &e.resistors.iter().find(|r|r.material.branch==1).unwrap().material.conductor;
            assert_eq!(conductor.source_basis(),ElectricalPropertyBasis::Conductivity);
            assert!(conductor.material().property("electrical_resistivity").is_none());
            supplied += r.mechanical.external_work_j+e.native.supplied+r.external_thermal_heat_j;
            outside += r.external_mechanical_loss_j+r.mechanical.numerical_dissipation_j;
            near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)+e.stored_after_j-initial+outside-supplied,0.0);
        }
        assert!(supplied>0.0);
    });
}

#[test]
fn late_electrothermal_refusal_preserves_charge_every_history_and_all_sources_for_retry() {
    fixture_powered(1.0,false,|m,plain,powered,mut heat,bases| {
        let gate = CancelGate::new_clock_free(); let (placement,material,native) = native(m);
        let mut owner = FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
        owner.bind_material_circuit(&heat,powered,circuit(powered,&heat,bases,None),bases,&gate).unwrap();
        let circuit = owner.electrical().unwrap().clone(); let cells = heat.cells().clone();
        let physical = placement.physical_state(owner.mechanical().components(),&gate).unwrap();
        let memory = owner.mechanical().material_memory().to_vec(); let source = owner.material().members()[0].source.clone();
        let force = vec![0.0;owner.mechanical().mode_count()];
        assert!(owner.step(&mut heat,&force,controls(),&[0.0;2],&gate).is_err());
        let mut exhausted = controls(); exhausted.maximum_material_trials = 1;
        assert!(owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
            exhausted,&[0.0;2],&gate).is_err());
        let cancelled = CancelGate::new_clock_free(); cancelled.request();
        assert!(owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
            controls(),&[0.0;2],&cancelled).is_err());
        assert_eq!(owner.electrical().unwrap(),&circuit); assert_eq!(heat.cells(),&cells);
        assert_eq!(heat.time_s(),0.0); assert_eq!(owner.mechanical().samples_rendered(),0);
        assert_eq!(owner.mechanical().material_memory(),memory);
        assert_eq!(owner.material().members()[0].source,source);
        assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
        owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
            controls(),&[0.0;2],&gate).unwrap();
        assert_eq!(owner.mechanical().samples_rendered(),1);
        assert_ne!(owner.electrical().unwrap().charges_c(),circuit.charges_c());
    });
}

#[test]
fn circuit_binding_refuses_separate_cells_bad_basis_and_late_rebinding_without_publication() {
    fixture_powered(1.0,false,|m,plain,powered,mut heat,bases| {
        let empty = MaterialNetwork {cells:&[],resistors:powered.resistors,links:powered.links};
        let gate = CancelGate::new_clock_free(); let (_,material,native) = native(m);
        let mut owner = FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
        let circuit = circuit(powered,&heat,bases,None);
        assert!(owner.bind_material_circuit(&heat,powered,circuit.clone(),&bases[..1],&gate).is_err());
        assert!(owner.electrical().is_none());
        assert!(owner.bind_material_circuit(&heat,&empty,circuit.clone(),bases,&gate).is_err());
        assert!(owner.electrical().is_none());
        owner.bind_material_circuit(&heat,powered,circuit.clone(),bases,&gate).unwrap();
        assert!(owner.bind_material_circuit(&heat,powered,circuit.clone(),bases,&gate).is_err());
        assert_eq!(owner.electrical().unwrap(),&circuit);
        let force = vec![0.0;owner.mechanical().mode_count()];
        owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {ports:&[],solve:electrical_controls()},
            controls(),&[0.0;2],&gate).unwrap();
        let accepted = owner.electrical().unwrap().clone();
        assert!(owner.bind_material_circuit(&heat,powered,circuit,bases,&gate).is_err());
        assert_eq!(owner.electrical().unwrap(),&accepted);
    });
}

#[test]
fn finite_flux_and_charge_follow_native_rlc_midpoint_with_same_endpoint_materials() {
    fixture_powered(1.0,false,|m,plain,powered,mut heat,bases| {
        let gate = CancelGate::new_clock_free(); let (_,material,native) = native(m);
        let inductance = 0.2; let flux = 0.00003;
        let initial = native.total_energy_j().unwrap()+caloric(&heat)
            +0.5*INITIAL_CHARGE.powi(2)/CAPACITANCE+0.5*flux*flux/inductance;
        let resistors = powered.prepare_with_bases(&heat,bases,||false).unwrap().resolve(&heat,||false).unwrap().resistors;
        let circuit = MaterialCircuit::new(CircuitGraph {node_count:4,branches:vec![
            (1,0,Branch::Capacitor {farads:CAPACITANCE}),
            (1,2,Branch::Resistor {ohms:999.0}), (2,3,Branch::Resistor {ohms:999.0}),
            (3,0,Branch::Inductor {henries:inductance})],transformers:vec![]},
            resistors,vec![flux],vec![INITIAL_CHARGE],&[]).unwrap();
        let mut owner = FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
        owner.bind_material_circuit(&heat,powered,circuit,bases,&gate).unwrap();
        let force = vec![0.0;owner.mechanical().mode_count()]; let mut outside = 0.0;
        for _ in 0..12 {
            let electrical = owner.electrical().unwrap();
            let (q0,f0) = (electrical.charges_c()[0],electrical.fluxes_wb()[0]);
            let report = owner.step_electrical(&mut heat,&force,RelaxationElectricalInput {
                ports:&[],solve:electrical_controls()},controls(),&[0.0;2],&gate).unwrap();
            let step = report.electrical.as_ref().unwrap();
            let resistance:f64 = step.resistors.iter().map(|r|r.material.conductor.resistance_ohm()).sum();
            let dt = owner.mechanical().sample_period_s();
            let (a,b,c) = (dt/(2.0*inductance),dt/(2.0*CAPACITANCE),dt*resistance/(2.0*inductance));
            let expected_flux = (2.0*b*q0+(1.0-c-a*b)*f0)/(1.0+c+a*b);
            let expected_charge = q0-a*(f0+expected_flux);
            let accepted = owner.electrical().unwrap();
            assert!((accepted.fluxes_wb()[0]-expected_flux).abs()<1e-12);
            assert!((accepted.charges_c()[0]-expected_charge).abs()<1e-12);
            let current = (f0+accepted.fluxes_wb()[0])/(2.0*inductance);
            for loss in &step.resistors {
                assert!((loss.heat_j-dt*current*current*loss.material.conductor.resistance_ohm()).abs()<1e-12);
            }
            outside += report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
            near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)+step.stored_after_j-initial+outside,0.0);
        }
        assert!(heat.temperature_k(CellId(0)).unwrap()>300.0);
    });
}

mod contacts;
mod transport;
