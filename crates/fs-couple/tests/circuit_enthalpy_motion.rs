//! G1/G3: actual field motion, material transport and finite enthalpy owners.
//! Synthetic constitutive fixtures test causality and physical accounting;
//! they do not establish a measured actuator's material or dynamic fidelity.
#[path = "circuit_enthalpy/common.rs"]
mod common;
#[path = "electrical_contacts/source.rs"]
mod source;
#[path = "circuit_enthalpy_motion/support.rs"]
mod support;

use common::*;
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve};
use fs_couple::circuit_enthalpy::material::*;
use fs_couple::circuit_enthalpy::contacts::material::*;
use fs_couple::circuit_enthalpy::motion::{self, FieldMotionEnthalpyInput};
use fs_couple::material_circuit::{MaterialCircuit, MaterialResistor};
use fs_couple::material_circuit::motion::{FieldMotionInput, MaterialFieldMotion};
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_ga::Vec3;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::field_storage::{ResolvedInductor, MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::{interface::ThermalInterfacePoint, source::ThermalMaterialPoint};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::Density;
use fs_thermal::enthalpy::{CellId, MaterialId};

fn cells<'a>(cards: &'a [MaterialCard; 2], curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint, key: &'a PropertyKey) -> [MaterialCell<'a>; 2] {
    [0, 1].map(|i| MaterialCell { cell: DESTINATIONS[i].cell, material_id: MaterialId(i as u64),
        curve: curves[i].curve(), point: ThermalMaterialPoint { material: &cards[i], context,
            temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)) })
}
fn resistors() -> [CellResistor; 2] {
    [0, 1].map(|i| CellResistor { destination: DESTINATIONS[i],
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: 1.0 })
}
fn links() -> [ConductionLink; 1] { [ConductionLink { link: 7, distance_a_m: 0.5, distance_b_m: 0.5 }] }
fn bases() -> [(usize, ElectricalPropertyBasis); 2] {
    [(1, ElectricalPropertyBasis::Resistivity), (2, ElectricalPropertyBasis::Resistivity)]
}
fn body(mass: f64) -> ResolvedBody {
    let card = source::scalar_card("density", Density::DIMS, mass);
    resolve_body(&[RegionBinding { key: "guided pole", material: &card, point: &QueryPoint::new(),
        selection: MaterialPropertySelection::SingleClaimOnly,
        geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap() }], 1, || false).unwrap()
}
fn gap_storage() -> MaterialStorage {
    let card = source::scalar_card(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, 1.0);
    MaterialStorage { branch: 0, element: MaterialStorageElement::Inductor(ResolvedInductor::resolve(
        &card, &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly, 1.0, 1.0, 1).unwrap()) }
}
fn finite_coil(resistors: Vec<MaterialResistor>, flux: f64) -> MaterialCircuit {
    MaterialCircuit::new_with_storage(CircuitGraph { node_count: 3, branches: vec![
        (1, 0, Branch::Inductor { henries: f64::NAN }),
        (1, 2, Branch::Resistor { ohms: f64::NAN }),
        (2, 0, Branch::Resistor { ohms: f64::NAN })], transformers: vec![] },
        resistors, vec![gap_storage()], vec![flux], vec![], &[]).unwrap()
}
fn input(model: &MaterialFieldMotion) -> FieldMotionEnthalpyInput<'_> {
    FieldMotionEnthalpyInput { motion: model, ports: &[], external_force_n: 0.0, destinations: &DESTINATIONS }
}
fn stored(model: &MaterialFieldMotion) -> f64 {
    model.inductor().unwrap().stored_energy_j(model.circuit().fluxes_wb()[0]).unwrap()
        +0.5*model.body().inertia().mass()*model.velocity_m_s().powi(2)
}

#[test]
fn heat_capacity_data_changes_resistivity_flux_and_motion_without_losing_energy() {
    let run = |cp| {
        let (a, ca) = material(1.0, cp, 0.005); let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key); let resistors = resistors(); let links = links();
        let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let mut thermal = initial(&curves, 999.0, None);
        let mut moving = MaterialFieldMotion::new(finite_coil(network.resolve(&thermal, || false).unwrap().resistors, 10.0),
            0, body(100.0), 0.0, [0.5, 1.5]).unwrap();
        let total = thermal.energy_j()+stored(&moving);
        let material = moving.inductor().unwrap().material().clone();
        let initial_resistors = moving.circuit().resistors().to_vec();
        let mut frozen = moving.clone();
        for n in 1..=8 {
            let target = f64::from(n)*0.1;
            let candidate = network.field_motion_candidate(&thermal, target, &moving, &[], 0.0,
                &bases(), config(), &[0.0; 2], || false).unwrap();
            let report = candidate.report();
            let temperature = candidate.thermal().temperature_k(CellId(30)).unwrap();
            let conductor = &report.electrical.resistors[0].material.conductor;
            near(conductor.resistance_ohm(), 1.0+0.005*(temperature-300.0));
            assert_eq!(conductor.source_property().answer().receipt.selected,
                cards[0].claims_for("electrical_resistivity")[0].0);
            near(report.source.links[0].conductance_w_k, 1.0);
            near(report.deposited_joule_j, report.electrical.heat_j);
            near(candidate.thermal().energy_j()+stored(candidate.motion()), total);
            near(report.total_energy_residual_j, 0.0);
            assert_eq!(candidate.motion().inductor().unwrap().material(), &material);
            frozen = frozen.step_candidate(initial_resistors.clone(), FieldMotionInput { electrical_ports: &[],
                external_force_n: 0.0, duration_s: target-thermal.time_s() }, config().circuit, || false).unwrap().into_parts().0;
            let (next_thermal, next_motion, _) = candidate.into_parts();
            thermal = next_thermal; moving = next_motion;
        }
        assert!((moving.gap_m()-frozen.gap_m()).abs()>1e-5,
            "real heating must change gap history compared with frozen resistance");
        (thermal.temperature_k(CellId(30)).unwrap(), moving.gap_m(), moving.circuit().fluxes_wb()[0])
    };
    let warm = run(1.0); let cooler = run(10.0);
    assert!(cooler.0<warm.0 && cooler.1<warm.1 && cooler.2>warm.2,
        "a larger caloric store suppresses resistance rise, retaining current and attraction");
}

#[test]
fn fixed_linkage_motion_and_endpoint_joule_feedback_match_independent_closed_forms() {
    let (a, ca) = material_with_conductivity(2.0, 2.0, 0.005, 0.0);
    let (b, cb) = material_with_conductivity(4.0, 3.0, 0.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key); let resistors = resistors(); let links = links();
    let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial(&curves, 999.0, None);
    let (flux, velocity, mass, force, dt, voltage) = (0.4, 0.02, 2.0, 0.03, 0.25, 3.0);
    // The winding's zero-voltage short conserves linkage. A second, independent
    // voltage port heats both fixed resistors. Its work must not drive the body.
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 3, branches: vec![
        (1, 0, Branch::Inductor { henries: f64::NAN }),
        (2, 0, Branch::Resistor { ohms: f64::NAN }),
        (2, 0, Branch::Resistor { ohms: f64::NAN }),
        (1, 0, Branch::VoltageSource { port: 0 }),
        (2, 0, Branch::VoltageSource { port: 1 })], transformers: vec![] },
        network.resolve(&thermal, || false).unwrap().resistors, vec![gap_storage()], vec![flux], vec![], &[0.0, voltage]).unwrap();
    let model = MaterialFieldMotion::new(circuit, 0, body(mass), velocity, [0.5, 1.5]).unwrap();
    let candidate = network.field_motion_candidate(&thermal, dt, &model, &[0.0, voltage], force,
        &bases(), config(), &[0.0; 2], || false).unwrap();
    let attraction = -0.5*flux*flux;
    let expected_gap = 1.0+velocity*dt+0.5*(force+attraction)/mass*dt*dt;
    let expected_velocity = velocity+(force+attraction)/mass*dt;
    let rhs = dt*voltage*voltage/(2.0*2.0);
    let temperature_rise = 2.0*rhs/(1.0+(1.0+4.0*0.005*rhs).sqrt());
    near(candidate.thermal().temperature_k(CellId(30)).unwrap(), 300.0+temperature_rise);
    near(candidate.thermal().temperature_k(CellId(10)).unwrap(), 300.0+dt*voltage*voltage/(3.0*4.0));
    near(candidate.motion().gap_m(), expected_gap);
    near(candidate.motion().velocity_m_s(), expected_velocity);
    near(candidate.circuit().fluxes_wb()[0], flux);
    let report = candidate.report();
    near(report.electrical.field_work_on_body_j, attraction*(expected_gap-1.0));
    near(report.electrical.mechanical_supply_j, force*(expected_gap-1.0));
    near(report.electrical.electrical_supply_j, report.deposited_joule_j);
    near(report.deposited_joule_j, 2.0*temperature_rise+dt*voltage*voltage/4.0);
    near(report.non_electrical_heat_j, 0.0);
    near(candidate.thermal().energy_j()-thermal.energy_j()+stored(candidate.motion())-stored(&model),
        report.electrical.electrical_supply_j+force*(expected_gap-1.0));
}

#[test]
fn contact_material_heat_partition_and_feedback_change_the_real_actuator() {
    let run = |resistance| {
        let (a, ca) = material(1.0, 1.0, 0.005); let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0); let interface_context = source::point();
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key); let resistors = resistors(); let links = links();
        let bulk = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let card = source::card(&cards[0], &cards[1], resistance, 0.005);
        let bindings = [MaterialElectricalContact { link: 7, branch: 1,
            point: ThermalInterfacePoint { interface: &card, materials: [&cards[0], &cards[1]],
                context: &interface_context, temperature_axes: ["Ta", "Tb"],
                liquid_fraction_axes: [Some("melt_a"), Some("melt_b")], mean_temperature_axis: None },
            selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: 0.25 }];
        let bases = bases();
        let network = ElectricalContactNetwork { bulk: &bulk, contacts: &bindings, electrical_bases: &bases };
        let mut thermal = initial(&curves, 999.0, None);
        let initial = network.resolve(&thermal, || false).unwrap();
        let circuit = finite_coil(initial.bulk.resistors, 10.0).rebind_contacts_candidate(initial.contacts, &[]).unwrap();
        let mut moving = MaterialFieldMotion::new(circuit, 0, body(100.0), 0.0, [0.5, 1.5]).unwrap();
        let total = thermal.energy_j()+stored(&moving);
        for n in 1..=4 {
            let target = f64::from(n)*0.1; let dt = target-thermal.time_s();
            let next = network.field_motion_candidate(&thermal, target, &moving, &[], 0.0, config(), &[0.0; 2], || false).unwrap();
            let report = next.report();
            let t_a = next.thermal().temperature_k(CellId(30)).unwrap();
            let t_b = next.thermal().temperature_k(CellId(10)).unwrap();
            let contact = &report.electrical.contacts[0];
            near(contact.contact.resistance_ohm(), resistance*(1.0+0.005*(t_a-300.0)));
            let exchanged = dt*report.source.material.links[0].conductance_w_k*(t_b-t_a);
            let heat_a = next.thermal().cells()[&CellId(30)].specific_h_j_kg-thermal.cells()[&CellId(30)].specific_h_j_kg;
            let heat_b = next.thermal().cells()[&CellId(10)].specific_h_j_kg-thermal.cells()[&CellId(10)].specific_h_j_kg;
            near(heat_a, report.electrical.resistors[0].heat_j+0.25*contact.heat_j+exchanged);
            near(heat_b, report.electrical.resistors[1].heat_j+0.75*contact.heat_j-exchanged);
            near(report.deposited_joule_j, report.electrical.heat_j);
            near(next.thermal().energy_j()+stored(next.motion()), total);
            assert_eq!(report.source.destinations[0].cells, [CellId(30), CellId(10)]);
            let (t, m, _) = next.into_parts(); thermal = t; moving = m;
        }
        (moving.gap_m(), moving.circuit().fluxes_wb()[0], thermal.temperature_k(CellId(30)).unwrap())
    };
    let low = run(0.2); let high = run(2.0);
    assert!(high.0>low.0 && high.1<low.1);
    assert!((high.2-low.2).abs()>0.01, "changing only interface data must alter real finite-body heat");
}

#[test]
fn material_budget_and_final_cancellation_refuse_all_three_histories() {
    let (a, ca) = material(1.0, 1.0, 0.005); let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key); let resistors = resistors(); let links = links();
    let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial(&curves, 999.0, None);
    let model = MaterialFieldMotion::new(finite_coil(network.resolve(&thermal, || false).unwrap().resistors, 10.0),
        0, body(100.0), 0.0, [0.5, 1.5]).unwrap();
    let saved_cells = thermal.cells().clone(); let saved_motion = model.clone();
    let mut polls = 0;
    let normal = network.field_motion_candidate(&thermal, 0.1, &model, &[], 0.0, &bases(), config(), &[0.0; 2],
        || { polls+=1; false }).unwrap();
    let mut actual = 0;
    assert!(matches!(network.field_motion_candidate(&thermal, 0.1, &model, &[], 0.0, &bases(), config(), &[0.0; 2],
        || { actual+=1; actual==polls }), Err(CircuitEnthalpyError::Cancelled)));
    let limited = CircuitEnthalpySolve { maximum_circuit_trials: normal.report().circuit_trials-1, ..config() };
    assert!(matches!(network.field_motion_candidate(&thermal, 0.1, &model, &[], 0.0, &bases(), limited, &[0.0; 2], || false),
        Err(CircuitEnthalpyError::CircuitBudget)));
    let prepared = network.prepare_with_bases(&thermal, &bases(), || false).unwrap();
    assert!(matches!(motion::step_candidate(&thermal, 0.1, input(&model), config(), &[0.0; 2], |state| {
        if state.temperature_k(CellId(30)).unwrap()>300.01 { return Err("winding temperature outside source support"); }
        prepared.resolve(state, || false).map_err(|_| "material source refused")
    }, || false), Err(CircuitEnthalpyError::Source("winding temperature outside source support"))));
    assert_eq!(thermal.cells(), &saved_cells); assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(model, saved_motion);
}
