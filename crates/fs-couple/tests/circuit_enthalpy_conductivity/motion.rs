//! G1/G3/G4: original field dynamics, full source bundles and finite receivers.
use super::*;
#[path = "../electrical_contacts/source.rs"]
mod scalar;
use fs_couple::circuit_enthalpy::motion::{FieldMotionEnthalpyInput, conductivity};
use fs_couple::material_circuit::motion::MaterialFieldMotion;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_ga::Vec3;
use fs_material::field_storage::{ResolvedCapacitor, ELECTRIC_PERMITTIVITY_DIMS,
    ELECTRIC_PERMITTIVITY_PROPERTY};
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};

fn moving(resistors: Vec<MaterialResistor>, driven: bool) -> MaterialFieldMotion {
    let dielectric = scalar::scalar_card(ELECTRIC_PERMITTIVITY_PROPERTY,ELECTRIC_PERMITTIVITY_DIMS,1.0);
    let capacitor = ResolvedCapacitor::resolve(&dielectric,&QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,1.0,1.0).unwrap();
    let material = scalar::scalar_card("density",Density::DIMS,100.0);
    let body = resolve_body(&[RegionBinding {key:"guided electrode",material:&material,
        point:&QueryPoint::new(),selection:MaterialPropertySelection::SingleClaimOnly,
        geometry:VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap()}],1,||false).unwrap();
    let branches = if driven {
        vec![(1,0,Branch::Capacitor {farads:f64::NAN}),
            (2,0,Branch::Resistor {ohms:f64::NAN}),
            (2,0,Branch::Resistor {ohms:f64::NAN}),
            (2,0,Branch::VoltageSource {port:0})]
    } else {
        vec![(1,0,Branch::Capacitor {farads:f64::NAN}),
            (1,2,Branch::Resistor {ohms:f64::NAN}),
            (2,0,Branch::Resistor {ohms:f64::NAN})]
    };
    let graph = CircuitGraph {node_count:3,branches,transformers:vec![]};
    let charge = if driven {0.8} else {10.0};
    let ports: &[f64] = if driven {&[3.0]} else {&[]};
    let circuit = MaterialCircuit::new_with_storage(graph,resistors,
        vec![MaterialStorage {branch:0,element:MaterialStorageElement::Capacitor(capacitor)}],
        vec![],vec![charge],ports).unwrap();
    MaterialFieldMotion::new(circuit,0,body,if driven {0.02} else {0.0},[0.5,1.5]).unwrap()
}
fn stored(motion: &MaterialFieldMotion) -> f64 {
    let q = motion.circuit().charges_c()[0];
    0.5*q*q*motion.gap_m()+0.5*100.0*motion.velocity_m_s().powi(2)
}

#[test]
fn independent_constant_charge_motion_and_nonlinear_heating_share_external_work() {
    let beta = 0.05;
    let (a,ca) = field_material(beta,2000.0);
    let (b,cb) = common::material(3.0,2.0,0.005);
    let cards = [a,b]; let curves = [ca,cb]; let thermal = initial(&curves,0.0,None);
    let model = moving(circuit(&cards).resistors().to_vec(),true);
    let dt = 0.2; let force = 0.1;
    let input = FieldMotionEnthalpyInput {motion:&model,ports:&[3.0],
        external_force_n:force,destinations:&DESTINATIONS};
    let result = conductivity::step_candidate(&thermal,dt,input,&responses(&cards),config(),
        conductivity_controls(),&[0.0;2],|s|coefficients(s,model.circuit(),&cards,&curves),||false).unwrap();
    let heat_a = dt*9.0/(1.0+beta*3.0);
    let rhs = dt*9.0/(2.0*3.0);
    let rise_b = 2.0*rhs/(1.0+(1.0+4.0*0.005*rhs).sqrt());
    let acceleration = (force-0.5*0.8*0.8)/100.0;
    let displacement = 0.02*dt+0.5*acceleration*dt*dt;
    near(result.thermal().temperature_k(CellId(30)).unwrap(),300.0+heat_a);
    near(result.thermal().temperature_k(CellId(10)).unwrap(),300.0+rise_b);
    near(result.motion().gap_m(),1.0+displacement);
    near(result.motion().velocity_m_s(),0.02+acceleration*dt);
    near(result.motion().circuit().charges_c()[0],0.8);
    near(result.report().electrical.evaluated[0].point.electric_field_v_m,3.0);
    near(result.report().deposited_joule_j,heat_a+2.0*rise_b);
    let field = &result.report().electrical.field;
    near(field.electrical_supply_j,heat_a+2.0*rise_b);
    near(field.mechanical_supply_j,force*displacement);
    near(result.thermal().energy_j()-thermal.energy_j()+stored(result.motion())-stored(&model),
        field.electrical_supply_j+force*displacement);
}

#[test]
fn prepared_material_swaps_change_heat_transport_charge_and_actual_gap_without_callbacks() {
    let mut endpoints = Vec::new();
    for beta in [0.02,0.08] {
        let (a,ca) = thermo_field_material(beta); let (b,cb) = common::material(3.0,2.0,0.005);
        let cards = [a,b]; let curves = [ca,cb]; let mut thermal = initial(&curves,1.0,None);
        let context = [query(300.0),point(300.0)];
        let key = PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cells = network_cells(&cards,&curves,&context,&key); let resistors = network_resistors();
        let network = MaterialNetwork {cells:&cells,resistors:&resistors,links:&NETWORK_LINKS};
        let responses = network_responses();
        let prepared = network.prepare_conductivity(&thermal,&responses,||false).unwrap();
        let seeds = prepared.resolve_conductivity_seeds(&thermal,&responses,||false).unwrap();
        let mut motion = moving(seeds.resistors,false); let total = thermal.energy_j()+stored(&motion);
        for interval in 1..=3 {
            let target = f64::from(interval)*0.1; let dt = target-thermal.time_s();
            let result = prepared.field_motion_conductivity_candidate(&thermal,target,&motion,&[],0.0,
                &responses,config(),conductivity_controls(),&[0.0;2],||false).unwrap();
            let ta = result.thermal().temperature_k(CellId(30)).unwrap();
            let tb = result.thermal().temperature_k(CellId(10)).unwrap();
            let electrical = &result.report().electrical;
            let operating = electrical.evaluated[0].point.electric_field_v_m.abs();
            let resistance = 1.0+0.004*(ta-300.0)+beta*operating;
            near(electrical.evaluated[0].conductor.resistance_ohm(),resistance);
            let k = 1.0+0.002*(ta-300.0)+0.3*operating;
            let conductance = 1.0/(0.5/k+0.5);
            near(result.report().source.links[0].conductance_w_k,conductance);
            let heat_a = result.thermal().cells()[&CellId(30)].specific_h_j_kg
                -thermal.cells()[&CellId(30)].specific_h_j_kg;
            near(heat_a,electrical.field.resistors[0].heat_j+dt*conductance*(tb-ta));
            near(result.report().deposited_joule_j,electrical.field.heat_j);
            near(result.thermal().energy_j()+stored(result.motion()),total);
            assert!(result.report().circuit_trials>electrical.iterations);
            let (t,m,_) = result.into_parts(); thermal=t; motion=m;
        }
        endpoints.push((motion.gap_m(),motion.circuit().charges_c()[0],
            thermal.temperature_k(CellId(30)).unwrap()));
    }
    assert!(endpoints[1].1>endpoints[0].1,"more field-dependent resistance retains charge");
    assert!(endpoints[1].0<endpoints[0].0,"retained charge strengthens actual attraction");
    assert!((endpoints[1].2-endpoints[0].2).abs()>0.01);
}

#[test]
fn all_native_trials_count_and_refusal_preserves_thermal_and_moving_histories() {
    let (a,ca) = field_material(0.05,2000.0); let (b,cb) = common::material(3.0,2.0,0.005);
    let cards = [a,b]; let curves = [ca,cb]; let thermal = initial(&curves,1.0,None);
    let model = moving(circuit(&cards).resistors().to_vec(),false);
    let input = FieldMotionEnthalpyInput {motion:&model,ports:&[],external_force_n:0.0,destinations:&DESTINATIONS};
    let accepted_model = model.clone(); let accepted_cells = thermal.cells().clone();
    let mut polls = 0;
    let result = conductivity::step_candidate(&thermal,0.1,input,&responses(&cards),config(),
        conductivity_controls(),&[0.0;2],|s|coefficients(s,model.circuit(),&cards,&curves),
        ||{polls+=1;false}).unwrap();
    let mut actual = 0;
    assert!(matches!(conductivity::step_candidate(&thermal,0.1,input,&responses(&cards),config(),
        conductivity_controls(),&[0.0;2],|s|coefficients(s,model.circuit(),&cards,&curves),
        ||{actual+=1;actual==polls}),Err(CircuitEnthalpyError::Cancelled)));
    let limit = fs_couple::circuit_enthalpy::CircuitEnthalpySolve {
        maximum_circuit_trials:result.report().circuit_trials-1,..config() };
    assert!(matches!(conductivity::step_candidate(&thermal,0.1,input,&responses(&cards),limit,
        conductivity_controls(),&[0.0;2],|s|coefficients(s,model.circuit(),&cards,&curves),||false),
        Err(CircuitEnthalpyError::CircuitBudget)));
    assert_eq!(model,accepted_model); assert_eq!(thermal.cells(),&accepted_cells); assert_eq!(thermal.time_s(),0.0);
}
