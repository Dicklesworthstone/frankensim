//! G1/G3: nonlinear material data reaches actual field/body trajectories.
use super::*;
use fs_couple::material_circuit::MaterialCircuitError;
use fs_couple::material_circuit::conductivity::{ConductionAxis, ConductionCoordinate,
    ConductivitySolve, ConductorResponse};
use fs_couple::material_circuit::conductivity::contact::{ContactAxis, ContactCoordinate,
    ContactDirection, ContactResponse};
use fs_qty::QuantitySpec;

fn axis(current: bool) -> ConductionAxis<'static> {
    ConductionAxis { name: "operating", coordinate: if current {
        ConductionCoordinate::CurrentDensityMagnitude
    } else { ConductionCoordinate::ElectricFieldMagnitude } }
}
fn nonlinear(beta: f64, current: bool) -> MaterialCard {
    let dims = if current { Dims([-2,0,0,0,1,0]) } else { Dims([1,1,-3,0,-1,0]) };
    let basis = ElectricalPropertyBasis::Conductivity;
    source(PropertyClaim { key: PropertyKey::new(basis.property_name(), basis.dimensions()),
        value: PropertyValue::Curve { abscissa: "operating".into(), abscissa_dims: dims,
            knots: vec![(0.0,1.0),(4.0,1.0+4.0*beta)], dims: basis.dimensions() },
        validity: ValidityDomain::unconstrained()
            .with_quantity("operating",QuantitySpec::dimensional(dims),0.0,4.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic nonlinear actuator conductor".into(),
            license: "test-only".into(), artifact: None } })
}
fn model(card: &MaterialCard, current: bool, polarity: f64) -> MaterialFieldMotion {
    model_with_storage(card,current,polarity,false)
}
fn model_with_storage(card: &MaterialCard, current: bool, polarity: f64, magnetic: bool)
    -> MaterialFieldMotion {
    let dims = if current { Dims([-2,0,0,0,1,0]) } else { Dims([1,1,-3,0,-1,0]) };
    let point = QueryPoint::new().with_quantity("operating",QuantitySpec::dimensional(dims),0.0).unwrap();
    let resistor = MaterialResistor { branch: 1, conductor: resolve_uniform_conductor_with_basis(
        card,&point,MaterialPropertySelection::SingleClaimOnly,ElectricalPropertyBasis::Conductivity,
        1.0,1.0).unwrap() };
    let element = if magnetic { Branch::Inductor {henries:f64::NAN} }
        else { Branch::Capacitor {farads:f64::NAN} };
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 2,
        branches: vec![(1,0,element),
            (1,0,Branch::Resistor { ohms:f64::NAN })], transformers:vec![] },
        vec![resistor],vec![if magnetic {inductor(0,1.0)} else {capacitor(0,2.0,0.5,1.0)}],
        if magnetic {vec![polarity*0.8]} else {vec![]},
        if magnetic {vec![]} else {vec![polarity*0.8]},&[]).unwrap();
    MaterialFieldMotion::new(circuit,0,plate(10_000.0),0.05,[0.5,1.5]).unwrap()
}
fn solve() -> ConductivitySolve {
    ConductivitySolve { maximum_iterations:64, absolute_resistance_tolerance_ohm:1e-13,
        relative_resistance_tolerance:1e-12 }
}
fn input(dt: f64) -> FieldMotionInput<'static> {
    FieldMotionInput { electrical_ports:&[],external_force_n:0.0,duration_s:dt }
}

#[test]
fn constant_response_is_exactly_the_original_joint_descriptor() {
    let card = nonlinear(0.0,false); let motion = model(&card,false,1.0);
    let axes = [axis(false)]; let bindings = [ConductorResponse {branch:1,card:&card,axes:&axes}];
    let original = motion.step_candidate(motion.circuit().resistors().to_vec(),input(0.1),controls(),||false).unwrap();
    let nonlinear = motion.conductivity_candidate(&bindings,input(0.1),controls(),solve(),||false).unwrap();
    assert_eq!(nonlinear.motion(),original.model());
    assert_eq!(nonlinear.step().field.native.x,original.step().native.x);
    assert_eq!(nonlinear.step().field.heat_j,original.step().heat_j);
    assert_eq!(nonlinear.step().iterations,1);
    assert_eq!(nonlinear.step().source_evaluations,1);

    let interface = contact_source::card(&card,&card,0.4,0.0);
    let circuit = motion.circuit().rebind_contacts_candidate(
        vec![contact_source::contact(7,1,&interface,0.5)],&[]).unwrap();
    let motion = MaterialFieldMotion::new(circuit,0,motion.body().clone(),0.05,[0.5,1.5]).unwrap();
    let original = motion.step_with_contacts_candidate(motion.circuit().resistors().to_vec(),
        motion.circuit().contacts().to_vec(),input(0.1),controls(),||false).unwrap();
    let result = motion.conductivity_with_contacts_candidate(motion.circuit().resistors().to_vec(),
        motion.circuit().contacts().to_vec(),&bindings,input(0.1),controls(),solve(),||false).unwrap();
    assert_eq!(result.motion(),original.model());
    assert_eq!(result.step().field.native.x,original.step().native.x);
    assert_eq!(result.step().field.heat_j,original.step().heat_j);
    assert!(result.step().evaluated_contacts.is_empty());
    assert_eq!(result.step().iterations,1);
}

fn trajectory(beta: f64, current: bool, polarity: f64, dt: f64, steps: usize) -> [f64;3] {
    let card = nonlinear(beta,current); let mut motion = model(&card,current,polarity);
    let axes = [axis(current)]; let bindings = [ConductorResponse {branch:1,card:&card,axes:&axes}];
    let initial = 0.5*0.8*0.8+0.5*0.05*0.05; let mut heat = 0.0;
    for _ in 0..steps {
        let next = motion.conductivity_candidate(&bindings,input(dt),controls(),solve(),||false).unwrap();
        let sample = &next.step().evaluated[0];
        let e = sample.point.electric_field_v_m;
        let j = sample.point.current_density_a_m2;
        let operating = if current { j.abs() } else { e.abs() };
        near(j,(1.0+beta*operating)*e,1e-10);
        near(sample.conductor.material().query().axes()["operating"],operating,1e-14);
        near(sample.power_defect_w,0.0,controls().power_tolerance_w);
        heat += next.step().field.heat_j;
        motion = next.into_parts().0;
    }
    let q = motion.circuit().charges_c()[0]; let d = motion.gap_m(); let v = motion.velocity_m_s();
    near(0.5*q*q*d+0.5*v*v+heat,initial,2e-9);
    [q,d,v]
}

// Independent continuous RC/Maxwell-force ODE. RK4 supplies a reference without
// native descriptors, discrete gradients, material iteration or solver calls.
fn reference(beta: f64, current: bool, polarity: f64) -> [f64;3] {
    let derivative = |x: [f64;3]| {
        let e = x[0]*x[1];
        let i = if current { e/(1.0-beta*e.abs()) } else { (1.0+beta*e.abs())*e };
        [-i,x[2],-0.5*x[0]*x[0]]
    };
    let add = |x: [f64;3], k: [f64;3], dt: f64| std::array::from_fn(|i| x[i]+dt*k[i]);
    let mut x = [polarity*0.8,1.0,0.05]; let dt = 0.0001;
    for _ in 0..2000 {
        let a = derivative(x); let b = derivative(add(x,a,dt*0.5));
        let c = derivative(add(x,b,dt*0.5)); let d = derivative(add(x,c,dt));
        x = std::array::from_fn(|i| x[i]+dt*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0);
    }
    x
}

#[test]
fn field_and_current_laws_change_real_motion_and_converge_to_independent_dynamics() {
    for current in [false,true] {
        for polarity in [-1.0,1.0] {
            let expected = reference(0.25,current,polarity);
            let coarse = trajectory(0.25,current,polarity,0.02,10);
            let fine = trajectory(0.25,current,polarity,0.01,20);
            let error = |x: [f64;3]| (0..3).map(|i|(x[i]-expected[i]).abs()).fold(0.0,f64::max);
            assert!(error(fine)<0.4*error(coarse));
            assert!(error(fine)<2e-5);
            let constant = trajectory(0.0,current,polarity,0.01,20);
            assert!(fine[0].abs()<constant[0].abs());
            assert!(fine[1]>constant[1],"stronger conduction reduces field attraction and closing motion");
        }
    }
}

#[test]
fn current_dependent_conductor_data_also_changes_magnetic_flux_and_attraction() {
    let mut endpoints = Vec::new();
    for beta in [0.0,0.25] {
        let card = nonlinear(beta,true); let mut motion = model_with_storage(&card,true,1.0,true);
        let axes = [axis(true)]; let bindings = [ConductorResponse {branch:1,card:&card,axes:&axes}];
        let mut heat = 0.0;
        for _ in 0..10 {
            let next = motion.conductivity_candidate(&bindings,input(0.02),controls(),solve(),||false).unwrap();
            let sample = &next.step().evaluated[0];
            let j = sample.point.current_density_a_m2;
            near(j,(1.0+beta*j.abs())*sample.point.electric_field_v_m,1e-10);
            heat += next.step().field.heat_j;
            motion = next.into_parts().0;
        }
        let flux = motion.circuit().fluxes_wb()[0];
        near(0.5*flux*flux*motion.gap_m()+0.5*motion.velocity_m_s().powi(2)+heat,
            0.5*0.8*0.8+0.5*0.05*0.05,2e-9);
        endpoints.push((flux,motion.gap_m()));
    }
    assert!(endpoints[1].0>endpoints[0].0,"reduced conductor loss retains magnetic linkage");
    assert!(endpoints[1].1<endpoints[0].1,"retained linkage strengthens magnetic attraction");
}

#[test]
fn nonlinear_failure_and_final_cancellation_preserve_the_complete_original_motion() {
    let card = nonlinear(0.25,false); let motion = model(&card,false,1.0); let saved = motion.clone();
    let axes = [axis(false)]; let bindings = [ConductorResponse {branch:1,card:&card,axes:&axes}];
    let limited = ConductivitySolve { maximum_iterations:1,..solve() };
    assert!(matches!(motion.conductivity_candidate(&bindings,input(0.1),controls(),limited,||false),
        Err(MaterialCircuitError::ConductivityConvergence { iterations:1,.. })));
    let mut polls = 0;
    let accepted = motion.conductivity_candidate(&bindings,input(0.1),controls(),solve(),||{polls+=1;false}).unwrap();
    let mut actual = 0;
    assert!(matches!(motion.conductivity_candidate(&bindings,input(0.1),controls(),solve(),||{
        actual+=1;actual==polls }),Err(MaterialCircuitError::Cancelled)));
    assert_eq!(motion,saved);
    let retry = motion.conductivity_candidate(&bindings,input(0.1),controls(),solve(),||false).unwrap();
    assert_eq!(retry.motion(),accepted.motion());
}

const CONTACT_V: Dims = Dims([2,1,-3,0,-1,0]);
const CONTACT_J: Dims = Dims([-2,0,0,0,1,0]);

// Independent continuous circuit/gap equations, with nonlinear bulk voltage
// V_b + 0.2 |V_b| V_b = I and two distinct source-backed series interfaces.
fn contact_voltages(current: f64, beta: f64) -> [f64;3] {
    let i = current.abs();
    [2.0*current/(1.0+(1.0+0.8*i).sqrt()),
        current*0.3/(1.0-0.3*beta*i),
        current*0.8*(1.0+beta*i/0.5)]
}
fn contact_reference(beta: f64, magnetic: bool, sign: f64) -> [f64;3] {
    let derivative = |x: [f64;3]| {
        let drive = x[0]*x[1];
        let decay = if magnetic { contact_voltages(drive,beta).iter().sum::<f64>() }
        else {
            let (mut lo,mut hi) = (0.0,4.0);
            for _ in 0..70 {
                let i = 0.5*(lo+hi);
                if contact_voltages(i,beta).iter().sum::<f64>()>drive.abs() { hi=i; }
                else { lo=i; }
            }
            drive.signum()*0.5*(lo+hi)
        };
        [-decay,x[2],-0.5*x[0]*x[0]]
    };
    let add = |x: [f64;3],k: [f64;3],h: f64| std::array::from_fn(|i|x[i]+h*k[i]);
    let mut x = [sign*0.8,1.0,0.05]; let h = 0.0001;
    for _ in 0..2000 {
        let a=derivative(x); let b=derivative(add(x,a,h*0.5));
        let c=derivative(add(x,b,h*0.5)); let d=derivative(add(x,c,h));
        x=std::array::from_fn(|i|x[i]+h*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0);
    }
    x
}
fn contact_trajectory(beta: f64, magnetic: bool, sign: f64, steps: usize) -> [f64;3] {
    let card = nonlinear(0.2,false);
    let interfaces = [contact_source::response_card(&card,&card,"V",CONTACT_V,0.3,0.0,beta,4.0),
        contact_source::response_card(&card,&card,"J",CONTACT_J,0.4,0.0,beta,4.0)];
    let seed = model_with_storage(&card,false,sign,magnetic);
    let circuit = seed.circuit().rebind_contacts_candidate(vec![
        contact_source::response_contact(9,1,&interfaces[1],0.5,"J",CONTACT_J),
        contact_source::response_contact(7,1,&interfaces[0],1.0,"V",CONTACT_V)],&[]).unwrap();
    let mut motion = MaterialFieldMotion::new(circuit,0,seed.body().clone(),0.05,[0.5,1.5]).unwrap();
    let axes = [axis(false)]; let bulk = [ConductorResponse {branch:1,card:&card,axes:&axes}];
    let voltage = [ContactAxis {name:"V",coordinate:ContactCoordinate::VoltageMagnitude}];
    let density = [ContactAxis {name:"J",coordinate:ContactCoordinate::CurrentDensityMagnitude}];
    let contacts = [ContactResponse {id:9,card:&interfaces[1],direction:ContactDirection::AlongBranch,axes:&density},
        ContactResponse {id:7,card:&interfaces[0],direction:ContactDirection::AlongBranch,axes:&voltage}];
    let dt = 0.2/steps as f64; let mut heat = 0.0;
    for _ in 0..steps {
        let next = motion.conductivity_with_contact_response_candidate(
            motion.circuit().resistors().to_vec(),motion.circuit().contacts().to_vec(),
            &bulk,Some(&contacts),input(dt),controls(),solve(),||false).unwrap();
        let report = next.step(); let branch = &report.field.resistors[0];
        let expected = contact_voltages(branch.current_a,beta);
        near(branch.voltage_v,expected[0],2e-10);
        for (slot,sample) in report.evaluated_contacts.iter().enumerate() {
            let native = &report.field.contacts[slot];
            near(native.voltage_v,expected[slot+1],2e-10);
            near(sample.point.current_a,branch.current_a,2e-12);
            near(sample.point.current_density_a_m2,branch.current_a/native.contact.area_m2(),2e-12);
            let (axis,coordinate) = if slot==0 { ("V",native.voltage_v.abs()) }
                else { ("J",(branch.current_a/0.5).abs()) };
            near(sample.contact.source().resolved().state_coordinate(axis).unwrap(),coordinate,2e-12);
            near(native.heat_j,dt*branch.current_a*expected[slot+1],2e-10);
            assert_eq!(sample.contact.source().resolved().card_identity(),interfaces[slot].content_hash());
        }
        near(report.field.heat_j,branch.heat_j+report.field.contacts.iter().map(|h|h.heat_j).sum::<f64>(),2e-10);
        assert_eq!(report.source_evaluations,3*report.iterations);
        heat += report.field.heat_j;
        motion = next.into_parts().0;
    }
    let coordinate = if magnetic {motion.circuit().fluxes_wb()[0]} else {motion.circuit().charges_c()[0]};
    let result = [coordinate,motion.gap_m(),motion.velocity_m_s()];
    near(0.5*coordinate*coordinate*result[1]+0.5*result[2]*result[2]+heat,
        0.5*0.8*0.8+0.5*0.05*0.05,3e-9);
    result
}

#[test]
fn g1_contact_material_changes_reciprocal_electric_and_magnetic_motion() {
    for magnetic in [false,true] {
        for sign in [-1.0,1.0] {
            let expected = contact_reference(0.08,magnetic,sign);
            let coarse = contact_trajectory(0.08,magnetic,sign,10);
            let fine = contact_trajectory(0.08,magnetic,sign,20);
            let error = |x: [f64;3]| (0..3).map(|i|(x[i]-expected[i]).abs()).fold(0.0,f64::max);
            assert!(error(fine)<0.4*error(coarse));
            assert!(error(fine)<3e-5);
            let fixed = contact_trajectory(0.0,magnetic,sign,20);
            assert!((fine[0]-fixed[0]).abs()>1e-4);
            assert!((fine[1]-fixed[1]).abs()>1e-6,
                "source-only contact changes reach actual gap motion, not only diagnostic current");
        }
    }
}

#[test]
fn g4_contact_motion_source_refusal_and_late_cancel_preserve_joint_history() {
    let card = nonlinear(0.2,false); let seed = model(&card,false,1.0);
    let axes = [axis(false)]; let bulk = [ConductorResponse {branch:1,card:&card,axes:&axes}];
    let voltage = [ContactAxis {name:"V",coordinate:ContactCoordinate::VoltageMagnitude}];
    for limit in [0.01,4.0] {
        let interface = contact_source::response_card(&card,&card,"V",CONTACT_V,0.4,0.0,0.08,limit);
        let circuit = seed.circuit().rebind_contacts_candidate(vec![
            contact_source::response_contact(7,1,&interface,1.0,"V",CONTACT_V)],&[]).unwrap();
        let motion = MaterialFieldMotion::new(circuit,0,seed.body().clone(),0.05,[0.5,1.5]).unwrap();
        let saved = motion.clone();
        let contacts = [ContactResponse {id:7,card:&interface,direction:ContactDirection::AlongBranch,axes:&voltage}];
        let run = |solve,check: &mut dyn FnMut()->bool| motion.conductivity_with_contact_response_candidate(
            motion.circuit().resistors().to_vec(),motion.circuit().contacts().to_vec(),
            &bulk,Some(&contacts),input(0.1),controls(),solve,check);
        if limit<1.0 { assert!(run(solve(),&mut ||false).is_err()); }
        else {
            assert!(motion.conductivity_with_contacts_candidate(motion.circuit().resistors().to_vec(),
                motion.circuit().contacts().to_vec(),&bulk,input(0.1),controls(),solve(),||false).is_err());
            assert!(matches!(run(ConductivitySolve {maximum_iterations:1,..solve()},&mut ||false),
                Err(MaterialCircuitError::ConductivityConvergence {..})));
            let mut polls=0; let accepted=run(solve(),&mut ||{polls+=1;false}).unwrap();
            let mut actual=0;
            assert!(matches!(run(solve(),&mut ||{actual+=1;actual==polls}),Err(MaterialCircuitError::Cancelled)));
            assert_eq!(run(solve(),&mut ||false).unwrap().motion(),accepted.motion());
        }
        assert_eq!(motion,saved);
    }
}
