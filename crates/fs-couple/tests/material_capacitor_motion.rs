//! Data-only material substitutions through reciprocal field/circuit/body motion.
#[path = "electrical_contacts/source.rs"]
mod contact_source;
use fs_couple::material_circuit::motion::{FieldMotionInput, MaterialFieldMotion};
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve, MaterialResistor};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::conductor::{ElectricalPropertyBasis, resolve_uniform_conductor_with_basis};
use fs_material::field_storage::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, Dims, Temperature};

fn card(property: &str, dims: Dims, value: f64) -> MaterialCard {
    source(PropertyClaim {
        key: PropertyKey::new(property, dims), value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(), uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
        provenance: Provenance { source: "synthetic electromechanical fixture".into(),
            license: "test-only".into(), artifact: None },
    })
}
fn source(claim: PropertyClaim) -> MaterialCard {
    let mut claims = ClaimSet::new(); claims.insert_claim(claim).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "authored scalar material".into(),
        phase: "solid".into(), process: "numerical fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn plate(density: f64) -> ResolvedBody {
    let material = card("density", Density::DIMS, density);
    resolve_body(&[RegionBinding { key: "electrode", material: &material,
        point: &QueryPoint::new(), selection: MaterialPropertySelection::SingleClaimOnly,
        geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.01)).unwrap(),
    }], 1, || false).unwrap()
}
fn capacitor(branch: usize, epsilon: f64, area: f64, gap: f64) -> MaterialStorage {
    let material = card(ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS, epsilon);
    MaterialStorage { branch, element: MaterialStorageElement::Capacitor(
        ResolvedCapacitor::resolve(&material, &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly, area, gap).unwrap()) }
}
fn inductor(branch: usize, inductance: f64) -> MaterialStorage {
    let material = card(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, inductance);
    MaterialStorage { branch, element: MaterialStorageElement::Inductor(
        ResolvedInductor::resolve(&material, &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly, 1.0, 1.0, 1).unwrap()) }
}
fn resistor(branch: usize, resistance: f64) -> MaterialResistor {
    // The actual consumer uses the shared conductivity basis, not a resistor literal.
    let basis = ElectricalPropertyBasis::Conductivity;
    let material = card(basis.property_name(), basis.dimensions(), 1.0/resistance);
    MaterialResistor { branch, conductor: resolve_uniform_conductor_with_basis(
        &material, &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly,
        basis, 1.0, 1.0).unwrap() }
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve { energy_tolerance_j: 1e-10, power_tolerance_w: 1e-10 }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a} != {b}, tolerance {tolerance}");
}
fn isolated(epsilon: f64, charge: f64, density: f64, limits: [f64; 2]) -> MaterialFieldMotion {
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 2,
        branches: vec![(1, 0, Branch::Capacitor { farads: f64::NAN })], transformers: vec![] },
        vec![], vec![capacitor(0, epsilon, 0.5, 1.0)], vec![], vec![charge], &[]).unwrap();
    MaterialFieldMotion::new(circuit, 0, plate(density), 0.02, limits).unwrap()
}

#[test]
fn dielectric_and_electrode_material_swaps_change_real_motion_by_the_energy_derivatives() {
    let (dt, external, initial_velocity) = (0.25, 0.01, 0.02);
    for q in [-0.2_f64, 0.0, 0.2] {
        for (epsilon, density) in [(1.0, 1000.0), (2.0, 1000.0), (1.0, 8000.0)] {
            let model = isolated(epsilon, q, density, [0.5, 1.5]);
            let before = model.clone();
            let mass = density*0.1*0.1*0.01;
            let attraction = -q*q/(2.0*epsilon*0.5);
            let acceleration = (external+attraction)/mass;
            near(model.body().inertia().mass(), mass, 1e-14);
            near(model.field_force_n().unwrap(), attraction, 1e-14);
            let candidate = model.step_candidate(vec![], FieldMotionInput {
                electrical_ports: &[], external_force_n: external, duration_s: dt,
            }, controls(), || false).unwrap();
            let next = candidate.model();
            let step = candidate.step();
            // Isolated charge gives constant force. The independent continuous
            // solution is exact for this midpoint step, including q=0 ballistic motion.
            let gap = 1.0+initial_velocity*dt+0.5*acceleration*dt*dt;
            let velocity = initial_velocity+acceleration*dt;
            near(next.gap_m(), gap, 2e-11);
            near(next.velocity_m_s(), velocity, 2e-11);
            near(next.circuit().charges_c()[0], q, 1e-13);
            near(next.capacitor().unwrap().capacitance_f(), epsilon*0.5/gap, 1e-12);
            near(step.mechanical_supply_j, external*(gap-1.0), 1e-12);
            near(step.field_work_on_body_j, attraction*(gap-1.0), 1e-12);
            near(step.electrical_energy_j[1], q*q*gap/(2.0*epsilon*0.5), 1e-12);
            near(step.kinetic_energy_j[1], 0.5*mass*velocity*velocity, 1e-12);
            near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
            assert_eq!(step.heat_j, 0.0);
            assert_eq!(step.electrical_supply_j, 0.0);
            assert_eq!(next.body(), model.body());
            assert_eq!(next.capacitor().unwrap().material(), model.capacitor().unwrap().material());
            assert_eq!(model, before);
        }
    }
}

fn network() -> MaterialCircuit {
    MaterialCircuit::new_with_storage(CircuitGraph { node_count: 3, branches: vec![
        (1, 2, Branch::Resistor { ohms: f64::NAN }),
        (1, 0, Branch::Capacitor { farads: f64::NAN }),
        (2, 0, Branch::Inductor { henries: f64::NAN }),
        (2, 0, Branch::Capacitor { farads: f64::NAN }),
        (1, 0, Branch::VoltageSource { port: 0 }),
        (2, 0, Branch::CurrentSource { port: 1 }),
    ], transformers: vec![] }, vec![resistor(0, 2.0)],
        // Unsorted bindings and the SECOND capacitor exercise original branch
        // ordering rather than a special circuit layout or active-index shortcut.
        vec![capacitor(3, 2.0, 0.5, 1.0), inductor(2, 0.5), capacitor(1, 0.5, 1.0, 1.0)],
        vec![0.1], vec![0.2, 0.5], &[0.4, 0.05]).unwrap()
}
fn sourced_electrical_energy(circuit: &MaterialCircuit) -> f64 {
    let (mut charge, mut flux) = (0, 0);
    circuit.storage().iter().map(|source| match &source.element {
        MaterialStorageElement::Capacitor(c) => {
            let energy = c.stored_energy_j(circuit.charges_c()[charge]).unwrap();
            charge += 1; energy
        }
        MaterialStorageElement::Inductor(l) => {
            let energy = l.stored_energy_j(circuit.fluxes_wb()[flux]).unwrap();
            flux += 1; energy
        }
    }).sum()
}

#[test]
fn moving_gap_feeds_back_through_the_driven_rlc_network_and_separate_contact_heat() {
    for contact_resistance in [None, Some(0.75)] {
        let contacts = if let Some(resistance) = contact_resistance {
            let basis = ElectricalPropertyBasis::Conductivity;
            let a = card(basis.property_name(), basis.dimensions(), 0.5);
            let b = card(basis.property_name(), basis.dimensions(), 1.0);
            let source = contact_source::card(&a, &b, resistance, 0.0);
            vec![contact_source::contact(7, 0, &source, 1.0)]
        } else { vec![] };
        let circuit = network().rebind_contacts_candidate(contacts.clone(), &[0.4, 0.05]).unwrap();
        let mut held = circuit.clone();
        let mut moving = MaterialFieldMotion::new(circuit, 3, plate(5000.0), 0.0, [0.5, 1.5]).unwrap();
        let body = moving.body().clone();
        let dielectric = moving.capacitor().unwrap().material().clone();
        let fixed_storage = moving.circuit().storage()[..2].to_vec();
        let initial_total = sourced_electrical_energy(moving.circuit());
        let mut supplied = 0.0;
        let mut heat = 0.0;
        for _ in 0..20 {
            let previous_gap = moving.gap_m();
            let candidate = moving.step_with_contacts_candidate(vec![resistor(0, 2.0)], contacts.clone(),
                FieldMotionInput { electrical_ports: &[0.4, 0.05], external_force_n: 0.03,
                    duration_s: 0.02 }, controls(), || false).unwrap();
            let (next, step) = candidate.into_parts();
            assert_eq!(step.native.y.len(), 3);
            near(step.electrical_supply_j, 0.02*(0.4*step.native.y[0]+0.05*step.native.y[1]), 1e-13);
            near(step.mechanical_supply_j, 0.03*(next.gap_m()-previous_gap), 1e-13);
            near(step.electrical_energy_j[1], sourced_electrical_energy(next.circuit()), 1e-12);
            near(step.kinetic_energy_j[1], 0.5*body.inertia().mass()*next.velocity_m_s().powi(2), 1e-12);
            assert_eq!(step.resistors.len(), 1);
            assert_eq!(step.resistors[0].material.branch, 0);
            let current = step.resistors[0].current_a;
            near(step.resistors[0].heat_j, 0.02*current*current*2.0, 1e-12);
            assert_eq!(step.contacts.len(), contacts.len());
            if let Some(resistance) = contact_resistance {
                near(step.contacts[0].heat_j, 0.02*current*current*resistance, 1e-12);
                near(step.contacts[0].current_a, current, 1e-12);
                assert_eq!(&step.contacts[0].contact, &contacts[0]);
            }
            near(step.heat_j, 0.02*current*current*(2.0+contact_resistance.unwrap_or(0.0)), 1e-12);
            assert!(step.heat_j > 0.0);
            supplied += step.electrical_supply_j+step.mechanical_supply_j;
            heat += step.heat_j;
            near(step.electrical_energy_j[1]+step.kinetic_energy_j[1]+heat,
                initial_total+supplied, 2e-9);
            near(step.field_work_on_body_j, step.electrical_supply_j-step.heat_j
                -(step.electrical_energy_j[1]-step.electrical_energy_j[0]), 1e-10);
            assert_eq!(next.body(), &body);
            assert_eq!(next.capacitor().unwrap().material(), &dielectric);
            assert_eq!(next.circuit().contacts(), contacts.as_slice());
            assert_eq!(&next.circuit().storage()[..2], fixed_storage.as_slice());
            assert_eq!(next.circuit().graph().branches[3].2,
                Branch::Capacitor { farads: next.capacitor().unwrap().capacitance_f() });
            near(next.circuit().charges_c()[0], 0.2, 1e-12);
            held = held.step_with_contacts_candidate(vec![resistor(0, 2.0)], contacts.clone(),
                &[0.4, 0.05], 0.02, controls(), || false).unwrap().into_parts().0;
            moving = next;
        }
        assert!(moving.gap_m() < 1.0);
        assert!(moving.capacitor().unwrap().capacitance_f() > 1.0);
        assert!((moving.circuit().charges_c()[1]-held.charges_c()[1]).abs() > 1e-6,
            "geometry must change actual circuit charge, not just the force report");
        assert!((moving.circuit().fluxes_wb()[0]-held.fluxes_wb()[0]).abs() > 1e-6);
    }
}

#[test]
fn small_si_gap_with_high_voltage_does_not_inherit_the_voltage_probe_scale() {
    let (epsilon, area, gap, voltage, dt) = (8.8541878128e-12, 0.01, 1e-4, 100.0, 1e-3);
    let charge = epsilon*area/gap*voltage;
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 3, branches: vec![
        (1, 0, Branch::VoltageSource { port: 0 }),
        (2, 0, Branch::Capacitor { farads: f64::NAN }),
    ], transformers: vec![] }, vec![], vec![capacitor(1, epsilon, area, gap)],
        vec![], vec![charge], &[voltage]).unwrap();
    // The separate voltage branch keeps the source in the same descriptor but
    // supplies no current. Charge stays isolated, giving an independent oracle.
    let moving = MaterialFieldMotion::new(circuit, 1, plate(1000.0), 0.0, [5e-5, 2e-4]).unwrap();
    let mass = moving.body().inertia().mass();
    let force = -charge*charge/(2.0*epsilon*area);
    let candidate = moving.step_candidate(vec![], FieldMotionInput {
        electrical_ports: &[voltage], external_force_n: 0.0, duration_s: dt,
    }, MaterialCircuitSolve { energy_tolerance_j: 1e-12, power_tolerance_w: 1e-12 }, || false).unwrap();
    near(candidate.model().gap_m(), gap+0.5*force/mass*dt*dt, 1e-11);
    near(candidate.model().momentum_kg_m_s(), force*dt, 1e-11);
    near(candidate.model().circuit().charges_c()[0], charge, 1e-15);
    near(candidate.step().electrical_supply_j, 0.0, 1e-12);
    near(candidate.step().field_work_on_body_j, candidate.step().kinetic_energy_j[1], 1e-12);
    near(candidate.step().energy_residual_j, 0.0, 1e-12);
}

#[test]
fn unsupported_gap_ideal_clamp_or_cancellation_cannot_partially_advance_the_model() {
    let bounded = isolated(1.0, 0.2, 1000.0, [0.99, 1.5]);
    let before = bounded.clone();
    assert!(bounded.step_candidate(vec![], FieldMotionInput {
        electrical_ports: &[], external_force_n: 0.0, duration_s: 1.0,
    }, controls(), || false).is_err());
    assert_eq!(bounded, before);
    let graph = CircuitGraph { node_count: 2, branches: vec![
        (1, 0, Branch::Capacitor { farads: f64::NAN }),
        (1, 0, Branch::VoltageSource { port: 0 }),
    ], transformers: vec![] };
    let circuit = MaterialCircuit::new_with_storage(graph, vec![], vec![capacitor(0, 1.0, 0.5, 1.0)],
        vec![], vec![0.2], &[0.4]).unwrap();
    let clamped = MaterialFieldMotion::new(circuit, 0, plate(1000.0), 0.0, [0.5, 1.5]).unwrap();
    let before = clamped.clone();
    assert!(clamped.step_candidate(vec![], FieldMotionInput {
        electrical_ports: &[0.4], external_force_n: 0.0, duration_s: 0.1,
    }, controls(), || false).is_err());
    assert_eq!(clamped, before);

    let model = isolated(1.0, 0.2, 1000.0, [0.5, 1.5]);
    let before = model.clone();
    let input = FieldMotionInput { electrical_ports: &[], external_force_n: 0.0, duration_s: 0.1 };
    let mut polls = 0;
    model.step_candidate(vec![], input, controls(), || { polls += 1; false }).unwrap();
    let mut count = 0;
    assert!(model.step_candidate(vec![], input, controls(), || {
        count += 1; count == polls
    }).is_err());
    assert_eq!(model, before);
}

#[test]
fn temperature_dependent_dielectric_cannot_invent_missing_reversible_heat() {
    let material = source(PropertyClaim {
        key: PropertyKey::new(ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 1.0), (400.0, 2.0)], dims: ELECTRIC_PERMITTIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic thermal dielectric".into(),
            license: "test-only".into(), artifact: None },
    });
    let cap = ResolvedCapacitor::resolve(&material, &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly, 0.5, 1.0).unwrap();
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 2,
        branches: vec![(1, 0, Branch::Capacitor { farads: f64::NAN })], transformers: vec![] },
        vec![], vec![MaterialStorage { branch: 0, element: MaterialStorageElement::Capacitor(cap) }],
        vec![], vec![0.2], &[]).unwrap();
    assert!(MaterialFieldMotion::new(circuit, 0, plate(1000.0), 0.0, [0.5, 1.5]).is_err());
}
