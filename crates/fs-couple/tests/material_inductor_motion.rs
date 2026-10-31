//! G1/G3: magnetic material and mass data drive the same field/body transaction.
//! The independent gap-pressure reference is H=lambda² d/(2 mu N² A),
//! F=-B² A/(2 mu), B=lambda/(N A), with an ideal magnetic return path.
//! MIT 6.013 Lecture 6, magnetic energy/pressure and fixed-linkage virtual work:
//! https://ocw.mit.edu/courses/6-013-electromagnetics-and-applications-spring-2009/89a8f6c53b4f4b62094ebbe821cfa436_MIT6_013S09_lec06.pdf

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
use fs_material::field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, ResolvedCapacitor, ResolvedInductor};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, Dims, Temperature};

fn card(property: &str, dims: Dims, value: f64) -> MaterialCard {
    source(PropertyClaim { key: PropertyKey::new(property, dims),
        value: PropertyValue::Scalar { value, dims }, validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "synthetic magnetic gap fixture".into(),
            license: "test-only".into(), artifact: None } })
}
fn source(claim: PropertyClaim) -> MaterialCard {
    let mut claims = ClaimSet::new(); claims.insert_claim(claim).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "numerical field/mass material".into(),
        phase: "declared fixture".into(), process: "numerical fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn body(density: f64) -> ResolvedBody {
    let material = card("density", Density::DIMS, density);
    resolve_body(&[RegionBinding { key: "guided pole", material: &material, point: &QueryPoint::new(),
        selection: MaterialPropertySelection::SingleClaimOnly,
        geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.01)).unwrap() }], 1, || false).unwrap()
}
fn magnetic(branch: usize, mu: f64, area: f64, gap: f64, turns: u32) -> MaterialStorage {
    let material = card(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, mu);
    MaterialStorage { branch, element: MaterialStorageElement::Inductor(ResolvedInductor::resolve(
        &material, &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly, area, gap, turns).unwrap()) }
}
fn dielectric(branch: usize, capacitance: f64) -> MaterialStorage {
    let material = card(ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS, capacitance);
    MaterialStorage { branch, element: MaterialStorageElement::Capacitor(ResolvedCapacitor::resolve(
        &material, &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly, 1.0, 1.0).unwrap()) }
}
fn resistor(branch: usize, resistance: f64) -> MaterialResistor {
    let basis = ElectricalPropertyBasis::Conductivity;
    let material = card(basis.property_name(), basis.dimensions(), 1.0/resistance);
    MaterialResistor { branch, conductor: resolve_uniform_conductor_with_basis(&material,
        &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly, basis, 1.0, 1.0).unwrap() }
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve { energy_tolerance_j: 1e-10, power_tolerance_w: 1e-10 }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual-expected).abs() <= tolerance, "{actual} != {expected}, tolerance {tolerance}");
}
fn shorted(source: MaterialStorage, flux: f64) -> MaterialCircuit {
    MaterialCircuit::new_with_storage(CircuitGraph { node_count: 2, branches: vec![
        (1, 0, Branch::Inductor { henries: f64::NAN }),
        (1, 0, Branch::VoltageSource { port: 0 })], transformers: vec![] },
        vec![], vec![source], vec![flux], vec![], &[0.0]).unwrap()
}

#[test]
fn moving_a_shared_winding_requires_a_shared_geometry_law() {
    use fs_couple::material_circuit::MaterialCircuitError;
    use fs_phs::circuit::MutualInductance;
    let circuit=MaterialCircuit::new_with_storage_and_mutual_inductance(
        CircuitGraph {node_count:3,branches:vec![
            (1,0,Branch::Inductor {henries:f64::NAN}),
            (2,0,Branch::Inductor {henries:f64::NAN}),
            (1,0,Branch::VoltageSource {port:0}),
            (2,0,Branch::VoltageSource {port:1}),
        ],transformers:vec![]},vec![],
        vec![magnetic(0,2.0,1.0,1.0,1),magnetic(1,3.0,1.0,1.0,1)],
        vec![MutualInductance {branches:(0,1),henries:1.0}],vec![0.2,0.1],vec![],&[0.0;2]).unwrap();
    assert!(matches!(MaterialFieldMotion::new(circuit,0,body(1000.0),0.0,[0.5,1.5]),
        Err(MaterialCircuitError::Invalid(
            "moving a mutually coupled winding requires a shared geometry-dependent inductance law"))));
}

#[test]
fn independent_gap_motion_retains_other_windings_mutual_energy() {
    use fs_phs::circuit::MutualInductance;
    let links=vec![MutualInductance {branches:(1,2),henries:1.0}];
    let circuit=MaterialCircuit::new_with_storage_and_mutual_inductance(
        CircuitGraph {node_count:4,branches:vec![
            (1,0,Branch::Inductor {henries:f64::NAN}),
            (2,0,Branch::Inductor {henries:f64::NAN}),
            (3,0,Branch::Inductor {henries:f64::NAN}),
            (1,0,Branch::VoltageSource {port:0}),
            (2,0,Branch::VoltageSource {port:1}),
            (3,0,Branch::VoltageSource {port:2}),
        ],transformers:vec![]},vec![],
        vec![magnetic(0,3.0,1.0,1.0,1),magnetic(1,2.0,1.0,1.0,1),
            magnetic(2,4.0,1.0,1.0,1)],links.clone(),vec![0.2,0.3,0.4],vec![],&[0.0;3]).unwrap();
    let moving=MaterialFieldMotion::new(circuit,0,body(1000.0),0.0,[0.5,1.5]).unwrap();
    let next=moving.step_candidate(vec![],FieldMotionInput {electrical_ports:&[0.0;3],
        external_force_n:0.0,duration_s:0.1},controls(),||false).unwrap();
    let gap=1.0-0.5*(0.5*0.2*0.2/3.0)/0.1*0.1*0.1;
    let shared=0.5*(4.0*0.3*0.3-2.0*0.3*0.4+2.0*0.4*0.4)/7.0;
    near(next.model().gap_m(),gap,1e-12);
    near(next.step().electrical_energy_j[0],shared+0.5*0.2*0.2/3.0,1e-12);
    near(next.step().electrical_energy_j[1],shared+0.5*0.2*0.2*gap/3.0,1e-12);
    assert_eq!(next.model().circuit().mutual_inductance(),links);
    for (actual,expected) in next.model().circuit().fluxes_wb().iter().zip([0.2,0.3,0.4]) {
        near(*actual,expected,1e-12);
    }
}

#[test]
fn permeability_turns_and_body_density_change_motion_at_fixed_flux_linkage() {
    let (area, initial_gap, initial_velocity, dt, external) = (0.5, 1.0, 0.02, 0.25, 0.01);
    for flux in [-0.2_f64, 0.0, 0.2] {
        for (mu, turns, density) in [(1.0, 2, 1000.0), (2.0, 2, 1000.0),
            (1.0, 4, 1000.0), (1.0, 2, 8000.0)] {
            // Zero terminal voltage conserves flux linkage, with no electrical work.
            let circuit = shorted(magnetic(0, mu, area, initial_gap, turns), flux);
            let moving = MaterialFieldMotion::new(circuit, 0, body(density), initial_velocity, [0.5, 1.5]).unwrap();
            let original = moving.clone();
            let b = flux/(f64::from(turns)*area);
            let force = -b*b*area/(2.0*mu);
            let mass = density*0.1*0.1*0.01;
            let acceleration = (force+external)/mass;
            let gap = initial_gap+initial_velocity*dt+0.5*acceleration*dt*dt;
            let velocity = initial_velocity+acceleration*dt;
            near(moving.field_force_n().unwrap(), force, 1e-14);
            let candidate = moving.step_candidate(vec![], FieldMotionInput { electrical_ports: &[0.0],
                external_force_n: external, duration_s: dt }, controls(), || false).unwrap();
            let next = candidate.model(); let step = candidate.step();
            near(next.gap_m(), gap, 2e-11);
            near(next.velocity_m_s(), velocity, 2e-11);
            near(next.circuit().fluxes_wb()[0], flux, 1e-13);
            near(next.inductor().unwrap().inductance_h(), mu*f64::from(turns).powi(2)*area/gap, 1e-12);
            near(step.field_work_on_body_j, force*(gap-initial_gap), 1e-12);
            near(step.mechanical_supply_j, external*(gap-initial_gap), 1e-12);
            near(step.electrical_energy_j[1], 0.5*b*b/mu*area*gap, 1e-12);
            near(step.kinetic_energy_j[1], 0.5*mass*velocity*velocity, 1e-12);
            near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
            assert_eq!(step.heat_j, 0.0); assert_eq!(step.electrical_supply_j, 0.0);
            assert!(next.capacitor().is_none());
            assert_eq!(next.inductor().unwrap().material(), moving.inductor().unwrap().material());
            assert_eq!(next.body(), moving.body()); assert_eq!(moving, original);
        }
    }
}

fn network() -> MaterialCircuit {
    MaterialCircuit::new_with_storage(CircuitGraph { node_count: 4, branches: vec![
        (1, 2, Branch::Resistor { ohms: f64::NAN }),
        (3, 0, Branch::Inductor { henries: f64::NAN }),
        (3, 0, Branch::Capacitor { farads: f64::NAN }),
        (2, 0, Branch::Inductor { henries: f64::NAN }),
        (1, 0, Branch::VoltageSource { port: 0 }),
        (2, 0, Branch::Capacitor { farads: f64::NAN })], transformers: vec![] },
        vec![resistor(0, 2.0)],
        // The moving branch is the second inductor; bindings arrive unsorted.
        vec![dielectric(5, 0.3), magnetic(3, 2.0, 0.5, 1.0, 1),
            dielectric(2, 0.4), magnetic(1, 0.5, 1.0, 1.0, 1)],
        vec![0.02, 0.3], vec![0.1, 0.1], &[0.4]).unwrap()
}
fn electrical_energy(circuit: &MaterialCircuit) -> f64 {
    let (mut q, mut flux) = (0, 0);
    circuit.storage().iter().map(|source| match &source.element {
        MaterialStorageElement::Capacitor(capacitor) => {
            let result = capacitor.stored_energy_j(circuit.charges_c()[q]).unwrap(); q += 1; result
        }
        MaterialStorageElement::Inductor(inductor) => {
            let result = inductor.stored_energy_j(circuit.fluxes_wb()[flux]).unwrap(); flux += 1; result
        }
    }).sum()
}

#[test]
fn magnetic_motion_feeds_back_through_source_bound_rlc_and_coil_contact_losses() {
    for contact_resistance in [None, Some(0.75)] {
        let contacts = if let Some(resistance) = contact_resistance {
            let basis = ElectricalPropertyBasis::Conductivity;
            let a = card(basis.property_name(), basis.dimensions(), 0.5);
            let b = card(basis.property_name(), basis.dimensions(), 1.0);
            let source = contact_source::card(&a, &b, resistance, 0.0);
            vec![contact_source::contact(7, 0, &source, 1.0)]
        } else { vec![] };
        let circuit = network().rebind_contacts_candidate(contacts.clone(), &[0.4]).unwrap();
        let initial = electrical_energy(&circuit);
        let mut held = circuit.clone();
        let mut moving = MaterialFieldMotion::new(circuit, 3, body(5000.0), 0.0, [0.5, 1.5]).unwrap();
        let source = moving.inductor().unwrap().material().clone();
        let unchanged_storage = moving.circuit().storage()[..2].to_vec();
        let (mut supplied, mut dissipated) = (0.0, 0.0);
        for _ in 0..20 {
            let gap = moving.gap_m();
            let (next, step) = moving.step_with_contacts_candidate(vec![resistor(0, 2.0)], contacts.clone(),
                FieldMotionInput { electrical_ports: &[0.4], external_force_n: 0.01, duration_s: 0.02 },
                controls(), || false).unwrap().into_parts();
            let current = step.resistors[0].current_a;
            near(step.heat_j, 0.02*current*current*(2.0+contact_resistance.unwrap_or(0.0)), 1e-12);
            near(step.resistors[0].heat_j, 0.02*current*current*2.0, 1e-12);
            if let Some(resistance) = contact_resistance {
                near(step.contacts[0].heat_j, 0.02*current*current*resistance, 1e-12);
                assert_eq!(&step.contacts[0].contact, &contacts[0]);
            }
            near(step.electrical_energy_j[1], electrical_energy(next.circuit()), 1e-12);
            near(step.mechanical_supply_j, 0.01*(next.gap_m()-gap), 1e-12);
            supplied += step.electrical_supply_j+step.mechanical_supply_j; dissipated += step.heat_j;
            near(step.electrical_energy_j[1]+step.kinetic_energy_j[1]+dissipated, initial+supplied, 2e-9);
            near(step.field_work_on_body_j, step.electrical_supply_j-step.heat_j
                -(step.electrical_energy_j[1]-step.electrical_energy_j[0]), 1e-10);
            assert_eq!(next.inductor().unwrap().material(), &source);
            assert_eq!(&next.circuit().storage()[..2], unchanged_storage.as_slice());
            assert_eq!(next.circuit().graph().branches[3].2,
                Branch::Inductor { henries: next.inductor().unwrap().inductance_h() });
            held = held.step_with_contacts_candidate(vec![resistor(0, 2.0)], contacts.clone(), &[0.4],
                0.02, controls(), || false).unwrap().into_parts().0;
            moving = next;
        }
        assert!(moving.gap_m() < 1.0); assert!(moving.inductor().unwrap().inductance_h() > 1.0);
        assert!((moving.circuit().fluxes_wb()[1]-held.fluxes_wb()[1]).abs() > 1e-6,
            "gap motion must alter actual circuit evolution");
        assert!((moving.circuit().charges_c()[1]-held.charges_c()[1]).abs() > 1e-6);
    }
}

#[test]
fn millimetre_airgap_and_flux_linkage_use_physical_si_scales() {
    let (mu, area, gap, turns, initial_current, dt) = (1.25663706212e-6, 1e-3, 1e-3, 100, 0.1, 0.01);
    let inductance = mu*f64::from(turns).powi(2)*area/gap;
    let flux = inductance*initial_current;
    let model = MaterialFieldMotion::new(shorted(magnetic(0, mu, area, gap, turns), flux),
        0, body(1000.0), 0.0, [5e-4, 2e-3]).unwrap();
    let b = mu*f64::from(turns)*initial_current/gap;
    let force = -b*b*area/(2.0*mu);
    let step = model.step_candidate(vec![], FieldMotionInput { electrical_ports: &[0.0],
        external_force_n: 0.0, duration_s: dt },
        MaterialCircuitSolve { energy_tolerance_j: 1e-12, power_tolerance_w: 1e-12 }, || false).unwrap();
    near(step.model().gap_m(), gap+0.5*force/model.body().inertia().mass()*dt*dt, 1e-11);
    near(step.model().momentum_kg_m_s(), force*dt, 1e-11);
    near(step.model().circuit().fluxes_wb()[0], flux, 1e-14);
    near(step.step().field_work_on_body_j, step.step().kinetic_energy_j[1], 1e-12);
    near(step.step().energy_residual_j, 0.0, 1e-12);
}

#[test]
fn magnetic_domain_cancellation_and_missing_thermal_field_energy_refuse_atomically() {
    let model = MaterialFieldMotion::new(shorted(magnetic(0, 1.0, 0.5, 1.0, 1), 0.2),
        0, body(1000.0), 0.0, [0.99, 1.5]).unwrap();
    let before = model.clone();
    assert!(model.step_candidate(vec![], FieldMotionInput { electrical_ports: &[0.0],
        external_force_n: 0.0, duration_s: 1.0 }, controls(), || false).is_err());
    assert_eq!(model, before);
    let input = FieldMotionInput { electrical_ports: &[0.0], external_force_n: 0.0, duration_s: 0.01 };
    let mut polls = 0;
    model.step_candidate(vec![], input, controls(), || { polls += 1; false }).unwrap();
    let mut count = 0;
    assert!(model.step_candidate(vec![], input, controls(), || { count += 1; count == polls }).is_err());
    assert_eq!(model, before);
    let material = source(PropertyClaim {
        key: PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 1.0), (400.0, 2.0)], dims: MAGNETIC_PERMEABILITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic thermal permeability".into(),
            license: "test-only".into(), artifact: None } });
    let magnetic = ResolvedInductor::resolve(&material, &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly, 0.5, 1.0, 2).unwrap();
    let circuit = shorted(MaterialStorage { branch: 0, element: MaterialStorageElement::Inductor(magnetic) }, 0.2);
    assert!(MaterialFieldMotion::new(circuit, 0, body(1000.0), 0.0, [0.5, 1.5]).is_err());
}
