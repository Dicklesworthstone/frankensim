//! Manufactured coupled SI equations, not an empirical material validation.
use super::*;
use crate::frame_vibration::{FrameVibrationBasis, FrameVibrationBudget};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use crate::circuit_enthalpy::{ResistorCell, material::{MaterialCell, CellResistor, ConductionLink}};
use crate::circuit_enthalpy::contacts::material::MaterialElectricalContact;
use crate::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_material::{
    field_storage::{ResolvedCapacitor, ResolvedInductor, thermal::ThermalFieldStorage,
        ELECTRIC_PERMITTIVITY_DIMS, MAGNETIC_PERMEABILITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_PROPERTY},
    conductor::ELECTRICAL_RESISTIVITY_DIMS,
    electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY},
    phase::{ResolvedHeatCapacityCurve, HeatCapacitySourceQuery, SolidLiquidPhase},
    section::{ConcentricElasticSection, ConcentricMaterialLayer},
    state_point::{MaterialPropertySelection, ScalarPropertyRequirement, ScalarAdmissibility,
        SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        resolve_material_state_point},
    thermal::{source::ThermalMaterialPoint, interface::ThermalInterfacePoint},
};
use fs_matdb::*;
use fs_qty::{Density, Pressure, DynViscosity, Temperature, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_evidence::ValidityDomain;
use fs_solid::frame2d::{Frame2d, Geometry, FrameMassFormulation, FrameViscousSection};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_thermal::enthalpy::*;
use fs_solver::NewtonKrylovConfig;

const RATE: u32 = 256;
const DT: f64 = 1.0 / RATE as f64;
const A: CellId = CellId(7);
const B: CellId = CellId(19);
const RADIUS: f64 = 0.1;
const K0: f64 = 0.1;
const Z0: f64 = 0.5;
const FRACTION: f64 = 0.25;
fn area() -> f64 { core::f64::consts::PI * RADIUS * RADIUS }
fn capacity(i: usize) -> f64 { 3.0 * area() * (i + 1) as f64 }
fn mass() -> f64 { area() } // rho*A*L/3, original consistent axial inertia
fn tq() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", tq(), t).unwrap().with("omega", 2048.0).unwrap()
}
#[derive(Clone, Copy)]
struct Spec { magnetic: bool, slope: f64, contact: bool, cached: bool, capacity_failure: bool }
fn spec() -> Spec { Spec { magnetic: false, slope: 0.001, contact: false, cached: true, capacity_failure: false } }
fn expansion(s: Spec) -> f64 { if s.capacity_failure { 0.0001 } else { 0.00002 } }
fn card(i: usize, s: Spec) -> MaterialCard {
    let (field, dims) = if s.magnetic { (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS) }
        else { (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS) };
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        ("density", Density::DIMS, 3.0, 0.0),
        ("young_modulus", Pressure::DIMS, 1000.0, -1.0),
        ("eta", DynViscosity::DIMS, 2.0, 0.0),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS, expansion(s), 0.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, (i + 1) as f64, 0.0),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, if s.contact { 0.02 } else { 0.0 }, 0.0),
        ("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, 5.0 * area(), 0.01 * area()),
        (field, dims, K0 / area(), K0 / area() * s.slope),
    ] {
        let mut validity = ValidityDomain::unconstrained().with_quantity("T", tq(), 200.0, 600.0);
        if name == "eta" { validity = validity.with("omega", 0.0, 4096.0); }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if slope == 0.0 { PropertyValue::Scalar { value, dims } }
                else { PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value - 100.0 * slope), (600.0, value + 300.0 * slope)], dims } },
            validity, uncertainty: UncertaintyModel::Unstated,
            interpolation: if slope == 0.0 { InterpolationPolicy::ConstantWithinValidity } else { InterpolationPolicy::LinearInside },
            observations: vec![], provenance: Provenance { source: "manufactured shared field/frame body".into(), license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: format!("body {i}"), phase: "solid".into(), process: "manufactured".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn controls() -> FrameThermalSolve {
    FrameThermalSolve { maximum_cells: 2, maximum_material_trials: 256, maximum_newton_iterations: 24,
        difference_step_j: 1e-7, cell_energy_tolerance_j: 1e-10, thermal_energy_tolerance_j: 1e-10,
        total_energy_tolerance_j: 5e-10, newton: NewtonKrylovConfig { absolute_tolerance: 1e-13,
            relative_tolerance: 1e-14, forcing_minimum: 1e-4, linear_restart: 2, ..NewtonKrylovConfig::default() } }
}
fn input() -> FrameElectricalInput<'static> {
    FrameElectricalInput { ports: &[], solve: MaterialCircuitSolve { energy_tolerance_j: 1e-10, power_tolerance_w: 1e-10 } }
}
fn fixture<R>(s: Spec, run: impl FnOnce(Result<FrameThermoelasticSystem<'_, '_>>, DepositionSimulation<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let cards = [card(0, s), card(1, s)];
    let context = point(300.0);
    let curves = cards.each_ref().map(|card| ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(200.0), upper: &point(600.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-9).unwrap());
    let mut events: Vec<_> = (0..2).map(|i| TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: [A, B][i], segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64),
        mass_kg: 3.0 * area(), specific_h_j_kg: 100.0 * (i + 1) as f64,
        geometry: CellGeometry { start_m: [0.0, i as f64, 0.0], end_m: [1.0, i as f64, 0.0],
            volume_m3: area(), surface_area_m2: 2.0 * area() + 2.0 * core::f64::consts::PI * RADIUS },
        boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
    }) }).collect();
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink { id: 3, a: A, b: B,
        area_m2: area(), conductance_w_k: 999.0 }) });
    let thermal = DepositionSimulation::new(events,
        &[(MaterialId(0), curves[0].curve() as &dyn ThermalLaw), (MaterialId(1), curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 128 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: DT, energy_tolerance_j: 1e-10 }, || Ok(())).unwrap().0;
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell { cell: [A, B][i], material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context: &context, temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)) });
    let links = [ConductionLink { link: 3, distance_a_m: 0.1, distance_b_m: 0.1 }];
    let thermal_source = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
    let resistors = [CellResistor { destination: ResistorCell { branch: 1, cell: A },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: area() }];
    let bulk = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let requirements = [("density", Density::DIMS), ("young_modulus", Pressure::DIMS), ("eta", DynViscosity::DIMS),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS)]
        .map(|(name, dims)| ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    let section = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer { outer_radius_m: RADIUS,
        material: resolve_material_state_point(&cards[0], &context, &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap() }]).unwrap();
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0); frame.add_node(1.0, 0.0); frame.fix(0, [true; 3]); frame.fix(1, [false, true, true]);
    frame.add_material_viscous_element(0, 1, 5, FrameViscousSection::try_new(section, "eta", "eta").unwrap()).unwrap();
    let basis = FrameVibrationBasis::from_frame(&frame, 1, FrameVibrationBudget { max_frame_dofs: 6, max_dense_entries: 36 }, || false).unwrap();
    let mut component = basis.coupled_component(RATE, &[0.0], ModalAcousticTimeBudget::audible_reference(), 18, || false).unwrap();
    component.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: 0.01 / basis.modes()[0].displacement[3], velocity_m_sqrt_kg_per_s: 0.0 }]).unwrap();
    let owner = component.into_thermoelastic_system(&thermal, &thermal_source, &[FrameLayerCell { element: 0, layer: 0, cell: A }],
        0.05, ModalCouplingConfig { max_modes: 1, max_connections: 6, max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0, maximum_abs_pressure_pa: 100.0,
            maximum_abs_connection_force_n: 1000.0, solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-10 },
        s.cached.then_some(PreparedModalBudget { max_modes: 1, max_dense_bytes: 8 * 1024 * 1024, max_multiply_terms: 1_000_000_000 }), &gate).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let (branch, element, law) = if s.magnetic {
        let value = ResolvedInductor::resolve(&cards[0], &context, selection, area(), 1.0, 1).unwrap();
        let law = ThermalFieldStorage::inductor(&cards[0], &value, "T", [200.0, 600.0]).unwrap();
        (Branch::Inductor { henries: f64::NAN }, MaterialStorageElement::Inductor(value), law)
    } else {
        let value = ResolvedCapacitor::resolve(&cards[0], &context, selection, area(), 1.0).unwrap();
        let law = ThermalFieldStorage::capacitor(&cards[0], &value, "T", [200.0, 600.0]).unwrap();
        (Branch::Capacitor { farads: f64::NAN }, MaterialStorageElement::Capacitor(value), law)
    };
    let maximum_abs_coordinate = if s.capacity_failure {
        let bound = (0.9999 * capacity(0) / -law.minimum_heat_capacity_j_k(1.0).unwrap()).sqrt();
        assert!(capacity(0) + law.minimum_heat_capacity_j_k(bound).unwrap() > 0.0);
        assert!(owner.members[0].capacity_floor_j_k[0] > 0.0);
        assert!(owner.members[0].capacity_floor_j_k[0] + law.minimum_heat_capacity_j_k(bound).unwrap() < 0.0);
        bound
    } else { Z0 };
    let fields = vec![FieldThermalBinding { branch: 0, cell: A, law, maximum_abs_coordinate }];
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, ELECTRICAL_CONTACT_RESISTANCE_DIMS),
        value: PropertyValue::Curve { abscissa: "Ta".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 1.6 * area()), (600.0, 3.2 * area())], dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS },
        validity: ValidityDomain::unconstrained().with("Ta", 200.0, 600.0).with("Tb", 200.0, 600.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside, observations: vec![],
        provenance: Provenance { source: "manufactured interface".into(), license: "test-only".into(), artifact: None } }).unwrap();
    let interface = InterfaceSystemCard::assemble(SurfaceSpec { material: cards[0].id().clone(), texture_frame: "A".into() },
        SurfaceSpec { material: cards[1].id().clone(), texture_frame: "B".into() },
        SystemContext { medium: "test".into(), third_body: None, environment: "test".into(), history: "fixed".into() }, claims, vec![]).unwrap();
    let interface_context = QueryPoint::new().with("Ta", 300.0).unwrap().with("Tb", 300.0).unwrap();
    let contacts = [MaterialElectricalContact { link: 3, branch: 1, selection: MaterialPropertySelection::SingleClaimOnly,
        heat_fraction_to_first: FRACTION, point: ThermalInterfacePoint { interface: &interface,
            materials: [&cards[0], &cards[1]], context: &interface_context,
            temperature_axes: ["Ta", "Tb"], liquid_fraction_axes: [None, None], mean_temperature_axis: None } }];
    let network = ElectricalContactNetwork { bulk: &bulk, contacts: if s.contact { &contacts } else { &[] }, electrical_bases: &bases };
    let prepared = network.prepare(&thermal, || false).unwrap();
    let initial = prepared.resolve(&thermal, || false).unwrap();
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 2,
        branches: vec![(1, 0, branch), (1, 0, Branch::Resistor { ohms: f64::NAN })], transformers: vec![] },
        initial.bulk.resistors, vec![MaterialStorage { branch: 0, element }],
        if s.magnetic { vec![Z0] } else { vec![] }, if s.magnetic { vec![] } else { vec![Z0] }, &[]).unwrap()
        .rebind_contacts_candidate(initial.contacts, &[]).unwrap();
    let owner = if s.contact {
        owner.with_material_contact_field_circuit(&thermal, &network, circuit, fields, &gate)
    } else { owner.with_material_field_circuit(&thermal, &bulk, circuit, &bases, fields, &gate) };
    run(owner, thermal)
}

fn near(a: f64, b: f64, tol: f64) { assert!((a - b).abs() <= tol, "{a:.17e} != {b:.17e}"); }
fn state(owner: &FrameThermoelasticSystem<'_, '_>) -> [f64; 2] {
    let scale = owner.basis().modes()[0].displacement[3];
    let q = owner.mechanical().components()[0].states()[0];
    [q.displacement_m_sqrt_kg * scale, q.velocity_m_sqrt_kg_per_s * scale]
}
fn coordinate(c: &MaterialCircuit, s: Spec) -> f64 {
    if s.magnetic { c.fluxes_wb()[0] } else { c.charges_c()[0] }
}
fn elastic(s: Spec, t: f64, x: f64) -> (f64, f64) {
    let e = 1000.0 - (t - 300.0);
    let a = expansion(s);
    let strain = x - a * (t - 300.0);
    let psi = 0.5 * area() * e * strain * strain;
    let entropy = area() * (e * a * strain + 0.5 * strain * strain);
    (psi, psi + t * entropy)
}
fn field_u(s: Spec, t: f64, z: f64) -> f64 {
    let f = 1.0 + s.slope * (t - 300.0);
    0.5 * z * z / (K0 * f) * (1.0 + t * s.slope / f)
}
// Exact frozen-T damped oscillator. Work/loss comes from independent Simpson
// integration of c*v^2, not any native port or combined-energy report.
fn motion(s: Spec, t: f64, initial: [f64; 2]) -> ([f64; 2], f64) {
    let k = area() * (1000.0 - (t - 300.0));
    let c = 2.0 * area();
    let decay = c / (2.0 * mass());
    let omega = (k / mass() - decay * decay).sqrt();
    let rest = expansion(s) * (t - 300.0);
    let x0 = initial[0] - rest;
    let b = (initial[1] + decay * x0) / omega;
    let at = |time: f64| {
        let (sin, cos) = (omega * time).sin_cos();
        let excess = (-decay * time).exp() * (x0 * cos + b * sin);
        [rest + excess, (-decay * time).exp() * omega * (-x0 * sin + b * cos) - decay * excess]
    };
    let mut loss = 0.0;
    for i in 0..=128 {
        let w = if i == 0 || i == 128 { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
        loss += w * c * at(DT * i as f64 / 128.0)[1].powi(2);
    }
    (at(DT), loss * DT / 384.0)
}
// Two physical caloric balances plus actual RC/RL midpoint evolution. Eliminate
// passive B analytically; solve A using independent total internal energy.
fn reference(s: Spec, old_t: [f64; 2], old_z: f64, old_x: [f64; 2]) -> ([f64; 2], f64, [f64; 2], [f64; 3]) {
    let (x, loss) = motion(s, old_t[0], old_x);
    let evaluate = |ta: f64| {
        let k = K0 * (1.0 + s.slope * (ta - 300.0));
        let bulk = 5.0 * (1.0 + 0.002 * (ta - 300.0));
        let contact = if s.contact { 2.0 * (1.0 + 0.002 * (ta - 300.0)) } else { 0.0 };
        let r = bulk + contact;
        let ratio = if s.magnetic { DT * r / (2.0 * k) } else { DT / (2.0 * k * r) };
        let z = old_z * (1.0 - ratio) / (1.0 + ratio);
        let current = 0.5 * (old_z + z) / k / if s.magnetic { 1.0 } else { r };
        let heat = [DT * current * current * bulk, DT * current * current * contact];
        let g = if s.contact { area() * 0.02 / 0.2 } else { 0.0 };
        let tb = (capacity(1) * old_t[1] + (1.0 - FRACTION) * heat[1] + DT * g * ta) / (capacity(1) + DT * g);
        let residual = capacity(0) * (ta - old_t[0]) + capacity(1) * (tb - old_t[1])
            + 0.5 * mass() * (x[1] * x[1] - old_x[1] * old_x[1])
            + elastic(s, ta, x[0]).1 - elastic(s, old_t[0], old_x[0]).1
            + field_u(s, ta, z) - field_u(s, old_t[0], old_z);
        (residual, [ta, tb], z, [loss, heat[0], heat[1]])
    };
    let (mut lo, mut hi) = (200.0, 550.0);
    assert!(evaluate(lo).0 < 0.0 && evaluate(hi).0 > 0.0);
    for _ in 0..90 { let mid = 0.5 * (lo + hi); if evaluate(mid).0 > 0.0 { hi = mid; } else { lo = mid; } }
    let (_, t, z, heat) = evaluate(0.5 * (lo + hi));
    (t, z, x, heat)
}

#[test]
fn field_contact_and_frame_materials_share_actual_internal_energy_and_motion() {
    for magnetic in [false, true] {
        for contact in [false, true] {
            for cached in [false, true] {
                let s = Spec { magnetic, contact, cached, ..spec() };
                fixture(s, |owner, mut thermal| {
                    let mut owner = owner.unwrap();
                    let gate = CancelGate::new_clock_free();
                    let total = thermal.energy_j() + owner.mechanical_internal_energy_j().unwrap() + field_u(s, 300.0, Z0);
                    for _ in 0..12 {
                        let t0 = [thermal.temperature_k(A).unwrap(), thermal.temperature_k(B).unwrap()];
                        let expected = reference(s, t0, coordinate(owner.electrical().unwrap(), s), state(&owner));
                        let report = owner.step_electrical(&mut thermal, &[0.0], input(), controls(), &[0.0; 2], &gate).unwrap();
                        near(thermal.temperature_k(A).unwrap(), expected.0[0], 2e-7);
                        near(thermal.temperature_k(B).unwrap(), expected.0[1], 2e-7);
                        near(coordinate(owner.electrical().unwrap(), s), expected.1, 2e-9);
                        near(state(&owner)[0], expected.2[0], 2e-10); near(state(&owner)[1], expected.2[1], 2e-9);
                        near(report.deposited_material_heat_j, expected.3[0], 2e-10);
                        near(report.deposited_joule_heat_j, expected.3[1] + expected.3[2], 2e-9);
                        let field = report.field_storage.as_ref().unwrap();
                        near(field.internal_energy_j[1], field_u(s, expected.0[0], expected.1), 2e-9);
                        near(thermal.energy_j() + owner.mechanical_internal_energy_j().unwrap() + field.internal_energy_j[1], total, 4e-9);
                        near(report.total_energy_residual_j, 0.0, 5e-10);
                        near(report.external_thermal_heat_j, 0.0, 2e-10);
                        near(report.entropy.accounting_residual_j_k, 0.0, 2e-11);
                        assert!(!report.entropy.has_unclassified_sources);
                        assert!(field.reversible_heat_j != 0.0);
                        if contact { near(report.electrical.as_ref().unwrap().contacts[0].heat_j, expected.3[2], 2e-9); }
                    }
                });
            }
        }
    }
}

#[test]
fn changing_only_the_field_material_law_changes_temperature_and_subsequent_motion() {
    let mut results = Vec::new();
    for slope in [-0.001, 0.0, 0.001] {
        let s = Spec { slope, contact: true, ..spec() };
        fixture(s, |owner, mut thermal| {
            let mut owner = owner.unwrap(); let gate = CancelGate::new_clock_free();
            for _ in 0..32 {
                let report = owner.step_electrical(&mut thermal, &[0.0], input(), controls(), &[0.0; 2], &gate).unwrap();
                if slope == 0.0 {
                    assert_eq!(report.field_storage.unwrap().reversible_heat_j, 0.0);
                    assert_eq!(report.entropy.field_storage_change_j_k, 0.0);
                }
            }
            results.push((thermal.temperature_k(A).unwrap(), state(&owner)[0], coordinate(owner.electrical().unwrap(), s)));
        });
    }
    assert!((results[0].0 - results[2].0).abs() > 0.01);
    assert!((results[0].1 - results[2].1).abs() > 1e-8);
    assert!((results[0].2 - results[2].2).abs() > 1e-5);
}

#[test]
fn shared_cell_capacity_cannot_be_spent_once_by_each_physics_subsystem() {
    fixture(Spec { capacity_failure: true, ..spec() }, |owner, thermal| {
        // Fixture independently checked that both pairwise capacity gates pass,
        // while the actual three-contribution envelope fails.
        assert!(matches!(owner, Err(CircuitEnthalpyError::Invalid(
            "combined caloric, elastic and all co-located field heat capacities must stay positive"))));
        assert_eq!(thermal.time_s(), 0.0);
        near(thermal.temperature_k(A).unwrap(), 300.0, 1e-12);
    });
}

#[test]
fn a_late_refusal_preserves_charge_motion_field_sources_and_both_thermal_histories() {
    let run = |refuse: bool| fixture(Spec { contact: true, ..spec() }, |owner, mut thermal| {
        let mut owner = owner.unwrap(); let gate = CancelGate::new_clock_free();
        for _ in 0..3 { owner.step_electrical(&mut thermal, &[0.0], input(), controls(), &[0.0; 2], &gate).unwrap(); }
        if refuse {
            let circuit = owner.electrical().unwrap().clone(); let q = state(&owner);
            let cells = thermal.cells().clone(); let time = thermal.time_s();
            let sample = owner.mechanical().samples_rendered();
            let old_frame = owner.mechanical().last_frame().unwrap().clone();
            let mut control = controls(); control.maximum_material_trials = 1;
            assert!(owner.step_electrical(&mut thermal, &[0.0], input(), control, &[0.0; 2], &gate).is_err());
            let cancelled = CancelGate::new_clock_free(); cancelled.request();
            assert!(matches!(owner.step_electrical(&mut thermal, &[0.0], input(), controls(), &[0.0; 2], &cancelled), Err(CircuitEnthalpyError::Cancelled)));
            assert_eq!(owner.electrical().unwrap(), &circuit); assert_eq!(state(&owner), q);
            assert_eq!(thermal.cells(), &cells); assert_eq!(thermal.time_s(), time);
            assert_eq!(owner.mechanical().samples_rendered(), sample);
            assert_eq!(owner.mechanical().last_frame(), Some(&old_frame));
        }
        let report = owner.step_electrical(&mut thermal, &[0.0], input(), controls(), &[0.0; 2], &gate).unwrap();
        (state(&owner), owner.electrical().unwrap().clone(), thermal.cells().clone(), report.field_storage, report.entropy)
    });
    assert_eq!(run(true), run(false));
}
