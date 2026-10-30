//! One material card drives structural motion, resistive heat and finite storage.
//! Manufactured SI fixtures; not experimental validation of a measured alloy.
#[path = "frame_electrical_thermal/contacts.rs"]
mod contacts;
#[path = "frame_electrical_thermal/heat.rs"]
mod heat;

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, ResistorCell,
    material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork},
};
use fs_couple::frame_vibration::{
    FrameVibrationBasis, FrameVibrationBudget, FrameVibrationComponent,
    thermal::{FrameLayerCell, FrameThermalSolve},
    thermoelastic::{FrameElectricalInput, FrameThermoelasticSystem},
};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig, PreparedModalBudget};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::conductor::{ElectricalPropertyBasis, ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_DIMS};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, FrameViscousSection, Geometry};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const RATE: u32 = 256;
const DT: f64 = 1.0 / RATE as f64;
const RADIUS: f64 = 0.1;
const RHO: f64 = 3.0;
const E0: f64 = 1000.0;
const ETA: f64 = 2.0;
const CAPACITANCE: f64 = 0.1;
const CHARGE: f64 = 0.5;
const CELL: CellId = CellId(7);
fn area() -> f64 { core::f64::consts::PI * RADIUS * RADIUS }
fn capacity() -> f64 { RHO * area() } // unit specific heat, length 1 m
fn mass() -> f64 { RHO * area() / 3.0 } // original consistent axial inertia
#[derive(Clone, Copy)]
struct Spec { resistance: f64, beta: f64, alpha: f64, e_slope: f64, conductivity: bool }
fn spec() -> Spec { Spec { resistance: 5.0, beta: 0.002, alpha: 0.0, e_slope: 0.0, conductivity: false } }
fn context(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
        .with("omega", 2048.0).unwrap()
}
fn card(s: Spec) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let electrical = if s.conductivity {
        assert_eq!(s.beta, 0.0); // equivalence case is constant, not reciprocal interpolation
        ("electrical_conductivity", 1.0 / (s.resistance * area()), ELECTRICAL_CONDUCTIVITY_DIMS, 0.0)
    } else { ("electrical_resistivity", s.resistance * area(), ELECTRICAL_RESISTIVITY_DIMS, s.resistance * area() * s.beta) };
    for (name, value, dims, slope) in [
        ("density", RHO, Density::DIMS, 0.0), ("young_modulus", E0, Pressure::DIMS, s.e_slope),
        ("eta", ETA, DynViscosity::DIMS, 0.0),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, s.alpha, INVERSE_TEMPERATURE_DIMS, 0.0),
        ("specific_heat_capacity", 1.0, SPECIFIC_HEAT_CAPACITY_DIMS, 0.0),
        ("thermal_conductivity", 0.0, fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS, 0.0), electrical,
    ] {
        let curved = slope != 0.0;
        let mut validity = ValidityDomain::unconstrained().with_quantity(
            "T", QuantitySpec::dimensional(Temperature::DIMS), 200.0, 600.0);
        if name == "eta" { validity = validity.with("omega", 0.0, 4096.0); }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if curved { PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, value - 100.0*slope), (600.0, value + 300.0*slope)], dims,
            }} else { PropertyValue::Scalar { value, dims } },
            validity, uncertainty: UncertaintyModel::Unstated,
            interpolation: if curved { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance { source: "manufactured frame/circuit thermal fixture".into(), license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "manufactured conductor".into(),
        phase: "solid".into(), process: "explicit linear properties".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn controls() -> FrameThermalSolve {
    FrameThermalSolve {
        maximum_cells: 1, maximum_material_trials: 128, maximum_newton_iterations: 16,
        difference_step_j: 1e-7, cell_energy_tolerance_j: 1e-10,
        thermal_energy_tolerance_j: 1e-10, total_energy_tolerance_j: 3e-10,
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-13, relative_tolerance: 1e-14,
            forcing_minimum: 1e-4, linear_restart: 1, ..NewtonKrylovConfig::default() },
    }
}
fn input(ports: &[f64]) -> FrameElectricalInput<'_> {
    FrameElectricalInput { ports, solve: MaterialCircuitSolve { energy_tolerance_j: 1e-10, power_tolerance_w: 1e-10 } }
}
fn frame_component(card: &MaterialCard, context: &QueryPoint, displacement: f64) -> FrameVibrationComponent {
    let requirements = [("density", Density::DIMS), ("young_modulus", Pressure::DIMS), ("eta", DynViscosity::DIMS),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS)]
        .map(|(name, dims)| ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    let section = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: RADIUS, material: resolve_material_state_point(card, context, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap(),
    }]).unwrap();
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0); frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]); frame.fix(1, [false, true, true]);
    frame.add_material_viscous_element(0, 1, 5, FrameViscousSection::try_new(section, "eta", "eta").unwrap()).unwrap();
    let basis = FrameVibrationBasis::from_frame(&frame, 1,
        FrameVibrationBudget { max_frame_dofs: 6, max_dense_entries: 36 }, || false).unwrap();
    let mut component = basis.coupled_component(RATE, &[0.0], ModalAcousticTimeBudget::audible_reference(), 18, || false).unwrap();
    component.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: displacement / basis.modes()[0].displacement[3], velocity_m_sqrt_kg_per_s: 0.0 }]).unwrap();
    component
}
fn fixture<R>(s: Spec, capacitor: bool, cached: bool,
    run: impl FnOnce(FrameThermoelasticSystem<'_, '_>, DepositionSimulation<'_>) -> R,
) -> R {
    let card = card(s);
    let context = context(300.0);
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &self::context(200.0), upper: &self::context(600.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-7).unwrap();
    let component = frame_component(&card, &context, if capacitor { 0.01 } else { 0.0 });
    let laws = [(MaterialId(0), curve.curve() as &dyn ThermalLaw)];
    let thermal = DepositionSimulation::new(vec![TimedEvent { time_s: 0.0,
        event: Event::Birth(DepositedCell { id: CELL, segment_id: 0, body_id: 0, material_id: MaterialId(0),
            mass_kg: RHO * area(), specific_h_j_kg: 100.0,
            geometry: CellGeometry { start_m: [0.0; 3], end_m: [1.0, 0.0, 0.0], volume_m3: area(), surface_area_m2: 2.0*area()+2.0*core::f64::consts::PI*RADIUS },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }),
    }], &laws, DepositionBudget { max_cells: 1, max_links: 0, max_events: 1, max_substeps: 1024 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: DT, energy_tolerance_j: 1e-10 }, || Ok(())).unwrap().0;
    let cells = [MaterialCell { cell: CELL, material_id: MaterialId(0), curve: curve.curve(),
        point: ThermalMaterialPoint { material: &card, context: &context, temperature_axis: "T", liquid_fraction_axis: None }, conductivity: None }];
    let links: [ConductionLink; 0] = [];
    let thermal_source = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
    let resistors = [CellResistor { destination: ResistorCell { branch: 0, cell: CELL }, selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0, area_m2: area() }];
    let circuit_source = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let bases = [(0, if s.conductivity { ElectricalPropertyBasis::Conductivity } else { ElectricalPropertyBasis::Resistivity })];
    let sources = circuit_source.resolve_with_bases(&thermal, &bases, || false).unwrap().resistors;
    let circuit = MaterialCircuit::new(CircuitGraph { node_count: 2,
        branches: if capacitor { vec![(1, 0, Branch::Resistor { ohms: 999.0 }), (1, 0, Branch::Capacitor { farads: CAPACITANCE })] }
            else { vec![(1, 0, Branch::Resistor { ohms: 999.0 }), (0, 1, Branch::CurrentSource { port: 0 })] }, transformers: vec![],
    }, sources, vec![], if capacitor { vec![CHARGE] } else { vec![] }, if capacitor { &[] } else { &[1.0] }).unwrap();
    let gate = CancelGate::new_clock_free();
    let owner = component.into_thermoelastic_system(&thermal, &thermal_source,
        &[FrameLayerCell { element: 0, layer: 0, cell: CELL }], 0.05,
        ModalCouplingConfig { max_modes: 1, max_connections: 6, max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0, maximum_abs_pressure_pa: 100.0,
            maximum_abs_connection_force_n: 1000.0, solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-10 },
        cached.then_some(PreparedModalBudget { max_modes: 1, max_dense_bytes: 8*1024*1024, max_multiply_terms: 1_000_000_000 }), &gate).unwrap();
    let owner = owner.with_material_circuit(&thermal, &circuit_source, circuit, &bases, &gate).unwrap();
    run(owner, thermal)
}
fn near(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a:e} != {b:e}, tolerance {tolerance:e}"); }
fn state(owner: &FrameThermoelasticSystem<'_, '_>) -> [f64; 2] {
    let shape = owner.basis().modes()[0].displacement[3];
    let state = owner.mechanical().components()[0].states()[0];
    [shape*state.displacement_m_sqrt_kg, shape*state.velocity_m_sqrt_kg_per_s]
}
fn potential(s: Spec, t: f64, x: f64) -> (f64, f64) {
    let e = E0+s.e_slope*(t-300.0);
    let strain = x-s.alpha*(t-300.0);
    let psi = 0.5*area()*e*strain*strain;
    let entropy = area()*(e*s.alpha*strain-0.5*s.e_slope*strain*strain);
    (psi, psi+t*entropy)
}
// Exact frozen-temperature mechanical trajectory and independent Simpson
// quadrature of eta*A*v^2, not the implementation's port work or heat report.
fn motion(s: Spec, t: f64, x: [f64; 2]) -> ([f64; 2], f64) {
    let k = area()*(E0+s.e_slope*(t-300.0));
    let c = ETA*area(); let a = c/(2.0*mass());
    let w = (k/mass()-a*a).sqrt();
    let rest = s.alpha*(t-300.0);
    let trajectory = |time: f64| {
        let (sin, cos) = (w*time).sin_cos();
        let b = (x[1]+a*(x[0]-rest))/w;
        let excess = (-a*time).exp()*((x[0]-rest)*cos+b*sin);
        [rest+excess, (-a*time).exp()*(-(x[0]-rest)*w*sin+b*w*cos)-a*excess]
    };
    let n = 128;
    let mut integral = 0.0;
    for i in 0..=n {
        let v = trajectory(DT*i as f64/n as f64)[1];
        let weight = if i==0 || i==n { 1.0 } else if i%2==0 { 2.0 } else { 4.0 };
        integral += weight*c*v*v;
    }
    (trajectory(DT), integral*DT/(3.0*n as f64))
}

#[test]
fn current_driven_source_resistance_and_frame_heat_solve_at_one_endpoint() {
    for cached in [false, true] {
        for current in [-1.0, 1.0] {
            let s = spec();
            fixture(s, false, cached, |mut owner, mut thermal| {
                let gate = CancelGate::new_clock_free();
                let mut theta = 0.0; let mut work = 0.0;
                let initial = thermal.energy_j();
                for sample in 1..=32 {
                    theta = (capacity()*theta + DT*current*current*s.resistance)
                        /(capacity()-DT*current*current*s.resistance*s.beta);
                    let report = owner.step_electrical(&mut thermal, &[0.0], input(&[current]), controls(), &[0.0], &gate).unwrap();
                    let electrical = report.electrical.as_ref().unwrap();
                    let r = s.resistance*(1.0+s.beta*theta);
                    near(thermal.temperature_k(CELL).unwrap(), 300.0+theta, 1e-7);
                    near(electrical.resistors[0].current_a.abs(), current.abs(), 1e-10);
                    near(electrical.resistors[0].material.conductor.resistance_ohm(), r, 1e-9);
                    near(report.deposited_joule_heat_j, DT*current*current*r, 1e-10);
                    near(report.deposited_material_heat_j, 0.0, 1e-15);
                    near(report.entropy.material_loss_generation_j_k, report.deposited_joule_heat_j/(300.0+theta), 1e-13);
                    assert!(!report.entropy.has_unclassified_sources);
                    work += electrical.native.supplied;
                    near(thermal.energy_j()-initial-work, 0.0, 3e-9);
                    assert_eq!(owner.mechanical().samples_rendered(), sample);
                    near(thermal.time_s(), sample as f64*DT, 0.0);
                    near(state(&owner)[0], 0.0, 1e-15);
                }
            });
        }
    }
}

#[test]
fn finite_charge_material_expansion_and_native_motion_close_without_repeatedly_draining_storage() {
    let s = Spec { alpha: 0.0001, e_slope: -1.0, ..spec() };
    for cached in [false, true] {
        fixture(s, true, cached, |mut owner, mut thermal| {
            let gate = CancelGate::new_clock_free();
            let initial = thermal.energy_j()+owner.mechanical_internal_energy_j().unwrap()+CHARGE*CHARGE/(2.0*CAPACITANCE);
            let mut t = 300.0; let mut x = [0.01, 0.0]; let mut charge = CHARGE;
            for _ in 0..32 {
                let (next_x, loss) = motion(s, t, x);
                let stress_work = potential(s, t, next_x[0]).0-potential(s, t, x[0]).0;
                let trial = |next_t: f64| {
                    let r = s.resistance*(1.0+s.beta*(next_t-300.0));
                    let a = DT/(2.0*CAPACITANCE*r);
                    let q = charge*(1.0-a)/(1.0+a);
                    let heat = (charge-q)*(charge+q)/(2.0*CAPACITANCE);
                    let residual = capacity()*(next_t-t)+potential(s, next_t, next_x[0]).1
                        -potential(s, t, x[0]).1-stress_work-loss-heat;
                    (residual, q, heat)
                };
                let (mut lo, mut hi) = (250.0, 450.0);
                assert!(trial(lo).0 < 0.0 && trial(hi).0 > 0.0);
                for _ in 0..70 { let mid = 0.5*(lo+hi); if trial(mid).0 > 0.0 { hi=mid; } else { lo=mid; } }
                let next_t = 0.5*(lo+hi); let (_, next_q, joule) = trial(next_t);
                let report = owner.step_electrical(&mut thermal, &[0.0], input(&[]), controls(), &[0.0], &gate).unwrap();
                near(thermal.temperature_k(CELL).unwrap(), next_t, 2e-7);
                near(owner.electrical().unwrap().charges_c()[0], next_q, 2e-9);
                near(state(&owner)[0], next_x[0], 2e-10); near(state(&owner)[1], next_x[1], 2e-9);
                near(report.deposited_material_heat_j, loss, 2e-11);
                near(report.deposited_joule_heat_j, joule, 2e-10);
                let electrical = report.electrical.as_ref().unwrap();
                near(electrical.native.supplied, 0.0, 1e-12);
                near(thermal.energy_j()+owner.mechanical_internal_energy_j().unwrap()+electrical.stored_after_j, initial, 5e-9);
                (t, x, charge) = (next_t, next_x, next_q);
            }
        });
    }
}

#[test]
fn data_only_conductor_substitution_changes_heating_expansion_and_following_motion() {
    let mut endpoints = Vec::new();
    for conductivity in [false, true] {
        for resistance in [3.0, 6.0] {
            let s = Spec { resistance, beta: 0.0, alpha: 0.0001, e_slope: -1.0, conductivity };
            fixture(s, false, true, |mut owner, mut thermal| {
                let gate = CancelGate::new_clock_free();
                for _ in 0..32 {
                    owner.step_electrical(&mut thermal, &[0.0], input(&[1.0]), controls(), &[0.0], &gate).unwrap();
                }
                endpoints.push((thermal.temperature_k(CELL).unwrap(), state(&owner)));
            });
        }
    }
    assert!(endpoints[1].0 > endpoints[0].0);
    assert!((endpoints[1].1[0]-endpoints[0].1[0]).abs() > 1e-6);
    for i in 0..2 {
        near(endpoints[i].0, endpoints[i+2].0, 1e-8);
        near(endpoints[i].1[0], endpoints[i+2].1[0], 1e-10);
        near(endpoints[i].1[1], endpoints[i+2].1[1], 1e-9);
    }
}

#[test]
fn refusal_keeps_motion_charge_heat_and_clock_together_and_retry_matches_uninterrupted_history() {
    let s = Spec { alpha: 0.0001, e_slope: -1.0, ..spec() };
    let run = |fail: bool| fixture(s, true, true, |mut owner, mut thermal| {
        let gate = CancelGate::new_clock_free();
        for _ in 0..3 { owner.step_electrical(&mut thermal, &[0.0], input(&[]), controls(), &[0.0], &gate).unwrap(); }
        if fail {
            let motion = state(&owner); let circuit = owner.electrical().unwrap().clone();
            let cells = thermal.cells().clone(); let time = thermal.time_s();
            let frame = owner.mechanical().last_frame().unwrap().clone();
            // Ordinary stepping must not freeze an attached electrical history.
            assert!(owner.step(&mut thermal, &[0.0], controls(), &[0.0], &gate).is_err());
            let mut exhausted = controls(); exhausted.maximum_material_trials=1;
            assert!(owner.step_electrical(&mut thermal, &[0.0], input(&[]), exhausted, &[0.0], &gate).is_err());
            assert!(owner.step_electrical(&mut thermal, &[0.0], input(&[f64::NAN]), controls(), &[0.0], &gate).is_err());
            let cancel = CancelGate::new_clock_free(); cancel.request();
            assert!(matches!(owner.step_electrical(&mut thermal, &[0.0], input(&[]), controls(), &[0.0], &cancel), Err(CircuitEnthalpyError::Cancelled)));
            assert_eq!(state(&owner), motion); assert_eq!(owner.electrical(), Some(&circuit));
            assert_eq!(thermal.cells(), &cells); assert_eq!(thermal.time_s(), time);
            assert_eq!(owner.mechanical().last_frame(), Some(&frame));
            assert_eq!(owner.mechanical().samples_rendered(), 3);
        }
        owner.step_electrical(&mut thermal, &[0.0], input(&[]), controls(), &[0.0], &gate).unwrap();
        (state(&owner), thermal.cells().clone(), owner.electrical().unwrap().clone(), thermal.time_s())
    });
    assert_eq!(run(false), run(true));
}
