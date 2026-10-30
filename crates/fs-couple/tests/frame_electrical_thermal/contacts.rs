//! Existing interface source laws feed both finite bodies and frame motion.
use super::*;
use fs_couple::circuit_enthalpy::contacts::material::{
    ElectricalContactNetwork, MaterialElectricalContact,
};
use fs_couple::frame_vibration::thermal::FrameThermalError;
use fs_couple::material_circuit::contacts::SeriesElectricalContact;
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};

const OTHER: CellId = CellId(17);
const LINK: u64 = 9;
#[derive(Clone, Copy)]
struct Interface {
    resistance: f64,
    slope: f64,
    fraction: f64,
    upper_k: f64,
    wrong_initial_area: bool,
}
fn interface() -> Interface {
    Interface { resistance: 2.0, slope: 0.002, fraction: 0.25,
        upper_k: 600.0, wrong_initial_area: false }
}
fn contact_card(cards: &[MaterialCard; 2], c: Interface) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    let value = c.resistance * area();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, ELECTRICAL_CONTACT_RESISTANCE_DIMS),
        value: PropertyValue::Curve {
            abscissa: "Ta".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, value*(1.0-100.0*c.slope)),
                (c.upper_k, value*(1.0+(c.upper_k-300.0)*c.slope))],
            dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("Ta", QuantitySpec::dimensional(Temperature::DIMS), 200.0, c.upper_k)
            .with_quantity("Tb", QuantitySpec::dimensional(Temperature::DIMS), 200.0, 600.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance {
            source: "manufactured frame electrical interface".into(), license: "test-only".into(), artifact: None,
        },
    }).unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec { material: cards[0].id().clone(), texture_frame: "first".into() },
        SurfaceSpec { material: cards[1].id().clone(), texture_frame: "second".into() },
        SystemContext { medium: "explicit contact".into(), third_body: None,
            environment: "manufactured reference".into(), history: "fixed geometry".into() },
        claims, vec![],
    ).unwrap()
}

fn fixture<R>(s: Spec, c: Interface, capacitor: bool, cached: bool,
    run: impl FnOnce(core::result::Result<FrameThermoelasticSystem<'_, '_>, FrameThermalError>, DepositionSimulation<'_>) -> R,
) -> R {
    let cards = [card(s), card(s)];
    let context = context(300.0);
    let curves = cards.each_ref().map(|card| ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &super::context(200.0), upper: &super::context(600.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-7).unwrap());
    let component = frame_component(&cards[0], &context, if capacitor { 0.01 } else { 0.0 });
    let mut events: Vec<_> = [CELL, OTHER].into_iter().enumerate().map(|(i, id)| TimedEvent {
        time_s: 0.0, event: Event::Birth(DepositedCell {
            id, segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64),
            mass_kg: RHO * area(), specific_h_j_kg: 100.0,
            geometry: CellGeometry { start_m: [0.0, i as f64, 0.0], end_m: [1.0, i as f64, 0.0],
                volume_m3: area(), surface_area_m2: 2.0*area()+2.0*core::f64::consts::PI*RADIUS },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }),
    }).collect();
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: LINK, a: CELL, b: OTHER, area_m2: area(), conductance_w_k: 999.0,
    }) });
    let laws = [(MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1), curves[1].curve() as &dyn ThermalLaw)];
    let thermal = DepositionSimulation::new(events, &laws,
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 1024 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: DT, energy_tolerance_j: 1e-10 }, || Ok(())).unwrap().0;
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0,1].map(|i| MaterialCell {
        cell: [CELL,OTHER][i], material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context: &context, temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let links = [ConductionLink { link: LINK, distance_a_m: 0.5, distance_b_m: 0.5 }];
    let thermal_source = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
    let resistors = [CellResistor { destination: ResistorCell { branch: 0, cell: CELL },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: area() }];
    let bulk = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let bases = [(0, ElectricalPropertyBasis::Resistivity)];
    let card = contact_card(&cards, c);
    let point = QueryPoint::new()
        .with_quantity("Ta", QuantitySpec::dimensional(Temperature::DIMS), 300.0).unwrap()
        .with_quantity("Tb", QuantitySpec::dimensional(Temperature::DIMS), 300.0).unwrap();
    let contacts = [MaterialElectricalContact { link: LINK, branch: 0,
        point: ThermalInterfacePoint { interface: &card, materials: [&cards[0],&cards[1]],
            context: &point, temperature_axes: ["Ta","Tb"], liquid_fraction_axes: [None,None], mean_temperature_axis: None },
        selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: c.fraction }];
    let source = ElectricalContactNetwork { bulk: &bulk, contacts: &contacts, electrical_bases: &bases };
    let mut initial = source.resolve(&thermal, || false).unwrap();
    if c.wrong_initial_area {
        let old = &initial.contacts[0];
        initial.contacts[0] = SeriesElectricalContact::new(old.id(), old.branch(), old.source().clone(), 2.0*old.area_m2()).unwrap();
    }
    let circuit = MaterialCircuit::new_with_contacts(CircuitGraph { node_count: 2,
        branches: if capacitor { vec![(1,0,Branch::Resistor { ohms: 999.0 }), (1,0,Branch::Capacitor { farads: CAPACITANCE })] }
            else { vec![(1,0,Branch::Resistor { ohms: 999.0 }), (0,1,Branch::CurrentSource { port: 0 })] }, transformers: vec![],
    }, initial.bulk.resistors, initial.contacts, vec![], if capacitor { vec![CHARGE] } else { vec![] },
        if capacitor { &[] } else { &[1.0] }).unwrap();
    let gate = CancelGate::new_clock_free();
    let owner = component.into_thermoelastic_system(&thermal, &thermal_source,
        &[FrameLayerCell { element: 0, layer: 0, cell: CELL }], 0.05,
        ModalCouplingConfig { max_modes: 1, max_connections: 6, max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0, maximum_abs_pressure_pa: 100.0,
            maximum_abs_connection_force_n: 1000.0, solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-10 },
        cached.then_some(PreparedModalBudget { max_modes: 1, max_dense_bytes: 8*1024*1024, max_multiply_terms: 1_000_000_000 }), &gate).unwrap();
    run(owner.with_material_contact_circuit(&thermal, &source, circuit, &gate), thermal)
}
fn controls() -> FrameThermalSolve {
    let mut controls = super::controls();
    controls.maximum_cells = 2;
    controls.newton.linear_restart = 2;
    controls
}

#[test]
fn two_sided_contact_heat_matches_independent_endpoint_equations() {
    let s = Spec { beta: 0.0, ..spec() };
    for cached in [false, true] {
        for fraction in [0.0, 0.25, 1.0] {
            for current in [-1.0, 1.0] {
                let c = Interface { fraction, ..interface() };
                fixture(s, c, false, cached, |owner, mut thermal| {
                    let mut owner = owner.unwrap();
                    let gate = CancelGate::new_clock_free();
                    let routes = owner.electrical_contact_heat_destinations();
                    assert_eq!(routes[0].cells, [CELL,OTHER]);
                    near(routes[0].heat_fraction_to_first, fraction, 0.0);
                    let mut excess = [0.0;2];
                    let initial = thermal.energy_j();
                    let mut work = 0.0;
                    for _ in 0..16 {
                        excess[0] = (capacity()*excess[0]+DT*current*current*(s.resistance+fraction*c.resistance))
                            /(capacity()-DT*current*current*fraction*c.resistance*c.slope);
                        let rc = c.resistance*(1.0+c.slope*excess[0]);
                        let contact_heat = DT*current*current*rc;
                        excess[1] += (1.0-fraction)*contact_heat/capacity();
                        let report = owner.step_electrical(&mut thermal, &[0.0], input(&[current]), controls(), &[0.0;2], &gate).unwrap();
                        let electrical = report.electrical.as_ref().unwrap();
                        near(thermal.temperature_k(CELL).unwrap(), 300.0+excess[0], 2e-7);
                        near(thermal.temperature_k(OTHER).unwrap(), 300.0+excess[1], 2e-7);
                        near(electrical.contacts[0].contact.resistance_ohm(), rc, 2e-10);
                        near(electrical.contacts[0].heat_j, contact_heat, 2e-10);
                        near(report.deposited_joule_heat_j, DT*current*current*s.resistance+contact_heat, 2e-10);
                        near(report.source.links[0].conductance_w_k, 0.0, 0.0);
                        assert!(!report.entropy.has_unclassified_sources);
                        work += electrical.native.supplied;
                        near(thermal.energy_j()-initial-work, 0.0, 5e-9);
                    }
                });
            }
        }
    }
}

#[test]
fn contact_data_changes_finite_charge_both_temperatures_and_thermoelastic_motion() {
    let s = Spec { alpha: 0.0001, e_slope: -1.0, ..spec() };
    let mut ends = Vec::new();
    for resistance in [1.0, 4.0] {
        let c = Interface { resistance, ..interface() };
        fixture(s, c, true, true, |owner, mut thermal| {
            let mut owner = owner.unwrap();
            let gate = CancelGate::new_clock_free();
            let initial = thermal.energy_j()+owner.mechanical_internal_energy_j().unwrap()+CHARGE*CHARGE/(2.0*CAPACITANCE);
            let (mut ta, mut tb, mut q, mut x) = (300.0, 300.0, CHARGE, [0.01,0.0]);
            for _ in 0..32 {
                let (next_x, loss) = motion(s, ta, x);
                let work = potential(s, ta, next_x[0]).0-potential(s, ta, x[0]).0;
                let trial = |t: f64| {
                    let rb = s.resistance*(1.0+s.beta*(t-300.0));
                    let rc = c.resistance*(1.0+c.slope*(t-300.0));
                    let a = DT/(2.0*CAPACITANCE*(rb+rc));
                    let nq = q*(1.0-a)/(1.0+a);
                    let i = (q-nq)/DT;
                    let hb = DT*i*i*rb;
                    let hc = DT*i*i*rc;
                    (capacity()*(t-ta)+potential(s,t,next_x[0]).1-potential(s,ta,x[0]).1
                        -work-loss-hb-c.fraction*hc, nq, hb, hc)
                };
                let (mut lo, mut hi) = (250.0,450.0);
                assert!(trial(lo).0<0.0 && trial(hi).0>0.0);
                for _ in 0..70 { let mid = 0.5*(lo+hi); if trial(mid).0>0.0 {hi=mid;} else {lo=mid;} }
                let next_t = 0.5*(lo+hi);
                let (_, next_q, hb, hc) = trial(next_t);
                tb += (1.0-c.fraction)*hc/capacity();
                let report = owner.step_electrical(&mut thermal, &[0.0], input(&[]), controls(), &[0.0;2], &gate).unwrap();
                near(thermal.temperature_k(CELL).unwrap(), next_t, 3e-7);
                near(thermal.temperature_k(OTHER).unwrap(), tb, 3e-7);
                near(owner.electrical().unwrap().charges_c()[0], next_q, 2e-9);
                near(state(&owner)[0], next_x[0], 2e-10);
                near(state(&owner)[1], next_x[1], 2e-9);
                let electrical = report.electrical.as_ref().unwrap();
                near(electrical.resistors[0].heat_j, hb, 2e-10);
                near(electrical.contacts[0].heat_j, hc, 2e-10);
                near(report.deposited_joule_heat_j, hb+hc, 2e-10);
                near(thermal.energy_j()+owner.mechanical_internal_energy_j().unwrap()+electrical.stored_after_j, initial, 6e-9);
                (ta,q,x) = (next_t,next_q,next_x);
            }
            ends.push((ta,tb,q,x));
        });
    }
    assert!((ends[0].0-ends[1].0).abs()>0.01);
    assert!((ends[0].1-ends[1].1).abs()>0.01);
    assert!(ends[1].2>ends[0].2);
    assert!((ends[0].3[0]-ends[1].3[0]).abs()>1e-7);
}

#[test]
fn interface_source_refusal_and_exhaustion_preserve_all_histories_and_retry() {
    let run = |fail| fixture(Spec { beta: 0.0, ..spec() },
        Interface { upper_k: 300.001, ..interface() }, false, true, |owner, mut thermal| {
            let mut owner = owner.unwrap();
            let gate = CancelGate::new_clock_free();
            for _ in 0..2 {
                owner.step_electrical(&mut thermal, &[0.0], input(&[0.001]), controls(), &[0.0;2], &gate).unwrap();
            }
            if fail {
                let cells = thermal.cells().clone(); let time = thermal.time_s();
                let circuit = owner.electrical().unwrap().clone(); let x = state(&owner);
                let frame = owner.mechanical().last_frame().unwrap().clone();
                // The full-current endpoint exceeds the original interface
                // source domain. No alternate source or clamped value is used.
                assert!(owner.step_electrical(&mut thermal, &[0.0], input(&[1.0]), controls(), &[0.0;2], &gate).is_err());
                let mut exhausted = controls(); exhausted.maximum_material_trials = 1;
                assert!(owner.step_electrical(&mut thermal, &[0.0], input(&[0.001]), exhausted, &[0.0;2], &gate).is_err());
                let cancel = CancelGate::new_clock_free(); cancel.request();
                assert!(matches!(owner.step_electrical(&mut thermal, &[0.0], input(&[0.001]), controls(), &[0.0;2], &cancel), Err(CircuitEnthalpyError::Cancelled)));
                assert_eq!(thermal.cells(), &cells); assert_eq!(thermal.time_s(), time);
                assert_eq!(owner.electrical(), Some(&circuit)); assert_eq!(state(&owner), x);
                assert_eq!(owner.mechanical().last_frame(), Some(&frame));
                assert_eq!(owner.mechanical().samples_rendered(), 2);
            }
            owner.step_electrical(&mut thermal, &[0.0], input(&[0.001]), controls(), &[0.0;2], &gate).unwrap();
            (thermal.cells().clone(), state(&owner), owner.electrical().unwrap().clone())
        });
    assert_eq!(run(false), run(true));
    fixture(spec(), Interface { wrong_initial_area: true, ..interface() }, false, false,
        |owner, thermal| { assert!(owner.is_err()); assert_eq!(thermal.time_s(), 0.0); });
}
