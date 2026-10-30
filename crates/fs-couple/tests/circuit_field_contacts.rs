//! Independent two-body caloric/field equations, not a measured material fit.
#[path = "electrical_contacts/source.rs"]
mod contact_source;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell,
    contacts::material::{ElectricalContactNetwork, MaterialElectricalContact, PreparedElectricalContactNetwork},
    field_storage::{self, FieldThermalBinding},
    material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork},
};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_couple::material_circuit::contacts::SeriesElectricalContact;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::{
    conductor::{ELECTRICAL_RESISTIVITY_DIMS, ElectricalPropertyBasis},
    field_storage::{ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        ResolvedCapacitor, ResolvedInductor, thermal::ThermalFieldStorage},
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS},
    thermal::{source::ThermalMaterialPoint, interface::ThermalInterfacePoint},
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const A: CellId = CellId(7);
const B: CellId = CellId(19);
const CP: [f64; 2] = [2.0, 3.0];
const BETA: f64 = 0.002;
const DT: f64 = 0.2;
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
}
fn card(i: usize, magnetic: bool, alpha: f64, g: f64) -> MaterialCard {
    let (name, dims) = if magnetic { (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS) }
        else { (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS) };
    let mut claims = ClaimSet::new();
    for (name, dims, low, high) in [
        (name, dims, 1.0 - 50.0 * alpha, 1.0 + 200.0 * alpha),
        ("density", Density::DIMS, 1.0, 1.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, CP[i], CP[i]),
        ("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, 1.0, 1.0),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, g, g),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, low), (500.0, high)], dims },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T", QuantitySpec::dimensional(Temperature::DIMS), 250.0, 500.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance {
                source: "manufactured two-body field/contact equations".into(), license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: format!("test {i}"), phase: "solid".into(),
        process: "manufactured".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(250.0), upper: &point(500.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-9).unwrap()
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2, maximum_circuit_trials: 500, maximum_newton_iterations: 24,
        difference_step_j: 1e-4, cell_energy_tolerance_j: 1e-7,
        thermal_energy_tolerance_j: 1e-7, total_energy_tolerance_j: 2e-7,
        circuit: MaterialCircuitSolve { energy_tolerance_j: 1e-9, power_tolerance_w: 1e-9 },
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-9, relative_tolerance: 1e-14,
            forcing_minimum: 1e-4, linear_restart: 4, ..NewtonKrylovConfig::default() },
    }
}
fn fixture<R>(magnetic: bool, alpha: f64, rc: Option<f64>, fraction: f64, g: f64,
    run: impl FnOnce(&PreparedElectricalContactNetwork<'_, '_>, &MaterialNetwork<'_>,
        DepositionSimulation<'_>, MaterialCircuit, FieldThermalBinding) -> R,
) -> R {
    let cards = [card(0, magnetic, alpha, g), card(1, magnetic, 0.0, g)];
    let curves = cards.each_ref().map(caloric);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell {
        cell: [A, B][i], material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context: &context,
            temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = [CellResistor { destination: ResistorCell { branch: 1, cell: A },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: 1.0 }];
    let links = [ConductionLink { link: 3, distance_a_m: 0.5, distance_b_m: 0.5 }];
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let mut events: Vec<_> = [0, 1].into_iter().map(|i| TimedEvent {
        time_s: 0.0, event: Event::Birth(DepositedCell {
            id: [A, B][i], segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64),
            mass_kg: 1.0, specific_h_j_kg: 50.0 * CP[i],
            geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 1.0, 0.0, 0.0],
                volume_m3: 1.0, surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }),
    }).collect();
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 3, a: A, b: B, area_m2: 1.0, conductance_w_k: 999.0,
    }) });
    let thermal = DepositionSimulation::new(events,
        &[(MaterialId(0), curves[0].curve() as &dyn ThermalLaw), (MaterialId(1), curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 100 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-8 }, || Ok(())).unwrap().0;
    let interface = contact_source::card(&cards[0], &cards[1], rc.unwrap_or(0.0), BETA);
    let contact_point = contact_source::point();
    let contacts: Vec<_> = rc.into_iter().map(|_| MaterialElectricalContact {
        link: 3, branch: 1, selection: MaterialPropertySelection::SingleClaimOnly,
        heat_fraction_to_first: fraction,
        point: ThermalInterfacePoint { interface: &interface, materials: [&cards[0], &cards[1]],
            context: &contact_point, temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [Some("melt_a"), Some("melt_b")], mean_temperature_axis: None },
    }).collect();
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let network = ElectricalContactNetwork { bulk: &material, contacts: &contacts, electrical_bases: &bases };
    let prepared = network.prepare(&thermal, || false).unwrap();
    let initial = prepared.resolve(&thermal, || false).unwrap();
    let (branch, element, law) = if magnetic {
        let value = ResolvedInductor::resolve(&cards[0], &context, MaterialPropertySelection::SingleClaimOnly, 1.0, 1.0, 1).unwrap();
        let law = ThermalFieldStorage::inductor(&cards[0], &value, "T", [250.0, 500.0]).unwrap();
        (Branch::Inductor { henries: f64::NAN }, MaterialStorageElement::Inductor(value), law)
    } else {
        let value = ResolvedCapacitor::resolve(&cards[0], &context, MaterialPropertySelection::SingleClaimOnly, 1.0, 1.0).unwrap();
        let law = ThermalFieldStorage::capacitor(&cards[0], &value, "T", [250.0, 500.0]).unwrap();
        (Branch::Capacitor { farads: f64::NAN }, MaterialStorageElement::Capacitor(value), law)
    };
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 2,
        branches: vec![(1, 0, branch), (1, 0, Branch::Resistor { ohms: f64::NAN })], transformers: vec![] },
        initial.bulk.resistors, vec![MaterialStorage { branch: 0, element }],
        if magnetic { vec![10.0] } else { vec![] }, if magnetic { vec![] } else { vec![10.0] }, &[]).unwrap()
        .rebind_contacts_candidate(initial.contacts, &[]).unwrap();
    run(&prepared, &material, thermal, circuit, FieldThermalBinding {
        branch: 0, cell: A, law, maximum_abs_coordinate: 10.0,
    })
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 4e-7, "{a:.16e} != {b:.16e}"); }
fn coordinate(c: &MaterialCircuit, magnetic: bool) -> f64 {
    if magnetic { c.fluxes_wb()[0] } else { c.charges_c()[0] }
}
fn field_u(alpha: f64, t: f64, z: f64) -> f64 {
    let k = 1.0 + alpha * (t - 300.0);
    0.5 * z * z / k * (1.0 + t * alpha / k)
}
// Independent equations eliminate B's backward-Euler temperature, then solve
// A's actual finite heat balance with original RC/RL midpoint storage decay.
fn reference(magnetic: bool, alpha: f64, rc: f64, fraction: f64, g: f64,
    before: [f64; 3], dt: f64) -> [f64; 6] {
    let [ta0, tb0, z0] = before;
    let k0 = 1.0 + alpha * (ta0 - 300.0);
    let evaluate = |ta: f64| {
        let k = 1.0 + alpha * (ta - 300.0);
        let contact = rc * (1.0 + BETA * (ta - 300.0));
        let r = 1.0 + contact;
        let ratio = if magnetic { dt * r / (2.0 * k) } else { dt / (2.0 * k * r) };
        let z = z0 * (1.0-ratio) / (1.0+ratio);
        let current = 0.5 * (z0+z) / k / if magnetic { 1.0 } else { r };
        let bulk_heat = dt * current * current;
        let contact_heat = bulk_heat * contact;
        let psi_change = 0.5*z0*z0*(1.0/k - 1.0/k0);
        let ts_change = 0.5*alpha*(ta*z*z/(k*k) - ta0*z0*z0/(k0*k0));
        let reversible = -psi_change-ts_change;
        let tb = (CP[1]*tb0 + (1.0-fraction)*contact_heat + dt*g*ta)/(CP[1]+dt*g);
        let residual = CP[0]*(ta-ta0)-bulk_heat-fraction*contact_heat-reversible-dt*g*(tb-ta);
        (residual, [ta, tb, z, bulk_heat, contact_heat, reversible])
    };
    let (mut lo, mut hi) = (250.0, 450.0);
    assert!(evaluate(lo).0 < 0.0 && evaluate(hi).0 > 0.0);
    for _ in 0..90 { let t = 0.5*(lo+hi); if evaluate(t).0 > 0.0 { hi=t; } else { lo=t; } }
    evaluate(0.5*(lo+hi)).1
}

#[test]
fn both_field_kinds_close_two_body_contact_heat_and_reversible_energy() {
    for magnetic in [false, true] {
        for alpha in [-0.001, 0.0, 0.001] {
            for fraction in [0.0, 0.3, 1.0] {
                fixture(magnetic, alpha, Some(2.0), fraction, 0.02, |prepared, _, thermal, circuit, field| {
                    let expected = reference(magnetic, alpha, 2.0, fraction, 0.02, [300.0, 300.0, 10.0], DT);
                    let next = prepared.field_storage_candidate(&thermal, DT, &circuit, &[], &[field], controls(), &[0.0; 2], || false).unwrap();
                    let report = next.report();
                    near(next.thermal().temperature_k(A).unwrap(), expected[0]);
                    near(next.thermal().temperature_k(B).unwrap(), expected[1]);
                    near(coordinate(next.circuit(), magnetic), expected[2]);
                    near(report.coupled.electrical.resistors[0].heat_j, expected[3]);
                    near(report.coupled.electrical.contacts[0].heat_j, expected[4]);
                    near(report.coupled.source.reversible_heat_j, expected[5]);
                    near(report.coupled.deposited_joule_j, expected[3]+expected[4]);
                    near(report.coupled.source.network.links[0].conductance_w_k, 0.02);
                    near(report.external_heat_j, 0.0); near(report.total_energy_residual_j, 0.0);
                    near(next.thermal().energy_j() + field_u(alpha, expected[0], expected[2]),
                        thermal.energy_j() + field_u(alpha, 300.0, 10.0));
                    assert_eq!(thermal.time_s(), 0.0); assert_eq!(coordinate(&circuit, magnetic), 10.0);
                    if alpha != 0.0 { assert!(expected[5].abs() > 0.01); }
                });
            }
        }
    }
}

#[test]
fn material_only_field_slope_changes_contact_heating_and_retained_electrical_history() {
    let mut ends = Vec::new();
    for alpha in [-0.001, 0.001] {
        fixture(false, alpha, Some(2.0), 0.3, 0.02, |prepared, _, mut thermal, mut circuit, field| {
            let fields = [field];
            let total = thermal.energy_j() + field_u(alpha, 300.0, 10.0);
            let mut expected = [300.0, 300.0, 10.0];
            for n in 1..=10 {
                let target = n as f64 * 0.1;
                let r = reference(false, alpha, 2.0, 0.3, 0.02, expected, target-thermal.time_s());
                let candidate = prepared.field_storage_candidate(&thermal, target, &circuit, &[], &fields, controls(), &[0.0; 2], || false).unwrap();
                let (t, c, report) = candidate.into_parts(); thermal=t; circuit=c;
                expected=[r[0], r[1], r[2]];
                near(thermal.temperature_k(A).unwrap(), r[0]); near(thermal.temperature_k(B).unwrap(), r[1]);
                near(coordinate(&circuit, false), r[2]);
                near(thermal.energy_j()+report.coupled.source.field_internal_energy_j[1], total);
                near(circuit.contacts()[0].source().resolved().query().axes()["Ta"], r[0]);
            }
            ends.push(expected);
        });
    }
    assert!((ends[0][0]-ends[1][0]).abs() > 0.01);
    assert!((ends[0][1]-ends[1][1]).abs() > 0.001);
    assert!((ends[0][2]-ends[1][2]).abs() > 0.001);
}

#[test]
fn refusal_cancellation_and_retry_preserve_all_three_energy_owners() {
    fixture(false, 0.001, Some(2.0), 0.3, 0.02, |prepared, _, thermal, circuit, field| {
        let saved = thermal.cells().clone(); let original = circuit.clone();
        let mut fields = [field];
        fields[0].maximum_abs_coordinate = 9.0;
        assert!(prepared.field_storage_candidate(&thermal, DT, &circuit, &[], &fields, controls(), &[0.0; 2], || false).is_err());
        fields[0].maximum_abs_coordinate = 10.0;
        let mut small = controls(); small.maximum_circuit_trials = 1;
        assert!(matches!(prepared.field_storage_candidate(&thermal, DT, &circuit, &[], &fields, small, &[0.0; 2], || false), Err(CircuitEnthalpyError::CircuitBudget)));
        let old_contact = &circuit.contacts()[0];
        let wrong = circuit.rebind_contacts_candidate(vec![SeriesElectricalContact::new(
            old_contact.id(), old_contact.branch(), old_contact.source().clone(), 2.0).unwrap()], &[]).unwrap();
        assert!(prepared.field_storage_candidate(&thermal, DT, &wrong, &[], &fields, controls(), &[0.0; 2], || false).is_err());
        let mut polls=0;
        let expected=prepared.field_storage_candidate(&thermal, DT, &circuit, &[], &fields, controls(), &[0.0; 2], || { polls+=1; false }).unwrap();
        let mut current=0;
        assert!(matches!(prepared.field_storage_candidate(&thermal, DT, &circuit, &[], &fields, controls(), &[0.0; 2], || { current+=1; current==polls }), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.cells(), &saved); assert_eq!(circuit, original); assert_eq!(thermal.time_s(), 0.0);
        let replay=prepared.field_storage_candidate(&thermal, DT, &circuit, &[], &fields, controls(), &[0.0; 2], || false).unwrap();
        assert_eq!(replay.circuit(), expected.circuit()); assert_eq!(replay.thermal().cells(), expected.thermal().cells());
        assert_eq!(replay.report().coupled.source, expected.report().coupled.source);
    });
}

#[test]
fn contact_free_composition_matches_original_field_storage_owner() {
    for magnetic in [false, true] {
        fixture(magnetic, 0.001, None, 0.3, 0.02, |prepared, material, thermal, circuit, field| {
            let fields = [field];
            let a = prepared.field_storage_candidate(&thermal, DT, &circuit, &[], &fields, controls(), &[0.0; 2], || false).unwrap();
            let b = field_storage::step_candidate(material, &thermal, DT, &circuit, &[], &fields,
                &[(1, ElectricalPropertyBasis::Resistivity)], controls(), &[0.0; 2], || false).unwrap();
            assert_eq!(a.circuit(), b.circuit()); assert_eq!(a.thermal().cells(), b.thermal().cells());
            assert_eq!(a.report().coupled.source, b.report().coupled.source);
        });
    }
}
