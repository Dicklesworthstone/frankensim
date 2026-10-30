//! Manufactured two-cell reference: induction, magnetic entropy, contact and photons.
use super::*;
use fs_conduction::{EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity};
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_couple::circuit_enthalpy::{ResistorCell, contacts::{heat::*, material::ElectricalContactNetwork},
    magnetic::MagneticThermalError, material::{CellResistor, ConductionLink, MaterialNetwork}};
use fs_exec::{CancelGate, Cx};
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::{ClaimSelection, InterpolationPolicy, PropertyClaim, PropertyValue, Provenance,
    SelectionPolicy, UncertaintyModel};

#[path = "../../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiation_geometry;
#[path = "../../../../fs-conduction/tests/contact_material_state/common.rs"]
mod contact_source;

const DT: f64 = 0.1;
const IDS: [CellId; 2] = [CellId(7), CellId(9)];

fn with_heat<R>(epsilon: f64, resistance: f64, graph: bool,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &MaterialHeatNetwork<'_>, DepositionSimulation<'_>, &MaterialCard) -> R,
) -> R {
    let cards = [epsilon, 0.6].map(|emissivity| {
        let source = field::absolute_affine_magnetic_card(1.0 / 18.0, 0.0, 1.0, 0.001,
            ElectricalPropertyBasis::Resistivity);
        let mut claims = source.claims().clone();
        for (name, dims, value) in [("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
            (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, emissivity)] {
            claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims }, validity: fs_evidence::ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![], provenance: Provenance { source: "manufactured induction/contact/radiation".into(),
                    license: "CC0-1.0".into(), artifact: None } }).unwrap();
        }
        MaterialCard::assemble(source.id().clone(), claims, vec![]).unwrap()
    });
    let curves = thermal_curves(&cards);
    let context = field::absolute_point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell { cell: IDS[i], material_id: MaterialId(i as u64),
        curve: curves[i].curve(), point: ThermalMaterialPoint { material: &cards[i], context: &context,
            temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)) });
    let mut events = [0, 1].map(|i| TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: IDS[i], segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64),
        mass_kg: 1.0, specific_h_j_kg: [400.0, 560.0][i], geometry: CellGeometry {
            start_m: [0.0; 3], end_m: [1.0, 0.0, 0.0], volume_m3: 1.0, surface_area_m2: 2.0 },
        boundary: common::boundary(0.0) }) }).to_vec();
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 7, a: IDS[0], b: IDS[1], area_m2: 1.0, conductance_w_k: 999.0 }) });
    let thermal = DepositionSimulation::new(events,
        &[(MaterialId(0), curves[0].curve() as &dyn ThermalLaw), (MaterialId(1), curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 1000 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: DT, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0;
    let links = [ConductionLink { link: 7, distance_a_m: 0.5, distance_b_m: 0.5 }];
    let resistors = [CellResistor { destination: ResistorCell { branch: 0, cell: IDS[0] },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: 1.0 }];
    let bulk = MaterialNetwork { cells: &cells, resistors: if graph { &resistors } else { &[] }, links: &links };
    let bases = [(0, ElectricalPropertyBasis::Resistivity)];
    let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[],
        electrical_bases: if graph { &bases } else { &[] } };
    let mut contact = contact_source::claim("fa", fs_qty::Dims::NONE, resistance);
    for axis in ["aT", "bT"] { contact.validity = contact.validity.with_quantity(axis,
        fs_qty::QuantitySpec::dimensional(fs_qty::Temperature::DIMS), 290.0, 360.0); }
    let mut claims = fs_matdb::ClaimSet::new();
    claims.insert_claim(contact).unwrap();
    let interface = contact_source::interface(&cards, claims);
    let point = contact_source::point();
    let contacts = [MaterialContactLink { link: 7, sides: [0, 1].map(|i| ContactCell {
        cell: IDS[i], material_id: cells[i].material_id, curve: curves[i].curve() }),
        point: contact_source::binding(&cards, &interface, &point),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) }];
    let surfaces = [0, 1].map(|i| RadiatingCell { surface: i, cell: IDS[i], material_id: cells[i].material_id,
        curve: curves[i].curve(), material: cells[i].point, selection: MaterialPropertySelection::SingleClaimOnly });
    let enclosure = radiation_geometry::enclosure([0, 1].map(|i| SurfaceEmissivity::from_thermal_state(
        "manufactured plate", cells[i].point, curves[i].curve().state_at_specific_enthalpy([400.0, 560.0][i]).unwrap(),
        MaterialPropertySelection::SingleClaimOnly).unwrap().0));
    let radiation = MaterialRadiation { enclosure: &enclosure, cells: &surfaces,
        power_tolerance_w: 1e-8, radiosity_tolerance_w_m2: 1e-8 };
    let assembly = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &contacts, radiation: Some(&radiation) };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, fs_exec::StreamKey { seed: 59, kernel_id: 73, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        run(&cx, &gate, &assembly, thermal, &cards[0])
    })
}

fn bracket(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..64 { let mid = (lo + hi) / 2.0;
        if residual(mid) > 0.0 { hi = mid; } else { lo = mid; } }
    (lo + hi) / 2.0
}
fn radiation_into_first(ta: f64, tb: f64, epsilon: f64) -> f64 {
    5.670_374_419e-8 * (tb.powi(4) - ta.powi(4)) / (1.0 / epsilon + 1.0 / 0.6 - 1.0)
}
fn reference(t0: [f64; 2], a0: f64, epsilon: f64, resistance: f64) -> ([f64; 2], f64) {
    // Independent single-P1-mode field equation, two finite capacities and the
    // analytic two-gray-plate exchange. No production residual is reused.
    let mu = |t| 1.0 + 0.001 * (t - 300.0);
    let amplitude = |ta| a0 * (10.0 - 2.0 / mu(t0[0])) / (10.0 + 2.0 / mu(ta));
    let internal = |t, a: f64| 2.0 * a * a / mu(t) * (1.0 + t * 0.001 / mu(t));
    let sink = |ta| bracket(290.0, 360.0, |tb|
        4.0 * (tb - t0[1]) + DT * ((tb - ta) / (1.0 + resistance) + radiation_into_first(ta, tb, epsilon)));
    let ta = bracket(290.0, 360.0, |ta| 4.0 * (ta - t0[0]) + 4.0 * (sink(ta) - t0[1])
        + internal(ta, amplitude(ta)) - internal(t0[0], a0));
    ([ta, sink(ta)], amplitude(ta))
}

#[test]
fn material_contact_and_enclosure_change_the_retained_field_at_one_caloric_endpoint() {
    let mut outputs = Vec::new();
    for (epsilon, resistance) in [(0.2, 0.2), (0.8, 0.2), (0.2, 2.0)] {
        outputs.push(with_heat(epsilon, resistance, false, |cx, _, assembly, mut thermal, card| {
            let problem = problem(card);
            let mut field = initial_field(&problem, 1.0);
            let bindings = laws(card, &problem);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            for step in 1..=2 {
                let t0 = IDS.map(|id| thermal.temperature_k(id).unwrap());
                let a0 = field.potential_wb_m()[4];
                let expected = reference(t0, a0, epsilon, resistance);
                let next = prepared.step_eddy_thermal_candidate(cx, &thermal, step as f64 * DT, &field,
                    &[0.0; 4], &single_routes(), 1e-12, &bindings, common::config(), field_controls(),
                    field::limits(), &[0.0; 2]).unwrap();
                let report = next.report();
                let source = &report.coupled.source.material.network;
                assert_eq!(source.thermal_contacts.len(), 1);
                common::near(source.thermal_contacts[0].conductance_w_k, 1.0 / (1.0 + resistance));
                for i in 0..2 {
                    common::near(next.thermal().temperature_k(IDS[i]).unwrap(), expected.0[i]);
                    common::near(source.thermal_contacts[0].phases[i].temperature_k(), expected.0[i]);
                }
                common::near(next.field().potential_wb_m()[4], expected.1);
                common::near(report.coupled.deposited_joule_j, 10.0 * (expected.1 - a0).powi(2));
                let radiation = source.radiation.as_ref().unwrap();
                common::near(radiation.sources_w()[&IDS[0]], radiation_into_first(expected.0[0], expected.0[1], epsilon));
                for surface in radiation.surfaces() {
                    let cell = source.material.cells.iter().find(|cell| cell.cell == surface.cell).unwrap();
                    assert_eq!(surface.phase, cell.phase);
                }
                common::near(report.external_heat_j, 0.0);
                common::near(report.total_energy_residual_j, 0.0);
                common::near(report.coupled.non_electrical_heat_j - report.coupled.source.reversible_heat_j
                    - DT * radiation.net_heat_w(), report.external_heat_j);
                common::near(next.thermal().energy_j() - thermal.energy_j()
                    + report.coupled.source.field_internal_energy_j[1] - report.coupled.source.field_internal_energy_j[0], 0.0);
                let (next_thermal, next_field, _) = next.into_parts();
                thermal = next_thermal; field = next_field;
            }
            [thermal.temperature_k(IDS[0]).unwrap(), thermal.temperature_k(IDS[1]).unwrap(), field.potential_wb_m()[4]]
        }));
    }
    for changed in &outputs[1..] { for i in 0..3 { assert!((outputs[0][i] - changed[i]).abs() > 1e-6); } }
}

#[test]
fn full_heat_induction_refuses_unowned_electrical_graphs_and_cancel_without_advancing() {
    for graph in [false, true] {
        with_heat(0.2, 0.2, graph, |cx, gate, assembly, thermal, card| {
            let problem = problem(card);
            let field = initial_field(&problem, 1.0);
            let bindings = laws(card, &problem);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let before = (thermal.cells().clone(), field.potential_wb_m().to_vec());
            let run = || prepared.step_eddy_thermal_candidate(cx, &thermal, DT, &field,
                &[0.0; 4], &single_routes(), 1e-12, &bindings, common::config(), field_controls(),
                field::limits(), &[0.0; 2]);
            if graph {
                assert!(matches!(run(), Err(CircuitEnthalpyError::Source(MagneticThermalError::Heat(HeatNetworkError::Invalid(_))))));
            } else { assert!(run().is_ok()); }
            gate.request();
            assert!(matches!(run(), Err(CircuitEnthalpyError::Cancelled)));
            assert_eq!(thermal.cells(), &before.0);
            assert_eq!(field.potential_wb_m(), before.1);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(field.time_s(), 0.0);
        });
    }
}
