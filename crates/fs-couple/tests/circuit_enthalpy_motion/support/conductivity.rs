//! Nonlinear source response through the original live support/enthalpy owners.
use super::*;
use fs_couple::circuit_enthalpy::contacts::ContactHeatCells;
use fs_couple::circuit_enthalpy::motion::support::step_conductivity_candidate;
use fs_couple::material_circuit::MaterialCircuitError;
use fs_couple::material_circuit::conductivity::{ConductionAxis, ConductionCoordinate,
    ConductivityIteration, ConductivitySolve, ConductorResponse};
use fs_couple::material_circuit::conductivity::contact::{ContactAxis, ContactCoordinate,
    ContactDirection, ContactResponse};
use fs_couple::material_circuit::contacts::SeriesElectricalContact;
use fs_qty::Dims;
use std::cell::Cell;

const FIELD: Dims = Dims([1, 1, -3, 0, -1, 0]);
const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);
const AXES: [ConductionAxis<'static>; 1] = [ConductionAxis {
    name: "E", coordinate: ConductionCoordinate::ElectricFieldMagnitude }];
const ROUTES: [ContactHeatCells; 1] = [ContactHeatCells { id: 7,
    cells: [DESTINATIONS[0].cell, DESTINATIONS[1].cell], heat_fraction_to_first: 0.25 }];

fn iteration() -> ConductivitySolve {
    ConductivitySolve { maximum_iterations: 96, absolute_resistance_tolerance_ohm: 1e-12,
        relative_resistance_tolerance: 1e-12 }
}

fn response_card(beta: f64, field_conductivity: Option<f64>, bounded: Option<&str>)
    -> (MaterialCard, ResolvedHeatCapacityCurve, Vec<ScalarPropertyRequirement>)
{
    let (original, _, requirements) = transport_card(1.0, 1e-3, 0.8);
    let mut claims = ClaimSet::new();
    for (_, old) in original.claims().claims_ordered() {
        let mut claim = old.clone();
        let slope = if claim.key.name() == ELECTRICAL_RESISTIVITY_PROPERTY { Some((0.004, beta)) }
            else if claim.key.name() == "thermal_conductivity" { field_conductivity.map(|k| (0.002, k)) }
            else { None };
        if let Some((thermal, field)) = slope {
            claim.value = PropertyValue::Surface2D { axes: ["T".into(), "E".into()],
                axis_dims: [Temperature::DIMS, FIELD],
                samples: [[280.0, 0.0], [400.0, 0.0], [400.0, 6.0], [280.0, 6.0]]
                    .map(|[t, e]| ([t, e], 1.0 + thermal * (t - 300.0) + field * e)).to_vec(),
                triangles: vec![[0, 1, 2], [0, 2, 3]], dims: claim.key.quantity().dims() };
            claim.validity = claim.validity.with_quantity("E", QuantitySpec::dimensional(FIELD), 0.0, 6.0);
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        if Some(claim.key.name()) == bounded {
            claim.validity = claim.validity.with_quantity("E", QuantitySpec::dimensional(FIELD), 0.0, 0.1);
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let query = |t| point(t).with_quantity("E", QuantitySpec::dimensional(FIELD), 0.0).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS), temperature_axis: "T",
        lower: &query(280.0), upper: &query(400.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, curve, requirements)
}

fn with_response_support<R>(beta: f64, field_conductivity: Option<f64>, bounded: Option<&str>,
    run: impl FnOnce(&MaterialNetwork<'_>, &PreparedMaterialNetwork<'_, '_>,
        &[CellConductionResponse<'_>], &[MaterialCard; 2], &SupportThermalBinding<'_>,
        DepositionSimulation<'_>, MaterialFieldMotion) -> R) -> R
{
    let (a, ca, ra) = response_card(beta, field_conductivity, bounded);
    let (b, cb, rb) = transport_card(3.0, 2e-3, 0.6);
    let cards = [a, b]; let curves = [ca, cb]; let requirements = [ra, rb];
    let context = [point(300.0).with_quantity("E", QuantitySpec::dimensional(FIELD), 0.0).unwrap(), point(300.0)];
    let layers = (0..2).map(|i| ConcentricMaterialLayer {
        outer_radius_m: ((i + 1) as f64 / core::f64::consts::PI).sqrt(),
        material: resolve_material_state_point(&cards[i], &context[i], &requirements[i],
            MaterialPropertySelection::SingleClaimOnly).unwrap(),
    }).collect();
    let law = ConcentricThermoelasticLaw::try_new(ConcentricSectionTemperatureLaw::try_new(
        ConcentricElasticSection::try_new(0.0, layers).unwrap(), &[&cards[0], &cards[1]], &["T", "T"],
    ).unwrap()).unwrap();
    let binding = SupportThermalBinding { law,
        layer_cells: DESTINATIONS.map(|d| d.cell).to_vec(), volume_tolerance_m3: 1e-12 };
    let thermal = thermal_with_contact(&curves, &binding.law, true);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell { cell: DESTINATIONS[i].cell,
        material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context: &context[i],
            temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)) });
    let resistors = [0, 1].map(|i| CellResistor {
        area_m2: binding.law.reference().layer_integrals()[i].area_m2(), length_m: 1.0,
        destination: DESTINATIONS[i], selection: MaterialPropertySelection::SingleClaimOnly });
    let links = links();
    let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let responses = [CellConductionResponse { branch: 1, basis: ElectricalPropertyBasis::Resistivity, axes: &AXES },
        CellConductionResponse { branch: 2, basis: ElectricalPropertyBasis::Resistivity, axes: &[] }];
    let prepared = network.prepare_conductivity(&thermal, &responses, || false).unwrap();
    let seeds = prepared.resolve_conductivity_seeds(&thermal, &responses, || false).unwrap();
    let support = AxialSupport::from_thermoelastic_law(&binding.law, &[300.0; 2], 1.0, 1.0, 0.2).unwrap();
    let motion = MaterialFieldMotion::new_with_support(driven_circuit(seeds.resistors), 0,
        body(100.0), support, 0.005, [0.9, 1.1]).unwrap();
    run(&network, &prepared, &responses, &cards, &binding, thermal, motion)
}

#[test]
fn prepared_nonlinear_support_uses_actual_field_for_resistance_and_heat_transport() {
    let run = |beta| with_response_support(beta, Some(0.3), None,
        |_, prepared, responses, _, binding, mut thermal, mut motion| {
        let initial = total_internal_layers(&thermal, &motion, [1e-3, 2e-3]);
        let mut supplied = 0.0;
        for n in 1..=4 {
            let target = f64::from(n) * 0.05; let dt = target - thermal.time_s();
            let next = prepared.field_motion_support_conductivity_candidate(&thermal, target,
                &motion, &[0.0, 3.0], 0.0, responses, binding, config(), iteration(), &[0.0; 2], || false).unwrap();
            let report = next.report();
            let electrical = &report.coupled.electrical;
            let temperatures = DESTINATIONS.map(|d| next.thermal().temperature_k(d.cell).unwrap());
            let sample = &electrical.evaluated[0];
            near(sample.point.electric_field_v_m, 3.0);
            near(sample.conductor.resistance_ohm(), 1.0 + 0.004 * (temperatures[0] - 300.0) + 3.0 * beta);
            let k = 1.0 + 0.002 * (temperatures[0] - 300.0) + 0.3 * sample.point.electric_field_v_m;
            let conductance = 1.0 / (0.5 / k + 0.5 / 4.0);
            near(report.coupled.source.network.links[0].conductance_w_k, conductance);
            let refreshed = report.coupled.source.network.cells.iter()
                .find(|s| s.cell == DESTINATIONS[0].cell).unwrap().conductivity.as_ref().unwrap();
            near(refreshed.query().axes()["E"], sample.point.electric_field_v_m);
            for (i, destination) in DESTINATIONS.iter().enumerate() {
                let old = &thermal.cells()[&destination.cell];
                let new = &next.thermal().cells()[&destination.cell];
                near(old.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg),
                    electrical.field.resistors[i].heat_j + report.coupled.source.layers[i].reversible_heat_j
                        + dt * conductance * (temperatures[1 - i] - temperatures[i]));
            }
            near(report.coupled.deposited_joule_j, electrical.field.heat_j);
            near(report.total_energy_residual_j, 0.0);
            assert!(electrical.iterations > 1);
            assert!(report.coupled.circuit_trials > electrical.iterations);
            supplied += electrical.field.electrical_supply_j;
            near(total_internal_layers(next.thermal(), next.motion(), [1e-3, 2e-3]), initial + supplied);
            let (t, m, _) = next.into_parts(); thermal = t; motion = m;
        }
        [thermal.temperature_k(DESTINATIONS[0].cell).unwrap(), motion.gap_m(), motion.velocity_m_s()]
    });
    let lower = run(0.02); let higher = run(0.08);
    assert!(lower[0] > higher[0]);
    assert!((lower[1] - higher[1]).abs() > 1e-8);
    assert!((lower[2] - higher[2]).abs() > 1e-7);
}

#[test]
fn nonlinear_contact_support_counts_native_trials_and_keeps_all_histories_atomic() {
    with_response_support(0.05, None, None, |network, prepared, responses, cards, binding, thermal, original| {
        let interface = source::response_card(&cards[0], &cards[1], "V", VOLTAGE, 0.5, 0.002, 0.1, 6.0);
        let circuit = original.circuit().rebind_contacts_candidate(vec![
            source::response_contact(7, 1, &interface, 1.0, "V", VOLTAGE)], &[0.0, 3.0]).unwrap();
        let motion = MaterialFieldMotion::new_with_support(circuit, 0, original.body().clone(),
            original.support().unwrap().clone(), original.velocity_m_s(), [0.9, 1.1]).unwrap();
        let source_responses = [ConductorResponse { branch: 1, card: &cards[0], axes: &AXES },
            ConductorResponse { branch: 2, card: &cards[1], axes: &[] }];
        let axes = [ContactAxis { name: "V", coordinate: ContactCoordinate::VoltageMagnitude }];
        let contact_responses = [ContactResponse { id: 7, card: &interface,
            direction: ContactDirection::AlongBranch, axes: &axes }];
        let native_trials = Cell::new(0usize); let resolutions = Cell::new(0usize);
        let refuse_at = Cell::new(None);
        let run = |controls, iteration: ConductivityIteration, cancelled: &mut dyn FnMut() -> bool| {
            step_conductivity_candidate(network, &thermal, 0.05, &motion, &[0.0, 3.0], 0.0,
                &source_responses, Some(&contact_responses), &ROUTES, binding, controls, iteration,
                &[0.0; 2], |state| {
                    let bulk = prepared.resolve_conductivity_seeds(state, responses, || false)
                        .map_err(|e| e.to_string())?;
                    let old = &motion.circuit().contacts()[0];
                    let point = old.source().resolved().query()
                        .with("Ta", state.temperature_k(DESTINATIONS[0].cell).map_err(|e| e.to_string())?)
                        .map_err(|e| e.to_string())?
                        .with("Tb", state.temperature_k(DESTINATIONS[1].cell).map_err(|e| e.to_string())?)
                        .map_err(|e| e.to_string())?;
                    let contact = SeriesElectricalContact::new(old.id(), old.branch(),
                        old.source().at_query_point(&interface, &point).map_err(|e| e.to_string())?, old.area_m2())
                        .map_err(|e| e.to_string())?;
                    Ok::<_, String>(ContactThermalCoefficients { bulk, contacts: vec![contact] })
                }, |_, electrical, mut transport| {
                    native_trials.set(native_trials.get() + electrical.iterations);
                    resolutions.set(resolutions.get() + 1);
                    if refuse_at.get() == Some(resolutions.get()) { return Err("late transport refusal".to_owned()); }
                    // Fourier laws are scalar constants in this fixture. Retain
                    // their original coefficients while forwarding solved-source observations.
                    for (i, destination) in DESTINATIONS.iter().enumerate() {
                        transport.payload.cells.iter_mut().find(|s| s.cell == destination.cell).unwrap()
                            .conductivity = Some(electrical.evaluated[i].conductor.material().clone());
                    }
                    Ok(transport)
                }, |_| 0.0, cancelled)
        };
        let accepted = motion.clone(); let cells = thermal.cells().clone();
        let policy = iteration().with_relaxation(0.6).unwrap();
        let polls = Cell::new(0usize);
        let normal = run(config(), policy, &mut || { polls.set(polls.get() + 1); false }).unwrap();
        let source_calls = resolutions.get();
        assert_eq!(normal.report().coupled.circuit_trials, native_trials.get());
        assert!(native_trials.get() > source_calls);
        let electrical = &normal.report().coupled.electrical;
        let field = &electrical.field;
        let temperatures = DESTINATIONS.map(|d| normal.thermal().temperature_k(d.cell).unwrap());
        let current = field.contacts[0].current_a.abs();
        // Independent implicit series laws: 3 = I[R_bulk(I,T) + R_contact(I,T)].
        let rb = (1.0 + 0.004 * (temperatures[0] - 300.0)) / (1.0 - 0.05 * current);
        let rc = 0.5 * (1.0 + 0.002 * (temperatures[0] - 300.0)) / (1.0 - 0.05 * current);
        near(3.0, current * (rb + rc));
        near(electrical.evaluated[0].conductor.resistance_ohm(), rb);
        near(electrical.evaluated_contacts[0].contact.resistance_ohm(), rc);
        for (i, destination) in DESTINATIONS.iter().enumerate() {
            let old = &thermal.cells()[&destination.cell]; let new = &normal.thermal().cells()[&destination.cell];
            let partition = if i == 0 { 0.25 } else { 0.75 };
            near(old.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg), field.resistors[i].heat_j
                + partition * field.contacts[0].heat_j + normal.report().coupled.source.layers[i].reversible_heat_j
                + 0.05 * 4.0 * (temperatures[1 - i] - temperatures[i]));
        }
        near(total_internal_layers(normal.thermal(), normal.motion(), [1e-3, 2e-3]),
            total_internal_layers(&thermal, &motion, [1e-3, 2e-3]) + field.electrical_supply_j);
        near(normal.report().coupled.deposited_joule_j, field.heat_j);
        near(normal.report().total_energy_residual_j, 0.0);
        native_trials.set(0); resolutions.set(0); refuse_at.set(Some(source_calls));
        assert!(matches!(run(config(), policy, &mut || false),
            Err(CircuitEnthalpyError::Source(SupportThermalError::Heat(message))) if message == "late transport refusal"));
        native_trials.set(0); resolutions.set(0); refuse_at.set(None);
        assert!(matches!(run(CircuitEnthalpySolve { maximum_circuit_trials:
            normal.report().coupled.circuit_trials - 1, ..config() }, policy, &mut || false),
            Err(CircuitEnthalpyError::CircuitBudget)));
        assert!(matches!(run(config(), ConductivitySolve { maximum_iterations: 1, ..iteration() }.into(),
            &mut || false), Err(CircuitEnthalpyError::Circuit(MaterialCircuitError::ConductivityConvergence { .. }))));
        let actual = Cell::new(0usize);
        assert!(matches!(run(config(), policy, &mut || { actual.set(actual.get() + 1); actual.get() == polls.get() }),
            Err(CircuitEnthalpyError::Cancelled)));
        let retry = run(config(), policy, &mut || false).unwrap();
        assert_eq!(retry.thermal().cells(), normal.thermal().cells());
        assert_eq!(retry.motion(), normal.motion());
        assert_eq!(thermal.cells(), &cells); assert_eq!(thermal.time_s(), 0.0); assert_eq!(motion, accepted);
    });
}

#[test]
fn electrically_independent_response_reduces_to_complete_heat_support_path() {
    with_heat_support(0.2, 0.2, |cx, _, assembly, binding, thermal, motion| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let bulk = assembly.electrical.bulk;
        let responses = bulk.resistors.iter().map(|resistor| ConductorResponse {
            branch: resistor.destination.branch,
            card: bulk.cells.iter().find(|cell| cell.cell == resistor.destination.cell).unwrap().point.material,
            axes: &[],
        }).collect::<Vec<_>>();
        let routes = prepared.heat_destinations(&thermal).unwrap();
        let actual = step_conductivity_candidate(bulk, &thermal, 0.05, &motion, &[0.0, 3.0], 0.0,
            &responses, None, &routes, binding, config(), iteration(), &[0.0; 2],
            |trial| prepared.resolve(cx, trial), |_, _, transport| Ok(transport),
            |sample| sample.radiation_net_power_w(), || cx.checkpoint().is_err()).unwrap();
        let original = prepared.field_motion_support_candidate(cx, &thermal, 0.05, &motion,
            &[0.0, 3.0], 0.0, binding, config(), &[0.0; 2]).unwrap();
        assert_eq!(actual.motion(), original.motion());
        assert_eq!(actual.thermal().cells(), original.thermal().cells());
        assert_eq!(actual.report().external_heat_j, original.report().external_heat_j);
        assert_eq!(actual.report().coupled.circuit_trials, original.report().coupled.circuit_trials);
        assert_eq!(actual.report().coupled.electrical.iterations, 1);
    });
}

#[test]
fn electrical_applicability_cannot_freeze_support_or_caloric_sources_at_seed_field() {
    for property in [YOUNG_MODULUS_PROPERTY, "specific_heat_capacity"] {
        with_response_support(0.05, None, Some(property),
            |_, prepared, responses, _, binding, thermal, motion| {
            let saved = motion.clone(); let cells = thermal.cells().clone();
            assert!(matches!(prepared.field_motion_support_conductivity_candidate(&thermal, 0.05,
                &motion, &[0.0, 3.0], 0.0, responses, binding, config(), iteration(), &[0.0; 2], || false),
                Err(CircuitEnthalpyError::Source(SupportThermalError::Invalid(_)))));
            assert_eq!(motion, saved); assert_eq!(thermal.cells(), &cells); assert_eq!(thermal.time_s(), 0.0);
        });
    }
}
