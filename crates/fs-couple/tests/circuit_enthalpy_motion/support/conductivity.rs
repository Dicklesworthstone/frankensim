//! Nonlinear source response through the original live support/enthalpy owners.
use super::*;
use fs_couple::circuit_enthalpy::contacts::ContactHeatCells;
use fs_couple::circuit_enthalpy::contacts::material::ContactConductionResponse;
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
            let limit = if Some(claim.key.name()) == bounded { 0.1 } else { 6.0 };
            claim.value = PropertyValue::Surface2D { axes: ["T".into(), "E".into()],
                axis_dims: [Temperature::DIMS, FIELD],
                samples: [[280.0, 0.0], [400.0, 0.0], [400.0, limit], [280.0, limit]]
                    .map(|[t, e]| ([t, e], 1.0 + thermal * (t - 300.0) + field * e)).to_vec(),
                triangles: vec![[0, 1, 2], [0, 2, 3]], dims: claim.key.quantity().dims() };
            claim.validity = claim.validity.with_quantity("E", QuantitySpec::dimensional(FIELD), 0.0, limit);
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
        let bulk_responses = bulk.resistors.iter().map(|resistor| CellConductionResponse {
            branch: resistor.destination.branch, basis: ElectricalPropertyBasis::Resistivity,
            axes: &[],
        }).collect::<Vec<_>>();
        let contacts = assembly.electrical.contacts.iter().map(|contact| ContactConductionResponse {
            id: contact.link, direction: ContactDirection::AlongBranch, axes: &[],
        }).collect::<Vec<_>>();
        let complete = assembly.prepare_conductivity(cx, &thermal, &bulk_responses, &contacts, &[]).unwrap();
        let actual = complete.field_motion_support_conductivity_candidate(cx, &thermal, 0.05,
            &motion, &[0.0, 3.0], 0.0, &bulk_responses, &contacts, binding, config(), iteration(),
            &[0.0; 2]).unwrap();
        assert_eq!(actual.motion(), original.motion());
        assert_eq!(actual.thermal().cells(), original.thermal().cells());
        assert_eq!(actual.report().external_heat_j, original.report().external_heat_j);
        assert_eq!(actual.report().coupled.circuit_trials, original.report().coupled.circuit_trials);
        assert_eq!(actual.report().coupled.electrical.iterations, 1);
    });
}

fn with_complete_response_support<R>(bounded: Option<&str>, contact_bound: bool, resistance: f64,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &MaterialHeatNetwork<'_>, &[CellConductionResponse<'_>],
        &[ContactConductionResponse<'_>], &SupportThermalBinding<'_>, DepositionSimulation<'_>,
        MaterialFieldMotion) -> R) -> R
{
    with_response_support(0.05, Some(0.3), bounded, |bulk, _, responses, cards, binding, thermal, original| {
        let response = source::response_card(&cards[0], &cards[1], "V", VOLTAGE, 0.5, 0.002, 0.1, 6.0);
        let thermal_card = support_interface(cards, resistance);
        let mut claims = response.claims().clone();
        for (_, claim) in thermal_card.claims().claims_ordered() {
            if claim.key.name() == CONTACT_RP {
                let mut claim = claim.clone();
                if contact_bound {
                    claim.validity = claim.validity.with_quantity("V", QuantitySpec::dimensional(VOLTAGE), 0.0, 0.1);
                }
                claims.insert_claim(claim).unwrap();
            }
        }
        let card = InterfaceSystemCard::assemble(response.surface_a().clone(), response.surface_b().clone(),
            response.context().clone(), claims, vec![]).unwrap();
        let context = source::point().with_quantity("V", QuantitySpec::dimensional(VOLTAGE), 0.0).unwrap();
        let point = ThermalInterfacePoint { interface: &card, materials: [&cards[0], &cards[1]],
            context: &context, temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [Some("melt_a"), Some("melt_b")], mean_temperature_axis: None };
        let contacts = [MaterialElectricalContact { link: 7, branch: 1, point,
            selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: 0.25 }];
        let axes = [ContactAxis { name: "V", coordinate: ContactCoordinate::Voltage }];
        let contact_responses = [ContactConductionResponse { id: 7,
            direction: ContactDirection::AgainstBranch, axes: &axes }];
        let bases = bases();
        let electrical = ElectricalContactNetwork { bulk, contacts: &contacts, electrical_bases: &bases };
        let thermal_contacts = [MaterialContactLink { link: 7, point,
            sides: [0, 1].map(|i| ContactCell { cell: DESTINATIONS[i].cell,
                material_id: bulk.cells[i].material_id, curve: bulk.cells[i].curve }),
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) }];
        let surfaces = [0, 1].map(|i| RadiatingCell { surface: i, cell: DESTINATIONS[i].cell,
            material_id: bulk.cells[i].material_id, curve: bulk.cells[i].curve, material: bulk.cells[i].point,
            selection: MaterialPropertySelection::SingleClaimOnly });
        let enclosure = radiation_common::enclosure([0, 1].map(|i|
            SurfaceEmissivity::from_thermal_state("nonlinear support enclosure", bulk.cells[i].point,
                bulk.cells[i].curve.state_at_specific_enthalpy(
                    thermal.cells()[&DESTINATIONS[i].cell].specific_h_j_kg).unwrap(),
                MaterialPropertySelection::SingleClaimOnly).unwrap().0));
        let radiation = MaterialRadiation { enclosure: &enclosure, cells: &surfaces,
            power_tolerance_w: 1e-8, radiosity_tolerance_w_m2: 1e-8 };
        let assembly = MaterialHeatNetwork { electrical: &electrical,
            thermal_contacts: &thermal_contacts, radiation: Some(&radiation) };
        let prepared = electrical.prepare_conductivity(&thermal, responses, &contact_responses, &[], || false).unwrap();
        let seeds = prepared.resolve_conductivity_seeds(&thermal, responses, &contact_responses, || false).unwrap();
        let circuit = driven_circuit(seeds.bulk.resistors).rebind_contacts_candidate(seeds.contacts, &[0.0, 3.0]).unwrap();
        let motion = MaterialFieldMotion::new_with_support(circuit, 0, original.body().clone(),
            original.support().unwrap().clone(), original.velocity_m_s(), [0.9, 1.1]).unwrap();
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(&gate, arena, fs_exec::StreamKey { seed: 73, kernel_id: 51, tile: 0, iteration: 0 },
                fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
            run(&cx, &gate, &assembly, responses, &contact_responses, binding, thermal, motion)
        })
    })
}

#[test]
fn complete_nonlinear_support_closes_interface_fourier_radiation_and_mechanical_energy() {
    let run = |resistance| with_complete_response_support(None, false, resistance,
        |cx, _, assembly, responses, contacts, binding, mut thermal, mut motion| {
        assert!(assembly.prepare(cx, &thermal).is_err(), "static source admission must stay static");
        let prepared = assembly.prepare_conductivity(cx, &thermal, responses, contacts, &[]).unwrap();
        let initial = total_internal_layers(&thermal, &motion, [1e-3, 2e-3]);
        let mut supplied = 0.0;
        for n in 1..=4 {
            let target = f64::from(n) * 0.05; let dt = target - thermal.time_s();
            let next = prepared.field_motion_support_conductivity_candidate(cx, &thermal, target,
                &motion, &[0.0, 3.0], 0.0, responses, contacts, binding, config(),
                iteration().with_relaxation(0.6).unwrap(), &[0.0; 2]).unwrap();
            let report = next.report(); let electrical = &report.coupled.electrical;
            let network = &report.coupled.source.network;
            let temperatures = DESTINATIONS.map(|d| next.thermal().temperature_k(d.cell).unwrap());
            let current = electrical.field.contacts[0].current_a.abs();
            // Independent series laws for bulk E=I*Rb and the oppositely
            // oriented interface V=-I*Rc. The contact never receives branch V.
            let rb = (1.0 + 0.004 * (temperatures[0] - 300.0)) / (1.0 - 0.05 * current);
            let rc = 0.5 * (1.0 + 0.002 * (temperatures[0] - 300.0)) / (1.0 + 0.05 * current);
            near(3.0, current * (rb + rc));
            near(electrical.evaluated[0].conductor.resistance_ohm(), rb);
            near(electrical.evaluated_contacts[0].contact.resistance_ohm(), rc);
            near(electrical.evaluated_contacts[0].point.voltage_v, -current * rc);
            let k = 1.0 + 0.002 * (temperatures[0] - 300.0) + 0.3 * current * rb;
            let bulk_resistance = 0.5 / k + 0.5 / 4.0;
            let interface_resistance = resistance * (1.0 + 0.001 * (temperatures[0] - 300.0));
            let g = 1.0 / (bulk_resistance + interface_resistance);
            near(network.material.links[0].conductance_w_k, 1.0 / bulk_resistance);
            near(network.thermal_contacts[0].conductance_w_k, g);
            near(network.thermal_contacts[0].interface.value_m2_k_per_w(), interface_resistance);
            let radiation = network.radiation.as_ref().unwrap();
            let plate_power = 5.670374419e-8 * (temperatures[0].powi(4) - temperatures[1].powi(4))
                / (1.0 / 0.8 + 1.0 / 0.6 - 1.0);
            for (i, destination) in DESTINATIONS.iter().enumerate() {
                let radiative = if i == 0 { -plate_power } else { plate_power };
                near(radiation.sources_w()[&destination.cell], radiative);
                let old = &thermal.cells()[&destination.cell];
                let new = &next.thermal().cells()[&destination.cell];
                let fraction = if i == 0 { 0.25 } else { 0.75 };
                near(old.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg),
                    electrical.field.resistors[i].heat_j + fraction * electrical.field.contacts[0].heat_j
                        + report.coupled.source.layers[i].reversible_heat_j
                        + dt * (g * (temperatures[1 - i] - temperatures[i]) + radiative));
            }
            near(report.coupled.deposited_joule_j, electrical.field.heat_j);
            near(report.external_heat_j, 0.0); near(report.total_energy_residual_j, 0.0);
            supplied += electrical.field.electrical_supply_j;
            near(total_internal_layers(next.thermal(), next.motion(), [1e-3, 2e-3]), initial + supplied);
            assert!(electrical.iterations > 1);
            let (t, m, _) = next.into_parts(); thermal = t; motion = m;
        }
        [thermal.temperature_k(DESTINATIONS[0].cell).unwrap(), motion.gap_m(), motion.velocity_m_s()]
    });
    let low = run(0.2); let high = run(2.0);
    for i in 0..3 { assert!((low[i] - high[i]).abs() > 1e-9); }
}

#[test]
fn complete_nonlinear_support_refuses_unowned_sources_and_preserves_all_histories() {
    for (bounded, contact_bound) in [(Some(SURFACE_EMISSIVITY_PROPERTY), false), (None, true)] {
        with_complete_response_support(bounded, contact_bound, 0.2,
            |cx, _, assembly, responses, contacts, _, thermal, _| {
            assert!(matches!(assembly.prepare_conductivity(cx, &thermal, responses, contacts, &[]),
                Err(HeatNetworkError::Invalid(_))));
        });
    }
    with_complete_response_support(None, false, 0.2,
        |cx, gate, assembly, responses, contacts, binding, thermal, motion| {
        let prepared = assembly.prepare_conductivity(cx, &thermal, responses, contacts, &[]).unwrap();
        let saved = motion.clone(); let cells = thermal.cells().clone();
        let run = |controls| prepared.field_motion_support_conductivity_candidate(cx, &thermal, 0.05,
            &motion, &[0.0, 3.0], 0.0, responses, contacts, binding, controls, iteration(), &[0.0; 2]);
        let accepted = run(config()).unwrap();
        assert!(matches!(run(CircuitEnthalpySolve { maximum_circuit_trials:
            accepted.report().coupled.circuit_trials - 1, ..config() }), Err(CircuitEnthalpyError::CircuitBudget)));
        let retry = run(config()).unwrap();
        assert_eq!(accepted.thermal().cells(), retry.thermal().cells());
        assert_eq!(accepted.motion(), retry.motion());
        let wrong = [ContactConductionResponse { id: 999, ..contacts[0] }];
        assert!(prepared.resolve_conductivity_seeds(cx, &thermal, responses, &wrong).is_err());
        let missing = [ContactConductionResponse { axes: &[], ..contacts[0] }];
        assert!(prepared.resolve_conductivity_seeds(cx, &thermal, responses, &missing).is_err());
        gate.request();
        assert!(matches!(run(config()), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.cells(), &cells); assert_eq!(thermal.time_s(), 0.0); assert_eq!(motion, saved);
    });
    with_complete_response_support(Some(ELECTRICAL_RESISTIVITY_PROPERTY), false, 0.2,
        |cx, _, assembly, responses, contacts, binding, thermal, motion| {
        let prepared = assembly.prepare_conductivity(cx, &thermal, responses, contacts, &[]).unwrap();
        let saved = motion.clone(); let cells = thermal.cells().clone();
        assert!(matches!(prepared.field_motion_support_conductivity_candidate(cx, &thermal, 0.05,
            &motion, &[0.0, 3.0], 0.0, responses, contacts, binding, config(), iteration(), &[0.0; 2]),
            Err(CircuitEnthalpyError::Circuit(_))));
        assert_eq!(thermal.cells(), &cells); assert_eq!(thermal.time_s(), 0.0); assert_eq!(motion, saved);
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
