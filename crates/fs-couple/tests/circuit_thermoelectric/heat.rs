//! G1/G3/G4: actual Kelvin conversion and complete shared heat transport.
use super::*;
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_conduction::{AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP, EMISSIVITY_DIMS,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity};
use fs_couple::circuit_enthalpy::contacts::{heat::{MaterialHeatNetwork, PreparedHeatNetwork}, material::*};
use fs_exec::{CancelGate, Cx};
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_DIMS as ED,
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as EP};
use fs_material::thermal::conduction::{LayerConductionConfig,
    network::{LayerLinkBinding, PreparedLayerConduction}, source::ResolvedConductionLayer};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_matdb::{ClaimSelection, ClaimSet, InterfaceSystemCard, SurfaceSpec, SystemContext};

#[path = "../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiation_geometry;

#[path = "heat/motion.rs"]
mod motion;

fn scoped<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| run(&Cx::new(&gate, arena,
        fs_exec::StreamKey { seed: 29, kernel_id: 61, tile: 0, iteration: 0 },
        fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic), &gate))
}

fn radiating_sources(e: f64) -> ([MaterialCard; 4], [ResolvedHeatCapacityCurve; 4]) {
    let (original, _) = fixture(0.02, -0.01, 3.0, 1.0);
    let cards = std::array::from_fn(|i| {
        let mut claims = ClaimSet::new();
        for (_, claim) in original[i].claims().claims_ordered() {
            let mut claim = claim.clone();
            if i < 2 && claim.key.name() == "thermal_conductivity" {
                claim.value = PropertyValue::Curve { abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, 0.5), (2000.0, 2.3)], dims: THERMAL_CONDUCTIVITY_DIMS };
            }
            claims.insert_claim(claim).unwrap();
        }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar { value: e, dims: EMISSIVITY_DIMS },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance {
                source: "manufactured coupled heat fixture, not measured thermoelectric material".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        MaterialCard::assemble(original[i].id().clone(), claims, vec![]).unwrap()
    });
    let curves = std::array::from_fn(|i| ResolvedHeatCapacityCurve::try_new(&cards[i],
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
            lower: &point(200.0), upper: &point(2000.0), selection: SelectionPolicy::SingleClaimOnly,
        }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap());
    (cards, curves)
}

fn interface(cards: &[MaterialCard; 4], rt: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [(RP, RD, rt), (EP, ED, 0.5)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("aT", 200.0, 2000.0)
                .with("bT", 200.0, 2000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance {
                source: "manufactured Ohmic/Fourier load interface".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec { material: cards[0].id().clone(), texture_frame: "a".into() },
        SurfaceSpec { material: cards[1].id().clone(), texture_frame: "b".into() },
        SystemContext { medium: "fixed".into(), third_body: None,
            environment: "manufactured".into(), history: "fixed".into() }, claims, vec![]).unwrap()
}

// Extra link 4 loads the two junctions; it does not replace either Kelvin leg.
fn with_heat<R>(e: f64, rt: f64, layers: bool, current: Option<f64>, thermal_link: u64,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &PreparedHeatNetwork<'_, '_>,
        &DepositionSimulation<'_>, &MaterialCircuit) -> R) -> R {
    let (cards, curves) = radiating_sources(e);
    let context = point(350.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let mut resistors = resistors().into_iter().collect::<Vec<_>>();
    resistors.push(CellResistor { destination: ResistorCell { branch: 5, cell: CellId(0) },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: 1.0 });
    let mut links = links().to_vec();
    links.push(ConductionLink { link: 4, distance_a_m: 0.5, distance_b_m: 0.5 });
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial_with_links(&curves, [400.0, 300.0, 350.0, 350.0], &[[0, 1]]);
    let interface = interface(&cards, rt);
    let ip = QueryPoint::new().with("aT", 400.0).unwrap().with("bT", 300.0).unwrap();
    let point = ThermalInterfacePoint { interface: &interface, materials: [&cards[0], &cards[1]],
        context: &ip, temperature_axes: ["aT", "bT"], liquid_fraction_axes: [None; 2],
        mean_temperature_axis: None };
    let contacts = [MaterialElectricalContact { link: 4, branch: 5, point,
        selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: 0.25 }];
    let bases = [BASES[0], BASES[1], (5, ElectricalPropertyBasis::Resistivity)];
    let electrical = ElectricalContactNetwork { bulk: &material, contacts: &contacts, electrical_bases: &bases };
    let other = if thermal_link == 0 { 2 } else { 1 };
    let thermal_contacts = [MaterialContactLink { link: thermal_link,
        point: ThermalInterfacePoint { materials: [&cards[0], &cards[other]], ..point },
        sides: [0, other].map(|i| ContactCell { cell: CellId(i as u64),
            material_id: MaterialId(i as u64), curve: curves[i].curve() }),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) }];
    let surfaces = [0, 1].map(|i| RadiatingCell { surface: i, cell: CellId(i as u64),
        material_id: MaterialId(i as u64), curve: curves[i].curve(), material: cells[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly });
    let emissivities = [0, 1].map(|i| SurfaceEmissivity::from_thermal_state("manufactured",
        cells[i].point, curves[i].curve().state_at_specific_enthalpy(if i == 0 { 2000.0 } else { 1000.0 }).unwrap(),
        MaterialPropertySelection::SingleClaimOnly).unwrap().0);
    let enclosure = radiation_geometry::enclosure(emissivities);
    let radiation = MaterialRadiation { enclosure: &enclosure, cells: &surfaces,
        power_tolerance_w: 1e-8, radiosity_tolerance_w_m2: 1e-8 };
    let assembly = MaterialHeatNetwork { electrical: &electrical,
        thermal_contacts: if layers { &[] } else { &thermal_contacts }, radiation: Some(&radiation) };
    let sources = electrical.resolve(&thermal, || false).unwrap();
    let circuit = MaterialCircuit::new_with_contacts(CircuitGraph { node_count: 6, branches: vec![
        (1, 2, Branch::VoltageSource { port: 0 }), (2, 3, Branch::Resistor { ohms: 999.0 }),
        (3, 4, Branch::VoltageSource { port: 1 }), (4, 0, Branch::Resistor { ohms: 999.0 }),
        (1, 5, if current.is_some() { Branch::CurrentSource { port: 2 } } else { Branch::VoltageSource { port: 2 } }),
        (5, 0, Branch::Resistor { ohms: 999.0 }),
    ], transformers: vec![] }, sources.bulk.resistors, sources.contacts,
        vec![], vec![], &[0.0, 0.0, current.unwrap_or(0.0)]).unwrap();
    let source_layers = [0, 1].map(|i| ResolvedConductionLayer::new(curves[i].curve(),
        cells[i].point, &key, ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid, [200.0, 2000.0], || Ok(())).unwrap());
    let layer_model = PreparedLayerConduction::new(&thermal, vec![LayerLinkBinding {
        link: 4, layers: [&source_layers[0], &source_layers[1]], distances_m: [0.5; 2],
        contact_area_resistance_m2_k_w: rt }], 1, || Ok(())).unwrap();
    scoped(|cx, gate| {
        let mut prepared = assembly.prepare(cx, &thermal).unwrap();
        if layers { prepared = prepared.with_integrated_conduction(cx, &thermal, &layer_model,
            LayerConductionConfig { max_iterations: 80, absolute_heat_tolerance_w: 1e-10,
                relative_heat_tolerance: 1e-12 }).unwrap(); }
        run(cx, gate, &prepared, &thermal, &circuit)
    })
}

#[test]
fn full_heat_changes_generator_current_and_preserves_original_interface_energy() {
    let run = |e, rt| with_heat(e, rt, false, None, 4, |cx, _, prepared, thermal, circuit| {
        let next = prepared.step_thermoelectric_candidate(cx, thermal, 0.01, circuit,
            &[0.0; 3], &legs(), config(), &[0.0; 4]).unwrap();
        let report = next.report();
        let [a, b] = [CellId(0), CellId(1)].map(|id| next.thermal().temperature_k(id).unwrap());
        let current = report.coupled.source.legs[0].transfer.current_a;
        near(current, 0.03 * (a - b) / 6.5);
        near(report.coupled.deposited_joule_j, 0.01 * current * current * 6.5);
        near(report.coupled.electrical.contacts[0].heat_j, 0.01 * current * current * 0.5);
        assert_eq!(report.coupled.source.network.thermal_contacts.len(), 1);
        for surface in report.coupled.source.network.radiation.as_ref().unwrap().surfaces() {
            let bulk = report.coupled.source.network.material.cells.iter().find(|x| x.cell == surface.cell).unwrap();
            assert_eq!(surface.phase, bulk.phase);
        }
        assert!((next.thermal().energy_j() - thermal.energy_j()).abs() < config().total_energy_tolerance_j);
        assert!(report.external_heat_j.abs() < 1e-9);
        assert!(report.total_energy_residual_j.abs() < config().total_energy_tolerance_j);
        (current, a - b)
    });
    let baseline = run(0.2, 0.5);
    for changed in [run(0.8, 0.5), run(0.2, 3.0)] {
        assert!((changed.0 - baseline.0).abs() > 1e-7);
        assert!((changed.1 - baseline.1).abs() > 1e-5);
    }
}

#[test]
fn powered_layers_retain_kelvin_reversal_and_deposit_contact_heat_once() {
    for current in [-0.5, 0.5] {
        with_heat(0.3, 0.0, true, Some(current), 4, |cx, _, prepared, thermal, circuit| {
            let next = prepared.step_thermoelectric_candidate(cx, thermal, 0.01, circuit,
                &[0.0, 0.0, current], &legs(), config(), &[0.0; 4]).unwrap();
            let report = next.report();
            let [a, b] = [CellId(0), CellId(1)].map(|id| next.thermal().temperature_k(id).unwrap());
            let kelvin = &report.coupled.source.legs;
            near(kelvin[0].transfer.current_a, current);
            near(kelvin[0].transfer.terminal_heat_w[0] + kelvin[1].transfer.terminal_heat_w[1], -0.03 * a * current);
            near(kelvin[0].transfer.terminal_heat_w[1] + kelvin[1].transfer.terminal_heat_w[0], 0.03 * b * current);
            let layer = &report.coupled.source.network.integrated_conduction[0];
            near(layer.applied_power_a_to_b_w, (0.3 + 0.001 * (a + b) * 0.5) * (a - b));
            near(report.coupled.deposited_joule_j, 0.01 * current * current * 6.5);
            near(report.external_electrical_work_j, 0.01 * current * (current * 6.5 + 0.03 * (b - a)));
            assert!((next.thermal().energy_j() - thermal.energy_j() - report.external_electrical_work_j).abs()
                < config().total_energy_tolerance_j);
        });
    }
}

#[test]
fn absent_additional_heat_matches_original_and_cancellation_is_atomic() {
    let (cards, curves) = fixture(0.02, -0.01, 3.0, 1.0);
    let context = point(350.0); let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key); let rs = resistors(); let links = links();
    let bulk = MaterialNetwork { cells: &cells, resistors: &rs, links: &links };
    let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
    let circuit = circuit(&bulk, &thermal, false, 0.0);
    let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[], electrical_bases: &BASES };
    let heat = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
    scoped(|cx, gate| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let original = thermoelectric::step_candidate(&bulk, &thermal, 0.1, &circuit,
            &[0.0; 3], &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
        let actual = prepared.step_thermoelectric_candidate(cx, &thermal, 0.1, &circuit,
            &[0.0; 3], &legs(), config(), &[0.0; 4]).unwrap();
        assert_eq!(actual.thermal().cells(), original.thermal().cells());
        assert_eq!(actual.circuit(), original.circuit());
        let saved = thermal.cells().clone(); let charge = circuit.clone();
        gate.request();
        assert!(matches!(prepared.step_thermoelectric_candidate(cx, &thermal, 0.1, &circuit,
            &[0.0; 3], &legs(), config(), &[0.0; 4]), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.cells(), &saved); assert_eq!(circuit, charge);
    });
}

#[test]
fn unresolved_contact_temperature_jump_on_a_kelvin_edge_refuses_atomically() {
    with_heat(0.2, 0.5, false, None, 0, |cx, _, prepared, thermal, circuit| {
        let saved = thermal.cells().clone(); let charge = circuit.clone();
        assert!(matches!(prepared.step_thermoelectric_candidate(cx, thermal, 0.01, circuit,
            &[0.0; 3], &legs(), config(), &[0.0; 4]), Err(CircuitEnthalpyError::Invalid(
                "a Kelvin edge requires its original junction-to-bulk temperature trace"))));
        assert_eq!(thermal.cells(), &saved); assert_eq!(*circuit, charge);
    });
}
