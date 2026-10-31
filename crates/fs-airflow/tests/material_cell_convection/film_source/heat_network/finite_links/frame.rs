//! Actual finite coolant, native C/R storage and live material frame motion.
//! Manufactured SI specimen; source substitution is not experimental validation.
use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::HeatAssembly;
use fs_couple::frame_vibration::{
    FrameVibrationBasis, FrameVibrationBudget, FrameVibrationComponent,
    thermal::{FrameLayerCell, FrameThermalSolve},
    thermoelastic::{FrameElectricalInput, FrameHeatError},
};
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_couple::render::schedule::force::coupled::ModalCouplingConfig;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{DynViscosity, Pressure};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, FrameViscousSection, Geometry};

const LENGTH: f64 = 0.1;
const AREA: f64 = 1e-5;
const ETA: f64 = 0.2;
const E_SLOPE: f64 = -20.0;

fn frame_point(t: f64) -> QueryPoint {
    body_point(t).with("omega", 100.0).unwrap()
}

fn frame_material() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = body_material();
    let mut claims = ClaimSet::new();
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        // Fixed original inventory requires source-declared scalar density.
        if claim.key.name() == DENSITY_PROPERTY {
            claim.value = PropertyValue::Scalar { value: 1000.0, dims: Density::DIMS };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, 1000.0),
        ("eta", DynViscosity::DIMS, ETA),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS, 0.0),
    ] {
        let curved = name == "young_modulus";
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity("T", temperature_quantity(), 280.0, 320.0);
        if name == "eta" { validity = validity.with("omega", 0.0, 1000.0); }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if curved { PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, 1600.0), (320.0, 800.0)], dims,
            } } else { PropertyValue::Scalar { value, dims } },
            validity,
            interpolation: if curved { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured finite-coolant material frame".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS), temperature_axis: "T",
        lower: &frame_point(280.0), upper: &frame_point(320.0),
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, curve)
}

fn frame_component(card: &MaterialCard, point: &QueryPoint) -> FrameVibrationComponent {
    let requirements = [
        (DENSITY_PROPERTY, Density::DIMS), ("young_modulus", Pressure::DIMS),
        ("eta", DynViscosity::DIMS),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS),
    ].map(|(name, dims)| ScalarPropertyRequirement::try_new(
        name, dims, ScalarAdmissibility::Finite).unwrap());
    let section = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: (AREA / core::f64::consts::PI).sqrt(),
        material: resolve_material_state_point(card, point, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap(),
    }]).unwrap();
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(LENGTH, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame.add_material_viscous_element(0, 1, 5,
        FrameViscousSection::try_new(section, "eta", "eta").unwrap()).unwrap();
    FrameVibrationBasis::from_frame(&frame, 1,
        FrameVibrationBudget { max_frame_dofs: 6, max_dense_entries: 36 }, || false)
        .unwrap().coupled_component(20, &[0.0], ModalAcousticTimeBudget::audible_reference(),
            18, || false).unwrap()
}

fn coolant(cp: f64, limited: bool) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, curve) = fluid_material(1.0, cp);
    if !limited { return (original, curve); }
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY {
            // The initial film is 300 K. Physical heat exchange subsequently
            // takes it below this source's original declared support.
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("T", temperature_quantity(), 300.0, 320.0)
                .with("p", 100000.0, 100000.0)
                .with("solute_mass_fraction", 0.2, 0.2);
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS), temperature_axis: "T",
        lower: &fluid_point(280.0), upper: &fluid_point(320.0),
        selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Liquid, 0.0, 1e-6).unwrap();
    (card, curve)
}

fn frame_controls() -> FrameThermalSolve {
    FrameThermalSolve {
        maximum_cells: 2, maximum_material_trials: 500, maximum_newton_iterations: 24,
        difference_step_j: 1e-5, cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8, total_energy_tolerance_j: 2e-8,
        newton: NewtonKrylovConfig { linear_restart: 2, ..controls().newton },
    }
}

fn coupling() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 1, max_connections: 6, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 100.0, maximum_abs_connection_force_n: 1000.0,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}

// Independent two-capacity, film-property, midpoint-RC equations. The first
// mechanical interval has zero displacement/velocity, so elastic work is zero.
fn reference(cp: f64, powered: bool) -> [f64; 3] {
    let evaluate = |wall: f64| {
        let resistance = 1.0 + 0.02 * (wall - 300.0);
        let ratio = DT / (2.0 * CAPACITANCE * resistance);
        let charge = INITIAL_CHARGE * (1.0 - ratio) / (1.0 + ratio);
        let joule = if powered { (INITIAL_CHARGE.powi(2) - charge.powi(2)) / (2.0 * CAPACITANCE) } else { 0.0 };
        let fluid = 290.0 + (joule - MASS * CP * (wall - 310.0)) / (MASS * cp);
        let film = 0.5 * (wall + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let mu = 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * LENGTH / mu).sqrt() * (mu * cp / k).cbrt() * k / LENGTH;
        (MASS * CP * (wall - 310.0) - joule + DT * SURFACE * h * (wall - fluid),
            [wall, fluid, charge])
    };
    let (mut low, mut high) = (290.0, 312.0);
    assert!(evaluate(low).0 < 0.0 && evaluate(high).0 > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if evaluate(middle).0 > 0.0 { high = middle; } else { low = middle; }
    }
    let (residual, expected) = evaluate(0.5 * (low + high));
    assert!(residual.abs() < 1e-10);
    expected
}

#[test]
fn finite_coolant_reaches_powered_and_unpowered_native_frame_histories() {
    for powered in [false, true] {
        let mut displacements = Vec::new();
        for (cp, limited) in [(2000.0, false), (8000.0, false), (2000.0, true)] {
            let (body_card, body_curve) = frame_material();
            let (fluid_card, fluid_curve) = coolant(cp, limited);
            let body_context = frame_point(310.0);
            let fluid_context = fluid_point(290.0);
            let fluid = liquid_law(&fluid_card);
            let mut thermal = initial_pair(&body_curve, &fluid_curve, cp, false, 0.0);
            let cells = [
                MaterialCell { cell: BODY, material_id: MaterialId(0), curve: body_curve.curve(),
                    point: ThermalMaterialPoint { material: &body_card, context: &body_context,
                        temperature_axis: "T", liquid_fraction_axis: None }, conductivity: None },
                MaterialCell { cell: FLUID, material_id: MaterialId(1), curve: fluid_curve.curve(),
                    point: ThermalMaterialPoint { material: &fluid_card, context: &fluid_context,
                        temperature_axis: "T", liquid_fraction_axis: None }, conductivity: None },
            ];
            let fourier = [];
            let resistors = [CellResistor { destination: ResistorCell { branch: 1, cell: BODY },
                selection: MaterialPropertySelection::SingleClaimOnly, length_m: LENGTH, area_m2: AREA }];
            let thermal_source = MaterialNetwork { cells: &cells, resistors: &[], links: &fourier };
            let source = MaterialNetwork { cells: &cells,
                resistors: if powered { &resistors } else { &[] }, links: &fourier };
            let bases = [(1, ElectricalPropertyBasis::Resistivity)];
            let electrical = ElectricalContactNetwork { bulk: &source, contacts: &[],
                electrical_bases: if powered { &bases } else { &[] } };
            let heat = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
            let links = [binding()];
            let initial = evaluate_source_link_convection(&thermal, &links, 1,
                |_, _, t| fluid.at_temperature(t), || false).unwrap();
            let gate = fs_exec::CancelGate::new_clock_free();
            // Ordinary cold preparation refuses a native link with no owner.
            assert!(frame_component(&body_card, &body_context).into_thermoelastic_system(
                &thermal, &thermal_source, &[FrameLayerCell { element: 0, layer: 0, cell: BODY }],
                0.02, coupling(), None, &gate).is_err());
            let mut owner = frame_component(&body_card, &body_context)
                .into_thermoelastic_system_with_link_conductances(&thermal, &thermal_source,
                    initial.conductances_w_k(), &[FrameLayerCell { element: 0, layer: 0, cell: BODY }],
                    0.02, coupling(), None, &gate).unwrap();
            let source_gate = fs_exec::CancelGate::new_clock_free();
            fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
                let cx = fs_exec::Cx::new(&source_gate, arena,
                    fs_exec::StreamKey { seed: 175, kernel_id: 0x4652414d, tile: 0, iteration: 0 },
                    fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
                let prepared = heat.prepare_with_link_conductances(&cx, &thermal,
                    initial.conductances_w_k()).unwrap();
                let network = SourceConvectionHeatNetwork { heat: &prepared, cells: &[], maximum_cells: 0,
                    resolve_fluid: |_: &SourceConvectiveCell, _: &DepositionSimulation<'_>, t: Temperature|
                        fluid.at_temperature(t) };
                let assembly = network.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
                if powered {
                    let circuit = MaterialCircuit::new(CircuitGraph { node_count: 2,
                        branches: vec![(1, 0, Branch::Capacitor { farads: CAPACITANCE }),
                            (1, 0, Branch::Resistor { ohms: 0.0 })], transformers: vec![] },
                        assembly.resolve_heat(&cx, &thermal).unwrap().bulk.resistors,
                        vec![], vec![INITIAL_CHARGE], &[]).unwrap();
                    owner = owner.with_material_circuit_and_link_conductances(&thermal, &source,
                        circuit, &bases, initial.conductances_w_k(), &gate).unwrap();
                }
                let input = powered.then_some(FrameElectricalInput { ports: &[], solve: controls().circuit });
                let original_cells = thermal.cells().clone();
                let original_links = thermal.links().clone();
                let ordinary = match input {
                    Some(input) => owner.step_electrical(&mut thermal, &[0.0], input,
                        frame_controls(), &[0.0; 2], &gate),
                    None => owner.step(&mut thermal, &[0.0], frame_controls(), &[0.0; 2], &gate),
                };
                assert!(ordinary.is_err());
                assert_eq!(thermal.cells(), &original_cells);
                assert_eq!(owner.mechanical().samples_rendered(), 0);
                let original_circuit = owner.electrical().cloned();
                let original_motion = owner.mechanical().components()[0].states()[0];
                let original_member = owner.members()[0].source.clone();
                let result = owner.step_heat_assembly(&mut thermal, &[0.0], &assembly, &cx,
                    input, frame_controls(), &[0.0; 2], &gate);
                if limited {
                    assert!(matches!(result,
                        Err(FrameHeatError::Heat(ConvectionHeatError::Convection(_)))));
                    assert_eq!(thermal.cells(), &original_cells);
                    assert_eq!(thermal.links(), &original_links);
                    assert_eq!(owner.electrical(), original_circuit.as_ref());
                    assert_eq!(owner.members()[0].source, original_member);
                    assert_eq!(owner.mechanical().components()[0].states()[0], original_motion);
                    assert_eq!(owner.mechanical().samples_rendered(), 0);
                    assert_eq!(thermal.time_s(), 0.0);
                    return;
                }
                let report = result.unwrap();
                let expected = reference(cp, powered);
                near(thermal.temperature_k(BODY).unwrap(), expected[0], 2e-8);
                near(thermal.temperature_k(FLUID).unwrap(), expected[1], 2e-8);
                let sample = &report.source.finite_convection.as_ref().unwrap().samples()[0];
                assert_eq!((sample.wall_cell, sample.fluid_cell), (BODY, FLUID));
                assert_eq!(sample.material.material().card_identity(), fluid_card.content_hash());
                near(sample.property_temperature_k, 0.5 * (expected[0] + expected[1]), 2e-8);
                near(MASS * cp * (expected[1] - 290.0), DT * sample.heat_into_fluid_w, 2e-8);
                for property in sample.material.material().properties() {
                    let original = fluid.reference_state().material()
                        .property(property.requirement().name()).unwrap();
                    assert_eq!(property.answer().receipt.selected, original.answer().receipt.selected);
                }
                assert!(report.heat_transfer.is_none());
                assert_eq!(thermal.links(), &original_links);
                assert_eq!(owner.mechanical().samples_rendered(), 1);
                near(report.external_thermal_heat_j, 0.0, 2e-8);
                near(report.total_energy_residual_j, 0.0, 2e-8);
                near(report.members[0].source.layers[0].young_modulus_pa,
                    1000.0 + E_SLOPE * (expected[0] - 310.0), 1e-6);
                if powered {
                    near(owner.electrical().unwrap().charges_c()[0], expected[2], 1e-9);
                    near(report.deposited_joule_heat_j,
                        (INITIAL_CHARGE.powi(2) - expected[2].powi(2)) / (2.0 * CAPACITANCE), 2e-8);
                } else {
                    assert!(report.electrical.is_none());
                    assert_eq!(report.deposited_joule_heat_j, 0.0);
                }
                // The accepted coolant endpoint changes the native oscillator
                // in the NEXT interval. Compare its exact forced damped motion.
                let shape = owner.basis().modes()[0].displacement[3];
                let force = 1e-5;
                let next = owner.step_heat_assembly(&mut thermal, &[shape * force], &assembly,
                    &cx, input, frame_controls(), &[0.0; 2], &gate).unwrap();
                let k = AREA / LENGTH * (1000.0 + E_SLOPE * (expected[0] - 310.0));
                let mass = MASS / 3.0;
                let decay = ETA * AREA / LENGTH / (2.0 * mass);
                let omega = (k / mass - decay * decay).sqrt();
                let expected_x = force / k * (1.0 - (-decay * DT).exp()
                    * ((omega * DT).cos() + decay / omega * (omega * DT).sin()));
                let motion = owner.mechanical().components()[0].states()[0];
                near(shape * motion.displacement_m_sqrt_kg, expected_x, 1e-10);
                near(next.total_energy_residual_j, 0.0, 2e-8);
                displacements.push(expected_x);

                let cells_before = thermal.cells().clone();
                let circuit_before = owner.electrical().cloned();
                let members_before = owner.members()[0].source.clone();
                // Distinct execution authorities must cancel together before
                // native mechanical or caloric publication.
                source_gate.request();
                let failure = owner.step_heat_assembly(&mut thermal, &[0.0], &assembly, &cx,
                    input, frame_controls(), &[0.0; 2], &gate);
                assert!(matches!(failure, Err(FrameHeatError::Frame(CircuitEnthalpyError::Cancelled))));
                assert!(gate.is_requested());
                assert_eq!(thermal.cells(), &cells_before);
                assert_eq!(owner.electrical(), circuit_before.as_ref());
                assert_eq!(owner.members()[0].source, members_before);
                assert_eq!(owner.mechanical().components()[0].states()[0], motion);
                assert_eq!(owner.mechanical().samples_rendered(), 2);
                near(thermal.time_s(), 2.0 * DT, 0.0);
            });
        }
        assert!((displacements[0] - displacements[1]).abs() > 1e-8);
    }
}
