//! The original conducting layer, exterior solve and thermal timeline own heat.
use super::*;

fn energy(temperature: f64) -> f64 {
    let delta = temperature - 300.0;
    4000.0 * delta + delta * delta
}

/// Independent analytic strip/half-cell balance, using only original gas values.
fn reference(gas: &IsobaricGasLaw, body: f64) -> (f64, [f64; 2]) {
    let powers = |surface: f64| {
        let state = gas
            .at_temperature(Temperature::new(0.5 * (surface + 300.0)))
            .unwrap();
        let state = state.state();
        let reynolds = state.density * 0.1 / state.dynamic_viscosity;
        let prandtl = state.dynamic_viscosity * state.specific_heat_cp / state.thermal_conductivity;
        // x1/L=1/4: the strip mean is the whole-plate mean / (1+sqrt(1/4)).
        let h = 0.664 * reynolds.sqrt() * prandtl.cbrt() * state.thermal_conductivity / 0.1 / 1.5;
        [
            AREA * h * (surface - 300.0),
            AREA * 0.7 * 5.670374419e-8 * (surface.powi(4) - 305.0_f64.powi(4)),
        ]
    };
    let mut lo = 300.0;
    let mut hi = body;
    for _ in 0..80 {
        let surface = 0.5 * (lo + hi);
        let conducted = AREA / 0.005
            * (0.5 * (body - surface)
                + 0.001 * ((body - 300.0).powi(2) - (surface - 300.0).powi(2)));
        let heat = powers(surface);
        if conducted > heat[0] + heat[1] {
            lo = surface;
        } else {
            hi = surface;
        }
    }
    let surface = 0.5 * (lo + hi);
    (surface, powers(surface))
}

fn advance(gas: &IsobaricGasLaw) -> f64 {
    advance_with(gas, settings(plate()), false, reference)
}

pub(super) fn advance_with(
    gas: &IsobaricGasLaw,
    configuration: ExteriorConvectionSettings,
    corrected_air: bool,
    reference: impl Fn(&IsobaricGasLaw, f64) -> (f64, [f64; 2]),
) -> f64 {
    let card = isobaric_liquid_card(0.0, 0.01, 320.0);
    let material = isobaric_liquid_law(&card);
    let point = |temperature| {
        material
            .reference_point()
            .clone()
            .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100000.0)
            .unwrap()
            .with("T", temperature)
            .unwrap()
    };
    let cp_key = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let k_key = PropertyKey::with_quantity(
        THERMAL_CONDUCTIVITY_PROPERTY,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    );
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let caloric = IsobaricHeatCapacityLaw::try_new(
        &card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &cp_key,
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &point(280.0),
            upper: &point(320.0),
            selection,
        },
        300.0,
    )
    .unwrap();
    let conduction = IsobaricConductionLaw::try_new(&card, &caloric, &k_key, selection).unwrap();
    let layer = CaloricConductionLayer::from_isobaric(&conduction).unwrap();
    let source = ImplicitGasConvection::try_new(
        layer.caloric_law(),
        layer.energy_bounds_j_kg(),
        gas,
        configuration,
    )
    .unwrap();
    let source = if corrected_air {
        source.with_gnielinski_air_wall_correction().unwrap()
    } else {
        source
    };
    let geometry = ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: 0.005,
        surface_resistance_m2_k_w: 0.0,
    };
    let boundary = IntegratedCaloricGrayBoundary::new(
        &layer,
        geometry,
        0.0,
        0.0,
        305.0,
        LayerConductionConfig {
            max_iterations: 64,
            absolute_heat_tolerance_w: 1e-10,
            relative_heat_tolerance: 1e-12,
        },
    )
    .unwrap()
    .with_exterior_law(&source)
    .unwrap();
    let initial_temperature = 315.0;
    let result = boundary
        .evaluate(energy(initial_temperature), 300.0, &mut || Ok(()))
        .unwrap();
    let (surface, powers) = reference(gas, initial_temperature);
    near(result.surface_temperature_k, surface, 2e-9);
    near(result.convection_power_w, powers[0], 1e-8);
    near(result.radiation_power_w, powers[1], 1e-8);
    near(result.surface_residual_w, 0.0, 1e-8);
    assert!(
        surface < initial_temperature - 0.5,
        "convection must use the resolved surface, not the body's bulk temperature"
    );
    let delta = 1e-3;
    let a = reference(gas, initial_temperature + delta).1;
    let b = reference(gas, initial_temperature - delta).1;
    let tangent = ((a[0] + a[1]) - (b[0] + b[1]))
        / (energy(initial_temperature + delta) - energy(initial_temperature - delta));
    near(
        result.boundary.power_enthalpy_derivative_w_kg_j,
        tangent,
        1e-9,
    );

    let cell = CellId(0);
    let material_id = MaterialId(0);
    let mass = 0.1; // rho=1000 times full thickness 2*0.005 times face area .01.
    let initial = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cell,
                segment_id: 0,
                body_id: 0,
                material_id,
                mass_kg: mass,
                specific_h_j_kg: energy(initial_temperature),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.01, 0.0, 0.0],
                    volume_m3: 1e-4,
                    surface_area_m2: AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: Some(Reservoir {
                        area_m2: AREA,
                        conductance_w_k: 0.0,
                        temperature_k: 300.0,
                    }),
                },
            }),
        }],
        &[(material_id, layer.caloric_law())],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let observer = ConductionObserver::new(
        &initial,
        &[(material_id, conduction.conductivity_law())],
        &[],
        &[ConductiveReservoir {
            cell_id: cell,
            half_distance_m: geometry.half_distance_m,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap();
    let solver = ImplicitConductionChain::new(
        &initial,
        &observer,
        &[(material_id, &conduction as &dyn ImplicitConductionLaw)],
        || Ok(()),
    )
    .unwrap()
    .with_boundary_exchanges(&initial, &[(cell, &boundary)], &[], || Ok(()))
    .unwrap();
    let dt = 100.0;
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = solver
        .advance_to(
            &initial,
            dt,
            ImplicitConductionControls {
                maximum_dt_s: dt,
                maximum_halvings: 0,
                cell_tolerance_j: 1e-7,
                total_tolerance_j: 1e-7,
                ..Default::default()
            },
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let temperature = next.temperature_k(cell).unwrap();
    let (_, heat) = reference(gas, temperature);
    near(
        mass * (energy(temperature) - energy(initial_temperature)) + dt * (heat[0] + heat[1]),
        0.0,
        1e-7,
    );
    near(
        report.energy.external_heat_j,
        -dt * (heat[0] + heat[1]),
        1e-7,
    );
    near(report.energy.residual_j, 0.0, 1e-7);
    near(
        report.entropy.reservoir_entropy_input_j_k,
        -dt * (heat[0] / 300.0 + heat[1] / 305.0),
        1e-9,
    );
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(initial.time_s(), 0.0);
    near(
        initial.temperature_k(cell).unwrap(),
        initial_temperature,
        1e-12,
    );
    assert_eq!(
        next.cells()[&cell].boundary,
        initial.cells()[&cell].boundary
    );
    temperature
}

#[test]
fn original_gas_changes_accepted_conducting_wall_history_with_separate_bath_accounts() {
    let dry = gas("300 100000 0");
    let mixture = gas("mixture 300 100000 mole 2 He 0.1 N2 0.9");
    let dry_temperature = advance(&dry);
    let mixed_temperature = advance(&mixture);
    assert!(mixed_temperature < dry_temperature - 1e-3);
    assert_eq!(dry.reference().state().temperature, 300.0);
    assert_eq!(mixture.reference().state().temperature, 300.0);
}

#[test]
fn original_thermochemical_calorics_drive_the_same_implicit_wall_and_bath_accounts() {
    let make = |slope| {
        IsobaricGasLaw::try_from_model(
            thermochemical_source::model_transport_with_heat_capacity_slope(0.8, slope),
            thermochemical_source::point(300.0, 100000.0),
        )
        .unwrap()
    };
    let constant_cp = make(0.0);
    let variable_cp = make(0.002);
    let saved = variable_cp.reference().clone();
    // The original transport claims, density, pressure and composition are
    // identical. A higher actual gas Cp increases Pr and the plate's heat flux;
    // the native finite solid therefore cools farther in this same interval.
    let constant_temperature = advance(&constant_cp);
    let variable_temperature = advance(&variable_cp);
    assert!(variable_temperature < constant_temperature - 1e-4);
    assert!(variable_temperature > 300.0);
    assert_eq!(variable_cp.reference(), &saved);
}

#[test]
fn original_pressure_temperature_surfaces_drive_native_wall_storage_and_bath_accounts() {
    let low = super::surface::law(65536.0, false);
    let high = super::surface::law(131072.0, false);
    let saved = high.reference().clone();
    let low_temperature = advance(&low);
    let high_temperature = advance(&high);
    assert!(high_temperature < low_temperature - 1e-4);
    assert!(high_temperature > 300.0);
    assert_eq!(high.reference(), &saved);
}
