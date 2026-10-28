//! Manufactured first-principles feedback fixtures, not identified rosin data.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::*;
use fs_matdb::SelectionPolicy;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{
    ResolvedShearLayer, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_qty::Temperature;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::collections::BTreeMap;
use std::convert::Infallible;

const AREA: f64 = 2e-6;
const THICKNESS: f64 = 1e-6;
const DENSITY: f64 = 1000.0;
const CP: f64 = 1000.0;
const FILM: CellId = CellId(10);

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("omega", 0.0)
        .unwrap()
}

fn coating(
    variable_loss: bool,
    variable_modulus: bool,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    coating_law(variable_loss, variable_modulus.then(||
        vec![(280.0, 8.4e6), (340.0, 7.2e6)]), None)
}

fn coating_law(
    variable_loss: bool,
    modulus_knots: Option<Vec<(f64, f64)>>,
    elastic_yield_pa: Option<f64>,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    coating_law_with_capacity(variable_loss, modulus_knots, elastic_yield_pa, CP)
}

fn coating_law_with_capacity(
    variable_loss: bool,
    modulus_knots: Option<Vec<(f64, f64)>>,
    elastic_yield_pa: Option<f64>,
    heat_capacity: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot, varies) in [
        ("density", Density::DIMS, DENSITY, DENSITY, false),
        (
            "shear_modulus",
            Pressure::DIMS,
            8e6,
            8e6,
            modulus_knots.is_some(),
        ),
        (
            "shear_yield_stress",
            Pressure::DIMS,
            elastic_yield_pa.unwrap_or(if variable_loss { 120.0 } else { 100.0 }),
            elastic_yield_pa.unwrap_or(if variable_loss { 60.0 } else { 100.0 }),
            variable_loss && elastic_yield_pa.is_none(),
        ),
        (
            "plastic_shear_viscosity",
            DynViscosity::DIMS,
            if variable_loss { 7.0 } else { 6.0 },
            if variable_loss { 4.0 } else { 6.0 },
            variable_loss,
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            heat_capacity,
            heat_capacity,
            false,
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            0.2,
            0.2,
            false,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if varies {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: if name == "shear_modulus" {
                            modulus_knots.clone().expect("variable modulus has its knots")
                        } else { vec![(280.0, cold), (340.0, hot)] },
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value: cold, dims }
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 340.0)
                    .with("omega", 0.0, 1e6),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if varies {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured thermal shear-layer fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured-film".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(340.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * heat_capacity,
        1e-8,
    )
    .unwrap();
    (card, caloric)
}

fn layer(card: &MaterialCard) -> ResolvedShearLayer {
    let requirements = [
        ("density", Density::DIMS),
        ("shear_modulus", Pressure::DIMS),
        ("shear_yield_stress", Pressure::DIMS),
        ("plastic_shear_viscosity", DynViscosity::DIMS),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    })
    .collect::<Vec<_>>();
    let source = resolve_material_state_point(
        card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedShearLayer::new(source, AREA, THICKNESS, 0.05).unwrap()
}

fn instrument(cx: &Cx<'_>, gate: &CancelGate, card: &MaterialCard) -> MultiBowStringBank {
    reacting_surface(joint(cx, gate, 7800.0), gate)
        .with_bow_shear_layers(vec![(0, layer(card))], gate)
        .unwrap()
}

fn nonlinear_instrument(cx: &Cx<'_>, gate: &CancelGate, card: &MaterialCard) -> MultiBowStringBank {
    use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;
    let bank = base(cx, gate, 7800.0).with_geometric_stretching(0.2,
        NonlinearModalConfig { max_iterations: 64, maximum_generalized_force: 1e5,
            force_absolute_tolerance: 1e-11, force_relative_tolerance: 1e-10 }, gate)
        .unwrap().with_additional_bows(vec![extra()], joint_limits(), gate).unwrap();
    reacting_surface(bank, gate).with_bow_shear_layers(vec![(0, layer(card))], gate).unwrap()
}

fn thermal<'a>(
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    neighbors: bool,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for i in 0..if neighbors { 3 } else { 1 } {
        let depth = if i == 0 { THICKNESS } else { 10.0 * THICKNESS };
        let volume = AREA * depth;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(10 + i),
                segment_id: i,
                body_id: i,
                material_id: MaterialId(0),
                mass_kg: DENSITY * volume,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64 * 1e-5],
                    end_m: [0.002, 0.0, i as f64 * 1e-5],
                    volume_m3: volume,
                    surface_area_m2: 2.0 * AREA + 0.006 * depth,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    if neighbors {
        for i in 1..3 {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: i,
                    a: FILM,
                    b: CellId(10 + i),
                    area_m2: AREA,
                    // Actual uniform source conductivity and half-cell distances.
                    conductance_w_k: series_conductance_w_k(
                        AREA,
                        [0.5 * THICKNESS, 5.0 * THICKNESS],
                        [0.2, 0.2],
                    )
                    .unwrap(),
                }),
            });
        }
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 8,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1e-5,
            energy_tolerance_j: 1e-14,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn binding<'a>(
    card: &'a MaterialCard,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
) -> [LayerThermalCell<'a>; 1] {
    [LayerThermalCell {
        port: 0,
        cell: FILM,
        material: card,
        caloric,
        temperature_axis: "T",
    }]
}
fn controls(cells: usize) -> LayerThermalSolve {
    LayerThermalSolve {
        maximum_cells: 3,
        maximum_mechanical_trials: 256,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-10,
        cell_energy_tolerance_j: 1e-13,
        thermal_energy_tolerance_j: 1e-13,
        total_energy_tolerance_j: 2e-12,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-15,
            relative_tolerance: 1e-12,
            linear_restart: cells,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn other(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, Infallible> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}
fn total_heat(thermal: &DepositionSimulation<'_>) -> f64 {
    thermal
        .cells()
        .values()
        .map(|c| c.mass_kg * c.specific_h_j_kg)
        .sum()
}

#[test]
fn layer_heat_requeries_same_material_and_changes_actual_instrument_pressure() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(true, false);
        let (frozen_card, frozen_caloric) = coating(false, false);
        let mut bank = instrument(cx, gate, &card);
        let mut frozen = instrument(cx, gate, &frozen_card);
        let mut heat = thermal(caloric.curve(), false);
        let mut frozen_heat = thermal(frozen_caloric.curve(), false);
        let bound = binding(&card, caloric.curve());
        let frozen_bound = binding(&frozen_card, frozen_caloric.curve());
        let original = bank.bow_shear_layer(0).unwrap().clone();
        let forces = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.001, 0.1), drive(0.0, 0.0)];
        let (mut plastic, mut pressure_difference) = (0.0, 0.0);
        for _ in 0..512 {
            let f = bank
                .step_with_layer_thermal_network(
                    &mut heat,
                    &bound,
                    controls(1),
                    &[0.0],
                    &forces,
                    &drives,
                    gate,
                    other,
                )
                .unwrap();
            let b = frozen
                .step_with_layer_thermal_network(
                    &mut frozen_heat,
                    &frozen_bound,
                    controls(1),
                    &[0.0],
                    &forces,
                    &drives,
                    gate,
                    other,
                )
                .unwrap();
            assert_eq!(f.mechanics.ports[0].dissipation_j, f.routed_plastic_heat_j);
            assert!(f.routed_plastic_heat_j >= 0.0);
            assert!((f.deposited_heat_j - f.routed_plastic_heat_j).abs() < 1e-15);
            assert_eq!(f.external_heat_j, 0.0);
            assert!(f.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
            assert!(f.mechanics.energy_residual_j.abs() <= f.mechanics.energy_tolerance_j);
            assert_eq!(
                bank.bow_shear_layer(0).unwrap().temperature_kelvin(),
                heat.temperature_k(FILM).unwrap()
            );
            plastic += f.routed_plastic_heat_j;
            pressure_difference +=
                (f.mechanics.observer_pressure_pa - b.mechanics.observer_pressure_pa).abs();
        }
        assert!(plastic > 0.0 && heat.temperature_k(FILM).unwrap() > 300.0);
        assert!((total_heat(&heat) - plastic).abs() < 1e-11);
        let warm = bank.bow_shear_layer(0).unwrap();
        assert_eq!(warm.stiffness_n_m(), original.stiffness_n_m());
        assert!(warm.plastic_resistance_n_s_m() < original.plastic_resistance_n_s_m());
        for (a, b) in warm
            .resolved()
            .properties()
            .iter()
            .zip(original.resolved().properties())
        {
            assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        }
        assert!(
            pressure_difference > 1e-14,
            "temperature must reach the original pressure observer"
        );
        assert_eq!(bank.system().samples_rendered(), 512);
        assert!(
            bank.system()
                .material_shear_state(0)
                .unwrap()
                .plastic_slip_m()[0]
                != 0.0
        );
    });
}

#[test]
fn nonlinear_string_and_film_internal_energy_close_while_heat_capacity_changes_pressure() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating_law_with_capacity(true,
            Some(vec![(280.0, 8.4e6), (340.0, 7.2e6)]), None, 10.0);
        let (slow_card, slow_caloric) = coating_law_with_capacity(true,
            Some(vec![(280.0, 8.4e6), (340.0, 7.2e6)]), None, CP);
        let mut bank = nonlinear_instrument(cx, gate, &card);
        let mut slow = nonlinear_instrument(cx, gate, &slow_card);
        let (mut heat, mut slow_heat) =
            (thermal(caloric.curve(), false), thermal(slow_caloric.curve(), false));
        let bound = binding(&card, caloric.curve());
        let slow_bound = binding(&slow_card, slow_caloric.curve());
        let original = bank.bow_shear_layer(0).unwrap().clone();
        let forces = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.001, 0.1), drive(0.0, 0.0)];
        let (mut work, mut external_loss, mut difference, mut strain) = (0.0, 0.0, 0.0, 0.0_f64);
        for _ in 0..256 {
            let f = bank.step_with_layer_thermal_network(&mut heat, &bound,
                controls(1), &[0.0], &forces, &drives, gate, other).unwrap();
            let reference = slow.step_with_layer_thermal_network(&mut slow_heat, &slow_bound,
                controls(1), &[0.0], &forces, &drives, gate, other).unwrap();
            let nonlinear = f.mechanics.nonlinear.unwrap();
            assert!(nonlinear.reaction_residual <= nonlinear.reaction_tolerance);
            assert!(f.mechanics.energy_residual_j.abs() <= f.mechanics.energy_tolerance_j);
            assert!(f.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
            assert!((f.routed_reversible_heat_j + f.coefficient_free_energy_change_j
                + f.mechanical_entropy_energy_change_j).abs() < 1e-15);
            work += f.mechanics.external_work_j + f.mechanics.drive_work_j;
            external_loss += f.mechanics.network_dissipation_j
                + f.mechanics.friction_dissipation_j - f.routed_plastic_heat_j;
            difference += (f.mechanics.observer_pressure_pa - reference.mechanics.observer_pressure_pa).abs();
            strain = strain.max(nonlinear.nonlinear_energy_j);

            // Independent internal-energy expression for affine G(T), plus
            // the original complete string/body/cavity/nonlinear storage.
            let e = bank.system().material_shear_state(0).unwrap().elastic_slip_m();
            let film_internal = 0.5 * AREA / THICKNESS * 14e6 * (e[0]*e[0] + e[1]*e[1]);
            let structural = bank.system().total_energy_j().unwrap() - f.mechanics.shear_layer_energy_j;
            assert!((structural + film_internal + total_heat(&heat) + external_loss - work).abs() < 2e-9);
        }
        assert!(strain > 0.0 && difference > 1e-14);
        assert!(heat.temperature_k(FILM).unwrap() > slow_heat.temperature_k(FILM).unwrap());
        assert!(slow_heat.temperature_k(FILM).unwrap() > 300.0);
        assert!(bank.bow_shear_layer(0).unwrap().stiffness_n_m() < original.stiffness_n_m());
        for (current, initial) in bank.bow_shear_layer(0).unwrap().resolved().properties().iter()
            .zip(original.resolved().properties()) {
            assert_eq!(current.answer().receipt.selected, initial.answer().receipt.selected);
        }
        assert_eq!(bank.system().samples_rendered(), 256);
        assert_eq!(heat.time_s(), 256.0 / 48_000.0);
    });
}

#[test]
fn original_network_conducts_film_heat_without_repeating_mechanical_history() {
  for nonlinear in [false, true] {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let make = || if nonlinear { nonlinear_instrument(cx, gate, &card) }
            else { instrument(cx, gate, &card) };
        let mut heated = make();
        let mut ordinary = make();
        let mut heat = thermal(caloric.curve(), true);
        let bound = binding(&card, caloric.curve());
        let forces = vec![0.0; heated.system().mode_count()];
        let drives = [drive(0.001, 0.1), drive(0.0, 0.0)];
        let mut plastic = 0.0;
        for _ in 0..128 {
            let f = heated
                .step_with_layer_thermal_network(
                    &mut heat,
                    &bound,
                    controls(3),
                    &[0.0; 3],
                    &forces,
                    &drives,
                    gate,
                    other,
                )
                .unwrap();
            let reference = ordinary.step(&forces, &drives).unwrap();
            assert_eq!(f.mechanics.nonlinear.is_some(), nonlinear);
            // Extra endpoint/Jacobian trials do not age mechanical/strain state.
            assert_eq!(
                f.mechanics.observer_pressure_pa,
                reference.observer_pressure_pa
            );
            assert_eq!(states(&heated), states(&ordinary));
            assert_eq!(
                heated
                    .system()
                    .material_shear_state(0)
                    .unwrap()
                    .elastic_slip_m(),
                ordinary
                    .system()
                    .material_shear_state(0)
                    .unwrap()
                    .elastic_slip_m()
            );
            assert_eq!(
                heated
                    .system()
                    .material_shear_state(0)
                    .unwrap()
                    .plastic_slip_m(),
                ordinary
                    .system()
                    .material_shear_state(0)
                    .unwrap()
                    .plastic_slip_m()
            );
            assert!(f.total_energy_residual_j.abs() <= controls(3).total_energy_tolerance_j);
            plastic += f.routed_plastic_heat_j;
        }
        assert!((total_heat(&heat) - plastic).abs() < 1e-11);
        assert!(heat.temperature_k(CellId(11)).unwrap() > 300.0);
        assert_eq!(
            heat.temperature_k(CellId(11)).unwrap(),
            heat.temperature_k(CellId(12)).unwrap()
        );
        assert!(heat.temperature_k(FILM).unwrap() > heat.temperature_k(CellId(11)).unwrap());
        assert_eq!(
            heated.system().samples_rendered(),
            ordinary.system().samples_rendered()
        );
    });
  }
}

#[test]
fn thermal_and_contact_failures_preserve_all_accepted_instrument_and_heat_state() {
  for nonlinear in [false, true] {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(true, true);
        let mut bank = if nonlinear { nonlinear_instrument(cx, gate, &card) }
            else { instrument(cx, gate, &card) };
        let mut heat = thermal(caloric.curve(), false);
        let bound = binding(&card, caloric.curve());
        let forces = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.001, 0.1), drive(0.0, 0.0)];
        for _ in 0..8 {
            bank.step_with_layer_thermal_network(
                &mut heat,
                &bound,
                controls(1),
                &[0.0],
                &forces,
                &drives,
                gate,
                other,
            )
            .unwrap();
        }
        let old_states = states(&bank);
        let old_source = bank.bow_shear_layer(0).unwrap().clone();
        let old_history = bank.system().material_shear_state(0).unwrap().clone();
        let old_frame = bank.system().last_frame().unwrap().clone();
        assert_eq!(old_frame.nonlinear.is_some(), nonlinear);
        let old_cells = heat.cells().clone();
        let old_time = heat.time_s();
        let assert_unchanged = |bank: &MultiBowStringBank, heat: &DepositionSimulation<'_>| {
            assert_eq!(states(bank), old_states);
            assert_eq!(bank.bow_shear_layer(0).unwrap(), &old_source);
            assert_eq!(bank.system().material_shear_state(0).unwrap(), &old_history);
            assert_eq!(bank.system().last_frame().unwrap(), &old_frame);
            assert_eq!(bank.system().samples_rendered(), 8);
            assert_eq!(heat.time_s(), old_time);
            assert_eq!(heat.cells(), &old_cells);
        };
        assert!(
            bank.step_with_layer_thermal_network(
                &mut heat,
                &bound,
                controls(1),
                &[0.0],
                &forces,
                &[drive(0.001, 0.0), drives[1]],
                gate,
                other
            )
            .is_err()
        );
        assert_unchanged(&bank, &heat);
        assert!(
            bank.step_with_layer_thermal_network(
                &mut heat,
                &bound,
                controls(1),
                &[0.0],
                &forces,
                &drives,
                gate,
                |_| Err::<StepCoefficients<()>, _>("transport source refused")
            )
            .is_err()
        );
        assert_unchanged(&bank, &heat);
        // Tiny physical work must not vanish beside an enormous other source.
        assert!(
            bank.step_with_layer_thermal_network(
                &mut heat,
                &bound,
                controls(1),
                &[0.0],
                &forces,
                &drives,
                gate,
                |_| Ok::<_, Infallible>(StepCoefficients {
                    sources_w: BTreeMap::from([(FILM, 1e100)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            )
            .is_err()
        );
        assert_unchanged(&bank, &heat);
        let mut exhausted = controls(1);
        exhausted.maximum_mechanical_trials = 1;
        assert!(
            bank.step_with_layer_thermal_network(
                &mut heat,
                &bound,
                exhausted,
                &[0.0],
                &forces,
                &drives,
                gate,
                other
            )
            .is_err()
        );
        assert_unchanged(&bank, &heat);
        assert!(
            bank.step_with_layer_thermal_network(
                &mut heat,
                &bound,
                controls(1),
                &[0.0],
                &forces,
                &drives,
                gate,
                |endpoint| {
                    gate.request();
                    other(endpoint)
                }
            )
            .is_err()
        );
        assert_unchanged(&bank, &heat);
    });
  }
}

#[test]
fn thermal_continuation_refuses_modulus_kinks_and_inconsistent_caloric_density() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating_law(false,
            Some(vec![(280.0, 8.4e6), (300.0, 8e6), (340.0, 7.6e6)]), None);
        let mut bank = instrument(cx, gate, &card);
        let mut heat = thermal(caloric.curve(), false);
        let forces = vec![0.0; bank.system().mode_count()];
        assert!(
            bank.step_with_layer_thermal_network(
                &mut heat,
                &binding(&card, caloric.curve()),
                controls(1),
                &[0.0],
                &forces,
                &[drive(0.001, 0.1), drive(0.0, 0.0)],
                gate,
                other
            )
            .is_err()
        );
        assert_eq!(bank.system().samples_rendered(), 0);
        assert_eq!(heat.time_s(), 0.0);

        let (card, _) = coating(false, false);
        let wrong = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            card.content_hash(),
            SolidLiquidPhase::Solid,
            vec![
                EnthalpyPhaseKnot {
                    temperature_k: 280.0,
                    specific_enthalpy_j_kg: -20.0 * CP,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: DENSITY,
                },
                EnthalpyPhaseKnot {
                    temperature_k: 300.0,
                    specific_enthalpy_j_kg: 0.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: DENSITY,
                },
                EnthalpyPhaseKnot {
                    temperature_k: 340.0,
                    specific_enthalpy_j_kg: 40.0 * CP,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1.1 * DENSITY,
                },
            ],
        )
        .unwrap();
        let mut bank = instrument(cx, gate, &card);
        let mut heat = thermal(&wrong, false);
        assert!(
            bank.step_with_layer_thermal_network(
                &mut heat,
                &binding(&card, &wrong),
                controls(1),
                &[0.0],
                &forces,
                &[drive(0.001, 0.1), drive(0.0, 0.0)],
                gate,
                other
            )
            .is_err()
        );
        assert_eq!(bank.system().samples_rendered(), 0);
        assert_eq!(heat.time_s(), 0.0);
    });
}

#[test]
fn affine_layer_softening_routes_reversible_exchange_through_the_same_instrument() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, true);
        let (constant_card, constant_caloric) = coating(false, false);
        let mut bank = instrument(cx, gate, &card);
        let mut reference = instrument(cx, gate, &constant_card);
        let mut heat = thermal(caloric.curve(), false);
        let mut reference_heat = thermal(constant_caloric.curve(), false);
        let bound = binding(&card, caloric.curve());
        let reference_bound = binding(&constant_card, constant_caloric.curve());
        let original = bank.bow_shear_layer(0).unwrap().clone();
        assert_eq!(original.stiffness_n_m(), reference.bow_shear_layer(0).unwrap().stiffness_n_m());
        let forces = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.001, 0.1), drive(0.0, 0.0)];
        let (mut plastic, mut reversible, mut stored_internal,
            mut pressure_difference) = (0.0, 0.0, 0.0, 0.0);
        for _ in 0..256 {
            let f = bank.step_with_layer_thermal_network(&mut heat, &bound,
                controls(1), &[0.0], &forces, &drives, gate, other).unwrap();
            let fixed = reference.step_with_layer_thermal_network(&mut reference_heat,
                &reference_bound, controls(1), &[0.0], &forces, &drives, gate, other).unwrap();
            plastic += f.routed_plastic_heat_j;
            reversible += f.routed_reversible_heat_j;
            stored_internal += f.layer_internal_energy_change_j;
            pressure_difference += (f.mechanics.observer_pressure_pa
                - fixed.mechanics.observer_pressure_pa).abs();
            assert!(f.routed_plastic_heat_j >= 0.0);
            assert!((f.routed_reversible_heat_j + f.coefficient_free_energy_change_j
                +f.mechanical_entropy_energy_change_j).abs() < 1e-15);
            assert!(f.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
            let layer = bank.bow_shear_layer(0).unwrap();
            let t = heat.temperature_k(FILM).unwrap();
            assert_eq!(layer.temperature_kelvin(), t);
            let expected = AREA/THICKNESS*(8.4e6-2e4*(t-280.0));
            assert!((layer.stiffness_n_m()-expected).abs() <= 1e-14*expected);
            for (a, b) in layer.resolved().properties().iter()
                .zip(original.resolved().properties()) {
                assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
            }
        }
        let e = bank.system().material_shear_state(0).unwrap().elastic_slip_m();
        // k(T)=(A/h)*(14e6-2e4*T), so the mechanical INTERNAL energy is
        // (A/h)*14e6*|e|^2/2, independently of the endpoint T.
        let expected_internal = 0.5*AREA/THICKNESS*14e6*(e[0]*e[0]+e[1]*e[1]);
        assert!((stored_internal-expected_internal).abs() < 1e-14);
        assert!((total_heat(&heat)-plastic-reversible).abs() < 1e-11);
        assert!(plastic > 0.0 && reversible != 0.0);
        assert!(bank.bow_shear_layer(0).unwrap().stiffness_n_m() < original.stiffness_n_m());
        assert!(pressure_difference > 1e-14);
    });
}

#[test]
fn reversible_shear_cooling_is_not_clipped_or_relabelled_as_plastic_dissipation() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating_law(false,
            Some(vec![(280.0, 8.4e6), (340.0, 7.2e6)]), Some(1e8));
        let mut bank = instrument(cx, gate, &card);
        let mut heat = thermal(caloric.curve(), false);
        let forces = vec![0.0; bank.system().mode_count()];
        let f = bank.step_with_layer_thermal_network(&mut heat,
            &binding(&card, caloric.curve()), controls(1), &[0.0], &forces,
            &[drive(0.001, 0.1), drive(0.0, 0.0)], gate, other).unwrap();
        let e = bank.system().material_shear_state(0).unwrap().elastic_slip_m();
        let e2 = e[0]*e[0]+e[1]*e[1];
        let stiffness_slope = -2e4*AREA/THICKNESS;
        let capacity = DENSITY*AREA*THICKNESS*CP;
        let expected_t = 300.0/(1.0-stiffness_slope*e2/(2.0*capacity));
        assert_eq!(f.routed_plastic_heat_j, 0.0);
        assert!(f.routed_reversible_heat_j < 0.0 && f.deposited_heat_j < 0.0);
        assert!(heat.temperature_k(FILM).unwrap() < 300.0);
        assert!((heat.temperature_k(FILM).unwrap()-expected_t).abs() < 1e-8);
        assert!((total_heat(&heat)-f.routed_reversible_heat_j).abs() < 1e-13);
        assert!(f.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
    });
}

// A geometry-derived, fixed-coefficient hair drive on the original instrument.
// The receiving caloric cells below are explicit manufactured lumped bodies;
// this does not identify a measured thermal law for horsehair or rosin.
fn compliant_instrument(cx: &Cx<'_>, gate: &CancelGate, distributed_damping: f64)
    -> (MultiBowStringBank, f64, f64) {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::compliant_ribbon::{
        CompliantBowRibbon, HairRibbonProperties,
    };
    let hair = HairRibbonProperties {
        youngs_modulus_pa: 7e9, hair_radius_m: 0.095e-3,
        hairs_in_contact: 50.0, free_length_m: 0.65, excitation_fraction: 0.5,
        distributed_damping_kg_s_m: distributed_damping,
        reference_width_m: 0.01, density_kg_m3: 1100.0,
    };
    let ribbon = CompliantBowRibbon::new(options().bow, 0.5, 0.01, 1, hair).unwrap();
    let bank = base(cx, gate, 7800.0)
        .with_compliant_bow_ribbon(0.01, 1, hair, joint_limits(), gate).unwrap();
    (bank, ribbon.lumped_tangential_stiffness_n_m(), ribbon.lumped_tangential_damping_kg_s())
}

#[test]
fn series_hair_heat_is_its_actual_dashpot_loss_and_preserves_original_motion() {
    with_cx(|cx, gate| {
        for damping in [0.0, 1000.0] {
            let (_, caloric) = coating(false, false);
            let (mut bank, k, c) = compliant_instrument(cx, gate, damping);
            let (mut ordinary, _, _) = compliant_instrument(cx, gate, damping);
            let mut heat = thermal(caloric.curve(), false);
            let routes = [SeriesDriveHeatDestination { port: 0, cell: FILM }];
            let forces = vec![0.0; bank.system().mode_count()];
            let drives = [drive(0.001, 0.1)];
            let dt = bank.system().sample_period_s();
            let (mut eta, mut generated) = (0.0, 0.0);
            for _ in 0..128 {
                let f = bank.step_with_contact_thermal_network(&mut heat, &[], &routes,
                    controls(1), &[0.0], &forces, &drives, gate, other).unwrap();
                let reference = ordinary.step(&forces, &drives).unwrap();
                assert_eq!(&f.mechanics, reference);
                assert_eq!(states(&bank), states(&ordinary));
                let traction = f.mechanics.ports[0].traction_n;
                let next = eta + (traction-k*eta)/(0.5*k+c/dt);
                // Independently subtract elastic storage from drive-end work.
                let work = traction*(next-eta);
                let stored = 0.5*k*(next*next-eta*eta);
                let expected = work-stored;
                assert!((f.routed_series_drive_heat_j-expected).abs()
                    <= 1e-12*work.abs().max(stored.abs())+1e-18);
                assert_eq!(f.routed_series_drive_heat_j,
                    f.mechanics.ports[0].series_drive_dissipation_j);
                assert_eq!(f.routed_plastic_heat_j, 0.0);
                assert_eq!(f.routed_reversible_heat_j, 0.0);
                assert!((f.deposited_heat_j-f.routed_series_drive_heat_j).abs() < 1e-15);
                assert!(f.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
                generated += f.routed_series_drive_heat_j;
                eta = next;
            }
            assert!((total_heat(&heat)-generated).abs() < 1e-11);
            if damping == 0.0 {
                assert_eq!(generated, 0.0);
                assert_eq!(heat.temperature_k(FILM).unwrap(), 300.0);
            } else {
                assert!(generated > 0.0 && heat.temperature_k(FILM).unwrap() > 300.0);
            }
        }
    });
}

#[test]
fn hair_heat_conducts_to_the_original_material_film_and_changes_instrument_pressure() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(true, true);
        let make = || compliant_instrument(cx, gate, 1000.0).0
            .with_bow_shear_layers(vec![(0, layer(&card))], gate).unwrap();
        let (mut bank, mut unrouted) = (make(), make());
        let (mut heat, mut reference_heat) =
            (thermal(caloric.curve(), true), thermal(caloric.curve(), true));
        let bound = binding(&card, caloric.curve());
        let routes = [SeriesDriveHeatDestination { port: 0, cell: CellId(11) }];
        let original = bank.bow_shear_layer(0).unwrap().clone();
        let forces = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.001, 0.1)];
        let (mut hair, mut plastic, mut reversible, mut pressure_difference) =
            (0.0, 0.0, 0.0, 0.0);
        for _ in 0..256 {
            let f = bank.step_with_contact_thermal_network(&mut heat, &bound, &routes,
                controls(3), &[0.0; 3], &forces, &drives, gate, other).unwrap();
            let reference = unrouted.step_with_layer_thermal_network(&mut reference_heat,
                &bound, controls(3), &[0.0; 3], &forces, &drives, gate, other).unwrap();
            hair += f.routed_series_drive_heat_j;
            plastic += f.routed_plastic_heat_j;
            reversible += f.routed_reversible_heat_j;
            pressure_difference += (f.mechanics.observer_pressure_pa
                -reference.mechanics.observer_pressure_pa).abs();
            assert!(f.routed_series_drive_heat_j >= 0.0 && f.routed_plastic_heat_j >= 0.0);
            assert!((f.deposited_heat_j-f.routed_plastic_heat_j
                -f.routed_reversible_heat_j-f.routed_series_drive_heat_j).abs() < 1e-15);
            assert_eq!(reference.routed_series_drive_heat_j, 0.0);
            assert!(f.total_energy_residual_j.abs() <= controls(3).total_energy_tolerance_j);
            assert_eq!(bank.bow_shear_layer(0).unwrap().temperature_kelvin(),
                heat.temperature_k(FILM).unwrap());
        }
        assert!(hair > 0.0 && plastic > 0.0);
        assert!((total_heat(&heat)-hair-plastic-reversible).abs() < 1e-11);
        assert!(heat.temperature_k(CellId(11)).unwrap() > heat.temperature_k(CellId(12)).unwrap());
        assert!((heat.temperature_k(FILM).unwrap()-reference_heat.temperature_k(FILM).unwrap()).abs() > 1e-10);
        assert!(pressure_difference > 1e-14);
        let current = bank.bow_shear_layer(0).unwrap();
        for (now, old) in current.resolved().properties().iter()
            .zip(original.resolved().properties()) {
            assert_eq!(now.answer().receipt.selected, old.answer().receipt.selected);
        }
        assert_eq!(bank.system().samples_rendered(), 256);
    });
}

#[test]
fn drive_heat_route_and_late_source_refusals_preserve_stretch_heat_and_pressure_history() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(true, true);
        let make = || compliant_instrument(cx, gate, 1000.0).0
            .with_bow_shear_layers(vec![(0, layer(&card))], gate).unwrap();
        let (mut bank, mut reference) = (make(), make());
        let (mut heat, mut reference_heat) =
            (thermal(caloric.curve(), true), thermal(caloric.curve(), true));
        let bound = binding(&card, caloric.curve());
        let route = SeriesDriveHeatDestination { port: 0, cell: CellId(11) };
        let forces = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.001, 0.1)];
        bank.step_with_contact_thermal_network(&mut heat, &bound, &[route],
            controls(3), &[0.0; 3], &forces, &drives, gate, other).unwrap();
        reference.step_with_contact_thermal_network(&mut reference_heat, &bound, &[route],
            controls(3), &[0.0; 3], &forces, &drives, gate, other).unwrap();
        let saved = (states(&bank), bank.system().last_frame().cloned(), heat.cells().clone(),
            bank.system().total_energy_j().unwrap(), heat.time_s(), bank.bow_shear_layer(0).cloned());
        for routes in [vec![route, route],
            vec![SeriesDriveHeatDestination { port: 1, cell: CellId(11) }],
            vec![SeriesDriveHeatDestination { port: 0, cell: CellId(999) }]] {
            assert!(bank.step_with_contact_thermal_network(&mut heat, &bound, &routes,
                controls(3), &[0.0; 3], &forces, &drives, gate, other).is_err());
        }
        assert!(matches!(bank.step_with_contact_thermal_network(&mut heat, &bound, &[route],
            controls(3), &[0.0; 3], &forces, &drives, gate,
            |_| -> Result<StepCoefficients<()>, &'static str> { Err("late source refusal") }),
            Err(LayerThermalError::Source("late source refusal"))));
        assert_eq!(states(&bank), saved.0);
        assert_eq!(bank.system().last_frame(), saved.1.as_ref());
        assert_eq!(heat.cells(), &saved.2);
        assert_eq!(bank.system().total_energy_j().unwrap(), saved.3);
        assert_eq!(heat.time_s(), saved.4);
        assert_eq!(bank.bow_shear_layer(0), saved.5.as_ref());
        let retry = bank.step_with_contact_thermal_network(&mut heat, &bound, &[route],
            controls(3), &[0.0; 3], &forces, &drives, gate, other).unwrap();
        let next = reference.step_with_contact_thermal_network(&mut reference_heat, &bound, &[route],
            controls(3), &[0.0; 3], &forces, &drives, gate, other).unwrap();
        assert_eq!(retry.mechanics, next.mechanics);
        assert_eq!(retry.routed_series_drive_heat_j, next.routed_series_drive_heat_j);
        assert_eq!(states(&bank), states(&reference));
        assert_eq!(heat.cells(), reference_heat.cells());
    });
}

#[path = "shear_thermal/performance.rs"]
mod performance;

#[path = "shear_thermal/dry.rs"]
mod dry;
