//! G1/G3/G4: actual circuit voltage, source-defined magnetic saturation and
//! finite enthalpy share the same endpoint. Synthetic data are not metal data.
#[path = "../../fs-em/tests/nonlinear_winding_circuit/common.rs"]
mod field;

use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell,
        material::MaterialCell,
        winding::{MagneticRegion, WindingThermalNetwork},
    },
    material_circuit::MaterialCircuitSolve,
};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint, SelectionPolicy};
use fs_material::{
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::SPECIFIC_HEAT_CAPACITY_DIMS,
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::Density;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const DEPTH: f64 = 0.25;
const VOLUME: f64 = 2.0 * DEPTH;
const CAPACITY: f64 =
    field::material::DENSITY_KG_M3 * field::material::SPECIFIC_HEAT_J_KG_K * VOLUME;
const DT: f64 = 0.25;
const DRIVE: f64 = 8.0;
const SERIES_OHM: f64 = 1.0;
const CELL: CellId = CellId(7);
const MATERIAL: MaterialId = MaterialId(9);
const GRAPH: [ResistorCell; 1] = [ResistorCell {
    branch: 2,
    cell: CELL,
}];

fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 1,
        maximum_circuit_trials: 1000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 1e-7,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-8,
            power_tolerance_w: 1e-8,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-13,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &field::material::point(200.0),
            upper: &field::material::point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}
fn thermal(curve: &ResolvedHeatCapacityCurve) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 0,
                body_id: 0,
                material_id: MATERIAL,
                mass_kg: field::material::DENSITY_KG_M3 * VOLUME,
                specific_h_j_kg: field::material::SPECIFIC_HEAT_J_KG_K * (300.0 - 200.0),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [2.0, 0.0, 0.0],
                    volume_m3: VOLUME,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MATERIAL, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn cells<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
) -> [MaterialCell<'a>; 1] {
    [MaterialCell {
        cell: CELL,
        material_id: MATERIAL,
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }]
}
fn routes() -> Vec<MagneticRegion> {
    (0..8)
        .map(|triangle| MagneticRegion::Cell {
            triangle,
            cell: CELL,
        })
        .collect()
}
fn network<'a>(
    cells: &'a [MaterialCell<'a>],
    regions: &'a [MagneticRegion],
) -> WindingThermalNetwork<'a> {
    WindingThermalNetwork {
        cells,
        links: &[],
        magnetic_regions: regions,
        graph_resistors: &GRAPH,
        volume_tolerance_m3: 1e-12,
    }
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-7 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

// Independent analytic reduction of two P1 square-star modes. No production
// field/constitutive routines or Newton routines are called by this oracle.
fn flux_density(current: f64, saturating: bool) -> f64 {
    if saturating && current > 1.0 {
        1.0 + (current - 1.0) / 4.0
    } else {
        current
    }
}
fn density_energy(b: f64, saturating: bool) -> f64 {
    if saturating && b > 1.0 {
        0.5 + (b - 1.0) + 2.0 * (b - 1.0).powi(2)
    } else {
        0.5 * b * b
    }
}
fn resistivity(t: f64) -> f64 {
    (1.0 + 0.005 * (t - 300.0)) / 18.0
}
fn winding_resistance(t: f64) -> f64 {
    72.0 * DEPTH * resistivity(t)
}
struct Reference {
    temperature: f64,
    current: f64,
    heat: f64,
    stored: f64,
    numerical: f64,
}
fn reference(i0: f64, t0: f64, saturating: bool) -> Reference {
    let b0 = flux_density(i0, saturating);
    let current = |t: f64| {
        let resistance = winding_resistance(t) + SERIES_OHM;
        let mut lo = 0.0;
        let mut hi = i0.max(DRIVE / resistance) + 1.0;
        for _ in 0..80 {
            let i = f64::midpoint(lo, hi);
            let residual =
                2.0 * DEPTH * (flux_density(i, saturating) - b0) + DT * resistance * i - DT * DRIVE;
            if residual > 0.0 {
                hi = i;
            } else {
                lo = i;
            }
        }
        f64::midpoint(lo, hi)
    };
    let mut lo = t0;
    let mut hi = t0 + 20.0;
    for _ in 0..80 {
        let t = f64::midpoint(lo, hi);
        let i = current(t);
        let heat = DT * (winding_resistance(t) + SERIES_OHM) * i * i;
        if CAPACITY * (t - t0) > heat {
            hi = t;
        } else {
            lo = t;
        }
    }
    let temperature = f64::midpoint(lo, hi);
    let current = current(temperature);
    let b = flux_density(current, saturating);
    let heat = DT * (winding_resistance(temperature) + SERIES_OHM) * current * current;
    let stored = VOLUME * density_energy(b, saturating);
    let numerical =
        current * 2.0 * DEPTH * (b - b0) - stored + VOLUME * density_energy(b0, saturating);
    Reference {
        temperature,
        current,
        heat,
        stored,
        numerical,
    }
}

#[test]
fn data_only_saturation_changes_actual_voltage_current_and_finite_thermal_feedback() {
    let mut outcomes = Vec::new();
    for saturating in [false, true] {
        let card = field::source(saturating, 18.0);
        let curve = caloric(&card);
        let mut thermal = thermal(&curve);
        let problem = field::problem(&card, DEPTH);
        let mut circuit = field::single(&problem, &card, SERIES_OHM, 0.0, DRIVE);
        let context = field::material::point(300.0);
        let cells = cells(&card, &curve, &context);
        let regions = routes();
        let layout = network(&cells, &regions);
        let mut cumulative_work = 0.0;
        let mut cumulative_numerical = 0.0;
        for step in 1..=2 {
            let old_temperature = thermal.temperature_k(CELL).unwrap();
            let old_current = circuit.currents_a()[0];
            let old_time = thermal.time_s();
            let expected = reference(old_current, old_temperature, saturating);
            let target = step as f64 * DT;
            let candidate = layout
                .step_nonlinear_candidate(
                    &thermal,
                    target,
                    &circuit,
                    &[DRIVE],
                    controls(),
                    field::controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let t = candidate.thermal().temperature_k(CELL).unwrap();
            let report = candidate.report();
            let electrical = &report.electrical;
            near(t, expected.temperature);
            near(candidate.circuit().currents_a()[0], expected.current);
            near(
                candidate.circuit().resistance_ohm()[0],
                winding_resistance(t),
            );
            near(electrical.stored_after_j, expected.stored);
            near(electrical.numerical_dissipation_j, expected.numerical);
            near(electrical.heat_j, expected.heat);
            near(report.deposited_joule_j, expected.heat);
            near(
                report.thermal.stored_energy_change_j,
                CAPACITY * (t - old_temperature),
            );
            near(report.thermal.stored_energy_change_j, expected.heat);
            near(electrical.supplied_work_j, DT * DRIVE * expected.current);
            near(report.non_electrical_heat_j, 0.0);
            assert!(report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
            assert!(report.maximum_cell_residual_j <= controls().cell_energy_tolerance_j);
            assert!(report.circuit_trials > 1);
            assert!(electrical.numerical_dissipation_j > 1e-4);
            // The raw native ledger includes internal winding-port work.
            assert!((electrical.native.supplied - electrical.supplied_work_j).abs() > 1e-3);
            assert_eq!(electrical.graph_resistors.len(), 1);
            near(
                electrical.graph_resistors[0].heat_j,
                DT * SERIES_OHM * expected.current.powi(2),
            );
            let mut coil_heat = 0.0;
            for loss in &electrical.cells {
                coil_heat += loss.heat_j;
                assert_eq!(
                    loss.source.conductor.material().card_identity(),
                    card.content_hash()
                );
                near(loss.source.conductor.resistivity_ohm_m(), resistivity(t));
                let source_t = loss
                    .source
                    .conductor
                    .material()
                    .query_point()
                    .iter()
                    .find(|(name, _)| name == "T")
                    .unwrap()
                    .1;
                near(source_t, t);
            }
            near(
                coil_heat,
                DT * winding_resistance(t) * expected.current.powi(2),
            );
            let frozen = circuit
                .step_to_candidate(
                    target,
                    circuit.field().conductors().to_vec(),
                    circuit.conductors().to_vec(),
                    &[DRIVE],
                    field::controls(),
                    &mut || false,
                )
                .unwrap();
            assert!(frozen.circuit().currents_a()[0] > expected.current + 1e-5);
            cumulative_work += electrical.supplied_work_j;
            cumulative_numerical += electrical.numerical_dissipation_j;
            near(
                candidate.circuit().field().energy_j()
                    + CAPACITY * (t - 300.0)
                    + cumulative_numerical,
                cumulative_work,
            );
            assert_eq!(circuit.time_s(), old_time);
            assert_eq!(thermal.time_s(), old_time);
            near(thermal.temperature_k(CELL).unwrap(), old_temperature);
            let (next_thermal, next_circuit, _) = candidate.into_parts();
            thermal = next_thermal;
            circuit = next_circuit;
        }
        outcomes.push((
            circuit.currents_a()[0],
            thermal.temperature_k(CELL).unwrap(),
            circuit.field().energy_j(),
        ));
    }
    assert!((outcomes[0].0 - outcomes[1].0).abs() > 0.2);
    assert!((outcomes[0].1 - outcomes[1].1).abs() > 0.1);
    assert!((outcomes[0].2 - outcomes[1].2).abs() > 0.05);
}

#[test]
fn incomplete_heat_routes_budget_and_late_cancel_preserve_all_three_histories() {
    let card = field::source(true, 18.0);
    let curve = caloric(&card);
    let thermal = thermal(&curve);
    let problem = field::problem(&card, DEPTH);
    let circuit = field::single(&problem, &card, SERIES_OHM, 0.0, DRIVE);
    let context = field::material::point(300.0);
    let cells = cells(&card, &curve, &context);
    let regions = routes();
    let layout = network(&cells, &regions);
    let missing_graph = WindingThermalNetwork {
        graph_resistors: &[],
        ..network(&cells, &regions)
    };
    assert!(
        missing_graph
            .step_nonlinear_candidate(
                &thermal,
                DT,
                &circuit,
                &[DRIVE],
                controls(),
                field::controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    assert!(
        network(&cells, &regions[..7])
            .step_nonlinear_candidate(
                &thermal,
                DT,
                &circuit,
                &[DRIVE],
                controls(),
                field::controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    let mut limited = controls();
    limited.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_nonlinear_candidate(
            &thermal,
            DT,
            &circuit,
            &[DRIVE],
            limited,
            field::controls(),
            &[0.0],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut checkpoints = 0;
    let success = layout
        .step_nonlinear_candidate(
            &thermal,
            DT,
            &circuit,
            &[DRIVE],
            controls(),
            field::controls(),
            &[0.0],
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut visits = 0;
    assert!(matches!(
        layout.step_nonlinear_candidate(
            &thermal,
            DT,
            &circuit,
            &[DRIVE],
            controls(),
            field::controls(),
            &[0.0],
            || {
                visits += 1;
                visits == checkpoints
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(circuit.currents_a(), &[0.0]);
    assert_eq!(circuit.time_s(), 0.0);
    assert_eq!(circuit.field().time_s(), 0.0);
    near(circuit.field().energy_j(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    near(thermal.temperature_k(CELL).unwrap(), 300.0);
    assert!(
        layout
            .step_nonlinear_candidate(
                &thermal,
                2.0 * DT,
                success.circuit(),
                &[DRIVE],
                controls(),
                field::controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    let retry = layout
        .step_nonlinear_candidate(
            &thermal,
            DT,
            &circuit,
            &[DRIVE],
            controls(),
            field::controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(retry.circuit().state(), success.circuit().state());
    assert_eq!(
        retry.circuit().field().potential_wb_m(),
        success.circuit().field().potential_wb_m()
    );
    assert_eq!(
        retry.thermal().cells()[&CELL].specific_h_j_kg.to_bits(),
        success.thermal().cells()[&CELL].specific_h_j_kg.to_bits()
    );
}
