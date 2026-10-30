//! G1/G3: the actual nonlinear field, source conductivity and finite enthalpy
//! receiver share one endpoint. Synthetic coefficients are not metal data.
#[path = "../../fs-em/tests/nonlinear_eddy/common.rs"]
mod field;

#[path = "nonlinear_eddy_enthalpy/thermal.rs"]
mod thermal_coupling;

use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve,
        magnetic::{MagneticRegion, MagneticThermalError, MagneticThermalNetwork},
        material::MaterialCell,
    },
    material_circuit::MaterialCircuitSolve,
};
use fs_em::planar::{
    eddy::nonlinear::PlanarNonlinearEddy, nonlinear::PlanarNonlinearMagnetostatics,
};
use fs_matdb::{
    ClaimSet, MaterialCard, MaterialStateId, PropertyKey, PropertyValue, SelectionPolicy,
};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    field_storage::{
        nonlinear::REVERSIBLE_FIELD_STRENGTH_PROPERTY,
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::SPECIFIC_HEAT_CAPACITY_DIMS,
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::{Density, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const VOLUME: f64 = 0.3;
const DT: f64 = 0.1;
const CELL: CellId = CellId(7);
const MATERIAL: MaterialId = MaterialId(9);

fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 1,
        maximum_circuit_trials: 1000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-9,
        thermal_energy_tolerance_j: 1e-9,
        total_energy_tolerance_j: 2e-9,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-9,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-14,
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
            lower: &field::point(200.0),
            upper: &field::point(500.0),
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
                mass_kg: field::DENSITY_KG_M3 * VOLUME,
                specific_h_j_kg: field::SPECIFIC_HEAT_J_KG_K * (300.0 - 200.0),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
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
fn state<'a>(
    card: &MaterialCard,
    problem: &'a PlanarNonlinearMagnetostatics,
) -> PlanarNonlinearEddy<'a> {
    PlanarNonlinearEddy::new(
        &field::initial(problem),
        field::conductors(card, 4, 300.0, ElectricalPropertyBasis::Resistivity),
        &[],
        &field::boundary(0.0),
        0.0,
        field::controls(),
        &mut || false,
    )
    .unwrap()
}
fn cells<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a fs_matdb::QueryPoint,
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
    (0..4)
        .map(|triangle| MagneticRegion::Cell {
            triangle,
            cell: CELL,
        })
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 5e-8 * (1.0 + b.abs()), "{a:e} != {b:e}");
}

fn h(b: f64, saturating: bool) -> f64 {
    if saturating && b > 1.0 {
        1.0 + 4.0 * (b - 1.0)
    } else {
        b
    }
}
fn energy(b: f64, saturating: bool) -> f64 {
    if saturating && b > 1.0 {
        0.5 + (b - 1.0) + 2.0 * (b - 1.0).powi(2)
    } else {
        0.5 * b * b
    }
}
fn resistivity(t: f64) -> f64 {
    (1.0 + 0.005 * (t - 300.0)) / 18.0
}

struct Reference {
    temperature: f64,
    amplitude: f64,
    heat: f64,
    numerical: f64,
    energy: f64,
}
// Independent nested scalar bisection, with no production constitutive/field
// evaluation or Newton routines. C44 follows analytic P1 moments 1/6-1/9.
fn reference(a0: f64, t0: f64, saturating: bool) -> Reference {
    let amplitude = |t: f64| {
        let c = VOLUME / (18.0 * resistivity(t));
        let mut lo = 0.0;
        let mut hi = a0;
        for _ in 0..80 {
            let a = f64::midpoint(lo, hi);
            let residual = c * (a - a0) / DT + 2.0 * VOLUME * h(2.0 * a, saturating);
            if residual > 0.0 {
                hi = a;
            } else {
                lo = a;
            }
        }
        f64::midpoint(lo, hi)
    };
    let capacity = field::DENSITY_KG_M3 * field::SPECIFIC_HEAT_J_KG_K * VOLUME;
    let mut lo = t0;
    let mut hi = t0 + VOLUME * energy(2.0 * a0, saturating) / capacity;
    for _ in 0..80 {
        let t = f64::midpoint(lo, hi);
        let a = amplitude(t);
        let c = VOLUME / (18.0 * resistivity(t));
        let q = c * (a - a0).powi(2) / DT;
        if capacity * (t - t0) > q {
            hi = t;
        } else {
            lo = t;
        }
    }
    let temperature = f64::midpoint(lo, hi);
    let amplitude = amplitude(temperature);
    let heat = VOLUME / (18.0 * resistivity(temperature)) * (amplitude - a0).powi(2) / DT;
    let stored = VOLUME * energy(2.0 * amplitude, saturating);
    let numerical = VOLUME * h(2.0 * amplitude, saturating) * 2.0 * (amplitude - a0) - stored
        + VOLUME * energy(2.0 * a0, saturating);
    Reference {
        temperature,
        amplitude,
        heat,
        numerical,
        energy: stored,
    }
}

#[test]
fn actual_nonlinear_induction_and_enthalpy_match_an_independent_coupled_oracle() {
    let mut outcomes = Vec::new();
    for saturating in [false, true] {
        let knots = if saturating {
            vec![(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]
        } else {
            vec![(0.0, 0.0), (4.0, 4.0)]
        };
        let card = field::source(&knots, 18.0);
        let curve = caloric(&card);
        let thermal = thermal(&curve);
        let problem = field::problem(&card, field::mesh(VOLUME));
        let initial = state(&card, &problem);
        let context = field::point(300.0);
        let cells = cells(&card, &curve, &context);
        let routes = routes();
        let layout = MagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &routes,
            volume_tolerance_m3: 1e-12,
        };
        let candidate = layout
            .step_nonlinear_eddy_candidate(
                &thermal,
                DT,
                &initial,
                &[0.0; 4],
                controls(),
                field::controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let expected = reference(if saturating { 0.625 } else { 1.0 }, 300.0, saturating);
        let t = candidate.thermal().temperature_k(CELL).unwrap();
        let report = candidate.report();
        let electromagnetic = &report.electrical;
        near(t, expected.temperature);
        near(candidate.field().potential_wb_m()[4], expected.amplitude);
        near(electromagnetic.heat_j, expected.heat);
        near(electromagnetic.numerical_dissipation_j, expected.numerical);
        near(electromagnetic.stored_after_j, expected.energy);
        let capacity = field::DENSITY_KG_M3 * field::SPECIFIC_HEAT_J_KG_K * VOLUME;
        near(
            report.thermal.stored_energy_change_j,
            capacity * (t - 300.0),
        );
        near(
            report.thermal.stored_energy_change_j,
            electromagnetic.heat_j,
        );
        near(report.deposited_joule_j, electromagnetic.heat_j);
        near(report.non_electrical_heat_j, 0.0);
        near(electromagnetic.supplied_work_j, 0.0);
        near(
            electromagnetic.stored_before_j - electromagnetic.stored_after_j,
            electromagnetic.heat_j + electromagnetic.numerical_dissipation_j,
        );
        assert!(electromagnetic.numerical_dissipation_j > 1e-3);
        assert!(
            (electromagnetic.stored_before_j
                - electromagnetic.stored_after_j
                - report.thermal.stored_energy_change_j)
                .abs()
                > 1e-3
        );
        assert!(report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
        assert!(report.maximum_cell_residual_j <= controls().cell_energy_tolerance_j);
        assert!(report.circuit_trials > 1);
        for loss in &electromagnetic.cells {
            near(loss.heat_j, electromagnetic.heat_j / 4.0);
            near(loss.current_density_a_m2.iter().sum(), 0.0);
            assert!(loss.heat_j > 0.0);
            assert_eq!(
                loss.source.conductor.material().card_identity(),
                card.content_hash()
            );
            near(loss.source.conductor.resistivity_ohm_m(), resistivity(t));
            assert!(loss.source.conductor.resistivity_ohm_m() > 1.0 / 18.0);
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
        let frozen = initial
            .step_to_candidate(
                DT,
                initial.conductors().to_vec(),
                &[0.0; 4],
                field::controls(),
                &mut || false,
            )
            .unwrap();
        assert!(candidate.field().potential_wb_m()[4] < frozen.state().potential_wb_m()[4] - 1e-6);
        let next_expected = reference(expected.amplitude, expected.temperature, saturating);
        let continued = layout
            .step_nonlinear_eddy_candidate(
                candidate.thermal(),
                2.0 * DT,
                candidate.field(),
                &[0.0; 4],
                controls(),
                field::controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        near(
            continued.thermal().temperature_k(CELL).unwrap(),
            next_expected.temperature,
        );
        near(
            continued.field().potential_wb_m()[4],
            next_expected.amplitude,
        );
        near(
            continued.report().electrical.stored_before_j,
            electromagnetic.stored_after_j,
        );
        near(continued.report().electrical.heat_j, next_expected.heat);
        near(
            continued.report().electrical.numerical_dissipation_j,
            next_expected.numerical,
        );
        near(
            continued.field().energy_j()
                + capacity * (next_expected.temperature - 300.0)
                + electromagnetic.numerical_dissipation_j
                + continued.report().electrical.numerical_dissipation_j,
            initial.energy_j(),
        );
        assert!(std::ptr::eq(continued.field().problem(), initial.problem()));
        assert_eq!(continued.field().time_s(), 2.0 * DT);
        assert_eq!(continued.thermal().time_s(), 2.0 * DT);
        assert_eq!(candidate.field().time_s(), DT);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        near(thermal.temperature_k(CELL).unwrap(), 300.0);
        outcomes.push((
            t,
            electromagnetic.heat_j,
            candidate.field().potential_wb_m()[4],
        ));
    }
    assert!((outcomes[0].0 - outcomes[1].0).abs() > 0.01);
    assert!((outcomes[0].1 - outcomes[1].1).abs() > 0.01);
    assert!((outcomes[0].2 - outcomes[1].2).abs() > 0.1);
}

#[test]
fn missing_heat_route_failed_field_solve_and_late_cancel_leave_both_histories_unchanged() {
    let card = field::source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)], 18.0);
    let curve = caloric(&card);
    let thermal = thermal(&curve);
    let problem = field::problem(&card, field::mesh(VOLUME));
    let initial = state(&card, &problem);
    let old_a = initial.potential_wb_m().to_vec();
    let old_h = thermal.cells()[&CELL].specific_h_j_kg;
    let context = field::point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes();
    let layout = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let missing = MagneticThermalNetwork {
        magnetic_regions: &routes[..3],
        ..layout
    };
    assert!(
        missing
            .step_nonlinear_eddy_candidate(
                &thermal,
                DT,
                &initial,
                &[0.0; 4],
                controls(),
                field::controls(),
                &[0.0],
                || false
            )
            .is_err()
    );
    let mut limited = field::controls();
    limited.nonlinear.max_newton_iterations = 0;
    assert!(
        layout
            .step_nonlinear_eddy_candidate(
                &thermal,
                DT,
                &initial,
                &[0.0; 4],
                controls(),
                limited,
                &[0.0],
                || false
            )
            .is_err()
    );
    let mut checkpoints = 0;
    let success = layout
        .step_nonlinear_eddy_candidate(
            &thermal,
            DT,
            &initial,
            &[0.0; 4],
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
    let cancelled = layout.step_nonlinear_eddy_candidate(
        &thermal,
        DT,
        &initial,
        &[0.0; 4],
        controls(),
        field::controls(),
        &[0.0],
        || {
            visits += 1;
            visits == checkpoints
        },
    );
    assert!(matches!(cancelled, Err(CircuitEnthalpyError::Cancelled)));
    assert_eq!(initial.potential_wb_m(), old_a);
    assert_eq!(
        thermal.cells()[&CELL].specific_h_j_kg.to_bits(),
        old_h.to_bits()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert!(
        layout
            .step_nonlinear_eddy_candidate(
                &thermal,
                2.0 * DT,
                success.field(),
                &[0.0; 4],
                controls(),
                field::controls(),
                &[0.0],
                || false
            )
            .is_err()
    );
    let retry = layout
        .step_nonlinear_eddy_candidate(
            &thermal,
            DT,
            &initial,
            &[0.0; 4],
            controls(),
            field::controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(
        retry.field().potential_wb_m(),
        success.field().potential_wb_m()
    );
    assert_eq!(
        retry.thermal().cells()[&CELL].specific_h_j_kg.to_bits(),
        success.thermal().cells()[&CELL].specific_h_j_kg.to_bits()
    );
}

#[test]
fn temperature_dependent_h_b_refuses_thermal_feedback_without_hidden_magnetic_work() {
    let original = field::source(&[(0.0, 0.0), (4.0, 4.0)], 18.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key == PropertyKey::new(REVERSIBLE_FIELD_STRENGTH_PROPERTY, FIELD_STRENGTH_DIMS) {
            claim.value = PropertyValue::Surface2D {
                axes: ["B".into(), "T".into()],
                axis_dims: [FLUX_DENSITY_DIMS, Temperature::DIMS],
                samples: vec![
                    ([0.0, 200.0], 0.0),
                    ([4.0, 200.0], 4.0),
                    ([0.0, 500.0], 0.0),
                    ([4.0, 500.0], 8.0),
                ],
                triangles: vec![[0, 1, 3], [0, 3, 2]],
                dims: FIELD_STRENGTH_DIMS,
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conductor".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = caloric(&card);
    let thermal = thermal(&curve);
    let problem = field::problem(&card, field::mesh(VOLUME));
    let initial = state(&card, &problem);
    let context = field::point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes();
    let layout = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let result = layout.step_nonlinear_eddy_candidate(
        &thermal,
        DT,
        &initial,
        &[0.0; 4],
        controls(),
        field::controls(),
        &[0.0],
        || false,
    );
    assert!(matches!(
        result,
        Err(CircuitEnthalpyError::Source(
            MagneticThermalError::Magnetic(_)
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    near(thermal.temperature_k(CELL).unwrap(), 300.0);
}
