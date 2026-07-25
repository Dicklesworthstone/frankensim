//! Duty-cycle drivers: the windowed energy audit, its exactness condition,
//! the tiling refusal, and cycle summary quantities (bead
//! `frankensim-extreal-program-f85xj.5.12`).

mod support;

use fs_conduction::ConductionError;
use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::duty::{DutyCycle, DutySegment, lower_envelope_duty_cycle, march_duty_cycle};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::LinearConfig;
use fs_conduction::transient::{TransientConfig, TransientProblem, VolumetricHeatCapacity};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::envelope::{
    AxisPoint, DutyPoint, DutyWeighting, EnvelopeDutyCycle, EnvelopePoint,
};
use support::{with_cancelled_cx, with_cx};

const K: f64 = 10.0;
const RHO_CP: f64 = 2.0e6;
const T_COLD: f64 = 300.0;
const FULL_LOAD_DENSITY: f64 = 2.0e5; // W/m^3 at scale 1

fn unit_mesh(n: usize) -> ConductionMesh {
    let (complex, positions) = box_grid([n, n, n], [1.0, 1.0, 1.0]);
    let complex = TetComplex::from_tets(positions.len(), complex.tets);
    ConductionMesh::new(complex, positions).expect("unit box mesh")
}

fn linear_config() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-13,
        max_iterations: 60_000,
        restart: 60,
    }
}

fn capacity() -> VolumetricHeatCapacity {
    VolumetricHeatCapacity::declared(RHO_CP).expect("capacity")
}

fn material() -> ConductivityModel {
    ConductivityModel::isotropic_declared(K).expect("material")
}

fn cold_wall(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "cold",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(T_COLD).expect("cold"),
        )
        .expect("cold region")
        .adiabatic_remainder()
        .finish()
        .expect("partition")
}

/// Idle, ramp up, burst, ramp down — a plausible small cycle with both
/// interpolation kinds and a total window of 4000 s.
fn mixed_cycle() -> DutyCycle {
    DutyCycle::new(vec![
        DutySegment::constant(1000.0, 0.2).expect("idle"),
        DutySegment::ramp(1000.0, 0.2, 1.0).expect("ramp up"),
        DutySegment::constant(1000.0, 1.0).expect("burst"),
        DutySegment::ramp(1000.0, 1.0, 0.2).expect("ramp down"),
    ])
    .expect("cycle admits")
}

// ---------------------------------------------------------------------------
// The declared schedule.
// ---------------------------------------------------------------------------

#[test]
fn the_declared_energy_integral_is_analytic_per_segment() {
    let cycle = mixed_cycle();
    assert!((cycle.window_s() - 4000.0).abs() < 1e-12);
    // 0.2*1000 + mean(0.2,1.0)*1000 + 1.0*1000 + mean(1.0,0.2)*1000
    let expected = 0.2f64.mul_add(1000.0, 0.6 * 1000.0) + 1.0f64.mul_add(1000.0, 0.6 * 1000.0);
    assert!(
        (cycle.energy_scale_seconds() - expected).abs() < 1e-12,
        "declared integral {} != {expected}",
        cycle.energy_scale_seconds()
    );
}

#[test]
fn the_schedule_interpolates_as_declared() {
    let cycle = mixed_cycle();
    // Constant segment holds.
    assert!((cycle.scale_at(0.0).expect("t0") - 0.2).abs() < 1e-12);
    assert!((cycle.scale_at(500.0).expect("mid idle") - 0.2).abs() < 1e-12);
    // Linear segment ramps.
    assert!((cycle.scale_at(1500.0).expect("mid ramp") - 0.6).abs() < 1e-12);
    assert!((cycle.scale_at(2000.0).expect("top") - 1.0).abs() < 1e-12);
    // And the far end.
    assert!((cycle.scale_at(4000.0).expect("end") - 0.2).abs() < 1e-12);
}

#[test]
fn a_time_past_the_window_refuses_instead_of_holding() {
    // Holding the last value is the tempting default and it is wrong: it
    // invents load history that nobody declared.
    let cycle = mixed_cycle();
    match cycle
        .scale_at(4000.1)
        .expect_err("past the window must refuse")
    {
        ConductionError::ScenarioRow { what, .. } => {
            assert!(what.contains("past the declared"), "what: {what}");
            assert!(
                what.contains("does not hold its last value"),
                "the refusal must say WHY: {what}"
            );
        }
        other => panic!("expected a scenario-row refusal, got {other:?}"),
    }
    assert!(cycle.scale_at(-1.0).is_err(), "negative time");
    assert!(cycle.scale_at(f64::NAN).is_err(), "non-finite time");
}

#[test]
fn segment_and_cycle_admission_refuse_degenerate_declarations() {
    assert!(DutySegment::constant(0.0, 1.0).is_err(), "zero duration");
    assert!(
        DutySegment::constant(-1.0, 1.0).is_err(),
        "negative duration"
    );
    assert!(
        DutySegment::constant(1.0, -0.1).is_err(),
        "a negative scale would make the source a sink"
    );
    assert!(DutySegment::ramp(1.0, 0.0, f64::NAN).is_err(), "non-finite");
    assert!(DutyCycle::new(Vec::new()).is_err(), "empty cycle");
    assert!(
        DutySegment::constant(1.0, 0.0).is_ok(),
        "an off state is legal"
    );
}

// ---------------------------------------------------------------------------
// The windowed energy audit.
// ---------------------------------------------------------------------------

fn march(
    cycle: &DutyCycle,
    config: &TransientConfig,
    limit_k: f64,
) -> fs_conduction::duty::DutyCycleSolution {
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", FULL_LOAD_DENSITY).expect("source");
    // Unit box, so the full-load total power is the density times 1 m^3.
    let base_power_w = FULL_LOAD_DENSITY;
    let initial = vec![T_COLD; mesh.vertex_count()];
    with_cx(|cx| {
        march_duty_cycle(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            config,
            cycle,
            base_power_w,
            &initial,
            limit_k,
        )
        .expect("duty march")
    })
}

#[test]
fn crank_nicolson_with_aligned_steps_delivers_the_declared_energy_exactly() {
    // The theta method's load weighting at theta = 0.5 IS the trapezoid rule,
    // and trapezoid is exact for a piecewise-linear profile when no step
    // straddles a segment boundary. So this is an identity, not a tolerance:
    // if it ever drifts, the weighting or the schedule sampling is wrong.
    let cycle = mixed_cycle();
    let dt = 250.0; // 4000 / 250 = 16 steps, every boundary on a step
    assert!(cycle.steps_align(dt, 1e-9), "the fixture must be aligned");
    let config = TransientConfig::crank_nicolson(dt, linear_config()).expect("config");

    let solution = march(&cycle, &config, f64::INFINITY);
    let energy = solution.energy;
    assert!(energy.steps_aligned());
    assert!(
        energy.residual_j().abs() < 1e-6,
        "declared {} J vs delivered {} J (residual {} J) must agree exactly",
        energy.declared_j(),
        energy.delivered_j(),
        energy.residual_j()
    );
    // And the declared value is the analytic one, not a re-quadrature.
    assert!((energy.declared_j() - 2400.0 * FULL_LOAD_DENSITY).abs() < 1e-6);
}

/// A cycle whose ramps do NOT mirror each other, so first-order integration
/// error cannot cancel between them.
fn asymmetric_cycle() -> DutyCycle {
    DutyCycle::new(vec![
        DutySegment::ramp(1000.0, 0.0, 1.0).expect("ramp up"),
        DutySegment::constant(3000.0, 1.0).expect("hold"),
    ])
    .expect("cycle admits")
}

#[test]
fn backward_euler_leaves_a_real_first_order_residual() {
    // The audit must not hide the scheme's own integration error. Backward
    // Euler is the right-endpoint rule, so on a ramp it is genuinely off —
    // and halving the step must halve the gap.
    let cycle = asymmetric_cycle();
    let residual_at = |dt: f64| {
        let config = TransientConfig::backward_euler(dt, linear_config()).expect("config");
        march(&cycle, &config, f64::INFINITY)
            .energy
            .residual_j()
            .abs()
    };
    let coarse = residual_at(500.0);
    let fine = residual_at(250.0);
    assert!(
        coarse > 1.0,
        "backward Euler on a ramp must show a real residual, got {coarse} J"
    );
    let ratio = coarse / fine;
    assert!(
        (1.6..2.4).contains(&ratio),
        "halving the step should roughly halve the first-order residual, ratio {ratio}"
    );
}

#[test]
fn a_symmetric_cycle_hides_the_first_order_error_by_cancellation() {
    // Found while writing the test above, and worth pinning: on a cycle whose
    // ramps mirror each other, backward Euler's right-endpoint rule
    // OVER-integrates the rising ramp by exactly what it UNDER-integrates the
    // falling one, so the windowed residual is zero even though the scheme is
    // first order and every individual step is wrong.
    //
    // The consequence matters for anyone reading an audit: a zero energy
    // residual is NOT evidence that the integration is exact. It is evidence
    // that the errors summed to zero over this particular window, which a
    // symmetric duty cycle — an extremely common shape — arranges for free.
    let symmetric = mixed_cycle();
    let config = TransientConfig::backward_euler(500.0, linear_config()).expect("config");
    let residual = march(&symmetric, &config, f64::INFINITY)
        .energy
        .residual_j();
    assert!(
        residual.abs() < 1e-6,
        "the symmetric fixture is expected to cancel, got {residual} J"
    );

    // Same scheme, same step, asymmetric profile: the error is exposed.
    let exposed = march(&asymmetric_cycle(), &config, f64::INFINITY)
        .energy
        .residual_j()
        .abs();
    assert!(
        exposed > 1.0,
        "the asymmetric fixture must expose what the symmetric one hides, got {exposed} J"
    );
}

#[test]
fn a_window_that_is_not_a_whole_number_of_steps_refuses() {
    // A step straddling the window end would inject load the schedule does
    // not declare, so it is refused rather than truncated or extended.
    let cycle = mixed_cycle();
    let config = TransientConfig::crank_nicolson(300.0, linear_config()).expect("config");
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", FULL_LOAD_DENSITY).expect("source");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let error = with_cx(|cx| {
        march_duty_cycle(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &cycle,
            FULL_LOAD_DENSITY,
            &initial,
            f64::INFINITY,
        )
        .expect_err("a ragged window must refuse")
    });
    match error {
        ConductionError::ScenarioRow { what, .. } => {
            assert!(what.contains("whole number"), "what: {what}");
        }
        other => panic!("expected a scenario-row refusal, got {other:?}"),
    }
}

#[test]
fn a_constant_full_load_cycle_matches_the_plain_transient_energy() {
    // Sanity against the simplest case: an all-ones cycle over the window
    // must deliver exactly base_power * window, for either scheme, because
    // both quadratures are exact on a constant.
    let cycle =
        DutyCycle::new(vec![DutySegment::constant(2000.0, 1.0).expect("flat")]).expect("cycle");
    for config in [
        TransientConfig::crank_nicolson(250.0, linear_config()).expect("cn"),
        TransientConfig::backward_euler(250.0, linear_config()).expect("be"),
    ] {
        let energy = march(&cycle, &config, f64::INFINITY).energy;
        let expected = 2000.0 * FULL_LOAD_DENSITY;
        assert!(
            (energy.delivered_j() - expected).abs() < 1e-6,
            "constant load must integrate exactly under either scheme: {} vs {expected}",
            energy.delivered_j()
        );
    }
}

// ---------------------------------------------------------------------------
// Cycle summary quantities.
// ---------------------------------------------------------------------------

#[test]
fn the_cycle_summary_reports_a_peak_that_lags_the_load() {
    // Thermal mass means the hottest moment trails the burst rather than
    // coinciding with it. If the peak ever landed at the load peak exactly,
    // the capacitance would not be doing anything.
    let cycle = mixed_cycle();
    let config = TransientConfig::crank_nicolson(125.0, linear_config()).expect("config");
    let solution = march(&cycle, &config, f64::INFINITY);
    let summary = &solution.summary;

    assert!(
        summary.peak_temperature_k() > T_COLD,
        "the body must warm above the cold wall"
    );
    assert!(
        summary.excursion_k() > 0.0,
        "excursion must be positive for a heated cycle"
    );
    // The burst ends at t = 3000; the peak must not precede it, and must fall
    // inside the window.
    assert!(
        summary.peak_time_s() >= 3000.0 - 1e-9,
        "peak at {} s arrived before the burst ended",
        summary.peak_time_s()
    );
    assert!(summary.peak_time_s() <= cycle.window_s() + 1e-9);
}

#[test]
fn time_above_limit_is_step_quantized_and_monotone_in_the_limit() {
    let cycle = mixed_cycle();
    let dt = 125.0;
    let config = TransientConfig::crank_nicolson(dt, linear_config()).expect("config");

    let unreachable = march(&cycle, &config, 1.0e6);
    assert_eq!(unreachable.summary.steps_above_limit(), 0);
    assert!(unreachable.summary.time_above_limit_s().abs() < 1e-12);

    let always = march(&cycle, &config, 0.0);
    assert!(always.summary.steps_above_limit() > 0);
    // Quantized to whole steps, by construction.
    let expected = dt * (always.summary.steps_above_limit() as f64);
    assert!((always.summary.time_above_limit_s() - expected).abs() < 1e-9);

    // A lower limit cannot be exceeded for less time than a higher one.
    let peak = unreachable.summary.peak_temperature_k();
    let low = march(&cycle, &config, T_COLD + (peak - T_COLD) * 0.25);
    let high = march(&cycle, &config, T_COLD + (peak - T_COLD) * 0.75);
    assert!(
        low.summary.time_above_limit_s() >= high.summary.time_above_limit_s(),
        "time above a LOWER limit must not be shorter"
    );
}

// ---------------------------------------------------------------------------
// Refusals, determinism, cancellation.
// ---------------------------------------------------------------------------

#[test]
fn a_negative_base_power_refuses() {
    let cycle =
        DutyCycle::new(vec![DutySegment::constant(1000.0, 1.0).expect("flat")]).expect("cycle");
    let config = TransientConfig::crank_nicolson(250.0, linear_config()).expect("config");
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", FULL_LOAD_DENSITY).expect("source");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let error = with_cx(|cx| {
        march_duty_cycle(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &cycle,
            -1.0,
            &initial,
            f64::INFINITY,
        )
        .expect_err("negative base power")
    });
    assert!(matches!(error, ConductionError::ScenarioRow { .. }));
}

#[test]
fn duty_marching_is_deterministic() {
    let cycle = mixed_cycle();
    let config = TransientConfig::crank_nicolson(250.0, linear_config()).expect("config");
    assert_eq!(
        march(&cycle, &config, T_COLD + 1.0),
        march(&cycle, &config, T_COLD + 1.0)
    );
}

#[test]
fn a_cancelled_duty_march_publishes_nothing() {
    let cycle = mixed_cycle();
    let config = TransientConfig::crank_nicolson(250.0, linear_config()).expect("config");
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", FULL_LOAD_DENSITY).expect("source");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let error = with_cancelled_cx(|cx| {
        march_duty_cycle(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &cycle,
            FULL_LOAD_DENSITY,
            &initial,
            f64::INFINITY,
        )
        .expect_err("a cancelled march publishes nothing")
    });
    assert!(matches!(error, ConductionError::Cancelled { .. }));
}

// ---------------------------------------------------------------------------
// The E17 seam: fs_scenario envelope duty cycles lowered onto the schedule
// ---------------------------------------------------------------------------

const AMBIENT: Dims = Dims([0, 0, 0, 1, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

/// A scenario point at a given power, holding ambient and fan state fixed.
fn scenario_point(power_w: f64, ambient_k: f64, fan: &str) -> EnvelopePoint {
    EnvelopePoint {
        coordinates: vec![
            (
                "ambient".to_string(),
                AxisPoint::Continuous(QtyAny::new(ambient_k, AMBIENT)),
            ),
            (
                "total-power".to_string(),
                AxisPoint::Continuous(QtyAny::new(power_w, POWER)),
            ),
            (
                "fan-state".to_string(),
                AxisPoint::Discrete(fan.to_string()),
            ),
        ],
    }
}

/// Idle at a quarter load for 3000 s, then burst at full load for 1000 s.
fn scenario_cycle() -> EnvelopeDutyCycle {
    EnvelopeDutyCycle {
        name: "office-day".to_string(),
        weighting: DutyWeighting::Dwells,
        points: vec![
            DutyPoint {
                name: "idle".to_string(),
                point: scenario_point(FULL_LOAD_DENSITY * 0.25, 300.0, "both-running"),
                weight: QtyAny::new(3000.0, TIME),
            },
            DutyPoint {
                name: "burst".to_string(),
                point: scenario_point(FULL_LOAD_DENSITY, 300.0, "both-running"),
                weight: QtyAny::new(1000.0, TIME),
            },
        ],
    }
}

#[test]
fn a_scenario_duty_cycle_lowers_onto_the_schedule_with_the_declared_scales() {
    let lowered = lower_envelope_duty_cycle(&scenario_cycle(), "total-power", FULL_LOAD_DENSITY)
        .expect("a dwell-weighted cycle over a power axis lowers");

    assert_eq!(lowered.base_power_w, FULL_LOAD_DENSITY);
    assert_eq!(lowered.cycle.segments().len(), 2);
    assert_eq!(lowered.cycle.boundaries_s(), &[0.0, 3000.0, 4000.0]);
    assert_eq!(lowered.cycle.scale_at(0.0).expect("in window"), 0.25);
    assert_eq!(lowered.cycle.scale_at(3500.0).expect("in window"), 1.0);

    // The dimensionless energy integral: 0.25*3000 + 1.0*1000 = 1750 s.
    assert_eq!(lowered.cycle.energy_scale_seconds(), 1750.0);
}

#[test]
fn lowering_records_the_axes_it_held_rather_than_applied() {
    let lowered = lower_envelope_duty_cycle(&scenario_cycle(), "total-power", FULL_LOAD_DENSITY)
        .expect("lowers");
    let names: Vec<&str> = lowered
        .held_axes
        .iter()
        .map(|(axis, _)| axis.as_str())
        .collect();
    assert_eq!(
        names,
        vec!["ambient", "fan-state"],
        "the caller can check that the boundary conditions they built match the conditions the \
         cycle was declared at"
    );
    assert_eq!(lowered.held_axes[1].1, "both-running");
}

#[test]
fn a_cycle_whose_ambient_varies_refuses_rather_than_marching_at_one_ambient() {
    // This is the failure the seam exists to prevent. A DutyCycle is a scale
    // on a volumetric source, and boundary conditions are assembled ONCE, so
    // applying the power axis while ignoring a varying ambient would march a
    // hot-afternoon cycle at morning ambient and report a confidently wrong
    // peak junction temperature. Nothing about that is visible downstream.
    let mut cycle = scenario_cycle();
    cycle.points[1].point = scenario_point(FULL_LOAD_DENSITY, 318.0, "both-running");

    let error = lower_envelope_duty_cycle(&cycle, "total-power", FULL_LOAD_DENSITY)
        .expect_err("a varying non-power axis must refuse");
    let ConductionError::ScenarioRow { region, what, .. } = &error else {
        panic!("expected a scenario-row refusal, got {error:?}");
    };
    assert_eq!(region, "ambient", "the refusal names the offending axis");
    assert!(
        what.contains("only the power axis is lowered"),
        "and says why it cannot be carried: {what}"
    );
    assert!(
        what.contains("idle") && what.contains("burst"),
        "and names both dwells that disagree: {what}"
    );
}

#[test]
fn a_cycle_whose_discrete_state_varies_also_refuses() {
    // Same hazard through the discrete axis: a fan failure changes the Robin
    // coefficient, which lives in the boundary, not the source.
    let mut cycle = scenario_cycle();
    cycle.points[1].point = scenario_point(FULL_LOAD_DENSITY, 300.0, "fan-a-failed");

    let error = lower_envelope_duty_cycle(&cycle, "total-power", FULL_LOAD_DENSITY)
        .expect_err("a varying fan state must refuse");
    let ConductionError::ScenarioRow { region, .. } = &error else {
        panic!("expected a scenario-row refusal, got {error:?}");
    };
    assert_eq!(region, "fan-state");
}

#[test]
fn a_fraction_weighted_scenario_cycle_refuses_because_it_has_no_absolute_time() {
    let mut cycle = scenario_cycle();
    cycle.weighting = DutyWeighting::Fractions;
    cycle.points[0].weight = QtyAny::dimensionless(0.75);
    cycle.points[1].weight = QtyAny::dimensionless(0.25);

    let error = lower_envelope_duty_cycle(&cycle, "total-power", FULL_LOAD_DENSITY)
        .expect_err("fractions carry no duration");
    let ConductionError::ScenarioRow { region, what, .. } = &error else {
        panic!("expected a scenario-row refusal, got {error:?}");
    };
    assert_eq!(region, "weighting");
    assert!(
        what.contains("no absolute dwell duration"),
        "the refusal explains what is missing, not merely that it failed: {what}"
    );
}

#[test]
fn lowering_refuses_a_missing_or_wrongly_dimensioned_power_axis() {
    let missing = lower_envelope_duty_cycle(&scenario_cycle(), "dissipation", FULL_LOAD_DENSITY)
        .expect_err("no such axis");
    let ConductionError::ScenarioRow { region, .. } = &missing else {
        panic!("expected a scenario-row refusal");
    };
    assert_eq!(region, "dissipation");

    // Pointing the power axis at the temperature axis: dimensions catch it.
    let wrong = lower_envelope_duty_cycle(&scenario_cycle(), "ambient", FULL_LOAD_DENSITY)
        .expect_err("kelvin is not watts");
    let ConductionError::ScenarioRow { what, .. } = &wrong else {
        panic!("expected a scenario-row refusal");
    };
    assert!(what.contains("is not power"), "got {what}");
}

#[test]
fn lowering_refuses_a_non_positive_reference_power() {
    for base in [0.0, -1.0, f64::NAN] {
        let error = lower_envelope_duty_cycle(&scenario_cycle(), "total-power", base)
            .expect_err("the reference power must be positive and finite");
        let ConductionError::ScenarioRow { region, .. } = &error else {
            panic!("expected a scenario-row refusal");
        };
        assert_eq!(region, "base_power_w");
    }
}

/// The seam EXERCISED, not assumed: a scenario-declared cycle drives a real
/// transient march, and the audit's declared energy is the SCENARIO's own.
#[test]
fn a_scenario_declared_cycle_drives_a_transient_solve() {
    let scenario = scenario_cycle();
    let lowered =
        lower_envelope_duty_cycle(&scenario, "total-power", FULL_LOAD_DENSITY).expect("lowers");

    let dt = 250.0;
    assert!(lowered.cycle.steps_align(dt, 1e-9));
    let config = TransientConfig::crank_nicolson(dt, linear_config()).expect("config");
    let solution = march(&lowered.cycle, &config, f64::INFINITY);
    assert!(solution.energy.steps_aligned());

    // The declared energy is the scenario's own arithmetic, computed here
    // independently of the lowering: sum of (point power x dwell).
    let from_scenario: f64 = scenario
        .points
        .iter()
        .map(|dwell| {
            dwell
                .point
                .continuous("total-power")
                .expect("power present")
                .value
                * dwell.weight.value
        })
        .sum();
    assert_eq!(from_scenario, 1750.0 * FULL_LOAD_DENSITY);
    assert!(
        (solution.energy.declared_j() - from_scenario).abs() < 1e-6,
        "the audit's declared energy is the scenario's own, so the seam is checked rather than \
         the lowering checked against itself: {} vs {from_scenario}",
        solution.energy.declared_j()
    );

    // And the physics ran: the body heated above its cold wall.
    assert!(
        solution.summary.peak_temperature_k() > T_COLD,
        "a heated body must rise above the wall it is clamped to"
    );
    assert!(
        solution.summary.peak_time_s() > 3000.0,
        "the peak is in the burst"
    );
}

#[test]
fn the_lowered_cycle_agrees_with_the_same_schedule_written_by_hand() {
    // The lowering must be a pure re-expression: no smoothing, no reordering,
    // no renormalisation. A scenario cycle declares dwells AT points, not
    // transitions BETWEEN them, so it lowers to piecewise-CONSTANT segments.
    let lowered = lower_envelope_duty_cycle(&scenario_cycle(), "total-power", FULL_LOAD_DENSITY)
        .expect("lowers");
    let by_hand = DutyCycle::new(vec![
        DutySegment::constant(3000.0, 0.25).expect("segment"),
        DutySegment::constant(1000.0, 1.0).expect("segment"),
    ])
    .expect("cycle");
    assert_eq!(lowered.cycle, by_hand);
}

/// A scenario-declared cycle CANNOT deliver its declared energy exactly, and
/// the gap is analytic rather than numerical noise.
///
/// This REVERSES the exactness property `mixed_cycle` demonstrates, and the
/// reason is structural rather than incidental: `mixed_cycle` is CONTINUOUS
/// (constant -> ramp -> constant -> ramp), while a scenario duty cycle
/// declares dwells AT points with no transitions between them, so it lowers to
/// a schedule with a JUMP. The theta weighting is the trapezoid rule at
/// theta = 0.5, exact for piecewise-LINEAR profiles — and a jump is not one.
///
/// Across a boundary where the scale jumps by `ds`, the endpoint weighting
/// misses exactly `(ds / 2) * dt` scale-seconds. Here `ds = 0.75` and
/// `dt = 250`, so the shortfall is `93.75` scale-seconds.
#[test]
fn the_discontinuity_a_scenario_cycle_introduces_costs_exactly_the_analytic_jump_error() {
    let lowered = lower_envelope_duty_cycle(&scenario_cycle(), "total-power", FULL_LOAD_DENSITY)
        .expect("lowers");

    let dt = 250.0;
    let config = TransientConfig::crank_nicolson(dt, linear_config()).expect("config");
    let solution = march(&lowered.cycle, &config, f64::INFINITY);

    let jump = 1.0 - 0.25;
    let predicted = -(jump / 2.0) * dt * FULL_LOAD_DENSITY;
    assert!(
        (solution.energy.residual_j() - predicted).abs() < 1e-6,
        "the residual must equal the analytic jump error {predicted} J, not merely be small: \
         got {}",
        solution.energy.residual_j()
    );
    assert!(
        solution.energy.residual_j() < 0.0,
        "the endpoint weighting UNDER-delivers across a rising jump, because the step containing \
         the jump is charged the pre-jump scale at its left end"
    );
}

#[test]
fn the_jump_error_is_first_order_in_the_step_so_it_is_resolution_not_a_lowering_bug() {
    let lowered = lower_envelope_duty_cycle(&scenario_cycle(), "total-power", FULL_LOAD_DENSITY)
        .expect("lowers");

    let coarse = march(
        &lowered.cycle,
        &TransientConfig::crank_nicolson(250.0, linear_config()).expect("config"),
        f64::INFINITY,
    );
    let fine = march(
        &lowered.cycle,
        &TransientConfig::crank_nicolson(125.0, linear_config()).expect("config"),
        f64::INFINITY,
    );

    let ratio = coarse.energy.residual_j() / fine.energy.residual_j();
    assert!(
        (ratio - 2.0).abs() < 1e-9,
        "halving the step must halve the gap exactly, which is what makes it a resolution \
         artifact rather than a defect in the lowering: ratio {ratio}"
    );
}

#[test]
fn declaring_the_transition_restores_exactness() {
    // The remedy available to a caller who needs an exact windowed balance: a
    // ramp is continuous, so the trapezoid rule integrates it exactly. Same
    // window and same total energy as the stepped cycle, zero residual.
    let ramped = DutyCycle::new(vec![
        DutySegment::constant(2875.0, 0.25).expect("idle"),
        DutySegment::ramp(250.0, 0.25, 1.0).expect("transition"),
        DutySegment::constant(875.0, 1.0).expect("burst"),
    ])
    .expect("cycle");
    // 0.25*2875 + 0.625*250 + 1.0*875 = 718.75 + 156.25 + 875 = 1750, the same
    // dimensionless energy as the stepped form.
    assert_eq!(ramped.energy_scale_seconds(), 1750.0);

    let dt = 125.0;
    assert!(ramped.steps_align(dt, 1e-9));
    let config = TransientConfig::crank_nicolson(dt, linear_config()).expect("config");
    let solution = march(&ramped, &config, f64::INFINITY);
    assert!(
        solution.energy.residual_j().abs() < 1e-6,
        "a continuous schedule balances exactly: residual {} J",
        solution.energy.residual_j()
    );
}
