//! Duty-cycle drivers: exact prescribed source energy across ramps and jumps,
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
use fs_conduction::transient::{
    TransientConfig, TransientProblem, VolumetricHeatCapacity, assemble_capacitance,
};
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
fn interval_energy_splits_at_switches_and_ramp_endpoints() {
    let cycle = DutyCycle::new(vec![
        DutySegment::constant(3.0, 0.25).expect("idle"),
        DutySegment::ramp(2.0, 0.0, 2.0).expect("ramp"),
        DutySegment::constant(5.0, 0.5).expect("hold"),
    ])
    .expect("cycle");
    // [2,3]: 0.25; [3,5]: area 2; [5,7]: area 1.
    assert_eq!(
        cycle.energy_scale_between(2.0, 7.0).expect("interval"),
        3.25
    );
    assert_eq!(
        cycle.energy_scale_between(3.5, 4.5).expect("ramp slice"),
        1.0
    );
    assert_eq!(cycle.energy_scale_between(3.0, 3.0).expect("empty"), 0.0);
    assert_eq!(cycle.energy_scale_between(10.0, 10.0).expect("end"), 0.0);
    assert_eq!(
        cycle.energy_scale_between(0.0, 10.0).expect("whole cycle"),
        cycle.energy_scale_seconds()
    );
    assert!(cycle.energy_scale_between(4.0, 3.0).is_err());
    assert!(cycle.energy_scale_between(-1.0, 3.0).is_err());
    assert!(cycle.energy_scale_between(9.0, 11.0).is_err());
    assert!(cycle.energy_scale_between(f64::NAN, 3.0).is_err());
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
        DutyCycle::new(vec![
            DutySegment::constant(1.0e20, 0.0).expect("long idle"),
            DutySegment::constant(1.0, 1.0).expect("pulse"),
        ])
        .is_err(),
        "a duration lost to the accumulated time must not erase a pulse"
    );
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
    // The declared segment integrals and delivered interval integrals agree.
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

#[test]
fn switching_off_power_preserves_the_robin_bath_equilibrium() {
    // The power map can turn off while a real 300 K environment remains.
    // Scaling the complete assembled load would instead make the body cool
    // toward zero kelvin by removing the bath's incoming Robin heat.
    let mesh = unit_mesh(3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "bath",
            |_| true,
            ThermalBc::robin(20.0, T_COLD).expect("Robin"),
        )
        .expect("bath region")
        .finish()
        .expect("boundary");
    let source = ScalarField::uniform("full-load source", FULL_LOAD_DENSITY).expect("source");
    let cycle =
        DutyCycle::new(vec![DutySegment::constant(4000.0, 0.0).expect("off")]).expect("cycle");
    let initial = vec![T_COLD; mesh.vertex_count()];

    for theta in [1.0, 0.75, 0.5] {
        let config = TransientConfig::new(theta, 1000.0, linear_config()).expect("config");
        let solution = with_cx(|cx| {
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
            .expect("power-off march")
        });
        for (vertex, &temperature) in solution.temperature.iter().enumerate() {
            assert!(
                (temperature - T_COLD).abs() < 2e-10,
                "theta={theta}: power-off vertex {vertex} left bath equilibrium: {temperature} K"
            );
        }
        assert_eq!(solution.energy.delivered_j(), 0.0);
    }
}

#[test]
fn prescribed_boundary_heat_survives_a_varying_power_schedule() {
    // Integrate the actual stored heat of a unit cube with a known inward
    // boundary flux. The affine P1 volume source has exact mean 4 W/m^3;
    // only that source follows the ramp. The six faces supply 12 W throughout.
    let mesh = unit_mesh(3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "inward heat",
            |_| true,
            ThermalBc::neumann(-2.0).expect("flux"),
        )
        .expect("flux region")
        .finish()
        .expect("boundary");
    let source = ScalarField::nodal(
        "affine full-load source",
        mesh.vertex_count(),
        mesh.positions().iter().map(|p| 2.0 + 4.0 * p[0]).collect(),
    )
    .expect("source");
    let cycle =
        DutyCycle::new(vec![DutySegment::ramp(6000.0, 0.0, 1.0).expect("ramp")]).expect("cycle");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let capacitance = with_cx(|cx| assemble_capacitance(cx, &mesh, capacity()).expect("capacity"));
    let dt_s = 2000.0;

    for theta in [1.0, 0.75, 0.5] {
        let config = TransientConfig::new(theta, dt_s, linear_config()).expect("config");
        let solution = with_cx(|cx| {
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
                4.0,
                &initial,
                f64::INFINITY,
            )
            .expect("heated boundary march")
        });
        // Prescribed ramp energy is independent of the transport integrator.
        let volume_heat_j = 4.0 * 3000.0;
        let boundary_heat_j = 12.0 * 6000.0;
        let change: Vec<_> = solution.temperature.iter().map(|t| t - T_COLD).collect();
        let mut nodal_energy = vec![0.0; mesh.vertex_count()];
        capacitance.spmv(&change, &mut nodal_energy);
        let stored_j: f64 = nodal_energy.iter().sum();
        let expected_j = volume_heat_j + boundary_heat_j;
        assert!(
            (stored_j - expected_j).abs() < 1e-5,
            "theta={theta}: stored {stored_j} J, expected {expected_j} J with the fixed boundary heat"
        );
        assert!((solution.energy.delivered_j() - volume_heat_j).abs() < 1e-9);
    }
}

#[test]
fn a_millikelvin_scheduled_pulse_is_resolved_against_its_temperature_rise() {
    let mesh = unit_mesh(4);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .expect("insulated boundary");
    let source_density = 2.0;
    let source = ScalarField::uniform("weak full-load source", source_density).expect("source");
    let cycle =
        DutyCycle::new(vec![DutySegment::ramp(6000.0, 0.0, 1.0).expect("ramp")]).expect("cycle");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let linear = LinearConfig {
        tolerance: 1e-7,
        ..LinearConfig::default()
    };
    let dt_s = 2000.0;

    for theta in [1.0, 0.75, 0.5] {
        let config = TransientConfig::new(theta, dt_s, linear).expect("config");
        let solution = with_cx(|cx| {
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
                source_density,
                &initial,
                f64::INFINITY,
            )
            .expect("weak scheduled pulse")
        });
        let expected_rise = source_density * 3000.0 / RHO_CP;
        for (vertex, &temperature) in solution.temperature.iter().enumerate() {
            let error = (temperature - T_COLD - expected_rise).abs();
            assert!(
                error < 20.0 * linear.tolerance * expected_rise,
                "theta={theta}: vertex {vertex} weak-pulse error {error} K exceeds the tolerance on the {expected_rise} K rise"
            );
        }
    }
}

/// An asymmetric ramp catches source-quadrature errors that cancel in a
/// symmetric up/down cycle.
fn asymmetric_cycle() -> DutyCycle {
    DutyCycle::new(vec![
        DutySegment::ramp(1000.0, 0.0, 1.0).expect("ramp up"),
        DutySegment::constant(3000.0, 1.0).expect("hold"),
    ])
    .expect("cycle admits")
}

#[test]
fn backward_euler_delivers_the_declared_ramp_energy_at_every_step_size() {
    let cycle = asymmetric_cycle();
    for dt in [500.0, 250.0, 800.0] {
        let config = TransientConfig::backward_euler(dt, linear_config()).expect("config");
        let energy = march(&cycle, &config, f64::INFINITY).energy;
        assert!(energy.residual_j().abs() < 1e-6, "dt={dt}: {energy:?}");
    }
}

#[test]
fn a_pulse_between_sample_times_heats_an_insulated_body_by_its_actual_joules() {
    let mesh = unit_mesh(2);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .expect("insulated boundary");
    let source = ScalarField::uniform("full-load source", FULL_LOAD_DENSITY).expect("source");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let cycle = DutyCycle::new(vec![
        DutySegment::constant(2.0, 0.0).expect("idle"),
        DutySegment::ramp(1.0, 0.0, 2.0).expect("ramp"),
        DutySegment::constant(0.25, 4.0).expect("burst"),
        DutySegment::constant(6.75, 0.0).expect("off"),
    ])
    .expect("cycle");
    // One scale-second in the ramp plus one in the short burst. In an
    // insulated homogeneous body, calorimetry gives ΔT = Q/(rho Cp V).
    let expected_j = 2.0 * FULL_LOAD_DENSITY;
    let expected_k = T_COLD + expected_j / RHO_CP;
    for theta in [1.0, 0.75, 0.5] {
        for dt in [10.0, 2.0] {
            let config = TransientConfig::new(theta, dt, linear_config()).expect("config");
            let solution = with_cx(|cx| {
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
                .expect("pulse march")
            });
            assert!(!solution.energy.steps_aligned());
            assert!((solution.energy.delivered_j() - expected_j).abs() < 1e-6);
            for &temperature in &solution.temperature {
                assert!(
                    (temperature - expected_k).abs() < 1e-9,
                    "theta={theta}, dt={dt}: {temperature} K instead of {expected_k} K"
                );
            }
        }
    }
}

#[test]
fn the_energy_audit_refuses_watts_unrelated_to_the_assembled_source() {
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("full-load source", FULL_LOAD_DENSITY).expect("source");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let cycle =
        DutyCycle::new(vec![DutySegment::constant(10.0, 1.0).expect("hold")]).expect("cycle");
    let config = TransientConfig::backward_euler(10.0, linear_config()).expect("config");
    for declared_w in [0.0, 0.5 * FULL_LOAD_DENSITY, 2.0 * FULL_LOAD_DENSITY] {
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
                declared_w,
                &initial,
                f64::INFINITY,
            )
            .expect_err("the audit cannot relabel the actual source power")
        });
        assert!(matches!(error, ConductionError::ScenarioRow { what, .. }
            if what.contains("assembled full-load source")));
    }
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

#[test]
fn scenario_dwells_deliver_their_authored_energy_without_smoothing_the_jump() {
    let lowered = lower_envelope_duty_cycle(&scenario_cycle(), "total-power", FULL_LOAD_DENSITY)
        .expect("lowers");

    let dt = 250.0;
    let config = TransientConfig::crank_nicolson(dt, linear_config()).expect("config");
    let solution = march(&lowered.cycle, &config, f64::INFINITY);

    let expected_j = (3000.0 * 0.25 + 1000.0) * FULL_LOAD_DENSITY;
    assert!(
        (solution.energy.delivered_j() - expected_j).abs() < 1e-6,
        "the full dwell energy must enter the march: {:?}",
        solution.energy
    );
}

#[test]
fn scenario_dwell_energy_is_preserved_when_a_step_crosses_the_switch() {
    let lowered = lower_envelope_duty_cycle(&scenario_cycle(), "total-power", FULL_LOAD_DENSITY)
        .expect("lowers");

    for theta in [1.0, 0.5] {
        let config = TransientConfig::new(theta, 800.0, linear_config()).expect("config");
        let energy = march(&lowered.cycle, &config, f64::INFINITY).energy;
        assert!(!energy.steps_aligned());
        assert!(
            energy.residual_j().abs() < 1e-6,
            "theta={theta}: {energy:?}"
        );
    }
}

#[test]
fn an_explicit_transition_preserves_its_own_authored_energy() {
    // This physically different profile has the same window and total energy
    // as the stepped cycle. Both deliver that energy without changing either
    // declaration's transition shape.
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
