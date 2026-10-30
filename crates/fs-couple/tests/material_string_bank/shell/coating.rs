//! Full film history reaches the existing reciprocal bow and pressure owner.
//! Published O1 coefficients are a numerical fixture on authored test geometry.
use super::*;
use fs_couple::render::schedule::force::coupled::PreparedModalBudget;
use fs_couple::shear_shell::coating::relaxation::CoatingRelaxationOptions;
use fs_couple::shear_shell::{CoatingSide, ShellCoating};
use fs_material::coating::ResolvedCoating;

fn film(strain_limit: f64) -> ResolvedCoating {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, 1097.0),
        ("varnish_e0", Pressure::DIMS, 20e6),
        ("varnish_einf", Pressure::DIMS, 1.37e9),
        ("varnish_alpha", Dims::NONE, 0.391),
        ("varnish_beta", Dims::NONE, 0.360),
        ("varnish_tau", Dims([0, 0, 1, 0, 0, 0]), 0.127),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 293.15, 293.15)
                    .with(
                        "omega",
                        20.0 * core::f64::consts::TAU,
                        20000.0 * core::f64::consts::TAU,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source:
                        "numerical fixture of Laemmlein2021 O1, doi:10.1016/j.rinma.2020.100137"
                            .into(),
                    license: "Attributed numerical facts".into(),
                    artifact: None,
                },
                observations: vec![],
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "modern copal-oil O1 numerical fixture".into(),
            phase: "film".into(),
            process: "authored test geometry".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", 293.15)
        .unwrap()
        .with("omega", 1000.0 * core::f64::consts::TAU)
        .unwrap();
    let source = resolve_material_state_point(
        &card,
        &point,
        &ResolvedCoating::requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedCoating::from_material(source, 0.3, strain_limit).unwrap()
}

fn prepared(cx: &Cx<'_>, gate: &CancelGate, coated: bool) -> ReciprocalStringBank {
    prepared_with_limit(cx, gate, coated, 0.0003)
}

fn prepared_with_limit(
    cx: &Cx<'_>,
    gate: &CancelGate,
    coated: bool,
    strain_limit: f64,
) -> ReciprocalStringBank {
    let mut body = shell(700.0, 0.02, 0.17, gate);
    body.projection_budget.max_projection_terms = 64_000_000;
    if coated {
        body.specimen = body
            .specimen
            .with_coatings(
                vec![ShellCoating {
                    region: "wood".into(),
                    material: film(strain_limit),
                    side: CoatingSide::Positive,
                    thickness_m: 20e-6,
                    reference_omega_rad_s: 1000.0 * core::f64::consts::TAU,
                    max_storage_relative_drift: 0.9,
                }],
                gate,
            )
            .unwrap()
            .with_coating_relaxation(
                CoatingRelaxationOptions {
                    terms: 16,
                    storage_relative_tolerance: 0.001,
                    loss_relative_tolerance: 0.001,
                    max_memory_coordinates: 64,
                },
                gate,
            )
            .unwrap();
    }
    let mut settings = options(0);
    settings.network.max_connections = 16;
    settings.network.max_setup_terms = 131_072;
    let bank = ReciprocalStringBank::from_material_shear_shell_prepared(
        cx,
        vec![member(7800.0).source, layered().source],
        gas(101_325.0),
        body,
        None,
        Some(PreparedModalBudget {
            max_modes: 32,
            max_dense_bytes: 128 * 1_048_576,
            max_multiply_terms: 100_000_000_000,
        }),
        law(),
        settings,
        gate,
    )
    .unwrap();
    if coated {
        bank.with_material_shell_coating_relaxation(None, gate)
            .unwrap()
    } else {
        bank
    }
}

#[test]
fn g0_loaded_coating_bridge_response_retains_strings_pressure_and_branch_heat() {
    with_cx(|cx, gate| {
        let bank = prepared(cx, gate, true);
        let before: Vec<_> = bank
            .system()
            .components()
            .iter()
            .map(|component| component.states().to_vec())
            .collect();
        let history = bank.system().material_memory().to_vec();
        let omega = core::f64::consts::TAU * 1000.0;
        let peak = fs_math::c64::C64::from_re(1e-3);
        let response = bank
            .material_shell_bridge_harmonic_response(omega, 0, peak, 100_000_000, gate)
            .unwrap();
        let dynamics = &response.released.dynamics;
        assert_eq!(response.bridge_velocity_m_s.len(), 2);
        assert!(response.coating_dissipation_w > 0.0);
        assert!(response.coating_mean_stored_energy_j > 0.0);
        assert!(
            (response.coating_dissipation_w + response.other_material_dissipation_w
                - dynamics.material_dissipation_w)
                .abs()
                < 1e-12 * dynamics.material_dissipation_w.max(1e-20)
        );
        assert!(
            dynamics.power_residual_w.abs() < 1e-8 * dynamics.external_power_w.abs().max(1e-20)
        );
        let string_modes: usize = bank.system().components()[..2]
            .iter()
            .map(|component| component.modes().len())
            .sum();
        assert!(
            dynamics.velocity_m_sqrt_kg_per_s[..string_modes]
                .iter()
                .any(|v| v.abs() > 1e-15),
            "bridge forcing lost sympathetic string backreaction"
        );
        let expected_pressure = bank
            .system()
            .acceleration_pressure_weights()
            .unwrap()
            .iter()
            .zip(&dynamics.displacement_m_sqrt_kg)
            .fold(fs_math::c64::C64::ZERO, |sum, (&weight, &q)| {
                sum + q.scale(-omega * omega * weight)
            });
        assert!(
            (response.released.observer_pressure_pa - expected_pressure).abs()
                < 1e-12 * expected_pressure.abs().max(1e-20)
        );
        assert_eq!(bank.system().material_memory(), history);
        assert_eq!(bank.system().samples_rendered(), 0);
        for (component, states) in bank.system().components().iter().zip(&before) {
            assert_eq!(component.states(), states);
        }

        let multi = bank.with_additional_bows(vec![],
            fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig {
                max_ports: 4, max_sweeps: 32, max_setup_terms: 1_000_000,
            }, gate).unwrap();
        let replay = multi
            .material_shell_bridge_harmonic_response(omega, 0, peak, 100_000_000, gate)
            .unwrap();
        assert_eq!(
            response, replay,
            "multiple-bow ownership changed the released loaded response"
        );
    });
}

#[test]
fn g0_loaded_coating_response_refuses_outside_sources_or_missing_bridges() {
    with_cx(|cx, gate| {
        let bank = prepared(cx, gate, true);
        let history = bank.system().material_memory().to_vec();
        let force = fs_math::c64::C64::from_re(1e-3);
        assert!(
            bank.material_shell_bridge_harmonic_response(
                10.0 * core::f64::consts::TAU,
                0,
                force,
                100_000_000,
                gate,
            )
            .is_err()
        );
        assert!(
            bank.material_shell_bridge_harmonic_response(
                1000.0 * core::f64::consts::TAU,
                2,
                force,
                100_000_000,
                gate,
            )
            .is_err()
        );
        assert!(
            bank.material_shell_bridge_harmonic_response(
                1000.0 * core::f64::consts::TAU,
                0,
                force,
                1,
                gate,
            )
            .is_err()
        );
        assert_eq!(bank.system().material_memory(), history);
        assert_eq!(bank.system().samples_rendered(), 0);
        let residual = bank.with_residual_flexibility(gate).unwrap();
        assert!(
            residual
                .material_shell_bridge_harmonic_response(
                    1000.0 * core::f64::consts::TAU,
                    0,
                    force,
                    100_000_000,
                    gate,
                )
                .is_err(),
            "retained-only bridge force silently omitted residual flexibility"
        );
    });
}

#[test]
fn g1_coating_memory_survives_cold_rebuild_and_changes_reciprocal_pressure() {
    with_cx(|cx, gate| {
        let bank = prepared(cx, gate, true);
        let image = bank.source().material_shell_coating_relaxation().unwrap();
        assert_eq!(image.mode_count(), 3);
        let rows = image.coupling_rows_per_s().to_vec();
        assert!(!rows.is_empty());
        let history = bank.system().material_memory().to_vec();
        assert_eq!(history.len(), rows.len());
        let source = bank.source().material_shell().unwrap().specimen.coatings()[0]
            .input()
            .material
            .source()
            .clone();
        let mut coated = bank.with_residual_flexibility(gate).unwrap();
        assert_eq!(coated.system().material_memory(), history);
        assert_eq!(
            coated
                .source()
                .material_shell_coating_relaxation()
                .unwrap()
                .coupling_rows_per_s(),
            rows
        );
        assert_eq!(
            coated
                .source()
                .material_shell()
                .unwrap()
                .specimen
                .coatings()[0]
                .input()
                .material
                .source(),
            &source
        );
        let mut bare = prepared(cx, gate, false)
            .with_residual_flexibility(gate)
            .unwrap();
        let zero = vec![0.0; coated.system().mode_count()];
        let mut difference = 0.0_f64;
        let mut peak = 0.0_f64;
        for _ in 0..96 {
            let a = coated.step(&zero, drive()).unwrap();
            let b = bare.step(&zero, drive()).unwrap();
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
            difference = difference.max((a.observer_pressure_pa - b.observer_pressure_pa).abs());
            peak = peak.max(a.observer_pressure_pa.abs());
        }
        assert!(peak > 1e-12);
        assert!(
            difference > 1e-12,
            "film did not change the shared pressure"
        );
        assert!(
            coated
                .system()
                .material_memory()
                .iter()
                .any(|z| z.abs() > 1e-14)
        );
    });
}

#[test]
fn g0_transient_coating_observation_reads_original_histories_without_advancing() {
    with_cx(|cx, gate| {
        let bare = prepared(cx, gate, false);
        assert_eq!(bare.material_shell_coating_observation_work(gate).unwrap(), None);
        assert_eq!(bare.material_shell_coating_observation(0, gate).unwrap(), None);

        let mut bank = prepared(cx, gate, true).with_residual_flexibility(gate).unwrap();
        let work = bank.material_shell_coating_observation_work(gate).unwrap().unwrap();
        let initial = bank.material_shell_coating_observation(work, gate).unwrap().unwrap();
        assert_eq!(initial.branch_count, bank.source().material_shell_coating_relaxation()
            .unwrap().relaxation_times_s().len());
        assert_eq!(initial.stored_energy_j, 0.0);
        assert_eq!(initial.dissipated_power_w, 0.0);
        assert!(bank.material_shell_coating_observation(work - 1, gate).is_err());
        let zero = vec![0.0; bank.system().mode_count()];
        for _ in 0..32 { bank.step(&zero, drive()).unwrap(); }
        let history = bank.system().material_memory().to_vec();
        let observed = bank.material_shell_coating_observation(work, gate).unwrap().unwrap();
        assert!(observed.stored_energy_j > 0.0);
        assert!(observed.dissipated_power_w > 0.0);
        // This fixture's only material histories are coating arms. Their
        // storage must equal the original owner's complete material account.
        assert_eq!(observed.stored_energy_j, bank.system().material_energy_j().unwrap());
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(bank.material_shell_coating_observation(work, &cancelled).is_err());
        assert_eq!(bank.system().material_memory(), history);
        assert_eq!(bank.system().samples_rendered(), 32);
    });
}

use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::{
    MultiBowPerformance, MultiBowPerformanceConfig,
    coating::{ShellCoatingMonitorError, ShellCoatingMonitorOptions},
    force::{StringForceCoupling, StringForcePort, StringForceProgram},
};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::schedule::force::{ForceInitialization, ForceRenderConfig};
use std::sync::Arc;

fn monitor_options() -> ShellCoatingMonitorOptions {
    ShellCoatingMonitorOptions {
        maximum_bytes: 32 * 1_048_576,
        maximum_work: usize::MAX / 8,
    }
}

fn monitor_performance(cx: &Cx<'_>, gate: &CancelGate, limit: f64) -> MultiBowPerformance {
    let bank = prepared_with_limit(cx, gate, true, limit).with_additional_bows(vec![],
        fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig {
            max_ports: 4, max_sweeps: 32, max_setup_terms: 1_000_000,
        }, gate).unwrap();
    let schedule = GestureSchedule::try_new(
        700,
        vec![GestureTrack {
            id: "bow".into(),
            target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow {
                velocity_m_per_s: 0.01,
                normal_force_n: 0.1,
                station: 0.2,
            },
            events: vec![],
        }],
    )
    .unwrap();
    bank.into_scheduled(
        schedule,
        vec!["bow".into()],
        MultiBowPerformanceConfig {
            samples: 32,
            max_block: 16,
            max_work: 4096,
            max_events: 64,
        },
    )
    .unwrap()
}

fn monitor_preload() -> StringForceProgram {
    StringForceProgram {
        ports: vec![StringForcePort {
            location: StringBowLocation {
                string: 0,
                station_fraction: 0.31,
            },
            coupling: StringForceCoupling::Centreline,
        }],
        initial_force_n: vec![0.001],
        events: vec![],
        initialization: ForceInitialization::StaticPreload,
        limits: ForceRenderConfig {
            sample_rate_hz: 48_000,
            max_block: 16,
            max_events: 64,
            max_controls: 1024,
            max_projection_terms: 32_768,
        },
    }
}

#[test]
fn g0_coating_monitor_checks_every_endpoint_without_changing_sound_or_mechanics() {
    with_cx(|cx, gate| {
        let mut plain = monitor_performance(cx, gate, 0.0003);
        let mut guarded = monitor_performance(cx, gate, 0.0003)
            .with_material_shell_coating_monitor(
                monitor_options(),
                Arc::new(CancelGate::new_clock_free()),
            )
            .unwrap();
        let report = guarded.coating_monitor_report().unwrap();
        assert_eq!(report.checked_endpoints, 1);
        assert_eq!(report.resources.evaluations, 33);
        assert_eq!(report.last_checked_sample, 0);
        assert_eq!(report.films().count(), 1);
        let memory = report.relaxation.unwrap();
        assert!(memory.current.branch_count > 0);
        assert_eq!(memory.current.stored_energy_j, 0.0);
        assert_eq!(memory.current.dissipated_power_w, 0.0);
        assert_eq!(memory.maximum_stored_energy_sample, 0);
        assert_eq!(memory.maximum_dissipated_power_sample, 0);
        assert_eq!(
            report.films().next().unwrap().maximum.frobenius_upper_bound,
            0.0
        );
        let mut a = [0.0; 32];
        let mut b = [0.0; 32];
        for block in a.chunks_mut(16) {
            plain.block(block).unwrap();
        }
        for block in b.chunks_mut(7) {
            guarded.block(block).unwrap();
        }
        assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
        assert!(a.iter().any(|pressure| pressure.abs() > 1e-12));
        assert_eq!(
            plain.bank().system().material_memory(),
            guarded.bank().system().material_memory()
        );
        for (a, b) in plain
            .bank()
            .system()
            .components()
            .iter()
            .zip(guarded.bank().system().components())
        {
            assert_eq!(a.states(), b.states());
        }
        let report = guarded.coating_monitor_report().unwrap();
        assert_eq!(report.checked_endpoints, 33);
        assert_eq!(report.last_checked_sample, 32);
        assert!(report.failure.is_none());
        let memory = report.relaxation.unwrap();
        let work = guarded.bank().material_shell_coating_observation_work(gate)
            .unwrap().unwrap();
        assert_eq!(memory.current, guarded.bank()
            .material_shell_coating_observation(work, gate).unwrap().unwrap());
        assert!(memory.current.stored_energy_j > 0.0);
        assert!(memory.current.dissipated_power_w > 0.0);
        assert!(memory.maximum_stored_energy_j >= memory.current.stored_energy_j);
        assert!(memory.maximum_dissipated_power_w >= memory.current.dissipated_power_w);
        assert!((1..=32).contains(&memory.maximum_stored_energy_sample));
        assert!((1..=32).contains(&memory.maximum_dissipated_power_sample));
        let film = report.films().next().unwrap();
        assert!(film.maximum.frobenius_upper_bound > 0.0);
        assert!(film.maximum.limit_ratio < 1.0);
        assert_eq!(
            film.maximum.covered_elements,
            guarded
                .bank()
                .source()
                .bank()
                .material_shell_image()
                .unwrap()
                .source()
                .mesh()
                .tris
                .len()
        );
        assert!(std::ptr::eq(
            film.source.source,
            &guarded
                .bank()
                .source()
                .bank()
                .material_shell_image()
                .unwrap()
                .source()
                .coatings()[0]
        ));
    });
}

#[test]
fn g4_coating_monitor_refuses_initial_and_accepted_over_envelope_states() {
    with_cx(|cx, gate| {
        let loaded = monitor_performance(cx, gate, 1e-30)
            .with_string_forces(monitor_preload(), gate)
            .unwrap();
        let error = match loaded.with_material_shell_coating_monitor(
            monitor_options(),
            Arc::new(CancelGate::new_clock_free()),
        ) {
            Ok(_) => panic!("preloaded film envelope was not checked at sample zero"),
            Err(error) => error,
        };
        assert!(
            matches!(error, ShellCoatingMonitorError::Envelope { sample: 0, ratio, .. } if ratio > 1.0)
        );

        let mut performance = monitor_performance(cx, gate, 1e-30)
            .with_material_shell_coating_monitor(
                monitor_options(),
                Arc::new(CancelGate::new_clock_free()),
            )
            .unwrap();
        assert!(performance.block(&mut [42.0; 16]).is_err());
        let clock = performance.samples_rendered();
        assert!((1..=16).contains(&clock));
        let report = performance.coating_monitor_report().unwrap();
        assert_eq!(report.last_checked_sample, clock);
        assert_eq!(report.checked_endpoints as u64, clock + 1);
        assert!(report.relaxation.unwrap().current.stored_energy_j > 0.0);
        assert!(matches!(report.failure,
            Some(ShellCoatingMonitorError::Envelope { sample, ratio, .. }) if *sample == clock && *ratio > 1.0));
        assert!(performance.validate_sample_count(0).is_err());
        assert!(performance.block(&mut [0.0]).is_err());
        assert_eq!(
            performance.samples_rendered(),
            clock,
            "refused output was resumed"
        );
    });
}

#[test]
fn g4_coating_monitor_budgets_cancellation_and_late_preload_are_explicit() {
    with_cx(|cx, gate| {
        let original = monitor_performance(cx, gate, 0.0003)
            .with_material_shell_coating_monitor(
                monitor_options(),
                Arc::new(CancelGate::new_clock_free()),
            )
            .unwrap();
        let resources = *original.coating_monitor_report().unwrap().resources;
        let exact = ShellCoatingMonitorOptions {
            maximum_bytes: resources.retained_bytes + resources.evaluation_bytes,
            maximum_work: resources.total_work,
        };
        let cancel = Arc::new(CancelGate::new_clock_free());
        let mut performance = monitor_performance(cx, gate, 0.0003)
            .with_material_shell_coating_monitor(exact, cancel.clone())
            .unwrap();
        for short in [
            ShellCoatingMonitorOptions {
                maximum_bytes: exact.maximum_bytes - 1,
                ..exact
            },
            ShellCoatingMonitorOptions {
                maximum_work: exact.maximum_work - 1,
                ..exact
            },
        ] {
            assert!(
                monitor_performance(cx, gate, 0.0003)
                    .with_material_shell_coating_monitor(
                        short,
                        Arc::new(CancelGate::new_clock_free())
                    )
                    .is_err()
            );
        }
        performance.block(&mut [0.0]).unwrap();
        cancel.request();
        assert!(performance.block(&mut [0.0]).is_err());
        assert_eq!(performance.samples_rendered(), 1);
        assert_eq!(
            performance
                .coating_monitor_report()
                .unwrap()
                .checked_endpoints,
            2
        );

        let mut changed = monitor_performance(cx, gate, 0.0003)
            .with_material_shell_coating_monitor(
                monitor_options(),
                Arc::new(CancelGate::new_clock_free()),
            )
            .unwrap()
            .with_string_forces(monitor_preload(), gate)
            .unwrap();
        assert!(
            changed.block(&mut [0.0]).is_err(),
            "late preload bypassed the initial observation"
        );
        assert_eq!(changed.samples_rendered(), 0);
    });
}
