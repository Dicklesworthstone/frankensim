//! Full film history reaches the existing reciprocal bow and pressure owner.
//! Published O1 coefficients are a numerical fixture on authored test geometry.
use super::*;
use fs_couple::render::schedule::force::coupled::PreparedModalBudget;
use fs_couple::shear_shell::coating::relaxation::CoatingRelaxationOptions;
use fs_couple::shear_shell::{CoatingSide, ShellCoating};
use fs_material::coating::ResolvedCoating;

fn film() -> ResolvedCoating {
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
    ResolvedCoating::from_material(source, 0.3, 0.0003).unwrap()
}

fn prepared(cx: &Cx<'_>, gate: &CancelGate, coated: bool) -> ReciprocalStringBank {
    let mut body = shell(700.0, 0.02, 0.17, gate);
    body.projection_budget.max_projection_terms = 64_000_000;
    if coated {
        body.specimen = body
            .specimen
            .with_coatings(
                vec![ShellCoating {
                    region: "wood".into(),
                    material: film(),
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
        let before: Vec<_> = bank.system().components().iter()
            .map(|component| component.states().to_vec()).collect();
        let history = bank.system().material_memory().to_vec();
        let omega = core::f64::consts::TAU * 1000.0;
        let peak = fs_math::c64::C64::from_re(1e-3);
        let response = bank.material_shell_bridge_harmonic_response(
            omega, 0, peak, 100_000_000, gate,
        ).unwrap();
        let dynamics = &response.released.dynamics;
        assert_eq!(response.bridge_velocity_m_s.len(), 2);
        assert!(response.coating_dissipation_w > 0.0);
        assert!(response.coating_mean_stored_energy_j > 0.0);
        assert!((response.coating_dissipation_w + response.other_material_dissipation_w
            - dynamics.material_dissipation_w).abs()
            < 1e-12 * dynamics.material_dissipation_w.max(1e-20));
        assert!(dynamics.power_residual_w.abs() < 1e-8 * dynamics.external_power_w.abs().max(1e-20));
        let string_modes: usize = bank.system().components()[..2].iter()
            .map(|component| component.modes().len()).sum();
        assert!(dynamics.velocity_m_sqrt_kg_per_s[..string_modes].iter().any(|v| v.abs() > 1e-15),
            "bridge forcing lost sympathetic string backreaction");
        let expected_pressure = bank.system().acceleration_pressure_weights().unwrap().iter()
            .zip(&dynamics.displacement_m_sqrt_kg)
            .fold(fs_math::c64::C64::ZERO, |sum, (&weight, &q)| sum + q.scale(-omega * omega * weight));
        assert!((response.released.observer_pressure_pa - expected_pressure).abs()
            < 1e-12 * expected_pressure.abs().max(1e-20));
        assert_eq!(bank.system().material_memory(), history);
        assert_eq!(bank.system().samples_rendered(), 0);
        for (component, states) in bank.system().components().iter().zip(&before) {
            assert_eq!(component.states(), states);
        }

        let multi = bank.with_additional_bows(vec![],
            fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig {
                max_ports: 4, max_sweeps: 32, max_setup_terms: 1_000_000,
            }, gate).unwrap();
        let replay = multi.material_shell_bridge_harmonic_response(
            omega, 0, peak, 100_000_000, gate,
        ).unwrap();
        assert_eq!(response, replay, "multiple-bow ownership changed the released loaded response");
    });
}

#[test]
fn g0_loaded_coating_response_refuses_outside_sources_or_missing_bridges() {
    with_cx(|cx, gate| {
        let bank = prepared(cx, gate, true);
        let history = bank.system().material_memory().to_vec();
        let force = fs_math::c64::C64::from_re(1e-3);
        assert!(bank.material_shell_bridge_harmonic_response(
            10.0 * core::f64::consts::TAU, 0, force, 100_000_000, gate,
        ).is_err());
        assert!(bank.material_shell_bridge_harmonic_response(
            1000.0 * core::f64::consts::TAU, 2, force, 100_000_000, gate,
        ).is_err());
        assert!(bank.material_shell_bridge_harmonic_response(
            1000.0 * core::f64::consts::TAU, 0, force, 1, gate,
        ).is_err());
        assert_eq!(bank.system().material_memory(), history);
        assert_eq!(bank.system().samples_rendered(), 0);
        let residual = bank.with_residual_flexibility(gate).unwrap();
        assert!(residual.material_shell_bridge_harmonic_response(
            1000.0 * core::f64::consts::TAU, 0, force, 100_000_000, gate,
        ).is_err(), "retained-only bridge force silently omitted residual flexibility");
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
