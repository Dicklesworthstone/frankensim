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
