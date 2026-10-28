//! G0/G3 source transport and refusal checks, not transient calibration.
use super::*;
use fs_matdb::EvaluationDecision;
use fs_qty::semantic::FrequencyConvention;

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);
const ACCELERATION: Dims = Dims([1, 0, -2, 0, 0, 0]);
const VISCOSITY: Dims = Dims([-1, 1, -1, 0, 0, 0]);

fn geometry(point: QueryPoint, height: f64) -> QueryPoint {
    point
        .with_quantity(
            "source-original-height",
            QuantitySpec::dimensional(LENGTH),
            height,
        )
        .unwrap()
}

fn deng_reference_point() -> QueryPoint {
    geometry(
        flags(&[
            ("source-sample-hongda-pu65-deng2025", 1.0),
            ("source-test-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
            ("source-fitted-spectrum-reference", 1.0),
            ("source-time-domain-validation-known", 0.0),
            ("source-precompression", 0.25),
            ("source-preconditioning-section2-2", 1.0),
        ]),
        0.08,
    )
    .with_quantity(
        "source-original-area",
        QuantitySpec::dimensional(AREA),
        0.0064,
    )
    .unwrap()
    .with_quantity(
        "source-acceleration-peak-to-peak",
        QuantitySpec::dimensional(ACCELERATION),
        0.25,
    )
    .unwrap()
    .with_quantity(
        "source-frequency",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Frequency(FrequencyConvention::Cyclic),
            ValueForm::Static,
        )),
        1.0,
    )
    .unwrap()
}

#[test]
fn g0_g3_deng_foam_preserves_axial_spectrum_and_original_geometry() {
    let pack = compile("hongda-pu65-deng2025-prony-reference", 13);
    let p = deng_reference_point();
    sample(&pack, "density", RHO, &p, 65.0);
    sample(&pack, "prony-branch-count", Dims::NONE, &p, 5.0);
    sample(&pack, "prony-equilibrium-modulus", PA, &p, 3050.0);
    for (branch, modulus, seconds) in [
        (1, 23900.0, 0.01),
        (2, 19700.0, 0.1),
        (3, 8990.0, 1.0),
        (4, 12400.0, 10.0),
        (5, 8140.0, 100.0),
    ] {
        let answer = sample(&pack, &format!("prony-modulus-{branch}"), PA, &p, modulus);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        sample(
            &pack,
            &format!("prony-relaxation-time-{branch}"),
            TIME,
            &p,
            seconds,
        );
    }
    for bad in [
        QueryPoint::new(),
        p.clone().with("source-original-height", 0.06).unwrap(),
        p.clone().with("source-original-area", 0.0048).unwrap(),
        p.clone().with("source-precompression", 0.0).unwrap(),
        p.clone().with("source-frequency", 20.01).unwrap(),
        p.clone()
            .with("source-fitted-spectrum-reference", 0.0)
            .unwrap(),
        p.clone()
            .with("source-time-domain-validation-known", 1.0)
            .unwrap(),
        p.clone()
            .with("source-test-temperature-known", 1.0)
            .unwrap(),
        p.clone()
            .with("source-preconditioning-section2-2", 0.0)
            .unwrap(),
    ] {
        refuses(&pack, "prony-equilibrium-modulus", PA, &bad);
    }
    refuses(&pack, "prony-modulus-1", VISCOSITY, &p);
    // The paper's coupled finite-strain model does not turn this axial
    // frequency fit into shear coefficients or a raw crushing curve.
    refuses(&pack, "shear-modulus", PA, &p);
    refuses(&pack, "compressive-stress", PA, &p);
    refuses(&pack, "specific-heat", CP, &p);
}

fn tpu_geometry(point: QueryPoint) -> QueryPoint {
    geometry(point, 0.02)
        .with_quantity(
            "source-original-diameter",
            QuantitySpec::dimensional(LENGTH),
            0.02,
        )
        .unwrap()
}

fn tpu_reference_point() -> QueryPoint {
    tpu_geometry(flags(&[
        ("source-sample-rec-tpu-easyflex-golovkin2026", 1.0),
        ("source-density-known", 0.0),
        ("source-infill-fraction", 1.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-fitted-spectrum-reference", 1.0),
        ("source-fit-method-creep", 1.0),
        ("source-canonical-parallel-maxwell-interpretation", 1.0),
        ("source-time-domain-validation-known", 0.0),
        ("source-time-domain-experiments-reported", 1.0),
        ("source-loading-level-definition-ambiguous", 1.0),
        ("source-test-temperature-known", 0.0),
    ]))
    .with_quantity(
        "source-test-duration",
        QuantitySpec::dimensional(TIME),
        6000.0,
    )
    .unwrap()
}

fn tpu_cycle_point(speed: f64) -> QueryPoint {
    temperature(
        tpu_geometry(flags(&[
            ("source-sample-rec-tpu-easyflex-golovkin2026", 1.0),
            ("source-density-known", 0.0),
            ("source-infill-fraction", 1.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
            ("source-fitted-spectrum-reference", 0.0),
            ("source-measured-curve-summary", 1.0),
            ("source-test-cycle", 1.0),
            ("source-peak-compressive-strain", 0.15),
            ("source-test-replicates", 3.0),
        ])),
        "temperature",
        296.15,
    )
    .with_quantity("test-speed", QuantitySpec::dimensional(SPEED), speed)
    .unwrap()
}

#[test]
fn g0_g3_tpu_creep_fit_and_measured_cycle_work_remain_separate() {
    let pack = compile("rec-tpu-easyflex-golovkin2026-creep-reference", 11);
    let p = tpu_reference_point();
    sample(&pack, "prony-branch-count", Dims::NONE, &p, 2.0);
    sample(&pack, "prony-equilibrium-modulus", PA, &p, 25.9e6);
    for (branch, modulus, viscosity) in [(1, 7.3e6, 1291e6), (2, 18.8e6, 53.9e6)] {
        sample(&pack, &format!("prony-modulus-{branch}"), PA, &p, modulus);
        sample(
            &pack,
            &format!("prony-branch-viscosity-{branch}"),
            VISCOSITY,
            &p,
            viscosity,
        );
        let tau = sample(
            &pack,
            &format!("prony-relaxation-time-{branch}"),
            TIME,
            &p,
            viscosity / modulus,
        );
        assert_eq!(tau.evidence.value.uncertainty, UncertaintyModel::Unstated);
    }
    refuses(&pack, "density", RHO, &p);
    refuses(&pack, "prony-modulus-3", PA, &p);
    for bad in [
        QueryPoint::new(),
        p.clone().with("source-original-height", 0.01).unwrap(),
        p.clone().with("source-infill-fraction", 0.5).unwrap(),
        p.clone().with("source-fit-method-creep", 0.0).unwrap(),
        p.clone()
            .with("source-time-domain-validation-known", 1.0)
            .unwrap(),
        p.clone()
            .with("source-loading-level-definition-ambiguous", 0.0)
            .unwrap(),
        p.clone()
            .with("source-canonical-parallel-maxwell-interpretation", 0.0)
            .unwrap(),
    ] {
        refuses(&pack, "prony-modulus-1", PA, &bad);
    }
    for (speed, work, tangent, sd) in [
        (0.000002, 119000.0, 33.4e6, 1.3e6),
        (0.00002, 134000.0, 37.5e6, 1.7e6),
        (0.0002, 159000.0, 41.7e6, 1.7e6),
        (0.001, 194000.0, 45.4e6, 0.6e6),
        (0.002, 218000.0, 49.6e6, 1.7e6),
    ] {
        let measured = tpu_cycle_point(speed);
        for (name, expected) in [
            ("first-cycle-work-density", work),
            ("compression-tangent-modulus", tangent),
            (
                "compression-tangent-modulus-reported-standard-deviation",
                sd,
            ),
        ] {
            let answer = sample(&pack, name, PA, &measured, expected);
            assert!(matches!(
                answer.receipt.decision,
                EvaluationDecision::ExactTabulated { .. }
            ));
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
        // A measured tangent is not another fitted branch.
        refuses(&pack, "prony-modulus-1", PA, &measured);
    }
    let measured = tpu_cycle_point(0.0002);
    for bad in [
        p,
        measured.clone().with("test-speed", 0.0003).unwrap(),
        measured.clone().with("test-speed", 0.003).unwrap(),
        measured.clone().with("source-test-cycle", 2.0).unwrap(),
        measured
            .clone()
            .with("source-peak-compressive-strain", 0.10)
            .unwrap(),
        measured.clone().with("temperature", 298.15).unwrap(),
        measured
            .clone()
            .with("source-original-diameter", 0.01)
            .unwrap(),
        measured
            .clone()
            .with("source-fitted-spectrum-reference", 1.0)
            .unwrap(),
    ] {
        refuses(&pack, "first-cycle-work-density", PA, &bad);
    }
}

fn prony_card(pack: &NormalizedPack) -> fs_matdb::NormalizedMaterialCardPack {
    fs_matdb::NormalizedMaterialCardPack::new(
        fs_matdb::MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "source specimen solid".into(),
            process:
                "source-fitted incremental Prony reference; no independent transient validation"
                    .into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap()
}

fn close_history(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 + 2e-9 * expected.abs(),
        "history response {actual} != independent reference {expected}"
    );
}

#[test]
fn g1_g3_source_foam_drives_relaxation_creep_and_zero_force_recovery() {
    use fs_material::visco::source::{
        PronyReferenceLimits, PronyReferenceRequest, ResolvedPronyReference,
    };
    let pack = compile("hongda-pu65-deng2025-prony-reference", 13);
    let card = prony_card(&pack);
    let reference = ResolvedPronyReference::resolve(
        card.card(),
        &deng_reference_point(),
        &PronyReferenceRequest {
            // Authored bounds for this numerical exercise, not measured limits.
            limits: PronyReferenceLimits {
                minimum_incremental_strain: -0.1,
                maximum_incremental_strain: 0.1,
                maximum_absolute_strain_rate_s: 1.0,
                maximum_elapsed_s: 2000.0,
            },
            original_length_axis: Some("source-original-height".into()),
            original_area_axis: Some("source-original-area".into()),
        },
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(reference.properties().len(), 13);
    assert_eq!(reference.card_identity(), card.card().content_hash());
    assert!(reference.coupon(0.06, 0.0064).is_err());
    let volume = 0.08 * 0.0064;
    let mut relaxed = reference.coupon(0.08, 0.0064).unwrap();
    close_history(relaxed.mass_kg().unwrap(), 65.0 * volume);
    let loaded = relaxed.advance_strain(0.01, 1.0).unwrap();
    let held = relaxed.advance_strain(0.01, 10.0).unwrap();
    // Independent analytic solution for a 1 s ramp followed by a 10 s hold,
    // using the published coefficients and the standard library exponential.
    let terms: [(f64, f64); 5] = [
        (23900.0, 0.01),
        (19700.0, 0.1),
        (8990.0, 1.0),
        (12400.0, 10.0),
        (8140.0, 100.0),
    ];
    let mut expected_stress = 3050.0 * 0.01;
    let mut expected_storage = 0.5 * 3050.0 * 0.01_f64.powi(2);
    for (e, tau) in terms {
        let q = -0.01 * tau * (-1.0 / tau).exp_m1() * (-10.0 / tau).exp();
        expected_stress += e * q;
        expected_storage += 0.5 * e * q * q;
    }
    close_history(held.incremental_stress_pa, expected_stress);
    close_history(held.incremental_force_n, expected_stress * 0.0064);
    close_history(held.stored_energy_j, expected_storage * volume);
    assert!(held.incremental_stress_pa < loaded.incremental_stress_pa);
    assert_eq!(held.incremental_model_work_j, 0.0);
    close_history(
        held.incremental_dissipated_energy_j,
        loaded.stored_energy_j - held.stored_energy_j,
    );

    let mut creep = reference.coupon(0.08, 0.0064).unwrap();
    let loaded = creep.advance_force(0.32, 1.0).unwrap();
    let held = creep.advance_force(0.32, 100.0).unwrap();
    assert!(held.incremental_strain > loaded.incremental_strain);
    close_history(held.incremental_force_n, 0.32);
    close_history(
        held.incremental_model_work_j,
        0.32 * (held.incremental_displacement_m - loaded.incremental_displacement_m),
    );
    let unloaded = creep.advance_force(0.0, 1.0).unwrap();
    let recovered = creep.advance_force(0.0, 1000.0).unwrap();
    assert!(unloaded.incremental_strain > 0.0);
    assert!(unloaded.incremental_strain < held.incremental_strain);
    assert!(recovered.incremental_strain > 0.0);
    assert!(recovered.incremental_strain < unloaded.incremental_strain);
    close_history(recovered.incremental_force_n, 0.0);
    close_history(recovered.incremental_model_work_j, 0.0);
    close_history(
        recovered.incremental_dissipated_energy_j,
        unloaded.stored_energy_j - recovered.stored_energy_j,
    );
    for response in [loaded, held, unloaded, recovered] {
        assert!(response.incremental_dissipated_energy_j >= 0.0);
        close_history(
            response.model_work_j,
            response.stored_energy_j + response.dissipated_energy_j,
        );
        close_history(response.energy_balance_residual_j, 0.0);
    }
}

#[test]
fn g3_tpu_triangular_reference_exposes_measured_cycle_work_mismatch() {
    use fs_material::visco::source::{
        PronyReferenceLimits, PronyReferenceRequest, ResolvedPronyReference,
    };
    let pack = compile("rec-tpu-easyflex-golovkin2026-creep-reference", 11);
    let card = prony_card(&pack);
    let reference = ResolvedPronyReference::resolve(
        card.card(),
        &tpu_reference_point(),
        &PronyReferenceRequest {
            limits: PronyReferenceLimits {
                minimum_incremental_strain: 0.0,
                maximum_incremental_strain: 0.15,
                maximum_absolute_strain_rate_s: 0.02,
                maximum_elapsed_s: 30.0,
            },
            original_length_axis: Some("source-original-height".into()),
            original_area_axis: None,
        },
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    let area = std::f64::consts::PI * 0.01_f64.powi(2);
    let volume = area * 0.02;
    let mut coupon = reference.coupon(0.02, area).unwrap();
    assert_eq!(coupon.mass_kg(), None);
    // 12 mm/min = 0.0002 m/s; 15% of 20 mm takes 15 s per ramp.
    coupon.advance_strain(0.15, 15.0).unwrap();
    let cycle = coupon.advance_strain(0.0, 15.0).unwrap();
    let measured = sample(
        &pack,
        "first-cycle-work-density",
        PA,
        &tpu_cycle_point(0.0002),
        159000.0,
    )
    .evidence
    .value
    .value;
    let predicted = cycle.model_work_j / volume;
    assert!((predicted - 124389.558).abs() < 0.001);
    let relative_error = (predicted - measured).abs() / measured;
    // This is an observed mismatch, deliberately not a physical-accuracy pass:
    // the paper's reconstructed stress-driven paths are a different experiment.
    assert!(relative_error > 0.20);
    assert!(cycle.stored_energy_j > 0.0);
    assert!(cycle.dissipated_energy_j > 0.0);
    close_history(
        cycle.model_work_j,
        cycle.stored_energy_j + cycle.dissipated_energy_j,
    );
    close_history(cycle.energy_balance_residual_j, 0.0);
}
