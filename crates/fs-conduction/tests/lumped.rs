//! The lumped reduced rung: its ODE, its validity gate, and — the test that
//! justifies the whole bead — a demonstration that the gate refuses exactly
//! the regime where the cheap model is wrong (bead
//! `frankensim-extreal-program-f85xj.5.13`).

mod support;

use fs_blake3::ContentHash;
use fs_conduction::ConductionError;
use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::box_grid;
use fs_conduction::lumped::{
    BiotGate, LUMPED_BIOT_CEILING, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig, LumpedNetwork,
    LumpedNode, LumpedThermalTransport, ValidityVerdict, extract_node_from_steady_rise,
    solve_gated, solve_lumped_enthalpy,
};
use fs_conduction::material::{CONDUCTIVITY_DIMS, ConductivityModel};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::radiation::SURFACE_EMISSIVITY_PROPERTY;
use fs_conduction::solve::LinearConfig;
use fs_conduction::transient::{TransientConfig, TransientProblem, VolumetricHeatCapacity, march};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use support::with_cx;

const K: f64 = 10.0;
const RHO_CP: f64 = 2.0e6;
const AMBIENT: f64 = 300.0;
const EXCESS: f64 = 100.0;
// Unit cube: V = 1 m^3, A = 6 m^2, so Lc = V/A = 1/6 and Bi = h/(6k).
const AREA: f64 = 6.0;
const CHAR_LENGTH: f64 = 1.0 / 6.0;

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

fn htc_for(biot: f64) -> f64 {
    biot * K / CHAR_LENGTH
}

/// The reduced model of the unit cube at a given Biot number.
fn node_at(biot: f64) -> LumpedNode {
    LumpedNode::new(
        "cube",
        RHO_CP, // C = rho c_p V, V = 1
        htc_for(biot) * AREA,
        CHAR_LENGTH,
        K,
        AREA,
    )
    .expect("node admits")
}

fn network_at(biot: f64) -> LumpedNetwork {
    LumpedNetwork::new(vec![node_at(biot)], AMBIENT).expect("network admits")
}

fn phase_curve() -> EquilibriumEnthalpyPhaseCurve {
    phase_curve_for(ContentHash([0x71; 32]))
}

fn phase_curve_for(material_card_identity: ContentHash) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        material_card_identity,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 11_300.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 30_000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 11_100.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 55_000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 10_600.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 95_000.0,
                temperature_k: 800.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 10_300.0,
            },
        ],
    )
    .expect("admitted synthetic phase curve")
}

fn thermal_material_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("thermal_conductivity", CONDUCTIVITY_DIMS, 35.0),
        (SURFACE_EMISSIVITY_PROPERTY, Dims::NONE, 0.5),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 800.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("synthetic card-backed {name}"),
                    license: "CC0-1.0".to_owned(),
                    artifact: None,
                },
            })
            .expect("thermal property claim");
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic phase-change metal".to_owned(),
            phase: "solid-liquid".to_owned(),
            process: "test".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("thermal material card")
}

fn enthalpy_body<'a>(curve: &'a EquilibriumEnthalpyPhaseCurve) -> LumpedEnthalpyBody<'a> {
    LumpedEnthalpyBody::try_new("phase body", 2.0, 0.01, 0.0, 0.0, 0.001, 35.0, curve)
        .expect("admitted phase body")
}

fn enthalpy_config() -> LumpedEnthalpyMarchConfig {
    LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: 20_000.0,
        ambient_temperature_k: 300.0,
        internal_power_w: 500.0,
        duration_s: 100.0,
        maximum_step_s: 1.0,
        maximum_steps: 100,
        enthalpy_tolerance_j_kg: 1.0e-8,
    }
}

/// Mean temperature of the FULL transient at one lumped time constant.
fn full_rung_mean_at_time_constant(biot: f64) -> f64 {
    let mesh = unit_mesh(3);
    let htc = htc_for(biot);
    let tau = RHO_CP / (htc * AREA);
    let boundary: ThermalBoundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "ambient",
            |_face| true,
            ThermalBc::robin(htc, AMBIENT).expect("robin"),
        )
        .expect("region")
        .adiabatic_remainder()
        .finish()
        .expect("partition");
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::uniform("volumetric source", 0.0).expect("no source");
    let capacity = VolumetricHeatCapacity::declared(RHO_CP).expect("capacity");
    let steps = 400;
    let config =
        TransientConfig::crank_nicolson(tau / f64::from(steps), linear_config()).expect("config");
    let initial = vec![AMBIENT + EXCESS; mesh.vertex_count()];

    let solution = with_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
                capacity,
            },
            &config,
            &initial,
            steps as usize,
        )
        .expect("full march")
    });
    volume_weighted_mean(&mesh, &solution.temperature)
}

/// Volume-weighted mean temperature, using the lumped nodal volumes
/// `w_a = (1/4) Σ_{e ∋ a} V_e`.
///
/// A plain VERTEX average would be wrong here and wrong in a way that looks
/// plausible: a coarse cube mesh carries most of its vertices on the boundary
/// (56 of 64 at n = 3), so an unweighted mean is dominated by the coldest
/// nodes and reports a body far cooler than it is. The lumped model's state
/// is a volumetric mean, so that is what it must be compared against.
fn volume_weighted_mean(mesh: &ConductionMesh, temperature: &[f64]) -> f64 {
    let mut weights = vec![0.0f64; mesh.vertex_count()];
    for element in 0..mesh.element_count() {
        let quarter = mesh.element_volume(element) / 4.0;
        for vertex in mesh.complex().tets[element] {
            weights[vertex as usize] += quarter;
        }
    }
    let total: f64 = weights.iter().sum();
    weights
        .iter()
        .zip(temperature.iter())
        .fold(0.0f64, |acc, (w, t)| w.mul_add(*t, acc))
        / total
}

// ---------------------------------------------------------------------------
// The reduced model's own behaviour.
// ---------------------------------------------------------------------------

#[test]
fn the_closed_form_response_satisfies_the_ode_it_claims_to_solve() {
    // Checking the analytic response against the exponential formula would be
    // circular — it IS that formula. The independent check is the DIFFERENTIAL
    // EQUATION: C dT/dt must equal P - hA (T - T_amb) at every instant.
    let network = network_at(0.05);
    let node = &network.nodes()[0];
    let power = [500.0f64];
    let initial = [AMBIENT + EXCESS];

    for t in [0.0f64, 0.25, 1.0, 3.0].map(|f| f * node.time_constant_s()) {
        let h = node.time_constant_s() * 1e-6;
        let now = network.response_at(&power, &initial, t).expect("t")[0];
        let later = network.response_at(&power, &initial, t + h).expect("t+h")[0];
        let numeric_rate = (later - now) / h;
        let ode_rate =
            (power[0] - node.conductance_w_per_k() * (now - AMBIENT)) / node.capacitance_j_per_k();
        assert!(
            (numeric_rate - ode_rate).abs() < 1e-6 * ode_rate.abs().max(1.0),
            "at t={t}: dT/dt {numeric_rate} does not satisfy the ODE ({ode_rate})"
        );
    }
}

#[test]
fn the_response_honours_its_endpoints_and_time_constant() {
    let network = network_at(0.05);
    let node = &network.nodes()[0];
    let power = [0.0f64];
    let initial = [AMBIENT + EXCESS];

    let at_zero = network.response_at(&power, &initial, 0.0).expect("t=0")[0];
    assert!((at_zero - (AMBIENT + EXCESS)).abs() < 1e-12);

    // One time constant is 1/e of the way, by definition.
    let at_tau = network
        .response_at(&power, &initial, node.time_constant_s())
        .expect("tau")[0];
    let expected = AMBIENT + EXCESS * (-1.0f64).exp();
    assert!((at_tau - expected).abs() < 1e-9, "{at_tau} != {expected}");

    // Far future settles on the steady state.
    let late = network
        .response_at(&power, &initial, node.time_constant_s() * 40.0)
        .expect("late")[0];
    assert!((late - AMBIENT).abs() < 1e-6);
}

#[test]
fn steady_temperature_is_ambient_plus_power_over_conductance() {
    let network = network_at(0.05);
    let node = &network.nodes()[0];
    let power = 750.0;
    let steady = network.steady(&[power]).expect("steady")[0];
    let expected = AMBIENT + power / node.conductance_w_per_k();
    assert!((steady - expected).abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// The gate, and the regime it protects.
// ---------------------------------------------------------------------------

#[test]
fn the_declared_ceiling_matches_the_corpus_row() {
    // The ceiling is duplicated here because fs-conduction must not depend on
    // the corpus registry. This test is what stops the duplicate drifting.
    let case = fs_vvreg::thermal_level_a::thermal_level_a_cases()
        .iter()
        .find(|case| case.id == "thermal-a-lumped-transient")
        .expect("Level-A lumped row");
    let corpus_ceiling = case
        .context
        .iter()
        .find(|entry| entry.name == "biot-number")
        .map(|entry| entry.hi)
        .expect("the row declares a Biot ceiling");
    assert!(
        (LUMPED_BIOT_CEILING - corpus_ceiling).abs() < 1e-12,
        "the crate ceiling {LUMPED_BIOT_CEILING} drifted from the corpus row's {corpus_ceiling}"
    );
}

#[test]
fn the_gate_admits_inside_the_regime_and_refuses_outside_it() {
    let gate = BiotGate::corpus_default();

    match gate.adjudicate(&network_at(0.05)) {
        ValidityVerdict::Admitted { worst_biot } => {
            assert!((worst_biot - 0.05).abs() < 1e-12);
        }
        other => panic!("small Biot must be admitted, got {other:?}"),
    }

    match gate.adjudicate(&network_at(2.0)) {
        ValidityVerdict::Refused {
            node,
            biot,
            ceiling,
        } => {
            assert_eq!(node, "cube");
            assert!((biot - 2.0).abs() < 1e-12);
            assert!((ceiling - LUMPED_BIOT_CEILING).abs() < 1e-12);
        }
        other => panic!("large Biot must be refused, got {other:?}"),
    }
}

#[test]
fn solving_outside_the_regime_refuses_rather_than_returning_a_cheap_wrong_number() {
    let error = solve_gated(
        &network_at(2.0),
        BiotGate::corpus_default(),
        &[0.0],
        &[AMBIENT + EXCESS],
        1.0,
    )
    .expect_err("outside the regime must refuse");
    match error {
        ConductionError::ScenarioRow { region, what, fix } => {
            assert_eq!(region, "cube");
            assert!(what.contains("not isothermal"), "what: {what}");
            assert!(
                what.contains("escalating"),
                "the refusal must point at the full rung, not at loosening the gate: {what}"
            );
            assert!(fix.contains("escalate"));
        }
        other => panic!("expected a scenario-row refusal, got {other:?}"),
    }
}

#[test]
fn a_loosened_gate_reports_the_ceiling_it_actually_applied() {
    // Loosening is the caller's decision to own, so it must be visible in the
    // verdict rather than indistinguishable from the corpus default.
    let loose = BiotGate::at(5.0).expect("loosened gate");
    assert!((loose.ceiling() - 5.0).abs() < 1e-12);
    assert!(loose.adjudicate(&network_at(2.0)).admitted());
    assert!(
        !BiotGate::corpus_default()
            .adjudicate(&network_at(2.0))
            .admitted(),
        "the corpus gate must still refuse what a loosened one admits"
    );
    assert!(BiotGate::at(0.0).is_err());
    assert!(BiotGate::at(f64::NAN).is_err());
}

// ---------------------------------------------------------------------------
// THE point of the bead: cost and authority are separate axes.
// ---------------------------------------------------------------------------

#[test]
fn the_reduced_rung_agrees_with_the_full_one_inside_the_regime() {
    // Cheap does not mean wrong — inside its declared regime the reduced rung
    // must reproduce what the expensive one says, or it is not a rung of the
    // same graph.
    let biot = 0.05;
    let network = network_at(biot);
    let node = &network.nodes()[0];
    let reduced = solve_gated(
        &network,
        BiotGate::corpus_default(),
        &[0.0],
        &[AMBIENT + EXCESS],
        node.time_constant_s(),
    )
    .expect("inside the regime")
    .temperature_k[0];

    let full = full_rung_mean_at_time_constant(biot);
    let relative = (reduced - full).abs() / (full - AMBIENT);
    assert!(
        relative < 0.05,
        "inside the regime the rungs must agree: reduced {reduced} vs full {full} ({:.2}% apart)",
        relative * 100.0
    );
}

#[test]
fn the_reduced_rung_is_visibly_wrong_outside_the_regime_which_is_why_the_gate_exists() {
    // This is the test that justifies the gate rather than asserting it. At
    // high Biot the body is NOT isothermal: the surface cools fast while the
    // core lags, so the real mean retains more heat than a lumped model
    // predicts. The reduced rung is therefore not merely less precise here —
    // it is wrong by a margin that would change a decision.
    let biot = 2.0;
    let network = network_at(biot);
    let node = &network.nodes()[0];

    // Compute it anyway, THROUGH A LOOSENED GATE, purely to measure the error
    // the corpus gate is protecting against.
    let reduced = solve_gated(
        &network,
        BiotGate::at(10.0).expect("deliberately loosened"),
        &[0.0],
        &[AMBIENT + EXCESS],
        node.time_constant_s(),
    )
    .expect("loosened gate admits")
    .temperature_k[0];

    let full = full_rung_mean_at_time_constant(biot);
    let relative = (reduced - full).abs() / (full - AMBIENT);
    assert!(
        relative > 0.10,
        "outside the regime the rungs must visibly DISAGREE, else the gate is protecting nothing: \
         reduced {reduced} vs full {full} ({:.2}% apart)",
        relative * 100.0
    );

    // And the real body retains more heat than the lumped model claims,
    // which is the physical direction: the lagging core keeps the mean up.
    assert!(
        full > reduced,
        "a non-isothermal body should retain MORE heat than lumping predicts: full {full}, reduced {reduced}"
    );

    // The corpus gate refuses exactly this case.
    assert!(
        !BiotGate::corpus_default().adjudicate(&network).admitted(),
        "the gate must refuse the case just shown to be wrong"
    );
}

// ---------------------------------------------------------------------------
// Extraction and admission.
// ---------------------------------------------------------------------------

#[test]
fn extraction_recovers_the_conductance_that_produced_the_rise() {
    // hA = P / dT is the extraction. Round-tripping it is the check that the
    // arithmetic is the stated one; it says nothing about whether the rise it
    // was extracted from was itself trustworthy.
    let truth = node_at(0.05);
    let power = 900.0;
    let rise = power / truth.conductance_w_per_k();
    let extracted =
        extract_node_from_steady_rise("cube", power, rise, RHO_CP, CHAR_LENGTH, K, AREA)
            .expect("extraction admits");
    assert!(
        (extracted.conductance_w_per_k() - truth.conductance_w_per_k()).abs()
            < 1e-9 * truth.conductance_w_per_k()
    );
    assert!((extracted.biot() - truth.biot()).abs() < 1e-12);
}

#[test]
fn admission_refuses_degenerate_declarations() {
    assert!(
        LumpedNode::new("  ", 1.0, 1.0, 1.0, 1.0, 1.0).is_err(),
        "blank name"
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            LumpedNode::new("n", bad, 1.0, 1.0, 1.0, 1.0).is_err(),
            "capacitance {bad}"
        );
        assert!(
            LumpedNode::new("n", 1.0, bad, 1.0, 1.0, 1.0).is_err(),
            "conductance {bad}"
        );
    }
    assert!(LumpedNetwork::new(Vec::new(), AMBIENT).is_err(), "no nodes");
    assert!(
        LumpedNetwork::new(vec![node_at(0.05), node_at(0.05)], AMBIENT).is_err(),
        "duplicate node name"
    );
    assert!(
        LumpedNetwork::new(vec![node_at(0.05)], f64::NAN).is_err(),
        "non-finite ambient"
    );

    let network = network_at(0.05);
    assert!(network.steady(&[]).is_err(), "power length mismatch");
    assert!(network.steady(&[-1.0]).is_err(), "negative node power");
    assert!(
        network.response_at(&[0.0], &[AMBIENT], -1.0).is_err(),
        "negative time"
    );
    assert!(
        network.response_at(&[0.0], &[], 1.0).is_err(),
        "initial length mismatch"
    );
    assert!(
        extract_node_from_steady_rise("n", 0.0, 1.0, 1.0, 1.0, 1.0, 1.0).is_err(),
        "zero extraction power"
    );
}

#[test]
fn the_reduced_rung_is_deterministic() {
    let network = network_at(0.05);
    let solve = || {
        solve_gated(
            &network,
            BiotGate::corpus_default(),
            &[400.0],
            &[AMBIENT + EXCESS],
            1000.0,
        )
        .expect("solve")
    };
    assert_eq!(solve(), solve());
}

#[test]
fn g1_enthalpy_march_conserves_energy_through_the_latent_heat_plateau() {
    let curve = phase_curve();
    let body = enthalpy_body(&curve);
    let march = with_cx(|cx| {
        solve_lumped_enthalpy(cx, &body, BiotGate::corpus_default(), enthalpy_config())
            .expect("Biot-admitted phase march")
    });

    let final_state = march.samples().last().expect("final sample").phase_state;
    assert_eq!(final_state.phase(), SolidLiquidPhase::SolidLiquid);
    assert!((final_state.temperature_k() - 600.0).abs() < 1.0e-10);
    assert!((final_state.specific_enthalpy_j_kg() - 45_000.0).abs() < 1.0e-6);
    assert!((final_state.liquid_mass_fraction() - 0.6).abs() < 1.0e-10);
    assert!(
        march.cumulative_absolute_energy_residual_j() < 1.0e-4,
        "discrete energy ledger drifted by {} J",
        march.cumulative_absolute_energy_residual_j()
    );
}

#[test]
fn g1_hot_environment_changes_phase_via_convection_and_radiation() {
    let curve = phase_curve();
    let body = LumpedEnthalpyBody::try_new(
        "radiatively heated body",
        0.1,
        0.02,
        20.0,
        0.5,
        0.001,
        100.0,
        &curve,
    )
    .expect("admitted body");
    let config = LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: 20_000.0,
        ambient_temperature_k: 1_000.0,
        internal_power_w: 0.0,
        // Long enough to enter the latent-heat interval, but intentionally
        // short of the phase curve's 95 kJ/kg evidence boundary. The solver
        // must refuse rather than extrapolate once that boundary is reached.
        duration_s: 10.0,
        maximum_step_s: 0.25,
        maximum_steps: 40,
        enthalpy_tolerance_j_kg: 1.0e-8,
    };
    let march = with_cx(|cx| {
        solve_lumped_enthalpy(cx, &body, BiotGate::corpus_default(), config)
            .expect("hot ambient march")
    });
    let initial = march.samples()[0].phase_state;
    let final_state = march.samples().last().expect("final sample").phase_state;
    assert_eq!(initial.phase(), SolidLiquidPhase::Solid);
    assert!(final_state.specific_enthalpy_j_kg() > initial.specific_enthalpy_j_kg());
    assert!(final_state.liquid_mass_fraction() > 0.0);
    assert!(
        march
            .samples()
            .iter()
            .skip(1)
            .all(|sample| sample.convection_into_body_w > 0.0
                && sample.radiation_into_body_w > 0.0)
    );
}

#[test]
fn g1_card_backed_transport_drives_the_same_phase_curve_without_extrapolation() {
    let card = thermal_material_card();
    let curve = phase_curve_for(card.content_hash());
    let transport = LumpedThermalTransport::from_material_card(
        &card,
        "thermal_conductivity",
        &[300.0, 600.0, 800.0],
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("card-backed conductivity and emissivity");
    assert_eq!(
        transport.material_card_identity(),
        Some(card.content_hash())
    );
    let body = LumpedEnthalpyBody::try_new_with_transport(
        "card-backed phase body",
        0.1,
        0.02,
        20.0,
        0.001,
        transport,
        &curve,
    )
    .expect("body and transport share one card");
    let config = LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: 20_000.0,
        ambient_temperature_k: 1_000.0,
        internal_power_w: 0.0,
        duration_s: 10.0,
        maximum_step_s: 0.25,
        maximum_steps: 40,
        enthalpy_tolerance_j_kg: 1.0e-8,
    };
    let march = with_cx(|cx| {
        solve_lumped_enthalpy(cx, &body, BiotGate::corpus_default(), config)
            .expect("card-backed hot ambient march")
    });
    assert!(
        march
            .samples()
            .last()
            .expect("final sample")
            .phase_state
            .liquid_mass_fraction()
            > 0.0
    );

    let mismatched_curve = phase_curve();
    assert!(
        LumpedEnthalpyBody::try_new_with_transport(
            "mismatched card body",
            0.1,
            0.02,
            20.0,
            0.001,
            body.transport().clone(),
            &mismatched_curve,
        )
        .is_err(),
        "transport from one material card must not drive another card's phase curve"
    );
}

#[test]
fn g0_enthalpy_march_refuses_invalid_fidelity_or_missing_phase_domain() {
    let curve = phase_curve();
    let nonuniform =
        LumpedEnthalpyBody::try_new("nonuniform body", 1.0, 1.0, 20.0, 0.8, 0.1, 1.0, &curve)
            .expect("body declaration itself is valid");
    assert!(
        with_cx(|cx| solve_lumped_enthalpy(
            cx,
            &nonuniform,
            BiotGate::corpus_default(),
            enthalpy_config(),
        ))
        .is_err()
    );

    let body = enthalpy_body(&curve);
    let mut outside = enthalpy_config();
    outside.initial_specific_enthalpy_j_kg = 90_000.0;
    outside.internal_power_w = 1.0e6;
    outside.duration_s = 1.0;
    outside.maximum_step_s = 1.0;
    outside.maximum_steps = 1;
    assert!(
        with_cx(|cx| solve_lumped_enthalpy(cx, &body, BiotGate::corpus_default(), outside,))
            .is_err()
    );

    let mut unbounded_steps = enthalpy_config();
    unbounded_steps.duration_s = f64::MAX;
    unbounded_steps.maximum_step_s = f64::MIN_POSITIVE;
    assert!(
        with_cx(
            |cx| solve_lumped_enthalpy(cx, &body, BiotGate::corpus_default(), unbounded_steps,)
        )
        .is_err()
    );
}
