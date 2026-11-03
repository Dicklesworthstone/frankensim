//! Exact P1 source moments, latent history and immutable material support.
use super::*;
use fs_conduction::bc::history::march::enthalpy::PhaseEnthalpyHistoryProblem;
use fs_conduction::bc::history::march::enthalpy::adaptive::PhaseErrorTolerance;
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::transient::enthalpy::material::phase::MaterialPhaseTransport;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard,
    MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{QuantitySpec, Temperature};

fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("k", THERMAL_CONDUCTIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "f".into(), abscissa_dims: Dims::NONE,
            knots: vec![(0.0, 2.0), (1.0, 8.0)], dims: THERMAL_CONDUCTIVITY_DIMS },
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 250.0, 1000.0)
            .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated, observations: Vec::new(),
        provenance: Provenance { source: "manufactured P1 history moment oracle".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic source".into(),
        phase: "equilibrium chart".into(), process: "numerical fixture".into(), revision: 0 },
        claims, Vec::new()).unwrap()
}

fn chart(card: &MaterialCard, datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), [
        (-100.0, 300.0, 0.0), (0.0, 400.0, 0.0),
        (100.0, 400.0, 1.0), (500.0, 800.0, 1.0),
    ].map(|(h, t, f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: datum + h,
        temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: 2000.0 }).to_vec()).unwrap()
}

fn source<'a>(card: &MaterialCard, curve: &'a EquilibriumEnthalpyPhaseCurve,
    bounds: [f64; 2]) -> PreparedPhaseConductivity<'a>
{
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 400.0).unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.5).unwrap();
    PreparedPhaseConductivity::new(curve, ThermalMaterialPoint { material: card,
        context: &context, temperature_axis: "T", liquid_fraction_axis: Some("f") },
        &PropertyKey::new("k", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly), bounds, || Ok(())).unwrap()
}

fn heated<'m>(cx: &fs_exec::Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .region("bath", |face| face.centroid[2] == 0.0, ThermalBc::dirichlet(410.0).unwrap()).unwrap()
        .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.0).unwrap()
}

fn histories<'h, 'm>(boundary: &'h ThermalBoundaryHistory<'m>, source: HeatSourceHistory)
    -> PhaseEnthalpyHistoryProblem<'h, 'm>
{
    PhaseEnthalpyHistoryProblem { boundary, interfaces: None, source }
}

// Three face nodes have h=110, T=410 K. The free apex is latent: 0<x<100,
// T=400 K. A tetra barycentric coordinate has density 3(1-lambda)^2.
// Integrating the ORIGINAL k=2+6f over h_h=110+(x-110)*lambda gives this
// mean, including the liquid/latent crossing inside the element.
fn mean_k(x: f64) -> f64 {
    8.0 - 0.015 * (100.0 - x).powi(4) / (110.0 - x).powi(3)
}

// Independent scalar backward-Euler equation, not the production moment leaf.
fn endpoint(old: f64, dt: f64, mass: f64, volumetric_source: f64) -> f64 {
    let residual = |x| mass * (x - old) - dt * ((10.0 / 6.0) * mean_k(x) + volumetric_source / 24.0);
    let (mut lo, mut hi) = (0.0, 100.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..70 {
        let middle = 0.5 * (lo + hi);
        if residual(middle) < 0.0 { lo = middle; } else { hi = middle; }
    }
    0.5 * (lo + hi)
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-8, "{actual:.16e} != {expected:.16e}");
}

#[test]
fn g1_original_k_of_phase_history_matches_exact_p1_moments_and_reference_mass() {
    let mesh = mesh(false);
    let card = card();
    for datum in [0.0, 4096.0] {
        let curve = chart(&card, datum);
        let source = source(&card, &curve, [datum - 100.0, datum + 500.0]);
        with_cx(|cx| {
            let materials = [MaterialPhaseTransport::try_new_isotropic(cx, &source, 24.0).unwrap()];
            let boundary = heated(cx, &mesh);
            let initial = [datum + 110.0, datum + 110.0, datum + 110.0, datum + 20.0];
            for nodal_density in [None, Some([12.0, 24.0, 36.0, 48.0])] {
                let p = histories(&boundary, HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)));
                let mut run = match nodal_density {
                    None => EnthalpyHistoryRun::from_phase_materials(cx, p, &materials,
                        &[0], &initial, budget(&mesh, 4)).unwrap(),
                    Some(ref density) => EnthalpyHistoryRun::from_phase_materials_with_nodal_reference_density(
                        cx, p, &materials, &[0], density, &initial, budget(&mesh, 4)).unwrap(),
                };
                let mass = if nodal_density.is_some() { 1.4 } else { 1.0 };
                close(run.storage().reference_nodal_masses_kg()[3], mass);
                assert!(std::ptr::eq(run.storage().phase_conductivity_source(0).unwrap(), &source));
                assert!(std::ptr::eq(run.storage().phase_curve_at(3).unwrap(), &curve));
                let mut old = 20.0;
                for step in 1..=4 {
                    let expected = endpoint(old, 0.5, mass, 0.0);
                    let report = run.advance_to(cx, 0.5 * f64::from(step), config()).unwrap();
                    close(report.specific_enthalpy_j_kg[3] - datum, expected);
                    close(report.liquid_mass_fraction[3], expected / 100.0);
                    assert_eq!(report.temperature, vec![410.0, 410.0, 410.0, 400.0]);
                    close(report.stored_energy_change_j, mass * (expected - old));
                    close(0.5 * report.dirichlet_in_w, report.stored_energy_change_j);
                    assert!(report.energy_residual_j.abs() <= config().energy_tolerance_j);
                    let flux = run.storage().phase_transport_heat_flux(cx, run.specific_enthalpy_j_kg()).unwrap();
                    close(flux[0][2], 10.0 * mean_k(expected));
                    // Averaging source answers at nodes is physically different.
                    assert!((mean_k(expected) - (24.0 + 2.0 + 0.06 * expected) / 4.0).abs() > 0.1);
                    old = expected;
                }
            }
        });
    }
}

#[test]
fn g1_source_histories_resolve_held_pulses_and_adapt_on_the_actual_latent_state() {
    let mesh = mesh(false); let card = card(); let curve = chart(&card, 0.0);
    let source = source(&card, &curve, [-100.0, 500.0]);
    with_cx(|cx| {
        let materials = [MaterialPhaseTransport::try_new_isotropic(cx, &source, 24.0).unwrap()];
        let boundary = heated(cx, &mesh);
        let create = || EnthalpyHistoryRun::from_phase_materials(cx, histories(&boundary,
            HeatSourceHistory::Uniform(TimeSignal::Table { times: vec![0.0, 0.5],
                values: vec![480.0, -480.0], dims: VOLUMETRIC_SOURCE_DIMS, interp: Interp::Hold })),
            &materials, &[0], &[110.0, 110.0, 110.0, 20.0], budget(&mesh, 20)).unwrap();
        let mut run = create();
        let first = run.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
        assert_eq!(first.end_time_s, 0.5);
        let h1 = endpoint(20.0, 0.5, 1.0, 480.0);
        close(first.report.specific_enthalpy_j_kg[3], h1);
        close(first.report.source_w, 80.0);
        let second = run.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
        assert_eq!(second.end_time_s, 1.0);
        close(second.report.specific_enthalpy_j_kg[3], endpoint(h1, 0.5, 1.0, -480.0));
        close(second.report.source_w, -80.0);
        close(first.report.stored_energy_change_j + second.report.stored_energy_change_j,
            0.5 * (first.report.dirichlet_in_w + second.report.dirichlet_in_w));

        let mut adaptive_run = create();
        let adaptive = AdaptiveHistoryConfig { absolute_temperature_k: 1.0,
            relative_temperature: 0.0, min_step_s: 1e-6, max_step_s: 0.5,
            max_attempts: 30 };
        let phase = PhaseErrorTolerance { absolute_enthalpy_j_kg: 1.0,
            relative_enthalpy: 0.0, absolute_liquid_fraction: 1e-5 };
        let accepted = adaptive_run.advance_adaptive_step(cx, 0.5, 0.5,
            config(), adaptive, phase).unwrap();
        assert!(accepted.end_time_s < 0.5);
        let half = accepted.end_time_s / 2.0;
        let hm = endpoint(20.0, half, 1.0, 480.0);
        let hf = endpoint(hm, half, 1.0, 480.0);
        close(accepted.first_half.specific_enthalpy_j_kg[3], hm);
        close(accepted.second_half.specific_enthalpy_j_kg[3], hf);
        close(adaptive_run.specific_enthalpy_j_kg()[3], hf);
        assert_eq!(adaptive_run.accepted_steps(), 2);
        assert_eq!(accepted.first_half.temperature[3], 400.0);
        assert_eq!(accepted.second_half.temperature[3], 400.0);
    });
}

#[test]
fn g0_original_conductivity_support_and_cancellation_preserve_accepted_history() {
    let mesh = mesh(false); let card = card(); let curve = chart(&card, 0.0);
    let source = source(&card, &curve, [40.0, 60.0]);
    with_cx(|cx| {
        let materials = [MaterialPhaseTransport::try_new_isotropic(cx, &source, 24.0).unwrap()];
        let boundary = natural(cx, &mesh);
        let create = |initial: f64| EnthalpyHistoryRun::from_phase_materials(cx,
            histories(&boundary, HeatSourceHistory::Fixed(ScalarField::Uniform(480.0))),
            &materials, &[0], &[initial; 4], budget(&mesh, 3));
        assert!(matches!(create(65.0), Err(HeterogeneousEnthalpyError::Enthalpy(
            EnthalpyError::SourceEnthalpyOutsideDomain { .. }))));
        let mut run = create(50.0).unwrap();
        with_cancelled_cx(|cancelled| assert!(matches!(run.advance_to(cancelled, 0.25, config()),
            Err(HeterogeneousEnthalpyError::Enthalpy(EnthalpyError::Conduction(
                ConductionError::Cancelled { .. }))))));
        assert!(run.advance_to(cx, 1.0, config()).is_err());
        assert_eq!(run.specific_enthalpy_j_kg(), &[50.0; 4]);
        assert_eq!(run.temperature(), &[400.0; 4]);
        assert_eq!(run.liquid_mass_fraction(), &[0.5; 4]);
        assert_eq!(run.time_s(), 0.0);
        assert_eq!(run.accepted_steps(), 0);
        let accepted = run.advance_to(cx, 0.25, config()).unwrap();
        for h in &accepted.specific_enthalpy_j_kg { close(*h, 55.0); }
        close(accepted.stored_energy_change_j, 20.0);
        assert_eq!(run.accepted_steps(), 1);
    });
}
