//! Closed-form latent endpoints distinguish scheduled source feedback from snapshots.
use super::*;
use fs_conduction::bc::history::march::enthalpy::adaptive::PhaseErrorTolerance;
use fs_conduction::bc::history::march::enthalpy::radiation::material::{
    MaterialRadiationRegionHistory, MaterialRadiativeEnthalpyHistoryRun,
};
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::transient::enthalpy::EnthalpyStepSolution;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::QueryPoint;

const LATENT: f64 = 1000.0;

fn source_card(upper_fraction: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "melt".into(), abscissa_dims: Dims::NONE,
            knots: vec![(0.0, 0.2), (upper_fraction, 0.2 + 0.6 * upper_fraction)],
            dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0)
            .with("melt", 0.0, upper_fraction).with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(), provenance: Provenance {
            source: "synthetic scheduled latent radiation oracle".into(),
            license: "test-only; not measured".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "chart".into(),
        process: "history oracle".into(), revision: 0 }, claims, Vec::new()).unwrap()
}

fn phase_chart(card: &MaterialCard, datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), [
        (-100.0, 300.0, 0.0), (0.0, 400.0, 0.0),
        (LATENT, 400.0, 1.0), (LATENT + 100.0, 500.0, 1.0),
    ].map(|(h, t, f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: datum + h,
        temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: 24.0 }).to_vec()).unwrap()
}

fn context() -> QueryPoint {
    QueryPoint::new().with("T", 400.0).unwrap().with("melt", 0.5).unwrap()
        .with("lot", 7.0).unwrap()
}

fn declaration<'a>(card: &'a MaterialCard, curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint, signal: TimeSignal) -> MaterialRadiationRegionHistory<'a>
{
    MaterialRadiationRegionHistory { region: "surface".into(), curve,
        material: ThermalMaterialPoint { material: card, context, temperature_axis: "T",
            liquid_fraction_axis: Some("melt") },
        selection: MaterialPropertySelection::SingleClaimOnly, ambient_temperature: signal }
}

fn face_history<'m>(cx: &fs_exec::Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .region("surface", |face| face.centroid[2] == 0.0, ThermalBc::adiabatic()).unwrap()
        .remainder("insulated", ThermalBc::adiabatic()).unwrap(),
        &Environment::earth_lab(), Vec::new(), 0.0).unwrap()
}

fn base<'h, 'm, 'c>(cx: &fs_exec::Cx<'_>, boundary: &'h ThermalBoundaryHistory<'m>,
    conductivity: &'h ConductivityModel, curve: &'c EquilibriumEnthalpyPhaseCurve,
    initial_h: f64, steps: usize) -> EnthalpyHistoryRun<'h, 'm, 'c>
{
    EnthalpyHistoryRun::new(cx, problem(boundary, conductivity,
        HeatSourceHistory::Fixed(ScalarField::Uniform(0.0))), EnthalpyHistoryMaterials {
            table: &[ReferenceEnthalpyMaterial { curve, reference_density_kg_m3: 24.0 }],
            of_element: &[0], phase_conductivity: None },
        &[initial_h; 4], budget(boundary.mesh(), steps)).unwrap()
}

fn constant(temperature: f64) -> TimeSignal {
    TimeSignal::Constant(QtyAny::new(temperature, TEMPERATURE_DIMS))
}

// On the z=0 face, three equal reference masses are each 1 kg. Their T stays
// 400 K, so conduction vanishes even though their latent h differs from the apex.
// 3(h1-h0) = -dt A sigma(Tm^4-Tb^4) (0.2 + 0.6 h1/L).
fn endpoint(old: f64, dt: f64, bath: f64) -> f64 {
    let c = 0.5 * STEFAN_BOLTZMANN_W_M2_K4 * (400.0_f64.powi(4) - bath.powi(4));
    (3.0 * old - dt * c * 0.2) / (3.0 + dt * c * 0.6 / LATENT)
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-8, "{actual:.16e} != {expected:.16e}");
}

fn admit_report(report: &EnthalpyStepSolution, expected: f64, datum: f64,
    old: f64, dt: f64, bath: f64)
{
    assert_eq!(report.temperature, vec![400.0; 4]);
    for vertex in 0..3 {
        close(report.specific_enthalpy_j_kg[vertex] - datum, expected);
        close(report.liquid_mass_fraction[vertex], expected / LATENT);
    }
    close(report.specific_enthalpy_j_kg[3] - datum, 500.0);
    let radiation = report.pointwise_radiation.as_ref().unwrap();
    let expected_power = 0.5 * STEFAN_BOLTZMANN_W_M2_K4
        * (400.0_f64.powi(4) - bath.powi(4)) * (0.2 + 0.6 * expected / LATENT);
    close(radiation.outward_heat_w, expected_power);
    close(report.stored_energy_change_j, 3.0 * (expected - old));
    close(report.stored_energy_change_j + dt * radiation.outward_heat_w, 0.0);
    assert!(report.energy_residual_j.abs() <= config().energy_tolerance_j);
    assert_eq!(report.robin_out_w, 0.0);
    assert_eq!(radiation.patches[0].material_states.len(), 3);
    for state in &radiation.patches[0].material_states {
        close(state.emissivity.value(), 0.2 + 0.6 * expected / LATENT);
    }
}

#[test]
fn direct_latent_heating_and_cooling_follow_original_epsilon_of_phase_and_enthalpy_datum() {
    let mesh = mesh(false); let card = source_card(1.0); let context = context();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    for datum in [0.0, 4096.0] {
        let curve = phase_chart(&card, datum);
        with_cx(|cx| {
            let boundary = face_history(cx, &mesh);
            let signal = TimeSignal::Table { times: vec![0.0, 0.2], values: vec![500.0, 300.0],
                dims: TEMPERATURE_DIMS, interp: Interp::Hold };
            let mut run = MaterialRadiativeEnthalpyHistoryRun::new(cx,
                base(cx, &boundary, &conductivity, &curve, datum + 500.0, 2),
                vec![declaration(&card, &curve, &context, signal)]).unwrap();
            let selected = run.reference_emissivities()[0].receipt().selected;
            let mut old = 500.0;
            for (time, bath) in [(0.1, 500.0), (0.2, 300.0)] {
                let expected = endpoint(old, 0.1, bath);
                let report = run.advance_to(cx, time, config()).unwrap();
                admit_report(&report, expected, datum, old, 0.1, bath);
                for state in &report.pointwise_radiation.as_ref().unwrap().patches[0].material_states {
                    assert_eq!(state.emissivity.receipt().selected, selected);
                    assert_eq!(state.emissivity.card_identity(), card.content_hash());
                }
                old = expected;
            }
            assert_eq!(run.conduction().accepted_steps(), 2);
        });
    }
}

#[test]
fn held_reservoir_pulse_clips_every_event_and_keeps_interval_end_left_heat() {
    let mesh = mesh(false); let card = source_card(1.0); let curve = phase_chart(&card, 0.0);
    let context = context(); let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    with_cx(|cx| {
        let boundary = face_history(cx, &mesh);
        let signal = TimeSignal::Table { times: vec![0.0, 0.125, 0.375],
            values: vec![400.0, 500.0, 300.0], dims: TEMPERATURE_DIMS, interp: Interp::Hold };
        let mut run = MaterialRadiativeEnthalpyHistoryRun::new(cx,
            base(cx, &boundary, &conductivity, &curve, 500.0, 3),
            vec![declaration(&card, &curve, &context, signal)]).unwrap();
        let mut old = 500.0; let mut previous = 0.0; let mut heat_in = 0.0;
        for (end, bath) in [(0.125, 400.0), (0.375, 500.0), (0.5, 300.0)] {
            let step = run.advance_event_step(cx, 0.5, 1.0, config()).unwrap();
            assert_eq!(step.start_time_s, previous); assert_eq!(step.end_time_s, end);
            let dt = end - previous; let expected = endpoint(old, dt, bath);
            admit_report(&step.report, expected, 0.0, old, dt, bath);
            heat_in -= dt * step.report.pointwise_radiation.as_ref().unwrap().outward_heat_w;
            old = expected; previous = end;
        }
        close(3.0 * (run.conduction().specific_enthalpy_j_kg()[0] - 500.0), heat_in);
        assert_eq!(run.conduction().time_s(), 0.5);
    });
}

#[test]
fn adaptive_phase_gate_rejects_constant_temperature_trials_and_publishes_both_half_steps() {
    let mesh = mesh(false); let card = source_card(1.0); let curve = phase_chart(&card, 0.0);
    let context = context(); let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    with_cx(|cx| {
        let boundary = face_history(cx, &mesh);
        let mut run = MaterialRadiativeEnthalpyHistoryRun::new(cx,
            base(cx, &boundary, &conductivity, &curve, 500.0, 2),
            vec![declaration(&card, &curve, &context, constant(500.0))]).unwrap();
        let result = run.advance_adaptive_step(cx, 0.25, 0.25, config(), AdaptiveHistoryConfig {
            absolute_temperature_k: 1e-10, relative_temperature: 0.0,
            min_step_s: 0.001, max_step_s: 0.25, max_attempts: 8,
        }, PhaseErrorTolerance { absolute_enthalpy_j_kg: 1.0, relative_enthalpy: 0.0,
            absolute_liquid_fraction: 5e-5 }).unwrap();
        assert!(result.rejected_attempts > 0); assert!(result.end_time_s < 0.25);
        assert!(result.normalized_error <= 1.0);
        let first_dt = result.midpoint_time_s - result.start_time_s;
        let second_dt = result.end_time_s - result.midpoint_time_s;
        let first = endpoint(500.0, first_dt, 500.0);
        let second = endpoint(first, second_dt, 500.0);
        admit_report(&result.first_half, first, 0.0, 500.0, first_dt, 500.0);
        admit_report(&result.second_half, second, 0.0, first, second_dt, 500.0);
        assert_eq!(run.conduction().accepted_steps(), 2);
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), result.second_half.specific_enthalpy_j_kg);
    });
}

#[test]
fn source_refusal_and_cancellation_preserve_resume_and_a_cloned_chart_is_not_original_storage() {
    let mesh = mesh(false); let card = source_card(0.6); let curve = phase_chart(&card, 0.0);
    let other_curve = curve.clone(); let context = context();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let boundary = with_cx(|cx| face_history(cx, &mesh));
    with_cx(|cx| {
        assert!(MaterialRadiativeEnthalpyHistoryRun::new(cx,
            base(cx, &boundary, &conductivity, &curve, 500.0, 4),
            vec![declaration(&card, &other_curve, &context, constant(500.0))]).is_err());
    });
    let mut run = with_cx(|cx| MaterialRadiativeEnthalpyHistoryRun::new(cx,
        base(cx, &boundary, &conductivity, &curve, 500.0, 4),
        vec![declaration(&card, &curve, &context, constant(500.0))]).unwrap());
    with_cx(|cx| { run.advance_to(cx, 0.1, config()).unwrap(); });
    let before = (run.conduction().specific_enthalpy_j_kg().to_vec(),
        run.conduction().temperature().to_vec(), run.conduction().liquid_mass_fraction().to_vec());
    with_cancelled_cx(|cx| assert!(run.advance_event_step(cx, 0.2, 0.1, config()).is_err()));
    with_cx(|cx| {
        // The endpoint would require f > 0.6, beyond the pinned source chart.
        assert!(run.advance_to(cx, 1.1, config()).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), before.0);
        assert_eq!(run.conduction().temperature(), before.1);
        assert_eq!(run.conduction().liquid_mass_fraction(), before.2);
        assert_eq!(run.conduction().time_s(), 0.1); assert_eq!(run.conduction().accepted_steps(), 1);
        let report = run.advance_to(cx, 0.2, config()).unwrap();
        let first = endpoint(500.0, 0.1, 500.0); let second = endpoint(first, 0.1, 500.0);
        admit_report(&report, second, 0.0, first, 0.1, 500.0);
    });
}
