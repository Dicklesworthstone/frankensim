//! G1 latent history against independent two-patch reflection and heat accounts.
use super::*;
use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{
    EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun,
    adaptive::PhaseErrorTolerance,
    radiation::enclosure::MaterialEnclosureEnthalpyHistoryRun,
};
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::transient::enthalpy::EnthalpyStepSolution;
use fs_conduction::{HTC_DIMS, TEMPERATURE_DIMS, VOLUMETRIC_SOURCE_DIMS};
use fs_qty::QtyAny;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};

fn constant(t: f64) -> TimeSignal { TimeSignal::Constant(QtyAny::new(t, TEMPERATURE_DIMS)) }

fn history<'m>(cx: &fs_exec::Cx<'_>, mesh: &'m ConductionMesh, robin: bool)
    -> ThermalBoundaryHistory<'m>
{
    let rows = if robin { vec![ThermalRegionHistory {
        row: BoundaryCondition { region: "right".into(), physics: Physics::Thermal,
            kind: BcKind::Robin, value: Some(BcValue::Uniform(QtyAny::new(5.0, HTC_DIMS))),
            compatibility: None, frame: 0 },
        reference_temperature: Some(TimeSignal::Table { times: vec![0.0, 0.3],
            values: vec![350.0, 400.0], dims: TEMPERATURE_DIMS, interp: Interp::Hold }),
    }] } else { Vec::new() };
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .region("left", |face| face.centroid[2] == 0.0 && face.vertices[0] < 4,
            ThermalBc::adiabatic()).unwrap()
        .region("right", |face| face.centroid[2] == 0.0 && face.vertices[0] >= 4,
            if robin { ThermalBc::robin(5.0, 350.0).unwrap() } else { ThermalBc::adiabatic() }).unwrap()
        .remainder("insulated", ThermalBc::adiabatic()).unwrap(),
        &Environment::earth_lab(), rows, 0.0).unwrap()
}

fn run<'h, 'm, 'c>(cx: &fs_exec::Cx<'_>, boundary: &'h ThermalBoundaryHistory<'m>,
    material: &'h ConductivityModel, curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    datum: f64, source: HeatSourceHistory, steps: usize) -> EnthalpyHistoryRun<'h, 'm, 'c>
{
    EnthalpyHistoryRun::new(cx, ThermalHistoryProblem { boundary, material, source,
        element_materials: None, interfaces: None }, EnthalpyHistoryMaterials {
        table: &[ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: RHO },
            ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: RHO }],
        of_element: &[0, 1], phase_conductivity: None }, &[datum + 500.0; 8],
        EnthalpyHistoryBudget { spatial: budget(), max_steps: steps }).unwrap()
}

fn zero() -> HeatSourceHistory { HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)) }

fn open(cx: &fs_exec::Cx<'_>, mesh: &ConductionMesh, cards: &[MaterialCard; 2],
    context: &QueryPoint) -> OpenDiffuseEnclosure
{
    OpenDiffuseEnclosure::new(cx, surfaces(mesh, cards, context),
        vec![vec![0.0, 0.4], vec![0.4, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "synthetic reciprocal history patches".into() },
        ViewFactorTolerance::default()).unwrap()
}

// Solve the two gray equations directly, without production factors/quadrature.
// Both traces remain on their own latent plateaus; epsilon depends on accepted h.
fn heat(h: [f64; 2], mutual: f64, bath: f64) -> [f64; 2] {
    let e = h.map(|h| 0.2 + 0.6 * h / LATENT);
    let black = [400.0_f64, 350.0].map(|t| SIGMA * t.powi(4));
    let ambient = SIGMA * bath.powi(4);
    let rhs: [f64; 2] = std::array::from_fn(|i|
        e[i] * black[i] + (1.0 - e[i]) * (1.0 - mutual) * ambient);
    let a = (1.0 - e[0]) * mutual; let b = (1.0 - e[1]) * mutual;
    let j = [(rhs[0] + a * rhs[1]) / (1.0 - a * b),
        (rhs[1] + b * rhs[0]) / (1.0 - a * b)];
    std::array::from_fn(|i| AREA * e[i]
        * (black[i] - mutual * j[1 - i] - (1.0 - mutual) * ambient))
}

fn endpoint(old: &[f64; 8], dt: f64, mutual: f64, bath: f64, q: f64, right_robin: f64)
    -> [f64; 8]
{
    let mut face = [old[0], old[4]];
    for _ in 0..128 {
        let power = heat(face, mutual, bath);
        face = std::array::from_fn(|i| old[4 * i] + dt
            * (q / RHO - (power[i] + if i == 1 { right_robin } else { 0.0 }) / 3.0));
    }
    std::array::from_fn(|i| if i % 4 == 3 { old[i] + dt * q / RHO } else { face[i / 4] })
}

fn admit(report: &EnthalpyStepSolution, expected: &[f64; 8], old: &[f64; 8],
    datum: f64, dt: f64, source: f64, right_robin: f64)
{
    assert_eq!(&report.temperature[..4], &[400.0; 4]);
    assert_eq!(&report.temperature[4..], &[350.0; 4]);
    for (i, &h) in expected.iter().enumerate() {
        close(report.specific_enthalpy_j_kg[i] - datum, h, 3e-8);
        close(report.liquid_mass_fraction[i], h / LATENT, 4e-11);
    }
    let radiation = report.pointwise_enclosure_radiation.as_ref().unwrap();
    let delta = expected.iter().zip(old).map(|(new, old)| new - old).sum::<f64>();
    close(report.stored_energy_change_j, delta, 1e-7);
    close(report.robin_out_w, right_robin, 1e-9);
    close(delta + dt * (radiation.outward_heat_w + right_robin - source / 3.0), 0.0, 1e-7);
    assert!(report.energy_residual_j.abs() <= config().energy_tolerance_j);
    close(radiation.system_energy_closure_w, 0.0, 1e-8);
    for (i, patch) in radiation.patches.iter().enumerate() {
        assert_eq!(patch.material_states.len(), 3);
        close(patch.material_mean_emissivity.unwrap(), 0.2 + 0.6 * expected[4 * i] / LATENT, 2e-10);
    }
}

#[test]
fn closed_history_evolves_distinct_latent_bodies_and_preserves_internal_energy() {
    let mesh = mesh(); let context = context();
    let cards = [material_card("left"), material_card("right")];
    let enclosure = closed(&mesh, &cards, &context);
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    for datum in [0.0, 1_048_576.0] {
        let curves = [phase_curve(&cards[0], 400.0, datum), phase_curve(&cards[1], 350.0, datum)];
        with_cx(|cx| {
            let boundary = history(cx, &mesh, false); let initial = boundary.at(cx, 0.0).unwrap();
            let base = PointwiseEnclosureRadiation::new(cx, &mesh, &initial, &enclosure,
                PointwiseEnclosureTolerance::default()).unwrap();
            let mut trajectory = MaterialEnclosureEnthalpyHistoryRun::new(cx,
                run(cx, &boundary, &conductivity, &curves, datum, zero(), 4), &base,
                Vec::from(bindings(&cards, &curves, &context)), None).unwrap();
            let selected = trajectory.reference_emissivities().iter()
                .map(|e| e.receipt().selected).collect::<Vec<_>>();
            let mut old = [500.0; 8];
            for i in 1..=4 {
                let expected = endpoint(&old, 0.2, 1.0, 300.0, 0.0, 0.0);
                let report = trajectory.advance_to(cx, f64::from(i) * 0.2, config()).unwrap();
                admit(&report, &expected, &old, datum, 0.2, 0.0, 0.0);
                let radiation = report.pointwise_enclosure_radiation.as_ref().unwrap();
                assert!(radiation.ambient.is_none());
                close(radiation.outward_heat_w, 0.0, 1e-8);
                for (p, patch) in radiation.patches.iter().enumerate() {
                    for state in &patch.material_states {
                        assert_eq!(state.emissivity.receipt().selected, selected[p]);
                        assert_eq!(state.emissivity.card_identity(), cards[p].content_hash());
                    }
                }
                old = expected;
            }
            assert!(old[0] < 500.0 && old[4] > 500.0);
            assert_eq!(trajectory.conduction().accepted_steps(), 4);
        });
    }
}

#[test]
fn open_history_merges_bath_source_and_robin_events_without_double_counting_heat() {
    let mesh = mesh(); let context = context();
    let cards = [material_card("left"), material_card("right")];
    let curves = [phase_curve(&cards[0], 400.0, 0.0), phase_curve(&cards[1], 350.0, 0.0)];
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    with_cx(|cx| {
        let boundary = history(cx, &mesh, true); let initial = boundary.at(cx, 0.0).unwrap();
        let enclosure = open(cx, &mesh, &cards, &context);
        let base = PointwiseEnclosureRadiation::new_open(cx, &mesh, &initial, &enclosure,
            350.0, PointwiseEnclosureTolerance::default()).unwrap();
        let source = HeatSourceHistory::Uniform(TimeSignal::Table { times: vec![0.0, 0.25],
            values: vec![0.0, 48.0], dims: VOLUMETRIC_SOURCE_DIMS, interp: Interp::Hold });
        let bath = TimeSignal::Table { times: vec![0.0, 0.125, 0.375],
            values: vec![350.0, 450.0, 300.0], dims: TEMPERATURE_DIMS, interp: Interp::Hold };
        let mut trajectory = MaterialEnclosureEnthalpyHistoryRun::new(cx,
            run(cx, &boundary, &conductivity, &curves, 0.0, source, 5), &base,
            Vec::from(bindings(&cards, &curves, &context)), Some(bath)).unwrap();
        let mut old = [500.0; 8]; let mut start = 0.0;
        for (end, bath, q, air) in [(0.125, 350.0, 0.0, 350.0), (0.25, 450.0, 0.0, 350.0),
            (0.3, 450.0, 48.0, 350.0), (0.375, 450.0, 48.0, 400.0), (0.5, 300.0, 48.0, 400.0)]
        {
            let dt = end - start; let robin = AREA * 5.0 * (350.0 - air);
            let expected = endpoint(&old, dt, 0.4, bath, q, robin);
            let step = trajectory.advance_event_step(cx, 0.5, 1.0, config()).unwrap();
            assert_eq!(step.start_time_s, start); assert_eq!(step.end_time_s, end);
            admit(&step.report, &expected, &old, 0.0, dt, q, robin);
            let radiation = step.report.pointwise_enclosure_radiation.as_ref().unwrap();
            let reservoir = radiation.ambient.as_ref().unwrap();
            assert_eq!(reservoir.temperature_k, bath);
            close(radiation.outward_heat_w + reservoir.outward_heat_w, 0.0, 1e-8);
            old = expected; start = end;
        }
        assert_eq!(trajectory.conduction().time_s(), 0.5);
    });
}

#[test]
fn latent_enclosure_adaptation_rejects_by_phase_and_publishes_both_half_step_accounts() {
    let mesh = mesh(); let context = context();
    let cards = [material_card("left"), material_card("right")];
    let curves = [phase_curve(&cards[0], 400.0, 0.0), phase_curve(&cards[1], 350.0, 0.0)];
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let enclosure = closed(&mesh, &cards, &context);
    with_cx(|cx| {
        let boundary = history(cx, &mesh, false); let initial = boundary.at(cx, 0.0).unwrap();
        let base = PointwiseEnclosureRadiation::new(cx, &mesh, &initial, &enclosure,
            PointwiseEnclosureTolerance::default()).unwrap();
        let mut trajectory = MaterialEnclosureEnthalpyHistoryRun::new(cx,
            run(cx, &boundary, &conductivity, &curves, 0.0, zero(), 2), &base,
            Vec::from(bindings(&cards, &curves, &context)), None).unwrap();
        let step = trajectory.advance_adaptive_step(cx, 1.0, 1.0, config(), AdaptiveHistoryConfig {
            absolute_temperature_k: 1e-9, relative_temperature: 0.0,
            min_step_s: 1e-4, max_step_s: 1.0, max_attempts: 12,
        }, PhaseErrorTolerance { absolute_enthalpy_j_kg: 1000.0, relative_enthalpy: 0.0,
            absolute_liquid_fraction: 1e-6 }).unwrap();
        assert!(step.rejected_attempts > 0); assert!(step.end_time_s < 1.0);
        assert!(step.normalized_error <= 1.0);
        let first_dt = step.midpoint_time_s - step.start_time_s;
        let second_dt = step.end_time_s - step.midpoint_time_s;
        let first = endpoint(&[500.0; 8], first_dt, 1.0, 300.0, 0.0, 0.0);
        let second = endpoint(&first, second_dt, 1.0, 300.0, 0.0, 0.0);
        admit(&step.first_half, &first, &[500.0; 8], 0.0, first_dt, 0.0, 0.0);
        admit(&step.second_half, &second, &first, 0.0, second_dt, 0.0, 0.0);
        assert_eq!(trajectory.conduction().accepted_steps(), 2);
        assert_eq!(trajectory.conduction().specific_enthalpy_j_kg(), step.second_half.specific_enthalpy_j_kg);
        assert!(trajectory.advance_to(cx, 1.0, config()).is_err());
        trajectory.extend_step_budget(3).unwrap();
        trajectory.advance_to(cx, step.end_time_s + 0.01, config()).unwrap();
    });
}

#[test]
fn enclosure_binding_and_invalid_reservoir_refuse_atomically_and_allow_resume() {
    let mesh = mesh(); let context = context();
    let cards = [material_card("left"), material_card("right")];
    let curves = [phase_curve(&cards[0], 400.0, 0.0), phase_curve(&cards[1], 350.0, 0.0)];
    let clones = curves.clone(); let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let boundary = with_cx(|cx| history(cx, &mesh, false));
    let initial = with_cx(|cx| boundary.at(cx, 0.0).unwrap());
    let enclosure = with_cx(|cx| open(cx, &mesh, &cards, &context));
    let base = with_cx(|cx| PointwiseEnclosureRadiation::new_open(cx, &mesh, &initial,
        &enclosure, 350.0, PointwiseEnclosureTolerance::default()).unwrap());
    with_cx(|cx| {
        assert!(MaterialEnclosureEnthalpyHistoryRun::new(cx,
            run(cx, &boundary, &conductivity, &curves, 0.0, zero(), 3), &base,
            Vec::from(bindings(&cards, &clones, &context)), Some(constant(350.0))).is_err());
        assert!(MaterialEnclosureEnthalpyHistoryRun::new(cx,
            run(cx, &boundary, &conductivity, &curves, 0.0, zero(), 3), &base,
            Vec::from(bindings(&cards, &curves, &context)), None).is_err());
    });
    let mut trajectory = with_cx(|cx| MaterialEnclosureEnthalpyHistoryRun::new(cx,
        run(cx, &boundary, &conductivity, &curves, 0.0, zero(), 3), &base,
        Vec::from(bindings(&cards, &curves, &context)), Some(TimeSignal::Table {
            times: vec![0.0, 0.2], values: vec![350.0, -10.0], dims: TEMPERATURE_DIMS,
            interp: Interp::Hold })).unwrap());
    with_cx(|cx| { trajectory.advance_to(cx, 0.1, config()).unwrap(); });
    let h = trajectory.conduction().specific_enthalpy_j_kg().to_vec();
    let t = trajectory.conduction().temperature().to_vec();
    let f = trajectory.conduction().liquid_mass_fraction().to_vec();
    with_cancelled_cx(|cx| assert!(trajectory.advance_to(cx, 0.15, config()).is_err()));
    with_cx(|cx| {
        assert!(trajectory.advance_to(cx, 0.3, config()).is_err());
        assert_eq!(trajectory.conduction().specific_enthalpy_j_kg(), h);
        assert_eq!(trajectory.conduction().temperature(), t);
        assert_eq!(trajectory.conduction().liquid_mass_fraction(), f);
        assert_eq!(trajectory.conduction().time_s(), 0.1);
        assert_eq!(trajectory.conduction().accepted_steps(), 1);
        trajectory.advance_to(cx, 0.15, config()).unwrap();
        assert_eq!(trajectory.conduction().time_s(), 0.15);
    });
}

#[test]
fn original_conductivity_and_emissivity_sources_compose_with_active_spatial_transport() {
    use fs_conduction::bc::history::march::enthalpy::PhaseEnthalpyHistoryProblem;
    use fs_conduction::transient::enthalpy::material::phase::MaterialPhaseTransport;
    use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
    use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
    use fs_matdb::{ClaimSelection, SelectionPolicy};
    use fs_qty::{QuantitySpec, Temperature};

    let mesh = mesh();
    let context = context().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 400.0).unwrap()
        .with_quantity("melt", QuantitySpec::dimensional(Dims::NONE), 0.5).unwrap();
    let cards = ["left", "right"].map(|name| {
        let original = material_card(name); let mut claims = ClaimSet::new();
        let validity = ValidityDomain::unconstrained()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0, 500.0)
            .with_quantity("melt", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
            .with("lot", 7.0, 7.0);
        let mut emissivity = original.claims_for(SURFACE_EMISSIVITY_PROPERTY)[0].1.clone();
        emissivity.validity = validity.clone();
        claims.insert_claim(emissivity).unwrap();
        let k = if name == "left" { [2.0, 8.0] } else { [3.0, 11.0] };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new("k", KD),
            value: PropertyValue::Curve { abscissa: "melt".into(), abscissa_dims: Dims::NONE,
                knots: vec![(0.0, k[0]), (1.0, k[1])], dims: KD },
            validity,
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(), provenance: Provenance {
                source: "manufactured coupled conduction and reflection fixture".into(),
                license: "internal-test-use".into(), artifact: None },
        }).unwrap();
        MaterialCard::assemble(original.id().clone(), claims, Vec::new()).unwrap()
    });
    let curves = [phase_curve(&cards[0], 400.0, 0.0), phase_curve(&cards[1], 350.0, 0.0)];
    let sources = [0, 1].map(|i| PreparedPhaseConductivity::new(&curves[i], binding(&cards[i], &context),
        &PropertyKey::new("k", KD), ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [-100.0, LATENT + 200.0], || Ok(())).unwrap());
    let old = [-30.0, 400.0, 1100.0, 200.0, -60.0, 250.0, 1050.0, 800.0];
    let dt = 0.02;
    with_cx(|cx| {
        let materials = [0, 1].map(|i| MaterialPhaseTransport::try_new_isotropic(cx, &sources[i], RHO).unwrap());
        let boundary = history(cx, &mesh, false); let initial = boundary.at(cx, 0.0).unwrap();
        let enclosure = closed(&mesh, &cards, &context);
        let base = PointwiseEnclosureRadiation::new(cx, &mesh, &initial, &enclosure,
            PointwiseEnclosureTolerance::default()).unwrap();
        let declarations = bindings(&cards, &curves, &context);
        let radiation = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
        let direct = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
            cx, &mesh, &materials, &[0, 1], budget()).unwrap();
        let flux = direct.phase_transport_heat_flux(cx, &old).unwrap();
        assert!(flux.iter().all(|v| v.iter().map(|x| x * x).sum::<f64>() > 1.0));
        let zero_source = ScalarField::Uniform(0.0);
        let expected = direct.advance_with_material_enclosure_radiation(cx,
            direct.phase_transport_problem(&initial, &zero_source).unwrap(),
            None, &old, dt, config(), &radiation).unwrap();
        let caloric = EnthalpyHistoryRun::from_phase_materials(cx,
            PhaseEnthalpyHistoryProblem { boundary: &boundary, interfaces: None, source: zero() },
            &materials, &[0, 1], &old, EnthalpyHistoryBudget { spatial: budget(), max_steps: 1 }).unwrap();
        let mut trajectory = MaterialEnclosureEnthalpyHistoryRun::new(cx, caloric, &base,
            Vec::from(declarations), None).unwrap();
        let actual = trajectory.advance_to(cx, dt, config()).unwrap();
        assert_eq!(actual.specific_enthalpy_j_kg, expected.specific_enthalpy_j_kg);
        assert_eq!(actual.temperature, expected.temperature);
        assert_eq!(actual.liquid_mass_fraction, expected.liquid_mass_fraction);
        assert_eq!(actual.energy_residual_j, expected.energy_residual_j);
        // The apex has no radiation trace. Its heat change proves bulk transport.
        assert!((actual.specific_enthalpy_j_kg[3] - old[3]).abs() > 1e-4);
        close(actual.stored_energy_change_j, 0.0, config().energy_tolerance_j);
        close(actual.pointwise_enclosure_radiation.as_ref().unwrap().outward_heat_w, 0.0, 1e-8);
        for i in 0..2 {
            assert!(std::ptr::eq(trajectory.conduction().storage().phase_conductivity_source(i).unwrap(), &sources[i]));
        }
    });
}
