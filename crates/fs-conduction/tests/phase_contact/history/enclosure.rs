//! Independent two-body balance with latent contact and spatial gray reflection.
use super::*;
use fs_conduction::radiation::OpenDiffuseEnclosure;
use fs_evidence::ValidityDomain;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::QtyAny;

const AREA: f64 = 0.5;
const RADIATING: [[usize; 3]; 2] = [[0, 1, 3], [4, 5, 7]];
const REGIONS: [&str; 2] = ["cold-radiation", "hot-radiation"];

fn material_card(name: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "melt".into(), abscissa_dims: Dims::NONE,
            knots: vec![(0.0, 0.2), (1.0, 0.8)], dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0, 800.0)
            .with_quantity("melt", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
            .with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(), provenance: Provenance {
            source: "manufactured latent contact and gray radiation".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "chart".into(),
        process: "combined heat fixture".into(), revision: 0 }, claims, Vec::new()).unwrap()
}
fn context() -> QueryPoint {
    point().with_quantity("melt", QuantitySpec::dimensional(Dims::NONE), 0.4).unwrap()
}
fn material<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> ThermalMaterialPoint<'a> {
    ThermalMaterialPoint { material: card, context, temperature_axis: "T", liquid_fraction_axis: Some("melt") }
}
fn patches<'a>(cards: &'a [MaterialCard; 2], curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    context: &'a QueryPoint) -> Vec<MaterialEnclosureRadiationPatch<'a>>
{
    [0, 1].map(|i| MaterialEnclosureRadiationPatch { region: REGIONS[i], curve: &curves[i],
        material: material(&cards[i], context), selection: MaterialPropertySelection::SingleClaimOnly }).into_iter().collect()
}
fn history<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .region(REGIONS[0], |face| face.centroid[1] == 0.0 && face.vertices[0] < 4,
            ThermalBc::adiabatic()).unwrap()
        .region(REGIONS[1], |face| face.centroid[1] == 0.0 && face.vertices[0] >= 4,
            ThermalBc::adiabatic()).unwrap().adiabatic_remainder(),
        &Environment::earth_lab(), Vec::new(), 0.0).unwrap()
}
fn surfaces(mesh: &ConductionMesh, cards: &[MaterialCard; 2], context: &QueryPoint) -> Vec<RadiationSurface> {
    [0, 1].map(|body| {
        let snapshot = SurfaceEmissivity::from_condition(REGIONS[body], material(&cards[body], context),
            MaterialPropertySelection::SingleClaimOnly).unwrap().0;
        RadiationSurface::new(mesh, REGIONS[body], |face| face.centroid[1] == 0.0
            && face.vertices[0] as usize / 4 == body, snapshot).unwrap()
    }).to_vec()
}
fn closed(mesh: &ConductionMesh, cards: &[MaterialCard; 2], context: &QueryPoint) -> GrayDiffuseEnclosure {
    // Manufactured reciprocal algebra: tetrahedra supply disjoint contact/radiation traces,
    // not a claim that these coplanar faces geometrically see one another.
    GrayDiffuseEnclosure::new(surfaces(mesh, cards, context),
        ViewFactorMatrix::infinite_parallel_plates(AREA).unwrap()).unwrap()
}
fn open(cx: &Cx<'_>, mesh: &ConductionMesh, cards: &[MaterialCard; 2], context: &QueryPoint) -> OpenDiffuseEnclosure {
    OpenDiffuseEnclosure::new(cx, surfaces(mesh, cards, context), vec![vec![0.0, 0.4], vec![0.4, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "manufactured reciprocal contact/radiation fixture".into() },
        ViewFactorTolerance::default()).unwrap()
}
fn bath(temperature_k: f64) -> TimeSignal {
    TimeSignal::Constant(QtyAny::new(temperature_k, TEMPERATURE_DIMS))
}

// Independent eight-coordinate BE equations. Each body stays on its latent plateau,
// so conduction vanishes while R(fa,fb) and nodal epsilon(f) evolve. Analytic triangle
// moments give integral(Ni*epsilon)=A*(sum(epsilon)+epsilon_i)/12. Two gray reflection
// equations are inverted explicitly, without production source queries or quadrature.
fn flux(h: &[f64; 8], mutual: f64, bath: f64) -> ([f64; 8], [f64; 2]) {
    let means: [f64; 2] = std::array::from_fn(|body|
        h[4 * body..4 * body + 3].iter().sum::<f64>() / 300.0);
    let contact = AREA * (400.0 - 500.0) / (1.0 + 2.0 * means[0] + means[1]);
    let epsilon = RADIATING.map(|nodes| nodes.map(|i| 0.2 + 0.006 * h[i]));
    let mean_e = epsilon.map(|e| e.iter().sum::<f64>() / 3.0);
    let black = [400.0_f64, 500.0].map(|t| STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4));
    let ambient = STEFAN_BOLTZMANN_W_M2_K4 * bath.powi(4);
    let rhs: [f64; 2] = std::array::from_fn(|i|
        mean_e[i] * black[i] + (1.0 - mean_e[i]) * (1.0 - mutual) * ambient);
    let a = (1.0 - mean_e[0]) * mutual; let b = (1.0 - mean_e[1]) * mutual;
    let j = [(rhs[0] + a * rhs[1]) / (1.0 - a * b), (rhs[1] + b * rhs[0]) / (1.0 - a * b)];
    let mut load = [0.0; 8]; let mut radiation = [0.0; 2];
    for body in 0..2 {
        for i in 0..3 { load[4 * body + i] = (if body == 0 { contact } else { -contact }) / 3.0; }
        let irradiation = mutual * j[1 - body] + (1.0 - mutual) * ambient;
        for (local, node) in RADIATING[body].into_iter().enumerate() {
            let q = AREA / 12.0 * (epsilon[body].iter().sum::<f64>() + epsilon[body][local])
                * (black[body] - irradiation);
            load[node] += q; radiation[body] += q;
        }
    }
    (load, radiation)
}
fn endpoint(old: &[f64; 8], dt: f64, mutual: f64, bath: f64, source: f64) -> [f64; 8] {
    let mut h = *old;
    for _ in 0..128 {
        let (outward, _) = flux(&h, mutual, bath);
        h = std::array::from_fn(|i| old[i] + dt * (source / 24.0 - outward[i])
            * if i < 4 { 24.0 } else { 12.0 });
    }
    assert!(h.iter().all(|h| *h > 0.0 && *h < 100.0));
    h
}
fn admit(report: &ContactEnthalpyStep, old: &[f64; 8], expected: &[f64; 8], dt: f64,
    mutual: f64, bath: f64, source: f64)
{
    for i in 0..8 {
        near(report.conduction.specific_enthalpy_j_kg[i], expected[i], 3e-8);
        near(report.conduction.liquid_mass_fraction[i], expected[i] / 100.0, 3e-10);
    }
    assert_eq!(report.conduction.temperature, [400.0, 400.0, 400.0, 400.0, 500.0, 500.0, 500.0, 500.0]);
    let radiation = report.conduction.pointwise_enclosure_radiation.as_ref().unwrap();
    let (_, heat) = flux(expected, mutual, bath);
    for (body, patch) in radiation.patches.iter().enumerate() {
        near(patch.outward_heat_w, heat[body], 2e-7);
        assert_eq!(patch.material_states.len(), 3);
        near(patch.material_mean_emissivity.unwrap(), RADIATING[body].iter()
            .map(|&i| 0.2 + 0.006 * expected[i]).sum::<f64>() / 3.0, 3e-10);
        for state in &patch.material_states {
            near(state.emissivity.value(), 0.2 + 0.006 * expected[state.vertex], 3e-10);
        }
    }
    near(radiation.outward_heat_w, heat.iter().sum(), 2e-7);
    near(radiation.system_energy_closure_w, 0.0, 1e-8);
    let stored: f64 = expected.iter().zip(old).enumerate()
        .map(|(i, (new, old))| (new - old) * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 }).sum();
    near(report.conduction.stored_energy_change_j, stored, 2e-9);
    near(stored + dt * (heat.iter().sum::<f64>() - source / 3.0), 0.0, 2e-9);
    near(report.physical.energy_residual_j, 0.0, config().step.energy_tolerance_j);
    near(report.physical.radiation_out_w, heat.iter().sum(), 2e-7);
    assert!(report.physical.residual_norm_j <= config().residual_tolerance_j);
    assert_eq!(report.iterations, 1); assert_eq!(report.radiation_iterations, 0);
    for trace in report.physical_contacts.material_trace_states() {
        for body in 0..2 { near(trace.liquid_mass_fractions[body],
            expected[4 * body..4 * body + 3].iter().sum::<f64>() / 300.0, 3e-10); }
    }
}

#[test]
fn coupled_closed_history_resolves_contact_and_spatial_reflection_at_the_same_enthalpy() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = history(cx, &mesh); let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let enclosure = closed(&mesh, &cards, &context);
        let radiation = PointwiseEnclosureRadiation::new(cx, &mesh, &initial, &enclosure,
            PointwiseEnclosureTolerance::default()).unwrap();
        let mut trajectory = ContactEnthalpyHistoryRun::new_newton(cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 4),
            contacts(cx, &mesh, &initial, &cards, &curves, &source, &context,
                sides(&mesh, 4, 2, 0.0), 0, false), config()).unwrap()
            .with_material_enclosure(cx, &radiation, patches(&cards, &curves, &context), None).unwrap();
        assert_eq!(trajectory.material_enclosure_patches().unwrap().len(), 2);
        assert!(trajectory.material_enclosure_ambient_temperature().is_none());
        let claims: Vec<_> = trajectory.material_enclosure_reference_emissivities().unwrap().iter()
            .map(|epsilon| epsilon.receipt().selected).collect();
        let mut old = OLD; let mut from = 0.0;
        for step in 1..=4 {
            let to = f64::from(step) * 0.0005; let dt = to - from;
            let expected = endpoint(&old, dt, 1.0, 300.0, 0.0);
            let report = trajectory.advance_to(cx, to, config().step).unwrap();
            admit(&report, &old, &expected, dt, 1.0, 300.0, 0.0);
            let actual = report.conduction.pointwise_enclosure_radiation.as_ref().unwrap();
            assert!(actual.ambient.is_none()); near(actual.outward_heat_w, 0.0, 1e-8);
            for (body, patch) in actual.patches.iter().enumerate() {
                for state in &patch.material_states {
                    assert_eq!(state.emissivity.receipt().selected, claims[body]);
                    assert_eq!(state.emissivity.card_identity(), cards[body].content_hash());
                }
            }
            old = expected; from = to;
        }
        assert!(old[0] > OLD[0] && old[4] < OLD[4]);
        assert_ne!(old[0], old[2]); // Radiation and contact act on different physical faces.
        assert_eq!(trajectory.conduction().accepted_steps(), 4);
    });
}

#[test]
fn coupled_open_history_merges_bath_and_heater_events_with_one_external_heat_account() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = history(cx, &mesh); let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let enclosure = open(cx, &mesh, &cards, &context);
        let radiation = PointwiseEnclosureRadiation::new_open(cx, &mesh, &initial, &enclosure,
            350.0, PointwiseEnclosureTolerance::default()).unwrap();
        let forcing = HeatSourceHistory::Uniform(TimeSignal::Table { times: vec![0.0, 0.0004],
            values: vec![0.0, 60.0], dims: VOLUMETRIC_SOURCE_DIMS, interp: Interp::Hold });
        let ambient = TimeSignal::Table { times: vec![0.0, 0.0002, 0.0004],
            values: vec![350.0, 600.0, 450.0], dims: TEMPERATURE_DIMS, interp: Interp::Hold };
        let mut trajectory = ContactEnthalpyHistoryRun::new_newton(cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, forcing, 3),
            contacts(cx, &mesh, &initial, &cards, &curves, &source, &context,
                sides(&mesh, 4, 2, 0.0), 0, false), config()).unwrap()
            .with_material_enclosure(cx, &radiation, patches(&cards, &curves, &context), Some(ambient)).unwrap();
        let mut old = OLD; let mut from = 0.0;
        for (to, bath, q) in [(0.0002, 350.0, 0.0), (0.0004, 600.0, 0.0), (0.0007, 450.0, 60.0)] {
            let dt = to - from; let expected = endpoint(&old, dt, 0.4, bath, q);
            let step = trajectory.advance_event_step(cx, 0.0007, 1.0, config().step).unwrap();
            assert_eq!(step.start_time_s, from); assert_eq!(step.end_time_s, to);
            admit(&step.report, &old, &expected, dt, 0.4, bath, q);
            let radiation = step.report.conduction.pointwise_enclosure_radiation.as_ref().unwrap();
            let reservoir = radiation.ambient.as_ref().unwrap();
            assert_eq!(reservoir.temperature_k, bath);
            near(radiation.outward_heat_w + reservoir.outward_heat_w, 0.0, 1e-8);
            old = expected; from = to;
        }
        assert_eq!(trajectory.conduction().time_s(), 0.0007);
        assert_eq!(trajectory.conduction().accepted_steps(), 3);
    });
}

#[test]
fn coupled_history_controls_latent_error_and_preserves_state_after_future_bath_failure() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = history(cx, &mesh); let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)]; let clones = curves.clone();
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let enclosure = open(cx, &mesh, &cards, &context);
        let radiation = PointwiseEnclosureRadiation::new_open(cx, &mesh, &initial, &enclosure,
            350.0, PointwiseEnclosureTolerance::default()).unwrap();
        let create = || ContactEnthalpyHistoryRun::new_newton(cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
            contacts(cx, &mesh, &initial, &cards, &curves, &source, &context,
                sides(&mesh, 4, 2, 0.0), 0, false), config()).unwrap();
        assert!(create().with_material_enclosure(cx, &radiation, patches(&cards, &curves, &context), None).is_err());
        assert!(create().with_material_enclosure(cx, &radiation, patches(&cards, &clones, &context), Some(bath(350.0))).is_err());
        let valid = || create().with_material_enclosure(cx, &radiation,
            patches(&cards, &curves, &context), Some(bath(350.0))).unwrap();
        assert!(valid().with_material_enclosure(cx, &radiation, patches(&cards, &curves, &context), Some(bath(350.0))).is_err());
        let mut adaptive = valid();
        let accepted = adaptive.advance_adaptive_step(cx, 0.001, 0.001, config().step,
            AdaptiveHistoryConfig { absolute_temperature_k: 1e-9, relative_temperature: 0.0,
                min_step_s: 1e-8, max_step_s: 0.001, max_attempts: 16 },
            PhaseErrorTolerance { absolute_enthalpy_j_kg: 100.0, relative_enthalpy: 0.0,
                absolute_liquid_fraction: 1e-7 }).unwrap();
        assert!(accepted.rejected_attempts > 0); assert!(accepted.end_time_s < 0.001);
        assert!(accepted.normalized_error <= 1.0);
        let first_dt = accepted.midpoint_time_s;
        let second_dt = accepted.end_time_s - accepted.midpoint_time_s;
        let first = endpoint(&OLD, first_dt, 0.4, 350.0, 0.0);
        let second = endpoint(&first, second_dt, 0.4, 350.0, 0.0);
        admit(&accepted.first_half, &OLD, &first, first_dt, 0.4, 350.0, 0.0);
        admit(&accepted.second_half, &first, &second, second_dt, 0.4, 350.0, 0.0);
        assert_eq!(adaptive.conduction().accepted_steps(), 2);
        assert_eq!(adaptive.conduction().specific_enthalpy_j_kg(), accepted.second_half.conduction.specific_enthalpy_j_kg);
        let mut invalid_bath = create().with_material_enclosure(cx, &radiation,
            patches(&cards, &curves, &context), Some(TimeSignal::Table {
                times: vec![0.0, 0.0004], values: vec![350.0, -10.0],
                dims: TEMPERATURE_DIMS, interp: Interp::Hold })).unwrap();
        invalid_bath.advance_to(cx, 0.0002, config().step).unwrap();
        let saved = invalid_bath.conduction().specific_enthalpy_j_kg().to_vec();
        assert!(invalid_bath.advance_to(cx, 0.0005, config().step).is_err());
        support::with_cancelled_cx(|cancelled|
            assert!(invalid_bath.advance_to(cancelled, 0.0003, config().step).is_err()));
        assert_eq!(invalid_bath.conduction().time_s(), 0.0002);
        assert_eq!(invalid_bath.conduction().accepted_steps(), 1);
        assert_eq!(invalid_bath.conduction().specific_enthalpy_j_kg(), saved);
        let resumed = invalid_bath.advance_to(cx, 0.0003, config().step).unwrap();
        let mut replay = valid(); replay.advance_to(cx, 0.0002, config().step).unwrap();
        let fresh = replay.advance_to(cx, 0.0003, config().step).unwrap();
        assert_eq!(resumed.conduction.specific_enthalpy_j_kg, fresh.conduction.specific_enthalpy_j_kg);
        assert_eq!(resumed.physical_contacts.material_trace_states(), fresh.physical_contacts.material_trace_states());
    });
}
