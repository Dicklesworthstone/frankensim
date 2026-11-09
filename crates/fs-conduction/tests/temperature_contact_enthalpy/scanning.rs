//! Moving absorbed power closes the original nonlinear thermal-contact model.
use super::*;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{
    EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun,
    adaptive::PhaseErrorTolerance,
    scanning::{ScanRadiation, ScanStep, ScanTransport, ScanningEnthalpyRun},
};
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::ScanSensitivityBudget;
use fs_conduction::bc::history::march::moving::sensitivity::MovingHeatGradient;
use fs_conduction::LinearConfig;
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyStep;
use fs_conduction::{AmbientRadiationConfig, SurfaceEmissivity, ThermalBc,
    COMPONENT_POWER_DIMS, TEMPERATURE_DIMS, STEFAN_BOLTZMANN_W_M2_K4};
use fs_qty::{Dims, QtyAny};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const OLD: [f64; 8] = [80., 80., 80., 80., 450., 450., 450., 450.];

fn constant(value: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Constant(QtyAny::new(value, dims))
}
fn heater(z: TimeSignal, power: TimeSignal) -> MovingHeatSource {
    heater_with_geometry(
        [constant(0.1, LENGTH), constant(0.1, LENGTH), z], [1.6, 1.6, 0.9], power,
    )
}
fn heater_with_geometry(center: [TimeSignal; 3], half_width: [f64; 3],
    power: TimeSignal) -> MovingHeatSource
{
    MovingHeatSource::new(center, half_width, power, 0.,
        MovingHeatBudget { max_vertices: 8, max_elements: 2,
            min_active_vertices: 1, power_relative_tolerance: 1e-12 },
    ).unwrap()
}
fn history<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh).adiabatic_remainder(),
        &Environment::earth_lab(), Vec::new(), 0.).unwrap()
}
fn base<'h, 'm, 'c>(cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>,
    material: &'h ConductivityModel, curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    max_steps: usize) -> EnthalpyHistoryRun<'h, 'm, 'c>
{
    EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
        boundary: history, material, element_materials: None, interfaces: None,
        source: HeatSourceHistory::Fixed(ScalarField::Uniform(2.)),
    }, EnthalpyHistoryMaterials {
        table: &[
            ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: 1. },
            ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: 2. },
        ], of_element: &[0, 1], phase_conductivity: None,
    }, &OLD, EnthalpyHistoryBudget {
        spatial: EnthalpyBudget { max_vertices: 8, max_elements: 2 }, max_steps,
    }).unwrap()
}
fn coupled(step: &ScanStep) -> &ContactEnthalpyStep {
    let ScanTransport::Contact(report) = &step.transport else {
        panic!("scan omitted its original source-backed contact owner");
    };
    report
}
fn stored(old: &[f64], new: &[f64]) -> f64 {
    old.iter().zip(new).enumerate()
        .map(|(i, (&before, &after))| (after - before) * if i < 4 { 1. / 24. } else { 1. / 12. })
        .sum()
}

#[test]
fn g1_scan_deposition_reaches_original_nonlinear_contacts_and_closes_total_energy() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let history = history(cx, &mesh); let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280., 1.6), (380., 3.6)]).unwrap());
        let initial = history.at(cx, 0.).unwrap();
        let source = heater(TimeSignal::Ramp { t_start: 0., t_end: 0.2,
            from: QtyAny::new(-0.65, LENGTH), to: QtyAny::new(0.65, LENGTH) },
            constant(3., COMPONENT_POWER_DIMS));
        let contacts = TemperatureContacts::new(cx, &mesh, &initial, vec![],
            vec![law(&mesh, &card)], 1).unwrap();
        let initial_resistance = contacts.evaluate(cx, &initial,
            &[308., 308., 308., 308., 345., 345., 345., 345.]).unwrap()
            .states()[0].resistance.value_m2_k_per_w();
        let mut scan = ScanningEnthalpyRun::new(cx, base(cx, &history, &material, &curves, 3),
            vec![source.clone()], 1, None).unwrap()
            .with_contacts(cx, contacts, coupling(), 1).unwrap();
        assert!(scan.contacts().is_some());
        assert!(scan.contact_coupling().is_some());
        assert!(scan.contact_newton_config().is_none());
        let mut old = OLD.to_vec();
        let mut from = 0.;
        let mut external = 0.;
        let mut final_resistance = initial_resistance;
        for end in [0.05, 0.15, 0.25] {
            let sample = source.at(cx, &mesh, end).unwrap();
            let combined = ScalarField::Nodal((0..8).map(|v| 2. + sample.source.at(v)).collect());
            let boundary = history.at(cx, end).unwrap();
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &combined };
            // An independent call to the existing contact owner verifies that
            // the scan supplied its actual P1 field and unchanged old history.
            let reference = scan.contacts().unwrap().advance_enthalpy(cx,
                scan.conduction().storage(), problem, &old, end - from, config(), coupling()).unwrap();
            let step = scan.advance_to(cx, end, config()).unwrap();
            let report = coupled(&step);
            for (&actual, &expected) in report.conduction.specific_enthalpy_j_kg.iter()
                .zip(&reference.conduction.specific_enthalpy_j_kg)
            {
                close(actual, expected, 1e-9);
            }
            close(step.deposition[0].center_m[2], sample.center_m[2], 0.);
            close(step.deposition[0].delivered_power_w, 3., 1e-12);
            let watts = 3. + 2. / 3.; // Two volume-1/6 tetrahedra plus one 3 W heater.
            close(report.physical.source_w, watts, 1e-11);
            close(report.physical.stored_energy_change_j, (end - from) * watts, 2e-8);
            close(stored(&old, &report.conduction.specific_enthalpy_j_kg),
                (end - from) * watts, 2e-8);
            close(report.physical.energy_residual_j, 0., config().energy_tolerance_j);
            assert!(report.physical.residual_norm_j <= coupling().residual_tolerance_j);
            let contact = &report.physical_contacts.states()[0];
            final_resistance = contact.resistance.value_m2_k_per_w();
            close(final_resistance, 0.1 + 0.006 * (contact.interface_temperature_k - 280.), 1e-14);
            assert_eq!(contact.resistance.card_identity(), card.content_hash());
            external += (end - from) * watts;
            old = report.conduction.specific_enthalpy_j_kg.clone();
            from = end;
        }
        assert!((final_resistance - initial_resistance).abs() > 1e-6);
        close(stored(&OLD, &old), external, 6e-8);
        assert_eq!(scan.conduction().specific_enthalpy_j_kg(), old.as_slice());
        assert_eq!(scan.conduction().accepted_steps(), 3);
    });
}

#[test]
fn g3_contact_scan_preserves_pulse_sides_and_both_adaptive_half_step_accounts() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let history = history(cx, &mesh); let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let initial = history.at(cx, 0.).unwrap();
        let source = heater(
            TimeSignal::Table { times: vec![0., 0.1, 0.2], values: vec![-0.65, 0.65, -0.65],
                dims: LENGTH, interp: Interp::Hold },
            TimeSignal::Table { times: vec![0., 0.1, 0.2], values: vec![4., 0., 2.],
                dims: COMPONENT_POWER_DIMS, interp: Interp::Hold });
        let create = || {
            let contacts = TemperatureContacts::new(cx, &mesh, &initial, vec![],
                vec![law(&mesh, &card)], 1).unwrap();
            ScanningEnthalpyRun::new(cx, base(cx, &history, &material, &curves, 4),
                vec![source.clone()], 1, None).unwrap()
                .with_contacts_newton(cx, contacts, coupled_newton()).unwrap()
        };
        let mut adaptive = create(); let mut manual = create();
        assert!(adaptive.contact_coupling().is_none());
        assert_eq!(adaptive.contact_newton_config().unwrap().step, coupled_newton().step);
        for (end, watts, z) in [(0.1, 4., -0.65), (0.2, 0., 0.65)] {
            let step = adaptive.advance_event_step(cx, 0.2, 1., coupled_newton().step).unwrap();
            let reference = manual.advance_event_step(cx, 0.2, 1., coupled_newton().step).unwrap();
            close(step.end_time_s, end, 0.);
            close(step.report.deposition[0].absorbed_power_w, watts, 0.);
            close(step.report.deposition[0].center_m[2], z, 0.);
            close(coupled(&step.report).physical.stored_energy_change_j, 0.1 * (watts + 2. / 3.), 2e-8);
            assert_eq!(step.report.solution().specific_enthalpy_j_kg,
                reference.report.solution().specific_enthalpy_j_kg);
        }
        let old = adaptive.conduction().specific_enthalpy_j_kg().to_vec();
        let pair = adaptive.advance_adaptive_step(cx, 0.25, 0.05, coupled_newton().step,
            AdaptiveHistoryConfig { absolute_temperature_k: 1000., relative_temperature: 0.,
                min_step_s: 1e-8, max_step_s: 0.05, max_attempts: 4 },
            PhaseErrorTolerance { absolute_enthalpy_j_kg: 1e6, relative_enthalpy: 0.,
                absolute_liquid_fraction: 1. }).unwrap();
        assert_eq!(pair.rejected_attempts, 0);
        for (end, half) in [(pair.midpoint_time_s, &pair.first_half), (pair.end_time_s, &pair.second_half)] {
            let reference = manual.advance_to(cx, end, coupled_newton().step).unwrap();
            assert_eq!(half.solution().specific_enthalpy_j_kg, reference.solution().specific_enthalpy_j_kg);
            close(coupled(half).physical.source_w, 2. + 2. / 3., 1e-11);
            close(coupled(half).physical.energy_residual_j, 0., coupled_newton().step.energy_tolerance_j);
        }
        close(stored(&old, adaptive.conduction().specific_enthalpy_j_kg()),
            coupled(&pair.first_half).physical.stored_energy_change_j
                + coupled(&pair.second_half).physical.stored_energy_change_j, 1e-10);
        close(stored(&old, adaptive.conduction().specific_enthalpy_j_kg()), 0.05 * (2. + 2. / 3.), 3e-8);
        assert_eq!(adaptive.conduction().accepted_steps(), 4);
        assert!(adaptive.advance_to(cx, 0.3, coupled_newton().step).is_err());
        adaptive.extend_step_budget(5).unwrap();
        assert!(matches!(adaptive.advance_to(cx, 0.3, coupled_newton().step).unwrap().transport,
            ScanTransport::Contact(_)));
    });
}

fn radiation() -> ScanRadiation {
    let dims = fs_conduction::EMISSIVITY_DIMS;
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::SURFACE_EMISSIVITY_PROPERTY, dims),
        value: PropertyValue::Scalar { value: 0.2, dims },
        validity: ValidityDomain::unconstrained().with("T", 280., 380.),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "synthetic scan/contact radiation".into(),
            license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    let card = fs_matdb::MaterialCard::assemble(MaterialStateId { chemistry: "scan surface".into(),
        phase: "solid".into(), process: "fixture".into(), revision: 0 }, claims, vec![]).unwrap();
    ScanRadiation { patches: vec![RadiationRegionHistory {
        region: "exterior".into(), emissivity: SurfaceEmissivity::from_card("exterior", &card,
            300., SelectionPolicy::SingleClaimOnly).unwrap(),
        ambient_temperature: constant(300., TEMPERATURE_DIMS),
    }], config: AmbientRadiationConfig::default() }
}

#[test]
fn g1_contact_scan_keeps_radiation_and_prescribed_reservoir_heat_distinct() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let card = card(); let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .region("bath", |face| face.centroid[2] > 0. && face.centroid[1] == 0.,
                ThermalBc::dirichlet(340.).unwrap()).unwrap()
            .region("exterior", |face| face.centroid[2] < 0. && face.centroid[1] == 0.,
                ThermalBc::neumann(0.).unwrap()).unwrap()
            .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.).unwrap();
        let initial = history.at(cx, 0.).unwrap();
        let make_scan = || ScanningEnthalpyRun::new(cx, base(cx, &history, &material, &curves, 2),
            vec![heater(constant(0.4, LENGTH), constant(1., COMPONENT_POWER_DIMS))],
            1, Some(radiation())).unwrap();
        let make_contacts = || TemperatureContacts::new(cx, &mesh, &initial, vec![],
            vec![law(&mesh, &card)], 1).unwrap();
        assert!(make_scan().with_contacts_newton(cx, make_contacts(), coupled_newton()).is_err());
        let mut scan = make_scan().with_contacts(cx, make_contacts(), coupling(), 1).unwrap();
        let step = scan.advance_to(cx, 0.05, config()).unwrap();
        let report = coupled(&step); let t = &report.conduction.temperature;
        let mean = (t[0] + t[1] + t[3]) / 3.;
        let radiation_w = 0.2 * STEFAN_BOLTZMANN_W_M2_K4 * 0.5 * (mean.powi(4) - 300_f64.powi(4));
        assert!(report.radiation.is_some() && report.combined_boundary.is_some());
        for vertex in [4, 5, 7] { close(t[vertex], 340., 1e-12); }
        close(report.physical.radiation_out_w, radiation_w, 1e-8);
        close(report.physical.convective_out_w, 0., 1e-12);
        close(report.physical.source_w, 1. + 2. / 3., 1e-11);
        assert!(report.physical.dirichlet_in_w.abs() > 1.);
        let expected = 0.05 * (1. + 2. / 3. + report.physical.dirichlet_in_w - radiation_w);
        close(stored(&OLD, &report.conduction.specific_enthalpy_j_kg), expected, 3e-8);
        close(report.physical.stored_energy_change_j, expected, 3e-8);
        close(report.physical.energy_residual_j, 0., config().energy_tolerance_j);
        assert!(report.physical.residual_norm_j <= coupling().residual_tolerance_j);
        assert_eq!(report.physical_contacts.states()[0].resistance.card_identity(), card.content_hash());
    });
}

#[test]
fn g4_contact_scan_policy_source_and_cancellation_refusals_preserve_history() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let history = history(cx, &mesh); let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let initial = history.at(cx, 0.).unwrap();
        let source = heater(constant(0.4, LENGTH), TimeSignal::Table {
            times: vec![0., 0.1], values: vec![4., 1e9], dims: COMPONENT_POWER_DIMS, interp: Interp::Hold,
        });
        let contacts = TemperatureContacts::new(cx, &mesh, &initial, vec![],
            vec![law(&mesh, &card)], 1).unwrap();
        let mut scan = ScanningEnthalpyRun::new(cx, base(cx, &history, &material, &curves, 3),
            vec![source], 1, None).unwrap()
            .with_contacts_newton(cx, contacts, coupled_newton()).unwrap();
        assert!(scan.advance_to(cx, 0.05, config()).is_err());
        let cancelled = CancelGate::new_clock_free(); cancelled.request();
        context(&cancelled, |cancelled_cx| {
            assert!(scan.advance_to(cancelled_cx, 0.05, coupled_newton().step).is_err());
        });
        assert_eq!(scan.conduction().specific_enthalpy_j_kg(), &OLD);
        assert_eq!(scan.conduction().time_s(), 0.);
        assert_eq!(scan.conduction().accepted_steps(), 0);
        scan.advance_to(cx, 0.05, coupled_newton().step).unwrap();
        let before = scan.conduction().specific_enthalpy_j_kg().to_vec();
        let temperature = scan.conduction().temperature().to_vec();
        let phase = scan.conduction().liquid_mass_fraction().to_vec();
        assert!(scan.advance_to(cx, 0.2, coupled_newton().step).is_err());
        assert_eq!(scan.conduction().specific_enthalpy_j_kg(), before.as_slice());
        assert_eq!(scan.conduction().temperature(), temperature.as_slice());
        assert_eq!(scan.conduction().liquid_mass_fraction(), phase.as_slice());
        assert_eq!(scan.conduction().time_s(), 0.05);
        assert_eq!(scan.conduction().accepted_steps(), 1);
        // A refused future right-side pulse did not poison the exact end-left
        // continuation to that switch or consume either remaining step slot.
        let resumed = scan.advance_event_step(cx, 0.1, 1., coupled_newton().step).unwrap();
        close(resumed.report.deposition[0].absorbed_power_w, 4., 0.);
        close(resumed.end_time_s, 0.1, 0.);
        assert_eq!(scan.conduction().accepted_steps(), 2);
    });
}

fn control_history<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .region("exterior", |face| face.centroid[2] < 0. && face.centroid[1] == 0.,
            ThermalBc::neumann(0.).unwrap()).unwrap()
        .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.).unwrap()
}
fn control_scan<'h, 'm, 'c>(cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>,
    material: &'h ConductivityModel, curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    card: &'c InterfaceSystemCard, p: [f64; 7], radiative: bool)
    -> ScanningEnthalpyRun<'h, 'm, 'c>
{
    let initial = history.at(cx, 0.).unwrap();
    let contacts = TemperatureContacts::new(cx, history.mesh(), &initial, vec![],
        vec![law(history.mesh(), card)], 1).unwrap();
    let source = heater_with_geometry(
        [p[0], p[1], p[2]].map(|v| constant(v, LENGTH)), [p[3], p[4], p[5]],
        constant(p[6], COMPONENT_POWER_DIMS));
    let scan = ScanningEnthalpyRun::new(cx, base(cx, history, material, curves, 2),
        vec![source], 1, radiative.then(radiation)).unwrap();
    if radiative {
        // Area-mean radiation uses the original coupled Picard owner. Keep its
        // complete physical residual well below the finite-difference scale.
        scan.with_contacts(cx, contacts, ContactIterationConfig {
            residual_tolerance_j: 1e-10, ..coupling()
        }, 1).unwrap()
    } else {
        scan.with_contacts_newton(cx, contacts, coupled_newton()).unwrap()
    }
}
fn control_step(radiative: bool) -> EnthalpyStepConfig {
    // Mean-radiation inner solves require positive relative tolerance; coupled
    // contact Newton instead retains its exact zero-relative-tolerance policy.
    if radiative { config() } else { coupled_newton().step }
}
fn control_budget() -> ScanSensitivityBudget {
    ScanSensitivityBudget { max_source_vertices: 8, max_radiation_feedback_entries: 1024 }
}
fn control_linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-10, max_iterations: 512, restart: 8 }
}
fn heater_controls(g: &MovingHeatGradient) -> [f64; 7] {
    [g.center_m[0], g.center_m[1], g.center_m[2],
        g.half_width_m[0], g.half_width_m[1], g.half_width_m[2], g.absorbed_power_w]
}
fn objective(seed: &[f64], state: &[f64]) -> f64 {
    seed.iter().zip(state).map(|(s, value)| s * value).sum()
}
fn close_control(actual: f64, expected: f64) {
    close(actual, expected, 5e-5 * (1. + actual.abs().max(expected.abs())));
}

#[test]
fn g3_contact_scan_controls_match_fresh_laws_with_and_without_mean_radiation() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let history = control_history(cx, &mesh); let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280., 1.6), (380., 3.6)]).unwrap());
        let p = [0.15, 0.2, 0.3, 1.6, 1.5, 1.2, 3.];
        let seed = [0.3, -0.1, 0.7, 0.5, -0.2, 0.4, 0.1, -0.3];
        for radiative in [false, true] {
            let step_config = control_step(radiative);
            let mut scan = control_scan(cx, &history, &material, &curves, &card, p, radiative);
            assert_eq!(scan.contact_coupling().is_some(), radiative);
            assert_eq!(scan.contact_newton_config().is_some(), !radiative);
            let endpoint = scan.advance_linearized_to(cx, 0.025, step_config,
                control_budget()).unwrap();
            let gradient = endpoint.pullback_temperature(cx, &seed, control_linear()).unwrap();
            assert_eq!(gradient.sources.len(), 1);
            assert_eq!(gradient.radiation.is_some(), radiative);
            assert_eq!(gradient.contact_resistance_multipliers.len(), 1);
            assert_eq!(gradient.contact_resistance_multipliers[0].interface, "joint");
            assert!(gradient.storage.transport.relative_residual <= control_linear().tolerance);
            assert!(gradient.sources[0].absorbed_power_w.abs() > 1e-5);
            assert!(gradient.contact_resistance_multipliers[0].derivative.abs() > 1e-3);
            assert_eq!(coupled(endpoint.report()).physical_contacts.states()[0]
                .resistance.card_identity(), card.content_hash());
            if let Some(radiation) = &gradient.radiation {
                assert_eq!(radiation.radiation_regions, vec!["exterior".to_string()]);
                assert_eq!(radiation.reservoir_temperatures.len(), 1);
                assert_eq!(radiation.emissivities.len(), 1);
            }
            // Binding the adjoint must preserve the accepted native endpoint,
            // not rerun a different frozen-contact or radiation approximation.
            let mut direct = control_scan(cx, &history, &material, &curves, &card, p, radiative);
            let primal = direct.advance_to(cx, 0.025, step_config).unwrap();
            assert_eq!(endpoint.report().solution().specific_enthalpy_j_kg,
                primal.solution().specific_enthalpy_j_kg);
            let evaluate = |controls: [f64; 7], resistance_scale: f64| {
                // Scale the entire original R(T) claim, including its slope.
                // Each finite difference re-queries live contact and radiation.
                let trial_card = card_with_curve(vec![
                    (280., 0.1 * resistance_scale), (380., 0.7 * resistance_scale)]);
                let mut trial = control_scan(cx, &history, &material, &curves,
                    &trial_card, controls, radiative);
                let step = trial.advance_to(cx, 0.025, step_config).unwrap();
                objective(&seed, &step.solution().temperature)
            };
            let actual = heater_controls(&gradient.sources[0]);
            let epsilon = 1e-4;
            for i in 0..7 {
                let (mut plus, mut minus) = (p, p);
                plus[i] += epsilon; minus[i] -= epsilon;
                close_control(actual[i],
                    (evaluate(plus, 1.) - evaluate(minus, 1.)) / (2. * epsilon));
            }
            close_control(gradient.contact_resistance_multipliers[0].derivative,
                (evaluate(p, 1. + epsilon) - evaluate(p, 1. - epsilon)) / (2. * epsilon));
        }
    });
}

#[test]
fn g3_contact_scan_reverse_sweep_carries_earlier_heat_and_resistance_through_melting() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let history = control_history(cx, &mesh); let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280., 1.6), (380., 3.6)]).unwrap());
        let p = [0.15, 0.2, 0.3, 1.6, 1.5, 1.2, 3.];
        let seed = [0.3, -0.1, 0.7, 0.5, -0.2, 0.4, 0.1, -0.3];
        let step_config = control_step(true);
        let mut scan = control_scan(cx, &history, &material, &curves, &card, p, true);
        let first = scan.advance_linearized_to(cx, 0.03, step_config,
            control_budget()).unwrap();
        let last = scan.advance_linearized_to(cx, 0.08, step_config,
            control_budget()).unwrap();
        assert!(first.report().solution().liquid_mass_fraction.iter().all(|&f| f == 0.));
        assert!(last.report().solution().liquid_mass_fraction[..3].iter()
            .any(|&f| f > 0.02 && f < 0.98));
        // The terminal objective observes h, so its earlier-energy dependence
        // survives even when the final temperature is on a latent plateau.
        let last_gradient = last.pullback(cx, &seed, control_linear()).unwrap();
        let first_gradient = first.pullback(cx,
            &last_gradient.storage.transport.previous_specific_enthalpy, control_linear()).unwrap();
        assert!(first_gradient.sources[0].absorbed_power_w.abs() > 1e-5);
        assert!(first_gradient.contact_resistance_multipliers[0].derivative.abs() > 1e-3);
        assert_eq!(first_gradient.contact_resistance_multipliers[0].interface, "joint");
        assert_eq!(last_gradient.contact_resistance_multipliers[0].interface, "joint");
        assert!(first_gradient.radiation.is_some() && last_gradient.radiation.is_some());
        let evaluate = |controls: [f64; 7], resistance_scale: f64| {
            let trial_card = card_with_curve(vec![
                (280., 0.1 * resistance_scale), (380., 0.7 * resistance_scale)]);
            let mut trial = control_scan(cx, &history, &material, &curves,
                &trial_card, controls, true);
            trial.advance_to(cx, 0.03, step_config).unwrap();
            let step = trial.advance_to(cx, 0.08, step_config).unwrap();
            objective(&seed, &step.solution().specific_enthalpy_j_kg)
        };
        let first_controls = heater_controls(&first_gradient.sources[0]);
        let last_controls = heater_controls(&last_gradient.sources[0]);
        let epsilon = 1e-4;
        // The authored constant center and power affect BOTH fixed intervals.
        for i in [2, 6] {
            let (mut plus, mut minus) = (p, p);
            plus[i] += epsilon; minus[i] -= epsilon;
            close_control(first_controls[i] + last_controls[i],
                (evaluate(plus, 1.) - evaluate(minus, 1.)) / (2. * epsilon));
        }
        close_control(first_gradient.contact_resistance_multipliers[0].derivative
                + last_gradient.contact_resistance_multipliers[0].derivative,
            (evaluate(p, 1. + epsilon) - evaluate(p, 1. - epsilon)) / (2. * epsilon));
    });
}
