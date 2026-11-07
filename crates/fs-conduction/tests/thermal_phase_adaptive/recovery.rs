//! Physical trial recovery on the original sharp-phase spatial residual.
use super::*;
use fs_conduction::bc::history::ThermalRegionHistory;
use fs_conduction::bc::history::march::enthalpy::adaptive::EntropyAdaptiveHistoryError;
use fs_conduction::transient::enthalpy::EnthalpyError;
use fs_conduction::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyError;
use fs_conduction::transient::enthalpy::sdirk::SDIRK2_GAMMA;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-8,
        "actual={actual:.17e}, expected={expected:.17e}");
}

fn limited() -> EnthalpyStepConfig {
    EnthalpyStepConfig { max_newton_iterations: 2, ..controls() }
}

fn recovery_controls() -> AdaptiveHistoryConfig {
    // Deliberately isolate nonlinear recovery from temporal-error rejection.
    AdaptiveHistoryConfig { absolute_temperature_k: 100.0, relative_temperature: 0.0,
        min_step_s: 0.015625, max_step_s: 8.0, max_attempts: 8 }
}

fn phase_controls() -> PhaseErrorTolerance {
    PhaseErrorTolerance { absolute_enthalpy_j_kg: 1000.0,
        relative_enthalpy: 0.0, absolute_liquid_fraction: 1.0 }
}

fn bath_history<'m>(
    cx: &Cx<'_>, mesh: &'m ConductionMesh, bath: f64,
    rows: Vec<ThermalRegionHistory>,
) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .region("bath", |face| face.centroid[2] == 0.0,
            ThermalBc::dirichlet(bath).unwrap()).unwrap()
        .adiabatic_remainder(), &Environment::earth_lab(), rows, 0.0).unwrap()
}

fn bath_h(bath: f64) -> f64 {
    if bath < 310.0 { 10.0 * (bath - 300.0) } else { 10.0 * (bath - 260.0) }
}

fn create<'h, 'm, 'c>(
    cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>,
    material: &'h ConductivityModel, chart: &'c EquilibriumEnthalpyPhaseCurve,
    bath: f64, old: f64,
) -> EnthalpyHistoryRun<'h, 'm, 'c> {
    let fixed = bath_h(bath);
    EnthalpyHistoryRun::new(cx,
        problem(history, material, HeatSourceHistory::Fixed(ScalarField::Uniform(0.0))),
        EnthalpyHistoryMaterials {
            table: &[ReferenceEnthalpyMaterial { curve: chart, reference_density_kg_m3: 1.0 }],
            of_element: &[0], phase_conductivity: None,
        }, &[fixed, fixed, fixed, old],
        EnthalpyHistoryBudget {
            spatial: EnthalpyBudget { max_vertices: 4, max_elements: 1 }, max_steps: 4,
        }).unwrap()
}

fn untouched(run: &EnthalpyHistoryRun<'_, '_, '_>, bath: f64, old: f64) {
    let fixed = bath_h(bath);
    assert_eq!(run.specific_enthalpy_j_kg(), &[fixed, fixed, fixed, old]);
    assert_eq!(run.time_s(), 0.0);
    assert_eq!(run.accepted_steps(), 0);
}

fn stalled(error: &HeterogeneousEnthalpyError) -> bool {
    matches!(error, HeterogeneousEnthalpyError::Enthalpy(EnthalpyError::NotConverged(_)))
}

// Independent one-free-node equations: the apex stores m=1/24 kg and
// exchanges G=1/3 W/K with the three prescribed face nodes. Thus
// dh/dt=8*(Tb-T(h)), with Cp=10 and a 400 J/kg sharp latent plateau.
fn temperature(h: f64) -> f64 {
    if h < 100.0 { 300.0 + h / 10.0 }
    else if h <= 500.0 { 310.0 }
    else { 260.0 + h / 10.0 }
}

fn implicit(old: f64, dt: f64, bath: f64) -> f64 {
    let solid = (old + 8.0 * dt * (bath - 300.0)) / (1.0 + 0.8 * dt);
    if solid <= 100.0 { return solid; }
    let latent = old + 8.0 * dt * (bath - 310.0);
    if latent <= 500.0 { return latent; }
    (old + 8.0 * dt * (bath - 260.0)) / (1.0 + 0.8 * dt)
}

fn sdirk(old: f64, dt: f64, bath: f64) -> [f64; 2] {
    let first = implicit(old, SDIRK2_GAMMA * dt, bath);
    // The algebraic carry is not queried as a physical material state.
    let carry = dt * (1.0 - SDIRK2_GAMMA) * 8.0 * (bath - temperature(first));
    [first, implicit(old + carry, SDIRK2_GAMMA * dt, bath)]
}

#[test]
fn g1_backward_euler_recovers_sharp_melting_and_freezing_without_extra_heat() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let chart = curve(0.0);
    with_cx(|cx| {
        for (bath, old, proposed) in [(340.0, 0.0, 4.0), (290.0, 1000.0, 6.0)] {
            let history = bath_history(cx, &mesh, bath, Vec::new());
            let mut run = create(cx, &history, &material, &chart, bath, old);
            let refusal = run.advance_to(cx, proposed, limited()).unwrap_err();
            assert!(stalled(&refusal), "{refusal:?}");
            untouched(&run, bath, old);

            let selected = run.advance_adaptive_step(cx, proposed, proposed, limited(),
                recovery_controls(), phase_controls()).unwrap();
            assert_eq!(selected.end_time_s, proposed / 2.0);
            assert_eq!(selected.rejected_attempts, 1);
            // One failed coarse solve, then the complete accepted comparison.
            assert_eq!(selected.endpoint_solves, 4);
            assert_eq!(run.accepted_steps(), 2);
            let half_dt = proposed / 4.0;
            let first = implicit(old, half_dt, bath);
            let second = implicit(first, half_dt, bath);
            close(selected.first_half.specific_enthalpy_j_kg[3], first);
            close(selected.second_half.specific_enthalpy_j_kg[3], second);
            close(run.liquid_mass_fraction()[3], (second - 100.0) / 400.0);
            assert_eq!(run.temperature()[3], 310.0);
            let heat = selected.first_half.dirichlet_in_w * half_dt
                + selected.second_half.dirichlet_in_w * half_dt;
            close(heat, (second - old) / 24.0);
            close(selected.first_half.stored_energy_change_j
                + selected.second_half.stored_energy_change_j, heat);

            let mut audited = create(cx, &history, &material, &chart, bath, old);
            let checked = audited.advance_adaptive_step_with_entropy(cx, proposed, proposed,
                limited(), recovery_controls(), phase_controls(), 1.0).unwrap();
            assert_eq!(checked.physical.endpoint_solves, 4);
            assert_eq!(checked.physical.rejected_attempts, 1);
            assert_eq!(audited.specific_enthalpy_j_kg(), run.specific_enthalpy_j_kg());
            let fixed = bath_h(bath);
            close(checked.entropy.thermal_storage_change_j_k,
                audited.storage().thermal_entropy_change_j_k(cx,
                    &[fixed, fixed, fixed, old], audited.specific_enthalpy_j_kg()).unwrap());
            checked.entropy.admit(1.0).unwrap();

            // Resume from the accepted prefix, without charging its heat twice.
            let resumed = run.advance_adaptive_step(cx, proposed, proposed, limited(),
                recovery_controls(), phase_controls()).unwrap();
            assert!(resumed.end_time_s > selected.end_time_s
                && resumed.end_time_s <= proposed);
            assert_eq!(run.accepted_steps(), 4);
            let resumed_half_dt = (resumed.end_time_s - selected.end_time_s) / 2.0;
            let middle = implicit(second, resumed_half_dt, bath);
            let final_h = implicit(middle, resumed_half_dt, bath);
            close(run.specific_enthalpy_j_kg()[3], final_h);
            close(heat + resumed.first_half.dirichlet_in_w * resumed_half_dt
                + resumed.second_half.dirichlet_in_w * resumed_half_dt,
                (final_h - old) / 24.0);
        }
    });
}

#[test]
fn g1_sdirk_recovers_with_authentic_fine_stages_and_entropy() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let chart = curve(0.0);
    with_cx(|cx| {
        let history = bath_history(cx, &mesh, 340.0, Vec::new());
        let mut run = create(cx, &history, &material, &chart, 340.0, 0.0).into_sdirk2();
        let refusal = run.advance_to(cx, 8.0, limited()).unwrap_err();
        assert!(stalled(&refusal), "{refusal:?}");
        untouched(run.conduction(), 340.0, 0.0);
        let selected = run.advance_adaptive_step(cx, 8.0, 8.0, limited(),
            recovery_controls(), phase_controls()).unwrap();
        assert_eq!(selected.end_time_s, 4.0);
        assert_eq!(selected.rejected_attempts, 1);
        assert_eq!(selected.endpoint_solves, 4);
        assert_eq!(run.conduction().accepted_steps(), 2);
        let first = sdirk(0.0, 2.0, 340.0);
        let second = sdirk(first[1], 2.0, 340.0);
        for (physical, expected) in [
            (&selected.first_half, first), (&selected.second_half, second),
        ] {
            close(physical.first_stage.specific_enthalpy_j_kg[3], expected[0]);
            close(physical.second_stage.specific_enthalpy_j_kg[3], expected[1]);
            assert!(physical.energy_residual_j.abs() <= limited().energy_tolerance_j);
            assert!(physical.nodal_residual_j.iter()
                .all(|r| r.abs() <= limited().energy_tolerance_j));
        }
        close(selected.first_half.dirichlet_in_j + selected.second_half.dirichlet_in_j,
            second[1] / 24.0);
        assert_eq!(selected.first_half.second_stage.temperature[3], 310.0);
        assert_eq!(run.conduction().liquid_mass_fraction()[3], 1.0);

        let mut audited = create(cx, &history, &material, &chart, 340.0, 0.0).into_sdirk2();
        let checked = audited.advance_adaptive_step_with_entropy(cx, 8.0, 8.0, limited(),
            recovery_controls(), phase_controls(), 1.0).unwrap();
        assert_eq!(checked.physical.endpoint_solves, 4);
        assert_eq!(checked.physical.rejected_attempts, 1);
        assert_eq!(audited.conduction().specific_enthalpy_j_kg(),
            run.conduction().specific_enthalpy_j_kg());
        close(checked.entropy.thermal_storage_change_j_k,
            audited.conduction().storage().thermal_entropy_change_j_k(cx,
                &[800.0, 800.0, 800.0, 0.0],
                audited.conduction().specific_enthalpy_j_kg()).unwrap());
        checked.entropy.admit(1.0).unwrap();

        let retained = run.conduction().specific_enthalpy_j_kg().to_vec();
        with_cancelled_cx(|cancelled| {
            assert!(run.advance_adaptive_step(cancelled, 8.0, 4.0, limited(),
                recovery_controls(), phase_controls()).is_err());
        });
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), retained);
        assert_eq!(run.conduction().time_s(), 4.0);
        assert_eq!(run.conduction().accepted_steps(), 2);

        let resumed = run.advance_adaptive_step(cx, 8.0, 4.0, limited(),
            recovery_controls(), phase_controls()).unwrap();
        assert_eq!(resumed.endpoint_solves, 3);
        assert_eq!(resumed.rejected_attempts, 0);
        let middle = sdirk(second[1], 2.0, 340.0);
        let final_h = sdirk(middle[1], 2.0, 340.0)[1];
        close(run.conduction().specific_enthalpy_j_kg()[3], final_h);
        assert_eq!(run.conduction().time_s(), 8.0);
        assert_eq!(run.conduction().accepted_steps(), 4);
    });
}

#[test]
fn g0_recovery_floor_and_attempt_cap_retain_the_original_typed_refusal() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let chart = curve(0.0);
    with_cx(|cx| {
        let history = bath_history(cx, &mesh, 340.0, Vec::new());
        for at_floor in [false, true] {
            let mut run = create(cx, &history, &material, &chart, 340.0, 0.0);
            let original = run.advance_to(cx, 4.0, limited()).unwrap_err();
            let policy = if at_floor {
                AdaptiveHistoryConfig { min_step_s: 4.0, ..recovery_controls() }
            } else {
                AdaptiveHistoryConfig { max_attempts: 1, ..recovery_controls() }
            };
            let error = run.advance_adaptive_step(cx, 4.0, 4.0, limited(),
                policy, phase_controls()).unwrap_err();
            assert_eq!(error, original);
            untouched(&run, 340.0, 0.0);
            let error = run.advance_adaptive_step_with_entropy(cx, 4.0, 4.0, limited(),
                policy, phase_controls(), 1.0).unwrap_err();
            assert_eq!(error, EntropyAdaptiveHistoryError::History(original));
            untouched(&run, 340.0, 0.0);

            let mut sdirk = run.into_sdirk2();
            let original = sdirk.advance_to(cx, 8.0, limited()).unwrap_err();
            let policy = if at_floor {
                AdaptiveHistoryConfig { min_step_s: 8.0, ..recovery_controls() }
            } else { policy };
            let error = sdirk.advance_adaptive_step(cx, 8.0, 8.0, limited(),
                policy, phase_controls()).unwrap_err();
            assert_eq!(error, original);
            untouched(sdirk.conduction(), 340.0, 0.0);
        }
        // Exhaust a real recovery after two different failed intervals.
        // The returned error belongs to the final eight-second trial.
        let mut run = create(cx, &history, &material, &chart, 340.0, 0.0).into_sdirk2();
        let last = run.advance_to(cx, 8.0, limited()).unwrap_err();
        let policy = AdaptiveHistoryConfig { max_step_s: 16.0, max_attempts: 2,
            ..recovery_controls() };
        let error = run.advance_adaptive_step(cx, 16.0, 16.0, limited(),
            policy, phase_controls()).unwrap_err();
        assert_eq!(error, last);
        untouched(run.conduction(), 340.0, 0.0);
    });
}

#[test]
fn g0_unsupported_scheduled_bath_is_not_hidden_by_a_shorter_interval() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let chart = curve(0.0);
    with_cx(|cx| {
        // The requested endpoint is outside the original chart; a bisected
        // endpoint would be inside it. This must refuse, not silently stop early.
        let row = BoundaryCondition {
            region: "bath".into(), physics: Physics::Thermal, kind: BcKind::Dirichlet,
            value: Some(BcValue::Signal(ramp(340.0, 410.0, 4.0, TEMPERATURE_DIMS))),
            compatibility: None, frame: 0,
        };
        let history = bath_history(cx, &mesh, 340.0,
            vec![ThermalRegionHistory { row, reference_temperature: None }]);
        let mut run = create(cx, &history, &material, &chart, 340.0, 0.0);
        let error = run.advance_adaptive_step(cx, 4.0, 4.0, limited(),
            recovery_controls(), phase_controls()).unwrap_err();
        assert!(matches!(error, HeterogeneousEnthalpyError::Enthalpy(
            EnthalpyError::PrescribedTemperatureOutsideDomain { temperature_k: 410.0, .. })));
        untouched(&run, 340.0, 0.0);
        let audited = run.advance_adaptive_step_with_entropy(cx, 4.0, 4.0, limited(),
            recovery_controls(), phase_controls(), 1.0).unwrap_err();
        assert_eq!(audited, EntropyAdaptiveHistoryError::History(error.clone()));
        untouched(&run, 340.0, 0.0);
        let mut sdirk = run.into_sdirk2();
        let sdirk_error = sdirk.advance_adaptive_step(cx, 4.0, 4.0, limited(),
            recovery_controls(), phase_controls()).unwrap_err();
        assert_eq!(sdirk_error, error);
        untouched(sdirk.conduction(), 340.0, 0.0);
        let audited = sdirk.advance_adaptive_step_with_entropy(cx, 4.0, 4.0, limited(),
            recovery_controls(), phase_controls(), 1.0).unwrap_err();
        assert_eq!(audited, EntropyAdaptiveHistoryError::History(error));
        untouched(sdirk.conduction(), 340.0, 0.0);
    });
}
