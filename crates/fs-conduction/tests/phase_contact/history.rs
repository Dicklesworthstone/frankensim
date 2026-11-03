//! Real phase-dependent contact trajectories on the existing h history owner.
#[path = "history/enclosure.rs"]
mod enclosure;
use super::*;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{
    EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun,
    adaptive::PhaseErrorTolerance, contact::ContactEnthalpyHistoryRun,
};
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyStep;
use fs_conduction::transient::enthalpy::{EnthalpyError, heterogeneous::HeterogeneousEnthalpyError};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};

const OLD: [f64; 8] = [25.0, 25.0, 25.0, 30.0, 60.0, 60.0, 60.0, 65.0];

fn boundary<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh).adiabatic_remainder(),
        &Environment::earth_lab(), Vec::new(), 0.0).unwrap()
}

fn base<'h, 'm, 'c>(cx: &Cx<'_>, boundary: &'h ThermalBoundaryHistory<'m>,
    conductivity: &'h ConductivityModel, curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    old: &[f64], source: HeatSourceHistory, max_steps: usize) -> EnthalpyHistoryRun<'h, 'm, 'c>
{
    EnthalpyHistoryRun::new(cx, ThermalHistoryProblem { boundary, material: conductivity,
        element_materials: None, interfaces: None, source }, EnthalpyHistoryMaterials {
        table: &[ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: 1.0 },
            ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: 2.0 }],
        of_element: &[0, 1], phase_conductivity: None }, old,
        EnthalpyHistoryBudget { spatial: EnthalpyBudget { max_vertices: 8, max_elements: 2 },
            max_steps }).unwrap()
}

fn zero() -> HeatSourceHistory { HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)) }

// R=1+2fa+fb, Q=A(Tb-Ta)/R. Face nodes have mA=1/24, mB=1/12.
// Their implicit source/contact increments give Q*(a+0.12*dt*Q)=50 W.
// This closed form is independent of the production source/tangent/Newton code.
fn endpoint(old: &[f64; 8], dt: f64, qa: f64, qb: f64) -> [f64; 8] {
    let a = 1.0 + 0.02 * (old[0] + dt * qa) + 0.01 * (old[4] + dt * qb / 2.0);
    let q = 100.0 / (a + (a * a + 24.0 * dt).sqrt());
    std::array::from_fn(|i| {
        let body_source = if i < 4 { qa } else { qb / 2.0 };
        let contact = if i % 4 == 3 { 0.0 } else if i < 4 { 8.0 * q } else { -4.0 * q };
        old[i] + dt * (body_source + contact)
    })
}

fn admit(report: &ContactEnthalpyStep, old: &[f64; 8], expected: &[f64; 8], dt: f64, qa: f64, qb: f64) {
    for i in 0..8 {
        near(report.conduction.specific_enthalpy_j_kg[i], expected[i], 2e-8);
        near(report.conduction.liquid_mass_fraction[i], expected[i] / 100.0, 2e-10);
    }
    assert_eq!(report.conduction.temperature, [400.0, 400.0, 400.0, 400.0, 500.0, 500.0, 500.0, 500.0]);
    let stored: f64 = expected.iter().zip(old).enumerate()
        .map(|(i, (new, old))| (new - old) * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 }).sum();
    near(stored, dt * (qa + qb) / 6.0, 1e-9);
    near(report.conduction.stored_energy_change_j, stored, 1e-9);
    near(report.physical.energy_residual_j, 0.0, config().step.energy_tolerance_j);
    assert!(report.physical.residual_norm_j <= config().residual_tolerance_j);
    assert_eq!(report.iterations, 1);
    assert!(!report.physical_contacts.material_trace_states().is_empty());
    for trace in report.physical_contacts.material_trace_states() {
        near(trace.liquid_mass_fractions[0], expected[0] / 100.0, 2e-10);
        near(trace.liquid_mass_fractions[1], expected[4] / 100.0, 2e-10);
    }
}

#[test]
fn g1_latent_contact_histories_follow_both_original_phase_coordinates_and_close_energy() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = boundary(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0); let point = point();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        for mode in 0..3 {
            let contacts = contacts(cx, &mesh, &initial, &cards, &curves, &source, &point,
                sides(&mesh, 4, 2, 0.0), mode, false);
            let mut trajectory = ContactEnthalpyHistoryRun::new_newton(cx,
                base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 4), contacts, config()).unwrap();
            assert!(trajectory.coupling().is_none());
            assert_eq!(trajectory.newton_config().unwrap().step, config().step);
            let mut old = OLD;
            for i in 1..=4 {
                let expected = endpoint(&old, 0.025, 0.0, 0.0);
                let report = trajectory.advance_to(cx, f64::from(i) * 0.025, config().step).unwrap();
                admit(&report, &old, &expected, 0.025, 0.0, 0.0);
                // Heat transfer changes phase while both absolute temperatures stay fixed.
                assert!(expected[0] > old[0] && expected[4] < old[4]);
                old = expected;
            }
            assert_eq!(trajectory.conduction().accepted_steps(), 4);
            assert_eq!(trajectory.conduction().time_s(), 0.1);
            if mode == 0 {
                trajectory.extend_step_budget(54).unwrap();
                for i in 1..=50 {
                    trajectory.advance_to(cx, 0.1 + f64::from(i), config().step).unwrap();
                }
                // Equilibrium lies above A's liquidus and below B's solidus:
                // hA=T-300, hB=T-500, with total masses 1/6 and 1/3 kg.
                let initial_energy: f64 = OLD.iter().enumerate()
                    .map(|(i, h)| h * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 }).sum();
                let equilibrium = (initial_energy + 300.0 / 6.0 + 500.0 / 3.0) / 0.5;
                for &temperature in trajectory.conduction().temperature() {
                    near(temperature, equilibrium, 1e-5);
                }
                assert_eq!(&trajectory.conduction().liquid_mass_fraction()[..4], &[1.0; 4]);
                assert_eq!(&trajectory.conduction().liquid_mass_fraction()[4..], &[0.0; 4]);
                let final_energy: f64 = trajectory.conduction().specific_enthalpy_j_kg().iter().enumerate()
                    .map(|(i, h)| h * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 }).sum();
                near(final_energy, initial_energy, 1e-7);
            }
        }
    });
}

#[test]
fn g1_phase_contact_histories_resolve_pulses_and_reject_temperature_invisible_time_error() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = boundary(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0); let point = point();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let create = || {
            let forcing = HeatSourceHistory::Scaled {
                field: ScalarField::Nodal(vec![60.0, 60.0, 60.0, 60.0, -60.0, -60.0, -60.0, -60.0]),
                scale: TimeSignal::Table { times: vec![0.0, 0.05], values: vec![1.0, -1.0],
                    dims: Dims::NONE, interp: Interp::Hold },
            };
            ContactEnthalpyHistoryRun::new_newton(cx,
                base(cx, &boundary, &conductivity, &curves, &OLD, forcing, 4),
                contacts(cx, &mesh, &initial, &cards, &curves, &source, &point,
                    sides(&mesh, 4, 2, 0.0), 1, false), config()).unwrap()
        };
        let mut events = create();
        let first = events.advance_event_step(cx, 0.15, 1.0, config().step).unwrap();
        assert_eq!(first.end_time_s, 0.05);
        let h1 = endpoint(&OLD, 0.05, 60.0, -60.0);
        admit(&first.report, &OLD, &h1, 0.05, 60.0, -60.0);
        let second = events.advance_event_step(cx, 0.15, 1.0, config().step).unwrap();
        let h2 = endpoint(&h1, 0.1, -60.0, 60.0);
        admit(&second.report, &h1, &h2, 0.1, -60.0, 60.0);
        assert_eq!(second.end_time_s, 0.15);

        let mut adaptive = create();
        let accepted = adaptive.advance_adaptive_step(cx, 0.05, 0.05, config().step,
            AdaptiveHistoryConfig { absolute_temperature_k: 1e-9, relative_temperature: 0.0,
                min_step_s: 1e-6, max_step_s: 0.05, max_attempts: 16 },
            PhaseErrorTolerance { absolute_enthalpy_j_kg: 1.0, relative_enthalpy: 0.0,
                absolute_liquid_fraction: 5e-6 }).unwrap();
        assert!(accepted.rejected_attempts > 0);
        assert!(accepted.end_time_s < 0.05);
        assert!(accepted.normalized_error <= 1.0);
        let dt = accepted.midpoint_time_s;
        let hm = endpoint(&OLD, dt, 60.0, -60.0);
        let hf = endpoint(&hm, accepted.end_time_s - dt, 60.0, -60.0);
        admit(&accepted.first_half, &OLD, &hm, dt, 60.0, -60.0);
        admit(&accepted.second_half, &hm, &hf, accepted.end_time_s - dt, 60.0, -60.0);
        assert_eq!(adaptive.conduction().accepted_steps(), 2);
        assert_eq!(adaptive.conduction().specific_enthalpy_j_kg(), accepted.second_half.conduction.specific_enthalpy_j_kg);
    });
}

#[test]
fn g0_phase_contact_source_cancellation_and_policy_refusals_preserve_resume() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = boundary(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let clones = curves.clone();
        let source = source(&cards, 1.0, false, false, 0.35); let point = point();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let contacts = || super::contacts(cx, &mesh, &initial, &cards, &curves, &source, &point,
            sides(&mesh, 4, 2, 0.0), 1, false);
        assert!(ContactEnthalpyHistoryRun::new_newton(cx,
            base(cx, &boundary, &conductivity, &clones, &OLD, zero(), 2), contacts(), config()).is_err());
        let create = || ContactEnthalpyHistoryRun::new_newton(cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2), contacts(), config()).unwrap();
        assert!(create().with_radiation(cx, Vec::new(), AmbientRadiationConfig::default()).is_err());
        let mut trajectory = create();
        let mut changed = config().step; changed.energy_tolerance_j *= 2.0;
        assert!(trajectory.advance_to(cx, 0.025, changed).is_err());
        support::with_cancelled_cx(|cancelled| assert!(matches!(
            trajectory.advance_to(cancelled, 0.025, config().step),
            Err(HeterogeneousEnthalpyError::Enthalpy(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))))));
        assert!(trajectory.advance_to(cx, 0.2, config().step).is_err());
        assert_eq!(trajectory.conduction().specific_enthalpy_j_kg(), &OLD);
        assert_eq!(trajectory.conduction().time_s(), 0.0);
        assert_eq!(trajectory.conduction().accepted_steps(), 0);
        let accepted = trajectory.advance_to(cx, 0.025, config().step).unwrap();
        admit(&accepted, &OLD, &endpoint(&OLD, 0.025, 0.0, 0.0), 0.025, 0.0, 0.0);
        let replay = create().advance_to(cx, 0.025, config().step).unwrap();
        assert_eq!(accepted.conduction.specific_enthalpy_j_kg, replay.conduction.specific_enthalpy_j_kg);
        assert_eq!(accepted.physical_contacts.material_node_states(), replay.physical_contacts.material_node_states());
    });
}
