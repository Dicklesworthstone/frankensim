//! Independent annular momentum and nonlinear equipment-source regressions.

use super::super::head::{FlowHeadLaw, FlowHeadResponse};
use super::*;

#[derive(Debug)]
struct Cubic {
    shutoff: f64,
    linear: f64,
    cubic: f64,
    support: [f64; 2],
}
impl Cubic {
    fn pressure(&self, q: f64) -> f64 {
        self.shutoff - self.linear * q - self.cubic * q * q * q
    }
    fn slope(&self, q: f64) -> f64 {
        -self.linear - 3.0 * self.cubic * q * q
    }
}
impl FlowHeadLaw for Cubic {
    fn flow_interval_m3_s(&self) -> [f64; 2] {
        self.support
    }
    fn evaluate(
        &self,
        flow: VolumetricFlowRate,
    ) -> Result<FlowHeadResponse, Box<dyn std::error::Error + Send + Sync>> {
        let q = flow.value();
        if !q.is_finite() || q < self.support[0] || q > self.support[1] {
            return Err(std::io::Error::other("original equipment flow support").into());
        }
        Ok(FlowHeadResponse {
            pressure: Pressure::new(self.pressure(q)),
            pressure_slope_pa_s_m3: self.slope(q),
        })
    }
}
fn equipment(shutoff: f64) -> Arc<Cubic> {
    Arc::new(Cubic {
        shutoff,
        linear: 2e6,
        cubic: 4e15,
        support: [-0.001, 0.001],
    })
}
fn single() -> TransientCarreauGraph {
    TransientCarreauGraph::new(
        CarreauPipeGraph::new(2, vec![branch("original", 0, 1, 0.5, 0.01, false)]).unwrap(),
        1000.0,
        &[2],
        vec![None; 2],
    )
    .unwrap()
}

// Direct two-annulus midpoint matrix, independent of the native radial solver.
fn radial(old: [f64; 2], dt: f64, pressure: f64) -> ([f64; 2], f64, f64, f64) {
    let area = std::f64::consts::PI * 0.01_f64.powi(2);
    let a = [area / 4.0, 3.0 * area / 4.0];
    let m = [500.0 * a[0], 500.0 * a[1]];
    let k = std::f64::consts::PI * 0.2;
    let aa = m[0] + dt * k / 2.0;
    let bb = -dt * k / 2.0;
    let dd = m[1] + dt * 5.0 * k / 2.0;
    let f0 = (m[0] - dt * k / 2.0) * old[0] + dt * k / 2.0 * old[1] + dt * a[0] * pressure;
    let f1 = dt * k / 2.0 * old[0] + (m[1] - dt * 5.0 * k / 2.0) * old[1] + dt * a[1] * pressure;
    let det = aa * dd - bb * bb;
    let next = [(dd * f0 - bb * f1) / det, (aa * f1 - bb * f0) / det];
    let mid = [0.5 * (old[0] + next[0]), 0.5 * (old[1] + next[1])];
    let q = a[0] * mid[0] + a[1] * mid[1];
    let heat = dt * k * ((mid[0] - mid[1]).powi(2) + 4.0 * mid[1].powi(2));
    let kinetic = (0..2).map(|i| m[i] * mid[i] * (next[i] - old[i])).sum();
    (next, q, heat, kinetic)
}
fn independent(old: [f64; 2], dt: f64, base: f64, law: &Cubic) -> ([f64; 2], f64, f64, f64, f64) {
    let offset = radial(old, dt, 0.0).1;
    let g = radial([0.0; 2], dt, 1.0).1;
    let residual = |q: f64| q - offset - g * (base + law.pressure(q));
    let [mut lo, mut hi] = law.support;
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..90 {
        let middle = f64::midpoint(lo, hi);
        if residual(middle) > 0.0 {
            hi = middle;
        } else {
            lo = middle;
        }
    }
    let q = f64::midpoint(lo, hi);
    let (v, flow, heat, kinetic) = radial(old, dt, base + law.pressure(q));
    (v, flow, heat, kinetic, g / (1.0 - g * law.slope(q)))
}

#[test]
fn g1_nonlinear_source_startup_reversal_and_recovering_work_match_independent_annuli() {
    with_cx(|cx| {
        let model = single();
        for (old, shutoff, base) in [
            ([0.0; 2], 800.0, 100.0),
            ([0.0; 2], -800.0, -100.0),
            ([0.6, 0.3], -200.0, 0.0),
        ] {
            let original = equipment(shutoff);
            let shared = FlowHeadSource::new(original.clone());
            let mut state = model
                .state_from_velocities(vec![old.to_vec()], &[None; 2], 0.0, cx)
                .unwrap();
            let dt = 0.02;
            let expected = independent(old, dt, base, &original);
            let step = model
                .step_with_head_laws(
                    &mut state,
                    &mut model.workspace().unwrap(),
                    &[boundary(0, base), boundary(1, 0.0)],
                    &[],
                    None,
                    &[Some(shared.clone())],
                    dt,
                    controls(),
                    cx,
                )
                .unwrap();
            balance(&step);
            for (&actual, expected) in state.branches()[0].velocities_m_s().iter().zip(expected.0) {
                near(actual, expected, 2e-9, 1e-12);
            }
            let accepted = step.flow_heads[0].as_ref().unwrap();
            assert_eq!(accepted.source, shared);
            near(accepted.midpoint_flow.value(), expected.1, 2e-9, 1e-16);
            near(
                accepted.effective_conductance_m3_s_pa,
                expected.4,
                2e-9,
                1e-18,
            );
            near(step.viscous_heat_j, expected.2, 2e-9, 1e-12);
            near(step.kinetic_energy_change_j, expected.3, 2e-9, 1e-12);
            near(
                step.branch_head_work_j,
                original.pressure(expected.1) * expected.1 * dt,
                2e-9,
                1e-12,
            );
            assert_eq!(step.branch_head_work_j.is_sign_negative(), old[0] > 0.0);
            assert!(step.head_evaluations >= 3 && step.branch_solves > 1);
            assert!(accepted.pressure_residual.value().abs() < 1e-7);
        }
    });
}

#[test]
fn g3_original_nonlinear_material_and_head_tangents_close_the_same_free_node_schur_system() {
    with_cx(|cx| {
        let model = TransientCarreauGraph::new(
            CarreauPipeGraph::new(
                3,
                vec![
                    branch("pump", 0, 1, 0.3, 0.01, true),
                    branch("return", 1, 2, 0.5, 0.012, true),
                ],
            )
            .unwrap(),
            1000.0,
            &[12, 12],
            vec![None, storage(2e-9), None],
        )
        .unwrap();
        let heads = [Some(FlowHeadSource::new(equipment(500.0))), None];
        let state = model
            .rest_state(&[None, Some(Pressure::new(100.0)), None], cx)
            .unwrap();
        let ports = [boundary(0, 200.0), boundary(2, 0.0)];
        let mut next = state.clone();
        let step = model
            .step_with_head_laws(
                &mut next,
                &mut model.workspace().unwrap(),
                &ports,
                &[],
                None,
                &heads,
                0.01,
                controls(),
                cx,
            )
            .unwrap();
        balance(&step);
        assert!(step.newton_iterations > 0);
        let layout = model
            .graph
            .boundary_layout_with_storage(&ports, &[], Some(&state.storage_pressures), cx)
            .unwrap();
        let p: Vec<_> = step.midpoint_pressures.iter().map(|p| p.value()).collect();
        let mut work = model.workspace().unwrap();
        let mut count = 0;
        let mut head_count = 0;
        model
            .evaluate(
                &state,
                &p,
                &layout,
                None,
                Some(&heads),
                0.01,
                controls(),
                &mut work.radial,
                &mut work.outputs,
                &mut work.current,
                &mut count,
                &mut head_count,
                cx,
            )
            .unwrap();
        let expected = work.current.conductance.iter().sum::<f64>() + 2.0 * 2e-9 / 0.01;
        let mut residuals = Vec::new();
        for offset in [0.01, -0.01] {
            let mut perturbed = p.clone();
            perturbed[1] += offset;
            model
                .evaluate(
                    &state,
                    &perturbed,
                    &layout,
                    None,
                    Some(&heads),
                    0.01,
                    controls(),
                    &mut work.radial,
                    &mut work.outputs,
                    &mut work.trial,
                    &mut count,
                    &mut head_count,
                    cx,
                )
                .unwrap();
            residuals.push(work.trial.residuals[1]);
        }
        near((residuals[0] - residuals[1]) / 0.02, expected, 2e-5, 1e-12);
    });
}

#[test]
fn g3_extreme_equipment_endpoints_do_not_require_unsupported_extreme_radial_profiles() {
    with_cx(|cx| {
        let mut branch = branch("bounded-original", 0, 1, 0.5, 0.01, false);
        branch.pipe = LaminarCarreauPipe::new(*branch.pipe.law(), 300.0, 0.01, 0.5, 100.0).unwrap();
        let model = TransientCarreauGraph::new(
            CarreauPipeGraph::new(2, vec![branch]).unwrap(),
            1000.0,
            &[2],
            vec![None; 2],
        )
        .unwrap();
        let old = model.rest_state(&[None; 2], cx).unwrap();
        let mut next = old.clone();
        let law = equipment(500.0);
        // Original curve endpoints imply multi-megapascal trial forcing, whose
        // shear is outside this original pipe's support. The true root is valid.
        assert!(
            model.pipes()[0]
                .step_from(
                    &old.branches()[0],
                    &mut next.branches[0],
                    &mut model.pipes()[0].workspace().unwrap(),
                    law.pressure(law.support[0]),
                    0.02,
                    controls().radial,
                    || false
                )
                .is_err()
        );
        next = old.clone();
        let step = model
            .step_with_head_laws(
                &mut next,
                &mut model.workspace().unwrap(),
                &[boundary(0, 0.0), boundary(1, 0.0)],
                &[],
                None,
                &[Some(FlowHeadSource::new(law))],
                0.02,
                controls(),
                cx,
            )
            .unwrap();
        assert!(step.branches[0].midpoint_volume_flow_m3_s > 0.0);
        balance(&step);
    });
}

#[test]
fn g1_final_equipment_support_is_checked_after_the_coupled_storage_pressure_solve() {
    with_cx(|cx| {
        let model = TransientCarreauGraph::new(
            CarreauPipeGraph::new(2, vec![branch("pump-to-store", 0, 1, 0.5, 0.01, false)])
                .unwrap(),
            1000.0,
            &[2],
            vec![None, storage(2e-11)],
        )
        .unwrap();
        let source = FlowHeadSource::new(Arc::new(Cubic {
            shutoff: 800.0,
            linear: 2e6,
            cubic: 4e15,
            support: [0.0, 6e-6],
        }));
        let laws = [Some(source)];
        let initial = model
            .rest_state(&[None, Some(Pressure::new(0.0))], cx)
            .unwrap();
        let ports = [boundary(0, 0.0)];
        let layout = model
            .graph
            .boundary_layout_with_storage(&ports, &[], Some(&initial.storage_pressures), cx)
            .unwrap();
        let mut workspace = model.workspace().unwrap();
        model
            .evaluate(
                &initial,
                &[0.0; 2],
                &layout,
                None,
                Some(&laws),
                0.02,
                controls(),
                &mut workspace.radial,
                &mut workspace.outputs,
                &mut workspace.current,
                &mut 0,
                &mut 0,
                cx,
            )
            .unwrap();
        let unbalanced = workspace.current.steps[0].unwrap();
        assert!(unbalanced.midpoint_volume_flow_m3_s < 6e-6);
        assert!(unbalanced.final_volume_flow_m3_s > 6e-6);
        let mut accepted = initial.clone();
        let step = model
            .step_with_head_laws(
                &mut accepted,
                &mut workspace,
                &ports,
                &[],
                None,
                &laws,
                0.02,
                controls(),
                cx,
            )
            .unwrap();
        assert!(step.branches[0].final_volume_flow_m3_s < 6e-6);
        assert!(step.midpoint_pressures[1].value() > 0.0);
        assert!(step.newton_iterations > 0);
        near(
            step.storage_volume_changes_m3[1],
            step.branches[0].midpoint_volume_flow_m3_s * 0.02,
            1e-8,
            1e-15,
        );
        balance(&step);
        assert_eq!(initial.accepted_time_s(), 0.0);
    });
}

#[test]
fn g4_equipment_support_slope_and_work_refusals_preserve_native_history_and_retry() {
    with_cx(|cx| {
        let model = single();
        let mut state = model.rest_state(&[None; 2], cx).unwrap();
        let before = snapshot(&state);
        let mut work = model.workspace().unwrap();
        let ports = [boundary(0, 0.0), boundary(1, 0.0)];
        for source in [
            Arc::new(Cubic {
                shutoff: 800.0,
                linear: -1e6,
                cubic: 0.0,
                support: [-0.001, 0.001],
            }),
            Arc::new(Cubic {
                shutoff: 800.0,
                linear: 2e6,
                cubic: 4e15,
                support: [0.0, 1e-6],
            }),
        ] {
            assert!(
                model
                    .step_with_head_laws(
                        &mut state,
                        &mut work,
                        &ports,
                        &[],
                        None,
                        &[Some(FlowHeadSource::new(source))],
                        0.02,
                        controls(),
                        cx
                    )
                    .is_err()
            );
            assert_eq!(snapshot(&state), before);
        }
        let source = FlowHeadSource::new(equipment(800.0));
        let mut limited = controls();
        limited.head.maximum_source_evaluations = 1;
        assert!(matches!(
            model.step_with_head_laws(
                &mut state,
                &mut work,
                &ports,
                &[],
                None,
                &[Some(source.clone())],
                0.02,
                limited,
                cx
            ),
            Err(TransientCarreauError::Head {
                source: FlowHeadError::SourceBudget,
                ..
            })
        ));
        limited = controls();
        limited.maximum_branch_solves = 1;
        assert!(matches!(
            model.step_with_head_laws(
                &mut state,
                &mut work,
                &ports,
                &[],
                None,
                &[Some(source.clone())],
                0.02,
                limited,
                cx
            ),
            Err(TransientCarreauError::BranchSolveBudget)
        ));
        assert_eq!(snapshot(&state), before);
        let accepted = model
            .step_with_head_laws(
                &mut state,
                &mut work,
                &ports,
                &[],
                None,
                &[Some(source)],
                0.02,
                controls(),
                cx,
            )
            .unwrap();
        balance(&accepted);
        assert_eq!(state.accepted_time_s(), 0.02);
    });
}
