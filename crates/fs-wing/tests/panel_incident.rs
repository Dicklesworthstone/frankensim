//! Actual control-point boundary loading and bound-point force (G0/G1/G3).
use fs_wing::nonlinear::{
    InfluenceOperator, SectionFrame, StripRegime, StripSpec, solve_nonlinear_with_panel_incident,
    solve_nonlinear_with_section_frames,
};
use fs_wing::{
    Panel, PanelIncidentFlow, SurfaceId, flat_surface, induced_velocity_free,
    solve_weissinger_linear, solve_weissinger_linear_with_incident,
};

const RHO: f64 = 1.2;
const WAKE: [f64; 3] = [10.0, 0.0, 0.5];
fn layout() -> (Vec<Panel>, Vec<StripSpec>) {
    let panels = flat_surface(SurfaceId::WingLower, 4.0, 1.0, 1.0, 0.0, 2, 3).unwrap();
    let strips = (0..2)
        .map(|s| StripSpec {
            panel_indices: vec![s, 2 + s, 4 + s],
            chord_m: 1.0,
            twist_rad: 0.0,
        })
        .collect();
    (panels, strips)
}
fn polar(_: usize, alpha: f64) -> (f64, StripRegime) {
    (core::f64::consts::TAU * alpha, StripRegime::Attached)
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn repeated_uniform_samples_keep_original_linear_and_nonlinear_arithmetic() {
    let (panels, strips) = layout();
    let samples = vec![WAKE; panels.len()];
    let flow = PanelIncidentFlow {
        control_mps: &samples,
        bound_mps: &samples,
    };
    let original = InfluenceOperator::build(&panels, WAKE, RHO).unwrap();
    let resolved = InfluenceOperator::build_with_panel_incident(&panels, WAKE, RHO, flow).unwrap();
    assert_eq!(original.linear(), resolved.linear());
    assert_eq!(
        solve_weissinger_linear(&panels, WAKE, RHO).unwrap(),
        solve_weissinger_linear_with_incident(&panels, WAKE, RHO, flow).unwrap()
    );
    let frames = [SectionFrame::BODY_XZ; 2];
    let old = solve_nonlinear_with_section_frames(
        &original, &panels, &strips, WAKE, RHO, &polar, &frames, &[WAKE; 2], None,
    )
    .unwrap();
    let actual = solve_nonlinear_with_panel_incident(
        &resolved, &panels, &strips, WAKE, RHO, &polar, &frames, &[WAKE; 2], flow, None,
    )
    .unwrap();
    assert_eq!(actual, old);
}

#[test]
fn manufactured_control_loading_drives_row_ratios_and_independent_bound_forces() {
    let (panels, strips) = layout();
    let specified = [-0.6, -0.8, -0.3, -0.2, -0.1, -0.15];
    // Manufacture only the normal component of the incident control velocity
    // from a specified physical vortex system. The native linear solve must
    // recover that circulation, irrespective of the separate force samples.
    let control: Vec<_> = panels
        .iter()
        .map(|p| {
            let w = induced_velocity_free(p.ctrl, &panels, &specified, WAKE);
            [10.0, 0.1, -w[2]]
        })
        .collect();
    let bound: Vec<_> = panels
        .iter()
        .enumerate()
        .map(|(j, p)| [9.0 + 0.2 * j as f64, 0.2 * p.a[0], 0.3 + 0.1 * p.a[0]])
        .collect();
    let flow = PanelIncidentFlow {
        control_mps: &control,
        bound_mps: &bound,
    };
    let op = InfluenceOperator::build_with_panel_incident(&panels, WAKE, RHO, flow).unwrap();
    for (&actual, expected) in op.linear().gamma.iter().zip(specified) {
        close(actual, expected, 2e-12);
    }
    for (j, p) in panels.iter().enumerate() {
        let induced = induced_velocity_free(p.ctrl, &panels, &op.linear().gamma, WAKE);
        close(control[j][2] + induced[2], 0.0, 2e-12);
    }
    let section = [bound[0], bound[1]];
    // A constant finite section law isolates the row distribution and force
    // integration from the independently tested nonlinear induction law.
    let closure = |s, _: f64| (0.3 + 0.05 * s as f64, StripRegime::Attached);
    let result = solve_nonlinear_with_panel_incident(
        &op,
        &panels,
        &strips,
        WAKE,
        RHO,
        &closure,
        &[SectionFrame::BODY_XZ; 2],
        &section,
        flow,
        None,
    )
    .unwrap();
    for (s, strip) in strips.iter().enumerate() {
        let target = -0.5 * section[s][0].hypot(section[s][2]) * (0.3 + 0.05 * s as f64);
        let original_total: f64 = strip.panel_indices.iter().map(|&j| specified[j]).sum();
        for &j in &strip.panel_indices {
            close(
                result.gamma[j],
                target * specified[j] / original_total,
                3e-10,
            );
        }
    }
    for (gamma, actual) in [
        (&op.linear().gamma, op.linear().total_force_n),
        (&result.gamma, result.total_force_n),
    ] {
        let mut force = [0.0; 3];
        for (j, p) in panels.iter().enumerate() {
            let segment = core::array::from_fn(|k| p.b[k] - p.a[k]);
            let row = cross(bound[j], segment);
            for k in 0..3 {
                force[k] += RHO * gamma[j] * row[k];
            }
            close(
                row.iter().zip(bound[j]).map(|(f, u)| f * u).sum(),
                0.0,
                2e-13,
            );
        }
        for k in 0..3 {
            close(actual[k], force[k], 2e-11);
        }
    }
}

#[test]
fn changed_trace_reverse_row_and_unresolved_opposing_rows_refuse() {
    let (panels, strips) = layout();
    let samples = vec![WAKE; panels.len()];
    let flow = PanelIncidentFlow {
        control_mps: &samples,
        bound_mps: &samples,
    };
    let op = InfluenceOperator::build_with_panel_incident(&panels, WAKE, RHO, flow).unwrap();
    for control_changes in [true, false] {
        let mut changed = samples.clone();
        changed[4][2] += 0.1;
        let changed = PanelIncidentFlow {
            control_mps: if control_changes { &changed } else { &samples },
            bound_mps: if control_changes { &samples } else { &changed },
        };
        assert_eq!(
            solve_nonlinear_with_panel_incident(
                &op,
                &panels,
                &strips,
                WAKE,
                RHO,
                &polar,
                &[SectionFrame::BODY_XZ; 2],
                &[WAKE; 2],
                changed,
                None,
            )
            .unwrap_err()
            .code,
            "influence-operator-stale"
        );
    }
    assert_eq!(
        solve_nonlinear_with_section_frames(
            &op,
            &panels,
            &strips,
            WAKE,
            RHO,
            &polar,
            &[SectionFrame::BODY_XZ; 2],
            &[WAKE; 2],
            None,
        )
        .unwrap_err()
        .code,
        "influence-operator-stale"
    );
    for invalid in [f64::NAN, -1.0] {
        let mut control = samples.clone();
        control[4][0] = invalid;
        let flow = PanelIncidentFlow {
            control_mps: &control,
            bound_mps: &samples,
        };
        if invalid.is_nan() {
            assert_eq!(
                InfluenceOperator::build_with_panel_incident(&panels, WAKE, RHO, flow)
                    .unwrap_err()
                    .code,
                "incident-flow-invalid"
            );
        } else {
            let op =
                InfluenceOperator::build_with_panel_incident(&panels, WAKE, RHO, flow).unwrap();
            assert_eq!(
                solve_nonlinear_with_panel_incident(
                    &op,
                    &panels,
                    &strips,
                    WAKE,
                    RHO,
                    &polar,
                    &[SectionFrame::BODY_XZ; 2],
                    &[WAKE; 2],
                    flow,
                    None,
                )
                .unwrap_err()
                .code,
                "incident-flow-invalid"
            );
        }
    }
    let short = PanelIncidentFlow {
        control_mps: &samples[..1],
        bound_mps: &samples,
    };
    assert_eq!(
        InfluenceOperator::build_with_panel_incident(&panels, WAKE, RHO, short)
            .unwrap_err()
            .code,
        "incident-flow-invalid"
    );
    let specified = [1.0, 1.0, -0.5, -0.5, -0.5, -0.5];
    let control: Vec<_> = panels
        .iter()
        .map(|p| {
            let w = induced_velocity_free(p.ctrl, &panels, &specified, WAKE);
            [10.0, 0.0, -w[2]]
        })
        .collect();
    let flow = PanelIncidentFlow {
        control_mps: &control,
        bound_mps: &samples,
    };
    let op = InfluenceOperator::build_with_panel_incident(&panels, WAKE, RHO, flow).unwrap();
    assert_eq!(
        solve_nonlinear_with_panel_incident(
            &op,
            &panels,
            &strips,
            WAKE,
            RHO,
            &polar,
            &[SectionFrame::BODY_XZ; 2],
            &[WAKE; 2],
            flow,
            None,
        )
        .unwrap_err()
        .code,
        "panel-row-loading-degenerate"
    );
}
