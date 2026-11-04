//! Original source-bound solid/time/Newton motion in prescribed hydrostatics.
use super::*;
use fs_solid::hyper3::HyperTetError;
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use std::sync::Arc;

fn surface(cx: &Cx<'_>) -> Arc<ClosedTetPressureSurface> {
    Arc::new(
        ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap(),
    )
}
fn field(clipped: bool) -> TetPressure<'static> {
    if clipped {
        TetPressure::FreeSurface {
            origin_m: [0.0; 3],
            pressure_pa: 0.35,
            gradient_pa_m: [0.0, 0.0, -0.7],
        }
    } else {
        TetPressure::Affine {
            origin_m: [0.0; 3],
            pressure_pa: 0.35,
            gradient_pa_m: [0.0, 0.0, -0.7],
        }
    }
}
fn nodal_pressure(body: &FiniteSolidBody<'_>) -> [f64; 4] {
    let u = body.displacement_m();
    core::array::from_fn(|a| 0.35 - 0.7 * (NODES[a][2] + u[a][2]))
}
// Independent determinant and truncated-power volume integral; fully wet/dry
// cases use their exact affine limits. The rotating fixture separates knots.
fn potential(body: &FiniteSolidBody<'_>, clipped: bool) -> f64 {
    let u = body.displacement_m();
    let f: [f64; 9] =
        core::array::from_fn(|ij| NODES[ij % 3 + 1][ij / 3] + u[ij % 3 + 1][ij / 3] - u[0][ij / 3]);
    let volume = (f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6]))
        / 6.0;
    let p = nodal_pressure(body);
    let mean = if !clipped || p.iter().all(|p| *p >= 0.0) {
        p.iter().sum::<f64>() / 4.0
    } else if p.iter().all(|p| *p <= 0.0) {
        0.0
    } else {
        (0..4)
            .map(|i| {
                p[i].max(0.0).powi(4)
                    / (4.0
                        * (0..4)
                            .filter(|&j| i != j)
                            .map(|j| p[i] - p[j])
                            .product::<f64>())
            })
            .sum()
    };
    volume * mean
}
fn submerged_motion() -> ([[f64; 3]; 4], [[f64; 3]; 4]) {
    let (mut u, mut v) = initial_motion();
    for a in 0..4 {
        u[a][2] += 0.55;
        // A uniform downward speed leaves the tumbling strain motion intact
        // while carrying initially dry vertices through the prescribed plane.
        v[a][2] -= 1.3;
    }
    (u, v)
}

#[test]
fn native_spatial_pressure_motion_closes_work_and_preserves_only_admissible_momenta() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let (u, v) = submerged_motion();
        for clipped in [false, true] {
            let loading = field(clipped);
            let mut body = FiniteSolidBody::new_with_pressure(
                problem(&material, &[], &[]),
                surface(cx),
                loading,
                &u,
                &v,
                0.0,
                Default::default(),
                cx,
            )
            .unwrap();
            let initial = invariants(&body);
            let initial_energy = initial.0 + potential(&body, clipped);
            let mut crossings = 0;
            let accepted = body.native_state().clone();
            for _ in 0..30 {
                let old_pressure = nodal_pressure(&body);
                let old_potential = potential(&body, clipped);
                let old_state = body.native_state().clone();
                let (next, report) = body
                    .step_energy_momentum_candidate_with_pressure(
                        0.02,
                        &[],
                        [0.0; 3],
                        loading,
                        controls(1e-10),
                        cx,
                    )
                    .unwrap();
                let measured = invariants(&next);
                let next_potential = potential(&next, clipped);
                near(measured.0 + next_potential, initial_energy, 3e-9);
                near(
                    report.pressure_work_j,
                    old_potential - next_potential,
                    1e-11,
                );
                near(report.source_work_j, report.pressure_work_j, 1e-14);
                near(report.integration_defect_j, 0.0, 1e-10);
                near(measured.1[0], initial.1[0], 1e-10);
                near(measured.1[1], initial.1[1], 1e-10);
                near(measured.2[2], initial.2[2], 1e-10);
                assert!(report.support_reactions_n.is_empty());
                assert_eq!(body.native_state(), &old_state);
                crossings += old_pressure
                    .into_iter()
                    .zip(nodal_pressure(&next))
                    .filter(|(old, new)| (*old > 0.0) != (*new > 0.0))
                    .count();
                body = next;
            }
            assert_eq!(accepted.steps, 0);
            assert_eq!(body.native_state().steps, 30);
            assert!(
                (invariants(&body).1[2] - initial.1[2]).abs() > 1e-4,
                "hydrostatic force must change vertical momentum"
            );
            if clipped {
                assert!(
                    crossings >= 2,
                    "fixture must cross the free plane during actual Newton motion; crossings={crossings}"
                );
            }
        }
    });
}

#[test]
fn spatial_pressure_refusal_and_cancellation_leave_the_native_candidate_unchanged() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let (u, v) = submerged_motion();
        let body = FiniteSolidBody::new_with_pressure(
            problem(&material, &[], &[]),
            surface(cx),
            field(true),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let accepted = body.native_state().clone();
        for invalid in [
            TetPressure::PerFace(&[0.0; 4]),
            TetPressure::FreeSurface {
                origin_m: [0.0; 3],
                pressure_pa: 0.35,
                gradient_pa_m: [0.0, 0.0, f64::NAN],
            },
        ] {
            assert!(
                body.step_energy_momentum_candidate_with_pressure(
                    0.02,
                    &[],
                    [0.0; 3],
                    invalid,
                    controls(1e-10),
                    cx
                )
                .is_err()
            );
            assert_eq!(body.native_state(), &accepted);
        }
        gate.request();
        assert!(matches!(
            body.step_energy_momentum_candidate_with_pressure(
                0.02,
                &[],
                [0.0; 3],
                field(true),
                controls(1e-10),
                cx
            ),
            Err(FiniteSolidError::Solid(HyperTetError::Cancelled))
        ));
        assert_eq!(body.native_state(), &accepted);
    });
}
