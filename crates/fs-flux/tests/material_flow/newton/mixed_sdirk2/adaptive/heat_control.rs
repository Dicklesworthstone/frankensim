//! G1/G4: local physical heat and signed SDIRK2 time-defect control.

use super::*;

fn heat_controls(tolerance: f64, defect_budget: f64) -> MaterialSdirk2HeatControls {
    let local = MaterialSdirk2HeatTolerance {
        absolute_j_per_m: tolerance,
        relative: 0.0,
    };
    MaterialSdirk2HeatControls {
        cell_bulk: local,
        boundary_slip: local,
        total_bulk: local,
        total_slip: local,
        absolute_temporal_defect_budget_j_per_m: defect_budget,
        relative_temporal_defect_budget: 0.0,
    }
}

#[test]
fn g1_heat_acceptance_refines_affine_shear_even_when_velocity_is_exact() {
    let gate = CancelGate::new_clock_free();
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 3.0;
    let mu = 0.6;
    let beta = 1.2;
    let reference = transport(density, mu);
    let owner = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            penalty: 15.0,
            linear_tolerance: 1e-12,
            ..MaterialFluxNumerics::default()
        },
    )
    .unwrap();
    let velocity = |time, p: [f64; 2]| [(1.0 + time) * p[1], 0.0];
    let wall = |time, p: [f64; 2]| {
        let sign = if p[1] < 0.5 { -1.0 } else { 1.0 };
        [(1.0 + time) * (p[1] + sign * mu / beta), 0.0]
    };
    let traction =
        |time, _: [f64; 2], n: [f64; 2]| [mu * (1.0 + time) * n[1], mu * (1.0 + time) * n[0]];
    let resistance = |_: f64| beta;
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[1].abs() > 0.5 {
                MaterialNavierStokesTimeBoundary::NavierSlip {
                    edge,
                    value: &wall,
                    tangential_resistance_pa_s_per_m: &resistance,
                }
            } else if geometry.normal[0] < -0.5 {
                MaterialNavierStokesTimeBoundary::OpenReservoir {
                    edge,
                    traction: &traction,
                    velocity: &velocity,
                }
            } else {
                MaterialNavierStokesTimeBoundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            }
        })
        .collect();
    let previous = moments(&mesh, &|p| velocity(0.0, p));
    let mut control = controls(1e-5);
    control.initial_step_s = 1.0;
    control.maximum_step_s = 1.0;
    let source = |_: f64, _: MaterialFluxRheologyPoint, _: &FluidTransportState| {
        Ok(MaterialFluxViscosityTangent {
            material: reference,
            derivative_pa_s2: 0.0,
        })
    };
    let force = |_: f64, p: [f64; 2]| [density * p[1], 0.0];
    let old = owner
        .advance_navier_stokes_mixed_rheology_sdirk2(
            0.0,
            1.0,
            &previous,
            &force,
            &boundaries,
            &control,
            source,
            &gate,
        )
        .unwrap();
    assert_eq!(old.report.accepted_comparisons, 1);
    let result = owner
        .advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control(
            0.0,
            1.0,
            &previous,
            &force,
            &boundaries,
            &control,
            &heat_controls(1e-6, 1e-4),
            source,
            &gate,
        )
        .unwrap();
    assert!(result.flow.report.accepted_comparisons > 4);
    assert!(result.flow.report.rejected_comparisons > 0);
    assert!(result.heat_control.maximum_local_heat_error_ratio <= 1.0);
    assert!(result.heat_control.maximum_temporal_defect_ratio < 0.01);
    let exact = 1.4;
    for (coarse, refined) in [
        (
            old.dissipated_heat.bulk_energy_j_per_m(),
            result.flow.dissipated_heat.bulk_energy_j_per_m(),
        ),
        (
            old.dissipated_heat.slip_energy_j_per_m(),
            result.flow.dissipated_heat.slip_energy_j_per_m(),
        ),
    ] {
        assert!((refined - exact).abs() < 0.05 * (coarse - exact).abs());
        near(refined, exact, 1.5e-4);
    }
    assert_eq!(result.flow.dissipated_heat.duration_s(), 1.0);
    for half in [&result.flow.first_half, &result.flow.second_half] {
        for stage in [&half.first_stage, &half.endpoint] {
            for sample in stage.flow.samples() {
                assert_eq!(sample.material(), &reference);
            }
        }
    }
}

// Curl of the continuous piecewise-P2 interpolant of sin(pi*x)*sin(pi*y)/pi.
// Edge moments are exact for that interpolant: this initial BDM1 field is
// divergence-free and has zero boundary-normal flux on every triangle.
fn vortex_moments(mesh: &TriMesh) -> Vec<f64> {
    let psi = |p: [f64; 2]| {
        if p[0] == 0.0 || p[0] == 1.0 || p[1] == 0.0 || p[1] == 1.0 {
            0.0
        } else {
            fs_math::det::sin(core::f64::consts::PI * p[0])
                * fs_math::det::sin(core::f64::consts::PI * p[1])
                / core::f64::consts::PI
        }
    };
    mesh.edges
        .iter()
        .flat_map(|edge| {
            let a = mesh.verts[edge.verts.0];
            let b = mesh.verts[edge.verts.1];
            let dx = b[0] - a[0];
            let dy = b[1] - a[1];
            let length = dx.hypot(dy);
            let sign = (-edge.normal[1] * dx + edge.normal[0] * dy) / length;
            let middle = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
            [
                sign * (psi(b) - psi(a)) / length,
                sign * (psi(a) - 2.0 * psi(middle) + psi(b)) / (3.0 * length),
            ]
        })
        .collect()
}

#[test]
fn g1_signed_time_defect_prevents_stiff_shear_heat_disappearing_at_small_endpoints() {
    let gate = CancelGate::new_clock_free();
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let density = 2.0;
    let reference = transport(density, density * 100.0);
    let owner = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            penalty: 15.0,
            linear_tolerance: 1e-12,
            ..MaterialFluxNumerics::default()
        },
    )
    .unwrap();
    let previous = vortex_moments(&mesh);
    let saved = previous.clone();
    let zero = |_: f64, _: [f64; 2]| [0.0; 2];
    let perfect_slip = |_: f64| 0.0;
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|edge| mesh.is_boundary(*edge))
        .map(|edge| MaterialNavierStokesTimeBoundary::NavierSlip {
            edge,
            value: &zero,
            tangential_resistance_pa_s_per_m: &perfect_slip,
        })
        .collect();
    let source = |_: f64, _: MaterialFluxRheologyPoint, _: &FluidTransportState| {
        Ok(MaterialFluxViscosityTangent {
            material: reference,
            derivative_pa_s2: 0.0,
        })
    };
    let mut control = controls(0.1);
    control.initial_step_s = 0.1;
    control.maximum_step_s = 0.1;
    control.maximum_attempts = 2000;
    let coarse = owner
        .advance_navier_stokes_mixed_rheology_sdirk2(
            0.0,
            0.1,
            &previous,
            &zero,
            &boundaries,
            &control,
            source,
            &gate,
        )
        .unwrap();
    assert_eq!(coarse.report.accepted_comparisons, 1);
    let result = owner
        .advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control(
            0.0,
            0.1,
            &previous,
            &zero,
            &boundaries,
            &control,
            &heat_controls(1.0, 0.01),
            source,
            &gate,
        )
        .unwrap();
    assert!(result.flow.report.rejected_comparisons > 0);
    assert!(result.flow.report.accepted_comparisons > 1);
    assert!(
        result.flow.dissipated_heat.bulk_energy_j_per_m()
            > 2.0 * coarse.dissipated_heat.bulk_energy_j_per_m()
    );
    assert_eq!(result.flow.dissipated_heat.slip_energy_j_per_m(), 0.0);
    let report = result.heat_control;
    assert!(report.maximum_temporal_defect_ratio <= 1.0);
    assert!(report.absolute_temporal_kinetic_defect_j_per_m <= 0.01);
    assert!(
        report.absolute_temporal_kinetic_defect_j_per_m
            >= report.signed_temporal_kinetic_defect_j_per_m.abs()
    );
    assert!(result.flow.report.stage_solves <= 6 * result.flow.report.attempts);
    assert!(result.flow.report.stage_solves >= 6 * result.flow.report.accepted_comparisons);
    assert_eq!(previous, saved);
    // Spatial SIP/upwind losses are separate. This test does not equate physical
    // bulk heat with the initial continuum kinetic energy or assert a continuum
    // eigenvalue on this coarse mesh.
}

#[test]
fn g4_heat_controls_and_late_work_refusal_leave_history_and_sources_unpublished() {
    let gate = CancelGate::new_clock_free();
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 0.4);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let previous = moments(&mesh, &|_| [1.0, 0.0]);
    let saved = previous.clone();
    let traction = |_: f64, _: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|edge| mesh.is_boundary(*edge))
        .map(|edge| MaterialNavierStokesTimeBoundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    let force = |t: f64, _: [f64; 2]| [6.0 * t * t, 0.0];
    let source = |_: f64, _: MaterialFluxRheologyPoint, _: &FluidTransportState| {
        Ok(MaterialFluxViscosityTangent {
            material: reference,
            derivative_pa_s2: 0.0,
        })
    };
    let mut invalid = heat_controls(1.0, 1.0);
    invalid.relative_temporal_defect_budget = f64::NAN;
    let bad = owner.advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control(
        0.0,
        0.4,
        &previous,
        &|_, _| panic!("invalid controls queried force"),
        &boundaries,
        &controls(1e-5),
        &invalid,
        |_, _, _| -> Result<MaterialFluxViscosityTangent<FluidTransportState>, MaterialError> {
            panic!("invalid controls queried material")
        },
        &gate,
    );
    assert!(matches!(bad, Err(MaterialSdirk2AdaptiveError::Invalid(_))));
    let mut limited = controls(1.0);
    limited.initial_step_s = 0.1;
    limited.maximum_step_s = 0.1;
    limited.maximum_attempts = 1;
    assert!(matches!(
        owner.advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control(
            0.0,
            0.4,
            &previous,
            &force,
            &boundaries,
            &limited,
            &heat_controls(1.0, 1.0),
            source,
            &gate,
        ),
        Err(MaterialSdirk2AdaptiveError::WorkLimit { attempts: 1, .. })
    ));
    let mut continuation = limited.clone();
    continuation.maximum_attempts = 10;
    let refusal = MaterialError::State {
        what: "late original heat-control source refusal".into(),
    };
    let refused = owner.advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control(
        0.0,
        0.4,
        &previous,
        &force,
        &boundaries,
        &continuation,
        &heat_controls(1.0, 1.0),
        |time, _, _| {
            if time > 0.25 {
                return Err(refusal.clone());
            }
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
        &gate,
    );
    assert!(
        matches!(refused, Err(MaterialSdirk2AdaptiveError::Stage(failure))
        if matches!(&failure.error.error, MaterialFluxError::Material(error) if error == &refusal))
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        owner.advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control(
            0.0,
            0.4,
            &previous,
            &force,
            &boundaries,
            &controls(1e-5),
            &heat_controls(1.0, 1.0),
            source,
            &cancelled,
        ),
        Err(MaterialSdirk2AdaptiveError::Cancelled)
    ));
    let retry = owner
        .advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control(
            0.0,
            0.4,
            &previous,
            &force,
            &boundaries,
            &controls(1e-5),
            &heat_controls(1.0, 1.0),
            source,
            &gate,
        )
        .unwrap();
    assert_eq!(retry.flow.report.end_time_s, 0.4);
    assert_eq!(previous, saved);
    // Accelerated rigid translation has essentially zero viscous heat but a
    // signed time-quadrature defect. The diagnostic is not labelled viscosity.
    assert!(retry.heat_control.signed_temporal_kinetic_defect_j_per_m > 0.0);
    assert!(retry.flow.dissipated_heat.bulk_energy_j_per_m() < 1e-12);
}
