//! Scalar equations, thermal budgets and original modeled-diffusivity coupling.

use super::*;
use fs_flux::ns::transport::scalar::{
    ScalarBoundary, ScalarSolveControls, ScalarTransportControls, ScalarTransportError,
};

fn values(mesh: &TriMesh, f: &dyn Fn([f64; 2]) -> f64) -> Vec<f64> {
    mesh.tris
        .iter()
        .flat_map(|tri| tri.map(|i| f(mesh.verts[i])))
        .collect()
}
fn prescribed<'a>(mesh: &TriMesh, value: &'a dyn Fn([f64; 2]) -> f64) -> Vec<ScalarBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(i, _)| mesh.is_boundary(*i))
        .map(|(edge, _)| ScalarBoundary::Value { edge, value })
        .collect()
}
fn insulated<'a>(
    mesh: &TriMesh,
    zero: &'a dyn Fn([f64; 2]) -> f64,
    inflow: Option<&'a dyn Fn([f64; 2]) -> f64>,
) -> Vec<ScalarBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(i, _)| mesh.is_boundary(*i))
        .map(|(edge, _)| ScalarBoundary::Flux {
            edge,
            outward: zero,
            inflow,
        })
        .collect()
}
fn controls() -> ScalarSolveControls {
    ScalarSolveControls {
        relative_tolerance: 2e-12,
        maximum_iterations: 1000,
    }
}

#[test]
fn all_sgs_models_evolve_temperature_with_their_original_endpoint_diffusivity_and_bdf_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let flow = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    let velocity = moments(&mesh, &flow);
    let nc = mesh.tris.len();
    let width = vec![0.4; nc];
    let viscosity = vec![0.01; nc];
    let molecular = vec![0.02; nc];
    let ratio = vec![0.7; nc];
    for model in [
        SgsModel::Smagorinsky { cs: 0.4 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.1 },
    ] {
        let closure = PlanarSgs::with_model(&system, model, &width, &viscosity).unwrap();
        let current_t: f64 = 0.4;
        let dt = 0.07;
        let old_dt = 0.1;
        let target = |p: [f64; 2]| 300.0 + p[0] + 2.0 * p[1] + (current_t + dt).powi(2);
        let previous = values(&mesh, &|p| 300.0 + p[0] + 2.0 * p[1] + current_t.powi(2));
        let older = values(&mesh, &|p| {
            300.0 + p[0] + 2.0 * p[1] + (current_t - old_dt).powi(2)
        });
        let snapshots = (previous.clone(), older.clone());
        let boundaries = prescribed(&mesh, &target);
        let prepared = closure
            .prepare_scalar_transport(
                &velocity,
                &molecular,
                &ratio,
                &boundaries,
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        for cell in 0..nc {
            assert!(prepared.turbulent_diffusivity_m2_s[cell] > 0.0);
            close(
                prepared.turbulent_diffusivity_m2_s[cell],
                prepared.sgs.cells[cell].eddy_viscosity_m2_s / ratio[cell],
                0.0,
            );
            close(
                prepared.transport.diffusivity_m2_s()[cell],
                molecular[cell] + prepared.turbulent_diffusivity_m2_s[cell],
                0.0,
            );
        }
        let source = |_: usize, p| {
            let u = flow(p);
            2.0 * (current_t + dt) + u[0] + 2.0 * u[1]
        };
        let step = prepared
            .transport
            .bdf2_variable_step(
                &previous,
                Some((&older, old_dt)),
                dt,
                &source,
                controls(),
                || false,
            )
            .unwrap();
        for (actual, exact) in step.values.iter().zip(values(&mesh, &target)) {
            close(*actual, exact, 2e-8);
        }
        close(step.storage_rate, 2.0 * (current_t + dt), 2e-8);
        close(step.balance_residual, 0.0, 2e-8);
        assert_eq!((previous, older), snapshots);
        // Different authored Pr_t changes heat mixing without changing the
        // original SGS momentum stress or silently converting it into heat.
        let wider_ratio = closure
            .prepare_scalar_transport(
                &velocity,
                &molecular,
                &vec![1.4; nc],
                &boundaries,
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        for cell in 0..nc {
            close(
                wider_ratio.turbulent_diffusivity_m2_s[cell],
                0.5 * prepared.turbulent_diffusivity_m2_s[cell],
                0.0,
            );
            close(
                wider_ratio.sgs.cells[cell].eddy_viscosity_m2_s,
                prepared.sgs.cells[cell].eddy_viscosity_m2_s,
                0.0,
            );
        }
        // A nonuniform thermal pulse actually mixes faster with the active
        // modeled diffusivity; coefficient-only diagnostics cannot pass this.
        let cold = |_: [f64; 2]| 0.0;
        let cold_boundaries = prescribed(&mesh, &cold);
        let pulse = values(&mesh, &|p| {
            (std::f64::consts::PI * p[0]).sin() * (std::f64::consts::PI * p[1]).sin()
        });
        let mixing = closure
            .prepare_scalar_transport(
                &velocity,
                &molecular,
                &ratio,
                &cold_boundaries,
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        let molecular_only = system
            .prepare_scalar_transport(
                &velocity,
                &molecular,
                &cold_boundaries,
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        let turbulent_step = mixing
            .transport
            .bdf1_step(&pulse, 0.15, &|_, _| 0.0, controls(), || false)
            .unwrap();
        let molecular_step = molecular_only
            .bdf1_step(&pulse, 0.15, &|_, _| 0.0, controls(), || false)
            .unwrap();
        assert!(molecular_step.integral - turbulent_step.integral > 1e-5);
        close(turbulent_step.balance_residual, 0.0, 2e-10);
        // Positive modeled momentum dissipation is not an implicit heat source.
        let uniform = |_: [f64; 2]| 300.0;
        let no_heating = closure
            .prepare_scalar_transport(
                &velocity,
                &molecular,
                &ratio,
                &prescribed(&mesh, &uniform),
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        assert!(no_heating.sgs.integrated_dissipation_m4_s3 > 0.0);
        let unchanged = no_heating
            .transport
            .bdf1_step(&vec![300.0; 3 * nc], 0.15, &|_, _| 0.0, controls(), || {
                false
            })
            .unwrap();
        for value in unchanged.values {
            close(value, 300.0, 2e-8);
        }
    }
}

#[test]
fn scalar_sip_diffusion_converges_on_skew_triangles_and_preserves_closed_inventory() {
    let pi = std::f64::consts::PI;
    let shear = 0.35;
    let mut errors = Vec::new();
    for resolution in [2, 4, 8] {
        let mut mesh2 = Mesh2::triangles(1.0, 1.0, resolution, resolution);
        for p in &mut mesh2.nodes {
            p[0] += shear * p[1];
        }
        let mesh = TriMesh::from_mesh2(&mesh2);
        let system = FluxSystem::new(&mesh);
        let zero = |_: [f64; 2]| 0.0;
        let exact = |p: [f64; 2]| (pi * (p[0] - shear * p[1])).sin() * (pi * p[1]).sin();
        let source = |_: usize, p: [f64; 2]| {
            let x = pi * (p[0] - shear * p[1]);
            let y = pi * p[1];
            (5.0 + (2.0 + shear * shear) * pi * pi) * x.sin() * y.sin()
                + 2.0 * shear * pi * pi * x.cos() * y.cos()
        };
        let native = system
            .prepare_scalar_transport(
                &vec![0.0; system.n_u],
                &vec![1.0; mesh.tris.len()],
                &prescribed(&mesh, &zero),
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        let step = native
            .bdf1_step(
                &vec![0.0; 3 * mesh.tris.len()],
                0.2,
                &source,
                controls(),
                || false,
            )
            .unwrap();
        let mut error = 0.0;
        for (cell, tri) in mesh.tris.iter().enumerate() {
            let points = tri.map(|v| mesh.verts[v]);
            for (point, weight) in fs_flux::bdm::tri_quad(points, mesh.areas[cell]) {
                let actual = p1_value(
                    points,
                    mesh.areas[cell],
                    &step.values[3 * cell..3 * cell + 3],
                    point,
                );
                error += weight * (actual - exact(point)).powi(2);
            }
        }
        errors.push(error.sqrt());
        assert!(step.balance_residual.abs() < 1e-9);
        let closed = system
            .prepare_scalar_transport(
                &vec![0.0; system.n_u],
                &vec![0.1; mesh.tris.len()],
                &insulated(&mesh, &zero, None),
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        let initial = values(&mesh, &|p| 2.0 + (pi * (p[0] - shear * p[1])).cos());
        let after = closed
            .bdf1_step(&initial, 0.5, &|_, _| 0.0, controls(), || false)
            .unwrap();
        close(after.integral, 2.0, 2e-11);
        close(after.storage_rate, 0.0, 2e-11);
        close(after.balance_residual, 0.0, 2e-11);
        assert!(
            after
                .values
                .iter()
                .map(|v| (v - 2.0).abs())
                .fold(0.0_f64, f64::max)
                < 0.9
        );
    }
    assert!(
        errors[0] / errors[1] > 2.7 && errors[1] / errors[2] > 2.7,
        "skew-grid L2 errors: {errors:?}"
    );
}

fn p1_value(p: [[f64; 2]; 3], area: f64, values: &[f64], x: [f64; 2]) -> f64 {
    let cross = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
    (0..3)
        .map(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            values[i]
                * cross(
                    [p[j][0] - x[0], p[j][1] - x[1]],
                    [p[k][0] - x[0], p[k][1] - x[1]],
                )
                / (2.0 * area)
        })
        .sum()
}

#[test]
fn reversing_open_faces_use_both_scalar_donors_and_close_the_physical_budget() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let flow = |p: [f64; 2]| [2.0 * p[1] - 1.0, 1.0 - 2.0 * p[0]];
    let velocity = moments(&mesh, &flow);
    let zero = |_: [f64; 2]| 0.0;
    let reservoir = |p: [f64; 2]| 5.0 + p[1];
    let boundaries = insulated(&mesh, &zero, Some(&reservoir));
    let plan = system
        .prepare_scalar_transport(
            &velocity,
            &[0.02, 0.02],
            &boundaries,
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    let previous = values(&mesh, &|p| 1.0 + 2.0 * p[0] + 0.5 * p[1]);
    let step = plan
        .bdf1_step(&previous, 0.1, &|_, _| 0.0, controls(), || false)
        .unwrap();
    for account in &step.boundaries {
        let edge = &mesh.edges[account.edge];
        let a = mesh.verts[edge.verts.0];
        let b = mesh.verts[edge.verts.1];
        let un = |p| {
            let u = flow(p);
            u[0] * edge.normal[0] + u[1] * edge.normal[1]
        };
        let u0 = un(a);
        let du = un(b) - u0;
        assert!(
            u0 * un(b) < 0.0,
            "each boundary reverses with zero net volume"
        );
        let split = -u0 / du;
        let cell = edge.tris.0;
        let points = mesh.tris[cell].map(|v| mesh.verts[v]);
        let local = &step.values[3 * cell..3 * cell + 3];
        let internal = [
            p1_value(points, mesh.areas[cell], local, a),
            p1_value(points, mesh.areas[cell], local, b),
        ];
        let external = [reservoir(a), reservoir(b)];
        let mut expected = 0.0;
        for (lo, hi) in [(0.0, split), (split, 1.0)] {
            let c = if u0 + du * (0.5 * (lo + hi)) > 0.0 {
                internal
            } else {
                external
            };
            let dc = c[1] - c[0];
            let antiderivative = |s: f64| {
                u0 * c[0] * s + (u0 * dc + du * c[0]) * s * s / 2.0 + du * dc * s * s * s / 3.0
            };
            expected += edge.len * (antiderivative(hi) - antiderivative(lo));
        }
        close(account.outward_advection, expected, 2e-12);
        close(account.outward_diffusion, 0.0, 0.0);
    }
    assert!(step.storage_rate > 1.0);
    close(step.balance_residual, 0.0, 2e-11);
    assert!(matches!(
        system.prepare_scalar_transport(
            &velocity,
            &[0.02, 0.02],
            &insulated(&mesh, &zero, None),
            ScalarTransportControls::default(),
            || false
        ),
        Err(ScalarTransportError::MissingInflow { .. })
    ));
}

#[test]
fn harmonic_diffusion_resolves_material_contrast_without_centroid_alignment() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 2));
    let system = FluxSystem::new(&mesh);
    let diffusivity: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| if p[1] < 0.5 { 0.1 } else { 10.0 })
        .collect();
    let target = |p: [f64; 2]| {
        if p[1] <= 0.5 {
            p[1] / 0.1
        } else {
            5.0 + (p[1] - 0.5) / 10.0
        }
    };
    let old = values(&mesh, &target);
    let plan = system
        .prepare_scalar_transport(
            &vec![0.0; system.n_u],
            &diffusivity,
            &prescribed(&mesh, &target),
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    let step = plan
        .bdf1_step(&old, 0.4, &|_, _| 0.0, controls(), || false)
        .unwrap();
    for (actual, exact) in step.values.iter().zip(old) {
        close(*actual, exact, 2e-9);
    }
    close(step.storage_rate, 0.0, 2e-9);
    close(step.balance_residual, 0.0, 2e-9);
}

#[test]
fn disabled_sgs_matches_native_scalar_and_refusal_never_changes_histories() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let flow = moments(&mesh, &|p| [0.2 + p[1], 0.0]);
    let field = |p: [f64; 2]| 2.0 + p[0] - p[1];
    let boundaries = prescribed(&mesh, &field);
    let old = values(&mesh, &field);
    let snapshot = old.clone();
    let nc = mesh.tris.len();
    let native = system
        .prepare_scalar_transport(
            &flow,
            &vec![0.03; nc],
            &boundaries,
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    let expected = native
        .bdf1_step(&old, 0.1, &|_, _| 0.0, controls(), || false)
        .unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let sgs = PlanarSgs::with_model(&system, model, &vec![0.2; nc], &vec![0.01; nc]).unwrap();
        let prepared = sgs
            .prepare_scalar_transport(
                &flow,
                &vec![0.03; nc],
                &vec![0.7; nc],
                &boundaries,
                ScalarTransportControls::default(),
                || false,
            )
            .unwrap();
        let actual = prepared
            .transport
            .bdf2_variable_step(&old, None, 0.1, &|_, _| 0.0, controls(), || false)
            .unwrap();
        assert_eq!(actual.values, expected.values);
        assert_eq!(
            actual.balance_residual.to_bits(),
            expected.balance_residual.to_bits()
        );
        assert!(
            prepared
                .turbulent_diffusivity_m2_s
                .iter()
                .all(|&v| v == 0.0)
        );
    }
    let mut polls = 0;
    assert!(matches!(
        native.bdf1_step(&old, 0.1, &|_, _| 0.0, controls(), || {
            polls += 1;
            polls == 120
        }),
        Err(ScalarTransportError::Cancelled)
    ));
    assert!(
        native
            .bdf2_variable_step(
                &old,
                Some((&old, 0.01)),
                0.1,
                &|_, _| 0.0,
                controls(),
                || false
            )
            .is_err()
    );
    assert!(
        native
            .bdf1_step(&old, 0.1, &|_, _| f64::NAN, controls(), || false)
            .is_err()
    );
    assert!(
        native
            .bdf1_step(
                &old,
                0.1,
                &|_, _| 0.0,
                ScalarSolveControls {
                    maximum_iterations: 0,
                    ..controls()
                },
                || false
            )
            .is_err()
    );
    assert_eq!(old, snapshot);
    let constant = |_p: [f64; 2]| 7.0;
    let no_heat = system
        .prepare_scalar_transport(
            &flow,
            &vec![0.0; nc],
            &prescribed(&mesh, &constant),
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap();
    let uniform = no_heat
        .bdf1_step(&vec![7.0; 3 * nc], 0.3, &|_, _| 0.0, controls(), || false)
        .unwrap();
    for value in uniform.values {
        close(value, 7.0, 2e-11);
    }
    close(uniform.storage_rate, 0.0, 2e-10);
}
