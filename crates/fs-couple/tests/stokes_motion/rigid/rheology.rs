//! G1/G3/G4: material rheology participates in the same three-momentum solve.

use super::*;
use fs_flux::MaterialFluxRheologyPoint;
use fs_material::MaterialError;
use fs_material::fluid_transport::FluidTransport;
use std::cell::Cell;

// Explicit numerical Carreau fixtures; these coefficients assert no sourced
// liquid/model authority. The reference snapshot is the zero-shear state.
fn viscosity(rate: f64, exponent: f64) -> f64 {
    0.06 + (0.8 - 0.06) * (1.0 + (1.8 * rate).powi(2)).powf(0.5 * (exponent - 1.0))
}

struct Response {
    point: MaterialFluxRheologyPoint,
    ordinal: usize,
    transport: FluidTransportState,
}
impl FluidTransport for Response {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}

/// Independently manufacture changing rigid spin/translation PLUS changing
/// affine shear. The finite-slip stress is mu(gamma)*gamma at the new field.
/// Constant initial viscosity or an old-shear viscosity gives the wrong load.
#[test]
fn g1_g3_rheological_rotating_body_resolves_actual_shear_force_torque_and_slip() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let (rho, dt) = (2.0, 0.03);
    let material = FluidTransportState::try_new(rho, 1000.0, 0.5, 0.8).unwrap();
    let pivot = [0.35, 0.4];
    let selected: Vec<_> = edges(&mesh)
        .into_iter()
        .filter(|&e| mesh.edges[e].normal[1] < -0.5)
        .collect();
    let (gamma0, gamma1) = (0.1, 0.65);
    let (omega0, omega1) = (0.1, 0.4);
    let angle: f64 = dt * omega1;
    let rotation = [angle.cos(), angle.sin()];
    let mut cfg = config();
    cfg.translation_stiffness_n_m = [0.8, 1.3];
    cfg.rotation_stiffness_n_m_rad = 0.4;
    cfg.translation_damping_n_s_m = [0.07, 0.09];
    cfg.rotation_damping_n_m_s_rad = 0.05;

    for (exponent, beta, boost) in [
        (0.55, None, [0.0; 2]),
        (0.55, Some(0.9), [0.6, -0.3]),
        (1.3, Some(0.9), [0.0; 2]),
    ] {
        let mu = viscosity(gamma1, exponent);
        let tau = mu * gamma1;
        let slip0 = beta.map_or(0.0, |b| viscosity(gamma0, exponent) * gamma0 / b);
        let slip1 = beta.map_or(0.0, |b| tau / b);
        let v0 = [0.05 + boost[0], -0.03 + boost[1]];
        let v1 = [-0.02 + boost[0], 0.04 + boost[1]];
        let end = [pivot[0] + dt * v1[0], pivot[1] + dt * v1[1]];
        let initial_field = |y: [f64; 2]| {
            let w = spin(omega0, [y[0] - pivot[0], y[1] - pivot[1]]);
            [v0[0] + w[0] + gamma0 * y[1] + slip0, v0[1] + w[1]]
        };
        let initial = history(&mesh, initial_field);
        let v1_ref = inverse(rotation, v1);
        let exact = |y: [f64; 2]| {
            let w = spin(omega1, [y[0] - pivot[0], y[1] - pivot[1]]);
            [v1_ref[0] + w[0] + gamma1 * y[1] + slip1, v1_ref[1] + w[1]]
        };
        let force = |x: [f64; 2]| {
            let r = inverse(rotation, [x[0] - end[0], x[1] - end[1]]);
            let y = [pivot[0] + r[0], pivot[1] + r[1]];
            let current = exact(y);
            let old = inverse(rotation, initial_field(y));
            // Uniform shear-dependent stress has zero divergence. This is the
            // independent covariant time derivative plus grad(u)*(u-w_mesh).
            rotate(
                rotation,
                [
                    rho * (current[0] - old[0]) / dt,
                    rho * ((current[1] - old[1]) / dt + omega1 * (gamma1 * y[1] + slip1)),
                ],
            )
        };
        let traction = |_: [f64; 2], n: [f64; 2]| {
            let n = inverse(rotation, n);
            rotate(rotation, [tau * n[1], tau * n[0]])
        };
        let exterior: Vec<_> = edges(&mesh)
            .into_iter()
            .filter(|e| !selected.contains(e))
            .map(|edge| MaterialNavierStokesBoundary::Open {
                edge,
                traction: &traction,
                backflow: OpenBoundaryPolicy::NaturalTraction,
            })
            .collect();
        let native = body(pivot, v0, omega0, 0.37);
        let mut owner = match beta {
            None => PlanarRigidMotion::new(
                &system,
                &material,
                MaterialFluxNumerics::default(),
                mass(),
                native,
                &initial,
                &selected,
                cfg,
            )
            .unwrap(),
            Some(b) => PlanarRigidMotion::new_with_slip(
                &system,
                &material,
                MaterialFluxNumerics::default(),
                mass(),
                native,
                &initial,
                &selected,
                cfg,
                b,
            )
            .unwrap(),
        };
        let physical_force = rotate(rotation, [cfg.depth_m * tau, 0.0]);
        let physical_torque = cfg.depth_m * pivot[1] * tau;
        let external: [f64; 2] = core::array::from_fn(|i| {
            mass().mass() * (v1[i] - v0[i]) / dt - physical_force[i]
                + cfg.translation_stiffness_n_m[i] * dt * v1[i]
                + cfg.translation_damping_n_s_m[i] * v1[i]
        });
        let external_torque = 0.3 * (omega1 - omega0) / dt - physical_torque
            + cfg.rotation_stiffness_n_m_rad * angle
            + cfg.rotation_damping_n_m_s_rad * omega1;
        let calls = Cell::new(0);
        let step = owner
            .step_with_rheology(
                dt,
                &force,
                &exterior,
                &[None],
                external,
                external_torque,
                fluid_controls(),
                controls(),
                |point, reference| {
                    calls.set(calls.get() + 1);
                    assert_eq!(*reference, material);
                    Ok(Response {
                        point,
                        ordinal: calls.get(),
                        transport: FluidTransportState::try_new(
                            reference.density_kg_m3(),
                            reference.specific_heat_j_kg_k(),
                            reference.thermal_conductivity_w_m_k(),
                            viscosity(point.shear_rate_per_s, exponent),
                        )?,
                    })
                },
                || false,
            )
            .unwrap();
        for i in 0..2 {
            near(step.velocity_m_s[i], v1[i], 4e-7);
            near(step.fluid_force_n[i], physical_force[i], 4e-7);
            near(step.flow.pivot_world_m[i], end[i], 3e-8);
        }
        near(step.angular_velocity_rad_s, omega1, 4e-7);
        near(step.fluid_torque_n_m, physical_torque, 4e-7);
        near(step.rotation_rad, angle, 3e-8);
        let (error, divergence) =
            system.velocity_error(step.flow.flow.solution.velocity_dofs(), &exact);
        assert!(error < 4e-7 && divergence < 4e-7, "{error}, {divergence}");
        near(
            step.interface_work_j,
            dt * (physical_force[0] * v1[0] + physical_force[1] * v1[1] + physical_torque * omega1),
            4e-8,
        );
        near(
            step.interface_slip_dissipation_j,
            beta.map_or(0.0, |b| dt * cfg.depth_m * tau * tau / b),
            4e-8,
        );
        near(step.mechanical_work_residual_j, 0.0, 4e-8);
        let flux = step.flow.boundary_transport_world(&selected, end).unwrap();
        near(flux.outward_mass_flux_kg_per_m_s, 0.0, 4e-8);

        let evaluated = calls.get();
        let samples = step.flow.rheology_samples();
        assert_eq!(samples.len(), mesh.tris.len());
        assert!(
            evaluated > samples.len(),
            "must resolve evolving candidate strain"
        );
        for (cell, sample) in samples.iter().enumerate() {
            assert_eq!(sample.point().cell, cell);
            assert_eq!(sample.material().point, sample.point());
            assert!(sample.material().ordinal > 0 && sample.material().ordinal <= evaluated);
            near(sample.point().shear_rate_per_s, gamma1, 4e-6);
            near(sample.transport().dynamic_viscosity_pa_s(), mu, 4e-6);
            let centroid = mesh.centroids[cell];
            let r = rotate(rotation, [centroid[0] - pivot[0], centroid[1] - pivot[1]]);
            for axis in 0..2 {
                near(sample.point().position_m[axis], end[axis] + r[axis], 5e-8);
            }
        }
        assert_eq!(
            calls.get(),
            evaluated,
            "retained owners must not re-query the source"
        );
    }
}

fn resolved(
    point: MaterialFluxRheologyPoint,
    reference: &FluidTransportState,
) -> Result<FluidTransportState, MaterialError> {
    FluidTransportState::try_new(
        reference.density_kg_m3(),
        reference.specific_heat_j_kg_k(),
        reference.thermal_conductivity_w_m_k(),
        viscosity(point.shear_rate_per_s, 0.55),
    )
}

/// A physical uniform-flow added-mass oracle is independent of constitutive
/// iteration. Both public paths must preserve the same native mechanical owner.
#[test]
fn g3_rheological_zero_strain_keeps_newtonian_added_mass_and_source_identity() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let rho = 2.0;
    let material = FluidTransportState::try_new(rho, 1000.0, 0.5, 0.8).unwrap();
    let selected = edges(&mesh);
    let v0 = [0.04, -0.02];
    let initial = history(&mesh, |_| v0);
    let native = body([0.5, 0.5], v0, 0.0, 0.0);
    let cfg = config();
    let make = || {
        PlanarRigidMotion::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            mass(),
            native,
            &initial,
            &selected,
            cfg,
        )
        .unwrap()
    };
    let (mut original, mut rheological) = (make(), make());
    let (dt, external) = (0.03, [0.5, -0.2]);
    let a = original
        .step(
            dt,
            &zero,
            &[],
            &[Some(0.0)],
            external,
            0.0,
            fluid_controls(),
            controls(),
            || false,
        )
        .unwrap();
    let b = rheological
        .step_with_rheology(
            dt,
            &zero,
            &[],
            &[Some(0.0)],
            external,
            0.0,
            fluid_controls(),
            controls(),
            resolved,
            || false,
        )
        .unwrap();
    for axis in 0..2 {
        let expected = v0[axis] + dt * external[axis] / (mass().mass() + rho * cfg.depth_m);
        near(b.velocity_m_s[axis], expected, 2e-7);
        near(a.velocity_m_s[axis], b.velocity_m_s[axis], 2e-8);
        near(a.fluid_force_n[axis], b.fluid_force_n[axis], 2e-8);
    }
    near(b.angular_velocity_rad_s, 0.0, 2e-7);
    near(b.fluid_torque_n_m, 0.0, 2e-7);
    assert_eq!(b.flow.rheology_samples().len(), mesh.tris.len());
    assert!(a.flow.rheology_samples().is_empty());
    assert_eq!(rheological.material(), &material);
    // A later constant-viscosity step uses the same accepted flow history.
    let a = original
        .step(
            dt,
            &zero,
            &[],
            &[Some(0.0)],
            external,
            0.0,
            fluid_controls(),
            controls(),
            || false,
        )
        .unwrap();
    let b = rheological
        .step(
            dt,
            &zero,
            &[],
            &[Some(0.0)],
            external,
            0.0,
            fluid_controls(),
            controls(),
            || false,
        )
        .unwrap();
    for axis in 0..2 {
        near(a.velocity_m_s[axis], b.velocity_m_s[axis], 4e-8);
        near(a.flow.pivot_world_m[axis], b.flow.pivot_world_m[axis], 4e-8);
    }
    near(b.flow.time_s, 2.0 * dt, 1e-14);
    assert!(b.flow.rheology_samples().is_empty());
}

/// Source refusal, changed density and terminal cancellation must not publish
/// any part of a coupled endpoint; a pure retry reproduces a fresh owner.
#[test]
fn g4_rheology_refusal_and_final_cancel_preserve_body_flow_and_retry() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let material = FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.8).unwrap();
    let selected = edges(&mesh);
    let initial = vec![0.0; system.n_u];
    let native = body([0.5, 0.5], [0.0; 2], 0.0, 0.0);
    let make = || {
        PlanarRigidMotion::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            mass(),
            native,
            &initial,
            &selected,
            config(),
        )
        .unwrap()
    };
    let (mut owner, mut reference) = (make(), make());
    let before = owner.flow_state().clone();
    for change_density in [false, true] {
        assert!(
            owner
                .step_with_rheology(
                    0.03,
                    &zero,
                    &[],
                    &[Some(0.0)],
                    [0.5, -0.2],
                    0.0,
                    fluid_controls(),
                    controls(),
                    |_, _| FluidTransportState::try_new(
                        if change_density { 3.0 } else { 2.0 },
                        1000.0,
                        0.5,
                        if change_density { 0.8 } else { -1.0 },
                    ),
                    || false,
                )
                .is_err()
        );
        assert_eq!(owner.body(), native);
        assert_eq!(owner.flow_state(), &before);
        assert_eq!(owner.rotation_rad(), 0.0);
    }
    let mut polls = 0;
    let expected = reference
        .step_with_rheology(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.0,
            fluid_controls(),
            controls(),
            resolved,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut observed = 0;
    assert!(matches!(
        owner.step_with_rheology(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.0,
            fluid_controls(),
            controls(),
            resolved,
            || {
                observed += 1;
                observed == polls
            },
        ),
        Err(PlanarRigidError::Cancelled)
    ));
    assert_eq!(observed, polls);
    assert_eq!(owner.body(), native);
    assert_eq!(owner.flow_state(), &before);
    assert_eq!(owner.rotation_rad(), 0.0);
    let actual = owner
        .step_with_rheology(
            0.03,
            &zero,
            &[],
            &[Some(0.0)],
            [0.5, -0.2],
            0.0,
            fluid_controls(),
            controls(),
            resolved,
            || false,
        )
        .unwrap();
    assert_eq!(owner.body(), reference.body());
    assert_eq!(owner.flow_state(), reference.flow_state());
    assert_eq!(owner.rotation_rad(), reference.rotation_rad());
    assert_eq!(actual.fluid_force_n, expected.fluid_force_n);
    assert_eq!(actual.fluid_torque_n_m, expected.fluid_torque_n_m);
    for (a, b) in actual
        .flow
        .rheology_samples()
        .iter()
        .zip(expected.flow.rheology_samples())
    {
        assert_eq!(a.point(), b.point());
        assert_eq!(a.transport(), b.transport());
        assert_eq!(a.material(), b.material());
    }
}
