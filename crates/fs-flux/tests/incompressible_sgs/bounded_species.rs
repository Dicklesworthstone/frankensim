//! Independent inventories, simplex preservation, signed reservoirs and SGS.

use super::*;
use fs_flux::ns::transport::scalar::bounded::species::{
    BoundedSpeciesError, BoundedSpeciesStep, PreparedBoundedSpeciesTransport, SPECIES_SUM_ROUNDOFF,
};
use fs_flux::ns::transport::scalar::bounded::{
    BoundedScalarControls, BoundedScalarError, ScalarBounds,
};
use fs_flux::ns::transport::scalar::{
    PreparedScalarTransport, ScalarBoundary, ScalarTransportControls,
};

fn prepare<'m>(
    system: &FluxSystem<'m>,
    mesh: &'m TriMesh,
    velocity: &[f64],
    diffusivity: &[f64],
    donor: &dyn Fn([f64; 2]) -> f64,
    value_boundary: bool,
) -> PreparedScalarTransport<'m> {
    let zero = |_: [f64; 2]| 0.0;
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| {
            if value_boundary {
                ScalarBoundary::Value { edge, value: donor }
            } else {
                ScalarBoundary::Flux {
                    edge,
                    outward: &zero,
                    inflow: Some(donor),
                }
            }
        })
        .collect();
    system
        .prepare_scalar_transport(
            velocity,
            diffusivity,
            &boundaries,
            ScalarTransportControls::default(),
            || false,
        )
        .unwrap()
}

fn pure(mesh: &TriMesh, species: usize, phase: impl Fn([f64; 2]) -> usize) -> Vec<Vec<f64>> {
    let mut values = vec![vec![0.0; 3 * mesh.tris.len()]; species];
    for (cell, &p) in mesh.centroids.iter().enumerate() {
        values[phase(p)][3 * cell..3 * cell + 3].fill(1.0);
    }
    values
}
fn references(values: &[Vec<f64>]) -> Vec<&[f64]> {
    values.iter().map(Vec::as_slice).collect()
}
fn audit(mesh: &TriMesh, previous: &[Vec<f64>], result: &BoundedSpeciesStep) {
    assert!(result.maximum_sum_residual <= SPECIES_SUM_ROUNDOFF);
    for (s, field) in result.species.iter().enumerate() {
        assert!(field.minimum >= -3e-14 && field.maximum <= 1.0 + 3e-14);
        // Recompute EVERY species inventory from original triangle areas and
        // accepted physical/numerical exchanges, independently of the receipt.
        let before: f64 = previous[s]
            .iter()
            .enumerate()
            .map(|(i, v)| v * mesh.areas[i / 3] / 3.0)
            .sum();
        let after: f64 = field
            .values
            .iter()
            .enumerate()
            .map(|(i, v)| v * mesh.areas[i / 3] / 3.0)
            .sum();
        let outward: f64 = field
            .boundaries
            .iter()
            .map(|b| b.advection + b.diffusion + b.limiter_diffusion)
            .sum();
        close(field.integral_before, before, 3e-13);
        close(field.integral_after, after, 3e-13);
        close(after - before + outward, 0.0, 3e-12);
        close(field.balance_residual, 0.0, 3e-12);
    }
    for dof in 0..previous[0].len() {
        close(
            result.species.iter().map(|s| s.values[dof]).sum(),
            1.0,
            3e-13,
        );
    }
}

#[test]
fn common_limiting_preserves_the_simplex_where_independent_scalars_fail() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 8, 4));
    let system = FluxSystem::new(&mesh);
    let velocity = moments(&mesh, &|_| [1.0, 0.35]);
    let diffusivity = vec![1e-8; mesh.tris.len()];
    let native: Vec<_> = [0.2, 0.3, 0.5]
        .into_iter()
        .map(|value| prepare(&system, &mesh, &velocity, &diffusivity, &|_| value, false))
        .collect();
    let previous = pure(&mesh, 3, |p| {
        if p[0] < 0.4 {
            0
        } else if p[1] < 0.5 {
            1
        } else {
            2
        }
    });
    let snapshot = previous.clone();
    let dt = 0.12;
    let independent: Vec<_> = native
        .iter()
        .zip(&previous)
        .map(|(op, values)| {
            op.prepare_bounded_transport(
                ScalarBounds {
                    lower: 0.0,
                    upper: 1.0,
                },
                || false,
            )
            .unwrap()
            .ssprk2_step(values, dt, BoundedScalarControls::default(), || false)
            .unwrap()
        })
        .collect();
    let defect = (0..previous[0].len())
        .map(|i| (independent.iter().map(|s| s.values[i]).sum::<f64>() - 1.0).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        defect > 1e-5,
        "independent-scalar negative control too weak: {defect:e}"
    );
    let bound =
        PreparedBoundedSpeciesTransport::new(&native.iter().collect::<Vec<_>>(), || false).unwrap();
    let result = bound
        .ssprk2_step(
            &references(&previous),
            dt,
            BoundedScalarControls::default(),
            || false,
        )
        .unwrap();
    audit(&mesh, &previous, &result);
    assert!(
        result
            .species
            .iter()
            .all(|s| s.substeps > 1 && s.limited_fluxes > 0)
    );
    assert!(
        result
            .species
            .iter()
            .any(|s| (s.integral_after - s.integral_before).abs() > 0.01)
    );
    assert_eq!(previous, snapshot);

    // All species are independent. Renaming or moving a component to the last
    // slot cannot change the accepted physical trajectory or its own budget.
    let order = [2, 0, 1];
    let permuted = PreparedBoundedSpeciesTransport::new(&order.map(|s| &native[s]), || false)
        .unwrap()
        .ssprk2_step(
            &order.map(|s| previous[s].as_slice()),
            dt,
            BoundedScalarControls::default(),
            || false,
        )
        .unwrap();
    for (slot, original) in order.into_iter().enumerate() {
        assert_eq!(
            permuted.species[slot].values,
            result.species[original].values
        );
        assert_eq!(
            permuted.species[slot].balance_residual,
            result.species[original].balance_residual
        );
    }
}

#[test]
fn closed_skew_grid_mixing_conserves_every_species_including_an_absent_component() {
    let mut geometry = Mesh2::triangles(1.0, 1.0, 4, 3);
    for p in &mut geometry.nodes {
        p[0] += 0.7 * p[1];
    }
    let mesh = TriMesh::from_mesh2(&geometry);
    let system = FluxSystem::new(&mesh);
    let diffusivity: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| if p[1] < 0.5 { 0.002 } else { 0.2 })
        .collect();
    let native = prepare(
        &system,
        &mesh,
        &vec![0.0; system.n_u],
        &diffusivity,
        &|_| 0.0,
        false,
    );
    let bound = PreparedBoundedSpeciesTransport::new(&[&native; 5], || false).unwrap();
    let previous = pure(&mesh, 5, |p| {
        usize::from(p[0] - 0.7 * p[1] > 0.5) + 2 * usize::from(p[1] > 0.5)
    });
    let result = bound
        .ssprk2_step(
            &references(&previous),
            0.2,
            BoundedScalarControls::default(),
            || false,
        )
        .unwrap();
    audit(&mesh, &previous, &result);
    for field in &result.species {
        close(field.integral_after, field.integral_before, 3e-12);
        assert!(
            field
                .boundaries
                .iter()
                .all(|b| b.advection == 0.0 && b.diffusion == 0.0 && b.limiter_diffusion == 0.0)
        );
    }
    assert!(result.species[..4].iter().all(|s| s.limited_fluxes > 0));
    assert!(
        result.species[0]
            .values
            .iter()
            .zip(&previous[0])
            .any(|(a, b)| (a - b).abs() > 0.05)
    );
    assert!(result.species[4].values.iter().all(|&v| v == 0.0));
}

#[test]
fn all_sgs_models_share_signed_reservoir_corrections_without_net_species_diffusion() {
    let mut geometry = Mesh2::triangles(1.0, 1.0, 2, 2);
    for p in &mut geometry.nodes {
        p[0] += 2.0 * p[1];
    }
    let mesh = TriMesh::from_mesh2(&geometry);
    let system = FluxSystem::new(&mesh);
    // Retain resolved reversing flow and positive SGS while exposing the
    // signed SIP load instead of masking it with dominant inflow advection.
    let velocity = moments(&mesh, &|p| [0.01 * (2.0 * p[1] - 0.7), 0.01 * (0.8 - p[0])]);
    let volume = system.edge_volume_transport(&velocity).unwrap();
    assert!(
        volume
            .iter()
            .enumerate()
            .any(|(edge, v)| mesh.is_boundary(edge)
                && v.along_normal_m2_s() > 0.0
                && v.against_normal_m2_s() > 0.0)
    );
    let donor: [fn([f64; 2]) -> f64; 3] = [
        |p| 0.05 + 0.2 * p[1],
        |p| 0.8 - 0.5 * p[1],
        |p| 0.15 + 0.3 * p[1],
    ];
    // A pure carrier makes signed weak reservoir loads try to remove an
    // initially absent component. The common reservoir limiter must act.
    let previous = pure(&mesh, 3, |_| 0);
    for model in [
        SgsModel::Smagorinsky { cs: 0.3 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.1 },
    ] {
        let closure = PlanarSgs::with_model(
            &system,
            model,
            &vec![0.2; mesh.tris.len()],
            &vec![0.01; mesh.tris.len()],
        )
        .unwrap();
        let coupled: Vec<_> = donor
            .iter()
            .map(|value| {
                let boundaries: Vec<_> = mesh
                    .edges
                    .iter()
                    .enumerate()
                    .filter(|(edge, _)| mesh.is_boundary(*edge))
                    .map(|(edge, _)| ScalarBoundary::Value { edge, value })
                    .collect();
                closure
                    .prepare_scalar_transport(
                        &velocity,
                        &vec![0.02; mesh.tris.len()],
                        &vec![0.7; mesh.tris.len()],
                        &boundaries,
                        ScalarTransportControls::default(),
                        || false,
                    )
                    .unwrap()
            })
            .collect();
        assert!(
            coupled[0]
                .transport
                .diffusivity_m2_s()
                .iter()
                .all(|&d| d > 0.02)
        );
        let bound = PreparedBoundedSpeciesTransport::new(
            &coupled.iter().map(|s| &s.transport).collect::<Vec<_>>(),
            || false,
        )
        .unwrap();
        let dt = 0.02;
        let result = bound
            .ssprk2_step(
                &references(&previous),
                dt,
                BoundedScalarControls::default(),
                || false,
            )
            .unwrap();
        audit(&mesh, &previous, &result);
        assert!(result.species.iter().any(|s| {
            s.boundaries
                .iter()
                .any(|b| b.limiter_diffusion.abs() > 1e-12)
        }));
        for i in 0..result.species[0].boundaries.len() {
            let edge = result.species[0].boundaries[i].edge;
            let advection: f64 = result
                .species
                .iter()
                .map(|s| s.boundaries[i].advection)
                .sum();
            let diffusion: f64 = result
                .species
                .iter()
                .map(|s| s.boundaries[i].diffusion)
                .sum();
            let numerical: f64 = result
                .species
                .iter()
                .map(|s| s.boundaries[i].limiter_diffusion)
                .sum();
            close(advection, dt * volume[edge].net_m2_s(), 3e-13);
            close(diffusion, 0.0, 3e-13);
            close(numerical, 0.0, 3e-13);
        }
    }
}

#[test]
fn smooth_unlimited_species_recover_the_original_consistent_mass_scalar_evolution() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = moments(&mesh, &|_| [0.2, -0.1]);
    let donor: [fn([f64; 2]) -> f64; 3] = [
        |p| 0.2 + 0.05 * p[0],
        |p| 0.3 + 0.04 * p[1],
        |p| 0.5 - 0.05 * p[0] - 0.04 * p[1],
    ];
    let native: Vec<_> = donor
        .iter()
        .map(|f| {
            prepare(
                &system,
                &mesh,
                &velocity,
                &vec![0.01; mesh.tris.len()],
                f,
                true,
            )
        })
        .collect();
    let previous: Vec<Vec<_>> = donor
        .iter()
        .map(|f| {
            mesh.tris
                .iter()
                .flat_map(|tri| tri.map(|v| f(mesh.verts[v])))
                .collect()
        })
        .collect();
    let dt = 1e-5;
    let result = PreparedBoundedSpeciesTransport::new(&native.iter().collect::<Vec<_>>(), || false)
        .unwrap()
        .ssprk2_step(
            &references(&previous),
            dt,
            BoundedScalarControls::default(),
            || false,
        )
        .unwrap();
    audit(&mesh, &previous, &result);
    for s in 0..3 {
        let single = native[s]
            .prepare_bounded_transport(
                ScalarBounds {
                    lower: 0.0,
                    upper: 1.0,
                },
                || false,
            )
            .unwrap()
            .ssprk2_step(&previous[s], dt, BoundedScalarControls::default(), || false)
            .unwrap();
        assert_eq!(result.species[s].limited_fluxes, 0);
        assert_eq!(result.species[s].substeps, 1);
        assert_eq!(single.values, result.species[s].values);
        assert_eq!(single.balance_residual, result.species[s].balance_residual);
    }
    // The original private scalar regression independently compares this
    // inactive path with original 3x3 mass-matrix LU, not a lumped target.
}

#[test]
fn mismatched_physics_invalid_compositions_and_late_cancellation_are_atomic() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = moments(&mesh, &|_| [1.0, 0.0]);
    let native: Vec<_> = [0.2, 0.3, 0.5]
        .into_iter()
        .map(|y| {
            prepare(
                &system,
                &mesh,
                &velocity,
                &vec![0.01; mesh.tris.len()],
                &|_| y,
                false,
            )
        })
        .collect();
    let wrong = prepare(
        &system,
        &mesh,
        &velocity,
        &vec![0.02; mesh.tris.len()],
        &|_| 0.3,
        false,
    );
    assert!(matches!(
        PreparedBoundedSpeciesTransport::new(&[&native[0], &wrong, &native[2]], || false),
        Err(BoundedSpeciesError::OperatorMismatch { species: 1 })
    ));
    let bad_sum = prepare(
        &system,
        &mesh,
        &velocity,
        &vec![0.01; mesh.tris.len()],
        &|_| 0.6,
        false,
    );
    assert!(matches!(
        PreparedBoundedSpeciesTransport::new(&[&native[0], &native[1], &bad_sum], || false),
        Err(BoundedSpeciesError::BoundarySum { .. })
    ));
    let bound =
        PreparedBoundedSpeciesTransport::new(&native.iter().collect::<Vec<_>>(), || false).unwrap();
    let previous = pure(&mesh, 3, |p| {
        if p[0] < 0.4 {
            0
        } else if p[1] < 0.5 {
            1
        } else {
            2
        }
    });
    let snapshot = previous.clone();
    let mut polls = 0;
    let complete = bound
        .ssprk2_step(
            &references(&previous),
            0.15,
            BoundedScalarControls::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut stopped = 0;
    assert!(matches!(
        bound.ssprk2_step(
            &references(&previous),
            0.15,
            BoundedScalarControls::default(),
            || {
                stopped += 1;
                stopped > polls * 3 / 4
            }
        ),
        Err(BoundedSpeciesError::Cancelled)
    ));
    let retry = bound
        .ssprk2_step(
            &references(&previous),
            0.15,
            BoundedScalarControls::default(),
            || false,
        )
        .unwrap();
    for s in 0..3 {
        assert_eq!(complete.species[s].values, retry.species[s].values);
    }
    assert_eq!(previous, snapshot);
    let mut bad_state = previous.clone();
    bad_state[0][0] = 0.7;
    assert!(matches!(
        bound.ssprk2_step(
            &references(&bad_state),
            0.15,
            BoundedScalarControls::default(),
            || false
        ),
        Err(BoundedSpeciesError::StateSum { .. })
    ));
    assert!(matches!(
        bound.ssprk2_step(
            &references(&previous),
            0.15,
            BoundedScalarControls {
                maximum_substeps: 1,
                ..BoundedScalarControls::default()
            },
            || false
        ),
        Err(BoundedSpeciesError::Scalar {
            error: BoundedScalarError::SubstepLimit { .. },
            ..
        })
    ));

    // Only the left edge has the adversarial affine trace. Its two Gauss
    // samples form valid simplex compositions, but its physical top endpoint
    // has Y_1=1.05 and Y_2=-0.05 and MUST be rejected before any transport.
    let coarse = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&coarse);
    let profiles: [fn([f64; 2]) -> f64; 3] = [
        |_| 0.0,
        |p| 0.9 + 0.3 * (p[1] - 0.5),
        |p| 0.1 - 0.3 * (p[1] - 0.5),
    ];
    let constants: [fn([f64; 2]) -> f64; 3] = [|_| 0.0, |_| 0.5, |_| 0.5];
    let adversarial: Vec<_> = (0..3)
        .map(|s| {
            let boundaries: Vec<_> = coarse
                .edges
                .iter()
                .enumerate()
                .filter(|(e, _)| coarse.is_boundary(*e))
                .map(|(edge, face)| ScalarBoundary::Value {
                    edge,
                    value: if face.mid[0] == 0.0 {
                        &profiles[s]
                    } else {
                        &constants[s]
                    },
                })
                .collect();
            system
                .prepare_scalar_transport(
                    &vec![0.0; system.n_u],
                    &vec![0.01; coarse.tris.len()],
                    &boundaries,
                    ScalarTransportControls::default(),
                    || false,
                )
                .unwrap()
        })
        .collect();
    assert!(matches!(
        PreparedBoundedSpeciesTransport::new(&adversarial.iter().collect::<Vec<_>>(), || false),
        Err(BoundedSpeciesError::Scalar {
            error: BoundedScalarError::BoundaryBounds { .. },
            ..
        })
    ));
}
