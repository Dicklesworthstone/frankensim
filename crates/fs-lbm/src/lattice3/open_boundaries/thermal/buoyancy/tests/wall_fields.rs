//! Spatial wall forcing through the original native scalar and buoyancy owners.

use super::*;
use crate::thermal::boundary::ThermalWallField;

fn wall_grid<P: BoussinesqForceLaw3>(
    axis: usize,
    initial: f64,
    law: P,
) -> BoussinesqOpenGrid3<D3q19, 19, P> {
    let mut native = [FaceBoundary3::Periodic; 6];
    let mut thermal = [ThermalOpenFace3::Periodic; 6];
    for side in 0..2 {
        native[2 * axis + side] = FaceBoundary3::stationary_wall();
        thermal[2 * axis + side] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
    }
    make_with_policy(
        CollisionModel3::Bgk { tau: 0.9 },
        [4, 3, 2],
        BoundarySpec3::new(native),
        thermal,
        initial,
        law,
    )
}

fn affine_grid(axis: usize) -> BoussinesqOpenGrid3<D3q19, 19> {
    wall_grid(
        axis,
        0.5,
        BoussinesqLaw3::try_new(REFERENCE_RHO, 0.5, BETA, GRAVITY, 0.02).unwrap(),
    )
}

// Independent inverse of the public plane ordering, without the native helper.
fn cell_on_face(dims: [usize; 3], axis: usize, side: usize, link: usize) -> usize {
    let xyz = match axis {
        0 => [side * (dims[0] - 1), link % dims[1], link / dims[1]],
        1 => [link % dims[0], side * (dims[1] - 1), link / dims[0]],
        2 => [link % dims[0], link / dims[0], side * (dims[2] - 1)],
        _ => unreachable!(),
    };
    xyz[0] + dims[0] * (xyz[1] + dims[1] * xyz[2])
}

#[test]
fn spatial_wall_laws_match_independent_reflection_and_storage_on_all_six_faces() {
    let dims = [4, 3, 2];
    let weight = 0.125;
    for axis in 0..3 {
        let mut owner = affine_grid(axis);
        let links = 24 / dims[axis];
        let mut expected_temperatures = [0.5; 24];
        let mut heat = [0.0; 2];
        let mut traces = [Vec::new(), Vec::new()];
        let population_address = owner.thermal().populations().as_ptr();
        for side in 0..2 {
            let face = 2 * axis + side;
            let fields: Vec<_> = (0..links)
                .map(|link| match (link + side) % 4 {
                    0 => ThermalWall::Temperature(0.55 + 0.002 * link as f64),
                    1 => ThermalWall::HeatFlux(-0.001 * (link + 1) as f64),
                    2 => ThermalWall::Convective {
                        transfer: 0.03 + 0.001 * link as f64,
                        ambient: 0.8,
                    },
                    _ => ThermalWall::Adiabatic,
                })
                .collect();
            let wall_address = owner.thermal().wall_states(face).unwrap().as_ptr();
            owner
                .set_wall_boundary_field(face, ThermalWallField::Spatial(&fields))
                .unwrap();
            assert_eq!(owner.thermal().wall_states(face).unwrap(), fields);
            assert_eq!(
                owner.thermal().wall_states(face).unwrap().as_ptr(),
                wall_address
            );
            // Original initial equilibrium, independent half-link heat relation.
            let signed_u = (2.0 * side as f64 - 1.0) * INITIAL_VELOCITY[axis];
            let outgoing = weight * 0.5 * (1.0 + signed_u / 0.25);
            let old_incoming = weight * 0.5 * (1.0 - signed_u / 0.25);
            for (link, wall) in fields.into_iter().enumerate() {
                let q = match wall {
                    ThermalWall::Temperature(t) => 2.0 * weight * t - 2.0 * outgoing,
                    ThermalWall::HeatFlux(q) => q,
                    ThermalWall::Convective { transfer, ambient } => {
                        (ambient - outgoing / weight) / (1.0 / transfer + 0.5 / weight)
                    }
                    ThermalWall::Adiabatic => 0.0,
                    ThermalWall::FiniteSolid(_) => unreachable!(),
                };
                heat[side] += q;
                traces[side].push((2.0 * outgoing + q) / (2.0 * weight));
                expected_temperatures[cell_on_face(dims, axis, side, link)] +=
                    outgoing + q - old_incoming;
            }
        }
        assert_eq!(owner.thermal().populations().as_ptr(), population_address);
        assert_eq!(owner.thermal().completed_steps(), 0);
        assert!(owner.last_step().is_none());
        let report = advance(&mut owner, &[0.0; 24]);
        for side in 0..2 {
            let wall = report.heat.walls[2 * axis + side].unwrap();
            near(wall.heat_in, heat[side], 5e-14);
            near(
                wall.mean_temperature,
                traces[side].iter().sum::<f64>() / links as f64,
                5e-14,
            );
            near(
                wall.min_temperature,
                traces[side].iter().copied().fold(f64::INFINITY, f64::min),
                5e-14,
            );
            near(
                wall.max_temperature,
                traces[side]
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max),
                5e-14,
            );
            assert!(wall.solid.is_none());
        }
        for (cell, expected) in expected_temperatures.into_iter().enumerate() {
            near(owner.thermal().temperature(cell).unwrap(), expected, 5e-14);
        }
        near(report.heat.stored_change, heat.into_iter().sum(), 2e-13);
        near(report.heat.combined_defect, 0.0, 2e-13);
    }
}

#[test]
fn uniform_and_spatial_wall_inputs_share_native_history_and_atomic_retry() {
    let mut actual = affine_grid(1);
    let mut reference = affine_grid(1);
    let fields = vec![ThermalWall::Temperature(0.6); 8];
    actual
        .set_wall_boundary_field(2, ThermalWallField::Spatial(&fields))
        .unwrap();
    reference
        .set_wall_boundary_field(2, ThermalWallField::Uniform(fields[0]))
        .unwrap();
    for _ in 0..3 {
        assert_eq!(
            advance(&mut actual, &[0.0; 24]),
            advance(&mut reference, &[0.0; 24])
        );
        assert_eq!(accepted(&actual), accepted(&reference));
    }
    let fields: Vec<_> = (0..8)
        .map(|link| ThermalWall::HeatFlux(0.0001 * (link as f64 - 3.5)))
        .collect();
    for owner in [&mut actual, &mut reference] {
        owner
            .set_wall_boundary_field(2, ThermalWallField::Spatial(&fields))
            .unwrap();
    }
    let before = accepted(&actual);
    let mut bad = fields.clone();
    bad[7] = ThermalWall::Temperature(f64::NAN);
    assert!(
        actual
            .set_wall_boundary_field(2, ThermalWallField::Spatial(&bad))
            .is_err()
    );
    assert!(
        actual
            .set_wall_boundary_field(2, ThermalWallField::Spatial(&fields[..7]))
            .is_err()
    );
    assert!(
        actual
            .set_wall_boundary_field(0, ThermalWallField::Spatial(&fields))
            .is_err()
    );
    bad[7] = ThermalWall::FiniteSolid(FiniteThermalWall {
        temperature: 0.5,
        heat_capacity: 1.0,
        resistance: 0.0,
        heat_per_step: 0.0,
    });
    assert!(
        actual
            .set_wall_boundary_field(2, ThermalWallField::Spatial(&bad))
            .is_err()
    );
    assert_eq!(accepted(&actual), before);

    let mut finite = channel(BETA);
    let old_finite = accepted(&finite);
    let replacement = vec![ThermalWall::Adiabatic; 18];
    assert!(
        finite
            .set_wall_boundary_field(2, ThermalWallField::Spatial(&replacement))
            .is_err()
    );
    assert_eq!(accepted(&finite), old_finite);
    assert!(finite.thermal().wall_conduction(1).is_some());

    let mut count = 0;
    let expected = reference
        .step_with_poll(None, &[0.0002; 24], None, [0.0; 3], &mut || {
            count += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    let failed = actual.step_with_poll(None, &[0.0002; 24], None, [0.0; 3], &mut || {
        calls += 1;
        if calls == count {
            Err(GridError3::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        failed,
        Err(BoussinesqOpenError3::Flow(GridError3::Cancelled))
    ));
    assert_eq!(calls, count);
    assert_eq!(accepted(&actual), before);
    assert_eq!(advance(&mut actual, &[0.0002; 24]), expected);
    assert_eq!(accepted(&actual), accepted(&reference));

    // The original uniform setter intentionally replaces the complete field.
    actual.set_wall_boundary(2, ThermalWall::Adiabatic).unwrap();
    assert!(
        actual
            .thermal()
            .wall_states(2)
            .unwrap()
            .iter()
            .all(|wall| *wall == ThermalWall::Adiabatic)
    );
}

#[derive(Clone, Copy)]
struct DensityPeak;
impl DensityDepartureLaw3 for DensityPeak {
    type Error = core::convert::Infallible;
    fn relative_density_departure(&self, theta: f64) -> Result<f64, Self::Error> {
        Ok(0.02 / (1.0 + 40000.0 * (theta - 0.5).powi(2)))
    }
}

fn nonlinear_grid() -> BoussinesqOpenGrid3<D3q19, 19, BoussinesqDensityLaw3<DensityPeak>> {
    wall_grid(
        1,
        0.2,
        BoussinesqDensityLaw3::try_new(REFERENCE_RHO, GRAVITY, 0.01, DensityPeak).unwrap(),
    )
}

#[test]
fn nonlinear_buoyancy_admits_actual_wall_fields_before_binding_and_every_resolved_trace() {
    let mut actual = nonlinear_grid();
    let mut prescribed = vec![ThermalWall::Temperature(0.4); 8];
    prescribed[7] = ThermalWall::Temperature(0.5);
    let initial = accepted(&actual);
    assert!(matches!(
        actual.set_wall_boundary_field(2, ThermalWallField::Spatial(&prescribed)),
        Err(BoussinesqOpenError3::BoundaryLaw { face: 2, .. })
    ));
    assert_eq!(accepted(&actual), initial);

    let BoussinesqOpenGrid3 {
        flow,
        mut thermal,
        law,
        background,
        limits,
        ..
    } = nonlinear_grid();
    // Source-free thermal input is legal, but attaching its actual source law
    // must inspect the installed per-link data, not its adiabatic declaration.
    thermal
        .set_wall_boundary_field(2, ThermalWallField::Spatial(&prescribed))
        .unwrap();
    let rejected = BoussinesqOpenGrid3::try_new_with_poll(
        flow.inner,
        thermal,
        law,
        background,
        limits,
        &mut || Ok(()),
    );
    assert!(matches!(
        rejected,
        Err(BoussinesqOpenError3::BoundaryLaw { face: 2, .. })
    ));

    // Inward heat fluxes reconstruct surface traces 0.4, 0.5 and 0.6.
    // Both extrema satisfy the nonlinear density envelope; the middle does not.
    assert!(actual.force_policy().sample(0.4, [0.0; 3]).is_ok());
    assert!(actual.force_policy().sample(0.6, [0.0; 3]).is_ok());
    let fields: Vec<_> = (0..8)
        .map(|link| ThermalWall::HeatFlux(0.05 + 0.025 * (link % 3) as f64))
        .collect();
    actual
        .set_wall_boundary_field(2, ThermalWallField::Spatial(&fields))
        .unwrap();
    let before = accepted(&actual);
    let rejected = actual.step_with_poll(None, &[0.0; 24], None, [0.0; 3], &mut || Ok(()));
    assert!(matches!(
        rejected,
        Err(BoussinesqOpenError3::BoundaryLaw {
            face: 2,
            source: BoussinesqDensityLawError3::Law(BoussinesqLawError3::DensityDeparture { .. }),
        })
    ));
    assert_eq!(accepted(&actual), before);
    let safe = vec![ThermalWall::HeatFlux(0.05); 8];
    actual
        .set_wall_boundary_field(2, ThermalWallField::Spatial(&safe))
        .unwrap();
    let report = advance(&mut actual, &[0.0; 24]);
    near(report.heat.walls[2].unwrap().heat_in, 8.0 * 0.05, 3e-14);
    near(report.heat.combined_defect, 0.0, 3e-14);
}

#[test]
fn inactive_wall_flux_and_late_invalid_data_cannot_partially_replace_inputs() {
    let mut solid = vec![false; 8];
    solid[0] = true;
    let faces = [
        ThermalOpenFace3::Periodic,
        ThermalOpenFace3::Periodic,
        ThermalOpenFace3::Wall(ThermalWall::Adiabatic),
        ThermalOpenFace3::Wall(ThermalWall::Adiabatic),
        ThermalOpenFace3::Periodic,
        ThermalOpenFace3::Periodic,
    ];
    let mut field = ThermalOpenField3::try_new::<ThermalOpenError3>(
        [2, 2, 2],
        &solid,
        &[0.5; 8],
        0.8,
        ThermalOpenBoundaries3 {
            faces,
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 {
            maximum_cells: 8,
            maximum_storage_bytes: 8192,
            temperature_range: [-1.0, 2.0],
            balance_tolerance: 1e-12,
        },
        0,
        |_| Ok([0.0; 3]),
        || Ok(()),
    )
    .unwrap();
    let mut laws = vec![ThermalWall::HeatFlux(0.0); 4];
    laws[3] = ThermalWall::Convective {
        transfer: 0.04,
        ambient: 3.0,
    };
    field
        .set_wall_boundary_field(2, ThermalWallField::Spatial(&laws))
        .unwrap();
    let previous = field.wall_states(2).unwrap().to_vec();
    let address = field.wall_states(2).unwrap().as_ptr();
    let heat = field.populations().to_vec();
    laws[0] = ThermalWall::HeatFlux(0.01);
    assert!(
        field
            .set_wall_boundary_field(2, ThermalWallField::Spatial(&laws))
            .is_err()
    );
    assert_eq!(field.wall_states(2).unwrap(), previous);
    laws[0] = ThermalWall::HeatFlux(0.0);
    laws[3] = ThermalWall::Convective {
        transfer: -1.0,
        ambient: 0.5,
    };
    assert!(
        field
            .set_wall_boundary_field(2, ThermalWallField::Spatial(&laws))
            .is_err()
    );
    assert_eq!(field.wall_states(2).unwrap(), previous);
    assert_eq!(field.wall_states(2).unwrap().as_ptr(), address);
    assert_eq!(field.populations(), heat);
    assert_eq!(field.completed_steps(), 0);
    assert!(field.last_balance().is_none());
}
