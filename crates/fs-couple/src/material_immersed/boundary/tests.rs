//! Conversion and native forwarding use original transport and fluid owners.
use super::*;
use fs_lbm::d3q19::{Face3, PartialSaturationBody3};
use fs_material::fluid_transport::FluidTransportState;

fn domain() -> ImmersedDomain {
    ImmersedDomain {
        cells: [8; 3],
        periodic: [true; 3],
        cell_size_m: 1.0 / 64.0,
        timestep_s: 1.0 / 8.0,
        maximum_speed_m_s: 1.0 / 512.0,
        maximum_cells: 512,
        samples_per_axis: 2,
    }
}

fn fluid() -> FluidTransportState {
    // Authored SI coefficients, not material-card or experimental authority.
    FluidTransportState::try_new(2.0, 1024.0, 0.25, 1.0 / 4096.0).unwrap()
}

fn faces() -> [FaceBoundarySi; 6] {
    [
        FaceBoundarySi::Periodic,
        FaceBoundarySi::Periodic,
        FaceBoundarySi::StationaryWall,
        FaceBoundarySi::StationaryWall,
        FaceBoundarySi::Velocity {
            velocity_m_s: [0.0, 0.0, 1.0 / 2048.0],
        },
        FaceBoundarySi::Pressure {
            perturbation_pa: 0.0,
        },
    ]
}

#[test]
fn si_pressure_and_velocity_reuse_original_units_and_legacy_defaults() {
    let fluid = fluid();
    for dx in [1.0 / 64.0, 1.0 / 32.0] {
        let domain = ImmersedDomain {
            cell_size_m: dx,
            ..domain()
        };
        let scaling =
            plan_scaling_from_material(&fluid, dx, domain.timestep_s, domain.maximum_speed_m_s)
                .unwrap();
        let speed_unit = dx / domain.timestep_s;
        let pressure_unit = fs_lbm::CS2 * 2.0 * speed_unit * speed_unit;
        for delta in [-0.5, -0.003, 0.0, 0.005, 0.5] {
            let pressure = pressure_unit * delta;
            let density = scaling
                .density_lattice_from_pressure_perturbation(pressure)
                .unwrap();
            assert!((density - (1.0 + delta)).abs() < 3e-16);
            assert!(
                (scaling.pressure_perturbation_si(density).unwrap() - pressure).abs()
                    < 3e-16 * pressure_unit
            );
        }
        assert_eq!(
            scaling
                .density_lattice_from_pressure_perturbation(0.0)
                .unwrap(),
            1.0
        );
        let mut faces = faces();
        faces[5] = FaceBoundarySi::Pressure {
            perturbation_pa: pressure_unit / 1024.0,
        };
        let config = domain.with_boundaries(faces);
        assert_eq!(config.faces(), &faces);
        let native = config.lattice_boundaries(&scaling).unwrap();
        assert_eq!(
            native.face(Face3::ZMax),
            FaceBoundary3::Pressure {
                density: 1.0 + 1.0 / 1024.0
            }
        );
        assert_eq!(
            native.face(Face3::ZMin),
            FaceBoundary3::Velocity {
                velocity: [0.0, 0.0, (1.0 / 2048.0) / speed_unit]
            }
        );
        assert_eq!(native.face(Face3::XMin), FaceBoundary3::Periodic);
        assert_eq!(native.face(Face3::YMin), FaceBoundary3::stationary_wall());
        let legacy = ImmersedBoundaryDomain::from(domain)
            .lattice_boundaries(&scaling)
            .unwrap();
        assert_eq!(legacy, BoundarySpec3::periodic());
        let mixed = ImmersedDomain {
            periodic: [true, false, false],
            ..domain
        };
        let legacy = ImmersedBoundaryDomain::from(mixed)
            .lattice_boundaries(&scaling)
            .unwrap();
        assert_eq!(legacy.face(Face3::XMin), FaceBoundary3::Periodic);
        assert_eq!(legacy.face(Face3::YMin), FaceBoundary3::stationary_wall());
    }
}

#[test]
fn si_faces_refuse_unrepresentable_pressure_speed_and_native_topology() {
    let fluid = fluid();
    let domain = domain();
    let scaling = plan_scaling_from_material(
        &fluid,
        domain.cell_size_m,
        domain.timestep_s,
        domain.maximum_speed_m_s,
    )
    .unwrap();
    let pressure_unit = fs_lbm::CS2 * 2.0 * (domain.cell_size_m / domain.timestep_s).powi(2);
    for pressure in [
        f64::NAN,
        f64::INFINITY,
        -pressure_unit,
        -2.0 * pressure_unit,
        f64::MIN_POSITIVE,
        pressure_unit * f64::EPSILON * 0.25,
    ] {
        let mut faces = faces();
        faces[5] = FaceBoundarySi::Pressure {
            perturbation_pa: pressure,
        };
        assert!(matches!(
            domain.with_boundaries(faces).lattice_boundaries(&scaling),
            Err(Error::Scaling(_))
        ));
    }
    let mut too_fast = faces();
    too_fast[4] = FaceBoundarySi::Velocity {
        velocity_m_s: [0.0, 0.0, 2.0 * domain.maximum_speed_m_s],
    };
    assert!(matches!(
        domain
            .with_boundaries(too_fast)
            .lattice_boundaries(&scaling),
        Err(Error::Invalid(_))
    ));
    let mut unpaired = faces();
    unpaired[0] = FaceBoundarySi::StationaryWall;
    assert!(matches!(
        domain
            .with_boundaries(unpaired)
            .lattice_boundaries(&scaling),
        Err(Error::Flow(_))
    ));
    let mut multiple_axes = faces();
    multiple_axes[0] = FaceBoundarySi::Pressure {
        perturbation_pa: 0.0,
    };
    multiple_axes[1] = FaceBoundarySi::StationaryWall;
    assert!(matches!(
        domain
            .with_boundaries(multiple_axes)
            .lattice_boundaries(&scaling),
        Err(Error::Flow(_))
    ));
    // Raw native admission is the same owner used by the legacy assert wrapper.
    for velocity in [[f64::NAN, 0.0, 0.0], [0.0, 0.01, 0.0], [0.2, 0.0, 0.0]] {
        assert!(BoundarySpec3::lid_cavity(velocity).try_validate().is_err());
    }
    assert!(
        BoundarySpec3::lid_cavity([0.01, 0.0, 0.0])
            .try_validate()
            .is_ok()
    );
}

#[test]
fn si_faces_drive_identical_native_fluid_body_and_reservoir_momentum() {
    let fluid = fluid();
    let domain = domain();
    let scaling = plan_scaling_from_material(
        &fluid,
        domain.cell_size_m,
        domain.timestep_s,
        domain.maximum_speed_m_s,
    )
    .unwrap();
    let from_si = domain
        .with_boundaries(faces())
        .lattice_boundaries(&scaling)
        .unwrap();
    let authored_native = BoundarySpec3::new([
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::stationary_wall(),
        FaceBoundary3::stationary_wall(),
        FaceBoundary3::Velocity {
            velocity: [0.0, 0.0, 1.0 / 256.0],
        },
        FaceBoundary3::Pressure { density: 1.0 },
    ]);
    assert_eq!(from_si, authored_native);
    let mut actual = BoundaryGrid3::new(8, 8, 8, scaling.flow().tau, [0.0; 3], from_si);
    let mut reference = BoundaryGrid3::new(8, 8, 8, scaling.flow().tau, [0.0; 3], authored_native);
    let mut body = PartialSaturationBody3 {
        center: [4.0; 3],
        mass: 2.0,
        inertia_world: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        velocity: [0.0; 3],
        angular_velocity: [0.0; 3],
    };
    let mut expected_body = body;
    let mut fractions = vec![0.0; 512];
    fractions[(3 * 8 + 3) * 8 + 3] = 0.4;
    let gate = CancelGate::new();
    let mut reservoir_impulse = 0.0;
    for _ in 0..6 {
        let actual_receipt = actual
            .step_with_partial_saturation_body(
                &fractions,
                &mut body,
                [0.0; 6],
                [-2.0, 1.0, 3.0],
                &gate,
            )
            .unwrap();
        let expected_receipt = reference
            .step_with_partial_saturation_body(
                &fractions,
                &mut expected_body,
                [0.0; 6],
                [-2.0, 1.0, 3.0],
                &gate,
            )
            .unwrap();
        assert_eq!(actual_receipt, expected_receipt);
        assert_eq!(body, expected_body);
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..8 {
                    assert_eq!(actual.populations(x, y, z), reference.populations(x, y, z));
                }
            }
        }
        assert!(actual_receipt.open.reconstructed_cells > 0);
        reservoir_impulse += actual_receipt.open.impulse[2].abs();
    }
    assert!(reservoir_impulse > 1e-3);
    assert!(
        body.velocity[2] > 0.0,
        "the inlet must accelerate the actual immersed body"
    );
}

fn profile_faces(axis: usize, upper: bool) -> [FaceBoundarySi; 6] {
    let mut faces = [FaceBoundarySi::Periodic; 6];
    faces[2 * axis + usize::from(upper)] = FaceBoundarySi::Velocity {
        velocity_m_s: [0.0; 3],
    };
    faces[2 * axis + usize::from(!upper)] = FaceBoundarySi::Pressure {
        perturbation_pa: 0.0,
    };
    faces
}

fn profile_velocity(position: [f64; 3], axis: usize, sign: f64) -> [f64; 3] {
    let mut value = [0.0; 3];
    value[axis] = -sign * (1.0 / 2048.0 + position[(axis + 1) % 3] / 256.0);
    value[(axis + 1) % 3] = position[(axis + 2) % 3] / 1024.0;
    value
}

fn populations(grid: &BoundaryGrid3) -> Vec<[f64; fs_lbm::d3q19::Q3]> {
    (0..512)
        .map(|i| grid.populations(i % 8, (i / 8) % 8, i / 64))
        .collect()
}

#[test]
fn si_profiles_use_physical_cell_centers_and_outward_mean_on_every_face() {
    let fluid = fluid();
    let domain = domain();
    let gate = CancelGate::new_clock_free();
    let scaling = plan_scaling_from_material(
        &fluid, domain.cell_size_m, domain.timestep_s, domain.maximum_speed_m_s,
    )
    .unwrap();
    let speed_unit = domain.cell_size_m / domain.timestep_s;
    let mean = -1.0 / 1024.0;
    for axis in 0..3 {
        for upper in [false, true] {
            let face = Face3::ALL[2 * axis + usize::from(upper)];
            let sign = f64::from(face.normal()[axis]);
            let boundary = domain.with_boundaries(profile_faces(axis, upper))
                .lattice_boundaries(&scaling).unwrap();
            let mut actual = BoundaryGrid3::new(8, 8, 8, scaling.flow().tau, [0.0; 3], boundary);
            let mut reference = BoundaryGrid3::new(8, 8, 8, scaling.flow().tau, [0.0; 3], boundary);
            let before = populations(&actual);
            let mut positions = Vec::new();
            let report = set_velocity_profile_si(
                &mut actual, &scaling, domain.maximum_speed_m_s, face,
                InletFluxSi::MeanNormalVelocity(mean), &gate,
                |position| {
                    positions.push(position);
                    Ok(profile_velocity(position, axis, sign))
                },
            )
            .unwrap();
            // Independent dimensions, without the material conversion helpers.
            let expected = reference.set_velocity_profile_with_speed_limit(
                face, InletFlux3::MeanNormal(mean / speed_unit),
                domain.maximum_speed_m_s / speed_unit, &gate,
                |cell| Ok(profile_velocity(
                    cell.map(|i| (i as f64 + 0.5) * domain.cell_size_m), axis, sign,
                ).map(|v| v / speed_unit)),
            )
            .unwrap();
            assert_eq!(report, expected);
            assert_eq!(actual.velocity_profile(face), reference.velocity_profile(face));
            assert_eq!(positions.len(), 64);
            let mut first = [0.5 * domain.cell_size_m; 3];
            let mut last = [7.5 * domain.cell_size_m; 3];
            first[axis] = if upper { 7.5 } else { 0.5 } * domain.cell_size_m;
            last[axis] = first[axis];
            assert_eq!(positions[0], first);
            assert_eq!(positions[63], last);
            assert!((sign * report.prescribed_mean_velocity_lu[axis] * speed_unit - mean).abs() < 1e-17);
            assert_eq!(populations(&actual), before);
            actual.step();
            reference.step();
            assert_eq!(populations(&actual), populations(&reference));
        }
    }
}

#[test]
fn si_profile_corrected_overspeed_source_failure_and_cancel_keep_previous_trace() {
    let fluid = fluid();
    let domain = domain();
    let gate = CancelGate::new_clock_free();
    let scaling = plan_scaling_from_material(
        &fluid, domain.cell_size_m, domain.timestep_s, domain.maximum_speed_m_s,
    )
    .unwrap();
    let boundary = domain.with_boundaries(profile_faces(2, false))
        .lattice_boundaries(&scaling).unwrap();
    let mut grid = BoundaryGrid3::new(8, 8, 8, scaling.flow().tau, [0.0; 3], boundary);
    let maximum = domain.maximum_speed_m_s;
    set_velocity_profile_si(
        &mut grid, &scaling, maximum, Face3::ZMin, InletFluxSi::Unconstrained,
        &gate, |_| Ok([0.0, 0.0, maximum / 4.0]),
    )
    .unwrap();
    let profile = grid.velocity_profile(Face3::ZMin).unwrap().to_vec();
    let before = populations(&grid);
    // The raw tangent and requested normal mean are each below the bound.
    // Their corrected vector is not: the native owner must refuse before swap.
    let error = set_velocity_profile_si(
        &mut grid, &scaling, maximum, Face3::ZMin,
        InletFluxSi::MeanNormalVelocity(-0.75 * maximum), &gate,
        |_| Ok([0.75 * maximum, 0.0, 0.0]),
    )
    .unwrap_err();
    assert!(matches!(error, Error::Profile(BoundaryProfileError3::InvalidVelocity { .. })));
    assert_eq!(grid.velocity_profile(Face3::ZMin).unwrap(), profile);
    let reason = BoundaryProfileError3::Source("coherent inlet field refused".into());
    let mut calls = 0;
    assert_eq!(
        set_velocity_profile_si(
            &mut grid, &scaling, maximum, Face3::ZMin, InletFluxSi::Unconstrained,
            &gate, |_| {
                calls += 1;
                if calls == 20 { Err(reason.clone()) } else { Ok([0.0; 3]) }
            },
        ),
        Err(Error::Profile(reason)),
    );
    assert!(matches!(
        set_velocity_profile_si(
            &mut grid, &scaling, maximum, Face3::ZMin, InletFluxSi::Unconstrained,
            &gate, |_| Ok([f64::NAN, 0.0, 0.0]),
        ),
        Err(Error::Profile(BoundaryProfileError3::Source(_))),
    ));
    let cancelled = CancelGate::new_clock_free();
    assert_eq!(
        set_velocity_profile_si(
            &mut grid, &scaling, maximum, Face3::ZMin, InletFluxSi::Unconstrained,
            &cancelled, |_| {
                cancelled.request();
                Ok([0.0; 3])
            },
        ),
        Err(Error::Profile(BoundaryProfileError3::Cancelled)),
    );
    assert_eq!(grid.velocity_profile(Face3::ZMin).unwrap(), profile);
    assert_eq!(populations(&grid), before);
    // Failed proposals left the same two-buffer owner reusable.
    set_velocity_profile_si(
        &mut grid, &scaling, maximum, Face3::ZMin,
        InletFluxSi::MeanNormalVelocity(-maximum / 4.0), &gate,
        |_| Ok([0.0, 0.0, maximum / 2.0]),
    )
    .unwrap();
    assert_eq!(grid.velocity_profile(Face3::ZMin).unwrap(), profile);
    assert_eq!(populations(&grid), before);
}
