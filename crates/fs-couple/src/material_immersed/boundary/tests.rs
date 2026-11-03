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
