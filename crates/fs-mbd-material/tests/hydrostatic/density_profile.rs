//! G1/G3/G4: sampled density drives the original forces, moments and motion.
use super::*;
use fs_mbd_material::hydrostatic::{
    DensityProfileHydrostaticMedium, HydrostaticDensitySample, HydrostaticStepControls,
};

fn samples<'a>(
    depths: &[f64],
    fluids: &'a [FluidTransportState],
) -> Vec<HydrostaticDensitySample<'a>> {
    depths
        .iter()
        .zip(fluids)
        .map(|(&depth_m, source)| HydrostaticDensitySample { depth_m, source })
        .collect()
}

// Independent Cartesian integrals through a rectangular column. Production
// instead clips tetrahedra and integrates positive barycentric moments.
fn column(a: f64, b: f64, profile: &[(f64, f64)]) -> (f64, f64, f64) {
    let (mut mass, mut first, mut pressure, mut head) = (0.0, 0.0, 0.0, 0.0);
    for (i, &(d, rho)) in profile.iter().enumerate() {
        let (end, slope) = profile
            .get(i + 1)
            .map_or((b, 0.0), |&(end, next)| (end, (next - rho) / (end - d)));
        let (lo, hi) = (a.max(d), b.min(end));
        if hi > lo {
            let intercept = rho - slope * d;
            mass += intercept * (hi - lo) + 0.5 * slope * (hi * hi - lo * lo);
            first +=
                0.5 * intercept * (hi * hi - lo * lo) + slope * (hi.powi(3) - lo.powi(3)) / 3.0;
            let (s, t) = (lo - d, hi - d);
            pressure += head * (hi - lo)
                + 0.5 * rho * (t * t - s * s)
                + slope * (t.powi(3) - s.powi(3)) / 6.0;
        }
        head += rho * (end - d) + 0.5 * slope * (end - d).powi(2);
    }
    (mass, first, pressure)
}

#[test]
fn multiple_thermoclines_match_independent_column_mass_center_and_pressure() {
    let depths = [0.0, 0.25, 1.0, 3.0];
    let densities = [1000.0, 1100.0, 1500.0, 1550.0];
    let fluids = densities.map(fluid);
    let medium = DensityProfileHydrostaticMedium::new(
        &samples(&depths, &fluids),
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let source = body(1300.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let profile: Vec<_> = depths.into_iter().zip(densities).collect();
    for z in [2.0_f64, 1.0, 0.6, 0.0, -0.9, -2.0, -4.0] {
        let (a, b) = ((-z - 1.0).max(0.0), (1.0 - z).max(0.0));
        let (mass, first, pressure) = column(a, b, &profile);
        let state = body.state_from_kinematics(pose(z), Twist::zero()).unwrap();
        let load = medium.load(&cube(), &body, state, || false).unwrap();
        close(load.displaced().volume_m3(), 4.0 * (b - a));
        close(load.displaced_mass_kg(), 4.0 * mass);
        close(load.wrench().force_world.z, 40.0 * mass);
        close(load.pressure_potential_j(), 40.0 * pressure);
        if mass > 0.0 {
            vec_close(
                load.center_of_buoyancy_world().unwrap(),
                Vec3::new(0.0, 0.0, -first / mass),
            );
        } else {
            assert_eq!(load.center_of_buoyancy_world(), None);
        }
    }
    for (depth, expected) in [
        (0.0, 1000.0),
        (0.125, 1050.0),
        (0.625, 1300.0),
        (2.0, 1525.0),
        (100.0, 1550.0),
    ] {
        close(medium.density_at_depth_kg_m3(depth).unwrap(), expected);
    }
}

#[test]
fn profile_force_and_torque_are_the_same_continuous_potential_derivatives() {
    let fluids = [fluid(800.0), fluid(1600.0), fluid(1700.0)];
    let medium = DensityProfileHydrostaticMedium::new(
        &samples(&[0.0, 0.4, 2.0], &fluids),
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let sides = Vec3::new(1.2, 1.8, 0.8);
    let source = body(1200.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.15, 0.3, 0.1))).unwrap();
    for z in [0.05, -0.4, -1.9] {
        let state = body
            .state_from_kinematics(
                Se3::from_parts(rotation, Vec3::new(0.1, -0.2, z)).unwrap(),
                Twist::zero(),
            )
            .unwrap();
        let load = medium.load(&mesh, &body, state, || false).unwrap();
        for unit in [
            fs_mbd::Vec3::new(1.0, 0.0, 0.0),
            fs_mbd::Vec3::new(0.0, 1.0, 0.0),
            fs_mbd::Vec3::new(0.0, 0.0, 1.0),
        ] {
            for rotational in [false, true] {
                let h = 1e-6;
                let energy = |sign: f64| {
                    let old = state.pose();
                    let position = if rotational {
                        old.position_world()
                    } else {
                        old.position_world().add(unit.scale(sign * h))
                    };
                    let orientation = if rotational {
                        old.orientation().right_exp(unit.scale(sign * h)).unwrap()
                    } else {
                        old.orientation()
                    };
                    let changed = RigidBodyState::new(
                        Pose::new(position, orientation).unwrap(),
                        state.linear_momentum_world(),
                        state.angular_momentum_body(),
                    )
                    .unwrap();
                    medium
                        .load(&mesh, &body, changed, || false)
                        .unwrap()
                        .pressure_potential_j()
                };
                let expected = if rotational {
                    load.wrench().torque_body
                } else {
                    load.wrench().force_world
                }
                .dot(unit);
                let measured = -(energy(1.0) - energy(-1.0)) / (2.0 * h);
                assert!(
                    (measured - expected).abs() < 3e-5,
                    "{z} {rotational}: {measured} != {expected}"
                );
            }
        }
    }
}

#[test]
fn source_refinement_preserves_loading_and_uniform_profiles_preserve_native_bits() {
    let coarse_fluids = [fluid(1000.0), fluid(1600.0), fluid(1700.0)];
    let refined_fluids = [
        fluid(1000.0),
        fluid(1300.0),
        fluid(1600.0),
        fluid(1650.0),
        fluid(1700.0),
    ];
    let coarse = DensityProfileHydrostaticMedium::new(
        &samples(&[0.0, 1.0, 3.0], &coarse_fluids),
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let refined = DensityProfileHydrostaticMedium::new(
        &samples(&[0.0, 0.5, 1.0, 2.0, 3.0], &refined_fluids),
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let uniform_fluids = [fluid(1000.0); 3];
    let uniform = DensityProfileHydrostaticMedium::new(
        &samples(&[0.0, 1.0, 3.0], &uniform_fluids),
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let original = HydrostaticMedium::new(&uniform_fluids[0], gravity(), origin()).unwrap();
    let source = body(1300.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    for z in [2.0, 0.3, -0.7, -2.0, -4.0] {
        let state = body.state_from_kinematics(pose(z), Twist::zero()).unwrap();
        let a = coarse.load(&cube(), &body, state, || false).unwrap();
        let b = refined.load(&cube(), &body, state, || false).unwrap();
        close(a.displaced_mass_kg(), b.displaced_mass_kg());
        close(a.pressure_potential_j(), b.pressure_potential_j());
        if let Some(center) = a.center_of_buoyancy_world() {
            vec_close(center, b.center_of_buoyancy_world().unwrap());
        }
        assert_eq!(
            uniform.load(&cube(), &body, state, || false).unwrap(),
            original.load(&cube(), &body, state, || false).unwrap()
        );
        let controls = HydrostaticStepControls {
            duration_s: 0.001,
            maximum_energy_defect_j: 0.1,
        };
        assert_eq!(
            uniform
                .step(&cube(), &body, state, controls, || false)
                .unwrap(),
            original
                .step(&cube(), &body, state, controls, || false)
                .unwrap()
        );
    }
}

#[test]
fn native_motion_restores_density_equilibrium_without_changing_body_mass() {
    let fluids = [fluid(1000.0), fluid(1500.0), fluid(1800.0)];
    let medium = DensityProfileHydrostaticMedium::new(
        &samples(&[0.0, 2.0, 4.0], &fluids),
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let sides = Vec3::new(0.5, 0.5, 0.5);
    let source = body(1250.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let exact = -1.0 + 0.05 * (2.0_f64.sqrt() * 1.2).cos();
    let mut errors = Vec::new();
    for (dt, steps) in [(0.02, 60), (0.01, 120), (0.005, 240)] {
        let mut state = body
            .state_from_kinematics(pose(-0.95), Twist::zero())
            .unwrap();
        for _ in 0..steps {
            state = medium
                .step(
                    &mesh,
                    &body,
                    state,
                    HydrostaticStepControls {
                        duration_s: dt,
                        maximum_energy_defect_j: 0.001,
                    },
                    || false,
                )
                .unwrap()
                .dynamics
                .state_after;
        }
        errors.push((state.pose().position_world().z - exact).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
    close(body.properties().mass(), 1250.0 * 0.125);
}

#[test]
fn source_snapshots_are_frozen_and_invalid_or_cancelled_profiles_refuse() {
    use std::cell::Cell;
    struct Counted {
        state: FluidTransportState,
        calls: Cell<usize>,
    }
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.calls.set(self.calls.get() + 1);
            Ok(self.state)
        }
    }
    let sources = [
        Counted {
            state: fluid(1000.0),
            calls: Cell::new(0),
        },
        Counted {
            state: fluid(1500.0),
            calls: Cell::new(0),
        },
    ];
    let samples = [
        HydrostaticDensitySample {
            depth_m: 0.0,
            source: &sources[0],
        },
        HydrostaticDensitySample {
            depth_m: 1.0,
            source: &sources[1],
        },
    ];
    let medium =
        DensityProfileHydrostaticMedium::new(&samples, gravity(), origin(), || false).unwrap();
    let source = body(1250.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(pose(-1.0), Twist::zero())
        .unwrap();
    for _ in 0..3 {
        medium.load(&cube(), &body, state, || false).unwrap();
    }
    assert_eq!(sources.each_ref().map(|s| s.calls.get()), [1, 1]);
    assert_eq!(medium.samples().count(), 2);
    assert!(matches!(
        medium.load(&cube(), &body, state, || true),
        Err(HydrostaticError::Cancelled)
    ));
    for depth in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(medium.density_at_depth_kg_m3(depth).is_err());
    }
    assert!(matches!(
        DensityProfileHydrostaticMedium::new(&samples, gravity(), origin(), || true),
        Err(HydrostaticError::Cancelled)
    ));
    for depths in [[0.0, 0.0], [0.0, -1.0], [0.0, f64::NAN], [1.0, 2.0]] {
        let invalid = [
            HydrostaticDensitySample {
                depth_m: depths[0],
                source: samples[0].source,
            },
            HydrostaticDensitySample {
                depth_m: depths[1],
                source: samples[1].source,
            },
        ];
        assert!(
            DensityProfileHydrostaticMedium::new(&invalid, gravity(), origin(), || false).is_err()
        );
    }
    let inverted = [
        HydrostaticDensitySample {
            depth_m: 0.0,
            source: samples[1].source,
        },
        HydrostaticDensitySample {
            depth_m: 1.0,
            source: samples[0].source,
        },
    ];
    assert!(
        DensityProfileHydrostaticMedium::new(&inverted, gravity(), origin(), || false).is_err()
    );
}
