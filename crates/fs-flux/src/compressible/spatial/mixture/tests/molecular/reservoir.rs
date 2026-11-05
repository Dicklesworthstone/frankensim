use super::*;
use crate::compressible::spatial::mixture::molecular::{
    MolecularEuler3d, MolecularStepControl, StabilizedDissipationControl,
};
use crate::compressible::spatial::mixture::{ReservoirBoundary, ReservoirFaceState};

const ALL_OPEN: Boundaries3d = Boundaries3d {
    x: OPEN,
    y: OPEN,
    z: OPEN,
};

fn bind(
    flow: &mut MixtureEuler3d,
    sides: &[BoundarySide],
    mut sample: impl FnMut(BoundarySide, [f64; 3]) -> ([f64; 3], [f64; 3], f64),
) {
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(sides)
        .unwrap();
    flow.prepare_reservoir_profiles(
        sides,
        bytes,
        |side, _, point, densities| {
            let (rho, velocity, pressure_pa) = sample(side, point);
            densities.copy_from_slice(&rho);
            Ok(ReservoirFaceState {
                velocity,
                pressure_pa,
            })
        },
        || false,
    )
    .unwrap();
}

#[test]
fn g1_six_bath_faces_use_half_cell_fourier_stress_and_prescribed_work() {
    let grid = CartesianGrid::new(2, 3, 4, 0.2, 0.3, 0.4).unwrap();
    let rho = [1.0, 0.0, 0.0];
    let initial = [0.07, -0.03, 0.05];
    let bath = [-0.04, 0.08, -0.02];
    let mu = 0.02;
    let bulk = 0.01;
    let conductivity = 0.03;
    let dt = 1e-9;
    for axis in 0..3 {
        let mut axes = [AxisBoundaries::Periodic; 3];
        axes[axis] = OPEN;
        let topology = Boundaries3d {
            x: axes[0],
            y: axes[1],
            z: axes[2],
        };
        let plain = field(equal_mass_gas(), grid, topology, |_| (rho, initial, 1.0));
        let sides = [BoundarySide::ALL[2 * axis], BoundarySide::ALL[2 * axis + 1]];
        let inputs = sides.map(|side| ReservoirBoundary {
            side,
            partial_densities: &rho,
            velocity: bath,
            pressure_pa: 2.0,
        });
        let bytes = MixtureEuler3d::required_workspace_bytes_with_reservoirs(
            plain.gas(),
            &grid,
            inputs.len(),
        )
        .unwrap();
        let mut flow = MixtureEuler3d::new_with_reservoirs(
            plain.gas().clone(),
            grid,
            plain.states().to_vec(),
            topology,
            &inputs,
            bytes,
            || false,
        )
        .unwrap();
        let mut operator = prepare(&flow, law(flow.gas(), 0.001, mu, bulk, conductivity));
        check(
            operator
                .advance(
                    &mut flow,
                    Time::new(dt),
                    DissipationControl::default(),
                    || false,
                )
                .unwrap(),
        );
        for side in sides {
            let gradient = |component: usize| {
                side.sign() * (bath[component] - initial[component]) / (0.5 * grid.spacing()[axis])
            };
            let traction: [f64; 3] = std::array::from_fn(|component| {
                -(if component == axis {
                    bulk + 4.0 * mu / 3.0
                } else {
                    mu
                }) * gradient(component)
            });
            // Equal molecular masses R give T=p/rho. Boundary energy uses
            // the prescribed bath velocity, including both shear components.
            let energy = -conductivity * side.sign() / (0.5 * grid.spacing()[axis])
                + (0..3).map(|d| bath[d] * traction[d]).sum::<f64>();
            for impulse in operator.boundary_impulses(side) {
                let scale = side.sign() * grid.face_area(side) * dt;
                for d in 0..3 {
                    close(impulse.species_masses[d], 0.0, 1e-15);
                    close(impulse.momentum[d] / scale, traction[d], 2e-7);
                }
                close(impulse.energy / scale, energy, 2e-7);
            }
        }
    }
}

#[test]
fn g1_profile_baths_preserve_three_axis_diffusion_and_original_enthalpy_in_both_integrators() {
    let grid = CartesianGrid::new(4, 3, 2, 0.2, 0.3, 0.4).unwrap();
    let slope = [0.1, 0.2, 0.15];
    let fraction = |p: [f64; 3]| 0.2 + (0..3).map(|d| slope[d] * p[d]).sum::<f64>();
    let diffusion = 0.02;
    let dt = 0.01;
    for stable in [false, true] {
        let mut flow = field(equal_mass_gas(), grid, ALL_OPEN, |p| {
            let a = fraction(p);
            ([a, 1.0 - a, 0.0], [0.0; 3], 1.0)
        });
        bind(&mut flow, &BoundarySide::ALL, |_, p| {
            let a = fraction(p);
            ([a, 1.0 - a, 0.0], [0.0; 3], 1.0)
        });
        let original = flow.states().to_vec();
        let source = SpatialTransport::new(
            flow.gas().clone(),
            law(flow.gas(), diffusion, 0.0, 0.0, 0.0),
            || false,
        )
        .unwrap();
        let bytes =
            MolecularTransport3d::required_stabilized_workspace_bytes(&flow, &source).unwrap();
        let mut operator = MolecularTransport3d::new(&flow, source, bytes, || false).unwrap();
        if stable {
            let report = operator
                .advance_stabilized(
                    &mut flow,
                    Time::new(dt),
                    StabilizedDissipationControl::default(),
                    || false,
                )
                .unwrap();
            check(report.transport);
        } else {
            check(
                operator
                    .advance(
                        &mut flow,
                        Time::new(dt),
                        DissipationControl::default(),
                        || false,
                    )
                    .unwrap(),
            );
        }
        for (actual, expected) in flow.states().iter().zip(&original) {
            close(*actual, *expected, 3e-11);
        }
        for side in BoundarySide::ALL {
            // Equal masses, rho=T=p=1 reduce Maxwell--Stefan to j_A=-D grad Y_A.
            // Original Cv_A-Cv_B=-1 and Tmin=0.001 give h_A-h_B=-0.999.
            let j = -diffusion * slope[side.axis()];
            let scale = side.sign() * grid.face_area(side) * dt;
            for impulse in operator.boundary_impulses(side) {
                close(impulse.species_masses[0] / scale, j, 2e-10);
                close(impulse.species_masses[1] / scale, -j, 2e-10);
                assert_eq!(impulse.species_masses[2], 0.0);
                close(impulse.energy / scale, -0.999 * j, 3e-10);
            }
        }
        assert_eq!(flow.time().value(), 0.0);
    }
}

#[test]
fn g1_tangential_profile_shear_uses_one_sided_and_periodic_face_derivatives() {
    let mu = 0.02;
    let dt = 1e-9;
    for periodic in [false, true] {
        let grid = CartesianGrid::new(2, 4, 3, 0.2, 0.3, 0.4).unwrap();
        let topology = Boundaries3d {
            x: OPEN,
            y: if periodic {
                AxisBoundaries::Periodic
            } else {
                WALLS
            },
            z: if periodic {
                AxisBoundaries::Periodic
            } else {
                WALLS
            },
        };
        let ky = 2.0 * std::f64::consts::PI / (4.0 * 0.3);
        let kz = 2.0 * std::f64::consts::PI / (3.0 * 0.4);
        let profile = |p: [f64; 3]| {
            if periodic {
                0.03 * (ky * p[1]).sin() + 0.02 * (kz * p[2]).cos()
            } else {
                0.02 * p[1] * p[1] + 0.03 * p[2]
            }
        };
        let mut flow = field(equal_mass_gas(), grid, topology, |_| {
            ([1.0, 0.0, 0.0], [0.0; 3], 1.0)
        });
        let sides = [BoundarySide::XMin, BoundarySide::XMax];
        bind(&mut flow, &sides, |_, p| {
            ([1.0, 0.0, 0.0], [profile(p), 0.0, 0.0], 1.0)
        });
        let mut operator = prepare(&flow, law(flow.gas(), 0.001, mu, 0.0, 0.0));
        check(
            operator
                .advance(
                    &mut flow,
                    Time::new(dt),
                    DissipationControl::default(),
                    || false,
                )
                .unwrap(),
        );
        for side in sides {
            for (face, impulse) in operator.boundary_impulses(side).enumerate() {
                let p = grid.face_center(side, face).unwrap();
                let [dy, dz] = if periodic {
                    [
                        0.03 * (ky * 0.3).sin() / 0.3 * (ky * p[1]).cos(),
                        -0.02 * (kz * 0.4).sin() / 0.4 * (kz * p[2]).sin(),
                    ]
                } else {
                    [0.04 * p[1], 0.03]
                };
                let scale = side.sign() * grid.face_area(side) * dt;
                // The cell trace is initially uniform: the source midpoint
                // tangential derivative is half the exterior-face derivative.
                close(impulse.momentum[1] / scale, -0.5 * mu * dy, 2e-8);
                close(impulse.momentum[2] / scale, -0.5 * mu * dz, 2e-8);
            }
        }
    }
}

fn coupled_profile() -> MolecularEuler3d<ConstantTransport> {
    let grid = CartesianGrid::new(2, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let mut flow = field(equal_mass_gas(), grid, ALL_OPEN, |_| {
        ([0.7, 0.3, 0.0], [0.0; 3], 1.0)
    });
    bind(&mut flow, &BoundarySide::ALL, |_, _| {
        ([0.7, 0.3, 0.0], [0.0; 3], 1.0)
    });
    let molecular = prepare(&flow, law(flow.gas(), 0.002, 0.003, 0.001, 0.004));
    let bytes = MolecularEuler3d::required_workspace_bytes(&flow, &molecular).unwrap();
    let result = MolecularEuler3d::new(flow, molecular, bytes, || false).unwrap();
    assert_eq!(result.workspace_bytes(), bytes);
    result
}

fn drive_profiles(owner: &mut MolecularEuler3d<ConstantTransport>) {
    owner
        .update_reservoir_profiles(
            |side, face, _, rho| {
                let a = 0.4 + 0.01 * face as f64;
                rho.copy_from_slice(&[a, 1.0 - a, 0.0]);
                Ok(ReservoirFaceState {
                    velocity: [0.03, -0.02, 0.01],
                    pressure_pa: 1.1 + 0.01 * side.axis() as f64,
                })
            },
            || false,
        )
        .unwrap();
}

#[test]
fn g4_coupled_profile_updates_and_late_step_cancellation_preserve_both_owners() {
    let mut reference = coupled_profile();
    let before = reference.flow().states().to_vec();
    let error = reference.update_reservoir_profiles(
        |side, face, _, rho| {
            rho.copy_from_slice(&[0.2, 0.8, 0.0]);
            if side == BoundarySide::ZMax && face == 3 {
                return Err(Error::Invalid("last face refused"));
            }
            Ok(ReservoirFaceState {
                velocity: [0.03; 3],
                pressure_pa: 1.2,
            })
        },
        || false,
    );
    assert!(error.is_err());
    assert_eq!(reference.flow().states(), before);
    for side in BoundarySide::ALL {
        for face in 0..reference.flow().grid().face_count(side) {
            let q = reference.flow().reservoir_face(side, face).unwrap();
            assert_eq!(q.partial_densities, &[0.7, 0.3, 0.0]);
            assert_eq!(q.velocity, [0.0; 3]);
            assert_eq!(q.pressure_pa, 1.0);
        }
    }
    drive_profiles(&mut reference);
    let mut calls = 0;
    let expected = reference
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            calls += 1;
            false
        })
        .unwrap();
    check(expected.first);
    check(expected.second);
    assert_ne!(reference.flow().states(), before);
    let mut cancelled = coupled_profile();
    drive_profiles(&mut cancelled);
    let mut count = 0;
    assert!(matches!(
        cancelled.step(Time::new(0.001), MolecularStepControl::default(), || {
            count += 1;
            count == calls
        }),
        Err(DissipationError::Cancelled)
    ));
    assert_eq!(cancelled.flow().states(), before);
    assert_eq!(cancelled.flow().time().value(), 0.0);
    assert!(cancelled.last_step().is_none());
    assert!(cancelled.molecular().last_step().is_none());
    let actual = cancelled
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(cancelled.flow().states(), reference.flow().states());
    for side in BoundarySide::ALL {
        for face in 0..cancelled.flow().grid().face_count(side) {
            assert_eq!(
                cancelled.flow().reservoir_face(side, face),
                reference.flow().reservoir_face(side, face)
            );
        }
    }
}
