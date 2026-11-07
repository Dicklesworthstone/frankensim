use super::*;
use fs_gas::multicomponent::dissipation::DissipationError;

#[derive(Clone, Copy)]
enum Profile {
    Constant,
    Linear,
    CellPeaks,
    FacePeak,
}

struct SampledViscosity {
    native: ConstantTransport,
    profile: Profile,
}

impl TransportLaw for SampledViscosity {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.native.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        let mut value =
            self.native
                .evaluate(temperature, pressure, fractions, binary, cancelled)?;
        let t = temperature.value();
        let midpoint = (1.25..1.75).contains(&t);
        value.dynamic_viscosity = DynViscosity::new(match self.profile {
            Profile::Constant => 1.0,
            Profile::Linear => 1.0 + 299.0 * (t - 1.0).max(0.0),
            Profile::CellPeaks => {
                if midpoint {
                    1.0
                } else {
                    300.0
                }
            }
            Profile::FacePeak => {
                if midpoint {
                    1e30
                } else {
                    1.0
                }
            }
        });
        Ok(value)
    }
}

fn source(gas: &FrozenGas, profile: Profile, conductivity: f64) -> SampledViscosity {
    SampledViscosity {
        native: law(gas, 1.0, 0.0, conductivity),
        profile,
    }
}

fn operator(
    flow: &Euler3d,
    profile: Profile,
    conductivity: f64,
    modeled: bool,
) -> MolecularTransport3d<SampledViscosity> {
    let native = source(flow.gas(), profile, conductivity);
    #[cfg(feature = "experimental-sgs")]
    if modeled {
        use crate::compressible::spatial::molecular::sgs::{SgsModel, SgsUse, SubgridConfig};
        return MolecularTransport3d::with_sgs(
            flow,
            native,
            SubgridConfig {
                model: SgsModel::Smagorinsky { cs: 0.18 },
                filter_width_factor: 1.0,
                turbulent_prandtl: 0.8,
                maximum_eddy_to_molecular: 20.0,
                usage: SgsUse::Exploration,
            },
            || false,
        )
        .unwrap();
    }
    assert!(!modeled);
    MolecularTransport3d::new(flow, native, || false).unwrap()
}

fn prepare<L: TransportLaw>(op: &mut MolecularTransport3d<L>, flow: &Euler3d) -> f64 {
    spatial(
        &op.gas,
        &op.grid,
        op.boundaries,
        &op.walls,
        None,
        flow.states(),
        &mut op.law,
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.viscosity,
        &mut op.work.flux,
        &mut op.work.binary,
        #[cfg(feature = "experimental-sgs")]
        op.sgs.as_mut(),
        &mut || false,
    )
    .unwrap()
}

fn vortex() -> Euler3d {
    vortex_in(PERIODIC)
}

fn vortex_in(boundaries: Boundaries3d) -> Euler3d {
    let a = [0.2, 0.3, 0.1, 0.3];
    let b = [0.0, 1.0, 0.0, -1.0];
    field(
        gas(0.0),
        CartesianGrid::new(4, 4, 1, 1.0, 1.0, 1.0).unwrap(),
        boundaries,
        |x| {
            let i = x[0] as usize;
            let j = x[1] as usize;
            [
                1.0,
                -a[i] * b[j],
                a[j] * b[i],
                0.0,
                if i == 0 && j == 0 { 2.0 } else { 1.0 },
            ]
        },
    )
}

#[test]
fn g1_variable_molecular_viscosity_dissipates_the_original_planar_and_spatial_vortex() {
    // The old midpoint-coefficient cross stress gives dK/dt=+3.613333... .
    // Original cell sampling and compatible mechanics give -2.366666... .
    // Periodic, no-slip, finite Navier, and natural free-slip domains use
    // the same molecular energy identity, including variable material data.
    for wall_mode in 0..4 {
        let boundaries = if wall_mode != 0 {
            Boundaries3d {
                x: WALL,
                y: WALL,
                z: AxisBoundaries::Periodic,
            }
        } else {
            PERIODIC
        };
        let mut space = vortex_in(boundaries);
        let source_gas = space.gas().clone();
        let cells = space
            .states()
            .iter()
            .map(|q| {
                let p = q.primitive(&source_gas, || false).unwrap();
                planar::State::from_primitive(
                    &source_gas,
                    p.density,
                    [p.velocity[0], p.velocity[1]],
                    p.pressure,
                )
                .unwrap()
            })
            .collect();
        let mut plane = planar::Euler2d::new(
            source_gas,
            planar::CartesianGrid::new(4, 4, 1.0, 1.0, 1.0).unwrap(),
            cells,
            planar::Boundaries2d {
                x: boundaries.x,
                y: boundaries.y,
            },
            || false,
        )
        .unwrap();
        let mut planar_op = planar::molecular::MolecularTransport2d::new(
            &plane,
            source(plane.gas(), Profile::Linear, 0.0),
            || false,
        )
        .unwrap();
        let mut spatial_op = operator(&space, Profile::Linear, 0.0, false);
        if wall_mode != 0 {
            for (plane_side, space_side) in [
                (planar::BoundarySide::XMin, BoundarySide::XMin),
                (planar::BoundarySide::XMax, BoundarySide::XMax),
                (planar::BoundarySide::YMin, BoundarySide::YMin),
                (planar::BoundarySide::YMax, BoundarySide::YMax),
            ] {
                let planar_wall = match wall_mode {
                    1 => planar::molecular::MolecularWall2d::NoSlip {
                        velocity: [0.0; 2],
                        thermal: WallThermalCondition::Adiabatic,
                    },
                    2 => planar::molecular::MolecularWall2d::NavierSlip {
                        velocity: [0.0; 2],
                        thermal: WallThermalCondition::Adiabatic,
                        slip_length_m: 0.3,
                        fluid_heat_fraction: 1.0,
                    },
                    _ => planar::molecular::MolecularWall2d::SlipAdiabatic,
                };
                let spatial_wall = match wall_mode {
                    1 => MolecularWall3d::NoSlip {
                        velocity: [0.0; 3],
                        thermal: WallThermalCondition::Adiabatic,
                    },
                    2 => MolecularWall3d::NavierSlip {
                        velocity: [0.0; 3],
                        thermal: WallThermalCondition::Adiabatic,
                        slip_length_m: 0.3,
                        fluid_heat_fraction: 1.0,
                    },
                    _ => MolecularWall3d::SlipAdiabatic,
                };
                planar_op = planar_op.with_wall(plane_side, planar_wall).unwrap();
                spatial_op = spatial_op.with_wall(space_side, spatial_wall).unwrap();
            }
        }
        assert!(planar_op.uses_compatible_molecular_stress());
        assert!(spatial_op.uses_compatible_molecular_stress());
        let before = partition(&space);
        let duration = Time::new(1e-6);
        let planar_report = planar_op
            .advance(&mut plane, duration, MolecularControl::default(), || false)
            .unwrap();
        let spatial_report = spatial_op
            .advance(&mut space, duration, MolecularControl::default(), || false)
            .unwrap();
        assert_eq!(planar_report.accepted_steps, 1);
        assert_eq!(spatial_report.accepted_steps, 1);
        let after = partition(&space);
        assert!(after.0 < before.0);
        assert!(after.1 > before.1);
        close(after.0 + after.1, before.0 + before.1, 3e-12);
        close(spatial_report.conservation_residual.energy, 0.0, 3e-12);
        close(planar_report.conservation_residual.energy, 0.0, 3e-12);
        for (two, three) in plane.states().iter().zip(space.states()) {
            close(two.momentum[0], three.momentum[0], 3e-12);
            close(two.momentum[1], three.momentum[1], 3e-12);
            close(two.energy, three.energy, 3e-12);
            assert_eq!(three.momentum[2], 0.0);
        }
    }
}

fn alternating() -> Euler3d {
    field(
        gas(0.0),
        CartesianGrid::new(2, 2, 2, 1.0, 1.0, 1.0).unwrap(),
        PERIODIC,
        |x| {
            let high = (x[0] as usize + x[1] as usize + x[2] as usize) % 2 == 0;
            let sign = if high { 1.0 } else { -1.0 };
            [
                1.0,
                0.1 * sign,
                -0.2 * sign,
                0.05 * sign,
                if high { 2.0 } else { 1.0 },
            ]
        },
    )
}

#[test]
fn g0_nonlinear_material_bounds_use_cells_and_extreme_face_peaks_keep_stress_heat_and_sgs() {
    let flow = alternating();
    let mut cells = operator(&flow, Profile::CellPeaks, 0.2, false);
    let rate = prepare(&mut cells, &flow);
    // All old midpoint samples are one, but every actual cell sample is300.
    close(rate, 4.0 * (20.0 / 3.0) * 300.0 * 3.0, 2e-10);
    for modeled in [false, true] {
        #[cfg(not(feature = "experimental-sgs"))]
        if modeled {
            continue;
        }
        let mut reference = operator(&flow, Profile::Constant, 0.2, modeled);
        let mut peaked = operator(&flow, Profile::FacePeak, 0.2, modeled);
        let reference_rate = prepare(&mut reference, &flow);
        let actual_rate = prepare(&mut peaked, &flow);
        close(actual_rate, reference_rate, 2e-12);
        for axis in 0..3 {
            for (actual, expected) in peaked.work.flux[axis]
                .iter()
                .zip(&reference.work.flux[axis])
            {
                close(actual.normal_momentum, expected.normal_momentum, 2e-13);
                for tangent in 0..2 {
                    close(
                        actual.transverse_momentum[tangent],
                        expected.transverse_momentum[tangent],
                        2e-13,
                    );
                }
                close(actual.energy, expected.energy, 2e-13);
            }
        }
    }
}

#[test]
fn g0_compatible_scope_admits_all_reflecting_walls_and_preserves_bath_and_sgs_fallback() {
    let flow = field(
        gas(0.0),
        CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap(),
        WALLS,
        |_| [1.0, 0.0, 0.0, 0.0, 1.0],
    );
    let free = operator(&flow, Profile::Linear, 0.0, false);
    assert!(free.uses_compatible_molecular_stress());
    let prescribed = free
        .with_wall(
            BoundarySide::XMin,
            MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Adiabatic,
            },
        )
        .unwrap();
    assert!(prescribed.uses_compatible_molecular_stress());
    for slip_length_m in [0.0, 0.1] {
        let navier = operator(&flow, Profile::Linear, 0.0, false)
            .with_wall(
                BoundarySide::XMin,
                MolecularWall3d::NavierSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m,
                    fluid_heat_fraction: 0.0,
                },
            )
            .unwrap();
        assert!(navier.uses_compatible_molecular_stress());
    }
    #[cfg(feature = "experimental-sgs")]
    {
        let modeled = operator(&flow, Profile::Linear, 0.0, true)
            .with_wall(
                BoundarySide::XMin,
                MolecularWall3d::NoSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                },
            )
            .unwrap();
        assert!(modeled.uses_compatible_molecular_stress());
        let modeled = modeled
            .with_equilibrium_sgs_wall(
                BoundarySide::XMin,
                crate::compressible::spatial::molecular::sgs::EquilibriumWallConfig {
                    maximum_wall_distance_plus: 1000.0,
                    maximum_mach_number: 0.3,
                },
            )
            .unwrap();
        assert!(!modeled.uses_compatible_molecular_stress());
    }
    assert!(
        !operator(&flow, Profile::Linear, 0.0, false)
            .with_reservoir_exchange()
            .uses_compatible_molecular_stress()
    );

    let mut compatible_flow = vortex();
    let mut native_flow = vortex();
    let mut compatible = operator(&compatible_flow, Profile::Constant, 0.02, false);
    // Empty bath exchange on a periodic grid selects the retained native
    // mechanical path without changing a physical face or its source data.
    let mut native =
        operator(&native_flow, Profile::Constant, 0.02, false).with_reservoir_exchange();
    let expected = native
        .advance(
            &mut native_flow,
            Time::new(1e-5),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    let actual = compatible
        .advance(
            &mut compatible_flow,
            Time::new(1e-5),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(compatible_flow.states(), native_flow.states());
}

#[test]
fn g0_prescribed_walls_use_cell_stress_and_keep_original_heat_at_extreme_wall_samples() {
    let grid = CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap();
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        let mut boundaries = PERIODIC;
        match axis {
            0 => boundaries.x = WALL,
            1 => boundaries.y = WALL,
            2 => boundaries.z = WALL,
            _ => unreachable!(),
        }
        let flow = field(gas(0.0), grid, boundaries, |x| {
            let mut velocity = [0.0; 3];
            velocity[axis] = 0.02 * x[a] + 0.01 * x[b];
            velocity[a] = 0.1 + 0.03 * x[a] - 0.02 * x[axis];
            velocity[b] = -0.1 + 0.04 * x[b] + 0.01 * x[axis];
            [1.0, velocity[0], velocity[1], velocity[2], 1.0]
        });
        let mut wall_velocity = [0.0; 3];
        wall_velocity[a] = 0.4;
        wall_velocity[b] = -0.3;
        for slip_length in [0.0, 0.7] {
            let thermal = WallThermalCondition::Isothermal(Temperature::new(2.0));
            let heat_fraction = if slip_length == 0.0 { 1.0 } else { 0.35 };
            let wall = if slip_length == 0.0 {
                MolecularWall3d::NoSlip {
                    velocity: wall_velocity,
                    thermal,
                }
            } else {
                MolecularWall3d::NavierSlip {
                    velocity: wall_velocity,
                    thermal,
                    slip_length_m: slip_length,
                    fluid_heat_fraction: heat_fraction,
                }
            };
            let mut reference = operator(&flow, Profile::Constant, 0.2, false)
                .with_wall(side, wall)
                .unwrap();
            let mut peaked = operator(&flow, Profile::FacePeak, 0.2, false)
                .with_wall(side, wall)
                .unwrap();
            // Cell mu=1, but the original thermal wall midpoint has mu=1e30.
            // Traction, physical slip heat, and rate use the cell material.
            let rate = prepare(&mut reference, &flow);
            close(prepare(&mut peaked, &flow), rate, 2e-12);
            for d in 0..3 {
                for (actual, expected) in peaked.work.flux[d].iter().zip(&reference.work.flux[d]) {
                    close(actual.normal_momentum, expected.normal_momentum, 2e-13);
                    for tangent in 0..2 {
                        close(
                            actual.transverse_momentum[tangent],
                            expected.transverse_momentum[tangent],
                            2e-13,
                        );
                    }
                    close(actual.energy, expected.energy, 2e-13);
                }
            }
            let mut face = [1; 3];
            face[axis] = if side.sign() > 0.0 { 3 } else { 0 };
            let flux = peaked.work.flux[axis][grid.cell_face(axis, face)];
            face[axis] = if side.sign() > 0.0 { 2 } else { 0 };
            let cell = grid.index(face);
            let p = peaked.work.primitive[cell];
            let g = peaked.work.gradient[cell];
            let h = 0.5 * grid.spacing[axis];
            let normal = -(4.0 / 3.0) * side.sign() * p.velocity[axis] / h
                - (2.0 / 3.0) * (g[a][a] + g[b][b]);
            let shear = [a, b].map(|t| {
                (side.sign() * (wall_velocity[t] - p.velocity[t]) + h * g[axis][t])
                    / (h + slip_length)
            });
            let relative_trace = [a, b].map(|t| {
                slip_length / (h + slip_length)
                    * (p.velocity[t] - wall_velocity[t] - side.sign() * h * g[axis][t])
            });
            let dissipation = -side.sign()
                * (shear[0] * relative_trace[0] + shear[1] * relative_trace[1]);
            assert!(dissipation >= 0.0);
            let heat = -0.2 * side.sign() * (2.0 - p.temperature) / h;
            close(flux.normal_momentum, -normal, 2e-13);
            for t in 0..2 {
                close(flux.transverse_momentum[t], -shear[t], 2e-13);
            }
            close(
                flux.energy,
                -shear[0] * wall_velocity[a] - shear[1] * wall_velocity[b]
                    + side.sign() * (1.0 - heat_fraction) * dissipation
                    + heat,
                2e-13,
            );
        }
    }
}

#[test]
fn g1_large_navier_slip_converges_to_the_same_natural_free_slip_gradient() {
    let flow = field(
        gas(0.0),
        CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap(),
        WALLS,
        |x| [
            1.0,
            0.03 * x[1] - 0.02 * x[2],
            0.04 * x[0] + 0.01 * x[2],
            -0.02 * x[0] + 0.03 * x[1],
            1.0 + 0.2 * x[0] + 0.1 * x[1],
        ],
    );
    let mut free = operator(&flow, Profile::Linear, 0.0, false);
    let mut large = operator(&flow, Profile::Linear, 0.0, false);
    for side in BoundarySide::ALL {
        large = large
            .with_wall(
                side,
                MolecularWall3d::NavierSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m: 1e12,
                    fluid_heat_fraction: 1.0,
                },
            )
            .unwrap();
    }
    assert!(free.uses_compatible_molecular_stress());
    assert!(large.uses_compatible_molecular_stress());
    let rate = prepare(&mut free, &flow);
    close(prepare(&mut large, &flow), rate, 2e-12);
    // The natural wall trace includes the normal-velocity cross derivative.
    // Matching only face traction while retaining an even ghost would fail
    // this comparison of the actual interior stress gradients.
    for (actual, expected) in large.work.gradient.iter().zip(&free.work.gradient) {
        for component in 0..3 {
            for derivative in 0..3 {
                close(
                    actual[component][derivative],
                    expected[component][derivative],
                    2e-11,
                );
            }
        }
    }
    for axis in 0..3 {
        for (actual, expected) in large.work.flux[axis].iter().zip(&free.work.flux[axis]) {
            close(actual.normal_momentum, expected.normal_momentum, 3e-9);
            for tangent in 0..2 {
                close(
                    actual.transverse_momentum[tangent],
                    expected.transverse_momentum[tangent],
                    3e-9,
                );
            }
            close(actual.energy, expected.energy, 3e-9);
        }
    }
}

fn moving_channel() -> Euler3d {
    field(
        gas(0.0),
        CartesianGrid::new(4, 1, 1, 0.25, 0.3, 0.4).unwrap(),
        Boundaries3d {
            x: WALL,
            y: AxisBoundaries::Periodic,
            z: AxisBoundaries::Periodic,
        },
        |x| [1.0, 0.0, 0.0, 0.0, if x[0] < 0.25 { 2.0 } else { 1.0 }],
    )
}

fn moving_operator(flow: &Euler3d, slip_length: f64) -> MolecularTransport3d<SampledViscosity> {
    let velocity = [0.0, 0.1, 0.0];
    let thermal = WallThermalCondition::Adiabatic;
    let wall = if slip_length == 0.0 {
        MolecularWall3d::NoSlip { velocity, thermal }
    } else {
        MolecularWall3d::NavierSlip {
            velocity,
            thermal,
            slip_length_m: slip_length,
            fluid_heat_fraction: 1.0,
        }
    };
    operator(flow, Profile::Linear, 0.01, false)
        .with_wall(BoundarySide::XMin, wall)
        .unwrap()
}

#[test]
fn g1_moving_wall_curl_uses_the_same_cell_stress_as_conservative_wall_power() {
    for slip_length in [0.0, 0.125] {
        let mut flow = moving_channel();
        let mut op = moving_operator(&flow, slip_length);
        assert!(op.uses_compatible_molecular_stress());
        let initial = flow.states().to_vec();
        let before = flow.totals();
        let time = flow.time();
        let mut work = op
            .vorticity_workspace(op.required_vorticity_workspace_bytes().unwrap())
            .unwrap();
        let observed = op.vorticity(&flow, &mut work, || false).unwrap();
        // With h=cell spacing, the first cell receives
        // mu*Uwall/[h*(h/2+b)] acceleration. Its neighbor's centered curl
        // is -mu*Uwall/[2*h^2*(h/2+b)], including the no-slip limit.
        close(
            observed[1].stress[2],
            -300.0 * 0.1 / (2.0 * 0.25_f64.powi(2) * (0.125 + slip_length)),
            3e-10,
        );
        close(observed[1].stress[0], 0.0, 2e-13);
        close(observed[1].stress[1], 0.0, 2e-13);
        assert_eq!(flow.states(), initial);
        assert_eq!(flow.time(), time);
        let report = op
            .advance(
                &mut flow,
                Time::new(1e-7),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        assert!(flow.totals().energy > before.energy);
        assert!(flow.totals().momentum[1] > before.momentum[1]);
        let mut outward_energy = 0.0;
        for impulse in op.boundary_impulses(BoundarySide::XMin) {
            // All interface heat returns to the fluid in this case.
            close(impulse.energy, 0.1 * impulse.momentum[1], 2e-13);
            outward_energy += impulse.energy;
        }
        close(flow.totals().energy - before.energy, -outward_energy, 3e-13);
        close(report.conservation_residual.energy, 0.0, 3e-13);
    }
}

#[test]
fn g4_variable_material_late_cancellation_preserves_flow_and_replays_native_sources() {
    for moving in [false, true] {
        let make_flow = || if moving { moving_channel() } else { vortex() };
        let make_op = |flow: &Euler3d| {
            if moving {
                moving_operator(flow, 0.125)
            } else {
                operator(flow, Profile::Linear, 0.01, false)
            }
        };
        let mut reference = make_flow();
        let mut reference_op = make_op(&reference);
        let mut calls = 0;
        let expected = reference_op
            .advance(
                &mut reference,
                Time::new(1e-6),
                MolecularControl::default(),
                || {
                    calls += 1;
                    false
                },
            )
            .unwrap();
        let mut flow = make_flow();
        let original = flow.states().to_vec();
        let time = flow.time();
        let mut op = make_op(&flow);
        let mut seen = 0;
        assert_eq!(
            op.advance(
                &mut flow,
                Time::new(1e-6),
                MolecularControl::default(),
                || {
                    seen += 1;
                    seen == calls
                }
            )
            .unwrap_err(),
            MolecularError::Cancelled
        );
        assert_eq!(flow.states(), original);
        assert_eq!(flow.time(), time);
        for side in BoundarySide::ALL {
            assert!(
                op.boundary_impulses(side)
                    .iter()
                    .all(|q| *q == Totals::default())
            );
        }
        let actual = op
            .advance(
                &mut flow,
                Time::new(1e-6),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(flow.states(), reference.states());
    }
}
