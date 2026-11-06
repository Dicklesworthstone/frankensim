use super::*;
use crate::compressible::mixture::molecular::{
    MolecularWall2d, StabilizedDissipationControl, TransportLaw, WallThermalCondition,
};

fn channel(axis: usize) -> Boundaries2d {
    if axis == 0 {
        Boundaries2d {
            x: WALLS,
            y: AxisBoundaries::Periodic,
        }
    } else {
        Boundaries2d {
            x: AxisBoundaries::Periodic,
            y: WALLS,
        }
    }
}

struct LocalTransport(ConstantTransport);

impl TransportLaw for LocalTransport {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.0.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        let mut value =
            self.0
                .evaluate(temperature, pressure, mole_fractions, binary, cancelled)?;
        value.dynamic_viscosity =
            DynViscosity::new(0.02 * temperature.value() * (1.0 + mole_fractions[0]));
        value.thermal_conductivity = ThermalConductivityQuantityV1::new(
            0.03 * temperature.value() * (1.0 + 0.5 * mole_fractions[0]),
        );
        Ok(value)
    }
}

#[test]
fn g1_planar_navier_walls_use_local_mole_composition_temperature_and_four_orientations() {
    let grid = CartesianGrid::new(3, 3, 0.2, 0.3, 0.4).unwrap();
    let (slip, fraction, duration, wall_temperature) = (0.17, 0.35, 1e-8, 2.5);
    for side in SIDES {
        let (axis, tangent) = (side.axis(), 1 - side.axis());
        for (first, cell_temperature) in [(0.2, 1.3), (0.7, 2.1)] {
            let mut velocity = [0.0; 2];
            velocity[tangent] = 0.15;
            let mut wall_velocity = [0.0; 2];
            wall_velocity[tangent] = 0.4;
            // A/B molar masses are R/2R: this exercises real mass-to-mole
            // conversion in the original source, not mass fractions as x_i.
            let concentration = first + 0.5 * (1.0 - first);
            let mut flow = field(standard(0.0), grid, channel(axis), |_| {
                (
                    [first, 1.0 - first, 0.0],
                    velocity,
                    concentration * cell_temperature,
                )
            });
            let source = PlanarTransport::new(
                flow.gas().clone(),
                LocalTransport(law(flow.gas(), 0.01, 0.02, 0.003, 0.03)),
                || false,
            )
            .unwrap();
            let bytes = MolecularTransport2d::required_workspace_bytes(&flow, &source).unwrap();
            let wall = MolecularWall2d::NavierSlip {
                velocity: wall_velocity,
                thermal: WallThermalCondition::Isothermal(Temperature::new(wall_temperature)),
                slip_length_m: slip,
                fluid_heat_fraction: fraction,
            };
            let mut op = MolecularTransport2d::new(&flow, source, bytes, || false)
                .unwrap()
                .with_wall(side, wall)
                .unwrap();
            assert_eq!(op.wall(side), wall);
            check(
                op.advance(
                    &mut flow,
                    Time::new(duration),
                    DissipationControl::default(),
                    || false,
                )
                .unwrap(),
            );
            let h = 0.5 * grid.spacing[axis];
            let x_a = first / concentration;
            let face_t = 0.5 * (cell_temperature + wall_temperature);
            let mu = 0.02 * face_t * (1.0 + x_a);
            let conductivity = 0.03 * face_t * (1.0 + 0.5 * x_a);
            let traction = mu * (wall_velocity[tangent] - velocity[tangent]) / (h + slip);
            let loss =
                mu * slip * ((wall_velocity[tangent] - velocity[tangent]) / (h + slip)).powi(2);
            let outward_power = -traction * wall_velocity[tangent]
                + (1.0 - fraction) * loss
                + conductivity * (cell_temperature - wall_temperature) / h;
            let area_time = grid.spacing[tangent] * grid.depth * duration;
            for impulse in op.boundary_impulses(side) {
                assert_eq!(impulse.species_masses, [0.0; 3]);
                close(impulse.momentum[tangent] / area_time, -traction, 3e-8);
                close(impulse.energy / area_time, outward_power, 3e-8);
            }
            assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
        }
    }
}

#[test]
fn g1_planar_navier_tangential_trace_enters_normal_stress() {
    let grid = CartesianGrid::new(3, 3, 0.25, 0.25, 0.6).unwrap();
    let (mu, bulk, slip, slope, duration) = (0.02, 0.003, 0.15, 0.05, 1e-8);
    for side in SIDES {
        let (axis, tangent) = (side.axis(), 1 - side.axis());
        let mut flow = field(equal_mass_gas(), grid, channel(axis), |x| {
            let distance = if side.sign() > 0.0 {
                0.75 - x[axis]
            } else {
                x[axis]
            };
            let mut velocity = [0.0; 2];
            velocity[tangent] = slope * x[tangent] * (distance + slip);
            ([0.4, 0.6, 0.0], velocity, 2.0)
        });
        let mut op = prepare(&flow, law(flow.gas(), 0.01, mu, bulk, 0.0))
            .with_wall(
                side,
                MolecularWall2d::NavierSlip {
                    velocity: [0.0; 2],
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m: slip,
                    fluid_heat_fraction: 1.0,
                },
            )
            .unwrap();
        check(
            op.advance(
                &mut flow,
                Time::new(duration),
                DissipationControl::default(),
                || false,
            )
            .unwrap(),
        );
        // The central face is away from periodic tangential end differences.
        // Its fluid trace is slope*b*x_t, so div_t u_f = slope*b exactly.
        let impulse = op.boundary_impulses(side).nth(1).unwrap();
        close(
            impulse.momentum[axis] / (duration * grid.spacing[tangent] * grid.depth),
            -side.sign() * (bulk - 2.0 * mu / 3.0) * slope * slip,
            3e-9,
        );
        assert_eq!(impulse.species_masses, [0.0; 3]);
    }
}

#[test]
fn g1_planar_slip_couette_has_exact_shear_and_partitioned_heat_in_ssp_and_rkl() {
    let grid = CartesianGrid::new(3, 6, 0.25, 0.2, 0.6).unwrap();
    let (mu, lower_slip, upper_slip, speed, duration) = (0.02, 0.1, 0.2, 0.3, 0.1);
    let slope = speed / (1.2 + lower_slip + upper_slip);
    for stabilized in [false, true] {
        for fraction in [0.0, 1.0] {
            // Equal molar masses eliminate pressure-driven separation when
            // interfacial heat makes the boundary cells warmer. The exact
            // independent shear/heat solution then has no species flux.
            let mut flow = field(equal_mass_gas(), grid, channel(1), |x| {
                ([0.3, 0.7, 0.0], [slope * (x[1] + lower_slip), 0.0], 2.0)
            });
            let initial = flow.states().to_vec();
            let source = PlanarTransport::new(
                flow.gas().clone(),
                law(flow.gas(), 0.01, mu, 0.0, 0.0),
                || false,
            )
            .unwrap();
            let bytes =
                MolecularTransport2d::required_stabilized_workspace_bytes(&flow, &source).unwrap();
            let mut op = MolecularTransport2d::new(&flow, source, bytes, || false).unwrap();
            for (side, velocity, slip) in [
                (BoundarySide::YMin, 0.0, lower_slip),
                (BoundarySide::YMax, speed, upper_slip),
            ] {
                op = op
                    .with_wall(
                        side,
                        MolecularWall2d::NavierSlip {
                            velocity: [velocity, 0.0],
                            thermal: WallThermalCondition::Adiabatic,
                            slip_length_m: slip,
                            fluid_heat_fraction: fraction,
                        },
                    )
                    .unwrap();
            }
            let report = if stabilized {
                let accepted = op
                    .advance_stabilized(
                        &mut flow,
                        Time::new(duration),
                        StabilizedDissipationControl::default(),
                        || false,
                    )
                    .unwrap();
                assert!(accepted.maximum_stages_used >= 3);
                accepted.transport
            } else {
                op.advance(
                    &mut flow,
                    Time::new(duration),
                    DissipationControl::default(),
                    || false,
                )
                .unwrap()
            };
            for (i, (state, before)) in flow
                .states()
                .chunks_exact(6)
                .zip(initial.chunks_exact(6))
                .enumerate()
            {
                for component in 0..5 {
                    close(state[component], before[component], 3e-12);
                }
                let wall_slip = if i / 3 == 0 {
                    lower_slip
                } else if i / 3 == 5 {
                    upper_slip
                } else {
                    0.0
                };
                let heat =
                    duration * mu * slope * slope * (1.0 + fraction * wall_slip / grid.spacing[1]);
                close(state[5] - before[5], heat, 5e-12);
            }
            close(report.absolute_energy_residual_j, 0.0, 3e-10);
            assert!(op.workspace_bytes() + flow.workspace_bytes() <= bytes);
        }
    }
}

#[test]
fn g1_planar_poiseuille_wall_traction_and_discrete_viscous_curvature_agree() {
    let (height, slip, mu, forcing) = (1.0, 0.15, 0.03, 0.02);
    for n in [8, 16] {
        let h = height / n as f64;
        let grid = CartesianGrid::new(3, n, 0.25, h, 0.6).unwrap();
        // The finite-volume steady profile differs from continuum Poiseuille
        // by exactly G*h^2/(8*mu), which converges away with grid refinement.
        let mut flow = field(standard(0.0), grid, channel(1), |x| {
            let u = forcing / (2.0 * mu) * (x[1] * (height - x[1]) + slip * height + 0.25 * h * h);
            ([0.3, 0.7, 0.0], [u, 0.0], 2.0)
        });
        let before = flow.states().to_vec();
        let mut op = prepare(&flow, law(flow.gas(), 0.01, mu, 0.0, 0.0));
        for side in [BoundarySide::YMin, BoundarySide::YMax] {
            op = op
                .with_wall(
                    side,
                    MolecularWall2d::NavierSlip {
                        velocity: [0.0; 2],
                        thermal: WallThermalCondition::Adiabatic,
                        slip_length_m: slip,
                        fluid_heat_fraction: 0.0,
                    },
                )
                .unwrap();
        }
        let dt = 1e-8;
        check(
            op.advance(
                &mut flow,
                Time::new(dt),
                DissipationControl::default(),
                || false,
            )
            .unwrap(),
        );
        for (state, old) in flow.states().chunks_exact(6).zip(before.chunks_exact(6)) {
            close((state[3] - old[3]) / dt, -forcing, 2e-8);
        }
        for side in [BoundarySide::YMin, BoundarySide::YMax] {
            for impulse in op.boundary_impulses(side) {
                close(
                    impulse.momentum[0] / (dt * grid.spacing[0] * grid.depth),
                    0.5 * forcing * height,
                    3e-10,
                );
                assert_eq!(impulse.species_masses, [0.0; 3]);
            }
        }
    }
}

#[test]
fn g0_planar_zero_slip_matches_no_slip_and_wall_admission_is_explicit() {
    let grid = CartesianGrid::new(3, 4, 0.25, 0.25, 0.6).unwrap();
    let create = || {
        field(standard(0.0), grid, channel(1), |x| {
            let a = 0.3 + 0.1 * x[1];
            ([a, 1.0 - a, 0.0], [0.1 * x[1], 0.0], 2.0 + 0.1 * x[1])
        })
    };
    let mut results = Vec::new();
    for navier in [false, true] {
        let mut flow = create();
        let mut op = prepare(&flow, law(flow.gas(), 0.01, 0.02, 0.003, 0.05));
        for (side, speed) in [(BoundarySide::YMin, -0.1), (BoundarySide::YMax, 0.3)] {
            let velocity = [speed, 0.0];
            let thermal = WallThermalCondition::Isothermal(Temperature::new(2.5));
            let wall = if navier {
                MolecularWall2d::NavierSlip {
                    velocity,
                    thermal,
                    slip_length_m: 0.0,
                    fluid_heat_fraction: 0.0,
                }
            } else {
                MolecularWall2d::NoSlip { velocity, thermal }
            };
            op = op.with_wall(side, wall).unwrap();
        }
        let report = op
            .advance(
                &mut flow,
                Time::new(0.002),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        results.push((
            flow.states()
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>(),
            report,
            impulses(&op),
        ));
    }
    assert_eq!(results[0], results[1]);
    for (side, velocity, slip, fraction) in [
        (BoundarySide::XMin, [0.0, 0.0], 0.1, 0.5),
        (BoundarySide::YMin, [0.0, 0.1], 0.1, 0.5),
        (BoundarySide::YMin, [0.0, 0.0], -0.1, 0.5),
        (BoundarySide::YMin, [0.0, 0.0], 0.1, 1.1),
    ] {
        let flow = create();
        let op = prepare(&flow, law(flow.gas(), 0.01, 0.02, 0.0, 0.0));
        assert!(matches!(
            op.with_wall(
                side,
                MolecularWall2d::NavierSlip {
                    velocity,
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m: slip,
                    fluid_heat_fraction: fraction,
                }
            ),
            Err(DissipationError::Invalid(_))
        ));
    }
}
