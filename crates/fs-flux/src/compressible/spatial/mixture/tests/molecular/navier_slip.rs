use super::*;
use crate::compressible::spatial::mixture::molecular::TransportLaw;
use fs_qty::Pressure;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;

fn channel(axis: usize) -> Boundaries3d {
    let mut boundaries = PERIODIC;
    match axis {
        0 => boundaries.x = WALLS,
        1 => boundaries.y = WALLS,
        _ => boundaries.z = WALLS,
    }
    boundaries
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
fn g1_navier_walls_use_actual_temperature_and_mole_composition_without_species_leakage() {
    let grid = CartesianGrid::new(3, 3, 3, 0.2, 0.3, 0.4).unwrap();
    let (slip, fraction, duration, wall_temperature) = (0.17, 0.35, 1e-8, 2.5);
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = [(axis + 1) % 3, (axis + 2) % 3];
        for (first, cell_temperature) in [(0.2, 1.3), (0.7, 2.1)] {
            let mut cell_velocity = [0.0; 3];
            cell_velocity[a] = 0.15;
            cell_velocity[b] = -0.1;
            let mut wall_velocity = [0.0; 3];
            wall_velocity[a] = 0.4;
            wall_velocity[b] = -0.3;
            // The original A and B molar masses are R and 2R. Their mass
            // fractions therefore differ from the mole fractions used by the law.
            let concentration = first + 0.5 * (1.0 - first);
            let mut flow = field(standard(0.0), grid, channel(axis), |_| {
                (
                    [first, 1.0 - first, 0.0],
                    cell_velocity,
                    concentration * cell_temperature,
                )
            });
            let law = LocalTransport(law(flow.gas(), 0.01, 0.02, 0.003, 0.03));
            let source = SpatialTransport::new(flow.gas().clone(), law, || false).unwrap();
            let bytes = MolecularTransport3d::required_workspace_bytes(&flow, &source).unwrap();
            let mut op = MolecularTransport3d::new(&flow, source, bytes, || false)
                .unwrap()
                .with_wall(
                    side,
                    MolecularWall3d::NavierSlip {
                        velocity: wall_velocity,
                        thermal: WallThermalCondition::Isothermal(Temperature::new(
                            wall_temperature,
                        )),
                        slip_length_m: slip,
                        fluid_heat_fraction: fraction,
                    },
                )
                .unwrap();
            let report = op
                .advance(
                    &mut flow,
                    Time::new(duration),
                    DissipationControl::default(),
                    || false,
                )
                .unwrap();
            check(report);
            let h = 0.5 * grid.spacing()[axis];
            let x_a = first / concentration;
            let face_temperature = 0.5 * (cell_temperature + wall_temperature);
            let mu = 0.02 * face_temperature * (1.0 + x_a);
            let conductivity = 0.03 * face_temperature * (1.0 + 0.5 * x_a);
            let traction = [a, b].map(|t| mu * (wall_velocity[t] - cell_velocity[t]) / (h + slip));
            let dissipation: f64 = [a, b]
                .into_iter()
                .map(|t| mu * slip * ((wall_velocity[t] - cell_velocity[t]) / (h + slip)).powi(2))
                .sum();
            let outward_power = -traction[0] * wall_velocity[a] - traction[1] * wall_velocity[b]
                + (1.0 - fraction) * dissipation
                + conductivity * (cell_temperature - wall_temperature) / h;
            let area_time = grid.face_area(side) * duration;
            for impulse in op.boundary_impulses(side) {
                assert_eq!(impulse.species_masses, [0.0; 3]);
                close(impulse.momentum[a] / area_time, -traction[0], 3e-8);
                close(impulse.momentum[b] / area_time, -traction[1], 3e-8);
                close(impulse.energy / area_time, outward_power, 3e-8);
            }
            assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
        }
    }
}

#[test]
fn g1_navier_trace_tangential_derivatives_enter_full_normal_stress() {
    let grid = CartesianGrid::new(3, 3, 3, 0.25, 0.25, 0.25).unwrap();
    let (mu, bulk, slip, slope, duration) = (0.02, 0.003, 0.15, 0.05, 1e-8);
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let tangent = (axis + 1) % 3;
        let mut flow = field(equal_mass_gas(), grid, channel(axis), |x| {
            let distance = if side.sign() > 0.0 {
                0.75 - x[axis]
            } else {
                x[axis]
            };
            let mut velocity = [0.0; 3];
            velocity[tangent] = slope * x[tangent] * (distance + slip);
            ([0.4, 0.6, 0.0], velocity, 2.0)
        });
        let transport = law(flow.gas(), 0.01, mu, bulk, 0.0);
        let mut op = prepare(&flow, transport)
            .with_wall(
                side,
                MolecularWall3d::NavierSlip {
                    velocity: [0.0; 3],
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
        // The center face avoids periodic end differences: u_t at the wall
        // equals slope*b*x_t, giving the exact tangential derivative slope*b.
        let impulse = op.boundary_impulses(side).nth(4).unwrap();
        close(
            impulse.momentum[axis] / (duration * grid.face_area(side)),
            -side.sign() * (bulk - 2.0 * mu / 3.0) * slope * slip,
            3e-9,
        );
        assert_eq!(impulse.species_masses, [0.0; 3]);
    }
}

#[test]
fn g1_mixture_navier_couette_preserves_species_and_turns_shear_work_into_heat() {
    let grid = CartesianGrid::new(3, 3, 6, 0.25, 0.25, 0.2).unwrap();
    let (mu, lower_slip, upper_slip, speed, duration) = (0.02, 0.1, 0.2, 0.3, 0.01);
    let slope = speed / (1.2 + lower_slip + upper_slip);
    let mut flow = field(standard(0.0), grid, channel(2), |x| {
        (
            [0.3, 0.7, 0.0],
            [slope * (x[2] + lower_slip), 0.0, 0.0],
            2.0,
        )
    });
    let initial = flow.states().to_vec();
    let mut op = prepare(&flow, law(flow.gas(), 0.01, mu, 0.0, 0.0));
    for (side, velocity, slip) in [
        (BoundarySide::ZMin, 0.0, lower_slip),
        (BoundarySide::ZMax, speed, upper_slip),
    ] {
        op = op
            .with_wall(
                side,
                MolecularWall3d::NavierSlip {
                    velocity: [velocity, 0.0, 0.0],
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m: slip,
                    fluid_heat_fraction: 0.0,
                },
            )
            .unwrap();
    }
    check(
        op.advance(
            &mut flow,
            Time::new(duration),
            DissipationControl::default(),
            || false,
        )
        .unwrap(),
    );
    for (state, before) in flow.states().chunks_exact(7).zip(initial.chunks_exact(7)) {
        for component in 0..6 {
            close(state[component], before[component], 3e-12);
        }
        close(state[6] - before[6], duration * mu * slope * slope, 3e-12);
    }
}

#[test]
fn g0_zero_slip_reuses_original_mixture_no_slip_arithmetic() {
    let grid = CartesianGrid::new(3, 3, 4, 0.25, 0.25, 0.25).unwrap();
    let mut results = Vec::new();
    for finite_slip in [false, true] {
        let mut flow = field(standard(0.0), grid, channel(2), |x| {
            let a = 0.3 + 0.1 * x[2];
            (
                [a, 1.0 - a, 0.0],
                [0.1 * x[2], -0.03, 0.0],
                2.0 + 0.1 * x[2],
            )
        });
        let mut op = prepare(&flow, law(flow.gas(), 0.01, 0.02, 0.003, 0.05));
        for (side, speed) in [(BoundarySide::ZMin, -0.1), (BoundarySide::ZMax, 0.3)] {
            let velocity = [speed, 0.0, 0.0];
            let thermal = WallThermalCondition::Isothermal(Temperature::new(2.5));
            let wall = if finite_slip {
                MolecularWall3d::NavierSlip {
                    velocity,
                    thermal,
                    slip_length_m: 0.0,
                    fluid_heat_fraction: 0.0,
                }
            } else {
                MolecularWall3d::NoSlip { velocity, thermal }
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
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            report,
            op.boundary_impulses(BoundarySide::ZMax)
                .map(|q| (q.species_masses.to_vec(), q.momentum, q.energy))
                .collect::<Vec<_>>(),
        ));
    }
    assert_eq!(results[0], results[1]);
}
