//! Source-owned variable viscosity through the actual planar time integrators.
use super::*;
use crate::compressible::mixture::molecular::{
    AdaptiveDissipationControl, MolecularWall2d, StabilizedDissipationControl, TransportLaw,
    WallThermalCondition,
};

#[derive(Clone, Copy)]
enum Material {
    Constant,
    HotSpot,
    Smooth,
    MidpointSpike,
}

struct VariableTransport {
    original: ConstantTransport,
    material: Material,
}
impl TransportLaw for VariableTransport {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.original.source()
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
            self.original
                .evaluate(temperature, pressure, fractions, binary, cancelled)?;
        let excess = (temperature.value() - 1.0).max(0.0);
        match self.material {
            Material::Constant => {}
            Material::HotSpot => {
                value.dynamic_viscosity = DynViscosity::new(1.0 + 299.0 * excess);
            }
            Material::Smooth => {
                value.dynamic_viscosity = DynViscosity::new(1.0 + 2.0 * excess);
                value.bulk_viscosity = DynViscosity::new(0.1 + 0.2 * excess * excess);
            }
            Material::MidpointSpike => {
                // Manufactured contrast: cell temperatures stay near 1 and 3,
                // while their shared face law samples the large value near 2.
                if (temperature.value() - 2.0).abs() < 0.1 {
                    value.dynamic_viscosity = DynViscosity::new(1e30);
                }
            }
        }
        Ok(value)
    }
}

#[derive(Clone, Copy)]
enum Method {
    Ssp,
    Rkl,
    Adaptive,
}

fn operator(
    flow: &MixtureEuler2d,
    material: Material,
    method: Method,
    conductivity: f64,
) -> MolecularTransport2d<VariableTransport> {
    let source = PlanarTransport::new(
        flow.gas().clone(),
        VariableTransport {
            original: law(flow.gas(), 0.02, 1.0, 0.0, conductivity),
            material,
        },
        || false,
    )
    .unwrap();
    let bytes = match method {
        Method::Ssp => MolecularTransport2d::required_workspace_bytes(flow, &source),
        Method::Rkl => MolecularTransport2d::required_stabilized_workspace_bytes(flow, &source),
        Method::Adaptive => MolecularTransport2d::required_adaptive_workspace_bytes(flow, &source),
    }
    .unwrap();
    let result = MolecularTransport2d::new(flow, source, bytes, || false).unwrap();
    assert!(result.workspace_bytes() + flow.workspace_bytes() <= bytes);
    result
}

fn advance(
    flow: &mut MixtureEuler2d,
    operator: &mut MolecularTransport2d<VariableTransport>,
    method: Method,
    duration: f64,
) -> MolecularReport {
    let transport = DissipationControl {
        maximum_step_s: duration,
        ..DissipationControl::default()
    };
    let stabilized = StabilizedDissipationControl {
        transport,
        ..StabilizedDissipationControl::default()
    };
    match method {
        Method::Ssp => operator
            .advance(flow, Time::new(duration), transport, || false)
            .unwrap(),
        Method::Rkl => {
            operator
                .advance_stabilized(flow, Time::new(duration), stabilized, || false)
                .unwrap()
                .transport
        }
        Method::Adaptive => {
            operator
                .advance_adaptive(
                    flow,
                    Time::new(duration),
                    AdaptiveDissipationControl {
                        stabilized,
                        absolute_partial_density_tolerance_kg_m3: 1e-7,
                        absolute_temperature_tolerance_k: 1e-7,
                        absolute_velocity_tolerance_m_s: 1e-7,
                        ..AdaptiveDissipationControl::default()
                    },
                    || false,
                )
                .unwrap()
                .stabilized
                .transport
        }
    }
}

fn vortex() -> MixtureEuler2d {
    let a = [0.2, 0.3, 0.1, 0.3];
    let b = [0.0, 1.0, 0.0, -1.0];
    field(
        equal_mass_gas(),
        CartesianGrid::new(4, 4, 1.0, 1.0, 1.0).unwrap(),
        PERIODIC,
        |p| {
            let (x, y) = (p[0] as usize, p[1] as usize);
            (
                [0.4, 0.6, 0.0],
                [-a[x] * b[y], a[y] * b[x]],
                if x == 0 && y == 0 { 2.0 } else { 1.0 },
            )
        },
    )
}

fn kinetic(flow: &MixtureEuler2d) -> f64 {
    (0..flow.cell_count())
        .map(|i| {
            let cell = flow.cell(i).unwrap();
            0.5 * cell.momentum.iter().map(|v| v * v).sum::<f64>()
                / cell.partial_densities.iter().sum::<f64>()
        })
        .sum::<f64>()
        * flow.grid.cell_volume()
}

#[test]
fn variable_source_vortex_decays_with_cell_bounds_in_ssp_rkl_and_adaptive_stages() {
    let duration = 1e-7;
    for method in [Method::Ssp, Method::Rkl, Method::Adaptive] {
        let mut flow = vortex();
        let mut op = operator(&flow, Material::HotSpot, method, 0.0);
        assert!(op.uses_compatible_molecular_stress());
        let before = kinetic(&flow);
        let report = advance(&mut flow, &mut op, method, duration);
        check(report);
        let rate = (kinetic(&flow) - before) / duration;
        // Independent compatible-stencil value; the former native stencil's
        // same material samples give +3.613333..., an unphysical energy gain.
        close(rate, -2.366666666666667, 0.005);
        assert_eq!(report.accepted_steps, 1);
        assert_eq!(report.rejected_steps, 0);
        // Cell mu=300 gives row sum1200 and geometric factor8. Using only
        // midpoint face mu=150.5 would understate this source/stencil bound.
        assert!(report.maximum_parabolic_number >= 0.49 * 9600.0 * duration);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(op.last_step(), Some(&report));
    }
}

#[test]
fn variable_viscosity_keeps_original_species_enthalpy_and_fourier_evolution() {
    let initial = || {
        field(
            equal_mass_gas(),
            CartesianGrid::new(6, 6, 1.0 / 6.0, 1.0 / 6.0, 0.6).unwrap(),
            PERIODIC,
            |p| {
                let angle = std::f64::consts::TAU * (p[0] + p[1]);
                let a = 0.5 + 0.1 * angle.sin();
                ([a, 1.0 - a, 0.0], [0.0; 2], 1.2 + 0.1 * angle.cos())
            },
        )
    };
    for method in [Method::Ssp, Method::Rkl, Method::Adaptive] {
        let mut actual = initial();
        let mut expected = initial();
        let before = actual.states().to_vec();
        let mut variable = operator(&actual, Material::Smooth, method, 0.05);
        let mut constant = operator(&expected, Material::Constant, method, 0.05);
        check(advance(&mut actual, &mut variable, method, 1e-6));
        check(advance(&mut expected, &mut constant, method, 1e-6));
        assert_ne!(actual.states(), before);
        assert_eq!(actual.states(), expected.states());
        for i in 0..actual.cell_count() {
            assert_eq!(actual.cell(i).unwrap().partial_densities[2], 0.0);
        }
    }
}

#[test]
fn huge_unused_midpoint_viscosity_preserves_small_cell_stress_heat_and_rate() {
    let initial = || {
        field(
            equal_mass_gas(),
            CartesianGrid::new(4, 2, 1.0, 1.0, 1.0).unwrap(),
            PERIODIC,
            |p| {
                let sign = if p[0] as usize % 2 == 0 { 1.0 } else { -1.0 };
                ([1.0, 0.0, 0.0], [0.0, 0.3 + 0.1 * sign], 2.0 - sign)
            },
        )
    };
    let mut actual = initial();
    let mut expected = initial();
    let mut spike = operator(&actual, Material::MidpointSpike, Method::Ssp, 0.2);
    let mut constant = operator(&expected, Material::Constant, Method::Ssp, 0.2);
    let report = advance(&mut actual, &mut spike, Method::Ssp, 1e-4);
    check(report);
    check(advance(&mut expected, &mut constant, Method::Ssp, 1e-4));
    assert_eq!(report.accepted_steps, 1);
    close(report.maximum_parabolic_number, 32.0e-4, 1e-12);
    for (&got, &want) in actual.states().iter().zip(expected.states()) {
        close(got, want, 2e-12);
    }
    assert!(
        actual
            .cell(0)
            .unwrap()
            .primitive(actual.gas(), || false)
            .unwrap()
            .temperature
            > 1.0
    );
}

fn loads(op: &MolecularTransport2d<VariableTransport>) -> Vec<([f64; 3], [f64; 2], f64)> {
    SIDES
        .into_iter()
        .flat_map(|side| op.boundary_impulses(side))
        .map(|q| {
            (
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.species_masses[2],
                ],
                q.momentum,
                q.energy,
            )
        })
        .collect()
}

fn reflecting_vortex() -> MixtureEuler2d {
    let a = [0.2, 0.3, 0.1, 0.3];
    let b = [0.0, 1.0, 0.0, -1.0];
    field(
        equal_mass_gas(),
        CartesianGrid::new(4, 4, 1.0, 1.0, 1.0).unwrap(),
        Boundaries2d { x: WALLS, y: WALLS },
        |p| {
            let (x, y) = (p[0] as usize, p[1] as usize);
            (
                [0.4, 0.6, 0.0],
                [-a[x] * b[y], a[y] * b[x]],
                if x == 0 && y == 0 { 2.0 } else { 1.0 },
            )
        },
    )
}

fn reflecting_operator(
    flow: &MixtureEuler2d,
    method: Method,
    slip: Option<f64>,
) -> MolecularTransport2d<VariableTransport> {
    let mut op = operator(flow, Material::HotSpot, method, 0.0);
    if let Some(slip_length_m) = slip {
        for side in SIDES {
            let mut velocity = [0.0; 2];
            velocity[1 - side.axis()] = 0.1 * side.sign();
            op = op
                .with_wall(
                    side,
                    MolecularWall2d::NavierSlip {
                        velocity,
                        thermal: WallThermalCondition::Adiabatic,
                        slip_length_m,
                        fluid_heat_fraction: 0.35,
                    },
                )
                .unwrap();
        }
    }
    assert!(op.uses_compatible_molecular_stress());
    op
}

#[test]
fn late_cancellation_preserves_the_accepted_variable_source_epoch_and_slip_wall_loads() {
    let mut reference = reflecting_vortex();
    let mut reference_op = reflecting_operator(&reference, Method::Ssp, Some(0.17));
    advance(&mut reference, &mut reference_op, Method::Ssp, 1e-7);
    let controls = DissipationControl {
        maximum_step_s: 1e-7,
        ..DissipationControl::default()
    };
    let mut polls = 0;
    let expected = reference_op
        .advance(&mut reference, Time::new(2e-7), controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut flow = reflecting_vortex();
    let mut op = reflecting_operator(&flow, Method::Ssp, Some(0.17));
    advance(&mut flow, &mut op, Method::Ssp, 1e-7);
    let before = flow.states().to_vec();
    let report = *op.last_step().unwrap();
    let before_loads = loads(&op);
    let mut attempted = 0;
    assert_eq!(
        op.advance(&mut flow, Time::new(2e-7), controls, || {
            attempted += 1;
            attempted >= polls
        })
        .unwrap_err(),
        DissipationError::Cancelled
    );
    assert_eq!(flow.states(), before);
    assert_eq!(op.last_step(), Some(&report));
    assert_eq!(loads(&op), before_loads);
    assert_eq!(
        op.advance(&mut flow, Time::new(2e-7), controls, || false)
            .unwrap(),
        expected
    );
    assert_eq!(flow.states(), reference.states());
}

#[test]
fn all_reflecting_slip_laws_are_compatible_while_baths_keep_the_original_operator() {
    let flow = wall_field();
    let op = operator(&flow, Material::Smooth, Method::Ssp, 0.05);
    assert!(op.uses_compatible_molecular_stress());
    for (wall, compatible) in [
        (
            MolecularWall2d::NoSlip {
                velocity: [0.0; 2],
                thermal: WallThermalCondition::Adiabatic,
            },
            true,
        ),
        (
            MolecularWall2d::NavierSlip {
                velocity: [0.0; 2],
                thermal: WallThermalCondition::Adiabatic,
                slip_length_m: 0.0,
                fluid_heat_fraction: 0.5,
            },
            true,
        ),
        (
            MolecularWall2d::NavierSlip {
                velocity: [0.0; 2],
                thermal: WallThermalCondition::Adiabatic,
                slip_length_m: 0.1,
                fluid_heat_fraction: 0.5,
            },
            true,
        ),
    ] {
        let op = operator(&flow, Material::Smooth, Method::Ssp, 0.05)
            .with_wall(BoundarySide::XMin, wall)
            .unwrap();
        assert_eq!(op.uses_compatible_molecular_stress(), compatible);
    }
    let bath = MixtureEuler2d::new_with_reservoirs(
        flow.gas().clone(),
        flow.grid,
        flow.states().to_vec(),
        Boundaries2d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::Transmissive,
                upper: BoundaryCondition::Reflecting,
            },
            y: WALLS,
        },
        &[ReservoirBoundary {
            side: BoundarySide::XMin,
            partial_densities: &[0.6, 0.4, 0.0],
            velocity: [0.0; 2],
            pressure_pa: 0.8,
        }],
        MAX_MIXTURE_WORKSPACE_BYTES,
        || false,
    )
    .unwrap();
    assert!(
        !operator(&bath, Material::Smooth, Method::Ssp, 0.05).uses_compatible_molecular_stress()
    );
}

#[test]
fn variable_source_no_slip_dissipates_relative_to_the_actual_moving_wall_work() {
    let a = [0.2, 0.3, 0.1, 0.3];
    let b = [0.0, 1.0, 0.0, -1.0];
    let duration = 1e-8;
    for method in [Method::Ssp, Method::Rkl, Method::Adaptive] {
        for (lid_speed, dissipation) in [(0.0, 10.131666666666667), (0.4, 9.971666666666667)] {
            let mut flow = field(
                equal_mass_gas(),
                CartesianGrid::new(4, 4, 1.0, 1.0, 1.0).unwrap(),
                Boundaries2d { x: WALLS, y: WALLS },
                |p| {
                    let (x, y) = (p[0] as usize, p[1] as usize);
                    (
                        [0.4, 0.6, 0.0],
                        [-a[x] * b[y], a[y] * b[x]],
                        if x == 0 && y == 0 { 2.0 } else { 1.0 },
                    )
                },
            );
            let mut op = operator(&flow, Material::HotSpot, method, 0.0);
            for side in SIDES {
                let velocity = if side == BoundarySide::YMax {
                    [lid_speed, 0.0]
                } else {
                    [0.0; 2]
                };
                op = op
                    .with_wall(
                        side,
                        MolecularWall2d::NoSlip {
                            velocity,
                            thermal: WallThermalCondition::Adiabatic,
                        },
                    )
                    .unwrap();
            }
            assert!(op.uses_compatible_molecular_stress());
            let before = kinetic(&flow);
            check(advance(&mut flow, &mut op, method, duration));
            let mut outward_work = 0.0;
            for side in SIDES {
                let speed = if side == BoundarySide::YMax {
                    lid_speed
                } else {
                    0.0
                };
                for load in op.boundary_impulses(side) {
                    assert_eq!(load.species_masses, [0.0; 3]);
                    close(load.energy, speed * load.momentum[0], 2e-17);
                    outward_work += load.energy;
                }
            }
            // Independent cell quadratic form, with odd wall traces and the
            // physical 3-D deviatoric trace. It includes boundary/corner cross
            // derivatives; the old stationary-wall operator gives -8.66666...
            // instead of -10.13166... for this original hot-spot source.
            close(
                (kinetic(&flow) - before + outward_work) / duration,
                -dissipation,
                0.003,
            );
            assert_eq!(flow.time().value(), 0.0);
        }
    }
}

#[test]
fn heated_no_slip_preserves_original_species_enthalpy_and_heat_without_unused_wall_stiffness() {
    let initial = || {
        field(
            equal_mass_gas(),
            CartesianGrid::new(4, 3, 1.0, 1.0, 1.0).unwrap(),
            Boundaries2d { x: WALLS, y: WALLS },
            |p| {
                let y = 0.4 + 0.03 * p[0];
                ([y, 1.0 - y, 0.0], [0.0; 2], 1.0)
            },
        )
    };
    for method in [Method::Ssp, Method::Rkl, Method::Adaptive] {
        let mut actual = initial();
        let mut expected = initial();
        let before = actual.states().to_vec();
        let mut spike = operator(&actual, Material::MidpointSpike, method, 0.2);
        let mut constant = operator(&expected, Material::Constant, method, 0.2);
        for side in SIDES {
            let wall = MolecularWall2d::NoSlip {
                velocity: [0.0; 2],
                thermal: WallThermalCondition::Isothermal(Temperature::new(3.0)),
            };
            spike = spike.with_wall(side, wall).unwrap();
            constant = constant.with_wall(side, wall).unwrap();
        }
        // Cell mu=1; the actual wall heat law samples mu=1e30 at mean T=2.
        // That unused mechanical face value cannot erase heat or dictate dt.
        let report = advance(&mut actual, &mut spike, method, 1e-5);
        check(report);
        check(advance(&mut expected, &mut constant, method, 1e-5));
        assert_eq!(actual.states(), expected.states());
        assert_eq!(loads(&spike), loads(&constant));
        assert_ne!(actual.states(), before);
        assert!(report.maximum_parabolic_number < 0.001);
    }
}

#[test]
fn large_navier_slip_converges_to_natural_free_slip_through_every_planar_integrator() {
    for method in [Method::Ssp, Method::Rkl, Method::Adaptive] {
        let mut free = reflecting_vortex();
        let initial_kinetic = kinetic(&free);
        let mut free_op = reflecting_operator(&free, method, None);
        let mut finite = reflecting_vortex();
        let mut finite_op = reflecting_operator(&finite, method, Some(1e200));
        check(advance(&mut free, &mut free_op, method, 1e-7));
        check(advance(&mut finite, &mut finite_op, method, 1e-7));
        assert!(kinetic(&free) < initial_kinetic);
        for (a, b) in free.states().iter().zip(finite.states()) {
            close(*a, *b, 1e-13);
        }
        for (a, b) in loads(&free_op).iter().zip(loads(&finite_op)) {
            assert_eq!(a.0, [0.0; 3]);
            assert_eq!(b.0, [0.0; 3]);
            for component in 0..2 {
                close(a.1[component], b.1[component], 1e-13);
            }
            close(a.2, b.2, 1e-13);
            assert_eq!(a.2, 0.0);
        }
    }
}
