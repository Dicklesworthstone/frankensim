use super::*;
use crate::d3q19::{BoundarySpec3, CollisionModel3, Face3, FaceBoundary3};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use crate::lattice3::{D3q19, Grid3, GridBudget3};
use crate::thermal::boundary::FiniteThermalWall;
use crate::thermal::boundary::conduction::WallCell;
use crate::thermal::open::ThermalOpenLimits3;
use crate::thermal::tensor::DiffusivityTensor;

const RHO: f64 = 1.2;
const REFERENCE_RHO: f64 = 1.7;
const INITIAL_VELOCITY: [f64; 3] = [0.025, 0.0, 0.0];
const BACKGROUND: [f64; 3] = [2e-6, 0.0, -1e-6];
const GRAVITY: [f64; 3] = [0.0, -0.02, 0.0];
const BETA: f64 = 0.02;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn make<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    dims: [usize; 3],
    native: BoundarySpec3,
    faces: [ThermalOpenFace3; 6],
    initial_temperature: f64,
    expansion: f64,
    departure_limit: f64,
) -> BoussinesqOpenGrid3<L, Q> {
    let cells: usize = dims.into_iter().product();
    let periodic = core::array::from_fn(|axis| {
        matches!(native.face(Face3::ALL[2 * axis]), FaceBoundary3::Periodic)
    });
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        dims,
        collision,
        [1e-5, 2e-5, -3e-5],
        periodic,
        GridBudget3 {
            maximum_cells: 128,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| (RHO, INITIAL_VELOCITY),
    )
    .unwrap();
    let flow = OpenBoundaryGrid3::try_new_with_poll(grid, native, &mut || Ok(())).unwrap();
    let thermal = ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &vec![false; cells],
        &vec![initial_temperature; cells],
        0.83,
        ThermalOpenBoundaries3 {
            faces,
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 {
            maximum_cells: 128,
            maximum_storage_bytes: 2 * 1024 * 1024,
            temperature_range: [-2.0, 3.0],
            balance_tolerance: 1e-9,
        },
        0,
        |cell| Ok(flow.macroscopic(cell).unwrap().1),
        || Ok(()),
    )
    .unwrap();
    let law =
        BoussinesqLaw3::try_new(REFERENCE_RHO, 0.5, expansion, GRAVITY, departure_limit).unwrap();
    BoussinesqOpenGrid3::try_new_with_poll(
        flow,
        thermal,
        law,
        vec![BACKGROUND; cells],
        BoussinesqOpenLimits3 {
            maximum_cells: 128,
            maximum_force_storage_bytes: 1024 * 1024,
            force_consistency_tolerance: 1e-12,
        },
        &mut || Ok(()),
    )
    .unwrap()
}

fn periodic<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    temperature: f64,
) -> BoussinesqOpenGrid3<L, Q> {
    make(
        collision,
        [4, 3, 2],
        BoundarySpec3::periodic(),
        [ThermalOpenFace3::Periodic; 6],
        temperature,
        BETA,
        0.02,
    )
}

fn native_channel() -> BoundarySpec3 {
    BoundarySpec3::new([
        FaceBoundary3::Velocity {
            velocity: INITIAL_VELOCITY,
        },
        FaceBoundary3::Pressure { density: RHO },
        FaceBoundary3::stationary_wall(),
        FaceBoundary3::stationary_wall(),
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
    ])
}

fn channel_faces() -> [ThermalOpenFace3; 6] {
    let wall = |temperature| {
        ThermalOpenFace3::Wall(ThermalWall::FiniteSolid(FiniteThermalWall {
            temperature,
            heat_capacity: 0.2,
            resistance: 0.3,
            heat_per_step: 0.0,
        }))
    };
    [
        ThermalOpenFace3::Reservoir { temperature: 0.5 },
        ThermalOpenFace3::Reservoir { temperature: 0.7 },
        wall(0.9),
        wall(0.4),
        ThermalOpenFace3::Periodic,
        ThermalOpenFace3::Periodic,
    ]
}

fn channel(expansion: f64) -> BoussinesqOpenGrid3<D3q19, 19> {
    let mut coupled = make(
        CollisionModel3::Bgk { tau: 0.9 },
        [6, 4, 3],
        native_channel(),
        channel_faces(),
        0.5,
        expansion,
        0.02,
    );
    coupled
        .set_wall_conduction(
            1,
            &[
                WallConductionLink {
                    a: WallCell { side: 0, link: 0 },
                    b: WallCell { side: 0, link: 17 },
                    conductance: 3.0,
                },
                WallConductionLink {
                    a: WallCell { side: 0, link: 1 },
                    b: WallCell { side: 1, link: 16 },
                    conductance: 1.7,
                },
            ],
            WallConductionControls {
                max_nodes: 4,
                max_links: 2,
                max_factor_bytes: 4 * 1024 * 1024,
                balance_tolerance: 1e-9,
            },
        )
        .unwrap();
    coupled
}

fn advance<L: EquilibriumLattice3<Q>, const Q: usize>(
    coupled: &mut BoussinesqOpenGrid3<L, Q>,
    heat: &[f64],
) -> BoussinesqOpenStep3 {
    coupled
        .step_with_poll(None, heat, None, [0.0; 3], &mut || Ok(()))
        .unwrap()
}

#[derive(Debug, PartialEq)]
struct Accepted<const Q: usize> {
    flow: Vec<[u64; Q]>,
    heat: Vec<[u64; 7]>,
    force: Vec<[u64; 3]>,
    walls: [Vec<ThermalWall>; 6],
    steps: [usize; 2],
    heat_receipt: Option<ThermalOpenBalance3>,
    receipt: Option<BoussinesqOpenStep3>,
}

fn accepted<L: EquilibriumLattice3<Q>, const Q: usize>(
    coupled: &BoussinesqOpenGrid3<L, Q>,
) -> Accepted<Q> {
    Accepted {
        flow: coupled
            .flow()
            .populations()
            .iter()
            .map(|f| f.map(f64::to_bits))
            .collect(),
        heat: coupled
            .thermal()
            .populations()
            .iter()
            .map(|f| f.map(f64::to_bits))
            .collect(),
        force: coupled
            .flow()
            .forces()
            .iter()
            .map(|f| f.map(f64::to_bits))
            .collect(),
        walls: core::array::from_fn(|face| coupled.thermal().wall_states(face).unwrap().to_vec()),
        steps: [
            coupled.flow().completed_steps(),
            coupled.thermal().completed_steps(),
        ],
        heat_receipt: coupled.thermal().last_balance(),
        receipt: coupled.last_step(),
    }
}

fn check_heated_periodic<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let mut coupled = periodic::<L, Q>(collision, 0.5);
    let cells = coupled.flow().populations().len();
    let source = 0.002;
    for cell in 0..cells {
        let (rho, velocity) = coupled.flow().macroscopic(cell).unwrap();
        near(rho, RHO, 2e-14);
        for axis in 0..3 {
            near(velocity[axis], INITIAL_VELOCITY[axis], 2e-14);
        }
    }
    let buoyancy_slope = -REFERENCE_RHO * BETA * GRAVITY[1];
    for step in 1..=20 {
        let report = advance(&mut coupled, &vec![source; cells]);
        let n = step as f64;
        let expected = [
            INITIAL_VELOCITY[0] + n * BACKGROUND[0] / RHO,
            0.5 * n * n * buoyancy_slope * source / RHO,
            n * BACKGROUND[2] / RHO,
        ];
        for cell in 0..cells {
            let (rho, velocity) = coupled.flow().macroscopic(cell).unwrap();
            near(rho, RHO, 2e-12);
            near(
                coupled.thermal().temperature(cell).unwrap(),
                0.5 + n * source,
                2e-13,
            );
            for axis in 0..3 {
                near(velocity[axis], expected[axis], 2e-12);
            }
        }
        near(
            report.buoyancy_impulse[1],
            cells as f64 * buoyancy_slope * source * (n - 0.5),
            2e-13,
        );
        near(report.body_impulse[0], cells as f64 * BACKGROUND[0], 2e-14);
        near(report.body_impulse[2], cells as f64 * BACKGROUND[2], 2e-14);
        near(report.heat.combined_defect, 0.0, 2e-12);
        assert!(report.maximum_force_defect <= 1e-12);
    }
}

#[test]
fn uniform_heating_drives_reference_density_buoyancy_on_both_native_lattices() {
    check_heated_periodic::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    check_heated_periodic::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

fn check_reversal<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let mut coupled = periodic::<L, Q>(collision, 0.6);
    let initial_force = coupled.flow().forces()[0][1];
    let cells = coupled.flow().populations().len();
    let receipt = advance(&mut coupled, &vec![-0.2; cells]);
    near(coupled.flow().forces()[0][1], -initial_force, 2e-16);
    for cell in 0..cells {
        near(coupled.thermal().temperature(cell).unwrap(), 0.4, 2e-14);
        near(coupled.flow().macroscopic(cell).unwrap().1[1], 0.0, 2e-14);
    }
    near(receipt.buoyancy_impulse[1], 0.0, 2e-15);
}

#[test]
fn thermal_reversal_cancels_the_buoyant_impulse_without_losing_background_force() {
    check_reversal::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    check_reversal::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

fn total_energy(coupled: &BoussinesqOpenGrid3<D3q19, 19>) -> f64 {
    let field = coupled.thermal();
    let fluid: f64 = (0..field.populations().len())
        .map(|cell| field.temperature(cell).unwrap())
        .sum();
    let solid: f64 = (0..6)
        .flat_map(|face| field.wall_states(face).unwrap())
        .map(|wall| {
            if let ThermalWall::FiniteSolid(wall) = wall {
                wall.heat_capacity * wall.temperature
            } else {
                0.0
            }
        })
        .sum();
    fluid + solid
}

#[test]
fn conducting_walls_and_open_reservoirs_drive_two_way_convection_and_close_energy() {
    let mut active = channel(BETA);
    let mut passive = channel(0.0);
    let heat: Vec<_> = (0..72)
        .map(|cell| 0.0002 + 0.00005 * (cell % 3) as f64)
        .collect();
    let initial_wall = active.thermal().wall_states(2).unwrap().to_vec();
    for _ in 0..12 {
        for coupled in [&mut active, &mut passive] {
            let before = total_energy(coupled);
            let report = advance(coupled, &heat);
            let reservoirs: f64 = report.heat.reservoirs.iter().map(|r| r.heat_in).sum();
            near(
                total_energy(coupled) - before,
                heat.iter().sum::<f64>() + reservoirs,
                2e-11,
            );
            near(report.heat.combined_defect, 0.0, 2e-12);
            assert!(coupled.thermal().wall_conduction(1).is_some());
        }
    }
    let mut velocity_difference: f64 = 0.0;
    let mut temperature_difference: f64 = 0.0;
    for cell in 0..72 {
        let ua = active.flow().macroscopic(cell).unwrap().1;
        let up = passive.flow().macroscopic(cell).unwrap().1;
        for axis in 0..3 {
            velocity_difference = velocity_difference.max((ua[axis] - up[axis]).abs());
        }
        temperature_difference = temperature_difference.max(
            (active.thermal().temperature(cell).unwrap()
                - passive.thermal().temperature(cell).unwrap())
            .abs(),
        );
    }
    assert!(velocity_difference > 1e-7, "{velocity_difference:e}");
    assert!(temperature_difference > 1e-9, "{temperature_difference:e}");
    assert_ne!(active.thermal().wall_states(2).unwrap(), initial_wall);
}

#[test]
fn final_poll_and_density_departure_refusals_preserve_the_joint_epoch_and_exact_retry() {
    let mut actual = channel(BETA);
    let mut reference = channel(BETA);
    let heat: Vec<_> = (0..72)
        .map(|cell| 0.0002 + 0.00005 * (cell % 3) as f64)
        .collect();
    advance(&mut actual, &heat);
    advance(&mut reference, &heat);
    assert_eq!(accepted(&actual), accepted(&reference));
    let before = accepted(&actual);
    let mut polls = 0;
    let expected = reference
        .step_with_poll(None, &heat, None, [0.0; 3], &mut || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(polls > 8);
    let mut calls = 0;
    let refused = actual.step_with_poll(None, &heat, None, [0.0; 3], &mut || {
        calls += 1;
        if calls == polls {
            Err(GridError3::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        refused,
        Err(BoussinesqOpenError3::Flow(GridError3::Cancelled))
    ));
    assert_eq!(calls, polls);
    assert_eq!(accepted(&actual), before);
    let refused = actual.step_with_poll(None, &vec![2.2; 72], None, [0.0; 3], &mut || Ok(()));
    assert!(
        matches!(
            refused,
            Err(BoussinesqOpenError3::Law {
                source: BoussinesqLawError3::DensityDeparture { .. },
                ..
            }) | Err(BoussinesqOpenError3::BoundaryLaw {
                source: BoussinesqLawError3::DensityDeparture { .. },
                ..
            })
        ),
        "{refused:?}"
    );
    assert_eq!(accepted(&actual), before);
    let retried = advance(&mut actual, &heat);
    assert_eq!(retried, expected);
    assert_eq!(accepted(&actual), accepted(&reference));
}

#[test]
fn law_rejects_a_wall_trace_even_when_every_mixed_cell_temperature_is_admissible() {
    let mut faces = channel_faces();
    faces[2] = ThermalOpenFace3::Wall(ThermalWall::HeatFlux(0.3));
    faces[3] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
    let mut coupled = make::<D3q19, 19>(
        CollisionModel3::Bgk { tau: 0.9 },
        [6, 4, 3],
        native_channel(),
        faces,
        0.5,
        BETA,
        0.02,
    );
    let before = accepted(&coupled);
    let boundaries = coupled.thermal().boundaries();
    assert!(matches!(
        coupled.set_reservoir_temperature(0, 1.6),
        Err(BoussinesqOpenError3::BoundaryLaw { face: 0, .. })
    ));
    assert_eq!(coupled.thermal().boundaries(), boundaries);
    let refused = coupled.step_with_poll(None, &[0.0; 72], None, [0.0; 3], &mut || Ok(()));
    assert!(
        matches!(
            refused,
            Err(BoussinesqOpenError3::BoundaryLaw {
                face: 2,
                source: BoussinesqLawError3::DensityDeparture { .. },
            })
        ),
        "{refused:?}"
    );
    // The normal population pair has T_face=0.5+0.3/(2*1/8)=1.7.
    // Its cell mixes to about 0.8, so a cell-only law gate misses this violation.
    for cell in 0..72 {
        let temperature = coupled.thermal.trial_temperature(cell).unwrap();
        assert!(coupled.law.sample(temperature, BACKGROUND).is_ok());
    }
    assert_eq!(accepted(&coupled), before);
}

fn tensor_advance<L: EquilibriumLattice3<Q>, const Q: usize>(
    coupled: &mut BoussinesqOpenGrid3<L, Q>,
    heat: &[f64],
    diffusion: DiffusivityField<'_, 3>,
) -> BoussinesqOpenStep3 {
    coupled
        .step_with_diffusion_poll(None, heat, None, Some(diffusion), [0.0; 3], &mut || Ok(()))
        .unwrap()
}

fn check_isotropic_tensor<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let mut scalar = periodic::<L, Q>(collision, 0.5);
    let mut tensor = periodic::<L, Q>(collision, 0.5);
    let cells = scalar.flow().populations().len();
    let alpha = 0.25 * (0.83 - 0.5);
    let coefficient =
        DiffusivityTensor::from_lattice([[alpha, 0.0, 0.0], [0.0, alpha, 0.0], [0.0, 0.0, alpha]])
            .unwrap();
    let field = vec![coefficient; cells];
    let heat: Vec<_> = (0..cells).map(|i| 0.001 * (1 + i % 4) as f64).collect();
    for step in 0..12 {
        let expected = advance(&mut scalar, &heat);
        let actual = tensor_advance(
            &mut tensor,
            &heat,
            if step % 2 == 0 {
                DiffusivityField::Uniform(&coefficient)
            } else {
                DiffusivityField::PerCell(&field)
            },
        );
        near(actual.heat.combined_defect, 0.0, 2e-12);
        for axis in 0..3 {
            near(
                actual.body_impulse[axis],
                expected.body_impulse[axis],
                5e-13,
            );
        }
    }
    // A tensor field is borrowed for its own step; it cannot persist into a
    // subsequent scalar step or detach the accepted force/temperature history.
    advance(&mut scalar, &heat);
    advance(&mut tensor, &heat);
    for cell in 0..cells {
        for (&actual, &expected) in tensor.flow().populations()[cell]
            .iter()
            .zip(&scalar.flow().populations()[cell])
        {
            near(actual, expected, 5e-13);
        }
        near(
            tensor.thermal().temperature(cell).unwrap(),
            scalar.thermal().temperature(cell).unwrap(),
            5e-13,
        );
        for axis in 0..3 {
            near(
                tensor.flow().forces()[cell][axis],
                scalar.flow().forces()[cell][axis],
                5e-13,
            );
        }
    }
}

#[test]
fn boussinesq_tensor_isotropic_limit_matches_both_native_lattices() {
    check_isotropic_tensor::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    check_isotropic_tensor::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn boussinesq_tensor_rotated_conduction_changes_buoyancy_and_closes_open_wall_energy() {
    let mut oblique = channel(BETA);
    let mut diagonal = channel(BETA);
    let rotated =
        DiffusivityTensor::from_lattice([[0.09, 0.035, 0.0], [0.035, 0.08, 0.0], [0.0, 0.0, 0.06]])
            .unwrap();
    let axis_aligned =
        DiffusivityTensor::from_lattice([[0.09, 0.0, 0.0], [0.0, 0.08, 0.0], [0.0, 0.0, 0.06]])
            .unwrap();
    let heat: Vec<_> = (0..72).map(|cell| 0.0002 * (1 + cell % 5) as f64).collect();
    for _ in 0..12 {
        for (coupled, coefficient) in [(&mut oblique, &rotated), (&mut diagonal, &axis_aligned)] {
            let energy_before = total_energy(coupled);
            let report = tensor_advance(coupled, &heat, DiffusivityField::Uniform(coefficient));
            let reservoir_heat: f64 = report.heat.reservoirs.iter().map(|r| r.heat_in).sum();
            near(
                total_energy(coupled) - energy_before,
                heat.iter().sum::<f64>() + reservoir_heat,
                2e-11,
            );
            near(report.heat.combined_defect, 0.0, 2e-12);
            for cell in 0..72 {
                let temperature = coupled.thermal().temperature(cell).unwrap();
                // Independent affine constitutive law, using the transported
                // endpoint rather than an initial or scalar-only temperature.
                for axis in 0..3 {
                    let expected = BACKGROUND[axis]
                        - REFERENCE_RHO * BETA * (temperature - 0.5) * GRAVITY[axis];
                    near(coupled.flow().forces()[cell][axis], expected, 2e-15);
                }
            }
        }
    }
    let largest_temperature_change = (0..72)
        .map(|cell| {
            (oblique.thermal().temperature(cell).unwrap()
                - diagonal.thermal().temperature(cell).unwrap())
            .abs()
        })
        .fold(0.0_f64, f64::max);
    let largest_force_change = (0..72)
        .map(|cell| (oblique.flow().forces()[cell][1] - diagonal.flow().forces()[cell][1]).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        largest_temperature_change > 1e-5,
        "{largest_temperature_change:e}"
    );
    assert!(largest_force_change > 1e-8, "{largest_force_change:e}");
}

#[test]
fn boussinesq_tensor_refusal_and_late_cancellation_preserve_force_wall_and_heat_epoch() {
    let mut actual = channel(BETA);
    let mut reference = channel(BETA);
    let coefficient =
        DiffusivityTensor::from_lattice([[0.09, 0.035, 0.0], [0.035, 0.08, 0.0], [0.0, 0.0, 0.06]])
            .unwrap();
    let field = vec![coefficient; 72];
    let heat = vec![0.0004; 72];
    let before = accepted(&actual);
    let refused = actual.step_with_diffusion_poll(
        None,
        &heat,
        None,
        Some(DiffusivityField::PerCell(&field[..71])),
        [0.0; 3],
        &mut || Ok(()),
    );
    assert!(matches!(
        refused,
        Err(BoussinesqOpenError3::Thermal(ThermalOpenError3::Invalid(_)))
    ));
    assert_eq!(accepted(&actual), before);
    let mut polls = 0;
    let expected = reference
        .step_with_diffusion_poll(
            None,
            &heat,
            None,
            Some(DiffusivityField::PerCell(&field)),
            [0.0; 3],
            &mut || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut calls = 0;
    let refused = actual.step_with_diffusion_poll(
        None,
        &heat,
        None,
        Some(DiffusivityField::PerCell(&field)),
        [0.0; 3],
        &mut || {
            calls += 1;
            if calls == polls {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        refused,
        Err(BoussinesqOpenError3::Flow(GridError3::Cancelled))
    ));
    assert_eq!(calls, polls);
    assert_eq!(accepted(&actual), before);
    let retried = tensor_advance(&mut actual, &heat, DiffusivityField::PerCell(&field));
    assert_eq!(retried, expected);
    assert_eq!(accepted(&actual), accepted(&reference));
}
