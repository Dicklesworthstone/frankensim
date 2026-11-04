//! G1/G3: source-temperature viscosity reaches native curved-wall mechanics.

use super::*;
use fs_evidence::ValidityDomain;
use fs_lbm::d3q19::BoundaryStepError3;
use fs_lbm::material::{boundary::MaterialBoundaryError3, plan_scaling_from_material};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
    liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis},
    state_point::MaterialPropertySelection,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use std::cell::Cell;

const DX: f64 = 0.01;
const DT: f64 = 0.01;
const LOWER: f64 = 1.2;
const UPPER: f64 = 11.3;
const INTERFACE: f64 = 6.0;
const CELLS: usize = 4 * 12 * 4;

fn snapshot(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, 1000.0, mu).unwrap()
}

fn retuned(model: CollisionModel3, tau: f64) -> CollisionModel3 {
    match model {
        CollisionModel3::Bgk { .. } => CollisionModel3::Bgk { tau },
        CollisionModel3::Trt { magic, .. } => CollisionModel3::Trt { tau, magic },
        CollisionModel3::CentralMoment {
            higher_order_rate, ..
        } => CollisionModel3::CentralMoment {
            second_order_rate: 1.0 / tau,
            higher_order_rate,
        },
        CollisionModel3::ReducedCumulant {
            third_order_rate,
            fourth_order_rate,
            ..
        } => CollisionModel3::ReducedCumulant {
            second_order_rate: 1.0 / tau,
            third_order_rate,
            fourth_order_rate,
        },
    }
}

fn fractional(model: CollisionModel3, force: [f64; 3]) -> (BoundaryGrid3, CurvedWall3) {
    let sdf = |p: [f64; 3]| (p[1] - LOWER).min(UPPER - p[1]);
    let mut grid = BoundaryGrid3::with_collision_model(4, 12, 4, model, force, channel([0.0; 3]));
    grid.voxelize_sdf(sdf);
    let wall =
        CurvedWall3::from_level_set(&grid, sdf, 48, 1e-12, &CancelGate::new_clock_free()).unwrap();
    (grid, wall)
}

fn motion(link: CurvedWallLink3) -> [f64; 3] {
    [
        if link.intersection()[1] < INTERFACE {
            -0.01
        } else {
            0.03
        },
        0.0,
        0.0,
    ]
}

#[derive(Debug)]
struct SingleRead {
    state: FluidTransportState,
    reads: Cell<usize>,
}
impl FluidTransport for SingleRead {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        assert_eq!(self.reads.get(), 1, "resolve the original transport once");
        Ok(self.state)
    }
}

#[test]
fn unchanged_viscosity_keeps_native_bits_and_every_collision_familys_other_rates() {
    let reference = snapshot(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    // Cp/k belong to the external thermal owner. This valid source has a heat
    // diffusivity whose thermal tau would round to 1/2; momentum still works.
    let current = FluidTransportState::try_new(1000.0, 1234.0, 1e-100, 1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    for family in MODELS {
        let model = retuned(family, scaling.flow().tau);
        let force = [2e-6, -1e-6, 3e-6];
        let (mut native, wall) = fractional(model, force);
        let (mut material, _) = fractional(model, force);
        seed(&mut native);
        seed(&mut material);
        for _ in 0..3 {
            let expected = native
                .step_with_curved_walls(&wall, motion, [2.0, 3.0, -1.0], &gate)
                .unwrap();
            let actual = scaling
                .step_curved_walls_with_material(
                    &mut material,
                    &wall,
                    motion,
                    [2.0, 3.0, -1.0],
                    CELLS,
                    &gate,
                    |_| {
                        Ok(SingleRead {
                            state: current,
                            reads: Cell::new(0),
                        })
                    },
                )
                .unwrap();
            assert_eq!(actual.momentum, expected);
            assert_eq!(snapshots(&material), snapshots(&native));
            assert_eq!(material.collision_model(), model);
            assert_eq!(actual.applied.len(), 4 * 10 * 4);
            for sample in actual.applied {
                assert_eq!(sample.material.reads.get(), 1);
                assert_eq!(sample.transport, current);
                assert_eq!(sample.shear_relaxation_time, scaling.flow().tau);
            }
        }
    }
}

#[test]
fn local_material_rates_use_original_cells_force_and_non_shear_controls() {
    // G3: x-fastest material coordinates differ from native 4^3 tile order.
    // Reconstruct collision and periodic pull independently for every model.
    let reference = snapshot(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let mu = |p: [usize; 3]| 0.5 + 0.05 * (p[0] + 3 * p[1] + 7 * p[2]) as f64;
    let gate = CancelGate::new_clock_free();
    let force = [1e-6, -2e-6, 3e-6];
    for model in MODELS {
        let mut grid = BoundaryGrid3::with_collision_model(8, 4, 4, model, force, periodic());
        seed(&mut grid);
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..8 {
                    let mut f = grid.populations(x, y, z);
                    f[0] -= 0.0006;
                    f[1] += 0.0002;
                    f[7] += 0.0004;
                    grid.set_populations(x, y, z, &f);
                }
            }
        }
        let mut post = Vec::new();
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..8 {
                    let tau = 0.5 + 3.0 * (mu([x, y, z]) / 1000.0) * (DT / DX) / DX;
                    post.push(
                        collide_cell3(grid.populations(x, y, z), retuned(model, tau), force)
                            .unwrap(),
                    );
                }
            }
        }
        let wall = CurvedWall3::new(&grid, &[], &gate).unwrap();
        let report = scaling
            .step_curved_walls_with_material(
                &mut grid,
                &wall,
                |_| panic!("no selected walls"),
                [0.0; 3],
                128,
                &gate,
                |cell| Ok(snapshot(mu(cell))),
            )
            .unwrap();
        assert_eq!(grid.collision_model(), model);
        for (index, sample) in report.applied.iter().enumerate() {
            let p = [index % 8, (index / 8) % 4, index / 32];
            assert_eq!(sample.cell, p);
            let f = grid.populations(p[0], p[1], p[2]);
            for (q, c) in E3.iter().enumerate() {
                let x = (p[0] as i32 - c.0).rem_euclid(8) as usize;
                let y = (p[1] as i32 - c.1).rem_euclid(4) as usize;
                let z = (p[2] as i32 - c.2).rem_euclid(4) as usize;
                near(f[q], post[(z * 4 + y) * 8 + x][q], 3e-15);
            }
        }
    }
}

fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Dims([-3, 1, 0, 0, 0, 0])),
        ("specific_heat_capacity", 1000.0, Dims([2, 0, -2, -1, 0, 0])),
        ("dynamic_viscosity", 2.0, Dims([-1, 1, -1, 0, 0, 0])),
        ("thermal_conductivity", 1000.0, Dims([1, 1, -3, -1, 0, 0])),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        let (value, interpolation) = if name == "dynamic_viscosity" {
            (
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, value), (320.0, 1.0)],
                    dims,
                },
                InterpolationPolicy::LinearInside,
            )
        } else {
            (
                PropertyValue::Scalar { value, dims },
                InterpolationPolicy::ConstantWithinValidity,
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value,
                interpolation,
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 320.0)
                    .with("p", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic two-temperature Couette oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic lubricant".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("p", 100000.0)
            .unwrap(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn two_temperature_source_couette_matches_integrated_fluidity_and_wall_work() {
    // G1: d/dy(mu du/dy)=0 gives constant traction and
    // u(y)=u_lower + delta_u * integral_lower^y(1/mu) / integral_lower^upper(1/mu).
    // The material jump is at a cell face; both BFL fractional branches apply.
    let card = card();
    let law = law(&card);
    let cold = law.at_temperature(Temperature::new(300.0)).unwrap();
    let hot = law.at_temperature(Temperature::new(320.0)).unwrap();
    let scaling = plan_scaling_from_material(&cold, DX, DT, 0.05).unwrap();
    let (mut grid, wall) = fractional(
        CollisionModel3::Bgk {
            tau: scaling.flow().tau,
        },
        [0.0; 3],
    );
    let gate = CancelGate::new_clock_free();
    let mut last = None;
    for _ in 0..3000 {
        last = Some(
            scaling
                .step_curved_walls_with_material(
                    &mut grid,
                    &wall,
                    motion,
                    [0.0; 3],
                    CELLS,
                    &gate,
                    |p| Ok(if p[1] < 6 { &cold } else { &hot }),
                )
                .unwrap(),
        );
    }
    let report = last.unwrap();
    let (nu_cold, nu_hot) = (0.2, 0.1);
    let resistance = (INTERFACE - LOWER) / nu_cold + (UPPER - INTERFACE) / nu_hot;
    for y in 1..11 {
        let position = y as f64 + 0.5;
        let partial =
            (position.min(INTERFACE) - LOWER) / nu_cold + (position - INTERFACE).max(0.0) / nu_hot;
        let expected = -0.01 + 0.04 * partial / resistance;
        near(grid.velocity(2, y, 2)[0], expected, 2e-8);
        near(grid.velocity(2, y, 2)[1], 0.0, 2e-11);
    }
    let expected_work = -16.0 * 0.04 * 0.04 / resistance;
    near(report.momentum.curved.wall_work, expected_work, 2e-9);
    let frozen_cold_work = 16.0 * 0.04 * 0.04 * nu_cold / (UPPER - LOWER);
    assert!(-report.momentum.curved.wall_work < 0.8 * frozen_cold_work);
    for sample in report.applied {
        let source = if sample.cell[1] < 6 { &cold } else { &hot };
        assert_eq!(sample.material, source);
        near(
            sample.shear_relaxation_time,
            if sample.cell[1] < 6 { 1.1 } else { 0.8 },
            3e-16,
        );
    }
}

#[test]
fn source_density_budget_and_late_wall_failures_preserve_the_accepted_flow() {
    let card = card();
    let law = law(&card);
    let initial_source = law.reference_state();
    let scaling = plan_scaling_from_material(initial_source, DX, DT, 0.05).unwrap();
    let model = CollisionModel3::Bgk {
        tau: scaling.flow().tau,
    };
    let (mut grid, wall) = fractional(model, [1e-6, 0.0, 0.0]);
    seed(&mut grid);
    let before = snapshots(&grid);
    let gate = CancelGate::new_clock_free();
    let error = scaling
        .step_curved_walls_with_material(
            &mut grid,
            &wall,
            motion,
            [0.0; 3],
            CELLS - 1,
            &gate,
            |_| -> Result<FluidTransportState, MaterialError> {
                panic!("budget precedes source queries")
            },
        )
        .unwrap_err();
    assert!(matches!(error, MaterialBoundaryError3::Invalid(_)));
    let last_cell = [3, 10, 3];
    let mut queries = 0;
    let error = scaling
        .step_curved_walls_with_material(&mut grid, &wall, motion, [0.0; 3], CELLS, &gate, |cell| {
            queries += 1;
            let t = if cell == last_cell { 321.0 } else { 300.0 };
            law.at_temperature(Temperature::new(t))
                .map_err(|e| MaterialError::State {
                    what: e.to_string(),
                })
        })
        .unwrap_err();
    assert_eq!(queries, 160);
    assert!(matches!(error, MaterialBoundaryError3::State { cell, .. } if cell == last_cell));
    assert_eq!(snapshots(&grid), before);
    for changed in [
        FluidTransportState::try_new(999.0, 1000.0, 1000.0, 2.0).unwrap(),
        snapshot(1e-20), // numerical tau cannot represent its positive viscosity
    ] {
        let error = scaling
            .step_curved_walls_with_material(
                &mut grid,
                &wall,
                motion,
                [0.0; 3],
                CELLS,
                &gate,
                |cell| {
                    Ok(if cell == last_cell {
                        changed
                    } else {
                        snapshot(2.0)
                    })
                },
            )
            .unwrap_err();
        assert!(matches!(error, MaterialBoundaryError3::State { cell, .. } if cell == last_cell));
        assert_eq!(snapshots(&grid), before);
    }
    for cancel in [false, true] {
        let late_gate = CancelGate::new_clock_free();
        let mut calls = 0;
        let error = scaling
            .step_curved_walls_with_material(
                &mut grid,
                &wall,
                |link| {
                    calls += 1;
                    if calls == wall.links().len() {
                        if cancel {
                            late_gate.request();
                        } else {
                            return [f64::NAN; 3];
                        }
                    }
                    motion(link)
                },
                [0.0; 3],
                CELLS,
                &late_gate,
                |_| Ok(initial_source),
            )
            .unwrap_err();
        assert_eq!(calls, wall.links().len());
        assert!(matches!(error, MaterialBoundaryError3::Flow(_)));
        assert_eq!(snapshots(&grid), before);
        assert_eq!(grid.collision_model(), model);
    }
    let (mut fresh, _) = fractional(model, [1e-6, 0.0, 0.0]);
    seed(&mut fresh);
    let expected = fresh
        .step_with_curved_walls(&wall, motion, [0.0; 3], &gate)
        .unwrap();
    let retried = scaling
        .step_curved_walls_with_material(&mut grid, &wall, motion, [0.0; 3], CELLS, &gate, |_| {
            Ok(initial_source)
        })
        .unwrap();
    assert_eq!(retried.momentum, expected);
    assert_eq!(snapshots(&grid), snapshots(&fresh));
}

#[test]
fn malformed_local_rate_fields_refuse_before_publication() {
    let (mut grid, wall) = fractional(CollisionModel3::Bgk { tau: 0.8 }, [0.0; 3]);
    seed(&mut grid);
    let before = snapshots(&grid);
    let gate = CancelGate::new_clock_free();
    for rates in [vec![0.8; CELLS - 1], vec![0.8; CELLS + 1]] {
        assert!(
            grid.step_with_curved_walls_and_relaxation(&wall, motion, [0.0; 3], &rates, &gate)
                .is_err()
        );
        assert_eq!(snapshots(&grid), before);
    }
    for invalid in [0.5, f64::NAN, f64::INFINITY] {
        let mut rates = vec![0.8; CELLS];
        rates[(3 * 12 + 10) * 4 + 3] = invalid;
        assert!(matches!(
            grid.step_with_curved_walls_and_relaxation(&wall, motion, [0.0; 3], &rates, &gate),
            Err(BoundaryStepError3::Collision {
                cell: [3, 10, 3],
                ..
            })
        ));
        assert_eq!(snapshots(&grid), before);
    }
}
