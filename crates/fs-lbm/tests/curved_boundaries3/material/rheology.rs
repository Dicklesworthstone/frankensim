//! G1/G3/G4: the original curved-wall owner consumes sourced nonlinear viscosity.
use super::*;
use fs_lbm::material::boundary::rheology::{MaterialBoundaryRheologyPoint3, RheologySolve};
use fs_material::liquid::resolve_liquid_state;
use fs_qty::Frequency;

fn controls() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 64,
        relaxation: 0.75,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-12,
    }
}

#[test]
fn constant_source_preserves_curved_populations_and_accounts_for_every_collision_family() {
    let reference = snapshot(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    // A valid transport response whose thermal tau would round to 1/2 must
    // still be usable by the original momentum-only admission.
    let source = FluidTransportState::try_new(1000.0, 1234.0, 1e-100, 1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    for family in MODELS {
        let model = retuned(family, scaling.flow().tau);
        let (mut native, wall) = fractional(model, [2e-6, -1e-6, 3e-6]);
        let (mut actual, _) = fractional(model, [2e-6, -1e-6, 3e-6]);
        seed(&mut native);
        seed(&mut actual);
        let mut queries = 0;
        for _ in 0..2 {
            let expected = native
                .step_with_curved_walls(&wall, motion, [2.0, 3.0, -1.0], &gate)
                .unwrap();
            let report = scaling
                .step_curved_walls_with_material_rheology(
                    &mut actual,
                    &wall,
                    motion,
                    [2.0, 3.0, -1.0],
                    CELLS,
                    controls(),
                    &gate,
                    |_| {
                        queries += 1;
                        Ok(SingleRead {
                            state: source,
                            reads: Cell::new(0),
                        })
                    },
                )
                .unwrap();
            assert_eq!(report.momentum, expected);
            assert_eq!(snapshots(&actual), snapshots(&native));
            assert_eq!(actual.collision_model(), model);
            for sample in report.applied {
                assert_eq!(sample.material.source.reads.get(), 1);
                assert_eq!(sample.material.iterations, 1);
                assert_eq!(sample.material.point.cell, sample.cell);
                assert_eq!(sample.transport, source);
                assert_eq!(sample.shear_relaxation_time, scaling.flow().tau);
            }
        }
        assert_eq!(queries, 2 * 160);
    }
}

fn query(temperature: f64, rate: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), temperature)
        .unwrap()
        .with_quantity("shear", QuantitySpec::dimensional(Frequency::DIMS), rate)
        .unwrap()
}

fn source_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), 1000.0),
        ("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0]), 1000.0),
        ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), 1000.0),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), 1.0),
    ] {
        let viscosity = name == "dynamic_viscosity";
        claims
            .insert_claim(PropertyClaim {
                key: if name == "thermal_conductivity" {
                    PropertyKey::with_quantity(
                        name,
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::ThermalConductivity,
                            ValueForm::Static,
                        )),
                    )
                } else {
                    PropertyKey::new(name, dims)
                },
                value: if viscosity {
                    PropertyValue::Surface2D {
                        axes: ["T".into(), "shear".into()],
                        axis_dims: [Temperature::DIMS, Frequency::DIMS],
                        samples: [[300.0, 0.0], [320.0, 0.0], [300.0, 100.0], [320.0, 100.0]]
                            .into_iter()
                            .map(|p| (p, 1.0 + 0.02 * (p[0] - 300.0) + 0.2 * p[1]))
                            .collect(),
                        triangles: vec![[0, 1, 2], [1, 3, 2]],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                interpolation: if viscosity {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.0,
                        320.0,
                    )
                    .with_quantity(
                        "shear",
                        QuantitySpec::dimensional(Frequency::DIMS),
                        0.0,
                        100.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured curved-wall temperature/shear law".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured generalized-Newtonian liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn temperature(cell: [usize; 3]) -> f64 {
    300.0 + 1.5 * cell[1] as f64
}

// Independently prescribe all six symmetric strain components and the raw
// Guo-frame moments with a Hermite projection. No production stress helper
// supplies the reference. The trace is zero, and rho/u/F vary independently.
fn tensor(cell: [usize; 3], factor: f64) -> [[f64; 3]; 3] {
    let a = factor * 0.0004 * (1 + cell[0] + 2 * cell[1] + 3 * cell[2]) as f64;
    [
        [0.5 * a, 0.4 * a, -0.15 * a],
        [0.4 * a, -0.2 * a, 0.25 * a],
        [-0.15 * a, 0.25 * a, -0.3 * a],
    ]
}

fn prescribed(grid: &mut BoundaryGrid3, factor: f64) {
    grid.set_force_field(|[x, y, z]| {
        [1e-4 * x as f64, -2e-5 * y as f64, 3e-5 * z as f64]
    });
    for z in 0..4 {
        for y in 1..11 {
            for x in 0..4 {
                let cell = [x, y, z];
                let rho = 0.99 + 0.004 * (x + z) as f64;
                let u = [0.01, -0.005, 0.003];
                let force = grid.force_density(x, y, z);
                let d = tensor(cell, factor);
                let equilibrium = equilibrium3(rho, u);
                let f = core::array::from_fn(|q| {
                    let (cx, cy, cz) = E3[q];
                    let c = [f64::from(cx), f64::from(cy), f64::from(cz)];
                    let momentum =
                        -0.5 * (c[0] * force[0] + c[1] * force[1] + c[2] * force[2]);
                    let mut contraction = 0.0;
                    for a in 0..3 {
                        for b in 0..3 {
                            let pi = -2.0 * rho * 0.8 * d[a][b] / 3.0
                                - 0.5 * (u[a] * force[b] + force[a] * u[b]);
                            contraction +=
                                (c[a] * c[b] - if a == b { 1.0 / 3.0 } else { 0.0 }) * pi;
                        }
                    }
                    equilibrium[q] + W3[q] * (3.0 * momentum + 4.5 * contraction)
                });
                grid.set_populations(x, y, z, &f);
            }
        }
    }
}

#[test]
fn sourced_thermal_shear_law_matches_independent_full_stress_root_and_curved_step() {
    let card = source_card();
    let reference = resolve_liquid_state(
        &card,
        &query(300.0, 0.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    for family in MODELS {
        for factor in [0.0, 1.0] {
            // A deliberately stale constructor tau must be only an initial guess.
            let model = retuned(family, 1.7);
            let force = [0.001, -0.0007, 0.0003];
            let (mut grid, wall) = fractional(model, force);
            let (mut manual, _) = fractional(model, force);
            prescribed(&mut grid, factor);
            prescribed(&mut manual, factor);
            let report = scaling
                .step_curved_walls_with_material_rheology(
                    &mut grid,
                    &wall,
                    motion,
                    [0.7, -1.3, 2.1],
                    CELLS,
                    controls(),
                    &gate,
                    |point| {
                        reference
                            .at_query_point(
                                &card,
                                &query(temperature(point.cell), point.shear_rate_per_s),
                            )
                            .map_err(|error| MaterialError::State {
                                what: error.to_string(),
                            })
                    },
                )
                .unwrap();
            let mut exact_tau = vec![1.0; CELLS];
            for sample in report.applied {
                let point = sample.material.point;
                let d = tensor(point.cell, factor);
                let gamma = (2.0 * d.iter().flatten().map(|x| x * x).sum::<f64>()).sqrt();
                let a = 0.5 + 0.3 * (1.0 + 0.02 * (temperature(point.cell) - 300.0));
                let b = 0.3 * 0.2 * 0.8 * gamma / DT;
                let tau = 0.5 * (a + (a * a + 4.0 * b).sqrt());
                near(sample.shear_relaxation_time, tau, 5e-12);
                near(point.shear_rate_per_s, 0.8 * gamma / (DT * tau), 8e-11);
                let tolerance = controls().absolute_relaxation_tolerance
                    + controls().relative_relaxation_tolerance
                        * (sample.shear_relaxation_time + sample.material.relaxation_defect);
                assert!(sample.material.relaxation_defect <= tolerance);
                let source = sample.material.source.material();
                assert_eq!(source.card_identity(), card.content_hash());
                assert_eq!(source.properties().len(), 4);
                assert_eq!(
                    source.query(),
                    query(temperature(point.cell), point.shear_rate_per_s)
                );
                assert_eq!(
                    source.property("dynamic_viscosity").unwrap().answer().receipt.selected,
                    reference
                        .material()
                        .property("dynamic_viscosity")
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
                let [x, y, z] = sample.cell;
                exact_tau[(z * 12 + y) * 4 + x] = tau;
            }
            let expected = manual
                .step_with_curved_walls_and_relaxation(
                    &wall,
                    motion,
                    [0.7, -1.3, 2.1],
                    &exact_tau,
                    &gate,
                )
                .unwrap();
            for (actual, expected) in snapshots(&grid)
                .iter()
                .flatten()
                .zip(snapshots(&manual).iter().flatten())
            {
                near(*actual, *expected, 2e-12);
            }
            near3(
                report.momentum.curved.wall_impulse,
                expected.curved.wall_impulse,
                2e-10,
            );
            near3(
                report.momentum.curved.wall_angular_impulse,
                expected.curved.wall_angular_impulse,
                2e-9,
            );
            near(report.momentum.curved.wall_work, expected.curved.wall_work, 2e-11);
            assert_eq!(report.momentum.body_force_impulse, expected.body_force_impulse);
            assert_eq!(grid.collision_model(), model);
        }
    }
}

#[test]
fn carreau_couette_uses_physical_shear_for_velocity_and_wall_work() {
    // Carreau has a finite zero-rate plateau; no artificial viscosity cap or
    // shear floor is needed at rest. The exact steady Couette rate is Δu/H.
    let viscosity = |rate: f64| 0.5 + 2.5 * (1.0 + (8.0 * rate).powi(2)).powf(-0.25);
    let reference = snapshot(viscosity(0.0));
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
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
                .step_curved_walls_with_material_rheology(
                    &mut grid,
                    &wall,
                    motion,
                    [0.0; 3],
                    CELLS,
                    controls(),
                    &gate,
                    |point| Ok(snapshot(viscosity(point.shear_rate_per_s))),
                )
                .unwrap(),
        );
    }
    let report = last.unwrap();
    let rate = 0.04 / ((UPPER - LOWER) * DT);
    let nu_lattice = viscosity(rate) / 1000.0 * DT / (DX * DX);
    for y in 1..11 {
        let expected = -0.01 + 0.04 * (y as f64 + 0.5 - LOWER) / (UPPER - LOWER);
        near(grid.velocity(2, y, 2)[0], expected, 5e-8);
    }
    for sample in report.applied {
        near(sample.material.point.shear_rate_per_s, rate, 3e-6);
        near(sample.transport.dynamic_viscosity_pa_s(), viscosity(rate), 3e-6);
    }
    let expected_work = -16.0 * nu_lattice * 0.04_f64.powi(2) / (UPPER - LOWER);
    near(report.momentum.curved.wall_work, expected_work, 2e-8);
    let frozen_work = 16.0 * 0.3 * 0.04_f64.powi(2) / (UPPER - LOWER);
    assert!(-report.momentum.curved.wall_work < 0.8 * frozen_work);
}

#[test]
fn source_density_budget_root_and_late_cancellation_refuse_atomically_then_retry() {
    let reference = snapshot(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let model = CollisionModel3::Bgk {
        tau: scaling.flow().tau,
    };
    let (mut grid, wall) = fractional(model, [0.001, -0.0007, 0.0003]);
    prescribed(&mut grid, 1.0);
    let saved = snapshots(&grid);
    let gate = CancelGate::new_clock_free();
    for (budget, solve) in [
        (CELLS - 1, controls()),
        (
            CELLS,
            RheologySolve {
                maximum_iterations: 0,
                ..controls()
            },
        ),
    ] {
        let error = scaling
            .step_curved_walls_with_material_rheology(
                &mut grid,
                &wall,
                motion,
                [0.0; 3],
                budget,
                solve,
                &gate,
                |_| -> Result<FluidTransportState, MaterialError> {
                    panic!("admission precedes source work")
                },
            )
            .unwrap_err();
        assert!(matches!(error, MaterialBoundaryError3::Invalid(_)));
        assert_eq!(snapshots(&grid), saved);
    }
    let error = scaling
        .step_curved_walls_with_material_rheology(
            &mut grid,
            &wall,
            motion,
            [0.0; 3],
            CELLS,
            RheologySolve {
                maximum_iterations: 1,
                ..controls()
            },
            &gate,
            |_| Ok(snapshot(2.0)),
        )
        .unwrap_err();
    assert!(matches!(
        error,
        MaterialBoundaryError3::RheologyConvergence {
            cell: [0, 1, 0],
            iterations: 1,
            ..
        }
    ));
    assert_eq!(snapshots(&grid), saved);
    let last = [3, 10, 3];
    for density in [false, true] {
        let error = scaling
            .step_curved_walls_with_material_rheology(
                &mut grid,
                &wall,
                motion,
                [0.0; 3],
                CELLS,
                controls(),
                &gate,
                |point| {
                    if point.cell == last {
                        if density {
                            return FluidTransportState::try_new(999.0, 1000.0, 1000.0, 1.0);
                        }
                        return Err(MaterialError::State {
                            what: "original source domain exhausted".into(),
                        });
                    }
                    Ok(reference)
                },
            )
            .unwrap_err();
        assert!(matches!(error, MaterialBoundaryError3::State { cell, .. } if cell == last));
        assert_eq!(snapshots(&grid), saved);
    }
    for cancel_in_source in [false, true] {
        let stopped = CancelGate::new_clock_free();
        let mut walls = 0;
        let error = scaling
            .step_curved_walls_with_material_rheology(
                &mut grid,
                &wall,
                |link| {
                    walls += 1;
                    if !cancel_in_source && walls == wall.links().len() {
                        stopped.request();
                    }
                    motion(link)
                },
                [0.0; 3],
                CELLS,
                controls(),
                &stopped,
                |point| {
                    if cancel_in_source && point.cell == last {
                        stopped.request();
                    }
                    Ok(reference)
                },
            )
            .unwrap_err();
        assert!(matches!(
            error,
            MaterialBoundaryError3::Flow(BoundaryStepError3::Cancelled)
        ));
        assert_eq!(snapshots(&grid), saved);
        assert_eq!(grid.collision_model(), model);
    }
    let (mut fresh, _) = fractional(model, [0.001, -0.0007, 0.0003]);
    prescribed(&mut fresh, 1.0);
    let resolve = |_: MaterialBoundaryRheologyPoint3| Ok(reference);
    let expected = scaling
        .step_curved_walls_with_material_rheology(
            &mut fresh,
            &wall,
            motion,
            [0.0; 3],
            CELLS,
            controls(),
            &gate,
            resolve,
        )
        .unwrap();
    let retried = scaling
        .step_curved_walls_with_material_rheology(
            &mut grid,
            &wall,
            motion,
            [0.0; 3],
            CELLS,
            controls(),
            &gate,
            resolve,
        )
        .unwrap();
    assert_eq!(expected, retried);
    assert_eq!(snapshots(&grid), snapshots(&fresh));
}
