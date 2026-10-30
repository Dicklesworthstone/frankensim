use super::*;
use crate::compressible as planar;
use crate::compressible::spatial::{
    AxisBoundaries, Boundaries3d, BoundaryCondition, BoundarySide, CartesianGrid, Euler3d,
    FrozenGas, State, StepControl, Totals,
};
use crate::compressible::tests::{close, gas};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::{DynViscosity, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const WALL: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const WALLS: Boundaries3d = Boundaries3d {
    x: WALL,
    y: WALL,
    z: WALL,
};

fn law(g: &FrozenGas, mu: f64, bulk: f64, conductivity: f64) -> ConstantTransport {
    let models = g.source().components().to_vec();
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap();
    ConstantTransport::new(
        source,
        "manufactured-spatial-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(mu),
            bulk_viscosity: DynViscosity::new(bulk),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        vec![0.0],
    )
    .unwrap()
}
fn field(
    g: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    initial: impl Fn([f64; 3]) -> [f64; 5],
) -> Euler3d {
    let count = grid.shape().iter().product();
    let states = (0..count)
        .map(|i| {
            let p = initial(grid.cell_center(i).unwrap());
            State::from_primitive(&g, p[0], [p[1], p[2], p[3]], p[4]).unwrap()
        })
        .collect();
    Euler3d::new(g, grid, states, boundaries, || false).unwrap()
}
fn cube(n: usize, initial: impl Fn([f64; 3]) -> [f64; 5]) -> Euler3d {
    let h = 1.0 / n as f64;
    field(
        gas(0.0),
        CartesianGrid::new(n, n, n, h, h, h).unwrap(),
        PERIODIC,
        initial,
    )
}
fn partition(flow: &Euler3d) -> (f64, f64) {
    let kinetic: f64 = flow
        .states()
        .iter()
        .map(|q| {
            0.5 * q.momentum.iter().map(|p| p * p).sum::<f64>() / q.density
                * flow.grid().cell_volume()
        })
        .sum();
    (kinetic, flow.totals().energy - kinetic)
}

#[test]
fn g1_spatial_fourier_heat_uses_all_physical_axes_and_refines_quadratically() {
    let lengths = [1.0, 1.5, 0.75];
    let k = lengths.map(|l| 2.0 * std::f64::consts::PI / l);
    let mode = |x: [f64; 3]| (k[0] * x[0]).sin() * (k[1] * x[1]).cos() * (k[2] * x[2]).sin();
    let end = 0.03;
    let mut errors = Vec::new();
    for n in [6, 12, 18] {
        let h = lengths.map(|l| l / n as f64);
        let mut flow = field(
            gas(0.0),
            CartesianGrid::new(n, n, n, h[0], h[1], h[2]).unwrap(),
            PERIODIC,
            |x| [1.0, 0.0, 0.0, 0.0, 2.0 + 0.2 * mode(x)],
        );
        let original = flow.totals();
        let mut op =
            MolecularTransport3d::new(&flow, law(flow.gas(), 0.0, 0.0, 0.05), || false).unwrap();
        let report = op
            .advance(
                &mut flow,
                Time::new(end),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        // Original source Cp/R=3.5 and R=1 give Cv=2.5; molecular diffusion is
        // kappa/(rho Cv), with all three physical reciprocal lengths present.
        let amplitude = 0.2 * (-0.05 / 2.5 * k.iter().map(|x| x * x).sum::<f64>() * end).exp();
        let error = (flow
            .states()
            .iter()
            .enumerate()
            .map(|(i, q)| {
                assert_eq!(q.density, 1.0);
                assert_eq!(q.momentum, [0.0; 3]);
                let expected = 2.0 + amplitude * mode(flow.grid().cell_center(i).unwrap());
                let actual = q.primitive(flow.gas(), || false).unwrap().temperature;
                (actual - expected).powi(2)
            })
            .sum::<f64>()
            / flow.states().len() as f64)
            .sqrt();
        errors.push(error);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(report.elapsed.value(), end);
        assert!(report.maximum_parabolic_number <= 0.40000000000001);
        close(flow.totals().energy, original.energy, 4e-13);
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 4e-13);
        }
    }
    eprintln!("spatial thermal mode L2 errors (6/12/18): {errors:?}");
    assert!(errors[0] < 0.002);
    for (i, ratio) in [2.0_f64, 1.5].into_iter().enumerate() {
        let order = (errors[i] / errors[i + 1]).ln() / ratio.ln();
        assert!(
            (1.8..2.2).contains(&order),
            "thermal refinement {errors:?}, order={order}"
        );
    }
}

#[test]
fn g1_spatial_oblique_shear_and_compression_include_every_cross_derivative() {
    let n = 8;
    let h = 1.0 / n as f64;
    let w = 2.0 * std::f64::consts::PI;
    let (mu, bulk, end) = (0.05, 0.03, 0.02);
    // Independent Fourier symbol of the centered Newtonian face stress.
    // All wave-vector components are nonzero; neither polarization suppresses
    // the third velocity. The off-diagonal symbol detects missing cross terms.
    let lambda = 4.0 * (std::f64::consts::PI / n as f64).sin().powi(2) / (h * h);
    let chi = (w / n as f64).sin().powi(2) / (h * h);
    for direction in [[1.0, 1.0, 1.0], [1.0, 1.0, -2.0]] {
        let longitudinal = direction[2] == 1.0;
        let mut flow = cube(n, |x| {
            let u = 0.1 * (w * (x[0] + x[1] + x[2])).sin();
            [
                1.0,
                direction[0] * u,
                direction[1] * u,
                direction[2] * u,
                2.0,
            ]
        });
        let initial = flow.totals();
        let (kinetic0, internal0) = partition(&flow);
        let mut op =
            MolecularTransport3d::new(&flow, law(flow.gas(), mu, bulk, 0.0), || false).unwrap();
        let controls = MolecularControl {
            maximum_step_s: 0.001,
            ..MolecularControl::default()
        };
        op.advance(&mut flow, Time::new(end), controls, || false)
            .unwrap();
        let rate = 3.0 * mu * lambda
            + (mu / 3.0 + bulk)
                * if longitudinal {
                    lambda + 2.0 * chi
                } else {
                    lambda - chi
                };
        let amplitude = 0.1 * (-rate * end).exp();
        for (i, q) in flow.states().iter().enumerate() {
            let x = flow.grid().cell_center(i).unwrap();
            let exact = amplitude * (w * (x[0] + x[1] + x[2])).sin();
            for axis in 0..3 {
                close(q.momentum[axis], direction[axis] * exact, 1e-6);
            }
        }
        let (kinetic1, internal1) = partition(&flow);
        assert!(kinetic1 < kinetic0 && internal1 > internal0);
        close(internal1 - internal0, kinetic0 - kinetic1, 3e-13);
        close(flow.totals().energy, initial.energy, 3e-13);
        for axis in 0..3 {
            close(flow.totals().momentum[axis], initial.momentum[axis], 3e-13);
        }
    }
}

#[test]
fn g3_spatial_molecular_transport_reduces_to_the_original_planar_extrusion() {
    let (nx, ny, nz) = (7, 5, 3);
    let (dx, dy, dz) = (0.12, 0.17, 0.3);
    let g = gas(0.4);
    let init = |x: [f64; 2]| {
        [
            1.0 + 0.05 * (3.0 * x[0]).sin(),
            0.12 * (2.0 * x[1]).cos(),
            0.2 * (4.0 * x[0]).sin(),
            2.0 + 0.15 * (x[0] + x[1]).sin(),
        ]
    };
    let grid2 = planar::CartesianGrid::new(nx, ny, dx, dy, nz as f64 * dz).unwrap();
    let states2 = (0..nx * ny)
        .map(|i| {
            let p = init(grid2.cell_center(i).unwrap());
            planar::State::from_primitive(&g, p[0], [p[1], p[2]], p[3]).unwrap()
        })
        .collect();
    let mut plane = planar::Euler2d::new(
        g.clone(),
        grid2,
        states2,
        planar::Boundaries2d {
            x: AxisBoundaries::Periodic,
            y: AxisBoundaries::Periodic,
        },
        || false,
    )
    .unwrap();
    let mut space = field(
        g.clone(),
        CartesianGrid::new(nx, ny, nz, dx, dy, dz).unwrap(),
        PERIODIC,
        |x| {
            let p = init([x[0], x[1]]);
            [p[0], p[1], p[2], 0.0, p[3]]
        },
    );
    let mut op2 =
        planar::molecular::MolecularTransport2d::new(&plane, law(&g, 0.02, 0.01, 0.03), || false)
            .unwrap();
    let mut op3 = MolecularTransport3d::new(&space, law(&g, 0.02, 0.01, 0.03), || false).unwrap();
    // Common small duration isolates the spatial reduction from the legitimately
    // different conservative 2D/3D parabolic timestep estimates.
    let controls = MolecularControl {
        maximum_step_s: 0.001,
        ..MolecularControl::default()
    };
    op2.advance(&mut plane, Time::new(0.001), controls, || false)
        .unwrap();
    op3.advance(&mut space, Time::new(0.001), controls, || false)
        .unwrap();
    for (i, q) in space.states().iter().enumerate() {
        let expected = plane.states()[i % (nx * ny)];
        close(q.density, expected.density, 2e-13);
        close(q.energy, expected.energy, 2e-12);
        close(q.momentum[0], expected.momentum[0], 2e-13);
        close(q.momentum[1], expected.momentum[1], 2e-13);
        assert_eq!(q.momentum[2], 0.0);
    }
    close(space.totals().energy, plane.totals().energy, 2e-12);
}

#[test]
fn g0_spatial_zero_law_preserves_euler_and_source_or_work_refusal_is_atomic() {
    let init = |x: [f64; 3]| [1.0, 0.1 * x[0], -0.1 * x[1], 0.2 * x[2], 2.0];
    let mut reference = cube(4, init);
    let initial = cube(4, init);
    let mut coupled =
        MolecularEuler3d::new(initial, law(reference.gas(), 0.0, 0.0, 0.0), || false).unwrap();
    let expected = reference
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let actual = coupled
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(actual.advection, expected);
    assert_eq!(coupled.flow().states(), reference.states());
    assert!(matches!(
        MolecularTransport3d::new(&reference, law(&gas(0.4), 0.1, 0.0, 0.1), || false),
        Err(MolecularError::Invalid(_))
    ));
    let mut flow = cube(4, init);
    let old = flow.states().to_vec();
    let mut op =
        MolecularTransport3d::new(&flow, law(flow.gas(), 0.1, 0.0, 0.1), || false).unwrap();
    let controls = MolecularControl {
        maximum_steps: 1,
        maximum_step_s: 1e-4,
        ..MolecularControl::default()
    };
    assert!(matches!(
        op.advance(&mut flow, Time::new(0.01), controls, || false),
        Err(MolecularError::WorkLimit { .. })
    ));
    assert_eq!(flow.states(), old);
    assert_eq!(flow.time().value(), 0.0);
    for side in BoundarySide::ALL {
        assert!(
            op.boundary_impulses(side)
                .iter()
                .all(|q| *q == Totals::default())
        );
    }
}

fn coupled_walls() -> MolecularEuler3d<ConstantTransport> {
    let flow = field(
        gas(0.4),
        CartesianGrid::new(4, 3, 2, 0.25, 0.2, 0.3).unwrap(),
        WALLS,
        |x| {
            [
                1.0 + 0.1 * x[0],
                0.3 + 0.1 * x[0],
                -0.2 + 0.2 * x[1],
                0.4 - 0.3 * x[2],
                2.0,
            ]
        },
    );
    let coefficients = law(flow.gas(), 0.02, 0.01, 0.05);
    MolecularEuler3d::new(flow, coefficients, || false).unwrap()
}

#[test]
fn g2_g4_all_six_walls_receive_complete_impulse_and_late_cancel_retains_history() {
    let mut probe = coupled_walls();
    probe
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    let before = probe.flow().totals();
    let mut calls = 0;
    let report = probe
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            calls += 1;
            false
        })
        .unwrap();
    assert_eq!(report.first.elapsed.value() * 2.0, report.elapsed.value());
    assert_eq!(report.second.elapsed, report.first.elapsed);
    assert_eq!(report.advection.elapsed, report.elapsed);
    assert_eq!(report.boundary_change.mass, 0.0);
    assert_eq!(report.boundary_change.energy, 0.0);
    close(probe.flow().totals().energy, before.energy, 3e-13);
    let mut impulse = [0.0; 3];
    for side in BoundarySide::ALL {
        for q in probe.boundary_impulses(side) {
            for axis in 0..3 {
                impulse[axis] += q.momentum[axis];
            }
        }
    }
    for axis in 0..3 {
        close(
            probe.flow().totals().momentum[axis] - before.momentum[axis] + impulse[axis],
            0.0,
            3e-13,
        );
        assert!(
            (report.first.boundary_change.momentum[axis]
                + report.second.boundary_change.momentum[axis])
                .abs()
                > 1e-9
        );
    }
    for cutoff in [calls / 2, calls] {
        let mut flow = coupled_walls();
        flow.step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let old = flow.flow().states().to_vec();
        let time = flow.flow().time();
        let previous = flow.last_step().copied();
        let impulses = BoundarySide::ALL.map(|s| flow.boundary_impulses(s).to_vec());
        let mut polls = 0;
        assert_eq!(
            flow.step(Time::new(0.001), MolecularStepControl::default(), || {
                polls += 1;
                polls == cutoff
            }),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(flow.flow().states(), old);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step().copied(), previous);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(flow.boundary_impulses(side), impulses[i]);
        }
        let next = flow
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(next, report);
        assert_eq!(flow.flow().states(), probe.flow().states());
    }
}
