use super::super::{BoundaryCondition, State, StepControl};
use super::*;
use crate::compressible::tests::{close, gas};
use fs_qty::Time;

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const INTERIOR: Boundaries3d = Boundaries3d {
    x: OPEN,
    y: OPEN,
    z: OPEN,
};

fn field(
    source: super::super::FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    sample: impl Fn([f64; 3]) -> (f64, [f64; 3], f64),
) -> Euler3d {
    let states = (0..grid.len())
        .map(|cell| {
            let (rho, velocity, pressure) = sample(grid.cell_center(cell).unwrap());
            State::from_primitive(&source, rho, velocity, pressure).unwrap()
        })
        .collect();
    Euler3d::new(source, grid, states, boundaries, || false).unwrap()
}

fn evaluate(flow: &Euler3d) -> Vec<VorticityBudget3d> {
    let bytes = flow.required_vorticity_workspace_bytes().unwrap();
    let mut work = flow.vorticity_workspace(bytes).unwrap();
    flow.vorticity(&mut work, || false).unwrap().to_vec()
}

#[test]
fn g0_native_rotation_expansion_and_boost_on_anisotropic_cells() {
    let grid = CartesianGrid::new(4, 5, 3, 0.3, 0.2, 0.4).unwrap();
    let expansion = 0.07;
    let rotation = 0.3;
    for boost in [[0.0; 3], [1.2, -0.8, 0.3]] {
        let flow = field(gas(0.4), grid, INTERIOR, |x| {
            (
                1.0,
                [
                    expansion * x[0] - rotation * x[1] + boost[0],
                    rotation * x[0] + expansion * x[1] + boost[1],
                    expansion * x[2] + boost[2],
                ],
                2.0,
            )
        });
        for sample in evaluate(&flow) {
            close(sample.divergence, 3.0 * expansion, 2e-12);
            for axis in 0..3 {
                let omega = if axis == 2 { 2.0 * rotation } else { 0.0 };
                close(sample.vorticity[axis], omega, 2e-12);
                close(sample.stretching[axis], expansion * omega, 2e-12);
                close(sample.dilation[axis], -3.0 * expansion * omega, 2e-12);
                close(sample.baroclinic[axis], 0.0, 2e-12);
                close(
                    sample.inviscid_material_rate()[axis],
                    -2.0 * expansion * omega,
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn g0_source_pressure_density_cross_product_all_orientations_and_barotropic_zero() {
    let grid = CartesianGrid::new(4, 3, 5, 0.17, 0.31, 0.23).unwrap();
    for density_axis in 0..3 {
        for pressure_axis in 0..3 {
            // Nonconstant Cp exercises original caloric pressure recovery.
            let source = gas(0.4);
            let flow = field(source.clone(), grid, INTERIOR, |x| {
                (
                    1.0 + 0.1 * x[density_axis],
                    [0.0; 3],
                    2.0 + 0.2 * x[pressure_axis],
                )
            });
            assert_eq!(flow.gas(), &source);
            for (cell, sample) in evaluate(&flow).iter().enumerate() {
                let rho = flow.states()[cell].density;
                let mut expected = [0.0; 3];
                if density_axis != pressure_axis {
                    let axis = 3 - density_axis - pressure_axis;
                    let sign = if (density_axis + 1) % 3 == pressure_axis {
                        1.0
                    } else {
                        -1.0
                    };
                    expected[axis] = sign * 0.02 / (rho * rho);
                }
                assert_eq!(sample.vorticity, [0.0; 3]);
                assert_eq!(sample.stretching, [0.0; 3]);
                for axis in 0..3 {
                    close(sample.baroclinic[axis], expected[axis], 2e-10);
                    close(sample.inviscid_material_rate()[axis], expected[axis], 2e-10);
                }
            }
        }
    }
}

#[test]
fn g1_one_sided_interior_stencils_reproduce_quadratics_at_every_boundary() {
    let grid = CartesianGrid::new(3, 4, 5, 0.13, 0.27, 0.19).unwrap();
    let flow = field(gas(0.4), grid, INTERIOR, |x| {
        (
            1.0,
            [
                0.0,
                0.0,
                0.4 * x[0] * x[0] + 0.2 * x[1] * x[1] - 0.1 * x[2] * x[2],
            ],
            2.0,
        )
    });
    for (cell, sample) in evaluate(&flow).iter().enumerate() {
        let x = grid.cell_center(cell).unwrap();
        close(sample.vorticity[0], 0.4 * x[1], 2e-12);
        close(sample.vorticity[1], -0.8 * x[0], 2e-12);
        close(sample.vorticity[2], 0.0, 2e-12);
        close(sample.divergence, -0.2 * x[2], 2e-12);
    }
}

#[test]
fn g1_periodic_curl_and_baroclinicity_refine_on_physical_axis_widths() {
    let mut errors = Vec::new();
    for n in [12, 24, 48] {
        let grid = CartesianGrid::new(n, n, 1, 2.0 / n as f64, 3.0 / n as f64, 0.7).unwrap();
        let kx = std::f64::consts::PI;
        let ky = 2.0 * std::f64::consts::PI / 3.0;
        let flow = field(gas(0.4), grid, PERIODIC, |x| {
            (
                1.0 + 0.1 * (kx * x[0]).sin(),
                [0.0, 0.3 * (kx * x[0]).sin(), 0.2 * (ky * x[1]).cos()],
                2.0 + 0.2 * (ky * x[1]).sin(),
            )
        });
        let mut error = 0.0;
        for (cell, sample) in evaluate(&flow).iter().enumerate() {
            let x = grid.cell_center(cell).unwrap();
            let rho = flow.states()[cell].density;
            let omega = [
                -0.2 * ky * (ky * x[1]).sin(),
                0.0,
                0.3 * kx * (kx * x[0]).cos(),
            ];
            let baroclinic =
                0.02 * kx * ky * (kx * x[0]).cos() * (ky * x[1]).cos() / (rho * rho);
            error += sample
                .vorticity
                .iter()
                .zip(omega)
                .map(|(actual, exact)| (actual - exact).abs())
                .sum::<f64>();
            error += (sample.baroclinic[2] - baroclinic).abs();
            close(sample.baroclinic[0], 0.0, 2e-12);
            close(sample.baroclinic[1], 0.0, 2e-12);
            close(sample.divergence, 0.0, 2e-12);
        }
        errors.push(error / grid.len() as f64);
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.8,
            "second-order smooth gradients: {errors:?}"
        );
    }
}

#[test]
fn g1_original_euler_step_generates_the_observed_baroclinic_vorticity() {
    let mut errors = Vec::new();
    for n in [12, 24, 48] {
        let dx = 2.0 / n as f64;
        let dy = 1.3 / n as f64;
        let grid = CartesianGrid::new(n, n, 1, dx, dy, 0.7).unwrap();
        let kx = std::f64::consts::PI;
        let ky = 2.0 * std::f64::consts::PI / 1.3;
        let mut flow = field(gas(0.0), grid, PERIODIC, |x| {
            (
                1.0 + 0.12 * (kx * x[0]).sin(),
                [0.0; 3],
                2.0 + 0.15 * (ky * x[1]).sin(),
            )
        });
        let bytes = flow.required_vorticity_workspace_bytes().unwrap();
        let mut work = flow.vorticity_workspace(bytes).unwrap();
        let before = flow.vorticity(&mut work, || false).unwrap().to_vec();
        // At rest D(omega)/Dt = partial_t(omega). A short native RK2 step
        // isolates baroclinic onset; it does not replace the conservative
        // pressure flux by an explicitly added vorticity source.
        let dt = 1e-4 * dx.min(dy);
        let report = flow
            .step(Time::new(dt), StepControl::default(), || false)
            .unwrap();
        close(report.elapsed.value(), dt, 1e-15);
        assert!(!report.first_order_retry);
        let after = flow.vorticity(&mut work, || false).unwrap();
        let mut error = 0.0;
        let mut norm = 0.0;
        let mut correlation = 0.0;
        for (initial, final_state) in before.iter().zip(after) {
            let observed =
                (final_state.vorticity[2] - initial.vorticity[2]) / report.elapsed.value();
            let expected = initial.inviscid_material_rate()[2];
            error += (observed - expected).abs();
            norm += expected.abs();
            correlation += observed * expected;
        }
        assert!(correlation > 0.0, "right-handed baroclinic production");
        errors.push(error / norm);
    }
    // MC reconstruction loses pointwise order at smooth extrema. Check
    // integrated refinement, allowing that known limiter behavior.
    assert!(
        errors[1] < 0.45 * errors[0] && errors[2] < 0.45 * errors[1],
        "{errors:?}"
    );
    assert!(
        errors[2] < 0.025,
        "native onset versus local material budget: {errors:?}"
    );
}

#[test]
fn g4_workspace_budget_cancellation_and_reuse_preserve_accepted_flow() {
    let grid = CartesianGrid::new(4, 3, 1, 0.2, 0.3, 0.7).unwrap();
    let flow = field(gas(0.4), grid, PERIODIC, |_| {
        (1.0, [0.2, -0.1, 0.3], 2.0)
    });
    let bytes = flow.required_vorticity_workspace_bytes().unwrap();
    assert!(flow.vorticity_workspace(bytes - 1).is_err());
    let mut work = flow.vorticity_workspace(bytes).unwrap();
    let capacity = work.workspace_bytes();
    let mut polls = 0;
    let expected = flow
        .vorticity(&mut work, || {
            polls += 1;
            false
        })
        .unwrap()
        .to_vec();
    let accepted = flow.states().to_vec();
    let time = flow.time();
    let mut seen = 0;
    assert_eq!(
        flow.vorticity(&mut work, || {
            seen += 1;
            seen == polls
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), accepted);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.vorticity(&mut work, || false).unwrap(), expected);
    assert_eq!(work.workspace_bytes(), capacity);
    let other_grid = CartesianGrid::new(3, 4, 1, 0.2, 0.3, 0.7).unwrap();
    let other = field(gas(0.4), other_grid, PERIODIC, |_| (1.0, [0.0; 3], 2.0));
    assert!(other.vorticity(&mut work, || false).is_err());
}

#[test]
fn g0_unresolved_axes_refuse_instead_of_manufacturing_zero_curl() {
    for (shape, boundaries) in [([2, 3, 1], PERIODIC), ([3, 3, 1], INTERIOR)] {
        let grid = CartesianGrid::new(shape[0], shape[1], shape[2], 0.2, 0.3, 0.4).unwrap();
        let flow = field(gas(0.4), grid, boundaries, |_| (1.0, [0.0; 3], 2.0));
        assert!(flow.vorticity_workspace(usize::MAX).is_err());
    }
}
