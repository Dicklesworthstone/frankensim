//! Weak applied loads must survive balanced resolved fluid tractions.
use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, PartialSaturationBody3, Q3};

const N: usize = 4;
const CELLS: [[usize; 3]; 4] = [[1, 1, 1], [2, 1, 1], [1, 2, 1], [2, 2, 1]];

fn body() -> PartialSaturationBody3 {
    PartialSaturationBody3 {
        center: [2.0, 2.0, 1.5],
        mass: 2.0,
        inertia_world: [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 2.0]],
        velocity: [0.0; 3],
        angular_velocity: [0.0; 3],
    }
}

fn field(fraction: f64) -> Vec<f64> {
    let mut result = vec![0.0; N * N * N];
    for [x, y, z] in CELLS {
        result[(z * N + y) * N + x] = fraction;
    }
    result
}

fn grid() -> BoundaryGrid3 {
    let mut result = BoundaryGrid3::new(N, N, N, 1.0, [0.0; 3], BoundarySpec3::periodic());
    for [x, y, z] in CELLS {
        // Exactly binary populations give rho=1 and ux=+/-1/64 without
        // introducing an equilibrium-rounding residual into the load oracle.
        let mut populations = [1.0 / 32.0; Q3];
        populations[0] = 7.0 / 16.0;
        let shift = if x == 1 { 1.0 / 128.0 } else { -1.0 / 128.0 };
        populations[1] += shift;
        populations[2] -= shift;
        result.set_populations(x, y, z, &populations);
    }
    result
}

fn close(actual: f64, expected: f64) {
    let tolerance = 16.0 * f64::EPSILON * expected.abs();
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn assert_response(body: PartialSaturationBody3, external: [f64; 6], share: f64) {
    // Four equal cells at r=(+/-1/2,+/-1/2,0): sum r=0,
    // Ktt=4B I, Krr=B diag(1,1,2). All resolved incoming forces
    // and couples cancel, so the remaining RHS is exactly the external load.
    let linear_inertia = 2.0 + 2.0 * share;
    let rotary_inertia = [3.0 + 0.5 * share, 3.0 + 0.5 * share, 2.0 + share];
    for axis in 0..3 {
        close(body.velocity[axis], external[axis] / linear_inertia);
        close(
            body.angular_velocity[axis],
            external[axis + 3] / rotary_inertia[axis],
        );
    }
}

fn loads() -> [[f64; 6]; 2] {
    let weak = 2.0_f64.powi(-62);
    [
        [weak, -2.0 * weak, 3.0 * weak, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 5.0 * weak, -7.0 * weak, 11.0 * weak],
    ]
}

#[test]
fn single_body_retains_small_force_and_pure_couple_under_opposing_fluid_loads() {
    for external in loads() {
        let mut flow = grid();
        let mut solid = body();
        let report = flow
            .step_with_partial_saturation_body(
                &field(1.0),
                &mut solid,
                external,
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        assert_response(solid, external, 1.0);
        for axis in 0..3 {
            assert_eq!(
                report.hydrodynamic.body_impulse[axis],
                -report.hydrodynamic.fluid_population_impulse[axis]
            );
            assert_eq!(
                report.hydrodynamic.body_angular_impulse[axis],
                -report.hydrodynamic.fluid_population_angular_impulse[axis]
            );
        }
        // Sub-ulp population reactions cannot be manufactured by a precise
        // body solve. Keep the actual represented mismatch visible to callers.
        assert!(
            report
                .momentum_residual
                .iter()
                .all(|value| value.is_finite())
        );
        assert!(report.momentum_residual.iter().any(|&value| value != 0.0));
    }
}

#[test]
fn simultaneous_bodies_use_the_same_small_force_and_couple_response() {
    for first in loads() {
        let external = [first, first.map(|value| -3.0 * value)];
        let fractions = field(0.5);
        let mut flow = grid();
        let mut solids = [body(), body()];
        let report = flow
            .step_with_partial_saturation_bodies(
                [&fractions, &fractions],
                &mut solids,
                external,
                [0.0; 3],
                &CancelGate::new(),
            )
            .unwrap();
        for a in 0..2 {
            // Total occupancy is one, so each member has B=1/2 at tau=1.
            assert_response(solids[a], external[a], 0.5);
            assert_eq!(report.bodies[a].hydrodynamic.coupled_cells, 4);
            for axis in 0..3 {
                assert_eq!(
                    report.bodies[a].hydrodynamic.body_impulse[axis],
                    -report.bodies[a].hydrodynamic.fluid_population_impulse[axis]
                );
                assert_eq!(
                    report.bodies[a].hydrodynamic.body_angular_impulse[axis],
                    -report.bodies[a]
                        .hydrodynamic
                        .fluid_population_angular_impulse[axis]
                );
            }
            assert!(
                report.bodies[a]
                    .momentum_residual
                    .iter()
                    .all(|value| value.is_finite())
            );
        }
    }
}
