//! G1/G3: actual re-equilibration and transpose of the coupled (u, physical S) map.

use super::*;
use fs_cutfem::elastic3::thermal::{
    ThermoelasticEntropyLimits3, ThermoelasticLinearResponseLimits3,
    ThermoelasticLinearSolveSettings3, ThermoelasticMotion3,
};
use fs_solver::krylov::ResidualClaim;

fn limits(n: usize) -> ThermoelasticLinearResponseLimits3 {
    ThermoelasticLinearResponseLimits3 {
        mechanical: ThermalJacobianLimits3 {
            max_points: POINTS,
            max_workspace_components: 1_000_000,
            max_action_terms: 2_000_000,
        },
        entropy: ThermoelasticEntropyLimits3 {
            max_points: POINTS,
            max_workspace_components: 1_000_000,
            max_action_terms: 2_000_000,
        },
        maximum_native_dofs: 8 * n,
        maximum_workspace_components: 2_000_000,
    }
}

fn solve_settings() -> ThermoelasticLinearSolveSettings3 {
    ThermoelasticLinearSolveSettings3 {
        relative_residual_tolerance: 2e-12,
        maximum_iterations: 30_000,
        poll_iterations: 32,
        maximum_workspace_components: 1_000_000,
    }
}

fn entropy_by_cell(points: &[ThermoelasticResponse3], cells: usize) -> Vec<f64> {
    let mut result = vec![0.0; cells];
    for point in points {
        result[point.mechanical.cell] += point.entropy_j_k;
    }
    result
}

fn motion_slope(position: [f64; 3], _: [f64; 3]) -> [f64; 3] {
    [2e-6, -1e-6, 0.7e-6].map(|value| value * position[0])
}

fn motion(position: [f64; 3], temperature: f64) -> [f64; 3] {
    let slope = motion_slope(position, [0.0; 3]);
    std::array::from_fn(|i| {
        [1e-4, -5e-5, 2e-5][i] * position[0] + slope[i] * (temperature - 300.0)
    })
}

fn owning_cell(bounds: &[HexCell], position: [f64; 3]) -> usize {
    // Native Gauss support points lie inside their owning x/y intervals.
    bounds.iter().position(|cell| (0..3).all(|axis| {
        position[axis] >= cell.lo()[axis] && position[axis] <= cell.hi()[axis]
    })).unwrap()
}

macro_rules! verify_response {
    ($op:ident, $type:ty, $law:ident, $bounds:ident) => {{
        let scales: Vec<_> = (0..$op.cells())
            .map(|i| 0.65 + 0.2 * i as f64 / $op.cells() as f64).collect();
        $op.set_scales(&scales).unwrap();
        let external = $op.body_load(&|_| [0.03, -0.01, 0.02], poll).unwrap();
        let direction: Vec<_> = (0..$op.cells())
            .map(|i| 0.4 - 0.17 * (i % 5) as f64).collect();
        let step = 0.002;
        let solve = |operator: &$type, sign: f64| {
            let eigenstrain = operator.prepare_thermoelastic(POINTS, poll).unwrap();
            let free = operator.eigenstrain_load(&eigenstrain, poll).unwrap();
            let prescribed = |position, _| {
                let cell = owning_cell(&$bounds, position);
                motion(position, temperature($bounds[cell]) + sign * step * direction[cell])
            };
            let supported = operator.prescribed_displacement_load(&prescribed, poll).unwrap();
            let rhs: Vec<_> = external.iter().zip(free).zip(supported)
                .map(|((external, free), supported)| external + free + supported).collect();
            operator.solve_controlled(&rhs, 2e-12, 30_000, 32, |_| poll()).unwrap()
        };
        let base = solve(&$op, 0.0);
        let prescribed = |position, _| {
            motion(position, temperature($bounds[owning_cell(&$bounds, position)]))
        };
        let prepared = $op.prepare_thermoelastic_linear_response(
            base.coefficients(),
            Some(ThermoelasticMotion3 {
                displacement: &prescribed,
                temperature_derivative: Some(&motion_slope),
            }),
            TOTAL_STRAIN, limits($op.n()), poll,
        ).unwrap();
        assert_eq!(prepared.n(), $op.n());
        assert_eq!(prepared.cells(), $op.cells());
        assert!(prepared.cell_potential(0).unwrap().heat_capacity_j_k.abs() > 1e-8);
        let response = $op.solve_thermoelastic_temperature_response(
            &prepared, &direction, solve_settings(), |_| poll(),
        ).unwrap();
        assert!(matches!(response.displacement.residual_claim(),
            ResidualClaim::TrueEuclidean(value) if value < 2e-12));
        assert!(response.displacement.coefficients().iter().any(|value| value.abs() > 1e-6));
        let (_, displacement_cotangent) = fields($op.nodes());
        let entropy_cotangent: Vec<_> = (0..$op.cells())
            .map(|i| 0.7 - 0.025 * i as f64).collect();
        let pullback = $op.solve_thermoelastic_temperature_pullback(
            &prepared, Some(&displacement_cotangent), &entropy_cotangent,
            solve_settings(), |_| poll(),
        ).unwrap();
        assert!(matches!(pullback.adjoint.residual_claim(),
            ResidualClaim::TrueEuclidean(value) if value < 2e-12));
        close(
            dot(response.displacement.coefficients(), &displacement_cotangent)
                + dot(&response.entropy_direction_j_k, &entropy_cotangent),
            dot(&direction, &pullback.temperature_derivative),
            3e-8, 3e-10,
        );
        for (node, position) in $op.nodes().iter().enumerate() {
            if position[0] == 0.0 {
                assert_eq!(&response.displacement.coefficients()[3 * node..3 * node + 3], &[0.0; 3]);
                assert_eq!(&pullback.adjoint.coefficients()[3 * node..3 * node + 3], &[0.0; 3]);
            }
        }
        // Independently update the original source laws AND local embedded
        // motion, solve complete native equilibrium again, and recover physical
        // source entropy. This is not a finite difference of the cached blocks.
        let values = [-1.0, 1.0].map(|sign| {
            let materials = |bounds: HexCell| {
                let cell = $bounds.iter().position(|cell| cell.lo() == bounds.lo()).unwrap();
                Ok(hot(&$law, temperature(bounds) + sign * step * direction[cell]))
            };
            $op.refresh_material_field(
                ElasticMaterialField3::ThermoelasticByCell(&materials),
                MaterialRefreshLimits3::default(), poll,
            ).unwrap();
            let perturbed = solve(&$op, sign);
            let points = $op.thermoelastic_response(
                perturbed.coefficients(), TOTAL_STRAIN, POINTS, poll,
            ).unwrap();
            (perturbed.coefficients().to_vec(), entropy_by_cell(&points, $op.cells()))
        });
        for (i, derivative) in response.displacement.coefficients().iter().enumerate() {
            close(*derivative, (values[1].0[i] - values[0].0[i]) / (2.0 * step),
                1e-5, 3e-9);
        }
        for (i, derivative) in response.entropy_direction_j_k.iter().enumerate() {
            close(*derivative, (values[1].1[i] - values[0].1[i]) / (2.0 * step),
                2e-5, 3e-10);
        }
        assert!($op.solve_thermoelastic_temperature_response(
            &prepared, &direction, solve_settings(), |_| poll(),
        ).is_err());
        assert!($op.solve_thermoelastic_temperature_pullback(
            &prepared, None, &entropy_cotangent, solve_settings(), |_| poll(),
        ).is_err());
    }};
}

#[test]
fn cartesian_response_matches_refreshed_oriented_equilibrium_and_physical_entropy() {
    let law = law(false, false);
    let mut op = cartesian(&law, ElasticityFormulation3::Displacement);
    let bounds: Vec<_> = op.cell_nodes().iter()
        .map(|nodes| HexCell::try_new(op.nodes()[nodes[0]], op.nodes()[nodes[7]]).unwrap())
        .collect();
    verify_response!(op, CutElasticity3, law, bounds);
}

#[test]
fn adaptive_bbar_response_keeps_original_hanging_map_and_temperature_dependent_supports() {
    let law = law(true, true);
    let mut op = adaptive(&law, ElasticityFormulation3::BBar);
    assert!(op.physical_nodes().len() > op.nodes().len());
    let bounds: Vec<_> = op.cell_nodes().iter()
        .map(|nodes| HexCell::try_new(op.physical_nodes()[nodes[0]], op.physical_nodes()[nodes[7]]).unwrap())
        .collect();
    let (u, _) = fields(op.nodes());
    assert!(op.prepare_thermoelastic_linear_response(
        &u, None, TOTAL_STRAIN,
        ThermoelasticLinearResponseLimits3 { maximum_native_dofs: op.n(), ..limits(op.n()) },
        poll,
    ).is_err());
    verify_response!(op, AdaptiveElasticity3, law, bounds);
}

#[test]
fn linear_response_bounds_and_late_cancellation_preserve_the_original_owner() {
    let law = law(true, true);
    let mut op = cartesian(&law, ElasticityFormulation3::BBar);
    let (u, _) = fields(op.nodes());
    let mut initial_action = vec![0.0; op.n()];
    op.apply(&u, &mut initial_action);
    let original = op.prepare_thermoelastic(POINTS, poll).unwrap();
    for denied in [
        ThermoelasticLinearResponseLimits3 {
            maximum_native_dofs: op.n() - 1, ..limits(op.n())
        },
        ThermoelasticLinearResponseLimits3 {
            maximum_workspace_components: 0, ..limits(op.n())
        },
    ] {
        assert!(op.prepare_thermoelastic_linear_response(
            &u, None, TOTAL_STRAIN, denied, poll,
        ).is_err());
    }
    let prepared = op.prepare_thermoelastic_linear_response(
        &u, None, TOTAL_STRAIN, limits(op.n()), poll,
    ).unwrap();
    let direction = vec![1.0; op.cells()];
    let mut calls = 0;
    let response = op.solve_thermoelastic_temperature_response(
        &prepared, &direction, solve_settings(), |_| { calls += 1; poll() },
    ).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut at = 0;
        assert!(matches!(op.solve_thermoelastic_temperature_response(
            &prepared, &direction, solve_settings(), |_| {
                at += 1;
                if at == stop { ControlFlow::Break(()) } else { poll() }
            },
        ), Err(ElasticityError3::Cancelled)));
    }
    assert!(op.solve_thermoelastic_temperature_response(
        &prepared, &direction,
        ThermoelasticLinearSolveSettings3 { maximum_workspace_components: 0, ..solve_settings() },
        |_| poll(),
    ).is_err());
    assert!(matches!(op.solve_thermoelastic_temperature_response(
        &prepared, &direction,
        ThermoelasticLinearSolveSettings3 { maximum_iterations: 0, ..solve_settings() },
        |_| poll(),
    ), Err(ElasticityError3::NotConverged { .. })));
    let mut invalid = direction.clone();
    invalid[0] = f64::NAN;
    assert!(op.solve_thermoelastic_temperature_response(
        &prepared, &invalid, solve_settings(), |_| poll(),
    ).is_err());
    assert!(op.solve_thermoelastic_temperature_pullback(
        &prepared, None, &direction[..direction.len() - 1], solve_settings(), |_| poll(),
    ).is_err());

    // Reuse treats a response as an algebraic direction, without imposing a
    // physical strain ceiling on its chosen amplitude.
    let large = vec![1e6; op.cells()];
    let algebraic = op.solve_thermoelastic_temperature_response(
        &prepared, &large, solve_settings(), |_| poll(),
    ).unwrap();
    assert!(algebraic.displacement.coefficients().iter().any(|value| value.abs() > 1.0));
    for (actual, base) in algebraic.displacement.coefficients().iter()
        .zip(response.displacement.coefficients())
    {
        close(*actual, 1e6 * base, 3e-8, 3e-6);
    }
    op.validate_eigenstrain(&original).unwrap();
    let mut after = vec![0.0; op.n()];
    op.apply(&u, &mut after);
    assert_eq!(after, initial_action);

    let other = cartesian(&law, ElasticityFormulation3::BBar);
    assert!(other.solve_thermoelastic_temperature_response(
        &prepared, &direction, solve_settings(), |_| poll(),
    ).is_err());
    op.set_scales(&vec![0.7; op.cells()]).unwrap();
    assert!(op.solve_thermoelastic_temperature_response(
        &prepared, &direction, solve_settings(), |_| poll(),
    ).is_err());
    assert!(op.solve_thermoelastic_temperature_pullback(
        &prepared, None, &direction, solve_settings(), |_| poll(),
    ).is_err());
}
