//! Source-potential Maxwell reciprocity and reversible mechanical heat exchange.

use super::*;
use fs_cutfem::elastic3::thermal::ThermoelasticEntropyLimits3;

fn limits() -> ThermoelasticEntropyLimits3 {
    ThermoelasticEntropyLimits3 {
        max_points: POINTS,
        max_workspace_components: 1_000_000,
        max_action_terms: 2_000_000,
    }
}

fn cell_entropy(points: &[ThermoelasticResponse3], cells: usize) -> Vec<f64> {
    let mut result = vec![0.0; cells];
    for point in points {
        result[point.mechanical.cell] += point.entropy_j_k;
    }
    result
}

macro_rules! verify_native_entropy {
    ($op:ident, $law:ident) => {{
        let scales: Vec<_> = (0..$op.cells())
            .map(|i| 0.6 + 0.3 * i as f64 / $op.cells() as f64).collect();
        $op.set_scales(&scales).unwrap();
        let (u, direction) = fields($op.nodes());
        let response = $op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll).unwrap();
        let prepared = $op.prepare_thermoelastic_entropy(&u, TOTAL_STRAIN, limits(), poll).unwrap();
        assert_eq!(prepared.n(), $op.n());
        assert_eq!(prepared.cells(), $op.cells());
        assert_eq!(prepared.points(), response.len());
        assert_eq!(prepared.workspace_components(), 30 * $op.cells());
        let mut totals = vec![[0.0; 4]; $op.cells()];
        for point in &response {
            for (sum, value) in totals[point.mechanical.cell].iter_mut().zip([
                point.free_energy_j, point.entropy_j_k,
                point.internal_energy_j, point.heat_capacity_j_k,
            ]) {
                *sum += value;
            }
        }
        for (id, expected) in totals.iter().enumerate() {
            let potential = prepared.cell_potential(id).unwrap();
            for (value, expected) in [
                potential.free_energy_j, potential.entropy_j_k,
                potential.internal_energy_j, potential.heat_capacity_j_k,
            ].into_iter().zip(expected) {
                close(value, *expected, 2e-13, 1e-13);
            }
            close(
                potential.internal_energy_j,
                potential.free_energy_j + potential.temperature_k * potential.entropy_j_k,
                2e-13, 1e-13,
            );
        }

        let mut entropy_direction = vec![0.0; $op.cells()];
        $op.thermoelastic_entropy_displacement_apply_into(
            &prepared, &direction, &mut entropy_direction, poll,
        ).unwrap();
        // The algebraic direction is deliberately much larger than the source
        // domain. Only each perturbed PHYSICAL state consumes the strain budget.
        let step = 1e-6;
        let values = [-1.0, 1.0].map(|sign| {
            let displaced: Vec<_> = u.iter().zip(&direction)
                .map(|(value, direction)| value + sign * step * direction).collect();
            cell_entropy(
                &$op.thermoelastic_response(&displaced, TOTAL_STRAIN, POINTS, poll).unwrap(),
                $op.cells(),
            )
        });
        for id in 0..$op.cells() {
            close(entropy_direction[id], (values[1][id] - values[0][id]) / (2.0 * step),
                3e-8, 3e-10);
        }
        let cotangent: Vec<_> = (0..$op.cells())
            .map(|i| 0.7 - 0.03 * i as f64).collect();
        let mut dual = vec![0.0; $op.n()];
        $op.thermoelastic_entropy_displacement_transpose_into(
            &prepared, &cotangent, &mut dual, poll,
        ).unwrap();
        close(dot(&entropy_direction, &cotangent), dot(&direction, &dual), 2e-13, 2e-12);
        for (node, position) in $op.nodes().iter().enumerate() {
            if position[0] == 0.0 {
                assert_eq!(&dual[3 * node..3 * node + 3], &[0.0; 3]);
            }
        }
        let mut heat = vec![0.0; $op.cells()];
        $op.thermoelastic_reversible_heat_rate_into(&prepared, &direction, &mut heat, poll).unwrap();
        for id in 0..$op.cells() {
            close(heat[id], -prepared.cell_potential(id).unwrap().temperature_k * entropy_direction[id],
                2e-13, 2e-12);
        }

        // The same cell capacity differentiates entropy at fixed TOTAL strain,
        // including changing C(T), integrated oriented expansion and B-bar.
        let step = 0.002;
        let values = [-1.0, 1.0].map(|sign| {
            let materials = |cell: HexCell| Ok(hot(&$law, temperature(cell) + sign * step));
            $op.refresh_material_field(
                ElasticMaterialField3::ThermoelasticByCell(&materials),
                MaterialRefreshLimits3::default(), poll,
            ).unwrap();
            cell_entropy(
                &$op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll).unwrap(),
                $op.cells(),
            )
        });
        for id in 0..$op.cells() {
            let potential = prepared.cell_potential(id).unwrap();
            close(
                potential.heat_capacity_j_k / potential.temperature_k,
                (values[1][id] - values[0][id]) / (2.0 * step),
                2e-6, 2e-11,
            );
        }
        // A published material change invalidates the old reusable block.
        let mut untouched = vec![97.0; $op.cells()];
        assert!($op.thermoelastic_entropy_displacement_apply_into(
            &prepared, &direction, &mut untouched, poll,
        ).is_err());
        assert!(untouched.iter().all(|value| *value == 97.0));
    }};
}

#[test]
fn physical_entropy_and_capacity_follow_native_cartesian_sources() {
    for (isotropic, formulation) in [
        (false, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::BBar),
    ] {
        let law = law(isotropic, formulation == ElasticityFormulation3::BBar);
        let mut op = cartesian(&law, formulation);
        verify_native_entropy!(op, law);
    }
}

#[test]
fn adaptive_entropy_reuses_hanging_reconstruction_and_exact_transpose() {
    for (isotropic, formulation) in [
        (false, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::BBar),
    ] {
        let law = law(isotropic, formulation == ElasticityFormulation3::BBar);
        let mut op = adaptive(&law, formulation);
        verify_native_entropy!(op, law);
    }
}

#[test]
fn bulk_maxwell_reciprocity_matches_unstabilized_mechanical_temperature_block() {
    for (isotropic, formulation) in [
        (false, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::BBar),
    ] {
        let law = law(isotropic, formulation == ElasticityFormulation3::BBar);
        let materials = |cell: HexCell| Ok(hot(&law, temperature(cell)));
        let op = control(|control| CutElasticity3::build_with_material_field(
            domain(), [2; 3], &Top,
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|position| position[0] == 0.0,
            ElasticityOptions3 { formulation, ghost_gamma: 0.0, ..Default::default() },
            control,
        ).unwrap());
        let (u, _) = fields(op.nodes());
        let entropy = op.prepare_thermoelastic_entropy(&u, TOTAL_STRAIN, limits(), poll).unwrap();
        let residual = op.prepare_thermal_residual_jacobian(
            &u, TOTAL_STRAIN,
            ThermalJacobianLimits3 {
                max_points: POINTS,
                max_workspace_components: 1_000_000,
                max_action_terms: 2_000_000,
            }, poll,
        ).unwrap();
        let temperature_direction: Vec<_> = (0..op.cells())
            .map(|i| 1.0 - 0.1 * i as f64).collect();
        let mut physical_transpose = vec![0.0; op.n()];
        let mut mechanical = vec![0.0; op.n()];
        op.thermoelastic_entropy_displacement_transpose_into(
            &entropy, &temperature_direction, &mut physical_transpose, poll,
        ).unwrap();
        op.thermal_residual_temperature_apply_into(
            &residual, &temperature_direction, &mut mechanical, poll,
        ).unwrap();
        for (entropy, force) in physical_transpose.iter().zip(&mechanical) {
            close(*entropy, -*force, 3e-11, 3e-11);
        }
    }
}

#[test]
fn physical_entropy_excludes_nitsche_and_ghost_thermal_work() {
    let law = law(false, false);
    let op = cartesian(&law, ElasticityFormulation3::Displacement);
    let (u, direction) = fields(op.nodes());
    let entropy = op.prepare_thermoelastic_entropy(&u, TOTAL_STRAIN, limits(), poll).unwrap();
    let residual = op.prepare_thermal_residual_jacobian(
        &u, TOTAL_STRAIN,
        ThermalJacobianLimits3 {
            max_points: POINTS,
            max_workspace_components: 1_000_000,
            max_action_terms: 2_000_000,
        }, poll,
    ).unwrap();
    let mut entropy_direction = vec![0.0; op.cells()];
    let mut full_residual_work = vec![0.0; op.cells()];
    op.thermoelastic_entropy_displacement_apply_into(
        &entropy, &direction, &mut entropy_direction, poll,
    ).unwrap();
    op.thermal_residual_temperature_transpose_into(
        &residual, &direction, &mut full_residual_work, poll,
    ).unwrap();
    let difference: f64 = entropy_direction.iter().zip(&full_residual_work)
        .map(|(entropy, residual)| (entropy + residual).abs()).sum();
    assert!(difference > 1e-6, "fixture must exercise nonphysical stabilization work");
    let step = 1e-6;
    let actual = [-1.0, 1.0].map(|sign| {
        let displaced: Vec<_> = u.iter().zip(&direction)
            .map(|(value, direction)| value + sign * step * direction).collect();
        cell_entropy(
            &op.thermoelastic_response(&displaced, TOTAL_STRAIN, POINTS, poll).unwrap(),
            op.cells(),
        )
    });
    for id in 0..op.cells() {
        close(entropy_direction[id], (actual[1][id] - actual[0][id]) / (2.0 * step),
            3e-8, 3e-10);
    }
}

#[test]
fn reversible_dilatation_cools_and_elastic_capacity_keeps_its_sign() {
    let law = law(true, false);
    let material = hot(&law, 300.0);
    let op = control(|control| CutElasticity3::build_with_material_field(
        domain(), [1; 3], &Top,
        ElasticMaterialField3::ThermoelasticUniform(&material),
        &|position| position == [0.0; 3],
        ElasticityOptions3 { ghost_gamma: 0.0, ..Default::default() },
        control,
    ).unwrap());
    let u = vec![0.0; op.n()];
    let velocity: Vec<_> = op.nodes().iter().flat_map(|position| *position).collect();
    let entropy = op.prepare_thermoelastic_entropy(&u, TOTAL_STRAIN, limits(), poll).unwrap();
    let mut heat = vec![0.0; op.cells()];
    op.thermoelastic_reversible_heat_rate_into(&entropy, &velocity, &mut heat, poll).unwrap();
    let point = material.source().state([0.0; 6], 1.0).unwrap();
    let volume: f64 = op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll).unwrap()
        .iter().map(|response| response.mechanical.weight).sum();
    let expected = 300.0 * volume
        * point.stress_temperature_derivative_mandel_pa_k[..3].iter().sum::<f64>();
    close(heat.iter().sum(), expected, 3e-13, 3e-12);
    assert!(heat.iter().sum::<f64>() < 0.0);
    assert!(entropy.cell_potential(0).unwrap().heat_capacity_j_k < 0.0);
    assert_eq!(entropy.cell_potential(0).unwrap().free_energy_j, 0.0);
}

#[test]
fn entropy_budgets_stale_scales_and_invalid_actions_refuse_before_writes() {
    let law = law(true, true);
    let mut op = cartesian(&law, ElasticityFormulation3::BBar);
    let (u, direction) = fields(op.nodes());
    for limits in [
        ThermoelasticEntropyLimits3 { max_points: 0, ..limits() },
        ThermoelasticEntropyLimits3 { max_workspace_components: 0, ..limits() },
        ThermoelasticEntropyLimits3 { max_action_terms: 0, ..limits() },
    ] {
        assert!(op.prepare_thermoelastic_entropy(&u, TOTAL_STRAIN, limits, poll).is_err());
    }
    let prepared = op.prepare_thermoelastic_entropy(&u, TOTAL_STRAIN, limits(), poll).unwrap();
    let mut output = vec![97.0; op.cells()];
    assert!(op.thermoelastic_entropy_displacement_apply_into(
        &prepared, &direction[..direction.len() - 1], &mut output, poll,
    ).is_err());
    assert!(output.iter().all(|value| *value == 97.0));
    let mut invalid = direction.clone();
    invalid[0] = f64::NAN;
    assert!(op.thermoelastic_entropy_displacement_apply_into(
        &prepared, &invalid, &mut output, poll,
    ).is_err());
    assert!(output.iter().all(|value| *value == 97.0));
    assert!(matches!(
        op.thermoelastic_entropy_displacement_apply_into(
            &prepared, &direction, &mut output, || ControlFlow::Break(()),
        ),
        Err(ElasticityError3::Cancelled),
    ));
    assert!(output.iter().all(|value| *value == 97.0));
    op.set_scales(&vec![0.7; op.cells()]).unwrap();
    assert!(op.thermoelastic_entropy_displacement_apply_into(
        &prepared, &direction, &mut output, poll,
    ).is_err());
    assert!(output.iter().all(|value| *value == 97.0));
}
