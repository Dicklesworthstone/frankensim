//! Finite source-owned physical exchange on native Cartesian/adaptive fields.

use super::*;
use fs_cutfem::elastic3::thermal::{
    ThermoelasticCellExchange3, ThermoelasticExchangeLimits3,
};
use fs_material::elastic::thermal::ThermoelasticExchange;
use fs_material::tensor::{contract, frobenius_norm, to_mandel};
use fs_math::dd::Dd;

fn limits(cells: usize, points: usize) -> ThermoelasticExchangeLimits3 {
    ThermoelasticExchangeLimits3 {
        maximum_cells: cells,
        maximum_points: points,
        maximum_workspace_bytes: cells * core::mem::size_of::<ThermoelasticCellExchange3>(),
    }
}

fn source_accounts(value: &ThermoelasticExchange) -> [Dd; 7] {
    [
        value.retained_mechanical_work_j(),
        value.retained_free_energy_change_j(),
        value.retained_entropy_change_j_k(),
        value.retained_internal_energy_change_j(),
        value.retained_reversible_heat_j(),
        value.retained_expansion_work_j(),
        value.retained_stiffness_work_j(),
    ]
}

fn native_accounts(value: &ThermoelasticCellExchange3) -> [Dd; 7] {
    [
        value.retained_mechanical_work_j(),
        value.retained_free_energy_change_j(),
        value.retained_entropy_change_j_k(),
        value.retained_internal_energy_change_j(),
        value.retained_reversible_heat_j(),
        value.retained_expansion_work_j(),
        value.retained_stiffness_work_j(),
    ]
}

fn reported_accounts(value: &ThermoelasticCellExchange3) -> [f64; 7] {
    [
        value.mechanical_work_j,
        value.free_energy_change_j,
        value.entropy_change_j_k,
        value.internal_energy_change_j,
        value.reversible_heat_j,
        value.expansion_work_j,
        value.stiffness_work_j,
    ]
}

fn total(point: &ThermoelasticResponse3) -> [f64; 6] {
    // Recovery publishes physical tensor components. The source potential
    // consumes Mandel total strain, including B-bar's constitutive projection.
    std::array::from_fn(|i| {
        point.mechanical.constitutive_elastic_strain[i] + point.mechanical.free_strain[i]
    })
}

fn independent_accounts(
    law: &AffineThermoelastic,
    before: &[ThermoelasticResponse3],
    after: &[ThermoelasticResponse3],
    scales: &[f64],
) -> Vec<[Dd; 7]> {
    assert_eq!(before.len(), after.len());
    let mut sums = vec![[Dd::ZERO; 7]; scales.len()];
    for (before, after) in before.iter().zip(after) {
        let p0 = &before.mechanical;
        let p1 = &after.mechanical;
        assert_eq!((p0.cell, p0.position, p0.weight), (p1.cell, p1.position, p1.weight));
        let interval = law.exchange(
            before.temperature_k, to_mandel(&total(before)),
            after.temperature_k, to_mandel(&total(after)), p0.weight,
        ).unwrap();
        close(scales[p0.cell] * interval.before.free_energy_j, before.free_energy_j, 2e-12, 2e-14);
        close(scales[p0.cell] * interval.after.free_energy_j, after.free_energy_j, 2e-12, 2e-14);
        for (sum, increment) in sums[p0.cell].iter_mut().zip(source_accounts(&interval)) {
            *sum = *sum + increment;
        }
    }
    for (sum, scale) in sums.iter_mut().zip(scales) {
        for value in sum {
            *value = *value * Dd::from_f64(*scale);
        }
    }
    sums
}

fn independent_work(
    before: &[ThermoelasticResponse3], after: &[ThermoelasticResponse3], cells: usize,
) -> Vec<f64> {
    let mut sums = vec![0.0; cells];
    for (before, after) in before.iter().zip(after) {
        let mean_stress = std::array::from_fn(|i| {
            0.5 * (before.mechanical.stress[i] + after.mechanical.stress[i])
        });
        let e0 = total(before);
        let e1 = total(after);
        let increment = std::array::from_fn(|i| e1[i] - e0[i]);
        sums[before.mechanical.cell] += before.mechanical.weight * contract(&mean_stress, &increment);
    }
    sums
}

macro_rules! native_exchange {
    ($before:ident, $after:ident, $law:ident, $projected:expr) => {{
        let cells = $after.cells();
        assert_eq!($before.nodes(), $after.nodes());
        assert!($after.has_embedded_supports());
        let bounds: Vec<_> = (0..cells).map(|cell| $after.reference_cell_bounds(cell).unwrap()).collect();
        assert!($after.reference_cell_bounds(cells).is_none());
        assert!($after.reference_cell_volume_m3(cells).is_none());
        let scales: Vec<_> = (0..cells).map(|i| 0.55 + 0.4 * i as f64 / cells as f64).collect();
        $before.set_scales(&scales).unwrap();
        $after.set_scales(&scales).unwrap();
        let temperatures0: Vec<_> = (0..cells).map(|i| 303.0 + 0.3 * i as f64).collect();
        let materials0 = |cell: HexCell| {
            let id = bounds.iter().position(|candidate| candidate.lo() == cell.lo()).unwrap();
            Ok(hot(&$law, temperatures0[id]))
        };
        $before.refresh_material_field(
            ElasticMaterialField3::ThermoelasticByCell(&materials0),
            MaterialRefreshLimits3::default(), poll,
        ).unwrap();
        let laws = vec![$law.clone(); cells];
        let (base, direction) = fields($after.nodes());
        let mut u0: Vec<_> = base.iter().zip(&direction).map(|(u, d)| u - 0.00009 * d).collect();
        let mut u1: Vec<_> = base.iter().zip(&direction).map(|(u, d)| u + 0.00017 * d).collect();
        // Finite coefficient entries at homogeneous supports must be masked,
        // including when their unused values are outside the physical domain.
        for (node, fixed) in $after.fixed().iter().enumerate() {
            if *fixed {
                for axis in 0..3 {
                    u0[3 * node + axis] = 17.0 + axis as f64;
                    u1[3 * node + axis] = -29.0 - axis as f64;
                }
            }
        }
        let p0 = $before.thermoelastic_response(&u0, TOTAL_STRAIN, POINTS, poll).unwrap();
        let p1 = $after.thermoelastic_response(&u1, TOTAL_STRAIN, POINTS, poll).unwrap();
        assert!(p1.iter().any(|p| p.temperature_k != p1[0].temperature_k));
        assert!(p1.iter().any(|p| p.mechanical.free_strain[5].abs() > 1e-5));
        if $projected {
            assert!(p1.iter().any(|p| p.mechanical.elastic_strain.iter()
                .zip(p.mechanical.constitutive_elastic_strain)
                .any(|(a, b)| (a - b).abs() > 1e-8)));
        }
        let expected = independent_accounts(&$law, &p0, &p1, &scales);
        let work = independent_work(&p0, &p1, cells);
        let exchange = $after.thermoelastic_exchange(
            &laws, &temperatures0, &u0, &u1, TOTAL_STRAIN, limits(cells, p1.len()), poll,
        ).unwrap();
        assert_eq!(exchange.len(), cells);
        let mut volume = vec![0.0; cells];
        let mut endpoint_changes = vec![[0.0; 3]; cells];
        for (before, after) in p0.iter().zip(&p1) {
            let id = after.mechanical.cell;
            volume[id] += after.mechanical.weight;
            for (sum, change) in endpoint_changes[id].iter_mut().zip([
                after.free_energy_j - before.free_energy_j,
                after.entropy_j_k - before.entropy_j_k,
                after.internal_energy_j - before.internal_energy_j,
            ]) {
                *sum += change;
            }
        }
        for (id, exchange) in exchange.iter().enumerate() {
            assert_eq!(exchange.cell, id);
            assert_eq!(exchange.stiffness_scale, scales[id]);
            assert_eq!(exchange.temperature_k, [temperatures0[id], temperature(bounds[id])]);
            close(exchange.reference_volume_m3, volume[id], 2e-13, 2e-14);
            assert_eq!(exchange.reference_volume_m3, $after.reference_cell_volume_m3(id).unwrap());
            for ((actual, reference), rounded) in native_accounts(exchange).into_iter()
                .zip(expected[id]).zip(reported_accounts(exchange))
            {
                // Reconstructing published total strain can round; the native
                // integration keeps the original unrounded kinematic owner.
                close(actual.to_f64(), reference.to_f64(), 3e-10, 3e-13);
                assert_eq!(rounded, actual.to_f64());
            }
            close(exchange.mechanical_work_j, work[id], 3e-10, 3e-13);
            for (actual, expected) in [
                exchange.free_energy_change_j,
                exchange.entropy_change_j_k,
                exchange.internal_energy_change_j,
            ].into_iter().zip(endpoint_changes[id]) {
                close(actual, expected, 3e-10, 3e-13);
            }
            close(
                (exchange.retained_mechanical_work_j() - exchange.retained_internal_energy_change_j())
                    .to_f64(),
                exchange.reversible_heat_j, 3e-14, 3e-15,
            );
        }
        assert_eq!($after.scales(), scales.as_slice());

        // At fixed temperature, linear-elastic trapezoidal stress work equals
        // the free-energy increment. Numerical ghost/Nitsche work is excluded.
        let temperatures1: Vec<_> = bounds.iter().copied().map(temperature).collect();
        let fixed0 = $after.thermoelastic_response(&u0, TOTAL_STRAIN, POINTS, poll).unwrap();
        let fixed = $after.thermoelastic_exchange(
            &laws, &temperatures1, &u0, &u1, TOTAL_STRAIN, limits(cells, p1.len()), poll,
        ).unwrap();
        let fixed_work = independent_work(&fixed0, &p1, cells);
        for (id, value) in fixed.iter().enumerate() {
            close(value.mechanical_work_j, fixed_work[id], 3e-10, 3e-13);
            close(value.mechanical_work_j, value.free_energy_change_j, 3e-13, 3e-14);
            close(value.reversible_heat_j, -temperatures1[id] * value.entropy_change_j_k, 3e-13, 3e-14);
            assert_eq!(value.expansion_work_j, 0.0);
            assert_eq!(value.stiffness_work_j, 0.0);
        }
    }};
}

#[test]
fn oriented_native_exchange_matches_original_source_and_physical_stress_work() {
    let law = law(false, false);
    let mut before = cartesian(&law, ElasticityFormulation3::Displacement);
    let mut after = cartesian(&law, ElasticityFormulation3::Displacement);
    native_exchange!(before, after, law, false);
}

#[test]
fn adaptive_bbar_exchange_uses_original_hanging_space_projection_and_fixed_scales() {
    let law = law(true, true);
    let mut before = adaptive(&law, ElasticityFormulation3::BBar);
    let mut after = adaptive(&law, ElasticityFormulation3::BBar);
    assert!(after.physical_nodes().len() > after.nodes().len());
    let nodes = after.physical_nodes().to_vec();
    native_exchange!(before, after, law, true);
    assert_eq!(after.physical_nodes(), nodes.as_slice());
}

#[test]
fn prepared_source_intervals_retain_tiny_work_when_rounded_endpoint_energy_is_identical() {
    for isotropic in [false, true] {
        let law = law(isotropic, false);
        let prepared = law.prepare_exchange(304.0, 329.0).unwrap();
        assert_eq!(prepared.before_source(), &law.at_temperature(304.0).unwrap());
        assert_eq!(prepared.after_source(), &law.at_temperature(329.0).unwrap());
        for volume in [0.03125, 0.17] {
            let e0 = [0.003, -0.001, 0.002, 0.0003, -0.0002, 0.0004];
            let e1 = [0.002, 0.0002, -0.001, -0.0001, 0.0004, 0.0002];
            assert_eq!(
                prepared.evaluate(e0, e1, volume).unwrap(),
                law.exchange(304.0, e0, 329.0, e1, volume).unwrap(),
            );
        }

        let prepared = law.prepare_exchange(319.0, 319.0).unwrap();
        let start: [f64; 6] = [1e-20, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut current = start;
        let mut sums = [Dd::ZERO; 7];
        for _ in 0..64 {
            let mut next = current;
            next[0] = next[0].next_up();
            let interval = prepared.evaluate(current, next, 0.3).unwrap();
            for (sum, increment) in sums.iter_mut().zip(source_accounts(&interval)) {
                *sum = *sum + increment;
            }
            current = next;
        }
        let entire = prepared.evaluate(start, current, 0.3).unwrap();
        assert_eq!(entire.before.elastic_strain_mandel, entire.after.elastic_strain_mandel);
        assert_eq!(entire.before.free_energy_j, entire.after.free_energy_j);
        assert_ne!(entire.free_energy_change_j, 0.0);
        assert_ne!(entire.mechanical_work_j, 0.0);
        for (sum, whole) in sums.into_iter().zip(source_accounts(&entire)) {
            close(sum.to_f64(), whole.to_f64(), 3e-13, 1e-47);
        }
        assert_eq!(
            entire, law.exchange(319.0, start, 319.0, current, 0.3).unwrap(),
        );
    }
}

#[test]
fn bbar_exchange_checks_actual_total_and_elastic_boxes_before_projected_energy() {
    let law = law(true, true);
    let op = cartesian(&law, ElasticityFormulation3::BBar);
    let laws = vec![law; op.cells()];
    let temperatures: Vec<_> = (0..op.cells())
        .map(|cell| temperature(op.reference_cell_bounds(cell).unwrap())).collect();
    for sign in [1.0, -1.0] {
        let safe: Vec<_> = op.nodes().iter().flat_map(|&[x, y, _]| {
            [sign * 0.001 * x * y, 0.0, 0.0]
        }).collect();
        let response = op.thermoelastic_response(&safe, TOTAL_STRAIN, POINTS, poll).unwrap();
        let mut actual_total_limit = f64::INFINITY;
        let mut actual_elastic_limit = f64::INFINITY;
        let mut projected_limit = f64::INFINITY;
        for point in &response {
            let actual = to_mandel(&point.mechanical.strain);
            let assumed = to_mandel(&total(point));
            let free = to_mandel(&point.mechanical.free_strain);
            for i in 0..6 {
                if actual[i] != 0.0 {
                    actual_total_limit = actual_total_limit.min(0.05 / actual[i].abs());
                    actual_elastic_limit = actual_elastic_limit.min(
                        (0.05 + actual[i].signum() * free[i]) / actual[i].abs(),
                    );
                }
                if assumed[i] != 0.0 {
                    projected_limit = projected_limit.min(0.05 / assumed[i].abs())
                        .min((0.05 + assumed[i].signum() * free[i]) / assumed[i].abs());
                }
            }
            projected_limit = projected_limit
                .min(TOTAL_STRAIN / frobenius_norm(&point.mechanical.strain))
                .min(TOTAL_STRAIN / frobenius_norm(&total(point)));
        }
        let forbidden = if sign > 0.0 { actual_total_limit } else { actual_elastic_limit };
        let allowed_projected = if sign > 0.0 {
            projected_limit
        } else {
            projected_limit.min(actual_total_limit)
        };
        assert!(allowed_projected > forbidden * 1.001);
        let factor = 0.5 * (forbidden + allowed_projected);
        let outside: Vec<_> = safe.iter().map(|u| u * factor).collect();
        // Both constitutive endpoint boxes still fit; an actual endpoint alone
        // violates a source box. The same refusal must apply at either end.
        for point in &response {
            let assumed = to_mandel(&total(point)).map(|e| factor * e);
            let free = to_mandel(&point.mechanical.free_strain);
            assert!(assumed.iter().all(|e| e.abs() < 0.05));
            assert!((0..6).all(|i| (assumed[i] - free[i]).abs() < 0.05));
        }
        for (before, after) in [(&safe, &outside), (&outside, &safe)] {
            assert!(op.thermoelastic_exchange(
                &laws, &temperatures, before, after, TOTAL_STRAIN,
                limits(op.cells(), response.len()), poll,
            ).is_err());
        }
    }
}

#[test]
fn native_exchange_refuses_other_sources_exhausted_work_and_cancellation_without_mutation() {
    let law = law(false, false);
    let op = cartesian(&law, ElasticityFormulation3::Displacement);
    let laws = vec![law; op.cells()];
    let temperatures: Vec<_> = (0..op.cells()).map(|i| 303.0 + 0.3 * i as f64).collect();
    let (u, _) = fields(op.nodes());
    let protected = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let original = op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll).unwrap();
    let controls = limits(op.cells(), original.len());
    let mut calls = 0;
    op.thermoelastic_exchange(
        &laws, &temperatures, &u, &u, TOTAL_STRAIN, controls,
        || { calls += 1; poll() },
    ).unwrap();
    let mut wrong = laws.clone();
    wrong[op.cells() / 2] = super::law(true, false);
    assert!(op.thermoelastic_exchange(
        &wrong, &temperatures, &u, &u, TOTAL_STRAIN, controls, poll,
    ).is_err());
    for exhausted in [
        ThermoelasticExchangeLimits3 { maximum_cells: op.cells() - 1, ..controls },
        ThermoelasticExchangeLimits3 { maximum_points: original.len() - 1, ..controls },
        ThermoelasticExchangeLimits3 { maximum_workspace_bytes: controls.maximum_workspace_bytes - 1, ..controls },
    ] {
        assert!(op.thermoelastic_exchange(
            &laws, &temperatures, &u, &u, TOTAL_STRAIN, exhausted, poll,
        ).is_err());
    }
    for stop in [1, calls / 2, calls] {
        let mut at = 0;
        assert!(matches!(op.thermoelastic_exchange(
            &laws, &temperatures, &u, &u, TOTAL_STRAIN, controls,
            || { at += 1; if at == stop { ControlFlow::Break(()) } else { poll() } },
        ), Err(ElasticityError3::Cancelled)));
    }
    op.validate_eigenstrain(&protected).unwrap();
    assert_eq!(op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll).unwrap(), original);
}
