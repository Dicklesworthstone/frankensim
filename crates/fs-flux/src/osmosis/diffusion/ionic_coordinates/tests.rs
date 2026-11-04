//! Kinematic checks independent of constitutive fluxes or nonlinear solvers.
use super::*;

fn policy(valences: &[i32]) -> ZeroCurrentDiffusion {
    ZeroCurrentDiffusion::new(valences.to_vec(), 1e-14, 1e-12).unwrap()
}

fn assert_constraints(
    original: &[f64],
    candidate: &[f64],
    volumes: &[f64],
    caps: &[f64],
    valences: &[i32],
    components: &[usize],
) {
    let n = valences.len();
    for cell in 0..volumes.len() {
        let old_charge: f64 = (0..n)
            .map(|i| f64::from(valences[i]) * original[cell * n + i])
            .sum();
        let charge: f64 = (0..n)
            .map(|i| f64::from(valences[i]) * candidate[cell * n + i])
            .sum();
        let scale: f64 = (0..n)
            .map(|i| (f64::from(valences[i]) * original[cell * n + i]).abs())
            .sum();
        assert!((charge - old_charge).abs() <= 32.0 * f64::EPSILON * (1.0 + scale));
        for i in 0..n {
            let value = candidate[cell * n + i];
            assert!(value >= 0.0 && value.is_finite());
            assert!(value / volumes[cell] <= caps[cell * n + i]);
        }
    }
    for &id in components {
        for i in 0..n {
            let old: f64 = components
                .iter()
                .enumerate()
                .filter(|(_, label)| **label == id)
                .map(|(cell, _)| original[cell * n + i])
                .sum();
            let new: f64 = components
                .iter()
                .enumerate()
                .filter(|(_, label)| **label == id)
                .map(|(cell, _)| candidate[cell * n + i])
                .sum();
            assert!((old - new).abs() <= 32.0 * f64::EPSILON * (1.0 + old.abs()));
        }
    }
}

#[test]
fn mixed_salts_preserve_every_charge_species_total_and_original_cap_on_large_probes() {
    let volumes = [1.0, 2.0, 0.5];
    let amounts = [1.0, 2.0, 3.0, 0.4, 3.0, 1.0, 4.0, 0.8, 2.0, 4.0, 6.0, 0.6];
    let caps = [20.0; 12];
    let ionic = policy(&[1, 1, -1, 0]);
    let components = [0; 3];
    let chart = IonicCoordinates::new(
        &volumes,
        &amounts,
        &caps,
        &ionic,
        &components,
        6,
        &mut || false,
    )
    .unwrap();
    assert_eq!(chart.dimension(), 6);
    assert_eq!(chart.independent_indices(), &[4, 5, 7, 8, 9, 11]);
    let mut output = [0.0; 12];
    chart
        .decode_into(chart.initial_guess(), &mut output, &mut || false)
        .unwrap();
    assert_eq!(output.map(f64::to_bits), amounts.map(f64::to_bits));
    for x in [
        [0.1, -0.2, 0.3, -0.4, 0.5, -0.6],
        [30.0, -20.0, 17.0, -12.0, 21.0, -36.0],
        [1e250, -2e250, 3e250, -4e250, 5e250, -6e250],
    ] {
        chart.decode_into(&x, &mut output, &mut || false).unwrap();
        assert_constraints(
            &amounts,
            &output,
            &volumes,
            &caps,
            ionic.valences(),
            &components,
        );
        assert!(output.iter().all(|&value| value > 0.0));
        let once = output.map(f64::to_bits);
        chart.decode_into(&x, &mut output, &mut || false).unwrap();
        assert_eq!(once, output.map(f64::to_bits));
    }
    assert!(
        chart.heap_bytes() + core::mem::size_of::<IonicCoordinates>()
            <= IonicCoordinates::payload_bound(3, 4, 6).unwrap()
    );
}

#[test]
fn boundary_reference_preserves_mixed_valences_and_exactly_absent_spectators() {
    let volumes = [1.0, 1.0];
    let amounts = [3.0, 2.0, 0.1, 0.0, 6.0, 4.0, 0.0, 0.0];
    let caps = [20.0; 8];
    let ionic = policy(&[2, -3, 0, 0]);
    let chart = IonicCoordinates::new(&volumes, &amounts, &caps, &ionic, &[0, 0], 2, &mut || false)
        .unwrap();
    assert_eq!(chart.dimension(), 2);
    let mut output = [0.0; 8];
    for x in [[0.0; 2], [0.25, -0.5], [-1e200, 1e200]] {
        chart.decode_into(&x, &mut output, &mut || false).unwrap();
        assert_constraints(
            &amounts,
            &output,
            &volumes,
            &caps,
            ionic.valences(),
            &[0, 0],
        );
        assert_eq!(output[3], 0.0);
        assert_eq!(output[7], 0.0);
        assert!(output[2] > 0.0 && output[6] > 0.0);
    }
}

#[test]
fn disconnected_components_neutral_only_pools_and_isolated_cells_remain_independent() {
    let volumes = [1.0; 5];
    let amounts = [
        1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 4.0, 4.0, 5.0,
    ];
    let caps = [20.0; 15];
    let ionic = policy(&[1, -1, 0]);
    let components = [0, 0, 2, 2, 4];
    let chart = IonicCoordinates::new(
        &volumes,
        &amounts,
        &caps,
        &ionic,
        &components,
        3,
        &mut || false,
    )
    .unwrap();
    assert_eq!(chart.dimension(), 3);
    let mut output = [0.0; 15];
    chart
        .decode_into(&[2.0, -3.0, 0.0], &mut output, &mut || false)
        .unwrap();
    assert_constraints(
        &amounts,
        &output,
        &volumes,
        &caps,
        ionic.valences(),
        &components,
    );
    assert_eq!(&output[6..], &amounts[6..]);
    chart
        .decode_into(&[0.0, 0.0, -100.0], &mut output, &mut || false)
        .unwrap();
    assert_constraints(
        &amounts,
        &output,
        &volumes,
        &caps,
        ionic.valences(),
        &components,
    );
    assert_eq!(&output[..6], &amounts[..6]);
    assert_eq!(&output[12..], &amounts[12..]);
    for index in [6, 7, 9, 10] {
        assert_eq!(output[index], 0.0);
    }
    assert_ne!(output[8], amounts[8]);
}

#[test]
fn saturated_pools_and_zero_dimensional_components_are_not_perturbed() {
    let ionic = policy(&[1, -1, 0]);
    let amounts = [1.0; 6];
    let caps = [10.0, 10.0, 1.0, 10.0, 10.0, 1.0];
    let chart = IonicCoordinates::new(&[1.0; 2], &amounts, &caps, &ionic, &[0, 0], 1, &mut || {
        false
    })
    .unwrap();
    assert_eq!(chart.dimension(), 1);
    let mut output = [0.0; 6];
    chart
        .decode_into(&[-1e100], &mut output, &mut || false)
        .unwrap();
    assert_eq!(output[2], 1.0);
    assert_eq!(output[5], 1.0);
    assert_constraints(
        &amounts,
        &output,
        &[1.0; 2],
        &caps,
        ionic.valences(),
        &[0, 0],
    );

    let isolated =
        IonicCoordinates::new(&[1.0; 2], &amounts, &caps, &ionic, &[0, 1], 0, &mut || {
            false
        })
        .unwrap();
    assert_eq!(isolated.dimension(), 0);
    isolated
        .decode_into(&[], &mut output, &mut || false)
        .unwrap();
    assert_eq!(output, amounts);
    assert!(IonicCoordinates::payload_bound(2, 3, 0).unwrap() > 0);
}

#[test]
fn reference_retains_original_admitted_charge_residuals_without_repair() {
    let amounts = [1.0 + 1e-10, 1.0, 0.0, 2.0, 2.0 + 2e-10, 1.0];
    let ionic = ZeroCurrentDiffusion::new(vec![1, -1, 0], 1e-14, 1e-9).unwrap();
    let caps = [10.0; 6];
    let chart = IonicCoordinates::new(&[1.0; 2], &amounts, &caps, &ionic, &[0, 0], 2, &mut || {
        false
    })
    .unwrap();
    let mut output = [0.0; 6];
    for x in [[0.0; 2], [2.0, -3.0]] {
        chart.decode_into(&x, &mut output, &mut || false).unwrap();
        assert_constraints(
            &amounts,
            &output,
            &[1.0; 2],
            &caps,
            ionic.valences(),
            &[0, 0],
        );
        assert!(output[0] - output[1] > 0.0);
        assert!(output[3] - output[4] < 0.0);
    }
}

#[test]
fn original_capacity_charge_dimension_and_finite_probe_refusals_remain_explicit() {
    let ionic = policy(&[1, -1]);
    let amounts = [1.0, 1.0, 2.0, 2.0];
    let caps = [10.0; 4];
    assert!(
        IonicCoordinates::new(&[1.0; 2], &amounts, &caps, &ionic, &[0, 0], 0, &mut || {
            false
        })
        .is_err()
    );
    assert!(
        IonicCoordinates::new(&[1.0; 2], &amounts, &caps, &ionic, &[0, 2], 1, &mut || {
            false
        })
        .is_err()
    );
    assert!(
        IonicCoordinates::new(
            &[1.0; 2],
            &[1.0, 2.0, 2.0, 2.0],
            &caps,
            &ionic,
            &[0, 0],
            1,
            &mut || false
        )
        .is_err()
    );
    assert!(
        IonicCoordinates::new(
            &[1.0; 2],
            &amounts,
            &[0.5; 4],
            &ionic,
            &[0, 0],
            1,
            &mut || false
        )
        .is_err()
    );
    // The conserved finite pool bounds a harmless overflowing cap*volume.
    let huge = IonicCoordinates::new(
        &[2.0; 2],
        &amounts,
        &[f64::MAX; 4],
        &ionic,
        &[0, 0],
        1,
        &mut || false,
    )
    .unwrap();
    let mut output = [0.0; 4];
    huge.decode_into(&[1.0], &mut output, &mut || false)
        .unwrap();
    assert_constraints(
        &amounts,
        &output,
        &[2.0; 2],
        &[f64::MAX; 4],
        ionic.valences(),
        &[0, 0],
    );
    assert!(
        huge.decode_into(&[f64::NAN], &mut output, &mut || false)
            .is_err()
    );
    assert!(huge.decode_into(&[], &mut output, &mut || false).is_err());
    assert!(IonicCoordinates::payload_bound(2, 2, usize::MAX).is_err());
}

#[test]
fn cancellation_never_changes_the_immutable_chart_and_retry_is_deterministic() {
    let ionic = policy(&[1, -1, 0]);
    let amounts = [1.0, 1.0, 0.5, 2.0, 2.0, 1.0, 3.0, 3.0, 1.5];
    let caps = [10.0; 9];
    assert!(matches!(
        IonicCoordinates::new(&[1.0; 3], &amounts, &caps, &ionic, &[0; 3], 4, &mut || true),
        Err(OsmoticNetworkError::Cancelled)
    ));
    let chart = IonicCoordinates::new(&[1.0; 3], &amounts, &caps, &ionic, &[0; 3], 4, &mut || {
        false
    })
    .unwrap();
    let mut output = [0.0; 9];
    let x = [1.0, -2.0, 3.0, -4.0];
    chart.decode_into(&x, &mut output, &mut || false).unwrap();
    let expected = output.map(f64::to_bits);
    let mut polls = 0;
    assert!(matches!(
        chart.decode_into(&x, &mut output, &mut || {
            polls += 1;
            polls == 3
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    chart.decode_into(&x, &mut output, &mut || false).unwrap();
    assert_eq!(output.map(f64::to_bits), expected);
}
