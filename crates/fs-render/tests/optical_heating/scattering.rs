use super::*;
use fs_render::volumes::heating::scattering::{ScatteringHeating, ScatteringHeatingTolerance};
use fs_render::volumes::transport::scattering::deposition::{
    ScatteringDepositionConfig, ScatteringPowerDeposition,
};

fn sampled(cx: &Cx<'_>, a: f64, s: f64, depth: u32) -> ScatteringPowerDeposition<1> {
    let medium = ScatteringGrid::try_new(
        cx,
        grid(
            cx,
            &[1.0, 2.0],
            MassAbsorptionSpectrum::constant(a).unwrap(),
        ),
        MassScatteringSpectrum::constant(s).unwrap(),
        PhaseFunction::Isotropic,
    )
    .unwrap();
    medium
        .deposit_multiple_scattering_power(
            cx,
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            &[550.0],
            &[100.0],
            ScatteringDepositionConfig {
                histories_per_lane: 2_000,
                max_collisions_per_history: depth,
                max_total_segments: 200_000,
                max_deposition_cells: 2,
                seed: 125,
                source_id: 7,
            },
        )
        .unwrap()
}
fn allowance() -> ScatteringHeatingTolerance {
    ScatteringHeatingTolerance {
        power_roundoff_w: 1e-10,
        max_unresolved_w: 100.0,
        max_total_standard_error_w: 10.0,
        max_cell_standard_error_w: 10.0,
    }
}

#[test]
fn multiple_scattering_heats_native_enthalpy_and_keeps_fixed_field_error_correlated_in_time() {
    with_cx(|cx| {
        let report = sampled(cx, 0.7, 2.0, 64);
        let absorbed = report.absorbed.mean_w[0];
        let error = report.absorbed.standard_error_w[0];
        let original = thermal(vec![cell(11, 2.0, 2.0, 2.0)], &SENSIBLE, vec![]);
        let heating = ScatteringHeating::try_new(
            report,
            &[(0, CellId(11)), (1, CellId(11))],
            &original,
            allowance(),
            || Ok(()),
        )
        .unwrap();
        let (whole, whole_heat, audit) = original
            .advance_observed(&heating, 2.0, config(0.17), || Ok(()))
            .unwrap();
        let (middle, first_heat, _) = original
            .advance_observed(&heating, 0.75, config(0.1), || Ok(()))
            .unwrap();
        let (split, split_heat, _) = middle
            .advance_observed(&first_heat, 2.0, config(0.1), || Ok(()))
            .unwrap();
        let expected_j = 2.0 * (absorbed + 2.0);
        close(audit.external_heat_j, expected_j);
        close(whole.cells()[&CellId(11)].specific_h_j_kg, expected_j / 2.0);
        close(split.cells()[&CellId(11)].specific_h_j_kg, expected_j / 2.0);
        close(
            SENSIBLE
                .temperature_k(whole.cells()[&CellId(11)].specific_h_j_kg)
                .unwrap(),
            300.0 + expected_j / 8.0,
        );
        close(whole_heat.energy().thermal_input_j, 2.0 * absorbed);
        close(
            split_heat.energy().absorbed_standard_error_j[0],
            2.0 * error,
        );
        close(
            whole_heat.energy().absorbed_standard_error_j[0],
            2.0 * error,
        );
        let account = whole_heat.energy();
        close(
            account.incident_j[0],
            account.ballistic_escape_j[0]
                + account.diffuse_escape_j[0]
                + account.absorbed_j[0]
                + account.unresolved_j[0]
                + account.residual_j[0],
        );
        assert!(account.diffuse_escape_j[0] > 0.0);
    });
}

#[test]
fn pure_scattering_and_unresolved_packets_receive_no_thermal_credit() {
    with_cx(|cx| {
        let report = sampled(cx, 0.0, 4.0, 1);
        assert!(report.unresolved.mean_w[0] > 50.0);
        assert_eq!(report.absorbed.mean_w, [0.0]);
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let heating =
            ScatteringHeating::try_new(report, &[], &original, allowance(), || Ok(())).unwrap();
        let (after, account, audit) = original
            .advance_observed(&heating, 2.0, config(0.1), || Ok(()))
            .unwrap();
        assert_eq!(after.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(account.energy().thermal_input_j, 0.0);
        assert_eq!(audit.external_heat_j, 0.0);
        assert!(account.energy().unresolved_j[0] > 100.0);
    });
}

#[test]
fn sampling_tail_and_aggregated_receiver_error_allowances_are_enforced() {
    with_cx(|cx| {
        let report = sampled(cx, 0.7, 2.0, 1);
        let original = thermal(vec![cell(11, 1.0, 2.0, 0.0)], &SENSIBLE, vec![]);
        let bindings = [(0, CellId(11)), (1, CellId(11))];
        let options = allowance();
        for tolerance in [
            ScatteringHeatingTolerance {
                max_unresolved_w: 0.0,
                ..options
            },
            ScatteringHeatingTolerance {
                max_total_standard_error_w: 0.0,
                ..options
            },
        ] {
            assert!(
                ScatteringHeating::try_new(report.clone(), &bindings, &original, tolerance, || Ok(
                    ()
                ))
                .is_err()
            );
        }
        let largest = report
            .cells
            .iter()
            .map(|c| c.absorbed.standard_error_w[0])
            .fold(0.0, f64::max);
        let combined: f64 = report
            .cells
            .iter()
            .map(|c| c.absorbed.standard_error_w[0])
            .sum();
        assert!(combined > largest);
        let tolerance = ScatteringHeatingTolerance {
            max_cell_standard_error_w: 0.5 * (largest + combined),
            ..options
        };
        assert!(
            ScatteringHeating::try_new(report, &bindings, &original, tolerance, || Ok(())).is_err()
        );
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn malformed_reports_and_failed_native_advances_leave_heat_and_radiation_unpublished() {
    with_cx(|cx| {
        let report = sampled(cx, 0.7, 2.0, 64);
        let original = thermal(vec![cell(11, 1.0, 2.0, 0.0)], &SENSIBLE, vec![]);
        let bindings = [(0, CellId(11)), (1, CellId(11))];
        let mut broken = report.clone();
        broken.cells[0].absorbed.mean_w[0] += 1.0;
        assert!(
            ScatteringHeating::try_new(broken, &bindings, &original, allowance(), || Ok(()))
                .is_err()
        );
        assert!(matches!(
            ScatteringHeating::try_new(report.clone(), &[], &original, allowance(), || Ok(())),
            Err(OpticalHeatingError::UnmappedCell(_))
        ));
        let heating =
            ScatteringHeating::try_new(report, &bindings, &original, allowance(), || Ok(()))
                .unwrap();
        let mut checkpoints = 0;
        let cancelled = original.advance_observed(&heating, 1.0, config(0.01), || {
            checkpoints += 1;
            if checkpoints > 20 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            cancelled,
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        assert!(matches!(
            original.advance_observed(&heating, 100.0, config(0.5), || Ok(())),
            Err(OpticalHeatingError::Thermal(
                DepositionError::MaterialDomain
            ))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(heating.time_s(), 0.0);
        assert_eq!(heating.energy().thermal_input_j, 0.0);
        let (_, accepted, audit) = original
            .advance_observed(&heating, 1.0, config(0.1), || Ok(()))
            .unwrap();
        close(accepted.energy().thermal_input_j, audit.external_heat_j);
    });
}
