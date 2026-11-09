//! Continuous heat, linear momentum and orbital momentum share native histories.

use super::*;

fn continuous_wrench(
    cx: &Cx<'_>,
    medium: &ScatteringGrid<'_>,
    origin: Point3,
    reference: Point3,
    options: ScatteringDepositionConfig,
) -> ScatteringWrenchDeposition<1> {
    medium
        .deposit_continuous_absorption_wrench(
            cx,
            origin,
            Vec3::new(1.0, 0.0, 0.0),
            &[550.0],
            &[LIGHT_M_S],
            reference,
            options,
        )
        .unwrap()
}

fn joint_balance(report: &ScatteringWrenchDeposition<1>) {
    balance(&ScatteringMomentumDeposition {
        power: report.power.clone(),
        momentum: report.momentum.clone(),
    });
    angular_balance(report);
    let energy = &report.power;
    let sum = energy.absorbed.mean_w[0]
        + energy.ballistic_escape.mean_w[0]
        + energy.diffuse_escape.mean_w[0]
        + energy.unresolved.mean_w[0];
    assert!((sum / energy.incident_power_w[0] - 1.0).abs() < 3e-13);
}

#[test]
fn g1_weak_continuous_absorption_retains_force_and_the_original_line_of_action() {
    with_cx(|cx, _| {
        let origin = Point3::new(-1.0, 0.2, 0.8);
        let reference = Point3::new(0.3, -0.4, 0.25);
        for scattering in [0.0, 1e-12] {
            let medium = medium(cx, &[1.0, 2.0], 1e-7, scattering, PhaseFunction::Isotropic);
            let options = config(32, 8);
            let report = continuous_wrench(cx, &medium, origin, reference, options);
            let heat = medium
                .deposit_continuous_absorption_power(
                    cx,
                    origin,
                    Vec3::new(1.0, 0.0, 0.0),
                    &[550.0],
                    &[LIGHT_M_S],
                    options,
                )
                .unwrap();
            assert_eq!(report.power, heat);
            let expected_force = -(-3e-7_f64).exp_m1();
            assert!((report.momentum.force.mean_n[0][0] / expected_force - 1.0).abs() < 2e-13);
            assert_eq!(report.momentum.force.standard_error_n, [[0.0; 3]]);
            assert_eq!(report.orbital.torque.standard_error_n_m, [[0.0; 3]]);
            for cell in &report.momentum.cells {
                let heat = report
                    .power
                    .cells
                    .iter()
                    .find(|h| h.cell_index == cell.cell_index)
                    .unwrap();
                let torque = report
                    .orbital
                    .cells
                    .iter()
                    .find(|t| t.cell_index == cell.cell_index)
                    .unwrap();
                let force = heat.absorbed.mean_w[0] / LIGHT_M_S;
                assert_eq!(cell.force.mean_n[0], [force, 0.0, 0.0]);
                for (actual, expected) in
                    torque.torque.mean_n_m[0]
                        .into_iter()
                        .zip([0.0, 0.55 * force, -0.6 * force])
                {
                    assert!((actual - expected).abs() <= 2e-13 * force);
                }
            }
            if scattering == 0.0 {
                let mut original = medium
                    .deposit_multiple_scattering_wrench(
                        cx,
                        origin,
                        Vec3::new(1.0, 0.0, 0.0),
                        &[550.0],
                        &[LIGHT_M_S],
                        reference,
                        options,
                    )
                    .unwrap();
                original.power.absorption_estimator = AbsorptionEstimator::Continuous;
                assert_eq!(report, original);
            }
            joint_balance(&report);
        }
    });
}

#[test]
fn g1_each_continuous_history_closes_energy_force_and_torque_at_the_cap() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1e6], 0.2, 0.8, PhaseFunction::Isotropic);
        let wavelengths = [550.0];
        let (absorption, scattering, extinction) = medium.coefficients(&wavelengths).unwrap();
        let prepared = UniformOptics {
            medium: &medium,
            wavelengths,
            absorption,
            scattering,
            extinction,
        };
        let origin = Point3::new(0.2, 0.3, 0.4);
        let reference = Point3::new(-0.4, 0.1, 0.2);
        for collisions in [1_u32, 2] {
            let mut tallies = MomentumTallies::<1>::with_orbital(1, reference).unwrap();
            let mut history = HistoryMomentum::new(&mut tallies, 0, 0);
            let mut segments = 0;
            let mut draws = [0.5, 0.2, 0.3, 0.4, 0.8, 0.7].into_iter();
            let result = super::super::super::continuous::trace_history(
                &prepared,
                &mut || Ok(()),
                (origin, Vec3::new(1.0, 0.0, 0.0), LIGHT_M_S),
                (collisions, u64::from(collisions)),
                &mut segments,
                &mut || draws.next().unwrap(),
                &mut |_, _| Ok(()),
                Some(&mut history),
            )
            .unwrap();
            assert_eq!(segments, u64::from(collisions));
            assert_eq!(draws.len(), 6 - 3 * collisions as usize);
            assert!(result.values[2] > 0.0 && result.values[3] > 0.0);
            assert!(((result.values[2] + result.values[3]) / LIGHT_M_S - 1.0).abs() < 2e-14);
            let orbital = history.orbital.as_ref().unwrap();
            let incident_torque = [0.0, origin.z - reference.z, -(origin.y - reference.y)];
            for axis in 0..3 {
                assert!(
                    (history.force[axis].value() + history.escape[2][axis]
                        - if axis == 0 { 1.0 } else { 0.0 })
                    .abs()
                        < 3e-14
                );
                assert!(
                    (orbital.torque[axis].value() + orbital.escape[2][axis]
                        - incident_torque[axis])
                        .abs()
                        < 3e-14
                );
            }
            if collisions == 1 {
                let absorbed = -(-0.25 * std::f64::consts::LN_2).exp_m1();
                assert!((result.values[2] / LIGHT_M_S - absorbed).abs() < 3e-14);
                let event = Point3::new(
                    origin.x + std::f64::consts::LN_2 / (1e6 * 0.8),
                    origin.y,
                    origin.z,
                );
                let force = history.force.map(PowerSum::value);
                let absorption_torque =
                    orbital::moment_about(reference, origin, [absorbed, 0.0, 0.0]).unwrap();
                let scattering_torque = orbital::moment_about(
                    reference,
                    event,
                    [force[0] - absorbed, force[1], force[2]],
                )
                .unwrap();
                for axis in 0..3 {
                    assert!(
                        (orbital.torque[axis].value()
                            - absorption_torque[axis]
                            - scattering_torque[axis])
                            .abs()
                            < 3e-14
                    );
                }
            }
            history.finish().unwrap();
        }
    });
}

#[test]
fn g3_continuous_escape_replay_and_reference_shift_keep_the_same_energy_histories() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.0, 2.0, 0.5], 0.15, 0.85, PhaseFunction::Rayleigh);
        let origin = Point3::new(-1.0, 0.2, 0.8);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let options = config(512, 3);
        let zero = Point3::new(0.0, 0.0, 0.0);
        let shift = Point3::new(0.3, -0.4, 0.2);
        let a = continuous_wrench(cx, &medium, origin, zero, options);
        let b = continuous_wrench(cx, &medium, origin, shift, options);
        assert_eq!(a, continuous_wrench(cx, &medium, origin, zero, options));
        assert_eq!(a.power, b.power);
        assert_eq!(a.momentum, b.momentum);
        assert!(a.power.diffuse_escape.mean_w[0] > 0.0);
        assert!(a.power.unresolved.mean_w[0] > 0.0);
        let heat = medium
            .deposit_continuous_absorption_power(
                cx,
                origin,
                direction,
                &[550.0],
                &[LIGHT_M_S],
                options,
            )
            .unwrap();
        let force = medium
            .deposit_continuous_absorption_momentum(
                cx,
                origin,
                direction,
                &[550.0],
                &[LIGHT_M_S],
                options,
            )
            .unwrap();
        assert_eq!(a.power, heat);
        assert_eq!(a.power, force.power);
        assert_eq!(a.momentum, force.momentum);
        for (before, after, rate) in [
            (
                a.orbital.torque.mean_n_m[0],
                b.orbital.torque.mean_n_m[0],
                a.momentum.force.mean_n[0],
            ),
            (
                a.orbital.ballistic_escape.mean_n_m[0],
                b.orbital.ballistic_escape.mean_n_m[0],
                a.momentum.ballistic_escape.mean_n[0],
            ),
            (
                a.orbital.diffuse_escape.mean_n_m[0],
                b.orbital.diffuse_escape.mean_n_m[0],
                a.momentum.diffuse_escape.mean_n[0],
            ),
            (
                a.orbital.unresolved.mean_n_m[0],
                b.orbital.unresolved.mean_n_m[0],
                a.momentum.unresolved.mean_n[0],
            ),
        ] {
            let shift_moment = orbital::moment_about(zero, shift, rate).unwrap();
            for axis in 0..3 {
                assert!((after[axis] - before[axis] + shift_moment[axis]).abs() < 3e-12);
            }
        }
        joint_balance(&a);
        joint_balance(&b);
    });
}

#[test]
fn g2_local_absorbers_and_elastic_spectral_scatterers_keep_distinct_heat_and_recoil() {
    with_cx(|cx, _| {
        let base = medium(cx, &[1.0; 2], 0.0, 0.0, PhaseFunction::Isotropic);
        let phases = [
            PhaseFunction::Isotropic,
            PhaseFunction::spectral_table(&[
                (380.0, PhaseFunction::henyey_greenstein(-0.8).unwrap()),
                (780.0, PhaseFunction::henyey_greenstein(0.8).unwrap()),
            ])
            .unwrap(),
        ];
        let field = CellScatteringGrid::try_new(
            cx,
            base.absorption_grid(),
            [380.0, 780.0],
            &[[0.2; 2], [0.0; 2]],
            &[[0.0; 2], [1e6; 2]],
            CellPhaseFunction::Indexed {
                palette: &phases,
                cell_indices: &[0, 1],
            },
            (380.0, 780.0),
            ContentHash([61; 32]),
        )
        .unwrap();
        let options = config(4_096, 1);
        let origin = Point3::new(-1.0, 0.2, 0.8);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let report = field
            .deposit_continuous_absorption_wrench(
                cx,
                origin,
                direction,
                &[LIGHT_M_S; 2],
                Point3::new(0.0, 0.0, 0.0),
                options,
            )
            .unwrap();
        let force = field
            .deposit_continuous_absorption_momentum(cx, origin, direction, &[LIGHT_M_S; 2], options)
            .unwrap();
        assert_eq!(report.power, force.power);
        assert_eq!(report.momentum, force.momentum);
        assert_eq!(report.power.cells.len(), 1);
        assert_eq!(report.power.cells[0].cell_index, 0);
        assert_eq!(report.momentum.cells.len(), 2);
        let absorbed = -(-0.2_f64).exp_m1();
        for (lane, g) in [-0.8, 0.8].into_iter().enumerate() {
            assert!((report.power.absorbed.mean_w[lane] / LIGHT_M_S - absorbed).abs() < 3e-14);
            let expected = 1.0 - (-0.2_f64).exp() * g;
            assert!(
                (report.momentum.force.mean_n[lane][0] - expected).abs()
                    <= 6.0 * report.momentum.force.standard_error_n[lane][0] + 3e-14
            );
        }
        balance(&force);
        angular_balance(&report);
    });
}

#[test]
fn g4_continuous_force_and_wrench_refusals_publish_no_partial_result() {
    with_cx(|cx, gate| {
        let medium = medium(cx, &[1.0], 0.0, 1e6, PhaseFunction::Isotropic);
        let origin = Point3::new(-1.0, 0.2, 0.8);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let options = config(2, 1);
        let identity = medium.content_identity();
        assert_eq!(
            medium.deposit_continuous_absorption_momentum(
                cx,
                origin,
                direction,
                &[550.0],
                &[1.0],
                ScatteringDepositionConfig {
                    max_deposition_cells: 0,
                    ..options
                },
            ),
            Err(VolumeTransportError::Budget)
        );
        assert_eq!(
            medium.deposit_continuous_absorption_wrench(
                cx,
                origin,
                direction,
                &[550.0],
                &[1.0],
                Point3::new(f64::NAN, 0.0, 0.0),
                options,
            ),
            Err(VolumeTransportError::Ray)
        );
        let phase = PhaseFunction::Isotropic;
        let field = CellScatteringGrid::try_new(
            cx,
            medium.absorption_grid(),
            [550.0],
            &[[0.0]],
            &[[1e6]],
            CellPhaseFunction::Uniform(&phase),
            (380.0, 780.0),
            identity,
        )
        .unwrap();
        assert_eq!(
            field.deposit_continuous_absorption_wrench_with_checkpoint(
                &mut || Err(VolumeTransportError::Cancelled),
                origin,
                direction,
                &[1.0],
                Point3::new(0.0, 0.0, 0.0),
                options,
            ),
            Err(VolumeTransportError::Cancelled)
        );
        assert_eq!(
            field.deposit_continuous_absorption_momentum_with_checkpoint(
                &mut || Err(VolumeTransportError::Cancelled),
                origin,
                direction,
                &[1.0],
                options,
            ),
            Err(VolumeTransportError::Cancelled)
        );
        assert_eq!(medium.content_identity(), identity);
        gate.request();
        assert_eq!(
            medium.deposit_continuous_absorption_wrench(
                cx,
                origin,
                direction,
                &[550.0],
                &[1.0],
                Point3::new(0.0, 0.0, 0.0),
                options,
            ),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
