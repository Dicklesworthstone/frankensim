//! Independent continuous-time limits around the original native endpoints.
//! All coefficients are authored verification data, not measured materials.
use super::*;
use fs_blake3::ContentHash;
use fs_couple::moisture::{adaptive::*, cartesian::*};

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}
fn time_controls(step: f64, tolerance: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: step,
        minimum_step_s: 1e-6,
        maximum_step_s: step,
        maximum_evaluations: 5000,
        moisture_absolute_tolerance: tolerance,
        activity_absolute_tolerance: tolerance,
        reservoir_water_absolute_tolerance_kg: tolerance,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}
fn grid(
    lengths_m: [f64; 3],
    cells: [usize; 3],
    diffusion: f64,
    initial: Vec<f64>,
    cx: &Cx<'_>,
) -> MoistureCartesian {
    MoistureCartesian::new(
        card([diffusion; 2], 0.3),
        point(0.1),
        ["moisture_diffusivity"; 3],
        MaterialPropertySelection::SingleClaimOnly,
        MoistureCartesianGeometry {
            lengths_m,
            cells,
            frame: ContentHash([4; 32]),
            maximum_moisture: 0.3,
            maximum_cells: 256,
        },
        initial,
        cx,
    )
    .unwrap()
}

#[test]
fn whole_duration_refinement_converges_to_discrete_diffusion_semigroup_and_replays() {
    with_cx(|cx| {
        let n = 8;
        let diffusion = 0.1;
        let wave = |i: usize| (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos();
        let profile: Vec<_> = (0..n).map(|i| 0.15 + 0.05 * wave(i)).collect();
        let original = slab(card([diffusion; 2], 0.3), 1.0, profile.clone(), cx)
            .with_implicit_diffusivity(nonlinear())
            .unwrap();
        let sealed = [MoistureFace::Sealed, MoistureFace::Sealed];
        let loose_time = time_controls(1.0, 1e-3);
        let tight_time = time_controls(1.0, 1e-5);
        let (loose, loose_report) = original
            .advance_adaptive_to(1.0, &sealed, config(), loose_time, cx)
            .unwrap();
        let (fine, report) = original
            .advance_adaptive_to(1.0, &sealed, config(), tight_time, cx)
            .unwrap();
        let dx = 1.0 / n as f64;
        let decay =
            4.0 * diffusion / (dx * dx) * (std::f64::consts::PI / (2.0 * n as f64)).sin().powi(2);
        // Exact continuous-time solution of the admitted spatial discretization;
        // no backward-Euler update is reused as the temporal oracle.
        let error = |state: &MoistureSlab| {
            state
                .moisture()
                .iter()
                .enumerate()
                .map(|(i, &m)| (m - (0.15 + 0.05 * (-decay).exp() * wave(i))).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(error(&fine) < 0.25 * error(&loose));
        assert!(error(&fine) < 5e-4);
        assert!(report.evaluations > loose_report.evaluations);
        assert!(report.rejected_intervals > 0);
        assert!(report.maximum_accepted_error_ratio <= 1.0);
        assert_eq!(report.interval_s, [0.0, 1.0]);
        assert_eq!(fine.time_s(), 1.0);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert_eq!(report.boundary_water_kg, [0.0; 2]);
        near(
            report.water_change_kg,
            0.0,
            tight_time.water_balance_tolerance_kg,
        );
        near(fine.water_mass_kg(), original.water_mass_kg(), 1e-12);
        let (repeat, repeated) = original
            .advance_adaptive_to(1.0, &sealed, config(), tight_time, cx)
            .unwrap();
        assert_eq!(repeat.moisture(), fine.moisture());
        assert_eq!(repeated.evaluations, report.evaluations);
        assert_eq!(repeated.boundary_water_kg, report.boundary_water_kg);
        // A real Cartesian grid with sealed transverse faces has the same mode.
        let spatial = grid([1.0, 0.02, 0.05], [n, 1, 1], diffusion, profile.clone(), cx);
        let (spatial_next, spatial_report) = spatial
            .advance_adaptive_to(
                1.0,
                &std::array::from_fn(|_| MoistureFace::Sealed),
                config(),
                nonlinear(),
                tight_time,
                cx,
            )
            .unwrap();
        for (&a, &b) in spatial_next.moisture().iter().zip(fine.moisture()) {
            near(a, b, 2e-10);
        }
        assert_eq!(spatial_report.boundary_water_kg, [0.0; 6]);
        assert_eq!(spatial_next.time_s(), 1.0);
        assert_eq!(original.moisture(), profile);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn cartesian_finite_vapor_approaches_two_store_decay_and_refuses_late_partial_work() {
    with_cx(|cx| {
        let initial_rh = 0.8;
        let slope = 0.25;
        let diffusion = 0.01;
        let target = 0.02;
        let original = grid([0.1; 3], [1; 3], diffusion, vec![0.1], cx);
        let air = atmosphere::atmosphere(initial_rh, 1.0);
        let original_gas = air.gas().clone();
        let carrier = atmosphere::carrier_mass(&air);
        let vapor = air.vapor_mass_kg().unwrap();
        let capacity = vapor / initial_rh;
        let dry = original.dry_mass_kg();
        let faces = std::array::from_fn(|_| Some(atmosphere::sorption(slope)));
        let controls = time_controls(target, 1e-5);
        let (next_air, next, report) = air
            .advance_cartesian_adaptive_to(
                &original,
                target,
                &faces,
                config(),
                nonlinear(),
                atmosphere::controls(),
                controls,
                cx,
            )
            .unwrap();
        // The exact two-inventory ODE is dM/dt = k*(s*RH-M),
        // C*dRH/dt = -m_dry*dM/dt, including all six half-cell faces.
        let k = 12.0 * diffusion / 0.1_f64.powi(2);
        let b = 1.0 + slope * dry / capacity;
        let equilibrium = (slope * initial_rh + slope * dry * 0.1 / capacity) / b;
        let expected = equilibrium + (0.1 - equilibrium) * (-k * b * target).exp();
        near(next.moisture()[0], expected, 5e-5);
        near(
            next_air.relative_humidity().unwrap(),
            initial_rh - dry * (expected - 0.1) / capacity,
            2e-3,
        );
        assert!(report.rejected_intervals > 0);
        assert_eq!(next.time_s(), target);
        assert!(next_air.relative_humidity().unwrap() < initial_rh);
        near(
            report.water_change_kg,
            dry * (next.moisture()[0] - 0.1),
            1e-14,
        );
        near(
            report.boundary_water_kg.iter().sum(),
            report.water_change_kg,
            controls.water_balance_tolerance_kg,
        );
        near(
            next_air.vapor_mass_kg().unwrap() + next.water_mass_kg(),
            vapor + original.water_mass_kg(),
            controls.water_balance_tolerance_kg,
        );
        near(atmosphere::carrier_mass(&next_air), carrier, 1e-12);
        assert!(
            report.maximum_carrier_mass_residual_kg.unwrap() <= controls.water_balance_tolerance_kg
        );
        assert!(
            report.total_mass_residual_kg.unwrap().abs() <= controls.water_balance_tolerance_kg
        );
        near(
            report.water_balance_residual_kg,
            0.0,
            controls.water_balance_tolerance_kg,
        );
        for source in report.last_substep.sorption.iter().flatten() {
            near(
                source.query().axes()["relative_humidity"],
                next_air.relative_humidity().unwrap(),
                1e-11,
            );
        }
        // Three endpoint calls admit the first private interval, but cannot
        // complete the second. No usable partial field or gas is returned.
        let mut short_work = time_controls(0.01, 1.0);
        short_work.maximum_evaluations = 3;
        assert!(
            matches!(air.advance_cartesian_adaptive_to(&original, target,
            &faces, config(), nonlinear(), atmosphere::controls(), short_work, cx),
            Err(MoistureAdvanceError::WorkLimit { time_s, accepted_intervals: 1,
                evaluations: 3 }) if time_s == 0.01)
        );
        assert_eq!(air.gas(), &original_gas);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn pure_solution_adaptive_accuracy_resolves_water_when_activity_is_constant() {
    with_cx(|cx| {
        let volume = 1e-3;
        let target = 10.0;
        let original = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let spatial = grid([1.0, 0.02, 0.05], [1; 3], 0.01, vec![0.1], cx);
        let bath = solution::reservoir(volume, vec![0.0]);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let controls = time_controls(target, 1e-5);
        let (next_bath, next, report) = bath
            .advance_adaptive_to(
                &original,
                target,
                &faces,
                config(),
                solution::controls(volume),
                controls,
                cx,
            )
            .unwrap();
        let expected = 0.2 - 0.1 * (-0.04 * target).exp();
        near(next.moisture()[0], expected, 5e-4);
        assert_eq!(report.activity, Some([1.0; 2]));
        assert!(report.rejected_intervals > 0);
        assert!(report.accepted_intervals > 1);
        near(
            bath.solvent_mass_kg() - next_bath.solvent_mass_kg(),
            original.dry_mass_kg() * (next.moisture()[0] - 0.1),
            controls.water_balance_tolerance_kg,
        );
        near(
            report.water_balance_residual_kg,
            0.0,
            controls.water_balance_tolerance_kg,
        );
        assert_eq!(next_bath.solution().solute_moles(), &[0.0]);
        // The same controller and original finite bath apply to genuine 3D
        // cells; the four transverse faces are sealed, without fake slab storage.
        let faces3 = [faces[0].clone(), faces[1].clone(), None, None, None, None];
        let (spatial_bath, spatial_next, spatial_report) = bath
            .advance_cartesian_adaptive_to(
                &spatial,
                target,
                &faces3,
                config(),
                nonlinear(),
                solution::controls(volume),
                controls,
                cx,
            )
            .unwrap();
        near(spatial_next.moisture()[0], next.moisture()[0], 1e-10);
        near(
            spatial_bath.solvent_mass_kg(),
            next_bath.solvent_mass_kg(),
            1e-11,
        );
        assert_eq!(&spatial_report.boundary_water_kg[2..], &[0.0; 4]);
        assert_eq!(bath.solution().solvent_volume_m3(), volume);
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn adaptive_retries_typed_increment_but_domain_and_cancellation_are_terminal() {
    with_cx(|cx| {
        let original = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let mut step = config();
        step.maximum_content_increment = 0.002;
        assert!(matches!(
            original.advance(1.0, &faces(0.2), step, cx),
            Err(MoistureError::ContentIncrement { cell: 0, .. })
        ));
        let (next, report) = original
            .advance_adaptive_to(2.0, &faces(0.2), step, time_controls(1.0, 1e-3), cx)
            .unwrap();
        assert!(report.rejected_intervals > 0);
        assert_eq!(next.time_s(), 2.0);
        assert!(next.moisture()[0] > 0.1 && next.moisture()[0] < 0.2);
        assert!(matches!(
            original.advance_adaptive_to(2.0, &faces(0.31), step, time_controls(1.0, 1e-3), cx),
            Err(MoistureAdvanceError::Step {
                evaluations: 1,
                error: MoistureError::Invalid(_),
                ..
            })
        ));
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let original = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], &cx);
        gate.request();
        assert!(matches!(
            original.advance_adaptive_to(1.0, &faces(0.2), config(), time_controls(1.0, 1e-3), &cx),
            Err(MoistureAdvanceError::Cancelled)
        ));
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}
