//! Pair events share the original finite inventory, photons and native heat.

use super::beams::{accounts, bind, inventory, limits, ray, sections};
use super::*;
use fs_render::volumes::heating::luminescence::saturation::annihilation::{
    AnnihilationLimits, PhosphorescenceAnnihilation, TripletAnnihilation,
};
use fs_render::volumes::heating::luminescence::saturation::beam::DepletingPumpBeam;

fn pair_limits() -> AnnihilationLimits {
    AnnihilationLimits {
        absolute_energy_tolerance_j: 1e-11,
        relative_tolerance: 2e-6,
        maximum_substeps: 4096,
    }
}

fn pair(grid_cell: usize, scaled_rate: f64, yield_: f64) -> PhosphorescenceAnnihilation {
    PhosphorescenceAnnihilation {
        grid_cell,
        law: TripletAnnihilation::try_new(scaled_rate / COUNT, 1.0, yield_).unwrap(),
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-4 * actual.abs().max(expected.abs()).max(1e-30),
        "{actual:.17e} != {expected:.17e}",
    );
}

#[test]
fn g1_finite_pairs_match_independent_rational_decay_and_native_heat() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, transparent_host());
        let (thermal, heating) = inventory(
            &medium,
            &[[0.1 * COUNT, 0.1 * COUNT, 0.8 * COUNT]],
            [1e12, 4.0],
            0.0,
        );
        let law = pair(0, 1.0, 0.4);
        let heating = heating
            .with_triplet_annihilation(&[law], pair_limits())
            .unwrap();
        assert_eq!(heating.triplet_annihilation_law(0), Some(law.law));
        let (warm, accepted, audit) = thermal
            .advance_observed(&heating, 2.0, config(2.0), || Ok(()))
            .unwrap();
        // Solve T'=-gamma*T-2*k*T²/V independently. The singlet lifetime
        // changes its initial/store count by less than 2e-12 during this test.
        let gamma = 0.25_f64;
        let duration = 2.0_f64;
        let initial_t = 0.8 * COUNT;
        let twice_k = 2.0 / COUNT;
        let denominator = 1.0 + twice_k * initial_t / gamma * -(-gamma * duration).exp_m1();
        let triplets = initial_t * (-gamma * duration).exp() / denominator;
        let residence_t = denominator.ln() / twice_k;
        let pairs = (initial_t - triplets - gamma * residence_t) / 2.0;
        let singlets = 0.1 * COUNT + 0.4 * pairs;
        let populations = accepted.populations(0).unwrap();
        near(populations[2], triplets);
        near(populations[1], singlets);
        near(populations[0], COUNT - triplets - singlets);
        near(
            accepted.energy().annihilation_heat_j,
            pairs * (2.0 * HC / 600.0 - 0.4 * HC / 500.0),
        );
        near(
            accepted.energy().escaped_j[1],
            gamma * residence_t * HC / 800.0,
        );
        assert_eq!(accepted.absorbed_pump_energy_j(), [0.0; 2]);
        assert_eq!(accepted.energy().reverse_crossing_heat_absorbed_j, 0.0);
        relative(warm.energy_j(), accepted.thermal_input_j());
        relative(audit.external_heat_j, accepted.thermal_input_j());
        balance(&accepted);
    });
}

#[test]
fn g1_finite_pairs_change_the_prescribed_pump_saturation_equilibrium() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, [0.0, 5.0], [COUNT, 0.0, 0.0], 0.0);
        let heating = heating
            .with_triplet_annihilation(&[pair(0, 2.0, 0.0)], pair_limits())
            .unwrap();
        let (warm, accepted, audit) = thermal
            .advance_observed(&heating, 4.0, config(0.25), || Ok(()))
            .unwrap();
        // No singlet channel: W*(N-T)=gamma*T+2*k*T²/V.
        // The transient decays by more than exp(-40) over these four seconds.
        let total_linear = 5.0_f64 + 0.2;
        let fraction = 10.0 / (total_linear + (total_linear * total_linear + 80.0).sqrt());
        near(accepted.populations(0).unwrap()[2], fraction * COUNT);
        near(
            accepted.populations(0).unwrap()[0],
            (1.0 - fraction) * COUNT,
        );
        assert_eq!(accepted.populations(0).unwrap()[1], 0.0);
        assert_eq!(accepted.energy().escaped_j[0], 0.0);
        assert!(accepted.energy().annihilation_heat_j > 0.0);
        relative(warm.energy_j(), accepted.thermal_input_j());
        relative(audit.external_heat_j, accepted.thermal_input_j());
        balance(&accepted);
    });
}

#[test]
fn g1_finite_pairs_recover_ground_for_shared_beams_and_emit_delayed_singlets() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, transparent_host());
        let (thermal, base) = inventory(&medium, &[[COUNT, 0.0, 0.0]], [0.05, 5.0], 0.0);
        let cross_sections = sections(1, [[0.0; 2], [0.0, 1.0 / COUNT]]);
        let watt = COUNT * HC / 450.0;
        let beams = [3.0, 2.0].map(|weight| DepletingPumpBeam {
            source_identity: ContentHash([186; 32]),
            beam: ray(1, weight == 2.0, [400.0, 450.0], [0.0, weight * watt]),
            cross_sections: &cross_sections,
            pump_window_s: None,
        });
        // Binding before or after the beam owner uses the same physical law.
        let combined = base
            .clone()
            .with_triplet_annihilation(&[pair(0, 1.0, 0.5)], pair_limits())
            .unwrap()
            .with_depleting_beam(
                cx,
                &medium,
                ContentHash([187; 32]),
                ray(1, false, [400.0, 450.0], [0.0, 5.0 * watt]),
                &cross_sections,
                &thermal,
                1,
                limits().column,
                || Ok(()),
            )
            .unwrap()
            .with_pump_window(0.0, 0.4)
            .unwrap();
        let split = bind(cx, &medium, &thermal, base.clone(), &beams)
            .with_triplet_annihilation(&[pair(0, 1.0, 0.5)], pair_limits())
            .unwrap()
            .with_pump_window(0.0, 0.4)
            .unwrap();
        let control = bind(cx, &medium, &thermal, base, &beams)
            .with_pump_window(0.0, 0.4)
            .unwrap();
        let (_, one, _) = thermal
            .advance_observed(&combined, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let (charged_thermal, charged, audit) = thermal
            .advance_observed(&split, 0.4, config(0.4), || Ok(()))
            .unwrap();
        let (linear_thermal, linear, _) = thermal
            .advance_observed(&control, 0.4, config(0.4), || Ok(()))
            .unwrap();
        for (x, y) in one
            .populations(0)
            .unwrap()
            .into_iter()
            .zip(charged.populations(0).unwrap())
        {
            near(x, y);
        }
        near(
            one.energy().annihilation_heat_j,
            charged.energy().annihilation_heat_j,
        );
        near(
            one.absorbed_pump_energy_j()[1],
            charged.absorbed_pump_energy_j()[1],
        );
        assert!(charged.populations(0).unwrap()[0] > linear.populations(0).unwrap()[0]);
        assert!(charged.absorbed_pump_energy_j()[1] > linear.absorbed_pump_energy_j()[1]);
        assert_eq!(charged.absorbed_pump_energy_j()[0], 0.0);
        assert_eq!(linear.populations(0).unwrap()[1], 0.0);
        assert!(charged.energy().annihilation_heat_j > 0.0);
        relative(audit.external_heat_j, charged.thermal_input_j());
        accounts(&charged);
        accounts(&one);

        let (_, dark, _) = charged_thermal
            .advance_observed(&charged, 1.0, config(0.6), || Ok(()))
            .unwrap();
        let (_, without_pairs, _) = linear_thermal
            .advance_observed(&linear, 1.0, config(0.6), || Ok(()))
            .unwrap();
        assert_eq!(
            dark.absorbed_pump_energy_j(),
            charged.absorbed_pump_energy_j()
        );
        assert_eq!(without_pairs.cell_emission_w()[0].1[0], 0.0);
        assert!(dark.cell_emission_w()[0].1[0] > 0.0);
        assert!(dark.energy().escaped_j[0] > charged.energy().escaped_j[0]);
        assert!(dark.energy().annihilation_heat_j > charged.energy().annihilation_heat_j);
        let view = LuminescenceGrid::try_new(
            cx,
            &medium,
            dark.emission_source_identity(),
            [650.0, 800.0],
            dark.cell_emission_w(),
            1,
        )
        .unwrap();
        let observed = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
            )
            .unwrap();
        relative(
            observed.emitted_radiance_w_m2_sr[0],
            dark.cell_emission_w()[0].1[0] / (4.0 * core::f64::consts::PI),
        );
        accounts(&dark);
    });
}

#[test]
fn g4_finite_pairs_include_unilluminated_cells_and_preserve_retry_and_zero_law_behavior() {
    with_cx(|cx| {
        let density = [1.0; 2];
        let medium = grid(cx, &density, transparent_host());
        let (thermal, base) = inventory(
            &medium,
            &[[COUNT, 0.0, 0.0], [0.2 * COUNT, 0.0, 0.8 * COUNT]],
            [1e12, 4.0],
            0.0,
        );
        let cross_sections = sections(1, [[1.0 / COUNT, 0.0], [0.0; 2]]);
        let mut beam = ray(2, false, [400.0, 450.0], [COUNT * HC / 400.0, 0.0]);
        beam.max_distance_m = 1.5; // Stops inside cell0; cell1 receives no pump photons.
        let beams = [DepletingPumpBeam {
            source_identity: ContentHash([188; 32]),
            beam,
            cross_sections: &cross_sections,
            pump_window_s: None,
        }];
        let base = bind(cx, &medium, &thermal, base, &beams);
        let zero = base
            .clone()
            .with_triplet_annihilation(&[pair(1, 0.0, 0.5)], pair_limits())
            .unwrap();
        let (_, no_pair, _) = thermal
            .advance_observed(&base, 0.2, config(0.2), || Ok(()))
            .unwrap();
        let (_, zero_pair, _) = thermal
            .advance_observed(&zero, 0.2, config(0.2), || Ok(()))
            .unwrap();
        assert_eq!(no_pair.energy(), zero_pair.energy());
        assert_eq!(
            no_pair.depleting_pump_report(),
            zero_pair.depleting_pump_report()
        );
        assert_eq!(no_pair.populations(1), zero_pair.populations(1));
        let owner = base
            .clone()
            .with_triplet_annihilation(&[pair(1, 1.0, 0.5)], pair_limits())
            .unwrap();
        let original_report = owner.depleting_pump_report().unwrap().clone();
        let mut candidate = owner.clone();
        candidate
            .step_sources(&thermal, 0.0, 0.2, 0.2, &mut || Ok(()))
            .unwrap();
        let retried = candidate
            .step_sources(&thermal, 0.0, 0.1, 0.1, &mut || Ok(()))
            .unwrap();
        let fresh = owner
            .clone()
            .step_sources(&thermal, 0.0, 0.1, 0.1, &mut || Ok(()))
            .unwrap();
        assert_eq!(retried, fresh);
        assert_eq!(candidate.energy(), owner.energy());
        assert_eq!(candidate.populations(1), owner.populations(1));
        assert_eq!(candidate.depleting_pump_report().unwrap(), &original_report);
        let mut polls = 0;
        assert!(matches!(
            thermal.advance_observed(&owner, 0.2, config(0.2), || {
                polls += 1;
                if polls > 30 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(OpticalHeatingError::Thermal(DepositionError::Cancelled))
        ));
        assert!(polls > 30);
        assert_eq!(owner.depleting_pump_report().unwrap(), &original_report);
        let (_, accepted, _) = thermal
            .advance_observed(&owner, 0.2, config(0.2), || Ok(()))
            .unwrap();
        assert_eq!(accepted.populations(0), no_pair.populations(0));
        assert_eq!(
            accepted.absorbed_pump_energy_j(),
            no_pair.absorbed_pump_energy_j()
        );
        assert!(accepted.populations(1).unwrap()[1] > 0.0);
        assert!(accepted.populations(1).unwrap()[2] < no_pair.populations(1).unwrap()[2]);
        relative(accepted.populations(1).unwrap().iter().sum(), COUNT);
        accounts(&accepted);

        let limited = base
            .clone()
            .with_triplet_annihilation(
                &[pair(1, 20.0, 0.5)],
                AnnihilationLimits {
                    absolute_energy_tolerance_j: 0.0,
                    relative_tolerance: 1e-12,
                    maximum_substeps: 2,
                },
            )
            .unwrap();
        assert!(
            thermal
                .advance_observed(&limited, 1.0, config(1.0), || Ok(()))
                .is_err()
        );
        assert_eq!(limited.time_s(), 0.0);
        assert_eq!(limited.energy(), base.energy());
        assert!(
            base.clone()
                .with_triplet_annihilation(&[pair(9, 1.0, 0.5)], pair_limits())
                .is_err()
        );
        assert!(
            base.with_triplet_annihilation(&[pair(1, 1.0, 0.5), pair(1, 1.0, 0.5)], pair_limits())
                .is_err()
        );
    });
}
