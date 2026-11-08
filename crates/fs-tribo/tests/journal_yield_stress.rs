use fs_matdb::Provenance;
use fs_material::slurry::rheology::{
    HerschelBulkleyCard, HerschelBulkleyParameters, ShearResponse, SlurryConditions, SlurryDomain,
    SlurryTestContext,
};
use fs_qty::{Dimensionless, Frequency, Stress, Temperature};
use fs_tribo::journal::JournalGeometry;
use fs_tribo::journal::yield_stress::{
    JournalYieldRegime, YieldStressJournal, YieldStressJournalControls, YieldStressJournalError,
};

fn conditions() -> SlurryConditions {
    SlurryConditions {
        temperature: Temperature::new(298.15),
        solids_mass_fraction: Dimensionless::new(0.4),
    }
}

fn source(yield_pa: f64, excess_pa: f64, index: f64) -> HerschelBulkleyCard {
    HerschelBulkleyCard::new(
        HerschelBulkleyParameters {
            yield_stress: Stress::new(yield_pa),
            reference_excess_stress: Stress::new(excess_pa),
            reference_shear_rate: Frequency::new(1.0),
            flow_index: Dimensionless::new(index),
        },
        SlurryDomain::new(
            [Frequency::new(0.0), Frequency::new(10_000.0)],
            [Temperature::new(293.0), Temperature::new(303.0)],
            [Dimensionless::new(0.35), Dimensionless::new(0.45)],
        )
        .unwrap(),
        SlurryTestContext {
            specimen: "authored analytic circular Couette fixture".into(),
            method: "declared HB parameters; no empirical calibration claimed".into(),
            provenance: Provenance {
                source: "project-owned analytic regression".into(),
                license: "project-owned synthetic fixture".into(),
                artifact: None,
            },
        },
    )
    .unwrap()
}

fn geometry() -> JournalGeometry {
    JournalGeometry {
        inner_radius_m: 0.01,
        outer_radius_m: 0.04,
        length_m: 0.5,
    }
}

fn factor(g: JournalGeometry) -> f64 {
    core::f64::consts::TAU * g.length_m * g.inner_radius_m.powi(2)
}

fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e} (relative tolerance {relative:e})"
    );
}

#[test]
fn zero_yield_newtonian_and_power_law_limits_are_analytic() {
    let g = geometry();
    for n in [0.5, 1.0, 1.5] {
        let card = source(0.0, 2.0, n);
        let journal = YieldStressJournal::new(&card, conditions(), g, 1000.0).unwrap();
        assert!(core::ptr::eq(journal.source(), &card));
        let inner_rate: f64 = 7.0;
        let stress = 2.0 * inner_rate.powf(n);
        let expected =
            0.5 * n * inner_rate * (1.0 - (g.inner_radius_m / g.outer_radius_m).powf(2.0 / n));
        let p = journal
            .at_outer_torque(factor(g) * stress, 0.25, Default::default())
            .unwrap();
        close(p.flow().relative_angular_velocity_rad_s, expected, 3e-14);
        assert_eq!(p.flow().quadrature_nodes, 0);
        assert_eq!(p.flow().regime, JournalYieldRegime::FullyYielded);
        let at_r = p.sample(0.02).unwrap();
        close(
            at_r.shear_rate_per_s,
            inner_rate * 0.5_f64.powf(2.0 / n),
            3e-14,
        );
        close(
            at_r.angular_velocity_rad_s,
            0.25 + 0.5 * n * inner_rate * (1.0 - 0.5_f64.powf(2.0 / n)),
            3e-14,
        );
        let inverse = journal
            .at_wall_speeds([0.0, expected], Default::default())
            .unwrap();
        close(
            inverse.flow().wall_torques_nm.unwrap()[1],
            factor(g) * stress,
            2e-10,
        );
        close(
            inverse.flow().dissipated_power_w,
            factor(g) * stress * expected,
            2e-10,
        );
    }
}

#[test]
fn bingham_partial_and_full_yield_match_circular_couette_solution() {
    let g = geometry();
    let (yield_pa, viscosity) = (10.0, 2.0);
    let card = source(yield_pa, viscosity, 1.0);
    let journal = YieldStressJournal::new(&card, conditions(), g, 1000.0).unwrap();
    for stress in [10.01_f64, 20.0, 200.0] {
        let yielded_radius = (g.inner_radius_m * (stress / yield_pa).sqrt()).min(g.outer_radius_m);
        let angle = |r: f64| {
            let r = r.min(yielded_radius);
            (stress * (1.0 - (g.inner_radius_m / r).powi(2))
                - 2.0 * yield_pa * (r / g.inner_radius_m).ln())
                / (2.0 * viscosity)
        };
        let expected = angle(g.outer_radius_m);
        let p = journal
            .at_outer_torque(factor(g) * stress, 0.0, Default::default())
            .unwrap();
        close(p.flow().relative_angular_velocity_rad_s, expected, 2e-9);
        close(p.flow().yielded_outer_radius_m, yielded_radius, 2e-14);
        assert_eq!(
            p.flow().regime,
            if stress < 160.0 {
                JournalYieldRegime::PartiallyYielded
            } else {
                JournalYieldRegime::FullyYielded
            }
        );
        for radius in [g.inner_radius_m, 0.012, 0.03, g.outer_radius_m] {
            let sample = p.sample(radius).unwrap();
            close(sample.angular_velocity_rad_s, angle(radius), 2e-9);
            let local_stress = stress * (g.inner_radius_m / radius).powi(2);
            close(sample.radial_shear_stress_pa.unwrap(), local_stress, 3e-14);
            close(
                sample.shear_rate_per_s,
                ((local_stress - yield_pa) / viscosity).max(0.0),
                3e-13,
            );
            if radius > yielded_radius {
                assert!(!sample.yielded);
                assert_eq!(sample.dissipation_w_m3, 0.0);
                assert_eq!(
                    sample.angular_velocity_rad_s,
                    p.flow().accepted_wall_angular_velocities_rad_s[1]
                );
            }
        }
        let recovered = journal
            .at_wall_speeds([0.0, expected], Default::default())
            .unwrap();
        close(
            recovered.flow().wall_torques_nm.unwrap()[1],
            factor(g) * stress,
            2e-9,
        );
    }
}

#[test]
fn general_hb_volume_heat_and_exact_tangent_match_independent_radial_integrals() {
    let g = geometry();
    for n in [0.4, 1.4] {
        let card = source(10.0, 2.0, n);
        let journal = YieldStressJournal::new(&card, conditions(), g, 1000.0).unwrap();
        let stress = 30.0;
        let torque = factor(g) * stress;
        let p = journal
            .at_outer_torque(torque, 0.0, Default::default())
            .unwrap();
        let c = g.inner_radius_m * 3.0_f64.sqrt();
        let dr = (c - g.inner_radius_m) / 65_536.0;
        let (mut speed, mut heat) = (0.0, 0.0);
        for i in 0..65_536 {
            let r = g.inner_radius_m + (f64::from(i) + 0.5) * dr;
            let tau = stress * (g.inner_radius_m / r).powi(2);
            let rate = ((tau - 10.0) / 2.0).powf(1.0 / n);
            speed += rate / r * dr;
            heat += core::f64::consts::TAU * g.length_m * r * tau * rate * dr;
        }
        close(p.flow().relative_angular_velocity_rad_s, speed, 2e-8);
        close(p.flow().dissipated_power_w, heat, 2e-8);
        let delta = torque * 1e-5;
        let lower = journal
            .at_outer_torque(torque - delta, 0.0, Default::default())
            .unwrap();
        let upper = journal
            .at_outer_torque(torque + delta, 0.0, Default::default())
            .unwrap();
        let finite_difference = 2.0 * delta
            / (upper.flow().relative_angular_velocity_rad_s
                - lower.flow().relative_angular_velocity_rad_s);
        close(
            p.flow().differential_resistance_nm_s.unwrap(),
            finite_difference,
            2e-8,
        );
        let reversed = journal
            .at_outer_torque(-torque, 8.0, Default::default())
            .unwrap();
        assert_eq!(
            reversed.flow().relative_angular_velocity_rad_s,
            -p.flow().relative_angular_velocity_rad_s
        );
        assert_eq!(
            reversed.flow().dissipated_power_w,
            p.flow().dissipated_power_w
        );
        close(
            reversed.sample(0.012).unwrap().angular_velocity_rad_s,
            8.0 - p.sample(0.012).unwrap().angular_velocity_rad_s,
            1e-12,
        );
    }
}

#[test]
fn static_torque_is_set_valued_until_load_is_specified_and_dissipates_no_heat() {
    let g = geometry();
    let card = source(10.0, 2.0, 1.0);
    let journal = YieldStressJournal::new(&card, conditions(), g, 1000.0).unwrap();
    let bounds = journal.static_outer_torque_bounds_nm();
    close(bounds[1], factor(g) * 10.0, 1e-15);
    assert_eq!(bounds[0], -bounds[1]);
    let unknown = journal
        .at_wall_speeds([3.0, 3.0], Default::default())
        .unwrap();
    assert_eq!(unknown.flow().wall_torques_nm, None);
    assert_eq!(unknown.flow().differential_resistance_nm_s, None);
    assert_eq!(unknown.sample(0.02).unwrap().radial_shear_stress_pa, None);
    for torque in [bounds[0], 0.5 * bounds[0], 0.0, 0.5 * bounds[1], bounds[1]] {
        let p = journal
            .at_outer_torque(torque, 3.0, Default::default())
            .unwrap();
        assert_eq!(p.flow().regime, JournalYieldRegime::Unyielded);
        assert_eq!(p.flow().accepted_wall_angular_velocities_rad_s, [3.0, 3.0]);
        assert_eq!(p.flow().dissipated_power_w, 0.0);
        assert_eq!(p.flow().wall_power_w, 0.0);
        let s = p.sample(0.02).unwrap();
        close(
            s.radial_shear_stress_pa.unwrap(),
            torque / factor(g) * 0.25,
            1e-14,
        );
        assert_eq!(s.shear_rate_per_s, 0.0);
        assert_eq!(s.dissipation_w_m3, 0.0);
        assert_eq!(s.angular_velocity_rad_s, 3.0);
    }
}

#[test]
fn narrow_gap_limit_and_original_source_ceiling_are_preserved() {
    let card = source(10.0, 2.0, 0.6);
    let g = JournalGeometry {
        inner_radius_m: 1.0,
        outer_radius_m: 1.0 + 1e-8,
        length_m: 0.5,
    };
    let ceiling = 47.3_f64;
    let journal = YieldStressJournal::new(&card, conditions(), g, ceiling).unwrap();
    let p = journal
        .at_outer_torque(factor(g) * 30.0, 0.0, Default::default())
        .unwrap();
    let planar =
        10.0_f64.powf(1.0 / 0.6) * (g.outer_radius_m - g.inner_radius_m) / g.inner_radius_m;
    close(p.flow().relative_angular_velocity_rad_s, planar, 1e-7);
    let ShearResponse::Yielded { stress, .. } = card
        .stress_at_shear_rate(Frequency::new(ceiling), conditions())
        .unwrap()
    else {
        panic!("positive ceiling must be yielded");
    };
    let maximum_torque = factor(g) * stress.value();
    let maximum = journal
        .at_outer_torque(maximum_torque, 0.0, Default::default())
        .unwrap();
    close(maximum.flow().wall_shear_rates_per_s[0], ceiling, 1e-14);
    let accepted = journal
        .at_wall_speeds(
            [0.0, maximum.flow().relative_angular_velocity_rad_s],
            Default::default(),
        )
        .unwrap();
    close(
        accepted.flow().wall_torques_nm.unwrap()[1],
        maximum_torque,
        1e-13,
    );
    assert!(matches!(
        journal.at_outer_torque(maximum_torque * 1.001, 0.0, Default::default()),
        Err(YieldStressJournalError::UnsupportedMotion)
    ));
    assert!(matches!(
        journal.at_wall_speeds(
            [0.0, maximum.flow().relative_angular_velocity_rad_s * 1.001],
            Default::default()
        ),
        Err(YieldStressJournalError::UnsupportedMotion)
    ));
    let unsupported_conditions = SlurryConditions {
        temperature: Temperature::new(400.0),
        ..conditions()
    };
    assert!(matches!(
        YieldStressJournal::new(&card, unsupported_conditions, g, ceiling),
        Err(YieldStressJournalError::Source(_))
    ));
    assert!(matches!(
        YieldStressJournal::new(&card, conditions(), g, 10_001.0),
        Err(YieldStressJournalError::Source(_))
    ));
}

#[test]
fn cancellation_and_work_limits_refuse_without_mutating_the_source() {
    let card = source(10.0, 2.0, 0.7);
    let journal = YieldStressJournal::new(&card, conditions(), geometry(), 1000.0).unwrap();
    let reference = journal
        .at_wall_speeds([0.0, 2.0], Default::default())
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        journal.with_cancel([0.0, 2.0], Default::default(), || {
            polls += 1;
            polls == 8
        }),
        Err(YieldStressJournalError::Cancelled)
    ));
    assert!(matches!(
        reference.sample_with_cancel(0.02, || true),
        Err(YieldStressJournalError::Cancelled)
    ));
    let one_node = YieldStressJournalControls {
        maximum_quadrature_nodes: 1,
        ..Default::default()
    };
    assert!(matches!(
        journal.at_wall_speeds([0.0, 2.0], one_node),
        Err(YieldStressJournalError::QuadratureBudget)
    ));
    let one_newton = YieldStressJournalControls {
        maximum_stress_iterations: 1,
        ..Default::default()
    };
    assert!(matches!(
        journal.at_wall_speeds([0.0, 2.0], one_newton),
        Err(YieldStressJournalError::StressSolveBudget)
    ));
    assert_eq!(
        journal
            .at_wall_speeds([0.0, 2.0], Default::default())
            .unwrap(),
        reference
    );
    assert!(reference.sample(0.05).is_err());
}
