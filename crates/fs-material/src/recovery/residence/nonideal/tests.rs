//! G0/G1/G3: independent physical answers, conservative transport and atomic refusal.
use super::*;
use crate::recovery::residence::StirredTankCascade;

fn source() -> Provenance {
    Provenance {
        source: "authored nonideal hydraulic known answers; no measured plant calibration".into(),
        license: "project fixture".into(),
        artifact: None,
    }
}
fn geometry(active: f64, stagnant: f64) -> TankGeometry {
    TankGeometry::new(Volume::new(active), Volume::new(stagnant)).unwrap()
}
fn state(active: f64, stagnant: f64) -> ZoneInventory {
    ZoneInventory {
        active: Mass::new(active),
        stagnant: Mass::new(stagnant),
    }
}
fn train(g: &[(f64, f64)], m: &[(f64, f64)]) -> NonidealTankCascade {
    NonidealTankCascade::new(
        "passive dissolved gold",
        g.iter().map(|&(a, s)| geometry(a, s)).collect(),
        m.iter().map(|&(a, s)| state(a, s)).collect(),
        source(),
    )
    .unwrap()
}
fn flows(q: f64, bypass: f64, exchange: &[f64]) -> ResidenceFlows {
    ResidenceFlows::new(
        VolumetricFlowRate::new(q),
        Dimensionless::new(bypass),
        exchange
            .iter()
            .copied()
            .map(VolumetricFlowRate::new)
            .collect(),
    )
    .unwrap()
}
fn advance(
    t: &mut NonidealTankCascade,
    dt: f64,
    f: &ResidenceFlows,
    feed: f64,
) -> NonidealInterval {
    t.advance(Time::new(dt), f, Density::new(feed)).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance={tolerance}"
    );
}

#[test]
fn equal_tanks_reduce_to_original_residence_owner() {
    let concentrations = [1.0, 3.0, 0.0, 0.5];
    let mut ideal = StirredTankCascade::new(
        Volume::new(8.0),
        concentrations.map(Density::new).to_vec(),
        source(),
    )
    .unwrap();
    let mut actual = train(&[(2.0, 0.0); 4], &concentrations.map(|c| (2.0 * c, 0.0)));
    for (dt, q, feed) in [(4.0, 0.5, 5.0), (19.0, 0.0, 100.0), (7.0, 1.3, 0.2)] {
        let reference = ideal
            .advance(
                Time::new(dt),
                VolumetricFlowRate::new(q),
                Density::new(feed),
            )
            .unwrap();
        let result = advance(&mut actual, dt, &flows(q, 0.0, &[0.0; 4]), feed);
        for (a, b) in actual.inventories().iter().zip(ideal.tank_concentrations()) {
            close(a.active.value() / 2.0, b.value(), 3e-14);
            assert_eq!(a.stagnant.value(), 0.0);
        }
        close(
            result.train_outlet_mass.value(),
            reference.outlet_mass.value(),
            2e-13,
        );
        close(
            result.final_inventory.active.value(),
            reference.final_inventory.value(),
            2e-13,
        );
        assert_eq!(
            result.outlet_mean.is_none(),
            reference.outlet_mean.is_none()
        );
        close(result.conservation_residual.value(), 0.0, 2e-13);
    }
    assert_eq!(actual.provenance(), &source());
    assert_eq!(actual.constituent(), "passive dissolved gold");
}

#[test]
fn bypass_uses_actual_train_flow_and_isolated_volume_keeps_its_mass() {
    let mut actual = train(&[(3.0, 2.0)], &[(6.0, 14.0)]);
    let f = flows(2.0, 0.25, &[0.0]);
    let result = advance(&mut actual, 2.0, &f, 5.0);
    // Independent single-CSTR solution: qtrain=1.5, tau=3/1.5=2 seconds.
    // The bypass has zero delay and the isolated volume has no outlet path.
    let e = (-1.0_f64).exp();
    close(
        actual.inventories()[0].active.value(),
        3.0 * (5.0 - 3.0 * e),
        8e-15,
    );
    assert_eq!(actual.inventories()[0].stagnant.value(), 14.0);
    close(
        result.train_outlet_mass.value(),
        15.0 - 9.0 * (1.0 - e),
        8e-15,
    );
    assert_eq!(result.bypass_mass.value(), 5.0);
    assert_eq!(result.total_inlet_volume.value(), 4.0);
    assert_eq!(result.train_volume.value(), 3.0);
    assert_eq!(result.bypass_volume.value(), 1.0);
    assert_eq!(result.inlet_mass.value(), 20.0);
    assert_eq!(result.flows, f);
    close(
        result.outlet_mean.unwrap().value(),
        (20.0 - 9.0 * (1.0 - e)) / 4.0,
        2e-15,
    );
    close(result.active_residual.value(), 0.0, 2e-14);
    assert_eq!(result.stagnant_residual.value(), 0.0);
}

#[test]
fn unequal_two_tank_step_matches_closed_form_and_integrated_outlet() {
    let mut actual = train(&[(1.0, 0.0), (2.0, 0.0)], &[(0.0, 0.0); 2]);
    let t = 4.0_f64;
    let result = advance(&mut actual, t, &flows(1.0, 0.0, &[0.0; 2]), 1.0);
    let c1 = 1.0 - (-t).exp();
    let c2 = 1.0 - 2.0 * (-t / 2.0).exp() + (-t).exp();
    let integrated = t - 4.0 * (1.0 - (-t / 2.0).exp()) + (1.0 - (-t).exp());
    close(actual.inventories()[0].active.value(), c1, 2e-15);
    close(actual.inventories()[1].active.value() / 2.0, c2, 2e-15);
    close(result.train_outlet_mass.value(), integrated, 4e-15);
    close(result.conservation_residual.value(), 0.0, 8e-15);
}

#[test]
fn exchange_continues_during_shutdown_and_integrates_both_directions() {
    let mut actual = train(&[(2.0, 3.0)], &[(4.0, 0.0)]);
    let t = 5.0_f64;
    let qex = 0.4;
    let k = qex * (1.0 / 2.0 + 1.0 / 3.0);
    let e = (-k * t).exp();
    let integral = (1.0 - e) / k;
    let result = advance(&mut actual, t, &flows(0.0, 0.0, &[qex]), 1e6);
    close(
        actual.inventories()[0].active.value(),
        2.0 * (0.8 + 1.2 * e),
        4e-15,
    );
    close(
        actual.inventories()[0].stagnant.value(),
        3.0 * (0.8 - 0.8 * e),
        4e-15,
    );
    close(
        result.transferred_to_stagnant.value(),
        qex * (0.8 * t + 1.2 * integral),
        8e-15,
    );
    close(
        result.transferred_from_stagnant.value(),
        qex * (0.8 * t - 0.8 * integral),
        4e-15,
    );
    assert_eq!(result.inlet_mass.value(), 0.0);
    assert_eq!(result.outlet_mass.value(), 0.0);
    assert_eq!(result.outlet_mean, None);
    close(result.conservation_residual.value(), 0.0, 8e-15);
    let held = actual.clone();
    let stopped = advance(&mut actual, 1e6, &flows(0.0, 0.0, &[0.0]), 1.0);
    assert_eq!(actual, held);
    assert_eq!(stopped.substeps, 0);
}

#[test]
fn full_bypass_preserves_train_and_equal_throughput_histories_agree() {
    let mut isolated = train(&[(2.0, 3.0)], &[(4.0, 7.0)]);
    let before = isolated.clone();
    let report = advance(&mut isolated, 4.0, &flows(2.0, 1.0, &[0.0]), 3.0);
    assert_eq!(isolated, before);
    assert_eq!(report.train_outlet_mass.value(), 0.0);
    assert_eq!(report.inlet_mass.value(), 24.0);
    assert_eq!(report.outlet_mass.value(), 24.0);
    assert_eq!(report.outlet_mean.unwrap().value(), 3.0);

    let mut varied = train(&[(2.0, 1.0), (3.0, 2.0)], &[(1.0, 0.3), (0.0, 0.1)]);
    let mut whole = varied.clone();
    let a = advance(&mut varied, 2.0, &flows(1.0, 0.2, &[0.1, 0.3]), 4.0);
    let b = advance(&mut varied, 3.0, &flows(2.0, 0.2, &[0.2, 0.6]), 4.0);
    let c = advance(&mut whole, 8.0, &flows(1.0, 0.2, &[0.1, 0.3]), 4.0);
    for (a, b) in varied.inventories().iter().zip(whole.inventories()) {
        close(a.active.value(), b.active.value(), 3e-14);
        close(a.stagnant.value(), b.stagnant.value(), 3e-14);
    }
    close(
        a.outlet_mass.value() + b.outlet_mass.value(),
        c.outlet_mass.value(),
        4e-14,
    );
    close(
        a.transferred_to_stagnant.value() + b.transferred_to_stagnant.value(),
        c.transferred_to_stagnant.value(),
        3e-14,
    );
    close(
        a.transferred_from_stagnant.value() + b.transferred_from_stagnant.value(),
        c.transferred_from_stagnant.value(),
        3e-14,
    );
}

#[test]
fn stagnant_storage_produces_a_long_tail_without_losing_the_pulse() {
    let mut actual = train(&[(1.0, 1.0)], &[(0.0, 0.0)]);
    let f = flows(1.0, 0.15, &[0.1]);
    let pulse = advance(&mut actual, 1.0, &f, 1.0);
    assert_eq!(pulse.bypass_mass.value(), 0.15);
    let mut out = pulse.outlet_mass.value();
    let mut late = 0.0;
    for _ in 0..20 {
        let r = advance(&mut actual, 2.0, &f, 0.0);
        out += r.outlet_mass.value();
        late = r.outlet_mass.value();
    }
    let m = actual.total_inventory().unwrap();
    close(m.active.value() + m.stagnant.value() + out, 1.0, 3e-14);
    assert!(m.stagnant.value() > m.active.value());
    assert!(late > 1e-5, "exchange retains a measurable late tail");
    assert!(m.stagnant.value() > 0.0);
}

#[test]
fn source_scaled_recurrence_preserves_small_representable_outlet() {
    let mut actual = train(&[(1e200, 0.0); 2], &[(0.0, 0.0); 2]);
    let r = advance(&mut actual, 1.0, &flows(1.0, 0.0, &[0.0; 2]), 1e200);
    close(actual.inventories()[0].active.value(), 1e200, 2e184);
    close(actual.inventories()[1].active.value(), 0.5, 2e-16);
    close(r.outlet_mass.value(), 1e-200 / 6.0, 2e-216);
    assert!(r.outlet_mass.value() > 0.0);
    let mut tiny = train(&[(1.0, 0.0); 2], &[(0.0, 0.0); 2]);
    let before = tiny.clone();
    assert!(matches!(
        tiny.advance(
            Time::new(1e-200),
            &flows(1.0, 0.0, &[0.0; 2]),
            Density::new(1.0)
        ),
        Err(NonidealError::Arithmetic(_))
    ));
    assert_eq!(tiny, before);
}

#[test]
fn invalid_sources_geometry_and_flows_refuse() {
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(TankGeometry::new(Volume::new(bad), Volume::new(0.0)).is_err());
        assert!(TankGeometry::new(Volume::new(1.0), Volume::new(bad)).is_err());
        assert!(
            ResidenceFlows::new(
                VolumetricFlowRate::new(bad),
                Dimensionless::new(0.0),
                vec![VolumetricFlowRate::new(0.0)]
            )
            .is_err()
        );
        assert!(
            ResidenceFlows::new(
                VolumetricFlowRate::new(1.0),
                Dimensionless::new(bad),
                vec![VolumetricFlowRate::new(0.0)]
            )
            .is_err()
        );
    }
    assert!(TankGeometry::new(Volume::new(0.0), Volume::new(1.0)).is_err());
    assert!(
        ResidenceFlows::new(
            VolumetricFlowRate::new(1.0),
            Dimensionless::new(1.1),
            vec![VolumetricFlowRate::new(0.0)]
        )
        .is_err()
    );
    assert!(
        NonidealTankCascade::new(
            "Au",
            vec![geometry(1.0, 0.0)],
            vec![state(0.0, 1.0)],
            source()
        )
        .is_err()
    );
    let mut missing = source();
    missing.source.clear();
    assert!(
        NonidealTankCascade::new(
            "Au",
            vec![geometry(1.0, 1.0)],
            vec![state(0.0, 1.0)],
            missing
        )
        .is_err()
    );
    assert!(
        NonidealTankCascade::new(
            "",
            vec![geometry(1.0, 0.0)],
            vec![state(0.0, 0.0)],
            source()
        )
        .is_err()
    );
}

#[test]
fn refusal_and_mid_step_cancellation_leave_the_original_inventory_unchanged() {
    let mut actual = train(&[(1.0, 1.0), (2.0, 0.0)], &[(1.0, 2.0), (3.0, 0.0)]);
    let before = actual.clone();
    for (dt, f, feed) in [
        (-1.0, flows(1.0, 0.0, &[0.1, 0.0]), 1.0),
        (1.0, flows(1.0, 0.0, &[0.1]), 1.0),
        (1.0, flows(1.0, 0.0, &[0.1, 0.2]), 1.0),
        (1.0, flows(1.0, 0.0, &[0.1, 0.0]), -1.0),
        (1e6, flows(1.0, 0.0, &[0.1, 0.0]), 1.0),
    ] {
        assert!(
            actual
                .advance(Time::new(dt), &f, Density::new(feed))
                .is_err()
        );
        assert_eq!(actual, before);
    }
    let mut calls = 0;
    let r = actual.advance_with_cancel(
        Time::new(20.0),
        &flows(1.0, 0.2, &[0.1, 0.0]),
        Density::new(4.0),
        || {
            calls += 1;
            calls == 20
        },
    );
    assert_eq!(r, Err(NonidealError::Cancelled));
    assert_eq!(actual, before);
    assert_eq!(calls, 20);
}

#[test]
fn long_interval_and_partitioned_interval_preserve_all_fluxes() {
    let mut full = train(
        &[(1.0, 0.5), (2.0, 1.0), (0.5, 0.25)],
        &[(2.0, 0.0), (0.1, 0.4), (0.0, 0.0)],
    );
    let mut parts = full.clone();
    let f = flows(0.7, 0.1, &[0.15, 0.1, 0.2]);
    let report = advance(&mut full, 100.0, &f, 0.3);
    let mut sums = [0.0; 3];
    for _ in 0..25 {
        let p = advance(&mut parts, 4.0, &f, 0.3);
        for (sum, x) in sums.iter_mut().zip([
            p.outlet_mass.value(),
            p.transferred_to_stagnant.value(),
            p.transferred_from_stagnant.value(),
        ]) {
            *sum += x;
        }
    }
    for (a, b) in full.inventories().iter().zip(parts.inventories()) {
        close(a.active.value(), b.active.value(), 1e-14);
        close(a.stagnant.value(), b.stagnant.value(), 1e-14);
    }
    for (a, b) in sums.into_iter().zip([
        report.outlet_mass.value(),
        report.transferred_to_stagnant.value(),
        report.transferred_from_stagnant.value(),
    ]) {
        close(a, b, 5e-13);
    }
    assert!(report.substeps > 1);
}
