use super::*;
use crate::dynamics::tests::gas;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_extensive_storage_retains_source_calorics_and_incoming_kinetic_energy() {
    for linear in [0.0, 0.4] {
        let mut reference = None;
        for formation in [0.0, -1e18, 1e18] {
            let gas = gas(linear, formation);
            let volume = 0.7;
            let tank = Plenum::new(volume, gas.conserved(1.0, 0.0, 2.0).unwrap()).unwrap();
            let donor = gas.conserved(1.0, 3.0, 4.0).unwrap();
            let mass = 0.2;
            let energy = mass * (donor.energy + 4.0);
            let candidate = tank
                .exchange_candidate(&gas, mass, energy, || false)
                .unwrap();
            let p = gas.primitive(candidate.state(), || false).unwrap();
            let total_mass = volume + mass;
            let total_energy = volume * tank.state().energy + energy;
            let sensible =
                total_energy / total_mass + 2.5 * 0.001 + 0.5 * linear * 0.001_f64.powi(2);
            let temperature = if linear == 0.0 {
                sensible / 2.5
            } else {
                2.0 * sensible / (2.5 + (6.25 + 2.0 * linear * sensible).sqrt())
            };
            close(candidate.volume_m3(), volume, 0.0);
            close(candidate.state().density * volume, total_mass, 2e-15);
            close(candidate.state().energy * volume, total_energy, 2e-15);
            close(p.temperature, temperature, 2e-10);
            close(p.pressure, total_mass * temperature / volume, 2e-10);
            assert_eq!(candidate.state().momentum, 0.0);
            assert_eq!(
                tank.exchange_candidate(&gas, 0.0, 0.0, || false).unwrap(),
                tank
            );
            if let Some(previous) = reference {
                assert_eq!(candidate.state(), previous);
            } else {
                reference = Some(candidate.state());
            }
            // The frozen inventory is a genuine Copy stage value, with no heap.
            let stages = [tank; 2];
            assert_eq!(stages[0], stages[1]);
        }
    }
}

fn blowdown(dt: f64, normal: f64) -> (f64, f64, f64) {
    let gas = gas(0.0, 0.0);
    let volume = 0.5;
    let mut tank = Plenum::new(volume, gas.conserved(1.0, 0.0, 2.0).unwrap()).unwrap();
    let ambient = gas.normal_conserved(0.01, 0.0, 0.0, 0.0001).unwrap();
    let area = 0.1;
    let mut delivered_mass = 0.0;
    let mut delivered_energy = 0.0;
    for _ in 0..(0.5 / dt).round() as usize {
        let a = gas
            .reservoir_normal_flux(ambient, tank.state(), normal, || false)
            .unwrap();
        let dm_a = normal * dt * area * a.mass;
        let de_a = normal * dt * area * a.energy;
        assert!(dm_a < 0.0);
        let stage = tank.exchange_candidate(&gas, dm_a, de_a, || false).unwrap();
        let b = gas
            .reservoir_normal_flux(ambient, stage.state(), normal, || false)
            .unwrap();
        let dm = 0.5 * (dm_a + normal * dt * area * b.mass);
        let de = 0.5 * (de_a + normal * dt * area * b.energy);
        tank = tank.exchange_candidate(&gas, dm, de, || false).unwrap();
        delivered_mass -= dm;
        delivered_energy -= de;
    }
    let p = gas.primitive(tank.state(), || false).unwrap();
    close(
        volume * tank.state().density + delivered_mass,
        volume,
        2e-14,
    );
    let initial = gas.conserved(1.0, 0.0, 2.0).unwrap();
    close(
        volume * tank.state().energy + delivered_energy,
        volume * initial.energy,
        2e-14,
    );
    let expected = 2.0 * tank.state().density.powf(1.4);
    (p.pressure, p.temperature, (p.pressure - expected).abs())
}

#[test]
fn g1_source_owned_blowdown_cools_with_pressure_feedback_and_rk2_refinement() {
    let coarse = blowdown(0.01, 1.0);
    let fine = blowdown(0.005, 1.0);
    let reversed = blowdown(0.005, -1.0);
    assert!(fine.0 < 2.0 && fine.1 < 2.0);
    assert!(coarse.2 > 3.5 * fine.2, "coarse={coarse:?}, fine={fine:?}");
    assert!(fine.2 < 2e-6, "fine={fine:?}");
    assert_eq!(fine, reversed);
}

#[test]
fn g4_finite_plenum_exhaustion_source_domain_and_cancellation_refuse() {
    let gas = gas(0.4, 0.0);
    let tank = Plenum::new(0.3, gas.conserved(1.0, 0.0, 1.0).unwrap()).unwrap();
    for (mass, energy) in [
        (f64::NAN, 0.0),
        (0.0, f64::INFINITY),
        (-0.3, 0.0),
        (-0.31, 0.0),
        (0.0, -10.0),
        (0.0, 1e6),
    ] {
        assert!(
            tank.exchange_candidate(&gas, mass, energy, || false)
                .is_err()
        );
    }
    let invalid = Plenum::new(
        0.3,
        Conserved {
            energy: -1.0,
            ..tank.state()
        },
    )
    .unwrap();
    assert!(
        invalid
            .exchange_candidate(&gas, 0.0, 3.0, || false)
            .is_err()
    );
    let mut count = 0;
    let candidate = tank
        .exchange_candidate(&gas, 0.01, 0.1, || {
            count += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        tank.exchange_candidate(&gas, 0.01, 0.1, || {
            calls += 1;
            calls == count
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(
        tank.exchange_candidate(&gas, 0.01, 0.1, || false).unwrap(),
        candidate
    );
}
