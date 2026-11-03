//! G0/G1/G3/G4: one shared stream, original calorics and finite vector moments.
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::drift::exchange::{
    DriftingEffusiveExchange, sample as sample_exchange,
};

const DT: f64 = 0.125;

fn pair(
    left: &ClosedGasParcel,
    right: &ClosedGasParcel,
    interface: &EffusiveDriftInterface,
    velocities: [[f64; 3]; 2],
    normal: [f64; 3],
) -> DriftingEffusiveExchange {
    sample_exchange(left, right, interface, AREA, velocities, normal, DT, || {
        false
    })
    .unwrap()
}

#[test]
fn rest_pair_preserves_original_beams_and_full_pressure_transfer_at_zero_net_mass() {
    let left = gas(0.0, 400.0);
    let right = left.clone();
    let interface = qualified(&left, [0.8, 0.3]);
    let normal = [0.6, 0.8, 0.0];
    let exchange = pair(&left, &right, &interface, [[0.0; 3]; 2], normal);
    let original = left
        .effusive_flux_to(
            &right,
            interface.momentum_interface().interface(),
            AREA,
            || false,
        )
        .unwrap();
    assert_eq!(exchange.duration_s, DT);
    assert_eq!(exchange.species.len(), 2);
    for ((row, original), (molar_mass, fraction, tau)) in exchange
        .species
        .iter()
        .zip(&original.species)
        .zip([(0.004, 0.25, 0.8), (0.028, 0.75, 0.3)])
    {
        let mass = AREA
            * DT
            * tau
            * fraction
            * 100_000.0
            * (molar_mass / (2.0 * core::f64::consts::PI * R * 400.0)).sqrt();
        let impulse = 0.5 * AREA * DT * tau * fraction * 100_000.0;
        near(row.mass_kg[0], mass);
        near(row.mass_kg[1], mass);
        near(row.mass_kg[0], original.left_to_right.mass_rate_kg_s * DT);
        near(row.mass_kg[1], original.right_to_left.mass_rate_kg_s * DT);
        assert_eq!(row.net_mass_left_to_right_kg, 0.0);
        assert_eq!(row.relative_energy_into_gas_j, [0.0; 2]);
        vector_near(row.momentum_relative_kg_m_s[0], normal.map(|n| impulse * n));
        vector_near(
            row.momentum_relative_kg_m_s[1],
            normal.map(|n| -impulse * n),
        );
        for d in 0..3 {
            // Both half-spaces transmit normal stress even when species mass
            // and heat balance. Remaining control-volume wall forces are separate.
            near(
                row.momentum_relative_kg_m_s[0][d] - row.momentum_relative_kg_m_s[1][d],
                2.0 * impulse * normal[d],
            );
        }
        near(row.energy_j[0], original.left_to_right.energy_rate_w * DT);
        near(row.energy_j[1], original.right_to_left.energy_rate_w * DT);
    }
}

#[test]
fn thermal_and_tangential_energy_crosswalk_retains_original_reference_independence() {
    let temperatures = [320.0, 900.0];
    let velocities = [[0.0, 150.0, -45.0], [0.0, -80.0, 95.0]];
    let left = gas(0.0, temperatures[0]);
    let right = gas(0.0, temperatures[1]);
    let interface = qualified(&left, [0.8, 0.3]);
    let exchange = pair(&left, &right, &interface, velocities, NORMAL);
    let shifted_left = gas(1e24, temperatures[0]);
    let shifted_right = gas(1e24, temperatures[1]);
    let shifted = pair(
        &shifted_left,
        &shifted_right,
        &qualified(&shifted_left, [0.8, 0.3]),
        velocities,
        NORMAL,
    );
    for ((row, shifted), (molar_mass, fraction, tau, cp_r, slope)) in exchange
        .species
        .iter()
        .zip(&shifted.species)
        .zip([(0.004, 0.25, 0.8, 2.5, 0.0), (0.028, 0.75, 0.3, 3.5, 0.001)])
    {
        let mass = temperatures.map(|t| {
            AREA * DT
                * tau
                * fraction
                * 100_000.0
                * (molar_mass / (2.0 * core::f64::consts::PI * R * t)).sqrt()
        });
        let h = temperatures.map(|t| R * (cp_r * t + 0.5 * slope * t * t) / molar_mass);
        let a: [f64; 2] = core::array::from_fn(|side| {
            -0.5 * R * temperatures[side] / molar_mass
                + 0.5 * (velocities[side][1].powi(2) + velocities[side][2].powi(2))
        });
        let dh = R / molar_mass
            * (temperatures[1] - temperatures[0])
            * (cp_r + 0.5 * slope * (temperatures[0] + temperatures[1]));
        for side in 0..2 {
            near(row.mass_kg[side], mass[side]);
            near(row.energy_j[side], mass[side] * (h[side] + a[side]));
            let direction = if side == 0 { 1.0 } else { -1.0 };
            vector_near(
                row.momentum_relative_kg_m_s[side],
                [
                    direction * 0.5 * AREA * DT * tau * fraction * 100_000.0,
                    mass[side] * velocities[side][1],
                    mass[side] * velocities[side][2],
                ],
            );
            near(
                shifted.energy_j[side],
                row.energy_j[side] + row.mass_kg[side] * 1e24,
            );
        }
        near(row.net_mass_left_to_right_kg, mass[0] - mass[1]);
        near(
            row.relative_energy_into_gas_j[0],
            -mass[0] * a[0] + mass[1] * (dh + a[1]),
        );
        near(
            row.relative_energy_into_gas_j[1],
            mass[0] * (-dh + a[0]) - mass[1] * a[1],
        );
        // The two ORIGINAL-temperature energy charts close after their source
        // enthalpy crosswalk, before either finite gas is reconstructed.
        let closure = row.relative_energy_into_gas_j[0]
            + row.relative_energy_into_gas_j[1]
            + row.net_mass_left_to_right_kg * dh;
        let scale = row.relative_energy_into_gas_j[0].abs()
            + row.relative_energy_into_gas_j[1].abs()
            + (row.net_mass_left_to_right_kg * dh).abs();
        assert!(closure.abs() <= 64.0 * f64::EPSILON * scale);
        assert_eq!(row.mass_kg, shifted.mass_kg);
        assert_eq!(
            row.net_mass_left_to_right_kg,
            shifted.net_mass_left_to_right_kg
        );
        assert_eq!(
            row.momentum_relative_kg_m_s,
            shifted.momentum_relative_kg_m_s
        );
        assert_eq!(
            row.relative_energy_into_gas_j,
            shifted.relative_energy_into_gas_j
        );
    }
}

#[test]
fn tiny_signed_transport_and_two_positive_rare_tail_beams_survive_cancellation() {
    let left = gas(0.0, 400.0);
    let right = left.clone();
    let interface = qualified(&left, [0.0, 0.6]);
    let sigma = (R * 400.0 / 0.028).sqrt();
    let density = 0.75 * 100_000.0 * 0.028 / (R * 400.0);
    for sign in [-1.0, 1.0] {
        let speed = sign * sigma * 1e-20;
        let exchange = pair(
            &left,
            &right,
            &interface,
            [[speed, 0.0, 0.0], [0.0; 3]],
            NORMAL,
        );
        let row = exchange
            .species
            .iter()
            .find(|row| row.species.as_str() == "B")
            .unwrap();
        near(
            row.net_mass_left_to_right_kg,
            0.5 * AREA * DT * 0.6 * density * speed,
        );
        assert!(row.net_mass_left_to_right_kg.abs() < f64::EPSILON * row.mass_kg[0]);
    }
    let tails = pair(
        &left,
        &right,
        &interface,
        [[-10.0 * sigma, 0.0, 0.0], [10.0 * sigma, 0.0, 0.0]],
        NORMAL,
    );
    let row = tails
        .species
        .iter()
        .find(|row| row.species.as_str() == "B")
        .unwrap();
    // Independent positive half-space integrals at a=-10; a baseline-minus-
    // baseline reconstruction cannot preserve this actual two-way exchange.
    let mass = AREA * DT * 0.6 * density * sigma * 7.474_560_254_589_328e-25;
    let impulse = mass * sigma * 0.194_383_033_412_553_3;
    assert!(row.mass_kg.iter().all(|mass| *mass > 0.0));
    near(row.mass_kg[0], mass);
    near(row.mass_kg[1], mass);
    near(row.momentum_relative_kg_m_s[0][0], impulse);
    near(row.momentum_relative_kg_m_s[1][0], -impulse);
    assert_eq!(row.net_mass_left_to_right_kg, 0.0);
    assert_eq!(row.relative_energy_into_gas_j, [0.0; 2]);
}

#[test]
fn pair_orientation_scaling_and_atomic_refusals_keep_original_sources() {
    let left = gas(0.0, 450.0);
    let right = gas(0.0, 800.0);
    let saved = [left.clone(), right.clone()];
    let interface = qualified(&left, [0.8, 0.3]);
    let velocities = [[120.0, -60.0, 80.0], [-170.0, 45.0, 15.0]];
    let normal = [0.6, 0.8, 0.0];
    let original = pair(&left, &right, &interface, velocities, normal);
    let swapped = pair(
        &right,
        &left,
        &interface,
        [velocities[1], velocities[0]],
        normal.map(|n| -n),
    );
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let rotated = pair(
        &left,
        &right,
        &interface,
        velocities.map(rotate),
        rotate(normal),
    );
    let scaled = sample_exchange(
        &left,
        &right,
        &interface,
        2.0 * AREA,
        velocities,
        normal,
        3.0 * DT,
        || false,
    )
    .unwrap();
    for (((row, swapped), rotated), scaled) in original
        .species
        .iter()
        .zip(&swapped.species)
        .zip(&rotated.species)
        .zip(&scaled.species)
    {
        near(
            swapped.net_mass_left_to_right_kg,
            -row.net_mass_left_to_right_kg,
        );
        near(
            rotated.net_mass_left_to_right_kg,
            row.net_mass_left_to_right_kg,
        );
        near(
            scaled.net_mass_left_to_right_kg,
            6.0 * row.net_mass_left_to_right_kg,
        );
        for side in 0..2 {
            near(swapped.mass_kg[side], row.mass_kg[1 - side]);
            vector_near(
                swapped.momentum_relative_kg_m_s[side],
                row.momentum_relative_kg_m_s[1 - side],
            );
            near(swapped.energy_j[side], row.energy_j[1 - side]);
            near(
                swapped.relative_energy_into_gas_j[side],
                row.relative_energy_into_gas_j[1 - side],
            );
            near(rotated.mass_kg[side], row.mass_kg[side]);
            vector_near(
                rotated.momentum_relative_kg_m_s[side],
                rotate(row.momentum_relative_kg_m_s[side]),
            );
            near(
                rotated.relative_energy_into_gas_j[side],
                row.relative_energy_into_gas_j[side],
            );
            near(scaled.mass_kg[side], 6.0 * row.mass_kg[side]);
            vector_near(
                scaled.momentum_relative_kg_m_s[side],
                row.momentum_relative_kg_m_s[side].map(|p| 6.0 * p),
            );
            near(scaled.energy_j[side], 6.0 * row.energy_j[side]);
            near(
                scaled.relative_energy_into_gas_j[side],
                6.0 * row.relative_energy_into_gas_j[side],
            );
        }
    }
    for (area, duration, velocity, normal) in [
        (-AREA, DT, velocities, normal),
        (AREA, 0.0, velocities, normal),
        (AREA, f64::NAN, velocities, normal),
        (AREA, DT, [[200_000.0, 0.0, 0.0], velocities[1]], normal),
        (AREA, DT, velocities, [2.0, 0.0, 0.0]),
    ] {
        assert!(
            sample_exchange(
                &left,
                &right,
                &interface,
                area,
                velocity,
                normal,
                duration,
                || false
            )
            .is_err()
        );
    }
    let foreign = gas(1.0, 800.0);
    assert!(
        sample_exchange(
            &left,
            &foreign,
            &interface,
            AREA,
            velocities,
            normal,
            DT,
            || false
        )
        .is_err()
    );
    let inactive = sample_exchange(
        &left,
        &foreign,
        &interface,
        0.0,
        [[f64::NAN; 3]; 2],
        [0.0; 3],
        DT,
        || false,
    )
    .unwrap();
    assert!(inactive.species.is_empty());
    let mut polls = 0;
    let accepted = sample_exchange(
        &left,
        &right,
        &interface,
        AREA,
        velocities,
        normal,
        DT,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        sample_exchange(
            &left,
            &right,
            &interface,
            AREA,
            velocities,
            normal,
            DT,
            || {
                count += 1;
                count == polls
            }
        ),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(accepted, original);
    assert_eq!([left, right], saved);
}
