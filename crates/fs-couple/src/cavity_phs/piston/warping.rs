//! G1 small-signal map of the actual source-backed gas/momentum owner.
use super::tests::{controls, drive, gas};
use super::*;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_math::det;

const ENERGY_BUDGET_J: f64 = 1e-10;

fn phase_controls(gas_mass_kg: f64) -> PistonSolve {
    let mut solve = controls();
    solve.momentum_tolerance_n_s = 1e-14;
    solve.energy_tolerance_j = ENERGY_BUDGET_J;
    // The material closure takes J/kg, while the outer piston gate takes J.
    // Reserve one eighth of that unchanged outer budget for material energy:
    // m*epsilon_u=epsilon_E/8. Midpoint's own quarter-budget inversions then
    // get epsilon_E/(32*m), rather than a dimensionally unrelated 2.5e-11 J/kg.
    solve.gas.energy_tolerance_j_kg = ENERGY_BUDGET_J / (8.0 * gas_mass_kg);
    solve
}

fn advance(omega: f64, dt: f64, coordinate: usize, signed_amplitude: f64) -> [f64; 2] {
    let mut parcel = gas(0.75, 0.0);
    let initial_volume = parcel.volume_m3();
    let area = 0.01;
    let medium = parcel.state().state();
    let mass = medium.density * medium.sound_speed * medium.sound_speed * area * area
        / (initial_volume * omega * omega);
    let exterior_pressure = medium.pressure;
    let q = if coordinate == 0 {
        signed_amplitude
    } else {
        0.0
    };
    let v = if coordinate == 1 {
        omega * signed_amplitude
    } else {
        0.0
    };
    if q != 0.0 {
        parcel
            .set_volume_isentropic(
                initial_volume + area * q,
                GasEntropySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    entropy_tolerance_j_kg_k: 1e-11,
                    maximum_iterations: 128,
                },
                || false,
            )
            .unwrap();
    }
    let mut model = GasPiston::new(parcel, mass, area, v).unwrap();
    let gas_mass = model.gas().mass_kg();
    let solve = phase_controls(gas_mass);
    let report = model
        .step(
            PistonInput {
                duration_s: dt,
                exterior_pressure_pa: exterior_pressure,
                ..drive(0.0)
            },
            solve,
            || false,
        )
        .unwrap();
    assert!(report.total_energy_residual_j.abs() <= ENERGY_BUDGET_J);
    assert!(report.gas.energy_residual_j.abs() <= ENERGY_BUDGET_J / 8.0);
    assert!(gas_mass * report.gas.midpoint_energy_residual_j_kg.abs() <= ENERGY_BUDGET_J / 8.0);
    [
        (model.gas().volume_m3() - initial_volume) / area,
        model.velocity_m_s() / omega,
    ]
}

#[test]
fn phase_energy_budget_resolves_material_ulps_in_joules() {
    let parcel = gas(0.75, 0.0);
    let u = parcel.state().specific_internal_energy_j_kg();
    let quantum = u.next_up() - u;
    let m = parcel.mass_kg();
    let solve = phase_controls(m);
    // The previous request became 2.5e-11 J/kg inside midpoint, below one
    // representable increment of the actual NASA-9 specific internal energy.
    // Its reported 2^-33 J/kg refusal is about two such increments, not a
    // 1e-10 J defect in this approximately 0.66 g parcel.
    assert!(quantum > 0.25e-10);
    assert!(2.0 * m * quantum < ENERGY_BUDGET_J / 8.0);
    assert!(0.25 * solve.gas.energy_tolerance_j_kg > quantum);
    assert_eq!(solve.energy_tolerance_j, ENERGY_BUDGET_J);
    // Exercise the real source-backed pressure/energy solve at all amplitudes;
    // the phase test and every physical residual gate still use their old bands.
    let omega = core::f64::consts::TAU * 1_000.0;
    for amplitude in [1e-4, 5e-5, 2.5e-5] {
        for coordinate in 0..2 {
            let x = advance(omega, 0.02 / omega, coordinate, amplitude);
            assert!(x.into_iter().all(f64::is_finite));
        }
    }
}

fn phase(omega: f64, dt: f64, amplitude: f64) -> f64 {
    let columns: Vec<_> = (0..2)
        .map(|coordinate| {
            let plus = advance(omega, dt, coordinate, amplitude);
            let minus = advance(omega, dt, coordinate, -amplitude);
            [
                (plus[0] - minus[0]) / (2.0 * amplitude),
                (plus[1] - minus[1]) / (2.0 * amplitude),
            ]
        })
        .collect();
    let re = 0.5 * (columns[0][0] + columns[1][1]);
    let determinant = columns[0][0] * columns[1][1] - columns[1][0] * columns[0][1];
    assert!(determinant > re * re);
    det::atan2(det::sqrt(determinant - re * re), re).abs()
}

#[test]
fn piston_small_signal_warping_matches_midpoint_on_amplitude_ladder() {
    let omega = core::f64::consts::TAU * 1_000.0;
    for z in [0.01, 0.02, 0.03] {
        let expected = 2.0 * det::atan(0.5 * z);
        let mut observed = Vec::new();
        // q/(V/A) ranges from 1e-3 to 2.5e-4: small enough for
        // linearization, large enough to resolve that axis above material ULPs.
        for amplitude in [1e-4, 5e-5, 2.5e-5] {
            let measured = phase(omega, z / omega, amplitude);
            let defect = 1200.0 * det::ln(measured / expected) / core::f64::consts::LN_2;
            observed.push(defect.abs());
            assert!((1200.0 * det::ln(measured / z) / core::f64::consts::LN_2).abs() <= 1.0);
            println!(
                "{{\"path\":\"cavity_phs::piston::GasPiston\",\"integrator\":\"midpoint gas/momentum\",\"omega_hz\":1000,\"sample_rate_hz\":{},\"omega_dt\":{},\"amplitude_m\":{},\"measured_cents\":{},\"predicted_cents\":{},\"linearization_error_cents\":{},\"tolerance_cents\":0.01,\"verdict\":\"resolved-band\"}}",
                omega / z,
                z,
                amplitude,
                1200.0 * det::ln(measured / z) / core::f64::consts::LN_2,
                1200.0 * det::ln(expected / z) / core::f64::consts::LN_2,
                defect
            );
        }
        assert!(
            observed[2] < 0.01,
            "source-backed linear limit differs by {} cents",
            observed[2]
        );
        assert!(
            observed[0] > observed[2],
            "amplitude axis must affect the linearization error"
        );
    }
}

#[test]
fn piston_phase_band_cap_and_cap_plus_one_are_transactional() {
    let parcel = gas(0.75, 0.0);
    let initial = parcel.state().state();
    let area = 0.01;
    let omega = core::f64::consts::TAU * 1_000.0;
    let mass = initial.density * initial.sound_speed.powi(2) * area * area
        / (parcel.volume_m3() * omega * omega);
    let pressure = initial.pressure;
    let mut model = GasPiston::new(parcel, mass, area, 0.0).unwrap();
    let cap = fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
    model
        .step(
            PistonInput {
                duration_s: cap / omega * (1.0 - 1e-12),
                exterior_pressure_pa: pressure,
                ..drive(0.0)
            },
            controls(),
            || false,
        )
        .unwrap();
    let before = model.clone();
    assert!(matches!(
        model.step(
            PistonInput {
                duration_s: cap / omega * (1.0 + 1e-10),
                exterior_pressure_pa: pressure,
                ..drive(0.0)
            },
            controls(),
            || false
        ),
        Err(PistonError::TimeResolution(_))
    ));
    assert_eq!(model, before);
}
