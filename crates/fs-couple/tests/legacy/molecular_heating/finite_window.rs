use super::*;
use fs_couple::molecular_heating::finite_window::{
    MolecularGasWindowHeating, MolecularWindowError, MolecularWindowLimits,
};

fn prepared<'a>(
    responses: &'a [PreparedMolecularResponse],
    quadrature: &'a [HarmonicVolumeQuadrature],
    gas: &'a IsochoricGasParcel,
    law: &'a RelaxationAbsorption,
    rates: &[&[C64]],
) -> MolecularGasWindowHeating<'a> {
    spectrum(responses, quadrature, gas, law)
        .prepare_coherent_window(
            rates,
            "original acoustic seconds",
            MolecularWindowLimits::default(),
            || false,
        )
        .unwrap()
}

fn phase(omega: f64, time: f64) -> C64 {
    C64::new((omega * time).cos(), -(omega * time).sin())
}

struct TemperatureOracle {
    source: fs_material::relaxation::RelaxationState,
    omega: Vec<f64>,
    pressures: Vec<[C64; 2]>,
    quadrature: Vec<HarmonicVolumeQuadrature>,
}
impl TemperatureOracle {
    fn new(responses: &[PreparedMolecularResponse], rates: &[&[C64]]) -> Self {
        let mut pressures = Vec::new();
        for (response, row) in responses.iter().zip(rates) {
            let mut fields = [HarmonicField::default(); 2];
            response
                .harmonic()
                .evaluate_coherent(row, &mut fields, || false)
                .unwrap();
            pressures.push(fields.map(|f| f.pressure_pa));
        }
        Self {
            source: responses[0].source().clone(),
            omega: responses
                .iter()
                .map(|r| r.harmonic().angular_frequency_rad_s())
                .collect(),
            pressures,
            quadrature: quadrature().to_vec(),
        }
    }

    // Independent caloric reservoir equations, with ORIGINAL molecular times,
    // not the eliminated Maxwell poles used by the product implementation:
    // Cv_inf*T' + sum C_j*T'_j = R*T0*epsilon',
    // T'_j=(T-T_j)/tau_j and p'=p0*epsilon+rho*R*T.
    // Their local balance is p*epsilon' = dE/dt + Q, where
    // E=p0*epsilon²/2 + rho*(Cv_inf*T²+sum C_j*T_j²)/(2*T0).
    fn at(&self, time: f64) -> [f64; 3] {
        let gas = self.source.gas();
        let mut result = [0.0; 3];
        for point in &self.quadrature {
            let mut p = 0.0;
            let mut strain = 0.0;
            let mut strain_rate = 0.0;
            let mut temperature = 0.0;
            let mut internal = vec![0.0; self.source.processes().len()];
            for (bin, &omega) in self.omega.iter().enumerate() {
                let pressure = self.pressures[bin][point.receiver];
                let caloric = self.source.caloric_response(omega).unwrap();
                let epsilon = pressure / caloric.bulk_modulus_pa;
                let theta = (pressure / caloric.cp_j_kg_k)
                    .scale(gas.temperature * gas.specific_gas_constant / gas.pressure);
                let e = phase(omega, time);
                p += (pressure * e).re;
                strain += (epsilon * e).re;
                strain_rate += (C64::new(0.0, -omega) * epsilon * e).re;
                temperature += (theta * e).re;
                for (j, reservoir) in self.source.processes().iter().enumerate() {
                    internal[j] +=
                        (theta / C64::new(1.0, -omega * reservoir.relaxation_time_s) * e).re;
                }
            }
            let mut heat = 0.0;
            let mut storage = self.source.frozen_cv_j_kg_k() * temperature * temperature;
            for (j, reservoir) in self.source.processes().iter().enumerate() {
                let mismatch = temperature - internal[j];
                heat += gas.density * reservoir.heat_capacity_j_kg_k
                    / (gas.temperature * reservoir.relaxation_time_s)
                    * mismatch
                    * mismatch;
                storage += reservoir.heat_capacity_j_kg_k * internal[j] * internal[j];
            }
            storage = 0.5 * gas.pressure * strain * strain
                + 0.5 * gas.density / gas.temperature * storage;
            for (sum, value) in result.iter_mut().zip([heat, p * strain_rate, storage]) {
                *sum += point.volume_m3 * value;
            }
        }
        result
    }

    fn integrate(&self, start: f64, duration: f64, intervals: usize) -> [f64; 2] {
        let mut result = [0.0; 2];
        for i in 0..=intervals {
            let values = self.at(start + duration * i as f64 / intervals as f64);
            let w = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for j in 0..2 {
                result[j] += w * values[j];
            }
        }
        result.map(|v| v * duration / (3.0 * intervals as f64))
    }
}

#[test]
fn g1_single_tone_keeps_carrier_phase_and_recovers_original_cycle_mean() {
    let gas = gas(300.0, VOLUME);
    let law = law(true, 1.0, 500.0);
    // At this lower frequency the original separated volume receivers retain
    // aligned carrier heating. At 1000 rad/s their physical phase difference
    // legitimately cancels much of the spatially integrated carrier ripple.
    let omega = 300.0;
    let responses = [response_at(&gas, &law, omega)];
    let q = quadrature();
    let rates = [C64::new(0.7, -0.2), C64::new(-0.1, 0.3)];
    let heating = prepared(&responses, &q, &gas, &law, &[&rates]);
    let oracle = TemperatureOracle::new(&responses, &[&rates]);
    let period = core::f64::consts::TAU / omega;
    let original_mean = direct_power(&responses[0], &q, &rates);
    close(heating.spectrum().dc_dd().to_f64(), original_mean, 2e-12);
    close(
        heating
            .heat_load(0.00023, period, || false)
            .unwrap()
            .mean_absorbed_power_w(),
        original_mean,
        2e-12,
    );
    let mut largest_carrier_effect: f64 = 0.0;
    for start in [0.0, period / 8.0, period / 4.0, 3.0 * period / 8.0] {
        let duration = period / 80.0;
        let load = heating.heat_load(start, duration, || false).unwrap();
        let expected = oracle.integrate(start, duration, 256)[0];
        close(load.absorbed_energy_j(), expected, 2e-11);
        largest_carrier_effect =
            largest_carrier_effect.max((load.mean_absorbed_power_w() / original_mean - 1.0).abs());
    }
    assert!(largest_carrier_effect > 0.5);
}

#[test]
fn g2_two_reservoir_beats_obey_independent_temperature_energy_balance() {
    let gas = gas(300.0, VOLUME);
    let data = "declared 1 A 2 vibration 1 100 1 power 101.325 300 0 200 500 vibration 1 1600 1 power 1013.25 300 0 200 500";
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let responses = [
        response_at(&gas, &law, 900.0),
        response_at(&gas, &law, 1080.0),
    ];
    let q = quadrature();
    let a = [C64::new(0.7, 0.2), C64::new(-0.1, 0.3)];
    let b = [C64::new(-0.5, 0.3), C64::new(0.2, -0.1)];
    let rates: [&[C64]; 2] = [&a, &b];
    let heating = prepared(&responses, &q, &gas, &law, &rates);
    let oracle = TemperatureOracle::new(&responses, &rates);
    let (start, duration) = (0.0012, 0.0031);
    let load = heating.heat_load(start, duration, || false).unwrap();
    let [heat, work] = oracle.integrate(start, duration, 2048);
    let stored = oracle.at(start + duration)[2] - oracle.at(start)[2];
    close(load.absorbed_energy_j(), heat, 2e-10);
    assert!((work - stored - heat).abs() < 2e-10 * (work.abs() + stored.abs() + heat));
    assert!(
        stored.abs() > 0.01 * heat,
        "reservoir storage must matter to this interval"
    );
    let dc = responses
        .iter()
        .zip(rates)
        .map(|(r, v)| direct_power(r, &q, v))
        .sum::<f64>();
    close(heating.spectrum().dc_dd().to_f64(), dc, 2e-11);
    assert!((load.mean_absorbed_power_w() / dc - 1.0).abs() > 0.1);
}

#[test]
fn g1_exact_source_cancellation_and_small_coherent_residual_survive() {
    let gas = gas(300.0, VOLUME);
    let law = law(true, 1.0, 500.0);
    let responses = [response(&gas, &law)];
    let q = quadrature();
    let cancelled = [C64::new(1e80, -2e80), C64::new(-1e80, 2e80)];
    let heating = prepared(&responses, &q, &gas, &law, &[&cancelled]);
    assert_eq!(
        heating
            .heat_load(0.00017, 0.00003, || false)
            .unwrap()
            .absorbed_energy_j(),
        0.0
    );
    assert!(heating.spectrum().tones().is_empty());
    let epsilon = 2.0_f64.powi(-26);
    let near = [C64::from_re(1.0 + epsilon), C64::from_re(-1.0)];
    let small = [C64::from_re(epsilon), C64::ZERO];
    let near = prepared(&responses, &q, &gas, &law, &[&near]);
    let small = prepared(&responses, &q, &gas, &law, &[&small]);
    let a = near
        .heat_load(0.00017, 0.00003, || false)
        .unwrap()
        .absorbed_energy_j();
    let b = small
        .heat_load(0.00017, 0.00003, || false)
        .unwrap()
        .absorbed_energy_j();
    assert!(a > 0.0);
    close(a, b, 2e-12);
}

#[test]
fn g1_clock_rotation_and_resumed_intervals_preserve_heat() {
    let gas = gas(300.0, VOLUME);
    let law = law(true, 1.0, 500.0);
    let responses = [
        response_at(&gas, &law, 700.0),
        response_at(&gas, &law, 1170.0),
    ];
    let q = quadrature();
    let rates = [
        [C64::new(0.8, 0.1), C64::new(0.2, -0.3)],
        [C64::new(-0.4, 0.2), C64::new(0.1, 0.3)],
    ];
    let heating = prepared(&responses, &q, &gas, &law, &[&rates[0], &rates[1]]);
    let shift = 0.217;
    let moved: Vec<Vec<C64>> = responses
        .iter()
        .zip(rates)
        .map(|(r, row)| {
            let rotation = phase(r.harmonic().angular_frequency_rad_s(), -shift);
            row.iter().map(|&v| v * rotation).collect()
        })
        .collect();
    let moved = prepared(&responses, &q, &gas, &law, &[&moved[0], &moved[1]]);
    let (start, duration, first) = (0.0013, 0.0071, 0.0029);
    let whole = heating
        .heat_load(start, duration, || false)
        .unwrap()
        .absorbed_energy_j();
    let parts = heating
        .heat_load(start, first, || false)
        .unwrap()
        .absorbed_energy_j()
        + heating
            .heat_load(start + first, duration - first, || false)
            .unwrap()
            .absorbed_energy_j();
    close(parts, whole, 2e-12);
    close(
        moved
            .heat_load(start + shift, duration, || false)
            .unwrap()
            .absorbed_energy_j(),
        whole,
        2e-12,
    );
}

#[test]
fn g1_finite_interval_heat_enters_original_nonlinear_gas_once() {
    let gas = gas(300.0, VOLUME);
    let original = gas.clone();
    let law = law(true, 1.0, 500.0);
    let responses = [
        response_at(&gas, &law, 700.0),
        response_at(&gas, &law, 1170.0),
    ];
    let q = quadrature();
    let rates = [C64::new(0.4, -0.2), C64::new(0.2, 0.1)];
    let unit = prepared(&responses, &q, &gas, &law, &[&rates, &rates]);
    let (start, duration) = (0.0017, 100.0031);
    let energy = energy_for_increment(&gas, 4.0);
    let factor = (energy
        / unit
            .heat_load(start, duration, || false)
            .unwrap()
            .absorbed_energy_j())
    .sqrt();
    let rates = rates.map(|z| z.scale(factor));
    let heating = prepared(&responses, &q, &gas, &law, &[&rates, &rates]);
    let load = heating.heat_load(start, duration, || false).unwrap();
    let step = load.heat_gas(controls(), 5.0, || false).unwrap();
    close(load.absorbed_energy_j(), energy, 2e-12);
    assert!((step.gas().state().state().temperature - 304.0).abs() < 1e-8);
    close(step.energy_balance().stored_energy_change_j, energy, 2e-10);
    assert_eq!(step.gas().mass_kg(), gas.mass_kg());
    assert_eq!(step.gas().volume_m3(), gas.volume_m3());
    assert_eq!(step.gas().model(), gas.model());
    assert_eq!(
        step.next_relaxation(),
        &law.evaluate(&RelaxingMedium::from_thermochemical(step.gas().state()).unwrap())
            .unwrap()
    );
    assert!(matches!(
        load.heat_gas(controls(), 2.0, || false),
        Err(MolecularWindowError::Heating(
            MolecularHeatingError::TemperatureIncrement { .. }
        ))
    ));
    assert!(load.heat_gas(controls(), 5.0, || true).is_err());
    assert_eq!(gas, original);
}

#[test]
fn g0_geometry_budget_clock_and_cancellation_refuse_atomically() {
    let gas = gas(300.0, VOLUME);
    let original = gas.clone();
    let law = law(true, 1.0, 500.0);
    let responses = [
        response_at(&gas, &law, 700.0),
        response_at(&gas, &law, 1170.0),
    ];
    let q = quadrature();
    let rates = [C64::ONE, C64::ZERO];
    let changed = [
        HarmonicVolumeQuadrature {
            receiver: 0,
            volume_m3: 0.4 * VOLUME,
        },
        HarmonicVolumeQuadrature {
            receiver: 1,
            volume_m3: 0.6 * VOLUME,
        },
    ];
    let banks = vec![
        responses[0]
            .prepare_absorption(&q, limits(), || false)
            .unwrap(),
        responses[1]
            .prepare_absorption(&changed, limits(), || false)
            .unwrap(),
    ];
    let mismatch = MolecularGasSpectrumHeating::new(banks, &gas, &law, 2, || false).unwrap();
    assert!(
        mismatch
            .prepare_coherent_window(
                &[&rates, &rates],
                "clock",
                MolecularWindowLimits::default(),
                || false
            )
            .is_err()
    );
    for (clock, budget) in [
        ("", MolecularWindowLimits::default()),
        (
            "clock",
            MolecularWindowLimits {
                maximum_products: 1,
                ..MolecularWindowLimits::default()
            },
        ),
        (
            "clock",
            MolecularWindowLimits {
                maximum_workspace_bytes: 1,
                ..MolecularWindowLimits::default()
            },
        ),
    ] {
        assert!(
            spectrum(&responses, &q, &gas, &law)
                .prepare_coherent_window(&[&rates, &rates], clock, budget, || false)
                .is_err()
        );
    }
    let mut calls = 0;
    assert!(
        spectrum(&responses, &q, &gas, &law)
            .prepare_coherent_window(
                &[&rates, &rates],
                "clock",
                MolecularWindowLimits::default(),
                || {
                    calls += 1;
                    calls == 12
                }
            )
            .is_err()
    );
    let heating = prepared(&responses, &q, &gas, &law, &[&rates, &rates]);
    assert!(heating.heat_load(0.0, 1.0, || true).is_err());
    let tiny = [C64::from_re(1e-200), C64::ZERO];
    assert!(
        spectrum(&responses, &q, &gas, &law)
            .prepare_coherent_window(
                &[&tiny, &tiny],
                "clock",
                MolecularWindowLimits::default(),
                || false
            )
            .is_err()
    );
    assert!(heating.heat_load(f64::NAN, 1.0, || false).is_err());
    assert!(heating.heat_load(0.0, 0.0, || false).is_err());
    assert!(
        heating
            .heat_load(0.0, 0.01, || false)
            .unwrap()
            .absorbed_energy_j()
            > 0.0
    );
    assert_eq!(gas, original);
}

#[test]
fn g0_no_relaxing_reservoirs_means_no_molecular_heat() {
    let gas = gas(300.0, VOLUME);
    let law = law(false, 1.0, 500.0);
    let responses = [response(&gas, &law)];
    let q = quadrature();
    let rates = [C64::new(10.0, -20.0), C64::new(3.0, 1.0)];
    let heating = prepared(&responses, &q, &gas, &law, &[&rates]);
    let load = heating.heat_load(0.00017, 0.00003, || false).unwrap();
    assert_eq!(load.absorbed_energy_j(), 0.0);
    assert!(heating.spectrum().tones().is_empty());
    assert_eq!(
        load.heat_gas(controls(), 1.0, || false).unwrap().gas(),
        &gas
    );
}
