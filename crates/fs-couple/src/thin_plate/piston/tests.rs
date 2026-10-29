use super::*;
use fs_material::gas::GasSpec;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};

fn air(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}

fn band() -> [f64; 2] {
    [40.0, 0.4 * core::f64::consts::PI * 48_000.0]
}

fn source(radius: f64) -> BaffledPistonRadiation {
    BaffledPistonRadiation::new(
        radius,
        &air(101_325.0),
        band(),
        257,
        128,
        &CancelGate::default(),
    )
    .unwrap()
}

#[test]
fn g1_continuous_piston_resolves_small_resistance_on_independent_frequency_grids() {
    for radius in [0.01, 0.05, 0.15] {
        let source = source(radius);
        let [lo, hi] = band();
        let (mut complex, mut resistance) = (0.0_f64, 0.0_f64);
        // Neither grid was used for fit coefficients or its stopping decision.
        // Log spacing explicitly exercises the weak low-frequency resistance
        // that a single full-band peak normalization used to hide.
        for i in 0..=512 {
            for omega in [
                lo + (hi - lo) * i as f64 / 512.0,
                (lo * (hi / lo).powf(i as f64 / 512.0)).min(hi),
            ] {
                let expected = fs_phs::baffled_piston_impedance(
                    source.gas().density,
                    source.gas().sound_speed,
                    radius,
                    omega,
                )
                .unwrap();
                let actual = source.impedance(omega).unwrap();
                let reference = C64::new(expected.0, expected.1);
                complex = complex.max((actual - reference).abs() / reference.abs());
                resistance = resistance.max((actual.re - reference.re).abs() / reference.re);
                assert!(actual.re > 0.0);
            }
        }
        assert!(complex < 0.01, "radius {radius}, complex error {complex}");
        assert!(
            resistance < 0.02,
            "radius {radius}, resistance error {resistance}"
        );
        assert!(source.report().complex_pointwise_relative_error <= 0.02);
        assert!(source.report().resistance_pointwise_relative_error <= 0.02);
        assert_eq!(source.report().poles, source.model().poles.len());
    }
}

fn physical_rhs(
    model: &Model,
    mass: f64,
    stiffness: f64,
    damping: f64,
    area: f64,
    force: f64,
    x: &[f64],
) -> Vec<f64> {
    let n = 2 + 2 * model.poles.len();
    let mut dx = vec![0.0; n + 2];
    let mut pressure = 0.0;
    let mut heat = damping * x[1] * x[1];
    for (i, pole) in model.poles.iter().enumerate() {
        let q = 2 + 2 * i;
        let rate = 2.0 * pole.zeta * pole.omega;
        pressure += pole.coupling[0] * x[q + 1];
        dx[q] = x[q + 1];
        dx[q + 1] =
            pole.coupling[0] * area * x[1] - pole.omega * pole.omega * x[q] - rate * x[q + 1];
        heat += rate * x[q + 1] * x[q + 1];
    }
    dx[0] = x[1];
    dx[1] = (force - stiffness * x[0] - damping * x[1] - area * pressure) / mass;
    dx[n] = force * x[1];
    dx[n + 1] = heat;
    dx
}

fn rk4_step(x: &mut Vec<f64>, dt: f64, rhs: impl Fn(&[f64]) -> Vec<f64>) {
    let k1 = rhs(x);
    let stage = |k: &[f64], fraction| {
        x.iter()
            .zip(k)
            .map(|(x, k)| x + fraction * dt * k)
            .collect::<Vec<_>>()
    };
    let k2 = rhs(&stage(&k1, 0.5));
    let k3 = rhs(&stage(&k2, 0.5));
    let k4 = rhs(&stage(&k3, 1.0));
    for i in 0..x.len() {
        x[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
}

#[test]
fn g1_reciprocal_piston_memory_matches_physical_velocity_equations_and_heat() {
    let source = source(0.01);
    let (mass, omega, damping, force, dt) = (2e-5, 600.0_f64, 4.8e-5, 2e-4, 1e-5);
    let stiffness = mass * omega * omega;
    let face = core::f64::consts::PI * source.radius_m().powi(2);
    let mut positive = None;
    for area in [face, -face, 0.0] {
        let host = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
        let system = host
            .with_parallel_load_columns(source.port_hamiltonian().unwrap(), &[0.0, area])
            .unwrap();
        let n = system.state_dim();
        let mut x = vec![0.0; n];
        x[0] = 1e-6;
        x[1] = mass * 0.004;
        let initial = system.hamiltonian(&x);
        let flow = PreparedExactZoh::new_owned(
            system,
            dt,
            LinearSplitBudget {
                max_modes: n,
                max_omega_dt: core::f64::consts::PI,
            },
        )
        .unwrap();
        let mut reference = vec![0.0; n + 2];
        reference[0] = x[0];
        reference[1] = x[1] / mass;
        let (mut work, mut heat) = (0.0, 0.0);
        for _ in 0..20 {
            let rec = flow.step(&x, &[force]).unwrap();
            work += rec.supplied;
            heat += rec.dissipated;
            assert!(rec.dissipated >= 0.0);
            x = rec.x;
            for _ in 0..200 {
                rk4_step(&mut reference, dt / 200.0, |state| {
                    physical_rhs(source.model(), mass, stiffness, damping, area, force, state)
                });
            }
        }
        let mut difference =
            mass * ((omega * (x[0] - reference[0])).powi(2) + (x[1] / mass - reference[1]).powi(2));
        for (i, pole) in source.model().poles.iter().enumerate() {
            let q = 2 + 2 * i;
            difference += (pole.omega * (x[q] - reference[q])).powi(2)
                + (x[q + 1] - reference[q + 1]).powi(2);
        }
        let scale = initial + work.abs();
        assert!(
            difference.sqrt() < 2e-8 * scale.sqrt(),
            "physical state error {difference:e}"
        );
        assert!((work - reference[n]).abs() < 2e-8 * scale, "force work");
        assert!(
            (heat - reference[n + 1]).abs() < 2e-8 * scale,
            "independent viscous/radiation heat"
        );
        assert!((flow.system().hamiltonian(&x) - initial + heat - work).abs() < 2e-11 * scale);
        if area > 0.0 {
            positive = Some(x.clone());
        }
        if area < 0.0 {
            let first = positive.as_ref().unwrap();
            assert!((x[0] - first[0]).abs() < 1e-13);
            assert!((x[1] - first[1]).abs() < 1e-13);
            for i in 2..n {
                assert!((x[i] + first[i]).abs() < 1e-12);
            }
        }
        if area == 0.0 {
            assert!(x[2..].iter().all(|&v| v == 0.0));
        }
    }
}

fn loaded_pole(
    model: &Model,
    mass: f64,
    stiffness: f64,
    damping: f64,
    area: f64,
    initial: C64,
) -> C64 {
    let mut s = initial;
    for _ in 0..16 {
        let (mut z, mut dz) = (C64::ZERO, C64::ZERO);
        for pole in &model.poles {
            let w2 = pole.omega * pole.omega;
            let den = s * s + s.scale(2.0 * pole.zeta * pole.omega) + C64::from_re(w2);
            let gain = pole.coupling[0] * pole.coupling[0];
            z = z + s.scale(gain) / den;
            dz = dz + (C64::from_re(w2) - s * s).scale(gain) / (den * den);
        }
        let residual = (s * s).scale(mass)
            + s.scale(damping)
            + C64::from_re(stiffness)
            + (s * z).scale(area * area);
        let derivative =
            s.scale(2.0 * mass) + C64::from_re(damping) + (z + s * dz).scale(area * area);
        let correction = residual / derivative;
        s = s - correction;
        if correction.abs() < 1e-12 * s.abs() {
            return s;
        }
    }
    panic!("loaded physical pole did not converge");
}

#[test]
fn g1_compact_added_mass_and_decay_remain_physical_after_closing_the_loop() {
    let source = source(0.01);
    let gas = source.gas();
    let radius = source.radius_m();
    let area = core::f64::consts::PI * radius * radius;
    let mass = 2e-5;
    let stiffness = mass * 600.0_f64.powi(2);
    let damping = 4.8e-5;
    // Independent small-ka expansion of the exact uniform piston:
    // m_air=8 rho a^3/3 and R_mech=rho*S*a^2*omega^2/(2c).
    let inertia = mass + 8.0 * gas.density * radius.powi(3) / 3.0;
    let omega = (stiffness / inertia).sqrt();
    let resistance = gas.density * area * radius.powi(2) * omega.powi(2) / (2.0 * gas.sound_speed);
    let decay = (damping + resistance) / (2.0 * inertia);
    let actual = loaded_pole(
        source.model(),
        mass,
        stiffness,
        damping,
        area,
        C64::new(-decay, omega),
    );
    assert!(omega * radius / gas.sound_speed < 0.02);
    assert!(
        actual.im < 600.0,
        "air loading must lower the mechanical frequency"
    );
    assert!(
        (actual.im / omega - 1.0).abs() < 0.01,
        "added-mass pole {actual:?}"
    );
    assert!(
        (-actual.re / decay - 1.0).abs() < 0.025,
        "loaded decay {actual:?} vs {decay}"
    );
}

#[test]
fn g0_source_band_capacity_density_and_held_out_power_errors_are_explicit() {
    let first = source(0.01);
    let gas = air(101_325.0);
    let count = first.model().poles.len();
    assert!(
        BaffledPistonRadiation::new(0.01, &gas, band(), 257, count - 1, &CancelGate::default())
            .is_err()
    );
    let second = BaffledPistonRadiation::new(
        0.01,
        &air(202_650.0),
        band(),
        257,
        count,
        &CancelGate::default(),
    )
    .unwrap();
    for omega in [40.0, 100.0, 600.0, 5000.0, band()[1]] {
        let a = first.impedance(omega).unwrap();
        let b = second.impedance(omega).unwrap();
        assert!((b - a.scale(2.0)).abs() < 1e-10 * a.abs());
    }
    for omega in [0.0, 39.0, band()[1] * 1.01, f64::NAN] {
        assert!(first.impedance(omega).is_err());
    }
    let mut bad = gas;
    bad.density *= 2.0;
    assert!(
        BaffledPistonRadiation::new(0.01, &bad, band(), 257, 128, &CancelGate::default()).is_err()
    );
    for points in [0, 32, 66, 259] {
        assert!(
            BaffledPistonRadiation::new(0.01, &gas, band(), points, 128, &CancelGate::default())
                .is_err()
        );
    }
    let gate = CancelGate::default();
    gate.request();
    assert!(matches!(
        BaffledPistonRadiation::new(0.01, &gas, band(), 257, 128, &gate),
        Err(AcousticRealizeError::Cancelled)
    ));

    let [lo, hi] = band();
    let omega: Vec<_> = (0..257)
        .map(|i| lo + (hi - lo) * (i as f64 / 256.0).powi(2))
        .collect();
    let mut data: Vec<_> = omega
        .iter()
        .map(|&w| {
            let (r, x) =
                fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, 0.01, w).unwrap();
            vec![C64::new(r, x)]
        })
        .collect();
    // The first held-out resistance is tiny compared with the full-band peak.
    // Doubling it must still refuse even though reactive/peak gates can hide it.
    data[1][0].re *= 2.0;
    assert!(fit::fit_resolved_relative(&omega, &data, 1, (hi - lo) / 8.0).is_err());
}
