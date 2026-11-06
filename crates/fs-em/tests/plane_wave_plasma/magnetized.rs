//! G1/G3/G0: circular Maxwell propagation retains the original Hall transport.
use super::*;
use fs_em::plane_wave::plasma::magnetized::MagnetizedIonizedGasInterface;
use fs_material::ionized_gas::IonizedGasMagneticLimits;

type Vector = [C64; 3];

fn magnetic_limits() -> IonizedGasMagneticLimits {
    IonizedGasMagneticLimits {
        maximum_cyclotron_to_thermal_energy_ratio: 0.01,
    }
}
fn vector_norm(v: Vector) -> f64 {
    v[0].abs().hypot(v[1].abs()).hypot(v[2].abs())
}
fn vector_near(a: Vector, b: Vector, tolerance: f64) {
    let error = vector_norm(std::array::from_fn(|d| a[d] - b[d]));
    assert!(
        error <= tolerance * vector_norm(b).max(1e-100),
        "{a:?} != {b:?}: {error:e}"
    );
}
fn add(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|d| a[d] + b[d])
}
fn times_i(z: C64) -> C64 {
    C64::new(-z.im, z.re)
}
fn norm_sq(v: Vector) -> f64 {
    v.into_iter().map(C64::norm_sq).sum()
}
fn peak(v: Vector) -> f64 {
    let bilinear = v.into_iter().fold(C64::ZERO, |s, z| s + z * z);
    ((norm_sq(v) + bilinear.abs()) * 0.5).sqrt()
}

fn drude_sigma(
    gas: &IonizedGasConductivity,
    omega: f64,
    field_along_normal: f64,
    sign: f64,
) -> C64 {
    gas.carriers().iter().fold(C64::ZERO, |s, carrier| {
        let q = carrier.charge_number() as f64 * fixture::E;
        let mass = carrier.mass_kg();
        let detuning = omega - sign * q * field_along_normal / mass;
        s + C64::from_re(carrier.number_density_per_m3() * q * q / mass)
            / C64::new(carrier.relaxation().rate_per_s(), detuning)
    })
}
fn reference_mode(
    gas: &IonizedGasConductivity,
    omega: f64,
    b: f64,
    sign: f64,
) -> (C64, C64, C64, C64) {
    let sigma = drude_sigma(gas, omega, b, sign);
    let epsilon = C64::new(EPSILON + sigma.im / omega, -sigma.re / omega);
    let k = epsilon.sqrt().scale(omega * MU.sqrt());
    let yi = C64::from_re((EPSILON / MU).sqrt());
    let yt = k.scale(1.0 / (omega * MU));
    let r = (yi - yt) / (yi + yt);
    let t = yi.scale(2.0) / (yi + yt);
    (epsilon, k, r, t)
}

#[test]
fn g1_signed_circular_dispersion_and_resonance_match_independent_drude_fresnel() {
    let gas = gas();
    let harmonic = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 1e12;
    for strength in [0.0_f64, 0.4, 1.0] {
        let b = strength * omega * fixture::ME / fixture::E;
        let transport = harmonic
            .prepare_magnetized(omega, [0.0, 0.0, b], magnetic_limits(), || false)
            .unwrap();
        for normal_sign in [1.0_f64, -1.0] {
            let interface = MagnetizedIonizedGasInterface::try_new(
                DielectricResponse::Lossless(&background),
                &transport,
                &background,
                [0.0, 0.0, normal_sign],
                || false,
            )
            .unwrap();
            assert!(std::ptr::eq(interface.transport(), &transport));
            assert!(std::ptr::eq(interface.transport().source(), &gas));
            assert!(std::ptr::eq(interface.background(), &background));
            for (m, sign) in [1.0_f64, -1.0].into_iter().enumerate() {
                let actual = interface.circular_modes()[m];
                let expected = reference_mode(&gas, omega, normal_sign * b, sign);
                complex_near(
                    interface.circular_conductivities_s_m()[m],
                    drude_sigma(&gas, omega, normal_sign * b, sign),
                    5e-13,
                );
                complex_near(actual.transmitted_permittivity_f_m, expected.0, 3e-12);
                complex_near(
                    actual.transmitted_normal_wavenumber_rad_m,
                    expected.1,
                    3e-12,
                );
                complex_near(actual.reflected_tangential_e, expected.2, 3e-12);
                complex_near(actual.transmitted_tangential_e, expected.3, 3e-12);
                assert!(actual.transmitted_normal_wavenumber_rad_m.re > 0.0);
                assert!(actual.transmitted_normal_wavenumber_rad_m.im < 0.0);
                near(
                    actual.reflected_power_fraction + actual.transmitted_power_fraction,
                    1.0,
                    5e-14,
                );
            }
            if strength == 1.0 {
                let sigma = interface.circular_conductivities_s_m();
                let resonant = if normal_sign > 0.0 { 1 } else { 0 };
                assert!(sigma[resonant].re > 1000.0 * sigma[1 - resonant].re);
            }
        }
    }
}

#[test]
fn g1_zero_bias_recovers_original_te_tm_fields_heat_and_boundary_continuity() {
    let gas = gas();
    let harmonic = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 2e12;
    let transport = harmonic
        .prepare_magnetized(omega, [0.0; 3], magnetic_limits(), || false)
        .unwrap();
    let interface = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &transport,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let original_law =
        IonizedGasWaveResponse::try_new(&harmonic, &background, omega, || false).unwrap();
    let original = PlanarDielectricInterface::try_new(
        &background,
        DielectricResponse::IonizedGas(original_law),
    )
    .unwrap();
    let incident = [C64::new(0.013, -0.009), C64::new(-0.005, 0.011), C64::ZERO];
    let wave = interface.prepare_wave(incident, || false).unwrap();
    let tm = original
        .prepare_ionized_gas_wave(1.0, Polarization::Tm, incident[0], || false)
        .unwrap();
    let te = original
        .prepare_ionized_gas_wave(1.0, Polarization::Te, incident[1], || false)
        .unwrap();
    let k = tm.interface_response().transmitted_normal_wavenumber_rad_m;
    let entrance = wave.field(0.0).unwrap();
    vector_near(
        add(wave.incident_electric_v_m(), wave.reflected_electric_v_m()),
        entrance.electric_v_m,
        3e-13,
    );
    vector_near(
        add(wave.incident_magnetic_a_m(), wave.reflected_magnetic_a_m()),
        entrance.magnetic_a_m,
        3e-13,
    );
    near(
        wave.incident_power_w_m2(),
        wave.reflected_power_w_m2() + wave.transmitted_power_w_m2(),
        3e-13,
    );
    for depth in [0.0_f64, 0.1, 2.0, 900.0].map(|d| d / -k.im) {
        let point = wave.field(depth).unwrap();
        let t = tm.field(depth).unwrap();
        let e = te.field(depth).unwrap();
        vector_near(
            point.electric_v_m,
            add(t.electric_v_m, e.electric_v_m),
            3e-12,
        );
        vector_near(
            point.magnetic_a_m,
            add(t.magnetic_a_m, e.magnetic_a_m),
            3e-12,
        );
        vector_near(
            point.current_density_a_m2,
            add(t.current_density_a_m2, e.current_density_a_m2),
            3e-12,
        );
        near(
            point.mean_heating_w_m3(),
            t.mean_heating_w_m3() + e.mean_heating_w_m3(),
            3e-12,
        );
        near(
            point.mean_drift_kinetic_energy_j_m3(),
            t.mean_drift_kinetic_energy_j_m3() + e.mean_drift_kinetic_energy_j_m3(),
            3e-12,
        );
        near(
            point.absorbed_power_to_depth_w_m2,
            t.absorbed_power_to_depth_w_m2 + e.absorbed_power_to_depth_w_m2,
            3e-12,
        );
    }
}

#[test]
fn g1_coupled_wave_closes_maxwell_flux_and_original_signed_thermal_accounts() {
    let original = gas();
    let gas = original
        .clone()
        .with_coulomb_drag(fixture::coulomb_model(original.state()), || false)
        .unwrap();
    let harmonic = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 1e12;
    let b = 0.4 * omega * fixture::ME / fixture::E;
    let transport = harmonic
        .prepare_magnetized(omega, [0.0, 0.0, b], magnetic_limits(), || false)
        .unwrap();
    let interface = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &transport,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let wave = interface
        .prepare_wave(
            [C64::new(0.011, -0.009), C64::new(-0.007, 0.004), C64::ZERO],
            || false,
        )
        .unwrap();
    let modes = interface.circular_modes();
    let fastest_decay = modes
        .iter()
        .map(|m| -m.transmitted_normal_wavenumber_rad_m.im)
        .fold(0.0_f64, f64::max);
    let fastest_phase = modes
        .iter()
        .map(|m| m.transmitted_normal_wavenumber_rad_m.abs())
        .fold(0.0_f64, f64::max);
    assert!(wave.power_balance_relative_residual() < 1e-10);
    for distance in [0.0_f64, 0.1, 0.7, 3.0] {
        let depth = distance / fastest_decay;
        let point = wave.field(depth).unwrap();
        let local = transport.response(point.electric_v_m, || false).unwrap();
        vector_near(
            point.current_density_a_m2,
            local.current_density_a_m2,
            3e-11,
        );
        vector_near(
            point.complex_poynting_w_m2,
            cross(point.electric_v_m, point.magnetic_a_m.map(C64::conj)).map(|s| s.scale(0.5)),
            3e-12,
        );
        near(
            point.mean_heating_w_m3(),
            local.mean_heating_w_m3.total_power_density_w_m3(),
            3e-11,
        );
        near(
            point.mean_electron_heating_w_m3(),
            local.mean_heating_w_m3.electron_power_density_w_m3(),
            3e-11,
        );
        near(
            point.mean_heavy_heating_w_m3(),
            local.mean_heating_w_m3.heavy_power_density_w_m3(),
            3e-11,
        );
        near(
            point.mean_exchange_correction_to_electrons_w_m3(),
            local
                .mean_heating_w_m3
                .thermal_exchange_correction_to_electrons_w_m3(),
            3e-11,
        );
        near(
            point.mean_drift_kinetic_energy_j_m3(),
            local.mean_drift_kinetic_energy_j_m3,
            3e-11,
        );
        near(
            point.reactive_carrier_power_var_m3(),
            local.reactive_power_var_m3,
            3e-11,
        );
        near(
            point.gyroscopic_reactive_power_var_m3(),
            local.gyroscopic_reactive_power_var_m3,
            3e-11,
        );
        near(
            point.reactive_carrier_power_var_m3(),
            2.0 * omega * point.mean_drift_kinetic_energy_j_m3()
                + point.gyroscopic_reactive_power_var_m3(),
            3e-11,
        );
        assert!(point.mean_magnetic_power_w_m3().abs() < 1e-10 * point.mean_heating_w_m3());
        near(
            point.mean_background_field_energy_j_m3(),
            0.25 * (EPSILON * norm_sq(point.electric_v_m) + MU * norm_sq(point.magnetic_a_m)),
            3e-12,
        );
        near(
            point.complex_poynting_w_m2[2].re + point.absorbed_power_to_depth_w_m2,
            wave.transmitted_power_w_m2(),
            3e-12,
        );
        for carrier in 0..gas.carriers().len() {
            vector_near(
                point.carrier_velocity_m_s(carrier).unwrap(),
                local.carrier_velocity_m_s[carrier],
                3e-11,
            );
        }
        assert!(point.carrier_velocity_m_s(gas.carriers().len()).is_none());
        if distance > 0.0 {
            // Independent centered spatial differences exercise the full
            // superposition; a single scalar effective k would fail here.
            let h = 2e-5 / fastest_phase;
            let left = wave.field(depth - h).unwrap();
            let right = wave.field(depth + h).unwrap();
            let derivative_e: Vector = std::array::from_fn(|d| {
                (right.electric_v_m[d] - left.electric_v_m[d]).scale(0.5 / h)
            });
            let derivative_h: Vector = std::array::from_fn(|d| {
                (right.magnetic_a_m[d] - left.magnetic_a_m[d]).scale(0.5 / h)
            });
            let normal = [C64::ZERO, C64::ZERO, C64::ONE];
            vector_near(
                cross(normal, derivative_e),
                point.magnetic_a_m.map(|h| times_i(h).scale(-omega * MU)),
                2e-8,
            );
            vector_near(
                cross(normal, derivative_h),
                add(
                    point.current_density_a_m2,
                    point
                        .electric_v_m
                        .map(|e| times_i(e).scale(omega * EPSILON)),
                ),
                2e-8,
            );
            near(
                -(right.complex_poynting_w_m2[2].re - left.complex_poynting_w_m2[2].re) / (2.0 * h),
                point.mean_heating_w_m3(),
                2e-8,
            );
        }
    }
    let slowest_decay = modes
        .iter()
        .map(|m| -m.transmitted_normal_wavenumber_rad_m.im)
        .fold(f64::INFINITY, f64::min);
    let opaque = wave.field(900.0 / slowest_decay).unwrap();
    assert_eq!(opaque.electric_v_m, [C64::ZERO; 3]);
    assert_eq!(opaque.mean_heating_w_m3(), 0.0);
    near(
        opaque.absorbed_power_to_depth_w_m2,
        wave.transmitted_power_w_m2(),
        3e-13,
    );
}

fn orientation(field: Vector) -> C64 {
    let q = field[0].norm_sq() - field[1].norm_sq();
    let u = 2.0 * (field[0] * field[1].conj()).re;
    C64::new(q, u).scale(1.0 / q.hypot(u))
}
fn circular_amplitudes(field: Vector) -> [C64; 2] {
    [
        (field[0] - times_i(field[1])).scale(core::f64::consts::FRAC_1_SQRT_2),
        (field[0] + times_i(field[1])).scale(core::f64::consts::FRAC_1_SQRT_2),
    ]
}

#[test]
fn g1_faraday_rotation_and_dichroism_reverse_with_bias_without_refitting() {
    let gas = gas();
    let harmonic = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 4e12;
    let b = 0.5 * omega * fixture::ME / fixture::E;
    let positive = harmonic
        .prepare_magnetized(omega, [0.0, 0.0, b], magnetic_limits(), || false)
        .unwrap();
    let negative = harmonic
        .prepare_magnetized(omega, [0.0, 0.0, -b], magnetic_limits(), || false)
        .unwrap();
    let plus = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &positive,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let minus = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &negative,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let amplitude = 0.02;
    let wave = plus
        .prepare_wave([C64::from_re(amplitude), C64::ZERO, C64::ZERO], || false)
        .unwrap();
    let reversed = minus
        .prepare_wave([C64::from_re(amplitude), C64::ZERO, C64::ZERO], || false)
        .unwrap();
    let first = reference_mode(&gas, omega, b, 1.0);
    let second = reference_mode(&gas, omega, b, -1.0);
    let delta_k = first.1.re - second.1.re;
    assert!(delta_k > 0.0);
    let depth = 0.6 / delta_k;
    let surface = wave.field(0.0).unwrap();
    let point = wave.field(depth).unwrap();
    let phase = |k: C64| {
        let (s, c) = (k.re * depth).sin_cos();
        C64::new(c, -s).scale((k.im * depth).exp())
    };
    let a = first.3 * phase(first.1);
    let b = second.3 * phase(second.1);
    let reference = [
        (a + b).scale(0.5 * amplitude),
        times_i(a - b).scale(0.5 * amplitude),
        C64::ZERO,
    ];
    vector_near(point.electric_v_m, reference, 4e-12);
    let rotated = orientation(point.electric_v_m) / orientation(surface.electric_v_m);
    complex_near(rotated, C64::new(0.6_f64.cos(), 0.6_f64.sin()), 4e-12);
    let surface_circular = circular_amplitudes(surface.electric_v_m);
    let point_circular = circular_amplitudes(point.electric_v_m);
    let ratio = (point_circular[0].abs() / point_circular[1].abs())
        / (surface_circular[0].abs() / surface_circular[1].abs());
    near(ratio, ((first.1.im - second.1.im) * depth).exp(), 4e-12);
    assert!((ratio - 1.0).abs() > 1e-4);
    let opposite = reversed.field(depth).unwrap();
    vector_near(
        opposite.electric_v_m,
        [point.electric_v_m[0], -point.electric_v_m[1], C64::ZERO],
        3e-12,
    );
    near(
        opposite.mean_heating_w_m3(),
        point.mean_heating_w_m3(),
        3e-12,
    );
}

#[test]
fn g3_arbitrary_proper_frame_rotation_preserves_fields_and_material_heat() {
    let gas = gas();
    let harmonic = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 1.7e12;
    let b = 0.3 * omega * fixture::ME / fixture::E;
    let n = [
        1.0 / 14.0_f64.sqrt(),
        2.0 / 14.0_f64.sqrt(),
        3.0 / 14.0_f64.sqrt(),
    ];
    let first = [2.0 / 5.0_f64.sqrt(), -1.0 / 5.0_f64.sqrt(), 0.0];
    let second = [
        3.0 / 70.0_f64.sqrt(),
        6.0 / 70.0_f64.sqrt(),
        -5.0 / 70.0_f64.sqrt(),
    ];
    let rotate = |v: Vector| -> Vector {
        std::array::from_fn(|d| v[0].scale(first[d]) + v[1].scale(second[d]) + v[2].scale(n[d]))
    };
    let original = harmonic
        .prepare_magnetized(omega, [0.0, 0.0, b], magnetic_limits(), || false)
        .unwrap();
    let rotated = harmonic
        .prepare_magnetized(omega, n.map(|v| v * b), magnetic_limits(), || false)
        .unwrap();
    let a = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &original,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let c = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &rotated,
        &background,
        n,
        || false,
    )
    .unwrap();
    let electric = [C64::new(0.013, 0.005), C64::new(-0.007, 0.002), C64::ZERO];
    let first_wave = a.prepare_wave(electric, || false).unwrap();
    let rotated_wave = c.prepare_wave(rotate(electric), || false).unwrap();
    let depth = 0.7 / -a.circular_modes()[0].transmitted_normal_wavenumber_rad_m.im;
    let p = first_wave.field(depth).unwrap();
    let q = rotated_wave.field(depth).unwrap();
    vector_near(q.electric_v_m, rotate(p.electric_v_m), 5e-12);
    vector_near(q.magnetic_a_m, rotate(p.magnetic_a_m), 5e-12);
    vector_near(
        q.current_density_a_m2,
        rotate(p.current_density_a_m2),
        5e-12,
    );
    vector_near(
        q.complex_poynting_w_m2,
        rotate(p.complex_poynting_w_m2),
        5e-12,
    );
    near(q.mean_heating_w_m3(), p.mean_heating_w_m3(), 5e-12);
    near(
        q.absorbed_power_to_depth_w_m2,
        p.absorbed_power_to_depth_w_m2,
        5e-12,
    );
}

#[test]
fn g0_combined_drift_conditions_geometry_and_cancellation_are_not_bypassed_by_modes() {
    let gas = gas();
    let harmonic = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 1e12;
    let zero_b = harmonic
        .prepare_magnetized(omega, [0.0; 3], magnetic_limits(), || false)
        .unwrap();
    let interface = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &zero_b,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let small = 1e-6;
    let reference = interface
        .prepare_wave([C64::from_re(small), C64::ZERO, C64::ZERO], || false)
        .unwrap();
    let maximum_ratio = gas
        .carriers()
        .iter()
        .zip(&reference.surface_response().carrier_velocity_m_s)
        .map(|(carrier, &v)| {
            let limit = (fixture::KB * carrier.temperature_k() / carrier.mass_kg()).sqrt()
                * fixture::limits().maximum_drift_thermal_ratio;
            peak(v) / limit
        })
        .fold(0.0_f64, f64::max);
    let excessive = 1.5 * small / maximum_ratio;
    let t = interface.circular_modes()[0]
        .transmitted_tangential_e
        .scale(0.5 * excessive);
    // Each circular mode individually passes, but their actual linear
    // superposition has twice the peak speed and must refuse.
    for sign in [1.0_f64, -1.0] {
        zero_b
            .response([t, times_i(t).scale(sign), C64::ZERO], || false)
            .unwrap();
    }
    assert!(
        interface
            .prepare_wave([C64::from_re(excessive), C64::ZERO, C64::ZERO], || false)
            .is_err()
    );
    assert!(interface.prepare_wave([C64::ONE; 3], || false).is_err());
    assert!(
        interface
            .prepare_wave([C64::new(f64::NAN, 0.0), C64::ZERO, C64::ZERO], || false)
            .is_err()
    );
    let zero = interface.prepare_wave([C64::ZERO; 3], || false).unwrap();
    let far = zero.field(f64::MAX).unwrap();
    assert_eq!(far.electric_v_m, [C64::ZERO; 3]);
    assert_eq!(far.mean_heating_w_m3(), 0.0);
    assert_eq!(far.absorbed_power_to_depth_w_m2, 0.0);
    for depth in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(zero.field(depth).is_err());
    }
    assert!(matches!(
        interface.prepare_wave([C64::ZERO; 3], || true),
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    ));
    let mut polls = 0;
    assert!(matches!(
        interface.prepare_wave([C64::ZERO; 3], || {
            polls += 1;
            polls >= 8
        }),
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    ));
    let biased = harmonic
        .prepare_magnetized(omega, [0.0, 0.0, 1.0], magnetic_limits(), || false)
        .unwrap();
    for normal in [
        [0.0; 3],
        [1.0, 0.0, 0.0],
        [0.01, 0.0, 1.0],
        [f64::NAN, 0.0, 1.0],
    ] {
        assert!(
            MagnetizedIonizedGasInterface::try_new(
                DielectricResponse::Lossless(&background),
                &biased,
                &background,
                normal,
                || false,
            )
            .is_err()
        );
    }
    let pressure = gas.state().pressure().value();
    for wrong in [
        background_at(6001.0, 12000.0, pressure, true),
        background_at(6000.0, 12001.0, pressure, true),
        background_at(6000.0, 12000.0, pressure * 1.01, true),
        background_at(6000.0, 12000.0, pressure, false),
    ] {
        assert!(
            MagnetizedIonizedGasInterface::try_new(
                DielectricResponse::Lossless(&background),
                &biased,
                &wrong,
                [0.0, 0.0, 1.0],
                || false,
            )
            .is_err()
        );
    }
    assert!(matches!(
        MagnetizedIonizedGasInterface::try_new(
            DielectricResponse::Lossless(&background),
            &biased,
            &background,
            [0.0, 0.0, 1.0],
            || true,
        ),
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    ));
}
