//! G1/G3: physical force transmission through a coherent finite fluid layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::{
    PlanarFluidLayer, PreparedLayeredFullWaveForce, PreparedLayeredFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{
    FullWaveForceLimits, FullWaveMonopoleLimits, PreparedFullWaveForce,
};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn layer(states: [FluidAcousticState; 3], h: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(
        states[0],
        states[1],
        states[2],
        [0.0; 3],
        [0.0, 0.0, 1.0],
        h,
    )
    .unwrap()
}
fn limits() -> FullWaveForceLimits {
    FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-12,
        absolute_velocity_tolerance_m_s_per_n: 1e-15,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn unit(axis: usize) -> [C64; 3] {
    core::array::from_fn(|i| if i == axis { C64::ONE } else { C64::ZERO })
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}

#[test]
fn all_tensor_channels_match_independent_full_spectrum_layer_references() {
    // Direct-q mpmath1.3.0 tanh-sinh at35/50digits, all critical points and
    // independently located narrow resonances split explicitly. No production
    // branch substitutions, denominator normalizers or GK nodes are reused.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
            700.0,
            0.025,
            0.018,
            0.023,
            0.04,
            [
                C64::new(10.857610346824266, 2.3839116809306973),
                C64::new(18.542617167535287, 2.346263024447859),
                C64::new(0.007089381850185072, -0.0102157582941132),
                C64::new(0.007528002537988684, -0.03428655469965853),
                C64::new(0.0001124268294976462, 0.03869171279617209),
                C64::new(-0.0016445457915486332, 0.040774202723601436),
                C64::new(0.006107285404237076, 0.03139600095473126),
            ],
        ),
        (
            [fluid(2.4, 600.0), fluid(1.4, 450.0), fluid(1.2, 300.0)],
            1500.0,
            0.04,
            0.07,
            0.06,
            0.12,
            [
                C64::new(-2.8711776978013483, -1.2817106329065895),
                C64::new(-2.393318580502536, -1.6422340434305886),
                C64::new(-0.00335301700050349, -0.001080800255017953),
                C64::new(-0.0009444008464839139, 0.002115565384723577),
                C64::new(-0.0014156896999471312, -0.00413132043369219),
                C64::new(-0.006100072765487454, -0.004157845975997523),
                C64::new(-0.005733623629066523, -0.003604170395914106),
            ],
        ),
        (
            [fluid(1.2, 300.0), fluid(1000.0, 600.0), fluid(1.2, 300.0)],
            1500.0,
            0.12,
            0.31,
            0.08,
            0.09,
            [
                C64::new(-0.0004953064688956653, -0.0005653188947342842),
                C64::new(-0.004461891532939172, -0.002757446473588368),
                C64::new(-6.376371125261489e-7, 2.175286683057699e-7),
                C64::new(-5.553907192605443e-7, 4.866078572230694e-7),
                C64::new(-1.1866322997503804e-6, -1.6236418769501898e-6),
                C64::new(-1.1866322997503804e-6, -1.6236418769501898e-6),
                C64::new(-1.1796336042834913e-5, -7.978655620422733e-6),
            ],
        ),
    ];
    for (states, hz, a, h, b, r, expected) in cases {
        let prepared = PreparedLayeredFullWaveForce::new(
            layer(states, h),
            hz,
            [0.0, 0.0, -a],
            [r, 0.0, h + b],
            limits(),
            &mut || false,
        )
        .unwrap();
        let x = prepared.evaluate(unit(0)).unwrap();
        let y = prepared.evaluate(unit(1)).unwrap();
        let z = prepared.evaluate(unit(2)).unwrap();
        let actual = [
            x.pressure_pa,
            z.pressure_pa,
            x.particle_velocity_m_s[0],
            y.particle_velocity_m_s[1],
            z.particle_velocity_m_s[0],
            x.particle_velocity_m_s[2],
            z.particle_velocity_m_s[2],
        ];
        for (a, e) in actual.into_iter().zip(expected) {
            close(a, e, 4e-8);
        }
        assert_eq!(y.pressure_pa, C64::ZERO);
        assert_eq!(x.particle_velocity_m_s[1], C64::ZERO);
        assert_eq!(z.particle_velocity_m_s[1], C64::ZERO);
    }
}

#[test]
fn collapsed_layer_matches_original_interface_and_equal_speed_images_keep_density() {
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.035, 0.0, 0.018];
    let first = fluid(1.2, 300.0);
    let last = fluid(2.4, 450.0);
    let f = [
        C64::new(0.2, -0.1),
        C64::new(-0.15, 0.3),
        C64::new(0.1, 0.25),
    ];
    // A zero-thickness slow layer has no trapped-mode problem or material path.
    let collapsed = PreparedLayeredFullWaveForce::new(
        layer([first, fluid(1000.0, 100.0), last], 0.0),
        700.0,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap()
    .evaluate(f)
    .unwrap();
    let interface = PlanarFluidInterface::new(first, last, [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
    let direct =
        PreparedFullWaveForce::new(interface, 700.0, source, receiver, limits(), &mut || false)
            .unwrap()
            .evaluate(f)
            .unwrap();
    close(collapsed.pressure_pa, direct.pressure_pa, 3e-14);
    for axis in 0..3 {
        close(
            collapsed.particle_velocity_m_s[axis],
            direct.particle_velocity_m_s[axis],
            3e-14,
        );
    }

    // When all sound speeds agree, R/T are angle-independent. The exact
    // transmitted Green function is a convergent series of spherical images;
    // force directions remain original because each return has two reflections.
    let states = [fluid(1.2, 300.0), fluid(5.0, 300.0), fluid(2.4, 300.0)];
    let h = 0.08;
    let receiver = [0.035, 0.01, h + 0.018];
    let prepared = PreparedLayeredFullWaveForce::new(
        layer(states, h),
        700.0,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = prepared.evaluate(f).unwrap();
    let transmission = (2.0 * 5.0 / (1.2 + 5.0)) * (2.0 * 2.4 / (5.0 + 2.4));
    let ratio = ((1.2 - 5.0) / (1.2 + 5.0)) * ((2.4 - 5.0) / (2.4 + 5.0));
    let mut expected = [C64::ZERO; 4];
    let mut weight = transmission;
    for image in 0..64 {
        let position = [source[0], source[1], source[2] - 2.0 * f64::from(image) * h];
        let sample = PreparedCompactSourceObserver::new(states[0], 700.0, receiver, &[position], 1)
            .unwrap()
            .evaluate(&[CompactSourcePhasor::force_dipole(f)])
            .unwrap();
        expected[0] = expected[0] + sample.pressure_pa.scale(weight);
        for axis in 0..3 {
            expected[axis + 1] =
                expected[axis + 1] + sample.particle_velocity_m_s[axis].scale(weight * 1.2 / 2.4);
        }
        weight *= ratio;
    }
    close(field.pressure_pa, expected[0], 3e-8);
    for axis in 0..3 {
        close(field.particle_velocity_m_s[axis], expected[axis + 1], 3e-8);
    }
}

#[test]
fn layer_force_uses_source_derivatives_and_the_receiving_fluid_momentum() {
    let states = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let slab = layer(states, 0.018);
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.0, 0.041];
    let hz = 700.0;
    let mono_limits = FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-9,
        absolute_velocity_tolerance_per_m2: 1e-12,
        relative_tolerance: 1e-12,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    };
    let prepared =
        PreparedLayeredFullWaveForce::new(slab, hz, source, receiver, limits(), &mut || false)
            .unwrap();
    let step = 2e-6;
    for axis in [0, 2] {
        let mut plus = source;
        plus[axis] += step;
        let mut minus = source;
        minus[axis] -= step;
        let q = |s| {
            PreparedLayeredFullWaveMonopole::new(slab, hz, s, receiver, mono_limits, &mut || false)
                .unwrap()
                .evaluate_harmonic_line(C64::ONE)
                .unwrap()
        };
        let a = q(plus);
        let b = q(minus);
        let conversion = C64::new(0.0, 1.0 / (TAU * hz * 1.2 * 2.0 * step));
        let f = prepared.evaluate(unit(axis)).unwrap();
        close(
            f.pressure_pa,
            (a.pressure_pa - b.pressure_pa) * conversion,
            2e-5,
        );
        for row in 0..3 {
            if row != 1 {
                close(
                    f.particle_velocity_m_s[row],
                    (a.particle_velocity_m_s[row] - b.particle_velocity_m_s[row]) * conversion,
                    2e-5,
                );
            }
        }
    }
    let f = [
        C64::new(0.2, 0.1),
        C64::new(-0.1, 0.3),
        C64::new(0.15, -0.2),
    ];
    let actual = prepared.evaluate(f).unwrap();
    for axis in 0..3 {
        let mut plus = receiver;
        plus[axis] += step;
        let mut minus = receiver;
        minus[axis] -= step;
        let p = |r| {
            PreparedLayeredFullWaveForce::new(slab, hz, source, r, limits(), &mut || false)
                .unwrap()
                .evaluate_harmonic_line(f)
                .unwrap()
                .pressure_pa
        };
        let gradient = (p(plus) - p(minus)).scale(1.0 / (2.0 * step));
        close(
            actual.particle_velocity_m_s[axis],
            gradient * C64::new(0.0, -1.0 / (TAU * hz * 1.8)),
            2e-5,
        );
    }
}

#[test]
fn axis_transverse_motion_clock_and_bounded_refusals_are_retained() {
    let states = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let slab = layer(states, 0.018);
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.0, 0.0, 0.041];
    let prepared =
        PreparedLayeredFullWaveForce::new(slab, 700.0, source, receiver, limits(), &mut || false)
            .unwrap();
    let field = prepared.evaluate(unit(0)).unwrap();
    assert_eq!(field.pressure_pa, C64::ZERO);
    assert!(field.particle_velocity_m_s[0].abs() > 0.0);
    assert_eq!(field.particle_velocity_m_s[1], C64::ZERO);
    assert_eq!(field.particle_velocity_m_s[2], C64::ZERO);
    let y = prepared.evaluate(unit(1)).unwrap();
    close(
        y.particle_velocity_m_s[1],
        field.particle_velocity_m_s[0],
        2e-14,
    );
    let line = prepared.evaluate_harmonic_line(unit(0)).unwrap();
    assert_eq!(line.frequency_hz, 700.0);
    assert_eq!(line.pressure_pa, field.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    assert_eq!(
        prepared.evaluate([C64::ZERO; 3]).unwrap().pressure_pa,
        C64::ZERO
    );
    assert!(prepared.evaluate([C64::new(f64::NAN, 0.0); 3]).is_err());
    assert!(
        PreparedLayeredFullWaveForce::new(
            slab,
            700.0,
            source,
            receiver,
            FullWaveForceLimits {
                maximum_evaluations: 30,
                ..limits()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedLayeredFullWaveForce::new(slab, 700.0, source, receiver, limits(), &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PreparedLayeredFullWaveForce::new(slab, 0.0, source, receiver, limits(), &mut || false)
            .is_err()
    );
    assert!(
        PreparedLayeredFullWaveForce::new(slab, 700.0, source, [0.0; 3], limits(), &mut || false)
            .is_err()
    );
    assert!(
        PlanarFluidLayer::new(
            states[0],
            fluid(1.0, 200.0),
            states[2],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.018
        )
        .is_err()
    );
}
