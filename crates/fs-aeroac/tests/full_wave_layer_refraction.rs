//! G1/G3: full spherical sources, coherent layer returns and critical spectra.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::{
    PlanarFluidLayer, PreparedLayeredFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveMonopoleLimits, PreparedFullWaveMonopole};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-10,
        absolute_velocity_tolerance_per_m2: 1e-13,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn slab(states: [FluidAcousticState; 3], h: f64) -> PlanarFluidLayer {
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
fn prepare(
    layer: PlanarFluidLayer,
    hz: f64,
    source: [f64; 3],
    receiver: [f64; 3],
) -> PreparedLayeredFullWaveMonopole {
    PreparedLayeredFullWaveMonopole::new(layer, hz, source, receiver, limits(), &mut || false)
        .unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}, relative tolerance={tolerance:e}"
    );
}

#[test]
fn absent_identical_and_matched_fluid_layers_keep_original_interface_limits() {
    let media = [fluid(1.2, 300.0), fluid(0.8, 200.0), fluid(2.4, 450.0)];
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.035, 0.01, 0.041];
    // An absent layer cannot trap, even if its unused state would be slower.
    let absent = prepare(slab(media, 0.0), 700.0, source, receiver);
    let interface =
        PlanarFluidInterface::new(media[0], media[2], [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
    let original =
        PreparedFullWaveMonopole::new(interface, 700.0, source, receiver, limits(), &mut || false)
            .unwrap();
    let q = C64::new(0.003, -0.001);
    let expected = original.evaluate(q).unwrap();
    let actual = absent.evaluate(q).unwrap();
    assert_eq!(actual.pressure_pa, expected.pressure_pa);
    assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
    assert!(absent.diagnostics().reduced_to_single_interface);
    assert_eq!(
        absent.diagnostics().evaluations,
        original.diagnostics().evaluations
    );
    let homogeneous = prepare(slab([media[0]; 3], 0.018), 700.0, source, receiver);
    let reference = PreparedCompactSourceObserver::new(media[0], 700.0, receiver, &[source], 1)
        .unwrap()
        .evaluate(&[CompactSourcePhasor::monopole(q)])
        .unwrap();
    let actual = homogeneous.evaluate(q).unwrap();
    close(actual.pressure_pa, reference.pressure_pa, 5e-14);
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            reference.particle_velocity_m_s[axis],
            5e-14,
        );
    }
    assert_eq!(homogeneous.diagnostics().evaluations, 0);
    assert!(homogeneous.diagnostics().reduced_to_single_interface);
    for (middle, interface_z) in [(media[0], 0.018), (media[2], 0.0)] {
        let actual = prepare(
            slab([media[0], middle, media[2]], 0.018),
            700.0,
            source,
            receiver,
        )
        .evaluate(q)
        .unwrap();
        let interface =
            PlanarFluidInterface::new(media[0], media[2], [0.0, 0.0, interface_z], [0.0, 0.0, 1.0])
                .unwrap();
        let expected = PreparedFullWaveMonopole::new(
            interface,
            700.0,
            source,
            receiver,
            limits(),
            &mut || false,
        )
        .unwrap()
        .evaluate(q)
        .unwrap();
        close(actual.pressure_pa, expected.pressure_pa, 2e-8);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
                2e-8,
            );
        }
    }
    let line = homogeneous.evaluate_harmonic_line(q).unwrap();
    assert_eq!(line.frequency_hz, 700.0);
    assert_eq!(
        line.pressure_pa,
        homogeneous.evaluate(q).unwrap().pressure_pa
    );
    assert_eq!(homogeneous.source_position_m(), source);
    assert_eq!(homogeneous.observer_position_m(), receiver);
    assert_eq!(homogeneous.layer().thickness_m(), 0.018);
}

#[test]
fn independent_full_spectrum_references_include_tunneling_and_narrow_resonances() {
    // Independent mpmath1.3.0 direct-q tanh-sinh at35/50 decimal digits.
    // All critical points and narrow resonances are independently split; the
    // production branch maps, scaled denominator and GK nodes are not reused.
    // Precision differences: p<1e-32, velocities<3e-35 for these parameters.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
            700.0,
            0.025,
            0.018,
            0.023,
            0.04,
            [
                C64::new(5165.123_941_891_204, -3398.782_564_943_966_5),
                C64::new(7.238_406_897_882_844, 1.589_274_453_953_798_2),
                C64::new(11.441_749_131_107_962, 3.529_384_156_172_694_4),
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
                C64::new(-5368.763_670_800_917, -633.890_403_890_175_9),
                C64::new(-5.742_355_395_602_697, -2.563_421_265_813_179),
                C64::new(-13.159_480_090_731_36, -2.370_123_092_281_613),
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
                C64::new(-1.683_421_558_150_960_2, -0.942_777_674_698_248_6),
                C64::new(-0.000_495_306_468_895_665_3, -0.000_565_318_894_734_284_2),
                C64::new(-0.004_461_891_532_939_172, -0.002_757_446_473_588_368),
            ],
        ),
    ];
    for (media, hz, a, h, b, r, expected) in cases {
        let observer = prepare(slab(media, h), hz, [0.0, 0.0, -a], [r, 0.0, h + b]);
        let field = observer.evaluate(C64::ONE).unwrap();
        close(field.pressure_pa, expected[0], 8e-9);
        close(field.particle_velocity_m_s[0], expected[1], 8e-9);
        close(field.particle_velocity_m_s[2], expected[2], 8e-9);
        assert_eq!(field.particle_velocity_m_s[1], C64::ZERO);
        assert!(!observer.diagnostics().reduced_to_single_interface);
        assert!(observer.diagnostics().evaluations > 0);
        assert!(observer.diagnostics().evaluations <= limits().maximum_evaluations);
        assert!(observer.diagnostics().maximum_bessel_argument < 3000.0);
    }
}

#[test]
fn equal_speeds_keep_coherent_density_layer_images_and_common_critical_circle() {
    let rho1 = 1.2;
    let rho3 = 2.4;
    let c = 350.0;
    let hz = 1100.0;
    let (a, h, b, r) = (0.03, 0.08, 0.04, 0.06);
    for rho2 in [5.0, 1.8] {
        let layer = slab([fluid(rho1, c), fluid(rho2, c), fluid(rho3, c)], h);
        let observer = prepare(layer, hz, [0.0, 0.0, -a], [r, 0.0, h + b]);
        let actual = observer.evaluate(C64::ONE).unwrap();
        // Here all angular Fresnel coefficients are constant, so the exact
        // coherent Green function is an independently constructed image series.
        let transmission = 2.0 * rho2 / (rho1 + rho2) * 2.0 * rho3 / (rho2 + rho3);
        let reflection = (rho1 - rho2) / (rho1 + rho2) * (rho3 - rho2) / (rho3 + rho2);
        let omega = TAU * hz;
        let mut weight = transmission;
        let mut expected = [C64::ZERO; 3];
        for order in 0..80 {
            let normal = a + b + (2.0 * f64::from(order) + 1.0) * h;
            let distance = f64::hypot(r, normal);
            let (sin, cos) = (omega * distance / c).sin_cos();
            let green = C64::new(cos, sin).scale(1.0 / (4.0 * PI * distance));
            expected[0] = expected[0] + C64::new(green.im, -green.re).scale(omega * rho1 * weight);
            let radial = green * C64::new(1.0 / distance, -omega / c).scale(weight * rho1 / rho3);
            expected[1] = expected[1] + radial.scale(r / distance);
            expected[2] = expected[2] + radial.scale(normal / distance);
            weight *= reflection;
        }
        assert!(weight.abs() < 1e-50);
        close(actual.pressure_pa, expected[0], 2e-8);
        close(actual.particle_velocity_m_s[0], expected[1], 2e-8);
        close(actual.particle_velocity_m_s[2], expected[2], 2e-8);
        assert!(!observer.diagnostics().reduced_to_single_interface);
    }
}

#[test]
fn original_receiving_momentum_reciprocity_and_rotation_survive_the_layer() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let h = 0.018;
    let hz = 700.0;
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.015, 0.041];
    let layer = slab(media, h);
    let actual = prepare(layer, hz, source, receiver)
        .evaluate(C64::ONE)
        .unwrap();
    let step = 2e-6;
    for axis in 0..3 {
        let mut plus = receiver;
        let mut minus = receiver;
        plus[axis] += step;
        minus[axis] -= step;
        let gradient = (prepare(layer, hz, source, plus)
            .evaluate(C64::ONE)
            .unwrap()
            .pressure_pa
            - prepare(layer, hz, source, minus)
                .evaluate(C64::ONE)
                .unwrap()
                .pressure_pa)
            .scale(0.5 / step);
        let velocity =
            C64::new(gradient.im, -gradient.re).scale(1.0 / (TAU * hz * media[2].density_kg_m3()));
        close(actual.particle_velocity_m_s[axis], velocity, 3e-6);
    }
    let reverse = PlanarFluidLayer::new(
        media[2],
        media[1],
        media[0],
        [0.0, 0.0, h],
        [0.0, 0.0, -1.0],
        h,
    )
    .unwrap();
    let reciprocal = prepare(reverse, hz, receiver, source)
        .evaluate(C64::ONE)
        .unwrap();
    close(actual.pressure_pa, reciprocal.pressure_pa, 3e-9);
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let shifted = |v: [f64; 3]| {
        let v = rotate(v);
        [v[0] + 0.5, v[1] - 0.25, v[2] + 0.75]
    };
    let rotated_layer = PlanarFluidLayer::new(
        media[0],
        media[1],
        media[2],
        shifted([0.0; 3]),
        rotate([0.0, 0.0, 1.0]),
        h,
    )
    .unwrap();
    let rotated = prepare(rotated_layer, hz, shifted(source), shifted(receiver))
        .evaluate(C64::ONE)
        .unwrap();
    close(rotated.pressure_pa, actual.pressure_pa, 3e-9);
    for (axis, original) in [2, 0, 1].into_iter().enumerate() {
        close(
            rotated.particle_velocity_m_s[axis],
            actual.particle_velocity_m_s[original],
            3e-9,
        );
    }
}

#[test]
fn near_zero_thickness_and_numerical_refinement_retain_original_limits() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.0, 0.041];
    let zero = prepare(slab(media, 0.0), 700.0, source, receiver)
        .evaluate(C64::ONE)
        .unwrap();
    let thin = prepare(slab(media, 1e-6), 700.0, source, receiver)
        .evaluate(C64::ONE)
        .unwrap();
    let thinner = prepare(slab(media, 5e-7), 700.0, source, receiver)
        .evaluate(C64::ONE)
        .unwrap();
    assert!(
        (thinner.pressure_pa - zero.pressure_pa).abs()
            < 0.51 * (thin.pressure_pa - zero.pressure_pa).abs()
    );
    let layer = slab(media, 0.018);
    let coarse = PreparedLayeredFullWaveMonopole::new(
        layer,
        700.0,
        source,
        receiver,
        FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 1e-5,
            absolute_velocity_tolerance_per_m2: 1e-8,
            relative_tolerance: 1e-7,
            ..limits()
        },
        &mut || false,
    )
    .unwrap();
    let fine = prepare(layer, 700.0, source, receiver);
    assert!(
        fine.diagnostics().pressure_tail_bound_pa_per_m3_s
            < coarse.diagnostics().pressure_tail_bound_pa_per_m3_s
    );
    assert!(
        fine.diagnostics().transverse_wavenumber_cutoff_rad_m
            > coarse.diagnostics().transverse_wavenumber_cutoff_rad_m
    );
    close(
        fine.evaluate(C64::ONE).unwrap().pressure_pa,
        coarse.evaluate(C64::ONE).unwrap().pressure_pa,
        3e-6,
    );
    assert!(
        fine.diagnostics().pressure_tail_bound_pa_per_m3_s
            <= limits().absolute_pressure_tolerance_pa_per_m3_s * 0.25
    );
    let velocity = fine.diagnostics().velocity_tail_bound_per_m2;
    assert!(velocity[0].hypot(velocity[1]) <= limits().absolute_velocity_tolerance_per_m2 * 0.25);
}

#[test]
fn trapping_geometry_unfinished_work_and_cancellation_never_publish_partial_fields() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    for h in [-0.1, f64::NAN, f64::INFINITY] {
        assert!(
            PlanarFluidLayer::new(media[0], media[1], media[2], [0.0; 3], [0.0, 0.0, 1.0], h)
                .is_err()
        );
    }
    assert!(
        PlanarFluidLayer::new(
            media[0],
            fluid(2.0, 200.0),
            media[2],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.01
        )
        .is_err()
    );
    let layer = slab(media, 0.018);
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.0, 0.041];
    for bad in [
        source,
        [0.0, 0.0, 0.01],
        [0.0, 0.0, 0.018],
        [f64::NAN, 0.0, 0.04],
    ] {
        assert!(
            PreparedLayeredFullWaveMonopole::new(layer, 700.0, source, bad, limits(), &mut || {
                false
            })
            .is_err()
        );
    }
    for limit in [
        FullWaveMonopoleLimits {
            maximum_evaluations: 44,
            ..limits()
        },
        FullWaveMonopoleLimits {
            maximum_depth: 1,
            ..limits()
        },
        FullWaveMonopoleLimits {
            relative_tolerance: 1.0,
            ..limits()
        },
    ] {
        assert!(
            PreparedLayeredFullWaveMonopole::new(
                layer,
                700.0,
                source,
                receiver,
                limit,
                &mut || false
            )
            .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        PreparedLayeredFullWaveMonopole::new(layer, 700.0, source, receiver, limits(), &mut || {
            polls += 1;
            polls == 100
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let observer = prepare(layer, 700.0, source, receiver);
    assert!(observer.evaluate(C64::new(f64::NAN, 0.0)).is_err());
    assert_eq!(observer.evaluate(C64::ZERO).unwrap().pressure_pa, C64::ZERO);
    assert_eq!(
        observer.evaluate(C64::ONE).unwrap().pressure_pa,
        prepare(layer, 700.0, source, receiver)
            .evaluate(C64::ONE)
            .unwrap()
            .pressure_pa
    );
}
