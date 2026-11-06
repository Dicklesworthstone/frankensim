//! G1/G3/G4: original volume stress on the incident side of a fluid interface.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::CoherentCompactField;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::FullWaveForceLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::{
    FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole,
};
use fs_aeroac::curle3d::refraction::full_wave::reflection::{
    PreparedReflectedFullWaveForce, PreparedReflectedFullWaveQuadrupole,
};
use fs_aeroac::observation::{
    HarmonicFieldLine, HarmonicObservationBudget, HarmonicSampleClock, PreparedHarmonicField,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};
use std::f64::consts::TAU;

const AXES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn fluid(rho: f64, speed: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, speed).unwrap()
}
fn interface(speed: [f64; 2]) -> PlanarFluidInterface {
    PlanarFluidInterface::new(
        fluid(1.2, speed[0]),
        fluid(2.4, speed[1]),
        [0.0; 3],
        [0.0, 0.0, 1.0],
    )
    .unwrap()
}
fn limits() -> FullWaveQuadrupoleLimits {
    FullWaveQuadrupoleLimits {
        absolute_pressure_tolerance_pa_per_n_m: 2e-6,
        absolute_velocity_tolerance_m_s_per_n_m: 2e-9,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn force_limits() -> FullWaveForceLimits {
    FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-9,
        absolute_velocity_tolerance_m_s_per_n: 1e-12,
        relative_tolerance: 1e-11,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn prepare(
    plane: PlanarFluidInterface,
    source: [f64; 3],
    observer: [f64; 3],
) -> PreparedReflectedFullWaveQuadrupole {
    PreparedReflectedFullWaveQuadrupole::new(plane, 700.0, source, observer, limits(), &mut || {
        false
    })
    .unwrap()
}
fn stress() -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(0.0007, 0.0001),
            C64::new(-0.0003, 0.0002),
            C64::new(0.0009, -0.0004),
            C64::new(0.00023, -0.00015),
            C64::new(-0.00011, 0.0003),
            C64::new(0.0004, 0.0002),
        ],
    }
}
fn basis(entry: usize) -> QuadrupolePhasor {
    let mut source = QuadrupolePhasor::default();
    source.stress_moment_n_m[entry] = C64::ONE;
    source
}
fn amplitudes(line: HarmonicFieldLine<3>) -> [C64; 4] {
    let [x, y, z] = line.particle_velocity_m_s;
    [line.pressure_pa, x, y, z]
}
fn close(actual: C64, expected: C64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:?}, expected={expected:?}, rel={relative:e}, abs={absolute:e}"
    );
}
fn compare_field(actual: &CoherentCompactField, expected: &CoherentCompactField) {
    close(actual.pressure_pa, expected.pressure_pa, 3e-12, 1e-14);
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
            3e-12,
            1e-16,
        );
        close(
            C64::new(
                actual.active_intensity_w_m2[axis],
                actual.reactive_intensity_w_m2[axis],
            ),
            expected.pressure_pa.scale(0.5) * expected.particle_velocity_m_s[axis].conj(),
            6e-12,
            1e-16,
        );
    }
}

#[test]
fn original_equal_speed_image_tensor_preserves_trace_shear_density_and_matched_zero() {
    let normal = [1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0];
    let tangent = [2.0 / 5.0_f64.sqrt(), -1.0 / 5.0_f64.sqrt(), 0.0];
    let origin = [0.2, -0.1, 0.4];
    let source = core::array::from_fn(|i| origin[i] - 0.08 * normal[i] + 0.03 * tangent[i]);
    let observer = core::array::from_fn(|i| origin[i] - 0.13 * normal[i] + 0.2 * tangent[i]);
    let image = core::array::from_fn(|i| source[i] + 0.16 * normal[i]);
    let original =
        PreparedQuadrupoleObserver::new(fluid(1.2, 343.0), 700.0, observer, &[source, image], 2)
            .unwrap();
    for rho2 in [0.6, 1.2, 2.4] {
        let plane =
            PlanarFluidInterface::new(fluid(1.2, 343.0), fluid(rho2, 343.0), origin, normal)
                .unwrap();
        let prepared = PreparedReflectedFullWaveQuadrupole::new(
            plane,
            700.0,
            source,
            observer,
            FullWaveQuadrupoleLimits {
                maximum_evaluations: 1,
                ..limits()
            },
            &mut || false,
        )
        .unwrap();
        assert_eq!(prepared.diagnostics().evaluations, 0);
        assert!(prepared.diagnostics().exact_equal_speed);
        assert_eq!(
            prepared.diagnostics().pressure_tail_bound_pa_per_n_m,
            [0.0; 4]
        );
        assert_eq!(
            prepared.diagnostics().velocity_tail_bound_m_s_per_n_m,
            [0.0; 10]
        );
        for moment in [stress(), basis(0), basis(3), basis(4), basis(5)] {
            // Independent full-matrix S*M*S image in the original Cartesian
            // frame: one normal source index reverses; two do not.
            let mut tensor = [[C64::ZERO; 3]; 3];
            for ((i, j), entry) in AXES.into_iter().zip(moment.stress_moment_n_m) {
                tensor[i][j] = entry;
                tensor[j][i] = entry;
            }
            let mirror =
                |i: usize, j: usize| (if i == j { 1.0 } else { 0.0 }) - 2.0 * normal[i] * normal[j];
            let reflected_moment = QuadrupolePhasor {
                stress_moment_n_m: AXES.map(|(i, j)| {
                    let mut value = C64::ZERO;
                    for a in 0..3 {
                        for b in 0..3 {
                            value = value + tensor[a][b].scale(mirror(i, a) * mirror(j, b));
                        }
                    }
                    value.scale((rho2 - 1.2) / (rho2 + 1.2))
                }),
            };
            let reflected = prepared.evaluate_reflected(moment).unwrap();
            compare_field(
                &reflected,
                &original
                    .evaluate(&[QuadrupolePhasor::default(), reflected_moment])
                    .unwrap(),
            );
            compare_field(
                &prepared.evaluate_total(moment).unwrap(),
                &original.evaluate(&[moment, reflected_moment]).unwrap(),
            );
            if rho2 == 1.2 {
                assert_eq!(reflected.pressure_pa, C64::ZERO);
                assert_eq!(reflected.particle_velocity_m_s, [C64::ZERO; 3]);
            }
        }
    }
}

#[test]
fn both_speed_orders_match_independent_spectra_original_force_derivatives_and_momentum() {
    // Direct-q mpmath tanh-sinh, original R, independent J2/J3 angular moments,
    // 35/50-digit refinement agreeing at all 30 printed digits. Tail t>8000:
    // p<2.64e-51 Pa/(N m), u<4.03e-51 (m/s)/(N m), per unweighted channel.
    let references = [
        [
            C64::new(498.93601771026204, 2.046949742356563),
            C64::new(-348.49423773064916, 3.172118820117301),
            C64::new(323.24374617173567, -11.494976059035309),
            C64::new(-285.08014350971081, -56.065567305010823),
            C64::new(-0.015579603256770408, 4.8169180286618278),
            C64::new(-0.0053296494208265053, -4.0140688714002541),
            C64::new(-0.050788154231183887, 5.4327289676557573),
            C64::new(0.015741953694547946, -1.370537198706824),
            C64::new(-0.0053296494208265053, -4.0140688714002541),
            C64::new(-0.054448876801144762, -1.5311261912730254),
            C64::new(0.050788154231183887, -5.4327289676557573),
            C64::new(0.054448876801144762, 1.5311261912730254),
            C64::new(-0.015741953694547946, 1.370537198706824),
            C64::new(0.02885849845298101, 4.1241015915888131),
        ],
        [
            C64::new(489.44027131726472, -25.015892764939988),
            C64::new(-377.51667671044709, -26.379071149172251),
            C64::new(404.16999153671321, 9.0978670090757827),
            C64::new(-135.27902702546939, 29.773631691241036),
            C64::new(0.019134511653560816, 4.6086984122805719),
            C64::new(0.0064570410168627215, -4.1065620156688611),
            C64::new(0.040788431389963152, 5.4741848947589258),
            C64::new(-0.014848847616435871, -0.65587760637334402),
            C64::new(0.0064570410168627215, -4.1065620156688611),
            C64::new(0.043094360300210958, -1.914453928645161),
            C64::new(-0.040788431389963152, -5.4741848947589258),
            C64::new(-0.043094360300210958, 1.914453928645161),
            C64::new(0.014848847616435871, 0.65587760637334402),
            C64::new(0.045997996229352356, 3.6338413792021642),
        ],
    ];
    let source = [0.0, 0.0, -0.006];
    let observer = [0.04, 0.0, -0.011];
    let step = 2e-5;
    for (speeds, expected) in [[300.0, 600.0], [600.0, 300.0]].into_iter().zip(references) {
        let plane = interface(speeds);
        let prepared = prepare(plane, source, observer);
        let values: [[C64; 4]; 6] = core::array::from_fn(|entry| {
            amplitudes(
                prepared
                    .evaluate_reflected_harmonic_line(basis(entry))
                    .unwrap(),
            )
        });
        let channels = [
            values[0][0],
            values[1][0],
            values[4][0].scale(0.5),
            values[2][0],
            values[0][1],
            values[1][1],
            values[4][1].scale(0.5),
            values[2][1],
            values[3][2].scale(0.5),
            values[5][2].scale(0.5),
            values[0][3],
            values[1][3],
            values[4][3].scale(0.5),
            values[2][3],
        ];
        for (index, (actual, expected)) in channels.into_iter().zip(expected).enumerate() {
            close(actual, expected, 3e-8, if index < 4 { 2e-6 } else { 2e-9 });
        }
        let diagnostics = prepared.diagnostics();
        assert!(!diagnostics.exact_equal_speed);
        assert!(
            diagnostics.evaluations > 0 && diagnostics.evaluations <= limits().maximum_evaluations
        );
        let norm = |v: &[f64]| v.iter().fold(0.0_f64, |sum, value| sum.hypot(*value));
        assert!(
            norm(&diagnostics.pressure_tail_bound_pa_per_n_m)
                <= 0.251 * limits().absolute_pressure_tolerance_pa_per_n_m
        );
        assert!(
            norm(&diagnostics.velocity_tail_bound_m_s_per_n_m)
                <= 0.251 * limits().absolute_velocity_tolerance_m_s_per_n_m
        );

        // Differentiate the independently existing reflected FORCE owner with
        // respect to each ORIGINAL source coordinate. Finite differences are
        // confined to this test; the product evaluates analytic derivatives.
        let mut gradient = [[[C64::ZERO; 3]; 4]; 3];
        for direction in 0..3 {
            for (offset, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut position = source;
                position[direction] += offset * step;
                let force = PreparedReflectedFullWaveForce::new(
                    plane,
                    700.0,
                    position,
                    observer,
                    force_limits(),
                    &mut || false,
                )
                .unwrap();
                for column in 0..3 {
                    let amplitude =
                        core::array::from_fn(
                            |axis| if axis == column { C64::ONE } else { C64::ZERO },
                        );
                    let field =
                        amplitudes(force.evaluate_reflected_harmonic_line(amplitude).unwrap());
                    for row in 0..4 {
                        gradient[direction][row][column] = gradient[direction][row][column]
                            + field[row].scale(weight / (12.0 * step));
                    }
                }
            }
        }
        for (entry, (i, j)) in AXES.into_iter().enumerate() {
            for row in 0..4 {
                let expected = gradient[j][row][i]
                    + if i == j {
                        C64::ZERO
                    } else {
                        gradient[i][row][j]
                    };
                close(
                    values[entry][row],
                    expected,
                    3e-6,
                    if row == 0 { 2e-5 } else { 2e-8 },
                );
            }
        }
        let actual = prepared.evaluate_reflected(stress()).unwrap();
        for axis in 0..3 {
            let mut derivative = C64::ZERO;
            for (offset, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut position = observer;
                position[axis] += offset * step;
                derivative = derivative
                    + prepare(plane, source, position)
                        .evaluate_reflected(stress())
                        .unwrap()
                        .pressure_pa
                        .scale(weight / (12.0 * step));
            }
            let velocity = derivative * C64::new(0.0, -1.0 / (TAU * 700.0 * 1.2));
            close(actual.particle_velocity_m_s[axis], velocity, 3e-6, 2e-10);
        }
    }
}

#[test]
fn original_interface_continuity_coherent_energy_and_source_clock_are_retained() {
    let source = [0.004, -0.006, -0.018];
    let observer = [0.043, 0.012, -0.027];
    for speeds in [[300.0, 600.0], [600.0, 300.0]] {
        let plane = interface(speeds);
        // One-sided linear extrapolation to the exact physical plane. Native
        // fields are prepared strictly on their own sides; both p and u_n
        // must agree for the same original stress and source clock.
        let mut pressure = [C64::ZERO; 2];
        let mut velocity = [C64::ZERO; 2];
        for (multiple, coefficient) in [(1.0, 2.0), (2.0, -1.0)] {
            let depth = multiple * 1e-6;
            let reflected = prepare(plane, source, [0.043, 0.012, -depth])
                .evaluate_total(stress())
                .unwrap();
            let transmitted = PreparedFullWaveQuadrupole::new(
                plane,
                700.0,
                source,
                [0.043, 0.012, depth],
                limits(),
                &mut || false,
            )
            .unwrap()
            .evaluate(stress())
            .unwrap();
            for (index, field) in [reflected, transmitted].into_iter().enumerate() {
                pressure[index] = pressure[index] + field.pressure_pa.scale(coefficient);
                velocity[index] =
                    velocity[index] + field.particle_velocity_m_s[2].scale(coefficient);
            }
        }
        close(pressure[0], pressure[1], 3e-6, 2e-9);
        close(velocity[0], velocity[1], 3e-6, 2e-11);

        let prepared = prepare(plane, source, observer);
        let direct =
            PreparedQuadrupoleObserver::new(plane.incident_state(), 700.0, observer, &[source], 1)
                .unwrap()
                .evaluate(&[stress()])
                .unwrap();
        let reflected = prepared.evaluate_reflected(stress()).unwrap();
        let total = prepared.evaluate_total(stress()).unwrap();
        close(
            total.pressure_pa,
            direct.pressure_pa + reflected.pressure_pa,
            3e-13,
            1e-14,
        );
        let mut largest_cross = 0.0_f64;
        let mut largest_flux = 0.0_f64;
        for axis in 0..3 {
            close(
                total.particle_velocity_m_s[axis],
                direct.particle_velocity_m_s[axis] + reflected.particle_velocity_m_s[axis],
                3e-13,
                1e-16,
            );
            let cross = (direct.pressure_pa * reflected.particle_velocity_m_s[axis].conj()
                + reflected.pressure_pa * direct.particle_velocity_m_s[axis].conj())
            .scale(0.5);
            let separate = C64::new(
                direct.active_intensity_w_m2[axis] + reflected.active_intensity_w_m2[axis],
                direct.reactive_intensity_w_m2[axis] + reflected.reactive_intensity_w_m2[axis],
            );
            let flux = C64::new(
                total.active_intensity_w_m2[axis],
                total.reactive_intensity_w_m2[axis],
            );
            close(flux, separate + cross, 5e-13, 1e-16);
            largest_cross = largest_cross.max(cross.abs());
            largest_flux = largest_flux.max(flux.abs());
        }
        assert!(largest_cross > 0.01 * largest_flux);
        for (field, line) in [
            (
                reflected,
                prepared.evaluate_reflected_harmonic_line(stress()).unwrap(),
            ),
            (
                total,
                prepared.evaluate_total_harmonic_line(stress()).unwrap(),
            ),
        ] {
            assert_eq!(line.frequency_hz, 700.0);
            assert_eq!(line.pressure_pa, field.pressure_pa);
            assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
            let observation = PreparedHarmonicField::new(
                &[line],
                HarmonicObservationBudget {
                    max_lines: 1,
                    max_line_pairs: 1,
                },
            )
            .unwrap();
            // 700*2^30 complete cycles, then a quarter period: propagated
            // phasors must select their imaginary parts with no added delay.
            let clock = HarmonicSampleClock::new(Dd::from_f64(1_073_741_824.0), 2800).unwrap();
            let sample = observation
                .field_at_clock_dd(clock.time_at(1).unwrap())
                .unwrap();
            close(
                C64::from_re(sample.pressure_pa),
                C64::from_re(field.pressure_pa.im),
                5e-12,
                1e-13,
            );
            for axis in 0..3 {
                close(
                    C64::from_re(sample.particle_velocity_m_s[axis]),
                    C64::from_re(field.particle_velocity_m_s[axis].im),
                    5e-12,
                    1e-15,
                );
            }
        }
    }
}

#[test]
fn axis_self_field_domains_and_bounded_cancellation_remain_explicit() {
    let plane = interface([300.0, 600.0]);
    let source = [0.0, 0.0, -0.04];
    let observer = [0.0, 0.0, -0.07];
    let prepared = prepare(plane, source, observer);
    let xz = prepared.evaluate_reflected(basis(4)).unwrap();
    let yz = prepared.evaluate_reflected(basis(5)).unwrap();
    assert_eq!(xz.pressure_pa, C64::ZERO);
    assert_eq!(yz.pressure_pa, C64::ZERO);
    assert!(xz.particle_velocity_m_s[0].abs() > 1e-6);
    close(
        xz.particle_velocity_m_s[0],
        yz.particle_velocity_m_s[1],
        2e-13,
        0.0,
    );
    assert!(prepared.has_total_field());
    let at_source = prepare(plane, source, source);
    assert!(!at_source.has_total_field());
    assert!(
        at_source
            .evaluate_reflected(stress())
            .unwrap()
            .pressure_pa
            .abs()
            .is_finite()
    );
    assert!(at_source.evaluate_total(stress()).is_err());
    assert!(at_source.evaluate_total_harmonic_line(stress()).is_err());
    let silent = prepared
        .evaluate_total(QuadrupolePhasor::default())
        .unwrap();
    assert_eq!(silent.pressure_pa, C64::ZERO);
    assert_eq!(silent.particle_velocity_m_s, [C64::ZERO; 3]);
    let invalid = QuadrupolePhasor {
        stress_moment_n_m: [C64::new(f64::NAN, 0.0); 6],
    };
    assert!(prepared.evaluate_reflected(invalid).is_err());
    assert!(prepared.evaluate_total(invalid).is_err());
    let new =
        |s, o, l| PreparedReflectedFullWaveQuadrupole::new(plane, 700.0, s, o, l, &mut || false);
    assert!(
        new(
            source,
            observer,
            FullWaveQuadrupoleLimits {
                maximum_evaluations: 1,
                ..limits()
            }
        )
        .is_err()
    );
    for control in [
        FullWaveQuadrupoleLimits {
            relative_tolerance: f64::NAN,
            ..limits()
        },
        FullWaveQuadrupoleLimits {
            absolute_velocity_tolerance_m_s_per_n_m: 0.0,
            ..limits()
        },
        FullWaveQuadrupoleLimits {
            maximum_depth: 0,
            ..limits()
        },
    ] {
        assert!(new(source, observer, control).is_err());
    }
    for position in [[0.0; 3], [0.0, 0.0, 0.04], [f64::INFINITY, 0.0, -0.04]] {
        assert!(new(position, observer, limits()).is_err());
        assert!(new(source, position, limits()).is_err());
    }
    for hz in [0.0, f64::NAN, f64::INFINITY] {
        assert!(
            PreparedReflectedFullWaveQuadrupole::new(
                plane,
                hz,
                source,
                observer,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(matches!(
        PreparedReflectedFullWaveQuadrupole::new(
            plane,
            700.0,
            source,
            observer,
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedReflectedFullWaveQuadrupole::new(
            plane,
            700.0,
            source,
            observer,
            limits(),
            &mut || {
                polls += 1;
                polls > 6
            },
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}
