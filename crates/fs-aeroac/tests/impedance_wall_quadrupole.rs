//! Original volume-stress trace, shear and near-field motion at a passive wall.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::CoherentCompactField;
use fs_aeroac::curle3d::quadrupole::QuadrupolePhasor;
use fs_aeroac::curle3d::refraction::full_wave::FullWaveForceLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::wall::{
    PlanarImpedanceWall, PreparedImpedanceWallForce, PreparedImpedanceWallQuadrupole,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

const HZ: f64 = 700.0;
const RHO: f64 = 1.2;
const SPEED: f64 = 343.0;
const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn limits() -> FullWaveQuadrupoleLimits {
    FullWaveQuadrupoleLimits {
        absolute_pressure_tolerance_pa_per_n_m: 1e-8,
        absolute_velocity_tolerance_m_s_per_n_m: 1e-11,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn wall(normalized: Option<C64>) -> PlanarImpedanceWall {
    PlanarImpedanceWall::new(
        FluidAcousticState::try_new(RHO, SPEED).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        normalized.map(|z| z.scale(RHO * SPEED)),
    )
    .unwrap()
}
fn basis(entry: usize) -> QuadrupolePhasor {
    let mut source = QuadrupolePhasor::default();
    source.stress_moment_n_m[entry] = C64::ONE;
    source
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
fn components(value: &CoherentCompactField) -> [C64; 4] {
    [
        value.pressure_pa,
        value.particle_velocity_m_s[0],
        value.particle_velocity_m_s[1],
        value.particle_velocity_m_s[2],
    ]
}
fn close(actual: C64, expected: C64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn matrix(source: QuadrupolePhasor) -> [[C64; 3]; 3] {
    let mut matrix = [[C64::ZERO; 3]; 3];
    for ((i, j), value) in ENTRIES.into_iter().zip(source.stress_moment_n_m) {
        matrix[i][j] = value;
        matrix[j][i] = value;
    }
    matrix
}

// Independent radial second/third derivatives of exp(i*k*r)/(4*pi*r).
// p = (G'/r) tr(M) + (G''-G'/r) n.M.n; u = grad(p)/(i*rho*omega).
fn quadrupole_green(source: [f64; 3], receiver: [f64; 3], stress: QuadrupolePhasor) -> [C64; 4] {
    let displacement: [f64; 3] = core::array::from_fn(|i| receiver[i] - source[i]);
    let radius = displacement.iter().map(|x| x * x).sum::<f64>().sqrt();
    let n = displacement.map(|x| x / radius);
    let moment = matrix(stress);
    let mn: [C64; 3] =
        core::array::from_fn(|i| (0..3).fold(C64::ZERO, |sum, j| sum + moment[i][j].scale(n[j])));
    let radial = (0..3).fold(C64::ZERO, |sum, i| sum + mn[i].scale(n[i]));
    let trace = moment[0][0] + moment[1][1] + moment[2][2];
    let omega = TAU * HZ;
    let k = omega / SPEED;
    let (sin, cos) = (k * radius).sin_cos();
    let green = C64::new(cos, sin).scale(1.0 / (4.0 * PI * radius));
    let a = C64::new(-1.0 / radius.powi(2), k / radius);
    let b = C64::new(3.0 / radius.powi(2) - k * k, -3.0 * k / radius);
    let third = C64::new(
        6.0 * k * k / radius - 15.0 / radius.powi(3),
        -k.powi(3) + 15.0 * k / radius.powi(2),
    );
    let velocity: [C64; 3] = core::array::from_fn(|i| {
        green
            * (b.scale(1.0 / radius) * (mn[i].scale(2.0) + trace.scale(n[i]))
                + third * radial.scale(n[i]))
            * C64::new(0.0, -1.0 / (RHO * omega))
    });
    [
        green * (a * trace + b * radial),
        velocity[0],
        velocity[1],
        velocity[2],
    ]
}

#[test]
fn ideal_images_transform_both_stress_indices_and_preserve_coherent_intensity() {
    let origin = [0.2, -0.1, 0.4];
    let normal = [1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0];
    let tangent = [2.0 / 5.0_f64.sqrt(), -1.0 / 5.0_f64.sqrt(), 0.0];
    let source = core::array::from_fn(|i| origin[i] - 0.08 * normal[i] + 0.03 * tangent[i]);
    let receiver = core::array::from_fn(|i| origin[i] - 0.13 * normal[i] + 0.2 * tangent[i]);
    let image = core::array::from_fn(|i| source[i] + 0.16 * normal[i]);
    let moment = matrix(stress());
    let mirror =
        |i: usize, j: usize| (if i == j { 1.0 } else { 0.0 }) - 2.0 * normal[i] * normal[j];
    let mirrored = QuadrupolePhasor {
        stress_moment_n_m: ENTRIES.map(|(i, j)| {
            let mut value = C64::ZERO;
            for (a, row) in moment.iter().enumerate() {
                for (b, entry) in row.iter().enumerate() {
                    value = value + entry.scale(mirror(i, a) * mirror(j, b));
                }
            }
            value
        }),
    };
    let direct = quadrupole_green(source, receiver, stress());
    for (impedance, sign) in [(None, 1.0), (Some(C64::ZERO), -1.0)] {
        let plane = PlanarImpedanceWall::new(
            FluidAcousticState::try_new(RHO, SPEED).unwrap(),
            origin,
            normal,
            impedance,
        )
        .unwrap();
        let owner = PreparedImpedanceWallQuadrupole::new(
            plane,
            HZ,
            source,
            receiver,
            FullWaveQuadrupoleLimits {
                maximum_evaluations: 1,
                ..limits()
            },
            &mut || false,
        )
        .unwrap();
        let reflected = quadrupole_green(image, receiver, mirrored).map(|z| z.scale(sign));
        let total = owner.evaluate_total(stress()).unwrap();
        for (actual, expected) in [
            (
                components(&owner.evaluate_reflected(stress()).unwrap()),
                reflected,
            ),
            (
                components(&total),
                core::array::from_fn(|i| direct[i] + reflected[i]),
            ),
        ] {
            for i in 0..4 {
                close(actual[i], expected[i], 5e-12, 1e-13);
            }
        }
        for i in 0..3 {
            close(
                C64::new(
                    total.active_intensity_w_m2[i],
                    total.reactive_intensity_w_m2[i],
                ),
                (direct[0] + reflected[0]).scale(0.5) * (direct[i + 1] + reflected[i + 1]).conj(),
                1e-11,
                1e-14,
            );
        }
        assert!(owner.diagnostics().exact_ideal_boundary);
        assert_eq!(owner.diagnostics().evaluations, 0);
    }
}

#[test]
fn finite_wall_stress_matches_all_original_source_derivatives_of_force() {
    let source = [0.025, -0.01, -0.08];
    let receiver = [0.17, 0.04, -0.11];
    let delta = 2e-4;
    for z in [C64::new(1.3, -0.4), C64::new(0.01, 1.0)] {
        let plane = wall(Some(z));
        let owner = PreparedImpedanceWallQuadrupole::new(
            plane,
            HZ,
            source,
            receiver,
            limits(),
            &mut || false,
        )
        .unwrap();
        let mut derivatives = [[[C64::ZERO; 4]; 3]; 3];
        for (j, derivative) in derivatives.iter_mut().enumerate() {
            for (offset, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut shifted = source;
                shifted[j] += offset * delta;
                let force = PreparedImpedanceWallForce::new(
                    plane,
                    HZ,
                    shifted,
                    receiver,
                    FullWaveForceLimits {
                        absolute_pressure_tolerance_pa_per_n: 1e-10,
                        absolute_velocity_tolerance_m_s_per_n: 1e-13,
                        relative_tolerance: 1e-11,
                        maximum_evaluations: 262_144,
                        maximum_depth: 40,
                    },
                    &mut || false,
                )
                .unwrap();
                for (i, derivative_i) in derivative.iter_mut().enumerate() {
                    let unit = core::array::from_fn(|a| if a == i { C64::ONE } else { C64::ZERO });
                    let value = components(&force.evaluate_reflected(unit).unwrap());
                    for (entry, value) in derivative_i.iter_mut().zip(value) {
                        *entry = *entry + value.scale(weight / (12.0 * delta));
                    }
                }
            }
        }
        for (entry, (i, j)) in ENTRIES.into_iter().enumerate() {
            let actual = components(&owner.evaluate_reflected(basis(entry)).unwrap());
            for (component, actual) in actual.into_iter().enumerate() {
                let expected = derivatives[j][i][component]
                    + if i == j {
                        C64::ZERO
                    } else {
                        derivatives[i][j][component]
                    };
                close(
                    actual,
                    expected,
                    2e-7,
                    if component == 0 { 2e-6 } else { 2e-9 },
                );
            }
        }
    }
}

#[test]
fn rotated_boundary_obeys_robin_law_and_absorbs_positive_energy_for_every_stress() {
    let origin = [0.2, -0.1, 0.4];
    let normal = [1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0];
    let tangent = [2.0 / 5.0_f64.sqrt(), -1.0 / 5.0_f64.sqrt(), 0.0];
    let source = core::array::from_fn(|i| origin[i] - 0.07 * normal[i] + 0.04 * tangent[i]);
    for normalized in [C64::new(1.7, -0.4), C64::new(0.03, 1.1)] {
        let impedance = normalized.scale(RHO * SPEED);
        let plane = PlanarImpedanceWall::new(
            FluidAcousticState::try_new(RHO, SPEED).unwrap(),
            origin,
            normal,
            Some(impedance),
        )
        .unwrap();
        let owner =
            PreparedImpedanceWallQuadrupole::new(plane, HZ, source, origin, limits(), &mut || {
                false
            })
            .unwrap();
        for source in (0..6).map(basis).chain([stress()]) {
            let value = owner.evaluate_total(source).unwrap();
            let un = (0..3).fold(C64::ZERO, |sum, i| {
                sum + value.particle_velocity_m_s[i].scale(normal[i])
            });
            close(value.pressure_pa, impedance * un, 3e-8, 2e-7);
            let flux = (0..3)
                .map(|i| value.active_intensity_w_m2[i] * normal[i])
                .sum::<f64>();
            let absorbed = 0.5 * impedance.re * un.abs().powi(2);
            assert!(flux >= -1e-8);
            assert!((flux - absorbed).abs() <= 1e-7 * absorbed.abs().max(1.0));
        }
    }
}

#[test]
fn narrow_surface_mode_retains_trace_and_shear_velocity_at_axial_pressure_nulls() {
    // Independent axial closed form G=F/(4*pi), H=-(zs+zo), s=i*k/z:
    // F=exp(i*k*H)/H + 2*s*exp(-s*H)*E1(-H*(s+i*k)). Principal-branch E1
    // uses its convergent complex power series, then analytic derivatives via
    // F'+s*F=g'-s*g, g=exp(i*k*H)/H. Axis Helmholtz supplies transverse
    // Hessians: G_xx=-(k*k*F+F'')/(8*pi), G_zz=F''/(4*pi).
    let budget = FullWaveQuadrupoleLimits {
        maximum_evaluations: 16_384,
        ..limits()
    };
    let owner = PreparedImpedanceWallQuadrupole::new(
        wall(Some(C64::new(0.001, 1.0))),
        HZ,
        [0.0, 0.0, -0.04],
        [0.0, 0.0, -0.02],
        budget,
        &mut || false,
    )
    .unwrap();
    let xx = owner.evaluate_reflected(basis(0)).unwrap();
    let yy = owner.evaluate_reflected(basis(1)).unwrap();
    let zz = owner.evaluate_reflected(basis(2)).unwrap();
    let xz = owner.evaluate_reflected(basis(4)).unwrap();
    let yz = owner.evaluate_reflected(basis(5)).unwrap();
    for (actual, expected) in [
        (
            xx.pressure_pa,
            C64::new(-761.1725508497339, -438.26520341309293),
        ),
        (
            zz.pressure_pa,
            C64::new(1577.79919345388, 531.9412133706572),
        ),
        (
            xx.particle_velocity_m_s[2],
            C64::new(-1.2157946130135568, 6.850711519552071),
        ),
        (
            zz.particle_velocity_m_s[2],
            C64::new(1.1270805771689134, -12.593481620516423),
        ),
        (
            xz.particle_velocity_m_s[0],
            C64::new(2.4315892260271137, -13.701423039104142),
        ),
    ] {
        close(actual, expected, 2e-8, 1e-8);
    }
    assert_eq!(components(&xx), components(&yy));
    assert_eq!(xz.pressure_pa, C64::ZERO);
    assert_eq!(yz.pressure_pa, C64::ZERO);
    assert_eq!(xz.particle_velocity_m_s[0], yz.particle_velocity_m_s[1]);
    assert_eq!(xz.particle_velocity_m_s[1..], [C64::ZERO; 2]);
    let isotropic = QuadrupolePhasor {
        stress_moment_n_m: [
            C64::ONE,
            C64::ONE,
            C64::ONE,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ],
    };
    close(
        owner.evaluate_reflected(isotropic).unwrap().pressure_pa,
        xx.pressure_pa.scale(2.0) + zz.pressure_pa,
        2e-13,
        1e-12,
    );
    assert!(owner.diagnostics().evaluations <= budget.maximum_evaluations);
    assert!(!owner.diagnostics().exact_ideal_boundary);
}

#[test]
fn reflected_self_field_is_finite_and_singular_total_or_invalid_sources_refuse() {
    let source = [0.02, -0.01, -0.08];
    let plane = wall(Some(C64::new(1.3, -0.4)));
    let owner =
        PreparedImpedanceWallQuadrupole::new(plane, HZ, source, source, limits(), &mut || false)
            .unwrap();
    assert!(!owner.has_total_field());
    assert!(owner.evaluate_reflected(stress()).is_ok());
    assert!(owner.evaluate_total(stress()).is_err());
    assert!(owner.evaluate_total_harmonic_line(stress()).is_err());
    assert_eq!(
        components(
            &owner
                .evaluate_reflected(QuadrupolePhasor::default())
                .unwrap()
        ),
        [C64::ZERO; 4]
    );
    assert!(
        owner
            .evaluate_reflected(QuadrupolePhasor {
                stress_moment_n_m: [C64::new(f64::NAN, 0.0); 6]
            })
            .is_err()
    );
    let line = owner.evaluate_reflected_harmonic_line(stress()).unwrap();
    assert_eq!(line.frequency_hz, HZ);
    assert_eq!(
        line.pressure_pa,
        owner.evaluate_reflected(stress()).unwrap().pressure_pa
    );
    assert_eq!(owner.source_position_m(), source);
    assert_eq!(owner.observer_position_m(), source);
    assert!(matches!(
        PreparedImpedanceWallQuadrupole::new(plane, HZ, source, source, limits(), &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PreparedImpedanceWallQuadrupole::new(
            plane,
            HZ,
            source,
            source,
            FullWaveQuadrupoleLimits {
                maximum_evaluations: 30,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
}
