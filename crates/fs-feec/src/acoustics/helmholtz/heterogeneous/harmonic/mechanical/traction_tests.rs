//! Physical force, source work and exact frequency tangents on original FE walls.

use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_phs::{WallPin, acoustic_outlet::AcousticWallRelaxation};
use fs_rep_mesh::TetComplex;

const POINTS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const MEDIUM: Medium = Medium {
    density_kg_m3: 2.0,
    sound_speed_m_s: 3.0,
};
const WALL: [u32; 3] = [0, 1, 3];
const ROBIN: [u32; 3] = [0, 2, 3];
const RELEASE: [u32; 3] = [0, 1, 2];
const SHAPE: [f64; 3] = [0.8, -0.4, 1.3];

fn uniform<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn wall() -> MechanicalBoundary {
    MechanicalBoundary {
        faces: vec![WALL],
        wall: WallPin {
            surface_density: 0.6,
            stiffness_per_area: 3.0,
            resistance: 0.4,
        },
        relaxation: vec![AcousticWallRelaxation {
            stiffness_per_area: 2.5,
            relaxation_time_s: 0.06,
        }],
    }
}

fn prepare<'m, 'f, 'c>(
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    omega: f64,
    wall: &MechanicalBoundary,
    release: &[[u32; 3]],
) -> PreparedMechanicalHarmonic<'m, 'f, 'c> {
    model
        .prepare_with_coupled_mechanical_walls(
            omega / core::f64::consts::TAU,
            &[ResistiveBoundary {
                faces: vec![ROBIN],
                admittance_m_pa_s: 0.12,
            }],
            core::slice::from_ref(wall),
            release,
            MechanicalHarmonicLimits::default(),
            &mut || false,
        )
        .unwrap()
}

fn profile() -> WallTractionProfile {
    WallTractionProfile {
        wall: 0,
        coefficients: vec![(0, SHAPE[0]), (1, SHAPE[1]), (3, SHAPE[2])],
    }
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 3e-10 * expected.abs().max(1.0),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn real_close(actual: f64, expected: f64) {
    close(C64::from_re(actual), C64::from_re(expected));
}

fn form(left: &[C64], right: &[C64]) -> f64 {
    let mut result = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            result += (left[i].conj() * right[j]).re / if i == j { 12.0 } else { 24.0 };
        }
    }
    result
}

fn dynamic(omega: f64) -> (C64, C64) {
    let rate = 1.0 / 0.06;
    let denominator = C64::new(rate, omega);
    (
        C64::new(3.0 - 0.6 * omega * omega, 0.4 * omega) + C64::new(0.0, 2.5 * omega) / denominator,
        C64::new(-1.2 * omega, 0.4) + C64::new(0.0, 2.5 * rate) / denominator / denominator,
    )
}

#[test]
fn independent_full_equations_and_analytic_spectra_preserve_signed_force_work() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POINTS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let amplitude = C64::new(0.7, -0.2);
    let amplitude_slope = C64::new(-0.03, 0.04);
    let speed = C64::new(0.06, 0.03);
    let speed_slope = C64::new(0.01, -0.02);
    // Both sides of the original-memory reconstruction's omega/rate branch.
    for omega in [1.7, 4.5, 30.0] {
        let prepared = prepare(&model, omega, &wall(), &[RELEASE]);
        let traction = prepared
            .prepare_wall_traction(&profile(), &mut || false)
            .unwrap();
        assert_eq!(prepared.wall_dofs(0).unwrap(), &[0, 1, 3]);
        let result = prepared
            .solve_with_wall_tractions_frequency_derivative(
                &[(WALL, [speed; 3])],
                Some(&[(WALL, [speed_slope; 3])]),
                &[],
                &[],
                &[(&traction, amplitude, amplitude_slope)],
                &mut || false,
            )
            .unwrap();
        // Independent unit-tetra moments: K33=1/12, M33=1/1080;
        // W_ii=1/12, W_ij=1/24, integral(phi3)=1/6. Only p3 is free.
        // D*w_j=p_j+shape_j*a, v=i*omega*w, z=c*rate*w/(rate+i*omega).
        let (d, dd) = dynamic(omega);
        let a = C64::new(1.0 / 12.0 - omega * omega / 1080.0, omega / 100.0);
        let da = C64::new(-omega / 540.0, 0.01);
        let b = speed.scale(1.0 / 6.0);
        let db = speed_slope.scale(1.0 / 6.0);
        let projected_shape = (SHAPE[0] + SHAPE[1]) / 24.0 + SHAPE[2] / 12.0;
        let denominator = a * d - C64::from_re(omega * omega / 12.0);
        let denominator_slope = da * d + a * dd - C64::from_re(omega / 6.0);
        let numerator =
            C64::new(0.0, omega) * b * d + amplitude.scale(omega * omega * projected_shape);
        let numerator_slope = C64::new(0.0, 1.0) * b * d
            + C64::new(0.0, omega) * (db * d + b * dd)
            + amplitude.scale(2.0 * omega * projected_shape)
            + amplitude_slope.scale(omega * omega * projected_shape);
        let pressure = numerator / denominator;
        let pressure_slope = (numerator_slope - pressure * denominator_slope) / denominator;
        close(result.response.acoustic.pressure_pa[3], pressure);
        close(result.pressure_derivative_pa_per_rad_s[3], pressure_slope);
        for global in [0, 1, 2] {
            assert_eq!(result.response.acoustic.pressure_pa[global], C64::ZERO);
            assert_eq!(result.pressure_derivative_pa_per_rad_s[global], C64::ZERO);
        }
        let phasors = &result.response.walls[0];
        let slopes = &result.wall_derivatives[0];
        let trace_pressure = [C64::ZERO, C64::ZERO, pressure];
        let rate = 1.0 / 0.06;
        let c = 2.5_f64.sqrt();
        for local in 0..3 {
            let p_slope = if local == 2 {
                pressure_slope
            } else {
                C64::ZERO
            };
            let w = (trace_pressure[local] + amplitude.scale(SHAPE[local])) / d;
            let dw = (p_slope + amplitude_slope.scale(SHAPE[local]) - w * dd) / d;
            let z_factor = C64::from_re(c * rate) / C64::new(rate, omega);
            let z = z_factor * w;
            close(phasors.displacement_m[local], w);
            close(phasors.normal_velocity_m_s[local], C64::new(0.0, omega) * w);
            close(phasors.relaxation_memory[0][local], z);
            close(slopes.displacement_m_per_rad_s[local], dw);
            close(
                slopes.normal_velocity_m_s_per_rad_s[local],
                C64::new(0.0, 1.0) * w + C64::new(0.0, omega) * dw,
            );
            close(
                slopes.relaxation_memory_per_rad_s[0][local],
                z_factor * dw - C64::new(0.0, 1.0) * z / C64::new(rate, omega),
            );
        }
        let force = SHAPE.map(|v| amplitude.scale(v));
        let extension: Vec<_> = phasors
            .displacement_m
            .iter()
            .zip(&phasors.relaxation_memory[0])
            .map(|(&w, &z)| w.scale(c) - z)
            .collect();
        let account = result.response.wall_accounting;
        real_close(
            account.wall_traction_power_w,
            0.5 * form(&force, &phasors.normal_velocity_m_s),
        );
        real_close(
            account.fluid_to_wall_power_w,
            0.5 * form(&trace_pressure, &phasors.normal_velocity_m_s),
        );
        real_close(
            account.kinetic_energy_j,
            0.15 * form(&phasors.normal_velocity_m_s, &phasors.normal_velocity_m_s),
        );
        real_close(
            account.elastic_energy_j,
            0.75 * form(&phasors.displacement_m, &phasors.displacement_m),
        );
        real_close(
            account.relaxation_energy_j,
            0.25 * form(&extension, &extension),
        );
        real_close(
            account.resistance_dissipation_w,
            0.2 * form(&phasors.normal_velocity_m_s, &phasors.normal_velocity_m_s),
        );
        real_close(
            account.relaxation_dissipation_w,
            0.5 * rate * form(&extension, &extension),
        );
        real_close(
            result.response.acoustic.input_power_w,
            0.5 * (pressure.conj() * b).re,
        );
        real_close(
            result.response.acoustic.boundary_dissipation_w,
            pressure.norm_sq() / 200.0,
        );
        real_close(result.response.power_balance_residual_w(), 0.0);
        real_close(result.response.fluid_power_balance_residual_w(), 0.0);
        real_close(account.power_balance_residual_w(), 0.0);
        assert!(result.response.acoustic.backward_error <= 5e-13);
        assert!(result.derivative_backward_error <= 5e-13);
    }
}

#[test]
fn p2_release_shared_tractions_retain_coherence_and_original_surface_work() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POINTS, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let omega = 3.4;
    let prepared = prepare(&model, omega, &wall(), &[RELEASE]);
    let dofs = prepared.wall_dofs(0).unwrap();
    assert_eq!(dofs.len(), 6);
    let profile = WallTractionProfile {
        wall: 0,
        coefficients: dofs
            .iter()
            .enumerate()
            .map(|(i, &global)| {
                (
                    global,
                    (0.2 + 0.3 * i as f64) * if i % 2 == 0 { 1.0 } else { -1.0 },
                )
            })
            .collect(),
    };
    let source = prepared
        .prepare_wall_traction(&profile, &mut || false)
        .unwrap();
    let amplitude = C64::new(0.7, -0.2);
    let slope = C64::new(-0.03, 0.04);
    let weak = prepared
        .solve_with_wall_tractions_frequency_derivative(
            &[],
            None,
            &[],
            &[],
            &[(&source, amplitude, slope)],
            &mut || false,
        )
        .unwrap();
    let huge = C64::new(1e16, -1e16);
    let coherent = prepared
        .solve_with_wall_tractions_frequency_derivative(
            &[],
            None,
            &[],
            &[],
            &[
                (&source, huge, huge),
                (&source, amplitude, slope),
                (&source, -huge, -huge),
            ],
            &mut || false,
        )
        .unwrap();
    for (a, b) in weak
        .response
        .acoustic
        .pressure_pa
        .iter()
        .zip(&coherent.response.acoustic.pressure_pa)
        .chain(
            weak.pressure_derivative_pa_per_rad_s
                .iter()
                .zip(&coherent.pressure_derivative_pa_per_rad_s),
        )
    {
        close(*a, *b);
    }
    real_close(
        coherent.response.wall_accounting.wall_traction_power_w,
        weak.response.wall_accounting.wall_traction_power_w,
    );
    let (d, _) = dynamic(omega);
    let mut velocity = vec![C64::ZERO; model.ndof()];
    let mut traction = vec![C64::ZERO; model.ndof()];
    for (local, &global) in dofs.iter().enumerate() {
        let value = amplitude.scale(source.coefficients()[local]);
        close(
            d * weak.response.walls[0].displacement_m[local]
                - weak.response.acoustic.pressure_pa[global],
            value,
        );
        velocity[global] = weak.response.walls[0].normal_velocity_m_s[local];
        traction[global] = value;
    }
    // Independent degree-four triangle integration of the physical work. No
    // mass matrix or cached source column participates in this observation.
    let mut work = 0.0;
    for (x, weight) in [
        (0.091576213509771, 0.109951743655322),
        (0.445948490915965, 0.223381589678011),
    ] {
        for triangle in [
            [1.0 - 2.0 * x, x, x],
            [x, 1.0 - 2.0 * x, x],
            [x, x, 1.0 - 2.0 * x],
        ] {
            let bary = [triangle[0], triangle[1], 0.0, triangle[2]];
            let force = model.pressure_at(&traction, 0, bary).unwrap();
            let speed = model.pressure_at(&velocity, 0, bary).unwrap();
            work += 0.25 * weight * (force.conj() * speed).re; // area=1/2, phasor factor=1/2
        }
    }
    real_close(weak.response.wall_accounting.wall_traction_power_w, work);
    assert_eq!(weak.response.acoustic.input_power_w, 0.0);
    assert!(weak.response.wall_accounting.wall_traction_power_w > 0.0);
    assert!(weak.response.wall_accounting.fluid_to_wall_power_w < 0.0);
    assert!(weak.response.walls[0].displacement_m[0].abs() > 1e-5);
    real_close(weak.response.power_balance_residual_w(), 0.0);
    real_close(weak.response.fluid_power_balance_residual_w(), 0.0);
    real_close(
        weak.response.wall_accounting.power_balance_residual_w(),
        0.0,
    );
}

#[test]
fn applied_force_crosses_a_local_lossless_impedance_zero_with_finite_tangent() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POINTS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let wall = MechanicalBoundary {
        faces: vec![WALL],
        wall: WallPin {
            surface_density: 1.0,
            stiffness_per_area: 1.0,
            resistance: 0.0,
        },
        relaxation: vec![],
    };
    let prepared = prepare(&model, 1.0, &wall, &[]);
    let source = prepared
        .prepare_wall_traction(&profile(), &mut || false)
        .unwrap();
    let amplitude = C64::new(0.7, -0.2);
    let slope = C64::new(-0.03, 0.04);
    let result = prepared
        .solve_with_wall_tractions_frequency_derivative(
            &[],
            None,
            &[],
            &[],
            &[(&source, amplitude, slope)],
            &mut || false,
        )
        .unwrap();
    for (local, &global) in prepared.wall_dofs(0).unwrap().iter().enumerate() {
        // At K-sigma*omega^2=0, original momentum fixes p=-traction;
        // differentiated momentum gives p'=-2*w-traction'.
        close(
            result.response.acoustic.pressure_pa[global],
            -amplitude.scale(SHAPE[local]),
        );
        close(
            result.pressure_derivative_pa_per_rad_s[global],
            -result.response.walls[0].displacement_m[local].scale(2.0) - slope.scale(SHAPE[local]),
        );
    }
    // Central difference: truncation ~ h^2 |f'''| / 6, roundoff ~ eps |f| / h (f = the solved field as a function of
    // omega). With h = 1e-5 the displacement block below (the pressure block, checked first, passed) missed the 2e-7
    // tolerance on the first run of this test. h = 3e-6 divides the truncation term by (1e-5 / 3e-6)^2 = 11 and
    // multiplies the roundoff term by 3.3; the tolerance is unchanged, and the planted slope error below must still
    // be refused by it.
    let h = 3e-6;
    let mut sides = Vec::new();
    for delta in [-h, h] {
        let at = prepare(&model, 1.0 + delta, &wall, &[]);
        let force = at.prepare_wall_traction(&profile(), &mut || false).unwrap();
        sides.push(
            at.solve_with_wall_tractions(
                &[],
                &[],
                &[],
                &[(&force, amplitude + slope.scale(delta))],
                &mut || false,
            )
            .unwrap(),
        );
    }
    for i in 0..model.ndof() {
        let difference =
            (sides[1].acoustic.pressure_pa[i] - sides[0].acoustic.pressure_pa[i]).scale(0.5 / h);
        assert!((difference - result.pressure_derivative_pa_per_rad_s[i]).abs() < 2e-7);
    }
    for i in 0..3 {
        let difference = (sides[1].walls[0].displacement_m[i]
            - sides[0].walls[0].displacement_m[i])
            .scale(0.5 / h);
        assert!((difference - result.wall_derivatives[0].displacement_m_per_rad_s[i]).abs() < 2e-7);
    }
    // Planted negative: the same comparison refuses a tangent whose traction slope is 0.1 % off.
    let off = prepared
        .solve_with_wall_tractions_frequency_derivative(
            &[],
            None,
            &[],
            &[],
            &[(&source, amplitude, slope.scale(1.001))],
            &mut || false,
        )
        .unwrap();
    assert!((0..model.ndof()).any(|i| {
        let difference =
            (sides[1].acoustic.pressure_pa[i] - sides[0].acoustic.pressure_pa[i]).scale(0.5 / h);
        (difference - off.pressure_derivative_pa_per_rad_s[i]).abs() >= 2e-7
    }));
    assert_eq!(result.response.wall_accounting.dissipation_w(), 0.0);
    real_close(result.response.power_balance_residual_w(), 0.0);
}

#[test]
fn preparation_binding_bad_spectra_and_cancellation_refuse_without_partial_fields() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POINTS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let prepared = prepare(&model, 3.0, &wall(), &[RELEASE]);
    let other = prepare(&model, 3.0, &wall(), &[RELEASE]);
    let source = prepared
        .prepare_wall_traction(&profile(), &mut || false)
        .unwrap();
    assert_eq!(source.wall(), 0);
    assert_eq!(source.coefficients(), SHAPE);
    for bad in [
        WallTractionProfile {
            wall: 1,
            coefficients: vec![],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(0, 1.0), (0, 2.0)],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(2, 0.0)],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(0, f64::NAN)],
        },
    ] {
        assert!(prepared.prepare_wall_traction(&bad, &mut || false).is_err());
    }
    assert!(
        other
            .solve_with_wall_tractions(&[], &[], &[], &[(&source, C64::ONE)], &mut || false)
            .is_err()
    );
    assert!(
        prepared
            .solve_with_wall_tractions(
                &[],
                &[],
                &[],
                &[(&source, C64::new(f64::INFINITY, 0.0))],
                &mut || false
            )
            .is_err()
    );
    assert!(
        prepared
            .solve_with_wall_tractions_frequency_derivative(
                &[],
                None,
                &[],
                &[],
                &[(&source, C64::ONE, C64::new(0.0, f64::NAN))],
                &mut || false
            )
            .is_err()
    );
    assert!(matches!(
        prepared.prepare_wall_traction(&profile(), &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert!(matches!(
        prepared.solve_with_wall_tractions(&[], &[], &[], &[(&source, C64::ONE)], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    let too_many = vec![(&source, C64::ONE); 129];
    assert!(matches!(
        prepared.solve_with_wall_tractions(&[], &[], &[], &too_many, &mut || false),
        Err(HeterogeneousError::Limit { .. })
    ));
    let zero = prepared
        .solve_with_wall_tractions(&[], &[], &[], &[], &mut || false)
        .unwrap();
    assert_eq!(zero.wall_accounting.wall_traction_power_w, 0.0);
    assert!(zero.acoustic.pressure_pa.iter().all(|&p| p == C64::ZERO));
}
