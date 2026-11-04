use super::*;
use fs_material::gas::GasSpec;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::WallPin;
use fs_thermochem::SpeciesId;

const WORK: usize = 100_000_000;
const PORT: PortId = PortId(7);

fn medium() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn duct(kind: usize) -> DgDuct {
    let mut duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: 0.05,
            }],
        },
        &medium(),
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Port(PORT),
        DgBoundary::Matched,
    )
    .unwrap();
    if kind == 1 || kind == 3 {
        duct = duct.with_viscothermal(8, [50.0, 2000.0], 0.02).unwrap();
    } else if kind == 2 || kind == 4 {
        duct = duct
            .with_spatial_viscothermal(DgWallResolution {
                poles_per_field: 8,
                quadrature_points: 2,
                band_rad_s: [50.0, 2000.0],
                maximum_relative_error: 0.02,
            })
            .unwrap();
    }
    if kind >= 3 {
        // Synthetic thermodynamic reservoirs exercise source equations; these
        // are not a fit or a claim of measured molecular relaxation parameters.
        let gas = medium();
        let species = SpeciesId::new("frequency-fixture").unwrap();
        let source = RelaxationAbsorption::new(vec![SpeciesRelaxation {
            species: species.clone(),
            processes: vec![RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: gas.pressure * 0.0003,
                    t_ref_k: gas.temperature,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            }],
        }])
        .unwrap()
        .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
        .unwrap();
        duct = duct
            .with_molecular_relaxation(source)
            .unwrap()
            .with_compliant_wall_relaxation(
                WallPin {
                    surface_density: 0.2,
                    stiffness_per_area: 8e7,
                    resistance: 30.0,
                },
                &GeneralizedMaxwell::new(2e9, vec![(8e8, 0.002)]).unwrap(),
            )
            .unwrap();
    }
    duct
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + 1e-22,
        "{actual:e} != {expected:e}"
    );
}

fn near_complex(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + 1e-22,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn analytic_frequency_tangents_preserve_all_native_material_channels() {
    let force = C64::new(1e-7, 4e-8);
    let force_slope = C64::new(3e-11, -2e-11);
    let step = 0.001;
    for kind in 0..5 {
        let duct = duct(kind);
        let unchanged = duct.zero_state();
        for omega in [650.0, 1050.0] {
            let result = duct
                .harmonic_response_with_frequency_derivative(
                    omega,
                    &[(PORT, force)],
                    Some(&[(PORT, force_slope)]),
                    WORK,
                    1e-11,
                    &mut || false,
                )
                .unwrap();
            assert_eq!(
                result.response,
                duct.harmonic_response(omega, &[(PORT, force)], WORK, 1e-11, &mut || false)
                    .unwrap()
            );
            let lo = duct
                .harmonic_response(
                    omega - step,
                    &[(PORT, force - force_slope.scale(step))],
                    WORK,
                    1e-11,
                    &mut || false,
                )
                .unwrap();
            let hi = duct
                .harmonic_response(
                    omega + step,
                    &[(PORT, force + force_slope.scale(step))],
                    WORK,
                    1e-11,
                    &mut || false,
                )
                .unwrap();
            let derivative = &result.derivative;
            for (actual, low, high) in [
                (
                    &derivative.energy_coordinates,
                    &lo.energy_coordinates,
                    &hi.energy_coordinates,
                ),
                (&derivative.pressure_pa, &lo.pressure_pa, &hi.pressure_pa),
                (
                    &derivative.volume_flow_m3_s,
                    &lo.volume_flow_m3_s,
                    &hi.volume_flow_m3_s,
                ),
            ] {
                assert_eq!(actual.len(), low.len());
                assert_eq!(actual.len(), high.len());
                for ((&value, &low), &high) in actual.iter().zip(low).zip(high) {
                    near_complex(value, (high - low).scale(0.5 / step), 3e-6);
                }
            }
            for (actual, low, high) in [
                (
                    derivative.mean_stored_energy_j,
                    lo.mean_stored_energy_j,
                    hi.mean_stored_energy_j,
                ),
                (
                    derivative.mean_wall_energy_j,
                    lo.mean_wall_energy_j,
                    hi.mean_wall_energy_j,
                ),
                (
                    derivative.wall_dissipation_w,
                    lo.wall_dissipation_w,
                    hi.wall_dissipation_w,
                ),
                (
                    derivative.molecular_dissipation_w,
                    lo.molecular_dissipation_w,
                    hi.molecular_dissipation_w,
                ),
                (
                    derivative.matched_absorption_w,
                    lo.matched_absorption_w,
                    hi.matched_absorption_w,
                ),
                (derivative.input_power_w, lo.input_power_w, hi.input_power_w),
            ] {
                near(actual, (high - low) * 0.5 / step, 3e-6);
            }
            near_complex(
                derivative.ports[0].pressure_pa,
                (hi.ports[0].pressure_pa - lo.ports[0].pressure_pa).scale(0.5 / step),
                3e-6,
            );
            assert_eq!(derivative.ports[0].inflow_m3_s, force_slope);
            near(
                derivative.ports[0].mean_input_power_w,
                derivative.input_power_w,
                1e-13,
            );
            near(
                derivative.wall_dissipation_w
                    + derivative.molecular_dissipation_w
                    + derivative.matched_absorption_w,
                derivative.input_power_w,
                2e-8,
            );
            assert!(derivative.solve_relative_residual <= 1e-11);
            if kind >= 3 {
                assert!(result.response.molecular_dissipation_w > 0.0);
                assert!(derivative.molecular_dissipation_w.abs() > 0.0);
            }
        }
        assert_eq!(duct.zero_state(), unchanged);
    }
}

#[test]
fn changing_flow_adds_its_exact_linear_tangent_and_zero_drive_keeps_zero_power_slope() {
    let duct = duct(3);
    let omega = 650.0;
    let flow = C64::new(1e-7, 4e-8);
    let slope = C64::new(3e-11, -2e-11);
    let fixed = duct
        .harmonic_response_with_frequency_derivative(
            omega,
            &[(PORT, flow)],
            None,
            WORK,
            1e-11,
            &mut || false,
        )
        .unwrap();
    let changing = duct
        .harmonic_response_with_frequency_derivative(
            omega,
            &[(PORT, flow)],
            Some(&[(PORT, slope)]),
            WORK,
            1e-11,
            &mut || false,
        )
        .unwrap();
    let linear = duct
        .harmonic_response(omega, &[(PORT, slope)], WORK, 1e-11, &mut || false)
        .unwrap();
    for ((&fixed, &changing), &linear) in fixed
        .derivative
        .energy_coordinates
        .iter()
        .zip(&changing.derivative.energy_coordinates)
        .zip(&linear.energy_coordinates)
    {
        near_complex(changing - fixed, linear, 2e-10);
    }
    let zero = duct
        .harmonic_response_with_frequency_derivative(
            omega,
            &[],
            Some(&[(PORT, slope)]),
            WORK,
            1e-11,
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        zero.derivative.energy_coordinates,
        linear.energy_coordinates
    );
    for value in [
        zero.derivative.mean_stored_energy_j,
        zero.derivative.mean_wall_energy_j,
        zero.derivative.wall_dissipation_w,
        zero.derivative.molecular_dissipation_w,
        zero.derivative.matched_absorption_w,
        zero.derivative.input_power_w,
    ] {
        assert_eq!(value, 0.0);
    }
}

#[test]
fn external_tangent_observation_checks_both_original_equations_and_source_admission() {
    let duct = duct(1);
    let flow = C64::new(1e-7, 4e-8);
    let result = duct
        .harmonic_response_with_frequency_derivative(
            650.0,
            &[(PORT, flow)],
            None,
            WORK,
            1e-11,
            &mut || false,
        )
        .unwrap();
    let view = duct.continuous_operator();
    assert_eq!(
        view.observe_harmonic_frequency_derivative(
            &result.response,
            &result.derivative.energy_coordinates,
            &[],
            WORK,
            1e-11,
            &mut || false
        )
        .unwrap(),
        result.derivative
    );
    let mut bad = result.response.clone();
    bad.energy_coordinates[0] = bad.energy_coordinates[0] + C64::from_re(1e-4);
    assert!(
        view.observe_harmonic_frequency_derivative(
            &bad,
            &result.derivative.energy_coordinates,
            &[],
            WORK,
            1e-11,
            &mut || false
        )
        .is_err()
    );
    assert!(
        view.observe_harmonic_frequency_derivative(
            &result.response,
            &vec![C64::ZERO; view.dimension()],
            &[],
            WORK,
            1e-11,
            &mut || false
        )
        .is_err()
    );
    assert!(
        duct.harmonic_response_with_frequency_derivative(
            650.0,
            &[(PORT, flow)],
            Some(&[(PortId(42), flow)]),
            WORK,
            1e-11,
            &mut || false
        )
        .is_err()
    );
    assert!(
        duct.harmonic_response_with_frequency_derivative(
            650.0,
            &[(PORT, flow)],
            Some(&[(PORT, C64::new(f64::NAN, 0.0))]),
            WORK,
            1e-11,
            &mut || false
        )
        .is_err()
    );
    assert!(
        duct.harmonic_response_with_frequency_derivative(
            650.0,
            &[(PORT, flow)],
            None,
            1,
            1e-11,
            &mut || false
        )
        .is_err()
    );
    assert!(
        duct.harmonic_response_with_frequency_derivative(
            1e9,
            &[(PORT, flow)],
            None,
            WORK,
            1e-11,
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        duct.harmonic_response_with_frequency_derivative(
            650.0,
            &[(PORT, flow)],
            None,
            WORK,
            1e-11,
            &mut || true
        ),
        Err(DuctError::Cancelled { .. })
    ));
}
