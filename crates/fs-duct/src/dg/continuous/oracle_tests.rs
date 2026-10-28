//! Independent cylinder fields and resonance references for the public solve.

use super::*;
use crate::dg::DgMesh;
use crate::{Duct, LossModel, Segment, segment_wave};
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;

const WORK: usize = 1_000_000_000;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn cylinder(elements: usize, order: usize, radius: f64, length: f64) -> DgDuct {
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        &air(),
        DgMesh { elements, order },
        DgBoundary::Port(PortId(7)),
        DgBoundary::IdealOpen,
    )
    .unwrap()
}

fn sin_cos(z: C64) -> (C64, C64) {
    let (sin, cos) = det::sin_cos(z.re);
    let positive = det::exp(z.im);
    let negative = 1.0 / positive;
    let sinh = 0.5 * (positive - negative);
    let cosh = 0.5 * (positive + negative);
    (
        C64::new(sin * cosh, cos * sinh),
        C64::new(cos * cosh, -sin * sinh),
    )
}

// Solve the uniform telegraph equations directly, with p(L)=0 and q(0)=U.
// This uses no DG matrix, mass map, eigensolver or diffusive realization.
// Under exp(-i omega t), p' = i k Zc q and q' = i k p/Zc.
fn fields(k: C64, zc: C64, length: f64, x: f64, inflow: C64) -> (C64, C64) {
    let denominator = sin_cos(k.scale(length)).1;
    let (sin, cos) = sin_cos(k.scale(length - x));
    (
        C64::new(0.0, -1.0) * zc * sin * inflow / denominator,
        cos * inflow / denominator,
    )
}

fn field_error(dg: &DgDuct, radius: f64, omega: f64, loss: LossModel) -> f64 {
    let inflow = C64::new(8e-8, -3e-8);
    let response = dg
        .harmonic_response(omega, &[(PortId(7), inflow)], WORK, 1e-10, &mut || false)
        .unwrap();
    let wave = segment_wave(dg.gas(), radius, omega, loss).unwrap();
    let denominator = sin_cos(wave.wavenumber.scale(dg.total_length())).1.abs();
    let flow_scale = inflow.abs() / denominator;
    let pressure_scale = wave.characteristic_impedance.abs() * flow_scale;
    let mut worst = 0.0_f64;
    for (i, x) in dg.node_positions().into_iter().enumerate() {
        let (pressure, flow) = fields(
            wave.wavenumber,
            wave.characteristic_impedance,
            dg.total_length(),
            x,
            inflow,
        );
        worst = worst.max((response.pressure_pa[i] - pressure).abs() / pressure_scale);
        worst = worst.max((response.volume_flow_m3_s[i] - flow).abs() / flow_scale);
    }
    let expected = fields(
        wave.wavenumber,
        wave.characteristic_impedance,
        dg.total_length(),
        0.0,
        inflow,
    )
    .0;
    worst = worst.max((response.ports[0].pressure_pa - expected).abs() / pressure_scale);
    assert_eq!(response.ports[0].id, PortId(7));
    assert_eq!(response.ports[0].inflow_m3_s, inflow);
    if loss == LossModel::Lossless {
        assert!(response.input_power_w.abs() < 1e-10 * pressure_scale * flow_scale);
    } else {
        assert!(response.input_power_w > 0.0);
        assert!(response.wall_dissipation_w > 0.0);
    }
    worst
}

#[test]
fn harmonic_cylinder_complex_fields_converge_with_spatial_resolution() {
    // G1: refine only the original DG space. Both complex field components
    // and the power-conjugate inlet pressure face the analytic cylinder.
    let length = 0.57;
    let first_pole = core::f64::consts::PI * air().sound_speed / (2.0 * length);
    let mut errors = Vec::new();
    for elements in [2, 4, 8] {
        let dg = cylinder(elements, 2, 0.0073, length);
        let error = [0.7, 1.4, 2.2]
            .into_iter()
            .map(|ratio| field_error(&dg, 0.0073, ratio * first_pole, LossModel::Lossless))
            .fold(0.0_f64, f64::max);
        errors.push(error);
    }
    assert!(
        errors[1] < 0.6 * errors[0],
        "spatial field errors: {errors:?}"
    );
    assert!(
        errors[2] < 0.6 * errors[1],
        "spatial field errors: {errors:?}"
    );
    assert!(errors[2] < 2e-3, "resolved complex field error: {errors:?}");
}

#[test]
fn harmonic_open_cylinder_pole_converges_to_the_analytic_quarter_wave() {
    // G2/G1: locate the zero of input admittance around the first resonance.
    // No native DG modal frequency is used as the reference or search center.
    let length = 0.57;
    let exact = core::f64::consts::PI * air().sound_speed / (2.0 * length);
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let dg = cylinder(elements, 2, 0.0073, length);
        let admittance = |omega| {
            let response = dg
                .harmonic_response(omega, &[(PortId(7), C64::ONE)], WORK, 1e-10, &mut || false)
                .unwrap();
            response.ports[0].pressure_pa.recip().im
        };
        let (mut lo, mut hi) = (0.8 * exact, 1.2 * exact);
        assert!(admittance(lo) > 0.0 && admittance(hi) < 0.0);
        for _ in 0..25 {
            let mid = 0.5 * (lo + hi);
            if admittance(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let numerical = 0.5 * (lo + hi);
        errors.push((1200.0 * det::ln(numerical / exact) / core::f64::consts::LN_2).abs());
    }
    assert!(
        errors[1] < errors[0],
        "quarter-wave pole errors in cents: {errors:?}"
    );
    assert!(
        errors[2] < errors[1],
        "quarter-wave pole errors in cents: {errors:?}"
    );
    assert!(
        errors[2] < 0.1,
        "quarter-wave pole errors in cents: {errors:?}"
    );
}

#[test]
fn harmonic_narrow_cylinder_complex_fields_converge_to_exact_bessel() {
    // G1: fixed high-order spatial space, independently refined constitutive
    // pole count. segment_wave uses the original exact Bessel coefficients,
    // while the production solve retains every actual DG memory coordinate.
    let mut errors = Vec::new();
    for poles in [4, 8, 16] {
        let dg = cylinder(1, 4, 0.001, 0.05)
            .with_viscothermal(poles, [80.0, 1200.0], 0.2)
            .unwrap();
        let error = [100.0, 500.0, 1000.0]
            .into_iter()
            .map(|frequency| {
                field_error(
                    &dg,
                    0.001,
                    core::f64::consts::TAU * frequency,
                    LossModel::Bessel,
                )
            })
            .fold(0.0_f64, f64::max);
        errors.push(error);
    }
    assert!(errors[1] < errors[0], "Bessel field errors: {errors:?}");
    assert!(errors[2] < errors[1], "Bessel field errors: {errors:?}");
    assert!(
        errors[2] < 2e-3,
        "resolved Bessel complex field error: {errors:?}"
    );
}
