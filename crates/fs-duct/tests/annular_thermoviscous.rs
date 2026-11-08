//! G0/G1/G3: actual annular diffusion, independent physical limits and h ladders.
//! Numerical material fixtures and manufactured radial references establish
//! implementation behavior, not measured microgap/entrance qualification.

use fs_duct::compose::{ComposeOptions, compose_chain};
use fs_duct::fluid::annular::{AnnularAdmission, AnnularSection, MAX_ANNULAR_ELEMENTS};
use fs_duct::fluid::{FluidDuctError, segment_wave_with_loss};
use fs_duct::{DuctError, LossModel};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use std::cell::Cell;
use std::f64::consts::PI;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn admission() -> AnnularAdmission {
    AnnularAdmission {
        maximum_transverse_phase: 0.8,
        maximum_cell_diffusion_phase: 0.75,
        maximum_relative_residual: 1e-8,
        maximum_payload_bytes: 2 << 20,
    }
}
fn options() -> ComposeOptions {
    ComposeOptions {
        preserve_well_conditioned_bits: false,
        ..ComposeOptions::default()
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-25),
        "{actual:?} != {expected:?} at tolerance {tolerance}"
    );
}
fn liquid(conductivity: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4200.0, conductivity, 0.001).unwrap(),
        4200.0,
    )
    .unwrap()
}

#[test]
fn both_annular_walls_recover_poiseuille_resistance_and_isothermal_compression() {
    let gas = air();
    let (inner, outer, length, omega) = (0.001_f64, 0.003_f64, 0.03, 0.001);
    let section = AnnularSection::new(inner, outer, 512).unwrap();
    assert_eq!(section.inner_radius_m(), inner);
    assert_eq!(section.outer_radius_m(), outer);
    assert_eq!(section.radius_m(0), Some(inner));
    assert_eq!(section.radius_m(512), Some(outer));
    assert_eq!(section.radius_m(513), None);
    let area = PI * (outer - inner) * (outer + inner);
    assert!((section.area_m2() / area - 1.0).abs() < 5e-16);
    let response = section.evaluate(&gas, omega, admission()).unwrap();
    assert_eq!(response.area_m2(), section.area_m2());
    let block = response
        .two_port(length, gas.characteristic_impedance / area)
        .unwrap();
    // Integrate the exact no-slip profile
    // u(r)=G/(4 mu) [b²-r²-(b²-a²) log(b/r)/log(b/a)] independently.
    let geometric = outer.powi(4)
        - inner.powi(4)
        - (outer * outer - inner * inner).powi(2) / (outer / inner).ln();
    let resistance = 8.0 * gas.dynamic_viscosity * length / (PI * geometric);
    let released = compose_chain(&[block], Some(C64::ZERO), options()).unwrap();
    assert!((released.z_in.re / resistance - 1.0).abs() < 7e-6);
    let capped = compose_chain(&[block], None, options()).unwrap();
    let isothermal_compliance = area * length * gas.gamma / (gas.density * gas.sound_speed.powi(2));
    assert!((capped.z_in.im * omega * isothermal_compliance - 1.0).abs() < 2e-7);
    let equal_area_circle = 8.0 * PI * gas.dynamic_viscosity * length / area.powi(2);
    assert!(resistance > 2.5 * equal_area_circle);
    for diffusion in [response.viscous, response.thermal.unwrap()] {
        assert!(diffusion.mean.im > 0.0 && diffusion.gradient_work > 0.0);
        assert!(diffusion.relative_work_defect < 1e-8);
        assert!(diffusion.relative_backward_residual < 1e-8);
        near(diffusion.mean + diffusion.wall, C64::ONE, 3e-16);
    }
}

// Independent shooting integration of h'' + gap/r h' = i*alpha*(h-1),
// h(0)=h(1)=0. A forced and a homogeneous primitive solution determine the
// unknown initial derivative. The mean is another integrated state. This
// uses uniform physical-radius RK4 panels, not production P1 matrices, radial
// mesh, coefficient assembly, tridiagonal elimination or mean quadrature.
fn radial_reference(a: f64, b: f64, omega: f64, diffusivity: f64, panels: usize) -> C64 {
    let gap = b - a;
    let mean_radius = 0.5 * (a + b);
    let alpha = omega * gap * gap / diffusivity;
    let derivative = |x: f64, y: [C64; 6]| {
        let r = a + gap * x;
        [
            y[1],
            C64::new(0.0, alpha) * (y[0] - C64::ONE) - y[1].scale(gap / r),
            y[0].scale(r / mean_radius),
            y[4],
            C64::new(0.0, alpha) * y[3] - y[4].scale(gap / r),
            y[3].scale(r / mean_radius),
        ]
    };
    let mut y = [C64::ZERO; 6];
    y[4] = C64::ONE;
    let step = 1.0 / panels as f64;
    for i in 0..panels {
        let x = i as f64 * step;
        let k1 = derivative(x, y);
        let k2 = derivative(
            x + 0.5 * step,
            core::array::from_fn(|j| y[j] + k1[j].scale(0.5 * step)),
        );
        let k3 = derivative(
            x + 0.5 * step,
            core::array::from_fn(|j| y[j] + k2[j].scale(0.5 * step)),
        );
        let k4 = derivative(x + step, core::array::from_fn(|j| y[j] + k3[j].scale(step)));
        for j in 0..6 {
            y[j] = y[j] + (k1[j] + k2[j].scale(2.0) + k3[j].scale(2.0) + k4[j]).scale(step / 6.0);
        }
    }
    y[2] - (y[0] / y[3]) * y[5]
}

#[test]
fn radial_h_refinement_converges_to_independent_primitive_diffusion() {
    let gas = air();
    let transport = gas.thermoviscous_state().unwrap().transport();
    let (a, b) = (0.0004, 0.0012);
    for omega in [100.0, 800.0, 3000.0] {
        let references = [
            transport.kinematic_viscosity_m2_s(),
            transport.thermal_diffusivity_m2_s(),
        ]
        .map(|d| {
            let reference = radial_reference(a, b, omega, d, 4096);
            near(reference, radial_reference(a, b, omega, d, 2048), 2e-10);
            reference
        });
        let mut previous = [f64::INFINITY; 2];
        for elements in [32, 64, 128] {
            let section = AnnularSection::new(a, b, elements).unwrap();
            assert_eq!(section.elements(), elements);
            let actual = section.evaluate(&gas, omega, admission()).unwrap();
            for (mechanism, response) in [actual.viscous, actual.thermal.unwrap()]
                .into_iter()
                .enumerate()
            {
                let error = (response.mean - references[mechanism]).abs();
                if previous[mechanism].is_finite() {
                    let ratio = previous[mechanism] / error;
                    assert!(
                        (3.8..4.2).contains(&ratio),
                        "omega={omega}, N={elements}, error ratio={ratio}"
                    );
                }
                previous[mechanism] = error;
                assert!(response.relative_work_defect < 1e-8);
                near(response.mean, references[mechanism], 0.003);
                near(response.wall, C64::ONE - references[mechanism], 0.02);
            }
        }
    }
}

fn slit_mean(gap: f64, omega: f64, diffusivity: f64) -> C64 {
    let z = C64::new(0.0, omega / diffusivity).sqrt().scale(0.5 * gap);
    let denominator = (2.0 * z.re).cosh() + (2.0 * z.im).cos();
    let tanh = C64::new((2.0 * z.re).sinh(), (2.0 * z.im).sin()).scale(1.0 / denominator);
    C64::ONE - tanh / z
}

#[test]
fn slit_and_disappearing_inner_wall_limits_keep_the_original_material() {
    let material = liquid(0.6);
    let (inner, outer) = (0.005, 0.00501);
    for omega in [10000.0, 100000.0] {
        let section = AnnularSection::new(inner, outer, 256).unwrap();
        let response = section.evaluate(&material, omega, admission()).unwrap();
        near(
            response.viscous.mean,
            slit_mean(outer - inner, omega, 1e-6),
            2e-5,
        );
        assert!(response.thermal.is_none());
        let changed = section.evaluate(&liquid(60.0), omega, admission()).unwrap();
        assert_eq!(response.wave.wavenumber, changed.wave.wavenumber);
        assert_eq!(
            response.wave.characteristic_impedance,
            changed.wave.characteristic_impedance
        );
    }
    // At low frequency an inner no-slip cylinder disappears only
    // logarithmically: do not silently substitute the circular response.
    let gas = air();
    let outer = 0.0012;
    let omega = 0.001;
    let circular = segment_wave_with_loss(&gas, outer, omega, LossModel::Bessel).unwrap();
    let circle_series =
        C64::new(0.0, -1.0) * circular.wavenumber * circular.characteristic_impedance;
    let mut previous = f64::INFINITY;
    for ratio in [0.1, 0.01, 0.0001] {
        let annular = AnnularSection::new(ratio * outer, outer, 512).unwrap();
        let response = annular.evaluate(&gas, omega, admission()).unwrap();
        let series =
            C64::new(0.0, -1.0) * response.wave.wavenumber * response.wave.characteristic_impedance;
        let error = (series - circle_series).abs() / circle_series.abs();
        assert!(error < previous && error > 0.05);
        previous = error;
    }
    assert!(previous < 0.13);
    assert!(AnnularSection::new(0.0, outer, 512).is_err());
}

#[test]
fn original_scattering_preserves_power_reciprocity_and_long_duct_composition() {
    let gas = air();
    let section = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    let response = section.evaluate(&gas, 3000.0, admission()).unwrap();
    let reference = gas.characteristic_impedance / section.area_m2();
    let whole = response.two_port(0.4, reference).unwrap();
    let quarter = response.two_port(0.1, reference).unwrap();
    let split = quarter
        .star(quarter, 1e12)
        .unwrap()
        .0
        .star(quarter, 1e12)
        .unwrap()
        .0
        .star(quarter, 1e12)
        .unwrap()
        .0;
    for (actual, expected) in whole.scattering.into_iter().zip(split.scattering) {
        near(actual, expected, 2e-12);
    }
    let [r, t, reverse, opposite] = whole.scattering;
    assert_eq!(r, opposite);
    assert_eq!(t, reverse);
    for phase in [0.0_f64, 0.7, 1.9, PI] {
        let other = C64::new(phase.cos(), phase.sin()).scale(0.7);
        let outgoing = (r + t * other).norm_sq() + (reverse + opposite * other).norm_sq();
        assert!(outgoing < 1.0 + other.norm_sq());
    }
    let matched = compose_chain(
        &[whole],
        Some(response.wave.characteristic_impedance),
        options(),
    )
    .unwrap();
    near(matched.z_in, response.wave.characteristic_impedance, 2e-12);
    let phase = response.wave.wavenumber.scale(0.4);
    near(
        matched.u_out_over_u_in,
        C64::new(phase.re.cos(), phase.re.sin()).scale((-phase.im).exp()),
        2e-12,
    );
    let opaque = response.two_port(100.0, reference).unwrap();
    assert!(
        opaque
            .scattering
            .iter()
            .all(|v| v.re.is_finite() && v.im.is_finite())
    );
    assert!(opaque.scattering[1].abs() < 1e-20);
}

struct OneSnapshot {
    state: ThermoviscousFluidState,
    queries: Cell<usize>,
}
impl FluidAcoustics for OneSnapshot {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        panic!("the consumer must use one complete original snapshot")
    }
}
impl FluidTransport for OneSnapshot {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("the consumer must not re-query separate transport")
    }
}
impl ThermoviscousFluid for OneSnapshot {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn one_snapshot_limits_actual_geometry_and_cancellation_remain_authoritative() {
    let gas = air();
    let material = OneSnapshot {
        state: gas.thermoviscous_state().unwrap(),
        queries: Cell::new(0),
    };
    let section = AnnularSection::new(0.001, 0.003, 128).unwrap();
    let exact_bytes = 192 * (section.elements() + 1) + 4096;
    let controls = AnnularAdmission {
        maximum_payload_bytes: exact_bytes,
        ..admission()
    };
    let response = section.evaluate(&material, 100.0, controls).unwrap();
    assert_eq!(material.queries.get(), 1);
    response.two_port(0.1, 1e6).unwrap();
    response.two_port(0.2, 1e6).unwrap();
    assert_eq!(material.queries.get(), 1);
    assert!(
        section
            .evaluate(
                &material,
                100.0,
                AnnularAdmission {
                    maximum_payload_bytes: exact_bytes - 1,
                    ..controls
                }
            )
            .is_err()
    );
    assert_eq!(material.queries.get(), 1);
    for (a, b, n) in [
        (0.0, 0.001, 2),
        (0.002, 0.001, 2),
        (f64::NAN, 0.003, 2),
        (0.001, f64::INFINITY, 2),
        (0.001, 0.003, 1),
        (0.001, 0.003, MAX_ANNULAR_ELEMENTS + 1),
    ] {
        assert!(AnnularSection::new(a, b, n).is_err());
    }
    assert!(AnnularSection::new(0.001, 0.003, MAX_ANNULAR_ELEMENTS).is_ok());
    assert!(section.evaluate(&gas, 1e6, admission()).is_err());
    assert!(
        AnnularSection::new(0.001, 0.003, 2)
            .unwrap()
            .evaluate(&gas, 3000.0, admission())
            .is_err()
    );
    for omega in [0.0, -1.0, f64::NAN] {
        assert!(section.evaluate(&gas, omega, admission()).is_err());
    }
    for invalid in [
        AnnularAdmission {
            maximum_transverse_phase: 1.0,
            ..admission()
        },
        AnnularAdmission {
            maximum_cell_diffusion_phase: 0.0,
            ..admission()
        },
        AnnularAdmission {
            maximum_relative_residual: f64::NAN,
            ..admission()
        },
    ] {
        assert!(section.evaluate(&gas, 100.0, invalid).is_err());
    }
    let mut invalid_gas = gas;
    invalid_gas.density = f64::NAN;
    assert!(matches!(
        section.evaluate(&invalid_gas, 100.0, admission()),
        Err(FluidDuctError::Material(_))
    ));
    assert!(matches!(
        AnnularSection::new_controlled(0.001, 0.003, 128, || true),
        Err(DuctError::Cancelled { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        section.evaluate_controlled(&material, 100.0, admission(), || {
            polls += 1;
            polls == 3
        }),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert_eq!(material.queries.get(), 2);
    // The smallest admissible mesh still uses both distinct boundary links.
    let smallest = AnnularSection::new(0.001, 0.003, 2).unwrap();
    let response = smallest.evaluate(&gas, 0.001, admission()).unwrap();
    assert!(response.viscous.gradient_work > 0.0 && response.viscous.relative_work_defect < 1e-8);
}
