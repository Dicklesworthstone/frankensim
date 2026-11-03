use super::*;
use crate::gradient_index::event::{BoundaryQuery, propagate_to_boundary_with_control};
use std::cell::Cell;

struct Uniform {
    absorption: f64,
    queries: Cell<u32>,
}
impl GradientIndexField for Uniform {
    fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        self.queries.set(self.queries.get() + 1);
        GradientIndexSample::try_new(1.3, Vec3::new(0.0, 0.0, 0.0), self.absorption, Some(2.0))
    }
    fn is_spatially_uniform(&self) -> bool {
        true
    }
}
struct Empty;
impl BoundaryOracle for Empty {
    fn query(&mut self, _: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        Ok(BoundaryQuery {
            clearance_m: f64::INFINITY,
            witness: None,
        })
    }
}
fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(1.0, 0.0, 0.0),
    }
}
fn options() -> GradientEventOptions {
    GradientEventOptions {
        integration: GradientIndexOptions {
            initial_step_m: 1.0,
            max_step_m: 1.0,
            relative_tolerance: 1e-11,
            ..GradientIndexOptions::default()
        },
        ..GradientEventOptions::default()
    }
}

// Invert SI Planck only to author a temperature with a known source function.
// The tested production path always uses the original Planck owner.
fn temperature_for_radiance(radiance: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0_f64;
    let k = 1.380_649e-23;
    let wavelength = 550e-9_f64;
    (h * c / (wavelength * k)) / (2.0 * h * c * c / wavelength.powi(5) * 1e-9 / radiance).ln_1p()
}

#[test]
fn g1_spatial_lte_source_controls_native_step_acceptance_and_reuses_prepared_optics() {
    let field = Uniform {
        absorption: 2.0,
        queries: Cell::new(0),
    };
    let temperatures = Cell::new(0_u32);
    let temperature = |p: Point3| {
        temperatures.set(temperatures.get() + 1);
        Ok(temperature_for_radiance(1.0 + 8.0 * p.x.powi(4)))
    };
    let run = |tolerance| {
        propagate_emission_to_boundary_with_control(
            &field,
            &temperature,
            tolerance,
            &mut Empty,
            ray(),
            Vec3::new(0.0, 1.0, 0.0),
            550.0,
            1.0,
            None,
            options(),
            || true,
        )
        .unwrap()
    };
    let (loose, approximate) = run(1e-2);
    field.queries.set(0);
    temperatures.set(0);
    let (tight, refined) = run(1e-12);
    let survival = (-2.0_f64).exp();
    let i0 = 1.0 - survival;
    let mut i4 = i0;
    for degree in 1..=4 {
        i4 = f64::from(degree) * i4 / 2.0 - survival;
    }
    let exact = i0 + 8.0 * i4;
    assert!((refined - exact).abs() < 2e-10, "{refined} != {exact}");
    assert!((refined - exact).abs() < (approximate - exact).abs());
    assert!(tight.propagation.accepted_steps > loose.propagation.accepted_steps);
    assert!(tight.propagation.rejected_steps > 0);
    assert_eq!(
        field.queries.get(),
        1,
        "uniform spectral optics prepare once"
    );
    assert_eq!(
        temperatures.get(),
        1 + 6 * (tight.propagation.accepted_steps + tight.propagation.rejected_steps)
    );
    assert!((tight.propagation.optical_depth - 2.0).abs() < 1e-14);
    assert!((tight.propagation.group_optical_path_length_m.unwrap() - 2.0).abs() < 1e-14);
}

#[test]
fn g0_transparent_uniform_source_preserves_original_native_path_without_temperature_queries() {
    let field = Uniform {
        absorption: 0.0,
        queries: Cell::new(0),
    };
    let forbidden = |_: Point3| -> Result<f64, GradientIndexError> {
        panic!("zero true absorption cannot require temperature");
    };
    let expected = propagate_to_boundary_with_control(
        &field,
        &mut Empty,
        ray(),
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        1.0,
        None,
        options(),
        || true,
    )
    .unwrap();
    let (actual, emitted) = propagate_emission_to_boundary_with_control(
        &field,
        &forbidden,
        1e-12,
        &mut Empty,
        ray(),
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        1.0,
        None,
        options(),
        || true,
    )
    .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(emitted.to_bits(), 0.0_f64.to_bits());
}

#[test]
fn g0_source_domain_cancellation_and_native_step_budget_are_not_partial_radiance() {
    let field = Uniform {
        absorption: 2.0,
        queries: Cell::new(0),
    };
    let bounded = |p: Point3| {
        if p.x > 0.2 {
            Err(GradientIndexError::OutsideFieldDomain)
        } else {
            Ok(2400.0)
        }
    };
    let result = propagate_emission_to_boundary_with_control(
        &field,
        &bounded,
        1e-12,
        &mut Empty,
        ray(),
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        1.0,
        None,
        options(),
        || true,
    );
    assert_eq!(result, Err(GradientIndexError::OutsideFieldDomain));
    let forbidden = |_: Point3| -> Result<f64, GradientIndexError> {
        panic!("cancelled source queried");
    };
    let result = propagate_emission_to_boundary_with_control(
        &field,
        &forbidden,
        1e-12,
        &mut Empty,
        ray(),
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        1.0,
        None,
        options(),
        || false,
    );
    assert_eq!(result, Err(GradientIndexError::Cancelled));
    let mut limited = options();
    limited.integration.max_steps = 1;
    let temperature = |p: Point3| Ok(temperature_for_radiance(1.0 + 8.0 * p.x.powi(4)));
    let result = propagate_emission_to_boundary_with_control(
        &field,
        &temperature,
        1e-12,
        &mut Empty,
        ray(),
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        1.0,
        None,
        limited,
        || true,
    );
    assert_eq!(result, Err(GradientIndexError::StepBudgetExceeded));
}
