use super::super::super::sources::{
    StructuralInteriorRelaxationSourceOrder as Order, StructuralInteriorRelaxationSourceSpectrum,
};
use super::super::transfer::{
    StructuralInteriorRelaxationFixedReceiver, StructuralInteriorRelaxationReceiverJet,
    StructuralInteriorRelaxationReceiverSpectrum, evaluate,
};
use super::*;
use crate::cavity_modal::interior::structure::resonance::{
    StructuralInteriorParameterDirection, StructuralInteriorResponseField as Field,
    StructuralInteriorSourceParameterJet, StructuralInteriorSources,
};
use core::cell::RefCell;
use fs_modal::harmonic::resonance::{self as query, HarmonicReceiverParameterJet};
use fs_phs::mechanical::MechanicalRelaxationDirection;

const FORCE: [f64; 6] = [1.0, 0.1, 0.08, 0.04, 0.2, -0.06];
const RECEIVER: [f64; 6] = [1.0, 0.07, 0.03, 0.02, 0.1, 0.02];
// Original authored law (a0+aw*w+at*theta+awt*w*theta)*exp(i*w*(t0+tt*theta)).
fn law(w: f64, [a0, aw, at, awt, t0, tt]: [f64; 6]) -> [C64; 5] {
    let (a, da) = (a0 + aw * w, at + awt * w);
    let e = C64::new((w * t0).cos(), (w * t0).sin());
    [
        C64::from_re(a),
        C64::new(aw, t0 * a),
        C64::new(-t0 * t0 * a, 2.0 * t0 * aw),
        C64::new(da, w * tt * a),
        C64::new(awt - t0 * w * tt * a, tt * a + w * tt * aw + t0 * da),
    ]
    .map(|z| z * e)
}
fn product(a: [C64; 5], b: [C64; 5]) -> [C64; 5] {
    [
        a[0] * b[0],
        a[1] * b[0] + a[0] * b[1],
        a[2] * b[0] + (a[1] * b[1]).scale(2.0) + a[0] * b[2],
        a[3] * b[0] + a[0] * b[3],
        a[4] * b[0] + a[1] * b[3] + a[3] * b[1] + a[0] * b[4],
    ]
}
// Independent physical response: uniform compression adds 0.3 to K_s=0.7.
// k_theta=.31,tau_theta=.16,K_theta=.12; no matrix or native response helper.
fn original_jet(w: f64) -> [C64; 5] {
    let s = C64::new(0.0, -w);
    let sw = C64::new(0.0, -1.0);
    let t = C64::ONE + s;
    let d = C64::from_re(1.0 - w * w) + s.scale(1.25) / t;
    let dw = C64::from_re(-2.0 * w) + sw.scale(1.25) / (t * t);
    let dww = C64::from_re(-2.0) + C64::from_re(2.5) / (t * t * t);
    let dt = C64::from_re(0.12) + s.scale(0.31) / t + s.scale(1.25 * 0.16) / (t * t);
    let dwt = sw.scale(0.31 + 1.25 * 0.16) / (t * t) - (sw * s).scale(2.5 * 0.16) / (t * t * t);
    let inv = C64::ONE / d;
    let inverse = [
        inv,
        -dw * inv * inv,
        (dw * dw * inv * inv * inv).scale(2.0) - dww * inv * inv,
        -dt * inv * inv,
        (dw * dt * inv * inv * inv).scale(2.0) - dwt * inv * inv,
    ];
    product(product(law(w, FORCE), law(w, RECEIVER)), inverse)
}
struct Spectrum {
    constants: [f64; 6],
    values: [[C64; 4]; 5],
    field: Field,
    calls: [usize; 2],
    mismatch: bool,
}
impl Spectrum {
    fn new(constants: [f64; 6], field: Field) -> Self {
        Self {
            constants,
            values: [[C64::ZERO; 4]; 5],
            field,
            calls: [0; 2],
            mismatch: false,
        }
    }
    fn update(&mut self, w: f64, order: Order) {
        self.calls[usize::from(order == Order::Parameter)] += 1;
        self.values = law(w, self.constants).map(|v| [v; 4]);
    }
}
impl StructuralInteriorRelaxationSourceSpectrum for Spectrum {
    fn jet(
        &mut self,
        w: f64,
        order: Order,
    ) -> Result<StructuralInteriorSourceParameterJet<'_>, Error> {
        self.update(w, order);
        let source = |i: usize| StructuralInteriorSources {
            structural_force: &self.values[i][..1],
            inward_flow_load: &[C64::ZERO; 4],
        };
        Ok(StructuralInteriorSourceParameterJet {
            values: source(0),
            angular_derivative: Some(source(1)),
            angular_second_derivative: Some(source(2)),
            parameter_derivative: Some(source(3)),
            angular_parameter_derivative: Some(source(4)),
        })
    }
}
impl StructuralInteriorRelaxationReceiverSpectrum for Spectrum {
    fn jet(
        &mut self,
        w: f64,
        order: Order,
    ) -> Result<StructuralInteriorRelaxationReceiverJet<'_>, Error> {
        self.update(w, order);
        let n = if self.field == Field::Pressure { 4 } else { 1 };
        let scale = if self.mismatch && order == Order::Parameter {
            2.0
        } else {
            1.0
        } / n as f64;
        for row in &mut self.values {
            for z in row {
                *z = z.scale(scale);
            }
        }
        Ok(StructuralInteriorRelaxationReceiverJet {
            field: self.field,
            receiver: HarmonicReceiverParameterJet {
                weights: &self.values[0][..n],
                angular_derivative: Some(&self.values[1][..n]),
                angular_second_derivative: Some(&self.values[2][..n]),
                parameter_derivative: Some(&self.values[3][..n]),
                angular_parameter_derivative: Some(&self.values[4][..n]),
            },
        })
    }
}
fn direction<'p, 'm>(p: &'p Prepared<'m>) -> PreparedDirection<'p, 'm> {
    p.prepare_parameter_direction(
        StructuralInteriorRelaxationParameterDirection {
            spatial: StructuralInteriorParameterDirection {
                structural_mass: &[0.0],
                structural_damping: &[0.0],
                structural_stiffness: &[0.12],
                acoustic_mass: &[0.0; 16],
                acoustic_stiffness: &[0.0; 16],
                wall_damping: &[0.0; 16],
                coupling: &[0.0; 4],
            },
            branches: Some(&[MechanicalRelaxationDirection {
                projection: None,
                stiffness: 0.31,
                relaxation_time_s: 0.16,
            }]),
        },
        &mut || false,
    )
    .unwrap()
}
fn controls() -> query::ResonanceControls {
    query::ResonanceControls {
        absolute_angular_tolerance_rad_s: 2e-12,
        relative_angular_tolerance: 2e-12,
        maximum_evaluations: 256,
    }
}

#[test]
fn authored_pressure_receiver_source_and_material_five_order_jet_matches_original_law() {
    let cavity = Cavity::new();
    let branches = [branch()];
    let p = Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let d = direction(&p);
    let mut source = Spectrum::new(FORCE, Field::StructuralDisplacement);
    let mut receiver = Spectrum::new(RECEIVER, Field::Pressure);
    let (s, r) = (RefCell::new(&mut source), RefCell::new(&mut receiver));
    for w in [0.3, 0.85, 1.4, 2.0] {
        let full = evaluate(&p, &s, &r, w, Order::Parameter, Some(&d), &mut || false).unwrap();
        let first = evaluate(&p, &s, &r, w, Order::Frequency, None, &mut || false).unwrap();
        // Original uniform p=-0.6*u, including physical pressure recovery.
        for (a, b) in [
            full.sample.amplitude,
            full.sample.angular_derivative,
            full.angular_second_derivative,
            full.parameter_derivative,
            full.angular_parameter_derivative,
        ]
        .into_iter()
        .zip(original_jet(w))
        {
            close(a, b.scale(-0.6), 3e-9);
        }
        close(first.sample.amplitude, full.sample.amplitude, 3e-12);
        close(
            first.sample.angular_derivative,
            full.sample.angular_derivative,
            3e-12,
        );
    }
}

// A separate finite-parameter physical law for moving-feature verification.
fn response(w: f64, theta: f64) -> (C64, C64) {
    let value = |[a0, aw, at, awt, t0, tt]: [f64; 6]| {
        let a = a0 + aw * w + theta * (at + awt * w);
        let delay = t0 + theta * tt;
        let e = C64::new((w * delay).cos(), (w * delay).sin());
        (e.scale(a), C64::new(aw + theta * awt, delay * a) * e)
    };
    let (f, fw) = value(FORCE);
    let (r, rw) = value(RECEIVER);
    let (k, tau) = (1.25 + 0.31 * theta, 1.0 + 0.16 * theta);
    let s = C64::new(0.0, -w);
    let t = C64::ONE + s.scale(tau);
    let d = C64::from_re(1.0 + 0.12 * theta - w * w) + s.scale(k * tau) / t;
    let dw = C64::from_re(-2.0 * w) + C64::new(0.0, -k * tau) / (t * t);
    (f * r / d, (fw * r + f * rw) / d - f * r * dw / (d * d))
}
fn bisect(mut a: f64, mut b: f64, f: impl Fn(f64) -> f64) -> f64 {
    let positive = f(a) > 0.0;
    assert_ne!(positive, f(b) > 0.0);
    for _ in 0..80 {
        let m = 0.5 * (a + b);
        if (f(m) > 0.0) == positive {
            a = m;
        } else {
            b = m;
        }
    }
    0.5 * (a + b)
}
fn features(theta: f64) -> [f64; 4] {
    let peak = bisect(0.8, 1.8, |w| {
        let (h, hw) = response(w, theta);
        (h.conj() * hw).re
    });
    let level = response(peak, theta).0.abs().powi(2) * 0.5;
    let lower = bisect(0.01, peak, |w| response(w, theta).0.abs().powi(2) - level);
    let upper = bisect(peak, 4.0, |w| response(w, theta).0.abs().powi(2) - level);
    [peak, lower, upper, peak / (upper - lower)]
}
#[test]
fn receiver_peak_and_linewidth_follow_the_original_material_and_source_law() {
    let cavity = Cavity::new();
    let branches = [branch()];
    let p = Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let d = direction(&p);
    let mut source = Spectrum::new(FORCE, Field::StructuralDisplacement);
    let mut receiver = Spectrum::new(RECEIVER, Field::StructuralDisplacement);
    let (s, r) = (RefCell::new(&mut source), RefCell::new(&mut receiver));
    let first =
        |w| evaluate(&p, &s, &r, w, Order::Frequency, None, &mut || false).map(|j| j.sample);
    let full = |w| evaluate(&p, &s, &r, w, Order::Parameter, Some(&d), &mut || false);
    let peak = query::refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (0.8, 1.8),
        controls(),
        || false,
        first,
        full,
    )
    .unwrap();
    assert_eq!(s.borrow().calls[1], 1);
    assert_eq!(r.borrow().calls[1], 1);
    let pw = peak.extremum.angular_frequency_rad_s;
    let width = query::refine_half_power_bandwidth_with_parameter_derivative_with_cancel(
        &peak,
        (0.01, pw),
        (pw, 4.0),
        controls(),
        || false,
        first,
        full,
    )
    .unwrap();
    assert_eq!(s.borrow().calls[1], 3);
    assert_eq!(r.borrow().calls[1], 3);
    assert_eq!(
        s.borrow().calls.iter().sum::<usize>(),
        peak.evaluations + width.evaluations
    );
    assert_eq!(s.borrow().calls, r.borrow().calls);
    let actual = [
        pw,
        width.bandwidth.lower_angular_frequency_rad_s,
        width.bandwidth.upper_angular_frequency_rad_s,
        width.bandwidth.response_quality_factor,
    ];
    for (a, b) in actual.into_iter().zip(features(0.0)) {
        near(a, b);
    }
    let derivatives = [
        peak.angular_frequency_derivative_rad_s,
        width.lower_angular_frequency_derivative_rad_s,
        width.upper_angular_frequency_derivative_rad_s,
        width.response_quality_factor_derivative,
    ];
    let h = 2e-5;
    let plus = features(h);
    let minus = features(-h);
    for i in 0..4 {
        let expected = (plus[i] - minus[i]) / (2.0 * h);
        assert!(
            (derivatives[i] - expected).abs() < 2e-7,
            "order{i}:{} !={expected}",
            derivatives[i]
        );
    }
    near(
        width.angular_bandwidth_derivative_rad_s,
        derivatives[2] - derivatives[1],
    );
}

#[test]
fn receiver_shape_cancellation_and_changed_final_observation_refuse() {
    let cavity = Cavity::new();
    let branches = [branch()];
    let p = Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let d = direction(&p);
    let mut source = Spectrum::new(FORCE, Field::StructuralDisplacement);
    let s = RefCell::new(&mut source);
    let mut malformed = StructuralInteriorRelaxationFixedReceiver {
        field: Field::Pressure,
        weights: &[C64::ONE],
    };
    let r = RefCell::new(&mut malformed);
    assert!(matches!(
        evaluate(&p, &s, &r, 1.0, Order::Frequency, None, &mut || false),
        Err(Error::Invalid(_))
    ));
    let before = s.borrow().calls;
    assert!(matches!(
        evaluate(&p, &s, &r, 1.0, Order::Frequency, None, &mut || true),
        Err(Error::Cancelled)
    ));
    assert_eq!(before, s.borrow().calls);
    let mut receiver = Spectrum::new(RECEIVER, Field::StructuralDisplacement);
    receiver.mismatch = true;
    let r = RefCell::new(&mut receiver);
    let result = query::refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (0.8, 1.8),
        controls(),
        || false,
        |w| evaluate(&p, &s, &r, w, Order::Frequency, None, &mut || false).map(|j| j.sample),
        |w| evaluate(&p, &s, &r, w, Order::Parameter, Some(&d), &mut || false),
    );
    assert!(matches!(result, Err(query::ResonanceError::Invalid(_))));
}
