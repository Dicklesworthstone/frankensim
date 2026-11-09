use super::super::super::{
    power::{StructuralInteriorRelaxationPowerChannel as Channel, evaluate},
    sources::{
        StructuralInteriorRelaxationSourceOrder as Order,
        StructuralInteriorRelaxationSourceSpectrum,
    },
};
use super::*;
use crate::cavity_modal::interior::structure::resonance::{
    StructuralInteriorParameterDirection, StructuralInteriorSourceParameterJet,
    StructuralInteriorSources,
};
use core::cell::RefCell;
use fs_modal::harmonic::resonance::{self as query};
use fs_phs::mechanical::MechanicalRelaxationDirection;

fn controls() -> query::ResonanceControls {
    query::ResonanceControls {
        absolute_angular_tolerance_rad_s: 2e-12,
        relative_angular_tolerance: 2e-12,
        maximum_evaluations: 256,
    }
}
fn stiffness_direction<'p, 'm>(p: &'p Prepared<'m>) -> PreparedDirection<'p, 'm> {
    p.prepare_parameter_direction(
        StructuralInteriorRelaxationParameterDirection {
            spatial: StructuralInteriorParameterDirection {
                structural_mass: &[0.0],
                structural_damping: &[0.0],
                structural_stiffness: &[0.0],
                acoustic_mass: &[0.0; 16],
                acoustic_stiffness: &[0.0; 16],
                wall_damping: &[0.0; 16],
                coupling: &[0.0; 4],
            },
            branches: Some(&[MechanicalRelaxationDirection {
                projection: None,
                stiffness: 1.0,
                relaxation_time_s: 0.0,
            }]),
        },
        &mut || false,
    )
    .unwrap()
}

// Independent scalar original law: the uniform cavity compression makes the
// relaxed stiffness1, and the one Maxwell arm has k=1.25,tau=1. No native
// state, work contraction or query derivative is used by this oracle.
fn original_heat(w: f64) -> (f64, f64, f64) {
    let k = 1.25;
    let x = w * w;
    let d = (1.0 - x).powi(2) + x * (x - 1.0 - k).powi(2);
    let dx = 3.0 * x * x - 2.0 * (1.0 + 2.0 * k) * x + k * k + 2.0 * k - 1.0;
    let dk = 2.0 * x * (k + 1.0 - x);
    (
        0.5 * k * x / d,
        k * w * (d - x * dx) / (d * d),
        0.5 * x / d - 0.5 * k * x * dk / (d * d),
    )
}
fn bisect(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    let sign = f(lo).is_sign_positive();
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid).is_sign_positive() == sign {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn original_maxwell_heat_peak_bandwidth_and_material_sensitivity() {
    let cavity = Cavity::new();
    let branches = [branch()];
    let p = Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let direction = stiffness_direction(&p);
    let mut sources = StructuralInteriorSources {
        structural_force: &[C64::ONE],
        inward_flow_load: &[C64::ZERO; 4],
    };
    let spectrum = RefCell::new(&mut sources);
    let sample = |w| {
        evaluate(
            &p,
            Channel::TotalDissipation,
            &spectrum,
            w,
            Order::Frequency,
            None,
            &mut || false,
        )
        .map(|j| j.sample)
    };
    let jet = |w| {
        evaluate(
            &p,
            Channel::TotalDissipation,
            &spectrum,
            w,
            Order::Parameter,
            Some(&direction),
            &mut || false,
        )
    };
    let peak = query::refine_power_extremum_with_parameter_derivative_with_cancel(
        (0.8, 1.8),
        controls(),
        || false,
        sample,
        jet,
    )
    .unwrap();
    let expected_w = bisect(0.8, 1.8, |w| original_heat(w).1);
    let (height, _, height_d) = original_heat(expected_w);
    let frequency_d = expected_w / (2.0 * (3.0 * expected_w * expected_w - 3.5));
    near(peak.extremum.angular_frequency_rad_s, expected_w);
    near(peak.extremum.sample.average_power_w, height);
    near(peak.angular_frequency_derivative_rad_s, frequency_d);
    near(peak.average_power_derivative, height_d);
    let width = query::refine_power_bandwidth_with_parameter_derivative_with_cancel(
        &peak,
        (0.1, 1.3),
        (1.45, 4.0),
        controls(),
        || false,
        sample,
        jet,
    )
    .unwrap();
    let lower = bisect(0.1, expected_w, |w| original_heat(w).0 - 0.5 * height);
    let upper = bisect(expected_w, 4.0, |w| original_heat(w).0 - 0.5 * height);
    let edge_d = |w| {
        let (_, pw, pk) = original_heat(w);
        (0.5 * height_d - pk) / pw
    };
    let lower_d = edge_d(lower);
    let upper_d = edge_d(upper);
    near(width.bandwidth.lower_angular_frequency_rad_s, lower);
    near(width.bandwidth.upper_angular_frequency_rad_s, upper);
    near(width.lower_angular_frequency_derivative_rad_s, lower_d);
    near(width.upper_angular_frequency_derivative_rad_s, upper_d);
    near(width.angular_bandwidth_derivative_rad_s, upper_d - lower_d);
    near(
        width.response_quality_factor_derivative,
        frequency_d / (upper - lower) - expected_w * (upper_d - lower_d) / (upper - lower).powi(2),
    );
    // The power line shape is not the displacement line shape. Search the
    // actual displacement and independently differentiate (1+w²)/D(w²).
    let motion = search_resonance(
        (0.8, 1.8),
        ResonanceControls {
            maximum_evaluations: 256,
            absolute_frequency_tolerance_rad_s: 2e-12,
            relative_frequency_tolerance: 2e-12,
            measure_half_power_bandwidth: false,
        },
        || false,
        |w, _, _| {
            let r = p
                .factor(w, &mut || false)?
                .solve_with_frequency_derivative(&[C64::ONE], &[C64::ZERO; 4], None, &mut || {
                    false
                })?;
            Ok::<_, Error>((
                r.response.structural_displacement[0],
                r.structural_displacement_derivative_per_rad_s[0],
                r.relative_derivative_residual,
            ))
        },
    )
    .unwrap();
    let expected_motion = bisect(0.8, 1.8, |w| {
        let x = w * w;
        let d = (1.0 - x).powi(2) + x * (x - 2.25).powi(2);
        let dx = 3.0 * x * x - 7.0 * x + 3.0625;
        d - (1.0 + x) * dx
    });
    near(motion.peak.angular_frequency_rad_s, expected_motion);
    assert!(expected_motion < expected_w - 1e-3);
}

struct CoherentSpectrum {
    force: [[C64; 1]; 5],
    flow: [[C64; 4]; 5],
    calls: [usize; 3],
}
impl CoherentSpectrum {
    fn new() -> Self {
        Self {
            force: [[C64::ZERO; 1]; 5],
            flow: [[C64::ZERO; 4]; 5],
            calls: [0; 3],
        }
    }
}
impl StructuralInteriorRelaxationSourceSpectrum for CoherentSpectrum {
    fn jet(
        &mut self,
        w: f64,
        order: Order,
    ) -> Result<StructuralInteriorSourceParameterJet<'_>, Error> {
        self.calls[match order {
            Order::Frequency => 0,
            Order::Curvature => 1,
            Order::Parameter => 2,
        }] += 1;
        // A coherent affine force/flow spectrum with a nonzero phase slope
        // and independently authored parameter/mixed source change.
        self.force = [
            [C64::new(0.8, 0.11 * w)],
            [C64::new(0.0, 0.11)],
            [C64::ZERO],
            [C64::new(0.03 * w, -0.04)],
            [C64::from_re(0.03)],
        ];
        self.flow = [
            [C64::new(-0.02 * w, 0.01); 4],
            [C64::from_re(-0.02); 4],
            [C64::ZERO; 4],
            [C64::new(0.006, -0.009 * w); 4],
            [C64::new(0.0, -0.009); 4],
        ];
        let source = |i: usize| StructuralInteriorSources {
            structural_force: &self.force[i][..],
            inward_flow_load: &self.flow[i][..],
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
#[test]
fn coherent_source_power_survey_reuses_first_order_and_keeps_original_channels() {
    let mut cavity = Cavity::new();
    cavity.damping[0] = 0.07;
    for i in 1..5 {
        cavity.damping[i * 5 + i] = 0.015;
    }
    let branches = [branch()];
    let p = Prepared::new(cavity.source(), &branches, Limits::default(), &mut || false).unwrap();
    let direction = stiffness_direction(&p);
    let mut sources = CoherentSpectrum::new();
    let spectrum = RefCell::new(&mut sources);
    for w in [0.4, 0.9, 1.3, 2.0] {
        for channel in [
            Channel::StructuralInput,
            Channel::AcousticInput,
            Channel::FluidToStructure,
            Channel::StructuralDissipation,
            Channel::WallDissipation,
            Channel::RelaxationDissipation,
            Channel::StructuralToRelaxation,
            Channel::AcousticToRelaxation,
            Channel::TotalDissipation,
            Channel::BranchDissipation(0),
            Channel::BranchStructuralWork(0),
            Channel::BranchAcousticWork(0),
        ] {
            let first = evaluate(
                &p,
                channel,
                &spectrum,
                w,
                Order::Frequency,
                None,
                &mut || false,
            )
            .unwrap();
            let full = evaluate(
                &p,
                channel,
                &spectrum,
                w,
                Order::Parameter,
                Some(&direction),
                &mut || false,
            )
            .unwrap();
            near(first.sample.average_power_w, full.sample.average_power_w);
            near(
                first.sample.angular_derivative,
                full.sample.angular_derivative,
            );
        }
    }
    let before = spectrum.borrow().calls;
    let grid = [0.3, 0.6, 0.9, 1.1, 1.3, 1.5, 1.8, 2.2, 3.0];
    let survey = query::survey_power_extrema_with_parameter_derivative_with_cancel(
        &grid,
        controls(),
        4,
        || false,
        |w| {
            evaluate(
                &p,
                Channel::TotalDissipation,
                &spectrum,
                w,
                Order::Frequency,
                None,
                &mut || false,
            )
            .map(|j| j.sample)
        },
        |w| {
            evaluate(
                &p,
                Channel::TotalDissipation,
                &spectrum,
                w,
                Order::Parameter,
                Some(&direction),
                &mut || false,
            )
        },
    )
    .unwrap();
    assert!(!survey.extrema.is_empty());
    let after = spectrum.borrow().calls;
    assert_eq!(after[1], before[1]);
    assert_eq!(after[2] - before[2], survey.extrema.len());
    assert_eq!(
        after[0] - before[0] + after[2] - before[2],
        survey.evaluations
    );
    assert!(matches!(
        evaluate(
            &p,
            Channel::BranchDissipation(1),
            &spectrum,
            1.0,
            Order::Frequency,
            None,
            &mut || false
        ),
        Err(Error::Invalid(_))
    ));
    let before = spectrum.borrow().calls;
    assert!(matches!(
        evaluate(
            &p,
            Channel::TotalDissipation,
            &spectrum,
            1.0,
            Order::Frequency,
            None,
            &mut || true
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(spectrum.borrow().calls, before);
}
