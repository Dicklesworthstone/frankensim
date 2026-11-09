use super::*;
use crate::d3q19::freesurface3::transport::MassRoute3;

fn nearly_empty_conversion(scale: f64) -> (MassTransport3, MassLocation3, MassLocation3, f64) {
    let source = MassLocation3::Cell(0, 0);
    let donor = MassLocation3::Cell(0, 1);
    let residual = 2.0f64.powi(-50) * scale;
    let transferred = 2.0f64.powi(-12) * scale;
    let mut trace = MassTransport3::new(BTreeMap::from([(source, residual), (donor, scale)]), 2);
    trace.after = BTreeMap::from([(source, 0.0), (donor, scale)]);
    trace.routes = vec![
        MassRoute3 {
            from: donor,
            to: source,
            mass: transferred,
        },
        MassRoute3 {
            from: source,
            to: donor,
            mass: transferred,
        },
    ];
    (trace, source, donor, residual)
}

#[test]
fn g1_conversion_density_roundoff_keeps_its_reference_energy_at_every_mass_scale() {
    let reference = 2.0f64.powi(40);
    let gate = CancelGate::new_clock_free();
    for scale in [2.0f64.powi(-60), 1.0, 2.0f64.powi(60)] {
        let (mut trace, source, donor, residual) = nearly_empty_conversion(scale);
        // The residual is below the native O(density) subtraction allowance,
        // but above roundoff measured against the tiny cell's route throughput.
        assert!(residual < 128.0 * f64::EPSILON * scale);
        assert!(residual > 2048.0 * f64::EPSILON * (trace.before[&source] + trace.routes[0].mass));
        trace.record_conversion_roundoff(source, residual).unwrap();
        let old_h = BTreeMap::from([(source, reference), (donor, reference)]);
        let (h, _) = advect(
            &trace,
            &old_h,
            FreeSurfaceEnthalpyControls::default(),
            &gate,
        )
        .unwrap();
        // The admitted mass defect remains in A*h = m_old*h_old. Removing it
        // from the centered RHS would return the unchanged reference instead.
        let expected_increment = residual / scale * reference;
        assert_eq!(h[&donor] - reference, expected_increment);
        assert_eq!(
            trace.after[&donor] * h[&donor],
            trace.before[&source] * old_h[&source] + trace.before[&donor] * old_h[&donor],
        );
    }
}

#[test]
fn g0_conversion_roundoff_is_local_finite_and_does_not_admit_unaccounted_defects() {
    let (mut trace, source, donor, residual) = nearly_empty_conversion(1.0);
    let old_h = BTreeMap::from([(source, 1.0), (donor, 1.0)]);
    let gate = CancelGate::new_clock_free();
    let controls = FreeSurfaceEnthalpyControls::default();
    let refuses = |trace: &MassTransport3| {
        assert!(matches!(
            advect(trace, &old_h, controls, &gate),
            Err(Error::Invalid(
                "native cell mass routes do not close to roundoff"
            ))
        ));
    };
    refuses(&trace);
    trace.record_conversion_roundoff(donor, residual).unwrap();
    refuses(&trace);
    trace
        .record_conversion_roundoff(source, residual * 0.5)
        .unwrap();
    refuses(&trace);
    // Independent conversion residuals accumulate by magnitude at their own
    // source. Opposite signs must not erase an already-admitted allowance.
    trace
        .record_conversion_roundoff(source, -residual * 0.5)
        .unwrap();
    assert_eq!(trace.conversion_roundoff_at(source), residual);
    advect(&trace, &old_h, controls, &gate).unwrap();
    *trace.before.get_mut(&source).unwrap() *= 2.0;
    refuses(&trace);
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(trace.record_conversion_roundoff(source, invalid).is_err());
        assert_eq!(trace.conversion_roundoff_at(source), residual);
    }
    let unrelated = MassLocation3::Cell(64, 0);
    trace
        .record_conversion_roundoff(unrelated, f64::MAX)
        .unwrap();
    assert!(
        trace
            .record_conversion_roundoff(unrelated, f64::MAX)
            .is_err()
    );
    assert_eq!(trace.conversion_roundoff_at(unrelated), f64::MAX);
    refuses(&trace);
}
