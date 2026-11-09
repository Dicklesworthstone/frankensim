//! Bounded mean-property closure; retains the actual accepted source trial.

#[derive(Debug)]
pub(crate) struct MeanState<T> {
    pub value: T,
    pub mean: f64,
    pub residual: f64,
    pub iterations: usize,
}

#[derive(Debug, PartialEq)]
pub(crate) enum MeanError<E> {
    Callback(E),
    InvalidOutlet,
    Unresolved { iterations: usize, residual: f64 },
}

/// A passive source trial returns its owned state and outlet temperature.
/// Its maximum principle brackets (T_in+T_out)/2 between T_in and
/// (T_in+T_wall)/2 in either heat direction. Source continuity belongs to the
/// material owner; no monotonicity or unique nonlinear root is assumed.
pub(crate) fn solve<T, E>(
    inlet: f64,
    wall: f64,
    tolerance: f64,
    iterations: usize,
    mut poll: impl FnMut() -> Result<(), E>,
    mut trial: impl FnMut(f64) -> Result<(T, f64), E>,
) -> Result<MeanState<T>, MeanError<E>> {
    let far_mean = 0.5 * inlet + 0.5 * wall;
    let mut low = inlet.min(far_mean);
    let mut high = inlet.max(far_mean);
    let rounding = 8.0 * f64::EPSILON * inlet.max(wall);
    let mut residual = f64::INFINITY;
    for iteration in 1..=iterations {
        poll().map_err(MeanError::Callback)?;
        let mean = 0.5 * low + 0.5 * high;
        let (value, outlet) = trial(mean).map_err(MeanError::Callback)?;
        if !outlet.is_finite()
            || outlet < inlet.min(wall) - rounding
            || outlet > inlet.max(wall) + rounding
        {
            return Err(MeanError::InvalidOutlet);
        }
        residual = mean - (0.5 * inlet + 0.5 * outlet);
        if residual.abs() <= tolerance {
            return Ok(MeanState {
                value,
                mean,
                residual,
                iterations: iteration,
            });
        }
        if residual < 0.0 {
            low = mean;
        } else {
            high = mean;
        }
    }
    Err(MeanError::Unresolved {
        iterations,
        residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_both_heat_directions_retain_the_actual_nonlinear_trial() {
        for (inlet, wall) in [(300.0, 400.0), (400.0, 300.0)] {
            let target = inlet + 0.371 * (wall - inlet);
            let result = solve(
                inlet,
                wall,
                1e-11,
                64,
                || Ok::<_, ()>(()),
                |mean| {
                    let offset = (mean - target) / (wall - inlet);
                    let outlet = 2.0 * target - inlet + 0.2 * (wall - inlet) * offset * offset;
                    Ok((mean.to_bits(), outlet))
                },
            )
            .unwrap();
            assert!((result.mean - target).abs() <= 1e-11);
            assert!(result.residual.abs() <= 1e-11);
            assert_eq!(result.value, result.mean.to_bits());
            assert!(result.iterations > 1 && result.iterations <= 64);
        }
        let equilibrium = solve(
            350.0,
            350.0,
            1e-12,
            1,
            || Ok::<_, ()>(()),
            |mean| Ok((mean.to_bits(), 350.0)),
        )
        .unwrap();
        assert_eq!(equilibrium.mean, 350.0);
        assert_eq!(equilibrium.residual, 0.0);
    }

    #[test]
    fn g4_iteration_cancellation_and_source_refusals_are_bounded() {
        let mut trials = 0;
        let exhausted = solve(
            300.0,
            400.0,
            1e-15,
            3,
            || Ok::<_, &str>(()),
            |mean| {
                trials += 1;
                Ok((mean, 374.2))
            },
        );
        assert!(matches!(
            exhausted,
            Err(MeanError::Unresolved { iterations: 3, .. })
        ));
        assert_eq!(trials, 3);
        let cancelled = solve(
            300.0,
            400.0,
            1e-9,
            64,
            || Err("cancelled"),
            |_| -> Result<((), f64), _> { panic!("cancelled step evaluated a source") },
        );
        assert!(matches!(cancelled, Err(MeanError::Callback("cancelled"))));
        let source = solve(
            300.0,
            400.0,
            1e-9,
            64,
            || Ok(()),
            |_| -> Result<((), f64), _> { Err("outside source") },
        );
        assert!(matches!(source, Err(MeanError::Callback("outside source"))));
    }

    #[test]
    fn g0_nonpassive_and_nonfinite_trial_outlets_refuse() {
        for outlet in [299.0, 401.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                solve(
                    300.0,
                    400.0,
                    1e-9,
                    64,
                    || Ok::<_, ()>(()),
                    |_| Ok(((), outlet))
                ),
                Err(MeanError::InvalidOutlet)
            ));
        }
    }
}
