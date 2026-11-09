//! Shared pressure search on one original, ordinary sonic branch.
//!
//! Thermochemical owners supply their own entropy path, density and acoustics.
//! This owner only selects a back-pressure state or brackets c^2-v^2=0. The
//! source callback owns a strict query/work ceiling and all source admission.
//! Sampled monotonicity is an engineering check, not a uniqueness certificate.

pub(crate) const MAX_PRESSURE_EVALUATIONS: usize = 256;

/// The same scalar coordinates for frozen-species phase or reacting gas flow.
pub(crate) trait Station {
    fn pressure(&self) -> f64;
    fn velocity_squared(&self) -> f64;
    fn sonic(&self) -> f64;
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct PressureSearch {
    pub minimum_pressure_pa: f64,
    pub pressure_step_ratio: f64,
    pub absolute_sonic_tolerance_m2_s2: f64,
    pub relative_sonic_tolerance: f64,
    pub absolute_pressure_tolerance_pa: f64,
    pub relative_pressure_tolerance: f64,
}
impl PressureSearch {
    pub fn validate(self) -> Result<(), &'static str> {
        if !self.minimum_pressure_pa.is_finite()
            || self.minimum_pressure_pa <= 0.0
            || !(0.5..=0.95).contains(&self.pressure_step_ratio)
            || !self.absolute_sonic_tolerance_m2_s2.is_finite()
            || self.absolute_sonic_tolerance_m2_s2 <= 0.0
            || !self.absolute_pressure_tolerance_pa.is_finite()
            || self.absolute_pressure_tolerance_pa <= 0.0
            || [
                self.relative_sonic_tolerance,
                self.relative_pressure_tolerance,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || !(0.0..1.0).contains(&v))
        {
            return Err(
                "finite positive sonic/pressure allowances and connected search bounds required",
            );
        }
        Ok(())
    }

    pub fn sonic_tolerance(self, sample: &impl Station) -> f64 {
        let c2 = sample.velocity_squared() + sample.sonic();
        self.absolute_sonic_tolerance_m2_s2
            + self.relative_sonic_tolerance * c2.abs().max(sample.velocity_squared())
    }

    fn check_descent<E>(
        self,
        lower: &impl Station,
        upper: &impl Station,
    ) -> Result<(), SearchError<E>> {
        if lower.sonic() > upper.sonic() + self.sonic_tolerance(lower) + self.sonic_tolerance(upper)
        {
            return Err(SearchError::NonMonotone {
                pressure_bracket_pa: [lower.pressure(), upper.pressure()],
                sonic_residuals_m2_s2: [lower.sonic(), upper.sonic()],
            });
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum SearchError<E> {
    Invalid(&'static str),
    Source(E),
    NonMonotone {
        pressure_bracket_pa: [f64; 2],
        sonic_residuals_m2_s2: [f64; 2],
    },
    OutsidePressureBracket {
        minimum_pressure_pa: f64,
        back_pressure_pa: f64,
        sonic_residual_m2_s2: f64,
    },
    NotConverged {
        pressure_bracket_pa: [f64; 2],
        sonic_residual_m2_s2: f64,
    },
    Cancelled,
}

/// Original descending scan and safeguarded bracket refinement. The source
/// callback must bound every attempted query, including final re-publication
/// if its numerical workspace no longer contains the selected endpoint.
pub(crate) fn resolve<S: Station, E, C: FnMut() -> bool>(
    initial: S,
    back_pressure: f64,
    controls: PressureSearch,
    mut sample: impl FnMut(f64, &mut C) -> Result<S, E>,
    cancelled: &mut C,
) -> Result<(S, Option<[f64; 2]>), SearchError<E>> {
    controls.validate().map_err(SearchError::Invalid)?;
    let p0 = initial.pressure();
    if !p0.is_finite()
        || p0 <= 0.0
        || !back_pressure.is_finite()
        || back_pressure < 0.0
        || back_pressure >= p0
        || controls.minimum_pressure_pa >= p0
    {
        return Err(SearchError::Invalid(
            "descending positive-pressure nozzle lane required",
        ));
    }
    let floor = back_pressure.max(controls.minimum_pressure_pa);
    let mut upper = initial;
    let (mut lower, mut upper) = loop {
        if cancelled() {
            return Err(SearchError::Cancelled);
        }
        let high = upper.pressure();
        let next = (high * controls.pressure_step_ratio).max(floor);
        if !(next > 0.0 && next < high) {
            return Err(SearchError::NotConverged {
                pressure_bracket_pa: [next, high],
                sonic_residual_m2_s2: upper.sonic(),
            });
        }
        let trial = sample(next, cancelled).map_err(SearchError::Source)?;
        controls.check_descent(&trial, &upper)?;
        if trial.sonic() <= 0.0 {
            break (trial, upper);
        }
        if next == back_pressure {
            return Ok((trial, None));
        }
        if next == floor {
            return Err(SearchError::OutsidePressureBracket {
                minimum_pressure_pa: floor,
                back_pressure_pa: back_pressure,
                sonic_residual_m2_s2: trial.sonic(),
            });
        }
        upper = trial;
    };
    loop {
        if cancelled() {
            return Err(SearchError::Cancelled);
        }
        let bracket = [lower.pressure(), upper.pressure()];
        let pressure_tolerance =
            controls.absolute_pressure_tolerance_pa + controls.relative_pressure_tolerance * p0;
        let best = if lower.sonic().abs() < upper.sonic().abs() {
            &lower
        } else {
            &upper
        };
        if bracket[1] - bracket[0] <= pressure_tolerance
            && best.sonic().abs() <= controls.sonic_tolerance(best)
        {
            let accepted = if lower.sonic().abs() < upper.sonic().abs() {
                lower
            } else {
                upper
            };
            return Ok((accepted, Some(bracket)));
        }
        let mid = bracket[0] + 0.5 * (bracket[1] - bracket[0]);
        if !(mid > bracket[0] && mid < bracket[1]) {
            return Err(SearchError::NotConverged {
                pressure_bracket_pa: bracket,
                sonic_residual_m2_s2: best.sonic(),
            });
        }
        let trial = sample(mid, cancelled).map_err(SearchError::Source)?;
        controls.check_descent(&lower, &trial)?;
        controls.check_descent(&trial, &upper)?;
        if trial.sonic() <= 0.0 {
            lower = trial;
        } else {
            upper = trial;
        }
    }
}
