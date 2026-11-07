//! Signed separation paths composed from the original analytic contact law.
//!
//! The positive part of an affine signed indentation defines at most one
//! active interval. Shear is affine over the WHOLE interval, including the gap.
//! No attraction is invented here: an adhesive normal owner supplies its own
//! conservative force and work. Only Hertz repulsion loads these springs.
use super::{
    HertzMindlinContact, HertzMindlinState, MindlinPathError, MindlinPathStep, MindlinPathTrial,
    TangentialLaw, TangentialStepInput, TriboError, checked_add, checked_mul, finite, invalid,
    poll, positive_finite,
};

/// Allocation-free whole-path response with exact contact timing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MindlinGapPathTrial {
    pub tangential: MindlinPathTrial,
    /// Fraction with positive real contact area. Do not reconstruct this by
    /// subtracting interval endpoints: a short interval near one can round away.
    pub contact_fraction: f64,
    /// Contact interval in normalized whole-step time; [0, 0] if inactive.
    /// Its endpoint instants may have zero area.
    pub contact_interval_fraction: [f64; 2],
}

/// Candidate history for the complete contact/gap path.
#[derive(Clone, Debug, PartialEq)]
pub struct MindlinGapPathStep {
    pub path: MindlinPathStep,
    pub contact_fraction: f64,
    pub contact_interval_fraction: [f64; 2],
}

impl HertzMindlinContact {
    /// Probe simultaneous affine SIGNED indentation and scalar shear.
    ///
    /// The two signed endpoints must match the positive Hertz indentation in
    /// the accepted state and `input.normal_n`, respectively. Negative values
    /// describe a gap, never tensile Hertz pressure. Geometry is checked before
    /// use. This is the same decoupled scalar approximation as the parent path.
    ///
    /// For active duration fraction f, u_active=f*u_total and Q_total=f*Q_active.
    /// Therefore the original active solve uses f*free_total and f²*compliance.
    /// The returned tangent differentiates Q_total with respect to u_total.
    /// Force limits/tolerances apply to the WHOLE path average. Heat and storage
    /// are physical changes over the active interval and are never scaled by f.
    pub fn trial_coupled_gap_path(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        signed_indentation_m: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MindlinGapPathTrial, MindlinPathError> {
        let window = self.gap_window(state, input, signed_indentation_m, cancelled)?;
        let tangential = if !window.crosses() {
            self.trial_coupled_path(state, input, cancelled)?
        } else {
            let active_input = window.active_input(input)?;
            let active = self.trial_coupled_path(state, &active_input, cancelled)?;
            window.whole_trial(state, input, active)?
        };
        poll(cancelled)?;
        Ok(MindlinGapPathTrial {
            tangential,
            contact_fraction: window.active,
            contact_interval_fraction: window.interval,
        })
    }

    /// Build one immutable history proposal for the same signed affine path.
    ///
    /// The original contact solver runs only on its actual active interval.
    /// Gap drift occurs in chronological order at exactly zero shear load;
    /// departure releases the original memory before that drift and reentry
    /// grows unstrained material after it. Only final history construction
    /// allocates, retaining the original caller-specified memory ceiling.
    pub fn advance_coupled_gap_path(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        signed_indentation_m: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MindlinGapPathStep, MindlinPathError> {
        let window = self.gap_window(state, input, signed_indentation_m, cancelled)?;
        let path = if !window.crosses() {
            self.advance_coupled_path(state, input, cancelled)?
        } else {
            let mut active_input = window.active_input(input)?;
            let active = self.trial_coupled_path(state, &active_input, cancelled)?;
            let whole = window.whole_trial(state, input, active)?;
            let mut gap_input = *input;
            gap_input.normal_n = 0.0;
            gap_input.compliance_m_per_n = 0.0;
            gap_input.dt_s = scaled(input.dt_s, window.gap, "mindlin_gap_duration_s")?;
            gap_input.slip_increment_m[0] =
                scaled(whole.slip_increment_m, window.gap, "mindlin_gap_slip_m")?;
            active_input.slip_increment_m[0] = active.slip_increment_m;
            active_input.compliance_m_per_n = 0.0;
            let entered;
            let contact_state = if window.gap_first {
                entered = self.advance(state, &gap_input)?.state;
                &entered
            } else {
                state
            };
            let mut contact =
                self.advance_prescribed_path(contact_state, &active_input, cancelled)?;
            if !window.gap_first {
                contact.tangential.state =
                    self.advance(&contact.tangential.state, &gap_input)?.state;
            }
            // Both real segments have advanced. Store the single macroscopic
            // endpoint without split-addition roundoff accumulating over cycles.
            // This changes no strain, force, energy or reversal interval.
            let slip = checked_add(
                state.slip_m,
                whole.slip_increment_m,
                "mindlin_accumulated_slip_m",
            )?;
            contact.tangential.state.slip_m = slip;
            if let Some(fixed) = &mut contact.tangential.state.fixed {
                fixed.slip_m = slip;
            }
            contact.tangential.traction_n[0] = scaled(
                contact.tangential.traction_n[0],
                window.active,
                "mindlin_gap_average_n",
            )?;
            check_force(contact.tangential.traction_n[0], input)?;
            contact.tangential.slip_increment_m[0] = whole.slip_increment_m;
            contact.tangential.compatibility_residual_m[0] = compatibility(
                input,
                whole.slip_increment_m,
                contact.tangential.traction_n[0],
            )?;
            contact.tangential.closure_residual_j = closure(
                contact.tangential.stored_change_j,
                contact.tangential.dissipated_j,
                contact.tangential.traction_n[0],
                whole.slip_increment_m,
            )?;
            contact.tangential.iterations = active.iterations;
            contact
        };
        poll(cancelled)?;
        Ok(MindlinGapPathStep {
            path,
            contact_fraction: window.active,
            contact_interval_fraction: window.interval,
        })
    }

    fn gap_window(
        &self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        signed: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Window, MindlinPathError> {
        self.admit_path(state, input, false, cancelled)?;
        for d in signed {
            finite(d, "mindlin_signed_indentation_m")?;
        }
        same_geometry(state.indentation_m, signed[0].max(0.0))?;
        let end = if input.normal_n == 0.0 {
            0.0
        } else {
            self.normal_geometry(input.normal_n)?.1
        };
        same_geometry(end, signed[1].max(0.0))?;
        Window::new(signed)
    }
}

#[derive(Clone, Copy)]
struct Window {
    active: f64,
    gap: f64,
    interval: [f64; 2],
    gap_first: bool,
}

impl Window {
    fn new(d: [f64; 2]) -> Result<Self, MindlinPathError> {
        if d[0] <= 0.0 && d[1] <= 0.0 {
            return Ok(Self {
                active: 0.0,
                gap: 1.0,
                interval: [0.0; 2],
                gap_first: false,
            });
        }
        if d[0] >= 0.0 && d[1] >= 0.0 {
            return Ok(Self {
                active: 1.0,
                gap: 0.0,
                interval: [0.0, 1.0],
                gap_first: false,
            });
        }
        let gap_first = d[0] < 0.0;
        let positive = d[0].max(d[1]);
        let negative = -d[0].min(d[1]);
        let scale = positive.max(negative);
        let denominator = positive / scale + negative / scale;
        // Compute BOTH durations directly, avoiding 1-f cancellation.
        let active = (positive / scale) / denominator;
        let gap = (negative / scale) / denominator;
        positive_finite(active, "mindlin_contact_fraction")?;
        positive_finite(gap, "mindlin_gap_fraction")?;
        Ok(Self {
            active,
            gap,
            gap_first,
            interval: if gap_first { [gap, 1.0] } else { [0.0, active] },
        })
    }

    fn crosses(self) -> bool {
        self.active > 0.0 && self.gap > 0.0
    }

    fn active_input(
        self,
        input: &TangentialStepInput,
    ) -> Result<TangentialStepInput, MindlinPathError> {
        let mut active = *input;
        active.slip_increment_m[0] = scaled(
            input.slip_increment_m[0],
            self.active,
            "mindlin_active_free_slip_m",
        )?;
        active.compliance_m_per_n = scaled(
            scaled(
                input.compliance_m_per_n,
                self.active,
                "mindlin_active_compliance_m_per_n",
            )?,
            self.active,
            "mindlin_active_compliance_m_per_n",
        )?;
        active.dt_s = scaled(input.dt_s, self.active, "mindlin_active_duration_s")?;
        active.solve.maximum_force_n = input.solve.maximum_force_n / self.active;
        active.solve.force_absolute_tolerance_n =
            input.solve.force_absolute_tolerance_n / self.active;
        active.solve.validate()?;
        Ok(active)
    }

    fn whole_trial(
        self,
        state: &HertzMindlinState,
        input: &TangentialStepInput,
        active: MindlinPathTrial,
    ) -> Result<MindlinPathTrial, MindlinPathError> {
        let increment = if input.compliance_m_per_n == 0.0 {
            input.slip_increment_m[0]
        } else {
            active.slip_increment_m / self.active
        };
        finite(increment, "mindlin_gap_total_slip_m")?;
        checked_add(state.slip_m, increment, "mindlin_accumulated_slip_m")?;
        let average = scaled(active.average_force_n, self.active, "mindlin_gap_average_n")?;
        check_force(average, input)?;
        let derivative = scaled(
            scaled(
                active.average_tangent_n_per_m,
                self.active,
                "mindlin_gap_tangent_n_per_m",
            )?,
            self.active,
            "mindlin_gap_tangent_n_per_m",
        )?;
        Ok(MindlinPathTrial {
            slip_increment_m: increment,
            average_force_n: average,
            endpoint_force_n: if self.gap_first {
                active.endpoint_force_n
            } else {
                0.0
            },
            average_tangent_n_per_m: derivative,
            stored_change_j: active.stored_change_j,
            dissipated_j: active.dissipated_j,
            closure_residual_j: closure(
                active.stored_change_j,
                active.dissipated_j,
                average,
                increment,
            )?,
            compatibility_residual_m: compatibility(input, increment, average)?,
            iterations: active.iterations,
        })
    }
}

fn same_geometry(actual: f64, expected: f64) -> Result<(), MindlinPathError> {
    let scale = actual.abs().max(expected.abs());
    if (actual == 0.0) != (expected == 0.0)
        || (actual - expected).abs() > 512.0 * f64::EPSILON * scale
    {
        return Err(invalid("mindlin_signed_path_geometry_mismatch").into());
    }
    Ok(())
}

// A nonzero physical segment or response may not silently disappear when
// transformed into the active interval's units.
fn scaled(value: f64, fraction: f64, field: &'static str) -> Result<f64, MindlinPathError> {
    let result = checked_mul(value, fraction, field)?;
    if value != 0.0 && result == 0.0 {
        return Err(invalid(field).into());
    }
    Ok(result)
}
fn check_force(force: f64, input: &TangentialStepInput) -> Result<(), MindlinPathError> {
    if force.abs() > input.solve.maximum_force_n {
        return Err(TriboError::OutsideApplicability {
            field: "maximum_force_n",
            value: force.abs(),
            minimum: 0.0,
            maximum: input.solve.maximum_force_n,
        }
        .into());
    }
    Ok(())
}
fn compatibility(
    input: &TangentialStepInput,
    slip: f64,
    force: f64,
) -> Result<f64, MindlinPathError> {
    let residual = slip
        - (input.slip_increment_m[0]
            - checked_mul(input.compliance_m_per_n, force, "mindlin_body_response_m")?);
    finite(residual, "mindlin_compatibility_residual_m")?;
    Ok(residual)
}
fn closure(storage: f64, heat: f64, force: f64, slip: f64) -> Result<f64, MindlinPathError> {
    let residual = storage + heat - checked_mul(force, slip, "mindlin_path_shear_work_j")?;
    finite(residual, "mindlin_path_closure_j")?;
    Ok(residual)
}
