//! Second-order time integration of the original radial Kirchhoff operator.
//!
//! The same Alexander tableau as the slab/Cartesian owners acts on physical
//! shell storage. The explicit load is assembled from the original first-stage
//! face rates, normalized by the common dry-mass scale, not from a difference
//! of rounded cell moisture or total water. Both implicit stages use the native
//! radial nonlinear solve, including its material and conservation admission.

use super::*;
use crate::moisture::sdirk::{GAMMA, StageLoad};

impl MoistureRadial {
    /// Propose a second-order SDIRK interval on the original curved shell mesh.
    ///
    /// Both stages retain the original source-integrated diffusivity, actual
    /// shell dry masses and curved-face surface resistance. The first-stage
    /// flux contributes with weight `(1-gamma)*dt`; the final implicit flux has
    /// weight `gamma*dt`. The existing report separates the diagnostic first
    /// stage from the complete original-to-final water account.
    ///
    /// Temperature and both reservoir conditions stay fixed over this call;
    /// cut intervals at environmental events. L-stability does not guarantee
    /// monotonicity: an inadmissible second stage returns `StageRange` and
    /// needs a shorter interval. No state is clipped or replaced by backward
    /// Euler. Ordinary [`Self::advance`] keeps its original first-order method.
    /// No global temporal-error bound is inferred from a converged stage.
    pub fn advance_sdirk2(
        &self,
        dt_s: f64,
        faces: &[MoistureFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureSdirk2Report), MoistureError> {
        poll(cx)?;
        nonlinear.validate()?;
        let implicit_dt = GAMMA * dt_s;
        let first_time = self.time_s + implicit_dt;
        let end = self.time_s + dt_s;
        if !positive(dt_s)
            || !positive(implicit_dt)
            || !end.is_finite()
            || first_time <= self.time_s
            || first_time >= end
        {
            return Err(MoistureError::Invalid(
                "radial SDIRK needs two representable forward stage times",
            ));
        }
        let temperature = self.temperature_kelvin();
        let boundaries = [
            faces[0].resolve(temperature)?,
            faces[1].resolve(temperature)?,
        ];
        let (first, first_stage) = self.advance(implicit_dt, faces, step, nonlinear, cx)?;
        let known = self.sdirk_load(dt_s, &first, &boundaries, cx)?;
        let (next, endpoint) = self.advance_stage(
            dt_s,
            implicit_dt,
            boundaries,
            step,
            nonlinear,
            Some(&known),
            cx,
        )?;
        poll(cx)?;
        Ok((
            next,
            MoistureSdirk2Report {
                endpoint,
                first_stage,
            },
        ))
    }

    /// Original-shell first-stage flux for fixed or source-owned reservoirs.
    /// The load uses normalized physical water, so unequal shell masses remain
    /// in the existing stage operator instead of becoming fictitious equal cells.
    pub(in crate::moisture) fn sdirk_load(
        &self,
        dt_s: f64,
        first: &Self,
        boundaries: &[ResolvedMoistureFace; 2],
        cx: &Cx<'_>,
    ) -> Result<StageLoad, MoistureError> {
        let reference = diffusion_reference(&self.moisture, boundaries);
        let upper = self
            .moisture
            .iter()
            .copied()
            .chain(
                boundaries
                    .iter()
                    .filter(|face| face.coefficient_m_s != Some(0.0))
                    .filter_map(|face| face.equilibrium),
            )
            .fold(0.0, f64::max);
        let profile = KirchhoffMoistureDiffusivity::from_binding(
            &self.card,
            &self.transport,
            &self.property,
            [reference, upper],
        )?;
        let (operator, _, coupling) = self.system(
            (1.0 - GAMMA) * dt_s,
            boundaries,
            &profile,
            &first.moisture,
            reference,
            cx,
        )?;
        let mut increment = vec![0.0; self.moisture.len()];
        for (i, &coefficient) in operator.off.iter().enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            let flux = coefficient * (first.moisture[i + 1] - first.moisture[i]);
            increment[i] += flux;
            increment[i + 1] -= flux;
        }
        let mut lower = self.moisture.iter().copied().fold(f64::INFINITY, f64::min);
        let mut upper = self.moisture.iter().copied().fold(0.0, f64::max);
        for (side, i) in [0, self.moisture.len() - 1].into_iter().enumerate() {
            if coupling[side].0 != 0.0 {
                let equilibrium = boundaries[side].equilibrium.expect("active radial face");
                increment[i] += coupling[side].0 * (equilibrium - first.moisture[i]);
                // A sealed surface or zero diffusivity cannot widen the range.
                lower = lower.min(equilibrium);
                upper = upper.max(equilibrium);
            }
        }
        if increment.iter().any(|value| !value.is_finite()) {
            return Err(MoistureError::Invalid("radial SDIRK face load overflow"));
        }
        let water = self.boundary_transfers(&first.moisture, boundaries, &coupling)?;
        poll(cx)?;
        Ok(StageLoad::new(increment, water, [lower, upper]))
    }
}
