//! Second-order time integration of the original three-dimensional Fick flux.
//!
//! Alexander's two-stage SDIRK method uses the same tableau and range admission
//! as the slab owner. Each stage retains the existing seven-point operator,
//! pinned principal material laws, nonlinear iteration and native CG solve.
//! Interior rates are computed as equal-and-opposite physical face transfers,
//! never inferred from subtraction of rounded stored-water inventories.

use super::*;
use crate::moisture::sdirk::{GAMMA, StageLoad};

/// One second-order spatial interval with both original implicit stage reports.
#[derive(Clone, Debug)]
pub struct MoistureCartesianSdirk2Report {
    /// Original-to-final interval and RK-weighted transfers through all six
    /// physical faces. Surface, source and solver diagnostics describe the
    /// final stage; `transport_start` describes the original accepted field.
    pub endpoint: MoistureCartesianReport,
    /// Internal stage at t0 + gamma*dt, reporting gamma*dt times its face rates.
    /// This is a diagnostic account, not an additional physical transfer.
    pub first_stage: MoistureCartesianReport,
}

impl MoistureCartesian {
    /// Propose a second-order SDIRK interval on all three principal directions.
    ///
    /// Both stages use the original nonlinear source laws, including Kirchhoff
    /// integrals for a degenerate dry front when selected. The explicit stage
    /// load uses actual first-stage face rates with weight (1-gamma)*dt; the
    /// final implicit flux has weight gamma*dt. No clipped field, extrapolated
    /// material state or first-order fallback is used.
    ///
    /// Temperature and the six faces remain constant; cut calls at environmental
    /// changes. Linear stability does not guarantee positivity: an inadmissible
    /// second stage returns `StageRange` and needs a shorter interval. Returned
    /// states remain directly usable by the existing solid/plate consumers.
    /// The atmosphere owner supplies the corresponding finite-vapor coupling.
    /// This method does not supply a global temporal-error bound. Ordinary
    /// `advance` remains backward Euler.
    pub fn advance_sdirk2(
        &self,
        dt_s: f64,
        faces: &[MoistureFace; 6],
        config: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureCartesianSdirk2Report), MoistureError> {
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
                "Cartesian SDIRK needs two representable forward stage times",
            ));
        }
        let temperature = self.temperature_kelvin();
        let boundaries = [
            faces[0].resolve(temperature)?,
            faces[1].resolve(temperature)?,
            faces[2].resolve(temperature)?,
            faces[3].resolve(temperature)?,
            faces[4].resolve(temperature)?,
            faces[5].resolve(temperature)?,
        ];
        let (first, first_stage) = self.advance(implicit_dt, faces, config, nonlinear, cx)?;
        let known = self.sdirk_load(dt_s, &first, &first_stage, &boundaries, cx)?;
        let (next, endpoint) = self.advance_stage(
            dt_s,
            implicit_dt,
            boundaries,
            config,
            nonlinear,
            Some(&known),
            cx,
        )?;
        poll(cx)?;
        Ok((
            next,
            MoistureCartesianSdirk2Report {
                endpoint,
                first_stage,
            },
        ))
    }
    /// Reuse actual first-stage face fluxes for fixed or evolving reservoirs.
    /// The original field defines storage and Kirchhoff source intervals.
    pub(in crate::moisture) fn sdirk_load(
        &self,
        dt_s: f64,
        first: &Self,
        first_stage: &MoistureCartesianReport,
        boundaries: &[ResolvedMoistureFace; 6],
        cx: &Cx<'_>,
    ) -> Result<StageLoad<6>, MoistureError> {
        let reservoir = std::array::from_fn(|k| boundaries[k].equilibrium);
        let profiles = if self.kirchhoff_diffusivity {
            Some(self.kirchhoff_profiles(reservoir, None, cx)?)
        } else {
            None
        };
        let system = self.system(
            (1.0 - GAMMA) * dt_s,
            boundaries,
            &first_stage.transport_end,
            profiles.as_ref(),
            cx,
        )?;
        let mut increment = vec![0.0; self.moisture.len()];
        for (i, coefficients) in system.operator.forward.iter().enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            for (axis, &coefficient) in coefficients.iter().enumerate() {
                if coefficient != 0.0 {
                    let j = i + system.operator.strides[axis];
                    let flux = coefficient * (first.moisture[j] - first.moisture[i]);
                    increment[i] += flux;
                    increment[j] -= flux;
                }
            }
        }
        let mut lower = self.moisture.iter().copied().fold(f64::INFINITY, f64::min);
        let mut upper = self.moisture.iter().copied().fold(0.0, f64::max);
        for face in &system.faces {
            poll(cx)?;
            if face.coefficient != 0.0 {
                let equilibrium = boundaries[face.side].equilibrium.expect("active boundary");
                increment[face.cell] +=
                    face.coefficient * (equilibrium - first.moisture[face.cell]);
                // Zero surface transfer and impermeable principal directions
                // cannot widen the physical range through a dormant reservoir.
                lower = lower.min(equilibrium);
                upper = upper.max(equilibrium);
            }
        }
        if increment.iter().any(|v| !v.is_finite()) {
            return Err(MoistureError::Invalid(
                "Cartesian SDIRK face increment overflow",
            ));
        }
        Ok(StageLoad::new(
            increment,
            self.transfers(&first.moisture, boundaries, &system),
            [lower, upper],
        ))
    }
}
