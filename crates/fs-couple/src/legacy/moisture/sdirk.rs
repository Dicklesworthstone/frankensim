//! Second-order stiff integration of the original conservative slab flux.
//!
//! Alexander's two-stage SDIRK method uses gamma = 1 - 1/sqrt(2):
//! M1-M0 = gamma*dt*F(M1),
//! M2-M0 = (1-gamma)*dt*F(M1) + gamma*dt*F(M2).
//! The first rate is evaluated from the actual faces, never reconstructed by
//! subtracting rounded stored moisture. Both nonlinear stages use the original
//! source owner and tridiagonal solver. Reference: R. Alexander, SIAM J. Numer.
//! Anal. 14(6), 1006-1021 (1977); the tableau is also documented at
//! <https://mooseframework.inl.gov/source/timeintegrators/LStableDirk2.html>.
//!
//! L-stability does not imply unconditional monotonicity. Every proposed second
//! stage must stay within the original field/reservoir range and source support;
//! an inadmissible step refuses rather than clipping or reverting to first order.

use super::*;

pub(super) const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

/// One second-order physical interval with both original implicit stage receipts.
#[derive(Clone, Debug)]
pub struct MoistureSdirk2Report {
    /// Complete original-to-final interval and water account. Boundary water is
    /// the RK-weighted sum of both actual stage fluxes. Surface/source and linear
    /// diagnostics describe the final stage, while `transport_start` is original.
    pub endpoint: MoistureStepReport,
    /// Internal stage at t0 + gamma*dt, containing gamma*dt times its flux.
    /// The endpoint uses (1-gamma)*dt times that same first rate plus gamma*dt
    /// times the final rate; do not add this diagnostic stage account again.
    pub first_stage: MoistureStepReport,
}

/// Already-integrated first-stage rate, separate from the physical old field.
pub(super) struct StageLoad<const FACES: usize = 2> {
    increment: Vec<f64>,
    pub(super) boundary_water_kg: [f64; FACES],
    bounds: [f64; 2],
}
impl<const FACES: usize> StageLoad<FACES> {
    pub(super) fn new(
        increment: Vec<f64>,
        boundary_water_kg: [f64; FACES],
        bounds: [f64; 2],
    ) -> Self {
        Self {
            increment,
            boundary_water_kg,
            bounds,
        }
    }

    pub(super) fn add_to(&self, rhs: &mut [f64], cx: &Cx<'_>) -> Result<(), MoistureError> {
        debug_assert_eq!(rhs.len(), self.increment.len());
        for (i, (value, &increment)) in rhs.iter_mut().zip(&self.increment).enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            *value += increment;
            if !value.is_finite() {
                return Err(MoistureError::Invalid("SDIRK moisture load overflow"));
            }
        }
        Ok(())
    }

    /// Source support already admitted by the original field and active
    /// first-stage reservoirs. A later reservoir probe can change its envelope
    /// without removing the physical first-stage load from the second stage.
    pub(super) const fn source_bounds(&self) -> [f64; 2] {
        self.bounds
    }

    pub(super) fn check_range(
        &self,
        candidate: &[f64],
        maximum: f64,
        cx: &Cx<'_>,
    ) -> Result<(), MoistureError> {
        self.check_range_with(candidate, maximum, std::iter::empty(), cx)
    }

    /// The second stage may see a different finite-reservoir activity. Only
    /// equilibria with an actual nonzero current face coefficient extend the
    /// retained first-stage envelope; a scalar search bracket is not a boundary.
    pub(super) fn check_range_with(
        &self,
        candidate: &[f64],
        maximum: f64,
        active_equilibria: impl IntoIterator<Item = f64>,
        cx: &Cx<'_>,
    ) -> Result<(), MoistureError> {
        let [mut lower, mut upper] = self.bounds;
        for (i, equilibrium) in active_equilibria.into_iter().enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            lower = lower.min(equilibrium);
            upper = upper.max(equilibrium);
        }
        let roundoff = 128.0 * f64::EPSILON * upper.max(1.0);
        for (cell, &moisture) in candidate.iter().enumerate() {
            if cell % 128 == 0 {
                poll(cx)?;
            }
            if !moisture.is_finite()
                || moisture < 0.0
                || moisture > maximum
                || moisture < lower - roundoff
                || moisture > upper + roundoff
            {
                return Err(MoistureError::StageRange {
                    cell,
                    moisture,
                    bounds: [lower, upper],
                });
            }
        }
        Ok(())
    }
}

impl MoistureSlab {
    /// Propose a second-order SDIRK interval with endpoint-dependent diffusivity.
    ///
    /// Both stages resolve the same pinned D(M), including the original
    /// Kirchhoff path when selected. Faces and temperature remain constant;
    /// cut intervals at environmental events. The accepted `Self` is consumed
    /// by the existing swelling/solid/plate owners without converting its field.
    /// Calling this method does not change the slab's ordinary advance mode.
    ///
    /// A large step may fail `StageRange` despite linear stability. Reduce dt
    /// rather than clipping water or interpreting a failed stage as physical.
    /// No temporal/global error bound is inferred. The finite-atmosphere owner
    /// supplies coupled RK stages with its actual evolving vapor inventory.
    pub fn advance_sdirk2(
        &self,
        dt_s: f64,
        faces: &[MoistureFace; 2],
        config: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureSdirk2Report), MoistureError> {
        poll(cx)?;
        nonlinear.validate()?;
        let end = self.time_s + dt_s;
        let implicit_dt = GAMMA * dt_s;
        let first_time = self.time_s + implicit_dt;
        if !positive(dt_s)
            || !positive(implicit_dt)
            || !end.is_finite()
            || first_time <= self.time_s
            || first_time >= end
        {
            return Err(MoistureError::Invalid(
                "SDIRK moisture needs two representable forward stage times",
            ));
        }
        let temperature = self.temperature_kelvin();
        let boundaries = [
            faces[0].resolve(temperature)?,
            faces[1].resolve(temperature)?,
        ];
        // Stage configuration is private. The physical old state and the
        // caller's method selection are unchanged by either failed stage.
        let mut staged = self.clone();
        staged.implicit_diffusivity = Some(nonlinear);
        if staged.kirchhoff_diffusivity.is_some() {
            staged.kirchhoff_diffusivity = Some(nonlinear);
        }
        let (first, first_stage) = staged.advance(implicit_dt, faces, config, cx)?;
        let known = staged.sdirk_load(dt_s, &first, &first_stage, &boundaries, cx)?;
        let (mut next, endpoint) =
            staged.advance_stage(dt_s, implicit_dt, boundaries, config, Some(&known), cx)?;
        next.implicit_diffusivity = self.implicit_diffusivity;
        next.kirchhoff_diffusivity = self.kirchhoff_diffusivity;
        poll(cx)?;
        Ok((
            next,
            MoistureSdirk2Report {
                endpoint,
                first_stage,
            },
        ))
    }

    /// Integrate the same first-stage face rates for fixed or evolving reservoirs.
    /// The physical storage reference remains the original accepted slab.
    pub(in crate::moisture) fn sdirk_load(
        &self,
        dt_s: f64,
        first: &Self,
        first_stage: &MoistureStepReport,
        boundaries: &[ResolvedMoistureFace; 2],
        cx: &Cx<'_>,
    ) -> Result<StageLoad, MoistureError> {
        let surface = [boundaries[0].equilibrium, boundaries[1].equilibrium];
        let profile = self
            .kirchhoff_diffusivity
            .map(|_| self.kirchhoff_profile(surface, None, cx))
            .transpose()?;
        let (operator, _, coupling) = self.transport_system(
            (1.0 - GAMMA) * dt_s,
            boundaries,
            &first_stage.transport_end,
            profile.as_ref(),
            cx,
        )?;
        let n = self.moisture.len();
        let mut increment = vec![0.0; n];
        for (i, &coefficient) in operator.off.iter().enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            let flux = coefficient * (first.moisture[i + 1] - first.moisture[i]);
            increment[i] += flux;
            increment[i + 1] -= flux;
        }
        for (side, i) in [0, n - 1].into_iter().enumerate() {
            if let Some(equilibrium) = surface[side] {
                increment[i] += coupling[side].0 * (equilibrium - first.moisture[i]);
            }
        }
        if increment.iter().any(|v| !v.is_finite()) {
            return Err(MoistureError::Invalid(
                "SDIRK moisture face increment overflow",
            ));
        }
        let active = || (0..2).filter(|&side| coupling[side].0 != 0.0)
            .filter_map(|side| surface[side]);
        Ok(StageLoad {
            increment,
            boundary_water_kg: self.boundary_transfer(&first.moisture, boundaries, &coupling),
            bounds: [
                self.moisture.iter().copied().chain(active()).fold(f64::INFINITY, f64::min),
                self.moisture.iter().copied().chain(active()).fold(0.0, f64::max),
            ],
        })
    }
}
