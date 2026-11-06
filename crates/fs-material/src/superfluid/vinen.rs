//! Homogeneous Vinen tangle evolution with an explicit effective line energy.
//!
//! With `W=|v_n-v_s|`, line length per volume `L [m^-2]`, and frozen material,
//! the authored closure is
//! `dL/dt = a W L^(3/2) - b L^2`, `dw/dt = -c L w`, `c = B kappa / 3`.
//! The production coefficient `a` is dimensionless and `b` includes circulation
//! (`b=beta_v*kappa` in the usual Vinen convention). This is a local unpolarized
//! tangle model, not vortex transport, nucleation, reconnections or a critical-
//! velocity prediction. A zero seed stays zero. Coefficients are not inferred
//! from a helium name, temperature, or channel geometry.
//!
//! Vinen-law convention: Sciacca, Mongiovi and Jou, *Alternative Vinen's equation
//! and its extension to rotating counterflow superfluid turbulence*, equation
//! (1.3), <https://cab.unime.it/mus/501/1/C1S0801017.pdf>. Coupling a changing
//! line population requires an energy account; see Sciacca, *Non-equilibrium
//! thermodynamics analysis of rotating counterflow superfluid turbulence*,
//! <https://arxiv.org/abs/1008.4390>. The restricted closure HERE is specified
//! by the following derivation, not claimed to reproduce that paper's model.
//!
//! The caller supplies a constant effective energy per length `e [J/m]`. Total
//! represented storage is `U=mu W^2/2 + e L`, `mu=rho_n*rho_s/rho`. Its heat
//! production is `-dU/dt=L(mu*c*W^2-e*a*W*sqrt(L)+e*b*L)`. We admit only
//! `e*a^2 <= 4*mu*c*b`, making this quadratic nonnegative for every state.
//! Growth draws from relative kinetic storage; decay releases line energy.
//! No logarithmic core/spacing dependence or universal vortex EOS is implied.
//!
//! The exposure `tau=integral(L dt)` makes the coupled dynamics linear:
//! `dw/dtau=-c*w`, `d(sqrt(L))/dtau=(a*W-b*sqrt(L))/2`. Its exact positive
//! exponential map preserves the momentum and dissipation properties for ANY
//! nonnegative exposure. A bounded exponential-midpoint march approximates the
//! physical clock to second order. Decay without production uses the exact
//! physical-time solution. There are no allocations or nonlinear iterations.

use super::{
    TwoFluidState, TwoFluidVelocity, finite_vector, invalid, nonnegative, positive, speed,
};
use crate::MaterialError;
use fs_math::{decay, det};

/// Authored coefficients at the same frozen reference as the two-fluid owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VinenParameters {
    /// Positive circulation magnitude [m²/s].
    pub circulation_quantum_m2_s: f64,
    /// Nonnegative dissipative Hall--Vinen B, in the parent's convention.
    pub dissipative_b: f64,
    /// Nonnegative dimensionless coefficient of `|w| L^(3/2)`.
    pub production_coefficient: f64,
    /// Positive coefficient of `L²` [m²/s], INCLUDING circulation.
    pub decay_coefficient_m2_s: f64,
    /// Positive, constant effective vortex energy per unit line length [J/m].
    pub line_energy_j_m: f64,
}

/// Explicit accuracy and bounded-work controls for the source's physical clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VinenStepControls {
    /// Maximum substep times any local production, decay or drag rate, in
    /// `(0, 0.25]`. Refining this controls time error; it is not an error bound.
    pub maximum_fractional_step: f64,
    /// Maximum source substeps per call, in `1..=4096`; exhaustion refuses.
    pub maximum_substeps: usize,
}

/// Immutable admitted local closure, including its original material reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VinenTangle {
    material: TwoFluidState,
    parameters: VinenParameters,
    controls: VinenStepControls,
    drag_per_line_m2_s: f64,
}

/// One private proposed source endpoint. The heat includes line decay and
/// excludes energy retained by line growth; it must be transferred only once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VinenTangleStep {
    pub velocity: TwoFluidVelocity,
    pub line_density_per_m2: f64,
    pub normal_impulse_kg_m2_s: [f64; 3],
    pub superfluid_impulse_kg_m2_s: [f64; 3],
    /// Relative kinetic storage lost to BOTH heat and line growth [J/m³].
    pub counterflow_energy_loss_j_m3: f64,
    pub vortex_energy_before_j_m3: f64,
    pub vortex_energy_after_j_m3: f64,
    pub dissipated_energy_j_m3: f64,
    /// Roundoff in `heat + (line_after-line_before) - kinetic_loss` [J/m³].
    pub energy_balance_defect_j_m3: f64,
    /// Integrated line density `integral(L dt)` [s/m²], with the declared
    /// physical-clock accuracy. This same exposure advances flow and density.
    pub integrated_line_density_s_m2: f64,
    /// Physical-clock substeps; zero for identity and one for exact pure decay.
    pub substeps: usize,
}

impl VinenTangle {
    /// Admit positive line energy/decay and a globally dissipative coefficient
    /// set. `e*a² <= 4*mu*c*b` is checked without squaring the dimensional
    /// coefficients. No calibration or parameter clipping is performed.
    pub fn try_new(
        material: TwoFluidState,
        parameters: VinenParameters,
        controls: VinenStepControls,
    ) -> Result<Self, MaterialError> {
        let p = parameters;
        if ![
            p.circulation_quantum_m2_s,
            p.decay_coefficient_m2_s,
            p.line_energy_j_m,
            controls.maximum_fractional_step,
        ]
        .into_iter()
        .all(positive)
            || ![p.dissipative_b, p.production_coefficient]
                .into_iter()
                .all(nonnegative)
            || controls.maximum_fractional_step > 0.25
            || !(1..=4096).contains(&controls.maximum_substeps)
        {
            return Err(invalid(
                "Vinen coefficients need positive circulation, decay and line energy; nonnegative production/B; fraction in (0,0.25] and 1..=4096 substeps",
            ));
        }
        let c = (p.dissipative_b / 3.0) * p.circulation_quantum_m2_s;
        if !nonnegative(c) || (p.dissipative_b > 0.0 && c == 0.0) {
            return Err(invalid("Vinen drag coefficient is unrepresentable"));
        }
        if p.production_coefficient > 0.0 {
            let permitted = 2.0
                * (material.reduced_density_kg_m3() / p.line_energy_j_m).sqrt()
                * c.sqrt()
                * p.decay_coefficient_m2_s.sqrt();
            if !positive(permitted) || p.production_coefficient > permitted {
                return Err(invalid(
                    "Vinen coefficients violate e*a^2 <= 4*mu*c*b for nonnegative heat production",
                ));
            }
        }
        Ok(Self {
            material,
            parameters,
            controls,
            drag_per_line_m2_s: c,
        })
    }

    pub const fn material(&self) -> &TwoFluidState {
        &self.material
    }
    pub const fn parameters(&self) -> VinenParameters {
        self.parameters
    }
    pub const fn controls(&self) -> VinenStepControls {
        self.controls
    }

    /// Effective stored vortex energy per total fluid volume [J/m³].
    pub fn line_energy_density_j_m3(&self, line_density_per_m2: f64) -> Result<f64, MaterialError> {
        let energy = self.parameters.line_energy_j_m * line_density_per_m2;
        if !nonnegative(line_density_per_m2)
            || !nonnegative(energy)
            || (line_density_per_m2 > 0.0 && energy == 0.0)
        {
            return Err(invalid(
                "Vinen line density or represented line energy is unrepresentable",
            ));
        }
        Ok(energy)
    }

    /// Advance a homogeneous seeded population and its relative flow. Any
    /// refusal leaves the caller's state untouched. Barycentric motion and
    /// relative direction are unchanged. Source work is capped by the authored
    /// substep budget, including when used inside a spatial cancellation tile.
    pub fn advance(
        &self,
        velocity: TwoFluidVelocity,
        line_density_per_m2: f64,
        dt_s: f64,
    ) -> Result<VinenTangleStep, MaterialError> {
        let relative = self.material.counterflow_velocity_m_s(velocity)?;
        let kinetic = self.material.kinetic_energy_density(velocity)?;
        let before = self.line_energy_density_j_m3(line_density_per_m2)?;
        if !nonnegative(dt_s) {
            return Err(invalid(
                "Vinen source duration must be finite and nonnegative",
            ));
        }
        if dt_s == 0.0 || line_density_per_m2 == 0.0 {
            return Ok(VinenTangleStep {
                velocity,
                line_density_per_m2,
                normal_impulse_kg_m2_s: [0.0; 3],
                superfluid_impulse_kg_m2_s: [0.0; 3],
                counterflow_energy_loss_j_m3: 0.0,
                vortex_energy_before_j_m3: before,
                vortex_energy_after_j_m3: before,
                dissipated_energy_j_m3: 0.0,
                energy_balance_defect_j_m3: 0.0,
                substeps: 0,
                integrated_line_density_s_m2: 0.0,
            });
        }
        let w = speed(relative);
        let root = line_density_per_m2.sqrt();
        let (exposure, substeps) = self.exposure(w, root, dt_s)?;
        let (next_root, root_change) = self.root_at_exposure(w, root, exposure)?;
        let next_density = next_root * next_root;
        let after = self.line_energy_density_j_m3(next_density)?;
        if next_density == 0.0 {
            return Err(invalid(
                "positive Vinen density underflowed during a finite source step",
            ));
        }
        let drag = self.drag_per_line_m2_s * exposure;
        if !nonnegative(drag) {
            return Err(invalid("Vinen drag exposure is unrepresentable"));
        }
        let velocity_loss = -det::expm1(-drag);
        let kinetic_loss = kinetic.counterflow_j_m3 * -det::expm1(-2.0 * drag);
        // Keep a weak line increment when next_density rounds to the old value.
        let line_change = self.parameters.line_energy_j_m * root_change * (next_root + root);
        let raw_heat = kinetic_loss - line_change;
        let roundoff = 128.0 * f64::EPSILON * (kinetic_loss + line_change.abs());
        if ![kinetic_loss, line_change, raw_heat, roundoff]
            .into_iter()
            .all(f64::is_finite)
            || raw_heat < -roundoff
        {
            return Err(invalid(
                "Vinen source lost its nonnegative heat or finite energy account",
            ));
        }
        // A negative value within this named floating-point bound is roundoff,
        // not an empirical correction. The returned residual retains the effect.
        let heat = raw_heat.max(0.0);
        let impulse = finite_vector(
            relative.map(|v| -self.material.reduced_density_kg_m3() * velocity_loss * v),
        )?;
        let next = if velocity_loss <= 0.5 {
            TwoFluidVelocity {
                normal_m_s: finite_vector(std::array::from_fn(|i| {
                    velocity.normal_m_s[i]
                        - self.material.superfluid_fraction() * velocity_loss * relative[i]
                }))?,
                superfluid_m_s: finite_vector(std::array::from_fn(|i| {
                    velocity.superfluid_m_s[i]
                        + self.material.normal_fraction() * velocity_loss * relative[i]
                }))?,
            }
        } else {
            self.material.velocities_from_barycentric(
                self.material.barycentric_velocity_m_s(velocity)?,
                relative.map(|v| decay::attenuate(v, drag)),
            )?
        };
        self.material.kinetic_energy_density(next)?;
        Ok(VinenTangleStep {
            velocity: next,
            line_density_per_m2: next_density,
            normal_impulse_kg_m2_s: impulse,
            superfluid_impulse_kg_m2_s: impulse.map(|v| -v),
            counterflow_energy_loss_j_m3: kinetic_loss,
            vortex_energy_before_j_m3: before,
            vortex_energy_after_j_m3: after,
            dissipated_energy_j_m3: heat,
            energy_balance_defect_j_m3: heat + (after - before) - kinetic_loss,
            integrated_line_density_s_m2: exposure,
            substeps,
        })
    }

    fn exposure(&self, w: f64, root: f64, duration: f64) -> Result<(f64, usize), MaterialError> {
        let b = self.parameters.decay_coefficient_m2_s;
        if self.parameters.production_coefficient == 0.0 || w == 0.0 {
            // L=L0/(1+b L0 t); integrate the drag using the SAME decaying L.
            let unattenuated_exposure = root * root * duration;
            let x = if unattenuated_exposure.is_finite() {
                b * unattenuated_exposure
            } else {
                b * root * root * duration
            };
            if !nonnegative(x) {
                return Err(invalid(
                    "Vinen pure-decay duration product is unrepresentable",
                ));
            }
            // Do not divide an underflowed b*L*t by b: L*t and the ensuing
            // friction impulse may still be representable when b is tiny.
            let exposure = if x == 0.0 {
                unattenuated_exposure
            } else if x < 1.0 {
                unattenuated_exposure * (log1p_nonnegative(x) / x)
            } else {
                log1p_nonnegative(x) / b
            };
            if !positive(exposure) {
                return Err(invalid("Vinen integrated line density is unrepresentable"));
            }
            return Ok((exposure, 1));
        }
        let a = self.parameters.production_coefficient;
        let mut elapsed = 0.0;
        let mut exposure = 0.0;
        let mut exposure_correction = 0.0;
        for substep in 1..=self.controls.maximum_substeps {
            let (current_root, _) = self.root_at_exposure(w, root, exposure)?;
            let current_w = decay::attenuate(w, self.drag_per_line_m2_s * exposure);
            let density = current_root * current_root;
            let rate = (self.drag_per_line_m2_s * density)
                .max(0.5 * b * density)
                .max(0.5 * a * current_w * current_root);
            if !positive(rate) || !positive(density) {
                return Err(invalid(
                    "Vinen source rate or positive density is unrepresentable",
                ));
            }
            let remaining = duration - elapsed;
            let h = remaining.min(self.controls.maximum_fractional_step / rate);
            let next_elapsed = if h == remaining {
                duration
            } else {
                elapsed + h
            };
            if !positive(h) || next_elapsed <= elapsed {
                return Err(invalid(
                    "Vinen source substep cannot advance the physical clock",
                ));
            }
            let prediction = exposure + 0.5 * h * density;
            let (mid_root, _) = self.root_at_exposure(w, root, prediction)?;
            let increment = h * mid_root * mid_root;
            let corrected = increment - exposure_correction;
            let next_exposure = exposure + corrected;
            exposure_correction = (next_exposure - exposure) - corrected;
            if !positive(increment) || !nonnegative(next_exposure) || next_exposure <= exposure {
                return Err(invalid("Vinen integrated line density cannot advance"));
            }
            exposure = next_exposure;
            elapsed = next_elapsed;
            if elapsed == duration {
                return Ok((exposure, substep));
            }
        }
        Err(invalid(
            "Vinen source exhausted its maximum_substeps; reduce duration or increase the explicit budget",
        ))
    }

    // Exact positive triangular exponential. The divided exponential handles
    // c=b/2 continuously, without subtracting almost equal exponentials.
    fn root_at_exposure(
        &self,
        w: f64,
        root: f64,
        exposure: f64,
    ) -> Result<(f64, f64), MaterialError> {
        let u = self.drag_per_line_m2_s * exposure;
        let v = 0.5 * self.parameters.decay_coefficient_m2_s * exposure;
        if ![exposure, u, v].into_iter().all(nonnegative) {
            return Err(invalid("Vinen exponential exposure is unrepresentable"));
        }
        let difference = (u - v).abs();
        let divided = if difference == 0.0 {
            1.0
        } else {
            -det::expm1(-difference) / difference
        };
        let production = 0.5 * self.parameters.production_coefficient * w * exposure;
        let produced = decay::attenuate(production, u.min(v)) * divided;
        let change = root * det::expm1(-v) + produced;
        let endpoint = if v <= 0.5 {
            root + change
        } else {
            decay::attenuate(root, v) + produced
        };
        if !positive(endpoint) || !change.is_finite() {
            return Err(invalid("Vinen line-density exponential is unrepresentable"));
        }
        Ok((endpoint, change))
    }
}

// Kahan's correction restores the low bits lost by forming 1+x. This retains
// the workspace's deterministic logarithm, including sub-ulp positive x.
fn log1p_nonnegative(x: f64) -> f64 {
    let sum = 1.0 + x;
    if sum == 1.0 {
        x
    } else {
        det::ln(sum) * (x / (sum - 1.0))
    }
}
