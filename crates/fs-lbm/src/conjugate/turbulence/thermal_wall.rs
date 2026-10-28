//! Equilibrium thermal wall resistance on the existing Spalding momentum law.
//!
//! For constant properties and constant wall-normal heat flux,
//! `dT+/dy+ = Pr / (1 + Pr/Pr_t * nu_t/nu)` and `nu_t/nu = dy+/du+ - 1`.
//! Integration in u+ resolves the viscous and buffer layers without imposing
//! a discontinuous linear/log switch. This is a low-Mach equilibrium model:
//! no separation, roughness, radiation, property variation or viscous heating
//! is inferred. A wall model is not a resolved-boundary-layer certificate.

use super::{KAPPA, SPALDING_E};
use crate::conjugate::{ChtError, domain::FluidProperties, poll};
use fs_exec::CancelGate;

mod coupling;
pub use coupling::{
    PreparedWallHeat, WallHeatFace, WallHeatFaceReport, WallHeatSolution,
    prepare_wall_heat, solve_energy_with_wall_heat,
};

/// Explicit quadrature and gradient-diffusion controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallHeatConfig {
    /// Positive turbulent Prandtl number; independent of molecular Pr.
    pub turbulent_prandtl: f64,
    /// Relative change required between successive composite Simpson rules.
    /// This is a convergence diagnostic, not a rigorous error enclosure.
    pub relative_tolerance: f64,
    /// Maximum mesh doublings (3..=20); at most 2^n + 1 integrand evaluations.
    pub max_refinements: u32,
}

impl WallHeatConfig {
    /// Admit the model and finite work budget before evaluating it.
    pub fn validate(self) -> Result<(), ChtError> {
        if !self.turbulent_prandtl.is_finite() || self.turbulent_prandtl <= 0.0
            || !self.relative_tolerance.is_finite()
            || !(1e-12..=1e-3).contains(&self.relative_tolerance)
            || !(3..=20).contains(&self.max_refinements)
        {
            return Err(invalid("positive Pr_t, tolerance in [1e-12,1e-3], and 3..=20 refinements required"));
        }
        Ok(())
    }
}

/// One admitted constant-property fluid-side wall law, independent of area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallHeatTransfer {
    /// Wall friction velocity inferred from the tangential speed, m/s.
    pub friction_velocity_m_s: f64,
    /// Tangential velocity in wall units (zero at rest).
    pub u_plus: f64,
    /// Distance in wall units (zero at rest).
    pub y_plus: f64,
    /// Temperature in wall units (zero at rest; dimensional resistance remains finite).
    pub temperature_plus: f64,
    /// Area-specific fluid-side thermal resistance, m^2 K/W.
    pub resistance_m2_k_w: f64,
    /// Reciprocal resistance, W/(m^2 K).
    pub heat_transfer_w_m2_k: f64,
    /// Successive-rule relative change; zero for analytic limits.
    pub quadrature_relative_change: f64,
    /// Integrand evaluations actually used; zero for analytic limits.
    pub quadrature_evaluations: usize,
}

pub(super) fn invalid(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput { field: "thermal.wall_law", reason: reason.into() }
}

/// Positive exponential remainder avoids cancellation in the viscous layer.
fn remainder(x: f64, order: u32) -> f64 {
    if x < 0.5 {
        let mut term = 1.0;
        for n in 1..=order { term *= x / f64::from(n); }
        let mut sum = term;
        for n in order + 1..=32 {
            term *= x / f64::from(n);
            sum += term;
        }
        sum
    } else {
        let r3 = x.exp() - 1.0 - x - 0.5 * x * x;
        if order == 3 { r3 } else { r3 - x * x * x / 6.0 }
    }
}

/// (y+, dy+/du+) on exactly the LVEL Spalding profile.
fn profile(u: f64) -> (f64, f64) {
    let x = KAPPA * u;
    (u + remainder(x, 4) / SPALDING_E,
        1.0 + KAPPA * remainder(x, 3) / SPALDING_E)
}

fn wall_coordinates(re: f64, gate: &CancelGate) -> Result<(f64, f64), ChtError> {
    // y+(u) >= u: sqrt(Re) is an upper bound. Grow a smaller bound
    // geometrically so very large Re does not waste bisections in overflow.
    let mut high = re.sqrt().min(1.0);
    for _ in 0..16 {
        if high * profile(high).0 >= re { break; }
        high *= 2.0;
    }
    let mut low = 0.0;
    for iteration in 0..96 {
        if iteration % 16 == 0 { poll(gate)?; }
        let mid = low + 0.5 * (high - low);
        if mid == low || mid == high { break; }
        if mid * profile(mid).0 > re { high = mid; } else { low = mid; }
    }
    let u = low + 0.5 * (high - low);
    let y = profile(u).0;
    if !u.is_finite() || !y.is_finite() || u <= 0.0 || y <= 0.0
        || (u * (y / re) - 1.0).abs() > 1e-11
    {
        return Err(invalid("Spalding wall coordinates are not representable"));
    }
    Ok((u, y))
}

fn integral(u: f64, ratio: f64, config: WallHeatConfig, gate: &CancelGate)
    -> Result<(f64, f64, usize), ChtError>
{
    // Integrate T+/Pr. This form avoids overflowing ratio*(slope-1).
    let value = |x| {
        let inverse = profile(x).1.recip();
        (inverse + ratio * (1.0 - inverse)).recip()
    };
    let mut trapezoid = 0.5 * u * (1.0 + value(u));
    let mut previous = f64::NAN;
    for level in 1..=config.max_refinements {
        poll(gate)?;
        let panels = 1usize << level;
        let h = u / panels as f64;
        let mut sum = 0.0;
        for j in 0..panels / 2 {
            if j % 64 == 0 { poll(gate)?; }
            sum += value((2 * j + 1) as f64 * h);
        }
        let next = 0.5 * trapezoid + h * sum;
        let simpson = next + (next - trapezoid) / 3.0;
        if !(simpson.is_finite() && simpson > 0.0) {
            return Err(invalid("thermal wall integral is not representable"));
        }
        let change = ((simpson - previous) / simpson).abs();
        if level >= 3 && change <= config.relative_tolerance {
            return Ok((simpson, change, panels + 1));
        }
        previous = simpson;
        trapezoid = next;
    }
    Err(invalid("thermal wall quadrature exhausted its refinement budget"))
}

/// Infer friction from the tangential wall-relative speed and integrate the
/// thermal resistance from the stationary wall to the sample point.
///
/// `distance_m` is the actual wall-normal distance, not a hydraulic diameter.
/// Only tangential velocity belongs in `tangential_speed_m_s`; normal velocity
/// is not a source of wall shear. At rest the result is exactly `distance/k`.
/// At `Pr == Pr_t`, `T+ = Pr*u+` gives the exact Reynolds-analogy limit.
///
/// # Errors
/// Nonfinite/nonpositive properties, negative speed, unrepresentable derived
/// values, exhausted quadrature, or cancellation. No silent molecular fallback.
pub fn equilibrium_wall_heat_transfer(
    tangential_speed_m_s: f64,
    distance_m: f64,
    fluid: &FluidProperties,
    config: WallHeatConfig,
    gate: &CancelGate,
) -> Result<WallHeatTransfer, ChtError> {
    poll(gate)?;
    config.validate()?;
    fluid.validate()?;
    if !tangential_speed_m_s.is_finite() || tangential_speed_m_s < 0.0
        || !distance_m.is_finite() || distance_m <= 0.0
    {
        return Err(invalid("finite nonnegative tangential speed and positive wall distance required"));
    }
    let molecular = distance_m / fluid.conductivity_w_m_k;
    let pr = fluid.prandtl();
    let ratio = pr / config.turbulent_prandtl;
    let re = tangential_speed_m_s * distance_m / fluid.kinematic_viscosity_m2_s;
    if ![molecular, pr, ratio].into_iter().all(|v| v.is_finite() && v > 0.0)
        || !re.is_finite() || (tangential_speed_m_s > 0.0 && re == 0.0)
    {
        return Err(invalid("derived molecular transport or Reynolds number is not representable"));
    }
    let (u, y, friction, temperature, resistance, change, evaluations) = if re == 0.0 {
        (0.0, 0.0, 0.0, 0.0, molecular, 0.0, 0)
    } else {
        let (u, y) = wall_coordinates(re, gate)?;
        let (thermal, change, count) = if ratio == 1.0 {
            (u, 0.0, 0)
        } else {
            integral(u, ratio, config, gate)?
        };
        (u, y, tangential_speed_m_s / u, pr * thermal,
            molecular * (thermal / y), change, count)
    };
    let heat_transfer = resistance.recip();
    if ![friction, temperature, resistance, heat_transfer].into_iter().all(f64::is_finite)
        || resistance <= 0.0 || heat_transfer <= 0.0
    {
        return Err(invalid("dimensional wall heat transport is not representable"));
    }
    poll(gate)?;
    Ok(WallHeatTransfer {
        friction_velocity_m_s: friction, u_plus: u, y_plus: y,
        temperature_plus: temperature, resistance_m2_k_w: resistance,
        heat_transfer_w_m2_k: heat_transfer, quadrature_relative_change: change,
        quadrature_evaluations: evaluations,
    })
}

#[cfg(test)]
mod tests;
