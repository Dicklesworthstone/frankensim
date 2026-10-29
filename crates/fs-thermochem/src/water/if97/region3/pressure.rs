//! Pressure-driven states from the original IF97 Region 3 Helmholtz equation.
//!
//! The density solve uses the original analytic isothermal stiffness and a
//! safeguarded Newton/bisection bracket on the requested stable outer branch.
//! Temperature powers are cached once per solve. No backward fit, extrapolated
//! remote density root, pressure replacement or caloric datum shift is used.
//! Final admission requires both the requested pressure residual and a numerical
//! density bracket width, and outward evaluation proves positive Cv/stiffness
//! across that final bracket. The reported bracket is a floating-point solve
//! diagnostic, not a certified enclosure of a real-arithmetic pressure root.
//! Near-critical loss of resolution therefore refuses explicitly.
//!
//! A saturation pair retains the original Region 4 temperature and both raw
//! Region 3 endpoints. Their pressure and Gibbs mismatch remain visible. The
//! pair defines no bulk mixture, latent-energy correction, or transport law.
//! Attribution: IAPWS R7-97(2012), equations (28), (30), and (31), Tables 30–33,
//! <https://iapws.org/technical-guidance/release/IF97-Rev.download>.

use super::super::{If97Error, saturation_temperature};
use super::{
    CRITICAL_DENSITY_KG_M3, CRITICAL_PRESSURE_PA, CRITICAL_TEMPERATURE_K, DENSITY_CEILING_KG_M3,
    If97DenseError, If97DenseState, LOG_COEFFICIENT, MAXIMUM_PRESSURE_PA, MAXIMUM_TEMPERATURE_K,
    MINIMUM_TEMPERATURE_K, R_SPECIFIC_J_KG_K, TERMS, boundary23_pressure, saturation_pressure,
};
use crate::{MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1};
use core::fmt;
use fs_qty::{Density, Pressure, Temperature};

/// Numerical solve allowances, distinct from the source fit's uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensePressureSolve {
    pub pressure_tolerance_pa: f64,
    pub density_tolerance_kg_m3: f64,
    /// Counts branch-search and root iterations together, in 1..=256.
    pub maximum_iterations: usize,
}
impl DensePressureSolve {
    fn validate(self) -> Result<(), DensePressureError> {
        if !self.pressure_tolerance_pa.is_finite()
            || self.pressure_tolerance_pa < 0.0
            || !self.density_tolerance_kg_m3.is_finite()
            || self.density_tolerance_kg_m3 <= 0.0
            || !(1..=256).contains(&self.maximum_iterations)
        {
            return Err(DensePressureError::Invalid(
                "finite nonnegative pressure allowance, positive density width, and 1..=256 iterations required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum DensePressureError {
    Invalid(&'static str),
    Source(If97DenseError),
    Saturation(If97Error),
    OutsideRegion {
        pressure_pa: f64,
        temperature_k: f64,
    },
    AmbiguousSaturation,
    CriticalSingularity,
    BranchNotBracketed {
        pressure_pa: f64,
        temperature_k: f64,
    },
    NotConverged {
        iterations: usize,
        evaluations: usize,
        pressure_residual_pa: f64,
        density_width_kg_m3: f64,
    },
    UnresolvedStability {
        stiffness_lower_bound_m2_s2: f64,
        cv_lower_bound_j_kg_k: f64,
    },
    Cancelled,
}
impl fmt::Display for DensePressureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IF97 dense pressure inversion: {self:?}")
    }
}
impl std::error::Error for DensePressureError {}
impl From<If97DenseError> for DensePressureError {
    fn from(error: If97DenseError) -> Self {
        Self::Source(error)
    }
}
impl From<If97Error> for DensePressureError {
    fn from(error: If97Error) -> Self {
        Self::Saturation(error)
    }
}

/// Unmodified source endpoint, achieved residual and bounded work diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensePressureSolution {
    pub state: If97DenseState,
    pub requested_pressure_pa: f64,
    pub pressure_residual_pa: f64,
    /// Numerical pressure-sign bracket; not a source-uncertainty interval.
    pub density_bracket_kg_m3: [f64; 2],
    pub iterations: usize,
    /// Original pressure/stiffness queries, excluding the final full state.
    pub evaluations: usize,
}

impl If97DenseState {
    /// Solve the stable original Region 3 state at prescribed pressure and T.
    /// Equality to Region 4 is ambiguous: request a saturation pair instead.
    pub fn at_pressure_temperature(
        pressure: Pressure,
        temperature: Temperature,
        controls: DensePressureSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensePressureSolution, DensePressureError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let (p, t) = (pressure.value(), temperature.value());
        validate_domain(p, t)?;
        let branch = if t < CRITICAL_TEMPERATURE_K {
            let ps = saturation_pressure(temperature)?.value();
            if p == ps {
                return Err(DensePressureError::AmbiguousSaturation);
            }
            if p > ps {
                Branch::Liquid
            } else {
                Branch::Vapor
            }
        } else {
            if t == CRITICAL_TEMPERATURE_K && p == CRITICAL_PRESSURE_PA {
                return Err(DensePressureError::CriticalSingularity);
            }
            Branch::Single
        };
        solve(p, t, branch, controls, &mut cancelled)
    }
}

/// Paired original liquid/vapor roots above the Region 1/2 saturation join.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97DenseSaturationPair {
    requested_pressure_pa: f64,
    temperature_k: f64,
    saturation_pressure_residual_pa: f64,
    liquid: DensePressureSolution,
    vapor: DensePressureSolution,
}
impl If97DenseSaturationPair {
    pub fn at_temperature(
        temperature: Temperature,
        controls: DensePressureSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DensePressureError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let t = temperature.value();
        validate_saturation_temperature(t)?;
        let p = saturation_pressure(temperature)?.value();
        Self::evaluate(p, t, p, controls, &mut cancelled)
    }
    pub fn at_pressure(
        pressure: Pressure,
        controls: DensePressureSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DensePressureError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let p = pressure.value();
        let t = saturation_temperature(pressure)?.value();
        validate_saturation_temperature(t)?;
        let ps = saturation_pressure(Temperature::new(t))?.value();
        Self::evaluate(p, t, ps, controls, &mut cancelled)
    }
    fn evaluate(
        requested_p: f64,
        t: f64,
        saturation_p: f64,
        controls: DensePressureSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DensePressureError> {
        // Each source state must remain on its admitted side of p_sat(T), even
        // when the inverse Region 4 equation does not round-trip exactly.
        let mut liquid = solve(saturation_p, t, Branch::Liquid, controls, cancelled)?;
        let mut vapor = solve(saturation_p, t, Branch::Vapor, controls, cancelled)?;
        for endpoint in [&mut liquid, &mut vapor] {
            endpoint.requested_pressure_pa = requested_p;
            endpoint.pressure_residual_pa = endpoint.state.pressure().value() - requested_p;
            if endpoint.pressure_residual_pa.abs() > controls.pressure_tolerance_pa {
                return Err(DensePressureError::NotConverged {
                    iterations: endpoint.iterations,
                    evaluations: endpoint.evaluations,
                    pressure_residual_pa: endpoint.pressure_residual_pa,
                    density_width_kg_m3: endpoint.density_bracket_kg_m3[1]
                        - endpoint.density_bracket_kg_m3[0],
                });
            }
        }
        if liquid.state.density().value() <= vapor.state.density().value()
            || liquid.state.specific_enthalpy().value() >= vapor.state.specific_enthalpy().value()
            || liquid.state.specific_entropy().value() >= vapor.state.specific_entropy().value()
        {
            return Err(DensePressureError::BranchNotBracketed {
                pressure_pa: requested_p,
                temperature_k: t,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            requested_pressure_pa: requested_p,
            temperature_k: t,
            saturation_pressure_residual_pa: saturation_p - requested_p,
            liquid,
            vapor,
        })
    }
    pub const fn pressure(self) -> Pressure {
        Pressure::new(self.requested_pressure_pa)
    }
    pub const fn temperature(self) -> Temperature {
        Temperature::new(self.temperature_k)
    }
    pub const fn liquid(&self) -> &If97DenseState {
        &self.liquid.state
    }
    pub const fn vapor(&self) -> &If97DenseState {
        &self.vapor.state
    }
    pub const fn liquid_solution(self) -> DensePressureSolution {
        self.liquid
    }
    pub const fn vapor_solution(self) -> DensePressureSolution {
        self.vapor
    }
    pub const fn saturation_pressure_residual_pa(self) -> f64 {
        self.saturation_pressure_residual_pa
    }
    pub fn latent_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(
            self.vapor.state.specific_enthalpy().value()
                - self.liquid.state.specific_enthalpy().value(),
        )
    }
    pub fn latent_entropy(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(
            self.vapor.state.specific_entropy().value()
                - self.liquid.state.specific_entropy().value(),
        )
    }
    pub fn gibbs_energy_gap(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(
            self.vapor.state.specific_gibbs_energy().value()
                - self.liquid.state.specific_gibbs_energy().value(),
        )
    }
}

fn validate_domain(p: f64, t: f64) -> Result<(), DensePressureError> {
    if !p.is_finite() || !t.is_finite() || p <= 0.0 {
        return Err(DensePressureError::Invalid(
            "positive finite pressure and finite temperature required",
        ));
    }
    if !(MINIMUM_TEMPERATURE_K..=MAXIMUM_TEMPERATURE_K).contains(&t)
        || p > MAXIMUM_PRESSURE_PA
        || p < boundary23_pressure(Temperature::new(t))?.value()
    {
        return Err(DensePressureError::OutsideRegion {
            pressure_pa: p,
            temperature_k: t,
        });
    }
    Ok(())
}
fn validate_saturation_temperature(t: f64) -> Result<(), DensePressureError> {
    if !t.is_finite() || t <= MINIMUM_TEMPERATURE_K || t >= CRITICAL_TEMPERATURE_K {
        return Err(DensePressureError::Invalid(
            "dense saturation requires 623.15 < T < 647.096 K",
        ));
    }
    Ok(())
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DensePressureError> {
    if cancelled() {
        Err(DensePressureError::Cancelled)
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Branch {
    Liquid,
    Vapor,
    Single,
}
#[derive(Clone, Copy)]
struct Point {
    rho: f64,
    pressure: f64,
    stiffness: f64,
}
struct Isotherm {
    powers: [f64; 27],
    rt: f64,
    evaluations: usize,
}
impl Isotherm {
    fn new(t: f64) -> Self {
        let mut powers = [1.0; 27];
        let tau = CRITICAL_TEMPERATURE_K / t;
        for j in 1..powers.len() {
            powers[j] = powers[j - 1] * tau;
        }
        Self {
            powers,
            rt: R_SPECIFIC_J_KG_K * t,
            evaluations: 0,
        }
    }
    fn evaluate(&mut self, rho: f64) -> Point {
        self.evaluations += 1;
        let mut powers = [1.0; 12];
        let delta = rho / CRITICAL_DENSITY_KG_M3;
        for i in 1..powers.len() {
            powers[i] = powers[i - 1] * delta;
        }
        // Preserve the original source d/dd operation order, while caching tau.
        let (mut d, mut dd) = (LOG_COEFFICIENT, -LOG_COEFFICIENT);
        for &(i, j, n) in &TERMS {
            let term = n * powers[i] * self.powers[j];
            let di = i as f64;
            d += di * term;
            dd += di * (di - 1.0) * term;
        }
        Point {
            rho,
            pressure: rho * self.rt * d,
            stiffness: self.rt * (2.0 * d + dd),
        }
    }
}

fn solve(
    p: f64,
    t: f64,
    branch: Branch,
    controls: DensePressureSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DensePressureSolution, DensePressureError> {
    validate_domain(p, t)?;
    poll(cancelled)?;
    let mut equation = Isotherm::new(t);
    let mut lo = equation.evaluate(0.0);
    let mut hi = equation.evaluate(DENSITY_CEILING_KG_M3);
    let middle = equation.evaluate(CRITICAL_DENSITY_KG_M3);
    let mut iterations = 0;
    match branch {
        Branch::Liquid => {
            // Move from the known high-density stable side toward the spinodal
            // only until a positive-stiffness point below the target is found.
            let mut unstable = middle;
            if middle.stiffness > 0.0 && middle.pressure < p {
                lo = middle;
            } else {
                while iterations < controls.maximum_iterations {
                    poll(cancelled)?;
                    iterations += 1;
                    let trial = equation.evaluate(unstable.rho + 0.5 * (hi.rho - unstable.rho));
                    if trial.stiffness > 0.0 {
                        if trial.pressure < p {
                            lo = trial;
                            break;
                        }
                        hi = trial;
                    } else {
                        unstable = trial;
                    }
                }
                if lo.rho == 0.0 {
                    return Err(DensePressureError::BranchNotBracketed {
                        pressure_pa: p,
                        temperature_k: t,
                    });
                }
            }
        }
        Branch::Vapor => {
            let mut unstable = middle;
            if middle.stiffness > 0.0 && middle.pressure > p {
                hi = middle;
            } else {
                while iterations < controls.maximum_iterations {
                    poll(cancelled)?;
                    iterations += 1;
                    let trial = equation.evaluate(lo.rho + 0.5 * (unstable.rho - lo.rho));
                    if trial.stiffness > 0.0 {
                        if trial.pressure > p {
                            hi = trial;
                            break;
                        }
                        lo = trial;
                    } else {
                        unstable = trial;
                    }
                }
                if hi.rho == DENSITY_CEILING_KG_M3 {
                    return Err(DensePressureError::BranchNotBracketed {
                        pressure_pa: p,
                        temperature_k: t,
                    });
                }
            }
        }
        Branch::Single => {
            // Avoid crossing the singular critical-density neighborhood when
            // the requested isotherm is exactly Tc and the target is off pc.
            if t == CRITICAL_TEMPERATURE_K {
                if p < middle.pressure {
                    hi = middle;
                } else {
                    lo = middle;
                }
            }
        }
    }
    if lo.pressure > p || hi.pressure < p || !hi.pressure.is_finite() {
        return Err(DensePressureError::BranchNotBracketed {
            pressure_pa: p,
            temperature_k: t,
        });
    }
    let mut current = equation.evaluate(lo.rho + 0.5 * (hi.rho - lo.rho));
    let mut best_residual = f64::INFINITY;
    while iterations < controls.maximum_iterations {
        poll(cancelled)?;
        iterations += 1;
        let residual = current.pressure - p;
        if !residual.is_finite() {
            return Err(DensePressureError::Source(If97DenseError::Unrepresentable(
                "pressure inverse",
            )));
        }
        if residual.abs() < best_residual.abs() {
            best_residual = residual;
        }
        if residual <= 0.0 {
            lo = current;
        } else {
            hi = current;
        }
        if residual.abs() <= controls.pressure_tolerance_pa {
            if let Ok(state) = If97DenseState::at_density_temperature(
                Density::new(current.rho),
                Temperature::new(t),
            ) {
                // Once Newton has found the pressure, test a narrow numerical
                // bracket directly instead of spending many one-sided Newton
                // steps shrinking a distant opposite endpoint.
                let radius = 0.25 * controls.density_tolerance_kg_m3;
                let left = (current.rho - radius).max(lo.rho);
                let right = (current.rho + radius).min(hi.rho);
                if left < right && right - left <= controls.density_tolerance_kg_m3 {
                    let a = equation.evaluate(left);
                    let b = equation.evaluate(right);
                    if a.pressure <= p && b.pressure >= p {
                        lo = a;
                        hi = b;
                    }
                }
                if hi.rho - lo.rho <= controls.density_tolerance_kg_m3 {
                    let (stiffness, cv) = stability_bounds(lo.rho, hi.rho, t);
                    if stiffness <= 0.0 || cv <= 0.0 || !stiffness.is_finite() || !cv.is_finite() {
                        return Err(DensePressureError::UnresolvedStability {
                            stiffness_lower_bound_m2_s2: stiffness,
                            cv_lower_bound_j_kg_k: cv,
                        });
                    }
                    poll(cancelled)?;
                    return Ok(DensePressureSolution {
                        state,
                        requested_pressure_pa: p,
                        pressure_residual_pa: state.pressure().value() - p,
                        density_bracket_kg_m3: [lo.rho, hi.rho],
                        iterations,
                        evaluations: equation.evaluations,
                    });
                }
            }
        }
        let newton = current.rho - residual / current.stiffness;
        let width = hi.rho - lo.rho;
        let candidate = if current.stiffness > 0.0
            && newton > lo.rho + 0.05 * width
            && newton < hi.rho - 0.05 * width
        {
            newton
        } else {
            lo.rho + 0.5 * width
        };
        if candidate == lo.rho || candidate == hi.rho {
            break;
        }
        current = equation.evaluate(candidate);
    }
    Err(DensePressureError::NotConverged {
        iterations,
        evaluations: equation.evaluations,
        pressure_residual_pa: best_residual,
        density_width_kg_m3: hi.rho - lo.rho,
    })
}

// A small outward polynomial enclosure is needed only for the final density
// bracket. It bounds the retained coefficients, not experimental uncertainty.
#[derive(Clone, Copy)]
struct Interval {
    lo: f64,
    hi: f64,
}
impl Interval {
    const fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn add(self, b: Self) -> Self {
        Self {
            lo: (self.lo + b.lo).next_down(),
            hi: (self.hi + b.hi).next_up(),
        }
    }
    fn mul(self, b: Self) -> Self {
        let p = [
            self.lo * b.lo,
            self.lo * b.hi,
            self.hi * b.lo,
            self.hi * b.hi,
        ];
        Self {
            lo: p.into_iter().fold(f64::INFINITY, f64::min).next_down(),
            hi: p.into_iter().fold(f64::NEG_INFINITY, f64::max).next_up(),
        }
    }
    fn scale(self, x: f64) -> Self {
        self.mul(Self::point(x))
    }
}
fn stability_bounds(lo: f64, hi: f64, t: f64) -> (f64, f64) {
    let delta = Interval {
        lo: (lo / CRITICAL_DENSITY_KG_M3).next_down(),
        hi: (hi / CRITICAL_DENSITY_KG_M3).next_up(),
    };
    let tau = Interval {
        lo: (CRITICAL_TEMPERATURE_K / t).next_down(),
        hi: (CRITICAL_TEMPERATURE_K / t).next_up(),
    };
    let mut density = [Interval::point(1.0); 12];
    let mut temperature = [Interval::point(1.0); 27];
    for i in 1..density.len() {
        density[i] = density[i - 1].mul(delta);
    }
    for j in 1..temperature.len() {
        temperature[j] = temperature[j - 1].mul(tau);
    }
    let mut stiffness = Interval::point(LOG_COEFFICIENT);
    let mut cv = Interval::point(0.0);
    for &(i, j, n) in &TERMS {
        let term = Interval::point(n).mul(density[i]).mul(temperature[j]);
        stiffness = stiffness.add(term.scale((i * (i + 1)) as f64));
        cv = cv.add(term.scale(-(j as f64) * (j as f64 - 1.0)));
    }
    (
        stiffness.scale(R_SPECIFIC_J_KG_K).scale(t).lo,
        cv.scale(R_SPECIFIC_J_KG_K).lo,
    )
}

#[cfg(test)]
mod tests;
