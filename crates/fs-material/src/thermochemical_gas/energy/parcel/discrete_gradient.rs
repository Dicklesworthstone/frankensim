//! Frozen ideal-gas heat/work steps with a discrete Gibbs identity.
//!
//! Write d = integral(Cv/T dT), e = integral(Cv dT), l = ln(V1/V0).
//! The original caloric source owns both integrals. Set theta=e/d (theta=T0
//! at equal temperatures) and p_bar=m R theta l/(V1-V0). Then
//! Delta U = theta Delta S - p_bar Delta V, with Delta S=m(d+R l).
//! Solving Delta U+p_bar Delta V=Q therefore gives Delta S=Q/theta:
//! adiabatic motion has no algorithmic entropy drift beyond the solve budget.
//!
//! This symmetric discrete-gradient law supplies a pressure conjugate to the
//! ACTUAL volume change and a temperature conjugate to entropy transfer. It is
//! not the midpoint pressure, a heat-transfer law, or a claim of physical
//! reversibility for arbitrary supplied heat. A thermal consumer must use
//! theta in its exchange law and include the other body's entropy. Mechanics
//! must use p_bar in both impulse and boundary work.
//!
//! Only the pressure sensitivities use quadrature: near equal temperatures,
//! four-point Gaussian integration of the original Cv/T moment avoids a
//! subtractive zero-over-zero derivative. Source-region seams are split. The
//! energy/entropy residuals always use the original integrated source laws.
//!
//! Derivation: MIT 16.050, notes 1.B.3--1.B.4 (Gibbs equation and ideal-gas
//! entropy); discrete-gradient framework: Anahory Simoes, Martin de Diego and
//! Maschke, https://arxiv.org/abs/2303.08034. This is the explicit frozen,
//! uniform ideal-gas specialization above, not that paper's general theorem.

use super::ClosedGasParcel;
use super::super::{
    GasEnergyError, GasEnergySolve, ResolvedThermochemicalGasState, ThermochemicalGasModel,
    map_root, poll, positive, root, updated_point,
};
use crate::thermochemical_gas::failure;
use fs_qty::Temperature;

/// Actual source-law balances and discrete pressure sensitivities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasDiscreteGradientStep {
    pub initial_volume_m3: f64,
    pub final_volume_m3: f64,
    /// Whole-step heat into the gas [J].
    pub supplied_heat_j: f64,
    /// Discrete pressure effort paired with the ACTUAL volume change [Pa].
    /// Use this for work and momentum. For a nonideal gas it can differ from
    /// the EOS pressure of the source-bearing sampling state.
    pub work_pressure_pa: f64,
    /// Work out of the gas, work_pressure_pa*(V1-V0) [J].
    pub work_by_gas_j: f64,
    /// Original-law m Delta u, independent of formation references [J].
    pub stored_energy_change_j: f64,
    /// Delta U + work - supplied heat [J].
    pub energy_residual_j: f64,
    /// Original-law entropy change at the two conserved states [J/K].
    pub entropy_change_j_k: f64,
    /// Supplied heat divided by the returned thermodynamic mean temperature.
    /// This is the gas-side entropy transfer, not the external bath's entropy.
    pub entropy_transfer_j_k: f64,
    /// Actual gas entropy change minus gas-side entropy transfer [J/K].
    pub entropy_residual_j_k: f64,
    /// Delta U + work - theta Delta S [J].
    pub gibbs_residual_j: f64,
    /// d(p_bar)/d(V1), holding the old state and supplied heat fixed [Pa/m3].
    /// Near-coincident source moments use the implementing owner's bounded quadrature.
    pub pressure_final_volume_derivative_pa_m3: f64,
    /// d(p_bar)/d(Q), holding both volumes and the old state fixed [Pa/J].
    pub pressure_heat_derivative_pa_j: f64,
    /// Interior iterations in the single bounded temperature inversion.
    pub iterations: usize,
}

/// An immutable trial; publish only after the complete coupled step accepts.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDiscreteGradientCandidate<P = ClosedGasParcel, S = ResolvedThermochemicalGasState> {
    pub(crate) endpoint: P,
    pub(crate) thermodynamic_mean: S,
    pub(crate) step: GasDiscreteGradientStep,
}
impl<P, S> GasDiscreteGradientCandidate<P, S> {
    #[must_use]
    pub const fn endpoint(&self) -> &P { &self.endpoint }
    /// Original material resolved at the entropy-conjugate temperature.
    /// Ideal gas uses logarithmic-mean volume; the PR implementation uses
    /// arithmetic mean specific volume. The latter state's EOS pressure need
    /// not be the discrete work effort: always use step().work_pressure_pa.
    #[must_use]
    pub const fn thermodynamic_mean(&self) -> &S {
        &self.thermodynamic_mean
    }
    #[must_use]
    pub const fn step(&self) -> &GasDiscreteGradientStep { &self.step }
    #[must_use]
    pub fn into_parts(self) -> (P, GasDiscreteGradientStep) {
        (self.endpoint, self.step)
    }
}

fn invalid(message: &'static str) -> GasEnergyError { GasEnergyError::Invalid(message) }

// The same ordered log1p correction as the original entropy source. Positive
// inputs avoid ratio overflow; reversal is exactly antisymmetric.
fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b { return 0.0; }
    let (lo, hi, sign) = if a > b { (b, a, 1.0) } else { (a, b, -1.0) };
    let delta = (hi - lo) / lo;
    let rounded = 1.0 + delta;
    let value = if rounded == 1.0 { delta }
        else if rounded.is_finite() {
            fs_math::det::ln(rounded) * (delta / (rounded - 1.0))
        } else { fs_math::det::ln(hi) - fs_math::det::ln(lo) };
    sign * value
}

// a=ln(V1/V0)/(V1-V0), and its derivative at fixed V0. The local
// series removes cancellation in d/V1-ln(V1/V0); remainder is O(x^10).
fn inverse_volume_mean(v0: f64, v1: f64, logarithm: f64) -> (f64, f64) {
    let delta = v1 - v0;
    if delta == 0.0 { return (1.0 / v0, (-0.5 / v0) / v0); }
    let x = delta / v0;
    let derivative = if x.abs() < 0.01 {
        let coefficients = [-0.5, 2.0/3.0, -3.0/4.0, 4.0/5.0, -5.0/6.0,
            6.0/7.0, -7.0/8.0, 8.0/9.0, -9.0/10.0, 10.0/11.0];
        let value = coefficients.into_iter().rev().fold(0.0, |sum, c| sum * x + c);
        (value / v0) / v0
    } else { ((delta / v1 - logarithm) / delta) / delta };
    (logarithm / delta, derivative)
}

struct CaloricMean {
    energy: f64,
    entropy: f64,
    temperature: f64,
    derivative: f64,
    cv: f64,
}

fn cv_at(model: &ThermochemicalGasModel, t: f64) -> Result<f64, GasEnergyError> {
    let pressure = model.thermodynamics().components()[0].convention().reference_pressure();
    let properties = model.thermodynamics().evaluate_properties(Temperature::new(t), pressure)
        .map_err(|error| failure("discrete Gibbs capacity", error))?;
    let cv = properties.mass_specific_cv().value();
    if !positive(cv) { return Err(invalid("positive original gas Cv required")); }
    Ok(cv)
}

// Integral_0^1 Cv(T0+z*dT)*(1-z)/(T0+z*dT) dz. Evaluating this positive
// moment, rather than T1*Delta s_cv-Delta u, retains tiny-step sensitivities.
// Four-point Gauss is exact for the polynomial part of NASA-9 here; its
// rational terms are smooth on the relative temperature span <=1e-3.
// Split every original source seam, with a hard bound on local work.
fn close_temperature_moment(
    model: &ThermochemicalGasModel, t0: f64, t1: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GasEnergyError> {
    let delta = t1 - t0;
    let nodes = [0.06943184420297371, 0.33000947820757187,
        0.6699905217924281, 0.9305681557970262];
    let weights = [0.17392742256872692, 0.32607257743127305,
        0.32607257743127305, 0.17392742256872692];
    let mut left = 0.0_f64;
    let mut moment = 0.0;
    for _ in 0..128 {
        poll(cancelled)?;
        let mut right = 1.0_f64;
        for species in model.thermodynamics().components() {
            for region in species.regions() {
                let z = (region.temperature_min().value() - t0) / delta;
                if z > left && z < right { right = z; }
            }
        }
        let a = t0 + left * delta;
        let b = t0 + right * delta;
        let lo = a.min(b);
        let hi = a.max(b);
        for (node, weight) in nodes.into_iter().zip(weights) {
            poll(cancelled)?;
            let z = left + (right - left) * node;
            let t = t0 + delta * z;
            // A one-ulp segment has no interior temperature. The inclusive
            // lower endpoint selects that segment's original source region.
            let sample = if lo.next_up() >= hi { lo }
                else { t.max(lo.next_up()).min(hi.next_down()) };
            moment += (right - left) * weight * (cv_at(model, sample)? / t) * (1.0 - z);
        }
        if right == 1.0 { return Ok(moment); }
        left = right;
    }
    Err(invalid("too many source seams in discrete Gibbs tangent; narrow temperature bracket"))
}

fn caloric_mean(
    model: &ThermochemicalGasModel, initial: f64, temperature: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<CaloricMean, GasEnergyError> {
    poll(cancelled)?;
    let cv = cv_at(model, temperature)?;
    if initial == temperature {
        return Ok(CaloricMean { energy: 0.0, entropy: 0.0,
            temperature, derivative: 0.5, cv });
    }
    let source = model.thermodynamics();
    let start = Temperature::new(initial);
    let end = Temperature::new(temperature);
    let energy = source.specific_internal_energy_change(start, end)
        .map_err(|error| failure("discrete Gibbs energy increment", error))?.value();
    poll(cancelled)?;
    let entropy = source.specific_isochoric_entropy_change(start, end)
        .map_err(|error| failure("discrete Gibbs entropy increment", error))?.value();
    let mean = energy / entropy;
    let delta = temperature - initial;
    let derivative = if delta.abs() <= 1e-3 * initial.min(temperature) {
        let divided_entropy = entropy / delta;
        (cv / temperature) * (close_temperature_moment(model, initial, temperature, cancelled)?
            / divided_entropy) / divided_entropy
    } else { (cv / entropy) * ((temperature - mean) / temperature) };
    if !energy.is_finite() || !entropy.is_finite() || energy == 0.0 || entropy == 0.0
        || !positive(mean) || !positive(derivative)
        || mean < initial.min(temperature) * (1.0 - 32.0 * f64::EPSILON)
        || mean > initial.max(temperature) * (1.0 + 32.0 * f64::EPSILON)
    {
        return Err(invalid("unrepresentable original-law thermodynamic mean"));
    }
    Ok(CaloricMean { energy, entropy, temperature: mean, derivative, cv })
}

impl ClosedGasParcel {
    /// Tentatively advance simultaneous heat and volume work using original
    /// internal-energy and entropy increments. Mass, composition and auxiliary
    /// query axes are unchanged; endpoint and mean transport are fully resolved.
    ///
    /// The caloric branch is prepared once. Its positive-Cv bound also admits
    /// an increasing discrete residual: compression requires
    /// -R*ln(V1/V0) < Cv_min*Tmin/T0. Refine excessive contraction steps rather
    /// than silently accepting a nonmonotone root. Expansion has no such limit.
    /// The source-temperature bracket still bounds every actual endpoint.
    /// Indeed theta'/Cv(T) <= T0/(Cv_min*Tmin), so the compression residual
    /// derivative is at least Cv_min + R*ln(V1/V0)*T0/Tmin > 0.
    ///
    /// The returned mean pressure MUST drive coupled work. Absolute energy,
    /// Gibbs and entropy-transfer residuals all pass the energy budget before
    /// return. No mutation occurs on source refusal, cancellation or budget
    /// exhaustion. At zero heat and zero motion the endpoint is bit-preserved.
    pub fn discrete_gradient_heat_volume(
        &self, volume_m3: f64, heat_j: f64, controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasDiscreteGradientCandidate, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let old = &self.parcel;
        let t0 = old.state.state().temperature;
        if !positive(volume_m3) || !heat_j.is_finite()
            || t0 < controls.minimum_temperature_k || t0 > controls.maximum_temperature_k
        {
            return Err(invalid("positive volume, finite heat and bracket containing old temperature required"));
        }
        let mass = old.mass_kg;
        let r = old.state.state().specific_gas_constant;
        let dv = volume_m3 - old.volume_m3;
        let logarithm = log_ratio(volume_m3, old.volume_m3);
        let (a, da) = inverse_volume_mean(old.volume_m3, volume_m3, logarithm);
        let density = mass / volume_m3;
        let q = heat_j / mass;
        let tolerance = mass * controls.energy_tolerance_j_kg;
        if !positive(a) || !da.is_finite() || da >= 0.0 || !positive(density)
            || !q.is_finite() || !positive(tolerance) || (heat_j != 0.0 && q == 0.0)
        {
            return Err(invalid("discrete gas geometry, heat or tolerance is unrepresentable"));
        }
        let k = r * logarithm;
        let unchanged = dv == 0.0 && heat_j == 0.0;
        let initial_residual = k * t0 - q;
        let mut iterations = 0;
        let t1 = if unchanged { t0 } else {
            let branch = old.model.prepare_caloric_branch(controls, &mut cancelled)?;
            let minimum_cv = branch.minimum_cv_j_kg_k();
            let compression_limit = (minimum_cv * (controls.minimum_temperature_k / t0)).next_down();
            if k < 0.0 && -k >= compression_limit {
                return Err(invalid("discrete gas compression exceeds positive-tangent branch; refine volume step"));
            }
            let minimum_tangent = if k < 0.0 {
                minimum_cv + k * (t0 / controls.minimum_temperature_k)
            } else { minimum_cv };
            if !positive(minimum_tangent) || !initial_residual.is_finite() {
                return Err(invalid("discrete gas thermal residual or tangent bound is unrepresentable"));
            }
            // Heat and boundary work may cancel to operation roundoff even
            // when their separately rounded per-mass values differ. This
            // allowance cannot swallow pure heating: with k=0 it is <|q|.
            let isothermal_tolerance = (8.0 * f64::EPSILON
                * q.abs().max((k * t0).abs())).min(0.25 * controls.energy_tolerance_j_kg);
            if initial_residual.abs() <= isothermal_tolerance { t0 } else {
                let span = (initial_residual.abs() / minimum_tangent).next_up();
                let (lo, hi) = if initial_residual < 0.0 {
                    (t0, controls.maximum_temperature_k.min((t0 + span).next_up()))
                } else {
                    (controls.minimum_temperature_k.max((t0 - span).next_down()), t0)
                };
                let solve_tolerance = (0.25 * controls.energy_tolerance_j_kg)
                    .min(0.5 * initial_residual.abs());
                if !positive(solve_tolerance) {
                    return Err(invalid("nonzero discrete gas increment is below represented solve resolution"));
                }
                let evaluate = |t| -> Result<root::Sample, GasEnergyError> {
                    let c = caloric_mean(&old.model, t0, t, &mut cancelled)?;
                    Ok(root::Sample { residual: c.energy + k * c.temperature - q,
                        tangent: c.cv + k * c.derivative, pressure_pa: density * r * t })
                };
                // The evaluator polls the real callback inside its bounded
                // source/seam work. Do not alias that same mutable callback
                // through the root driver's separate between-sample poll.
                let solved = root::solve(lo, hi, solve_tolerance, controls.maximum_iterations,
                    evaluate, &mut || false).map_err(map_root)?;
                iterations = solved.iterations;
                solved.temperature_k
            }
        };
        poll(&mut cancelled)?;
        let c = caloric_mean(&old.model, t0, t1, &mut cancelled)?;
        let mean_density = mass * a;
        if !positive(mean_density) {
            return Err(invalid("thermodynamic mean density is unrepresentable"));
        }
        let (endpoint, mean) = if unchanged { (self.clone(), old.state.clone()) } else {
            let point = updated_point(&old.point, t1, density * r * t1)?;
            let state = old.model.evaluate(&point)?;
            poll(&mut cancelled)?;
            let mean_point = updated_point(&old.point, c.temperature, mean_density * r * c.temperature)?;
            let mean = old.model.evaluate(&mean_point)?;
            if (state.state().density / density - 1.0).abs() > 64.0 * f64::EPSILON
                || (mean.state().density / mean_density - 1.0).abs() > 64.0 * f64::EPSILON
            {
                return Err(invalid("discrete gas material states did not retain their conserved density"));
            }
            let mut parcel = old.clone();
            parcel.point = point;
            parcel.state = state;
            parcel.volume_m3 = volume_m3;
            (Self::from_isochoric(parcel), mean)
        };
        let work = mean.state().pressure * dv;
        let stored = mass * c.energy;
        let entropy = mass * (c.entropy + k);
        let transferred = heat_j / c.temperature;
        let energy_residual = (stored + work) - heat_j;
        let entropy_residual = entropy - transferred;
        let gibbs_residual = (stored + work) - c.temperature * entropy;
        let tangent = c.cv + k * c.derivative;
        let dpdv = (mass * r * c.temperature)
            * (da - a * c.derivative * r / (volume_m3 * tangent));
        let dpdq = r * a * c.derivative / tangent;
        if !work.is_finite() || (dv != 0.0 && work == 0.0)
            || !energy_residual.is_finite() || energy_residual.abs() > tolerance
            || !gibbs_residual.is_finite() || gibbs_residual.abs() > tolerance
            || !entropy_residual.is_finite()
            || entropy_residual.abs() > tolerance / c.temperature
            || !positive(tangent) || !dpdv.is_finite() || dpdv >= 0.0 || !positive(dpdq)
        {
            return Err(invalid("discrete gas Gibbs, energy, entropy or sensitivity budget refused"));
        }
        poll(&mut cancelled)?;
        let work_pressure_pa = mean.state().pressure;
        Ok(GasDiscreteGradientCandidate { endpoint, thermodynamic_mean: mean,
            step: GasDiscreteGradientStep {
                initial_volume_m3: old.volume_m3, final_volume_m3: volume_m3,
                supplied_heat_j: heat_j, work_pressure_pa,
                work_by_gas_j: work, stored_energy_change_j: stored,
                energy_residual_j: energy_residual, entropy_change_j_k: entropy,
                entropy_transfer_j_k: transferred, entropy_residual_j_k: entropy_residual,
                gibbs_residual_j: gibbs_residual,
                pressure_final_volume_derivative_pa_m3: dpdv,
                pressure_heat_derivative_pa_j: dpdq, iterations,
            } })
    }
}

#[cfg(test)]
mod tests;
