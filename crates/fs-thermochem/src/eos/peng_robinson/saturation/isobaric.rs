//! Original pure-fluid sensible/latent calorics at one coexistence pressure.
//!
//! Saturated-liquid enthalpy is the zero coordinate. Between zero and the
//! original latent enthalpy, conserved H determines phase amount at constant T.
//! Outside that interval the original stable PR roots and NASA integrals supply
//! sensible calorics. No apparent heat capacity, new EOS or heat integrator is
//! introduced. The phase model excludes solids, kinetics and critical merging.

use super::super::{
    MixingRecord, PengRobinsonErrorV1, PengRobinsonMixtureV1, PengRobinsonPhaseRootV1,
    PengRobinsonPhaseV1,
};
use super::{PengRobinsonSaturatedMixtureV1, PengRobinsonSaturationV1};
use crate::{
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStatePhaseV1, ThermochemErrorV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    kinetics::{KineticsErrorV1, admit_caloric_interval},
};
use fs_math::dd::Dd;
use fs_qty::{Density, Pressure, Temperature};

/// Original-source inversion accuracy and bounded new PT work per query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PureFluidIsobaricSolveV1 {
    /// Absolute specific enthalpy residual [J/kg], independent of formation data.
    pub enthalpy_tolerance_j_kg: f64,
    /// Final Newton correction or bracket width [K]; H must also converge.
    pub temperature_tolerance_k: f64,
    /// Maximum original PT root evaluations per inversion, from 1 through 256.
    pub maximum_evaluations: usize,
}
impl Default for PureFluidIsobaricSolveV1 {
    fn default() -> Self {
        Self {
            enthalpy_tolerance_j_kg: 1e-6,
            temperature_tolerance_k: 1e-9,
            maximum_evaluations: 128,
        }
    }
}
impl PureFluidIsobaricSolveV1 {
    /// Validate finite physical tolerances and the bounded work policy.
    pub fn validate(self) -> Result<(), PureFluidIsobaricErrorV1> {
        if !positive(self.enthalpy_tolerance_j_kg)
            || !positive(self.temperature_tolerance_k)
            || !(1..=256).contains(&self.maximum_evaluations)
        {
            return Err(E::Invalid(
                "positive isobaric tolerances and bounded PT work required",
            ));
        }
        Ok(())
    }
}

/// Refusal never changes the immutable source or publishes a partial result.
#[derive(Clone, Debug, PartialEq)]
pub enum PureFluidIsobaricErrorV1 {
    /// Invalid source identity, physical domain or unresolved arithmetic.
    Invalid(&'static str),
    /// Original EOS refusal.
    Eos(PengRobinsonErrorV1),
    /// Original NASA evaluation refusal.
    Source(ThermochemErrorV1),
    /// Original continuous caloric source admission refusal.
    Caloric(KineticsErrorV1),
    /// T alone leaves the phase amount and H undetermined at coexistence.
    CoexistenceAmountRequired,
    /// A coordinate lies outside the prepared physical interval.
    OutsideDomain,
    /// The admitted PT budget or representable bracket was exhausted.
    NotConverged {
        evaluations: usize,
        residual_j_kg: f64,
    },
    /// Cancellation before a bounded evaluation or final publication.
    Cancelled,
}
impl core::fmt::Display for PureFluidIsobaricErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pure PR isobaric calorics: {self:?}")
    }
}
impl std::error::Error for PureFluidIsobaricErrorV1 {}
impl From<PengRobinsonErrorV1> for PureFluidIsobaricErrorV1 {
    fn from(error: PengRobinsonErrorV1) -> Self {
        Self::Eos(error)
    }
}
impl From<ThermochemErrorV1> for PureFluidIsobaricErrorV1 {
    fn from(error: ThermochemErrorV1) -> Self {
        Self::Source(error)
    }
}
type E = PureFluidIsobaricErrorV1;

/// Original homogeneous root or original coexistence with an H-selected amount.
#[derive(Clone, Debug, PartialEq)]
pub enum PureFluidIsobaricPhaseV1 {
    /// Stable liquid at or below the transition.
    Liquid(PengRobinsonPhaseRootV1),
    /// Original saturated roots; mass and mole quality coincide for a pure fluid.
    Coexistence(PengRobinsonSaturatedMixtureV1),
    /// Stable vapor at or above the transition.
    Vapor(PengRobinsonPhaseRootV1),
}

/// Source-relative physical state and explicitly one-sided caloric derivatives.
#[derive(Clone, Debug, PartialEq)]
pub struct PureFluidIsobaricStateV1 {
    /// Original absolute temperature.
    pub temperature: Temperature,
    /// Source enthalpy relative to the saturated liquid [J/kg].
    pub specific_enthalpy_j_kg: f64,
    /// Actual original source entropy relative to the saturated liquid [J/(kg K)].
    pub source_entropy_change_j_kg_k: f64,
    /// Conserved mass divided by the actual original phase volume [kg/m3].
    pub density_kg_m3: f64,
    /// dT/dh [kg K/J]; exactly zero on the inclusive latent interval.
    pub temperature_slope_kg_k_j: f64,
    /// d(1/rho)/dh [m3/J], including latent expansion on the plateau.
    pub specific_volume_slope_m3_j: f64,
    /// Complete original phase roots or coexistence allocation.
    pub phase: PureFluidIsobaricPhaseV1,
}
impl PureFluidIsobaricStateV1 {
    /// Actual vapor amount, without interpreting density as phase fraction.
    pub fn vapor_fraction(&self) -> f64 {
        match &self.phase {
            PureFluidIsobaricPhaseV1::Liquid(_) => 0.0,
            PureFluidIsobaricPhaseV1::Vapor(_) => 1.0,
            PureFluidIsobaricPhaseV1::Coexistence(state) => state.vapor_fraction(),
        }
    }
}

/// A bounded H inverse retains its actual source residual and new EOS work.
#[derive(Clone, Debug, PartialEq)]
pub struct PureFluidIsobaricSolutionV1 {
    pub state: PureFluidIsobaricStateV1,
    pub residual_j_kg: f64,
    /// New PT root evaluations; cached endpoints and latent queries use zero.
    pub evaluations: usize,
}

/// Direct original-source changes between two canonical inversions.
#[derive(Clone, Debug, PartialEq)]
pub struct PureFluidIsobaricChangeV1 {
    pub initial: PureFluidIsobaricSolutionV1,
    pub final_state: PureFluidIsobaricSolutionV1,
    pub requested_enthalpy_change_j_kg: f64,
    pub source_enthalpy_change_j_kg: f64,
    pub source_entropy_change_j_kg_k: f64,
    /// Sensible source entropy plus latent Delta h/T. The visible coexistence
    /// fugacity defect is not silently fitted into the original source entropy.
    pub caloric_entropy_change_j_kg_k: f64,
}

/// Prepared pure-fluid path retaining the identical original EOS/NASA sources.
///
/// The admitted PR alpha-curvature sign proves residual Cv is nonnegative on
/// this subcritical interval. Original continuous NASA admission proves ideal
/// Cv positive; stable-root Cp >= Cv. A consumer needing a numerical global
/// inverse-slope bound must additionally bound that original ideal Cv chart.
pub struct PengRobinsonIsobaricCaloricV1<'a> {
    model: &'a PengRobinsonMixtureV1,
    source: &'a Nasa9StandardStateModelV1,
    saturation: PengRobinsonSaturationV1,
    controls: PureFluidIsobaricSolveV1,
    domain: [Temperature; 2],
    endpoints: [PureFluidIsobaricStateV1; 2],
    latent_h: f64,
    latent_s: f64,
}
impl<'a> PengRobinsonIsobaricCaloricV1<'a> {
    /// Bind an already accepted original coexistence point to adjacent sensible
    /// branches. The transition must be strictly inside a subcritical continuous
    /// NASA interval. The pressure is the ORIGINAL accepted coexistence pressure.
    pub fn new(
        model: &'a PengRobinsonMixtureV1,
        source: &'a Nasa9StandardStateModelV1,
        saturation: PengRobinsonSaturationV1,
        domain: [Temperature; 2],
        controls: PureFluidIsobaricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, E> {
        poll(&mut cancelled)?;
        controls.validate()?;
        let ts = saturation.temperature().value();
        if model.components().len() != 1
            || model.mole_fractions() != [1.0]
            || model.components()[0] != *saturation.component()
            || model.molar_mass() != source.molar_mass()
            || source.convention().phase() != StandardStatePhaseV1::Gas
            || source.convention().eos() != ReferenceEquationOfStateV1::IdealGas
            || !positive(domain[0].value())
            || !(domain[0].value() < ts && ts < domain[1].value())
            || domain[1].value() >= saturation.component().critical_temperature.value()
        {
            return Err(E::Invalid(
                "same pure source and a subcritical bracket straddling coexistence required",
            ));
        }
        admit_caloric_interval(source, domain[0], domain[1], &mut cancelled).map_err(|error| {
            if matches!(error, KineticsErrorV1::Cancelled) {
                E::Cancelled
            } else {
                E::Caloric(error)
            }
        })?;
        // The original component traversal owns kappa. For its pure alpha,
        // a''(T) is proportional to kappa*(1+kappa), hence residual Cv >= 0
        // when kappa >= 0. Below Tc the alpha branch is strictly positive.
        let mut nonnegative = false;
        model.mixing_caloric(domain[0].value(), |record| {
            if let MixingRecord::Component { kappa, .. } = record {
                nonnegative = kappa >= 0.0;
            }
        })?;
        if !nonnegative {
            return Err(E::Invalid(
                "original PR alpha curvature does not establish a nonnegative residual Cv",
            ));
        }
        let mass = source.molar_mass().value();
        let latent_h = saturation.latent_enthalpy_j_mol() / mass;
        let latent_s = (saturation.vapor().residual_entropy_j_mol_k()
            - saturation.liquid().residual_entropy_j_mol_k())
            / mass;
        if !positive(latent_h) || !positive(latent_s) {
            return Err(E::Invalid(
                "positive represented original latent H and S required",
            ));
        }
        let placeholder = plateau(&saturation, latent_h, latent_s, 0.0)?;
        let mut path = Self {
            model,
            source,
            saturation,
            controls,
            domain,
            endpoints: [placeholder.clone(), placeholder],
            latent_h,
            latent_s,
        };
        path.endpoints = [
            path.sensible(domain[0], false, &mut cancelled)?,
            path.sensible(domain[1], true, &mut cancelled)?,
        ];
        if path.endpoints[0].specific_enthalpy_j_kg >= 0.0
            || path.endpoints[1].specific_enthalpy_j_kg <= latent_h
        {
            return Err(E::Invalid(
                "sensible phase spans must be represented on both sides of latent storage",
            ));
        }
        poll(&mut cancelled)?;
        Ok(path)
    }
    pub const fn source(&self) -> &'a Nasa9StandardStateModelV1 {
        self.source
    }
    pub const fn constitutive_law(&self) -> &'a PengRobinsonMixtureV1 {
        self.model
    }
    pub const fn saturation(&self) -> &PengRobinsonSaturationV1 {
        &self.saturation
    }
    pub const fn controls(&self) -> PureFluidIsobaricSolveV1 {
        self.controls
    }
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.domain
    }
    pub fn pressure(&self) -> Pressure {
        self.saturation.pressure()
    }
    pub fn specific_enthalpy_domain_j_kg(&self) -> [f64; 2] {
        self.endpoints
            .each_ref()
            .map(|state| state.specific_enthalpy_j_kg)
    }
    pub const fn latent_enthalpy_j_kg(&self) -> f64 {
        self.latent_h
    }
    /// Original latent H-T*S closure [J/kg], unchanged by the thermal coordinate.
    pub fn latent_closure_residual_j_kg(&self) -> f64 {
        (Dd::from_f64(self.latent_h)
            - Dd::from_f64(self.saturation.temperature().value()) * Dd::from_f64(self.latent_s))
        .to_f64()
    }

    /// T alone admits sensible branches; the transition itself is set-valued.
    pub fn sample_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidIsobaricStateV1, E> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        if !t.is_finite() || t < self.domain[0].value() || t > self.domain[1].value() {
            return Err(E::OutsideDomain);
        }
        if temperature == self.saturation.temperature() {
            return Err(E::CoexistenceAmountRequired);
        }
        for endpoint in &self.endpoints {
            if temperature == endpoint.temperature {
                poll(&mut cancelled)?;
                return Ok(endpoint.clone());
            }
        }
        self.sensible(
            temperature,
            t > self.saturation.temperature().value(),
            &mut cancelled,
        )
    }

    /// H selects quality exactly on the plateau; otherwise bounded original
    /// Newton/bisection queries invert the stable sensible branch. Every query
    /// begins from immutable prepared endpoints, independent of prior calls.
    pub fn sample_specific_enthalpy(
        &self,
        h: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidIsobaricSolutionV1, E> {
        poll(&mut cancelled)?;
        let [lo, hi] = self.specific_enthalpy_domain_j_kg();
        if !h.is_finite() || h < lo || h > hi {
            return Err(E::OutsideDomain);
        }
        let publish = |state: PureFluidIsobaricStateV1, evaluations| {
            let residual_j_kg = state.specific_enthalpy_j_kg - h;
            PureFluidIsobaricSolutionV1 {
                state,
                residual_j_kg,
                evaluations,
            }
        };
        if (0.0..=self.latent_h).contains(&h) {
            let quality = h / self.latent_h;
            if (h != 0.0 && quality == 0.0) || (h != self.latent_h && quality == 1.0) {
                return Err(E::Invalid("nonzero latent phase amount is unrepresentable"));
            }
            let state = plateau(&self.saturation, self.latent_h, self.latent_s, quality)?;
            let answer = publish(state, 0);
            if answer.residual_j_kg.abs() > self.controls.enthalpy_tolerance_j_kg {
                return Err(E::Invalid(
                    "latent fraction cannot resolve requested enthalpy",
                ));
            }
            poll(&mut cancelled)?;
            return Ok(answer);
        }
        let vapor = h > self.latent_h;
        let endpoint = &self.endpoints[usize::from(vapor)];
        if h == endpoint.specific_enthalpy_j_kg {
            poll(&mut cancelled)?;
            return Ok(publish(endpoint.clone(), 0));
        }
        let root = if vapor {
            self.saturation.vapor()
        } else {
            self.saturation.liquid()
        };
        let boundary = self.root_state(self.saturation.temperature(), root.clone(), vapor)?;
        let (mut lower, mut upper) = if vapor {
            (boundary, endpoint.clone())
        } else {
            (endpoint.clone(), boundary)
        };
        let mut residual = f64::INFINITY;
        let mut evaluations = 0;
        while evaluations < self.controls.maximum_evaluations {
            let tl = lower.temperature.value();
            let th = upper.temperature.value();
            let best = if (lower.specific_enthalpy_j_kg - h).abs()
                <= (upper.specific_enthalpy_j_kg - h).abs()
            {
                &lower
            } else {
                &upper
            };
            let correction = (best.specific_enthalpy_j_kg - h) * best.temperature_slope_kg_k_j;
            let newton = best.temperature.value() - correction;
            let t = if evaluations % 4 != 3 && newton > tl && newton < th {
                newton
            } else {
                tl + 0.5 * (th - tl)
            };
            if t <= tl || t >= th {
                break;
            }
            evaluations += 1;
            let state = self.sensible(Temperature::new(t), vapor, &mut cancelled)?;
            residual = state.specific_enthalpy_j_kg - h;
            if residual.abs() <= self.controls.enthalpy_tolerance_j_kg
                && (residual.abs() * state.temperature_slope_kg_k_j
                    <= self.controls.temperature_tolerance_k
                    || th - tl <= self.controls.temperature_tolerance_k)
            {
                poll(&mut cancelled)?;
                return Ok(publish(state, evaluations));
            }
            if residual < 0.0 {
                lower = state;
            } else {
                upper = state;
            }
        }
        Err(E::NotConverged {
            evaluations,
            residual_j_kg: residual,
        })
    }

    /// Direct source differences preserve small sensible heat and the original
    /// latent source entropy. No formation-loaded absolute total is subtracted.
    pub fn change_between_specific_enthalpies(
        &self,
        a: f64,
        b: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidIsobaricChangeV1, E> {
        let initial = self.sample_specific_enthalpy(a, &mut cancelled)?;
        let final_state = self.sample_specific_enthalpy(b, &mut cancelled)?;
        let (h, source_s, caloric_s) = self.change(&initial.state, &final_state.state)?;
        poll(&mut cancelled)?;
        Ok(PureFluidIsobaricChangeV1 {
            initial,
            final_state,
            requested_enthalpy_change_j_kg: (Dd::from_f64(b) - Dd::from_f64(a)).to_f64(),
            source_enthalpy_change_j_kg: h,
            source_entropy_change_j_kg_k: source_s,
            caloric_entropy_change_j_kg_k: caloric_s,
        })
    }

    fn sensible(
        &self,
        temperature: Temperature,
        vapor: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PureFluidIsobaricStateV1, E> {
        poll(cancelled)?;
        let evaluated = self.model.evaluate_pt(temperature, self.pressure())?;
        let root = evaluated.select(if vapor {
            PengRobinsonPhaseV1::Vapor
        } else {
            PengRobinsonPhaseV1::Liquid
        });
        let stable = evaluated.select(PengRobinsonPhaseV1::Stable);
        let gap = root.residual_gibbs_over_rt() - stable.residual_gibbs_over_rt();
        let allowance = self.saturation.log_fugacity_residual().abs()
            + 64.0 * f64::EPSILON * root.residual_gibbs_over_rt().abs().max(1.0);
        if !gap.is_finite() || gap > allowance {
            return Err(E::Invalid(
                "sensible root is not the original stable fluid branch",
            ));
        }
        let state = self.root_state(temperature, root.clone(), vapor)?;
        poll(cancelled)?;
        Ok(state)
    }

    fn root_state(
        &self,
        t: Temperature,
        root: PengRobinsonPhaseRootV1,
        vapor: bool,
    ) -> Result<PureFluidIsobaricStateV1, E> {
        let reference = self.source.evaluate_properties(t)?;
        let mass = self.source.molar_mass().value();
        let response = root
            .caloric_mechanics(reference.cp().value() - R, mass)
            .ok_or(E::Invalid(
                "original sensible phase has no stable caloric response",
            ))?;
        let anchor = if vapor {
            self.saturation.vapor()
        } else {
            self.saturation.liquid()
        };
        let (h, s) = self.homogeneous_change(self.saturation.temperature(), anchor, t, &root)?;
        let specific_enthalpy_j_kg =
            (Dd::from_f64(h) + Dd::from_f64(if vapor { self.latent_h } else { 0.0 })).to_f64();
        let source_entropy_change_j_kg_k =
            (Dd::from_f64(s) + Dd::from_f64(if vapor { self.latent_s } else { 0.0 })).to_f64();
        let slope = mass / response.molar_cp_j_mol_k;
        let dv = -(response.dp_dt_v_pa_k / response.dp_dv_t_pa_mol_m3) / response.molar_cp_j_mol_k;
        if !specific_enthalpy_j_kg.is_finite()
            || !source_entropy_change_j_kg_k.is_finite()
            || !positive(slope)
            || !dv.is_finite()
            || (response.dp_dt_v_pa_k != 0.0 && dv == 0.0)
        {
            return Err(E::Invalid(
                "unrepresentable original isobaric calorics or derivative",
            ));
        }
        Ok(PureFluidIsobaricStateV1 {
            temperature: t,
            specific_enthalpy_j_kg,
            source_entropy_change_j_kg_k,
            density_kg_m3: root.mass_density_kg_per_m3(),
            temperature_slope_kg_k_j: slope,
            specific_volume_slope_m3_j: dv,
            phase: if vapor {
                PureFluidIsobaricPhaseV1::Vapor(root)
            } else {
                PureFluidIsobaricPhaseV1::Liquid(root)
            },
        })
    }

    fn homogeneous_change(
        &self,
        t0: Temperature,
        a: &PengRobinsonPhaseRootV1,
        t1: Temperature,
        b: &PengRobinsonPhaseRootV1,
    ) -> Result<(f64, f64), E> {
        if t0 == t1 && a == b {
            return Ok((0.0, 0.0));
        }
        let change = self
            .model
            .prepare_density_caloric_change(
                t0,
                Density::new(a.mass_density_kg_per_m3()),
                Density::new(b.mass_density_kg_per_m3()),
            )?
            .evaluate(t1)?;
        let mass = Dd::from_f64(self.source.molar_mass().value());
        let h = ((Dd::from_f64(self.source.molar_enthalpy_change(t0, t1)?.value())
            + Dd::from_f64(change.residual_enthalpy_change_j_mol()))
            / mass)
            .to_f64();
        let s = ((Dd::from_f64(self.source.molar_entropy_change(t0, t1)?.value())
            + Dd::from_f64(change.residual_entropy_at_density_change_j_mol_k())
            + Dd::from_f64(change.ideal_entropy_density_correction_j_mol_k()))
            / mass)
            .to_f64();
        if !h.is_finite() || !s.is_finite() {
            return Err(E::Invalid("unrepresentable original source increment"));
        }
        Ok((h, s))
    }

    fn change(
        &self,
        a: &PureFluidIsobaricStateV1,
        b: &PureFluidIsobaricStateV1,
    ) -> Result<(f64, f64, f64), E> {
        match (&a.phase, &b.phase) {
            (PureFluidIsobaricPhaseV1::Liquid(ar), PureFluidIsobaricPhaseV1::Liquid(br))
            | (PureFluidIsobaricPhaseV1::Vapor(ar), PureFluidIsobaricPhaseV1::Vapor(br)) => {
                let (h, s) = self.homogeneous_change(a.temperature, ar, b.temperature, br)?;
                return Ok((h, s, s));
            }
            _ => {}
        }
        let leg = |state: &PureFluidIsobaricStateV1| -> Result<(f64, f64), E> {
            match &state.phase {
                PureFluidIsobaricPhaseV1::Liquid(root) => self.homogeneous_change(
                    self.saturation.temperature(),
                    self.saturation.liquid(),
                    state.temperature,
                    root,
                ),
                PureFluidIsobaricPhaseV1::Vapor(root) => self.homogeneous_change(
                    self.saturation.temperature(),
                    self.saturation.vapor(),
                    state.temperature,
                    root,
                ),
                PureFluidIsobaricPhaseV1::Coexistence(_) => Ok((0.0, 0.0)),
            }
        };
        let (ha, sa) = leg(a)?;
        let (hb, sb) = leg(b)?;
        let q = Dd::from_f64(b.vapor_fraction()) - Dd::from_f64(a.vapor_fraction());
        let dh = Dd::from_f64(hb) - Dd::from_f64(ha) + q * Dd::from_f64(self.latent_h);
        let sensible_s = Dd::from_f64(sb) - Dd::from_f64(sa);
        let source_s = sensible_s + q * Dd::from_f64(self.latent_s);
        let caloric_s = sensible_s
            + q * Dd::from_f64(self.latent_h) / Dd::from_f64(self.saturation.temperature().value());
        let result = (dh.to_f64(), source_s.to_f64(), caloric_s.to_f64());
        if ![result.0, result.1, result.2]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(E::Invalid("unrepresentable sensible/latent interval"));
        }
        Ok(result)
    }
}

fn plateau(
    saturation: &PengRobinsonSaturationV1,
    latent_h: f64,
    latent_s: f64,
    quality: f64,
) -> Result<PureFluidIsobaricStateV1, E> {
    let phase = saturation
        .at_vapor_fraction(quality)
        .map_err(|_| E::Invalid("latent enthalpy cannot resolve a physical phase amount"))?;
    let specific_volume_slope_m3_j = (saturation.vapor().molar_volume_m3_per_mol()
        - saturation.liquid().molar_volume_m3_per_mol())
        / saturation.latent_enthalpy_j_mol();
    let density_kg_m3 = phase.mass_density_kg_m3();
    if !positive(specific_volume_slope_m3_j) || !positive(density_kg_m3) {
        return Err(E::Invalid("unrepresentable latent expansion"));
    }
    Ok(PureFluidIsobaricStateV1 {
        temperature: saturation.temperature(),
        specific_enthalpy_j_kg: quality * latent_h,
        source_entropy_change_j_kg_k: quality * latent_s,
        density_kg_m3,
        temperature_slope_kg_k_j: 0.0,
        specific_volume_slope_m3_j,
        phase: PureFluidIsobaricPhaseV1::Coexistence(phase),
    })
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), E> {
    if cancelled() {
        Err(E::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
