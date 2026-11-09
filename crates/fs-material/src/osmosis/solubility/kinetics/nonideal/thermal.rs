//! Temperature-qualified neutral activities and their calibrated aqueous heat.
//!
//! The immutable source is B(T) = B_ref + C (1/T - 1/T_ref), with C in
//! K m³/mol. The original solid caloric datum is retained. The AQUEOUS excess
//! potential is RT [nᵀ B n/(2V) - Σ n_i B_ii csat_i(T)]; its temperature
//! derivative supplies both excess mixing heat and the full diagonal reference
//! correction. Original dissolved Cp curves are baseline reference curves,
//! not independently measured effective mixture heat capacities.
//!
//! Endpoint convexity admits the complete inverse-temperature family on one
//! declared concentration box. Conservative interval bounds additionally
//! require each effective dissolved Cp to stay positive. These mathematical
//! admissions do not establish empirical applicability of caller coefficients.
//! Thermal trials reuse the original scalar phase kernel and heat controller;
//! they allocate nothing and never rebuild a dense thermodynamic model.

use super::super::{KineticSolubilityError, invalid};
use super::conditional_root;
use crate::osmosis::{MAX_OSMOTIC_SOLUTES, Sum, nonideal::VirialSolutionModel};
use fs_math::det::{exp, expm1};
use fs_qty::{AmountConcentration, Temperature, Volume};
use fs_thermochem::{
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    solution::solubility::{
        caloric::{
            SolubilityCaloricComponentV1, SolubilityCaloricErrorV1,
            SolubilityCaloricModelV1,
            kinetics::SolubilityThermalInteractionV1,
        },
        kinetics::{SolubilityKineticStateV1, SolubilityKineticsV1},
    },
};
use std::sync::Arc;

type Result<T> = core::result::Result<T, SolubilityCaloricErrorV1>;

/// Explicit, bounded inverse-temperature interaction source for the SAME
/// original ordered phase laws, reference heat capacities and solvent basis.
///
/// Pure-pair saturation remains the original law at every admitted temperature.
/// Its preservation adds an aqueous standard-state enthalpy and heat-capacity
/// correction; measured pure-pair dissolution enthalpy is not asserted to remain
/// unchanged by this newly selected nonideal thermal constitutive source.
#[derive(Clone, Debug, PartialEq)]
pub struct VirialSolubilityThermalModel {
    caloric: SolubilityCaloricModelV1,
    reference: Arc<VirialSolutionModel>,
    inverse_temperature_m3_k_mol: Box<[f64]>,
    temperature_domain: [Temperature; 2],
}

impl VirialSolubilityThermalModel {
    /// Prepare a symmetric family in the original source's species order.
    ///
    /// The reference temperature must lie inside the declared domain, which
    /// must fit every original phase and caloric source. Both endpoint matrices
    /// must admit the SAME original concentration caps. A bounded conservative
    /// Cp check may refuse a family whose positivity cannot be established;
    /// it never silently changes its coefficients or support.
    pub fn new(
        caloric: SolubilityCaloricModelV1,
        reference: Arc<VirialSolutionModel>,
        inverse_temperature_m3_k_mol: Vec<f64>,
        temperature_domain: [Temperature; 2],
    ) -> Result<Self> {
        let count = caloric.components().len();
        let [low, high] = temperature_domain.map(|t| t.value());
        let support = caloric.temperature_domain();
        let reference_t = reference.temperature_k();
        if count != reference.species_count()
            || inverse_temperature_m3_k_mol.len() != count * count
            || caloric.solvent_molar_volume_m3_mol().to_bits()
                != reference.solvent_molar_volume_m3_mol().to_bits()
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || low >= high
            || low < support[0].value()
            || high > support[1].value()
            || reference_t < low
            || reference_t > high
            || inverse_temperature_m3_k_mol.iter().any(|x| !x.is_finite())
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal virial source layout, solvent convention or temperature support",
            ));
        }
        for i in 0..count {
            for j in 0..i {
                if inverse_temperature_m3_k_mol[i * count + j]
                    != inverse_temperature_m3_k_mol[j * count + i]
                {
                    return Err(SolubilityCaloricErrorV1::Invalid(
                        "inverse-temperature interaction matrix must be exactly symmetric",
                    ));
                }
            }
        }
        let result = Self {
            caloric,
            reference,
            inverse_temperature_m3_k_mol: inverse_temperature_m3_k_mol.into_boxed_slice(),
            temperature_domain,
        };
        // Hessians are affine in 1/T. The cone of positive-definite matrices
        // is convex, so the two original box admissions cover every interior T.
        result.make_slice(temperature_domain[0])?;
        result.make_slice(temperature_domain[1])?;
        result.admit_caloric_domain()?;
        Ok(result)
    }

    pub const fn caloric(&self) -> &SolubilityCaloricModelV1 {
        &self.caloric
    }
    pub const fn reference_model(&self) -> &Arc<VirialSolutionModel> {
        &self.reference
    }
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.temperature_domain
    }
    pub fn inverse_temperature_m3_k_mol(&self) -> &[f64] {
        &self.inverse_temperature_m3_k_mol
    }

    /// Conservative retained-source charge, including shared source payloads.
    pub fn heap_bytes(&self) -> usize {
        core::mem::size_of::<Self>()
            + 2 * core::mem::size_of::<usize>()
            + core::mem::align_of::<Self>()
            + core::mem::size_of_val(self.inverse_temperature_m3_k_mol.as_ref())
            + core::mem::size_of::<VirialSolutionModel>()
            + 2 * core::mem::size_of::<usize>()
            + core::mem::align_of::<VirialSolutionModel>()
            + self.reference.heap_bytes()
            + self.caloric.heap_bytes()
    }

    /// Prepare a fixed-temperature view for the original material/transport
    /// owners. This preparation allocates; thermal trial callbacks do not.
    pub fn at_temperature(
        &self,
        temperature: Temperature,
    ) -> core::result::Result<Arc<VirialSolutionModel>, KineticSolubilityError> {
        self.check_temperature(temperature)
            .map_err(|_| invalid("temperature outside the original thermal virial source"))?;
        if temperature.value().to_bits() == self.reference.temperature_k().to_bits() {
            return Ok(Arc::clone(&self.reference));
        }
        self.make_slice(temperature)
            .map(Arc::new)
            .map_err(|_| invalid("thermal virial slice is not representable on its original support"))
    }

    fn make_slice(&self, temperature: Temperature) -> Result<VirialSolutionModel> {
        let mut matrix = Vec::with_capacity(self.inverse_temperature_m3_k_mol.len());
        for i in 0..self.inverse_temperature_m3_k_mol.len() {
            matrix.push(self.interaction(temperature.value(), i)?);
        }
        VirialSolutionModel::new(
            temperature.value(),
            self.reference.solvent_molar_volume_m3_mol(),
            self.reference.reference_concentration_mol_m3(),
            self.reference.maximum_concentrations_mol_m3().to_vec(),
            matrix,
        )
        .map_err(|_| {
            SolubilityCaloricErrorV1::Invalid(
                "thermal virial endpoint is outside its original convex concentration support",
            )
        })
    }

    fn interaction(&self, temperature: f64, index: usize) -> Result<f64> {
        let reference_t = self.reference.temperature_k();
        let change = inverse_difference(reference_t, temperature)?;
        finite(
            self.reference.interaction_matrix_m3_mol()[index]
                + self.inverse_temperature_m3_k_mol[index] * change,
        )
    }

    fn check_temperature(&self, temperature: Temperature) -> Result<()> {
        let t = temperature.value();
        if !t.is_finite()
            || t < self.temperature_domain[0].value()
            || t > self.temperature_domain[1].value()
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "temperature outside original nonideal thermal support",
            ));
        }
        Ok(())
    }

    fn concentrations(
        &self,
        temperature: Temperature,
        volume: Volume,
        phases: &[SolubilityKineticStateV1],
    ) -> Result<[f64; MAX_OSMOTIC_SOLUTES]> {
        self.check_temperature(temperature)?;
        if !volume.value().is_finite()
            || volume.value() <= 0.0
            || phases.len() != self.caloric.components().len()
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal virial phase inventory layout or solvent volume",
            ));
        }
        let mut concentrations = [0.0; MAX_OSMOTIC_SOLUTES];
        for ((out, phase), &cap) in concentrations
            .iter_mut()
            .zip(phases)
            .zip(self.reference.maximum_concentrations_mol_m3())
        {
            let amount = phase.dissolved().value();
            *out = finite(amount / volume.value())?;
            if *out < 0.0 || *out > cap || (amount > 0.0 && *out == 0.0) {
                return Err(SolubilityCaloricErrorV1::Invalid(
                    "dissolved inventory outside original thermal concentration support",
                ));
            }
        }
        Ok(concentrations)
    }

    fn calibration(&self, temperature: Temperature, index: usize) -> Result<f64> {
        let law = self.caloric.components()[index].law();
        let saturation = law.saturation_concentration(temperature)?.value();
        let enthalpy = law.dissolution_enthalpy_at(temperature)?;
        let diagonal = index * self.caloric.components().len() + index;
        let mut coefficient = Sum::default();
        coefficient.add(product(self.interaction(temperature.value(), diagonal)?, enthalpy)?);
        coefficient.add(product(-R, self.inverse_temperature_m3_k_mol[diagonal])?);
        product(saturation, sum(coefficient)?)
    }

    fn calibration_change(
        &self,
        from: Temperature,
        to: Temperature,
        index: usize,
    ) -> Result<f64> {
        if from == to {
            return Ok(0.0);
        }
        let law = self.caloric.components()[index].law();
        let s0 = law.saturation_concentration(from)?.value();
        let s1 = law.saturation_concentration(to)?.value();
        let logarithm = law.log_saturation_ratio(from, to)?;
        let ds_direct = s0 * expm1(logarithm);
        let ds = if ds_direct.is_finite() {
            ds_direct
        } else {
            s1 - s0
        };
        let h0 = law.dissolution_enthalpy_at(from)?;
        let dh = law.dissolution_enthalpy_change_j_mol(from, to)?;
        let diagonal = index * self.caloric.components().len() + index;
        let c = self.inverse_temperature_m3_k_mol[diagonal];
        let b0 = self.interaction(from.value(), diagonal)?;
        let b1 = self.interaction(to.value(), diagonal)?;
        let db = product(c, inverse_difference(from.value(), to.value())?)?;
        // Exact telescoping identity avoids subtracting two large standards.
        let mut change = Sum::default();
        change.add(product3(ds, b0, h0)?);
        change.add(product3(-R, c, ds)?);
        change.add(product3(s1, db, h0)?);
        change.add(product3(s1, b1, dh)?);
        sum(change)
    }

    fn admit_caloric_domain(&self) -> Result<()> {
        // Fixed bounded source preparation, never repeated in thermal trials.
        const INTERVALS: usize = 32;
        let lo = self.temperature_domain[0].value();
        let hi = self.temperature_domain[1].value();
        let count = self.caloric.components().len();
        for (i, component) in self.caloric.components().iter().enumerate() {
            let diagonal = i * count + i;
            for part in 0..INTERVALS {
                let a = if part == 0 { lo } else { lo + (hi - lo) * part as f64 / INTERVALS as f64 };
                let b = if part + 1 == INTERVALS { hi } else { lo + (hi - lo) * (part + 1) as f64 / INTERVALS as f64 };
                if a == b {
                    continue;
                }
                let domain = [Temperature::new(a), Temperature::new(b)];
                let (cp, dh, h, saturation) = component_bounds(component, domain)?;
                if saturation.1 > self.reference.maximum_concentrations_mol_m3()[i] {
                    return Err(SolubilityCaloricErrorV1::Invalid(
                        "pure-pair saturation leaves original thermal concentration support",
                    ));
                }
                let b0 = self.interaction(a, diagonal)?;
                let b1 = self.interaction(b, diagonal)?;
                let interaction = Span(b0.min(b1), b0.max(b1));
                let inv_t2 = Span::new((1.0 / b / b).next_down(), (1.0 / a / a).next_up())?;
                let h_squared = h.square();
                let first = interaction.mul(h_squared.mul(inv_t2).scale(1.0 / R).add(dh));
                let second = h.mul(inv_t2).scale(-2.0 * self.inverse_temperature_m3_k_mol[diagonal]);
                let correction = saturation.mul(first.add(second));
                let effective = cp.add(correction);
                if !effective.0.is_finite() || !effective.1.is_finite() || effective.0 <= 0.0 {
                    return Err(SolubilityCaloricErrorV1::Invalid(
                        "positive effective dissolved heat capacity cannot be certified over thermal source",
                    ));
                }
            }
        }
        Ok(())
    }
}

impl SolubilityThermalInteractionV1 for VirialSolubilityThermalModel {
    fn temperature_domain(&self) -> [Temperature; 2] {
        self.temperature_domain
    }

    fn admit_binding(
        &self,
        caloric: &SolubilityCaloricModelV1,
        kinetics: &[SolubilityKineticsV1],
    ) -> Result<()> {
        if caloric != &self.caloric
            || kinetics.len() != self.caloric.components().len()
            || kinetics.iter().zip(self.caloric.components()).any(|(kinetic, component)| {
                kinetic.law() != component.law()
            })
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal virial source must retain exact original calorics and ordered phase laws",
            ));
        }
        // Complete-domain capacity and convexity were established once in new().
        Ok(())
    }

    fn conditional_saturation(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        phases: &[SolubilityKineticStateV1],
        component: usize,
    ) -> Result<AmountConcentration> {
        let concentration = self.concentrations(temperature, solvent_volume, phases)?;
        let count = self.caloric.components().len();
        if component >= count {
            return Err(SolubilityCaloricErrorV1::Invalid("thermal saturation component layout"));
        }
        let pure = self.caloric.components()[component]
            .law()
            .saturation_concentration(temperature)?
            .value();
        if pure > self.reference.maximum_concentrations_mol_m3()[component] {
            return Err(SolubilityCaloricErrorV1::Invalid("pure saturation outside thermal source support"));
        }
        let mut other = Sum::default();
        for (j, &c) in concentration[..count].iter().enumerate() {
            if j != component {
                other.add(product(self.interaction(temperature.value(), component * count + j)?, c)?);
            }
        }
        conditional_root(
            pure,
            self.interaction(temperature.value(), component * count + component)?,
            sum(other)?,
            self.reference.maximum_concentrations_mol_m3()[component],
        )
        .map(AmountConcentration::new)
        .map_err(|_| SolubilityCaloricErrorV1::Invalid("conditional thermal saturation outside representable source support"))
    }

    fn enthalpy_increment_j(
        &self,
        solvent_volume: Volume,
        from_temperature: Temperature,
        from_phases: &[SolubilityKineticStateV1],
        to_temperature: Temperature,
        to_phases: &[SolubilityKineticStateV1],
    ) -> Result<f64> {
        self.concentrations(from_temperature, solvent_volume, from_phases)?;
        self.concentrations(to_temperature, solvent_volume, to_phases)?;
        let count = self.caloric.components().len();
        let mut delta = [0.0; MAX_OSMOTIC_SOLUTES];
        let mut mean = [0.0; MAX_OSMOTIC_SOLUTES];
        for (i, (from, to)) in from_phases.iter().zip(to_phases).enumerate() {
            if from.total_solute() != to.total_solute() {
                return Err(SolubilityCaloricErrorV1::Invalid(
                    "thermal interaction increment requires conserved component totals",
                ));
            }
            // Preserve a resolvable trace exchange when the larger phase rounds
            // to the same f64 stock at both ends.
            delta[i] = if from.dissolved().value().max(to.dissolved().value())
                <= from.solid().value().max(to.solid().value())
            {
                to.dissolved().value() - from.dissolved().value()
            } else {
                from.solid().value() - to.solid().value()
            };
            mean[i] = from.dissolved().value() + 0.5 * delta[i];
        }
        let mut increment = Sum::default();
        for i in 0..count {
            let mut row = Sum::default();
            for (j, &mean_amount) in mean[..count].iter().enumerate() {
                row.add(product(self.inverse_temperature_m3_k_mol[i * count + j], mean_amount / solvent_volume.value())?);
            }
            increment.add(product3(R, delta[i], sum(row)?)?);
            let q0 = self.calibration(from_temperature, i)?;
            let dq = self.calibration_change(from_temperature, to_temperature, i)?;
            increment.add(product(delta[i], q0 + 0.5 * dq)?);
            increment.add(product(mean[i], dq)?);
        }
        sum(increment)
    }
}

fn component_bounds(
    component: &SolubilityCaloricComponentV1,
    domain: [Temperature; 2],
) -> Result<(Span, Span, Span, Span)> {
    let aq = component.dissolved_heat_capacity().heat_capacity_bounds(domain)?;
    let solid = component.solid_heat_capacity().heat_capacity_bounds(domain)?;
    let cp = Span::new(aq[0].value(), aq[1].value())?;
    let dcp = if component.common_phase_heat_capacity().is_some() {
        Span(0.0, 0.0)
    } else {
        Span::new((aq[0].value() - solid[1].value()).next_down(), (aq[1].value() - solid[0].value()).next_up())?
    };
    let law = component.law();
    let h0 = law.dissolution_enthalpy_at(domain[0])?;
    let h = Span(h0, h0).add(dcp.mul(Span(0.0, domain[1].value() - domain[0].value())));
    let s0 = law.saturation_concentration(domain[0])?.value();
    let s1 = law.saturation_concentration(domain[1])?.value();
    let saturation = if h.0 >= 0.0 || h.1 <= 0.0 {
        Span(s0.min(s1), s0.max(s1))
    } else {
        let inverse_span = -inverse_difference(domain[0].value(), domain[1].value())?;
        let logs = h.mul(Span(0.0, inverse_span)).scale(1.0 / R);
        Span::new((s0 * exp(logs.0)).next_down().min(s0.min(s1)), (s0 * exp(logs.1)).next_up().max(s0.max(s1)))?
    };
    Ok((cp, dcp, h, saturation))
}

// Tiny local interval arithmetic for bounded SOURCE admission. It has no
// runtime solver role. Outward rounding avoids treating sampled Cp as a proof.
#[derive(Clone, Copy)]
struct Span(f64, f64);
impl Span {
    fn new(low: f64, high: f64) -> Result<Self> {
        if low.is_finite() && high.is_finite() && low <= high {
            Ok(Self(low, high))
        } else {
            Err(SolubilityCaloricErrorV1::Arithmetic("thermal source interval is not representable"))
        }
    }
    fn add(self, rhs: Self) -> Self {
        if rhs.0 == 0.0 && rhs.1 == 0.0 { return self; }
        if self.0 == 0.0 && self.1 == 0.0 { return rhs; }
        Self((self.0 + rhs.0).next_down(), (self.1 + rhs.1).next_up())
    }
    fn mul(self, rhs: Self) -> Self {
        if (self.0 == 0.0 && self.1 == 0.0) || (rhs.0 == 0.0 && rhs.1 == 0.0) { return Self(0.0, 0.0); }
        let values = [self.0 * rhs.0, self.0 * rhs.1, self.1 * rhs.0, self.1 * rhs.1];
        if values.iter().any(|x| !x.is_finite()) { return Self(f64::NEG_INFINITY, f64::INFINITY); }
        Self(values.iter().copied().fold(f64::INFINITY, f64::min).next_down(),
             values.iter().copied().fold(f64::NEG_INFINITY, f64::max).next_up())
    }
    fn scale(self, factor: f64) -> Self { self.mul(Self(factor, factor)) }
    fn square(self) -> Self {
        if self.0 == 0.0 && self.1 == 0.0 { return self; }
        let lo = if self.0 <= 0.0 && self.1 >= 0.0 { 0.0 } else { (self.0.abs().min(self.1.abs()) * self.0.abs().min(self.1.abs())).next_down() };
        Self(lo, (self.0.abs().max(self.1.abs()) * self.0.abs().max(self.1.abs())).next_up())
    }
}

fn inverse_difference(from: f64, to: f64) -> Result<f64> {
    let stable = ((from - to) / from) / to;
    finite(if stable.is_finite() { stable } else { 1.0 / to - 1.0 / from })
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else {
        Err(SolubilityCaloricErrorV1::Arithmetic("thermal virial quantity is not representable"))
    }
}
fn product(a: f64, b: f64) -> Result<f64> { finite(a * b) }
fn product3(a: f64, b: f64, c: f64) -> Result<f64> {
    let ordinary = (a * b) * c;
    if ordinary.is_finite() && (ordinary != 0.0 || a == 0.0 || b == 0.0 || c == 0.0) {
        return Ok(ordinary);
    }
    for alternative in [(a * c) * b, (b * c) * a] {
        if alternative.is_finite() && alternative != 0.0 { return Ok(alternative); }
    }
    finite(ordinary)
}
fn sum(value: Sum) -> Result<f64> {
    value.finish().map_err(|_| SolubilityCaloricErrorV1::Arithmetic("thermal virial sum is not representable"))
}

#[cfg(test)]
mod tests;
