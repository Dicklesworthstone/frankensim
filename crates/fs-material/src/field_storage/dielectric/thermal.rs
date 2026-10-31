//! Temperature dependence of the original Debye polarization state.
//!
//! Each independent charge stores Helmholtz energy psi=q²/(2 C(T)). For an
//! admitted affine C(T), S=psi C'/C and U=psi+T S; its fixed-charge contribution
//! to heat capacity is -2 T psi (C'/C)². Sum every instantaneous/polarization
//! contribution with the original zero-field caloric chart. Relaxation time
//! changes affect dissipation, not these potentials. The circuit remains the
//! owner of charge history and irreversible polarization heat.
//!
//! This reuses the existing field thermodynamic potential and exact source-path
//! admission. It supplies no new heat store, time integrator or fitted pole data.

use fs_matdb::MaterialCard;

use super::{ResolvedDebyeCapacitor, finite, invalid};
use crate::{
    field_storage::{
        FieldStorageError, admit_source_temperature_refresh, energy,
        thermal::{FieldThermalState, affine_source_path, capacity_bound, thermal_state},
    },
    state_point::{ResolvedMaterialStatePoint, pinned_temperature_line},
};

/// One original material and fixed geometry over an explicit temperature range.
/// Every electric storage coefficient must have one positive affine source path.
/// Relaxation times may use positive continuous piecewise-linear paths, since
/// their derivatives do not enter stored energy or entropy.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalDebyeCapacitor {
    card: MaterialCard,
    original: ResolvedDebyeCapacitor,
    temperature_axis: String,
    temperature_range_k: [f64; 2],
    original_temperature_k: f64,
    slopes_f_k: Vec<f64>,
    minimum_capacitances_f: Vec<f64>,
}

/// The complete refreshed source and the same admitted potential derivatives.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalDebyeCapacitor {
    source: ResolvedDebyeCapacitor,
    temperature_k: f64,
    slopes_f_k: Vec<f64>,
}

fn coefficient(source: &ResolvedDebyeCapacitor, index: usize) -> f64 {
    if index == 0 {
        source.instantaneous_capacitance_f()
    } else {
        source.poles()[index - 1].capacitance_f
    }
}

fn add(sum: &mut f64, correction: &mut f64, term: f64) -> Result<(), FieldStorageError> {
    let next = finite(*sum + term)?;
    *correction = finite(
        *correction
            + if sum.abs() >= term.abs() {
                (*sum - next) + term
            } else {
                (term - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}

impl ThermalDebyeCapacitor {
    /// Admit all original epsilon/relaxation paths, including interior knots,
    /// while retaining every other property in the complete source bundle.
    /// The temperature axis needs explicit absolute-temperature semantics.
    pub fn new(
        card: &MaterialCard,
        original: &ResolvedDebyeCapacitor,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        let material = original.material();
        let scale = original.area_m2() / original.gap_m();
        let mut slopes_f_k = Vec::with_capacity(original.poles().len() + 1);
        let mut minimum_capacitances_f = Vec::with_capacity(original.poles().len() + 1);
        let mut original_temperature_k = 0.0;
        for name in std::iter::once(original.high_frequency_property()).chain(
            original
                .pole_properties()
                .iter()
                .map(|p| p.permittivity_increment.as_str()),
        ) {
            let (reference, slope, endpoints) = affine_source_path(
                card,
                material,
                name,
                temperature_axis,
                temperature_range_k,
                true,
            )?;
            original_temperature_k = reference;
            let slope_f_k = finite(slope * scale)?;
            if slope != 0.0 && slope_f_k == 0.0 {
                return Err(invalid(
                    "dielectric capacitance temperature derivative underflows",
                ));
            }
            slopes_f_k.push(slope_f_k);
            minimum_capacitances_f.push(finite(endpoints[0].min(endpoints[1]) * scale)?);
        }
        let point = |t| {
            material
                .query()
                .with(temperature_axis, t)
                .map_err(|_| invalid("dielectric temperature query is unrepresentable"))
        };
        let lower = point(temperature_range_k[0])?;
        let upper = point(temperature_range_k[1])?;
        // These endpoint queries validate all original properties, not only
        // the coefficients used in the thermal potential.
        original.at_query_point(card, &lower)?;
        original.at_query_point(card, &upper)?;
        for pole in original.pole_properties() {
            let property = material.property(&pole.relaxation_time).ok_or(invalid(
                "dielectric relaxation source is missing or ambiguous",
            ))?;
            if !property.has_only_fixed_external_conditions(Some(temperature_axis), false) {
                return Err(invalid(
                    "dielectric relaxation requires one owned temperature and fixed external coordinates",
                ));
            }
            let selected = property.answer().receipt.selected;
            let claim = card.claims().claim(selected).ok_or(invalid(
                "dielectric relaxation claim is absent from its original card",
            ))?;
            let knots = pinned_temperature_line(
                card,
                &claim.key,
                temperature_axis,
                &lower,
                &upper,
                selected,
            )?;
            if knots.first().map(|v| v.0) != Some(temperature_range_k[0])
                || knots.last().map(|v| v.0) != Some(temperature_range_k[1])
                || knots
                    .iter()
                    .any(|(_, value)| !value.is_finite() || *value <= 0.0)
            {
                return Err(invalid(
                    "positive continuous relaxation times must cover the full dielectric temperature range",
                ));
            }
        }
        Ok(Self {
            card: card.clone(),
            original: original.clone(),
            temperature_axis: temperature_axis.to_owned(),
            temperature_range_k,
            original_temperature_k,
            slopes_f_k,
            minimum_capacitances_f,
        })
    }

    pub fn card(&self) -> &MaterialCard {
        &self.card
    }
    pub fn original(&self) -> &ResolvedDebyeCapacitor {
        &self.original
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.original.material()
    }
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    pub fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }
    /// dC/dT in instantaneous-then-polarization order [F/K].
    pub fn capacitance_slopes_f_k(&self) -> &[f64] {
        &self.slopes_f_k
    }

    /// Refresh all pinned source claims, with every nonthermal coordinate fixed.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedThermalDebyeCapacitor, FieldStorageError> {
        if !temperature_k.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&temperature_k)
        {
            return Err(invalid(
                "dielectric temperature is outside its admitted source interval",
            ));
        }
        let point = self
            .material()
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|_| invalid("dielectric temperature query is unrepresentable"))?;
        let source = self.original.at_query_point(&self.card, &point)?;
        for (index, &slope) in self.slopes_f_k.iter().enumerate() {
            let expected = finite(slope.mul_add(
                temperature_k - self.original_temperature_k,
                coefficient(&self.original, index),
            ))?;
            let actual = coefficient(&source, index);
            if (actual - expected).abs() > 512.0 * f64::EPSILON * actual.abs().max(expected.abs()) {
                return Err(invalid(
                    "dielectric source left its admitted affine storage law",
                ));
            }
        }
        Ok(ResolvedThermalDebyeCapacitor {
            source,
            temperature_k,
            slopes_f_k: self.slopes_f_k.clone(),
        })
    }

    /// Match a native source refresh to this exact thermal derivative law.
    pub fn at_source(
        &self,
        source: &ResolvedDebyeCapacitor,
    ) -> Result<ResolvedThermalDebyeCapacitor, FieldStorageError> {
        admit_source_temperature_refresh(self.material(), source.material())?;
        if source.high_frequency_property() != self.original.high_frequency_property()
            || source.pole_properties() != self.original.pole_properties()
            || source.area_m2().to_bits() != self.original.area_m2().to_bits()
            || source.gap_m().to_bits() != self.original.gap_m().to_bits()
        {
            return Err(invalid(
                "native dielectric and thermal law must retain the same geometry and ordered pole roles",
            ));
        }
        let temperature = source
            .material()
            .state_coordinate(&self.temperature_axis)
            .ok_or(invalid("dielectric temperature coordinate is missing"))?;
        let resolved = self.at_temperature(temperature)?;
        if source
            .material()
            .properties()
            .iter()
            .zip(resolved.source.material().properties())
            .any(|(a, b)| a.value_si().to_bits() != b.value_si().to_bits())
        {
            return Err(invalid(
                "native dielectric and thermal law must retain the same complete material snapshot",
            ));
        }
        Ok(resolved)
    }

    /// Conservative lower bound on fixed-charge dielectric dU/dT [J/K] for the
    /// specified absolute charge envelopes over the complete source interval.
    /// A coupled owner must add this to the original caloric capacity, require a
    /// positive sum and enforce every envelope. This is a model bound with
    /// roundoff padding, not a time-discretization or experimental error bound.
    pub fn minimum_heat_capacity_j_k(
        &self,
        maximum_instantaneous_charge_c: f64,
        maximum_polarization_charges_c: &[f64],
    ) -> Result<f64, FieldStorageError> {
        if maximum_polarization_charges_c.len() != self.original.poles().len()
            || std::iter::once(&maximum_instantaneous_charge_c)
                .chain(maximum_polarization_charges_c)
                .any(|q| !q.is_finite() || *q < 0.0)
        {
            return Err(invalid(
                "complete finite nonnegative dielectric charge envelopes required",
            ));
        }
        let mut sum = 0.0;
        let mut correction = 0.0;
        for (index, charge) in std::iter::once(maximum_instantaneous_charge_c)
            .chain(maximum_polarization_charges_c.iter().copied())
            .enumerate()
        {
            let slope = self.slopes_f_k[index];
            if charge == 0.0 || slope == 0.0 {
                continue;
            }
            let minimum = self.minimum_capacitances_f[index];
            let ratio = finite(slope / minimum)?;
            if ratio == 0.0 {
                return Err(invalid("dielectric heat capacity derivative underflows"));
            }
            let bound =
                capacity_bound(energy(charge, minimum)?, self.temperature_range_k[1], ratio)?;
            add(&mut sum, &mut correction, bound)?;
        }
        finite((sum + correction) * (1.0 + 64.0 * f64::EPSILON))
    }
}

impl ResolvedThermalDebyeCapacitor {
    pub fn source(&self) -> &ResolvedDebyeCapacitor {
        &self.source
    }
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub fn capacitance_slopes_f_k(&self) -> &[f64] {
        &self.slopes_f_k
    }

    /// Cancellation-stable change in the thermal offset T*S between complete
    /// charge/temperature endpoints of this same affine source law. A coupled
    /// owner combines it with the separately resolved Helmholtz change; it is
    /// signed reversible exchange, not additional polarization dissipation.
    ///
    /// Split at the original charges and final temperature. The charge change
    /// uses its quadratic midpoint increment, including midpoint roundoff. At
    /// fixed charge S is proportional to 1/C², so its change follows the bounded
    /// relative capacitance difference without subtracting large entropies.
    /// In d(T*S), the lower temperature multiplies dS; this keeps both products
    /// bounded by the original endpoint offsets. All pole terms are compensated.
    pub fn thermal_offset_change_j(
        &self,
        after: &Self,
        before_instantaneous_charge_c: f64,
        before_polarization_charges_c: &[f64],
        after_instantaneous_charge_c: f64,
        after_polarization_charges_c: &[f64],
    ) -> Result<f64, FieldStorageError> {
        if self.slopes_f_k != after.slopes_f_k {
            return Err(invalid(
                "dielectric thermal offset requires the same admitted affine derivatives",
            ));
        }
        self.state(before_instantaneous_charge_c, before_polarization_charges_c)?;
        after.state(after_instantaneous_charge_c, after_polarization_charges_c)?;
        self.source.fixed_charge_free_energy_change_j(
            &after.source,
            before_instantaneous_charge_c,
            before_polarization_charges_c,
        )?;
        let temperature_change = finite(after.temperature_k - self.temperature_k)?;
        let mut sum = 0.0;
        let mut correction = 0.0;
        let before_charges = std::iter::once(before_instantaneous_charge_c)
            .chain(before_polarization_charges_c.iter().copied());
        let after_charges = std::iter::once(after_instantaneous_charge_c)
            .chain(after_polarization_charges_c.iter().copied());
        for (index, (q0, q1)) in before_charges.zip(after_charges).enumerate() {
            let slope = self.slopes_f_k[index];
            if slope == 0.0 {
                continue;
            }
            let c0 = coefficient(&self.source, index);
            let c1 = coefficient(&after.source, index);
            let r0 = finite(slope / c0)?;
            let r1 = finite(slope / c1)?;
            let initial = thermal_state(energy(q0, c0)?, self.temperature_k, r0)?;
            let bridge = thermal_state(energy(q0, c1)?, after.temperature_k, r1)?;
            let charge_change = finite(q1 - q0)?;
            let midpoint = f64::midpoint(q0, q1);
            let midpoint_offset = 0.5 * (q0 - midpoint) + 0.5 * (q1 - midpoint);
            for effort in [midpoint / c1, midpoint_offset / c1] {
                let energy_change = finite(charge_change * effort)?;
                let entropy_change = finite(energy_change * r1)?;
                let offset_change = finite(entropy_change * after.temperature_k)?;
                if energy_change != 0.0 && (entropy_change == 0.0 || offset_change == 0.0) {
                    return Err(invalid(
                        "nonzero dielectric charge thermal offset underflows",
                    ));
                }
                add(&mut sum, &mut correction, offset_change)?;
            }
            let largest_entropy = if c0 <= c1 {
                initial.entropy_j_k
            } else {
                bridge.entropy_j_k
            };
            let entropy_change = finite(
                largest_entropy * ((c0 - c1) / c0.max(c1)) * (1.0 + c0.min(c1) / c0.max(c1)),
            )?;
            if largest_entropy != 0.0 && c0 != c1 && entropy_change == 0.0 {
                return Err(invalid(
                    "nonzero dielectric parameter entropy change underflows",
                ));
            }
            let endpoint_entropy = if self.temperature_k <= after.temperature_k {
                bridge.entropy_j_k
            } else {
                initial.entropy_j_k
            };
            add(
                &mut sum,
                &mut correction,
                finite(endpoint_entropy * temperature_change)?,
            )?;
            add(
                &mut sum,
                &mut correction,
                finite(self.temperature_k.min(after.temperature_k) * entropy_change)?,
            )?;
        }
        finite(sum + correction)
    }

    /// Helmholtz energy, entropy, internal energy and fixed-charge heat capacity
    /// of ALL retained dielectric charges. The caloric owner supplies its own
    /// zero-field energy/entropy. This bounded evaluation does not allocate.
    pub fn state(
        &self,
        instantaneous_charge_c: f64,
        polarization_charges_c: &[f64],
    ) -> Result<FieldThermalState, FieldStorageError> {
        if polarization_charges_c.len() != self.source.poles().len() {
            return Err(invalid(
                "one original charge per dielectric polarization pole required",
            ));
        }
        let mut sums = [0.0; 4];
        let mut corrections = [0.0; 4];
        for (index, charge) in std::iter::once(instantaneous_charge_c)
            .chain(polarization_charges_c.iter().copied())
            .enumerate()
        {
            let c = coefficient(&self.source, index);
            let ratio = finite(self.slopes_f_k[index] / c)?;
            if self.slopes_f_k[index] != 0.0 && ratio == 0.0 {
                return Err(invalid(
                    "dielectric thermal potential derivative underflows",
                ));
            }
            let state = thermal_state(energy(charge, c)?, self.temperature_k, ratio)?;
            for (slot, value) in [
                state.free_energy_j,
                state.entropy_j_k,
                state.internal_energy_j,
                state.heat_capacity_j_k,
            ]
            .into_iter()
            .enumerate()
            {
                add(&mut sums[slot], &mut corrections[slot], value)?;
            }
        }
        for i in 0..4 {
            sums[i] = finite(sums[i] + corrections[i])?;
        }
        Ok(FieldThermalState {
            free_energy_j: sums[0],
            entropy_j_k: sums[1],
            internal_energy_j: sums[2],
            heat_capacity_j_k: sums[3],
        })
    }
}
