//! Original bounded material sources in explicitly declared tensor directions.
//!
//! This borrows the prepared source owners. Querying a tensor neither selects
//! a new claim nor rebuilds a conductivity table.

use super::{
    ConductivityTensor, TensorConductivityLaw, invalid, proper_rotation, rotated_diagonal,
};
use crate::conjugate::{ChtError, poll};
use fs_exec::CancelGate;
use fs_matdb::PropertyLineSlice;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::thermal::conduction::source::ResolvedConductionLayer;

/// Three original directional material laws in one declared proper frame.
///
/// Columns of the rotation are the material principal directions in grid
/// coordinates. The caller declares that these scalar source keys describe
/// those directions; no isotropy or orientation is inferred from their names.
/// All three sources must retain the same original caloric curve, pure phase,
/// temperature-axis convention and fixed pressure/composition/etc. context.
///
/// The common temperature and energy support is the INTERSECTION of the
/// original source domains. Point and interval queries never extrapolate.
/// Every crossed original conductivity knot contributes to an interval mean.
/// Preparation and evaluation preserve the borrowed source lines and receipts.
///
/// With enthalpy-coupling, this same value implements TensorCaloricLaw. Supply
/// it to both the conductivity and storage fields to preserve the restricted
/// original energy domain, including the solid-side tangent at melting onset.
/// This is a pure-phase source view, not a new phase-mixture law, a density-driven
/// geometry update, or a newly identified thermodynamic chart.
#[derive(Clone, Copy, Debug)]
pub struct RotatedSourceConductivity<'a, 'm> {
    layers: [&'a ResolvedConductionLayer<'m>; 3],
    rotation: [[f64; 3]; 3],
    temperature_range_k: [f64; 2],
    specific_energy_range_j_kg: [f64; 2],
    energy_boundary_owner: [usize; 2],
}

impl<'a, 'm> RotatedSourceConductivity<'a, 'm> {
    /// Borrow three original principal-direction sources once.
    ///
    /// Cancellation is checked during fixed-context comparisons. No source
    /// interpolation, material selection, normalization or table copying occurs.
    pub fn new(
        layers: [&'a ResolvedConductionLayer<'m>; 3],
        rotation: [[f64; 3]; 3],
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        if !proper_rotation(rotation) {
            return Err(invalid(
                "source tensor directions require a proper orthogonal frame",
            ));
        }
        let mut temperature_range_k = layers[0].conductivity_curve().temperature_range_k();
        let mut specific_energy_range_j_kg = layers[0].enthalpy_range_j_kg();
        let mut energy_boundary_owner = [0; 2];
        for axis in 1..3 {
            poll(gate)?;
            if layers[axis].caloric_curve().identity() != layers[0].caloric_curve().identity()
                || layers[axis].phase() != layers[0].phase()
                || !same_context(layers[0].source_line(), layers[axis].source_line(), gate)?
            {
                return Err(invalid(
                    "directional sources must share the original caloric curve, phase and fixed context",
                ));
            }
            let temperatures = layers[axis].conductivity_curve().temperature_range_k();
            temperature_range_k[0] = temperature_range_k[0].max(temperatures[0]);
            temperature_range_k[1] = temperature_range_k[1].min(temperatures[1]);
            let energies = layers[axis].enthalpy_range_j_kg();
            if energies[0] > specific_energy_range_j_kg[0] {
                specific_energy_range_j_kg[0] = energies[0];
                energy_boundary_owner[0] = axis;
            }
            if energies[1] < specific_energy_range_j_kg[1] {
                specific_energy_range_j_kg[1] = energies[1];
                energy_boundary_owner[1] = axis;
            }
        }
        if temperature_range_k[0] >= temperature_range_k[1]
            || specific_energy_range_j_kg[0] >= specific_energy_range_j_kg[1]
        {
            return Err(invalid(
                "directional conductivity sources have no common finite interval",
            ));
        }
        poll(gate)?;
        Ok(Self {
            layers,
            rotation,
            temperature_range_k,
            specific_energy_range_j_kg,
            energy_boundary_owner,
        })
    }

    /// The original owners, in declared principal-direction order.
    pub const fn layers(&self) -> [&'a ResolvedConductionLayer<'m>; 3] {
        self.layers
    }

    /// The unchanged original thermodynamic chart; the source view is narrower.
    pub fn caloric_curve(&self) -> &'m EquilibriumEnthalpyPhaseCurve {
        self.layers[0].caloric_curve()
    }

    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }

    pub const fn specific_energy_range_j_kg(&self) -> [f64; 2] {
        self.specific_energy_range_j_kg
    }

    fn admit_temperature(&self, temperature_k: f64) -> Result<(), ChtError> {
        let [lo, hi] = self.temperature_range_k;
        if !temperature_k.is_finite() || !(lo..=hi).contains(&temperature_k) {
            return Err(invalid(
                "temperature leaves the original directional source intersection",
            ));
        }
        Ok(())
    }

    fn tensor(&self, principal: [f64; 3]) -> Result<ConductivityTensor, ChtError> {
        ConductivityTensor::new(
            rotated_diagonal(principal, self.rotation)
                .ok_or_else(|| invalid("original directional tensor is unrepresentable"))?,
        )
    }
}

fn original_error(error: impl std::fmt::Display) -> ChtError {
    ChtError::InvalidInput {
        field: "thermal.tensor_conduction.original_source",
        reason: error.to_string(),
    }
}

fn same_context(
    a: &PropertyLineSlice,
    b: &PropertyLineSlice,
    gate: &CancelGate,
) -> Result<bool, ChtError> {
    if a.axis != b.axis {
        return Ok(false);
    }
    let ar = &a
        .knots
        .first()
        .ok_or_else(|| invalid("empty original conductivity source"))?
        .1
        .receipt;
    let br = &b
        .knots
        .first()
        .ok_or_else(|| invalid("empty original conductivity source"))?
        .1
        .receipt;
    if ar.query_point.len() != br.query_point.len()
        || ar.axis_quantities.len() != br.axis_quantities.len()
    {
        return Ok(false);
    }
    for (index, ((an, av), (bn, bv))) in ar.query_point.iter().zip(&br.query_point).enumerate() {
        if index.is_multiple_of(256) {
            poll(gate)?;
        }
        if an != bn || (an != &a.axis && av.to_bits() != bv.to_bits()) {
            return Ok(false);
        }
    }
    for (index, (a, b)) in ar
        .axis_quantities
        .iter()
        .zip(&br.axis_quantities)
        .enumerate()
    {
        if index.is_multiple_of(256) {
            poll(gate)?;
        }
        if a != b {
            return Ok(false);
        }
    }
    Ok(true)
}

impl TensorConductivityLaw for RotatedSourceConductivity<'_, '_> {
    fn conductivity_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ConductivityTensor, ChtError> {
        self.admit_temperature(temperature_k)?;
        let mut principal = [0.0; 3];
        for (axis, value) in principal.iter_mut().enumerate() {
            *value = self.layers[axis]
                .conductivity_curve()
                .conductivity_w_m_k(temperature_k)
                .map_err(original_error)?;
        }
        self.tensor(principal)
    }

    fn mean_conductivity(&self, from_k: f64, to_k: f64) -> Result<ConductivityTensor, ChtError> {
        self.admit_temperature(from_k)?;
        self.admit_temperature(to_k)?;
        let mut principal = [0.0; 3];
        for (axis, value) in principal.iter_mut().enumerate() {
            *value = self.layers[axis]
                .conductivity_curve()
                .mean_conductivity_w_m_k(from_k, to_k)
                .map_err(original_error)?;
        }
        self.tensor(principal)
    }
}

#[cfg(feature = "enthalpy-coupling")]
mod caloric {
    use super::super::super::TensorCaloricLaw;
    use super::*;
    use fs_conduction::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};

    impl RotatedSourceConductivity<'_, '_> {
        fn energy_owner(
            &self,
            energy: f64,
        ) -> Result<&ResolvedConductionLayer<'_>, DepositionError> {
            let [lo, hi] = self.specific_energy_range_j_kg;
            if !energy.is_finite() || !(lo..=hi).contains(&energy) {
                return Err(DepositionError::MaterialDomain);
            }
            // Use the original owner of a binding endpoint. It retains exact
            // source-temperature pinning and its admitted one-sided tangent.
            let owner = if energy == lo {
                self.energy_boundary_owner[0]
            } else if energy == hi {
                self.energy_boundary_owner[1]
            } else {
                0
            };
            Ok(self.layers[owner])
        }
    }

    impl ThermalLaw for RotatedSourceConductivity<'_, '_> {
        fn storage_basis(&self) -> ThermalStorageBasis {
            self.caloric_curve().storage_basis()
        }

        fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
            self.energy_owner(energy)?
                .temperature_at_enthalpy(energy)
                .map_err(|_| DepositionError::MaterialDomain)
        }

        fn bulk_density_kg_m3(&self, energy: f64) -> Result<Option<f64>, DepositionError> {
            self.energy_owner(energy)?;
            self.caloric_curve().bulk_density_kg_m3(energy)
        }

        fn maximum_temperature_slope(&self) -> f64 {
            // A bound on the full original chart also bounds the restricted view.
            self.caloric_curve().maximum_temperature_slope()
        }

        fn specific_entropy_change_j_kg_k(
            &self,
            initial: f64,
            final_energy: f64,
        ) -> Result<Option<f64>, DepositionError> {
            self.energy_owner(initial)?;
            self.energy_owner(final_energy)?;
            ThermalLaw::specific_entropy_change_j_kg_k(self.caloric_curve(), initial, final_energy)
        }
    }

    impl TensorCaloricLaw for RotatedSourceConductivity<'_, '_> {
        fn specific_energy_bounds_j_kg(&self) -> [f64; 2] {
            self.specific_energy_range_j_kg
        }

        fn temperature_slope_kg_k_j(&self, energy: f64) -> Result<f64, DepositionError> {
            Ok(self.temperature_and_slope_at_specific_energy(energy)?[1])
        }

        fn temperature_and_slope_at_specific_energy(
            &self,
            energy: f64,
        ) -> Result<[f64; 2], DepositionError> {
            let state = self
                .energy_owner(energy)?
                .conduction_state_at_enthalpy(energy)?;
            Ok([state.temperature_k, state.temperature_slope_kg_k_j])
        }
    }
}

#[cfg(test)]
mod tests;
