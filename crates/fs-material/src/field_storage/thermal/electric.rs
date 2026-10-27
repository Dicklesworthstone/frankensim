//! Source-pinned affine dielectric thermodynamics on fixed material axes.
//!
//! At fixed displacement D, psi=V D.epsilon(T)^-1.D/2. The same elementary
//! Helmholtz derivative used by C(T), L(T) and spatial magnetics is applied to
//! each actual principal component: S=-partial_T psi and U=psi+T*S. The
//! nonpositive fixed-D field heat capacity must be combined with the original
//! caloric chart. An epsilon(T) snapshot alone does not define thermal work.

use super::{
    FieldStorageError, FieldThermalState, affine_path, capacity_bound, energy, finite, positive,
    thermal_state,
};
use crate::{
    field_storage::planar::PlanarDielectricLaw,
    state_point::ResolvedMaterialStatePoint,
};
use fs_matdb::MaterialCard;

/// One complete material bundle, fixed dielectric axes and supported affine
/// principal permittivities. Other source coordinates and selections stay fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalPlanarDielectricMaterial {
    card: MaterialCard,
    material: ResolvedMaterialStatePoint,
    law: PlanarDielectricLaw,
    temperature_axis: String,
    temperature_range_k: [f64; 2],
    reference_temperature_k: f64,
    slopes_f_m_k: [f64; 2],
    minimum_f_m: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalPlanarDielectricMaterial {
    material: ResolvedMaterialStatePoint,
    law: PlanarDielectricLaw,
    temperature_k: f64,
    principal_f_m: [f64; 2],
    slopes_f_m_k: [f64; 2],
}

fn principal(law: PlanarDielectricLaw, material: &ResolvedMaterialStatePoint) -> [f64; 2] {
    let names = law.source_property_names();
    let first = material.property(names[0]).expect("admitted dielectric source").value_si();
    let second = material.property(*names.last().expect("dielectric law has a source"))
        .expect("admitted dielectric source").value_si();
    [first, second]
}

impl ThermalPlanarDielectricMaterial {
    /// Admit the exact selected source paths, including their interior knots.
    /// Scalar and orthotropic choices share the native dielectric admission;
    /// no normal response, mixture, dispersion or rotating axes are inferred.
    pub fn new(
        card: &MaterialCard,
        material: &ResolvedMaterialStatePoint,
        law: PlanarDielectricLaw,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        law.resolve(material)?;
        let names = law.source_property_names();
        let mut slopes_f_m_k = [0.0; 2];
        let mut minimum_f_m = [0.0; 2];
        let mut reference_temperature_k = 0.0;
        for axis in 0..2 {
            let name = names[axis.min(names.len() - 1)];
            let (reference, slope, minimum) = affine_path(
                card, material, name, temperature_axis, temperature_range_k,
            )?;
            reference_temperature_k = reference;
            slopes_f_m_k[axis] = slope;
            minimum_f_m[axis] = minimum;
        }
        Ok(Self {
            card: card.clone(), material: material.clone(), law,
            temperature_axis: temperature_axis.to_owned(), temperature_range_k,
            reference_temperature_k, slopes_f_m_k, minimum_f_m,
        })
    }
    pub fn card(&self) -> &MaterialCard { &self.card }
    pub fn material(&self) -> &ResolvedMaterialStatePoint { &self.material }
    pub fn law(&self) -> PlanarDielectricLaw { self.law }
    pub fn temperature_axis(&self) -> &str { &self.temperature_axis }
    pub fn temperature_range_k(&self) -> [f64; 2] { self.temperature_range_k }
    pub fn principal_slopes_f_m_k(&self) -> [f64; 2] { self.slopes_f_m_k }

    /// Re-evaluate the entire original bundle, retaining all nonthermal axes,
    /// source selections and validity gates, not just its permittivity entries.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedThermalPlanarDielectricMaterial, FieldStorageError> {
        if !temperature_k.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&temperature_k)
        {
            return Err(FieldStorageError::Invalid("dielectric temperature is outside its supported source interval"));
        }
        let query = self.material.query().with(&self.temperature_axis, temperature_k)
            .map_err(|_| FieldStorageError::Invalid("dielectric thermal query is unrepresentable"))?;
        let material = self.material.at_query_point(&self.card, &query)?;
        self.law.resolve(&material)?;
        let values = principal(self.law, &material);
        let original = principal(self.law, &self.material);
        for axis in 0..2 {
            let expected = self.slopes_f_m_k[axis].mul_add(
                temperature_k - self.reference_temperature_k, original[axis],
            );
            if !expected.is_finite()
                || (expected - values[axis]).abs()
                    > 512.0 * f64::EPSILON * expected.abs().max(values[axis].abs())
            {
                return Err(FieldStorageError::Invalid("dielectric source left its admitted affine law"));
            }
        }
        Ok(ResolvedThermalPlanarDielectricMaterial {
            material, law: self.law, temperature_k, principal_f_m: values,
            slopes_f_m_k: self.slopes_f_m_k,
        })
    }

    /// Match a native field refresh to this same derivative law and full source
    /// bundle. A changed phase/composition/pressure or claim is not a T refresh.
    pub fn at_material(
        &self,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<ResolvedThermalPlanarDielectricMaterial, FieldStorageError> {
        // Re-evaluation pins the selected claims, so its selection-policy
        // receipt can differ from the original single-claim admission. Retain
        // the physical/source identity rather than requiring receipt-byte equality.
        super::super::admit_source_temperature_refresh(&self.material, material)?;
        let temperature_k = material.query_point().iter()
            .find_map(|(name, value)| (name == &self.temperature_axis).then_some(*value))
            .ok_or(FieldStorageError::Invalid("dielectric temperature coordinate is missing"))?;
        let resolved = self.at_temperature(temperature_k)?;
        if resolved.material.properties().iter().zip(material.properties())
            .any(|(a,b)| a.value_si().to_bits() != b.value_si().to_bits()) {
            return Err(FieldStorageError::Invalid("native dielectric and thermal derivative must retain the same complete original source bundle"));
        }
        Ok(resolved)
    }

    /// Conservative fixed-|D| capacity bound throughout the complete source
    /// interval. The worst principal contribution bounds their quadratic sum;
    /// the thermal owner adds ALL triangles sharing each caloric receiver.
    pub fn minimum_heat_capacity_j_k(
        &self,
        maximum_displacement_c_m2: f64,
        volume_m3: f64,
    ) -> Result<f64, FieldStorageError> {
        positive(volume_m3, "dielectric volume must be finite positive")?;
        if !maximum_displacement_c_m2.is_finite() || maximum_displacement_c_m2 < 0.0 {
            return Err(FieldStorageError::Invalid("dielectric displacement envelope must be finite nonnegative"));
        }
        let mut minimum: f64 = 0.0;
        for axis in 0..2 {
            let worst = finite(energy(maximum_displacement_c_m2, self.minimum_f_m[axis])? * volume_m3,
                "dielectric capacity energy bound is unrepresentable")?;
            if maximum_displacement_c_m2 != 0.0 && worst == 0.0 {
                return Err(FieldStorageError::Invalid("dielectric capacity energy bound underflows"));
            }
            let ratio = finite(self.slopes_f_m_k[axis] / self.minimum_f_m[axis],
                "dielectric capacity derivative is unrepresentable")?;
            if self.slopes_f_m_k[axis] != 0.0 && ratio == 0.0 {
                return Err(FieldStorageError::Invalid("dielectric capacity derivative underflows"));
            }
            minimum = minimum.min(capacity_bound(worst, self.temperature_range_k[1], ratio)?);
        }
        Ok(minimum)
    }
}

impl ResolvedThermalPlanarDielectricMaterial {
    pub fn material(&self) -> &ResolvedMaterialStatePoint { &self.material }
    pub fn temperature_k(&self) -> f64 { self.temperature_k }

    /// Integrated Helmholtz, entropy, internal energy and fixed-D heat capacity
    /// on the actual mesh-coordinate displacement. Signed field contributions
    /// use the zero-field caloric chart as reference, not another heat store.
    pub fn state(
        &self,
        displacement_c_m2: [f64; 2],
        volume_m3: f64,
    ) -> Result<FieldThermalState, FieldStorageError> {
        positive(volume_m3, "dielectric volume must be finite positive")?;
        if displacement_c_m2.iter().any(|d| !d.is_finite()) {
            return Err(FieldStorageError::Invalid("dielectric displacement must be finite"));
        }
        let principal_d = match self.law.frame() {
            None => displacement_c_m2,
            Some(frame) => {
                let q = crate::tensor::rotation([0.0, 0.0, 1.0], frame.angle_rad());
                [q[0][0].mul_add(displacement_c_m2[0], q[1][0] * displacement_c_m2[1]),
                 q[0][1].mul_add(displacement_c_m2[0], q[1][1] * displacement_c_m2[1])]
            }
        };
        let mut result = FieldThermalState {
            free_energy_j: 0.0, entropy_j_k: 0.0, internal_energy_j: 0.0, heat_capacity_j_k: 0.0,
        };
        for (axis, d) in principal_d.into_iter().enumerate() {
            let psi = finite(energy(d, self.principal_f_m[axis])? * volume_m3,
                "dielectric free energy is unrepresentable")?;
            if d != 0.0 && psi == 0.0 {
                return Err(FieldStorageError::Invalid("dielectric free energy underflows"));
            }
            let ratio = finite(self.slopes_f_m_k[axis] / self.principal_f_m[axis],
                "dielectric derivative is unrepresentable")?;
            if self.slopes_f_m_k[axis] != 0.0 && ratio == 0.0 {
                return Err(FieldStorageError::Invalid("dielectric derivative underflows"));
            }
            let value = thermal_state(psi, self.temperature_k, ratio)?;
            result.free_energy_j = finite(result.free_energy_j + value.free_energy_j, "dielectric free energy is unrepresentable")?;
            result.entropy_j_k = finite(result.entropy_j_k + value.entropy_j_k, "dielectric entropy is unrepresentable")?;
            result.internal_energy_j = finite(result.internal_energy_j + value.internal_energy_j, "dielectric internal energy is unrepresentable")?;
            result.heat_capacity_j_k = finite(result.heat_capacity_j_k + value.heat_capacity_j_k, "dielectric heat capacity is unrepresentable")?;
        }
        Ok(result)
    }
}
