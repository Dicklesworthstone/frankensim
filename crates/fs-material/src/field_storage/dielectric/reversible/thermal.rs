//! Explicit separable dielectric thermodynamics: D(E,T)=s(T)D_ref(E).
//!
//! The original dimensionless scale is positive and affine, and s(T_ref)=1.
//! With z=D/s, psi=s w_ref(z), S=s' g_ref(z), U=psi+TS and
//! C_D=-T (s')² z²/(s epsilon_ref,diff). The negative fixed-charge field
//! capacity must be combined with the original zero-field caloric capacity.
//! General D(E,T) surfaces, latent heat, and dielectric losses are not inferred.

use super::*;
use crate::field_storage::{
    admit_source_temperature_refresh,
    thermal::{FieldThermalState, affine_source_path},
};
use fs_matdb::PropertyKey;
use fs_qty::Temperature;

#[derive(Debug, PartialEq)]
struct ScaleLaw {
    material: ResolvedMaterialStatePoint,
    property: String,
    temperature_axis: String,
    temperature_range_k: [f64; 2],
    reference_temperature_k: f64,
    slope_per_k: f64,
    minimum_scale: f64,
}

#[derive(Debug, PartialEq)]
pub(super) struct ResolvedScale {
    law: Arc<ScaleLaw>,
    pub material: ResolvedMaterialStatePoint,
    pub value: f64,
    temperature_k: f64,
}

/// Both original claims and the complete material bundle, pinned on one
/// supported affine-temperature strip. Clones share immutable prepared data.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalReversibleDielectricMaterial {
    reference: ReversibleDielectricMaterial,
    scale: Arc<ScaleLaw>,
}

/// One physical scaled D(E) snapshot and its source-exact thermal derivatives.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalReversibleDielectricMaterial {
    dielectric: ReversibleDielectricMaterial,
}

fn volume(volume_m3: f64) -> Result<f64, FieldStorageError> {
    if volume_m3.is_finite() && volume_m3 > 0.0 {
        Ok(volume_m3)
    } else {
        Err(invalid(
            "nonlinear dielectric volume must be finite and positive",
        ))
    }
}
fn nonzero(value: f64) -> Result<f64, FieldStorageError> {
    let value = finite(value)?;
    if value == 0.0 {
        Err(invalid(
            "nonzero nonlinear dielectric thermal contribution underflows",
        ))
    } else {
        Ok(value)
    }
}
fn add_product(
    sum: &mut f64,
    correction: &mut f64,
    left: f64,
    right: f64,
) -> Result<(), FieldStorageError> {
    let product = finite(left * right)?;
    if left != 0.0 && right != 0.0 && product == 0.0 {
        return Err(invalid("nonzero dielectric entropy work underflows"));
    }
    add(sum, correction, product)?;
    add(sum, correction, left.mul_add(right, -product))
}

// Prove the original field path still exists at the new applicability point.
// Its absent temperature value axis plus rectangular validity then proves the
// entire E/T strip from the two boundary paths, without a sampled surface fit.
fn admit_reference_path(
    reference: &ReversibleDielectricMaterial,
    material: &ResolvedMaterialStatePoint,
) -> Result<(), FieldStorageError> {
    let original = reference
        .material()
        .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
        .expect("prepared reversible source");
    let actual = material
        .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
        .ok_or(invalid(
            "thermal bundle lost its original dielectric source",
        ))?;
    if reference.material().card_identity() != material.card_identity()
        || reference.material().material() != material.material()
        || original.requirement() != actual.requirement()
        || original.source_axes() != actual.source_axes()
        || original.answer().receipt.selected != actual.answer().receipt.selected
    {
        return Err(invalid(
            "thermal dielectric changed its original material or source",
        ));
    }
    let query = material.query();
    let checked = reference
        .material()
        .at_query_point(reference.card(), &query)?;
    admit_source_temperature_refresh(reference.material(), &checked)?;
    let selected = original.answer().receipt.selected;
    let claim = reference
        .card()
        .claims()
        .claim(selected)
        .expect("pinned D source");
    let upper = query
        .clone()
        .with(
            &reference.source_path().axis,
            reference.maximum_field_strength_v_m(),
        )
        .map_err(|_| invalid("thermal dielectric field range is unrepresentable"))?;
    let path = reference
        .card()
        .claims()
        .query_line_typed(
            &claim.key,
            &reference.source_path().axis,
            &query,
            &upper,
            ClaimSelection::Pinned(selected),
        )
        .map_err(|source| MaterialStatePointError::Support {
            property: REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY.to_owned(),
            source: Box::new(source),
        })?;
    if path.knots.len() != reference.source_path().knots.len()
        || path
            .knots
            .iter()
            .zip(&reference.source_path().knots)
            .any(|((e, d), (e0, d0))| {
                e.to_bits() != e0.to_bits()
                    || d.evidence.value.value.to_bits() != d0.evidence.value.value.to_bits()
                    || d.receipt.selected != selected
            })
    {
        return Err(invalid(
            "thermal dielectric reference path changed its original D(E) law",
        ));
    }
    Ok(())
}

impl ThermalReversibleDielectricMaterial {
    /// Explicitly choose D=s(T)D_ref(E), retaining BOTH original source claims.
    /// The reference must have no temperature value axis; the separate positive
    /// dimensionless scale must be exactly one at the common reference point.
    /// Every original requirement, source pin and nonthermal coordinate remains
    /// fixed. No separability is inferred from a two-coordinate D(E,T) surface.
    pub fn new(
        reference: &ReversibleDielectricMaterial,
        material: &ResolvedMaterialStatePoint,
        scale_property: &PropertyKey,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        let source = reference
            .material()
            .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
            .expect("prepared reversible source");
        if reference.is_temperature_scaled()
            || source
                .source_axes()
                .iter()
                .any(|(axis, dims)| axis == temperature_axis || *dims == Temperature::DIMS)
            || scale_property.name() == REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY
            || scale_property.quantity() != QuantitySpec::dimensional(Dims([0; 6]))
            || reference.material().query() != material.query()
            || reference.material().card_identity() != material.card_identity()
            || reference.material().material() != material.material()
            || reference.material().properties().iter().any(|original| {
                material
                    .property(original.requirement().name())
                    .is_none_or(|actual| {
                        original.requirement() != actual.requirement()
                            || original.source_axes() != actual.source_axes()
                            || original.answer().receipt.selected
                                != actual.answer().receipt.selected
                            || original.value_si().to_bits() != actual.value_si().to_bits()
                    })
            })
        {
            return Err(invalid(
                "thermal scaling needs the original temperature-independent D(E), complete common context and distinct dimensionless scale",
            ));
        }
        let checked = material.at_query_point(reference.card(), &material.query())?;
        admit_source_temperature_refresh(material, &checked)?;
        if material
            .properties()
            .iter()
            .zip(checked.properties())
            .any(|(a, b)| a.value_si().to_bits() != b.value_si().to_bits())
        {
            return Err(invalid(
                "thermal dielectric bundle is not its original source evaluation",
            ));
        }
        let property = material.property(scale_property.name()).ok_or(invalid(
            "complete thermal dielectric bundle lacks its scale source",
        ))?;
        if property.requirement().quantity() != scale_property.quantity()
            || property.value_si() != 1.0
        {
            return Err(invalid(
                "reference dielectric scale must be dimensionless and exactly one",
            ));
        }
        let (reference_temperature_k, slope_per_k, endpoints) = affine_source_path(
            reference.card(),
            material,
            scale_property.name(),
            temperature_axis,
            temperature_range_k,
            true,
        )?;
        let law = Self {
            reference: reference.clone(),
            scale: Arc::new(ScaleLaw {
                material: material.clone(),
                property: scale_property.name().to_owned(),
                temperature_axis: temperature_axis.to_owned(),
                temperature_range_k,
                reference_temperature_k,
                slope_per_k,
                minimum_scale: endpoints[0].min(endpoints[1]),
            }),
        };
        for temperature in temperature_range_k {
            let snapshot = law.at_temperature(temperature)?;
            admit_reference_path(reference, snapshot.material())?;
        }
        Ok(law)
    }

    pub fn reference(&self) -> &ReversibleDielectricMaterial {
        &self.reference
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.scale.material
    }
    pub fn card(&self) -> &MaterialCard {
        self.reference.card()
    }
    pub fn temperature_axis(&self) -> &str {
        &self.scale.temperature_axis
    }
    pub fn temperature_range_k(&self) -> [f64; 2] {
        self.scale.temperature_range_k
    }
    pub fn scale_slope_per_k(&self) -> f64 {
        self.scale.slope_per_k
    }

    /// Requery every original source role at T, sharing the immutable D_ref
    /// table. Source receipts stay reference evidence; the physical snapshot
    /// evaluates the actual scaled D, tangent, inverse and both potentials.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedThermalReversibleDielectricMaterial, FieldStorageError> {
        if !temperature_k.is_finite()
            || !(self.scale.temperature_range_k[0]..=self.scale.temperature_range_k[1])
                .contains(&temperature_k)
        {
            return Err(invalid(
                "nonlinear dielectric temperature left its complete source interval",
            ));
        }
        let query = self
            .material()
            .query()
            .with(self.temperature_axis(), temperature_k)
            .map_err(|_| invalid("nonlinear dielectric thermal query is unrepresentable"))?;
        let material = self.material().at_query_point(self.card(), &query)?;
        admit_source_temperature_refresh(self.material(), &material)?;
        let value = material
            .property(&self.scale.property)
            .expect("pinned scale")
            .value_si();
        let expected = finite(
            self.scale
                .slope_per_k
                .mul_add(temperature_k - self.scale.reference_temperature_k, 1.0),
        )?;
        if !value.is_finite()
            || value <= 0.0
            || expected <= 0.0
            || (value - expected).abs() > 512.0 * f64::EPSILON * value.abs().max(expected.abs())
        {
            return Err(invalid(
                "dielectric scale left its positive affine source path",
            ));
        }
        let dielectric = ReversibleDielectricMaterial {
            prepared: self.reference.prepared.clone(),
            temperature_scale: Some(Arc::new(ResolvedScale {
                law: self.scale.clone(),
                material,
                value,
                temperature_k,
            })),
        };
        // Multiplying raw receipts is not a new source fit. Admission ensures
        // that every scaled interval and both potentials remain representable.
        for i in 0..dielectric.prepared.slopes_f_m.len() {
            let d0 = dielectric.displacement_knot(i);
            let d1 = dielectric.displacement_knot(i + 1);
            if !d0.is_finite()
                || !d1.is_finite()
                || d1 <= d0
                || !dielectric.slope(i).is_finite()
                || dielectric.slope(i) <= 0.0
                || !dielectric.slope(i).recip().is_finite()
                || dielectric.slope(i).recip() <= 0.0
                || nonzero(dielectric.energy_knot(i + 1)).is_err()
                || nonzero(dielectric.coenergy_knot(i + 1)).is_err()
            {
                return Err(invalid(
                    "scaled dielectric source support or energy is unrepresentable",
                ));
            }
        }
        Ok(ResolvedThermalReversibleDielectricMaterial { dielectric })
    }

    /// Conservative lower bound on fixed-D field capacity over the full
    /// temperature interval and 0<=|D|<=maximum. Each original inverse segment
    /// contributes its worst endpoint, including both sides of tangent jumps.
    /// A finite thermal owner must add its caloric capacity and require >0.
    pub fn minimum_heat_capacity_j_k(
        &self,
        maximum_displacement_c_m2: f64,
        volume_m3: f64,
    ) -> Result<f64, FieldStorageError> {
        let volume = volume(volume_m3)?;
        let s = self.scale.minimum_scale;
        if !maximum_displacement_c_m2.is_finite()
            || maximum_displacement_c_m2 < 0.0
            || maximum_displacement_c_m2 > s * self.reference.maximum_displacement_c_m2()
        {
            return Err(invalid(
                "dielectric charge envelope must fit the full thermal source strip",
            ));
        }
        let a = self.scale.slope_per_k;
        if a == 0.0 || maximum_displacement_c_m2 == 0.0 {
            return Ok(0.0);
        }
        let maximum = maximum_displacement_c_m2 / s;
        let mut magnitude = 0.0_f64;
        for (i, &epsilon) in self.reference.prepared.slopes_f_m.iter().enumerate() {
            if self.reference.displacement_knot(i) > maximum {
                break;
            }
            let z = maximum.min(self.reference.displacement_knot(i + 1));
            let rate = nonzero(a * z)?;
            let bound = nonzero(
                (self.scale.temperature_range_k[1] * rate) * (rate / (s * epsilon)) * volume,
            )?;
            magnitude = magnitude.max(bound);
        }
        // The analytic extrema are exact; outward roundoff allowance keeps
        // this an actual lower bound when state() groups products differently.
        Ok(-nonzero(magnitude * (1.0 + 32.0 * f64::EPSILON))?)
    }
}

impl ResolvedThermalReversibleDielectricMaterial {
    fn scale(&self) -> &ResolvedScale {
        self.dielectric
            .temperature_scale
            .as_ref()
            .expect("resolved thermal scale")
    }
    pub fn dielectric(&self) -> &ReversibleDielectricMaterial {
        &self.dielectric
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.dielectric.material()
    }
    pub fn card(&self) -> &MaterialCard {
        self.dielectric.card()
    }
    pub fn temperature_k(&self) -> f64 {
        self.scale().temperature_k
    }
    pub fn displacement_scale(&self) -> f64 {
        self.scale().value
    }
    pub fn scale_slope_per_k(&self) -> f64 {
        self.scale().law.slope_per_k
    }

    /// Same original reference law, scale claim, full material bundle and
    /// nonthermal context; only the admitted temperature may change. This is
    /// source admission, not permission to omit parameter work or field entropy.
    pub fn admit_source_refresh(&self, after: &Self) -> Result<(), FieldStorageError> {
        if self.dielectric.prepared != after.dielectric.prepared
            || self.scale().law != after.scale().law
        {
            return Err(invalid(
                "thermal dielectric refresh changed the reference or affine scale source",
            ));
        }
        admit_source_temperature_refresh(self.material(), after.material())
    }

    /// Field thermodynamics at a NONNEGATIVE actual D magnitude. The zero-field
    /// caloric chart owns the thermal datum; S and U may be signed contributions.
    /// These are reversible terms, never an additional Joule-heat source.
    pub fn state(
        &self,
        displacement_c_m2: f64,
        volume_m3: f64,
    ) -> Result<FieldThermalState, FieldStorageError> {
        let volume = volume(volume_m3)?;
        let response = self.dielectric.evaluate_displacement(displacement_c_m2)?;
        let free_energy_j = finite(volume * response.energy_density_j_m3)?;
        if displacement_c_m2 != 0.0 && free_energy_j == 0.0 {
            return Err(invalid(
                "nonzero nonlinear dielectric free energy underflows",
            ));
        }
        if self.scale_slope_per_k() == 0.0 || displacement_c_m2 == 0.0 {
            return Ok(FieldThermalState {
                free_energy_j,
                entropy_j_k: 0.0,
                internal_energy_j: free_energy_j,
                heat_capacity_j_k: 0.0,
            });
        }
        let ratio = nonzero(self.scale_slope_per_k() / self.displacement_scale())?;
        let entropy_j_k = nonzero((ratio * response.coenergy_density_j_m3) * volume)?;
        let internal_energy_j = finite(self.temperature_k().mul_add(entropy_j_k, free_energy_j))?;
        let rate = nonzero(ratio * displacement_c_m2)?;
        let heat_capacity_j_k = -nonzero(
            (self.temperature_k() * rate)
                * (rate / response.differential_permittivity_f_m)
                * volume,
        )?;
        Ok(FieldThermalState {
            free_energy_j,
            entropy_j_k,
            internal_energy_j,
            heat_capacity_j_k,
        })
    }

    /// Stable change of T*S [J] for SIGNED scalar D endpoints (entropy is even).
    /// Original inverse-source segments integrate the coenergy change directly.
    /// The affine scale supplies ds=s' dT before endpoint rounding, retaining
    /// reversible exchange beneath large unchanged field/entropy biases.
    pub fn thermal_offset_change_j(
        &self,
        after: &Self,
        before_displacement_c_m2: f64,
        after_displacement_c_m2: f64,
        volume_m3: f64,
    ) -> Result<f64, FieldStorageError> {
        self.admit_source_refresh(after)?;
        let volume = volume(volume_m3)?;
        let before_d = before_displacement_c_m2.abs();
        let after_d = after_displacement_c_m2.abs();
        let before = self.dielectric.evaluate_displacement(before_d)?;
        after.dielectric.evaluate_displacement(after_d)?;
        let a = self.scale_slope_per_k();
        if a == 0.0 || before_d == 0.0 && after_d == 0.0 {
            return Ok(0.0);
        }
        let dt = after.temperature_k() - self.temperature_k();
        let ds = finite(a * dt)?;
        if dt != 0.0 && ds == 0.0 {
            return Err(invalid("nonzero dielectric scale work underflows"));
        }
        // At the final scale, map the original D to the SAME reference state:
        // D_mapped = D_before*(1+ds/s_before). Integrate from actual D_after
        // using the explicit increment, including sub-ulp changes in D_mapped.
        // No fictitious intermediate D at the final T must fit source support.
        let scale_increment = ds / self.displacement_scale();
        let scale_increment_error = (-scale_increment).mul_add(self.displacement_scale(), ds)
            / self.displacement_scale()
            + a.mul_add(dt, -ds) / self.displacement_scale();
        let increment = finite(
            before_d.mul_add(scale_increment, before_d - after_d)
                + before_d * scale_increment_error,
        )?;
        let (coenergy_change, coenergy_error) = after.dielectric.coenergy_displacement_increment(
            after_d,
            increment,
            before_d == self.dielectric.maximum_displacement_c_m2(),
        )?;
        let (mut sum, mut correction) = (0.0, 0.0);
        add_product(
            &mut sum,
            &mut correction,
            scale_increment,
            before.coenergy_density_j_m3,
        )?;
        add_product(
            &mut sum,
            &mut correction,
            scale_increment_error,
            before.coenergy_density_j_m3,
        )?;
        let weighted_temperature = finite(a * after.temperature_k())?;
        let coefficient = weighted_temperature / after.displacement_scale();
        let coefficient_error = (-coefficient)
            .mul_add(after.displacement_scale(), weighted_temperature)
            / after.displacement_scale()
            + a.mul_add(after.temperature_k(), -weighted_temperature) / after.displacement_scale();
        add_product(&mut sum, &mut correction, -coefficient, coenergy_change)?;
        add_product(&mut sum, &mut correction, -coefficient, coenergy_error)?;
        add_product(
            &mut sum,
            &mut correction,
            -coefficient_error,
            coenergy_change,
        )?;
        let result = finite((sum + correction) * volume)?;
        if sum + correction != 0.0 && result == 0.0 {
            return Err(invalid("nonzero dielectric thermal offset underflows"));
        }
        Ok(result)
    }
}
