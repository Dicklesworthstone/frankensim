//! Explicit thermal eigenstrain and CLAMPED pyroelectric polarization.
//!
//! p is measured at fixed total strain and electric field. The stress-free,
//! zero-field response also contains secondary piezoelectricity e*alpha; a total
//! pyroelectric coefficient must be converted before it is supplied under these
//! clamped source roles. No spontaneous polarization or reference strain is
//! inferred. The caller explicitly identifies the temperature at which both
//! offsets vanish in its reference geometry and charge convention.
//!
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.mems/mems_ug_sme.08.3.html>

use super::*;
use crate::piezoelectric::tensor::TensorPiezoelectricResponse;
use crate::state_point::{
    INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_qty::Dims;

/// Clamped pyroelectric coefficient dimensions [C/(m² K)].
pub const PYROELECTRIC_CLAMPED_DIMS: Dims = Dims([-2, 0, 1, -1, 1, 0]);
const EXPANSION: &str = "piezoelectric_thermal_expansion_mandel";
const PYROELECTRIC: &str = "piezoelectric_pyroelectric_clamped";

fn norm<const N: usize>(values: &[f64; N]) -> Result<f64> {
    finite(
        values
            .iter()
            .fold(0.0_f64, |norm, value| norm.hypot(*value)),
    )
}

#[derive(Clone, Debug)]
pub(super) struct ThermalBiasSource {
    source: ResolvedMaterialStatePoint,
    reference_temperature_k: f64,
    expansion: [f64; 6],
    pyroelectric: [f64; 3],
    pub(super) expansion_norm: f64,
    pub(super) polarization_norm: f64,
    pub(super) maximum_temperature_change_k: f64,
    pub(super) maximum_expansion_curvature: f64,
    pub(super) expansion_slope_norm: f64,
    pub(super) maximum_pyroelectric_inverse_norm: f64,
}

#[derive(Clone, Debug)]
pub(super) struct ResolvedThermalBias {
    source: ResolvedMaterialStatePoint,
    reference_temperature_k: f64,
    pub(super) expansion: [f64; 6],
    pub(super) pyroelectric: [f64; 3],
    pub(super) expansion_curvature: f64,
    thermal_strain: [f64; 6],
    polarization: [f64; 3],
    stress_offset: [f64; 6],
    identity: ContentHash,
}

impl ThermalTensorPiezoelectric {
    /// Six Mandel linear expansion coefficients [1/K], followed by three
    /// Cartesian CLAMPED pyroelectric coefficients [C/(m² K)]. All nine named
    /// components, including zeros and Mandel shear factors, require explicit
    /// source claims in the original tensor's SOURCE frame.
    pub fn thermal_bias_requirements() -> Result<Vec<ScalarPropertyRequirement>> {
        [
            (EXPANSION, 6, INVERSE_TEMPERATURE_DIMS),
            (PYROELECTRIC, 3, PYROELECTRIC_CLAMPED_DIMS),
        ]
        .into_iter()
        .flat_map(|(prefix, count, dims)| {
            (1..=count).map(move |i| {
                ScalarPropertyRequirement::try_new(
                    format!("{prefix}_{i}"),
                    dims,
                    ScalarAdmissibility::Finite,
                )
                .map_err(PiezoelectricError::Material)
            })
        })
        .collect()
    }

    /// Attach constant, explicitly sourced thermal expansion and clamped
    /// pyroelectric coefficients over this law's complete temperature interval.
    /// The same immutable card, physical context and original frame are used.
    /// Reference temperature defines ZERO thermal strain and polarization; it
    /// must lie inside the interval. The native tensor itself remains unchanged.
    ///
    /// S=alpha*(T-Tref), E=0 has zero field free energy. A consuming caloric
    /// chart therefore represents stress-free, zero-field heat capacity.
    /// Supplying an unconverted clamped heat capacity would count constraint
    /// effects twice. Total pyroelectric coefficients are not clamped sources:
    /// with constant e the free shorted response has slope p_clamped+e*alpha.
    pub fn with_thermal_bias(
        mut self,
        reference_temperature_k: f64,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        if self.bias.is_some()
            || !reference_temperature_k.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1])
                .contains(&reference_temperature_k)
        {
            return Err(PiezoelectricError::Invalid(
                "thermal bias needs one in-range reference temperature and cannot replace an attached bias",
            ));
        }
        let requirements = Self::thermal_bias_requirements()?;
        let source = resolve_material_state_point(
            &self.card,
            &self.original.material().query(),
            &requirements,
            selection,
        )?;
        for requirement in &requirements {
            let property = source
                .property(requirement.name())
                .expect("all bias roles resolved");
            if let Some(original) = self.original.material().property(requirement.name()) {
                if original.requirement().quantity() != property.requirement().quantity()
                    || original.answer().receipt.selected != property.answer().receipt.selected
                    || original.value_si().to_bits() != property.value_si().to_bits()
                {
                    return Err(PiezoelectricError::Invalid(
                        "thermal bias cannot replace an already selected original source",
                    ));
                }
            }
            let (_, slope, _) = affine_source_path(
                &self.card,
                &source,
                requirement.name(),
                &self.temperature_axis,
                self.temperature_range_k,
                false,
            )
            .map_err(field_error)?;
            if slope != 0.0 {
                return Err(PiezoelectricError::Invalid(
                    "thermal expansion and clamped pyroelectric coefficients must be constant on the admitted interval",
                ));
            }
        }
        let value = |prefix: &str, i: usize| {
            source
                .property(&format!("{prefix}_{}", i + 1))
                .expect("complete thermal bias source")
                .value_si()
        };
        let source_expansion = core::array::from_fn(|i| value(EXPANSION, i));
        let source_pyroelectric = core::array::from_fn(|i| value(PYROELECTRIC, i));
        let q = *self.original.source_to_frame_rotation();
        let expansion = apply(&mandel_rotation(q), &source_expansion);
        let pyroelectric = apply(&q, &source_pyroelectric);
        let up = 1.0 + 4096.0 * f64::EPSILON;
        let mut bias = ThermalBiasSource {
            source,
            reference_temperature_k,
            expansion,
            pyroelectric,
            expansion_norm: finite(norm(&expansion)? * up)?,
            polarization_norm: finite(norm(&pyroelectric)? * up)?,
            maximum_temperature_change_k: self
                .temperature_range_k
                .iter()
                .fold(0.0_f64, |bound, t| {
                    bound.max((t - reference_temperature_k).abs())
                }),
            maximum_expansion_curvature: 0.0,
            expansion_slope_norm: finite(norm(&apply(&self.slopes.stiffness, &expansion))? * up)?,
            maximum_pyroelectric_inverse_norm: 0.0,
        };
        if norm(&source_expansion)? != 0.0 && bias.expansion_norm == 0.0
            || norm(&source_pyroelectric)? != 0.0 && bias.polarization_norm == 0.0
        {
            return Err(PiezoelectricError::Invalid(
                "rotated thermal bias coefficients underflow",
            ));
        }
        for temperature in self.temperature_range_k {
            let at = self.at_temperature(temperature)?;
            let material = at.material();
            let elastic = finite(
                2.0 * material
                    .field_response(expansion, [0.0; 3])?
                    .stored_energy_j_m3,
            )?;
            let electric = norm(&material.dielectric_inverse_factor_action(pyroelectric)?)?;
            bias.maximum_expansion_curvature = bias.maximum_expansion_curvature.max(elastic);
            bias.maximum_pyroelectric_inverse_norm =
                bias.maximum_pyroelectric_inverse_norm.max(electric);
        }
        // alpha.C.alpha is affine. p.epsilon^-1.p is convex on SPD matrices,
        // so the endpoint maxima bound the entire affine temperature interval.
        bias.maximum_expansion_curvature = finite(bias.maximum_expansion_curvature * up)?;
        bias.maximum_pyroelectric_inverse_norm =
            finite(bias.maximum_pyroelectric_inverse_norm * up)?;
        self.bias = Some(bias);
        Ok(self)
    }
}

impl ThermalBiasSource {
    pub(super) fn at_temperature(
        &self,
        card: &MaterialCard,
        material: &TensorPiezoelectric,
        temperature_k: f64,
    ) -> Result<ResolvedThermalBias> {
        let point = material.material().query();
        let refreshed = self.source.at_query_point(card, &point)?;
        if self
            .source
            .properties()
            .iter()
            .zip(refreshed.properties())
            .any(|(a, b)| a.value_si() != b.value_si())
        {
            return Err(PiezoelectricError::Invalid(
                "thermal bias refresh changed an admitted constant source coefficient",
            ));
        }
        let source = if point == self.source.query() {
            self.source.clone()
        } else {
            refreshed
        };
        let delta = temperature_k - self.reference_temperature_k;
        let thermal_strain = self.expansion.map(|value| value * delta);
        let polarization = self.pyroelectric.map(|value| value * delta);
        if thermal_strain
            .iter()
            .chain(&polarization)
            .any(|value| !value.is_finite())
            || delta != 0.0
                && (self
                    .expansion
                    .iter()
                    .zip(thermal_strain)
                    .any(|(&a, value)| a != 0.0 && value == 0.0)
                    || self
                        .pyroelectric
                        .iter()
                        .zip(polarization)
                        .any(|(&p, value)| p != 0.0 && value == 0.0))
        {
            return Err(PiezoelectricError::Invalid(
                "thermal bias offsets are unrepresentable",
            ));
        }
        let stress_offset = apply(material.stiffness_mandel_pa(), &thermal_strain).map(|v| -v);
        if stress_offset.iter().any(|value| !value.is_finite()) {
            return Err(PiezoelectricError::Invalid(
                "thermal stress offset is unrepresentable",
            ));
        }
        let expansion_curvature = finite(
            2.0 * material
                .field_response(self.expansion, [0.0; 3])?
                .stored_energy_j_m3,
        )?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.tensor-piezoelectric-thermal-bias.v1");
        identity.update(material.identity().as_bytes());
        identity.update(source.identity().as_bytes());
        identity.update(&self.reference_temperature_k.to_bits().to_le_bytes());
        Ok(ResolvedThermalBias {
            source,
            reference_temperature_k: self.reference_temperature_k,
            expansion: self.expansion,
            pyroelectric: self.pyroelectric,
            expansion_curvature,
            thermal_strain,
            polarization,
            stress_offset,
            identity: identity.finalize(),
        })
    }
}

impl ResolvedThermalTensorPiezoelectric {
    /// Original selected bias claims at this sample's actual physical query.
    pub fn thermal_bias_source(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.bias.as_ref().map(|bias| &bias.source)
    }
    /// Explicit reference geometry/charge temperature; absent for an unbiased law.
    pub fn thermal_bias_reference_temperature_k(&self) -> Option<f64> {
        self.bias.as_ref().map(|bias| bias.reference_temperature_k)
    }
    /// Source, current tensor/frame, actual temperature and reference-anchor identity.
    pub fn thermal_bias_identity(&self) -> Option<ContentHash> {
        self.bias.as_ref().map(|bias| bias.identity)
    }
    /// Free thermal strain in current-frame Mandel coordinates, alpha*(T-Tref).
    pub fn thermal_strain_mandel(&self) -> [f64; 6] {
        self.bias
            .as_ref()
            .map_or([0.0; 6], |bias| bias.thermal_strain)
    }
    /// Direct CLAMPED polarization offset p*(T-Tref) [C/m²]. The secondary
    /// contribution e*alpha*(T-Tref) is already produced by actual free strain.
    pub fn clamped_polarization_c_m2(&self) -> [f64; 3] {
        self.bias
            .as_ref()
            .map_or([0.0; 3], |bias| bias.polarization)
    }
    /// Additive stress at S=E=0, -C*alpha*(T-Tref), in current Mandel axes [Pa].
    pub fn stress_offset_mandel_pa(&self) -> [f64; 6] {
        self.bias
            .as_ref()
            .map_or([0.0; 6], |bias| bias.stress_offset)
    }

    /// Stress, displacement and the positive strain/charge potential at total
    /// strain S and electric field E. The tangent tensors are unchanged.
    /// Native positive energy factors act on S-alpha*(T-Tref) and E; coupling
    /// still acts on TOTAL strain, preserving clamped versus total pyroelectricity.
    pub fn field_response(
        &self,
        strain_mandel: [f64; 6],
        electric_field_v_m: [f64; 3],
    ) -> Result<TensorPiezoelectricResponse> {
        let Some(bias) = &self.bias else {
            return self
                .material
                .field_response(strain_mandel, electric_field_v_m);
        };
        let elastic_strain = core::array::from_fn(|i| strain_mandel[i] - bias.thermal_strain[i]);
        let mut response = self
            .material
            .field_response(elastic_strain, electric_field_v_m)?;
        let piezoelectric = apply(self.material.coupling_c_m2(), &strain_mandel);
        let dielectric = apply(self.material.permittivity_f_m(), &electric_field_v_m);
        response.displacement_c_m2 =
            core::array::from_fn(|i| piezoelectric[i] + dielectric[i] + bias.polarization[i]);
        if response
            .displacement_c_m2
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(PiezoelectricError::Invalid(
                "biased piezoelectric displacement is unrepresentable",
            ));
        }
        Ok(response)
    }

    /// Evaluate the shared Helmholtz potential and its stress/field gradients
    /// at fixed TOTAL strain and electric displacement, including thermal bias.
    pub fn charge_response(
        &self,
        strain_mandel: [f64; 6],
        displacement_c_m2: [f64; 3],
    ) -> Result<TensorPiezoelectricResponse> {
        let Some(bias) = &self.bias else {
            return self
                .material
                .charge_response(strain_mandel, displacement_c_m2);
        };
        let piezoelectric = apply(self.material.coupling_c_m2(), &strain_mandel);
        let field = self.material.solve_dielectric(core::array::from_fn(|i| {
            displacement_c_m2[i] - piezoelectric[i] - bias.polarization[i]
        }))?;
        let mut response = self.field_response(strain_mandel, field)?;
        response.displacement_c_m2 = displacement_c_m2;
        Ok(response)
    }
}
