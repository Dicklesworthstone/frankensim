//! Fixed-condition electroquasistatic dielectric storage with Ohmic leakage.
//!
//! The Debye law owns D and polarization storage. The independently selected
//! DC conductor owns J = sigma E. Combining their integrated total current
//! gives D1 + dt J_mid = (epsilon_mid + dt sigma/2) E1 + history,
//! with history = Debye_history + dt sigma E0/2. The spatial owner must enforce
//! charge continuity using this total-current moment rather than conserving
//! dielectric free charge in a leaking cell.
//!
//! Polarization heat and Joule heat remain separate. Conductivity adds no
//! storage and must not include the loss already represented by Debye poles.
//! Both source roles retain one identical complete material snapshot; frequency-
//! or field-dependent conductivity cannot be frozen into this Ohmic law.
//! Midpoint is passive and second order, but stiff charge relaxation can change
//! sign between steps. It is not a monotone exponential relaxation method.
//! This material law supplies neither induction nor a full-wave validity test.

use fs_matdb::PropertyKey;
use fs_qty::QuantitySpec;

use crate::conductor::ELECTRICAL_CONDUCTIVITY_DIMS;

use super::{
    DebyePoleProperties, ELECTRIC_PERMITTIVITY_DIMS, FieldStorageError, MaterialCard,
    MaterialPropertySelection, QueryPoint, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, Time, admit_names,
    bulk::{DebyeMidpointStep, ResolvedDebyeMaterial},
    finite, invalid, resolve_material_state_point,
};

/// Scalar bulk leakage combined with isotropic multipole polarization.
/// No geometry, material defaults, source refit or parameter refresh is implied.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLeakyDielectric {
    dielectric: ResolvedDebyeMaterial,
    conductivity_property: String,
    conductivity_s_m: f64,
}

/// Local recovery after a total-current field solve, in the original pole order.
/// Densities must be integrated using the spatial owner's physical cell volume.
#[derive(Clone, Debug, PartialEq)]
pub struct LeakyDielectricMidpointStep<const N: usize> {
    /// Original Debye recovery, including storage and polarization heat only.
    pub dielectric: DebyeMidpointStep<N>,
    /// Midpoint DC conduction current density [A/m2], excluding polarization.
    pub ohmic_current_density_a_m2: [f64; N],
    /// Integrated midpoint conduction current [C/m2].
    pub integrated_ohmic_current_c_m2: [f64; N],
    /// Stable Delta D + dt J_mid, using actual retained polarization increments.
    pub total_current_increment_c_m2: [f64; N],
    /// Nonnegative dt sigma |E_mid|^2; excludes every Debye pole loss [J/m3].
    pub joule_heat_j_m3: f64,
    /// Polarization heat plus Joule heat, each counted exactly once [J/m3].
    pub total_heat_j_m3: f64,
    /// Midpoint total-current work, including the Ohmic contribution [J/m3].
    pub supplied_work_j_m3: f64,
    /// Debye storage change + total heat - supplied work [J/m3].
    pub energy_residual_j_m3: f64,
}

/// Admit an independent DC source in the original complete storage bundle.
/// Both relaxation and inertial polarization use this same conductivity gate.
pub(super) fn dc_conductivity<'a>(
    material: &ResolvedMaterialStatePoint,
    conductivity_property: &str,
    storage_properties: impl IntoIterator<Item = &'a str>,
) -> Result<f64, FieldStorageError> {
    if conductivity_property.trim().is_empty()
        || storage_properties
            .into_iter()
            .any(|name| name == conductivity_property)
    {
        return Err(invalid(
            "DC leakage requires its own nonempty source property name",
        ));
    }
    let key = PropertyKey::new(conductivity_property, ELECTRICAL_CONDUCTIVITY_DIMS);
    let source = material
        .property_by_key(&key)
        .ok_or_else(|| invalid("resolved DC leakage conductivity is missing"))?;
    let conductivity_s_m = source.value_si();
    if source.requirement().quantity() != QuantitySpec::dimensional(ELECTRICAL_CONDUCTIVITY_DIMS)
        || !source.is_static_at_fixed_external_conditions()
        || !conductivity_s_m.is_finite()
        || conductivity_s_m < 0.0
    {
        return Err(invalid(
            "DC leakage requires nonnegative conductivity at fixed external conditions",
        ));
    }
    Ok(conductivity_s_m)
}

impl ResolvedLeakyDielectric {
    /// Adopt a DC conductivity from the complete source bundle already retained
    /// by the Debye law. No second resolution or duplicate source bundle is made.
    /// Zero conductivity is an explicitly sourced perfect-insulator limit.
    ///
    /// Conductivity must mean DC Ohmic leakage. A dielectric loss tangent or
    /// frequency-dependent total loss cannot be reused as DC conductivity here.
    ///
    /// # Errors
    /// Missing/ambiguous names, incorrect dimensions, nonfinite or negative
    /// conductivity, or unowned frequency/field/strain/stress dependence refuse.
    pub fn try_from_dielectric(
        dielectric: ResolvedDebyeMaterial,
        conductivity_property: &str,
    ) -> Result<Self, FieldStorageError> {
        let conductivity_s_m = dc_conductivity(
            dielectric.material(),
            conductivity_property,
            core::iter::once(dielectric.high_frequency_property()).chain(
                dielectric.pole_properties().iter().flat_map(|pole| {
                    [
                        pole.permittivity_increment.as_str(),
                        pole.relaxation_time.as_str(),
                    ]
                }),
            ),
        )?;
        Ok(Self {
            dielectric,
            conductivity_property: conductivity_property.to_owned(),
            conductivity_s_m,
        })
    }

    /// Bind polarization and DC leakage roles in one existing multiphysics
    /// snapshot, preserving every original property and source-use receipt.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        high_frequency_property: &str,
        pole_properties: &[DebyePoleProperties],
        conductivity_property: &str,
    ) -> Result<Self, FieldStorageError> {
        Self::try_from_dielectric(
            ResolvedDebyeMaterial::try_from_material(
                material,
                high_frequency_property,
                pole_properties,
            )?,
            conductivity_property,
        )
    }

    /// Resolve all storage, memory and leakage requirements together. DC
    /// conductivity alone permits zero; every Debye storage/time remains positive.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        high_frequency_property: &str,
        pole_properties: &[DebyePoleProperties],
        conductivity_property: &str,
    ) -> Result<Self, FieldStorageError> {
        admit_names(high_frequency_property, pole_properties)?;
        let positive = |name: &str, dims| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        };
        let mut requirements = vec![positive(
            high_frequency_property,
            ELECTRIC_PERMITTIVITY_DIMS,
        )?];
        for pole in pole_properties {
            requirements.push(positive(
                &pole.permittivity_increment,
                ELECTRIC_PERMITTIVITY_DIMS,
            )?);
            requirements.push(positive(&pole.relaxation_time, Time::DIMS)?);
        }
        requirements.push(ScalarPropertyRequirement::try_new(
            conductivity_property,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::NonNegative,
        )?);
        Self::try_from_material(
            resolve_material_state_point(card, point, &requirements, selection)?,
            high_frequency_property,
            pole_properties,
            conductivity_property,
        )
    }

    /// Complete original source bundle; no synthetic permittivity or fitted
    /// conductivity is inserted into the material database.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.dielectric.material()
    }

    /// Retained polarization law and all its original parameter roles.
    pub const fn dielectric(&self) -> &ResolvedDebyeMaterial {
        &self.dielectric
    }

    /// Exact source property selected as DC conductivity.
    pub fn conductivity_property(&self) -> &str {
        &self.conductivity_property
    }

    /// Absolute source-backed scalar DC conductivity [S/m], including zero.
    pub const fn conductivity_s_m(&self) -> f64 {
        self.conductivity_s_m
    }

    /// Move the original source bundle once into a field owner's existing Debye
    /// storage, with the admitted conductivity role and value alongside it.
    /// Callers must not treat the algorithmic midpoint coefficient as storage.
    pub fn into_parts(self) -> (ResolvedDebyeMaterial, String, f64) {
        (
            self.dielectric,
            self.conductivity_property,
            self.conductivity_s_m,
        )
    }

    /// Positive local coefficient for the total-current midpoint field solve.
    /// This algorithmic coefficient is not stored electric permittivity.
    pub fn midpoint_permittivity_f_m(&self, dt_s: f64) -> Result<f64, FieldStorageError> {
        let dielectric = self.dielectric.midpoint_permittivity_f_m(dt_s)?;
        finite(dielectric + self.half_interval_conductivity(dt_s)?)
    }

    /// Endpoint affine history satisfying
    /// D1 + dt J_mid = midpoint_permittivity E1 + midpoint_history.
    /// The plus sign on the old-field leakage moment is physical: the weak
    /// scalar-potential RHS obtains the corresponding minus K_sigma phi0 term.
    pub fn midpoint_history_c_m2<const N: usize>(
        &self,
        dt_s: f64,
        electric_before_v_m: [f64; N],
        polarization_before_c_m2: &[[f64; N]],
    ) -> Result<[f64; N], FieldStorageError> {
        let mut history = self.dielectric.midpoint_history_c_m2(
            dt_s,
            electric_before_v_m,
            polarization_before_c_m2,
        )?;
        let half = self.half_interval_conductivity(dt_s)?;
        for (value, electric) in history.iter_mut().zip(electric_before_v_m) {
            *value = finite(half.mul_add(electric, *value))?;
        }
        Ok(history)
    }

    /// Recover polarization and both irreversible heat contributions for one
    /// immutable interval. Only the existing Debye state output allocates.
    /// The current moment uses actual state increments, avoiding subtraction
    /// of large absolute D values under a small perturbation or DC bias.
    pub fn step_midpoint<const N: usize>(
        &self,
        dt_s: f64,
        electric_before_v_m: [f64; N],
        electric_after_v_m: [f64; N],
        polarization_before_c_m2: &[[f64; N]],
    ) -> Result<LeakyDielectricMidpointStep<N>, FieldStorageError> {
        let dielectric = self.dielectric.step_midpoint(
            dt_s,
            electric_before_v_m,
            electric_after_v_m,
            polarization_before_c_m2,
        )?;
        // Use precisely the same interval admission as operator preparation.
        self.midpoint_permittivity_f_m(dt_s)?;
        let mut ohmic_current_density_a_m2 = [0.0; N];
        let mut integrated_ohmic_current_c_m2 = [0.0; N];
        let mut total_current_increment_c_m2 = [0.0; N];
        let mut joule_terms = [0.0; N];
        for i in 0..N {
            let midpoint = f64::midpoint(electric_before_v_m[i], electric_after_v_m[i]);
            ohmic_current_density_a_m2[i] = finite(self.conductivity_s_m * midpoint)?;
            integrated_ohmic_current_c_m2[i] = product3(dt_s, self.conductivity_s_m, midpoint)?;
            joule_terms[i] = finite(midpoint * integrated_ohmic_current_c_m2[i])?;
            if self.conductivity_s_m > 0.0
                && midpoint != 0.0
                && (ohmic_current_density_a_m2[i] == 0.0
                    || integrated_ohmic_current_c_m2[i] == 0.0
                    || joule_terms[i] == 0.0)
            {
                return Err(invalid("nonzero Ohmic current or Joule heat underflows"));
            }
            let field_increment = finite(electric_after_v_m[i] - electric_before_v_m[i])?;
            let instantaneous =
                finite(self.dielectric.high_frequency_permittivity_f_m() * field_increment)?;
            total_current_increment_c_m2[i] = compensated_sum(
                [instantaneous, integrated_ohmic_current_c_m2[i]]
                    .into_iter()
                    .chain(
                        dielectric
                            .polarization_after_c_m2
                            .iter()
                            .zip(polarization_before_c_m2)
                            .map(|(after, before)| after[i] - before[i]),
                    ),
            )?;
        }
        let joule_heat_j_m3 = compensated_sum(joule_terms)?;
        let total_heat_j_m3 = finite(dielectric.heat_j_m3 + joule_heat_j_m3)?;
        let supplied_work_j_m3 = finite(dielectric.supplied_work_j_m3 + joule_heat_j_m3)?;
        // Ohmic work and heat are the same computed quantity. Cancel that
        // identity before reducing a possibly much smaller dielectric residual.
        let energy_residual_j_m3 = dielectric.energy_residual_j_m3;
        Ok(LeakyDielectricMidpointStep {
            dielectric,
            ohmic_current_density_a_m2,
            integrated_ohmic_current_c_m2,
            total_current_increment_c_m2,
            joule_heat_j_m3,
            total_heat_j_m3,
            supplied_work_j_m3,
            energy_residual_j_m3,
        })
    }

    fn half_interval_conductivity(&self, dt_s: f64) -> Result<f64, FieldStorageError> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid(
                "leaky dielectric interval must be finite and positive",
            ));
        }
        if self.conductivity_s_m == 0.0 {
            return Ok(0.0);
        }
        // Halve the larger factor first: this avoids both dt*sigma overflow
        // and double rounding a representable subnormal half-product to zero.
        let half = if dt_s >= self.conductivity_s_m {
            (0.5 * dt_s) * self.conductivity_s_m
        } else {
            dt_s * (0.5 * self.conductivity_s_m)
        };
        if half.is_finite() && half > 0.0 {
            Ok(half)
        } else {
            Err(invalid("midpoint Ohmic coefficient is unrepresentable"))
        }
    }
}

// Multiplying the smallest and largest magnitudes first avoids avoidable
// intermediate overflow/underflow for three finite factors.
fn product3(a: f64, b: f64, c: f64) -> Result<f64, FieldStorageError> {
    let mut factors = [a, b, c];
    factors.sort_by(|x, y| x.abs().total_cmp(&y.abs()));
    finite((factors[0] * factors[2]) * factors[1])
}

fn compensated_sum(values: impl IntoIterator<Item = f64>) -> Result<f64, FieldStorageError> {
    let mut sum = 0.0_f64;
    let mut correction = 0.0;
    for value in values {
        let value = finite(value)?;
        let next = finite(sum + value)?;
        correction = finite(
            correction
                + if sum.abs() >= value.abs() {
                    (sum - next) + value
                } else {
                    (value - next) + sum
                },
        )?;
        sum = next;
    }
    finite(sum + correction)
}
