//! Reversible longitudinal magnetostriction about an explicitly held DC bias.
//!
//! Incremental source convention: `strain = s_H stress + d H`,
//! `B = d stress + mu_T H`, where every field/strain/stress here is a perturbation
//! of the prepared, stress-free biased reference. The signed piezomagnetic d is
//! BOTH direct and converse response. See Deng and Dapino, *Magnetic flux biasing
//! of magnetostrictive sensors*, Smart Mater. Struct. 26 (2017) 055027, Eq. (1):
//! <https://doi.org/10.1088/1361-665X/aa688b>.
//!
//! This is a local reversible tangent, not a global law through H=0. The
//! original complete material query retains the nonzero bias field, temperature
//! and all fixed conditions. Independent variables strain and incremental B
//! give the positive perturbation energy
//! `w = strain^2/(2*s_H) + (B - e*strain)^2/(2*mu_axial)`,
//! `e=d/s_H`, `mu_axial=mu_T-d*d/s_H > 0`. Its derivatives supply BOTH efforts.
//! Absolute bias preparation energy, hysteresis, eddy loss and thermal changes
//! require their own physical owners and are not included in this energy.

use crate::field_storage::MAGNETIC_PERMEABILITY_DIMS;
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, QuantitySpec};

/// Incremental axial compliance at fixed H and free transverse stress [1/Pa].
pub const COMPLIANCE_33_H_PROPERTY: &str =
    "magnetostrictive_incremental_compliance_33_at_constant_field";
/// Signed reciprocal incremental d33 [m/A = T/Pa], about the stated H bias.
pub const PIEZOMAGNETIC_33_PROPERTY: &str = "magnetostrictive_incremental_piezomagnetic_33";
/// Incremental absolute mu33 at fixed stress [H/m], not relative permeability.
pub const PERMEABILITY_33_T_PROPERTY: &str =
    "magnetostrictive_incremental_permeability_33_at_constant_stress";
/// Axial compliance dimensions.
pub const COMPLIANCE_DIMS: Dims = Dims([1, -1, 2, 0, 0, 0]);
/// Reciprocal piezomagnetic coefficient dimensions [m/A].
pub const PIEZOMAGNETIC_DIMS: Dims = Dims([1, 0, 0, 0, -1, 0]);
/// Magnetic field strength dimensions [A/m]; bias is distinct from flux density.
pub const MAGNETIC_FIELD_DIMS: Dims = Dims([-1, 0, 0, 0, 1, 0]);

/// Original source refusal or unsupported/unstable incremental law.
#[derive(Debug)]
pub enum MagnetostrictiveError {
    /// Original source resolver rejected the complete material point.
    Material(MaterialStatePointError),
    /// Unsupported source convention or unrepresentable coupled storage.
    Invalid(&'static str),
}
impl From<MaterialStatePointError> for MagnetostrictiveError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}
impl core::fmt::Display for MagnetostrictiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "incremental magnetostrictive material: {self:?}")
    }
}
impl std::error::Error for MagnetostrictiveError {}

/// One complete source snapshot at a specified, fixed magnetic bias. Axis 3 is
/// the longitudinal material direction; the rod's reference shape is already
/// equilibrated at this bias and zero incremental stress. This is not a full
/// tensor and supplies no off-axis, rotating-domain or magnetization history law.
#[derive(Clone, Debug, PartialEq)]
pub struct LongitudinalMagnetostrictive {
    material: ResolvedMaterialStatePoint,
    bias_axis: String,
    bias_field_a_m: f64,
    compliance_h: f64,
    d: f64,
    mu_t: f64,
    mu_axial: f64,
    e: f64,
    density: f64,
}

/// Perturbation energy and work-conjugate derivatives of the SAME source law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MagnetostrictiveResponse {
    /// Axial tensile stress increment, partial w / partial strain [Pa].
    pub stress_increment_pa: f64,
    /// Axial H increment, partial w / partial B [A/m].
    pub field_increment_a_m: f64,
    /// Reversible perturbation energy per reference volume [J/m^3].
    pub incremental_energy_j_m3: f64,
}
impl LongitudinalMagnetostrictive {
    /// Resolve the complete four-role material at an explicitly typed nonzero
    /// H-bias coordinate. A source curve in that coordinate selects a local
    /// tangent; it is not queried at the evolving perturbation field. The caller
    /// must separately declare the admitted small-signal amplitude envelope.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        bias_axis: &str,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MagnetostrictiveError> {
        let roles = [
            (
                COMPLIANCE_33_H_PROPERTY,
                COMPLIANCE_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                PIEZOMAGNETIC_33_PROPERTY,
                PIEZOMAGNETIC_DIMS,
                ScalarAdmissibility::Finite,
            ),
            (
                PERMEABILITY_33_T_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                DENSITY_PROPERTY,
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
        ];
        let required = roles
            .into_iter()
            .map(|(name, dims, domain)| ScalarPropertyRequirement::try_new(name, dims, domain))
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_material_state(
            resolve_material_state_point(card, point, &required, selection)?,
            bias_axis,
        )
    }

    /// Preserve an existing complete source bundle. This owner explicitly fixes
    /// its H-bias while all other axes use the shared fixed-external-condition
    /// policy. Unowned stress, frequency, rate, or dynamic field support refuses;
    /// a point-only scalar cannot silently become a finite-amplitude tangent law.
    pub fn from_material_state(
        material: ResolvedMaterialStatePoint,
        bias_axis: &str,
    ) -> Result<Self, MagnetostrictiveError> {
        let invalid = MagnetostrictiveError::Invalid;
        let bias = material
            .query_point()
            .iter()
            .find(|(axis, _)| axis == bias_axis)
            .map(|(_, value)| *value)
            .ok_or(invalid("explicit magnetic H bias is missing"))?;
        if !bias.is_finite()
            || bias == 0.0
            || bias_axis.is_empty()
            || material.axis_quantities().get(bias_axis)
                != Some(&QuantitySpec::dimensional(MAGNETIC_FIELD_DIMS))
        {
            return Err(invalid(
                "nonzero finite H bias with explicit A/m quantity is required",
            ));
        }
        let coefficient = |name, dims| {
            let p = material
                .property(name)
                .ok_or(invalid("missing incremental magnetostrictive role"))?;
            if p.requirement().quantity() != QuantitySpec::dimensional(dims)
                || !p.has_only_owned_and_fixed_external_conditions(
                    bias_axis,
                    QuantitySpec::dimensional(MAGNETIC_FIELD_DIMS),
                )
            {
                return Err(invalid(
                    "incremental coefficients require their declared SI roles, held bias and fixed external conditions",
                ));
            }
            Ok(p.value_si())
        };
        let compliance_h = coefficient(COMPLIANCE_33_H_PROPERTY, COMPLIANCE_DIMS)?;
        let d = coefficient(PIEZOMAGNETIC_33_PROPERTY, PIEZOMAGNETIC_DIMS)?;
        let mu_t = coefficient(PERMEABILITY_33_T_PROPERTY, MAGNETIC_PERMEABILITY_DIMS)?;
        let density = coefficient(DENSITY_PROPERTY, Density::DIMS)?;
        let e = d / compliance_h;
        let reduction = d * e;
        let mu_axial = mu_t - reduction;
        if [
            compliance_h,
            mu_t,
            mu_axial,
            density,
            1.0 / compliance_h,
            1.0 / mu_axial,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
            || !d.is_finite()
            || !e.is_finite()
            || !reduction.is_finite()
            || d != 0.0 && (e == 0.0 || reduction == 0.0)
        {
            return Err(invalid(
                "stable representable reciprocal block requires s_H>0 and mu_T-d^2/s_H>0",
            ));
        }
        let result = Self {
            material,
            bias_axis: bias_axis.to_owned(),
            bias_field_a_m: bias,
            compliance_h,
            d,
            mu_t,
            mu_axial,
            e,
            density,
        };
        if !result.modulus_at_fixed_flux_pa().is_finite() {
            return Err(invalid("fixed-flux modulus is unrepresentable"));
        }
        Ok(result)
    }
    /// Complete original source bundle, including its fixed bias and conditions.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Caller-selected, explicitly typed H-bias coordinate name.
    pub fn bias_axis(&self) -> &str {
        &self.bias_axis
    }
    /// Signed fixed DC H bias [A/m], distinct from the evolving field increment.
    pub const fn bias_field_a_m(&self) -> f64 {
        self.bias_field_a_m
    }
    /// Longitudinal compliance s_H at fixed H [1/Pa].
    pub const fn compliance_at_fixed_field_pa_inv(&self) -> f64 {
        self.compliance_h
    }
    /// Signed direct/converse incremental coefficient d33 [m/A].
    pub const fn piezomagnetic_coefficient_m_a(&self) -> f64 {
        self.d
    }
    /// Absolute incremental permeability at fixed axial stress [H/m].
    pub const fn permeability_at_fixed_stress_h_m(&self) -> f64 {
        self.mu_t
    }
    /// Permeability at fixed axial strain and free transverse stress [H/m].
    pub const fn axially_clamped_permeability_h_m(&self) -> f64 {
        self.mu_axial
    }
    /// Reciprocal coefficient e=d/s_H [T] in B=mu_axial*H+e*strain.
    pub const fn stress_coefficient_t(&self) -> f64 {
        self.e
    }
    /// Density of the same prepared reference material [kg/m^3].
    pub const fn density_kg_m3(&self) -> f64 {
        self.density
    }
    /// Incremental axial modulus at fixed H [Pa].
    pub fn modulus_at_fixed_field_pa(&self) -> f64 {
        1.0 / self.compliance_h
    }
    /// Incremental axial modulus at fixed B [Pa].
    pub fn modulus_at_fixed_flux_pa(&self) -> f64 {
        1.0 / self.compliance_h + self.e * (self.e / self.mu_axial)
    }

    /// Evaluate perturbation energy at independent axial strain and flux-density
    /// increment. No scalar gain is allowed to differ between actuation and sensing.
    pub fn response(
        &self,
        strain: f64,
        flux_increment_t: f64,
    ) -> Result<MagnetostrictiveResponse, MagnetostrictiveError> {
        let elastic_stress = strain / self.compliance_h;
        let residual = flux_increment_t - self.e * strain;
        let field_increment_a_m = residual / self.mu_axial;
        let stress_increment_pa = elastic_stress - self.e * field_increment_a_m;
        let elastic = 0.5 * strain * elastic_stress;
        let magnetic = 0.5 * residual * field_increment_a_m;
        let incremental_energy_j_m3 = elastic + magnetic;
        if [
            strain,
            flux_increment_t,
            field_increment_a_m,
            stress_increment_pa,
            incremental_energy_j_m3,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || strain != 0.0 && elastic == 0.0
            || residual != 0.0 && magnetic == 0.0
        {
            return Err(MagnetostrictiveError::Invalid(
                "incremental magnetostrictive response is unrepresentable",
            ));
        }
        Ok(MagnetostrictiveResponse {
            stress_increment_pa,
            field_increment_a_m,
            incremental_energy_j_m3,
        })
    }
}
