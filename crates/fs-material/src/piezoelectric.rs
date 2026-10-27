//! Source-bound, reversible longitudinal piezoelectricity at a fixed material state.
//!
//! The declared material direction is axis 3, with transverse stress free and
//! electric field and displacement along that axis. The source convention is
//! `strain = s_E stress + d E`, `D = d stress + epsilon_T E`. These are the
//! small-signal strain/charge equations, e.g. PI's *Fundamentals of Piezo
//! Technology*, Electromechanics and Assignment of Axes:
//! <https://www.pi-usa.us/en/expertise/technology/piezo-technology/fundamentals>.
//!
//! Eliminating stress gives `e = d/s_E`, `epsilon_axial = epsilon_T-d*d/s_E`.
//! The latter holds AXIAL strain fixed while transverse stress remains free;
//! it is not the fully clamped three-dimensional permittivity. Positive
//! `s_E` and `epsilon_axial` admit the reciprocal energy
//! `w(strain,D) = strain^2/(2 s_E) + (D-e strain)^2/(2 epsilon_axial)`.
//! Its derivatives, rather than two independently tuned gains, supply stress
//! and electric field. No material name selects a coefficient.

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, QuantitySpec};

use crate::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};

/// Reversible caloric response of the original coupled compliance block.
pub mod thermal;

/// Axial `s33^E` [1/Pa], with other stresses zero and electric field fixed.
pub const COMPLIANCE_33_E_PROPERTY: &str = "piezoelectric_compliance_33_at_constant_field";
/// Signed `d33` [m/V = C/N] in the same oriented material axis as the other roles.
pub const STRAIN_CHARGE_33_PROPERTY: &str = "piezoelectric_strain_charge_33";
/// Absolute `epsilon33^T` [F/m], with mechanical stresses fixed. Not relative epsilon.
pub const PERMITTIVITY_33_T_PROPERTY: &str = "piezoelectric_permittivity_33_at_constant_stress";
/// Inverse-pressure dimensions of elastic compliance.
pub const COMPLIANCE_DIMS: Dims = Dims([1, -1, 2, 0, 0, 0]);
/// Dimensions shared by the reciprocal longitudinal strain and charge coefficients.
pub const STRAIN_CHARGE_DIMS: Dims = Dims([-1, -1, 3, 0, 1, 0]);

/// Source selection, unsupported constitutive state or coupled-stability refusal.
#[derive(Debug)]
pub enum PiezoelectricError {
    /// Original material-card query failure.
    Material(MaterialStatePointError),
    /// Unsupported source convention or unrepresentable coupled response.
    Invalid(&'static str),
}
impl From<MaterialStatePointError> for PiezoelectricError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for PiezoelectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "longitudinal piezoelectric material: {self:?}")
    }
}
impl std::error::Error for PiezoelectricError {}

/// One complete fixed source snapshot; coefficients are never refreshed silently
/// during a trajectory. The property names declare the common longitudinal
/// direction and stress/field conventions. This is a one-dimensional reduction,
/// not an inference of a complete tensor or its off-axis response.
#[derive(Clone, Debug, PartialEq)]
pub struct LongitudinalPiezoelectric {
    material: ResolvedMaterialStatePoint,
    compliance_e: f64,
    d: f64,
    epsilon_t: f64,
    epsilon_axial: f64,
    e: f64,
    density: f64,
}

/// Derivatives and value of the SAME electric-displacement/strain energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PiezoelectricResponse {
    /// Axial tensile stress, `dw/d(strain)` [Pa].
    pub stress_pa: f64,
    /// Field along material axis 3, `dw/dD` [V/m].
    pub electric_field_v_m: f64,
    /// Reversible energy per reference volume [J/m^3].
    pub stored_energy_j_m3: f64,
}

impl LongitudinalPiezoelectric {
    /// Resolve all four roles from one immutable card at the same full query.
    /// Signed or zero d is supported. Missing directional conventions, relative
    /// permittivity, ambiguous claims and unsupported source domains refuse.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<Self, PiezoelectricError> {
        let roles = [
            (
                COMPLIANCE_33_E_PROPERTY,
                COMPLIANCE_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                STRAIN_CHARGE_33_PROPERTY,
                STRAIN_CHARGE_DIMS,
                ScalarAdmissibility::Finite,
            ),
            (
                PERMITTIVITY_33_T_PROPERTY,
                ELECTRIC_PERMITTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                DENSITY_PROPERTY,
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
        ];
        let requirements = roles
            .into_iter()
            .map(|(name, dims, admissibility)| {
                ScalarPropertyRequirement::try_new(name, dims, admissibility)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_material_state(resolve_material_state_point(
            card,
            point,
            &requirements,
            selection,
        )?)
    }

    /// Consume an existing shared source bundle without dropping its other roles.
    /// Static temperature/pressure/composition conditions remain attached. Field,
    /// frequency, strain or stress dependence AND applicability bounds cannot be
    /// frozen into this constant-coefficient reversible model.
    pub fn from_material_state(
        material: ResolvedMaterialStatePoint,
    ) -> Result<Self, PiezoelectricError> {
        let invalid = PiezoelectricError::Invalid;
        let coefficient = |name, dims| {
            let property = material
                .property(name)
                .ok_or(invalid("missing longitudinal piezoelectric role"))?;
            if property.requirement().quantity() != QuantitySpec::dimensional(dims)
                || !property.is_static_at_fixed_external_conditions()
            {
                return Err(invalid(
                    "piezoelectric coefficients require the declared SI convention at fixed external conditions",
                ));
            }
            Ok(property.value_si())
        };
        let compliance_e = coefficient(COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS)?;
        let d = coefficient(STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS)?;
        let epsilon_t = coefficient(PERMITTIVITY_33_T_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS)?;
        let density = coefficient(DENSITY_PROPERTY, Density::DIMS)?;
        let e = d / compliance_e;
        let reduction = d * e;
        let epsilon_axial = epsilon_t - reduction;
        if [
            compliance_e,
            epsilon_t,
            epsilon_axial,
            density,
            1.0 / compliance_e,
            1.0 / epsilon_axial,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
            || !d.is_finite()
            || !e.is_finite()
            || !reduction.is_finite()
            || d != 0.0 && (e == 0.0 || reduction == 0.0)
        {
            return Err(invalid(
                "positive compliance, density and epsilon_T-d^2/s_E are required for stable representable coupled storage",
            ));
        }
        let result = Self {
            material,
            compliance_e,
            d,
            epsilon_t,
            epsilon_axial,
            e,
            density,
        };
        if !result.open_circuit_modulus_pa().is_finite() {
            return Err(invalid("open-circuit stiffness is not representable"));
        }
        Ok(result)
    }

    /// Original complete source bundle, quantity descriptors and selected claims.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Constant-field axial compliance `s33^E` [1/Pa].
    pub const fn compliance_at_fixed_field_pa_inv(&self) -> f64 {
        self.compliance_e
    }
    /// Signed reciprocal `d33` [m/V = C/N].
    pub const fn strain_charge_coefficient_m_v(&self) -> f64 {
        self.d
    }
    /// Stress-free absolute longitudinal permittivity `epsilon33^T` [F/m].
    pub const fn permittivity_at_fixed_stress_f_m(&self) -> f64 {
        self.epsilon_t
    }
    /// Permittivity at fixed axial strain and free transverse stress [F/m].
    pub const fn axially_clamped_permittivity_f_m(&self) -> f64 {
        self.epsilon_axial
    }
    /// Signed reduced stress coefficient `d33/s33^E` [C/m^2].
    pub const fn stress_coefficient_c_m2(&self) -> f64 {
        self.e
    }
    /// Density at the same original material conditions [kg/m^3].
    pub const fn density_kg_m3(&self) -> f64 {
        self.density
    }
    /// Axial modulus at fixed field, including a short circuit [Pa].
    pub fn short_circuit_modulus_pa(&self) -> f64 {
        1.0 / self.compliance_e
    }
    /// Axial modulus at fixed electric displacement, including open electrodes [Pa].
    pub fn open_circuit_modulus_pa(&self) -> f64 {
        1.0 / self.compliance_e + self.e * (self.e / self.epsilon_axial)
    }

    /// Evaluate the energy and its work-conjugate derivatives. Strain and D are
    /// independent variables. The physical consumer must additionally enforce
    /// its small-signal strain/field envelope; no saturation, depoling, leakage,
    /// dielectric loss, hysteresis, heating or variable-coefficient work is added.
    pub fn response(
        &self,
        strain: f64,
        electric_displacement_c_m2: f64,
    ) -> Result<PiezoelectricResponse, PiezoelectricError> {
        let elastic_stress = strain / self.compliance_e;
        let residual_displacement = electric_displacement_c_m2 - self.e * strain;
        let electric_field_v_m = residual_displacement / self.epsilon_axial;
        let stress_pa = elastic_stress - self.e * electric_field_v_m;
        let elastic = 0.5 * strain * elastic_stress;
        let electric = 0.5 * residual_displacement * electric_field_v_m;
        let stored_energy_j_m3 = elastic + electric;
        if [
            strain,
            electric_displacement_c_m2,
            electric_field_v_m,
            stress_pa,
            stored_energy_j_m3,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || strain != 0.0 && elastic == 0.0
            || residual_displacement != 0.0 && electric == 0.0
        {
            return Err(PiezoelectricError::Invalid(
                "piezoelectric response is not representable",
            ));
        }
        Ok(PiezoelectricResponse {
            stress_pa,
            electric_field_v_m,
            stored_energy_j_m3,
        })
    }
}
