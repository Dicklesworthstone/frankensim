//! Source-bound finite membrane thermodynamics in original-area units.
//!
//! Transverse equilibrium is eliminated before taking the thermal response.
//! Its Schur correction changes both thermal stress and heat capacity. The
//! existing Green-strain map then supplies work-conjugate sheet resultants;
//! multiplying by the ORIGINAL thickness preserves the original inventory.
//! These are local derivatives of a held equilibrium free energy, not a heat
//! integrator, a finite thermal reference-work approximation, or a new caloric
//! source. Elastic coefficients remain those of the selected native model.

use super::expansion::{SheetExpansionFrame, SheetExpansionReference};
use super::{
    HyperelasticMembrane, VolumetricTreatment, admitted_representative, representative,
    response_from_plane_stress,
};
use crate::hyper::thermoelastic::ThermalReferenceJet;
use crate::hyper::thermoelastic::plane_stress::FiniteThermoelasticPlaneStressState;
use crate::sheet::{SheetError, SheetResponse};
use crate::visco_hyper::membrane::ViscoMembraneState;
use crate::visco_hyper::{ViscoHyperelastic, ViscoState};
use fs_ad::Dual64;

impl SheetExpansionReference {
    /// Exact local jet of the retained original engineering-length path in
    /// fixed sheet axes. Endpoints and slope-changing source knots refuse;
    /// isolated mean expansion evidence cannot manufacture this derivative.
    pub fn thermal_reference_jet(&self) -> Result<ThermalReferenceJet, SheetError> {
        let slopes = self
            .source()
            .at_temperature_with_smooth_tangent(self.temperature_k())?
            .principal_stretch_derivatives_per_k();
        let q = match self.frame() {
            SheetExpansionFrame::Isotropic => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            SheetExpansionFrame::FixedPrincipal {
                principal_to_sheet, ..
            } => principal_to_sheet,
        };
        let first = core::array::from_fn(|ij| {
            (0..3)
                .map(|k| q[ij / 3][k] * slopes[k] * q[ij % 3][k])
                .sum()
        });
        Ok(ThermalReferenceJet::try_new(
            *self.reference_deformation(),
            first,
            [0.0; 9],
        )?)
    }
}

/// One source-retaining, locally equilibrated finite membrane response.
/// Engineering strain is `[E11,E22,2E12]`. All storage and entropy are per
/// original reference area; current thickness/area do not rescale them.
#[derive(Clone, Debug)]
pub struct ThermoelasticMembraneState<'a> {
    source: &'a HyperelasticMembrane,
    strain: [f64; 3],
    native: FiniteThermoelasticPlaneStressState,
    response: SheetResponse,
    resultant_temperature_derivative_n_m_k: [f64; 3],
    viscoelastic_history: Option<(&'a ViscoHyperelastic, &'a ViscoState)>,
}

impl HyperelasticMembrane {
    /// Reduce the original thermal natural map and full compressible energy
    /// at its actual source temperature. The same native plane-stress owner
    /// relaxes thickness and both transverse shears once. Exact volume needs
    /// a separate constrained thermal potential and is explicitly refused.
    ///
    /// Source and model temperature support must be two-sided. Native elastic
    /// activation conventions are retained; this does not certify a two-sided
    /// mechanical derivative at a fibre activation surface. Temperature-varying
    /// elastic coefficients, viscous history and caloric storage are not inferred.
    pub fn thermoelastic_state(
        &self,
        strain: [f64; 3],
    ) -> Result<ThermoelasticMembraneState<'_>, SheetError> {
        let (in_plane, temperature, jet) = self.thermal_data(strain)?;
        let native = FiniteThermoelasticPlaneStressState::try_new(
            self.model().law(),
            &in_plane,
            temperature,
            jet,
            self.settings(),
        )?;
        ThermoelasticMembraneState::from_native(self, strain, native, None)
    }

    /// The same original sheet source at fixed accepted material-frame viscous
    /// metrics. The explicit spectrum must retain this exact equilibrium law
    /// and natural reference; it acquires no additional source qualification.
    /// The full native physical potential supplies every branch contribution
    /// to entropy, internal energy, capacity and thermal stress. No relaxation
    /// interval passes, and no accepted metric or energy ledger is changed.
    pub fn frozen_visco_thermoelastic_state<'a>(
        &'a self,
        law: &'a ViscoHyperelastic,
        accepted: &'a ViscoMembraneState,
        strain: [f64; 3],
    ) -> Result<ThermoelasticMembraneState<'a>, SheetError> {
        let (in_plane, temperature, jet) = self.thermal_data(strain)?;
        let ViscoMembraneState::Full(history) = accepted else {
            return Err(SheetError::Invalid(
                "exact-volume viscous membrane thermodynamics requires its constrained thermal owner",
            ));
        };
        if &law.equilibrium != self.model().law() || law.reference_deformation() != jet.reference()
        {
            return Err(SheetError::Invalid(
                "viscous membrane thermodynamics must retain its original sheet equilibrium and natural reference",
            ));
        }
        let native = law.frozen_thermoelastic_plane_stress_state(
            history,
            &in_plane,
            temperature,
            jet,
            None,
            self.settings(),
        )?;
        ThermoelasticMembraneState::from_native(self, strain, native, Some((law, history)))
    }

    fn thermal_data(
        &self,
        strain: [f64; 3],
    ) -> Result<([f64; 4], f64, ThermalReferenceJet), SheetError> {
        if self.volumetric_treatment() != VolumetricTreatment::OriginalMaterial {
            return Err(SheetError::Invalid(
                "exact-volume membrane thermodynamics requires its constrained thermal owner",
            ));
        }
        let expansion = self.free_expansion_reference().ok_or(SheetError::Invalid(
            "membrane thermodynamics requires the retained original expansion path",
        ))?;
        let temperature = expansion.temperature_k();
        if self
            .model()
            .model_card()
            .validity
            .bound("T")
            .is_some_and(|(low, high)| temperature <= low || temperature >= high)
        {
            return Err(SheetError::Invalid(
                "membrane thermodynamics requires a two-sided model temperature neighborhood",
            ));
        }
        Ok((
            admitted_representative(strain)?,
            temperature,
            expansion.thermal_reference_jet()?,
        ))
    }
}

impl<'a> ThermoelasticMembraneState<'a> {
    fn from_native(
        source: &'a HyperelasticMembrane,
        strain: [f64; 3],
        native: FiniteThermoelasticPlaneStressState,
        viscoelastic_history: Option<(&'a ViscoHyperelastic, &'a ViscoState)>,
    ) -> Result<Self, SheetError> {
        let thickness = source.reference_thickness_m();
        let response = response_from_plane_stress(strain, *native.response(), thickness)?;
        // This is the SAME objective kinematic map as the mechanical sheet.
        // Temperature does not change this map at held total Green strain.
        let variables = core::array::from_fn(|i| Dual64::<3>::variable(strain[i], i));
        let representative = representative(variables);
        let resultant_temperature_derivative_n_m_k = core::array::from_fn(|i| {
            thickness
                * (0..4)
                    .map(|a| {
                        representative[a].eps[i] * native.stress_temperature_derivative_pa_k()[a]
                    })
                    .sum::<f64>()
        });
        if [
            thickness * native.entropy_j_m3_k(),
            thickness * native.internal_energy_j_m3(),
            thickness * native.mechanical_heat_capacity_j_m3_k(),
        ]
        .iter()
        .chain(&resultant_temperature_derivative_n_m_k)
        .any(|value| !value.is_finite())
        {
            return Err(SheetError::Invalid(
                "unrepresentable original-area membrane thermodynamics",
            ));
        }
        Ok(Self {
            source,
            strain,
            native,
            response,
            resultant_temperature_derivative_n_m_k,
            viscoelastic_history,
        })
    }

    pub const fn source(&self) -> &'a HyperelasticMembrane {
        self.source
    }
    pub const fn strain(&self) -> &[f64; 3] {
        &self.strain
    }
    pub const fn native(&self) -> &FiniteThermoelasticPlaneStressState {
        &self.native
    }
    pub const fn response(&self) -> &SheetResponse {
        &self.response
    }
    /// Original explicit law and immutable accepted history, when present.
    pub const fn viscoelastic_history(&self) -> Option<(&'a ViscoHyperelastic, &'a ViscoState)> {
        self.viscoelastic_history
    }
    pub fn entropy_j_m2_k(&self) -> f64 {
        self.source.reference_thickness_m() * self.native.entropy_j_m3_k()
    }
    pub fn internal_energy_j_m2(&self) -> f64 {
        self.source.reference_thickness_m() * self.native.internal_energy_j_m3()
    }
    /// Mechanical contribution only. Add the explicit caloric reference once.
    pub fn mechanical_heat_capacity_j_m2_k(&self) -> f64 {
        self.source.reference_thickness_m() * self.native.mechanical_heat_capacity_j_m3_k()
    }
    pub const fn resultant_temperature_derivative_n_m_k(&self) -> &[f64; 3] {
        &self.resultant_temperature_derivative_n_m_k
    }

    /// Eliminate temperature at constant TOTAL entropy, using an explicitly
    /// supplied caloric capacity per original area. The mechanical correction
    /// is added once, after transverse relaxation. The caller owns the caloric
    /// reference meaning and its source/temperature admission.
    pub fn adiabatic_response(
        &self,
        caloric_heat_capacity_j_m2_k: f64,
    ) -> Result<AdiabaticMembraneResponse, SheetError> {
        let total = caloric_heat_capacity_j_m2_k + self.mechanical_heat_capacity_j_m2_k();
        if !caloric_heat_capacity_j_m2_k.is_finite()
            || caloric_heat_capacity_j_m2_k <= 0.0
            || !total.is_finite()
            || total <= 0.0
        {
            return Err(SheetError::Invalid(
                "positive finite caloric and total fixed-strain membrane capacities required",
            ));
        }
        let scale = self.native.full_state().temperature_k() / total;
        let beta = self.resultant_temperature_derivative_n_m_k;
        let temperature_gradient_k = beta.map(|value| scale * value);
        let tangent_n_m = core::array::from_fn(|ij| {
            temperature_gradient_k[ij / 3].mul_add(beta[ij % 3], self.response.tangent_n_m[ij])
        });
        if temperature_gradient_k
            .iter()
            .chain(&tangent_n_m)
            .any(|value| !value.is_finite())
        {
            return Err(SheetError::Invalid(
                "unrepresentable entropy-constrained membrane response",
            ));
        }
        Ok(AdiabaticMembraneResponse {
            fixed_strain_heat_capacity_j_m2_k: total,
            temperature_gradient_k,
            tangent_n_m,
        })
    }
}

/// Local entropy-constrained constitutive response, excluding spatial geometric
/// stiffness, inertia, heat transport and irreversible relaxation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdiabaticMembraneResponse {
    pub fixed_strain_heat_capacity_j_m2_k: f64,
    /// `dT/d[E11,E22,2E12]` at fixed total entropy.
    pub temperature_gradient_k: [f64; 3],
    /// Original-area material tangent in engineering Green coordinates.
    pub tangent_n_m: [f64; 9],
}
