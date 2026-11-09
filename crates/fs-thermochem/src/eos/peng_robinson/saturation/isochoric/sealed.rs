//! One fixed-density source through stable liquid, coexistence and vapor.
//!
//! Phase selection is the original convex saturation-volume comparison followed
//! by at most one original homogeneous density evaluation. Direct U/S use the
//! parent phase integrals, also when a phase appears or disappears. There is no
//! new flash, UV inversion, latent regularization or caloric source.

use super::super::super::PengRobinsonDensityEvaluationV1;
use super::super::PengRobinsonSaturationV1;
use super::{
    Dd, Density, E, Nasa9StandardStateModelV1, PengRobinsonMixtureV1,
    PengRobinsonSaturatedMixtureV1, PhaseCalorics, R, Temperature, admit_source, poll, positive,
    source_changes,
};
use crate::eos::IdealGasMolarVolumeQuantityV1 as MolarVolume;
use fs_qty::Pressure;

/// Original phase allocation. Homogeneous states are outside the resolved
/// saturation envelope; mechanical stability alone never selects them.
#[derive(Clone, Debug, PartialEq)]
pub enum PengRobinsonIsochoricPhaseV1 {
    Liquid(PengRobinsonDensityEvaluationV1),
    Coexistence(PengRobinsonSaturatedMixtureV1),
    Vapor(PengRobinsonDensityEvaluationV1),
}

/// Source-bound original density equilibrium, including its saturation witness.
/// Only the original EOS constructor below can create this value.
#[derive(Clone, Debug, PartialEq)]
pub struct PengRobinsonEquilibriumDensityV1 {
    saturation: PengRobinsonSaturationV1,
    density: Density,
    phase: PengRobinsonIsochoricPhaseV1,
}
impl PengRobinsonEquilibriumDensityV1 {
    pub const fn saturation(&self) -> &PengRobinsonSaturationV1 {
        &self.saturation
    }
    pub const fn density(&self) -> Density {
        self.density
    }
    pub const fn phase(&self) -> &PengRobinsonIsochoricPhaseV1 {
        &self.phase
    }
    pub fn temperature(&self) -> Temperature {
        self.saturation.temperature()
    }
    pub fn pressure(&self) -> Pressure {
        match &self.phase {
            PengRobinsonIsochoricPhaseV1::Coexistence(_) => self.saturation.pressure(),
            PengRobinsonIsochoricPhaseV1::Liquid(s) | PengRobinsonIsochoricPhaseV1::Vapor(s) => {
                s.pressure()
            }
        }
    }
    /// Additional EOS evaluations after the supplied saturation solve.
    pub fn homogeneous_evaluations(&self) -> usize {
        usize::from(!matches!(
            self.phase,
            PengRobinsonIsochoricPhaseV1::Coexistence(_)
        ))
    }
    pub fn vapor_fraction(&self) -> f64 {
        match &self.phase {
            PengRobinsonIsochoricPhaseV1::Liquid(_) => 0.0,
            PengRobinsonIsochoricPhaseV1::Coexistence(s) => s.vapor_fraction(),
            PengRobinsonIsochoricPhaseV1::Vapor(_) => 1.0,
        }
    }
    /// Finite stepping retains the one-sided coexistence tangent here. A
    /// two-sided harmonic/linear-response derivative does not exist.
    pub fn is_phase_boundary(&self) -> bool {
        match &self.phase {
            PengRobinsonIsochoricPhaseV1::Coexistence(s) => {
                let v = s.molar_volume().value();
                s.vapor_fraction() == 0.0
                    || s.vapor_fraction() == 1.0
                    || v <= self.saturation.liquid().molar_volume_m3_per_mol()
                    || v >= self.saturation.vapor().molar_volume_m3_per_mol()
            }
            _ => false,
        }
    }
    fn calorics(&self) -> PhaseCalorics<'_> {
        let (root, quality) = match &self.phase {
            PengRobinsonIsochoricPhaseV1::Coexistence(s) => return PhaseCalorics::coexistence(s),
            PengRobinsonIsochoricPhaseV1::Liquid(s) => (s.root(), 0.0),
            PengRobinsonIsochoricPhaseV1::Vapor(s) => (s.root(), 1.0),
        };
        // The absent phase contributes no distinct latent datum. Duplicating
        // this one root makes the same symmetric amount identity exact.
        PhaseCalorics {
            temperature: self.temperature(),
            roots: [root, root],
            quality,
        }
    }
    fn molar_volume(&self) -> f64 {
        match &self.phase {
            PengRobinsonIsochoricPhaseV1::Coexistence(s) => s.molar_volume().value(),
            PengRobinsonIsochoricPhaseV1::Liquid(s) | PengRobinsonIsochoricPhaseV1::Vapor(s) => {
                s.root().molar_volume_m3_per_mol()
            }
        }
    }
}

impl PengRobinsonMixtureV1 {
    /// Reuse an ORIGINAL accepted saturation state to resolve fixed-density
    /// equilibrium. Its source component must match exactly. No additional
    /// saturation solve occurs; outside the dome one original homogeneous
    /// density evaluation supplies pressure and calorics. Critical/saturation
    /// failures must be resolved by that owner before this method is called.
    pub fn pure_equilibrium_density_at_saturation(
        &self,
        saturation: PengRobinsonSaturationV1,
        density: Density,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PengRobinsonEquilibriumDensityV1, E> {
        poll(&mut cancelled)?;
        if self.components().len() != 1
            || self.mole_fractions() != [1.0]
            || self.components()[0] != *saturation.component()
            || !positive(density.value())
        {
            return Err(E::Invalid(
                "same pure PR source and positive conserved density required",
            ));
        }
        let volume = self.molar_mass().value() / density.value();
        if !positive(volume) {
            return Err(E::Invalid("fixed molar volume is not representable"));
        }
        let vl = saturation.liquid().molar_volume_m3_per_mol();
        let vv = saturation.vapor().molar_volume_m3_per_mol();
        let phase = if volume >= vl && volume <= vv {
            PengRobinsonIsochoricPhaseV1::Coexistence(
                saturation.at_molar_volume(MolarVolume::new(volume))?,
            )
        } else {
            poll(&mut cancelled)?;
            let homogeneous =
                self.evaluate_density_temperature(saturation.temperature(), density)?;
            if volume < vl {
                PengRobinsonIsochoricPhaseV1::Liquid(homogeneous)
            } else {
                PengRobinsonIsochoricPhaseV1::Vapor(homogeneous)
            }
        };
        poll(&mut cancelled)?;
        Ok(PengRobinsonEquilibriumDensityV1 {
            saturation,
            density,
            phase,
        })
    }
}

/// Actual fixed-volume caloric/pressure tangent. At a phase boundary it is the
/// one-sided coexistence tangent, never a fabricated common acoustic response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PureFluidSealedResponseV1 {
    pub molar_cv_j_mol_k: f64,
    pub dp_dt_v_pa_k: f64,
    pub ds_dt_v_j_mol_k2: f64,
    pub phase_boundary: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PureFluidSealedChangeV1 {
    pub specific_internal_energy_change_j_kg: f64,
    pub source_entropy_change_j_kg_k: f64,
    pub response: PureFluidSealedResponseV1,
    pub relative_volume_residual: f64,
    /// Zero for the identical state; two for one active homogeneous phase;
    /// otherwise four original scalar phase caloric evaluations.
    pub caloric_evaluations: usize,
}

/// Immutable original NASA/PR calorics at fixed density across the dome.
/// Every endpoint is an admitted native density equilibrium. The same analytic
/// ideal-Cv floor applies to homogeneous residual Cv and moving-phase Cv on
/// the original kappa>=0 subcritical branch. Saturation and UV search budgets
/// remain with their existing owners.
pub struct PengRobinsonSealedCaloricReferenceV1<'a> {
    model: &'a PengRobinsonMixtureV1,
    source: &'a Nasa9StandardStateModelV1,
    initial: PengRobinsonEquilibriumDensityV1,
    volume: f64,
    domain: [Temperature; 2],
}
impl<'a> PengRobinsonSealedCaloricReferenceV1<'a> {
    pub fn new(
        model: &'a PengRobinsonMixtureV1,
        source: &'a Nasa9StandardStateModelV1,
        initial: PengRobinsonEquilibriumDensityV1,
        domain: [Temperature; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, E> {
        poll(&mut cancelled)?;
        let volume = admit_source(model, source, initial.density(), domain, &mut cancelled)?;
        let result = Self {
            model,
            source,
            initial,
            volume,
            domain,
        };
        result.admit(&result.initial)?;
        poll(&mut cancelled)?;
        Ok(result)
    }
    pub const fn source(&self) -> &'a Nasa9StandardStateModelV1 {
        self.source
    }
    pub const fn constitutive_law(&self) -> &'a PengRobinsonMixtureV1 {
        self.model
    }
    pub const fn initial_state(&self) -> &PengRobinsonEquilibriumDensityV1 {
        &self.initial
    }
    pub fn density(&self) -> Density {
        self.initial.density()
    }
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.domain
    }
    pub fn response_at(
        &self,
        state: &PengRobinsonEquilibriumDensityV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidSealedResponseV1, E> {
        poll(&mut cancelled)?;
        let (_, response) = self.admit(state)?;
        poll(&mut cancelled)?;
        Ok(response)
    }
    /// Exact bounded algebra cost used by an outer whole-call EOS work budget.
    /// It does not admit an endpoint; `change_to` still checks the full binding.
    pub fn caloric_evaluations_to(&self, state: &PengRobinsonEquilibriumDensityV1) -> usize {
        if state == &self.initial {
            0
        } else {
            self.initial.calorics().work_to(&state.calorics())
        }
    }
    pub fn change_to(
        &self,
        state: &PengRobinsonEquilibriumDensityV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidSealedChangeV1, E> {
        poll(&mut cancelled)?;
        let (relative_volume_residual, response) = self.admit(state)?;
        let (
            specific_internal_energy_change_j_kg,
            source_entropy_change_j_kg_k,
            caloric_evaluations,
        ) = if state == &self.initial {
            (0.0, 0.0, 0)
        } else {
            source_changes(
                self.model,
                self.source,
                self.initial.calorics(),
                state.calorics(),
                &mut cancelled,
            )?
        };
        poll(&mut cancelled)?;
        Ok(PureFluidSealedChangeV1 {
            specific_internal_energy_change_j_kg,
            source_entropy_change_j_kg_k,
            response,
            relative_volume_residual,
            caloric_evaluations,
        })
    }
    fn admit(
        &self,
        state: &PengRobinsonEquilibriumDensityV1,
    ) -> Result<(f64, PureFluidSealedResponseV1), E> {
        let t = state.temperature().value();
        if state.saturation().component() != &self.model.components()[0]
            || state.density() != self.density()
            || t < self.domain[0].value()
            || t > self.domain[1].value()
        {
            return Err(E::Invalid(
                "same original component, exact fixed density and source interval required",
            ));
        }
        let relative = ((Dd::from_f64(state.molar_volume()) - Dd::from_f64(self.volume))
            / Dd::from_f64(self.volume))
        .to_f64();
        if !relative.is_finite() || relative.abs() > 128.0 * f64::EPSILON {
            return Err(E::Invalid(
                "original phase volume does not close the conserved density",
            ));
        }
        let ideal = self
            .source
            .evaluate_properties(state.temperature())?
            .cp()
            .value()
            - R;
        let (residual, pt) = match state.phase() {
            PengRobinsonIsochoricPhaseV1::Coexistence(s) => {
                for root in [s.saturation().liquid(), s.saturation().vapor()] {
                    if !positive(ideal + root.residual_cv_j_mol_k()) {
                        return Err(E::Invalid(
                            "both original phase capacities must be positive",
                        ));
                    }
                }
                (
                    s.equilibrium_isochoric_residual_cv_j_mol_k()?,
                    s.saturation().saturation_pressure_slope_pa_k(),
                )
            }
            PengRobinsonIsochoricPhaseV1::Liquid(s) | PengRobinsonIsochoricPhaseV1::Vapor(s) => (
                s.root().residual_cv_j_mol_k(),
                s.root().pressure_derivatives().0,
            ),
        };
        let cv = ideal + residual;
        let dsdt = cv / t;
        if !positive(ideal) || !positive(cv) || cv < ideal || !positive(pt) || !positive(dsdt) {
            return Err(E::Invalid(
                "original fixed-volume response violates its stable Cv floor",
            ));
        }
        Ok((
            relative,
            PureFluidSealedResponseV1 {
                molar_cv_j_mol_k: cv,
                dp_dt_v_pa_k: pt,
                ds_dt_v_j_mol_k2: dsdt,
                phase_boundary: state.is_phase_boundary(),
            },
        ))
    }
}
