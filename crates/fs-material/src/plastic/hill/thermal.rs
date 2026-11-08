//! Source-bound heating of the original fixed-axis Hill return map.
//!
//! The anisotropic metric and material frame stay fixed. The reference yield
//! radius may vary with temperature; elastic and energetic isotropic/Prager
//! hardening moduli are affine, with fixed Poisson ratio and reference density.
//! Explicitly selected Hill Voce saturation has affine Q(T) and fixed b, with
//! g(alpha)=alpha-(1-exp(-b*alpha))/b. Its additional free energy is Q*g,
//! entropy is -Q'*g, and no frozen-history heat-capacity term is introduced.
//! Optional sourced isotropic expansion is constant. These restrictions keep the retained Hill
//! multiplier meaningful and the mechanical entropy differentiable.
//!
//! With e = eps - eps_p - a(T-T_ref)I, the free energy is
//! psi = e:C_el(T):e/2 + H(T)alpha²/2 + C_k(T)eps_p:eps_p/2 + Q(T)g(alpha).
//! Its mechanical entropy is
//! s = a tr(sigma) - (E'/E)psi_elastic - H'alpha²/2
//!     - C_k'eps_p:eps_p/2 - Q'g(alpha).
//! Prager hardening requires an explicit independent source selection. Its
//! retained coordinate is plastic strain, so changing C_k, including reaching
//! zero, never divides by the modulus or discards a backstress history.
//! Optional sourced linear overstress adds eta(T)*delta_alpha/dt to the
//! endpoint flow resistance. Its yield and viscous work are physical heat;
//! viscosity stores no energy and adds no mechanical entropy. Rate histories
//! may remain outside the static yield surface and relax during a held strain.
//! Native endpoint yield work is physical heat. Coefficient storage and delta(Ts)
//! give the separate reversible caloric transfer; numerical loss is never heat.
//! The existing thermal section/shell owners perform reduction, equilibrium,
//! dynamics, contact and heat transport without projecting Hill history to J2.

mod source;
pub use source::ThermalHillSourceQuery;

use super::rate::{self, LinearOverstressHill48};
use super::source::{HILL_OVERSTRESS_VISCOSITY_PROPERTY, HillSourceError, ResolvedHill48};
use super::{HillError, HillState};
use crate::phase::{EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve};
use crate::plastic::energy::EnergyTolerances;
use crate::plastic::section::thermal::source::{
    ThermalSectionIncrement, ThermalSectionSourceError, ThermalSectionSourceFamily,
    ThermalSectionStep,
};
use crate::plastic::voce::hardening;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ScalarPropertyRequirement,
};
use crate::tensor::{contract, frobenius_norm, trace};
use crate::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;
use fs_matdb::{MaterialCard, QueryPoint};

/// Original source, caloric support, history or thermodynamic refusal.
#[derive(Debug)]
pub enum ThermalHillError {
    Invalid(&'static str),
    Source(HillSourceError),
    Caloric(PhaseStateError),
    Transport(ThermalMaterialError),
    Mechanical(HillError),
}
impl core::fmt::Display for ThermalHillError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "thermal Hill48: {s}"),
            Self::Source(e) => e.fmt(f),
            Self::Caloric(e) => e.fmt(f),
            Self::Transport(e) => e.fmt(f),
            Self::Mechanical(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for ThermalHillError {}
impl From<HillSourceError> for ThermalHillError {
    fn from(e: HillSourceError) -> Self {
        Self::Source(e)
    }
}
impl From<MaterialStatePointError> for ThermalHillError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Source(e.into())
    }
}
impl From<PhaseStateError> for ThermalHillError {
    fn from(e: PhaseStateError) -> Self {
        Self::Caloric(e)
    }
}
impl From<ThermalMaterialError> for ThermalHillError {
    fn from(e: ThermalMaterialError) -> Self {
        Self::Transport(e)
    }
}
impl From<HillError> for ThermalHillError {
    fn from(e: HillError) -> Self {
        Self::Mechanical(e)
    }
}

/// One immutable source binding; every layer resolves its endpoint law once.
/// Histories remain local-axis Hill states at every native quadrature point.
#[derive(Clone, Debug)]
pub struct ThermalHillSource<'a> {
    card: &'a MaterialCard,
    caloric: &'a ResolvedHeatCapacityCurve,
    context: QueryPoint,
    temperature_axis: String,
    requirements: Vec<ScalarPropertyRequirement>,
    selection: MaterialPropertySelection,
    reference: ResolvedHill48,
    reference_temperature_k: f64,
    density_kg_m3: f64,
    expansion_per_k: f64,
    youngs_slope_pa_k: f64,
    hardening_slope_pa_k: f64,
    /// None means no kinematic source was selected; Some(0) is an explicitly
    /// selected constant modulus. The original source bundle owns its value.
    kinematic_slope_pa_k: Option<f64>,
    /// Explicit saturation selection, including Q=0 and Q'=0. The reference
    /// law retains b per Hill work-conjugate alpha and its original source.
    voce_slope_pa_k: Option<f64>,
    /// Explicitly pinned dissipative viscosity; never a stored-energy modulus.
    linear_overstress_selected: bool,
}

/// The original oriented Hill law at one endpoint, with stress-free expansion.
/// No thermal transport, extra constitutive step or per-point source query.
#[derive(Clone, Debug)]
pub struct ThermalHill {
    source: ResolvedHill48,
    expansion_strain: f64,
    rate: Option<LinearOverstressHill48>,
}
impl ThermalHill {
    /// Original claim identities, conditions, orientation and endpoint values.
    pub const fn resolved(&self) -> &ResolvedHill48 {
        &self.source
    }

    /// Prepared oriented overstress law at the actual physical interval.
    /// None retains the original rate-independent Hill response exactly.
    pub const fn linear_overstress(&self) -> Option<&LinearOverstressHill48> {
        self.rate.as_ref()
    }

    fn mechanical_strain(&self, total: &Voigt) -> Result<Voigt, HillError> {
        if total.iter().any(|v| !v.is_finite())
            || frobenius_norm(total) > self.source.local_law().elastic().strain_limit
        {
            return Err(HillError::State(
                "total thermal strain exceeds the small-strain domain",
            ));
        }
        Ok(core::array::from_fn(|i| {
            total[i] - if i < 3 { self.expansion_strain } else { 0.0 }
        }))
    }
}
impl SmallStrainLaw for ThermalHill {
    type State = HillState;
    fn initial_state(&self) -> HillState {
        self.source.initial_state()
    }
    fn stress(&self, strain: &Voigt, old: &HillState) -> Voigt {
        self.mechanical_strain(strain).map_or([f64::NAN; 6], |e| {
            self.rate
                .as_ref()
                .map_or_else(|| self.source.stress(&e, old), |law| law.stress(&e, old))
        })
    }
    fn tangent(&self, strain: &Voigt, old: &HillState) -> Tangent6 {
        self.mechanical_strain(strain)
            .map_or([[f64::NAN; 6]; 6], |e| {
                self.rate
                    .as_ref()
                    .map_or_else(|| self.source.tangent(&e, old), |law| law.tangent(&e, old))
            })
    }
    fn update_state(&self, strain: &Voigt, old: &HillState) -> HillState {
        self.mechanical_strain(strain).map_or_else(
            |_| HillState {
                plastic_strain: [f64::NAN; 6],
                alpha: f64::NAN,
            },
            |e| {
                self.rate.as_ref().map_or_else(
                    || self.source.update_state(&e, old),
                    |law| law.update_state(&e, old),
                )
            },
        )
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.rate
            .as_ref()
            .map_or_else(|| self.source.admissibility(), |law| law.admissibility())
    }
    fn card(&self) -> ModelCard {
        self.rate
            .as_ref()
            .map_or_else(|| self.source.card(), |law| law.card())
    }
}

impl ThermalHillSource<'_> {
    fn phase(&self, h: f64) -> Result<EquilibriumPhaseState, ThermalHillError> {
        Ok(self.caloric.curve().state_at_specific_enthalpy(h)?)
    }

    fn endpoint(&self, h: f64) -> Result<(EquilibriumPhaseState, ThermalHill), ThermalHillError> {
        let phase = self.phase(h)?;
        let source = self.resolve_temperature(phase.temperature_k())?;
        let expansion_strain = product(
            self.expansion_per_k,
            phase.temperature_k() - self.reference_temperature_k,
        )?;
        Ok((
            phase,
            ThermalHill {
                source,
                expansion_strain,
                rate: None,
            },
        ))
    }

    fn endpoint_for_step(
        &self,
        h: f64,
        duration_s: f64,
    ) -> Result<(EquilibriumPhaseState, ThermalHill), ThermalHillError> {
        duration(duration_s)?;
        let (phase, mut endpoint) = self.endpoint(h)?;
        if self.linear_overstress_selected {
            let viscosity = endpoint
                .source
                .material()
                .property(HILL_OVERSTRESS_VISCOSITY_PROPERTY)
                .ok_or(ThermalHillError::Invalid(
                    "selected Hill viscosity is absent",
                ))?
                .value_si();
            endpoint.rate = Some(LinearOverstressHill48::new_oriented(
                endpoint.source.law().clone(),
                viscosity,
                duration_s,
            )?);
        }
        Ok((phase, endpoint))
    }

    #[allow(clippy::too_many_arguments)]
    fn increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &HillState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<ThermalSectionIncrement<HillState>, ThermalHillError> {
        duration(duration_s)?;
        let (old_phase, previous) = self.endpoint(old_h)?;
        let (phase, current) = self.endpoint_for_step(new_h, duration_s)?;
        let old_mechanical = previous.mechanical_strain(old_strain)?;
        let mechanical = current.mechanical_strain(strain)?;
        // Validate the old history at its actual temperature. A selected rate
        // law admits physical overstress; it must not project the old endpoint
        // onto the static surface or invent a second constitutive interval.
        let native = if let Some(rate) = &current.rate {
            rate.advance_with_energy_from(
                previous.source.law(),
                &old_mechanical,
                old,
                &mechanical,
                tolerances,
            )?
        } else {
            current.source.law().advance_with_energy_from(
                previous.source.law(),
                &old_mechanical,
                old,
                &mechanical,
                tolerances,
            )?
        };
        let energy = &native.energy;
        let old_law = previous.source.local_law();
        let law = current.source.local_law();
        let old_plastic = previous.source.law().plastic_strain_spatial(old);
        let old_elastic = core::array::from_fn(|i| old_mechanical[i] - old_plastic[i]);
        // Elastic reconstruction inspects accepted history; a second return
        // would erase the thermal release that this interval must account for.
        let old_stress = old_law.elastic().stress(&old_elastic, &());
        let old_square = product(old.alpha, old.alpha)?;
        let new_square = product(native.response.state.alpha, native.response.state.alpha)?;
        let elastic_derivative = self.youngs_slope_pa_k / law.elastic().youngs;
        let entropy = |stress: &Voigt, elastic_energy: f64, alpha_square: f64| {
            finite(
                product(self.expansion_per_k, trace(stress))?
                    - product(elastic_derivative, elastic_energy)?
                    - product(0.5 * self.hardening_slope_pa_k, alpha_square)?,
            )
        };
        // Both native energies use E(T1); E'/E1 times native initial energy
        // equals E'/E0 times actual initial energy, including E'(T)=0.
        let mut old_entropy = entropy(&old_stress, energy.initial.elastic_j_m3, old_square)?;
        let mut new_entropy = entropy(
            &native.response.stress_pa,
            energy.final_state.elastic_j_m3,
            new_square,
        )?;
        let mut coefficient_storage_change_j_m3 = finite(
            product(
                (law.elastic().youngs - old_law.elastic().youngs) / law.elastic().youngs,
                energy.initial.elastic_j_m3,
            )? + product(
                0.5 * (law.hardening_modulus_pa() - old_law.hardening_modulus_pa()),
                old_square,
            )?,
        )?;
        let kinematic_slope = self.kinematic_slope_pa_k.unwrap_or(0.0);
        if kinematic_slope != 0.0 {
            // Plastic strain is already the retained Prager coordinate in
            // fixed material axes. Differentiate the native C*eps_p:eps_p/2
            // directly; neither C=0 nor a stress-free endpoint loses history.
            let old_shape = product(0.5, contract(&old.plastic_strain, &old.plastic_strain))?;
            let new_shape = product(
                0.5,
                contract(
                    &native.response.state.plastic_strain,
                    &native.response.state.plastic_strain,
                ),
            )?;
            old_entropy = finite(old_entropy - product(kinematic_slope, old_shape)?)?;
            new_entropy = finite(new_entropy - product(kinematic_slope, new_shape)?)?;
            coefficient_storage_change_j_m3 = finite(
                coefficient_storage_change_j_m3
                    + product(
                        law.kinematic_modulus_pa() - old_law.kinematic_modulus_pa(),
                        old_shape,
                    )?,
            )?;
        }
        let voce_slope = self.voce_slope_pa_k.unwrap_or(0.0);
        if voce_slope != 0.0 {
            let rate = law.voce_rate_per_strain();
            // Evaluate the signed coefficient inside the shared stable
            // primitive. A unit-Q shape could underflow for tiny b while
            // Q'*g or delta_Q*g remains representable. Q=0 does not erase
            // entropy from a nonzero admitted temperature derivative.
            old_entropy = finite(old_entropy + hardening::storage(-voce_slope, rate, old.alpha))?;
            new_entropy = finite(
                new_entropy + hardening::storage(-voce_slope, rate, native.response.state.alpha),
            )?;
            coefficient_storage_change_j_m3 = finite(
                coefficient_storage_change_j_m3
                    + hardening::storage(
                        law.voce_saturation_pa() - old_law.voce_saturation_pa(),
                        rate,
                        old.alpha,
                    ),
            )?;
        }
        let thermoelastic_storage_change_j_m3 = finite(
            product(phase.temperature_k(), new_entropy)?
                - product(old_phase.temperature_k(), old_entropy)?,
        )?;
        let reversible_heat_j_m3 = if self.youngs_slope_pa_k == 0.0
            && self.hardening_slope_pa_k == 0.0
            && kinematic_slope == 0.0
            && voce_slope == 0.0
        {
            // Exact constant-modulus reduction avoids subtracting large T*s terms.
            product(
                -old_phase.temperature_k(),
                finite(new_entropy - old_entropy)?,
            )?
        } else {
            finite(
                product(
                    trace(&native.response.stress_pa),
                    product(
                        self.expansion_per_k,
                        phase.temperature_k() - old_phase.temperature_k(),
                    )?,
                )? - coefficient_storage_change_j_m3
                    - thermoelastic_storage_change_j_m3,
            )?
        };
        let stored_change_j_m3 =
            finite(energy.stored_change_j_m3() + coefficient_storage_change_j_m3)?;
        let delta = core::array::from_fn(|i| strain[i] - old_strain[i]);
        let mechanical_work_j_m3 = finite(contract(&native.response.stress_pa, &delta))?;
        let residual = finite(
            mechanical_work_j_m3
                - stored_change_j_m3
                - thermoelastic_storage_change_j_m3
                - energy.physical_dissipation_j_m3
                - reversible_heat_j_m3
                - energy.numerical_dissipation_j_m3,
        )?;
        if residual.abs() > tolerances.energy_j_m3 {
            return Err(HillError::EnergyResidual {
                residual_j_m3: residual,
                tolerance_j_m3: tolerances.energy_j_m3,
            }
            .into());
        }
        Ok(ThermalSectionIncrement {
            state: native.response.state,
            stress_pa: native.response.stress_pa,
            physical_heat_j_m3: energy.physical_dissipation_j_m3,
            creep_heat_j_m3: 0.0,
            reversible_heat_j_m3,
            stored_change_j_m3,
            coefficient_storage_change_j_m3,
            thermoelastic_storage_change_j_m3,
            numerical_loss_j_m3: energy.numerical_dissipation_j_m3,
            mechanical_work_j_m3,
        })
    }
}

impl<'a> ThermalSectionSourceFamily<'a> for ThermalHillSource<'a> {
    type Material = ThermalHill;
    fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.caloric
    }
    fn reference_density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    fn thermal_material_point<'s>(&'s self) -> ThermalMaterialPoint<'s>
    where
        'a: 's,
    {
        ThermalMaterialPoint {
            material: self.card,
            context: &self.context,
            temperature_axis: &self.temperature_axis,
            liquid_fraction_axis: None,
        }
    }
    fn admit_section_enthalpy(
        &self,
        h: f64,
    ) -> Result<EquilibriumPhaseState, ThermalSectionSourceError> {
        Ok(self.endpoint(h)?.0)
    }
    fn resolve_section_step(
        &self,
        old_h: f64,
        new_h: f64,
        duration_s: f64,
    ) -> Result<ThermalSectionStep<ThermalHill>, ThermalSectionSourceError> {
        duration(duration_s)?;
        self.phase(old_h)?;
        let (phase, law) = self.endpoint_for_step(new_h, duration_s)?;
        Ok(ThermalSectionStep {
            temperature_k: phase.temperature_k(),
            law,
        })
    }
    fn section_increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &HillState,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<ThermalSectionIncrement<HillState>, ThermalSectionSourceError> {
        self.increment(
            old_h, new_h, old_strain, old, new_strain, tolerances, duration_s,
        )
        .map_err(Into::into)
    }
    fn validate_section_history(
        &self,
        h: f64,
        strain: &Voigt,
        state: &HillState,
        stress_tolerance_pa: f64,
    ) -> Result<(), ThermalSectionSourceError> {
        let (_, law) = self.endpoint(h)?;
        let mechanical = law
            .mechanical_strain(strain)
            .map_err(ThermalHillError::from)?;
        if self.linear_overstress_selected {
            if !stress_tolerance_pa.is_finite() || stress_tolerance_pa <= 0.0 {
                return Err(
                    ThermalHillError::Invalid("positive finite stress allowance required").into(),
                );
            }
            // No artificial duration is needed to inspect an accepted history.
            rate::validate_history(law.source.law(), &mechanical, state)
                .map_err(ThermalHillError::from)?;
        } else {
            law.source
                .law()
                .validate_committed_history(&mechanical, state, stress_tolerance_pa)
                .map_err(ThermalHillError::from)?;
        }
        Ok(())
    }
}

fn duration(value: f64) -> Result<(), ThermalHillError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(ThermalHillError::Invalid(
            "positive finite physical duration required",
        ))
    }
}
fn finite(value: f64) -> Result<f64, ThermalHillError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ThermalHillError::Invalid(
            "unrepresentable thermal Hill arithmetic",
        ))
    }
}
fn product(a: f64, b: f64) -> Result<f64, ThermalHillError> {
    let value = finite(a * b)?;
    if value == 0.0 && a != 0.0 && b != 0.0 {
        Err(ThermalHillError::Invalid(
            "nonzero thermal work underflowed",
        ))
    } else {
        Ok(value)
    }
}
