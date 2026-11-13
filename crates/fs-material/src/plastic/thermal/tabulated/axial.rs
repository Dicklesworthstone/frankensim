//! Owned thermal intervals on the exact free-sided tabulated J2 subspace.
//!
//! A prepared interval retains its original source choices, actual endpoint
//! temperatures and physical duration. Newton trials only evaluate its cached
//! native axial law against the same accepted history. Thermal expansion shifts
//! all three normal strains; neither transverse plastic flow nor Poisson
//! contraction is suppressed. The shared tensor thermal owner audits storage,
//! entropy, physical heat and numerical loss without another energy model.

use super::{ThermalTabulatedJ2, ThermalTabulatedSource, source_error};
use crate::fiber::Uniaxial;
use crate::plastic::energy::{EnergyError, EnergyTolerances, StoredEnergy};
use crate::plastic::reduction::{UniaxialJ2Response, UniaxialJ2State, UniaxialTabulated};
use crate::plastic::source::HardeningPropertyBasis;
use crate::plastic::tabulated::rate::axial::UniaxialTabulatedOverstress;
use crate::plastic::tabulated::thermal_source::ResolvedThermalCauchyTabulatedHardening;
use crate::plastic::thermal::expansion::{PreparedThermoplasticEnergy, ThermalHistoryContribution};
use crate::plastic::thermal::{
    FeedbackError, OverstressSelection, ThermoplasticIncrement, product,
};
use crate::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use crate::tensor::contract;
use crate::{MaterialAdmissibility, SmallStrainLaw, Voigt};
use fs_blake3::ContentHash;
use fs_evidence::ModelCard;
use fs_matdb::QueryPoint;
use std::sync::Arc;

/// Accepted permanent strain, accumulated plastic flow and actual temperature.
/// Preparing a new physical interval never changes this retained history.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialThermalTabulatedState {
    pub mechanical: UniaxialJ2State,
    temperature_k: f64,
}

impl UniaxialThermalTabulatedState {
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
}

/// Uncommitted scalar response with its full compatible TOTAL tensor strain.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialThermalTabulatedResponse {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialThermalTabulatedState,
}

/// Exact native tensor thermal account paired with its retained axial state.
#[derive(Clone, Debug)]
pub struct UniaxialThermalTabulatedIncrement {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialThermalTabulatedState,
    pub thermal: ThermoplasticIncrement,
}

#[derive(Clone, Debug)]
enum Mechanical {
    Independent(UniaxialTabulated),
    Overstress(UniaxialTabulatedOverstress),
}

impl Mechanical {
    fn response(
        &self,
        strain: f64,
        state: &UniaxialJ2State,
    ) -> Result<UniaxialJ2Response, FeedbackError> {
        match self {
            Self::Independent(law) => law.response(strain, state),
            Self::Overstress(law) => law.response(strain, state),
        }
        .map_err(source_error)
    }

    fn initial_state(&self) -> UniaxialJ2State {
        match self {
            Self::Independent(law) => law.initial_state(),
            Self::Overstress(law) => law.initial_state(),
        }
    }
}

// These immutable choices describe the original calibration, independent of
// the trial temperatures. Native source refresh has already pinned every
// mechanical, kinetic, expansion, Cp and density claim in `selection`.
#[derive(Clone, Debug, PartialEq)]
struct Family {
    context: QueryPoint,
    temperature_axis: String,
    selection: MaterialPropertySelection,
    hardening_basis: HardeningPropertyBasis,
    kinetics: Option<OverstressSelection>,
    reference_temperature_k: f64,
    strain_limit: f64,
    density_kg_m3: f64,
    // The caloric owner hashes every original Cp/density receipt, in addition
    // to the complete chart. Equal numerical knots from other claims differ.
    caloric_curve_identity: ContentHash,
}

#[derive(Clone, Debug)]
struct Prepared {
    family: Family,
    hardening: Arc<ResolvedThermalCauchyTabulatedHardening>,
    old: ThermalTabulatedJ2,
    next: ThermalTabulatedJ2,
    mechanical: Mechanical,
    energy: PreparedThermoplasticEnergy,
    old_temperature_k: f64,
    temperature_k: f64,
    duration_s: f64,
}

/// One owned source-bound physical interval, reusable across fibers and trials.
/// Cloning shares all preparation, including original observations and tables.
#[derive(Clone, Debug)]
pub struct UniaxialThermalTabulated {
    prepared: Arc<Prepared>,
}

impl ThermalTabulatedSource<'_> {
    /// Prepare an actual thermal/mechanical interval without advancing history.
    /// Both temperatures must belong to the original caloric/source span. The
    /// same positive duration is retained even for rate-independent selection,
    /// because its structural owner is advancing one physical thermal interval.
    pub fn at_uniaxial_temperatures(
        &self,
        old_temperature_k: f64,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<UniaxialThermalTabulated, FeedbackError> {
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(FeedbackError::Invalid(
                "thermal axial response requires the actual positive physical duration",
            ));
        }
        let old = self.resolve_temperature(old_temperature_k, Some(duration_s))?;
        let next = self.resolve_temperature(temperature_k, Some(duration_s))?;
        let mechanical = match &next.rate {
            Some(rate) => Mechanical::Overstress(
                UniaxialTabulatedOverstress::new(rate.clone())
                    .map_err(|error| FeedbackError::Source(error.into()))?,
            ),
            None => Mechanical::Independent(
                UniaxialTabulated::new(next.law.clone())
                    .map_err(|error| FeedbackError::Source(error.into()))?,
            ),
        };
        let variable_hardening = self
            .hardening
            .hardening()
            .maximum_absolute_storage_slope_j_m3_k()
            != 0.0;
        let energy = self.base.prepare_energy_at_temperatures(
            old_temperature_k,
            temperature_k,
            &next.source,
            variable_hardening,
        )?;
        let family = Family {
            context: self.base.context.clone(),
            temperature_axis: self.base.temperature_axis.clone(),
            selection: self.base.selection.clone(),
            hardening_basis: self.base.source.hardening_basis(),
            kinetics: self.base.overstress,
            reference_temperature_k: self.base.reference_temperature_k,
            strain_limit: self.base.source.law().elastic.strain_limit,
            density_kg_m3: self.base.density_kg_m3,
            caloric_curve_identity: self.base.caloric.curve().identity(),
        };
        Ok(UniaxialThermalTabulated {
            prepared: Arc::new(Prepared {
                family,
                hardening: self.hardening.clone(),
                old,
                next,
                mechanical,
                energy,
                old_temperature_k,
                temperature_k,
                duration_s,
            }),
        })
    }
}

impl UniaxialThermalTabulated {
    pub fn old_temperature_k(&self) -> f64 {
        self.prepared.old_temperature_k
    }
    pub fn temperature_k(&self) -> f64 {
        self.prepared.temperature_k
    }
    pub fn duration_s(&self) -> f64 {
        self.prepared.duration_s
    }
    pub fn density_kg_m3(&self) -> f64 {
        self.prepared.family.density_kg_m3
    }
    pub fn caloric_curve_identity(&self) -> ContentHash {
        self.prepared.family.caloric_curve_identity
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.prepared.next.source.material()
    }
    /// Native tensor law at the new temperature; the full source is untouched.
    pub fn tensor_law(&self) -> &ThermalTabulatedJ2 {
        &self.prepared.next
    }
    pub fn hardening_source(&self) -> &ResolvedThermalCauchyTabulatedHardening {
        &self.prepared.hardening
    }

    /// Admit a succeeding interval once, outside constitutive iterations.
    /// Temperatures, duration and their legitimately refreshed coefficients may
    /// change; original source choices, context and caloric reference cannot.
    pub fn admit_successor(&self, next: &Self) -> Result<(), FeedbackError> {
        let a = &self.prepared;
        let b = &next.prepared;
        if a.temperature_k != b.old_temperature_k
            || a.family != b.family
            || !(Arc::ptr_eq(&a.hardening, &b.hardening) || a.hardening == b.hardening)
        {
            return Err(FeedbackError::Invalid(
                "thermal axial successor must preserve original calibration and continue the accepted temperature",
            ));
        }
        Ok(())
    }

    fn endpoint(
        &self,
        state: &UniaxialThermalTabulatedState,
    ) -> Result<&ThermalTabulatedJ2, FeedbackError> {
        if state.temperature_k == self.old_temperature_k() {
            Ok(&self.prepared.old)
        } else if state.temperature_k == self.temperature_k() {
            Ok(&self.prepared.next)
        } else {
            Err(FeedbackError::Invalid(
                "retained axial temperature is not an endpoint of this physical interval",
            ))
        }
    }

    fn admit_total(&self, total: &Voigt, mechanical: &Voigt) -> Result<(), FeedbackError> {
        let limit = self.prepared.family.strain_limit;
        for strain in [total, mechanical] {
            let norm = contract(strain, strain);
            if strain.iter().any(|value| !value.is_finite())
                || !norm.is_finite()
                || norm > limit * limit
            {
                return Err(FeedbackError::Invalid(
                    "free-sided total or mechanical tensor exceeds its original small-strain domain",
                ));
            }
        }
        Ok(())
    }

    /// Reconstruct the actual free-sided total strain at an accepted endpoint.
    /// This inspection does not perform a return map or advance material time.
    pub fn tensor_strain(
        &self,
        axial_strain: f64,
        state: &UniaxialThermalTabulatedState,
    ) -> Result<Voigt, FeedbackError> {
        let endpoint = self.endpoint(state)?;
        let free = endpoint.expansion_strain;
        let axial = axial_strain - free;
        let plastic = state.mechanical.plastic_strain();
        let poisson = endpoint.law.backbone().elastic.poisson;
        let transverse = -poisson * (axial - plastic) - 0.5 * plastic;
        let mechanical = [axial, transverse, transverse, 0.0, 0.0, 0.0];
        let total = add_free(mechanical, free);
        self.admit_total(&total, &mechanical)?;
        endpoint
            .law
            .backbone()
            .validate_overstress_history(&mechanical, &state.mechanical.tensor_history())
            .map_err(FeedbackError::Mechanical)?;
        endpoint
            .law
            .hardening()
            .yield_stress_pa(state.mechanical.equivalent_plastic_strain())
            .map_err(source_error)?;
        Ok(total)
    }

    /// Check accepted history at its own temperature, preserving valid viscous
    /// overstress when kinetics were explicitly selected for this source.
    pub fn validate_history(
        &self,
        axial_strain: f64,
        state: &UniaxialThermalTabulatedState,
        yield_tolerance_pa: f64,
    ) -> Result<(), FeedbackError> {
        if !yield_tolerance_pa.is_finite() || yield_tolerance_pa <= 0.0 {
            return Err(FeedbackError::Invalid(
                "positive finite thermal axial yield allowance required",
            ));
        }
        let total = self.tensor_strain(axial_strain, state)?;
        let endpoint = self.endpoint(state)?;
        let mechanical = endpoint.mechanical_strain(&total);
        let history = state.mechanical.tensor_history();
        match &endpoint.rate {
            Some(rate) => rate.validate_history(&mechanical, &history),
            None => {
                endpoint
                    .law
                    .validate_committed_history(&mechanical, &history, yield_tolerance_pa)
            }
        }
        .map_err(FeedbackError::Mechanical)
    }

    /// Exact Helmholtz storage at the retained state's actual temperature.
    /// The separately reported entropy contribution is not added a second time.
    pub fn stored_energy(
        &self,
        axial_strain: f64,
        state: &UniaxialThermalTabulatedState,
    ) -> Result<StoredEnergy, FeedbackError> {
        let total = self.tensor_strain(axial_strain, state)?;
        let endpoint = self.endpoint(state)?;
        endpoint
            .law
            .stored_energy(
                &endpoint.mechanical_strain(&total),
                &state.mechanical.tensor_history(),
            )
            .map_err(FeedbackError::Mechanical)
    }

    /// Pure free-sided candidate using the endpoint's original scalar return.
    /// The old temperature must match this interval; reuse of a newly accepted
    /// hot state as another old state cannot silently create a second step.
    pub fn response(
        &self,
        axial_strain: f64,
        old: &UniaxialThermalTabulatedState,
    ) -> Result<UniaxialThermalTabulatedResponse, FeedbackError> {
        if old.temperature_k != self.old_temperature_k() {
            return Err(FeedbackError::Invalid(
                "thermal axial trial must retain the physical interval's old temperature",
            ));
        }
        let endpoint = &self.prepared.next;
        let candidate = self
            .prepared
            .mechanical
            .response(axial_strain - endpoint.expansion_strain, &old.mechanical)?;
        let total = add_free(candidate.strain, endpoint.expansion_strain);
        self.admit_total(&total, &candidate.strain)?;
        // Admit the full native predictor as well as the final history. The
        // scalar reduction must not evade a tensor elastic/history envelope.
        endpoint
            .law
            .backbone()
            .validate_overstress_history(&candidate.strain, &old.mechanical.tensor_history())
            .map_err(FeedbackError::Mechanical)?;
        Ok(UniaxialThermalTabulatedResponse {
            strain: total,
            stress_pa: candidate.stress_pa,
            tangent_pa: candidate.tangent_pa,
            state: UniaxialThermalTabulatedState {
                mechanical: candidate.state,
                temperature_k: self.temperature_k(),
            },
        })
    }

    /// Audit through the same tensor return and prepared thermal first-law
    /// account as the source family. All histories remain immutable on refusal.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old: &UniaxialThermalTabulatedState,
        strain: f64,
        tolerances: EnergyTolerances,
    ) -> Result<UniaxialThermalTabulatedIncrement, FeedbackError> {
        self.validate_history(old_strain, old, tolerances.yield_pa)?;
        let old_total = self.tensor_strain(old_strain, old)?;
        let candidate = self.response(strain, old)?;
        let old_state = old.mechanical.tensor_history();
        let old_mechanical = self.prepared.old.mechanical_strain(&old_total);
        let endpoint = &self.prepared.next;
        let next_mechanical = endpoint.mechanical_strain(&candidate.strain);
        let mechanical = match &endpoint.rate {
            Some(rate) => {
                rate.advance_with_energy(&old_mechanical, &old_state, &next_mechanical, tolerances)
            }
            None => endpoint.law.advance_with_energy_from_admitted_history(
                &old_mechanical,
                &old_state,
                &next_mechanical,
                tolerances,
            ),
        }
        .map_err(FeedbackError::Mechanical)?;
        let expected = candidate.state.mechanical.tensor_history();
        let modulus = endpoint.law.backbone().elastic.youngs;
        let mut residual = (mechanical.stress_pa[0] - candidate.stress_pa).abs();
        for value in &mechanical.stress_pa[1..] {
            residual = residual.max(value.abs());
        }
        for (actual, expected) in mechanical
            .state
            .plastic_strain
            .iter()
            .zip(expected.plastic_strain)
        {
            residual = residual.max(modulus * (actual - expected).abs());
        }
        residual = residual.max(modulus * (mechanical.state.alpha - expected.alpha).abs());
        if !residual.is_finite() || residual > tolerances.yield_pa {
            return Err(FeedbackError::Mechanical(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.yield_pa,
            }));
        }
        let hardening = self.prepared.hardening.hardening();
        let mut history = ThermalHistoryContribution::default();
        if hardening.maximum_absolute_storage_slope_j_m3_k() != 0.0 {
            let old_derivative = hardening
                .hardening_energy_temperature_derivative_j_m3_k(old_state.alpha)
                .map_err(source_error)?;
            let new_derivative = hardening
                .hardening_energy_temperature_derivative_j_m3_k(mechanical.state.alpha)
                .map_err(source_error)?;
            history.variable_storage = true;
            history.old_entropy_j_m3_k = -old_derivative;
            history.new_entropy_j_m3_k = -new_derivative;
            history.coefficient_storage_change_j_m3 = product(
                self.temperature_k() - self.old_temperature_k(),
                old_derivative,
            )?;
        }
        // A free-sided member's lateral faces carry no traction, so its stress
        // is uniaxial. The solve's lateral residual, refused above beyond the
        // yield tolerance, does no work on the member and enters no account.
        let mut mechanical = mechanical;
        mechanical.stress_pa = [mechanical.stress_pa[0], 0.0, 0.0, 0.0, 0.0, 0.0];
        let thermal = self.prepared.energy.finish(
            &old_total,
            &old_state,
            &candidate.strain,
            mechanical,
            history,
            tolerances,
        )?;
        Ok(UniaxialThermalTabulatedIncrement {
            strain: candidate.strain,
            stress_pa: candidate.stress_pa,
            tangent_pa: candidate.tangent_pa,
            state: candidate.state,
            thermal,
        })
    }
}

fn add_free(mut mechanical: Voigt, free: f64) -> Voigt {
    for value in &mut mechanical[..3] {
        *value += free;
    }
    mechanical
}

impl Uniaxial for UniaxialThermalTabulated {
    type State = UniaxialThermalTabulatedState;

    fn initial_state(&self) -> Self::State {
        Self::State {
            mechanical: self.prepared.mechanical.initial_state(),
            temperature_k: self.old_temperature_k(),
        }
    }
    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |response| response.stress_pa)
    }
    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |response| response.tangent_pa)
    }
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state).map_or(
            Self::State {
                mechanical: state.mechanical.clone(),
                temperature_k: f64::NAN,
            },
            |response| response.state,
        )
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.prepared.next.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.prepared.next.card()
    }
}
