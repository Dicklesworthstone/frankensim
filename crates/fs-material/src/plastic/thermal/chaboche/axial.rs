//! Source-bound thermal intervals for exact free-sided cyclic plastic fibers.
//!
//! Every Chaboche branch retains its native tensor history. Changing temperature
//! transports its strain-like coordinate X/C before the original implicit axial
//! return; no branch averaging or fitted effective backstress is introduced.
//! The prepared interval and native thermal work account are shared by clones.
//! Source queries and expansion integration stay outside constitutive trials.

use super::{
    ThermalChaboche, ThermalChabocheError, ThermalChabocheSource, branch_history_contribution,
};
use crate::fiber::Uniaxial;
use crate::plastic::chaboche::axial::{
    UniaxialChaboche, UniaxialChabocheOverstress, UniaxialChabocheResponse, UniaxialChabocheState,
};
use crate::plastic::chaboche::source::BackstressProperties;
use crate::plastic::chaboche::{ChabocheError, ChabocheJ2};
use crate::plastic::energy::{EnergyError, EnergyIncrement, EnergyTolerances, StoredEnergy};
use crate::plastic::source::ResolvedJ2Plasticity;
use crate::plastic::thermal::expansion::PreparedThermoplasticEnergy;
use crate::plastic::thermal::{
    FeedbackError, OverstressSelection, ThermoplasticIncrement, product,
};
use crate::state_point::ResolvedMaterialStatePoint;
use crate::tensor::contract;
use crate::{MaterialAdmissibility, SmallStrainLaw, Voigt};
use fs_blake3::ContentHash;
use fs_evidence::ModelCard;
use fs_matdb::QueryPoint;
use std::sync::Arc;

/// Accepted independent backstresses, permanent strain and actual temperature.
/// No interval preparation advances the retained history or constitutive time.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialThermalChabocheState<const N: usize = 1> {
    pub mechanical: UniaxialChabocheState<N>,
    temperature_k: f64,
}
impl<const N: usize> UniaxialThermalChabocheState<N> {
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
}

/// Pure axial candidate with full free-sided TOTAL strain and native history.
#[derive(Clone, Debug)]
pub struct UniaxialThermalChabocheResponse<const N: usize = 1> {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialThermalChabocheState<N>,
}

/// Native audited state and heat/work account for this physical interval.
#[derive(Clone, Debug)]
pub struct UniaxialThermalChabocheIncrement<const N: usize = 1> {
    pub strain: Voigt,
    pub stress_pa: f64,
    pub tangent_pa: f64,
    pub state: UniaxialThermalChabocheState<N>,
    pub thermal: ThermoplasticIncrement,
    /// Already included in the physical dissipation; never heat it twice.
    pub recovery_dissipation_j_m3: f64,
    /// Already included in the physical dissipation; not numerical loss.
    pub overstress_dissipation_j_m3: f64,
}

#[derive(Clone, Debug)]
enum Mechanical<const N: usize> {
    Independent(UniaxialChaboche<N>),
    Overstress(UniaxialChabocheOverstress<N>),
    Recovery(crate::plastic::chaboche::recovery::ChabocheStaticRecoveryJ2<N>),
}
impl<const N: usize> Mechanical<N> {
    fn response(
        &self,
        strain: f64,
        state: &UniaxialChabocheState<N>,
    ) -> Result<UniaxialChabocheResponse<N>, ChabocheError> {
        match self {
            Self::Independent(law) => law.response(strain, state),
            Self::Overstress(law) => law.response(strain, state),
            Self::Recovery(law) => law.uniaxial_response(strain, state),
        }
    }
}

/// Complete original reference bundle binds the immutable card, all selected
/// claims, quantity descriptors and coefficients; order binds each branch role.
#[derive(Clone, Debug, PartialEq)]
struct Family<const N: usize> {
    reference: ResolvedJ2Plasticity,
    context: QueryPoint,
    temperature_axis: String,
    kinetics: Option<OverstressSelection>,
    properties: [BackstressProperties; N],
    modulus_slopes_pa_k: [f64; N],
    static_recovery_properties: Option<[String; N]>,
    density_kg_m3: f64,
    caloric_curve_identity: ContentHash,
}

#[derive(Clone, Debug)]
struct Prepared<const N: usize> {
    family: Family<N>,
    next: ThermalChaboche<N>,
    mechanical: Mechanical<N>,
    energy: PreparedThermoplasticEnergy,
    old_temperature_k: f64,
    temperature_k: f64,
    old_expansion_strain: f64,
    duration_s: f64,
}

/// One immutable source-bound thermal interval for N independent backstresses.
/// Clones share every prepared law, coefficient account and original receipt.
#[derive(Clone, Debug)]
pub struct UniaxialThermalChaboche<const N: usize = 1> {
    prepared: Arc<Prepared<N>>,
}

impl<const N: usize> ThermalChabocheSource<'_, N> {
    /// Prepare the original selected family for one actual physical interval.
    /// The caloric owner and every source must cover both temperatures; source
    /// resolution, branch transport and Newton trials never advance the clock.
    pub fn at_uniaxial_temperatures(
        &self,
        old_temperature_k: f64,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<UniaxialThermalChaboche<N>, FeedbackError> {
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(FeedbackError::Invalid(
                "thermal cyclic fibers require the actual positive physical duration",
            ));
        }
        let next =
            self.resolve_at_temperatures(old_temperature_k, temperature_k, Some(duration_s))?;
        let mechanical = if let Some(recovery) = &next.recovery {
            Mechanical::Recovery(recovery.clone())
        } else {
            match &next.rate {
                Some(rate) => Mechanical::Overstress(
                    UniaxialChabocheOverstress::new(rate.clone())
                        .map_err(|error| FeedbackError::Source(error.into()))?,
                ),
                None => Mechanical::Independent(
                    UniaxialChaboche::new(next.law.clone())
                        .map_err(|error| FeedbackError::Source(error.into()))?,
                ),
            }
        };
        let old_expansion_strain = self
            .base
            .expansion_strain_between(self.base.reference_temperature_k, old_temperature_k)?;
        let energy = self.base.prepare_energy_at_temperatures(
            old_temperature_k,
            temperature_k,
            &next.source,
            self.modulus_slopes_pa_k.iter().any(|slope| *slope != 0.0),
        )?;
        let family = Family {
            reference: self.base.source.clone(),
            context: self.base.context.clone(),
            temperature_axis: self.base.temperature_axis.clone(),
            kinetics: self.base.overstress,
            properties: self.properties.clone(),
            modulus_slopes_pa_k: self.modulus_slopes_pa_k,
            static_recovery_properties: self.static_recovery_properties.clone(),
            density_kg_m3: self.base.density_kg_m3,
            // Includes original Cp and density receipts, not only equal knots.
            caloric_curve_identity: self.base.caloric.curve().identity(),
        };
        Ok(UniaxialThermalChaboche {
            prepared: Arc::new(Prepared {
                family,
                next,
                mechanical,
                energy,
                old_temperature_k,
                temperature_k,
                old_expansion_strain,
                duration_s,
            }),
        })
    }
}

impl<const N: usize> UniaxialThermalChaboche<N> {
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
    /// Original tensor transition, including its exact old/new branch transport.
    pub fn tensor_law(&self) -> &ThermalChaboche<N> {
        &self.prepared.next
    }

    /// Admit a successor once, outside material iterations. Temperature, dt and
    /// legitimately refreshed values may change; original calibration cannot.
    pub fn admit_successor(&self, next: &Self) -> Result<(), FeedbackError> {
        if self.temperature_k() != next.old_temperature_k()
            || self.prepared.family != next.prepared.family
        {
            return Err(FeedbackError::Invalid(
                "cyclic thermal successor must retain original ordered branches, source and caloric identity and continue the accepted temperature",
            ));
        }
        Ok(())
    }

    fn endpoint(
        &self,
        state: &UniaxialThermalChabocheState<N>,
    ) -> Result<(&ChabocheJ2<N>, f64), ThermalChabocheError> {
        if state.temperature_k == self.old_temperature_k() {
            Ok((
                &self.prepared.next.previous_law,
                self.prepared.old_expansion_strain,
            ))
        } else if state.temperature_k == self.temperature_k() {
            Ok((&self.prepared.next.law, self.prepared.next.expansion_strain))
        } else {
            Err(FeedbackError::Invalid(
                "retained cyclic temperature is not an endpoint of this physical interval",
            )
            .into())
        }
    }

    fn admit_total(&self, total: &Voigt, mechanical: &Voigt) -> Result<(), ThermalChabocheError> {
        let limit = self.prepared.family.reference.law().elastic.strain_limit;
        for strain in [total, mechanical] {
            let norm = contract(strain, strain);
            if strain.iter().any(|value| !value.is_finite())
                || !norm.is_finite()
                || norm > limit * limit
            {
                return Err(FeedbackError::Invalid(
                    "free-sided total or mechanical cyclic strain exceeds its original small-strain domain",
                )
                .into());
            }
        }
        Ok(())
    }

    /// Reconstruct accepted free-sided TOTAL strain without another return map.
    pub fn tensor_strain(
        &self,
        strain: f64,
        state: &UniaxialThermalChabocheState<N>,
    ) -> Result<Voigt, ThermalChabocheError> {
        let (law, free) = self.endpoint(state)?;
        let axial = strain - free;
        let plastic = state.mechanical.plastic_strain();
        let poisson = law.backbone().backbone().elastic.poisson;
        let transverse = -poisson * (axial - plastic) - 0.5 * plastic;
        let mechanical = [axial, transverse, transverse, 0.0, 0.0, 0.0];
        let total = add_free(mechanical, free);
        self.admit_total(&total, &mechanical)?;
        law.validate_history(&mechanical, state.mechanical.tensor_history())?;
        Ok(total)
    }

    /// Inspect the actual accepted temperature and native yield interval.
    /// A valid viscous state is not returned again or spuriously relaxed.
    pub fn validate_history(
        &self,
        strain: f64,
        state: &UniaxialThermalChabocheState<N>,
        yield_tolerance_pa: f64,
    ) -> Result<(), ThermalChabocheError> {
        if !yield_tolerance_pa.is_finite() || yield_tolerance_pa <= 0.0 {
            return Err(
                FeedbackError::Invalid("positive finite cyclic yield allowance required").into(),
            );
        }
        let total = self.tensor_strain(strain, state)?;
        let (law, free) = self.endpoint(state)?;
        let mechanical = subtract_free(total, free);
        let history = state.mechanical.tensor_history();
        let elastic = core::array::from_fn(|i| mechanical[i] - history.plastic.plastic_strain[i]);
        let stress = law.backbone().backbone().elastic.stress(&elastic, &());
        let residual = law.yield_function_with_state(&stress, history);
        if !residual.is_finite() {
            return Err(
                ChabocheError::Invalid("accepted cyclic yield is not representable").into(),
            );
        }
        if self.prepared.family.kinetics.is_none() && residual > yield_tolerance_pa {
            return Err(ChabocheError::Energy(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: yield_tolerance_pa,
            })
            .into());
        }
        Ok(())
    }

    /// Native Helmholtz storage at the retained state's actual temperature.
    /// Entropy and numerical dissipation are separate accounts.
    pub fn stored_energy(
        &self,
        strain: f64,
        state: &UniaxialThermalChabocheState<N>,
    ) -> Result<StoredEnergy, ThermalChabocheError> {
        let total = self.tensor_strain(strain, state)?;
        let (law, free) = self.endpoint(state)?;
        Ok(law.stored_energy(
            &subtract_free(total, free),
            state.mechanical.tensor_history(),
        )?)
    }

    /// Native scalar return with thermal transport of every independent branch.
    pub fn response(
        &self,
        strain: f64,
        old: &UniaxialThermalChabocheState<N>,
    ) -> Result<UniaxialThermalChabocheResponse<N>, ThermalChabocheError> {
        if old.temperature_k != self.old_temperature_k() {
            return Err(FeedbackError::Invalid(
                "cyclic thermal trial must retain the physical interval's old temperature",
            )
            .into());
        }
        let endpoint = &self.prepared.next;
        let transported = UniaxialChabocheState::from_transported_history(
            endpoint.transport_history(old.mechanical.tensor_history())?,
        );
        let candidate = self
            .prepared
            .mechanical
            .response(strain - endpoint.expansion_strain, &transported)?;
        let total = add_free(candidate.strain, endpoint.expansion_strain);
        self.admit_total(&total, &candidate.strain)?;
        endpoint
            .law
            .validate_history(&candidate.strain, transported.tensor_history())?;
        Ok(UniaxialThermalChabocheResponse {
            strain: total,
            stress_pa: candidate.stress_pa,
            tangent_pa: candidate.tangent_pa,
            state: UniaxialThermalChabocheState {
                mechanical: candidate.state,
                temperature_k: self.temperature_k(),
            },
        })
    }

    /// Audit through the original tensor account and prepared thermal owner.
    /// The audited native state is authoritative; scalar/tensor roundoff is
    /// checked in stress units for all plastic and backstress components.
    pub fn advance_with_energy(
        &self,
        old_strain: f64,
        old: &UniaxialThermalChabocheState<N>,
        strain: f64,
        tolerances: EnergyTolerances,
    ) -> Result<UniaxialThermalChabocheIncrement<N>, ThermalChabocheError> {
        self.validate_history(old_strain, old, tolerances.yield_pa)?;
        let old_total = self.tensor_strain(old_strain, old)?;
        let candidate = self.response(strain, old)?;
        let endpoint = &self.prepared.next;
        let old_state = old.mechanical.tensor_history();
        let transported = endpoint.transport_history(old_state)?;
        let old_mechanical = subtract_free(old_total, self.prepared.old_expansion_strain);
        let next_mechanical = subtract_free(candidate.strain, endpoint.expansion_strain);
        let mechanical = if let Some(recovery) = &endpoint.recovery {
            recovery.advance_with_energy(
                &old_mechanical,
                &transported,
                &next_mechanical,
                tolerances,
            )
        } else {
            match &endpoint.rate {
                Some(rate) => rate.advance_with_energy(
                    &old_mechanical,
                    &transported,
                    &next_mechanical,
                    tolerances,
                ),
                None => endpoint.law.advance_with_energy_from_admitted_history(
                    &old_mechanical,
                    &transported,
                    &next_mechanical,
                    tolerances,
                ),
            }
        }?;
        let expected = candidate.state.mechanical.tensor_history();
        let modulus = endpoint.law.backbone().backbone().elastic.youngs;
        let mut residual = (mechanical.stress_pa[0] - candidate.stress_pa).abs();
        for value in &mechanical.stress_pa[1..] {
            residual = residual.max(value.abs());
        }
        for (actual, expected) in mechanical
            .state
            .plastic
            .plastic_strain
            .iter()
            .zip(expected.plastic.plastic_strain)
        {
            residual = residual.max(modulus * (actual - expected).abs());
        }
        residual =
            residual.max(modulus * (mechanical.state.plastic.alpha - expected.plastic.alpha).abs());
        for (actual, expected) in mechanical
            .state
            .backstress_pa
            .iter()
            .flatten()
            .zip(expected.backstress_pa.iter().flatten())
        {
            residual = residual.max((actual - expected).abs());
        }
        if !residual.is_finite() || residual > tolerances.yield_pa {
            return Err(ChabocheError::Energy(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.yield_pa,
            })
            .into());
        }
        let contribution = branch_history_contribution(
            &self.prepared.family.modulus_slopes_pa_k,
            endpoint,
            old_state,
            &mechanical.state,
        )?;
        let delta = mechanical.state.plastic.alpha - old_state.plastic.alpha;
        let overstress_dissipation_j_m3 = match &endpoint.rate {
            Some(rate) => product(rate.endpoint_overstress_pa(delta)?, delta)?,
            None => 0.0,
        };
        let state = UniaxialThermalChabocheState {
            mechanical: UniaxialChabocheState::from_transported_history(mechanical.state.clone()),
            temperature_k: self.temperature_k(),
        };
        let recovery_dissipation_j_m3 = mechanical.recovery_dissipation_j_m3;
        let native = EnergyIncrement {
            strain: mechanical.strain,
            stress_pa: mechanical.stress_pa,
            state: mechanical.state.plastic,
            energy: mechanical.energy,
        };
        let thermal = self.prepared.energy.finish(
            &old_total,
            &old_state.plastic,
            &candidate.strain,
            native,
            contribution,
            tolerances,
        )?;
        Ok(UniaxialThermalChabocheIncrement {
            strain: candidate.strain,
            stress_pa: thermal.mechanical.stress_pa[0],
            tangent_pa: candidate.tangent_pa,
            state,
            thermal,
            recovery_dissipation_j_m3,
            overstress_dissipation_j_m3,
        })
    }
}

fn add_free(mut mechanical: Voigt, free: f64) -> Voigt {
    for value in &mut mechanical[..3] {
        *value += free;
    }
    mechanical
}
fn subtract_free(mut total: Voigt, free: f64) -> Voigt {
    for value in &mut total[..3] {
        *value -= free;
    }
    total
}

impl<const N: usize> Uniaxial for UniaxialThermalChaboche<N> {
    type State = UniaxialThermalChabocheState<N>;

    fn initial_state(&self) -> Self::State {
        Self::State {
            mechanical: UniaxialChabocheState::from_transported_history(
                self.prepared.next.previous_law.initial_state(),
            ),
            temperature_k: self.old_temperature_k(),
        }
    }
    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |r| r.stress_pa)
    }
    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |r| r.tangent_pa)
    }
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state).map_or(
            Self::State {
                mechanical: state.mechanical.clone(),
                temperature_k: f64::NAN,
            },
            |r| r.state,
        )
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            tangent_symmetric: true,
            ..self.prepared.next.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        self.prepared.next.card()
    }
}
