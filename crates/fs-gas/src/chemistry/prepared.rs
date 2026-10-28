use super::*;

/// Source-owning, prepared reactor for repeated independent fixed-U,V solves.
///
/// Construction admits the source interval once and reserves all concentration,
/// reaction-extent, Newton, stage and endpoint storage. Successful advances use
/// this storage without allocating. Inputs use the immutable mechanism's
/// canonical species order. The owner can be reused after every refusal;
/// unsuccessful scratch state is never returned as an accepted endpoint.
pub struct IsochoricReactor {
    pub(super) reactor: solve::Reactor,
    pub(super) current: Vec<f64>,
    pub(super) full: Vec<f64>,
    pub(super) half: Vec<f64>,
    pub(super) fine: Vec<f64>,
    pub(super) current_concentrations: Vec<f64>,
    pub(super) full_concentrations: Vec<f64>,
    pub(super) amounts: Vec<Amount>,
}

impl IsochoricReactor {
    /// Bind one complete original mechanism and admit its caloric interval.
    pub fn new(
        mechanism: IdealGasKineticsV1,
        controls: IsochoricKineticsSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IsochoricKineticsError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        mechanism.admit_temperature_interval(
            Temperature::new(controls.minimum_temperature_k),
            Temperature::new(controls.maximum_temperature_k),
            &mut cancelled,
        )?;
        let ns = mechanism.thermodynamics().models().len();
        let nr = mechanism.reactions().len();
        let mut amounts = Vec::new();
        amounts
            .try_reserve_exact(ns)
            .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
        amounts.resize(ns, Amount::new(0.0));
        let prepared = Self {
            reactor: solve::Reactor::new(mechanism, controls)?,
            current: zeros(nr)?,
            full: zeros(nr)?,
            half: zeros(nr)?,
            fine: zeros(nr)?,
            current_concentrations: zeros(ns)?,
            full_concentrations: zeros(ns)?,
            amounts,
        };
        poll(&mut cancelled)?;
        Ok(prepared)
    }

    /// Exact original source bank, reaction laws and canonical species order.
    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        &self.reactor.mechanism
    }

    /// Unchanged physical tolerances and work bounds admitted at construction.
    pub const fn controls(&self) -> IsochoricKineticsSolve {
        self.reactor.controls
    }

    /// Advance dense canonical amounts from their original source energy at T.
    ///
    /// The endpoint borrows prepared storage until dropped. A subsequent call
    /// starts from its supplied amounts, never from this previous endpoint.
    pub fn advance_adiabatic(
        &mut self,
        amounts: &[Amount],
        volume: Volume,
        initial_temperature: Temperature,
        duration: Time,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsochoricKineticsView<'_>, IsochoricKineticsError> {
        advance_prepared(
            self,
            amounts,
            None,
            volume,
            initial_temperature,
            duration,
            cancelled,
        )
    }

    /// Advance dense canonical amounts at one explicit total internal energy.
    ///
    /// Formation energy is included exactly once in `target`; the supplied
    /// temperature is a root hint. No caller buffer is mutated on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adiabatic_energy(
        &mut self,
        amounts: &[Amount],
        target: Energy,
        volume: Volume,
        initial_temperature: Temperature,
        duration: Time,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsochoricKineticsView<'_>, IsochoricKineticsError> {
        advance_prepared(
            self,
            amounts,
            Some(target),
            volume,
            initial_temperature,
            duration,
            cancelled,
        )
    }
}

/// Accepted endpoint borrowed from a prepared reactor, with no per-call allocation.
///
/// Original thermodynamic identity is retained by `mechanism()`, and
/// `temperature()` identifies the accepted source evaluation. Call `to_owned`
/// when independently owned arrays and individual NASA-9 receipts are needed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsochoricKineticsView<'a> {
    pub(super) mechanism: &'a IdealGasKineticsV1,
    pub(super) amounts: &'a [Amount],
    pub(super) temperature: Temperature,
    pub(super) pressure: Pressure,
    pub(super) internal_energy: Energy,
    pub(super) target_internal_energy: Energy,
    pub(super) mass: Mass,
    pub(super) volume: Volume,
    pub(super) elapsed_time: Time,
    pub(super) balance: IsochoricKineticsBalance,
}

impl IsochoricKineticsView<'_> {
    /// Exact immutable source bank and mechanism used for this endpoint.
    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.mechanism
    }
    /// Accepted amounts in canonical mechanism species order.
    pub const fn amounts(&self) -> &[Amount] {
        self.amounts
    }
    /// Actual accepted temperature.
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Pressure from actual molecular count and temperature.
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Actual endpoint total internal energy, J.
    pub const fn internal_energy(&self) -> Energy {
        self.internal_energy
    }
    /// Single original conserved energy target, J.
    pub const fn target_internal_energy(&self) -> Energy {
        self.target_internal_energy
    }
    /// Endpoint source molecular mass, kg.
    pub const fn mass(&self) -> Mass {
        self.mass
    }
    /// Fixed volume, m3.
    pub const fn volume(&self) -> Volume {
        self.volume
    }
    /// Complete requested duration, s.
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed_time
    }
    /// Actual conservation, entropy and bounded work report.
    pub const fn balance(&self) -> IsochoricKineticsBalance {
        self.balance
    }

    /// Materialize an independent endpoint and the original source receipts.
    /// This explicit conversion allocates; the prepared advance itself does not.
    pub fn to_owned(
        self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsochoricKineticsResult, IsochoricKineticsError> {
        poll(&mut cancelled)?;
        let mut amounts = Vec::new();
        amounts
            .try_reserve_exact(self.amounts.len())
            .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
        amounts.extend_from_slice(self.amounts);
        let mut sources = Vec::new();
        sources
            .try_reserve_exact(amounts.len())
            .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
        for model in self.mechanism.thermodynamics().models() {
            poll(&mut cancelled)?;
            let source = model.evaluate(self.temperature).map_err(|source| {
                KineticsErrorV1::SpeciesEvaluation {
                    species: model.species().clone(),
                    source,
                }
            })?;
            sources.push(source.receipt().clone());
        }
        let result = IsochoricKineticsResult {
            mechanism: self.mechanism.clone(),
            amounts,
            temperature: self.temperature,
            pressure: self.pressure,
            internal_energy: self.internal_energy,
            target_internal_energy: self.target_internal_energy,
            mass: self.mass,
            volume: self.volume,
            elapsed_time: self.elapsed_time,
            balance: self.balance,
            sources,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}
