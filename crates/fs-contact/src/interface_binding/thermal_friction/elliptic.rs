//! Original directional coefficient claims through the shared elliptical law.
//!
//! Each principal direction is an already admitted `SourceThermalCoulomb`.
//! Both must describe the same ordered interface, fixed query and explicit
//! temperature/pressure axes. Their selected static and kinetic claims remain
//! pinned through every update. The effective interface texture direction is
//! explicit geometry; it is not inferred from opaque texture-frame names.

use fs_tribo::tangential::elliptic::{EllipticCoulomb, EllipticCoulombState};
use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};

use super::*;

/// The mechanical owner's admission of the interval's friction branch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EllipticFrictionMode {
    /// Admit static reaction first, as required for accepted rest.
    StaticOrSliding,
    /// Continue accepted motion on the kinetic ellipse. The mechanics owner
    /// must separately admit any zero-slip result as a physical arrest event.
    Kinetic,
}

/// Two sourced principal coefficient pairs with one retained material direction.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceThermalEllipticCoulomb {
    principal: [SourceThermalCoulomb; 2],
    initial_law: EllipticCoulomb,
    identity: ContentHash,
}

/// Accepted source receipts, actual slip and objectively transported texture.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceThermalEllipticState {
    owner: ContentHash,
    law: EllipticCoulomb,
    constitutive: EllipticCoulombState,
    sources: [ResolvedInterfaceStatePoint; 2],
}

impl SourceThermalEllipticState {
    /// Original native checkpoint in its accepted tangent frame.
    #[must_use]
    pub const fn constitutive_state(&self) -> &EllipticCoulombState {
        &self.constitutive
    }

    /// Actual selected claims, values, uncertainty and use receipts, by axis.
    #[must_use]
    pub const fn friction_sources(&self) -> &[ResolvedInterfaceStatePoint; 2] {
        &self.sources
    }

    /// Current source-derived coefficients; no mutable access to retained history.
    #[must_use]
    pub const fn constitutive_law(&self) -> EllipticCoulomb {
        self.law
    }
}

impl SourceThermalEllipticCoulomb {
    /// Compose two existing source owners, without querying a new material name
    /// or losing their selection policy, authority, conditions or source pins.
    ///
    /// The supplied nonzero tangent vector is the effective first material
    /// principal direction. Combining independently oriented surface textures
    /// requires an explicit external model; this constructor does not invent one.
    pub fn from_principal_sources(
        principal: [SourceThermalCoulomb; 2],
        principal_axis: [f64; 2],
    ) -> Result<Self, SourceThermalCoulombError> {
        if principal[0].interface != principal[1].interface
            || principal[0].initial_source.query() != principal[1].initial_source.query()
            || principal[0].axes != principal[1].axes
            || principal[0].query_axes != principal[1].query_axes
        {
            return Err(SourceThermalCoulombError::Invalid(
                "directional coefficients require the same ordered interface and source conditions",
            ));
        }
        let (static_mu, kinetic_mu) = principal_coefficients(
            &principal,
            [&principal[0].initial_source, &principal[1].initial_source],
        )?;
        let initial_law = EllipticCoulomb::new(static_mu, kinetic_mu, principal_axis)?;
        let mut hash = DomainHasher::new("org.frankensim.source-thermal-elliptic-coulomb.v1");
        for source in &principal {
            hash.update(source.identity.as_bytes());
        }
        for component in initial_law.initial_state().principal_axis() {
            hash.update(&component.to_bits().to_le_bytes());
        }
        Ok(Self {
            principal,
            initial_law,
            identity: hash.finalize(),
        })
    }

    /// Original ordered source owners. Axis order has physical meaning.
    #[must_use]
    pub const fn principal_sources(&self) -> &[SourceThermalCoulomb; 2] {
        &self.principal
    }

    /// Complete source-role, query, authority and initial-texture identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    fn check(&self, state: &SourceThermalEllipticState) -> Result<(), SourceThermalCoulombError> {
        if state.owner != self.identity {
            return Err(SourceThermalCoulombError::Invalid(
                "directional Coulomb checkpoint belongs to another original source",
            ));
        }
        Ok(())
    }

    /// Resolve both original coefficient pairs at the supplied temperature and
    /// pressure, then solve one associated ellipse against the actual body motion.
    /// Zero native storage permits coefficient changes without parameter work;
    /// accepted slip and transported texture survive the refresh unchanged.
    /// A failure on either source publishes neither source nor mechanical state.
    pub fn advance_sourced(
        &self,
        state: &SourceThermalEllipticState,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<SourceThermalEllipticState>, SourceThermalCoulombError> {
        self.check(state)?;
        let sources = [
            self.principal[0].source_at_input(&state.sources[0], input)?,
            self.principal[1].source_at_input(&state.sources[1], input)?,
        ];
        self.advance_with_sources(state, input, sources, EllipticFrictionMode::StaticOrSliding)
    }

    /// Advance at the exact query produced by the thermal/mechanical owner.
    /// Both principal source bundles refresh atomically, including explicitly
    /// owned surface temperatures, speed and path. Every unowned coordinate
    /// stays fixed, and the native checkpoint retains its full static limits
    /// even when the mechanics owner selects an already-moving kinetic branch.
    pub fn advance_at_query(
        &self,
        state: &SourceThermalEllipticState,
        input: &TangentialStepInput,
        point: &QueryPoint,
        mode: EllipticFrictionMode,
    ) -> Result<TangentialStep<SourceThermalEllipticState>, SourceThermalCoulombError> {
        self.check(state)?;
        let sources = [
            self.principal[0].source_at_query(&state.sources[0], input, point)?,
            self.principal[1].source_at_query(&state.sources[1], input, point)?,
        ];
        self.advance_with_sources(state, input, sources, mode)
    }

    fn advance_with_sources(
        &self,
        state: &SourceThermalEllipticState,
        input: &TangentialStepInput,
        sources: [ResolvedInterfaceStatePoint; 2],
        mode: EllipticFrictionMode,
    ) -> Result<TangentialStep<SourceThermalEllipticState>, SourceThermalCoulombError> {
        let (static_mu, kinetic_mu) =
            principal_coefficients(&self.principal, [&sources[0], &sources[1]])?;
        let (law, history) =
            state
                .law
                .continue_coefficients(&state.constitutive, static_mu, kinetic_mu)?;
        let step = match mode {
            EllipticFrictionMode::StaticOrSliding => law.advance(&history, input)?,
            EllipticFrictionMode::Kinetic => law.advance_kinetic(&history, input)?,
        };
        Ok(TangentialStep {
            traction_n: step.traction_n,
            slip_increment_m: step.slip_increment_m,
            state: SourceThermalEllipticState {
                owner: self.identity,
                law,
                constitutive: step.state,
                sources,
            },
            regime: step.regime,
            stored_change_j: step.stored_change_j,
            dissipated_j: step.dissipated_j,
            closure_residual_j: step.closure_residual_j,
            compatibility_residual_m: step.compatibility_residual_m,
            iterations: step.iterations,
        })
    }
}

impl TangentialLaw for SourceThermalEllipticCoulomb {
    type State = SourceThermalEllipticState;

    fn model_id(&self) -> &str {
        "source-temperature-dependent-elliptic-coulomb"
    }

    fn initial_state(&self) -> Self::State {
        Self::State {
            owner: self.identity,
            law: self.initial_law,
            constitutive: self.initial_law.initial_state(),
            sources: self.principal.each_ref().map(|p| p.initial_source.clone()),
        }
    }

    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        if self.check(state).is_ok() {
            0.0
        } else {
            f64::NAN
        }
    }

    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        self.advance_sourced(state, input).map_err(as_tribo)
    }
}

impl ObjectiveTangentialLaw for SourceThermalEllipticCoulomb {
    fn transport_state(
        &self,
        state: &Self::State,
        rotation: TangentialHistoryRotation,
    ) -> Result<Self::State, TriboError> {
        self.check(state).map_err(as_tribo)?;
        let mut next = state.clone();
        next.constitutive = state.law.transport_state(&state.constitutive, rotation)?;
        Ok(next)
    }
}

fn principal_coefficients(
    principal: &[SourceThermalCoulomb; 2],
    sources: [&ResolvedInterfaceStatePoint; 2],
) -> Result<([f64; 2], [f64; 2]), SourceThermalCoulombError> {
    let mut static_mu = [0.0; 2];
    let mut kinetic_mu = [0.0; 2];
    for i in 0..2 {
        let law = coefficients(
            sources[i],
            &principal[i].static_property,
            &principal[i].kinetic_property,
        )?;
        static_mu[i] = law.static_coefficient()?;
        kinetic_mu[i] = law.zero_speed_kinetic_coefficient()?;
    }
    Ok((static_mu, kinetic_mu))
}

fn as_tribo(error: SourceThermalCoulombError) -> TriboError {
    match error {
        SourceThermalCoulombError::Tribo(error) => error,
        SourceThermalCoulombError::Invalid(field) => TriboError::InvalidInput { field },
        _ => TriboError::InvalidInput {
            field: "original directional Coulomb coefficient source",
        },
    }
}
