//! Temperature-dependent Coulomb friction on one immutable ordered interface.
//!
//! Each trial resolves the explicitly selected static and kinetic coefficients
//! at its contact temperature and, optionally, supplied nominal pressure. The
//! original claims remain pinned. Explicit physical owners may additionally
//! supply surface temperatures, speed and accumulated path through a complete
//! query; every undeclared external coordinate stays fixed.
//! The native set-valued Coulomb law owns the actual traction/compliance solve.
//! Its zero elastic potential permits coefficient evolution without silently
//! changing stored energy. This is not a transfer of an elastic contact history.

/// Directional principal coefficients from the same original ordered interface.
pub mod elliptic;
/// Source-preserving joint Hertz sliding/spin friction at actual thermal conditions.
pub mod hertz;
mod query;
pub use query::{FrictionQueryAxis, FrictionQueryAxisKind};

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MatDbError, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarPropertyRequirement, is_fixed_external_coordinate, resolve_interface_state_point,
};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};
use fs_tribo::{
    DryTangentialState, FrictionLaw, InputAuthority, TangentialLaw, TangentialStep,
    TangentialStepInput, TriboError,
};

use super::BoundDryInterfaceSystem;

/// Explicit correspondence between mechanical inputs and source coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ThermalFrictionAxes {
    /// Absolute temperature [K]; typed aliases or the legacy exact `T` axis.
    pub temperature: String,
    /// Optional nominal pressure [Pa], with static pressure quantity semantics.
    /// This value is supplied by the caller, never inferred from force or area.
    pub nominal_pressure: Option<String>,
}

/// Source-preserving dissipative friction with coefficients evaluated per trial.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceThermalCoulomb {
    interface: BoundDryInterfaceSystem,
    initial_source: ResolvedInterfaceStatePoint,
    initial_law: FrictionLaw,
    static_property: String,
    kinetic_property: String,
    axes: ThermalFrictionAxes,
    query_axes: Vec<FrictionQueryAxis>,
    identity: ContentHash,
}

/// Accumulated slip and latest property receipts, bound to the original source.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceThermalCoulombState {
    owner: ContentHash,
    constitutive: DryTangentialState,
    source: ResolvedInterfaceStatePoint,
}

impl SourceThermalCoulombState {
    /// Native accepted slip and contact temperature, exposed read-only.
    #[must_use]
    pub const fn constitutive_state(&self) -> &DryTangentialState {
        &self.constitutive
    }

    /// Actual last coefficient queries, values, uncertainty and usage receipts.
    #[must_use]
    pub const fn friction_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }
}

/// Complete source or constitutive refusal, before publishing any candidate.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceThermalCoulombError {
    /// The original selected source does not support the requested state.
    Material(MaterialStatePointError),
    /// An external coordinate is malformed or has the wrong quantity convention.
    Query(MatDbError),
    /// The native Coulomb law refuses its coefficients or mechanical inputs.
    Tribo(TriboError),
    /// A request lies outside this explicit temperature/pressure friction model.
    Invalid(&'static str),
}

impl fmt::Display for SourceThermalCoulombError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "source-backed thermal Coulomb friction: {self:?}")
    }
}
impl std::error::Error for SourceThermalCoulombError {}
impl From<MaterialStatePointError> for SourceThermalCoulombError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<MatDbError> for SourceThermalCoulombError {
    fn from(error: MatDbError) -> Self {
        Self::Query(error)
    }
}
impl From<TriboError> for SourceThermalCoulombError {
    fn from(error: TriboError) -> Self {
        Self::Tribo(error)
    }
}

impl SourceThermalCoulomb {
    /// Admit two explicit dimensionless properties from the original dry card.
    ///
    /// Equal requirements explicitly select one shared coefficient. Otherwise
    /// the property names must differ. Each selected claim must be a static
    /// scalar at external conditions: rate, slip, strain, area and history bounds
    /// refuse, even on constant-valued claims. No material-name defaults apply.
    /// The initial selection policy and receipts bind the checkpoint identity;
    /// all refreshes pin those same claims instead of reselecting another source.
    pub fn resolve(
        interface: BoundDryInterfaceSystem,
        initial_query: QueryPoint,
        static_requirement: ScalarPropertyRequirement,
        kinetic_requirement: ScalarPropertyRequirement,
        selection: MaterialPropertySelection,
        axes: ThermalFrictionAxes,
    ) -> Result<Self, SourceThermalCoulombError> {
        Self::resolve_with_query_axes(
            interface,
            initial_query,
            static_requirement,
            kinetic_requirement,
            selection,
            axes,
            Vec::new(),
        )
    }

    /// Admit additional coordinates supplied by an explicit physical owner.
    ///
    /// Surface temperatures, speed and accumulated path may evolve together.
    /// Every response axis AND applicability bound must be owned here or be a
    /// fixed external condition. Runtime callers must use `advance_at_query`;
    /// the temperature-only `TangentialLaw` entry refuses these extra axes.
    /// This owns source selection, not the thermal or mechanical evolution.
    #[allow(clippy::too_many_arguments)]
    pub fn resolve_with_query_axes(
        interface: BoundDryInterfaceSystem,
        initial_query: QueryPoint,
        static_requirement: ScalarPropertyRequirement,
        kinetic_requirement: ScalarPropertyRequirement,
        selection: MaterialPropertySelection,
        axes: ThermalFrictionAxes,
        mut query_axes: Vec<FrictionQueryAxis>,
    ) -> Result<Self, SourceThermalCoulombError> {
        query::validate_additional_axes(&initial_query, &axes, &query_axes)?;
        validate_axes(&initial_query, &axes, &query_axes)?;
        query_axes.sort_by(|a, b| a.name.cmp(&b.name));
        if static_requirement.quantity() != QuantitySpec::dimensional(Dims::NONE)
            || kinetic_requirement.quantity() != QuantitySpec::dimensional(Dims::NONE)
        {
            return Err(SourceThermalCoulombError::Invalid(
                "static and kinetic coefficients require dimensionless scalar properties",
            ));
        }
        let static_property = static_requirement.name().to_owned();
        let kinetic_property = kinetic_requirement.name().to_owned();
        let requirements = if static_requirement == kinetic_requirement {
            vec![static_requirement]
        } else if static_property == kinetic_property {
            return Err(SourceThermalCoulombError::Invalid(
                "distinct static and kinetic requirements must have distinct property names",
            ));
        } else {
            vec![static_requirement, kinetic_requirement]
        };
        let initial_source = resolve_interface_state_point(
            interface.card(),
            &initial_query,
            &requirements,
            selection,
        )?;
        let owned: Vec<_> = query_axes
            .iter()
            .map(|axis| {
                (
                    axis.name.as_str(),
                    initial_query.axis_quantities()[&axis.name],
                )
            })
            .collect();
        if initial_source.properties().iter().any(|p| {
            if owned.is_empty() {
                !p.is_static_at_fixed_external_conditions()
            } else {
                !p.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
            }
        }) {
            return Err(SourceThermalCoulombError::Invalid(
                "Coulomb coefficient sources require owned or fixed external conditions",
            ));
        }
        let initial_law = coefficients(&initial_source, &static_property, &kinetic_property)?;
        let mut hash = DomainHasher::new("org.frankensim.source-thermal-coulomb.v1");
        hash.update(initial_source.identity().as_bytes());
        hash.update(&[match interface.interface().provenance().authority() {
            InputAuthority::CallerDeclared => 0,
            InputAuthority::SyntheticFixture => 1,
            InputAuthority::Estimated => 2,
        }]);
        hash.update(&[u8::from(axes.nominal_pressure.is_some())]);
        for name in [
            static_property.as_str(),
            kinetic_property.as_str(),
            axes.temperature.as_str(),
            axes.nominal_pressure.as_deref().unwrap_or(""),
        ] {
            hash.update(&(name.len() as u64).to_le_bytes());
            hash.update(name.as_bytes());
        }
        if !query_axes.is_empty() {
            hash.update(b"\0explicit-friction-query-axes.v1");
            for axis in &query_axes {
                hash.update(&[axis.kind.identity_tag()]);
                hash.update(&(axis.name.len() as u64).to_le_bytes());
                hash.update(axis.name.as_bytes());
            }
        }
        Ok(Self {
            interface,
            initial_source,
            initial_law,
            static_property,
            kinetic_property,
            axes,
            query_axes,
            identity: hash.finalize(),
        })
    }

    /// Complete immutable ordered interface, texture/history and authority ceiling.
    #[must_use]
    pub const fn interface(&self) -> &BoundDryInterfaceSystem {
        &self.interface
    }

    /// Original selected coefficient claims and fixed external query coordinates.
    #[must_use]
    pub const fn initial_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.initial_source
    }

    /// Explicit temperature and pressure coordinate ownership.
    #[must_use]
    pub const fn axes(&self) -> &ThermalFrictionAxes {
        &self.axes
    }

    /// Additional evolving coordinates in canonical name order. The primary
    /// temperature and nominal pressure remain in `axes()`.
    #[must_use]
    pub fn query_axes(&self) -> &[FrictionQueryAxis] {
        &self.query_axes
    }

    /// Original source, selection, coefficient roles, axes and authority identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Resolve the original selected coefficients without advancing friction.
    ///
    /// Coupled physical owners may need the law before solving their complete
    /// mechanical system. This uses the same complete-query admission and
    /// pinned claims as `advance_at_query`; it invents no slip increment or
    /// constitutive checkpoint. The caller must close every evolving query
    /// coordinate against that mechanical/thermal candidate before publication.
    pub fn coefficients_at_query(
        &self,
        point: &QueryPoint,
    ) -> Result<(ResolvedInterfaceStatePoint, FrictionLaw), SourceThermalCoulombError> {
        let temperature = point.axes().get(&self.axes.temperature).copied();
        let pressure = self.axes.nominal_pressure.as_ref()
            .and_then(|axis| point.axes().get(axis).copied());
        let source = self.source_at_query_conditions(
            &self.initial_source, temperature, pressure, point,
        )?;
        let law = coefficients(&source, &self.static_property, &self.kinetic_property)?;
        Ok((source, law))
    }

    // Shared source refresh for isotropic and directional Coulomb owners. This
    // performs no mechanical step and never reselects an alternative claim.
    fn source_at_input(
        &self,
        previous: &ResolvedInterfaceStatePoint,
        input: &TangentialStepInput,
    ) -> Result<ResolvedInterfaceStatePoint, SourceThermalCoulombError> {
        let temperature = input
            .contact_temperature_k
            .filter(|t| t.is_finite() && *t > 0.0)
            .ok_or(SourceThermalCoulombError::Invalid(
                "thermal Coulomb requires a positive absolute contact temperature",
            ))?;
        self.source_at_thermal_conditions(previous, temperature, input.nominal_pressure_pa)
    }

    // Shared T/p-only admission for original Coulomb and Hertz consumers.
    // Additional physical axes require the complete-query owner; they cannot
    // silently remain at their initial values during a Hertz thermal trial.
    fn source_at_thermal_conditions(
        &self,
        previous: &ResolvedInterfaceStatePoint,
        temperature: f64,
        nominal_pressure_pa: Option<f64>,
    ) -> Result<ResolvedInterfaceStatePoint, SourceThermalCoulombError> {
        if !self.query_axes.is_empty() {
            return Err(SourceThermalCoulombError::Invalid(
                "owned friction coordinates require the complete physical query",
            ));
        }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(SourceThermalCoulombError::Invalid(
                "thermal Coulomb requires a positive absolute contact temperature",
            ));
        }
        let mut point = self
            .initial_source
            .query()
            .with(&self.axes.temperature, temperature)?;
        match (&self.axes.nominal_pressure, nominal_pressure_pa) {
            (Some(axis), Some(pressure)) if pressure.is_finite() && pressure >= 0.0 => {
                point = point.with(axis, pressure)?;
            }
            (None, None) => {}
            _ => {
                return Err(SourceThermalCoulombError::Invalid(
                    "nominal pressure must match the explicitly owned source axis",
                ));
            }
        }
        if point == previous.query() {
            Ok(previous.clone())
        } else {
            Ok(self
                .initial_source
                .at_query_point(self.interface.card(), &point)?)
        }
    }

    fn source_at_query(
        &self,
        previous: &ResolvedInterfaceStatePoint,
        input: &TangentialStepInput,
        point: &QueryPoint,
    ) -> Result<ResolvedInterfaceStatePoint, SourceThermalCoulombError> {
        self.source_at_query_conditions(
            previous,
            input.contact_temperature_k,
            input.nominal_pressure_pa,
            point,
        )
    }

    fn source_at_query_conditions(
        &self,
        previous: &ResolvedInterfaceStatePoint,
        temperature_k: Option<f64>,
        pressure_pa: Option<f64>,
        point: &QueryPoint,
    ) -> Result<ResolvedInterfaceStatePoint, SourceThermalCoulombError> {
        query::validate_refresh(
            &self.initial_source.query(),
            point,
            &self.axes,
            &self.query_axes,
            temperature_k,
            pressure_pa,
        )?;
        let previous_point = previous.query();
        if self.query_axes.iter().any(|axis| {
            axis.kind == FrictionQueryAxisKind::AccumulatedDistance
                && point.axes()[&axis.name] < previous_point.axes()[&axis.name]
        }) {
            return Err(SourceThermalCoulombError::Invalid(
                "accumulated friction path cannot retreat from its accepted checkpoint",
            ));
        }
        if point == &previous_point {
            Ok(previous.clone())
        } else {
            Ok(self
                .initial_source
                .at_query_point(self.interface.card(), point)?)
        }
    }

    /// Advance at the physical owner's complete endpoint query, retaining every
    /// original selected claim. Undeclared coordinate changes and mismatched
    /// mechanical temperature/pressure inputs refuse before a candidate exists.
    pub fn advance_at_query(
        &self,
        state: &SourceThermalCoulombState,
        input: &TangentialStepInput,
        point: &QueryPoint,
    ) -> Result<TangentialStep<SourceThermalCoulombState>, SourceThermalCoulombError> {
        if state.owner != self.identity {
            return Err(SourceThermalCoulombError::Invalid(
                "thermal Coulomb checkpoint belongs to another original source",
            ));
        }
        let source = self.source_at_query(&state.source, input, point)?;
        self.advance_with_source(state, input, source)
    }

    /// Prepare one complete candidate with the original material error preserved.
    ///
    /// The supplied temperature and pressure are held for this constitutive
    /// interval. This operation only borrows accepted history; source-domain,
    /// coefficient-ordering and mechanical-solve failures leave it unchanged.
    /// The returned native dissipated work is allocated by the thermal consumer.
    pub fn advance_sourced(
        &self,
        state: &SourceThermalCoulombState,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<SourceThermalCoulombState>, SourceThermalCoulombError> {
        if state.owner != self.identity {
            return Err(SourceThermalCoulombError::Invalid(
                "thermal Coulomb checkpoint belongs to another original source",
            ));
        }
        let source = self.source_at_input(&state.source, input)?;
        self.advance_with_source(state, input, source)
    }

    fn advance_with_source(
        &self,
        state: &SourceThermalCoulombState,
        input: &TangentialStepInput,
        source: ResolvedInterfaceStatePoint,
    ) -> Result<TangentialStep<SourceThermalCoulombState>, SourceThermalCoulombError> {
        let law = coefficients(&source, &self.static_property, &self.kinetic_property)?;
        let step = law.advance(&state.constitutive, input)?;
        Ok(TangentialStep {
            traction_n: step.traction_n,
            slip_increment_m: step.slip_increment_m,
            state: SourceThermalCoulombState {
                owner: self.identity,
                constitutive: step.state,
                source,
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

impl TangentialLaw for SourceThermalCoulomb {
    type State = SourceThermalCoulombState;

    fn model_id(&self) -> &str {
        "source-temperature-dependent-coulomb"
    }

    fn initial_state(&self) -> Self::State {
        SourceThermalCoulombState {
            owner: self.identity,
            constitutive: self.initial_law.initial_state(),
            source: self.initial_source.clone(),
        }
    }

    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        if state.owner == self.identity {
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
        self.advance_sourced(state, input)
            .map_err(|error| match error {
                SourceThermalCoulombError::Tribo(error) => error,
                SourceThermalCoulombError::Invalid(field) => TriboError::InvalidInput { field },
                _ => TriboError::InvalidInput {
                    field: "original thermal Coulomb coefficient source",
                },
            })
    }
}

impl ObjectiveTangentialLaw for SourceThermalCoulomb {
    fn transport_state(
        &self,
        state: &SourceThermalCoulombState,
        rotation: TangentialHistoryRotation,
    ) -> Result<SourceThermalCoulombState, TriboError> {
        if state.owner != self.identity {
            return Err(TriboError::InvalidInput {
                field: "thermal Coulomb checkpoint belongs to another original source",
            });
        }
        // Orientation changes only the native slip vector. The original owner,
        // last temperature/pressure and selected coefficient receipts survive;
        // the following physical step performs any actual thermal refresh.
        let constitutive = self
            .initial_law
            .transport_state(&state.constitutive, rotation)?;
        Ok(SourceThermalCoulombState {
            owner: state.owner,
            constitutive,
            source: state.source.clone(),
        })
    }
}

fn coefficients(
    source: &ResolvedInterfaceStatePoint,
    static_property: &str,
    kinetic_property: &str,
) -> Result<FrictionLaw, SourceThermalCoulombError> {
    let value = |name| {
        source
            .property(name)
            .map(|p| p.value_si())
            .ok_or(SourceThermalCoulombError::Invalid(
                "selected Coulomb coefficient is ambiguous or missing",
            ))
    };
    let law = FrictionLaw::Coulomb {
        static_mu: value(static_property)?,
        kinetic_mu: value(kinetic_property)?,
    };
    // The original owner enforces finite nonnegative coefficients and mu_k <= mu_s.
    law.static_coefficient()?;
    Ok(law)
}

fn validate_axes(
    point: &QueryPoint,
    axes: &ThermalFrictionAxes,
    query_axes: &[FrictionQueryAxis],
) -> Result<(), SourceThermalCoulombError> {
    let temperature_quantity = point.axis_quantities().get(&axes.temperature).copied();
    let semantic_temperature = temperature_quantity
        .and_then(QuantitySpec::semantic_type)
        .is_some_and(|q| {
            q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
        });
    let legacy_temperature = axes.temperature == "T"
        && (temperature_quantity.is_none()
            || temperature_quantity == Some(QuantitySpec::dimensional(Temperature::DIMS)));
    if !(semantic_temperature || legacy_temperature)
        || !point
            .axes()
            .get(&axes.temperature)
            .is_some_and(|t| t.is_finite() && *t > 0.0)
        || point.axes().iter().any(|(axis, &value)| {
            !query_axes.iter().any(|owned| &owned.name == axis)
                && !is_fixed_external_coordinate(
                    axis,
                    point.axis_quantities().get(axis).copied(),
                    Some(value),
                )
        })
    {
        return Err(SourceThermalCoulombError::Invalid(
            "source coordinates must be static external conditions with absolute temperature",
        ));
    }
    if let Some(axis) = &axes.nominal_pressure {
        let pressure = point
            .axis_quantities()
            .get(axis)
            .copied()
            .and_then(QuantitySpec::semantic_type)
            .is_some_and(|q| q.kind() == QuantityKind::Pressure && q.form() == ValueForm::Static);
        if axis == &axes.temperature
            || !pressure
            || !point
                .axes()
                .get(axis)
                .is_some_and(|p| p.is_finite() && *p >= 0.0)
        {
            return Err(SourceThermalCoulombError::Invalid(
                "nominal pressure requires a distinct nonnegative static pressure coordinate",
            ));
        }
    }
    Ok(())
}
