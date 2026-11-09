//! Bounded dispatch of native progressive axial ductility.
//!
//! These variants retain the material owner's return map, consistent tangent,
//! irreversible history and separate energy accounts. A prepared law owns one
//! physical duration and temperature; ordinary untimed frame routes reject it.
//! Cloning and duration rebinding share the immutable native model. There is no
//! virtual material call, source query or second damage calculation in dispatch.
//!
//! The axial model degrades both tensile and compressive axial response. It does
//! not inherit the three-dimensional model's intact bulk response or represent
//! crack closure/contact. A sourced initiation fit alone does not make authored
//! mechanical coefficients or equivalent-work capacity into measured data; use
//! `FiberMaterial::with_ductile_evolution` to bind the original mechanical source.

use crate::ductile::DuctileMaterialBinding;
use fs_material::ductile::evolution::{
    DuctileEvolutionError, DuctileJ2Conditions, DuctileJ2Settings, DuctileWorkEvolution,
};
use fs_material::ductile::source::{ResolvedDuctileLocus, ResolvedJohnsonCookInitiation};
use fs_material::ductile::uniaxial::{
    UniaxialDuctileJ2Law, UniaxialDuctileJ2State, UniaxialDuctileJ2Step,
};
use fs_material::ductile::{JohnsonCookInitiation, TabulatedDuctileLocus};
use fs_material::plastic::reduction::UniaxialJ2;
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Prepared native axial damage models, with one common owned history family.
#[derive(Clone, Debug)]
pub enum DuctileFiberLaw {
    /// Explicitly authored Johnson--Cook initiation and equivalent-work capacity.
    JohnsonCook(UniaxialDuctileJ2Law<JohnsonCookInitiation>),
    /// Explicitly authored tabulated initiation and equivalent-work capacity.
    Tabulated(UniaxialDuctileJ2Law<TabulatedDuctileLocus>),
    /// Original pinned Johnson--Cook fit; mechanical/capacity source admission
    /// belongs to the material-binding constructor, not this variant tag.
    SourcedJohnsonCook(UniaxialDuctileJ2Law<ResolvedJohnsonCookInitiation>),
    /// Original initiation observations, including a retained Lode surface.
    SourcedLocus(UniaxialDuctileJ2Law<ResolvedDuctileLocus>),
}

macro_rules! with_law {
    ($value:expr, $law:ident, $body:expr) => {
        match $value {
            DuctileFiberLaw::JohnsonCook($law) => $body,
            DuctileFiberLaw::Tabulated($law) => $body,
            DuctileFiberLaw::SourcedJohnsonCook($law) => $body,
            DuctileFiberLaw::SourcedLocus($law) => $body,
        }
    };
}

impl DuctileFiberLaw {
    /// Prepared physical duration [s] and actual fixed temperature [K].
    pub fn conditions(&self) -> DuctileJ2Conditions {
        with_law!(self, law, law.conditions())
    }

    /// Original characteristic length, equivalent-work capacity and local gates.
    pub fn settings(&self) -> DuctileJ2Settings {
        with_law!(self, law, law.model().settings())
    }

    /// Original normalized attenuation profile, never an invented residual floor.
    pub fn evolution(&self) -> DuctileWorkEvolution {
        with_law!(self, law, law.model().evolution())
    }

    /// Original free-sided axial J2 law, before degradation.
    pub fn native(&self) -> &UniaxialJ2 {
        with_law!(self, law, law.model().native())
    }

    /// Start an original model-owned virgin plastic/damage history.
    pub fn initial_state(&self) -> UniaxialDuctileJ2State {
        with_law!(self, law, law.initial_state())
    }

    /// Admit history ownership without comparing coefficients or rebuilding data.
    pub fn owns_state(&self, state: &UniaxialDuctileJ2State) -> bool {
        with_law!(self, law, law.model().owns_state(state))
    }

    /// Pure native response from the original accepted history.
    pub fn response(
        &self,
        strain: f64,
        state: &UniaxialDuctileJ2State,
    ) -> Result<UniaxialDuctileJ2Step, DuctileEvolutionError> {
        with_law!(self, law, law.response(strain, state))
    }

    /// Reprepare a physical interval while retaining the same immutable model.
    /// A thermal evolution owner must separately qualify changing temperature;
    /// ordinary frame damage paths use [`Self::with_duration`] at fixed temperature.
    pub fn with_conditions(
        &self,
        conditions: DuctileJ2Conditions,
    ) -> Result<Self, DuctileEvolutionError> {
        Ok(match self {
            Self::JohnsonCook(law) => Self::JohnsonCook(law.with_conditions(conditions)?),
            Self::Tabulated(law) => Self::Tabulated(law.with_conditions(conditions)?),
            Self::SourcedJohnsonCook(law) => {
                Self::SourcedJohnsonCook(law.with_conditions(conditions)?)
            }
            Self::SourcedLocus(law) => Self::SourcedLocus(law.with_conditions(conditions)?),
        })
    }

    /// Change only physical duration, preserving original actual temperature.
    pub fn with_duration(&self, duration_s: f64) -> Result<Self, DuctileEvolutionError> {
        self.with_conditions(DuctileJ2Conditions {
            duration_s,
            ..self.conditions()
        })
    }

    /// Whether two prepared intervals share the original model and histories.
    pub fn same_model(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::JohnsonCook(a), Self::JohnsonCook(b)) => a.model().same_model(b.model()),
            (Self::Tabulated(a), Self::Tabulated(b)) => a.model().same_model(b.model()),
            (Self::SourcedJohnsonCook(a), Self::SourcedJohnsonCook(b)) => {
                a.model().same_model(b.model())
            }
            (Self::SourcedLocus(a), Self::SourcedLocus(b)) => a.model().same_model(b.model()),
            _ => false,
        }
    }

    /// Authored numerical initiation carries no original material-card authority.
    pub fn permits_unbound_material(&self) -> bool {
        with_law!(
            self,
            law,
            law.model().initiation_law().permits_unbound_material()
        )
    }

    /// Match original initiation source and actual prepared temperature to the
    /// native mechanical material. Numeric laws cannot claim this source binding.
    pub fn matches_material(&self, material: &ResolvedMaterialStatePoint) -> bool {
        !self.permits_unbound_material()
            && with_law!(
                self,
                law,
                law.model()
                    .initiation_law()
                    .matches_material_temperature(material, law.conditions().temperature_k,)
            )
    }
}

impl From<UniaxialDuctileJ2Law<JohnsonCookInitiation>> for DuctileFiberLaw {
    fn from(law: UniaxialDuctileJ2Law<JohnsonCookInitiation>) -> Self {
        Self::JohnsonCook(law)
    }
}
impl From<UniaxialDuctileJ2Law<TabulatedDuctileLocus>> for DuctileFiberLaw {
    fn from(law: UniaxialDuctileJ2Law<TabulatedDuctileLocus>) -> Self {
        Self::Tabulated(law)
    }
}
impl From<UniaxialDuctileJ2Law<ResolvedJohnsonCookInitiation>> for DuctileFiberLaw {
    fn from(law: UniaxialDuctileJ2Law<ResolvedJohnsonCookInitiation>) -> Self {
        Self::SourcedJohnsonCook(law)
    }
}
impl From<UniaxialDuctileJ2Law<ResolvedDuctileLocus>> for DuctileFiberLaw {
    fn from(law: UniaxialDuctileJ2Law<ResolvedDuctileLocus>) -> Self {
        Self::SourcedLocus(law)
    }
}
