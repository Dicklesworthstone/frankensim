//! Source-preserving temperature endpoints for finite viscoelastic solids.
//!
//! The original engineering-length path supplies every natural deformation
//! `G(T)` relative to its ORIGINAL datum and fixed material axes. Optional
//! kinetic shifts always act on the original relaxation spectrum. Repeated
//! updates therefore neither compound expansion ratios nor rescale an already
//! shifted relaxation time.
//!
//! An explicit temperature/reference change keeps total deformation and viscous
//! metrics fixed and returns the signed change of mechanical storage. This is
//! prescribed-reference work, not dashpot heat or a caloric/entropy law. The
//! caller's thermal balance must account for that exchange. A subsequent native
//! mechanical interval uses the prepared law and produces its own physical
//! relaxation loss. Changing-temperature trajectories require time refinement.

use std::sync::Arc;

use super::ExpansionThermalReference;
use crate::finite_solid::FiniteSolidError;
use fs_matdb::MaterialCard;
use fs_material::MaterialError;
use fs_material::hyper::{Hyperelastic, reference::ReferenceDeformation};
use fs_material::visco::temperature::{HyperelasticTemperatureError, TemperatureShift};
use fs_material::visco_hyper::{
    ViscoHyperelastic, ViscoReferenceChange, ViscoRemainderReferenceChange, ViscoRemainderState,
    ViscoState,
};

/// Original source, material or temperature admission failed. Accepted owners
/// and histories are borrowed and remain unchanged on every refusal.
#[derive(Clone, Debug)]
pub enum ExpansionViscoelasticError {
    Thermal(FiniteSolidError),
    Kinetics(HyperelasticTemperatureError),
    Material(MaterialError),
    Incompatible(&'static str),
}

impl core::fmt::Display for ExpansionViscoelasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "finite viscoelastic expansion: {error}"),
            Self::Kinetics(error) => error.fmt(f),
            Self::Material(error) => error.fmt(f),
            Self::Incompatible(message) => write!(f, "finite viscoelastic expansion: {message}"),
        }
    }
}

impl std::error::Error for ExpansionViscoelasticError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Thermal(error) => Some(error),
            Self::Kinetics(error) => Some(error),
            Self::Material(error) => Some(error),
            Self::Incompatible(_) => None,
        }
    }
}

/// One held-temperature endpoint retaining the same expansion source and
/// native finite-strain material. Its borrowed [`Self::law`] drives the
/// existing full or mixed native solids without a second constitutive solver.
///
/// Expansion data are source-bound. The supplied equilibrium/branch energies
/// and optional kinetic calibration remain caller declarations; this owner
/// does not infer those measurements from the expansion card. Energies are
/// temperature-independent in the relaxed configuration. Mass and original
/// reference geometry remain with the native spatial owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpansionViscoelastic {
    original_law: Arc<ViscoHyperelastic>,
    thermal: ExpansionThermalReference,
    kinetics: Option<TemperatureShift>,
    law: ViscoHyperelastic,
}

impl ExpansionViscoelastic {
    /// Bind the original unshifted material and one source-admitted thermal
    /// endpoint. Relaxation times refer to the optional shift's declared
    /// calibration temperature, which need not equal the expansion datum.
    /// Without a shift declaration, kinetics are explicitly held constant.
    ///
    /// The supplied law must have the identity natural reference; an existing
    /// prestrain is refused rather than silently overwritten or composed in an
    /// inferred order. Preserve manufactured prestrain in its own explicit
    /// constitutive/reference owner when extending this thermal model.
    pub fn try_new(
        original_law: ViscoHyperelastic,
        thermal: ExpansionThermalReference,
        kinetics: Option<TemperatureShift>,
    ) -> Result<Self, ExpansionViscoelasticError> {
        if original_law.reference_deformation() != &ReferenceDeformation::identity() {
            return Err(ExpansionViscoelasticError::Incompatible(
                "the original thermal material must have an identity natural reference",
            ));
        }
        Self::prepare(Arc::new(original_law), thermal, kinetics)
    }

    fn prepare(
        original_law: Arc<ViscoHyperelastic>,
        thermal: ExpansionThermalReference,
        kinetics: Option<TemperatureShift>,
    ) -> Result<Self, ExpansionViscoelasticError> {
        let temperature_k = thermal.geometry().temperature_k();
        let law = if let Some(shift) = kinetics {
            shift
                .shift_hyperelastic(&original_law, temperature_k)
                .map_err(ExpansionViscoelasticError::Kinetics)?
                .into_law()
        } else {
            // Readmit public parameters through their original owners. No
            // return map or spectral solve is needed to prepare an endpoint.
            let equilibrium = Hyperelastic::new(
                original_law.equilibrium.model.clone(),
                original_law.equilibrium.stretch_limit,
            )
            .map_err(ExpansionViscoelasticError::Material)?;
            ViscoHyperelastic::new(equilibrium, original_law.branches.clone())
                .and_then(|law| {
                    law.with_local_newton(
                        original_law.local_newton_tolerance(),
                        original_law.local_newton_max_iterations(),
                    )
                })
                .map_err(ExpansionViscoelasticError::Material)?
        }
        .with_reference_deformation(*thermal.reference_deformation());
        Ok(Self {
            original_law,
            thermal,
            kinetics,
            law,
        })
    }

    /// Re-evaluate the same source path at an admitted absolute temperature.
    /// Original card identity, datum, fixed external conditions and material
    /// frame are retained. Source and kinetic validity are both enforced.
    /// The old endpoint and every retained history remain untouched.
    pub fn at_temperature(
        &self,
        card: &MaterialCard,
        temperature_k: f64,
    ) -> Result<Self, ExpansionViscoelasticError> {
        let thermal = ExpansionThermalReference::try_new(
            card,
            self.thermal.source().clone(),
            temperature_k,
            self.thermal.frame(),
        )
        .map_err(ExpansionViscoelasticError::Thermal)?;
        Self::prepare(Arc::clone(&self.original_law), thermal, self.kinetics)
    }

    /// Native equilibrium, viscous branches and current reference deformation.
    pub const fn law(&self) -> &ViscoHyperelastic {
        &self.law
    }

    /// Source path, original datum, fixed conditions and material axes.
    pub const fn thermal_reference(&self) -> &ExpansionThermalReference {
        &self.thermal
    }

    /// Independent original kinetic calibration, if explicitly supplied.
    pub const fn kinetics(&self) -> Option<TemperatureShift> {
        self.kinetics
    }

    /// Actual held temperature [K] of the expansion and kinetic endpoint.
    pub const fn temperature_k(&self) -> f64 {
        self.thermal.geometry().temperature_k()
    }

    fn check_transition(&self, previous: &Self) -> Result<(), ExpansionViscoelasticError> {
        if self.original_law != previous.original_law
            || self.kinetics != previous.kinetics
            || self.thermal.source().identity() != previous.thermal.source().identity()
            || self.thermal.source().material_card_identity()
                != previous.thermal.source().material_card_identity()
            || self.thermal.frame() != previous.thermal.frame()
            || self.thermal.reference_context() != previous.thermal.reference_context()
        {
            return Err(ExpansionViscoelasticError::Incompatible(
                "a thermal history must retain its original material, expansion path, datum, fixed frame and kinetic calibration",
            ));
        }
        Ok(())
    }

    /// Propose the full-history temperature/reference change at fixed total F
    /// and viscous metrics. Returned signed reference work must be retained in
    /// the external thermal/parameter energy account; it is not viscous heat.
    pub fn change_reference_from(
        &self,
        previous: &Self,
        accepted: &ViscoState,
    ) -> Result<ViscoReferenceChange, ExpansionViscoelasticError> {
        self.check_transition(previous)?;
        self.law
            .change_reference_from(&previous.law, accepted)
            .map_err(ExpansionViscoelasticError::Material)
    }

    /// Propose the same change for a bulk-free history. A native mixed solid
    /// must separately restore/admit its weak pressure constraint and account
    /// for bulk storage before accepting a body-level thermal endpoint.
    pub fn change_remainder_reference_from(
        &self,
        previous: &Self,
        accepted: &ViscoRemainderState,
    ) -> Result<ViscoRemainderReferenceChange, ExpansionViscoelasticError> {
        self.check_transition(previous)?;
        self.law
            .change_remainder_reference_from(&previous.law, accepted)
            .map_err(ExpansionViscoelasticError::Material)
    }
}
