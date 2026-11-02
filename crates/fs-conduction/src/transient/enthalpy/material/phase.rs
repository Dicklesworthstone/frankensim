//! Original multiaxis phase conductivity integrated over P1 enthalpy fields.
//!
//! The caller explicitly declares the scalar source isotropic. Its original
//! caloric object, source paths, fixed coordinates and selected claim survive.
//! h_h is the P1 interpolant of nodal specific enthalpy. The existing scalar
//! moment leaf integrates k(h_h) and lambda_i dk/dh across every crossed source
//! knot; this is not interpolation of nodal k or a sample at mean phase.
//! The existing nodal T(h) interpolation, reference masses, sparse thermal
//! residual, Newton, radiation, adaptive stepping and adjoints remain owners.
//! This fixed-reference model introduces no phase mixture or latent account.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;

use super::super::{ASSEMBLY_TILE, EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, poll};
use crate::{ConductionMesh, ConductionProblem, ConductivityModel, ScalarField, ThermalBoundary};

mod field;
pub(in crate::transient::enthalpy) use field::PhaseFieldTransport;

/// One original phase source interpreted explicitly as isotropic conductivity.
/// The independently declared rho0 is fixed reference mass per mesh volume.
/// A directional scalar source does not establish physical isotropy by itself.
#[derive(Debug)]
pub struct MaterialPhaseTransport<'a> {
    source: &'a PreparedPhaseConductivity<'a>,
    // The original source knots are converted once into the shared leaf's
    // coordinate/value view. Their coordinate remains h [J/kg], never T [K].
    knots: Vec<(f64, f64)>,
    reference_density_kg_m3: f64,
    numerical_basis: ConductivityModel,
    identity: ContentHash,
}

impl<'a> MaterialPhaseTransport<'a> {
    /// Retain one complete original scalar k(T,f) source on its caloric h path.
    /// Native FEM conduction requires strictly positive k; insulating zero-k
    /// charts remain available to the separate finite-volume conduction owner.
    /// No conductivity or density is inferred from a material name.
    pub fn try_new_isotropic(
        cx: &Cx<'_>,
        source: &'a PreparedPhaseConductivity<'a>,
        reference_density_kg_m3: f64,
    ) -> Result<Self, EnthalpyError> {
        poll(cx, 0)?;
        if !reference_density_kg_m3.is_finite() || reference_density_kg_m3 <= 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "positive fixed reference density required for phase field transport",
            ));
        }
        let mut knots = Vec::with_capacity(source.knots().len());
        for (index, &[h, k]) in source.knots().iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            if !k.is_finite() || k <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "spatial phase conductivity must be strictly positive",
                ));
            }
            knots.push((h, k));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-conduction.phase-field-transport.v1");
        hash.update(source.identity().as_bytes());
        hash.update(&reference_density_kg_m3.to_bits().to_le_bytes());
        poll(cx, knots.len())?;
        Ok(Self {
            source,
            knots,
            reference_density_kg_m3,
            numerical_basis: ConductivityModel::isotropic_declared(1.0)?,
            identity: hash.finalize(),
        })
    }

    /// Identity of the original source, isotropic interpretation and fixed rho0.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact original source, including its borrowed caloric chart and paths.
    #[must_use]
    pub const fn conductivity_source(&self) -> &'a PreparedPhaseConductivity<'a> {
        self.source
    }

    /// Fixed reference mass per mesh volume, independent of current bulk density.
    #[must_use]
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density_kg_m3
    }

    /// Bind the SAME caloric chart to existing P1 reference-mass storage.
    /// The ordinary forward, endpoint, radiation, adaptive and adjoint APIs
    /// then evaluate the source-integrated conductivity through that owner.
    pub fn prepare_uniform<'m, 's>(
        &'s self,
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        budget: EnthalpyBudget,
    ) -> Result<EnthalpyBackwardEuler<'m, 's>, EnthalpyError> {
        let mut storage = EnthalpyBackwardEuler::uniform(
            cx,
            mesh,
            self.source.caloric_curve(),
            self.reference_density_kg_m3,
            budget,
        )?;
        storage.sourced_phase = Some(PhaseFieldTransport::Uniform(self));
        Ok(storage)
    }

    /// Retain this law's unit numerical basis. An unrelated material tensor,
    /// element override or second phase multiplier cannot multiply physical k.
    pub fn problem<'s>(
        &'s self,
        mesh: &'s ConductionMesh,
        boundary: &'s ThermalBoundary,
        source: &'s ScalarField,
    ) -> ConductionProblem<'s> {
        ConductionProblem {
            mesh,
            boundary,
            material: &self.numerical_basis,
            element_materials: None,
            source,
        }
    }

    /// Element-average physical heat flux, -<k(h_h)> grad(T_h), in W/m².
    /// T_h remains the native P1 interpolant of nodal chart temperatures.
    /// This uses the complete original conductivity, not the unit numerical
    /// basis retained by the problem. It is not a face-conservative recovery.
    pub fn element_heat_flux(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        h: &[f64],
    ) -> Result<Vec<[f64; 3]>, EnthalpyError> {
        PhaseFieldTransport::Uniform(self).element_heat_flux(cx, mesh, h)
    }
}
