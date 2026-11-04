//! Original multiaxis phase conductivity integrated over P1 enthalpy fields.
//!
//! The caller explicitly declares one scalar source isotropic or binds three
//! principal sources to a fixed orthonormal frame. Original caloric objects,
//! source paths, fixed coordinates and selected claims survive.
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
use crate::{
    ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable, ScalarField,
    ThermalBoundary,
};

mod field;
pub(in crate::transient::enthalpy) use field::PhaseFieldTransport;

/// Original phase conductivity with an explicit isotropic or orthotropic frame.
/// The independently declared rho0 is fixed reference mass per mesh volume.
/// A directional scalar source does not establish physical isotropy by itself.
#[derive(Debug)]
pub struct MaterialPhaseTransport<'a> {
    source: &'a PreparedPhaseConductivity<'a>,
    sources: [&'a PreparedPhaseConductivity<'a>; 3],
    // The original source knots are converted once into the shared leaf's
    // coordinate/value view. Their coordinate remains h [J/kg], never T [K].
    // Isotropic transport uses only the first table; the other two are empty.
    knots: [Vec<(f64, f64)>; 3],
    principal_axes: Option<[[f64; 3]; 3]>,
    enthalpy_bounds_j_kg: [f64; 2],
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
        let bounds = source.enthalpy_bounds_j_kg();
        let knots = source_knots(cx, source, bounds)?;
        let mut hash = DomainHasher::new("org.frankensim.fs-conduction.phase-field-transport.v1");
        hash.update(source.identity().as_bytes());
        hash.update(&reference_density_kg_m3.to_bits().to_le_bytes());
        poll(cx, knots.len())?;
        Ok(Self {
            source,
            sources: [source; 3],
            knots: [knots, Vec::new(), Vec::new()],
            principal_axes: None,
            enthalpy_bounds_j_kg: bounds,
            reference_density_kg_m3,
            numerical_basis: ConductivityModel::isotropic_declared(1.0)?,
            identity: hash.finalize(),
        })
    }

    /// Bind three original principal k(T,f) sources to one caloric history.
    ///
    /// Each row of `principal_axes` is the physical-frame direction of the
    /// corresponding source. The frame is fixed and must be orthonormal. All
    /// sources must borrow the exact same caloric object and retain compatible
    /// complete typed state coordinates at every admitted caloric segment.
    /// Their source domains may differ; only their nonempty common enthalpy
    /// interval is admitted. Each principal conductivity must be positive on
    /// that interval, including its endpoints and original interior knots.
    ///
    /// No mixture, phase-dependent reorientation or density update is inferred.
    /// Unlike a scalar multiplier of a fixed tensor, the principal ratios may
    /// change independently during a latent plateau where dT/dh is zero.
    pub fn try_new_orthotropic(
        cx: &Cx<'_>,
        sources: [&'a PreparedPhaseConductivity<'a>; 3],
        principal_axes: [[f64; 3]; 3],
        reference_density_kg_m3: f64,
    ) -> Result<Self, EnthalpyError> {
        poll(cx, 0)?;
        if !reference_density_kg_m3.is_finite() || reference_density_kg_m3 <= 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "positive fixed reference density required for phase field transport",
            ));
        }
        // Reuse the original frame admission. The unit tensor is only the
        // problem's identity-bound numerical basis, never physical conductivity.
        let numerical_basis = ConductivityModel::orthotropic(
            principal_axes,
            [
                ConductivityTable::declared(1.0)?,
                ConductivityTable::declared(1.0)?,
                ConductivityTable::declared(1.0)?,
            ],
        )?;
        let mut bounds = sources[0].enthalpy_bounds_j_kg();
        for source in &sources[1..] {
            if !std::ptr::eq(sources[0].caloric_curve(), source.caloric_curve()) {
                return Err(EnthalpyError::InvalidInput(
                    "orthotropic phase sources must borrow the same original caloric object",
                ));
            }
            let other = source.enthalpy_bounds_j_kg();
            bounds[0] = bounds[0].max(other[0]);
            bounds[1] = bounds[1].min(other[1]);
        }
        if bounds[0] >= bounds[1] {
            return Err(EnthalpyError::InvalidInput(
                "orthotropic phase sources require a common enthalpy interval",
            ));
        }
        admit_source_contexts(cx, &sources)?;
        let knots = [
            source_knots(cx, sources[0], bounds)?,
            source_knots(cx, sources[1], bounds)?,
            source_knots(cx, sources[2], bounds)?,
        ];
        let mut hash =
            DomainHasher::new("org.frankensim.fs-conduction.orthotropic-phase-field-transport.v1");
        for source in sources {
            hash.update(source.identity().as_bytes());
        }
        for axis in principal_axes {
            for component in axis {
                hash.update(&component.to_bits().to_le_bytes());
            }
        }
        hash.update(&reference_density_kg_m3.to_bits().to_le_bytes());
        poll(cx, 0)?;
        Ok(Self {
            source: sources[0],
            sources,
            knots,
            principal_axes: Some(principal_axes),
            enthalpy_bounds_j_kg: bounds,
            reference_density_kg_m3,
            numerical_basis,
            identity: hash.finalize(),
        })
    }

    /// Identity of the original sources, interpretation/frame and fixed rho0.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact isotropic source, or principal source zero for orthotropic laws.
    /// All principal sources borrow this same original caloric chart; this
    /// single-source accessor does not assert that an orthotropic law is scalar.
    #[must_use]
    pub const fn conductivity_source(&self) -> &'a PreparedPhaseConductivity<'a> {
        self.source
    }

    /// One isotropic source or three orthotropic sources in principal-axis order.
    #[must_use]
    pub fn conductivity_sources(&self) -> &[&'a PreparedPhaseConductivity<'a>] {
        if self.principal_axes.is_some() {
            &self.sources
        } else {
            &self.sources[..1]
        }
    }

    /// Fixed principal direction vectors as rows, absent for explicit isotropy.
    #[must_use]
    pub const fn principal_axes(&self) -> Option<[[f64; 3]; 3]> {
        self.principal_axes
    }

    /// Common inclusive source support in the original caloric enthalpy datum.
    #[must_use]
    pub const fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.enthalpy_bounds_j_kg
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

    /// Element-average physical heat flux, -<K(h_h)> grad(T_h), in W/m².
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

fn source_knots(
    cx: &Cx<'_>,
    source: &PreparedPhaseConductivity<'_>,
    bounds: [f64; 2],
) -> Result<Vec<(f64, f64)>, EnthalpyError> {
    let positive = |k: f64| {
        if !k.is_finite() || k <= 0.0 {
            Err(EnthalpyError::InvalidInput(
                "spatial phase conductivity must be strictly positive",
            ))
        } else {
            Ok(())
        }
    };
    for h in bounds {
        let sample = source.state_at_specific_enthalpy(h).map_err(|_| {
            EnthalpyError::InvalidInput("phase conductivity endpoint must be source-supported")
        })?;
        positive(sample.conductivity_w_m_k())?;
    }
    let mut knots = Vec::with_capacity(source.knots().len());
    for (index, &[h, k]) in source.knots().iter().enumerate() {
        if index % ASSEMBLY_TILE == 0 {
            poll(cx, index)?;
        }
        if h >= bounds[0] && h <= bounds[1] {
            positive(k)?;
        }
        knots.push((h, k));
    }
    poll(cx, knots.len())?;
    Ok(knots)
}

fn admit_source_contexts(
    cx: &Cx<'_>,
    sources: &[&PreparedPhaseConductivity<'_>; 3],
) -> Result<(), EnthalpyError> {
    // The original source owner guarantees that ONLY its declared T/f axes
    // evolve. Fixed values and quantities must agree exactly; source ranges
    // may be clipped differently along the shared caloric chart.
    let a = &sources[0]
        .source_paths()
        .first()
        .ok_or(EnthalpyError::InvalidInput(
            "phase source must retain its original query context",
        ))?
        .start;
    for (index, source) in sources[1..].iter().enumerate() {
        poll(cx, index)?;
        let b = &source
            .source_paths()
            .first()
            .ok_or(EnthalpyError::InvalidInput(
                "phase source must retain its original query context",
            ))?
            .start;
        if source.temperature_axis() != sources[0].temperature_axis()
            || source.liquid_fraction_axis() != sources[0].liquid_fraction_axis()
            || a.axes().keys().ne(b.axes().keys())
            || a.axis_quantities() != b.axis_quantities()
        {
            return Err(EnthalpyError::InvalidInput(
                "orthotropic phase sources require the same typed state coordinates",
            ));
        }
        if a.axes().iter().any(|(axis, value)| {
            axis != sources[0].temperature_axis()
                && Some(axis.as_str()) != sources[0].liquid_fraction_axis()
                && b.axes().get(axis) != Some(value)
        }) {
            return Err(EnthalpyError::InvalidInput(
                "orthotropic phase sources require the same fixed external conditions",
            ));
        }
    }
    Ok(())
}
