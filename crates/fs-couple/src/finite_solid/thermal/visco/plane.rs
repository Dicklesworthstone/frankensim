//! Source-bound thermal reference changes for native thin viscoelastic solids.
//!
//! Prepare immutable temperature endpoints once per element, then borrow the
//! field for an accepted body. The native owner retains independent quadrature
//! histories, original area and thickness; no source query runs in a Newton
//! or local material iteration. Expansion and kinetic shifts always start from
//! the original datum and declared calibration.
//!
//! A reference event and the following physical interval have separate energy
//! reports. Reference input includes both fixed-deformation work and the exact
//! volume-constraint actuation. Neither is relaxation heat. This quasistatic
//! adapter introduces no temperature-dependent mass or caloric closure.

use super::{ExpansionViscoelastic, ExpansionViscoelasticError};
use fs_exec::Cx;
use fs_matdb::MaterialCard;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::visco_hyper::ViscoHyperelastic;
use fs_solid::hyper2d::visco::{
    FiniteViscoPlaneEnergy, FiniteViscoPlaneError, FiniteViscoPlaneMode,
    FiniteViscoPlaneReferenceEnergy, FiniteViscoPlaneSolid, FiniteViscoStepControls,
    ViscoelasticLawField,
};
use fs_solid::{Mesh2, NewtonReport, NewtonSettings, Patch, SolidError};

/// Failed source admission or native mechanics. The borrowed accepted field
/// and solid remain unchanged, including on a late physical-interval refusal.
#[derive(Debug)]
pub enum ExpansionViscoPlaneError {
    Invalid(&'static str),
    Source {
        element: usize,
        source: ExpansionViscoelasticError,
    },
    Native(FiniteViscoPlaneError),
}

impl core::fmt::Display for ExpansionViscoPlaneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "thermal thin viscoelastic solid: {message}"),
            Self::Source { element, source } => {
                write!(f, "thermal thin element {element}: {source}")
            }
            Self::Native(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for ExpansionViscoPlaneError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Invalid(_) => None,
            Self::Source { source, .. } => Some(source),
            Self::Native(source) => Some(source),
        }
    }
}

impl From<FiniteViscoPlaneError> for ExpansionViscoPlaneError {
    fn from(error: FiniteViscoPlaneError) -> Self {
        Self::Native(error)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ExpansionViscoPlaneError> {
    cx.checkpoint()
        .map_err(|_| FiniteViscoPlaneError::Native(SolidError::Cancelled).into())
}

#[derive(Clone, Debug)]
enum SourceField {
    Uniform(ExpansionViscoelastic),
    PerElement {
        endpoints: Vec<ExpansionViscoelastic>,
        // A single preparation-time copy supplies the native contiguous view.
        // It cannot diverge from the immutable source endpoints.
        laws: Vec<ViscoHyperelastic>,
    },
}

/// Prepared uniform or original-element thermal endpoints. Each native law
/// remains tied to its source, fixed frame, original datum, manufactured
/// prestrain, composition order and caller-declared branch calibration.
#[derive(Clone, Debug)]
pub struct ExpansionViscoPlaneField {
    data: SourceField,
}

impl ExpansionViscoPlaneField {
    pub fn uniform(endpoint: ExpansionViscoelastic) -> Self {
        Self {
            data: SourceField::Uniform(endpoint),
        }
    }

    /// Prepare one endpoint for every original element. Exact mesh coverage is
    /// checked when the field is bound to a body or used for a reference event.
    pub fn per_element(
        endpoints: Vec<ExpansionViscoelastic>,
        cx: &Cx<'_>,
    ) -> Result<Self, ExpansionViscoPlaneError> {
        checkpoint(cx)?;
        if endpoints.is_empty() {
            return Err(ExpansionViscoPlaneError::Invalid(
                "a per-element thermal field cannot be empty",
            ));
        }
        let mut laws = Vec::with_capacity(endpoints.len());
        for endpoint in &endpoints {
            checkpoint(cx)?;
            laws.push(endpoint.law().clone());
        }
        Ok(Self {
            data: SourceField::PerElement { endpoints, laws },
        })
    }

    /// One endpoint for a uniform field, otherwise original-element order.
    pub fn endpoints(&self) -> &[ExpansionViscoelastic] {
        match &self.data {
            SourceField::Uniform(endpoint) => std::slice::from_ref(endpoint),
            SourceField::PerElement { endpoints, .. } => endpoints,
        }
    }

    /// Uniform fields have no prescribed element count.
    pub fn element_count(&self) -> Option<usize> {
        match &self.data {
            SourceField::Uniform(_) => None,
            SourceField::PerElement { endpoints, .. } => Some(endpoints.len()),
        }
    }

    pub fn endpoint_at(&self, element: usize) -> Option<&ExpansionViscoelastic> {
        match &self.data {
            SourceField::Uniform(endpoint) => Some(endpoint),
            SourceField::PerElement { endpoints, .. } => endpoints.get(element),
        }
    }

    /// Immutable native view, with no source resolution or per-point cloning.
    pub fn native_laws(&self) -> ViscoelasticLawField<'_> {
        match &self.data {
            SourceField::Uniform(endpoint) => ViscoelasticLawField::Uniform(endpoint.law()),
            SourceField::PerElement { laws, .. } => ViscoelasticLawField::PerElement(laws),
        }
    }

    /// Re-evaluate the retained original expansion paths and kinetic spectra.
    /// Supply one card/temperature for a uniform field, otherwise exactly one
    /// pair per original element. This preserves the field's coverage shape.
    pub fn at_temperatures(
        &self,
        cards: &[&MaterialCard],
        temperatures_k: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Self, ExpansionViscoPlaneError> {
        checkpoint(cx)?;
        let endpoints = self.endpoints();
        if cards.len() != endpoints.len() || temperatures_k.len() != endpoints.len() {
            return Err(ExpansionViscoPlaneError::Invalid(
                "one source card and temperature per stored endpoint required",
            ));
        }
        let mut prepared = Vec::with_capacity(endpoints.len());
        for (element, ((endpoint, card), temperature)) in
            endpoints.iter().zip(cards).zip(temperatures_k).enumerate()
        {
            checkpoint(cx)?;
            prepared.push(
                endpoint
                    .at_temperature(card, *temperature)
                    .map_err(|source| ExpansionViscoPlaneError::Source { element, source })?,
            );
        }
        match &self.data {
            SourceField::Uniform(_) => Ok(Self::uniform(prepared.remove(0))),
            SourceField::PerElement { .. } => Self::per_element(prepared, cx),
        }
    }

    fn check_transition(
        &self,
        previous: &Self,
        element_count: usize,
        cx: &Cx<'_>,
    ) -> Result<(), ExpansionViscoPlaneError> {
        if self
            .element_count()
            .is_some_and(|count| count != element_count)
        {
            return Err(ExpansionViscoPlaneError::Invalid(
                "one target thermal endpoint per original element required",
            ));
        }
        for element in 0..element_count {
            checkpoint(cx)?;
            let target = self
                .endpoint_at(element)
                .ok_or(ExpansionViscoPlaneError::Invalid(
                    "missing target thermal endpoint",
                ))?;
            let original =
                previous
                    .endpoint_at(element)
                    .ok_or(ExpansionViscoPlaneError::Invalid(
                        "missing accepted thermal endpoint",
                    ))?;
            target
                .check_transition(original)
                .map_err(|source| ExpansionViscoPlaneError::Source { element, source })?;
        }
        Ok(())
    }
}

/// Source identity and native quadrature histories accepted together. The
/// immutable prepared field outlives this body; a successor borrows its own
/// prepared target field while retaining the same original mesh and thickness.
#[derive(Debug)]
pub struct ExpansionViscoPlaneSolid<'a> {
    source: &'a ExpansionViscoPlaneField,
    native: FiniteViscoPlaneSolid<'a>,
}

#[derive(Debug)]
pub struct ExpansionViscoPlaneStep<'a> {
    pub solid: ExpansionViscoPlaneSolid<'a>,
    pub equilibrium: NewtonReport,
    pub energy: FiniteViscoPlaneEnergy,
    pub reactions_n: Vec<(usize, f64)>,
}

#[derive(Debug)]
pub struct ExpansionViscoPlaneReferenceChange<'a> {
    pub solid: ExpansionViscoPlaneSolid<'a>,
    pub energy: FiniteViscoPlaneReferenceEnergy,
}

/// Atomic source/reference event followed by a physical mechanical interval.
/// The energy reports retain their distinct physical meanings and allowances.
#[derive(Debug)]
pub struct ExpansionViscoPlaneThermalStep<'a> {
    pub solid: ExpansionViscoPlaneSolid<'a>,
    pub reference: FiniteViscoPlaneReferenceEnergy,
    pub mechanical: FiniteViscoPlaneEnergy,
    pub equilibrium: NewtonReport,
    pub reactions_n: Vec<(usize, f64)>,
}

impl ExpansionViscoPlaneThermalStep<'_> {
    /// Total prescribed reference input, including exact-volume actuation once.
    pub fn reference_input_j(&self) -> f64 {
        self.reference.reference_work_j + self.reference.constraint_reference_work_j
    }
}

impl<'a> ExpansionViscoPlaneSolid<'a> {
    pub fn try_new(
        mesh: &'a Mesh2,
        source: &'a ExpansionViscoPlaneField,
        mode: FiniteViscoPlaneMode,
        reference_thickness_m: f64,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, ExpansionViscoPlaneError> {
        let native = FiniteViscoPlaneSolid::try_new(
            mesh,
            source.native_laws(),
            mode,
            reference_thickness_m,
            local,
            cx,
        )?;
        Ok(Self { source, native })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn try_new_from(
        mesh: &'a Mesh2,
        source: &'a ExpansionViscoPlaneField,
        mode: FiniteViscoPlaneMode,
        reference_thickness_m: f64,
        initial_displacement_m: &[[f64; 2]],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, ExpansionViscoPlaneError> {
        let native = FiniteViscoPlaneSolid::try_new_from(
            mesh,
            source.native_laws(),
            mode,
            reference_thickness_m,
            initial_displacement_m,
            local,
            cx,
        )?;
        Ok(Self { source, native })
    }

    pub fn native(&self) -> &FiniteViscoPlaneSolid<'a> {
        &self.native
    }

    pub fn source_field(&self) -> &'a ExpansionViscoPlaneField {
        self.source
    }

    pub fn advance(
        &self,
        dt_s: f64,
        dirichlet: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        traction: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        settings: NewtonSettings,
        controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<ExpansionViscoPlaneStep<'a>, ExpansionViscoPlaneError> {
        let step = self
            .native
            .advance(dt_s, dirichlet, traction, settings, controls, cx)?;
        Ok(ExpansionViscoPlaneStep {
            solid: Self {
                source: self.source,
                native: step.solid,
            },
            equilibrium: step.equilibrium,
            energy: step.energy,
            reactions_n: step.reactions_n,
        })
    }

    /// Preserve original source identity and all material memory. Uniform and
    /// per-element fields may be interchanged only when each original element
    /// retains the same source, frame, datum, prestrain order and calibration.
    pub fn change_reference<'b>(
        &self,
        target: &'b ExpansionViscoPlaneField,
        controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<ExpansionViscoPlaneReferenceChange<'b>, ExpansionViscoPlaneError>
    where
        'a: 'b,
    {
        checkpoint(cx)?;
        target.check_transition(self.source, self.native.reference_mesh().elems.len(), cx)?;
        let changed = self
            .native
            .change_reference(target.native_laws(), controls, cx)?;
        Ok(ExpansionViscoPlaneReferenceChange {
            solid: ExpansionViscoPlaneSolid {
                source: target,
                native: changed.solid,
            },
            energy: changed.energy,
        })
    }

    /// Propose the source event and physical interval as one candidate. Every
    /// source, local, global and energy gate must pass before a successor is
    /// returned. Refusal never partially accepts a target field or history.
    /// The caller refines the prescribed temperature path and physical time;
    /// this endpoint split does not infer an unseen heating trajectory.
    #[allow(clippy::too_many_arguments)]
    pub fn change_reference_and_advance<'b>(
        &self,
        target: &'b ExpansionViscoPlaneField,
        dt_s: f64,
        dirichlet: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        traction: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        settings: NewtonSettings,
        reference_controls: FiniteViscoStepControls,
        step_controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<ExpansionViscoPlaneThermalStep<'b>, ExpansionViscoPlaneError>
    where
        'a: 'b,
    {
        let changed = self.change_reference(target, reference_controls, cx)?;
        let step = changed
            .solid
            .advance(dt_s, dirichlet, traction, settings, step_controls, cx)?;
        Ok(ExpansionViscoPlaneThermalStep {
            solid: step.solid,
            reference: changed.energy,
            mechanical: step.energy,
            equilibrium: step.equilibrium,
            reactions_n: step.reactions_n,
        })
    }
}
