//! Conservative molecular transport on a stationary constant-area gas mesh.
//!
//! Maxwell--Stefan species fluxes use the mass-averaged frame. Their original
//! partial enthalpies accompany diffusion; Fourier heat conduction and Newtonian
//! longitudinal stress share the same total-energy face flux. No independent
//! heat-release or viscous-heating source is added to conservative total energy.
//! Soret/Dufour cross effects and viscous stresses in variable-area geometry
//! require additional closures and are not represented by this operator.

pub(super) mod evolve;
mod law;

pub use law::{
    ConstantTransport, KineticTransport, TransportCoefficients, TransportDomain, TransportLaw,
};
#[cfg(feature = "material")]
pub(crate) use law::validate_binary;

use super::{
    Boundaries, Boundary, Error as FlowError, MAX_COMPONENTS, Mesh1d, MixtureEuler1d, MixtureGas,
    Primitive, Totals, allocated, positive,
};
use fs_qty::Time;
use fs_thermochem::transport::kinetic::KineticTransportErrorV1;
use fs_thermochem::transport::maxwell_stefan::{
    MaxwellStefanErrorV1, MaxwellStefanSolveV1, MaxwellStefanV1, MaxwellStefanWorkspaceV1,
};

/// Explicit diffusion work, step-size, and constitutive residual budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissipationControl {
    /// Maximum estimated parabolic number, in (0,0.45]. Checked at both stages.
    pub parabolic_cfl: f64,
    /// Largest internal molecular-transport substep [s].
    pub maximum_step_s: f64,
    /// Smallest permitted retry step [s], except an exact final remainder.
    pub minimum_step_s: f64,
    /// Maximum accepted plus rejected internal attempts, at most 1,000,000.
    pub maximum_steps: usize,
    /// Maximum halvings for one substep, at most 32.
    pub maximum_halvings: usize,
    /// Checks on the original Maxwell--Stefan equations and zero net mass flux.
    pub diffusion: MaxwellStefanSolveV1,
}
impl Default for DissipationControl {
    fn default() -> Self {
        Self {
            parabolic_cfl: 0.4,
            maximum_step_s: 1.0,
            minimum_step_s: 1e-16,
            maximum_steps: 10000,
            maximum_halvings: 20,
            diffusion: MaxwellStefanSolveV1::default(),
        }
    }
}
impl DissipationControl {
    pub(super) fn validate(self) -> Result<(), DissipationError> {
        if !positive(self.parabolic_cfl)
            || self.parabolic_cfl > 0.45
            || !positive(self.maximum_step_s)
            || !positive(self.minimum_step_s)
            || self.minimum_step_s > self.maximum_step_s
            || !(1..=1_000_000).contains(&self.maximum_steps)
            || self.maximum_halvings > 32
            || !positive(self.diffusion.relative_tolerance)
            || self.diffusion.relative_tolerance > 1e-3
            || !positive(self.diffusion.absolute_gradient_tolerance_per_m)
            || !positive(self.diffusion.absolute_mass_flux_tolerance_kg_m2_s)
        {
            return Err(DissipationError::Invalid("molecular transport controls"));
        }
        Ok(())
    }
}

/// A complete requested-duration molecular update; the flow clock is unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissipationReport {
    /// Complete physical duration represented by the molecular operator.
    pub elapsed: Time,
    /// Accepted internal SSP RK2 steps.
    pub accepted_steps: usize,
    /// Rejected and halved internal attempts.
    pub rejected_steps: usize,
    /// Largest estimated parabolic number over accepted stages.
    pub maximum_parabolic_number: f64,
    /// Largest original Maxwell--Stefan equation residual [1/m].
    pub maximum_gradient_residual_per_m: f64,
    /// Largest absolute net diffusive mass flux [kg/(m² s)].
    pub maximum_mass_frame_residual_kg_m2_s: f64,
    /// Initial physical integrals.
    pub before: Totals,
    /// Final physical integrals.
    pub after: Totals,
    /// Time/area-integrated left inflow minus right outflow.
    pub boundary_change: Totals,
    /// `after - before - boundary_change` for every conservative component.
    pub conservation_residual: Totals,
    /// Physical total-energy balance including original species references [J].
    pub absolute_energy_residual_j: f64,
}

/// Refusal before any caller field or clock is changed.
#[derive(Clone, Debug, PartialEq)]
pub enum DissipationError {
    /// Invalid source identity, geometry, state, law or controls.
    Invalid(&'static str),
    /// Original conservative/caloric flow owner refused.
    Flow(FlowError),
    /// Original Maxwell--Stefan constitutive owner refused.
    MaxwellStefan(MaxwellStefanErrorV1),
    /// Original kinetic transport-property owner refused.
    Kinetic(KineticTransportErrorV1),
    /// Original source-card material resolver refused the actual face state.
    #[cfg(feature = "material")]
    Material(fs_material::MaterialError),
    /// An accepted-stage parabolic estimate exceeds its declared bound.
    ParabolicNumber { observed: f64, maximum: f64 },
    /// Complete requested evolution exceeded its bounded attempt count.
    WorkLimit {
        elapsed_time_s: f64,
        accepted_steps: usize,
        rejected_steps: usize,
    },
    /// No admissible complete internal substep remained.
    StepFailure {
        elapsed_time_s: f64,
        step_s: f64,
        last: Box<DissipationError>,
    },
    /// Cancellation before complete publication.
    Cancelled,
}
impl core::fmt::Display for DissipationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "molecular gas transport: {self:?}")
    }
}
impl std::error::Error for DissipationError {}
impl From<FlowError> for DissipationError {
    fn from(error: FlowError) -> Self {
        match error {
            FlowError::Cancelled => Self::Cancelled,
            other => Self::Flow(other),
        }
    }
}
impl From<MaxwellStefanErrorV1> for DissipationError {
    fn from(error: MaxwellStefanErrorV1) -> Self {
        match error {
            MaxwellStefanErrorV1::Cancelled => Self::Cancelled,
            other => Self::MaxwellStefan(other),
        }
    }
}
impl From<KineticTransportErrorV1> for DissipationError {
    fn from(error: KineticTransportErrorV1) -> Self {
        match error {
            KineticTransportErrorV1::Cancelled => Self::Cancelled,
            other => Self::Kinetic(other),
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipationError> {
    if cancelled() {
        Err(DissipationError::Cancelled)
    } else {
        Ok(())
    }
}

struct Workspace {
    current: Vec<f64>,
    stage: Vec<f64>,
    candidate: Vec<f64>,
    primitive: Vec<Primitive>,
    heat_capacity_v: Vec<f64>,
    faces: Vec<f64>,
    face_rates: Vec<f64>,
    binary_diffusivities: Vec<f64>,
    diffusion: MaxwellStefanWorkspaceV1,
}

/// Prepared conservative molecular, heat and viscous field operator.
///
/// The original flow source, constant-area mesh, boundaries and transport law
/// are bound at construction. Complete requested-duration updates reuse one
/// flat candidate, SSP-stage, face and molecular workspace; physical state is
/// published once, after all internal substeps pass. The operator never moves
/// the flow's clock, allowing a composing integrator to own physical time.
pub struct DissipationOperator<L: TransportLaw> {
    gas: MixtureGas,
    mesh: Mesh1d,
    boundaries: Boundaries,
    law: L,
    maxwell_stefan: MaxwellStefanV1,
    widths: Vec<f64>,
    mass_ratio: f64,
    work: Workspace,
}
impl<L: TransportLaw> DissipationOperator<L> {
    /// Bind one stationary, constant-area, source-consistent flow.
    pub fn new(
        flow: &MixtureEuler1d,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if law.source() != flow.gas().source() {
            return Err(DissipationError::Invalid(
                "flow and molecular law have different original species banks",
            ));
        }
        if matches!(
            flow.boundaries,
            Boundaries::Ends {
                left: Boundary::Reservoir(_),
                ..
            } | Boundaries::Ends {
                right: Boundary::Reservoir(_),
                ..
            }
        ) {
            return Err(DissipationError::Invalid(
                "molecular reservoir requires species-bearing state and energy references",
            ));
        }
        let area = flow
            .mesh()
            .uniform_cross_section()
            .ok_or(DissipationError::Invalid(
                "molecular operator requires constant cross-sectional area",
            ))?;
        let count = flow.cell_count();
        let n = flow.gas().species_count();
        let m = n + 2;
        let mut widths = allocated(count, 0.0)?;
        for i in 0..count {
            poll(&mut cancelled)?;
            widths[i] = flow.mesh().uniform_cell_width().map_or_else(
                || flow.mesh().positions()[i + 1] - flow.mesh().positions()[i],
                |w| w.value(),
            );
            let expected = area * widths[i];
            if !positive(expected)
                || (flow.mesh().volumes()[i] - expected).abs()
                    > 64.0 * f64::EPSILON * expected.max(flow.mesh().volumes()[i])
            {
                return Err(DissipationError::Invalid(
                    "constant-area volume differs from area times cell width",
                ));
            }
        }
        let species = flow
            .gas()
            .source()
            .models()
            .iter()
            .map(|s| (s.species().clone(), s.molar_mass()))
            .collect();
        let maxwell_stefan = MaxwellStefanV1::new(species)?;
        let minimum = flow
            .gas()
            .source()
            .models()
            .iter()
            .map(|s| s.molar_mass().value())
            .fold(f64::INFINITY, f64::min);
        let maximum = flow
            .gas()
            .source()
            .models()
            .iter()
            .map(|s| s.molar_mass().value())
            .fold(0.0, f64::max);
        let mass_ratio = maximum / minimum;
        if !positive(mass_ratio) {
            return Err(DissipationError::Invalid(
                "molecular mass ratio is unrepresentable",
            ));
        }
        let work = Workspace {
            current: allocated(count * m, 0.0)?,
            stage: allocated(count * m, 0.0)?,
            candidate: allocated(count * m, 0.0)?,
            primitive: allocated(count, Primitive::default())?,
            heat_capacity_v: allocated(count, 0.0)?,
            faces: allocated((count + 1) * m, 0.0)?,
            face_rates: allocated(count + 1, 0.0)?,
            binary_diffusivities: allocated(n * n, 0.0)?,
            diffusion: maxwell_stefan.workspace()?,
        };
        poll(&mut cancelled)?;
        Ok(Self {
            gas: flow.gas().clone(),
            mesh: flow.mesh().clone(),
            boundaries: flow.boundaries,
            law,
            maxwell_stefan,
            widths,
            mass_ratio,
            work,
        })
    }
    /// Original immutable caloric/source bank admitted at construction.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Original transport-law owner, including its authored data and validity.
    pub const fn law(&self) -> &L {
        &self.law
    }

    /// Advance the entire requested duration without advancing the flow clock.
    ///
    /// SSP RK2 uses actual-state molecular coefficients and a parabolic scale
    /// estimate at both stages. Failed stages halve the step instead of clipping
    /// density, composition or temperature. This is an Estimated engineering
    /// discretization, not a fully discrete entropy or invariant-domain proof.
    pub fn advance_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DissipationReport, DissipationError> {
        evolve::advance(self, flow, duration, controls, &mut cancelled)
    }
}

#[cfg(test)]
mod tests;
