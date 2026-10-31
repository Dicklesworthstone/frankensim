//! Finite solution activity coupled to the original bound-water diffusion owner.
//!
//! The aqueous solution and specimens exchange only water through an equilibrated
//! vapor path at equal prescribed pressure and temperature, using the same
//! pure-water reference pressure as the sorption source. Nonvolatile solutes
//! remain in the solution. The shared aqueous thermodynamics supplies `a_w`;
//! ideal-vapor equilibrium gives `RH = a_w` for the original source-backed
//! sorption boundary. The original material's diffusivity and optional vapor
//! surface-exchange resistance still control uptake and desorption.
//!
//! One bounded scalar solve closes endpoint solution solvent plus endpoint
//! bound water against their original inventories. Solute amounts never change,
//! so uptake concentrates the solution and desorption dilutes it. The accepted
//! source queries, native slab or Cartesian fields and solution publish together. Incompressible
//! solvent density is the explicit molar mass divided by partial molar volume.
//!
//! This is the equal-pressure, negligible-vapor-storage sorption limit. It does
//! not turn a vapor surface coefficient into a liquid-immersion permeability,
//! resolve an intervening vapor space, admit volatile solutes, transport ions
//! into the solid, or represent sorption heat and stress-assisted transport.
//! A separate selective-membrane flow law is needed for pressure-driven liquid
//! osmosis. Existing solid/plate consumers can consume the native moisture field.
//! Mechanical response is one way; elastic work does not feed solvent chemical
//! potential, and no closed chemo-mechanical energy balance is claimed.

mod shared;

// These are views of the original physical fields and their native receipts;
// they carry no vapor inventory or atmospheric solver state.
pub use super::atmosphere::{
    CartesianMoistureParticipant, MoistureParticipant, SharedMoistureSolid,
};
pub use shared::{
    SharedCartesianSolutionMoistureCandidate, SharedSolutionMoistureCandidate,
    SharedSolutionMoistureStep,
};

use super::{
    MoistureError, MoistureNonlinearConfig, MoistureSlab, MoistureStepConfig, MoistureStepReport,
    atmosphere::DryBasisSorption,
    cartesian::{MoistureCartesian, MoistureCartesianReport},
    poll, positive,
    reservoir::{ReservoirField, ReservoirParticipant},
};
use fs_exec::Cx;
use fs_material::osmosis::{AqueousSolution, OsmosisError};

/// Original thermodynamic/source/diffusion refusal; no endpoint is published.
#[derive(Debug)]
pub enum SolutionMoistureError {
    /// Original moisture/source/solid owner refused.
    Moisture(MoistureError),
    /// Original aqueous thermodynamic or concentration-support owner refused.
    Solution(OsmosisError),
    /// Malformed, unsupported or unrepresentable data.
    Invalid(&'static str),
    /// The supplied source-supported solvent-volume interval contains no sign change.
    Unbracketed { residual_kg: [f64; 2] },
    /// The bounded scalar solve did not close the requested water balance.
    NotConverged { residual_kg: f64, iterations: usize },
    /// Accepted endpoints do not close the original finite inventories.
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
}
impl core::fmt::Display for SolutionMoistureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite solution moisture exchange: {self:?}")
    }
}
impl std::error::Error for SolutionMoistureError {}
impl From<MoistureError> for SolutionMoistureError {
    fn from(value: MoistureError) -> Self {
        Self::Moisture(value)
    }
}
impl From<OsmosisError> for SolutionMoistureError {
    fn from(value: OsmosisError) -> Self {
        Self::Solution(value)
    }
}

/// Explicit scalar work, source-domain and finite-inventory tolerances.
#[derive(Clone, Copy, Debug)]
pub struct SolutionMoistureSolve {
    /// Positive ordered solvent volumes [m³] containing the accepted volume.
    /// All source/diffusion proposals in this interval must be admissible.
    pub solvent_volume_bracket_m3: [f64; 2],
    /// Maximum safeguarded scalar iterations, in 1..=256. Initial and two
    /// bracket-endpoint trials are additional; equilibrium and pure solvent
    /// each need only the initial trial across all participating native fields.
    pub maximum_iterations: usize,
    /// Absolute solution-plus-solid water conservation tolerance [kg].
    pub mass_tolerance_kg: f64,
}

/// One original finite solution inventory and an explicit solvent mass scale.
#[derive(Clone, Debug)]
pub struct FiniteMoistureSolution {
    solution: AqueousSolution,
    solvent_molar_mass_kg_mol: f64,
    solvent_density_kg_m3: f64,
}

/// Accepted native diffusion and finite-reservoir exchange over one interval.
#[derive(Clone, Debug)]
pub struct SolutionMoistureStep<Report = MoistureStepReport> {
    /// Original material/diffusion receipt, including endpoint source queries.
    pub transport: Report,
    /// Initial and endpoint solvent activity, dimensionless.
    pub water_activity: [f64; 2],
    /// Initial and endpoint osmotic pressure [Pa], as a thermodynamic diagnostic.
    /// It is not silently applied as a mechanical surface traction.
    pub osmotic_pressure_pa: [f64; 2],
    /// Initial and endpoint solvent inventory [kg].
    pub solvent_mass_kg: [f64; 2],
    /// Positive water transfer enters the original solid [kg].
    pub water_to_solid_kg: f64,
    /// Solution solvent change plus solid water change [kg].
    pub water_balance_residual_kg: f64,
    /// Number of interior scalar trials, excluding initial/endpoints.
    pub iterations: usize,
}

/// Solution and native moisture candidates admitted as one physical endpoint.
#[derive(Clone, Debug)]
pub struct SolutionMoistureCandidate<Field = MoistureSlab, Report = MoistureStepReport> {
    reservoir: FiniteMoistureSolution,
    moisture: Field,
    step: SolutionMoistureStep<Report>,
}
impl<Field, Report> SolutionMoistureCandidate<Field, Report> {
    /// Accepted solution reservoir candidate.
    pub const fn reservoir(&self) -> &FiniteMoistureSolution {
        &self.reservoir
    }
    /// Original native geometry with the accepted bound-water field.
    pub const fn moisture(&self) -> &Field {
        &self.moisture
    }
    /// Accepted source, transport and conservation report.
    pub const fn step(&self) -> &SolutionMoistureStep<Report> {
        &self.step
    }
    /// Publish both original inventories and their report together.
    pub fn into_parts(self) -> (FiniteMoistureSolution, Field, SolutionMoistureStep<Report>) {
        (self.reservoir, self.moisture, self.step)
    }
}

/// Six-face diffusion plus original finite-solution activity and water accounts.
pub type CartesianSolutionMoistureStep = SolutionMoistureStep<MoistureCartesianReport>;
/// Original solution and complete 3D moisture field admitted together.
pub type CartesianSolutionMoistureCandidate =
    SolutionMoistureCandidate<MoistureCartesian, MoistureCartesianReport>;

impl FiniteMoistureSolution {
    /// Retain the shared ideal or nonideal solution and an explicit water molar mass [kg/mol].
    /// Density follows the declared incompressible solvent volume convention;
    /// it is not inferred from a mixture name or a room-temperature constant.
    pub fn new(
        solution: impl Into<AqueousSolution>,
        solvent_molar_mass_kg_mol: f64,
    ) -> Result<Self, SolutionMoistureError> {
        let solution = solution.into();
        solution.water_activity()?;
        let density = solvent_molar_mass_kg_mol / solution.solvent_molar_volume_m3_mol();
        let mass = solution.solvent_volume_m3() * density;
        if !positive(solvent_molar_mass_kg_mol) || !positive(density) || !positive(mass) {
            return Err(SolutionMoistureError::Invalid(
                "finite positive solvent molar mass, density and water inventory required",
            ));
        }
        Ok(Self {
            solution,
            solvent_molar_mass_kg_mol,
            solvent_density_kg_m3: density,
        })
    }

    /// Shared thermodynamic owner and its original solute amounts.
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    /// Explicit solvent molar mass [kg/mol].
    pub const fn solvent_molar_mass_kg_mol(&self) -> f64 {
        self.solvent_molar_mass_kg_mol
    }
    /// Actual solvent inventory [kg], excluding dissolved solutes.
    pub fn solvent_mass_kg(&self) -> f64 {
        self.solution.solvent_volume_m3() * self.solvent_density_kg_m3
    }

    /// Advance native diffusion against this finite solution's endpoint activity.
    /// `None` explicitly seals a face; each `Some` retains the original vapor
    /// equilibrium sorption and optional vapor-exchange source. Both faces share
    /// this single inventory. Accepted inputs remain immutable on every refusal.
    ///
    /// The equal-pressure ideal-vapor bridge requires the solution's fixed
    /// temperature to equal the source/slab temperature. A safeguarded secant
    /// solve refines one finite solvent-volume bracket. Native interval, source,
    /// content-increment, linear-iteration and water-balance controls remain active.
    pub fn advance(
        &self,
        slab: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SolutionMoistureCandidate, SolutionMoistureError> {
        self.advance_field(slab, dt_s, faces, transport, controls, cx)
    }

    /// Couple all six original Cartesian faces to the finite solution's endpoint
    /// water activity. Solutes stay in the original thermodynamic owner while
    /// local uptake concentrates, and desorption dilutes, the one water inventory.
    /// This retains the same equal-pressure, negligible-vapor-storage limit;
    /// the sorption and surface coefficients remain vapor-path quantities.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_cartesian(
        &self,
        moisture: &MoistureCartesian,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<CartesianSolutionMoistureCandidate, SolutionMoistureError> {
        self.advance_field(moisture, dt_s, faces, (transport, nonlinear), controls, cx)
    }

    fn advance_field<F: ReservoirField<N>, const N: usize>(
        &self,
        field: &F,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; N],
        transport: F::Controls,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SolutionMoistureCandidate<F, F::Report>, SolutionMoistureError> {
        let participants = [ReservoirParticipant {
            id: 0,
            field,
            faces,
            transport,
        }];
        let candidate = self.advance_fields(
            &participants,
            dt_s,
            field.inventory().moisture.len(),
            controls,
            cx,
        )?;
        let (reservoir, mut solids, shared) = candidate.into_parts();
        let solid = solids.pop().expect("one admitted solution participant");
        let step = SolutionMoistureStep {
            transport: solid.transport,
            water_activity: shared.water_activity,
            osmotic_pressure_pa: shared.osmotic_pressure_pa,
            solvent_mass_kg: shared.solvent_mass_kg,
            water_to_solid_kg: shared.water_to_solids_kg,
            water_balance_residual_kg: shared.water_balance_residual_kg,
            iterations: shared.iterations,
        };
        Ok(SolutionMoistureCandidate {
            reservoir,
            moisture: solid.moisture,
            step,
        })
    }
}
