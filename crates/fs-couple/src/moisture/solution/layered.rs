//! Original heterogeneous bodies exchanging solvent with one finite solution.
//!
//! Every solvent-volume trial advances the whole original stack or coated body.
//! Interior interfaces retain their own partition/resistance laws and never
//! become reservoir ports. Only exposed physical end faces see solution activity
//! through the existing equal-pressure, negligible-vapor-storage bridge. All
//! material-specific dry inventories debit the shared solvent exactly once.

use super::field::SolutionParticipant;
use super::*;
use crate::moisture::{
    layered::{
        MoistureActivityTransfer, MoistureLayeredFace, MoistureLayeredReport, MoistureLayeredSlab,
    },
    radial::layered::MoistureLayeredRadial,
    reservoir::water_sum,
    surface::admit_environment,
};

// These types contain native specimen/source references, not a gas inventory.
pub use crate::moisture::layered::atmosphere::LayeredMoistureParticipant;
pub use crate::moisture::radial::layered::atmosphere::LayeredRadialMoistureParticipant;

mod sdirk;
pub use sdirk::{
    LayeredRadialSolutionSdirk2Candidate, LayeredSolutionSdirk2Candidate,
    LayeredSolutionSdirk2Step, SharedLayeredRadialSolutionSdirk2Candidate,
    SharedLayeredSolutionSdirk2Candidate,
};

/// One original layered slab and the finite solution endpoint.
pub type LayeredSolutionMoistureCandidate =
    SolutionMoistureCandidate<MoistureLayeredSlab, MoistureLayeredReport>;
/// Independently sourced heterogeneous slabs sharing actual finite solvent.
pub type SharedLayeredSolutionMoistureCandidate =
    SharedSolutionMoistureCandidate<MoistureLayeredSlab, MoistureLayeredReport>;
/// One original coated cylinder/sphere and the finite solution endpoint.
pub type LayeredRadialSolutionMoistureCandidate =
    SolutionMoistureCandidate<MoistureLayeredRadial, MoistureLayeredReport>;
/// Independently sourced coated radial bodies sharing actual finite solvent.
pub type SharedLayeredRadialSolutionMoistureCandidate =
    SharedSolutionMoistureCandidate<MoistureLayeredRadial, MoistureLayeredReport>;

fn exposure(requested: [bool; 2], transfer: &[MoistureActivityTransfer; 2]) -> [bool; 2] {
    std::array::from_fn(|side| {
        requested[side] && transfer[side] != MoistureActivityTransfer::Sealed
    })
}

fn boundary(exposed: [bool; 2], activity: f64) -> [MoistureLayeredFace; 2] {
    exposed.map(|exposed| {
        if exposed {
            MoistureLayeredFace::PrescribedRelativeHumidity(activity)
        } else {
            MoistureLayeredFace::Sealed
        }
    })
}

impl SolutionParticipant for LayeredMoistureParticipant<'_> {
    type Field = MoistureLayeredSlab;
    type Report = MoistureLayeredReport;

    fn id(&self) -> u64 {
        self.id
    }
    fn time_s(&self) -> f64 {
        self.moisture.time_s()
    }
    fn temperature_kelvin(&self) -> f64 {
        self.moisture.layers()[0].sorption().temperature_kelvin()
    }
    fn cells(&self) -> usize {
        self.moisture.activity().len()
    }
    fn sealed(&self) -> bool {
        exposure(self.exposed, self.moisture.exterior_transfer()) == [false; 2]
    }
    fn admit(&self, cx: &Cx<'_>) -> Result<(), SolutionMoistureError> {
        if !self.sealed() {
            for layer in self.moisture.layers() {
                poll(cx)?;
                // This bridge supplies temperature and water activity only.
                // Extra pressure/composition coordinates need an explicit owner.
                admit_environment(layer.sorption().source_point())?;
            }
        }
        Ok(())
    }
    fn advance(
        &self,
        dt_s: f64,
        _temperature: f64,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, Self::Report), SolutionMoistureError> {
        let faces = boundary(
            exposure(self.exposed, self.moisture.exterior_transfer()),
            activity,
        );
        Ok(self
            .moisture
            .advance(dt_s, &faces, self.transport, self.nonlinear, cx)?)
    }
    fn water_change(&self, next: &Self::Field) -> f64 {
        crate::moisture::layered::atmosphere::water_change(self.moisture, next)
    }
    fn water_mass_kg(&self, endpoint: Option<&Self::Field>) -> f64 {
        let moisture = endpoint.unwrap_or(self.moisture);
        water_sum(moisture.layers().iter().map(|layer| {
            let field = layer.slab().field();
            water_sum(
                field
                    .moisture
                    .iter()
                    .enumerate()
                    .map(|(i, m)| m * field.dry_mass_at(i)),
            )
        }))
    }
}

impl SolutionParticipant for LayeredRadialMoistureParticipant<'_> {
    type Field = MoistureLayeredRadial;
    type Report = MoistureLayeredReport;

    fn id(&self) -> u64 {
        self.id
    }
    fn time_s(&self) -> f64 {
        self.moisture.time_s()
    }
    fn temperature_kelvin(&self) -> f64 {
        self.moisture.layers()[0].sorption().temperature_kelvin()
    }
    fn cells(&self) -> usize {
        self.moisture.activity().len()
    }
    fn sealed(&self) -> bool {
        exposure(self.exposed, self.moisture.exterior_transfer()) == [false; 2]
    }
    fn admit(&self, cx: &Cx<'_>) -> Result<(), SolutionMoistureError> {
        if self.exposed[0] && self.moisture.layers()[0].radial().geometry().faces_m[0] == 0.0 {
            return Err(SolutionMoistureError::Invalid(
                "a solid radial center has zero area and cannot contact a finite solution vapor path",
            ));
        }
        if !self.sealed() {
            for layer in self.moisture.layers() {
                poll(cx)?;
                admit_environment(layer.sorption().source_point())?;
            }
        }
        Ok(())
    }
    fn advance(
        &self,
        dt_s: f64,
        _temperature: f64,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, Self::Report), SolutionMoistureError> {
        let faces = boundary(
            exposure(self.exposed, self.moisture.exterior_transfer()),
            activity,
        );
        Ok(self
            .moisture
            .advance(dt_s, &faces, self.transport, self.nonlinear, cx)?)
    }
    fn water_change(&self, next: &Self::Field) -> f64 {
        crate::moisture::radial::layered::atmosphere::water_change(self.moisture, next)
    }
    fn water_mass_kg(&self, endpoint: Option<&Self::Field>) -> f64 {
        let moisture = endpoint.unwrap_or(self.moisture);
        water_sum(moisture.layers().iter().flat_map(|layer| {
            layer
                .radial()
                .moisture()
                .iter()
                .zip(layer.radial().dry_masses_kg())
                .map(|(&m, &mass)| m * mass)
        }))
    }
}

fn single<Field, Report>(
    candidate: SharedSolutionMoistureCandidate<Field, Report>,
    cx: &Cx<'_>,
) -> Result<SolutionMoistureCandidate<Field, Report>, SolutionMoistureError> {
    let (reservoir, mut solids, shared) = candidate.into_parts();
    let solid = solids
        .pop()
        .expect("one admitted solution-coupled layered body");
    let step = SolutionMoistureStep {
        transport: solid.transport,
        water_activity: shared.water_activity,
        osmotic_pressure_pa: shared.osmotic_pressure_pa,
        solvent_mass_kg: shared.solvent_mass_kg,
        water_to_solid_kg: shared.water_to_solids_kg,
        water_balance_residual_kg: shared.water_balance_residual_kg,
        iterations: shared.iterations,
    };
    poll(cx)?;
    Ok(SolutionMoistureCandidate {
        reservoir,
        moisture: solid.moisture,
        step,
    })
}

impl FiniteMoistureSolution {
    /// Exchange water with one entire heterogeneous slab. Original nonlinear
    /// sorption, diffusivity and explicit exterior/interface resistances govern
    /// every volume trial. No additional isotherm overrides the end materials.
    /// Uptake concentrates and release dilutes the actual finite solution.
    /// Both original inventories publish only after the common balance checks.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered(
        &self,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<LayeredSolutionMoistureCandidate, SolutionMoistureError> {
        let participant = LayeredMoistureParticipant {
            id: 0,
            moisture,
            exposed,
            transport,
            nonlinear,
        };
        single(
            self.advance_shared_layered(
                &[participant],
                dt_s,
                moisture.activity().len(),
                controls,
                cx,
            )?,
            cx,
        )
    }

    /// Simultaneous uptake/desorption of heterogeneous slabs sharing one
    /// original solution. Stable IDs define result and compensated summation
    /// order. Every trial restarts all native stacks and the original solution
    /// phases. Accepted state closes both local-increment and directly stored
    /// water accounts, retaining each layer's actual dry mass and source.
    pub fn advance_shared_layered(
        &self,
        participants: &[LayeredMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredSolutionMoistureCandidate, SolutionMoistureError> {
        self.advance_participants(participants, dt_s, maximum_total_cells, controls, cx)
    }

    /// The same finite aqueous feedback for a coated cylinder, annulus, sphere
    /// or spherical shell. Real unequal shell storage and curved resistance
    /// remain native. Inner/outer exposure is explicit; a solid center cannot
    /// be exposed. Interior material joins never become solution boundaries.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_radial(
        &self,
        moisture: &MoistureLayeredRadial,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<LayeredRadialSolutionMoistureCandidate, SolutionMoistureError> {
        let participant = LayeredRadialMoistureParticipant {
            id: 0,
            moisture,
            exposed,
            transport,
            nonlinear,
        };
        single(
            self.advance_shared_layered_radial(
                &[participant],
                dt_s,
                moisture.activity().len(),
                controls,
                cx,
            )?,
            cx,
        )
    }

    /// Share one actual finite solution among independently sourced coated
    /// radial bodies. All original temperature, phase, source, cell-count and
    /// scalar-work admissions remain in force. Exterior seals permit internal
    /// redistribution and reservoir phase evolution without imposing bath
    /// activity on sealed materials. Failure returns no body or bath candidate.
    ///
    /// These are backward-Euler water endpoints. Selected kinetic phases retain
    /// their original first-order endpoint-volume approximation; no kinetic
    /// SDIRK accuracy is inferred. Temperature/reference geometry stay fixed.
    /// Nonvolatile solute remains in the solution; the equilibrated vapor bridge
    /// does not supply liquid immersion, sorption heat or osmotic traction.
    pub fn advance_shared_layered_radial(
        &self,
        participants: &[LayeredRadialMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredRadialSolutionMoistureCandidate, SolutionMoistureError> {
        self.advance_participants(participants, dt_s, maximum_total_cells, controls, cx)
    }
}
