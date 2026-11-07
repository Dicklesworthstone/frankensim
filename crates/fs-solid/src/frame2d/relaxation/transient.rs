//! Original longitudinal/shear Maxwell histories on the full native trajectory.
//! The shared coupled exact-flow owner advances mechanics, every source memory
//! and every ordinary boundary state together. These observations never replay
//! a constitutive update or turn the integration defect into material heat.

/// Tensile normal stress from the actual accepted equilibrium and Maxwell state.
pub mod strength;

use super::{
    FrameRelaxationCoordinate, FrameRelaxationError, PreparedFrameRelaxation, invalid, poll,
};
use crate::frame2d::coupled::transient::{
    FrameTransientControls, FrameTransientReport, PreparedFrameTransient, finite, sum,
};
use fs_exec::Cx;

type Result<T> = core::result::Result<T, FrameRelaxationError>;

/// Explicit original source history at the initial nodal displacement.
#[derive(Clone, Copy, Debug)]
pub enum FrameRelaxationInitialHistory<'a> {
    /// Every arm is relaxed: zero elastic extension and zero arm storage.
    Relaxed,
    /// Zero viscous extension; each original strain is initially stored in its arm.
    Unrelaxed,
    /// One elastic extension [m] per `PreparedFrameRelaxation::coordinates()`.
    /// Positive arms invisible to free motion still require their own history.
    ElasticExtensions(&'a [f64]),
}

/// One original layer/arm/strain coordinate at an accepted endpoint.
#[derive(Clone, Debug)]
pub struct FrameRelaxationMemory {
    pub coordinate: FrameRelaxationCoordinate,
    pub extension_m: f64,
    pub elastic_extension_m: f64,
    pub viscous_extension_m: f64,
    pub opposing_force_n: f64,
    pub stored_energy_j: f64,
    /// Exact interval integral of this source's z^2/tau, not a storage difference.
    pub dissipated_energy_j: f64,
}

/// Borrowed accepted reports; mechanism heat is included in the complete native
/// motion report and must not be added a second time to its energy balance.
pub struct FrameRelaxationTransientReport<'a> {
    pub motion: &'a FrameTransientReport,
    pub memory: &'a [FrameRelaxationMemory],
    pub material_memory_energy_j: f64,
    pub material_dissipation_j: f64,
    /// Original boundary's own resistance, separately from constitutive memory.
    pub boundary_dissipation_j: f64,
}

/// Fixed-source, fixed-clock native history with atomic memory observations.
/// No mode truncation, fitted damping ratio, split held-velocity forcing or
/// parallel integrator replaces the original longitudinal or shear Maxwell law.
pub struct PreparedFrameRelaxationTransient<'a> {
    source: &'a PreparedFrameRelaxation,
    native: PreparedFrameTransient<'a>,
    memory: Vec<FrameRelaxationMemory>,
    staged_memory: Vec<FrameRelaxationMemory>,
    totals: [f64; 3],
    staged_totals: [f64; 3],
}

impl PreparedFrameRelaxation {
    /// Prepare exact native ringing/relaxation with all original layer arms,
    /// nonproportional frame viscosity and optional ordinary boundary memory.
    /// Explicit broadband material support is required. Descriptor constraints
    /// remain available to AC/poles but refuse this ordinary exact time path.
    /// Starts at zero motion and zero memory; use `initialize` or `restore_state`
    /// for an explicit nonzero history. Temperature and geometry remain fixed.
    pub fn prepare_transient(
        &self,
        duration_s: f64,
        controls: FrameTransientControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedFrameRelaxationTransient<'_>> {
        poll(cx)?;
        let nf = self.coupled.free_dofs().len();
        let mut groups = Vec::with_capacity(self.coordinates.len() + 2);
        groups.push(nf..2 * nf);
        groups.extend(
            self.coordinates
                .iter()
                .map(|c| c.native_state..c.native_state + 1),
        );
        groups.push(self.boundary_state_range());
        let native =
            PreparedFrameTransient::prepare(&self.coupled, duration_s, controls, groups, cx)?;
        let memory: Vec<_> = self
            .coordinates
            .iter()
            .map(|&coordinate| FrameRelaxationMemory {
                coordinate,
                extension_m: 0.0,
                elastic_extension_m: 0.0,
                viscous_extension_m: 0.0,
                opposing_force_n: 0.0,
                stored_energy_j: 0.0,
                dissipated_energy_j: 0.0,
            })
            .collect();
        let staged_memory = memory.clone();
        poll(cx)?;
        Ok(PreparedFrameRelaxationTransient {
            source: self,
            native,
            memory,
            staged_memory,
            totals: [0.0; 3],
            staged_totals: [0.0; 3],
        })
    }
}

impl PreparedFrameRelaxationTransient<'_> {
    pub fn source(&self) -> &PreparedFrameRelaxation {
        self.source
    }
    pub fn duration_s(&self) -> f64 {
        self.native.duration_s()
    }
    pub fn canonical_state(&self) -> &[f64] {
        self.native.canonical_state()
    }
    pub fn report(&self) -> FrameRelaxationTransientReport<'_> {
        FrameRelaxationTransientReport {
            motion: self.native.report(),
            memory: &self.memory,
            material_memory_energy_j: self.totals[0],
            material_dissipation_j: self.totals[1],
            boundary_dissipation_j: self.totals[2],
        }
    }

    /// Original full nodal motion and explicit original arm history. Boundary
    /// states use `boundary_state_range()` order; no existing memory is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn initialize(
        &mut self,
        time_s: f64,
        displacement: &[f64],
        velocity: &[f64],
        history: FrameRelaxationInitialHistory<'_>,
        boundary_state: &[f64],
        cx: &Cx<'_>,
    ) -> Result<()> {
        poll(cx)?;
        let count = self.source.coordinates.len();
        if displacement.len() != self.source.coupled.ndof()
            || displacement.iter().any(|x| !x.is_finite())
            || boundary_state.len() != self.source.boundary_state_range().len()
            || boundary_state.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "initial Maxwell frame and boundary must retain every original coordinate",
            ));
        }
        if let FrameRelaxationInitialHistory::ElasticExtensions(values) = history {
            if values.len() != count || values.iter().any(|x| !x.is_finite()) {
                return Err(invalid(
                    "one finite elastic extension per original Maxwell coordinate required",
                ));
            }
        }
        let mut attached = Vec::with_capacity(count + boundary_state.len());
        for (index, coordinate) in self.source.coordinates.iter().enumerate() {
            poll(cx)?;
            let extension = match history {
                FrameRelaxationInitialHistory::Relaxed => 0.0,
                FrameRelaxationInitialHistory::Unrelaxed => {
                    extension(self.source, coordinate, displacement)?
                }
                FrameRelaxationInitialHistory::ElasticExtensions(values) => values[index],
            };
            attached.push(finite(coordinate.coefficient.root * extension)?);
        }
        attached.extend_from_slice(boundary_state);
        self.native
            .stage_initialize(time_s, displacement, velocity, &attached, cx)?;
        self.observe_staged(cx)?;
        poll(cx)?;
        self.publish_staged();
        Ok(())
    }

    /// Complete exact checkpoint, including dark Maxwell and original boundary
    /// states. A late observation/source/cancellation refusal changes no history.
    pub fn restore_state(&mut self, time_s: f64, state: &[f64], cx: &Cx<'_>) -> Result<()> {
        self.native.stage_restore(time_s, state, cx)?;
        self.observe_staged(cx)?;
        poll(cx)?;
        self.publish_staged();
        Ok(())
    }

    /// One original held nodal force/moment interval. Both native motion and
    /// source-specific storage/heat are staged before either is published.
    pub fn step(
        &mut self,
        nodal_force: &[f64],
        cx: &Cx<'_>,
    ) -> Result<FrameRelaxationTransientReport<'_>> {
        self.native.stage_step(nodal_force, cx)?;
        self.observe_staged(cx)?;
        poll(cx)?;
        self.publish_staged();
        Ok(self.report())
    }

    fn observe_staged(&mut self, cx: &Cx<'_>) -> Result<()> {
        let state = self.native.staged_state();
        let displacement = &self.native.staged_report().displacement;
        let heat = self.native.staged_heat();
        for (index, output) in self.staged_memory.iter_mut().enumerate() {
            poll(cx)?;
            let coordinate = &output.coordinate;
            let value = state[coordinate.native_state];
            output.extension_m = extension(self.source, coordinate, displacement)?;
            output.elastic_extension_m = finite(value / coordinate.coefficient.root)?;
            output.viscous_extension_m = finite(output.extension_m - output.elastic_extension_m)?;
            output.opposing_force_n = finite(coordinate.coefficient.root * value)?;
            output.stored_energy_j = finite((0.5 * value) * value)?;
            output.dissipated_energy_j = heat[index + 1];
        }
        self.staged_totals = [
            sum(self.staged_memory.iter().map(|m| m.stored_energy_j))?,
            sum(self.staged_memory.iter().map(|m| m.dissipated_energy_j))?,
            heat[heat.len() - 1],
        ];
        poll(cx)
    }

    fn publish_staged(&mut self) {
        self.native.publish_staged();
        core::mem::swap(&mut self.memory, &mut self.staged_memory);
        core::mem::swap(&mut self.totals, &mut self.staged_totals);
    }
}

fn extension(
    source: &PreparedFrameRelaxation,
    coordinate: &FrameRelaxationCoordinate,
    displacement: &[f64],
) -> Result<f64> {
    let member = &source.members[coordinate.member];
    Ok(sum(member
        .dofs
        .iter()
        .zip(&member.weights[coordinate.channel.index()])
        .map(|(&dof, &weight)| weight * displacement[dof]))?)
}
