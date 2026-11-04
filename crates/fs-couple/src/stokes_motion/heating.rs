//! Fixed-reference Stokes viscous heating of the original finite thermal cells.
//!
//! The native accepted endpoint supplies `2 mu D:D`; its per-triangle rates
//! are held over the same motion interval while the existing thermal owner
//! subcycles its authored sources, contacts, reservoirs and events. This is a
//! first-order, one-way split with frozen transport. It adds no advection,
//! pressure work, temperature-dependent viscosity or global passivity claim.
//! Numerical flow damping, body dashpots and wall-slip losses are not deposited.

use std::collections::BTreeMap;

use fs_material::fluid_transport::FluidTransport;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig,
};

use super::{
    GuidedStokesControls, GuidedStokesError, GuidedStokesMotion, GuidedStokesStep,
    MaterialStokesBoundary, nonnegative, poll, positive,
};

/// Explicit correspondence to the original x/y mesh extruded along world +z.
/// Each triangle owns exactly one distinct active thermal cell. Additional
/// thermal cells may carry supports and heat sinks through original links.
#[derive(Clone, Copy, Debug)]
pub struct StokesThermalCells<'a> {
    /// Thermal cell IDs in original `TriMesh::tris` order.
    pub triangle_cells: &'a [CellId],
    /// Extrusion start plane [m]; the end plane is this plus the motion depth.
    pub plane_z_m: f64,
}

/// A refused coupled step changes neither accepted owner.
#[derive(Clone, Debug, PartialEq)]
pub enum GuidedStokesHeatError {
    /// Incomplete, ambiguous, geometrically inconsistent or unsynchronized map.
    Binding(&'static str),
    /// Original flow/body refusal, including cancellation between flow solves.
    Motion(GuidedStokesError),
    /// Original thermal material, budget, energy or cancellation refusal.
    Thermal(DepositionError),
}

impl std::fmt::Display for GuidedStokesHeatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Stokes thermal coupling refusal: {self:?}")
    }
}
impl std::error::Error for GuidedStokesHeatError {}
impl From<GuidedStokesError> for GuidedStokesHeatError {
    fn from(error: GuidedStokesError) -> Self {
        Self::Motion(error)
    }
}
impl From<DepositionError> for GuidedStokesHeatError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// Receipts from the two original owners, committed at one identical time.
pub struct GuidedStokesHeatStep {
    /// Original body/fluid endpoint and its separate work accounts.
    pub motion: GuidedStokesStep,
    /// Native thermal audit, including viscous heat and all authored sources.
    pub thermal: EnergyAudit,
    /// Physical bulk heat actually supplied during native thermal substeps [J].
    /// Already included in `thermal.external_heat_j`; never add it twice.
    pub deposited_viscous_heat_j: f64,
}

impl<F: FluidTransport + ?Sized> GuidedStokesMotion<'_, '_, '_, F> {
    /// Advance native guided flow/body motion and its finite thermal receivers.
    ///
    /// A receiver's centreline must run from the original triangle centroid at
    /// `plane_z_m` to that centroid plus the original depth along +z. Its volume
    /// and total surface area must describe that same triangular prism; mass
    /// must equal the frozen fluid density times that volume. Admission allows
    /// only floating-point roundoff, not a configurable geometric tolerance.
    /// The caller retains responsibility for choosing the same physical fluid's
    /// caloric law; a transport snapshot alone supplies no caloric identity.
    /// No thermal law, boundary source or link is replaced by this composition.
    ///
    /// Heat uses the original endpoint strain field throughout this interval.
    /// Refine both motion and thermal timesteps for temporal accuracy. Both
    /// clocks must initially agree exactly; the native thermal timeline advances
    /// to the exact proposed motion endpoint, including its scheduled events.
    ///
    /// # Errors
    /// Mapping, flow, mechanics, thermal-domain, energy, budget and cancellation
    /// failures preserve both owners, even after successful thermal substeps.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_heating(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        thermal_config: StepConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GuidedStokesHeatStep, GuidedStokesHeatError> {
        poll(&mut cancelled)?;
        self.check_thermal_cells(thermal, cells, &mut cancelled)?;
        let (state, motion) = self.propose_step(
            timestep_s,
            force_density,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            controls,
            &mut cancelled,
        )?;
        let mut source = ViscousHeat {
            watts: BTreeMap::new(),
            heat_j: 0.0,
        };
        let heating = motion.flow.viscous_dissipation().cell_heating_w_m3();
        if heating.len() != cells.triangle_cells.len() {
            return Err(GuidedStokesHeatError::Binding(
                "incomplete native strain field",
            ));
        }
        for (index, (&id, &rate)) in cells.triangle_cells.iter().zip(heating).enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            let watts = rate * self.system.mesh.areas[index] * self.config.depth_m;
            if !nonnegative(watts) {
                return Err(DepositionError::Arithmetic.into());
            }
            source.watts.insert(id, watts);
        }
        let (candidate, source, audit) =
            thermal.advance_observed(&source, state.time_s, thermal_config, || {
                if cancelled() {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })?;
        // This checks conversion and substep integration, not a new fluid
        // energy theorem. Thermal storage has its own original energy gate.
        if candidate.time_s() != state.time_s
            || (source.heat_j - motion.bulk_viscous_dissipation_j).abs()
                > thermal_config.energy_tolerance_j
        {
            return Err(DepositionError::EnergyBalance(
                source.heat_j - motion.bulk_viscous_dissipation_j,
            )
            .into());
        }
        poll(&mut cancelled)?;
        self.state = state;
        *thermal = candidate;
        Ok(GuidedStokesHeatStep {
            motion,
            thermal: audit,
            deposited_viscous_heat_j: source.heat_j,
        })
    }

    fn check_thermal_cells(
        &self,
        thermal: &DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GuidedStokesHeatError> {
        let mesh = self.system.mesh;
        let depth = self.config.depth_m;
        let end_z = cells.plane_z_m + depth;
        if thermal.time_s() != self.state.time_s
            || cells.triangle_cells.len() != mesh.tris.len()
            || !cells.plane_z_m.is_finite()
            || !end_z.is_finite()
            || !same(end_z - cells.plane_z_m, depth)
        {
            return Err(GuidedStokesHeatError::Binding(
                "matching clocks, complete triangle map and representable extrusion required",
            ));
        }
        let mut seen = std::collections::BTreeSet::new();
        for (index, &id) in cells.triangle_cells.iter().enumerate() {
            if index % 256 == 0 {
                poll(cancelled)?;
            }
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(GuidedStokesHeatError::Binding(
                    "inactive thermal triangle receiver",
                ))?;
            if !seen.insert(id) {
                return Err(GuidedStokesHeatError::Binding(
                    "duplicate thermal triangle receiver",
                ));
            }
            let [x, y] = mesh.centroids[index];
            let start = [x, y, cells.plane_z_m];
            let end = [x, y, end_z];
            let volume = mesh.areas[index] * depth;
            let mass = self.fluid.transport().density_kg_m3() * volume;
            let perimeter: f64 = mesh.tri_edges[index]
                .iter()
                .map(|&(edge, _)| mesh.edges[edge].len)
                .sum();
            let surface = 2.0 * mesh.areas[index] + perimeter * depth;
            if !positive(volume)
                || !positive(mass)
                || !positive(surface)
                || !same(cell.geometry.volume_m3, volume)
                || !same(cell.geometry.surface_area_m2, surface)
                || !same(cell.mass_kg, mass)
                || !(0..3).all(|i| {
                    same(cell.geometry.start_m[i], start[i]) && same(cell.geometry.end_m[i], end[i])
                })
            {
                return Err(GuidedStokesHeatError::Binding(
                    "thermal geometry and mass must match the original fluid triangle prism",
                ));
            }
        }
        Ok(())
    }
}

fn same(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 64.0 * f64::EPSILON * a.abs().max(b.abs())
}

#[derive(Clone)]
struct ViscousHeat {
    watts: BTreeMap<CellId, f64>,
    heat_j: f64,
}

impl DepositionObserver for ViscousHeat {
    type Error = DepositionError;

    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        Ok(self.watts.clone())
    }

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        // Native events can change boundaries/links and add other cells, but
        // cannot replace an existing receiver's mass, material or geometry.
        Ok(())
    }

    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        for (index, watts) in self.watts.values().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            self.heat_j += dt * watts;
            if !self.heat_j.is_finite() {
                return Err(DepositionError::Arithmetic);
            }
        }
        Ok(())
    }
}
