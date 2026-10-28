//! Imposed magnetic fields couple live velocity, electric current, heat and force.
//!
//! The original entry point is explicit; [`midpoint`] additionally closes
//! velocity/temperature source feedback on the same fixed-grid Boussinesq solvers. It is NOT an induction solver or a total-energy time integrator.
//! The instantaneous electrical/Joule/Lorentz work identity is checked, while
//! the actual kinetic-energy update remains the existing flow owner's result.

pub mod midpoint;
pub mod transport;
mod isolated;

use super::{
    ChtError, ElectricHeatError, ElectricHeatInput, FluidTransportState,
    HeatBalance, HeatSourceLimits, MaterialSlab, Slab, ThermalLbm, ThermalLbm3,
    ThermalSourceError, VoxelDomain, admit, candidate, deposited_power,
    invalid, poll, positive, temperatures,
};
use crate::conjugate::{ElectricMotionInput, MovingConductorSolution};
use fs_exec::CancelGate;

/// Accepted fields from one explicit material-bound magnetic/thermal flow step.
#[derive(Debug, Clone, PartialEq)]
pub struct MagneticHeatStep {
    /// Actual current, conductor-frame heat, force and instantaneous work used.
    pub motion: MovingConductorSolution,
    /// Actual source, storage and independent wall exchange from the heat owner.
    pub heat: HeatBalance,
    /// Electrical source power at the accepted initial state times dt [J].
    pub electrode_energy_j: f64,
    /// Sum of cell Joule powers times dt [J]; the only electrical heat source.
    pub joule_energy_j: f64,
    /// Initial-state Lorentz power times dt [J], NOT measured kinetic-energy
    /// change. Negative in magnetic braking; source-lag error is not hidden.
    pub lagged_mechanical_work_j: f64,
    /// Thermal source minus deposited Joule energy [J].
    pub transfer_defect_j: f64,
    /// Stored heat minus wall input plus lagged mechanical work minus electrode
    /// energy [J]. A coupling-accounting check, NOT total-energy conservation.
    pub coupling_defect_j: f64,
}

// Use the existing Guo-corrected velocity convention, including the force field
// retained by the accepted previous step. Do not introduce a second flow state.
trait MovingSlab<const D: usize, const Q: usize>: Slab<D, Q> {
    fn lattice_velocity(&self, cell: usize) -> [f64; 3];
}
impl MovingSlab<2, 5> for ThermalLbm {
    fn lattice_velocity(&self, cell: usize) -> [f64; 3] {
        let u = self.grid.moments(cell).u;
        [u[0], u[1], 0.0]
    }
}
impl MovingSlab<3, 7> for ThermalLbm3 {
    fn lattice_velocity(&self, cell: usize) -> [f64; 3] { self.grid.moments3(cell).1 }
}

fn velocities<S: MovingSlab<D, Q>, const D: usize, const Q: usize>(
    slab: &MaterialSlab<S>, cells: usize, plane: usize, gate: &CancelGate,
) -> Result<Vec<[f64; 3]>, ElectricHeatError> {
    let mut velocity = vec![[0.0; 3]; cells];
    for (cell, target) in velocity.iter_mut().enumerate().take(cells-plane).skip(plane) {
        if cell % 256 == 0 { poll(gate)?; }
        let u = slab.lattice.lattice_velocity(cell);
        let mach = u.iter().map(|v| v*v).sum::<f64>().sqrt() / crate::CS2.sqrt();
        if !mach.is_finite() || mach >= crate::MACH_LIMIT {
            return Err(invalid("magnetic coupling requires the existing low-lattice-Mach range"));
        }
        for a in 0..3 {
            target[a] = u[a] * slab.force_units.velocity_scale;
            if !target[a].is_finite() { return Err(invalid("physical flow velocity is unrepresentable")); }
        }
    }
    Ok(velocity)
}

// Both explicit and midpoint coupling use this exact source and energy owner.
struct MagneticSource {
    motion: MovingConductorSolution,
    power: Vec<f64>,
    electrode_energy_j: f64,
    joule_energy_j: f64,
    mechanical_work_j: f64,
}

fn admit_field<const D: usize>(field_t: &[[f64; 3]], cells: usize) -> Result<(), ElectricHeatError> {
    if field_t.len() != cells || !field_t.iter().flatten().all(|b| b.is_finite()) {
        return Err(invalid("one finite magnetic flux-density vector per cell is required"));
    }
    // A 2-D flow does not evolve out-of-plane momentum or current. Refuse the
    // general vector instead of silently discarding its EMF/force components.
    if D == 2 && field_t.iter().any(|b| b[0] != 0.0 || b[1] != 0.0) {
        return Err(invalid("2D magnetic flow requires an out-of-plane magnetic field"));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn source_at(
    domain: &VoxelDomain, plane: usize, input: ElectricHeatInput<'_>,
    field_t: &[[f64; 3]], temperature: &[f64], velocity: &[[f64; 3]],
    transport: &FluidTransportState, dt: f64, isolated_current_tolerance_a: Option<f64>, gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
) -> Result<MagneticSource, ElectricHeatError> {
    source_at_with_open_faces(domain, plane, input, &[], field_t, temperature, velocity,
        transport, dt, isolated_current_tolerance_a, gate, resolve)
}

#[allow(clippy::too_many_arguments)]
fn source_at_with_open_faces(
    domain: &VoxelDomain, plane: usize, input: ElectricHeatInput<'_>, open_faces: &[[usize; 2]],
    field_t: &[[f64; 3]], temperature: &[f64], velocity: &[[f64; 3]],
    transport: &FluidTransportState, dt: f64, isolated_current_tolerance_a: Option<f64>, gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
) -> Result<MagneticSource, ElectricHeatError> {
    poll(gate)?;
    let cells = domain.cell_count();
    let tolerance_w = input.limits.coupling_tolerance_j / dt;
    if !tolerance_w.is_finite() || (input.limits.coupling_tolerance_j > 0.0 && tolerance_w == 0.0) {
        return Err(invalid("magnetic power tolerance is unrepresentable"));
    }
    let motion_input = ElectricMotionInput {
        velocity_m_s: velocity, magnetic_flux_density_t: field_t,
        max_cells: input.limits.thermal.max_cells, power_tolerance_w: tolerance_w,
    };
    let law = |cell, _, t| {
        if cell < plane || cell >= cells-plane { Ok(None) }
        else { resolve(cell, t, transport) }
    };
    let motion = if let Some(tolerance_a) = isolated_current_tolerance_a {
        input.setup.solve_isolated_cell_motion_topology(
            domain, temperature, input.contacts, open_faces, motion_input, tolerance_a, law, gate,
        )?
    } else {
        input.setup.solve_with_cell_motion_topology(
            domain, temperature, input.contacts, open_faces, motion_input, law, gate,
        )?
    };
    let (power, joule_energy_j) = deposited_power(domain, plane, &motion.electric.joule_w, dt, gate)?;
    let electrode_energy_j: f64 = motion.electric.electrodes.iter()
        .map(|e| (e.voltage_v * e.current_a) * dt).sum();
    let mechanical_work_j = motion.mechanical_power_w * dt;
    let electrical_defect = electrode_energy_j - joule_energy_j - mechanical_work_j;
    let tol = input.limits.coupling_tolerance_j;
    if ![electrode_energy_j, joule_energy_j, mechanical_work_j, electrical_defect]
        .iter().all(|v| v.is_finite()) || electrical_defect.abs() > tol
        || (motion.electric.total_joule_w * dt - joule_energy_j).abs() > tol
        || (motion.electric.delivered_w * dt - electrode_energy_j).abs() > tol
    { return Err(invalid("magnetic electrical/heat/mechanical source energy does not close")); }
    Ok(MagneticSource { motion, power, electrode_energy_j, joule_energy_j, mechanical_work_j })
}

fn combined_forces<const D: usize>(
    external: &[[f64; D]], motion: &MovingConductorSolution, gate: &CancelGate,
) -> Result<Vec<[f64; D]>, ElectricHeatError> {
    let mut forces = external.to_vec();
    for (cell, target) in forces.iter_mut().enumerate() {
        if cell % 256 == 0 { poll(gate)?; }
        for a in 0..D {
            target[a] += motion.lorentz_force_density_n_m3[cell][a];
            if !target[a].is_finite() { return Err(invalid("combined magnetic/external force overflows")); }
        }
    }
    Ok(forces)
}

fn balance(source: MagneticSource, heat: HeatBalance, input: ElectricHeatInput<'_>) -> Result<MagneticHeatStep, ElectricHeatError> {
    let MagneticSource { motion, electrode_energy_j, joule_energy_j,
        mechanical_work_j: lagged_mechanical_work_j, .. } = source;
    let tol = input.limits.coupling_tolerance_j;
    let transfer_defect_j = heat.source_j - joule_energy_j;
    let coupling_defect_j = heat.stored_change_j - heat.bottom_wall_j - heat.top_wall_j
        + lagged_mechanical_work_j - electrode_energy_j;
    if !transfer_defect_j.is_finite() || transfer_defect_j.abs() > tol
        || !coupling_defect_j.is_finite()
        || coupling_defect_j.abs() > input.limits.thermal.balance_tolerance_j + 2.0 * tol
    { return Err(invalid("magnetic heat transfer or coupling accounting failed")); }
    Ok(MagneticHeatStep { motion, heat, electrode_energy_j, joule_energy_j,
        lagged_mechanical_work_j, transfer_defect_j, coupling_defect_j })
}

#[allow(clippy::too_many_arguments)]
fn advance<S: MovingSlab<D, Q>, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, domain: VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, field_t: &[[f64; 3]], external: &[[f64; D]],
    isolated_current_tolerance_a: Option<f64>,
    gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    step: impl FnOnce(&mut MaterialSlab<S>, &[f64], &[[f64; D]], HeatSourceLimits) -> Result<HeatBalance, ThermalSourceError>,
) -> Result<MagneticHeatStep, ElectricHeatError> {
    poll(gate)?;
    let cells = domain.cell_count();
    admit_field::<D>(field_t, cells)?;
    let temperature = temperatures::<S, D, Q>(slab, input.limits.thermal, gate)?;
    let velocity = velocities::<S, D, Q>(slab, cells, plane, gate)?;
    let source = source_at(&domain, plane, input, field_t, &temperature, &velocity,
        &slab.force_units.transport, slab.force_units.timestep_s, isolated_current_tolerance_a, gate, &mut resolve)?;
    let forces = combined_forces(external, &source.motion, gate)?;
    poll(gate)?;
    let mut next = candidate::<S, D, Q>(slab);
    let heat = step(&mut next, &source.power, &forces, input.limits.thermal)?;
    velocities::<S, D, Q>(&next, cells, plane, gate)?;
    let report = balance(source, heat, input)?;
    poll(gate)?;
    *slab = next;
    Ok(report)
}

impl MaterialSlab<ThermalLbm> {
    /// Advance 2-D flow with self-consistent instantaneous motional current,
    /// its Lorentz force, Joule heating, external force [N/m3] and buoyancy.
    /// `field_t` has one xyz magnetic flux-density [T] vector per cell; only
    /// out-of-plane Bz is admitted. Depth [m] controls the electrical geometry
    /// and joule accounting, not the 2-D dynamics. Material queries have the
    /// same source semantics as `step_electrothermal`.
    ///
    /// Temperature and Guo-corrected velocity are read from the accepted state.
    /// Electric sources are held over one thermal step and recomputed on the
    /// next call. Strong magnetic damping needs timestep refinement; successful
    /// admission is not a stability, time-accuracy or total-energy certificate.
    /// No implicit force persists and no trial resets accepted field history.
    ///
    /// The original state survives material, field, balance, lattice-Mach and
    /// cancellation refusals. Electrical outer faces remain insulating even
    /// where lateral flow is periodic. Ideal electrodes are not obstacles;
    /// every conducting component needs a voltage reference. B is prescribed:
    /// no induction, magnetic source dynamics, Hall or magnetization forces.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 2]], depth_m: f64, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticHeatStep, ElectricHeatError> {
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, forces_n_m3, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) { return Err(invalid("invalid represented 2D cell volume")); }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ])?;
        advance::<ThermalLbm, 2, 5>(self, domain, plane, input, field_t, forces_n_m3, None, gate, resolve,
            |s, p, f, limits| s.step_with_power_density(p, f, depth_m, limits))
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional counterpart of `MaterialSlab::step_magnetothermal`:
    /// all three components of u, B and Lorentz force participate. Uses the
    /// same retained material snapshot, current solver, source conversion,
    /// transactional thermal step and explicit-coupling limitations.
    pub fn step_magnetothermal(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 3]], gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticHeatStep, ElectricHeatError> {
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, forces_n_m3, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume")); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx)?;
        advance::<ThermalLbm3, 3, 7>(self, domain, plane, input, field_t, forces_n_m3, None, gate, resolve,
            |s, p, f, limits| s.step_with_power_density(p, f, limits))
    }
}
