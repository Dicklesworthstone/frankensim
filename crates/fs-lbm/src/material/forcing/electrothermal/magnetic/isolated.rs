//! Opt-in electrode-free magnetic braking on the same thermal-flow owners.

use super::*;
use super::midpoint::{MagneticMidpointError, MagneticMidpointOptions, MagneticMidpointStep};

fn admit_isolated(input: ElectricHeatInput<'_>, tolerance_a: f64) -> Result<(), ElectricHeatError> {
    if !input.setup.electrodes.is_empty() || !tolerance_a.is_finite() || tolerance_a < 0.0 {
        return Err(invalid("isolated magnetic flow requires no electrodes and a finite nonnegative current tolerance"));
    }
    Ok(())
}

impl MaterialSlab<ThermalLbm> {
    /// Electrode-free counterpart of `step_magnetothermal`, using a numerical
    /// potential gauge per conducting component. `input.setup.electrodes` must
    /// be empty. Current continuity at every cell is checked against the explicit
    /// `current_tolerance_a` [A]; heat and force use the same retained material.
    /// Only out-of-plane Bz is admitted. Depth controls represented volume.
    /// Electrical walls insulate and fluid walls retain their thermal conditions.
    /// Sources are lagged, not an exact kinetic-energy or induction evolution.
    /// All failures preserve the original slab; query receipts remain caller-owned.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_isolated(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 2]], depth_m: f64, current_tolerance_a: f64,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticHeatStep, ElectricHeatError> {
        poll(gate)?;
        admit_isolated(input, current_tolerance_a)?;
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, external_n_m3, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) { return Err(invalid("invalid represented 2D cell volume")); }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ])?;
        advance::<ThermalLbm, 2, 5>(self, domain, plane, input, field_t, external_n_m3,
            Some(current_tolerance_a), gate, resolve, |s,p,f,limits| s.step_with_power_density(p,f,depth_m,limits))
    }

    /// Electrode-free counterpart of `step_magnetothermal_midpoint`. All current,
    /// material, force, work and heat consistency checks and transactional failure
    /// semantics are shared; numerical gauges do not create physical electrodes.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_isolated(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 2]], depth_m: f64, current_tolerance_a: f64,
        options: MagneticMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticMidpointStep, MagneticMidpointError> {
        self.step_magnetothermal_midpoint_isolated_with_open_faces(
            input, &[], field_t, external_n_m3, depth_m, current_tolerance_a, options, gate, resolve,
        )
    }

    /// Isolated midpoint feedback across a fixed electrically partitioned slab.
    /// Open faces split potential-gauge components without deleting material,
    /// force or heat storage. Each trial and final current recheck uses the
    /// same insulating faces; current continuity includes every gauge cell.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_isolated_with_open_faces(
        &mut self, input: ElectricHeatInput<'_>, open_faces: &[[usize; 2]], field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 2]], depth_m: f64, current_tolerance_a: f64,
        options: MagneticMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticMidpointStep, MagneticMidpointError> {
        poll(gate)?;
        admit_isolated(input, current_tolerance_a)?;
        options.admit()?;
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, external_n_m3, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) { return Err(invalid("invalid represented 2D cell volume").into()); }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ]).map_err(ElectricHeatError::from)?;
        midpoint::solve::<ThermalLbm, 2, 5>(self, domain, plane, input, open_faces, field_t, external_n_m3,
            options, Some(current_tolerance_a), gate, resolve, |s,p,f,limits| s.step_with_power_density(p,f,depth_m,limits))
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional isolated conducting fluid: all vector components,
    /// existing cubic geometry, per-component gauges, shared source and flow
    /// solvers. No external electrode power, induction or total-energy claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_isolated(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 3]], current_tolerance_a: f64, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticHeatStep, ElectricHeatError> {
        poll(gate)?;
        admit_isolated(input, current_tolerance_a)?;
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, external_n_m3, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume")); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx)?;
        advance::<ThermalLbm3, 3, 7>(self, domain, plane, input, field_t, external_n_m3,
            Some(current_tolerance_a), gate, resolve, |s,p,f,limits| s.step_with_power_density(p,f,limits))
    }

    /// Coupled velocity/temperature midpoint for isolated 3-D conductors. Uses
    /// the same nonlinear equation and local source rechecks as the wired path.
    /// Both initial and trial states retain the existing Guo velocity convention.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_isolated(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 3]], current_tolerance_a: f64,
        options: MagneticMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticMidpointStep, MagneticMidpointError> {
        self.step_magnetothermal_midpoint_isolated_with_open_faces(
            input, &[], field_t, external_n_m3, current_tolerance_a, options, gate, resolve,
        )
    }

    /// Three-dimensional isolated counterpart with open electrical faces.
    /// All velocity/temperature feedback and source-work checks use the actual
    /// disconnected conductor graph. Mechanical and thermal walls are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_isolated_with_open_faces(
        &mut self, input: ElectricHeatInput<'_>, open_faces: &[[usize; 2]], field_t: &[[f64; 3]],
        external_n_m3: &[[f64; 3]], current_tolerance_a: f64,
        options: MagneticMidpointOptions, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<MagneticMidpointStep, MagneticMidpointError> {
        poll(gate)?;
        admit_isolated(input, current_tolerance_a)?;
        options.admit()?;
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, external_n_m3, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume").into()); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx).map_err(ElectricHeatError::from)?;
        midpoint::solve::<ThermalLbm3, 3, 7>(self, domain, plane, input, open_faces, field_t, external_n_m3,
            options, Some(current_tolerance_a), gate, resolve, |s,p,f,limits| s.step_with_power_density(p,f,limits))
    }
}
