//! Imposed low-induction magnetic field on the original liquid/current owners.

use super::*;
use crate::conjugate::ElectricMotionInput;
use crate::material::free_surface::MaterialFreeSurfaceScaling;

/// Prescribed stationary B on the electrical core's complete Cartesian cells.
/// The caller owns field geometry, solenoidal validity and the low-induction
/// approximation. This does not infer permeability, magnetization or a field
/// source from the liquid name; no induction, Hall or Maxwell-traction term.
#[derive(Clone, Debug)]
pub struct FreeSurfaceMagneticInput {
    /// Core-local voxel order, in the original Cartesian velocity frame [T].
    pub magnetic_flux_density_t: Vec<[f64; 3]>,
    /// Largest admitted endpoint-minus-initial core speed-vector change [m/s].
    /// Controls lag in the explicit frozen-velocity current/flow split. This
    /// is not a rigorous time-error bound or an energy-balancing correction.
    pub maximum_velocity_change_m_s: f64,
}

/// Native moving-conductor power/force observation at its actual sampled u.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceMagneticState {
    pub velocity_m_s: Vec<[f64; 3]>,
    pub lorentz_force_density_n_m3: Vec<[f64; 3]>,
    /// Signed F dot u; negative magnetic braking becomes positive I^2 R.
    pub mechanical_power_w: f64,
    /// Native transpose-work residual, independently summed from motional EMFs.
    pub transfer_defect_w: f64,
}

/// Actual endpoint motion, separate from frozen-velocity electrical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceMagneticMotion {
    pub maximum_velocity_change_m_s: f64,
    /// Actual core control-volume kinetic storage change [J]. Includes native
    /// advection, pressure, viscosity and other forces; it is NOT Lorentz work.
    /// Neither this value nor its difference from sampled work becomes heat.
    pub core_kinetic_energy_change_j: f64,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Enable native motional EMF and reciprocal Lorentz force on the retained
    /// electrical core. The original current interpolation owns BOTH u x B and
    /// its transpose force; no second J x B reconstruction is added.
    ///
    /// Current samples actual accepted Guo-corrected velocity, held through the
    /// source-resolved thermal substeps. Their I^2 R heats the liquid once and
    /// their integrated Lorentz force supplies one time-average native Guo load.
    /// Gravity and buoyancy remain present. Endpoint force-frame restoration
    /// preserves physical velocity/stress and adds no second half-force impulse.
    ///
    /// This is an explicit first-order thermal/current/flow split. The electrical
    /// identity is electrode power = Joule power + sampled Lorentz power; the
    /// native LBM kinetic-energy change is a separate observation, not an exact
    /// global energy claim. Strong braking requires time refinement; the stated
    /// velocity-change gate refuses an excessive lag before any participant
    /// commits. Fixed complete current geometry and all existing source limits
    /// remain mandatory.
    pub fn with_imposed_magnetic_field(mut self, input: FreeSurfaceMagneticInput)
        -> Result<Self, Error>
    {
        let core = self.electrical.as_mut()
            .ok_or(Error::Invalid("imposed B requires an existing electrical core"))?;
        if input.magnetic_flux_density_t.len() != core.coordinates.len()
            || input.magnetic_flux_density_t.iter().flatten().any(|x| !x.is_finite())
            || !input.maximum_velocity_change_m_s.is_finite()
            || input.maximum_velocity_change_m_s < 0.0
        {
            return Err(Error::Invalid("invalid imposed magnetic field or velocity-change bound"));
        }
        core.magnetic = Some(input);
        Ok(self)
    }
}

impl ElectricalCore {
    pub(super) fn velocities(&self, lattice: &FreeSurface3, scaling: &MaterialFreeSurfaceScaling<'_>)
        -> Result<Vec<[f64; 3]>, Error>
    {
        if self.magnetic.is_none() { return Ok(Vec::new()); }
        self.coordinates.iter().map(|&p| {
            let sample = lattice.sample_at(p[0], p[1], p[2])
                .ok_or(Error::Invalid("magnetic core has no native velocity"))?;
            let [x, y, z] = sample.velocity;
            Ok([scaling.velocity_si(x)?, scaling.velocity_si(y)?, scaling.velocity_si(z)?])
        }).collect()
    }

    pub(crate) fn admit_motion(&self, before: &FreeSurface3, after: &FreeSurface3,
        scaling: &MaterialFreeSurfaceScaling<'_>, gate: &CancelGate)
        -> Result<Option<FreeSurfaceMagneticMotion>, Error>
    {
        let Some(input) = &self.magnetic else { return Ok(None); };
        let initial = self.velocities(before, scaling)?;
        let endpoint = self.velocities(after, scaling)?;
        let mut maximum = 0.0_f64;
        let mut kinetic = 0.0;
        for (cell, &p) in self.coordinates.iter().enumerate() {
            poll(gate)?;
            let u = initial[cell];
            let v = endpoint[cell];
            maximum = maximum.max((v[0] - u[0]).hypot(v[1] - u[1]).hypot(v[2] - u[2]));
            let rho0 = before.sample_at(p[0], p[1], p[2]).expect("sampled initial core").density;
            let rho1 = after.sample_at(p[0], p[1], p[2]).expect("admitted endpoint core").density;
            let old = u.iter().map(|x| x * x).sum::<f64>();
            let new = v.iter().map(|x| x * x).sum::<f64>();
            kinetic = finite(kinetic + 0.5 * scaling.mass_unit_kg() * (rho1 * new - rho0 * old))?;
        }
        if !maximum.is_finite() || maximum > input.maximum_velocity_change_m_s {
            return Err(Error::Invalid("magnetic frozen-velocity split exceeded its velocity-change bound"));
        }
        Ok(Some(FreeSurfaceMagneticMotion {
            maximum_velocity_change_m_s: maximum, core_kinetic_energy_change_j: kinetic,
        }))
    }
}

impl ElectricalTimeline {
    pub(super) fn solve_current(&self, temperature: &[f64], materials: &[ResolvedConductorState],
        gate: &CancelGate) -> Result<(ElectricSolution, Option<FreeSurfaceMagneticState>), Error>
    {
        let core = &self.core;
        if let Some(input) = &core.magnetic {
            let result = core.input.setup.solve_with_cell_motion(
                &core.domain, temperature, &core.input.contacts,
                ElectricMotionInput {
                    velocity_m_s: &self.velocity_m_s,
                    magnetic_flux_density_t: &input.magnetic_flux_density_t,
                    max_cells: core.input.maximum_cells,
                    power_tolerance_w: core.input.power_tolerance_w,
                }, |cell, _, _| Ok(Some(materials[cell].resistivity_ohm_m())), gate,
            )?;
            Ok((result.electric, Some(FreeSurfaceMagneticState {
                velocity_m_s: self.velocity_m_s.clone(),
                lorentz_force_density_n_m3: result.lorentz_force_density_n_m3,
                mechanical_power_w: result.mechanical_power_w,
                transfer_defect_w: result.transfer_defect_w,
            })))
        } else {
            Ok((core.input.setup.solve_with_cell_contacts(
                &core.domain, temperature, &core.input.contacts,
                |cell, _, _| Ok(Some(materials[cell].resistivity_ohm_m())), gate,
            )?, None))
        }
    }

    /// Convert the physical interval impulse using original rho0, dx and dt.
    /// The original native force-frame owner composes and applies this once.
    pub(in crate::material::free_surface::enthalpy::conduction) fn average_forces(
        &self, scaling: &MaterialFreeSurfaceScaling<'_>,
    ) -> Result<Option<BTreeMap<MassLocation3, [f64; 3]>>, Error> {
        if self.core.magnetic.is_none() { return Ok(None); }
        let mut forces = BTreeMap::new();
        for (cell, (&p, impulse)) in self.core.coordinates.iter()
            .zip(&self.lorentz_impulse_n_s).enumerate()
        {
            let mut force = [0.0; 3];
            for axis in 0..3 {
                let acceleration = finite(impulse[axis] / scaling.timestep_s()
                    / self.core.domain.volume(cell) / scaling.material().density_kg_m3())?;
                force[axis] = scaling.gravity_lattice(acceleration)?;
                if impulse[axis] != 0.0 && force[axis] == 0.0 {
                    return Err(Error::Invalid("magnetic force underflows native lattice units"));
                }
            }
            let location = FreeSurface3::transport_location(p)
                .ok_or(Error::Invalid("magnetic core location is unrepresentable"))?;
            forces.insert(location, force);
        }
        Ok(Some(forces))
    }
}
