//! Captured-SI plate motion for the original thermal and rheological flow owner.

use super::{MaterialLbmError, MaterialSlab, ThermalLbm3};
use crate::d3q19::{PlateStepError3, WallMomentumExchange3};

/// Native mechanical exchange on one prescribed plate over one accepted step.
/// These SI values do not include thermal wall exchange or molecular heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPlateExchange {
    /// Fluid impulse delivered to the plate [N s].
    pub wall_impulse_n_s: [f64; 3],
    /// Change of fluid population momentum across its reflected links [N s].
    pub fluid_population_impulse_n_s: [f64; 3],
    /// Moving-wall mass correction [N s]. Together with the two impulses above,
    /// `wall + fluid_population = wall_velocity_mass` to summation roundoff.
    pub wall_velocity_mass_impulse_n_s: [f64; 3],
    /// Wall angular impulse about the caller's physical origin [N m s].
    pub wall_angular_impulse_n_m_s: [f64; 3],
    /// Fluid population angular impulse about that same origin [N m s].
    pub fluid_population_angular_impulse_n_m_s: [f64; 3],
    /// Moving-wall mass angular correction [N m s].
    pub wall_velocity_mass_angular_impulse_n_m_s: [f64; 3],
    /// Signed mass change on this plate's reflected links [kg]. A complete
    /// uniform tangential plane is impermeable apart from summation roundoff.
    pub fluid_mass_change_kg: f64,
    /// Mechanical work delivered to the plate [J]. Negative means work supplied
    /// to the fluid by its prescribed actuator; it is not deposited as heat here.
    pub wall_work_j: f64,
    /// Number of actual reflected native links included in this receipt.
    pub measured_links: usize,
}

fn invalid(reason: &'static str) -> MaterialLbmError {
    MaterialLbmError::Invalid(reason)
}

fn component(value: f64, scale: f64) -> Result<f64, MaterialLbmError> {
    let result = value * scale;
    if !value.is_finite() || !result.is_finite() || (value != 0.0 && result == 0.0) {
        return Err(invalid("plate momentum/work conversion is unrepresentable"));
    }
    Ok(result)
}

fn vector(value: [f64; 3], scale: f64) -> Result<[f64; 3], MaterialLbmError> {
    let mut result = [0.0; 3];
    for (output, input) in result.iter_mut().zip(value) {
        *output = component(input, scale)?;
    }
    Ok(result)
}

fn native_error(error: PlateStepError3) -> MaterialLbmError {
    match error {
        PlateStepError3::WallMotion(reason) => invalid(reason),
        PlateStepError3::Layout => invalid("plate motion requires valid native grid storage"),
        PlateStepError3::WallState { .. } | PlateStepError3::Collision { .. } => {
            invalid("native plate flow refused its physical state")
        }
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// Prescribe bottom/top tangential plate velocities [m/s] using this slab's
    /// original dx/dt. `moment_origin_m` uses the grid's physical cell-index
    /// coordinates: the halfway plates lie at dx/2 and (nz-3/2)*dx.
    ///
    /// Every conversion, native tangency check and Mach bound precedes mutation.
    /// A changed drive affects the next native interval without restarting
    /// populations, force frames, viscosity or thermal history. The same drive
    /// survives the checked material/source candidate snapshots.
    ///
    /// Use [`Self::step_with_material_rheology_and_viscous_heating`] to resolve
    /// the actual temperature/shear material and thermalize its molecular stress
    /// power once. The independently measured wall work is a mechanical input,
    /// so adding it again to the prescribed volumetric heat would double count it.
    /// These are fixed, impermeable plates with externally prescribed motion;
    /// this method does not evolve a finite actuator or moving geometry.
    pub fn set_plate_motion_si(
        &mut self,
        velocities_m_s: [[f64; 3]; 2],
        moment_origin_m: [f64; 3],
    ) -> Result<(), MaterialLbmError> {
        let convert = |value: f64, unit: f64| {
            let result = value / unit;
            if !value.is_finite() || !result.is_finite() || (value != 0.0 && result == 0.0) {
                return Err(invalid(
                    "prescribed plate motion is unrepresentable in lattice units",
                ));
            }
            Ok(result)
        };
        let mut velocities = [[0.0; 3]; 2];
        let mut origin = [0.0; 3];
        for (output, input) in velocities
            .iter_mut()
            .flatten()
            .zip(velocities_m_s.into_iter().flatten())
        {
            *output = convert(input, self.force_units.velocity_scale)?;
        }
        for (output, input) in origin.iter_mut().zip(moment_origin_m) {
            *output = convert(input, self.force_units.cell_size_m)?;
        }
        self.lattice
            .grid
            .set_plate_motion(velocities, origin)
            .map_err(native_error)
    }

    /// Return to stationary reflection without erasing fluid or thermal history.
    /// The previous mechanical plate receipt is invalidated.
    pub fn reset_plate_motion(&mut self) {
        self.lattice.grid.reset_plate_motion();
    }

    /// Last native plate exchange in SI, bottom then top. `None` means no driven
    /// interval has completed since initialization or the last drive change.
    /// Impulses are integrated over one retained dt; divide by that dt for force
    /// or torque. This read neither advances the solver nor deposits heat.
    pub fn last_plate_exchange_si(
        &self,
    ) -> Result<Option<[MaterialPlateExchange; 2]>, MaterialLbmError> {
        let Some(exchange) = self.lattice.grid.last_plate_exchange() else {
            return Ok(None);
        };
        let units = self.force_units;
        let dx = units.cell_size_m;
        let mass = ((units.transport.density_kg_m3() * dx) * dx) * dx;
        let impulse = mass * units.velocity_scale;
        let angular = impulse * dx;
        let work = impulse * units.velocity_scale;
        if [mass, impulse, angular, work]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "physical plate mass, impulse or work unit is unrepresentable",
            ));
        }
        let convert =
            |native: WallMomentumExchange3| -> Result<MaterialPlateExchange, MaterialLbmError> {
                Ok(MaterialPlateExchange {
                    wall_impulse_n_s: vector(native.wall_impulse, impulse)?,
                    fluid_population_impulse_n_s: vector(native.fluid_population_impulse, impulse)?,
                    wall_velocity_mass_impulse_n_s: vector(
                        native.wall_velocity_mass_impulse,
                        impulse,
                    )?,
                    wall_angular_impulse_n_m_s: vector(native.wall_angular_impulse, angular)?,
                    fluid_population_angular_impulse_n_m_s: vector(
                        native.fluid_population_angular_impulse,
                        angular,
                    )?,
                    wall_velocity_mass_angular_impulse_n_m_s: vector(
                        native.wall_velocity_mass_angular_impulse,
                        angular,
                    )?,
                    fluid_mass_change_kg: component(native.fluid_mass_change, mass)?,
                    wall_work_j: component(native.wall_work, work)?,
                    measured_links: native.measured_links,
                })
            };
        Ok(Some([convert(exchange[0])?, convert(exchange[1])?]))
    }
}
