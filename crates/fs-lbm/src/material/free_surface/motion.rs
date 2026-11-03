//! SI wall motion and reaction through the original frozen liquid trajectory.

use super::*;
use crate::d3q19::freesurface3::pressure::GasPressureField3;
use crate::d3q19::freesurface3::{FreeSurfaceError3, PlanarWallMotionReport3};

mod rigid;
pub use rigid::{MaterialPlanarRigidWallMotion3, MaterialPlanarRigidWallReport3};
mod field;
pub use field::{
    MaterialPlanarWallFaceExchange3, MaterialPlanarWallFieldReport3,
    MaterialPlanarWallVelocityField3, MaterialWallWrenchExchange3,
};

/// One accepted isothermal wall-driven free-surface timestep in SI units.
///
/// The paired raw impulses retain the numerical LBM pressure datum for momentum
/// accounting. Physical wall loads are reported separately about the declared
/// reference; an absolute load is available only with with_pressure_reference.
/// Wall work is tangential and therefore independent of the pressure datum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPlanarWallMotionReport3 {
    /// Domain plane whose material moves tangentially.
    pub wall: PlanarWall3,
    /// Prescribed velocity in the existing liquid inertial frame [m/s].
    pub velocity_m_s: [f64; 3],
    /// Actual native population impulse into liquid [kg m/s], including datum.
    pub fluid_impulse_kg_m_s: [f64; 3],
    /// Exactly reciprocal raw impulse on the wall [kg m/s], including datum.
    pub wall_impulse_kg_m_s: [f64; 3],
    /// Numerical density-one pressure impulse into liquid [kg m/s].
    pub reference_pressure_fluid_impulse_kg_m_s: [f64; 3],
    /// Wall force about the declared reference pressure [N].
    pub pressure_perturbation_wall_force_n: [f64; 3],
    /// Wall force including the explicitly declared absolute reference [N].
    /// None means no absolute pressure was supplied; no pressure is invented.
    pub absolute_wall_force_n: Option<[f64; 3]>,
    /// Mechanical energy delivered by the prescribed wall to liquid [J].
    pub fluid_work_j: f64,
    /// Equal and opposite energy received by the external wall actuator [J].
    pub wall_work_j: f64,
    /// Complete flooded plane area [m²].
    pub wetted_area_m2: f64,
    /// Native total liquid mass before and after the accepted step [kg].
    pub mass_before_kg: f64,
    /// Native total liquid mass after the accepted step [kg].
    pub mass_after_kg: f64,
    /// Absolute relative native liquid mass defect.
    pub relative_mass_defect: f64,
    /// Maximum physical liquid speed at the two admitted endpoints [m/s].
    pub maximum_speed_m_s: f64,
}

fn vector_si(
    value: [f64; 3],
    unit: f64,
    reason: &'static str,
) -> Result<[f64; 3], MaterialLbmError> {
    let mut result = [0.0; 3];
    for (out, value) in result.iter_mut().zip(value) {
        *out = finite(value * unit, reason)?;
    }
    Ok(result)
}

impl MaterialFreeSurface3<'_> {
    /// Propose one original wall-driven physical step without changing this liquid.
    ///
    /// Uses the identical native candidate and complete SI/source/speed/mass
    /// admission as [`Self::step_with_planar_wall_motion_si`]. Coupled mechanics
    /// can therefore test different wall velocities from one accepted state
    /// without first copying its full population history a second time.
    /// The returned candidate retains the original immutable material/scaling.
    ///
    /// # Errors
    /// All original fully flooded tangent-wall, source-domain, native-step,
    /// representability, mass and cancellation refusals remain operative.
    pub fn planar_wall_motion_candidate_si(
        &self,
        wall: PlanarWall3,
        velocity_m_s: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<(Self, MaterialPlanarWallMotionReport3), MaterialLbmError> {
        self.planar_wall_motion_candidate_inner(
            wall,
            velocity_m_s,
            maximum_relative_mass_defect,
            gate,
            None,
        )
    }

    /// Advance the original liquid with one fully flooded moving planar wall.
    ///
    /// Velocity is prescribed tangentially in the same inertial frame as
    /// sample_at_si, constant over the original physical timestep. A moving
    /// bottom below a free surface is supported. Normal motion, cell
    /// cover/vacate, gas or partial fill on the selected plane and a moving
    /// contact line refuse. The caller owns the substrate/actuator dynamics.
    ///
    /// Original source density, viscosity, capillarity, force and native VOF
    /// transport remain authoritative. Actual link reflection supplies
    /// reciprocal impulse and mechanical wall work; no second fluid integrator
    /// or thermal dissipation model is introduced. This frozen isothermal
    /// method does not deposit wall work as caloric heat.
    ///
    /// The same private native candidate is committed only after source-scale,
    /// endpoint speed, SI representability, mass and cancellation checks.
    ///
    /// # Errors
    /// Refuses wall or liquid speed above the original scenario bound,
    /// unsupported native wall state, excessive mass defect, native failure,
    /// unrepresentable output, or observed cancellation without publication.
    pub fn step_with_planar_wall_motion_si(
        &mut self,
        wall: PlanarWall3,
        velocity_m_s: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<MaterialPlanarWallMotionReport3, MaterialLbmError> {
        self.step_with_planar_wall_motion_inner(
            wall,
            velocity_m_s,
            maximum_relative_mass_defect,
            gate,
            None,
        )
    }

    /// Advance gas pressure [Pa] and tangential wall motion [m/s] together.
    ///
    /// Requires an explicit fixed [`MaterialFreeSurfaceScaling::with_pressure_reference`].
    /// The callback supplies positive absolute gas pressure at exposed liquid
    /// cell centers [m], frozen for this physical timestep. Original capillary
    /// traction adds to this load under the supplied density-excursion bound.
    /// Wall reflection and gas reconstruction use ONE native candidate; the
    /// reported reaction and mechanical work come from its actual populations.
    /// Pressure disturbances propagate through that original liquid evolution.
    ///
    /// The source coefficients remain frozen and wall motion remains tangential
    /// on a fully flooded, fixed domain plane. This neither evolves a gas
    /// reservoir nor deposits boundary work as thermal energy. The caller owns
    /// the prescribed gas field and wall actuator. A common shift of the gas
    /// pressure and reference leaves the trajectory and wall work unchanged;
    /// the reported absolute normal wall force includes the shifted datum.
    ///
    /// # Errors
    /// The original wall, source-scale, speed, mass and SI-output gates apply.
    /// Missing reference, invalid pressure, excess combined gas/capillary
    /// density excursion, native failure or cancellation preserves the entire
    /// accepted liquid state.
    pub fn step_with_planar_wall_motion_and_gas_pressure_si(
        &mut self,
        wall: PlanarWall3,
        velocity_m_s: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 3]) -> f64,
    ) -> Result<MaterialPlanarWallMotionReport3, MaterialLbmError> {
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        let reference = self
            .scaling
            .reference_pressure_pa
            .ok_or(MaterialLbmError::Invalid(
                "gas-loaded moving wall requires an explicit fixed pressure reference",
            ))?;
        let pressure = super::pressure::lattice_pressure_field(
            reference,
            self.scaling.cell_size_m,
            self.scaling.pressure_unit_pa,
            pressure_pa,
        );
        self.step_with_planar_wall_motion_inner(
            wall,
            velocity_m_s,
            maximum_relative_mass_defect,
            gate,
            Some((
                maximum_density_excursion,
                GasPressureField3::Cell(&pressure),
            )),
        )
    }

    /// Advance one wall-driven physical step with distinct pressure per gas link.
    ///
    /// Positive absolute pressure [Pa] uses [`MaterialGasPressureLink3`]'s cell
    /// coordinates and midpoint, with the fixed construction-time reference.
    /// The original physical scaling, flooded tangent-only wall, reciprocal
    /// wall impulse/work, endpoint speed admission and atomic publication of
    /// [`Self::step_with_planar_wall_motion_and_gas_pressure_si`] all apply.
    /// No moving gas geometry or gas pressure-work owner is introduced.
    pub fn step_with_planar_wall_motion_and_gas_link_pressure_si(
        &mut self,
        wall: PlanarWall3,
        velocity_m_s: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(MaterialGasPressureLink3) -> f64,
    ) -> Result<MaterialPlanarWallMotionReport3, MaterialLbmError> {
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        let reference = self
            .scaling
            .reference_pressure_pa
            .ok_or(MaterialLbmError::Invalid(
                "gas-loaded moving wall requires an explicit fixed pressure reference",
            ))?;
        let field = super::pressure::lattice_link_pressure_field(
            reference,
            self.scaling.cell_size_m,
            self.scaling.pressure_unit_pa,
            pressure_pa,
        );
        self.step_with_planar_wall_motion_inner(
            wall,
            velocity_m_s,
            maximum_relative_mass_defect,
            gate,
            Some((maximum_density_excursion, GasPressureField3::Link(&field))),
        )
    }

    fn step_with_planar_wall_motion_inner(
        &mut self,
        wall: PlanarWall3,
        velocity_m_s: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
        gas_pressure: Option<(f64, GasPressureField3<'_>)>,
    ) -> Result<MaterialPlanarWallMotionReport3, MaterialLbmError> {
        let (candidate, report) = self.planar_wall_motion_candidate_inner(
            wall,
            velocity_m_s,
            maximum_relative_mass_defect,
            gate,
            gas_pressure,
        )?;
        *self = candidate;
        Ok(report)
    }

    fn planar_wall_motion_candidate_inner(
        &self,
        wall: PlanarWall3,
        velocity_m_s: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
        gas_pressure: Option<(f64, GasPressureField3<'_>)>,
    ) -> Result<(Self, MaterialPlanarWallMotionReport3), MaterialLbmError> {
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        let scale = &self.scaling;
        let velocity = velocity_m_s.map(|v| v / scale.velocity_scale);
        let wall_speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
        if !wall_speed.is_finite() || wall_speed > scale.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "moving-wall velocity exceeds the original material scenario speed",
            ));
        }
        let units = self.planar_wall_units()?;
        let (candidate, raw): (_, PlanarWallMotionReport3) =
            self.lattice.planar_wall_motion_candidate(
                wall,
                velocity,
                maximum_relative_mass_defect,
                gate,
                gas_pressure,
            )?;
        if raw.maximum_speed > scale.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "wall-driven liquid exceeds the original material scenario speed",
            ));
        }
        let report = self.planar_wall_report_si(raw, velocity_m_s, units)?;
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        Ok((
            Self {
                lattice: candidate,
                scaling: self.scaling,
            },
            report,
        ))
    }

    // One conversion and pressure-reference owner for both uniform motion and
    // the linear part of an actual rigid-wall wrench.
    fn planar_wall_units(&self) -> Result<[f64; 3], MaterialLbmError> {
        let scale = &self.scaling;
        let impulse_unit = positive(
            scale.mass_unit_kg * scale.velocity_scale,
            "moving-wall SI impulse unit is unrepresentable",
        )?;
        let force_unit = positive(
            impulse_unit / scale.timestep_s,
            "moving-wall SI force unit is unrepresentable",
        )?;
        let work_unit = positive(
            impulse_unit * scale.velocity_scale,
            "moving-wall SI work unit is unrepresentable",
        )?;
        Ok([impulse_unit, force_unit, work_unit])
    }

    fn planar_wall_report_si(
        &self,
        raw: PlanarWallMotionReport3,
        velocity_m_s: [f64; 3],
        [impulse_unit, force_unit, work_unit]: [f64; 3],
    ) -> Result<MaterialPlanarWallMotionReport3, MaterialLbmError> {
        let scale = &self.scaling;
        let fluid_impulse_kg_m_s = vector_si(
            raw.fluid_impulse,
            impulse_unit,
            "moving-wall SI fluid impulse is unrepresentable",
        )?;
        let wall_impulse_kg_m_s = fluid_impulse_kg_m_s.map(|v| -v);
        let reference_pressure_fluid_impulse_kg_m_s = vector_si(
            raw.reference_pressure_fluid_impulse,
            impulse_unit,
            "moving-wall SI reference impulse is unrepresentable",
        )?;
        let pressure_perturbation_wall_force_n = vector_si(
            core::array::from_fn(|d| raw.wall_impulse[d] + raw.reference_pressure_fluid_impulse[d]),
            force_unit,
            "moving-wall SI pressure-perturbation force is unrepresentable",
        )?;
        let wetted_area_m2 = positive(
            (raw.wetted_area * scale.cell_size_m) * scale.cell_size_m,
            "moving-wall SI wetted area is unrepresentable",
        )?;
        let absolute_wall_force_n = scale
            .reference_pressure_pa
            .map(|pressure| {
                let reference_force = positive(
                    pressure * wetted_area_m2,
                    "moving-wall absolute reference force is unrepresentable",
                )?;
                let normal = raw.reference_pressure_fluid_impulse.map(|v| v.signum());
                let mut absolute = pressure_perturbation_wall_force_n;
                for d in 0..3 {
                    // signum(0.0) is +1, so zero reference components must remain zero.
                    if raw.reference_pressure_fluid_impulse[d] != 0.0 {
                        absolute[d] = finite(
                            absolute[d] - normal[d] * reference_force,
                            "moving-wall absolute force is unrepresentable",
                        )?;
                    }
                }
                Ok::<_, MaterialLbmError>(absolute)
            })
            .transpose()?;
        let fluid_work_j = finite(
            raw.fluid_work * work_unit,
            "moving-wall SI actuator work is unrepresentable",
        )?;
        let report = MaterialPlanarWallMotionReport3 {
            wall: raw.wall,
            velocity_m_s,
            fluid_impulse_kg_m_s,
            wall_impulse_kg_m_s,
            reference_pressure_fluid_impulse_kg_m_s,
            pressure_perturbation_wall_force_n,
            absolute_wall_force_n,
            fluid_work_j,
            wall_work_j: -fluid_work_j,
            wetted_area_m2,
            mass_before_kg: positive(
                raw.mass_before * scale.mass_unit_kg,
                "moving-wall SI initial mass is unrepresentable",
            )?,
            mass_after_kg: positive(
                raw.mass_after * scale.mass_unit_kg,
                "moving-wall SI final mass is unrepresentable",
            )?,
            relative_mass_defect: raw.relative_mass_defect,
            maximum_speed_m_s: finite(
                raw.maximum_speed * scale.velocity_scale,
                "moving-wall SI liquid speed is unrepresentable",
            )?,
        };
        Ok(report)
    }
}
