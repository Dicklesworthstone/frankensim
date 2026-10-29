//! Complete three-dimensional material valve images through the SAME aperture.
//! A single explicit change of coordinates makes q = sqrt(M_total) x, where q
//! is the source's mass-normalized coordinate and x is metre-valued. Thus
//! M_x=M_total I, K_x=M_total diag(lambda), C_x=M_total C, and every pressure,
//! force and displacement map receives the same sqrt(M_total) factor.
//! No mode is normalized by its mean opening and no cross-mode loss is dropped.

/// Original hereditary tensor-solid images through the same aperture owner.
pub mod hereditary;
pub use hereditary::{HereditaryModalSolid, HereditarySolidOptions};
/// Original exterior-face contact and slit geometry on the same solid modes.
pub mod surface;
pub use surface::SolidSurfaceClosureSpec;

use super::{ApertureMode, ModalAperture, ModalApertureSpec, invalid, nonlinear};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::paired::{BladePropagator, BladeState, PairedApertureLimits};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_solid::reduce::visco::MaterialModalValve;

/// Explicit clock, load and existing runtime limits for the source image.
/// ExactZoh retains the full frozen material K/C; it does not make the nonlinear
/// slit/contact or sampled pressure/force histories exact.
#[derive(Clone, Copy, Debug)]
pub struct MaterialSolidOptions {
    pub impedance_pa_s_m3: f64,
    pub time_step_s: f64,
    pub max_steps: u64,
    pub limits: PairedApertureLimits,
}

/// Source ownership and the declared physical coordinate conversion. The
/// material state is frozen exactly as admitted by the reduction, not evolved
/// or relabeled as a full hereditary solid.
pub struct MaterialModalSolid {
    source: MaterialModalValve,
    root_mass: f64,
}
impl MaterialModalSolid {
    /// Same mass conversion and original nodal strain map as the solid owner.
    /// The aperture checks one bounded physical endpoint before publication;
    /// this does not certify strain extrema between samples.
    pub(super) fn validate_state(&self, state: &[BladeState]) -> Result<(), AcousticRealizeError> {
        let coordinates = state.iter().map(|value|
            scaled(value.displacement_m, self.root_mass)).collect::<Result<Vec<_>, _>>()?;
        self.source.validate_modal_displacement(&coordinates, || false).map_err(nonlinear)
    }

    pub fn source(&self) -> &MaterialModalValve {
        &self.source
    }
    /// Every physical modal coordinate uses this same positive mass [kg].
    pub fn coordinate_mass_kg(&self) -> f64 {
        self.source.total_mass_kg()
    }
    /// Physical coordinates [m] -> original mass-normalized q [m sqrt(kg)].
    /// The same conversion applies to velocities and preserves kinetic energy.
    pub fn mass_normalized_coordinates(
        &self,
        coordinates: &[f64],
        gate: &CancelGate,
    ) -> Result<Vec<f64>, AcousticRealizeError> {
        poll(gate)?;
        if coordinates.len() != self.source.modes().len()
            || coordinates.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "solid modal observation needs every finite physical coordinate",
            ));
        }
        let result = coordinates
            .iter()
            .map(|v| scaled(*v, self.root_mass))
            .collect::<Result<Vec<_>, _>>()?;
        poll(gate)?;
        Ok(result)
    }
    /// Observe nodal displacement or velocity on the original support map.
    pub fn expand_nodal(
        &self,
        coordinates: &[f64],
        gate: &CancelGate,
    ) -> Result<Vec<[f64; 3]>, AcousticRealizeError> {
        let q = self.mass_normalized_coordinates(coordinates, gate)?;
        self.source
            .expand_nodal(&q, gate)
            .map_err(|error| source_error(error, gate))
    }
}

fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}
fn source_error(error: impl core::fmt::Display, gate: &CancelGate) -> AcousticRealizeError {
    if gate.is_requested() {
        AcousticRealizeError::Cancelled
    } else {
        nonlinear(error)
    }
}
fn scaled(value: f64, scale: f64) -> Result<f64, AcousticRealizeError> {
    let result = value * scale;
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        Err(invalid(
            "complete material modal coordinate conversion is unrepresentable",
        ))
    } else {
        Ok(result)
    }
}

impl ModalAperture {
    /// Consume every mode of an admitted three-dimensional material image,
    /// including zero-opening modes, in the existing ExactZoh/full-C pressure,
    /// swept-volume and contact solve. Geometry and uniform pressure projection
    /// come from that source; no independent scalar damping is synthesized.
    ///
    /// States are metre-valued x with q=sqrt(total_mass_kg)*x. Each optional
    /// force pattern is one dimensionless vector per original node, defining
    /// nodal force per applied newton at the corresponding existing force port.
    /// None means exactly zero force participation. All supplied projections
    /// are pulled back through the same source basis, including dark modes.
    ///
    /// The source's frozen-frequency dispersion report remains the authority;
    /// this adapter does not establish omitted-mode, spatial-contact, broadband
    /// material or large-strain accuracy. The existing angular/energy/flow/
    /// pressure/contact admission and sample-refinement obligations all remain.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_solid(
        source: MaterialModalValve,
        force_patterns_per_n: [Option<&[[f64; 3]]>; 2],
        gas: GasState,
        options: MaterialSolidOptions,
        state: Vec<BladeState>,
        contact: Obstacle,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        gas.validate_coherence().map_err(nonlinear)?;
        let n = source.modes().len();
        if n == 0 || n > super::MAX_MODES || state.len() != n {
            return Err(invalid(
                "complete solid image must fit the aperture without dropping modes",
            ));
        }
        let mass = source.total_mass_kg();
        let root_mass = mass.sqrt();
        if !root_mass.is_finite() || root_mass <= 0.0 {
            return Err(invalid("complete solid modal mass is not representable"));
        }
        let mut force = [vec![0.0; n], vec![0.0; n]];
        for (port, pattern) in force_patterns_per_n.into_iter().enumerate() {
            if let Some(pattern) = pattern {
                force[port] = source
                    .project_nodal_force(pattern, gate)
                    .map_err(|error| source_error(error, gate))?;
            }
        }
        let mut modes = Vec::with_capacity(n);
        for (i, mode) in source.modes().iter().enumerate() {
            poll(gate)?;
            modes.push(ApertureMode {
                mass_kg: mass,
                stiffness_n_m: scaled(mode.lambda, mass)?,
                damping_n_s_m: 0.0,
                pressure_area_m2: scaled(
                    source.pressure_projection_m2_per_sqrt_kg()[i],
                    root_mass,
                )?,
                gap_projection: scaled(source.opening_projection_per_sqrt_kg()[i], root_mass)?,
                force_projection: [
                    scaled(force[0][i], root_mass)?,
                    scaled(force[1][i], root_mass)?,
                ],
            });
        }
        let viscosity = source
            .approximation()
            .modal_viscosity_per_s
            .iter()
            .map(|c| scaled(*c, mass))
            .collect::<Result<Vec<_>, _>>()?;
        let spec = ModalApertureSpec {
            rest_opening_m: source.rest_opening_m(),
            width_m: source.width_m(),
            modes,
            density_kg_m3: gas.density,
            impedance_pa_s_m3: options.impedance_pa_s_m3,
            time_step_s: options.time_step_s,
            max_steps: options.max_steps,
            limits: options.limits,
        };
        let mut valve = Self::new_with_propagator(spec, state, contact, BladePropagator::ExactZoh)?
            .with_viscosity(viscosity, gate)?;
        valve.material_solid = Some(MaterialModalSolid { source, root_mass });
        valve.material_gas = Some(gas);
        valve.validate_state(&valve.state)?;
        poll(gate)?;
        Ok(valve)
    }

    /// Original complete solid source, modal shapes and physical coordinate map.
    pub fn material_solid(&self) -> Option<&MaterialModalSolid> {
        self.material_solid.as_ref()
    }

    /// Accepted nodal displacement and velocity from this original solid image.
    /// Uses the same signs and scaling as pressure and the two force ports.
    pub fn solid_nodal_motion(
        &self,
        gate: &CancelGate,
    ) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), AcousticRealizeError> {
        let displacement: Vec<_> = self.state.iter().map(|s| s.displacement_m).collect();
        let velocity: Vec<_> = self.state.iter().map(|s| s.velocity_m_s).collect();
        if let Some(solid) = &self.hereditary_solid {
            return Ok((
                solid.expand_nodal(&displacement, gate)?,
                solid.expand_nodal(&velocity, gate)?,
            ));
        }
        let solid = self.material_solid.as_ref()
            .ok_or_else(|| invalid("modal aperture has no original solid image"))?;
        Ok((
            solid.expand_nodal(&displacement, gate)?,
            solid.expand_nodal(&velocity, gate)?,
        ))
    }
}
