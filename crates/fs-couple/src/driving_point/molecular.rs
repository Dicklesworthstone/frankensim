//! Molecular gas storage and wall diffusion in the ordinary characteristic port.
//!
//! The prepared `fs-duct` owner supplies every positive-frequency reflection;
//! the existing FIR ladder, passivity check and reed runtime remain unchanged.
//! Reservoir times therefore change pressure feedback, not an output filter.
//! No second sound-speed model, extra ISO absorption or copied wall loss is
//! added. The source gas stays at local equilibrium for the inlet pressure/flow
//! reference, while the original frozen sound speed sets the wavefront delay.
//!
//! This retains the owner's frozen plane-wave model and explicit
//! gas shear/thermal diffusion choice. Compact finite tone-hole chimneys use
//! those SAME reservoirs and wall law. Their independent flanged mouth loads
//! and the terminal radiation use the explicitly bound exterior gas. An
//! optional source-resolved material wall supplies its original local membrane
//! impedance to that same gas operator, including every finite chimney. The
//! feedback and outlet observation therefore share the wall's storage and
//! hereditary loss. Nonlocal shell motion, external wall radiation, nonlinear
//! vortices and mutual exterior radiation remain outside this consumer.

use super::regional::material_wall::LocalMaterialWall;
use super::{
    DrivingPointError, LineDc, LineLadder, LineRealization, bulk_delay_pole_state,
    duct_axial_length, realize_characteristic_response,
};
use fs_duct::compose::{ComposeOptions, compose_chain_reflection};
use fs_duct::molecular::{MolecularDuct, MolecularHoleSite, MolecularWallResolution};
use fs_duct::{ConeSlices, Duct, Segment, Termination};
use fs_material::gas::GasState;
use fs_material::relaxation::RelaxationState;
use fs_math::c64::C64;
use fs_vfit::discretize::DelayedFilter;
use std::sync::Arc;

mod outlet;
pub use outlet::MAX_MOLECULAR_OUTLET_BINS;
mod field;

/// Original molecular gas, geometry and exterior load carried into time playback.
#[derive(Debug, Clone)]
pub struct MolecularDrivingPoint {
    geometry: Duct,
    material: MolecularDuct,
    exterior: GasState,
    composition: ComposeOptions,
    one_way_wavefront_s: f64,
    open_hole_sites: Vec<MolecularHoleSite>,
    open_hole_wavefront_s: Vec<f64>,
    material_wall: Option<Arc<LocalMaterialWall>>,
}

impl MolecularDrivingPoint {
    /// Bind the existing molecular constitutive owner to an explicit exterior.
    /// `wall` and `cones` have exactly `MolecularDuct::new`'s meaning; all
    /// original diffusion/source/geometry and composition limits remain binding.
    /// No material or gas coefficient is inferred from an instrument name.
    pub fn new(
        geometry: &Duct,
        state: RelaxationState,
        wall: Option<MolecularWallResolution>,
        cones: ConeSlices,
        composition: ComposeOptions,
        exterior: GasState,
    ) -> Result<Self, DrivingPointError> {
        Self::build(geometry, state, wall, None, cones, composition, exterior)
    }

    /// Bind molecular gas and a source-resolved deformable wall to the same
    /// pressure feedback and physical mouth-flow observations.
    ///
    /// The existing wall owner retains the material identity, thickness,
    /// axial constraint, anisotropy and supported elastic/viscous/Maxwell law.
    /// Each actual bore, profile and chimney radius must satisfy that owner's
    /// local membrane admission. This is one material over the complete duct;
    /// no regional shell coupling or radiation from the solid wall is inferred.
    ///
    /// `wall` is the separate gas shear/thermal resolution. Its coupled
    /// molecular heat law remains unchanged: mechanical wall motion adds only
    /// the original pressure/normal-velocity shunt. The material's positive
    /// surface mass leaves the gas's frozen high-frequency wavefront intact.
    /// Every evaluated frequency, including the near-DC realization probe,
    /// must lie in both source domains. Missing low-frequency material data
    /// refuses instead of silently substituting a rigid wall.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_material_wall(
        geometry: &Duct,
        state: RelaxationState,
        wall: Option<MolecularWallResolution>,
        material_wall: LocalMaterialWall,
        cones: ConeSlices,
        composition: ComposeOptions,
        exterior: GasState,
    ) -> Result<Self, DrivingPointError> {
        Self::build(
            geometry,
            state,
            wall,
            Some(Arc::new(material_wall)),
            cones,
            composition,
            exterior,
        )
    }

    fn build(
        geometry: &Duct,
        state: RelaxationState,
        wall: Option<MolecularWallResolution>,
        material_wall: Option<Arc<LocalMaterialWall>>,
        cones: ConeSlices,
        composition: ComposeOptions,
        exterior: GasState,
    ) -> Result<Self, DrivingPointError> {
        Self::admit_geometry_exterior(geometry, state.gas(), &exterior)?;
        let material = match &material_wall {
            Some(material) => MolecularDuct::new_with_compliant_wall(
                geometry,
                state,
                wall,
                material.clone(),
                cones,
                composition.maximum_blocks,
            )?,
            None => MolecularDuct::new(geometry, state, wall, cones, composition.maximum_blocks)?,
        };
        Self::from_prepared(geometry, material, material_wall, composition, exterior)
    }

    fn admit_geometry_exterior(
        geometry: &Duct,
        inlet: &GasState,
        exterior: &GasState,
    ) -> Result<(), DrivingPointError> {
        exterior
            .validate_coherence()
            .map_err(DrivingPointError::Material)?;
        if exterior.pressure != inlet.pressure {
            return Err(DrivingPointError::Invalid {
                what: "exterior and molecular duct need a common static pressure; mean flow is not modeled",
            });
        }
        if matches!(geometry.segments.first(), Some(Segment::ToneHole { .. }))
            || matches!(geometry.segments.last(), Some(Segment::ToneHole { .. }))
        {
            return Err(DrivingPointError::Invalid {
                what: "molecular characteristic playback requires axial inlet and terminal segments",
            });
        }
        Ok(())
    }

    fn from_prepared(
        geometry: &Duct,
        material: MolecularDuct,
        material_wall: Option<Arc<LocalMaterialWall>>,
        composition: ComposeOptions,
        exterior: GasState,
    ) -> Result<Self, DrivingPointError> {
        let open_hole_sites = material.open_hole_sites();
        let (one_way_wavefront_s, open_hole_wavefront_s) =
            if let Some(times) = material.field_travel_times_s() {
                (times.terminal_s, times.open_mouths_s)
            } else {
                let speed = material.state().frozen_sound_speed_m_s();
                (
                    duct_axial_length(geometry) / speed,
                    open_hole_sites
                        .iter()
                        .map(|site| (site.axial_position_m + site.open_chimney_length_m) / speed)
                        .collect(),
                )
            };
        if !one_way_wavefront_s.is_finite()
            || one_way_wavefront_s <= 0.0
            || open_hole_wavefront_s.len() != open_hole_sites.len()
            || open_hole_wavefront_s
                .iter()
                .any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(DrivingPointError::Invalid {
                what: "molecular line and every open mouth need finite positive wavefront times",
            });
        }
        Ok(Self {
            geometry: geometry.clone(),
            material,
            exterior,
            composition,
            one_way_wavefront_s,
            open_hole_sites,
            open_hole_wavefront_s,
            material_wall,
        })
    }

    /// Original geometry, without a replacement equivalent cylinder.
    #[must_use]
    pub const fn geometry(&self) -> &Duct {
        &self.geometry
    }
    /// Original material reservoirs and prepared wall-diffusion operators.
    #[must_use]
    pub const fn material(&self) -> &MolecularDuct {
        &self.material
    }
    /// Original wall source used by feedback, terminal flow and every chimney.
    #[must_use]
    pub fn material_wall(&self) -> Option<&LocalMaterialWall> {
        self.material_wall.as_deref()
    }
    /// Equilibrium inlet state shared by the pressure/flow and exciter owners.
    #[must_use]
    pub fn inlet_gas(&self) -> &GasState {
        self.material.state().gas()
    }
    /// Original surrounding gas used by the terminal radiation law.
    #[must_use]
    pub const fn exterior_gas(&self) -> &GasState {
        &self.exterior
    }
    /// One-way axial inlet-to-terminal wavefront time with reservoirs frozen [s].
    /// Side branches can return pressure before the full terminal round trip.
    #[must_use]
    pub const fn one_way_wavefront_s(&self) -> f64 {
        self.one_way_wavefront_s
    }
    /// Actual inlet radius [m], retained from the admitted geometry.
    #[must_use]
    pub fn inlet_radius_m(&self) -> f64 {
        self.geometry.segments[0].inlet_radius()
    }

    /// Open or partially vented mouths in authored order. Each uses its full
    /// radius for the flanged exterior load; its original admittance mixture
    /// applies `open_fraction` exactly once to the solved mouth volume flow.
    #[must_use]
    pub fn open_hole_sites(&self) -> &[MolecularHoleSite] {
        &self.open_hole_sites
    }

    /// Inlet-to-mouth gas wavefront times, in `open_hole_sites()` order [s].
    /// Each finite chimney uses its local material; exterior delay is separate.
    #[must_use]
    pub fn open_hole_wavefront_times_s(&self) -> &[f64] {
        &self.open_hole_wavefront_s
    }

    fn hole_mouth_loads(&self, omega: f64) -> Result<Vec<C64>, DrivingPointError> {
        self.open_hole_sites
            .iter()
            .map(|site| {
                fs_duct::compact_radiation_impedance(
                    Termination::FlangedOpen,
                    &self.exterior,
                    site.hole_radius_m,
                    omega,
                )
                .map_err(DrivingPointError::from)
            })
            .collect()
    }

    fn terminal_load(
        &self,
        omega: f64,
        termination: Termination,
    ) -> Result<Option<C64>, DrivingPointError> {
        Ok(match termination {
            Termination::Closed => None,
            Termination::IdealOpen => Some(C64::ZERO),
            other => Some(fs_duct::compact_radiation_impedance(
                other,
                &self.exterior,
                self.geometry
                    .segments
                    .last()
                    .expect("admitted nonempty geometry")
                    .outlet_radius(),
                omega,
            )?),
        })
    }

    /// Actual source-bound reflection at a positive frequency, exp(-i omega t).
    /// Scattering composition avoids forming an infinite impedance at a
    /// lossless resonance. `zc` only selects the real pressure-wave reference.
    pub fn reflectance(
        &self,
        omega_rad_s: f64,
        termination: Termination,
        zc: f64,
    ) -> Result<C64, DrivingPointError> {
        let loads = self.hole_mouth_loads(omega_rad_s)?;
        let chain = self.material.prepared_chain_with_hole_loads_controlled(
            omega_rad_s,
            &loads,
            self.composition,
            &mut || false,
        )?;
        let load = self.terminal_load(omega_rad_s, termination)?;
        Ok(compose_chain_reflection(&chain.blocks, load, self.composition)?.referenced_to(zc)?)
    }

    /// Realize this constitutive response through the shared converged FIR owner.
    /// Every sampled bin, including Nyquist, must lie in the wall/source band.
    /// Open mouths also require the existing 1e-3 rad/s near-DC probe; a refusal
    /// propagates instead of replacing the low-frequency material response.
    /// The source's diffusion and molecular storage are counted only once.
    pub fn realize_line(
        &self,
        termination: Termination,
        sample_rate_hz: u32,
        zc: f64,
        ladder: &LineLadder,
    ) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
        let dc = if termination == Termination::Closed && self.open_hole_sites.is_empty() {
            1.0
        } else {
            self.reflectance(1e-3, termination, zc)?.re
        };
        let delay = 2.0 * self.one_way_wavefront_s * f64::from(sample_rate_hz);
        realize_characteristic_response(
            delay,
            LineDc::Explicit(dc),
            sample_rate_hz,
            zc,
            vec![bulk_delay_pole_state(delay)],
            ladder,
            |omega| self.reflectance(omega, termination, zc),
        )
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod material_wall_tests;
