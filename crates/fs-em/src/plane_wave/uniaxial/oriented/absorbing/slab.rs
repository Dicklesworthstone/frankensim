//! Finite absorbing uniaxial films with coherent internal reflection.
//!
//! Both crystal eigenpolarizations retain their original complex spectral
//! indices. Tangential E/H continuity is solved simultaneously at both faces;
//! dielectric heating uses the total field, including counterpropagating
//! interference. Modes are anchored at their incoming face so attenuation
//! never requires an exponentially growing transfer matrix.
//!
//! The source is the existing finite-slab owner. The surrounding half-spaces
//! must be lossless and isotropic; the crystal is homogeneous and uniaxial,
//! with real scalar permeability and a fixed arbitrary optic axis. Coalesced or
//! unresolved eigenmodes refuse. This is a monochromatic coherent solution,
//! not a causal time-domain constitutive model. See Berreman, JOSA 62, 502
//! (1972), <https://doi.org/10.1364/JOSA.62.000502> for tangential Maxwell matching.

use super::super::super::kernel::unit_axis;
use super::super::slab::{OrientedUniaxialSlab, SlabRegion};
use super::super::{
    Incidence, OrientedCrystalField, ReflectionError,
    kernel::{FieldSum, WaveGeometry, unit},
    positive,
};
use super::{AbsorbingCrystalField, kernel as crystal};
use fs_material::field_storage::linear_magnetic_permeability;
use fs_math::{c64::C64, det};

mod absorption;
mod kernel;
pub use absorption::AbsorbingSlabVolumeAbsorption;
mod illumination;
pub use illumination::{AbsorbingSlabCoherentIllumination, AbsorbingSlabIllumination};
const C0: f64 = 299_792_458.0;
type Result<T> = core::result::Result<T, ReflectionError>;
type Vector = [C64; 3];

/// A physical absorbing eigenwave, without a lossless bulk-ray classification.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingSlabWave {
    pub normal_wavenumber_rad_m: C64,
    pub electric_v_m: Vector,
    pub magnetic_a_m: Vector,
}

/// Coherent fields and independently evaluated exterior normal power.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingSlabResponse {
    pub thickness_m: f64,
    pub tangential_wavenumber_rad_m: f64,
    pub incident_normal_wavenumber_rad_m: f64,
    pub exit_normal_wavenumber_rad_m: C64,
    pub incident_electric_v_m: Vector,
    pub incident_magnetic_a_m: Vector,
    pub reflected_electric_v_m: Vector,
    pub reflected_magnetic_a_m: Vector,
    /// Exit-face amplitudes, with no vacuum-phase removal.
    pub transmitted_electric_v_m: Vector,
    pub transmitted_magnetic_a_m: Vector,
    /// Ordinary then extraordinary, anchored at z=0.
    pub forward: [AbsorbingSlabWave; 2],
    /// Ordinary then extraordinary, anchored at z=thickness_m.
    pub backward: [AbsorbingSlabWave; 2],
    pub incident_normal_power_w_m2: f64,
    pub reflected_power_fraction: f64,
    pub transmitted_power_fraction: f64,
    /// Net power dissipated in the crystal, from exterior Poynting balance.
    /// Like a difference of measured fluxes, this can round to zero for losses
    /// below floating-point resolution; local tensor heating remains available.
    pub absorbed_power_fraction: f64,
    omega: f64,
    axis: [f64; 3],
    loss_f_m: [f64; 2],
    plane: [f64; 2],
}

impl AbsorbingSlabResponse {
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn incidence_plane_axis(&self) -> [f64; 2] {
        self.plane
    }
    pub const fn principal_loss_permittivities_f_m(&self) -> [f64; 2] {
        self.loss_f_m
    }

    /// Coherent peak-phasor fields, Poynting flux and nonnegative local tensor
    /// work at a physical point. Repeated observations reuse the solved modes.
    pub fn field_at(&self, region: SlabRegion, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        admit_point(region, point, self.thickness_m)?;
        let mut sum = FieldSum::at_point(point).map_err(ReflectionError::Invalid)?;
        self.add_waves(region, &mut sum)?;
        self.with_heating(region, sum.finish().map_err(ReflectionError::Invalid)?)
    }

    // Single fields and coherent groups share the original anchored waves.
    // A group finalizes E/H only after all components have interfered.
    fn add_waves(&self, region: SlabRegion, sum: &mut FieldSum) -> Result<()> {
        let mut add = |e, h, normal, anchor| {
            sum.add_wave(
                e,
                h,
                WaveGeometry {
                    tangential_wavenumber_rad_m: self.tangential_wavenumber_rad_m,
                    incidence_plane_axis: self.plane,
                    normal_wavenumber_rad_m: normal,
                    anchor_z_m: anchor,
                },
            )
            .map_err(ReflectionError::Invalid)
        };
        if matches!(region, SlabRegion::Incident | SlabRegion::TotalIncidentSide) {
            add(
                &self.incident_electric_v_m,
                &self.incident_magnetic_a_m,
                C64::from_re(self.incident_normal_wavenumber_rad_m),
                0.0,
            )?;
        }
        if matches!(
            region,
            SlabRegion::Reflected | SlabRegion::TotalIncidentSide
        ) {
            add(
                &self.reflected_electric_v_m,
                &self.reflected_magnetic_a_m,
                C64::from_re(-self.incident_normal_wavenumber_rad_m),
                0.0,
            )?;
        }
        if region == SlabRegion::Crystal {
            for w in &self.forward {
                add(
                    &w.electric_v_m,
                    &w.magnetic_a_m,
                    w.normal_wavenumber_rad_m,
                    0.0,
                )?;
            }
            for w in &self.backward {
                add(
                    &w.electric_v_m,
                    &w.magnetic_a_m,
                    w.normal_wavenumber_rad_m,
                    self.thickness_m,
                )?;
            }
        }
        if region == SlabRegion::Transmitted {
            add(
                &self.transmitted_electric_v_m,
                &self.transmitted_magnetic_a_m,
                self.exit_normal_wavenumber_rad_m,
                self.thickness_m,
            )?;
        }
        Ok(())
    }

    fn with_heating(
        &self,
        region: SlabRegion,
        field: OrientedCrystalField,
    ) -> Result<AbsorbingCrystalField> {
        let absorbed_power_density_w_m3 = if region == SlabRegion::Crystal {
            crystal::heating(field.electric_v_m, self.axis, self.loss_f_m, self.omega)
                .map_err(ReflectionError::Invalid)?
        } else {
            0.0
        };
        Ok(AbsorbingCrystalField {
            field,
            absorbed_power_density_w_m3,
        })
    }
}

fn admit_point(region: SlabRegion, point: [f64; 3], thickness_m: f64) -> Result<()> {
    let admitted = match region {
        SlabRegion::Crystal => (0.0..=thickness_m).contains(&point[2]),
        SlabRegion::Transmitted => point[2] >= thickness_m,
        _ => point[2] <= 0.0,
    };
    if !admitted || point.iter().any(|v| !v.is_finite()) {
        return Err(ReflectionError::Invalid(
            "absorbing film point is outside its selected region",
        ));
    }
    Ok(())
}

/// Solved finite film retaining the exact entrance, crystal and exit sources.
pub struct PreparedAbsorbingSlabField<'f, 's, 'm> {
    source: &'f OrientedUniaxialSlab<'s, 'm>,
    response: AbsorbingSlabResponse,
}

impl<'f, 's, 'm> PreparedAbsorbingSlabField<'f, 's, 'm> {
    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
        self.source
    }
    pub const fn response(&self) -> &AbsorbingSlabResponse {
        &self.response
    }
    pub fn field_at(&self, region: SlabRegion, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        self.response.field_at(region, point)
    }
}

impl<'s, 'm> OrientedUniaxialSlab<'s, 'm> {
    /// Original global Ex/Ey peak amplitudes, with the incident wave in x-z.
    pub fn prepare_absorbing_field(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        electric_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'_, 's, 'm>> {
        self.prepare_absorbing(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            [1.0, 0.0],
            electric_v_m,
        )
    }

    /// Signed conserved transverse phase index; preserves near-critical input
    /// without converting it through a rounded incidence angle.
    pub fn prepare_absorbing_field_at_tangential_index(
        &self,
        omega_rad_s: f64,
        index: f64,
        electric_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'_, 's, 'm>> {
        self.prepare_absorbing(
            omega_rad_s,
            Incidence::TangentialIndex(index),
            [1.0, 0.0],
            electric_v_m,
        )
    }

    /// Arbitrary global transverse incidence direction and Cartesian Ex/Ey.
    /// The original optic axis, fields and conserved transverse wavevector
    /// rotate together; the crystal is never silently reoriented.
    pub fn prepare_absorbing_field_at_tangential_index_in_plane(
        &self,
        omega_rad_s: f64,
        index: f64,
        plane_axis: [f64; 2],
        electric_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'_, 's, 'm>> {
        self.prepare_absorbing(
            omega_rad_s,
            Incidence::TangentialIndex(index),
            plane_axis,
            electric_v_m,
        )
    }

    fn prepare_absorbing(
        &self,
        omega: f64,
        incidence: Incidence,
        plane_axis: [f64; 2],
        electric: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'_, 's, 'm>> {
        let omega = positive(omega)?;
        let source = self.entrance_interface();
        let ni = positive(
            (C0 * det::sqrt(source.incident_mu_h_m))
                * det::sqrt(source.incident.lossless_permittivity(omega)?),
        )?;
        let exit = self.exit_response();
        let exit_mu = linear_magnetic_permeability(exit.material())?;
        let exit_index =
            positive((C0 * det::sqrt(exit_mu)) * det::sqrt(exit.lossless_permittivity(omega)?))?;
        let axis = unit(source.optic_axis).map_err(ReflectionError::Invalid)?;
        let plane = unit_axis(plane_axis)
            .map_err(|_| ReflectionError::Invalid("finite nonzero incidence plane required"))?;
        let local_axis = [
            axis[0] * plane[0] + axis[1] * plane[1],
            axis[1] * plane[0] - axis[0] * plane[1],
            axis[2],
        ];
        let local_e = [
            electric[0].scale(plane[0]) + electric[1].scale(plane[1]),
            electric[1].scale(plane[0]) - electric[0].scale(plane[1]),
        ];
        let mut response = kernel::solve(
            crystal::Input {
                omega,
                ni,
                indices: source.crystal.complex_principal_indices(omega)?,
                mui: source.incident_mu_h_m,
                mut_: source.crystal.permeability_h_m(),
                axis: local_axis,
                incidence,
                electric: local_e,
                loss_f_m: source
                    .crystal
                    .principal_permittivities_f_m(omega)?
                    .map(|e| -e.im),
            },
            exit_index,
            exit_mu,
            self.thickness_m(),
        )
        .map_err(ReflectionError::Invalid)?;
        response
            .rotate(plane, axis)
            .map_err(ReflectionError::Invalid)?;
        Ok(PreparedAbsorbingSlabField {
            source: self,
            response,
        })
    }
}
