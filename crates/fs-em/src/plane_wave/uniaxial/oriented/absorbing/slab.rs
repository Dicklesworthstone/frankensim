//! Finite absorbing uniaxial films with coherent internal reflection.
//!
//! Both crystal eigenpolarizations retain their original complex spectral
//! indices. Tangential E/H continuity is solved simultaneously at both faces;
//! dielectric heating uses the total field, including counterpropagating
//! interference. Modes are anchored at their incoming face so attenuation
//! never requires an exponentially growing transfer matrix.
//! Either exterior port can supply a coherent incident wave. The two drives
//! share frequency and transverse wavevector; their relative phase is authored
//! at their respective crystal faces, including any external optical path.
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
use fs_math::{c64::C64, dd::Dd, det};

mod absorption;
mod kernel;
pub use absorption::AbsorbingSlabVolumeAbsorption;
mod illumination;
pub use illumination::{AbsorbingSlabCoherentIllumination, AbsorbingSlabIllumination};
mod polarization;
pub use polarization::PreparedAbsorbingSlabPolarization;
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
    /// Incoming exit-side amplitudes anchored at z=thickness_m, with kz<0.
    pub right_incident_electric_v_m: Vector,
    pub right_incident_magnetic_a_m: Vector,
    pub reflected_electric_v_m: Vector,
    pub reflected_magnetic_a_m: Vector,
    /// Exit-face amplitudes, with no vacuum-phase removal.
    pub transmitted_electric_v_m: Vector,
    pub transmitted_magnetic_a_m: Vector,
    /// Ordinary then extraordinary, anchored at z=0.
    pub forward: [AbsorbingSlabWave; 2],
    /// Ordinary then extraordinary, anchored at z=thickness_m.
    pub backward: [AbsorbingSlabWave; 2],
    /// Positive inward power at the entrance port; zero for right-only driving.
    pub incident_normal_power_w_m2: f64,
    /// Positive inward power at the exit port; zero for entrance-only driving.
    pub right_incident_normal_power_w_m2: f64,
    /// Outgoing entrance-port power divided by the sum of both incident powers.
    pub reflected_power_fraction: f64,
    /// Outgoing exit-port power divided by the sum of both incident powers.
    pub transmitted_power_fraction: f64,
    /// Net power dissipated in the crystal, divided by the sum of incident powers.
    /// Evaluated from exterior Poynting balance.
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

    /// Sum of the positive incoming powers at the two exterior ports [W/m²].
    pub fn total_incident_normal_power_w_m2(&self) -> f64 {
        (Dd::from_f64(self.incident_normal_power_w_m2)
            + Dd::from_f64(self.right_incident_normal_power_w_m2))
        .to_f64()
    }

    /// Coherent peak-phasor fields, Poynting flux and nonnegative local tensor
    /// work at a physical point. Repeated observations reuse the solved modes.
    pub fn field_at(&self, region: SlabRegion, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        admit_point(region, point, self.thickness_m)?;
        let mut sum = FieldSum::at_point(point).map_err(ReflectionError::Invalid)?;
        self.add_waves(region, &mut sum)?;
        self.with_heating(region, sum.finish().map_err(ReflectionError::Invalid)?)
    }

    /// Incoming exit-side field for z>=thickness_m, with negative normal flux.
    /// A dark exit input returns a zero field, including an evanescent exit.
    pub fn exit_incident_field_at(&self, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        self.exit_field_at(point, false)
    }

    /// Coherent sum of the incoming and outgoing exit-side waves. The existing
    /// `field_at(SlabRegion::Transmitted, ...)` continues to return outgoing
    /// light only. This total field matches the crystal at the exit face.
    pub fn total_exit_field_at(&self, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        self.exit_field_at(point, true)
    }

    fn exit_field_at(
        &self,
        point: [f64; 3],
        include_outgoing: bool,
    ) -> Result<AbsorbingCrystalField> {
        admit_point(SlabRegion::Transmitted, point, self.thickness_m)?;
        let mut sum = FieldSum::at_point(point).map_err(ReflectionError::Invalid)?;
        if include_outgoing {
            self.add_waves(SlabRegion::Transmitted, &mut sum)?;
        }
        // An undriven evanescent exit has no incoming wave. Do not evaluate a
        // growing evanescent propagation factor merely to multiply it by zero.
        if self
            .right_incident_electric_v_m
            .iter()
            .any(|&v| v != C64::ZERO)
        {
            sum.add_wave(
                &self.right_incident_electric_v_m,
                &self.right_incident_magnetic_a_m,
                WaveGeometry {
                    tangential_wavenumber_rad_m: self.tangential_wavenumber_rad_m,
                    incidence_plane_axis: self.plane,
                    normal_wavenumber_rad_m: -self.exit_normal_wavenumber_rad_m,
                    anchor_z_m: self.thickness_m,
                },
            )
            .map_err(ReflectionError::Invalid)?;
        }
        self.with_heating(
            SlabRegion::Transmitted,
            sum.finish().map_err(ReflectionError::Invalid)?,
        )
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
    pub fn exit_incident_field_at(&self, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        self.response.exit_incident_field_at(point)
    }
    pub fn total_exit_field_at(&self, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        self.response.total_exit_field_at(point)
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

    /// Coherent illumination from both faces, sharing the transverse wavevector
    /// set by the entrance-side incidence cosine. Both arrays are global
    /// Cartesian Ex/Ey peak phasors: the entrance input is anchored at z=0 and
    /// the exit input at z=thickness_m. Either input may be dark; both dark
    /// inputs refuse because power fractions would be undefined. The entrance
    /// incidence must remain a propagating channel even when its input is dark;
    /// a nonzero exit input also requires a propagating exit channel. Thus this
    /// two-port drive uses the shared propagating transverse-wavevector range.
    pub fn prepare_absorbing_field_from_both_sides(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        electric_v_m: [C64; 2],
        right_electric_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'_, 's, 'm>> {
        self.prepare_absorbing_from_both_sides(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            [1.0, 0.0],
            electric_v_m,
            right_electric_v_m,
        )
    }

    /// Two-face coherent drive with a signed conserved transverse phase index.
    /// Port amplitudes and dark-input admission match
    /// [`Self::prepare_absorbing_field_from_both_sides`].
    pub fn prepare_absorbing_field_from_both_sides_at_tangential_index(
        &self,
        omega_rad_s: f64,
        index: f64,
        electric_v_m: [C64; 2],
        right_electric_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'_, 's, 'm>> {
        self.prepare_absorbing_from_both_sides(
            omega_rad_s,
            Incidence::TangentialIndex(index),
            [1.0, 0.0],
            electric_v_m,
            right_electric_v_m,
        )
    }

    /// Two-face coherent drive in an arbitrary transverse incidence direction.
    /// Both Cartesian Jones inputs rotate with the same original crystal axis.
    pub fn prepare_absorbing_field_from_both_sides_at_tangential_index_in_plane(
        &self,
        omega_rad_s: f64,
        index: f64,
        plane_axis: [f64; 2],
        electric_v_m: [C64; 2],
        right_electric_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'_, 's, 'm>> {
        self.prepare_absorbing_from_both_sides(
            omega_rad_s,
            Incidence::TangentialIndex(index),
            plane_axis,
            electric_v_m,
            right_electric_v_m,
        )
    }

    fn prepare_absorbing(
        &self,
        omega: f64,
        incidence: Incidence,
        plane_axis: [f64; 2],
        electric: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'_, 's, 'm>> {
        self.prepare_absorbing_from_both_sides(
            omega,
            incidence,
            plane_axis,
            electric,
            [C64::ZERO; 2],
        )
    }

    fn prepare_absorbing_from_both_sides(
        &self,
        omega: f64,
        incidence: Incidence,
        plane_axis: [f64; 2],
        electric: [C64; 2],
        right_electric: [C64; 2],
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
        let local = |e: [C64; 2]| {
            [
                e[0].scale(plane[0]) + e[1].scale(plane[1]),
                e[1].scale(plane[0]) - e[0].scale(plane[1]),
            ]
        };
        let mut response = kernel::solve_from_both_sides(
            crystal::Input {
                omega,
                ni,
                indices: source.crystal.complex_principal_indices(omega)?,
                mui: source.incident_mu_h_m,
                mut_: source.crystal.permeability_h_m(),
                axis: local_axis,
                incidence,
                electric: local(electric),
                loss_f_m: source
                    .crystal
                    .principal_permittivities_f_m(omega)?
                    .map(|e| -e.im),
            },
            exit_index,
            exit_mu,
            self.thickness_m(),
            local(right_electric),
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
