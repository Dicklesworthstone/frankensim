//! Coherent absorbing stacks of the original independently oriented crystals.
//!
//! The transparent stack's source owner and authored placements are retained.
//! Every finite layer keeps its own complex principal indices, scalar real
//! permeability and optic axis. A bounded reflection sweep matches tangential
//! Maxwell fields without growing transfer factors. Preparation and retained
//! storage are O(layers); a point or layer-volume query has fixed work.
//!
//! Exterior half-spaces must be lossless and isotropic. This frequency-domain
//! path admits reciprocal passive uniaxial layers, not gain, biaxial/Jordan
//! modes, causal pulse reservoirs or an implicit thermal feedback update.

pub use super::super::stack::StackRegion;
use super::super::{
    Incidence, OrientedCrystalField,
    kernel::{FieldSum, WaveGeometry},
    positive,
    slab::kernel::finite,
    spectrum::rotate,
    stack::OrientedUniaxialStack,
};
use super::{
    AbsorbingCrystalField,
    kernel::{heating, slab_volume_m3},
    slab::{
        AbsorbingSlabWave,
        absorption::{LayerWaves, integrated_layer_heat},
    },
};
use crate::plane_wave::{angular_spectrum::AngularSpectrumError, uniaxial::kernel::unit_axis};
use fs_math::{c64::C64, dd::Dd, det};

mod kernel;
mod polarization;
pub use polarization::PreparedAbsorbingStackPolarization;

const C0: f64 = 299_792_458.0;
type Result<T> = core::result::Result<T, AngularSpectrumError>;
type Vector = [C64; 3];

/// Four physical eigenwaves with local face anchors and original loss axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingStackLayerField {
    pub entrance_z_m: f64,
    pub exit_z_m: f64,
    pub thickness_m: f64,
    /// Ordinary then extraordinary, anchored at the layer's entrance.
    pub forward: [AbsorbingSlabWave; 2],
    /// Ordinary then extraordinary, anchored at the layer's exit.
    pub backward: [AbsorbingSlabWave; 2],
    omega: f64,
    axis: [f64; 3],
    loss_f_m: [f64; 2],
}

impl AbsorbingStackLayerField {
    pub const fn optic_axis(&self) -> [f64; 3] {
        self.axis
    }
    pub const fn principal_loss_permittivities_f_m(&self) -> [f64; 2] {
        self.loss_f_m
    }

    /// Positive tensor work in a prism of this layer. Depths are LOCAL to its
    /// entrance, so a thin layer remains representable deep in a large stack.
    /// All ordinary/extraordinary and forward/backward cross terms remain.
    /// No flux subtraction, quadrature, allocation or source query is used.
    pub fn absorbed_power_w(
        &self,
        area_m2: f64,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<f64> {
        integrated_layer_heat(
            LayerWaves {
                thickness_m: self.thickness_m,
                omega: self.omega,
                axis: self.axis,
                loss_f_m: self.loss_f_m,
                forward: &self.forward,
                backward: &self.backward,
            },
            area_m2,
            start_depth_m,
            end_depth_m,
        )
        .map_err(bad)
    }

    fn field_sum(&self, point: [f64; 3]) -> Result<FieldSum> {
        if !(self.entrance_z_m..=self.exit_z_m).contains(&point[2]) {
            return Err(bad("absorbing stack point is outside its selected layer"));
        }
        let depth = if point[2] == self.exit_z_m {
            Dd::from_f64(self.thickness_m)
        } else {
            Dd::from_f64(point[2]) - Dd::from_f64(self.entrance_z_m)
        };
        FieldSum::at_local_point([point[0], point[1]], depth).map_err(bad)
    }

    fn add_waves(&self, sum: &mut FieldSum, tangential: f64, plane: [f64; 2]) -> Result<()> {
        for (waves, anchor) in [(&self.forward, 0.0), (&self.backward, self.thickness_m)] {
            for wave in waves {
                sum.add_wave(
                    &wave.electric_v_m,
                    &wave.magnetic_a_m,
                    WaveGeometry {
                        tangential_wavenumber_rad_m: tangential,
                        incidence_plane_axis: plane,
                        normal_wavenumber_rad_m: wave.normal_wavenumber_rad_m,
                        anchor_z_m: anchor,
                    },
                )
                .map_err(bad)?;
            }
        }
        Ok(())
    }

    fn face(&self, end: bool) -> Result<OrientedCrystalField> {
        let z = if end {
            self.exit_z_m
        } else {
            self.entrance_z_m
        };
        let mut sum = self.field_sum([0.0, 0.0, z])?;
        self.add_waves(&mut sum, 0.0, [1.0, 0.0])?;
        sum.finish().map_err(bad)
    }
}

/// Complete coherent snapshot, including the physical field in every layer.
#[derive(Clone, Debug, PartialEq)]
pub struct AbsorbingStackResponse {
    pub thickness_m: f64,
    /// Signed magnitude along `incidence_plane_axis()`.
    pub tangential_wavenumber_rad_m: f64,
    pub incident_normal_wavenumber_rad_m: f64,
    pub exit_normal_wavenumber_rad_m: C64,
    pub exit_wave_kind: crate::plane_wave::uniaxial::UniaxialWaveKind,
    pub incident_electric_v_m: Vector,
    pub incident_magnetic_a_m: Vector,
    pub reflected_electric_v_m: Vector,
    pub reflected_magnetic_a_m: Vector,
    pub transmitted_electric_v_m: Vector,
    pub transmitted_magnetic_a_m: Vector,
    pub layers: Vec<AbsorbingStackLayerField>,
    pub incident_normal_power_w_m2: f64,
    pub reflected_power_fraction: f64,
    pub transmitted_power_fraction: f64,
    /// Exterior power balance; weak absorption may round to zero here.
    /// Use direct layer tensor work to retain a sub-resolution loss channel.
    pub absorbed_power_fraction: f64,
    omega: f64,
    plane: [f64; 2],
}

impl AbsorbingStackResponse {
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn incidence_plane_axis(&self) -> [f64; 2] {
        self.plane
    }

    /// Global Cartesian peak E/H, cycle-mean Poynting and local dielectric heat.
    /// Adjacent layers both admit their shared face and retain distinct normal
    /// fields. Repeated queries never solve or reevaluate a material law.
    pub fn field_at(&self, region: StackRegion, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        let field = self.plain_field_at(region, point)?;
        let absorbed_power_density_w_m3 = if let StackRegion::Layer(index) = region {
            let layer = &self.layers[index];
            heating(field.electric_v_m, layer.axis, layer.loss_f_m, layer.omega).map_err(bad)?
        } else {
            0.0
        };
        Ok(AbsorbingCrystalField {
            field,
            absorbed_power_density_w_m3,
        })
    }

    fn plain_field_at(&self, region: StackRegion, point: [f64; 3]) -> Result<OrientedCrystalField> {
        if point.iter().any(|v| !v.is_finite()) {
            return Err(bad("finite absorbing stack field point required"));
        }
        if let StackRegion::Layer(index) = region {
            let layer = self
                .layers
                .get(index)
                .ok_or_else(|| bad("unknown absorbing stack layer"))?;
            let mut sum = layer.field_sum(point)?;
            layer.add_waves(&mut sum, self.tangential_wavenumber_rad_m, self.plane)?;
            return sum.finish().map_err(bad);
        }
        if if region == StackRegion::Transmitted {
            point[2] < self.thickness_m
        } else {
            point[2] > 0.0
        } {
            return Err(bad(
                "absorbing stack point is outside its selected exterior",
            ));
        }
        let mut sum = FieldSum::at_point(point).map_err(bad)?;
        let mut add = |e: &Vector, h: &Vector, normal, anchor| {
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
            .map_err(bad)
        };
        if matches!(
            region,
            StackRegion::Incident | StackRegion::TotalIncidentSide
        ) {
            add(
                &self.incident_electric_v_m,
                &self.incident_magnetic_a_m,
                C64::from_re(self.incident_normal_wavenumber_rad_m),
                0.0,
            )?;
        }
        if matches!(
            region,
            StackRegion::Reflected | StackRegion::TotalIncidentSide
        ) {
            add(
                &self.reflected_electric_v_m,
                &self.reflected_magnetic_a_m,
                C64::from_re(-self.incident_normal_wavenumber_rad_m),
                0.0,
            )?;
        }
        if region == StackRegion::Transmitted {
            add(
                &self.transmitted_electric_v_m,
                &self.transmitted_magnetic_a_m,
                self.exit_normal_wavenumber_rad_m,
                self.thickness_m,
            )?;
        }
        sum.finish().map_err(bad)
    }

    // Mode matching happens in the local incidence plane. Source axes already
    // remain global; only solved vectors and the transverse direction rotate.
    fn rotate(&mut self, plane: [f64; 2]) -> Result<()> {
        for vector in [
            &mut self.incident_electric_v_m,
            &mut self.incident_magnetic_a_m,
            &mut self.reflected_electric_v_m,
            &mut self.reflected_magnetic_a_m,
            &mut self.transmitted_electric_v_m,
            &mut self.transmitted_magnetic_a_m,
        ] {
            *vector = rotate(*vector, plane);
            if vector.iter().any(|&v| !finite(v)) {
                return Err(bad("absorbing stack rotated field is unrepresentable"));
            }
        }
        for layer in &mut self.layers {
            for wave in layer.forward.iter_mut().chain(&mut layer.backward) {
                wave.electric_v_m = rotate(wave.electric_v_m, plane);
                wave.magnetic_a_m = rotate(wave.magnetic_a_m, plane);
                if wave
                    .electric_v_m
                    .iter()
                    .chain(&wave.magnetic_a_m)
                    .any(|&v| !finite(v))
                {
                    return Err(bad("absorbing stack rotated field is unrepresentable"));
                }
            }
        }
        self.plane = plane;
        Ok(())
    }
}

/// The solved field borrows the exact original ordered material sources.
pub struct PreparedAbsorbingStackField<'f, 's, 'm> {
    source: &'f OrientedUniaxialStack<'s, 'm>,
    response: AbsorbingStackResponse,
}

impl<'f, 's, 'm> PreparedAbsorbingStackField<'f, 's, 'm> {
    pub const fn source(&self) -> &'f OrientedUniaxialStack<'s, 'm> {
        self.source
    }
    pub const fn response(&self) -> &AbsorbingStackResponse {
        &self.response
    }
    pub fn field_at(&self, region: StackRegion, point: [f64; 3]) -> Result<AbsorbingCrystalField> {
        self.response.field_at(region, point)
    }

    /// Exact work in a finite layer prism. The returned observation retains
    /// this prepared illumination and the original layer/material identity.
    pub fn layer_absorption(
        &self,
        layer_index: usize,
        area_m2: f64,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<AbsorbingStackVolumeAbsorption<'_, 'f, 's, 'm>> {
        let layer = self
            .response
            .layers
            .get(layer_index)
            .ok_or_else(|| bad("unknown absorbing stack layer"))?;
        let absorbed_power_w = layer.absorbed_power_w(area_m2, start_depth_m, end_depth_m)?;
        let volume_m3 = slab_volume_m3(area_m2, start_depth_m, end_depth_m).map_err(bad)?;
        Ok(AbsorbingStackVolumeAbsorption {
            prepared_field: self,
            layer_index,
            area_m2,
            start_depth_m,
            end_depth_m,
            volume_m3,
            absorbed_power_w,
        })
    }
}

pub struct AbsorbingStackVolumeAbsorption<'w, 'f, 's, 'm> {
    prepared_field: &'w PreparedAbsorbingStackField<'f, 's, 'm>,
    layer_index: usize,
    area_m2: f64,
    start_depth_m: f64,
    end_depth_m: f64,
    volume_m3: f64,
    absorbed_power_w: f64,
}

impl<'w, 'f, 's, 'm> AbsorbingStackVolumeAbsorption<'w, 'f, 's, 'm> {
    pub const fn prepared_field(&self) -> &'w PreparedAbsorbingStackField<'f, 's, 'm> {
        self.prepared_field
    }
    pub const fn source(&self) -> &'f OrientedUniaxialStack<'s, 'm> {
        self.prepared_field.source
    }
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn start_depth_m(&self) -> f64 {
        self.start_depth_m
    }
    pub const fn end_depth_m(&self) -> f64 {
        self.end_depth_m
    }
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    pub const fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }
}

impl<'s, 'm> OrientedUniaxialStack<'s, 'm> {
    /// One-sided coherent global Ex/Ey illumination in the x-z plane.
    pub fn prepare_absorbing_field(
        &self,
        omega: f64,
        cosine: f64,
        electric: [C64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAbsorbingStackField<'_, 's, 'm>> {
        self.prepare_absorbing_impl(
            omega,
            Incidence::Cosine(cosine),
            [1.0, 0.0],
            electric,
            &mut cancelled,
        )
    }

    /// Arbitrary global transverse direction with the original positive
    /// incidence cosine preserved at grazing incidence.
    pub fn prepare_absorbing_field_in_plane(
        &self,
        omega: f64,
        cosine: f64,
        plane: [f64; 2],
        electric: [C64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAbsorbingStackField<'_, 's, 'm>> {
        self.prepare_absorbing_impl(
            omega,
            Incidence::Cosine(cosine),
            plane,
            electric,
            &mut cancelled,
        )
    }

    /// Signed conserved transverse index, retaining near-critical input
    /// without converting it through a rounded incidence angle.
    pub fn prepare_absorbing_field_at_tangential_index_in_plane(
        &self,
        omega: f64,
        index: f64,
        plane: [f64; 2],
        electric: [C64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAbsorbingStackField<'_, 's, 'm>> {
        self.prepare_absorbing_impl(
            omega,
            Incidence::TangentialIndex(index),
            plane,
            electric,
            &mut cancelled,
        )
    }

    fn prepare_absorbing_impl(
        &self,
        omega: f64,
        incidence: Incidence,
        plane: [f64; 2],
        electric: [C64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedAbsorbingStackField<'_, 's, 'm>> {
        poll(cancelled)?;
        if kernel::workspace_bytes(self.layers().len(), 1)?
            > self.preparation_limits().maximum_workspace_bytes
        {
            return Err(bad(
                "absorbing stack exceeds its original preparation budget",
            ));
        }
        let plane = unit_axis(plane).map_err(bad)?;
        let data = self.resolve_absorbing_frequency(omega, cancelled)?;
        let response = kernel::solve(&data, incidence, electric, plane, cancelled)?;
        Ok(PreparedAbsorbingStackField {
            source: self,
            response,
        })
    }

    fn resolve_absorbing_frequency(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Frequency> {
        poll(cancelled)?;
        let omega = positive(omega)?;
        let [mui, mut_] = self.exterior_permeabilities_h_m();
        let ni = positive(
            (C0 * det::sqrt(mui))
                * det::sqrt(self.incident_response().lossless_permittivity(omega)?),
        )?;
        let nt = positive(
            (C0 * det::sqrt(mut_)) * det::sqrt(self.exit_response().lossless_permittivity(omega)?),
        )?;
        let mut layers = Vec::new();
        layers
            .try_reserve_exact(self.layers().len())
            .map_err(|_| bad("absorbing stack source allocation failed"))?;
        let mut entrance = 0.0;
        for source in self.layers() {
            poll(cancelled)?;
            let crystal = source.crystal_response();
            let thickness = source.thickness_m();
            let exit = entrance + thickness;
            layers.push(ResolvedLayer {
                indices: crystal.complex_principal_indices(omega)?,
                mu: crystal.permeability_h_m(),
                axis: source.optic_axis(),
                loss_f_m: crystal.principal_permittivities_f_m(omega)?.map(|v| -v.im),
                thickness,
                entrance,
                exit,
            });
            entrance = exit;
        }
        Ok(Frequency {
            omega,
            ni,
            nt,
            mui,
            mut_,
            layers,
        })
    }
}

#[derive(Clone, Copy)]
struct ResolvedLayer {
    indices: [C64; 2],
    mu: f64,
    axis: [f64; 3],
    loss_f_m: [f64; 2],
    thickness: f64,
    entrance: f64,
    exit: f64,
}

struct Frequency {
    omega: f64,
    ni: f64,
    nt: f64,
    mui: f64,
    mut_: f64,
    layers: Vec<ResolvedLayer>,
}

fn bad(message: &'static str) -> AngularSpectrumError {
    AngularSpectrumError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(AngularSpectrumError::Cancelled)
    } else {
        Ok(())
    }
}
