//! Source-owned shell motion in the resolved pressure/velocity/temperature solve.
//!
//! The original shell's signed GLOBAL P1 midsurface translations are represented
//! exactly by the fluid's hierarchical P2 trace: vertex coefficients are copied,
//! edge bubbles are zero. Rotations retain their original mechanical influence
//! through the source modes; they are not added to midsurface translation.
//! Native no-slip congruence supplies the transpose pressure AND shear reaction.
//! No compact volume spring, exterior radiation load, refit or new mode is added.
//!
//! This bounded adapter requires coincident conforming planar boundary facets.
//! The source's cavity selection and orientation are authoritative. It does not
//! weld vertices, extrapolate a skin offset, invent missing faces, or transfer
//! nonmatching meshes. Modes with no resolved boundary translation refuse under
//! the native owner; none are silently removed. The first law envelope is an
//! instantaneous Kelvin--Voigt substrate and its original linear interfaces.
//! Frequency-local coating/constant-loss slices and separate material memories
//! require another causal owner and refuse instead of losing their storage.

use super::{ThermoviscousInterior, ThermoviscousSourceError};
use crate::shear_shell::MaterialShellImage;
use fs_exec::CancelGate;
use fs_feec::acoustics::helmholtz::HelmholtzFe;
use fs_feec::acoustics::thermoviscous::moving_wall::{MAX_MOVING_WALL_MODES, MovingWall};
use fs_feec::acoustics::thermoviscous::{
    ThermoviscousBoundary, ThermoviscousError, ThermoviscousFe, ThermoviscousLimits,
    ThermoviscousMedium, ThermoviscousThermodynamics,
};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_material::{gas::GasState, relaxation::RelaxationState};
use fs_math::{c64::C64, dd::Dd};
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeError, ModalAcousticTimeModel,
};
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};

#[cfg(test)]
mod tests;

/// Explicit resources for the complete original source and its trace transfer.
#[derive(Clone, Copy, Debug)]
pub struct ShellWallControls {
    pub maximum_source_nodes: usize,
    pub maximum_source_triangles: usize,
    pub maximum_source_interfaces: usize,
    pub maximum_force_ports: usize,
    /// Complete original mode count, additionally capped by the native 64.
    pub maximum_modes: usize,
    /// Material-factor columns PLUS nonzero interface dashpots, at most 64.
    /// Excess channels refuse; no cross term is dropped or refactored to fit.
    pub maximum_damping_channels: usize,
    /// Peak additional numeric/mapping payload, excluding borrowed sources,
    /// the supplied geometry, allocator bookkeeping and native FE assembly.
    /// Admission includes all mode traces, channels and their dense Gram,
    /// vertex maps, selected facets, modal metadata and bounded scratch.
    pub maximum_setup_bytes: usize,
    /// Conservative deterministic scalar-work score for the complete mapping,
    /// source validation, volume witness and channel Gram. Native assembly has
    /// its independent `fluid` extent; no source truncation satisfies a budget.
    pub maximum_setup_work: u64,
    /// Source mode admission metadata; no source clock or state is advanced.
    /// The native causal cache still uses its explicitly supplied duration.
    pub sample_rate_hz: u32,
    pub modal: ModalAcousticTimeBudget,
    pub fluid: ThermoviscousLimits,
}

/// Source/geometry admission remains distinct from native numerical refusals.
#[derive(Debug)]
pub enum ShellWallError {
    Invalid(&'static str),
    Budget {
        quantity: &'static str,
        required: u128,
        maximum: u128,
    },
    Modal(ModalAcousticTimeError),
    Fluid(ThermoviscousSourceError),
    Allocation,
    Cancelled,
}
impl core::fmt::Display for ShellWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-COUPLE-SHELL-THERMOVISCOUS: {self:?}")
    }
}
impl std::error::Error for ShellWallError {}
impl From<ThermoviscousSourceError> for ShellWallError {
    fn from(value: ThermoviscousSourceError) -> Self {
        Self::Fluid(value)
    }
}
impl From<ThermoviscousError> for ShellWallError {
    fn from(value: ThermoviscousError) -> Self {
        Self::Fluid(value.into())
    }
}
impl From<ModalAcousticTimeError> for ShellWallError {
    fn from(value: ModalAcousticTimeError) -> Self {
        Self::Modal(value)
    }
}
fn poll(gate: &CancelGate) -> Result<(), ShellWallError> {
    if gate.is_requested() {
        Err(ShellWallError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit(quantity: &'static str, required: u128, maximum: u128) -> Result<(), ShellWallError> {
    if required > maximum {
        Err(ShellWallError::Budget {
            quantity,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn reserved<T>(count: usize) -> Result<Vec<T>, ShellWallError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| ShellWallError::Allocation)?;
    Ok(result)
}
fn canonical(mut face: [u32; 3]) -> [u32; 3] {
    face.sort_unstable();
    face
}

struct ShellBinding<'s> {
    image: &'s MaterialShellImage,
    vertex_map: Vec<Option<u32>>,
}

/// One complete material shell and native fluid operator before fluid binding.
/// Consume through `bind_gas`, `bind_relaxation` or `bind_material`; all reuse
/// the original source checks and retain the same immutable shell reference.
pub struct ShellWallAssembly<'s, 'c> {
    binding: ShellBinding<'s>,
    fe: ThermoviscousFe<'c>,
}

/// Original shell image, exact vertex correspondence and source-owned fluid.
/// Borrow `interior()` for existing harmonic, causal, field and duct operations.
pub struct ShellFluidInterior<'s, 'c, S = GasState> {
    binding: ShellBinding<'s>,
    interior: ThermoviscousInterior<'c, S>,
}

impl<'s, 'c> ShellWallAssembly<'s, 'c> {
    /// Assemble the exact source cavity selection on a conforming fluid mesh.
    ///
    /// `vertex_map` has one entry per ORIGINAL shell node. Every node of the
    /// selected cavity facets must map once to its coincident fluid vertex;
    /// other nodes must be None. Cartesian coordinates must agree exactly as
    /// represented. Relabeling is allowed; welding, offset skins and tolerance-
    /// based nearest-point matching are not. The source's default cavity
    /// selection means all its original triangles with their positive winding.
    /// Each mapped face must be a native exterior no-slip face, and its source
    /// side must point OUT of the fluid. A mismatch refuses without sign repair.
    ///
    /// Every source mode, material channel and interface dashpot is retained.
    /// Interface stiffness is already in the eigenbasis and is not added again.
    /// Source modes are admitted at the caller's clock with zero scalar loss
    /// and observer transfer before attaching the complete original damping.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        image: &'s MaterialShellImage,
        geometry: HelmholtzFe<'c>,
        vertex_map: &[Option<u32>],
        medium: ThermoviscousMedium,
        thermodynamics: ThermoviscousThermodynamics,
        boundary: &ThermoviscousBoundary,
        controls: ShellWallControls,
        gate: &CancelGate,
    ) -> Result<Self, ShellWallError> {
        poll(gate)?;
        let source = image.source();
        let mesh = source.mesh();
        let modes = image.physical_modes().len();
        let interfaces = source.interfaces();
        let selected = source.pressure_surfaces().map(|s| s.cavity.as_slice());
        let facet_count = selected.map_or(mesh.tris.len(), <[_]>::len);
        if geometry.order() != 2
            || modes == 0
            || facet_count == 0
            || facet_count > boundary.no_slip_faces.len()
            || vertex_map.len() != mesh.nodes.len()
            || image.radiators().len() != modes
            || image.damping().mode_count() != modes
            || image.interface_columns().len() != interfaces.len()
            || image.cavity_volume_shapes_m2_per_sqrt_kg().len() != modes
        {
            return Err(ShellWallError::Invalid(
                "complete nonempty original shell basis, cavity and P2 fluid geometry required",
            ));
        }
        admit(
            "source nodes",
            mesh.nodes.len() as u128,
            controls.maximum_source_nodes as u128,
        )?;
        admit(
            "source triangles",
            mesh.tris.len() as u128,
            controls.maximum_source_triangles as u128,
        )?;
        admit(
            "source interfaces",
            interfaces.len() as u128,
            controls.maximum_source_interfaces as u128,
        )?;
        admit(
            "original force ports",
            image.force_ports().len() as u128,
            controls.maximum_force_ports as u128,
        )?;
        admit(
            "original modes",
            modes as u128,
            controls.maximum_modes.min(MAX_MOVING_WALL_MODES) as u128,
        )?;
        if image.coating_relaxation().is_some()
            || image.substrate_relaxation().is_some()
            || !source.coatings().is_empty()
            || source
                .regions()
                .iter()
                .any(|region| region.constant_loss_reference_omega_rad_s().is_some())
        {
            return Err(ShellWallError::Invalid(
                "frequency-local coating/constant-loss slices or separate material memory need their original causal owner",
            ));
        }
        let channel_count = image
            .damping()
            .channel_count()
            .checked_add(
                interfaces
                    .iter()
                    .filter(|source| source.law().damping_n_s_m != 0.0)
                    .count(),
            )
            .ok_or(ShellWallError::Invalid("channel extent overflow"))?;
        admit(
            "complete damping channels",
            channel_count as u128,
            controls.maximum_damping_channels.min(MAX_MOVING_WALL_MODES) as u128,
        )?;
        let n = modes as u128;
        let s = mesh.nodes.len() as u128;
        let t = mesh.tris.len() as u128;
        let f = facet_count as u128;
        let d = geometry.ndof() as u128;
        let v = geometry.positions().len() as u128;
        let p = channel_count as u128;
        let bytes = 24 * n * d
            + 8 * p * n
            + 8 * n * n
            + 24 * (n + p)
            + 512 * n
            + (core::mem::size_of::<Option<u32>>() as u128 + 1) * s
            + core::mem::size_of::<Option<usize>>() as u128 * v
            + (core::mem::size_of::<OrientedShellFacet>() as u128 + 12) * f
            + 4096;
        let work = 64 * (s + t + interfaces.len() as u128 + v)
            + 64 * n * (s + d + f + image.force_ports().len() as u128)
            + 64 * f
                * (boundary.no_slip_faces.len() as u128 + geometry.tetrahedra().len() as u128 + 1)
            + 32 * p * n * n;
        admit(
            "shell mapping bytes",
            bytes,
            controls.maximum_setup_bytes as u128,
        )?;
        admit(
            "shell mapping work",
            work,
            controls.maximum_setup_work as u128,
        )?;
        admit(
            "fluid elements",
            geometry.tetrahedra().len() as u128,
            controls.fluid.max_elements as u128,
        )?;
        let mut facets = reserved(facet_count)?;
        if let Some(selected) = selected {
            facets.extend_from_slice(selected);
        } else {
            facets.extend((0..mesh.tris.len()).map(|triangle| OrientedShellFacet {
                triangle,
                side: ShellFace::Positive,
            }));
        }
        let mut active = reserved(mesh.nodes.len())?;
        active.resize(mesh.nodes.len(), false);
        for selected in &facets {
            poll(gate)?;
            let triangle = mesh
                .tris
                .get(selected.triangle)
                .ok_or(ShellWallError::Invalid(
                    "source cavity facet is not in original mesh",
                ))?;
            for &node in triangle {
                *active
                    .get_mut(node)
                    .ok_or(ShellWallError::Invalid("source cavity node out of range"))? = true;
            }
        }
        let mut reverse = reserved(geometry.positions().len())?;
        reverse.resize(geometry.positions().len(), None);
        for (node, (&mapped, &active)) in vertex_map.iter().zip(&active).enumerate() {
            poll(gate)?;
            if mapped.is_some() != active {
                return Err(ShellWallError::Invalid(
                    "vertex map must cover exactly the original cavity nodes",
                ));
            }
            if let Some(vertex) = mapped {
                let fluid = geometry
                    .positions()
                    .get(vertex as usize)
                    .ok_or(ShellWallError::Invalid("mapped fluid vertex out of range"))?;
                if *fluid != mesh.nodes[node] || reverse[vertex as usize].replace(node).is_some() {
                    return Err(ShellWallError::Invalid(
                        "source and fluid vertices must coincide one-to-one without welding",
                    ));
                }
            }
        }
        let mut faces = reserved(facet_count)?;
        let mut volume = vec![Dd::ZERO; modes];
        let mut volume_scale = vec![0.0_f64; modes];
        for mode in image.physical_modes() {
            poll(gate)?;
            if mode.len() != 6 * mesh.nodes.len() {
                return Err(ShellWallError::Invalid(
                    "finite complete original global shell modes required",
                ));
            }
            for values in mode.chunks(1024) {
                poll(gate)?;
                if values.iter().any(|v| !v.is_finite()) {
                    return Err(ShellWallError::Invalid(
                        "finite complete original global shell modes required",
                    ));
                }
            }
        }
        for selected in &facets {
            poll(gate)?;
            let triangle = mesh.tris[selected.triangle];
            let face =
                canonical(triangle.map(|node| vertex_map[node].expect("admitted cavity vertex")));
            if faces.contains(&face)
                || !boundary
                    .no_slip_faces
                    .iter()
                    .any(|&wall| canonical(wall) == face)
            {
                return Err(ShellWallError::Invalid(
                    "mapped cavity faces must be distinct original no-slip boundary faces",
                ));
            }
            let (_, normal, area) = geometry
                .boundary_geometry(face)
                .map_err(ThermoviscousError::from)?;
            let original = mesh
                .facet(selected.triangle)
                .map_err(|_| ShellWallError::Invalid("invalid original shell facet geometry"))?;
            let sign = match selected.side {
                ShellFace::Positive => 1.0,
                ShellFace::Negative => -1.0,
            };
            let cosine: f64 = (0..3)
                .map(|c| sign * original.frame[2][c] * normal[c])
                .sum();
            if !cosine.is_finite()
                || cosine < 1.0 - 1024.0 * f64::EPSILON
                || (area - original.area_m2).abs() > 1024.0 * f64::EPSILON * area
            {
                return Err(ShellWallError::Invalid(
                    "source cavity side must match the actual fluid outward face without an offset",
                ));
            }
            for (r, mode) in image.physical_modes().iter().enumerate() {
                for node in triangle {
                    for c in 0..3 {
                        let term = Dd::from_f64(area / 3.0)
                            * Dd::from_f64(normal[c])
                            * Dd::from_f64(mode[6 * node + c]);
                        volume[r] = volume[r] + term;
                        volume_scale[r] += term.to_f64().abs();
                    }
                }
            }
            faces.push(face);
        }
        for r in 0..modes {
            poll(gate)?;
            let expected = image.cavity_volume_shapes_m2_per_sqrt_kg()[r];
            let value = volume[r].to_f64();
            let allowance = 64.0
                * f64::EPSILON
                * (facet_count as f64 + 1.0)
                * volume_scale[r].max(expected.abs());
            if !value.is_finite()
                || !expected.is_finite()
                || !allowance.is_finite()
                || (value - expected).abs() > allowance
            {
                return Err(ShellWallError::Invalid(
                    "mapped fluid trace disagrees with original signed cavity volume work",
                ));
            }
        }
        let mut shapes = reserved(modes)?;
        let mut modal = reserved(modes)?;
        for (r, original) in image.physical_modes().iter().enumerate() {
            poll(gate)?;
            let mut shape = reserved(geometry.ndof())?;
            shape.resize(geometry.ndof(), [0.0; 3]);
            for (node, &vertex) in vertex_map.iter().enumerate() {
                if node % 1024 == 0 {
                    poll(gate)?;
                }
                if let Some(vertex) = vertex {
                    shape[vertex as usize].copy_from_slice(&original[6 * node..6 * node + 3]);
                }
            }
            shapes.push(shape);
            modal.push(ModalAcousticMode {
                angular_frequency_rad_s: image.radiators()[r].omega,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            });
        }
        for row in image.force_ports() {
            poll(gate)?;
            if row.len() != modes || row.iter().any(|v| !v.is_finite()) {
                return Err(ShellWallError::Invalid(
                    "original force ports must cover the same finite modal basis",
                ));
            }
        }
        let mut channels = reserved(channel_count)?;
        let resistance = image.damping().reference_resistance_n_s_m().sqrt();
        for original in image.damping().channel_shapes() {
            poll(gate)?;
            channels.push(original.iter().map(|&v| resistance * v).collect());
        }
        for (interface, row) in interfaces.iter().zip(image.interface_columns()) {
            poll(gate)?;
            let damping = interface.law().damping_n_s_m;
            if row.len() != modes
                || row.iter().any(|v| !v.is_finite())
                || !damping.is_finite()
                || damping < 0.0
            {
                return Err(ShellWallError::Invalid(
                    "finite complete original interface damping required",
                ));
            }
            if damping > 0.0 {
                let root = damping.sqrt();
                channels.push(row.iter().map(|&v| root * v).collect());
            }
        }
        let structure =
            ModalAcousticTimeModel::try_new(controls.sample_rate_hz, modal, controls.modal)?;
        let fe = ThermoviscousFe::from_geometry_with_moving_wall(
            geometry,
            medium,
            thermodynamics,
            boundary,
            MovingWall {
                faces,
                velocity_shapes: shapes,
                structure,
            },
            controls.fluid,
            &mut || gate.is_requested(),
        )?
        .with_moving_wall_damping_channels(channels, &mut || gate.is_requested())?;
        let mut owned_map = reserved(vertex_map.len())?;
        owned_map.extend_from_slice(vertex_map);
        poll(gate)?;
        Ok(Self {
            binding: ShellBinding {
                image,
                vertex_map: owned_map,
            },
            fe,
        })
    }

    pub fn native(&self) -> &ThermoviscousFe<'c> {
        &self.fe
    }
    pub fn image(&self) -> &'s MaterialShellImage {
        self.binding.image
    }

    /// Bind the same original gas; requested drive band must also be admitted
    /// by every original shell material and interface, not just its eigenmodes.
    pub fn bind_gas(
        self,
        gas: GasState,
        band_hz: [f64; 2],
        gate: &CancelGate,
    ) -> Result<ShellFluidInterior<'s, 'c>, ShellWallError> {
        self.binding.admit_band(band_hz, gate)?;
        let interior = ThermoviscousInterior::from_gas(self.fe, gas, band_hz)?;
        poll(gate)?;
        Ok(ShellFluidInterior {
            binding: self.binding,
            interior,
        })
    }

    /// Retain original gas molecular reservoirs with the same mechanical source.
    pub fn bind_relaxation(
        self,
        source: RelaxationState,
        band_hz: [f64; 2],
        gate: &CancelGate,
    ) -> Result<ShellFluidInterior<'s, 'c>, ShellWallError> {
        self.binding.admit_band(band_hz, gate)?;
        let interior = ThermoviscousInterior::from_relaxation(self.fe, source, band_hz)?;
        poll(gate)?;
        Ok(ShellFluidInterior {
            binding: self.binding,
            interior,
        })
    }

    /// Bind a generic original material fluid through its existing exact EOS,
    /// signed expansion, caloric and transport-coherence checks.
    pub fn bind_material<S: ThermoviscousFluid + ThermalFluidTransport>(
        self,
        source: S,
        band_hz: [f64; 2],
        gate: &CancelGate,
    ) -> Result<ShellFluidInterior<'s, 'c, S>, ShellWallError> {
        self.binding.admit_band(band_hz, gate)?;
        let interior = ThermoviscousInterior::from_material(self.fe, source, band_hz)?;
        poll(gate)?;
        Ok(ShellFluidInterior {
            binding: self.binding,
            interior,
        })
    }
}

impl ShellBinding<'_> {
    fn admit_band(&self, band_hz: [f64; 2], gate: &CancelGate) -> Result<(), ShellWallError> {
        poll(gate)?;
        let omega = band_hz.map(|hz| std::f64::consts::TAU * hz);
        if omega.iter().any(|v| !v.is_finite() || *v <= 0.0) || omega[0] >= omega[1] {
            return Err(ShellWallError::Invalid(
                "positive finite increasing shell/fluid band required",
            ));
        }
        for band in self
            .image
            .source()
            .regions()
            .iter()
            .map(|r| r.omega_band_rad_s())
            .chain(
                self.image
                    .source()
                    .interfaces()
                    .iter()
                    .map(|i| i.omega_band_rad_s()),
            )
        {
            poll(gate)?;
            if omega[0] < band.0 || omega[1] > band.1 {
                return Err(ShellWallError::Invalid(
                    "drive band exceeds original shell material or interface support",
                ));
            }
        }
        Ok(())
    }
}

impl<'s, 'c, S> ShellFluidInterior<'s, 'c, S> {
    pub fn image(&self) -> &'s MaterialShellImage {
        self.binding.image
    }
    pub fn vertex_map(&self) -> &[Option<u32>] {
        &self.binding.vertex_map
    }
    pub fn interior(&self) -> &ThermoviscousInterior<'c, S> {
        &self.interior
    }

    /// Continue existing consuming fluid-source attachments, including finite
    /// conductive walls, on this exact assembled interior. All native shell
    /// modes, traces and damping remain. This releases the wrapper's borrowed
    /// shell metadata; the caller retains its original image and force rows.
    pub fn into_interior(self) -> ThermoviscousInterior<'c, S> {
        self.interior
    }

    /// Project original physical force-port amplitudes [N] into generalized
    /// wall forces [N/sqrt(kg)] without allocation or basis/sign changes.
    /// Complex values use the native exp(+i omega t) convention. Real inputs
    /// give the same signed map for causal loads. Ports may include original
    /// moment/N lever arms; no extra lever or area is applied here.
    /// DD coherent sums retain cancellations. Invalid input or cancellation
    /// leaves `output` unchanged; source force-port count was bounded at setup.
    pub fn project_forces(
        &self,
        amplitudes: &[C64],
        output: &mut [C64],
        gate: &CancelGate,
    ) -> Result<(), ShellWallError> {
        poll(gate)?;
        let rows = self.binding.image.force_ports();
        let modes = self.binding.image.physical_modes().len();
        if amplitudes.len() != rows.len()
            || output.len() != modes
            || amplitudes
                .iter()
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(ShellWallError::Invalid(
                "finite amplitudes for original shell force ports and complete output required",
            ));
        }
        let mut candidate = [C64::ZERO; MAX_MOVING_WALL_MODES];
        for r in 0..modes {
            poll(gate)?;
            let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
            for (row, amplitude) in rows.iter().zip(amplitudes) {
                re = re + Dd::from_f64(row[r]) * Dd::from_f64(amplitude.re);
                im = im + Dd::from_f64(row[r]) * Dd::from_f64(amplitude.im);
            }
            candidate[r] = C64::new(re.to_f64(), im.to_f64());
            if !candidate[r].re.is_finite() || !candidate[r].im.is_finite() {
                return Err(ShellWallError::Invalid(
                    "original force projection is unrepresentable",
                ));
            }
        }
        poll(gate)?;
        output.copy_from_slice(&candidate[..modes]);
        Ok(())
    }
}
