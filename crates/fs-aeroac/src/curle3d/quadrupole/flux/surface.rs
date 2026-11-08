//! Stationary permeable FW-H sources from accepted compressible flow.
//!
//! A closed control box supplies, with its normal pointing OUT of the box,
//! `Q = rho (u.n) dA / rho0` and
//! `F = [(p-p0) n + rho u (u.n)] dA`.
//! These equivalent volume/force sources use the existing compact Green
//! operators and source-major `[Q,Fx,Fy,Fz]` Welch channels. They preserve
//! source phase and full near/far pressure and vector velocity. A permeable
//! surface's F includes convected momentum; it is not just wall pressure.
//!
//! The two original adjacent cell states supply their physical fluxes BEFORE
//! centered interpolation to each face. Face-midpoint area quadrature and
//! this reconstruction require spatial refinement; they are not the native
//! Riemann solver's dissipative numerical flux. No whole-box compactness is
//! assumed, but each face's wavelength/range resolution remains the caller's
//! responsibility. Sampling is O(faces), allocates no adapter storage, and
//! publishes all channels at the same actual accepted flow time.
//!
//! Scope is INVISCID, stationary Cartesian geometry and a quiescent homogeneous
//! exterior reference. The surface represents sources INSIDE the box. The
//! FW-H volume term OUTSIDE the box is still required when non-negligible;
//! hydrodynamic disturbances crossing the surface can otherwise create false
//! sound. Never add enclosed Lighthill volume stresses a second time. Only an
//! exterior acoustic observer may identify c0^2 rho' with microphone pressure.
//! No viscosity, SGS stress, moving geometry, wake correction, resampling,
//! source derivatives, radiation feedback or measured-SPL claim is inferred.
//!
//! Reference: Lockard & Casper (2005), equations (1)-(6),
//! <https://ntrs.nasa.gov/api/citations/20050182105/downloads/20050182105.pdf>.

use super::{
    Boundaries3d, CartesianGrid, CompressibleVolumeConfig3, CompressibleVolumeError3, Euler3d,
    FrozenGas, MixtureEuler3d, MixtureGas, Primitive, Source, admit_bytes, allocated, kernel, poll,
};
use fs_flux::compressible::spatial::mixture::{Cell, state_into};
use fs_flux::compressible::spatial::{BoundarySide, State};
use fs_math::dd::Dd;

/// Native joint recording of this surface and selected exterior volume stresses.
pub mod exterior;

/// Half-open original cell box. Its six bounding planes form the source surface.
/// A real cell must exist on BOTH sides of every face: periodic wrapping,
/// physical-wall ghosts and extrapolated exterior states are never inferred.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PermeableBox3 {
    pub lower: [usize; 3],
    pub upper: [usize; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct PermeableSurfaceControls3 {
    /// Bounds retained geometry and the two native EOS decodes per face/sample.
    pub maximum_faces: usize,
    /// Numerical geometry, source, trial and ambient-composition storage.
    /// Original flow/EOS owners and allocator metadata are separate.
    pub maximum_workspace_bytes: usize,
    /// Numerical agreement of native ambient rho/p/c with the supplied
    /// propagation reference, strictly positive and at most 1e-3. This admits
    /// EOS inversion roundoff; the supplied medium is never adjusted to fit.
    pub reference_relative_tolerance: f64,
}

/// One actual grid face in the fixed source order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PermeableFace3 {
    /// Original negative-axis and positive-axis neighboring cell indices.
    pub adjacent_cells: [usize; 2],
    /// Side of the CONTROL BOX, not a physical boundary of the flow domain.
    pub side: BoundarySide,
    pub area_m2: f64,
}

/// Complete instantaneous surface data, bound to the original gas and grid.
///
/// Ordering is XMin, XMax, YMin, YMax, ZMin, ZMax; within each side the lower
/// numbered tangential axis varies fastest. All four source components share
/// the same position, physical clock and source order. Sampling failure leaves
/// the preceding accepted sample intact. This owner does not advance the flow.
#[derive(Debug)]
pub struct PreparedPermeableSurface3 {
    source: Source,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    config: CompressibleVolumeConfig3,
    control_box: PermeableBox3,
    reference_partial_densities: Vec<f64>,
    faces: Vec<PermeableFace3>,
    positions: Vec<[f64; 3]>,
    sources: Vec<[f64; 4]>,
    trial: Vec<[f64; 4]>,
    time_s: Option<f64>,
}

impl PermeableBox3 {
    /// Complete closed surface count, after admitting a one-cell native halo.
    pub fn face_count(self, grid: &CartesianGrid) -> Result<usize, CompressibleVolumeError3> {
        let shape = grid.shape();
        if (0..3).any(|a| {
            self.lower[a] == 0 || self.lower[a] >= self.upper[a] || self.upper[a] >= shape[a]
        }) {
            return Err(invalid(
                "permeable box requires nonempty extents and a native cell on both sides",
            ));
        }
        let d: [usize; 3] = std::array::from_fn(|a| self.upper[a] - self.lower[a]);
        d[0].checked_mul(d[1])
            .and_then(|xy| d[0].checked_mul(d[2]).and_then(|xz| xy.checked_add(xz)))
            .and_then(|sum| d[1].checked_mul(d[2]).and_then(|yz| sum.checked_add(yz)))
            .and_then(|sum| sum.checked_mul(2))
            .ok_or_else(|| invalid("permeable surface extent overflow"))
    }
}

impl PreparedPermeableSurface3 {
    /// Conservative peak numerical payload. The final 32 bytes also cover the
    /// four mechanical entries of the temporary native mixture reference.
    pub fn required_workspace_bytes(
        faces: usize,
        reference_species: usize,
    ) -> Result<usize, CompressibleVolumeError3> {
        faces
            .checked_mul(
                size_of::<PermeableFace3>() + size_of::<[f64; 3]>() + 2 * size_of::<[f64; 4]>(),
            )
            .and_then(|bytes| {
                reference_species
                    .checked_mul(size_of::<f64>())
                    .and_then(|r| bytes.checked_add(r))
            })
            .and_then(|bytes| bytes.checked_add(32))
            .ok_or_else(|| invalid("permeable surface workspace overflow"))
    }

    /// Retain the exact frozen source. Native ambient reconstruction checks
    /// the explicitly supplied rho0, p0 and c0 before publishing any sample.
    pub fn from_frozen(
        flow: &Euler3d,
        control_box: PermeableBox3,
        config: CompressibleVolumeConfig3,
        controls: PermeableSurfaceControls3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        let count = Self::admit(
            flow.grid(),
            control_box,
            config,
            controls,
            0,
            &mut cancelled,
        )?;
        let reference = State::from_primitive(
            flow.gas(),
            config.acoustic.density_kg_m3(),
            [0.0; 3],
            config.reference_pressure_pa,
        )
        .map_err(frozen_error)?
        .primitive(flow.gas(), &mut cancelled)
        .map_err(frozen_error)?;
        admit_reference(reference, config, controls.reference_relative_tolerance)?;
        let mut result = Self::prepare(
            Source::Frozen(flow.gas().clone()),
            *flow.grid(),
            flow.boundaries(),
            control_box,
            config,
            &[],
            count,
            controls.maximum_workspace_bytes,
            &mut cancelled,
        )?;
        result.sample_frozen(flow, &mut cancelled)?;
        Ok(result)
    }

    /// Ambient composition is explicit in the original source's species order,
    /// including exact zeros. Its sum must match rho0 and its native sound speed
    /// at p0 must match c0. Current surface cells retain their OWN compositions.
    pub fn from_mixture(
        flow: &MixtureEuler3d,
        reference_partial_densities_kg_m3: &[f64],
        control_box: PermeableBox3,
        config: CompressibleVolumeConfig3,
        controls: PermeableSurfaceControls3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        let species = flow.gas().species_count();
        let count = Self::admit(
            flow.grid(),
            control_box,
            config,
            controls,
            species,
            &mut cancelled,
        )?;
        if reference_partial_densities_kg_m3.len() != species {
            return Err(invalid(
                "complete original ambient species-density axis required",
            ));
        }
        let mut state = allocated(species + 4, 0.0)?;
        state_into(
            flow.gas(),
            reference_partial_densities_kg_m3,
            [0.0; 3],
            config.reference_pressure_pa,
            &mut state,
        )
        .map_err(mixture_error)?;
        let reference = Cell {
            partial_densities: &state[..species],
            momentum: [state[species], state[species + 1], state[species + 2]],
            energy: state[species + 3],
        }
        .primitive(flow.gas(), &mut cancelled)
        .map_err(mixture_error)?;
        admit_reference(reference, config, controls.reference_relative_tolerance)?;
        drop(state);
        let mut result = Self::prepare(
            Source::Mixture(flow.gas().clone()),
            *flow.grid(),
            flow.boundaries(),
            control_box,
            config,
            reference_partial_densities_kg_m3,
            count,
            controls.maximum_workspace_bytes,
            &mut cancelled,
        )?;
        result.sample_mixture(flow, &mut cancelled)?;
        Ok(result)
    }

    fn admit(
        grid: &CartesianGrid,
        control_box: PermeableBox3,
        config: CompressibleVolumeConfig3,
        controls: PermeableSurfaceControls3,
        species: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize, CompressibleVolumeError3> {
        poll(cancelled)?;
        let count = control_box.face_count(grid)?;
        if count > controls.maximum_faces {
            return Err(invalid("complete permeable surface exceeds face budget"));
        }
        admit_bytes(
            Self::required_workspace_bytes(count, species)?,
            controls.maximum_workspace_bytes,
        )?;
        if !config.reference_pressure_pa.is_finite()
            || config.reference_pressure_pa <= 0.0
            || config.origin_m.iter().any(|x| !x.is_finite())
            || !kernel::proper_rotation(config.grid_to_acoustic)
            || !controls.reference_relative_tolerance.is_finite()
            || controls.reference_relative_tolerance <= 0.0
            || controls.reference_relative_tolerance > 1e-3
        {
            return Err(invalid(
                "finite reference pressure/frame and tight native reference tolerance required",
            ));
        }
        Ok(count)
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        source: Source,
        grid: CartesianGrid,
        boundaries: Boundaries3d,
        control_box: PermeableBox3,
        config: CompressibleVolumeConfig3,
        reference_partial_densities: &[f64],
        count: usize,
        maximum: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        let mut result = Self {
            source,
            grid,
            boundaries,
            config,
            control_box,
            reference_partial_densities: allocated(reference_partial_densities.len(), 0.0)?,
            faces: allocated(
                count,
                PermeableFace3 {
                    adjacent_cells: [0; 2],
                    side: BoundarySide::XMin,
                    area_m2: 0.0,
                },
            )?,
            positions: allocated(count, [0.0; 3])?,
            sources: allocated(count, [0.0; 4])?,
            trial: allocated(count, [0.0; 4])?,
            time_s: None,
        };
        admit_bytes(result.workspace_bytes(), maximum)?;
        result
            .reference_partial_densities
            .copy_from_slice(reference_partial_densities);
        let shape = grid.shape();
        let widths = grid.spacing();
        let index = |c: [usize; 3]| (c[2] * shape[1] + c[1]) * shape[0] + c[0];
        let mut next = 0;
        for side in BoundarySide::ALL {
            let axis = side.axis();
            let [a, b] = match axis {
                0 => [1, 2],
                1 => [0, 2],
                _ => [0, 1],
            };
            let plane = if side.sign() < 0.0 {
                control_box.lower[axis]
            } else {
                control_box.upper[axis]
            };
            for j in control_box.lower[b]..control_box.upper[b] {
                for i in control_box.lower[a]..control_box.upper[a] {
                    if next % 64 == 0 {
                        poll(cancelled)?;
                    }
                    let mut left = [0; 3];
                    left[axis] = plane - 1;
                    left[a] = i;
                    left[b] = j;
                    let mut right = left;
                    right[axis] += 1;
                    let mut local = grid
                        .cell_center(index(left))
                        .ok_or(CompressibleVolumeError3::GeometryMismatch)?;
                    local[axis] = plane as f64 * widths[axis];
                    let rotated = kernel::rotate(config.grid_to_acoustic, local);
                    let position = std::array::from_fn(|d| config.origin_m[d] + rotated[d]);
                    if position.iter().any(|x| !x.is_finite()) {
                        return Err(invalid("unrepresentable permeable source position"));
                    }
                    result.faces[next] = PermeableFace3 {
                        adjacent_cells: [index(left), index(right)],
                        side,
                        area_m2: grid.face_area(side),
                    };
                    result.positions[next] = position;
                    next += 1;
                }
            }
        }
        poll(cancelled)?;
        Ok(result)
    }

    pub fn faces(&self) -> &[PermeableFace3] {
        &self.faces
    }
    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.positions
    }
    /// Instantaneous source-major `[Q,Fx,Fy,Fz]` in m^3/s and N. These are
    /// physical histories, NOT phasors. Transform every channel on the original
    /// clock before using the native harmonic or coherent Welch field owners.
    /// Keep all faces in one record; per-face intensity addition loses work.
    pub fn source_components(&self) -> &[[f64; 4]] {
        &self.sources
    }
    pub fn sample_time_s(&self) -> f64 {
        self.time_s
            .expect("constructor publishes initial surface sample")
    }
    pub const fn config(&self) -> CompressibleVolumeConfig3 {
        self.config
    }
    /// Begin one joint spectrum with the CURRENT accepted physical sample.
    /// The supplied clock must place that sample at index zero. All original
    /// face Q/F channels share one window/transform; no independent-face power
    /// approximation or extra source derivative is introduced.
    ///
    /// This preserves the stationary inviscid FW-H source scope above. It
    /// does not add omitted exterior volume stresses, antialias adaptive CFD,
    /// or validate disturbances crossing the permeable surface as sound.
    pub fn begin_recording(
        &self,
        clock: crate::curle3d::array::recorded::SourceSampleClock,
        expected_samples: usize,
        config: fs_fft::welch::WelchConfig,
        budget: crate::curle3d::array::recorded::CompactSourceRecordingBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<crate::curle3d::array::recorded::CompactSourceWelchStream, crate::AeroacError> {
        let mut recording = crate::curle3d::array::recorded::CompactSourceWelchStream::new(
            self.config.acoustic,
            &self.positions,
            clock,
            expected_samples,
            config,
            budget,
            &mut cancelled,
        )?;
        self.record_sample(&mut recording, &mut cancelled)?;
        Ok(recording)
    }

    /// Record the actual accepted mass/momentum sample on its actual clock.
    /// A different acoustic reference, geometry/order, repeated timestamp or
    /// cancellation leaves the previously accepted recording unchanged.
    pub fn record_sample(
        &self,
        recording: &mut crate::curle3d::array::recorded::CompactSourceWelchStream,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), crate::AeroacError> {
        if recording.acoustic() != self.config.acoustic {
            return Err(crate::AeroacError::InvalidParameter {
                what: "permeable surface and source recording acoustic references differ",
            });
        }
        recording.push_sample(
            &self.positions,
            self.sample_time_s(),
            &self.sources,
            cancelled,
        )
    }
    pub const fn control_box(&self) -> PermeableBox3 {
        self.control_box
    }
    /// Explicit original ambient composition; empty for a frozen source.
    pub fn reference_partial_densities_kg_m3(&self) -> &[f64] {
        &self.reference_partial_densities
    }
    pub fn workspace_bytes(&self) -> usize {
        self.faces.capacity() * size_of::<PermeableFace3>()
            + self.positions.capacity() * size_of::<[f64; 3]>()
            + (self.sources.capacity() + self.trial.capacity()) * size_of::<[f64; 4]>()
            + self.reference_partial_densities.capacity() * size_of::<f64>()
    }
    pub fn frozen_source(&self) -> Option<&FrozenGas> {
        match &self.source {
            Source::Frozen(source) => Some(source),
            Source::Mixture(_) => None,
        }
    }
    pub fn mixture_source(&self) -> Option<&MixtureGas> {
        match &self.source {
            Source::Mixture(source) => Some(source),
            Source::Frozen(_) => None,
        }
    }

    /// Publish only a strictly later, complete accepted frozen-gas sample.
    pub fn sample_frozen(
        &mut self,
        flow: &Euler3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(&mut cancelled)?;
        if self.frozen_source() != Some(flow.gas()) {
            return Err(CompressibleVolumeError3::SourceMismatch);
        }
        self.admit_sample(flow.grid(), flow.boundaries(), flow.time().value())?;
        self.sample_cells(flow.time().value(), &mut cancelled, |cell, cancel| {
            let state = flow.states()[cell];
            Ok((
                state.primitive(flow.gas(), cancel).map_err(frozen_error)?,
                state.momentum,
            ))
        })
    }

    /// Publish current compositions through the same bound original species laws.
    pub fn sample_mixture(
        &mut self,
        flow: &MixtureEuler3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(&mut cancelled)?;
        if self.mixture_source() != Some(flow.gas()) {
            return Err(CompressibleVolumeError3::SourceMismatch);
        }
        self.admit_sample(flow.grid(), flow.boundaries(), flow.time().value())?;
        self.sample_cells(flow.time().value(), &mut cancelled, |cell, cancel| {
            let state = flow
                .cell(cell)
                .ok_or(CompressibleVolumeError3::GeometryMismatch)?;
            Ok((
                state.primitive(flow.gas(), cancel).map_err(mixture_error)?,
                state.momentum,
            ))
        })
    }

    fn admit_sample(
        &self,
        grid: &CartesianGrid,
        boundaries: Boundaries3d,
        time: f64,
    ) -> Result<(), CompressibleVolumeError3> {
        if *grid != self.grid || boundaries != self.boundaries {
            return Err(CompressibleVolumeError3::GeometryMismatch);
        }
        if !time.is_finite() || time < 0.0 {
            return Err(invalid("invalid native accepted surface time"));
        }
        if let Some(previous_s) = self.time_s
            && time <= previous_s
        {
            return Err(CompressibleVolumeError3::SampleClock {
                previous_s,
                candidate_s: time,
            });
        }
        Ok(())
    }

    fn sample_cells<C: FnMut() -> bool>(
        &mut self,
        time: f64,
        cancelled: &mut C,
        mut cell: impl FnMut(usize, &mut C) -> Result<(Primitive, [f64; 3]), CompressibleVolumeError3>,
    ) -> Result<(), CompressibleVolumeError3> {
        for (i, face) in self.faces.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let states = [
                cell(face.adjacent_cells[0], cancelled)?,
                cell(face.adjacent_cells[1], cancelled)?,
            ];
            self.trial[i] = surface_flux(face, self.config, states)?;
        }
        poll(cancelled)?;
        std::mem::swap(&mut self.sources, &mut self.trial);
        self.time_s = Some(time);
        Ok(())
    }
}

fn surface_flux(
    face: &PermeableFace3,
    config: CompressibleVolumeConfig3,
    states: [(Primitive, [f64; 3]); 2],
) -> Result<[f64; 4], CompressibleVolumeError3> {
    let axis = face.side.axis();
    let half_area = Dd::from_f64(face.area_m2) * Dd::from_f64(0.5 * face.side.sign());
    let mut mass = Dd::from_f64(0.0);
    let mut loading = [Dd::from_f64(0.0); 3];
    for (state, momentum) in states {
        mass = mass + Dd::from_f64(momentum[axis]);
        for d in 0..3 {
            loading[d] =
                loading[d] + Dd::from_f64(momentum[d]) * Dd::from_f64(state.velocity[axis]);
        }
        // Subtract ambient before combining the two pressure traces. Physical
        // conserved momentum supplies rho*u directly; no mass flux is inferred
        // from the lattice or from a face-averaged velocity/density product.
        loading[axis] = loading[axis]
            + (Dd::from_f64(state.pressure) - Dd::from_f64(config.reference_pressure_pa));
    }
    let mut out = [0.0; 4];
    out[0] = ((mass * half_area) / Dd::from_f64(config.acoustic.density_kg_m3())).to_f64();
    for (d, row) in config.grid_to_acoustic.iter().enumerate() {
        let sum = row
            .iter()
            .zip(loading)
            .fold(Dd::from_f64(0.0), |sum, (&r, f)| sum + Dd::from_f64(r) * f);
        out[d + 1] = (sum * half_area).to_f64();
    }
    if out.iter().any(|v| !v.is_finite()) {
        return Err(invalid("unrepresentable permeable mass/momentum source"));
    }
    Ok(out)
}

fn admit_reference(
    reference: Primitive,
    config: CompressibleVolumeConfig3,
    tolerance: f64,
) -> Result<(), CompressibleVolumeError3> {
    for (actual, expected) in [
        (reference.density, config.acoustic.density_kg_m3()),
        (reference.pressure, config.reference_pressure_pa),
        (reference.sound_speed, config.acoustic.sound_speed_m_s()),
    ] {
        if !actual.is_finite()
            || actual <= 0.0
            || (actual - expected).abs() > tolerance * actual.abs().max(expected.abs())
        {
            return Err(invalid(
                "native ambient source and propagation rho0/p0/c0 differ",
            ));
        }
    }
    Ok(())
}
fn frozen_error(error: fs_flux::compressible::Error) -> CompressibleVolumeError3 {
    if error == fs_flux::compressible::Error::Cancelled {
        CompressibleVolumeError3::Cancelled
    } else {
        CompressibleVolumeError3::Frozen(error)
    }
}
fn mixture_error(
    error: fs_flux::compressible::spatial::mixture::Error,
) -> CompressibleVolumeError3 {
    if error == fs_flux::compressible::spatial::mixture::Error::Cancelled {
        CompressibleVolumeError3::Cancelled
    } else {
        CompressibleVolumeError3::Mixture(error)
    }
}
fn invalid(what: &'static str) -> CompressibleVolumeError3 {
    CompressibleVolumeError3::Invalid(what)
}

#[cfg(test)]
mod tests;
