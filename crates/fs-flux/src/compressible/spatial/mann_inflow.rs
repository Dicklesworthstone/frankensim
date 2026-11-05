//! Coherent Mann velocity traces for physical compressible reservoir faces.
//!
//! The original atmosphere owns the full spectral tensor and stochastic path
//! (Mann 1998, doi:10.1016/S0266-8920(97)00036-2). This adapter observes that
//! same realization on the native Cartesian faces. Original gas sources still
//! own calorics and Riemann fluxes; the prescribed baseline owns density,
//! composition, pressure and mean velocity. No per-face random draws, RMS
//! restoration, invented thermodynamic fluctuations or flux correction occur.
//!
//! `apply_*` publishes a trace held through the existing SSP stages and retries.
//! `step_*` instead observes a frozen advected realization at the actual two
//! SSP stage times, including every retry, and publishes the accepted endpoint
//! with the fluid. `step_*_bridge` observes a prepared conditional OU path
//! between the original fixed-clock endpoints at those same physical times.
//! These are complete-state
//! exterior Riemann data, not nonreflecting or characteristic inlet conditions.

use super::mixture::{
    Error as MixtureError, MixtureEuler3d, ReservoirFaceState as MixtureReservoirFaceState,
};
use super::{
    AxisBoundaries, Boundaries3d, BoundarySide, CartesianGrid, Error, Euler3d, ReservoirFaceState,
};
use fs_atmo::mann_field::filter::MannBoxFilter;
use fs_atmo::mann_field::fixed_step::{FixedStepMannField, MannOuBridge, MannOuBridgeFrame};
use fs_atmo::mann_field::{HomogeneousMannField, MannSpectralModel};

mod geometry;
use geometry::{half_edges, rotation_is_proper, speed, to_field, to_grid};
mod stages;
pub use stages::MannReservoirStepReport3;
mod bridge;
pub use bridge::MannMixtureForce3;
pub use bridge::MannReservoirBridgeReport3;

/// The Cartesian grid uses metres already. Its origin is a cell CORNER, unlike
/// the on-site lattice-node origin used by the D3Q19 adapter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannReservoirConfig3 {
    /// Field-frame position of Cartesian corner [0,0,0], in metres.
    pub grid_origin_in_field_m: [f64; 3],
    /// Proper orthonormal rotation: field vectors -> grid vectors.
    pub field_to_grid: [[f64; 3]; 3],
    /// Maximum face observations times represented wavevectors per apply or
    /// entire `step_*`, including every native stage and rejected attempt.
    pub maximum_mode_samples: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MannReservoirObservation3 {
    Point,
    /// Exact rectangular spatial average of velocity over each native face.
    /// This is not an average of nonlinear momentum or energy fluxes.
    FaceAverage,
}

/// Borrow one existing source. Neither observation nor publication advances it.
#[derive(Clone, Copy, Debug)]
pub enum MannReservoirSource3<'a> {
    /// Freeze the supplied CURRENT snapshot at the authored reference time,
    /// then translate it. The source's already accumulated phase is retained.
    Frozen {
        source: &'a HomogeneousMannField,
        reference_time_s: f64,
        /// Convection velocity in grid Cartesian axes [m/s].
        convection_velocity_m_s: [f64; 3],
    },
    /// Observe an accepted OU endpoint at exactly the native fluid time.
    /// Caller coordinates its fixed clock with accepted fluid timesteps.
    Evolving(&'a FixedStepMannField),
    /// One accepted observation of a prepared conditional OU bridge. This is
    /// a held trace; use `step_*_bridge` to sample changing integration times.
    BridgeFrame(&'a MannOuBridgeFrame),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MannReservoirAction3 {
    /// Allocate the native owner's active/trial profiles once, within this
    /// additional storage bound. Existing prescribed sides cannot overlap.
    Prepare { maximum_profile_bytes: usize },
    /// Refresh all and only the already prepared sides, without numerical
    /// allocation. Failed admission/cancellation leaves accepted traces intact.
    Update,
}

#[derive(Debug)]
pub enum MannReservoirError3 {
    Input(&'static str),
    SpatialNyquist {
        axis: usize,
        phase_bound_rad: f64,
    },
    /// Bound on advected spatial phases only. OU noise is not bandlimited.
    AdvectiveNyquist {
        phase_bound_rad: f64,
    },
    Work {
        required: usize,
        maximum: usize,
    },
    Atmosphere(fs_atmo::Refusal),
    Flow(Error),
    Mixture(MixtureError),
    Molecular(super::molecular::MolecularError),
    MixtureMolecular(fs_gas::multicomponent::dissipation::DissipationError),
    Reactive(fs_gas::multicomponent::reactive::ReactiveError),
    /// Conditional OU normal variates across every stage and rejected trial.
    BridgeWork {
        required: usize,
        maximum: usize,
    },
}
impl core::fmt::Display for MannReservoirError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-FLUX-MANN-RESERVOIR3: {self:?}")
    }
}
impl std::error::Error for MannReservoirError3 {}
impl From<fs_atmo::Refusal> for MannReservoirError3 {
    fn from(error: fs_atmo::Refusal) -> Self {
        Self::Atmosphere(error)
    }
}
impl From<Error> for MannReservoirError3 {
    fn from(error: Error) -> Self {
        Self::Flow(error)
    }
}
impl From<MixtureError> for MannReservoirError3 {
    fn from(error: MixtureError) -> Self {
        Self::Mixture(error)
    }
}
impl From<super::molecular::MolecularError> for MannReservoirError3 {
    fn from(error: super::molecular::MolecularError) -> Self {
        Self::Molecular(error)
    }
}
impl From<fs_gas::multicomponent::dissipation::DissipationError> for MannReservoirError3 {
    fn from(error: fs_gas::multicomponent::dissipation::DissipationError) -> Self {
        Self::MixtureMolecular(error)
    }
}
impl From<fs_gas::multicomponent::reactive::ReactiveError> for MannReservoirError3 {
    fn from(error: fs_gas::multicomponent::reactive::ReactiveError) -> Self {
        Self::Reactive(error)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannReservoirReport3 {
    pub sampled_time_s: f64,
    /// Caller must refresh no later than this bound; no fluid step is taken here.
    pub maximum_hold_s: f64,
    pub face_count: usize,
    pub mode_samples: usize,
    pub spatial_phase_bound_rad: [f64; 3],
    pub advective_phase_bound_rad: f64,
}

/// Reusable velocity observation shared by frozen and evolving-composition gas.
/// Opposite faces share one prepared filter for their common normal axis.
pub struct MannReservoirInflow3 {
    grid: CartesianGrid,
    config: MannReservoirConfig3,
    sides: [BoundarySide; 6],
    side_count: usize,
    selected: [bool; 6],
    filters: [Option<MannBoxFilter>; 3],
    k_max_per_m: f64,
    faces: usize,
    mode_samples: usize,
    spatial_phase_bound_rad: [f64; 3],
}

/// An immutable observation at one admitted fluid time. It can also feed the
/// native molecular wrappers' profile callbacks, with exactly the same samples.
pub struct MannReservoirTrace3<'a> {
    inlet: &'a MannReservoirInflow3,
    source: MannReservoirSource3<'a>,
    displacement_m: [f64; 3],
    report: MannReservoirReport3,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MannReservoirError3> {
    if cancelled() {
        Err(Error::Cancelled.into())
    } else {
        Ok(())
    }
}

impl MannReservoirInflow3 {
    pub fn new(
        model: &MannSpectralModel,
        grid: CartesianGrid,
        sides: &[BoundarySide],
        config: MannReservoirConfig3,
        observation: MannReservoirObservation3,
        maximum_filter_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MannReservoirError3> {
        poll(&mut cancelled)?;
        if sides.is_empty()
            || sides.len() > 6
            || config.maximum_mode_samples == 0
            || !config.grid_origin_in_field_m.iter().all(|x| x.is_finite())
            || !rotation_is_proper(config.field_to_grid)
        {
            return Err(MannReservoirError3::Input(
                "finite proper frame, nonempty unique sides and work bound required",
            ));
        }
        let mut selected = [false; 6];
        let mut faces = 0usize;
        for &side in sides {
            let slot = 2 * side.axis() + usize::from(side.sign() > 0.0);
            if selected[slot] {
                return Err(MannReservoirError3::Input("duplicate inlet side"));
            }
            selected[slot] = true;
            faces = faces
                .checked_add(grid.face_count(side))
                .ok_or(MannReservoirError3::Input("face count overflow"))?;
        }
        let mode_samples = faces
            .checked_mul(model.mode_count())
            .ok_or(MannReservoirError3::Input("inlet work overflow"))?;
        if mode_samples > config.maximum_mode_samples {
            return Err(MannReservoirError3::Work {
                required: mode_samples,
                maximum: config.maximum_mode_samples,
            });
        }
        let k_max_per_m = model.quadrature().max_kl / model.parameters().length_m;
        let spatial_phase_bound_rad = grid.spacing().map(|dx| k_max_per_m * dx);
        for (axis, phase) in spatial_phase_bound_rad.into_iter().enumerate() {
            if !phase.is_finite() || phase > core::f64::consts::PI {
                return Err(MannReservoirError3::SpatialNyquist {
                    axis,
                    phase_bound_rad: phase,
                });
            }
        }
        let mut filters: [Option<MannBoxFilter>; 3] = core::array::from_fn(|_| None);
        let mut remaining = maximum_filter_bytes;
        for axis in 0..3 {
            if selected[2 * axis] || selected[2 * axis + 1] {
                let filter = MannBoxFilter::prepare(
                    model,
                    half_edges(
                        config.field_to_grid,
                        grid.spacing(),
                        axis,
                        observation == MannReservoirObservation3::Point,
                    ),
                    remaining,
                    &mut cancelled,
                )?;
                remaining = remaining
                    .checked_sub(filter.workspace_bytes())
                    .ok_or(MannReservoirError3::Input("filter capacity exceeds budget"))?;
                filters[axis] = Some(filter);
            }
        }
        let mut canonical = BoundarySide::ALL;
        let mut side_count = 0;
        for (slot, side) in BoundarySide::ALL.into_iter().enumerate() {
            if selected[slot] {
                canonical[side_count] = side;
                side_count += 1;
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            grid,
            config,
            sides: canonical,
            side_count,
            selected,
            filters,
            k_max_per_m,
            faces,
            mode_samples,
            spatial_phase_bound_rad,
        })
    }

    pub fn sides(&self) -> &[BoundarySide] {
        &self.sides[..self.side_count]
    }
    pub const fn config(&self) -> MannReservoirConfig3 {
        self.config
    }
    pub fn workspace_bytes(&self) -> usize {
        self.filters
            .iter()
            .flatten()
            .map(MannBoxFilter::workspace_bytes)
            .sum()
    }
    /// The original filter also supports the atmosphere's exact covariance API.
    pub fn spatial_filter(&self, side: BoundarySide) -> Option<&MannBoxFilter> {
        let slot = 2 * side.axis() + usize::from(side.sign() > 0.0);
        self.selected[slot]
            .then(|| self.filters[side.axis()].as_ref())
            .flatten()
    }

    /// Bind the actual native grid, topology and physical time before sampling.
    /// Arbitrary Mann modes are not periodic on a CFD cross-section; transverse
    /// periodic axes refuse. Wall-touching mean profiles remain caller-authored.
    pub fn trace<'a>(
        &'a self,
        grid: &CartesianGrid,
        boundaries: Boundaries3d,
        fluid_time_s: f64,
        maximum_hold_s: f64,
        source: MannReservoirSource3<'a>,
    ) -> Result<MannReservoirTrace3<'a>, MannReservoirError3> {
        self.validate_geometry(grid, boundaries)?;
        if !fluid_time_s.is_finite()
            || !maximum_hold_s.is_finite()
            || maximum_hold_s <= 0.0
            || !(fluid_time_s + maximum_hold_s).is_finite()
            || fluid_time_s + maximum_hold_s <= fluid_time_s
        {
            return Err(MannReservoirError3::Input(
                "original grid and representable fluid time/hold bound required",
            ));
        }
        let (displacement_m, advection) = match source {
            MannReservoirSource3::Frozen {
                reference_time_s,
                convection_velocity_m_s,
                ..
            } => {
                let elapsed = fluid_time_s - reference_time_s;
                if !reference_time_s.is_finite()
                    || !elapsed.is_finite()
                    || !convection_velocity_m_s.iter().all(|x| x.is_finite())
                {
                    return Err(MannReservoirError3::Input(
                        "finite frozen epoch and convection velocity required",
                    ));
                }
                (
                    convection_velocity_m_s.map(|v| -v * elapsed),
                    convection_velocity_m_s,
                )
            }
            MannReservoirSource3::Evolving(field) => {
                if field.time_s() != fluid_time_s {
                    return Err(MannReservoirError3::Input(
                        "accepted Mann endpoint differs from native fluid time",
                    ));
                }
                // Norm is frame invariant, so the source's field-frame vector
                // supplies the same conservative advective phase bound.
                ([0.0; 3], field.advection_mps())
            }
            MannReservoirSource3::BridgeFrame(frame) => {
                if frame.time_s() != fluid_time_s {
                    return Err(MannReservoirError3::Input(
                        "Mann bridge frame differs from native fluid time",
                    ));
                }
                ([0.0; 3], frame.advection_mps())
            }
        };
        if !displacement_m.iter().all(|x| x.is_finite()) {
            return Err(MannReservoirError3::Input(
                "frozen displacement is not representable",
            ));
        }
        let advective_phase_bound_rad = self.k_max_per_m * speed(advection) * maximum_hold_s;
        if !advective_phase_bound_rad.is_finite()
            || advective_phase_bound_rad > core::f64::consts::PI
        {
            return Err(MannReservoirError3::AdvectiveNyquist {
                phase_bound_rad: advective_phase_bound_rad,
            });
        }
        Ok(MannReservoirTrace3 {
            inlet: self,
            source,
            displacement_m,
            report: MannReservoirReport3 {
                sampled_time_s: fluid_time_s,
                maximum_hold_s,
                face_count: self.faces,
                mode_samples: self.mode_samples,
                spatial_phase_bound_rad: self.spatial_phase_bound_rad,
                advective_phase_bound_rad,
            },
        })
    }

    /// Observe and atomically publish source-admitted frozen-composition gas.
    /// `baseline` supplies actual density, pressure and mean grid velocity at
    /// every physical face; the shared fluctuation is added only to velocity.
    pub fn apply_frozen(
        &self,
        flow: &mut Euler3d,
        source: MannReservoirSource3<'_>,
        maximum_hold_s: f64,
        action: MannReservoirAction3,
        mut baseline: impl FnMut(BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirReport3, MannReservoirError3> {
        poll(&mut cancelled)?;
        let trace = self.trace(
            flow.grid(),
            flow.boundaries(),
            flow.time().value(),
            maximum_hold_s,
            source,
        )?;
        if action == MannReservoirAction3::Update
            && flow.prepared_reservoir_profile_sides() != self.selected
        {
            return Err(MannReservoirError3::Input(
                "update requires exactly the prepared Mann sides",
            ));
        }
        let mut sample_error = None;
        let sample = |side, face, point| {
            let mean = baseline(side, face, point)?;
            trace.frozen_state(side, face, mean).map_err(|error| {
                sample_error = Some(error);
                Error::Invalid("Mann reservoir observation refused")
            })
        };
        let result = match action {
            MannReservoirAction3::Prepare {
                maximum_profile_bytes,
            } => flow.prepare_reservoir_profiles(
                self.sides(),
                maximum_profile_bytes,
                sample,
                &mut cancelled,
            ),
            MannReservoirAction3::Update => flow.update_reservoir_profiles(sample, &mut cancelled),
        };
        match result {
            Ok(()) => Ok(trace.report()),
            Err(error) => Err(sample_error.unwrap_or_else(|| error.into())),
        }
    }

    /// Same observation with every original species supplied by `baseline` in
    /// the canonical partial-density slice, including explicit absent species.
    pub fn apply_mixture(
        &self,
        flow: &mut MixtureEuler3d,
        source: MannReservoirSource3<'_>,
        maximum_hold_s: f64,
        action: MannReservoirAction3,
        mut baseline: impl FnMut(
            BoundarySide,
            usize,
            [f64; 3],
            &mut [f64],
        ) -> Result<MixtureReservoirFaceState, MixtureError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MannReservoirReport3, MannReservoirError3> {
        poll(&mut cancelled)?;
        let trace = self.trace(
            flow.grid(),
            flow.boundaries(),
            flow.time().value(),
            maximum_hold_s,
            source,
        )?;
        if action == MannReservoirAction3::Update
            && flow.prepared_reservoir_profile_sides() != self.selected
        {
            return Err(MannReservoirError3::Input(
                "update requires exactly the prepared Mann sides",
            ));
        }
        let mut sample_error = None;
        let sample = |side, face, point, densities: &mut [f64]| {
            let mean = baseline(side, face, point, densities)?;
            trace.mixture_state(side, face, mean).map_err(|error| {
                sample_error = Some(error);
                MixtureError::Invalid("Mann reservoir observation refused")
            })
        };
        let result = match action {
            MannReservoirAction3::Prepare {
                maximum_profile_bytes,
            } => flow.prepare_reservoir_profiles(
                self.sides(),
                maximum_profile_bytes,
                sample,
                &mut cancelled,
            ),
            MannReservoirAction3::Update => flow.update_reservoir_profiles(sample, &mut cancelled),
        };
        match result {
            Ok(()) => Ok(trace.report()),
            Err(error) => Err(sample_error.unwrap_or_else(|| error.into())),
        }
    }
}

impl MannReservoirInflow3 {
    fn validate_geometry(
        &self,
        grid: &CartesianGrid,
        boundaries: Boundaries3d,
    ) -> Result<(), MannReservoirError3> {
        if grid != &self.grid {
            return Err(MannReservoirError3::Input(
                "original prepared Mann grid required",
            ));
        }
        for &side in self.sides() {
            for axis in 0..3 {
                if axis != side.axis() && boundaries.axis(axis) == AxisBoundaries::Periodic {
                    return Err(MannReservoirError3::Input(
                        "Mann cross-sections require nonperiodic transverse axes",
                    ));
                }
            }
        }
        Ok(())
    }

    // Share the physical point, exact spatial filter and vector rotation with
    // held traces and mutable stage-frame consumers, without self-borrows or
    // reconstructing inlet geometry for every face.
    fn fluctuation_at(
        &self,
        side: BoundarySide,
        face: usize,
        displacement_m: [f64; 3],
        sample: impl FnOnce([f64; 3], &MannBoxFilter) -> Result<[f64; 3], fs_atmo::Refusal>,
    ) -> Result<[f64; 3], MannReservoirError3> {
        let filter = self
            .spatial_filter(side)
            .ok_or(MannReservoirError3::Input("unprepared Mann inlet side"))?;
        let point = self
            .grid
            .face_center(side, face)
            .ok_or(MannReservoirError3::Input(
                "Mann inlet face outside prepared grid",
            ))?;
        let field_point = to_field(
            self.config.field_to_grid,
            self.config.grid_origin_in_field_m,
            point,
            displacement_m,
        );
        let velocity = to_grid(self.config.field_to_grid, sample(field_point, filter)?);
        if !velocity.iter().all(|x| x.is_finite()) {
            return Err(MannReservoirError3::Input(
                "rotated Mann velocity is not representable",
            ));
        }
        Ok(velocity)
    }
}

impl MannReservoirTrace3<'_> {
    pub const fn report(&self) -> MannReservoirReport3 {
        self.report
    }

    /// Read the same retained realization. No random draw, source advance or
    /// numerical allocation occurs on successful sampling.
    pub fn fluctuation(
        &self,
        side: BoundarySide,
        face: usize,
    ) -> Result<[f64; 3], MannReservoirError3> {
        self.inlet
            .fluctuation_at(side, face, self.displacement_m, |point, filter| {
                let sample = match self.source {
                    MannReservoirSource3::Frozen { source, .. } => {
                        source.sample_filtered(point, filter)
                    }
                    MannReservoirSource3::Evolving(source) => source.sample_filtered(point, filter),
                    MannReservoirSource3::BridgeFrame(frame) => {
                        frame.sample_filtered(point, filter)
                    }
                }?;
                Ok(sample.u)
            })
    }

    pub fn frozen_state(
        &self,
        side: BoundarySide,
        face: usize,
        mut baseline: ReservoirFaceState,
    ) -> Result<ReservoirFaceState, MannReservoirError3> {
        let u = self.fluctuation(side, face)?;
        for i in 0..3 {
            baseline.velocity[i] += u[i];
        }
        Ok(baseline)
    }
    pub fn mixture_state(
        &self,
        side: BoundarySide,
        face: usize,
        mut baseline: MixtureReservoirFaceState,
    ) -> Result<MixtureReservoirFaceState, MannReservoirError3> {
        let u = self.fluctuation(side, face)?;
        for i in 0..3 {
            baseline.velocity[i] += u[i];
        }
        Ok(baseline)
    }
}

#[cfg(test)]
mod tests;
