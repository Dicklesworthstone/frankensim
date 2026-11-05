//! Full-tensor Mann observations through the native D3Q19 velocity boundary.
//!
//! The existing atmosphere owns the spectrum and stochastic evolution. This
//! adapter prepares only a physical frame and a point/face observation, then
//! publishes through the original atomic velocity-profile owner. The same
//! profile is consumed by ordinary, pooled and SGS stepping. No RMS restoration,
//! periodicization, viscous wall layer or new boundary reconstruction is added.
//! The spatial model is Mann (1998), doi:10.1016/S0266-8920(97)00036-2.

use super::{
    BoundaryGrid3, BoundaryProfileError3, Face3, FaceBoundary3, InletFlux3, VelocityProfileReport3,
    planar_face_cells,
};
use fs_atmo::mann_field::filter::MannBoxFilter;
use fs_atmo::mann_field::fixed_step::{FixedStepMannField, MannStepAverage};
use fs_atmo::mann_field::{HomogeneousMannField, MannSpectralModel};
use fs_exec::CancelGate;
use fs_math::det;

/// Explicit SI/lattice and field/grid mapping; no mean or length is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryMannConfig3 {
    pub face: Face3,
    pub dx_m: f64,
    pub dt_s: f64,
    /// Mean velocity in grid Cartesian axes [m/s].
    pub mean_velocity_m_s: [f64; 3],
    /// Field-frame position of lattice node [0, 0, 0], not a cell corner [m].
    pub grid_origin_in_field_m: [f64; 3],
    /// Proper orthonormal rotation: field velocity -> grid velocity.
    pub field_to_grid: [[f64; 3]; 3],
    /// Maximum active face nodes times represented wavevectors per apply.
    pub maximum_mode_samples: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryMannObservation3 {
    Point,
    /// Exact dx-by-dx spatial average centered on each boundary node.
    /// This filters velocity, not nonlinear momentum or energy fluxes.
    FaceAverage,
}

/// Borrow an existing accepted source; applying never advances or redraws it.
#[derive(Clone, Copy, Debug)]
pub enum BoundaryMannSource3<'a> {
    /// Translate a retained snapshot under Taylor's frozen-field hypothesis.
    /// At reference_time_s the snapshot occupies its original coordinates.
    Frozen {
        source: &'a HomogeneousMannField,
        time_s: f64,
        reference_time_s: f64,
        /// Convection velocity in grid Cartesian axes [m/s].
        convection_velocity_m_s: [f64; 3],
    },
    /// Current accepted OU endpoint, including the source's own advection.
    Evolving(&'a FixedStepMannField),
    /// Exact stationary-observer mean of the source's last accepted interval.
    /// The interval must end at source.step() and use the declared CFD dt.
    StepAverage {
        source: &'a FixedStepMannField,
        average: &'a MannStepAverage,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub enum BoundaryMannError3 {
    Input(&'static str),
    SpatialNyquist {
        phase_bound_rad: f64,
    },
    /// Advection only: OU innovations are not temporally bandlimited.
    TemporalNyquist {
        phase_bound_rad: f64,
    },
    Work {
        required: usize,
        maximum: usize,
    },
    Atmosphere(fs_atmo::Refusal),
    Boundary(BoundaryProfileError3),
}

impl core::fmt::Display for BoundaryMannError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-BOUNDARY-MANN3: {self:?}")
    }
}
impl std::error::Error for BoundaryMannError3 {}
impl From<fs_atmo::Refusal> for BoundaryMannError3 {
    fn from(error: fs_atmo::Refusal) -> Self {
        Self::Atmosphere(error)
    }
}
impl From<BoundaryProfileError3> for BoundaryMannError3 {
    fn from(error: BoundaryProfileError3) -> Self {
        Self::Boundary(error)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryMannReport3 {
    pub profile: VelocityProfileReport3,
    /// Equal endpoints for an instantaneous trace; otherwise the exact mean's interval.
    pub sampled_interval_s: [f64; 2],
    pub mode_samples: usize,
    pub spatial_phase_bound_rad: f64,
    pub advective_phase_bound_rad: f64,
}

/// Reusable observation bound to the complete immutable Mann model.
/// Point sampling stores no numerical gains; face sampling stores one gain
/// per wavevector. The existing boundary owner retains its two plane buffers.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryMannInflow3 {
    config: BoundaryMannConfig3,
    filter: MannBoxFilter,
    mode_count: usize,
    velocity_scale: f64,
    spatial_phase_bound_rad: f64,
}

fn poll(gate: &CancelGate) -> Result<(), BoundaryMannError3> {
    if gate.is_requested() {
        Err(BoundaryProfileError3::Cancelled.into())
    } else {
        Ok(())
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

fn speed(v: [f64; 3]) -> f64 {
    let scale = v.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if scale == 0.0 {
        0.0
    } else {
        let normalized = v.map(|x| x / scale);
        scale * det::sqrt(dot(normalized, normalized))
    }
}

impl BoundaryMannInflow3 {
    /// Prepare a finite-band observation without copying modes or source state.
    /// The byte bound covers only the original MannBoxFilter's numerical gains.
    pub fn new(
        model: &MannSpectralModel,
        config: BoundaryMannConfig3,
        observation: BoundaryMannObservation3,
        maximum_filter_bytes: usize,
        gate: &CancelGate,
    ) -> Result<Self, BoundaryMannError3> {
        poll(gate)?;
        let velocity_scale = config.dt_s / config.dx_m;
        if !config.dx_m.is_finite()
            || config.dx_m <= 0.0
            || !config.dt_s.is_finite()
            || config.dt_s <= 0.0
            || !velocity_scale.is_finite()
            || velocity_scale <= 0.0
            || config.maximum_mode_samples == 0
            || !config
                .mean_velocity_m_s
                .iter()
                .chain(&config.grid_origin_in_field_m)
                .chain(config.field_to_grid.iter().flatten())
                .all(|x| x.is_finite())
        {
            return Err(BoundaryMannError3::Input(
                "invalid physical scale, frame, mean or work bound",
            ));
        }
        let rotation = config.field_to_grid;
        for i in 0..3 {
            for j in 0..3 {
                let error = dot(rotation[i], rotation[j]) - if i == j { 1.0 } else { 0.0 };
                if !error.is_finite() || error.abs() > 256.0 * f64::EPSILON {
                    return Err(BoundaryMannError3::Input(
                        "field-to-grid rotation must be orthonormal",
                    ));
                }
            }
        }
        let determinant = rotation[0][0]
            * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
            - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
            + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0]);
        if (determinant - 1.0).abs() > 512.0 * f64::EPSILON {
            return Err(BoundaryMannError3::Input(
                "field-to-grid rotation must preserve handedness",
            ));
        }
        let spatial_phase_bound_rad =
            model.quadrature().max_kl / model.parameters().length_m * config.dx_m;
        if !spatial_phase_bound_rad.is_finite() || spatial_phase_bound_rad > core::f64::consts::PI {
            return Err(BoundaryMannError3::SpatialNyquist {
                phase_bound_rad: spatial_phase_bound_rad,
            });
        }
        let half_edges = core::array::from_fn(|axis| {
            if observation == BoundaryMannObservation3::Point || axis == config.face.axis() {
                [0.0; 3]
            } else {
                rotation[axis].map(|v| v * (0.5 * config.dx_m))
            }
        });
        let filter = MannBoxFilter::prepare(model, half_edges, maximum_filter_bytes, || {
            gate.is_requested()
        })?;
        Ok(Self {
            config,
            filter,
            mode_count: model.mode_count(),
            velocity_scale,
            spatial_phase_bound_rad,
        })
    }

    pub const fn config(&self) -> BoundaryMannConfig3 {
        self.config
    }

    /// Original point/face observation, also usable by the source's covariance API.
    pub const fn spatial_filter(&self) -> &MannBoxFilter {
        &self.filter
    }

    /// Publish a resolved trace at the caller-selected source endpoint/interval.
    /// No source, fluid clock or populations advance. A subsequent refused fluid
    /// step retains this accepted control, exactly as set_velocity_profile does.
    ///
    /// The arbitrary Mann wavevectors are not periodic on a CFD cross-section;
    /// transverse periodic faces refuse. Transverse walls retain the original
    /// wall/open ownership and mean-normal control's mixed-seam refusal.
    /// All sampling uses on-site lattice nodes, not NSF exterior ghost centers.
    /// The interval mean averages velocity only, not nonlinear fluxes, and does
    /// not claim a higher-order time-dependent boundary reconstruction.
    pub fn apply(
        &self,
        grid: &mut BoundaryGrid3,
        source: BoundaryMannSource3<'_>,
        flux: InletFlux3,
        gate: &CancelGate,
    ) -> Result<BoundaryMannReport3, BoundaryMannError3> {
        poll(gate)?;
        let (interval, displacement, advection) = match source {
            BoundaryMannSource3::Frozen {
                time_s,
                reference_time_s,
                convection_velocity_m_s,
                ..
            } => {
                let elapsed = time_s - reference_time_s;
                if !time_s.is_finite()
                    || !reference_time_s.is_finite()
                    || !elapsed.is_finite()
                    || !convection_velocity_m_s.iter().all(|v| v.is_finite())
                {
                    return Err(BoundaryMannError3::Input(
                        "invalid frozen time or convection velocity",
                    ));
                }
                (
                    [time_s; 2],
                    convection_velocity_m_s.map(|v| -v * elapsed),
                    convection_velocity_m_s,
                )
            }
            BoundaryMannSource3::Evolving(field) => {
                if field.dt_s() != self.config.dt_s {
                    return Err(BoundaryMannError3::Input(
                        "evolving source must use the declared CFD dt",
                    ));
                }
                ([field.time_s(); 2], [0.0; 3], field.advection_mps())
            }
            BoundaryMannSource3::StepAverage {
                source: field,
                average,
            } => {
                average.validate_source(field)?;
                if field.dt_s() != self.config.dt_s || average.end_step() != Some(field.step()) {
                    return Err(BoundaryMannError3::Input(
                        "mean must end at the source's accepted CFD endpoint",
                    ));
                }
                (
                    average
                        .interval()
                        .ok_or(BoundaryMannError3::Input("no accepted interval mean"))?,
                    [0.0; 3],
                    field.advection_mps(),
                )
            }
        };
        let next_time = interval[1] + self.config.dt_s;
        if !interval.iter().chain(&displacement).all(|v| v.is_finite())
            || !next_time.is_finite()
            || next_time <= interval[1]
        {
            return Err(BoundaryMannError3::Input(
                "source time or displacement is not representable",
            ));
        }
        let advective_phase_bound_rad =
            self.spatial_phase_bound_rad * (speed(advection) * self.velocity_scale);
        if !advective_phase_bound_rad.is_finite()
            || advective_phase_bound_rad > core::f64::consts::PI
        {
            return Err(BoundaryMannError3::TemporalNyquist {
                phase_bound_rad: advective_phase_bound_rad,
            });
        }
        let face = self.config.face;
        if !matches!(grid.boundaries().face(face), FaceBoundary3::Velocity { .. }) {
            return Err(BoundaryProfileError3::NotVelocityFace.into());
        }
        for axis in 0..3 {
            if axis != face.axis()
                && matches!(
                    grid.boundaries().face(Face3::ALL[2 * axis]),
                    FaceBoundary3::Periodic
                )
            {
                return Err(BoundaryMannError3::Input(
                    "Mann wavevectors require nonperiodic transverse faces",
                ));
            }
        }
        let mut cells = 0_usize;
        for (index, p) in planar_face_cells(grid.dimensions(), face).enumerate() {
            if index.is_multiple_of(64) {
                poll(gate)?;
            }
            if !grid.is_solid(p[0], p[1], p[2]) && grid.open_link_mask(p[0], p[1], p[2]) != 0 {
                cells += 1;
            }
        }
        let mode_samples = cells
            .checked_mul(self.mode_count)
            .ok_or(BoundaryMannError3::Input("inlet sampling work overflows"))?;
        if mode_samples > self.config.maximum_mode_samples {
            return Err(BoundaryMannError3::Work {
                required: mode_samples,
                maximum: self.config.maximum_mode_samples,
            });
        }
        let mut source_error = None;
        let result = grid.set_velocity_profile(face, flux, gate, |node| {
            let point_grid: [f64; 3] =
                core::array::from_fn(|a| node[a] as f64 * self.config.dx_m + displacement[a]);
            let point_field = core::array::from_fn(|a| {
                self.config.grid_origin_in_field_m[a]
                    + (0..3)
                        .map(|i| self.config.field_to_grid[i][a] * point_grid[i])
                        .sum::<f64>()
            });
            let sampled = match source {
                BoundaryMannSource3::Frozen { source, .. } => {
                    source.sample_filtered(point_field, &self.filter)
                }
                BoundaryMannSource3::Evolving(source) => {
                    source.sample_filtered(point_field, &self.filter)
                }
                BoundaryMannSource3::StepAverage { average, .. } => {
                    average.sample_filtered(point_field, &self.filter)
                }
            };
            match sampled {
                Ok(field) => Ok(core::array::from_fn(|a| {
                    (self.config.mean_velocity_m_s[a] + dot(self.config.field_to_grid[a], field.u))
                        * self.velocity_scale
                })),
                Err(error) => {
                    source_error = Some(error);
                    Err(BoundaryProfileError3::Source(
                        "Mann observation refused".into(),
                    ))
                }
            }
        });
        let profile = match result {
            Ok(profile) => profile,
            Err(error) => {
                return Err(match source_error {
                    Some(error) => error.into(),
                    None => error.into(),
                });
            }
        };
        Ok(BoundaryMannReport3 {
            profile,
            sampled_interval_s: interval,
            mode_samples,
            spatial_phase_bound_rad: self.spatial_phase_bound_rad,
            advective_phase_bound_rad,
        })
    }
}
