//! Source-qualified longitudinal mechanical dispersion on actual hydraulic faces.
//!
//! The local seepage speed is |q|/theta and D_mech=alpha_L |q|/theta.
//! Thus the bulk coefficient consumed by the original finite volumes is
//! theta*D_mech=alpha_L |q|, separate from theta*D_molecular/tau.
//! This is the longitudinal specialization of USGS MODFLOW 6 GWT, section 4,
//! equation 4-19: <https://pubs.usgs.gov/tm/06/a61/tm6a61.pdf>.
//! Calibration scale remains an explicit physical source coordinate, not a
//! numerical cell width. The reference supplies no material-specific alpha.

use super::langmuir::LangmuirSoluteError;
use super::transport::{
    self, DiluteSoluteControls, DiluteSoluteError, finite, integrated_flux, nonnegative, poll,
    positive,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_qty::{Dims, Length, Qty, QuantitySpec};
use std::sync::Arc;

pub const LONGITUDINAL_DISPERSION_LAW: &str = "porous-liquid-longitudinal-mechanical-dispersion";
pub const DISPERSIVITY_PARAMETER: &str = "longitudinal-dispersivity";
pub const PORE_SPEED_AXIS: &str = "pore-speed";
pub const LIQUID_SATURATION_AXIS: &str = "liquid-saturation";
pub const POROSITY_AXIS: &str = "total-porosity";
pub const TRANSPORT_SCALE_AXIS: &str = "dispersion-transport-scale";
type PoreSpeed = Qty<1, 0, -1, 0, 0, 0>;

/// Original constant dispersivity and declared applicability of this local law.
#[derive(Clone, Debug, PartialEq)]
pub struct LongitudinalDispersivity {
    source: Arc<ConstitutiveModelCard>,
    alpha_m: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SoluteDispersionError {
    Source(MatDbError),
    Invalid(&'static str),
    OutsideValidity {
        cell: usize,
        side: usize,
        axis: String,
        value: f64,
        minimum: f64,
        maximum: f64,
    },
    SourceWork {
        required_points: usize,
        maximum_points: usize,
    },
    Transport(DiluteSoluteError),
    Langmuir(LangmuirSoluteError),
}
impl From<MatDbError> for SoluteDispersionError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e)
    }
}
impl From<DiluteSoluteError> for SoluteDispersionError {
    fn from(e: DiluteSoluteError) -> Self {
        Self::Transport(e)
    }
}
impl core::fmt::Display for SoluteDispersionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "porous mechanical dispersion: {self:?}")
    }
}
impl std::error::Error for SoluteDispersionError {}
type Result<T> = core::result::Result<T, SoluteDispersionError>;

fn owned_axes() -> [(&'static str, QuantitySpec); 4] {
    [
        (PORE_SPEED_AXIS, QuantitySpec::dimensional(PoreSpeed::DIMS)),
        (
            LIQUID_SATURATION_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
        ),
        (POROSITY_AXIS, QuantitySpec::dimensional(Dims::NONE)),
        (
            TRANSPORT_SCALE_AXIS,
            QuantitySpec::dimensional(Length::DIMS),
        ),
    ]
}

impl LongitudinalDispersivity {
    /// Admit a stateless version-one, single-parameter original source. Alpha
    /// has length units; it is neither molecular D nor a tortuosity divisor.
    /// Source support must cover local speed, saturation, porosity and physical
    /// transport scale, with any extra conditions explicitly typed and bounded.
    pub fn from_card(source: ConstitutiveModelCard) -> Result<Self> {
        source.validate()?;
        if source.law.0 != LONGITUDINAL_DISPERSION_LAW
            || source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.parameters.len() != 1
            || source.sources.is_empty()
            || source.provenance.source.trim().is_empty()
            || source.provenance.license.trim().is_empty()
        {
            return Err(SoluteDispersionError::Invalid(
                "original sourced stateless longitudinal dispersion card required",
            ));
        }
        let alpha =
            source
                .parameters
                .get(DISPERSIVITY_PARAMETER)
                .ok_or(SoluteDispersionError::Invalid(
                    "missing longitudinal dispersivity",
                ))?;
        if alpha.dims != Length::DIMS || !nonnegative(alpha.value) {
            return Err(SoluteDispersionError::Invalid(
                "finite nonnegative longitudinal dispersivity in metres required",
            ));
        }
        for (axis, &(low, high)) in source.validity.bounds() {
            if !low.is_finite()
                || !high.is_finite()
                || low > high
                || !source.validity.axis_quantities().contains_key(axis)
            {
                return Err(SoluteDispersionError::Invalid(
                    "every dispersion source condition needs finite ordered typed support",
                ));
            }
        }
        for (axis, quantity) in owned_axes() {
            let (low, high) = source
                .validity
                .bound(axis)
                .ok_or(SoluteDispersionError::Invalid(
                    "missing owned dispersion source axis",
                ))?;
            if source.validity.axis_quantities().get(axis) != Some(&quantity)
                || low < 0.0
                || ([POROSITY_AXIS, TRANSPORT_SCALE_AXIS].contains(&axis) && low == 0.0)
                || ([LIQUID_SATURATION_AXIS, POROSITY_AXIS].contains(&axis) && high > 1.0)
            {
                return Err(SoluteDispersionError::Invalid(
                    "dispersion source axis has invalid bounds or quantity",
                ));
            }
        }
        let alpha_m = alpha.value;
        Ok(Self {
            source: Arc::new(source),
            alpha_m,
        })
    }
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    pub fn dispersivity(&self) -> Length {
        Length::new(self.alpha_m)
    }

    fn admit(
        &self,
        cell: usize,
        side: usize,
        pore_speed: f64,
        geometry: DispersionCellState,
        scale: f64,
        context: &QueryPoint,
    ) -> Result<()> {
        // Borrow the existing source and explicit context. No source-card or
        // query cloning is needed for these four owner-supplied coordinates.
        for (axis, &(minimum, maximum)) in self.source.validity.bounds() {
            let value = match axis.as_str() {
                PORE_SPEED_AXIS => pore_speed,
                LIQUID_SATURATION_AXIS => geometry.liquid_saturation,
                POROSITY_AXIS => geometry.porosity,
                TRANSPORT_SCALE_AXIS => scale,
                _ => {
                    if context.axis_quantities().get(axis)
                        != self.source.validity.axis_quantities().get(axis)
                    {
                        return Err(SoluteDispersionError::Invalid(
                            "missing or mistyped additional original dispersion condition",
                        ));
                    }
                    *context
                        .axes()
                        .get(axis)
                        .ok_or(SoluteDispersionError::Invalid(
                            "missing additional original dispersion condition",
                        ))?
                }
            };
            if !value.is_finite() || value < minimum || value > maximum {
                return Err(SoluteDispersionError::OutsideValidity {
                    cell,
                    side,
                    axis: axis.clone(),
                    value,
                    minimum,
                    maximum,
                });
            }
        }
        Ok(())
    }
}

/// One immutable calibration binding per original cell, shared by all passive
/// constituent slots. Its source must apply to those constituents and context.
#[derive(Clone, Copy, Debug)]
pub struct SoluteDispersionCell<'a> {
    pub model: &'a LongitudinalDispersivity,
    /// Physical scale covered by the calibration, unchanged by mesh refinement.
    pub transport_scale: Length,
    /// Additional original coordinates only. Their consistency with any prior
    /// hydraulic temperature/chemistry is caller-owned; that history is not
    /// retained by the hydraulic candidate and is not inferred here.
    pub conditions: &'a QueryPoint,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SoluteDispersionControls {
    /// Both faces of every cell, including zero-speed/zero-alpha points.
    pub maximum_source_points: usize,
    /// Preparation payload cap. The subsequent transport cap ALSO includes
    /// this retained payload; it is not an extra unbudgeted allocation.
    pub maximum_storage_bytes: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DispersionCellState {
    pub porosity: f64,
    pub liquid_saturation: f64,
    /// Actual accepted water capacity divided by physical cell width.
    pub liquid_volume_fraction: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SoluteDispersionReport {
    pub interval_s: [f64; 2],
    pub cells: Vec<DispersionCellState>,
    /// Left/right face speed magnitudes in each cell's own liquid volume.
    pub pore_speed_m_s: Vec<[f64; 2]>,
    /// Bulk half-cell contribution alpha_L*|q_face| [m²/s], before molecular
    /// addition and series assembly. These coefficients are not extra fluxes.
    pub mechanical_bulk_m2_s: Vec<[f64; 2]>,
    pub source_points: usize,
    pub preparation_storage_bytes: usize,
}

/// A native conservative transport result with original dispersion evidence.
/// Its nested `diffusive_face_amount_mol_m2` holds TOTAL molecular plus
/// mechanical gradient transport. A separate mechanism split is not inferred
/// at an interface where the two adjacent material fractions differ.
#[derive(Clone, Debug)]
pub struct DispersiveSoluteStep<'a, T> {
    pub transport: T,
    pub dispersion: SoluteDispersionReport,
    bindings: &'a [SoluteDispersionCell<'a>],
}
impl<'a, T> DispersiveSoluteStep<'a, T> {
    pub fn models(&self) -> &'a [SoluteDispersionCell<'a>] {
        self.bindings
    }
}

pub(super) struct PreparedDispersion<'a> {
    pub(super) report: SoluteDispersionReport,
    bindings: &'a [SoluteDispersionCell<'a>],
}
impl<'a> PreparedDispersion<'a> {
    /// Reserve retained preparation storage inside the WHOLE transport cap.
    pub(super) fn transport_controls<const N: usize>(
        &self,
        mut controls: DiluteSoluteControls<N>,
    ) -> Result<DiluteSoluteControls<N>> {
        controls.maximum_storage_bytes = controls
            .maximum_storage_bytes
            .checked_sub(self.report.preparation_storage_bytes)
            .ok_or_else(|| {
                SoluteDispersionError::Transport(DiluteSoluteError::Workspace {
                    required_bytes: self.report.preparation_storage_bytes,
                    maximum_bytes: controls.maximum_storage_bytes,
                })
            })?;
        Ok(controls)
    }
    pub(super) fn transport_error(&self, error: DiluteSoluteError) -> SoluteDispersionError {
        SoluteDispersionError::Transport(self.restore_budget(error))
    }
    pub(super) fn langmuir_error(&self, error: LangmuirSoluteError) -> SoluteDispersionError {
        SoluteDispersionError::Langmuir(match error {
            LangmuirSoluteError::Transport(e) => {
                LangmuirSoluteError::Transport(self.restore_budget(e))
            }
            other => other,
        })
    }
    fn restore_budget(&self, error: DiluteSoluteError) -> DiluteSoluteError {
        match error {
            DiluteSoluteError::Workspace {
                required_bytes,
                maximum_bytes,
            } => DiluteSoluteError::Workspace {
                required_bytes: required_bytes
                    .saturating_add(self.report.preparation_storage_bytes),
                maximum_bytes: maximum_bytes.saturating_add(self.report.preparation_storage_bytes),
            },
            other => other,
        }
    }
    pub(super) fn finish<T>(self, transport: T) -> DispersiveSoluteStep<'a, T> {
        DispersiveSoluteStep {
            transport,
            dispersion: self.report,
            bindings: self.bindings,
        }
    }
}

/// Both coordinates and integrated face volumes come from the same admitted
/// original hydraulic interval. The private numerical seam is also the test
/// seam; the public interval adapter owns the exact state/time binding.
pub(super) fn prepare<'a>(
    interval_s: [f64; 2],
    duration_s: f64,
    bindings: &'a [SoluteDispersionCell<'a>],
    controls: SoluteDispersionControls,
    cell: impl Fn(usize) -> DispersionCellState,
    face_volume: impl Fn(usize) -> transport::Result<f64>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedDispersion<'a>> {
    poll(&mut cancelled)?;
    let n = bindings.len();
    if n == 0
        || n > 4096
        || !positive(duration_s)
        || !interval_s.iter().all(|v| v.is_finite())
        || interval_s[1] <= interval_s[0]
    {
        return Err(SoluteDispersionError::Invalid(
            "bounded cells and advancing physical dispersion interval required",
        ));
    }
    if 2 * n > controls.maximum_source_points {
        return Err(SoluteDispersionError::SourceWork {
            required_points: 2 * n,
            maximum_points: controls.maximum_source_points,
        });
    }
    let bytes = 7 * n * core::mem::size_of::<f64>() + 512;
    if bytes > controls.maximum_storage_bytes {
        return Err(DiluteSoluteError::Workspace {
            required_bytes: bytes,
            maximum_bytes: controls.maximum_storage_bytes,
        }
        .into());
    }
    let mut report = SoluteDispersionReport {
        interval_s,
        cells: Vec::with_capacity(n),
        pore_speed_m_s: Vec::with_capacity(n),
        mechanical_bulk_m2_s: Vec::with_capacity(n),
        source_points: 0,
        preparation_storage_bytes: bytes,
    };
    for (i, binding) in bindings.iter().enumerate() {
        poll(&mut cancelled)?;
        let geometry = cell(i);
        if !positive(geometry.porosity)
            || geometry.porosity > 1.0
            || !positive(geometry.liquid_saturation)
            || geometry.liquid_saturation > 1.0
            || !positive(geometry.liquid_volume_fraction)
            || geometry.liquid_volume_fraction > 1.0
            || (geometry.liquid_volume_fraction - geometry.porosity * geometry.liquid_saturation)
                .abs()
                > 16.0 * f64::EPSILON * geometry.liquid_volume_fraction
            || !positive(binding.transport_scale.value())
            || owned_axes()
                .iter()
                .any(|(axis, _)| binding.conditions.axes().contains_key(*axis))
        {
            return Err(SoluteDispersionError::Invalid(
                "strictly wet actual cell coordinates, physical transport scale and unmodified owned axes required",
            ));
        }
        let (mut speed, mut bulk) = ([0.0; 2], [0.0; 2]);
        for side in 0..2 {
            poll(&mut cancelled)?;
            let volume = face_volume(i + side)?;
            let q = finite(volume.abs() / duration_s)?;
            if volume != 0.0 && q == 0.0 {
                return Err(SoluteDispersionError::Invalid(
                    "nonzero hydraulic interval speed is unrepresentable",
                ));
            }
            speed[side] = finite(q / geometry.liquid_volume_fraction)?;
            binding.model.admit(
                i,
                side,
                speed[side],
                geometry,
                binding.transport_scale.value(),
                binding.conditions,
            )?;
            bulk[side] = integrated_flux(binding.model.alpha_m, q)?;
            report.source_points += 1;
        }
        report.cells.push(geometry);
        report.pore_speed_m_s.push(speed);
        report.mechanical_bulk_m2_s.push(bulk);
    }
    poll(&mut cancelled)?;
    Ok(PreparedDispersion { report, bindings })
}

#[cfg(test)]
mod tests;
