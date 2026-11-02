//! Frequency-independent geometry for repeated evaluations of one modal duct.

use super::*;

pub(super) struct PreparedJunction {
    pub(super) expansion: bool,
    pub(super) projection: Vec<C64>,
    pub(super) transpose: Vec<C64>,
}

/// Original staircase geometry and its radial basis, prepared once for a sweep.
///
/// Preparation uses exactly the one-frequency API's station discretization,
/// J1 root search and junction projection arithmetic. The borrowed source cannot
/// be changed while this preparation is used, and the retained mode count and
/// geometry have no mutating or replacement API. A different geometry or mode
/// rung requires a new preparation.
///
/// No gas state, frequency, loss model or mouth impedance is cached. Every
/// evaluation recomputes propagation and the physical load from its supplied
/// material state, so pressure, temperature and composition substitutions remain
/// live. This is the same staircase model, not a different approximation.
pub struct PreparedModalDuct<'a> {
    source: &'a Duct,
    pub(super) n_modes: usize,
    pub(super) stations: Vec<ModalStation>,
    pub(super) gammas: Vec<f64>,
    pub(super) junctions: Vec<Option<PreparedJunction>>,
}

impl<'a> PreparedModalDuct<'a> {
    /// Prepare an admitted, hole-free cylindrical/conical staircase.
    ///
    /// # Errors
    /// The same geometry, retained-mode and unsupported-segment refusals as
    /// [`mm_input_impedance`].
    pub fn new(duct: &'a Duct, n_modes: usize, extra_slices: usize) -> Result<Self, DuctError> {
        Self::new_controlled(duct, n_modes, extra_slices, || false)
    }

    /// [`Self::new`] with cancellation before preparation and at each junction.
    /// A cancellation reports the number of fully prepared station junctions.
    ///
    /// # Errors
    /// The same admission failures as [`Self::new`], or [`DuctError::Cancelled`].
    pub fn new_controlled(
        duct: &'a Duct,
        n_modes: usize,
        extra_slices: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        if cancelled() {
            return Err(DuctError::Cancelled { completed_steps: 0 });
        }
        if n_modes == 0 {
            return Err(DuctError::BadParameter {
                what: "at least one mode (the plane wave)",
            });
        }
        if n_modes > MAX_MODES {
            return Err(DuctError::BadParameter {
                what: "mode count exceeds MAX_MODES",
            });
        }
        for segment in &duct.segments {
            segment.validate()?;
        }
        let stations = stations_of(duct, extra_slices)?;
        let gammas = j1_roots(n_modes - 1)?;
        let mut junctions = Vec::with_capacity(stations.len() - 1);
        for pair in stations.windows(2) {
            if cancelled() {
                return Err(DuctError::Cancelled { completed_steps: junctions.len() });
            }
            let upstream = pair[0].radius;
            let downstream = pair[1].radius;
            let junction = if (upstream - downstream).abs() > 1e-12 * upstream.max(1e-300) {
                let expansion = upstream < downstream;
                let projection = if expansion {
                    junction_projection(upstream, downstream, &gammas, n_modes)
                } else {
                    junction_projection(downstream, upstream, &gammas, n_modes)
                };
                Some(PreparedJunction {
                    expansion,
                    transpose: transpose(&projection, n_modes),
                    projection,
                })
            } else {
                None
            };
            junctions.push(junction);
        }
        if cancelled() {
            return Err(DuctError::Cancelled { completed_steps: junctions.len() });
        }
        Ok(Self { source: duct, n_modes, stations, gammas, junctions })
    }

    /// The immutable source geometry that owns this preparation.
    #[must_use]
    pub const fn source(&self) -> &'a Duct { self.source }

    /// Retained radial mode count, including the plane mode.
    #[must_use]
    pub const fn n_modes(&self) -> usize { self.n_modes }

    /// Actual last discrete station's mouth radius [m]. For a sliced cone this
    /// is the last midpoint radius, where this staircase attaches its load.
    #[must_use]
    pub fn mouth_radius_m(&self) -> f64 {
        self.stations.last().expect("admitted nonempty modal duct").radius
    }

    /// [`mm_input_impedance`] using this exact prepared geometry.
    pub fn input_impedance(
        &self, gas: &GasState, omega: f64, loss: LossModel, termination: Termination,
    ) -> Result<ModalResponse, DuctError> {
        mm_evaluate_prepared(self, gas, omega, loss, MmLoad::Analytic(termination), false)
            .map(|response| response.input)
    }

    /// [`mm_input_impedance_tabulated`] with the original table admission.
    pub fn input_impedance_tabulated(
        &self, gas: &GasState, omega: f64, loss: LossModel, table: &crate::TabulatedLoad,
    ) -> Result<ModalResponse, DuctError> {
        mm_evaluate_prepared(self, gas, omega, loss, MmLoad::Tabulated(table), false)
            .map(|response| response.input)
    }

    /// [`mm_input_impedance_load`] with a finite passive plane-mouth load.
    pub fn input_impedance_load(
        &self, gas: &GasState, omega: f64, loss: LossModel, impedance: C64,
    ) -> Result<ModalResponse, DuctError> {
        if !impedance.re.is_finite() || !impedance.im.is_finite() || impedance.re < 0.0 {
            return Err(DuctError::BadParameter {
                what: "modal mouth impedance must be finite with nonnegative resistance",
            });
        }
        mm_evaluate_prepared(self, gas, omega, loss, MmLoad::Impedance(impedance), false)
            .map(|response| response.input)
    }

    /// [`mm_input_impedance_flanged`] with live interior and exterior states.
    pub fn input_impedance_flanged(
        &self, gas: &GasState, exterior: &GasState, omega: f64, loss: LossModel,
    ) -> Result<ModalResponse, DuctError> {
        mm_evaluate_prepared(self, gas, omega, loss, MmLoad::Flanged(exterior), false)
            .map(|response| response.input)
    }

    /// [`mm_input_impedance_matrix`] on this exact prepared geometry and basis.
    pub fn input_impedance_matrix(
        &self, gas: &GasState, omega: f64, loss: LossModel, load: &ModalMouthImpedance,
    ) -> Result<ModalResponse, DuctError> {
        mm_evaluate_prepared(self, gas, omega, loss, MmLoad::Matrix(load), false)
            .map(|response| response.input)
    }

    /// Native input reflection before the input-impedance conversion.
    ///
    /// The same propagation and load equations as [`Self::input_impedance`]
    /// apply. Flow nodes remain finite reflected waves; no impedance inverse
    /// is required. Use [`ModalReflectionResponse::plane_reflection`] for a
    /// plane source or [`ModalReflectionResponse::at_references`] for all ports.
    ///
    /// # Errors
    /// The same wave/load refusals as [`Self::input_impedance`], without its
    /// final input-impedance-pole refusal.
    pub fn input_reflection(
        &self, gas: &GasState, omega: f64, loss: LossModel, termination: Termination,
    ) -> Result<ModalReflectionResponse, DuctError> {
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Analytic(termination), false)
            .map(|response| response.input)
    }

    /// [`Self::input_reflection`] with the same tabulated mouth and matched
    /// higher-mode closure as [`Self::input_impedance_tabulated`].
    ///
    /// # Errors
    /// The table's support and wave admission refusals apply at this frequency.
    pub fn input_reflection_tabulated(
        &self, gas: &GasState, omega: f64, loss: LossModel, table: &crate::TabulatedLoad,
    ) -> Result<ModalReflectionResponse, DuctError> {
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Tabulated(table), false)
            .map(|response| response.input)
    }

    /// [`Self::input_reflection`] with a finite passive plane-mouth load,
    /// including an independently declared static load for quasi-static use.
    /// Higher mouth modes retain the matched closure.
    ///
    /// # Errors
    /// Invalid waves or a nonfinite/negative-resistance mouth impedance.
    pub fn input_reflection_load(
        &self, gas: &GasState, omega: f64, loss: LossModel, impedance: C64,
    ) -> Result<ModalReflectionResponse, DuctError> {
        if !impedance.re.is_finite() || !impedance.im.is_finite() || impedance.re < 0.0 {
            return Err(DuctError::BadParameter {
                what: "modal mouth impedance must be finite with nonnegative resistance",
            });
        }
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Impedance(impedance), false)
            .map(|response| response.input)
    }

    /// [`Self::input_reflection`] with the complete flanged radiation matrix.
    /// The interior/exterior material states and all mutual loads stay live.
    ///
    /// # Errors
    /// The same physical admissions as [`Self::input_impedance_flanged`].
    pub fn input_reflection_flanged(
        &self, gas: &GasState, exterior: &GasState, omega: f64, loss: LossModel,
    ) -> Result<ModalReflectionResponse, DuctError> {
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Flanged(exterior), false)
            .map(|response| response.input)
    }

    /// Native reflected waves with the same full physical matrix as
    /// [`Self::input_impedance_matrix`]. An inlet impedance pole alone does not
    /// make these wave coordinates singular.
    ///
    /// # Errors
    /// Matrix declaration or native wave refusals, without final conversion to
    /// input impedance.
    pub fn input_reflection_matrix(
        &self, gas: &GasState, omega: f64, loss: LossModel, load: &ModalMouthImpedance,
    ) -> Result<ModalReflectionResponse, DuctError> {
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Matrix(load), false)
            .map(|response| response.input)
    }

    /// Native input reflection and mouth-flow observation from one loaded
    /// solve. Geometry, roots and junction projections are reused across bins.
    ///
    /// # Errors
    /// Original wave/load admissions; no final input-impedance inverse.
    pub fn mouth_waves(
        &self, gas: &GasState, omega: f64, loss: LossModel, termination: Termination,
    ) -> Result<outlet::MouthWaveResponse, DuctError> {
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Analytic(termination), true)
            .and_then(outlet::wave_response)
    }

    /// [`Self::mouth_waves`] with a strict table bound to the original mouth
    /// plane, radius and gas. Higher mouth modes retain the matched closure.
    ///
    /// # Errors
    /// Table reference/support mismatch or the original wave admissions.
    pub fn mouth_waves_tabulated(
        &self, gas: &GasState, omega: f64, loss: LossModel, table: &crate::TabulatedLoad,
    ) -> Result<outlet::MouthWaveResponse, DuctError> {
        let radius = self.source.segments.last().ok_or(DuctError::EmptyDuct)?.outlet_radius();
        table.validate_mouth(self.source.segments.iter().map(Segment::axial_length).sum(), radius, gas)?;
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Tabulated(table), true)
            .and_then(outlet::wave_response)
    }

    /// [`Self::mouth_waves`] with a finite passive plane-mouth impedance,
    /// including an independently declared static load at a quasi-static probe.
    ///
    /// # Errors
    /// Invalid waves or a nonfinite/negative-resistance mouth impedance.
    pub fn mouth_waves_load(
        &self, gas: &GasState, omega: f64, loss: LossModel, impedance: C64,
    ) -> Result<outlet::MouthWaveResponse, DuctError> {
        if !impedance.re.is_finite() || !impedance.im.is_finite() || impedance.re < 0.0 {
            return Err(DuctError::BadParameter {
                what: "modal mouth impedance must be finite with nonnegative resistance",
            });
        }
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Impedance(impedance), true)
            .and_then(outlet::wave_response)
    }

    /// [`Self::mouth_waves`] with the full coupled flanged radiation matrix
    /// and original interior/exterior gases at their common static pressure.
    ///
    /// # Errors
    /// The same physical admissions as [`Self::input_reflection_flanged`].
    pub fn mouth_waves_flanged(
        &self, gas: &GasState, exterior: &GasState, omega: f64, loss: LossModel,
    ) -> Result<outlet::MouthWaveResponse, DuctError> {
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Flanged(exterior), true)
            .and_then(outlet::wave_response)
    }

    /// [`Self::mouth_waves`] with an original complete physical mouth matrix,
    /// including geometry-derived exterior loads on every retained mode.
    ///
    /// # Errors
    /// The same physical admissions as [`Self::input_reflection_matrix`].
    pub fn mouth_waves_matrix(
        &self, gas: &GasState, omega: f64, loss: LossModel, load: &ModalMouthImpedance,
    ) -> Result<outlet::MouthWaveResponse, DuctError> {
        mm_reflection_prepared(self, gas, omega, loss, MmLoad::Matrix(load), true)
            .and_then(outlet::wave_response)
    }
}

#[cfg(test)]
mod tests;
