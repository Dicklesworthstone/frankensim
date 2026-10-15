//! Prescribed, piecewise-uniform gases connected at TMM pressure/volume-flow ports.
//!
//! Each constant-medium run is evaluated by `fs-duct`, including its cones,
//! wall law and tone-hole junctions. Its downstream load is the next run's input
//! impedance. Pressure and volume flow are continuous at interfaces: no averaged
//! gas, artificial reflection coefficient or additional propagation law is used.
//! This is frozen, no-mean-flow acoustics, not thermal/composition transport.

/// Physical outlet observation of the existing regional transfer.
pub mod outlet;

use super::{DrivingPointError, characteristic_line_response};
use fs_duct::{Duct, DuctError, DuctResponse, HoleState, LossModel, Segment, Termination,
    input_impedance_load, input_impedance_wall};
use fs_material::gas::GasState;
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

/// Cold-path cap, including axial segments and zero-length tone-hole records.
pub const MAX_REGIONAL_SEGMENTS: usize = 64;

#[derive(Debug, Clone)]
struct Region {
    duct: Duct,
    gas: GasState,
}

/// Static regional duct composed from the existing uniform TMM owners.
///
/// Adjacent identical gases coalesce, preserving the exact uniform calculation.
/// Every open/vented tone hole must lie in the SAME constant-medium run so that
/// `fs-duct` still evaluates every mutual-hole correction. Cross-region mutual
/// hole laws are not available and refuse rather than silently disappearing.
#[derive(Debug, Clone)]
pub struct RegionalDuct {
    regions: Vec<Region>,
    gases: Vec<GasState>,
    one_way_time_s: f64,
}

fn invalid(what: &'static str) -> DuctError { DuctError::BadParameter { what } }

impl RegionalDuct {
    /// Bind exactly one resolved state to every authored geometry record.
    ///
    /// States must be coherent and at exactly the same static pressure. Species
    /// and temperature may vary. A tone-hole record uses its assigned gas for
    /// both its compact bore junction and chimney; no exterior medium jump is
    /// inferred. Geometry/radiation/loss validity remains checked by `fs-duct`
    /// at each requested frequency. Construction does not advance sound state.
    ///
    /// # Errors
    /// Empty/oversized or incomplete binding, incoherent gases, unequal pressure,
    /// nonpositive axial lengths, nonaxial ends, or unsupported mutual-hole scope.
    pub fn new(duct: &Duct, gases: &[GasState]) -> Result<Self, DuctError> {
        if duct.segments.is_empty() { return Err(DuctError::EmptyDuct); }
        if duct.segments.len() > MAX_REGIONAL_SEGMENTS || gases.len() != duct.segments.len() {
            return Err(invalid("regional duct requires 1..=64 segments and one gas per segment"));
        }
        if matches!(duct.segments.first(), Some(Segment::ToneHole { .. }))
            || matches!(duct.segments.last(), Some(Segment::ToneHole { .. })) {
            return Err(invalid("regional duct ends must be explicit axial segments"));
        }
        let mut regions: Vec<Region> = Vec::new();
        let mut open_hole_region = None;
        for (segment, gas) in duct.segments.iter().zip(gases) {
            gas.validate_coherence().map_err(|_| invalid("regional gas is not a coherent material state"))?;
            if gas.pressure != gases[0].pressure {
                return Err(invalid("regional duct needs a common static pressure; mean flow is not modeled"));
            }
            if regions.last().is_none_or(|r| r.gas != *gas) {
                regions.push(Region { duct: Duct { segments: Vec::new() }, gas: *gas });
            }
            match *segment {
                Segment::Cylinder { length, .. } | Segment::Cone { length, .. } => {
                    if !length.is_finite() || length <= 0.0 {
                        return Err(invalid("regional axial lengths must be positive and finite"));
                    }
                }
                Segment::ToneHole { state, .. } => {
                    let open = match state {
                        HoleState::Open => true,
                        HoleState::Closed => false,
                        HoleState::Vent(s) => {
                            if !s.is_finite() || !(0.0..=1.0).contains(&s) {
                                return Err(invalid("regional vent fraction must be finite and in [0,1]"));
                            }
                            s > 0.0
                        }
                    };
                    if open {
                        let index = regions.len() - 1;
                        if open_hole_region.is_some_and(|previous| previous != index) {
                            return Err(invalid("open tone holes in different gas regions need an unavailable mutual-hole law"));
                        }
                        open_hole_region = Some(index);
                    }
                }
            }
            regions.last_mut().expect("run just admitted").duct.segments.push(*segment);
        }
        let one_way_time_s = regions.iter().map(|region| {
            let length: f64 = region.duct.segments.iter().map(|segment| match *segment {
                Segment::Cylinder { length, .. } | Segment::Cone { length, .. } => length,
                Segment::ToneHole { .. } => 0.0,
            }).sum();
            length / region.gas.sound_speed
        }).sum::<f64>();
        if !one_way_time_s.is_finite() || one_way_time_s <= 0.0 {
            return Err(invalid("regional travel time must be positive and finite"));
        }
        Ok(Self { regions, gases: gases.to_vec(), one_way_time_s })
    }

    /// Resolved media in authored segment order, including hole records.
    #[must_use]
    pub fn segment_gases(&self) -> &[GasState] { &self.gases }

    /// Actual inlet medium for Bernoulli density and the characteristic port.
    #[must_use]
    pub fn inlet_gas(&self) -> &GasState { &self.gases[0] }

    /// Actual terminal medium used by the shared analytic radiation load.
    #[must_use]
    pub fn outlet_gas(&self) -> &GasState { self.gases.last().expect("nonempty binding") }

    /// Sum of each constant-medium axial length divided by its local sound speed.
    /// This excludes side-branch chimney length, as does the uniform line owner.
    #[must_use]
    pub const fn one_way_time_s(&self) -> f64 { self.one_way_time_s }

    /// Evaluate one frequency with the original TMM loss, wall and load laws.
    ///
    /// # Errors
    /// The existing duct owner refuses geometry, frequency, loss or radiation;
    /// or an intermediate port is nonfinite (for example an exact lossless pole).
    pub fn response(&self, omega: f64, loss: LossModel, termination: Termination,
        wall: Option<&WallPin>) -> Result<DuctResponse, DuctError> {
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid("frequency must be positive and finite"));
        }
        let mut runs = self.regions.iter().rev();
        let terminal = runs.next().expect("nonempty binding");
        let mut response = input_impedance_wall(&terminal.duct, &terminal.gas,
            omega, loss, termination, wall)?;
        for region in runs {
            let upstream = input_impedance_load(&region.duct, &region.gas, omega,
                loss, response.impedance, wall)?;
            // A downstream one-port receives the upstream run's output flow.
            // Keep the REAL mouth's ka and load; this interface is not another
            // radiating termination. Units are pressure and VOLUME flow on both
            // sides, so no density, area or power-wave renormalization is needed.
            response.u_mouth_over_u_in = response.u_mouth_over_u_in * upstream.u_mouth_over_u_in;
            response.p_mouth_over_u_in = response.p_mouth_over_u_in * upstream.u_mouth_over_u_in;
            response.impedance = upstream.impedance;
            response.min_shear_number = response.min_shear_number.min(upstream.min_shear_number);
        }
        Ok(response)
    }

    /// Realize the regional driving-point reflectance using the SAME FIR owner.
    ///
    /// The physical return time uses local sound speeds. Port impedance uses
    /// the inlet gas and inlet area. `n_fft_override` has exactly the uniform
    /// line's meaning, budgets and finite-history approximation.
    ///
    /// # Errors
    /// As [`Self::response`], plus original history/Fourier/filter admission.
    pub fn characteristic_line(&self, termination: Termination, sample_rate_hz: u32,
        n: usize, wall: Option<&WallPin>, n_fft_override: Option<usize>)
        -> Result<DelayedFilter, DrivingPointError> {
        // Delegate the uniform case verbatim, including the order of travel-time
        // arithmetic at the FFT-size thresholds.
        if self.regions.len() == 1 {
            let r = &self.regions[0];
            return super::characteristic_line_dense(&r.duct, &r.gas, termination,
                sample_rate_hz, n, self.inlet_impedance(), wall, n_fft_override);
        }
        let dt = 1.0 / f64::from(sample_rate_hz);
        characteristic_line_response(2.0 * self.one_way_time_s / dt, termination,
            sample_rate_hz, n, self.inlet_impedance(), n_fft_override,
            |omega| self.response(omega, LossModel::Bessel, termination, wall))
    }

    /// Plane-wave impedance rho*c/S at the original inlet geometry.
    #[must_use]
    pub fn inlet_impedance(&self) -> f64 {
        let radius = self.regions[0].duct.segments[0].inlet_radius();
        let area = core::f64::consts::PI * radius * radius;
        self.inlet_gas().density * self.inlet_gas().sound_speed / area
    }
}

#[cfg(test)]
mod tests;
