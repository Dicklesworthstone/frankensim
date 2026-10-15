//! Quasi-steady separated flow through a fixed opening, using the existing
//! continuous Bernoulli owner. Effective jet area is Cc times opening area:
//! dp = rho Q |Q| / (2 (Cc A)^2). This is a symmetric, memoryless loss model,
//! not resolved vortex dynamics, radiation, or compressible mean flow.

use super::moving;
use crate::acoustic_realize::AcousticRealizeError;
use fs_material::gas::GasState;

/// Explicit empirical closure and its use limit; no instrument preset.
#[derive(Clone, Debug, PartialEq)]
pub struct ContractionSpec {
    /// Contracted jet area / physical opening area, in (0,1]. The same value
    /// applies in both flow directions. Unity still has Bernoulli separation
    /// loss; omit this model entirely for a lossless comparison.
    pub area_ratio: f64,
    /// Maximum |Q|/(Cc A c), in (0,0.3]. Checked on the solved flow, not on
    /// blowing pressure or a previous sample. Not a full validity certificate.
    pub maximum_contracted_mach: f64,
    /// Source/condition identity or an explicit authored-model label.
    pub provenance: String,
}

/// One consistent flow/drop observation; positive power is energy lost to the
/// unresolved separated flow, NOT far-field acoustic radiation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ContractionFrame {
    /// Signed through-flow [m^3/s].
    pub volume_flow_m3_s: f64,
    /// Signed opposing pressure drop [Pa].
    pub pressure_drop_pa: f64,
    /// dp * Q [W], nonnegative in either direction.
    pub dissipated_power_w: f64,
}

/// Cold binding to physical geometry and one complete coherent material gas.
/// No mutable coefficients or independent density overrides escape this object.
#[derive(Clone, Debug, PartialEq)]
pub struct ContractionLoss {
    spec: ContractionSpec,
    gas: GasState,
    width_m: f64,
    opening_m: f64,
    effective_area_m2: f64,
    unit_pressure_flow: f64,
}
fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
impl ContractionLoss {
    /// Bind an explicitly supplied area contraction to an actual fixed window.
    /// Density and sound speed come from the SAME gas state. No source identity
    /// implies experimental accuracy; the quasi-steady closure neglects vortex
    /// memory, pressure recovery and direction-dependent edge geometry.
    ///
    /// # Errors
    /// Invalid geometry, coefficient, source, gas coherence or representability.
    pub fn new(width_m: f64, opening_m: f64, gas: GasState, spec: ContractionSpec)
        -> Result<Self, AcousticRealizeError>
    {
        gas.validate_coherence().map_err(|_| invalid("contraction needs one coherent material gas"))?;
        if ![width_m, opening_m, spec.area_ratio, spec.maximum_contracted_mach]
            .iter().all(|x| x.is_finite() && *x > 0.0)
            || spec.area_ratio > 1.0 || spec.maximum_contracted_mach > 0.3
            || spec.provenance.trim().is_empty()
        { return Err(invalid("contraction needs positive geometry, explicit source, Cc<=1 and Mach<=0.3")); }
        let effective_width = spec.area_ratio * width_m;
        let effective_area_m2 = effective_width * opening_m;
        // The original Bernoulli law supplies all dimensional scaling. Its
        // continuous extension below the old pressure cutoff remains in moving.
        let unit_pressure_flow = moving::volume_flow(effective_width, opening_m, 1.0, gas.density);
        if ![effective_width, effective_area_m2, unit_pressure_flow]
            .iter().all(|x| x.is_finite() && *x > 0.0)
        { return Err(invalid("contraction geometry or Bernoulli scale is not representable")); }
        Ok(Self { spec, gas, width_m, opening_m, effective_area_m2, unit_pressure_flow })
    }
    /// Original coefficient, source and explicit contracted-flow limit.
    #[must_use] pub const fn spec(&self) -> &ContractionSpec { &self.spec }
    /// Complete material gas used by the pressure drop and speed admission.
    #[must_use] pub const fn gas(&self) -> &GasState { &self.gas }
    /// Actual opening width [m].
    #[must_use] pub const fn width_m(&self) -> f64 { self.width_m }
    /// Actual opening height [m], not the thickness of an exciting jet.
    #[must_use] pub const fn opening_m(&self) -> f64 { self.opening_m }
    /// Cc A [m^2], retained from the original geometry.
    #[must_use] pub const fn effective_area_m2(&self) -> f64 { self.effective_area_m2 }

    /// Invert the original Bernoulli flow law without an invented low-flow cut.
    ///
    /// # Errors
    /// Nonfinite flow, solved contracted Mach outside the supplied domain, or
    /// unrepresentable pressure/power (including nonzero loss underflow).
    pub fn at_flow(&self, flow: f64) -> Result<ContractionFrame, AcousticRealizeError> {
        let mach = (flow.abs() / self.effective_area_m2) / self.gas.sound_speed;
        if !flow.is_finite() || !mach.is_finite() || mach > self.spec.maximum_contracted_mach {
            return Err(invalid("solved contraction flow exceeds its finite low-Mach domain"));
        }
        let root = flow / self.unit_pressure_flow;
        let pressure = root * root.abs();
        let power = pressure * flow;
        if ![pressure, power].iter().all(|x| x.is_finite())
            || (flow != 0.0 && (pressure == 0.0 || power == 0.0)) {
            return Err(invalid("contraction pressure or dissipated power is not representable"));
        }
        Ok(ContractionFrame { volume_flow_m3_s: flow, pressure_drop_pa: pressure,
            dissipated_power_w: power })
    }

    /// Simultaneously solve drive = R Q + dp(Q), including R=0. R is the
    /// caller's discrete passive response [Pa s/m^3], not another flow law.
    /// Uses the SAME rationalized pressure partition as moving reed valves;
    /// no explicit old-flow loss, iterations, clipping or audio normalization.
    ///
    /// # Errors
    /// Invalid response/drive, original Bernoulli refusal, unresolved pressure
    /// equation or an at_flow refusal at the candidate (never the old flow).
    pub fn solve(&self, resistance: f64, drive_pa: f64)
        -> Result<ContractionFrame, AcousticRealizeError>
    {
        if !resistance.is_finite() || resistance < 0.0 || !drive_pa.is_finite() {
            return Err(invalid("contraction load must be finite and nonnegative"));
        }
        let flow = if resistance == 0.0 {
            moving::volume_flow(self.spec.area_ratio*self.width_m, self.opening_m,
                drive_pa, self.gas.density)
        } else {
            let (_, linear_pressure) = moving::pressure_partition(drive_pa,
                resistance*self.unit_pressure_flow)?;
            // Retain the linear component directly: squaring the much smaller
            // nonlinear drop first can underflow before Q itself does.
            linear_pressure / resistance
        };
        let frame = self.at_flow(flow)?;
        let linear = resistance*flow;
        let residual = linear + frame.pressure_drop_pa - drive_pa;
        let scale = linear.abs().max(frame.pressure_drop_pa.abs()).max(drive_pa.abs());
        if ![linear, residual].iter().all(|x| x.is_finite())
            || residual.abs() > 512.0*f64::EPSILON*scale {
            return Err(invalid("contraction pressure equation is unresolved"));
        }
        Ok(frame)
    }
}
