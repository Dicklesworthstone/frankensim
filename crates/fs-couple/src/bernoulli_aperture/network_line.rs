//! Material valves on the original regional duct-network driving point.
//!
//! For the admitted FIR, `p_minus = r0*p_plus + h` and
//! `p = Zc*Q + 2*p_minus`. Eliminating the SAME-sample coefficient gives
//! `Zeff = Zc*(1+r0)/(1-r0)` and `a_eff = h/(1-r0)`. The original single,
//! paired or material-modal aperture solves this load before any owner advances.
//! Its accepted flow is then inserted into the original `DelayedFilter`; there
//! is no added feedback delay or alternative convolution/valve solver.
//!
//! This is a static small-signal network medium with a nonlinear moving source.
//! It retains the existing line ladder, Nyquist taper and precursor folding;
//! exact algebra for that finite sampled load is not a claim of exact continuous
//! propagation. Material, contact and slit-domain refusals remain the valve's.
//! Source-plane pressure is internal diagnostic pressure, not a microphone.
//! The local mechanical work account is retained; no FIR storage functional or
//! complete network energy certificate is invented here.

/// Canonical finite pressure gestures on this retained physical composition.
pub mod performance;

use super::dynamic::{ApertureDrive, ApertureFrame};
use super::modal::pressure::MaterialPressureFrame;
use super::tube::TubeDrive;
use super::two_sided::LoadedAperture;
use crate::acoustic_realize::AcousticRealizeError;
use crate::driving_point::network::{CoupledSource, NetworkDrivingPoint};
use crate::driving_point::{DrivingPointError, LineLadder, LineRealization};
use crate::wind_line::{WindLineError, history_line};
use fs_duct::network::NetworkMedia;
use fs_duct::{DuctError, LossModel};
use fs_material::gas_species::input::GasRecord;
use fs_vfit::discretize::{DelayedFilter, DiscretizeError};

/// Original owner refusals, including the complete network ladder diagnostic.
#[derive(Debug)]
pub enum NetworkLineError {
    /// A composition does not identify a finite, matching physical load.
    Invalid {
        /// Failed requirement.
        what: &'static str,
    },
    /// Geometry, regional material or coupled boundary refusal.
    Duct(DuctError),
    /// The original sampled-response ladder refused.
    Realize(DrivingPointError),
    /// The shared characteristic-history owner refused.
    History(WindLineError),
    /// An atomic filter push refused.
    Discrete(DiscretizeError),
    /// The original mechanical/material/contact owner refused.
    Aperture(AcousticRealizeError),
}

impl core::fmt::Display for NetworkLineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "network material valve: {what}"),
            Self::Duct(e) => write!(f, "network material valve geometry/medium: {e}"),
            Self::Realize(e) => write!(f, "network material valve realization: {e}"),
            Self::History(e) => write!(f, "network material valve history: {e}"),
            Self::Discrete(e) => write!(f, "network material valve filter: {e}"),
            Self::Aperture(e) => write!(f, "network material valve source: {e}"),
        }
    }
}
impl core::error::Error for NetworkLineError {}

fn invalid(what: &'static str) -> NetworkLineError {
    NetworkLineError::Invalid { what }
}

/// One realized network retaining the complete geometry, regional source records,
/// coupled-boundary identities, and original admission/convergence report.
/// Construction resolves the physical source plane before constructing a valve.
pub struct NetworkLineLoad {
    driving_point: NetworkDrivingPoint,
    media: NetworkMedia,
    realization: LineRealization,
    coupled_source_ids: Vec<String>,
    sample_rate_hz: u32,
    density_kg_m3: f64,
    characteristic_impedance_pa_s_m3: f64,
    effective_impedance_pa_s_m3: f64,
    direct_reflection: f64,
    tap_count: usize,
    tail: DelayedFilter,
}

impl NetworkLineLoad {
    /// Realize every edge with its own source gas and the explicit exterior.
    /// The active plane must be one edge END; its actual radius, density and
    /// sound speed set `Zc = rho*c/area`. A coupled source is required exactly
    /// for declared coupled ports and all distinct supplied source IDs survive.
    /// The complete FIR is read directly, including taps after sample 4096.
    pub fn new(
        driving_point: NetworkDrivingPoint,
        media: NetworkMedia,
        loss: LossModel,
        sample_rate_hz: u32,
        ladder: &LineLadder,
        coupled_at: Option<CoupledSource<'_>>,
    ) -> Result<Self, NetworkLineError> {
        let gas = driving_point
            .active_port_record(&media)
            .map_err(NetworkLineError::Duct)?
            .state();
        let radius = driving_point
            .active_port_radius()
            .map_err(NetworkLineError::Duct)?;
        let density = gas.density;
        let zc = density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        if sample_rate_hz == 0 || !zc.is_finite() || zc <= 0.0 {
            return Err(invalid(
                "source plane needs finite positive impedance and sample rate",
            ));
        }
        let mut ids = Vec::new();
        let (line, realization) = match coupled_at {
            Some(source) => {
                let mut retained = |omega| {
                    let boundary = source(omega)?;
                    if !ids.contains(&boundary.source_id) {
                        ids.push(boundary.source_id.clone());
                    }
                    Ok(boundary)
                };
                driving_point.realize_line_regional(
                    &media,
                    loss,
                    sample_rate_hz,
                    zc,
                    ladder,
                    Some(&mut retained),
                )
            }
            None => {
                driving_point.realize_line_regional(&media, loss, sample_rate_hz, zc, ladder, None)
            }
        }
        .map_err(NetworkLineError::Realize)?;
        let (dt, taps) = line
            .quiescent_fir()
            .ok_or_else(|| invalid("network realization must retain a quiescent plain FIR"))?;
        let direct = taps[0];
        let effective = zc * (1.0 + direct) / (1.0 - direct);
        if !(direct > -1.0 && direct < 1.0 && effective.is_finite() && effective > 0.0)
            || dt.to_bits() != (1.0 / f64::from(sample_rate_hz)).to_bits()
        {
            return Err(invalid(
                "direct reflection must yield a finite positive same-sample valve load",
            ));
        }
        let tail = history_line(dt, taps, &[]).map_err(NetworkLineError::History)?;
        Ok(Self {
            driving_point,
            media,
            realization,
            coupled_source_ids: ids,
            sample_rate_hz,
            density_kg_m3: density,
            characteristic_impedance_pa_s_m3: zc,
            effective_impedance_pa_s_m3: effective,
            direct_reflection: direct,
            tap_count: taps.len(),
            tail,
        })
    }

    /// Original complete graph, loads and port schema.
    #[must_use]
    pub fn driving_point(&self) -> &NetworkDrivingPoint {
        &self.driving_point
    }
    /// Full edge and exterior material records; no global air default.
    #[must_use]
    pub fn media(&self) -> &NetworkMedia {
        &self.media
    }
    /// Actual source gas record, selected by edge identity and endpoint.
    #[must_use]
    pub fn source_record(&self) -> &GasRecord {
        self.driving_point
            .active_port_record(&self.media)
            .expect("immutable admitted source plane")
    }
    /// Every rung and original convergence diagnostic.
    #[must_use]
    pub fn realization(&self) -> &LineRealization {
        &self.realization
    }
    /// Distinct coupled boundary source identities in first-evaluation order.
    #[must_use]
    pub fn coupled_source_ids(&self) -> &[String] {
        &self.coupled_source_ids
    }
    /// Mechanical and line sample clock [Hz].
    #[must_use]
    pub fn sample_rate_hz(&self) -> u32 {
        self.sample_rate_hz
    }
    /// Actual source gas density required by the retained valve [kg/m³].
    #[must_use]
    pub fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    /// Source-plane characteristic pressure/volume-flow impedance [Pa s/m³].
    #[must_use]
    pub fn characteristic_impedance_pa_s_m3(&self) -> f64 {
        self.characteristic_impedance_pa_s_m3
    }
    /// Same-sample eliminated impedance to declare when building the valve.
    #[must_use]
    pub fn effective_impedance_pa_s_m3(&self) -> f64 {
        self.effective_impedance_pa_s_m3
    }
    /// Original sampled reflection's coefficient at zero delay.
    #[must_use]
    pub fn direct_reflection(&self) -> f64 {
        self.direct_reflection
    }
    /// Complete original FIR length; no fixed probing horizon is used.
    #[must_use]
    pub fn tap_count(&self) -> usize {
        self.tap_count
    }
    /// Strictly past-dependent reflected pressure for the next solve [Pa].
    #[must_use]
    pub fn history_pressure_pa(&self) -> f64 {
        self.tail.incoming()
    }
    /// Physical outgoing-wave history, oldest first, for inspection only.
    #[must_use]
    pub fn outgoing_history(&self) -> Vec<f64> {
        self.tail.history()
    }
}

/// Accepted physical source-plane state, preserving the original valve account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkLineFrame {
    /// Mechanical work/contact/material results; outgoing is mapped to physical
    /// `p_plus`, while bore pressure and volume flow are the unchanged solution.
    pub aperture: ApertureFrame,
    /// Same-sample physical `p_minus` [Pa], including the instantaneous load.
    pub incoming_pressure_pa: f64,
    /// Past-dependent `h` used by this solve [Pa].
    pub history_pressure_pa: f64,
    /// Original modal spatial pressure results, in original source port order.
    pub spatial_pressure: [MaterialPressureFrame; 2],
}

/// Original valve and network histories share one atomic acceptance boundary.
pub struct LoadedApertureNetworkLine {
    aperture: LoadedAperture,
    load: NetworkLineLoad,
    last: Option<NetworkLineFrame>,
}

impl LoadedApertureNetworkLine {
    /// Bind an unadvanced physical valve, using the compiled load's actual gas
    /// density, effective impedance and clock. Any retained material/slit gas
    /// must also match in full, including transport properties. Never rewrite
    /// source coefficients.
    pub fn new(aperture: LoadedAperture, load: NetworkLineLoad) -> Result<Self, NetworkLineError> {
        let (rho, impedance) = aperture.load();
        let gas = load.source_record().state();
        let mismatched_gas = match &aperture {
            LoadedAperture::Single(a) => a
                .viscous_slit()
                .is_some_and(|law| law.gas().is_some_and(|retained| retained != gas)),
            LoadedAperture::Paired(a) => a.gas().is_some_and(|retained| retained != gas),
            LoadedAperture::Modal(a) => a.gas().is_some_and(|retained| retained != gas),
        };
        if aperture.accepted_steps() != 0
            || mismatched_gas
            || rho.to_bits() != load.density_kg_m3.to_bits()
            || impedance.to_bits() != load.effective_impedance_pa_s_m3.to_bits()
            || aperture.time_step_s().to_bits() != (1.0 / f64::from(load.sample_rate_hz)).to_bits()
        {
            return Err(invalid(
                "cold valve must declare the retained source gas, eliminated impedance and clock",
            ));
        }
        Ok(Self {
            aperture,
            load,
            last: None,
        })
    }
    /// Original single, paired or material-modal owner, including all memory.
    #[must_use]
    pub fn aperture(&self) -> &LoadedAperture {
        &self.aperture
    }
    /// Retained physical graph, material records and original ladder.
    #[must_use]
    pub fn load(&self) -> &NetworkLineLoad {
        &self.load
    }
    /// Most recent accepted source-plane observation.
    #[must_use]
    pub fn last_frame(&self) -> Option<&NetworkLineFrame> {
        self.last.as_ref()
    }
    /// Advance without independent structural loading.
    pub fn step(&mut self, drive: TubeDrive) -> Result<NetworkLineFrame, NetworkLineError> {
        self.step_with_forces_and_pressures(drive, [0.0; 2], [0.0; 2])
    }
    /// Advance with the original signed force ports and optional material-mesh
    /// pressure ports. An absent port refuses instead of discarding a drive.
    /// A failed preview, mapping or filter push changes neither physical owner.
    pub fn step_with_forces_and_pressures(
        &mut self,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
    ) -> Result<NetworkLineFrame, NetworkLineError> {
        let h = self.load.tail.incoming();
        let denominator = 1.0 - self.load.direct_reflection;
        let trial = self
            .aperture
            .preview(
                ApertureDrive {
                    upstream_pressure_pa: drive.upstream_pressure_pa,
                    incoming_pressure_pa: h / denominator,
                    body_flow_m3_s: drive.body_flow_m3_s,
                },
                forces_n,
                pressures_pa,
            )
            .map_err(NetworkLineError::Aperture)?;
        let mut aperture = trial.frame();
        let outgoing = (self.load.characteristic_impedance_pa_s_m3 * aperture.bore_flow_m3_s + h)
            / denominator;
        let incoming = self.load.direct_reflection * outgoing + h;
        if !outgoing.is_finite() || !incoming.is_finite() {
            return Err(invalid("physical source-plane wave mapping is nonfinite"));
        }
        aperture.outgoing_pressure_pa = outgoing;
        let frame = NetworkLineFrame {
            aperture,
            incoming_pressure_pa: incoming,
            history_pressure_pa: h,
            spatial_pressure: trial.spatial_pressure(),
        };
        // DelayedFilter validates its complete next tail before mutating memory.
        // After that single atomic push, accepting the exclusive valve trial is
        // infallible. No second convolution or state clone is needed per sample.
        self.load
            .tail
            .push(outgoing)
            .map_err(NetworkLineError::Discrete)?;
        self.aperture.accept(trial);
        self.last = Some(frame);
        Ok(frame)
    }
}
