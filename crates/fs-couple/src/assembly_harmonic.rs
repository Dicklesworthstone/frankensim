//! Original admitted MIMO loads on the accepted loaded shell assembly.
//!
//! The existing rational source supplies Z(-i*omega), and the original physical
//! attachment owner supplies its work-conjugate map. fs-modal owns the loaded
//! dynamic solve. There is no new fit, state realization or frequency interpolation.
//! The source's declared band and two admission-report classes remain visible;
//! a pointwise harmonic solve does not promote them into a broadband certificate.

use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicPortControls, HarmonicResponseControls};
use fs_plate::shell::preload::assembly::vibration::harmonic::ports::{
    AssemblyHarmonicPorts, AssemblyPortHarmonicResponse, PreparedAssemblyPortHarmonicResponse,
};
use fs_vfit::mimo::MimoPassivityReport;
use fs_vfit::mimo::load::{MimoLoadError, PreparedMimoLoad};

/// Source-band/admission, structural solve or cancellation refusal.
#[derive(Debug)]
pub enum AssemblyMimoHarmonicError {
    /// Cancellation before the prepared response was published.
    Cancelled,
    /// The original continuous load refused this frequency.
    Source(MimoLoadError),
    /// Attachment, resource, dynamic solve or physical work refusal.
    Structure(String),
}

impl core::fmt::Display for AssemblyMimoHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "assembled MIMO harmonic response refused: {self:?}")
    }
}
impl std::error::Error for AssemblyMimoHarmonicError {}

/// One source-bound continuous frequency, reusable for physical forcing cases.
pub struct PreparedAssemblyMimoHarmonicResponse<'l, 'p, 'm, 'v, 'body> {
    source: &'l PreparedMimoLoad,
    angular_frequency_rad_s: f64,
    response: PreparedAssemblyPortHarmonicResponse<'p, 'm, 'v, 'body>,
}

impl<'l, 'p, 'm, 'v, 'body> PreparedAssemblyMimoHarmonicResponse<'l, 'p, 'm, 'v, 'body> {
    /// Use the already admitted proper continuous operator in its original port
    /// order. The supplied attachment units/order must match that source. The
    /// load's clock, pole history and passivity reports are borrowed unchanged.
    /// In particular, a nonzero runtime history never becomes periodic forcing.
    pub fn new(
        ports: &'p AssemblyHarmonicPorts<'m, 'v, 'body>,
        source: &'l PreparedMimoLoad,
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        port_controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<Self, AssemblyMimoHarmonicError> {
        let impedance = source_impedance(
            source,
            ports.port_count(),
            angular_frequency_rad_s,
            port_controls,
            cx,
        )?;
        let response = ports
            .prepare(
                angular_frequency_rad_s,
                &impedance,
                controls,
                port_controls,
                cx,
            )
            .map_err(AssemblyMimoHarmonicError::Structure)?;
        poll(cx)?;
        Ok(Self {
            source,
            angular_frequency_rad_s,
            response,
        })
    }

    /// Admitted continuous driving frequency [rad/s].
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.angular_frequency_rad_s
    }

    /// Original admitted source, with its band, port order and unchanged history.
    pub fn source(&self) -> &PreparedMimoLoad {
        self.source
    }

    /// Continuous and actual sampled reports, in that order. GridOnly remains
    /// GridOnly; neither report supplies a load storage/heat decomposition.
    pub fn source_reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        self.source.reports()
    }

    /// The returned external-load account is physical boundary work. It must not
    /// be relabeled material heat or far-field radiation without that source's
    /// independent physical owner. Structural glue loss is counted separately.
    pub fn solve(
        &self,
        component_loads: &[&[C64]],
        cx: &Cx<'_>,
    ) -> Result<AssemblyPortHarmonicResponse, AssemblyMimoHarmonicError> {
        self.response
            .solve(component_loads, cx)
            .map_err(AssemblyMimoHarmonicError::Structure)
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), AssemblyMimoHarmonicError> {
    cx.checkpoint()
        .map_err(|_| AssemblyMimoHarmonicError::Cancelled)
}

fn source_impedance(
    source: &PreparedMimoLoad,
    ports: usize,
    omega: f64,
    controls: HarmonicPortControls,
    cx: &Cx<'_>,
) -> Result<Vec<C64>, AssemblyMimoHarmonicError> {
    poll(cx)?;
    if ports != source.ports() {
        return Err(AssemblyMimoHarmonicError::Structure(
            "physical attachment and original MIMO port counts disagree".into(),
        ));
    }
    // Preflight before the original evaluator allocates its one full matrix.
    // The shared loaded solve applies its stronger complete payload/work bounds.
    let bytes = (ports as u128) * (ports as u128) * core::mem::size_of::<C64>() as u128;
    if ports > controls.maximum_ports || bytes > controls.maximum_payload_bytes as u128 {
        return Err(AssemblyMimoHarmonicError::Structure(
            "original MIMO impedance exceeds the declared port or payload limit".into(),
        ));
    }
    source
        .continuous_operator()
        .admit_frequency(omega)
        .map_err(AssemblyMimoHarmonicError::Source)?;
    // fs-vfit's original Laplace model uses s=+i*w for exp(+i*w*t).
    // Evaluate the negative axis directly for the structural exp(-i*w*t) chart.
    let impedance = source.source_model().eval(C64::new(0.0, -omega));
    poll(cx)?;
    Ok(impedance)
}

#[cfg(test)]
mod tests;
