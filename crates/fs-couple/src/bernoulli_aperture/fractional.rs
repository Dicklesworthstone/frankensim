//! Geometry/gas-derived sub-sample transits using the existing passive network.
//!
//! The original integer image stays available. This selected image preserves
//! low-frequency L/c rather than rounding length to the nearest acoustic cell.
//! It adds no fitted loss, changes no impedance, and retains actual allpass
//! storage. Finite-band phase error and numerical memory are explicit, not a
//! broadband-exact propagation or measured instrument-fidelity claim.
use super::network::{SectionRealization, TubeNetworkSpec};
use super::tube::UniformTubeSpec;
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::waveguide::fractional::FractionalTransit;
use fs_vfit::waveguide::network::{NetworkSegment, WaveguideNetwork};

/// Approximation admission, independent of gas/material source authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FractionalTubePolicy {
    /// Highest modeled physical frequency [Hz], strictly below Nyquist.
    /// This is not an output filter or a claim that nonlinear forcing is bandlimited.
    pub maximum_frequency_hz: f64,
    /// Sum of the one-way section phase-error bounds at that frequency [rad].
    /// This bounds a traversal visiting each section at most once, not resonant
    /// amplification, cycles, multiple returns or an entire transfer function.
    pub maximum_summed_phase_error_rad: f64,
}

/// Original geometry and actual selected phase image, in realized section order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FractionalSection {
    /// Original requested length [m].
    pub requested_length_m: f64,
    /// Rounded-coefficient zero-frequency effective length [m].
    pub realized_dc_length_m: f64,
    /// Original dimensionless request L/(c dt), before splitting.
    pub requested_samples: f64,
    /// Complete positive integer buffer length.
    pub integer_samples: usize,
    /// Explicit residual model, including coefficient-rounding information.
    pub transit: FractionalTransit,
    /// One-way phase discrepancy allowance in the policy's physical band [rad].
    pub phase_error_bound_rad: f64,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn wave(error: fs_vfit::waveguide::WaveguideError) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(error.to_string())
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {Err(AcousticRealizeError::Cancelled)} else {Ok(())}
}

impl TubeNetworkSpec {
    /// Lower each section's original length and coherent gas into integer plus
    /// normalized-allpass travel. Reuse the original impedance, topology, load
    /// and memory owners; there is no material-name branch or pitch input.
    ///
    /// SectionRealization's integer count describes only the buffer on this
    /// image. Its represented length includes fractional DC transit; the third
    /// return value retains all residuals and their band errors explicitly.
    /// "Original" here means this graph's supplied section geometry. If another
    /// reduction already rounded that geometry (such as viscothermal subdivision),
    /// this operation does not recover its earlier source lengths or certify that
    /// earlier approximation; retain and inspect the reduction's own records.
    /// The section's original length tolerance still applies to DC representation.
    /// A section shorter than one complete sample refuses; it needs a finer
    /// clock or an explicitly different lumped model, not an invented delay.
    ///
    /// Payload budget covers returned geometry and transit records, temporary
    /// segment/transit arrays and the original network, not caller data or RSS.
    pub fn realize_fractional_waveguide(&self, gases: &[GasState], dt: f64,
        policy: FractionalTubePolicy, gate: &CancelGate)
        -> Result<(WaveguideNetwork, Vec<SectionRealization>, Vec<FractionalSection>), AcousticRealizeError>
    {
        poll(gate)?;
        self.validate_section_gases(gases)?;
        let omega=core::f64::consts::TAU*policy.maximum_frequency_hz*dt;
        if !dt.is_finite() || dt<=0.0 || !policy.maximum_frequency_hz.is_finite()
            || policy.maximum_frequency_hz<=0.0 || !omega.is_finite()
            || omega<=0.0 || omega>=core::f64::consts::PI
            || !policy.maximum_summed_phase_error_rad.is_finite()
            || policy.maximum_summed_phase_error_rad<0.0 {
            return Err(invalid("fractional tube requires a resolved clock and explicit finite phase budget"));
        }
        let bytes=self.sections.len().checked_mul(core::mem::size_of::<NetworkSegment>()
            +core::mem::size_of::<FractionalTransit>()+core::mem::size_of::<FractionalSection>()
            +core::mem::size_of::<SectionRealization>())
            .ok_or_else(||invalid("fractional tube geometry size overflow"))?;
        let budget=self.max_wave_memory_bytes.checked_sub(bytes)
            .ok_or_else(||invalid("fractional tube geometry exceeds its memory allowance"))?;
        let mut segments=Vec::new();let mut transits=Vec::new();
        let mut represented=Vec::new();let mut details=Vec::new();
        segments.try_reserve_exact(self.sections.len()).map_err(|_|invalid("fractional geometry allocation failed"))?;
        transits.try_reserve_exact(self.sections.len()).map_err(|_|invalid("fractional geometry allocation failed"))?;
        represented.try_reserve_exact(self.sections.len()).map_err(|_|invalid("fractional geometry allocation failed"))?;
        details.try_reserve_exact(self.sections.len()).map_err(|_|invalid("fractional geometry allocation failed"))?;
        let mut phase_sum=0.0;
        for (section,gas) in self.sections.iter().zip(gases) {
            poll(gate)?;
            // The original uniform-tube owner validates physical geometry and Z.
            let z=UniformTubeSpec {length_m:section.length_m,radius_m:section.radius_m,
                sound_speed_m_s:gas.sound_speed,terminal_reflection:0.0,
                max_length_error_m:section.max_length_error_m,max_wave_memory_bytes:0}
                .characteristic_impedance(gas.density)?;
            let cell=gas.sound_speed*dt;
            let requested=section.length_m/cell;
            if !cell.is_finite() || cell<=0.0 || !requested.is_finite()
                || requested<1.0 || requested>=usize::MAX as f64
                || !section.max_length_error_m.is_finite() || section.max_length_error_m<0.0 {
                return Err(invalid("fractional section needs finite geometry and at least one full transit sample"));
            }
            let integer=requested.floor() as usize;
            let transit=FractionalTransit::new(requested-integer as f64).map_err(wave)?;
            let actual=(integer as f64+transit.realized_dc_samples())*cell;
            if !actual.is_finite() || (actual-section.length_m).abs()>section.max_length_error_m {
                return Err(invalid("fractional DC geometry exceeds the original length tolerance"));
            }
            let phase=transit.phase_error_bound_rad(omega).map_err(wave)?;
            phase_sum+=phase;
            if !phase_sum.is_finite() || phase_sum>policy.maximum_summed_phase_error_rad {
                return Err(invalid("fractional tube exceeds its declared accumulated phase-error budget"));
            }
            segments.push(NetworkSegment {nodes:section.nodes,one_way_samples:integer,impedance_pa_s_m3:z});
            transits.push(transit);
            represented.push(SectionRealization {represented_length_m:actual,one_way_samples:integer,impedance_pa_s_m3:z});
            details.push(FractionalSection {requested_length_m:section.length_m,realized_dc_length_m:actual,
                requested_samples:requested,integer_samples:integer,transit,phase_error_bound_rad:phase});
        }
        let network=WaveguideNetwork::with_fractional_transits(&self.nodes,&segments,&transits,dt,budget).map_err(wave)?;
        poll(gate)?;
        Ok((network,represented,details))
    }
}
