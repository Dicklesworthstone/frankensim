//! A lip-valve performance file: a measured lip, a held bore, a pressure
//! phrase and a physical outlet microphone.
//!
//! Nothing here names an instrument, a note or a pitch. The file supplies
//! SI primitives and their provenance; the note that sounds is the lock the
//! loop finds between the lip and the bore ([`super::BrassVoice`]). The
//! radiated sound is the retarded derivative of the bore's MOUTH VOLUME FLOW
//! at a baffled point, never the mouthpiece pressure:
//!
//! `p+` (wave sent into the bore) -> [`plane_flow_line`] (the same multimodal
//! bore's total `U_mouth / p+`) -> [`BaffledPressure`] (Rayleigh integral of
//! the actual discrete outlet disk).
//!
//! ```text
//! frankensim-lip-performance-v1
//! # comment lines and blank lines are ignored
//! audio <rate_hz> <samples>                  output clock (also the mechanics clock of the lines)
//! full_scale_pa <pa>
//! ambient moist_air <temperature_k> <pressure_pa> <relative_humidity>
//! ambient mixture <temperature_k> <pressure_pa> <mole|mass> <count> <species fraction>...
//! ambient thermochemical <temperature_k> <pressure_pa> <mole|mass> <count> <id fraction>...
//! bore_gas <the same legacy, mixture or pack record as ambient>  optional; the gas the bore, cup and lip hold (needs modal_baffle, or on bore nonlinear a thermochemical record and mouth_load)
//! lip_response <f_hz> <q> <inverse_mu_m2_per_kg> <width_m> <rest_gap_m> <face_area_m2>
//! lip_response_mode2 <f_hz> <q> <inverse_mu_m2_per_kg>      optional; the lip's second, inward-striking mode (a negative 1/mu)
//! lip_area_exponent <q>                                    optional; flow area w h (h/H0)^(q-1), 1 <= q <= 4 (q = 1 is the plain w h)
//! lip_rim_stop inelastic|elastic                           optional; the lip gap may not exceed the first bore cone's diameter (two limit cases of an unmeasured restitution)
//! bore_wall rigid|thin-brass <thickness_m> <quality>        optional; on bore nonlinear the bell's breathing wall: every cell a radial oscillator of cartridge brass (rigid is the row absent)
//! bore_wall shell-brass <thickness_m> <loss_factor> <max_frequency_hz> <azimuths>  optional alternative to thin-brass: the bell as a MITC3+ shell of brass (its bending modes up to the cut-off, on a mesh of that many points round, even, 16 to 128), throat ring clamped
//! bore_wall_source <free text to end of line>               required with bore_wall thin-brass: whose instrument and whose measurement the thickness and the quality factor are
//! player_tract mouth <tau_s>                               optional (or `player_tract lung`); a player's tract upstream of the lip, the card's pressure measured in the mouth or at the lung
//! player_tract_section <length_m> <area_m2>                1 to 256 rows beside `player_tract`, lip side first, each a cylinder of that area (the glottis a section of its own)
//! player_tract_wall <surface_density> <stiffness> <resistance>  optional; every tract tube's locally reacting wall (rigid without it)
//! player_airway weibel-1963 saab-2021-table-1-B <generations>   optional; the first generations of Weibel's tree beyond the last section
//! lip <mass_kg> <stiffness_n_m> <damping_n_s_m> <width_m> <rest_gap_m> <face_area_m2>
//! lip_source <free text to end of line>      required provenance sentence
//! cup_volume_m3 <volume>
//! mouth_reservoir <volume_m3> <feed_area_m2> <initial_pressure_pa> <initial_temperature_k>  optional; finite upstream gas on bore nonlinear
//! termination <flanged|unflanged|ideal_open>
//! modes <n_modes> <extra_slices>
//! mouth_flow_tail_tol <relative_energy>      optional; default 1e-10 (the line ladder's own)
//! mouth_flow_max_taps <taps>                 optional; power of two in 512..=262144; default 65536 (the ladder's own ceiling)
//! modal_baffle <relative_kernel_tol> <maximum_kernel_terms>  optional; full modal load and observation
//! mouth_load modal <n_modes> <extra_slices> <band_hz> <bin_hz> <tail_tolerance>  optional; on bore nonlinear: the bore's load is the plane-mode equivalent of the complete flanged multimodal mouth, not the baffled piston
//! jet_noise <seed> <gain> <st_lo> <st_hi> <re_on>            optional; constriction noise of the lip jet, Verge's fraction
//! constriction_noise <seed> <duct_radius_m>                  optional; Oldham and Ukpoho's printed constriction noise (both bores)
//! cylinder <radius_m> <length_m>             bore records, mouthpiece to bell
//! cone <inlet_radius_m> <outlet_radius_m> <length_m>
//! bessel <inlet_radius_m> <outlet_radius_m> <length_m> <flare> <pieces>
//! bore_source <free text to end of line>     required provenance sentence
//! bore nonlinear <max_cell_m> <mechanics_per_output> <gas_steps_per_mechanics> <wall_poles> <wall_band_lo_hz> <wall_band_hi_hz> <wall_tolerance>
//! microphone <x_m> <y_m> <z_m> <radial_rings> <angular_points> <max_frequency_hz>
//! observer_gas <stokes-kirchhoff|lossless> <max_paths> <max_filter_taps> <max_kernel_terms>  optional
//! observer_gas molecular <max_paths> <max_filter_taps> <max_kernel_terms> <shared relaxation law>  optional
//! schedule
//! <canonical GestureSchedule bytes: blowing-pressure, optional pre-stress>
//! ```
//!
//! `observer_gas` applies the shared material propagation law to every original
//! bell-surface-to-microphone path, using the complete `ambient` gas and the
//! caller's work limits. Molecular laws bind that original gas record, including
//! its composition and source packs; species properties are never inferred from
//! a name. Both the linear and nonlinear total-flow receivers carry this law,
//! retaining every mechanics substep and the common causal filter latency.
//! It changes observation only, not lip motion, bore feedback or bell loading.
//! The signed `modal_baffle` receiver uses the same per-path law and checks its
//! actual material kernels on the complete spatial quadrature ladder. Without
//! the row, observation is unchanged.
//!
//! `mouth_reservoir` attaches a rigid adiabatic upstream gas volume to
//! `bore nonlinear`. Its pressure and temperature evolve with the bore's
//! signed mass and enthalpy exchanges. The blowing-pressure track controls
//! the held feed reservoir through the authored effective opening; it never
//! resets storage. Initial pressure is absolute, whereas the track is relative
//! to room pressure. Zero feed area seals the supply. Feed and stored gas use
//! the bore's original thermochemical source. The existing jet/constriction
//! noise rows stand beside finite storage in either row order: their summed
//! pressure is a series offset between the stored gas and the lip face, offered
//! to the lip at the stored gas's own temperature (an isothermal modulation,
//! as the original reservoir applies it at the room temperature); the feed and
//! the stored inventory never see it, and the work such a source would do is
//! reported by the finite-mouth step and rejected as heat, not drawn from the
//! storage. The constriction noise keeps its amplitude on the held blowing
//! pressure, the jet noise reads the stored gas's gauge pressure. This compact
//! gas volume supplies no anatomical lung, neck-inertia or moving-wall model.
//! Without this row the held mouth is unchanged.
//!
//! `modal_baffle` explicitly selects the complete flanged modal radiation load
//! and observes every signed mouth mode at the declared finite-distance point.
//! It requires `termination flanged` and a linear bore. The same original gas
//! and actual staircase aperture enter both load and receiver. Its microphone
//! radial count is Gauss-Legendre nodes; radial and angular counts must be
//! multiples of four for the quarter/half/full kernel convergence checks.
//! Files without that record retain the original total-flow disk observation.
//!
//! `bore_gas` separates the gas inside the instrument from the room. Without it one `ambient` gas serves both, as it
//! always did. With it the bore's propagation and wall loss, the line references, the cup compliance, and the lip's and
//! the constriction noise's density and viscosity use this gas (a played bore holds warm, humid, CO2-rich breath), while the
//! flanged modal radiation load, the receiver and the observer stay in the `ambient` room gas. It takes the grammar of
//! `ambient` (so `pack` names a sourced seed gas), needs `modal_baffle` (whose exterior load names the room gas explicitly;
//! the analytic mouth loads would radiate into the bore gas), and refuses a different static pressure (mean flow is not modeled).
//! On `bore nonlinear` it stands beside `mouth_load modal` instead (MB 24783): the bore evolves the `bore_gas` row's own NASA-9 laws, so the record must be `thermochemical`
//! with its packs; the kernel of the mouth load names the `ambient` room as its exterior (the two-gas flanged law, with the BORE gas's characteristic impedance as its
//! reference, which is the bore end's own), and the receiver and the observer stay in the room; the `ambient` row may then be any gas. A baffled-piston end would radiate
//! into the bore gas, so a nonlinear `bore_gas` without `mouth_load` is refused. Without the row nothing changes, bit for bit.
//! It is one uniform gas: no axial temperature or composition profile, no condensation. A file without the row is
//! bit-identical to one written before it existed.
//!
//! `jet_noise` is the opt-in constriction noise of the lip jet, the voice's row on the lip channel
//! ([`crate::turbulence_noise::LipJetNoise`]): a series
//! pressure source added to the lip's driving pressure, of rms `gain * rho U_j^2 / 2` over `St_h = f h / U_j` in
//! `[st_lo, st_hi]`, with `U_j = Q / (w h)` the Bernoulli speed of the lip channel and `h` the lip gap, shut while
//! `Re_h = Q / (w nu)` is below `re_on` or the jet flow is not forward. Every field is required. The form is Stevens' (1971),
//! the constant, band and gate are Verge's (1995, a flue jet at a labium, slit 1 mm), and putting them on a lip is an
//! extrapolation that any page showing a render must say: "constriction noise from a flue-jet measurement (Verge 1995: flue
//! jet at a labium, slit 1 mm), applied to the lip jet as an extrapolation; Stevens 1971 gives the form (a pressure source
//! proportional to the pressure drop), not this constant". Both bores carry it: the linear loop adds it to the lip's driving
//! pressure each sample, the nonlinear bore each gas step ([`super::nonlinear::NonlinearLip::with_jet_noise`]; there
//! `Re_h = mdot / (w mu)` with the wall binding's own viscosity). A file without the record is bit-identical to one written
//! before it existed, and so is one whose gate never opens.
//!
//! `constriction_noise` is the second opt-in source on the lip channel, the printed duct-constriction law of Oldham and
//! Ukpoho (J. Sound Vib. 140(2), 1990; Guedel, AIVC Contributed Report 02, eq. 3 and Fig. 10) in place of Verge's flue-jet
//! fraction ([`crate::turbulence_noise::ConstrictionNoise`]): a series pressure source added to the lip's driving pressure, in
//! each third octave below the duct cutoff `100 / r` whose Strouhal number `f d / U_c` lies on the printed line, of rms
//! `K(St) P_blow` with `K` the printed trend, `U_c = sqrt(2 P_blow / rho)` and `d = pi r / 2`. Two fields and no gain, so it
//! has nothing to tune: `r` is the duct the constriction opens into, the cup rim radius the bore record's first station
//! gives (the card says which). The law was measured on dampers and orifice plates in a 300 mm duct at 10-25 m/s, so putting
//! it on a lip is an extrapolation in size and Mach number that any page showing a render says. It is silent while the jet
//! flow is not forward, and a file without the record is bit-identical to one written before it existed. Both bores carry
//! it, beside (not instead of) `jet_noise` when both are given: the linear loop adds it to the lip's driving pressure each
//! sample, the nonlinear bore each gas step ([`super::nonlinear::NonlinearLip::with_constriction_noise`]; its bank runs at the
//! gas step, with the wall binding's own gas density for `U_c`).
//!
//! Exactly one of `lip_response` and `lip` is given. `lip_response` carries
//! the quantities lip-reed measurements publish: the outward-striking
//! resonance `f`, its quality factor `Q`, and `1/mu`, the pressure-collecting
//! area per effective mass (`mu = m / A`). The island follows
//! `m = A / (1/mu)`, `k = m (2 pi f)^2`, `c = m (2 pi f) / Q`; the face area
//! `A` cancels in every motion and only fixes the force normalization.
//!
//! `lip_response_mode2` is the opt-in second, inward-striking mode of the same lip ([`super::LipMode2`]). Cullen (2000, PhD,
//! University of Edinburgh, Appendix A Table A.2) prints it beside the first for twenty embouchures of an artificial lip pair, with a
//! NEGATIVE `1/mu` because it moves against the force, and Newton, Campbell and Gilbert (JASA 123(1) EL14-EL20, 2008) report the same
//! outward/inward pair in human and artificial lips, the played note between them, and conclude that at least two degrees of
//! freedom are needed to simulate the lip realistically. The aperture is the sum of the two modes' displacements under the same
//! force. It needs `lip_response` (the face area is shared), not `lip`. Its published numbers carry no pre-stress law, so the
//! linear loop and `bore nonlinear` ([`super::nonlinear::NonlinearLip::with_lip_mode2`], the same integrator) both refuse any lip tension
//! but 1.0 by name. A file without the record is bit-identical to one written before it existed.
//!
//! `lip_area_exponent` is the opt-in exponent `q` of the flow aperture, `w h (h/H0)^(q-1)` with `H0` the lip record's rest gap
//! ([`super::BrassVoice::with_lip_area_exponent`]). Copley and Strong (1996) describe the lip area against height as almost linear and
//! print no exponent, so `q` is a stated bracket (`1 <= q <= 4`) and the normalisation at `H0`, which leaves the aperture at the rest
//! gap unchanged, is ours. Only the flow area changes, not the lip's motion; `bore nonlinear` takes the same record
//! ([`super::nonlinear::NonlinearLip::with_lip_area_exponent`]). A file without the record, or with `q` exactly 1, is bit-identical
//! to one written before it existed.
//!
//! `lip_rim_stop` is the opt-in upper bound of the lip gap ([`super::BrassVoice::with_lip_rim_stop`]): the gap cannot open wider than the
//! mouthpiece rim, and the bound is not an argument but twice the entry radius of the first bore cone, a geometry the card already
//! holds. At the bound the first lip mode's outward momentum is removed (`inelastic`) or reversed (`elastic`): the two limit cases of a
//! restitution no source measures, reported side by side and neither called measured. Linear loop only; refused with `bore nonlinear`
//! and with `lip_response_mode2`. A file without the record is bit-identical to one written before it existed.
//!
//! `player_tract` puts a player's tract upstream of the lip ([`super::BrassVoice::with_player_tract`], and [`super::player_tract`]
//! for what it is): the rows are the reed side's, word for word. The tract is a Thevenin source in series with the aperture; with
//! `mouth` the airway's DC resistance comes out at DC only through a one-pole of time constant `tau_s` on the flow (a numerical
//! regulation constant with no default, not a measurement) and the voice refuses a tract whose mouth-pressure load is active on
//! any bin; with `lung` the card's pressure is the source at the far end. The airway takes the bore gas and its own default
//! convergence ladder. Linear loop only: refused with `bore nonlinear` and `jet_noise`; the `constriction_noise` joins the tract as a series source. The rows may stand in
//! any order; `player_tract_section` rows keep their file order. A file without the rows is bit-identical to one written before.
//!
//! A `bessel` record is a Bessel horn `r(s) = b (u - s)^(-flare)` through the
//! stated end radii over the stated length, expanded into `pieces` cones whose
//! end radii lie exactly on that curve (the piecewise-conic image every
//! transfer-matrix brass model uses; `pieces` is the refinement axis).
//!
//! The optional `bore nonlinear` record replaces the lines with the shared
//! nonlinear gas bore ([`super::nonlinear`], the `nonlinear-brass` feature):
//! the ambient gas's own NASA-9 laws evolve in cells no longer than
//! `max_cell_m`, with Zwikker-Kosten wall memory admitted on the stated band
//! and tolerance and a baffled-piston bell. It needs an `ambient
//! thermochemical` record and its packs. The mechanics clock is
//! `mechanics_per_output` (1..=16) times the output clock, each mechanics
//! sample taking `gas_steps_per_mechanics` fixed gas steps; the microphone
//! receives every mechanics sample's bell volume flow.
//!
//! Limits: one held bore, a stationary infinite rigid baffle, and the original
//! fs-duct radiation-band ceiling at the declared clock's Nyquist frequency.
//! The legacy receiver sees total flow only; `modal_baffle` resolves the retained
//! axisymmetric mouth-mode directivity. Neither resolves omitted transverse
//! modes, finite bell/body scattering or a realistic instrument; wall vibration
//! is resolved only by the opt-in `bore_wall thin-brass` row on `bore nonlinear`
//! (its breathing mode, [`super::bore_wall`]), and by nothing else.

use fs_blake3::{ContentHash, hash_domain};
use fs_duct::{Duct, LossModel, Segment, Termination};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{
    GasRecord, ThermochemicalSources, read_frozen_thermochemical_model, read_gas_record,
};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_math::det;
use fs_scenario::gesture::{GestureSchedule, GestureTarget};
use fs_vfit::discretize::DelayedFilter;

use super::bore_wall::{BELL_MAX_AZIMUTHS, BELL_MIN_AZIMUTHS, BellShellRecord, BoreWallRecord};
use super::player_tract::TractPressure;
use super::{BrassBlockDiag, BrassControl, BrassVoice, LipIsland, LipMode2};
use crate::acoustic_realize::AcousticRealizeError;
use crate::upstream_airway::{AirwayChain, AirwaySection, WEIBEL_SAAB_2021_TABLE_1_B_MM};
use crate::brass_sequence::mouth::{MAX_MOUTH_FLOW_TAPS, MouthFlowConfig, modal_flow_lines, plane_flow_line};
use crate::driving_point::{LineLadder, PassivityGain};
use crate::mm_line::{MmComboReceipt, MmLineConfig, MmLoad};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::pcm_wav::baffled::modal::{MAX_MODAL_KERNEL_TERMS, ModalOutletOptions, ModalOutletReport};
use crate::pcm_wav::baffled::openings::CircularOpening;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::turbulence_noise::{ConstrictionNoiseCard, JetNoiseCard};

/// First line of every lip performance file.
pub const LIP_PERFORMANCE_SCHEMA: &str = "frankensim-lip-performance-v1";
/// Largest accepted file.
pub const MAX_LIP_PERFORMANCE_BYTES: usize = 1 << 20;
const MAX_SEGMENTS: usize = 512;
const MAX_BESSEL_PIECES: usize = 256;
const MAX_MODES: u64 = 16;
const MAX_CONTROL_TICKS: u64 = 262_144;
const MAX_EVENTS: usize = 16_384;
const MAX_SECONDS: u64 = 600;
const MAX_BLOCK: usize = 65_536;
const MIN_RATE_HZ: u32 = 8_000;
const MAX_RATE_HZ: u32 = 192_000;

/// The lip's lumped numbers before provenance is attached.
#[derive(Debug, Clone, Copy)]
struct LipCore {
    mass_kg: f64,
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    width_m: f64,
    rest_gap_m: f64,
    face_area_m2: f64,
}

/// A `bore nonlinear` record: how the shared nonlinear gas bore is cut,
/// clocked and walled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearBoreRecord {
    /// Largest axial gas cell [m].
    pub maximum_cell_m: f64,
    /// Mechanics (microphone) samples per output sample, 1..=16.
    pub mechanics_per_output: u32,
    /// Fixed gas steps per mechanics sample.
    pub gas_steps_per_mechanics: u32,
    /// Zwikker-Kosten diffusive poles per wall field.
    pub wall_poles: usize,
    /// Band on which the wall memory is admitted [Hz].
    pub wall_band_hz: [f64; 2],
    /// Largest sampled relative error of either wall operator on that band.
    pub wall_tolerance: f64,
}

/// A finite upstream gas volume and its feed opening, before source admission.
///
/// The blowing-pressure track controls the feed, not the stored pressure.
/// This record requires the original thermochemical nonlinear bore.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteMouthRecord {
    /// Fixed, acoustically compact storage volume [m^3].
    pub volume_m3: f64,
    /// Effective source-owned nozzle feed area [m^2]; zero seals the supply.
    pub feed_area_m2: f64,
    /// Initial absolute storage pressure [Pa].
    pub initial_pressure_pa: f64,
    /// Initial storage temperature [K], admitted by the bore's source at build.
    pub initial_temperature_k: f64,
}

/// A `mouth_load modal` record: the nonlinear bore's radiating end is the plane-mode equivalent of the complete flanged
/// multimodal mouth ([`fs_duct::modal::plane_mouth_kernel`]) in place of the baffled piston. The bore's gas is the ambient gas,
/// the exterior the same room, the kernel's step the bore's own gas step. Lossless modal propagation: the nonlinear bore
/// carries its own wall memory, and the equivalent load needs only the modal flow ratios at the mouth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MouthLoadRecord {
    /// Modes of the multimodal solve, `1..=`[`fs_duct::modal::MAX_MODES`] (8 modes are cut on up to ka = 24, the kernel's evaluation limit).
    pub n_modes: usize,
    /// Extra staircase slices per cone of that solve.
    pub extra_slices: usize,
    /// Top of the band over which the realized load is claimed to match the equivalent load [Hz].
    pub band_hz: f64,
    /// Requested frequency spacing of the kernel's evaluation grid [Hz].
    pub bin_hz: f64,
    /// Largest fraction of the kernel's energy dropped from its tail.
    pub tail_tolerance: f64,
}

/// What a built `mouth_load modal` kernel is, for the receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MouthLoadInfo {
    /// The declared record.
    pub record: MouthLoadRecord,
    /// Taps of the kept kernel.
    pub taps: usize,
    /// The multimodal solve's staircase mouth radius [m].
    pub mouth_radius_m: f64,
    /// `rho c / (pi a^2)` at that radius [Pa s/m^3].
    pub z_inf: f64,
    /// Fraction of the kernel's energy in its dropped tail.
    pub dropped_tail_energy: f64,
    /// Sum of the dropped weights.
    pub dropped_sum: f64,
    /// Frequency spacing of the evaluation grid actually used [Hz].
    pub bin_hz: f64,
    /// Top of the evaluated band [Hz].
    pub evaluated_hz: f64,
    /// Smallest `Re Z / Z_inf` of the REALIZED load on that grid: a value below zero is gain at the mouth.
    pub minimum_resistance_ratio: f64,
}

/// The smallest realized resistance over the kernel's own evaluation grid, as a fraction of `Z_inf`.
#[cfg(feature = "nonlinear-brass")]
fn minimum_resistance_ratio(kernel: &fs_duct::modal::MouthKernel) -> f64 {
    let top = (kernel.evaluated_hz / kernel.bin_hz).round() as usize;
    (1..=top).fold(f64::INFINITY, |worst, k| {
        let omega = core::f64::consts::TAU * k as f64 * kernel.bin_hz;
        worst.min(kernel.impedance_at(omega).re / kernel.z_inf)
    })
}

/// Smallest realized resistance a `mouth_load modal` kernel may have, as a fraction of `Z_inf` (the equivalent load is not
/// certified passive; the commit that added it measured +3.1e-3 as the smallest on the tenor trombone to 3 kHz).
#[cfg(feature = "nonlinear-brass")]
const MOUTH_LOAD_MINIMUM_RESISTANCE_RATIO: f64 = -1.0e-3;

/// `player_tract_section` rows: a measured area function, its glottis among them, with room to spare (the reed side's bound).
const MAX_PLAYER_TRACT_SECTIONS: usize = 256;

const NONLINEAR_FEATURE: &str = "bore nonlinear needs fs-couple's nonlinear-brass feature (build music_render with --features nonlinear-brass)";

/// A parsed, validated, not yet realized performance.
#[derive(Debug, Clone)]
pub struct LipDescription {
    input_hash: ContentHash,
    rate_hz: u32,
    samples: u64,
    full_scale_pa: f64,
    gas: GasRecord,
    /// The optional `bore_gas` record: the gas the bore, cup and lip hold. `None` is the ambient gas.
    bore_gas: Option<GasRecord>,
    lip: LipIsland,
    lip_source: String,
    bore_source: String,
    cup_volume_m3: f64,
    termination: Termination,
    modes: usize,
    extra_slices: usize,
    duct: Duct,
    microphone: CircularOutletReceiver,
    observer_gas: Option<GasPropagation>,
    mouth_flow_tail_tol: f64,
    mouth_flow_max_taps: usize,
    modal_baffle: Option<ModalOutletOptions>,
    schedule: GestureSchedule,
    blow_track: String,
    tension_track: Option<String>,
    nonlinear: Option<(NonlinearBoreRecord, ThermochemicalGasModel)>,
    mouth_reservoir: Option<FiniteMouthRecord>,
    mouth_load: Option<MouthLoadRecord>,
    jet_noise: Option<JetNoiseCard>,
    constriction_noise: Option<ConstrictionNoiseCard>,
    lip_mode2: Option<LipMode2>,
    lip_area_exponent: Option<f64>,
    lip_rim_stop: Option<super::RimStopKind>,
    /// The optional breathing bore wall (`bore_wall thin-brass`) and the provenance sentence that names whose thickness and quality factor it carries.
    bore_wall: Option<(BoreWallRecord, String)>,
    /// The optional modal bell wall (`bore_wall shell-brass`) and its provenance sentence.
    bell_shell: Option<(BellShellRecord, String)>,
    /// The optional player's tract upstream of the lip: its airway chain (lip side first) and where the card's pressure acts.
    player_tract: Option<(AirwayChain, TractPressure)>,
}

struct Fields<'a> {
    words: core::str::SplitAsciiWhitespace<'a>,
    line: usize,
    key: &'a str,
}

impl<'a> Fields<'a> {
    fn err<T>(&self, what: &str) -> Result<T, String> {
        Err(format!("line {}: {}: {what}", self.line, self.key))
    }

    fn word(&mut self) -> Result<&'a str, String> {
        match self.words.next() {
            Some(word) => Ok(word),
            None => self.err("missing field"),
        }
    }

    fn parse<T: core::str::FromStr>(&mut self) -> Result<T, String> {
        let word = self.word()?;
        match word.parse::<T>() {
            Ok(value) => Ok(value),
            Err(_) => self.err("malformed field"),
        }
    }

    fn number(&mut self) -> Result<f64, String> {
        let x: f64 = self.parse()?;
        if x.is_finite() {
            Ok(x)
        } else {
            self.err("physical scalars must be finite")
        }
    }

    fn positive(&mut self) -> Result<f64, String> {
        let x = self.number()?;
        if x > 0.0 {
            Ok(x)
        } else {
            self.err("expected a positive number")
        }
    }

    fn count(&mut self, max: u64) -> Result<u64, String> {
        let n: u64 = self.parse()?;
        if (1..=max).contains(&n) {
            Ok(n)
        } else {
            self.err("count outside its admitted range")
        }
    }

    fn finish(mut self) -> Result<(), String> {
        if self.words.next().is_some() {
            return self.err("unexpected extra field");
        }
        Ok(())
    }
}

fn once<T>(slot: &mut Option<T>, value: T, f: &Fields<'_>) -> Result<(), String> {
    if slot.is_some() {
        return f.err("may only appear once");
    }
    *slot = Some(value);
    Ok(())
}

fn text_of(f: &Fields<'_>, rest: &str) -> Result<String, String> {
    if rest.is_empty() {
        return f.err("needs a provenance sentence");
    }
    Ok(rest.to_string())
}

/// Expand a Bessel horn `r(s) = r_in (u / (u - s))^flare` through `r_out`
/// at `s = length` into `pieces` cones with end radii on that curve.
fn bessel_horn(
    r_in: f64,
    r_out: f64,
    length: f64,
    flare: f64,
    pieces: usize,
    out: &mut Vec<Segment>,
) -> Result<(), String> {
    if !(r_out > r_in) || !(flare > 0.0) || !(1..=MAX_BESSEL_PIECES).contains(&pieces) {
        return Err(
            "a Bessel horn needs a growing radius, a positive flare and 1..=256 pieces".to_string(),
        );
    }
    // r_out / r_in = (u / (u - length))^flare  =>  u / (u - length) = q.
    let q = det::pow(r_out / r_in, 1.0 / flare);
    if !(q.is_finite() && q > 1.0) {
        return Err("Bessel horn end radii and flare admit no horn".to_string());
    }
    let u = q * length / (q - 1.0);
    let step = length / pieces as f64;
    let mut previous = r_in;
    for k in 1..=pieces {
        let radius = if k == pieces {
            r_out
        } else {
            r_in * det::pow(u / (u - step * k as f64), flare)
        };
        if !(radius.is_finite() && radius > previous) {
            return Err("Bessel horn pieces must have strictly growing radius".to_string());
        }
        out.push(Segment::Cone {
            inlet_radius: previous,
            outlet_radius: radius,
            length: step,
        });
        previous = radius;
    }
    Ok(())
}

fn decode_schedule(text: &str) -> Result<GestureSchedule, String> {
    // The shared decoder reserves from declared counts, including trailing ones.
    let lines = text.lines().count();
    for line in text.lines() {
        if let Some(value) = line.strip_prefix("tracks\t") {
            if !matches!(value.parse::<usize>(), Ok(1 | 2)) {
                return Err("one blowing-pressure track and an optional pre-stress track".into());
            }
        }
        if let Some(value) = line.strip_prefix("events\t") {
            let count = value
                .parse::<usize>()
                .map_err(|_| "invalid gesture event count".to_string())?;
            if count > MAX_EVENTS || count > lines {
                return Err("gesture count exceeds the source/16384-event budget".into());
            }
        }
    }
    let schedule = GestureSchedule::from_canonical_bytes(text.as_bytes())
        .map_err(|e| format!("gesture schedule refused: {e:?}"))?;
    if schedule.to_canonical_bytes() != text.as_bytes() {
        return Err(
            "gesture bytes must be canonical without ignored fields or trailing records".into(),
        );
    }
    Ok(schedule)
}

/// The control tick whose value holds at mechanical sample `sample`: tick `k`
/// applies at `ceil(k * rate / control_rate)`, the clock the pressure
/// scheduler of the reed files uses.
fn tick_at(sample: u64, control_rate_hz: u64, rate_hz: u64) -> u64 {
    sample * control_rate_hz / rate_hz
}

/// First mechanical sample of control tick `tick`.
fn tick_start(tick: u64, control_rate_hz: u64, rate_hz: u64) -> u64 {
    (tick * rate_hz).div_ceil(control_rate_hz)
}

impl LipDescription {
    /// Decode and admit every record before any physics is realized.
    ///
    /// # Errors
    /// A named refusal for malformed, repeated, missing, unphysical or
    /// out-of-budget records.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        Self::from_bytes_with_gas_sources(bytes, None)
    }

    /// Admit the same shared gas records as ducts, reeds and flue sources.
    /// NASA-9 records require explicitly associated species/transport packs;
    /// their receipts stay attached to the exact state used by every consumer.
    /// The gas is frozen throughout this performance, including the exterior.
    ///
    /// # Errors
    /// Ordinary file refusals plus missing, ambiguous or out-of-domain gas
    /// sources. No material properties are inferred from a species identifier.
    #[allow(clippy::too_many_lines)] // one record switch
    pub fn from_bytes_with_gas_sources(
        bytes: &[u8],
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, String> {
        if bytes.len() > MAX_LIP_PERFORMANCE_BYTES {
            return Err("source exceeds 1 MiB".into());
        }
        let text = core::str::from_utf8(bytes).map_err(|_| "source must be UTF-8".to_string())?;
        let (header, schedule_text) = text.split_once("\nschedule\n").ok_or_else(|| {
            "expected an LF-delimited schedule record and canonical gesture bytes".to_string()
        })?;
        let mut lines = header.lines().enumerate();
        match lines.next() {
            Some((_, first)) if first.trim_end() == LIP_PERFORMANCE_SCHEMA => {}
            _ => return Err(format!("expected the {LIP_PERFORMANCE_SCHEMA} header")),
        }
        let mut audio: Option<(u32, u64)> = None;
        let mut full_scale: Option<f64> = None;
        let mut gas: Option<GasRecord> = None;
        let mut bore_gas: Option<GasRecord> = None;
        let mut lip_core: Option<LipCore> = None;
        let mut lip_from_response = false;
        let mut lip_mode2_raw: Option<(f64, f64, f64)> = None;
        let mut lip_area_exponent_raw: Option<f64> = None;
        let mut lip_rim_stop_raw: Option<super::RimStopKind> = None;
        let mut bore_wall_raw: Option<Option<BoreWallRecord>> = None;
        let mut bell_shell_raw: Option<BellShellRecord> = None;
        let mut bore_wall_source: Option<String> = None;
        let mut player_tract_raw: Option<TractPressure> = None;
        let mut tract_sections: Vec<(f64, f64)> = Vec::new();
        let mut tract_wall: Option<fs_phs::WallPin> = None;
        let mut tract_generations: Option<usize> = None;
        let mut lip_source: Option<String> = None;
        let mut bore_source: Option<String> = None;
        let mut cup: Option<f64> = None;
        let mut termination: Option<Termination> = None;
        let mut modes: Option<(usize, usize)> = None;
        let mut microphone: Option<CircularOutletReceiver> = None;
        let mut observer_gas_record = None;
        let mut mouth_flow_tail_tol: Option<f64> = None;
        let mut mouth_flow_max_taps: Option<usize> = None;
        let mut modal_baffle: Option<ModalOutletOptions> = None;
        let mut nonlinear: Option<NonlinearBoreRecord> = None;
        let mut mouth_reservoir: Option<FiniteMouthRecord> = None;
        let mut mouth_load: Option<MouthLoadRecord> = None;
        let mut jet_noise: Option<JetNoiseCard> = None;
        let mut constriction_noise: Option<ConstrictionNoiseCard> = None;
        let mut ambient_words = None;
        let mut bore_gas_words = None;
        let mut segments: Vec<Segment> = Vec::new();
        for (index, raw) in lines {
            let line = index + 1;
            let trimmed = raw.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let (key, rest) = match trimmed.split_once(char::is_whitespace) {
                Some((key, rest)) => (key, rest.trim()),
                None => (trimmed, ""),
            };
            let mut f = Fields {
                words: rest.split_ascii_whitespace(),
                line,
                key,
            };
            match key {
                "audio" => {
                    let rate: u32 = f.parse()?;
                    let samples: u64 = f.parse()?;
                    if !(MIN_RATE_HZ..=MAX_RATE_HZ).contains(&rate)
                        || samples == 0
                        || samples > u64::from(rate) * MAX_SECONDS
                    {
                        return f.err("needs a rate in 8000..=192000 Hz and 1..=600 s of samples");
                    }
                    once(&mut audio, (rate, samples), &f)?;
                    f.finish()?;
                }
                "full_scale_pa" => {
                    let value = f.positive()?;
                    once(&mut full_scale, value, &f)?;
                    f.finish()?;
                }
                "ambient" => {
                    if f.words.clone().next() == Some("moist_air") {
                        f.words.next();
                    }
                    ambient_words = Some(f.words.clone());
                    let record = read_gas_record(&mut f.words, None, sources)
                        .map_err(|e| format!("line {line}: ambient gas refused: {e}"))?;
                    once(&mut gas, record, &f)?;
                    f.finish()?;
                }
                "bore_gas" => {
                    // The grammar of `ambient`, including its optional legacy keyword. A thermochemical bore gas is the record a nonlinear bore evolves
                    // (beside `mouth_load`, below); its words are kept for that, as the ambient's are.
                    if f.words.clone().next() == Some("moist_air") {
                        f.words.next();
                    }
                    bore_gas_words = Some(f.words.clone());
                    let record = read_gas_record(&mut f.words, None, sources)
                        .map_err(|e| format!("line {line}: bore_gas refused: {e}"))?;
                    once(&mut bore_gas, record, &f)?;
                    f.finish()?;
                }
                "lip_response" => {
                    let (freq, q, inverse_mu) = (f.positive()?, f.positive()?, f.positive()?);
                    let (width, gap, area) = (f.positive()?, f.positive()?, f.positive()?);
                    let mass = area / inverse_mu;
                    let omega = core::f64::consts::TAU * freq;
                    let core = LipCore {
                        mass_kg: mass,
                        stiffness_n_m: mass * omega * omega,
                        damping_n_s_m: mass * omega / q,
                        width_m: width,
                        rest_gap_m: gap,
                        face_area_m2: area,
                    };
                    once(&mut lip_core, core, &f)?;
                    lip_from_response = true;
                    f.finish()?;
                }
                "lip_response_mode2" => {
                    let (freq, q) = (f.positive()?, f.positive()?);
                    let inverse_mu = f.number()?;
                    if inverse_mu == 0.0 {
                        return f.err("inverse mu must be nonzero (negative for an inward-striking mode)");
                    }
                    once(&mut lip_mode2_raw, (freq, q, inverse_mu), &f)?;
                    f.finish()?;
                }
                "lip_area_exponent" => {
                    let q = f.positive()?;
                    once(&mut lip_area_exponent_raw, q, &f)?;
                    f.finish()?;
                }
                "lip_rim_stop" => {
                    let kind = match f.word()? {
                        "inelastic" => super::RimStopKind::Inelastic,
                        "elastic" => super::RimStopKind::Elastic,
                        _ => return f.err("expected inelastic or elastic (the two limit cases of an unmeasured restitution)"),
                    };
                    once(&mut lip_rim_stop_raw, kind, &f)?;
                    f.finish()?;
                }
                "bore_wall" => {
                    let kind = f.word()?;
                    if kind == "shell-brass" {
                        let thickness_m = f.positive()?;
                        let loss_factor = f.positive()?;
                        let max_frequency_hz = f.positive()?;
                        let azimuths = f.count(BELL_MAX_AZIMUTHS as u64)? as usize;
                        if azimuths < BELL_MIN_AZIMUTHS || azimuths % 2 != 0 {
                            return f.err(&format!(
                                "the bell mesh needs an even azimuth count from {BELL_MIN_AZIMUTHS} to {BELL_MAX_AZIMUTHS}"
                            ));
                        }
                        if bore_wall_raw.is_some() {
                            return f.err("may only appear once");
                        }
                        once(
                            &mut bell_shell_raw,
                            BellShellRecord { thickness_m, loss_factor, max_frequency_hz, azimuths },
                            &f,
                        )?;
                    } else {
                        let record = match kind {
                            "rigid" => None,
                            "thin-brass" => {
                                let thickness_m = f.positive()?;
                                let quality = f.positive()?;
                                Some(BoreWallRecord { thickness_m, quality })
                            }
                            _ => {
                                return f.err(
                                    "expected `rigid`, `thin-brass <thickness_m> <quality>` or `shell-brass <thickness_m> <loss_factor> <max_frequency_hz> <azimuths>`",
                                );
                            }
                        };
                        if bell_shell_raw.is_some() {
                            return f.err("may only appear once");
                        }
                        once(&mut bore_wall_raw, record, &f)?;
                    }
                    f.finish()?;
                }
                "bore_wall_source" => {
                    let value = text_of(&f, rest)?;
                    once(&mut bore_wall_source, value, &f)?;
                }
                "player_tract" => {
                    let pressure = match f.word()? {
                        "mouth" => TractPressure::Mouth { tau_s: f.positive()? },
                        "lung" => TractPressure::Lung,
                        _ => return f.err("expected `mouth <tau_s>` or `lung`"),
                    };
                    once(&mut player_tract_raw, pressure, &f)?;
                    f.finish()?;
                }
                "player_tract_section" => {
                    let (length_m, area_m2) = (f.positive()?, f.positive()?);
                    if tract_sections.len() == MAX_PLAYER_TRACT_SECTIONS {
                        return f.err("a player tract takes at most 256 sections");
                    }
                    tract_sections.push((length_m, area_m2));
                    f.finish()?;
                }
                "player_tract_wall" => {
                    let surface_density = f.positive()?;
                    let (stiffness_per_area, resistance) = (f.number()?, f.number()?);
                    if stiffness_per_area < 0.0 || resistance < 0.0 {
                        return f.err("needs a nonnegative stiffness and resistance");
                    }
                    let pin = fs_phs::WallPin { surface_density, stiffness_per_area, resistance };
                    once(&mut tract_wall, pin, &f)?;
                    f.finish()?;
                }
                "player_airway" => {
                    if f.word()? != "weibel-1963" || f.word()? != "saab-2021-table-1-B" {
                        return f.err("takes `weibel-1963 saab-2021-table-1-B <generations>`");
                    }
                    let generations = f.count(WEIBEL_SAAB_2021_TABLE_1_B_MM.len() as u64)? as usize;
                    once(&mut tract_generations, generations, &f)?;
                    f.finish()?;
                }
                "lip" => {
                    let (mass, stiffness) = (f.positive()?, f.positive()?);
                    let damping = f.number()?;
                    let (width, gap, area) = (f.positive()?, f.positive()?, f.positive()?);
                    if damping < 0.0 {
                        return f.err("damping must be nonnegative");
                    }
                    let core = LipCore {
                        mass_kg: mass,
                        stiffness_n_m: stiffness,
                        damping_n_s_m: damping,
                        width_m: width,
                        rest_gap_m: gap,
                        face_area_m2: area,
                    };
                    once(&mut lip_core, core, &f)?;
                    f.finish()?;
                }
                "lip_source" => {
                    let value = text_of(&f, rest)?;
                    once(&mut lip_source, value, &f)?;
                }
                "bore_source" => {
                    let value = text_of(&f, rest)?;
                    once(&mut bore_source, value, &f)?;
                }
                "cup_volume_m3" => {
                    let value = f.positive()?;
                    once(&mut cup, value, &f)?;
                    f.finish()?;
                }
                "mouth_reservoir" => {
                    let volume_m3 = f.positive()?;
                    let feed_area_m2 = f.number()?;
                    if feed_area_m2 < 0.0 {
                        return f.err("feed area must be non-negative; zero seals the supply");
                    }
                    let initial_pressure_pa = f.positive()?;
                    let initial_temperature_k = f.positive()?;
                    once(&mut mouth_reservoir, FiniteMouthRecord {
                        volume_m3,
                        feed_area_m2,
                        initial_pressure_pa,
                        initial_temperature_k,
                    }, &f)?;
                    f.finish()?;
                }
                "termination" => {
                    let kind = match f.word()? {
                        "flanged" => Termination::FlangedOpen,
                        "unflanged" => Termination::UnflangedOpen,
                        "ideal_open" => Termination::IdealOpen,
                        _ => return f.err("expected flanged, unflanged or ideal_open"),
                    };
                    once(&mut termination, kind, &f)?;
                    f.finish()?;
                }
                "modes" => {
                    let n_modes = f.count(MAX_MODES)? as usize;
                    let extra = f.count(MAX_MODES)? as usize;
                    once(&mut modes, (n_modes, extra), &f)?;
                    f.finish()?;
                }
                "cylinder" => {
                    let (radius, length) = (f.positive()?, f.positive()?);
                    segments.push(Segment::Cylinder { radius, length });
                    f.finish()?;
                }
                "cone" => {
                    let (inlet, outlet, length) = (f.positive()?, f.positive()?, f.positive()?);
                    segments.push(if inlet == outlet {
                        Segment::Cylinder {
                            radius: inlet,
                            length,
                        }
                    } else {
                        Segment::Cone {
                            inlet_radius: inlet,
                            outlet_radius: outlet,
                            length,
                        }
                    });
                    f.finish()?;
                }
                "bessel" => {
                    let (inlet, outlet, length) = (f.positive()?, f.positive()?, f.positive()?);
                    let flare = f.positive()?;
                    let pieces = f.count(MAX_BESSEL_PIECES as u64)? as usize;
                    bessel_horn(inlet, outlet, length, flare, pieces, &mut segments)
                        .map_err(|e| format!("line {line}: bessel: {e}"))?;
                    f.finish()?;
                }
                "bore" => {
                    if f.word()? != "nonlinear" {
                        return f.err("expected nonlinear");
                    }
                    let maximum_cell_m = f.positive()?;
                    let mechanics_per_output = f.count(16)? as u32;
                    let gas_steps_per_mechanics = f.count(64)? as u32;
                    let wall_poles = f.count(128)? as usize;
                    let wall_band_hz = [f.positive()?, f.positive()?];
                    let wall_tolerance = f.positive()?;
                    if wall_band_hz[0] >= wall_band_hz[1] || wall_tolerance >= 1.0 {
                        return f
                            .err("needs a rising wall band and a relative wall tolerance below 1");
                    }
                    let record = NonlinearBoreRecord {
                        maximum_cell_m,
                        mechanics_per_output,
                        gas_steps_per_mechanics,
                        wall_poles,
                        wall_band_hz,
                        wall_tolerance,
                    };
                    once(&mut nonlinear, record, &f)?;
                    f.finish()?;
                }
                "mouth_flow_tail_tol" => {
                    let value = f.positive()?;
                    if value >= 1.0 {
                        return f.err("a relative tail energy must be below 1");
                    }
                    once(&mut mouth_flow_tail_tol, value, &f)?;
                    f.finish()?;
                }
                "mouth_flow_max_taps" => {
                    let taps = f.count(MAX_MOUTH_FLOW_TAPS as u64)? as usize;
                    if taps < 512 || !taps.is_power_of_two() {
                        return f.err("mouth-flow taps must be a power of two in 512..=262144");
                    }
                    once(&mut mouth_flow_max_taps, taps, &f)?;
                    f.finish()?;
                }
                "modal_baffle" => {
                    let relative_kernel_tolerance = f.positive()?;
                    if relative_kernel_tolerance >= 1.0 {
                        return f.err("modal quadrature tolerance must be below 1");
                    }
                    let maximum_kernel_terms = f.count(MAX_MODAL_KERNEL_TERMS as u64)? as usize;
                    once(&mut modal_baffle, ModalOutletOptions {
                        relative_kernel_tolerance, maximum_kernel_terms,
                    }, &f)?;
                    f.finish()?;
                }
                "mouth_load" => {
                    if f.word()? != "modal" {
                        return f.err("expected modal");
                    }
                    let n_modes = f.count(fs_duct::modal::MAX_MODES as u64)? as usize;
                    let extra_slices = f.count(64)? as usize;
                    let (band_hz, bin_hz, tail_tolerance) = (f.positive()?, f.positive()?, f.positive()?);
                    if !(band_hz > 2.0 * bin_hz) || tail_tolerance >= 1.0 {
                        return f.err("needs a band of at least two bins and a tail tolerance below 1");
                    }
                    let record = MouthLoadRecord { n_modes, extra_slices, band_hz, bin_hz, tail_tolerance };
                    once(&mut mouth_load, record, &f)?;
                    f.finish()?;
                }
                "jet_noise" => {
                    // Every field is required: the gain, the band and the onset are a source's measured numbers (Verge 1995,
                    // a flue jet), and the card, not this code, says which it uses and where they come from.
                    let seed: u64 = f.parse()?;
                    let card = JetNoiseCard {
                        seed,
                        gain: f.number()?,
                        strouhal_lo: f.number()?,
                        strouhal_hi: f.number()?,
                        reynolds_onset: f.number()?,
                    };
                    if let Err(why) = card.validate() {
                        return f.err(why);
                    }
                    once(&mut jet_noise, card, &f)?;
                    f.finish()?;
                }
                "constriction_noise" => {
                    // Two fields and no gain: the law is Oldham and Ukpoho's printed trend, with no free number to declare.
                    // The radius is the cup rim radius the bore record's first station gives; the card, not this code, says it.
                    let seed: u64 = f.parse()?;
                    let card = ConstrictionNoiseCard { seed, duct_radius_m: f.positive()? };
                    if let Err(why) = card.validate() {
                        return f.err(why);
                    }
                    once(&mut constriction_noise, card, &f)?;
                    f.finish()?;
                }
                "microphone" => {
                    let position = [f.number()?, f.number()?, f.number()?];
                    let rings = f.count(128)? as usize;
                    let points = f.count(512)? as usize;
                    let band = f.positive()?;
                    if position[2] <= 0.0 {
                        return f.err("the microphone must stand in front of the baffle (z > 0)");
                    }
                    let receiver = CircularOutletReceiver {
                        position_m: position,
                        radial_rings: rings,
                        angular_points: points,
                        maximum_frequency_hz: band,
                    };
                    once(&mut microphone, receiver, &f)?;
                    f.finish()?;
                }
                "observer_gas" => {
                    // Bind after all records so molecular laws use the actual
                    // ambient composition regardless of declaration order.
                    once(&mut observer_gas_record, (line, rest), &f)?;
                }
                _ => return Err(format!("line {line}: unknown key {key:?}")),
            }
            if segments.len() > MAX_SEGMENTS {
                return Err(format!("line {line}: bore exceeds {MAX_SEGMENTS} segments"));
            }
        }
        let missing = |key: &str| format!("missing required record {key}");
        let (rate_hz, samples) = audio.ok_or_else(|| missing("audio"))?;
        let core = lip_core.ok_or_else(|| missing("lip_response or lip"))?;
        let lip_source = lip_source.ok_or_else(|| missing("lip_source"))?;
        let bore_source = bore_source.ok_or_else(|| missing("bore_source"))?;
        if segments.is_empty() {
            return Err(missing("cylinder, cone or bessel"));
        }
        let lip = LipIsland {
            mass_kg: core.mass_kg,
            stiffness_n_m: core.stiffness_n_m,
            damping_n_s_m: core.damping_n_s_m,
            width_m: core.width_m,
            rest_gap_m: core.rest_gap_m,
            face_area_m2: core.face_area_m2,
            provenance: lip_source.clone(),
        };
        let (modes_n, extra_slices) = modes.ok_or_else(|| missing("modes"))?;
        let schedule = decode_schedule(schedule_text)?;
        let mut blow_track = None;
        let mut tension_track = None;
        for track in schedule.tracks() {
            let slot = match &track.target {
                GestureTarget::BlowingPressure => &mut blow_track,
                GestureTarget::PreStress => &mut tension_track,
                _ => {
                    return Err(
                        "a lip performance accepts only blowing-pressure and pre-stress tracks"
                            .into(),
                    );
                }
            };
            if slot.replace(track.id.clone()).is_some() {
                return Err("each control track may appear only once".into());
            }
        }
        let blow_track =
            blow_track.ok_or_else(|| "a blowing-pressure track is required".to_string())?;
        if mouth_reservoir.is_some() && nonlinear.is_none() {
            return Err("mouth_reservoir requires bore nonlinear and its original thermochemical source".into());
        }
        if mouth_load.is_some() && nonlinear.is_none() {
            return Err("mouth_load requires bore nonlinear (a linear bore takes modal_baffle)".into());
        }
        if modal_baffle.is_some()
            && (termination != Some(Termination::FlangedOpen) || nonlinear.is_some())
        {
            return Err("modal_baffle requires a linear bore and termination flanged; no foreign exterior load is substituted".into());
        }
        let gas = gas.ok_or_else(|| missing("ambient"))?;
        if let Some(bore) = &bore_gas {
            // A linear bore names the room in `modal_baffle`'s flanged load, a nonlinear bore in `mouth_load`'s kernel (the two-gas flanged law).
            if modal_baffle.is_none() && !(nonlinear.is_some() && mouth_load.is_some()) {
                return Err("bore_gas needs modal_baffle (or, on bore nonlinear, mouth_load): its flanged exterior load names the room gas, whereas the analytic mouth loads would radiate into the bore gas".into());
            }
            if bore.state().pressure != gas.state().pressure {
                return Err("bore_gas and ambient need one static pressure; mean flow is not modeled".into());
            }
        }
        let observer_gas = observer_gas_record.map(|(line, record)| {
            GasPropagation::read_record_for_gas(&mut record.split_ascii_whitespace(), &gas)
                .map_err(|error| format!("line {line}: observer_gas refused: {error}"))
        }).transpose()?;
        let lip_mode2 = match lip_mode2_raw {
            None => None,
            Some(_) if !lip_from_response => {
                return Err("lip_response_mode2 needs the lip_response record (its first mode and shared face area), not lip".into());
            }
            Some((freq, q, inverse_mu)) => Some(
                LipMode2::from_response(freq, q, inverse_mu, core.face_area_m2)
                    .map_err(|e| format!("lip_response_mode2 refused: {e}"))?,
            ),
        };
        let lip_area_exponent = match lip_area_exponent_raw {
            None => None,
            Some(q) if !(1.0..=4.0).contains(&q) => {
                return Err(format!("lip_area_exponent must be in [1, 4] (a stated bracket), got {q}"));
            }
            Some(q) => Some(q),
        };
        let lip_rim_stop = match lip_rim_stop_raw {
            None => None,
            Some(_) if lip_mode2.is_some() => {
                return Err("lip_rim_stop acts on the first lip mode only and is refused with lip_response_mode2".into());
            }
            Some(kind) => Some(kind),
        };
        // `bore_wall` (module doc): rigid is the row absent; thin-brass needs the nonlinear bore and a sentence naming whose numbers.
        let (bore_wall, bell_shell) = match (bore_wall_raw, bell_shell_raw, bore_wall_source) {
            (None, None, None) | (Some(None), None, None) => (None, None),
            (None, None, Some(_)) => return Err("bore_wall_source needs a bore_wall row".into()),
            (Some(None), None, Some(_)) => {
                return Err("bore_wall_source needs bore_wall thin-brass or shell-brass; bore_wall rigid carries no source".into());
            }
            (Some(Some(_)), None, None) => {
                return Err("bore_wall thin-brass needs a bore_wall_source row naming whose instrument and whose measurement the thickness and the quality factor are".into());
            }
            (None, Some(_), None) => {
                return Err("bore_wall shell-brass needs a bore_wall_source row naming whose instrument and whose measurement the thickness and the loss factor are".into());
            }
            (Some(Some(_)), None, Some(_)) if nonlinear.is_none() => {
                return Err("bore_wall thin-brass is on the nonlinear bore only: the linear mode-matching bore has rigid walls and refuses the row (add a bore nonlinear row)".into());
            }
            (None, Some(_), Some(_)) if nonlinear.is_none() => {
                return Err("bore_wall shell-brass is on the nonlinear bore only: the linear mode-matching bore has rigid walls and refuses the row (add a bore nonlinear row)".into());
            }
            (Some(Some(record)), None, Some(source)) => (Some((record, source)), None),
            (None, Some(record), Some(source)) => (None, Some((record, source))),
            (Some(_), Some(_), _) => return Err("bore_wall may only appear once".into()),
        };
        // `player_tract` and its rows (module doc): the chain is built here, so a malformed tract refuses at the parse.
        let player_tract = match player_tract_raw {
            None => {
                if !tract_sections.is_empty() || tract_wall.is_some() || tract_generations.is_some() {
                    return Err("player_tract_section, player_tract_wall and player_airway need a player_tract row".into());
                }
                None
            }
            Some(pressure) => {
                if tract_sections.is_empty() {
                    return Err("player_tract needs at least one player_tract_section row".into());
                }
                if nonlinear.is_some() {
                    return Err("player_tract is on the linear bore only and is refused with bore nonlinear".into());
                }
                if jet_noise.is_some() {
                    return Err("player_tract is not combined with jet_noise".into());
                }
                let tubes = tract_sections
                    .iter()
                    .map(|&(length_m, area_m2)| AirwaySection::from_area(length_m, area_m2, 1, tract_wall))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|_| "player_tract sections need positive finite lengths and areas".to_string())?;
                let chain = AirwayChain::new(tubes, LossModel::Bessel)
                    .map_err(|_| "player_tract sections did not form an airway".to_string())?;
                let chain = match tract_generations {
                    Some(generations) => AirwayChain::weibel_saab_2021_b(generations, tract_wall)
                        .and_then(|tree| chain.ahead_of(&tree))
                        .map_err(|_| "player_airway did not join the tract".to_string())?,
                    None => chain,
                };
                Some((chain, pressure))
            }
        };
        let nonlinear = match nonlinear {
            None => None,
            Some(_) if !cfg!(feature = "nonlinear-brass") => return Err(NONLINEAR_FEATURE.into()),
            // Its bell is a rigid piston in an infinite baffle, nothing else.
            Some(_) if termination != Some(Termination::FlangedOpen) => {
                return Err(
                    "bore nonlinear radiates as a baffled piston: it needs termination flanged"
                        .into(),
                );
            }
            Some(record) => {
                // The bore evolves its own gas's laws, never a snapshot: the `bore_gas` row's when the file declares one (the room is then only the exterior
                // of the mouth load, the receiver and the observer), else the ambient's.
                let (words, bore_row) = match bore_gas_words {
                    Some(words) => (Some(words), true),
                    None => (ambient_words, false),
                };
                let model = match (words, sources) {
                    (Some(mut words), Some(sources))
                        if words.clone().next() == Some("thermochemical") =>
                    {
                        let which = if bore_row { "bore_gas" } else { "ambient" };
                        read_frozen_thermochemical_model(&mut words, sources)
                            .map_err(|e| {
                                format!("bore nonlinear: {which} gas source refused: {e}")
                            })?
                            .0
                    }
                    _ if bore_row => {
                        return Err(
                            "bore nonlinear evolves the bore gas's own NASA-9 laws: its bore_gas row needs a thermochemical record and its packs"
                                .into(),
                        );
                    }
                    _ => {
                        return Err(
                            "bore nonlinear evolves the ambient gas's own NASA-9 laws: it needs an ambient thermochemical record and its packs"
                                .into(),
                        );
                    }
                };
                Some((record, model))
            }
        };
        let control_rate = u64::from(schedule.control_rate_hz);
        if control_rate == 0 || control_rate > u64::from(rate_hz) {
            return Err(
                "the control rate must be positive and no greater than the mechanics clock".into(),
            );
        }
        let last_tick = tick_at(samples - 1, control_rate, u64::from(rate_hz));
        if last_tick >= MAX_CONTROL_TICKS {
            return Err("control compilation exceeds the 262144 control-tick budget".into());
        }
        let last_time = last_tick as f64 / control_rate as f64;
        if schedule
            .tracks()
            .iter()
            .any(|track| track.events.iter().any(|event| event.time_s > last_time))
        {
            return Err("an authored command starts after the last observed control tick".into());
        }
        Ok(Self {
            input_hash: hash_domain("org.frankensim.fs-couple.lip-performance.v1", bytes),
            rate_hz,
            samples,
            full_scale_pa: full_scale.ok_or_else(|| missing("full_scale_pa"))?,
            gas,
            bore_gas,
            lip,
            lip_source,
            bore_source,
            cup_volume_m3: cup.ok_or_else(|| missing("cup_volume_m3"))?,
            termination: termination.ok_or_else(|| missing("termination"))?,
            modes: modes_n,
            extra_slices,
            duct: Duct { segments },
            microphone: microphone.ok_or_else(|| missing("microphone"))?,
            observer_gas,
            mouth_flow_tail_tol: mouth_flow_tail_tol.unwrap_or(LineLadder::default().tail_tol),
            mouth_flow_max_taps: mouth_flow_max_taps.unwrap_or(MouthFlowConfig::default().max_taps),
            modal_baffle,
            schedule,
            blow_track,
            tension_track,
            nonlinear,
            mouth_reservoir,
            mouth_load,
            jet_noise,
            constriction_noise,
            lip_mode2,
            lip_area_exponent,
            lip_rim_stop,
            bore_wall,
            bell_shell,
            player_tract,
        })
    }

    /// Where the opt-in player tract's card pressure acts, if the file declared one ([`super::BrassVoice::with_player_tract`]).
    #[must_use]
    pub fn player_tract_pressure(&self) -> Option<TractPressure> {
        self.player_tract.as_ref().map(|(_, pressure)| *pressure)
    }

    /// The opt-in player tract's airway chain (lip side first), if the file declared one.
    #[must_use]
    pub fn player_tract_chain(&self) -> Option<&AirwayChain> {
        self.player_tract.as_ref().map(|(chain, _)| chain)
    }

    /// The opt-in lip jet-noise record, if the file declared one ([`crate::turbulence_noise::LipJetNoise`]).
    #[must_use]
    pub const fn jet_noise(&self) -> Option<&JetNoiseCard> {
        self.jet_noise.as_ref()
    }

    /// The opt-in second lip mode, if the file declared one ([`LipMode2`]).
    #[must_use]
    pub const fn lip_mode2(&self) -> Option<&LipMode2> {
        self.lip_mode2.as_ref()
    }

    /// The opt-in flow-area exponent `q`, if the file declared one ([`super::BrassVoice::with_lip_area_exponent`]).
    #[must_use]
    pub const fn lip_area_exponent(&self) -> Option<f64> {
        self.lip_area_exponent
    }

    /// The opt-in breathing bore wall and its provenance sentence, if the file declared `bore_wall thin-brass` ([`super::bore_wall`]).
    #[must_use]
    pub fn bore_wall(&self) -> Option<&(BoreWallRecord, String)> {
        self.bore_wall.as_ref()
    }

    /// The opt-in modal bell wall and its provenance sentence, if the file declared `bore_wall shell-brass` ([`super::bell_shell`]).
    #[must_use]
    pub fn bell_shell(&self) -> Option<&(BellShellRecord, String)> {
        self.bell_shell.as_ref()
    }

    /// The opt-in rim stop's limit case, if the file declared one ([`super::BrassVoice::with_lip_rim_stop`]).
    #[must_use]
    pub const fn lip_rim_stop(&self) -> Option<super::RimStopKind> {
        self.lip_rim_stop
    }

    /// The largest lip gap the opt-in rim stop allows: the mouthpiece rim opening, twice the entry radius of the first bore cone.
    fn rim_max_gap_m(&self) -> Result<f64, String> {
        self.duct
            .segments
            .first()
            .map(|segment| 2.0 * segment.inlet_radius())
            .ok_or_else(|| "the bore has no first segment".to_string())
    }

    /// The opt-in lip constriction-noise record, if the file declared one ([`crate::turbulence_noise::ConstrictionNoise`]).
    #[must_use]
    pub const fn constriction_noise(&self) -> Option<&ConstrictionNoiseCard> {
        self.constriction_noise.as_ref()
    }

    /// Finite upstream storage and its feed, if the file declares it.
    #[must_use]
    pub const fn mouth_reservoir(&self) -> Option<&FiniteMouthRecord> {
        self.mouth_reservoir.as_ref()
    }

    /// The `mouth_load modal` record, if the file declares one: the nonlinear bore's load is the plane-mode equivalent of
    /// the complete flanged multimodal mouth.
    #[must_use]
    pub const fn mouth_load(&self) -> Option<&MouthLoadRecord> {
        self.mouth_load.as_ref()
    }

    /// The `bore nonlinear` record and the ambient source model it evolves,
    /// when the file declares one.
    #[must_use]
    pub const fn nonlinear_bore(&self) -> Option<&(NonlinearBoreRecord, ThermochemicalGasModel)> {
        self.nonlinear.as_ref()
    }

    /// The exact shared state used by the lip, bore, mouth-flow and receiver.
    #[must_use]
    pub fn ambient_gas(&self) -> &GasState {
        self.gas.state()
    }

    /// Original gas source receipts, retained without constitutive conversion.
    #[must_use]
    pub const fn gas_record(&self) -> &GasRecord {
        &self.gas
    }

    /// The gas the bore, cup and lip hold: the `bore_gas` row's record when the file declared one, otherwise the ambient
    /// gas (the room), as before the row existed.
    #[must_use]
    pub fn bore_gas(&self) -> &GasState {
        self.bore_gas.as_ref().map_or_else(|| self.gas.state(), GasRecord::state)
    }

    /// The `bore_gas` row's original gas record, retained without constitutive conversion; `None` without the row.
    #[must_use]
    pub const fn bore_gas_record(&self) -> Option<&GasRecord> {
        self.bore_gas.as_ref()
    }

    /// Explicit exterior propagation on the original bell surface paths.
    #[must_use]
    pub const fn observer_gas(&self) -> Option<&GasPropagation> {
        self.observer_gas.as_ref()
    }

    /// Retained transverse modes in each section of the bore's acoustic image.
    #[must_use]
    pub const fn transverse_mode_count(&self) -> usize {
        self.modes
    }

    /// Retained source control items, including each track's initial value.
    /// Ramps remain owned and evaluated by this source on its control clock.
    #[must_use]
    pub fn control_event_count(&self) -> usize {
        self.schedule.tracks().iter().map(|track| 1 + track.events.len()).sum()
    }

    /// The passive bore (no lip, cup shunt or observer) at `frequencies_hz`.
    /// The modal impedance uses the render's actual load, including a declared
    /// full modal baffle. Scalar impedance/flow remain comparison diagnostics.
    ///
    /// # Errors
    /// The duct owner's refusal at any frequency (loss-model or radiation validity).
    pub fn bore_scan(&self, frequencies_hz: &[f64]) -> Result<BoreScan, String> {
        // `room` is the exterior (the radiation load and the receiver); `gas` is what the bore holds (the same state
        // unless a `bore_gas` row separates them).
        let room = self.gas.state();
        let gas = self.bore_gas();
        let radius = self
            .duct
            .segments
            .first()
            .map(Segment::inlet_radius)
            .ok_or_else(|| "empty bore".to_string())?;
        let [x, y, z] = self.microphone.position_m;
        let mut rows = Vec::with_capacity(frequencies_hz.len());
        for &frequency_hz in frequencies_hz {
            let omega = core::f64::consts::TAU * frequency_hz;
            let multimodal = if self.modal_baffle.is_some() {
                fs_duct::modal::mm_input_impedance_flanged(&self.duct, gas, room, omega,
                    LossModel::Bessel, self.modes, self.extra_slices)
            } else { fs_duct::modal::mm_input_impedance(
                &self.duct,
                gas,
                omega,
                LossModel::Bessel,
                self.termination,
                self.modes,
                self.extra_slices,
            ) }
            .map_err(|e| format!("multimodal bore at {frequency_hz} Hz refused: {e}"))?;
            let scalar = fs_duct::input_impedance(
                &self.duct,
                gas,
                omega,
                LossModel::Bessel,
                self.termination,
            )
            .map_err(|e| format!("plane-wave bore at {frequency_hz} Hz refused: {e}"))?;
            rows.push(BoreResponse {
                frequency_hz,
                impedance_multimodal: [
                    multimodal.plane_impedance.re,
                    multimodal.plane_impedance.im,
                ],
                impedance_scalar: [scalar.impedance.re, scalar.impedance.im],
                mouth_over_input_flow: [scalar.u_mouth_over_u_in.re, scalar.u_mouth_over_u_in.im],
            });
        }
        Ok(BoreScan {
            input_characteristic_impedance: gas.density * gas.sound_speed
                / (core::f64::consts::PI * radius * radius),
            // The receiver's medium: the room, which the pressure at the microphone is read in.
            density_kg_m3: room.density,
            microphone_distance_m: x.hypot(y).hypot(z),
            rows,
        })
    }

    /// Realize the bore (the lines, or the nonlinear gas bore when the file
    /// declares one) and the microphone.
    ///
    /// # Errors
    /// Any owner's refusal: the multimodal bank (passivity, ladder
    /// convergence, radiation `ka` ceiling), the mouth-flow ladder, the
    /// nonlinear bore (gas source, mesh, wall binding, clock) or the
    /// receiver's resolved-band admission.
    pub fn build(&self) -> Result<LipPerformance, String> {
        let gas = self.ambient_gas();
        #[cfg(feature = "nonlinear-brass")]
        if let Some((record, model)) = &self.nonlinear {
            let mouth_radius_m = self.duct.segments.last().map(Segment::outlet_radius)
                .ok_or_else(|| "empty bore".to_string())?;
            return self.build_nonlinear(*record, model, mouth_radius_m);
        }
        let config = MmLineConfig {
            sample_rate_hz: self.rate_hz,
            n_modes: self.modes,
            extra_slices: self.extra_slices,
            loss: LossModel::Bessel,
        };
        let load = if self.modal_baffle.is_some() { MmLoad::FlangedModal(gas) }
            else { MmLoad::Analytic(self.termination) };
        // Reject an unresolved or over-budget receiver before the frequency
        // sweeps. Its aperture comes from the same immutable geometry owner.
        let modal_receiver = if let Some(options) = self.modal_baffle {
            let geometry = fs_duct::modal::PreparedModalDuct::new(&self.duct,
                self.modes, self.extra_slices).map_err(|e| e.to_string())?;
            let radius = geometry.mouth_radius_m();
            let (receiver, report, propagation) = if let Some(propagation) = &self.observer_gas {
                let (receiver, report, gas_report) =
                    BaffledPressure::circular_modal_outlet_multirate_gas_differentiating(
                        radius, fs_duct::modal::ModeCount(self.modes), self.rate_hz,
                        self.rate_hz, self.microphone, gas, options, propagation, || false,
                    ).map_err(|error| format!("observer_gas refused: {error}"))?;
                (receiver, report, Some(gas_report))
            } else {
                let (receiver, report) = BaffledPressure::circular_modal_outlet_multirate_differentiating(
                    radius, fs_duct::modal::ModeCount(self.modes), self.rate_hz,
                    self.rate_hz, self.microphone, gas, options, || false)?;
                (receiver, report, None)
            };
            Some((radius, receiver, report, propagation))
        } else { None };
        // What the bore, cup and lip hold; `gas` stays the exterior (the room) for the radiation load and the receiver.
        // Without a `bore_gas` row the two are the same state, so every operand below is the one it always was.
        let bore_state = self.bore_gas();
        let voice = BrassVoice::new(
            core::slice::from_ref(&self.duct),
            &["held"],
            bore_state,
            &load,
            &config,
            self.lip.clone(),
            self.cup_volume_m3,
        )
        .map_err(|e| format!("brass voice refused: {e}"))?;
        let voice = match self.jet_noise {
            Some(card) => voice
                .with_jet_noise(card)
                .map_err(|e| format!("jet_noise refused: {e}"))?,
            None => voice,
        };
        let voice = match self.lip_mode2 {
            Some(mode) => voice
                .with_lip_mode2(mode)
                .map_err(|e| format!("lip_response_mode2 refused: {e}"))?,
            None => voice,
        };
        let voice = match self.lip_area_exponent {
            Some(q) => voice
                .with_lip_area_exponent(q)
                .map_err(|e| format!("lip_area_exponent refused: {e}"))?,
            None => voice,
        };
        let voice = match self.lip_rim_stop {
            Some(kind) => voice
                .with_lip_rim_stop(self.rim_max_gap_m()?, kind)
                .map_err(|e| format!("lip_rim_stop refused: {e}"))?,
            None => voice,
        };
        let voice = match self.constriction_noise {
            Some(card) => voice
                .with_constriction_noise(card)
                .map_err(|e| format!("constriction_noise refused: {e}"))?,
            None => voice,
        };
        let voice = match &self.player_tract {
            Some((chain, pressure)) => voice
                .with_player_tract(chain, *pressure, bore_state, &LineLadder::default())
                .map_err(|e| format!("player_tract refused: {e}"))?,
            None => voice,
        };
        let line = voice
            .receipts()
            .first()
            .cloned()
            .ok_or_else(|| "the bank realized no line".to_string())?;
        let output = MouthFlowConfig {
            tail_tolerance: self.mouth_flow_tail_tol,
            tap_change_tolerance: self.mouth_flow_tail_tol.sqrt(),
            max_taps: self.mouth_flow_max_taps,
            ..MouthFlowConfig::default()
        };
        let (mouth, mouth_receipt) = if self.modal_baffle.is_some() {
            modal_flow_lines(&self.duct, bore_state, load, &config, voice.port_impedance(), output)
        } else {
            plane_flow_line(&self.duct, gas, load, &config, voice.port_impedance(), output)
                .map(|(filter, receipt)| (vec![filter], receipt))
        }
        .map_err(|e| format!("mouth-flow line refused: {e}"))?;
        let mouth_radius_m = mouth_receipt.mouth_radius_m;
        let (receiver, modal_radiation, observer_gas) = if let Some((radius, receiver, report, propagation)) = modal_receiver {
            if radius != mouth_radius_m {
                return Err("modal load and receiver must retain the same actual mouth radius".into());
            }
            (receiver, Some(report), propagation)
        } else {
            let (receiver, observer_gas) = self.scalar_receiver(mouth_radius_m, self.rate_hz)?;
            (receiver, None, observer_gas)
        };
        let bore = LipBoreInfo::Linear(LinearBoreInfo {
            line,
            mouth_flow_modes: mouth.len(),
            modal_radiation,
            mouth_flow_taps: mouth_receipt.n_fft,
            mouth_flow_rungs: mouth_receipt.rungs.len(),
            mouth_flow_tail_tol: self.mouth_flow_tail_tol,
            mouth_flow_ladder: mouth_receipt
                .rungs
                .iter()
                .map(|r| (r.n_fft, r.tail_energy.iter().copied().fold(0.0, f64::max),
                    r.tap_change.iter().copied().fold(0.0, f64::max)))
                .collect(),
            player_tract_passivity_gain: voice
                .player_tract()
                .and_then(|tract| tract.airway().receipt().reflectance.passivity_gain),
        });
        Ok(self.performance(
            mouth_radius_m,
            bore,
            Engine::Linear { voice, mouth },
            receiver,
            observer_gas,
            1,
        ))
    }

    #[cfg(feature = "nonlinear-brass")]
    fn build_nonlinear(
        &self,
        record: NonlinearBoreRecord,
        model: &ThermochemicalGasModel,
        mouth_radius_m: f64,
    ) -> Result<LipPerformance, String> {
        use super::nonlinear::{FiniteMouthSpec, NonlinearClock, NonlinearLip};
        use fs_duct::nonlinear::bore::{WallGas, WallMemory};
        use fs_duct::nonlinear::{CaloricDomain, FrozenGas};
        use fs_qty::Temperature;
        // What the bore, cup and lip hold. The room (the exterior of the mouth load; `scalar_receiver` reads it for the receiver and the observer) is the same state
        // unless a `bore_gas` row separates them, so every operand below is then the one it always was.
        let gas = self.bore_gas();
        let room = self.ambient_gas();
        if self.lip_rim_stop.is_some() {
            return Err("lip_rim_stop acts on the linear loop only and is refused with bore nonlinear".to_string());
        }
        let [lowest, highest] = room_region(model.thermodynamics(), gas.temperature)?;
        let source = FrozenGas::new(
            model.thermodynamics().clone(),
            CaloricDomain::new(Temperature::new(lowest), Temperature::new(highest)),
            || false,
        )
        .map_err(|e| format!("nonlinear bore gas refused: {e}"))?;
        let clock = NonlinearClock {
            maximum_cell_m: record.maximum_cell_m,
            output_rate_hz: self.rate_hz,
            mechanics_per_output: record.mechanics_per_output,
            gas_steps_per_mechanics: record.gas_steps_per_mechanics,
            wall: WallMemory {
                poles: record.wall_poles,
                band_hz: record.wall_band_hz,
                maximum_relative_error: record.wall_tolerance,
            },
        };
        let wall_gas = WallGas {
            density_kg_m3: gas.density,
            dynamic_viscosity_pa_s: gas.dynamic_viscosity,
            gamma: gas.gamma,
            prandtl: gas.prandtl,
        };
        // The bore's load: the plane-mode equivalent of the flanged multimodal mouth (the bore gas inside, the room outside), on the bore's own gas step, or the baffled piston.
        let mouth = match self.mouth_load {
            None => None,
            Some(load) => Some(
                fs_duct::modal::plane_mouth_kernel(
                    &self.duct,
                    gas,
                    room,
                    LossModel::Lossless,
                    load.n_modes,
                    load.extra_slices,
                    fs_duct::modal::MouthKernelSpec {
                        step_s: clock.gas_step_s(),
                        bin_hz: load.bin_hz,
                        band_hz: load.band_hz,
                        tail_tolerance: load.tail_tolerance,
                    },
                )
                .map_err(|e| format!("mouth_load refused: {e}"))?,
            ),
        };
        let mouth_info = match &mouth {
            None => None,
            Some(kernel) => {
                let minimum = minimum_resistance_ratio(kernel);
                if !(minimum >= MOUTH_LOAD_MINIMUM_RESISTANCE_RATIO) {
                    return Err(format!(
                        "mouth_load refused: the realized load's smallest Re Z / Z_inf is {minimum:.3e}, below {MOUTH_LOAD_MINIMUM_RESISTANCE_RATIO:e}: a negative resistance is gain at the mouth"
                    ));
                }
                self.mouth_load.map(|record| MouthLoadInfo {
                    record,
                    taps: kernel.weights.len(),
                    mouth_radius_m: kernel.mouth_radius_m,
                    z_inf: kernel.z_inf,
                    dropped_tail_energy: kernel.dropped_tail_energy,
                    dropped_sum: kernel.dropped_sum,
                    bin_hz: kernel.bin_hz,
                    evaluated_hz: kernel.evaluated_hz,
                    minimum_resistance_ratio: minimum,
                })
            }
        };
        let lip = NonlinearLip::new_with_mouth(
            &self.duct,
            source,
            gas.temperature,
            gas.pressure,
            wall_gas,
            self.lip.clone(),
            self.cup_volume_m3,
            clock,
            mouth.as_ref(),
        )?;
        let lip = match self.mouth_reservoir {
            Some(record) => lip.with_finite_mouth(FiniteMouthSpec {
                volume_m3: record.volume_m3,
                feed_area_m2: record.feed_area_m2,
                initial_pressure_pa: record.initial_pressure_pa,
                initial_temperature_k: record.initial_temperature_k,
            }).map_err(|e| format!("mouth_reservoir refused: {e}"))?,
            None => lip,
        };
        let lip = match self.jet_noise {
            Some(card) => lip.with_jet_noise(card).map_err(|e| format!("jet_noise refused: {e}"))?,
            None => lip,
        };
        let lip = match self.lip_mode2 {
            Some(mode) => lip.with_lip_mode2(mode).map_err(|e| format!("lip_response_mode2 refused: {e}"))?,
            None => lip,
        };
        let lip = match self.lip_area_exponent {
            Some(q) => lip.with_lip_area_exponent(q).map_err(|e| format!("lip_area_exponent refused: {e}"))?,
            None => lip,
        };
        let lip = match self.constriction_noise {
            Some(card) => lip
                .with_constriction_noise(card)
                .map_err(|e| format!("constriction_noise refused: {e}"))?,
            None => lip,
        };
        let lip = match &self.bore_wall {
            Some((record, _)) => lip
                .with_bore_wall(&self.duct, *record)
                .map_err(|e| format!("bore_wall refused: {e}"))?,
            None => lip,
        };
        let lip = match &self.bell_shell {
            Some((record, _)) => lip
                .with_bell_shell(&self.duct, *record)
                .map_err(|e| format!("bore_wall shell-brass refused: {e}"))?,
            None => lip,
        };
        let bell_shell_info = lip.bell_shell_info();
        let bore_wall_info = match (&self.bore_wall, lip.bore().shell_receipt()) {
            (Some((record, _)), Some(receipt)) => Some(BoreWallInfo {
                record: *record,
                cells: receipt.cells,
                membrane_cells: receipt.membrane_cells,
                thick_cells: receipt.thick_cells,
                largest_half_angle_rad: receipt.largest_half_angle_rad,
                largest_thick_half_angle_rad: receipt.largest_thick_half_angle_rad,
                lowest_resonance_hz: receipt.lowest_resonance_hz,
                highest_resonance_hz: receipt.highest_resonance_hz,
            }),
            _ => None,
        };
        let (receiver, observer_gas) =
            self.scalar_receiver(mouth_radius_m, clock.mechanics_rate_hz())?;
        let wall = lip
            .wall_receipt()
            .ok_or_else(|| "the nonlinear bore bound no wall memory".to_string())?;
        let bore = LipBoreInfo::Nonlinear(NonlinearBoreInfo {
            record,
            mechanics_rate_hz: clock.mechanics_rate_hz(),
            gas_step_s: clock.gas_step_s(),
            cells: lip.bore().centers_m().len(),
            radiation_taps: lip.bore().radiation_taps(),
            mouth_load: mouth_info,
            wall_distinct_radii: wall.distinct_radii,
            wall_maximum_arms: wall.maximum_arms,
            wall_maximum_complex_error: wall.maximum_complex_error,
            wall_maximum_loss_error: wall.maximum_loss_error,
            bore_wall: bore_wall_info,
            bell_shell: bell_shell_info,
        });
        Ok(self.performance(
            mouth_radius_m,
            bore,
            Engine::Nonlinear(Box::new(lip)),
            receiver,
            observer_gas,
            record.mechanics_per_output as usize,
        ))
    }

    // Both bore engines deliver actual bell volume flow to the same aperture
    // owner; only their mechanics/output clock ratio differs.
    fn scalar_receiver(
        &self,
        mouth_radius_m: f64,
        mechanics_rate_hz: u32,
    ) -> Result<(BaffledPressure, Option<GasPropagationReport>), String> {
        let gas = self.ambient_gas();
        if let Some(propagation) = &self.observer_gas {
            let opening = CircularOpening { center_m: [0.0; 2], radius_m: mouth_radius_m };
            let (receiver, report) = BaffledPressure::circular_openings_multirate_gas_differentiating(
                &[opening], mechanics_rate_hz, self.rate_hz, self.microphone,
                gas, propagation, &mut || false,
            ).map_err(|error| format!("observer_gas refused: {error}"))?;
            Ok((receiver, Some(report)))
        } else {
            let receiver = BaffledPressure::circular_outlet_multirate_differentiating(
                mouth_radius_m, mechanics_rate_hz, self.rate_hz, self.microphone,
                RayleighMedium::from(gas),
            )?;
            Ok((receiver, None))
        }
    }

    fn performance(
        &self,
        mouth_radius_m: f64,
        bore: LipBoreInfo,
        engine: Engine,
        receiver: BaffledPressure,
        observer_gas: Option<GasPropagationReport>,
        ratio: usize,
    ) -> LipPerformance {
        let gas = self.bore_gas();
        let info = LipInfo {
            input_hash: self.input_hash,
            sample_rate_hz: self.rate_hz,
            samples: self.samples,
            full_scale_pa: self.full_scale_pa,
            segments: self.duct.segments.len(),
            bore_length_m: self.duct.segments.iter().map(Segment::axial_length).sum(),
            mouth_radius_m,
            lip: self.lip.clone(),
            lip_mode2: self.lip_mode2,
            lip_area_exponent: self.lip_area_exponent,
            rim_stop: self
                .lip_rim_stop
                .zip(self.rim_max_gap_m().ok())
                .map(|(kind, max_gap_m)| RimStopInfo { kind, max_gap_m }),
            bore_wall: self.bore_wall.clone(),
            bell_shell: self.bell_shell.clone(),
            lip_source: self.lip_source.clone(),
            bore_source: self.bore_source.clone(),
            cup_volume_m3: self.cup_volume_m3,
            gas_temperature_k: gas.temperature,
            sound_speed_m_s: gas.sound_speed,
            microphone: self.microphone,
            observer_gas,
            bore,
            control_rate_hz: self.schedule.control_rate_hz,
            has_lip_tension_track: self.tension_track.is_some(),
        };
        let inputs = match &engine {
            Engine::Linear { mouth, .. } => mouth.len(),
            #[cfg(feature = "nonlinear-brass")]
            Engine::Nonlinear(_) => 1,
        };
        LipPerformance {
            info,
            engine,
            receiver,
            schedule: self.schedule.clone(),
            blow_track: self.blow_track.clone(),
            tension_track: self.tension_track.clone(),
            completed: 0,
            applied_pressure: f64::NAN,
            applied_tension: f64::NAN,
            poisoned: false,
            mouthpiece: vec![0.0; MAX_BLOCK],
            outgoing: vec![0.0; MAX_BLOCK],
            trace: vec![0.0; MAX_BLOCK * ratio * inputs],
            ratio,
            window_len: u64::from(self.rate_hz / 10).max(1),
            acc: Vec::new(),
            tap: None,
        }
    }
}

/// The temperatures the nonlinear field may reach: for every species, the
/// NASA-9 region that holds the room. The packs state no caloric continuity
/// across their region boundaries, so a state that would cross one refuses
/// instead. A bound shared with the next region belongs to that region.
#[cfg(feature = "nonlinear-brass")]
fn room_region(
    source: &fs_thermochem::mixture::FrozenIdealGasMixtureModelV1,
    room_k: f64,
) -> Result<[f64; 2], String> {
    let mut domain = [0.0_f64, f64::INFINITY];
    for component in source.components() {
        let regions = component.regions();
        let shared = |index: usize, maximum: f64| {
            regions
                .get(index + 1)
                .is_some_and(|next| next.temperature_min().value().to_bits() == maximum.to_bits())
        };
        let (index, region) = regions
            .iter()
            .enumerate()
            .find(|(index, region)| {
                let (minimum, maximum) = (
                    region.temperature_min().value(),
                    region.temperature_max().value(),
                );
                room_k >= minimum
                    && (room_k < maximum || (!shared(*index, maximum) && room_k <= maximum))
            })
            .ok_or_else(|| {
                format!(
                    "bore nonlinear: the room ({room_k} K) lies outside every NASA-9 region of {}",
                    component.species()
                )
            })?;
        let maximum = region.temperature_max().value();
        domain[0] = domain[0].max(region.temperature_min().value());
        domain[1] = domain[1].min(if shared(index, maximum) {
            maximum.next_down()
        } else {
            maximum
        });
    }
    Ok(domain)
}

/// The opt-in rim stop as realized: its limit case and the largest gap it allows (twice the entry radius of the first bore cone).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RimStopInfo {
    /// The limit case of the unmeasured restitution.
    pub kind: super::RimStopKind,
    /// The largest lip gap [m].
    pub max_gap_m: f64,
}

/// What a realized performance reports about itself.
#[derive(Debug, Clone)]
pub struct LipInfo {
    /// Domain-separated hash of the exact input bytes.
    pub input_hash: ContentHash,
    /// Mechanics and output clock [Hz].
    pub sample_rate_hz: u32,
    /// Declared samples.
    pub samples: u64,
    /// Pascals at PCM full scale.
    pub full_scale_pa: f64,
    /// Bore segments after Bessel expansion.
    pub segments: usize,
    /// Axial length of the bore [m].
    pub bore_length_m: f64,
    /// Radius of the radiating outlet [m].
    pub mouth_radius_m: f64,
    /// The lip as the loop uses it, with its provenance sentence.
    pub lip: LipIsland,
    /// The opt-in second lip mode, when the file declared one.
    pub lip_mode2: Option<LipMode2>,
    /// The opt-in flow-area exponent, when the file declared one.
    pub lip_area_exponent: Option<f64>,
    /// The opt-in rim stop, when the file declared one.
    pub rim_stop: Option<RimStopInfo>,
    /// The opt-in breathing bore wall with its provenance sentence, when the file declared `bore_wall thin-brass`.
    pub bore_wall: Option<(BoreWallRecord, String)>,
    /// The opt-in modal bell wall with its provenance sentence, when the file declared `bore_wall shell-brass`.
    pub bell_shell: Option<(BellShellRecord, String)>,
    /// Provenance sentence of the lip numbers.
    pub lip_source: String,
    /// Provenance sentence of the bore geometry.
    pub bore_source: String,
    /// Lumped mouthpiece-cup volume [m^3].
    pub cup_volume_m3: f64,
    /// Bore gas temperature [K] (the `bore_gas` row's gas when the file declared one, else the ambient gas).
    pub gas_temperature_k: f64,
    /// Bore sound speed [m/s] (as above).
    pub sound_speed_m_s: f64,
    /// The declared microphone.
    pub microphone: CircularOutletReceiver,
    /// Original per-path gas propagation, absent for the legacy lossless receiver.
    /// Its extra causal padding is included once in the receiver's filter delay.
    pub observer_gas: Option<GasPropagationReport>,
    /// How the bore was realized, with its receipts.
    pub bore: LipBoreInfo,
    /// Control clock of the pressure phrase [Hz].
    pub control_rate_hz: u32,
    /// Whether a pre-stress (lip tension) track drives the lip.
    pub has_lip_tension_track: bool,
}

/// The realized bore of a performance.
#[derive(Debug, Clone)]
pub enum LipBoreInfo {
    /// The multimodal reflectance and its total-mouth-flow observation.
    Linear(LinearBoreInfo),
    /// The shared nonlinear gas bore.
    Nonlinear(NonlinearBoreInfo),
}

/// Receipts of the linear lines.
#[derive(Debug, Clone)]
pub struct LinearBoreInfo {
    /// Realization receipt of the multimodal reflectance line.
    pub line: MmComboReceipt,
    /// Observed signed mouth coordinates; one for the legacy total-flow route.
    pub mouth_flow_modes: usize,
    /// Present only for the explicitly selected full modal baffled exterior.
    pub modal_radiation: Option<ModalOutletReport>,
    /// Taps of the mouth-flow line.
    pub mouth_flow_taps: usize,
    /// Ladder rungs the mouth-flow line traversed.
    pub mouth_flow_rungs: usize,
    /// Relative tail energy the mouth-flow ladder was asked to reach (the
    /// owner's default unless the file declares `mouth_flow_tail_tol`).
    pub mouth_flow_tail_tol: f64,
    /// Every mouth-flow rung: `(taps, worst modal tail energy, worst modal tap change)`.
    pub mouth_flow_ladder: Vec<(usize, f64, f64)>,
    /// The passivity gain the player tract's airway line took ([`PassivityGain`]); `None` without a tract or when its
    /// line was admitted without one.
    pub player_tract_passivity_gain: Option<PassivityGain>,
}

/// Receipts of the nonlinear gas bore.
#[derive(Debug, Clone, Copy)]
pub struct NonlinearBoreInfo {
    /// The declared record.
    pub record: NonlinearBoreRecord,
    /// Mechanics (microphone) clock [Hz].
    pub mechanics_rate_hz: u32,
    /// Fixed gas step [s].
    pub gas_step_s: f64,
    /// Gas cells.
    pub cells: usize,
    /// Past bell flows the radiating end convolves (the piston's, or the `mouth_load modal` kernel's).
    pub radiation_taps: usize,
    /// The `mouth_load modal` kernel's receipts, or `None` for the baffled piston.
    pub mouth_load: Option<MouthLoadInfo>,
    /// Distinct wall radii bound.
    pub wall_distinct_radii: usize,
    /// Largest viscous or thermal arm count in any cell.
    pub wall_maximum_arms: usize,
    /// Worst sampled complex error of the wall operators on the band.
    pub wall_maximum_complex_error: f64,
    /// Worst sampled error of their dissipative parts.
    pub wall_maximum_loss_error: f64,
    /// The breathing wall's admission receipt, or `None` for the rigid-walled bore (the row absent or `rigid`).
    pub bore_wall: Option<BoreWallInfo>,
    /// The modal bell wall's admission receipt, or `None` without a `bore_wall shell-brass` row.
    pub bell_shell: Option<BellShellInfo>,
}

/// What binding `bore_wall shell-brass` admitted.
#[derive(Debug, Clone, Copy)]
pub struct BellShellInfo {
    /// The declared thickness, loss factor and cut-off.
    pub record: BellShellRecord,
    /// Bore cells the shell covers (the bell's).
    pub cells_covered: usize,
    /// Modes that couple to the bore (axisymmetric, below the cut-off).
    pub coupling_modes: usize,
    /// Modes the cut-off admitted, of every circumferential harmonic.
    pub modes_below_cutoff: usize,
    /// The first bore cell with a wall band: the clamp is at the bore face at or after the throat that starts it.
    pub first_covered_cell: usize,
    /// Lowest and highest coupling mode frequency [Hz].
    pub lowest_mode_hz: f64,
    /// Highest coupling mode frequency [Hz].
    pub highest_mode_hz: f64,
}

/// What binding `bore_wall thin-brass` admitted: the declared numbers and the wall's cells.
#[derive(Debug, Clone, Copy)]
pub struct BoreWallInfo {
    /// The declared thickness and quality factor.
    pub record: BoreWallRecord,
    /// Cells of the bore.
    pub cells: usize,
    /// Cells whose wall ran the thin-membrane expression (the classical cone factor applies).
    pub membrane_cells: usize,
    /// Cells whose wall ran Lame's thick cylinder (the cylinder stiffness is carried unchanged).
    pub thick_cells: usize,
    /// Largest cell half-angle [rad].
    pub largest_half_angle_rad: f64,
    /// Largest half-angle among the thick-wall cells [rad]; 0 with none.
    pub largest_thick_half_angle_rad: f64,
    /// Lowest breathing resonance of any cell [Hz].
    pub lowest_resonance_hz: f64,
    /// Highest breathing resonance of any cell [Hz].
    pub highest_resonance_hz: f64,
}

/// The passive bore's response at one frequency, at its input plane (the mouthpiece cup rim).
#[derive(Debug, Clone, Copy)]
pub struct BoreResponse {
    /// Frequency [Hz].
    pub frequency_hz: f64,
    /// Multimodal plane-mode input impedance `[re, im]` [Pa s/m^3]: what the playing loop's reflectance is built from.
    pub impedance_multimodal: [f64; 2],
    /// Plane-wave (scalar transfer-matrix) input impedance `[re, im]` [Pa s/m^3], retained as a comparison to the modal bore.
    pub impedance_scalar: [f64; 2],
    /// Mouth volume flow over input volume flow `[re, im]` (scalar chain).
    pub mouth_over_input_flow: [f64; 2],
}

/// A scan of the passive bore with the quantities needed to read it.
#[derive(Debug, Clone)]
pub struct BoreScan {
    /// Plane-wave characteristic impedance of the input plane `rho c / (pi r^2)` [Pa s/m^3].
    pub input_characteristic_impedance: f64,
    /// Gas density [kg/m^3].
    pub density_kg_m3: f64,
    /// Distance from the outlet plane to the microphone [m].
    pub microphone_distance_m: f64,
    /// One row per requested frequency.
    pub rows: Vec<BoreResponse>,
}

/// Lip and bore statistics over one tenth of a second of a render: whether the lip closes and how large
/// the pressure at the mouthpiece is. A diagnostic of the loop, not a graded quantity.
#[derive(Debug, Clone, Copy)]
pub struct LipWindow {
    /// Start of the window [s].
    pub start_s: f64,
    /// Mean pressure at the mouthpiece (the mouth's DC operating point) [Pa].
    pub bore_mean_pa: f64,
    /// RMS of the alternating mouthpiece pressure [Pa].
    pub bore_ac_rms_pa: f64,
    /// Half the peak-to-peak mouthpiece pressure [Pa].
    pub bore_ac_peak_pa: f64,
    /// Smallest lip gap in the window [m] (0 = the lip closed).
    pub gap_min_m: f64,
    /// Mean lip gap [m] (mean of the per-chunk means).
    pub gap_mean_m: f64,
    /// Largest lip gap [m].
    pub gap_max_m: f64,
    /// Fraction of the control-tick chunks in which the lip closed at least once.
    pub closed_chunk_fraction: f64,
    /// Signed tension-actuator work [J] by lip mode. Each instantaneous
    /// control change belongs to the window containing its control sample;
    /// it is neither blowing-pressure work nor dissipated heat.
    pub lip_tension_work_j: [f64; 2],
}

#[derive(Debug, Clone, Copy)]
struct WindowAcc {
    n: u64,
    sum: f64,
    sum_sq: f64,
    p_min: f64,
    p_max: f64,
    chunks: u64,
    closed_chunks: u64,
    gap_min: f64,
    gap_max: f64,
    gap_mean_sum: f64,
    tension_work: super::TensionWork,
}

impl WindowAcc {
    fn new() -> Self {
        Self {
            n: 0,
            sum: 0.0,
            sum_sq: 0.0,
            p_min: f64::INFINITY,
            p_max: f64::NEG_INFINITY,
            chunks: 0,
            closed_chunks: 0,
            gap_min: f64::INFINITY,
            gap_max: 0.0,
            gap_mean_sum: 0.0,
            tension_work: super::TensionWork::default(),
        }
    }
}

enum Engine {
    Linear {
        voice: BrassVoice,
        mouth: Vec<DelayedFilter>,
    },
    #[cfg(feature = "nonlinear-brass")]
    Nonlinear(Box<super::nonlinear::NonlinearLip>),
}

impl Engine {
    fn apply(&mut self, control: BrassControl) -> Result<(), String> {
        match self {
            Self::Linear { voice, .. } => voice.apply(control).map_err(|e| e.to_string()),
            #[cfg(feature = "nonlinear-brass")]
            Self::Nonlinear(lip) => lip.apply(control),
        }
    }

    fn diagnostics(&self) -> &[BrassBlockDiag] {
        match self {
            Self::Linear { voice, .. } => voice.diagnostics(),
            #[cfg(feature = "nonlinear-brass")]
            Self::Nonlinear(lip) => lip.diagnostics(),
        }
    }

    fn lip_tension_work_j(&self) -> [f64; 2] {
        match self {
            Self::Linear { voice, .. } => voice.lip_tension_work_j(),
            #[cfg(feature = "nonlinear-brass")]
            Self::Nonlinear(lip) => lip.lip_tension_work_j(),
        }
    }

    fn first_order_fallbacks(&self) -> Option<u64> {
        match self {
            Self::Linear { .. } => None,
            #[cfg(feature = "nonlinear-brass")]
            Self::Nonlinear(lip) => Some(lip.bore().first_order_fallbacks()),
        }
    }

    /// The breathing walls' energy account on the nonlinear bore, if one is bound.
    #[cfg(feature = "nonlinear-brass")]
    fn wall_ledger(&self) -> Option<fs_duct::nonlinear::bore::ShellLedger> {
        match self {
            Self::Linear { .. } => None,
            Self::Nonlinear(lip) => lip.bore().shell_ledger(),
        }
    }

    fn rim_stop_hits(&self) -> Option<u64> {
        match self {
            Self::Linear { voice, .. } => voice.rim_stop_hits(),
            #[cfg(feature = "nonlinear-brass")]
            Self::Nonlinear(_) => None,
        }
    }

    fn jet_noise_coverage(&self) -> Option<crate::turbulence_noise::JetNoiseCoverage> {
        match self {
            Self::Linear { voice, .. } => voice.jet_noise_coverage(),
            #[cfg(feature = "nonlinear-brass")]
            Self::Nonlinear(lip) => lip.jet_noise_coverage(),
        }
    }
}

/// A realized performance: the lip, the bore and the microphone, advanced
/// together on one physical clock.
pub struct LipPerformance {
    info: LipInfo,
    engine: Engine,
    receiver: BaffledPressure,
    schedule: GestureSchedule,
    blow_track: String,
    tension_track: Option<String>,
    completed: u64,
    applied_pressure: f64,
    applied_tension: f64,
    poisoned: bool,
    mouthpiece: Vec<f64>,
    outgoing: Vec<f64>,
    trace: Vec<f64>,
    ratio: usize,
    window_len: u64,
    acc: Vec<WindowAcc>,
    tap: Option<Vec<f64>>,
}

fn refusal(what: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(what))
}

impl LipPerformance {
    /// The realized description.
    #[must_use]
    pub const fn info(&self) -> &LipInfo {
        &self.info
    }

    /// Existing microphone filter latency, in output samples. Physical
    /// source-to-listener travel remains in the pressure and is not included.
    #[must_use]
    pub fn filter_delay_samples(&self) -> f64 {
        self.receiver.filter_delay_samples()
    }

    /// The loop's per-chunk diagnostics (lock estimate, gap and power ledger).
    #[must_use]
    pub fn diagnostics(&self) -> &[BrassBlockDiag] {
        self.engine.diagnostics()
    }

    /// Nonlinear bore only (`None` for the linear bore): the gas steps accepted after a
    /// first-order retry of the failed second-order update, so far.
    #[must_use]
    pub fn first_order_fallbacks(&self) -> Option<u64> {
        self.engine.first_order_fallbacks()
    }

    /// The breathing walls' energy account so far (gas work, work against the room, dissipated and stored energy): `None` for the
    /// linear bore and for a nonlinear bore without a `bore_wall` row.
    #[cfg(feature = "nonlinear-brass")]
    #[must_use]
    pub fn wall_ledger(&self) -> Option<fs_duct::nonlinear::bore::ShellLedger> {
        self.engine.wall_ledger()
    }

    /// Samples at which the opt-in rim stop acted so far; `None` without the `lip_rim_stop` row (and for `bore nonlinear`).
    #[must_use]
    pub fn rim_stop_hits(&self) -> Option<u64> {
        self.engine.rim_stop_hits()
    }

    /// Actual accepted finite upstream pressure, temperature and gas inventory.
    /// Absent for a linear bore or the nonlinear held-mouth default.
    #[cfg(feature = "nonlinear-brass")]
    #[must_use]
    pub fn finite_mouth_state(&self) -> Option<super::nonlinear::FiniteMouthState> {
        match &self.engine {
            Engine::Linear { .. } => None,
            Engine::Nonlinear(lip) => lip.finite_mouth_state(),
        }
    }

    /// Last accepted finite-mouth exchange on the native gas clock.
    #[cfg(feature = "nonlinear-brass")]
    #[must_use]
    pub fn finite_mouth_report(&self) -> Option<&super::nonlinear::FiniteMouthReport> {
        match &self.engine {
            Engine::Linear { .. } => None,
            Engine::Nonlinear(lip) => lip.finite_mouth_report(),
        }
    }

    /// The opt-in lip jet noise's open-gate mechanics steps so far by how much of the band each could represent at the
    /// lip's clock; `None` without the `jet_noise` row.
    #[must_use]
    pub fn jet_noise_coverage(&self) -> Option<crate::turbulence_noise::JetNoiseCoverage> {
        self.engine.jet_noise_coverage()
    }

    /// Start recording the pressure at the mouthpiece (the internal pressure the
    /// lip sees, DC included) for every sample rendered from now on.
    pub fn record_mouthpiece(&mut self) {
        self.tap.get_or_insert_with(Vec::new);
    }

    /// The recorded mouthpiece pressure [Pa] (empty unless [`Self::record_mouthpiece`] ran).
    #[must_use]
    pub fn mouthpiece_pressure(&self) -> &[f64] {
        self.tap.as_deref().unwrap_or(&[])
    }

    /// Cumulative signed work [J] supplied by the tension actuator to each
    /// lip mode. Uses accepted controls in sample order, independent of audio
    /// block partitioning. It does not assert a closed lip/bore energy balance.
    #[must_use]
    pub fn lip_tension_work_j(&self) -> [f64; 2] {
        self.engine.lip_tension_work_j()
    }

    /// Lip and bore statistics per tenth of a second rendered so far.
    #[must_use]
    pub fn windows(&self) -> Vec<LipWindow> {
        let rate = f64::from(self.info.sample_rate_hz);
        self.acc
            .iter()
            .enumerate()
            .filter(|(_, a)| a.n > 0 && a.chunks > 0)
            .map(|(k, a)| {
                let n = a.n as f64;
                let mean = a.sum / n;
                LipWindow {
                    start_s: (k as u64 * self.window_len) as f64 / rate,
                    bore_mean_pa: mean,
                    bore_ac_rms_pa: (a.sum_sq / n - mean * mean).max(0.0).sqrt(),
                    bore_ac_peak_pa: 0.5 * (a.p_max - a.p_min),
                    gap_min_m: a.gap_min,
                    gap_mean_m: a.gap_mean_sum / a.chunks as f64,
                    gap_max_m: a.gap_max,
                    closed_chunk_fraction: a.closed_chunks as f64 / a.chunks as f64,
                    lip_tension_work_j: a.tension_work.total_j(),
                }
            })
            .collect()
    }

    fn apply_controls(&mut self, tick: u64) -> Result<(), String> {
        let pressure = self
            .schedule
            .sample(&self.blow_track, tick)
            .map_err(|e| format!("pressure schedule refused: {e:?}"))?;
        if pressure != self.applied_pressure {
            self.engine
                .apply(BrassControl::SetBlowingPressure(pressure))
                .map_err(|e| format!("blowing pressure refused: {e}"))?;
            self.applied_pressure = pressure;
        }
        if let Some(track) = &self.tension_track {
            let tension = self
                .schedule
                .sample(track, tick)
                .map_err(|e| format!("lip tension schedule refused: {e:?}"))?;
            if tension != self.applied_tension {
                self.engine
                    .apply(BrassControl::SetLipTension(tension))
                    .map_err(|e| format!("lip tension refused: {e}"))?;
                self.applied_tension = tension;
            }
        }
        Ok(())
    }

    fn advance(&mut self, output: &mut [f64]) -> Result<(), String> {
        let rate = u64::from(self.info.sample_rate_hz);
        let control = u64::from(self.info.control_rate_hz);
        let mut done = 0usize;
        while done < output.len() {
            let n = self.completed + done as u64;
            let tick = tick_at(n, control, rate);
            self.apply_controls(tick)?;
            let until_next = tick_start(tick + 1, control, rate) - n;
            let chunk = usize::try_from(until_next)
                .unwrap_or(usize::MAX)
                .min(output.len() - done)
                .min(MAX_BLOCK)
                .max(1);
            match &mut self.engine {
                Engine::Linear { voice, .. } => voice
                    .step_block_outgoing(&mut self.mouthpiece[..chunk], &mut self.outgoing[..chunk])
                    .map_err(|e| format!("brass loop refused: {e}"))?,
                #[cfg(feature = "nonlinear-brass")]
                Engine::Nonlinear(lip) => lip.step_block(
                    &mut self.mouthpiece[..chunk],
                    &mut self.trace[..chunk * self.ratio],
                )?,
            }
            let diag = self.engine.diagnostics().last().cloned();
            let window = usize::try_from(n / self.window_len).unwrap_or(usize::MAX);
            while self.acc.len() <= window {
                self.acc.push(WindowAcc::new());
            }
            if let Some(tap) = self.tap.as_mut() {
                tap.extend_from_slice(&self.mouthpiece[..chunk]);
            }
            let acc = &mut self.acc[window];
            for &p in &self.mouthpiece[..chunk] {
                acc.n += 1;
                acc.sum += p;
                acc.sum_sq += p * p;
                acc.p_min = acc.p_min.min(p);
                acc.p_max = acc.p_max.max(p);
            }
            if let Some(d) = diag {
                acc.tension_work = acc.tension_work.with_increment(d.lip_tension_work_j)?;
                acc.chunks += 1;
                acc.closed_chunks += u64::from(d.gap_min_m <= 0.0);
                acc.gap_min = acc.gap_min.min(d.gap_min_m);
                acc.gap_max = acc.gap_max.max(d.gap_max_m);
                acc.gap_mean_sum += d.gap_mean_m;
            }
            for i in 0..chunk {
                let frame = match &mut self.engine {
                    Engine::Linear { mouth, .. } => {
                        for (flow, filter) in self.trace.iter_mut().zip(mouth.iter_mut()) {
                            *flow = filter.push(self.outgoing[i])
                                .map_err(|e| format!("mouth-flow line refused: {e:?}"))?;
                        }
                        &self.trace[..mouth.len()]
                    }
                    #[cfg(feature = "nonlinear-brass")]
                    Engine::Nonlinear(_) => &self.trace[i * self.ratio..(i + 1) * self.ratio],
                };
                output[done + i] = self.receiver.step_trace(frame)?;
            }
            done += chunk;
        }
        self.completed += output.len() as u64;
        Ok(())
    }
}

impl PressureRenderer for LipPerformance {
    fn samples_rendered(&self) -> u64 {
        self.completed
    }

    fn max_block_len(&self) -> usize {
        MAX_BLOCK
    }

    fn validate_sample_rate(&self, sample_rate_hz: u32) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if sample_rate_hz != self.info.sample_rate_hz {
            return Err(RenderError::Sizing {
                what: "lip performance output clock differs from the requested clock",
            });
        }
        Ok(())
    }

    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if samples > self.info.samples - self.completed {
            return Err(RenderError::Sizing {
                what: "request runs past the samples the lip performance declares",
            });
        }
        Ok(())
    }

    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        if output.len() > MAX_BLOCK {
            return Err(RenderError::Sizing {
                what: "lip callback exceeds its admitted block length",
            });
        }
        self.validate_sample_count(output.len() as u64)?;
        match self.advance(output) {
            Ok(()) => Ok(()),
            Err(what) => {
                // A failed physical callback leaves a prefix advanced.
                self.poisoned = true;
                Err(refusal(what))
            }
        }
    }
}

#[cfg(test)]
#[path = "performance_bank_limit.rs"]
mod bank_limit;

#[cfg(test)]
#[path = "performance_observer_gas.rs"]
mod observer_gas_tests;

#[cfg(test)]
#[path = "../render/schedule/reed/field_input/tests/support.rs"]
mod gas_test_support;

#[cfg(test)]
mod tests {
    use super::*;

    // The canonical pressure phrase of examples/reed-duct.performance.
    const SCHEDULE: &str = "frankensim-gesture-schedule-v1\n\
        control_rate_hz\t700\n\
        tracks\t1\n\
        track\tblow\tblowing-pressure\n\
        initial\tpressure-pa\t0e0\n\
        events\t5\n\
        event\t0e0\t1e-2\tpressure-pa\t2.8e3\n\
        event\t5e-3\t1e-2\tpressure-pa\t3.5e3\n\
        event\t2e-2\t0e0\tpressure-pa\t0e0\n\
        event\t4e-2\t5e-3\tpressure-pa\t2.8e3\n\
        event\t6e-2\t0e0\tpressure-pa\t0e0\n";

    // Test fixture, not data: a short brass-like bore and a lip-shaped
    // response so every owner in the chain is exercised.
    fn fixture(extra: &str) -> String {
        format!(
            "frankensim-lip-performance-v1\n\
             # test fixture, not data\n\
             audio 24000 2400\n\
             full_scale_pa 30000\n\
             ambient moist_air 293.15 101325 0\n\
             lip_response 200 7 0.1 0.012 0.0005 0.0001\n\
             lip_source test fixture, not a measurement\n\
             bore_source test fixture, not a measurement\n\
             cup_volume_m3 1.5e-6\n\
             termination flanged\n\
             modes 3 1\n\
             cylinder 0.006 0.30\n\
             cone 0.006 0.012 0.25\n\
             cone 0.012 0.035 0.12\n\
             microphone 0 0 1 32 128 3000\n\
             {extra}schedule\n{SCHEDULE}"
        )
    }

    #[test]
    fn jet_noise_parses_and_is_refused_by_name() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        let ok = parse("jet_noise 7 1.3e-3 0.045 0.115 1700\n").expect("parses");
        assert_eq!(
            ok.jet_noise(),
            Some(&JetNoiseCard { seed: 7, gain: 1.3e-3, strouhal_lo: 0.045, strouhal_hi: 0.115, reynolds_onset: 1700.0 })
        );
        assert_eq!(parse("").expect("parses").jet_noise(), None);
        for (extra, needle) in [
            ("jet_noise\n", "missing field"),
            ("jet_noise 7 1.3e-3 0.045 0.115\n", "missing field"),
            ("jet_noise 7 1.3e-3 0.045 0.115 1700 9\n", "unexpected extra field"),
            ("jet_noise x 1.3e-3 0.045 0.115 1700\n", "malformed field"),
            ("jet_noise 7 nan 0.045 0.115 1700\n", "must be finite"),
            ("jet_noise 7 0 0.045 0.115 1700\n", "gain must be finite and positive"),
            ("jet_noise 7 1.3e-3 0.115 0.045 1700\n", "Strouhal band"),
            ("jet_noise 7 1.3e-3 0.045 0.115 -1\n", "reynolds onset"),
            ("jet_noise 7 1.3e-3 0.045 0.115 1700\njet_noise 8 1.3e-3 0.045 0.115 1700\n", "may only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    #[test]
    fn lip_response_mode2_parses_and_is_refused_by_name() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        let ok = parse("lip_response_mode2 244.5 11.0 -0.073\n").expect("parses");
        let mode = ok.lip_mode2().expect("declared");
        let omega = core::f64::consts::TAU * 244.5;
        assert!((mode.mass_kg - 0.0001 / 0.073).abs() <= 1e-15 * mode.mass_kg);
        assert!((mode.stiffness_n_m - mode.mass_kg * omega * omega).abs() <= 1e-12 * mode.stiffness_n_m);
        assert!((mode.damping_n_s_m - mode.mass_kg * omega / 11.0).abs() <= 1e-12 * mode.damping_n_s_m);
        assert_eq!(mode.sign, -1.0);
        assert!(parse("").expect("parses").lip_mode2().is_none());
        for (extra, needle) in [
            ("lip_response_mode2\n", "missing field"),
            ("lip_response_mode2 244.5 11.0\n", "missing field"),
            ("lip_response_mode2 244.5 11.0 -0.073 1\n", "unexpected extra field"),
            ("lip_response_mode2 x 11.0 -0.073\n", "malformed field"),
            ("lip_response_mode2 244.5 0 -0.073\n", "expected a positive number"),
            ("lip_response_mode2 244.5 11.0 0\n", "nonzero"),
            ("lip_response_mode2 244.5 11.0 nan\n", "must be finite"),
            ("lip_response_mode2 244.5 11.0 -0.073\nlip_response_mode2 250 11.0 -0.07\n", "may only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        // The second mode shares the first's face area, so it needs lip_response and not the lumped lip record.
        let lumped = fixture("lip_response_mode2 244.5 11.0 -0.073\n").replace(
            "lip_response 200 7 0.1 0.012 0.0005 0.0001\n",
            "lip 1.0e-4 160 0.05 0.012 0.0005 0.0001\n",
        );
        let error = LipDescription::from_bytes(lumped.as_bytes()).expect_err("lumped lip");
        assert!(error.contains("needs the lip_response record"), "{error}");
    }

    #[test]
    fn lip_response_mode2_reaches_the_radiated_pressure_and_a_file_without_it_is_unchanged() {
        let build = |extra: &str| pressure(&LipDescription::from_bytes(fixture(extra).as_bytes()).expect("parses"));
        let plain = build("");
        let again = build("");
        assert!(
            plain.iter().zip(&again).all(|(a, b)| a.to_bits() == b.to_bits()),
            "the same file must render the same bits"
        );
        let two = build("lip_response_mode2 244.5 11.0 -0.073\n");
        assert_ne!(plain, two, "a second lip mode must reach the radiated pressure");
        assert!(two.iter().all(|x| x.is_finite()));
        // The same bytes with the record are a different input (the hash covers the file); the first mode is untouched.
        let (a, b) = (
            LipDescription::from_bytes(fixture("").as_bytes()).expect("parses"),
            LipDescription::from_bytes(fixture("lip_response_mode2 244.5 11.0 -0.073\n").as_bytes()).expect("parses"),
        );
        assert_ne!(a.input_hash, b.input_hash);
        assert_eq!(a.lip.mass_kg.to_bits(), b.lip.mass_kg.to_bits());
        assert_eq!(a.lip.stiffness_n_m.to_bits(), b.lip.stiffness_n_m.to_bits());
    }

    #[test]
    fn lip_area_exponent_parses_is_refused_by_name_and_reaches_the_radiated_pressure() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        assert_eq!(parse("lip_area_exponent 2\n").expect("parses").lip_area_exponent(), Some(2.0));
        assert_eq!(parse("").expect("parses").lip_area_exponent(), None);
        for (extra, needle) in [
            ("lip_area_exponent\n", "missing field"),
            ("lip_area_exponent 2 3\n", "unexpected extra field"),
            ("lip_area_exponent x\n", "malformed field"),
            ("lip_area_exponent 0\n", "expected a positive number"),
            ("lip_area_exponent 0.5\n", "lip_area_exponent must be in [1, 4]"),
            ("lip_area_exponent 5\n", "lip_area_exponent must be in [1, 4]"),
            ("lip_area_exponent 2\nlip_area_exponent 3\n", "may only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        let build = |extra: &str| pressure(&LipDescription::from_bytes(fixture(extra).as_bytes()).expect("parses"));
        let plain = build("");
        let unit = build("lip_area_exponent 1\n");
        assert!(plain.iter().zip(&unit).all(|(a, b)| a.to_bits() == b.to_bits()), "q = 1 must be bit-identical to no row");
        let steeper = build("lip_area_exponent 2\n");
        assert_ne!(plain, steeper, "q = 2 must reach the radiated pressure");
        assert!(steeper.iter().all(|x| x.is_finite()));
        let (a, b) = (
            LipDescription::from_bytes(fixture("").as_bytes()).expect("parses"),
            LipDescription::from_bytes(fixture("lip_area_exponent 2\n").as_bytes()).expect("parses"),
        );
        assert_ne!(a.input_hash, b.input_hash);
        assert_eq!(a.lip.mass_kg.to_bits(), b.lip.mass_kg.to_bits(), "the lip's own motion is untouched");
    }

    // The digitized LOW tract of a trumpeter (Kaburagi et al. 2011, Fig. 3, Eb4) as card rows: 40 equal sections from the lips, then the
    // printed glottis. A tract invented for this test (five sections) was refused by the airway's own passivity check, so the card test
    // takes numbers the airway realizes.
    const TRACT_ROWS: &str = "player_tract_section 0.00435 1.42e-05\n\
        player_tract_section 0.00435 6.79e-05\n\
        player_tract_section 0.00435 0.0004216\n\
        player_tract_section 0.00435 0.0005289\n\
        player_tract_section 0.00435 0.0004784\n\
        player_tract_section 0.00435 0.0003205\n\
        player_tract_section 0.00435 0.0003016\n\
        player_tract_section 0.00435 0.0002763\n\
        player_tract_section 0.00435 0.0002511\n\
        player_tract_section 0.00435 0.0002574\n\
        player_tract_section 0.00435 0.0001658\n\
        player_tract_section 0.00435 0.0001089\n\
        player_tract_section 0.00435 5.84e-05\n\
        player_tract_section 0.00435 4.58e-05\n\
        player_tract_section 0.00435 3.63e-05\n\
        player_tract_section 0.00435 6.79e-05\n\
        player_tract_section 0.00435 5.84e-05\n\
        player_tract_section 0.00435 4.89e-05\n\
        player_tract_section 0.00435 4.26e-05\n\
        player_tract_section 0.00435 0.0001279\n\
        player_tract_section 0.00435 9.63e-05\n\
        player_tract_section 0.00435 9.95e-05\n\
        player_tract_section 0.00435 0.0001153\n\
        player_tract_section 0.00435 0.0001626\n\
        player_tract_section 0.00435 0.0002068\n\
        player_tract_section 0.00435 0.0003458\n\
        player_tract_section 0.00435 0.0003047\n\
        player_tract_section 0.00435 0.0002858\n\
        player_tract_section 0.00435 0.0002732\n\
        player_tract_section 0.00435 0.0002984\n\
        player_tract_section 0.00435 0.0004626\n\
        player_tract_section 0.00435 0.0004342\n\
        player_tract_section 0.00435 0.0003995\n\
        player_tract_section 0.00435 0.0002037\n\
        player_tract_section 0.00435 0.0001247\n\
        player_tract_section 0.00435 9.95e-05\n\
        player_tract_section 0.00435 7.74e-05\n\
        player_tract_section 0.00435 9.32e-05\n\
        player_tract_section 0.00435 9.95e-05\n\
        player_tract_section 0.00435 0.0001216\n\
        player_tract_section 1.2e-3 1.4e-5\n";

    #[test]
    fn player_tract_parses_is_refused_by_name_and_reaches_the_radiated_pressure() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        assert_eq!(parse("").expect("parses").player_tract_pressure(), None);
        // The rows stand in any order: the sections here come before the `player_tract` row.
        let ok = parse(&format!("{TRACT_ROWS}player_tract mouth 0.02\n")).expect("parses");
        assert_eq!(ok.player_tract_pressure(), Some(TractPressure::Mouth { tau_s: 0.02 }));
        assert_eq!(ok.player_tract_chain().expect("a chain").sections().len(), 41);
        assert_eq!(
            parse(&format!("player_tract lung\n{TRACT_ROWS}")).expect("parses").player_tract_pressure(),
            Some(TractPressure::Lung)
        );
        let one = "player_tract lung\nplayer_tract_section 0.01 0.0001\n";
        for (extra, needle) in [
            ("player_tract\n".to_string(), "missing field"),
            ("player_tract mouth\n".to_string(), "missing field"),
            ("player_tract mouth 0\n".to_string(), "expected a positive number"),
            ("player_tract mouth 0.02 1\n".to_string(), "unexpected extra field"),
            ("player_tract sideways\n".to_string(), "expected `mouth <tau_s>` or `lung`"),
            ("player_tract lung\nplayer_tract lung\n".to_string(), "may only appear once"),
            ("player_tract lung\n".to_string(), "at least one player_tract_section"),
            ("player_tract_section 0.01 0.0001\n".to_string(), "need a player_tract row"),
            ("player_tract_wall 1 1 1\n".to_string(), "need a player_tract row"),
            ("player_airway weibel-1963 saab-2021-table-1-B 3\n".to_string(), "need a player_tract row"),
            ("player_tract lung\nplayer_tract_section 0.01\n".to_string(), "missing field"),
            ("player_tract lung\nplayer_tract_section 0.01 -1\n".to_string(), "expected a positive number"),
            ("player_tract lung\nplayer_tract_section 0.01 0.0001 5\n".to_string(), "unexpected extra field"),
            (format!("{one}player_tract_wall 1 -1 0\n"), "nonnegative"),
            (format!("{one}player_tract_wall 1 1 1\nplayer_tract_wall 1 1 1\n"), "may only appear once"),
            (format!("{one}player_airway weibel-1963 other 3\n"), "takes `weibel-1963"),
            (format!("{one}player_airway weibel-1963 saab-2021-table-1-B 0\n"), "count outside its admitted range"),
            (format!("{one}player_airway weibel-1963 saab-2021-table-1-B 25\n"), "count outside its admitted range"),
            (format!("{one}jet_noise 7 1.3e-3 0.045 0.115 1700\n"), "not combined with jet_noise"),
        ] {
            let error = parse(&extra).expect_err(&extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        // The rendered pressure: the same file renders the same bits, the tract reaches the radiated pressure, and a file
        // without the rows is the original input (the hash covers the file; the lip's own numbers are untouched).
        let build = |extra: &str| pressure(&LipDescription::from_bytes(fixture(extra).as_bytes()).expect("parses"));
        let plain = build("");
        assert!(plain.iter().zip(&build("")).all(|(a, b)| a.to_bits() == b.to_bits()), "the same file must render the same bits");
        let lung = build(&format!("player_tract lung\n{TRACT_ROWS}"));
        assert!(lung.iter().all(|x| x.is_finite()));
        assert_ne!(plain, lung, "the tract must reach the radiated pressure");
        // The constriction noise rides on the tract as a series source: the pair renders, and differs from the tract alone and from the noise alone.
        let noisy = build(&format!("player_tract lung\n{TRACT_ROWS}constriction_noise 7 8.4e-3\n"));
        assert!(noisy.iter().all(|x| x.is_finite()));
        assert_ne!(noisy, lung, "the noise must reach the pressure through the tract");
        assert_ne!(noisy, build("constriction_noise 7 8.4e-3\n"), "the tract must reach the pressure through the noise");
        let (a, b) = (
            LipDescription::from_bytes(fixture("").as_bytes()).expect("parses"),
            LipDescription::from_bytes(fixture(&format!("player_tract lung\n{TRACT_ROWS}")).as_bytes()).expect("parses"),
        );
        assert_ne!(a.input_hash, b.input_hash);
        assert_eq!(a.lip.mass_kg.to_bits(), b.lip.mass_kg.to_bits(), "the lip's own motion numbers are untouched");
    }

    #[test]
    fn lip_rim_stop_parses_is_refused_by_name_reads_its_bound_from_the_first_cone_and_changes_nothing_it_does_not_reach() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        assert_eq!(parse("lip_rim_stop inelastic\n").expect("parses").lip_rim_stop(), Some(crate::brass_loop::RimStopKind::Inelastic));
        assert_eq!(parse("lip_rim_stop elastic\n").expect("parses").lip_rim_stop(), Some(crate::brass_loop::RimStopKind::Elastic));
        assert_eq!(parse("").expect("parses").lip_rim_stop(), None);
        for (extra, needle) in [
            ("lip_rim_stop\n", "missing field"),
            ("lip_rim_stop inelastic 2\n", "unexpected extra field"),
            ("lip_rim_stop sticky\n", "expected inelastic or elastic"),
            ("lip_rim_stop inelastic\nlip_rim_stop elastic\n", "may only appear once"),
            ("lip_response_mode2 244.5 11.0 -0.073\nlip_rim_stop inelastic\n", "refused with lip_response_mode2"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        // The bound is read from the card's own bore: twice the entry radius of the first cone (the fixture's cylinder, 6 mm).
        let performance = parse("lip_rim_stop inelastic\n").expect("parses").build().expect("builds");
        let info = performance.info().rim_stop.expect("the stop is reported");
        assert_eq!(info.kind, crate::brass_loop::RimStopKind::Inelastic);
        assert!((info.max_gap_m - 0.012).abs() < 1e-15, "{}", info.max_gap_m);
        assert_eq!(performance.rim_stop_hits(), Some(0));
        // A stop the lip never reaches changes nothing and counts nothing: the fixture's lip stays far under a 12 mm gap.
        let plain = pressure(&parse("").expect("parses"));
        let mut run = parse("lip_rim_stop elastic\n").expect("parses").build().expect("builds");
        let mut stopped = vec![0.0; 1200];
        for block in stopped.chunks_mut(37) {
            run.block(block).expect("physical pressure");
        }
        assert!(plain.iter().zip(&stopped).all(|(a, b)| a.to_bits() == b.to_bits()), "a stop that never acts must be bit-identical");
        assert_eq!(run.rim_stop_hits(), Some(0), "the fixture's lip never reaches a 12 mm rim");
        let (a, b) = (parse("").expect("parses"), parse("lip_rim_stop elastic\n").expect("parses"));
        assert_ne!(a.input_hash, b.input_hash);
    }

    /// Gauntlet G0 (a gate that never opens adds exactly zero), G3 (an open gate reaches the radiated sound, is seeded and
    /// finite) and a planted negative (an ungated source speaks earlier than the gated one).
    #[test]
    fn jet_noise_is_bit_identical_when_its_gate_never_opens_and_otherwise_gated_seeded_and_audible() {
        let build = |extra: &str| pressure(&LipDescription::from_bytes(fixture(extra).as_bytes()).expect("parses"));
        let plain = build("");
        let shut = build("jet_noise 7 1.3e-3 0.045 0.115 1e12\n");
        assert!(
            plain.iter().zip(&shut).all(|(a, b)| a.to_bits() == b.to_bits()),
            "a gate that never opens must leave the radiated pressure bit-identical to no row"
        );
        let open = build("jet_noise 7 1.3e-3 0.045 0.115 1700\n");
        assert_ne!(plain, open, "an open gate must reach the radiated pressure");
        let again = build("jet_noise 7 1.3e-3 0.045 0.115 1700\n");
        assert!(open.iter().zip(&again).all(|(a, b)| a.to_bits() == b.to_bits()), "a declared seed repeats bit for bit");
        assert_ne!(open, build("jet_noise 8 1.3e-3 0.045 0.115 1700\n"), "another seed gives another noise");
        // The fixture's phrase ramps the mouth pressure up from 5 ms, so the jet crosses the 1700 gate partway up the ramp,
        // while a source with no gate speaks from the first forward flow.
        let first = |x: &[f64]| plain.iter().zip(x).position(|(a, b)| a.to_bits() != b.to_bits());
        let ungated = build("jet_noise 7 1.3e-3 0.045 0.115 0\n");
        let (at_gate, at_flow) = (first(&open).expect("gated source speaks"), first(&ungated).expect("ungated source speaks"));
        assert!(at_flow + 20 < at_gate, "the gated source must wait for the gate: ungated at {at_flow}, gated at {at_gate}");
    }

    #[test]
    fn constriction_noise_parses_and_is_refused_by_name() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        let ok = parse("constriction_noise 7 8.4e-3\n").expect("parses");
        assert_eq!(ok.constriction_noise(), Some(&ConstrictionNoiseCard { seed: 7, duct_radius_m: 8.4e-3 }));
        assert_eq!(parse("").expect("parses").constriction_noise(), None);
        for (extra, needle) in [
            ("constriction_noise\n", "missing field"),
            ("constriction_noise 7\n", "missing field"),
            ("constriction_noise 7 8.4e-3 9\n", "unexpected extra field"),
            ("constriction_noise x 8.4e-3\n", "malformed field"),
            ("constriction_noise 7 nan\n", "must be finite"),
            ("constriction_noise 7 0\n", "expected a positive number"),
            ("constriction_noise 7 -8.4e-3\n", "expected a positive number"),
            ("constriction_noise 7 8.4e-3\nconstriction_noise 8 8.4e-3\n", "may only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    /// Gauntlet G0 (a cutoff under every driven band leaves nothing to add: bit-identical to no row), G3 (the printed law
    /// reaches the radiated sound, is seeded and finite) and the planted negative that keeps the first from being vacuous (the
    /// same row with the cup radius of a real cup speaks).
    #[test]
    fn constriction_noise_is_bit_identical_with_no_band_below_the_cutoff_and_otherwise_seeded_and_audible() {
        let build = |extra: &str| pressure(&LipDescription::from_bytes(fixture(extra).as_bytes()).expect("parses"));
        let plain = build("");
        // A 12 m "duct" has its cutoff `100 / r` at 8.3 Hz, under the bank's lowest band (10 Hz), so nothing is driven at any
        // pressure. (A 3 m radius does not do it: the Strouhal length grows with the radius and the 10-31.5 Hz bands are driven.)
        let cut_off = build("constriction_noise 7 12\n");
        assert!(
            plain.iter().zip(&cut_off).all(|(a, b)| a.to_bits() == b.to_bits()),
            "a row whose cutoff leaves no driven band must leave the radiated pressure bit-identical to no row"
        );
        let open = build("constriction_noise 7 8.4e-3\n");
        assert!(open.iter().all(|p| p.is_finite()));
        assert_ne!(plain, open, "the printed law drives bands on a real cup, so it must reach the radiated pressure");
        let again = build("constriction_noise 7 8.4e-3\n");
        assert!(open.iter().zip(&again).all(|(a, b)| a.to_bits() == b.to_bits()), "a declared seed repeats bit for bit");
        assert_ne!(open, build("constriction_noise 8 8.4e-3\n"), "another seed gives another noise");
    }

    fn with_gas(record: &str) -> String {
        fixture("").replace("moist_air 293.15 101325 0", record)
    }

    fn pressure(description: &LipDescription) -> Vec<f64> {
        let mut performance = description.build().expect("gas reaches all physics owners");
        assert_eq!(
            performance.info().sound_speed_m_s,
            description.ambient_gas().sound_speed
        );
        let mut result = vec![0.0; 1200];
        for block in result.chunks_mut(37) {
            performance.block(block).expect("physical pressure");
        }
        assert!(result.iter().all(|value| value.is_finite()));
        assert!(result.iter().any(|value| value.abs() > 1e-12));
        result
    }

    #[test]
    fn g3_shared_mixture_composition_and_pressure_change_the_radiated_lip_sound() {
        let legacy = LipDescription::from_bytes(fixture("").as_bytes()).unwrap();
        assert_eq!(
            *legacy.ambient_gas(),
            GasState::try_new_moist_air(293.15, 101325.0, 0.0).unwrap()
        );
        let records = [
            "mixture 293.15 101325 mole 2 N2 0.8 O2 0.2",
            "mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2",
            "mixture 293.15 80000 mole 2 N2 0.8 O2 0.2",
        ];
        let mut rendered = Vec::new();
        for record in records {
            let description = LipDescription::from_bytes(with_gas(record).as_bytes()).unwrap();
            let expected =
                read_gas_record(&mut record.split_ascii_whitespace(), None, None).unwrap();
            assert_eq!(description.gas_record(), &expected);
            rendered.push(pressure(&description));
        }
        assert_ne!(
            rendered[0], rendered[1],
            "composition must reach the audible loop"
        );
        assert_ne!(
            rendered[0], rendered[2],
            "pressure must reach flow and radiation"
        );
    }

    #[test]
    fn g3_nasa9_source_swap_changes_lip_pressure_and_retains_the_original_packs() {
        use fs_material::gas_species::input::ThermochemicalPhase;
        let record = "thermochemical 293.15 101325 mole 1 specimen-gas 1";
        let text = with_gas(record);
        assert!(
            LipDescription::from_bytes(text.as_bytes())
                .unwrap_err()
                .contains("explicit source registry")
        );
        let mut rendered = Vec::new();
        for cp_r in [3.5, 4.5] {
            let packs = gas_test_support::packs("specimen-gas", 0.028, cp_r, 200.0);
            let bindings = [packs.binding()];
            let sources = ThermochemicalSources {
                species: &bindings,
                water_species: None,
                phase: ThermochemicalPhase::FrozenGas,
            };
            let description =
                LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources))
                    .unwrap();
            let expected =
                read_gas_record(&mut record.split_ascii_whitespace(), None, Some(&sources))
                    .unwrap();
            assert_eq!(description.gas_record(), &expected);
            let GasRecord::Thermochemical(gas) = description.gas_record() else {
                panic!("retained source");
            };
            assert_eq!(gas.source_packs()[0].nasa9_pack, packs.nasa9.content_hash());
            rendered.push(pressure(&description));
        }
        assert_ne!(
            rendered[0], rendered[1],
            "data-only heat-capacity swap must reach sound"
        );
    }

    #[test]
    fn gas_refusals_keep_source_admission_strict() {
        for record in [
            "mixture 293.15 101325 mole 2 N2 0.7 O2 0.2",
            "mixture 293.15 101325 mole 1 unavailable 1",
            "mixture 293.15 101325 mole 1 N2 1 extra",
            "moist_air 293.15 101325 0 extra",
            "thermochemical 293.15 101325 mole 1 specimen-gas 1",
        ] {
            assert!(
                LipDescription::from_bytes(with_gas(record).as_bytes()).is_err(),
                "{record}"
            );
        }
    }

    #[test]
    fn bessel_pieces_lie_on_the_horn_and_reach_both_ends() {
        let (r_in, r_out, length, flare) = (0.01105, 0.01507, 0.15042, 1.26);
        let mut segments = Vec::new();
        bessel_horn(r_in, r_out, length, flare, 16, &mut segments).expect("horn");
        assert_eq!(segments.len(), 16);
        let mut radii = vec![r_in];
        let mut total = 0.0;
        for segment in &segments {
            let Segment::Cone {
                inlet_radius,
                outlet_radius,
                length,
            } = *segment
            else {
                panic!("a Bessel piece is a cone");
            };
            assert_eq!(inlet_radius, *radii.last().expect("previous"));
            assert!(outlet_radius > inlet_radius);
            radii.push(outlet_radius);
            total += length;
        }
        assert!((total - length).abs() < 1e-15);
        assert_eq!(*radii.last().expect("end"), r_out);
        // r = b (u - s)^-flare, so r^(-1/flare) is affine in s: equal steps
        // in s give equal steps in r^(-1/flare).
        let affine: Vec<f64> = radii.iter().map(|r| det::pow(*r, -1.0 / flare)).collect();
        let first = affine[1] - affine[0];
        for pair in affine.windows(2) {
            assert!(((pair[1] - pair[0]) - first).abs() < 1e-9 * first.abs());
        }
        assert!(bessel_horn(r_out, r_in, length, flare, 4, &mut segments).is_err());
        assert!(bessel_horn(r_in, r_out, length, 0.0, 4, &mut segments).is_err());
    }

    #[test]
    fn lip_response_reproduces_the_published_quantities() {
        let (f, q, inverse_mu, area) = (204.5, 7.0, 0.063, 1.0e-4);
        let text = fixture("").replace(
            "lip_response 200 7 0.1 0.012 0.0005 0.0001",
            &format!("lip_response {f} {q} {inverse_mu} 0.012 0.00053 {area}"),
        );
        let description = LipDescription::from_bytes(text.as_bytes()).expect("parses");
        let lip = &description.lip;
        let resonance = (lip.stiffness_n_m / lip.mass_kg).sqrt() / core::f64::consts::TAU;
        let quality = lip.mass_kg * core::f64::consts::TAU * f / lip.damping_n_s_m;
        assert!((resonance - f).abs() < 1e-9 * f);
        assert!((quality - q).abs() < 1e-9 * q);
        assert!((lip.face_area_m2 / lip.mass_kg - inverse_mu).abs() < 1e-12);
        assert_eq!(lip.provenance, "test fixture, not a measurement");
    }

    #[test]
    fn the_mouth_flow_tail_tolerance_defaults_to_the_ladders_own_and_is_declarable() {
        let default = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        assert_eq!(default.mouth_flow_tail_tol, LineLadder::default().tail_tol);
        let declared = LipDescription::from_bytes(fixture("mouth_flow_tail_tol 1e-6\n").as_bytes())
            .expect("parses");
        assert_eq!(declared.mouth_flow_tail_tol, 1.0e-6);
    }

    #[test]
    fn the_mouth_flow_tap_ceiling_defaults_to_the_ladders_own_and_is_declarable() {
        let default = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        assert_eq!(default.mouth_flow_max_taps, MouthFlowConfig::default().max_taps);
        assert_eq!(default.mouth_flow_max_taps, 65_536);
        for taps in [512_usize, 65_536, 131_072, 262_144] {
            let declared = LipDescription::from_bytes(fixture(&format!("mouth_flow_max_taps {taps}\n")).as_bytes())
                .expect("parses");
            assert_eq!(declared.mouth_flow_max_taps, taps);
        }
        for (extra, needle) in [
            ("mouth_flow_max_taps 100000\n", "power of two"),
            ("mouth_flow_max_taps 256\n", "power of two"),
            ("mouth_flow_max_taps 524288\n", "outside its admitted range"),
            ("mouth_flow_max_taps 0\n", "outside its admitted range"),
            ("mouth_flow_max_taps 131072\nmouth_flow_max_taps 131072\n", "only appear once"),
        ] {
            let error = LipDescription::from_bytes(fixture(extra).as_bytes()).unwrap_err();
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    #[test]
    fn modal_baffle_is_explicit_and_cannot_replace_a_foreign_exterior_or_nonlinear_bore() {
        assert!(LipDescription::from_bytes(fixture("").as_bytes()).unwrap().modal_baffle.is_none());
        let selected = "modal_baffle 1e-5 8000000\n";
        let source = fixture(selected);
        assert!(LipDescription::from_bytes(source.as_bytes()).unwrap().modal_baffle.is_some());
        for foreign in ["unflanged", "ideal_open"] {
            let source = source.replace("termination flanged", &format!("termination {foreign}"));
            assert!(LipDescription::from_bytes(source.as_bytes()).unwrap_err().contains("modal_baffle requires"));
        }
        let nonlinear = fixture(&format!("{selected}bore nonlinear 0.0024 4 1 24 20 8000 1e-2\n"));
        assert!(LipDescription::from_bytes(nonlinear.as_bytes()).unwrap_err().contains("modal_baffle requires"));
        for declaration in ["modal_baffle 0 8000000\n", "modal_baffle 1 8000000\n",
            "modal_baffle 1e-5 0\n", "modal_baffle 1e-5 134217729\n",
            "modal_baffle 1e-5 8000000\nmodal_baffle 1e-5 8000000\n"] {
            assert!(LipDescription::from_bytes(fixture(declaration).as_bytes()).is_err());
        }
    }

    // The woodwind exhaled-breath bore gas the clarinet and oboe cards cite (Coltman 1966, Nederveen 1969, Fletcher and
    // Rossing 1998): 27 C, 101325 Pa, saturated, 4.5 % dry CO2, 16.5 % dry O2. Here as a test input, not as a claim for a brass bore.
    const BREATH: &str = "bore_gas mixture 300.15 101325 mole 5 N2 0.7534 O2 0.1591 Ar 0.0090 CO2 0.0434 H2O 0.0351\n";
    const MODAL: &str = "modal_baffle 1e-5 8000000\n";
    // The fixture's modal mouth output does not meet the default tail tolerance (1e-10) even at 65536 taps (tail energy 5e-7 on the
    // second and third stage), so a render declares one; the Bb3 page card declares 2e-6, as here.
    const TAIL: &str = "mouth_flow_tail_tol 2e-6\n";

    #[test]
    fn bore_gas_defaults_to_the_room_and_is_refused_by_name() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        let plain = parse(MODAL).expect("parses");
        assert!(plain.bore_gas_record().is_none());
        assert_eq!(plain.bore_gas(), plain.ambient_gas(), "without the row the bore holds the room gas");
        let separated = parse(&format!("{MODAL}{BREATH}")).expect("parses");
        assert_eq!(separated.bore_gas().temperature, 300.15);
        assert_eq!(separated.ambient_gas().temperature, 293.15, "the room stays what ambient says");
        assert!(separated.bore_gas().sound_speed != separated.ambient_gas().sound_speed);
        assert!(separated.bore_gas_record().is_some());
        // The optional legacy keyword and the legacy triple read as they do on `ambient`.
        let legacy = parse(&format!("{MODAL}bore_gas moist_air 300.15 101325 1\n")).expect("parses");
        assert_eq!(legacy.bore_gas().temperature, 300.15);
        for (extra, needle) in [
            (BREATH.to_string(), "bore_gas needs modal_baffle"),
            (format!("{MODAL}{BREATH}{BREATH}"), "may only appear once"),
            (format!("{MODAL}bore_gas mixture 300.15 90000 mole 1 N2 1\n"), "one static pressure"),
            (format!("{MODAL}bore_gas\n"), "bore_gas refused"),
            (format!("{MODAL}bore_gas mixture 300.15 101325 mole 1 unobtainium 1\n"), "bore_gas refused"),
            (format!("{MODAL}bore_gas moist_air 300.15 101325 1 7\n"), "unexpected extra field"),
        ] {
            let error = parse(&extra).expect_err(&extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    #[test]
    fn bore_gas_equal_to_the_room_changes_no_bit_and_a_breath_bore_scales_the_resonances_by_the_sound_speed() {
        let scan = |extra: &str, frequencies: &[f64]| {
            LipDescription::from_bytes(fixture(extra).as_bytes()).unwrap().bore_scan(frequencies).unwrap()
        };
        let probe = [95.0, 131.0, 262.5, 411.0];
        let plain = scan(MODAL, &probe);
        let same = scan(&format!("{MODAL}bore_gas moist_air 293.15 101325 0\n"), &probe);
        assert_eq!(plain.input_characteristic_impedance.to_bits(), same.input_characteristic_impedance.to_bits());
        assert_eq!(plain.density_kg_m3.to_bits(), same.density_kg_m3.to_bits());
        for (a, b) in plain.rows.iter().zip(&same.rows) {
            for (x, y) in a.impedance_multimodal.iter().chain(&a.impedance_scalar).zip(b.impedance_multimodal.iter().chain(&b.impedance_scalar)) {
                assert_eq!(x.to_bits(), y.to_bits());
            }
        }
        // A breath bore: the input reference is the BORE gas's rho c / S; the pressure at the microphone is read in the ROOM's density.
        let breath = scan(&format!("{MODAL}{BREATH}"), &probe);
        let room = LipDescription::from_bytes(fixture(MODAL).as_bytes()).unwrap();
        let bore = LipDescription::from_bytes(fixture(&format!("{MODAL}{BREATH}")).as_bytes()).unwrap();
        let (r, c_room, c_bore) = (0.006_f64, room.ambient_gas().sound_speed, bore.bore_gas().sound_speed);
        assert_eq!(breath.input_characteristic_impedance, bore.bore_gas().density * c_bore / (core::f64::consts::PI * r * r));
        assert_eq!(breath.density_kg_m3.to_bits(), plain.density_kg_m3.to_bits());
        // The first input-impedance maximum of the passive bore moves by the sound-speed ratio (the end correction and the
        // room's radiation load do not depend on the bore gas; the wall loss shifts it a little).
        let peak = |extra: &str| {
            let description = LipDescription::from_bytes(fixture(extra).as_bytes()).unwrap();
            let coarse: Vec<f64> = (0..400).map(|i| 60.0 + 0.5 * f64::from(i)).collect();
            let magnitude = |scan: &BoreScan| scan.rows.iter().map(|row| row.impedance_multimodal[0].hypot(row.impedance_multimodal[1])).collect::<Vec<_>>();
            // The first local maximum above the scan's start: a narrow bore's |Z| may also be large toward DC.
            let first_peak = |values: &[f64]| (1..values.len() - 1).find(|&i| values[i - 1] < values[i] && values[i] >= values[i + 1]).expect("a peak");
            let scan = description.bore_scan(&coarse).unwrap();
            let centre = coarse[first_peak(&magnitude(&scan))];
            let fine: Vec<f64> = (0..=240).map(|i| centre - 3.0 + 0.025 * f64::from(i)).collect();
            let scan = description.bore_scan(&fine).unwrap();
            fine[first_peak(&magnitude(&scan))]
        };
        let (f_room, f_bore) = (peak(MODAL), peak(&format!("{MODAL}{BREATH}")));
        let ratio = f_bore / f_room;
        let expected = c_bore / c_room;
        assert!((ratio / expected - 1.0).abs() < 2e-3, "peak ratio {ratio} against sound-speed ratio {expected} ({f_room} Hz -> {f_bore} Hz)");
        assert!(ratio > 1.0 + 0.5 * (expected - 1.0), "a breath bore must move the peak, not leave it");
    }

    #[test]
    fn a_bore_gas_render_runs_the_bank_the_mouth_flow_and_the_lip_in_the_bore_gas_and_the_receiver_in_the_room() {
        let separated = LipDescription::from_bytes(fixture(&format!("{MODAL}{TAIL}{BREATH}")).as_bytes()).unwrap();
        let mut performance = separated.build().unwrap();
        let info = performance.info();
        assert_eq!(info.gas_temperature_k, 300.15);
        assert_eq!(info.sound_speed_m_s, separated.bore_gas().sound_speed);
        assert!(info.sound_speed_m_s != separated.ambient_gas().sound_speed);
        let mut output = vec![0.0; 600];
        for block in output.chunks_mut(61) {
            performance.block(block).unwrap();
        }
        assert!(output.iter().all(|p| p.is_finite()));
        assert!(output.iter().any(|p| p.abs() > 1e-9), "the loop must sound");
        // Without the row the same file is the room-gas render it always was.
        let plain = LipDescription::from_bytes(fixture(&format!("{MODAL}{TAIL}")).as_bytes()).unwrap().build().unwrap();
        assert_eq!(plain.info().gas_temperature_k, 293.15);
    }

    #[test]
    fn modal_baffle_example_renders_all_mouth_modes_from_the_declared_flanged_bore() {
        let description = LipDescription::from_bytes(
            include_bytes!("../../examples/modal-baffled-lip.performance")).unwrap();
        let gas = description.ambient_gas();
        let frequency = 2107.0;
        let scan = description.bore_scan(&[frequency]).unwrap();
        let expected = fs_duct::modal::mm_input_impedance_flanged(&description.duct, gas, gas,
            core::f64::consts::TAU*frequency, LossModel::Bessel, description.modes,
            description.extra_slices).unwrap();
        assert_eq!(scan.rows[0].impedance_multimodal,
            [expected.plane_impedance.re, expected.plane_impedance.im]);
        let mut performance = description.build().unwrap();
        let LipBoreInfo::Linear(lines) = &performance.info().bore else { panic!("linear fixture"); };
        assert_eq!(lines.mouth_flow_modes, 3);
        assert_eq!(lines.modal_radiation.as_ref().unwrap().modes, 3);
        assert_eq!(lines.modal_radiation.as_ref().unwrap().radial_nodes, [16,32,64]);
        let mouth = match &performance.engine {
            Engine::Linear { mouth, .. } => mouth,
            #[cfg(feature = "nonlinear-brass")]
            Engine::Nonlinear(_) => panic!("linear fixture"),
        };
        assert_eq!(mouth.len(), 3);
        assert!(mouth[1..].iter().all(|filter|
            filter.quiescent_fir().unwrap().1.iter().any(|tap| *tap != 0.0)));
        let mut output = vec![0.0; description.samples as usize];
        let mut contributions = [0.0; 3];
        let mut higher_peak = 0.0_f64;
        for block in output.chunks_mut(61) {
            performance.block(block).unwrap();
            performance.receiver.mode_pressures(&mut contributions).unwrap();
            higher_peak = higher_peak.max(contributions[1].abs()).max(contributions[2].abs());
        }
        assert!(output.iter().all(|p| p.is_finite()));
        assert!(output.iter().any(|p| p.abs() > 1e-9));
        assert!(higher_peak > 1e-12, "retained higher mouth modes must reach emitted pressure");
        assert_eq!(performance.samples_rendered(), description.samples);
    }

    #[test]
    fn control_ticks_apply_at_the_ceiling_of_their_start() {
        let (control, rate) = (700, 19_200);
        assert_eq!(tick_at(0, control, rate), 0);
        assert_eq!(tick_start(1, control, rate), 28);
        assert_eq!(tick_at(27, control, rate), 0);
        assert_eq!(tick_at(28, control, rate), 1);
        for tick in 0..2_000u64 {
            let start = tick_start(tick, control, rate);
            assert_eq!(tick_at(start, control, rate), tick);
            if start > 0 {
                assert_eq!(tick_at(start - 1, control, rate), tick - 1);
            }
        }
    }

    #[test]
    fn malformed_and_non_physical_files_are_refused_by_name() {
        for (extra, needle) in [
            ("bogus 1\n", "unknown key"),
            ("audio 24000 10\n", "only appear once"),
            ("cylinder 0.006\n", "missing field"),
            ("cone 0.006 0.012 0.1 0.1\n", "extra field"),
            ("microphone 0 0 -1 4 16 1000\n", "in front of the baffle"),
            ("lip 1e-4 100 0.1 0.01 0.0005 0.0001\n", "only appear once"),
            (
                "lip 1e-4 100 -1 0.01 0.0005 0.0001\n",
                "damping must be nonnegative",
            ),
            ("modes 0 1\n", "outside its admitted range"),
            ("mouth_flow_tail_tol 1\n", "below 1"),
            ("mouth_flow_tail_tol 0\n", "positive"),
            (
                "mouth_flow_tail_tol 1e-6\nmouth_flow_tail_tol 1e-6\n",
                "only appear once",
            ),
            ("bessel 0.01 0.02 0.1 1 0\n", "outside its admitted range"),
            ("bore linear 0.0024\n", "expected nonlinear"),
            (
                "bore nonlinear 0.0024 17 1 24 20 8000 1e-2\n",
                "outside its admitted range",
            ),
            (
                "bore nonlinear 0.0024 16 1 129 20 8000 1e-2\n",
                "outside its admitted range",
            ),
            (
                "bore nonlinear 0.0024 16 1 24 8000 20 1e-2\n",
                "rising wall band",
            ),
            ("bore nonlinear 0.0024 16 1 24 20 8000 1\n", "below 1"),
        ] {
            let error = LipDescription::from_bytes(fixture(extra).as_bytes()).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        let no_bore_source =
            fixture("").replace("bore_source test fixture, not a measurement\n", "");
        let error = LipDescription::from_bytes(no_bore_source.as_bytes()).expect_err("source");
        assert!(error.contains("bore_source"), "{error}");
        let no_source_text =
            fixture("").replace("lip_source test fixture, not a measurement", "lip_source");
        let error = LipDescription::from_bytes(no_source_text.as_bytes()).expect_err("text");
        assert!(error.contains("provenance"), "{error}");
        let noncanonical = fixture("").replace("1e-2", "0.01");
        assert!(LipDescription::from_bytes(noncanonical.as_bytes()).is_err());
        assert!(LipDescription::from_bytes(b"not a lip file\n").is_err());
    }

    #[test]
    fn window_diagnostics_cover_the_render_and_describe_a_possible_lip() {
        let description = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        assert!(performance.windows().is_empty());
        let mut block = vec![0.0; 600];
        for _ in 0..4 {
            performance.block(&mut block).expect("block");
        }
        let windows = performance.windows();
        // 2400 samples at 24 kHz are exactly one tenth of a second.
        assert_eq!(windows.len(), 1);
        let w = windows[0];
        assert_eq!(w.start_s, 0.0);
        assert!(w.gap_min_m >= 0.0 && w.gap_min_m <= w.gap_mean_m && w.gap_mean_m <= w.gap_max_m);
        assert!((0.0..=1.0).contains(&w.closed_chunk_fraction));
        assert!(w.bore_ac_rms_pa.is_finite() && w.bore_ac_rms_pa <= w.bore_ac_peak_pa + 1e-9);
        assert!(w.bore_mean_pa.is_finite());
        assert_eq!(w.lip_tension_work_j, [0.0; 2]);
        assert_eq!(performance.lip_tension_work_j(), [0.0; 2]);
    }

    #[test]
    fn scheduled_lip_tension_work_reaches_windows_once_independent_of_block_partition() {
        use fs_scenario::gesture::{GestureEvent, GestureTrack, GestureValue};
        let schedule = GestureSchedule::try_new(
            700,
            vec![
                GestureTrack {
                    id: "blow".into(),
                    target: GestureTarget::BlowingPressure,
                    initial: GestureValue::PressurePa(2800.0),
                    events: Vec::new(),
                },
                GestureTrack {
                    id: "tension".into(),
                    target: GestureTarget::PreStress,
                    initial: GestureValue::Multiplier(1.0),
                    events: [(0.02, 1.2), (0.10, 0.8), (0.16, 1.1)]
                        .into_iter()
                        .map(|(time_s, value)| GestureEvent {
                            time_s,
                            transition_s: 0.0,
                            value: GestureValue::Multiplier(value),
                        })
                        .collect(),
                },
            ],
        )
        .unwrap();
        let schedule = String::from_utf8(schedule.to_canonical_bytes()).unwrap();
        let text = fixture("")
            .replace("audio 24000 2400", "audio 24000 4800")
            .replace(SCHEDULE, &schedule);
        let description = LipDescription::from_bytes(text.as_bytes()).unwrap();
        let render = |chunk: usize| {
            let mut performance = description.build().unwrap();
            let mut pressure = vec![0.0; 4800];
            for block in pressure.chunks_mut(chunk) {
                performance.block(block).unwrap();
            }
            let changes: Vec<_> = performance
                .engine
                .diagnostics()
                .iter()
                .map(|d| d.lip_tension_work_j)
                .filter(|w| *w != [0.0; 2])
                .collect();
            (
                pressure,
                performance.lip_tension_work_j(),
                performance.windows(),
                changes,
            )
        };
        let (a, work_a, windows_a, changes_a) = render(37);
        let (b, work_b, windows_b, changes_b) = render(600);
        assert_eq!(
            a, b,
            "accounting must preserve the original sample trajectory"
        );
        assert_eq!(work_a, work_b);
        assert_eq!(changes_a, changes_b);
        assert_eq!(changes_a.len(), 3, "each control event is debited once");
        assert!(changes_a[0][0] > 0.0 && changes_a[1][0] < 0.0 && changes_a[2][0] > 0.0);
        assert!(changes_a.iter().all(|w| w[1] == 0.0));
        assert_eq!(windows_a.len(), 2);
        assert_eq!(windows_b.len(), 2);
        for (a, b) in windows_a.iter().zip(&windows_b) {
            assert_eq!(a.start_s, b.start_s);
            assert_eq!(a.lip_tension_work_j, b.lip_tension_work_j);
        }
        assert_eq!(windows_a[0].lip_tension_work_j, changes_a[0]);
        let in_second = changes_a[1][0] + changes_a[2][0];
        assert_eq!(windows_a[1].lip_tension_work_j, [in_second, 0.0]);
        let from_windows = windows_a
            .iter()
            .map(|w| w.lip_tension_work_j[0])
            .sum::<f64>();
        let scale = changes_a.iter().map(|w| w[0].abs()).sum::<f64>();
        assert!((from_windows - work_a[0]).abs() <= 4.0 * f64::EPSILON * scale);
    }

    #[test]
    fn the_mouthpiece_tap_records_every_rendered_sample_only_when_asked() {
        let description = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        assert!(performance.mouthpiece_pressure().is_empty());
        let mut block = vec![0.0; 600];
        performance.block(&mut block).expect("block");
        assert!(performance.mouthpiece_pressure().is_empty());
        performance.record_mouthpiece();
        for _ in 0..3 {
            performance.block(&mut block).expect("block");
        }
        let tap = performance.mouthpiece_pressure();
        assert_eq!(tap.len(), 1800);
        assert!(tap.iter().all(|p| p.is_finite()));
    }

    #[test]
    fn the_bore_scan_reads_the_passive_bore_without_building_the_loop() {
        let description = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let scan = description
            .bore_scan(&[150.0, 300.0, 450.0])
            .expect("scans");
        assert_eq!(scan.rows.len(), 3);
        assert!(scan.input_characteristic_impedance > 0.0 && scan.microphone_distance_m > 0.0);
        for row in &scan.rows {
            let z = row.impedance_multimodal[0].hypot(row.impedance_multimodal[1]);
            assert!(z.is_finite() && z > 0.0);
            // a passive load: the real part of both input impedances is not negative
            assert!(row.impedance_multimodal[0] >= -1e-6 * z);
            assert!(row.impedance_scalar[0] >= -1e-6 * z);
            assert!(
                row.mouth_over_input_flow[0]
                    .hypot(row.mouth_over_input_flow[1])
                    .is_finite()
            );
        }
    }

    #[test]
    fn the_chain_builds_and_renders_finite_pressure_on_its_own_clock() {
        let description = LipDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        assert_eq!(performance.info().samples, 2400);
        let LipBoreInfo::Linear(lines) = &performance.info().bore else {
            panic!("a file without a bore record realizes the lines");
        };
        assert!(lines.mouth_flow_taps >= 256);
        // The microphone is fed by the same retained modal bore as the lip,
        // not a separately propagated scalar flow through this flare.
        {
            use fs_math::c64::C64;
            let (voice, mouth) = match &performance.engine {
                Engine::Linear { voice, mouth } => (voice, mouth),
                #[cfg(feature = "nonlinear-brass")]
                Engine::Nonlinear(_) => panic!("linear fixture"),
            };
            let reference = voice.port_impedance();
            let gas = description.ambient_gas();
            let geometry = fs_duct::modal::PreparedModalDuct::new(
                &description.duct, description.modes, description.extra_slices,
            ).unwrap();
            let (_, taps) = mouth[0].quiescent_fir().unwrap();
            let (mut scale, mut error, mut scalar_difference) = (0.0_f64, 0.0_f64, 0.0_f64);
            for frequency in [310.0, 1013.0, 2207.0, 5119.0] {
                let omega = core::f64::consts::TAU*frequency;
                let waves = geometry.mouth_waves(gas, omega, LossModel::Bessel,
                    description.termination).unwrap();
                let expected = waves.per_incident_plane_pressure(reference).unwrap()[0];
                assert_eq!(performance.info().mouth_radius_m, waves.mouth_radius_m);
                let realized = taps.iter().enumerate().fold(C64::ZERO, |sum, (sample, &h)| {
                    let phase = omega*sample as f64/f64::from(description.rate_hz);
                    sum+C64::new(phase.cos(), phase.sin()).scale(h)
                });
                let scalar = crate::driving_point::outlet::mouth_flow_per_incident_pressure(
                    &description.duct, gas, omega, description.termination, None, reference,
                ).unwrap();
                scale = scale.max(expected.abs());
                error = error.max((realized-expected).abs());
                scalar_difference = scalar_difference.max((scalar-expected).abs());
            }
            assert!(error < 1e-3*scale, "modal mouth FIR error {error:e}, scale {scale:e}");
            assert!(scalar_difference > 20.0*error,
                "modal/scalar distinction {scalar_difference:e} must exceed realization error {error:e}");
        }
        assert!(performance.validate_sample_rate(48_000).is_err());
        assert!(performance.validate_sample_rate(24_000).is_ok());
        assert!(performance.validate_sample_count(2401).is_err());
        let mut block = vec![0.0; 600];
        for _ in 0..4 {
            performance.block(&mut block).expect("block");
            assert!(block.iter().all(|p| p.is_finite()));
        }
        assert_eq!(performance.samples_rendered(), 2400);
        assert!(performance.block(&mut block).is_err());
    }

    const NONLINEAR_RECORD: &str = "bore nonlinear 0.0024 16 1 64 20 8000 1e-2\n";

    const FINITE_MOUTH_RECORD: &str = "mouth_reservoir 2e-5 0 106325 293.15\n";

    #[test]
    fn finite_mouth_record_requires_a_physical_nonlinear_owner_and_admitted_scalars() {
        assert_eq!(LipDescription::from_bytes(fixture("").as_bytes()).unwrap().mouth_reservoir(), None);
        let error = LipDescription::from_bytes(fixture(FINITE_MOUTH_RECORD).as_bytes()).unwrap_err();
        assert!(error.contains("mouth_reservoir") && error.contains("bore nonlinear"), "{error}");
        for (row, reason) in [
            ("mouth_reservoir 0 0 106325 293.15\n", "positive"),
            ("mouth_reservoir 2e-5 -1 106325 293.15\n", "non-negative"),
            ("mouth_reservoir 2e-5 0 0 293.15\n", "positive"),
            ("mouth_reservoir 2e-5 0 106325 NaN\n", "finite"),
            ("mouth_reservoir 2e-5 0 106325\n", "missing field"),
            ("mouth_reservoir 2e-5 0 106325 293.15 7\n", "extra field"),
        ] {
            let error = LipDescription::from_bytes(fixture(row).as_bytes()).unwrap_err();
            assert!(error.contains("mouth_reservoir") && error.contains(reason), "{row}: {error}");
        }
        let duplicate = format!("{FINITE_MOUTH_RECORD}{FINITE_MOUTH_RECORD}");
        let error = LipDescription::from_bytes(fixture(&duplicate).as_bytes()).unwrap_err();
        assert!(error.contains("may only appear once"), "{error}");
    }

    // The series noise rows stand beside the finite mouth in either row order. Without the nonlinear engine the description
    // is still refused for its engine, never for the noise.
    #[cfg(not(feature = "nonlinear-brass"))]
    #[test]
    fn finite_mouth_record_beside_series_noise_is_refused_only_for_its_engine_without_the_feature() {
        for noise in [
            "jet_noise 7 1.3e-3 0.045 0.115 1700\n",
            "constriction_noise 7 0.006\n",
        ] {
            for extra in [
                format!("{NONLINEAR_RECORD}{FINITE_MOUTH_RECORD}{noise}"),
                format!("{noise}{FINITE_MOUTH_RECORD}{NONLINEAR_RECORD}"),
            ] {
                let error = LipDescription::from_bytes(fixture(&extra).as_bytes()).unwrap_err();
                assert!(error.contains("nonlinear-brass") && !error.contains("work"), "{error}");
            }
        }
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn finite_mouth_record_admits_series_noise_in_either_row_order_and_the_lip_renders() {
        for noise in [
            "jet_noise 7 1.3e-3 0.045 0.115 1700\n",
            "constriction_noise 7 0.006\n",
        ] {
            for row in [format!("{FINITE_MOUTH_RECORD}{noise}"), format!("{noise}{FINITE_MOUTH_RECORD}")] {
                let description = finite_description(&row, "3e3");
                assert!(description.mouth_reservoir().is_some());
                let mut performance = description.build().unwrap();
                assert!(performance.finite_mouth_state().is_some());
                performance.apply_controls(0).unwrap();
                performance.record_mouthpiece();
                let mut radiated = [0.0; 32];
                performance.block(&mut radiated).unwrap();
                assert!(radiated.iter().all(|p| p.is_finite()));
            }
        }
    }

    const MOUTH_LOAD_RECORD: &str = "mouth_load modal 3 1 3000 100 1e-8\n";

    #[test]
    fn mouth_load_record_requires_bore_nonlinear_and_admitted_scalars() {
        assert_eq!(LipDescription::from_bytes(fixture("").as_bytes()).unwrap().mouth_load(), None);
        let error = LipDescription::from_bytes(fixture(MOUTH_LOAD_RECORD).as_bytes()).unwrap_err();
        assert!(error.contains("mouth_load") && error.contains("bore nonlinear"), "{error}");
        for (row, reason) in [
            ("mouth_load flat 3 1 3000 100 1e-8\n", "expected modal"),
            ("mouth_load modal 0 1 3000 100 1e-8\n", "count outside"),
            ("mouth_load modal 9 1 3000 100 1e-8\n", "count outside"),
            ("mouth_load modal 3 0 3000 100 1e-8\n", "count outside"),
            ("mouth_load modal 3 1 150 100 1e-8\n", "two bins"),
            ("mouth_load modal 3 1 3000 100 1\n", "below 1"),
            ("mouth_load modal 3 1 3000 -100 1e-8\n", "positive"),
            ("mouth_load modal 3 1 3000 100\n", "missing field"),
            ("mouth_load modal 3 1 3000 100 1e-8 7\n", "extra field"),
        ] {
            let error = LipDescription::from_bytes(fixture(&format!("{NONLINEAR_RECORD}{row}")).as_bytes()).unwrap_err();
            assert!(error.contains(reason), "{row}: {error}");
        }
        let duplicate = format!("{NONLINEAR_RECORD}{MOUTH_LOAD_RECORD}{MOUTH_LOAD_RECORD}");
        let error = LipDescription::from_bytes(fixture(&duplicate).as_bytes()).unwrap_err();
        assert!(error.contains("may only appear once"), "{error}");
    }

    // The flanged-equivalent load is wired into the nonlinear bore (bounds fixed before the first run, stamp W1): with the row the
    // lip speaks on a different bore end, so the radiated sound differs from the baffled piston's by more than roundoff and less
    // than its own size, and the bore says it convolves a kernel of at least two taps.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn the_flanged_equivalent_mouth_load_reaches_the_nonlinear_bore() {
        use fs_material::gas_species::input::ThermochemicalPhase;
        const STEADY: &str = "frankensim-gesture-schedule-v1\n\
            control_rate_hz\t700\n\
            tracks\t1\n\
            track\tblow\tblowing-pressure\n\
            initial\tpressure-pa\t0e0\n\
            events\t1\n\
            event\t0e0\t1e-2\tpressure-pa\t3e3\n";
        let packs = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
        let bindings = [packs.binding()];
        let sources = ThermochemicalSources {
            species: &bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        };
        let render = |load_expected: &str| {
            let load = load_expected;
            let text = fixture(&format!("{NONLINEAR_RECORD}{load}"))
                .replace("moist_air 293.15 101325 0", SPECIMEN)
                .replace(SCHEDULE, STEADY);
            let description = LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources)).unwrap();
            let mut performance = description.build().unwrap();
            let LipBoreInfo::Nonlinear(info) = &performance.info().bore else {
                panic!("nonlinear fixture");
            };
            let taps = info.radiation_taps;
            if let Some(load) = &info.mouth_load {
                println!(
                    "W1: kernel {} taps, mouth radius {:.6e} m, z_inf {:.6e}, dropped tail {:.2e}, dropped sum {:.2e}, evaluated to {:.0} Hz, minimum Re Z / Z_inf {:.3e}",
                    load.taps, load.mouth_radius_m, load.z_inf, load.dropped_tail_energy, load.dropped_sum, load.evaluated_hz, load.minimum_resistance_ratio
                );
                assert!(load.minimum_resistance_ratio >= -1.0e-3);
            } else {
                assert!(load_expected.is_empty(), "the row was declared but the bore reports no kernel");
            }
            performance.record_mouthpiece();
            let mut radiated = vec![0.0; 2400];
            for block in radiated.chunks_mut(600) {
                performance.block(block).unwrap();
            }
            (radiated, taps)
        };
        let (piston, piston_taps) = render("");
        let (modal, modal_taps) = render(MOUTH_LOAD_RECORD);
        fn rms(x: &[f64]) -> f64 {
            (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt()
        }
        let difference: Vec<f64> = piston.iter().zip(&modal).map(|(a, b)| a - b).collect();
        let relative = rms(&difference) / rms(&piston);
        println!(
            "W1: piston {piston_taps} taps rms {:.4e}, modal kernel {modal_taps} taps rms {:.4e}, relative rms difference {relative:.3e}",
            rms(&piston),
            rms(&modal)
        );
        assert!(piston.iter().chain(&modal).all(|p| p.is_finite()));
        assert!(rms(&piston) > 0.0 && rms(&modal) > 0.0);
        assert!(modal_taps >= 2);
        assert!(relative > 1.0e-6 && relative < 1.0, "relative difference {relative}");
    }

    // B2 (MB 24783): the nonlinear bore holds its own gas inside and the room outside, through `bore_gas` beside `mouth_load`.
    #[test]
    fn the_breath_bore_differs_from_the_room_in_characteristic_impedance_by_the_card_figures() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        let description = parse(&format!("{MODAL}{BREATH}")).expect("parses");
        let (room, bore) = (description.ambient_gas(), description.bore_gas());
        let ratio = bore.density * bore.sound_speed / (room.density * room.sound_speed);
        let reflection = (1.0 - ratio) / (1.0 + ratio);
        println!(
            "B2: rho c breath {:.2} room {:.2} Pa s/m, ratio {ratio:.5}, reflection {reflection:.5}, sound speed ratio {:.5}",
            bore.density * bore.sound_speed,
            room.density * room.sound_speed,
            bore.sound_speed / room.sound_speed
        );
        assert!((0.985..=0.989).contains(&ratio), "{ratio}");
        assert!((0.005..=0.008).contains(&reflection), "{reflection}");
    }

    #[cfg(feature = "nonlinear-brass")]
    const BORE_GAS_ROW: &str = "bore_gas thermochemical 300.15 101325 mole 1 bore-gas 1\n";

    #[cfg(feature = "nonlinear-brass")]
    const STEADY_3KPA: &str = "frankensim-gesture-schedule-v1\n\
        control_rate_hz\t700\n\
        tracks\t1\n\
        track\tblow\tblowing-pressure\n\
        initial\tpressure-pa\t0e0\n\
        events\t1\n\
        event\t0e0\t1e-2\tpressure-pa\t3e3\n";

    #[cfg(feature = "nonlinear-brass")]
    fn two_gas_description(rows: &str, bore_molar_mass: f64) -> Result<LipDescription, String> {
        use fs_material::gas_species::input::ThermochemicalPhase;
        let room = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
        let bore = gas_test_support::packs("bore-gas", bore_molar_mass, 3.5, 200.0);
        let bindings = [room.binding(), bore.binding()];
        let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
        let text = fixture(&format!("{NONLINEAR_RECORD}{rows}")).replace("moist_air 293.15 101325 0", SPECIMEN).replace(SCHEDULE, STEADY_3KPA);
        LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources))
    }

    #[cfg(feature = "nonlinear-brass")]
    fn short_render(description: &LipDescription) -> Vec<f64> {
        let mut performance = description.build().unwrap();
        performance.record_mouthpiece();
        let mut radiated = vec![0.0; 2400];
        for block in radiated.chunks_mut(600) {
            performance.block(block).unwrap();
        }
        radiated
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn a_bore_gas_equal_to_the_room_changes_no_bit_on_the_nonlinear_bore() {
        let plain = two_gas_description(MOUTH_LOAD_RECORD, 0.029).unwrap();
        let equal = two_gas_description(&format!("{MOUTH_LOAD_RECORD}bore_gas thermochemical 293.15 101325 mole 1 specimen-gas 1\n"), 0.029).unwrap();
        assert!(plain.bore_gas_record().is_none() && equal.bore_gas_record().is_some());
        let (a, b) = (short_render(&plain), short_render(&equal));
        assert!(a.iter().any(|x| *x != 0.0), "the fixture must sound");
        assert_eq!(a.len(), b.len());
        assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()), "equal gases must be bit-identical");
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn a_distinct_bore_gas_reaches_the_nonlinear_bore_and_the_kernel_names_the_room_as_its_exterior() {
        let equal = two_gas_description(MOUTH_LOAD_RECORD, 0.029).unwrap();
        let distinct = two_gas_description(&format!("{MOUTH_LOAD_RECORD}{BORE_GAS_ROW}"), 0.0265).unwrap();
        let (room, bore) = (distinct.ambient_gas(), distinct.bore_gas());
        assert_eq!(room.temperature, 293.15);
        assert_eq!(bore.temperature, 300.15);
        assert!(bore.sound_speed > 1.04 * room.sound_speed, "{} against {}", bore.sound_speed, room.sound_speed);
        let performance = distinct.build().unwrap();
        let LipBoreInfo::Nonlinear(info) = &performance.info().bore else {
            panic!("nonlinear fixture");
        };
        let load = info.mouth_load.as_ref().expect("the kernel");
        let area = core::f64::consts::PI * load.mouth_radius_m * load.mouth_radius_m;
        let rho_c = |gas: &GasState| gas.density * gas.sound_speed;
        println!("B2: kernel z_inf x area {:.4} Pa s/m, bore rho c {:.4}, room rho c {:.4}", load.z_inf * area, rho_c(bore), rho_c(room));
        assert!((load.z_inf * area - rho_c(bore)).abs() <= 1e-12 * rho_c(bore), "the kernel's reference is not the bore gas's");
        assert!((load.z_inf * area - rho_c(room)).abs() > 1e-3 * rho_c(room), "the kernel's reference is the room's");
        let (a, b) = (short_render(&equal), short_render(&distinct));
        let rms = |x: &[f64]| (x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64).sqrt();
        let difference: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x - y).collect();
        let relative = rms(&difference) / rms(&a);
        println!("B2: equal-gas render rms {:.4e}, distinct-gas render rms {:.4e}, relative rms difference {relative:.3e}", rms(&a), rms(&b));
        assert!(a.iter().chain(&b).all(|p| p.is_finite()) && rms(&a) > 0.0 && rms(&b) > 0.0);
        assert!(relative > 1e-3, "the bore gas did not reach the render: {relative}");
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn a_nonlinear_bore_gas_is_refused_by_name_lacking_mouth_load_on_another_pressure_or_as_a_mixture() {
        for (rows, needle) in [
            (BORE_GAS_ROW.to_string(), "bore_gas needs modal_baffle (or, on bore nonlinear, mouth_load)"),
            (format!("{MOUTH_LOAD_RECORD}bore_gas thermochemical 300.15 90000 mole 1 bore-gas 1\n"), "one static pressure"),
            (format!("{MOUTH_LOAD_RECORD}bore_gas mixture 300.15 101325 mole 1 N2 1\n"), "its bore_gas row needs a thermochemical record"),
            (format!("{MOUTH_LOAD_RECORD}{MODAL}{BORE_GAS_ROW}"), "modal_baffle requires a linear bore"),
        ] {
            let error = two_gas_description(&rows, 0.0265).unwrap_err();
            assert!(error.contains(needle), "{needle}: {error}");
        }
    }

    #[cfg(not(feature = "nonlinear-brass"))]
    #[test]
    fn finite_mouth_record_cannot_silently_select_a_linear_engine_without_the_feature() {
        let rows = format!("{NONLINEAR_RECORD}{FINITE_MOUTH_RECORD}");
        let error = LipDescription::from_bytes(fixture(&rows).as_bytes()).unwrap_err();
        assert!(error.contains("nonlinear-brass"), "{error}");
    }

    #[cfg(feature = "nonlinear-brass")]
    fn finite_description(row: &str, pressure: &str) -> LipDescription {
        use fs_material::gas_species::input::ThermochemicalPhase;
        let packs = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
        let bindings = [packs.binding()];
        let sources = ThermochemicalSources {
            species: &bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        };
        let steady = format!(
            "frankensim-gesture-schedule-v1\n\
             control_rate_hz\t700\n\
             tracks\t1\n\
             track\tblow\tblowing-pressure\n\
             initial\tpressure-pa\t{pressure}\n\
             events\t0\n"
        );
        let rows = format!("{NONLINEAR_RECORD}{row}");
        let text = fixture(&rows)
            .replace("moist_air 293.15 101325 0", SPECIMEN)
            .replace("audio 24000 2400", "audio 24000 32")
            .replace(SCHEDULE, &steady);
        LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources)).unwrap()
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn finite_mouth_row_reaches_the_real_lip_and_a_sealed_supply_ignores_feed_gestures() {
        let render = |pressure: &str| {
            let description = finite_description(FINITE_MOUTH_RECORD, pressure);
            assert_eq!(description.mouth_reservoir(), Some(&FiniteMouthRecord {
                volume_m3: 2e-5,
                feed_area_m2: 0.0,
                initial_pressure_pa: 106_325.0,
                initial_temperature_k: 293.15,
            }));
            let mut performance = description.build().unwrap();
            assert!(matches!(performance.info().bore, LipBoreInfo::Nonlinear(_)));
            let before = performance.finite_mouth_state().unwrap();
            performance.apply_controls(0).unwrap();
            assert_eq!(performance.finite_mouth_state().unwrap(), before);
            performance.record_mouthpiece();
            let mut radiated = [0.0; 32];
            performance.block(&mut radiated).unwrap();
            assert!(radiated.iter().all(|p| p.is_finite()));
            let after = performance.finite_mouth_state().unwrap();
            assert!(after.mass_kg < before.mass_kg);
            assert!(after.pressure_pa < before.pressure_pa);
            assert!(after.temperature_k < before.temperature_k);
            let report = performance.finite_mouth_report().unwrap();
            assert_eq!(report.after, after);
            assert_eq!(report.feed_mass_kg, 0.0);
            assert_eq!(report.feed_energy_j, 0.0);
            let mouthpiece = performance.mouthpiece_pressure().to_vec();
            assert!(mouthpiece.iter().all(|p| p.is_finite()));
            assert!(mouthpiece.iter().any(|p| *p > 1.0), "the stored charge must drive the actual bore");
            (after, mouthpiece)
        };
        let (quiet_feed, quiet_trace) = render("0e0");
        let (pressurized_feed, pressurized_trace) = render("5e4");
        assert_eq!(quiet_feed, pressurized_feed, "zero feed area isolates stored gas from the gesture");
        assert!(quiet_trace.iter().zip(pressurized_trace)
            .all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn finite_mouth_feed_gesture_fills_storage_in_the_rendered_nonlinear_loop() {
        let description = finite_description("mouth_reservoir 2e-5 4e-6 101325 293.15\n", "3e3");
        let mut performance = description.build().unwrap();
        let before = performance.finite_mouth_state().unwrap();
        performance.apply_controls(0).unwrap();
        assert_eq!(performance.finite_mouth_state().unwrap(), before);
        performance.record_mouthpiece();
        performance.block(&mut [0.0; 32]).unwrap();
        let after = performance.finite_mouth_state().unwrap();
        let report = performance.finite_mouth_report().unwrap();
        assert!(after.pressure_pa > before.pressure_pa + 10.0);
        assert!(after.mass_kg > before.mass_kg);
        assert!(report.feed_mass_kg > 0.0 && report.feed_energy_j > 0.0);
        assert!(performance.mouthpiece_pressure().iter().any(|p| *p > 1.0));
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn finite_mouth_initial_temperature_keeps_the_original_source_domain() {
        let description = finite_description("mouth_reservoir 2e-5 0 106325 1200\n", "0e0");
        let error = description.build().err().expect("initial gas leaves the source's caloric domain");
        assert!(error.contains("mouth_reservoir"), "{error}");
    }

    // MB 3915: the two lip rows run on bore nonlinear (the parser refused them there before). Thirty-two samples of a 3 kPa blow through the
    // real lip and bore: q = 1 is bit-identical to no row, q = 1.4 and the second mode each reach the mouthpiece pressure, and the rows
    // reach the parsed description.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn lip_area_exponent_and_second_mode_run_on_bore_nonlinear() {
        let render = |rows: &str| {
            let mut performance = finite_description(rows, "3e3").build().unwrap();
            performance.apply_controls(0).unwrap();
            performance.record_mouthpiece();
            performance.block(&mut [0.0; 32]).unwrap();
            performance.mouthpiece_pressure().to_vec()
        };
        let same = |a: &[f64], b: &[f64]| a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits());
        let plain = render("");
        assert!(plain.iter().any(|p| p.abs() > 1.0), "the blow must reach the cup for the comparison to mean anything");
        assert!(same(&plain, &render("lip_area_exponent 1\n")), "q = 1 must be bit-identical to no row on the nonlinear bore");
        let steeper = render("lip_area_exponent 1.4\n");
        assert!(steeper.iter().all(|p| p.is_finite()));
        assert!(!same(&plain, &steeper), "q = 1.4 must reach the mouthpiece pressure on the nonlinear bore");
        let two = render("lip_response_mode2 244.5 11.0 -0.073\n");
        assert!(two.iter().all(|p| p.is_finite()));
        assert!(!same(&plain, &two), "the second mode must reach the mouthpiece pressure on the nonlinear bore");
        assert!(same(&two, &render("lip_response_mode2 244.5 11.0 -0.073\n")), "the same card repeats bit for bit");
        let both = finite_description("lip_response_mode2 244.5 11.0 -0.073\nlip_area_exponent 1.4\n", "3e3");
        assert_eq!(both.lip_area_exponent(), Some(1.4));
        assert!(both.lip_mode2().is_some());
        assert!(both.build().is_ok());
    }

    #[test]
    fn bore_wall_parses_and_is_refused_by_name() {
        let parse = |extra: &str| LipDescription::from_bytes(fixture(extra).as_bytes());
        assert_eq!(parse("").expect("parses").bore_wall(), None);
        assert_eq!(parse("bore_wall rigid\n").expect("parses").bore_wall(), None);
        for (extra, needle) in [
            ("bore_wall\n", "missing field"),
            ("bore_wall soft\n", "expected `rigid`, `thin-brass"),
            ("bore_wall thin-brass 0.0005\n", "missing field"),
            ("bore_wall thin-brass 0 100\nbore_wall_source x\n", "expected a positive number"),
            ("bore_wall thin-brass 0.0005 0\nbore_wall_source x\n", "expected a positive number"),
            ("bore_wall rigid\nbore_wall rigid\n", "may only appear once"),
            ("bore_wall_source x\n", "bore_wall_source needs a bore_wall row"),
            ("bore_wall rigid\nbore_wall_source x\n", "bore_wall rigid carries no source"),
            ("bore_wall thin-brass 0.0005 100\n", "needs a bore_wall_source row"),
            ("bore_wall thin-brass 0.0005 100\nbore_wall_source x\n", "nonlinear bore only"),
            ("bore_wall shell-brass 0.0005\n", "missing field"),
            ("bore_wall shell-brass 0.0005 1e-4\n", "missing field"),
            ("bore_wall shell-brass 0.0005 1e-4 4000\n", "missing field"),
            ("bore_wall shell-brass 0.0005 1e-4 4000 15\nbore_wall_source x\n", "even azimuth count"),
            ("bore_wall shell-brass 0.0005 1e-4 4000 12\nbore_wall_source x\n", "even azimuth count"),
            ("bore_wall shell-brass 0.0005 1e-4 4000 256\nbore_wall_source x\n", "count outside its admitted range"),
            ("bore_wall shell-brass 0 1e-4 4000 16\nbore_wall_source x\n", "expected a positive number"),
            ("bore_wall shell-brass 0.0005 1e-4 4000 16\n", "needs a bore_wall_source row"),
            ("bore_wall shell-brass 0.0005 1e-4 4000 16\nbore_wall_source x\n", "nonlinear bore only"),
            ("bore_wall thin-brass 0.0005 100\nbore_wall shell-brass 0.0005 1e-4 4000 16\nbore_wall_source x\n", "may only appear once"),
            ("bore_wall shell-brass 0.0005 1e-4 4000 16\nbore_wall rigid\n", "may only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    // MB 23820 (item 4): `bore_wall rigid` is the row absent bit for bit on bore nonlinear, and the brass wall reaches the
    // mouthpiece pressure by a small, finite amount (the static compliance of a thin brass wall is about 1e-4 of the air's).
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn bore_wall_rigid_is_the_row_absent_bit_for_bit_and_thin_brass_reaches_the_mouthpiece_pressure() {
        let render = |rows: &str| {
            let mut performance = finite_description(rows, "3e3").build().unwrap();
            performance.apply_controls(0).unwrap();
            performance.record_mouthpiece();
            performance.block(&mut [0.0; 32]).unwrap();
            performance.mouthpiece_pressure().to_vec()
        };
        let same = |a: &[f64], b: &[f64]| a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits());
        let plain = render("");
        assert!(plain.iter().any(|p| p.abs() > 1.0), "the blow must reach the cup for the comparison to mean anything");
        assert!(same(&plain, &render("bore_wall rigid\n")), "bore_wall rigid must be bit-identical to no row");
        let rows = "bore_wall thin-brass 0.0005 100\nbore_wall_source test fixture, not a measurement\n";
        let walled = render(rows);
        assert!(walled.iter().all(|p| p.is_finite()));
        assert!(!same(&plain, &walled), "the brass wall must reach the mouthpiece pressure");
        assert!(same(&walled, &render(rows)), "the same card repeats bit for bit");
        let scale = plain.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        let gap = plain.iter().zip(&walled).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(gap <= 0.05 * scale, "{gap} of {scale}");
        let description = finite_description(rows, "3e3");
        let (record, source) = description.bore_wall().expect("declared");
        assert_eq!(*record, BoreWallRecord { thickness_m: 0.0005, quality: 100.0 });
        assert_eq!(source, "test fixture, not a measurement");
        let performance = description.build().unwrap();
        let LipBoreInfo::Nonlinear(info) = &performance.info().bore else {
            panic!("a nonlinear bore");
        };
        let wall = info.bore_wall.expect("the wall's receipt");
        // The 6 mm lead pipe (thickness / radius 0.083) takes Lame's cylinder, the flared bore the membrane.
        assert!(wall.thick_cells > 0 && wall.membrane_cells > 0, "{wall:?}");
        assert_eq!(wall.thick_cells + wall.membrane_cells, wall.cells);
        assert!(wall.largest_thick_half_angle_rad < wall.largest_half_angle_rad);
        assert!(performance.info().bore_wall.is_some());
    }

    // The nonlinear fixture with 10 mm gas cells (the 37 cells of the bell stay cheap for the bell's shell solve) and `samples`
    // output samples of a steady blow.
    #[cfg(feature = "nonlinear-brass")]
    fn coarse_description(row: &str, pressure: &str, samples: u32) -> LipDescription {
        use fs_material::gas_species::input::ThermochemicalPhase;
        let packs = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
        let bindings = [packs.binding()];
        let sources = ThermochemicalSources {
            species: &bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        };
        let steady = format!(
            "frankensim-gesture-schedule-v1\n\
             control_rate_hz\t700\n\
             tracks\t1\n\
             track\tblow\tblowing-pressure\n\
             initial\tpressure-pa\t{pressure}\n\
             events\t0\n"
        );
        let rows = format!("bore nonlinear 0.01 16 1 64 20 8000 1e-2\n{row}");
        let text = fixture(&rows)
            .replace("moist_air 293.15 101325 0", SPECIMEN)
            .replace("audio 24000 2400", &format!("audio 24000 {samples}"))
            .replace(SCHEDULE, &steady);
        LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources)).unwrap()
    }

    // MB 24332 (the bending branch): the bell of the fixture bore (the two widening cones, 0.37 m) as a brass shell, throat clamped.
    // Bounds fixed before the first run (satwall/STAMP-bell-tests.txt): the shell reaches the mouthpiece pressure (a nonzero, finite
    // difference of at most 5 % of the peak); the wall's effect changes by under 25 % when the cut-off doubles (3 kHz to 6 kHz);
    // the same card repeats bit for bit; the receipt names the coupling modes inside the cut-off.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn bore_wall_shell_brass_reaches_the_mouthpiece_pressure_and_converges_in_its_cut_off() {
        let samples = 600;
        let render = |rows: &str| {
            let mut performance = coarse_description(rows, "3e3", samples).build().unwrap();
            performance.apply_controls(0).unwrap();
            performance.record_mouthpiece();
            performance.block(&mut [0.0; 600]).unwrap();
            performance.mouthpiece_pressure().to_vec()
        };
        let source = "bore_wall_source test fixture, not a measurement\n";
        let card = |cut: &str, azimuths: usize| format!("bore_wall shell-brass 4.826e-4 1e-4 {cut} {azimuths}\n{source}");
        let plain = render("");
        let low = render(&card("3000", 16));
        let high = render(&card("6000", 16));
        let fine = render(&card("3000", 32));
        assert!(plain.iter().chain(&low).chain(&high).chain(&fine).all(|p| p.is_finite()));
        let peak = plain.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        assert!(peak > 1.0, "the blow must reach the cup");
        let gap = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
        let (effect, doubling, mesh) = (gap(&plain, &low), gap(&low, &high), gap(&low, &fine));
        println!("bell shell: effect of the wall {effect:.3e} Pa of a {peak:.3e} Pa peak; doubling the cut-off moves it by {doubling:.3e} Pa, doubling the azimuth count by {mesh:.3e} Pa");
        assert!(effect > 0.0, "the shell must reach the mouthpiece pressure");
        assert!(effect <= 0.05 * peak, "{effect} of {peak}");
        assert!(doubling <= 0.25 * effect, "the cut-off is not converged: {doubling} against {effect}");
        assert!(mesh <= 0.25 * effect, "the azimuth count is not converged: {mesh} against {effect}");
        assert!(render(&card("3000", 16)).iter().zip(&low).all(|(a, b)| a.to_bits() == b.to_bits()), "the same card repeats bit for bit");
        let description = coarse_description(&card("3000", 16), "3e3", 32);
        let (record, text) = description.bell_shell().expect("declared");
        assert_eq!(*record, BellShellRecord { thickness_m: 4.826e-4, loss_factor: 1e-4, max_frequency_hz: 3000.0, azimuths: 16 });
        assert_eq!(text, "test fixture, not a measurement");
        let performance = description.build().unwrap();
        let LipBoreInfo::Nonlinear(info) = &performance.info().bore else {
            panic!("a nonlinear bore");
        };
        let bell = info.bell_shell.expect("the shell's receipt");
        assert!(bell.coupling_modes >= 1 && bell.cells_covered >= 20, "{bell:?}");
        assert!(bell.lowest_mode_hz > 50.0 && bell.highest_mode_hz <= 3000.0, "{bell:?}");
        assert!(performance.info().bell_shell.is_some() && info.bore_wall.is_none());
    }

    // The walls' energy account on the card bell (MB 24332: "energy account closes"). Bounds fixed before the run: gas work equals
    // stored + dissipated + room work within 1e-9 of the largest of the four (the account of the bore's own step); the wall holds
    // energy and has dissipated some.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn bore_wall_shell_brass_account_closes_on_the_card_bell() {
        let source = "bore_wall_source test fixture, not a measurement\n";
        let rows = format!("bore_wall shell-brass 4.826e-4 1e-4 3000 16\n{source}");
        let mut performance = coarse_description(&rows, "3e3", 600).build().unwrap();
        performance.apply_controls(0).unwrap();
        performance.record_mouthpiece();
        performance.block(&mut [0.0; 600]).unwrap();
        let ledger = performance.wall_ledger().expect("the shell's account");
        println!("bell shell account {ledger:?}: closing error {:.3e}", ledger.closing_error());
        assert!(ledger.closing_error() <= 1.0e-9, "{ledger:?}");
        assert!(ledger.stored_j > 0.0 && ledger.dissipated_j > 0.0 && ledger.gas_work_j != 0.0, "{ledger:?}");
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn lip_rim_stop_is_refused_with_bore_nonlinear_by_name() {
        let error = finite_description("lip_rim_stop inelastic\n", "3e3").build().err().expect("refused");
        assert!(format!("{error:?}").contains("lip_rim_stop"), "{error:?}");
    }

    #[cfg(feature = "nonlinear-brass")]
    const SPECIMEN: &str = "thermochemical 293.15 101325 mole 1 specimen-gas 1";

    #[cfg(not(feature = "nonlinear-brass"))]
    #[test]
    fn a_nonlinear_bore_is_refused_by_name_without_its_feature() {
        let error = LipDescription::from_bytes(fixture(NONLINEAR_RECORD).as_bytes())
            .expect_err("no nonlinear bore in this build");
        assert!(error.contains("nonlinear-brass"), "{error}");
    }

    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn a_nonlinear_bore_evolves_only_a_source_backed_ambient_gas_through_a_flanged_bell() {
        let error = LipDescription::from_bytes(fixture(NONLINEAR_RECORD).as_bytes())
            .expect_err("a snapshot gas cannot evolve");
        assert!(error.contains("ambient thermochemical"), "{error}");
        let unflanged =
            fixture(NONLINEAR_RECORD).replace("termination flanged", "termination unflanged");
        let error = LipDescription::from_bytes(unflanged.as_bytes())
            .expect_err("the nonlinear bell is a baffled piston");
        assert!(error.contains("termination flanged"), "{error}");
    }

    /// One NASA-9 species: each region is (lowest K, highest K, a0..=a8).
    #[cfg(feature = "nonlinear-brass")]
    fn nasa9_species(
        name: &str,
        molar_mass: f64,
        reference_pressure: f64,
        elemental_reference: &str,
        regions: &[(f64, f64, [f64; 9])],
        source: &str,
    ) -> fs_thermochem::Nasa9StandardStateModelV1 {
        use fs_evidence::ValidityDomain;
        use fs_matdb::{
            ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance,
        };
        use fs_qty::{MolarMass, Pressure};
        use fs_thermochem::{
            ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
            NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
            ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
            required_nasa9_parameter_dims,
        };
        let regions = regions
            .iter()
            .map(|&(low, high, a)| {
                let mut parameters = std::collections::BTreeMap::new();
                for ((name, dims), value) in
                    required_nasa9_parameter_dims().into_iter().take(9).zip(a)
                {
                    parameters.insert(name.to_owned(), LawParameter { value, dims });
                }
                parameters.insert(
                    "reference_pressure".into(),
                    LawParameter {
                        value: reference_pressure,
                        dims: Pressure::DIMS,
                    },
                );
                Nasa9RegionV1::from_card(ConstitutiveModelCard {
                    law: LawId(NASA9_LAW_ID_V1.into()),
                    law_version: NASA9_LAW_VERSION_V1,
                    parameters,
                    state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                    initial_state: InitialStatePolicy::ZeroInternalState,
                    validity: ValidityDomain::unconstrained().with("T", low, high),
                    sources: Vec::new(),
                    provenance: Provenance {
                        source: source.into(),
                        license: "test fixture".into(),
                        artifact: None,
                    },
                })
                .unwrap()
            })
            .collect();
        Nasa9StandardStateModelV1::new(
            SpeciesId::new(name).unwrap(),
            MolarMass::new(molar_mass),
            StandardStateConventionV1::new(
                StandardStatePhaseV1::Gas,
                ReferenceEquationOfStateV1::IdealGas,
                Pressure::new(reference_pressure),
                ElementalReferenceIdV1::new(elemental_reference).unwrap(),
            )
            .unwrap(),
            regions,
        )
        .unwrap()
    }

    /// A data/matdb/seed-v1 species: molar mass, reference pressure and
    /// elemental reference from species.tsv; every region's bounds and
    /// coefficients from the McBride nasa9.tsv.
    #[cfg(feature = "nonlinear-brass")]
    fn seed_species(species: &str, nasa9: &str) -> fs_thermochem::Nasa9StandardStateModelV1 {
        let row: Vec<&str> = species
            .lines()
            .find(|line| line.starts_with("species\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(
            (row[3], row[4], row[5], row[7]),
            ("g/mol", "gas", "ideal-gas", "kPa"),
            "{row:?}"
        );
        let reference_pressure = row[6].parse::<f64>().unwrap() * 1e3;
        let mut regions: Vec<(String, f64, f64, [f64; 9])> = Vec::new();
        for line in nasa9.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            match fields.first() {
                Some(&"region") => {
                    assert_eq!((fields[5], fields[7]), ("K", "Pa"), "{line}");
                    assert_eq!(fields[6].parse::<f64>().unwrap(), reference_pressure);
                    regions.push((
                        fields[2].to_owned(),
                        fields[3].parse().unwrap(),
                        fields[4].parse().unwrap(),
                        [f64::NAN; 9],
                    ));
                }
                Some(&"coefficient") => {
                    let region = regions.iter_mut().find(|r| r.0 == fields[2]).unwrap();
                    let index: usize = fields[3].strip_prefix('a').unwrap().parse().unwrap();
                    region.3[index] = fields[4].parse().unwrap();
                }
                _ => {}
            }
        }
        let regions: Vec<_> = regions
            .into_iter()
            .map(|(_, low, high, a)| {
                assert!(a.iter().all(|c| c.is_finite()), "{a:?}");
                (low, high, a)
            })
            .collect();
        nasa9_species(
            row[1],
            row[2].parse::<f64>().unwrap() * 1e-3,
            reference_pressure,
            row[8],
            &regions,
            "McBride, Zehe and Gordon (2002), NASA/TP-2002-211556 (data/matdb/seed-v1)",
        )
    }

    // A species split at 1000 K, like every McBride (2002) air pack: the field
    // keeps to the room's region, which the source admits; the old fixed
    // 200-1000 K domain crossed the shared bound and refused.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn the_nonlinear_field_keeps_to_the_nasa9_region_that_holds_the_room() {
        use fs_duct::nonlinear::{CaloricDomain, FrozenGas};
        use fs_qty::Temperature;
        use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
        use fs_thermochem::{Composition, CompositionBasis};
        let mut a = [0.0; 9];
        a[2] = 3.5;
        let source = FrozenIdealGasMixtureModelV1::new(
            vec![nasa9_species(
                "split-gas",
                0.029,
                1.0e5,
                "split-fixture",
                &[(200.0, 1000.0, a), (1000.0, 6000.0, a)],
                "two-region test fixture, not data",
            )],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap();
        let crossing = FrozenGas::new(
            source.clone(),
            CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
            || false,
        )
        .expect_err("1000 K belongs to the upper region");
        assert!(
            format!("{crossing}").contains("caloric continuity"),
            "{crossing}"
        );
        let [low, high] = room_region(&source, 293.15).unwrap();
        assert_eq!(low, 200.0);
        assert!(high < 1000.0 && high == 1000.0_f64.next_down(), "{high}");
        FrozenGas::new(
            source.clone(),
            CaloricDomain::new(Temperature::new(low), Temperature::new(high)),
            || false,
        )
        .unwrap();
        // A hot room is held in the upper region; a room off the source refuses.
        assert_eq!(room_region(&source, 1200.0).unwrap(), [1000.0, 6000.0]);
        let error = room_region(&source, 150.0).unwrap_err();
        assert!(error.contains("outside every NASA-9 region"), "{error}");
    }

    // The air the trumpet A/B renders breathe: the McBride (2002) N2, O2, Ar
    // and CO2 tables of data/matdb/seed-v1, each split at 1000 K with no
    // stated caloric continuity. The old fixed 200-1000 K domain refused this
    // air by name; the room's region admits it, a room state round-trips at
    // the ideal-gas sound speed of dry air, and a state above the region
    // refuses rather than extrapolating. The sound-speed bound (0.5 m/s about
    // Cramer's 343.2 m/s for dry air at 20 C, JASA 93 (1993) 2510) was fixed
    // before the first run.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn the_nonlinear_field_admits_the_mcbride_air_in_the_rooms_region() {
        use fs_duct::nonlinear::{CaloricDomain, FrozenGas};
        use fs_qty::Temperature;
        use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
        use fs_thermochem::{Composition, CompositionBasis};
        macro_rules! seed {
            ($name:literal) => {
                seed_species(
                    include_str!(concat!(
                        env!("CARGO_MANIFEST_DIR"),
                        "/../../data/matdb/seed-v1/",
                        $name,
                        "/species.tsv"
                    )),
                    include_str!(concat!(
                        env!("CARGO_MANIFEST_DIR"),
                        "/../../data/matdb/seed-v1/",
                        $name,
                        "-nasa9-mcbride2002/nasa9.tsv"
                    )),
                )
            };
        }
        let air = vec![
            seed!("nitrogen"),
            seed!("oxygen"),
            seed!("argon"),
            seed!("carbon-dioxide"),
        ];
        for species in &air {
            let regions = species.regions();
            assert_eq!(regions.len(), 2, "{:?}", species.species());
            assert_eq!(regions[0].temperature_max().value(), 1000.0);
            assert_eq!(regions[1].temperature_min().value(), 1000.0);
            assert!(species.caloric_reference().is_none());
        }
        // Dry air by mole: N2 78.084, O2 20.946, Ar 0.934, CO2 0.036 %.
        let source = FrozenIdealGasMixtureModelV1::new(
            air,
            Composition::new(
                CompositionBasis::MoleFraction,
                vec![0.78084, 0.20946, 0.00934, 0.00036],
            )
            .unwrap(),
        )
        .unwrap();
        let crossing = FrozenGas::new(
            source.clone(),
            CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
            || false,
        )
        .expect_err("the McBride air states no continuity at 1000 K");
        assert!(
            format!("{crossing}").contains("caloric continuity"),
            "{crossing}"
        );
        let [low, high] = room_region(&source, 293.15).unwrap();
        assert_eq!([low, high], [200.0, 1000.0_f64.next_down()]);
        let gas = FrozenGas::new(
            source,
            CaloricDomain::new(Temperature::new(low), Temperature::new(high)),
            || false,
        )
        .unwrap();
        let r = gas.gas_constant();
        let room = gas
            .conserved(101_325.0 / (r * 293.15), 0.0, 101_325.0)
            .unwrap();
        let state = gas.primitive(room, || false).unwrap();
        println!(
            "mcbride air: R {r:.4} J/(kg K), T {:.9} K, c {:.4} m/s",
            state.temperature, state.sound_speed
        );
        assert!((state.temperature - 293.15).abs() <= 1e-9, "{state:?}");
        assert!((state.sound_speed - 343.2).abs() <= 0.5, "{state:?}");
        let hot = gas
            .conserved(101_325.0 / (r * 1200.0), 0.0, 101_325.0)
            .expect_err("1200 K lies above the room's region");
        assert!(
            format!("{hot}").contains("outside declared interval"),
            "{hot}"
        );
    }

    // The same lip, bore, gas and steady blow through both bore engines: the
    // nonlinear bore must speak in the regime the linear loop finds. The
    // tolerance (3 % on the playing frequency) was fixed before the first run.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn the_lip_speaks_on_the_nonlinear_bore_in_the_linear_loops_regime() {
        use fs_material::gas_species::input::ThermochemicalPhase;
        const STEADY: &str = "frankensim-gesture-schedule-v1\n\
            control_rate_hz\t700\n\
            tracks\t1\n\
            track\tblow\tblowing-pressure\n\
            initial\tpressure-pa\t0e0\n\
            events\t1\n\
            event\t0e0\t1e-2\tpressure-pa\t3e3\n";
        let packs = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
        let bindings = [packs.binding()];
        let sources = ThermochemicalSources {
            species: &bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        };
        let text = |extra: &str| {
            fixture(extra)
                .replace("moist_air 293.15 101325 0", SPECIMEN)
                .replace("audio 24000 2400", "audio 24000 4800")
                .replace(SCHEDULE, STEADY)
        };
        let render = |text: String| {
            let description =
                LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources))
                    .unwrap();
            let mut performance = description.build().unwrap();
            performance.record_mouthpiece();
            let started = std::time::Instant::now();
            let mut radiated = vec![0.0; 4800];
            for block in radiated.chunks_mut(600) {
                performance.block(block).unwrap();
            }
            let seconds = started.elapsed().as_secs_f64();
            // The last 0.1 s: AC level and the playing frequency from the
            // interpolated first and last rising mean crossings.
            let tail = &performance.mouthpiece_pressure()[2400..];
            let mean = tail.iter().sum::<f64>() / tail.len() as f64;
            let rms = (tail.iter().map(|p| (p - mean) * (p - mean)).sum::<f64>()
                / tail.len() as f64)
                .sqrt();
            let rising: Vec<f64> = tail
                .windows(2)
                .enumerate()
                .filter(|(_, w)| w[0] < mean && w[1] >= mean)
                .map(|(i, w)| i as f64 + (mean - w[0]) / (w[1] - w[0]))
                .collect();
            let frequency = match rising.as_slice() {
                [first, .., last] => (rising.len() - 1) as f64 * 24_000.0 / (last - first),
                _ => 0.0,
            };
            let level = radiated.iter().map(|p| p * p).sum::<f64>().sqrt();
            assert!(radiated.iter().all(|p| p.is_finite()) && level > 0.0);
            (performance.info().bore.clone(), rms, frequency, seconds)
        };
        let (_, rms_linear, f_linear, s_linear) = render(text(""));
        let (bore, rms_nonlinear, f_nonlinear, s_nonlinear) = render(text(NONLINEAR_RECORD));
        let LipBoreInfo::Nonlinear(receipt) = bore else {
            panic!("the record realizes the nonlinear bore");
        };
        println!(
            "lip regime: linear {f_linear:.2} Hz rms {rms_linear:.1} Pa ({s_linear:.1} s) | nonlinear {f_nonlinear:.2} Hz rms {rms_nonlinear:.1} Pa ({s_nonlinear:.1} s) | {receipt:?}"
        );
        assert!(
            rms_linear > 150.0,
            "the linear loop must speak: {rms_linear} Pa"
        );
        assert!(
            rms_nonlinear > 150.0,
            "the nonlinear loop must speak: {rms_nonlinear} Pa"
        );
        assert!((f_nonlinear / f_linear - 1.0).abs() <= 0.03);
    }

    // `--mouthpiece-wav` on the nonlinear route (MB 16445, 16793). `advance` appends the block's mouthpiece buffer to
    // the tap after the engine match, so both engines fill it from the very buffer the window diagnostics read.
    // Gauntlet G3 (metamorphic): the tap is a pure observer, so a tapped render is bit-identical to an untapped one;
    // it is empty unless asked, holds one entry per rendered sample, and reproduces the engine's own window statistics.
    #[cfg(feature = "nonlinear-brass")]
    #[test]
    fn the_mouthpiece_tap_on_the_nonlinear_bore_is_the_diagnostics_buffer_and_a_pure_observer() {
        use fs_material::gas_species::input::ThermochemicalPhase;
        const STEADY: &str = "frankensim-gesture-schedule-v1\n\
            control_rate_hz\t700\n\
            tracks\t1\n\
            track\tblow\tblowing-pressure\n\
            initial\tpressure-pa\t0e0\n\
            events\t1\n\
            event\t0e0\t1e-2\tpressure-pa\t3e3\n";
        let packs = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
        let bindings = [packs.binding()];
        let sources = ThermochemicalSources {
            species: &bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        };
        let render = |tap: bool| {
            let text = fixture(NONLINEAR_RECORD)
                .replace("moist_air 293.15 101325 0", SPECIMEN)
                .replace("audio 24000 2400", "audio 24000 4800")
                .replace(SCHEDULE, STEADY);
            let description =
                LipDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources))
                    .unwrap();
            let mut performance = description.build().unwrap();
            assert!(matches!(performance.info().bore, LipBoreInfo::Nonlinear(_)));
            if tap {
                performance.record_mouthpiece();
            }
            let mut radiated = vec![0.0; 4800];
            for block in radiated.chunks_mut(600) {
                performance.block(block).unwrap();
            }
            (
                radiated,
                performance.mouthpiece_pressure().to_vec(),
                performance.windows(),
            )
        };
        let (plain, untapped, _) = render(false);
        let (tapped_audio, tap, windows) = render(true);
        assert!(untapped.is_empty(), "no tap unless asked");
        assert_eq!(tap.len(), 4800, "one tap entry per rendered sample");
        assert!(tap.iter().all(|p| p.is_finite()));
        assert!(
            plain
                .iter()
                .zip(&tapped_audio)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "the tap must change nothing else"
        );
        assert_eq!(windows.len(), 2, "4800 samples at 24 kHz are two tenth-second windows");
        for (k, window) in windows.iter().enumerate() {
            let slice = &tap[k * 2400..(k + 1) * 2400];
            let mean = slice.iter().sum::<f64>() / slice.len() as f64;
            let peak = 0.5
                * (slice.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                    - slice.iter().copied().fold(f64::INFINITY, f64::min));
            assert!(
                (window.bore_mean_pa - mean).abs() <= 1e-9 * mean.abs().max(1.0),
                "window {k}: diagnostics mean {} against the tap's {mean}",
                window.bore_mean_pa
            );
            assert!(
                (window.bore_ac_peak_pa - peak).abs() <= 1e-9 * peak.max(1.0),
                "window {k}: diagnostics peak {} against the tap's {peak}",
                window.bore_ac_peak_pa
            );
        }
        assert!(
            windows[1].bore_ac_peak_pa > 150.0,
            "the lip speaks on the nonlinear bore, so the tap holds a real signal: {} Pa",
            windows[1].bore_ac_peak_pa
        );
    }
}
