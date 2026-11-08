//! A glottal source driving a tract area function, observed at a microphone.
//!
//! The file declares physical inputs only: the lumped fold card, an optional
//! tract wall, the mouth termination, the area function as ordered sections,
//! the gas, the subglottal pressure the player applies, and where the
//! microphone stands. Nothing here names a vowel, a register or a singer.
//!
//! Composition, all existing owners: [`GlottalIsland`] (fold pHS, Bernoulli
//! aperture, fold collision) against the tract's characteristic line; the
//! mouth volume flow of the SAME transfer-matrix duct
//! ([`crate::driving_point::outlet`], driven by the island's incident wave);
//! and the baffled circular-opening receiver, whose pressure is the retarded
//! derivative of that flow. No source waveform or formant filter exists, the one
//! gain is the opt-in `pulse_shimmer` below, and the supraglottal pressure is never
//! the output.
//!
//! Limits, not hidden: infinite-baffle radiation (no head diffraction) unless the
//! opt-in `head_sphere` row puts the mouth on a rigid sphere, a plane-wave series
//! tract (no side branches or nasal coupling), a lumped two-mass, body-cover or
//! one-mass fold, a rigid or single-pin yielding wall, and an ideal pressure
//! source under the folds unless the opt-in `subglottal_airway` row puts a
//! symmetric tree there.
//!
//! ```text
//! frankensim-voice-performance-v1
//! # comment lines and blank lines are ignored
//! audio <output_rate_hz> <samples> <mechanics_per_output>
//! full_scale_pa <pa>
//! ambient moist_air <temperature_k> <pressure_pa> <relative_humidity>
//! ambient mixture <temperature_k> <pressure_pa> <mole|mass> <count> <species fraction>...
//! ambient thermochemical <temperature_k> <pressure_pa> <mole|mass> <count> <id fraction>...
//! observer_law <none|iso9613|stokes-kirchhoff|relaxation LAW> (optional; exterior propagation)
//! fold <two_mass|body_cover|one_dof> <m_lower_kg> <m_upper_kg> <k_lower_n_m> <k_upper_n_m>
//!      <k_couple_n_m> <damping_ratio> <rest_gap_m> <fold_length_m> <thickness_m>
//!      [<thickness_upper_m>]      (the second viscous section; default 0.05 cm)
//! fold_damping_upper <open_ratio>                    (required on two_mass and body_cover, refused on one_dof)
//! fold_body <m_body_kg> <k_body_n_m> <damping_ratio_body> <rest_gap_upper_m>   (required on body_cover only)
//! tract_wall <surface_density_kg_m2> <stiffness_n_m3> <resistance_pa_s_m>   (optional)
//! termination <flanged|unflanged|ideal_open>
//! tract_source <identifier without spaces>
//! tract_licence <free text to end of line>
//! section <area_m2> <length_m>                       (repeat, glottis first)
//! microphone <x_m> <y_m> <z_m> <radial_rings> <angular_points> <max_frequency_hz>
//! subglottal <time_s> <pressure_pa>                  (repeat, strictly increasing time)
//! subglottal_airway weibel-1963 saab-2021-table-1-B <generations>   (optional; the trachea and bronchial tree below the folds)
//! tension_modulation <rate_hz> <relative_amplitude>  (optional; vibrato)
//! tension_flutter <seed> <centre_hz> <bandwidth_hz> <sd_rel>  (optional; pitch movement, adds to any vibrato)
//! tension_flutter_second <seed> <centre_hz> <bandwidth_hz> <sd_rel>  (optional, needs tension_flutter; a second filter)
//! fold_mapping bilateral                             (optional; the fold line is one fold of a pair)
//! fold_shape <sine|rigid> <uniform|tapered>          (optional, two_mass only; the gap varies along the fold)
//! aspiration_noise <seed> <xi_pa_s3_per_m6> [<re_crit>]   (optional; turbulence noise at the glottis)
//! jet_noise <seed> <gain> <st_lo> <st_hi> <re_on>         (optional; constriction noise of the glottal jet, Verge's fraction)
//! constriction_noise <seed> <duct_radius_m>               (optional; Oldham and Ukpoho's printed constriction noise)
//! pulse_shimmer <seed> <sd_rel>                           (optional; a per-cycle gain on the wave the tract receives)
//! head_sphere <radius_m> <radius source without spaces>   (optional; the mouth as a cap on a rigid sphere)
//! head_load cap                                           (optional, needs head_sphere and termination flanged; the lips see the sphere's load)
//! torso_sphere <radius_m> <neck_gap_m> <source without spaces>   (optional, needs head_sphere; a torso sphere under the head)
//! tract_gas <any `ambient` record>                        (optional, needs head_load; the gas inside the tract)
//! keyframe <name> source <identifier without spaces>      (optional; a further measured shape of a moving tract)
//! keyframe <name> licence <free text to end of line>
//! keyframe_section <name> <area_m2> <length_m>            (repeat, glottis first)
//! tract_move <start_s> <end_s> <name|card> <steps>        (repeat, in time order; the tract moves to that shape)
//! ```
//!
//! A moving tract (pexfp.327.1): each `tract_move` takes the tract from the shape it holds (the card's own `section`
//! rows until the first move) to keyframe `name`, or back to the card's shape with `card`, through `steps` shapes at
//! fractions `k / steps` of [`TractChart::log_area_morph`] (log area interpolated, length linear). Shape `k`'s
//! crossfade starts at `start + (k - 1) (end - start) / steps` and lasts `(end - start) / steps`, so the tract is
//! exactly the target from `end`. Every shape's island line and mouth-flow line are realized at build time by the same
//! ladder, load, wall and gases as the card's own shape, and the render switches them with the outgoing history carried
//! ([`crate::keyframe_line`]). Three limits, not hidden: the glottal section is held constant, so a keyframe whose
//! first-section area is not the card's is refused by name and every intermediate shape takes the card's first area
//! (the island's `zc` and fold law read it); the microphone, the receiver and the head sphere keep the card's mouth
//! radius; and the shapes are frozen linear tracts crossfaded, coefficient interpolation with no volume source from the
//! moving walls a real tract has. A keyframe that no move names is refused. Without the rows the build and the render
//! are bit for bit what they were.
//!
//! `tract_gas` fills the tract with its own gas (breath) while `ambient` stays the exterior (the room). The island's
//! line, the mouth-flow line, the passive formants and the glottal flow's density and viscosity take the tract gas; the
//! lips' load table, the receiver and the head and torso spheres take `ambient`. The table is an absolute impedance
//! (Pa s/m^3) computed in the exterior and checked against it, so nothing is converted at the lips: pressure and volume
//! flow are continuous there. It needs `head_load`, because a declared termination would be evaluated in the tract's gas,
//! and the two gases share one static pressure (mean flow is not modeled). The medium is stationary: the breath does not
//! mix with the room or replace the tract's air over a note. Without the row the one `ambient` gas fills the tract and the
//! exterior exactly as before.
//!
//! `subglottal_airway` puts an airway below the folds ([`crate::upstream_airway`]): the first `<generations>` of
//! Weibel's symmetric tree as Saab, Issa, Zeineddine, Tartakovsky and Lakkis (2021, arXiv:2111.07627) print it in Table
//! 1 column B, generation `n` being `2^(n-1)` cylinders with Bessel losses in the tract's gas, and every tube carrying
//! the card's `tract_wall` pin on its own perimeter (rigid without one). The `subglottal` schedule then becomes the lung
//! pressure at the tree's far end, and every step the folds see the airway's port ([`GlottalIsland::with_upstream`]).
//! The aspiration, jet and constriction sources stay in series at the glottis, and their laws still read the schedule
//! as the blowing pressure. How many generations a card takes, and by which rule, is the card's to state; nothing here
//! picks them. No measured per-generation airway wall is carried, so a page showing a render must say which wall the
//! airway took. Without the row the schedule is the pressure under the folds and a performance is bit-identical to one
//! written before the row existed.
//!
//! `observer_law` binds the existing causal gas law to that exterior record on
//! every mouth-quadrature path. Absorption, dispersion and physical travel
//! reach the microphone before the paths are summed; the receiver reports any
//! common causal padding separately. The law never changes the tract, its load
//! or the fold mechanics. `relaxation LAW` uses the shared material grammar
//! (`iso9613-terrestrial-air` or explicit species data). No law is inferred
//! for arbitrary mixtures. Without this row the original lossless receiver
//! remains; `none` is an explicit lossless comparison. Preparation uses the
//! finite `GasPropagationLimits::surface_reference()` work ceilings. An
//! optional head/torso FIR retains its existing lossless scattering model;
//! this does not solve scattering in a dissipative, dispersive medium.
//!
//! `head_sphere` observes the mouth as a cap of its own area on a rigid sphere of that radius, centred `radius` behind
//! the mouth plane, instead of a disk in an infinite baffle ([`crate::pcm_wav::baffled::head::HeadSphere`]: the
//! baffled receiver's pressure through a causal FIR, `S / B`, blended to 1 above the microphone's resolved band). Only
//! the observed pressure changes: the lips' radiation load is still the `termination` row's, which a page using it must
//! state as a limit, unless `head_load` is also present. The second field names where the radius is printed. The FIR's
//! latency joins [`VoicePerformance::filter_delay_samples`].
//!
//! `head_load cap` puts the same cap's radiation load at the lips: the declared flanged load plus the cap's change on the
//! sphere, `Z_flanged + (conj(Z_cap) - Z_piston)` ([`crate::pcm_wav::baffled::head::head_mouth_load`]), tabulated on the
//! bins of the lines' finest rung with DC 0. The island's line, the mouth-flow line and the passive formants all read
//! that one table. It needs the `head_sphere` row (the radius is that row's) and `termination flanged` (the table adds to
//! that load). The receiver stays the head's `S / B` on the mouth flow, which now sees the sphere's load. Without the row
//! every line takes exactly its old path.
//!
//! `torso_sphere` puts a rigid torso sphere of that radius under the head sphere, the neck gap between their surfaces
//! (Algazi, Duda and Thompson's snowman, [`crate::pcm_wav::baffled::head::TorsoSphere`]; the mouth on the head's
//! horizontal equator facing front and the card's `+y` up are our choices). The torso's change joins the head's `S`
//! below 1.25 times the microphone's band, and with `head_load` the lips' load gains the cap's change from the torso over
//! the same band ([`crate::pcm_wav::baffled::head::head_torso_mouth_load`]), so the lips and the microphone see one
//! geometry. It needs the `head_sphere` row; the third field names where the sizes are printed. Without the row nothing
//! changes.
//!
//! `constriction_noise` is an opt-in glottal source on the printed duct-constriction law of Oldham and Ukpoho (J. Sound Vib.
//! 140(2), 1990; Guedel, AIVC Contributed Report 02, eq. 3 and Fig. 10; [`crate::turbulence_noise::ConstrictionNoise`]): a
//! series pressure source added to the step's drive, in each third octave below the duct cutoff `100 / r` whose Strouhal
//! number `f d / U_c` lies on the printed line, of rms `K(St) P` with `K` the printed trend, `P` the subglottal gesture,
//! `U_c = sqrt(2 P / rho)` and `d = pi r / 2`. Two fields and no gain, so it has nothing to tune: `r` is the duct the
//! constriction opens into, and the card says which radius it took. The law was measured on dampers and orifice plates in a
//! 300 mm duct at 10-25 m/s, so putting it on a glottis is an extrapolation in size and Mach number that any page showing a
//! render must say: "constriction noise from a duct-constriction law (Oldham and Ukpoho 1990, via Guedel AIVC CR02: dampers
//! and orifice plates in a 300 mm duct at 10-25 m/s), applied to the glottal jet as an extrapolation". It is silent while the
//! glottal flow is not forward. It may stand on a card with `jet_noise`, but those are two laws for one jet, so a page
//! carrying both must say so. Without the key nothing is added and a performance is bit-identical to one written before the
//! key existed.
//!
//! `jet_noise` is a second opt-in glottal source, beside `aspiration_noise` and not a replacement for it, and the two may
//! be on one card: a series pressure source of rms `gain * rho U_j^2 / 2` over `St_h = f h / U_j` in `[st_lo, st_hi]`,
//! shut while `Re_h = U_j h / nu` (which is `Q / (L nu)`) is below `re_on`, with `U_j = Q / (L h)` the glottal jet speed
//! from the island's own previous-step flow and gap, bounded by the Bernoulli speed of the driving pressure. Every field
//! is required. The form is Stevens' (1971), the constant, band and gate are Verge's (1995, a flue jet at a labium), and
//! applying them to the glottis is an extrapolation that any page showing a render must say (see
//! [`crate::turbulence_noise::JetNoiseCard`] for the law, the label and the measured shortfall: at this constant the noise
//! it adds is about 20-30 dB too weak to close the flux or hnr gap on the F4 notes). It is added to the step's drive in series with the
//! glottal channel exactly as the aspiration source is. Without the key nothing is added and a performance is
//! bit-identical to one written before the key existed; a gate that never opens adds exactly zero.
//!
//! `aspiration_noise` adds a seeded turbulence pressure source in series with the glottal channel, driven by the
//! island's own flow and aperture (the law, its sources and the four readings it takes of them are in
//! [`crate::turbulence_noise`]). The gain `xi` has no default: no source gives one for glottal aspiration, so a card
//! without it is refused by name. `re_crit` is optional and defaults to the 1700 that Elie and Laprie take from
//! Sondhi and Schroeter 1987 (who call it arbitrary). The source is added to the step's subglottal drive, which for a
//! series source downstream of the folds is the same loop equation (`p_sub + p_n`) and the same fold force
//! (`p_sub + p_n - p_supra`). Without the key nothing is added and a performance is bit-identical to one written
//! before the key existed; with the key but a gate that never opens (flow below the onset), the added pressure is
//! exactly zero and the output is bit-identical as well.
//!
//! `tension_modulation` is a cricothyroid gesture: the fold's three stiffnesses are scaled by
//! `1 + a sin(2 pi r t)` ([`GlottalIsland::set_tension`]), so the emergent f0 swings by about
//! `1731 e a` cents peak (`e` ~ 0.5, f0 ~ sqrt(stiffness)). Without the key the tension is never touched
//! and a performance is bit-identical to one written before the key existed.
//!
//! `tension_flutter` is the same gesture driven by noise instead of a sinusoid: the slow, unrepeating F0 movement of
//! a tone sung as straight as a singer can, in the form Ternstrom and Friberg print (STL-QPSR 30(3):1-14, 1989):
//! white noise through a resonant second-order filter of centre frequency and bandwidth. The stiffnesses scale by
//! `1 + y`, with `y` of stationary standard deviation `sd_rel`, set once per output sample. It is not per-cycle
//! jitter: the filter's bandwidth is a few hertz, so `y` has no steps. The row moves the tension, not the pitch
//! directly: the card, not this code, says how `sd_rel` was derived from a pitch extent, which depends on the note.
//! On a card with `tension_modulation` the two add: the stiffnesses scale by `1 + a sin(2 pi r t) + y`, one tension set
//! once per output sample. Without the key the tension is never touched and a
//! performance is bit-identical to one written before the key existed; `sd_rel` 0 is bit-identical as well.
//!
//! `tension_flutter_second` is a second filter on that same gesture, for a card that models two bands of movement: the
//! same law and the same four fields, its `y` added to the first row's. It draws from tile 1 of the row's key, so it is
//! independent of the first row even on the same seed. A card needs `tension_flutter` to take it. Without the key a
//! performance is bit-identical to one written before the key existed; `sd_rel` 0 is bit-identical as well.
//!
//! `pulse_shimmer` is the one gain in this module, and it is opt-in. At each opening of the glottis (a step that ends
//! open after one that ended shut) it draws `g` from a seeded normal stream of standard deviation `sd_rel`, and until the
//! next opening the island's incident wave reaches the mouth-flow line multiplied by `1 + g`. It stands for the
//! cycle-to-cycle amplitude variation printed for real voices (shimmer), whose cause this code does not model. The
//! island, its own load, the tension and every noise source never see the gain, so the fold's motion and phase are bit
//! for bit those of the card without the key; only the radiated pressure changes. The gain multiplies the whole incident
//! wave of a cycle, which includes the tract's own wave returning from the glottis, so a cycle's gain also lands on the
//! end of the cycle before. For independent normal gains, local shimmer (the mean absolute difference of consecutive
//! cycle amplitudes over their mean) is `2 sd_rel / sqrt(pi)`; the card, not this code, says which print it took and
//! whether that print's tracker floor was taken out. A fold that never shuts draws nothing and is left as it was.
//! Without the key nothing is multiplied and a performance is bit-identical to one written before the key existed;
//! `sd_rel` 0 is bit-identical as well.
//!
//! `fold_mapping bilateral` reads the fold line's masses and stiffnesses as ONE fold of Ishizaka and
//! Flanagan's mirror pair and halves them for the island's full-gap coordinate
//! ([`FoldCard::bilateral`]); without it they enter the gap coordinate as written, and a performance is
//! bit-identical to one written before the key existed.
//!
//! `fold_shape` lets the gap vary along a `two_mass` fold, so closure spreads along it (the law is in
//! [`crate::glottis`]). `sine` holds the ends, `rigid` moves the whole mass; `tapered` narrows the rest gap to nothing
//! at the anterior commissure, `uniform` keeps it. Held ends are printed for prescribed motion only (Birkholz 2019,
//! after Titze); the self-oscillating form and its modal mass are ours, and any page showing a render must say so.
//! `rigid tapered` is the printed self-oscillating fold of Birkholz, Kroger and Neuschaefer-Rube 2011. `rigid uniform`
//! is the plain island reached through the strips. Without the key a performance is bit-identical to one written
//! before the key existed.
//!
//! On a `two_mass` fold the fold line's damping ratio is the lower mass's open ratio and
//! `fold_damping_upper` the upper mass's, and the law adds `1.0` to each while its own section is shut
//! (Ishizaka and Flanagan 1972, eq 3 and 4; [`FoldCard::damping_ratio`]). The key has no default: a
//! `two_mass` card without it is refused by name. It is a key and not a field of the fold line because
//! that line is positional, and a second ratio there would shift the fields after it.
//!
//! `body_cover` is Story and Titze's three-mass fold with its own pressure and flow law (NCVS Status and Progress
//! Report 6, 1994; [`crate::glottis`], Body-cover fold). The fold line's masses, springs, coupling, damping ratio, rest
//! gap, length and thicknesses are the cover's, lower first, with the rest gap the lower section's; `fold_damping_upper`
//! is the upper cover mass's ratio, and `fold_body` gives the body's mass, its spring to the thyroid, its damping ratio
//! (it does not step at closure) and the upper section's rest gap. All are for one fold, so `fold_mapping bilateral`
//! halves the body with the cover. A `fold_shape` is refused on it. Tension and its gestures scale the body's spring
//! with the cover's.

use crate::acoustic_realize::AcousticRealizeError;
use crate::driving_point::LineLadder;
use crate::driving_point::outlet::{mouth_flow_line, mouth_flow_line_exterior};
use crate::glottis::{
    FoldBody, FoldCard, FoldCardSource, FoldMode, FoldRestGap, FoldShape, GLOTTAL_LINE_MAX_TAPS, GlottalIsland,
    glottal_tract_line_taps,
};
use crate::keyframe_line::ShapeLine;
use crate::pcm_wav::baffled::head::{HeadSphere, TorsoSphere, head_mouth_load, head_torso_mouth_load};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use crate::air_path::input::{ObserverLawInput, resolve_observer_law};
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationLimits, GasPropagationReport};
use crate::pcm_wav::baffled::openings::CircularOpening;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::tract::{TractChart, TractSection};
use crate::turbulence_noise::{
    AspirationCard, ConstrictionNoise, ConstrictionNoiseCard, GlottalAspiration, GlottalJetNoise, JetNoiseCard,
    RE_CRIT_SONDHI_SCHROETER,
};
use crate::upstream_airway::{AirwayChain, AirwayReceipt, UpstreamAirway, WEIBEL_SAAB_2021_TABLE_1_B_MM};
use fs_duct::{MouthLoad, Termination};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources, read_gas_record};
use fs_math::det;
use fs_phs::WallPin;
use fs_rand::{Stream, StreamKey};

/// First line of every voice performance file.
pub const VOICE_PERFORMANCE_SCHEMA: &str = "frankensim-voice-performance-v1";
const HEAD_LOAD_NEEDS_SPHERE: &str = "head_load needs the head_sphere row: the load is that sphere's";
const HEAD_LOAD_NEEDS_FLANGED: &str = "head_load needs `termination flanged`: it adds the sphere's change to the flanged load";
const TORSO_NEEDS_SPHERE: &str = "torso_sphere needs the head_sphere row: the torso sits under that head";
const TRACT_GAS_NEEDS_HEAD_LOAD: &str =
    "tract_gas needs the head_load row: the lips' load must be computed in the exterior, not the tract's gas";
/// The `subglottal_airway` row's tree token: Weibel's symmetric model.
pub const SUBGLOTTAL_AIRWAY_TREE: &str = "weibel-1963";
/// The `subglottal_airway` row's geometry token: the dimensions as Saab et al. (2021) print them, Table 1 column B.
pub const SUBGLOTTAL_AIRWAY_GEOMETRY: &str = "saab-2021-table-1-B";
/// Largest accepted file.
pub const MAX_VOICE_PERFORMANCE_BYTES: usize = 1 << 20;
const MAX_SECTIONS: usize = 1024;
const MAX_KEYFRAMES: usize = 64;
const MAX_TRACT_MOVES: usize = 64;
const MAX_MOVE_STEPS: u64 = 256;
/// A moving tract's whole menu of shapes, the card's own included.
const MAX_TRACT_SHAPES: usize = 512;
/// Threads that build a moving tract's further shapes at once. Each shape is built alone by the same serial code, so the
/// lines do not depend on the count (`build_shapes`).
const MAX_SHAPE_WORKERS: usize = 16;
const MAX_SCHEDULE_POINTS: usize = 4096;
const MAX_SECONDS: u64 = 120;
const MAX_BLOCK: usize = 65_536;

/// A parsed, validated, not yet realized performance.
#[derive(Debug, Clone)]
pub struct VoiceDescription {
    /// Output (microphone) clock [Hz].
    pub output_rate_hz: u32,
    /// Output samples to render.
    pub samples: u64,
    /// Integer physical substeps per output sample, 1..=16.
    pub mechanics_per_output: u32,
    /// Pascals mapped to PCM full scale by the command output.
    pub full_scale_pa: f64,
    /// One immutable gas record fills the exterior, and the tract too unless `tract_gas` is present.
    gas: GasRecord,
    /// Opt-in `tract_gas`: the gas inside the tract (module doc); `None` fills the tract with `gas`.
    tract_gas: Option<GasRecord>,
    /// Explicit exterior absorption/dispersion; absent keeps the lossless microphone.
    pub observer_law: Option<ObserverLawInput>,
    /// `true` for the two-mass fold, `false` for the one-mass image.
    pub two_mass: bool,
    /// Lumped fold parameters.
    pub fold: FoldCard,
    /// Opt-in shape of a two-mass fold along its length; `None` is the printed island.
    pub fold_shape: Option<FoldShape>,
    /// The body of a `body_cover` fold, in the island's coordinate; `None` on the other folds.
    pub fold_body: Option<FoldBody>,
    /// Optional single-pin yielding tract wall.
    pub wall: Option<WallPin>,
    /// Mouth radiation load.
    pub termination: Termination,
    /// Area function, glottis first.
    pub sections: Vec<TractSection>,
    /// Identifier of the area-function source.
    pub tract_source: String,
    /// Licence or terms the data were admitted under.
    pub tract_licence: String,
    /// Microphone in the mouth-plane frame (origin at the mouth centre,
    /// +z out of the baffle) and its spatial resolution.
    pub microphone: CircularOutletReceiver,
    /// Subglottal pressure breakpoints `(time_s, pressure_pa)`: the pressure under the folds, or with
    /// `subglottal_airway` the lung pressure at the airway's far end.
    pub subglottal: Vec<(f64, f64)>,
    /// Opt-in `subglottal_airway`: this many generations of Weibel's tree as Saab et al. (2021) print it in Table 1
    /// column B, below the glottis (module doc); `None` applies the schedule under the folds.
    pub subglottal_airway: Option<usize>,
    /// Vibrato as a fold-tension gesture `(rate_hz, relative_amplitude)`: the stiffnesses scale by
    /// `1 + a sin(2 pi r t)`. `None` leaves the tension at the card's value for the whole performance.
    pub tension_modulation: Option<(f64, f64)>,
    /// Slow random movement of the fold tension (filtered noise); `None` leaves the tension at the card's value.
    pub tension_flutter: Option<TensionFlutterCard>,
    /// A second filter on the same gesture, its `y` added to the first's; `None` adds nothing.
    pub tension_flutter_second: Option<TensionFlutterCard>,
    /// Opt-in glottal aspiration noise ([`crate::turbulence_noise`]); `None` adds nothing.
    pub aspiration_noise: Option<AspirationCard>,
    /// Opt-in glottal jet noise (Verge's fraction on the glottal jet, an extrapolation; [`crate::turbulence_noise::JetNoiseCard`]);
    /// `None` adds nothing.
    pub jet_noise: Option<JetNoiseCard>,
    /// Opt-in glottal constriction noise (Oldham and Ukpoho's printed law on the glottal jet, an extrapolation;
    /// [`crate::turbulence_noise::ConstrictionNoise`]); `None` adds nothing.
    pub constriction_noise: Option<ConstrictionNoiseCard>,
    /// Opt-in per-cycle gain on the wave the tract receives ([`PulseShimmerCard`]); `None` multiplies nothing.
    pub pulse_shimmer: Option<PulseShimmerCard>,
    /// Opt-in head: the mouth as a cap on a rigid sphere ([`HeadSphereCard`]); `None` keeps the infinite baffle.
    pub head_sphere: Option<HeadSphereCard>,
    /// Opt-in `head_load cap`: the lips see that sphere's radiation load (module doc); `false` keeps the declared one.
    pub head_load: bool,
    /// Opt-in torso under the head sphere ([`TorsoSphereCard`]); `None` keeps the head alone.
    pub torso_sphere: Option<TorsoSphereCard>,
    /// Further measured shapes of a moving tract, in the order first named (module doc); empty for a fixed tract.
    pub keyframes: Vec<TractKeyframe>,
    /// The tract's movements, in time order; empty for a fixed tract.
    pub tract_moves: Vec<TractMove>,
}

/// A further measured shape of a moving tract: the `keyframe` and `keyframe_section` rows of one name.
#[derive(Debug, Clone, PartialEq)]
pub struct TractKeyframe {
    /// The name the `tract_move` rows use.
    pub name: String,
    /// Identifier of the area-function source.
    pub source: String,
    /// Licence or terms the data were admitted under.
    pub licence: String,
    /// Area function, glottis first; the first area is the card's own (the glottal section is held).
    pub sections: Vec<TractSection>,
}

/// A `tract_move` row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TractMove {
    /// When the tract starts to move [s].
    pub start_s: f64,
    /// When it holds the target [s].
    pub end_s: f64,
    /// The target: an index into [`VoiceDescription::keyframes`], or `None` for the card's own shape.
    pub target: Option<usize>,
    /// Shapes on the way, the target included.
    pub steps: u32,
}

/// A keyframe being read: its rows may come in any order.
#[derive(Default)]
struct KeyframeDraft {
    name: String,
    source: Option<String>,
    licence: Option<String>,
    sections: Vec<TractSection>,
}

fn keyframe_draft<'a>(drafts: &'a mut Vec<KeyframeDraft>, name: &str) -> Result<&'a mut KeyframeDraft, String> {
    if name == "card" {
        return Err("`card` names the card's own shape and cannot be a keyframe".into());
    }
    let at = match drafts.iter().position(|d| d.name == name) {
        Some(at) => at,
        None => {
            if drafts.len() >= MAX_KEYFRAMES {
                return Err("too many keyframes".into());
            }
            drafts.push(KeyframeDraft {
                name: name.to_string(),
                ..KeyframeDraft::default()
            });
            drafts.len() - 1
        }
    };
    Ok(&mut drafts[at])
}

/// Check the keyframes and movements against the card's own tract (module doc).
fn admit_motion(
    card_sections: &[TractSection],
    drafts: Vec<KeyframeDraft>,
    moves: Vec<(f64, f64, String, u64)>,
    duration_s: f64,
) -> Result<(Vec<TractKeyframe>, Vec<TractMove>), String> {
    let mut keyframes = Vec::with_capacity(drafts.len());
    for draft in drafts {
        let name = draft.name;
        let source = draft.source.ok_or_else(|| format!("keyframe {name:?} needs a `source` row"))?;
        let licence = draft.licence.ok_or_else(|| format!("keyframe {name:?} needs a `licence` row"))?;
        if draft.sections.is_empty() {
            return Err(format!("keyframe {name:?} needs `keyframe_section` rows"));
        }
        TractChart::try_new(draft.sections.clone(), &source, &licence)
            .map_err(|e| format!("keyframe {name:?} chart refused: {e:?}"))?;
        if draft.sections[0].area_m2.to_bits() != card_sections[0].area_m2.to_bits() {
            return Err(format!(
                "keyframe {name:?}: the first section's area must be the card's own, the glottal section is held constant"
            ));
        }
        keyframes.push(TractKeyframe {
            name,
            source,
            licence,
            sections: draft.sections,
        });
    }
    let mut admitted = Vec::with_capacity(moves.len());
    let mut used = vec![false; keyframes.len()];
    let mut holding: Option<usize> = None;
    let mut free_from = 0.0;
    for (start_s, end_s, name, steps) in moves {
        let target = if name == "card" {
            None
        } else {
            let at = keyframes
                .iter()
                .position(|k| k.name == name)
                .ok_or_else(|| format!("tract_move names an unknown keyframe {name:?}"))?;
            used[at] = true;
            Some(at)
        };
        if !(start_s >= free_from && end_s > start_s && end_s <= duration_s) {
            return Err("tract_move times must increase, not overlap, and end inside the render".into());
        }
        if !(1..=MAX_MOVE_STEPS).contains(&steps) {
            return Err(format!("tract_move steps must be in 1..={MAX_MOVE_STEPS}"));
        }
        if target == holding {
            return Err("a tract_move must change the shape the tract holds".into());
        }
        admitted.push(TractMove {
            start_s,
            end_s,
            target,
            steps: steps as u32,
        });
        holding = target;
        free_from = end_s;
    }
    if let Some(at) = used.iter().position(|u| !u) {
        return Err(format!("keyframe {:?} is named by no tract_move", keyframes[at].name));
    }
    let shapes: u64 = admitted.iter().map(|m| u64::from(m.steps)).sum();
    if shapes >= MAX_TRACT_SHAPES as u64 {
        return Err(format!("a moving tract may hold at most {MAX_TRACT_SHAPES} shapes"));
    }
    Ok((keyframes, admitted))
}

/// A `head_sphere` row: the sphere's radius and where it is printed.
#[derive(Debug, Clone, PartialEq)]
pub struct HeadSphereCard {
    /// Head radius [m].
    pub radius_m: f64,
    /// Where the radius is printed (an identifier without spaces).
    pub source: String,
}

/// A `torso_sphere` row: the torso sphere under the head and where its sizes are printed.
#[derive(Debug, Clone, PartialEq)]
pub struct TorsoSphereCard {
    /// The torso's radius and the neck gap between the two surfaces [m].
    pub torso: TorsoSphere,
    /// Where the sizes are printed (an identifier without spaces).
    pub source: String,
}

fn field(rest: &[&str], i: usize, what: &str) -> Result<f64, String> {
    let token = rest.get(i).ok_or_else(|| format!("missing {what}"))?;
    let value: f64 = token
        .parse()
        .map_err(|_| format!("{what} is not a number: {token:?}"))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(format!("{what} must be finite"))
    }
}

fn whole(rest: &[&str], i: usize, what: &str) -> Result<u64, String> {
    let token = rest.get(i).ok_or_else(|| format!("missing {what}"))?;
    token
        .parse::<u64>()
        .map_err(|_| format!("{what} is not a non-negative integer: {token:?}"))
}

fn arity(rest: &[&str], n: usize, key: &str) -> Result<(), String> {
    if rest.len() == n {
        Ok(())
    } else {
        Err(format!("`{key}` takes {n} fields, got {}", rest.len()))
    }
}

fn once<T>(slot: &mut Option<T>, value: T, key: &str) -> Result<(), String> {
    if slot.replace(value).is_some() {
        Err(format!("`{key}` may only appear once"))
    } else {
        Ok(())
    }
}

impl VoiceDescription {
    /// Parse and validate a performance file.
    ///
    /// # Errors
    /// Oversize or non-UTF-8 input, an unknown or repeated key, a malformed or
    /// non-finite field, or any non-physical value, each named.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        Self::from_bytes_with_gas_sources(bytes, None)
    }

    /// Admit legacy air, named mixtures or source-backed NASA-9 gas records.
    /// All tract, fold-flow and radiation consumers borrow the same state,
    /// unless a `tract_gas` row fills the tract with its own (module doc);
    /// source receipts are retained without changing its coefficients.
    /// The medium is stationary; this does not model species transport.
    ///
    /// # Errors
    /// Ordinary file refusals plus missing, ambiguous or out-of-domain gas
    /// sources. NASA-9 records never fall back to fixed heat capacities.
    pub fn from_bytes_with_gas_sources(
        bytes: &[u8],
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, String> {
        if bytes.len() > MAX_VOICE_PERFORMANCE_BYTES {
            return Err("voice performance file exceeds its size limit".into());
        }
        let text = std::str::from_utf8(bytes).map_err(|_| "voice performance must be UTF-8")?;
        let mut lines = text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty() && !l.starts_with('#'));
        if lines.next() != Some(VOICE_PERFORMANCE_SCHEMA) {
            return Err(format!("first line must be {VOICE_PERFORMANCE_SCHEMA}"));
        }
        let mut audio = None;
        let mut full_scale = None;
        let mut gas = None;
        let mut tract_gas = None;
        let mut observer_law = None;
        let mut fold = None;
        let mut wall = None;
        let mut termination = None;
        let mut source = None;
        let mut licence = None;
        let mut microphone = None;
        let mut sections: Vec<TractSection> = Vec::new();
        let mut subglottal: Vec<(f64, f64)> = Vec::new();
        let mut subglottal_airway = None;
        let mut tension_modulation = None;
        let mut tension_flutter = None;
        let mut tension_flutter_second = None;
        let mut aspiration_noise = None;
        let mut jet_noise = None;
        let mut constriction_noise = None;
        let mut pulse_shimmer = None;
        let mut head_sphere = None;
        let mut head_load = None;
        let mut torso_sphere = None;
        let mut fold_mapping = None;
        let mut fold_shape = None;
        let mut fold_damping_upper = None;
        let mut fold_body = None;
        let mut keyframe_drafts: Vec<KeyframeDraft> = Vec::new();
        let mut tract_moves: Vec<(f64, f64, String, u64)> = Vec::new();
        for line in lines {
            let mut tokens = line.split_whitespace();
            let key = tokens.next().unwrap_or_default();
            let rest: Vec<&str> = tokens.collect();
            match key {
                "audio" => {
                    arity(&rest, 3, key)?;
                    let rate = whole(&rest, 0, "audio rate")?;
                    let samples = whole(&rest, 1, "audio samples")?;
                    let ratio = whole(&rest, 2, "mechanics per output")?;
                    if !(8_000..=192_000).contains(&rate) || !(1..=16).contains(&ratio) {
                        return Err("audio needs a rate in 8000..=192000 and a 1..=16 ratio".into());
                    }
                    if samples == 0 || samples > MAX_SECONDS * rate {
                        return Err(format!("audio samples must be in 1..={MAX_SECONDS} s"));
                    }
                    once(&mut audio, (rate as u32, samples, ratio as u32), key)?;
                }
                "full_scale_pa" => {
                    arity(&rest, 1, key)?;
                    let pa = field(&rest, 0, "full scale")?;
                    if pa <= 0.0 {
                        return Err("full_scale_pa must be positive".into());
                    }
                    once(&mut full_scale, pa, key)?;
                }
                "ambient" => {
                    let record_text = rest.join(" ");
                    let mut fields = record_text.split_ascii_whitespace();
                    if fields.clone().next() == Some("moist_air") {
                        fields.next();
                    }
                    let record = read_gas_record(&mut fields, None, sources)
                        .map_err(|e| format!("ambient gas refused: {e}"))?;
                    if fields.next().is_some() {
                        return Err("ambient gas contains extra fields".into());
                    }
                    once(&mut gas, record, key)?;
                }
                "observer_law" => {
                    let law = ObserverLawInput::parse(&rest)
                        .map_err(|why| format!("observer_law refused: {why}"))?;
                    once(&mut observer_law, law, key)?;
                }
                "tract_gas" => {
                    let record_text = rest.join(" ");
                    let mut fields = record_text.split_ascii_whitespace();
                    if fields.clone().next() == Some("moist_air") {
                        fields.next();
                    }
                    let record = read_gas_record(&mut fields, None, sources)
                        .map_err(|e| format!("tract gas refused: {e}"))?;
                    if fields.next().is_some() {
                        return Err("tract gas contains extra fields".into());
                    }
                    once(&mut tract_gas, record, key)?;
                }
                "fold" => {
                    // Ten fields; an optional eleventh sets the upper-section thickness of the
                    // viscous channel (default: Ishizaka-Flanagan's 0.05 cm, `FoldCard`).
                    if rest.len() != 11 {
                        arity(&rest, 10, key)?;
                    }
                    let (two_mass, body_cover) = match rest[0] {
                        "two_mass" => (true, false),
                        "body_cover" => (true, true),
                        "one_dof" => (false, false),
                        _ => return Err("fold model must be `two_mass`, `body_cover` or `one_dof`".into()),
                    };
                    let card = FoldCard {
                        mass_lower_kg: field(&rest, 1, "lower mass")?,
                        mass_upper_kg: field(&rest, 2, "upper mass")?,
                        stiffness_lower_n_m: field(&rest, 3, "lower stiffness")?,
                        stiffness_upper_n_m: field(&rest, 4, "upper stiffness")?,
                        coupling_n_m: field(&rest, 5, "coupling stiffness")?,
                        damping_ratio: field(&rest, 6, "damping ratio")?,
                        damping_ratio_upper: 0.0, // from `fold_damping_upper`, checked after the loop
                        rest_gap_m: field(&rest, 7, "rest gap")?,
                        fold_length_m: field(&rest, 8, "fold length")?,
                        thickness_m: field(&rest, 9, "thickness")?,
                        thickness_upper_m: if rest.len() == 11 {
                            field(&rest, 10, "upper thickness")?
                        } else {
                            FoldCard::two_mass_standard().thickness_upper_m
                        },
                        // Per-field provenance lives in the file's comments.
                        source: FoldCardSource::AuthoredEstimate,
                    };
                    once(&mut fold, (two_mass, body_cover, card), key)?;
                }
                "tract_wall" => {
                    arity(&rest, 3, key)?;
                    let pin = WallPin {
                        surface_density: field(&rest, 0, "wall surface density")?,
                        stiffness_per_area: field(&rest, 1, "wall stiffness")?,
                        resistance: field(&rest, 2, "wall resistance")?,
                    };
                    once(&mut wall, pin, key)?;
                }
                "termination" => {
                    arity(&rest, 1, key)?;
                    let load = match rest[0] {
                        "flanged" => Termination::FlangedOpen,
                        "unflanged" => Termination::UnflangedOpen,
                        "ideal_open" => Termination::IdealOpen,
                        _ => return Err("termination must be flanged, unflanged or ideal_open".into()),
                    };
                    once(&mut termination, load, key)?;
                }
                "tract_source" => {
                    arity(&rest, 1, key)?;
                    once(&mut source, rest[0].to_string(), key)?;
                }
                "tract_licence" => {
                    let text = rest.join(" ");
                    if text.is_empty() {
                        return Err("tract_licence needs text".into());
                    }
                    once(&mut licence, text, key)?;
                }
                "section" => {
                    arity(&rest, 2, key)?;
                    if sections.len() >= MAX_SECTIONS {
                        return Err("too many tract sections".into());
                    }
                    sections.push(TractSection {
                        area_m2: field(&rest, 0, "section area")?,
                        length_m: field(&rest, 1, "section length")?,
                    });
                }
                "microphone" => {
                    arity(&rest, 6, key)?;
                    let receiver = CircularOutletReceiver {
                        position_m: [
                            field(&rest, 0, "microphone x")?,
                            field(&rest, 1, "microphone y")?,
                            field(&rest, 2, "microphone z")?,
                        ],
                        radial_rings: whole(&rest, 3, "radial rings")? as usize,
                        angular_points: whole(&rest, 4, "angular points")? as usize,
                        maximum_frequency_hz: field(&rest, 5, "resolved band")?,
                    };
                    if receiver.position_m[2] <= 0.0 {
                        return Err("the microphone must stand in front of the baffle (z > 0)".into());
                    }
                    once(&mut microphone, receiver, key)?;
                }
                "subglottal" => {
                    arity(&rest, 2, key)?;
                    if subglottal.len() >= MAX_SCHEDULE_POINTS {
                        return Err("too many subglottal pressure points".into());
                    }
                    let point = (field(&rest, 0, "time")?, field(&rest, 1, "pressure")?);
                    if point.0 < 0.0 || point.1 < 0.0 {
                        return Err("subglottal time and pressure must be non-negative".into());
                    }
                    if subglottal.last().is_some_and(|last| point.0 <= last.0) {
                        return Err("subglottal times must be strictly increasing".into());
                    }
                    subglottal.push(point);
                }
                "subglottal_airway" => {
                    arity(&rest, 3, key)?;
                    if rest[0] != SUBGLOTTAL_AIRWAY_TREE || rest[1] != SUBGLOTTAL_AIRWAY_GEOMETRY {
                        return Err(format!(
                            "subglottal_airway takes `{SUBGLOTTAL_AIRWAY_TREE} {SUBGLOTTAL_AIRWAY_GEOMETRY} <generations>`"
                        ));
                    }
                    let generations = whole(&rest, 2, "subglottal_airway generations")?;
                    if !(1..=WEIBEL_SAAB_2021_TABLE_1_B_MM.len() as u64).contains(&generations) {
                        return Err("subglottal_airway generations must be in 1..=24".into());
                    }
                    once(&mut subglottal_airway, generations as usize, key)?;
                }
                "tension_modulation" => {
                    arity(&rest, 2, key)?;
                    let rate = field(&rest, 0, "tension modulation rate")?;
                    let amplitude = field(&rest, 1, "tension modulation amplitude")?;
                    if !(rate > 0.0 && rate <= 30.0) {
                        return Err("tension modulation rate must be in (0, 30] Hz".into());
                    }
                    // 1 -+ a must stay inside the island's admitted tension range (1/4..=4).
                    if !(0.0..=0.5).contains(&amplitude) {
                        return Err("tension modulation amplitude must be in 0..=0.5".into());
                    }
                    once(&mut tension_modulation, (rate, amplitude), key)?;
                }
                "tension_flutter" => {
                    arity(&rest, 4, key)?;
                    let seed = whole(&rest, 0, "tension_flutter seed")?;
                    let centre_hz = field(&rest, 1, "tension_flutter centre")?;
                    let bandwidth_hz = field(&rest, 2, "tension_flutter bandwidth")?;
                    let sd_rel = field(&rest, 3, "tension_flutter sd")?;
                    let card = TensionFlutterCard { seed, centre_hz, bandwidth_hz, sd_rel };
                    card.validate()?;
                    once(&mut tension_flutter, card, key)?;
                }
                "tension_flutter_second" => {
                    arity(&rest, 4, key)?;
                    let seed = whole(&rest, 0, "tension_flutter_second seed")?;
                    let centre_hz = field(&rest, 1, "tension_flutter_second centre")?;
                    let bandwidth_hz = field(&rest, 2, "tension_flutter_second bandwidth")?;
                    let sd_rel = field(&rest, 3, "tension_flutter_second sd")?;
                    let card = TensionFlutterCard { seed, centre_hz, bandwidth_hz, sd_rel };
                    card.validate()?;
                    once(&mut tension_flutter_second, card, key)?;
                }
                "aspiration_noise" => {
                    // The gain has no default (no source gives one for glottal aspiration); the threshold has the
                    // Sondhi-Schroeter value via Elie and Laprie.
                    if !(rest.len() == 2 || rest.len() == 3) {
                        return Err(format!(
                            "`{key}` takes 2 or 3 fields (seed, xi, optional re_crit), got {}; xi has no default",
                            rest.len()
                        ));
                    }
                    let seed = whole(&rest, 0, "aspiration seed")?;
                    let xi = field(&rest, 1, "aspiration gain xi")?;
                    let re_crit = if rest.len() == 3 {
                        field(&rest, 2, "aspiration re_crit")?
                    } else {
                        RE_CRIT_SONDHI_SCHROETER
                    };
                    let card = AspirationCard { seed, xi, re_crit };
                    card.validate().map_err(str::to_string)?;
                    once(&mut aspiration_noise, card, key)?;
                }
                "jet_noise" => {
                    // Every field is required: the gain, the band and the onset are a source's measured numbers (Verge 1995 for a flue
                    // jet), and the card, not this code, says which it uses and where they come from.
                    arity(&rest, 5, key)?;
                    let seed = whole(&rest, 0, "jet_noise seed")?;
                    let gain = field(&rest, 1, "jet_noise gain")?;
                    let strouhal_lo = field(&rest, 2, "jet_noise ST_LO")?;
                    let strouhal_hi = field(&rest, 3, "jet_noise ST_HI")?;
                    let reynolds_onset = field(&rest, 4, "jet_noise RE_ON")?;
                    let card = JetNoiseCard { seed, gain, strouhal_lo, strouhal_hi, reynolds_onset };
                    card.validate().map_err(str::to_string)?;
                    once(&mut jet_noise, card, key)?;
                }
                "constriction_noise" => {
                    // No gain and no gate: the law has no free number. The radius is the duct the constriction opens into, and
                    // the card, not this code, says which it took.
                    arity(&rest, 2, key)?;
                    let seed = whole(&rest, 0, "constriction_noise seed")?;
                    let duct_radius_m = field(&rest, 1, "constriction_noise duct radius")?;
                    let card = ConstrictionNoiseCard { seed, duct_radius_m };
                    card.validate().map_err(str::to_string)?;
                    once(&mut constriction_noise, card, key)?;
                }
                "pulse_shimmer" => {
                    arity(&rest, 2, key)?;
                    let seed = whole(&rest, 0, "pulse_shimmer seed")?;
                    let sd_rel = field(&rest, 1, "pulse_shimmer sd")?;
                    let card = PulseShimmerCard { seed, sd_rel };
                    card.validate()?;
                    once(&mut pulse_shimmer, card, key)?;
                }
                "head_sphere" => {
                    arity(&rest, 2, key)?;
                    let radius_m = field(&rest, 0, "head radius")?;
                    if radius_m <= 0.0 {
                        return Err("the head radius must be positive".into());
                    }
                    let card = HeadSphereCard { radius_m, source: rest[1].to_string() };
                    once(&mut head_sphere, card, key)?;
                }
                "head_load" => {
                    arity(&rest, 1, key)?;
                    if rest[0] != "cap" {
                        return Err("head_load must be `cap`".into());
                    }
                    once(&mut head_load, (), key)?;
                }
                "torso_sphere" => {
                    arity(&rest, 3, key)?;
                    let torso =
                        TorsoSphere { radius_m: field(&rest, 0, "torso radius")?, neck_gap_m: field(&rest, 1, "neck gap")? };
                    torso.validate()?;
                    once(&mut torso_sphere, TorsoSphereCard { torso, source: rest[2].to_string() }, key)?;
                }
                "fold_mapping" => {
                    arity(&rest, 1, key)?;
                    if rest[0] != "bilateral" {
                        return Err("fold_mapping must be `bilateral`".into());
                    }
                    once(&mut fold_mapping, (), key)?;
                }
                "fold_shape" => {
                    arity(&rest, 2, key)?;
                    let mode = match rest[0] {
                        "sine" => FoldMode::Sine,
                        "rigid" => FoldMode::Rigid,
                        other => return Err(format!("fold_shape mode must be `sine` or `rigid`, got {other:?}")),
                    };
                    let rest_gap = match rest[1] {
                        "uniform" => FoldRestGap::Uniform,
                        "tapered" => FoldRestGap::Tapered,
                        other => return Err(format!("fold_shape rest gap must be `uniform` or `tapered`, got {other:?}")),
                    };
                    once(&mut fold_shape, FoldShape { mode, rest: rest_gap }, key)?;
                }
                "fold_damping_upper" => {
                    arity(&rest, 1, key)?;
                    let ratio = field(&rest, 0, "upper damping ratio")?;
                    if ratio < 0.0 {
                        return Err("fold_damping_upper must be non-negative".into());
                    }
                    once(&mut fold_damping_upper, ratio, key)?;
                }
                "fold_body" => {
                    arity(&rest, 4, key)?;
                    let body = FoldBody {
                        mass_kg: field(&rest, 0, "body mass")?,
                        stiffness_n_m: field(&rest, 1, "body stiffness")?,
                        damping_ratio: field(&rest, 2, "body damping ratio")?,
                        rest_gap_upper_m: field(&rest, 3, "upper rest gap")?,
                    };
                    once(&mut fold_body, body, key)?;
                }
                "keyframe" => {
                    if rest.len() < 3 {
                        return Err("keyframe rows are `keyframe <name> source <id>` and `keyframe <name> licence <text>`".into());
                    }
                    let draft = keyframe_draft(&mut keyframe_drafts, rest[0])?;
                    match rest[1] {
                        "source" => {
                            arity(&rest, 3, "keyframe source")?;
                            once(&mut draft.source, rest[2].to_string(), "keyframe source")?;
                        }
                        "licence" => once(&mut draft.licence, rest[2..].join(" "), "keyframe licence")?,
                        _ => {
                            return Err(
                                "keyframe rows are `keyframe <name> source <id>` and `keyframe <name> licence <text>`".into(),
                            );
                        }
                    }
                }
                "keyframe_section" => {
                    arity(&rest, 3, key)?;
                    let draft = keyframe_draft(&mut keyframe_drafts, rest[0])?;
                    if draft.sections.len() >= MAX_SECTIONS {
                        return Err("too many keyframe sections".into());
                    }
                    draft.sections.push(TractSection {
                        area_m2: field(&rest, 1, "keyframe section area")?,
                        length_m: field(&rest, 2, "keyframe section length")?,
                    });
                }
                "tract_move" => {
                    arity(&rest, 4, key)?;
                    if tract_moves.len() >= MAX_TRACT_MOVES {
                        return Err("too many tract_move rows".into());
                    }
                    tract_moves.push((
                        field(&rest, 0, "tract_move start")?,
                        field(&rest, 1, "tract_move end")?,
                        rest[2].to_string(),
                        whole(&rest, 3, "tract_move steps")?,
                    ));
                }
                other => return Err(format!("unknown key {other:?}")),
            }
        }
        let missing = |what: &str| format!("voice performance needs `{what}`");
        let (output_rate_hz, samples, mechanics_per_output) = audio.ok_or_else(|| missing("audio"))?;
        let (two_mass, body_cover, mut fold) = fold.ok_or_else(|| missing("fold"))?;
        if fold_mapping.is_some() {
            fold = fold.bilateral();
        }
        let fold_body = match (body_cover, fold_body) {
            (true, Some(body)) if fold_mapping.is_some() => Some(FoldBody::bilateral(body)),
            (true, Some(body)) => Some(body),
            (true, None) => return Err(missing("fold_body")),
            (false, Some(_)) => return Err("fold_body needs a body_cover fold".into()),
            (false, None) => None,
        };
        if body_cover && fold_shape.is_some() {
            return Err("fold_shape needs a two_mass fold, not body_cover".into());
        }
        match (two_mass, fold_damping_upper) {
            (true, Some(ratio)) => fold.damping_ratio_upper = ratio,
            (true, None) => return Err(missing("fold_damping_upper")),
            (false, Some(_)) => return Err("fold_damping_upper needs a two_mass fold".into()),
            (false, None) => {}
        }
        if !two_mass && fold_shape.is_some() {
            return Err("fold_shape needs a two_mass fold".into());
        }
        if sections.is_empty() {
            return Err(missing("section"));
        }
        if subglottal.is_empty() {
            return Err(missing("subglottal"));
        }
        if tension_flutter_second.is_some() && tension_flutter.is_none() {
            return Err("tension_flutter_second needs tension_flutter: it adds to that row".into());
        }
        if head_load.is_some() {
            if head_sphere.is_none() {
                return Err(HEAD_LOAD_NEEDS_SPHERE.into());
            }
            if !matches!(termination, Some(Termination::FlangedOpen)) {
                return Err(HEAD_LOAD_NEEDS_FLANGED.into());
            }
        }
        if torso_sphere.is_some() && head_sphere.is_none() {
            return Err(TORSO_NEEDS_SPHERE.into());
        }
        // In the order the fields were always refused: full scale, then the gas.
        let full_scale_pa = full_scale.ok_or_else(|| missing("full_scale_pa"))?;
        let gas = gas.ok_or_else(|| missing("ambient"))?;
        if let Some(tract) = &tract_gas {
            if head_load.is_none() {
                return Err(TRACT_GAS_NEEDS_HEAD_LOAD.into());
            }
            if tract.state().pressure != gas.state().pressure {
                return Err("tract_gas and ambient need one static pressure: mean flow is not modeled".into());
            }
        }
        let (keyframes, tract_moves) = admit_motion(
            &sections,
            keyframe_drafts,
            tract_moves,
            samples as f64 / f64::from(output_rate_hz),
        )?;
        Ok(Self {
            output_rate_hz,
            samples,
            mechanics_per_output,
            full_scale_pa,
            gas,
            tract_gas,
            observer_law,
            two_mass,
            fold,
            fold_shape,
            fold_body,
            wall,
            termination: termination.ok_or_else(|| missing("termination"))?,
            sections,
            tract_source: source.ok_or_else(|| missing("tract_source"))?,
            tract_licence: licence.ok_or_else(|| missing("tract_licence"))?,
            microphone: microphone.ok_or_else(|| missing("microphone"))?,
            subglottal,
            subglottal_airway,
            tension_modulation,
            tension_flutter,
            tension_flutter_second,
            aspiration_noise,
            jet_noise,
            constriction_noise,
            pulse_shimmer,
            head_sphere,
            head_load: head_load.is_some(),
            torso_sphere,
            keyframes,
            tract_moves,
        })
    }

    /// The exact state of the exterior, which also fills the tract unless a `tract_gas` row is present.
    #[must_use]
    pub fn ambient_gas(&self) -> &GasState {
        self.gas.state()
    }

    /// Original source receipts remain attached to the numerical gas view.
    #[must_use]
    pub const fn gas_record(&self) -> &GasRecord {
        &self.gas
    }

    /// The exact state filling the tract: the `tract_gas` row's, else [`Self::ambient_gas`].
    #[must_use]
    pub fn tract_gas(&self) -> &GasState {
        self.tract_gas.as_ref().map_or_else(|| self.gas.state(), GasRecord::state)
    }

    /// The `tract_gas` row's record with its source receipts; `None` when the tract holds the exterior's gas.
    #[must_use]
    pub const fn tract_gas_record(&self) -> Option<&GasRecord> {
        self.tract_gas.as_ref()
    }

    /// Realize the island, the tract's mouth-flow line and the microphone.
    ///
    /// # Errors
    /// Any owner's refusal: the chart licence check, the fold card, a line
    /// that does not converge or is not passive, or an under-resolved band.
    pub fn build(&self) -> Result<VoicePerformance, String> {
        let workers = std::thread::available_parallelism().map_or(1, |p| p.get()).clamp(1, MAX_SHAPE_WORKERS);
        self.build_with_shape_workers(workers)
    }

    /// [`Self::build`] with a moving tract's further shapes built on `workers` threads; every count gives the same
    /// performance.
    fn build_with_shape_workers(&self, workers: usize) -> Result<VoicePerformance, String> {
        // `tract_gas`: the tract's lines, formants and glottal flow read `gas`; the lips' load, the receiver and the
        // spheres read `exterior`. Without the row the two are one state.
        let gas = self.tract_gas();
        let exterior = self.ambient_gas();
        let observer = match &self.observer_law {
            Some(declared) => {
                let (law, relaxation) = resolve_observer_law(
                    declared, exterior, core::slice::from_ref(&self.gas),
                ).map_err(|error| format!("voice observer law refused: {error}"))?;
                Some(GasPropagation {
                    law, relaxation, limits: GasPropagationLimits::surface_reference(),
                })
            }
            None => None,
        };
        let chart = TractChart::try_new(
            self.sections.clone(),
            &self.tract_source,
            &self.tract_licence,
        )
        .map_err(|e| format!("tract chart refused: {e:?}"))?;
        let duct = chart.to_duct();
        let mechanics_rate = self.output_rate_hz * self.mechanics_per_output;
        if self.torso_sphere.is_some() && self.head_sphere.is_none() {
            return Err(TORSO_NEEDS_SPHERE.into());
        }
        // `head_load`: one table on the finest rung's bins, read by the island's line, the mouth-flow line and the
        // passive formants alike, with the torso's change when there is a torso (one geometry for the lips and the
        // microphone). Without the row each takes the declared termination exactly as before.
        let head_load = match (self.head_load, &self.head_sphere) {
            (false, _) => None,
            (true, None) => return Err(HEAD_LOAD_NEEDS_SPHERE.into()),
            (true, Some(_)) if !matches!(self.termination, Termination::FlangedOpen) => {
                return Err(HEAD_LOAD_NEEDS_FLANGED.into());
            }
            (true, Some(card)) => Some(
                match &self.torso_sphere {
                    None => head_mouth_load(card.radius_m, &duct, exterior, mechanics_rate, GLOTTAL_LINE_MAX_TAPS),
                    Some(t) => head_torso_mouth_load(card.radius_m, t.torso, self.microphone.maximum_frequency_hz, &duct,
                        exterior, mechanics_rate, GLOTTAL_LINE_MAX_TAPS),
                }
                .map_err(|e| format!("head load refused: {e}"))?,
            ),
        };
        let mouth_load = head_load.as_ref().map_or(MouthLoad::Analytic(self.termination), MouthLoad::Tabulated);
        // With `tract_gas` both tract lines check the table against the exterior it was computed in.
        let split = match (&self.tract_gas, head_load.as_ref()) {
            (None, _) => None,
            (Some(_), Some(table)) => Some(table),
            (Some(_), None) => return Err(TRACT_GAS_NEEDS_HEAD_LOAD.into()),
        };
        // `subglottal_airway`: Weibel's tree below the folds in the tract's gas, every tube carrying the card's
        // `tract_wall` (rigid without one), on the mechanics clock. Without the row the island takes no airway.
        let airway = match self.subglottal_airway {
            Some(generations) => Some(
                AirwayChain::weibel_saab_2021_b(generations, self.wall)
                    .and_then(|chain| {
                        let ladder = LineLadder { max_taps: GLOTTAL_LINE_MAX_TAPS, ..LineLadder::default() };
                        UpstreamAirway::realize(&chain, gas, mechanics_rate, &ladder)
                    })
                    .map_err(|e| format!("subglottal airway refused: {e}"))?,
            ),
            None => None,
        };
        let island = match split {
            None => GlottalIsland::new_with_wall(
                self.fold,
                self.two_mass,
                &duct,
                gas,
                mouth_load,
                mechanics_rate,
                self.wall.as_ref(),
            ),
            Some(table) => GlottalIsland::new_with_wall_exterior(
                self.fold,
                self.two_mass,
                &duct,
                gas,
                exterior,
                table,
                mechanics_rate,
                self.wall.as_ref(),
            ),
        }
        .and_then(|island| match self.fold_shape {
            Some(shape) => island.with_fold_shape(shape),
            None => Ok(island),
        })
        .and_then(|island| match self.fold_body {
            Some(body) => island.with_body(body),
            None => Ok(island),
        })
        .and_then(|island| match airway {
            Some(airway) => island.with_upstream(airway),
            None => Ok(island),
        })
        .map_err(|e| format!("glottal island refused: {e:?}"))?;
        let (mouth, receipt) = match split {
            None => mouth_flow_line(
                &duct,
                gas,
                mouth_load,
                mechanics_rate,
                island.port_impedance(),
                self.wall.as_ref(),
                &LineLadder::default(),
            ),
            Some(table) => mouth_flow_line_exterior(
                &duct,
                gas,
                exterior,
                table,
                mechanics_rate,
                island.port_impedance(),
                self.wall.as_ref(),
                &LineLadder::default(),
            ),
        }
        .map_err(|e| format!("mouth-flow line refused: {e}"))?;
        // A moving tract: every further shape's two lines, and the steps at which the render switches to them. Without
        // `tract_move` rows nothing here runs and both lines stay the fixed ones.
        let mut mouth = ShapeLine::Fixed(mouth);
        let (island, motion) = if self.tract_moves.is_empty() {
            (island, TractMotion::default())
        } else {
            let mut motion = self.realize_motion(&chart, mechanics_rate, island.port_impedance(), workers)?;
            mouth = mouth
                .to_moving(core::mem::take(&mut motion.mouth_taps))
                .map_err(|e| format!("moving mouth-flow line refused: {e}"))?;
            let island = island
                .with_tract_keyframes(core::mem::take(&mut motion.island_taps))
                .map_err(|e| format!("moving tract line refused: {e:?}"))?;
            (island, motion)
        };
        let mouth_area = self.sections.last().map_or(0.0, |s| s.area_m2);
        let mouth_radius_m = (mouth_area / core::f64::consts::PI).sqrt();
        let medium = RayleighMedium {
            density: exterior.density,
            sound_speed: exterior.sound_speed,
        };
        let (receiver, observer_paths) = match observer.as_ref() {
            Some(propagation) => {
                let opening = CircularOpening {
                    center_m: [0.0, 0.0], radius_m: mouth_radius_m,
                };
                let (receiver, report) = BaffledPressure::circular_openings_multirate_gas_differentiating(
                    &[opening], mechanics_rate, self.output_rate_hz, self.microphone,
                    exterior, propagation, &mut || false,
                ).map_err(|error| format!("voice observer paths refused: {error}"))?;
                (receiver, Some(report))
            }
            None => (
                BaffledPressure::circular_outlet_multirate_differentiating(
                    mouth_radius_m, mechanics_rate, self.output_rate_hz,
                    self.microphone, medium,
                )?,
                None,
            ),
        };
        let head = match &self.head_sphere {
            Some(card) => Some(
                match &self.torso_sphere {
                    None => HeadSphere::new(card.radius_m, mouth_radius_m, self.microphone, medium, self.output_rate_hz),
                    Some(t) => HeadSphere::with_torso(card.radius_m, mouth_radius_m, self.microphone, medium,
                        self.output_rate_hz, t.torso),
                }
                .map_err(|e| format!("head sphere refused: {e}"))?,
            ),
            None => None,
        };
        let passive_formants_hz = match &head_load {
            Some(table) => chart.tmm_formants_tabulated(gas, self.wall.as_ref(), table, 4),
            None => chart.tmm_formants_terminated(gas, self.wall.as_ref(), self.termination, 4),
        }
        .map_err(|e| format!("{e:?}"));
        let info = VoiceInfo {
            output_rate_hz: self.output_rate_hz,
            mechanics_rate_hz: mechanics_rate,
            samples: self.samples,
            full_scale_pa: self.full_scale_pa,
            two_mass: self.two_mass,
            sections: self.sections.len(),
            tract_length_m: self.sections.iter().map(|s| s.length_m).sum(),
            glottal_radius_m: (self.sections[0].area_m2 / core::f64::consts::PI).sqrt(),
            mouth_radius_m,
            sound_speed_m_s: gas.sound_speed,
            has_wall: self.wall.is_some(),
            mouth_flow_taps: receipt.n_fft,
            mouth_flow_rungs: receipt.rungs.len(),
            passive_formants_hz,
            head_sphere_radius_m: self.head_sphere.as_ref().map(|card| card.radius_m),
            head_load: head_load.is_some(),
            torso_sphere: self.torso_sphere.as_ref().map(|card| card.torso),
            tract_shapes: motion.switches.len() + 1,
            keyframe_formants_hz: motion.keyframe_formants_hz,
            subglottal_airway: island.upstream().map(|airway| airway.receipt().clone()),
        };
        let aspiration = match self.aspiration_noise {
            Some(card) => Some(GlottalAspiration::new(card).map_err(str::to_string)?),
            None => None,
        };
        let jet_noise = match self.jet_noise {
            Some(card) => Some(GlottalJetNoise::new(card).map_err(str::to_string)?),
            None => None,
        };
        let constriction_noise = match self.constriction_noise {
            Some(card) => Some(ConstrictionNoise::new(card, 1.0 / f64::from(mechanics_rate)).map_err(str::to_string)?),
            None => None,
        };
        Ok(VoicePerformance {
            island,
            mouth,
            receiver,
            observer_paths,
            head,
            aspiration,
            jet_noise,
            constriction_noise,
            pulse_shimmer: self.pulse_shimmer.map(PulseShimmer::new),
            gas_density: gas.density,
            gas_viscosity: gas.dynamic_viscosity,
            fold_length_m: self.fold.fold_length_m,
            schedule: Schedule {
                points: self.subglottal.clone(),
                cursor: 0,
            },
            tension_modulation: self.tension_modulation,
            tension_flutter: self.tension_flutter.map(|card| TensionFlutter::new(card, self.output_rate_hz)),
            tension_flutter_second: self
                .tension_flutter_second
                .map(|card| TensionFlutter::on_tile(card, self.output_rate_hz, 1)),
            trace: vec![0.0; self.mechanics_per_output as usize],
            ratio: self.mechanics_per_output as usize,
            mechanics_rate,
            info,
            steps: 0,
            completed: 0,
            poisoned: false,
            peak_flow_m3_s: 0.0,
            min_flow_m3_s: 0.0,
            closed_steps: 0,
            tract_switches: motion.switches,
            next_switch: 0,
        })
    }

    /// Every further shape of the moving tract, in the order the render reaches them, with the step each starts at
    /// (module doc). A shape is the morph at `k / steps` from the shape the tract holds to the move's target, with the
    /// card's own first-section area (the glottal section is held); `k = steps` is the target as written.
    ///
    /// The shapes and switch steps come first, in one serial pass; each shape's lines are then built on up to `workers`
    /// threads ([`build_shapes`]). A refusal is the one a shape-by-shape build would meet first: a build refusal of an
    /// earlier shape before a later shape's morph or step refusal.
    fn realize_motion(&self, chart: &TractChart, mechanics_rate: u32, zc: f64, workers: usize) -> Result<TractMotion, String> {
        let glottal_area = chart.sections()[0].area_m2;
        let held = |shape: TractChart| -> Result<TractChart, String> {
            let mut sections = shape.sections().to_vec();
            sections[0].area_m2 = glottal_area;
            TractChart::try_new(sections, shape.source_id(), shape.license())
                .map_err(|e| format!("moving tract shape refused: {e:?}"))
        };
        let keyframes = self
            .keyframes
            .iter()
            .map(|k| TractChart::try_new(k.sections.clone(), &k.source, &k.licence))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format!("keyframe chart refused: {e:?}"))?;
        let rate = f64::from(mechanics_rate);
        let mut motion = TractMotion::default();
        let mut shapes = Vec::new();
        let mut targets = Vec::new();
        let mut cut = None;
        let mut holding = chart.clone();
        'moves: for m in &self.tract_moves {
            let (target, name) = match m.target {
                None => (chart.clone(), "card"),
                Some(at) => (keyframes[at].clone(), self.keyframes[at].name.as_str()),
            };
            let start = (m.start_s * rate).round() as u64;
            let total = ((m.end_s * rate).round() as u64).saturating_sub(start);
            let steps = u64::from(m.steps);
            if total < steps {
                cut = Some("a tract_move needs at least one mechanics step per shape".to_string());
                break;
            }
            for k in 1..=steps {
                let shape = holding
                    .log_area_morph(&target, k as f64 / steps as f64)
                    .map_err(|e| format!("moving tract morph refused: {e:?}"))
                    .and_then(&held);
                match shape {
                    Ok(shape) => shapes.push(shape),
                    Err(e) => {
                        cut = Some(e);
                        break 'moves;
                    }
                }
                let from = start + total * (k - 1) / steps;
                let to = start + total * k / steps;
                motion.switches.push(TractSwitch {
                    step: from,
                    shape: shapes.len(),
                    fade: (to - from) as usize,
                });
            }
            targets.push((shapes.len() - 1, name.to_string()));
            holding = target;
        }
        let built = build_shapes(shapes.len(), workers, |i| self.shape_taps(&shapes[i], mechanics_rate, zc))?;
        if let Some(e) = cut {
            return Err(e);
        }
        let mut formants = Vec::with_capacity(built.len());
        for (island_taps, mouth_taps, f) in built {
            motion.island_taps.push(island_taps);
            motion.mouth_taps.push(mouth_taps);
            formants.push(Some(f));
        }
        for (at, name) in targets {
            let f = formants[at].take().ok_or("a moving tract target was reached twice")?;
            motion.keyframe_formants_hz.push((name, f));
        }
        Ok(motion)
    }

    /// The island's line taps, the mouth-flow line taps and the passive formants of one further shape, realized with
    /// the card's mouth load, wall and gases exactly as `build` realizes the card's own shape. A shape whose `zc` is
    /// not the card's is refused: one glottal section for every shape.
    #[allow(clippy::type_complexity)]
    fn shape_taps(
        &self,
        chart: &TractChart,
        mechanics_rate: u32,
        zc: f64,
    ) -> Result<(Vec<f64>, Vec<f64>, Result<Vec<f64>, String>), String> {
        let gas = self.tract_gas();
        let exterior = self.ambient_gas();
        let duct = chart.to_duct();
        // `build` has already refused every other combination of these rows.
        let head_load = match (self.head_load, &self.head_sphere) {
            (true, Some(card)) => Some(
                match &self.torso_sphere {
                    None => head_mouth_load(card.radius_m, &duct, exterior, mechanics_rate, GLOTTAL_LINE_MAX_TAPS),
                    Some(t) => head_torso_mouth_load(card.radius_m, t.torso, self.microphone.maximum_frequency_hz,
                        &duct, exterior, mechanics_rate, GLOTTAL_LINE_MAX_TAPS),
                }
                .map_err(|e| format!("moving tract head load refused: {e}"))?,
            ),
            _ => None,
        };
        let mouth_load = head_load.as_ref().map_or(MouthLoad::Analytic(self.termination), MouthLoad::Tabulated);
        let split = if self.tract_gas.is_some() { head_load.as_ref() } else { None };
        let wall = self.wall.as_ref();
        let (island_taps, shape_zc) = match split {
            None => glottal_tract_line_taps(&duct, gas, None, mouth_load, mechanics_rate, wall),
            Some(table) => {
                glottal_tract_line_taps(&duct, gas, Some(exterior), MouthLoad::Tabulated(table), mechanics_rate, wall)
            }
        }
        .map_err(|e| format!("moving tract line refused: {e:?}"))?;
        if shape_zc.to_bits() != zc.to_bits() {
            return Err("every shape of a moving tract needs the card's glottal section: one zc for the island".into());
        }
        let (mouth, _) = match split {
            None => mouth_flow_line(&duct, gas, mouth_load, mechanics_rate, zc, wall, &LineLadder::default()),
            Some(table) => mouth_flow_line_exterior(
                &duct, gas, exterior, table, mechanics_rate, zc, wall, &LineLadder::default(),
            ),
        }
        .map_err(|e| format!("moving mouth-flow line refused: {e}"))?;
        let (_, mouth_taps) = mouth
            .quiescent_fir()
            .ok_or("a moving tract needs every shape's mouth-flow line to be a plain FIR")?;
        let formants = match &head_load {
            Some(table) => chart.tmm_formants_tabulated(gas, wall, table, 4),
            None => chart.tmm_formants_terminated(gas, wall, self.termination, 4),
        }
        .map_err(|e| format!("{e:?}"));
        Ok((island_taps, mouth_taps.to_vec(), formants))
    }
}

/// One further shape's island taps, mouth-flow taps and passive formants (`VoiceDescription::shape_taps`).
type ShapeLines = (Vec<f64>, Vec<f64>, Result<Vec<f64>, String>);

/// `build(0..n)` on up to `workers` scoped threads, each on one contiguous block of shapes, in the manner of the snowman's
/// `par_map` (`baffled_pressure_head_torso`): every shape is built alone by the same serial code and kept at its index, so
/// the result is the same for any worker count. The first refusal in shape order is returned; a block stops at its own
/// first refusal.
fn build_shapes<F>(n: usize, workers: usize, build: F) -> Result<Vec<ShapeLines>, String>
where
    F: Fn(usize) -> Result<ShapeLines, String> + Sync,
{
    let chunk = n.div_ceil(workers.clamp(1, n.max(1))).max(1);
    let mut out: Vec<Option<Result<ShapeLines, String>>> = (0..n).map(|_| None).collect();
    std::thread::scope(|scope| {
        for (c, slots) in out.chunks_mut(chunk).enumerate() {
            let build = &build;
            scope.spawn(move || {
                for (i, slot) in slots.iter_mut().enumerate() {
                    let lines = build(c * chunk + i);
                    let refused = lines.is_err();
                    *slot = Some(lines);
                    if refused {
                        break;
                    }
                }
            });
        }
    });
    // A slot is empty only after its block's own refusal, which comes first in this scan.
    out.into_iter()
        .map(|slot| slot.unwrap_or_else(|| Err("a moving tract shape was not built".into())))
        .collect()
}

/// A moving tract as realized at build time.
#[derive(Default)]
struct TractMotion {
    /// The island's line taps of shapes 1.. (shape 0 is the island's own).
    island_taps: Vec<Vec<f64>>,
    /// The mouth-flow line taps of shapes 1...
    mouth_taps: Vec<Vec<f64>>,
    /// When the render starts each shape's crossfade, in step order.
    switches: Vec<TractSwitch>,
    /// Passive formants of each move's target, by name.
    keyframe_formants_hz: Vec<(String, Result<Vec<f64>, String>)>,
}

/// The physical step at which the render starts crossfading to `shape`, and over how many steps.
#[derive(Debug, Clone, Copy)]
struct TractSwitch {
    step: u64,
    shape: usize,
    fade: usize,
}

/// What a realized performance reports about itself.
#[derive(Debug, Clone, PartialEq)]
pub struct VoiceInfo {
    /// Output clock [Hz].
    pub output_rate_hz: u32,
    /// Physical clock [Hz].
    pub mechanics_rate_hz: u32,
    /// Output samples declared.
    pub samples: u64,
    /// Pascals at PCM full scale.
    pub full_scale_pa: f64,
    /// Two-mass or body-cover (`true`), or one-mass fold.
    pub two_mass: bool,
    /// Tract section count.
    pub sections: usize,
    /// Sum of section lengths [m].
    pub tract_length_m: f64,
    /// Radius of the glottal-end section [m].
    pub glottal_radius_m: f64,
    /// Radius of the mouth-end section, also the radiating opening [m].
    pub mouth_radius_m: f64,
    /// Speed of sound in the gas filling the tract (the `tract_gas` row's, else `ambient`'s) [m/s].
    pub sound_speed_m_s: f64,
    /// Whether a yielding wall was declared.
    pub has_wall: bool,
    /// Converged length of the mouth-flow line.
    pub mouth_flow_taps: usize,
    /// Rungs the mouth-flow ladder traversed.
    pub mouth_flow_rungs: usize,
    /// Passive resonances of exactly this tract, mouth load and wall: the
    /// first four `|Z_in|` peaks seen from a closed glottis
    /// ([`TractChart::tmm_formants_terminated`]), or why it refused. The
    /// glottis opens and closes, so this is the rigid-end limit that the
    /// formants of the radiated sound are compared with, not a prediction of
    /// them.
    pub passive_formants_hz: Result<Vec<f64>, String>,
    /// The `head_sphere` row's radius [m]; `None` on an infinite baffle.
    pub head_sphere_radius_m: Option<f64>,
    /// Whether the lips see the head sphere's radiation load (`head_load`) instead of the declared termination.
    pub head_load: bool,
    /// The `torso_sphere` row's torso; `None` for the head alone.
    pub torso_sphere: Option<TorsoSphere>,
    /// Shapes of the tract's lines, the card's own included: 1 for a fixed tract.
    pub tract_shapes: usize,
    /// Passive formants of each `tract_move` target, by keyframe name (`card` for the card's own shape), as
    /// [`VoiceInfo::passive_formants_hz`] reads them; empty for a fixed tract.
    pub keyframe_formants_hz: Vec<(String, Result<Vec<f64>, String>)>,
    /// The `subglottal_airway` row's realized airway; `None` without the row.
    pub subglottal_airway: Option<AirwayReceipt>,
}

/// Source-side summaries over the physical substeps rendered so far.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VoiceStats {
    /// Largest glottal volume flow [m^3/s].
    pub peak_glottal_flow_m3_s: f64,
    /// Smallest glottal volume flow [m^3/s], negative when the flow reverses.
    pub min_glottal_flow_m3_s: f64,
    /// Fraction of physical substeps with the folds in contact: with a `fold_shape`, any strip of either section.
    pub closed_fraction: f64,
    /// Most passes of the two-mass fold step's Newton solve in any one physical substep
    /// ([`GlottalIsland::force_passes`]).
    pub max_force_passes: u32,
    /// Mean passes of that solve per physical substep.
    pub mean_force_passes: f64,
    /// Gains `pulse_shimmer` has drawn, one per opening of the glottis; 0 without the key.
    pub pulse_shimmer_draws: u64,
}

struct Schedule {
    points: Vec<(f64, f64)>,
    cursor: usize,
}

impl Schedule {
    fn at(&mut self, t: f64) -> f64 {
        let last = self.points.len() - 1;
        if t <= self.points[0].0 {
            return self.points[0].1;
        }
        while self.cursor < last && t >= self.points[self.cursor + 1].0 {
            self.cursor += 1;
        }
        if self.cursor >= last {
            return self.points[last].1;
        }
        let (t0, y0) = self.points[self.cursor];
        let (t1, y1) = self.points[self.cursor + 1];
        y0 + (y1 - y0) * (t - t0) / (t1 - t0)
    }
}

/// Kernel id of the `tension_flutter` noise stream in `fs_rand`'s `(seed, kernel, tile)` key: "flut".
pub const TENSION_FLUTTER_KERNEL: u32 = 0x666c_7574;

/// A `tension_flutter` row: the filtered-noise model of Ternstrom and Friberg (1989), its two parameters and the
/// standard deviation it is scaled to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensionFlutterCard {
    /// Seed of the noise stream.
    pub seed: u64,
    /// Natural frequency of the resonant filter [Hz].
    pub centre_hz: f64,
    /// Bandwidth of the resonant filter [Hz]; its poles decay at `pi * bandwidth_hz`.
    pub bandwidth_hz: f64,
    /// Stationary standard deviation of the relative stiffness change.
    pub sd_rel: f64,
}

impl TensionFlutterCard {
    fn validate(&self) -> Result<(), String> {
        if !(self.centre_hz > 0.0 && self.centre_hz <= 30.0) {
            return Err("tension_flutter centre must be in (0, 30] Hz".into());
        }
        // At twice the centre the poles turn real and the filter no longer resonates.
        if !(self.bandwidth_hz > 0.0 && self.bandwidth_hz < 2.0 * self.centre_hz) {
            return Err("tension_flutter bandwidth must be positive and below twice the centre".into());
        }
        // A margin, not a measurement: 0.1 keeps the island's admitted tension range (1/4..=4)
        // 7.5 standard deviations away.
        if !(0.0..=0.1).contains(&self.sd_rel) {
            return Err("tension_flutter sd must be in 0..=0.1".into());
        }
        Ok(())
    }
}

/// White noise through a two-pole resonator, `y[n] = a1 y[n-1] + a2 y[n-2] + gain w[n]`, one step per output sample.
struct TensionFlutter {
    white: Stream,
    a1: f64,
    a2: f64,
    gain: f64,
    y1: f64,
    y2: f64,
}

impl TensionFlutter {
    fn new(card: TensionFlutterCard, output_rate_hz: u32) -> Self {
        Self::on_tile(card, output_rate_hz, 0)
    }

    /// `tile` is the stream's place in the row's key: 0 for `tension_flutter`, 1 for `tension_flutter_second`.
    fn on_tile(card: TensionFlutterCard, output_rate_hz: u32, tile: u32) -> Self {
        let rate = f64::from(output_rate_hz);
        let half_bandwidth = card.bandwidth_hz / 2.0;
        let damped_hz = det::sqrt(card.centre_hz * card.centre_hz - half_bandwidth * half_bandwidth);
        let r = det::exp(-core::f64::consts::PI * card.bandwidth_hz / rate);
        let theta = core::f64::consts::TAU * damped_hz / rate;
        // The stationary variance of the unit-innovation recursion is
        // (1 - a2) / ((1 + a2) (1 - a1 - a2) (1 + a1 - a2)); these forms of the last two factors keep their digits with
        // both poles this close to 1.
        let half_sin = det::sin(theta / 2.0);
        let chord = 4.0 * r * half_sin * half_sin;
        let low = (1.0 - r) * (1.0 - r) + chord;
        let high = (1.0 + r) * (1.0 + r) - chord;
        let variance = (1.0 + r * r) / ((1.0 - r * r) * low * high);
        let mut flutter = Self {
            white: StreamKey { seed: card.seed, kernel: TENSION_FLUTTER_KERNEL, tile }.stream(),
            a1: 2.0 * r * det::cos(theta),
            a2: -r * r,
            gain: card.sd_rel / det::sqrt(variance),
            y1: 0.0,
            y2: 0.0,
        };
        // Ten decay times from rest, so the performance starts in the stationary state and not at y = 0.
        let warm_up = (10.0 * rate / (core::f64::consts::PI * card.bandwidth_hz)).ceil() as u64;
        for _ in 0..warm_up {
            flutter.next_relative();
        }
        flutter
    }

    fn next_relative(&mut self) -> f64 {
        let y = self.a1 * self.y1 + self.a2 * self.y2 + self.gain * self.white.next_normal();
        self.y2 = self.y1;
        self.y1 = y;
        y
    }
}

/// Seed-stream kernel of `pulse_shimmer` ("shim").
pub const PULSE_SHIMMER_KERNEL: u32 = 0x7368_696d;

/// A `pulse_shimmer` row: the seed and the standard deviation of the per-cycle relative gain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PulseShimmerCard {
    /// Seed of the gain stream.
    pub seed: u64,
    /// Standard deviation of `g`, drawn once per glottal cycle; the wave is multiplied by `1 + g`.
    pub sd_rel: f64,
}

impl PulseShimmerCard {
    fn validate(&self) -> Result<(), String> {
        // 0.05 puts a gain of zero 20 standard deviations away.
        if !(0.0..=0.05).contains(&self.sd_rel) {
            return Err("pulse_shimmer sd must be in 0..=0.05".into());
        }
        Ok(())
    }
}

/// The gain of the current cycle: 1 until the first opening, then `1 + g` with a new `g` at each opening.
struct PulseShimmer {
    white: Stream,
    sd_rel: f64,
    gain: f64,
    was_shut: bool,
    draws: u64,
}

impl PulseShimmer {
    fn new(card: PulseShimmerCard) -> Self {
        Self {
            white: StreamKey { seed: card.seed, kernel: PULSE_SHIMMER_KERNEL, tile: 0 }.stream(),
            sd_rel: card.sd_rel,
            gain: 1.0,
            was_shut: false,
            draws: 0,
        }
    }

    /// One mechanics step's frame: a step that ends open after one that ended shut starts a cycle.
    fn observe(&mut self, shut: bool) {
        if self.was_shut && !shut {
            self.gain = 1.0 + self.sd_rel * self.white.next_normal();
            self.draws += 1;
        }
        self.was_shut = shut;
    }
}

/// A realized performance: the fold, the tract and the microphone, advanced
/// together on one physical clock.
pub struct VoicePerformance {
    island: GlottalIsland,
    mouth: ShapeLine,
    /// A moving tract's switches in step order, and the next one due; empty for a fixed tract.
    tract_switches: Vec<TractSwitch>,
    next_switch: usize,
    receiver: BaffledPressure,
    observer_paths: Option<GasPropagationReport>,
    /// Opt-in head: the baffled pressure through the rigid sphere's causal FIR.
    head: Option<HeadSphere>,
    schedule: Schedule,
    tension_modulation: Option<(f64, f64)>,
    tension_flutter: Option<TensionFlutter>,
    tension_flutter_second: Option<TensionFlutter>,
    /// Opt-in aspiration source, with the gas properties and fold length its law needs.
    aspiration: Option<GlottalAspiration>,
    /// Opt-in glottal jet-noise source (Verge's fraction on the glottal jet, an extrapolation).
    jet_noise: Option<GlottalJetNoise>,
    /// Opt-in glottal constriction noise (Oldham and Ukpoho's printed law, an extrapolation).
    constriction_noise: Option<ConstrictionNoise>,
    /// Opt-in per-cycle gain on the wave the tract receives.
    pulse_shimmer: Option<PulseShimmer>,
    gas_density: f64,
    gas_viscosity: f64,
    fold_length_m: f64,
    trace: Vec<f64>,
    ratio: usize,
    mechanics_rate: u32,
    info: VoiceInfo,
    steps: u64,
    completed: u64,
    poisoned: bool,
    peak_flow_m3_s: f64,
    min_flow_m3_s: f64,
    closed_steps: u64,
}

fn refusal(what: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(what))
}

impl VoicePerformance {
    /// The realized description.
    #[must_use]
    pub const fn info(&self) -> &VoiceInfo {
        &self.info
    }

    /// Complete exterior gas and original per-path filter receipts; no report
    /// is synthesized for the legacy lossless microphone.
    #[must_use]
    pub fn observer_gas_report(&self) -> Option<&GasPropagationReport> {
        self.observer_paths.as_ref()
    }

    /// Existing microphone filter latency, in output samples, with any common
    /// material-path padding and the head FIR's when present. Physical
    /// source-to-listener travel remains in the pressure and is not included.
    #[must_use]
    pub fn filter_delay_samples(&self) -> f64 {
        self.receiver.filter_delay_samples() + self.head.as_ref().map_or(0.0, |head| head.latency_samples() as f64)
    }

    /// The opt-in glottal jet noise's open-gate mechanics steps so far by how much of the band each could represent at the
    /// mechanics clock; `None` without the `jet_noise` row.
    #[must_use]
    pub fn jet_noise_coverage(&self) -> Option<crate::turbulence_noise::JetNoiseCoverage> {
        self.jet_noise.as_ref().map(GlottalJetNoise::coverage)
    }

    /// Source-side summaries so far.
    #[must_use]
    pub fn stats(&self) -> VoiceStats {
        let (max_force_passes, all_passes) = self.island.force_passes();
        VoiceStats {
            peak_glottal_flow_m3_s: self.peak_flow_m3_s,
            min_glottal_flow_m3_s: self.min_flow_m3_s,
            closed_fraction: if self.steps == 0 {
                0.0
            } else {
                self.closed_steps as f64 / self.steps as f64
            },
            max_force_passes,
            mean_force_passes: if self.steps == 0 {
                0.0
            } else {
                all_passes as f64 / self.steps as f64
            },
            pulse_shimmer_draws: self.pulse_shimmer.as_ref().map_or(0, |row| row.draws),
        }
    }

    fn sample(&mut self) -> Result<f64, String> {
        for j in 0..self.ratio {
            let t = self.steps as f64 / f64::from(self.mechanics_rate);
            // The gesture is far slower than the audio clock, so the tension is set once per output sample.
            if j == 0 {
                let vibrato = self
                    .tension_modulation
                    .map(|(rate_hz, amplitude)| amplitude * det::sin(core::f64::consts::TAU * rate_hz * t));
                let flutter = self.tension_flutter.as_mut().map(|flutter| {
                    let first = flutter.next_relative();
                    // No sum with zero when the row is absent: such a card keeps the first row's value bit for bit.
                    match self.tension_flutter_second.as_mut() {
                        Some(second) => first + second.next_relative(),
                        None => first,
                    }
                });
                // The same: only a card with both gestures takes the sum, and only such a card was refused before it.
                let moved = match (vibrato, flutter) {
                    (Some(v), Some(y)) => Some(v + y),
                    (one, None) | (None, one) => one,
                };
                if let Some(moved) = moved {
                    self.island
                        .set_tension(1.0 + moved)
                        .map_err(|e| format!("fold tension refused: {e:?}"))?;
                }
            }
            // A moving tract switches both lines on the same physical step, before the island scatters into them.
            if let Some(switch) = self.tract_switches.get(self.next_switch)
                && switch.step == self.steps
            {
                self.island
                    .switch_tract(switch.shape, switch.fade)
                    .map_err(|e| format!("moving tract switch refused: {e:?}"))?;
                self.mouth
                    .switch(switch.shape, switch.fade)
                    .map_err(|e| format!("moving mouth-flow switch refused: {e}"))?;
                self.next_switch += 1;
            }
            let gesture = self.schedule.at(t);
            let mut drive = gesture;
            // The same sources summed alone, for an airway below the folds: there the gesture is the lung pressure and
            // the sources stay in series at the glottis. Without the row this sum is never read.
            let mut series = 0.0;
            if let Some(source) = self.constriction_noise.as_mut() {
                // K(St) times the gesture's own pressure, not the noisy drive: the law's blowing pressure is a held value. Series
                // like the sources below; zero is not added, so a silent source leaves the drive bit for bit what it was.
                let noise = source.pressure_pa(gesture, self.gas_density);
                if noise != 0.0 {
                    drive += noise;
                    series += noise;
                }
            }
            if let Some(source) = self.jet_noise.as_mut() {
                // The jet cannot be faster than the Bernoulli speed of the whole driving pressure. A series source like the
                // aspiration one (the loop sees p_sub + p_n); zero is not added, so a shut gate leaves the drive bit for bit
                // what it was.
                let max_speed = if drive > 0.0 { det::sqrt(2.0 * drive / self.gas_density) } else { 0.0 };
                let noise = source.pressure_pa(
                    self.fold_length_m,
                    max_speed,
                    self.gas_density,
                    self.gas_viscosity,
                    1.0 / f64::from(self.mechanics_rate),
                );
                if noise != 0.0 {
                    drive += noise;
                    series += noise;
                }
            }
            if let Some(source) = self.aspiration.as_mut() {
                // A series pressure source downstream of the folds: the loop sees p_sub + p_n. Zero is not added, so a
                // shut gate leaves the drive bit for bit what it was.
                let noise = source.pressure_pa(self.gas_density, self.gas_viscosity, 1.0 / f64::from(self.mechanics_rate));
                if noise != 0.0 {
                    drive += noise;
                    series += noise;
                }
            }
            let frame = if self.island.has_upstream() {
                self.island.step_lung(gesture, series)
            } else {
                self.island.step(drive)
            }
            .map_err(|e| format!("glottal island refused: {e:?}"))?;
            if let Some(source) = self.aspiration.as_mut() {
                source.record(frame.flow_m3_s, self.fold_length_m * frame.gap_m.max(0.0));
            }
            if let Some(source) = self.jet_noise.as_mut() {
                source.record(frame.flow_m3_s, frame.gap_m.max(0.0));
            }
            if let Some(source) = self.constriction_noise.as_mut() {
                source.record(frame.flow_m3_s);
            }
            let mut wave = self.island.incident_wave_pa();
            if let Some(shimmer) = self.pulse_shimmer.as_mut() {
                // After the island's step, so nothing the island sees changes; a gain of 1 multiplies nothing, so a card
                // with `sd_rel` 0 or a fold that never shuts sends the wave bit for bit.
                shimmer.observe(frame.shut);
                if shimmer.gain != 1.0 {
                    wave *= shimmer.gain;
                }
            }
            let flow = self
                .mouth
                .push(wave)
                .map_err(|e| format!("mouth-flow line refused: {e:?}"))?;
            self.trace[j] = flow;
            self.peak_flow_m3_s = self.peak_flow_m3_s.max(frame.flow_m3_s);
            self.min_flow_m3_s = self.min_flow_m3_s.min(frame.flow_m3_s);
            if frame.shut {
                self.closed_steps += 1;
            }
            self.steps += 1;
        }
        let pressure = self.receiver.step_trace(&self.trace)?;
        // The head filters the baffle's pressure into the sphere's; without the row the pressure is what it was.
        let pressure = match self.head.as_mut() {
            Some(head) => head.step(pressure),
            None => pressure,
        };
        self.completed += 1;
        Ok(pressure)
    }
}

impl PressureRenderer for VoicePerformance {
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
        if sample_rate_hz != self.info.output_rate_hz {
            return Err(RenderError::Sizing {
                what: "voice performance output clock differs from the requested clock",
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
                what: "request runs past the samples the voice performance declares",
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
                what: "voice callback exceeds its admitted block length",
            });
        }
        self.validate_sample_count(output.len() as u64)?;
        for slot in output.iter_mut() {
            match self.sample() {
                Ok(pressure) => *slot = pressure,
                Err(what) => {
                    // A failed physical callback leaves a prefix advanced.
                    self.poisoned = true;
                    return Err(refusal(what));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "render/schedule/reed/field_input/tests/support.rs"]
mod gas_test_support;

#[cfg(test)]
mod tests {
    use super::*;

    // Test fixture, not data: a short uniform-ish tract and fold values scaled
    // from the two-mass pack so the build exercises every owner.
    fn fixture(extra: &str) -> String {
        let mut text = String::from(
            "frankensim-voice-performance-v1\n\
             audio 48000 2400 2\n\
             full_scale_pa 2.0\n\
             ambient moist_air 310.15 101325 1.0\n\
             fold two_mass 1.25e-4 2.5e-5 128 12.8 40 0.1 3.57e-4 1.4e-2 2.5e-3\n\
             fold_damping_upper 0.6\n\
             termination flanged\n\
             tract_source test-fixture\n\
             tract_licence test fixture, not data\n\
             microphone 0 0 0.3 4 16 4000\n\
             subglottal 0 0\n\
             subglottal 0.02 900\n",
        );
        for area in [0.5e-4, 0.6e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4] {
            text.push_str(&format!("section {area:e} 0.02\n"));
        }
        text.push_str(extra);
        text
    }

    fn with_gas(record: &str) -> String {
        fixture("").replace("moist_air 310.15 101325 1.0", record)
    }

    fn pressure(description: &VoiceDescription) -> Vec<f64> {
        let mut performance = description.build().expect("gas reaches all physics owners");
        assert_eq!(
            performance.info().sound_speed_m_s,
            description.tract_gas().sound_speed
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
    fn observer_law_rows_are_explicit_and_refuse_invalid_source_bindings() {
        for row in ["none", "iso9613", "stokes-kirchhoff"] {
            let description = VoiceDescription::from_bytes(
                fixture(&format!("observer_law {row}\n")).as_bytes(),
            ).unwrap();
            assert!(description.observer_law.is_some());
        }
        for rows in [
            "observer_law\n",
            "observer_law guessed-air\n",
            "observer_law iso9613 ignored\n",
            "observer_law none\nobserver_law stokes-kirchhoff\n",
            "observer_law relaxation declared 0\n",
        ] {
            let error = VoiceDescription::from_bytes(fixture(rows).as_bytes()).unwrap_err();
            assert!(error.contains("observer_law"), "{error}");
        }
        let text = with_gas("mixture 310.15 101325 mole 2 N2 0.8 CO2 0.2")
            + "observer_law iso9613\n";
        let description = VoiceDescription::from_bytes(text.as_bytes()).unwrap();
        let error = description.build().err().expect("a named mixture cannot borrow ISO air data");
        assert!(error.contains("observer law refused"), "{error}");
    }

    /// G3: explicitly lossless material paths keep the original microphone's
    /// pressure, physical travel and numerical latency, up to summation roundoff.
    #[test]
    fn declared_lossless_voice_matches_the_original_microphone() {
        let plain = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        let explicit = VoiceDescription::from_bytes(fixture("observer_law none\n").as_bytes()).unwrap();
        let (mut a, mut b) = (plain.build().unwrap(), explicit.build().unwrap());
        assert!(a.observer_gas_report().is_none());
        let report = b.observer_gas_report().unwrap();
        assert_eq!(report.law, crate::air_path::GasAbsorptionModel::Lossless);
        assert_eq!(report.additional_latency_samples, 0);
        assert_eq!(a.filter_delay_samples(), b.filter_delay_samples());
        let (mut before, mut after) = (vec![0.0; 1200], vec![0.0; 1200]);
        a.block(&mut before).unwrap();
        b.block(&mut after).unwrap();
        assert!(before.iter().any(|value| value.abs() > 1e-12));
        for (x, y) in before.iter().zip(&after) {
            assert!((x - y).abs() < 1e-11 * x.abs().max(1.0));
        }
        assert_eq!(a.stats(), b.stats(), "observation must not change the fold mechanics");
    }

    /// G3/G5: the same exterior law reaches every physical path and real
    /// pressure, with a callback-independent state and unchanged source motion.
    #[test]
    fn voice_absorption_uses_original_path_filters_and_is_partition_independent() {
        let plain = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        let selected = VoiceDescription::from_bytes(fixture("observer_law iso9613\n").as_bytes()).unwrap();
        let (mut before, mut whole, mut split) = (
            plain.build().unwrap(), selected.build().unwrap(), selected.build().unwrap(),
        );
        let report = whole.observer_gas_report().unwrap();
        assert_eq!(report.gas, *selected.ambient_gas());
        assert!(report.filters.len() > 1, "the aperture must retain distinct physical paths");
        for filter in &report.filters {
            let expected = crate::air_path::ObserverFilter::for_law(
                report.law, &report.gas, None, filter.range_m, selected.output_rate_hz,
            ).unwrap();
            assert_eq!(filter, expected.receipt());
        }
        let (mut lossless, mut a, mut b) = (vec![0.0; 1200], vec![0.0; 1200], vec![0.0; 1200]);
        before.block(&mut lossless).unwrap();
        whole.block(&mut a).unwrap();
        for block in b.chunks_mut(37) { split.block(block).unwrap(); }
        assert_eq!(a, b);
        assert_ne!(a, lossless, "declared absorption must reach the actual microphone");
        assert!(a.iter().all(|value| value.is_finite()));
        assert_eq!(before.stats(), whole.stats());
        assert_eq!(whole.stats(), split.stats());
    }

    /// G3: a separate tract gas cannot replace the exterior composition at the
    /// observer. Its molecular reservoirs reach each existing filter exactly once.
    #[test]
    fn voice_molecular_observer_binds_the_exterior_beside_a_separate_tract() {
        let rows = format!(
            "{HEAD_LOAD}tract_gas moist_air 293.15 101325 0.5\n\
             observer_law relaxation declared 3 N2 0 O2 0 CO2 1 vibration 1 960 2 power 10 293.15 0 200 400\n"
        );
        let text = fixture(&rows).replace(
            "ambient moist_air 310.15 101325 1.0",
            "ambient mixture 310.15 101325 mole 3 N2 0.4 O2 0.2 CO2 0.4",
        );
        let description = VoiceDescription::from_bytes(text.as_bytes()).unwrap();
        let performance = description.build().unwrap();
        let report = performance.observer_gas_report().unwrap();
        assert_eq!(report.gas, *description.ambient_gas());
        assert_ne!(report.gas, *description.tract_gas());
        assert_eq!(report.law, crate::air_path::GasAbsorptionModel::MolecularRelaxation);
        let ObserverLawInput::MolecularRelaxation(law) = description.observer_law.as_ref().unwrap() else {
            panic!("the authored relaxation law must be retained");
        };
        let state = law.evaluate(
            &fs_material::relaxation::input::medium_from_record(description.gas_record()).unwrap(),
        ).unwrap();
        for filter in &report.filters {
            let expected = crate::air_path::ObserverFilter::molecular_relaxation(
                &state, filter.range_m, description.output_rate_hz,
            ).unwrap();
            assert_eq!(filter, expected.receipt());
            assert!(filter.absorption.relaxation_omitted.is_none());
        }
    }

    #[test]
    fn g3_shared_mixture_composition_and_pressure_change_the_radiated_voice() {
        let legacy = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        assert_eq!(
            *legacy.ambient_gas(),
            GasState::try_new_moist_air(310.15, 101325.0, 1.0).unwrap()
        );
        let records = [
            "mixture 310.15 101325 mole 2 N2 0.8 O2 0.2",
            "mixture 310.15 101325 mole 2 N2 0.8 CO2 0.2",
            "mixture 310.15 80000 mole 2 N2 0.8 O2 0.2",
        ];
        let mut rendered = Vec::new();
        for record in records {
            let description = VoiceDescription::from_bytes(with_gas(record).as_bytes()).unwrap();
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
    fn g3_nasa9_source_swap_changes_voice_pressure_and_retains_the_original_packs() {
        use fs_material::gas_species::input::ThermochemicalPhase;
        let record = "thermochemical 310.15 101325 mole 1 specimen-gas 1";
        let text = with_gas(record);
        assert!(
            VoiceDescription::from_bytes(text.as_bytes())
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
                VoiceDescription::from_bytes_with_gas_sources(text.as_bytes(), Some(&sources))
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
            "mixture 310.15 101325 mole 2 N2 0.7 O2 0.2",
            "mixture 310.15 101325 mole 1 unavailable 1",
            "mixture 310.15 101325 mole 1 N2 1 extra",
            "moist_air 310.15 101325 1 extra",
            "thermochemical 310.15 101325 mole 1 specimen-gas 1",
        ] {
            assert!(
                VoiceDescription::from_bytes(with_gas(record).as_bytes()).is_err(),
                "{record}"
            );
        }
    }

    #[test]
    fn malformed_and_non_physical_files_are_refused_by_name() {
        for (extra, needle) in [
            ("bogus 1\n", "unknown key"),
            ("audio 48000 10 1\n", "only appear once"),
            ("section 1e-4\n", "takes 2 fields"),
            ("subglottal 0.01 5\n", "strictly increasing"),
            ("microphone 0 0 -1 4 16 1000\n", "in front of the baffle"),
        ] {
            let error = VoiceDescription::from_bytes(fixture(extra).as_bytes()).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        let no_licence = fixture("").replace("tract_licence test fixture, not data\n", "");
        let error = VoiceDescription::from_bytes(no_licence.as_bytes()).expect_err("licence");
        assert!(error.contains("tract_licence"), "{error}");
        assert!(VoiceDescription::from_bytes(b"not a voice file\n").is_err());
    }

    #[test]
    fn unlicensed_data_is_refused_at_build() {
        let text = fixture("").replace("test fixture, not data", "unknown");
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let error = description.build().err().expect("licensing-first refusal");
        assert!(error.contains("tract chart refused"), "{error}");
    }

    #[test]
    fn passive_formants_follow_the_declared_mouth_load() {
        // Eight equal sections are a uniform tube of length L and radius a,
        // closed at the glottis. Ideal-open mouth: F1 = c / (4 L). A flanged
        // mouth lengthens the tube by its end correction, about 0.8216 a.
        let gas = GasState::try_new_moist_air(310.15, 101_325.0, 1.0).expect("gas");
        let area = 2.5e-4;
        let sections = vec![TractSection { area_m2: area, length_m: 0.02 }; 8];
        let chart = TractChart::try_new(sections, "uniform-test-fixture", "test fixture, not data")
            .expect("chart");
        let length = 0.16;
        let radius = (area / core::f64::consts::PI).sqrt();
        let ideal = chart
            .tmm_formants_terminated(&gas, None, Termination::IdealOpen, 1)
            .expect("ideal open")[0];
        let flanged = chart
            .tmm_formants_terminated(&gas, None, Termination::FlangedOpen, 1)
            .expect("flanged")[0];
        let ideal_expected = gas.sound_speed / (4.0 * length);
        let flanged_expected = gas.sound_speed / (4.0 * (length + 0.8216 * radius));
        assert!((ideal / ideal_expected - 1.0).abs() < 0.02, "ideal F1 {ideal} vs {ideal_expected}");
        assert!(
            (flanged / flanged_expected - 1.0).abs() < 0.02,
            "flanged F1 {flanged} vs {flanged_expected}"
        );
        assert!(flanged < ideal, "an end correction must lower the resonance");
        assert_eq!(
            chart.tmm_formants(&gas, None, 1).expect("legacy entry point"),
            vec![ideal],
            "the original oracle is the ideal-open case, unchanged"
        );
    }

    #[test]
    fn subglottal_schedule_interpolates_and_holds() {
        let mut schedule = Schedule {
            points: vec![(0.1, 10.0), (0.2, 30.0)],
            cursor: 0,
        };
        assert_eq!(schedule.at(0.0), 10.0);
        assert!((schedule.at(0.15) - 20.0).abs() < 1.0e-9);
        assert_eq!(schedule.at(0.5), 30.0);
    }

    #[test]
    fn rendering_is_independent_of_callback_partition() {
        let description = VoiceDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        let mut whole = description.build().expect("builds");
        let mut split = description.build().expect("builds");
        let mut a = vec![0.0; 600];
        whole.block(&mut a).expect("whole block");
        let mut b = vec![0.0; 600];
        for chunk in b.chunks_mut(37) {
            split.block(chunk).expect("split block");
        }
        assert!(a.iter().all(|v| v.is_finite()));
        assert!(
            a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()),
            "block partition changed the pressure"
        );
        assert_eq!(whole.samples_rendered(), 600);
        let past = split.validate_sample_count(2400);
        assert!(past.is_err(), "a request past the declared end must refuse");
    }

    // Test fixture, not data: the licensed standard two-mass card (80 / 8 / 25 N/m) at 1400 Pa on a uniform
    // 17.6 cm tube of 1 cm radius with a flanged mouth, 1.3 s long. The tube is the one the island's own
    // tests phonate against; it stands for no vowel.
    fn modulation_text(extra: &str) -> String {
        let mut text = String::from(
            "frankensim-voice-performance-v1\n\
             audio 48000 62400 2\n\
             full_scale_pa 2.0\n\
             ambient moist_air 310.15 101325 1.0\n\
             fold two_mass 1.25e-4 2.5e-5 80 8 25 0.1 3.57e-4 1.4e-2 2.5e-3\n\
             fold_damping_upper 0.6\n\
             termination flanged\n\
             tract_source test-fixture\n\
             tract_licence test fixture, not data\n\
             microphone 0 0 0.3 4 16 4000\n\
             subglottal 0 0\n\
             subglottal 0.03 1400\n",
        );
        for _ in 0..8 {
            text.push_str("section 3.14e-4 0.022\n");
        }
        text.push_str(extra);
        text
    }

    fn render(text: &str) -> Vec<f64> {
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        let mut out = vec![0.0; description.samples as usize];
        for block in out.chunks_mut(4096) {
            performance.block(block).expect("renders");
        }
        out
    }

    /// The period of a steady tone: the shortest lag within 70..400 Hz whose autocorrelation reaches 95 % of the
    /// best, so a subharmonic peak cannot win. Only used as the prior of [`f0_track`].
    fn steady_f0_hz(pcm: &[f64], rate: f64) -> f64 {
        let (lo, hi) = ((rate / 400.0) as usize, (rate / 70.0) as usize);
        let n = pcm.len() - hi;
        let r: Vec<f64> = (lo..=hi)
            .map(|lag| (0..n).map(|i| pcm[i] * pcm[i + lag]).sum::<f64>())
            .collect();
        let best = r.iter().copied().fold(f64::MIN, f64::max);
        let lag = lo + r.iter().position(|&v| v >= 0.95 * best).expect("a peak");
        rate / lag as f64
    }

    /// f0 track `(time, Hz)`: autocorrelation of 40 ms windows every 10 ms, the lag searched within +-25 % of
    /// `prior_hz`, parabolic peak interpolation.
    fn f0_track(pcm: &[f64], rate: f64, prior_hz: f64) -> Vec<(f64, f64)> {
        let (win, hop) = ((0.040 * rate) as usize, (0.010 * rate) as usize);
        let lo = (rate / (1.25 * prior_hz)) as usize;
        let hi = (rate / (0.75 * prior_hz)) as usize + 1;
        let mut track = Vec::new();
        let mut start = 0;
        while start + win + hi + 2 <= pcm.len() {
            let seg = &pcm[start..start + win + hi + 2];
            let r: Vec<f64> = (lo - 1..=hi + 1)
                .map(|lag| (0..win).map(|i| seg[i] * seg[i + lag]).sum::<f64>())
                .collect();
            let k = (1..r.len() - 1)
                .max_by(|&a, &b| r[a].total_cmp(&r[b]))
                .expect("lags");
            let curvature = r[k - 1] - 2.0 * r[k] + r[k + 1];
            let shift = if curvature < 0.0 {
                0.5 * (r[k - 1] - r[k + 1]) / curvature
            } else {
                0.0
            };
            let lag = (lo - 1 + k) as f64 + shift;
            track.push(((start + win / 2) as f64 / rate, rate / lag));
            start += hop;
        }
        track
    }

    /// Amplitude of the component of `y(t)` at `freq` (mean-removed projection onto one sinusoid; exact for whole
    /// periods and an even sampling).
    fn amplitude_at(points: &[(f64, f64)], freq: f64) -> f64 {
        let n = points.len() as f64;
        let mean = points.iter().map(|p| p.1).sum::<f64>() / n;
        let (mut re, mut im) = (0.0_f64, 0.0_f64);
        for &(t, y) in points {
            let w = core::f64::consts::TAU * freq * t;
            re += (y - mean) * det::cos(w);
            im += (y - mean) * det::sin(w);
        }
        2.0 * re.hypot(im) / n
    }

    #[test]
    fn tension_modulation_parses_and_is_refused_by_name() {
        let ok = VoiceDescription::from_bytes(modulation_text("tension_modulation 5.9 0.082\n").as_bytes())
            .expect("parses");
        assert_eq!(ok.tension_modulation, Some((5.9, 0.082)));
        let none = VoiceDescription::from_bytes(modulation_text("").as_bytes()).expect("parses");
        assert_eq!(none.tension_modulation, None);
        for (extra, needle) in [
            ("tension_modulation 5.9\n", "takes 2 fields"),
            ("tension_modulation x 0.08\n", "not a number"),
            ("tension_modulation 0 0.08\n", "rate"),
            ("tension_modulation 31 0.08\n", "rate"),
            ("tension_modulation 5.9 0.6\n", "amplitude"),
            ("tension_modulation 5.9 -0.1\n", "amplitude"),
            (
                "tension_modulation 5.9 0.08\ntension_modulation 6 0.08\n",
                "only appear once",
            ),
        ] {
            let error = VoiceDescription::from_bytes(modulation_text(extra).as_bytes()).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    #[test]
    fn fold_mapping_and_fold_damping_upper_parse_and_are_refused_by_name() {
        let parse = |text: &str| VoiceDescription::from_bytes(text.as_bytes());
        let plain = parse(&fixture("")).expect("parses");
        assert_eq!(plain.fold.damping_ratio, 0.1, "the fold line's ratio is the lower mass's");
        assert_eq!(plain.fold.damping_ratio_upper, 0.6);
        // The mapping halves masses and stiffnesses, and the damping ratios stay.
        let mapped = parse(&fixture("fold_mapping bilateral\n")).expect("parses");
        assert_eq!(mapped.fold, plain.fold.bilateral());
        let upper = |line: &str| fixture("").replace("fold_damping_upper 0.6\n", line);
        for (text, needle) in [
            (fixture("fold_mapping\n"), "takes 1 fields"),
            (fixture("fold_mapping unilateral\n"), "must be `bilateral`"),
            (
                fixture("fold_mapping bilateral\nfold_mapping bilateral\n"),
                "only appear once",
            ),
            (upper(""), "needs `fold_damping_upper`"),
            (upper("fold_damping_upper\n"), "takes 1 fields"),
            (upper("fold_damping_upper x\n"), "not a number"),
            (upper("fold_damping_upper -0.6\n"), "non-negative"),
            (upper("fold_damping_upper nan\n"), "must be finite"),
            (fixture("fold_damping_upper 0.6\n"), "only appear once"),
            (
                fixture("").replace("fold two_mass", "fold one_dof"),
                "needs a two_mass fold",
            ),
        ] {
            let error = parse(&text).expect_err(needle);
            assert!(error.contains(needle), "{needle:?} -> {error}");
        }
        let one_dof = upper("").replace("fold two_mass", "fold one_dof");
        parse(&one_dof).expect("a one_dof card has one ratio");
    }

    #[test]
    fn the_upper_damping_ratio_and_the_mapping_reach_the_radiated_pressure() {
        let plain = render(&fixture(""));
        let other = fixture("").replace("fold_damping_upper 0.6\n", "fold_damping_upper 0.3\n");
        assert_ne!(plain, render(&other), "the upper ratio must reach the radiated pressure");
        let mapped = render(&fixture("fold_mapping bilateral\n"));
        assert_ne!(plain, mapped, "the mapping must reach the radiated pressure");
    }

    /// The `body_cover` row: the fold line is the cover and `fold_body` the body, both halved by the mapping; every
    /// missing, extra or misplaced part is refused by name, a body the island refuses is refused at build, and the
    /// body's spring reaches the radiated pressure.
    #[test]
    fn the_body_cover_row_parses_maps_its_body_and_is_refused_by_name() {
        let parse = |text: &str| VoiceDescription::from_bytes(text.as_bytes());
        let body_cover = |extra: &str| fixture(extra).replace("fold two_mass", "fold body_cover");
        let row = "fold_body 5.0e-5 100 0.2 3.598e-4\n";
        let body = FoldBody {
            mass_kg: 5.0e-5,
            stiffness_n_m: 100.0,
            damping_ratio: 0.2,
            rest_gap_upper_m: 3.598e-4,
        };
        let plain = parse(&body_cover(row)).expect("parses");
        assert!(plain.two_mass, "the cover is a two-mass fold");
        assert_eq!(plain.fold, parse(&fixture("")).expect("parses").fold, "the fold line reads as before");
        assert_eq!(plain.fold_body, Some(body));
        assert_eq!(parse(&fixture("")).expect("parses").fold_body, None);
        let mapped = parse(&body_cover(&format!("fold_mapping bilateral\n{row}"))).expect("parses");
        assert_eq!(mapped.fold, plain.fold.bilateral());
        assert_eq!(mapped.fold_body, Some(body.bilateral()));
        assert_eq!(body.bilateral().mass_kg, 2.5e-5);
        assert_eq!(body.bilateral().stiffness_n_m, 50.0);
        assert_eq!(body.bilateral().rest_gap_upper_m, body.rest_gap_upper_m);
        for (text, needle) in [
            (body_cover(""), "needs `fold_body`"),
            (fixture(row), "fold_body needs a body_cover fold"),
            (body_cover("fold_body 5.0e-5 100 0.2\n"), "takes 4 fields"),
            (body_cover("fold_body 5.0e-5 x 0.2 3.598e-4\n"), "not a number"),
            (body_cover(&format!("{row}{row}")), "only appear once"),
            (body_cover(&format!("{row}fold_shape sine uniform\n")), "not body_cover"),
            (
                body_cover(row).replace("fold_damping_upper 0.6\n", ""),
                "needs `fold_damping_upper`",
            ),
            (fixture("").replace("fold two_mass", "fold three_mass"), "`body_cover`"),
        ] {
            let error = parse(&text).expect_err(needle);
            assert!(error.contains(needle), "{needle:?} -> {error}");
        }
        let refused = parse(&body_cover("fold_body 5.0e-5 0 0.2 3.598e-4\n")).expect("parses");
        let error = refused.build().err().expect("a zero body spring is refused");
        assert!(error.contains("body stiffness"), "{error}");
        let stiffer = render(&body_cover("fold_body 5.0e-5 300 0.2 3.598e-4\n"));
        assert_ne!(render(&body_cover(row)), stiffer, "the body must reach the radiated pressure");
    }

    #[test]
    fn zero_amplitude_modulation_is_bit_identical_to_no_key_and_a_real_one_is_not() {
        let plain = render(&fixture(""));
        let zero = render(&fixture("tension_modulation 5.9 0\n"));
        assert!(
            plain.iter().zip(&zero).all(|(a, b)| a.to_bits() == b.to_bits()),
            "a zero-amplitude gesture must leave the tension untouched"
        );
        let swung = render(&fixture("tension_modulation 5.9 0.08\n"));
        assert!(swung.iter().all(|v| v.is_finite()));
        assert_ne!(plain, swung, "a real gesture must reach the radiated pressure");
    }

    #[test]
    fn tension_flutter_parses_and_is_refused_by_name() {
        let ok = VoiceDescription::from_bytes(modulation_text("tension_flutter 7 4 1.3 0.02\n").as_bytes())
            .expect("parses");
        let card = TensionFlutterCard { seed: 7, centre_hz: 4.0, bandwidth_hz: 1.3, sd_rel: 0.02 };
        assert_eq!(ok.tension_flutter, Some(card));
        let none = VoiceDescription::from_bytes(modulation_text("").as_bytes()).expect("parses");
        assert_eq!(none.tension_flutter, None);
        let wide = modulation_text("tension_flutter 7 4 1.3 0.1\n");
        let wide = VoiceDescription::from_bytes(wide.as_bytes()).expect("0.1 is admitted");
        assert_eq!(wide.tension_flutter.map(|card| card.sd_rel), Some(0.1));
        for (extra, needle) in [
            ("tension_flutter 7 4 1.3\n", "takes 4 fields"),
            ("tension_flutter x 4 1.3 0.02\n", "seed"),
            ("tension_flutter 7 x 1.3 0.02\n", "not a number"),
            ("tension_flutter 7 0 1.3 0.02\n", "centre"),
            ("tension_flutter 7 31 1.3 0.02\n", "centre"),
            ("tension_flutter 7 4 0 0.02\n", "bandwidth"),
            ("tension_flutter 7 4 8 0.02\n", "bandwidth"),
            ("tension_flutter 7 4 1.3 -0.01\n", "sd"),
            ("tension_flutter 7 4 1.3 0.1001\n", "sd"),
            ("tension_flutter 7 4 1.3 0.02\ntension_flutter 8 4 1.3 0.02\n", "only appear once"),
        ] {
            let error = VoiceDescription::from_bytes(modulation_text(extra).as_bytes()).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    #[test]
    fn tension_flutter_second_parses_and_is_refused_by_name() {
        let first = "tension_flutter 7 4 1.3 0.02\n";
        let parse = |extra: &str| VoiceDescription::from_bytes(modulation_text(extra).as_bytes());
        let ok = parse(&format!("{first}tension_flutter_second 7 20 20 0.005\n")).expect("parses");
        let card = TensionFlutterCard { seed: 7, centre_hz: 20.0, bandwidth_hz: 20.0, sd_rel: 0.005 };
        assert_eq!(ok.tension_flutter_second, Some(card));
        assert_eq!(parse(first).expect("parses").tension_flutter_second, None);
        for (row, needle) in [
            ("tension_flutter_second 7 20 20\n", "takes 4 fields"),
            ("tension_flutter_second x 20 20 0.005\n", "seed"),
            ("tension_flutter_second 7 x 20 0.005\n", "not a number"),
            ("tension_flutter_second 7 31 20 0.005\n", "centre"),
            ("tension_flutter_second 7 20 40 0.005\n", "bandwidth"),
            ("tension_flutter_second 7 20 20 0.1001\n", "sd"),
            ("tension_flutter_second 7 20 20 0.005\ntension_flutter_second 8 20 20 0.005\n", "only appear once"),
        ] {
            let error = parse(&format!("{first}{row}")).expect_err(row);
            assert!(error.contains(needle), "{row:?} -> {error}");
        }
        for extra in [
            "tension_flutter_second 7 20 20 0.005\n",
            "tension_modulation 5.9 0.08\ntension_flutter_second 7 20 20 0.005\n",
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains("needs tension_flutter"), "{extra:?} -> {error}");
        }
    }

    /// The row's law written out from its key, so a change to `TensionFlutter` can't move both sides of a comparison:
    /// the recursion, the ten decay times of warm-up, and the key's tile (0 is what the first row has always drawn).
    fn flutter_written_out(card: TensionFlutterCard, tile: u32, rate: u32, samples: usize) -> Vec<f64> {
        let filter = TensionFlutter::new(card, rate);
        let mut white = StreamKey { seed: card.seed, kernel: TENSION_FLUTTER_KERNEL, tile }.stream();
        let warm_up = (10.0 * f64::from(rate) / (core::f64::consts::PI * card.bandwidth_hz)).ceil() as usize;
        let (mut y1, mut y2) = (0.0_f64, 0.0_f64);
        let mut out = Vec::with_capacity(samples);
        for n in 0..warm_up + samples {
            let y = filter.a1 * y1 + filter.a2 * y2 + filter.gain * white.next_normal();
            y2 = y1;
            y1 = y;
            if n >= warm_up {
                out.push(y);
            }
        }
        out
    }

    /// A card with no flutter row, its tension set by hand before each output sample, where the row sets it.
    fn render_by_hand(text: &str, y: &[f64]) -> Vec<f64> {
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        let mut out = vec![0.0; description.samples as usize];
        for (n, sample) in out.iter_mut().enumerate() {
            performance.island.set_tension(1.0 + y[n]).expect("tension");
            performance.block(core::slice::from_mut(sample)).expect("renders");
        }
        out
    }

    /// With the second row in the engine, a card with one row or none plays what it played before the key existed; the
    /// old law is written out above. A silent second row changes nothing, and a sounding one adds to the first.
    #[test]
    fn one_flutter_row_or_none_stays_bit_identical_beside_the_second_row() {
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let plain = modulation_text("");
        let n = VoiceDescription::from_bytes(plain.as_bytes()).expect("parses").samples as usize;
        let first = TensionFlutterCard { seed: 7, centre_hz: 4.0, bandwidth_hz: 1.3, sd_rel: 0.02 };
        let second = TensionFlutterCard { seed: 7, centre_hz: 20.0, bandwidth_hz: 20.0, sd_rel: 0.005 };
        let y1 = flutter_written_out(first, 0, 48_000, n);
        let y2 = flutter_written_out(second, 1, 48_000, n);
        let one = render(&modulation_text("tension_flutter 7 4 1.3 0.02\n"));
        assert_eq!(bits(&render(&plain)), bits(&render_by_hand(&plain, &vec![0.0; n])), "no row");
        assert_eq!(bits(&one), bits(&render_by_hand(&plain, &y1)), "one row");
        let silent = render(&modulation_text("tension_flutter 7 4 1.3 0.02\ntension_flutter_second 7 20 20 0\n"));
        assert_eq!(bits(&one), bits(&silent), "a silent second row");
        let both = render(&modulation_text("tension_flutter 7 4 1.3 0.02\ntension_flutter_second 7 20 20 0.005\n"));
        let sum: Vec<f64> = y1.iter().zip(&y2).map(|(a, b)| a + b).collect();
        assert_eq!(bits(&both), bits(&render_by_hand(&plain, &sum)), "two rows add");
        assert_ne!(bits(&both), bits(&one), "the second row must reach the pressure");
    }

    /// A plain card with its tension set by hand before each output sample: `vibrato` is the old law written out, at the
    /// time the engine reads it, and `y` adds to it. With no vibrato this is `render_by_hand`.
    fn render_by_hand_with_vibrato(text: &str, vibrato: (f64, f64), y: &[f64]) -> Vec<f64> {
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        let mut out = vec![0.0; description.samples as usize];
        let (rate_hz, amplitude) = vibrato;
        for (n, sample) in out.iter_mut().enumerate() {
            let t = performance.steps as f64 / f64::from(performance.mechanics_rate);
            let v = amplitude * det::sin(core::f64::consts::TAU * rate_hz * t);
            performance.island.set_tension(1.0 + (v + y[n])).expect("tension");
            performance.block(core::slice::from_mut(sample)).expect("renders");
        }
        out
    }

    /// Vibrato alone still plays its old law. A card with both takes one tension, the sum, in either order of its lines:
    /// the row does not overwrite the vibrato, and a silent row or a silent vibrato leaves the other bit for bit.
    #[test]
    fn vibrato_and_flutter_add_on_one_tension_and_alone_stay_bit_identical() {
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let plain = modulation_text("");
        let n = VoiceDescription::from_bytes(plain.as_bytes()).expect("parses").samples as usize;
        let row = TensionFlutterCard { seed: 7, centre_hz: 20.0, bandwidth_hz: 20.0, sd_rel: 0.005 };
        let y = flutter_written_out(row, 0, 48_000, n);
        let alone = render(&modulation_text("tension_modulation 5.9 0.11\n"));
        let by_hand = render_by_hand_with_vibrato(&plain, (5.9, 0.11), &vec![0.0; n]);
        assert_eq!(bits(&alone), bits(&by_hand), "vibrato alone");
        let both = render(&modulation_text("tension_modulation 5.9 0.11\ntension_flutter 7 20 20 0.005\n"));
        let swapped = render(&modulation_text("tension_flutter 7 20 20 0.005\ntension_modulation 5.9 0.11\n"));
        assert_eq!(bits(&both), bits(&swapped), "line order");
        assert_eq!(bits(&both), bits(&render_by_hand_with_vibrato(&plain, (5.9, 0.11), &y)), "the two add");
        assert_ne!(bits(&both), bits(&alone), "the row must reach the pressure");
        let row_alone = render(&modulation_text("tension_flutter 7 20 20 0.005\n"));
        assert_ne!(bits(&both), bits(&row_alone), "the vibrato must survive the row");
        let silent_row = render(&modulation_text("tension_modulation 5.9 0.11\ntension_flutter 7 20 20 0\n"));
        assert_eq!(bits(&alone), bits(&silent_row), "a silent row");
        let silent_vibrato = render(&modulation_text("tension_modulation 5.9 0\ntension_flutter 7 20 20 0.005\n"));
        assert_eq!(bits(&row_alone), bits(&silent_vibrato), "a silent vibrato");
    }

    /// Same seed, same filter, and the two rows still don't correlate; each keeps its stated spread.
    #[test]
    fn tension_flutter_second_draws_its_own_stream() {
        let sd_rel = 0.02;
        let card = TensionFlutterCard { seed: 1, centre_hz: 20.0, bandwidth_hz: 20.0, sd_rel };
        let rate = 4800_u32;
        let mut a = TensionFlutter::on_tile(card, rate, 0);
        let mut b = TensionFlutter::on_tile(card, rate, 1);
        let n = rate as usize * 500;
        let (mut ab, mut aa, mut bb) = (0.0_f64, 0.0_f64, 0.0_f64);
        for _ in 0..n {
            let (x, y) = (a.next_relative(), b.next_relative());
            ab += x * y;
            aa += x * x;
            bb += y * y;
        }
        let r = ab / det::sqrt(aa * bb);
        assert!(r.abs() < 0.03, "correlation {r}");
        for sum in [aa, bb] {
            let sd = det::sqrt(sum / n as f64);
            assert!((sd / sd_rel - 1.0).abs() < 0.05, "sd {sd:e}");
        }
    }

    #[test]
    fn fold_shape_parses_and_is_refused_by_name() {
        let parse = |extra: &str| VoiceDescription::from_bytes(modulation_text(extra).as_bytes());
        let sine = parse("fold_shape sine uniform\n").expect("parses");
        assert_eq!(sine.fold_shape, Some(FoldShape { mode: FoldMode::Sine, rest: FoldRestGap::Uniform }));
        let rigid = parse("fold_shape rigid tapered\n").expect("parses");
        assert_eq!(rigid.fold_shape, Some(FoldShape { mode: FoldMode::Rigid, rest: FoldRestGap::Tapered }));
        assert_eq!(parse("").expect("parses").fold_shape, None);
        for (extra, needle) in [
            ("fold_shape sine\n", "takes 2 fields"),
            ("fold_shape wobble uniform\n", "mode must be"),
            ("fold_shape sine flat\n", "rest gap must be"),
            ("fold_shape sine uniform\nfold_shape sine uniform\n", "only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        let one_dof = modulation_text("fold_shape sine uniform\n")
            .replace("fold two_mass", "fold one_dof")
            .replace("fold_damping_upper 0.6\n", "");
        let error = VoiceDescription::from_bytes(one_dof.as_bytes()).expect_err("one_dof");
        assert!(error.contains("fold_shape needs a two_mass fold"), "{error}");
    }

    /// The row reaches the radiated pressure and its contact count, the partition stays free, and the tension gesture
    /// runs on a shaped fold. `min_glottal_flow_m3_s` is what a page's flow column reads.
    #[test]
    fn fold_shape_renders_and_reports_contact_and_flow() {
        let text = modulation_text("fold_shape sine uniform\n");
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let mut whole = description.build().expect("builds");
        let mut split = description.build().expect("builds");
        let (mut a, mut b) = (vec![0.0_f64; 9600], vec![0.0_f64; 9600]);
        whole.block(&mut a).expect("whole");
        for chunk in b.chunks_mut(37) {
            split.block(chunk).expect("split");
        }
        for (x, y) in a.iter().zip(&b) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        let plain = render(&modulation_text(""));
        let shaped = render(&text);
        assert!(shaped.iter().all(|v| v.is_finite()));
        assert!(plain.iter().zip(&shaped).any(|(x, y)| x != y), "the row must reach the pressure");
        let mut performance = description.build().expect("builds");
        let mut out = vec![0.0; description.samples as usize];
        for block in out.chunks_mut(4096) {
            performance.block(block).expect("renders");
        }
        let stats = performance.stats();
        assert!(stats.peak_glottal_flow_m3_s > 0.0);
        assert!(stats.min_glottal_flow_m3_s <= stats.peak_glottal_flow_m3_s);
        assert!(stats.closed_fraction > 0.0 && stats.closed_fraction < 1.0, "contact {}", stats.closed_fraction);
        let gesture = modulation_text("fold_shape sine uniform\ntension_flutter 7 4 1.3 0.02\n");
        let flutter = render(&gesture);
        assert!(flutter.iter().all(|v| v.is_finite()));
        assert!(flutter.iter().zip(&shaped).any(|(x, y)| x != y), "the gesture must reach the shaped fold");
    }

    /// The five filters of Ternstrom and Friberg's Fig. 9, at a rate low enough for a long run. Standard deviation
    /// against the closed form the gain is scaled from, mean, and the lag correlation at half the damped period,
    /// `-exp(-pi B tau)`, which a filter of the wrong centre or bandwidth would miss.
    #[test]
    fn tension_flutter_matches_its_closed_form_on_the_printed_filters() {
        for (centre_hz, bandwidth_hz) in [(6.0, 1.5), (6.0, 3.0), (4.0, 1.3), (4.0, 5.0), (3.0, 3.0)] {
            let sd_rel = 0.02;
            let card = TensionFlutterCard { seed: 1, centre_hz, bandwidth_hz, sd_rel };
            let rate = 4800_u32;
            let mut flutter = TensionFlutter::new(card, rate);
            let y: Vec<f64> = (0..rate as usize * 1000).map(|_| flutter.next_relative()).collect();
            let mean = y.iter().sum::<f64>() / y.len() as f64;
            let sd = det::sqrt(y.iter().map(|v| v * v).sum::<f64>() / y.len() as f64);
            assert!(mean.abs() < 0.1 * sd_rel, "({centre_hz}, {bandwidth_hz}): mean {mean:e}");
            assert!((sd / sd_rel - 1.0).abs() < 0.05, "({centre_hz}, {bandwidth_hz}): sd {sd:e}");
            let half_bandwidth = bandwidth_hz / 2.0;
            let tau = 0.5 / det::sqrt(centre_hz * centre_hz - half_bandwidth * half_bandwidth);
            let lag = (tau * f64::from(rate)).round() as usize;
            let covariance = y.iter().zip(&y[lag..]).map(|(a, b)| a * b).sum::<f64>() / (y.len() - lag) as f64;
            let read = covariance / (sd * sd);
            let closed = -det::exp(-core::f64::consts::PI * bandwidth_hz * tau);
            assert!(
                (read - closed).abs() < 0.05,
                "({centre_hz}, {bandwidth_hz}): lag correlation {read} against {closed}"
            );
        }
    }

    /// At the audio rate the row must move smoothly. For (4, 1.3) at 48 kHz the closed-form rms step is `sqrt(2 (1 -
    /// rho_1))` = 5.2e-4 of the extent, so a limit of 5e-3 is about ten times that: a process that stepped once per
    /// cycle would pass it by far. The standard deviation holds at that rate too, and from the first sample: over many
    /// seeds the first value already has the full spread.
    #[test]
    fn tension_flutter_is_smooth_at_the_audio_rate_and_starts_stationary() {
        let sd_rel = 0.02;
        let card = TensionFlutterCard { seed: 1, centre_hz: 4.0, bandwidth_hz: 1.3, sd_rel };
        let mut flutter = TensionFlutter::new(card, 48_000);
        let y: Vec<f64> = (0..48_000 * 200).map(|_| flutter.next_relative()).collect();
        let sd = det::sqrt(y.iter().map(|v| v * v).sum::<f64>() / y.len() as f64);
        assert!((sd / sd_rel - 1.0).abs() < 0.1, "sd {sd:e} at 48 kHz");
        let largest = y.windows(2).map(|w| (w[1] - w[0]).abs()).fold(0.0, f64::max);
        assert!(largest < 5e-3 * sd_rel, "a step of {largest:e} against an extent of {sd_rel:e}");
        let first: Vec<f64> = (0..400)
            .map(|seed| TensionFlutter::new(TensionFlutterCard { seed, ..card }, 4800).next_relative())
            .collect();
        let spread = det::sqrt(first.iter().map(|v| v * v).sum::<f64>() / first.len() as f64);
        assert!((spread / sd_rel - 1.0).abs() < 0.15, "first-sample spread {spread:e}");
    }

    #[test]
    fn zero_sd_flutter_is_bit_identical_to_no_key_and_a_real_one_is_seeded() {
        let plain = render(&fixture(""));
        let zero = render(&fixture("tension_flutter 7 4 1.3 0\n"));
        assert!(first_difference(&plain, &zero).is_none(), "a zero extent must leave the tension untouched");
        let moved = render(&fixture("tension_flutter 7 4 1.3 0.03\n"));
        assert!(moved.iter().all(|v| v.is_finite()));
        assert!(first_difference(&plain, &moved).is_some(), "a real row must reach the radiated pressure");
        assert!(first_difference(&moved, &render(&fixture("tension_flutter 7 4 1.3 0.03\n"))).is_none());
        assert!(first_difference(&moved, &render(&fixture("tension_flutter 8 4 1.3 0.03\n"))).is_some());
    }

    #[test]
    fn tension_flutter_rendering_is_independent_of_callback_partition() {
        let text = fixture("tension_flutter 7 4 1.3 0.03\n");
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let mut whole = description.build().expect("builds");
        let mut split = description.build().expect("builds");
        let mut a = vec![0.0; 600];
        whole.block(&mut a).expect("whole block");
        let mut b = vec![0.0; 600];
        for chunk in b.chunks_mut(37) {
            split.block(chunk).expect("split block");
        }
        assert!(first_difference(&a, &b).is_none(), "block partition changed the pressure");
    }

    /// With a planted negative: the radiated f0 reads steady with no key and moves with the row, and doubling `sd_rel`
    /// on one seed about doubles the movement, so the noise reaches the pitch. Measured on the fixture's own tube: the
    /// result is no note's extent. The broad filter (4, 5) puts several independent values in this 1 s hold.
    #[test]
    fn tension_flutter_moves_the_radiated_pitch_in_proportion_to_its_sd() {
        let rate = 48_000.0_f64;
        let steady = render(&modulation_text(""));
        let prior = steady_f0_hz(&steady[(0.4 * rate) as usize..(0.9 * rate) as usize], rate);
        let sd_cents = |pcm: &[f64]| -> f64 {
            let track = f0_track(&pcm[(0.30 * rate) as usize..], rate, prior);
            let count = track.len() as f64;
            let mean = track.iter().map(|p| p.1).sum::<f64>() / count;
            let cents = track.iter().map(|p| 1200.0 * det::ln(p.1 / mean) / core::f64::consts::LN_2);
            det::sqrt(cents.map(|c| c * c).sum::<f64>() / count)
        };
        let none = sd_cents(&steady);
        let small = sd_cents(&render(&modulation_text("tension_flutter 7 4 5 0.01\n")));
        let full = sd_cents(&render(&modulation_text("tension_flutter 7 4 5 0.02\n")));
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-tension-flutter\",\"prior_hz\":{prior:.2},\"sd_cents_none\":{none:.2},\"sd_cents_0p01\":{small:.1},\"sd_cents_0p02\":{full:.1}}}"
        );
        assert!(none < 0.2 * full, "no key, yet {none:.2} cents against {full:.1}");
        assert!((1.6..=2.4).contains(&(full / small)), "doubling the sd gave {small:.1} -> {full:.1} cents");
    }

    fn first_difference(a: &[f64], b: &[f64]) -> Option<usize> {
        a.iter().zip(b).position(|(x, y)| x.to_bits() != y.to_bits())
    }

    #[test]
    fn aspiration_noise_parses_and_is_refused_by_name() {
        let parse = |extra: &str| VoiceDescription::from_bytes(fixture(extra).as_bytes());
        let ok = parse("aspiration_noise 7 1e-4\n").expect("parses");
        assert_eq!(ok.aspiration_noise, Some(AspirationCard { seed: 7, xi: 1.0e-4, re_crit: 1700.0 }));
        let own = parse("aspiration_noise 7 1e-4 2500\n").expect("parses");
        assert_eq!(own.aspiration_noise.map(|card| card.re_crit), Some(2500.0));
        assert_eq!(parse("").expect("parses").aspiration_noise, None);
        for (extra, needle) in [
            ("aspiration_noise\n", "takes 2 or 3 fields"),
            ("aspiration_noise 7\n", "xi has no default"),
            ("aspiration_noise 7 1e-4 1700 9\n", "takes 2 or 3 fields"),
            ("aspiration_noise x 1e-4\n", "not a non-negative integer"),
            ("aspiration_noise 7 x\n", "not a number"),
            ("aspiration_noise 7 nan\n", "must be finite"),
            ("aspiration_noise 7 0\n", "xi must be finite and positive"),
            ("aspiration_noise 7 -1e-4\n", "xi must be finite and positive"),
            ("aspiration_noise 7 1e-4 -1\n", "re_crit"),
            ("aspiration_noise 7 1e-4\naspiration_noise 8 1e-4\n", "only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    #[test]
    fn aspiration_is_silent_below_the_onset_flow_noisy_above_it_and_seeded() {
        let plain = render(&fixture(""));
        // A threshold no glottal jet reaches: the gate never opens, nothing is added, and every bit stays.
        let shut = render(&fixture("aspiration_noise 7 1e-4 1e12\n"));
        assert!(
            plain.iter().zip(&shut).all(|(a, b)| a.to_bits() == b.to_bits()),
            "a gate that never opens must leave the output bit-identical to no key"
        );
        // The default threshold (1700): the fixture's jet reaches Re of order 1e4, so the gate opens and the noise reaches the output.
        let gated = render(&fixture("aspiration_noise 7 1e-4\n"));
        assert!(gated.iter().all(|v| v.is_finite()));
        let first_gated = first_difference(&plain, &gated).expect("noise above the onset flow reaches the output");
        // Planted negative: a gate with no threshold (Re_c = 0) opens at the first flow instead of at Re_c, so its first
        // difference comes far earlier. The test above would pass a source that ignored its threshold; this one would not.
        let ungated = render(&fixture("aspiration_noise 7 1e-4 0\n"));
        let first_ungated = first_difference(&plain, &ungated).expect("an ungated source speaks");
        assert!(
            first_ungated + 100 < first_gated,
            "the gated source must stay silent while Re < Re_c: ungated first difference at {first_ungated}, gated at {first_gated}"
        );
        // A declared seed gives one result, another seed another.
        let again = render(&fixture("aspiration_noise 7 1e-4\n"));
        assert!(gated.iter().zip(&again).all(|(a, b)| a.to_bits() == b.to_bits()), "the same seed must repeat bit for bit");
        let other = render(&fixture("aspiration_noise 8 1e-4\n"));
        assert_ne!(gated, other, "another seed must give another noise");
        // The noise is broadband and the gain scales it: a bigger gain moves the output further from the noiseless render.
        let rms_from_plain = |x: &[f64]| (x.iter().zip(&plain).map(|(a, b)| (a - b) * (a - b)).sum::<f64>() / x.len() as f64).sqrt();
        let loud = render(&fixture("aspiration_noise 7 4e-4\n"));
        assert!(rms_from_plain(&loud) > 1.5 * rms_from_plain(&gated), "a 4x gain must move the output further");
    }

    #[test]
    fn jet_noise_parses_and_is_refused_by_name() {
        let parse = |extra: &str| VoiceDescription::from_bytes(fixture(extra).as_bytes());
        let ok = parse("jet_noise 7 1.3e-3 0.045 0.115 1700\n").expect("parses");
        assert_eq!(
            ok.jet_noise,
            Some(JetNoiseCard { seed: 7, gain: 1.3e-3, strouhal_lo: 0.045, strouhal_hi: 0.115, reynolds_onset: 1700.0 })
        );
        assert_eq!(parse("").expect("parses").jet_noise, None);
        for (extra, needle) in [
            ("jet_noise\n", "takes 5 fields"),
            ("jet_noise 7 1.3e-3 0.045 0.115\n", "takes 5 fields"),
            ("jet_noise 7 1.3e-3 0.045 0.115 1700 9\n", "takes 5 fields"),
            ("jet_noise x 1.3e-3 0.045 0.115 1700\n", "not a non-negative integer"),
            ("jet_noise 7 x 0.045 0.115 1700\n", "not a number"),
            ("jet_noise 7 nan 0.045 0.115 1700\n", "must be finite"),
            ("jet_noise 7 0 0.045 0.115 1700\n", "gain must be finite and positive"),
            ("jet_noise 7 -1e-3 0.045 0.115 1700\n", "gain must be finite and positive"),
            ("jet_noise 7 1.3e-3 0.115 0.045 1700\n", "Strouhal band"),
            ("jet_noise 7 1.3e-3 0 0.115 1700\n", "Strouhal band"),
            ("jet_noise 7 1.3e-3 0.045 0.115 -1\n", "reynolds onset"),
            ("jet_noise 7 1.3e-3 0.045 0.115 1700\njet_noise 8 1.3e-3 0.045 0.115 1700\n", "only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    /// Gauntlet G0 (a card without the row, or with a gate that never opens, is bit-identical to one written before the
    /// row existed), G3 (the gate, the seed, the gain) and a planted negative (a source that ignored its gate would speak
    /// from the first flow).
    #[test]
    fn jet_noise_is_silent_when_its_gate_never_opens_and_otherwise_gated_seeded_and_scaled_by_its_gain() {
        let plain = render(&fixture(""));
        let shut = render(&fixture("jet_noise 7 1.3e-3 0.045 0.115 1e12\n"));
        assert!(
            plain.iter().zip(&shut).all(|(a, b)| a.to_bits() == b.to_bits()),
            "a gate that never opens must leave the output bit-identical to no row"
        );
        let gated = render(&fixture("jet_noise 7 1.3e-3 0.045 0.115 300\n"));
        assert!(gated.iter().all(|v| v.is_finite()));
        let first_gated = first_difference(&plain, &gated).expect("a jet above the onset flow reaches the output");
        // Planted negative: an onset of zero opens at the first flow, so its first difference comes far earlier.
        let ungated = render(&fixture("jet_noise 7 1.3e-3 0.045 0.115 0\n"));
        let first_ungated = first_difference(&plain, &ungated).expect("an ungated source speaks");
        assert!(
            first_ungated + 50 < first_gated,
            "the gated source must stay silent below its onset flow: ungated first difference at {first_ungated}, gated at {first_gated}"
        );
        let again = render(&fixture("jet_noise 7 1.3e-3 0.045 0.115 300\n"));
        assert!(gated.iter().zip(&again).all(|(a, b)| a.to_bits() == b.to_bits()), "the same seed must repeat bit for bit");
        assert_ne!(gated, render(&fixture("jet_noise 8 1.3e-3 0.045 0.115 300\n")), "another seed must give another noise");
        let rms_from_plain = |x: &[f64]| (x.iter().zip(&plain).map(|(a, b)| (a - b) * (a - b)).sum::<f64>() / x.len() as f64).sqrt();
        let loud = render(&fixture("jet_noise 7 5.2e-3 0.045 0.115 300\n"));
        assert!(rms_from_plain(&loud) > 1.5 * rms_from_plain(&gated), "a 4x gain must move the output further");
    }

    #[test]
    fn jet_noise_and_aspiration_noise_stand_on_one_card_with_streams_of_their_own() {
        let both = render(&fixture("aspiration_noise 7 1e-4\njet_noise 7 1.3e-3 0.045 0.115 300\n"));
        assert!(both.iter().all(|v| v.is_finite()));
        let aspiration = render(&fixture("aspiration_noise 7 1e-4\n"));
        let jet = render(&fixture("jet_noise 7 1.3e-3 0.045 0.115 300\n"));
        assert_ne!(both, aspiration, "the jet row must add to a card that has aspiration");
        assert_ne!(both, jet, "the aspiration row must add to a card that has the jet row");
    }

    /// The head is the baffled pressure through the sphere's FIR and nothing else: the tract, the fold and the
    /// receiver are the plain card's, and the head's latency joins the reported filter latency.
    #[test]
    fn head_sphere_filters_the_baffled_pressure_and_nothing_else() {
        let plain = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        let headed = VoiceDescription::from_bytes(fixture("head_sphere 0.0875 brown-duda-1998\n").as_bytes()).unwrap();
        assert_eq!(headed.head_sphere, Some(HeadSphereCard { radius_m: 0.0875, source: "brown-duda-1998".into() }));
        let (a, b) = (plain.build().unwrap(), headed.build().unwrap());
        assert_eq!(b.filter_delay_samples(), a.filter_delay_samples() + 64.0);
        assert_eq!(a.info().head_sphere_radius_m, None);
        assert_eq!(b.info().head_sphere_radius_m, Some(0.0875));
        let gas = plain.ambient_gas();
        let medium = RayleighMedium { density: gas.density, sound_speed: gas.sound_speed };
        let mut head = HeadSphere::new(0.0875, a.info().mouth_radius_m, plain.microphone, medium, 48_000).unwrap();
        let by_hand: Vec<f64> = pressure(&plain).into_iter().map(|p| head.step(p)).collect();
        assert_eq!(pressure(&headed), by_hand);
    }

    #[test]
    fn head_sphere_rows_are_refused_by_name() {
        for (row, why) in [
            ("head_sphere 0.0875\n", "takes 2 fields"),
            ("head_sphere 0 src\n", "must be positive"),
            ("head_sphere -0.1 src\n", "must be positive"),
            ("head_sphere nan src\n", "must be finite"),
            ("head_sphere 0.0875 a\nhead_sphere 0.0875 b\n", "only appear once"),
        ] {
            let err = VoiceDescription::from_bytes(fixture(row).as_bytes()).unwrap_err();
            assert!(err.contains(why), "{row}: {err}");
        }
        // A head no larger than the mouth (radius 8 mm here) is refused when the performance is built.
        let small = VoiceDescription::from_bytes(fixture("head_sphere 0.005 src\n").as_bytes()).unwrap();
        let err = small.build().err().expect("refused");
        assert!(err.contains("smaller than the head"), "{err}");
    }

    #[test]
    fn head_load_rows_are_refused_by_name() {
        for (row, why) in [
            ("head_load cap\n", "needs the head_sphere row"),
            ("head_sphere 0.0875 src\nhead_load\n", "takes 1 fields"),
            ("head_sphere 0.0875 src\nhead_load disk\n", "must be `cap`"),
            ("head_sphere 0.0875 src\nhead_load cap\nhead_load cap\n", "only appear once"),
        ] {
            let err = VoiceDescription::from_bytes(fixture(row).as_bytes()).unwrap_err();
            assert!(err.contains(why), "{row}: {err}");
        }
        let unflanged = fixture("head_sphere 0.0875 src\nhead_load cap\n").replace("termination flanged", "termination unflanged");
        let err = VoiceDescription::from_bytes(unflanged.as_bytes()).unwrap_err();
        assert!(err.contains("termination flanged"), "{err}");
        // The fields are public, so `build` refuses the same combinations again.
        let mut bare = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        bare.head_load = true;
        let err = bare.build().err().expect("refused");
        assert!(err.contains("needs the head_sphere row"), "{err}");
    }

    /// The head's load at the lips: the island's line, the mouth flow and the passive formants read it. The end
    /// correction falls (mail 7212's prediction), so the first two formants rise by under 2 %, and the pressure moves.
    #[test]
    fn head_load_moves_the_lips_load() {
        let headed = VoiceDescription::from_bytes(fixture("head_sphere 0.0875 brown-duda-1998\n").as_bytes()).unwrap();
        let loaded =
            VoiceDescription::from_bytes(fixture("head_sphere 0.0875 brown-duda-1998\nhead_load cap\n").as_bytes()).unwrap();
        assert!(!headed.head_load && loaded.head_load);
        let (a, b) = (headed.build().unwrap(), loaded.build().unwrap());
        assert!(!a.info().head_load && b.info().head_load);
        assert_eq!(b.filter_delay_samples(), a.filter_delay_samples());
        let (fa, fb) = (a.info().passive_formants_hz.clone().unwrap(), b.info().passive_formants_hz.clone().unwrap());
        assert_eq!(fa.len(), fb.len());
        for (x, y) in fa.iter().zip(&fb).take(2) {
            assert!(*y > *x && *y < 1.02 * *x, "formant {x} -> {y} Hz");
        }
        for (x, y) in fa.iter().zip(&fb) {
            assert!((y / x - 1.0).abs() < 0.02, "formant {x} -> {y} Hz");
        }
        assert_ne!(pressure(&headed), pressure(&loaded));
    }

    const SNOWMAN: &str = "head_sphere 0.087 algazi-duda-thompson-2002\ntorso_sphere 0.169 0.053 algazi-duda-thompson-2002\n";

    /// The torso without `head_load` moves the microphone and nothing else: the pressure is the plain card's baffled
    /// pressure through the snowman's FIR, by hand. With `head_load` the lips' load gains the torso's change too, so the
    /// passive formants move from the head's own, by well under 1 % (the change is 0.2-0.8 % of the cap's load).
    #[test]
    fn torso_sphere_moves_the_microphone_and_with_head_load_the_lips() {
        let plain = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        let headed = VoiceDescription::from_bytes(fixture("head_sphere 0.087 algazi-duda-thompson-2002\n").as_bytes()).unwrap();
        let snow = VoiceDescription::from_bytes(fixture(SNOWMAN).as_bytes()).unwrap();
        let kemar = TorsoSphere { radius_m: 0.169, neck_gap_m: 0.053 };
        assert_eq!(snow.torso_sphere, Some(TorsoSphereCard { torso: kemar, source: "algazi-duda-thompson-2002".into() }));
        assert_eq!(headed.torso_sphere, None);
        let (a, b) = (headed.build().unwrap(), snow.build().unwrap());
        assert_eq!(b.filter_delay_samples(), a.filter_delay_samples());
        assert_eq!((a.info().torso_sphere, b.info().torso_sphere), (None, Some(kemar)));
        assert_eq!(a.info().passive_formants_hz, b.info().passive_formants_hz, "no head_load: the lips' load is the card's");
        let gas = plain.ambient_gas();
        let medium = RayleighMedium { density: gas.density, sound_speed: gas.sound_speed };
        let mouth = plain.build().unwrap().info().mouth_radius_m;
        let mut head = HeadSphere::with_torso(0.087, mouth, plain.microphone, medium, 48_000, kemar).unwrap();
        let by_hand: Vec<f64> = pressure(&plain).into_iter().map(|p| head.step(p)).collect();
        assert_eq!(pressure(&snow), by_hand);
        assert_ne!(pressure(&snow), pressure(&headed));
        let head_loaded =
            VoiceDescription::from_bytes(fixture("head_sphere 0.087 algazi-duda-thompson-2002\nhead_load cap\n").as_bytes()).unwrap();
        let snow_loaded = VoiceDescription::from_bytes(fixture(&format!("{SNOWMAN}head_load cap\n")).as_bytes()).unwrap();
        let (c, d) = (head_loaded.build().unwrap(), snow_loaded.build().unwrap());
        assert!(c.info().head_load && d.info().head_load);
        let (fc, fd) = (c.info().passive_formants_hz.clone().unwrap(), d.info().passive_formants_hz.clone().unwrap());
        assert_eq!(fc.len(), fd.len());
        assert_ne!(fc, fd, "the torso must reach the lips' load");
        for (x, y) in fc.iter().zip(&fd) {
            assert!((y / x - 1.0).abs() < 0.01, "formant {x} -> {y} Hz");
        }
        assert_ne!(pressure(&head_loaded), pressure(&snow_loaded));
    }

    #[test]
    fn torso_sphere_rows_are_refused_by_name() {
        for (row, why) in [
            ("torso_sphere 0.169 0.053 src\n", "needs the head_sphere row"),
            ("head_sphere 0.087 s\ntorso_sphere 0.169 0.053\n", "takes 3 fields"),
            ("head_sphere 0.087 s\ntorso_sphere 0 0.053 src\n", "torso radius must be finite and positive"),
            ("head_sphere 0.087 s\ntorso_sphere 0.169 0 src\n", "must not touch"),
            ("head_sphere 0.087 s\ntorso_sphere 0.169 nan src\n", "neck gap must be finite"),
            ("head_sphere 0.087 s\ntorso_sphere 0.169 0.053 a\ntorso_sphere 0.169 0.053 b\n", "only appear once"),
        ] {
            let err = VoiceDescription::from_bytes(fixture(row).as_bytes()).unwrap_err();
            assert!(err.contains(why), "{row}: {err}");
        }
        // The fields are public, so `build` refuses a torso without a head again.
        let mut bare = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        bare.torso_sphere = Some(TorsoSphereCard { torso: TorsoSphere { radius_m: 0.169, neck_gap_m: 0.053 }, source: "s".into() });
        let err = bare.build().err().expect("refused");
        assert!(err.contains("needs the head_sphere row"), "{err}");
    }

    const HEAD_LOAD: &str = "head_sphere 0.087 s\nhead_load cap\n";

    /// A `tract_gas` row holding the exterior's own record takes the exterior-checked lines with one gas, which are the
    /// old lines exactly: the receipts, the formants and the pressure are bit for bit the card's without the row.
    #[test]
    fn tract_gas_equal_to_ambient_is_the_card_without_the_row() {
        let base = VoiceDescription::from_bytes(fixture(HEAD_LOAD).as_bytes()).unwrap();
        let same = VoiceDescription::from_bytes(
            fixture(&format!("{HEAD_LOAD}tract_gas moist_air 310.15 101325 1.0\n")).as_bytes(),
        )
        .unwrap();
        assert!(base.tract_gas_record().is_none() && same.tract_gas_record().is_some());
        assert_eq!(same.tract_gas(), same.ambient_gas());
        assert_eq!(base.tract_gas(), base.ambient_gas());
        let (a, b) = (base.build().unwrap(), same.build().unwrap());
        assert_eq!(a.info(), b.info());
        assert_eq!(pressure(&base), pressure(&same));
    }

    /// A tract gas apart from the exterior: the tract's lines and formants take its sound speed, so each passive formant
    /// moves by the ratio of the two sound speeds to within 1 % (the lips' load stays the exterior's absolute impedance,
    /// so the end correction seen from the tract scales with the density ratio instead), and the receiver and the head
    /// stay the exterior's.
    #[test]
    fn tract_gas_moves_the_formants_by_the_sound_speed_ratio() {
        let base = VoiceDescription::from_bytes(fixture(HEAD_LOAD).as_bytes()).unwrap();
        let breath = VoiceDescription::from_bytes(
            fixture(&format!("{HEAD_LOAD}tract_gas moist_air 290.15 101325 0.5\n")).as_bytes(),
        )
        .unwrap();
        assert_eq!(breath.ambient_gas(), base.ambient_gas());
        let ratio = breath.tract_gas().sound_speed / base.ambient_gas().sound_speed;
        assert!(ratio < 0.97, "the test gas must move the sound speed: {ratio}");
        let (a, b) = (base.build().unwrap(), breath.build().unwrap());
        assert_eq!(b.info().sound_speed_m_s, breath.tract_gas().sound_speed);
        assert_eq!(b.filter_delay_samples(), a.filter_delay_samples(), "the exterior's head is unchanged");
        let (fa, fb) = (a.info().passive_formants_hz.clone().unwrap(), b.info().passive_formants_hz.clone().unwrap());
        assert_eq!(fa.len(), fb.len());
        for (x, y) in fa.iter().zip(&fb) {
            assert!((y / x / ratio - 1.0).abs() < 0.01, "formant {x} -> {y} Hz against the ratio {ratio}");
        }
        assert_ne!(pressure(&base), pressure(&breath));
    }

    #[test]
    fn tract_gas_rows_are_refused_by_name() {
        let rows: [(String, &str); 5] = [
            ("tract_gas moist_air 310.15 101325 1.0\n".into(), "needs the head_load row"),
            ("head_sphere 0.087 s\ntract_gas moist_air 310.15 101325 1.0\n".into(), "needs the head_load row"),
            (format!("{HEAD_LOAD}tract_gas moist_air 310.15 90000 1.0\n"), "one static pressure"),
            (format!("{HEAD_LOAD}tract_gas moist_air 310.15 101325 1.0 7\n"), "contains extra fields"),
            (
                format!("{HEAD_LOAD}tract_gas moist_air 310.15 101325 1.0\ntract_gas moist_air 300 101325 1.0\n"),
                "only appear once",
            ),
        ];
        for (row, why) in rows {
            let err = VoiceDescription::from_bytes(fixture(&row).as_bytes()).unwrap_err();
            assert!(err.contains(why), "{row}: {err}");
        }
        // The fields are public, so `build` refuses a tract gas without the head's load again.
        let mut flipped = VoiceDescription::from_bytes(
            fixture(&format!("{HEAD_LOAD}tract_gas moist_air 290.15 101325 0.5\n")).as_bytes(),
        )
        .unwrap();
        flipped.head_load = false;
        let err = flipped.build().err().expect("refused");
        assert!(err.contains("needs the head_load row"), "{err}");
    }

    #[test]
    fn pulse_shimmer_parses_and_is_refused_by_name() {
        let ok = VoiceDescription::from_bytes(modulation_text("pulse_shimmer 7 0.0113\n").as_bytes()).expect("parses");
        assert_eq!(ok.pulse_shimmer, Some(PulseShimmerCard { seed: 7, sd_rel: 0.0113 }));
        let none = VoiceDescription::from_bytes(modulation_text("").as_bytes()).expect("parses");
        assert_eq!(none.pulse_shimmer, None);
        for (extra, needle) in [
            ("pulse_shimmer 7\n", "takes 2 fields"),
            ("pulse_shimmer x 0.01\n", "seed"),
            ("pulse_shimmer 7 x\n", "not a number"),
            ("pulse_shimmer 7 -0.01\n", "sd"),
            ("pulse_shimmer 7 0.06\n", "sd"),
            ("pulse_shimmer 7 0.01\npulse_shimmer 8 0.01\n", "only appear once"),
        ] {
            let error = VoiceDescription::from_bytes(modulation_text(extra).as_bytes()).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    /// The row on its own: no draw without an opening (a step open after a shut one), one draw per opening, and the
    /// gains have mean 1 and the stated spread.
    #[test]
    fn pulse_shimmer_draws_once_per_opening_with_its_stated_spread() {
        let card = PulseShimmerCard { seed: 3, sd_rel: 0.02 };
        let mut row = PulseShimmer::new(card);
        for shut in [false, false, false] {
            row.observe(shut);
        }
        assert_eq!((row.draws, row.gain.to_bits()), (0, 1.0_f64.to_bits()), "open from the start draws nothing");
        for _ in 0..5 {
            row.observe(true);
        }
        assert_eq!(row.draws, 0, "a closing draws nothing");
        let mut gains = Vec::new();
        for _ in 0..20_000 {
            row.observe(false);
            gains.push(row.gain);
            row.observe(false);
            assert_eq!(row.gain.to_bits(), gains[gains.len() - 1].to_bits(), "the gain holds through the cycle");
            row.observe(true);
        }
        assert_eq!(row.draws, 20_000);
        let mean = gains.iter().sum::<f64>() / gains.len() as f64;
        let sd = det::sqrt(gains.iter().map(|g| (g - mean) * (g - mean)).sum::<f64>() / gains.len() as f64);
        assert!((mean - 1.0).abs() < 0.05 * card.sd_rel, "mean {mean}");
        assert!((sd / card.sd_rel - 1.0).abs() < 0.03, "sd {sd:e}");
        let lag = gains.iter().zip(&gains[1..]).map(|(a, b)| (a - mean) * (b - mean)).sum::<f64>() / (gains.len() - 1) as f64;
        assert!((lag / (sd * sd)).abs() < 0.03, "consecutive gains correlate {}", lag / (sd * sd));
    }

    /// Gauntlet G0: no key or `sd_rel` 0 leaves the pressure bit for bit. A real gain reaches the pressure, repeats on its
    /// seed, draws about once per cycle, keeps the callback partition free, and leaves the island alone: the island's
    /// incident wave and gap are those of the card without the key at every output sample.
    #[test]
    fn pulse_shimmer_is_bit_identical_when_silent_and_otherwise_leaves_the_island_alone() {
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let plain = render(&modulation_text(""));
        assert_eq!(bits(&plain), bits(&render(&modulation_text("pulse_shimmer 7 0\n"))), "sd 0");
        let text = modulation_text("pulse_shimmer 7 0.0113\n");
        let shimmer = render(&text);
        assert!(shimmer.iter().all(|v| v.is_finite()));
        assert_ne!(bits(&plain), bits(&shimmer), "a real gain must reach the pressure");
        assert_eq!(bits(&shimmer), bits(&render(&text)), "the same seed repeats");
        assert_ne!(bits(&shimmer), bits(&render(&modulation_text("pulse_shimmer 8 0.0113\n"))), "another seed");
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let mut split = description.build().expect("builds");
        let mut chunked = vec![0.0; description.samples as usize];
        for chunk in chunked.chunks_mut(37) {
            split.block(chunk).expect("split");
        }
        assert_eq!(bits(&shimmer), bits(&chunked), "block partition changed the pressure");
        let build = |t: &str| VoiceDescription::from_bytes(t.as_bytes()).expect("parses").build().expect("builds");
        let (mut a, mut b) = (build(&modulation_text("")), build(&text));
        let (mut x, mut y) = (0.0_f64, 0.0_f64);
        for n in 0..description.samples {
            a.block(core::slice::from_mut(&mut x)).expect("renders");
            b.block(core::slice::from_mut(&mut y)).expect("renders");
            assert_eq!(a.island.incident_wave_pa().to_bits(), b.island.incident_wave_pa().to_bits(), "wave at {n}");
            assert_eq!(a.island.gap_m().to_bits(), b.island.gap_m().to_bits(), "gap at {n}");
        }
        assert_eq!(a.closed_steps, b.closed_steps);
        let rate = 48_000.0_f64;
        let prior = steady_f0_hz(&plain[(0.4 * rate) as usize..(0.9 * rate) as usize], rate);
        assert_eq!(a.stats().pulse_shimmer_draws, 0, "no key, no draws");
        let draws = b.stats().pulse_shimmer_draws as f64;
        let seconds = description.samples as f64 / rate;
        assert!(
            draws > 0.6 * prior * seconds && draws <= 1.05 * prior * seconds,
            "{draws} draws in {seconds} s at {prior:.1} Hz"
        );
    }

    #[test]
    fn constriction_noise_parses_and_is_refused_by_name() {
        let parse = |extra: &str| VoiceDescription::from_bytes(fixture(extra).as_bytes());
        let ok = parse("constriction_noise 7 5.3e-3\n").expect("parses");
        assert_eq!(ok.constriction_noise, Some(ConstrictionNoiseCard { seed: 7, duct_radius_m: 5.3e-3 }));
        assert_eq!(parse("").expect("parses").constriction_noise, None);
        for (extra, needle) in [
            ("constriction_noise\n", "takes 2 fields"),
            ("constriction_noise 7\n", "takes 2 fields"),
            ("constriction_noise 7 5.3e-3 9\n", "takes 2 fields"),
            ("constriction_noise x 5.3e-3\n", "not a non-negative integer"),
            ("constriction_noise 7 x\n", "not a number"),
            ("constriction_noise 7 nan\n", "must be finite"),
            ("constriction_noise 7 0\n", "duct radius must be finite and positive"),
            ("constriction_noise 7 -5.3e-3\n", "duct radius must be finite and positive"),
            ("constriction_noise 7 5.3e-3\nconstriction_noise 8 5.3e-3\n", "only appear once"),
        ] {
            let error = parse(extra).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
    }

    /// Gauntlet G0 (a cutoff under every driven band leaves nothing to add: bit-identical to no row), G3 (the printed law
    /// reaches the radiated sound, finite and seeded) and the planted negative that keeps the first from being vacuous (the
    /// same row at a real radius speaks).
    #[test]
    fn constriction_noise_is_bit_identical_with_no_band_below_its_cutoff_and_otherwise_seeded_and_audible() {
        let plain = render(&fixture(""));
        // A 12 m "duct" has its cutoff `100 / r` at 8.3 Hz, under the bank's lowest band (10 Hz), so nothing is driven.
        let cut_off = render(&fixture("constriction_noise 7 12\n"));
        assert!(
            plain.iter().zip(&cut_off).all(|(a, b)| a.to_bits() == b.to_bits()),
            "a row that drives no band must leave the output bit-identical to no row"
        );
        let spoken = render(&fixture("constriction_noise 7 5.3e-3\n"));
        assert!(spoken.iter().all(|v| v.is_finite()));
        assert!(first_difference(&plain, &spoken).is_some(), "the row at a real radius must reach the output");
        let again = render(&fixture("constriction_noise 7 5.3e-3\n"));
        assert!(spoken.iter().zip(&again).all(|(a, b)| a.to_bits() == b.to_bits()), "the same seed must repeat bit for bit");
        assert_ne!(spoken, render(&fixture("constriction_noise 8 5.3e-3\n")), "another seed must give another noise");
    }

    #[test]
    fn constriction_noise_and_jet_noise_stand_on_one_card_with_streams_of_their_own() {
        let both = render(&fixture("jet_noise 7 1.3e-3 0.045 0.115 300\nconstriction_noise 7 5.3e-3\n"));
        assert!(both.iter().all(|v| v.is_finite()));
        assert_ne!(both, render(&fixture("jet_noise 7 1.3e-3 0.045 0.115 300\n")), "the row must add to a card that has jet_noise");
        assert_ne!(both, render(&fixture("constriction_noise 7 5.3e-3\n")), "jet_noise must add to a card that has the row");
    }

    /// The radiated f0's peak swing at the gesture rate with no key, at `a = 0.04` and at `a = 0.08`, and the largest
    /// swing off the rate at `a = 0.08`, in cents. Measured on the radiated pressure of a 1.3 s render.
    fn radiated_tension_swing() -> (f64, f64, f64, f64) {
        let rate = 48_000.0_f64;
        let gesture_hz = 6.0_f64;
        let steady = render(&modulation_text(""));
        let prior = steady_f0_hz(&steady[(0.4 * rate) as usize..(0.9 * rate) as usize], rate);
        // Five whole periods of the gesture, after the onset transient.
        let (i0, i1) = (
            (0.30 * rate) as usize,
            (0.30 * rate + 5.0 / gesture_hz * rate) as usize + 4096,
        );
        let measure = |pcm: &[f64]| -> (f64, f64, f64) {
            // 83 hops of 10 ms: five periods of the gesture to within one hop, so the 3.6 and 8.4 Hz
            // projections are (nearly) whole-period nulls.
            let track: Vec<(f64, f64)> = f0_track(&pcm[i0..i1], rate, prior).into_iter().take(83).collect();
            assert_eq!(track.len(), 83, "the render is long enough for five periods");
            let mean = track.iter().map(|p| p.1).sum::<f64>() / track.len() as f64;
            let cents: Vec<(f64, f64)> = track
                .iter()
                .map(|&(t, f)| (t, 1200.0 * det::ln(f / mean) / core::f64::consts::LN_2))
                .collect();
            let at = amplitude_at(&cents, gesture_hz);
            let off = amplitude_at(&cents, gesture_hz - 2.4).max(amplitude_at(&cents, gesture_hz + 2.4));
            (at, off, mean)
        };
        let (none, _, mean_none) = measure(&steady);
        let (small, _, _) = measure(&render(&modulation_text("tension_modulation 6 0.04\n")));
        let (full, off, mean_full) = measure(&render(&modulation_text("tension_modulation 6 0.08\n")));
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-tension-modulation\",\"prior_hz\":{prior:.2},\"steady_hz\":{mean_none:.2},\"mean_hz_modulated\":{mean_full:.2},\"peak_cents_none\":{none:.2},\"peak_cents_a0p04\":{small:.1},\"peak_cents_a0p08\":{full:.1},\"max_cents_off_rate\":{off:.1}}}"
        );
        (none, small, full, off)
    }

    /// Gauntlet G1 with a planted negative. The radiated f0 swings at the gesture rate: nothing without the key, and
    /// DOUBLING the amplitude about doubles the extent.
    #[test]
    fn tension_modulation_swings_the_radiated_pitch_at_the_gesture_rate() {
        let (none, small, full, off) = radiated_tension_swing();
        assert!(none < 0.2 * full, "no key, yet {none:.1} cents at the gesture rate");
        assert!(
            (1.6..=2.4).contains(&(full / small)),
            "doubling the amplitude gave {small:.1} -> {full:.1} cents"
        );
        assert!(full > 4.0 * off, "{full:.1} cents at the rate against {off:.1} cents off it");
    }

    /// The same extent against `1731 e a` cents peak (`e` the square-root exponent, 0.3..0.8).
    #[test]
    #[ignore = "reads 41.0 cents for a = 0.08 against the floor of 41.5, a 0.5 cent miss; a pitch step on the a-C4 q grid reaches 3.93 cents, which is wider; bead frankensim-irp-v5-instrument-realism-pexfp.325.2"]
    fn tension_modulation_extent_stands_in_the_square_root_band() {
        let (_, _, full, _) = radiated_tension_swing();
        let per_a = 1200.0 / core::f64::consts::LN_2;
        assert!(
            (0.3 * 0.08 * per_a..=0.8 * 0.08 * per_a).contains(&full),
            "f0 swing {full:.1} cents peak for a = 0.08"
        );
    }

    // Moving tract (pexfp.327.1). Test charts, not data: an /a/-like card shape, narrow behind and wide in front, and
    // an /i/-like keyframe of ten sections, wide behind and narrow in front, on one glottal section.
    const SHAPE_A: [f64; 8] = [0.5e-4, 0.5e-4, 0.6e-4, 0.8e-4, 3.0e-4, 5.0e-4, 5.0e-4, 4.0e-4];
    const SHAPE_I: [f64; 10] = [0.5e-4, 3.0e-4, 4.0e-4, 4.0e-4, 3.0e-4, 1.5e-4, 0.4e-4, 0.3e-4, 0.5e-4, 1.0e-4];

    fn moving_card(sections: &[f64], section_m: f64, seconds: f64, extra: &str) -> String {
        let mut text = format!(
            "frankensim-voice-performance-v1\n\
             audio 48000 {} 1\n\
             full_scale_pa 2.0\n\
             ambient moist_air 310.15 101325 1.0\n\
             fold two_mass 1.25e-4 2.5e-5 128 12.8 40 0.1 3.57e-4 1.4e-2 2.5e-3\n\
             fold_damping_upper 0.6\n\
             termination flanged\n\
             tract_wall 15 1.48e+06 1.6e4\n\
             tract_source test-shape\n\
             tract_licence test fixture, not data\n\
             microphone 0 0 0.3 4 16 4000\n\
             subglottal 0 0\n\
             subglottal 0.02 900\n",
            (48_000.0 * seconds).round() as u64
        );
        for area in sections {
            text.push_str(&format!("section {area:e} {section_m:e}\n"));
        }
        text.push_str(extra);
        text
    }

    fn keyframe_i() -> String {
        let mut rows = String::from("keyframe i source test-i\nkeyframe i licence test fixture, not data\n");
        for area in SHAPE_I {
            rows.push_str(&format!("keyframe_section i {area:e} 1.7e-2\n"));
        }
        rows
    }

    fn render_performance(mut performance: VoicePerformance, samples: usize) -> Vec<f64> {
        let mut out = vec![0.0; samples];
        for block in out.chunks_mut(4096) {
            performance.block(block).expect("renders");
        }
        out
    }

    /// Largest second difference of `p` over `range`.
    fn largest_curvature(p: &[f64], range: core::ops::Range<usize>) -> f64 {
        range.map(|n| (p[n + 1] - 2.0 * p[n] + p[n - 1]).abs()).fold(0.0, f64::max)
    }

    /// F1 and F2 of a window: autocorrelation LPC of order 12 on the window decimated to 12 kHz (a 101-tap
    /// Hann-windowed sinc at 5 kHz first), pre-emphasized and Hann-windowed; the first two peaks of the envelope
    /// above 200 Hz on a 5 Hz grid.
    fn lpc_formants(x: &[f64]) -> (f64, f64) {
        let (rate, half) = (48_000.0, 50usize);
        let cutoff = 5_000.0 / rate;
        let taps: Vec<f64> = (0..=2 * half)
            .map(|k| {
                let m = k as f64 - half as f64;
                let sinc = if m == 0.0 { 2.0 * cutoff } else { (core::f64::consts::TAU * cutoff * m).sin() / (core::f64::consts::PI * m) };
                sinc * 0.5 * (1.0 - (core::f64::consts::TAU * k as f64 / (2 * half) as f64).cos())
            })
            .collect();
        let decimated: Vec<f64> = (half..x.len() - half)
            .step_by(4)
            .map(|n| taps.iter().enumerate().map(|(k, h)| h * x[n + half - k]).sum())
            .collect();
        let len = decimated.len();
        let y: Vec<f64> = (0..len)
            .map(|n| {
                let pre = decimated[n] - if n > 0 { 0.9 * decimated[n - 1] } else { 0.0 };
                pre * 0.5 * (1.0 - (core::f64::consts::TAU * n as f64 / (len - 1) as f64).cos())
            })
            .collect();
        let order = 12;
        let r: Vec<f64> = (0..=order).map(|lag| (lag..len).map(|n| y[n] * y[n - lag]).sum()).collect();
        // Levinson-Durbin: A(z) = 1 + a_1 z^-1 + ... + a_p z^-p.
        let mut a = vec![0.0; order + 1];
        a[0] = 1.0;
        let mut error = r[0];
        for i in 1..=order {
            let acc: f64 = (1..i).map(|j| a[j] * r[i - j]).sum();
            let k = -(r[i] + acc) / error;
            let previous = a.clone();
            for j in 1..i {
                a[j] = previous[j] + k * previous[i - j];
            }
            a[i] = k;
            error *= 1.0 - k * k;
        }
        let envelope = |f: f64| {
            let w = core::f64::consts::TAU * f / 12_000.0;
            let (re, im) = a.iter().enumerate().fold((0.0, 0.0), |(re, im), (k, c)| {
                (re + c * (w * k as f64).cos(), im - c * (w * k as f64).sin())
            });
            1.0 / (re * re + im * im)
        };
        let mut peaks = Vec::new();
        let mut f = 205.0;
        while f < 5_000.0 && peaks.len() < 2 {
            if envelope(f) > envelope(f - 5.0) && envelope(f) > envelope(f + 5.0) {
                peaks.push(f);
            }
            f += 5.0;
        }
        assert_eq!(peaks.len(), 2, "two envelope peaks below 5 kHz");
        (peaks[0], peaks[1])
    }

    #[test]
    fn tract_move_rows_parse_and_are_refused_by_name() {
        let ok = moving_card(&SHAPE_A, 0.02, 0.2, &(keyframe_i() + "tract_move 0.05 0.1 i 8\ntract_move 0.12 0.15 card 4\n"));
        let description = VoiceDescription::from_bytes(ok.as_bytes()).expect("parses");
        assert_eq!(description.keyframes.len(), 1);
        assert_eq!(description.keyframes[0].sections.len(), 10);
        assert_eq!(
            description.tract_moves,
            vec![
                TractMove { start_s: 0.05, end_s: 0.1, target: Some(0), steps: 8 },
                TractMove { start_s: 0.12, end_s: 0.15, target: None, steps: 4 },
            ]
        );
        let refused = |extra: String, why: &str| {
            let error = VoiceDescription::from_bytes(moving_card(&SHAPE_A, 0.02, 0.2, &extra).as_bytes())
                .expect_err("refused");
            assert!(error.contains(why), "{error}");
        };
        refused(keyframe_i(), "named by no tract_move");
        refused(keyframe_i() + "tract_move 0.05 0.1 u 8\n", "unknown keyframe");
        refused(keyframe_i().replacen("keyframe_section i 5e-5", "keyframe_section i 6e-5", 1) + "tract_move 0.05 0.1 i 8\n", "glottal section is held constant");
        refused(keyframe_i().replace("keyframe i source test-i\n", "") + "tract_move 0.05 0.1 i 8\n", "needs a `source` row");
        refused(keyframe_i().replace("keyframe i licence test fixture, not data\n", "") + "tract_move 0.05 0.1 i 8\n", "needs a `licence` row");
        refused("keyframe card source x\n".into(), "cannot be a keyframe");
        refused(keyframe_i() + "tract_move 0.05 0.1 card 8\ntract_move 0.12 0.15 i 4\n", "must change the shape");
        refused(keyframe_i() + "tract_move 0.05 0.1 i 8\ntract_move 0.09 0.15 card 4\n", "not overlap");
        refused(keyframe_i() + "tract_move 0.05 0.3 i 8\n", "end inside the render");
        refused(keyframe_i() + "tract_move 0.05 0.1 i 0\n", "steps must be in");
        refused(keyframe_i() + "keyframe i source again\ntract_move 0.05 0.1 i 8\n", "may only appear once");
        // At least one mechanics step per shape: 10 steps of 1/48000 s for 16 shapes.
        let crowded = moving_card(&SHAPE_A, 0.02, 0.2, &(keyframe_i() + "tract_move 0.05 0.0502 i 16\n"));
        let error = VoiceDescription::from_bytes(crowded.as_bytes()).expect("parses").build().err().expect("refused");
        assert!(error.contains("at least one mechanics step per shape"), "{error}");
    }

    /// An island given no further shapes plays bit for bit as the plain island: the moving path with one shape is the
    /// fixed one.
    #[test]
    fn an_island_with_one_tract_shape_is_bit_identical_to_the_plain_island() {
        let description = VoiceDescription::from_bytes(moving_card(&SHAPE_A, 0.02, 0.1, "").as_bytes()).expect("parses");
        let chart = TractChart::try_new(description.sections.clone(), &description.tract_source, &description.tract_licence)
            .expect("chart");
        let duct = chart.to_duct();
        let island = || {
            GlottalIsland::new_with_wall(
                description.fold, description.two_mass, &duct, description.tract_gas(), description.termination, 48_000,
                description.wall.as_ref(),
            )
            .expect("island")
        };
        let mut plain = island();
        let mut one = island().with_tract_keyframes(Vec::new()).expect("one shape");
        assert!(one.switch_tract(1, 0).is_err(), "no shape 1");
        for n in 0..20_000 {
            let p = 900.0 * (n as f64 / 960.0).min(1.0);
            let (a, b) = (plain.step(p).expect("plain"), one.step(p).expect("one shape"));
            assert_eq!(a.flow_m3_s.to_bits(), b.flow_m3_s.to_bits(), "step {n}");
            assert_eq!(plain.incident_wave_pa().to_bits(), one.incident_wave_pa().to_bits(), "step {n}");
        }
        assert!(plain.with_tract_keyframes(Vec::new()).is_err(), "keyframes come before the first step");
    }

    /// Test group 4 of the plan line: a tract moving from the /a/-like card shape to the /i/-like keyframe over 0.25 to
    /// 0.35 s is bit for bit the fixed card before the move, has no click, and after the move has the keyframe's
    /// formants. The same move as one hard switch (fade 0) must click, or the click test could not fail.
    #[test]
    fn a_moving_tract_has_no_click_and_reaches_the_keyframe_formants() {
        let (seconds, samples) = (0.6, 28_800usize);
        let moving_text = moving_card(&SHAPE_A, 0.02, seconds, &(keyframe_i() + "tract_move 0.25 0.35 i 20\n"));
        let moving = VoiceDescription::from_bytes(moving_text.as_bytes()).expect("parses");
        let fixed = VoiceDescription::from_bytes(moving_card(&SHAPE_A, 0.02, seconds, "").as_bytes()).expect("parses");
        let static_i = VoiceDescription::from_bytes(moving_card(&SHAPE_I, 0.017, seconds, "").as_bytes()).expect("parses");
        let started = std::time::Instant::now();
        let fixed_build = fixed.build().expect("fixed builds");
        let fixed_s = started.elapsed().as_secs_f64();
        let started = std::time::Instant::now();
        let moving_build = moving.build().expect("moving builds");
        let moving_s = started.elapsed().as_secs_f64();
        assert_eq!(moving_build.info().tract_shapes, 21);
        assert_eq!(moving_build.info().keyframe_formants_hz.len(), 1);
        let fixed_p = render_performance(fixed_build, samples);
        let moving_p = render_performance(moving_build, samples);
        let static_p = render_performance(static_i.build().expect("static builds"), samples);
        let mut hard = VoiceDescription::from_bytes(
            moving_card(&SHAPE_A, 0.02, seconds, &(keyframe_i() + "tract_move 0.25 0.35 i 1\n")).as_bytes(),
        )
        .expect("parses")
        .build()
        .expect("hard builds");
        assert_eq!(hard.tract_switches.len(), 1);
        hard.tract_switches[0].fade = 0;
        let hard_p = render_performance(hard, samples);
        // Nothing switches before 0.25 s: those samples are the fixed card's, bit for bit.
        let first = 12_000;
        assert!(moving_p[..first].iter().zip(&fixed_p[..first]).all(|(a, b)| a.to_bits() == b.to_bits()));
        assert!(moving_p[first..].iter().zip(&fixed_p[first..]).any(|(a, b)| a != b), "the move changes the sound");
        // No click: the move's largest curvature stays within 1.5 times the larger of the two held shapes'.
        let held = largest_curvature(&moving_p, 4_800..12_000).max(largest_curvature(&moving_p, 21_600..28_799));
        let moved = largest_curvature(&moving_p, 12_000..16_800);
        let held_hard = largest_curvature(&hard_p, 4_800..12_000).max(largest_curvature(&hard_p, 21_600..28_799));
        let clicked = largest_curvature(&hard_p, 12_000..16_800);
        // Formants after the move against the static keyframe, and the two shapes' own, which must differ.
        let after = lpc_formants(&moving_p[21_600..28_800]);
        let keyframe = lpc_formants(&static_p[21_600..28_800]);
        let card = lpc_formants(&fixed_p[21_600..28_800]);
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-moving-tract\",\"moved_over_held\":{:.3},\"hard_over_held\":{:.3},\
             \"after_hz\":[{:.0},{:.0}],\"keyframe_hz\":[{:.0},{:.0}],\"card_hz\":[{:.0},{:.0}],\
             \"further_shapes_wall_s\":{:.3}}}",
            moved / held, clicked / held_hard, after.0, after.1, keyframe.0, keyframe.1, card.0, card.1,
            moving_s - fixed_s
        );
        assert!(moved <= 1.5 * held, "a click: {moved:e} against held {held:e}");
        assert!(clicked > 1.5 * held_hard, "the hard switch did not click ({clicked:e} against {held_hard:e}): the test cannot fail");
        for (x, y) in [(after.0, keyframe.0), (after.1, keyframe.1)] {
            assert!((x - y).abs() <= 0.03 * y, "after the move {x:.0} Hz against the keyframe's {y:.0} Hz");
        }
        assert!(
            (card.0 - keyframe.0).abs() > 0.1 * keyframe.0 || (card.1 - keyframe.1).abs() > 0.1 * keyframe.1,
            "the two test shapes must differ in F1 or F2: {card:?} against {keyframe:?}"
        );
    }

    /// The further shapes built on one thread and on three (blocks of 3, 3 and 3 over two moves) give the same switches,
    /// the same targets' formants and the same samples, bit for bit.
    #[test]
    fn a_moving_tract_builds_the_same_on_any_number_of_threads() {
        let (seconds, samples) = (0.3, 14_400usize);
        let text = moving_card(&SHAPE_A, 0.02, seconds, &(keyframe_i() + "tract_move 0.05 0.1 i 5\ntract_move 0.15 0.2 card 4\n"));
        let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
        let (one, three) = (
            description.build_with_shape_workers(1).expect("one thread builds"),
            description.build_with_shape_workers(3).expect("three threads build"),
        );
        assert_eq!(one.info().tract_shapes, 10);
        let switches = |p: &VoicePerformance| p.tract_switches.iter().map(|s| (s.step, s.shape, s.fade)).collect::<Vec<_>>();
        assert_eq!(switches(&one), switches(&three));
        let formants = |p: &VoicePerformance| {
            p.info()
                .keyframe_formants_hz
                .iter()
                .map(|(name, f)| (name.clone(), f.as_ref().map(|v| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>()).ok()))
                .collect::<Vec<_>>()
        };
        assert_eq!(formants(&one).len(), 2);
        assert_eq!(formants(&one), formants(&three));
        let (a, b) = (render_performance(one, samples), render_performance(three, samples));
        assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()), "the thread count changed a sample");
    }

    /// Build time per shape at the live cards' clock and rows (head sphere and load, torso, tract gas, 192 kHz), which a
    /// render's stop line reads before any moving card renders. The four further shapes build at once on the host's
    /// threads, so their wall is one round of up to `MAX_SHAPE_WORKERS` shapes, and a move of `n` shapes builds in about
    /// `fixed_build_s + ceil(n / shape_workers) * further_shapes_wall_s`. Test shapes of 90 sections, not data; nothing
    /// is rendered. The keyframe's
    /// mouth stays at 3.0e-4 m^2 (9.8 mm), inside the live cards' 8.6-14.0 mm: a 1.0e-4 m^2 mouth (5.6 mm) on the
    /// 8.7 cm head needs more than the head load's 65536 orders, and the build refuses it by name.
    #[test]
    #[ignore = "timing probe for the build host: run with --ignored --nocapture"]
    fn moving_tract_build_time_per_shape_at_the_live_clock() {
        let a: Vec<f64> = (0..90).map(|k| if k == 0 { 1.4e-4 } else { 1.0e-4 + 4.0e-4 * (k as f64 / 89.0) }).collect();
        let mut extra = String::from("keyframe i source test-i\nkeyframe i licence test fixture, not data\n");
        for k in 0..98 {
            let area = if k == 0 { 1.4e-4 } else { 5.0e-4 - 2.0e-4 * (k as f64 / 97.0) };
            extra.push_str(&format!("keyframe_section i {area:e} 2.0e-3\n"));
        }
        let rows = "head_sphere 0.087 test\nhead_load cap\ntorso_sphere 0.169 0.053 test\n\
                    tract_gas moist_air 310.15 101325 1.0\n";
        let card = |motion: &str| {
            moving_card(&a, 2.0e-3, 0.02, &format!("{rows}{motion}"))
                .replace("audio 48000 960 1", "audio 48000 960 4")
                .replace("microphone 0 0 0.3 4 16 4000", "microphone 0.212132 0 0.212132 16 64 8000")
        };
        let fixed = VoiceDescription::from_bytes(card("").as_bytes()).expect("parses");
        let moving = VoiceDescription::from_bytes(card(&(extra + "tract_move 0.005 0.015 i 4\n")).as_bytes()).expect("parses");
        let started = std::time::Instant::now();
        fixed.build().expect("fixed builds");
        let fixed_s = started.elapsed().as_secs_f64();
        let started = std::time::Instant::now();
        let built = moving.build().expect("moving builds");
        let moving_s = started.elapsed().as_secs_f64();
        assert_eq!(built.info().tract_shapes, 5);
        let workers = std::thread::available_parallelism().map_or(1, |p| p.get()).clamp(1, MAX_SHAPE_WORKERS);
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-moving-tract-build-time\",\"fixed_build_s\":{fixed_s:.2},\
             \"moving_build_s\":{moving_s:.2},\"shape_workers\":{workers},\"further_shapes_wall_s\":{:.2}}}",
            moving_s - fixed_s
        );
    }

    /// The `subglottal_airway` row parses, is refused by name, and reaches the folds: the airway is realized on the
    /// mechanics clock in the tract's gas with the card's wall, its receipt is reported, and the radiated pressure
    /// moves; without the row nothing is reported and nothing moves.
    #[test]
    fn subglottal_airway_parses_is_refused_by_name_and_loads_the_folds() {
        let row = "subglottal_airway weibel-1963 saab-2021-table-1-B 17\n";
        let description = VoiceDescription::from_bytes(fixture(row).as_bytes()).unwrap();
        assert_eq!(description.subglottal_airway, Some(17));
        assert_eq!(VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap().subglottal_airway, None);
        for rows in [
            "subglottal_airway weibel-1963 saab-2021-table-1-B\n",
            "subglottal_airway horsfield-1971 saab-2021-table-1-B 17\n",
            "subglottal_airway weibel-1963 saab-2021-table-1-A 17\n",
            "subglottal_airway weibel-1963 saab-2021-table-1-B 0\n",
            "subglottal_airway weibel-1963 saab-2021-table-1-B 25\n",
            "subglottal_airway weibel-1963 saab-2021-table-1-B 17.5\n",
            "subglottal_airway weibel-1963 saab-2021-table-1-B 17\nsubglottal_airway weibel-1963 saab-2021-table-1-B 16\n",
        ] {
            let error = VoiceDescription::from_bytes(fixture(rows).as_bytes()).unwrap_err();
            assert!(error.contains("subglottal_airway"), "{rows:?}: {error}");
        }
        let wall = "tract_wall 15 1.48e+06 1.6e4\n";
        let plain = VoiceDescription::from_bytes(fixture(wall).as_bytes()).unwrap();
        let walled = VoiceDescription::from_bytes(fixture(&format!("{wall}{row}")).as_bytes()).unwrap();
        let performance = walled.build().expect("the airway realizes");
        let receipt = performance.info().subglottal_airway.clone().expect("receipt");
        assert_eq!(receipt.sections, 17);
        assert!(plain.build().unwrap().info().subglottal_airway.is_none());
        let (without, with) = (pressure(&plain), pressure(&walled));
        let moved = without.iter().zip(&with).any(|(a, b)| a.to_bits() != b.to_bits());
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-subglottal-airway-row\",\"zs\":{:e},\"r_dc\":{:e},\"r0\":{:e},\"reflectance_taps\":{},\"source_taps\":{}}}",
            receipt.port_impedance,
            receipt.dc_resistance,
            receipt.same_sample_reflection,
            receipt.reflectance.n_fft,
            receipt.source_taps
        );
        assert!(moved, "the airway must reach the radiated pressure");
    }
}
