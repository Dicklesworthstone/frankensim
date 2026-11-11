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
//! tract (no side branches or nasal coupling) unless the opt-in nasal rows open a
//! nasal tract beside it (`nasal`, a graph of plane-wave edges), a lumped two-mass, body-cover or
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
//! termination <flanged|unflanged|ideal_open>    or    termination unflanged_thin_wall <quadrature_nodes>
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
//! tract_gas <any `ambient` record>                        (optional; without head_load needs termination flanged and no nasal rows; the gas inside the tract)
//! keyframe <name> source <identifier without spaces>      (optional; a further measured shape of a moving tract)
//! keyframe <name> licence <free text to end of line>
//! keyframe_section <name> <area_m2> <length_m>            (repeat, glottis first)
//! tract_move <start_s> <end_s> <name|card> <steps>        (repeat, in time order; the tract moves to that shape)
//! glottal_orifice <area_m2>                               (optional; a fixed hole in place of the folds: a flow rig)
//! subglottal_tube <area_m2> <length_m>                    (repeat, optional; a tube below the glottis, glottis end first)
//! fricative_noise <seed> <k> <j> [stevens [<m>]]          (optional; turbulence after section k, or the obstacle m)
//! fricative_throat <area_m2> <source without spaces>      (optional, needs fricative_noise; the constriction's measured throat)
//! receiver free_field_monopole                            (optional; the mouth flow as a monopole in free space)
//! rest_gap <time_s> <gap_m>                               (repeat, optional, two_mass only; the fold's rest gap over time)
//! fold_tension <time_s> <stiffness_scale> <vibrato_scale> <flutter_scale>
//!                                                         (repeat, optional; the fold's tension over time)
//! nasal_source <identifier without spaces>                (optional; with any nasal row: the nasal tables)
//! nasal_licence <free text to end of line>
//! velum <after_section> <port_area_m2> <attach_source>    (the port opens after 0-based oral section after_section)
//! nasal_section <area_m2> <length_m>                      (repeat, velum first; the common nasal passage)
//! nasal_pair <left_area_m2> <right_area_m2> <length_m>    (repeat, choana first; the two passages to the nostrils)
//! sinus <left|right> <from_nostril_m> <ostium_area_m2> <ostium_length_m> <volume_m3>   (optional, one per side)
//! nostrils <x_m> <y_m> <source without spaces>            (the nostrils beside the mouth on the baffle)
//! ```
//!
//! The nasal rows (N3, MB 15368) open the nasal tract, and with `sinus` rows the maxillary sinuses, beside the card's
//! oral tract as one network; the lips and both nostrils radiate into one receiver. The graph, the port, the cavity
//! loads, the one-plane receiver and the rows refused with them are in the `nasal` child module's doc. Without the rows
//! nothing here changes.
//!
//! The flow rig and the supraglottal source (MB 15368; rulings 15773, 15828, 15902, 16025; interface NavyTower 15491,
//! 15765). `glottal_orifice` replaces the folds by a fixed hole of that area ([`GlottalIsland::with_orifice`]): the jet
//! with the two-mass exit's recovery into the first section (Ishizaka and Flanagan eq 6; a hole wider than the first
//! section is refused), in series with the tract, as a stationary flow rig blows through a printed model. The fold rows stay
//! (the grammar needs them) and reach no number of the render; the rows that move or read a fold are refused with it
//! (`fold_shape`, `tension_modulation`, `tension_flutter`, `pulse_shimmer`, and the glottal sources `aspiration_noise`,
//! `jet_noise` and `constriction_noise`, whose laws read a fold's gap or take the whole drive as their drop). A rig's
//! stats carry its mean glottal flow over the render's last 2 s ([`VoiceStats::rig_mean_glottal_flow_m3_s`]).
//! `subglottal_tube` rows are a tube below the glottis, one cylinder of each area and length, glottis end first, with
//! Bessel losses and the card's `tract_wall` (rigid without one), through the same seam as `subglottal_airway`, which is
//! refused with them: the `subglottal` schedule becomes the pressure at the tube's far end (a rig's pressure tap).
//! `fricative_noise` is Oldham and Ukpoho's printed law in its supraglottal form
//! ([`crate::turbulence_noise::ConstrictionNoise::new_fricative`]) at the constriction `k` (0-based section) opening into
//! the section `j`, which the card names (MB 15828: the largest section between `k` and the lips): a series pressure at
//! the junction after section `k`, reaching the mouth flow and the island's returned wave through the junction lines
//! ([`crate::driving_point::outlet::junction_source_lines`]), with the constriction's Borda-Carnot drop and Shadle's
//! viscous term `K3 U` (MB 18431; [`crate::glottis::shadle_constriction_k3`]) in the island's junction solve
//! ([`GlottalIsland::with_constriction`]), K3 less the Poiseuille resistance the line already carries over the sections
//! the constriction spans, every shape its own (MB 18484; [`GlottalIsland::with_constriction_spans`]). On the flow rig
//! its drop and jet speed are those of the held flow through the orifice and the constriction at the schedule's
//! pressure ([`GlottalIsland::orifice_held_flow_m3_s`]);
//! on the folds, which have no held flow, they are those of the island's own flow of the step before (the source reads
//! the flow it does not yet change, as the glottal sources do), and the source is silent while that flow is not positive
//! (no jet leaves `k` toward `j`). That flow is the glottal flow, the lumped law of the island's own drop at `k`, not the
//! volume velocity at `k`: in a voiced frame the noise follows the glottal pulse without the tract's filtering between
//! the glottis and `k` (NavyTower 16796). Its lines end in the island's own mouth load: the declared termination, or with
//! `head_load` the same table, read in the exterior with `tract_gas` (the bore, the law's density and the lines stay in
//! the tract's gas; [`crate::driving_point::outlet::junction_source_lines_exterior`]). On a moving tract (`tract_move`;
//! MB 16735 piece 2) every shape takes its own `k` by [`TractChart::log_area_morph`]'s own correspondence (MB 16891,
//! ruling (c) on NavyTower 16887): at the card's section count `k` itself, at any other count the shape's section at the
//! card's section `k`'s centre as a fraction of the tract's length from the glottis; its own junction lines after that
//! section, switched with the island's and crossfaded the same
//! way, and its own constriction, `A_k` that section and `A_j` the shape's largest section after it: the law's `sigma`,
//! `d` and cutoff from the switch's step on, and the island's Borda-Carnot factor and K3 with them
//! ([`GlottalIsland::with_constriction_keyframes`]; neither is faded). A shape with no expansion after its `k`, the
//! card's own included (a word's card is its vowel), has `sigma >= 1`, where the printed law drives no band, and the
//! island takes no constriction on it: the source forms and turns off through the law as the tract closes and opens,
//! with no fade added. With `tract_move` the card's `j` must be the card shape's largest section after `k`; without it
//! the card's one shape must be an expansion. A fourth field `stevens` (MB 17504, 17645) takes Stevens' printed source at
//! an obstacle from Shadle's measurements in place of Oldham and Ukpoho's
//! ([`crate::turbulence_noise::ConstrictionNoise::pressure_pa_stevens`]): the same junction, lines, bank, streams, held flow
//! and shapes, its level read from the flow and `A_k` alone, silent on a shape with `sigma >= 1` as the printed law is
//! there; without the field a card is bit-identical to one written before it existed. Stevens' `p_s` is defined at the
//! obstacle (Shadle's source placement), so after `stevens` a fifth field `m` (MB 17712) puts the junction lines after
//! section `m` instead of `k`, `k <= m` and a section after `m`: the obstacle the card names (the incisors, measured on
//! the tract's own surface), with the law's level still read from the flow and `A_k` and the Borda-Carnot drop and K3
//! still at `k`; without it `m = k`, bit-identical. On a moving tract (MB 25469) each shape's `m` is the card's `m`
//! by the same correspondence as its `k` (`shape_constriction`: `m` itself at the card's section count, else the
//! shape's section at the card's `m`'s centre fraction), so it is never before the shape's `k`; a shape whose `m`
//! falls on its last section is refused by name, as one whose `k` does. With `m = k` every shape's lines and `A_m` are
//! the ones it took before `m` was allowed to move. The printed source is the source in Shadle's 5.067 cm^2 tube; the
//! law carries Shadle's eq. 3.3 duct-area factor (MB 18252; [`crate::turbulence_noise::shadle_duct_area_factor`]): at
//! the junction after section `m` the source is the printed one times `A_tube / A_m`, `A_m` the area of section `m`
//! (on a moving tract, each shape's own section `m`).
//! `fricative_throat <area_m2> <source>` (MB 25440 ruling A, 25482) gives the constriction laws
//! the throat measured on the tract's own surface by a flat-section route, the route and its read named in the source
//! field: the tract's sections are an acoustic-equivalent area function (a level cut keeps each section's volume and
//! energy), which is not the geometric throat the jet passes, and one number cannot serve both. The row moves `A_k` in
//! the constriction laws only: the island's Borda-Carnot factor `(1/A_k - 1/A_j)^2` and Shadle's K3 at `A_k`
//! ([`GlottalIsland::with_constriction`]), so the rig's held flow with them, and the source's `A_k` and `sigma = A_k /
//! A_j` (Oldham and Ukpoho's jet speed, Stevens' level). The line's sections, section `k` included, the Poiseuille
//! resistance taken out of K3 over the line's span, `A_j` and `A_m` stay the tract's. It is refused without
//! `fricative_noise`, with `tract_move` (a moving shape has no measured throat), and for an area that is not finite and
//! positive, not narrower than section `j`, or more than 1.5 times section `k` (a stated sanity bound on a measured
//! value, not a fit); without the row a card is bit-identical to one written before it existed.
//! `receiver free_field_monopole` observes
//! the mouth flow as a point monopole in free space, `rho dQ/dt / (4 pi R)` at the microphone's distance, through the
//! baffled receiver's own differentiator and travel-time kernel (half a point source on the rigid baffle, whose image
//! doubles it); the microphone row's rings and points are not read by it, the mouth's size and the directivity of its
//! opening are not modelled, and `head_sphere` and `observer_law` are refused with it. `termination unflanged_thin_wall
//! <nodes>` is Levine and Schwinger's exact thin-wall unflanged load on that many quadrature nodes, refused by fs-duct
//! above its first cut-on (`ka` 3.83), so the mechanics clock must keep the lip section below it. Without these rows
//! every line takes exactly its old path.
//!
//! `rest_gap` rows move a `two_mass` fold's rest gap over time, as an abducting or adducting gesture of the
//! arytenoids does (the first word, MB 16735): linear between rows in strictly increasing time, the first row's gap
//! held before it and the last row's after it, each gap the whole opening as the `fold` row's own `rest_gap_m` (with a
//! `fold_shape`, the mean over the fold's length). The gap is set once per output sample before the mechanics steps,
//! as the tension gestures are ([`GlottalIsland::set_rest_gap`]): the displacement from rest is kept, and the springs'
//! energy change is reported as work on the fold, so the ledger stays closed. The rows say nothing about how the gap
//! moves between them; a card states where its gaps and times come from. `body_cover` and `one_dof` folds and
//! `glottal_orifice` refuse them by name. Without the rows the fold keeps the `fold` row's gap and a performance is
//! bit-identical to one written before the rows existed.
//!
//! `fold_tension` rows move the fold's tension over time, as a singer does from one note or vowel to the next within
//! one render (the lyric phrase, MB 23138): linear between rows in strictly increasing time, the first row held before
//! it and the last row after it, as `rest_gap` is. At time `t` the three factors `s`, `g_v` and `g_f` are read, and
//! the fold's stiffnesses take `s (1 + g_v a sin(2 pi r t) + g_f y)` in place of `1 + a sin(2 pi r t) + y`, where
//! `a` and `r` are the `tension_modulation` row's and `y` is the flutter rows' sum ([`GlottalIsland::set_tension`],
//! once per output sample, so the springs' energy change is reported as work on the fold as the gestures' already is).
//! A card states where its factors come from; nothing here picks them. A vibrato scale other than 1 needs
//! `tension_modulation`, a flutter scale other than 1 needs `tension_flutter`, the scaled vibrato amplitude and the
//! scaled flutter sds keep their rows' own bounds (0.5 and 0.1), the stiffness factor lies in the fold's admitted
//! tension range (1/4 to 4), and `glottal_orifice` refuses the rows by name. The product with the gestures must stay
//! inside that range too; a render that leaves it is refused there. Without the rows nothing is read or
//! called and a performance is bit-identical to one written before the rows existed; rows at `1 1 1` are bit-identical
//! too, since a product with 1 is exact.
//!
//! A moving tract (pexfp.327.1): each `tract_move` takes the tract from the shape it holds (the card's own `section`
//! rows until the first move) to keyframe `name`, or back to the card's shape with `card`, through `steps` shapes at
//! fractions `k / steps` of [`TractChart::log_area_morph`] (log area interpolated, length linear). Shape `k`'s
//! crossfade starts at `start + (k - 1) (end - start) / steps` and lasts `(end - start) / steps`, so the tract is
//! exactly the target from `end`. Every shape's island line and mouth-flow line are realized at build time by the same
//! ladder, load, wall and gases as the card's own shape, and the render switches them with the outgoing history carried
//! ([`crate::keyframe_line`]). The glottal (first) section moves with the rest (MB 25468): a keyframe's first area is
//! its own, every shape's is the morph's (with unequal section counts the length-weighted mean over the first of the
//! common spans, which can take in part of the second section), and each shape's island line, mouth-flow line and
//! fricative junction lines are referenced to that shape's own `rho c / A_1`. Three limits, not hidden: the island's `zc` and the two-mass exit's
//! `A_1` (eq 6) step to the new shape's at its switch while the reflected wave crossfades, so they are not faded; the
//! microphone, the receiver and the head sphere keep the card's mouth radius; and the shapes are frozen linear tracts
//! crossfaded, coefficient interpolation with no volume source from the moving walls a real tract has. A keyframe that
//! no move names is refused. Without the rows the build and the render are bit for bit what they were.
//!
//! `tract_gas` fills the tract with its own gas (breath) while `ambient` stays the exterior (the room). The island's
//! line, the mouth-flow line, the passive formants and the glottal flow's density and viscosity take the tract gas; the
//! receiver and the head and torso spheres take `ambient`. The lips' load table is an absolute impedance (Pa s/m^3)
//! bound to the exterior and checked against it, so nothing is converted at the lips: pressure and volume flow are
//! continuous there. Its flanged part is fs-duct's two-gas law, whose radiation is in the exterior and whose eliminated
//! higher duct modes are in the tract's gas. Without `head_load` the table is that flanged load alone, on the same bins
//! ([`crate::pcm_wav::baffled::head::flanged_mouth_load`]; MB 23833, 24424); with it the head's cap, piston and torso
//! terms, all in the exterior, are added to that flanged load, and each nostril's table takes the tract's gas the same
//! way ([`crate::pcm_wav::baffled::head::head_mouth_load`]; MB 25423). Without `head_load` the row needs
//! `termination flanged`; a nasal card with a tract gas and no head is untried and refused. The two gases share one
//! static pressure (mean flow is not modeled). The medium is stationary: the breath does not
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
use crate::driving_point::outlet::{
    junction_source_lines, junction_source_lines_exterior, mouth_flow_line, mouth_flow_line_exterior,
};
use crate::glottis::{
    FoldBody, FoldCard, FoldCardSource, FoldMode, FoldRestGap, FoldShape, GLOTTAL_LINE_MAX_TAPS, GlottalIsland,
    constriction_span, glottal_tract_line_taps,
};
use crate::keyframe_line::ShapeLine;
use crate::pcm_wav::baffled::head::{
    HeadSphere, TorsoSphere, flanged_mouth_load, head_mouth_load, head_opening_load, head_torso_mouth_load,
};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium, SurfaceSample};
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
use crate::upstream_airway::{AirwayChain, AirwayReceipt, AirwaySection, UpstreamAirway, WEIBEL_SAAB_2021_TABLE_1_B_MM};
use fs_duct::levine_schwinger::LsQuadrature;
use fs_duct::{LossModel, MouthLoad, Termination};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources, read_gas_record};
use fs_math::det;
use fs_phs::WallPin;
use fs_rand::{Stream, StreamKey};
use fs_vfit::discretize::DelayedFilter;

mod nasal;
pub use nasal::{NasalCard, NasalPair, NasalSinus};

/// First line of every voice performance file.
pub const VOICE_PERFORMANCE_SCHEMA: &str = "frankensim-voice-performance-v1";
const HEAD_LOAD_NEEDS_SPHERE: &str = "head_load needs the head_sphere row: the load is that sphere's";
const HEAD_LOAD_NEEDS_FLANGED: &str = "head_load needs `termination flanged`: it adds the sphere's change to the flanged load";
const TORSO_NEEDS_SPHERE: &str = "torso_sphere needs the head_sphere row: the torso sits under that head";
const TRACT_GAS_NEEDS_FLANGED: &str =
    "tract_gas without head_load needs `termination flanged`: the lips' flanged load is computed in the exterior";
const TRACT_GAS_NASAL_NEEDS_HEAD_LOAD: &str =
    "tract_gas with the nasal rows needs the head_load row: a nasal card with a tract gas and no head is untried";
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
    /// Opt-in `glottal_orifice` [m^2]: a fixed hole in place of the folds (module doc); `None` is the fold.
    pub glottal_orifice: Option<f64>,
    /// Opt-in `subglottal_tube` rows, glottis end first (module doc); empty for none.
    pub subglottal_tube: Vec<TractSection>,
    /// Opt-in `fricative_noise` (module doc); `None` adds nothing.
    pub fricative_noise: Option<FricativeNoiseCard>,
    /// Opt-in `fricative_throat` (module doc; MB 25440 ruling A): the constriction laws' measured `A_k`; `None` keeps
    /// section `k`'s own area.
    pub fricative_throat: Option<FricativeThroat>,
    /// Opt-in `receiver free_field_monopole` (module doc); `false` keeps the infinite baffle.
    pub free_field_receiver: bool,
    /// Opt-in `rest_gap` rows `(time_s, gap_m)` in strictly increasing time (module doc); empty keeps the `fold` row's gap.
    pub rest_gap: Vec<(f64, f64)>,
    /// Opt-in `fold_tension` rows in strictly increasing time (module doc); empty keeps the card's own tension.
    pub fold_tension: Vec<FoldTensionPoint>,
    /// Opt-in nasal rows (module doc; the `nasal` child module); `None` keeps the series tract and the one opening.
    pub nasal: Option<NasalCard>,
}

/// One `fold_tension <time_s> <stiffness_scale> <vibrato_scale> <flutter_scale>` row (module doc): at `time_s` the
/// fold's stiffnesses are scaled by `stiffness`, the `tension_modulation` amplitude by `vibrato` and the flutter rows'
/// output by `flutter`; linear between rows. The card names all four numbers; this code chooses none.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FoldTensionPoint {
    /// Time of the row [s].
    pub time_s: f64,
    /// Common factor on the fold's stiffnesses.
    pub stiffness: f64,
    /// Factor on the `tension_modulation` amplitude.
    pub vibrato: f64,
    /// Factor on the flutter rows' summed output.
    pub flutter: f64,
}

/// The `fricative_noise <seed> <k> <j> [stevens [<m>]]` row: the supraglottal constriction `k` (0-based section index)
/// opening into the section `j > k`, the printed law, and with `stevens` the obstacle section `m >= k`. The card names all
/// four; this code chooses none.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FricativeNoiseCard {
    /// Declared seed of the noise streams.
    pub seed: u64,
    /// The constriction's section.
    pub constriction: usize,
    /// The section it opens into.
    pub expansion: usize,
    /// The printed source law; a row without the fourth field is Oldham and Ukpoho's.
    pub law: FricativeLaw,
    /// The section after which the source's junction lines sit: the fifth field, else `constriction`.
    pub obstacle: usize,
}

/// The `fricative_throat <area_m2> <source>` row (module doc; MB 25440 ruling A, 25482): the geometric throat of the
/// `fricative_noise` constriction, measured on the tract's own surface, and the route and read that measured it. The card
/// names both; this code chooses neither.
#[derive(Debug, Clone, PartialEq)]
pub struct FricativeThroat {
    /// The throat's area [m^2]: the constriction laws' `A_k`.
    pub area_m2: f64,
    /// The route and its read (an identifier without spaces).
    pub source: String,
}

/// The printed law of the `fricative_noise` source (module doc).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FricativeLaw {
    /// Oldham and Ukpoho (1990), the duct-constriction law in its supraglottal form: no fourth field.
    OldhamUkpoho,
    /// Stevens, *Acoustic Phonetics*, Fig. 2.33, on Shadle's (1985a) data: the fourth field `stevens`.
    StevensShadle,
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
    /// Area function, glottis first; the first area is the keyframe's own (the glottal section moves with the tract).
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

/// The area [m^2] of the largest section after section `k` (`fricative_noise`'s `A_j` on a moving tract, MB 15828: the
/// largest section between `k` and the lips); 0 when none follows.
fn widest_after(sections: &[TractSection], k: usize) -> f64 {
    sections.iter().skip(k + 1).map(|s| s.area_m2).fold(0.0, f64::max)
}

/// The line sections `fricative_noise`'s constriction at section `k` spans (MB 18484), as `(area [m^2], length [m])`:
/// the tract's own sections, which are the line's cylinders (`TractChart::to_duct`), cut by [`constriction_span`]. The
/// line already carries their Poiseuille resistance, so it comes out of the constriction's K3
/// ([`GlottalIsland::with_constriction_spans`]).
fn line_span(sections: &[TractSection], k: usize) -> Vec<(f64, f64)> {
    let pairs: Vec<(f64, f64)> = sections.iter().map(|s| (s.area_m2, s.length_m)).collect();
    constriction_span(&pairs, k).to_vec()
}

/// `fricative_noise`'s `k` on a moving tract's shape `shape` (MB 16891, ruling (c) on NavyTower 16887): at the card's
/// section count, `k` itself (the morph's correspondence for equal counts is section by section, so every equal-count
/// card stays bit for bit); at any other count, the shape's section at the card's section `k`'s centre as a fraction of
/// the tract's length from the glottis, the correspondence [`TractChart::log_area_morph`] reads unequal counts by.
fn shape_constriction(card: &[TractSection], k: usize, shape: &[TractSection]) -> usize {
    if shape.len() == card.len() {
        k
    } else {
        section_at_fraction(shape, centre_fraction(card, k))
    }
}

/// The distance of section `k`'s centre from the glottis as a fraction of the tract's length.
fn centre_fraction(sections: &[TractSection], k: usize) -> f64 {
    let total: f64 = sections.iter().map(|s| s.length_m).sum();
    let before: f64 = sections[..k].iter().map(|s| s.length_m).sum();
    (before + 0.5 * sections[k].length_m) / total
}

/// The section whose span holds the distance `fraction` of the tract's length from the glottis: the correspondence
/// [`TractChart::log_area_morph`] reads unequal section counts by (`fricative_noise`'s `k` on a moving tract, NavyTower
/// 16887 (c)). The last section if `fraction` is at or past the lips.
fn section_at_fraction(sections: &[TractSection], fraction: f64) -> usize {
    let at = fraction * sections.iter().map(|s| s.length_m).sum::<f64>();
    let mut end = 0.0;
    for (i, s) in sections.iter().enumerate() {
        end += s.length_m;
        if at < end {
            return i;
        }
    }
    sections.len() - 1
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

/// Check the keyframes and movements (module doc).
fn admit_motion(
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
        let mut glottal_orifice = None;
        let mut subglottal_tube: Vec<TractSection> = Vec::new();
        let mut fricative_noise = None;
        let mut fricative_throat = None;
        let mut free_field_receiver = None;
        let mut rest_gap: Vec<(f64, f64)> = Vec::new();
        let mut fold_tension: Vec<FoldTensionPoint> = Vec::new();
        let mut nasal_draft = nasal::NasalDraft::default();
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
                "termination" if rest.first() == Some(&"unflanged_thin_wall") => {
                    // Levine and Schwinger's exact load; the quadrature's node count is the card's, never a default.
                    arity(&rest, 2, key)?;
                    let nodes = whole(&rest, 1, "unflanged_thin_wall quadrature nodes")?;
                    let nodes = usize::try_from(nodes).map_err(|_| "unflanged_thin_wall quadrature nodes out of range")?;
                    if !(fs_duct::levine_schwinger::LS_MIN_NODES..=fs_duct::levine_schwinger::LS_MAX_NODES).contains(&nodes) {
                        return Err(format!(
                            "unflanged_thin_wall quadrature nodes must be in {}..={}",
                            fs_duct::levine_schwinger::LS_MIN_NODES,
                            fs_duct::levine_schwinger::LS_MAX_NODES
                        ));
                    }
                    once(&mut termination, Termination::UnflangedThinWall { quadrature: LsQuadrature { nodes } }, key)?;
                }
                "termination" => {
                    arity(&rest, 1, key)?;
                    let load = match rest[0] {
                        "flanged" => Termination::FlangedOpen,
                        "unflanged" => Termination::UnflangedOpen,
                        "ideal_open" => Termination::IdealOpen,
                        _ => {
                            return Err(
                                "termination must be flanged, unflanged, ideal_open or unflanged_thin_wall <nodes>".into(),
                            );
                        }
                    };
                    once(&mut termination, load, key)?;
                }
                "glottal_orifice" => {
                    arity(&rest, 1, key)?;
                    let area = field(&rest, 0, "glottal_orifice area")?;
                    if area <= 0.0 {
                        return Err("glottal_orifice area must be positive".into());
                    }
                    once(&mut glottal_orifice, area, key)?;
                }
                "subglottal_tube" => {
                    arity(&rest, 2, key)?;
                    if subglottal_tube.len() >= MAX_SECTIONS {
                        return Err("too many subglottal_tube rows".into());
                    }
                    let section = TractSection {
                        area_m2: field(&rest, 0, "subglottal_tube area")?,
                        length_m: field(&rest, 1, "subglottal_tube length")?,
                    };
                    if !(section.area_m2 > 0.0 && section.length_m > 0.0) {
                        return Err("subglottal_tube area and length must be positive".into());
                    }
                    subglottal_tube.push(section);
                }
                "fricative_noise" => {
                    // The card names the constriction, the section it opens into and, optionally, the printed law and
                    // its obstacle; no gain, no gate.
                    let law = match rest.len() {
                        3 => FricativeLaw::OldhamUkpoho,
                        4 | 5 if rest[3] == "stevens" => FricativeLaw::StevensShadle,
                        4 | 5 => return Err(format!("fricative_noise law must be `stevens`, got `{}`", rest[3])),
                        n => {
                            return Err(format!(
                                "`{key}` takes 3 fields, or 4 with the law `stevens`, or 5 with its obstacle section m, got {n}"
                            ))
                        }
                    };
                    let seed = whole(&rest, 0, "fricative_noise seed")?;
                    let index = |i: usize, what: &str| -> Result<usize, String> {
                        usize::try_from(whole(&rest, i, what)?).map_err(|_| format!("{what} out of range"))
                    };
                    let constriction = index(1, "fricative_noise k")?;
                    let card = FricativeNoiseCard {
                        seed,
                        constriction,
                        expansion: index(2, "fricative_noise j")?,
                        law,
                        obstacle: if rest.len() == 5 { index(4, "fricative_noise m")? } else { constriction },
                    };
                    once(&mut fricative_noise, card, key)?;
                }
                "fricative_throat" => {
                    arity(&rest, 2, key)?;
                    let throat =
                        FricativeThroat { area_m2: field(&rest, 0, "fricative_throat area")?, source: rest[1].to_string() };
                    once(&mut fricative_throat, throat, key)?;
                }
                "receiver" => {
                    arity(&rest, 1, key)?;
                    if rest[0] != "free_field_monopole" {
                        return Err("receiver must be `free_field_monopole`".into());
                    }
                    once(&mut free_field_receiver, (), key)?;
                }
                "rest_gap" => {
                    arity(&rest, 2, key)?;
                    if rest_gap.len() >= MAX_SCHEDULE_POINTS {
                        return Err("too many rest_gap points".into());
                    }
                    let point = (field(&rest, 0, "rest_gap time")?, field(&rest, 1, "rest_gap gap")?);
                    if point.0 < 0.0 || point.1 <= 0.0 {
                        return Err("rest_gap time must be non-negative and its gap positive".into());
                    }
                    if rest_gap.last().is_some_and(|last: &(f64, f64)| point.0 <= last.0) {
                        return Err("rest_gap times must be strictly increasing".into());
                    }
                    rest_gap.push(point);
                }
                "fold_tension" => {
                    arity(&rest, 4, key)?;
                    if fold_tension.len() >= MAX_SCHEDULE_POINTS {
                        return Err("too many fold_tension points".into());
                    }
                    let point = FoldTensionPoint {
                        time_s: field(&rest, 0, "fold_tension time")?,
                        stiffness: field(&rest, 1, "fold_tension stiffness scale")?,
                        vibrato: field(&rest, 2, "fold_tension vibrato scale")?,
                        flutter: field(&rest, 3, "fold_tension flutter scale")?,
                    };
                    if point.time_s < 0.0 {
                        return Err("fold_tension time must be non-negative".into());
                    }
                    // The fold's admitted tension range (`GlottalIsland::set_tension`, 1/4..=4): a factor outside it
                    // could never be set, so it is refused here, by the row's name.
                    if !(0.25..=4.0).contains(&point.stiffness) {
                        return Err("fold_tension stiffness scale must be in 0.25..=4, the fold's tension range".into());
                    }
                    if point.vibrato < 0.0 || point.flutter < 0.0 {
                        return Err("fold_tension vibrato and flutter scales must be non-negative".into());
                    }
                    if fold_tension.last().is_some_and(|last: &FoldTensionPoint| point.time_s <= last.time_s) {
                        return Err("fold_tension times must be strictly increasing".into());
                    }
                    fold_tension.push(point);
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
                nasal_key if nasal::NASAL_KEYS.contains(&nasal_key) => nasal_draft.row(nasal_key, &rest)?,
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
        if body_cover && !rest_gap.is_empty() {
            return Err("rest_gap needs a two_mass fold, not body_cover".into());
        }
        if !two_mass && !rest_gap.is_empty() {
            return Err("rest_gap needs a two_mass fold".into());
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
        // `fold_tension` (module doc): a factor on a gesture needs that gesture; the scaled vibrato keeps its bound.
        if tension_modulation.is_none() && fold_tension.iter().any(|p| p.vibrato != 1.0) {
            return Err("fold_tension's vibrato scale needs tension_modulation: it scales that row's amplitude".into());
        }
        if tension_flutter.is_none() && fold_tension.iter().any(|p| p.flutter != 1.0) {
            return Err("fold_tension's flutter scale needs tension_flutter: it scales that row's output".into());
        }
        if let Some((_, amplitude)) = tension_modulation
            && fold_tension.iter().any(|p| p.vibrato * amplitude > 0.5)
        {
            return Err("fold_tension's vibrato scale takes the tension modulation amplitude above 0.5".into());
        }
        // The flutter rows' own bound (sd 0.1, `TensionFlutterCard::validate`) holds for the scaled rows too.
        for card in [tension_flutter, tension_flutter_second].into_iter().flatten() {
            if fold_tension.iter().any(|p| p.flutter * card.sd_rel > 0.1) {
                return Err("fold_tension's flutter scale takes a tension_flutter sd above 0.1".into());
            }
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
            if head_load.is_none() && !matches!(termination, Some(Termination::FlangedOpen)) {
                return Err(TRACT_GAS_NEEDS_FLANGED.into());
            }
            if tract.state().pressure != gas.state().pressure {
                return Err("tract_gas and ambient need one static pressure: mean flow is not modeled".into());
            }
        }
        // The flow rig's rows (module doc): each refusal names the row it meets.
        if glottal_orifice.is_some() {
            let fold_rows = [
                ("fold_shape", fold_shape.is_some()),
                ("rest_gap", !rest_gap.is_empty()),
                ("fold_tension", !fold_tension.is_empty()),
                ("tension_modulation", tension_modulation.is_some()),
                ("tension_flutter", tension_flutter.is_some()),
                ("pulse_shimmer", pulse_shimmer.is_some()),
                ("aspiration_noise", aspiration_noise.is_some()),
                ("jet_noise", jet_noise.is_some()),
                ("constriction_noise", constriction_noise.is_some()),
            ];
            if let Some((row, _)) = fold_rows.iter().find(|(_, present)| *present) {
                return Err(format!("glottal_orifice replaces the folds: `{row}` moves or reads a fold and is refused with it"));
            }
        }
        if !subglottal_tube.is_empty() && subglottal_airway.is_some() {
            return Err("subglottal_tube and subglottal_airway are two airways below one glottis: refused together".into());
        }
        if let Some(card) = fricative_noise {
            let (k, j) = (card.constriction, card.expansion);
            if !(k < j && j < sections.len()) {
                return Err(format!(
                    "fricative_noise needs k < j < the section count {}: got k {k}, j {j}",
                    sections.len()
                ));
            }
            // A fixed tract has one shape, and the Borda-Carnot loss needs it to be an expansion after `k`.
            if tract_moves.is_empty() && !(sections[j].area_m2 > sections[k].area_m2) {
                return Err("fricative_noise: the Borda-Carnot loss needs section j wider than section k".into());
            }
            // A moving tract (MB 16735 piece 2; NavyTower 16887): every shape takes its own `A_j`, its largest section
            // after its own `k` (MB 15828), so the card's `j` must be that section of the card's own shape. A shape with
            // no expansion after `k`, the card's own included (the word's card is its vowel), takes no constriction and
            // its law drives no band.
            if !tract_moves.is_empty() && sections[j].area_m2.to_bits() != widest_after(&sections, k).to_bits() {
                return Err("fricative_noise with tract_move needs j to be the largest section between k and the lips: \
                            every shape takes its own"
                    .into());
            }
            // The obstacle (`stevens m`, MB 17712): at or after the constriction, with a section after it for the
            // junction. On a moving tract each further shape takes its own `m` (`shape_taps`, MB 25469).
            let m = card.obstacle;
            if !(k <= m && m + 1 < sections.len()) {
                return Err(format!(
                    "fricative_noise needs k <= m < the last section {}: got k {k}, m {m}",
                    sections.len() - 1
                ));
            }
        }
        // `fricative_throat` (MB 25440 ruling A, 25482): the measured throat of the card's constriction, for its laws only.
        if let Some(throat) = &fricative_throat {
            let Some(card) = fricative_noise else {
                return Err("fricative_throat needs fricative_noise: it is that row's constriction's throat".into());
            };
            if !tract_moves.is_empty() {
                return Err("fricative_throat is refused with tract_move: a moving shape has no measured throat".into());
            }
            let (area, a_k, a_j) =
                (throat.area_m2, sections[card.constriction].area_m2, sections[card.expansion].area_m2);
            if !(area.is_finite() && area > 0.0) {
                return Err("fricative_throat area must be finite and positive".into());
            }
            if !(area < a_j) {
                return Err("fricative_throat: the Borda-Carnot loss needs the throat narrower than section j".into());
            }
            if !(area <= 1.5 * a_k) {
                return Err(format!(
                    "fricative_throat {area:e} m^2 is more than 1.5 times section k ({a_k:e} m^2): the sanity bound on a \
                     measured throat"
                ));
            }
        }
        if free_field_receiver.is_some() && (head_sphere.is_some() || observer_law.is_some()) {
            return Err("receiver free_field_monopole is refused with head_sphere and observer_law: both observe the baffle".into());
        }
        // The nasal rows (the `nasal` child module): admitted against the oral sections, then refused by name with each
        // row whose lines or load are one series mouth's. `head_sphere`, `head_load`, `torso_sphere` and `tract_gas` are
        // taken per opening (N3b; the `nasal` module doc).
        let nasal = nasal_draft.admit(&sections)?;
        if nasal.is_some() && tract_gas.is_some() && head_load.is_none() {
            return Err(TRACT_GAS_NASAL_NEEDS_HEAD_LOAD.into());
        }
        if nasal.is_some() {
            let series_rows = [
                ("tract_move", !tract_moves.is_empty()),
                ("fricative_noise", fricative_noise.is_some()),
                ("receiver free_field_monopole", free_field_receiver.is_some()),
            ];
            if let Some((row, _)) = series_rows.iter().find(|(_, present)| *present) {
                return Err(format!("the nasal rows are refused with `{row}`: its lines or its load are one series mouth's"));
            }
        }
        let (keyframes, tract_moves) = admit_motion(
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
            glottal_orifice,
            subglottal_tube,
            fricative_noise,
            fricative_throat,
            free_field_receiver: free_field_receiver.is_some(),
            rest_gap,
            fold_tension,
            nasal,
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

    /// The constriction laws' `A_k` [m^2] for `card` (MB 25440 ruling A): the `fricative_throat` row's measured throat,
    /// else section `k`'s own area, the value every card took before the row existed.
    fn constriction_law_area_m2(&self, card: FricativeNoiseCard) -> f64 {
        self.fricative_throat.as_ref().map_or(self.sections[card.constriction].area_m2, |throat| throat.area_m2)
    }

    /// The `tract_gas` row's record with its source receipts; `None` when the tract holds the exterior's gas.
    #[must_use]
    pub const fn tract_gas_record(&self) -> Option<&GasRecord> {
        self.tract_gas.as_ref()
    }

    /// `tract_gas` without `head_load` on the series tract with `termination flanged`: the lips' declared flanged load
    /// across the gas interface (radiation in the exterior, the eliminated duct modes in the tract's gas), tabulated on
    /// the head table's bins (module doc; MB 23833, 24424). `None` for every other card, which keeps its load exactly
    /// as before; `build` refuses the combinations that need a table and have none.
    fn exterior_flanged_load(&self, duct: &fs_duct::Duct, mechanics_rate: u32)
        -> Result<Option<fs_duct::TabulatedLoad>, String> {
        if self.tract_gas.is_none() || self.head_load || self.nasal.is_some()
            || !matches!(self.termination, Termination::FlangedOpen)
        {
            return Ok(None);
        }
        flanged_mouth_load(duct, self.tract_gas(), self.ambient_gas(), mechanics_rate, GLOTTAL_LINE_MAX_TAPS)
            .map(Some)
            .map_err(|e| format!("flanged load refused: {e}"))
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
        // microphone), its flanged term's duct modes in the tract's gas. Without the row each takes the declared
        // termination exactly as before.
        let head_load = match (self.head_load, &self.head_sphere) {
            (false, _) => None,
            (true, None) => return Err(HEAD_LOAD_NEEDS_SPHERE.into()),
            (true, Some(_)) if !matches!(self.termination, Termination::FlangedOpen) => {
                return Err(HEAD_LOAD_NEEDS_FLANGED.into());
            }
            // The nasal rows: each opening's own table is made with the network's lines (the `nasal` child module).
            (true, Some(_)) if self.nasal.is_some() => None,
            (true, Some(card)) => Some(
                match &self.torso_sphere {
                    None => head_mouth_load(card.radius_m, &duct, gas, exterior, mechanics_rate, GLOTTAL_LINE_MAX_TAPS),
                    Some(t) => head_torso_mouth_load(card.radius_m, t.torso, self.microphone.maximum_frequency_hz, &duct,
                        gas, exterior, mechanics_rate, GLOTTAL_LINE_MAX_TAPS),
                }
                .map_err(|e| format!("head load refused: {e}"))?,
            ),
        };
        // `tract_gas` without `head_load`: the declared flanged load across the gas interface (module doc). Without
        // `tract_gas` this is `None` and the lips' table is the head's, as before.
        let flanged = self.exterior_flanged_load(&duct, mechanics_rate)?;
        let lips_table = head_load.as_ref().or(flanged.as_ref());
        let mouth_load = lips_table.map_or(MouthLoad::Analytic(self.termination), MouthLoad::Tabulated);
        // With `tract_gas` both tract lines check the table against the exterior it was computed in.
        let split = match (&self.tract_gas, lips_table) {
            (None, _) => None,
            (Some(_), Some(table)) => Some(table),
            // The nasal rows with `head_load`: the network's ports read their own tables in the exterior.
            (Some(_), None) if self.head_load && self.nasal.is_some() => None,
            (Some(_), None) if self.nasal.is_some() => return Err(TRACT_GAS_NASAL_NEEDS_HEAD_LOAD.into()),
            (Some(_), None) => return Err(TRACT_GAS_NEEDS_FLANGED.into()),
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
            // `subglottal_tube` (refused with `subglottal_airway` at parse): the card's cylinders through the same seam.
            None if !self.subglottal_tube.is_empty() => {
                let sections = self
                    .subglottal_tube
                    .iter()
                    .map(|s| AirwaySection {
                        length_m: s.length_m,
                        radius_m: (s.area_m2 / core::f64::consts::PI).sqrt(),
                        count: 1,
                        wall: self.wall,
                    })
                    .collect();
                Some(
                    AirwayChain::new(sections, LossModel::Bessel)
                        .and_then(|chain| {
                            let ladder = LineLadder { max_taps: GLOTTAL_LINE_MAX_TAPS, ..LineLadder::default() };
                            UpstreamAirway::realize(&chain, gas, mechanics_rate, &ladder)
                        })
                        .map_err(|e| format!("subglottal tube refused: {e}"))?,
                )
            }
            None => None,
        };
        // The nasal rows: the network's island line and lip flow replace the series duct's, and the nostrils' flow lines
        // join them (the `nasal` child module; every edge in the tract's gas, the exterior the ambient, and with
        // `head_load` each opening's own table, so `split` is `None` there). Without the rows every match below takes
        // the arm it always took.
        let (nasal_line, nasal_lips, nasal_nostrils) = match &self.nasal {
            Some(card) => {
                let head = match (self.head_load, &self.head_sphere) {
                    (true, Some(sphere)) => Some(nasal::NasalHead {
                        radius_m: sphere.radius_m,
                        torso: self.torso_sphere.as_ref().map(|t| (t.torso, self.microphone.maximum_frequency_hz)),
                    }),
                    _ => None,
                };
                let lines = nasal::nasal_lines(
                    &self.sections,
                    card,
                    gas,
                    self.wall.as_ref(),
                    nasal::Exterior { gas: exterior, termination: self.termination, head },
                    mechanics_rate,
                    nasal::Velum::Open,
                )?;
                (Some((lines.line, lines.inlet_radius_m)), Some(lines.lips), lines.nostrils.zip(Some(lines.nostril_radius_m)))
            }
            None => (None, None, None),
        };
        let island = match (split, nasal_line) {
            (_, Some((line, inlet_radius_m))) => {
                GlottalIsland::new_with_tract_line(self.fold, self.two_mass, line, inlet_radius_m, gas, mechanics_rate)
            }
            (None, None) => GlottalIsland::new_with_wall(
                self.fold,
                self.two_mass,
                &duct,
                gas,
                mouth_load,
                mechanics_rate,
                self.wall.as_ref(),
            ),
            (Some(table), None) => GlottalIsland::new_with_wall_exterior(
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
        .and_then(|island| match self.glottal_orifice {
            Some(area) => island.with_orifice(area),
            None => Ok(island),
        })
        // A moving tract takes every shape's constriction below (with_constriction_keyframes), shape 0's among them.
        // K3 is less the line's own Poiseuille resistance over the sections the constriction spans (MB 18484). `A_k` is
        // the `fricative_throat` row's measured throat when the card has one (ruling A); the span stays the line's.
        .and_then(|island| match self.fricative_noise {
            Some(card) if self.tract_moves.is_empty() => {
                let span = line_span(chart.sections(), card.constriction);
                island
                    .with_constriction(self.constriction_law_area_m2(card), self.sections[card.expansion].area_m2)
                    .and_then(|island| island.with_constriction_spans(&[span.as_slice()]))
            }
            _ => Ok(island),
        })
        .map_err(|e| format!("glottal island refused: {e:?}"))?;
        let (mouth, receipt) = match (split, nasal_lips) {
            (_, Some(lips)) => Ok(lips),
            (None, None) => mouth_flow_line(
                &duct,
                gas,
                mouth_load,
                mechanics_rate,
                island.port_impedance(),
                self.wall.as_ref(),
                &LineLadder::default(),
            ),
            (Some(table), None) => mouth_flow_line_exterior(
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
        let (island, mut motion) = if self.tract_moves.is_empty() {
            (island, TractMotion::default())
        } else {
            let mut motion = self.realize_motion(&chart, mechanics_rate, workers)?;
            mouth = mouth
                .to_moving(core::mem::take(&mut motion.mouth_taps))
                .map_err(|e| format!("moving mouth-flow line refused: {e}"))?;
            let island = island
                .with_tract_keyframes(core::mem::take(&mut motion.island_taps))
                .map_err(|e| format!("moving tract line refused: {e:?}"))?;
            // `fricative_noise` (MB 16735 piece 2): every shape's constriction, the card's own first, switched with the
            // tract line; a shape with no expansion after its `k` takes none.
            let island = match self.fricative_noise {
                Some(card) => {
                    // Shape 0's `A_k` by the same rule as the source's shape 0 (`fricative_throat` is refused with
                    // `tract_move` at parse, so this is section `k`'s own area).
                    let mut shapes =
                        vec![(self.constriction_law_area_m2(card), self.sections[card.expansion].area_m2)];
                    shapes.extend(motion.fricative.iter().map(|f| (f.area_k_m2, f.area_j_m2)));
                    // Each shape's K3 less the line's own Poiseuille resistance over that shape's span (MB 18484).
                    let first = line_span(chart.sections(), card.constriction);
                    let mut spans = vec![first.as_slice()];
                    spans.extend(motion.fricative.iter().map(|f| f.span.as_slice()));
                    island
                        .with_constriction_keyframes(&shapes)
                        .and_then(|island| island.with_constriction_spans(&spans))
                }
                None => Ok(island),
            }
            .map_err(|e| format!("moving constriction refused: {e:?}"))?;
            (island, motion)
        };
        let mouth_area = self.sections.last().map_or(0.0, |s| s.area_m2);
        let mouth_radius_m = (mouth_area / core::f64::consts::PI).sqrt();
        let medium = RayleighMedium {
            density: exterior.density,
            sound_speed: exterior.sound_speed,
        };
        // The openings on the baffle: the mouth at the origin, and with the nasal rows the left nostril at (x, y) and the
        // right at (-x, y), in the order the trace carries their flows.
        let mut openings = vec![CircularOpening {
            center_m: [0.0, 0.0], radius_m: mouth_radius_m,
        }];
        if let (Some(card), Some((_, radius))) = (&self.nasal, &nasal_nostrils) {
            openings.push(CircularOpening { center_m: [card.nostril_x_m, card.nostril_y_m], radius_m: radius[0] });
            openings.push(CircularOpening { center_m: [-card.nostril_x_m, card.nostril_y_m], radius_m: radius[1] });
        }
        // The nasal rows with `head_sphere` (N3b): each opening's own receiver and head FIR, summed after the FIRs (the
        // `nasal` module doc); the mouth's pair takes the fields every card has. Otherwise each arm is the one it was.
        let (receiver, observer_paths, head, nostril_heads) = match (&self.head_sphere, openings.len() > 1) {
            (Some(card), true) => {
                let heads = nasal::opening_heads(
                    &openings,
                    card.radius_m,
                    self.torso_sphere.as_ref().map(|t| t.torso),
                    mechanics_rate,
                    self.output_rate_hz,
                    self.microphone,
                    exterior,
                    observer.as_ref(),
                )?;
                let (receiver, head) = heads.mouth;
                (receiver, heads.observer_paths, Some(head), Some(heads.nostrils))
            }
            _ => {
                let (receiver, observer_paths) = match observer.as_ref() {
                    Some(propagation) => {
                        let (receiver, report) = BaffledPressure::circular_openings_multirate_gas_differentiating(
                            &openings, mechanics_rate, self.output_rate_hz, self.microphone,
                            exterior, propagation, &mut || false,
                        ).map_err(|error| format!("voice observer paths refused: {error}"))?;
                        (receiver, Some(report))
                    }
                    // `receiver free_field_monopole` (refused with observer_law at parse).
                    None if self.free_field_receiver => (
                        free_field_monopole(mouth_area, mechanics_rate, self.output_rate_hz, self.microphone, medium)?,
                        None,
                    ),
                    None if openings.len() > 1 => (
                        BaffledPressure::circular_openings_multirate_differentiating(
                            &openings, mechanics_rate, self.output_rate_hz,
                            self.microphone, medium,
                        )?,
                        None,
                    ),
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
                            None => HeadSphere::new(card.radius_m, mouth_radius_m, self.microphone, medium,
                                self.output_rate_hz),
                            Some(t) => HeadSphere::with_torso(card.radius_m, mouth_radius_m, self.microphone, medium,
                                self.output_rate_hz, t.torso),
                        }
                        .map_err(|e| format!("head sphere refused: {e}"))?,
                    ),
                    None => None,
                };
                (receiver, observer_paths, head, None)
            }
        };
        let passive_formants_hz = if self.nasal.is_some() {
            Err("the passive formants are the series tract's: a nasal card does not read them".to_string())
        } else {
            match lips_table {
                Some(table) => chart.tmm_formants_tabulated(gas, self.wall.as_ref(), table, 4),
                None => chart.tmm_formants_terminated(gas, self.wall.as_ref(), self.termination, 4),
            }
            .map_err(|e| format!("{e:?}"))
        };
        let nostril_flow_taps = nasal_nostrils.as_ref().map(|(lines, _)| [lines[0].1.n_fft, lines[1].1.n_fft]);
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
            nostril_flow_taps,
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
        // `fricative_noise`: the source after section k, referenced to the island's port, with the island's own mouth
        // load: the declared termination, or with `head_load` the same table, read in the exterior with `tract_gas`
        // (module doc). Shape 0 is the card's own tract; a moving tract adds its further shapes' lines and constrictions
        // (`shape_taps`), switched with the island's.
        let fricative = match self.fricative_noise {
            Some(card) => {
                // `A_k` the `fricative_throat` row's throat when the card has one (ruling A); `A_j` and `A_m` the tract's.
                let (area_k, area_j) = (self.constriction_law_area_m2(card), self.sections[card.expansion].area_m2);
                let noise = ConstrictionNoise::new_fricative(
                    ConstrictionNoiseCard { seed: card.seed, duct_radius_m: (area_j / core::f64::consts::PI).sqrt() },
                    1.0 / f64::from(mechanics_rate),
                )
                .map_err(str::to_string)?;
                // The junction lines after the obstacle section, `k` itself unless `stevens m` names another.
                let lines = match split {
                    None => junction_source_lines(
                        &duct,
                        card.obstacle,
                        gas,
                        mouth_load,
                        mechanics_rate,
                        island.port_impedance(),
                        self.wall.as_ref(),
                        &LineLadder::default(),
                    ),
                    Some(table) => junction_source_lines_exterior(
                        &duct,
                        card.obstacle,
                        gas,
                        exterior,
                        table,
                        mechanics_rate,
                        island.port_impedance(),
                        self.wall.as_ref(),
                        &LineLadder::default(),
                    ),
                }
                .map_err(|e| format!("fricative junction lines refused: {e}"))?;
                let mut to_mouth_flow = ShapeLine::Fixed(lines.to_mouth_flow);
                let mut to_inlet_wave = ShapeLine::Fixed(lines.to_inlet_wave);
                let mut shapes = vec![FricativeShape {
                    section: card.constriction,
                    area_k_m2: area_k,
                    area_m_m2: self.sections[card.obstacle].area_m2,
                    open_area_ratio: area_k / area_j,
                    duct_radius_m: noise.card().duct_radius_m,
                }];
                if !self.tract_moves.is_empty() {
                    let further = core::mem::take(&mut motion.fricative);
                    shapes.extend(further.iter().map(|f| FricativeShape {
                        section: f.section,
                        area_k_m2: f.area_k_m2,
                        area_m_m2: f.area_m_m2,
                        open_area_ratio: f.area_k_m2 / f.area_j_m2,
                        duct_radius_m: (f.area_j_m2 / core::f64::consts::PI).sqrt(),
                    }));
                    let (mouth_taps, inlet_taps): (Vec<_>, Vec<_>) =
                        further.into_iter().map(|f| (f.mouth_taps, f.inlet_taps)).unzip();
                    to_mouth_flow = to_mouth_flow
                        .to_moving(mouth_taps)
                        .map_err(|e| format!("moving fricative mouth line refused: {e}"))?;
                    to_inlet_wave = to_inlet_wave
                        .to_moving(inlet_taps)
                        .map_err(|e| format!("moving fricative inlet line refused: {e}"))?;
                }
                Some(FricativeSource {
                    noise,
                    to_mouth_flow,
                    to_inlet_wave,
                    shapes,
                    shape: 0,
                    rig: self.glottal_orifice.is_some(),
                    law: card.law,
                    last_flow_m3_s: 0.0,
                })
            }
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
            fricative,
            pulse_shimmer: self.pulse_shimmer.map(PulseShimmer::new),
            glottal_trace: None,
            gas_density: gas.density,
            gas_viscosity: gas.dynamic_viscosity,
            fold_length_m: self.fold.fold_length_m,
            schedule: Schedule {
                points: self.subglottal.clone(),
                cursor: 0,
            },
            rest_gap: (!self.rest_gap.is_empty()).then(|| Schedule {
                points: self.rest_gap.clone(),
                cursor: 0,
            }),
            fold_tension: (!self.fold_tension.is_empty()).then(|| {
                let factor = |read: fn(&FoldTensionPoint) -> f64| Schedule {
                    points: self.fold_tension.iter().map(|p| (p.time_s, read(p))).collect(),
                    cursor: 0,
                };
                [factor(|p| p.stiffness), factor(|p| p.vibrato), factor(|p| p.flutter)]
            }),
            tension_modulation: self.tension_modulation,
            tension_flutter: self.tension_flutter.map(|card| TensionFlutter::new(card, self.output_rate_hz)),
            tension_flutter_second: self
                .tension_flutter_second
                .map(|card| TensionFlutter::on_tile(card, self.output_rate_hz, 1)),
            trace: vec![0.0; self.mechanics_per_output as usize * openings.len()],
            nostrils: nasal_nostrils.map(|([(left, _), (right, _)], _)| [left, right]),
            nostril_heads,
            ratio: self.mechanics_per_output as usize,
            mechanics_rate,
            info,
            steps: 0,
            completed: 0,
            poisoned: false,
            peak_flow_m3_s: 0.0,
            min_flow_m3_s: 0.0,
            closed_steps: 0,
            rig_tail: self.glottal_orifice.map(|_| {
                let total = self.samples * u64::from(self.mechanics_per_output);
                (total.saturating_sub(2 * u64::from(mechanics_rate)), 0.0, 0)
            }),
            tract_switches: motion.switches,
            next_switch: 0,
        })
    }

    /// Every further shape of the moving tract, in the order the render reaches them, with the step each starts at
    /// (module doc). A shape is the morph at `k / steps` from the shape the tract holds to the move's target, its first
    /// section among the rest (MB 25468); `k = steps` is the target as written.
    ///
    /// The shapes and switch steps come first, in one serial pass; each shape's lines are then built on up to `workers`
    /// threads ([`build_shapes`]). A refusal is the one a shape-by-shape build would meet first: a build refusal of an
    /// earlier shape before a later shape's morph or step refusal.
    fn realize_motion(&self, chart: &TractChart, mechanics_rate: u32, workers: usize) -> Result<TractMotion, String> {
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
                    .map_err(|e| format!("moving tract morph refused: {e:?}"));
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
        let built = build_shapes(shapes.len(), workers, |i| self.shape_taps(&shapes[i], mechanics_rate))?;
        if let Some(e) = cut {
            return Err(e);
        }
        let mut formants = Vec::with_capacity(built.len());
        for (island_taps, mouth_taps, f, fricative) in built {
            motion.island_taps.push(island_taps);
            motion.mouth_taps.push(mouth_taps);
            formants.push(Some(f));
            motion.fricative.extend(fricative);
        }
        for (at, name) in targets {
            let f = formants[at].take().ok_or("a moving tract target was reached twice")?;
            motion.keyframe_formants_hz.push((name, f));
        }
        Ok(motion)
    }

    /// The island's line taps with the shape's own `rho c / S` and first section, the mouth-flow line taps and the
    /// passive formants of one further shape, realized with the card's mouth load, wall and gases exactly as `build`
    /// realizes the card's own shape. Every line of the shape is referenced to the shape's own first section (MB 25468).
    #[allow(clippy::type_complexity)]
    fn shape_taps(
        &self,
        chart: &TractChart,
        mechanics_rate: u32,
    ) -> Result<ShapeLines, String> {
        let gas = self.tract_gas();
        let exterior = self.ambient_gas();
        let duct = chart.to_duct();
        // `build` has already refused every other combination of these rows.
        let head_load = match (self.head_load, &self.head_sphere) {
            (true, Some(card)) => Some(
                match &self.torso_sphere {
                    None => head_mouth_load(card.radius_m, &duct, gas, exterior, mechanics_rate, GLOTTAL_LINE_MAX_TAPS),
                    Some(t) => head_torso_mouth_load(card.radius_m, t.torso, self.microphone.maximum_frequency_hz,
                        &duct, gas, exterior, mechanics_rate, GLOTTAL_LINE_MAX_TAPS),
                }
                .map_err(|e| format!("moving tract head load refused: {e}"))?,
            ),
            _ => None,
        };
        let flanged = self.exterior_flanged_load(&duct, mechanics_rate).map_err(|e| format!("moving tract {e}"))?;
        let lips_table = head_load.as_ref().or(flanged.as_ref());
        let mouth_load = lips_table.map_or(MouthLoad::Analytic(self.termination), MouthLoad::Tabulated);
        let split = if self.tract_gas.is_some() { lips_table } else { None };
        let wall = self.wall.as_ref();
        // This shape's own `rho c / S` and first section: the island steps to them at the shape's switch, and its
        // mouth-flow and fricative junction lines below are referenced to this `zc`.
        let (island_taps, zc, inlet_area_m2) = match split {
            None => glottal_tract_line_taps(&duct, gas, None, mouth_load, mechanics_rate, wall),
            Some(table) => {
                glottal_tract_line_taps(&duct, gas, Some(exterior), MouthLoad::Tabulated(table), mechanics_rate, wall)
            }
        }
        .map_err(|e| format!("moving tract line refused: {e:?}"))?;
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
        let formants = match lips_table {
            Some(table) => chart.tmm_formants_tabulated(gas, wall, table, 4),
            None => chart.tmm_formants_terminated(gas, wall, self.termination, 4),
        }
        .map_err(|e| format!("{e:?}"));
        // `fricative_noise` (MB 16735 piece 2): this shape's `k` by the morph's own correspondence (`shape_constriction`;
        // MB 16891 (c)), and its obstacle `m` by the same correspondence (MB 25469; `m = k` gives `k` itself); its
        // junction lines after `m` with this shape's own mouth load, exactly as `build` makes the card's; and its
        // constriction, `A_j` the shape's largest section after its `k` (MB 15828). The correspondence is monotone in
        // the centre fraction, so the card's `k <= m` gives the shape's `k <= m`.
        let fricative = match self.fricative_noise {
            Some(card) => {
                let sections = chart.sections();
                let k = shape_constriction(&self.sections, card.constriction, sections);
                if k + 1 >= sections.len() {
                    return Err("fricative_noise: the card's section k falls on a moving shape's last section, with \
                                nothing after it"
                        .into());
                }
                let m = shape_constriction(&self.sections, card.obstacle, sections);
                if m + 1 >= sections.len() {
                    return Err("fricative_noise: the card's obstacle section m falls on a moving shape's last section, \
                                with nothing after it for the junction"
                        .into());
                }
                let lines = match split {
                    None => junction_source_lines(&duct, m, gas, mouth_load, mechanics_rate, zc, wall, &LineLadder::default()),
                    Some(table) => junction_source_lines_exterior(
                        &duct, m, gas, exterior, table, mechanics_rate, zc, wall, &LineLadder::default(),
                    ),
                }
                .map_err(|e| format!("moving fricative junction lines refused: {e}"))?;
                let plain = |line: &DelayedFilter| {
                    line.quiescent_fir()
                        .map(|(_, taps)| taps.to_vec())
                        .ok_or("a moving tract needs every shape's fricative lines to be plain FIRs")
                };
                Some(FricativeShapeLines {
                    mouth_taps: plain(&lines.to_mouth_flow)?,
                    inlet_taps: plain(&lines.to_inlet_wave)?,
                    section: k,
                    area_k_m2: sections[k].area_m2,
                    area_m_m2: sections[m].area_m2,
                    area_j_m2: widest_after(sections, k),
                    span: line_span(sections, k),
                })
            }
            None => None,
        };
        Ok(((island_taps, zc, inlet_area_m2), mouth_taps.to_vec(), formants, fricative))
    }
}

/// One further shape's island taps with its `rho c / S` and first section, mouth-flow taps, passive formants and, with
/// `fricative_noise`, its fricative lines (`VoiceDescription::shape_taps`).
type ShapeLines = ((Vec<f64>, f64, f64), Vec<f64>, Result<Vec<f64>, String>, Option<FricativeShapeLines>);

/// One further shape's `fricative_noise` (MB 16735 piece 2): its junction lines' taps and its constriction.
struct FricativeShapeLines {
    mouth_taps: Vec<f64>,
    inlet_taps: Vec<f64>,
    /// The shape's `k`: its section at the card's section `k`'s fractional distance from the glottis.
    section: usize,
    area_k_m2: f64,
    /// The area of the shape's obstacle section `m` (the card's `m` by the same correspondence; MB 25469), which the
    /// junction lines sit after and the Stevens form's duct-area factor reads.
    area_m_m2: f64,
    /// The shape's largest section after `k` (MB 15828).
    area_j_m2: f64,
    /// The shape's line sections its constriction spans (`line_span`): their Poiseuille resistance comes out of the
    /// shape's K3 (MB 18484).
    span: Vec<(f64, f64)>,
}

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
    /// The island's line taps, `rho c / S` and first section of shapes 1.. (shape 0 is the island's own).
    island_taps: Vec<(Vec<f64>, f64, f64)>,
    /// The mouth-flow line taps of shapes 1...
    mouth_taps: Vec<Vec<f64>>,
    /// With `fricative_noise`, the fricative lines and constriction of shapes 1.. (empty without the row).
    fricative: Vec<FricativeShapeLines>,
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
    /// Converged lengths of the left and right nostrils' flow lines; `None` without the nasal rows.
    pub nostril_flow_taps: Option<[usize; 2]>,
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
    /// With `glottal_orifice` only: the mean glottal flow [m^3/s] over the render's last 2 s (the flow rig's read
    /// window; the whole render when it is shorter). `None` without the row, or before the window has begun.
    pub rig_mean_glottal_flow_m3_s: Option<f64>,
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

/// The columns of [`VoicePerformance::glottal_trace`], in row order, each an `f32`:
/// - `gap_m`: [`GlottalFrame::gap_m`](crate::glottis::GlottalFrame::gap_m) at the end of the step [m]: the narrower
///   section's gap, negative in contact, or with a `fold_shape` the narrower section's mean open gap (never
///   negative). Both [`FoldMode`]s are non-negative at every strip and the two sections share the strips' rest gaps,
///   so one section is the narrower at every strip and `l max(gap_m, 0)` is the glottal area a view from above sees;
/// - `flow_m3_s`: [`GlottalFrame::flow_m3_s`](crate::glottis::GlottalFrame::flow_m3_s), the glottal volume flow
///   [m^3/s];
/// - `x_lower_m`, `x_upper_m`: the lower and upper masses' displacements [m]
///   ([`GlottalIsland::mass_displacements_m`]); with a `fold_shape`, each section's modal coordinate `q`, strip `j`
///   sitting at `rest[j] + mode[j] q` (the sine mode's middle moves `(pi/2) q`); NaN for the one-mass fold.
pub const GLOTTAL_TRACE_COLUMNS: [&str; 4] = ["gap_m", "flow_m3_s", "x_lower_m", "x_upper_m"];

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
    /// Opt-in `rest_gap` rows: the fold's rest gap over time [m] (module doc); `None` keeps the `fold` row's gap.
    rest_gap: Option<Schedule>,
    /// Opt-in `fold_tension` rows: the stiffness, vibrato and flutter factors over time (module doc); `None` keeps the
    /// card's own tension.
    fold_tension: Option<[Schedule; 3]>,
    tension_modulation: Option<(f64, f64)>,
    tension_flutter: Option<TensionFlutter>,
    tension_flutter_second: Option<TensionFlutter>,
    /// Opt-in aspiration source, with the gas properties and fold length its law needs.
    aspiration: Option<GlottalAspiration>,
    /// Opt-in glottal jet-noise source (Verge's fraction on the glottal jet, an extrapolation).
    jet_noise: Option<GlottalJetNoise>,
    /// Opt-in glottal constriction noise (Oldham and Ukpoho's printed law, an extrapolation).
    constriction_noise: Option<ConstrictionNoise>,
    /// Opt-in supraglottal (fricative) source and its junction lines.
    fricative: Option<FricativeSource>,
    /// Opt-in per-cycle gain on the wave the tract receives.
    pulse_shimmer: Option<PulseShimmer>,
    /// The render option `--glottal-trace` ([`Self::keep_glottal_trace`]): per physical substep the
    /// [`GLOTTAL_TRACE_COLUMNS`], appended after the step, so nothing the island, the tract or the receiver sees
    /// changes. `None` keeps nothing.
    glottal_trace: Option<Vec<f32>>,
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
    /// With `glottal_orifice` only, the flow rig's read window: the first physical substep of the render's last 2 s,
    /// the sum of the glottal flow from it [m^3/s] and the substeps summed. `None` without the row.
    rig_tail: Option<(u64, f64, u64)>,
    /// With the nasal rows, the left and right nostrils' flow lines, driven by the island's wave as the mouth's is; the
    /// trace then carries three flows per substep, time-major (mouth, left, right). `None` keeps the one opening.
    nostrils: Option<[DelayedFilter; 2]>,
    /// With the nasal rows and `head_sphere` (N3b): the left and right nostrils' own receivers and head FIRs (the mouth's
    /// are `receiver` and `head`); the trace is then opening-major (the mouth's substeps, the left's, the right's) and
    /// the three FIRs' outputs are summed. `None` keeps the one receiver.
    nostril_heads: Option<[(BaffledPressure, HeadSphere); 2]>,
}

fn refusal(what: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(what))
}

/// The `fricative_noise` row realized: the supraglottal source, its two junction lines, and what its drop needs.
struct FricativeSource {
    noise: ConstrictionNoise,
    /// `push(p_s)`: the source's own mouth flow [m^3/s] (`JunctionSourceLines::to_mouth_flow`), moving with the tract.
    to_mouth_flow: ShapeLine,
    /// `push(p_s)`: the wave [Pa] it returns to the island (`JunctionSourceLines::to_inlet_wave`), moving with the tract.
    to_inlet_wave: ShapeLine,
    /// The constriction of every shape, shape 0 the card's own (one shape on a fixed tract).
    shapes: Vec<FricativeShape>,
    /// The present shape: the last switch's.
    shape: usize,
    /// The card has `glottal_orifice`: the source reads the held flow, not [`Self::last_flow_m3_s`].
    rig: bool,
    /// The card's printed law.
    law: FricativeLaw,
    /// On the folds: the island's flow [m^3/s] of the step before, the flow the source reads (module doc).
    last_flow_m3_s: f64,
}

/// One shape's constriction for the fricative source: its section `k`, `A_k`, `sigma = A_k / A_j` and the radius of
/// `A_j` (`sigma >= 1` on a shape with no expansion after `k`, where the law drives no band).
#[derive(Debug, Clone, Copy)]
struct FricativeShape {
    /// The (c) rule's section `k`, read by the tests; the junction lines were built after the shape's obstacle section
    /// `m`, which is `k` unless `stevens m` names another, and the step reads the rest.
    #[cfg_attr(not(test), allow(dead_code))]
    section: usize,
    area_k_m2: f64,
    /// The area of the duct the source sits in, section `m` (the card's obstacle; on a moving shape, its own `m` by the
    /// card's correspondence, MB 25469): the Stevens form's duct-area factor reads it (MB 18252).
    area_m_m2: f64,
    open_area_ratio: f64,
    duct_radius_m: f64,
}

/// The `receiver free_field_monopole` row (module doc): a point monopole at the mouth centre in free space, `rho dQ/dt /
/// (4 pi R)`, on the baffled receiver's own band-limited differentiator and travel-time kernel. One surface sample's
/// Rayleigh gain is `rho A / (2 pi R)` times its velocity per unit flow; a rigid baffle's image doubles a point source, so
/// the free field is that sample with half the flow basis, `1 / (2 A)`. The band rule is the disk receiver's.
fn free_field_monopole(
    mouth_area_m2: f64,
    mechanics_rate: u32,
    output_rate: u32,
    microphone: CircularOutletReceiver,
    medium: RayleighMedium,
) -> Result<BaffledPressure, String> {
    if !(microphone.maximum_frequency_hz > 0.0 && microphone.maximum_frequency_hz <= 0.45 * f64::from(output_rate)) {
        return Err("free_field_monopole band must be positive and at most 0.45 of its output clock".into());
    }
    if !(mouth_area_m2.is_finite() && mouth_area_m2 > 0.0) {
        return Err("free_field_monopole needs a finite positive mouth area".into());
    }
    let surface = [SurfaceSample { position_m: [0.0; 3], area_m2: mouth_area_m2, mode_shape: vec![0.5 / mouth_area_m2] }];
    BaffledPressure::new_multirate_differentiating(&surface, 1, mechanics_rate, output_rate, microphone.position_m, medium)
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

    /// Keep the glottal trace from the next physical substep on ([`Self::glottal_trace`]), with room reserved for
    /// `substeps` rows. The render option `--glottal-trace` calls it; the card has no key for it, so a card's text and
    /// its pressure are those of a render without it, sample for sample.
    pub fn keep_glottal_trace(&mut self, substeps: usize) {
        self.glottal_trace = Some(Vec::with_capacity(substeps.saturating_mul(GLOTTAL_TRACE_COLUMNS.len())));
    }

    /// The glottal trace kept so far: time-major, one row of [`GLOTTAL_TRACE_COLUMNS`] per physical substep, on the
    /// mechanics clock ([`VoiceInfo::mechanics_rate_hz`]). `None` unless [`Self::keep_glottal_trace`] was called.
    #[must_use]
    pub fn glottal_trace(&self) -> Option<&[f32]> {
        self.glottal_trace.as_deref()
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
            rig_mean_glottal_flow_m3_s: self.rig_tail.and_then(|(_, sum, n)| (n > 0).then_some(sum / n as f64)),
        }
    }

    fn sample(&mut self) -> Result<f64, String> {
        for j in 0..self.ratio {
            let t = self.steps as f64 / f64::from(self.mechanics_rate);
            // The gesture is far slower than the audio clock, so the tension is set once per output sample.
            if j == 0 {
                // `fold_tension`: the stiffness, vibrato and flutter factors at `t` (module doc). Without the rows
                // nothing is read, and every expression below is the one written before the rows existed.
                let scales = self
                    .fold_tension
                    .as_mut()
                    .map(|[stiffness, vibrato, flutter]| [stiffness.at(t), vibrato.at(t), flutter.at(t)]);
                let vibrato = self.tension_modulation.map(|(rate_hz, amplitude)| {
                    let amplitude = scales.map_or(amplitude, |[_, g_v, _]| g_v * amplitude);
                    amplitude * det::sin(core::f64::consts::TAU * rate_hz * t)
                });
                let flutter = self.tension_flutter.as_mut().map(|flutter| {
                    let first = flutter.next_relative();
                    // No sum with zero when the row is absent: such a card keeps the first row's value bit for bit.
                    let y = match self.tension_flutter_second.as_mut() {
                        Some(second) => first + second.next_relative(),
                        None => first,
                    };
                    scales.map_or(y, |[_, _, g_f]| g_f * y)
                });
                // The same: only a card with both gestures takes the sum, and only such a card was refused before it.
                let moved = match (vibrato, flutter) {
                    (Some(v), Some(y)) => Some(v + y),
                    (one, None) | (None, one) => one,
                };
                let tension = match (scales, moved) {
                    (None, moved) => moved.map(|moved| 1.0 + moved),
                    (Some([s, _, _]), Some(moved)) => Some(s * (1.0 + moved)),
                    (Some([s, _, _]), None) => Some(s),
                };
                if let Some(tension) = tension {
                    self.island
                        .set_tension(tension)
                        .map_err(|e| format!("fold tension refused: {e:?}"))?;
                }
                // `rest_gap`: the abducting gesture, on the same clock (module doc). Without the rows nothing is called.
                if let Some(rest_gap) = self.rest_gap.as_mut() {
                    self.island
                        .set_rest_gap(rest_gap.at(t))
                        .map_err(|e| format!("fold rest gap refused: {e:?}"))?;
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
                // The fricative's lines take the same crossfade and its law the new shape's constriction (MB 16735
                // piece 2); the island has already taken that shape's Borda-Carnot factor and K3 in `switch_tract`.
                if let Some(source) = self.fricative.as_mut() {
                    if switch.shape >= source.shapes.len() {
                        return Err("a moving tract switch names a shape the fricative source does not have".into());
                    }
                    source
                        .to_mouth_flow
                        .switch(switch.shape, switch.fade)
                        .map_err(|e| format!("moving fricative mouth line switch refused: {e}"))?;
                    source
                        .to_inlet_wave
                        .switch(switch.shape, switch.fade)
                        .map_err(|e| format!("moving fricative inlet line switch refused: {e}"))?;
                    source.shape = switch.shape;
                }
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
            // `fricative_noise`: the series pressure after section k, from the last step's flow, at the held flow through
            // the orifice and the constriction at the gesture's pressure, or on the folds at the island's flow of the step
            // before (module doc). It is not in the glottal drive: it reaches the island only through its own returned
            // wave, below.
            let fricative_pa = match self.fricative.as_mut() {
                Some(source) => {
                    let flow = if source.rig {
                        self.island.orifice_held_flow_m3_s(gesture).unwrap_or(0.0)
                    } else {
                        source.last_flow_m3_s.max(0.0)
                    };
                    // The present shape's constriction: past `sigma = 1` (section k as wide as every section after it)
                    // the printed law drives no band, so the source turns off through the law as the tract opens.
                    let shape = source.shapes[source.shape];
                    match source.law {
                        FricativeLaw::OldhamUkpoho => {
                            let jet_speed = flow / shape.area_k_m2;
                            let open = 1.0 - shape.open_area_ratio;
                            let dp = 0.5 * self.gas_density * jet_speed * jet_speed * open * open;
                            source.noise.pressure_pa_jet_in(dp, jet_speed, shape.open_area_ratio, shape.duct_radius_m)
                        }
                        // Stevens' law reads the flow and `A_k`, and its duct-area factor the obstacle section's `A_m`; on
                        // a shape with no constriction (`sigma >= 1`) it is silent as the printed law is, every stream
                        // still advancing.
                        FricativeLaw::StevensShadle => {
                            let constricted = shape.open_area_ratio < 1.0;
                            source.noise.pressure_pa_stevens(
                                if constricted { flow } else { 0.0 },
                                shape.area_k_m2,
                                shape.area_m_m2,
                            )
                        }
                    }
                }
                None => 0.0,
            };
            let frame = if self.island.has_upstream() {
                self.island.step_lung(gesture, series)
            } else {
                self.island.step(drive)
            }
            .map_err(|e| format!("glottal island refused: {e:?}"))?;
            if let Some(source) = self.fricative.as_mut() {
                source.noise.record(frame.flow_m3_s);
                source.last_flow_m3_s = frame.flow_m3_s;
                let returned = source
                    .to_inlet_wave
                    .push(fricative_pa)
                    .map_err(|e| format!("fricative inlet line refused: {e:?}"))?;
                self.island
                    .add_returning_wave_pa(returned)
                    .map_err(|e| format!("glottal island refused: {e:?}"))?;
            }
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
            let mut flow = self
                .mouth
                .push(wave)
                .map_err(|e| format!("mouth-flow line refused: {e:?}"))?;
            if let Some(source) = self.fricative.as_mut() {
                // Superposition on the linear tract: the source's own mouth flow joins the island's.
                flow += source
                    .to_mouth_flow
                    .push(fricative_pa)
                    .map_err(|e| format!("fricative mouth line refused: {e:?}"))?;
            }
            match self.nostrils.as_mut() {
                None => self.trace[j] = flow,
                // The nasal rows: the same wave into both nostrils' lines; the receiver takes the flows time-major.
                Some([left, right]) => {
                    let left_flow = left.push(wave).map_err(|e| format!("nostril flow line refused: {e:?}"))?;
                    let right_flow = right.push(wave).map_err(|e| format!("nostril flow line refused: {e:?}"))?;
                    if self.nostril_heads.is_some() {
                        // One receiver per opening: each takes its own substeps (opening-major).
                        let r = self.ratio;
                        self.trace[j] = flow;
                        self.trace[r + j] = left_flow;
                        self.trace[2 * r + j] = right_flow;
                    } else {
                        self.trace[3 * j] = flow;
                        self.trace[3 * j + 1] = left_flow;
                        self.trace[3 * j + 2] = right_flow;
                    }
                }
            }
            self.peak_flow_m3_s = self.peak_flow_m3_s.max(frame.flow_m3_s);
            self.min_flow_m3_s = self.min_flow_m3_s.min(frame.flow_m3_s);
            if frame.shut {
                self.closed_steps += 1;
            }
            if let Some((from, sum, count)) = self.rig_tail.as_mut()
                && self.steps >= *from
            {
                *sum += frame.flow_m3_s;
                *count += 1;
            }
            if let Some(trace) = self.glottal_trace.as_mut() {
                // Read after the step: the island's state is the one the frame reports, and nothing here feeds back.
                let [lower, upper] = self.island.mass_displacements_m().unwrap_or([f64::NAN; 2]);
                trace.extend_from_slice(&[frame.gap_m as f32, frame.flow_m3_s as f32, lower as f32, upper as f32]);
            }
            self.steps += 1;
        }
        let pressure = match self.nostril_heads.as_mut() {
            None => {
                let pressure = self.receiver.step_trace(&self.trace)?;
                // The head filters the baffle's pressure into the sphere's; without the row the pressure is what it was.
                match self.head.as_mut() {
                    Some(head) => head.step(pressure),
                    None => pressure,
                }
            }
            // The nasal rows with `head_sphere` (N3b): each opening's pressure through its own head FIR, then the sum.
            Some([left, right]) => {
                let r = self.ratio;
                let head = self.head.as_mut().ok_or("the nasal head has no mouth FIR")?;
                let mouth = head.step(self.receiver.step_trace(&self.trace[..r])?);
                let left = left.1.step(left.0.step_trace(&self.trace[r..2 * r])?);
                let right = right.1.step(right.0.step_trace(&self.trace[2 * r..3 * r])?);
                mouth + left + right
            }
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

    /// The `rest_gap` rows (MB 16735, the first word's piece 4) parse in time order, and each refusal names what it
    /// meets: the fields, a gap that is not positive, times that do not increase, and folds that have no two-mass rest
    /// gap to move (`one_dof`, `body_cover`, `glottal_orifice`).
    #[test]
    fn rest_gap_parses_and_is_refused_by_name() {
        let parse = |text: &str| VoiceDescription::from_bytes(text.as_bytes());
        let rows = "rest_gap 0 2.79e-3\nrest_gap 0.2 2.79e-3\nrest_gap 0.35 3.57e-4\n";
        let ok = parse(&modulation_text(rows)).expect("parses");
        assert_eq!(ok.rest_gap, vec![(0.0, 2.79e-3), (0.2, 2.79e-3), (0.35, 3.57e-4)]);
        assert!(parse(&modulation_text("")).expect("parses").rest_gap.is_empty());
        for (extra, needle) in [
            ("rest_gap 0\n", "takes 2 fields"),
            ("rest_gap 0 x\n", "not a number"),
            ("rest_gap 0 inf\n", "must be finite"),
            ("rest_gap -0.1 3e-4\n", "non-negative"),
            ("rest_gap 0 0\n", "gap positive"),
            ("rest_gap 0 -3e-4\n", "gap positive"),
            ("rest_gap 0.1 3e-4\nrest_gap 0.1 4e-4\n", "strictly increasing"),
            ("rest_gap 0.2 3e-4\nrest_gap 0.1 4e-4\n", "strictly increasing"),
        ] {
            let error = parse(&modulation_text(extra)).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        let one_dof = modulation_text(rows)
            .replace("fold two_mass", "fold one_dof")
            .replace("fold_damping_upper 0.6\n", "");
        let error = parse(&one_dof).expect_err("one_dof");
        assert!(error.contains("rest_gap needs a two_mass fold"), "{error}");
        let body_cover =
            modulation_text(&format!("fold_body 5.0e-5 100 0.2 3.598e-4\n{rows}")).replace("fold two_mass", "fold body_cover");
        let error = parse(&body_cover).expect_err("body_cover");
        assert!(error.contains("rest_gap needs a two_mass fold, not body_cover"), "{error}");
        let error = parse(&fixture(&format!("glottal_orifice 7.854e-5\n{rows}"))).expect_err("orifice");
        assert!(error.contains("`rest_gap` moves or reads a fold"), "{error}");
    }

    /// The rows reach the fold and nothing else. A schedule at the `fold` row's own gap is bit for bit the card without
    /// it, held from the start or from a later first row. A gesture that opens the fold and brings it back to that gap
    /// reaches the pressure and keeps the partition free. Held at Yoshioka's /s/ voicing-offset area (39 mm^2 on the
    /// 14 mm fold, a 2.79 mm uniform gap) the fold shuts less often than at the card's own gap.
    #[test]
    fn rest_gap_at_the_cards_own_gap_is_bit_identical_and_a_gesture_reaches_the_fold() {
        let plain = render(&fixture(""));
        for rows in ["rest_gap 0 3.57e-4\n", "rest_gap 0.01 3.57e-4\nrest_gap 0.04 3.57e-4\n"] {
            let own = render(&fixture(rows));
            assert_eq!(first_difference(&plain, &own), None, "{rows:?}: the card's own gap must change nothing");
        }
        let gesture = fixture("rest_gap 0 1.0e-3\nrest_gap 0.03 3.57e-4\n");
        let moved = render(&gesture);
        assert!(moved.iter().all(|v| v.is_finite()));
        assert!(first_difference(&plain, &moved).is_some(), "the gesture must reach the pressure");
        let description = VoiceDescription::from_bytes(gesture.as_bytes()).expect("parses");
        let (mut whole, mut split) = (description.build().expect("builds"), description.build().expect("builds"));
        let (mut a, mut b) = (vec![0.0_f64; 2400], vec![0.0_f64; 2400]);
        whole.block(&mut a).expect("whole");
        for chunk in b.chunks_mut(37) {
            split.block(chunk).expect("split");
        }
        assert_eq!(first_difference(&a, &b), None, "the gesture must not depend on the partition");
        let closed = |extra: &str| {
            let description = VoiceDescription::from_bytes(modulation_text(extra).as_bytes()).expect("parses");
            let mut performance = description.build().expect("builds");
            let mut out = vec![0.0; description.samples as usize];
            for block in out.chunks_mut(4096) {
                performance.block(block).expect("renders");
            }
            assert!(out.iter().all(|v| v.is_finite()));
            performance.stats().closed_fraction
        };
        let (card, abducted) = (closed(""), closed("rest_gap 0 2.79e-3\n"));
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-rest-gap-abducted\",\"closed_fraction_card\":{card},\"closed_fraction_39mm2\":{abducted}}}"
        );
        assert!(card > 0.0, "the card's own fold must shut ({card})");
        assert!(abducted < card, "an abducted fold shuts less often: {abducted} against {card}");
    }

    /// The `fold_tension` rows (MB 23138, the lyric phrase) parse in time order, and each refusal names what it meets:
    /// the fields, a stiffness scale outside the fold's tension range, a negative gesture scale, times that do not
    /// increase, a factor on a gesture the card lacks, a scaled vibrato or flutter past its row's own bound, and the
    /// flow rig (no fold).
    #[test]
    fn fold_tension_parses_and_is_refused_by_name() {
        let parse = |text: &str| VoiceDescription::from_bytes(text.as_bytes());
        let gestures = "tension_modulation 5.9 0.11\ntension_flutter 7 4 1.3 0.02\n";
        let rows = "fold_tension 0 1 1 1\nfold_tension 0.9 1.08 1.02 0.87\nfold_tension 1.8 1.58 3.8 2.9\n";
        let ok = parse(&modulation_text(&format!("{gestures}{rows}"))).expect("parses");
        let point = |time_s, stiffness, vibrato, flutter| FoldTensionPoint { time_s, stiffness, vibrato, flutter };
        assert_eq!(
            ok.fold_tension,
            vec![point(0.0, 1.0, 1.0, 1.0), point(0.9, 1.08, 1.02, 0.87), point(1.8, 1.58, 3.8, 2.9)]
        );
        assert!(parse(&modulation_text(gestures)).expect("parses").fold_tension.is_empty());
        parse(&modulation_text("fold_tension 0 1.2 1 1\n")).expect("a stiffness factor alone needs no gesture");
        for (extra, needle) in [
            ("fold_tension 0 1 1\n", "takes 4 fields"),
            ("fold_tension 0 x 1 1\n", "not a number"),
            ("fold_tension 0 inf 1 1\n", "must be finite"),
            ("fold_tension -0.1 1 1 1\n", "time must be non-negative"),
            ("fold_tension 0 0 1 1\n", "stiffness scale must be in 0.25..=4"),
            ("fold_tension 0 -1 1 1\n", "stiffness scale must be in 0.25..=4"),
            ("fold_tension 0 0.24 1 1\n", "stiffness scale must be in 0.25..=4"),
            ("fold_tension 0 4.01 1 1\n", "stiffness scale must be in 0.25..=4"),
            ("fold_tension 0 1 -1 1\n", "scales must be non-negative"),
            ("fold_tension 0 1 1 -1\n", "scales must be non-negative"),
            ("fold_tension 0.1 1 1 1\nfold_tension 0.1 1.1 1 1\n", "strictly increasing"),
            ("fold_tension 0.2 1 1 1\nfold_tension 0.1 1.1 1 1\n", "strictly increasing"),
            ("fold_tension 0 1 4.6 1\n", "above 0.5"),
            ("fold_tension 0 1 1 5.1\n", "sd above 0.1"),
        ] {
            let error = parse(&modulation_text(&format!("{gestures}{extra}"))).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        for (extra, needle) in [
            ("fold_tension 0 1 1.1 1\n", "vibrato scale needs tension_modulation"),
            ("fold_tension 0 1 1 1.1\n", "flutter scale needs tension_flutter"),
        ] {
            let error = parse(&modulation_text(extra)).expect_err(extra);
            assert!(error.contains(needle), "{extra:?} -> {error}");
        }
        let error = parse(&fixture("glottal_orifice 7.854e-5\nfold_tension 0 1 1 1\n")).expect_err("orifice");
        assert!(error.contains("`fold_tension` moves or reads a fold"), "{error}");
    }

    /// Rows at `1 1 1` change nothing, held from the start or from a later first row, on a card with both gestures and
    /// on a card with neither: a product with 1 is exact, and a tension of 1 is the island's own (a no-op).
    #[test]
    fn fold_tension_at_one_is_bit_identical() {
        let gestures =
            "tension_modulation 5.9 0.11\ntension_flutter 7 20 20 0.005\ntension_flutter_second 7 20 20 0.002\n";
        for card in ["", gestures] {
            let plain = render(&fixture(card));
            for rows in ["fold_tension 0 1 1 1\n", "fold_tension 0.01 1 1 1\nfold_tension 0.04 1 1 1\n"] {
                let own = render(&fixture(&format!("{card}{rows}")));
                let at = first_difference(&plain, &own);
                assert_eq!(at, None, "{card:?} {rows:?}: rows at 1 1 1 must change nothing");
            }
        }
    }

    /// The rows reach the fold exactly as the module doc writes them: a card with the rows plays the card without them,
    /// its tension set by hand before each output sample to `s (1 + (g_v a sin(2 pi r t) + g_f (y1 + y2)))`, with the
    /// factors interpolated as written out here and the flutter streams written out from their keys. The factors reach
    /// the pressure, and the result does not depend on the partition.
    #[test]
    fn fold_tension_reaches_the_fold_as_written() {
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let (rate_hz, amplitude) = (5.9, 0.11);
        let gestures = format!(
            "tension_modulation {rate_hz} {amplitude}\ntension_flutter 7 20 20 0.005\n\
             tension_flutter_second 7 20 20 0.002\n"
        );
        // Two rows: the first held before 0.01 s, the second after 0.03 s, linear between.
        let rows = [(0.01, 1.0, 1.0, 1.0), (0.03, 1.2, 1.5, 0.5)];
        let card = fixture(&format!(
            "{gestures}fold_tension {} {} {} {}\nfold_tension {} {} {} {}\n",
            rows[0].0, rows[0].1, rows[0].2, rows[0].3, rows[1].0, rows[1].1, rows[1].2, rows[1].3
        ));
        let engine = render(&card);
        let plain = fixture("");
        let n = VoiceDescription::from_bytes(plain.as_bytes()).expect("parses").samples as usize;
        let first = TensionFlutterCard { seed: 7, centre_hz: 20.0, bandwidth_hz: 20.0, sd_rel: 0.005 };
        let second = TensionFlutterCard { seed: 7, centre_hz: 20.0, bandwidth_hz: 20.0, sd_rel: 0.002 };
        let (y1, y2) = (flutter_written_out(first, 0, 48_000, n), flutter_written_out(second, 1, 48_000, n));
        let factor = |t: f64, (a, b): (f64, f64)| -> f64 {
            let (t0, t1) = (rows[0].0, rows[1].0);
            if t <= t0 {
                a
            } else if t >= t1 {
                b
            } else {
                a + (b - a) * (t - t0) / (t1 - t0)
            }
        };
        let description = VoiceDescription::from_bytes(plain.as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        let mut by_hand = vec![0.0; n];
        for (k, sample) in by_hand.iter_mut().enumerate() {
            let t = performance.steps as f64 / f64::from(performance.mechanics_rate);
            let s = factor(t, (rows[0].1, rows[1].1));
            let v = (factor(t, (rows[0].2, rows[1].2)) * amplitude) * det::sin(core::f64::consts::TAU * rate_hz * t);
            let y = factor(t, (rows[0].3, rows[1].3)) * (y1[k] + y2[k]);
            performance.island.set_tension(s * (1.0 + (v + y))).expect("tension");
            performance.block(core::slice::from_mut(sample)).expect("renders");
        }
        assert_eq!(bits(&engine), bits(&by_hand), "the rows must reach the fold as the module doc writes them");
        let unscaled = render(&fixture(&gestures));
        assert!(first_difference(&unscaled, &engine).is_some(), "the factors must reach the pressure");
        let description = VoiceDescription::from_bytes(card.as_bytes()).expect("parses");
        let mut split = description.build().expect("builds");
        let mut chunked = vec![0.0_f64; n];
        for chunk in chunked.chunks_mut(37) {
            split.block(chunk).expect("split");
        }
        assert_eq!(first_difference(&engine, &chunked), None, "the rows must not depend on the partition");
    }

    /// A held stiffness factor is the card whose `fold` row carries the scaled stiffnesses: every stiffness-derived
    /// quantity goes through the fold's one tensioning path, and 80, 8 and 25 times 1.21 are exactly 96.8, 9.68 and
    /// 30.25 in binary64, so the two render bit for bit the same. A step to the factor then settles on that note: the
    /// read window starts 0.25 s after the step ends, and 1e-3 (1.7 cents) is a stated bound on what is left of the
    /// self-oscillation's return to its cycle plus the f0 track's window-phase sensitivity, not a derived one (the
    /// value is printed). The factor must move f0 by more than 20 times that bound.
    #[test]
    fn a_stiffness_factor_plays_the_scaled_cards_note() {
        let rate = 48_000.0;
        let mean_f0 = |pcm: &[f64]| {
            let tail = &pcm[(0.9 * rate) as usize..];
            let track = f0_track(tail, rate, steady_f0_hz(tail, rate));
            assert!(!track.is_empty(), "the read window holds no f0 frame");
            track.iter().map(|(_, f)| f).sum::<f64>() / track.len() as f64
        };
        let card = modulation_text("");
        let scaled_text = card.replace("2.5e-5 80 8 25 0.1", "2.5e-5 96.8 9.68 30.25 0.1");
        assert_ne!(scaled_text, card, "the scaled card must differ from the card");
        let at_one = mean_f0(&render(&card));
        let (row_pcm, scaled_pcm) = (render(&modulation_text("fold_tension 0 1.21 1 1\n")), render(&scaled_text));
        assert_eq!(first_difference(&row_pcm, &scaled_pcm), None, "the held factor must be the scaled card bit for bit");
        let row = mean_f0(&row_pcm);
        let stepped = mean_f0(&render(&modulation_text("fold_tension 0.6 1 1 1\nfold_tension 0.65 1.21 1 1\n")));
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-fold-tension-note\",\"f0_card_hz\":{at_one},\"f0_row_hz\":{row},\"f0_stepped_hz\":{stepped}}}"
        );
        assert!((stepped / row - 1.0).abs() < 1e-3, "stepped {stepped} Hz against the held row {row} Hz");
        assert!(row / at_one - 1.0 > 2e-2, "the factor must raise f0: {row} Hz against {at_one} Hz");
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
    /// moves by the ratio of the two sound speeds to within 1 % (the lips' load is an absolute impedance whose radiation
    /// stays in the exterior gas, only its flanged term's eliminated duct modes taking the tract's, so the end correction
    /// seen from the tract does not follow the sound speed ratio), and the receiver and the head's cap, piston and torso
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
        let unflanged =
            fixture("tract_gas moist_air 310.15 101325 1.0\n").replace("termination flanged", "termination unflanged");
        let err = VoiceDescription::from_bytes(unflanged.as_bytes()).unwrap_err();
        assert!(err.contains("needs `termination flanged`"), "{err}");
        let rows: [(String, &str); 3] = [
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
        // The fields are public, so `build` refuses a tract gas with neither the head's load nor a flanged mouth.
        let mut flipped = VoiceDescription::from_bytes(
            fixture(&format!("{HEAD_LOAD}tract_gas moist_air 290.15 101325 0.5\n")).as_bytes(),
        )
        .unwrap();
        flipped.head_load = false;
        flipped.termination = Termination::UnflangedOpen;
        let err = flipped.build().err().expect("refused");
        assert!(err.contains("needs `termination flanged`"), "{err}");
    }

    /// `tract_gas` with no head (MB 23833, 24424): the lips' load is the declared flanged load across the gas interface
    /// (radiation in the exterior, the eliminated duct modes in the tract's gas). The formants are, bit for bit, the
    /// tract gas's TMM against that table made by hand, and not the TMM with the
    /// flanged load evaluated in the tract's gas; they move by the sound-speed ratio to within 1 % (as with the head's
    /// load), the receiver stays the exterior's, and a tract gas equal to the ambient gives the card's own formants to
    /// within the table's linear interpolation between rows (1e-6).
    #[test]
    fn tract_gas_without_the_head_load_reads_the_flanged_load_in_the_exterior() {
        let base = VoiceDescription::from_bytes(fixture("").as_bytes()).unwrap();
        let breath =
            VoiceDescription::from_bytes(fixture("tract_gas moist_air 290.15 101325 0.5\n").as_bytes()).unwrap();
        let same =
            VoiceDescription::from_bytes(fixture("tract_gas moist_air 310.15 101325 1.0\n").as_bytes()).unwrap();
        assert_eq!(same.tract_gas(), same.ambient_gas());
        let (a, b, s) = (base.build().unwrap(), breath.build().unwrap(), same.build().unwrap());
        assert!(!a.info().head_load && !b.info().head_load && !s.info().head_load);
        let chart = TractChart::try_new(breath.sections.clone(), &breath.tract_source, &breath.tract_licence).unwrap();
        let rate = breath.output_rate_hz * breath.mechanics_per_output;
        let (tract, exterior) = (breath.tract_gas(), breath.ambient_gas());
        let table = flanged_mouth_load(&chart.to_duct(), tract, exterior, rate, GLOTTAL_LINE_MAX_TAPS).unwrap();
        let by_hand = chart.tmm_formants_tabulated(breath.tract_gas(), breath.wall.as_ref(), &table, 4).unwrap();
        let fb = b.info().passive_formants_hz.clone().unwrap();
        assert_eq!(fb, by_hand);
        let in_tract = chart
            .tmm_formants_terminated(breath.tract_gas(), breath.wall.as_ref(), Termination::FlangedOpen, 4)
            .unwrap();
        assert_ne!(fb, in_tract, "the lips' load must not be evaluated in the tract's gas");
        let ratio = breath.tract_gas().sound_speed / base.ambient_gas().sound_speed;
        assert!(ratio < 0.97, "the test gas must move the sound speed: {ratio}");
        let fa = a.info().passive_formants_hz.clone().unwrap();
        assert_eq!(fa.len(), fb.len());
        for (x, y) in fa.iter().zip(&fb) {
            assert!((y / x / ratio - 1.0).abs() < 0.01, "formant {x} -> {y} Hz against the ratio {ratio}");
        }
        assert_eq!(b.filter_delay_samples(), a.filter_delay_samples(), "the exterior's receiver is unchanged");
        let fs = s.info().passive_formants_hz.clone().unwrap();
        assert_eq!(fa.len(), fs.len());
        for (x, y) in fa.iter().zip(&fs) {
            assert!((y / x - 1.0).abs() < 1e-6, "equal gases: formant {x} -> {y} Hz");
        }
        assert_ne!(pressure(&base), pressure(&breath));
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

    /// The render option `--glottal-trace` (MB 22209): keeping the trace leaves the pressure bit for bit, on the
    /// uniform fold and on a shaped one. The trace holds one row per physical substep; its last row is the island's own
    /// gap and displacements at the end of the render; its flow column's extremes (from 0) are the stats' peak and
    /// minimum; on the uniform fold its rows with no open gap are the steps `closed_fraction` counts, and on the shaped
    /// fold (contact on any strip) they are no more than those.
    #[test]
    fn glottal_trace_leaves_the_pressure_bit_for_bit_and_reads_the_island() {
        const N: usize = GLOTTAL_TRACE_COLUMNS.len();
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        for (extra, uniform) in [("", true), ("fold_shape sine uniform\n", false)] {
            let text = modulation_text(extra);
            let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
            let plain = render(&text);
            let mut kept = description.build().expect("builds");
            assert!(kept.glottal_trace().is_none(), "nothing kept before the call");
            let substeps = description.samples as usize * description.mechanics_per_output as usize;
            kept.keep_glottal_trace(substeps);
            let mut out = vec![0.0; description.samples as usize];
            for block in out.chunks_mut(4096) {
                kept.block(block).expect("renders");
            }
            assert_eq!(bits(&plain), bits(&out), "the trace changed the pressure ({extra:?})");
            assert_eq!(kept.info().mechanics_rate_hz, 48_000 * description.mechanics_per_output);
            let (rows, rest) = kept.glottal_trace().expect("kept").as_chunks::<N>();
            assert!(rest.is_empty());
            assert_eq!(rows.len(), substeps, "one row per substep");
            assert_eq!(kept.steps as usize, substeps);
            let last = rows[rows.len() - 1];
            assert_eq!(last[0], kept.island.gap_m() as f32, "gap column");
            let [lower, upper] = kept.island.mass_displacements_m().expect("two masses");
            assert_eq!((last[2], last[3]), (lower as f32, upper as f32), "displacement columns");
            let stats = kept.stats();
            // The stats' peak and minimum start from 0, so the folds over the column do too.
            let flows = rows.iter().map(|row| row[1]);
            assert_eq!(flows.clone().fold(0.0_f32, f32::max), stats.peak_glottal_flow_m3_s as f32, "peak flow");
            assert_eq!(flows.fold(0.0_f32, f32::min), stats.min_glottal_flow_m3_s as f32, "minimum flow");
            let shut_rows = rows.iter().filter(|row| row[0] <= 0.0).count() as u64;
            assert!(kept.closed_steps > 0, "the fixture's fold makes contact ({extra:?})");
            if uniform {
                assert_eq!(shut_rows, kept.closed_steps, "uniform: no open gap is contact");
            } else {
                // The sine mode holds the ends, so the whole length may never shut while some strip touches.
                let contact = kept.closed_steps;
                assert!(shut_rows <= contact, "shaped: {shut_rows} rows shut, {contact} steps in contact");
            }
        }
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

    /// The ratio of the `i2` keyframe's first section to the card's in
    /// `a_moving_tract_plays_the_keyframes_own_note_on_its_own_glottal_section`.
    const I2_FIRST_RATIO: f64 = 1.5;

    /// A moving tract's glottal section moves with it (MB 25412 step 1, MB 25468). The /a/-like card moves over 0.10 to
    /// 0.18 s to `i2`, the /i/-like keyframe with its first section `I2_FIRST_RATIO` times the card's, and from 0.30
    /// to 0.60 s it plays the static `i2` card's note: the mean f0 (`f0_track` on the prior of `steady_f0_hz`) within 3
    /// cents. The /i/-like shape with the card's first section, which is the tract a held glottal section left after
    /// the move, must play a note more than 6 cents from `i2`'s, twice the tolerance, so that a tract that held the
    /// card's section could not pass. The move's largest curvature against the larger of the held windows' (before
    /// the move from 0.06 s, and the read window) is printed, not asserted: the island's `zc` and `A_1` step at each
    /// switch, a stated limit, and this is its measure on the fixture.
    #[test]
    fn a_moving_tract_plays_the_keyframes_own_note_on_its_own_glottal_section() {
        let (seconds, rate) = (0.6, 48_000.0);
        let i2: Vec<f64> =
            SHAPE_I.iter().enumerate().map(|(k, &area)| if k == 0 { I2_FIRST_RATIO * area } else { area }).collect();
        let mut rows = String::from("keyframe i2 source test-i2\nkeyframe i2 licence test fixture, not data\n");
        for area in &i2 {
            rows.push_str(&format!("keyframe_section i2 {area:e} 1.7e-2\n"));
        }
        rows.push_str("tract_move 0.10 0.18 i2 8\n");
        let render = |text: String| {
            let description = VoiceDescription::from_bytes(text.as_bytes()).expect("parses");
            render_performance(description.build().expect("builds"), (seconds * rate) as usize)
        };
        let mean_f0 = |pcm: &[f64]| {
            let window = &pcm[(0.30 * rate) as usize..(0.60 * rate) as usize];
            let track = f0_track(window, rate, steady_f0_hz(window, rate));
            assert!(!track.is_empty(), "the read window holds no f0 frame");
            track.iter().map(|(_, f)| f).sum::<f64>() / track.len() as f64
        };
        let cents = |a: f64, b: f64| 1200.0 * (a / b).log2();
        let moving_p = render(moving_card(&SHAPE_A, 0.02, seconds, &rows));
        let held_curvature = largest_curvature(&moving_p, 2_880..4_800).max(largest_curvature(&moving_p, 14_400..28_799));
        let moved_over_held = largest_curvature(&moving_p, 4_800..8_640) / held_curvature;
        let moved = mean_f0(&moving_p);
        let own = mean_f0(&render(moving_card(&i2, 0.017, seconds, "")));
        let held = mean_f0(&render(moving_card(&SHAPE_I, 0.017, seconds, "")));
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-moving-tract-glottal-section\",\"first_ratio\":{I2_FIRST_RATIO},\
             \"f0_moved_hz\":{moved:.4},\"f0_own_hz\":{own:.4},\"f0_held_hz\":{held:.4},\"moved_c\":{:.3},\"held_c\":{:.3},\
             \"moved_over_held_curvature\":{moved_over_held:.3}}}",
            cents(moved, own),
            cents(held, own)
        );
        assert!(
            cents(held, own).abs() > 6.0,
            "the card's first section moves the keyframe's note by {:.2} c, not more than 6: the test cannot fail",
            cents(held, own)
        );
        assert!(
            cents(moved, own).abs() <= 3.0,
            "after the move {moved:.3} Hz against the keyframe's own {own:.3} Hz ({:.2} c)",
            cents(moved, own)
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

    /// `tract_gas` with no head on a moving tract (MB 24424): `shape_taps` makes each further shape's own flanged table
    /// across the gas interface, from that shape's duct. The keyframe target's passive formants are, bit for bit, the
    /// tract gas's TMM against that table made by hand on the shape `build` makes (the log-area morph at fraction 1, the
    /// keyframe as written), and not the TMM with the flanged load evaluated in the tract's gas.
    #[test]
    fn a_moving_tract_with_tract_gas_and_no_head_reads_each_shapes_own_flanged_table() {
        let rows = keyframe_i() + "tract_move 0.05 0.1 i 3\ntract_gas moist_air 290.15 101325 0.5\n";
        let description =
            VoiceDescription::from_bytes(moving_card(&SHAPE_A, 0.02, 0.15, &rows).as_bytes()).expect("parses");
        assert!(!description.head_load && description.tract_gas_record().is_some());
        let built = description.build().expect("builds");
        let formants = &built.info().keyframe_formants_hz;
        assert_eq!(formants.len(), 1);
        assert_eq!(formants[0].0, "i");
        let card = TractChart::try_new(description.sections.clone(), &description.tract_source,
            &description.tract_licence).expect("card chart");
        let k = &description.keyframes[0];
        let keyframe = TractChart::try_new(k.sections.clone(), &k.source, &k.licence).expect("keyframe chart");
        // The move's last step: `k / steps` with k = steps = 3 is exactly 1.
        let shape = card.log_area_morph(&keyframe, 1.0).expect("morph");
        let (tract, exterior, wall) = (description.tract_gas(), description.ambient_gas(), description.wall.as_ref());
        let rate = description.output_rate_hz * description.mechanics_per_output;
        let table = flanged_mouth_load(&shape.to_duct(), tract, exterior, rate, GLOTTAL_LINE_MAX_TAPS).expect("table");
        let by_hand = shape.tmm_formants_tabulated(tract, wall, &table, 4).expect("formants");
        assert_eq!(formants[0].1.as_ref().expect("keyframe formants"), &by_hand);
        let in_tract = shape.tmm_formants_terminated(tract, wall, Termination::FlangedOpen, 4).expect("formants");
        assert_ne!(by_hand, in_tract, "the keyframe's lips' load must not be evaluated in the tract's gas");
    }

    /// Build time per shape at the live cards' clock and rows (head sphere and load, torso, tract gas, 192 kHz), which a
    /// render's stop line reads before any moving card renders. The four further shapes build at once on the host's
    /// threads, so their wall is one round of up to `MAX_SHAPE_WORKERS` shapes, and a move of `n` shapes builds in about
    /// `fixed_build_s + ceil(n / shape_workers) * further_shapes_wall_s`. Test shapes of 90 sections, not data; nothing
    /// is rendered. The keyframe's
    /// mouth stays at 3.0e-4 m^2 (9.8 mm), inside the live cards' 8.6-14.0 mm: a 1.0e-4 m^2 mouth (5.6 mm) on the
    /// 8.7 cm head needs about 76,000-78,000 of the head load's orders at 192 kHz, which the old cap of 65536 refused by
    /// name and `LOAD_MAX_ORDER` (2^17 since N3b) admits.
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

    /// The flow rig's rows (MB 15368; 15773, 15828, 15902; NavyTower 15765) parse, and each refusal names what it meets:
    /// the fold rows with the orifice, the two airways together, the fricative source with an unusable `k`/`j`, on a fixed
    /// tract with no expansion, on a moving tract with a card `j` that is not the card shape's largest section after `k`,
    /// the monopole with the baffle's own observers, and the thin-wall load outside its quadrature range. The fricative
    /// source on the folds (the word, MB 16735 piece 1) and on a moving tract (piece 2) parses, with a keyframe of another
    /// section count and with a card shape that is no expansion (NavyTower 16887), and a fixed tract keeps any wider `j`
    /// the card names.
    #[test]
    fn flow_rig_rows_parse_and_are_refused_by_name() {
        let rig = "glottal_orifice 7.854e-5\nsubglottal_tube 2.5e-4 0.08\nsubglottal_tube 1.2e-3 0.038\nfricative_noise 5 1 4\nreceiver free_field_monopole\n";
        let ok = VoiceDescription::from_bytes(fixture(rig).as_bytes()).expect("the rig parses");
        assert_eq!(ok.glottal_orifice, Some(7.854e-5));
        assert_eq!(ok.subglottal_tube.len(), 2);
        assert_eq!(
            ok.fricative_noise,
            Some(FricativeNoiseCard {
                seed: 5,
                constriction: 1,
                expansion: 4,
                law: FricativeLaw::OldhamUkpoho,
                obstacle: 1,
            })
        );
        assert!(ok.free_field_receiver);
        let plain = VoiceDescription::from_bytes(fixture("").as_bytes()).expect("parses");
        assert!(plain.glottal_orifice.is_none() && plain.subglottal_tube.is_empty() && plain.fricative_noise.is_none());
        assert!(!plain.free_field_receiver);
        let thin = fixture("").replace("termination flanged", "termination unflanged_thin_wall 256");
        assert_eq!(
            VoiceDescription::from_bytes(thin.as_bytes()).expect("parses").termination,
            Termination::UnflangedThinWall { quadrature: LsQuadrature { nodes: 256 } }
        );
        let orifice = "glottal_orifice 7.854e-5\n";
        let cases: Vec<(String, &str)> = vec![
            (format!("{orifice}fold_shape sine uniform\n"), "`fold_shape` moves or reads a fold"),
            (format!("{orifice}tension_modulation 5 0.02\n"), "`tension_modulation`"),
            (format!("{orifice}pulse_shimmer 3 0.05\n"), "`pulse_shimmer`"),
            (format!("{orifice}aspiration_noise 3 1e7\n"), "`aspiration_noise`"),
            (format!("{orifice}constriction_noise 7 5.3e-3\n"), "`constriction_noise`"),
            ("glottal_orifice 0\n".into(), "must be positive"),
            (format!("{orifice}glottal_orifice 1e-4\n"), "only appear once"),
            ("subglottal_tube 2.5e-4 0.08\nsubglottal_airway weibel-1963 saab-2021-table-1-B 3\n".into(), "refused together"),
            ("subglottal_tube 2.5e-4 0\n".into(), "must be positive"),
            (format!("{orifice}fricative_noise 5 4 1\n"), "k < j"),
            (format!("{orifice}fricative_noise 5 1 8\n"), "k < j"),
            (format!("{orifice}fricative_noise 5 4 5\n"), "wider than section k"),
            (format!("{orifice}fricative_noise 5 1\n"), "takes 3 fields"),
            ("receiver free_field_monopole\nhead_sphere 0.0875 src\n".into(), "refused with head_sphere and observer_law"),
            ("receiver baffle\n".into(), "free_field_monopole"),
        ];
        for (rows, word) in &cases {
            let error = VoiceDescription::from_bytes(fixture(rows).as_bytes()).unwrap_err();
            assert!(error.contains(word), "{rows:?}: {error}");
        }
        let moving = fixture(&format!(
            "{orifice}fricative_noise 5 1 4\nkeyframe b source src\nkeyframe b licence test\n{}tract_move 0.01 0.02 b 2\n",
            [0.5e-4, 0.7e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4]
                .iter()
                .map(|a| format!("keyframe_section b {a:e} 0.02\n"))
                .collect::<String>()
        ));
        VoiceDescription::from_bytes(moving.as_bytes()).expect("the source on a moving tract parses (MB 16735 piece 2)");
        let not_widest = moving.replace("fricative_noise 5 1 4", "fricative_noise 5 1 3");
        assert_ne!(not_widest, moving, "the replacement must change the card");
        let error = VoiceDescription::from_bytes(not_widest.as_bytes()).unwrap_err();
        assert!(error.contains("largest section between k and the lips"), "{error}");
        VoiceDescription::from_bytes(fixture(&format!("{orifice}fricative_noise 5 1 3\n")).as_bytes())
            .expect("a fixed tract keeps the card's j");
        let recounted = format!("{moving}keyframe_section b 2.0e-4 0.02\n");
        VoiceDescription::from_bytes(recounted.as_bytes()).expect("a keyframe of another section count parses (16887 (c))");
        let vowel = moving.replace("fricative_noise 5 1 4", "fricative_noise 5 5 6");
        assert_ne!(vowel, moving, "the replacement must change the card");
        VoiceDescription::from_bytes(vowel.as_bytes()).expect("a moving card's own shape need not be an expansion");
        let error = VoiceDescription::from_bytes(fixture(&format!("{orifice}fricative_noise 5 5 6\n")).as_bytes()).unwrap_err();
        assert!(error.contains("wider than section k"), "a fixed card's one shape must be an expansion: {error}");
        let folds = VoiceDescription::from_bytes(fixture("fricative_noise 5 1 4\n").as_bytes()).expect("folds parse");
        assert!(folds.glottal_orifice.is_none());
        assert_eq!(
            folds.fricative_noise,
            Some(FricativeNoiseCard {
                seed: 5,
                constriction: 1,
                expansion: 4,
                law: FricativeLaw::OldhamUkpoho,
                obstacle: 1,
            })
        );
        for nodes in ["8", "4096", "x"] {
            let text = fixture("").replace("termination flanged", &format!("termination unflanged_thin_wall {nodes}"));
            let error = VoiceDescription::from_bytes(text.as_bytes()).unwrap_err();
            assert!(error.contains("unflanged_thin_wall quadrature nodes"), "{nodes}: {error}");
        }
    }

    /// The rig renders: finite, the fricative row reaches the radiated pressure, one seed gives one noise and another seed
    /// another, and the orifice card's fold numbers reach nothing (the same rig on a stiffer fold renders bit for bit). A
    /// hole wider than the tract's first section (the fixture's 0.5 cm^2) has no eq 6 expansion and is refused at build.
    #[test]
    fn flow_rig_renders_its_fricative_seeded_and_the_fold_reaches_nothing() {
        let wide = fixture("glottal_orifice 7.854e-5\nreceiver free_field_monopole\n");
        let error = VoiceDescription::from_bytes(wide.as_bytes()).expect("parses").build().err().expect("refused");
        assert!(error.contains("no expansion"), "{error}");
        let base = "glottal_orifice 4.0e-5\nsubglottal_tube 2.5e-4 0.08\nreceiver free_field_monopole\n";
        let quiet = render(&fixture(base));
        let spoken = render(&fixture(&format!("{base}fricative_noise 5 1 4\n")));
        assert!(quiet.iter().chain(&spoken).all(|x| x.is_finite()));
        assert!(quiet.iter().zip(&spoken).any(|(a, b)| a.to_bits() != b.to_bits()), "the fricative must reach the pressure");
        assert_eq!(spoken, render(&fixture(&format!("{base}fricative_noise 5 1 4\n"))), "one seed, one noise");
        assert_ne!(spoken, render(&fixture(&format!("{base}fricative_noise 6 1 4\n"))), "another seed, another noise");
        let stiff = fixture(&format!("{base}fricative_noise 5 1 4\n")).replace(
            "fold two_mass 1.25e-4 2.5e-5 128 12.8 40",
            "fold two_mass 1.25e-4 2.5e-5 256 25.6 80",
        );
        assert_ne!(stiff, fixture(&format!("{base}fricative_noise 5 1 4\n")), "the replacement must change the card");
        let bits = |xs: &[f64]| xs.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
        assert_eq!(bits(&spoken), bits(&render(&stiff)), "through an orifice the fold's numbers reach nothing");
        // The rig's read window: this card is shorter than 2 s, so the mean runs over the whole render.
        let rig = VoiceDescription::from_bytes(fixture(base).as_bytes()).expect("parses");
        let mut performance = rig.build().expect("builds");
        let mut out = vec![0.0; rig.samples as usize];
        performance.block(&mut out).expect("renders");
        let stats = performance.stats();
        let mean = stats.rig_mean_glottal_flow_m3_s.expect("the rig's mean flow");
        assert!(mean > 0.0 && mean <= stats.peak_glottal_flow_m3_s, "rig mean flow {mean} against the peak {}", stats.peak_glottal_flow_m3_s);
        let plain = VoiceDescription::from_bytes(fixture("").as_bytes()).expect("parses").build().expect("builds");
        assert!(plain.stats().rig_mean_glottal_flow_m3_s.is_none(), "a fold card has no rig window");
    }

    /// The Stevens form (MB 17504, 17645, 17712), test 4 of 4: the row's fourth field names the printed law and nothing
    /// else is read as one; a row without it is Oldham and Ukpoho's, so a card written before the field renders as it did
    /// (the rendered byte controls of the change carry that across the whole engine). The fifth field is the obstacle
    /// section `m`: without it `m = k`, and `m = k` written out renders bit for bit the same; `m < k`, `m` on the last
    /// section and an unreadable `m` are refused by name, and an `m` other than `k` parses on a moving tract too (MB
    /// 25469; `fricative_obstacle_on_a_moving_tract_is_each_shapes_own_m`). On the rig the named law
    /// reaches the radiated pressure, renders finite, one seed gives one noise, it is not the other law's noise, and the
    /// source at an obstacle downstream is not the source at `k`.
    #[test]
    fn fricative_noise_names_its_printed_law_and_keeps_oldham_ukpoho_without_it() {
        let parse = |row: &str| VoiceDescription::from_bytes(fixture(row).as_bytes());
        let card = |row: &str| parse(row).expect("parses").fricative_noise.expect("the row");
        assert_eq!(card("fricative_noise 5 1 4\n").law, FricativeLaw::OldhamUkpoho);
        assert_eq!(card("fricative_noise 5 1 4\n").obstacle, 1);
        assert_eq!(card("fricative_noise 5 1 4 stevens\n").law, FricativeLaw::StevensShadle);
        assert_eq!(card("fricative_noise 5 1 4 stevens\n").obstacle, 1, "no fifth field: m = k");
        assert_eq!(card("fricative_noise 5 1 4 stevens 3\n").obstacle, 3);
        assert_eq!(card("fricative_noise 5 1 4 stevens 6\n").obstacle, 6, "the section before the last");
        for (row, word) in [
            ("fricative_noise 5 1 4 shadle\n", "law must be `stevens`"),
            ("fricative_noise 5 1 4 Stevens\n", "law must be `stevens`"),
            ("fricative_noise 5 1 4 shadle 3\n", "law must be `stevens`"),
            ("fricative_noise 5 1 4 stevens 3 1\n", "takes 3 fields, or 4"),
            ("fricative_noise 5 1\n", "takes 3 fields"),
            ("fricative_noise 5 1 4 stevens x\n", "fricative_noise m"),
            ("fricative_noise 5 1 4 stevens 0\n", "k <= m < the last section 7"),
            ("fricative_noise 5 1 4 stevens 7\n", "k <= m < the last section 7"),
        ] {
            let error = parse(row).unwrap_err();
            assert!(error.contains(word), "{row:?}: {error}");
        }
        let keyframes = [0.5e-4, 0.7e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4]
            .iter()
            .map(|a| format!("keyframe_section b {a:e} 0.02\n"))
            .collect::<String>();
        let moving = |tail: &str| {
            fixture(&format!(
                "glottal_orifice 7.854e-5\nfricative_noise 5 1 4{tail}\nkeyframe b source src\nkeyframe b licence test\n\
                 {keyframes}tract_move 0.01 0.02 b 2\n"
            ))
        };
        VoiceDescription::from_bytes(moving(" stevens").as_bytes()).expect("the law on a moving tract parses");
        VoiceDescription::from_bytes(moving(" stevens 1").as_bytes()).expect("m = k on a moving tract parses");
        VoiceDescription::from_bytes(moving(" stevens 3").as_bytes()).expect("m after k on a moving tract parses");
        let base = "glottal_orifice 4.0e-5\nsubglottal_tube 2.5e-4 0.08\nreceiver free_field_monopole\n";
        let quiet = render(&fixture(base));
        let oldham = render(&fixture(&format!("{base}fricative_noise 5 1 4\n")));
        let stevens = render(&fixture(&format!("{base}fricative_noise 5 1 4 stevens\n")));
        assert!(stevens.iter().all(|x| x.is_finite()));
        assert!(quiet.iter().zip(&stevens).any(|(a, b)| a.to_bits() != b.to_bits()), "the named law must reach the pressure");
        assert_ne!(stevens, oldham, "the named law is not the other law's noise");
        assert_eq!(stevens, render(&fixture(&format!("{base}fricative_noise 5 1 4 stevens\n"))), "one seed, one noise");
        let bits = |xs: &[f64]| xs.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
        let at_k = render(&fixture(&format!("{base}fricative_noise 5 1 4 stevens 1\n")));
        assert_eq!(bits(&stevens), bits(&at_k), "m = k written out is the row without m, bit for bit");
        let downstream = render(&fixture(&format!("{base}fricative_noise 5 1 4 stevens 3\n")));
        assert!(downstream.iter().all(|x| x.is_finite()));
        assert_ne!(bits(&downstream), bits(&stevens), "the source at an obstacle downstream is not the source at k");
    }

    /// `fricative_throat` (MB 25440 ruling A, 25482) parses with its source, and each refusal names what it meets: no
    /// `fricative_noise`, `tract_move`, an area that is not a finite positive number, a throat not narrower than section
    /// `j` or more than 1.5 times section `k`, a field missing or extra, and the row twice.
    #[test]
    fn fricative_throat_parses_and_is_refused_by_name() {
        let rig = "glottal_orifice 4.0e-5\nfricative_noise 5 1 4 stevens 3\n";
        let ok = VoiceDescription::from_bytes(fixture(&format!("{rig}fricative_throat 6.5e-5 route-read\n")).as_bytes())
            .expect("the row parses");
        assert_eq!(ok.fricative_throat, Some(FricativeThroat { area_m2: 6.5e-5, source: "route-read".into() }));
        assert!(VoiceDescription::from_bytes(fixture(rig).as_bytes()).expect("parses").fricative_throat.is_none());
        VoiceDescription::from_bytes(fixture(&format!("{rig}fricative_throat 8.9e-5 src\n")).as_bytes())
            .expect("a throat inside 1.5 times section k (6e-5) parses");
        let cases: Vec<(String, &str)> = vec![
            ("fricative_throat 6.5e-5 src\n".into(), "needs fricative_noise"),
            (format!("{rig}fricative_throat 0 src\n"), "finite and positive"),
            (format!("{rig}fricative_throat -6.5e-5 src\n"), "finite and positive"),
            (format!("{rig}fricative_throat nan src\n"), "must be finite"),
            (format!("{rig}fricative_throat inf src\n"), "must be finite"),
            (format!("{rig}fricative_throat x src\n"), "not a number"),
            (format!("{rig}fricative_throat 3e-4 src\n"), "narrower than section j"),
            (format!("{rig}fricative_throat 9.1e-5 src\n"), "more than 1.5 times section k"),
            (format!("{rig}fricative_throat 6.5e-5\n"), "takes 2 fields"),
            (format!("{rig}fricative_throat 6.5e-5 src more\n"), "takes 2 fields"),
            (format!("{rig}fricative_throat 6.5e-5 src\nfricative_throat 6.5e-5 src\n"), "only appear once"),
        ];
        for (rows, word) in &cases {
            let error = VoiceDescription::from_bytes(fixture(rows).as_bytes()).unwrap_err();
            assert!(error.contains(word), "{rows:?}: {error}");
        }
        let keyframes = [0.5e-4, 0.7e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4]
            .iter()
            .map(|a| format!("keyframe_section b {a:e} 0.02\n"))
            .collect::<String>();
        let moving = |rows: &str| {
            fixture(&format!(
                "glottal_orifice 4.0e-5\nfricative_noise 5 1 4\n{rows}keyframe b source src\nkeyframe b licence test\n\
                 {keyframes}tract_move 0.01 0.02 b 2\n"
            ))
        };
        VoiceDescription::from_bytes(moving("").as_bytes()).expect("the moving card without the row parses");
        let error = VoiceDescription::from_bytes(moving("fricative_throat 6.5e-5 src\n").as_bytes()).unwrap_err();
        assert!(error.contains("refused with tract_move"), "{error}");
    }

    /// `fricative_throat` moves the constriction laws' `A_k` and nothing else (MB 25440 ruling A, 25482). The row at section
    /// `k`'s own area renders bit for bit the card without it. With a measured throat, the rig's held flow through the
    /// orifice and the constriction is the closed form at the throat to 1e-12: the root of `P = rho K U^2 / 2 + K3_eff U`,
    /// `K = K_o + (1/A_t - 1/A_j)^2` with the orifice's eq 6 factor `K_o`, `K3_eff` Shadle's K3 at the throat less the
    /// line's own Poiseuille resistance over section `k`'s span (a quantity the row does not give), and it is not the flow
    /// at section `k`'s area (the same closed form at section `k` is the card without the row, to 1e-12). The source's
    /// shape takes the throat's `A_k` and `sigma` and keeps the tract's obstacle `A_m` and `A_j`'s radius, and the throat
    /// reaches the radiated pressure.
    #[test]
    fn fricative_throat_moves_only_the_constriction_laws_area() {
        let bits = |xs: &[f64]| xs.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
        // Test chart, not data: section k = 1 narrowed to 0.2 cm^2, inside Shadle's printed rows, so K3 is read between two
        // rows and the line's span through k is section k alone.
        let narrow = |rows: &str| {
            let text = fixture(rows).replace("section 6e-5 0.02\n", "section 2e-5 0.02\n");
            assert_ne!(text, fixture(rows), "the replacement must change the card");
            text
        };
        let rig = "glottal_orifice 4.0e-5\nsubglottal_tube 2.5e-4 0.08\nreceiver free_field_monopole\n\
                   fricative_noise 5 1 4 stevens 3\n";
        let without = render(&narrow(rig));
        assert!(without.iter().all(|x| x.is_finite()));
        let at_k = render(&narrow(&format!("{rig}fricative_throat 2e-5 the-level-cut-itself\n")));
        assert_eq!(bits(&without), bits(&at_k), "the row at section k's own area is the card without it, bit for bit");
        let throat = 2.5e-5;
        let card = narrow(&format!("{rig}fricative_throat {throat:e} flat-section-test-chart\n"));
        let description = VoiceDescription::from_bytes(card.as_bytes()).expect("parses");
        let measured = description.build().expect("builds");
        let plain = VoiceDescription::from_bytes(narrow(rig).as_bytes()).expect("parses").build().expect("builds");
        let gas = description.tract_gas();
        let (rho, mu) = (gas.density, gas.dynamic_viscosity);
        let span = line_span(&description.sections, 1);
        assert_eq!(span, vec![(2.0e-5, 0.02)], "the span through k is section k alone");
        let k3_line = crate::glottis::line_poiseuille_resistance(mu, &span);
        let held = |a_t: f64, p: f64| {
            let (a_o, a_1, a_j) = (4.0e-5, 0.5e-4, 3.0e-4);
            let n = a_o / a_1;
            let gap = 1.0 / a_t - 1.0 / a_j;
            let kinetic = (1.0 - 2.0 * n * (1.0 - n)) / (a_o * a_o) + gap * gap;
            let k3 = (crate::glottis::shadle_constriction_k3(a_t).expect("K3") - k3_line).max(0.0);
            // The root of rho K U^2 / 2 + K3 U = P, rationalized: 2 P / (K3 + sqrt(K3^2 + 2 rho K P)).
            2.0 * p / (k3 + (k3 * k3 + 2.0 * rho * kinetic * p).sqrt())
        };
        for p in [100.0, 900.0] {
            let at_throat = measured.island.orifice_held_flow_m3_s(p).expect("the rig's held flow");
            let closed = held(throat, p);
            assert!((at_throat / closed - 1.0).abs() < 1.0e-12, "{p} Pa: island {at_throat:e}, closed form {closed:e}");
            let at_section = plain.island.orifice_held_flow_m3_s(p).expect("the rig's held flow");
            let closed_k = held(2.0e-5, p);
            assert!((at_section / closed_k - 1.0).abs() < 1.0e-12, "{p} Pa: island {at_section:e}, closed {closed_k:e}");
            assert!((at_throat / at_section - 1.0).abs() > 1.0e-3, "{p} Pa: the throat, not section k, sets the flow");
        }
        let (shape, own) = (
            measured.fricative.as_ref().expect("the source").shapes[0],
            plain.fricative.as_ref().expect("the source").shapes[0],
        );
        assert_eq!(shape.section, 1);
        assert_eq!(shape.area_k_m2.to_bits(), throat.to_bits(), "A_k: the throat");
        assert_eq!(shape.open_area_ratio.to_bits(), (throat / 3.0e-4).to_bits(), "sigma: the throat over A_j");
        assert_eq!(own.area_k_m2.to_bits(), 2.0e-5_f64.to_bits(), "without the row A_k is section k's");
        assert_eq!(
            [shape.area_m_m2, shape.duct_radius_m].map(f64::to_bits),
            [own.area_m_m2, own.duct_radius_m].map(f64::to_bits),
            "A_m and A_j's radius stay the tract's"
        );
        assert_eq!(shape.area_m_m2.to_bits(), 2.5e-4_f64.to_bits(), "A_m: the obstacle section 3");
        let spoken = render(&card);
        assert!(spoken.iter().all(|x| x.is_finite()));
        assert_ne!(bits(&spoken), bits(&without), "the throat reaches the radiated pressure");
    }

    /// The fricative source on the folds (MB 16735 piece 1): it renders finite, reaches the radiated pressure, one seed
    /// gives one noise and another seed another, and it reads the island's own flow: with no lung pressure there is no
    /// flow, no jet and no noise, so the card renders as the same card without the row's noise would (the constriction's
    /// Borda-Carnot drop and K3 term stay in the island's solve either way, and at zero flow they take nothing).
    #[test]
    fn fricative_on_the_folds_reads_the_islands_own_flow() {
        let voiced = render(&fixture(""));
        let spoken = render(&fixture("fricative_noise 5 1 4\n"));
        assert!(spoken.iter().all(|x| x.is_finite()));
        assert!(voiced.iter().zip(&spoken).any(|(a, b)| a.to_bits() != b.to_bits()), "the fricative must reach the pressure");
        assert_eq!(spoken, render(&fixture("fricative_noise 5 1 4\n")), "one seed, one noise");
        assert_ne!(spoken, render(&fixture("fricative_noise 6 1 4\n")), "another seed, another noise");
        let still = |rows: &str| fixture(rows).replace("subglottal 0.02 900", "subglottal 0.02 0");
        assert_ne!(still(""), fixture(""), "the replacement must change the card");
        let silent = render(&still("fricative_noise 5 1 4\n"));
        assert!(silent.iter().all(|x| *x == 0.0), "no lung pressure: no flow, no jet, no noise");
        assert_eq!(silent, render(&still("")), "a still fold card renders the same with the row");
    }

    /// The fricative source's lines end in the island's own load (MB 16735 piece 3): with `head_load` the card renders
    /// finite and the noise reaches the pressure through the head-load lines (a seed pair: the constriction's drop in the
    /// island's solve is the same for both seeds, so only the noise can tell them apart; NavyTower 16796); a `tract_gas`
    /// row holding the exterior's own record takes the exterior-checked lines with one gas, which are the head-load lines
    /// exactly, so the render is bit for bit the card's without the row; a different tract gas renders finite and its
    /// noise reaches the pressure through the exterior-checked lines (a seed pair there too).
    #[test]
    fn fricative_lines_take_the_islands_head_load_and_tract_gas() {
        let bits = |xs: &[f64]| xs.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
        let rows = |seed: u32| format!("{HEAD_LOAD}fricative_noise {seed} 1 4\n");
        let spoken = render(&fixture(&rows(5)));
        assert!(spoken.iter().all(|x| x.is_finite()));
        assert_ne!(bits(&spoken), bits(&render(&fixture(&rows(6)))), "the noise reaches the pressure through the head's load");
        let same = render(&fixture(&format!("{}tract_gas moist_air 310.15 101325 1.0\n", rows(5))));
        assert_eq!(bits(&spoken), bits(&same), "one gas: the exterior-checked lines are the head-load lines");
        let warm = |seed: u32| render(&fixture(&format!("{}tract_gas moist_air 290.15 101325 0.5\n", rows(seed))));
        let warm_5 = warm(5);
        assert!(warm_5.iter().all(|x| x.is_finite()));
        assert_ne!(bits(&warm_5), bits(&spoken), "another tract gas reaches the pressure");
        assert_ne!(bits(&warm_5), bits(&warm(6)), "the noise reaches the pressure through the exterior-checked lines");
    }

    /// The fricative source on a moving tract (MB 16735 piece 2). Every shape takes its own junction lines and its own
    /// constriction, `A_k` that shape's section `k` and `A_j` its largest section after `k` (MB 15828). Shape 0 is the
    /// card's own, bit for bit the fixed card's. On a move toward a tract whose section `k` is wider than every later
    /// section, the middle shape is still an expansion (`sigma < 1`) and the target is past `sigma = 1`, where the printed
    /// law drives no band: the source turns off through the law as the tract opens. The card renders finite, the noise
    /// reaches the pressure on the moving lines (a seed pair), and the render ends with the source on the target.
    #[test]
    fn fricative_on_a_moving_tract_takes_each_shapes_lines_and_constriction() {
        let bits = |xs: &[f64]| xs.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
        // Test chart, not data: section k = 1 opens from 0.6 cm^2 to 4.0 cm^2, past every later section (3.0 at most).
        let opened = [0.5e-4, 4.0e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4];
        let card = |seed: u32| {
            let mut rows = format!("fricative_noise {seed} 1 4\nkeyframe o source src\nkeyframe o licence test\n");
            for area in opened {
                rows.push_str(&format!("keyframe_section o {area:e} 0.02\n"));
            }
            rows.push_str("tract_move 0.01 0.02 o 2\n");
            fixture(&rows)
        };
        let description = VoiceDescription::from_bytes(card(5).as_bytes()).expect("parses");
        let mut performance = description.build().expect("builds");
        let fixed = VoiceDescription::from_bytes(fixture("fricative_noise 5 1 4\n").as_bytes())
            .expect("parses")
            .build()
            .expect("builds");
        let own = fixed.fricative.as_ref().expect("the fixed card's source").shapes[0];
        let source = performance.fricative.as_ref().expect("the source");
        assert_eq!(source.shapes.len(), performance.tract_switches.len() + 1, "one constriction per shape");
        assert_eq!(source.shapes.len(), 3, "the card and the move's two steps");
        let (first, middle, last) = (source.shapes[0], source.shapes[1], source.shapes[2]);
        let of = |s: FricativeShape| [s.area_k_m2, s.open_area_ratio, s.duct_radius_m].map(f64::to_bits);
        assert_eq!(of(first), of(own), "shape 0: the fixed card's own");
        assert_eq!([first.section, middle.section, last.section], [1, 1, 1], "equal counts and lengths: the card's k");
        assert!(first.open_area_ratio < middle.open_area_ratio && middle.open_area_ratio < 1.0, "{middle:?}");
        assert!((last.area_k_m2 / 4.0e-4 - 1.0).abs() < 1.0e-12, "{last:?}");
        assert!((last.duct_radius_m / (3.0e-4 / core::f64::consts::PI).sqrt() - 1.0).abs() < 1.0e-12, "{last:?}");
        assert!(last.open_area_ratio > 1.0, "the target is past the expansion: {last:?}");
        assert!(source.noise.driven_bands_jet_in(100.0, 30.0, last.open_area_ratio, last.duct_radius_m).is_empty());
        let mut spoken = vec![0.0; description.samples as usize];
        for block in spoken.chunks_mut(4096) {
            performance.block(block).expect("renders");
        }
        assert!(spoken.iter().all(|x| x.is_finite()));
        assert_eq!(performance.fricative.as_ref().expect("the source").shape, 2, "the render ends on the target");
        assert_eq!(bits(&spoken), bits(&render(&card(5))), "one seed, one noise");
        assert_ne!(bits(&spoken), bits(&render(&card(6))), "the noise reaches the pressure on the moving lines");
        // A word's card is its vowel (NavyTower 16887): the card's own shape need not be an expansion, and a keyframe of
        // another section count takes its section at the card's k's centre. Test chart, not data: card k = 5 (3.0 cm^2,
        // nothing wider after it; j = 6, the widest after it), toward a 10-section keyframe of the card's length that
        // narrows to 0.2 cm^2 at section 6, where the card's k's centre (0.6875 of the length) falls.
        let narrowing = [0.5e-4, 0.6e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 0.2e-4, 3.0e-4, 2.5e-4, 2.0e-4];
        let closing = |seed: u32| {
            let mut rows = format!("fricative_noise {seed} 5 6\nkeyframe n source src\nkeyframe n licence test\n");
            for area in narrowing {
                rows.push_str(&format!("keyframe_section n {area:e} 0.016\n"));
            }
            rows.push_str("tract_move 0.01 0.02 n 2\n");
            fixture(&rows)
        };
        let closed = VoiceDescription::from_bytes(closing(5).as_bytes()).expect("parses").build().expect("builds");
        let source = closed.fricative.as_ref().expect("the source");
        let (vowel, target) = (source.shapes[0], source.shapes[2]);
        assert_eq!(vowel.section, 5);
        assert!(vowel.open_area_ratio > 1.0, "the card's vowel is no expansion after k: {vowel:?}");
        assert_eq!(target.section, 6, "the card's k's centre falls on the keyframe's section 6: {target:?}");
        assert!((target.open_area_ratio / (0.2 / 3.0) - 1.0).abs() < 1.0e-12, "{target:?}");
        let closing_5 = render(&closing(5));
        assert!(closing_5.iter().all(|x| x.is_finite()));
        assert_ne!(bits(&closing_5), bits(&render(&closing(6))), "the noise forms as the tract closes");
        // MB 16891's condition on (c): at the card's section count the rule gives k itself, so every equal-count card
        // stays bit for bit. Test chart, not data: the opened areas on unequal lengths of the card's total, where the
        // card's k's centre (0.1875 of the length) falls in the keyframe's section 0; every shape still takes k = 1.
        let uneven_length = |i: usize| if i == 0 { 0.05 } else { 0.11 / 7.0 };
        let uneven = |seed: u32| {
            let mut rows = format!("fricative_noise {seed} 1 4\nkeyframe u source src\nkeyframe u licence test\n");
            for (i, area) in opened.iter().enumerate() {
                rows.push_str(&format!("keyframe_section u {area:e} {:e}\n", uneven_length(i)));
            }
            rows.push_str("tract_move 0.01 0.02 u 2\n");
            fixture(&rows)
        };
        let keyframe: Vec<TractSection> = opened
            .iter()
            .enumerate()
            .map(|(i, &area_m2)| TractSection { area_m2, length_m: uneven_length(i) })
            .collect();
        let at_equal_count = VoiceDescription::from_bytes(uneven(5).as_bytes()).expect("parses");
        let card_sections = &at_equal_count.sections;
        assert_eq!(card_sections.len(), keyframe.len());
        assert_eq!(section_at_fraction(&keyframe, centre_fraction(card_sections, 1)), 0, "the fraction alone moves k");
        assert_eq!(shape_constriction(card_sections, 1, &keyframe), 1, "an equal count gives k itself");
        let uneven_source = at_equal_count.build().expect("builds").fricative.expect("the source");
        assert_eq!(uneven_source.shapes.iter().map(|s| s.section).collect::<Vec<_>>(), [1, 1, 1], "k on every shape");
        assert_eq!(of(uneven_source.shapes[0]), of(own), "shape 0: the fixed card's own");
        assert!(render(&uneven(5)).iter().all(|x| x.is_finite()));
    }

    /// The obstacle section `m` on a moving tract (MB 25469): every shape takes the card's `m` by the correspondence
    /// it takes `k` by, its junction lines after that section and its `A_m` that section's area. At the card's
    /// section count that is `m` itself; at another count, the shape's section at the card's `m`'s centre fraction.
    /// With `m = k` every shape's `A_m` is its `A_k` bit for bit, as before `m` could move. The source after `m` is
    /// not the source after `k`, and a shape whose `m` falls on its last section is refused by name at build.
    #[test]
    fn fricative_obstacle_on_a_moving_tract_is_each_shapes_own_m() {
        let bits = |xs: &[f64]| xs.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
        // Test chart, not data: the card's k = 1, j = 4, m = 3 (2.5 cm^2); the keyframe narrows section 3 to 2.2 cm^2.
        let card = |seed: u32, tail: &str, keyframe: &[(f64, f64)]| {
            let mut rows = format!("fricative_noise {seed} 1 4 stevens{tail}\n");
            rows.push_str("keyframe b source src\nkeyframe b licence test\n");
            for (area, length) in keyframe {
                rows.push_str(&format!("keyframe_section b {area:e} {length:e}\n"));
            }
            rows.push_str("tract_move 0.01 0.02 b 2\n");
            fixture(&rows)
        };
        let equal: Vec<(f64, f64)> =
            [0.5e-4, 0.7e-4, 1.5e-4, 2.2e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4].iter().map(|&a| (a, 0.02)).collect();
        let shapes = |text: &str| {
            VoiceDescription::from_bytes(text.as_bytes())
                .expect("parses")
                .build()
                .expect("builds")
                .fricative
                .expect("the source")
                .shapes
        };
        let at_m = shapes(&card(5, " 3", &equal));
        assert_eq!(at_m.len(), 3, "the card and the move's two steps");
        assert_eq!(at_m[0].area_m_m2.to_bits(), 2.5e-4_f64.to_bits(), "shape 0: the card's own section m");
        assert!((at_m[2].area_m_m2 / 2.2e-4 - 1.0).abs() < 1.0e-12, "the target's own section m: {:?}", at_m[2]);
        assert!((at_m[2].area_k_m2 / 0.7e-4 - 1.0).abs() < 1.0e-12, "k stays the target's section 1: {:?}", at_m[2]);
        assert_eq!(at_m.iter().map(|s| s.section).collect::<Vec<_>>(), [1, 1, 1], "k on every shape");
        let at_k = shapes(&card(5, "", &equal));
        for (s, (with_m, without)) in at_m.iter().zip(&at_k).enumerate() {
            assert_eq!(without.area_m_m2.to_bits(), without.area_k_m2.to_bits(), "m = k: shape {s}'s A_m is its A_k");
            let law = |x: &FricativeShape| [x.area_k_m2, x.open_area_ratio, x.duct_radius_m].map(f64::to_bits);
            assert_eq!(law(with_m), law(without), "shape {s}: m moves neither A_k, sigma nor A_j's radius");
        }
        let spoken = render(&card(5, " 3", &equal));
        assert!(spoken.iter().all(|x| x.is_finite()));
        assert_ne!(bits(&spoken), bits(&render(&card(6, " 3", &equal))), "the noise reaches the pressure from m");
        assert_ne!(bits(&spoken), bits(&render(&card(5, "", &equal))), "the source after m is not the source after k");
        // Another section count, the card's length: ten sections of 1.6 cm. The card's m's centre (0.4375 of the
        // length, 7.0 cm) falls in the keyframe's section 4 and its k's (0.1875, 3.0 cm) in section 1.
        let recounted: Vec<(f64, f64)> =
            [0.5e-4, 0.7e-4, 1.5e-4, 2.0e-4, 2.2e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4, 1.8e-4]
                .iter()
                .map(|&a| (a, 0.016))
                .collect();
        let description = VoiceDescription::from_bytes(card(5, " 3", &recounted).as_bytes()).expect("parses");
        let keyframe: Vec<TractSection> =
            recounted.iter().map(|&(area_m2, length_m)| TractSection { area_m2, length_m }).collect();
        assert_eq!(shape_constriction(&description.sections, 3, &keyframe), 4, "the card's m's centre fraction");
        assert_eq!(shape_constriction(&description.sections, 1, &keyframe), 1, "the card's k's centre fraction");
        let target = description.build().expect("builds").fricative.expect("the source").shapes[2];
        assert!((target.area_m_m2 / 2.2e-4 - 1.0).abs() < 1.0e-12, "the keyframe's section 4: {target:?}");
        assert!((target.area_k_m2 / 0.7e-4 - 1.0).abs() < 1.0e-12, "the keyframe's section 1: {target:?}");
        // Seven sections, six of 2 cm and a last of 4 cm: the card's m = 6's centre (0.8125 of the length, 13 cm)
        // falls on the keyframe's last section, while its k's (3 cm) falls on section 1.
        let short: Vec<(f64, f64)> = [0.5e-4, 0.7e-4, 1.5e-4, 2.5e-4, 3.0e-4, 2.5e-4, 2.0e-4]
            .iter()
            .enumerate()
            .map(|(i, &a)| (a, if i == 6 { 0.04 } else { 0.02 }))
            .collect();
        let error = VoiceDescription::from_bytes(card(5, " 6", &short).as_bytes())
            .expect("parses: the card's own m is before its last section")
            .build()
            .err()
            .expect("refused at build");
        assert!(error.contains("obstacle section m falls on a moving shape's last section"), "{error}");
    }

    /// The free-field monopole: a low-frequency sinusoidal flow `Q0 sin(w t)` reaches the microphone at distance `R` as
    /// `rho w Q0 / (4 pi R)` once the filter latency and the travel time have passed, half the baffled point source.
    #[test]
    fn free_field_monopole_radiates_rho_dq_dt_over_four_pi_r() {
        let (rate, f, q0) = (40_000_u32, 500.0, 1.0e-4);
        let medium = RayleighMedium { density: 1.1886, sound_speed: 344.4 };
        let microphone = CircularOutletReceiver {
            position_m: [0.3, 0.0, 0.3],
            radial_rings: 4,
            angular_points: 16,
            maximum_frequency_hz: 16_000.0,
        };
        let mut receiver = free_field_monopole(3.06e-4, rate, rate, microphone, medium).expect("monopole");
        let w = core::f64::consts::TAU * f;
        let out: Vec<f64> = (0..8000)
            .map(|n| receiver.step_trace(&[q0 * (w * f64::from(n) / f64::from(rate)).sin()]).expect("step"))
            .collect();
        let peak = out[4000..].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let r = (0.3_f64 * 0.3 + 0.3 * 0.3).sqrt();
        let expected = medium.density * w * q0 / (4.0 * core::f64::consts::PI * r);
        assert!((peak / expected - 1.0).abs() < 0.01, "peak {peak} against rho w Q0 / (4 pi R) = {expected}");
        let wide = CircularOutletReceiver { maximum_frequency_hz: 18_001.0, ..microphone };
        assert!(free_field_monopole(3.06e-4, rate, rate, wide, medium).is_err(), "the disk receiver's band rule");
    }
}
