//! One struck bar as a route: card -> TET10 modes -> Burton-Miller receiver transfers -> cited-mallet strike -> WAV.
//!
//! cargo run --release -p fs-couple --example bar_strike -- CARD.bar OUT_DIR   (from the repository root)
//!
//! The bar is data (see `card.rs`); the physics are the existing owners the bar_s1 test composes. Nothing
//! here is fitted to a recording. What the card does not give is reported as absent: no resonator unless
//! the card names one (see `resonator.rs` for what that model omits), no supports, no mallet shaft, and
//! the wood's loss is the declared scenario named in the card.
//! Gauntlet: G2 (benchmarks) for the uniform-prism and density-scaling tests; no G3/G4/G5 claim.
use fs_bem::helmholtz::Medium;
use fs_couple::modal_acoustic_time::broadband::BroadbandModalObserver;
use fs_matdb::ContactLawCard;
use fs_material::gas_species::input::read_gas_state;
use fs_modal::ModePair;
use std::path::Path;
use std::time::Instant;

mod bar;
mod card;
mod observer;
mod radiation;
mod render;
mod resonator;
#[cfg(test)]
mod tests;
mod tube;

use card::Card;

type Error = Box<dyn std::error::Error>;

/// The free body at one mesh rung: band modes and the six rigid coordinates.
pub struct Structure {
    /// The structured P1 grid the TET10 elements were promoted from (its vertices come first).
    pub mesh: bar::Mesh,
    /// TET10 element count.
    pub n_tets: usize,
    /// Free degrees of freedom of the assembly (no constraints: the bar is free-free).
    pub free_dofs: usize,
    /// Every certified elastic mode from 20 Hz to the card's band, mass-normalized.
    pub modes: Vec<ModePair>,
    /// The six rigid coordinates' (negated) displacement at the strike node along the thickness axis.
    pub rigid: Vec<f64>,
    /// The same coordinates' mass-normalized fields, laid out as a mode's `phi`.
    pub rigid_fields: Vec<Vec<f64>>,
    /// Wall time of assembly, slicing and the rigid projections [s].
    pub wall_s: f64,
}

/// Assemble the bar at mesh rung `rung` with the given density [kg/m3] and slice its band modes.
/// The density is explicit so a test can plant a wrong one; the route passes the card's.
pub fn structure(card: &Card, rung: usize, density: f64) -> Result<Structure, String> {
    let started = Instant::now();
    bar::with_cx(|cx| {
        let mesh = bar::mesh(card, bar::rung_counts(card, rung));
        let mat = bar::material(card, density)?;
        let (nodes, n_tets, a) = bar::assemble(card, &mesh, &mat, cx)?;
        let modes = bar::band_modes(card, &a, cx)?;
        let (rigid, rigid_fields) = bar::rigid_projectors(card, &mesh, &nodes, &a)?;
        Ok(Structure {
            n_tets,
            free_dofs: a.free_dofs.len(),
            mesh,
            modes,
            rigid,
            rigid_fields,
            wall_s: started.elapsed().as_secs_f64(),
        })
    })
}

fn main() -> Result<(), Error> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage = "usage: bar_strike CARD.bar OUT_DIR | bar_strike CARD.bar --modes-only | bar_strike CARD.bar --observer-cost | bar_strike CARD.bar --surface-check X,Y,Z ... (run from the repository root)";
    let surface_check = args.len() > 2 && args[1] == "--surface-check";
    if args.len() != 2 && !surface_check {
        return Err(usage.into());
    }
    let root = std::env::current_dir()?;
    let card = Card::load(Path::new(&args[0]), &root)?;
    if args[1] == "--modes-only" {
        // Geometry checks: the structural rungs only, no radiation and no strike.
        for &rung in &card.mesh_rungs {
            report_structure(&card, rung, &structure(&card, rung, card.density_kg_m3)?);
        }
        return Ok(());
    }
    if args[1] == "--observer-cost" {
        return observer_cost(&card);
    }
    if surface_check {
        return surface_check_run(&card, &args[2..]);
    }
    run(&card, Path::new(&args[1]))
}

/// How the receiver transfers move with the radiating surface: every mode at its own frequency, on the card's
/// finest rung's modes, for each named stride of that rung's grid along length, width and thickness. Nothing
/// is struck and nothing is written.
fn surface_check_run(card: &Card, strides: &[String]) -> Result<(), Error> {
    let mut finest = None;
    for &rung in &card.mesh_rungs {
        let s = structure(card, rung, card.density_kg_m3)?;
        report_structure(card, rung, &s);
        finest = Some(s);
    }
    let s = finest.ok_or("the card lists no mesh rung")?;
    let gas_text = std::fs::read_to_string(&card.gas_path)
        .map_err(|e| format!("{}: {e}", card.gas_path.display()))?;
    let gas = read_gas_state(&mut gas_text.split_ascii_whitespace(), None)
        .map_err(|e| format!("gas record: {e:?}"))?;
    let medium = Medium::from(&gas);
    for text in strides {
        let parsed = text
            .split(',')
            .map(str::parse::<usize>)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format!("stride `{text}`: {e}"))?;
        let stride: [usize; 3] = parsed
            .try_into()
            .map_err(|_| format!("stride `{text}` needs three counts"))?;
        if (0..3).any(|a| stride[a] == 0 || s.mesh.counts[a] % stride[a] != 0) {
            return Err(format!(
                "stride `{text}` does not divide the grid {:?}",
                s.mesh.counts
            )
            .into());
        }
        let started = Instant::now();
        let (surface, shapes) = radiation::surface_and_shapes_at(card, &s.mesh, &s.modes, stride)?;
        let t = radiation::modal_transfers(card, medium, &surface, &shapes, &s.modes)?;
        let row: Vec<[f64; 2]> = t.pressure.iter().map(|p| [p.re, p.im]).collect();
        println!(
            "{{\"suite\":\"bar-strike\",\"case\":\"surface-check\",\"bar\":\"{}\",\"grid\":{:?},\"stride\":{stride:?},\"panels\":{},\"min_panels_per_wavelength\":{},\"p_pa_per_unit_modal_velocity_re_im\":{row:?},\"wall_s\":{}}}",
            card.name,
            s.mesh.counts,
            t.panels,
            t.min_panels_per_wavelength,
            started.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

/// What the broadband observer would cost on this card and how its fit behaves on this surface: the band's
/// BEM, every fit order and the time of a step. Nothing is struck and nothing is written.
fn observer_cost(card: &Card) -> Result<(), Error> {
    let mut finest = None;
    for &rung in &card.mesh_rungs {
        let s = structure(card, rung, card.density_kg_m3)?;
        report_structure(card, rung, &s);
        finest = Some(s);
    }
    let s = finest.ok_or("the card lists no mesh rung")?;
    let gas_text = std::fs::read_to_string(&card.gas_path)
        .map_err(|e| format!("{}: {e}", card.gas_path.display()))?;
    let gas = read_gas_state(&mut gas_text.split_ascii_whitespace(), None)
        .map_err(|e| format!("gas record: {e:?}"))?;
    let medium = Medium::from(&gas);
    let (surface, shapes) = radiation::surface_and_shapes(card, &s.mesh, &s.modes)?;
    let started = Instant::now();
    let transfers = radiation::modal_transfers(card, medium, &surface, &shapes, &s.modes)?;
    println!(
        "{{\"suite\":\"bar-strike\",\"case\":\"observer-narrow-baseline\",\"bar\":\"{}\",\"panels\":{},\"modes\":{},\"min_panels_per_wavelength\":{},\"wall_s\":{}}}",
        card.name,
        transfers.panels,
        s.modes.len(),
        transfers.min_panels_per_wavelength,
        started.elapsed().as_secs_f64()
    );
    let band = observer::band(card, medium, &surface, &shapes, &s, &transfers)?;
    observer::report_band(card, &band);
    let modes = render::bar_modes(&s.modes, None, None, &transfers.pressure, s.rigid.len());
    let direct = render::DirectBody::new(card, modes)?;
    let own_hz: Vec<f64> = s.modes.iter().map(bar::hz).collect();
    observer::report_fits(card, direct.model(), &band, &transfers.pressure, &own_hz);
    Ok(())
}

fn run(card: &Card, out: &Path) -> Result<(), Error> {
    // `create_dir` refuses an existing directory: a render never overwrites an earlier one.
    std::fs::create_dir(out)?;
    let mut finest = None;
    for &rung in &card.mesh_rungs {
        let s = structure(card, rung, card.density_kg_m3)?;
        report_structure(card, rung, &s);
        finest = Some(s);
    }
    let s = finest.ok_or("the card lists no mesh rung")?;
    let r = strike(card, &s)?;
    let transfer_row: Vec<[f64; 2]> = r.transfers.pressure.iter().map(|p| [p.re, p.im]).collect();
    println!(
        "{{\"suite\":\"bar-strike\",\"case\":\"burton-miller-modal-receiver\",\"bar\":\"{}\",\"panels\":{},\"min_panels_per_wavelength\":{},\"receiver_m\":{:?},\"p_pa_per_unit_modal_velocity_re_im\":{transfer_row:?},\"gas_sound_speed_m_s\":{},\"elapsed_s\":{}}}",
        card.name,
        r.transfers.panels,
        r.transfers.min_panels_per_wavelength,
        card.receiver_m,
        r.sound_speed_m_s,
        r.radiation_wall_s
    );
    let resonator = r
        .resonator_isolated_hz
        .map_or("null".to_string(), |f| f.to_string());
    let (stroke, mallet, loss, validity) = (&r.stroke, &r.mallet, &r.loss, &r.validity);
    let rms =
        (stroke.samples.iter().map(|p| p * p).sum::<f64>() / stroke.samples.len() as f64).sqrt();
    let peak = stroke.samples.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
    let (wav, clipped) =
        fs_couple::pcm_wav::encode_pcm16_wav(&stroke.samples, card.rate_hz, card.full_scale_pa)
            .map_err(|e| format!("pcm16 wav: {e:?}"))?;
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav).to_hex();
    std::fs::write(out.join(format!("{}.wav", card.name)), wav)?;
    let f32_wav = fs_couple::pcm_wav::high_resolution::encode_pressure_wav(
        &stroke.samples,
        card.rate_hz,
        1,
        card.full_scale_pa,
        fs_couple::pcm_wav::high_resolution::Encoding::Float32,
    )
    .map_err(|e| format!("float32 wav: {e:?}"))?;
    std::fs::write(out.join(format!("{}-f32.wav", card.name)), f32_wav)?;
    if let Some(o) = &r.observer {
        // The same run's narrow-band receiver, kept beside the broadband one: it must replay a render
        // without the observer bit for bit, which is the check that the bar did not change.
        let narrow = fs_couple::pcm_wav::high_resolution::encode_pressure_wav(
            &r.narrow_samples,
            card.rate_hz,
            1,
            card.full_scale_pa,
            fs_couple::pcm_wav::high_resolution::Encoding::Float32,
        )
        .map_err(|e| format!("float32 wav: {e:?}"))?;
        std::fs::write(out.join(format!("{}-narrow-f32.wav", card.name)), narrow)?;
        println!(
            "{{\"suite\":\"bar-strike\",\"case\":\"broadband-observer\",\"bar\":\"{}\",\"order\":{},\"held_out_max\":{},\"held_out_budget\":{},\"producer_estimate\":{},\"delay_samples\":{},\"top_hz\":{},\"training\":{},\"held_out\":{},\"band_wall_s\":{},\"fit_wall_s\":{},\"source\":\"{}\"}}",
            card.name,
            o.report.controls.fit_order,
            o.report.held_out_max_rel_error,
            o.report.controls.max_held_out_rel_error,
            o.report.producer_relative_error_estimate,
            o.report.delay_samples,
            o.top_hz,
            o.training,
            o.held_out,
            o.band_wall_s,
            o.fit_wall_s,
            o.report.source_id
        );
    }
    // An outside card's refusal is its Debug text, which quotes its fields: escaped for the JSON string
    // (backslash first, then the double quote) the row still parses. `inside` has neither and is unchanged.
    let validity = validity.replace('\\', "\\\\").replace('"', "\\\"");
    // Under `bar_self_radiation on`: every band mode's W in 1/s, negatives included, and the floor they set
    // (the largest negative |W|; a positive W of that size is no better resolved). Absent, the record is
    // as before, byte for byte.
    let self_radiation = if card.bar_self_radiation {
        let w = &r.transfers.self_radiation;
        let floor = w.iter().fold(0.0_f64, |m, &x| m.max(-x));
        format!(",\"bar_self_radiation_per_s\":{w:?},\"bar_self_radiation_floor_per_s\":{floor}")
    } else {
        String::new()
    };
    let record = format!(
        "{{\"suite\":\"bar-strike\",\"case\":\"strike-render\",\"bar\":\"{}\",\"mallet_mass_kg\":{},\"mallet_speed_m_s\":{},\"k_n_per_m_alpha\":{},\"alpha\":{},\"contact_time_s\":{},\"peak_force_n\":{},\"contact_card_validity\":\"{validity}\",\"strike_x_m\":{},\"loss\":\"{}\"{self_radiation},\"rms_pa\":{rms},\"peak_pa\":{peak},\"clipped_samples\":{clipped},\"wav_blake3\":\"{hash}\",\"sample_rate_hz\":{},\"full_scale_pa\":{},\"resonator_isolated_hz\":{resonator},\"supports\":false,\"mallet_shaft\":false,\"acceptance\":false,\"elapsed_s\":{}}}",
        card.name,
        mallet.mass_kg,
        mallet.speed_m_s,
        mallet.k_n_per_m_alpha,
        mallet.alpha,
        stroke.contact_time_s,
        stroke.peak_force_n,
        bar::strike_x(card, &s.mesh),
        if loss.is_some() {
            "declared scenario from the card"
        } else {
            "none (lossless diagnostic)"
        },
        card.rate_hz,
        card.full_scale_pa,
        r.render_wall_s
    );
    println!("{record}");
    std::fs::write(out.join(format!("{}.json", card.name)), record + "\n")?;
    Ok(())
}

/// Everything the strike produces for one bar: the receiver transfers, the stroke, the cited mallet and
/// the contact card's verdict on the stroke's peak force and speed.
pub struct Strike {
    /// Burton-Miller receiver pressure per unit modal velocity, one per elastic mode.
    pub transfers: radiation::Transfers,
    /// The rendered pressure history and its contact summary.
    pub stroke: render::Stroke,
    /// The cited mallet that struck it.
    pub mallet: render::Mallet,
    /// The wood's declared decay law, when the card names one.
    pub loss: Option<render::Loss>,
    /// `inside`, or the contact card's refusal text for the stroke's peak force and speed.
    pub validity: String,
    /// Isolated resonance of the card's resonator, when it has one.
    pub resonator_isolated_hz: Option<f64>,
    /// The broadband observer's fit and cost, when the card names one.
    pub observer: Option<ObserverOutcome>,
    /// The narrow-band receiver of the same run when the broadband observer was the one heard.
    pub narrow_samples: Vec<f64>,
    /// Sound speed of the card's gas record [m/s].
    pub sound_speed_m_s: f64,
    /// Wall time of the surface and the modal transfers [s].
    pub radiation_wall_s: f64,
    /// Wall time of the strike render [s].
    pub render_wall_s: f64,
}

/// What the observer row cost and what its fit reported.
pub struct ObserverOutcome {
    /// The broadband fit's own report.
    pub report: fs_couple::modal_acoustic_time::broadband::BroadbandModalReport,
    /// Top of the band the observer was fitted on [Hz].
    pub top_hz: f64,
    /// Number of frequencies the fit was trained on.
    pub training: usize,
    /// Number of frequencies held out to check the fit.
    pub held_out: usize,
    /// Wall time of the whole band [s].
    pub band_wall_s: f64,
    /// Wall time of the fit [s].
    pub fit_wall_s: f64,
    /// The whole bar's translation against a compact dipole, when the rigid rows are on.
    pub rigid_check: Option<observer::RigidCheck>,
}

/// The broadband observer fitted on `model`, the bar's modal model as the render steps it, with what its band
/// and fit cost. A bar alone and a bar over its tube build it the same way.
fn broadband_observer(
    card: &Card,
    row: card::ObserverRow,
    medium: Medium,
    (surface, shapes): (&fs_bem::panel3d::SpherePanels, &[Vec<f64>]),
    s: &Structure,
    transfers: &radiation::Transfers,
    model: &fs_couple::modal_acoustic_time::ModalAcousticTimeModel,
) -> Result<(BroadbandModalObserver, ObserverOutcome), String> {
    let band = observer::band(card, medium, surface, shapes, s, transfers)?;
    observer::report_band(card, &band);
    let fit_started = Instant::now();
    let (heard, report) = BroadbandModalObserver::build(
        model,
        &band.samples,
        observer::controls(row.fit_order, row.held_out_budget),
        card.rate_hz * card.substeps,
    )
    .map_err(|e| format!("broadband observer: {e:?}"))?;
    let outcome = ObserverOutcome {
        report,
        top_hz: band.top_hz,
        training: band.samples.omega_rad_s.len(),
        held_out: band.samples.held_out.len(),
        band_wall_s: band.wall_s,
        fit_wall_s: fit_started.elapsed().as_secs_f64(),
        rigid_check: band.rigid_check,
    };
    Ok((heard, outcome))
}

/// Radiate every mode of `s` to the card's receiver and render the cited mallet's stroke.
#[allow(clippy::too_many_lines)] // one coherent realization dispatcher
pub fn strike(card: &Card, s: &Structure) -> Result<Strike, String> {
    let gas_text = std::fs::read_to_string(&card.gas_path)
        .map_err(|e| format!("{}: {e}", card.gas_path.display()))?;
    let gas = read_gas_state(&mut gas_text.split_ascii_whitespace(), None)
        .map_err(|e| format!("gas record: {e:?}"))?;
    let medium = Medium::from(&gas);
    let read = |p: &Path| std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()));
    let loss = card
        .loss_path
        .as_ref()
        .map(|p| read(p).and_then(|text| render::Loss::read(&text)))
        .transpose()?;
    if card.bar_self_radiation {
        render::self_radiation_loss(loss.as_ref())?;
    }

    let started = Instant::now();
    let (surface, shapes) = radiation::surface_and_shapes(card, &s.mesh, &s.modes)?;
    let transfers = radiation::modal_transfers(card, medium, &surface, &shapes, &s.modes)?;
    let radiation_wall_s = started.elapsed().as_secs_f64();

    let contact =
        ContactLawCard::parse(&card.mallet_contact_pack, &read(&card.mallet_contact_path)?)
            .map_err(|e| format!("mallet contact card: {e:?}"))?;
    let mallet = render::Mallet {
        mass_kg: render::scalar(
            &read(&card.mallet_properties_path)?,
            "mallet_equivalent_mass",
        )?,
        speed_m_s: card.mallet_speed_m_s,
        k_n_per_m_alpha: contact.k_n_per_m_alpha,
        alpha: contact.alpha,
    };
    let node = bar::strike_node(card, &s.mesh);
    let mut projection: Vec<f64> = s.modes.iter().map(|m| -m.phi[3 * node + 2]).collect();
    projection.extend_from_slice(&s.rigid);
    let self_radiation = if card.bar_self_radiation {
        if !transfers.self_radiation.iter().all(|w| w.is_finite()) {
            return Err("non-finite bar self-radiation".into());
        }
        Some(transfers.self_radiation.as_slice())
    } else {
        None
    };
    let modes = render::bar_modes(
        &s.modes,
        loss.as_ref(),
        self_radiation,
        &transfers.pressure,
        s.rigid.len(),
    );

    let mut started = Instant::now();
    let mut resonator_isolated_hz = None;
    let mut observer_outcome = None;
    let mut body: Box<dyn render::Body> = match (&card.resonator, &card.tube) {
        (None, None) => {
            let direct = render::DirectBody::new(card, modes)?;
            match card.observer {
                None => Box::new(direct),
                Some(row) => {
                    let (heard, outcome) = broadband_observer(
                        card,
                        row,
                        medium,
                        (&surface, &shapes),
                        s,
                        &transfers,
                        direct.model(),
                    )?;
                    observer_outcome = Some(outcome);
                    started = Instant::now();
                    Box::new(render::ObservedBody::new(direct, heard))
                }
            }
        }
        (Some(res), _) => {
            resonator_isolated_hz = Some(resonator::isolated_hz(&resonator::cavity(res, &gas)?));
            let shapes = resonator::wet_volume_shapes(&s.mesh, &s.modes, res.wet_x_m);
            Box::new(resonator::CavityBody::new(card, res, &gas, modes, shapes)?)
        }
        (None, Some(t)) => {
            let isolated = tube::isolated_modes(t, &gas)?;
            let fields = radiation::mouth_fields(card, t, medium, &surface, &shapes, &isolated)?;
            let bar_omega: Vec<f64> = modes[..s.modes.len()]
                .iter()
                .map(|m| m.angular_frequency_rad_s)
                .collect();
            let network = tube::Network::new(&gas, t, &isolated, &fields, &bar_omega)?;
            tube::report(card, t, &isolated, &fields, &network);
            let (monopole, delay) = tube::monopole_to_receiver(card, t, &gas);
            let mut over_tube = tube::TubeBody::new(
                card.rate_hz * card.substeps,
                &modes,
                s.modes.len(),
                &network,
                monopole,
                delay,
            )?;
            if let Some(row) = card.observer {
                let (heard, outcome) = broadband_observer(
                    card,
                    row,
                    medium,
                    (&surface, &shapes),
                    s,
                    &transfers,
                    over_tube.model(),
                )?;
                observer_outcome = Some(outcome);
                started = Instant::now();
                over_tube.observe(heard);
            }
            Box::new(over_tube)
        }
    };
    let stroke = render::render(card, body.as_mut(), &projection, &mallet)?;
    let narrow_samples = body.take_narrow();
    let validity = contact
        .lookup(stroke.peak_force_n, mallet.speed_m_s)
        .map(|_| "inside".to_string())
        .unwrap_or_else(|e| format!("{e:?}"));
    Ok(Strike {
        transfers,
        stroke,
        mallet,
        loss,
        validity,
        resonator_isolated_hz,
        observer: observer_outcome,
        narrow_samples,
        sound_speed_m_s: medium.sound_speed,
        radiation_wall_s,
        render_wall_s: started.elapsed().as_secs_f64(),
    })
}

/// One JSON row per rung: every band mode, its certified uncertainty, and the vertical bending modes
/// against the printed frequencies of the same bar (a report, never an input).
fn report_structure(card: &Card, rung: usize, s: &Structure) {
    let f: Vec<f64> = s.modes.iter().map(bar::hz).collect();
    let e: Vec<f64> = s.modes.iter().map(bar::hz_error).collect();
    let vertical: Vec<f64> = bar::vertical_modes(&s.mesh, &s.modes)
        .iter()
        .map(|&i| f[i])
        .collect();
    let misses: Vec<f64> = vertical
        .iter()
        .zip(&card.reference_vertical_hz)
        .map(|(a, b)| 100.0 * (a / b - 1.0))
        .collect();
    println!(
        "{{\"suite\":\"bar-strike\",\"case\":\"structural-rung\",\"bar\":\"{}\",\"element\":\"TET10\",\"rung\":{rung},\"counts\":{:?},\"n_tets\":{},\"free_dofs\":{},\"band_hz\":[20,{}],\"mode_count\":{},\"f_modes_hz\":{f:?},\"eigen_uncertainty_hz\":{e:?},\"vertical_hz\":{vertical:?},\"reference_vertical_hz\":{:?},\"vertical_miss_percent\":{misses:?},\"wall_s\":{}}}",
        card.name,
        s.mesh.counts,
        s.n_tets,
        s.free_dofs,
        card.band_hz,
        f.len(),
        card.reference_vertical_hz,
        s.wall_s
    );
}
