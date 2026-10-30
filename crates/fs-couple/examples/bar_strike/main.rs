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
use fs_matdb::ContactLawCard;
use fs_material::gas_species::input::read_gas_state;
use fs_modal::ModePair;
use std::path::Path;
use std::time::Instant;

mod bar;
mod card;
mod radiation;
mod render;
mod resonator;
#[cfg(test)]
mod tests;

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
        let rigid = bar::rigid_projectors(card, &mesh, &nodes, &a)?;
        Ok(Structure {
            n_tets,
            free_dofs: a.free_dofs.len(),
            mesh,
            modes,
            rigid,
            wall_s: started.elapsed().as_secs_f64(),
        })
    })
}

fn main() -> Result<(), Error> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage = "usage: bar_strike CARD.bar OUT_DIR | bar_strike CARD.bar --modes-only (run from the repository root)";
    if args.len() != 2 {
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
    run(&card, Path::new(&args[1]))
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
    let (stroke, mallet, loss, validity) = (r.stroke, r.mallet, r.loss, r.validity);
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
    let record = format!(
        "{{\"suite\":\"bar-strike\",\"case\":\"strike-render\",\"bar\":\"{}\",\"mallet_mass_kg\":{},\"mallet_speed_m_s\":{},\"k_n_per_m_alpha\":{},\"alpha\":{},\"contact_time_s\":{},\"peak_force_n\":{},\"contact_card_validity\":\"{validity}\",\"strike_x_m\":{},\"loss\":\"{}\",\"rms_pa\":{rms},\"peak_pa\":{peak},\"clipped_samples\":{clipped},\"wav_blake3\":\"{hash}\",\"sample_rate_hz\":{},\"full_scale_pa\":{},\"resonator_isolated_hz\":{resonator},\"supports\":false,\"mallet_shaft\":false,\"acceptance\":false,\"elapsed_s\":{}}}",
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
    /// Sound speed of the card's gas record [m/s].
    pub sound_speed_m_s: f64,
    /// Wall time of the surface and the modal transfers [s].
    pub radiation_wall_s: f64,
    /// Wall time of the strike render [s].
    pub render_wall_s: f64,
}

/// Radiate every mode of `s` to the card's receiver and render the cited mallet's stroke.
pub fn strike(card: &Card, s: &Structure) -> Result<Strike, String> {
    let gas_text = std::fs::read_to_string(&card.gas_path)
        .map_err(|e| format!("{}: {e}", card.gas_path.display()))?;
    let gas = read_gas_state(&mut gas_text.split_ascii_whitespace(), None)
        .map_err(|e| format!("gas record: {e:?}"))?;
    let medium = Medium::from(&gas);
    let read = |p: &Path| std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()));

    let started = Instant::now();
    let (surface, shapes) = radiation::surface_and_shapes(card, &s.mesh, &s.modes)?;
    let transfers = radiation::modal_transfers(card, medium, &surface, &shapes, &s.modes)?;
    let radiation_wall_s = started.elapsed().as_secs_f64();

    let loss = card
        .loss_path
        .as_ref()
        .map(|p| read(p).and_then(|text| render::Loss::read(&text)))
        .transpose()?;
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
    let modes = render::bar_modes(&s.modes, loss.as_ref(), &transfers.pressure, s.rigid.len());

    let started = Instant::now();
    let mut resonator_isolated_hz = None;
    let mut body: Box<dyn render::Body> = match &card.resonator {
        None => Box::new(render::DirectBody::new(card, modes)?),
        Some(res) => {
            resonator_isolated_hz = Some(resonator::isolated_hz(&resonator::cavity(res, &gas)?));
            let shapes = resonator::wet_volume_shapes(&s.mesh, &s.modes, res.wet_x_m);
            Box::new(resonator::CavityBody::new(card, res, &gas, modes, shapes)?)
        }
    };
    let stroke = render::render(card, body.as_mut(), &projection, &mallet)?;
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
