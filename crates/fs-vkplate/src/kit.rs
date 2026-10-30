//! A standard drum kit built from one solver: every head is a tensioned Mylar film (membrane
//! tension plus the film's own bending stiffness) and every cymbal a free-edged von Karman shallow
//! shell. All numbers below are declared ESTIMATES of ordinary instruments (dimensions from common
//! sizes, film and bronze from handbook values), not measurements of a particular specimen.

use crate::*;
use crate::coil::{CANOPUS_VINTAGE, CoilGeometry, CoilWire};
use crate::materials::KitMaterials;
use fs_material::state_point::IsotropicElasticStatePoint;

pub const OUT_RATE: f64 = 48_000.0;

/// Effective tip mass of a hickory stick in a short impact: FrankenSim's own Z5A stick estimate
/// (crates/fs-couple/examples/percussion: contact_effective_mass_kg = 0.00724314).
pub const STICK_TIP_KG: f64 = 0.00724314;

/// Default hi-hat stick-tip radius used by native and live impacts [m].
/// A declared KitMaterials striker keeps its own sourced impact law.
pub const HIHAT_STICK_TIP_RADIUS_M: f64 = 0.004;

/// Polyester (Mylar) film.
const MYLAR: (f64, f64, f64) = (4.0e9, 0.38, 1390.0);
/// B20 cymbal bronze (E Pa, nu, rho kg/m^3): handbook values, not a measured specimen.
pub const BRONZE: (f64, f64, f64) = (1.05e11, 0.33, 8600.0);

pub struct Piece {
    pub instrument: Instrument,
    pub oversample: usize,
    /// Seconds of simulation to discard at the start (snare-wire preload settling, hi-hat closing).
    pub settle: f64,
    pub notes: Vec<String>,
}

#[derive(Clone, Copy)]
pub struct Head {
    pub thickness: f64,
    pub tension: f64,
    pub t60: [(f64, f64); 2],
}

/// Snare wires: coiled steel strands under the resonant head, pressed up into it by the strainer.
/// Strand count, spacing, length, tension, preload and contact are declared, not sourced. Without
/// a sourced wire (`KitMaterials::with_snare_wire`) each strand is a plain string of the fitted
/// `mass_per_m` and `t60` (FITTED_WIRES_RECEIPT); with one it is the coil `coil` of that steel
/// (`crate::coil`), and `mass_per_m` and `t60` are unused.
#[derive(Clone, Copy, Debug)]
pub struct Wires {
    pub strands: i32,
    pub spacing: f64,
    pub length: f64,
    pub mass_per_m: f64,
    pub tension: f64,
    pub preload_per_m: f64,
    pub t60: f64,
    /// The wire-head contact law; a coil of a bound steel keeps only its loss and computes the
    /// stiffness (`coil::turn_contact`).
    pub contact: Law,
    pub coil: CoilGeometry,
}

impl Wires {
    pub fn estimated() -> Wires {
        // a loosely set strainer: slack strands slap the head instead of tracking it (measured: the most buzz)
        Wires { strands: 18, spacing: 0.004, length: 0.36, mass_per_m: 0.004, tension: 0.8, preload_per_m: 0.22, t60: 0.35,
                contact: Law::Power { k: 3.0e7, alpha: 1.5, mu: 0.35 }, coil: CANOPUS_VINTAGE }
    }
}

/// Carried by every drum whose snare wires are the fitted strings.
pub const FITTED_WIRES_RECEIPT: &str = "snare wires: fitted to a recording (96ccf8f1b), not physics: plain strings with a declared mass per length and one T60, no coil bending and no material or air loss law";

#[derive(Clone, Copy)]
pub struct DrumSpec {
    pub name: &'static str,
    pub radius: f64,
    pub depth: f64,
    pub batter: Head,
    pub reso: Head,
    pub grid: f64,
    pub oversample: usize,
    pub strike_xy: (f64, f64),
    /// (mass kg, law, footprint radius m)
    pub striker: (f64, Law, f64),
    pub muffle: Option<(f64, f64, f64, f64)>,
    pub wires: Option<Wires>,
    /// Front-head port radius (kick drums), m.
    pub port: Option<f64>,
    pub mic: [f64; 3],
    /// The floor, in the drum's frame (the batter head at the origin; None: free field).
    pub floor: Option<Floor>,
    /// Where the tensions, strike point, stick and (kick) beater, pillow and microphone come
    /// from. It joins the notes of every drum built from the spec.
    pub provenance: &'static str,
}

/// The presets' provenance (ruling 7486). Their head tensions, strike points and stick loss, and
/// the kick's felt crown and pillow, were adopted in 96ccf8f1b because they graded closer to the
/// Salamander overheads. `card_drum` replaces them. The kick's microphone, also adopted there, is
/// the declared overhead since ruling 8318.
pub const FITTED_DRUM_RECEIPT: &str = "drum: fitted to a recording (96ccf8f1b), not physics: head tensions, strike point and stick Hunt-Crossley loss 0.05 s/m; the kick's felt crown and pillow";

/// Carried by a drum whose microphone, and its image in the floor, both lie on the batter's side
/// of the resonant head (the kick's overhead, ruling 8318): the baffled microphone never hears a
/// surface facing away from it.
pub const OVERHEAD_HEARS_BATTER_ONLY_RECEIPT: &str = "overhead hears the batter only: the resonant head faces away and the baffled (Rayleigh) microphone has no path around the shell; its radiation, and the together pair's dipole cancellation, are not modelled";

/// How far the kick's hoops stand beyond its shell radius [m], a declared estimate (ruling 8318):
/// the drum lies on them, so the floor is tangent to their lowest point.
pub const KICK_HOOP_RISE_M: f64 = 0.02;

/// The presets' stick Hunt-Crossley loss [s/m], adopted by grade in 96ccf8f1b (it was 0.2).
const FITTED_STICK_LOSS_S_M: f64 = 0.05;

/// The stick's Hunt-Crossley loss [s/m] on the card route: 88cccbb78's declared estimate. No
/// source gives the loss of a hickory tip on PET film or bronze; the 3-6 ms contact Wagner (2006)
/// measured is a falsifier, not an input.
pub const DECLARED_STICK_LOSS_S_M: f64 = 0.2;

/// The presets' hickory tip (E Pa, nu): a handbook value, not the hickory card.
const PRESET_HICKORY: (f64, f64) = (12.0e9, 0.3);

pub(crate) fn hertz_stick(tip_radius: f64, e_target: f64, nu_target: f64) -> Law {
    // hickory tip on the target: K = 4/3 E* sqrt(R)
    hertz_pair(tip_radius, PRESET_HICKORY, (e_target, nu_target), FITTED_STICK_LOSS_S_M)
}

/// The existing two-elastic-half-space contact, shared by presets and sources.
pub(crate) fn hertz_pair(tip_radius: f64, (e_w, nu_w): (f64, f64), (e_target, nu_target): (f64, f64), loss_s_m: f64) -> Law {
    let es = 1.0 / ((1.0 - nu_w * nu_w) / e_w + (1.0 - nu_target * nu_target) / e_target);
    Law::Power { k: 4.0 / 3.0 * es * tip_radius.sqrt(), alpha: 1.5, mu: loss_s_m }
}

fn elastic_values(material: &IsotropicElasticStatePoint) -> (f64, f64, f64) {
    (material.young_modulus_pa(), material.poisson_ratio(), material.density_kg_m3())
}

fn instrument(rate: f64, gas: GasState, materials: Option<&KitMaterials>) -> Result<Instrument, String> {
    match materials {
        Some(materials) => Instrument::with_materials(rate, gas, materials.clone()),
        None => Instrument::with_gas(rate, gas),
    }
}

/// Where a drum head's own loss comes from: the preset's recording-fitted T60
/// (no imported materials), a sourced bending viscosity, or a stated unknown
/// entered as zero.
#[derive(Clone, Copy)]
enum FilmLoss<'a> {
    FittedT60,
    Sourced(&'a PlaneStressViscosity, f64),
    /// No surface loss on the plate: the film loss is a stated unknown, or
    /// the cards' film loss curve, which enters through the head radiation.
    NoSurfaceLoss,
}

fn head_body(name: &str, r: f64, grid: f64, hd: Head, muffle: Option<(f64, f64, f64, f64)>, elastic: (f64, f64, f64), loss: FilmLoss) -> Result<Body, String> {
    let t = hd.thickness;
    let th = move |_r: f64| t;
    let z = |_r: f64| 0.0;
    let spec = PlateSpec {
        radius: r, hole: 0.0, grid, edge: Edge::Clamped, thickness: &th, rest_height: &z,
        young: elastic.0, poisson: elastic.1, density: elastic.2, tension: hd.tension, bending: true,
        nonlinear: Nonlinear::Berger, t60: hd.t60, support: None, muffle,
    };
    match loss {
        FilmLoss::FittedT60 => Body::plate(name, &spec),
        FilmLoss::Sourced(law, angle) => Body::plate_with_viscosity(name, &spec, law, angle),
        FilmLoss::NoSurfaceLoss => Body::plate_without_surface_loss(name, &spec),
    }
}

/// Lowest (0,1) membrane frequency of a head, ignoring air (for the radiation resistance).
fn f01(r: f64, hd: Head, density: f64) -> f64 {
    2.405 * (hd.tension / (density * hd.thickness)).sqrt() / (2.0 * PI * r)
}

pub fn drum(s: &DrumSpec, hits: &[(f64, f64)]) -> Result<Piece, String> {
    drum_with_gas(s, hits, standard_gas())
}

/// The same structural data in one material-resolved cavity/exterior gas.
pub fn drum_with_gas(s: &DrumSpec, hits: &[(f64, f64)], gas: GasState) -> Result<Piece, String> {
    build_drum(s, hits, gas, None, false)
}

/// Fixed geometry and installed tension, with both head materials and a
/// spherical striker resolved from source data. The head loss is the selected
/// material bending loss, the head cards' film loss curve (its receipt and the
/// edge loss's stated unknown join the notes), or a declared stated unknown
/// entered as zero (its receipt joins the notes); with none the build refuses with
/// `materials::FILM_LOSS_REFUSAL` instead of keeping the preset's fitted T60.
/// Support/wire laws stay as declared in `s`. The supplied impactor replaces
/// `s.striker` completely.
pub fn drum_with_materials(s: &DrumSpec, hits: &[(f64, f64)], gas: GasState, materials: &KitMaterials) -> Result<Piece, String> {
    build_drum(s, hits, gas, Some(materials), false)
}

/// The kick's receipt on the card route: its beater is the spec's (the card drum's is 88cccbb78's
/// declared felt; see its provenance).
pub const KICK_BEATER_PRESET_RECEIPT: &str = "kick beater: the spec's felt beater (declared estimate, not a material card)";

/// `drum_with_materials` for the kick: heads from the cards, while the drum
/// keeps its declared `s.striker`, a felt beater that no material card
/// describes yet. KICK_BEATER_PRESET_RECEIPT replaces the card striker's tip
/// note, and the card striker is unused.
pub fn kick_with_materials(s: &DrumSpec, hits: &[(f64, f64)], gas: GasState, materials: &KitMaterials) -> Result<Piece, String> {
    let mut piece = build_drum(s, hits, gas, Some(materials), true)?;
    piece.notes.push(kick_pillow_note(s));
    Ok(piece)
}

/// How the kick's receipt on the card route begins: the preset names a pillow,
/// and nothing in the render models one.
pub const KICK_PILLOW_RECEIPT: &str = "kick: pillow not modelled";

/// KICK_PILLOW_RECEIPT with what stands in for the pillow: the spec's declared
/// muffle dashpot, if any, which is an estimate and not a pillow model.
pub fn kick_pillow_note(s: &DrumSpec) -> String {
    let stand_in = s.muffle.map_or_else(
        || "no muffle dashpot is declared".to_string(),
        |(_, _, radius, c)| format!("only the preset's declared muffle dashpot ({c} N s/m^3 over a {:.0} mm radius on each head, an estimate) stands in for it", radius * 1e3),
    );
    format!("{KICK_PILLOW_RECEIPT}; {stand_in}; a porous absorber in the heads' near field, from a measured flow resistivity, is a separate unit")
}

fn build_drum(s: &DrumSpec, hits: &[(f64, f64)], gas: GasState, materials: Option<&KitMaterials>, preset_striker: bool) -> Result<Piece, String> {
    let (batter_loss, reso_loss) = match materials {
        None => (FilmLoss::FittedT60, FilmLoss::FittedT60),
        Some(m) => match (m.primary_bending_loss(), m.secondary_bending_loss()) {
            (Some((a, angle_a)), Some((b, angle_b))) => (FilmLoss::Sourced(a, angle_a), FilmLoss::Sourced(b, angle_b)),
            _ if m.film_loss_is_stated_unknown() || m.has_film_loss_curve() => (FilmLoss::NoSurfaceLoss, FilmLoss::NoSurfaceLoss),
            _ => return Err(crate::materials::FILM_LOSS_REFUSAL.into()),
        },
    };
    let rate = OUT_RATE * s.oversample as f64;
    let mut ins = instrument(rate, gas, materials)?;
    let primary = materials.map_or(MYLAR, |m| elastic_values(m.primary()));
    let secondary = materials.map_or(MYLAR, |m| elastic_values(m.secondary()));
    let mut notes = vec![s.provenance.to_string()];
    match materials.and_then(KitMaterials::film_loss_note) {
        Some(film) => notes.extend([film, crate::materials::EDGE_LOSS_STATED_UNKNOWN_RECEIPT.into()]),
        None if matches!(batter_loss, FilmLoss::NoSurfaceLoss) => notes.push(crate::materials::FILM_LOSS_STATED_UNKNOWN_RECEIPT.into()),
        None => {}
    }
    if let Some(m) = materials {
        notes.extend(m.head_reduction_notes());
        notes.push(if preset_striker { KICK_BEATER_PRESET_RECEIPT.into() } else { m.striker_tip_note() });
    }
    let mut batter = head_body(&format!("{} batter", s.name), s.radius, s.grid, s.batter, s.muffle, primary, batter_loss)?;
    batter.origin = [0.0, 0.0, 0.0];
    batter.facing = 1.0;
    let mut reso = head_body(&format!("{} resonant", s.name), s.radius, s.grid, s.reso, s.muffle, secondary, reso_loss)?;
    reso.origin = [0.0, 0.0, -s.depth];
    reso.facing = -1.0;
    notes.push(format!("heads: {} nodes each at {:.1} mm; (0,1) in vacuo {:.0} / {:.0} Hz", batter.n, s.grid * 1e3, f01(s.radius, s.batter, primary.2), f01(s.radius, s.reso, secondary.2)));
    // outside air on each head's volume motion: a baffled piston's radiation mass, 8 rho a^3 / 3
    batter.air_mass = 8.0 * gas.density * s.radius.powi(3) / 3.0;
    reso.air_mass = batter.air_mass;
    // Without the fitted T60 (which lumped it in), the non-volume shapes take
    // their exterior resistance and air mass from the baked table, in this
    // gas, whatever the film-loss source.
    if !matches!(batter_loss, FilmLoss::FittedT60) {
        let table = crate::radiation::RadiationTable::shipped()?;
        // The volume shapes stay on the compact monopole below; the receipt
        // states its error against the baked baffled self-resistance.
        let compact = RadDamp::monopole(vec![], f01(s.radius, s.batter, primary.2), &gas).resistance;
        let film = materials.filter(|m| m.has_film_loss_curve());
        for (index, head, label, hd, elastic) in [
            (0, &mut batter, "batter", s.batter, primary),
            (1, &mut reso, "resonant", s.reso, secondary),
        ] {
            let coincidence_hz = gas.sound_speed.powi(2) / (2.0 * PI * hd.thickness)
                * (12.0 * elastic.2 * (1.0 - elastic.1 * elastic.1) / elastic.0).sqrt();
            // The card's film loss is the head's own generalized-Maxwell law on its bending
            // operator (every grid mode); the radiation columns keep air mass, radiation and
            // the viscous layer on every shape they retain.
            let mut lines = match film {
                Some(m) => {
                    let law = crate::film::fit(
                        &|hz| m.film_tan_delta(index, hz),
                        m.film_band_hz(index)?,
                        crate::film::matsuzawa_anchor_hz(),
                    )?;
                    head.attach_film_maxwell(law)?
                }
                None => vec![],
            };
            lines.extend(head.attach_head_radiation_with_viscous(
                &table,
                s.radius,
                &gas,
                Some(compact),
                coincidence_hz,
            )?);
            for line in lines {
                notes.push(format!("{label} {line}"));
            }
        }
    }
    let bi = ins.add(batter);
    let ri = ins.add(reso);
    let port = s.port.map(|pr| {
        let mut q = Port::new(pr, s.reso.thickness, [0.0, -0.6 * s.radius, -s.depth], -1.0);
        // small-ka radiation resistance at the Helmholtz frequency of cavity plus port
        let v = PI * s.radius * s.radius * s.depth;
        let fh = gas.sound_speed / (2.0 * PI) * (q.area / (v * q.length)).sqrt();
        q.resistance = gas.density * (2.0 * PI * fh).powi(2) / (2.0 * PI * gas.sound_speed);
        notes.push(format!("port {:.0} mm: Helmholtz {fh:.0} Hz with the cavity alone", 2e3 * pr));
        q
    });
    ins.cavities.push(Cavity { heads: vec![bi, ri], volume: PI * s.radius * s.radius * s.depth, port });
    ins.rads.push(RadDamp::monopole(vec![bi, ri], f01(s.radius, s.batter, primary.2), &gas));
    let mut settle = 0.0;
    if let Some(wr) = s.wires {
        let (mass_per_m, tension, preload_per_m, len) = (wr.mass_per_m, wr.tension, wr.preload_per_m, wr.length);
        // a coil of a bound steel presses the film with its computed Hertz law (the declared loss
        // kept); fitted strings keep the spec's law
        let (coil, contact) = match materials.and_then(|m| m.snare_wire().map(|steel| (m, steel))) {
            Some((m, steel)) => {
                let c = CoilWire::new(&wr.coil, steel, &gas)?;
                notes.extend(c.notes(tension)?);
                let Law::Power { mu, .. } = wr.contact;
                let film = (m.secondary().young_modulus_pa(), m.secondary().poisson_ratio());
                let (law, note) = crate::coil::turn_contact(&wr.coil, steel, film, s.reso.thickness, s.grid, preload_per_m, mu)?;
                notes.push(note);
                (Some(c), law)
            }
            None => {
                notes.push(FITTED_WIRES_RECEIPT.into());
                (None, wr.contact)
            }
        };
        let h = s.grid;
        let n = (len / h).round() as usize;
        let x0 = -(n as f64 - 1.0) * h / 2.0;
        let mut count = 0;
        let mut rows: Vec<i64> = vec![];
        for k in 0..wr.strands {
            let y = (k as f64 - (wr.strands as f64 - 1.0) / 2.0) * wr.spacing;
            let row = (y / h).round() as i64;
            if rows.contains(&row) {
                continue; // two strands on one grid row would share contact nodes
            }
            rows.push(row);
            let y = row as f64 * h;
            let name = format!("wire {row}");
            let mut wb = match &coil {
                Some(c) => Body::coil(&name, n, h, c.mass_per_m(), tension, c.bending_n_m2(), &|omega, share| c.decay_rate(omega, share), x0, y)?,
                None => Body::string(&name, n, h, mass_per_m, tension, wr.t60, x0, y),
            };
            wb.origin = [0.0, 0.0, -s.depth];
            for k in 0..n {
                wb.bias[k] = -preload_per_m * h; // toward the head (outward is +, the head is above)
            }
            let wi = ins.add(wb);
            for k in 0..n {
                let x = x0 + k as f64 * h;
                if let Some(hn) = ins.bodies[ri].node_at(x, y) {
                    let mut c = Contact::new(vec![(ri, hn, 1.0)], vec![(wi, k, 1.0)], 0.0, contact);
                    c.fixed_point_exit = true;
                    ins.contacts.push(c);
                }
            }
            count += 1;
        }
        notes.push(format!("snare: {count} wires of {n} nodes, {:.0} contact points", ins.contacts.len()));
        settle = 0.25;
    }
    let (m, law, foot) = match materials {
        Some(m) if !preset_striker => (m.striker_mass_kg(), m.impact_law(), 0.0),
        _ => s.striker,
    };
    let st = ins.add(Body::point("striker", m));
    let fp = ins.bodies[bi].footprint(s.strike_xy.0, s.strike_xy.1, foot);
    notes.push(format!("striker {m} kg over {} nodes", fp.len()));
    let mut c = Contact::new(fp.into_iter().map(|(k, w)| (bi, k, w)).collect(), vec![(st, 0, 1.0)], 0.0, law);
    c.active = false;
    ins.contacts.push(c);
    let ci = ins.contacts.len() - 1;
    for &(t, v) in hits {
        ins.strikes.push(Strike { time: t + settle, striker: st, contact: ci, speed: v, clear: 2e-4, done: false });
    }
    let mut mic = Mic::new(s.mic);
    mic.floor = s.floor;
    ins.mics.push(mic);
    // the resonant head, at z = -depth facing -z, reaches the microphone only from beyond it
    let beyond_reso = |p: [f64; 3]| p[2] < -s.depth;
    if !beyond_reso(s.mic) && !s.floor.is_some_and(|f| beyond_reso(f.image(s.mic))) {
        notes.push(OVERHEAD_HEARS_BATTER_ONLY_RECEIPT.into());
    }
    let oversample = if !matches!(batter_loss, FilmLoss::FittedT60) {
        pick_oversample(&mut ins, 0.6)?
    } else { s.oversample };
    ins.validate_bending_loss_clock()?;
    ins.validate_head_radiation()?;
    Ok(Piece { instrument: ins, oversample, settle, notes })
}

// Sizes and heads follow the reference kit we grade against (Salamander Drumkit, A. Holm, CC BY-SA
// 3.0): 18x18 kick with a 2-ply batter and a coated single-ply front head, 14x6 snare, 12x7 rack tom,
// 14x14 floor tom. Tensions are the drummer's tuning: chosen so each fundamental lands near the
// recording's. Loss points come from the recordings' own band decay times (see README).

const PLY1: f64 = 2.54e-4; // single-ply 10 mil
const PLY2: f64 = 3.56e-4; // two plies of 7 mil
const SNARE_SIDE: f64 = 7.62e-5; // 3 mil
const RESO: f64 = 1.9e-4; // 7.5 mil clear resonant head

pub fn snare() -> DrumSpec {
    DrumSpec {
        name: "snare 14x6",
        radius: 0.1703,
        depth: 0.1524,
        // recording: ~0.55 s T60 across the mid bands, with off-centre strike exciting (1,1) at ~235 Hz and wire rattle
        batter: Head { thickness: PLY1, tension: 1450.0, t60: [(150.0, 0.50), (3000.0, 0.45)] },
        reso: Head { thickness: SNARE_SIDE, tension: 2250.0, t60: [(150.0, 0.50), (3000.0, 0.45)] },
        grid: 3.5e-3,
        oversample: 2,
        strike_xy: (0.072, 0.025),
        striker: (0.0060, hertz_stick(0.005, MYLAR.0, MYLAR.1), 0.0),
        muffle: None,
        wires: Some(Wires::estimated()),
        port: None,
        mic: [0.35, 0.25, 0.9], // an overhead about 1 m away, 25 degrees off the drum's axis
        floor: Some(floor_z(-0.65)),
        provenance: FITTED_DRUM_RECEIPT,
    }
}

pub fn snare_wires_off() -> DrumSpec {
    DrumSpec { wires: None, name: "snare 14x6, wires off", ..snare() }
}

pub fn kick() -> DrumSpec {
    DrumSpec {
        name: "kick 18x18",
        radius: 0.221,
        depth: 0.4572,
        // a 2-ply batter and a pillow: fundamental on the recording's 44.5 Hz, the beater crown's
        // click decaying within ~60-100 ms
        batter: Head { thickness: PLY2, tension: 358.0, t60: [(45.0, 0.48), (4000.0, 0.18)] },
        reso: Head { thickness: PLY1, tension: 253.0, t60: [(45.0, 0.48), (4000.0, 0.18)] },
        grid: 4.5e-3,
        oversample: 1,
        strike_xy: (0.055, 0.015),
        // convex felt beater crown: point contact at initial impact launches the 0.5-4 kHz beater click
        striker: (0.032, Law::Power { k: 5.5e9, alpha: 2.25, mu: 0.04 }, 0.0),
        // a pillow resting against the lower part of the heads
        muffle: Some((0.0, -0.16, 0.07, 22.0)),
        wires: None,
        port: None,
        // the drum lies on its side, y up: a declared overhead 1 m above the shell's top edge,
        // 0.1 m toward the drummer so that the batter (z = 0, facing +z) passes the baffled
        // microphone's facing test, and the floor under the hoops (ruling 8318)
        mic: [0.0, 0.221 + 1.0, 0.1],
        floor: Some(floor_y(-(0.221 + KICK_HOOP_RISE_M))),
        provenance: FITTED_DRUM_RECEIPT,
    }
}

pub fn tom(size_in: f64) -> DrumSpec {
    let (radius, depth, batter, reso, strike_xy, grid, name) = if size_in <= 12.5 {
        (0.146, 0.1778,
         Head { thickness: PLY1, tension: 1300.0, t60: [(139.0, 1.55), (2500.0, 0.25)] },
         Head { thickness: RESO, tension: 975.0, t60: [(139.0, 1.55), (2500.0, 0.25)] },
         (0.275 * 0.146, 0.0),
         4.0e-3,
         "rack tom 12x7")
    } else {
        (0.1703, 0.3556,
         Head { thickness: PLY2, tension: 410.0, t60: [(69.0, 1.5), (3500.0, 0.28)] },
         Head { thickness: RESO, tension: 215.0, t60: [(69.0, 1.5), (3500.0, 0.28)] },
         (0.43 * 0.1703, 0.0),
         4.2e-3,
         "floor tom 14x14")
    };
    DrumSpec {
        name,
        radius,
        depth,
        batter,
        reso,
        grid,
        oversample: 1,
        strike_xy,
        striker: (STICK_TIP_KG, hertz_stick(0.005, MYLAR.0, MYLAR.1), 0.0),
        muffle: None,
        wires: None,
        port: None,
        mic: [0.35, 0.25, 0.9],
        floor: Some(floor_z(if size_in <= 12.5 { -0.75 } else { -0.55 })),
        provenance: FITTED_DRUM_RECEIPT,
    }
}

/// Scale both heads' tension (a drummer's tuning turn), keeping everything else.
pub fn retune(mut s: DrumSpec, factor: f64) -> DrumSpec {
    s.batter.tension *= factor;
    s.reso.tension *= factor;
    s
}

/// One card drum's tuning and playing position (ruling 7486).
struct CardDrum {
    name: &'static str,
    /// Batter and resonant tension [N/m]: the batter tuned so the card drum's own lowest mode sits
    /// at the late pitch of the Salamander strokes (rule 1, below), the resonant head at
    /// 88cccbb78's ratio to it.
    tensions: (f64, f64),
    /// 88cccbb78's declared playing position [m] (rule 2).
    strike_xy: (f64, f64),
    provenance: &'static str,
}

// Rule 1 of ruling 7486: a drummer tunes the heads, so the tension may put the drum's own computed
// lowest mode at the pitch measured on the reference stroke. The pitch is the lowest spectral peak
// within 20 dB of the strongest, 30-400 Hz (Hann, 2^19 points), in the small-amplitude tail of the
// Salamander overhead strokes the kit is compared with, where two dynamics agree (mail 7506):
// kick 40.6 Hz (kick_OH_F_1, _P_3, 0.3-0.6 s), snare 134.4 Hz (snare2OFF_OH_F_1, _P_2,
// 0.6-1.0 s), rack tom 133.9 Hz (hiTom_OH_F_4, _P_3, 1.0-1.6 s), floor tom 60.4 Hz
// (loTom_OH_FF_5, _MP_3, 1.0-1.6 s). The card drum's mode is read the same way from a 0.1 m/s
// stroke (the Berger stretch negligible), 0.3-1.0 s after onset; both heads keep 88cccbb78's
// ratio. Found by scaling both preset heads together (f proportional to the square root of the
// tension) until that mode was within 0.3 Hz of the target, rounding the batter to four figures
// and setting the resonant head at 88cccbb78's ratio to it; the snare, 0.37 Hz high at that ratio,
// took one more step. These drums' own modes: snare 134.40, rack tom 133.85, floor tom 60.42
// and kick 40.56 Hz (trj, 2026-10-26). The early pitch of a hard stroke is up to 15 Hz sharper
// (the stretch); nothing grades balance or decay.
const CARD_DRUMS: [CardDrum; 4] = [
    CardDrum {
        name: "snare 14x6",
        tensions: (1072.0, 1072.0 * 2080.0 / 1330.0),
        strike_xy: (0.06, 0.02),
        provenance: "drum (card route, ruling 7486): tensions tuned to the recording's pitch, lowest mode 134.4 Hz (Salamander snare2OFF_OH_F_1/_P_2, late tail); strike point (60, 20) mm, a declared playing-position estimate (88cccbb78); stick Hunt-Crossley loss 0.2 s/m and the wires' tension 2.0 N, preload 0.25 N/m and contact loss 2.0 s/m are declared estimates (88cccbb78); the wires' contact stiffness is computed from the cards (Hertz, elliptical, helix; ruling 7520)",
    },
    CardDrum {
        name: "kick 18x18",
        tensions: (375.5, 375.5 * 253.0 / 358.0),
        strike_xy: (0.0, 0.02),
        provenance: "drum (card route, ruling 7486): tensions tuned to the recording's pitch, lowest mode 40.6 Hz (Salamander kick_OH_F_1/_P_3, late tail); strike point (0, 20) mm, a declared playing-position estimate (88cccbb78); beater 0.12 kg with a 3 cm felt face (k 1e8 N/m^2.5, mu 0.3 s/m) and the muffle dashpot 25 N s/m^3 standing in for the pillow are declared estimates (88cccbb78); the overhead microphone is a declared estimate 1 m above the shell's top edge, (0, 1.221, 0.1) m with y up, its 0.1 m toward the drummer there to clear the batter's facing test, and the floor is tangent to the hoops at a declared 0.02 m beyond the shell radius (ruling 8318); the nominal 18x18 cylinder overestimates the air volume, the shells being very thick and unlathed (Salamander README); unknowns: whether the front head is ported (no source says, and none is modelled; a trial 102 mm port raised the first together line from 46.9 to 54.9 Hz and added a together line at 21-25 Hz, which the linear lumped model does not predict, perhaps by the strike's amplitude, untested; ruling 8334), and the one-of-a-kind rope tuning, without lugs, so the heads' tension need not be uniform around the rim",
    },
    CardDrum {
        name: "rack tom 12x7",
        tensions: (1264.0, 1264.0 * 950.0 / 1270.0),
        strike_xy: (0.35 * 0.146, 0.0),
        provenance: "drum (card route, ruling 7486): tensions tuned to the recording's pitch, lowest mode 133.9 Hz (Salamander hiTom_OH_F_4/_P_3, late tail); strike point 0.35 R, a declared playing-position estimate (88cccbb78); stick Hunt-Crossley loss 0.2 s/m, a declared estimate (88cccbb78)",
    },
    CardDrum {
        name: "floor tom 14x14",
        tensions: (478.0, 478.0 * 175.0 / 330.0),
        strike_xy: (0.35 * 0.1703, 0.0),
        provenance: "drum (card route, ruling 7486): tensions tuned to the recording's pitch, lowest mode 60.4 Hz (Salamander loTom_OH_FF_5/_MP_3, late tail); strike point 0.35 R, a declared playing-position estimate (88cccbb78); stick Hunt-Crossley loss 0.2 s/m, a declared estimate (88cccbb78)",
    },
];

/// The card route's drum (ruling 7486): `preset` with every value 96ccf8f1b adopted by grade
/// replaced. The tensions are tuned to the recording's pitch (rule 1); the strike point, the
/// wires' tension, preload and contact loss, and the kick's beater and muffle are 88cccbb78's
/// declared estimates (rules 2 and 3). Geometry, heads, grids and the microphone are kept. The
/// stick's loss is `DECLARED_STICK_LOSS_S_M`, which the card route passes to the hickory card.
/// The wires' contact stiffness is computed from the steel and film cards once the steel is
/// bound (`coil::turn_contact`); 88cccbb78's 2e7 N/m^1.5 here only reaches strings without one.
/// The wires-off snare is the snare without wires. Refuses a spec it has no card values for.
pub fn card_drum(preset: &DrumSpec) -> Result<DrumSpec, String> {
    let card_name = if preset.name == "snare 14x6, wires off" { "snare 14x6" } else { preset.name };
    let card = CARD_DRUMS.iter().find(|c| c.name == card_name)
        .ok_or_else(|| format!("card_drum: no card values for {:?}", preset.name))?;
    let mut s = DrumSpec { strike_xy: card.strike_xy, provenance: card.provenance, ..*preset };
    s.batter.tension = card.tensions.0;
    s.reso.tension = card.tensions.1;
    s.wires = preset.wires.map(|w| Wires { tension: 2.0, preload_per_m: 0.25, contact: Law::Power { k: 2.0e7, alpha: 1.5, mu: 2.0 }, ..w });
    if card_name == "kick 18x18" {
        // a felt beater face of 3 cm radius, power-law felt stiffening with compression (about
        // 5 mm of felt compression at 180 N), and the muffle dashpot against the lower part of the
        // heads, standing in for the pillow (not modelled; KICK_PILLOW_RECEIPT)
        s.striker = (0.12, Law::Power { k: 1.0e8, alpha: 2.5, mu: 0.3 }, 0.03);
        s.muffle = preset.muffle.map(|(x, y, r, _)| (x, y, r, 25.0));
    }
    Ok(s)
}

pub struct CymbalSpec {
    pub name: &'static str,
    pub radius: f64,
    pub bow_thickness: f64,
    pub bell_thickness: f64,
    pub bell_radius: f64,
    pub bow_rise: f64,
    pub bell_rise: f64,
    pub grid: f64,
    pub t60: [(f64, f64); 2],
    /// Where the shell's dimensions or loss points come from, when a receipt is owed; it joins the
    /// notes of every cymbal built from the spec.
    pub provenance: Option<&'static str>,
}

/// The stick and rim constants a cymbal or hi-hat is played with (ruling 7486).
#[derive(Clone, Copy, Debug)]
pub struct CymbalInterface {
    /// The stick's Hunt-Crossley loss on bronze without a stick card [s/m].
    pub stick_loss_s_m: f64,
    /// The hi-hat stick's tip radius without a stick card [m].
    pub hat_stick_tip_m: f64,
    /// Two-lobe waviness of the hi-hat rims: the rim gap grows by `w (1 - cos 2 theta)` [m].
    pub hat_rim_waviness_m: f64,
    /// Dry friction coefficient between the hi-hat's rims (`hihat_with_interface`).
    pub hat_rim_friction: f64,
    /// Joins a cymbal's notes when no stick card is bound.
    pub cymbal_receipt: &'static str,
    /// Joins every hi-hat's notes.
    pub hat_receipt: &'static str,
}

/// The presets' interface, each value adopted in 96ccf8f1b because it graded closer to the
/// Salamander overheads (88cccbb78 had a 0.2 s/m stick loss, a 5 mm tip and flat rims). Its rim
/// law is HIHAT_RIM_LAW.
pub const FITTED_CYMBAL_INTERFACE: CymbalInterface = CymbalInterface {
    stick_loss_s_m: FITTED_STICK_LOSS_S_M,
    hat_stick_tip_m: HIHAT_STICK_TIP_RADIUS_M,
    hat_rim_waviness_m: 7.5e-5,
    hat_rim_friction: RIM_FRICTION,
    cymbal_receipt: "cymbal stick: fitted to a recording (96ccf8f1b), not physics: Hunt-Crossley loss 0.05 s/m",
    hat_receipt: "hi-hat: fitted to a recording (96ccf8f1b), not physics: rim loss 0.08 s/m (HIHAT_RIM_LAW), a two-lobe rim waviness of up to 0.15 mm, a 4 mm stick tip and, without a stick card, stick loss 0.05 s/m; rim friction RIM_FRICTION 0.5, picked by grade (963d1e278)",
};

/// The card route's interface: 88cccbb78's declared estimates, and the measured BRONZE_RIM_FRICTION.
/// Its rim law is DECLARED_HIHAT_RIM_LAW. No source gives a hi-hat's rim waviness, so the rims
/// are flat.
pub const DECLARED_CYMBAL_INTERFACE: CymbalInterface = CymbalInterface {
    stick_loss_s_m: DECLARED_STICK_LOSS_S_M,
    hat_stick_tip_m: 0.005,
    hat_rim_waviness_m: 0.0,
    hat_rim_friction: BRONZE_RIM_FRICTION,
    cymbal_receipt: "cymbal stick (ruling 7486): Hunt-Crossley loss 0.2 s/m, a declared estimate (88cccbb78)",
    hat_receipt: "hi-hat (ruling 7486): flat rims, a 5 mm stick tip and stick loss 0.2 s/m, declared estimates (88cccbb78); no rim contact loss (ruling 7865): 88cccbb78's 3.0 s/m stood in for rim friction and no source gives a separate loss at closure, so rim dissipation is the friction alone; rim friction 0.6 (Buckley, Cu-10Sn self-mated; ruling 7520): dry in air, a polished hemispherical pin at 2.5 N and 0.05 m/s (NASA TN D-6912, 1972, fig. 2), not B20 rims; primary values span 0.2-0.6 across contact geometry (Morin 1834: 0.20, a loaded bronze sled on bronze)",
};

/// Dry kinetic friction of tin bronze on itself, for the card route's hi-hat rims (ruling 7520
/// G1): Cu-10 at.% Sn (17 wt% Sn; B20 is 20 wt%) sliding on itself in air, 0.6, read from the
/// bar of Buckley 1972 (NASA TN D-6912, fig. 2: a 0.48 cm hemispherical rider on a disk of the
/// same alloy, polished, 250 g, 300 cm/min, 25 C, averaged over an hour of adhesive wear). The
/// other primary measurement, Morin 1834 (Tableau LXXVII), is 0.20 for an unnamed bronze sled on
/// a bronze bed at about 1.2 m/s; neither is a lathed B20 rim at a hi-hat's load and slip.
pub const BRONZE_RIM_FRICTION: f64 = 0.6;

/// 88cccbb78's declared hi-hat rim stiffness (reference cell as HIHAT_RIM_LAW), with no
/// Hunt-Crossley loss (ruling 7865). 88cccbb78's 3.0 s/m stood in for rim friction, which
/// BRONZE_RIM_FRICTION now carries from data, so keeping it would count that loss twice. No
/// source gives a separate loss at closure, so the rims dissipate through the friction alone.
pub const DECLARED_HIHAT_RIM_LAW: Law = Law::Power { k: 2.0e9, alpha: 1.5, mu: 0.0 };

/// Grid points (spacing h) with r in [r0, r1]. A grid too coarse to put four of them in the
/// annulus (a real-time grid) widens r1 to the first ring that does; the caller then spreads the
/// annulus's total stiffness and damping over those nodes. Offline grids never widen.
fn ring_nodes(h: f64, r0: f64, r1: f64) -> (f64, usize) {
    let count = |r1: f64| {
        let m = (r1 / h).ceil() as i64 + 1;
        let mut n = 0;
        for j in -m..=m {
            for i in -m..=m {
                let r = h * ((i * i + j * j) as f64).sqrt();
                if r >= r0 && r <= r1 + 1e-12 {
                    n += 1;
                }
            }
        }
        n
    };
    let mut r = r1;
    loop {
        let n = count(r);
        if n >= 4 {
            return (r, n);
        }
        r += 0.25 * h;
    }
}

fn cymbal_body(name: &str, c: &CymbalSpec, elastic: (f64, f64, f64), viscosity: Option<(&PlaneStressViscosity, f64)>) -> Result<Body, String> {
    cymbal_shell(name, c, true, elastic, viscosity)
}

/// A cymbal's thickness at radius `x` outside the bell (the bow law of `cymbal_shell`).
fn bow_thickness(c: &CymbalSpec, x: f64) -> f64 {
    c.bow_thickness * (1.15 - 0.15 * x / c.radius)
}

/// The felt washers' support of a stand-mounted cymbal: (r0, r1, stiffness N/m per m^2,
/// damping N s/m per m^2) over the annulus round the 13 mm centre hole. Authored estimates
/// (~1e5 N/m and ~3 N s/m in all), not felt data.
pub fn felt_support(c: &CymbalSpec) -> (f64, f64, f64, f64) {
    // felt washers clamp the annulus round the 13 mm centre hole
    let (r1, n) = ring_nodes(c.grid, 0.0065, 0.02);
    let mount_scale = if c.radius < 0.15 { 0.25 } else { 1.0 };
    let (ks, cs) = if r1 > 0.02 {
        let per = PI * (0.02f64.powi(2) - 0.0065f64.powi(2)) / (n as f64 * c.grid * c.grid);
        (9.0e7 * per * mount_scale, 2.7e3 * per * mount_scale)
    } else {
        (9.0e7 * mount_scale, 2.7e3 * mount_scale)
    };
    (0.0065, r1, ks, cs)
}

/// A cymbal shell; `felt` puts it on a stand's felt washers (a hi-hat's top shell hangs from the
/// clutch on the pull rod instead, so all the pedal force reaches the rims).
fn cymbal_shell(name: &str, c: &CymbalSpec, felt: bool, elastic: (f64, f64, f64), viscosity: Option<(&PlaneStressViscosity, f64)>) -> Result<Body, String> {
    let (rb, tb, tw, r) = (c.bell_radius, c.bell_thickness, c.bow_thickness, c.radius);
    // thickness: bell, a blend over 2 cm, then the bow thinning slightly toward the rim
    let th = move |x: f64| {
        let bow = tw * (1.15 - 0.15 * x / r);
        if x <= rb { tb } else if x < rb + 0.02 { let u = (x - rb) / 0.02; tb + (bow - tb) * (3.0 * u * u - 2.0 * u * u * u) } else { bow }
    };
    let (ab, al) = (c.bow_rise, c.bell_rise);
    // rest height above the rim: a slightly convex bow plus a raised bell (cosine cap)
    let z = move |x: f64| {
        let u = (1.0 - x / r).max(0.0);
        let bow = ab * (0.7 * u + 0.3 * u * u);
        let bell = if x < rb { al * 0.5 * (1.0 + (PI * x / rb).cos()) } else { 0.0 };
        bow + bell
    };
    let (r0, r1, ks, cs) = felt_support(c);
    let spec = PlateSpec {
        radius: r, hole: 0.0065, grid: c.grid, edge: Edge::Free, thickness: &th, rest_height: &z,
        young: elastic.0, poisson: elastic.1, density: elastic.2, tension: 0.0, bending: true,
        nonlinear: Nonlinear::VonKarman, t60: c.t60,
        // felt washers: ~1e5 N/m and ~3 N s/m in all, spread over the clamped annulus
        support: if felt { Some((r0, r1, ks, cs)) } else { None }, muffle: None,
    };
    match viscosity {
        Some((law, angle)) => Body::plate_with_viscosity(name, &spec, law, angle),
        None => Body::plate(name, &spec),
    }
}

pub fn splash() -> CymbalSpec {
    CymbalSpec { name: "splash 8", radius: 0.1016, bow_thickness: 0.68e-3, bell_thickness: 1.5e-3, bell_radius: 0.039, bow_rise: 0.008, bell_rise: 0.014, grid: 4.5e-3, t60: [(500.0, 6.8), (10000.0, 3.1)], provenance: None }
}

pub fn crash() -> CymbalSpec {
    CymbalSpec { name: "crash 18", radius: 0.2286, bow_thickness: 0.95e-3, bell_thickness: 1.6e-3, bell_radius: 0.075, bow_rise: 0.025, bell_rise: 0.015, grid: 6.0e-3, t60: [(500.0, 10.0), (10000.0, 2.0)], provenance: None }
}

pub fn ride() -> CymbalSpec {
    CymbalSpec { name: "ride 20", radius: 0.254, bow_thickness: 1.65e-3, bell_thickness: 2.25e-3, bell_radius: 0.085, bow_rise: 0.022, bell_rise: 0.020, grid: 6.2e-3, t60: [(500.0, 14.0), (10000.0, 3.2)],
                 provenance: Some("ride: fitted to a recording (96ccf8f1b), not physics: bow 1.65 mm, bell 2.25 mm, 10 kHz T60 3.2 s") }
}

/// The card route's ride (ruling 7486): the preset with 88cccbb78's declared shell, bow 1.60 mm
/// and bell 2.20 mm, and its 10 kHz loss point, 3.0 s. The grid (6.2 mm) is numerics and stays.
pub fn card_ride() -> CymbalSpec {
    CymbalSpec { bow_thickness: 1.6e-3, bell_thickness: 2.2e-3, t60: [(500.0, 14.0), (10000.0, 3.0)],
                 provenance: Some("ride (ruling 7486): bow 1.60 mm, bell 2.20 mm and the 10 kHz T60 3.0 s, declared estimates (88cccbb78)"), ..ride() }
}

fn pick_oversample(ins: &mut Instrument, max_ratio: f64) -> Result<usize, String> {
    let (lam, who) = ins.admission_rate2();
    let structural = (lam - ins.gas_admission_rate2(ins.rate)
        - ins.viscosity_admission_rate2(ins.rate) - ins.link_admission_rate2(ins.rate)).max(0.0);
    // Finer spatial meshes need smaller time steps (plate bending scales as
    // h^-4). Extend the bounded clock search, retaining the same CFL margin
    // and every previously admitted rate; never coarsen the physical grid.
    for os in [1usize, 2, 3, 4, 6, 8, 12, 16, 24, 32] {
        let rate = OUT_RATE * os as f64;
        if rate < ins.rate { continue; }
        let k = 1.0 / rate;
        // Explicit radiation's contribution is 2 D / dt, so it must use
        // this candidate's clock rather than the instrument's original rate.
        let candidate = structural + ins.gas_admission_rate2(rate)
            + ins.viscosity_admission_rate2(rate) + ins.link_admission_rate2(rate);
        if k * k * candidate / 4.0 <= max_ratio {
            ins.rate = rate;
            ins.validate_bending_loss_clock()?;
            return Ok(os);
        }
    }
    Err(format!("{who}: leapfrog unstable even at 32x oversampling"))
}

#[test]
fn g0_clock_refinement_preserves_the_cfl_margin_and_bounded_refusal() {
    // Two equal string masses have eigenvalues T and 3T when h=rho=1.
    // Choose physical tensions whose exact highest frequency falls between
    // ladder rungs; the oracle does not reuse the power-iteration estimate.
    for (minimum, expected) in [(0.5_f64, Some(1)), (7.0, Some(8)),
        (10.0, Some(12)), (14.0, Some(16)), (20.0, Some(24)),
        (28.0, Some(32)), (40.0, None)] {
        let lambda = 4.0 * 0.6 * (OUT_RATE * minimum).powi(2);
        let mut ins = Instrument::new(OUT_RATE);
        ins.add(Body::string("clock test", 2, 1.0, 1.0, lambda / 3.0,
            f64::INFINITY, 0.0, 0.0));
        let result = pick_oversample(&mut ins, 0.6);
        if let Some(os) = expected {
            assert_eq!(result.unwrap(), os);
            assert_eq!(ins.rate, OUT_RATE * os as f64);
            assert!(lambda / (4.0 * ins.rate.powi(2)) <= 0.6);
        } else {
            assert!(result.unwrap_err().contains("32x oversampling"));
            assert_eq!(ins.rate, OUT_RATE);
        }
    }
}

#[test]
fn g0_link_clock_covers_stiff_springs_and_retimes_dashpots() {
    // Relative motion of two unit masses has stiffness 2K and damping 2C.
    // K dt^2 + 2 C dt < 2 is the analytical leapfrog condition; the kit
    // retains its stricter 0.6 margin, independent of the row-bound code.
    for (spring, damping, expected) in [(1e12, 0.0, 24), (0.0, 2.0 * OUT_RATE, 4)] {
        let mut ins = Instrument::new(OUT_RATE);
        ins.add(Body::point("a", 1.0));
        ins.add(Body::point("b", 1.0));
        ins.links.push(Link { a: (0, 0), b: Some((1, 0)), k: spring, c: damping });
        assert_eq!(pick_oversample(&mut ins, 0.6).unwrap(), expected);
        let ratio = spring / (2.0 * ins.rate.powi(2)) + damping / ins.rate;
        assert!(ratio <= 0.6);
    }
}

#[test]
fn g0_link_clock_refuses_beyond_existing_limit() {
    let mut ins = Instrument::new(OUT_RATE);
    ins.add(Body::point("a", 1.0));
    ins.links.push(Link { a: (0, 0), b: None, k: 1e15, c: 0.0 });
    assert!(pick_oversample(&mut ins, 0.6).unwrap_err().contains("32x oversampling"));
    assert_eq!(ins.rate, OUT_RATE);
}

/// A suspended cymbal struck once per hit: (time s, speed m/s, radius fraction, stick tip radius m).
pub fn cymbal(c: &CymbalSpec, hits: &[(f64, f64, f64, f64)]) -> Result<Piece, String> {
    cymbal_with_gas(c, hits, standard_gas())
}

/// The same plate and contact mechanics with material-resolved propagation.
pub fn cymbal_with_gas(c: &CymbalSpec, hits: &[(f64, f64, f64, f64)], gas: GasState) -> Result<Piece, String> {
    build_cymbal(c, hits, gas, None, &FITTED_CYMBAL_INTERFACE, None)
}

/// The same with the stick's loss from `iface` (the card route passes DECLARED_CYMBAL_INTERFACE).
pub fn cymbal_with_interface(c: &CymbalSpec, hits: &[(f64, f64, f64, f64)], iface: &CymbalInterface) -> Result<Piece, String> {
    build_cymbal(c, hits, standard_gas(), None, iface, None)
}

/// Source-bound shell stiffness/mass and spherical-striker mass/contact.
/// With materials supplied, hit tuple radii remain strike positions; the
/// material binding's tip radius declares the common physical impactor.
pub fn cymbal_with_materials(c: &CymbalSpec, hits: &[(f64, f64, f64, f64)], gas: GasState, materials: &KitMaterials) -> Result<Piece, String> {
    build_cymbal(c, hits, gas, Some(materials), &FITTED_CYMBAL_INTERFACE, None)
}

/// A declared striker: its moving mass [kg], its contact law and, if it crushes, its crush.
pub type Striker = (f64, Law, Option<Crush>);

/// A cymbal of preset bronze struck by a declared striker, the same for every hit (time s, speed
/// m/s, radius fraction). `notes`, the striker's receipts, join the piece's in place of a stick
/// receipt. The shaft strike (`shaft::ShaftStrike`) enters here.
pub fn cymbal_struck_by(c: &CymbalSpec, hits: &[(f64, f64, f64)], gas: GasState, striker: Striker, notes: &[String]) -> Result<Piece, String> {
    // the hit tuple's tip radius is read only by the presets' own stick
    let hits: Vec<_> = hits.iter().map(|&(t, v, rf)| (t, v, rf, f64::NAN)).collect();
    let mut piece = build_cymbal(c, &hits, gas, None, &DECLARED_CYMBAL_INTERFACE, Some(striker))?;
    piece.notes.extend(notes.iter().cloned());
    Ok(piece)
}

/// A cymbal's thickness at its free edge [m].
pub fn edge_thickness(c: &CymbalSpec) -> f64 {
    bow_thickness(c, c.radius)
}

fn build_cymbal(c: &CymbalSpec, hits: &[(f64, f64, f64, f64)], gas: GasState, materials: Option<&KitMaterials>, iface: &CymbalInterface, striker: Option<Striker>) -> Result<Piece, String> {
    let mut ins = instrument(OUT_RATE, gas, materials)?;
    let elastic = materials.map_or(BRONZE, |m| elastic_values(m.primary()));
    let mut body = cymbal_body(c.name, c, elastic, materials.and_then(KitMaterials::primary_bending_loss))?;
    body.origin = [0.0, 0.0, 0.0];
    let n = body.n;
    let ci = ins.add(body);
    let mass = striker.map_or_else(|| materials.map_or(STICK_TIP_KG, KitMaterials::striker_mass_kg), |(m, _, _)| m);
    if !(mass.is_finite() && mass > 0.0) {
        return Err(format!("cymbal striker mass must be positive, got {mass} kg"));
    }
    let st = ins.add(Body::point("stick", mass));
    let mut notes = vec![format!("{}: {n} nodes at {:.1} mm", c.name, c.grid * 1e3)];
    notes.extend(c.provenance.map(String::from));
    if materials.is_none() && striker.is_none() {
        notes.push(iface.cymbal_receipt.into());
    }
    for (k, &(t, v, rf, tip)) in hits.iter().enumerate() {
        let fp = ins.bodies[ci].footprint(rf * c.radius, 0.0, 0.0);
        // the stick comes down onto the upper surface: plate on side a
        let law = striker.map_or_else(
            || materials.map_or_else(|| hertz_pair(tip, PRESET_HICKORY, (elastic.0, elastic.1), iface.stick_loss_s_m), KitMaterials::impact_law),
            |(_, l, _)| l,
        );
        let mut con = Contact::new(fp.into_iter().map(|(n_, w)| (ci, n_, w)).collect(), vec![(st, 0, 1.0)], 0.0, law);
        if let Some((_, _, Some(crush))) = striker {
            con = con.crushing(crush)?;
        }
        con.active = false;
        ins.contacts.push(con);
        ins.strikes.push(Strike { time: t, striker: st, contact: ins.contacts.len() - 1, speed: v, clear: 2e-4, done: false });
        let _ = k;
    }
    let mut mic = Mic::new([0.3, 0.2, 1.0]);
    mic.unbaffled = Some(c.radius);
    ins.mics.push(mic);
    let os = pick_oversample(&mut ins, 0.6)?;
    if let Some((loss, angle)) = ins.bodies[ci].bending_viscosity() {
        notes.push(format!("internal rate {} kHz; sourced Kelvin-Voigt bending loss {}, orientation {angle} rad; independent supports retained",
            OUT_RATE * os as f64 / 1000.0, loss.source().identity()));
    } else {
        let (s0, s1) = ins.bodies[ci].loss();
        notes.push(format!("internal rate {} kHz; loss s0 {s0:.3}/s s1 {s1:.3e} m^2/s", OUT_RATE * os as f64 / 1000.0));
    }
    Ok(Piece { instrument: ins, oversample: os, settle: 0.0, notes })
}

pub fn hihat_spec() -> (CymbalSpec, CymbalSpec) {
    let top = CymbalSpec { name: "hi-hat top 14", radius: 0.1778, bow_thickness: 0.85e-3, bell_thickness: 1.4e-3, bell_radius: 0.055, bow_rise: 0.015, bell_rise: 0.01, grid: 5.5e-3, t60: [(500.0, 8.0), (10000.0, 1.5)], provenance: None };
    let bottom = CymbalSpec { name: "hi-hat bottom 14", bow_thickness: 1.1e-3, ..top };
    (top, bottom)
}

/// Reference-cell normal contact estimate for bronze hi-hat rims, shared by
/// offline and live construction. Grid assembly scales K with cell area;
/// Hunt-Crossley loss mu [s/m] and the Hertz exponent remain unchanged.
pub const HIHAT_RIM_LAW: Law = Law::Power { k: 2.0e9, alpha: 1.5, mu: 0.08 };

/// Hi-hat: opposing shells on one rod. `gap` is the undeformed face-to-face rim opening [m].
/// Bow clearance follows the supplied profiles; pedal loading deforms them from this geometry.
pub fn hihat(gap: f64, pedal_n: f64, hits: &[(f64, f64)]) -> Result<Piece, String> {
    hihat_scaled(gap, pedal_n, hits, 1.0)
}

/// The same physical hi-hat and rim loss on an explicitly coarsened grid.
pub fn hihat_scaled(gap: f64, pedal_n: f64, hits: &[(f64, f64)], scale: f64) -> Result<Piece, String> {
    // bronze rims pressed together: stiff and lightly damped so rim chatter radiates 4-16 kHz sizzle
    hihat_grid(gap, pedal_n, hits, HIHAT_RIM_LAW, 0.8, scale)
}

/// The same with an explicit rim contact law and the inner edge of the contact annulus (fraction of the radius).
pub fn hihat_with(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64) -> Result<Piece, String> {
    hihat_grid(gap, pedal_n, hits, rim, from, 1.0)
}

/// Dry friction coefficient between the hi-hat's bronze rims. Unlubricated copper alloys on themselves
/// span about 0.3-1; with 0.3 and 0.5 the closed hat's octave-band decay grades inside the real-vs-real
/// ruler (0.40 and 0.32 octaves from the Salamander stroke, two real strokes 0.47; without friction 0.95),
/// and 0.5 was the closer of the two. Those grades predate opposing geometry and
/// finite washer rocking; they do not qualify the current preset. Picked by grade: the presets
/// (and so live, the browser and GPU) keep it and say so; the card route takes
/// BRONZE_RIM_FRICTION (ruling 7520 G1).
pub const RIM_FRICTION: f64 = 0.5;

/// The same with grid spacing multiplied by `scale`.
/// The reference-cell rim force law scales with area to preserve traction under refinement.
pub fn hihat_grid(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64, scale: f64) -> Result<Piece, String> {
    hihat_grid_mu(gap, pedal_n, hits, rim, from, scale, RIM_FRICTION)
}

/// The same with an explicit rim friction coefficient (0: the rims press but do not rub).
pub fn hihat_grid_mu(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64, scale: f64, friction: f64) -> Result<Piece, String> {
    hihat_grid_mu_with_gas(gap, pedal_n, hits, rim, from, scale, friction, standard_gas())
}

/// Opposed plates, their shared contacts and microphones in the same gas state.
pub fn hihat_grid_mu_with_gas(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64, scale: f64, friction: f64, gas: GasState) -> Result<Piece, String> {
    build_hihat(gap, pedal_n, hits, rim, from, scale, friction, gas, None, &FITTED_CYMBAL_INTERFACE)
}

/// The kit's own hi-hat (contact annulus from 0.8 R, offline grid) with the rim law `rim` and the
/// stick, rim waviness and rim friction of `iface`. The card route passes DECLARED_HIHAT_RIM_LAW
/// and DECLARED_CYMBAL_INTERFACE.
pub fn hihat_with_interface(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, iface: &CymbalInterface) -> Result<Piece, String> {
    build_hihat(gap, pedal_n, hits, rim, 0.8, 1.0, iface.hat_rim_friction, standard_gas(), None, iface)
}

/// Material-resolved opposed shells and impactor. Selected material bending
/// loss replaces surface T60; rim contact, rim friction and support remain
/// explicitly separate input data.
pub fn hihat_grid_mu_with_materials(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64, scale: f64, friction: f64, gas: GasState, materials: &KitMaterials) -> Result<Piece, String> {
    build_hihat(gap, pedal_n, hits, rim, from, scale, friction, gas, Some(materials), &FITTED_CYMBAL_INTERFACE)
}

fn build_hihat(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64, scale: f64, friction: f64, gas: GasState, materials: Option<&KitMaterials>, iface: &CymbalInterface) -> Result<Piece, String> {
    if !gap.is_finite() || gap < 0.0 {
        return Err(format!("hi-hat: rim opening must be finite and nonnegative, got {gap} m"));
    }
    let (mut ts, mut bs) = hihat_spec();
    let Law::Power { k, alpha, mu } = rim;
    let rim = Law::Power { k: k * scale * scale, alpha, mu };
    ts.grid *= scale;
    bs.grid *= scale;
    // Displacements remain upward; the bottom's bell and bow face away from the top.
    bs.bow_rise = -bs.bow_rise;
    bs.bell_rise = -bs.bell_rise;
    let mut ins = instrument(OUT_RATE, gas, materials)?;
    let primary = materials.map_or(BRONZE, |m| elastic_values(m.primary()));
    let secondary = materials.map_or(BRONZE, |m| elastic_values(m.secondary()));
    let mut top = cymbal_shell(ts.name, &ts, false, primary, materials.and_then(KitMaterials::primary_bending_loss))?;
    top.origin = [0.0, 0.0, 0.0];
    let mut bot = cymbal_shell(bs.name, &bs, false, secondary, materials.and_then(KitMaterials::secondary_bending_loss))?;
    // Profiles describe midplanes, while the requested opening is between their faces.
    let rim_half_thickness = 0.5 * (bow_thickness(&ts, ts.radius) + bow_thickness(&bs, bs.radius));
    bot.origin = [0.0, 0.0, -rim_half_thickness - gap];
    let ti = ins.add(top);
    let bi = ins.add(bot);
    // the top hangs from the clutch on the pull rod; the bottom sits on the stand's felt cup
    let rod = ins.add(Body::point("rod", 0.3));
    ins.links.push(Link { a: (rod, 0), b: None, k: 2000.0, c: 20.0 });
    // Preserve the physical totals on both coarse and refined discretizations.
    let (r1, _) = ring_nodes(ins.bodies[ti].h, 0.0065, 0.02);
    let clutch: Vec<usize> = (0..ins.bodies[ti].n).filter(|&k| {
        let p = ins.bodies[ti].xy[k];
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        r >= 0.0065 && r <= r1 + 1e-12
    }).collect();
    if clutch.is_empty() { return Err("hi-hat: grid leaves no clutch nodes".into()); }
    let nodes: Vec<_> = clutch.iter().map(|&i| (i,1.0)).collect();
    // Two felt bearings at the estimated mid-radius of the physical annulus.
    // Samejima (2021), doi:10.1250/ast.42.314, measured line stiffness
    // 1.65e5 N/m/m on a crash washer, NOT on this hi-hat. Its transfer here is
    // an estimate, not a fit: each ring gives K_theta = integral(K_line*x^2 ds)
    // = pi*K_line*r^3. Angular damping is unmeasured and remains zero.
    // Keep the existing axial clutch grip, but do not let coarse grid lever
    // arms silently turn it into an increasingly rigid rotational clamp.
    let bearing_radius = 0.5_f64*(0.0065+0.02);
    let angular_k = 2.0*PI*1.65e5*bearing_radius.powi(3);
    ins.links.extend(Link::affine_mount(ti,&ins.bodies[ti],&nodes,Some((rod,0)),
        [8.0e6,80.0],[angular_k,0.0,angular_k],[0.0;3])?);
    // The lower cup has one effective felt bearing, using the same declared
    // line-stiffness estimate. Integrate its axial foundation over the physical
    // annulus, independent of how many grid nodes represent that annulus.
    let felt_area=PI*(0.02_f64.powi(2)-0.0065_f64.powi(2));
    ins.links.extend(Link::affine_mount(bi,&ins.bodies[bi],&nodes,None,
        [9.0e7*felt_area,2.7e3*felt_area],
        [0.5*angular_k,0.0,0.5*angular_k],[0.0;3])?);
    // rim contact between the shells: matching nodes in the outer 20 % of the radius, with the
    // interface's 2-lobe rim waviness (the presets' ~0.15 mm rests the pedal on high spots)
    let mut pairs = 0;
    for k in 0..ins.bodies[ti].n {
        let p = ins.bodies[ti].xy[k];
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        if r < from * ts.radius {
            continue;
        }
        if let Some(kb) = ins.bodies[bi].node_at(p[0], p[1]) {
            let cos2 = if r > 0.0 { (p[0] * p[0] - p[1] * p[1]) / (r * r) } else { 1.0 };
            let (tb, tt) = (bow_thickness(&bs, r), bow_thickness(&ts, r));
            let top_z = ins.bodies[ti].origin[2] + ins.bodies[ti].rest_z[k];
            let bottom_z = ins.bodies[bi].origin[2] + ins.bodies[bi].rest_z[kb];
            let local_gap = top_z - bottom_z - 0.5 * (tt + tb) + iface.hat_rim_waviness_m * (1.0 - cos2);
            ins.contacts.push(Contact::new(vec![(bi, kb, 1.0)], vec![(ti, k, 1.0)], local_gap, rim));
            if friction > 0.0 {
                // the faces the contact presses together rub as the shells bend (see Friction)
                let (tb, tt) = (bow_thickness(&bs, r), bow_thickness(&ts, r));
                let fr = Friction::between(&ins.bodies, ins.contacts.len() - 1, bi, kb, tb, ti, k, tt, friction);
                ins.frictions.push(fr);
            }
            pairs += 1;
        }
    }
    // the pedal pulls the hanging top shell ~1 cm down onto the bottom; let the suspension settle
    let settle = if pedal_n != 0.0 { 0.5 } else { 0.0 };
    if pedal_n != 0.0 {
        ins.drives.push(Drive { body: rod, node: 0, knots: vec![(0.0, 0.0), (0.12, -pedal_n), (1e9, -pedal_n)] });
    }
    let st = ins.add(Body::point("stick", materials.map_or(STICK_TIP_KG, KitMaterials::striker_mass_kg)));
    for &(t, v) in hits {
        let fp = ins.bodies[ti].footprint(0.85 * ts.radius, 0.0, 0.0);
        let tip = iface.hat_stick_tip_m;
        let law = materials.map_or_else(|| hertz_pair(tip, PRESET_HICKORY, (primary.0, primary.1), iface.stick_loss_s_m), KitMaterials::impact_law);
        let mut con = Contact::new(fp.into_iter().map(|(n_, w)| (ti, n_, w)).collect(), vec![(st, 0, 1.0)], 0.0, law);
        con.active = false;
        ins.contacts.push(con);
        ins.strikes.push(Strike { time: t + settle, striker: st, contact: ins.contacts.len() - 1, speed: v, clear: 2e-4, done: false });
    }
    let mut mic = Mic::new([0.3, 0.2, 1.0]);
    mic.unbaffled = Some(ts.radius);
    ins.mics.push(mic);
    let os = pick_oversample(&mut ins, 0.6)?;
    let notes = vec![format!("hi-hat: {} + {} nodes, {pairs} rim contact pairs (friction {friction}), rim opening {:.2} mm, pedal {pedal_n} N, internal {} kHz", ins.bodies[ti].n, ins.bodies[bi].n, gap * 1e3, OUT_RATE * os as f64 / 1000.0),
                     iface.hat_receipt.into()];
    Ok(Piece { instrument: ins, oversample: os, settle, notes })
}

/// Run a piece for `seconds` of output after its settling time; returns 48 kHz pressure [Pa].
pub fn render(p: &mut Piece, seconds: f64, mut progress: impl FnMut(f64)) -> Vec<Vec<f64>> {
    let ins = &mut p.instrument;
    let nm = ins.mics.len();
    let total = ((seconds + p.settle) * ins.rate).round() as usize;
    let skip = (p.settle * ins.rate).round() as usize;
    // Keep the FIR's left support across the settling boundary. No additional
    // physical steps or output latency are introduced.
    let retained_from = skip.saturating_sub(DECIMATION_RADIUS_FRAMES * p.oversample);
    let mut raw = vec![Vec::with_capacity(total - retained_from); nm];
    let mut buf = vec![0.0; nm];
    for s in 0..total {
        ins.step(&mut buf);
        if s >= retained_from {
            for (c, v) in raw.iter_mut().zip(&buf) {
                c.push(*v);
            }
        }
        if s % 4800 == 0 {
            progress(s as f64 / total as f64);
        }
    }
    raw.iter().map(|c| decimate_from(c, p.oversample, skip - retained_from)).collect()
}

#[cfg(test)]
mod settling_render_tests {
    use super::*;

    #[test]
    fn output_is_a_crop_of_the_continuously_filtered_physical_trajectory() {
        let mut p = hihat_scaled(0.002, 0.0, &[(0.0, 2.5)], 8.0).unwrap();
        // The small test mesh normally needs no oversampling. Refine its
        // clock before the first step so this test actually exercises the FIR.
        p.oversample *= 2;
        p.instrument.rate = OUT_RATE * p.oversample as f64;
        let skip_frames = 1024;
        p.settle = skip_frames as f64 / OUT_RATE;
        let seconds = 512.0 / OUT_RATE;
        let mut control = p.instrument.clone();
        let total = ((seconds + p.settle) * control.rate).round() as usize;
        let mut raw = vec![Vec::new(); control.mics.len()];
        let mut buf = vec![0.0; raw.len()];
        for _ in 0..total {
            control.step(&mut buf);
            for (channel, value) in raw.iter_mut().zip(&buf) { channel.push(*value); }
        }
        let actual = render(&mut p, seconds, |_| {});
        assert_eq!(p.instrument.step_index, control.step_index);
        assert!(actual[0][..16].iter().any(|x| x.abs() > 1e-12));
        for (channel, output) in raw.iter().zip(&actual) {
            let continuous = decimate(channel, p.oversample);
            let skip = skip_frames * p.oversample;
            assert!(channel[skip - 16 * p.oversample..skip].iter().any(|x| x.abs() > 1e-12));
            assert_ne!(decimate(&channel[skip..], p.oversample), continuous[skip_frames..]);
            assert_eq!(output, &continuous[skip_frames..]);
        }
    }
}
