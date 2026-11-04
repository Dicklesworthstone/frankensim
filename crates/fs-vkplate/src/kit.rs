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
    /// Seconds of simulation to discard at the start (the snare wires' load coming on, hi-hat closing).
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
    /// Opt-in: the heads' compact dipole radiation (`HeadDipole`). None leaves the render untouched.
    pub head_dipole: Option<HeadDipole>,
}

/// What the build needs for the heads' compact dipole radiation (`crate::RadDipole`).
#[derive(Clone, Copy, Debug)]
pub struct HeadDipole {
    /// The dipole's length [m], from the piece's closed-body table (`head_dipole_from_table`).
    pub length_m: f64,
    /// The air modes' `f_max` [Hz] when they are attached, so the law's frequency carries their mass.
    pub air_modes_f_max: Option<f64>,
}

/// `HeadDipole` from the piece's closed-body table.
pub fn head_dipole_from_table(text: &str, air_modes_f_max: Option<f64>) -> Result<HeadDipole, String> {
    let length_m = crate::closed_body::ClosedBodyTable::parse(text)?.dipole_length()?;
    Ok(HeadDipole { length_m, air_modes_f_max })
}

/// The lower hybrid [Hz] at small amplitude, the frequency the dipole's law is taken at: both heads' (0,1) shapes,
/// each with its own tension, film mass and piston air mass, on the cavity's stiffness, and with `air_modes_f_max` the
/// mass of the axisymmetric air modes below it. Each mode adds `rho c^2 g^2 / (N (w_m^2 - w^2))` to the two heads'
/// mass matrix, with `g` its overlap with the shape, `N` its volume norm, and sign `(-1)^p` on the far head.
/// This is the membrane's law: it leaves out bending, the wires and the hit, and a port (the build refuses one).
pub fn lower_hybrid_hz(s: &DrumSpec, densities: (f64, f64), gas: &GasState, air_modes_f_max: Option<f64>) -> Result<f64, String> {
    use crate::radiation::{ModeShape, bessel_j};
    let j = ModeShape::new(0, 1).zero;
    let (a, depth) = (s.radius, s.depth);
    let (rho, c) = (gas.density, gas.sound_speed);
    let area = PI * a * a;
    let shape_volume = 2.0 * PI * a * a / j;
    let air = 8.0 * rho * a.powi(3) / 3.0 * (shape_volume / area).powi(2);
    let cavity = rho * c * c / (area * depth) * shape_volume * shape_volume;
    let stiffness = |hd: Head| hd.tension * (j / a).powi(2) * area;
    let k = [stiffness(s.batter) + cavity, cavity, stiffness(s.reso) + cavity];
    let m = [densities.0 * s.batter.thickness * area + air, densities.1 * s.reso.thickness * area + air];
    let modes: Vec<(f64, f64, f64)> = air_modes_f_max
        .map_or_else(Vec::new, |f_max| cylinder_air_modes(a, depth, c, f_max))
        .into_iter()
        .filter(|md| md.m == 0)
        .map(|md| {
            let g = 2.0 * PI * a * a * j * bessel_j(0, md.alpha) / (j * j - md.alpha * md.alpha);
            (md.omega, rho * c * c * g * g / (md.lambda * area * depth), if md.p % 2 == 1 { -1.0 } else { 1.0 })
        })
        .collect();
    let lowest = |w: f64| {
        let (mut m11, mut m12, mut m22) = (m[0], 0.0, m[1]);
        for &(wm, weight, sign) in &modes {
            let added = weight / (wm * wm - w * w);
            m11 += added;
            m12 += sign * added;
            m22 += added;
        }
        // det(K - lambda M) = 0 for two coordinates; the smaller root, in the form that keeps its digits.
        let det_m = m11 * m22 - m12 * m12;
        let det_k = k[0] * k[2] - k[1] * k[1];
        let b = m11 * k[2] + m22 * k[0] - 2.0 * m12 * k[1];
        (2.0 * det_k / (b + (b * b - 4.0 * det_m * det_k).sqrt())).sqrt() - w
    };
    let first_mode = modes.iter().map(|md| md.0).fold(f64::INFINITY, f64::min);
    let (mut lo, mut hi) = (2.0 * PI, (2.0 * PI * 5000.0_f64).min(first_mode * (1.0 - 1e-6)));
    if !(lowest(lo) > 0.0 && lowest(hi) < 0.0) {
        return Err(format!("head dipole: no lower hybrid of {} between {} and {} Hz", s.name, lo / (2.0 * PI), hi / (2.0 * PI)));
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if lowest(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi) / (2.0 * PI))
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
    let head_dipole = match s.head_dipole {
        None => None,
        Some(_) if matches!(batter_loss, FilmLoss::FittedT60) => {
            return Err(format!("head dipole: {}'s fitted T60 already lumps its radiation in; it needs the card route", s.name));
        }
        Some(_) if s.port.is_some() => {
            return Err(format!("head dipole: {} has a port, which the law's lower hybrid does not carry", s.name));
        }
        Some(d) => {
            let f_ref = lower_hybrid_hz(s, (primary.2, secondary.2), &gas, d.air_modes_f_max)?;
            Some((d, f_ref, RadDipole::compact_resistance(d.length_m, f_ref, &gas)))
        }
    };
    // Without the fitted T60 (which lumped it in), the heads take their exterior
    // resistance, air mass and viscous layer from the baked table, in this gas,
    // whatever the film-loss source. The viscous layer is always attached, so
    // the load retains every shape the grid resolves.
    if !matches!(batter_loss, FilmLoss::FittedT60) {
        let table = crate::radiation::RadiationTable::shipped()?;
        // (0,1) stays on the compact monopole below; the volume shapes (0,n>=2)
        // take the baffled load less the monopole's share (each receipt says so), and less the
        // dipole's own-head share when it is attached, or the shape would pay that part twice.
        let compact = RadDamp::monopole(vec![], f01(s.radius, s.batter, primary.2), &gas).resistance
            + head_dipole.map_or(0.0, |(_, _, resistance)| RadDipole::OWN_HEAD_SHARE * resistance);
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
    if let Some((d, f_ref, resistance)) = head_dipole {
        ins.rad_dipole = Some(RadDipole { heads: (bi, ri), resistance, radiated_j: 0.0 });
        notes.push(format!(
            "head dipole (opt-in, an estimate): the heads' same-way motion radiates R_d Q^2, Q half the difference of their outward volume velocities, R_d = rho w^4 D^2 / (12 pi c^3) = {resistance:.4e} kg/(m^4 s). D {:.4} m from the closed-body table's low end; w at the membrane-law lower hybrid {f_ref:.2} Hz ({}), fixed, not the rendered line. Each head's own share R_d / 4 comes off the volume shapes' baffled load. Free-field law, no floor",
            d.length_m,
            d.air_modes_f_max.map_or_else(|| "no air-mode mass".to_string(), |f_max| format!("with the axisymmetric air modes below {f_max} Hz")),
        ));
    }
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
        // an odd node count puts every wire node on a head node; with an even one the nodes sit on
        // half cells, and rounding them onto the head grid pairs some with one head node and skips
        // the next
        let n = 2 * ((len / h - 1.0) / 2.0).round() as usize + 1;
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
    // last, on the final clock: the wires reach their load at rest, so only the stick starts the drum
    if settle > 0.0 {
        ins.load_from_rest(settle);
    }
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
        head_dipole: None,
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
        head_dipole: None,
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
        head_dipole: None,
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
// kick 40.6 Hz (kick_OH_F_1, _P_3, 0.3-0.6 s), rack tom 133.9 Hz (hiTom_OH_F_4, _P_3, 1.0-1.6 s),
// floor tom 60.4 Hz (loTom_OH_FF_5, _MP_3, 1.0-1.6 s). The snare is the exception: that read gives
// 134.3 Hz at about 30 dB/s, the same with the wires on or off, which is the rack tom's own pitch
// and rate (133.85-134.22 Hz, 29.8-30.1 dB/s). It is probably the rack tom ringing; that is
// inferred, since no stroke exists with the tom muted. The target is the snare's own line, 150.4 Hz
// at 78-109 dB/s (snare2OFF_OH_F_1, _P_2, 0.10-0.50 s). The card drum's mode is read as above from a
// 0.1 m/s stroke (the Berger stretch negligible), 0.3-1.0 s after onset; both heads keep 88cccbb78's
// ratio. Found by scaling both preset heads together (f proportional to the square root of the
// tension) until that mode was within 0.3 Hz of the target, rounding the batter to four figures and
// setting the resonant head at 88cccbb78's ratio to it; the snare came from three rendered rungs,
// 148.76, 150.60 and 152.42 Hz at 1608, 1651 and 1694 N/m. The rack tom came from two rendered rungs
// and one secant step, with the air modes and the head dipole on: 121.84 and 135.04 Hz at 1264 and
// 1554.7 N/m, then 133.97 Hz at 1530 N/m, against 133.95 Hz on its F and P strokes at 0.80-1.20 s.
// These drums' own modes: snare 150.40 (rendered at 1646), rack tom 133.97 (rendered at 1530), floor
// tom 60.42 and kick 40.56 Hz (trj, 2026-10-26). The early pitch of a hard stroke is up to 15 Hz
// sharper (the stretch); nothing grades balance or decay.
const CARD_DRUMS: [CardDrum; 4] = [
    CardDrum {
        name: "snare 14x6",
        tensions: (1646.0, 1646.0 * 2080.0 / 1330.0),
        strike_xy: (0.06, 0.02),
        provenance: "drum (card route, ruling 7486): tensions tuned to the recording's pitch, lowest mode 150.4 Hz (Salamander snare2OFF_OH_F_1/_P_2, the snare's own line, 0.10-0.50 s; the late-tail 134.3 Hz line is probably the rack tom ringing, inferred, since no stroke exists with the tom muted); strike point (60, 20) mm, a declared playing-position estimate (88cccbb78); stick Hunt-Crossley loss 0.2 s/m and the wires' tension 2.0 N, preload 0.25 N/m and contact loss 2.0 s/m are declared estimates (88cccbb78); the wires' contact stiffness is computed from the cards (Hertz, elliptical, helix; ruling 7520)",
    },
    CardDrum {
        name: "kick 18x18",
        tensions: (375.5, 375.5 * 253.0 / 358.0),
        strike_xy: (0.0, 0.02),
        provenance: "drum (card route, ruling 7486): tensions tuned to the recording's pitch, lowest mode 40.6 Hz (Salamander kick_OH_F_1/_P_3, late tail); strike point (0, 20) mm, a declared playing-position estimate (88cccbb78); beater 0.12 kg with a 3 cm felt face (k 1e8 N/m^2.5, mu 0.3 s/m) and the muffle dashpot 25 N s/m^3 standing in for the pillow are declared estimates (88cccbb78); the overhead microphone is a declared estimate 1 m above the shell's top edge, (0, 1.221, 0.1) m with y up, its 0.1 m toward the drummer there to clear the batter's facing test, and the floor is tangent to the hoops at a declared 0.02 m beyond the shell radius (ruling 8318); the nominal 18x18 cylinder overestimates the air volume, the shells being very thick and unlathed (Salamander README); unknowns: whether the front head is ported (no source says, and none is modelled; a trial 102 mm port raised the first together line from 46.9 to 54.9 Hz and added a together line at 21-25 Hz, which the linear lumped model does not predict, perhaps by the strike's amplitude, untested; ruling 8334), and the one-of-a-kind rope tuning, without lugs, so the heads' tension need not be uniform around the rim",
    },
    CardDrum {
        name: "rack tom 12x7",
        tensions: (1530.0, 1530.0 * 950.0 / 1270.0),
        strike_xy: (0.35 * 0.146, 0.0),
        provenance: "drum (card route, ruling 7486): tensions tuned to the recording's pitch, lowest mode 133.95 Hz with the air modes and head dipole on (Salamander hiTom_OH F and P strokes, 0.80-1.20 s); strike point 0.35 R, a declared playing-position estimate (88cccbb78); stick Hunt-Crossley loss 0.2 s/m, a declared estimate (88cccbb78)",
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
/// law is HIHAT_RIM_LAW. Its rim friction is not fitted: it is the cited BRONZE_RIM_FRICTION
/// (ruling 10862), which replaced 0.5, picked by grade (963d1e278).
pub const FITTED_CYMBAL_INTERFACE: CymbalInterface = CymbalInterface {
    stick_loss_s_m: FITTED_STICK_LOSS_S_M,
    hat_stick_tip_m: HIHAT_STICK_TIP_RADIUS_M,
    hat_rim_waviness_m: 7.5e-5,
    hat_rim_friction: BRONZE_RIM_FRICTION,
    cymbal_receipt: "cymbal stick: fitted to a recording (96ccf8f1b), not physics: Hunt-Crossley loss 0.05 s/m",
    hat_receipt: "hi-hat: fitted to a recording (96ccf8f1b), not physics: rim loss 0.08 s/m (HIHAT_RIM_LAW), a two-lobe rim waviness of up to 0.15 mm, a 4 mm stick tip and, without a stick card, stick loss 0.05 s/m; rim friction 0.6 (Buckley, Cu-10Sn self-mated; ruling 10862): dry in air, a polished hemispherical pin at 2.5 N and 0.05 m/s (NASA TN D-6912, 1972, fig. 2), a kinetic average only, with no static limit or micro-slip data; the alloy, contact geometry and motion differ from B20 rims",
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

/// Dry kinetic friction of tin bronze on itself, for every hi-hat's rims: the card route's (ruling
/// 7520 G1) and, since ruling 10862, the presets' and so live, the browser and GPU, which took
/// 0.5, picked by grade (963d1e278), before. Cu-10 at.% Sn (17 wt% Sn; B20 is 20 wt%) sliding on
/// itself in air, 0.6, read from the bar of Buckley 1972 (NASA TN D-6912, fig. 2: a 0.48 cm
/// hemispherical rider on a disk of the same alloy, polished, 250 g, 300 cm/min, 25 C, averaged
/// over an hour of adhesive wear).
/// It is a kinetic average only: the source gives no static limit and no micro-slip data,
/// yet `Friction` uses this one value as both its stick limit and its sliding force; the
/// alloy, the contact geometry and the motion all differ from a lathed B20 rim. The
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

/// The same with grid spacing multiplied by `scale`.
/// The reference-cell rim force law scales with area to preserve traction under refinement.
pub fn hihat_grid(gap: f64, pedal_n: f64, hits: &[(f64, f64)], rim: Law, from: f64, scale: f64) -> Result<Piece, String> {
    hihat_grid_mu(gap, pedal_n, hits, rim, from, scale, BRONZE_RIM_FRICTION)
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

/// Hear a drum piece's overhead as one closed body on its floor (pexfp.378), from a table baked by
/// fs-bem's `examples/drum_observer.rs` for this drum: refused unless the table's radius, depth and
/// microphone are the spec's bit for bit and `Instrument::attach_closed_body` admits it (clock, gas,
/// floor, head placement). The baffled microphone's batter-only receipt no longer holds; the table's
/// route receipt replaces it.
pub fn attach_closed_body(p: &mut Piece, spec: &DrumSpec, text: &str) -> Result<(), String> {
    let table = crate::closed_body::ClosedBodyTable::parse(text)?;
    if table.radius.to_bits() != spec.radius.to_bits() || table.depth.to_bits() != spec.depth.to_bits() || table.mic != spec.mic {
        return Err(format!(
            "closed-body table {} was baked for radius {} m, depth {} m, microphone {:?}; {} is {} m, {} m, {:?}",
            table.piece, table.radius, table.depth, table.mic, spec.name, spec.radius, spec.depth, spec.mic
        ));
    }
    if p.instrument.mics.len() != 1 {
        return Err(format!("closed-body table {}: a drum piece has one microphone, this one has {}", table.piece, p.instrument.mics.len()));
    }
    let find = |head: &str| p.instrument.bodies.iter().position(|b| b.name == format!("{} {head}", spec.name));
    let (Some(batter), Some(resonant)) = (find("batter"), find("resonant")) else {
        return Err(format!("closed-body table {}: {} has no batter and resonant head bodies", table.piece, spec.name));
    };
    let receipt = p.instrument.attach_closed_body(0, &table, batter, resonant)?;
    p.notes.retain(|n| n != OVERHEAD_HEARS_BATTER_ONLY_RECEIPT);
    p.notes.push(receipt);
    Ok(())
}

/// One acoustic mode of a rigid closed cylinder of radius `a` and depth `L`:
/// `psi = J_m(alpha r / a) {cos, sin}(m theta) cos(p pi z / L)`, `J_m'(alpha) = 0`, at angular frequency
/// `c sqrt((alpha / a)^2 + (p pi / L)^2)`; `lambda` is the volume mean of `psi^2`.
#[derive(Clone, Copy, Debug)]
pub struct CylinderAirMode {
    pub m: usize,
    pub n: usize,
    pub p: usize,
    pub sine: bool,
    pub alpha: f64,
    pub omega: f64,
    pub lambda: f64,
}

/// The positive zeros of `J_m'` below `x_max` (for `m = 0`, those of `J_1`), by a 0.005 scan and bisection.
pub fn bessel_jprime_zeros(m: usize, x_max: f64) -> Vec<f64> {
    use crate::radiation::bessel_j;
    let d = |x: f64| if m == 0 { -bessel_j(1, x) } else { 0.5 * (bessel_j(m - 1, x) - bessel_j(m + 1, x)) };
    let mut out = vec![];
    let start = if m == 0 { 1e-3 } else { m as f64 }; // J_m' > 0 on (0, j'_m1) and j'_m1 > m
    let (mut x0, mut f0) = (start, d(start));
    let mut x = x0;
    while x < x_max {
        x += 5e-3;
        let f1 = d(x);
        if f0 == 0.0 || f0.signum() != f1.signum() {
            let (mut lo, mut hi, mut flo) = (x0, x, f0);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi {
                    break;
                }
                let fm = d(mid);
                if fm.signum() == flo.signum() {
                    lo = mid;
                    flo = fm;
                } else {
                    hi = mid;
                }
            }
            let z = 0.5 * (lo + hi);
            if z < x_max {
                out.push(z);
            }
        }
        x0 = x;
        f0 = f1;
    }
    out
}

/// Every mode of the rigid closed cylinder below `f_max` except the uniform (0, 0, 0) one, which a `Cavity`
/// carries as its lumped term; cosine and sine partners both for `m >= 1`. Sorted by frequency.
pub fn cylinder_air_modes(radius: f64, depth: f64, sound_speed: f64, f_max: f64) -> Vec<CylinderAirMode> {
    use crate::radiation::bessel_j;
    let kmax = 2.0 * PI * f_max / sound_speed;
    let mut out = vec![];
    for m in 0.. {
        let mut alphas: Vec<(usize, f64)> = bessel_jprime_zeros(m, kmax * radius).into_iter().enumerate().map(|(i, a)| (i + 1, a)).collect();
        if m == 0 {
            alphas.insert(0, (0, 0.0));
        }
        if alphas.is_empty() || (m > 0 && alphas[0].1 > kmax * radius) {
            break;
        }
        for &(n, alpha) in &alphas {
            let radial = if m == 0 {
                if n == 0 { 1.0 } else { bessel_j(0, alpha).powi(2) }
            } else {
                0.5 * (1.0 - (m * m) as f64 / (alpha * alpha)) * bessel_j(m, alpha).powi(2)
            };
            for p in 0.. {
                let kz = p as f64 * PI / depth;
                let k2 = (alpha / radius).powi(2) + kz * kz;
                if k2 > kmax * kmax {
                    break;
                }
                if m == 0 && n == 0 && p == 0 {
                    continue;
                }
                let axial = if p == 0 { 1.0 } else { 0.5 };
                let omega = sound_speed * k2.sqrt();
                for sine in if m == 0 { vec![false] } else { vec![false, true] } {
                    out.push(CylinderAirMode { m, n, p, sine, alpha, omega, lambda: radial * axial });
                }
            }
        }
    }
    out.sort_by(|a, b| a.omega.total_cmp(&b.omega).then(a.m.cmp(&b.m)).then(a.p.cmp(&b.p)).then(a.sine.cmp(&b.sine)));
    out
}

/// The viscothermal boundary-layer loss of one rigid-cylinder air mode (`cylinder_air_modes`), as the free
/// mode's amplitude decay rate sigma [1/s] (Q = omega / (2 sigma)).
///
/// The law: A. D. Pierce, Acoustics: An Introduction to Its Physical Principles and Applications (McGraw-Hill
/// 1981; ASA reprint 1989), Sec. 10-4, Eq. (10-4.16). It gives the power the acoustic mode loses per unit
/// area of a rigid wall held at ambient temperature:
/// (omega rho mu / 2)^(1/2) <v_T^2> + (gamma - 1) (omega rho kappa / (2 c_p))^(1/2) <p^2> / (rho c)^2.
/// Here v_T is the tangential velocity and p the pressure just outside the layers. The page is 594 of the
/// retypeset PDF read (jontalle.web.engr.illinois.edu/Public/BOOKS/Pierce-Acoustics.81.pdf); the equation
/// number is the book's.
///
/// Take p = A psi cos(omega t), so v = A grad(psi) sin(omega t) / (rho omega), and the stored energy
/// A^2 int psi^2 dV / (2 rho c^2). The energy balance then gives
/// sigma = [ (c / omega)^2 (omega rho mu / 2)^(1/2) int_S |grad_T psi|^2 dS
///         + (gamma - 1) (omega rho kappa / (2 c_p))^(1/2) int_S psi^2 dS ] / (2 rho int psi^2 dV).
///
/// For psi = J_m(alpha r / a) {cos, sin}(m theta) cos(p pi z / L) with J_m'(alpha) = 0, every integral is
/// closed:
/// - On the shell (r = a): psi^2, and its azimuthal and axial gradients.
/// - On both heads (z = 0 and z = L, where cos^2 = 1): int psi^2 dA = R Theta, and
///   int |grad psi|^2 dA = (alpha / a)^2 R Theta, by Green's identity with the Neumann rim. Here
///   R = int_0^a J_m^2 r dr = a^2 (1 - m^2 / alpha^2) J_m(alpha)^2 / 2, or a^2 / 2 for the axial modes.
///
/// The heads count as rigid walls at ambient temperature. Pierce's criteria (10-4.4a, b) hold for the film:
/// its thermal effusivity is about 10^4 times the air's, and a 0.19 mm film is about 10^3 times
/// (10-4.4b)'s bound at 600 Hz. The layers' share of the heads' own motion is not modelled.
///
/// Values come from the gas record: viscosity, conductivity, c_p, gamma, density and sound speed. Nothing is
/// fitted. A non-finite or negative property is refused; zero viscosity and conductivity give sigma = 0,
/// so Q is infinite.
pub fn cylinder_mode_wall_sigma(md: &CylinderAirMode, radius: f64, depth: f64, gas: &crate::GasState) -> Result<f64, String> {
    let (viscous, thermal, scale) = cylinder_mode_wall_terms(md, radius, depth, gas)?;
    Ok((viscous + thermal) / scale)
}

/// `cylinder_mode_wall_sigma`'s two terms as decay rates [1/s]: (viscous, thermal), the shares of the free
/// mode's sigma from the viscous layer (the tangential velocity, so the air's kinetic energy) and from the
/// thermal layer (the pressure, so its potential energy). Their sum is the sigma to rounding. `AirLayers`
/// steps each share with its half-order law.
pub fn cylinder_mode_wall_sigma_parts(md: &CylinderAirMode, radius: f64, depth: f64, gas: &crate::GasState) -> Result<(f64, f64), String> {
    let (viscous, thermal, scale) = cylinder_mode_wall_terms(md, radius, depth, gas)?;
    Ok((viscous / scale, thermal / scale))
}

/// The two wall-loss integrals of Eq. (10-4.16) for one mode and their common divisor `2 rho int psi^2 dV`.
fn cylinder_mode_wall_terms(md: &CylinderAirMode, radius: f64, depth: f64, gas: &crate::GasState) -> Result<(f64, f64, f64), String> {
    use crate::radiation::bessel_j;
    let (rho, c) = (gas.density, gas.sound_speed);
    let (mu, kap, cp, gamma) = (gas.dynamic_viscosity, gas.thermal_conductivity, gas.specific_heat_cp, gas.gamma);
    for (name, v) in [("dynamic viscosity", mu), ("thermal conductivity", kap)] {
        if !(v.is_finite() && v >= 0.0) {
            return Err(format!("cavity wall loss: the gas record's {name} {v} must be finite and not negative"));
        }
    }
    if !(rho.is_finite() && rho > 0.0 && c.is_finite() && c > 0.0 && cp.is_finite() && cp > 0.0 && gamma.is_finite() && gamma >= 1.0) {
        return Err(format!("cavity wall loss: gas density {rho}, sound speed {c}, c_p {cp}, gamma {gamma} out of range"));
    }
    if !(radius > 0.0 && depth > 0.0 && md.omega.is_finite() && md.omega > 0.0) {
        return Err(format!("cavity wall loss: mode {md:?} in a {radius} x {depth} m cylinder is out of range"));
    }
    let m = md.m as f64;
    let axial_only = md.m == 0 && md.n == 0;
    let jm2 = if axial_only { 1.0 } else { bessel_j(md.m, md.alpha).powi(2) };
    // int over theta of cos^2 (or sin^2) m theta, and of the derivative's sin^2 (or cos^2): equal for m >= 1
    let theta = if md.m == 0 { 2.0 * PI } else { PI };
    let (zc, zs) = if md.p == 0 { (depth, 0.0) } else { (0.5 * depth, 0.5 * depth) };
    let rr = if axial_only {
        0.5 * radius * radius
    } else {
        0.5 * radius * radius * (1.0 - m * m / (md.alpha * md.alpha)) * jm2
    };
    let (kr, kz) = (md.alpha / radius, md.p as f64 * PI / depth);
    let volume_sq = rr * theta * zc;
    let surface_sq = jm2 * radius * theta * zc + 2.0 * rr * theta;
    let surface_grad = jm2 * radius * theta * ((m / radius).powi(2) * zc + kz * kz * zs) + 2.0 * kr * kr * rr * theta;
    let w = md.omega;
    let viscous = (c / w).powi(2) * (w * rho * mu / 2.0).sqrt() * surface_grad;
    let thermal = (gamma - 1.0) * (w * rho * kap / (2.0 * cp)).sqrt() * surface_sq;
    Ok((viscous, thermal, 2.0 * rho * volume_sq))
}

/// Opt-in (ALL-IN, MB 14741): the drum's enclosed air as the rigid cylinder's acoustic modes below `f_max`
/// (`cylinder_air_modes`), beside the lumped cavity, coupling both heads (`AirModes`). Values: the gas record's
/// density and sound speed and the drum's radius and depth; `f_max` is a declared truncation. The shell is
/// rigid and the air lossless here. Returns the receipt line.
pub fn attach_cavity_modes(p: &mut Piece, spec: &DrumSpec, f_max: f64) -> Result<String, String> {
    attach_cavity_modes_with(p, spec, f_max, false)
}

/// `attach_cavity_modes`, and with `wall_loss` (opt-in, MB 15788) each mode also carries the viscothermal
/// boundary-layer loss on shell and heads (`cylinder_mode_wall_sigma_parts`), stepped as Pierce's two layers
/// with their half-order law (`AirLayers`): the viscous layer on the air's kinetic energy, the thermal layer on
/// its potential energy, each loss growing as `sqrt(omega)`. This replaces two earlier forms:
/// - Form B (79988cd3e), a dashpot on the stretch. Below the air modes it charged only about
///   `sigma (omega / omega_m)^4` of the air's kinetic share, a lower bound; above them it overcharged the loss by
///   about `sigma (omega / omega_m)^2` of the air spring's share.
/// - The form landed in 9d3b26242, which damped the air coordinate's own velocity and overcharged it below them.
///
/// Without `wall_loss` the result is
/// `attach_cavity_modes`'s, bit for bit.
pub fn attach_cavity_modes_with(p: &mut Piece, spec: &DrumSpec, f_max: f64, wall_loss: bool) -> Result<String, String> {
    use crate::radiation::bessel_j;
    if !(f_max.is_finite() && f_max > 0.0) {
        return Err(format!("cavity modes: f_max {f_max} Hz must be positive"));
    }
    let find = |head: &str| p.instrument.bodies.iter().position(|b| b.name == format!("{} {head}", spec.name));
    let (Some(batter), Some(resonant)) = (find("batter"), find("resonant")) else {
        return Err(format!("cavity modes: {} has no batter and resonant head bodies", spec.name));
    };
    let (rho, c) = (p.instrument.gas.density, p.instrument.gas.sound_speed);
    let (a, l) = (spec.radius, spec.depth);
    let volume = PI * a * a * l;
    let modes = cylinder_air_modes(a, l, c, f_max);
    if let Some(bad) = modes.iter().find(|md| !(md.lambda > 0.0 && md.omega.is_finite())) {
        return Err(format!("cavity modes: mode {bad:?} has no positive mean square"));
    }
    if modes.is_empty() {
        return Err(format!("cavity modes: none below {f_max} Hz"));
    }
    let k = 1.0 / p.instrument.rate;
    // modes that differ only in their axial order p share one shape on a head: J_m(alpha r / a) {cos, sin}(m th)
    let mut keys: Vec<(usize, usize, bool, f64)> = vec![];
    let mut shape = Vec::with_capacity(modes.len());
    for md in &modes {
        let s = match keys.iter().position(|q| (q.0, q.1, q.2) == (md.m, md.n, md.sine)) {
            Some(s) => s,
            None => {
                keys.push((md.m, md.n, md.sine, md.alpha));
                keys.len() - 1
            }
        };
        shape.push(s);
    }
    let ns = keys.len();
    let mut heads = vec![];
    for (body, at_far_end) in [(batter, false), (resonant, true)] {
        let b = &p.instrument.bodies[body];
        let mut g = vec![0.0; ns * b.n];
        for (s, &(m, _, sine, alpha)) in keys.iter().enumerate() {
            for i in 0..b.n {
                let (x, y) = (b.origin[0] + b.xy[i][0], b.origin[1] + b.xy[i][1]);
                let r = (x * x + y * y).sqrt().min(a);
                let th = y.atan2(x);
                let ang = if m == 0 { 1.0 } else if sine { (m as f64 * th).sin() } else { (m as f64 * th).cos() };
                g[s * b.n + i] = b.h * b.h * bessel_j(m, alpha * r / a) * ang;
            }
        }
        let sign: Vec<f64> = modes.iter().map(|md| if at_far_end && md.p % 2 == 1 { -1.0 } else { 1.0 }).collect();
        heads.push((body, g, sign));
    }
    let kappa: Vec<f64> = modes.iter().map(|md| rho * c * c / (volume * md.lambda)).collect();
    let wt2: Vec<f64> = modes.iter().map(|md| (2.0 - 2.0 * (md.omega * k).cos()) / (k * k)).collect();
    let nm = modes.len();
    let gas = p.instrument.gas.clone();
    let sigma: Vec<f64> = if wall_loss {
        modes.iter().map(|md| cylinder_mode_wall_sigma(md, a, l, &gas)).collect::<Result<_, _>>()?
    } else {
        vec![]
    };
    let layers = if wall_loss {
        let parts: Vec<(f64, f64)> = modes.iter().map(|md| cylinder_mode_wall_sigma_parts(md, a, l, &gas)).collect::<Result<_, _>>()?;
        let (sv, st): (Vec<f64>, Vec<f64>) = parts.into_iter().unzip();
        Some(AirLayers::new(k, &kappa, &wt2, &sv, &st).map_err(|e| format!("cavity wall loss: {e}"))?)
    } else {
        None
    };
    let loss = if let Some(layers) = &layers {
        let q = |i: usize| modes[i].omega / (2.0 * sigma[i]);
        let (lo, hi) = (0..nm).fold((0, 0), |(lo, hi), i| (if q(i) < q(lo) { i } else { lo }, if q(i) > q(hi) { i } else { hi }));
        let band: Vec<usize> = (0..nm).filter(|&i| (2000.0..=8000.0).contains(&(modes[i].omega / (2.0 * PI)))).collect();
        let t60 = |i: usize| 3.0 * std::f64::consts::LN_10 / sigma[i];
        let (t_lo, t_hi) = band.iter().fold((f64::INFINITY, 0.0f64), |(x, y), &i| (x.min(t60(i)), y.max(t60(i))));
        let name = |i: usize| format!("({},{},{}) {:.0} Hz", modes[i].m, modes[i].n, modes[i].p, modes[i].omega / (2.0 * PI));
        let share = |i: usize| layers.sigma_viscous[i] / (layers.sigma_viscous[i] + layers.sigma_thermal[i]);
        let (v_lo, v_hi) = (0..nm).fold((1.0f64, 0.0f64), |(x, y), i| (x.min(share(i)), y.max(share(i))));
        format!(
            "shell rigid; viscothermal wall loss (opt-in, MB 15788): Pierce's boundary-layer loss (Acoustics 1981, Eq. (10-4.16)) on the shell and both heads as rigid walls at ambient temperature, from the gas record (mu {:.4e} Pa s, kappa {:.4e} W/(m K), c_p {:.1} J/(kg K), gamma {:.4}): Q {:.0} at {} to {:.0} at {}, T60 {t_lo:.2}-{t_hi:.2} s for the {} modes at 2-8 kHz. Each mode carries Pierce's two layers with their half-order law (sec. 10-4): the viscous layer on the air's kinetic energy, the thermal layer on its potential energy, each loss growing as sqrt(f) away from the mode's own frequency (viscous share {:.0}-{:.0} %). The half-order operators are {} poles fitted over {:.0}-{:.0} Hz on this clock, within {:.1e} (viscous) and {:.1e} (thermal) of the law; the thermal layer softens the air springs at rest by at most {:.2} %. The layers' reactance lowers each mode by less than the law's sigma / (2 pi) Hz (about 0.7-0.9 of it), since only the loss is fitted. No loss inside the air, none from the heads' own motion in the layers, and no coupling between modes through the layers. This replaces form B (79988cd3e), a dashpot on the stretch that was a lower bound below the air modes and overcharged the loss above them, and 9d3b26242's form, which overcharged it below them",
            gas.dynamic_viscosity, gas.thermal_conductivity, gas.specific_heat_cp, gas.gamma, q(lo), name(lo), q(hi), name(hi), band.len(),
            100.0 * v_lo, 100.0 * v_hi, layers.xi.len(), layers.band.0, layers.band.1, layers.fit_error.0, layers.fit_error.1, 100.0 * layers.softening
        )
    } else {
        "shell rigid, air lossless (no boundary-layer loss: a lower bound)".to_string()
    };
    p.instrument.air_modes.push(AirModes {
        kappa,
        wt2,
        x: vec![0.0; nm],
        xp: vec![0.0; nm],
        shape,
        heads,
        y: vec![0.0; nm],
        q: vec![0.0; ns],
        sigma,
        layers,
    });
    let first: Vec<String> = modes.iter().filter(|md| !md.sine).take(10)
        .map(|md| format!("({},{},{}) {:.0} Hz", md.m, md.n, md.p, md.omega / (2.0 * PI))).collect();
    Ok(format!(
        "cavity modes (opt-in, ALL-IN): the enclosed air as {nm} acoustic modes on {ns} head shapes of a rigid closed cylinder (radius {a} m, depth {l} m; gas rho {rho:.4} kg/m^3, c {c:.2} m/s) from {:.0} to {:.0} Hz (declared truncation {f_max} Hz), coupling both heads beside the lumped (0,0,0) term; first: {}; {loss}",
        modes[0].omega / (2.0 * PI), modes[nm - 1].omega / (2.0 * PI), first.join(", ")
    ))
}

/// The declared reference impact speed of the wire's read (MB 14513 ruling 1, (b)): the log-centre of
/// 0.01-1 m/s; a decade in it moves A by 1.58x.
pub const WIRE_LOSS_REFERENCE_SPEED_M_S: f64 = 0.1;

/// The viscoelastic Hertz retardation time A of one coil turn on the film (.372, MB 14513): the contact's
/// loss tangent (the film card's, times the film's share of the contact compliance) read at the half-period
/// frequency of the elastic Hertz contact of one turn of mass `m_turn` [kg] and stiffness `k_turn`
/// [N/m^1.5] on a rigid film at `v_ref`: t_c = 2.9432 xi0 / v, xi0 = (5 m v^2 / 4K)^(2/5), f = 1 / 2 t_c;
/// A = tan_delta / (2 pi f). Returns (f_ref [Hz], A [s]).
pub fn wire_retardation(m_turn: f64, k_turn: f64, v_ref: f64, tan_delta_contact: f64) -> (f64, f64) {
    let xi0 = (5.0 * m_turn * v_ref * v_ref / (4.0 * k_turn)).powf(0.4);
    let f_ref = 0.5 / (2.9432 * xi0 / v_ref);
    (f_ref, tan_delta_contact / (2.0 * PI * f_ref))
}

/// Opt-in (ALL-IN, MB 14741; .372 rulings 14513): the snare wires' contacts take the viscoelastic Hertz
/// loss (`Contact::retardation_s`) instead of the declared Hunt-Crossley 2 s/m. A comes from the cards
/// (`wire_retardation`): the coil's turn mass and Hertz stiffness, the film card's loss tangent at its
/// band's upper edge (the read lies above it, ruling (b)), the film's share of the contact compliance.
/// One retardation time for shear and bulk is an assumption (with a lossless bulk, Brilliantov's Eq. 4 gives
/// 0.19 tau_E and the E* correspondence 0.89 tau_E). Returns the receipt line.
pub fn attach_viscoelastic_wires(p: &mut Piece, spec: &DrumSpec) -> Result<String, String> {
    let wr = spec.wires.ok_or_else(|| format!("viscoelastic wire contact: {} has no wires", spec.name))?;
    let m = p.instrument.materials.clone().ok_or("viscoelastic wire contact: the drum was built without card materials")?;
    let steel = m.snare_wire().ok_or("viscoelastic wire contact: no sourced snare wire is bound (fitted strings)")?;
    let gas = p.instrument.gas.clone();
    let coil = crate::coil::CoilWire::new(&wr.coil, steel, &gas)?;
    let film = (m.secondary().young_modulus_pa(), m.secondary().poisson_ratio());
    let card = card_drum(spec)?;
    let wires = card.wires.ok_or("viscoelastic wire contact: the card drum has no wires")?;
    let (law, _) = crate::coil::turn_contact(&wr.coil, steel, film, spec.reso.thickness, wr.coil.pitch_m, wires.preload_per_m, 0.0)?;
    let Law::Power { k: k_turn, .. } = law;
    let m_turn = coil.mass_per_m() * wr.coil.pitch_m;
    let (_, band_hi) = m.film_band_hz(1)?;
    let tan_film = m.film_tan_delta(1, band_hi)?.ok_or("viscoelastic wire contact: no film loss tangent at the band edge")?;
    let nu_s = crate::coil::STEEL_POISSON_RATIO_ASSUMED;
    let c_steel = (1.0 - nu_s * nu_s) / steel.young_modulus_pa();
    let c_film = (1.0 - film.1 * film.1) / film.0;
    let share = c_film / (c_steel + c_film);
    let tan_contact = share * tan_film;
    let (f_ref, a) = wire_retardation(m_turn, k_turn, WIRE_LOSS_REFERENCE_SPEED_M_S, tan_contact);
    let wire_bodies: Vec<usize> = p.instrument.bodies.iter().enumerate().filter(|(_, b)| b.name.starts_with("wire ")).map(|(i, _)| i).collect();
    let mut count = 0;
    for c in p.instrument.contacts.iter_mut() {
        if c.b.iter().any(|&(b, _, _)| wire_bodies.contains(&b)) {
            if c.crush.is_some() {
                return Err("viscoelastic wire contact: a wire contact crushes".into());
            }
            let Law::Power { k, alpha, .. } = c.law;
            if alpha != 1.5 {
                return Err(format!("viscoelastic wire contact: law exponent {alpha}, not Hertz"));
            }
            c.law = Law::Power { k, alpha, mu: 0.0 };
            c.retardation_s = Some(a);
            count += 1;
        }
    }
    if count == 0 {
        return Err("viscoelastic wire contact: no wire contacts".into());
    }
    Ok(format!(
        "snare wire contact loss (opt-in, ALL-IN; .372 rulings 14513): viscoelastic Hertz F = K (d^1.5 + 1.5 A d^0.5 d') on all {count} wire contacts, replacing the declared Hunt-Crossley 2 s/m; A = {a:.4e} s = tan delta {tan_contact:.5} / (2 pi {f_ref:.1} Hz): the film card's {tan_film:.4} at its band edge {band_hi:.0} Hz (the read lies above the band) times the film's {share:.4} share of the contact compliance (steel E {:.4e} Pa, nu {nu_s} assumed; film E {:.4e} Pa, nu {}), read at the half-period of one turn ({m_turn:.4e} kg, K {k_turn:.4e} N/m^1.5) on a rigid film at the declared v_ref {WIRE_LOSS_REFERENCE_SPEED_M_S} m/s (the log-centre of 0.01-1 m/s; a decade in v_ref moves A by 1.58x); one retardation time for shear and bulk is an assumption (lossless bulk: Brilliantov Eq. 4 gives 0.19 tau_E, the E* correspondence 0.89 tau_E); the model's contacts (a turn against one free head node) last about half as long, so they carry about twice the loss a read at their own duration would give",
        steel.young_modulus_pa(), film.0, film.1
    ))
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
    // A closed-body observer (pexfp.378) delays its microphone by its latency: run that many more
    // steps and keep each channel's own window, so every channel stays on the physical clock.
    // Without one every latency is 0 and this is the plain loop.
    let latency: Vec<usize> = ins.mics.iter().map(Mic::output_latency).collect();
    let extra = latency.iter().copied().max().unwrap_or(0);
    let mut raw = vec![Vec::with_capacity(total - retained_from); nm];
    let mut buf = vec![0.0; nm];
    for s in 0..total + extra {
        ins.step(&mut buf);
        for ((c, v), &l) in raw.iter_mut().zip(&buf).zip(&latency) {
            if s >= retained_from + l && s < total + l {
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

    /// Wires at rest on a head at rest make no sound until the stick lands. The bound comes from the
    /// card, not from a render. A load F on a film of areal mass sigma gives it a volume acceleration
    /// of at most F / sigma. That reaches a microphone at r as at most rho F / (2 pi r sigma), and the
    /// floor's image can double it. Coming on as a raised cosine over T, the load leaves a mode of
    /// frequency f at most 1 / ((2 f T)^2 - 1) of that, a share that only falls with f, so the
    /// assembly's slowest motion, the wire's own fundamental, is the worst case.
    #[test]
    fn a_snare_under_its_wire_load_is_quiet_before_the_strike() {
        let mut spec = card_drum(&snare()).unwrap();
        spec.grid *= 6.0;
        let wires = spec.wires.unwrap();
        let mut p = drum(&spec, &[]).unwrap();
        let load: f64 = p.instrument.bodies.iter().flat_map(|b| &b.bias).map(|b| b.abs()).sum();
        let sigma = MYLAR.2 * spec.reso.thickness;
        let r = spec.mic.iter().map(|x| x * x).sum::<f64>().sqrt();
        let step_ring = p.instrument.gas().density * load / (PI * r * sigma);
        let slowest = (wires.tension / wires.mass_per_m).sqrt() / (2.0 * wires.length);
        let bound = step_ring / ((2.0 * slowest * p.settle).powi(2) - 1.0);
        let output = render(&mut p, 0.02, |_| {});
        let loudest = output[0].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        eprintln!("{loudest:.3e} Pa before any strike, bound {bound:.3e} Pa");
        assert!(loudest <= bound, "{loudest} Pa before any strike, against {bound} Pa");
    }

    /// The browser kit plays clones, so the ramp's clock travels with the state: a clone taken while
    /// the load comes on finishes it, and both end on exactly the full load.
    #[test]
    fn a_clone_taken_while_the_load_comes_on_finishes_it() {
        let mut spec = snare();
        spec.grid *= 8.0;
        let p = drum(&spec, &[]).unwrap();
        let steps = (p.settle * p.instrument.rate).round() as usize;
        let full: Vec<Vec<f64>> = p.instrument.bodies.iter().map(|b| b.bias.clone()).collect();
        let mut original = p.instrument;
        let mut out = vec![0.0; original.mics.len()];
        for _ in 0..steps / 2 {
            original.step(&mut out);
        }
        assert!(original.load_is_coming_on());
        let mut clone = original.clone();
        for _ in steps / 2..=steps {
            original.step(&mut out);
            clone.step(&mut out);
        }
        assert!(!original.load_is_coming_on() && !clone.load_is_coming_on());
        for ((a, b), f) in original.bodies.iter().zip(&clone.bodies).zip(&full) {
            assert_eq!(&a.bias, f);
            assert_eq!(&b.bias, f);
            assert_eq!(a.w, b.w);
        }
    }
}

#[cfg(test)]
mod wire_alignment_tests {
    use super::*;

    #[test]
    fn every_wire_node_presses_the_head_node_under_it_whatever_the_grid() {
        // The 0.36 m wires span 17.1, 25.7 and 34.3 cells at 6, 4 and 3 times the kit's grid: the
        // rounded counts 17, 26 and 34 put the two even ones on half cells. The kit's own 3.5 mm
        // grid rounds to 103 nodes and keeps them.
        for (scale, nodes) in [(6.0, 17), (4.0, 25), (3.0, 35)] {
            let mut spec = snare();
            spec.grid *= scale;
            let piece = drum(&spec, &[]).unwrap();
            let bodies = &piece.instrument.bodies;
            let wires: Vec<usize> = (0..bodies.len()).filter(|&i| bodies[i].name.starts_with("wire ")).collect();
            assert!(!wires.is_empty());
            for &w in &wires {
                let wire = &bodies[w];
                let mut pressed = vec![];
                for c in piece.instrument.contacts.iter().filter(|c| c.b.iter().any(|&(b, _, _)| b == w)) {
                    let ((head, node, _), (_, k, _)) = (c.a[0], c.b[0]);
                    let (at, under) = (wire.xy[k], bodies[head].xy[node]);
                    assert!(
                        (at[0] - under[0]).abs() < 1e-9 * spec.grid && (at[1] - under[1]).abs() < 1e-9 * spec.grid,
                        "{} node {k} at {at:?} presses the head node at {under:?}",
                        wire.name
                    );
                    pressed.push(node);
                }
                let count = pressed.len();
                pressed.sort_unstable();
                pressed.dedup();
                assert_eq!(pressed.len(), count, "two nodes of {} press one head node", wire.name);
                assert!(2 * count > wire.n, "{} presses the head at {count} of {} nodes", wire.name, wire.n);
                assert_eq!(wire.n, nodes, "{} at {} mm", wire.name, 1e3 * spec.grid);
            }
        }
        let h = snare().grid;
        assert_eq!(2 * ((snare().wires.unwrap().length / h - 1.0) / 2.0).round() as usize + 1, 103);
    }
}

#[cfg(test)]
mod air_mode_tests {
    use super::*;
    use crate::radiation::bessel_j;

    #[test]
    fn jprime_zeros_match_the_tables() {
        // Abramowitz and Stegun table 9.5: j'_{0,1} = 3.8317060, j'_{1,1} = 1.8411838, j'_{1,2} = 5.3314428,
        // j'_{2,1} = 3.0542369, j'_{3,1} = 4.2011889
        let z0 = bessel_jprime_zeros(0, 8.0);
        let z1 = bessel_jprime_zeros(1, 6.0);
        assert!((z0[0] - 3.8317060).abs() < 1e-6 && (z0[1] - 7.0155867).abs() < 1e-6, "{z0:?}");
        assert!((z1[0] - 1.8411838).abs() < 1e-6 && (z1[1] - 5.3314428).abs() < 1e-6, "{z1:?}");
        assert!((bessel_jprime_zeros(2, 4.0)[0] - 3.0542369).abs() < 1e-6);
        assert!((bessel_jprime_zeros(3, 5.0)[0] - 4.2011889).abs() < 1e-6);
    }

    #[test]
    fn mean_squares_match_quadrature() {
        let (a, l) = (0.1703, 0.1524);
        let modes = cylinder_air_modes(a, l, 343.0, 4000.0);
        assert!(modes.len() > 20);
        for md in modes.iter().take(40) {
            // volume mean of psi^2 by a midpoint rule in r, theta, z
            let (nr, nt, nz) = (400, 256, 64);
            let mut s = 0.0;
            for ir in 0..nr {
                let r = a * (ir as f64 + 0.5) / nr as f64;
                let jr = bessel_j(md.m, md.alpha * r / a);
                for it in 0..nt {
                    let th = 2.0 * PI * (it as f64 + 0.5) / nt as f64;
                    let ang = if md.m == 0 { 1.0 } else if md.sine { (md.m as f64 * th).sin() } else { (md.m as f64 * th).cos() };
                    for iz in 0..nz {
                        let z = l * (iz as f64 + 0.5) / nz as f64;
                        let v = jr * ang * (md.p as f64 * PI * z / l).cos();
                        s += v * v * r;
                    }
                }
            }
            let mean = s * (a / nr as f64) * (2.0 * PI / nt as f64) * (l / nz as f64) / (PI * a * a * l);
            assert!((mean / md.lambda - 1.0).abs() < 2e-3, "{md:?}: quadrature {mean}");
        }
    }

    /// A point "head" of mass M on two air modes: the leapfrog of heads and air conserves
    /// E = 0.5 M v^2 + sum 0.5 M_air,m vx_m^2 + sum 0.5 kappa_m (Y+ - X+)(Y - X), M_air = kappa / wt2. The planted
    /// negative steps the same system by hand with the air's coupling sign flipped; its E must drift.
    #[test]
    fn air_modes_conserve_the_discrete_energy_and_a_flipped_coupling_does_not() {
        let (mass, k) = (4.0e-3, 1.0 / 144000.0);
        let kappa = vec![1.0e7, 3.0e7];
        let omega: [f64; 2] = [2.0 * PI * 600.0, 2.0 * PI * 1900.0];
        let wt2: Vec<f64> = omega.iter().map(|w| (2.0 - 2.0 * (w * k).cos()) / (k * k)).collect();
        let g = vec![1.2e-3, -0.7e-3];
        let mut ins = Instrument::new(144000.0);
        let head = ins.add(Body::point("head", mass));
        ins.bodies[head].w[0] = 1e-4;
        ins.bodies[head].wp[0] = 1e-4;
        ins.air_modes.push(AirModes {
            kappa: kappa.clone(),
            wt2: wt2.clone(),
            x: vec![0.0; 2],
            xp: vec![0.0; 2],
            shape: vec![0, 1],
            heads: vec![(head, g.clone(), vec![1.0, 1.0])],
            y: vec![0.0; 2],
            q: vec![0.0; 2],
            sigma: vec![],
            layers: None,
        });
        let energy = |w: f64, wp: f64, x: &[f64], xp: &[f64]| -> f64 {
            let mut e = 0.5 * mass * ((w - wp) / k).powi(2);
            for m in 0..2 {
                e += 0.5 * kappa[m] / wt2[m] * ((x[m] - xp[m]) / k).powi(2) + 0.5 * kappa[m] * (g[m] * w - x[m]) * (g[m] * wp - xp[m]);
            }
            e
        };
        let mut out: [f64; 0] = [];
        ins.step(&mut out);
        let e0 = energy(ins.bodies[head].w[0], ins.bodies[head].wp[0], &ins.air_modes[0].x, &ins.air_modes[0].xp);
        let mut drift: f64 = 0.0;
        for _ in 0..20000 {
            ins.step(&mut out);
            let e = energy(ins.bodies[head].w[0], ins.bodies[head].wp[0], &ins.air_modes[0].x, &ins.air_modes[0].xp);
            drift = drift.max((e / e0 - 1.0).abs());
        }
        eprintln!("air modes: energy drift {drift:.2e} over 20000 steps (e0 {e0:.3e} J)");
        assert!(e0 > 0.0 && drift < 1e-9, "{drift}");
        // planted negative: the same leapfrog by hand, the air driven by -(Y - X)
        let (mut w, mut wp) = (1e-4, 1e-4);
        let (mut x, mut xp) = (vec![0.0; 2], vec![0.0; 2]);
        let mut e_first = None;
        let mut flipped: f64 = 0.0;
        for _ in 0..20000 {
            let mut f = 0.0;
            let mut xn = vec![0.0; 2];
            for m in 0..2 {
                let s = g[m] * w - x[m];
                f += -kappa[m] * s * g[m];
                xn[m] = 2.0 * x[m] - xp[m] - k * k * wt2[m] * s;
            }
            let wn = 2.0 * w - wp + k * k * f / mass;
            (wp, w) = (w, wn);
            (xp, x) = (x, xn);
            let e = energy(w, wp, &x, &xp);
            let e1 = *e_first.get_or_insert(e);
            flipped = flipped.max((e / e1 - 1.0).abs());
        }
        eprintln!("planted negative (air coupling sign flipped): energy drift {flipped:.2e}");
        assert!(flipped > 1e-3, "the energy check cannot see a flipped coupling: {flipped}");
    }

    /// Two point "heads" on three air modes, the first two sharing one shape with opposite signs on the
    /// second head (an odd axial order). Every step, from the instrument's own state, a dense hand step of the
    /// same system (each mode with its own weights) must give the instrument's air coordinates and its head
    /// increments k^2 f / m to rounding (a few ulps of the terms; the increment also to 1e-12 of itself), and a
    /// planted wrong sign on the shared shape must miss by orders. Energy alone cannot see a wrong shape or
    /// sign (those couplings conserve too); a free-running trajectory comparison also carries rounding amplified
    /// by the heads' slow motion over thousands of steps (1.2e-9 of scale here), so the comparison is per step.
    #[test]
    fn the_shape_factored_air_step_matches_a_dense_hand_step() {
        let k = 1.0 / 144000.0;
        let mass = [4.0e-3, 1.5e-3];
        let kappa = [1.0e7, 2.0e7, 3.0e7];
        let wt2: Vec<f64> = [600.0, 1100.0, 1900.0].iter().map(|f: &f64| (2.0 - 2.0 * (2.0 * PI * f * k).cos()) / (k * k)).collect();
        let shape = [0usize, 0, 1];
        let gs = [[1.2e-3, -0.7e-3], [0.9e-3, 0.4e-3]];
        let sign = [[1.0, 1.0, 1.0], [1.0, -1.0, 1.0]];
        let w0 = [1e-4, -0.5e-4];
        let eps = f64::EPSILON;
        // the largest discrepancy of the instrument's step from the dense hand step, in units of its rounding
        // allowance, when the instrument carries `sign_b` on the second head
        let run = |sign_b: [f64; 3]| -> f64 {
            let mut ins = Instrument::new(144000.0);
            let heads = [ins.add(Body::point("a", mass[0])), ins.add(Body::point("b", mass[1]))];
            for j in 0..2 {
                ins.bodies[heads[j]].w[0] = w0[j];
                ins.bodies[heads[j]].wp[0] = w0[j];
            }
            let signs = [sign[0], sign_b];
            ins.air_modes.push(AirModes {
                kappa: kappa.to_vec(),
                wt2: wt2.clone(),
                x: vec![0.0; 3],
                xp: vec![0.0; 3],
                shape: shape.to_vec(),
                heads: (0..2).map(|j| (heads[j], gs[j].to_vec(), signs[j].to_vec())).collect(),
                y: vec![0.0; 3],
                q: vec![0.0; 2],
                sigma: vec![],
                layers: None,
            });
            let g = |j: usize, m: usize| sign[j][m] * gs[j][shape[m]];
            let mut out: [f64; 0] = [];
            let mut worst: f64 = 0.0;
            for _ in 0..5000 {
                let w = [0, 1].map(|j| ins.bodies[heads[j]].w[0]);
                let wp = [0, 1].map(|j| ins.bodies[heads[j]].wp[0]);
                let (x, xp) = (ins.air_modes[0].x.clone(), ins.air_modes[0].xp.clone());
                let mut f = [0.0; 2];
                let mut xn = [0.0; 3];
                for m in 0..3 {
                    let s = g(0, m) * w[0] + g(1, m) * w[1] - x[m];
                    for j in 0..2 {
                        f[j] += -kappa[m] * s * g(j, m);
                    }
                    xn[m] = 2.0 * x[m] - xp[m] + k * k * wt2[m] * s;
                }
                ins.step(&mut out);
                for m in 0..3 {
                    let allow = 8.0 * eps * (2.0 * x[m].abs() + xp[m].abs() + (xn[m] - 2.0 * x[m] + xp[m]).abs()) + f64::MIN_POSITIVE;
                    worst = worst.max((ins.air_modes[0].x[m] - xn[m]).abs() / allow);
                }
                for j in 0..2 {
                    let hand = k * k * f[j] / mass[j];
                    let got = ins.bodies[heads[j]].w[0] - 2.0 * w[j] + wp[j];
                    let allow = 1e-12 * hand.abs() + 8.0 * eps * (2.0 * w[j].abs() + wp[j].abs()) + f64::MIN_POSITIVE;
                    worst = worst.max((got - hand).abs() / allow);
                }
            }
            worst
        };
        let correct = run(sign[1]);
        // planted: the shared shape's odd-order mode carries the even one's sign on the far head
        let planted = run([1.0, 1.0, 1.0]);
        eprintln!("factored vs dense air step, per step in rounding allowances: worst {correct:.3}; planted wrong sign {planted:.3e}");
        assert!(correct <= 1.0, "the factored step leaves the dense step by {correct} rounding allowances");
        assert!(planted >= 1e3, "the per-step check cannot see a wrong sign: {planted}");
    }

    /// The air modes' admission term on the snare (890 modes to 8 kHz) is an upper bound on the largest
    /// eigenvalue of their Gram matrix (built densely, mode by mode, and power-iterated: a lower estimate of
    /// it), and no larger than the trace bound it replaced.
    #[test]
    fn the_air_bound_lies_between_the_gram_eigenvalue_and_the_trace_bound() {
        let spec = snare();
        let mut p = drum(&spec, &[(0.005, 4.5)]).unwrap();
        let rate = p.instrument.rate;
        let before = p.instrument.gas_admission_rate2(rate);
        attach_cavity_modes(&mut p, &spec, 8000.0).unwrap();
        let term = p.instrument.gas_admission_rate2(rate) - before;
        let ins = &p.instrument;
        let air = &ins.air_modes[0];
        let nm = air.kappa.len();
        let mut v: Vec<Vec<f64>> = vec![vec![]; nm];
        for (body, weights, sign) in &air.heads {
            let b = &ins.bodies[*body];
            for kk in 0..nm {
                let s = air.shape[kk];
                let scale = air.kappa[kk].sqrt() * sign[kk];
                v[kk].extend((0..b.n).map(|i| scale * weights[s * b.n + i] / b.mass[i].sqrt()));
            }
        }
        let mut gram = vec![0.0; nm * nm];
        for a in 0..nm {
            for b in a..nm {
                let d: f64 = v[a].iter().zip(&v[b]).map(|(x, y)| x * y).sum();
                gram[a * nm + b] = d;
                gram[b * nm + a] = d;
            }
        }
        let trace_bound = (0..nm).map(|a| gram[a * nm + a]).sum::<f64>() + air.wt2.iter().fold(0.0f64, |x, &y| x.max(y));
        for a in 0..nm {
            gram[a * nm + a] += air.wt2[a];
        }
        let mut u = vec![1.0; nm];
        let mut lam = 0.0;
        for _ in 0..400 {
            let z: Vec<f64> = (0..nm).map(|a| (0..nm).map(|b| gram[a * nm + b] * u[b]).sum()).collect();
            let norm = z.iter().map(|x| x * x).sum::<f64>().sqrt();
            lam = z.iter().zip(&u).map(|(x, y)| x * y).sum::<f64>() / u.iter().map(|x| x * x).sum::<f64>();
            u = z.iter().map(|x| x / norm).collect();
        }
        eprintln!("air bound on the snare: power-iterated lambda {lam:.4e}, admitted term {term:.4e}, trace bound {trace_bound:.4e} ({:.3} of it)", term / trace_bound);
        assert!(lam > 0.0 && term >= lam * (1.0 - 1e-9), "the bound {term} is below the Gram eigenvalue {lam}");
        assert!(term <= trace_bound * (1.0 + 1e-12), "{term} > {trace_bound}");
    }

    /// The snare with air modes to 8 kHz stays admitted (the off path's byte identity is a separate
    /// check: this binary against the pin's with no new flags, `cmp` of the WAVs).
    #[test]
    fn the_snare_with_air_modes_to_8_khz_is_admitted() {
        let spec = snare();
        let mut p = drum(&spec, &[(0.005, 4.5)]).unwrap();
        let (before, _) = p.instrument.stability();
        let receipt = attach_cavity_modes(&mut p, &spec, 8000.0).unwrap();
        let (after, _) = p.instrument.stability();
        assert!(after >= before && after < 0.95, "{before} -> {after}: {receipt}");
    }
}

#[cfg(test)]
mod wire_retardation_tests {
    use super::*;

    #[test]
    fn kat5_a_recomputed_from_the_receipt() {
        // refresh-2 snare receipt: coil 1.9474e-3 kg/m x 3.4 mm pitch, K 1.9094e8 per turn; steel E 2.12e11 (nu 0.3),
        // film E 4.9033e9 (nu 0.38); film card tan delta 0.0372 at its 20 kHz edge. Design (14510): f_ref 23.73 kHz,
        // A 2.44e-7 s.
        let share = ((1.0 - 0.38f64.powi(2)) / 4.9033e9) / ((1.0 - 0.09) / 2.12e11 + (1.0 - 0.38f64.powi(2)) / 4.9033e9);
        let (f, a) = wire_retardation(1.9474e-3 * 3.4e-3, 1.9094e8, WIRE_LOSS_REFERENCE_SPEED_M_S, share * 0.0372);
        eprintln!("KAT5 share {share:.4}, f_ref {f:.1} Hz, A {a:.4e} s");
        assert!((f / 23730.0 - 1.0).abs() < 2e-3, "{f}");
        assert!((a / 2.44e-7 - 1.0).abs() < 5e-3, "{a}");
    }
}

#[cfg(test)]
mod cavity_wall_loss_tests {
    use super::*;
    use crate::radiation::bessel_j;

    const A: f64 = 0.1703;
    const L: f64 = 0.1524;

    /// Pierce's book route, independent of the surface integrals in `cylinder_mode_wall_sigma`, for the
    /// purely axial modes (0,0,p) of the closed cylinder.
    /// - The shell: a wide duct's wall attenuation, Eq. (10-5.8),
    ///   alpha_walls = 2^(-3/2) eta_mu (1 + (gamma - 1) / Pr^(1/2)) L_P / A_c with L_P / A_c = 2 / a.
    ///   The standing wave's two grazing components decay at c alpha_walls.
    /// - The two heads, at normal incidence: Eq. (10-4.21)'s absorption 2^(3/2) eta_kappa takes
    ///   2^(3/2) eta_kappa (A/2)^2 / (2 rho c) per unit area from each.
    ///   That is sigma = 2^(1/2) eta_kappa c / L, with eta from Eq. (10-4.19).
    #[test]
    fn axial_modes_match_pierces_duct_and_normal_incidence_q() {
        let g = crate::standard_gas();
        let (rho, c, mu, kap, cp, gamma) = (g.density, g.sound_speed, g.dynamic_viscosity, g.thermal_conductivity, g.specific_heat_cp, g.gamma);
        let pr = mu * cp / kap;
        let mut worst: f64 = 0.0;
        let mut count = 0;
        for md in cylinder_air_modes(A, L, c, 8000.0).iter().filter(|md| md.m == 0 && md.n == 0) {
            let w = md.omega;
            let eta_mu = (w * mu / (rho * c * c)).sqrt();
            let eta_kappa = (gamma - 1.0) * (w * kap / (rho * c * c * cp)).sqrt();
            let alpha_walls = 2f64.powf(-1.5) * eta_mu * (1.0 + (gamma - 1.0) / pr.sqrt()) * 2.0 / A;
            let book = c * alpha_walls + 2f64.sqrt() * eta_kappa * c / L;
            let code = cylinder_mode_wall_sigma(md, A, L, &g).unwrap();
            eprintln!("axial (0,0,{}) {:.0} Hz: sigma {code:.4} 1/s, book {book:.4}, Q {:.0}", md.p, w / (2.0 * PI), w / (2.0 * code));
            worst = worst.max((code / book - 1.0).abs());
            count += 1;
        }
        eprintln!("axial modes: {count} checked, worst relative difference {worst:.2e}");
        assert!(count >= 7 && worst < 1e-12, "{count} {worst}");
    }

    /// Every other mode, against the same law integrated numerically.
    /// - Midpoint sums over the volume, the shell and both heads.
    /// - Gradients by central differences of psi itself.
    /// This uses none of the closed forms: not R, not Green's identity, not the angular and axial factors.
    /// Bound 1e-3.
    #[test]
    fn general_modes_match_surface_quadrature() {
        let g = crate::standard_gas();
        let (rho, c, mu, kap, cp, gamma) = (g.density, g.sound_speed, g.dynamic_viscosity, g.thermal_conductivity, g.specific_heat_cp, g.gamma);
        let modes = cylinder_air_modes(A, L, c, 6000.0);
        let general: Vec<&CylinderAirMode> = modes.iter().filter(|md| !(md.m == 0 && md.n == 0)).collect();
        // The pick covers every class by construction: m = 0, m >= 1 cosine, m >= 1 sine, each with p = 0 and p >= 1.
        // It takes the lowest and the highest mode of each class below 6 kHz.
        // Then it adds the highest azimuthal order, the highest radial order, and the mode nearest the 5754 Hz line.
        let classes: [fn(&CylinderAirMode) -> bool; 6] = [
            |md| md.m == 0 && md.p == 0,
            |md| md.m == 0 && md.p >= 1,
            |md| md.m >= 1 && !md.sine && md.p == 0,
            |md| md.m >= 1 && md.sine && md.p == 0,
            |md| md.m >= 1 && !md.sine && md.p >= 1,
            |md| md.m >= 1 && md.sine && md.p >= 1,
        ];
        let mut picked: Vec<&CylinderAirMode> = vec![];
        let mut extra: Vec<&CylinderAirMode> = vec![];
        for class in classes {
            let members: Vec<&CylinderAirMode> = general.iter().copied().filter(|md| class(md)).collect();
            assert!(!members.is_empty(), "no mode of a class below 6 kHz");
            extra.push(members[0]);
            extra.push(members[members.len() - 1]);
        }
        extra.push(general.iter().copied().max_by(|x, y| x.m.cmp(&y.m).then(y.omega.total_cmp(&x.omega))).unwrap());
        extra.push(general.iter().copied().max_by(|x, y| x.n.cmp(&y.n).then(y.omega.total_cmp(&x.omega))).unwrap());
        let target = 2.0 * PI * 5754.0;
        extra.push(general.iter().copied().min_by(|x, y| (x.omega - target).abs().total_cmp(&(y.omega - target).abs())).unwrap());
        for md in extra {
            if !picked.iter().any(|q| std::ptr::eq(*q, md)) {
                picked.push(md);
            }
        }
        eprintln!("{} general modes below 6 kHz; {} picked", general.len(), picked.len());
        assert!(picked.len() >= 12, "{}", picked.len());
        assert!(picked.iter().any(|md| md.sine) && picked.iter().any(|md| md.m == 0) && picked.iter().any(|md| md.p >= 1));
        let mut worst: f64 = 0.0;
        for md in picked {
            let kz = md.p as f64 * PI / L;
            let ang = |th: f64| if md.m == 0 { 1.0 } else if md.sine { (md.m as f64 * th).sin() } else { (md.m as f64 * th).cos() };
            let psi = |r: f64, th: f64, z: f64| bessel_j(md.m, md.alpha * r / A) * ang(th) * (kz * z).cos();
            let (nr, nt, nz) = (400usize, 256usize, 256usize);
            let (dr, dt, dz) = (A / nr as f64, 2.0 * PI / nt as f64, L / nz as f64);
            let (hr, ht, hz) = (1e-6 * A, 1e-6, 1e-6 * L);
            // volume, by its separable factors summed numerically
            let sr: f64 = (0..nr).map(|i| { let r = (i as f64 + 0.5) * dr; bessel_j(md.m, md.alpha * r / A).powi(2) * r * dr }).sum();
            let st: f64 = (0..nt).map(|j| ang((j as f64 + 0.5) * dt).powi(2) * dt).sum();
            let sz: f64 = (0..nz).map(|q| ((kz * (q as f64 + 0.5) * dz).cos()).powi(2) * dz).sum();
            let volume = sr * st * sz;
            // the shell, r = a
            let (mut shell_sq, mut shell_grad) = (0.0, 0.0);
            for j in 0..nt {
                let th = (j as f64 + 0.5) * dt;
                for q in 0..nz {
                    let z = (q as f64 + 0.5) * dz;
                    let v = psi(A, th, z);
                    let d_th = (psi(A, th + ht, z) - psi(A, th - ht, z)) / (2.0 * ht) / A;
                    let d_z = (psi(A, th, z + hz) - psi(A, th, z - hz)) / (2.0 * hz);
                    shell_sq += v * v * A * dt * dz;
                    shell_grad += (d_th * d_th + d_z * d_z) * A * dt * dz;
                }
            }
            // both heads, z = 0 and z = L
            let (mut heads_sq, mut heads_grad) = (0.0, 0.0);
            for z in [0.0, L] {
                for i in 0..nr {
                    let r = (i as f64 + 0.5) * dr;
                    for j in 0..nt {
                        let th = (j as f64 + 0.5) * dt;
                        let v = psi(r, th, z);
                        let d_r = (psi(r + hr, th, z) - psi(r - hr, th, z)) / (2.0 * hr);
                        let d_th = (psi(r, th + ht, z) - psi(r, th - ht, z)) / (2.0 * ht) / r;
                        heads_sq += v * v * r * dr * dt;
                        heads_grad += (d_r * d_r + d_th * d_th) * r * dr * dt;
                    }
                }
            }
            let w = md.omega;
            let quad = ((c / w).powi(2) * (w * rho * mu / 2.0).sqrt() * (shell_grad + heads_grad)
                + (gamma - 1.0) * (w * rho * kap / (2.0 * cp)).sqrt() * (shell_sq + heads_sq)) / (2.0 * rho * volume);
            let code = cylinder_mode_wall_sigma(md, A, L, &g).unwrap();
            let rel = (code / quad - 1.0).abs();
            eprintln!("({},{},{}) {} {:.0} Hz: sigma {code:.4} 1/s, quadrature {quad:.4}, rel {rel:.1e}, Q {:.0}",
                md.m, md.n, md.p, if md.sine { "sin" } else { "cos" }, w / (2.0 * PI), w / (2.0 * code));
            worst = worst.max(rel);
        }
        eprintln!("general modes: worst relative difference to quadrature {worst:.2e}");
        assert!(worst < 1e-3, "{worst}");
    }

    /// The stepped air coordinate, heads held: one mode near 5754 Hz rings down at its sigma.
    /// - The mode carries the half-order layers (`AirLayers`) with its viscous and thermal shares.
    /// - The discrete energy E = (kappa / wt2) ((X - X-) / k)^2 / 2 + kappa X X- / 2 falls as exp(-2 sigma t).
    ///   Checked at 144 kHz over 1 / sigma, while the oscillation's energy falls by e^2; bound 1e-3.
    ///   Over longer windows the layers' slow memory (their low poles) is a larger share of what is left.
    /// - Also printed: the ringing frequency against the mode's, read from interpolated zero crossings.
    /// - The same mode without sigma or layers keeps E to 1e-9.
    #[test]
    fn a_held_head_air_mode_rings_down_at_its_sigma() {
        let g = crate::standard_gas();
        let modes = cylinder_air_modes(A, L, g.sound_speed, 8000.0);
        let target = 2.0 * PI * 5754.0;
        let md = modes.iter().min_by(|x, y| (x.omega - target).abs().total_cmp(&(y.omega - target).abs())).unwrap();
        let sigma = cylinder_mode_wall_sigma(md, A, L, &g).unwrap();
        let (sv, st) = cylinder_mode_wall_sigma_parts(md, A, L, &g).unwrap();
        let rate = 144000.0;
        let k = 1.0 / rate;
        let wt2 = (2.0 - 2.0 * (md.omega * k).cos()) / (k * k);
        let kappa = 1.0e7;
        let steps = (rate / sigma).round() as usize;
        let run = |lossy: bool| -> (f64, f64) {
            let mut ins = Instrument::new(rate);
            ins.add(Body::point("held", 1.0));
            let layers = lossy.then(|| AirLayers::new(k, &[kappa], &[wt2], &[sv], &[st]).unwrap());
            ins.air_modes.push(AirModes {
                kappa: vec![kappa], wt2: vec![wt2], x: vec![1.0], xp: vec![1.0], shape: vec![0], heads: vec![],
                y: vec![0.0], q: vec![0.0], sigma: if lossy { vec![sigma] } else { vec![] }, layers,
            });
            let energy = |a: &AirModes| 0.5 * kappa / wt2 * ((a.x[0] - a.xp[0]) / k).powi(2) + 0.5 * kappa * a.x[0] * a.xp[0];
            let mut out: [f64; 0] = [];
            ins.step(&mut out);
            let e0 = energy(&ins.air_modes[0]);
            let mut crossings: Vec<f64> = vec![];
            for i in 0..steps {
                let before = ins.air_modes[0].x[0];
                ins.step(&mut out);
                let after = ins.air_modes[0].x[0];
                if before != 0.0 && (before < 0.0) != (after < 0.0) {
                    crossings.push((i as f64 + before / (before - after)) * k);
                }
            }
            let freq = (crossings.len() - 1) as f64 / (2.0 * (crossings[crossings.len() - 1] - crossings[0]));
            (energy(&ins.air_modes[0]) / e0, freq)
        };
        let (ratio, f_lossy) = run(true);
        let measured = -ratio.ln() / (2.0 * steps as f64 * k);
        let (lossless, f_lossless) = run(false);
        eprintln!(
            "({},{},{}) {:.0} Hz: sigma {sigma:.4} 1/s (viscous {sv:.4}, thermal {st:.4}), stepped {measured:.4} 1/s over {steps} steps (T60 {:.2} s, Q {:.0}); rings at {f_lossy:.3} Hz against {f_lossless:.3} Hz lossless ({:+.3} Hz; sigma / 2 pi = {:.3}); lossless energy ratio {lossless:.12}",
            md.m, md.n, md.p, md.omega / (2.0 * PI), 3.0 * std::f64::consts::LN_10 / sigma, md.omega / (2.0 * sigma), f_lossy - f_lossless, sigma / (2.0 * PI)
        );
        assert!((measured / sigma - 1.0).abs() < 1e-3, "{measured} vs {sigma}");
        assert!((lossless - 1.0).abs() < 1e-9, "{lossless}");
    }

    /// (2) The loss is opt-in.
    /// - Without it, the snare-off's air modes carry no sigma and no layers, and the receipt is the lossless one.
    /// - With it, only sigma and the layers differ: kappa, wt2, the shapes, the weights and the signs are the
    ///   same bits. The layers carry one share pair per mode, summing to its sigma.
    /// - The WAV-level check, this binary against build8's with no new flag, is the build's `cmp`.
    #[test]
    fn the_wall_loss_is_opt_in_and_changes_only_sigma() {
        let spec = snare_wires_off();
        let mut off = drum(&spec, &[(0.005, 4.5)]).unwrap();
        let mut on = drum(&spec, &[(0.005, 4.5)]).unwrap();
        let r_off = attach_cavity_modes(&mut off, &spec, 8000.0).unwrap();
        let r_on = attach_cavity_modes_with(&mut on, &spec, 8000.0, true).unwrap();
        let (a, b) = (&off.instrument.air_modes[0], &on.instrument.air_modes[0]);
        assert!(a.sigma.is_empty() && r_off.ends_with("shell rigid, air lossless (no boundary-layer loss: a lower bound)"), "{r_off}");
        assert!(a.layers.is_none());
        assert_eq!(b.sigma.len(), b.kappa.len());
        assert!(b.sigma.iter().all(|s| s.is_finite() && *s > 0.0));
        let layers = b.layers.as_ref().expect("the wall loss carries its layers");
        assert_eq!((layers.sigma_viscous.len(), layers.sigma_thermal.len()), (b.kappa.len(), b.kappa.len()));
        for m in 0..b.kappa.len() {
            let sum = layers.sigma_viscous[m] + layers.sigma_thermal[m];
            assert!((sum / b.sigma[m] - 1.0).abs() < 1e-14, "mode {m}: shares {sum} vs sigma {}", b.sigma[m]);
        }
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&a.kappa), bits(&b.kappa));
        assert_eq!(bits(&a.wt2), bits(&b.wt2));
        assert_eq!(a.shape, b.shape);
        for (ha, hb) in a.heads.iter().zip(&b.heads) {
            assert_eq!(ha.0, hb.0);
            assert_eq!(bits(&ha.1), bits(&hb.1));
            assert_eq!(bits(&ha.2), bits(&hb.2));
        }
        eprintln!("{r_on}");
    }

    /// (3) Planted negatives.
    /// - A gas with no viscosity and no conductivity gives sigma = 0 exactly, so Q is infinite, not a default.
    /// - Viscosity alone at zero leaves exactly the thermal share.
    /// - A missing (NaN) or negative property is refused.
    #[test]
    fn an_inviscid_nonconducting_gas_gives_infinite_q_and_a_bad_record_is_refused() {
        let air = crate::standard_gas();
        let modes = cylinder_air_modes(A, L, air.sound_speed, 8000.0);
        let mut ideal = air.clone();
        ideal.dynamic_viscosity = 0.0;
        ideal.thermal_conductivity = 0.0;
        for md in &modes {
            let s = cylinder_mode_wall_sigma(md, A, L, &ideal).unwrap();
            assert!(s == 0.0 && (md.omega / (2.0 * s)).is_infinite(), "{md:?}: {s}");
        }
        let mut inviscid = air.clone();
        inviscid.dynamic_viscosity = 0.0;
        for md in modes.iter().filter(|md| md.m == 0 && md.n == 0) {
            let w = md.omega;
            let therm = (air.gamma - 1.0) * (w * air.thermal_conductivity / (2.0 * air.density * air.specific_heat_cp)).sqrt();
            let s = cylinder_mode_wall_sigma(md, A, L, &inviscid).unwrap();
            assert!((s / (therm / A + 2.0 * therm / L) - 1.0).abs() < 1e-12, "{md:?}: {s}");
            let (sv, st) = cylinder_mode_wall_sigma_parts(md, A, L, &inviscid).unwrap();
            assert!(sv == 0.0 && (st / s - 1.0).abs() < 1e-15, "{md:?}: shares {sv} {st} against {s}");
        }
        for (field, v) in [("mu", f64::NAN), ("mu", -1e-5), ("kappa", f64::NAN), ("kappa", -1e-3)] {
            let mut bad = air.clone();
            if field == "mu" { bad.dynamic_viscosity = v } else { bad.thermal_conductivity = v }
            assert!(cylinder_mode_wall_sigma(&modes[0], A, L, &bad).is_err(), "{field} {v} was not refused");
        }
    }

    /// The discrete energy at the half step just stepped to, of a point head (mass, ground spring `kh`) on air modes
    /// with their layers:
    /// E = m v^2 / 2 + kh w w- / 2 + sum_m [ (kappa / wt2) vx^2 / 2 + kappa_rest s s- / 2
    ///     + sum_j ( (G_j / 2) (s- - r_j) (s - r_j) + (c_j / xi_j) u_j^2 / 2 ) ],
    /// with s = g w - x and s- = g w- - x-, and r, u the layers' half-step states. Without layers kappa_rest is
    /// kappa and the sum over j is empty.
    #[allow(clippy::too_many_arguments)]
    fn layered_energy(mass: f64, kh: f64, g: &[f64], w: f64, wp: f64, kappa: &[f64], wt2: &[f64], x: &[f64], xp: &[f64], layers: Option<&AirLayers>, k: f64) -> f64 {
        let mut e = 0.5 * mass * ((w - wp) / k).powi(2) + 0.5 * kh * w * wp;
        for m in 0..kappa.len() {
            let (s, sp) = (g[m] * w - x[m], g[m] * wp - xp[m]);
            e += 0.5 * kappa[m] / wt2[m] * ((x[m] - xp[m]) / k).powi(2);
            match layers {
                None => e += 0.5 * kappa[m] * s * sp,
                Some(l) => {
                    let nj = l.xi.len();
                    e += 0.5 * l.kappa_rest[m] * s * sp;
                    for j in 0..nj {
                        let den = 1.0 + 0.5 * k * l.xi[j];
                        let (gj, cj) = (l.arm[m * nj + j] * den, l.visc[m * nj + j] * den);
                        let (r, u) = (l.r[m * nj + j], l.u[m * nj + j]);
                        e += 0.5 * gj * (sp - r) * (s - r) + 0.5 * cj / l.xi[j] * u * u;
                    }
                }
            }
        }
        e
    }

    /// One air step of `AirModes` with layers, by hand, as `Instrument::step` takes it, for one point head at
    /// displacement `w` (projections `g`). Returns the force on the head, with or without the thermal arms'
    /// share of it (`reaction`).
    #[allow(clippy::too_many_arguments)]
    fn layered_air_step_by_hand(l: &mut AirLayers, kappa: &[f64], wt2: &[f64], g: &[f64], w: f64, x: &mut [f64], xp: &mut [f64], k: f64, reaction: bool) -> f64 {
        let nj = l.xi.len();
        let mut f = 0.0;
        for m in 0..kappa.len() {
            let s = g[m] * w - x[m];
            let base = m * nj;
            let mut p = l.kappa_rest[m] * s;
            for j in 0..nj {
                let a = s - l.r[base + j];
                p += l.arm[base + j] * a;
                l.r[base + j] += l.gain[j] * a;
            }
            let u_sum: f64 = (0..nj).map(|j| l.visc[base + j] * l.u[base + j]).sum();
            let inv_mass = wt2[m] / kappa[m];
            let sk = 0.5 * k * l.visc_sum[m] * inv_mass;
            let next = (2.0 * x[m] - (1.0 - sk) * xp[m] + k * k * inv_mass * (p + u_sum)) / (1.0 + sk);
            let v = (next - xp[m]) / (2.0 * k);
            for j in 0..nj {
                l.u[base + j] += l.gain[j] * (v - l.u[base + j]);
            }
            xp[m] = x[m];
            x[m] = next;
            f -= g[m] * if reaction { p } else { l.kappa_rest[m] * s };
        }
        f
    }

    /// (4) The layers on a moving head. A point head on a ground spring (near 150 Hz) carries two lossy air modes,
    /// at 600 and 1900 Hz, each with its half-order layers.
    /// - The discrete energy is `layered_energy`'s: the lossless one plus each thermal arm's spring and each viscous
    ///   element's mass.
    /// - The viscous elements are exactly passive. Each thermal arm changes the energy in a step by
    ///   -G k xi (a^2 + a (s+ - 2s + s-) / 4), where a = s - r_bar. So it can rise by at most
    ///   sum G k xi (s+ - 2s + s-)^2 / 64, the arms' share of the stretch's second difference.
    /// - Checked every step for 40000 steps at 144 kHz, with an allowance of 1e-12 of E0. The energy must fall.
    /// - Planted negative: the same step by hand, but the head takes only the spring at rest, not the arms' force.
    ///   It must break the bound.
    #[test]
    fn the_layers_on_a_moving_head_rise_at_most_their_second_difference_share_and_a_missing_reaction_does_not() {
        let (mass, rate) = (4.0e-3, 144000.0);
        let k = 1.0 / rate;
        let kappa = [1.0e7, 3.0e7];
        let wt2: Vec<f64> = [600.0, 1900.0].iter().map(|f: &f64| (2.0 - 2.0 * (2.0 * PI * f * k).cos()) / (k * k)).collect();
        let g = [1.2e-3, -0.7e-3];
        let (sv, st) = ([1.2, 2.4], [0.8, 1.6]);
        let kh = mass * (2.0 * PI * 150.0).powi(2);
        let layers = AirLayers::new(k, &kappa, &wt2, &sv, &st).unwrap();
        let nj = layers.xi.len();
        let allowed_rise = |d2: &[f64; 2]| -> f64 {
            (0..2)
                .map(|m| (0..nj).map(|j| layers.arm[m * nj + j] * (1.0 + 0.5 * k * layers.xi[j]) * k * layers.xi[j] * d2[m] * d2[m] / 64.0).sum::<f64>())
                .sum()
        };
        let second = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| [a[0] - 2.0 * b[0] + c[0], a[1] - 2.0 * b[1] + c[1]];
        let mut ins = Instrument::new(rate);
        let head = ins.add(Body::point("head", mass));
        ins.bodies[head].w[0] = 1e-4;
        ins.bodies[head].wp[0] = 1e-4;
        ins.links.push(Link { a: (head, 0), b: None, k: kh, c: 0.0 });
        ins.air_modes.push(AirModes {
            kappa: kappa.to_vec(), wt2: wt2.clone(), x: vec![0.0; 2], xp: vec![0.0; 2], shape: vec![0, 1],
            heads: vec![(head, g.to_vec(), vec![1.0, 1.0])], y: vec![0.0; 2], q: vec![0.0; 2],
            sigma: vec![sv[0] + st[0], sv[1] + st[1]], layers: Some(layers.clone()),
        });
        let mut out: [f64; 0] = [];
        let state = |ins: &Instrument| {
            let (b, a) = (&ins.bodies[head], &ins.air_modes[0]);
            layered_energy(mass, kh, &g, b.w[0], b.wp[0], &a.kappa, &a.wt2, &a.x, &a.xp, a.layers.as_ref(), k)
        };
        let stretch = |ins: &Instrument| {
            let (w, x) = (ins.bodies[head].w[0], &ins.air_modes[0].x);
            [g[0] * w - x[0], g[1] * w - x[1]]
        };
        let mut s_prev = [g[0] * 1e-4, g[1] * 1e-4];
        ins.step(&mut out);
        let mut s_now = stretch(&ins);
        let e0 = state(&ins);
        let (mut e, mut worst) = (e0, f64::NEG_INFINITY);
        for _ in 0..40000 {
            ins.step(&mut out);
            let s_next = stretch(&ins);
            let en = state(&ins);
            worst = worst.max((en - e - allowed_rise(&second(s_next, s_now, s_prev))) / e0);
            e = en;
            (s_prev, s_now) = (s_now, s_next);
        }
        // planted negative: the same step by hand, the head without the thermal arms' force
        let mut l = layers.clone();
        let (mut w, mut wp) = (1e-4, 1e-4);
        let (mut x, mut xp) = (vec![0.0; 2], vec![0.0; 2]);
        let (mut e_prev, mut e_first) = (None, 0.0);
        let mut s_hist = vec![[g[0] * w, g[1] * w]];
        let mut broken = f64::NEG_INFINITY;
        for _ in 0..40001 {
            let f = -kh * w + layered_air_step_by_hand(&mut l, &kappa, &wt2, &g, w, &mut x, &mut xp, k, false);
            let wn = 2.0 * w - wp + k * k * f / mass;
            (wp, w) = (w, wn);
            s_hist.push([g[0] * w - x[0], g[1] * w - x[1]]);
            let en = layered_energy(mass, kh, &g, w, wp, &kappa, &wt2, &x, &xp, Some(&l), k);
            let n = s_hist.len();
            match e_prev {
                Some(ep) => broken = broken.max((en - ep - allowed_rise(&second(s_hist[n - 1], s_hist[n - 2], s_hist[n - 3]))) / e_first),
                None => e_first = en,
            }
            e_prev = Some(en);
        }
        eprintln!(
            "layers, moving head ({} poles): worst step (dE - allowed rise) / E0 {worst:.3e}, E end / E0 {:.6}; planted (no arms' force on the head): {broken:.3e}",
            nj, e / e0
        );
        assert!(worst <= 1e-12, "the energy rose past the arms' second-difference share by {worst} of E0");
        assert!(e < e0, "the loss took no energy: {e} vs {e0}");
        assert!(broken > 1e-9, "the per-step bound cannot see a missing reaction: {broken}");
    }

    /// The energy log-slope decay [1/s] of a run's samples (time, energy) between `t_a` and `t_b`.
    fn energy_decay(samples: &[(f64, f64)], t_a: f64, t_b: f64) -> f64 {
        let pts: Vec<(f64, f64)> = samples.iter().filter(|(t, _)| *t >= t_a && *t <= t_b).map(|(t, e)| (*t, e.ln())).collect();
        let n = pts.len() as f64;
        let (st, se) = (pts.iter().map(|p| p.0).sum::<f64>() / n, pts.iter().map(|p| p.1).sum::<f64>() / n);
        let (mut num, mut den) = (0.0, 0.0);
        for (t, e) in &pts {
            num += (t - st) * (e - se);
            den += (t - st) * (t - st);
        }
        -0.5 * num / den
    }

    /// (5) Below the air modes (MB 17080). The test cylinder's lowest air mode is carried by a point head on a
    /// ground spring. The head sits at a fifth of that mode's frequency, and the air is a fifth of the moving mass.
    /// - Bound: the head's decay must be at or under the Pierce sqrt(omega) bound. That is the mode's whole sigma,
    ///   scaled to the head's frequency by sqrt(omega_h / omega_m), times the air's share of the moving mass.
    /// - Both runs start in the low mode, the air carried quasi-statically, and are read from the discrete
    ///   energy between 0.1 s and 2 s at 144 kHz.
    /// - The half-order layers are stepped by the instrument.
    /// - Planted negative: 9d3b26242's form (2 sigma dX/dt on the air coordinate, the head with the spring alone),
    ///   stepped by hand. It must exceed the bound.
    #[test]
    fn below_the_air_modes_the_head_decays_at_most_at_the_pierce_sqrt_omega_bound_and_the_9d3b26242_form_does_not() {
        let gas = crate::standard_gas();
        let modes = cylinder_air_modes(A, L, gas.sound_speed, 8000.0);
        let md = modes.iter().min_by(|x, y| x.omega.total_cmp(&y.omega)).unwrap();
        let sigma = cylinder_mode_wall_sigma(md, A, L, &gas).unwrap();
        let (sv, st) = cylinder_mode_wall_sigma_parts(md, A, L, &gas).unwrap();
        let rate = 144000.0;
        let k = 1.0 / rate;
        let wt2 = (2.0 - 2.0 * (md.omega * k).cos()) / (k * k);
        let (mass, g) = (4.0e-3, 1.0e-3);
        let air_mass = 0.25 * mass;
        let kappa = air_mass * wt2 / (g * g);
        let wh = md.omega / 5.0;
        let moving = mass + air_mass;
        let kh = wh * wh * moving;
        let share = air_mass / moving;
        let bound = sigma * (wh / md.omega).sqrt() * share;
        let (w0, t1, t2) = (1e-4, 14400usize, 288000usize);
        let x0 = g * w0 * wt2 / (wt2 - wh * wh);
        let energy = |w: f64, wp: f64, x: f64, xp: f64| {
            0.5 * mass * ((w - wp) / k).powi(2) + 0.5 * kh * w * wp + 0.5 * kappa / wt2 * ((x - xp) / k).powi(2) + 0.5 * kappa * (g * w - x) * (g * wp - xp)
        };
        let decay = |e1: f64, e2: f64| -(e2 / e1).ln() / (2.0 * (t2 - t1) as f64 * k);
        let mut ins = Instrument::new(rate);
        let head = ins.add(Body::point("head", mass));
        ins.bodies[head].w[0] = w0;
        ins.bodies[head].wp[0] = w0;
        ins.links.push(Link { a: (head, 0), b: None, k: kh, c: 0.0 });
        ins.air_modes.push(AirModes {
            kappa: vec![kappa], wt2: vec![wt2], x: vec![x0], xp: vec![x0], shape: vec![0],
            heads: vec![(head, vec![g], vec![1.0])], y: vec![0.0], q: vec![0.0], sigma: vec![sigma],
            layers: Some(AirLayers::new(k, &[kappa], &[wt2], &[sv], &[st]).unwrap()),
        });
        let mut out: [f64; 0] = [];
        let state = |ins: &Instrument| {
            let (b, a) = (&ins.bodies[head], &ins.air_modes[0]);
            layered_energy(mass, kh, &[g], b.w[0], b.wp[0], &a.kappa, &a.wt2, &a.x, &a.xp, a.layers.as_ref(), k)
        };
        let mut e1 = 0.0;
        for i in 1..=t2 {
            ins.step(&mut out);
            if i == t1 {
                e1 = state(&ins);
            }
        }
        let law = decay(e1, state(&ins));
        // planted negative: 9d3b26242's form, by hand
        let (mut w, mut wp, mut x, mut xp) = (w0, w0, x0, x0);
        let sk = sigma * k;
        let mut f1 = 0.0;
        for i in 1..=t2 {
            let s = g * w - x;
            let xn = (2.0 * x - (1.0 - sk) * xp + k * k * wt2 * s) / (1.0 + sk);
            let wn = 2.0 * w - wp + k * k * (-kh * w - kappa * s * g) / mass;
            (wp, w, xp, x) = (w, wn, x, xn);
            if i == t1 {
                f1 = energy(w, wp, x, xp);
            }
        }
        let old_form = decay(f1, energy(w, wp, x, xp));
        eprintln!(
            "lowest air mode ({},{},{}) {:.0} Hz, sigma {sigma:.4} 1/s (viscous {sv:.4}, thermal {st:.4}); head at {:.0} Hz, air share {share:.3}: the half-order layers decay at {law:.3e} 1/s, Pierce sqrt(omega) bound {bound:.3e}, 9d3b26242's form {old_form:.3e} (whole sigma x share {:.3e})",
            md.m, md.n, md.p, md.omega / (2.0 * PI), wh / (2.0 * PI), sigma * share
        );
        assert!(law.is_finite() && law > 0.0 && law <= bound, "the layers charge the head {law} 1/s, over the bound {bound}");
        assert!(old_form > bound, "the below-mode check cannot see 9d3b26242's form: {old_form} <= {bound}");
    }

    /// (6) Off resonance, on both sides of a mode. The test cylinder's lowest air mode, with its viscous
    /// and thermal shares, is carried by a point head on a ground spring. The air is a fifth of the head's mass.
    /// The head mode sits at a fifth and at a half of the air mode's frequency, and at twice and four times it.
    /// - The half-order law: the head mode decays at
    ///   sqrt(omega / omega_m) M (sigma_v X^2 + sigma_t (omega_m / omega)^2 (Y - X)^2) / (m w^2 + M X^2).
    ///   Here (w, Y = g w, X) is the lossless eigenvector at the head mode's discrete eigenvalue and M = kappa / wt2.
    ///   This is the viscous layer on the air's kinetic energy and the thermal layer on its potential energy, each
    ///   scaled from the mode's own frequency by sqrt(omega / omega_m).
    /// - Read from `layered_energy`'s log-slope from 5 ms to the smaller of 2 s and 1 / (the law's decay), at
    ///   144 kHz. It must match the law within the fit's worst error plus 2e-3.
    /// - Planted negatives, stepped by hand: 9d3b26242's form below the mode and form B above it each decay at
    ///   more than twice the law.
    #[test]
    fn off_resonance_the_head_follows_the_half_order_law_on_both_sides_and_the_two_earlier_forms_do_not() {
        let gas = crate::standard_gas();
        let modes = cylinder_air_modes(A, L, gas.sound_speed, 8000.0);
        let md = modes.iter().min_by(|x, y| x.omega.total_cmp(&y.omega)).unwrap();
        let sigma = cylinder_mode_wall_sigma(md, A, L, &gas).unwrap();
        let (sv, st) = cylinder_mode_wall_sigma_parts(md, A, L, &gas).unwrap();
        let rate = 144000.0;
        let k = 1.0 / rate;
        let wt2 = (2.0 - 2.0 * (md.omega * k).cos()) / (k * k);
        let (mass, g) = (4.0e-3, 1.0e-3);
        let air_mass = 0.2 * mass;
        let kappa = air_mass * wt2 / (g * g);
        let big_m = kappa / wt2;
        let layers = AirLayers::new(k, &[kappa], &[wt2], &[sv], &[st]).unwrap();
        let tolerance = layers.fit_error.0.max(layers.fit_error.1) + 2e-3;
        // one run: (form, head frequency ratio) -> (measured decay, the law's decay)
        let run = |form: &str, ratio: f64| -> (f64, f64) {
            let wh = md.omega * ratio;
            let lam = (2.0 - 2.0 * (wh * k).cos()) / (k * k);
            let carried = kappa / (kappa - big_m * lam); // X / Y in the lossless eigenvector
            let kh = mass * lam - kappa * g * g * (1.0 - carried);
            let (y, x0) = (g, carried * g);
            let law = (wh / md.omega).sqrt() * big_m * (sv * x0 * x0 + st * (md.omega / wh).powi(2) * (y - x0).powi(2)) / (mass + big_m * x0 * x0);
            let (w0, steps) = (1e-4, (rate * (2.0f64).min(1.0 / law)).round() as usize);
            let (mut w, mut wp) = (w0, w0);
            let (mut x, mut xp) = (vec![x0 * w0], vec![x0 * w0]);
            let mut l = layers.clone();
            let mut samples = vec![];
            for i in 1..=steps {
                let s = g * w - x[0];
                let f = match form {
                    "layers" => -kh * w + layered_air_step_by_hand(&mut l, &[kappa], &[wt2], &[g], w, &mut x, &mut xp, k, true),
                    "9d3b26242" => {
                        let sk = sigma * k;
                        let xn = (2.0 * x[0] - (1.0 - sk) * xp[0] + k * k * wt2 * s) / (1.0 + sk);
                        (xp[0], x[0]) = (x[0], xn);
                        -kh * w - kappa * s * g
                    }
                    _ => {
                        // form B (79988cd3e): a dashpot 2 sigma kappa / wt2 on the stretch, the head taking its reaction
                        let (sk, dy) = (sigma * k, g * (w - wp));
                        let xn = (2.0 * x[0] - (1.0 - sk) * xp[0] + 2.0 * sk * dy + k * k * wt2 * s) / (1.0 + sk);
                        let relative = (0.5 * (xn - xp[0]) - dy) / k;
                        (xp[0], x[0]) = (x[0], xn);
                        -kh * w + g * (-kappa * s + 2.0 * sigma * kappa / wt2 * relative)
                    }
                };
                let wn = 2.0 * w - wp + k * k * f / mass;
                (wp, w) = (w, wn);
                if i % 100 == 0 {
                    let lay = if form == "layers" { Some(&l) } else { None };
                    samples.push((i as f64 * k, layered_energy(mass, kh, &[g], w, wp, &[kappa], &[wt2], &x, &xp, lay, k)));
                }
            }
            (energy_decay(&samples, 0.005, steps as f64 * k), law)
        };
        eprintln!(
            "lowest air mode ({},{},{}) {:.0} Hz, sigma {sigma:.4} 1/s (viscous {sv:.4}, thermal {st:.4}); fit error {:.2e} / {:.2e}; tolerance {tolerance:.2e}",
            md.m, md.n, md.p, md.omega / (2.0 * PI), layers.fit_error.0, layers.fit_error.1
        );
        for ratio in [0.2, 0.5, 2.0, 4.0] {
            let (measured, law) = run("layers", ratio);
            eprintln!("  head at {ratio} x the mode ({:.0} Hz): the layers decay at {measured:.5e} 1/s, the law {law:.5e} ({:+.2e})", ratio * md.omega / (2.0 * PI), measured / law - 1.0);
            assert!((measured / law - 1.0).abs() <= tolerance, "head at {ratio} x the mode: {measured} against the law's {law}");
        }
        let (old_below, law_below) = run("9d3b26242", 0.2);
        let (b_above, law_above) = run("B", 2.0);
        eprintln!("  planted: 9d3b26242's form at 0.2 x decays at {old_below:.4e} ({:.2} x the law); form B at 2 x at {b_above:.4e} ({:.2} x the law)", old_below / law_below, b_above / law_above);
        assert!(old_below > 2.0 * law_below, "the law cannot tell 9d3b26242's form below the mode: {old_below} vs {law_below}");
        assert!(b_above > 2.0 * law_above, "the law cannot tell form B above the mode: {b_above} vs {law_above}");
    }

    /// (7) The half-order fit, and its error over the band.
    /// - On the drums' clocks (96 and 144 kHz) and at 192 kHz, the fit's worst relative error in each layer's loss
    ///   is at most `LAYER_FIT_TOLERANCE` (1 %). The band is 20 Hz to the smaller of 20 kHz and a fifth of the
    ///   rate, and the error is read on a grid ten times finer than the fit's.
    /// - At each mode's own frequency, its stepped coefficients give back its viscous and thermal shares to 1e-12,
    ///   by the scheme's first-order decay. So at resonance the loss is the mode's sigma.
    /// - Refused: a mode below the band, a negative share, mismatched lengths.
    #[test]
    fn the_half_order_fit_is_within_one_percent_over_the_band_and_exact_at_each_mode() {
        for rate in [96000.0, 144000.0, 192000.0] {
            let k = 1.0 / rate;
            let freqs = [60.0, 591.0, 2400.0, 7990.0];
            let wt2: Vec<f64> = freqs.iter().map(|f| (2.0 - 2.0 * (2.0 * PI * f * k).cos()) / (k * k)).collect();
            let kappa = vec![1.0e7; 4];
            let (sv, st) = (vec![0.5, 1.6, 3.0, 6.0], vec![0.3, 1.1, 2.0, 4.0]);
            let l = AirLayers::new(k, &kappa, &wt2, &sv, &st).unwrap();
            eprintln!(
                "{rate} Hz: {} poles over {:.0}-{:.0} Hz; worst relative error in the loss: viscous {:.2e}, thermal {:.2e}; air springs softened at rest by at most {:.3} %",
                l.xi.len(), l.band.0, l.band.1, l.fit_error.0, l.fit_error.1, 100.0 * l.softening
            );
            assert!(l.fit_error.0 <= LAYER_FIT_TOLERANCE && l.fit_error.1 <= LAYER_FIT_TOLERANCE, "{:?}", l.fit_error);
            assert!(l.band.0 == 20.0 && (l.band.1 / (0.2 * rate).min(20000.0) - 1.0).abs() < 1e-12, "band {:?}", l.band);
            let nj = l.xi.len();
            for m in 0..4 {
                let mass = kappa[m] / wt2[m];
                let om = 2.0 / k * (0.5 * k * wt2[m].sqrt()).asin();
                let w = 2.0 / k * (0.5 * om * k).tan();
                let (mut re, mut im) = (0.0, 0.0);
                for j in 0..nj {
                    let (den, d) = (1.0 + 0.5 * k * l.xi[j], w * w + l.xi[j] * l.xi[j]);
                    re += l.visc[m * nj + j] * den * w * w / d;
                    im += l.arm[m * nj + j] * den * w * l.xi[j] / d;
                }
                let (back_v, back_t) = (re / (2.0 * mass), im * k / (2.0 * mass * (om * k).sin()));
                assert!(
                    (back_v / sv[m] - 1.0).abs() < 1e-12 && (back_t / st[m] - 1.0).abs() < 1e-12,
                    "{rate} Hz, mode at {} Hz: shares back {back_v} {back_t} against {} {}",
                    freqs[m], sv[m], st[m]
                );
            }
        }
        let k = 1.0 / 96000.0;
        let wt2 = |f: f64| (2.0 - 2.0 * (2.0 * PI * f * k).cos()) / (k * k);
        assert!(AirLayers::new(k, &[1e7], &[wt2(10.0)], &[1.0], &[1.0]).is_err(), "a mode below the band was not refused");
        assert!(AirLayers::new(k, &[1e7], &[wt2(1000.0)], &[-1.0], &[1.0]).is_err(), "a negative share was not refused");
        assert!(AirLayers::new(k, &[1e7, 1e7], &[wt2(1000.0)], &[1.0], &[1.0]).is_err(), "mismatched lengths were not refused");
    }
}

#[cfg(test)]
mod head_dipole_tests {
    use super::*;
    use crate::radiation::{ModeShape, RadiationTable};

    const RATE: f64 = 48_000.0;
    const RADIUS: f64 = 0.06;
    const DEPTH: f64 = 0.05;
    const GRID: f64 = 0.006;

    fn k() -> f64 {
        1.0 / RATE
    }

    /// Two identical heads of the card's build (Berger stretch, no surface loss) on a small drum, so these run in seconds.
    fn head(name: &str) -> Body {
        head_body(name, RADIUS, GRID, tom(12.0).batter, None, MYLAR, FilmLoss::NoSurfaceLoss).unwrap()
    }

    /// Start a head at rest in the (0,1) shape, `outward` metres per unit of the shape, outward positive.
    fn start(b: &mut Body, outward: f64) {
        let shape = ModeShape::new(0, 1);
        for i in 0..b.n {
            let [x, y] = b.xy[i];
            b.w[i] = outward * shape.value(x / RADIUS, y / RADIUS, false);
            b.wp[i] = b.w[i];
        }
        b.sync();
    }

    fn volume(b: &Body, x: &[f64]) -> f64 {
        x.iter().sum::<f64>() * b.h * b.h
    }

    /// A head's outward volume velocity at the half step.
    fn rate(b: &Body) -> f64 {
        (volume(b, &b.w) - volume(b, &b.wp)) / k()
    }

    fn pair(a_outward: f64, b_outward: f64) -> Instrument {
        let mut ins = Instrument::with_gas(RATE, standard_gas()).unwrap();
        for (name, outward) in [("a", a_outward), ("b", b_outward)] {
            let mut b = head(name);
            start(&mut b, outward);
            ins.add(b);
        }
        ins
    }

    /// The test MB 1067 asks for: the rendered heads lose exactly the energy the term says it removes. Everything a card
    /// head carries that stores or takes energy is on: the Berger stretch, piston air mass, cavity, baffled load, monopole.
    #[test]
    fn the_rendered_heads_lose_exactly_the_energy_the_dipole_reports() {
        let gas = standard_gas();
        let table = RadiationTable::shipped().unwrap();
        let hd = tom(12.0).batter;
        let spec = DrumSpec { radius: RADIUS, depth: DEPTH, batter: hd, reso: hd, ..tom(12.0) };
        let f_ref = lower_hybrid_hz(&spec, (MYLAR.2, MYLAR.2), &gas, None).unwrap();
        let r_d = RadDipole::compact_resistance(0.15, f_ref, &gas);
        let r_m = RadDamp::monopole(vec![], f01(RADIUS, hd, MYLAR.2), &gas).resistance;
        let coincidence_hz = gas.sound_speed.powi(2) / (2.0 * PI * hd.thickness)
            * (12.0 * MYLAR.2 * (1.0 - MYLAR.1 * MYLAR.1) / MYLAR.0).sqrt();
        let mut ins = pair(1.0e-4, -1.0e-4);
        for b in ins.bodies.iter_mut() {
            b.air_mass = 8.0 * gas.density * RADIUS.powi(3) / 3.0;
            b.attach_head_radiation_with_viscous(&table, RADIUS, &gas, Some(r_m + RadDipole::OWN_HEAD_SHARE * r_d), coincidence_hz).unwrap();
        }
        let cavity_volume = PI * RADIUS * RADIUS * DEPTH;
        ins.cavities.push(Cavity { heads: vec![0, 1], volume: cavity_volume, port: None });
        ins.rads.push(RadDamp { bodies: vec![0, 1], resistance: r_m });
        ins.rad_dipole = Some(RadDipole { heads: (0, 1), resistance: r_d, radiated_j: 0.0 });
        assert!(ins.stability().0 < 0.5);
        let net = |ins: &Instrument| rate(&ins.bodies[0]) + rate(&ins.bodies[1]);
        let difference = |ins: &Instrument| 0.5 * (rate(&ins.bodies[0]) - rate(&ins.bodies[1]));
        // Stored energy at the half step, less the backward-velocity correction of the two radiation terms.
        let stored = |ins: &Instrument| -> f64 {
            let heads: f64 = ins.bodies.iter().map(|b| b.discrete_energy(k()) + 0.5 * b.air_mass / (b.area * b.area) * rate(b).powi(2)).sum();
            let (now, before): (f64, f64) = (ins.bodies.iter().map(|b| volume(b, &b.w)).sum(), ins.bodies.iter().map(|b| volume(b, &b.wp)).sum());
            let air = 0.5 * gas.density * gas.sound_speed.powi(2) / cavity_volume * now * before;
            heads + air - 0.25 * k() * (r_m * net(ins).powi(2) + r_d * difference(ins).powi(2))
        };
        let e0 = stored(&ins);
        let (mut net_before, mut difference_before) = (net(&ins), difference(&ins));
        let (mut monopole_lost, mut dipole_lost) = (0.0, 0.0);
        for _ in 0..2000 {
            ins.step(&mut []);
            let (n, d) = (net(&ins), difference(&ins));
            monopole_lost += k() * r_m * (0.5 * (n + net_before)).powi(2);
            dipole_lost += k() * r_d * (0.5 * (d + difference_before)).powi(2);
            (net_before, difference_before) = (n, d);
        }
        let fell = e0 - stored(&ins);
        let load: f64 = ins.bodies.iter().map(Body::head_radiated_j).sum();
        let reported = ins.rad_dipole.as_ref().unwrap().radiated_j;
        assert!((reported / dipole_lost - 1.0).abs() < 1e-12, "the term reports {reported} J, its own definition gives {dipole_lost} J");
        assert!(reported > 1e-4 * e0, "a vacuous run: the dipole took {reported} J of {e0} J");
        let residual = fell - (monopole_lost + load + reported);
        assert!(residual.abs() < 1e-9 * e0, "the heads lost {fell} J; monopole {monopole_lost}, load {load}, dipole {reported}");
        // Planted negative: a term that took a quarter more than it reports would not close.
        assert!((fell - (monopole_lost + load + 0.8 * reported)).abs() > 1e4 * residual.abs());
        eprintln!("energy {e0:.4e} J, fell {fell:.4e} J: monopole {monopole_lost:.3e}, load {load:.3e}, dipole {reported:.3e}, residual {residual:.1e}");
    }

    /// Same-way motion under a dipole of `R_d = 4 R` is, for each head, breathing motion under a monopole of `R`: bit for bit.
    /// Each kind of radiation pays nothing on the other's motion.
    #[test]
    fn the_dipole_pays_on_same_way_motion_and_the_monopole_on_breathing_motion() {
        let r = 3.0e3;
        let run = |a: f64, b: f64, monopole: Option<f64>, dipole: Option<f64>| {
            let mut ins = pair(a, b);
            if let Some(resistance) = monopole {
                ins.rads.push(RadDamp { bodies: vec![0, 1], resistance });
            }
            ins.rad_dipole = dipole.map(|resistance| RadDipole { heads: (0, 1), resistance, radiated_j: 0.0 });
            for _ in 0..1500 {
                ins.step(&mut []);
            }
            ins
        };
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let amp = 1.0e-4;
        let bare_same = run(amp, -amp, None, None);
        let bare_breathing = run(amp, amp, None, None);
        let dipole_same = run(amp, -amp, None, Some(4.0 * r));
        let monopole_breathing = run(amp, amp, Some(r), None);
        assert_ne!(bits(&dipole_same.bodies[0].w), bits(&bare_same.bodies[0].w), "the dipole did nothing to same-way motion");
        assert_eq!(bits(&dipole_same.bodies[0].w), bits(&monopole_breathing.bodies[0].w));
        assert!(dipole_same.rad_dipole.as_ref().unwrap().radiated_j > 0.0);
        let dipole_breathing = run(amp, amp, None, Some(4.0 * r));
        assert_eq!(bits(&dipole_breathing.bodies[0].w), bits(&bare_breathing.bodies[0].w), "breathing motion paid the dipole");
        assert_eq!(dipole_breathing.rad_dipole.as_ref().unwrap().radiated_j, 0.0);
        let monopole_same = run(amp, -amp, Some(r), None);
        assert_eq!(bits(&monopole_same.bodies[0].w), bits(&bare_same.bodies[0].w), "same-way motion paid the monopole");
    }

    /// The grid's own lowest mode: power iteration on `8 I - L`, L the plate's graph Laplacian, from the sampled Bessel shape.
    /// The stepped circle's mode is what decays, and at this grid its `V^2 / M` is 4.4 % above the sampled shape's.
    fn grid_mode(b: &Body) -> Vec<f64> {
        let shape = ModeShape::new(0, 1);
        let node = |x: f64| (x / b.h).round() as i64;
        let index: std::collections::BTreeMap<(i64, i64), usize> = b.xy.iter().enumerate().map(|(n, &[x, y])| ((node(x), node(y)), n)).collect();
        let mut phi: Vec<f64> = b.xy.iter().map(|&[x, y]| shape.value(x / RADIUS, y / RADIUS, false)).collect();
        for _ in 0..4000 {
            let mut next = vec![0.0; b.n];
            for (&(i, j), &a) in &index {
                next[a] = 4.0 * phi[a] + [(1, 0), (-1, 0), (0, 1), (0, -1)].iter().filter_map(|&(di, dj)| index.get(&(i + di, j + dj))).map(|&c| phi[c]).sum::<f64>();
            }
            let norm = next.iter().map(|x| x * x).sum::<f64>().sqrt();
            phi = next.iter().map(|x| x / norm).collect();
        }
        phi
    }

    /// Same-way motion's amplitude falls at `sigma = R_d V^2 / (4 M)`: `Q = V qdot` and each head takes half the force.
    #[test]
    fn same_way_motion_decays_at_the_laws_rate() {
        let sample = head("sample");
        let phi = grid_mode(&sample);
        let shape_volume = volume(&sample, &phi);
        let modal_mass: f64 = phi.iter().zip(&sample.mass).map(|(p, m)| m * p * p).sum();
        let sigma = 25.0;
        let mut ins = pair(0.0, 0.0);
        let amplitude = 1.0e-6 / phi.iter().fold(0.0_f64, |top, p| top.max(p.abs()));
        for (head, sign) in [(0, 1.0), (1, -1.0)] {
            let b = &mut ins.bodies[head];
            for i in 0..b.n {
                b.w[i] = sign * amplitude * phi[i];
                b.wp[i] = b.w[i];
            }
            b.sync();
        }
        ins.rad_dipole = Some(RadDipole { heads: (0, 1), resistance: 4.0 * modal_mass * sigma / (shape_volume * shape_volume), radiated_j: 0.0 });
        let mut peaks = vec![];
        let (mut before, mut now) = (f64::NAN, f64::NAN);
        for step in 0..(0.12 * RATE) as usize {
            ins.step(&mut []);
            let next = volume(&ins.bodies[0], &ins.bodies[0].w).abs();
            if before < now && now >= next {
                peaks.push((step as f64 * k(), now.ln()));
            }
            (before, now) = (now, next);
        }
        let n = peaks.len() as f64;
        let (mt, my) = (peaks.iter().map(|p| p.0).sum::<f64>() / n, peaks.iter().map(|p| p.1).sum::<f64>() / n);
        let slope = peaks.iter().map(|p| (p.0 - mt) * (p.1 - my)).sum::<f64>() / peaks.iter().map(|p| (p.0 - mt).powi(2)).sum::<f64>();
        eprintln!("same-way decay {} 1/s, the law's {sigma}", -slope);
        assert!((-slope / sigma - 1.0).abs() < 0.01, "decay {} vs the law's {sigma}", -slope);
    }

    /// Against an independent script's two-coordinate solve (brentq, numpy eigenvalues), both with and without the air modes' mass.
    #[test]
    fn the_lower_hybrid_matches_its_independent_solve() {
        let gas = standard_gas();
        let card = card_drum(&snare_wires_off()).unwrap();
        let film = (MYLAR.2, MYLAR.2);
        let read = |s: &DrumSpec, f_max| lower_hybrid_hz(s, film, &gas, f_max).unwrap();
        assert!((read(&card, None) - 166.2546).abs() < 2e-3, "{}", read(&card, None));
        assert!((read(&card, Some(8000.0)) - 152.3994).abs() < 2e-3, "{}", read(&card, Some(8000.0)));
        let tuned = retune(card, 1.5);
        assert!((read(&tuned, Some(8000.0)) - 183.5781).abs() < 2e-3, "{}", read(&tuned, Some(8000.0)));
    }

    #[test]
    fn a_fitted_preset_refuses_the_dipole() {
        let mut spec = tom(12.0);
        spec.grid *= 6.0;
        spec.head_dipole = Some(HeadDipole { length_m: 0.37, air_modes_f_max: None });
        let error = drum(&spec, &[]).err().expect("the preset's fitted T60 already lumps its radiation in");
        assert!(error.contains("fitted T60"), "{error}");
    }
}
