//! A closed-open tube under the bar, coupled both ways. The air column's odd plane-wave modes are oscillators
//! on the mouth's volume displacement; each hears the bar through the bar's own Burton-Miller field averaged
//! over the mouth, and by reciprocity the same number is the force of the mouth's flow on each bar mode, so
//! the coupling has no constant of its own.
//!
//! The tube's own frequency and loss are Doutaut, Matignon and Chaigne, JASA 104 (1998): Eq. 7 for the
//! viscothermal wall, Eq. 12 (Levine and Schwinger, unflanged) for the open end. The rigid bar's reflection of
//! the mouth's flow back onto the mouth is solved by the same BEM, which is what lowers a tube close under a bar.
//!
//! With the card's `tube_bore_m` row the bore is the one drawn, and each mode is the same fixed point with the column's
//! length replaced by its plane-wave length in the mouth's units: a walk up the bore, each cylinder's phase at its own
//! Eq. 7 speed, and at each change of radius the pressure and the volume velocity continuous. The open end's correction
//! is Eq. 12's at the mouth, as above. Each mode's mass per unit volume velocity, its wall loss (Eq. 7's viscous term on
//! the mode's kinetic energy, its thermal term on the potential energy, each cylinder at its own radius) and its open-end
//! loss come from the mode's energy along the bore. Drawn as one cylinder, the bore gives the law above bit for bit.
//! Declared omissions of the drawn bore: at a change of radius only the plane wave is carried, so the inertance of the
//! step's evanescent field and the loss on the step's annular face are left out; sloped pieces are cylinders
//! `BORE_CUT_M` long.
//!
//! Declared omissions, none corrected by a constant: the tube's wall is not a scatterer of the bar's field; the
//! mouth radiates to the receiver as a free compact source, with no reflection off the bar; the mouth's flow is a
//! piston; each tube mode's reflection is taken at that mode's frequency and modes aren't coupled through it; the
//! neighbouring bars and the frame are absent.
//!
//! Without the card's observer row the receiver hears the bar itself through the narrow-band law, as a bar alone
//! does. With it, the broadband observer hears the bar: its one-way copy of the bar steps under the same held
//! forces the bar takes, the contact and the tube's reaction on each mode, so it hears the coupled bar and moves
//! nothing, and the mouth's term is added to it as it is to the narrow-band one.
//!
//! Each bar mode exchanges with its nearest tube mode alone, resistance and inertance alike. A transfer held at
//! one tube mode's frequency is only true near that frequency, and with every pair kept A3's tube modes at 2027
//! and 2500 Hz against the bar's 220 Hz mode broke passivity. This rule was chosen after seeing that refusal. A
//! broadband fit per pair is the real cure and hasn't been done.
use crate::bar::TAU;
use crate::card::{Card, Tube};
use crate::render::{Body, strike_budget};
use fs_couple::modal_acoustic_time::broadband::BroadbandModalObserver;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_math::det;

/// Eq. 7's characteristic lengths of viscous and heat-propagation effects [m], as printed.
const VISCOUS_LENGTH_M: f64 = 4e-8;
const THERMAL_LENGTH_M: f64 = 5.6e-8;
/// Eq. 12 holds below this, and the plane-wave column well below its cut-off (Eq. 5: ka = 1.84).
const LARGEST_KA: f64 = 1.5;
/// A sloped piece of a drawn bore is cut into cylinders no longer than this [m], each at its mid radius. On the
/// vibraphone C4 and E4 tubes' rings this cut is within 8e-7 Hz of a 0.002 mm cut on every mode, and a 0.1 mm cut
/// within 7e-5 Hz.
const BORE_CUT_M: f64 = 1e-5;

/// One plane-wave mode of the isolated column.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsolatedMode {
    pub omega_rad_s: f64,
    /// Amplitude decay rate through the wall [1/s], Eq. 7: alpha c.
    pub wall_per_s: f64,
    /// Amplitude decay rate through the open end [1/s]: c Re(Z)/L, with Z normalized by Z0.
    pub radiation_per_s: f64,
    pub ka: f64,
    /// Column length with the open end's correction [m]. Over a drawn bore, the length of a column of the mouth's
    /// radius with this mode's mass per unit volume velocity at the open end's plane.
    pub acoustic_length_m: f64,
}

/// Eq. 12, normalized radiation impedance of an unflanged open end, real and imaginary parts, for ka < 1.5.
pub fn open_end(ka: f64) -> (f64, f64) {
    let (x2, x3, x4) = (ka * ka, ka * ka * ka, ka * ka * ka * ka);
    let log = det::ln(ka);
    (
        x2 / 4.0 + 0.0127 * x4 + 0.082 * x4 * log - 0.023 * x3 * x3,
        0.6133 * ka - 0.036 * x3 + 0.034 * x3 * log - 0.0187 * x4 * ka,
    )
}

/// The column's first `tube.modes` odd modes: each frequency is a fixed point, since the wall slows the wave
/// and the open end's correction both depend on it. With a drawn bore, `drawn_modes`.
pub fn isolated_modes(tube: &Tube, gas: &GasState) -> Result<Vec<IsolatedMode>, String> {
    if tube.bore.is_some() {
        return drawn_modes(tube, gas);
    }
    let c0 = gas.sound_speed;
    let a = tube.radius_m;
    let eps = VISCOUS_LENGTH_M.sqrt() + (gas.gamma - 1.0) * THERMAL_LENGTH_M.sqrt();
    let wall = eps / (a * core::f64::consts::SQRT_2);
    (0..tube.modes)
        .map(|m| {
            let quarter = (2 * m + 1) as f64 * 0.25 * TAU;
            let mut w = quarter * c0 / tube.air_length_m;
            let mut out = None;
            for _ in 0..64 {
                let ka = w * a / c0;
                if ka >= LARGEST_KA {
                    return Err(format!(
                        "tube mode {} reaches ka {ka:.3}: the open-end law holds below {LARGEST_KA}; lower tube_modes",
                        m + 1
                    ));
                }
                let (re, im) = open_end(ka);
                let length = tube.air_length_m + a * im / ka;
                let c = c0 * (1.0 - wall * (c0 / w).sqrt());
                let next = quarter * c / length;
                out = Some(IsolatedMode {
                    omega_rad_s: next,
                    wall_per_s: wall * (next / c0).sqrt() * c,
                    radiation_per_s: c * re / length,
                    ka,
                    acoustic_length_m: length,
                });
                if (next - w).abs() <= 1e-12 * w {
                    return Ok(out.expect("just set"));
                }
                w = next;
            }
            out.ok_or_else(|| "tube mode has no fixed point".to_string())
        })
        .collect()
}

/// The drawn bore as cylinders from the floor face up, [length, radius] [m]: the card's points joined by straight
/// lines, each sloped piece cut into cylinders `BORE_CUT_M` long or less at their mid radii, neighbours of one radius
/// made one. Without the row, the one cylinder of the card's radius.
pub fn cylinders(tube: &Tube) -> Vec<[f64; 2]> {
    let Some(points) = &tube.bore else {
        return vec![[tube.air_length_m, tube.radius_m]];
    };
    let mut out: Vec<[f64; 2]> = Vec::new();
    for pair in points.windows(2) {
        let ([za, ra], [zb, rb]) = (pair[0], pair[1]);
        if zb == za {
            continue;
        }
        let pieces = if ra == rb {
            1
        } else {
            ((zb - za) / BORE_CUT_M).ceil() as usize
        };
        for j in 0..pieces {
            let piece = if ra == rb {
                [zb - za, ra]
            } else {
                let n = pieces as f64;
                [(zb - za) / n, ra + (rb - ra) * (j as f64 + 0.5) / n]
            };
            match out.last_mut() {
                Some(last) if last[1] == piece[1] => last[0] += piece[0],
                _ => out.push(piece),
            }
        }
    }
    out
}

/// One cylinder of the walk up a drawn bore at one frequency.
struct Run {
    /// [m]; the last carries the open end's correction.
    length: f64,
    radius: f64,
    /// Eq. 7's slowed speed at this radius [m/s].
    speed: f64,
    /// Eq. 7's wall factor at this radius.
    wall: f64,
    /// Pressure amplitude, the floor's being 1: p = A cos(phase), Z0 U = A sin(phase).
    amplitude: f64,
    /// Phase at the cylinder's lower end [rad].
    phase: f64,
}

/// `isolated_modes` over the drawn bore (see the header). Each mode is the law's fixed point with the column's length
/// replaced by the walk's length in the mouth's units. The mass, the wall loss and the open-end loss come from the
/// mode's energy over the walk, with the open end's correction carried at the mouth's radius to the plane where the
/// pressure is zero; the mass is per unit volume velocity there, as the law's is.
fn drawn_modes(tube: &Tube, gas: &GasState) -> Result<Vec<IsolatedMode>, String> {
    let c0 = gas.sound_speed;
    let a = tube.radius_m;
    let viscous = VISCOUS_LENGTH_M.sqrt();
    let thermal = (gas.gamma - 1.0) * THERMAL_LENGTH_M.sqrt();
    let eps = viscous + thermal;
    // Eq. 7's viscous share less its thermal share.
    let split = (viscous - thermal) / eps;
    let bore = cylinders(tube);
    let walls: Vec<f64> = bore
        .iter()
        .map(|c| eps / (c[1] * core::f64::consts::SQRT_2))
        .collect();
    let pi = core::f64::consts::PI;
    (0..tube.modes)
        .map(|m| {
            let quarter = (2 * m + 1) as f64 * 0.25 * TAU;
            let mut w = quarter * c0 / tube.air_length_m;
            let mut runs: Vec<Run> = Vec::with_capacity(bore.len());
            for _ in 0..64 {
                let ka = w * a / c0;
                if ka >= LARGEST_KA {
                    return Err(format!(
                        "tube mode {} reaches ka {ka:.3}: the open-end law holds below {LARGEST_KA}; lower tube_modes",
                        m + 1
                    ));
                }
                let (re, im) = open_end(ka);
                let root = (c0 / w).sqrt();
                let mouth = c0 * (1.0 - walls[walls.len() - 1] * root);
                // The walk, its length in the mouth's units.
                runs.clear();
                let (mut length, mut amplitude) = (0.0, 1.0);
                for (i, (c, wall)) in bore.iter().zip(&walls).enumerate() {
                    let speed = c0 * (1.0 - wall * root);
                    if i > 0 {
                        // A change of radius: p and U continuous.
                        let phase = w * length / mouth;
                        let turns = (phase / pi).floor();
                        let (s, co) = det::sin_cos(phase - turns * pi);
                        let below = bore[i - 1][1];
                        let ratio = (below * below) / (c[1] * c[1]);
                        amplitude *= det::hypot(co, ratio * s);
                        length = (turns * pi + det::atan2(ratio * s, co)) * mouth / w;
                    }
                    runs.push(Run {
                        length: c[0],
                        radius: c[1],
                        speed,
                        wall: *wall,
                        amplitude,
                        phase: w * length / mouth,
                    });
                    length += c[0] * (mouth / speed);
                }
                let correction = a * im / ka;
                let column = length + correction;
                let next = quarter * mouth / column;
                if (next - w).abs() <= 1e-12 * w {
                    let last = runs.len() - 1;
                    runs[last].length += correction;
                    let (area, peak) = (a * a, runs[last].amplitude);
                    let (mut kinetic, mut potential) = (0.0, 0.0);
                    let mut parts = Vec::with_capacity(runs.len());
                    for (i, r) in runs.iter().enumerate() {
                        let k = w / r.speed;
                        // The last run ends where the pressure is zero.
                        let top = if i == last {
                            0.0
                        } else {
                            det::sin(2.0 * (r.phase + k * r.length))
                        };
                        let swing = (top - det::sin(2.0 * r.phase)) / (4.0 * k);
                        let weight = (r.radius * r.radius / area)
                            * (r.amplitude / peak)
                            * (r.amplitude / peak);
                        let (ke, pe) = (
                            weight * (r.length / 2.0 - swing),
                            weight * (r.length / 2.0 + swing),
                        );
                        kinetic += ke;
                        potential += pe;
                        parts.push((ke, pe));
                    }
                    let total = kinetic + potential;
                    let mass_length = 2.0 * kinetic;
                    let rise = (next / c0).sqrt();
                    let wall_per_s = runs.iter().zip(&parts).fold(0.0, |sum, (r, (ke, pe))| {
                        sum + r.wall
                            * rise
                            * r.speed
                            * ((ke + pe) / total + split * (ke - pe) / total)
                    });
                    return Ok(IsolatedMode {
                        omega_rad_s: next,
                        wall_per_s,
                        radiation_per_s: mouth * re / mass_length,
                        ka,
                        acoustic_length_m: mass_length,
                    });
                }
                w = next;
            }
            Err(format!(
                "tube mode {} has no fixed point over the drawn bore",
                m + 1
            ))
        })
        .collect()
}

/// A mode's acoustic mass per unit mouth volume velocity [kg/m^4]: rho L / (2 S), with the mode's own acoustic length.
pub fn acoustic_mass(gas: &GasState, tube: &Tube, mode: &IsolatedMode) -> f64 {
    let area = 0.5 * TAU * tube.radius_m * tube.radius_m;
    gas.density * mode.acoustic_length_m / (2.0 * area)
}

/// The mouth's centre: under the bar's middle, mid-width, `gap_m` below the plane of the uncut lower face.
pub fn mouth_centre(card: &Card, tube: &Tube) -> [f64; 3] {
    [
        0.5 * card.length_m,
        0.0,
        -0.5 * card.thickness_m - tube.gap_m,
    ]
}

/// Nineteen equal-area points across the mouth: 1, then rings of 6 and 12.
pub fn mouth_points(centre: [f64; 3], radius: f64) -> Vec<[f64; 3]> {
    let ring = |n: usize, r: f64, start: f64| -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let t = start + TAU * i as f64 / n as f64;
                [
                    centre[0] + r * det::cos(t),
                    centre[1] + r * det::sin(t),
                    centre[2],
                ]
            })
            .collect()
    };
    let mut points = vec![centre];
    points.extend(ring(6, radius * (4.0_f64 / 19.0).sqrt(), 0.0));
    points.extend(ring(12, radius * (13.0_f64 / 19.0).sqrt(), TAU / 24.0));
    points
}

/// The numbers the BEM gives one tube mode at that mode's frequency, exp(-iwt).
#[derive(Clone, Debug, PartialEq)]
pub struct MouthField {
    /// Mouth-mean pressure per unit modal velocity of each bar mode.
    pub per_bar_mode: Vec<C64>,
    /// Mouth-mean pressure reflected by the rigid bar per unit outward mouth volume velocity [kg/(m^4 s)].
    pub reflection: C64,
}

/// One closed-open tube mode as a mass-normalized oscillator with its exact held-force step.
#[derive(Clone, Copy, Debug)]
struct Oscillator {
    omega: f64,
    decay: f64,
    /// 1 / sqrt(acoustic mass including the bar's reflection).
    inverse_root_mass: f64,
    displacement: f64,
    velocity: f64,
    held: f64,
    keep_q: f64,
    keep_v: f64,
    swap_q: f64,
    swap_v: f64,
}

impl Oscillator {
    fn new(omega: f64, decay: f64, inverse_root_mass: f64, h: f64) -> Result<Oscillator, String> {
        if !(decay > 0.0 && decay < omega) {
            return Err(format!(
                "tube mode at {:.2} Hz has decay {decay:.3} /s: it needs positive and underdamped",
                omega / TAU
            ));
        }
        let wd = (omega * omega - decay * decay).sqrt();
        let e = det::exp(-decay * h);
        let (s, c) = (det::sin(wd * h), det::cos(wd * h));
        Ok(Oscillator {
            omega,
            decay,
            inverse_root_mass,
            displacement: 0.0,
            velocity: 0.0,
            held: 0.0,
            keep_q: e * (c + decay / wd * s),
            keep_v: e * s / wd,
            swap_q: -e * omega * omega * s / wd,
            swap_v: e * (c - decay / wd * s),
        })
    }

    fn acceleration(&self) -> f64 {
        self.held - 2.0 * self.decay * self.velocity - self.omega * self.omega * self.displacement
    }

    fn step(&mut self, force: f64) {
        let rest = force / (self.omega * self.omega);
        let q = self.displacement - rest;
        let v = self.velocity;
        self.displacement = rest + self.keep_q * q + self.keep_v * v;
        self.velocity = self.swap_q * q + self.swap_v * v;
        self.held = force;
    }
}

/// True when the symmetric matrix `m` (row-major, n by n) is positive definite.
fn positive_definite(m: &[f64], n: usize) -> bool {
    let mut l = m.to_vec();
    for i in 0..n {
        for j in 0..=i {
            let mut sum = l[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if !(sum > 0.0) {
                    return false;
                }
                l[i * n + i] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    true
}

/// What a card's tube becomes: the tube modes with the bar's reflection in them and the coupling to each bar mode.
#[derive(Clone, Debug, PartialEq)]
pub struct Network {
    /// omega [rad/s], decay [1/s], 1/sqrt(mass) per tube mode.
    pub tube: Vec<[f64; 3]>,
    /// Velocity coupling (mutual resistance) per tube mode, then per bar mode.
    pub resistance: Vec<Vec<f64>>,
    /// Acceleration coupling (mutual inertance), dimensionless, same shape.
    pub inertance: Vec<Vec<f64>>,
}

impl Network {
    /// Build from the isolated modes and the BEM fields. `bar_omega` are the bar's elastic modes [rad/s]; each keeps
    /// its coupling to the tube mode nearest it, with the bar's reflection in, and to no other. See the header.
    pub fn new(
        gas: &GasState,
        tube: &Tube,
        modes: &[IsolatedMode],
        fields: &[MouthField],
        bar_omega: &[f64],
    ) -> Result<Network, String> {
        let mut out = Network {
            tube: Vec::with_capacity(modes.len()),
            resistance: Vec::with_capacity(modes.len()),
            inertance: Vec::with_capacity(modes.len()),
        };
        for (mode, field) in modes.iter().zip(fields) {
            let w0 = mode.omega_rad_s;
            let mass0 = acoustic_mass(gas, tube, mode);
            // exp(-iwt): i dU/dt-free reactance X in p = R U + X dU/dt is -Im(Z)/w.
            let mass = mass0 - field.reflection.im / w0;
            if !(mass > 0.0) {
                return Err("the bar's reflection leaves a tube mode no mass".into());
            }
            let w = w0 * (mass0 / mass).sqrt();
            let decay = ((mode.wall_per_s + mode.radiation_per_s) * 2.0 * mass0
                + field.reflection.re)
                / (2.0 * mass);
            let root = mass.sqrt();
            out.tube.push([w, decay, 1.0 / root]);
            out.resistance
                .push(field.per_bar_mode.iter().map(|g| g.re / root).collect());
            out.inertance.push(
                field
                    .per_bar_mode
                    .iter()
                    .map(|g| -g.im / (w0 * root))
                    .collect(),
            );
        }
        if out
            .resistance
            .iter()
            .any(|row| row.len() != bar_omega.len())
        {
            return Err(
                "the mouth fields and the bar modes disagree on how many bar modes there are"
                    .into(),
            );
        }
        let nearest: Vec<Option<usize>> = bar_omega
            .iter()
            .map(|&w| {
                (0..out.tube.len()).min_by(|&a, &b| {
                    (out.tube[a][0] - w)
                        .abs()
                        .total_cmp(&(out.tube[b][0] - w).abs())
                })
            })
            .collect();
        for (m, (r, x)) in out
            .resistance
            .iter_mut()
            .zip(&mut out.inertance)
            .enumerate()
        {
            for (k, kept) in nearest.iter().enumerate() {
                if *kept != Some(m) {
                    r[k] = 0.0;
                    x[k] = 0.0;
                }
            }
        }
        Ok(out)
    }
}

/// The bar over its tube: the bar's existing modal body, the tube modes beside it, and the coupling between them
/// held for one substep from the states at its start.
pub struct TubeBody {
    body: ModalAcousticTimeModel,
    workspace: ModalAcousticWorkspace,
    tube: Vec<Oscillator>,
    resistance: Vec<Vec<f64>>,
    inertance: Vec<Vec<f64>>,
    bar_omega: Vec<f64>,
    bar_decay: Vec<f64>,
    bar_held: Vec<f64>,
    total: Vec<f64>,
    tube_force: Vec<f64>,
    bar_acceleration: Vec<f64>,
    tube_acceleration: Vec<f64>,
    bar_pressure: f64,
    /// rho / (4 pi r): pressure of a free compact source per unit volume acceleration.
    monopole: f64,
    delay_samples: f64,
    history: Vec<f64>,
    /// The broadband observer, when the card names one: it hears the bar in place of the narrow-band law.
    observer: Option<BroadbandModalObserver>,
    /// The observer's pressure after the substeps just taken [Pa].
    observed: f64,
    /// The narrow-band receiver plus the mouth, kept apart when the observer is the one heard.
    narrow: Vec<f64>,
}

impl TubeBody {
    /// `bar` are the bar's modal oscillators (`elastic` band modes, then rigid coordinates). Refuses a
    /// coupling that isn't passive at this loss: mass and velocity matrices must both be positive definite.
    pub fn new(
        rate_hz: u32,
        bar: &[ModalAcousticMode],
        elastic: usize,
        network: &Network,
        monopole: f64,
        delay_samples: f64,
    ) -> Result<TubeBody, String> {
        let h = 1.0 / f64::from(rate_hz);
        let body = ModalAcousticTimeModel::try_new_with_free_coordinates(
            rate_hz,
            bar.to_vec(),
            strike_budget(),
        )
        .map_err(|e| format!("modal body: {e:?}"))?;
        let workspace = ModalAcousticWorkspace::new(&body);
        let bar_omega: Vec<f64> = bar[..elastic]
            .iter()
            .map(|m| m.angular_frequency_rad_s)
            .collect();
        let bar_decay: Vec<f64> = bar[..elastic]
            .iter()
            .map(|m| m.damping_ratio * m.angular_frequency_rad_s)
            .collect();
        let tube: Vec<Oscillator> = network
            .tube
            .iter()
            .map(|t| Oscillator::new(t[0], t[1], t[2], h))
            .collect::<Result<_, _>>()?;
        let (t, n) = (tube.len(), elastic);
        let size = t + n;
        let mut mass = vec![0.0; size * size];
        let mut damping = vec![0.0; size * size];
        for i in 0..size {
            mass[i * size + i] = 1.0;
        }
        for (k, a) in bar_decay.iter().enumerate() {
            damping[k * size + k] = 2.0 * a;
        }
        for (m, o) in tube.iter().enumerate() {
            damping[(n + m) * size + n + m] = 2.0 * o.decay;
            for k in 0..n {
                let (i, j) = (n + m, k);
                mass[i * size + j] = network.inertance[m][k];
                mass[j * size + i] = network.inertance[m][k];
                damping[i * size + j] = network.resistance[m][k];
                damping[j * size + i] = network.resistance[m][k];
            }
        }
        if !positive_definite(&mass, size) {
            return Err("tube coupling refused: the mutual inertance isn't passive".into());
        }
        if !positive_definite(&damping, size) {
            return Err(
                "tube coupling refused: the mutual resistance isn't passive at this loss".into(),
            );
        }
        Ok(TubeBody {
            body,
            workspace,
            tube_force: vec![0.0; t],
            tube_acceleration: vec![0.0; t],
            tube,
            resistance: network.resistance.clone(),
            inertance: network.inertance.clone(),
            bar_omega,
            bar_decay,
            bar_held: vec![0.0; n],
            total: vec![0.0; bar.len()],
            bar_acceleration: vec![0.0; n],
            bar_pressure: 0.0,
            monopole,
            delay_samples,
            history: Vec::new(),
            observer: None,
            observed: 0.0,
            narrow: Vec::new(),
        })
    }

    /// The bar's own modal model, built as a bar alone builds it: the broadband observer is fitted on it.
    pub fn model(&self) -> &ModalAcousticTimeModel {
        &self.body
    }

    /// Hear the bar through `observer` in place of the narrow-band law. The narrow-band receiver, plus the
    /// mouth, is still computed and kept apart as the check that the bar and the tube did not change.
    pub fn observe(&mut self, observer: BroadbandModalObserver) {
        self.observer = Some(observer);
    }

    /// Mechanical energy: both sets of modes plus the mutual kinetic term.
    #[cfg(test)]
    pub fn energy(&self) -> f64 {
        let states = self.body.states();
        let mut e = 0.0;
        for (k, w) in self.bar_omega.iter().enumerate() {
            let s = &states[k];
            e += 0.5 * s.velocity_m_sqrt_kg_per_s * s.velocity_m_sqrt_kg_per_s;
            e += 0.5 * w * w * s.displacement_m_sqrt_kg * s.displacement_m_sqrt_kg;
        }
        for (m, o) in self.tube.iter().enumerate() {
            e += 0.5 * o.velocity * o.velocity;
            e += 0.5 * o.omega * o.omega * o.displacement * o.displacement;
            for (k, x) in self.inertance[m].iter().enumerate() {
                e += x * states[k].velocity_m_sqrt_kg_per_s * o.velocity;
            }
        }
        e
    }

    /// Power lost to the two loss laws and the mutual resistance at this instant [W per unit modal mass].
    #[cfg(test)]
    pub fn dissipation(&self) -> f64 {
        let states = self.body.states();
        let mut p = 0.0;
        for (k, a) in self.bar_decay.iter().enumerate() {
            let v = states[k].velocity_m_sqrt_kg_per_s;
            p += 2.0 * a * v * v;
        }
        for (m, o) in self.tube.iter().enumerate() {
            p += 2.0 * o.decay * o.velocity * o.velocity;
            for (k, r) in self.resistance[m].iter().enumerate() {
                p += 2.0 * r * states[k].velocity_m_sqrt_kg_per_s * o.velocity;
            }
        }
        p
    }

    /// The tube modes' velocities, for a test.
    #[cfg(test)]
    pub fn tube_velocities(&self) -> Vec<f64> {
        self.tube.iter().map(|o| o.velocity).collect()
    }

    /// The bar modes' velocities, for a test.
    #[cfg(test)]
    pub fn bar_velocities(&self) -> Vec<f64> {
        self.body
            .states()
            .iter()
            .take(self.bar_omega.len())
            .map(|s| s.velocity_m_sqrt_kg_per_s)
            .collect()
    }

    /// Volume acceleration of the mouth, `delay_samples` ago, linearly interpolated.
    fn delayed_acceleration(&self) -> f64 {
        let whole = self.delay_samples.floor();
        let frac = self.delay_samples - whole;
        let at = |back: usize| {
            self.history
                .len()
                .checked_sub(1 + back)
                .map_or(0.0, |i| self.history[i])
        };
        let n = whole as usize;
        (1.0 - frac) * at(n) + frac * at(n + 1)
    }
}

impl Body for TubeBody {
    fn strike_displacement(&self, projection: &[f64]) -> f64 {
        self.body
            .states()
            .iter()
            .zip(projection)
            .map(|(s, p)| s.displacement_m_sqrt_kg * p)
            .sum()
    }

    fn step(&mut self, forces: &[f64]) -> Result<(), String> {
        {
            let states = self.body.states();
            for (k, a) in self.bar_acceleration.iter_mut().enumerate() {
                let (w, s) = (self.bar_omega[k], &states[k]);
                *a = self.bar_held[k]
                    - 2.0 * self.bar_decay[k] * s.velocity_m_sqrt_kg_per_s
                    - w * w * s.displacement_m_sqrt_kg;
            }
            for (a, o) in self.tube_acceleration.iter_mut().zip(&self.tube) {
                *a = o.acceleration();
            }
            self.total.copy_from_slice(forces);
            for m in 0..self.tube.len() {
                let mut on_tube = 0.0;
                let (r, x) = (&self.resistance[m], &self.inertance[m]);
                let (tube_v, tube_a) = (self.tube[m].velocity, self.tube_acceleration[m]);
                for k in 0..self.bar_omega.len() {
                    on_tube -=
                        r[k] * states[k].velocity_m_sqrt_kg_per_s + x[k] * self.bar_acceleration[k];
                    self.total[k] -= r[k] * tube_v + x[k] * tube_a;
                }
                self.tube_force[m] = on_tube;
            }
        }
        let frame = self
            .body
            .step_into(&self.total, &mut self.workspace)
            .map_err(|e| format!("modal step: {e:?}"))?;
        self.bar_pressure = frame.observer_pressure_pa;
        if let Some(o) = self.observer.as_mut() {
            // The same held forces the bar just took: the contact and the tube's reaction on each mode.
            self.observed = o
                .step(&self.total)
                .map_err(|e| format!("broadband observer step: {e:?}"))?;
        }
        let n = self.bar_held.len();
        self.bar_held.copy_from_slice(&self.total[..n]);
        for (o, f) in self.tube.iter_mut().zip(&self.tube_force) {
            o.step(*f);
        }
        Ok(())
    }

    fn sample(&mut self) -> f64 {
        let volume_acceleration: f64 = self
            .tube
            .iter()
            .map(|o| o.acceleration() * o.inverse_root_mass)
            .sum();
        self.history.push(volume_acceleration);
        let mouth = self.monopole * self.delayed_acceleration();
        let narrow = self.bar_pressure + mouth;
        if self.observer.is_none() {
            return narrow;
        }
        self.narrow.push(narrow);
        self.observed + mouth
    }

    fn take_narrow(&mut self) -> Vec<f64> {
        std::mem::take(&mut self.narrow)
    }
}

/// Mouth to receiver as a free compact source: rho / (4 pi r), and the delay in output samples.
pub fn monopole_to_receiver(card: &Card, tube: &Tube, gas: &GasState) -> (f64, f64) {
    let mouth = mouth_centre(card, tube);
    let r = (0..3)
        .map(|c| (card.receiver_m[c] - mouth[c]).powi(2))
        .sum::<f64>()
        .sqrt();
    (
        gas.density / (2.0 * TAU * r),
        r / gas.sound_speed * f64::from(card.rate_hz),
    )
}

/// One JSON row: the isolated column, then what the bar's reflection and field make of it.
pub fn report(
    card: &Card,
    tube: &Tube,
    modes: &[IsolatedMode],
    fields: &[MouthField],
    network: &Network,
) {
    let hz: Vec<f64> = modes.iter().map(|m| m.omega_rad_s / TAU).collect();
    let coupled: Vec<f64> = network.tube.iter().map(|t| t[0] / TAU).collect();
    let decay: Vec<f64> = network.tube.iter().map(|t| t[1]).collect();
    let wall: Vec<f64> = modes.iter().map(|m| m.wall_per_s).collect();
    let open: Vec<f64> = modes.iter().map(|m| m.radiation_per_s).collect();
    let ka: Vec<f64> = modes.iter().map(|m| m.ka).collect();
    let reflection: Vec<[f64; 2]> = fields
        .iter()
        .map(|f| [f.reflection.re, f.reflection.im])
        .collect();
    // Only a drawn bore adds to the row, so a card without one writes the row it always wrote.
    let bore = tube.bore.as_ref().map_or(String::new(), |b| {
        format!(
            ",\"bore_m\":{b:?},\"bore_cylinders\":{},\"acoustic_length_m\":{:?}",
            cylinders(tube).len(),
            modes
                .iter()
                .map(|m| m.acoustic_length_m)
                .collect::<Vec<_>>()
        )
    });
    println!(
        "{{\"suite\":\"bar-strike\",\"case\":\"tube\",\"bar\":\"{}\",\"radius_m\":{},\"air_length_m\":{},\"gap_m\":{},\"mouth_m\":{:?},\"isolated_hz\":{hz:?},\"with_bar_reflection_hz\":{coupled:?},\"decay_per_s\":{decay:?},\"wall_per_s\":{wall:?},\"open_end_per_s\":{open:?},\"ka\":{ka:?},\"reflection_re_im\":{reflection:?},\"resistance\":{:?},\"inertance\":{:?}{bore}}}",
        card.name,
        tube.radius_m,
        tube.air_length_m,
        tube.gap_m,
        mouth_centre(card, tube),
        network.resistance,
        network.inertance,
    );
}
