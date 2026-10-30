//! The strike: a point-mass mallet with a cited power-law contact on the bar's modal body, heard through
//! the per-mode receiver transfers. Mirrors the bar_s1 test's `render_stroke` on card data. The body is
//! either the bar alone ([`DirectBody`]) or the bar over a Helmholtz resonator (`resonator::CavityBody`);
//! both run the same mallet loop.
use crate::bar::{TAU, hz};
use crate::card::Card;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use fs_math::{c64::C64, det};
use fs_modal::ModePair;

/// One `scalar` row of a matdb properties file.
pub fn scalar(text: &str, property: &str) -> Result<f64, String> {
    text.lines()
        .find_map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            (fields.first() == Some(&"scalar") && fields.get(3) == Some(&property))
                .then(|| fields.get(4).and_then(|v| v.parse::<f64>().ok()))
                .flatten()
        })
        .ok_or_else(|| format!("required source property missing: {property}"))
}

/// Amplitude-decay law 1/tau = a0 + a2 f^2 as a modal damping ratio (Chaigne and Doutaut 1997,
/// Eqs. 18-19: a0 = gamma_B / 2, a2 = 2 pi^2 eta); `None` when the card declares no loss.
pub struct Loss {
    a0: f64,
    a2: f64,
}

impl Loss {
    pub fn read(text: &str) -> Result<Loss, String> {
        let pi = core::f64::consts::PI;
        Ok(Loss {
            a0: scalar(text, "fluid_damping_gamma_b")? / 2.0,
            a2: 2.0 * pi * pi * scalar(text, "viscoelastic_constant_eta")?,
        })
    }

    pub fn ratio(&self, f_hz: f64) -> f64 {
        (self.a0 + self.a2 * f_hz * f_hz) / (TAU * f_hz)
    }
}

pub struct Mallet {
    pub mass_kg: f64,
    pub speed_m_s: f64,
    pub k_n_per_m_alpha: f64,
    pub alpha: f64,
}

impl Mallet {
    pub fn force(&self, overlap: f64) -> f64 {
        if overlap > 0.0 {
            self.k_n_per_m_alpha * det::pow(overlap, self.alpha)
        } else {
            0.0
        }
    }
}

pub struct Stroke {
    pub samples: Vec<f64>,
    pub contact_time_s: f64,
    pub peak_force_n: f64,
}

/// The audible reference ceilings except the displacement one: the bar has no supports, so after a hard
/// stroke its six free coordinates drift without bound. The drift does not radiate (rigid transfers are
/// zero) and the contact has ended.
pub fn strike_budget() -> ModalAcousticTimeBudget {
    ModalAcousticTimeBudget {
        maximum_abs_displacement_m_sqrt_kg: 100.0,
        ..ModalAcousticTimeBudget::audible_reference()
    }
}

/// The bar's modal oscillators: every band mode with its loss and receiver transfer, then the rigid
/// free-body coordinates (zero frequency, no direct radiation).
pub fn bar_modes(
    modes: &[ModePair],
    loss: Option<&Loss>,
    transfers: &[C64],
    rigid: usize,
) -> Vec<ModalAcousticMode> {
    let mut out: Vec<_> = modes
        .iter()
        .zip(transfers)
        .map(|(m, p)| ModalAcousticMode {
            angular_frequency_rad_s: TAU * hz(m),
            damping_ratio: loss.map_or(0.0, |l| l.ratio(hz(m))),
            pressure_per_modal_velocity: *p,
        })
        .collect();
    out.extend((0..rigid).map(|_| ModalAcousticMode {
        angular_frequency_rad_s: 0.0,
        damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::ZERO,
    }));
    out
}

/// What the mallet loop needs of a body: its displacement at the strike node and a held-force step.
pub trait Body {
    /// Displacement at the strike node along the strike direction, from the retained modal states.
    fn strike_displacement(&self, projection: &[f64]) -> f64;
    /// One substep with held modal forces (flattened bar-mode order, rigid coordinates last).
    fn step(&mut self, forces: &[f64]) -> Result<(), String>;
    /// Receiver pressure after the substeps just taken; called once per output sample.
    fn sample(&mut self) -> f64;
}

/// The bar alone in free air.
pub struct DirectBody {
    body: ModalAcousticTimeModel,
    workspace: ModalAcousticWorkspace,
    pressure: f64,
}

impl DirectBody {
    pub fn new(card: &Card, modes: Vec<ModalAcousticMode>) -> Result<DirectBody, String> {
        let rate = card.rate_hz * card.substeps;
        let body =
            ModalAcousticTimeModel::try_new_with_free_coordinates(rate, modes, strike_budget())
                .map_err(|e| format!("modal body: {e:?}"))?;
        let workspace = ModalAcousticWorkspace::new(&body);
        Ok(DirectBody {
            body,
            workspace,
            pressure: 0.0,
        })
    }
}

impl Body for DirectBody {
    fn strike_displacement(&self, projection: &[f64]) -> f64 {
        self.body
            .states()
            .iter()
            .zip(projection)
            .map(|(s, p)| s.displacement_m_sqrt_kg * p)
            .sum()
    }

    fn step(&mut self, forces: &[f64]) -> Result<(), String> {
        let frame = self
            .body
            .step_into(forces, &mut self.workspace)
            .map_err(|e| format!("modal step: {e:?}"))?;
        self.pressure = frame.observer_pressure_pa;
        Ok(())
    }

    fn sample(&mut self) -> f64 {
        self.pressure
    }
}

/// Render one stroke. `projection` is each coordinate's (negated) displacement at the strike node:
/// every elastic mode, then the rigid coordinates.
pub fn render(
    card: &Card,
    body: &mut dyn Body,
    projection: &[f64],
    mallet: &Mallet,
) -> Result<Stroke, String> {
    let substeps = card.substeps;
    let h = 1.0 / f64::from(card.rate_hz * substeps);
    let mut forces = vec![0.0; projection.len()];
    let (mut position, mut speed) = (0.0, mallet.speed_m_s);
    let (mut detached, mut contact_time, mut peak_force) = (false, 0.0, 0.0_f64);
    let frames = (card.duration_s * f64::from(card.rate_hz)).round() as usize;
    let mut samples = Vec::with_capacity(frames);
    for sample in 0..frames {
        for substep in 0..substeps {
            let overlap = position - body.strike_displacement(projection);
            if !detached && overlap < 0.0 && speed < 0.0 {
                detached = true;
                contact_time = (sample as f64 * f64::from(substeps) + f64::from(substep)) * h;
            }
            let load = if detached { 0.0 } else { mallet.force(overlap) };
            peak_force = peak_force.max(load);
            for (f, p) in forces.iter_mut().zip(projection) {
                *f = p * load;
            }
            body.step(&forces)
                .map_err(|e| format!("{e} at sample {sample}"))?;
            speed -= h * load / mallet.mass_kg;
            position += h * speed;
        }
        samples.push(body.sample());
    }
    if !(detached && contact_time > 0.0) {
        return Err("the mallet never separated within the render".into());
    }
    Ok(Stroke {
        samples,
        contact_time_s: contact_time,
        peak_force_n: peak_force,
    })
}
