//! A Helmholtz resonator under the bar, from the existing owners: `HelmholtzModalCavity` (fs-phs
//! compliance, neck end corrections and resonance-pinned radiation resistance) coupled to the bar's modal
//! body through `SharedCavity` (one compression coordinate, two-way reactions). The neck's outflow is
//! heard as a compact flanged monopole. Declared omissions, none of them corrected by a constant:
//! the box and the neighbouring bars do not screen the bar's own radiation (the bar is heard in free
//! space); the neck's mouth is a point under the bar's lower face; bars on the same box are absent, so
//! the measured detuning by the bars (Rucz Table IV, 6-15 %) is not reproduced; the six rigid
//! coordinates are not coupled to the cavity (the bar's supports are absent).
use crate::bar::{Mesh, TAU};
use crate::card::{Card, Resonator};
use crate::render::{Body, strike_budget};
use fs_couple::cavity_modal::shared::{CavitySurface, SharedCavity};
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::ModalCouplingConfig;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_modal::ModePair;
use fs_scenario::{HelmholtzCavity, NeckEnd};

/// The card's resonator as the library's lumped cavity: the opening's circle of equal area as the neck,
/// flanged at both ends (the box's plate), loss from the owner's radiation at resonance plus the exact
/// viscothermal neck.
pub fn cavity(res: &Resonator, gas: &GasState) -> Result<HelmholtzModalCavity, String> {
    let radius = (res.opening_area_m2 / core::f64::consts::PI).sqrt();
    HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: res.volume_m3,
            neck_radius_m: radius,
            neck_length_m: res.neck_length_m,
            inner_end: NeckEnd::Flanged,
            outer_end: NeckEnd::Flanged,
        },
        *gas,
        CavityNeckLoss::RadiationAtResonance,
    )
    .map_err(|e| format!("Helmholtz cavity: {e:?}"))
}

/// Isolated resonance [Hz] of the lumped cavity, before any coupling to a bar.
pub fn isolated_hz(cavity: &HelmholtzModalCavity) -> f64 {
    cavity.isolated_omega_rad_s() / TAU
}

/// Outward (away from the cavity, +z) volume projection of each mode over the wet patch of the lower
/// face [m^2/sqrt(kg)]: the lower face's area vectors, weighted by the share of each grid cell inside the
/// wet interval, dotted with the mode's mean vertex displacement on that cell.
pub fn wet_volume_shapes(mesh: &Mesh, modes: &[ModePair], wet_x_m: [f64; 2]) -> Vec<f64> {
    let [nx, ny, nz] = mesh.counts;
    let node = |i: usize, j: usize| (i * (ny + 1) + j) * (nz + 1);
    let mut shapes = vec![0.0; modes.len()];
    for i in 0..nx {
        for j in 0..ny {
            let corners = [
                node(i, j),
                node(i + 1, j),
                node(i + 1, j + 1),
                node(i, j + 1),
            ];
            let [p00, p10, p11, p01] = corners.map(|n| mesh.nodes[n]);
            let (xa, xb) = (p00[0], p10[0]);
            let overlap = (xb.min(wet_x_m[1]) - xa.max(wet_x_m[0])).max(0.0) / (xb - xa);
            if overlap == 0.0 {
                continue;
            }
            // Area vector of the quad from its diagonals; the lower face is a graph over (x, y), so +z.
            let d1: [f64; 3] = core::array::from_fn(|c| p11[c] - p00[c]);
            let d2: [f64; 3] = core::array::from_fn(|c| p01[c] - p10[c]);
            let mut area = [
                0.5 * (d1[1] * d2[2] - d1[2] * d2[1]),
                0.5 * (d1[2] * d2[0] - d1[0] * d2[2]),
                0.5 * (d1[0] * d2[1] - d1[1] * d2[0]),
            ];
            if area[2] < 0.0 {
                area = area.map(|a| -a);
            }
            for (shape, mode) in shapes.iter_mut().zip(modes) {
                let mean: f64 = (0..3)
                    .map(|c| {
                        area[c] * corners.iter().map(|&n| mode.phi[3 * n + c]).sum::<f64>() / 4.0
                    })
                    .sum();
                *shape += overlap * mean;
            }
        }
    }
    shapes
}

/// The bar over its resonator, stepped as one reciprocal modal network.
pub struct CavityBody {
    shared: SharedCavity,
    gate: CancelGate,
    bar_pressure: f64,
    q: f64,
    q_prev: f64,
    substep_s: f64,
    /// rho / (2 pi r): pressure of a flanged compact source per unit volume acceleration.
    monopole: f64,
    delay_samples: f64,
    history: Vec<f64>,
}

impl CavityBody {
    /// `bar` are the bar's modal oscillators (band modes then rigid coordinates); `shapes` the cavity
    /// volume projection of each elastic mode (the rigid coordinates get zero).
    pub fn new(
        card: &Card,
        res: &Resonator,
        gas: &GasState,
        bar: Vec<ModalAcousticMode>,
        mut shapes: Vec<f64>,
    ) -> Result<CavityBody, String> {
        shapes.resize(bar.len(), 0.0);
        let rate = card.rate_hz * card.substeps;
        let model =
            ModalAcousticTimeModel::try_new_with_free_coordinates(rate, bar, strike_budget())
                .map_err(|e| format!("modal body: {e:?}"))?;
        let modes = shapes.len();
        let config = ModalCouplingConfig {
            max_modes: modes + 1,
            max_connections: 2,
            max_setup_terms: 1 << 24,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 10.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8,
        };
        let gate = CancelGate::new_clock_free();
        let shared = SharedCavity::new(
            cavity(res, gas)?,
            vec![CavitySurface {
                model,
                outward_volume_shapes: shapes,
            }],
            rate,
            ModalAcousticState::default(),
            strike_budget(),
            config,
            &gate,
        )
        .map_err(|e| format!("shared cavity: {e:?}"))?;
        // The mouth: a point on the bar's lower face at the middle of the wet interval, mid-width.
        let x = 0.5 * (res.wet_x_m[0] + res.wet_x_m[1]);
        let z = 0.5 * card.thickness_m - card.profile.thickness(card.length_m, card.thickness_m, x);
        let mouth = [x, 0.0, z];
        let r = (0..3)
            .map(|c| (card.receiver_m[c] - mouth[c]).powi(2))
            .sum::<f64>()
            .sqrt();
        Ok(CavityBody {
            shared,
            gate,
            bar_pressure: 0.0,
            q: 0.0,
            q_prev: 0.0,
            substep_s: 1.0 / f64::from(rate),
            monopole: gas.density / (2.0 * core::f64::consts::PI * r),
            delay_samples: r / gas.sound_speed * f64::from(card.rate_hz),
            history: Vec::new(),
        })
    }

    /// Volume acceleration of the neck's outflow, `delay_samples` ago, linearly interpolated.
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

impl Body for CavityBody {
    fn strike_displacement(&self, projection: &[f64]) -> f64 {
        self.shared.network().components()[0]
            .states()
            .iter()
            .zip(projection)
            .map(|(s, p)| s.displacement_m_sqrt_kg * p)
            .sum()
    }

    fn step(&mut self, forces: &[f64]) -> Result<(), String> {
        let frame = self
            .shared
            .step(forces, 0.0, &self.gate)
            .map_err(|e| format!("cavity step: {e:?}"))?;
        self.bar_pressure = frame.observer_pressure_pa;
        let observation = self
            .shared
            .observe()
            .map_err(|e| format!("cavity observation: {e:?}"))?;
        self.q_prev = self.q;
        self.q = observation.neck_outflow_m3_s;
        Ok(())
    }

    fn sample(&mut self) -> f64 {
        self.history.push((self.q - self.q_prev) / self.substep_s);
        self.bar_pressure + self.monopole * self.delayed_acceleration()
    }
}
