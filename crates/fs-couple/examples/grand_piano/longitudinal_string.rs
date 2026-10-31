//! Retained longitudinal inertia on the existing string/contact clock.
//!
//! The generic fs-nlmodal owner supplies strain, storage and the discrete
//! gradient. This adapter pulls both transverse slope forces and axial bridge
//! tension through the moving bridge and advances axial coordinates using the
//! bank's existing exact-ZOH transitions.

use fs_nlmodal::{
    KcStringParams,
    longitudinal::{LongitudinalString, MAX_LONGITUDINAL_MODES},
};

use super::super::{ModalAcousticMode, Transition, transitions};
use super::{Channel, Observation, StringMode, StringPort, WaveMaterial};
use fs_math::{c64::C64, det};

const MAX_SLOPES: usize = 2 * (super::super::MAX_STRING_MODES + 1);

/// Axial boundary extension of the speaking termination is `sum_n a_n x_n` over the board's loaded
/// coordinates `x_n`: `a_n` is the axial participation row of the board's modes at the string's
/// attachment (`t . (U_n + Theta_n x r)`, from the board's retained nodal motion; see
/// `BoardGeometry::bridge_attachments`). It is a per-mode vector, not a scalar, and it has no default:
/// a `waves` key without its row refuses before this type is built.
pub(super) struct Prepared {
    /// Axial participation row in the bank's loaded board coordinates (length = board mode count).
    axial: Vec<f64>,
    model: LongitudinalString,
    string: usize,
    second: Option<usize>,
    length: f64,
    ea: f64,
    rest_tension: f64,
    material: WaveMaterial,
    transverse_scales: Vec<f64>,
    axial_betas: Vec<f64>,
    transition: Vec<Transition>,
    pub q: Vec<f64>,
    pub v: Vec<f64>,
    next_q: Vec<f64>,
    next_v: Vec<f64>,
    saved_q: Vec<f64>,
    saved_v: Vec<f64>,
    force: Vec<f64>,
    required: Vec<f64>,
    start_slopes: Vec<f64>,
    end_slopes: Vec<f64>,
    slope_gradient: Vec<f64>,
    pub loss_j: f64,
}

impl Prepared {
    pub fn new(
        channel: &Channel,
        density: f64,
        rate: u32,
        material: WaveMaterial,
        axial: Vec<f64>,
    ) -> Result<Self, String> {
        if axial.is_empty() || axial.iter().any(|a| !a.is_finite()) {
            return Err("longitudinal string needs a finite, nonempty axial participation row".into());
        }
        if material.bandwidth_hz > 0.45 * f64::from(rate) {
            return Err("longitudinal band exceeds the mechanics clock's 0.45*rate guard; supply more substeps or a lower explicit band".into());
        }
        let fundamental =
            0.5 * det::sqrt(channel.material.axial_rigidity_n / density) / channel.length;
        if !fundamental.is_finite() || fundamental <= 0. {
            return Err("longitudinal string wave speed is unrepresentable".into());
        }
        let mut indices = Vec::new();
        for mode in 1..=material.mode_cap + 1 {
            if mode as f64 * fundamental > material.bandwidth_hz {
                break;
            }
            if mode > material.mode_cap {
                return Err("longitudinal mode cap cannot cover the supplied axial band; increase the explicit cap or lower the band".into());
            }
            indices.push(mode);
        }
        let model = LongitudinalString::new(
            KcStringParams {
                length: channel.length,
                tension: channel.rest_tension,
                lin_density: density,
                ea: channel.material.axial_rigidity_n,
            },
            channel.diagonal.len(),
            1 + usize::from(channel.second.is_some()),
            &indices,
        )
        .map_err(|e| e.to_string())?;
        let modes: Vec<_> = model
            .axial_omegas()
            .iter()
            .map(|&omega| ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: material.damping_ratio,
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect();
        let transition = transitions(rate, &modes)?;
        let n = indices.len();
        let slopes = model.slope_count();
        let axial_betas = model
            .axial_modes()
            .iter()
            .zip(model.axial_slope_scales())
            .map(|(&mode, &scale)| {
                let sign = if mode % 2 == 0 { -1.0 } else { 1.0 };
                2.0 * sign / (channel.length * scale)
            })
            .collect();
        Ok(Self {
            axial,
            model,
            string: channel.string,
            second: channel.second,
            length: channel.length,
            ea: channel.material.axial_rigidity_n,
            rest_tension: channel.rest_tension,
            material,
            transverse_scales: channel.diagonal.iter().map(|x| det::sqrt(*x)).collect(),
            axial_betas,
            transition,
            q: vec![0.; n],
            v: vec![0.; n],
            next_q: vec![0.; n],
            next_v: vec![0.; n],
            saved_q: vec![0.; n],
            saved_v: vec![0.; n],
            force: vec![0.; n],
            required: vec![0.; n],
            start_slopes: vec![0.; slopes],
            end_slopes: vec![0.; slopes],
            slope_gradient: vec![0.; slopes],
            loss_j: 0.,
        })
    }

    pub fn mode_count(&self) -> usize {
        self.q.len()
    }
    pub fn coupling_terms(&self) -> usize {
        self.model.slope_count() * (1 + self.q.len())
    }

    // A channel retaining no axial mode still has the quasi-static boundary
    // extension `sum a_n x_n`: `evaluate` books its bridge reaction
    // unconditionally, so the potential must include it too (otherwise that
    // force does work the energy does not see).
    fn bridge_axial(&self, q: &[f64], modes: &[StringMode]) -> f64 {
        self.axial.iter().zip(&q[modes.len()..]).map(|(a, x)| a * x).sum()
    }

    fn fill_relative_axial(&self, b_axial: f64, u: &[f64], out: &mut [f64]) {
        for ((slot, &coord), &beta) in out.iter_mut().zip(u).zip(&self.axial_betas) {
            *slot = coord - beta * b_axial;
        }
    }

    fn fill_slopes(
        string: usize,
        second: Option<usize>,
        length: f64,
        scales: &[f64],
        q: &[f64],
        strings: &[StringPort],
        modes: &[StringMode],
        out: &mut [f64],
    ) {
        for (plane, si) in std::iter::once(string).chain(second).enumerate() {
            let s = &strings[si];
            let row = &mut out[plane * (scales.len() + 1)..(plane + 1) * (scales.len() + 1)];
            let bridge = s
                .bridge
                .iter()
                .zip(&q[modes.len()..])
                .map(|(g, x)| g * x)
                .sum::<f64>();
            row[0] = bridge / length;
            for ((slot, &scale), k) in row[1..].iter_mut().zip(scales).zip(s.modes.clone()) {
                *slot = scale * (q[k] - modes[k].beta * bridge);
            }
        }
    }

    fn slopes(&self, q: &[f64], strings: &[StringPort], modes: &[StringMode], out: &mut [f64]) {
        Self::fill_slopes(
            self.string,
            self.second,
            self.length,
            &self.transverse_scales,
            q,
            strings,
            modes,
            out,
        );
    }

    fn potential_and_bounds(
        &self,
        q: &[f64],
        u: &[f64],
        strings: &[StringPort],
        modes: &[StringMode],
    ) -> Option<(f64, f64, f64)> {
        let mut scratch = [0.; MAX_SLOPES];
        let slopes = &mut scratch[..self.model.slope_count()];
        self.slopes(q, strings, modes, slopes);
        let b_axial = self.bridge_axial(q, modes);
        let eps_bridge = b_axial / self.length;
        let mut u_rel_buf = [0.; MAX_LONGITUDINAL_MODES];
        let u_rel = &mut u_rel_buf[..u.len()];
        self.fill_relative_axial(b_axial, u, u_rel);
        let (s0, spread) = self.model.strain_bounds(slopes, u_rel).ok()?;
        let base = self.model.energy(slopes, u_rel).ok()?;
        let extra_mean = 0.5 * self.ea * self.length * eps_bridge * (2.0 * s0 + eps_bridge);
        let potential = base + extra_mean;
        let mean = s0 + eps_bridge;
        if potential.is_finite() && mean.is_finite() && spread.is_finite() {
            Some((potential, mean, spread))
        } else {
            None
        }
    }

    pub fn energy(
        &self,
        q: &[f64],
        strings: &[StringPort],
        modes: &[StringMode],
        candidate: bool,
    ) -> f64 {
        let (u, v) = if candidate {
            (&self.next_q, &self.next_v)
        } else {
            (&self.q, &self.v)
        };
        // Invalid energy propagates to the original finite/energy gate.
        let potential = self
            .potential_and_bounds(q, u, strings, modes)
            .map_or(f64::NAN, |(p, _, _)| p);
        potential + 0.5 * v.iter().map(|x| x * x).sum::<f64>()
    }

    pub fn observe(
        &self,
        q: &[f64],
        strings: &[StringPort],
        modes: &[StringMode],
        out: &mut Observation,
    ) {
        let (potential, mean, spread) = self
            .potential_and_bounds(q, &self.q, strings, modes)
            .unwrap_or((f64::NAN, f64::NAN, f64::NAN));
        out.additional_strain = mean;
        out.tension_n = self.rest_tension + self.ea * mean;
        out.minimum_tension_n = self.rest_tension + self.ea * (mean - spread);
        out.maximum_tension_n = self.rest_tension + self.ea * (mean + spread);
        out.stretching_energy_j = potential;
        out.longitudinal_modes = self.q.len();
        out.longitudinal_kinetic_energy_j = 0.5 * self.v.iter().map(|x| x * x).sum::<f64>();
    }

    pub fn begin(&mut self) {
        self.force.fill(0.);
    }
    pub fn predict(&mut self) {
        self.loss_j = 0.;
        for k in 0..self.q.len() {
            let t = self.transition[k];
            self.next_q[k] = t.qq * self.q[k] + t.qv * self.v[k] + t.bq * self.force[k];
            self.next_v[k] = t.vq * self.q[k] + t.vv * self.v[k] + t.bv * self.force[k];
            let omega = self.model.axial_omegas()[k];
            self.loss_j += self.force[k] * (self.next_q[k] - self.q[k])
                - 0.5
                    * (self.next_v[k].powi(2) - self.v[k].powi(2)
                        + omega * omega * (self.next_q[k].powi(2) - self.q[k].powi(2)));
        }
    }

    pub fn evaluate(
        &mut self,
        a: &[f64],
        b: &[f64],
        strings: &[StringPort],
        modes: &[StringMode],
        mechanical_force: &mut [f64],
    ) -> Result<(), &'static str> {
        Self::fill_slopes(
            self.string,
            self.second,
            self.length,
            &self.transverse_scales,
            a,
            strings,
            modes,
            &mut self.start_slopes,
        );
        Self::fill_slopes(
            self.string,
            self.second,
            self.length,
            &self.transverse_scales,
            b,
            strings,
            modes,
            &mut self.end_slopes,
        );
        let b_axial_a = self.bridge_axial(a, modes);
        let b_axial_b = self.bridge_axial(b, modes);
        let n = self.q.len();
        let mut u_rel_a_buf = [0.; MAX_LONGITUDINAL_MODES];
        let mut u_rel_b_buf = [0.; MAX_LONGITUDINAL_MODES];
        let u_rel_a = &mut u_rel_a_buf[..n];
        let u_rel_b = &mut u_rel_b_buf[..n];
        self.fill_relative_axial(b_axial_a, &self.q, u_rel_a);
        self.fill_relative_axial(b_axial_b, &self.next_q, u_rel_b);
        self.model
            .discrete_gradient(
                &self.start_slopes,
                &self.end_slopes,
                u_rel_a,
                u_rel_b,
                &mut self.slope_gradient,
                &mut self.required,
            )
            .map_err(|_| "longitudinal string discrete gradient failed finite admission")?;
        let s0_a = self
            .model
            .strain_bounds(&self.start_slopes, u_rel_a)
            .map_err(|_| "longitudinal string strain is unrepresentable")?
            .0;
        let s0_b = self
            .model
            .strain_bounds(&self.end_slopes, u_rel_b)
            .map_err(|_| "longitudinal string strain is unrepresentable")?
            .0;
        let s0_mid = f64::midpoint(s0_a, s0_b);
        let eps_bridge_mid = f64::midpoint(b_axial_a, b_axial_b) / self.length;
        let extra_mean_factor = self.ea * self.length * eps_bridge_mid;
        let width = self.transverse_scales.len() + 1;
        for ((sa, sb), out) in self
            .start_slopes
            .chunks_exact(width)
            .zip(self.end_slopes.chunks_exact(width))
            .zip(self.slope_gradient.chunks_exact_mut(width))
        {
            out[0] += extra_mean_factor * f64::midpoint(sa[0], sb[0]);
            for k in 1..width {
                out[k] += 0.5 * extra_mean_factor * f64::midpoint(sa[k], sb[k]);
            }
        }
        // Force on the termination along the string axis; it reaches the board through the axial row.
        let mut axial_bridge_force = self.ea * (s0_mid + eps_bridge_mid);
        for (j, ((value, &omega), &beta)) in self
            .required
            .iter_mut()
            .zip(self.model.axial_omegas())
            .zip(&self.axial_betas)
            .enumerate()
        {
            let g_axial = *value;
            axial_bridge_force -= beta * g_axial;
            let u_mid = f64::midpoint(self.q[j], self.next_q[j]);
            *value = -(g_axial - omega * omega * u_mid);
        }
        for (plane, si) in std::iter::once(self.string).chain(self.second).enumerate() {
            let s = &strings[si];
            let row = &self.slope_gradient[plane * width..(plane + 1) * width];
            let mut reaction = row[0] / self.length;
            for ((&gradient, &scale), k) in row[1..]
                .iter()
                .zip(&self.transverse_scales)
                .zip(s.modes.clone())
            {
                let gradient = gradient * scale;
                mechanical_force[k] -= gradient;
                reaction -= modes[k].beta * gradient;
            }
            for (&shape, value) in s.bridge.iter().zip(&mut mechanical_force[modes.len()..]) {
                *value -= shape * reaction;
            }
            // One physical string, one axial termination: its reaction is booked once, on the first plane.
            if plane == 0 {
                for (&a, value) in self.axial.iter().zip(&mut mechanical_force[modes.len()..]) {
                    *value -= a * axial_bridge_force;
                }
            }
        }
        Ok(())
    }

    pub fn validate(
        &self,
        q: &[f64],
        strings: &[StringPort],
        modes: &[StringMode],
        candidate: bool,
    ) -> Result<(), &'static str> {
        let (u, v) = if candidate {
            (&self.next_q, &self.next_v)
        } else {
            (&self.q, &self.v)
        };
        if u.iter().chain(v).any(|x| !x.is_finite() || x.abs() > 1e5) {
            return Err("longitudinal string state is outside its finite motion domain");
        }
        let (_, mean, spread) = self
            .potential_and_bounds(q, u, strings, modes)
            .ok_or("longitudinal string strain is unrepresentable")?;
        if mean.abs() + spread > self.material.maximum_strain
            || self.rest_tension + self.ea * (mean - spread) <= 0.
        {
            return Err(
                "longitudinal string exceeds its supplied strain bound or cannot guarantee positive tension across the span",
            );
        }
        Ok(())
    }

    pub fn residual(&self, residual: &mut f64, scale: &mut f64) {
        for (&want, &got) in self.required.iter().zip(&self.force) {
            *residual = residual.max((want - got).abs());
            *scale = scale.max(want.abs());
        }
    }
    pub fn correct(&mut self, relaxation: f64) {
        for (got, &want) in self.force.iter_mut().zip(&self.required) {
            *got += relaxation * (want - *got);
        }
    }
    pub fn commit(&mut self) {
        self.q.copy_from_slice(&self.next_q);
        self.v.copy_from_slice(&self.next_v);
    }
    pub fn checkpoint(&mut self) {
        self.saved_q.copy_from_slice(&self.q);
        self.saved_v.copy_from_slice(&self.v);
    }
    pub fn restore(&mut self) {
        self.q.copy_from_slice(&self.saved_q);
        self.v.copy_from_slice(&self.saved_v);
    }
}
