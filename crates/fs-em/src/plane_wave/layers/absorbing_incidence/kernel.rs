//! Stable real-reference E/H composition with a complex physical entrance.
use super::super::{
    MAX_LAYERS, Slab, admittance_state_from_normal, bad, divided, maximum, propagate_state,
};
use crate::plane_wave::{
    Polarization, ReflectionError, WaveMedium,
    absorbing_incidence::kernel::{
        InterfaceBoundaryField, Medium, entrance_power, exit_power, quadratic, reference_scales,
    },
    finite,
    layers::fields::kernel::LayerFieldKernel,
    normalized_interface_pair, positive,
};
use fs_math::{c64::C64, dd::Dd, det};

type Result<T> = core::result::Result<T, ReflectionError>;

/// Actual peak fields at the entrance and final exit, and cycle-mean normal
/// powers [W/m²]. Entrance interference is signed; film heat is nonnegative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingLayerReflection {
    pub incident_normal_wavenumber_rad_m: C64,
    pub exit_normal_wavenumber_rad_m: C64,
    /// Fixed tangential Ex (TM), Ey (TE); both gains reference their own face.
    pub reflected_tangential_e_gain: [C64; 2],
    pub transmitted_tangential_e_gain: [C64; 2],
    pub incident: InterfaceBoundaryField,
    pub reflected: InterfaceBoundaryField,
    pub transmitted: InterfaceBoundaryField,
    pub incident_normal_power_w_m2: f64,
    /// Signed along +z; outgoing reflected power is nonpositive.
    pub reflected_normal_power_w_m2: f64,
    /// Incident/reflected cross flux at the first face, not absorbed heat.
    pub interference_normal_power_w_m2: f64,
    /// Flux entering the exit half-space, before any subsequent bulk heating.
    pub transmitted_normal_power_w_m2: f64,
    /// Positive loss in the finite films, independently evaluated from their
    /// dissipative modes, never the remainder of exterior power subtraction.
    pub absorbed_power_w_m2: f64,
}

/// The total coherent field inside one original film at the requested depth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingLayerPoint {
    pub field: InterfaceBoundaryField,
    pub normal_power_w_m2: f64,
    /// Original dielectric loss, -omega Im(epsilon) |E|²/2 [W/m³].
    pub heating_w_m3: f64,
}

#[derive(Clone, Copy)]
pub(super) struct Input {
    pub(super) omega: f64,
    pub(super) kx: f64,
    pub(super) epsilon_i: C64,
    pub(super) mu_i: f64,
    pub(super) epsilon_t: C64,
    pub(super) mu_t: f64,
}
#[derive(Clone, Copy)]
pub(super) struct FilmInput {
    pub(super) epsilon: C64,
    pub(super) mu: f64,
    pub(super) thickness: f64,
}
struct Film {
    medium: Medium,
    fields: [LayerFieldKernel; 2],
    // Dimensionless unit-incident real-reference loss modes, TM then TE.
    loss_modes: [[C64; 2]; 2],
}
pub(super) struct Prepared {
    pub(super) input: Input,
    incident: Medium,
    exit: Medium,
    yi: [C64; 2],
    admittance_ref: f64,
    reflected: [C64; 2],
    transmitted_e: [C64; 2],
    transmitted_h: [C64; 2],
    films: Vec<Film>,
}
pub(super) struct Illumination {
    pub(super) response: AbsorbingLayerReflection,
    pub(super) absorbed: [f64; MAX_LAYERS],
}
struct Polarized {
    reflected: C64,
    exit: [C64; 2],
    fields: Vec<LayerFieldKernel>,
    loss_modes: Vec<[C64; 2]>,
}

impl Prepared {
    pub(super) fn new(input: Input, films: &[FilmInput]) -> Result<Self> {
        positive(input.omega)?;
        if !input.kx.is_finite()
            || films.len() > MAX_LAYERS
            || films
                .iter()
                .any(|f| !f.thickness.is_finite() || f.thickness < 0.0)
        {
            return Err(bad(
                "absorbing-incidence stack requires real kx and <=64 finite nonnegative films",
            ));
        }
        let incident = Medium::new(input.omega, input.kx, input.epsilon_i, input.mu_i)?;
        let exit = Medium::new(input.omega, input.kx, input.epsilon_t, input.mu_t)?;
        if incident.tm_impedance == C64::ZERO {
            return Err(bad("stack incident wave must carry positive normal power"));
        }
        let yi = [C64::ONE / incident.tm_impedance, incident.te_admittance];
        if yi.iter().any(|&y| !finite(y) || y.re <= 0.0)
            || exit.te_admittance.re < 0.0
            || exit.tm_impedance.re < 0.0
        {
            return Err(bad(
                "stack incident wave must carry positive power into passive films",
            ));
        }
        let [epsilon_ref, mu_ref, wave_ref, admittance_ref] =
            reference_scales(input.omega, input.epsilon_i, input.mu_i)?;
        let prepare_wave = |epsilon: C64, mu: f64| {
            let wave = WaveMedium::new(epsilon, mu, epsilon_ref, mu_ref)?;
            if (epsilon.re != 0.0
                && (wave.epsilon_ratio.re == 0.0 || wave.index_squared_ratio.re == 0.0))
                || (epsilon.im != 0.0
                    && (wave.epsilon_ratio.im == 0.0 || wave.index_squared_ratio.im == 0.0))
            {
                return Err(bad(
                    "film constitutive contrast underflows the shared reference",
                ));
            }
            Ok(wave)
        };
        let exit_wave = prepare_wave(input.epsilon_t, input.mu_t)?;
        let mut media = Vec::with_capacity(films.len());
        let mut waves = Vec::with_capacity(films.len());
        for film in films {
            media.push(Medium::new(input.omega, input.kx, film.epsilon, film.mu)?);
            waves.push(prepare_wave(film.epsilon, film.mu)?);
        }
        let prepare = |pol, y| {
            prepare_polarization(
                input.kx,
                films,
                &media,
                &waves,
                exit,
                exit_wave,
                wave_ref,
                admittance_ref,
                pol,
                y,
            )
        };
        let tm = prepare(Polarization::Tm, yi[0])?;
        let te = prepare(Polarization::Te, yi[1])?;
        let prepared_films = (0..films.len())
            .map(|i| Film {
                medium: media[i],
                fields: [tm.fields[i], te.fields[i]],
                loss_modes: [tm.loss_modes[i], te.loss_modes[i]],
            })
            .collect();
        Ok(Self {
            input,
            incident,
            exit,
            yi,
            admittance_ref,
            reflected: [tm.reflected, te.reflected],
            transmitted_e: [tm.exit[0], te.exit[0]],
            transmitted_h: [tm.exit[1], te.exit[1]],
            films: prepared_films,
        })
    }

    pub(super) fn illuminate(&self, electric: [C64; 2]) -> Result<Illumination> {
        if electric.iter().any(|&v| !finite(v)) {
            return Err(bad(
                "finite incident peak tangential electric field required",
            ));
        }
        let hi = core::array::from_fn(|c| self.yi[c] * electric[c]);
        let er = core::array::from_fn(|c| self.reflected[c] * electric[c]);
        let hr = core::array::from_fn(|c| -self.yi[c] * er[c]);
        let et = core::array::from_fn(|c| self.transmitted_e[c] * electric[c]);
        let ht = core::array::from_fn(|c| self.transmitted_h[c] * electric[c]);
        let incident = self.incident.field(self.input.kx, electric, hi)?;
        let reflected = self.incident.field(self.input.kx, er, hr)?;
        let transmitted = self.exit.field(self.input.kx, et, ht)?;
        let [pi, pr, interference] = entrance_power(self.yi, self.reflected, electric, er);
        let pt = exit_power(self.exit, et, ht);
        let weight = Dd::from_f64(self.admittance_ref) * Dd::from_f64(0.5);
        let mut absorbed = [0.0; MAX_LAYERS];
        let mut heat = Dd::ZERO;
        for (i, film) in self.films.iter().enumerate() {
            let mut loss = Dd::ZERO;
            for (c, amplitude) in electric.iter().enumerate() {
                for mode in film.loss_modes[c] {
                    loss = loss + quadratic(weight, mode * *amplitude);
                }
            }
            absorbed[i] = loss.to_f64();
            if !absorbed[i].is_finite() || absorbed[i] < 0.0 {
                return Err(bad("finite-film absorbed power is unrepresentable"));
            }
            heat = heat + loss;
        }
        let power = [pi, pr, interference, pt, heat].map(Dd::to_f64);
        if power.iter().any(|p| !p.is_finite())
            || power[0] < 0.0
            || power[1] > 0.0
            || power[3] < 0.0
            || power[4] < 0.0
        {
            return Err(bad("stack exterior power is unrepresentable or nonpassive"));
        }
        if electric.iter().any(|&e| e != C64::ZERO) && power[0] == 0.0 {
            return Err(bad("stack incident power underflows"));
        }
        let scale = power.iter().fold(0.0_f64, |s, p| s.max(p.abs()));
        if scale > 0.0 {
            let scale = Dd::from_f64(scale);
            let defect =
                (pi / scale + pr / scale + interference / scale - pt / scale - heat / scale)
                    .to_f64();
            if !defect.is_finite()
                || defect.abs() > 512.0 * f64::EPSILON * (self.films.len() + 1) as f64
            {
                return Err(bad(
                    "stack Poynting balance with incident interference failed",
                ));
            }
        }
        Ok(Illumination {
            response: AbsorbingLayerReflection {
                incident_normal_wavenumber_rad_m: self.incident.normal,
                exit_normal_wavenumber_rad_m: self.exit.normal,
                reflected_tangential_e_gain: self.reflected,
                transmitted_tangential_e_gain: self.transmitted_e,
                incident,
                reflected,
                transmitted,
                incident_normal_power_w_m2: power[0],
                reflected_normal_power_w_m2: power[1],
                interference_normal_power_w_m2: power[2],
                transmitted_normal_power_w_m2: power[3],
                absorbed_power_w_m2: power[4],
            },
            absorbed,
        })
    }

    pub(super) fn field(
        &self,
        electric: [C64; 2],
        layer: usize,
        depth: f64,
    ) -> Result<AbsorbingLayerPoint> {
        if electric.iter().any(|&e| !finite(e)) {
            return Err(bad(
                "finite incident peak tangential electric field required",
            ));
        }
        let film = self
            .films
            .get(layer)
            .ok_or_else(|| bad("field layer is outside the original stack"))?;
        let mut e = [C64::ZERO; 2];
        let mut h = e;
        for c in 0..2 {
            let trace = film.fields[c]
                .field(depth)
                .ok_or_else(|| bad("field depth is outside the film or is unrepresentable"))?;
            e[c] = trace[0] * electric[c];
            h[c] = trace[1] * electric[c];
        }
        let field = film.medium.field(self.input.kx, e, h)?;
        let mut flux = Dd::ZERO;
        for c in 0..2 {
            flux = flux
                + (Dd::from_f64(e[c].re) * Dd::from_f64(h[c].re)
                    + Dd::from_f64(e[c].im) * Dd::from_f64(h[c].im))
                    * Dd::from_f64(0.5);
        }
        let weight = Dd::from_f64(-film.medium.omega_epsilon.im) * Dd::from_f64(0.5);
        let heat = field
            .electric_v_m
            .iter()
            .fold(Dd::ZERO, |sum, &v| sum + quadratic(weight, v))
            .to_f64();
        let normal_power_w_m2 = flux.to_f64();
        if !normal_power_w_m2.is_finite() || !heat.is_finite() || heat < 0.0 {
            return Err(bad("interior power or heat is unrepresentable"));
        }
        Ok(AbsorbingLayerPoint {
            field,
            normal_power_w_m2,
            heating_w_m3: heat,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare_polarization(
    kx: f64,
    films: &[FilmInput],
    media: &[Medium],
    waves: &[WaveMedium],
    exit: Medium,
    exit_wave: WaveMedium,
    wave_ref: f64,
    admittance_ref: f64,
    pol: Polarization,
    yi: C64,
) -> Result<Polarized> {
    let n = films.len();
    let mut states = [[C64::ZERO; 2]; MAX_LAYERS + 1];
    let mut gains = [C64::ZERO; MAX_LAYERS];
    let mut slabs = [Slab::IDENTITY; MAX_LAYERS];
    let exit_normal = divided(exit.normal, wave_ref);
    if (exit.normal.re != 0.0 && exit_normal.re == 0.0)
        || (exit.normal.im != 0.0 && exit_normal.im == 0.0)
    {
        return Err(bad("exit wave number underflows the shared reference"));
    }
    states[n] = admittance_state_from_normal(exit_wave, 1.0, pol, exit_normal)?;
    for i in (0..n).rev() {
        slabs[i] = Slab::new_at_transverse(
            waves[i],
            wave_ref,
            kx,
            media[i].normal,
            pol,
            films[i].thickness,
        )?;
        (states[i], gains[i]) = propagate_state(slabs[i], states[i + 1])?;
    }
    // The reference port is a real numerical E/H coordinate. The entrance
    // alone is matched against its original complex physical admittance.
    let y = divided(yi, admittance_ref);
    if !finite(y) || y.re <= 0.0 || (yi.im != 0.0 && y.im == 0.0) {
        return Err(bad(
            "incident admittance is unrepresentable in the shared reference",
        ));
    }
    let admittance_scale = positive(maximum(y).max(1.0))?;
    let a = divided(y, admittance_scale) * states[0][0];
    let b = divided(states[0][1], admittance_scale);
    let (a, b, denominator, _) = normalized_interface_pair(a, b)?;
    let reflected = (a - b) / denominator;
    // E+H/Y_ref, evaluated from the independent E and H numerators.
    // Their separate source factors avoid reconstructing either as 1 +/- r.
    let input_e = a.scale(2.0) / denominator;
    let input_h = (b.scale(2.0) / denominator) * y;
    let mut incoming = (input_e + input_h).scale(0.5);
    if !finite(reflected) || !finite(incoming) {
        return Err(bad("absorbing stack entrance amplitude is unrepresentable"));
    }
    let trace = |u: C64, state: [C64; 2]| -> Result<[C64; 2]> {
        let sum = state[0] + state[1];
        if sum == C64::ZERO || !finite(sum) {
            return Err(bad("film field trace is unrepresentable"));
        }
        let e = u * (state[0] / sum).scale(2.0);
        let h = (u * (state[1] / sum).scale(2.0)).scale(admittance_ref);
        if !finite(e) || !finite(h) {
            return Err(bad("film dimensional field is unrepresentable"));
        }
        Ok([e, h])
    };
    let mut left = trace(incoming, states[0])?;
    let mut fields = Vec::with_capacity(n);
    let mut loss_modes = Vec::with_capacity(n);
    for i in 0..n {
        let outgoing = incoming * gains[i];
        let right = trace(outgoing, states[i + 1])?;
        let tail_r = (states[i + 1][0] - states[i + 1][1]) / (states[i + 1][0] + states[i + 1][1]);
        let opposite = tail_r * outgoing;
        let modes = [
            (incoming + opposite)
                .scale(det::sqrt(slabs[i].even_loss) * core::f64::consts::FRAC_1_SQRT_2),
            (incoming - opposite)
                .scale(det::sqrt(slabs[i].odd_loss) * core::f64::consts::FRAC_1_SQRT_2),
        ];
        if modes.iter().any(|&m| !finite(m)) {
            return Err(bad("film loss mode is unrepresentable"));
        }
        loss_modes.push(modes);
        let medium = media[i];
        let (a, b) = match pol {
            Polarization::Te => (
                C64::from_re(medium.omega_mu),
                medium.normal * medium.te_admittance,
            ),
            Polarization::Tm => (medium.normal * medium.tm_impedance, medium.omega_epsilon),
        };
        fields.push(
            LayerFieldKernel::new(
                films[i].thickness,
                medium.normal,
                Dd::from_f64(medium.normal.re),
                a,
                b,
                left,
                right,
            )
            .ok_or_else(|| bad("finite-film depth field is unrepresentable"))?,
        );
        incoming = outgoing;
        left = right;
    }
    Ok(Polarized {
        reflected,
        exit: left,
        fields,
        loss_modes,
    })
}

#[cfg(test)]
mod tests;
