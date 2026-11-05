//! Settled harmonic carrier inertia on the original neutral and Coulomb drag.
//!
//! For peak phasors with exp(+i omega t), (K+i omega M)u = b E,
//! M_ss=n_s m_s, b_s=n_s q_s. K retains the original neutral anchors and
//! symmetric pair friction. Independent carriers reduce to
//! sigma_s=n_s q_s^2/[m_s(nu_s+i omega)], without converting collision rates
//! into angular frequencies. See Fitzpatrick, equations 1201--1207:
//! <https://farside.ph.utexas.edu/teaching/em/lectures/node103.html>.
//!
//! The mass-normalized symmetric drag is diagonalized once; frequency queries
//! reuse its real modes. Heat comes from the original friction, and reversible
//! drift kinetic storage is reported separately. This local weak-drift fluid
//! closure freezes populations, temperatures and Markovian collision rates.
//! It supplies no Hall/kinetic/nonlocal response, transient initial condition,
//! dielectric background, wave solve or measured broadband qualification.

use super::{
    BOLTZMANN_J_K, ELEMENTARY_CHARGE_C, IonizedGasConductivity, IonizedGasError,
    IonizedGasJouleHeating, poll, positive,
};
use fs_math::{c64::C64, det};

const MAX_CARRIERS: usize = 128;
const BALANCE_TOLERANCE: f64 = 1e-10;

/// Original source plus one prepared mass-normalized relaxation spectrum.
/// Preparation is bounded by 128 carriers and one dense eigendecomposition.
/// Independent carriers keep the diagonal O(N) path; coupled frequency queries
/// use O(N²) work without refactoring a complex system at each frequency.
#[derive(Debug)]
pub struct IonizedGasHarmonic<'a> {
    source: &'a IonizedGasConductivity,
    mass_density: Vec<f64>,
    mass_sqrt: Vec<f64>,
    forcing: Vec<f64>,
    rates: Vec<f64>,
    modes: Option<Vec<f64>>,
    projected_forcing: Vec<f64>,
    pair_carriers: Vec<[usize; 2]>,
}

/// Actual three-dimensional peak current/velocity phasors and cycle averages.
/// Carrier rows preserve source order. Signed thermal-exchange corrections
/// have the same meaning as in the original DC owner; total friction is heat.
#[derive(Debug)]
pub struct IonizedGasHarmonicResponse<'a> {
    /// Exact original frozen populations, temperatures and collision sources.
    pub source: &'a IonizedGasConductivity,
    /// Strictly positive angular frequency [rad/s]; exactly DC is separate.
    pub omega_rad_s: f64,
    /// Isotropic complex free-carrier conductivity [S/m] under exp(+i omega t).
    /// Its real part dissipates; its negative imaginary part stores drift energy.
    pub conductivity_s_m: C64,
    /// Authored three-dimensional peak electric-field phasor [V/m].
    pub electric_field_v_m: [C64; 3],
    /// Total free-carrier peak current density [A/m²], excluding displacement current.
    pub current_density_a_m2: [C64; 3],
    /// Three-dimensional peak drift velocity [m/s] for each original carrier,
    /// in the same order as source.carriers(), relative to the fixed background.
    pub carrier_velocity_m_s: Vec<[C64; 3]>,
    /// Cycle-mean electron/heavy thermal powers [W/m³]. Includes the original
    /// drift-dependent exchange correction once, without zero-drift thermal exchange.
    pub mean_heating_w_m3: IonizedGasJouleHeating,
    /// Cycle-mean directed kinetic storage [J/m³], excluding equilibrium thermal energy.
    pub mean_drift_kinetic_energy_j_m3: f64,
    /// Cycle-mean electrical work Re(sum(E conj(J)))/2 [W/m³], equal to total heat.
    pub mean_electrical_power_w_m3: f64,
    /// Im(sum(E conj(J)))/2 [var/m³], positive for inductive carrier inertia.
    pub reactive_power_var_m3: f64,
    /// Maximum per-carrier complex momentum residual divided by gross force:
    /// electric plus background, inertial and individual pair-force magnitudes.
    pub force_balance_relative_residual: f64,
    /// Maximum relative mismatch of real work versus heat and reactive power
    /// versus 2*omega*kinetic storage, each normalized by its larger magnitude.
    pub power_balance_relative_residual: f64,
}

impl IonizedGasConductivity {
    /// Retain this exact state, all original rates and every admitted mutual
    /// drag pair. No independent-carrier approximation replaces Coulomb drag.
    /// A cancelled or numerically unrepresentable preparation publishes nothing.
    pub fn prepare_harmonic(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IonizedGasHarmonic<'_>, IonizedGasError> {
        poll(&mut cancelled)?;
        let n = self.carriers.len();
        if !(1..=MAX_CARRIERS).contains(&n) {
            return Err(invalid("harmonic carrier count must lie in 1..=128"));
        }
        let mut mass_density = Vec::with_capacity(n);
        let mut mass_sqrt = Vec::with_capacity(n);
        let mut forcing = Vec::with_capacity(n);
        let mut rates = Vec::with_capacity(n);
        let mut projected_forcing = Vec::with_capacity(n);
        for carrier in &self.carriers {
            poll(&mut cancelled)?;
            let mass = carrier.number_density_per_m3 * carrier.mass_kg;
            let root = det::sqrt(mass);
            let force = carrier.number_density_per_m3
                * (carrier.charge_number as f64 * ELEMENTARY_CHARGE_C);
            let weighted = force / root;
            if !positive(mass) || !positive(root) || !weighted.is_finite() || weighted == 0.0 {
                return Err(invalid(
                    "harmonic carrier inertia or forcing is unrepresentable",
                ));
            }
            mass_density.push(mass);
            mass_sqrt.push(root);
            forcing.push(force);
            rates.push(carrier.relaxation.rate_per_s());
            projected_forcing.push(weighted);
        }
        let mut pair_carriers = Vec::with_capacity(self.coulomb_pairs().len());
        let modes = if self.coulomb_pairs().is_empty() {
            None
        } else {
            let mut matrix = vec![0.0; n * n];
            for i in 0..n {
                matrix[i * n + i] = rates[i];
            }
            for pair in self.coulomb_pairs() {
                poll(&mut cancelled)?;
                let index = |original| {
                    self.state
                        .thermodynamics()
                        .models()
                        .get(original)
                        .and_then(|model| {
                            self.carriers
                                .binary_search_by(|carrier| carrier.species().cmp(model.species()))
                                .ok()
                        })
                        .ok_or(invalid("original Coulomb carrier is absent"))
                };
                let i = index(pair.first_index())?;
                let j = index(pair.second_index())?;
                let k = pair.friction_coefficient_kg_m3_s();
                let off = k / mass_sqrt[i] / mass_sqrt[j];
                if !positive(off) {
                    return Err(invalid(
                        "mass-normalized pair friction is unrepresentable",
                    ));
                }
                matrix[i * n + i] += k / mass_density[i];
                matrix[j * n + j] += k / mass_density[j];
                matrix[i * n + j] = -off;
                matrix[j * n + i] = -off;
                pair_carriers.push([i, j]);
            }
            let scale = (0..n)
                .map(|i| matrix[i * n + i])
                .fold(0.0_f64, f64::max);
            if !positive(scale) {
                return Err(invalid("harmonic relaxation scale is unrepresentable"));
            }
            for value in &mut matrix {
                let scaled = *value / scale;
                if !scaled.is_finite() || (*value != 0.0 && scaled == 0.0) {
                    return Err(invalid("harmonic relaxation matrix loses range"));
                }
                *value = scaled;
            }
            poll(&mut cancelled)?;
            let (eigenvalues, eigenvectors) = fs_la::eigen::jacobi_eigh(&matrix, n);
            poll(&mut cancelled)?;
            rates = eigenvalues.into_iter().map(|x| x * scale).collect();
            if rates.iter().any(|&rate| !positive(rate))
                || eigenvectors.iter().any(|x| !x.is_finite())
            {
                return Err(invalid("harmonic relaxation lost positive definiteness"));
            }
            let original = projected_forcing;
            projected_forcing = (0..n)
                .map(|j| sum((0..n).map(|i| eigenvectors[i * n + j] * original[i])))
                .collect();
            if projected_forcing.iter().any(|x| !x.is_finite()) {
                return Err(invalid("harmonic forcing projection is unrepresentable"));
            }
            Some(eigenvectors)
        };
        poll(&mut cancelled)?;
        Ok(IonizedGasHarmonic {
            source: self,
            mass_density,
            mass_sqrt,
            forcing,
            rates,
            modes,
            projected_forcing,
            pair_carriers,
        })
    }
}

impl<'a> IonizedGasHarmonic<'a> {
    /// Original immutable thermochemical state and complete relaxation sources.
    pub const fn source(&self) -> &'a IonizedGasConductivity {
        self.source
    }

    /// Evaluate a positive angular frequency [rad/s] and arbitrary vector peak
    /// phasor [V/m]. Exactly DC uses the existing static owner, avoiding a false
    /// factor of one half at zero frequency. Peak instantaneous drift (including
    /// elliptical polarization and every pair-relative drift) must remain
    /// within the original source limits. Frequency does not enlarge the source
    /// temperature, composition or collision-model applicability.
    pub fn response(
        &self,
        omega_rad_s: f64,
        electric_field_v_m: [C64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IonizedGasHarmonicResponse<'a>, IonizedGasError> {
        poll(&mut cancelled)?;
        if !positive(omega_rad_s) || electric_field_v_m.iter().any(|&e| !finite(e)) {
            return Err(invalid(
                "harmonic response needs positive frequency and finite peak fields",
            ));
        }
        let peak_field = peak(electric_field_v_m);
        let mean_field_squared = mean_square(electric_field_v_m);
        if !peak_field.is_finite() || !mean_field_squared.is_finite() {
            return Err(invalid("harmonic field magnitude is unrepresentable"));
        }
        let n = self.source.carriers.len();
        let modal: Vec<_> = self
            .projected_forcing
            .iter()
            .zip(&self.rates)
            .map(|(&force, &rate)| C64::from_re(force) / C64::new(rate, omega_rad_s))
            .collect();
        let mut mobility = Vec::with_capacity(n);
        let mut velocities = Vec::with_capacity(n);
        let mut forces = Vec::with_capacity(n);
        let mut gross_forces = Vec::with_capacity(n);
        let mut thermal = [0.0; 3]; // heavy friction, electron friction, correction to electrons
        let mut kinetic = 0.0;
        for (i, carrier) in self.source.carriers.iter().enumerate() {
            poll(&mut cancelled)?;
            let weighted = if let Some(modes) = &self.modes {
                C64::new(
                    sum((0..n).map(|j| modes[i * n + j] * modal[j].re)),
                    sum((0..n).map(|j| modes[i * n + j] * modal[j].im)),
                )
            } else {
                modal[i]
            };
            let mu = C64::new(
                weighted.re / self.mass_sqrt[i],
                weighted.im / self.mass_sqrt[i],
            );
            let limit = det::sqrt(BOLTZMANN_J_K * carrier.temperature_k / carrier.mass_kg)
                * self.source.limits.maximum_drift_thermal_ratio;
            if !finite(mu) || mu.abs() * peak_field > limit {
                return Err(invalid("harmonic carrier drift exceeds its source domain"));
            }
            let velocity = electric_field_v_m.map(|e| mu * e);
            let speed_squared = mean_square(velocity);
            let energy = 0.5 * self.mass_density[i] * speed_squared;
            let drag = mu.scale(self.mass_density[i] * carrier.relaxation.rate_per_s());
            let inertia = C64::new(-mu.im, mu.re).scale(self.mass_density[i] * omega_rad_s);
            kinetic += energy;
            thermal[usize::from(carrier.is_electron)] +=
                self.mass_density[i] * carrier.relaxation.rate_per_s() * speed_squared;
            forces.push(drag + inertia);
            gross_forces.push(self.forcing[i].abs() + drag.abs() + inertia.abs());
            mobility.push(mu);
            velocities.push(velocity);
        }
        for (pair, &[i, j]) in self
            .source
            .coulomb_pairs()
            .iter()
            .zip(&self.pair_carriers)
        {
            poll(&mut cancelled)?;
            let relative = mobility[i] - mobility[j];
            let force = relative.scale(pair.friction_coefficient_kg_m3_s());
            forces[i] = forces[i] + force;
            forces[j] = forces[j] - force;
            gross_forces[i] += force.abs();
            gross_forces[j] += force.abs();
            if relative.abs() * peak_field
                > det::sqrt(pair.maximum_relative_speed_squared_m2_s2())
            {
                return Err(invalid("harmonic relative drift exceeds its Coulomb source domain"));
            }
            let speed_squared = mean_square(electric_field_v_m.map(|e| relative * e));
            if speed_squared > 0.0 {
                let heat = pair
                    .drift_heating(speed_squared)
                    .map_err(IonizedGasError::Coulomb)?;
                thermal[usize::from(self.source.carriers[i].is_electron)] +=
                    heat.first_friction_power_density_w_m3();
                thermal[usize::from(self.source.carriers[j].is_electron)] +=
                    heat.second_friction_power_density_w_m3();
                let correction = heat.thermal_exchange_correction_to_first_w_m3();
                if self.source.carriers[i].is_electron {
                    thermal[2] += correction;
                } else if self.source.carriers[j].is_electron {
                    thermal[2] -= correction;
                }
            }
        }
        let mut force_residual = 0.0_f64;
        for i in 0..n {
            let residual = (forces[i] - C64::from_re(self.forcing[i])).abs() / gross_forces[i];
            if !finite(forces[i])
                || !positive(gross_forces[i])
                || !residual.is_finite()
                || residual > BALANCE_TOLERANCE
            {
                return Err(invalid("original harmonic momentum balance failed"));
            }
            force_residual = force_residual.max(residual);
        }
        let conductivity = C64::new(
            sum((0..n).map(|i| self.forcing[i] * mobility[i].re)),
            sum((0..n).map(|i| self.forcing[i] * mobility[i].im)),
        );
        let current = electric_field_v_m.map(|e| conductivity * e);
        let work = conductivity.re * mean_field_squared;
        let reactive = -conductivity.im * mean_field_squared;
        let heat = IonizedGasJouleHeating {
            electron_power_density_w_m3: thermal[1] + thermal[2],
            heavy_power_density_w_m3: thermal[0] - thermal[2],
            electron_friction_power_density_w_m3: thermal[1],
            heavy_friction_power_density_w_m3: thermal[0],
            thermal_exchange_correction_to_electrons_w_m3: thermal[2],
        };
        if !finite(conductivity)
            || conductivity.re <= 0.0
            || conductivity.im >= 0.0
            || current.iter().any(|&j| !finite(j))
            || velocities.iter().flatten().any(|&v| !finite(v))
            || thermal.iter().any(|x| !x.is_finite())
            || thermal[0] < 0.0
            || thermal[1] < 0.0
            || !kinetic.is_finite()
            || kinetic < 0.0
            || !heat.electron_power_density_w_m3.is_finite()
            || !heat.heavy_power_density_w_m3.is_finite()
            || (peak_field > 0.0 && (work <= 0.0 || kinetic <= 0.0))
        {
            return Err(invalid("harmonic current, storage or heat is unrepresentable"));
        }
        let power_residual = balance(work, heat.total_power_density_w_m3())?
            .max(balance(reactive, 2.0 * omega_rad_s * kinetic)?);
        poll(&mut cancelled)?;
        Ok(IonizedGasHarmonicResponse {
            source: self.source,
            omega_rad_s,
            conductivity_s_m: conductivity,
            electric_field_v_m,
            current_density_a_m2: current,
            carrier_velocity_m_s: velocities,
            mean_heating_w_m3: heat,
            mean_drift_kinetic_energy_j_m3: kinetic,
            mean_electrical_power_w_m3: work,
            reactive_power_var_m3: reactive,
            force_balance_relative_residual: force_residual,
            power_balance_relative_residual: power_residual,
        })
    }
}

fn invalid(message: &'static str) -> IonizedGasError {
    IonizedGasError::Invalid(message)
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn sum(values: impl Iterator<Item = f64>) -> f64 {
    let mut result = 0.0;
    let mut correction = 0.0;
    for value in values {
        let increment = value - correction;
        let next = result + increment;
        correction = (next - result) - increment;
        result = next;
    }
    result
}
fn mean_square(vector: [C64; 3]) -> f64 {
    0.5 * sum(vector.into_iter().map(C64::norm_sq))
}
fn peak(vector: [C64; 3]) -> f64 {
    let scale = vector
        .iter()
        .fold(0.0_f64, |s, v| s.max(v.re.abs()).max(v.im.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    let unit = vector.map(|v| C64::new(v.re / scale, v.im / scale));
    let square = sum(unit.iter().map(|z| z.norm_sq()));
    let phase = C64::new(
        sum(unit.iter().map(|z| (z.re - z.im) * (z.re + z.im))),
        sum(unit.iter().map(|z| 2.0 * z.re * z.im)),
    );
    scale * det::sqrt(0.5 * (square + phase.abs()))
}
fn balance(first: f64, second: f64) -> Result<f64, IonizedGasError> {
    let scale = first.abs().max(second.abs());
    let residual = if scale == 0.0 {
        0.0
    } else {
        (first - second).abs() / scale
    };
    if !first.is_finite()
        || !second.is_finite()
        || first < 0.0
        || second < 0.0
        || !residual.is_finite()
        || residual > BALANCE_TOLERANCE
    {
        return Err(invalid(
            "harmonic electrical work, friction or kinetic storage disagree",
        ));
    }
    Ok(residual)
}
