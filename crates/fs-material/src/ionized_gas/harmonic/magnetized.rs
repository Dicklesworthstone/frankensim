//! Prescribed-field gyromotion on the original collisional carrier system.
//!
//! With exp(+i omega t), `(K+i omega M)u + nq B cross u = nq E`.
//! The signed cyclotron frequency is Omega_s=q_s |B|/m_s. Parallel and
//! circular projectors reduce this to three N-carrier systems, with transverse
//! diagonals `K_ss+i M_ss(omega-/+Omega_s)`. Mutual drag stays in K; it is
//! never replaced by independent effective carrier frequencies.
//! See Fitzpatrick, *Cold-Plasma Dielectric Permittivity*, equations 5.18--5.24:
//! <https://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node68.html>.
//! That reference uses exp(-i omega t); the convention above fixes our signs.
//!
//! B does no real work. Its skew term CAN contribute imaginary complex power:
//! `Im(E.conj(J))/2 = 2 omega U + Im(sum(nq (B cross u).conj(u)))/2`.
//! In particular, cyclotron resonance does not erase the drift kinetic storage.
//! Only the original background/pair friction enters the thermal receivers.
//!
//! Populations, temperatures, the uniform prescribed B and isotropic Markovian
//! collision rates stay fixed. This classical local weak-drift fluid model
//! does not provide magnetized collision kinetics, finite-Larmor-radius terms,
//! quantum response, evolving background momentum, self-consistent Maxwell
//! fields, or a scalar conductivity for the existing unmagnetized wave solver.

use super::super::PLANCK_J_S;
use super::{
    BALANCE_TOLERANCE, BOLTZMANN_J_K, C64, ELEMENTARY_CHARGE_C, IonizedGasConductivity,
    IonizedGasError, IonizedGasHarmonic, IonizedGasJouleHeating, balance, det, finite, invalid,
    mean_square, peak, poll, positive, sum,
};
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};

type Result<T> = core::result::Result<T, IonizedGasError>;
type Vector = [C64; 3];

/// Caller-declared classical magnetic domain, in addition to all original
/// population, Coulomb and weak-drift limits. There is no empirical default.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonizedGasMagneticLimits {
    /// Maximum hbar*|q_s B/m_s|/(k*T_s), positive and strictly below one.
    /// A small supplied ceiling excludes Landau-quantized thermal populations;
    /// it is an applicability screen, not an error or validation certificate.
    pub maximum_cyclotron_to_thermal_energy_ratio: f64,
}

/// Fixed-frequency, fixed-B gyrotropic response retaining the original sources.
/// At most 128 carriers are inherited from the harmonic preparation. Independent
/// carriers use O(N) work. With mutual drag, two bounded N-by-N complex solves
/// are prepared once; the existing parallel spectrum is reused. Each subsequent
/// field evaluation costs O(N + number of original charged pairs).
#[derive(Debug)]
pub struct IonizedGasMagnetizedHarmonic<'a> {
    source: &'a IonizedGasConductivity,
    omega_rad_s: f64,
    magnetic_field_t: [f64; 3],
    axis: [f64; 3],
    cyclotron_rad_s: Vec<f64>,
    maximum_cyclotron_to_thermal_energy_ratio: f64,
    mass_density: Vec<f64>,
    forcing: Vec<f64>,
    pair_carriers: Vec<[usize; 2]>,
    // Keep original circular factors: forming Pedersen/Hall coefficients first
    // can lose the weaker handedness near a strongly resolved resonance.
    mobility: Vec<[C64; 3]>,
    conductivity: [C64; 3],
    tensor: [[C64; 3]; 3],
    preparation_relative_residual: f64,
}

/// Actual vector peak phasors and cycle means under exp(+i omega t).
#[derive(Debug)]
pub struct IonizedGasMagnetizedHarmonicResponse<'a> {
    /// Exact original thermochemical and collision sources.
    pub source: &'a IonizedGasConductivity,
    /// Positive drive angular frequency [rad/s].
    pub omega_rad_s: f64,
    /// Uniform prescribed static magnetic flux density [T].
    pub magnetic_field_t: [f64; 3],
    /// Authored Cartesian peak electric field [V/m].
    pub electric_field_v_m: Vector,
    /// Total free-carrier current; excludes dielectric/displacement current.
    pub current_density_a_m2: Vector,
    /// Original carrier order and velocity relative to the fixed background.
    pub carrier_velocity_m_s: Vec<Vector>,
    /// Same source-owned friction and signed thermal-exchange correction as DC.
    pub mean_heating_w_m3: IonizedGasJouleHeating,
    /// Cycle-mean directed kinetic storage, excluding equilibrium thermal energy.
    pub mean_drift_kinetic_energy_j_m3: f64,
    /// Re(sum(E*conj(J)))/2, equal to total physical friction heat.
    pub mean_electrical_power_w_m3: f64,
    /// Im(sum(E*conj(J)))/2, computed per circular mode as
    /// -sum_m Im(sigma_m)*|P_m E|^2/2. Can have either sign under gyromotion.
    pub reactive_power_var_m3: f64,
    /// Signed imaginary power of nq*(B cross u) on the LEFT of momentum balance.
    /// Reactive electrical power = 2*omega*kinetic storage + this term.
    pub gyroscopic_reactive_power_var_m3: f64,
    /// Computed mean work of the actual Lorentz force nq*(u cross B).
    /// Zero to the reported numerical balance tolerance; never deposited as heat.
    pub mean_magnetic_power_w_m3: f64,
    /// Maximum original per-carrier momentum residual normalized by gross forces.
    pub force_balance_relative_residual: f64,
    /// Maximum friction, gyroscopic/reactive and magnetic-no-work residual.
    /// Reactive closure is scaled by its gross terms, retaining resonant cancellation.
    pub power_balance_relative_residual: f64,
}

impl<'a> IonizedGasHarmonic<'a> {
    /// Prepare Hall/cyclotron response for a prescribed static B [T] and a
    /// positive angular frequency [rad/s]. B and subsequent E share one fixed
    /// Cartesian frame. Zero B reproduces the original isotropic response.
    /// Exactly DC remains a different averaging convention and is refused.
    pub fn prepare_magnetized(
        &self,
        omega_rad_s: f64,
        magnetic_field_t: [f64; 3],
        limits: IonizedGasMagneticLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IonizedGasMagnetizedHarmonic<'a>> {
        poll(&mut cancelled)?;
        let ceiling = limits.maximum_cyclotron_to_thermal_energy_ratio;
        if !positive(omega_rad_s)
            || !positive(ceiling)
            || ceiling >= 1.0
            || magnetic_field_t.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "magnetized response needs finite B, positive frequency and a classical-domain limit in (0,1)",
            ));
        }
        let magnitude = det::hypot(
            det::hypot(magnetic_field_t[0], magnetic_field_t[1]),
            magnetic_field_t[2],
        );
        if !magnitude.is_finite() {
            return Err(invalid("magnetic field magnitude is unrepresentable"));
        }
        let axis = if magnitude == 0.0 {
            [0.0, 0.0, 1.0]
        } else {
            magnetic_field_t.map(|b| b / magnitude)
        };
        let mut cyclotron = Vec::with_capacity(self.source.carriers.len());
        let mut quantum_ratio = 0.0_f64;
        for carrier in &self.source.carriers {
            poll(&mut cancelled)?;
            let frequency =
                carrier.charge_number as f64 * ELEMENTARY_CHARGE_C / carrier.mass_kg * magnitude;
            let ratio = (PLANCK_J_S / core::f64::consts::TAU) * frequency.abs()
                / (BOLTZMANN_J_K * carrier.temperature_k);
            if !frequency.is_finite() || !ratio.is_finite() || ratio > ceiling {
                return Err(invalid(
                    "magnetized carrier exceeds the declared classical cyclotron domain",
                ));
            }
            quantum_ratio = quantum_ratio.max(ratio);
            cyclotron.push(frequency);
        }
        let parallel = self.mobility(omega_rad_s, &mut cancelled)?;
        let (plus, minus) = if magnitude == 0.0 {
            (parallel.clone(), parallel.clone())
        } else {
            (
                self.circular_mobility(omega_rad_s, &cyclotron, 1.0, &mut cancelled)?,
                self.circular_mobility(omega_rad_s, &cyclotron, -1.0, &mut cancelled)?,
            )
        };
        let mobility: Vec<_> = (0..parallel.len())
            .map(|i| [parallel[i], plus[i], minus[i]])
            .collect();
        let conductivity = std::array::from_fn(|mode| {
            csum(
                mobility
                    .iter()
                    .zip(&self.forcing)
                    .map(|(u, &b)| u[mode].scale(b)),
            )
        });
        if conductivity.iter().any(|&s| !finite(s) || s.re <= 0.0) {
            return Err(invalid(
                "magnetized conductivity lost representable passivity",
            ));
        }
        let mut tensor = [[C64::ZERO; 3]; 3];
        for column in 0..3 {
            let mut unit = [C64::ZERO; 3];
            unit[column] = C64::ONE;
            let current = combine(conductivity, project(axis, unit));
            for row in 0..3 {
                if !finite(current[row]) {
                    return Err(invalid("magnetized conductivity tensor is unrepresentable"));
                }
                tensor[row][column] = current[row];
            }
        }
        let mut owner = IonizedGasMagnetizedHarmonic {
            source: self.source,
            omega_rad_s,
            magnetic_field_t,
            axis,
            cyclotron_rad_s: cyclotron,
            maximum_cyclotron_to_thermal_energy_ratio: quantum_ratio,
            mass_density: self.mass_density.clone(),
            forcing: self.forcing.clone(),
            pair_carriers: self.pair_carriers.clone(),
            mobility,
            conductivity,
            tensor,
            preparation_relative_residual: 0.0,
        };
        owner.preparation_relative_residual = owner.check_mobility(&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(owner)
    }

    fn circular_mobility(
        &self,
        omega: f64,
        cyclotron: &[f64],
        handedness: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>> {
        let n = self.source.carriers.len();
        let shifts: Vec<_> = cyclotron.iter().map(|&o| omega - handedness * o).collect();
        if shifts.iter().any(|x| !x.is_finite()) {
            return Err(invalid("cyclotron detuning is unrepresentable"));
        }
        if self.pair_carriers.is_empty() {
            let mut result = Vec::with_capacity(n);
            for (i, carrier) in self.source.carriers.iter().enumerate() {
                poll(cancelled)?;
                let charge_mass =
                    carrier.charge_number as f64 * ELEMENTARY_CHARGE_C / carrier.mass_kg;
                let u = C64::from_re(charge_mass)
                    / C64::new(carrier.relaxation.rate_per_s(), shifts[i]);
                if !finite(u) {
                    return Err(invalid("circular carrier mobility is unrepresentable"));
                }
                result.push(u);
            }
            return Ok(result);
        }
        let mut matrix = vec![C64::ZERO; n * n];
        for i in 0..n {
            matrix[i * n + i] =
                C64::new(self.source.carriers[i].relaxation.rate_per_s(), shifts[i]);
        }
        for (pair, &[i, j]) in self.source.coulomb_pairs().iter().zip(&self.pair_carriers) {
            poll(cancelled)?;
            let k = pair.friction_coefficient_kg_m3_s();
            let off = k / self.mass_sqrt[i] / self.mass_sqrt[j];
            matrix[i * n + i].re += k / self.mass_density[i];
            matrix[j * n + j].re += k / self.mass_density[j];
            matrix[i * n + j].re -= off;
            matrix[j * n + i].re -= off;
        }
        // Symmetric diagonal scaling keeps unlike electron/ion inertia and
        // detuning on comparable scales without dropping the mutual forces.
        let scales: Vec<_> = (0..n).map(|i| det::sqrt(matrix[i * n + i].abs())).collect();
        let mut rhs = Vec::with_capacity(n);
        for i in 0..n {
            poll(cancelled)?;
            if !positive(scales[i]) {
                return Err(invalid("magnetized momentum scaling is unrepresentable"));
            }
            for j in 0..n {
                let old = matrix[i * n + j];
                let value = C64::new(
                    old.re / scales[i] / scales[j],
                    old.im / scales[i] / scales[j],
                );
                if !finite(value)
                    || (old.re != 0.0 && value.re == 0.0)
                    || (old.im != 0.0 && value.im == 0.0)
                {
                    return Err(invalid("magnetized momentum matrix loses range"));
                }
                matrix[i * n + j] = value;
            }
            let force = self.forcing[i] / self.mass_sqrt[i] / scales[i];
            if !force.is_finite() || force == 0.0 {
                return Err(invalid("magnetized forcing loses range"));
            }
            rhs.push(C64::from_re(force));
        }
        let factor = lu_complex_with_cancel(&matrix, n, &mut *cancelled).map_err(lu_error)?;
        factor
            .solve_with_cancel(&mut rhs, &mut *cancelled)
            .map_err(lu_error)?;
        for i in 0..n {
            poll(cancelled)?;
            rhs[i] = C64::new(
                rhs[i].re / scales[i] / self.mass_sqrt[i],
                rhs[i].im / scales[i] / self.mass_sqrt[i],
            );
            if !finite(rhs[i]) {
                return Err(invalid("coupled circular velocity is unrepresentable"));
            }
        }
        Ok(rhs)
    }
}

impl<'a> IonizedGasMagnetizedHarmonic<'a> {
    /// Exact immutable population and collision owner retained by preparation.
    pub const fn source(&self) -> &'a IonizedGasConductivity {
        self.source
    }
    /// Positive prepared drive angular frequency [rad/s].
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    /// Original prescribed Cartesian flux density [T].
    pub const fn magnetic_field_t(&self) -> [f64; 3] {
        self.magnetic_field_t
    }
    /// Signed q_s*|B|/m_s [rad/s], in the unchanged source carrier order.
    pub fn cyclotron_frequencies_rad_s(&self) -> &[f64] {
        &self.cyclotron_rad_s
    }
    /// Largest actual hbar*|Omega_s|/(k*T_s), admitted against the supplied ceiling.
    pub const fn maximum_cyclotron_to_thermal_energy_ratio(&self) -> f64 {
        self.maximum_cyclotron_to_thermal_energy_ratio
    }
    /// Row-major J_i = sum_j sigma_ij E_j in the authored Cartesian frame.
    /// Generally nonsymmetric; it must not enter a symmetric-Ohmic field owner.
    pub const fn conductivity_tensor_s_m(&self) -> &[[C64; 3]; 3] {
        &self.tensor
    }
    /// [parallel, e_plus, e_minus] conductivity, with e_+/−=(e1 +/− i e2)/sqrt(2)
    /// for any right-handed (e1,e2,B/|B|). Reversal of B swaps handedness in the
    /// fixed laboratory frame. Positive real parts dissipate physical friction.
    pub const fn circular_conductivities_s_m(&self) -> [C64; 3] {
        self.conductivity
    }

    fn check_mobility(&self, cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
        let n = self.mobility.len();
        let mut maximum = 0.0_f64;
        for (mode, handedness) in [0.0, 1.0, -1.0].into_iter().enumerate() {
            let mut forces = Vec::with_capacity(n);
            let mut gross = Vec::with_capacity(n);
            for i in 0..n {
                poll(cancelled)?;
                let u = self.mobility[i][mode];
                let m = self.mass_density[i];
                let drag = u.scale(m * self.source.carriers[i].relaxation.rate_per_s());
                // Check the detuned equation before separating its two large
                // reactive terms. Their cancellation must not hide a bad solve
                // when collision damping is much smaller than cyclotron rate.
                let inertia = times_i(u)
                    .scale(m * (self.omega_rad_s - handedness * self.cyclotron_rad_s[i]));
                forces.push(drag + inertia - C64::from_re(self.forcing[i]));
                gross.push(drag.abs() + inertia.abs() + self.forcing[i].abs());
            }
            for (pair, &[i, j]) in self.source.coulomb_pairs().iter().zip(&self.pair_carriers) {
                poll(cancelled)?;
                let f = (self.mobility[i][mode] - self.mobility[j][mode])
                    .scale(pair.friction_coefficient_kg_m3_s());
                forces[i] = forces[i] + f;
                forces[j] = forces[j] - f;
                gross[i] += f.abs();
                gross[j] += f.abs();
            }
            for i in 0..n {
                maximum = maximum.max(check_residual(forces[i].abs(), gross[i])?);
            }
        }
        Ok(maximum)
    }

    /// Reuse prepared circular mobilities for an arbitrary peak field phasor.
    /// Actual instantaneous peak carrier AND pair-relative speeds are screened,
    /// including elliptically polarized excitation near cyclotron resonance.
    /// Source failure or cancellation returns no partially accepted response.
    pub fn response(
        &self,
        electric_field_v_m: Vector,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IonizedGasMagnetizedHarmonicResponse<'a>> {
        poll(&mut cancelled)?;
        if electric_field_v_m.iter().any(|&e| !finite(e)) {
            return Err(invalid(
                "magnetized response needs finite peak electric fields",
            ));
        }
        let projection = project(self.axis, electric_field_v_m);
        let n = self.mobility.len();
        let mut velocities = Vec::with_capacity(n);
        let mut forces = Vec::with_capacity(n);
        let mut gross = Vec::with_capacity(n);
        let mut thermal = [0.0; 3];
        let mut kinetic = 0.0;
        let mut gyro_powers = Vec::with_capacity(n);
        let mut gyro_gross = 0.0;
        for (i, carrier) in self.source.carriers.iter().enumerate() {
            poll(&mut cancelled)?;
            let u = combine(self.mobility[i], projection);
            let speed = peak(u);
            let limit = det::sqrt(BOLTZMANN_J_K * carrier.temperature_k / carrier.mass_kg)
                * self.source.limits.maximum_drift_thermal_ratio;
            if u.iter().any(|&v| !finite(v)) || !speed.is_finite() || speed > limit {
                return Err(invalid(
                    "magnetized carrier drift exceeds its source domain",
                ));
            }
            let a = self.mass_density[i] * carrier.relaxation.rate_per_s();
            let drag = u.map(|v| v.scale(a));
            let inertia = u.map(|v| times_i(v).scale(self.mass_density[i] * self.omega_rad_s));
            let gyro = cross(self.magnetic_field_t, u).map(|v| v.scale(self.forcing[i]));
            let electric = electric_field_v_m.map(|e| e.scale(self.forcing[i]));
            forces.push(std::array::from_fn(|j| {
                drag[j] + inertia[j] + gyro[j] - electric[j]
            }));
            gross.push(norm(drag) + norm(inertia) + norm(gyro) + norm(electric));
            let speed2 = mean_square(u);
            kinetic += 0.5 * self.mass_density[i] * speed2;
            thermal[usize::from(carrier.is_electron)] += a * speed2;
            gyro_powers.push(csum((0..3).map(|j| gyro[j] * u[j].conj())).scale(0.5));
            gyro_gross += 0.5 * sum((0..3).map(|j| gyro[j].abs() * u[j].abs()));
            velocities.push(u);
        }
        for (pair, &[i, j]) in self.source.coulomb_pairs().iter().zip(&self.pair_carriers) {
            poll(&mut cancelled)?;
            let relative = std::array::from_fn(|d| velocities[i][d] - velocities[j][d]);
            let speed = peak(relative);
            if !speed.is_finite() || speed > det::sqrt(pair.maximum_relative_speed_squared_m2_s2())
            {
                return Err(invalid(
                    "magnetized relative drift exceeds its Coulomb source domain",
                ));
            }
            let pair_force = relative.map(|v| v.scale(pair.friction_coefficient_kg_m3_s()));
            for d in 0..3 {
                forces[i][d] = forces[i][d] + pair_force[d];
                forces[j][d] = forces[j][d] - pair_force[d];
            }
            gross[i] += norm(pair_force);
            gross[j] += norm(pair_force);
            let speed2 = mean_square(relative);
            if speed2 > 0.0 {
                let heat = pair
                    .drift_heating(speed2)
                    .map_err(IonizedGasError::Coulomb)?;
                thermal[usize::from(self.source.carriers[i].is_electron)] +=
                    heat.first_friction_power_density_w_m3();
                thermal[usize::from(self.source.carriers[j].is_electron)] +=
                    heat.second_friction_power_density_w_m3();
                if self.source.carriers[i].is_electron {
                    thermal[2] += heat.thermal_exchange_correction_to_first_w_m3();
                } else if self.source.carriers[j].is_electron {
                    thermal[2] -= heat.thermal_exchange_correction_to_first_w_m3();
                }
            }
        }
        let mut force_residual = self.preparation_relative_residual;
        for i in 0..n {
            force_residual = force_residual.max(check_residual(norm(forces[i]), gross[i])?);
        }
        let current =
            std::array::from_fn(|d| csum((0..n).map(|i| velocities[i][d].scale(self.forcing[i]))));
        let electrical = csum((0..3).map(|d| electric_field_v_m[d] * current[d].conj())).scale(0.5);
        let gyro = csum(gyro_powers.into_iter());
        let heat = IonizedGasJouleHeating {
            electron_power_density_w_m3: thermal[1] + thermal[2],
            heavy_power_density_w_m3: thermal[0] - thermal[2],
            electron_friction_power_density_w_m3: thermal[1],
            heavy_friction_power_density_w_m3: thermal[0],
            thermal_exchange_correction_to_electrons_w_m3: thermal[2],
        };
        if current.iter().any(|&j| !finite(j))
            || !finite(electrical)
            || !finite(gyro)
            || thermal.iter().any(|v| !v.is_finite())
            || thermal[0] < 0.0
            || thermal[1] < 0.0
            || !kinetic.is_finite()
            || kinetic < 0.0
            || !heat.electron_power_density_w_m3.is_finite()
            || !heat.heavy_power_density_w_m3.is_finite()
            || (electric_field_v_m.iter().any(|&e| e != C64::ZERO)
                && (electrical.re <= 0.0 || kinetic <= 0.0))
        {
            return Err(invalid(
                "magnetized current, heat or storage is unrepresentable",
            ));
        }
        let inertial = 2.0 * self.omega_rad_s * kinetic;
        // Reactive power from the circular modes. The projections P_m E of
        // `project` are Hermitian-orthogonal and J = sum_m sigma_m P_m E, so
        // Im(E.conj(J))/2 = -sum_m Im(sigma_m)*|P_m E|^2/2. Each term is one
        // mode's own conductivity times its own field power. The Cartesian
        // E.conj(J) instead forms each Im(E_d conj(J_d)) as a difference of
        // products as large as |E_d||J_d|; when collisions dominate
        // (|Im sigma| << Re sigma, as at omega = 1e-6 rad/s) the rounding of that
        // difference is larger than 2*omega*U itself.
        let modal_reactive = -0.5
            * sum((0..3)
                .map(|m| self.conductivity[m].im * sum(projection[m].iter().map(|v| v.norm_sq()))));
        let reactive_scale = modal_reactive.abs().max(inertial.abs()).max(gyro.im.abs());
        let reactive_residual = if reactive_scale == 0.0 {
            0.0
        } else {
            (modal_reactive / reactive_scale - inertial / reactive_scale - gyro.im / reactive_scale)
                .abs()
        };
        // Second, independent check of the current built from the carrier
        // velocities: E.conj(J) closes the same balance to the tolerance plus its
        // own rounding floor. Each Im(E_d conj(J_d)) = Im(E_d)Re(J_d) -
        // Re(E_d)Im(J_d) is a difference of two products bounded by |E_d||J_d|;
        // the products, the difference and the compensated sum each round by at
        // most eps of that bound, and J_d carries a few eps of |J_d| from its mode
        // and carrier sums when those sums do not cancel. 8 eps*sum_d |E_d||J_d|
        // bounds the floor of the halved sum with margin. Under cancellation that
        // floor, not the relative tolerance, sets the allowance.
        let electrical_scale = electrical.im.abs().max(inertial.abs()).max(gyro.im.abs());
        let electrical_floor = 8.0
            * f64::EPSILON
            * sum((0..3).map(|d| electric_field_v_m[d].abs() * current[d].abs()));
        let electrical_closes = electrical_scale == 0.0
            || (electrical.im / electrical_scale
                - inertial / electrical_scale
                - gyro.im / electrical_scale)
                .abs()
                <= BALANCE_TOLERANCE + electrical_floor / electrical_scale;
        if !inertial.is_finite()
            || !reactive_residual.is_finite()
            || reactive_residual > BALANCE_TOLERANCE
            || !electrical_closes
        {
            return Err(invalid(
                "magnetized reactive power, gyromotion and kinetic storage disagree",
            ));
        }
        let power_residual = balance(electrical.re, heat.total_power_density_w_m3())?
            .max(reactive_residual)
            .max(check_residual(gyro.re.abs(), gyro_gross)?);
        poll(&mut cancelled)?;
        Ok(IonizedGasMagnetizedHarmonicResponse {
            source: self.source,
            omega_rad_s: self.omega_rad_s,
            magnetic_field_t: self.magnetic_field_t,
            electric_field_v_m,
            current_density_a_m2: current,
            carrier_velocity_m_s: velocities,
            mean_heating_w_m3: heat,
            mean_drift_kinetic_energy_j_m3: kinetic,
            mean_electrical_power_w_m3: electrical.re,
            reactive_power_var_m3: modal_reactive,
            gyroscopic_reactive_power_var_m3: gyro.im,
            mean_magnetic_power_w_m3: -gyro.re,
            force_balance_relative_residual: force_residual,
            power_balance_relative_residual: power_residual,
        })
    }
}

fn lu_error(error: LuComplexError) -> IonizedGasError {
    match error {
        LuComplexError::Cancelled => IonizedGasError::Cancelled,
        LuComplexError::Singular { .. } => {
            invalid("magnetized momentum factorization is numerically singular")
        }
    }
}
fn times_i(z: C64) -> C64 {
    C64::new(-z.im, z.re)
}
fn csum(values: impl Iterator<Item = C64>) -> C64 {
    let mut total = C64::ZERO;
    let mut correction = C64::ZERO;
    for value in values {
        let increment = value - correction;
        let next = total + increment;
        correction = (next - total) - increment;
        total = next;
    }
    total
}
fn cross(a: [f64; 3], b: Vector) -> Vector {
    [
        b[2].scale(a[1]) - b[1].scale(a[2]),
        b[0].scale(a[2]) - b[2].scale(a[0]),
        b[1].scale(a[0]) - b[0].scale(a[1]),
    ]
}
fn project(axis: [f64; 3], field: Vector) -> [Vector; 3] {
    let parallel = csum((0..3).map(|d| field[d].scale(axis[d])));
    let along = axis.map(|b| parallel.scale(b));
    let across = std::array::from_fn(|d| field[d] - along[d]);
    let rotated = cross(axis, across).map(times_i);
    [
        along,
        std::array::from_fn(|d| (across[d] + rotated[d]).scale(0.5)),
        std::array::from_fn(|d| (across[d] - rotated[d]).scale(0.5)),
    ]
}
fn combine(mobility: [C64; 3], projection: [Vector; 3]) -> Vector {
    std::array::from_fn(|d| csum((0..3).map(|m| mobility[m] * projection[m][d])))
}
fn norm(vector: Vector) -> f64 {
    det::hypot(
        det::hypot(vector[0].abs(), vector[1].abs()),
        vector[2].abs(),
    )
}
fn check_residual(residual: f64, scale: f64) -> Result<f64> {
    let relative = if scale == 0.0 && residual == 0.0 {
        0.0
    } else {
        residual / scale
    };
    if !residual.is_finite()
        || !scale.is_finite()
        || !relative.is_finite()
        || relative > BALANCE_TOLERANCE
    {
        return Err(invalid(
            "original magnetized momentum or magnetic no-work balance failed",
        ));
    }
    Ok(relative)
}
