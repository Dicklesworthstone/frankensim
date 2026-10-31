//! Ideal-gas equilibrium at fixed volume, with optional fixed energy.
//!
//! TV closes p V = N(T,p) R T around the original element/charge-constrained Gibbs
//! minimizer. The pressure derivative is the same analytic equilibrium
//! compressibility used by the sound-speed owner. UV then closes the actual
//! species internal energy with equilibrium Cv, including composition response.
//! Formation energy enters once. This is an equilibrium endpoint, not a rate,
//! ignition criterion, real-gas closure or phase-stability determination.
//!
//! Thermodynamic basis: Gordon and McBride, NASA RP-1311, Part I (1994),
//! sections 2.3–2.4 and Table 2.2, <https://ntrs.nasa.gov/citations/19950013764>.

use super::*;
use fs_qty::Volume;

/// Work and residual limits for a fixed-temperature, fixed-volume equilibrium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumVolumeSolveV1 {
    /// Each pressure trial uses the original composition solver and tolerances.
    pub equilibrium: EquilibriumSolveV1,
    /// Maximum absolute `(N R T / p - V) / V` at publication.
    pub relative_volume_tolerance: f64,
    /// At most 256 pressure trials for each temperature.
    pub maximum_pressure_iterations: usize,
}

impl Default for EquilibriumVolumeSolveV1 {
    fn default() -> Self {
        Self {
            equilibrium: EquilibriumSolveV1::default(),
            relative_volume_tolerance: 1e-10,
            maximum_pressure_iterations: 64,
        }
    }
}

/// Temperature path and total internal-energy tolerance for UV equilibrium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumInternalEnergySolveV1 {
    /// Every temperature trial closes volume and chemical equilibrium.
    pub volume: EquilibriumVolumeSolveV1,
    /// Source-supported lower absolute temperature [K].
    pub minimum_temperature_k: f64,
    /// Source-supported upper absolute temperature [K].
    pub maximum_temperature_k: f64,
    /// Absolute total internal-energy tolerance [J] for this feed amount.
    pub energy_tolerance_j: f64,
    /// At most 256 temperature trials, excluding bracket endpoints.
    pub maximum_temperature_iterations: usize,
}

/// A physical or numerical refusal, with no partial equilibrium publication.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumVolumeErrorV1 {
    /// Invalid volume, temperature, work budget or tolerance.
    Invalid(&'static str),
    /// The original source or chemical-equilibrium owner refused the state.
    Equilibrium(CompositionEquilibriumErrorV1),
    /// A source path or reacting thermodynamic derivative was not admissible.
    Response(EquilibriumExpansionErrorV1),
    /// Pressure work was exhausted before the ideal-gas volume equation closed.
    VolumeNoConvergence {
        /// Actual number of pressure trials.
        iterations: usize,
        /// Last `(N R T / p - V) / V` residual.
        relative_residual: f64,
    },
    /// The temperature bracket does not contain the requested total energy.
    InternalEnergyNotBracketed {
        /// Equilibrated internal energy at the lower temperature [J].
        lower_energy_j: f64,
        /// Equilibrated internal energy at the upper temperature [J].
        upper_energy_j: f64,
    },
    /// Temperature work was exhausted without internal-energy closure.
    InternalEnergyNoConvergence {
        /// Actual number of interior temperature trials.
        iterations: usize,
        /// Last product-minus-target internal-energy residual [J].
        residual_j: f64,
    },
}

impl fmt::Display for EquilibriumVolumeErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "fixed-volume ideal-gas equilibrium: {self:?}")
    }
}
impl std::error::Error for EquilibriumVolumeErrorV1 {}
impl From<CompositionEquilibriumErrorV1> for EquilibriumVolumeErrorV1 {
    fn from(value: CompositionEquilibriumErrorV1) -> Self {
        Self::Equilibrium(value)
    }
}
impl From<EquilibriumExpansionErrorV1> for EquilibriumVolumeErrorV1 {
    fn from(value: EquilibriumExpansionErrorV1) -> Self {
        Self::Response(value)
    }
}

/// Complete fixed-volume equilibrium with original species sources and response.
#[derive(Clone, Debug, PartialEq)]
pub struct IdealGasIsochoricEquilibriumV1 {
    thermodynamic: EquilibriumThermodynamicStateV1,
    volume: Volume,
    internal_energy: Energy,
    target_internal_energy: Option<Energy>,
    relative_volume_residual: f64,
    pressure_iterations: usize,
}

impl IdealGasIsochoricEquilibriumV1 {
    /// Accepted chemical composition, temperature, pressure and source receipts.
    pub const fn equilibrium(&self) -> &IdealGasEquilibriumSolutionV1 {
        self.thermodynamic.equilibrium()
    }
    /// Consistent equilibrium entropy, Cp, Cv, compressibility and sound speed.
    pub const fn thermodynamic(&self) -> &EquilibriumThermodynamicStateV1 {
        &self.thermodynamic
    }
    /// Caller-declared fixed total volume [m³].
    pub const fn volume(&self) -> Volume {
        self.volume
    }
    /// Sum of original species internal energies at the accepted composition [J].
    pub const fn internal_energy(&self) -> Energy {
        self.internal_energy
    }
    /// UV target, absent for a temperature-and-volume solve.
    pub const fn target_internal_energy(&self) -> Option<Energy> {
        self.target_internal_energy
    }
    /// Actual product-minus-target internal energy [J], when UV was requested.
    pub fn energy_residual(&self) -> Option<Energy> {
        self.target_internal_energy
            .map(|target| Energy::new(self.internal_energy.value() - target.value()))
    }
    /// Actual `(N R T / p - V) / V` at the accepted state.
    pub const fn relative_volume_residual(&self) -> f64 {
        self.relative_volume_residual
    }
    /// Total pressure trials across all temperature evaluations.
    pub const fn pressure_iterations(&self) -> usize {
        self.pressure_iterations
    }
}

impl EquilibriumVolumeSolveV1 {
    fn validate(self) -> Result<(), EquilibriumVolumeErrorV1> {
        self.equilibrium.validate()?;
        if !self.relative_volume_tolerance.is_finite()
            || self.relative_volume_tolerance <= 0.0
            || self.relative_volume_tolerance >= 0.01
            || !(1..=256).contains(&self.maximum_pressure_iterations)
        {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "invalid volume tolerance or pressure budget",
            ));
        }
        Ok(())
    }
}
impl EquilibriumInternalEnergySolveV1 {
    fn validate(self) -> Result<(), EquilibriumVolumeErrorV1> {
        self.volume.validate()?;
        if !self.minimum_temperature_k.is_finite()
            || self.minimum_temperature_k <= 0.0
            || !self.maximum_temperature_k.is_finite()
            || self.maximum_temperature_k <= self.minimum_temperature_k
            || !self.energy_tolerance_j.is_finite()
            || self.energy_tolerance_j <= 0.0
            || !(1..=256).contains(&self.maximum_temperature_iterations)
        {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "invalid UV temperature bracket, energy tolerance or budget",
            ));
        }
        Ok(())
    }
}

impl IdealGasEquilibriumV1 {
    /// Original feed internal energy at one temperature, without equilibration.
    /// Sum separate streams if their initial temperatures differ. Formation
    /// energy belongs to the original u=h-RT, with no extra heat-release term.
    pub fn feed_internal_energy(
        &self,
        temperature: Temperature,
        feed: &[(SpeciesId, Amount)],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Energy, CompositionEquilibriumErrorV1> {
        poll(&mut cancelled)?;
        let mut energy = Sum::default();
        for (i, amount) in self.feed_amounts(feed)?.into_iter().enumerate() {
            if amount == 0.0 {
                continue;
            }
            poll(&mut cancelled)?;
            let properties = self.models[i]
                .evaluate_properties(temperature)
                .map_err(|source| self.species_error(i, source))?;
            energy.add(amount * properties.u().value());
        }
        poll(&mut cancelled)?;
        if !energy.value().is_finite() {
            return Err(E::Arithmetic("feed internal energy"));
        }
        Ok(Energy::new(energy.value()))
    }

    /// Equilibrate at fixed temperature, total volume and elemental inventory.
    /// Pressure is solved, not held at its feed value. Reuses one bounded Gibbs
    /// workspace through all pressure trials; only the endpoint gets receipts.
    pub fn equilibrate_tv(
        &self,
        temperature: Temperature,
        volume: Volume,
        feed: &[(SpeciesId, Amount)],
        controls: EquilibriumVolumeSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasIsochoricEquilibriumV1, EquilibriumVolumeErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        let mut work = self.volume_workspace(volume, feed)?;
        let mut counts = Counts::default();
        let residual = work.solve_tv(
            self,
            temperature,
            volume,
            controls,
            &mut counts,
            &mut cancelled,
        )?;
        work.publish_volume(
            self,
            temperature,
            volume,
            controls,
            counts,
            0,
            None,
            residual,
            &mut cancelled,
        )
    }

    /// Adiabatic rigid-volume equilibrium for feed at one initial temperature.
    /// The gas can change molecular count, temperature and pressure; volume,
    /// elements and total internal energy remain constrained.
    pub fn equilibrate_uv(
        &self,
        volume: Volume,
        feed: &[(SpeciesId, Amount)],
        inlet_temperature: Temperature,
        controls: EquilibriumInternalEnergySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasIsochoricEquilibriumV1, EquilibriumVolumeErrorV1> {
        let target = self.feed_internal_energy(inlet_temperature, feed, &mut cancelled)?;
        self.equilibrate_uv_energy(volume, feed, target, controls, cancelled)
    }

    /// Fixed-volume equilibrium for an explicit total internal-energy target [J].
    /// The complete temperature bracket must have contiguous source support;
    /// crossing NASA regions requires their explicit caloric-continuity policy.
    /// Analytic reacting Cv supplies the safeguarded Newton slope. Frozen
    /// Cp-R, feed pressure and constant gamma are not substituted.
    pub fn equilibrate_uv_energy(
        &self,
        volume: Volume,
        feed: &[(SpeciesId, Amount)],
        target: Energy,
        controls: EquilibriumInternalEnergySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasIsochoricEquilibriumV1, EquilibriumVolumeErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !target.value().is_finite() {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "UV target must be finite",
            ));
        }
        let mut work = self.volume_workspace(volume, feed)?;
        work.admit_temperature_range(
            self,
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        )?;
        let mut counts = Counts::default();
        let mut low = controls.minimum_temperature_k;
        let mut high = controls.maximum_temperature_k;
        let mut endpoints = [0.0; 2];
        for (index, t) in [low, high].into_iter().enumerate() {
            let temperature = Temperature::new(t);
            let residual = work.solve_tv(
                self,
                temperature,
                volume,
                controls.volume,
                &mut counts,
                &mut cancelled,
            )?;
            endpoints[index] = work.internal_energy(temperature)?;
            if (endpoints[index] - target.value()).abs() <= controls.energy_tolerance_j {
                return work.publish_volume(
                    self,
                    temperature,
                    volume,
                    controls.volume,
                    counts,
                    0,
                    Some(target),
                    residual,
                    &mut cancelled,
                );
            }
        }
        if !(endpoints[0] < target.value() && target.value() < endpoints[1]) {
            return Err(EquilibriumVolumeErrorV1::InternalEnergyNotBracketed {
                lower_energy_j: endpoints[0],
                upper_energy_j: endpoints[1],
            });
        }
        let mut t = low + 0.5 * (high - low);
        let mut residual = endpoints[1] - target.value();
        let mut steps = 0;
        for step in 1..=controls.maximum_temperature_iterations {
            steps = step;
            let temperature = Temperature::new(t);
            let volume_residual = work.solve_tv(
                self,
                temperature,
                volume,
                controls.volume,
                &mut counts,
                &mut cancelled,
            )?;
            residual = work.internal_energy(temperature)? - target.value();
            if residual.abs() <= controls.energy_tolerance_j {
                return work.publish_volume(
                    self,
                    temperature,
                    volume,
                    controls.volume,
                    counts,
                    step,
                    Some(target),
                    volume_residual,
                    &mut cancelled,
                );
            }
            if residual < 0.0 {
                low = t;
            } else {
                high = t;
            }
            let response = work.response(
                self,
                temperature,
                controls.volume.equilibrium,
                &mut cancelled,
            )?;
            let next = t - residual / response.heat_capacity_v_j_k;
            let margin = 0.01 * (high - low);
            t = if next.is_finite() && next > low + margin && next < high - margin {
                next
            } else {
                low + 0.5 * (high - low)
            };
            if t == low || t == high {
                break;
            }
        }
        Err(EquilibriumVolumeErrorV1::InternalEnergyNoConvergence {
            iterations: steps,
            residual_j: residual,
        })
    }

    fn volume_workspace(
        &self,
        volume: Volume,
        feed: &[(SpeciesId, Amount)],
    ) -> Result<Workspace, EquilibriumVolumeErrorV1> {
        if !volume.value().is_finite() || volume.value() <= 0.0 {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "equilibrium volume must be positive and finite",
            ));
        }
        // Initial pressure is only a numerical guess; solve_tv replaces it
        // with an EOS pressure before the first source/composition evaluation.
        Ok(Workspace::new(
            self,
            self.convention.reference_pressure(),
            feed,
        )?)
    }
}

#[derive(Default)]
struct Counts {
    composition: usize,
    pressure: usize,
}

impl Workspace {
    fn internal_energy(&self, temperature: Temperature) -> Result<f64, E> {
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature.value();
        let mut energy = Sum::default();
        for &i in &self.active {
            energy.add(self.n[i] * (self.h[i] - rt));
        }
        let value = energy.value() * self.scale;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(E::Arithmetic("equilibrium internal energy"))
        }
    }

    fn solve_tv(
        &mut self,
        owner: &IdealGasEquilibriumV1,
        temperature: Temperature,
        volume: Volume,
        controls: EquilibriumVolumeSolveV1,
        counts: &mut Counts,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, EquilibriumVolumeErrorV1> {
        poll(cancelled)?;
        if !temperature.value().is_finite() || temperature.value() <= 0.0 {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "equilibrium temperature must be positive and finite",
            ));
        }
        let base = (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature.value())
            * (self.scale / volume.value());
        // Conserved nuclei bound heavy-particle count. Free electrons add at
        // most the positive ionic charge: N <= sum (1 + max(z_i,0)) n_i.
        // Electrons have zero nuclei; including them in a minimum atom count
        // would divide by zero, while ignoring their pressure underbounds N.
        let mut atoms = Sum::default();
        for &b in &self.inventory[..owner.elemental.row_count()] {
            atoms.add(b);
        }
        let mut minimum_atoms = f64::INFINITY;
        let mut maximum_atoms = 0.0_f64;
        let mut maximum_particles_per_atom = 0.0_f64;
        let mut has_electron = false;
        let mut previous_total = Sum::default();
        for &i in &self.active {
            let mut count = Sum::default();
            for k in 0..owner.elemental.row_count() {
                count.add(owner.coefficients[k * self.n.len() + i]);
            }
            let count = count.value();
            if count > 0.0 {
                minimum_atoms = minimum_atoms.min(count);
                maximum_atoms = maximum_atoms.max(count);
                maximum_particles_per_atom = maximum_particles_per_atom
                    .max((1.0 + self.charges[i].max(0.0)) / count);
            } else {
                has_electron = true;
            }
            previous_total.add(self.n[i]);
        }
        let padding = 4.0 * controls.equilibrium.relative_balance_tolerance + 64.0 * f64::EPSILON;
        let lower = base * (atoms.value() / maximum_atoms) * (1.0 - padding);
        let maximum_particles = if has_electron {
            atoms.value() * maximum_particles_per_atom
        } else {
            atoms.value() / minimum_atoms
        };
        let upper = base * maximum_particles * (1.0 + padding);
        if !lower.is_finite() || lower <= 0.0 || !upper.is_finite() || upper < lower {
            return Err(E::Arithmetic("fixed-volume pressure bounds").into());
        }
        let (mut low, mut high) = (ln(lower), ln(upper));
        let mut p = (base * previous_total.value()).clamp(lower, upper);
        let mut residual = f64::INFINITY;
        let mut steps = 0;
        for step in 1..=controls.maximum_pressure_iterations {
            poll(cancelled)?;
            steps = step;
            counts.pressure += 1;
            self.pressure = Pressure::new(p);
            self.log_pressure_ratio =
                crate::entropy::log_ratio(p, owner.convention.reference_pressure().value());
            counts.composition +=
                self.solve(owner, temperature, controls.equilibrium, cancelled)?;
            let mut total = Sum::default();
            for &i in &self.active {
                total.add(self.n[i]);
            }
            let ideal_pressure = base * total.value();
            if !ideal_pressure.is_finite() || ideal_pressure <= 0.0 {
                return Err(E::Arithmetic("fixed-volume ideal pressure").into());
            }
            let log_residual = crate::entropy::log_ratio(ideal_pressure, p);
            residual = fs_math::det::expm1(log_residual);
            if residual.abs() <= controls.relative_volume_tolerance {
                return Ok(residual);
            }
            let log_p = ln(p);
            if residual > 0.0 {
                low = log_p;
            } else {
                high = log_p;
            }
            self.pressure_sensitivity(cancelled)?;
            let compressibility =
                1.0 - self.admit_sensitivity(owner, None, controls.equilibrium)?;
            if !compressibility.is_finite() || compressibility <= 0.0 {
                return Err(EquilibriumVolumeErrorV1::Invalid(
                    "nonpositive equilibrium compressibility",
                ));
            }
            let next = log_p + log_residual / compressibility;
            let candidate = if next.is_finite() && next > low && next < high {
                exp(next)
            } else {
                exp(low + 0.5 * (high - low))
            };
            if !candidate.is_finite() || candidate <= 0.0 || candidate == p {
                break;
            }
            p = candidate;
        }
        Err(EquilibriumVolumeErrorV1::VolumeNoConvergence {
            iterations: steps,
            relative_residual: residual,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn publish_volume(
        &mut self,
        owner: &IdealGasEquilibriumV1,
        temperature: Temperature,
        volume: Volume,
        controls: EquilibriumVolumeSolveV1,
        counts: Counts,
        temperature_iterations: usize,
        target: Option<Energy>,
        relative_volume_residual: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IdealGasIsochoricEquilibriumV1, EquilibriumVolumeErrorV1> {
        let internal_energy = Energy::new(self.internal_energy(temperature)?);
        let thermodynamic = self.publish_response(
            owner,
            temperature,
            controls.equilibrium,
            counts.composition,
            temperature_iterations,
            None,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(IdealGasIsochoricEquilibriumV1 {
            thermodynamic,
            volume,
            internal_energy,
            target_internal_energy: target,
            relative_volume_residual,
            pressure_iterations: counts.pressure,
        })
    }
}
