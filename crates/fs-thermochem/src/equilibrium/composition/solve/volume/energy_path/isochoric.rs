//! Repeatable reacting internal-energy calorics for sealed rigid storage.
//!
//! The original TV/UV owner fixes elemental inventory and total volume. Every
//! query restores one numerical anchor, so temperature, pressure and reacting
//! derivatives do not depend on earlier thermal trials or failed queries.
//! This supplies constitutive data, not a second heat or chemistry integrator.

use super::*;

// Original assembly overwrites its linear/trial scratch before using it.
// Retain the actual solved variables and source evaluations needed to restart
// TV/UV and to publish a cached cold state without stale hot caloric values.
struct NumericalAnchor {
    n: Vec<f64>,
    log_n: Vec<f64>,
    lambda: Vec<Dd>,
    g: Vec<f64>,
    h: Vec<f64>,
    cp: Vec<f64>,
    log_total: f64,
    log_pressure_ratio: f64,
    sample: EquilibriumEnergyVolumeSampleV1,
}
impl NumericalAnchor {
    fn new(path: &EquilibriumEnergyVolumePathV1<'_>) -> Result<Self, E> {
        fn copy<T: Copy>(source: &[T], zero: T) -> Result<Vec<T>, E> {
            let mut values = allocated(source.len(), zero)?;
            values.copy_from_slice(source);
            Ok(values)
        }
        let w = &path.work;
        Ok(Self {
            n: copy(&w.n, 0.0)?,
            log_n: copy(&w.log_n, 0.0)?,
            lambda: copy(&w.lambda, Dd::ZERO)?,
            g: copy(&w.g, 0.0)?,
            h: copy(&w.h, 0.0)?,
            cp: copy(&w.cp, 0.0)?,
            log_total: w.log_total,
            log_pressure_ratio: w.log_pressure_ratio,
            sample: path
                .current
                .ok_or(E::Invalid("missing original isochoric reference"))?,
        })
    }

    fn restore(&self, path: &mut EquilibriumEnergyVolumePathV1<'_>) {
        let w = &mut path.work;
        w.n.copy_from_slice(&self.n);
        w.log_n.copy_from_slice(&self.log_n);
        w.lambda.copy_from_slice(&self.lambda);
        w.g.copy_from_slice(&self.g);
        w.h.copy_from_slice(&self.h);
        w.cp.copy_from_slice(&self.cp);
        w.log_total = self.log_total;
        w.pressure = self.sample.pressure;
        w.log_pressure_ratio = self.log_pressure_ratio;
        path.current = Some(self.sample);
    }
}

/// Direct original-source changes between two fixed-volume storage endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumIsochoricCaloricChangeV1 {
    pub initial: EquilibriumEnergyVolumeSampleV1,
    pub final_state: EquilibriumEnergyVolumeSampleV1,
    /// Actual requested difference in stored specific energy times feed mass [J].
    pub requested_internal_energy_change_j: f64,
    /// Direct local original-source U_final-U_initial [J].
    pub internal_energy_change_j: f64,
    /// Direct local original-source S_final-S_initial, including mixing [J/K].
    pub entropy_change_j_k: f64,
}

/// Canonical reacting u(T) and T(u) at fixed conserved density.
///
/// The coordinate is (U-U_cold)/original feed mass [J/kg]. It applies to sealed
/// cells with these elemental ratios and this same mass/volume ratio. Pressure
/// and molecular count evolve; reacting Cv is not obtained from frozen Cp-R.
/// The closed-storage datum does not identify a compatible flow enthalpy datum.
pub struct EquilibriumIsochoricCaloricV1<'a> {
    path: EquilibriumEnergyVolumePathV1<'a>,
    anchor: NumericalAnchor,
    density_kg_m3: f64,
    upper_specific_energy_j_kg: f64,
}
impl<'a> EquilibriumIsochoricCaloricV1<'a> {
    /// Admit the complete continuous source bracket, hot TV equilibrium and
    /// exact inclusive upper specific-energy endpoint through the original UV
    /// inverse. Unsupported coordinate rounding refuses without widening U
    /// tolerance. Original composition/pressure/temperature budgets remain.
    pub fn new(
        owner: &'a IdealGasEquilibriumV1,
        volume: Volume,
        feed: &[(SpeciesId, Amount)],
        controls: EquilibriumInternalEnergySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, EquilibriumVolumeErrorV1> {
        let mut path = EquilibriumEnergyVolumePathV1::new(
            owner,
            Temperature::new(controls.minimum_temperature_k),
            volume,
            feed,
            controls,
            &mut cancelled,
        )?;
        let density_kg_m3 = path.feed_mass_kg() / volume.value();
        if !density_kg_m3.is_finite() || density_kg_m3 <= 0.0 {
            return Err(E::Arithmetic("conserved isochoric caloric density").into());
        }
        let anchor = NumericalAnchor::new(&path)?;
        let hot = path.sample_temperature(
            Temperature::new(controls.maximum_temperature_k),
            volume,
            &mut cancelled,
        )?;
        let upper_specific_energy_j_kg = hot.internal_energy_change_j / path.feed_mass_kg();
        let upper_target = upper_specific_energy_j_kg * path.feed_mass_kg();
        if !upper_specific_energy_j_kg.is_finite()
            || upper_specific_energy_j_kg <= 0.0
            || !upper_target.is_finite()
            || upper_target <= 0.0
        {
            return Err(E::Arithmetic("positive representable isochoric energy span").into());
        }
        anchor.restore(&mut path);
        path.sample(volume, Energy::new(upper_target), &mut cancelled)?;
        anchor.restore(&mut path);
        poll(&mut cancelled)?;
        Ok(Self {
            path,
            anchor,
            density_kg_m3,
            upper_specific_energy_j_kg,
        })
    }

    /// Original cold equilibrium, native residuals and complete source receipts.
    pub const fn reference_state(&self) -> &IdealGasIsochoricEquilibriumV1 {
        self.path.reference_state()
    }
    pub fn volume(&self) -> Volume {
        self.anchor.sample.volume
    }
    pub fn feed_mass_kg(&self) -> f64 {
        self.path.feed_mass_kg()
    }
    /// Original feed mass divided by the exact declared sealed volume.
    /// The native response density separately retains its admitted residuals.
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        [
            Temperature::new(self.path.controls.minimum_temperature_k),
            Temperature::new(self.path.controls.maximum_temperature_k),
        ]
    }
    /// Inclusive cold-relative specific internal-energy interval [J/kg].
    pub fn specific_internal_energy_domain_j_kg(&self) -> [f64; 2] {
        [0.0, self.upper_specific_energy_j_kg]
    }

    /// Original UV inverse from the same numerical anchor on every call.
    /// Any refusal invalidates publication but leaves that anchor unchanged.
    pub fn sample_specific_internal_energy(
        &mut self,
        u_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumEnergyVolumeSampleV1, EquilibriumVolumeErrorV1> {
        self.path.current = None;
        poll(&mut cancelled)?;
        if !u_j_kg.is_finite() || !(0.0..=self.upper_specific_energy_j_kg).contains(&u_j_kg) {
            return Err(EquilibriumVolumeErrorV1::Invalid(
                "specific internal energy outside the reacting caloric domain",
            ));
        }
        let target = u_j_kg * self.feed_mass_kg();
        if !target.is_finite() || (u_j_kg != 0.0 && target == 0.0) {
            return Err(E::Arithmetic("reacting specific internal-energy target").into());
        }
        self.anchor.restore(&mut self.path);
        self.path
            .sample(self.volume(), Energy::new(target), &mut cancelled)
    }

    /// Direct original TV evaluation; pressure follows chemistry and volume.
    pub fn sample_temperature(
        &mut self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumEnergyVolumeSampleV1, EquilibriumVolumeErrorV1> {
        self.path.current = None;
        poll(&mut cancelled)?;
        self.anchor.restore(&mut self.path);
        self.path
            .sample_temperature(temperature, self.volume(), &mut cancelled)
    }

    /// Evaluate local source differences before subtracting large cold-relative
    /// coordinates. Actual and requested energy changes remain distinct because
    /// the endpoints retain their original absolute U and volume tolerances.
    pub fn change_between_specific_internal_energies(
        &mut self,
        initial_u_j_kg: f64,
        final_u_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumIsochoricCaloricChangeV1, EquilibriumVolumeErrorV1> {
        let result = (|| {
            let initial = self.sample_specific_internal_energy(initial_u_j_kg, &mut cancelled)?;
            if initial_u_j_kg == final_u_j_kg {
                poll(&mut cancelled)?;
                return Ok(EquilibriumIsochoricCaloricChangeV1 {
                    initial,
                    final_state: initial,
                    requested_internal_energy_change_j: 0.0,
                    internal_energy_change_j: 0.0,
                    entropy_change_j_k: 0.0,
                });
            }
            let local = Reference::new(
                self.path.owner,
                &self.path.work,
                initial.temperature,
                &mut cancelled,
            )?;
            let final_state = self.sample_specific_internal_energy(final_u_j_kg, &mut cancelled)?;
            let internal_energy_change_j = local.energy_change(
                self.path.owner,
                &self.path.work,
                final_state.temperature,
                &mut cancelled,
            )?;
            let entropy_change_j_k = local.entropy_change(
                self.path.owner,
                &self.path.work,
                final_state.temperature,
                &mut cancelled,
            )?;
            let requested_internal_energy_change_j = ((Dd::from_f64(final_u_j_kg)
                - Dd::from_f64(initial_u_j_kg))
                * Dd::from_f64(self.feed_mass_kg()))
            .to_f64();
            if !requested_internal_energy_change_j.is_finite()
                || requested_internal_energy_change_j == 0.0
            {
                return Err(E::Arithmetic("reacting isochoric interval energy").into());
            }
            poll(&mut cancelled)?;
            Ok(EquilibriumIsochoricCaloricChangeV1 {
                initial,
                final_state,
                requested_internal_energy_change_j,
                internal_energy_change_j,
                entropy_change_j_k,
            })
        })();
        if result.is_err() {
            self.path.current = None;
        }
        result
    }

    /// Publish only the last complete query using the original source owner.
    pub fn publish_current(
        &mut self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasIsochoricEquilibriumV1, EquilibriumVolumeErrorV1> {
        self.path.publish_current(cancelled)
    }
}
