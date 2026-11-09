//! Fixed-pressure calorics on one original reacting Gibbs/HP workspace.
//!
//! The cold equilibrium defines a source-relative specific enthalpy coordinate.
//! Every query restores the same bounded numerical anchor before using the
//! original inverse or forward law. Constitutive values therefore do not depend
//! on earlier Newton trials, failed queries, or caller query order. No sampled
//! caloric curve, frozen-composition Cp or separate heat integrator is introduced.

use super::*;

// Only state needed to restart the original composition equations. Linear
// matrices, trial arrays and any allocated full-KKT scratch are reused; each
// original assembly overwrites them before reading. Original source models,
// constraints and conserved inventory are never duplicated or modified.
struct NumericalAnchor {
    n: Vec<f64>,
    log_n: Vec<f64>,
    lambda: Vec<Dd>,
    g: Vec<f64>,
    h: Vec<f64>,
    cp: Vec<f64>,
    log_total: f64,
    log_pressure_ratio: f64,
    sample: EquilibriumEnthalpyPressureSampleV1,
}
impl NumericalAnchor {
    fn new(path: &EquilibriumEnthalpyPressurePathV1<'_>) -> Result<Self, E> {
        fn copy<T: Copy>(source: &[T], zero: T) -> Result<Vec<T>, E> {
            let mut values = allocated(source.len(), zero)?;
            values.copy_from_slice(source);
            Ok(values)
        }
        let w = &path.work;
        Ok(Self {
            n: copy(&w.n, 0.0)?, log_n: copy(&w.log_n, 0.0)?,
            lambda: copy(&w.lambda, Dd::ZERO)?,
            g: copy(&w.g, 0.0)?, h: copy(&w.h, 0.0)?, cp: copy(&w.cp, 0.0)?,
            log_total: w.log_total, log_pressure_ratio: w.log_pressure_ratio,
            sample: path.current.ok_or(E::Invalid("missing original caloric reference"))?,
        })
    }

    fn restore(&self, path: &mut EquilibriumEnthalpyPressurePathV1<'_>) {
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

/// Direct original-source change between two fixed-pressure caloric endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumIsobaricCaloricChangeV1 {
    /// Accepted initial endpoint, on the immutable cold enthalpy reference.
    pub initial: EquilibriumEnthalpyPressureSampleV1,
    /// Accepted final endpoint, on that same reference.
    pub final_state: EquilibriumEnthalpyPressureSampleV1,
    /// Difference of the requested stored specific energies times feed mass [J].
    pub requested_enthalpy_change_j: f64,
    /// Direct source-integrated H_final-H_initial [J].
    pub enthalpy_change_j: f64,
    /// Direct source-integrated S_final-S_initial including redistribution [J/K].
    pub entropy_change_j_k: f64,
}

/// Repeatable reacting caloric law at fixed pressure and elemental inventory.
///
/// The specific coordinate is (H-H_cold)/original feed mass [J/kg], with zero
/// at the admitted lower temperature. The prepared mass is a normalization:
/// the same intensive law applies to any positive mass of these element ratios.
/// This closed caloric datum is not an independently compatible stream datum.
pub struct EquilibriumIsobaricCaloricV1<'a> {
    path: EquilibriumEnthalpyPressurePathV1<'a>,
    anchor: NumericalAnchor,
    upper_specific_enthalpy_j_kg: f64,
}
impl<'a> EquilibriumIsobaricCaloricV1<'a> {
    /// Admit the original source bracket, hot equilibrium and bounded inverse
    /// at the exact advertised inclusive upper specific-enthalpy endpoint.
    ///
    /// All queries retain the original chemistry, source and HP work limits.
    /// This constructor does not assert material transport support or provide
    /// a thermal timestep bound; those remain with the original source adapter.
    pub fn new(
        owner: &'a IdealGasEquilibriumV1,
        pressure: Pressure,
        feed: &[(SpeciesId, Amount)],
        controls: EquilibriumEnthalpySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, EquilibriumExpansionErrorV1> {
        let mut path = EquilibriumEnthalpyPressurePathV1::new(
            owner, Temperature::new(controls.minimum_temperature_k), pressure,
            feed, controls, &mut cancelled,
        )?;
        let anchor = NumericalAnchor::new(&path)?;
        let hot = path.sample_temperature(
            Temperature::new(controls.maximum_temperature_k), pressure, &mut cancelled,
        )?;
        let upper_specific_enthalpy_j_kg = hot.enthalpy_change_j / path.feed_mass_kg();
        if !upper_specific_enthalpy_j_kg.is_finite() || upper_specific_enthalpy_j_kg <= 0.0 {
            return Err(E::Arithmetic("positive reacting caloric enthalpy span").into());
        }
        // Admit the exact extensive target represented by the inclusive specific
        // upper bound. Dividing by mass and multiplying back need not recover
        // the hot source enthalpy within the original absolute HP tolerance.
        // Refuse an unrepresentable endpoint here, without relaxing that gate.
        let upper_target_j = upper_specific_enthalpy_j_kg * path.feed_mass_kg();
        if !upper_target_j.is_finite() || upper_target_j <= 0.0 {
            return Err(E::Arithmetic("representable reacting caloric endpoint target").into());
        }
        anchor.restore(&mut path);
        path.sample(pressure, Energy::new(upper_target_j), &mut cancelled)?;
        anchor.restore(&mut path);
        poll(&mut cancelled)?;
        Ok(Self { path, anchor, upper_specific_enthalpy_j_kg })
    }

    /// Original source-bearing cold equilibrium and its reacting response.
    pub const fn reference_state(&self) -> &EquilibriumThermodynamicStateV1 {
        self.path.reference_state()
    }
    /// Constant pressure [Pa], with its usual pressure-work enthalpy convention.
    pub fn pressure(&self) -> Pressure { self.anchor.sample.pressure }
    /// Original normalization feed mass [kg].
    pub fn feed_mass_kg(&self) -> f64 { self.path.feed_mass_kg() }
    /// Fully admitted original caloric source interval [K].
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        [Temperature::new(self.path.controls.minimum_temperature_k),
            Temperature::new(self.path.controls.maximum_temperature_k)]
    }
    /// Inclusive source-relative specific enthalpy interval [J/kg].
    pub fn specific_enthalpy_domain_j_kg(&self) -> [f64; 2] {
        [0.0, self.upper_specific_enthalpy_j_kg]
    }

    /// Original reacting inverse from the same numerical anchor on every call.
    /// Failed queries invalidate publication without changing that anchor.
    pub fn sample_specific_enthalpy(
        &mut self, h_j_kg: f64, mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumEnthalpyPressureSampleV1, EquilibriumExpansionErrorV1> {
        self.path.current = None;
        poll(&mut cancelled)?;
        if !h_j_kg.is_finite() || !(0.0..=self.upper_specific_enthalpy_j_kg).contains(&h_j_kg) {
            return Err(E::Invalid("specific enthalpy outside the reacting caloric domain").into());
        }
        let target = h_j_kg * self.feed_mass_kg();
        if !target.is_finite() || (h_j_kg != 0.0 && target == 0.0) {
            return Err(E::Arithmetic("reacting specific enthalpy target").into());
        }
        self.anchor.restore(&mut self.path);
        let pressure = self.pressure();
        self.path.sample(pressure, Energy::new(target), &mut cancelled)
    }

    /// Direct forward source enthalpy and reacting Cp at one admitted temperature.
    pub fn sample_temperature(
        &mut self, temperature: Temperature, mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumEnthalpyPressureSampleV1, EquilibriumExpansionErrorV1> {
        self.path.current = None;
        poll(&mut cancelled)?;
        self.anchor.restore(&mut self.path);
        let pressure = self.pressure();
        self.path.sample_temperature(temperature, pressure, &mut cancelled)
    }

    /// Direct source changes between two stored-energy endpoints.
    ///
    /// This uses the initial endpoint as the local source reference before
    /// summing, retaining small changes far from the cold storage anchor.
    /// Actual and requested H changes remain separate because each bounded
    /// inverse retains its original absolute enthalpy tolerance.
    pub fn change_between_specific_enthalpies(
        &mut self, initial_h_j_kg: f64, final_h_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumIsobaricCaloricChangeV1, EquilibriumExpansionErrorV1> {
        let result = (|| {
            let initial = self.sample_specific_enthalpy(initial_h_j_kg, &mut cancelled)?;
            if initial_h_j_kg == final_h_j_kg {
                poll(&mut cancelled)?;
                return Ok(EquilibriumIsobaricCaloricChangeV1 {
                    initial, final_state: initial, requested_enthalpy_change_j: 0.0,
                    enthalpy_change_j: 0.0, entropy_change_j_k: 0.0,
                });
            }
            let local = Reference::new(
                self.path.owner, &self.path.work, initial.temperature, &mut cancelled,
            )?;
            let final_state = self.sample_specific_enthalpy(final_h_j_kg, &mut cancelled)?;
            let enthalpy_change_j = local.enthalpy_change(
                self.path.owner, &self.path.work, final_state.temperature, &mut cancelled,
            )?;
            let entropy_change_j_k = local.entropy_change(
                self.path.owner, &self.path.work, final_state.temperature, &mut cancelled,
            )?;
            let requested_enthalpy_change_j = (final_h_j_kg - initial_h_j_kg) * self.feed_mass_kg();
            if !requested_enthalpy_change_j.is_finite() || requested_enthalpy_change_j == 0.0 {
                return Err(E::Arithmetic("reacting caloric interval energy").into());
            }
            poll(&mut cancelled)?;
            Ok(EquilibriumIsobaricCaloricChangeV1 {
                initial, final_state, requested_enthalpy_change_j,
                enthalpy_change_j, entropy_change_j_k,
            })
        })();
        if result.is_err() { self.path.current = None; }
        result
    }

    /// Publish only the last complete query with its original species receipts.
    pub fn publish_current(
        &mut self, cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumThermodynamicStateV1, EquilibriumExpansionErrorV1> {
        self.path.publish_current(cancelled)
    }
}
