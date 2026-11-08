//! Finite, frozen-composition NASA-mixture heat storage coupled to chemistry.
//!
//! Both gases retain their own sources, amount, mechanical constraint and
//! temperature domain. Only the reactor changes composition. Exchange carries
//! heat, not species; original formation references cancel before arithmetic.

use super::*;
use crate::chemistry::isobaric::{IsobaricKineticsView, IsobaricReactor};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::kinetics::admit_caloric_interval;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

mod solve;

/// Independently prescribed mechanical constraint of the finite partner gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GasHeatConstraint {
    /// Rigid storage: heat changes U and pressure, with no boundary work.
    ConstantVolume(Volume),
    /// Held pressure: heat changes H; expansion work is reported separately.
    ConstantPressure(Pressure),
}

/// One original immutable mixture and its finite, uniformly mixed inventory.
/// Composition and declared exchange geometry remain fixed; gray coupling is explicit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasHeatBody<'m> {
    model: &'m FrozenIdealGasMixtureModelV1,
    amount: Amount,
    mass: Mass,
    constraint: GasHeatConstraint,
    bounds: [Temperature; 2],
    bath: HeatBath,
}

impl<'m> FiniteGasHeatBody<'m> {
    /// Admit the whole partner source path before evolution. Its bounds are
    /// independent of the reacting gas's temperature bounds. No source fit,
    /// constant-Cp replacement, equilibrium projection or mass transfer occurs.
    pub fn new(
        model: &'m FrozenIdealGasMixtureModelV1,
        amount: Amount,
        constraint: GasHeatConstraint,
        temperature: Temperature,
        bounds: [Temperature; 2],
        conductance: ThermalConductance,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IsochoricKineticsError> {
        poll(&mut cancelled)?;
        if !positive(amount.value())
            || !positive(bounds[0].value())
            || !positive(bounds[1].value())
            || bounds[0].value() > bounds[1].value()
            || !temperature.value().is_finite()
            || temperature.value() < bounds[0].value()
            || temperature.value() > bounds[1].value()
            || !positive(match constraint {
                GasHeatConstraint::ConstantVolume(v) => v.value(),
                GasHeatConstraint::ConstantPressure(p) => p.value(),
            })
        {
            return Err(IsochoricKineticsError::Invalid(
                "invalid finite gas heat inventory or source bounds",
            ));
        }
        for source in model.components() {
            admit_caloric_interval(source, bounds[0], bounds[1], &mut cancelled)?;
        }
        let mut body = Self {
            model,
            amount,
            mass: Mass::new(0.0),
            constraint,
            bounds,
            bath: HeatBath::new(conductance, temperature)?,
        };
        for t in bounds {
            let (p, v) = body.mechanics(t);
            if !positive(p.value()) || !positive(v.value()) {
                return Err(IsochoricKineticsError::Invalid(
                    "unrepresentable finite gas pressure or volume",
                ));
            }
        }
        let properties = model.evaluate_properties(temperature, body.pressure())?;
        body.mass = Mass::new(amount.value() * properties.molar_mass().value());
        if !positive(body.mass.value()) {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite gas mass",
            ));
        }
        poll(&mut cancelled)?;
        Ok(body)
    }

    pub const fn model(self) -> &'m FrozenIdealGasMixtureModelV1 {
        self.model
    }
    pub const fn amount(self) -> Amount {
        self.amount
    }
    pub const fn mass(self) -> Mass {
        self.mass
    }
    pub const fn constraint(self) -> GasHeatConstraint {
        self.constraint
    }
    pub const fn temperature(self) -> Temperature {
        self.bath.temperature()
    }
    pub const fn temperature_bounds(self) -> [Temperature; 2] {
        self.bounds
    }
    pub const fn conductance(self) -> ThermalConductance {
        self.bath.conductance()
    }
    pub const fn effective_radiating_area(self) -> Area {
        self.bath.effective_radiating_area()
    }
    pub fn pressure(self) -> Pressure {
        self.mechanics(self.temperature()).0
    }
    pub fn volume(self) -> Volume {
        self.mechanics(self.temperature()).1
    }

    pub(super) fn radiation_capacity_density(self, step: f64, volume: f64) -> f64 {
        self.bath.radiation_capacity_density(step, volume)
    }

    /// Retain caller-declared epsilon*A under the existing gray-bath convention.
    /// This does not derive a participating-gas opacity from the species list.
    pub fn with_gray_radiation(
        mut self,
        area: Area,
        effective_emissivity: f64,
    ) -> Result<Self, IsochoricKineticsError> {
        self.bath = self.bath.with_gray_radiation(area, effective_emissivity)?;
        Ok(self)
    }

    fn mechanics(self, temperature: Temperature) -> (Pressure, Volume) {
        match self.constraint {
            GasHeatConstraint::ConstantVolume(v) => (
                Pressure::new((self.amount.value() / v.value()) * R * temperature.value()),
                v,
            ),
            GasHeatConstraint::ConstantPressure(p) => (
                p,
                Volume::new((self.amount.value() / p.value()) * R * temperature.value()),
            ),
        }
    }

    fn caloric(self, temperature: Temperature) -> Result<(f64, f64), IsochoricKineticsError> {
        let cp = self.model.evaluate_heat_capacities(temperature)?;
        let (change, capacity) = match self.constraint {
            GasHeatConstraint::ConstantVolume(_) => (
                self.model
                    .molar_internal_energy_change(self.temperature(), temperature)?
                    .value(),
                cp.cv().value(),
            ),
            GasHeatConstraint::ConstantPressure(_) => (
                self.model
                    .molar_enthalpy_change(self.temperature(), temperature)?
                    .value(),
                cp.cp().value(),
            ),
        };
        let change = self.amount.value() * change;
        let capacity = self.amount.value() * capacity;
        if !change.is_finite() || !positive(capacity) {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite gas caloric increment or capacity",
            ));
        }
        Ok((change, capacity))
    }

    pub(super) fn temperature_at_heat(
        self,
        gas_heat: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Temperature, IsochoricKineticsError> {
        poll(cancelled)?;
        if gas_heat == 0.0 {
            return Ok(self.temperature());
        }
        solve::temperature(
            self.bounds.map(Temperature::value),
            self.temperature().value(),
            controls,
            solve::Side::Either,
            |t| {
                let (change, capacity) = self.caloric(Temperature::new(t))?;
                let mut residual = Sum::default();
                residual.add(change);
                residual.add(gas_heat);
                Ok((residual.value(), capacity))
            },
            cancelled,
        )
        .map(Temperature::new)
    }

    fn entropy_change(
        self,
        initial: Temperature,
        terminal: Temperature,
    ) -> Result<f64, IsochoricKineticsError> {
        let change = match self.constraint {
            GasHeatConstraint::ConstantVolume(_) => self
                .model
                .molar_isochoric_entropy_change(initial, terminal)?,
            GasHeatConstraint::ConstantPressure(p) => {
                self.model.molar_entropy_change(initial, p, terminal, p)?
            }
        };
        let value = self.amount.value() * change.value();
        if !value.is_finite() {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable finite gas entropy change",
            ));
        }
        Ok(value)
    }

    pub(super) fn entropy_transfer(
        self,
        initial: Sum,
        terminal: Sum,
        volume: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, IsochoricKineticsError> {
        let initial = self.temperature_at_heat(initial.value() * volume, controls, cancelled)?;
        let terminal = self.temperature_at_heat(terminal.value() * volume, controls, cancelled)?;
        self.entropy_change(initial, terminal).map(|s| -s)
    }

    pub(super) fn additional_error(
        self,
        base: Sum,
        full: Sum,
        fine: Sum,
        volume: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, IsochoricKineticsError> {
        let base = self
            .temperature_at_heat(base.value() * volume, controls, cancelled)?
            .value();
        let full = self
            .temperature_at_heat(full.value() * volume, controls, cancelled)?
            .value();
        let fine = self
            .temperature_at_heat(fine.value() * volume, controls, cancelled)?
            .value();
        Ok((fine - full).abs() / (controls.relative_tolerance * fine.max(base)))
    }

    pub(super) fn increment_density(
        self,
        step: f64,
        temperature: f64,
        heat: Sum,
        volume: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, IsochoricKineticsError> {
        let mut bath = self.bath;
        bath.temperature = self.temperature_at_heat(heat.value() * volume, controls, cancelled)?;
        bath.increment_density(step, temperature, volume)
    }

    pub(super) fn admit_heat(
        self,
        heat_j: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), IsochoricKineticsError> {
        let temperature = self.temperature_at_heat(heat_j, controls, cancelled)?;
        if (self.caloric(temperature)?.0 + heat_j).abs() > controls.energy_tolerance_j {
            return Err(IsochoricKineticsError::Invalid(
                "finite gas caloric heat residual exceeds tolerance",
            ));
        }
        Ok(())
    }

    fn finish<G>(
        self,
        gas: G,
        heat_j: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FiniteGasHeatBodyStep<'m, G>, IsochoricKineticsError> {
        let temperature = self.temperature_at_heat(heat_j, controls, cancelled)?;
        let n = self.amount.value();
        let enthalpy_change_j = n * self
            .model
            .molar_enthalpy_change(self.temperature(), temperature)?
            .value();
        let internal_energy_change_j = n * self
            .model
            .molar_internal_energy_change(self.temperature(), temperature)?
            .value();
        let boundary_work_j = match self.constraint {
            GasHeatConstraint::ConstantVolume(_) => 0.0,
            GasHeatConstraint::ConstantPressure(_) => {
                (n * R) * (temperature.value() - self.temperature().value())
            }
        };
        let heat_residual_j = match self.constraint {
            GasHeatConstraint::ConstantVolume(_) => internal_energy_change_j + heat_j,
            GasHeatConstraint::ConstantPressure(_) => enthalpy_change_j + heat_j,
        };
        let entropy_change_j_k = self.entropy_change(self.temperature(), temperature)?;
        if [
            enthalpy_change_j,
            internal_energy_change_j,
            boundary_work_j,
            heat_residual_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || heat_residual_j.abs() > controls.energy_tolerance_j
            || (internal_energy_change_j + boundary_work_j + heat_j).abs()
                > controls.energy_tolerance_j
        {
            return Err(IsochoricKineticsError::Invalid(
                "finite gas source heat and pressure work do not close",
            ));
        }
        let mut body = self;
        body.bath.temperature = temperature;
        poll(cancelled)?;
        Ok(FiniteGasHeatBodyStep {
            gas,
            body,
            report: FiniteGasHeatReport {
                heat_input_j: -heat_j,
                enthalpy_change_j,
                internal_energy_change_j,
                boundary_work_j,
                entropy_change_j_k,
                heat_residual_j,
            },
        })
    }
}

/// Original frozen-mixture source increments; positive work is outward.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasHeatReport {
    pub heat_input_j: f64,
    pub enthalpy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub boundary_work_j: f64,
    pub entropy_change_j_k: f64,
    pub heat_residual_j: f64,
}

/// Both original-source gas endpoints publish together after all admissions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasHeatBodyStep<'m, G> {
    gas: G,
    body: FiniteGasHeatBody<'m>,
    report: FiniteGasHeatReport,
}
impl<'m, G> FiniteGasHeatBodyStep<'m, G> {
    pub const fn gas(&self) -> G
    where
        G: Copy,
    {
        self.gas
    }
    pub const fn body(&self) -> FiniteGasHeatBody<'m> {
        self.body
    }
    pub const fn body_report(&self) -> &FiniteGasHeatReport {
        &self.report
    }
    pub fn into_parts(self) -> (G, FiniteGasHeatBody<'m>, FiniteGasHeatReport) {
        (self.gas, self.body, self.report)
    }
}

impl IsochoricReactor {
    /// Evolve fixed-volume chemistry and finite original gas heat storage.
    /// Each gas retains its own mechanical constraint and formation convention.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_gas_body<'a, 'm>(
        &'a mut self,
        amounts: &[Amount],
        volume: Volume,
        initial_temperature: Temperature,
        duration: Time,
        body: FiniteGasHeatBody<'m>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteGasHeatBodyStep<'m, IsochoricKineticsView<'a>>, IsochoricKineticsError> {
        let controls = self.controls();
        let gas = advance_prepared_constraint(
            self,
            amounts,
            None,
            volume,
            None,
            initial_temperature,
            duration,
            Some(HeatExchange::Gas(&body)),
            &mut cancelled,
        )?;
        let heat = gas.balance().heat_input_j;
        body.finish(gas, heat, controls, &mut cancelled)
    }
}

impl IsobaricReactor {
    /// Evolve held-pressure chemistry and a finite frozen-composition gas.
    /// Gas enthalpy includes its pressure work; each partner reports its own.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_gas_body<'a, 'm>(
        &'a mut self,
        amounts: &[Amount],
        pressure: Pressure,
        initial_temperature: Temperature,
        duration: Time,
        body: FiniteGasHeatBody<'m>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteGasHeatBodyStep<'m, IsobaricKineticsView<'a>>, IsochoricKineticsError> {
        let controls = self.controls();
        let gas = self.advance(
            amounts,
            None,
            pressure,
            initial_temperature,
            duration,
            Some(HeatExchange::Gas(&body)),
            &mut cancelled,
        )?;
        let heat = gas.balance().heat_input_j;
        body.finish(gas, heat, controls, &mut cancelled)
    }
}
