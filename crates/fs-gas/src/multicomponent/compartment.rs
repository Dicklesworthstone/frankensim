//! Conservative two-layer control volumes with a moving pressure-equilibrium interface.
//!
//! Each layer retains the original species masses and NASA sensible internal
//! energy. The common pressure and interface follow from the actual composition
//! and variable heat capacities, not a prescribed constant gamma. Advected
//! energy is enthalpy. The interface exchanges equal and opposite p dV work.
//!
//! This is an estimated, well-mixed, low-Mach zone model. It does not resolve a
//! flame, plume, ceiling jet, acoustic wave, or a stratified layer's internal
//! temperature profile. Fire, vent and wall owners supply the extensive rates.
//! The conservation equations are those of NIST TN 1889v1, chapter 2; their
//! constant-Cp reduction is not imposed on the retained caloric sources.

use super::{Error as GasError, MAX_COMPONENTS, MAX_SPECIES, MixtureGas, sum};
use fs_qty::{Area, Energy, Length, MassFlowRate, Power, Pressure, Temperature, Time, Volume};

/// Original-source, reactant-inventory-limited global burning increments.
pub mod burning;

#[cfg(test)]
mod tests;

/// Lower or upper well-mixed gas layer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layer {
    /// Gas below the interface.
    Lower,
    /// Gas above the interface.
    Upper,
}
impl Layer {
    /// Canonical lower/upper index used by all layer arrays.
    pub const fn index(self) -> usize {
        match self {
            Self::Lower => 0,
            Self::Upper => 1,
        }
    }
}

/// One admitted extensive gas inventory. An absent layer has exactly zero mass
/// and volume; it is not populated with artificial species or thermal seeds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerState {
    masses: [f64; MAX_SPECIES],
    species: usize,
    energy: f64,
    volume: f64,
    temperature: f64,
}
impl LayerState {
    /// Canonical species masses [kg].
    pub fn species_masses(&self) -> &[f64] {
        &self.masses[..self.species]
    }
    /// Total gas mass [kg].
    pub fn mass_kg(&self) -> f64 {
        sum(self.species_masses().iter().copied())
    }
    /// Energy relative to the original per-species internal energies at Tmin.
    pub const fn sensible_internal_energy(&self) -> Energy {
        Energy::new(self.energy)
    }
    /// Physical layer volume.
    pub const fn volume(&self) -> Volume {
        Volume::new(self.volume)
    }
    /// Temperature; an absent layer retains its explicit initialization value.
    pub const fn temperature(&self) -> Temperature {
        Temperature::new(self.temperature)
    }
    /// Density [kg/m³], or None for an absent layer.
    pub fn density_kg_m3(&self) -> Option<f64> {
        (self.volume > 0.0).then(|| self.mass_kg() / self.volume)
    }
}

/// Extensive source rates in the same canonical species order and energy gauge
/// as the room. Positive values enter the indicated layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerRates {
    /// Species mass rates [kg/s]. Entries beyond the room's axis must be zero.
    pub species_kg_s: [f64; MAX_SPECIES],
    /// Heat plus advected enthalpy less the associated species references [W].
    pub referenced_energy_w: f64,
}
impl Default for LayerRates {
    fn default() -> Self {
        Self {
            species_kg_s: [0.0; MAX_SPECIES],
            referenced_energy_w: 0.0,
        }
    }
}
impl LayerRates {
    /// Add signed thermal power without changing mass.
    pub fn add_heat(&mut self, heat: Power) -> Result<(), ZoneError> {
        let next = self.referenced_energy_w + heat.value();
        if !next.is_finite() {
            return Err(ZoneError::Invalid("non-finite layer heat rate"));
        }
        self.referenced_energy_w = next;
        Ok(())
    }
    /// Add another already source-matched extensive rate.
    pub fn add(&mut self, other: &Self) -> Result<(), ZoneError> {
        let mut next = *self;
        for (a, b) in next.species_kg_s.iter_mut().zip(other.species_kg_s) {
            *a += b;
        }
        next.referenced_energy_w += other.referenced_energy_w;
        if !next.species_kg_s.iter().all(|v| v.is_finite()) || !next.referenced_energy_w.is_finite()
        {
            return Err(ZoneError::Invalid("non-finite combined layer rates"));
        }
        *self = next;
        Ok(())
    }
    /// Opposite side of the identical conservative transfer.
    pub fn opposite(self) -> Self {
        Self {
            species_kg_s: self.species_kg_s.map(|v| -v),
            referenced_energy_w: -self.referenced_energy_w,
        }
    }
}

/// Bounded pressure/interface solve. Source rates are held over the requested
/// interval; a caller supplies time quadrature and network coupling separately.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZoneSolve {
    /// Maximum safeguarded iterations of the scalar interface-energy equation.
    pub maximum_iterations: usize,
    /// Absolute layer energy-balance tolerance [J].
    pub absolute_energy_tolerance_j: f64,
    /// Relative layer energy-balance tolerance.
    pub relative_energy_tolerance: f64,
}
impl Default for ZoneSolve {
    fn default() -> Self {
        Self {
            maximum_iterations: 80,
            absolute_energy_tolerance_j: 1e-8,
            relative_energy_tolerance: 1e-10,
        }
    }
}

/// One fully published two-zone increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZoneStepReport {
    /// Exact accepted duration.
    pub elapsed: Time,
    /// Iterations of the pressure-equilibrating interface equation.
    pub iterations: usize,
    /// Work done by each layer on the other [J]; these sum to zero exactly.
    pub interface_work_j: [f64; 2],
    /// Layer sensible-energy residuals including interface work [J].
    pub layer_energy_residual_j: [f64; 2],
    /// Whole-room mass residual against all supplied source rates [kg].
    pub mass_residual_kg: f64,
    /// Whole-room physical internal-energy residual against all sources [J].
    pub physical_energy_residual_j: f64,
    /// Final common thermodynamic pressure.
    pub pressure: Pressure,
    /// Final layer interface above the compartment floor.
    pub interface_height: Length,
}

/// A refusal leaves all public layer inventories and the physical clock intact.
#[derive(Clone, Debug, PartialEq)]
pub enum ZoneError {
    /// Bad geometry, composition, source identity, numerical controls or rate.
    Invalid(&'static str),
    /// Original source owner refused the state or caloric inversion.
    Gas(GasError),
    /// Negative species inventory would result from the proposed source duration.
    SpeciesDepletion { layer: Layer, species: usize },
    /// Source support cannot contain the requested thermodynamic endpoint.
    CaloricDomain,
    /// No interface inside the room satisfies both original layer energy balances.
    InterfaceDomain,
    /// The bounded scalar solve did not meet the declared balance tolerance.
    IterationLimit,
    /// Cancellation occurred before publication.
    Cancelled,
}
impl core::fmt::Display for ZoneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "two-zone gas: {self:?}")
    }
}
impl std::error::Error for ZoneError {}
impl From<GasError> for ZoneError {
    fn from(value: GasError) -> Self {
        match value {
            GasError::Cancelled => Self::Cancelled,
            other => Self::Gas(other),
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ZoneError> {
    if cancelled() {
        Err(ZoneError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}

/// A fixed-volume compartment with two original-source gas inventories and a
/// mechanically equilibrated moving interface. No mixture heat capacity, gas
/// constant, formation energy or gamma is authored a second time.
#[derive(Clone, Debug, PartialEq)]
pub struct TwoZone {
    gas: MixtureGas,
    floor_area: f64,
    height: f64,
    pressure: f64,
    layers: [LayerState; 2],
    time: f64,
}
impl TwoZone {
    /// Initialize both layers from explicit temperatures/compositions at one
    /// pressure. Interface endpoints admit an absent upper or lower layer.
    pub fn new(
        gas: MixtureGas,
        floor_area: Area,
        height: Length,
        interface: Length,
        pressure: Pressure,
        temperatures: [Temperature; 2],
        mass_fractions: [Vec<f64>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ZoneError> {
        poll(&mut cancelled)?;
        let (a, h, z, p) = (
            floor_area.value(),
            height.value(),
            interface.value(),
            pressure.value(),
        );
        if !positive(a)
            || !positive(h)
            || !positive(a * h)
            || !positive(p)
            || !(p * (a * h)).is_finite()
            || !z.is_finite()
            || z < 0.0
            || z > h
        {
            return Err(ZoneError::Invalid(
                "room area, height, interface or pressure",
            ));
        }
        let n = gas.species_count();
        let volumes = [a * z, a * (h - z)];
        let mut layers = [LayerState {
            masses: [0.0; MAX_SPECIES],
            species: n,
            energy: 0.0,
            volume: 0.0,
            temperature: 0.0,
        }; 2];
        for k in 0..2 {
            poll(&mut cancelled)?;
            let y = &mass_fractions[k];
            let t = temperatures[k].value();
            let d = gas.domain();
            if y.len() != n
                || y.iter().any(|v| !v.is_finite() || *v < 0.0)
                || (sum(y.iter().copied()) - 1.0).abs() > 1e-12
                || !(d.minimum_temperature..=d.maximum_temperature).contains(&t)
            {
                return Err(ZoneError::Invalid(
                    "initial layer composition or temperature",
                ));
            }
            let inv_mass = sum(y.iter().zip(gas.inverse_molar_masses()).map(|(y, m)| y * m));
            let rho = p / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t * inv_mass);
            let mut q = [0.0; MAX_COMPONENTS];
            let mut densities = [0.0; MAX_SPECIES];
            for i in 0..n {
                densities[i] = rho * y[i];
                layers[k].masses[i] = densities[i] * volumes[k];
            }
            gas.state_at(&densities[..n], 0.0, p, &mut q[..n + 2])?;
            layers[k].energy = q[n + 1] * volumes[k];
            layers[k].volume = volumes[k];
            layers[k].temperature = t;
            if !layers[k].masses[..n].iter().all(|m| m.is_finite())
                || !layers[k].mass_kg().is_finite()
                || !layers[k].energy.is_finite()
            {
                return Err(ZoneError::Invalid(
                    "unrepresentable initial extensive inventory",
                ));
            }
        }
        if !(layers[0].energy + layers[1].energy).is_finite()
            || !(layers[0].mass_kg() + layers[1].mass_kg()).is_finite()
        {
            return Err(ZoneError::Invalid("unrepresentable room inventory"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            gas,
            floor_area: a,
            height: h,
            pressure: p,
            layers,
            time: 0.0,
        })
    }
    /// Original complete caloric/species owner.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Read one admitted layer.
    pub const fn layer(&self, layer: Layer) -> &LayerState {
        &self.layers[layer.index()]
    }
    /// Fixed total room volume.
    pub fn volume(&self) -> Volume {
        Volume::new(self.floor_area * self.height)
    }
    /// Floor area used for the moving layer height.
    pub const fn floor_area(&self) -> Area {
        Area::new(self.floor_area)
    }
    /// Room height.
    pub const fn height(&self) -> Length {
        Length::new(self.height)
    }
    /// Common thermodynamic pressure. Hydrostatic vent profiles are separate.
    pub const fn pressure(&self) -> Pressure {
        Pressure::new(self.pressure)
    }
    /// Lower-layer depth above the floor.
    pub fn interface_height(&self) -> Length {
        Length::new(self.layers[0].volume / self.floor_area)
    }
    /// Physical accepted clock.
    pub const fn time(&self) -> Time {
        Time::new(self.time)
    }
    /// Reuse a prepared compartment's storage for another state of the identical
    /// source and geometry. This supports atomic multi-room nonlinear trials.
    pub fn copy_state_from(&mut self, source: &Self) -> Result<(), ZoneError> {
        if self.gas != source.gas
            || self.floor_area != source.floor_area
            || self.height != source.height
        {
            return Err(ZoneError::Invalid(
                "compartment workspace source or geometry mismatch",
            ));
        }
        self.layers = source.layers;
        self.pressure = source.pressure;
        self.time = source.time;
        Ok(())
    }
    /// Whole-room sensible internal energy in the retained species gauge.
    pub fn sensible_internal_energy(&self) -> Energy {
        Energy::new(self.layers[0].energy + self.layers[1].energy)
    }
    /// Whole-room physical internal energy, retaining species formation energies.
    pub fn physical_internal_energy(&self) -> Energy {
        Energy::new(
            self.sensible_internal_energy().value()
                + sum((0..self.gas.species_count()).map(|i| {
                    (self.layers[0].masses[i] + self.layers[1].masses[i])
                        * self.gas.energy_references()[i]
                })),
        )
    }
    /// Build an incoming rate from this exact source bank, mass composition and
    /// donor temperature. Its energy is partial enthalpy, not internal energy.
    pub fn inflow(
        &self,
        mass_rate: MassFlowRate,
        mass_fractions: &[f64],
        temperature: Temperature,
    ) -> Result<LayerRates, ZoneError> {
        let n = self.gas.species_count();
        let mdot = mass_rate.value();
        if !mdot.is_finite()
            || mdot < 0.0
            || mass_fractions.len() != n
            || mass_fractions.iter().any(|v| !v.is_finite() || *v < 0.0)
            || (sum(mass_fractions.iter().copied()) - 1.0).abs() > 1e-12
        {
            return Err(ZoneError::Invalid("inflow rate or mass composition"));
        }
        let domain = self.gas.domain();
        if !(domain.minimum_temperature..=domain.maximum_temperature).contains(&temperature.value())
        {
            return Err(ZoneError::CaloricDomain);
        }
        let mut out = LayerRates::default();
        let t0 = Temperature::new(domain.minimum_temperature);
        let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        for (i, model) in self.gas.source().models().iter().enumerate() {
            out.species_kg_s[i] = mdot * mass_fractions[i];
            let h = (model
                .molar_internal_energy_change(t0, temperature)
                .map_err(GasError::from)?
                .value()
                + r * temperature.value())
                * self.gas.inverse_molar_masses()[i];
            out.referenced_energy_w += out.species_kg_s[i] * h;
        }
        if !out.referenced_energy_w.is_finite() {
            return Err(ZoneError::Invalid("inflow enthalpy overflow"));
        }
        Ok(out)
    }
    /// Positive transfer carried out of one layer, ready to be added to the
    /// receiver and negated for the donor using the identical values.
    pub fn carried_rate(
        &self,
        layer: Layer,
        mass_rate: MassFlowRate,
    ) -> Result<LayerRates, ZoneError> {
        let state = self.layer(layer);
        let mass = state.mass_kg();
        if !positive(mass) {
            return Err(ZoneError::Invalid("cannot advect from an absent layer"));
        }
        let mut fractions = [0.0; MAX_SPECIES];
        for (i, y) in fractions[..self.gas.species_count()].iter_mut().enumerate() {
            *y = state.masses[i] / mass;
        }
        self.inflow(
            mass_rate,
            &fractions[..self.gas.species_count()],
            state.temperature(),
        )
    }
    /// Advance one exact supplied-source duration atomically. An external driver
    /// must recompute its rates when retrying a shorter interval.
    pub fn advance(
        &mut self,
        duration: Time,
        rates: [LayerRates; 2],
        controls: ZoneSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ZoneStepReport, ZoneError> {
        poll(&mut cancelled)?;
        let dt = duration.value();
        let n = self.gas.species_count();
        let volume = self.volume().value();
        if !positive(dt)
            || !positive(self.time + dt - self.time)
            || !(self.time + dt).is_finite()
            || !(1..=256).contains(&controls.maximum_iterations)
            || !positive(controls.absolute_energy_tolerance_j)
            || !positive(controls.relative_energy_tolerance)
            || controls.relative_energy_tolerance > 1e-3
        {
            return Err(ZoneError::Invalid("zone duration or solve controls"));
        }
        let mut next = self.layers;
        let mut targets = [0.0; 2];
        let mut maximum = [0.0; 2];
        for k in 0..2 {
            if !rates[k].referenced_energy_w.is_finite()
                || rates[k].species_kg_s.iter().any(|v| !v.is_finite())
                || rates[k].species_kg_s[n..].iter().any(|v| *v != 0.0)
            {
                return Err(ZoneError::Invalid("layer source axis or rate"));
            }
            for i in 0..n {
                let m = self.layers[k].masses[i] + dt * rates[k].species_kg_s[i];
                if !m.is_finite() || m < 0.0 {
                    return Err(ZoneError::SpeciesDepletion {
                        layer: if k == 0 { Layer::Lower } else { Layer::Upper },
                        species: i,
                    });
                }
                next[k].masses[i] = m;
            }
            targets[k] = self.layers[k].energy + dt * rates[k].referenced_energy_w;
            maximum[k] =
                sum((0..n).map(|i| next[k].masses[i] * self.gas.maximum_sensible_energies()[i]));
        }
        let total = targets[0] + targets[1];
        if !total.is_finite()
            || total < 0.0
            || !maximum.iter().all(|v| v.is_finite())
            || total > maximum[0] + maximum[1]
        {
            return Err(ZoneError::CaloricDomain);
        }
        let scale = total.abs() + self.pressure * volume;
        let tolerance =
            controls.absolute_energy_tolerance_j + controls.relative_energy_tolerance * scale;
        let mut lower = 0.0_f64.max(total - maximum[1]);
        let mut upper = total.min(maximum[0]);
        if lower > upper {
            return Err(ZoneError::CaloricDomain);
        }
        let left = evaluate_partition(self, &next, total, lower, targets[0], &mut cancelled)?;
        let right = evaluate_partition(self, &next, total, upper, targets[0], &mut cancelled)?;
        if left.residual > tolerance || right.residual < -tolerance {
            return Err(ZoneError::InterfaceDomain);
        }
        let mut energy = targets[0].max(lower).min(upper);
        let mut accepted = None;
        for iteration in 0..controls.maximum_iterations {
            poll(&mut cancelled)?;
            let evaluation =
                evaluate_partition(self, &next, total, energy, targets[0], &mut cancelled)?;
            if evaluation.residual.abs() <= tolerance {
                accepted = Some((evaluation, energy, iteration + 1));
                break;
            }
            if evaluation.residual < 0.0 {
                lower = energy;
            } else {
                upper = energy;
            }
            let trial = energy - evaluation.residual / evaluation.derivative;
            energy = if trial.is_finite() && trial > lower && trial < upper {
                trial
            } else {
                lower + 0.5 * (upper - lower)
            };
        }
        let (evaluation, energy, iterations) = accepted.ok_or(ZoneError::IterationLimit)?;
        let energies = [energy, total - energy];
        for k in 0..2 {
            next[k].energy = energies[k];
            next[k].volume = evaluation.volumes[k];
            next[k].temperature = evaluation.temperatures[k];
        }
        let work =
            0.5 * (self.pressure + evaluation.pressure) * (next[0].volume - self.layers[0].volume);
        let residuals = [
            next[0].energy - targets[0] + work,
            next[1].energy - targets[1] - work,
        ];
        let mass_residual = sum((0..n).map(|i| {
            sum((0..2).map(|k| {
                next[k].masses[i] - self.layers[k].masses[i] - dt * rates[k].species_kg_s[i]
            }))
        }));
        let physical_residual =
            (next[0].energy + next[1].energy - self.layers[0].energy - self.layers[1].energy)
                - dt * (rates[0].referenced_energy_w + rates[1].referenced_energy_w)
                + sum((0..n).map(|i| {
                    (next[0].masses[i] + next[1].masses[i]
                        - self.layers[0].masses[i]
                        - self.layers[1].masses[i]
                        - dt * (rates[0].species_kg_s[i] + rates[1].species_kg_s[i]))
                        * self.gas.energy_references()[i]
                }));
        if !work.is_finite()
            || !mass_residual.is_finite()
            || !physical_residual.is_finite()
            || residuals
                .iter()
                .any(|r| !r.is_finite() || r.abs() > 2.0 * tolerance)
        {
            return Err(ZoneError::Invalid("unrepresentable zone balances"));
        }
        poll(&mut cancelled)?;
        self.layers = next;
        self.pressure = evaluation.pressure;
        self.time += dt;
        Ok(ZoneStepReport {
            elapsed: duration,
            iterations,
            interface_work_j: [work, -work],
            layer_energy_residual_j: residuals,
            mass_residual_kg: mass_residual,
            physical_energy_residual_j: physical_residual,
            pressure: self.pressure(),
            interface_height: self.interface_height(),
        })
    }
}

struct Partition {
    residual: f64,
    derivative: f64,
    pressure: f64,
    volumes: [f64; 2],
    temperatures: [f64; 2],
}
fn evaluate_partition(
    room: &TwoZone,
    next: &[LayerState; 2],
    total: f64,
    lower_energy: f64,
    lower_target: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Partition, ZoneError> {
    poll(cancelled)?;
    let n = room.gas.species_count();
    let volume = room.volume().value();
    let energies = [lower_energy, total - lower_energy];
    let mut primitives = [None; 2];
    let mut pressures = [0.0; 2];
    let mut temperatures = [0.0; 2];
    let mut gamma_minus_one = [0.0; 2];
    for k in 0..2 {
        temperatures[k] = room.layers[k].temperature;
        if next[k].mass_kg() == 0.0 {
            if energies[k] != 0.0 {
                return Err(ZoneError::CaloricDomain);
            }
            continue;
        }
        // Evaluate the existing EOS on a common reference volume. Its partial
        // pressures then sum to the actual common room pressure; no new caloric
        // inversion, frozen gamma or per-species formation subtraction is used.
        let mut q = [0.0; MAX_COMPONENTS];
        for i in 0..n {
            q[i] = next[k].masses[i] / volume;
        }
        q[n + 1] = energies[k] / volume;
        let p = room
            .gas
            .decode(&q[..n + 2], Some(room.layers[k].temperature), cancelled)?;
        pressures[k] = p.pressure;
        temperatures[k] = p.temperature;
        gamma_minus_one[k] = p.sound_speed * p.sound_speed * p.density / p.pressure - 1.0;
        primitives[k] = Some(p);
    }
    let pressure = pressures[0] + pressures[1];
    if !positive(pressure) {
        return Err(ZoneError::Invalid("empty or unrepresentable room pressure"));
    }
    let lower_volume = volume * (pressures[0] / pressure);
    let volumes = [lower_volume, volume - lower_volume];
    if primitives[1].is_none() && volumes[1] != 0.0 {
        return Err(ZoneError::InterfaceDomain);
    }
    let work = 0.5 * (room.pressure + pressure) * (lower_volume - room.layers[0].volume);
    // On positive-Cv source paths this derivative is positive: it is one plus
    // a nonnegative weighted combination of the two actual gamma-1 values.
    let ratio = room.pressure / pressure;
    let weight = ratio * (lower_volume / volume) + room.layers[0].volume / volume;
    let derivative = 1.0 + 0.5 * (1.0 + ratio) * gamma_minus_one[0]
        - 0.5 * weight * (gamma_minus_one[0] - gamma_minus_one[1]);
    Ok(Partition {
        residual: lower_energy - lower_target + work,
        derivative,
        pressure,
        volumes,
        temperatures,
    })
}
