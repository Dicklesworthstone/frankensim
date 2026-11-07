use super::*;
use fs_math::det::expm1;
use fs_qty::Amount;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateErrorV1;
use fs_thermochem::kinetics::atomic_ionization::AtomicIonizationErrorV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    PhotoionizationErrorV1, PhotoionizationEvaluationV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::RadiativeRecombinationErrorV1;
use fs_thermochem::kinetics::{KineticsErrorV1, admit_caloric_interval};

const CV: f64 = 1.5 * R;

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let increment = value - self.correction;
        let next = self.value + increment;
        self.correction = (next - self.value) - increment;
        self.value = next;
    }
}

// Net, photo, net collisional, capture, absorbed energy, captured kinetic energy.
#[derive(Clone, Copy, Default)]
pub(super) struct Events([Sum; 6]);

pub(super) struct Point {
    pub amounts: Vec<f64>,
    pub events: Vec<Events>,
    te: f64,
}
impl Point {
    pub fn new(species: usize, channels: usize) -> Self {
        Self {
            amounts: vec![0.0; species],
            events: vec![Events::default(); channels],
            te: 0.0,
        }
    }
    pub fn copy_from(&mut self, source: &Self) {
        self.amounts.copy_from_slice(&source.amounts);
        self.events.copy_from_slice(&source.events);
        self.te = source.te;
    }
}

struct Edge<'a> {
    source: &'a AtomicRadiativeRecombinationV1,
    lower: usize,
    upper: usize,
    photon: PhotoionizationEvaluationV1,
}

pub(super) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    edges: Vec<Edge<'a>>,
    electron: Option<usize>,
    ladders: Vec<Vec<usize>>,
    level_energies: [f64; MAX_EQUILIBRIUM_SPECIES_V1],
    controls: IonizationNetworkSolve,
}
impl<'a> Problem<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        sources: &'a [AtomicRadiativeRecombinationV1],
        photons: &[PhotonNumberFluxV1],
        controls: IonizationNetworkSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, IonizationNetworkError> {
        let electron = initial.electron_index();
        if initial.relative_charge_residual() > controls.local.relative_charge_tolerance {
            return Err(IonizationNetworkError::Invalid(
                "initial network charge residual",
            ));
        }
        let mut result = Self {
            initial,
            edges: Vec::with_capacity(sources.len()),
            electron,
            ladders: Vec::new(),
            level_energies: [0.0; MAX_EQUILIBRIUM_SPECIES_V1],
            controls,
        };
        if !sources.is_empty() {
            result.admit_temperature(initial.electron_temperature().value())?;
        }
        let mut incoming = [None; MAX_EQUILIBRIUM_SPECIES_V1];
        let mut outgoing = [None; MAX_EQUILIBRIUM_SPECIES_V1];
        for (ordinal, source) in sources.iter().enumerate() {
            poll(cancelled)?;
            let atom = source.source().channel();
            if atom.thermodynamics() != initial.thermodynamics()
                || electron != Some(atom.electron_index())
            {
                return Err(IonizationNetworkError::Invalid(
                    "network requires one original bank and electron",
                ));
            }
            let lower = atom.lower_index();
            let upper = atom.upper_index();
            let charge = initial.thermodynamics().charge_vector();
            if atom.lower_charge() < 0
                || charge.get(upper) != atom.lower_charge().checked_add(1)
                || incoming[upper].is_some()
                || outgoing[lower].is_some()
                || sources[..ordinal].iter().any(|other| {
                    other
                        .source()
                        .channel()
                        .common_temperature_mechanism()
                        .reactions()[0]
                        .reaction()
                        == atom.common_temperature_mechanism().reactions()[0].reaction()
                })
            {
                return Err(IonizationNetworkError::Invalid(
                    "duplicate/branching edge, reaction identity or adjacent charge",
                ));
            }
            incoming[upper] = Some(ordinal);
            outgoing[lower] = Some(ordinal);
            let photon = source
                .source()
                .evaluate(photons, &mut *cancelled)
                .map_err(photo_error)?;
            result.edges.push(Edge {
                source,
                lower,
                upper,
                photon,
            });
        }
        result.edges.sort_by_key(|edge| (edge.lower, edge.upper));
        incoming.fill(None);
        outgoing.fill(None);
        for (index, edge) in result.edges.iter().enumerate() {
            incoming[edge.upper] = Some(index);
            outgoing[edge.lower] = Some(index);
        }
        // A ladder's original chemical level relative to its root is
        // F_q-F_root+(q-q_root)*F_e = sum(edge gaps). Conserved nuclei and
        // charge remove elemental/charge reference offsets before subtraction.
        // Endpoint energy will use these levels with actual population changes,
        // independently of the accumulated reaction extents.
        let mut visited = 0;
        for root in 0..initial.amounts().len() {
            if incoming[root].is_some() || outgoing[root].is_none() {
                continue;
            }
            let mut ladder = vec![root];
            let mut node = root;
            let mut energy = Sum::default();
            while let Some(index) = outgoing[node] {
                poll(cancelled)?;
                visited += 1;
                if visited > result.edges.len() {
                    return Err(IonizationNetworkError::Invalid("cyclic atomic network"));
                }
                let edge = &result.edges[index];
                node = edge.upper;
                energy.add(edge.source.source().channel().ionization_energy_j_mol());
                if !positive(energy.value) {
                    return Err(IonizationNetworkError::Arithmetic("network level gap"));
                }
                result.level_energies[node] = energy.value;
                ladder.push(node);
            }
            result.ladders.push(ladder);
        }
        if visited != result.edges.len() {
            return Err(IonizationNetworkError::Invalid("cyclic atomic network"));
        }
        for (index, model) in initial.thermodynamics().models().iter().enumerate() {
            let (lo, hi) = if Some(index) == electron && !sources.is_empty() {
                (
                    Temperature::new(controls.local.minimum_electron_temperature_k),
                    Temperature::new(controls.local.maximum_electron_temperature_k),
                )
            } else if incoming[index].is_some() || outgoing[index].is_some() {
                (initial.heavy_temperature(), initial.heavy_temperature())
            } else {
                continue;
            };
            admit_caloric_interval(model, lo, hi, &mut *cancelled).map_err(
                |error| match error {
                    KineticsErrorV1::Cancelled => IonizationNetworkError::Cancelled,
                    error => IonizationNetworkError::Admission(error),
                },
            )?;
        }
        Ok(result)
    }

    pub fn sources(&self) -> impl Iterator<Item = &AtomicRadiativeRecombinationV1> {
        self.edges.iter().map(|edge| edge.source)
    }

    pub fn initial_point(&self) -> Point {
        Point {
            amounts: self.initial.amounts().iter().map(|n| n.value()).collect(),
            events: vec![Events::default(); self.edges.len()],
            te: self.initial.electron_temperature().value(),
        }
    }

    pub fn inactive(&self, point: &Point) -> bool {
        self.edges.is_empty()
            || (self.electron.is_none_or(|i| point.amounts[i] == 0.0)
                && self.edges.iter().all(|edge| {
                    point.amounts[edge.lower] == 0.0 || edge.photon.rate_per_neutral_per_s == 0.0
                }))
    }

    pub fn compose(
        &self,
        point: &mut Point,
        dt: f64,
        balance: &mut IonizationNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), IonizationNetworkError> {
        let last = self.edges.len() - 1;
        for index in 0..last {
            self.midpoint(index, point, 0.5 * dt, balance, cancelled)?;
        }
        self.midpoint(last, point, dt, balance, cancelled)?;
        for index in (0..last).rev() {
            self.midpoint(index, point, 0.5 * dt, balance, cancelled)?;
        }
        Ok(())
    }

    fn midpoint(
        &self,
        index: usize,
        point: &mut Point,
        dt: f64,
        balance: &mut IonizationNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), IonizationNetworkError> {
        poll(cancelled)?;
        let edge = &self.edges[index];
        let electron = self.electron.expect("admitted structural electron");
        let pair = Pair {
            lower: point.amounts[edge.lower],
            upper: point.amounts[edge.upper],
            electron: point.amounts[electron],
            te: point.te,
        };
        if pair.lower == 0.0 && pair.upper == 0.0 {
            return Ok(());
        }
        let rates = self.rates(edge, pair, balance, cancelled)?;
        let middle = self.frozen(edge, pair, rates, 0.5 * dt)?.0;
        let middle_rates = self.rates(edge, middle, balance, cancelled)?;
        let (next, increments) = self.frozen(edge, pair, middle_rates, dt)?;
        point.amounts[edge.lower] = next.lower;
        point.amounts[edge.upper] = next.upper;
        point.amounts[electron] = next.electron;
        point.te = next.te;
        for (sum, increment) in point.events[index].0.iter_mut().zip(increments) {
            sum.add(increment);
            if !sum.value.is_finite() {
                return Err(IonizationNetworkError::Arithmetic("network event integral"));
            }
        }
        Ok(())
    }

    fn rates(
        &self,
        edge: &Edge<'_>,
        pair: Pair,
        balance: &mut IonizationNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Rates, IonizationNetworkError> {
        poll(cancelled)?;
        if balance.rate_evaluations >= self.controls.maximum_rate_evaluations {
            return Err(IonizationNetworkError::StepBudget);
        }
        balance.rate_evaluations += 1;
        let photon = edge.photon;
        let te = if pair.electron == 0.0 && photon.rate_per_neutral_per_s > 0.0 {
            (NA / CV) * (photon.electron_heating_per_neutral_w / photon.rate_per_neutral_per_s)
        } else {
            pair.te
        };
        self.admit_temperature(te)?;
        let volume = self.initial.volume().value();
        let atom = edge
            .source
            .source()
            .channel()
            .evaluate(
                self.initial.heavy_temperature(),
                Temperature::new(te),
                [
                    pair.lower / volume,
                    pair.upper / volume,
                    pair.electron / volume,
                ],
                &mut *cancelled,
            )
            .map_err(|error| match error {
                AtomicIonizationErrorV1::Cancelled
                | AtomicIonizationErrorV1::Kinetics(KineticsErrorV1::Cancelled) => {
                    IonizationNetworkError::Cancelled
                }
                error => IonizationNetworkError::Atomic(error),
            })?;
        let radiative = edge
            .source
            .evaluate(
                Temperature::new(te),
                self.controls.local.quadrature,
                &mut *cancelled,
            )
            .map_err(|error| match error {
                RadiativeRecombinationErrorV1::Cancelled => IonizationNetworkError::Cancelled,
                error => IonizationNetworkError::Recombination(error),
            })?;
        balance.maximum_relative_quadrature_error = balance
            .maximum_relative_quadrature_error
            .max(radiative.relative_quadrature_error_estimate);
        balance.maximum_relative_tail_bound = balance
            .maximum_relative_tail_bound
            .max(radiative.relative_tail_bound);
        let minor = pair.upper.min(pair.electron);
        let three_body = if pair.upper <= pair.electron {
            -atom.concentration_derivatives_per_s[1]
        } else if minor > 0.0 {
            atom.reverse_rate_mol_m3_s / (minor / volume)
        } else {
            0.0
        };
        let capture = NA * radiative.coefficient_m3_s * (pair.upper.max(pair.electron) / volume);
        let rates = Rates {
            birth: photon.rate_per_neutral_per_s + atom.concentration_derivatives_per_s[0],
            death: three_body + capture,
            capture,
            captured_energy: NA * radiative.mean_captured_energy_j().unwrap_or(0.0),
        };
        if [
            rates.birth,
            rates.death,
            rates.capture,
            rates.captured_energy,
        ]
        .iter()
        .any(|v| !nonnegative(*v))
        {
            return Err(IonizationNetworkError::Arithmetic(
                "network production/capture rates",
            ));
        }
        Ok(rates)
    }

    fn frozen(
        &self,
        edge: &Edge<'_>,
        pair: Pair,
        rates: Rates,
        dt: f64,
    ) -> Result<(Pair, [f64; 6]), IonizationNetworkError> {
        let lambda = rates.birth + rates.death;
        if !nonnegative(lambda) {
            return Err(IonizationNetworkError::Arithmetic("network turnover rate"));
        }
        if lambda == 0.0 || dt == 0.0 {
            return Ok((pair, [0.0; 6]));
        }
        let minor = pair.upper.min(pair.electron);
        let total = pair.lower + minor;
        if !nonnegative(total) {
            return Err(IonizationNetworkError::Arithmetic("network pair inventory"));
        }
        let dose = lambda * dt;
        let depleted = if dose.is_finite() { -expm1(-dose) } else { 1.0 };
        let surviving = if dose.is_finite() { exp(-dose) } else { 0.0 };
        let relaxation = depleted / lambda;
        let remainder = if dose < 1e-3 {
            dt * dose
                * (0.5
                    + dose
                        * (-1.0 / 6.0 + dose * (1.0 / 24.0 + dose * (-1.0 / 120.0 + dose / 720.0))))
        } else {
            dt - relaxation
        };
        let birth = rates.birth / lambda;
        let death = rates.death / lambda;
        let lower_time = pair.lower * relaxation + total * death * remainder;
        let minor_time = minor * relaxation + total * birth * remainder;
        let mut lower = pair.lower * surviving + total * death * depleted;
        let mut next_minor = minor * surviving + total * birth * depleted;
        if next_minor <= lower {
            lower = total - next_minor;
        } else {
            next_minor = total - lower;
        }
        let extent = (birth * pair.lower - death * minor) * depleted;
        let photo = edge.photon.rate_per_neutral_per_s * lower_time;
        let capture = rates.capture * minor_time;
        let collisional = (extent - photo) + capture;
        let absorbed = NA * edge.photon.absorbed_power_per_neutral_w * lower_time;
        let captured = rates.captured_energy * capture;
        let sensible = NA * edge.photon.electron_heating_per_neutral_w * lower_time
            - edge.source.source().channel().ionization_energy_j_mol() * collisional
            - captured;
        let electron = next_minor + (pair.electron - pair.upper).max(0.0);
        let upper = next_minor + (pair.upper - pair.electron).max(0.0);
        if [
            lower, upper, electron, lower_time, minor_time, photo, capture, absorbed, captured,
        ]
        .iter()
        .any(|v| !nonnegative(*v))
            || !sensible.is_finite()
            || !extent.is_finite()
            || !collisional.is_finite()
            || (electron == 0.0 && (pair.electron != 0.0 || sensible != 0.0))
        {
            return Err(IonizationNetworkError::OutsideDomain);
        }
        let te = if electron == 0.0 {
            pair.te
        } else if pair.electron == 0.0 {
            sensible / (CV * electron)
        } else {
            let increment = (sensible - CV * pair.te * extent) / (CV * electron);
            if increment.is_finite() && increment.abs() < 0.5 * pair.te {
                pair.te + increment
            } else {
                (CV * pair.electron * pair.te + sensible) / (CV * electron)
            }
        };
        self.admit_temperature(te)?;
        Ok((
            Pair {
                lower,
                upper,
                electron,
                te,
            },
            [extent, photo, collisional, capture, absorbed, captured],
        ))
    }

    fn admit_temperature(&self, temperature: f64) -> Result<(), IonizationNetworkError> {
        if !positive(temperature)
            || !(self.controls.local.minimum_electron_temperature_k
                ..=self.controls.local.maximum_electron_temperature_k)
                .contains(&temperature)
        {
            Err(IonizationNetworkError::OutsideDomain)
        } else {
            Ok(())
        }
    }

    pub fn error(
        &self,
        start: &Point,
        coarse: &Point,
        fine: &Point,
    ) -> Result<f64, IonizationNetworkError> {
        let controls = self.controls.local;
        let mut error = (coarse.te - fine.te).abs()
            / (controls.temperature_tolerance_k
                + controls.relative_tolerance * coarse.te.max(fine.te));
        for ((initial, full), next) in start.amounts.iter().zip(&coarse.amounts).zip(&fine.amounts)
        {
            error = error.max(
                (full - next).abs()
                    / (controls.absolute_amount_tolerance_mol
                        + controls.relative_tolerance * initial.max(*full).max(*next)),
            );
        }
        for ((old, full), next) in start.events.iter().zip(&coarse.events).zip(&fine.events) {
            for slot in 0..6 {
                let absolute = if slot < 4 {
                    controls.absolute_amount_tolerance_mol
                } else {
                    controls.energy_tolerance_j
                };
                error = error.max(
                    (full.0[slot].value - next.0[slot].value).abs()
                        / (absolute
                            + controls.relative_tolerance
                                * old.0[slot]
                                    .value
                                    .abs()
                                    .max(full.0[slot].value.abs())
                                    .max(next.0[slot].value.abs())),
                );
            }
        }
        if !nonnegative(error) {
            return Err(IonizationNetworkError::Arithmetic(
                "network local error estimate",
            ));
        }
        Ok(error / 3.0)
    }

    pub fn publish(
        &self,
        point: &Point,
        balance: &mut IonizationNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(FrozenIonizedGasStateV1, Vec<IonizationChannelBalance>), IonizationNetworkError>
    {
        let mut expected = [Sum::default(); MAX_EQUILIBRIUM_SPECIES_V1];
        let mut totals = [Sum::default(); 4];
        let mut channels = Vec::with_capacity(self.edges.len());
        for (edge, events) in self.edges.iter().zip(&point.events) {
            poll(cancelled)?;
            let [net, photo, collisional, capture, absorbed, captured] =
                events.0.map(|sum| sum.value);
            expected[edge.lower].add(-net);
            expected[edge.upper].add(net);
            expected[self.electron.expect("admitted electron")].add(net);
            let emitted =
                edge.source.source().channel().ionization_energy_j_mol() * capture + captured;
            for (sum, value) in totals.iter_mut().zip([photo, capture, absorbed, emitted]) {
                sum.add(value);
            }
            channels.push(IonizationChannelBalance {
                net_extent_mol: net,
                collisional_extent_mol: collisional,
                photoionization_extent_mol: photo,
                radiative_recombination_extent_mol: capture,
                absorbed_radiation_energy_j: absorbed,
                emitted_radiation_energy_j: emitted,
            });
        }
        let mut chemical = Sum::default();
        let mut heavy_count = Sum::default();
        for (index, (before, after)) in self
            .initial
            .amounts()
            .iter()
            .zip(&point.amounts)
            .enumerate()
        {
            poll(cancelled)?;
            let change = after - before.value();
            let residual = (change - expected[index].value).abs();
            let allowance = self.controls.local.absolute_amount_tolerance_mol
                + 128.0 * f64::EPSILON * before.value().max(*after);
            if !nonnegative(*after) || residual > allowance {
                return Err(IonizationNetworkError::Arithmetic(
                    "network endpoint stoichiometry",
                ));
            }
            balance.maximum_stoichiometric_residual_mol =
                balance.maximum_stoichiometric_residual_mol.max(residual);
            chemical.add(self.level_energies[index] * change);
            if Some(index) != self.electron {
                heavy_count.add(change);
            }
        }
        for ladder in &self.ladders {
            let mut change = Sum::default();
            let mut initial = Sum::default();
            for &index in ladder {
                change.add(point.amounts[index] - self.initial.amounts()[index].value());
                initial.add(self.initial.amounts()[index].value());
            }
            balance.maximum_nucleus_residual_mol =
                balance.maximum_nucleus_residual_mol.max(change.value.abs());
            if change.value.abs()
                > self.controls.local.absolute_amount_tolerance_mol
                    + 128.0 * f64::EPSILON * initial.value
            {
                return Err(IonizationNetworkError::Arithmetic(
                    "network nucleus inventory",
                ));
            }
        }
        let sensible = self.electron.map_or(0.0, |index| {
            let before = self.initial.amounts()[index].value();
            CV * (before * (point.te - self.initial.electron_temperature().value())
                + (point.amounts[index] - before) * point.te)
        });
        let heavy = CV * self.initial.heavy_temperature().value() * heavy_count.value;
        let [photo, capture, absorbed, emitted] = totals.map(|sum| sum.value);
        let mut residual = Sum::default();
        for value in [sensible, chemical.value, heavy, emitted, -absorbed] {
            residual.add(value);
        }
        if [photo, capture, absorbed, emitted]
            .iter()
            .any(|v| !nonnegative(*v))
            || !nonnegative(NA * photo)
            || !nonnegative(NA * capture)
            || !residual.value.is_finite()
        {
            return Err(IonizationNetworkError::Arithmetic(
                "network photon/gas energy",
            ));
        }
        if residual.value.abs() > self.controls.local.energy_tolerance_j {
            return Err(IonizationNetworkError::EnergyBalance {
                residual_j: residual.value,
            });
        }
        let amounts: Vec<_> = point.amounts.iter().map(|&n| Amount::new(n)).collect();
        let state = if amounts == self.initial.amounts()
            && point.te == self.initial.electron_temperature().value()
        {
            self.initial.clone()
        } else {
            FrozenIonizedGasStateV1::new(
                self.initial.thermodynamics(),
                &amounts,
                self.initial.volume(),
                self.initial.heavy_temperature(),
                Temperature::new(point.te),
                self.controls.local.relative_charge_tolerance,
                &mut *cancelled,
            )
            .map_err(|error| match error {
                FrozenIonizedGasStateErrorV1::Cancelled => IonizationNetworkError::Cancelled,
                error => IonizationNetworkError::State(error),
            })?
        };
        balance.absorbed_photons = NA * photo;
        balance.emitted_photons = NA * capture;
        balance.absorbed_radiation_energy_j = absorbed;
        balance.emitted_radiation_energy_j = emitted;
        balance.chemical_storage_change_j = chemical.value;
        balance.electron_sensible_energy_change_j = sensible;
        balance.heavy_sensible_energy_change_j = heavy;
        balance.energy_residual_j = residual.value;
        balance.relative_charge_residual = state.relative_charge_residual();
        Ok((state, channels))
    }
}

#[derive(Clone, Copy)]
struct Pair {
    lower: f64,
    upper: f64,
    electron: f64,
    te: f64,
}
#[derive(Clone, Copy)]
struct Rates {
    birth: f64,
    death: f64,
    capture: f64,
    captured_energy: f64,
}

fn photo_error(error: PhotoionizationErrorV1) -> IonizationNetworkError {
    match error {
        PhotoionizationErrorV1::Cancelled => IonizationNetworkError::Cancelled,
        error => IonizationNetworkError::Photoionization(error),
    }
}
