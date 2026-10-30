//! Finite rigid thermochemical storage at one original NSF reservoir opening.
//!
//! The accepted chamber pressure and temperature are held through one complete
//! fluid/acoustic interval. Its actual streamed mass and conservative total
//! energy then update the same source-owned gas inventory. This is an explicit
//! first-order partition: interval refinement remains necessary. No nozzle law,
//! mean transport into the acoustic bore, or separate acoustic work is added.

use super::*;
use fs_lbm::nsf3::{NsfGasMedium, NsfReservoirExchange};
use fs_material::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, parcel::ClosedGasParcel,
};
use fs_qty::Temperature;

/// Explicit inventory-movement and conservation limits for the held-source step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfFiniteReservoirLimits {
    /// Sum of absolute original substep mass transfers divided by starting mass.
    /// Strictly between zero and one; opposing transfers cannot hide movement.
    pub maximum_relative_mass_exchange: f64,
    /// Sum of absolute original substep energy transfers divided by the
    /// magnitude of starting sensible-chart energy. Strictly between zero and one.
    pub maximum_relative_energy_exchange: f64,
    /// Absolute allowance for actual inventory change minus streamed mass [kg].
    pub mass_absolute_tolerance_kg: f64,
    /// Absolute allowance for original-law sensible storage minus NSF energy [J].
    pub energy_absolute_tolerance_j: f64,
    /// Relative allowance on the actual and requested changes, in (0, 1e-3].
    /// Formation energy and the complete starting inventory do not scale residuals.
    pub relative_tolerance: f64,
}

impl NsfFiniteReservoirLimits {
    pub fn validate(&self) -> Result<(), NsfDgError> {
        if ![
            self.maximum_relative_mass_exchange,
            self.maximum_relative_energy_exchange,
            self.mass_absolute_tolerance_kg,
            self.energy_absolute_tolerance_j,
            self.relative_tolerance,
        ].into_iter().all(f64::is_finite)
            || !(0.0 < self.maximum_relative_mass_exchange
                && self.maximum_relative_mass_exchange < 1.0)
            || !(0.0 < self.maximum_relative_energy_exchange
                && self.maximum_relative_energy_exchange < 1.0)
            || self.mass_absolute_tolerance_kg <= 0.0
            || self.energy_absolute_tolerance_j <= 0.0
            || !(0.0 < self.relative_tolerance && self.relative_tolerance <= 1e-3)
        {
            return Err(NsfDgError::Input(
                "finite reservoir requires relative exchange bounds in (0,1), positive absolute balances and relative tolerance in (0,1e-3]",
            ));
        }
        Ok(())
    }
}

/// Accepted chamber storage and the original source exchange that changed it.
/// Positive source mass/energy leaves NSF and enters this rigid chamber.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfFiniteReservoirFrame {
    pub aperture: usize,
    pub volume_m3: f64,
    pub initial_mass_kg: f64,
    pub final_mass_kg: f64,
    pub mass_change_kg: f64,
    pub mass_residual_kg: f64,
    pub initial_pressure_pa: f64,
    pub final_pressure_pa: f64,
    pub initial_temperature_k: f64,
    pub final_temperature_k: f64,
    /// Storage in the same fixed sensible chart as the NSF total energy.
    pub initial_sensible_energy_j: f64,
    pub final_sensible_energy_j: f64,
    /// Original NASA internal-energy increment plus actual mass change times
    /// original sensible specific energy; no absolute formation subtraction.
    pub sensible_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub energy_reference_j_kg: f64,
    /// Includes pressure work, advection, conduction and viscous work already.
    /// Acoustic pressure work is a diagnostic subset, never a second debit.
    pub source_exchange: NsfReservoirExchange,
    pub caloric_iterations: usize,
}

pub(super) struct FiniteReservoir {
    pub(super) aperture: usize,
    gas: ClosedGasParcel,
    limits: NsfFiniteReservoirLimits,
    energy_solve: GasEnergySolve,
    energy_reference_j_kg: f64,
}

pub(super) struct Candidate {
    gas: ClosedGasParcel,
    frame: NsfFiniteReservoirFrame,
}

fn material(error: GasEnergyError) -> NsfDgError {
    if error == GasEnergyError::Cancelled {
        NsfDgError::Cancelled
    } else {
        NsfDgError::FiniteReservoir(error)
    }
}

impl<L: Lattice3<27>> NsfDg<L> {
    /// Attach one finite rigid source to an existing quiescent transverse opening.
    ///
    /// Its initial density and temperature come from that opening. The SAME
    /// thermochemical medium supplies pressure, composition, caloric reference
    /// and all transport properties. No legacy constant-cv replacement is used.
    /// Only FiniteReservoir drives are admitted while this owner is attached.
    /// A distinct opening may retain passive radiation in either binding order.
    ///
    /// Incoming resolved kinetic energy thermalizes in this well-mixed chamber.
    /// There is no stored chamber momentum, moving wall, pump or extra p*Q debit.
    pub fn with_finite_reservoir(
        mut self,
        aperture: usize,
        volume_m3: f64,
        limits: NsfFiniteReservoirLimits,
        gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        limits.validate()?;
        if self.steps != 0 || self.finite_reservoir.is_some()
            || !volume_m3.is_finite() || volume_m3 <= 0.0
            || self.reservoir_radiation_aperture() == Some(aperture)
        {
            return Err(NsfDgError::Input(
                "finite reservoir requires one cold positive volume and a source distinct from passive radiation",
            ));
        }
        let opening = self.zone.reservoir_aperture_at(aperture).ok_or(NsfDgError::Input(
            "finite reservoir requires an existing selected transverse opening",
        ))?;
        if opening.state.velocity_m_s != [0.0; 3] {
            return Err(NsfDgError::Input(
                "a quiescent finite rigid reservoir cannot prescribe ghost momentum",
            ));
        }
        let NsfGasMedium::Thermochemical(source) = self.zone.medium() else {
            return Err(NsfDgError::Input(
                "finite reservoir requires the original thermochemical NSF medium",
            ));
        };
        let point = source.context().clone().with("temperature", opening.state.temperature_k)
            .map_err(|_| NsfDgError::Input("finite reservoir temperature coordinate is invalid"))?;
        // A zero original-law energy increment reconstructs pressure from the
        // authored density and temperature without introducing another EOS.
        let solved = source.model().reconstruct_density_internal_energy_increment(
            &point, opening.state.density_kg_m3, 0.0, source.energy_solve(),
            || gate.is_requested(),
        ).map_err(material)?;
        let point = point.with("temperature", solved.state().state().temperature)
            .and_then(|p| p.with("pressure", solved.state().state().pressure))
            .map_err(|_| NsfDgError::Input("finite reservoir source coordinates are invalid"))?;
        let gas = ClosedGasParcel::new(source.model().clone(), point, volume_m3)
            .map_err(material)?;
        let authored_mass = opening.state.density_kg_m3 * volume_m3;
        let mass_residual = gas.mass_kg() - authored_mass;
        let sensible = gas.state().specific_internal_energy_j_kg()
            - source.energy_reference_j_kg();
        let storage = gas.mass_kg() * sensible;
        if ![authored_mass, mass_residual, sensible, storage].into_iter().all(f64::is_finite)
            || authored_mass <= 0.0 || storage == 0.0
            || mass_residual.abs() > limits.mass_absolute_tolerance_kg
        {
            return Err(NsfDgError::Input(
                "finite reservoir initial source inventory or sensible energy is unrepresentable",
            ));
        }
        poll(gate)?;
        self.finite_reservoir = Some(FiniteReservoir {
            aperture, gas, limits, energy_solve: source.energy_solve(),
            energy_reference_j_kg: source.energy_reference_j_kg(),
        });
        Ok(self)
    }

    /// Accepted original material inventory; failed trials leave it unchanged.
    pub fn finite_reservoir(&self) -> Option<&ClosedGasParcel> {
        self.finite_reservoir.as_ref().map(|owner| &owner.gas)
    }

    pub fn finite_reservoir_aperture(&self) -> Option<usize> {
        self.finite_reservoir.as_ref().map(|owner| owner.aperture)
    }

    pub fn finite_reservoir_limits(&self) -> Option<&NsfFiniteReservoirLimits> {
        self.finite_reservoir.as_ref().map(|owner| &owner.limits)
    }

    pub(super) fn admit_finite_reservoir_drive(
        &self, drive: NsfDgDrive,
    ) -> Result<(), NsfDgError> {
        match (self.finite_reservoir.as_ref(), drive) {
            (Some(owner), NsfDgDrive::FiniteReservoir { aperture })
                if owner.aperture == aperture => Ok(()),
            (None, NsfDgDrive::FiniteReservoir { .. }) => Err(NsfDgError::Input(
                "finite reservoir drive requires its attached inventory",
            )),
            (Some(_), _) => Err(NsfDgError::Input(
                "an attached finite reservoir exclusively owns its selected source control",
            )),
            (None, _) => Ok(()),
        }
    }

    pub(super) fn finite_reservoir_state(
        &self, aperture: usize,
    ) -> Result<NsfOpenState, NsfDgError> {
        self.admit_finite_reservoir_drive(NsfDgDrive::FiniteReservoir { aperture })?;
        let gas = self.finite_reservoir.as_ref()
            .expect("admitted finite reservoir").gas.state().state();
        self.zone.reservoir_aperture_state_for(aperture, gas.pressure, gas.temperature)
            .map_err(NsfDgError::Flow)
    }

    pub(super) fn step_with_finite_reservoir(
        &mut self, aperture: usize, gate: &CancelGate,
    ) -> Result<NsfDgFrame, NsfDgError> {
        poll(gate)?;
        let state = self.finite_reservoir_state(aperture)?;
        let drive = NsfDgDrive::FiniteReservoir { aperture };
        if self.upstream.is_some() {
            self.step_with_upstream_conditions(0.0, &[(aperture, state)], drive, gate)
        } else {
            self.step_with_boundary_conditions(
                NsfBoundaryDrive::IncomingCharacteristic(0.0), &[(aperture, state)], drive, gate,
            )
        }
    }
}

/// Refuse a held-source interval if ANY prefix permitted by the original total
/// variations can leave the caloric bracket or consume its admitted inventory.
/// All four energy/mass corners are checked; sensible energy may have either sign.
fn admit_envelope(
    owner: &FiniteReservoir, exchange: &NsfReservoirExchange, gate: &CancelGate,
) -> Result<(f64, f64), NsfDgError> {
    poll(gate)?;
    let old = &owner.gas;
    let mass = old.mass_kg();
    let sensible = old.state().specific_internal_energy_j_kg() - owner.energy_reference_j_kg;
    let energy = mass * sensible;
    let dm = exchange.mass_out_kg;
    let de = exchange.total_energy_out_j;
    let mv = exchange.mass_exchange_variation_kg;
    let ev = exchange.energy_exchange_variation_j;
    if ![mass, sensible, energy, dm, de, mv, ev].into_iter().all(f64::is_finite)
        || mass <= 0.0 || energy == 0.0 || mv < 0.0 || ev < 0.0
        || mv < dm.abs() * (1.0 - 64.0 * f64::EPSILON)
        || ev < de.abs() * (1.0 - 64.0 * f64::EPSILON)
        || mv > owner.limits.maximum_relative_mass_exchange * mass
        || ev > owner.limits.maximum_relative_energy_exchange * energy.abs()
    {
        return Err(NsfDgError::Input(
            "finite reservoir interval exceeds its original mass or sensible-energy movement bound",
        ));
    }
    if mv == 0.0 && ev == 0.0 && dm == 0.0 && de == 0.0 {
        return Ok((sensible, energy));
    }
    // Outward rounding encloses signed-reduction roundoff as well as every
    // chronological prefix. These are whole original substep transfers.
    let mass_variation = mv.max(dm.abs()).next_up();
    let energy_variation = ev.max(de.abs()).next_up();
    let low_mass = (mass - mass_variation).next_down();
    let high_mass = (mass + mass_variation).next_up();
    let low_energy = (energy - energy_variation).next_down();
    let high_energy = (energy + energy_variation).next_up();
    let temperature = old.state().state().temperature;
    let lower_change = old.model().thermodynamics().specific_internal_energy_change(
        Temperature::new(temperature),
        Temperature::new(owner.energy_solve.minimum_temperature_k),
    ).map_err(|_| NsfDgError::Input(
        "original caloric source refused the finite reservoir lower energy envelope",
    ))?.value();
    let upper_change = old.model().thermodynamics().specific_internal_energy_change(
        Temperature::new(temperature),
        Temperature::new(owner.energy_solve.maximum_temperature_k),
    ).map_err(|_| NsfDgError::Input(
        "original caloric source refused the finite reservoir upper energy envelope",
    ))?.value();
    let lower = sensible + lower_change;
    let upper = sensible + upper_change;
    if ![low_mass, high_mass, low_energy, high_energy, lower, upper]
        .into_iter().all(f64::is_finite)
        || low_mass <= 0.0 || lower >= upper
        || [low_energy, high_energy].into_iter().any(|e| {
            [low_mass, high_mass].into_iter().any(|m| {
                let value = e / m;
                !value.is_finite() || value.next_down() < lower || value.next_up() > upper
            })
        })
    {
        return Err(NsfDgError::Input(
            "finite reservoir substep envelope can leave its original caloric branch; refine the interval",
        ));
    }
    poll(gate)?;
    Ok((sensible, energy))
}

/// Prepare only after the original fluid/acoustic solution is admitted.
/// No accepted inventory, ghost, clock or acoustic history is changed here.
pub(super) fn prepare<L: Lattice3<27>>(
    owner: Option<&FiniteReservoir>,
    fluid: &NsfCharacteristicTrial<'_, L>,
    gate: &CancelGate,
) -> Result<Option<Candidate>, NsfDgError> {
    let Some(owner) = owner else { return Ok(None); };
    poll(gate)?;
    let NsfGasMedium::Thermochemical(source) = fluid.grid().medium() else {
        return Err(NsfDgError::Input("finite reservoir lost its original thermochemical medium"));
    };
    if owner.gas.model() != source.model()
        || owner.energy_reference_j_kg.to_bits() != source.energy_reference_j_kg().to_bits()
        || owner.energy_solve != source.energy_solve()
    {
        return Err(NsfDgError::Input("finite reservoir and NSF must retain identical source laws and energy reference"));
    }
    let exchange = *fluid.reservoir_exchange_at(owner.aperture).ok_or(NsfDgError::Input(
        "finite reservoir lost its original indexed conservative exchange",
    ))?;
    if exchange.aperture != owner.aperture {
        return Err(NsfDgError::Input("finite reservoir exchange changed its original owner"));
    }
    let (sensible, energy) = admit_envelope(owner, &exchange, gate)?;
    let old = &owner.gas;
    let dm = exchange.mass_out_kg;
    let de = exchange.total_energy_out_j;
    let target_mass = old.mass_kg() + dm;
    if !target_mass.is_finite() || target_mass <= 0.0
        || (dm != 0.0 && target_mass == old.mass_kg())
    {
        return Err(NsfDgError::Input("finite reservoir cannot exhaust or lose a streamed mass increment"));
    }
    let (gas, iterations) = if dm == 0.0 && de == 0.0 {
        (old.clone(), 0)
    } else {
        // Total FV energy already contains pressure/viscous work, advection and
        // heat. Only the mass-dependent chart conversion belongs in this solve.
        let thermal = de - dm * sensible;
        let change = thermal / target_mass;
        if !thermal.is_finite() || !change.is_finite()
            || (thermal != 0.0 && change == 0.0)
        {
            return Err(NsfDgError::Input("finite reservoir thermal increment is unrepresentable"));
        }
        let solved = old.model().reconstruct_density_internal_energy_increment(
            old.point(), target_mass / old.volume_m3(), change, owner.energy_solve,
            || gate.is_requested(),
        ).map_err(material)?;
        let point = old.point().clone().with("temperature", solved.state().state().temperature)
            .and_then(|p| p.with("pressure", solved.state().state().pressure))
            .map_err(|_| NsfDgError::Input("finite reservoir endpoint coordinates are invalid"))?;
        let next = ClosedGasParcel::new(old.model().clone(), point, old.volume_m3())
            .map_err(material)?;
        (next, solved.iterations())
    };
    let mass_change = gas.mass_kg() - old.mass_kg();
    let mass_residual = mass_change - dm;
    let du = old.model().thermodynamics().specific_internal_energy_change(
        Temperature::new(old.state().state().temperature),
        Temperature::new(gas.state().state().temperature),
    ).map_err(|_| NsfDgError::Input(
        "original caloric source refused the finite reservoir storage audit",
    ))?.value();
    let stored = gas.mass_kg() * du + mass_change * sensible;
    let final_energy = gas.mass_kg() * (sensible + du);
    let energy_residual = stored - de;
    let mass_allowance = owner.limits.mass_absolute_tolerance_kg
        + owner.limits.relative_tolerance * (mass_change.abs() + dm.abs());
    let energy_allowance = owner.limits.energy_absolute_tolerance_j
        + owner.limits.relative_tolerance * (stored.abs() + de.abs());
    if ![mass_change, mass_residual, du, stored, final_energy, energy_residual,
        mass_allowance, energy_allowance].into_iter().all(f64::is_finite)
        || mass_residual.abs() > mass_allowance
        || energy_residual.abs() > energy_allowance
        || (dm != 0.0 && (mass_change == 0.0
            || mass_change.is_sign_positive() != dm.is_sign_positive()))
        || (de != 0.0 && (stored == 0.0
            || stored.is_sign_positive() != de.is_sign_positive()))
    {
        return Err(NsfDgError::Input(
            "finite reservoir actual source inventory failed mass or sensible-energy conservation",
        ));
    }
    let frame = NsfFiniteReservoirFrame {
        aperture: owner.aperture,
        volume_m3: old.volume_m3(),
        initial_mass_kg: old.mass_kg(),
        final_mass_kg: gas.mass_kg(),
        mass_change_kg: mass_change,
        mass_residual_kg: mass_residual,
        initial_pressure_pa: old.state().state().pressure,
        final_pressure_pa: gas.state().state().pressure,
        initial_temperature_k: old.state().state().temperature,
        final_temperature_k: gas.state().state().temperature,
        initial_sensible_energy_j: energy,
        final_sensible_energy_j: final_energy,
        sensible_energy_change_j: stored,
        energy_residual_j: energy_residual,
        energy_reference_j_kg: owner.energy_reference_j_kg,
        source_exchange: exchange,
        caloric_iterations: iterations,
    };
    poll(gate)?;
    Ok(Some(Candidate { gas, frame }))
}

/// Infallible publication after the caller's final cancellation poll and audits.
pub(super) fn publish(
    owner: &mut Option<FiniteReservoir>, candidate: Option<Candidate>,
) -> Option<NsfFiniteReservoirFrame> {
    candidate.map(|candidate| {
        owner.as_mut().expect("prepared finite reservoir owner").gas = candidate.gas;
        candidate.frame
    })
}
