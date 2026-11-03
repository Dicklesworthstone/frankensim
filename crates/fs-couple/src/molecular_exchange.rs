//! Conservative molecular transfer between two ORIGINAL finite moving gases.
//!
//! One thin aperture joins two homogeneous, fixed-volume gas owners in a
//! declared common inertial frame. Both drifting Maxwellians are sampled in
//! the aperture frame; each actual particle stream is debited and credited
//! with the SAME species mass, vector momentum and complete crossing energy.
//! Original bulk kinetic changes are removed before one original caloric
//! reconstruction per gas. No second EOS or receiver-mixing heat is introduced.
//!
//! This is a frozen-rate, first-order control-volume exchange, not particle
//! flight between separated vessels, a moving-container model or a resolved
//! gas-flow solver. Other boundary pressure/support forces and moving-volume
//! work require their own owners. Equal counterflow still transfers normal
//! momentum; an aperture alone does not supply a whole vessel's wall balance.

use crate::liquid_gas::evaporation::momentum::{
    MovingGasParcel,
    boundary::{MovingGasBoundaryError, MovingGasBoundarySolve, MovingGasBoundaryStep},
};
use fs_material::thermochemical_gas::energy::{
    GasEnergyError,
    parcel::transfer::species::{
        boundary::SpeciesBoundaryRelativeTransfer,
        effusion::{
            MAX_EFFUSIVE_SPECIES,
            drift::{
                EffusiveDriftInterface,
                exchange::{self, DriftingEffusiveExchange},
            },
        },
    },
};
use fs_math::dd::Dd;
use fs_thermochem::SpeciesId;

/// One shared aperture and energy frame, in the SAME axes as both gas momenta.
/// The original gas state has no position or independently transformable frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularExchangeSurface {
    pub area_m2: f64,
    /// Unit direction from the left gas into the right gas.
    pub normal_left_to_right: [f64; 3],
    /// Inertial aperture/control-frame velocity, held over this interval [m/s].
    /// Zero area leaves this velocity and the source distributions inactive.
    pub velocity_world_m_s: [f64; 3],
}

/// Original per-gas solves and independent budgets for the joint transaction.
#[derive(Clone, Copy, Debug)]
pub struct MolecularExchangeControls {
    pub gas: [MovingGasBoundarySolve; 2],
    pub maximum_duration_s: f64,
    /// Maximum gross outgoing mass divided by each ORIGINAL species inventory.
    /// Both directions are admitted before either gas is reconstructed.
    pub maximum_species_fraction: f64,
    pub mass_absolute_kg: f64,
    /// Per-component allowance, independent of momentum on other axes.
    pub momentum_absolute_n_s: f64,
    pub energy_absolute_j: f64,
    /// Scales exchanged mass, each requested net impulse component, and the
    /// formation-independent relative energy/kinetic changes, never large
    /// original inventories, common drift or absolute formation energy.
    pub relative_tolerance: f64,
}

#[derive(Debug)]
pub enum MolecularExchangeError {
    Invalid(&'static str),
    Material(GasEnergyError),
    MovingGas {
        side: usize,
        error: MovingGasBoundaryError,
    },
    /// A bounded adaptive owner may retry this interval at a smaller duration.
    Duration {
        duration_s: f64,
        maximum_duration_s: f64,
    },
    /// Gross withdrawal cannot be funded by simultaneous incoming material.
    Withdrawal {
        side: usize,
        species: SpeciesId,
        fraction: f64,
        maximum: f64,
    },
    Balance {
        account: &'static str,
        residual: f64,
        allowance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MolecularExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite moving-gas molecular exchange refused: {self:?}")
    }
}
impl std::error::Error for MolecularExchangeError {}
impl From<GasEnergyError> for MolecularExchangeError {
    fn from(error: GasEnergyError) -> Self {
        if error == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Material(error)
        }
    }
}

/// The common streams and the TWO original finite-state receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularExchangeStep {
    pub duration_s: f64,
    pub surface: MolecularExchangeSurface,
    pub exchange: DriftingEffusiveExchange,
    pub gas: [MovingGasBoundaryStep; 2],
    /// Positive signed mass/impulse goes from left to right.
    pub net_mass_left_to_right_kg: f64,
    pub net_momentum_left_to_right_n_s: [f64; 3],
    /// Difference between rounded gross masses and the retained small-drift
    /// net coordinate. Gated on BOTH original gas mass budgets before solving.
    pub stream_mass_rounding_kg: f64,
    pub maximum_stream_mass_rounding_kg: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub momentum_residual_n_s: [f64; 3],
    /// Canonical reciprocal-stream cancellation, with formation terms canceled.
    pub stream_energy_residual_j: f64,
    /// Joint internal-plus-bulk energy closure in the aperture frame. Original
    /// source-reference residuals are used; inventory rounding is bounded by
    /// the separate actual total and individual-species mass checks.
    pub energy_residual_j: f64,
    /// The same account in the retained momenta's inertial frame, including
    /// V_aperture dot momentum_residual - |V_aperture|² mass_residual/2.
    pub laboratory_energy_residual_j: f64,
    pub mass_allowance_kg: f64,
    pub momentum_allowance_n_s: [f64; 3],
    pub energy_allowance_j: f64,
}

/// Both finite endpoints escape together, after original and joint admission.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularExchangeCandidate {
    gases: [MovingGasParcel; 2],
    step: MolecularExchangeStep,
}
impl MolecularExchangeCandidate {
    pub const fn gases(&self) -> &[MovingGasParcel; 2] {
        &self.gases
    }
    pub const fn step(&self) -> &MolecularExchangeStep {
        &self.step
    }
    pub fn into_parts(self) -> ([MovingGasParcel; 2], MolecularExchangeStep) {
        (self.gases, self.step)
    }
}

type Error = MolecularExchangeError;

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: Dd) -> Result<f64, Error> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value.to_f64())
    } else {
        Err(Error::Invalid(
            "molecular exchange coordinate is unrepresentable",
        ))
    }
}
fn allowance(absolute: f64, relative: f64, scale: Dd) -> Result<f64, Error> {
    let result = finite(Dd::from_f64(absolute) + Dd::from_f64(relative) * scale.abs())?;
    if positive(result) {
        Ok(result)
    } else {
        Err(Error::Invalid(
            "joint molecular allowance is unrepresentable",
        ))
    }
}
fn balance(account: &'static str, residual: Dd, allowance: f64) -> Result<f64, Error> {
    let value = finite(residual)?;
    if Dd::from_f64(allowance).lt(residual.abs()) {
        Err(Error::Balance {
            account,
            residual: value,
            allowance,
        })
    } else {
        Ok(value)
    }
}
fn moving_error(side: usize, error: MovingGasBoundaryError) -> Error {
    match error {
        MovingGasBoundaryError::Cancelled
        | MovingGasBoundaryError::Gas(GasEnergyError::Cancelled) => Error::Cancelled,
        error => Error::MovingGas { side, error },
    }
}

/// Transfer the SAME original molecular beams between two finite gas owners.
///
/// Frozen gas temperatures, drift distributions and aperture data make this
/// first order in time; accepted endpoint source/drift admission is mandatory
/// but is not a continuous trajectory enclosure. Resample on subsequent steps.
/// Only duration and gross-withdrawal refusals explicitly request refinement;
/// source, caloric and conservation refusals retain their physical meaning.
///
/// Fixed gas volumes and all original auxiliary source coordinates survive.
/// No wall-pressure force, container trajectory, mixing-heat patch or finite
/// flight path is inferred. The scalar and vector accounts describe this one
/// shared molecular interface; all other loads need explicit composing owners.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    gases: [&MovingGasParcel; 2],
    interface: &EffusiveDriftInterface,
    surface: MolecularExchangeSurface,
    duration_s: f64,
    controls: MolecularExchangeControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MolecularExchangeCandidate, Error> {
    poll(&mut cancelled)?;
    if !positive(duration_s)
        || !positive(controls.maximum_duration_s)
        || !positive(controls.maximum_species_fraction)
        || controls.maximum_species_fraction >= 1.0
        || !positive(controls.mass_absolute_kg)
        || !positive(controls.momentum_absolute_n_s)
        || !positive(controls.energy_absolute_j)
        || !controls.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_tolerance)
        || !surface.area_m2.is_finite()
        || surface.area_m2 < 0.0
    {
        return Err(Error::Invalid(
            "positive physical interval/budgets, bounded original withdrawal and nonnegative area required",
        ));
    }
    for solve in controls.gas {
        solve.gas.gas.validate()?;
        if !positive(solve.gas.mass_tolerance_kg)
            || !positive(solve.gas.energy_tolerance_j)
            || !positive(solve.momentum_absolute_n_s)
            || !solve.momentum_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&solve.momentum_relative_tolerance)
        {
            return Err(Error::Invalid(
                "positive original moving-gas budgets required",
            ));
        }
    }
    if duration_s > controls.maximum_duration_s {
        return Err(Error::Duration {
            duration_s,
            maximum_duration_s: controls.maximum_duration_s,
        });
    }
    let active = surface.area_m2 > 0.0;
    let frame = if active {
        surface.velocity_world_m_s
    } else {
        [0.0; 3]
    };
    if frame.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid(
            "active aperture needs a finite velocity in the common gas frame",
        ));
    }
    let relative = if active {
        gases.map(|gas| {
            let velocity = gas.velocity_m_s();
            core::array::from_fn(|d| velocity[d] - frame[d])
        })
    } else {
        [[0.0; 3]; 2]
    };
    let exchange = exchange::sample(
        gases[0].parcel(),
        gases[1].parcel(),
        interface,
        surface.area_m2,
        relative,
        surface.normal_left_to_right,
        duration_s,
        &mut cancelled,
    )?;
    if exchange.species.len() > MAX_EFFUSIVE_SPECIES {
        return Err(Error::Invalid(
            "paired molecular source exceeds the original species bound",
        ));
    }
    let catalog = interface.momentum_interface().interface().gas_catalog();
    let mass_rounding_allowance = controls.gas[0]
        .gas
        .mass_tolerance_kg
        .min(controls.gas[1].gas.mass_tolerance_kg);
    let mut rows = [
        Vec::with_capacity(exchange.species.len()),
        Vec::with_capacity(exchange.species.len()),
    ];
    let mut net_mass = Dd::ZERO;
    let mut gross_mass = Dd::ZERO;
    let mut stream_mass_rounding = Dd::ZERO;
    let mut maximum_stream_mass_rounding = 0.0_f64;
    let mut impulse = [Dd::ZERO; 3];
    let mut stream_energy_residual = Dd::ZERO;
    let mut energy_scale = Dd::ZERO;
    let temperatures = gases.map(|gas| gas.parcel().state().state().temperature);
    for row in &exchange.species {
        poll(&mut cancelled)?;
        for side in 0..2 {
            let mass = row.mass_kg[side];
            if mass > 0.0 {
                let available = gases[side].parcel().species_mass_kg(&row.species)?;
                if !positive(available) {
                    return Err(Error::Invalid(
                        "outgoing molecular species lost its original inventory",
                    ));
                }
                if (Dd::from_f64(available) * Dd::from_f64(controls.maximum_species_fraction))
                    .lt(Dd::from_f64(mass))
                {
                    return Err(Error::Withdrawal {
                        side,
                        species: row.species.clone(),
                        fraction: mass / available,
                        maximum: controls.maximum_species_fraction,
                    });
                }
            }
        }
        let transferred = Dd::from_f64(row.net_mass_left_to_right_kg);
        let rounded = Dd::from_f64(row.mass_kg[0]) - Dd::from_f64(row.mass_kg[1]) - transferred;
        maximum_stream_mass_rounding = maximum_stream_mass_rounding.max(
            balance(
                "individual gross/net molecular mass",
                rounded,
                mass_rounding_allowance,
            )?
            .abs(),
        );
        stream_mass_rounding = stream_mass_rounding + rounded;
        net_mass = net_mass + transferred;
        gross_mass = gross_mass + Dd::from_f64(row.mass_kg[0]) + Dd::from_f64(row.mass_kg[1]);
        for (d, component) in impulse.iter_mut().enumerate() {
            *component = *component + Dd::from_f64(row.momentum_relative_kg_m_s[0][d])
                - Dd::from_f64(row.momentum_relative_kg_m_s[1][d]);
        }
        for side in 0..2 {
            let relative_energy = row.relative_energy_into_gas_j[side];
            rows[side].push(SpeciesBoundaryRelativeTransfer {
                species: row.species.clone(),
                mass_into_gas_kg: if side == 0 {
                    -row.net_mass_left_to_right_kg
                } else {
                    row.net_mass_left_to_right_kg
                },
                relative_energy_into_gas_j: relative_energy,
            });
            energy_scale = energy_scale + Dd::from_f64(relative_energy).abs();
        }
        // Both reference enthalpies use the SAME original source. This check
        // removes formation constants analytically rather than subtracting
        // the displayed absolute energies of the two gross streams.
        stream_energy_residual = stream_energy_residual
            + Dd::from_f64(row.relative_energy_into_gas_j[0])
            + Dd::from_f64(row.relative_energy_into_gas_j[1]);
        if row.net_mass_left_to_right_kg != 0.0 {
            let dh = catalog.species_specific_enthalpy_change_j_kg(
                &row.species,
                temperatures[0],
                temperatures[1],
            )?;
            stream_energy_residual = stream_energy_residual + transferred * Dd::from_f64(dh);
        }
    }
    let stream_mass_rounding_kg = balance(
        "total gross/net molecular mass",
        stream_mass_rounding,
        mass_rounding_allowance,
    )?;
    let net_mass_left_to_right_kg = finite(net_mass)?;
    // Transform the already compensated NET species mass once. Large opposite
    // species streams need not form and subtract their common drift impulses.
    for (d, component) in impulse.iter_mut().enumerate() {
        *component = *component + Dd::from_f64(frame[d]) * net_mass;
    }
    let net_momentum_left_to_right_n_s = [
        finite(impulse[0])?,
        finite(impulse[1])?,
        finite(impulse[2])?,
    ];
    let left = gases[0]
        .exchange_species_momentum_relative(
            if active {
                catalog
            } else {
                gases[0].parcel().model()
            },
            &rows[0],
            net_momentum_left_to_right_n_s.map(|p| -p),
            frame,
            controls.gas[0],
            &mut cancelled,
        )
        .map_err(|error| moving_error(0, error))?;
    poll(&mut cancelled)?;
    let right = gases[1]
        .exchange_species_momentum_relative(
            if active {
                catalog
            } else {
                gases[1].parcel().model()
            },
            &rows[1],
            net_momentum_left_to_right_n_s,
            frame,
            controls.gas[1],
            &mut cancelled,
        )
        .map_err(|error| moving_error(1, error))?;
    let (left, _, left_step) = left.into_parts();
    let (right, _, right_step) = right.into_parts();
    let next = [left, right];
    let gas = [left_step, right_step];
    let final_relative = if active {
        core::array::from_fn(|side| {
            let velocity = next[side].velocity_m_s();
            core::array::from_fn(|d| velocity[d] - frame[d])
        })
    } else {
        [[0.0; 3]; 2]
    };
    exchange.admit_endpoints(
        [next[0].parcel(), next[1].parcel()],
        final_relative,
        &mut cancelled,
    )?;
    let mut mass_residual = Dd::ZERO;
    let mut momentum_residual = [Dd::ZERO; 3];
    let mut energy_residual = stream_energy_residual;
    for side in 0..2 {
        poll(&mut cancelled)?;
        mass_residual = mass_residual + Dd::from_f64(next[side].parcel().mass_kg())
            - Dd::from_f64(gases[side].parcel().mass_kg());
        for (d, residual) in momentum_residual.iter_mut().enumerate() {
            *residual = *residual + Dd::from_f64(next[side].momentum_kg_m_s()[d])
                - Dd::from_f64(gases[side].momentum_kg_m_s()[d]);
        }
        energy_residual = energy_residual + Dd::from_f64(gas[side].total_energy_residual_j);
        energy_scale = energy_scale + Dd::from_f64(gas[side].kinetic_energy_change_j).abs();
    }
    let mass_allowance_kg = allowance(
        controls.mass_absolute_kg,
        controls.relative_tolerance,
        gross_mass,
    )?;
    let mass_residual_kg = balance("joint actual gas mass", mass_residual, mass_allowance_kg)?;
    let mut maximum_species_mass_residual_kg = 0.0_f64;
    for row in &exchange.species {
        poll(&mut cancelled)?;
        let mut residual = Dd::ZERO;
        for side in 0..2 {
            residual = residual + Dd::from_f64(next[side].parcel().species_mass_kg(&row.species)?)
                - Dd::from_f64(gases[side].parcel().species_mass_kg(&row.species)?);
        }
        let limit = allowance(
            controls.mass_absolute_kg,
            controls.relative_tolerance,
            Dd::from_f64(row.mass_kg[0]) + Dd::from_f64(row.mass_kg[1]),
        )?;
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg
            .max(balance("joint actual species mass", residual, limit)?.abs());
    }
    let mut momentum_residual_n_s = [0.0; 3];
    let mut momentum_allowance_n_s = [0.0; 3];
    let mut frame_conversion = Dd::ZERO;
    for d in 0..3 {
        momentum_allowance_n_s[d] = allowance(
            controls.momentum_absolute_n_s,
            controls.relative_tolerance,
            impulse[d],
        )?;
        momentum_residual_n_s[d] = balance(
            "joint actual gas momentum",
            momentum_residual[d],
            momentum_allowance_n_s[d],
        )?;
        let velocity = Dd::from_f64(frame[d]);
        frame_conversion = frame_conversion + velocity * momentum_residual[d]
            - (mass_residual * Dd::from_f64(0.5) * velocity) * velocity;
    }
    let energy_allowance_j = allowance(
        controls.energy_absolute_j,
        controls.relative_tolerance,
        energy_scale,
    )?;
    let stream_energy_residual_j = balance(
        "reciprocal molecular energy",
        stream_energy_residual,
        energy_allowance_j,
    )?;
    let energy_residual_j = balance(
        "joint aperture-frame total energy",
        energy_residual,
        energy_allowance_j,
    )?;
    let laboratory_energy_residual_j = balance(
        "joint inertial-frame total energy",
        energy_residual + frame_conversion,
        energy_allowance_j,
    )?;
    poll(&mut cancelled)?;
    Ok(MolecularExchangeCandidate {
        gases: next,
        step: MolecularExchangeStep {
            duration_s,
            surface,
            exchange,
            gas,
            net_mass_left_to_right_kg,
            net_momentum_left_to_right_n_s,
            stream_mass_rounding_kg,
            maximum_stream_mass_rounding_kg: maximum_stream_mass_rounding,
            mass_residual_kg,
            maximum_species_mass_residual_kg,
            momentum_residual_n_s,
            stream_energy_residual_j,
            energy_residual_j,
            laboratory_energy_residual_j,
            mass_allowance_kg,
            momentum_allowance_n_s,
            energy_allowance_j,
        },
    })
}

#[cfg(test)]
mod tests;
