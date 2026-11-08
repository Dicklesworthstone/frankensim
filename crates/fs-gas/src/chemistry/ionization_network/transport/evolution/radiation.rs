//! Radiation modes share the original network time controller and storage audit.

use super::super::diffuse::field::PreparedDiffuseField;
use super::super::diffuse::{NetworkDiffuseCellExposure, NetworkDiffuseSolve};
use super::*;
use crate::chemistry::diffuse_recombination::{
    DiffuseRecombinationFlux, DiffuseRecombinationSolve,
};

pub(super) enum RadiationField<'a> {
    Beams(field::PreparedField),
    Diffuse {
        field: PreparedDiffuseField<'a>,
        controls: NetworkDiffuseSolve,
        work: usize,
        evaluations: usize,
        maxima: [f64; 3],
    },
}

impl<'a> RadiationField<'a> {
    pub fn new(
        cells: &[NetworkRadiationCell<'a>],
        problems: &[Problem<'a>],
        area: Area,
        beams: &[DiffuseRecombinationBeam],
        controls: NetworkRadiationSolve,
        diffuse: Option<DiffuseRecombinationSolve>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if let Some(radiation) = diffuse {
            let controls = NetworkDiffuseSolve::new(controls, radiation);
            let field =
                PreparedDiffuseField::new(cells, problems, area, beams, controls, cancelled)?;
            let work = field.preparation_work();
            Ok(Self::Diffuse {
                field,
                controls,
                work,
                evaluations: 0,
                maxima: [0.0; 3],
            })
        } else {
            Ok(Self::Beams(field::PreparedField::new(
                cells, problems, area, beams, controls, cancelled,
            )?))
        }
    }

    pub fn is_diffuse(&self) -> bool {
        matches!(self, Self::Diffuse { .. })
    }

    pub fn illuminated(&self) -> bool {
        match self {
            Self::Beams(field) => field.illuminated(),
            Self::Diffuse { field, .. } => field.illuminated(),
        }
    }

    pub fn radiation_enabled(&self) -> bool {
        self.is_diffuse() || self.illuminated()
    }

    pub fn evaluate(
        &mut self,
        points: &[Point],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        match self {
            Self::Beams(field) => field.evaluate(points, cancelled),
            Self::Diffuse {
                field,
                controls,
                work,
                evaluations,
                maxima,
            } => {
                if *evaluations >= controls.network.maximum_radiation_evaluations {
                    return Err(Error::WorkBudget);
                }
                let remaining = controls
                    .network
                    .maximum_radiation_work
                    .checked_sub(*work)
                    .ok_or(Error::WorkBudget)?;
                *evaluations += 1;
                let result = field.evaluate(points, remaining, cancelled);
                *work = work
                    .checked_add(field.work_used())
                    .ok_or(Error::WorkBudget)?;
                result?;
                if *work > controls.network.maximum_radiation_work {
                    return Err(Error::WorkBudget);
                }
                for (maximum, value) in maxima.iter_mut().zip([
                    field.balance.maximum_normalized_spectral_estimate,
                    field.balance.maximum_normalized_angular_estimate,
                    field.balance.maximum_normalized_tail_bound,
                ]) {
                    *maximum = maximum.max(value);
                }
                Ok(())
            }
        }
    }

    pub fn work(&self) -> usize {
        match self {
            Self::Beams(field) => field.work,
            Self::Diffuse { work, .. } => *work,
        }
    }

    pub fn evaluations(&self) -> usize {
        match self {
            Self::Beams(field) => field.evaluations,
            Self::Diffuse { evaluations, .. } => *evaluations,
        }
    }

    pub fn diffuse_controls(&self) -> Option<NetworkDiffuseSolve> {
        match self {
            Self::Beams(_) => None,
            Self::Diffuse { controls, .. } => Some(*controls),
        }
    }

    pub fn maxima(&self) -> [f64; 3] {
        match self {
            Self::Beams(_) => [0.0; 3],
            Self::Diffuse { maxima, .. } => *maxima,
        }
    }

    pub fn beam_rates(&self) -> &[NetworkRadiationBeamRate] {
        match self {
            Self::Beams(field) => &field.beam_rates,
            Self::Diffuse { field, .. } => &field.beam_rates,
        }
    }

    pub fn cell_beam_rates(&self, cell: usize) -> &[NetworkRadiationBeamRate] {
        match self {
            Self::Beams(field) => field.cell_beam_rates(cell),
            Self::Diffuse { field, .. } => field.cell_beam_rates(cell),
        }
    }

    pub fn interfaces(&self) -> &[DiffuseRecombinationFlux] {
        match self {
            Self::Beams(_) => &[],
            Self::Diffuse { field, .. } => &field.interfaces,
        }
    }

    /// Fill the existing reusable edge increments and return independent shared
    /// cell totals when channel rounding cannot preserve small trapped escape.
    pub fn increments(
        &self,
        cell: usize,
        dt: f64,
        increments: &mut [DiffuseIncrement],
    ) -> Result<Option<[Dd; 3]>, Error> {
        match self {
            Self::Beams(field) => {
                for (increment, rate) in increments.iter_mut().zip(field.channel_rates(cell)) {
                    let extent = product(rate.ionization_rate_mol_per_s, dt)?;
                    let energy = product(rate.absorbed_power_w, dt)?;
                    *increment = DiffuseIncrement {
                        net_extent_mol: extent,
                        net_electron_energy_j: product(rate.electron_heating_power_w, dt)?,
                        net_gas_energy_j: energy,
                        absorbed_mol: extent,
                        captured_mol: 0.0,
                        absorbed_energy_j: energy,
                        captured_energy_j: 0.0,
                    };
                }
                Ok(None)
            }
            Self::Diffuse { field, .. } => {
                for (increment, rate) in increments.iter_mut().zip(field.channel_rates(cell)) {
                    let mut absorbed = Dd::ZERO;
                    for value in [
                        rate.self_absorbed_photons_per_s,
                        rate.cross_channel_absorbed_photons_per_s,
                        rate.nonlocal_absorbed_photons_per_s,
                        rate.incident_absorbed_photons_per_s,
                    ] {
                        absorbed = sum(absorbed, Dd::from_f64(value))?;
                    }
                    *increment = DiffuseIncrement {
                        net_extent_mol: signed_product(rate.net_ionization_rate_mol_per_s, dt)?,
                        net_electron_energy_j: signed_product(rate.net_electron_power_w, dt)?,
                        net_gas_energy_j: signed_product(rate.net_gas_power_w, dt)?,
                        absorbed_mol: product(absorbed.to_f64() / NA, dt)?,
                        captured_mol: product(rate.emitted_photons_per_s / NA, dt)?,
                        absorbed_energy_j: product(rate.absorbed_power_w, dt)?,
                        captured_energy_j: product(rate.captured_electron_power_w, dt)?,
                    };
                }
                let totals = field
                    .cell_net_rates(cell)
                    .map(|rate| rate * Dd::from_f64(dt));
                if totals
                    .iter()
                    .any(|v| !v.hi.is_finite() || !v.lo.is_finite())
                {
                    return Err(Error::Arithmetic("shared radiation increment"));
                }
                Ok(Some(totals))
            }
        }
    }

    pub fn diffuse_cell_rates(&self, cell: usize) -> Result<[Dd; 10], Error> {
        let Self::Diffuse { field, .. } = self else {
            return Err(Error::Arithmetic(
                "diffuse accounts without diffuse radiation",
            ));
        };
        let b = field.cell_balance(cell);
        let mut rates = [Dd::ZERO; 10];
        rates[0] = Dd::from_f64(b.emitted_photons_per_s);
        rates[5] = Dd::from_f64(b.emitted_power_w);
        rates[6] = Dd::from_f64(b.captured_electron_power_w);
        for channel in field.channel_rates(cell) {
            for (target, value) in rates[1..5].iter_mut().zip([
                channel.self_absorbed_photons_per_s,
                channel.cross_channel_absorbed_photons_per_s,
                channel.nonlocal_absorbed_photons_per_s,
                channel.diffuse_absorbed_power_w,
            ]) {
                *target = sum(*target, Dd::from_f64(value))?;
            }
        }
        rates[7..10].copy_from_slice(&field.cell_net_rates(cell));
        Ok(rates)
    }
}

pub(super) fn exposure(values: [Dd; 10]) -> NetworkDiffuseCellExposure {
    let a = values.map(Dd::to_f64);
    NetworkDiffuseCellExposure {
        emitted_photons: a[0],
        self_absorbed_photons: a[1],
        cross_channel_absorbed_photons: a[2],
        nonlocal_absorbed_photons: a[3],
        diffuse_absorbed_energy_j: a[4],
        emitted_energy_j: a[5],
        captured_electron_energy_j: a[6],
        net_ionization_extent_mol: a[7],
        net_electron_energy_j: a[8],
        net_gas_energy_j: a[9],
    }
}

fn signed_product(a: f64, b: f64) -> Result<f64, Error> {
    let value = a * b;
    if !value.is_finite() || (a != 0.0 && b != 0.0 && value == 0.0) {
        Err(Error::Arithmetic("signed radiation exposure increment"))
    } else {
        Ok(value)
    }
}
