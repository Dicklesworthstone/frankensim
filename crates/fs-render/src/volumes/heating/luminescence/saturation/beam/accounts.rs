//! Accepted spectral accounts; shared kinetic work is counted once.

use super::{BeamPath, CheckedSum, ContentHash, OpticalHeatingError, invalid};
use fs_thermal::luminescence::transfer::pump_beams::{PumpBeamStep, PumpBeamsStep};
use fs_thermal::luminescence::transfer::pump_column::PumpColumnStep;

/// Original single-beam account, retained for the one-column interface.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DepletingPumpEnergy {
    pub duration_s: f64,
    pub incident_j: [f64; 2],
    pub transmitted_j: [f64; 2],
    pub absorbed_j: [f64; 2],
    pub passive_absorbed_j: [f64; 2],
    pub residual_j: [f64; 2],
    pub accepted_kinetic_substeps: u64,
    pub cell_evaluations: u64,
    /// Temporal convergence indicator, not a certified continuum error bound.
    pub maximum_energy_discrepancy_j: f64,
}

/// One declared beam's accepted photon-energy budget, in its own wavelength
/// order. Different beams' lane indices need not represent the same wavelength.
/// Repeated source identities are permitted for independently declared rays.
#[derive(Clone, Debug, PartialEq)]
pub struct DepletingPumpBeamEnergy {
    pub source_identity: ContentHash,
    pub wavelengths_nm: [f64; 2],
    pub incident_j: [f64; 2],
    pub transmitted_j: [f64; 2],
    pub absorbed_j: [f64; 2],
    pub passive_absorbed_j: [f64; 2],
    /// Incident minus transmitted, molecular and passive absorption.
    pub residual_j: [f64; 2],
}

/// Accepted accounts for the simultaneous physical pump system. A shared
/// cell's kinetic work is counted here once, never once per illuminating beam.
#[derive(Clone, Debug, PartialEq)]
pub struct DepletingPumpReport {
    pub duration_s: f64,
    /// Original declaration order, including inactive and missed beams.
    pub beams: Vec<DepletingPumpBeamEnergy>,
    pub accepted_kinetic_substeps: u64,
    pub cell_evaluations: u64,
    pub coupling_iterations: u64,
    /// Largest accepted whole-system temporal-refinement discrepancy, joules.
    /// This is an estimator, not a certified spatial or experimental bound.
    pub maximum_energy_discrepancy_j: f64,
}

#[derive(Clone, Debug, Default)]
struct SpectralAccounts {
    incident: [CheckedSum; 2],
    transmitted: [CheckedSum; 2],
    absorbed: [CheckedSum; 2],
    passive: [CheckedSum; 2],
}

#[derive(Clone, Debug)]
pub(in crate::volumes::heating::luminescence::saturation) struct BeamAccounts {
    duration: CheckedSum,
    spectra: Vec<SpectralAccounts>,
    report: DepletingPumpReport,
    single: Option<DepletingPumpEnergy>,
}

impl BeamAccounts {
    pub(super) fn new(paths: &[BeamPath]) -> Self {
        Self {
            duration: CheckedSum::default(),
            spectra: vec![SpectralAccounts::default(); paths.len()],
            report: DepletingPumpReport {
                duration_s: 0.0,
                beams: paths
                    .iter()
                    .map(|path| DepletingPumpBeamEnergy {
                        source_identity: path.source_identity,
                        wavelengths_nm: path.wavelengths_nm,
                        incident_j: [0.0; 2],
                        transmitted_j: [0.0; 2],
                        absorbed_j: [0.0; 2],
                        passive_absorbed_j: [0.0; 2],
                        residual_j: [0.0; 2],
                    })
                    .collect(),
                accepted_kinetic_substeps: 0,
                cell_evaluations: 0,
                coupling_iterations: 0,
                maximum_energy_discrepancy_j: 0.0,
            },
            single: (paths.len() == 1).then(DepletingPumpEnergy::default),
        }
    }

    pub(super) fn report(&self) -> &DepletingPumpReport {
        &self.report
    }
    pub(super) fn single(&self) -> Option<&DepletingPumpEnergy> {
        self.single.as_ref()
    }

    pub(super) fn add_column(
        &mut self,
        step: &PumpColumnStep<2>,
        duration_s: f64,
        energy_tolerance_j: f64,
    ) -> Result<(), OpticalHeatingError> {
        if self.spectra.len() != 1 {
            return Err(invalid("single column differs from declared pump system"));
        }
        self.add_spectrum(
            0,
            [
                step.incident_energy_j,
                step.transmitted_energy_j,
                step.absorbed_energy_j,
                step.passive_absorbed_energy_j,
            ],
            energy_tolerance_j,
        )?;
        self.add_work(
            duration_s,
            step.substeps,
            step.cell_evaluations,
            0,
            step.max_energy_discrepancy_j,
        )
    }

    pub(super) fn add_beams(
        &mut self,
        step: &PumpBeamsStep<2>,
        duration_s: f64,
        energy_tolerance_j: f64,
    ) -> Result<(), OpticalHeatingError> {
        self.add_paths(
            &step.beams,
            duration_s,
            step.substeps,
            step.cell_evaluations,
            step.coupling_iterations,
            step.max_energy_discrepancy_j,
            energy_tolerance_j,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn add_paths(
        &mut self,
        paths: &[PumpBeamStep<2>],
        duration_s: f64,
        substeps: usize,
        cell_evaluations: usize,
        coupling_iterations: usize,
        discrepancy_j: f64,
        energy_tolerance_j: f64,
    ) -> Result<(), OpticalHeatingError> {
        if paths.len() != self.spectra.len() {
            return Err(invalid("coupled result differs from declared pump beams"));
        }
        for (index, beam) in paths.iter().enumerate() {
            self.add_spectrum(
                index,
                [
                    beam.incident_energy_j,
                    beam.transmitted_energy_j,
                    beam.absorbed_energy_j,
                    beam.passive_absorbed_energy_j,
                ],
                energy_tolerance_j,
            )?;
        }
        self.add_work(
            duration_s,
            substeps,
            cell_evaluations,
            coupling_iterations,
            discrepancy_j,
        )
    }

    fn add_spectrum(
        &mut self,
        index: usize,
        energies: [[f64; 2]; 4],
        tolerance: f64,
    ) -> Result<(), OpticalHeatingError> {
        let account = &mut self.spectra[index];
        let report = &mut self.report.beams[index];
        for lane in 0..2 {
            account.incident[lane].add(energies[0][lane])?;
            account.transmitted[lane].add(energies[1][lane])?;
            account.absorbed[lane].add(energies[2][lane])?;
            account.passive[lane].add(energies[3][lane])?;
            report.incident_j[lane] = account.incident[lane].value()?;
            report.transmitted_j[lane] = account.transmitted[lane].value()?;
            report.absorbed_j[lane] = account.absorbed[lane].value()?;
            report.passive_absorbed_j[lane] = account.passive[lane].value()?;
            let mut balance = CheckedSum::default();
            for value in [
                report.incident_j[lane],
                -report.transmitted_j[lane],
                -report.absorbed_j[lane],
                -report.passive_absorbed_j[lane],
            ] {
                balance.add(value)?;
            }
            report.residual_j[lane] = balance.value()?;
            if report.residual_j[lane].abs() > tolerance {
                return Err(invalid(
                    "one pump beam's incident, transmitted and absorbed energies do not balance",
                ));
            }
        }
        Ok(())
    }

    fn add_work(
        &mut self,
        duration_s: f64,
        substeps: usize,
        evaluations: usize,
        iterations: usize,
        discrepancy: f64,
    ) -> Result<(), OpticalHeatingError> {
        self.duration.add(duration_s)?;
        self.report.duration_s = self.duration.value()?;
        for (counter, increment) in [
            (&mut self.report.accepted_kinetic_substeps, substeps),
            (&mut self.report.cell_evaluations, evaluations),
            (&mut self.report.coupling_iterations, iterations),
        ] {
            let increment = u64::try_from(increment)
                .map_err(|_| invalid("depleting pump work counter overflow"))?;
            *counter = counter
                .checked_add(increment)
                .ok_or_else(|| invalid("depleting pump work counter overflow"))?;
        }
        self.report.maximum_energy_discrepancy_j =
            self.report.maximum_energy_discrepancy_j.max(discrepancy);
        if let Some(single) = &mut self.single {
            let beam = &self.report.beams[0];
            *single = DepletingPumpEnergy {
                duration_s: self.report.duration_s,
                incident_j: beam.incident_j,
                transmitted_j: beam.transmitted_j,
                absorbed_j: beam.absorbed_j,
                passive_absorbed_j: beam.passive_absorbed_j,
                residual_j: beam.residual_j,
                accepted_kinetic_substeps: self.report.accepted_kinetic_substeps,
                cell_evaluations: self.report.cell_evaluations,
                maximum_energy_discrepancy_j: self.report.maximum_energy_discrepancy_j,
            };
        }
        Ok(())
    }
}
