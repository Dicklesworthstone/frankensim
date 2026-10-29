//! Molecular compression memory on the existing reciprocal cavity network.
//! The material owner supplies the exact positive bulk spectrum. Geometry maps
//! pressure/strain into force/displacement; the original coupled propagator
//! advances all surfaces, neck, radiation and molecular histories together.
//! Gas temperature/composition remain frozen, and the cavity wall is adiabatic.
//! A wall-thermal shunt cannot be added independently to this caloric law.

use super::*;
use fs_material::relaxation::{RelaxationState, bulk::BulkRelaxation};

/// Source reservoirs and the exact volume/area chart used by the cavity.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularCavityCompression {
    bulk: BulkRelaxation,
    stiffness_per_bulk_m: f64,
}

impl MolecularCavityCompression {
    fn new(
        source: &HelmholtzModalCavity,
        relaxation: RelaxationState,
    ) -> Result<Self, ModalCouplingError> {
        if relaxation.gas() != source.gas() {
            return Err(network_invalid(
                "molecular cavity law must retain the exact original gas snapshot",
            ));
        }
        let area = source.neck_area_m2();
        let scale = area * (area / source.geometry().volume_m3);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(network_invalid(
                "molecular cavity compression chart is unrepresentable",
            ));
        }
        let bulk = relaxation.bulk_relaxation().map_err(|_| {
            network_invalid("original molecular bulk realization refused its source spectrum")
        })?;
        let law = Self {
            bulk,
            stiffness_per_bulk_m: scale,
        };
        if !law.equilibrium_stiffness_n_m().is_finite()
            || law.equilibrium_stiffness_n_m() <= 0.0
            || !law.frozen_stiffness_n_m().is_finite()
            || law.frozen_stiffness_n_m() <= 0.0
        {
            return Err(network_invalid(
                "molecular cavity compression storage is unrepresentable",
            ));
        }
        Ok(law)
    }

    /// Source-owned bulk modulus with species, original rates and heat capacities.
    pub const fn bulk(&self) -> &BulkRelaxation {
        &self.bulk
    }
    /// Locally equilibrated adiabatic compression stiffness [N/m].
    pub fn equilibrium_stiffness_n_m(&self) -> f64 {
        self.bulk.equilibrium_bulk_modulus_pa() * self.stiffness_per_bulk_m
    }
    /// Instantaneous compression with every internal reservoir frozen [N/m].
    pub fn frozen_stiffness_n_m(&self) -> f64 {
        self.bulk.frozen_bulk_modulus_pa() * self.stiffness_per_bulk_m
    }

    pub(super) fn branches(
        &self,
        column: &[f64],
    ) -> Result<Vec<fs_phs::RelaxationBranch>, ModalCouplingError> {
        let mut projection = vec![0.0; 2 * column.len()];
        for (i, &value) in column.iter().enumerate() {
            projection[2 * i] = value;
        }
        self.bulk
            .arms()
            .iter()
            .map(|arm| {
                let stiffness = arm.delta_bulk_modulus_pa * self.stiffness_per_bulk_m;
                if !stiffness.is_finite() || stiffness <= 0.0 {
                    return Err(network_invalid(
                        "molecular cavity memory stiffness is unrepresentable",
                    ));
                }
                Ok(fs_phs::RelaxationBranch {
                    projection: projection.clone(),
                    stiffness,
                    relaxation_time_s: arm.relaxation_time_s,
                })
            })
            .collect()
    }

    pub(super) fn observation(
        &self,
        extension: f64,
        memory: &[f64],
    ) -> Result<(f64, f64), ModalCouplingError> {
        if memory.len() != self.bulk.arms().len() {
            return Err(network_invalid(
                "molecular cavity observation requires every accepted reservoir history",
            ));
        }
        let equilibrium = self.equilibrium_stiffness_n_m();
        let mut force = equilibrium * extension;
        let mut energy = 0.5 * equilibrium * extension * extension;
        for (arm, &z) in self.bulk.arms().iter().zip(memory) {
            let root = det::sqrt(arm.delta_bulk_modulus_pa * self.stiffness_per_bulk_m);
            let slip = root * extension - z;
            force += root * slip;
            energy += 0.5 * slip * slip;
        }
        Ok((force, energy))
    }

    pub(super) fn harmonic_observation(
        &self,
        omega: f64,
        extension: C64,
    ) -> Result<(C64, f64), ModalCouplingError> {
        let stiffness = self
            .bulk
            .response(omega)
            .map_err(|_| {
                network_invalid("molecular cavity harmonic bulk response is unrepresentable")
            })?
            .scale(self.stiffness_per_bulk_m);
        let mut energy = 0.25 * self.equilibrium_stiffness_n_m() * extension.norm_sq();
        let s = C64::new(0.0, -omega);
        for arm in self.bulk.arms() {
            let slip = s / (s + C64::from_re(arm.relaxation_time_s.recip())) * extension;
            energy += 0.25 * arm.delta_bulk_modulus_pa * self.stiffness_per_bulk_m * slip.norm_sq();
        }
        Ok((stiffness * extension, energy))
    }
}

impl SharedCavity {
    /// Use the original gas's molecular reservoirs in the SAME time solve.
    ///
    /// Geometry and equilibrium gas must match the supplied material law exactly.
    /// Changing species, collision partners, temperature or pressure reevaluates
    /// the upstream material source; no gas-specific branch or modal-Q input is
    /// needed here. Every moving surface sees the same reciprocal pressure.
    /// Optional viscous-neck memory and subsequent exterior/DG attachment retain
    /// their original owners, histories, budgets and work accounts.
    ///
    /// `None` initializes molecular history relaxed at the supplied displacement.
    /// `Some(z)` sets explicit energy-normalized histories in bulk().arms() order;
    /// zeros represent an instantaneous displacement with frozen reservoirs.
    /// These internal states use the existing material-memory dissipation account.
    /// They do not deposit the same acoustic loss again as gas heat.
    ///
    /// The cavity remains spatially uniform and wall-adiabatic. For a conducting
    /// cylindrical wall, use [`Self::new_with_cylindrical_molecular_wall`], whose
    /// joint caloric law replaces this adiabatic specialization. The two owners
    /// cannot be added as independent compression shunts. This does not evolve
    /// finite gas temperature/species or apply propagation loss to uniform pressure.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_molecular_relaxation(
        source: HelmholtzModalCavity,
        surfaces: Vec<CavitySurface>,
        sample_rate_hz: u32,
        neck_initial: ModalAcousticState,
        neck_budget: ModalAcousticTimeBudget,
        config: ModalCouplingConfig,
        neck_options: Option<ViscousNeckOptions>,
        relaxation: RelaxationState,
        initial_molecular_memory: Option<Vec<f64>>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled);
        }
        let molecular = MolecularCavityCompression::new(&source, relaxation)?;
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled);
        }
        let viscous = neck_options
            .map(|options| {
                Self::prepare_viscous_neck(&source, sample_rate_hz, config, options, gate)
            })
            .transpose()?;
        Self::new_inner(
            source,
            surfaces,
            sample_rate_hz,
            neck_initial,
            neck_budget,
            config,
            viscous,
            None,
            Some((molecular, initial_molecular_memory)),
            gate,
        )
    }
}

#[cfg(test)]
mod tests;
