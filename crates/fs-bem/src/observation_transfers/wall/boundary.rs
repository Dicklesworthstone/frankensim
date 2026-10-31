//! Local surface laws evaluated by their original mechanical/material owners.
use super::*;
use crate::helmholtz::bc::PorousRigidBoundary;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::porous::PorousFluidLaw;
use fs_phs::WallPin;

/// A passive locally reacting surface, in specific-impedance units [Pa s/m].
/// These are boundary models, not a new shell or moving porous-frame solver.
#[derive(Clone, Debug)]
pub enum LocalWall {
    Rigid,
    /// Original surface mass, stiffness and resistance, with exp(-i omega t).
    Mechanical(WallPin),
    /// Original JCA(L) source and one coherent saturating-fluid snapshot.
    /// The BEM adapter performs its existing time-convention conversion once.
    PorousRigid {
        law: PorousFluidLaw,
        fluid: ThermoviscousFluidState,
        thickness_m: f64,
    },
}

/// An immutable material table and one table index per physical boundary panel.
/// Each distinct law is evaluated once per frequency, irrespective of panel or
/// source count. Every declared law must be used; at most 256 are admitted.
#[derive(Clone, Debug)]
pub struct LocalBoundary {
    medium: Medium,
    laws: Vec<LocalWall>,
    panel_laws: Vec<usize>,
    supported_hz: (f64, f64),
    hash: ContentHash,
}

impl LocalBoundary {
    /// Bind physical panel assignments without defaults, interpolation or
    /// extrapolation. Porous exterior density and speed must match the actual
    /// saturating fluid; the material's positive-frequency support is retained.
    pub fn try_new(
        medium: Medium,
        laws: Vec<LocalWall>,
        panel_laws: Vec<usize>,
    ) -> Result<Self, TransferError> {
        if laws.is_empty()
            || laws.len() > 256
            || panel_laws.is_empty()
            || panel_laws.len() > helmholtz::MAX_DENSE_PANELS
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
            || !(medium.density * medium.sound_speed).is_finite()
        {
            return Err(TransferError::Invalid(
                "invalid local boundary dimensions or medium",
            ));
        }
        let mut used = vec![false; laws.len()];
        for &index in &panel_laws {
            let Some(flag) = used.get_mut(index) else {
                return Err(TransferError::Invalid(
                    "local boundary panel names an absent law",
                ));
            };
            *flag = true;
        }
        if used.iter().any(|used| !used) {
            return Err(TransferError::Invalid(
                "local boundary contains an unused law",
            ));
        }
        let mut supported_hz = (0.0_f64, f64::INFINITY);
        for law in &laws {
            match law {
                LocalWall::Rigid => {}
                LocalWall::Mechanical(pin) => {
                    // This original owner admits the physical coefficients.
                    fs_phs::wall_specific_impedance(pin, 1.0)?;
                }
                LocalWall::PorousRigid {
                    law,
                    fluid,
                    thickness_m,
                } => {
                    let band = law.reference_card().valid_hz;
                    let prepared = PorousRigidBoundary::new(
                        law,
                        fluid,
                        core::f64::consts::TAU * band.0,
                        *thickness_m,
                    )?;
                    if prepared.medium() != medium {
                        return Err(TransferError::Invalid(
                            "porous and exterior fluid coefficients differ",
                        ));
                    }
                    supported_hz.0 = supported_hz.0.max(band.0);
                    supported_hz.1 = supported_hz.1.min(band.1);
                }
            }
        }
        if supported_hz.0 >= supported_hz.1 {
            return Err(TransferError::Invalid(
                "local wall materials have no common frequency band",
            ));
        }
        let mut boundary = Self {
            medium,
            laws,
            panel_laws,
            supported_hz,
            hash: ContentHash([0; 32]),
        };
        boundary.hash = boundary.compute_hash();
        Ok(boundary)
    }

    pub const fn medium(&self) -> Medium {
        self.medium
    }
    pub fn panel_count(&self) -> usize {
        self.panel_laws.len()
    }
    /// Common source support [Hz]; infinity means no finite constitutive upper
    /// cutoff, not that mesh resolution or floating-point range is unlimited.
    pub const fn supported_band_hz(&self) -> (f64, f64) {
        self.supported_hz
    }
    pub const fn content_hash(&self) -> ContentHash {
        self.hash
    }

    pub(super) fn rigid_at(&self, panel: usize) -> bool {
        matches!(self.laws[self.panel_laws[panel]], LocalWall::Rigid)
    }

    pub(super) fn admit_band(
        &self,
        highest_omega: f64,
        lower: Option<LowerBandWindow>,
    ) -> Result<(), TransferError> {
        let tau = core::f64::consts::TAU;
        let lowest = lower.map_or(0.0, |window| window.stop_hz);
        if lowest < self.supported_hz.0
            || highest_omega > tau * self.supported_hz.1
            || lower.is_some_and(|window| tau * window.pass_hz >= highest_omega)
        {
            return Err(TransferError::Invalid(
                "observation band exceeds local-wall source support",
            ));
        }
        for law in &self.laws {
            if let LocalWall::Mechanical(pin) = law
                && pin.resistance == 0.0
            {
                let resonance = det::sqrt(pin.stiffness_per_area) / det::sqrt(pin.surface_density);
                if resonance > tau * lowest && resonance <= highest_omega {
                    return Err(TransferError::Invalid(
                        "lossless wall resonance requires a pressure-release boundary owner",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Dimensionless beta=rho*c/Z in exp(-i omega t), in original panel order.
    /// Native constitutive refusals propagate. A zero-impedance lossless
    /// resonance refuses; admittance is never clipped or replaced with rigid.
    pub fn normalized_admittance(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        if !omega.is_finite() || omega <= 0.0 {
            return Err(TransferError::Invalid(
                "local wall frequency must be positive and finite",
            ));
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let mut coefficients = Vec::with_capacity(self.laws.len());
        for law in &self.laws {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let beta = match law {
                LocalWall::Rigid => C64::ZERO,
                LocalWall::Mechanical(pin) => {
                    let impedance = fs_phs::wall_specific_impedance(pin, omega)?;
                    if !finite(impedance) || impedance == C64::ZERO {
                        return Err(TransferError::Invalid(
                            "wall impedance is singular or unrepresentable",
                        ));
                    }
                    C64::from_re(rho_c) / impedance
                }
                LocalWall::PorousRigid {
                    law,
                    fluid,
                    thickness_m,
                } => PorousRigidBoundary::new(law, fluid, omega, *thickness_m)?
                    .normalized_admittance(),
            };
            if !finite(beta) || beta.re < 0.0 {
                return Err(TransferError::Invalid(
                    "local wall admittance is nonfinite or active",
                ));
            }
            coefficients.push(beta);
        }
        Ok(self
            .panel_laws
            .iter()
            .map(|&index| coefficients[index])
            .collect())
    }

    fn compute_hash(&self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-bem:local-boundary:v1");
        hash.update(&(self.laws.len() as u64).to_le_bytes());
        hash.update(&(self.panel_laws.len() as u64).to_le_bytes());
        hash_scalars(&mut hash, &[self.medium.density, self.medium.sound_speed]);
        for law in &self.laws {
            match law {
                LocalWall::Rigid => hash.update(&[0]),
                LocalWall::Mechanical(pin) => {
                    hash.update(&[1]);
                    hash_scalars(
                        &mut hash,
                        &[pin.surface_density, pin.stiffness_per_area, pin.resistance],
                    );
                }
                LocalWall::PorousRigid {
                    law,
                    fluid,
                    thickness_m,
                } => {
                    hash.update(&[2]);
                    let card = law.reference_card();
                    for text in [card.name, card.provenance] {
                        hash.update(&(text.len() as u64).to_le_bytes());
                        hash.update(text.as_bytes());
                    }
                    hash_scalars(
                        &mut hash,
                        &[
                            card.porosity,
                            card.flow_resistivity,
                            card.tortuosity,
                            card.viscous_length_m,
                            card.thermal_length_m,
                            card.valid_hz.0,
                            card.valid_hz.1,
                            law.reference_dynamic_viscosity_pa_s(),
                            *thickness_m,
                            fluid.acoustic().density_kg_m3(),
                            fluid.acoustic().sound_speed_m_s(),
                            fluid.acoustic().adiabatic_bulk_modulus_pa(),
                            fluid.transport().density_kg_m3(),
                            fluid.transport().specific_heat_j_kg_k(),
                            fluid.transport().thermal_conductivity_w_m_k(),
                            fluid.transport().dynamic_viscosity_pa_s(),
                            fluid.heat_capacity_ratio(),
                            fluid.specific_heat_cv_j_kg_k(),
                        ],
                    );
                    match card.thermal_permeability_m2 {
                        None => hash.update(&[0]),
                        Some(value) => {
                            hash.update(&[1]);
                            hash_scalars(&mut hash, &[value]);
                        }
                    }
                }
            }
        }
        for &index in &self.panel_laws {
            hash.update(&(index as u64).to_le_bytes());
        }
        hash.finalize()
    }
}

fn hash_scalars(hash: &mut DomainHasher, values: &[f64]) {
    for value in values {
        hash.update(&value.to_bits().to_le_bytes());
    }
}
