//! Circular thermoviscous diffusion in the existing DG exponential owner.
//!
//! Uniform circular geometry and gas make each normalized wall coefficient
//! constant in space. It therefore commutes with the exact DG mass maps; no
//! nodal mass lumping, collocated gas averaging, or second spatial solver is
//! needed. With normalized main coordinate x and memory y, each arm adds
//! `[-g, sqrt(g*p); sqrt(g*p), -p]`, the negative outer product of
//! `[sqrt(g),-sqrt(p)]`. Storage is `(x²+y²)/2`; dissipation is
//! `(sqrt(g)*x-sqrt(p)*y)²`. Every state advances in ONE augmented exponential.

use super::*;

mod spatial;
pub use spatial::{DgWallAdmission, DgWallResolution};
mod heat;
pub use heat::DgViscousHeat;

#[derive(Debug, Clone)]
pub(super) enum ViscothermalMemory {
    Uniform(fs_phs::ZkDiffusive),
    Spatial(spatial::SpatialViscothermal),
}

pub(super) fn uniform_radius(duct: &Duct) -> Option<f64> {
    let Some(Segment::Cylinder { radius, .. }) = duct.segments.first() else {
        return None;
    };
    duct.segments
        .iter()
        .all(|segment| matches!(segment, Segment::Cylinder { radius: r, .. } if r == radius))
        .then_some(*radius)
}

impl DgDuct {
    /// Attach exact-Bessel-derived viscous and thermal memory before creating a
    /// state or propagator. Gas density, viscosity, heat-capacity ratio and
    /// Prandtl number all come from the already bound material state.
    ///
    /// `poles` is a refinable diffusion discretization for EACH wall operator.
    /// The same count covers any number of adjacent cylinders of identical
    /// radius. This first attachment requires uniform gas and radius; a taper,
    /// profile or nonuniform gas field refuses instead of averaging its data.
    /// Existing lossless variable-area and gas-field images remain available.
    ///
    /// The pole producer admits both complex wall operators and real losses on
    /// a 129-point logarithmic grid over the declared physical frequency band.
    /// The DG spatial and time-propagator errors are separate axes. The wall
    /// memories use the EXISTING certified DenseSeries exponential, with no
    /// bilinear relaxation, changed clock or extra inviscid storage.
    ///
    /// The existing dense ledger still infers dissipation from energy deficit
    /// and reports `balance_checked=false`; it is not an independent resistor
    /// integral. This attachment does not close the full M2.2-W3b ledger gate.
    ///
    /// # Errors
    /// Unsupported geometry/gas variation, repeated attachment, unresolved pole
    /// band, or total acoustic-plus-memory DOFs above the existing dense cap.
    pub fn with_viscothermal(
        self,
        poles: usize,
        band_hz: [f64; 2],
        maximum_relative_error: f64,
    ) -> Result<Self, DuctError> {
        if self.compliant_wall.is_some() || self.molecular.is_some() {
            return Err(DuctError::BadParameter {
                what: "attach DG viscothermal memory before compliant wall mechanics",
            });
        }
        if self.viscothermal.is_some() {
            return Err(DuctError::BadParameter {
                what: "DG viscothermal memory is already attached",
            });
        }
        let Some(radius) = self.uniform_radius_m else {
            return Err(DuctError::BadParameter {
                what: "DG viscothermal memory currently requires identical circular cylinders",
            });
        };
        if self
            .gas_field
            .as_ref()
            .is_some_and(|field| field.uniform_state().is_none())
        {
            return Err(DuctError::BadParameter {
                what: "DG viscothermal memory currently requires uniform material gas",
            });
        }
        let gas = &self.gas;
        let model = fs_phs::ZkDiffusive::new(
            radius,
            gas.density,
            gas.dynamic_viscosity,
            gas.gamma,
            gas.prandtl,
            poles,
            band_hz.map(|f| core::f64::consts::TAU * f),
            maximum_relative_error,
        )
        .map_err(|_| DuctError::BadParameter {
            what: "DG viscothermal Bessel pole data or declared band is unresolved",
        })?;
        self.with_uniform_viscothermal_model(model)
    }
}

impl<M: DgMedium> DgDuct<M> {
    // Stamp the SAME native uniform Bessel realization. Live material refresh
    // uses this owner with the original angular band, radius and pole count.
    pub(super) fn with_uniform_viscothermal_model(
        mut self,
        model: fs_phs::ZkDiffusive,
    ) -> Result<Self, DuctError> {
        if self.viscothermal.is_some()
            || self.compliant_wall.is_some()
            || self.molecular.is_some()
            || self.gas.porous_diffusive_model().is_some()
        {
            return Err(DuctError::BadParameter {
                what: "uniform viscothermal stamping requires the bare acoustic owner",
            });
        }
        self.stamp_uniform_diffusion(
            model.series_drag_per_s(),
            model
                .series_terms()
                .iter()
                .map(|t| (t.gain_per_s, t.rate_per_s)),
            model
                .shunt_terms()
                .iter()
                .map(|t| (t.gain_per_s, t.rate_per_s)),
        )?;
        self.viscothermal = Some(ViscothermalMemory::Uniform(model));
        Ok(self)
    }

    // One energy-normalized owner for circular wall and intrinsic porous
    // positive arms. Their original producers retain material/source admission.
    // The primary acoustic masses already own the infinite-frequency storage.
    pub(super) fn stamp_uniform_diffusion(
        &mut self,
        series_drag_per_s: f64,
        series: impl ExactSizeIterator<Item = (f64, f64)>,
        shunt: impl ExactSizeIterator<Item = (f64, f64)>,
    ) -> Result<(), DuctError> {
        let m = self.field_dofs;
        let original = self.total_dofs();
        if original != 2 * m
            || self.memory_dofs != 0
            || !series_drag_per_s.is_finite()
            || series_drag_per_s < 0.0
        {
            return Err(DuctError::BadParameter {
                what: "uniform diffusion stamping needs bare acoustic storage and finite nonnegative drag",
            });
        }
        let branches = series
            .len()
            .checked_add(shunt.len())
            .ok_or(DuctError::BadParameter {
                what: "DG uniform diffusion branch count overflow",
            })?;
        let memory_dofs = m.checked_mul(branches).ok_or(DuctError::BadParameter {
            what: "DG uniform diffusion memory size overflow",
        })?;
        let total = original
            .checked_add(memory_dofs)
            .ok_or(DuctError::BadParameter {
                what: "DG uniform diffusion total size overflow",
            })?;
        if total > MAX_DG_PROPAGATOR_DOFS {
            return Err(DuctError::BadParameter {
                what: "DG acoustic plus diffusion memory exceeds the dense DOF cap",
            });
        }
        let mut operator = vec![0.0; total * total];
        for i in 0..original {
            operator[i * total..i * total + original]
                .copy_from_slice(&self.l_tilde[i * original..(i + 1) * original]);
        }
        for i in 0..m {
            operator[(m + i) * total + m + i] -= series_drag_per_s;
        }
        for (branch, (field_offset, gain, rate)) in series
            .map(|(gain, rate)| (m, gain, rate))
            .chain(shunt.map(|(gain, rate)| (0, gain, rate)))
            .enumerate()
        {
            if !gain.is_finite() || gain <= 0.0 || !rate.is_finite() || rate <= 0.0 {
                return Err(DuctError::BadParameter {
                    what: "DG uniform diffusion needs positive finite original poles and weights",
                });
            }
            let coupling = det::sqrt(gain) * det::sqrt(rate);
            for i in 0..m {
                let main = field_offset + i;
                let memory = original + branch * m + i;
                operator[main * total + main] -= gain;
                operator[main * total + memory] = coupling;
                operator[memory * total + main] = coupling;
                operator[memory * total + memory] = -rate;
            }
        }
        if operator.iter().any(|value| !value.is_finite()) {
            return Err(DuctError::BadParameter {
                what: "DG uniform diffusion augmented operator is not finite",
            });
        }
        for (_, port) in &mut self.ports {
            port.resize(total, 0.0);
        }
        self.l_tilde = operator;
        self.g_block = None;
        self.memory_dofs = memory_dofs;
        self.operator_identity = hash_operator(
            self.mesh,
            self.inlet,
            self.outlet,
            &self.elements,
            &self.l_tilde,
            &self.ports,
        );
        Ok(())
    }

    /// Admitted angular-frequency band of circular wall diffusion or intrinsic
    /// porous equivalent-fluid memory. None means no declared diffusion band;
    /// other physical losses can still be present.
    #[must_use]
    pub fn viscothermal_band_rad_s(&self) -> Option<[f64; 2]> {
        match &self.viscothermal {
            Some(ViscothermalMemory::Uniform(model)) => Some(model.band_rad_s()),
            Some(ViscothermalMemory::Spatial(wall)) => Some(wall.admission.resolution.band_rad_s),
            None => self
                .gas
                .porous_diffusive_model()
                .map(|model| model.band_rad_s()),
        }
    }

    /// Uniform circular-wall pole model and its admitted band/error. A spatial
    /// wall reports spatial_viscothermal_admission. None does not imply a
    /// lossless duct: intrinsic porous memory has its own porous_diffusive_model.
    #[must_use]
    pub fn viscothermal_model(&self) -> Option<&fs_phs::ZkDiffusive> {
        match &self.viscothermal {
            Some(ViscothermalMemory::Uniform(model)) => Some(model),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests;
