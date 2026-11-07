//! Original solid resistance on a shared wall's gas-facing ports.

use super::*;
use crate::multicomponent::network::wall_exchange::{
    MaterialWallConduction, WallConductionEvaluation,
};

pub(super) struct SharedConduction {
    layers: Vec<Option<MaterialWallConduction>>,
    initial_conductance: Vec<f64>,
    // Only point-sampled state-dependent conductivity lacks an exact tangent.
    pub(super) implicit: bool,
    pub(super) integrated: bool,
}

impl Group<'_> {
    pub(crate) fn bind_conduction(
        &mut self,
        layers: Vec<Option<MaterialWallConduction>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        if self.conduction.is_some() || layers.len() != self.accepted.ports.len() {
            return Err(WallExchangeError::Invalid(
                "shared wall conduction port count or state",
            ));
        }
        let caloric = self.accepted.storage.heat_storage.as_ref();
        let mut implicit = false;
        let mut integrated = false;
        for layer in layers.iter().flatten() {
            poll(cancelled)?;
            integrated |= layer.has_integrated_profile();
            if layer.evaluation() == WallConductionEvaluation::ImplicitEndpoint {
                if caloric.is_none() {
                    return Err(WallExchangeError::Invalid(
                        "implicit shared wall conductivity requires original caloric storage",
                    ));
                }
                implicit |= !layer.has_integrated_profile();
            }
            // This retains source identity and the phase chart's original h
            // coordinate. The storage's own film remains exactly zero.
            layer.conductance(self.accepted.storage.thermal, caloric)?;
        }
        if layers.iter().any(Option::is_some) {
            let initial_conductance = super::super::super::allocated(layers.len(), 0.0)?;
            self.conduction = Some(SharedConduction {
                layers,
                initial_conductance,
                implicit,
                integrated,
            });
        }
        poll(cancelled)?;
        Ok(())
    }

    pub(crate) fn conduction_layers(&self) -> Option<&[Option<MaterialWallConduction>]> {
        self.conduction.as_ref().map(|c| c.layers.as_slice())
    }

    /// A later wall link may move the retained wall beyond the supported FACE
    /// of an earlier gas-port solve. Admit the complete final profile against the
    /// unchanged original gas energies, without applying another heat transfer.
    pub(crate) fn admit_integrated_profile_endpoints<G: WallGasCaloric>(
        &self,
        gas: &impl Fn(SharedWallGas, f64) -> G,
        thermal: WallThermal,
        duration: Time,
        control: WallControl,
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        let Some(binding) = &self.conduction else {
            return Ok(());
        };
        if !binding.integrated {
            return Ok(());
        }
        let wall = self.after_state(thermal, 0.0, work)?;
        for (i, port) in self.accepted.ports.iter().enumerate() {
            let Some(layer) = binding.profile(i) else {
                continue;
            };
            poll(cancelled)?;
            let original = MeteredGas {
                original: gas(port.gas, wall.temperature),
                work,
            };
            let state = original.primitive(self.initial[i][0], self.energies[i], cancelled)?;
            let bound = layer
                .profile_gas_sensitivity_bound(port.conductance_w_per_k)
                .expect("selected integrated profile");
            let temperature = refine_caloric_temperature_with_allowance(
                &original,
                state.temperature,
                self.energies[i] - 0.5 * self.initial[i][0] * state.velocity,
                self.volumes[i],
                bound * duration.value(),
                control.tolerance(self.heat[i].abs()),
                control,
                cancelled,
            )?;
            work.charge()?;
            layer.profile_response(
                [wall.temperature, temperature],
                port.conductance_w_per_k,
                duration,
                control,
                cancelled,
            )?;
        }
        Ok(())
    }

    /// Lagged coefficients still require a supported final material state.
    /// Later wall-link stages request all layers, because they can move a wall
    /// beyond the constitutive endpoint of its preceding gas-contact solve.
    pub(crate) fn admit_conduction_endpoint(
        &self,
        thermal: WallThermal,
        include_implicit: bool,
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        if let Some(binding) = &self.conduction {
            for layer in binding.layers.iter().flatten() {
                poll(cancelled)?;
                if include_implicit || layer.evaluation() == WallConductionEvaluation::InitialState
                {
                    work.charge()?;
                    layer.conductance(thermal, self.accepted.storage.heat_storage.as_ref())?;
                }
            }
        }
        Ok(())
    }
}

impl SharedConduction {
    pub(super) fn profile(&self, index: usize) -> Option<&MaterialWallConduction> {
        self.layers[index]
            .as_ref()
            .filter(|layer| layer.has_integrated_profile())
    }

    pub(super) fn prepare(
        &mut self,
        thermal: WallThermal,
        caloric: Option<&MaterialWallCaloric<'_>>,
        ports: &[SharedWallPort],
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        for (i, layer) in self.layers.iter().enumerate() {
            poll(cancelled)?;
            if let Some(layer) = layer {
                if layer.evaluation() == WallConductionEvaluation::InitialState {
                    work.charge()?;
                    let solid = layer.conductance(thermal, caloric)?.0;
                    self.initial_conductance[i] = series(ports[i].conductance_w_per_k, solid)?;
                }
            }
        }
        Ok(())
    }

    pub(super) fn conductance(
        &self,
        index: usize,
        film: f64,
        thermal: WallThermal,
        caloric: Option<&MaterialWallCaloric<'_>>,
        work: &Work,
    ) -> Result<f64, WallExchangeError> {
        match &self.layers[index] {
            None => Ok(film),
            Some(layer) if layer.evaluation() == WallConductionEvaluation::InitialState => {
                Ok(self.initial_conductance[index])
            }
            Some(layer) => {
                work.charge()?;
                let solid = layer.conductance(thermal, caloric)?.0;
                series(film, solid)
            }
        }
    }
}

// Stable resistance addition, preserving an explicitly disabled film without
// reciprocal overflow. The authored film is immutable; the combined value is
// never fed back as a second film in a later source half.
pub(super) fn series(film: f64, solid: f64) -> Result<f64, WallExchangeError> {
    let smaller = film.min(solid);
    let combined = smaller / (1.0 + smaller / film.max(solid));
    if !combined.is_finite() || combined < 0.0 || film > 0.0 && combined == 0.0 {
        return Err(WallExchangeError::Invalid(
            "shared wall series conductance is unrepresentable",
        ));
    }
    Ok(combined)
}
