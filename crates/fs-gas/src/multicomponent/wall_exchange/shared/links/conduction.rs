//! Two original bulk resistances around a finite authored wall interface.

use super::*;
use crate::multicomponent::network::wall_exchange::{
    MaterialWallConduction, WallConductionEvaluation,
};
use crate::multicomponent::wall_exchange::shared::conduction::series;

pub(crate) struct LinkConduction {
    layers: [Option<MaterialWallConduction>; 2],
}

impl WallLinks {
    pub(crate) fn bind_conduction<'a, 'law: 'a>(
        &mut self,
        layers: Vec<[Option<MaterialWallConduction>; 2]>,
        group: impl Fn(usize) -> &'a Group<'law>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        if self.conduction.is_some()
            || layers.len() != self.ports.len()
            || !layers.iter().flatten().any(Option::is_some)
        {
            return Err(WallExchangeError::Invalid(
                "wall link material binding count or state",
            ));
        }
        let mut bindings = Vec::new();
        bindings
            .try_reserve_exact(layers.len())
            .map_err(|_| WallExchangeError::Invalid("wall link material allocation"))?;
        for (link, layers) in self.ports.iter().zip(layers) {
            poll(cancelled)?;
            for (side, layer) in layers.iter().enumerate() {
                if let Some(layer) = layer {
                    let wall = &group(link.walls[side]).accepted.storage;
                    if layer.evaluation() == WallConductionEvaluation::ImplicitEndpoint
                        && wall.heat_storage.is_none()
                    {
                        return Err(WallExchangeError::Invalid(
                            "implicit link conductivity requires original caloric storage",
                        ));
                    }
                    layer.conductance(wall.thermal, wall.heat_storage.as_ref())?;
                }
            }
            bindings.push(if layers.iter().any(Option::is_some) {
                Some(LinkConduction { layers })
            } else {
                None
            });
        }
        poll(cancelled)?;
        self.conduction = Some(bindings);
        Ok(())
    }

    pub(crate) fn conduction_layers(
        &self,
        index: usize,
    ) -> Option<&[Option<MaterialWallConduction>; 2]> {
        Some(&self.conduction.as_ref()?.get(index)?.as_ref()?.layers)
    }

    /// Later stages can change an earlier link's material state. Admit the
    /// final candidate against EVERY retained source before full publication.
    pub(crate) fn admit_conduction_endpoints<'a, 'law: 'a>(
        &self,
        group: impl Fn(usize) -> &'a Group<'law>,
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        if let Some(bindings) = &self.conduction {
            for (link, binding) in self.ports.iter().zip(bindings) {
                if let Some(binding) = binding {
                    let groups = link.walls.map(&group);
                    binding.admit_endpoint(
                        groups,
                        groups.map(|g| g.candidate),
                        false,
                        work,
                        cancelled,
                    )?;
                }
            }
        }
        Ok(())
    }
}

impl LinkConduction {
    pub(super) fn prepare<'a>(
        &'a self,
        groups: [&Group<'_>; 2],
        states: [WallThermal; 2],
        interface: f64,
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Prepared<'a>, WallExchangeError> {
        let mut result = Prepared {
            binding: self,
            solid: [None; 2],
            initial: interface,
            implicit: false,
        };
        for (i, layer) in self.layers.iter().enumerate() {
            if let Some(layer) = layer {
                poll(cancelled)?;
                work.charge()?;
                let solid = layer
                    .conductance(states[i], groups[i].accepted.storage.heat_storage.as_ref())?
                    .0;
                result.solid[i] = Some(solid);
                result.initial = series(result.initial, solid)?;
                result.implicit |= layer.evaluation() == WallConductionEvaluation::ImplicitEndpoint;
            }
        }
        Ok(result)
    }

    fn admit_endpoint(
        &self,
        groups: [&Group<'_>; 2],
        states: [WallThermal; 2],
        only_lagged: bool,
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        for (i, layer) in self.layers.iter().enumerate() {
            if let Some(layer) = layer {
                if !only_lagged || layer.evaluation() == WallConductionEvaluation::InitialState {
                    poll(cancelled)?;
                    work.charge()?;
                    layer
                        .conductance(states[i], groups[i].accepted.storage.heat_storage.as_ref())?;
                }
            }
        }
        Ok(())
    }
}

pub(super) struct Prepared<'a> {
    binding: &'a LinkConduction,
    solid: [Option<f64>; 2],
    pub(super) initial: f64,
    pub(super) implicit: bool,
}

impl Prepared<'_> {
    pub(super) fn at(
        &self,
        groups: [&Group<'_>; 2],
        states: [WallThermal; 2],
        mut interface: f64,
        work: &Work,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, WallExchangeError> {
        // With only lagged layers, retain the prepared coefficient verbatim.
        if !self.implicit {
            return Ok(self.initial);
        }
        for (i, layer) in self.binding.layers.iter().enumerate() {
            if let Some(layer) = layer {
                let solid = if layer.evaluation() == WallConductionEvaluation::ImplicitEndpoint {
                    if cancelled() {
                        return Err(WallExchangeError::Cancelled);
                    }
                    work.charge()?;
                    layer
                        .conductance(states[i], groups[i].accepted.storage.heat_storage.as_ref())?
                        .0
                } else {
                    self.solid[i].expect("prepared lagged material")
                };
                interface = series(interface, solid)?;
            }
        }
        Ok(interface)
    }

    pub(super) fn admit_lagged_endpoint(
        &self,
        groups: [&Group<'_>; 2],
        states: [WallThermal; 2],
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        self.binding
            .admit_endpoint(groups, states, true, work, cancelled)
    }
}
