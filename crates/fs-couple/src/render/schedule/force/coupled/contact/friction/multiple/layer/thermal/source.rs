//! Trial-only source coefficient changes on the original simultaneous solver.
use super::*;
use interface::DryFrictionBinding;

pub(super) struct DryTemperatureUnknown {
    pub port: usize,
    original: DryFrictionBinding,
    scale_j_per_k: f64,
}

impl DryTemperatureUnknown {
    pub(super) fn prepare<E>(
        system: &MultiFrictionModalSystem, dry: &[DryContactHeat<'_>],
        thermal: &DepositionSimulation<'_>, dt: f64,
    ) -> Result<Vec<Self>, LayerThermalError<E>> {
        let mut unknowns = Vec::new();
        for route in dry {
            let original = system.interface(route.port).expect("admitted dry source");
            if original.temperature_source().is_some() {
                unknowns.push(Self {
                    port: route.port, original: original.clone(),
                    scale_j_per_k: route.temperature_scale(thermal, dt)?,
                });
            }
        }
        Ok(unknowns)
    }

    pub(super) fn temperature<E>(&self, coordinate_j: f64) -> Result<f64, LayerThermalError<E>> {
        number(self.original.temperature_kelvin() + coordinate_j / self.scale_j_per_k)
    }

    pub(super) fn binding<E>(&self, coordinate_j: f64)
        -> Result<(usize, DryFrictionBinding), LayerThermalError<E>>
    {
        Ok((self.port, self.original.at_temperature(self.temperature(coordinate_j)?)?))
    }
}

impl MultiFrictionModalSystem {
    // Keep the scalar port-zero owner and every multiport law in agreement.
    // Cached nonmonotone load lines do not identify their law in matches().
    pub(super) fn install_dry_source(&mut self, port: usize, source: DryFrictionBinding) {
        self.points[port].port.law = source.card().friction_law().clone();
        self.load_lines[port] = None;
        if port == 0 {
            self.core.port.law = source.card().friction_law().clone();
            self.core.interface = Some(source);
            self.core.load_line = None;
        } else {
            self.interfaces[port - 1] = Some(source);
        }
    }

    pub(super) fn stage_temperature_sources(
        &mut self, candidates: &[(usize, DryFrictionBinding)],
        external: &[f64], drives: &[FrictionDrive], obstacle_offsets: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<MultiFrictionFrame, ModalFrictionError> {
        if candidates.is_empty() {
            return self.stage_constitutive(external, drives, obstacle_offsets, Some(gate));
        }
        let original: Vec<_> = candidates.iter().map(|(port, _)| (*port,
            self.interface(*port).expect("admitted original source").clone(),
            self.points[*port].port.law.clone())).collect();
        let load_lines = self.load_lines.clone();
        let core_line = self.core.load_line.clone();
        let core_law = self.core.port.law.clone();
        for (port, candidate) in candidates {
            self.install_dry_source(*port, candidate.clone());
        }
        let staged = self.stage_constitutive(external, drives, obstacle_offsets, Some(gate));
        // Restore accepted sources, laws and caches BEFORE propagating any
        // mechanical/source/cancellation refusal. Only candidates are staged.
        for (port, source, law) in original {
            self.install_dry_source(port, source);
            self.points[port].port.law = law;
        }
        self.core.port.law = core_law;
        self.core.load_line = core_line;
        self.load_lines = load_lines;
        staged
    }
}
