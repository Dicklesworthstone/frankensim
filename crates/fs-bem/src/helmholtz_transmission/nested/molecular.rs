//! Original molecular compression in all three regions of the nested solver.
use super::*;
use crate::ComplexSum;
use crate::helmholtz::laplace::DampedWaveNumber;
use crate::helmholtz::transmission::molecular::{admit_products, kernels as molecular_kernels};
use fs_material::relaxation::RelaxationState;

#[path = "molecular/fields.rs"]
mod receivers;

/// One stationary shell and core with source-owned molecular phase and loss.
/// The source order is exterior, shell, core. All mean pressures must agree;
/// momentum remains inviscid with real rho*omega. No extra bulk attenuation,
/// shear viscosity, mean flow, mixing, or heat transport is introduced.
pub struct MolecularCoreShellFluidInterface<'a> {
    surfaces: [&'a SpherePanels; 2],
    sources: [&'a RelaxationState; 3],
    omega: f64,
    waves: [DampedWaveNumber; 3],
    impedance: f64,
    ppw: [f64; 2],
    offsets: [Option<usize>; 2],
    dimension: usize,
    drive: Vec<C64>,
    factor: Option<LuComplex>,
    condition: f64,
}

impl<'a> MolecularCoreShellFluidInterface<'a> {
    /// Reuse native geometric, wavelength, dense-memory and worker admission.
    /// Evaluate each distinct original molecular owner once at this frequency.
    /// Exact equality of density AND complex k removes an interface; matching
    /// equilibrium sound speed alone does not remove molecular refraction.
    pub fn new(
        outer: &'a SpherePanels,
        core: &'a SpherePanels,
        sources: [&'a RelaxationState; 3],
        omega_rad_s: f64,
        minimum_interface_clearance_m: f64,
        options: Options,
        gate: &CancelGate,
    ) -> Result<Self, TransmissionError> {
        poll(gate)?;
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(bad("molecular core-shell frequency must be positive and finite").into());
        }
        let workers = geometry::admit(outer, core, minimum_interface_clearance_m, options, gate)?;
        for source in sources {
            source.gas().validate_coherence()?;
            if source.gas().pressure != sources[0].gas().pressure {
                return Err(
                    bad("stationary molecular core-shell regions need equal mean pressures").into(),
                );
            }
        }
        let mut waves = [DampedWaveNumber {
            k: C64::ZERO,
            phase_origin_m: 0.0,
        }; 3];
        for region in 0..3 {
            poll(gate)?;
            waves[region] = if let Some(previous) =
                (0..region).find(|&index| core::ptr::eq(sources[index], sources[region]))
            {
                waves[previous]
            } else {
                DampedWaveNumber {
                    k: sources[region].relaxing_wavenumber(omega_rad_s)?,
                    phase_origin_m: 0.0,
                }
            };
            let k = waves[region].k;
            if !finite(k) || k.re <= 0.0 || k.im < 0.0 {
                return Err(
                    bad("molecular core-shell waves must be finite, outgoing and passive").into(),
                );
            }
        }
        let density = sources.map(|source| source.gas().density);
        let k_reference = omega_rad_s / sources[0].gas().sound_speed;
        let beta = [k_reference, k_reference * (density[1] / density[0])];
        let impedance = density[0] * sources[0].gas().sound_speed;
        if beta
            .into_iter()
            .chain([impedance])
            .chain(density.map(|rho| omega_rad_s * rho))
            .any(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(bad("unrepresentable molecular core-shell momentum scales").into());
        }
        let surfaces = [outer, core];
        let mut ppw = [
            admit_resolution(outer, waves[0].k.abs())?
                .min(admit_resolution(outer, waves[1].k.abs())?),
            admit_resolution(core, waves[1].k.abs())?
                .min(admit_resolution(core, waves[2].k.abs())?),
        ];
        let mut offsets = [None; 2];
        let mut dimension = 0;
        for boundary in 0..2 {
            if density[boundary] != density[boundary + 1]
                || waves[boundary].k != waves[boundary + 1].k
            {
                offsets[boundary] = Some(dimension);
                dimension += 2 * surfaces[boundary].areas().len();
            }
        }
        if offsets.iter().all(Option::is_some) {
            ppw[1] = ppw[1].min(admit_resolution(core, waves[0].k.abs())?);
        }
        let (drive, factor, condition) = if dimension == 0 {
            (Vec::new(), None, 1.0)
        } else {
            let (operator, drive) = assembly::NestedAssembly {
                surfaces,
                density,
                waves,
                beta,
                offsets,
                dimension,
            }
            .assemble(gate, molecular_kernels::difference, molecular_kernels::self_contrast)?;
            let (operator, factor, condition) =
                factor_contrast(operator, dimension, workers, gate)?;
            drop(operator);
            (drive, Some(factor), condition)
        };
        poll(gate)?;
        Ok(Self {
            surfaces,
            sources,
            omega: omega_rad_s,
            waves,
            impedance,
            ppw,
            offsets,
            dimension,
            drive,
            factor,
            condition,
        })
    }

    pub fn surface(&self, interface: CoreShellInterface) -> &'a SpherePanels {
        self.surfaces[interface.index()]
    }
    pub fn material_sources(&self) -> [&'a RelaxationState; 3] {
        self.sources
    }
    pub fn wavenumbers_per_m(&self) -> [C64; 3] {
        self.waves.map(|wave| wave.k)
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn panels_per_wavelength(&self, interface: CoreShellInterface) -> f64 {
        self.ppw[interface.index()]
    }
    pub fn active_interface_count(&self) -> usize {
        self.offsets.iter().filter(|offset| offset.is_some()).count()
    }
    pub fn condition_lower_bound(&self) -> f64 {
        self.condition
    }
    fn radiating_interface(&self) -> usize {
        if self.offsets[0].is_some() || self.offsets[1].is_none() {
            0
        } else {
            1
        }
    }

    /// Solve one original coherent source field, regular throughout the outer
    /// volume, with all shell/core reflections. The original law/factorization
    /// is reused, without material resampling or a post-hoc attenuation factor.
    /// Source cardinality is bounded across BOTH surfaces and each later query.
    pub fn scatter_field<'p>(
        &'p self,
        incident: IncidentField<'p>,
        minimum_source_clearance_m: f64,
        maximum_source_point_products: u64,
        gate: &CancelGate,
    ) -> Result<MolecularCoreShellScatteredField<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        incident.validate()?;
        admit_products(
            incident.source_count(),
            self.surfaces.iter().map(|s| s.areas().len()).sum(),
            maximum_source_point_products,
        )?;
        if !minimum_source_clearance_m.is_finite() || minimum_source_clearance_m <= 0.0 {
            return Err(bad("molecular core-shell source clearance must be positive and finite"));
        }
        let positions = incident.source_positions();
        if !positions.is_empty() {
            near_field::Geometry::new_iterative(
                self.surfaces[0], &positions, minimum_source_clearance_m, gate,
            )?;
        }
        let mut reference = [(Vec::new(), Vec::new()), (Vec::new(), Vec::new())];
        let mut incoming = zero_matrix(self.dimension, "molecular core-shell incident traces")?;
        for (boundary, (pressure, velocity)) in reference.iter_mut().enumerate() {
            let surface = self.surfaces[boundary];
            let n = surface.areas().len();
            for (&point, &normal) in surface.centroids().iter().zip(surface.normals()) {
                poll(gate)?;
                let sample = incident.sample_wave(
                    self.waves[0],
                    self.sources[0].gas().density * self.omega,
                    None,
                    point,
                    &mut || gate.is_requested(),
                )?;
                pressure.push(sample.pressure_pa);
                let mut normal_velocity = ComplexSum::default();
                for (component, direction) in sample.velocity_m_s.into_iter().zip(normal) {
                    normal_velocity.add(component.scale(direction));
                }
                velocity.push(normal_velocity.finish());
            }
            if let Some(offset) = self.offsets[boundary] {
                incoming[offset..offset + n].copy_from_slice(pressure);
                for (j, &value) in velocity.iter().enumerate() {
                    incoming[offset + n + j] = value * C64::new(0.0, self.impedance);
                }
            }
        }
        if incoming.iter().any(|&value| !finite(value)) {
            return Err(bad("nonfinite molecular core-shell incident traces"));
        }
        let correction =
            assembly::solve_correction(&self.drive, self.factor.as_ref(), &incoming, gate)?;
        let mut traces = [None, None];
        for boundary in 0..2 {
            let Some(offset) = self.offsets[boundary] else {
                continue;
            };
            let n = self.surfaces[boundary].areas().len();
            let pressure: Vec<_> = (0..n)
                .map(|i| reference[boundary].0[i] + correction[offset + i])
                .collect();
            let velocity: Vec<_> = (0..n)
                .map(|i| {
                    reference[boundary].1[i]
                        + physical_velocity(correction[offset + n + i], self.impedance)
                })
                .collect();
            let (_, work_interval, work) = evaluate_radiated_power(
                &pressure,
                &velocity,
                self.surfaces[boundary].areas(),
                Medium::from(self.sources[boundary].gas()),
            );
            if pressure.iter().chain(&velocity).any(|&value| !finite(value))
                || !finite(work.surface_power)
                || !work.apparent_power.is_finite()
            {
                return Err(bad("nonfinite molecular core-shell total traces or work"));
            }
            traces[boundary] = Some(CoreShellTrace {
                pressure, velocity, work, work_interval,
            });
        }
        let boundary = self.radiating_interface();
        let n = self.surfaces[boundary].areas().len();
        let (scattered_pressure, scattered_velocity) = if let Some(offset) = self.offsets[boundary] {
            (
                correction[offset..offset + n].to_vec(),
                correction[offset + n..offset + 2 * n]
                    .iter()
                    .map(|&value| physical_velocity(value, self.impedance))
                    .collect(),
            )
        } else {
            (vec![C64::ZERO; n], vec![C64::ZERO; n])
        };
        poll(gate)?;
        Ok(MolecularCoreShellScatteredField {
            interface: self,
            incident,
            maximum_source_point_products,
            traces,
            scattered_pressure,
            scattered_velocity,
        })
    }
}

/// Coherent refracted field retaining the three original molecular sources.
/// No lossless far-field power is inferred for an absorbing exterior.
pub struct MolecularCoreShellScatteredField<'p, 'a> {
    interface: &'p MolecularCoreShellFluidInterface<'a>,
    incident: IncidentField<'p>,
    maximum_source_point_products: u64,
    traces: [Option<CoreShellTrace>; 2],
    scattered_pressure: Vec<C64>,
    scattered_velocity: Vec<C64>,
}

impl MolecularCoreShellScatteredField<'_, '_> {
    /// Matched interfaces have no artificial unknowns. Their surrounding
    /// physical fields are still observable with strict region admission.
    pub fn interface_trace(&self, boundary: CoreShellInterface) -> Option<&CoreShellTrace> {
        self.traces[boundary.index()].as_ref()
    }
    pub fn radiating_surface(&self) -> &SpherePanels {
        self.interface.surfaces[self.interface.radiating_interface()]
    }
    pub fn scattered_pressure_pa(&self) -> &[C64] {
        &self.scattered_pressure
    }
    pub fn scattered_normal_velocity_m_s(&self) -> &[C64] {
        &self.scattered_velocity
    }
    /// Signed -Re(work) enclosed by an active boundary: core alone or shell
    /// plus core. Includes spatial error and is never clipped to passivity.
    /// A removed interface returns None rather than inventing a heat partition.
    pub fn enclosed_absorbed_power_w(&self, boundary: CoreShellInterface) -> Option<f64> {
        self.interface_trace(boundary)
            .map(|trace| -trace.work.surface_power.re)
    }
    /// Shell dissipation is Re(work_core - work_outer), because its outward
    /// normal on the core is reversed. Requires both physical boundary traces;
    /// subtracting two independent isolated-body absorptions is not equivalent.
    pub fn shell_absorbed_power_w(&self) -> Option<f64> {
        Some(
            self.traces[1].as_ref()?.work.surface_power.re
                - self.traces[0].as_ref()?.work.surface_power.re,
        )
    }
}
