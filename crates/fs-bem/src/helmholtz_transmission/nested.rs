//! Coupled refraction through one curved fluid shell and its enclosed core.
//!
//! Both retained surfaces point outward from their enclosed volumes. The three
//! original homogeneous fluids are exterior (0), shell (1), and core (2).
//! Pressure and physical normal velocity are shared across each interface;
//! shell propagation uses its OWN density and wavenumber on both surfaces.
//!
//! With exp(-i omega t), use one pressure coordinate w=i*Z0*v_n, Z0=rho0*c0,
//! and beta_r=k0*rho_r/rho0. The original local transmission operator is
//! M_j=I+C_j, assembled with k_reference=beta_parent(j). For a cross-surface
//! block E_ij^r=[D_r,-beta_r*S_r;N_r/beta_parent(i),
//! -(beta_r/beta_parent(i))*D'_r], the complete system is
//! [M_outer,-E_outer,core^1;E_core,outer^1,M_core] x=[incident_outer;0].
//! The opposite signs follow from the shell representation -T_outer+T_core,
//! T=D*p-S*q. One factorization retains every reflection between the surfaces.
//!
//! We solve corrections to the original exterior incident field, analytically
//! cancelling its regular Green identities before assembling the forcing.
//! Cross-kernel differences reuse the original weak-contrast kernel owner.
//! Exactly matched adjacent media remove that interface from the linear
//! system; no numerical representation of a fictitious surface can scatter.
//!
//! The scope is stationary, inviscid, nondispersive positive-density fluids,
//! one convex closed outer surface, and one source-free enclosed core with
//! positive separation. Constant panels, equivalent-disc self terms, mesh
//! refinement needs and source authority remain those of the original BEM.

use super::*;

#[path = "nested/fields.rs"]
mod fields;
#[path = "nested/geometry.rs"]
mod geometry;

/// The original physical interface, with its normal pointing outward from
/// the volume it encloses. The shell's outward normal on `Core` is the opposite.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoreShellInterface {
    Outer,
    Core,
}
impl CoreShellInterface {
    const fn index(self) -> usize {
        match self {
            Self::Outer => 0,
            Self::Core => 1,
        }
    }
}

/// A physical observation region, never selected from an approximate radius.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoreShellRegion {
    Exterior,
    Shell,
    Core,
}

/// One frequency and reusable, fully coupled factorization for two nested
/// interfaces. Original sources and meshes remain borrowed. Source snapshots
/// use the same exact nonzero-sized trait-object-view reuse as `FluidInterface`.
pub struct CoreShellFluidInterface<'a> {
    surfaces: [&'a SpherePanels; 2],
    sources: [&'a dyn FluidAcoustics; 3],
    states: [FluidAcousticState; 3],
    omega: f64,
    k: [f64; 3],
    beta: [f64; 2],
    impedance: f64,
    ppw: [f64; 2],
    offsets: [Option<usize>; 2],
    dimension: usize,
    drive: Vec<C64>,
    factor: Option<LuComplex>,
    condition: f64,
}

impl<'a> CoreShellFluidInterface<'a> {
    /// `materials` are ordered exterior, shell, core. The outer mesh must be
    /// convex; every core vertex must lie strictly inside all outer face
    /// planes by the requested metric clearance and floating-point guard.
    /// Both meshes retain the original closed/outward/non-self-intersecting
    /// input duty. No sphere fit, mesh repair or nesting guess is introduced.
    ///
    /// The dense budget conservatively admits both interfaces before any
    /// quadratic allocation, including the original worker-dependent LU
    /// workspace. Exact transparency may subsequently reduce the actual solve.
    pub fn new(
        outer: &'a SpherePanels,
        core: &'a SpherePanels,
        materials: [&'a dyn FluidAcoustics; 3],
        omega_rad_s: f64,
        minimum_interface_clearance_m: f64,
        options: Options,
        gate: &CancelGate,
    ) -> Result<Self, TransmissionError> {
        poll(gate)?;
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(bad("core-shell acoustics requires positive finite frequency").into());
        }
        let workers = geometry::admit(outer, core, minimum_interface_clearance_m, options, gate)?;
        let exterior = materials[0].acoustic_state()?;
        poll(gate)?;
        let shell = acoustic_snapshot_for_view(materials[1], [(materials[0], exterior)])?;
        poll(gate)?;
        let core_state = acoustic_snapshot_for_view(
            materials[2],
            [(materials[0], exterior), (materials[1], shell)],
        )?;
        let states = [exterior, shell, core_state];
        let k = states.map(|state| omega_rad_s / state.sound_speed_m_s());
        let impedance = exterior.characteristic_impedance_pa_s_m();
        let density_ratio = shell.density_kg_m3() / exterior.density_kg_m3();
        let beta = [k[0], k[0] * density_ratio];
        if k.iter().chain(&beta).any(|v| !v.is_finite() || *v <= 0.0)
            || states.iter().any(|state| {
                let momentum = omega_rad_s * state.density_kg_m3();
                !momentum.is_finite() || momentum <= 0.0
            })
        {
            return Err(bad("unrepresentable core-shell acoustic scales").into());
        }
        let surfaces = [outer, core];
        let mut ppw = [
            admit_resolution(outer, k[0])?.min(admit_resolution(outer, k[1])?),
            admit_resolution(core, k[1])?.min(admit_resolution(core, k[2])?),
        ];
        let mut offsets = [None; 2];
        let mut dimension = 0;
        for boundary in 0..2 {
            if medium(states[boundary]) != medium(states[boundary + 1]) {
                offsets[boundary] = Some(dimension);
                dimension += 2 * surfaces[boundary].areas().len();
            }
        }
        if offsets.iter().all(Option::is_some) {
            // Contrast forcing samples the exterior reference wave on the
            // inner surface as well. Resolve that reference wavelength, even
            // when both actual adjacent media have longer wavelengths.
            ppw[1] = ppw[1].min(admit_resolution(core, k[0])?);
        }
        let mut owner = Self {
            surfaces,
            sources: materials,
            states,
            omega: omega_rad_s,
            k,
            beta,
            impedance,
            ppw,
            offsets,
            dimension,
            drive: Vec::new(),
            factor: None,
            condition: 1.0,
        };
        if dimension != 0 {
            let (operator, drive) = owner.assemble(gate)?;
            let (operator, factor, condition) =
                factor_contrast(operator, dimension, workers, gate)?;
            drop(operator);
            owner.drive = drive;
            owner.factor = Some(factor);
            owner.condition = condition;
        }
        poll(gate)?;
        Ok(owner)
    }

    pub fn surface(&self, interface: CoreShellInterface) -> &'a SpherePanels {
        self.surfaces[interface.index()]
    }
    /// Original exterior, shell and core owners, in that order.
    pub fn material_sources(&self) -> [&'a dyn FluidAcoustics; 3] {
        self.sources
    }
    /// The coherent snapshots actually used by this frequency preparation.
    pub fn material_states(&self) -> [FluidAcousticState; 3] {
        self.states
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Matched interfaces remain part of region/receiver geometry but have no
    /// artificial boundary unknowns or outgoing representation.
    pub fn active_interface_count(&self) -> usize {
        self.offsets
            .iter()
            .filter(|offset| offset.is_some())
            .count()
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

    fn assemble(&self, gate: &CancelGate) -> Result<(Vec<C64>, Vec<C64>), HelmholtzError> {
        let m = self.dimension;
        let mut operator = zero_matrix(m * m, "core-shell coupled matrix")?;
        let mut drive = zero_matrix(m * m, "core-shell contrast forcing")?;
        for boundary in 0..2 {
            let Some(offset) = self.offsets[boundary] else {
                continue;
            };
            let outer_density = self.states[boundary].density_kg_m3();
            let inner_density = self.states[boundary + 1].density_kg_m3();
            let ratio = inner_density / outer_density;
            let contrast = (inner_density - outer_density) / outer_density;
            if !ratio.is_finite() || ratio <= 0.0 || !contrast.is_finite() {
                return Err(bad("unrepresentable core-shell density contrast"));
            }
            let local = assemble_wave(
                self.surfaces[boundary],
                self.k[boundary],
                self.k[boundary + 1],
                self.beta[boundary],
                ratio,
                contrast,
                gate,
                kernel_difference,
                self_contrast,
            )?;
            let size = 2 * self.surfaces[boundary].areas().len();
            for row in 0..size {
                poll(gate)?;
                for col in 0..size {
                    let target = (offset + row) * m + offset + col;
                    operator[target] = local[row * size + col];
                    drive[target] = -local[row * size + col];
                }
            }
        }
        if self.offsets.iter().all(Option::is_some) {
            self.assemble_shell_coupling(&mut operator, &mut drive, gate)?;
        }
        if operator.iter().chain(&drive).any(|&value| !finite(value)) {
            return Err(bad("nonfinite core-shell coupled operator or forcing"));
        }
        Ok((operator, drive))
    }

    fn assemble_shell_coupling(
        &self,
        operator: &mut [C64],
        drive: &mut [C64],
        gate: &CancelGate,
    ) -> Result<(), HelmholtzError> {
        let m = self.dimension;
        let rho0 = self.states[0].density_kg_m3();
        let rho1 = self.states[1].density_kg_m3();
        let density_change = (rho1 - rho0) / rho0;
        for target in 0..2 {
            let source = 1 - target;
            let target_offset = self.offsets[target].expect("two active interfaces");
            let source_offset = self.offsets[source].expect("two active interfaces");
            let a = self.surfaces[target];
            let b = self.surfaces[source];
            let (na, nb) = (a.areas().len(), b.areas().len());
            let sign = if target == 0 { -1.0 } else { 1.0 };
            let beta = self.beta[target];
            for i in 0..na {
                poll(gate)?;
                for j in 0..nb {
                    if j % 64 == 0 {
                        poll(gate)?;
                    }
                    let (x, nx) = (a.centroids()[i], a.normals()[i]);
                    let (y, ny) = (b.centroids()[j], b.normals()[j]);
                    let (s, d, dp, n) = kernels(self.k[1], x, nx, y, ny);
                    let (s0, _, dp0, _) = kernels(self.k[0], x, nx, y, ny);
                    let [ds, dd, ddp, dn] = kernel_difference(self.k[0], self.k[1], x, nx, y, ny);
                    let area = b.areas()[j];
                    let transfer = [
                        d,
                        s.scale(-self.beta[1]),
                        n.scale(1.0 / beta),
                        dp.scale(-self.beta[1] / beta),
                    ];
                    // E_shell-E_exterior, with the TARGET derivative row's
                    // scale unchanged. Original differences cancel static
                    // singular terms and retain representable weak contrasts.
                    let difference = [
                        dd,
                        -s0.scale(self.beta[0] * density_change) - ds.scale(self.beta[1]),
                        dn.scale(1.0 / beta),
                        -dp0.scale((self.beta[0] / beta) * density_change)
                            - ddp.scale(self.beta[1] / beta),
                    ];
                    let row = target_offset + i;
                    let col = source_offset + j;
                    let entries = [
                        row * m + col,
                        row * m + col + nb,
                        (row + na) * m + col,
                        (row + na) * m + col + nb,
                    ];
                    for entry in 0..4 {
                        operator[entries[entry]] = transfer[entry].scale(sign * area);
                        drive[entries[entry]] = difference[entry].scale(-sign * area);
                    }
                }
            }
        }
        // The reference wave obeys exterior momentum. At the inner target its
        // derivative Green identity therefore returns (rho0/rho1)*w, not w.
        // Factor the missing density contribution instead of subtracting one.
        let core_offset = self.offsets[1].expect("two active interfaces");
        let n = self.surfaces[1].areas().len();
        let density_correction = (rho1 - rho0) / rho1;
        for i in 0..n {
            let diagonal = (core_offset + n + i) * m + core_offset + n + i;
            drive[diagonal] = drive[diagonal] - C64::from_re(density_correction);
        }
        Ok(())
    }

    /// Retain one original coherent incident field, regular throughout the
    /// entire outer volume. Analytic incident traces on the core are only the
    /// contrast-forcing reference; the actual core receives the coupled field.
    /// Repeated illuminations reuse the factorization and sample no materials.
    pub fn scatter_field<'p>(
        &'p self,
        incident: IncidentField<'p>,
        minimum_source_clearance_m: f64,
        gate: &CancelGate,
    ) -> Result<CoreShellScatteredField<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        let exterior = medium(self.states[0]);
        if !minimum_source_clearance_m.is_finite() || minimum_source_clearance_m <= 0.0 {
            return Err(bad("incident source clearance must be positive and finite"));
        }
        incident.validate()?;
        let positions = incident.source_positions();
        if !positions.is_empty() {
            near_field::Geometry::new_iterative(
                self.surfaces[0],
                &positions,
                minimum_source_clearance_m,
                gate,
            )?;
        }
        let reference = [
            IncidentTraces::from_admitted_field(
                self.surfaces[0],
                self.k[0],
                exterior,
                incident,
                Some(gate),
            )?,
            IncidentTraces::from_admitted_field(
                self.surfaces[1],
                self.k[0],
                exterior,
                incident,
                Some(gate),
            )?,
        ];
        let m = self.dimension;
        let mut incoming = zero_matrix(m, "core-shell incident traces")?;
        for (boundary, traces) in reference.iter().enumerate() {
            if let Some(offset) = self.offsets[boundary] {
                let n = self.surfaces[boundary].areas().len();
                incoming[offset..offset + n].copy_from_slice(traces.pressure_pa());
                for (j, &velocity) in traces.normal_velocity_m_s().iter().enumerate() {
                    incoming[offset + n + j] = velocity * C64::new(0.0, self.impedance);
                }
            }
        }
        if incoming.iter().any(|&value| !finite(value)) {
            return Err(bad("nonfinite pressure-scaled core-shell incident field"));
        }
        let mut correction = zero_matrix(m, "core-shell trace correction")?;
        if let Some(factor) = &self.factor {
            for (row, value) in correction.iter_mut().enumerate() {
                poll(gate)?;
                for (col, &input) in incoming.iter().enumerate() {
                    *value = *value + self.drive[row * m + col] * input;
                }
            }
            if correction.iter().any(|&value| !finite(value)) {
                return Err(bad("nonfinite core-shell contrast drive"));
            }
            factor
                .solve_with_cancel(&mut correction, || gate.is_requested())
                .map_err(lu_error)?;
        }
        if correction.iter().any(|&value| !finite(value)) {
            return Err(bad("nonfinite core-shell scattered solution"));
        }
        let mut traces = [None, None];
        for boundary in 0..2 {
            let Some(offset) = self.offsets[boundary] else {
                continue;
            };
            let n = self.surfaces[boundary].areas().len();
            let pressure: Vec<_> = (0..n)
                .map(|i| reference[boundary].pressure_pa()[i] + correction[offset + i])
                .collect();
            let velocity: Vec<_> = (0..n)
                .map(|i| {
                    reference[boundary].normal_velocity_m_s()[i]
                        + physical_velocity(correction[offset + n + i], self.impedance)
                })
                .collect();
            let (_, interval, work) = evaluate_radiated_power(
                &pressure,
                &velocity,
                self.surfaces[boundary].areas(),
                medium(self.states[boundary]),
            );
            if pressure
                .iter()
                .chain(&velocity)
                .any(|&value| !finite(value))
                || !finite(work.surface_power)
                || !work.apparent_power.is_finite()
            {
                return Err(bad("nonfinite core-shell total interface traces"));
            }
            traces[boundary] = Some(CoreShellTrace {
                pressure,
                velocity,
                work,
                work_interval: interval,
            });
        }
        let boundary = self.radiating_interface();
        let n = self.surfaces[boundary].areas().len();
        let (pressure, velocity) = if let Some(offset) = self.offsets[boundary] {
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
        let scattered = self.radiation(boundary, 0, pressure, velocity)?;
        poll(gate)?;
        Ok(CoreShellScatteredField {
            interface: self,
            incident,
            traces,
            scattered,
        })
    }

    fn radiation(
        &self,
        boundary: usize,
        region: usize,
        pressure: Vec<C64>,
        velocity: Vec<C64>,
    ) -> Result<RadiationSolution, HelmholtzError> {
        let surface = self.surfaces[boundary];
        let fluid = medium(self.states[region]);
        let (power, interval, diagnostics) =
            evaluate_radiated_power(&pressure, &velocity, surface.areas(), fluid);
        if !power.is_finite()
            || !finite(diagnostics.surface_power)
            || !diagnostics.apparent_power.is_finite()
        {
            return Err(bad("nonfinite core-shell representation work"));
        }
        Ok(RadiationSolution {
            pressure,
            velocity,
            k: self.k[region],
            medium: fluid,
            surface_fingerprint: surface_fingerprint(surface),
            panels_per_wavelength: self.ppw[boundary],
            radiated_power: power,
            radiated_power_roundoff_interval: interval,
            power_diagnostics: diagnostics,
            condition_lower_bound: self.condition,
            dense_cap_utilization: self.dimension as f64 / MAX_DENSE_PANELS as f64,
        })
    }
}

/// Total physical traces on one active interface. Both neighboring fluids use
/// this same pressure and this same outward normal velocity. Real work defects
/// are retained as discretization diagnostics, never relabeled absorption.
pub struct CoreShellTrace {
    pressure: Vec<C64>,
    velocity: Vec<C64>,
    work: RadiationPowerDiagnostics,
    work_interval: (f64, f64),
}
impl CoreShellTrace {
    pub fn pressure_pa(&self) -> &[C64] {
        &self.pressure
    }
    pub fn normal_velocity_m_s(&self) -> &[C64] {
        &self.velocity
    }
    pub fn interface_work(&self) -> &RadiationPowerDiagnostics {
        &self.work
    }
    pub fn work_roundoff_interval(&self) -> (f64, f64) {
        self.work_interval
    }
}

/// Complete nested response to the retained coherent incident field.
pub struct CoreShellScatteredField<'p, 'a> {
    interface: &'p CoreShellFluidInterface<'a>,
    incident: IncidentField<'p>,
    traces: [Option<CoreShellTrace>; 2],
    scattered: RadiationSolution,
}
impl CoreShellScatteredField<'_, '_> {
    /// A matched interface is eliminated, so it has no artificial trace
    /// unknowns. Its physical surrounding field remains observable by region.
    pub fn interface_trace(&self, boundary: CoreShellInterface) -> Option<&CoreShellTrace> {
        self.traces[boundary.index()].as_ref()
    }
    /// Original surface supporting the exterior outgoing representation. This
    /// is the core only when the outer interface is exactly transparent.
    pub fn radiating_surface(&self) -> &SpherePanels {
        self.interface.surfaces[self.interface.radiating_interface()]
    }
    pub fn scattered(&self) -> &RadiationSolution {
        &self.scattered
    }

    /// Original retained-triangle far-field observer of the complete nested
    /// solution. No independent-body far fields are summed.
    pub fn far_field(
        &self,
        directions: &[[f64; 3]],
        maximum_kernel_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<Vec<C64>, HelmholtzError> {
        super::super::far_field_integrated_cancellable(
            self.radiating_surface(),
            &self.scattered,
            medium(self.interface.states[0]),
            directions,
            maximum_kernel_evaluations,
            || gate.is_requested(),
        )
    }
}

fn physical_velocity(value: C64, impedance: f64) -> C64 {
    C64::new(value.im / impedance, -value.re / impedance)
}
fn zero_matrix(entries: usize, what: &'static str) -> Result<Vec<C64>, HelmholtzError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(entries)
        .map_err(|_| HelmholtzError::AllocationFailed { what })?;
    values.resize(entries, C64::ZERO);
    Ok(values)
}
