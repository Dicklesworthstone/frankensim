//! Resolved radial heat diffusion through an actual finite-thickness solid wall.
//!
//! A source-bound, homogeneous solid at T0 obeys C_vol dT/dt = k Laplacian(T).
//! The original annular P1 geometry and native O(N) tridiagonal solver retain
//! every solid temperature coordinate. A unit gas-interface temperature gives
//! a reusable Dirichlet-to-Neumann response at one positive angular frequency.
//! The other face exchanges heat with the same background bath through an
//! explicitly declared isothermal, insulated or finite-conductance boundary.
//!
//! Gas/solid contact is perfect. The exact discrete solid admittance is
//! condensed into the EXISTING gas thermal-boundary solve, including its
//! original molecular Cp. Its same-frequency G/C representation is private:
//! public results distinguish solid Fourier loss, external bath exchange and
//! solid thermal storage. No uniform-through-thickness wall is substituted.
//!
//! This is a linear tangent about uniform equilibrium. Source density,
//! selected specific heat and radial conductivity are evaluated once at T0;
//! no isotropy, heat-capacity constraint, phase change, thermal expansion,
//! contact resistance or mechanical wall law is inferred. The solid is a
//! hollow cylindrical shell: an axis, axial heat diffusion and moving geometry
//! require their own owners. Positive losses are acoustic availability losses,
//! not an automatic assignment of mean caloric heating.

use super::{
    AnnularAdmission, AnnularSection, AnnularThermalWall, AnnularThermalWallField, DuctError,
    FluidDuctError, PreparedAnnularThermalWalls, SegmentWave, bad, filled, finite, nodal_load,
    poll, positive, radial_mass, relative,
};
use crate::compose::TwoPort;
use fs_la::eigen_complex::solve_tridiagonal_complex;
use fs_material::gas::GasState;
use fs_material::relaxation::RelaxationState;
use fs_material::thermal::transport::{
    PreparedThermalTransport, ThermalTransportError, ThermalTransportSample,
};
use fs_math::{c64::C64, dd::Dd, det};

/// Face touching the gas; the opposite face touches the declared bath.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RadialWallFace {
    Inner,
    Outer,
}
impl RadialWallFace {
    fn side(self) -> usize {
        usize::from(self == Self::Outer)
    }
}

/// External boundary at the same background temperature as the solid and gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RadialWallBath {
    Isothermal,
    Adiabatic,
    /// Explicit conductance per area of the BATH face, not the gas face.
    Conductance {
        conductance_w_m2_k: f64,
    },
}
impl RadialWallBath {
    fn conductance(self) -> Result<f64, DuctError> {
        let value = match self {
            Self::Isothermal | Self::Adiabatic => 0.0,
            Self::Conductance { conductance_w_m2_k } => conductance_w_m2_k,
        };
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(bad(
                "radial solid bath conductance must be finite and nonnegative",
            ))
        }
    }
}

/// Explicit radial resolution, algebraic and resource limits.
#[derive(Clone, Copy, Debug)]
pub struct RadialThermalWallAdmission {
    /// max cell width * sqrt(omega / diffusivity), in (0,1].
    /// This is a resolution ceiling, not a continuum-error certificate.
    pub maximum_cell_diffusion_phase: f64,
    /// Componentwise solve error, weak reaction/positive-form agreement and
    /// complex heat-balance defect, in (0,1). Refine the mesh separately.
    pub maximum_relative_residual: f64,
    /// Conservative numerical peak: 320*(elements+1)+4096 bytes, including
    /// borrowed radial geometry and native solve scratch. Excludes allocator
    /// metadata and the independently prepared source-material owner.
    pub maximum_payload_bytes: usize,
}

#[derive(Debug)]
pub enum RadialThermalWallError {
    Transport(ThermalTransportError),
    Fluid(FluidDuctError),
    Duct(DuctError),
}
impl From<ThermalTransportError> for RadialThermalWallError {
    fn from(error: ThermalTransportError) -> Self {
        Self::Transport(error)
    }
}
impl From<FluidDuctError> for RadialThermalWallError {
    fn from(error: FluidDuctError) -> Self {
        Self::Fluid(error)
    }
}
impl From<DuctError> for RadialThermalWallError {
    fn from(error: DuctError) -> Self {
        Self::Duct(error)
    }
}
impl core::fmt::Display for RadialThermalWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radial thermal wall: {self:?}")
    }
}
impl std::error::Error for RadialThermalWallError {}

/// Solid P1 fields, with exp(-i omega t) peak phasors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialThermalWallField {
    pub radius_m: f64,
    pub temperature_k: C64,
    pub radial_temperature_gradient_k_m: C64,
    /// Positive in the direction of increasing radius, at either contact orientation.
    pub radial_heat_flux_w_m2: C64,
    pub conduction_dissipation_w_m3: f64,
    /// C_vol |T'|²/(4 T0), distinct from irreversible work.
    pub thermal_storage_j_m3: f64,
}

/// Physical integrated and face quantities for a prescribed interface temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialThermalWallObservation {
    pub interface_radius_m: f64,
    pub bath_radius_m: f64,
    pub interface_temperature_k: C64,
    pub bath_temperature_k: C64,
    /// Directed FROM gas INTO solid, independent of which shell face is wetted.
    pub fluid_to_wall_heat_flux_w_m2: C64,
    /// Directed FROM solid INTO bath, per bath-face area.
    pub wall_to_bath_heat_flux_w_m2: C64,
    pub solid_conduction_w_m: f64,
    pub bath_exchange_w_m: f64,
    pub thermal_storage_j_m: f64,
    /// Exact P1 integral of the temperature phasor over solid cross-section.
    pub temperature_integral_k_m2: C64,
}

/// Source and actual shell retained with their one-frequency radial response.
/// Preparation is O(N), all subsequent cell/face/integral queries O(1), with
/// no allocation or source reevaluation.
#[derive(Debug)]
pub struct PreparedRadialThermalWall<'a> {
    section: &'a AnnularSection,
    source: &'a PreparedThermalTransport,
    sample: ThermalTransportSample,
    omega: f64,
    face: RadialWallFace,
    bath: RadialWallBath,
    // Native -i omega fields for unit interface temperature.
    temperature_response: Vec<C64>,
    admittance: C64,
    bath_flux: C64,
    temperature_integral: C64,
    loss_roots: [f64; 2],
    storage_root: f64,
    local_conduction_root: f64,
    local_storage_root: f64,
    relative_backward_residual: f64,
    relative_reaction_defect: f64,
    relative_heat_balance_defect: f64,
}

impl AnnularSection {
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_radial_thermal_wall<'a>(
        &'a self,
        source: &'a PreparedThermalTransport,
        temperature_k: f64,
        omega: f64,
        face: RadialWallFace,
        bath: RadialWallBath,
        admission: RadialThermalWallAdmission,
    ) -> Result<PreparedRadialThermalWall<'a>, RadialThermalWallError> {
        self.prepare_radial_thermal_wall_controlled(
            source,
            temperature_k,
            omega,
            face,
            bath,
            admission,
            || false,
        )
    }

    /// Polls before source evaluation/allocation, every 64 nodes and immediately
    /// before/after the bounded native tridiagonal tile (at most 4097 nodes).
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_radial_thermal_wall_controlled<'a>(
        &'a self,
        source: &'a PreparedThermalTransport,
        temperature_k: f64,
        omega: f64,
        face: RadialWallFace,
        bath: RadialWallBath,
        admission: RadialThermalWallAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRadialThermalWall<'a>, RadialThermalWallError> {
        poll(&mut cancelled, 0)?;
        if !positive(omega)
            || !positive(temperature_k)
            || !positive(admission.maximum_cell_diffusion_phase)
            || admission.maximum_cell_diffusion_phase > 1.0
            || !positive(admission.maximum_relative_residual)
            || admission.maximum_relative_residual >= 1.0
            || 320 * (self.elements() + 1) + 4096 > admission.maximum_payload_bytes
        {
            return Err(
                bad("invalid radial solid frequency, numerical limits or payload budget").into(),
            );
        }
        let conductance = bath.conductance()?;
        let sample = source.at_temperature(temperature_k)?;
        let conductivity = sample.conductivity_w_m_k();
        let capacity = sample.volumetric_heat_capacity_j_m3_k();
        let phase =
            self.gap_m * self.maximum_cell_fraction * det::sqrt(omega / sample.diffusivity_m2_s());
        let argument = omega / sample.diffusivity_m2_s() * self.gap_m * self.gap_m;
        let mean_radius = self.inner_m + 0.5 * self.gap_m;
        let radii = [self.inner_m, self.outer_m];
        let contact = face.side();
        let bath_side = 1 - contact;
        let boundary = conductance * self.gap_m / conductivity * (radii[bath_side] / mean_radius);
        if !positive(argument)
            || !phase.is_finite()
            || phase > admission.maximum_cell_diffusion_phase
            || !boundary.is_finite()
        {
            return Err(bad(
                "radial solid mesh or scaled diffusion is unresolved or unrepresentable",
            )
            .into());
        }
        let (mut nodes, backward) =
            solve_solid(self, argument, boundary, contact, bath, &mut cancelled)?;
        let n = self.elements();
        let (mut stiffness, mut mass, mut integral_re, mut integral_im) =
            (Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for i in 0..n {
            if i % 64 == 0 {
                poll(&mut cancelled, i)?;
            }
            stiffness = stiffness
                + Dd::from_f64(self.stiffness_edges[i])
                    * Dd::from_f64((nodes[i + 1] - nodes[i]).norm_sq());
            let x = self.fractions[i];
            let width = self.fractions[i + 1] - x;
            let left = self.inner_m / mean_radius + self.gap_m / mean_radius * x;
            let right = self.inner_m / mean_radius + self.gap_m / mean_radius * (x + width);
            let (diagonal, off) = radial_mass(width, left, right);
            let ratio = off / diagonal[0];
            let form = (nodes[i] + nodes[i + 1].scale(ratio)).norm_sq() * diagonal[0]
                + nodes[i + 1].norm_sq() * (diagonal[1] - off * ratio);
            mass = mass + Dd::from_f64(form);
        }
        for (i, value) in nodes.iter().enumerate() {
            if i % 64 == 0 {
                poll(&mut cancelled, i)?;
            }
            let weight = Dd::from_f64(nodal_load(self, i));
            integral_re = integral_re + weight * Dd::from_f64(value.re);
            integral_im = integral_im + weight * Dd::from_f64(value.im);
        }
        let contact_perimeter = core::f64::consts::TAU * radii[contact];
        let bath_perimeter = core::f64::consts::TAU * radii[bath_side];
        let fourier_form =
            conductivity * (self.area_m2 / self.gap_m / self.gap_m) * stiffness.to_f64();
        let bath_form = conductance * bath_perimeter * nodes[bath_side * n].norm_sq();
        let storage_form = capacity * self.area_m2 * mass.to_f64();
        // Green's identity gives the exact discrete admittance through positive
        // forms. In particular a tiny real loss is not obtained by subtracting
        // almost equal endpoint flux terms in the quasistatic insulated limit.
        let admittance = C64::new(
            (fourier_form + bath_form) / contact_perimeter,
            -omega * (storage_form / contact_perimeter),
        );
        let reaction = endpoint_reaction(self, argument, contact, &nodes)
            .conj()
            .scale(conductivity / self.gap_m * mean_radius / radii[contact]);
        let reaction_defect = relative(
            (admittance - reaction).abs(),
            admittance.abs() + reaction.abs(),
        );
        let bath_flux = match bath {
            RadialWallBath::Isothermal => -endpoint_reaction(self, argument, bath_side, &nodes)
                .conj()
                .scale(conductivity / self.gap_m * mean_radius / radii[bath_side]),
            _ => nodes[bath_side * n].conj().scale(conductance),
        };
        let integral = C64::new(integral_re.to_f64(), -integral_im.to_f64()).scale(self.area_m2);
        let input = admittance.scale(contact_perimeter);
        let output = bath_flux.scale(bath_perimeter);
        let stored_rate = C64::new(0.0, omega * capacity) * integral;
        let heat_defect = relative(
            (input - output + stored_rate).abs(),
            input.abs() + output.abs() + stored_rate.abs(),
        );
        let loss_roots = [
            det::sqrt(fourier_form / (2.0 * temperature_k)),
            det::sqrt(bath_form / (2.0 * temperature_k)),
        ];
        let storage_root = det::sqrt(storage_form / (4.0 * temperature_k));
        let local_conduction_root = det::sqrt(conductivity / (2.0 * temperature_k));
        let local_storage_root = det::sqrt(capacity / (4.0 * temperature_k));
        if ![admittance, reaction, bath_flux, integral]
            .into_iter()
            .all(finite)
            || admittance.re < 0.0
            || admittance.im >= 0.0
            || !loss_roots
                .into_iter()
                .chain([storage_root, local_conduction_root, local_storage_root])
                .all(|v| v.is_finite() && v >= 0.0)
            || backward > admission.maximum_relative_residual
            || reaction_defect > admission.maximum_relative_residual
            || heat_defect > admission.maximum_relative_residual
        {
            return Err(bad(
                "radial solid response or independent heat/work residual exceeds admission",
            )
            .into());
        }
        for (i, value) in nodes.iter_mut().enumerate() {
            if i % 64 == 0 {
                poll(&mut cancelled, i)?;
            }
            *value = value.conj();
        }
        poll(&mut cancelled, n)?;
        Ok(PreparedRadialThermalWall {
            section: self,
            source,
            sample,
            omega,
            face,
            bath,
            temperature_response: nodes,
            admittance,
            bath_flux,
            temperature_integral: integral,
            loss_roots,
            storage_root,
            local_conduction_root,
            local_storage_root,
            relative_backward_residual: backward,
            relative_reaction_defect: reaction_defect,
            relative_heat_balance_defect: heat_defect,
        })
    }
}

// Original area-normalized K+i*omega*C/k M; only boundary conditions change.
fn endpoint_reaction(section: &AnnularSection, argument: f64, side: usize, values: &[C64]) -> C64 {
    let n = section.elements();
    let edge = if side == 0 { 0 } else { n - 1 };
    let node = side * n;
    let other = if side == 0 { 1 } else { n - 1 };
    let diagonal = C64::new(
        section.stiffness_edges[edge],
        argument * (section.boundary_weights[side] - section.mass_edges[edge]),
    );
    let off = C64::new(
        -section.stiffness_edges[edge],
        argument * section.mass_edges[edge],
    );
    diagonal * values[node] + off * values[other]
}

fn solve_solid(
    section: &AnnularSection,
    argument: f64,
    boundary: f64,
    contact: usize,
    bath: RadialWallBath,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, f64), DuctError> {
    let n = section.elements();
    let scale = 1.0 / argument.max(boundary).max(1.0);
    let z = argument * scale;
    let mut diagonal = filled(n + 1, C64::ZERO)?;
    let mut off = filled(n, C64::ZERO)?;
    let mut values = filled(n + 1, C64::ZERO)?;
    for i in 0..=n {
        if i % 64 == 0 {
            poll(cancelled, i)?;
        }
        let (stiffness, mass) = if i == 0 {
            (
                section.stiffness_edges[0],
                section.boundary_weights[0] - section.mass_edges[0],
            )
        } else if i == n {
            (
                section.stiffness_edges[n - 1],
                section.boundary_weights[1] - section.mass_edges[n - 1],
            )
        } else {
            (
                section.stiffness_edges[i - 1] + section.stiffness_edges[i],
                section.mass_diagonal[i - 1],
            )
        };
        diagonal[i] = C64::new(scale * stiffness, z * mass);
        if i < n {
            off[i] = C64::new(
                -scale * section.stiffness_edges[i],
                z * section.mass_edges[i],
            );
        }
    }
    let contact_node = contact * n;
    let contact_edge = if contact == 0 { 0 } else { n - 1 };
    let neighbor = if contact == 0 { 1 } else { n - 1 };
    let eliminated_drive = -off[contact_edge];
    values[neighbor] = eliminated_drive;
    diagonal[contact_node] = C64::ONE;
    off[contact_edge] = C64::ZERO;
    values[contact_node] = C64::ONE;
    let bath_node = (1 - contact) * n;
    if bath == RadialWallBath::Isothermal {
        diagonal[bath_node] = C64::ONE;
        off[if bath_node == 0 { 0 } else { n - 1 }] = C64::ZERO;
    } else {
        diagonal[bath_node] = diagonal[bath_node] + C64::from_re(scale * boundary);
    }
    if diagonal
        .iter()
        .chain(&off)
        .chain(&values)
        .any(|&v| !finite(v))
    {
        return Err(bad("radial solid thermal matrix is unrepresentable"));
    }
    poll(cancelled, n)?;
    solve_tridiagonal_complex(&off, &diagonal, &off, &mut values, 1)
        .map_err(|_| DuctError::Singular)?;
    poll(cancelled, n)?;
    let mut residual = 0.0_f64;
    for i in 0..=n {
        if i % 64 == 0 {
            poll(cancelled, i)?;
        }
        let drive = if i == contact_node {
            C64::ONE
        } else if i == neighbor {
            eliminated_drive
        } else {
            C64::ZERO
        };
        let mut error = diagonal[i] * values[i] - drive;
        let mut magnitude = diagonal[i].abs() * values[i].abs() + drive.abs();
        if i > 0 {
            error = error + off[i - 1] * values[i - 1];
            magnitude += off[i - 1].abs() * values[i - 1].abs();
        }
        if i < n {
            error = error + off[i] * values[i + 1];
            magnitude += off[i].abs() * values[i + 1].abs();
        }
        residual = residual.max(relative(error.abs(), magnitude));
    }
    Ok((values, residual))
}

impl<'a> PreparedRadialThermalWall<'a> {
    pub const fn source(&self) -> &'a PreparedThermalTransport {
        self.source
    }
    pub const fn section(&self) -> &'a AnnularSection {
        self.section
    }
    pub const fn sample(&self) -> &ThermalTransportSample {
        &self.sample
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn gas_face(&self) -> RadialWallFace {
        self.face
    }
    pub const fn bath(&self) -> RadialWallBath {
        self.bath
    }
    /// Native heat-flux/temperature response at the original gas interface.
    pub const fn interface_admittance_w_m2_k(&self) -> C64 {
        self.admittance
    }
    pub const fn relative_backward_residual(&self) -> f64 {
        self.relative_backward_residual
    }
    pub const fn relative_reaction_defect(&self) -> f64 {
        self.relative_reaction_defect
    }
    pub const fn relative_heat_balance_defect(&self) -> f64 {
        self.relative_heat_balance_defect
    }

    pub fn field_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        interface_temperature_k: C64,
    ) -> Result<RadialThermalWallField, DuctError> {
        if !finite(interface_temperature_k)
            || cell >= self.section.elements()
            || !fraction.is_finite()
            || !(0.0..=1.0).contains(&fraction)
        {
            return Err(bad(
                "invalid radial solid field cell, fraction or interface temperature",
            ));
        }
        let x = self.section.fractions[cell];
        let width = self.section.fractions[cell + 1] - x;
        let a = self.temperature_response[cell];
        let b = self.temperature_response[cell + 1];
        let temperature = (a.scale(1.0 - fraction) + b.scale(fraction)) * interface_temperature_k;
        let gradient = (b - a) / C64::from_re(self.section.gap_m * width) * interface_temperature_k;
        let flux = gradient.scale(-self.sample.conductivity_w_m_k());
        let loss = gradient.scale(self.local_conduction_root).norm_sq();
        let storage = temperature.scale(self.local_storage_root).norm_sq();
        if ![temperature, gradient, flux].into_iter().all(finite)
            || !loss.is_finite()
            || !storage.is_finite()
        {
            return Err(bad("radial solid local field is unrepresentable"));
        }
        Ok(RadialThermalWallField {
            radius_m: self.section.inner_m + self.section.gap_m * (x + fraction * width),
            temperature_k: temperature,
            radial_temperature_gradient_k_m: gradient,
            radial_heat_flux_w_m2: flux,
            conduction_dissipation_w_m3: loss,
            thermal_storage_j_m3: storage,
        })
    }

    pub fn observation(
        &self,
        interface_temperature_k: C64,
    ) -> Result<RadialThermalWallObservation, DuctError> {
        if !finite(interface_temperature_k) {
            return Err(bad("radial solid interface temperature must be finite"));
        }
        let side = self.face.side();
        let radii = [self.section.inner_m, self.section.outer_m];
        let bath_temperature = self.temperature_response[(1 - side) * self.section.elements()]
            * interface_temperature_k;
        let input = self.admittance * interface_temperature_k;
        let output = self.bath_flux * interface_temperature_k;
        let integral = self.temperature_integral * interface_temperature_k;
        let losses = self
            .loss_roots
            .map(|r| interface_temperature_k.scale(r).norm_sq());
        let storage = interface_temperature_k.scale(self.storage_root).norm_sq();
        if ![bath_temperature, input, output, integral]
            .into_iter()
            .all(finite)
            || !losses.into_iter().chain([storage]).all(f64::is_finite)
        {
            return Err(bad("radial solid observation is unrepresentable"));
        }
        Ok(RadialThermalWallObservation {
            interface_radius_m: radii[side],
            bath_radius_m: radii[1 - side],
            interface_temperature_k,
            bath_temperature_k: bath_temperature,
            fluid_to_wall_heat_flux_w_m2: input,
            wall_to_bath_heat_flux_w_m2: output,
            solid_conduction_w_m: losses[0],
            bath_exchange_w_m: losses[1],
            thermal_storage_j_m: storage,
            temperature_integral_k_m2: integral,
        })
    }
}

/// Either an existing explicit boundary or a resolved original solid shell.
/// Immutable response lifetimes may be shortened to the common preparation borrow.
#[derive(Clone, Copy, Debug)]
pub enum AnnularConjugateWall<'a> {
    Boundary(AnnularThermalWall),
    Solid(&'a PreparedRadialThermalWall<'a>),
}

/// Face and integrated work for either an explicit boundary or a resolved solid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnularConjugateWallObservation {
    pub interface_radius_m: f64,
    pub bath_radius_m: f64,
    pub interface_temperature_k: C64,
    pub bath_face_temperature_k: C64,
    pub fluid_to_wall_heat_flux_w_m2: C64,
    /// Per physical bath-face area; bath and interface radii can differ.
    pub wall_to_bath_heat_flux_w_m2: C64,
    pub solid_conduction_w_m: f64,
    pub bath_exchange_w_m: f64,
    pub thermal_storage_j_m: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AnnularConjugateWallDissipation {
    pub gas_viscous_w_m: f64,
    pub gas_conduction_w_m: f64,
    pub gas_molecular_w_m: f64,
    pub solid_conduction_w_m: [f64; 2],
    pub bath_exchange_w_m: [f64; 2],
    pub total_w_m: f64,
}

/// The existing gas solve with exact discrete solid boundary elimination.
/// The internally equivalent boundary coefficients are deliberately private.
#[derive(Debug)]
pub struct PreparedConjugateAnnularWalls<'g, 'm, 'w, M: ?Sized> {
    gas: PreparedAnnularThermalWalls<'g, 'm, M>,
    walls: [AnnularConjugateWall<'w>; 2],
}

fn condensed_walls(
    section: &AnnularSection,
    gas: &GasState,
    omega: f64,
    walls: [AnnularConjugateWall<'_>; 2],
) -> Result<[AnnularThermalWall; 2], DuctError> {
    let mut result = [AnnularThermalWall::Adiabatic; 2];
    for (side, wall) in walls.into_iter().enumerate() {
        result[side] = match wall {
            AnnularConjugateWall::Boundary(boundary) => boundary,
            AnnularConjugateWall::Solid(solid) => {
                let expected_face = if side == 0 {
                    RadialWallFace::Outer
                } else {
                    RadialWallFace::Inner
                };
                let solid_radius = if solid.face == RadialWallFace::Inner {
                    solid.section.inner_m
                } else {
                    solid.section.outer_m
                };
                let gas_radius = if side == 0 {
                    section.inner_m
                } else {
                    section.outer_m
                };
                if solid.face != expected_face
                    || solid_radius != gas_radius
                    || solid.omega != omega
                    || solid.sample.temperature_k() != gas.temperature
                {
                    return Err(bad(
                        "resolved solid face, radius, frequency or equilibrium temperature does not match gas",
                    ));
                }
                let capacity = -solid.admittance.im / omega;
                if !positive(capacity) {
                    return Err(bad("condensed solid storage is unrepresentable"));
                }
                AnnularThermalWall::Finite {
                    conductance_w_m2_k: solid.admittance.re,
                    areal_heat_capacity_j_m2_k: capacity,
                }
            }
        };
    }
    Ok(result)
}

impl AnnularSection {
    pub fn prepare_conjugate_wall_fields<'g, 'm, 'w>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        walls: [AnnularConjugateWall<'w>; 2],
        admission: AnnularAdmission,
    ) -> Result<PreparedConjugateAnnularWalls<'g, 'm, 'w, GasState>, RadialThermalWallError> {
        self.prepare_conjugate_wall_fields_controlled(gas, omega, walls, admission, || false)
    }

    /// The gas budget excludes separately admitted retained solid/material owners.
    pub fn prepare_conjugate_wall_fields_controlled<'g, 'm, 'w>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        walls: [AnnularConjugateWall<'w>; 2],
        admission: AnnularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedConjugateAnnularWalls<'g, 'm, 'w, GasState>, RadialThermalWallError> {
        poll(&mut cancelled, 0)?;
        let boundaries = condensed_walls(self, gas, omega, walls)?;
        let gas = self
            .prepare_thermal_wall_fields_controlled(gas, omega, boundaries, admission, cancelled)?;
        Ok(PreparedConjugateAnnularWalls { gas, walls })
    }

    pub fn prepare_molecular_conjugate_wall_fields<'g, 'm, 'w>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        walls: [AnnularConjugateWall<'w>; 2],
        admission: AnnularAdmission,
    ) -> Result<PreparedConjugateAnnularWalls<'g, 'm, 'w, RelaxationState>, RadialThermalWallError>
    {
        self.prepare_molecular_conjugate_wall_fields_controlled(
            material,
            omega,
            walls,
            admission,
            || false,
        )
    }

    pub fn prepare_molecular_conjugate_wall_fields_controlled<'g, 'm, 'w>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        walls: [AnnularConjugateWall<'w>; 2],
        admission: AnnularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedConjugateAnnularWalls<'g, 'm, 'w, RelaxationState>, RadialThermalWallError>
    {
        poll(&mut cancelled, 0)?;
        let boundaries = condensed_walls(self, material.gas(), omega, walls)?;
        let gas = self.prepare_molecular_thermal_wall_fields_controlled(
            material, omega, boundaries, admission, cancelled,
        )?;
        Ok(PreparedConjugateAnnularWalls { gas, walls })
    }
}

impl<'g, 'm, 'w, M: ?Sized> PreparedConjugateAnnularWalls<'g, 'm, 'w, M> {
    pub fn source(&self) -> &'m M {
        self.gas.source()
    }
    pub fn section(&self) -> &'g AnnularSection {
        self.gas.section()
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.gas.angular_frequency_rad_s()
    }
    pub const fn walls(&self) -> &[AnnularConjugateWall<'w>; 2] {
        &self.walls
    }
    pub fn wave(&self) -> &SegmentWave {
        self.gas.wave()
    }
    pub fn mean_temperature_response(&self) -> C64 {
        self.gas.mean_temperature_response()
    }
    pub fn relative_backward_residual(&self) -> f64 {
        self.gas.relative_backward_residual()
    }
    pub fn relative_thermal_work_defect(&self) -> f64 {
        self.gas.relative_thermal_work_defect()
    }
    pub fn two_port(&self, length_m: f64, reference: f64) -> Result<TwoPort, DuctError> {
        self.gas.two_port(length_m, reference)
    }
    pub fn field_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<AnnularThermalWallField, DuctError> {
        self.gas
            .field_in_cell(cell, fraction, pressure_pa, volume_flow_m3_s)
    }
    pub fn solid_field_in_cell(
        &self,
        side: usize,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
    ) -> Result<RadialThermalWallField, DuctError> {
        let interface = self.gas.wall_observation(side, pressure_pa)?;
        match self.walls[side] {
            AnnularConjugateWall::Solid(solid) => {
                solid.field_in_cell(cell, fraction, interface.temperature_k)
            }
            AnnularConjugateWall::Boundary(_) => {
                Err(bad("this annular wall has no resolved solid field"))
            }
        }
    }
    pub fn wall_observation(
        &self,
        side: usize,
        pressure_pa: C64,
    ) -> Result<AnnularConjugateWallObservation, DuctError> {
        let fluid = self.gas.wall_observation(side, pressure_pa)?;
        let result = match self.walls[side] {
            AnnularConjugateWall::Boundary(_) => AnnularConjugateWallObservation {
                interface_radius_m: fluid.radius_m,
                bath_radius_m: fluid.radius_m,
                interface_temperature_k: fluid.temperature_k,
                bath_face_temperature_k: fluid.temperature_k,
                fluid_to_wall_heat_flux_w_m2: fluid.fluid_to_wall_heat_flux_w_m2,
                wall_to_bath_heat_flux_w_m2: fluid.wall_to_bath_heat_flux_w_m2,
                solid_conduction_w_m: 0.0,
                bath_exchange_w_m: fluid.bath_dissipation_w_m,
                thermal_storage_j_m: fluid.thermal_storage_j_m,
            },
            AnnularConjugateWall::Solid(solid) => {
                let response = solid.observation(fluid.temperature_k)?;
                AnnularConjugateWallObservation {
                    interface_radius_m: response.interface_radius_m,
                    bath_radius_m: response.bath_radius_m,
                    interface_temperature_k: response.interface_temperature_k,
                    bath_face_temperature_k: response.bath_temperature_k,
                    fluid_to_wall_heat_flux_w_m2: response.fluid_to_wall_heat_flux_w_m2,
                    wall_to_bath_heat_flux_w_m2: response.wall_to_bath_heat_flux_w_m2,
                    solid_conduction_w_m: response.solid_conduction_w_m,
                    bath_exchange_w_m: response.bath_exchange_w_m,
                    thermal_storage_j_m: response.thermal_storage_j_m,
                }
            }
        };
        Ok(result)
    }
    pub fn dissipation_per_axial_length(
        &self,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<AnnularConjugateWallDissipation, DuctError> {
        let gas = self
            .gas
            .dissipation_per_axial_length(pressure_pa, volume_flow_m3_s)?;
        let mut solid = [0.0; 2];
        let mut bath = [0.0; 2];
        for side in 0..2 {
            let wall = self.wall_observation(side, pressure_pa)?;
            solid[side] = wall.solid_conduction_w_m;
            bath[side] = wall.bath_exchange_w_m;
        }
        let total = gas.viscous_w_m
            + gas.conduction_w_m
            + gas.molecular_w_m
            + solid.into_iter().sum::<f64>()
            + bath.into_iter().sum::<f64>();
        if !total.is_finite() {
            return Err(bad(
                "conjugate annular integrated acoustic loss is unrepresentable",
            ));
        }
        Ok(AnnularConjugateWallDissipation {
            gas_viscous_w_m: gas.viscous_w_m,
            gas_conduction_w_m: gas.conduction_w_m,
            gas_molecular_w_m: gas.molecular_w_m,
            solid_conduction_w_m: solid,
            bath_exchange_w_m: bath,
            total_w_m: total,
        })
    }
}
