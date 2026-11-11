//! Full pressure fields in passive, spatially varying equivalent fluids.
//!
//! In exp(+i omega t), alpha = 1/rho_eff and beta = 1/K_eff give
//! div(alpha grad p) + omega² beta p = 0. Both coefficients remain INSIDE
//! their original volume forms. This is the same inhomogeneous pressure
//! equation as the lossless owner, with complex constitutive coefficients:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.164.html>.
//! JCA/JCAL material laws belong to fs-material, not this numerical owner.
//!
//! Loss is evaluated from its two positive local forms, independently of the
//! power balance: Im(alpha) |grad p|²/(2 omega) and
//! -omega Im(beta) |p|²/2. These are dynamic-density and compression losses;
//! their microscopic mechanism and destination remain material-owner data.
//! A sampled dispersive law alone does not determine stored energy or a causal
//! transient realization. Neither is inferred from the real coefficients.

use super::super::super::HelmholtzFe;
use super::super::super::pressure_release::PressureReleaseConstraints;
use super::super::MAX_QUADRATURE_ORDER;
use super::super::source::PressureSourceOwner;
use super::super::time::mechanical::{MechanicalBoundary, traction::WallTractionProfile};
use super::super::{BoundaryVelocitySource, VolumeVelocitySource};
use super::frequency::coupled::{
    CoupledPencil, MechanicalHarmonicLimits, MechanicalWallAccounting, MechanicalWallPhasors,
    MechanicalWallPower, WallTraction,
};
use super::mechanical::mechanical_walls_for;
use super::{
    BoundaryAdmittance, HarmonicLimits, HeterogeneousError, HeterogeneousHelmholtzFe,
    MAX_HARMONIC_DOFS, PreparedHeterogeneous, VelocityTrace, finite, invalid, poll,
    quadratic_nonnegative,
};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use crate::highorder::simplex::duffy_quadrature;
use fs_math::{c64::C64, dd::Dd};
use fs_sparse::{Coo, Csr};

mod sources;

/// Coherent effective material response at the exact preparation frequency.
///
/// Density includes the material owner's exposed-area/porosity convention.
/// For an equivalent porous fluid its velocity is superficial (Darcy) velocity,
/// so multiplying by porosity again would violate interface continuity.
/// This lane admits positive-real passive coefficients, including their
/// lossless limit; negative-parameter metamaterials need a separate model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquivalentFluidCoefficients {
    /// Effective density [kg/m³], with positive real and nonpositive imaginary part.
    pub dynamic_density_kg_m3: C64,
    /// Effective bulk modulus [Pa], with positive real and nonnegative imaginary part.
    pub dynamic_bulk_modulus_pa: C64,
}

impl EquivalentFluidCoefficients {
    fn inverse(self) -> Result<(C64, C64), HeterogeneousError> {
        let rho = self.dynamic_density_kg_m3;
        let bulk = self.dynamic_bulk_modulus_pa;
        if !finite(rho)
            || rho.re <= 0.0
            || rho.im > 0.0
            || !finite(bulk)
            || bulk.re <= 0.0
            || bulk.im < 0.0
        {
            return Err(invalid(
                "positive-real passive equivalent-fluid density and modulus required",
            )
            .into());
        }
        let alpha = rho.recip();
        let beta = bulk.recip();
        if !finite(alpha)
            || alpha.re <= 0.0
            || alpha.im < 0.0
            || !finite(beta)
            || beta.re <= 0.0
            || beta.im > 0.0
        {
            return Err(invalid("unrepresentable equivalent-fluid inverse coefficients").into());
        }
        Ok((alpha, beta))
    }
}

/// Independent spatial sampling, retained local loss and full-space solve limits.
#[derive(Clone, Copy, Debug)]
pub struct EquivalentFluidControls {
    /// Duffy points per axis, from pressure order + 2 through 16. The original
    /// curved map is evaluated at these same physical quadrature locations.
    pub quadrature_order: usize,
    /// Total volume samples, checked before any material callback.
    pub max_volume_points: usize,
    /// Retained per-element loss matrices, indices and descriptors in bytes.
    /// Sparse assembly storage and one element's temporary matrices are additional.
    pub max_element_loss_bytes: usize,
    /// Original dense factorization and componentwise physical residual limits.
    /// With mechanical walls, `max_dofs` bounds ALL free pressure, displacement,
    /// velocity and relaxation-memory coordinates, not just the pressure block.
    pub harmonic: HarmonicLimits,
}

impl Default for EquivalentFluidControls {
    fn default() -> Self {
        Self {
            quadrature_order: 4,
            max_volume_points: 2_000_000,
            max_element_loss_bytes: 64 * 1024 * 1024,
            harmonic: HarmonicLimits::default(),
        }
    }
}

/// Original tetrahedron's volume and independently integrated absorbed powers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquivalentFluidElementPower {
    pub volume_m3: f64,
    /// integral(Im(1/rho_eff) |grad p|²)/(2 omega) [W].
    pub density_dissipation_w: f64,
    /// -omega integral(Im(1/K_eff) |p|²)/2 [W].
    pub compressional_dissipation_w: f64,
}

impl EquivalentFluidElementPower {
    #[must_use]
    pub fn total_dissipation_w(self) -> f64 {
        self.density_dissipation_w + self.compressional_dissipation_w
    }
}

/// Complete original pressure coefficients and separate irreversible work.
/// No dispersive stored energy is inferred from a single-frequency sample.
#[derive(Clone, Debug)]
pub struct EquivalentFluidPressureField {
    pub pressure_pa: Vec<C64>,
    /// Total signed work of imposed boundary velocity and volume injection [W],
    /// separate from mechanical wall traction. Includes the two subtotals below.
    pub input_power_w: f64,
    /// Signed work of prepared quadratic boundary sources [W], already included
    /// in input_power_w. Coherent sources may individually absorb power.
    pub boundary_source_power_w: f64,
    /// Signed work of distributed volume sources [W], already included in
    /// input_power_w. This source work is distinct from material absorption.
    pub volume_source_power_w: f64,
    pub boundary_dissipation_w: f64,
    pub density_dissipation_w: f64,
    pub compressional_dissipation_w: f64,
    /// Same order as the original geometry().tetrahedra().
    pub element_power: Vec<EquivalentFluidElementPower>,
    /// Every original wall displacement, velocity and Maxwell memory.
    pub walls: Vec<MechanicalWallPhasors>,
    /// Original wall storage, irreversible heat and signed pressure/traction work.
    /// Bulk dispersive stored energy remains unspecified.
    pub wall_accounting: MechanicalWallAccounting,
    /// Complete fluid/wall equation error, not spatial or material-model error.
    pub backward_error: f64,
}

impl EquivalentFluidPressureField {
    /// Fluid and mechanical source work minus independent boundary, volume and wall heat.
    #[must_use]
    pub fn power_balance_residual_w(&self) -> f64 {
        (Dd::from_f64(self.input_power_w)
            + Dd::from_f64(self.wall_accounting.wall_traction_power_w)
            - Dd::from_f64(self.boundary_dissipation_w)
            - Dd::from_f64(self.density_dissipation_w)
            - Dd::from_f64(self.compressional_dissipation_w)
            - Dd::from_f64(self.wall_accounting.resistance_dissipation_w)
            - Dd::from_f64(self.wall_accounting.relaxation_dissipation_w))
        .to_f64()
    }

    /// Fluid source work minus boundary/volume heat and signed pressure work on walls.
    #[must_use]
    pub fn fluid_power_balance_residual_w(&self) -> f64 {
        (Dd::from_f64(self.input_power_w)
            - Dd::from_f64(self.boundary_dissipation_w)
            - Dd::from_f64(self.density_dissipation_w)
            - Dd::from_f64(self.compressional_dissipation_w)
            - Dd::from_f64(self.wall_accounting.fluid_to_wall_power_w))
        .to_f64()
    }
}

struct ElementLoss {
    dofs: Vec<usize>,
    // Packed upper triangles of the two positive damping forms, BEFORE the
    // peak-to-mean factor 1/2. No dense original-coordinate matrix is retained.
    density: Vec<f64>,
    compression: Vec<f64>,
    volume_m3: f64,
}

/// Full equivalent-fluid pressure operator at one frozen physical frequency.
///
/// The original FE geometry, pressure order, exact pressure-release constraints,
/// boundary forms, equilibrated LU and residual correction are reused. Every
/// material callback runs once per volume point during preparation; subsequent
/// coherent drives and element heat queries do not resample material or geometry.
/// Jumps must conform to element faces. Smooth fields and curved geometry need
/// independent material quadrature, geometry and h/p convergence.
pub struct PreparedEquivalentFluid<'f, 'c> {
    geometry: &'f HelmholtzFe<'c>,
    frequency_hz: f64,
    pencil: CoupledPencil,
    pressure_release: PressureReleaseConstraints,
    boundary_conductance: Csr,
    losses: Vec<ElementLoss>,
    boundaries: Vec<MechanicalBoundary>,
}

/// One applied traction shape bound to this exact fluid/wall preparation.
/// Its amplitude is peak pressure [Pa], positive outward like interior pressure.
/// The cached force uses the original consistent wall trace matrix.
pub struct PreparedEquivalentFluidTraction<'p, 'f, 'c> {
    owner: &'p PreparedEquivalentFluid<'f, 'c>,
    wall: usize,
    coefficients: Vec<f64>,
    column: WallTraction,
}

impl PreparedEquivalentFluidTraction<'_, '_, '_> {
    #[must_use]
    pub const fn wall(&self) -> usize {
        self.wall
    }

    /// Dimensionless coefficients in the original `wall_dofs` order.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }
}

impl<'f, 'c> PreparedEquivalentFluid<'f, 'c> {
    /// Assemble all original pressure coordinates with passive volume material.
    ///
    /// The callback receives the original tet index and SI position. It must
    /// return rho/K from ONE material state at frequency_hz; its errors retain
    /// their type. The geometry owner's homogeneous medium is never sampled.
    /// Robin admittances retain specific outward velocity/pressure units; ideal
    /// pressure release is exact elimination, not an approximate radiation law.
    /// Work/storage admission precedes callbacks. Cancellation publishes nothing.
    pub fn new<E: From<HeterogeneousError>>(
        geometry: &'f HelmholtzFe<'c>,
        frequency_hz: f64,
        admittances: &[BoundaryAdmittance],
        pressure_release_faces: &[[u32; 3]],
        controls: EquivalentFluidControls,
        material_at: &mut impl FnMut(usize, [f64; 3]) -> Result<EquivalentFluidCoefficients, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, E> {
        Self::new_with_coupled_mechanical_walls(
            geometry,
            frequency_hz,
            admittances,
            &[],
            pressure_release_faces,
            controls,
            material_at,
            cancelled,
        )
    }

    /// Retain reciprocal motion of original local mechanical walls in the SAME
    /// complex-material pressure solve, with every Maxwell memory coordinate.
    ///
    /// Robin, mechanical and pressure-release faces must be disjoint. Adjacent
    /// pressure release eliminates no wall coordinates. The shared original
    /// face mass transfers pressure force and normal wall motion reciprocally;
    /// no inverse wall impedance or added damping is used, even at local Z=0.
    /// `controls.harmonic.max_dofs` bounds all free fluid and wall states, and
    /// dense admission uses their complete extent BEFORE material callbacks.
    ///
    /// Material coefficients are frozen at this exact frequency. This method
    /// infers neither a dispersive frequency tangent nor a causal bulk time law.
    /// Walls retain their original locally reacting law, without lateral bending,
    /// moving porous frames, or an inferred exterior radiation impedance.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_coupled_mechanical_walls<E: From<HeterogeneousError>>(
        geometry: &'f HelmholtzFe<'c>,
        frequency_hz: f64,
        admittances: &[BoundaryAdmittance],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        controls: EquivalentFluidControls,
        material_at: &mut impl FnMut(usize, [f64; 3]) -> Result<EquivalentFluidCoefficients, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, E> {
        poll(cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        let limits = controls.harmonic;
        if !omega.is_finite()
            || omega <= 0.0
            || !(omega * omega).is_finite()
            || omega * omega == 0.0
            || geometry.ndof() == 0
            || geometry.tetrahedra().is_empty()
            || limits.max_dofs == 0
            || limits.max_dense_bytes == 0
            || !limits.residual_tolerance.is_finite()
            || limits.residual_tolerance <= 0.0
            || controls.quadrature_order < geometry.order + 2
            || controls.quadrature_order > MAX_QUADRATURE_ORDER
        {
            return Err(HeterogeneousError::from(invalid(
                "invalid equivalent-fluid frequency, geometry or controls",
            ))
            .into());
        }
        let elements = geometry.tetrahedra().len();
        let points = elements
            .checked_mul(controls.quadrature_order.pow(3))
            .ok_or_else(|| {
                HeterogeneousError::from(invalid("equivalent-fluid quadrature extent overflow"))
            })?;
        admit_limit(
            "equivalent-fluid volume quadrature points",
            points,
            controls.max_volume_points,
        )?;
        let order = geometry.order;
        let local_dofs = (order + 1) * (order + 2) * (order + 3) / 6;
        let triangle = local_dofs * (local_dofs + 1) / 2;
        let loss_bytes = elements
            .checked_mul(
                2 * triangle * core::mem::size_of::<f64>()
                    + local_dofs * core::mem::size_of::<usize>()
                    + core::mem::size_of::<ElementLoss>(),
            )
            .ok_or_else(|| {
                HeterogeneousError::from(invalid("equivalent-fluid loss storage overflow"))
            })?;
        admit_limit(
            "equivalent-fluid element loss bytes",
            loss_bytes,
            controls.max_element_loss_bytes,
        )?;
        let cap = limits.max_dofs.min(MAX_HARMONIC_DOFS);
        if pressure_release_faces.is_empty() {
            admit_limit(
                "unconstrained pressure degrees of freedom",
                geometry.ndof(),
                cap,
            )?;
        }
        let pressure_release = geometry
            .pressure_release_constraints(
                pressure_release_faces,
                EmptyFreePressure::Refuse,
                cancelled,
            )
            .map_err(HeterogeneousError::from)?;
        let n = pressure_release.free_dofs.len();
        admit_limit("unconstrained pressure degrees of freedom", n, cap)?;
        let bytes = n
            .checked_mul(n)
            .and_then(|value| value.checked_mul(32))
            .ok_or_else(|| {
                HeterogeneousError::from(invalid("equivalent-fluid dense extent overflow"))
            })?;
        admit_limit(
            "full pressure dense matrix bytes",
            bytes,
            limits.max_dense_bytes,
        )?;
        let mut seen = pressure_release.faces.clone();
        for faces in admittances
            .iter()
            .map(|b| &b.faces)
            .chain(boundaries.iter().map(|b| &b.faces))
        {
            for &face in faces {
                poll(cancelled)?;
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(HeterogeneousError::from(invalid(
                        "mechanical, admittance and pressure-release faces must be disjoint",
                    ))
                    .into());
                }
            }
        }
        let (boundary_conductance, susceptance) =
            HeterogeneousHelmholtzFe::boundary_coefficient_matrices_for(
                geometry,
                admittances
                    .iter()
                    .map(|b| (b.faces.as_slice(), b.admittance_m_pa_s)),
                true,
                cancelled,
            )?;
        let (walls, states) = mechanical_walls_for(geometry, boundaries, n, cap, cancelled)?;
        let complete_bytes = states
            .checked_mul(states)
            .and_then(|value| value.checked_mul(32))
            .ok_or_else(|| {
                HeterogeneousError::from(invalid("equivalent-fluid complete dense extent overflow"))
            })?;
        admit_limit(
            "complete harmonic dense matrix bytes",
            complete_bytes,
            limits.max_dense_bytes,
        )?;
        let quadrature = duffy_quadrature(controls.quadrature_order);
        let mut stiffness = Coo::new(geometry.ndof(), geometry.ndof());
        let mut mass = Coo::new(geometry.ndof(), geometry.ndof());
        let mut conductance = Coo::new(geometry.ndof(), geometry.ndof());
        let mut losses = Vec::with_capacity(elements);
        for element in 0..elements {
            poll(cancelled)?;
            let dofs = geometry.space.element_dofs(element);
            let count = dofs.len();
            let triangle = count * (count + 1) / 2;
            let mut local_k = vec![0.0; triangle];
            let mut local_m = vec![0.0; triangle];
            let mut density = vec![0.0; triangle];
            let mut compression = vec![0.0; triangle];
            let mut values = vec![0.0; count];
            let mut gradients = vec![[0.0; 3]; count];
            let mut volume = Dd::ZERO;
            for &(barycentric, weight) in &quadrature {
                poll(cancelled)?;
                let point = geometry
                    .geometry_point(element, barycentric)
                    .map_err(HeterogeneousError::from)?;
                let (alpha, beta) = material_at(element, point.position_m)?.inverse()?;
                let measure = weight * point.determinant_m3;
                let coefficients = [
                    measure * alpha.re,
                    measure * beta.re,
                    measure * (alpha.im / omega),
                    measure * (-beta.im * omega),
                ];
                if coefficients.iter().any(|v| !v.is_finite())
                    || coefficients[0] <= 0.0
                    || coefficients[1] <= 0.0
                {
                    return Err(HeterogeneousError::from(invalid(
                        "equivalent-fluid quadrature coefficient overflow",
                    ))
                    .into());
                }
                volume = volume + Dd::from_f64(measure);
                for (i, (_, shape)) in dofs.iter().enumerate() {
                    values[i] = shape.eval(barycentric, order);
                    gradients[i] = [0.0; 3];
                    for (a, derivative) in shape.d_lambda(barycentric, order).iter().enumerate() {
                        for axis in 0..3 {
                            gradients[i][axis] += derivative * point.barycentric_gradients[a][axis];
                        }
                    }
                }
                let mut index = 0;
                for i in 0..count {
                    for j in i..count {
                        let gradient: f64 = gradients[i]
                            .iter()
                            .zip(&gradients[j])
                            .map(|(a, b)| a * b)
                            .sum();
                        let shape = values[i] * values[j];
                        local_k[index] += coefficients[0] * gradient;
                        local_m[index] += coefficients[1] * shape;
                        density[index] += coefficients[2] * gradient;
                        compression[index] += coefficients[3] * shape;
                        index += 1;
                    }
                }
            }
            let mut index = 0;
            for (i, &(row, _)) in dofs.iter().enumerate() {
                for (j, &(col, _)) in dofs.iter().enumerate().skip(i) {
                    let g = density[index] + compression[index];
                    for (matrix, value) in [
                        (&mut stiffness, local_k[index]),
                        (&mut mass, local_m[index]),
                        (&mut conductance, g),
                    ] {
                        if !value.is_finite() {
                            return Err(HeterogeneousError::from(invalid(
                                "equivalent-fluid element assembly overflow",
                            ))
                            .into());
                        }
                        matrix.push(row, col, value);
                        if i != j {
                            matrix.push(col, row, value);
                        }
                    }
                    index += 1;
                }
            }
            let volume_m3 = volume.to_f64();
            if !volume_m3.is_finite() || volume_m3 <= 0.0 {
                return Err(
                    HeterogeneousError::from(invalid("equivalent-fluid volume overflow")).into(),
                );
            }
            losses.push(ElementLoss {
                dofs: dofs.iter().map(|&(dof, _)| dof).collect(),
                density,
                compression,
                volume_m3,
            });
        }
        for row in 0..geometry.ndof() {
            poll(cancelled)?;
            let (cols, values) = boundary_conductance.row(row);
            for (&col, &value) in cols.iter().zip(values) {
                conductance.push(row, col, value);
            }
        }
        let (stiffness, mass, conductance) = (
            stiffness.assemble(),
            mass.assemble(),
            conductance.assemble(),
        );
        let pencil = CoupledPencil::prepare(
            omega,
            [&stiffness, &mass, &conductance, &susceptance],
            &pressure_release.free_dofs,
            walls,
            MechanicalHarmonicLimits {
                max_states: cap,
                max_dense_bytes: limits.max_dense_bytes,
                residual_tolerance: limits.residual_tolerance,
            },
            cancelled,
        )
        .map_err(HeterogeneousError::from)?;
        poll(cancelled)?;
        Ok(Self {
            geometry,
            frequency_hz,
            pencil,
            pressure_release,
            boundary_conductance,
            losses,
            boundaries: boundaries.to_vec(),
        })
    }

    #[must_use]
    pub const fn geometry(&self) -> &'f HelmholtzFe<'c> {
        self.geometry
    }

    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }

    #[must_use]
    pub fn wall_count(&self) -> usize {
        self.boundaries.len()
    }

    #[must_use]
    pub fn wall(&self, wall: usize) -> Option<&MechanicalBoundary> {
        self.boundaries.get(wall)
    }

    /// Original sorted trace coordinates, including pressure-constrained vertices.
    #[must_use]
    pub fn wall_dofs(&self, wall: usize) -> Option<&[usize]> {
        self.pencil.wall_dofs(wall)
    }

    /// Recompute separate positive resistance and Maxwell absorption in original
    /// wall order, with each patch's complete physical area and surface mass.
    ///
    /// The coherent wall phasors must belong to THIS geometry, wall order and
    /// frequency; the caller owns that binding, as for `element_dissipation`.
    /// Native relaxed extensions are retained directly so tiny physical heat
    /// survives cancellation in `sqrt(k)*displacement - memory`. Complete state
    /// dimensions and finite values are checked. Pressure exchange, boundary
    /// Robin loss, bulk heat and aggregate power receipts never enter these forms.
    pub fn wall_dissipation(
        &self,
        walls: &[MechanicalWallPhasors],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<MechanicalWallPower>, HeterogeneousError> {
        Ok(self.pencil.wall_dissipation(walls, cancelled)?)
    }

    /// Project an applied mechanical traction shape once on this exact owner.
    /// The force drives wall momentum; it does not prescribe wall velocity.
    pub fn prepare_wall_traction(
        &self,
        profile: &WallTractionProfile,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedEquivalentFluidTraction<'_, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let dofs = self
            .wall_dofs(profile.wall)
            .ok_or_else(|| invalid("traction wall index"))?;
        let coefficients = profile.local_coefficients(dofs, cancelled)?;
        let column = self
            .pencil
            .prepare_traction(profile.wall, &coefficients, cancelled)?;
        Ok(PreparedEquivalentFluidTraction {
            owner: self,
            wall: profile.wall,
            coefficients,
            column,
        })
    }

    /// Solve simultaneous original linear normal-velocity traces [m/s].
    /// Sources interfere in the same pressure field before work is evaluated.
    /// The retained factor and original residual corrections are reused.
    pub fn solve(
        &self,
        inward_velocity: &[VelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EquivalentFluidPressureField, HeterogeneousError> {
        self.solve_with_wall_tractions(inward_velocity, &[], cancelled)
    }

    /// Solve coherent fluid velocity [m/s] and independently applied wall
    /// traction [Pa] through the same complete pressure/motion/memory factor.
    /// The returned fluid, wall and total balances retain their signed work.
    pub fn solve_with_wall_tractions(
        &self,
        inward_velocity: &[VelocityTrace],
        tractions: &[(&PreparedEquivalentFluidTraction<'_, 'f, 'c>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EquivalentFluidPressureField, HeterogeneousError> {
        self.solve_with_sources(inward_velocity, &[], &[], tractions, cancelled)
    }

    /// Resolve simultaneous linear boundary velocity, original prepared
    /// quadratic boundary motion, distributed volume injection and independent
    /// mechanical wall traction through the SAME complete fluid/wall factor.
    ///
    /// Prepared boundary amplitudes are inward velocity [m/s], volume amplitudes
    /// are volume flow [m^3/s], and wall traction is outward pressure [Pa]. All
    /// are complex peak phasors in exp(+i omega t). Their original load columns
    /// add before solving, retaining coherent interference and signed source
    /// work without normalization, clipping or material/source resampling.
    /// Every prepared source must belong to THIS exact owner, even at zero
    /// amplitude. Boundary sources cannot overlap pressure-release faces.
    /// The original source-count, complete-state and cancellation limits apply.
    pub fn solve_with_sources(
        &self,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        tractions: &[(&PreparedEquivalentFluidTraction<'_, 'f, 'c>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EquivalentFluidPressureField, HeterogeneousError> {
        poll(cancelled)?;
        super::super::source::admit_count(tractions.len())?;
        let mut columns = Vec::with_capacity(tractions.len());
        for &(source, amplitude) in tractions {
            poll(cancelled)?;
            if !core::ptr::eq(source.owner, self) || !finite(amplitude) {
                return Err(invalid(
                    "finite traction amplitude and original prepared equivalent-fluid owner required",
                ).into());
            }
            columns.push((&source.column, amplitude));
        }
        let load = if boundary_sources.is_empty() && volume_sources.is_empty() {
            // Keep the original unprepared-drive path and its arithmetic.
            PreparedHeterogeneous::admit_velocity_faces_for(
                &self.pressure_release,
                inward_velocity,
                cancelled,
            )?;
            HeterogeneousHelmholtzFe::velocity_load_for(self.geometry, inward_velocity, cancelled)?
        } else {
            PreparedHeterogeneous::project_sources_for_owner(
                PressureSourceOwner::Equivalent(self),
                &self.pressure_release,
                inward_velocity,
                boundary_sources,
                volume_sources,
                cancelled,
            )?
        };
        let (mut solved, traction_power) = self
            .pencil
            .solve_with_tractions(&load, &columns, cancelled)?;
        let (walls, mut wall_accounting) = self.pencil.wall_fields(&solved.pressure, cancelled)?;
        wall_accounting.wall_traction_power_w = traction_power;
        solved.pressure.truncate(self.geometry.ndof());
        let mut work = Dd::ZERO;
        for (pressure, source) in solved.pressure.iter().zip(&load) {
            work = work
                + Dd::from_f64(pressure.re) * Dd::from_f64(source.re)
                + Dd::from_f64(pressure.im) * Dd::from_f64(source.im);
        }
        let input_power_w = 0.5 * work.to_f64();
        let (boundary_source_power_w, volume_source_power_w) = PreparedHeterogeneous::source_power(
            &solved.pressure,
            boundary_sources,
            volume_sources,
            cancelled,
        )?;
        let boundary_dissipation_w =
            0.5 * quadratic_nonnegative(&self.boundary_conductance, &solved.pressure, cancelled)?;
        let element_power = self.element_dissipation(&solved.pressure, cancelled)?;
        let mut density = Dd::ZERO;
        let mut compression = Dd::ZERO;
        for power in &element_power {
            density = density + Dd::from_f64(power.density_dissipation_w);
            compression = compression + Dd::from_f64(power.compressional_dissipation_w);
        }
        let density_dissipation_w = density.to_f64();
        let compressional_dissipation_w = compression.to_f64();
        if [
            input_power_w,
            boundary_dissipation_w,
            density_dissipation_w,
            compressional_dissipation_w,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("equivalent-fluid work overflow").into());
        }
        poll(cancelled)?;
        Ok(EquivalentFluidPressureField {
            pressure_pa: solved.pressure,
            input_power_w,
            boundary_source_power_w,
            volume_source_power_w,
            boundary_dissipation_w,
            density_dissipation_w,
            compressional_dissipation_w,
            element_power,
            walls,
            wall_accounting,
            backward_error: solved.backward_error,
        })
    }

    /// Re-evaluate original positive volume-loss forms from pressure coefficients.
    ///
    /// This also accepts externally composed pressure fields on THIS original
    /// space and frequency; the caller owns that binding. Surface losses never
    /// enter these rows. No source/material query, numerical gradient, net-power
    /// subtraction or heat deposition is performed. Cancellation returns no rows.
    pub fn element_dissipation(
        &self,
        pressure_pa: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<EquivalentFluidElementPower>, HeterogeneousError> {
        poll(cancelled)?;
        if pressure_pa.len() != self.geometry.ndof() || pressure_pa.iter().any(|&p| !finite(p)) {
            return Err(invalid(
                "equivalent-fluid heat requires original finite pressure coefficients",
            )
            .into());
        }
        let mut result = Vec::with_capacity(self.losses.len());
        for loss in &self.losses {
            poll(cancelled)?;
            result.push(EquivalentFluidElementPower {
                volume_m3: loss.volume_m3,
                density_dissipation_w: 0.5 * local_form(&loss.density, &loss.dofs, pressure_pa)?,
                compressional_dissipation_w: 0.5
                    * local_form(&loss.compression, &loss.dofs, pressure_pa)?,
            });
        }
        poll(cancelled)?;
        Ok(result)
    }
}

fn admit_limit(
    what: &'static str,
    requested: usize,
    maximum: usize,
) -> Result<(), HeterogeneousError> {
    if requested > maximum {
        Err(HeterogeneousError::Limit {
            what,
            requested,
            maximum,
        })
    } else {
        Ok(())
    }
}

fn local_form(matrix: &[f64], dofs: &[usize], pressure: &[C64]) -> Result<f64, HeterogeneousError> {
    let mut sum = Dd::ZERO;
    let mut absolute = 0.0;
    let mut index = 0;
    for (i, &row) in dofs.iter().enumerate() {
        let p = pressure[row];
        for (j, &col) in dofs.iter().enumerate().skip(i) {
            let q = pressure[col];
            let coefficient = matrix[index] * if i == j { 1.0 } else { 2.0 };
            sum = sum
                + Dd::from_f64(coefficient)
                    * (Dd::from_f64(p.re) * Dd::from_f64(q.re)
                        + Dd::from_f64(p.im) * Dd::from_f64(q.im));
            absolute += coefficient.abs() * p.abs() * q.abs();
            index += 1;
        }
    }
    let value = sum.to_f64();
    if !value.is_finite() || !absolute.is_finite() || value < -64.0 * f64::EPSILON * absolute {
        return Err(
            invalid("equivalent-fluid local dissipation form is nonfinite or negative").into(),
        );
    }
    Ok(value.max(0.0))
}
