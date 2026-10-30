//! Harmonic duct propagation with the material owner's molecular reservoirs.
//!
//! Molecular relaxation changes both storage and loss. With `s = -i omega`,
//! the original reservoir sum supplies `cp(s)`. The circular wall's original
//! positive thermal measure `(w_n, lambda_n)` then gives
//!
//! `kappa(s) = [1 - R sum_n w_n s/(s cp(s) + cp_eq lambda_n)]/(rho R T)`.
//!
//! Without a wall the sum is `1/cp(s)`. This is the SAME thermodynamic law
//! used by [`crate::dg::DgDuct::with_molecular_relaxation`], evaluated directly
//! in frequency instead of eliminating a time-domain state matrix. Independent
//! molecular and equilibrium thermal shunts would double-count storage and
//! miss the isothermal limit. No ISO or Stokes--Kirchhoff loss is added.
//!
//! Geometry, positive circular diffusion spectra, spherical-cone wave bases,
//! Magnus integration and scattering composition retain their existing owners.
//! Preparation caches each distinct local radius's diffusion eigensolve;
//! frequency sweeps only evaluate the retained rational coefficients.
//!
//! This is frozen, linear plane-wave acoustics with rigid walls. The homogeneous
//! constructor uses the existing local-radius spherical-cone approximation;
//! profiles and prescribed gas fields use the original Webster/Magnus owner.
//! Molecular data, diffusion poles, cone slices and profile/field integration
//! steps are independent resolution/model choices. The diffusion admission is
//! against equilibrium Bessel response, not an error bound for the joint
//! molecular law. Compact tone-hole chimneys consume the same reservoirs and
//! wall spectra; each open mouth requires an explicit pressure/flow load.
//! Their original inner/series/mutual corrections and vent-admittance mixture
//! remain a compact model without resolved pads. Scalar mouth loads or a
//! complete exterior impedance matrix close the same propagation operator.
//! Exterior radiation is explicit, never inferred from the bore gas.

use crate::compose::{self, ChainResult, ComposeOptions, PreparedChain, TwoPort};
use crate::{ConeSlices, Duct, DuctError, LossModel, Segment, SegmentWave};
use fs_material::MaterialError;
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, det};
use fs_phs::ZkDiffusive;
use std::collections::BTreeMap;

mod holes;
pub use holes::{MolecularHoleFlow, MolecularHoleSite, MolecularResponse};
mod field;
mod exterior;
pub use exterior::{MolecularExteriorResponse, MolecularPortScattering};

/// Original constitutive or propagation refusal, without a fallback gas law.
#[derive(Debug)]
pub enum MolecularDuctError {
    Material(MaterialError),
    Duct(DuctError),
}
impl core::fmt::Display for MolecularDuctError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "molecular duct: {error}"),
            Self::Duct(error) => write!(f, "molecular duct: {error}"),
        }
    }
}
impl std::error::Error for MolecularDuctError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Duct(error) => Some(error),
        }
    }
}
impl From<MaterialError> for MolecularDuctError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<DuctError> for MolecularDuctError {
    fn from(error: DuctError) -> Self {
        Self::Duct(error)
    }
}

/// Frequency-independent wall-diffusion resolution at every physical radius.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MolecularWallResolution {
    /// Positive continued-fraction arms for EACH of shear and heat diffusion.
    pub poles_per_field: usize,
    /// Admitted physical angular-frequency band [rad/s].
    pub band_rad_s: [f64; 2],
    /// Sampled equilibrium complex AND dissipative wall-error ceiling.
    pub maximum_relative_error: f64,
}

#[derive(Debug, Clone)]
struct MaterialPoint {
    radius: f64,
    x_m: f64,
    wall: Option<ZkDiffusive>,
    /// Absent on the original homogeneous path; otherwise the exact source
    /// state at this geometry/material quadrature point.
    state: Option<usize>,
}
#[derive(Debug, Clone)]
enum Piece {
    Wave { geometry: Segment, point: usize },
    Magnus { h: f64, points: [usize; 2] },
    ToneHole(holes::ToneHole),
}

/// Prepared harmonic counterpart of the existing molecular DG duct.
#[derive(Debug, Clone)]
pub struct MolecularDuct {
    state: RelaxationState,
    field_binding: Option<field::FieldBinding>,
    points: Vec<MaterialPoint>,
    pieces: Vec<Piece>,
    reference_impedance: f64,
}

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(cancelled: &mut impl FnMut() -> bool, completed: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled {
            completed_steps: completed,
        })
    } else {
        Ok(())
    }
}

impl MolecularDuct {
    /// Bind the original gas and reservoirs to duct and tone-hole geometry.
    /// `None` selects adiabatic, slip walls; `Some` retains coupled shear and
    /// thermal diffusion. Cone slices and each profile's own RiccatiSteps are
    /// preserved, and `maximum_blocks` is an admission cap, never a clamp.
    pub fn new(
        duct: &Duct,
        state: RelaxationState,
        wall: Option<MolecularWallResolution>,
        cones: ConeSlices,
        maximum_blocks: usize,
    ) -> Result<Self, MolecularDuctError> {
        Self::new_controlled(duct, state, wall, cones, maximum_blocks, &mut || false)
    }

    /// Cancel between original diffusion solves; no partial operator escapes.
    pub fn new_controlled(
        duct: &Duct,
        state: RelaxationState,
        wall: Option<MolecularWallResolution>,
        cones: ConeSlices,
        maximum_blocks: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MolecularDuctError> {
        poll(cancelled, 0)?;
        crate::validate_gas(state.gas())?;
        let first = duct.segments.first().ok_or(DuctError::EmptyDuct)?;
        if maximum_blocks == 0 || maximum_blocks > compose::MAX_COMPOSITION_BLOCKS {
            return Err(bad("molecular duct composition block cap is invalid").into());
        }
        // Admit the entire geometry/work request before any eigensolve.
        let mut count = 0usize;
        for segment in &duct.segments {
            poll(cancelled, 0)?;
            segment.validate()?;
            let n = match segment {
                Segment::Cylinder { .. } => 1,
                Segment::Cone {
                    inlet_radius,
                    outlet_radius,
                    ..
                } => {
                    let n = cones.count(*inlet_radius, *outlet_radius)?;
                    if wall.is_some() { n } else { 1 }
                }
                Segment::Profile { steps, .. } => steps.0,
                Segment::ToneHole { state, .. } => {
                    holes::opening(*state)?;
                    3 // two series halves surrounding the projective shunt
                }
            };
            count = count
                .checked_add(n)
                .ok_or_else(|| bad("molecular duct block count overflow"))?;
            if count > maximum_blocks {
                return Err(bad("molecular duct composition block work cap exceeded").into());
            }
        }
        let gas = state.gas();
        let r = first.inlet_radius();
        let reference = gas.characteristic_impedance / (core::f64::consts::PI * r * r);
        if !positive(reference) {
            return Err(bad("molecular reference impedance is unrepresentable").into());
        }
        let mut result = Self {
            state,
            field_binding: None,
            points: Vec::new(),
            pieces: Vec::with_capacity(count),
            reference_impedance: reference,
        };
        let mut cache = BTreeMap::new();
        let open_sites = crate::collect_open_hole_sites(duct);
        let mut origin = 0.0;
        let root3 = det::sqrt(3.0);
        for segment in &duct.segments {
            match segment {
                Segment::Cylinder { radius, .. } => {
                    let point = result.point(*radius, origin, wall, &mut cache, cancelled)?;
                    result.pieces.push(Piece::Wave {
                        geometry: segment.clone(),
                        point,
                    });
                }
                Segment::Cone {
                    inlet_radius: a,
                    outlet_radius: b,
                    length,
                } => {
                    let n = if wall.is_some() {
                        cones.count(*a, *b)?
                    } else {
                        1
                    };
                    for j in 0..n {
                        let rin = a + (b - a) * j as f64 / n as f64;
                        let rout = a + (b - a) * (j + 1) as f64 / n as f64;
                        let radius = f64::midpoint(rin, rout);
                        let point = result.point(
                            radius,
                            origin + length * (j as f64 + 0.5) / n as f64,
                            wall,
                            &mut cache,
                            cancelled,
                        )?;
                        result.pieces.push(Piece::Wave {
                            geometry: Segment::Cone {
                                inlet_radius: rin,
                                outlet_radius: rout,
                                length: length / n as f64,
                            },
                            point,
                        });
                    }
                }
                Segment::Profile { profile, steps } => {
                    let h = -profile.length() / steps.0 as f64;
                    for j in 0..steps.0 {
                        let right = -((j + 1) as f64) * h;
                        let mut points = [0; 2];
                        for (i, c) in [0.5 - root3 / 6.0, 0.5 + root3 / 6.0]
                            .into_iter()
                            .enumerate()
                        {
                            let x = right + c * h;
                            points[i] = result.point(
                                profile.radius(x)?,
                                origin + x,
                                wall,
                                &mut cache,
                                cancelled,
                            )?;
                        }
                        result.pieces.push(Piece::Magnus { h, points });
                    }
                }
                Segment::ToneHole { hole_radius, .. } => {
                    let point = result.point(*hole_radius, origin, wall, &mut cache, cancelled)?;
                    result.pieces.push(Piece::ToneHole(holes::ToneHole::new(
                        segment,
                        point,
                        origin,
                        &open_sites,
                    )?));
                }
            }
            origin += segment.axial_length();
        }
        Ok(result)
    }

    fn point(
        &mut self,
        radius: f64,
        x: f64,
        resolution: Option<MolecularWallResolution>,
        cache: &mut BTreeMap<u64, usize>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize, MolecularDuctError> {
        poll(cancelled, self.pieces.len())?;
        if let Some(&index) = cache.get(&radius.to_bits()) {
            return Ok(index);
        }
        let gas = self.state.gas();
        let wall = resolution
            .map(|r| {
                ZkDiffusive::new(
                    radius,
                    gas.density,
                    gas.dynamic_viscosity,
                    gas.gamma,
                    gas.prandtl,
                    r.poles_per_field,
                    r.band_rad_s,
                    r.maximum_relative_error,
                )
                .map_err(|reason| DuctError::WallDiffusionUnresolved {
                    reason,
                    radius_m: radius,
                    x_m: x,
                })
            })
            .transpose()?;
        let index = self.points.len();
        self.points.push(MaterialPoint {
            radius,
            x_m: x,
            wall,
            state: None,
        });
        cache.insert(radius.to_bits(), index);
        Ok(index)
    }

    /// Original evaluated inlet material, never an averaged/effective state.
    /// A field binding retains its local states at each material point.
    pub fn state(&self) -> &RelaxationState {
        &self.state
    }
    /// Exact number of retained scattering blocks, including three per hole.
    pub fn blocks(&self) -> usize {
        self.pieces
            .iter()
            .map(|piece| match piece {
                Piece::ToneHole(_) => 3,
                _ => 1,
            })
            .sum()
    }
    /// Distinct radius/material pairs prepared once and reused.
    pub fn material_points(&self) -> usize {
        self.points.len()
    }

    fn state_at_point(&self, index: usize) -> &RelaxationState {
        match self.points[index].state {
            Some(state) => &self.field_binding.as_ref().expect("field material point").states[state],
            None => &self.state,
        }
    }

    fn wave(
        &self,
        index: usize,
        omega: f64,
        cp: C64,
        relaxing_k: C64,
    ) -> Result<SegmentWave, MolecularDuctError> {
        let point = &self.points[index];
        let state = self.state_at_point(index);
        let gas = state.gas();
        let area = core::f64::consts::PI * point.radius * point.radius;
        if point.wall.is_none() {
            if state.processes().is_empty() {
                return Ok(crate::segment_wave(
                    gas,
                    point.radius,
                    omega,
                    LossModel::Lossless,
                )?);
            }
            let k = relaxing_k;
            let z = C64::from_re(omega * gas.density) / k;
            return Ok(SegmentWave {
                wavenumber: k,
                characteristic_impedance: z.scale(1.0 / area),
                specific_impedance: z,
                shear_number: f64::INFINITY,
            });
        }
        let model = point.wall.as_ref().expect("checked above");
        let (series, equilibrium_shunt) =
            model
                .response_at(omega)
                .map_err(|reason| DuctError::WallDiffusionUnresolved {
                    reason,
                    radius_m: point.radius,
                    x_m: point.x_m,
                })?;
        let z_plus = series.scale(gas.density / area);
        let y_plus = if state.processes().is_empty() {
            equilibrium_shunt.scale(area / (gas.density * gas.sound_speed * gas.sound_speed))
        } else {
            let mut temperature = C64::ZERO;
            for term in model.shunt_terms() {
                let weight = term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0));
                // s/(s cp+d) = 1/(cp+i d/omega), avoiding large s*cp.
                let denominator =
                    cp + C64::new(0.0, gas.specific_heat_cp * (term.rate_per_s / omega));
                temperature = temperature + denominator.recip().scale(weight);
            }
            let kappa = (C64::ONE - temperature.scale(gas.specific_gas_constant))
                .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature));
            (C64::new(0.0, -omega) * kappa.scale(area)).conj()
        };
        let (k, zc) = crate::helmholtz_wave_from_telegraph(z_plus, y_plus);
        if !finite(k) || k.re <= 0.0 || k.im < 0.0 || !finite(zc) || zc.re <= 0.0 {
            return Err(bad("molecular duct wave is not finite, forward and passive").into());
        }
        Ok(SegmentWave {
            wavenumber: k,
            characteristic_impedance: zc,
            specific_impedance: zc.scale(area),
            shear_number: point.radius * det::sqrt(gas.density * omega / gas.dynamic_viscosity),
        })
    }

    /// Retain the solved two-ports for existing graph or scattering consumers.
    /// Cancellation is checked at each point and each axial block.
    pub fn prepared_chain_controlled(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedChain, MolecularDuctError> {
        self.prepare_with_holes_controlled(omega, &[], ComposeOptions::default(), cancelled)
            .map(|(chain, _)| chain)
    }

    fn prepare_with_holes_controlled(
        &self,
        omega: f64,
        mouth_loads: &[C64],
        options: ComposeOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(PreparedChain, Vec<holes::PreparedHole>), MolecularDuctError> {
        poll(cancelled, 0)?;
        options.validate()?;
        if !positive(omega) {
            return Err(bad("molecular duct frequency must be positive and finite").into());
        }
        if self.blocks() > options.maximum_blocks {
            return Err(bad("molecular duct composition block work cap exceeded").into());
        }
        let open_count = self
            .pieces
            .iter()
            .filter(|piece| matches!(piece, Piece::ToneHole(hole) if hole.site.open_fraction > 0.0))
            .count();
        if mouth_loads.len() != open_count {
            return Err(bad("each open molecular chimney requires one explicit mouth load").into());
        }
        if mouth_loads
            .iter()
            .any(|&load| !finite(load) || load.re < 0.0)
        {
            return Err(bad("molecular chimney loads must be finite and passive").into());
        }
        let (waves, min_shear) = self.waves_controlled(omega, cancelled)?;
        let reference = self.reference_impedance;
        let mut blocks = Vec::with_capacity(self.blocks());
        let mut observations = Vec::new();
        let mut loads = mouth_loads.iter();
        for (i, piece) in self.pieces.iter().enumerate() {
            poll(cancelled, i)?;
            let block = match piece {
                Piece::ToneHole(hole) => {
                    let load = if hole.site.open_fraction > 0.0 {
                        Some(*loads.next().expect("admitted load count"))
                    } else {
                        None
                    };
                    let (local, observation) = hole.prepare(
                        &waves[hole.point],
                        self.state_at_point(hole.point).gas().density,
                        omega,
                        load,
                        reference,
                        options,
                        blocks.len(),
                    )?;
                    blocks.extend(local);
                    observations.push(observation);
                    continue;
                }
                _ => self.piece_wave_block(piece, &waves)?,
            };
            blocks.push(block);
        }
        Ok((
            PreparedChain {
                blocks,
                reference_impedance: reference,
                min_shear_number: min_shear,
            },
            observations,
        ))
    }

    fn waves_controlled(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<SegmentWave>, f64), MolecularDuctError> {
        poll(cancelled, 0)?;
        if !positive(omega) {
            return Err(bad("molecular duct frequency must be positive and finite").into());
        }
        let caloric = |state: &RelaxationState| -> Result<(C64, C64), MolecularDuctError> {
            let response = state.caloric_response(omega)?;
            let k = response.bulk_modulus_pa.scale(1.0 / state.gas().density)
                .sqrt().recip().scale(omega);
            if !finite(k) || k.re <= 0.0 || k.im < 0.0 {
                return Err(bad("molecular duct material wavenumber is unrepresentable").into());
            }
            Ok((response.cp_j_kg_k, k))
        };
        // Preserve one evaluation and the original operation order for a
        // homogeneous duct. Varying fields evaluate each cached local state.
        let uniform = caloric(&self.state)?;
        let mut local = Vec::new();
        if let Some(binding) = &self.field_binding {
            local.reserve(binding.states.len());
            for (index, state) in binding.states.iter().enumerate() {
                poll(cancelled, 0)?;
                local.push(if index == 0 { uniform } else { caloric(state)? });
            }
        }
        let mut waves = Vec::with_capacity(self.points.len());
        let mut min_shear = f64::INFINITY;
        for (index, point) in self.points.iter().enumerate() {
            poll(cancelled, 0)?;
            let (cp, k) = match point.state {
                Some(state) => local[state],
                None => uniform,
            };
            let wave = self.wave(index, omega, cp, k)?;
            min_shear = min_shear.min(wave.shear_number);
            waves.push(wave);
        }
        Ok((waves, min_shear))
    }

    fn piece_wave_block(
        &self,
        piece: &Piece,
        waves: &[SegmentWave],
    ) -> Result<TwoPort, MolecularDuctError> {
        let reference = self.reference_impedance;
        Ok(match piece {
            Piece::Wave { geometry, point } => {
                TwoPort::from_wave(geometry, &waves[*point], reference)?
            }
            Piece::Magnus { h, points } => {
                let generator = |point: usize| {
                    let wave = waves[point];
                    let ik = C64::new(0.0, 1.0) * wave.wavenumber;
                    (
                        ik * wave.characteristic_impedance.scale(1.0 / reference),
                        ik * wave.characteristic_impedance.recip().scale(reference),
                    )
                };
                let mut matrix = crate::profile::magnus_step(
                    generator(points[0]), generator(points[1]), *h, det::sqrt(3.0),
                )?;
                matrix[1] = matrix[1].scale(reference);
                matrix[2] = matrix[2].scale(1.0 / reference);
                TwoPort::from_abcd(matrix, reference)?
            }
            Piece::ToneHole(_) => return Err(bad("tone holes require their branch preparation").into()),
        })
    }

    /// Actual terminated pressure/volume-flow response; `None` is a rigid cap,
    /// `Some(0)` pressure release. A radiating exterior supplies its own load
    /// in Pa s/m^3; it is never inferred from the interior relaxation model.
    pub fn response(
        &self,
        omega: f64,
        load: Option<C64>,
        options: ComposeOptions,
    ) -> Result<ChainResult, MolecularDuctError> {
        self.response_controlled(omega, load, options, &mut || false)
    }

    /// Original composition cancellation and conditioning limits apply.
    pub fn response_controlled(
        &self,
        omega: f64,
        load: Option<C64>,
        options: ComposeOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChainResult, MolecularDuctError> {
        if self.blocks() > options.maximum_blocks {
            return Err(bad("molecular duct composition block work cap exceeded").into());
        }
        let (chain, _) = self.prepare_with_holes_controlled(omega, &[], options, cancelled)?;
        Ok(compose::compose_chain_controlled(
            &chain.blocks,
            load,
            options,
            cancelled,
        )?)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod holes_tests;
